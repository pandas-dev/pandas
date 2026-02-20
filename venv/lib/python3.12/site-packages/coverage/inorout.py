# Licensed under the Apache License: http://www.apache.org/licenses/LICENSE-2.0
# For details: https://github.com/coveragepy/coveragepy/blob/main/NOTICE.txt

"""Determining whether files are being measured/reported or not."""

from __future__ import annotations

import importlib.util
import inspect
import itertools
import os
import os.path
import sys
import sysconfig
import traceback
from collections.abc import Iterable
from dataclasses import dataclass
from types import FrameType, ModuleType
from typing import TYPE_CHECKING, Any, cast

from coverage import env
from coverage.disposition import FileDisposition, disposition_init
from coverage.exceptions import ConfigError, CoverageException, PluginError
from coverage.files import (
    GlobMatcher,
    ModuleMatcher,
    TreeMatcher,
    canonical_filename,
    find_python_files,
    prep_patterns,
)
from coverage.misc import isolate_module, sys_modules_saved
from coverage.python import source_for_file, source_for_morf
from coverage.types import TDebugCtl, TFileDisposition, TMorf, TWarnFn

if TYPE_CHECKING:
    from coverage.config import CoverageConfig
    from coverage.plugin_support import Plugins


os = isolate_module(os)


def canonical_path(morf: TMorf, directory: bool = False) -> str:
    """Return the canonical path of the module or file `morf`.

    If the module is a package, then return its directory. If it is a
    module, then return its file, unless `directory` is True, in which
    case return its enclosing directory.

    """
    morf_path = canonical_filename(source_for_morf(morf))
    if morf_path.endswith("__init__.py") or directory:
        morf_path = os.path.split(morf_path)[0]
    return morf_path


def name_for_module(filename: str, frame: FrameType | None) -> str | None:
    """Get the name of the module for a filename and frame.

    For configurability's sake, we allow __main__ modules to be matched by
    their importable name.

    If loaded via runpy (aka -m), we can usually recover the "original"
    full dotted module name, otherwise, we resort to interpreting the
    file name to get the module's name.  In the case that the module name
    can't be determined, None is returned.

    """
    module_globals = frame.f_globals if frame is not None else {}
    dunder_name: str | None = module_globals.get("__name__", None)

    if isinstance(dunder_name, str) and dunder_name != "__main__":
        # This is the usual case: an imported module.
        return dunder_name

    spec = module_globals.get("__spec__", None)
    if spec:
        fullname = spec.name
        if isinstance(fullname, str) and fullname != "__main__":
            # Module loaded via: runpy -m
            return fullname

    # Script as first argument to Python command line.
    inspectedname = inspect.getmodulename(filename)
    if inspectedname is not None:
        return inspectedname
    else:
        return dunder_name


def module_is_namespace(mod: ModuleType) -> bool:
    """Is the module object `mod` a PEP420 namespace module?"""
    return hasattr(mod, "__path__") and getattr(mod, "__file__", None) is None


def module_has_file(mod: ModuleType) -> bool:
    """Does the module object `mod` have an existing __file__ ?"""
    mod__file__ = getattr(mod, "__file__", None)
    if mod__file__ is None:
        return False
    return os.path.exists(mod__file__)


def file_and_path_for_module(modulename: str) -> tuple[str | None, list[str]]:
    """Find the file and search path for `modulename`.

    Returns:
        filename: The filename of the module, or None.
        path: A list (possibly empty) of directories to find submodules in.

    """
    filename = None
    path = []
    try:
        spec = importlib.util.find_spec(modulename)
    except Exception:
        pass
    else:
        if spec is not None:
            filename = spec.origin
            path = list(spec.submodule_search_locations or ())
    return filename, path


def _add_sysconfig_paths(paths: set[str], path_names: list[str]) -> None:
    """Get paths from `sysconfig.get_paths`"""
    scheme_names = set(sysconfig.get_scheme_names())

    for scheme in scheme_names:
        config_paths = sysconfig.get_paths(scheme)
        for path_name in path_names:
            if path_name in config_paths:
                paths.add(config_paths[path_name])


def _add_stdlib_paths(paths: set[str]) -> None:
    """Add paths where the stdlib can be found to the set `paths`."""
    _add_sysconfig_paths(paths, ["stdlib", "platstdlib"])


def _add_third_party_paths(paths: set[str]) -> None:
    """Add locations for third-party packages to the set `paths`."""

    # These sysconfig locations are where third-party packages are installed.
    _add_sysconfig_paths(paths, ["platlib", "purelib", "scripts"])

    # Any importable directory that is a venv is also a third-party location.
    for d in sys.path:
        detail = _analyze_directory(d)
        if detail.exists and detail.venv is not None:
            paths.add(d)


def _add_coverage_paths(paths: set[str]) -> None:
    """Add paths where coverage.py code can be found to the set `paths`."""
    cover_path = canonical_path(__file__, directory=True)
    paths.add(cover_path)
    if env.TESTING:
        # Don't include our own test code.
        paths.add(os.path.join(cover_path, "tests"))


@dataclass
class DirectoryDetail:
    """Details about a directory."""

    exists: bool
    venv: str | None


def _analyze_directory(d: str) -> DirectoryDetail:
    """Analyze the directory `d` for existence and venv status."""
    detail = DirectoryDetail(exists=os.path.exists(d), venv=None)
    if detail.exists:
        while True:
            d = os.path.dirname(d)
            if d == os.path.dirname(d):
                break
            pyvenv = os.path.join(d, "pyvenv.cfg")
            if os.path.exists(pyvenv):
                detail.venv = d
                break
    return detail


def _dir_detail(d: str) -> str:
    """Get a string describing the directory `d` for debugging."""
    detail = _analyze_directory(d)
    if not detail.exists:
        describe = "does not exist"
    elif detail.venv is not None:
        describe = f"venv at {detail.venv}"
    else:
        describe = "not a venv"
    return f"{d!r} ({describe})"


class InOrOut:
    """Machinery for determining what files to measure."""

    def __init__(
        self,
        config: CoverageConfig,
        warn: TWarnFn,
        debug: TDebugCtl | None,
        include_namespace_packages: bool,
    ) -> None:
        self.warn = warn
        self.debug = debug
        self.include_namespace_packages = include_namespace_packages

        self.plugins: Plugins
        self.disp_class: type[TFileDisposition] = FileDisposition

        self.source_pkgs: list[str] = list(config.source_pkgs)
        self.source_dirs: list[str] = list(config.source_dirs)
        for src in config.source or []:
            if os.path.isdir(src):
                self.source_dirs.append(src)
            else:
                self.source_pkgs.append(src)

        # Canonicalize everything in `source_dirs`.
        # Also confirm that they actually are directories.
        for i, src in enumerate(self.source_dirs):
            if not os.path.isdir(src):
                raise ConfigError(f"Source dir is not a directory: {src!r}")
            self.source_dirs[i] = canonical_filename(src)

        self.source_pkgs_unmatched = self.source_pkgs[:]

        self.include = prep_patterns(config.run_include)
        self.omit = prep_patterns(config.run_omit)

        # The directories for files considered "installed with the interpreter".
        self.pylib_paths: set[str] = set()
        if not config.cover_pylib:
            _add_stdlib_paths(self.pylib_paths)

        # To avoid tracing the coverage.py code itself, we skip anything
        # located where we are.
        self.cover_paths: set[str] = set()
        _add_coverage_paths(self.cover_paths)

        # Find where third-party packages are installed.
        self.third_paths: set[str] = set()
        _add_third_party_paths(self.third_paths)

        # Generally useful information
        if self.debug:
            self._debug("sysconfig paths:")
            for scheme in sorted(sysconfig.get_scheme_names()):
                self._debug(f"    {scheme}:")
                for k, v in sysconfig.get_paths(scheme).items():
                    self._debug(f"        {k}: {_dir_detail(v)}")

        # Create the matchers we need for should_trace
        self.source_match = None
        self.source_pkgs_match = None
        self.pylib_match = None
        self.include_match = None
        self.omit_match = None

        if self.source_dirs or self.source_pkgs:
            if self.source_dirs:
                self.source_match = TreeMatcher(
                    self.source_dirs, "source", "Source directory", self._debug
                )
            if self.source_pkgs:
                self.source_pkgs_match = ModuleMatcher(
                    self.source_pkgs, "source_pkgs", "Source imports", self._debug
                )
        else:
            if self.pylib_paths:
                self.pylib_match = TreeMatcher(
                    self.pylib_paths, "pylib", "Python stdlib", self._debug
                )
        if self.include:
            self.include_match = GlobMatcher(self.include, "include", "Include", self._debug)
        if self.omit:
            self.omit_match = GlobMatcher(self.omit, "omit", "Omit", self._debug)

        self.coverage_match = TreeMatcher(
            self.cover_paths, "coverage", "Coverage code", self._debug
        )

        self.last_sys_path = list(sys.path)
        self.set_matchers_depending_on_syspath()

    def _debug(self, msg: str) -> None:
        """A more convenient way to write debug messages."""
        if self.debug:
            self.debug.write(msg)

    def set_matchers_depending_on_syspath(self) -> None:
        """Set up matchers that depend on sys.path.

        This is called at initialization time, and later if sys.path changes,
        which can happen when test runners like pytest manipulate sys.path.

        """
        self._debug("sys.path:" + "".join(f"\n    {_dir_detail(d)}" for d in sys.path))

        self.third_paths = set()
        _add_third_party_paths(self.third_paths)
        self.third_match = TreeMatcher(self.third_paths, "third", "Third-party lib", self._debug)

        # Check if the source we want to measure has been installed as a
        # third-party package.
        # Is the source inside a third-party area?
        self.source_in_third_paths = set()
        with sys_modules_saved():
            for pkg in self.source_pkgs:
                try:
                    modfile, path = file_and_path_for_module(pkg)
                    self._debug(f"Imported source package {pkg!r} as {modfile!r}")
                except CoverageException as exc:
                    self._debug(f"Couldn't import source package {pkg!r}: {exc}")
                    continue
                if modfile:
                    if self.third_match.match(modfile):
                        self._debug(
                            f"Source in third-party: source_pkg {pkg!r} at {modfile!r}",
                        )
                        self.source_in_third_paths.add(canonical_path(source_for_file(modfile)))
                else:
                    for pathdir in path:
                        if self.third_match.match(pathdir):
                            self._debug(
                                f"Source in third-party: {pkg!r} path directory at {pathdir!r}",
                            )
                            self.source_in_third_paths.add(pathdir)

        for src in self.source_dirs:
            if self.third_match.match(src):
                self._debug(f"Source in third-party: source directory {src!r}")
                self.source_in_third_paths.add(src)
        self.source_in_third_match = TreeMatcher(
            self.source_in_third_paths, "source_in_third", "Source in third-party", self._debug
        )

    def should_trace(self, filename: str, frame: FrameType | None = None) -> TFileDisposition:
        """Decide whether to trace execution in `filename`, with a reason.

        This function is called from the trace function.  As each new file name
        is encountered, this function determines whether it is traced or not.

        Returns a FileDisposition object.

        """
        if sys.path != self.last_sys_path:
            self.set_matchers_depending_on_syspath()
            self.last_sys_path = list(sys.path)

        original_filename = filename
        disp = disposition_init(self.disp_class, filename)

        def nope(disp: TFileDisposition, reason: str) -> TFileDisposition:
            """Simple helper to make it easy to return NO."""
            disp.trace = False
            disp.reason = reason
            return disp

        if original_filename.startswith("<"):
            return nope(disp, "original file name is not real")

        if frame is not None:
            # Compiled Python files have two file names: frame.f_code.co_filename is
            # the file name at the time the .pyc was compiled.  The second name is
            # __file__, which is where the .pyc was actually loaded from.  Since
            # .pyc files can be moved after compilation (for example, by being
            # installed), we look for __file__ in the frame and prefer it to the
            # co_filename value.
            dunder_file = frame.f_globals and frame.f_globals.get("__file__")
            if dunder_file:
                # Danger: __file__ can (rarely?) be of type Path.
                filename = source_for_file(str(dunder_file))
                if original_filename and not original_filename.startswith("<"):
                    orig = os.path.basename(original_filename)
                    if orig != os.path.basename(filename):
                        # Files shouldn't be renamed when moved. This happens when
                        # exec'ing code.  If it seems like something is wrong with
                        # the frame's file name, then just use the original.
                        filename = original_filename

        if not filename:
            # Empty string is pretty useless.
            return nope(disp, "empty string isn't a file name")

        if filename.startswith("memory:"):
            return nope(disp, "memory isn't traceable")

        if filename.startswith("<"):
            # Lots of non-file execution is represented with artificial
            # file names like "<string>", "<doctest readme.txt[0]>", or
            # "<exec_function>".  Don't ever trace these executions, since we
            # can't do anything with the data later anyway.
            return nope(disp, "file name is not real")

        canonical = canonical_filename(filename)
        disp.canonical_filename = canonical

        # Try the plugins, see if they have an opinion about the file.
        plugin = None
        for plugin in self.plugins.file_tracers:
            if not plugin._coverage_enabled:
                continue

            try:
                file_tracer = plugin.file_tracer(canonical)
                if file_tracer is not None:
                    file_tracer._coverage_plugin = plugin
                    disp.trace = True
                    disp.file_tracer = file_tracer
                    if file_tracer.has_dynamic_source_filename():
                        disp.has_dynamic_filename = True
                    else:
                        disp.source_filename = canonical_filename(
                            file_tracer.source_filename(),
                        )
                    break
            except Exception:
                plugin_name = plugin._coverage_plugin_name
                tb = traceback.format_exc()
                self.warn(f"Disabling plug-in {plugin_name!r} due to an exception:\n{tb}")
                plugin._coverage_enabled = False
                continue
        else:
            # No plugin wanted it: it's Python.
            disp.trace = True
            disp.source_filename = canonical

        if not disp.has_dynamic_filename:
            if not disp.source_filename:
                raise PluginError(
                    f"Plugin {plugin!r} didn't set source_filename for '{disp.original_filename}'",
                )
            reason = self.check_include_omit_etc(disp.source_filename, frame)
            if reason:
                nope(disp, reason)

        return disp

    def check_include_omit_etc(self, filename: str, frame: FrameType | None) -> str | None:
        """Check a file name against the include, omit, etc, rules.

        Returns a string or None.  String means, don't trace, and is the reason
        why.  None means no reason found to not trace.

        """
        modulename = name_for_module(filename, frame)

        # If the user specified source or include, then that's authoritative
        # about the outer bound of what to measure and we don't have to apply
        # any canned exclusions. If they didn't, then we have to exclude the
        # stdlib and coverage.py directories.
        if self.source_match or self.source_pkgs_match:
            extra = ""
            ok = False
            if self.source_pkgs_match:
                if isinstance(modulename, str) and self.source_pkgs_match.match(modulename):
                    ok = True
                    if modulename in self.source_pkgs_unmatched:
                        self.source_pkgs_unmatched.remove(modulename)
                else:
                    extra = f"module {modulename!r} "
            if not ok and self.source_match:
                if self.source_match.match(filename):
                    ok = True
            if not ok:
                return extra + "falls outside the --source spec"
            if self.third_match.match(filename) and not self.source_in_third_match.match(filename):
                return "inside --source, but is third-party"
        elif self.include_match:
            if not self.include_match.match(filename):
                return "falls outside the --include trees"
        else:
            # We exclude the coverage.py code itself, since a little of it
            # will be measured otherwise.
            if self.coverage_match.match(filename):
                return "is part of coverage.py"

            # Exclude anything in the third-party installation areas. Check this before
            # the stdlib, since site-packages is nested inside the stdlib area. If we
            # do it the other way around, third-party code will be labeled as stdlib
            # in the debug output.
            if self.third_match.match(filename):
                return "is a third-party module"

            # If we aren't supposed to trace installed code, then check if this
            # is in the Python standard library and skip it if so.
            if self.pylib_match and self.pylib_match.match(filename):
                return "is in the stdlib"

        # Check the file against the omit pattern.
        if self.omit_match and self.omit_match.match(filename):
            return "is inside an --omit pattern"

        # No point tracing a file we can't later write to SQLite.
        try:
            filename.encode("utf-8")
        except UnicodeEncodeError:
            return "non-encodable filename"

        # No reason found to skip this file.
        return None

    def warn_conflicting_settings(self) -> None:
        """Warn if there are settings that conflict."""
        if self.include:
            if self.source_dirs or self.source_pkgs:
                self.warn("--include is ignored because --source is set", slug="include-ignored")

    def warn_already_imported_files(self) -> None:
        """Warn if files have already been imported that we will be measuring."""
        if self.include or self.source_dirs or self.source_pkgs:
            warned = set()
            for mod in list(sys.modules.values()):
                filename = getattr(mod, "__file__", None)
                if filename is None:
                    continue
                if filename in warned:
                    continue

                if len(getattr(mod, "__path__", ())) > 1:
                    # A namespace package, which confuses this code, so ignore it.
                    continue

                disp = self.should_trace(filename)
                if disp.has_dynamic_filename:
                    # A plugin with dynamic filenames: the Python file
                    # shouldn't cause a warning, since it won't be the subject
                    # of tracing anyway.
                    continue
                if disp.trace:
                    msg = f"Already imported a file that will be measured: {filename}"
                    self.warn(msg, slug="already-imported")
                    warned.add(filename)
                elif self.debug and self.debug.should("trace"):
                    self.debug.write(
                        "Didn't trace already imported file {!r}: {}".format(
                            disp.original_filename,
                            disp.reason,
                        ),
                    )

    def warn_unimported_source(self) -> None:
        """Warn about source packages that were of interest, but never traced."""
        for pkg in self.source_pkgs_unmatched:
            self._warn_about_unmeasured_code(pkg)

    def _warn_about_unmeasured_code(self, pkg: str) -> None:
        """Warn about a package or module that we never traced.

        `pkg` is a string, the name of the package or module.

        """
        mod = sys.modules.get(pkg)
        if mod is None:
            self.warn(f"Module {pkg} was never imported.", slug="module-not-imported")
            return

        if module_is_namespace(mod):
            # A namespace package. It's OK for this not to have been traced,
            # since there is no code directly in it.
            return

        if not module_has_file(mod):
            self.warn(f"Module {pkg} has no Python source.", slug="module-not-python")
            return

        # The module was in sys.modules, and seems like a module with code, but
        # we never measured it. I guess that means it was imported before
        # coverage even started.
        msg = f"Module {pkg} was previously imported, but not measured"
        self.warn(msg, slug="module-not-measured")

    def find_possibly_unexecuted_files(self) -> Iterable[tuple[str, str | None]]:
        """Find files in the areas of interest that might be untraced.

        Yields pairs: file path, and responsible plug-in name.
        """
        for pkg in self.source_pkgs:
            if pkg not in sys.modules or not module_has_file(sys.modules[pkg]):
                continue
            pkg_file = source_for_file(cast(str, sys.modules[pkg].__file__))
            yield from self._find_executable_files(canonical_path(pkg_file))

        for src in self.source_dirs:
            yield from self._find_executable_files(src)

    def _find_plugin_files(self, src_dir: str) -> Iterable[tuple[str, str]]:
        """Get executable files from the plugins."""
        for plugin in self.plugins.file_tracers:
            for x_file in plugin.find_executable_files(src_dir):
                yield x_file, plugin._coverage_plugin_name

    def _find_executable_files(self, src_dir: str) -> Iterable[tuple[str, str | None]]:
        """Find executable files in `src_dir`.

        Search for files in `src_dir` that can be executed because they
        are probably importable. Don't include ones that have been omitted
        by the configuration.

        Yield the file path, and the plugin name that handles the file.

        """
        py_files = (
            (py_file, None)
            for py_file in find_python_files(src_dir, self.include_namespace_packages)
        )
        plugin_files = self._find_plugin_files(src_dir)

        for file_path, plugin_name in itertools.chain(py_files, plugin_files):
            file_path = canonical_filename(file_path)
            if self.omit_match and self.omit_match.match(file_path):
                # Turns out this file was omitted, so don't pull it back
                # in as un-executed.
                continue
            yield file_path, plugin_name

    def sys_info(self) -> Iterable[tuple[str, Any]]:
        """Our information for Coverage.sys_info.

        Returns a list of (key, value) pairs.
        """
        info = [
            ("coverage_paths", self.cover_paths),
            ("stdlib_paths", self.pylib_paths),
            ("third_party_paths", self.third_paths),
            ("source_in_third_party_paths", self.source_in_third_paths),
        ]

        matcher_names = [
            "source_match",
            "source_pkgs_match",
            "include_match",
            "omit_match",
            "coverage_match",
            "pylib_match",
            "third_match",
            "source_in_third_match",
        ]

        for matcher_name in matcher_names:
            matcher = getattr(self, matcher_name)
            if matcher:
                matcher_info = matcher.info()
            else:
                matcher_info = "-none-"
            info.append((matcher_name, matcher_info))

        return info
