"""Concrete Python interpreter information, also used as subprocess interrogation script (stdlib only)."""

from __future__ import annotations

import json
import logging
import os
import platform
import re
import struct
import sys
import sysconfig
import warnings
from collections import OrderedDict
from itertools import product
from string import digits
from typing import TYPE_CHECKING, ClassVar, Final, NamedTuple

if TYPE_CHECKING:
    import tkinter as tk
    from collections.abc import Generator, Mapping

    from ._cache import PyInfoCache
    from ._py_spec import PythonSpec


class VersionInfo(NamedTuple):
    major: int
    minor: int
    micro: int
    releaselevel: str
    serial: int


_LOGGER: Final[logging.Logger] = logging.getLogger(__name__)


def _get_path_extensions() -> list[str]:
    return list(OrderedDict.fromkeys(["", *os.environ.get("PATHEXT", "").lower().split(os.pathsep)]))


EXTENSIONS: Final[list[str]] = _get_path_extensions()
_32BIT_POINTER_SIZE: Final[int] = 4
_CONF_VAR_RE: Final[re.Pattern[str]] = re.compile(
    r"""
    \{ \w+  }   # sysconfig variable placeholder like {base}
    """,
    re.VERBOSE,
)


class PythonInfo:  # ruff:ignore[too-many-public-methods]
    """Contains information for a Python interpreter."""

    def __init__(self) -> None:
        self._init_identity()
        self._init_prefixes()
        self._init_schemes()
        self._init_sysconfig()

    def _init_identity(self) -> None:
        self.platform = sys.platform
        self.implementation = platform.python_implementation()
        if self.implementation == "GraalVM":
            self.implementation = "GraalPy"
        if self.implementation == "PyPy":
            self.pypy_version_info = tuple(sys.pypy_version_info)  # ty: ignore[unresolved-attribute] # pypy only

        self.version_info = VersionInfo(*sys.version_info)
        # same as stdlib platform.architecture to account for pointer size != max int
        self.architecture = 32 if struct.calcsize("P") == _32BIT_POINTER_SIZE else 64
        self.sysconfig_platform = sysconfig.get_platform()
        self.version_nodot = sysconfig.get_config_var("py_version_nodot")
        self.version = sys.version
        self.os = os.name
        self.free_threaded = sysconfig.get_config_var("Py_GIL_DISABLED") == 1
        self.debug_build = bool(sysconfig.get_config_var("Py_DEBUG"))

    def _init_prefixes(self) -> None:
        def abs_path(value: str | None) -> str | None:
            return None if value is None else os.path.abspath(value)

        self.prefix = abs_path(getattr(sys, "prefix", None))
        self.base_prefix = abs_path(getattr(sys, "base_prefix", None))
        self.real_prefix = abs_path(getattr(sys, "real_prefix", None))
        self.base_exec_prefix = abs_path(getattr(sys, "base_exec_prefix", None))
        self.exec_prefix = abs_path(getattr(sys, "exec_prefix", None))

        self.executable = abs_path(sys.executable)
        self.original_executable = abs_path(self.executable)
        self.system_executable = self._fast_get_system_executable()

        try:
            __import__("venv")
            has = True
        except ImportError:  # pragma: no cover # venv is always available in standard CPython
            has = False
        self.has_venv = has
        self.path = sys.path
        self.file_system_encoding = sys.getfilesystemencoding()
        self.stdout_encoding = getattr(sys.stdout, "encoding", None)

    def _init_schemes(self) -> None:
        scheme_names = sysconfig.get_scheme_names()

        if "venv" in scheme_names:  # pragma: >=3.11 cover
            self.sysconfig_scheme = "venv"
            self.sysconfig_paths = {
                i: sysconfig.get_path(i, expand=False, scheme=self.sysconfig_scheme) for i in sysconfig.get_path_names()
            }
            self.distutils_install = {}
        # debian / ubuntu python 3.10 without `python3-distutils` will report mangled `local/bin` / etc. names
        elif sys.version_info[:2] == (3, 10) and "deb_system" in scheme_names:  # pragma: no cover # Debian/Ubuntu 3.10
            self.sysconfig_scheme = "posix_prefix"
            self.sysconfig_paths = {
                i: sysconfig.get_path(i, expand=False, scheme=self.sysconfig_scheme) for i in sysconfig.get_path_names()
            }
            self.distutils_install = {}
        else:  # pragma: no cover # "venv" scheme always present on Python 3.12+
            self.sysconfig_scheme = None
            self.sysconfig_paths = {i: sysconfig.get_path(i, expand=False) for i in sysconfig.get_path_names()}
            self.distutils_install = self._distutils_install().copy()

    def _init_sysconfig(self) -> None:
        makefile = getattr(sysconfig, "get_makefile_filename", getattr(sysconfig, "_get_makefile_filename", None))
        self.sysconfig = {
            k: v
            for k, v in [
                ("makefile_filename", makefile() if makefile is not None else None),
            ]
            if k is not None
        }

        config_var_keys = set()
        for element in self.sysconfig_paths.values():
            config_var_keys.update(k[1:-1] for k in _CONF_VAR_RE.findall(element))
        config_var_keys.add("PYTHONFRAMEWORK")
        config_var_keys.update(("Py_ENABLE_SHARED", "INSTSONAME", "LIBDIR"))

        self.sysconfig_vars = {i: sysconfig.get_config_var(i or "") for i in config_var_keys}

        if "TCL_LIBRARY" in os.environ:
            self.tcl_lib, self.tk_lib = self._get_tcl_tk_libs()
        else:
            self.tcl_lib, self.tk_lib = None, None

        confs = {
            k: (self.system_prefix if isinstance(v, str) and v.startswith(self.prefix) else v)
            for k, v in self.sysconfig_vars.items()
        }
        self.system_stdlib = self.sysconfig_path("stdlib", confs)
        self.system_stdlib_platform = self.sysconfig_path("platstdlib", confs)
        self.max_size = getattr(sys, "maxsize", getattr(sys, "maxint", None))
        self._creators = None  # virtualenv-specific, set via monkey-patch

    @staticmethod
    def _get_tcl_tk_libs() -> tuple[
        str | None,
        str | None,
    ]:  # pragma: no cover # tkinter availability varies; tested indirectly via __init__
        """Detect the tcl and tk libraries using tkinter."""
        tcl_lib, tk_lib = None, None
        try:
            import tkinter as tk  # ruff:ignore[import-outside-top-level]
        except ImportError:
            pass
        else:
            try:
                tcl = tk.Tcl()
                tcl_lib = tcl.eval("info library")
                tk_lib = PythonInfo._resolve_tk_lib(tcl, tcl_lib)
            except tk.TclError:
                pass

        return tcl_lib, tk_lib

    @staticmethod
    def _query_tk_library(tcl: tk.Tk) -> str | None:  # pragma: no cover
        """Try to get the TK library path directly from Tcl."""
        import tkinter as tk  # ruff:ignore[import-outside-top-level]

        try:
            if (tk_lib := tcl.eval("set tk_library")) and os.path.isdir(tk_lib):
                return tk_lib
        except tk.TclError:
            pass
        return None

    @staticmethod
    def _resolve_tk_lib(tcl: tk.Tk, tcl_lib: str) -> str | None:  # pragma: no cover
        """Resolve the TK library path by direct query or path construction."""
        if (tk_lib := PythonInfo._query_tk_library(tcl)) is not None:
            return tk_lib
        tk_version = tcl.eval("package require Tk")
        tcl_parent = os.path.dirname(tcl_lib)
        for version in (tk_version, ".".join(tk_version.split(".")[:2]), tk_version.split(".")[0]):
            tk_lib_path = os.path.join(tcl_parent, f"tk{version}")
            if os.path.isdir(tk_lib_path) and os.path.exists(os.path.join(tk_lib_path, "tk.tcl")):
                return tk_lib_path
        return None

    def _fast_get_system_executable(self) -> str | None:
        """Try to get the system executable by just looking at properties."""
        # if we're not in a virtual environment, this is already a system python, so return the original executable
        # note we must choose the original and not the pure executable as shim scripts might throw us off
        if not (self.real_prefix or (self.base_prefix is not None and self.base_prefix != self.prefix)):
            return self._resolve_executable_symlink(self.original_executable)

        # if this is NOT a virtual environment, can't determine easily, bail out
        if self.real_prefix is not None:
            return None

        base_executable = getattr(sys, "_base_executable", None)  # some platforms may set this to help us
        if base_executable is None:  # use the saved system executable if present
            return None

        # we know we're in a virtual environment, can not be us
        if sys.executable == base_executable:
            return None

        # We're not in a venv and base_executable exists; use it directly
        if os.path.exists(base_executable):  # pragma: >=3.11 cover
            return self._resolve_executable_symlink(base_executable)

        # Try fallback for POSIX virtual environments
        return self._try_posix_fallback_executable(base_executable)  # pragma: >=3.11 cover

    def _resolve_executable_symlink(self, path: str, *, framework: bool | None = None) -> str:
        """
        Resolve symlinks of the executable itself, but never of its parent directories.

        Mirrors CPython's ``getpath.realpath`` (and ``venv`` in python/cpython#115237): an executable-only symlink
        resolves to the real interpreter so its home can be located, while a fully symlinked interpreter tree is
        kept as-is. Like ``getpath``, resolution stops as soon as the stdlib landmark is reachable from the current
        directory - an alias such as Debian's ``/usr/bin/python3`` is a usable home and stays untouched.
        """
        result = os.path.abspath(path)
        if self.os != "posix":  # CPython only does this where HAVE_READLINK
            return result
        if framework is None:
            framework = bool(sysconfig.get_config_var("PYTHONFRAMEWORK"))
        if framework:  # macOS framework builds self-locate via dyld from the real binary; e.g. for Homebrew
            return result  # resolving would pin the versioned Cellar path into the recorded home
        real_path = os.path.realpath(result)
        if not os.path.exists(real_path):  # symlink loop or broken symlink
            return result
        while os.path.islink(result):
            if self._stdlib_landmark_exists(os.path.dirname(result)):
                return result
            link = os.readlink(result)
            candidate = link if os.path.isabs(link) else os.path.normpath(os.path.join(os.path.dirname(result), link))
            # normpath through a symlinked directory may point at a different file - stop resolving there
            if not (os.path.exists(candidate) and os.path.samefile(real_path, candidate)):
                return result
            result = candidate
        return result

    @staticmethod
    def _stdlib_landmark_exists(dir_path: str) -> bool:
        lib_name = os.path.basename(os.path.dirname(os.__file__))
        return any(
            os.path.exists(os.path.join(dir_path, os.pardir, lib, lib_name, "os.py")) for lib in ("lib", "lib64")
        )

    def _try_posix_fallback_executable(self, base_executable: str) -> str | None:
        """Find a versioned Python binary as fallback for POSIX virtual environments."""
        major, minor = self.version_info.major, self.version_info.minor
        if self.os != "posix" or (major, minor) < (3, 11):
            return None

        # search relative to the directory of sys._base_executable
        base_dir = os.path.dirname(base_executable)
        candidates = [f"python{major}", f"python{major}.{minor}"]
        if self.implementation == "PyPy":
            candidates.extend(["pypy", "pypy3", f"pypy{major}", f"pypy{major}.{minor}"])

        for candidate in candidates:
            full_path = os.path.join(base_dir, candidate)
            if os.path.exists(full_path):
                return full_path

        return None  # in this case we just can't tell easily without poking around FS and calling them, bail

    def install_path(self, key: str) -> str:
        """
        Return the relative installation path for a given installation scheme *key*.

        :param key: sysconfig installation scheme key (e.g. ``"scripts"``, ``"purelib"``).
        """
        result = self.distutils_install.get(key)
        if result is None:  # pragma: >=3.11 cover # distutils is empty when "venv" scheme is available
            # set prefixes to empty => result is relative from cwd
            prefixes = self.prefix, self.exec_prefix, self.base_prefix, self.base_exec_prefix
            config_var = {k: "" if v in prefixes else v for k, v in self.sysconfig_vars.items()}
            result = self.sysconfig_path(key, config_var=config_var).lstrip(os.sep)
        return result

    @staticmethod
    def _distutils_install() -> dict[str, str]:
        # use distutils primarily because that's what pip does
        # https://github.com/pypa/pip/blob/main/src/pip/_internal/locations.py#L95
        # note here we don't import Distribution directly to allow setuptools to patch it
        with warnings.catch_warnings():  # disable warning for PEP-632
            warnings.simplefilter("ignore")
            try:
                # ruff:ignore[import-outside-top-level]
                from distutils import dist  # ty: ignore[unresolved-import]

                # ruff:ignore[import-outside-top-level]
                from distutils.command.install import SCHEME_KEYS  # ty: ignore[unresolved-import]
            except ImportError:  # pragma: no cover # if removed or not installed ignore
                return {}

        distribution = dist.Distribution({
            "script_args": "--no-user-cfg",
        })  # conf files not parsed so they do not hijack paths
        if hasattr(sys, "_framework"):  # pragma: no cover # macOS framework builds only
            sys._framework = None  # ruff:ignore[private-member-access]  # disable macOS static paths for framework

        with warnings.catch_warnings():  # disable warning for PEP-632
            warnings.simplefilter("ignore")
            install = distribution.get_command_obj("install", create=True)

        install.prefix = os.sep  # paths generated are relative to prefix that contains the path sep
        install.finalize_options()
        return {key: (getattr(install, f"install_{key}")[1:]).lstrip(os.sep) for key in SCHEME_KEYS}

    @property
    def version_str(self) -> str:
        """The full version as ``major.minor.micro`` string (e.g. ``3.13.2``)."""
        return ".".join(str(i) for i in self.version_info[0:3])

    @property
    def version_release_str(self) -> str:
        """The release version as ``major.minor`` string (e.g. ``3.13``)."""
        return ".".join(str(i) for i in self.version_info[0:2])

    @property
    def python_name(self) -> str:
        """The python executable name as ``pythonX.Y`` (e.g. ``python3.13``)."""
        version_info = self.version_info
        return f"python{version_info.major}.{version_info.minor}"

    @property
    def is_old_virtualenv(self) -> bool:
        """``True`` if this interpreter runs inside an old-style virtualenv (has ``real_prefix``)."""
        return self.real_prefix is not None

    @property
    def is_venv(self) -> bool:
        """``True`` if this interpreter runs inside a PEP 405 venv (has ``base_prefix``)."""
        return self.base_prefix is not None

    def sysconfig_path(self, key: str, config_var: dict[str, str] | None = None, sep: str = os.sep) -> str:
        """
        Return the sysconfig install path for a scheme *key*, optionally substituting config variables.

        :param key: sysconfig path key (e.g. ``"purelib"``, ``"include"``).
        :param config_var: replacement mapping for sysconfig variables; when ``None`` uses the interpreter's own values.
        :param sep: path separator to use in the result.
        """
        pattern = self.sysconfig_paths.get(key)
        if pattern is None:
            return ""
        if config_var is None:
            config_var = self.sysconfig_vars
        else:
            base = self.sysconfig_vars.copy()
            base.update(config_var)
            config_var = base
        return pattern.format(**config_var).replace("/", sep)

    @property
    def system_include(self) -> str:
        """The path to the system include directory for C headers."""
        path = self.sysconfig_path(
            "include",
            {
                k: (self.system_prefix if isinstance(v, str) and v.startswith(self.prefix) else v)
                for k, v in self.sysconfig_vars.items()
            },
        )
        if not os.path.exists(path):  # pragma: no cover # broken packaging fallback
            fallback = os.path.join(self.prefix, os.path.dirname(self.install_path("headers")))
            if os.path.exists(fallback):
                path = fallback
        return path

    @property
    def system_prefix(self) -> str:
        """The prefix of the system Python this interpreter is based on."""
        return self.real_prefix or self.base_prefix or self.prefix

    @property
    def system_exec_prefix(self) -> str:
        """The exec prefix of the system Python this interpreter is based on."""
        return self.real_prefix or self.base_exec_prefix or self.exec_prefix

    def __repr__(self) -> str:
        return "{}({!r})".format(
            self.__class__.__name__,
            {k: v for k, v in self.__dict__.items() if not k.startswith("_")},
        )

    def __str__(self) -> str:
        return "{}({})".format(
            self.__class__.__name__,
            ", ".join(
                f"{k}={v}"
                for k, v in (
                    ("spec", self.spec),
                    (
                        "system"
                        if self.system_executable is not None and self.system_executable != self.executable
                        else None,
                        self.system_executable,
                    ),
                    (
                        "original"
                        if self.original_executable not in {self.system_executable, self.executable}
                        else None,
                        self.original_executable,
                    ),
                    ("exe", self.executable),
                    ("platform", self.platform),
                    ("version", repr(self.version)),
                    ("encoding_fs_io", f"{self.file_system_encoding}-{self.stdout_encoding}"),
                )
                if k is not None
            ),
        )

    @property
    def machine(self) -> str:
        """The instruction set architecture (ISA) derived from :func:`sysconfig.get_platform`."""
        plat = self.sysconfig_platform
        if plat is None:
            return "unknown"
        if plat == "win32":
            return "x86"
        isa = plat.rsplit("-", 1)[-1]
        if isa == "universal2":
            isa = platform.machine().lower()
        return normalize_isa(isa)

    @property
    def spec(self) -> str:
        """A specification string identifying this interpreter (e.g. ``CPython3.13.2-64-arm64``)."""
        return "{}{}{}{}-{}-{}".format(
            self.implementation,
            ".".join(str(i) for i in self.version_info),
            "t" if self.free_threaded else "",
            "d" if self.debug_build else "",
            self.architecture,
            self.machine,
        )

    @classmethod
    def clear_cache(cls, cache: PyInfoCache) -> None:
        """
        Clear all cached interpreter information from *cache*.

        :param cache: the cache store to clear.
        """
        from ._cached_py_info import clear  # ruff:ignore[import-outside-top-level]

        clear(cache)
        cls._cache_exe_discovery.clear()

    def satisfies(self, spec: PythonSpec, *, impl_must_match: bool) -> bool:  # ruff:ignore[too-many-return-statements]
        """
        Check if a given specification can be satisfied by this python interpreter instance.

        :param spec: the specification to check against.
        :param impl_must_match: when ``True``, the implementation name must match exactly.
        """
        if spec.path and not self._satisfies_path(spec):
            return False
        if impl_must_match and not self._satisfies_implementation(spec):
            return False
        if spec.architecture is not None and spec.architecture != self.architecture:
            return False
        if spec.machine is not None and spec.machine != self.machine:
            return False
        if spec.free_threaded is not None and spec.free_threaded != self.free_threaded:
            return False
        if spec.debug is not None and spec.debug != self.debug_build:
            return False
        if spec.version_specifier is not None and not self._satisfies_version_specifier(spec):
            return False
        return all(
            req is None or our is None or our == req
            for our, req in zip(self.version_info[0:3], (spec.major, spec.minor, spec.micro))
        )

    def _satisfies_path(self, spec: PythonSpec) -> bool:
        if self.executable == os.path.abspath(spec.path):
            return True
        if spec.is_abs:
            return True
        basename = os.path.basename(self.original_executable)
        spec_path = spec.path
        if sys.platform == "win32":
            basename, suffix = os.path.splitext(basename)
            spec_path = spec_path[: -len(suffix)] if suffix and spec_path.endswith(suffix) else spec_path
        return basename == spec_path

    def _satisfies_implementation(self, spec: PythonSpec) -> bool:
        return spec.implementation is None or spec.implementation.lower() == self.implementation.lower()

    def _satisfies_version_specifier(self, spec: PythonSpec) -> bool:
        if spec.version_specifier is None:  # pragma: no cover
            return True
        version_info = self.version_info
        for specifier in spec.version_specifier:
            assert specifier.version is not None  # ruff:ignore[assert]
            numeric_version = specifier.version_str
            for prefix in ("rc", "b", "a"):
                if prefix in numeric_version:
                    numeric_version = numeric_version.split(prefix)[0]
                    break
            precision = numeric_version.count(".") + 1
            release = ".".join(str(c) for c in [version_info.major, version_info.minor, version_info.micro][:precision])
            if (
                version_info.releaselevel != "final"
                and (precision == 3 or specifier.version.pre_type is not None)  # ruff:ignore[magic-value-comparison]
                and (suffix := {"alpha": "a", "beta": "b", "candidate": "rc"}.get(version_info.releaselevel))
            ):
                release = f"{release}{suffix}{version_info.serial}"
            if not specifier.contains(release):
                return False
        return True

    _current_system = None
    _current = None

    @classmethod
    def current(cls, cache: PyInfoCache | None = None) -> PythonInfo:
        """
        Locate the current host interpreter information.

        :param cache: interpreter metadata cache; when ``None`` results are not cached.
        """
        if cls._current is None:
            result = cls.from_exe(sys.executable, cache, raise_on_error=True, resolve_to_host=False)
            if result is None:
                msg = "failed to query current Python interpreter"
                raise RuntimeError(msg)
            cls._current = result
        return cls._current

    @classmethod
    def current_system(cls, cache: PyInfoCache | None = None) -> PythonInfo:
        """
        Locate the current system interpreter information, resolving through any virtualenv layers.

        :param cache: interpreter metadata cache; when ``None`` results are not cached.
        """
        if cls._current_system is None:
            result = cls.from_exe(sys.executable, cache, raise_on_error=True, resolve_to_host=True)
            if result is None:
                msg = "failed to query current system Python interpreter"
                raise RuntimeError(msg)
            cls._current_system = result
        return cls._current_system

    def to_json(self) -> str:
        """Serialize this interpreter information to a JSON string."""
        return json.dumps(self.to_dict(), indent=2)

    def to_dict(self) -> dict[str, object]:
        """Convert this interpreter information to a plain dictionary."""
        data = {var: (getattr(self, var) if var != "_creators" else None) for var in vars(self)}
        version_info = data["version_info"]
        data["version_info"] = version_info._asdict() if hasattr(version_info, "_asdict") else version_info
        return data

    @classmethod
    def from_exe(  # ruff:ignore[too-many-arguments]
        cls,
        exe: str,
        cache: PyInfoCache | None = None,
        *,
        raise_on_error: bool = True,
        ignore_cache: bool = False,
        resolve_to_host: bool = True,
        env: Mapping[str, str] | None = None,
    ) -> PythonInfo | None:
        """
        Get the python information for a given executable path.

        :param exe: path to the Python executable.
        :param cache: interpreter metadata cache; when ``None`` results are not cached.
        :param raise_on_error: raise on failure instead of returning ``None``.
        :param ignore_cache: bypass the cache and re-query the interpreter.
        :param resolve_to_host: resolve through virtualenv layers to the system interpreter.
        :param env: environment mapping; defaults to :data:`os.environ`.
        """
        from ._cached_py_info import from_exe  # ruff:ignore[import-outside-top-level]

        env = os.environ if env is None else env
        proposed = from_exe(cls, cache, exe, env=env, raise_on_error=raise_on_error, ignore_cache=ignore_cache)

        if isinstance(proposed, PythonInfo) and resolve_to_host:
            try:
                proposed = proposed.resolve_to_system(cache, proposed)
            except Exception as exception:
                if raise_on_error:
                    raise
                _LOGGER.info("ignore %s due cannot resolve system due to %r", proposed.original_executable, exception)
                proposed = None
        return proposed

    @classmethod
    def from_json(cls, payload: str) -> PythonInfo:
        """
        Deserialize interpreter information from a JSON string.

        :param payload: JSON produced by :meth:`to_json`.
        """
        raw = json.loads(payload)
        return cls.from_dict(raw.copy())

    @classmethod
    def from_dict(cls, data: dict[str, object]) -> PythonInfo:
        """
        Reconstruct a :class:`PythonInfo` from a plain dictionary.

        :param data: dictionary produced by :meth:`to_dict`.
        """
        data["version_info"] = VersionInfo(**data["version_info"])  # restore this to a named tuple structure
        result = cls()
        result.__dict__ = data.copy()
        return result

    @classmethod
    def resolve_to_system(cls, cache: PyInfoCache | None, target: PythonInfo) -> PythonInfo:
        """
        Walk virtualenv/venv prefix chains to find the underlying system interpreter.

        :param cache: interpreter metadata cache; when ``None`` results are not cached.
        :param target: the interpreter to resolve.
        """
        start_executable = target.executable
        prefixes = OrderedDict()
        while target.system_executable is None:
            prefix = target.real_prefix or target.base_prefix or target.prefix
            if prefix in prefixes:
                if len(prefixes) == 1:
                    _LOGGER.info("%r links back to itself via prefixes", target)
                    target.system_executable = target.executable
                    break
                for at, (p, t) in enumerate(prefixes.items(), start=1):
                    _LOGGER.error("%d: prefix=%s, info=%r", at, p, t)
                _LOGGER.error("%d: prefix=%s, info=%r", len(prefixes) + 1, prefix, target)
                msg = "prefixes are causing a circle {}".format("|".join(prefixes.keys()))
                raise RuntimeError(msg)
            prefixes[prefix] = target
            target = target.discover_exe(cache, prefix=prefix, exact=False)
        if target.executable != target.system_executable:
            resolved = cls.from_exe(target.system_executable, cache)
            if resolved is not None:
                target = resolved
        target.executable = start_executable
        return target

    _cache_exe_discovery: ClassVar[dict[tuple[str, bool], PythonInfo]] = {}

    def discover_exe(
        self,
        cache: PyInfoCache,
        prefix: str,
        *,
        exact: bool = True,
        env: Mapping[str, str] | None = None,
    ) -> PythonInfo:
        """
        Discover a matching Python executable under a given *prefix* directory.

        :param cache: interpreter metadata cache.
        :param prefix: directory prefix to search under.
        :param exact: when ``True``, require an exact version match.
        :param env: environment mapping; defaults to :data:`os.environ`.
        """
        key = prefix, exact
        if key in self._cache_exe_discovery and prefix:
            _LOGGER.debug("discover exe from cache %s - exact %s: %r", prefix, exact, self._cache_exe_discovery[key])
            return self._cache_exe_discovery[key]
        _LOGGER.debug("discover exe for %s in %s", self, prefix)
        possible_names = self._find_possible_exe_names()
        possible_folders = self._find_possible_folders(prefix)
        discovered = []
        env = os.environ if env is None else env
        for folder in possible_folders:
            for name in possible_names:
                info = self._check_exe(cache, folder, name, discovered, env, exact=exact)
                if info is not None:
                    self._cache_exe_discovery[key] = info
                    return info
        if exact is False and discovered:
            info = self._select_most_likely(discovered, self)
            folders = os.pathsep.join(possible_folders)
            self._cache_exe_discovery[key] = info
            _LOGGER.debug("no exact match found, chosen most similar of %s within base folders %s", info, folders)
            return info
        msg = "failed to detect {} in {}".format("|".join(possible_names), os.pathsep.join(possible_folders))
        raise RuntimeError(msg)

    def _check_exe(  # ruff:ignore[too-many-arguments]
        self,
        cache: PyInfoCache | None,
        folder: str,
        name: str,
        discovered: list[PythonInfo],
        env: Mapping[str, str],
        *,
        exact: bool,
    ) -> PythonInfo | None:
        exe_path = os.path.join(folder, name)
        if not os.path.exists(exe_path):
            return None
        info = self.from_exe(exe_path, cache, resolve_to_host=False, raise_on_error=False, env=env)
        if info is None:  # ignore if for some reason we can't query
            return None
        for item in ["implementation", "architecture", "machine", "version_info", "free_threaded", "debug_build"]:
            found = getattr(info, item)
            searched = getattr(self, item)
            if found != searched:
                if item == "version_info":
                    found, searched = ".".join(str(i) for i in found), ".".join(str(i) for i in searched)
                executable = info.executable
                _LOGGER.debug("refused interpreter %s because %s differs %s != %s", executable, item, found, searched)
                if exact is False:
                    discovered.append(info)
                break
        else:
            return info
        return None

    @staticmethod
    def _select_most_likely(discovered: list[PythonInfo], target: PythonInfo) -> PythonInfo:
        def sort_by(info: PythonInfo) -> int:
            # we need to setup some priority of traits, this is as follows:
            # implementation, major, minor, architecture, machine, micro, tag, serial
            matches = [
                info.implementation == target.implementation,
                info.version_info.major == target.version_info.major,
                info.version_info.minor == target.version_info.minor,
                info.architecture == target.architecture,
                info.machine == target.machine,
                info.version_info.micro == target.version_info.micro,
                info.version_info.releaselevel == target.version_info.releaselevel,
                info.version_info.serial == target.version_info.serial,
            ]
            return sum((1 << pos if match else 0) for pos, match in enumerate(reversed(matches)))

        sorted_discovered = sorted(discovered, key=sort_by, reverse=True)  # sort by priority in decreasing order
        return sorted_discovered[0]

    def _find_possible_folders(self, inside_folder: str) -> list[str]:
        candidate_folder = OrderedDict()
        executables = OrderedDict()
        executables[os.path.realpath(self.executable)] = None
        executables[self.executable] = None
        executables[os.path.realpath(self.original_executable)] = None
        executables[self.original_executable] = None
        for exe in executables:
            base = os.path.dirname(exe)
            if base.startswith(self.prefix):
                relative = base[len(self.prefix) :]
                candidate_folder[f"{inside_folder}{relative}"] = None

        # or at root level
        candidate_folder[inside_folder] = None
        return [i for i in candidate_folder if os.path.exists(i)]

    def _find_possible_exe_names(self) -> list[str]:
        name_candidate = OrderedDict()
        mods = ["", "t"] if self.free_threaded else [""]
        # _d (Windows python_d.exe), -dbg (Debian's own name), d (the abiflags form used by
        # upstream/Fedora and also shipped by Debian, so pythonX.Yd / pythonX.Ytd resolve too)
        debug_suffixes = ["_d", "-dbg", "d", ""] if self.debug_build else [""]
        archs = [f"-{self.architecture}", ""]
        for name in self._possible_base():
            for at in (3, 2, 1, 0):
                version = ".".join(str(i) for i in self.version_info[:at])
                for mod, debug, arch, ext in product(mods, debug_suffixes, archs, EXTENSIONS):
                    name_candidate[f"{name}{version}{mod}{debug}{arch}{ext}"] = None
        return list(name_candidate.keys())

    def _possible_base(self) -> Generator[str, None, None]:
        possible_base = OrderedDict()
        basename = os.path.splitext(os.path.basename(self.executable))[0].rstrip(digits)
        possible_base[basename] = None
        possible_base[self.implementation] = None
        # python is always the final option as in practice is used by multiple implementation as exe name
        if "python" in possible_base:
            del possible_base["python"]
        possible_base["python"] = None
        for base in possible_base:
            lower = base.lower()
            yield lower
            from ._compat import fs_is_case_sensitive  # ruff:ignore[import-outside-top-level]

            if fs_is_case_sensitive():  # pragma: no branch
                if base != lower:
                    yield base
                upper = base.upper()
                if upper != base:
                    yield upper


KNOWN_ARCHITECTURES: frozenset[str] = frozenset({
    "arm64",
    "loongarch64",
    "ppc",
    "ppc64",
    "ppc64le",
    "riscv64",
    "s390x",
    "sparc64",
    "x86",
    "x86_64",
})
"""Known CPU architecture (ISA) values after normalization.

.. deprecated::
    Use :func:`normalize_isa` instead, which handles both known and unknown architectures.
"""


def normalize_isa(isa: str) -> str:
    """
    Normalize an ISA (instruction set architecture) string to a canonical form.

    Known aliases are mapped (e.g. ``amd64`` → ``x86_64``, ``aarch64`` → ``arm64``).
    Unrecognized values are lowercased and returned as-is.
    """
    low = isa.lower()
    return {
        "amd64": "x86_64",
        "aarch64": "arm64",
        "i386": "x86",
        "i486": "x86",
        "i586": "x86",
        "i686": "x86",
        "powerpc": "ppc",
        "powerpc64": "ppc64",
        "powerpc64le": "ppc64le",
        "sparcv9": "sparc64",
    }.get(low, low)


def _main() -> None:  # pragma: no cover
    argv = sys.argv[1:]

    if len(argv) >= 1:
        start_cookie = argv[0]
        argv = argv[1:]
    else:
        start_cookie = ""

    if len(argv) >= 1:
        end_cookie = argv[0]
        argv = argv[1:]
    else:
        end_cookie = ""

    sys.argv = sys.argv[:1] + argv

    result = PythonInfo().to_json()
    sys.stdout.write("".join((start_cookie[::-1], result, end_cookie[::-1])))
    sys.stdout.flush()


if __name__ == "__main__":
    _main()
