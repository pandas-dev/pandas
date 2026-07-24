"""Utilities for installing Javascript extensions for the notebook."""

# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.

from __future__ import annotations

import importlib
import json
import os
import platform
import shutil
import subprocess
import sys
from pathlib import Path
from typing import TYPE_CHECKING, Any, cast

if TYPE_CHECKING:
    import logging

from importlib.metadata import PackageNotFoundError, version

from jupyter_core.paths import ENV_JUPYTER_PATH, SYSTEM_JUPYTER_PATH, jupyter_data_dir
from jupyter_core.utils import ensure_dir_exists

from .federated_extensions_requirements import get_federated_extensions

if sys.version_info >= (3, 11):
    from tomllib import load
else:
    from tomli import load

from .commands import _test_overlap
from .core_path import get_core_meta
from .jlpm import _which_node_js
from .jupyterlab_semver import clean, satisfies

DEPRECATED_ARGUMENT = object()

HERE = str(Path(__file__).resolve().parent)


class ArgumentConflict(ValueError):  # noqa: N818
    """Raised when conflicting arguments are passed to a labextension function."""


# ------------------------------------------------------------------------------
# Public API
# ------------------------------------------------------------------------------


def develop_labextension(  # noqa: PLR0913, C901, PLR0912
    path: str | os.PathLike[str],
    symlink: bool = True,
    overwrite: bool = False,
    user: bool = False,
    labextensions_dir: str | None = None,
    destination: str | None = None,
    logger: logging.Logger | None = None,
    sys_prefix: bool = False,
) -> str:
    """Install a prebuilt extension for JupyterLab.

    Stages files and/or directories into the labextensions directory.
    By default, this compares modification time, and only stages files that need updating.
    If `overwrite` is specified, matching files are purged before proceeding.

    Parameters
    ----------
    path : path to file, directory, zip or tarball archive, or URL to install
        By default, the file will be installed with its base name, so '/path/to/foo'
        will install to 'labextensions/foo'. See the destination argument below to change this.
        Archives (zip or tarballs) will be extracted into the labextensions directory.
    user : bool [default: False]
        Whether to install to the user's labextensions directory.
        Otherwise do a system-wide install (e.g. /usr/local/share/jupyter/labextensions).
    overwrite : bool [default: False]
        If True, always install the files, regardless of what may already be installed.
    symlink : bool [default: True]
        If True, create a symlink in labextensions, rather than copying files.
        Windows support for symlinks requires a permission bit which only admin users
        have by default, so don't rely on it.
    labextensions_dir : str [optional]
        Specify absolute path of labextensions directory explicitly.
    destination : str [optional]
        name the labextension is installed to.  For example, if destination is 'foo', then
        the source file will be installed to 'labextensions/foo', regardless of the source name.
    logger : Jupyter logger [optional]
        Logger instance to use
    sys_prefix : bool [default: False]
        If True, install to sys.prefix (i.e. the virtualenv/conda env prefix).

    """
    # the actual path to which we eventually installed
    full_dest = None

    labext = _get_labextension_dir(
        user=user,
        sys_prefix=sys_prefix,
        labextensions_dir=labextensions_dir,
    )
    # make sure labextensions dir exists
    ensure_dir_exists(labext)

    if isinstance(path, list | tuple):
        msg = "path must be a string pointing to a single extension to install; call this function multiple times to install multiple extensions"  # noqa: E501
        raise TypeError(msg)

    if not destination:
        destination = Path(path).name

    full_dest = str(Path(labext) / destination)
    if overwrite and os.path.lexists(full_dest):
        if logger:
            logger.info("Removing: %s", full_dest)
        if Path(full_dest).is_dir() and not Path(full_dest).is_symlink():
            shutil.rmtree(full_dest)
        else:
            Path(full_dest).unlink()

    # Make sure the parent directory exists
    Path(full_dest).parent.mkdir(parents=True, exist_ok=True)

    if symlink:
        path = str(Path(path).resolve())
        if not Path(full_dest).exists():
            if logger:
                logger.info("Symlinking: %s -> %s", full_dest, path)
            try:
                Path(full_dest).symlink_to(path)
            except OSError as e:
                if platform.platform().startswith("Windows"):
                    msg = (
                        "Symlinks can be activated on Windows 10 for Python version 3.8 or higher"
                        " by activating the 'Developer Mode'. That may not be allowed by your administrators.\n"  # noqa: E501
                        "See https://docs.microsoft.com/en-us/windows/apps/get-started/enable-your-device-for-development"
                    )
                    raise OSError(msg) from e
                raise

        elif not Path(full_dest).is_symlink():
            msg = f"{full_dest} exists and is not a symlink"
            raise ValueError(msg)

    elif Path(path).is_dir():
        src_root = Path(path).resolve()
        for parent_str, _, files in os.walk(src_root):
            parent_path = Path(parent_str)
            dest_dir = Path(full_dest) / parent_path.relative_to(src_root)
            if not dest_dir.exists():
                if logger:
                    logger.info("Making directory: %s", dest_dir)
                dest_dir.mkdir(parents=True)
            for file_name in files:
                src = str(parent_path / file_name)
                dest_file = str(dest_dir / file_name)
                _maybe_copy(src, dest_file, logger=logger)
    else:
        src = str(path)
        _maybe_copy(src, full_dest, logger=logger)

    return full_dest


def develop_labextension_py(  # noqa: PLR0913
    module: str,
    user: bool = False,
    sys_prefix: bool = False,
    overwrite: bool = True,
    symlink: bool = True,
    labextensions_dir: str | None = None,
    logger: logging.Logger | None = None,
) -> list[str]:
    """Develop a labextension bundled in a Python package.

    Returns a list of installed/updated directories.

    See develop_labextension for parameter information.
    """
    m, labexts = _get_labextension_metadata(module)
    base_path = str(Path(m.__file__).parent)

    full_dests = []

    for labext in labexts:
        src = str(Path(base_path) / labext["src"])
        dest = labext["dest"]
        if logger:
            logger.info("Installing %s -> %s", src, dest)

        if not Path(src).exists():
            build_labextension(base_path, logger=logger)

        full_dest = develop_labextension(
            src,
            overwrite=overwrite,
            symlink=symlink,
            user=user,
            sys_prefix=sys_prefix,
            labextensions_dir=labextensions_dir,
            destination=dest,
            logger=logger,
        )
        full_dests.append(full_dest)

    return full_dests


def build_labextension(  # noqa: PLR0913
    path: str | os.PathLike[str],
    logger: logging.Logger | None = None,
    development: bool = False,
    static_url: str | None = None,
    source_map: bool = False,
    core_version: str | None = None,
    core_package_file: str | None = None,
    core_path: str | None = None,
) -> None:
    """Build a labextension in the given path."""
    ext_path = str(Path(path).resolve())

    if core_path is not None:
        if logger:
            logger.warning(
                "\033[33m(Deprecated) `core_path` is deprecated and will be removed "
                "in a future release. Use `core_package_file` instead.\n \033[0m",
            )
        core_path_package = Path(core_path).resolve() / "package.json"
        if core_path_package.exists():
            core_package_file = str(core_path_package)

    core_package_file = core_package_file or get_core_meta(
        core_version,
        ext_path=ext_path,
        logger=logger,
    )
    core_package_file = str(Path(core_package_file).resolve())

    if logger:
        logger.info("Building extension in %s", path)

    builder, marker_pkg = _ensure_builder(ext_path, core_package_file)
    _check_node_version(builder, ext_path, logger=logger)

    if marker_pkg == "@jupyterlab/builder":
        core_flag = ["--core-path", _resolve_core_path_for_jupyterlab_builder(core_package_file)]
    else:
        core_flag = ["--core-package-file", core_package_file]

    node = _which_node_js()
    arguments = [node, builder, *core_flag, ext_path]
    if static_url is not None:
        arguments.extend(["--static-url", static_url])
    if development:
        arguments.append("--development")
    if source_map:
        arguments.append("--source-map")

    subprocess.check_call(arguments, cwd=ext_path)  # noqa: S603


def watch_labextension(  # noqa: PLR0913
    path: str | os.PathLike[str],
    labextensions_path: list[str],
    logger: logging.Logger | None = None,
    development: bool = False,
    source_map: bool = False,
    core_version: str | None = None,
    core_package_file: str | None = None,
) -> None:
    """Watch a labextension in a given path."""
    ext_path = str(Path(path).resolve())
    core_package_file = core_package_file or get_core_meta(
        core_version,
        ext_path=ext_path,
        logger=logger,
    )
    core_package_file = str(Path(core_package_file).resolve())

    if logger:
        logger.info("Building extension in %s", path)

    # Check to see if we need to create a symlink
    federated_extensions = get_federated_extensions(labextensions_path)

    with (Path(ext_path) / "package.json").open() as fid:
        ext_data = json.load(fid)

    if ext_data["name"] not in federated_extensions:
        develop_labextension_py(ext_path, sys_prefix=True)
    else:
        full_dest = str(
            Path(federated_extensions[ext_data["name"]]["ext_dir"]) / ext_data["name"],
        )
        output_dir = str(
            Path(ext_path) / ext_data["jupyterlab"].get("outputDir", "static"),
        )
        if not Path(full_dest).is_symlink():
            shutil.rmtree(full_dest)
            Path(full_dest).symlink_to(output_dir)

    builder, marker_pkg = _ensure_builder(ext_path, core_package_file)
    _check_node_version(builder, ext_path, logger=logger)

    if marker_pkg == "@jupyterlab/builder":
        core_flag = ["--core-path", _resolve_core_path_for_jupyterlab_builder(core_package_file)]
    else:
        core_flag = ["--core-package-file", core_package_file]

    node = _which_node_js()
    arguments = [node, builder, *core_flag, "--watch", ext_path]
    if development:
        arguments.append("--development")
    if source_map:
        arguments.append("--source-map")

    subprocess.check_call(arguments, cwd=ext_path)  # noqa: S603


# ------------------------------------------------------------------------------
# Private API
# ------------------------------------------------------------------------------


# Marker packages an extension may declare to identify its builder, in order
# of preference.
_BUILDER_MARKER_CANDIDATES = ("@jupyter/builder", "@jupyterlab/builder")

# Minimum Node.js range required by `@rspack/core` when its own `engines.node`
# field cannot be read. `@rspack/core` is a pure ES module that older Node.js
# versions cannot `require()`.
_FALLBACK_NODE_RANGE = "^20.19.0 || >=22.12.0"


def _read_rspack_node_range(builder: str, ext_path: str) -> str:
    """Return the ``engines.node`` range declared by ``@rspack/core``."""
    for root in (Path(builder).parent, Path(ext_path)):
        target = root
        while True:
            pkg = target / "node_modules" / "@rspack" / "core" / "package.json"
            if pkg.exists():
                try:
                    with pkg.open() as fid:
                        node_range = json.load(fid).get("engines", {}).get("node")
                except (OSError, ValueError):
                    node_range = None
                if node_range:
                    return str(node_range)
                break
            if target.parent == target:
                break
            target = target.parent
    return _FALLBACK_NODE_RANGE


def _check_node_version(
    builder: str,
    ext_path: str,
    logger: logging.Logger | None = None,
) -> None:
    """Fail early with a clear message when Node.js is too old to load the builder."""
    node = _which_node_js()
    node_range = _read_rspack_node_range(builder, ext_path)
    try:
        raw = subprocess.check_output([node, "--version"]).decode("utf8").strip()  # noqa: S603
    except (OSError, subprocess.CalledProcessError):
        return
    current = clean(raw, loose=True)  # type: ignore[no-untyped-call]
    if current is None or satisfies(current, node_range, loose=True):  # type: ignore[no-untyped-call]
        return
    msg = (
        f"Building this extension requires Node.js {node_range} (found {raw}). "
        "Please upgrade Node.js."
    )
    if logger:
        logger.error(msg)
    raise RuntimeError(msg)


def _resolve_core_path_for_jupyterlab_builder(core_package_file: str) -> str:
    """Return the core path directory for @jupyterlab/builder.

    @jupyterlab/builder's downstream script expects a file named package.json
    inside the core path directory. If core_package_file is not named
    package.json, a copy named package.json is created in the same directory.
    """
    core_file = Path(core_package_file)
    core_dir = core_file.parent

    if core_file.name != "package.json":
        target = core_dir / "package.json"
        if not target.exists():
            shutil.copy2(core_file, target)

    return str(core_dir)


def _select_builder_marker(ext_data: dict[str, Any]) -> tuple[str | None, str | None]:
    """Return (marker_pkg, dep_spec) for the builder marker the extension declares.

    Prefers `@jupyter/builder`. Returns (None, None) if neither marker is present.
    """
    for candidate in _BUILDER_MARKER_CANDIDATES:
        v = ext_data.get("devDependencies", {}).get(candidate)
        v = v or ext_data.get("dependencies", {}).get(candidate)
        if v is not None:
            return candidate, v
    return None, None


def _ensure_builder(ext_path: str, core_package_file: str) -> tuple[str, str]:
    """Ensure that we can build the extension and return ``(script, marker_pkg)``."""
    with Path(core_package_file).open() as fid:
        core_data = json.load(fid)
    with (Path(ext_path) / "package.json").open() as fid:
        ext_data = json.load(fid)

    marker_pkg, dep_version2 = _select_builder_marker(ext_data)
    if marker_pkg is None:
        msg = f"Extensions require a devDependency on {' or '.join(_BUILDER_MARKER_CANDIDATES)}"
        raise ValueError(msg)
    dep_version2 = cast("str", dep_version2)

    # if we have installed from disk (version is a path), assume we know what
    # we are doing and do not check versions.
    if "/" in dep_version2:
        with (Path(ext_path) / dep_version2 / "package.json").open() as fid:
            dep_version2 = str(json.load(fid).get("version", ""))
    if not (Path(ext_path) / "node_modules").exists():
        subprocess.check_call(["jlpm"], cwd=ext_path)  # noqa: S607

    # Find the marker package in node_modules using node module resolution.
    # We cannot use a script because the script path is a shell script on Windows.
    marker_parts = marker_pkg.split("/", 1)
    target = ext_path
    while not Path(target).joinpath("node_modules", *marker_parts).exists():
        if Path(target).parent == Path(target):
            msg = f"Could not find {marker_pkg}"
            raise ValueError(msg)
        target = str(Path(target).parent)

    # Check for compatible versions
    dep_version1 = core_data.get("devDependencies", {}).get(marker_pkg)
    if dep_version1 is not None:
        overlap = _test_overlap(
            dep_version1,
            dep_version2,
            drop_prerelease1=True,
            drop_prerelease2=True,
        )
        if not overlap:
            with Path(target).joinpath("node_modules", *marker_parts, "package.json").open() as fid:
                dep_version2 = str(json.load(fid).get("version", ""))
            overlap = _test_overlap(
                dep_version1,
                dep_version2,
                drop_prerelease1=True,
                drop_prerelease2=True,
            )

        if not overlap:
            msg = (
                f"Extensions require a devDependency on {marker_pkg}@{dep_version1}, "
                f"you have a dependency on {dep_version2}"
            )
            raise ValueError(msg)

    script = str(
        Path(target).joinpath("node_modules", *marker_parts, "lib", "build-labextension.js"),
    )
    return script, marker_pkg


def _should_copy(src: str, dest: str, logger: logging.Logger | None = None) -> bool:
    """Check whether a file should be copied because it is missing or the source is newer.

    Returns whether the file needs to be updated.

    Parameters
    ----------
    src : string
        A path that should exist from which to copy a file
    dest : string
        A path that might exist to which to copy a file
    logger : Jupyter logger [optional]
        Logger instance to use

    """
    if not Path(dest).exists():
        return True
    if Path(src).stat().st_mtime - Path(dest).stat().st_mtime > 1e-6:  # noqa: PLR2004
        # we add a fudge factor to work around a bug in python 2.x
        # that was fixed in python 3.x: https://bugs.python.org/issue12904
        if logger:
            logger.warning("Out of date: %s", dest)
        return True
    if logger:
        logger.info("Up to date: %s", dest)
    return False


def _maybe_copy(src: str, dest: str, logger: logging.Logger | None = None) -> None:
    """Copy a file if it needs updating.

    Parameters
    ----------
    src : string
        A path that should exist from which to copy a file
    dest : string
        A path that might exist to which to copy a file
    logger : Jupyter logger [optional]
        Logger instance to use

    """
    if _should_copy(src, dest, logger=logger):
        if logger:
            logger.info("Copying: %s -> %s", src, dest)
        shutil.copy2(src, dest)


def _get_labextension_dir(
    user: bool = False,
    sys_prefix: bool = False,
    prefix: str | None = None,
    labextensions_dir: str | None = None,
) -> str:
    """Return the labextension directory specified.

    Parameters
    ----------
    user : bool [default: False]
        Get the user's .jupyter/labextensions directory
    sys_prefix : bool [default: False]
        Get sys.prefix, i.e. ~/.envs/my-env/share/jupyter/labextensions
    prefix : str [optional]
        Get custom prefix
    labextensions_dir : str [optional]
        Get what you put in

    """
    conflicting = [
        ("user", user),
        ("prefix", prefix),
        ("labextensions_dir", labextensions_dir),
        ("sys_prefix", sys_prefix),
    ]
    conflicting_set = [f"{n}={v!r}" for n, v in conflicting if v]
    if len(conflicting_set) > 1:
        conflict = ", ".join(conflicting_set)
        msg = f"cannot specify more than one of user, sys_prefix, prefix, or labextensions_dir, but got: {conflict}"  # noqa: E501
        raise ArgumentConflict(msg)
    if user:
        labext = str(Path(jupyter_data_dir()) / "labextensions")
    elif sys_prefix:
        labext = str(Path(ENV_JUPYTER_PATH[0]) / "labextensions")
    elif prefix:
        labext = str(Path(prefix) / "share" / "jupyter" / "labextensions")
    elif labextensions_dir:
        labext = labextensions_dir
    else:
        labext = str(Path(SYSTEM_JUPYTER_PATH[0]) / "labextensions")
    return labext


# Directory names that are valid Python identifiers but never contain
# importable extension source we care about.
_SKIP_DIRS = frozenset({"__pycache__", "node_modules", "venv", "env"})


def _valid_package_dirs(dirs: list[str]) -> list[str]:
    """Filter out dirs that can never be Python package components."""
    return [d for d in dirs if d.isidentifier() and d not in _SKIP_DIRS]


def _find_packages(path: str) -> list[str]:
    """Find importable regular packages (dirs with ``__init__.py``) under *path*.

    Recursion only continues into directories that are themselves regular packages.
    Also prunes non-identifier names and common non-source directories.
    """
    path_obj = Path(path)
    packages: list[str] = []
    for root, dirs, files in os.walk(str(path_obj), followlinks=True):
        # Only keep descending into subdirectories that are themselves
        # packages; prune everything else.
        dirs[:] = [
            d for d in _valid_package_dirs(dirs) if (Path(root) / d / "__init__.py").is_file()
        ]
        rel = Path(root).relative_to(path_obj)
        if rel.parts and "__init__.py" in files:
            packages.append(".".join(rel.parts))
    return packages


def _find_namespace_packages(path: str) -> list[str]:
    """Find namespace packages (dirs with .py files, no __init__.py required)."""
    path_obj = Path(path)
    found: set[str] = set()
    for root, dirs, files in os.walk(str(path_obj), followlinks=True):
        dirs[:] = _valid_package_dirs(dirs)
        if any(f.endswith(".py") for f in files):
            rel = Path(root).relative_to(path_obj)
            if rel.parts:
                for i in range(len(rel.parts)):
                    found.add(".".join(rel.parts[: i + 1]))
    return list(found)


def _get_labextension_metadata(module: str) -> tuple[Any, list[dict[str, str]]]:  # noqa: C901
    """Get the list of labextension paths associated with a Python module.

    Returns a tuple of (the module path,             [{
        'src': 'mockextension',
        'dest': '_mockdestination'
    }])

    Parameters
    ----------
    module : str
        Importable Python module exposing the
        magic-named `_jupyter_labextension_paths` function

    """
    mod_path = str(Path(module).resolve())
    if not Path(mod_path).exists():
        msg = f"The path `{mod_path}` does not exist."
        raise FileNotFoundError(msg)

    errors = []

    # Check if the path is a valid labextension
    try:
        m = importlib.import_module(module)
        if hasattr(m, "_jupyter_labextension_paths"):
            return m, m._jupyter_labextension_paths()  # noqa: SLF001
    except Exception as exc:  # noqa: BLE001
        errors.append(exc)

    # Try to get the package name
    package = None

    # Try getting the package name from pyproject.toml
    if (Path(mod_path) / "pyproject.toml").exists():
        with (Path(mod_path) / "pyproject.toml").open("rb") as fid:
            data = load(fid)
        package = data.get("project", {}).get("name")

    # Try getting the package name from setup.py
    if not package:
        try:
            package = (
                subprocess.check_output(
                    [sys.executable, "setup.py", "--name"],
                    cwd=mod_path,
                )
                .decode("utf8")
                .strip()
            )
        except subprocess.CalledProcessError:
            msg = (
                f"The Python package `{module}` is not a valid package, "
                "it does not specify a `name` in `pyproject.toml` nor has a legacy `setup.py` file."
            )
            raise FileNotFoundError(msg) from None

    # Make sure the package is installed
    try:
        version(package)
    except PackageNotFoundError:
        subprocess.check_call([sys.executable, "-m", "pip", "install", "-e", mod_path])  # noqa: S603
        sys.path.insert(0, mod_path)

    package_candidates = [
        package.replace("-", "_"),  # Module with the same name as package
    ]
    package_candidates.extend(_find_packages(mod_path))  # Packages in the module path
    package_candidates.extend(
        _find_namespace_packages(mod_path),
    )  # Namespace packages in the module path

    for package in package_candidates:
        try:
            m = importlib.import_module(package)
            if hasattr(m, "_jupyter_labextension_paths"):
                return m, m._jupyter_labextension_paths()  # noqa: SLF001
        except Exception as exc:  # noqa: BLE001, PERF203
            errors.append(exc)

    msg = f"There is no labextension at {module}. Errors encountered: {errors}"
    raise ModuleNotFoundError(msg)
