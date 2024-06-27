"""Utilities for installing Javascript extensions for the notebook"""

# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.

import importlib
import json
import os
import os.path as osp
import platform
import shutil
import subprocess
import sys
from pathlib import Path

try:
    from importlib.metadata import PackageNotFoundError, version
except ImportError:
    from importlib_metadata import PackageNotFoundError, version

from os.path import basename, normpath
from os.path import join as pjoin

from jupyter_core.paths import ENV_JUPYTER_PATH, SYSTEM_JUPYTER_PATH, jupyter_data_dir
from jupyter_core.utils import ensure_dir_exists
from jupyter_server.extension.serverextension import ArgumentConflict
from jupyterlab_server.config import get_federated_extensions

try:
    from tomllib import load  # Python 3.11+
except ImportError:
    from tomli import load

from .commands import _test_overlap

DEPRECATED_ARGUMENT = object()

HERE = osp.abspath(osp.dirname(__file__))


# ------------------------------------------------------------------------------
# Public API
# ------------------------------------------------------------------------------


def develop_labextension(  # noqa
    path,
    symlink=True,
    overwrite=False,
    user=False,
    labextensions_dir=None,
    destination=None,
    logger=None,
    sys_prefix=False,
):
    """Install a prebuilt extension for JupyterLab

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
    """
    # the actual path to which we eventually installed
    full_dest = None

    labext = _get_labextension_dir(
        user=user, sys_prefix=sys_prefix, labextensions_dir=labextensions_dir
    )
    # make sure labextensions dir exists
    ensure_dir_exists(labext)

    if isinstance(path, (list, tuple)):
        msg = "path must be a string pointing to a single extension to install; call this function multiple times to install multiple extensions"
        raise TypeError(msg)

    if not destination:
        destination = basename(normpath(path))

    full_dest = normpath(pjoin(labext, destination))
    if overwrite and os.path.lexists(full_dest):
        if logger:
            logger.info("Removing: %s" % full_dest)
        if os.path.isdir(full_dest) and not os.path.islink(full_dest):
            shutil.rmtree(full_dest)
        else:
            os.remove(full_dest)

    # Make sure the parent directory exists
    os.makedirs(os.path.dirname(full_dest), exist_ok=True)

    if symlink:
        path = os.path.abspath(path)
        if not os.path.exists(full_dest):
            if logger:
                logger.info(f"Symlinking: {full_dest} -> {path}")
            try:
                os.symlink(path, full_dest)
            except OSError as e:
                if platform.platform().startswith("Windows"):
                    msg = (
                        "Symlinks can be activated on Windows 10 for Python version 3.8 or higher"
                        " by activating the 'Developer Mode'. That may not be allowed by your administrators.\n"
                        "See https://docs.microsoft.com/en-us/windows/apps/get-started/enable-your-device-for-development"
                    )
                    raise OSError(msg) from e
                raise

        elif not os.path.islink(full_dest):
            raise ValueError("%s exists and is not a symlink" % full_dest)

    elif os.path.isdir(path):
        path = pjoin(os.path.abspath(path), "")  # end in path separator
        for parent, _, files in os.walk(path):
            dest_dir = pjoin(full_dest, parent[len(path) :])
            if not os.path.exists(dest_dir):
                if logger:
                    logger.info("Making directory: %s" % dest_dir)
                os.makedirs(dest_dir)
            for file_name in files:
                src = pjoin(parent, file_name)
                dest_file = pjoin(dest_dir, file_name)
                _maybe_copy(src, dest_file, logger=logger)
    else:
        src = path
        _maybe_copy(src, full_dest, logger=logger)

    return full_dest


def develop_labextension_py(
    module,
    user=False,
    sys_prefix=False,
    overwrite=True,
    symlink=True,
    labextensions_dir=None,
    logger=None,
):
    """Develop a labextension bundled in a Python package.

    Returns a list of installed/updated directories.

    See develop_labextension for parameter information."""
    m, labexts = _get_labextension_metadata(module)
    base_path = os.path.split(m.__file__)[0]

    full_dests = []

    for labext in labexts:
        src = os.path.join(base_path, labext["src"])
        dest = labext["dest"]
        if logger:
            logger.info(f"Installing {src} -> {dest}")

        if not os.path.exists(src):
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


def build_labextension(
    path, logger=None, development=False, static_url=None, source_map=False, core_path=None
):
    """Build a labextension in the given path"""
    core_path = osp.join(HERE, "staging") if core_path is None else str(Path(core_path).resolve())

    ext_path = str(Path(path).resolve())

    if logger:
        logger.info("Building extension in %s" % path)

    builder = _ensure_builder(ext_path, core_path)

    arguments = ["node", builder, "--core-path", core_path, ext_path]
    if static_url is not None:
        arguments.extend(["--static-url", static_url])
    if development:
        arguments.append("--development")
    if source_map:
        arguments.append("--source-map")

    subprocess.check_call(arguments, cwd=ext_path)  # noqa S603


def watch_labextension(
    path, labextensions_path, logger=None, development=False, source_map=False, core_path=None
):
    """Watch a labextension in a given path"""
    core_path = osp.join(HERE, "staging") if core_path is None else str(Path(core_path).resolve())
    ext_path = str(Path(path).resolve())

    if logger:
        logger.info("Building extension in %s" % path)

    # Check to see if we need to create a symlink
    federated_extensions = get_federated_extensions(labextensions_path)

    with open(pjoin(ext_path, "package.json")) as fid:
        ext_data = json.load(fid)

    if ext_data["name"] not in federated_extensions:
        develop_labextension_py(ext_path, sys_prefix=True)
    else:
        full_dest = pjoin(federated_extensions[ext_data["name"]]["ext_dir"], ext_data["name"])
        output_dir = pjoin(ext_path, ext_data["jupyterlab"].get("outputDir", "static"))
        if not osp.islink(full_dest):
            shutil.rmtree(full_dest)
            os.symlink(output_dir, full_dest)

    builder = _ensure_builder(ext_path, core_path)
    arguments = ["node", builder, "--core-path", core_path, "--watch", ext_path]
    if development:
        arguments.append("--development")
    if source_map:
        arguments.append("--source-map")

    subprocess.check_call(arguments, cwd=ext_path)  # noqa S603


# ------------------------------------------------------------------------------
# Private API
# ------------------------------------------------------------------------------


def _ensure_builder(ext_path, core_path):
    """Ensure that we can build the extension and return the builder script path"""
    # Test for compatible dependency on @jupyterlab/builder
    with open(osp.join(core_path, "package.json")) as fid:
        core_data = json.load(fid)
    with open(osp.join(ext_path, "package.json")) as fid:
        ext_data = json.load(fid)
    dep_version1 = core_data["devDependencies"]["@jupyterlab/builder"]
    dep_version2 = ext_data.get("devDependencies", {}).get("@jupyterlab/builder")
    dep_version2 = dep_version2 or ext_data.get("dependencies", {}).get("@jupyterlab/builder")
    if dep_version2 is None:
        raise ValueError(
            "Extensions require a devDependency on @jupyterlab/builder@%s" % dep_version1
        )

    # if we have installed from disk (version is a path), assume we know what
    # we are doing and do not check versions.
    if "/" in dep_version2:
        with open(osp.join(ext_path, dep_version2, "package.json")) as fid:
            dep_version2 = json.load(fid).get("version")
    if not osp.exists(osp.join(ext_path, "node_modules")):
        subprocess.check_call(["jlpm"], cwd=ext_path)  # noqa S603 S607

    # Find @jupyterlab/builder using node module resolution
    # We cannot use a script because the script path is a shell script on Windows
    target = ext_path
    while not osp.exists(osp.join(target, "node_modules", "@jupyterlab", "builder")):
        if osp.dirname(target) == target:
            msg = "Could not find @jupyterlab/builder"
            raise ValueError(msg)
        target = osp.dirname(target)

    overlap = _test_overlap(
        dep_version1, dep_version2, drop_prerelease1=True, drop_prerelease2=True
    )
    if not overlap:
        with open(
            osp.join(target, "node_modules", "@jupyterlab", "builder", "package.json")
        ) as fid:
            dep_version2 = json.load(fid).get("version")
        overlap = _test_overlap(
            dep_version1, dep_version2, drop_prerelease1=True, drop_prerelease2=True
        )

    if not overlap:
        msg = f"Extensions require a devDependency on @jupyterlab/builder@{dep_version1}, you have a dependency on {dep_version2}"
        raise ValueError(msg)

    return osp.join(
        target, "node_modules", "@jupyterlab", "builder", "lib", "build-labextension.js"
    )


def _should_copy(src, dest, logger=None):
    """Should a file be copied, if it doesn't exist, or is newer?

    Returns whether the file needs to be updated.

    Parameters
    ----------

    src : string
        A path that should exist from which to copy a file
    src : string
        A path that might exist to which to copy a file
    logger : Jupyter logger [optional]
        Logger instance to use
    """
    if not os.path.exists(dest):
        return True
    if os.stat(src).st_mtime - os.stat(dest).st_mtime > 1e-6:  # noqa
        # we add a fudge factor to work around a bug in python 2.x
        # that was fixed in python 3.x: https://bugs.python.org/issue12904
        if logger:
            logger.warning("Out of date: %s" % dest)
        return True
    if logger:
        logger.info("Up to date: %s" % dest)
    return False


def _maybe_copy(src, dest, logger=None):
    """Copy a file if it needs updating.

    Parameters
    ----------

    src : string
        A path that should exist from which to copy a file
    src : string
        A path that might exist to which to copy a file
    logger : Jupyter logger [optional]
        Logger instance to use
    """
    if _should_copy(src, dest, logger=logger):
        if logger:
            logger.info(f"Copying: {src} -> {dest}")
        shutil.copy2(src, dest)


def _get_labextension_dir(user=False, sys_prefix=False, prefix=None, labextensions_dir=None):
    """Return the labextension directory specified

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
        msg = "cannot specify more than one of user, sys_prefix, prefix, or labextensions_dir, but got: {}".format(
            ", ".join(conflicting_set)
        )
        raise ArgumentConflict(msg)
    if user:
        labext = pjoin(jupyter_data_dir(), "labextensions")
    elif sys_prefix:
        labext = pjoin(ENV_JUPYTER_PATH[0], "labextensions")
    elif prefix:
        labext = pjoin(prefix, "share", "jupyter", "labextensions")
    elif labextensions_dir:
        labext = labextensions_dir
    else:
        labext = pjoin(SYSTEM_JUPYTER_PATH[0], "labextensions")
    return labext


def _get_labextension_metadata(module):  # noqa
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
    mod_path = osp.abspath(module)
    if not osp.exists(mod_path):
        msg = f"The path `{mod_path}` does not exist."
        raise FileNotFoundError(msg)

    errors = []

    # Check if the path is a valid labextension
    try:
        m = importlib.import_module(module)
        if hasattr(m, "_jupyter_labextension_paths"):
            return m, m._jupyter_labextension_paths()
    except Exception as exc:
        errors.append(exc)

    # Try to get the package name
    package = None

    # Try getting the package name from pyproject.toml
    if os.path.exists(os.path.join(mod_path, "pyproject.toml")):
        with open(os.path.join(mod_path, "pyproject.toml"), "rb") as fid:
            data = load(fid)
        package = data.get("project", {}).get("name")

    # Try getting the package name from setup.py
    if not package:
        try:
            package = (
                subprocess.check_output(
                    [sys.executable, "setup.py", "--name"],  # noqa S603
                    cwd=mod_path,
                )
                .decode("utf8")
                .strip()
            )
        except subprocess.CalledProcessError:
            msg = (
                f"The Python package `{module}` is not a valid package, "
                "it is missing the `setup.py` file."
            )
            raise FileNotFoundError(msg) from None

    # Make sure the package is installed
    try:
        version(package)
    except PackageNotFoundError:
        subprocess.check_call([sys.executable, "-m", "pip", "install", "-e", mod_path])  # noqa S603
        sys.path.insert(0, mod_path)

    from setuptools import find_namespace_packages, find_packages

    package_candidates = [
        package.replace("-", "_"),  # Module with the same name as package
    ]
    package_candidates.extend(find_packages(mod_path))  # Packages in the module path
    package_candidates.extend(
        find_namespace_packages(mod_path)
    )  # Namespace packages in the module path

    for package in package_candidates:
        try:
            m = importlib.import_module(package)
            if hasattr(m, "_jupyter_labextension_paths"):
                return m, m._jupyter_labextension_paths()
        except Exception as exc:
            errors.append(exc)

    msg = f"There is no labextension at {module}. Errors encountered: {errors}"
    raise ModuleNotFoundError(msg)
