"""JupyterLab command handler"""

# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.
import contextlib
import errno
import hashlib
import itertools
import json
import logging
import os
import os.path as osp
import re
import shutil
import site
import stat
import subprocess
import sys
import tarfile
from copy import deepcopy
from dataclasses import dataclass
from glob import glob
from pathlib import Path
from tempfile import TemporaryDirectory
from threading import Event
from typing import Optional
from urllib.error import URLError
from urllib.request import Request, quote, urljoin, urlopen

from jupyter_core.paths import jupyter_config_dir
from jupyter_server.extension.serverextension import GREEN_ENABLED, GREEN_OK, RED_DISABLED, RED_X
from jupyterlab_server.config import (
    get_allowed_levels,
    get_federated_extensions,
    get_package_url,
    get_page_config,
    get_static_page_config,
    write_page_config,
)
from jupyterlab_server.process import Process, WatchHelper, list2cmdline, which
from packaging.version import Version
from traitlets import Bool, HasTraits, Instance, List, Unicode, default

from jupyterlab._version import __version__
from jupyterlab.coreconfig import CoreConfig
from jupyterlab.jlpmapp import HERE, YARN_PATH
from jupyterlab.semver import Range, gt, gte, lt, lte, make_semver

# The regex for expecting the webpack output.
WEBPACK_EXPECT = re.compile(r".*theme-light-extension/style/theme.css")


# The repo root directory
REPO_ROOT = osp.abspath(osp.join(HERE, ".."))


# The dev mode directory.
DEV_DIR = osp.join(REPO_ROOT, "dev_mode")


# If we are pinning the package, rename it `pin@<alias>`
PIN_PREFIX = "pin@"


# Default Yarn registry used in default yarn.lock
YARN_DEFAULT_REGISTRY = "https://registry.yarnpkg.com"


class ProgressProcess(Process):
    def __init__(self, cmd, logger=None, cwd=None, kill_event=None, env=None):
        """Start a subprocess that can be run asynchronously.

        Parameters
        ----------
        cmd: list
            The command to run.
        logger: :class:`~logger.Logger`, optional
            The logger instance.
        cwd: string, optional
            The cwd of the process.
        kill_event: :class:`~threading.Event`, optional
            An event used to kill the process operation.
        env: dict, optional
            The environment for the process.
        """
        if not isinstance(cmd, (list, tuple)):
            msg = "Command must be given as a list"
            raise ValueError(msg)

        if kill_event and kill_event.is_set():
            msg = "Process aborted"
            raise ValueError(msg)

        self.logger = _ensure_logger(logger)
        self._last_line = ""
        self.cmd = cmd
        self.logger.debug(f"> {list2cmdline(cmd)}")

        self.proc = self._create_process(
            cwd=cwd,
            env=env,
            stderr=subprocess.STDOUT,
            stdout=subprocess.PIPE,
            universal_newlines=True,
            encoding="utf-8",
        )
        self._kill_event = kill_event or Event()

        Process._procs.add(self)

    def wait(self):
        cache = []
        proc = self.proc
        kill_event = self._kill_event
        spinner = itertools.cycle(["-", "\\", "|", "/"])
        while proc.poll() is None:
            sys.stdout.write(next(spinner))  # write the next character
            sys.stdout.flush()  # flush stdout buffer (actual character display)
            sys.stdout.write("\b")
            if kill_event.is_set():
                self.terminate()
                msg = "Process was aborted"
                raise ValueError(msg)
            try:
                out, _ = proc.communicate(timeout=0.1)
                cache.append(out)
            except subprocess.TimeoutExpired:
                continue
        self.logger.debug("\n".join(cache))
        sys.stdout.flush()
        return self.terminate()


def pjoin(*args):
    """Join paths to create a real path."""
    return osp.abspath(osp.join(*args))


def get_user_settings_dir():
    """Get the configured JupyterLab user settings directory."""
    settings_dir = os.environ.get("JUPYTERLAB_SETTINGS_DIR")
    settings_dir = settings_dir or pjoin(jupyter_config_dir(), "lab", "user-settings")
    return osp.abspath(settings_dir)


def get_workspaces_dir():
    """Get the configured JupyterLab workspaces directory."""
    workspaces_dir = os.environ.get("JUPYTERLAB_WORKSPACES_DIR")
    workspaces_dir = workspaces_dir or pjoin(jupyter_config_dir(), "lab", "workspaces")
    return osp.abspath(workspaces_dir)


def get_app_dir():
    """Get the configured JupyterLab app directory."""
    # Default to the override environment variable.
    if os.environ.get("JUPYTERLAB_DIR"):
        # We must resolve the path to get the canonical case of the path for
        # case-sensitive systems
        return str(Path(os.environ["JUPYTERLAB_DIR"]).resolve())

    # Use the default locations for data_files.
    app_dir = pjoin(sys.prefix, "share", "jupyter", "lab")

    # Check for a user level install.
    # Ensure that USER_BASE is defined
    if hasattr(site, "getuserbase"):
        site.getuserbase()
    userbase = getattr(site, "USER_BASE", None)
    if HERE.startswith(userbase) and not app_dir.startswith(userbase):
        app_dir = pjoin(userbase, "share", "jupyter", "lab")

    # Check for a system install in '/usr/local/share'.
    elif (
        sys.prefix.startswith("/usr")
        and not osp.exists(app_dir)
        and osp.exists("/usr/local/share/jupyter/lab")
    ):
        app_dir = "/usr/local/share/jupyter/lab"

    # We must resolve the path to get the canonical case of the path for
    # case-sensitive systems
    return str(Path(app_dir).resolve())


def dedupe_yarn(path, logger=None):
    """`yarn-deduplicate` with the `fewer` strategy to minimize total
    packages installed in a given staging directory

    This means a extension (or dependency) _could_ cause a downgrade of an
    version expected at publication time, but core should aggressively set
    pins above, for example, known-bad versions
    """
    had_dupes = (
        ProgressProcess(
            [
                "node",
                YARN_PATH,
                "dlx",
                "yarn-berry-deduplicate",
                "-s",
                "fewerHighest",
                "--fail",
            ],
            cwd=path,
            logger=logger,
        ).wait()
        != 0
    )

    if had_dupes:
        yarn_proc = ProgressProcess(["node", YARN_PATH], cwd=path, logger=logger)
        yarn_proc.wait()


def ensure_node_modules(cwd, logger=None):
    """Ensure that node_modules is up to date.

    Returns true if the node_modules was updated.
    """
    logger = _ensure_logger(logger)
    yarn_proc = ProgressProcess(
        ["node", YARN_PATH, "--immutable", "--immutable-cache"], cwd=cwd, logger=logger
    )
    ret = yarn_proc.wait()

    # Update node_modules if needed.
    if ret != 0:
        yarn_proc = ProgressProcess(["node", YARN_PATH], cwd=cwd, logger=logger)
        yarn_proc.wait()
        dedupe_yarn(REPO_ROOT, logger)

    return ret != 0


def ensure_dev(logger=None):
    """Ensure that the dev assets are available."""
    logger = _ensure_logger(logger)
    target = pjoin(DEV_DIR, "static")

    # Determine whether to build.
    if ensure_node_modules(REPO_ROOT, logger) or not osp.exists(target):
        yarn_proc = ProgressProcess(["node", YARN_PATH, "build"], cwd=REPO_ROOT, logger=logger)
        yarn_proc.wait()


def ensure_core(logger=None):
    """Ensure that the core assets are available."""
    staging = pjoin(HERE, "staging")
    logger = _ensure_logger(logger)

    # Determine whether to build.
    target = pjoin(HERE, "static", "index.html")
    if not osp.exists(target):
        ensure_node_modules(staging, logger)
        yarn_proc = ProgressProcess(["node", YARN_PATH, "build"], cwd=staging, logger=logger)
        yarn_proc.wait()


def ensure_app(app_dir):
    """Ensure that an application directory is available.

    If it does not exist, return a list of messages to prompt the user.
    """
    if osp.exists(pjoin(app_dir, "static", "index.html")):
        return

    msgs = [
        f'JupyterLab application assets not found in "{app_dir}"',
        "Please run `jlpm run build:core` then `jupyter lab build` ",
        "or use a different app directory",
    ]
    return msgs


def watch_packages(logger=None):
    """Run watch mode for the source packages.

    Parameters
    ----------
    logger: :class:`~logger.Logger`, optional
        The logger instance.

    Returns
    -------
    A list of `WatchHelper` objects.
    """
    logger = _ensure_logger(logger)
    ensure_node_modules(REPO_ROOT, logger)

    ts_dir = osp.abspath(osp.join(REPO_ROOT, "packages", "metapackage"))

    # Run typescript watch and wait for the string indicating it is done.
    ts_regex = r".* Found 0 errors\. Watching for file changes\."
    ts_proc = WatchHelper(
        ["node", YARN_PATH, "run", "watch"], cwd=ts_dir, logger=logger, startup_regex=ts_regex
    )

    return [ts_proc]


def watch_dev(logger=None):
    """Run watch mode in a given directory.

    Parameters
    ----------
    logger: :class:`~logger.Logger`, optional
        The logger instance.

    Returns
    -------
    A list of `WatchHelper` objects.
    """
    logger = _ensure_logger(logger)

    package_procs = watch_packages(logger)

    # Run webpack watch and wait for compilation.
    wp_proc = WatchHelper(
        ["node", YARN_PATH, "run", "watch"],
        cwd=DEV_DIR,
        logger=logger,
        startup_regex=WEBPACK_EXPECT,
    )

    return [*package_procs, wp_proc]


class AppOptions(HasTraits):
    """Options object for build system"""

    def __init__(self, logger=None, core_config=None, **kwargs):
        if core_config is not None:
            kwargs["core_config"] = core_config
        if logger is not None:
            kwargs["logger"] = logger

        # use the default if app_dir is empty
        if "app_dir" in kwargs and not kwargs["app_dir"]:
            kwargs.pop("app_dir")

        super().__init__(**kwargs)

    app_dir = Unicode(help="The application directory")

    use_sys_dir = Bool(
        True,
        help=("Whether to shadow the default app_dir if that is set to a non-default value"),
    )

    logger = Instance(logging.Logger, help="The logger to use")

    core_config = Instance(CoreConfig, help="Configuration for core data")

    kill_event = Instance(Event, args=(), help="Event for aborting call")

    labextensions_path = List(
        Unicode(), help="The paths to look in for prebuilt JupyterLab extensions"
    )

    registry = Unicode(help="NPM packages registry URL")

    splice_source = Bool(False, help="Splice source packages into app directory.")

    skip_full_build_check = Bool(
        False,
        help=(
            "If true, perform only a quick check that the lab build is up to date."
            " If false, perform a thorough check, which verifies extension contents."
        ),
    )

    verbose = Bool(False, help="Increase verbosity level.")

    @default("logger")
    def _default_logger(self):
        return logging.getLogger("jupyterlab")

    # These defaults need to be federated to pick up
    # any changes to env vars:
    @default("app_dir")
    def _default_app_dir(self):
        return get_app_dir()

    @default("core_config")
    def _default_core_config(self):
        return CoreConfig()

    @default("registry")
    def _default_registry(self):
        config = _yarn_config(self.logger)["yarn config"]
        return config.get("registry", YARN_DEFAULT_REGISTRY)


def _ensure_options(options):
    """Helper to use deprecated kwargs for AppOption"""
    if options is None:
        return AppOptions()
    elif issubclass(options.__class__, AppOptions):
        return options
    else:
        return AppOptions(**options)


def watch(app_options=None):
    """Watch the application.

    Parameters
    ----------
    app_options: :class:`AppOptions`, optional
        The application options.

    Returns
    -------
    A list of processes to run asynchronously.
    """
    app_options = _ensure_options(app_options)
    _node_check(app_options.logger)
    handler = _AppHandler(app_options)

    package_procs = watch_packages(app_options.logger) if app_options.splice_source else []

    return package_procs + handler.watch()


def install_extension(extension, app_options=None, pin=None):
    """Install an extension package into JupyterLab.

    The extension is first validated.

    Returns `True` if a rebuild is recommended, `False` otherwise.
    """
    app_options = _ensure_options(app_options)
    _node_check(app_options.logger)
    handler = _AppHandler(app_options)
    return handler.install_extension(extension, pin=pin)


def uninstall_extension(name=None, app_options=None, all_=False):
    """Uninstall an extension by name or path.

    Returns `True` if a rebuild is recommended, `False` otherwise.
    """
    app_options = _ensure_options(app_options)
    _node_check(app_options.logger)
    handler = _AppHandler(app_options)
    if all_ is True:
        return handler.uninstall_all_extensions()
    return handler.uninstall_extension(name)


def update_extension(name=None, all_=False, app_dir=None, app_options=None):
    """Update an extension by name, or all extensions.
    Either `name` must be given as a string, or `all_` must be `True`.
    If `all_` is `True`, the value of `name` is ignored.
    Returns `True` if a rebuild is recommended, `False` otherwise.
    """
    app_options = _ensure_options(app_options)
    _node_check(app_options.logger)
    handler = _AppHandler(app_options)
    if all_ is True:
        return handler.update_all_extensions()
    return handler.update_extension(name)


def clean(app_options=None):
    """Clean the JupyterLab application directory."""
    app_options = _ensure_options(app_options)
    logger = app_options.logger
    app_dir = app_options.app_dir

    logger.info(f"Cleaning {app_dir}...")
    if app_dir == pjoin(HERE, "dev"):
        msg = "Cannot clean the dev app"
        raise ValueError(msg)
    if app_dir == pjoin(HERE, "core"):
        msg = "Cannot clean the core app"
        raise ValueError(msg)

    if getattr(app_options, "all", False):
        logger.info(f"Removing everything in {app_dir}...")
        _rmtree_star(app_dir, logger)
    else:
        possible_targets = ["extensions", "settings", "staging", "static"]
        targets = [t for t in possible_targets if getattr(app_options, t)]

        for name in targets:
            target = pjoin(app_dir, name)
            if osp.exists(target):
                logger.info(f"Removing {name}...")
                _rmtree(target, logger)
            else:
                logger.info(f"{name} not present, skipping...")

    logger.info("Success!")
    if getattr(app_options, "all", False) or getattr(app_options, "extensions", False):
        logger.info("All of your extensions have been removed, and will need to be reinstalled")


def build(
    name=None,
    version=None,
    static_url=None,
    kill_event=None,
    clean_staging=False,
    app_options=None,
    production=True,
    minimize=True,
):
    """Build the JupyterLab application."""
    app_options = _ensure_options(app_options)
    _node_check(app_options.logger)
    handler = _AppHandler(app_options)
    return handler.build(
        name=name,
        version=version,
        static_url=static_url,
        production=production,
        minimize=minimize,
        clean_staging=clean_staging,
    )


def get_app_info(app_options=None):
    """Get a dictionary of information about the app."""
    handler = _AppHandler(app_options)
    handler._ensure_disabled_info()
    return handler.info


def enable_extension(extension, app_options=None, level="sys_prefix"):
    """Enable a JupyterLab extension/plugin.

    Returns `True` if a rebuild is recommended, `False` otherwise.
    """
    handler = _AppHandler(app_options)
    return handler.toggle_extension(extension, False, level=level)


def disable_extension(extension, app_options=None, level="sys_prefix"):
    """Disable a JupyterLab extension/plugin.

    Returns `True` if a rebuild is recommended, `False` otherwise.
    """
    handler = _AppHandler(app_options)
    return handler.toggle_extension(extension, True, level=level)


def check_extension(extension, installed=False, app_options=None):
    """Check if a JupyterLab extension is enabled or disabled."""
    handler = _AppHandler(app_options)
    return handler.check_extension(extension, installed)


def lock_extension(extension, app_options=None, level="sys_prefix"):
    """Lock a JupyterLab extension/plugin."""
    handler = _AppHandler(app_options)
    return handler.toggle_extension_lock(extension, True, level=level)


def unlock_extension(extension, app_options=None, level="sys_prefix"):
    """Unlock a JupyterLab extension/plugin."""
    handler = _AppHandler(app_options)
    return handler.toggle_extension_lock(extension, False, level=level)


def build_check(app_options=None):
    """Determine whether JupyterLab should be built.

    Returns a list of messages.
    """
    app_options = _ensure_options(app_options)
    _node_check(app_options.logger)
    handler = _AppHandler(app_options)
    return handler.build_check()


def list_extensions(app_options=None):
    """List the extensions."""
    handler = _AppHandler(app_options)
    return handler.list_extensions()


def link_package(path, app_options=None):
    """Link a package against the JupyterLab build.

    Returns `True` if a rebuild is recommended, `False` otherwise.
    """
    handler = _AppHandler(app_options)
    return handler.link_package(path)


def unlink_package(package, app_options=None):
    """Unlink a package from JupyterLab by path or name.

    Returns `True` if a rebuild is recommended, `False` otherwise.
    """
    handler = _AppHandler(app_options)
    return handler.unlink_package(package)


def get_app_version(app_options=None):
    """Get the application version."""
    handler = _AppHandler(app_options)
    return handler.info["version"]


def get_latest_compatible_package_versions(names, app_options=None):
    """Get the latest compatible version of a list of packages."""
    handler = _AppHandler(app_options)
    return handler.latest_compatible_package_versions(names)


def read_package(target):
    """Read the package data in a given target tarball."""
    with tarfile.open(target, "r") as tar:
        with tar.extractfile("package/package.json") as f:
            data = json.loads(f.read().decode("utf8"))
        data["jupyterlab_extracted_files"] = [f.path[len("package/") :] for f in tar.getmembers()]
    return data


# ----------------------------------------------------------------------
# Implementation details
# ----------------------------------------------------------------------


class _AppHandler:
    def __init__(self, options):
        """Create a new _AppHandler object"""
        options = _ensure_options(options)
        self._options = options
        self.app_dir = options.app_dir
        self.sys_dir = get_app_dir() if options.use_sys_dir else self.app_dir
        self.logger = options.logger
        # Make a deep copy of the core data so we don't influence the original copy
        self.core_data = deepcopy(options.core_config._data)
        self.labextensions_path = options.labextensions_path
        self.verbose = options.verbose
        self.kill_event = options.kill_event
        self.registry = options.registry
        self.skip_full_build_check = options.skip_full_build_check

        # Do this last since it relies on other attributes
        self.info = self._get_app_info()
        # Migrate from 4.0 which did not have "locked" status
        try:
            self._maybe_mirror_disabled_in_locked(level="sys_prefix")
        except (PermissionError, OSError):
            try:
                self.logger.info(
                    "`sys_prefix` level settings are read-only, using `user` level for migration to `lockedExtensions`"
                )
                self._maybe_mirror_disabled_in_locked(level="user")
            except (PermissionError, OSError):
                self.logger.warning(
                    "Both `sys_prefix` and `user` level settings are read-only, cannot auto-migrate `disabledExtensions` to `lockedExtensions`"
                )

    def install_extension(self, extension, existing=None, pin=None):
        """Install an extension package into JupyterLab.

        The extension is first validated.

        Returns `True` if a rebuild is recommended, `False` otherwise.
        """
        extension = _normalize_path(extension)
        extensions = self.info["extensions"]

        # Check for a core extensions.
        if extension in self.info["core_extensions"]:
            config = self._read_build_config()
            uninstalled = config.get("uninstalled_core_extensions", [])
            if extension in uninstalled:
                self.logger.info(f"Installing core extension {extension}")
                uninstalled.remove(extension)
                config["uninstalled_core_extensions"] = uninstalled
                self._write_build_config(config)
                return True
            return False

        # Create the app dirs if needed.
        self._ensure_app_dirs()

        # Install the package using a temporary directory.
        with TemporaryDirectory() as tempdir:
            info = self._install_extension(extension, tempdir, pin=pin)

        name = info["name"]

        # Local directories get name mangled and stored in metadata.
        if info["is_dir"]:
            config = self._read_build_config()
            local = config.setdefault("local_extensions", {})
            local[name] = info["source"]
            self._write_build_config(config)

        # Remove an existing extension with the same name and different path
        if name in extensions:
            other = extensions[name]
            if other["path"] != info["path"] and other["location"] == "app":
                os.remove(other["path"])

        return True

    def build(
        self,
        name=None,
        version=None,
        static_url=None,
        clean_staging=False,
        production=True,
        minimize=True,
    ):
        """Build the application."""
        if production is None:
            production = not (self.info["linked_packages"] or self.info["local_extensions"])

        if not production:
            minimize = False

        # If splicing, make sure the source packages are built
        if self._options.splice_source:
            ensure_node_modules(REPO_ROOT, logger=self.logger)
            self._run(["node", YARN_PATH, "build:packages"], cwd=REPO_ROOT)

        info = ["production" if production else "development"]
        if production:
            info.append("minimized" if minimize else "not minimized")
        self.logger.info(f"Building jupyterlab assets ({', '.join(info)})")

        # Set up the build directory.
        app_dir = self.app_dir

        self._populate_staging(
            name=name, version=version, static_url=static_url, clean=clean_staging
        )

        staging = pjoin(app_dir, "staging")

        # Make sure packages are installed.
        ret = self._run(["node", YARN_PATH, "install"], cwd=staging)
        if ret != 0:
            msg = "npm dependencies failed to install"
            self.logger.debug(msg)
            raise RuntimeError(msg)

        # Build the app.
        dedupe_yarn(staging, self.logger)
        command = f"build:{'prod' if production else 'dev'}{':minimize' if minimize else ''}"
        ret = self._run(["node", YARN_PATH, "run", command], cwd=staging)
        if ret != 0:
            msg = "JupyterLab failed to build"
            self.logger.debug(msg)
            raise RuntimeError(msg)

    def watch(self):
        """Start the application watcher and then run the watch in
        the background.
        """
        staging = pjoin(self.app_dir, "staging")

        self._populate_staging()

        # Make sure packages are installed.
        self._run(["node", YARN_PATH, "install"], cwd=staging)
        dedupe_yarn(staging, self.logger)

        proc = WatchHelper(
            ["node", YARN_PATH, "run", "watch"],
            cwd=pjoin(self.app_dir, "staging"),
            startup_regex=WEBPACK_EXPECT,
            logger=self.logger,
        )
        return [proc]

    def list_extensions(self):  # noqa
        """Print an output of the extensions."""
        self._ensure_disabled_info()
        logger = self.logger
        info = self.info
        version = info["version"]
        logger.info(f"JupyterLab v{version}")

        if info["federated_extensions"] or info["extensions"]:
            info["compat_errors"] = self._get_extension_compat()

        if info["federated_extensions"]:
            self._list_federated_extensions()

        if info["extensions"]:
            logger.info("Other labextensions (built into JupyterLab)")
            self._list_extensions(info, "app")
            self._list_extensions(info, "sys")

        local = info["local_extensions"]
        if local:
            logger.info("\n   local extensions:")
            for name in sorted(local):
                logger.info(f"        {name}: {local[name]}")

        linked_packages = info["linked_packages"]
        if linked_packages:
            logger.info("\n   linked packages:")
            for key in sorted(linked_packages):
                source = linked_packages[key]["source"]
                logger.info(f"        {key}: {source}")

        uninstalled_core = info["uninstalled_core"]
        if uninstalled_core:
            logger.info("\nUninstalled core extensions:")
            [logger.info(f"    {item}") for item in sorted(uninstalled_core)]

        all_exts = (
            list(info["federated_extensions"])
            + list(info["extensions"])
            + list(info["core_extensions"])
        )
        # Ignore disabled extensions that are not installed
        disabled = [i for i in info["disabled"] if i.partition(":")[0] in all_exts]
        if disabled:
            logger.info("\nDisabled extensions:")
            for item in sorted(disabled):
                # Show that all plugins will be disabled if the whole extension matches
                if item in all_exts:
                    item += " (all plugins)"  # noqa PLW2901
                logger.info(f"    {item}")

        # Here check if modules are improperly shadowed
        improper_shadowed = []
        for ext_name in self.info["shadowed_exts"]:
            source_version = self.info["extensions"][ext_name]["version"]
            prebuilt_version = self.info["federated_extensions"][ext_name]["version"]
            if not gte(prebuilt_version, source_version, True):
                improper_shadowed.append(ext_name)

        if improper_shadowed:
            logger.info(
                "\nThe following source extensions are overshadowed by older prebuilt extensions:"
            )
            [logger.info(f"    {name}") for name in sorted(improper_shadowed)]

        messages = self.build_check(fast=True)
        if messages:
            logger.info("\nBuild recommended, please run `jupyter lab build`:")
            [logger.info(f"    {item}") for item in messages]

    def build_check(self, fast=None):  # noqa
        """Determine whether JupyterLab should be built.

        Returns a list of messages.
        """
        if fast is None:
            fast = self.skip_full_build_check
        app_dir = self.app_dir
        local = self.info["local_extensions"]
        linked = self.info["linked_packages"]
        messages = []

        # Check for no application.
        pkg_path = pjoin(app_dir, "static", "package.json")
        if not osp.exists(pkg_path):
            return ["No built application"]

        static_data = self.info["static_data"]
        old_jlab = static_data["jupyterlab"]
        old_deps = static_data.get("dependencies", {})

        # Look for mismatched version.
        static_version = old_jlab.get("version", "")
        if not static_version.endswith("-spliced"):
            core_version = old_jlab["version"]
            if Version(static_version) != Version(core_version):
                msg = f"Version mismatch: {static_version} (built), {core_version} (current)"
                return [msg]

        shadowed_exts = self.info["shadowed_exts"]

        # Look for mismatched extensions.
        new_package = self._get_package_template(silent=fast)
        new_jlab = new_package["jupyterlab"]
        new_deps = new_package.get("dependencies", {})

        for ext_type in ["extensions", "mimeExtensions"]:
            # Extensions that were added.
            for ext in new_jlab[ext_type]:
                if ext in shadowed_exts:
                    continue
                if ext not in old_jlab[ext_type]:
                    messages.append(f"{ext} needs to be included in build")

            # Extensions that were removed.
            for ext in old_jlab[ext_type]:
                if ext in shadowed_exts:
                    continue
                if ext not in new_jlab[ext_type]:
                    messages.append(f"{ext} needs to be removed from build")

        # Look for mismatched dependencies
        src_pkg_dir = pjoin(REPO_ROOT, "packages")
        for pkg, dep in new_deps.items():
            if old_deps.get(pkg, "").startswith(src_pkg_dir):
                continue
            if pkg not in old_deps:
                continue
            # Skip local and linked since we pick them up separately.
            if pkg in local or pkg in linked:
                continue
            if old_deps[pkg] != dep:
                msg = f"{pkg} changed from {old_deps[pkg]} to {new_deps[pkg]}"
                messages.append(msg)

        # Look for updated local extensions.
        for name, source in local.items():
            if fast or name in shadowed_exts:
                continue
            dname = pjoin(app_dir, "extensions")
            if self._check_local(name, source, dname):
                messages.append(f"{name} content changed")

        # Look for updated linked packages.
        for name, item in linked.items():
            if fast or name in shadowed_exts:
                continue
            dname = pjoin(app_dir, "staging", "linked_packages")
            if self._check_local(name, item["source"], dname):
                messages.append(f"{name} content changed")

        return messages

    def uninstall_extension(self, name):
        """Uninstall an extension by name.

        Returns `True` if a rebuild is recommended, `False` otherwise.
        """
        info = self.info
        logger = self.logger

        if name in info["federated_extensions"]:
            if (
                info["federated_extensions"][name]
                .get("install", {})
                .get("uninstallInstructions", None)
            ):
                instructions = info["federated_extensions"][name]["install"][
                    "uninstallInstructions"
                ]
                logger.error(f"JupyterLab cannot uninstall this extension. {instructions}")
            else:
                logger.error(
                    f"JupyterLab cannot uninstall {name} since it was installed outside of JupyterLab. Use the same method used to install this extension to uninstall this extension."
                )
            return False

        # Allow for uninstalled core extensions.
        if name in info["core_extensions"]:
            config = self._read_build_config()
            uninstalled = config.get("uninstalled_core_extensions", [])
            if name not in uninstalled:
                logger.info(f"Uninstalling core extension {name}")
                uninstalled.append(name)
                config["uninstalled_core_extensions"] = uninstalled
                self._write_build_config(config)
                return True
            return False

        local = info["local_extensions"]

        for extname, data in info["extensions"].items():
            path = data["path"]
            if extname == name:
                msg = f"Uninstalling {name} from {osp.dirname(path)}"
                logger.info(msg)
                os.remove(path)
                # Handle local extensions.
                if extname in local:
                    config = self._read_build_config()
                    data = config.setdefault("local_extensions", {})  # noqa PLW2901
                    del data[extname]
                    self._write_build_config(config)
                return True

        logger.warning(f'No labextension named "{name}" installed')
        return False

    def uninstall_all_extensions(self):
        """Uninstalls all extensions

        Returns `True` if a rebuild is recommended, `False` otherwise
        """
        should_rebuild = False
        for extname, _ in self.info["extensions"].items():
            uninstalled = self.uninstall_extension(extname)
            should_rebuild = should_rebuild or uninstalled
        return should_rebuild

    def update_all_extensions(self):
        """Update all non-local extensions.

        Returns `True` if a rebuild is recommended, `False` otherwise.
        """
        should_rebuild = False
        for extname, _ in self.info["extensions"].items():
            if extname in self.info["local_extensions"]:
                continue
            updated = self._update_extension(extname)
            # Rebuild if at least one update happens:
            should_rebuild = should_rebuild or updated
        return should_rebuild

    def update_extension(self, name):
        """Update an extension by name.

        Returns `True` if a rebuild is recommended, `False` otherwise.
        """
        if name not in self.info["extensions"]:
            self.logger.warning(f'No labextension named "{name}" installed')
            return False
        return self._update_extension(name)

    def _update_extension(self, name):
        """Update an extension by name.

        Returns `True` if a rebuild is recommended, `False` otherwise.
        """
        data = self.info["extensions"][name]
        if data["alias_package_source"]:
            self.logger.warning(f"Skipping updating pinned extension '{name}'.")
            return False
        try:
            latest = self._latest_compatible_package_version(name)
        except URLError:
            return False
        if latest is None:
            self.logger.warning(f"No compatible version found for {name}!")
            return False
        if latest == data["version"]:
            self.logger.info(f"Extension {name!r} already up to date")
            return False
        self.logger.info(f"Updating {name} to version {latest}")
        return self.install_extension(f"{name}@{latest}")

    def link_package(self, path):
        """Link a package at the given path.

        Returns `True` if a rebuild is recommended, `False` otherwise.
        """
        path = _normalize_path(path)
        if not osp.exists(path) or not osp.isdir(path):
            msg = f'Cannot install "{path}" only link local directories'
            raise ValueError(msg)

        with TemporaryDirectory() as tempdir:
            info = self._extract_package(path, tempdir)

        messages = _validate_extension(info["data"])
        if not messages:
            return self.install_extension(path)

        # Warn that it is a linked package.
        self.logger.warning(
            f"Installing {path} as a linked package because it does not have extension metadata:"
        )
        [self.logger.warning(f"   {m}") for m in messages]

        # Add to metadata.
        config = self._read_build_config()
        linked = config.setdefault("linked_packages", {})
        linked[info["name"]] = info["source"]
        self._write_build_config(config)

        return True

    def unlink_package(self, path):
        """Unlink a package by name or at the given path.

        A ValueError is raised if the path is not an unlinkable package.

        Returns `True` if a rebuild is recommended, `False` otherwise.
        """
        path = _normalize_path(path)
        config = self._read_build_config()
        linked = config.setdefault("linked_packages", {})

        found = None
        for name, source in linked.items():
            if path in {name, source}:
                found = name

        if found:
            del linked[found]
        else:
            local = config.setdefault("local_extensions", {})
            for name, source in local.items():
                if path in {name, source}:
                    found = name
            if found:
                del local[found]
                path = self.info["extensions"][found]["path"]
                os.remove(path)

        if not found:
            msg = f"No linked package for {path}"
            raise ValueError(msg)

        self._write_build_config(config)
        return True

    def _is_extension_locked(self, extension, level="sys_prefix", include_higher_levels=True):
        app_settings_dir = osp.join(self.app_dir, "settings")
        page_config = get_static_page_config(
            app_settings_dir=app_settings_dir,
            logger=self.logger,
            level=level,
            include_higher_levels=True,
        )

        locked = page_config.get("lockedExtensions", {})
        return locked.get(extension, False)

    def toggle_extension(self, extension, value, level="sys_prefix"):
        """Enable or disable a lab extension.

        Returns `True` if a rebuild is recommended, `False` otherwise.
        """
        app_settings_dir = osp.join(self.app_dir, "settings")

        # If extension is locked at a higher level, we don't toggle it.
        # The highest level at which an extension can be locked is system,
        # so we do not need to check levels above that one.
        if level != "system":
            allowed = get_allowed_levels()
            if self._is_extension_locked(
                extension, level=allowed[allowed.index(level) + 1], include_higher_levels=True
            ):
                self.logger.info("Extension locked at a higher level, cannot toggle status")
                return False

        complete_page_config = get_static_page_config(
            app_settings_dir=app_settings_dir, logger=self.logger, level="all"
        )

        level_page_config = get_static_page_config(
            app_settings_dir=app_settings_dir, logger=self.logger, level=level
        )

        disabled = complete_page_config.get("disabledExtensions", {})
        disabled_at_level = level_page_config.get("disabledExtensions", {})
        did_something = False
        is_disabled = disabled.get(extension, False)

        if value and not is_disabled:
            disabled_at_level[extension] = True
            did_something = True
        elif not value and is_disabled:
            disabled_at_level[extension] = False
            did_something = True

        if did_something:
            level_page_config["disabledExtensions"] = disabled_at_level
            write_page_config(level_page_config, level=level)
        return did_something

    def _maybe_mirror_disabled_in_locked(self, level="sys_prefix"):
        """Lock all extensions that were previously disabled.

        This exists to facilitate migration from 4.0 (which did not include lock
        function) to 4.1 which exposes the plugin management to users in UI.

        Returns `True` if migration happened, `False` otherwise.
        """
        app_settings_dir = osp.join(self.app_dir, "settings")

        page_config = get_static_page_config(
            app_settings_dir=app_settings_dir, logger=self.logger, level=level
        )
        if "lockedExtensions" in page_config:
            # short-circuit if migration already happened
            return False

        # copy disabled onto lockedExtensions, ensuring the mapping format
        disabled = page_config.get("disabledExtensions", {})
        if isinstance(disabled, list):
            disabled = dict.fromkeys(disabled, True)
        page_config["lockedExtensions"] = disabled
        write_page_config(page_config, level=level)
        return True

    def toggle_extension_lock(self, extension, value, level="sys_prefix"):
        """Lock or unlock a lab extension (/plugin)."""
        app_settings_dir = osp.join(self.app_dir, "settings")

        # The highest level at which an extension can be locked is system,
        # so we do not need to check levels above that one.
        if level != "system":
            allowed = get_allowed_levels()
            if self._is_extension_locked(
                extension, level=allowed[allowed.index(level) + 1], include_higher_levels=True
            ):
                self.logger.info("Extension locked at a higher level, cannot toggle")
                return False

        page_config = get_static_page_config(
            app_settings_dir=app_settings_dir, logger=self.logger, level=level
        )

        locked = page_config.get("lockedExtensions", {})
        locked[extension] = value
        page_config["lockedExtensions"] = locked
        write_page_config(page_config, level=level)

    def check_extension(self, extension, check_installed_only=False):
        """Check if a lab extension is enabled or disabled"""
        self._ensure_disabled_info()
        info = self.info

        if extension in info["core_extensions"]:
            return self._check_core_extension(extension, info, check_installed_only)

        if extension in info["linked_packages"]:
            self.logger.info(f"{extension}:{GREEN_ENABLED}")
            return True

        return self._check_common_extension(extension, info, check_installed_only)

    def _check_core_extension(self, extension, info, check_installed_only):
        """Check if a core extension is enabled or disabled"""
        if extension in info["uninstalled_core"]:
            self.logger.info(f"{extension}:{RED_X}")
            return False
        if check_installed_only:
            self.logger.info(f"{extension}: {GREEN_OK}")
            return True
        if extension in info["disabled_core"]:
            self.logger.info(f"{extension}: {RED_DISABLED}")
            return False
        self.logger.info(f"{extension}:{GREEN_ENABLED}")
        return True

    def _check_common_extension(self, extension, info, check_installed_only):
        """Check if a common (non-core) extension is enabled or disabled"""
        if extension not in info["extensions"]:
            self.logger.info(f"{extension}:{RED_X}")
            return False

        errors = self._get_extension_compat()[extension]
        if errors:
            self.logger.info(f"{extension}:{RED_X} (compatibility errors)")
            return False

        if check_installed_only:
            self.logger.info(f"{extension}: {GREEN_OK}")
            return True

        if _is_disabled(extension, info["disabled"]):
            self.logger.info(f"{extension}: {RED_DISABLED}")
            return False

        self.logger.info(f"{extension}:{GREEN_ENABLED}")
        return True

    def _get_app_info(self):
        """Get information about the app."""

        info = {}
        info["core_data"] = core_data = self.core_data
        info["extensions"] = extensions = self._get_extensions(core_data)

        info["local_extensions"] = self._get_local_extensions()
        info["linked_packages"] = self._get_linked_packages()
        info["app_extensions"] = app = []
        info["sys_extensions"] = sys = []
        for name, data in extensions.items():
            data["is_local"] = name in info["local_extensions"]
            if data["location"] == "app":
                app.append(name)
            else:
                sys.append(name)

        info["uninstalled_core"] = self._get_uninstalled_core_extensions()

        info["static_data"] = _get_static_data(self.app_dir)
        app_data = info["static_data"] or core_data
        info["version"] = app_data["jupyterlab"]["version"]
        info["staticUrl"] = app_data["jupyterlab"].get("staticUrl", "")

        info["sys_dir"] = self.sys_dir
        info["app_dir"] = self.app_dir

        info["core_extensions"] = _get_core_extensions(self.core_data)

        info["federated_extensions"] = get_federated_extensions(self.labextensions_path)
        info["shadowed_exts"] = [
            ext for ext in info["extensions"] if ext in info["federated_extensions"]
        ]
        return info

    def _ensure_disabled_info(self):
        info = self.info
        if "disabled" in info:
            return

        labextensions_path = self.labextensions_path
        app_settings_dir = osp.join(self.app_dir, "settings")

        page_config = get_page_config(
            labextensions_path, app_settings_dir=app_settings_dir, logger=self.logger
        )

        disabled = page_config.get("disabledExtensions", {})
        # handle disabledExtensions specified as a list (jupyterlab_server < 2.10)
        # see https://github.com/jupyterlab/jupyterlab_server/pull/192 for more info
        if isinstance(disabled, list):
            disabled = dict.fromkeys(disabled, True)

        info["disabled"] = disabled

        locked = page_config.get("lockedExtensions", {})
        if isinstance(locked, list):
            locked = dict.fromkeys(locked, True)
        info["locked"] = locked

        disabled_core = []
        for key in info["core_extensions"]:
            if key in info["disabled"]:
                disabled_core.append(key)

        info["disabled_core"] = disabled_core

    def _populate_staging(self, name=None, version=None, static_url=None, clean=False):  # noqa
        """Set up the assets in the staging directory."""
        app_dir = self.app_dir
        staging = pjoin(app_dir, "staging")
        if clean and osp.exists(staging):
            self.logger.info(f"Cleaning {staging}")
            _rmtree(staging, self.logger)

        self._ensure_app_dirs()
        if not version:
            version = self.info["core_data"]["jupyterlab"]["version"]

        splice_source = self._options.splice_source
        if splice_source:
            self.logger.debug("Splicing dev packages into app directory.")
            source_dir = DEV_DIR
            version = __version__ + "-spliced"
        else:
            source_dir = pjoin(HERE, "staging")

        # Look for mismatched version.
        pkg_path = pjoin(staging, "package.json")

        if osp.exists(pkg_path):
            with open(pkg_path) as fid:
                data = json.load(fid)
            if data["jupyterlab"].get("version", "") != version:
                _rmtree(staging, self.logger)
                os.makedirs(staging)

        for fname in [
            "index.js",
            "bootstrap.js",
            "publicpath.js",
            "webpack.config.js",
            "webpack.prod.config.js",
            "webpack.prod.minimize.config.js",
        ]:
            target = pjoin(staging, fname)
            shutil.copy(pjoin(source_dir, fname), target)

        for fname in [".yarnrc.yml", "yarn.js"]:
            target = pjoin(staging, fname)
            shutil.copy(pjoin(HERE, "staging", fname), target)

        # Ensure a clean templates directory
        templates = pjoin(staging, "templates")
        if osp.exists(templates):
            _rmtree(templates, self.logger)

        try:
            shutil.copytree(pjoin(source_dir, "templates"), templates)
        except shutil.Error as error:
            # `copytree` throws an error if copying to + from NFS even though
            # the copy is successful (see https://bugs.python.org/issue24564
            # and https://github.com/jupyterlab/jupyterlab/issues/5233)

            real_error = "[Errno 22]" not in str(error) and "[Errno 5]" not in str(error)
            if real_error or not osp.exists(templates):
                raise

        # Ensure a clean linked packages directory.
        linked_dir = pjoin(staging, "linked_packages")
        if osp.exists(linked_dir):
            _rmtree(linked_dir, self.logger)
        os.makedirs(linked_dir)

        # Template the package.json file.
        # Update the local extensions.
        extensions = self.info["extensions"]
        removed = False
        for key, source in self.info["local_extensions"].items():
            # Handle a local extension that was removed.
            if key not in extensions:
                config = self._read_build_config()
                data = config.setdefault("local_extensions", {})
                del data[key]
                self._write_build_config(config)
                removed = True
                continue
            dname = pjoin(app_dir, "extensions")
            self._update_local(key, source, dname, extensions[key], "local_extensions")

        # Update the list of local extensions if any were removed.
        if removed:
            self.info["local_extensions"] = self._get_local_extensions()

        # Update the linked packages.
        linked = self.info["linked_packages"]
        for key, item in linked.items():
            dname = pjoin(staging, "linked_packages")
            self._update_local(key, item["source"], dname, item, "linked_packages")

        # Then get the package template.
        data = self._get_package_template()
        jlab = data["jupyterlab"]

        if version:
            jlab["version"] = version

        if name:
            jlab["name"] = name

        if static_url:
            jlab["staticUrl"] = static_url

        # Handle splicing of packages
        if splice_source:
            # Splice workspace tree as linked dependencies
            for path in glob(pjoin(REPO_ROOT, "packages", "*", "package.json")):
                local_path = osp.dirname(osp.abspath(path))
                pkg_data = json.loads(Path(path).read_text(encoding="utf-8"))
                name = pkg_data["name"]
                if name in data["dependencies"]:
                    data["dependencies"][name] = local_path
                    jlab["linkedPackages"][name] = local_path
                if name in data["resolutions"]:
                    data["resolutions"][name] = local_path

            # splice the builder as well
            local_path = osp.abspath(pjoin(REPO_ROOT, "builder"))
            data["devDependencies"]["@jupyterlab/builder"] = local_path
            target = osp.join(staging, "node_modules", "@jupyterlab", "builder")

            # Remove node_modules so it gets re-populated
            node_modules = pjoin(staging, "node_modules")
            if osp.exists(node_modules):
                shutil.rmtree(node_modules, ignore_errors=True)

        # Write the package file
        pkg_path = pjoin(staging, "package.json")
        with open(pkg_path, "w") as fid:
            json.dump(data, fid, indent=4)

        # copy known-good yarn.lock if missing
        lock_path = pjoin(staging, "yarn.lock")
        lock_template = pjoin(HERE, "staging", "yarn.lock")
        if not osp.exists(lock_path):
            shutil.copy(lock_template, lock_path)
            os.chmod(lock_path, stat.S_IWRITE | stat.S_IREAD)

    def _get_package_template(self, silent=False):  # noqa
        """Get the template the for staging package.json file."""
        logger = self.logger
        # make a deep copy of the data so we don't influence the core data
        data = deepcopy(self.info["core_data"])
        local = self.info["local_extensions"]
        linked = self.info["linked_packages"]
        extensions = self.info["extensions"]
        shadowed_exts = self.info["shadowed_exts"]
        jlab = data["jupyterlab"]

        def format_path(path):
            path = osp.relpath(path, osp.abspath(osp.realpath(pjoin(self.app_dir, "staging"))))
            path = "file:" + path.replace(os.sep, "/")
            if os.name == "nt":
                path = path.lower()
            return path

        jlab["linkedPackages"] = {}

        # Handle local extensions.
        for key, source in local.items():
            if key in shadowed_exts:
                continue
            jlab["linkedPackages"][key] = source
            data["resolutions"][key] = "file:" + self.info["extensions"][key]["path"]

        # Handle linked packages.
        for key, item in linked.items():
            if key in shadowed_exts:
                continue
            path = pjoin(self.app_dir, "staging", "linked_packages")
            path = pjoin(path, item["filename"])
            data["dependencies"][key] = format_path(path)
            jlab["linkedPackages"][key] = item["source"]
            data["resolutions"][key] = format_path(path)

        data["jupyterlab"]["extensionMetadata"] = {}

        # Handle extensions
        compat_errors = self._get_extension_compat()
        for key, value in extensions.items():
            # Reject incompatible extensions with a message.
            errors = compat_errors[key]
            if errors:
                if not silent:
                    _log_single_compat_errors(logger, key, value["version"], errors)
                continue

            data["dependencies"][key] = format_path(value["path"])

            jlab_data = value["jupyterlab"]
            for item in ["extension", "mimeExtension"]:
                ext = jlab_data.get(item, False)
                if not ext:
                    continue
                if ext is True:
                    ext = ""
                jlab[item + "s"][key] = ext

                # Add metadata for the extension
                data["jupyterlab"]["extensionMetadata"][key] = jlab_data

        # Handle uninstalled core extensions.
        for item in self.info["uninstalled_core"]:
            if item in jlab["extensions"]:
                data["jupyterlab"]["extensions"].pop(item)
            elif item in jlab["mimeExtensions"]:
                data["jupyterlab"]["mimeExtensions"].pop(item)
            # Remove from dependencies as well.
            if item in data["dependencies"]:
                data["dependencies"].pop(item)

        return data

    def _check_local(self, name, source, dname):
        """Check if a local package has changed.

        `dname` is the directory name of existing package tar archives.
        """
        # Extract the package in a temporary directory.
        with TemporaryDirectory() as tempdir:
            info = self._extract_package(source, tempdir)
            # Test if the file content has changed.
            # This relies on `_extract_package` adding the hashsum
            # to the filename, allowing a simple exist check to
            # compare the hash to the "cache" in dname.
            target = pjoin(dname, info["filename"])
            return not osp.exists(target)

    def _update_local(self, name, source, dname, data, dtype):
        """Update a local dependency.  Return `True` if changed."""
        # Extract the package in a temporary directory.
        existing = data["filename"]
        if not osp.exists(pjoin(dname, existing)):
            existing = ""

        with TemporaryDirectory() as tempdir:
            info = self._extract_package(source, tempdir)

            # Bail if the file content has not changed.
            if info["filename"] == existing:
                return existing

            shutil.move(info["path"], pjoin(dname, info["filename"]))

        # Remove the previous tarball and return the new file name.
        if existing:
            os.remove(pjoin(dname, existing))

        data["filename"] = info["filename"]
        data["path"] = pjoin(data["tar_dir"], data["filename"])
        return info["filename"]

    def _get_extensions(self, core_data):
        """Get the extensions for the application."""
        app_dir = self.app_dir
        extensions = {}

        # Get system level packages.
        sys_path = pjoin(self.sys_dir, "extensions")
        app_path = pjoin(self.app_dir, "extensions")

        extensions = self._get_extensions_in_dir(self.sys_dir, core_data)

        # Look in app_dir if different.
        app_path = pjoin(app_dir, "extensions")
        if app_path == sys_path or not osp.exists(app_path):
            return extensions

        extensions.update(self._get_extensions_in_dir(app_dir, core_data))

        return extensions

    def _get_extensions_in_dir(self, dname, core_data):
        """Get the extensions in a given directory."""
        extensions = {}
        location = "app" if dname == self.app_dir else "sys"
        for target in glob(pjoin(dname, "extensions", "*.tgz")):
            data = read_package(target)
            deps = data.get("dependencies", {})
            name = data["name"]
            jlab = data.get("jupyterlab", {})
            path = osp.abspath(target)

            filename = osp.basename(target)
            if filename.startswith(PIN_PREFIX):
                alias = filename[len(PIN_PREFIX) : -len(".tgz")]
            else:
                alias = None
            url = get_package_url(data)
            extensions[alias or name] = {
                "description": data.get("description", ""),
                "path": path,
                "filename": osp.basename(path),
                "url": url,
                "version": data["version"],
                # Only save the package name if the extension name is an alias
                "alias_package_source": name if alias else None,
                "jupyterlab": jlab,
                "dependencies": deps,
                "tar_dir": osp.dirname(path),
                "location": location,
            }
        return extensions

    def _get_extension_compat(self):
        """Get the extension compatibility info."""
        compat = {}
        core_data = self.info["core_data"]
        seen = set()
        for name, data in self.info["federated_extensions"].items():
            deps = data["dependencies"]
            compat[name] = _validate_compatibility(name, deps, core_data)
            seen.add(name)
        for name, data in self.info["extensions"].items():
            if name in seen:
                continue
            deps = data["dependencies"]
            compat[name] = _validate_compatibility(name, deps, core_data)
        return compat

    def _get_local_extensions(self):
        """Get the locally installed extensions."""
        return self._get_local_data("local_extensions")

    def _get_linked_packages(self):
        """Get the linked packages."""
        info = self._get_local_data("linked_packages")
        dname = pjoin(self.app_dir, "staging", "linked_packages")
        for name, source in info.items():
            info[name] = {"source": source, "filename": "", "tar_dir": dname}

        if not osp.exists(dname):
            return info

        for path in glob(pjoin(dname, "*.tgz")):
            path = osp.abspath(path)  # noqa PLW2901
            data = read_package(path)
            name = data["name"]
            if name not in info:
                self.logger.warning(f"Removing orphaned linked package {name}")
                os.remove(path)
                continue
            item = info[name]
            item["filename"] = osp.basename(path)
            item["path"] = path
            item["version"] = data["version"]
            item["data"] = data
        return info

    def _get_uninstalled_core_extensions(self):
        """Get the uninstalled core extensions."""
        config = self._read_build_config()
        return config.get("uninstalled_core_extensions", [])

    def _ensure_app_dirs(self):
        """Ensure that the application directories exist"""
        dirs = ["extensions", "settings", "staging", "schemas", "themes"]
        for dname in dirs:
            path = pjoin(self.app_dir, dname)
            if not osp.exists(path):
                try:
                    os.makedirs(path)
                except OSError as e:
                    if e.errno != errno.EEXIST:
                        raise

    def _list_extensions(self, info, ext_type):
        """List the extensions of a given type."""
        self._ensure_disabled_info()
        logger = self.logger
        names = info[f"{ext_type}_extensions"]
        if not names:
            return

        dname = info[f"{ext_type}_dir"]

        error_accumulator = {}

        logger.info(f"   {ext_type} dir: {dname}")
        for name in sorted(names):
            if name in info["federated_extensions"]:
                continue
            data = info["extensions"][name]
            version = data["version"]
            errors = info["compat_errors"][name]
            extra = self._compose_extra_status(name, info, data, errors)

            # If we have the package name in the data, this means this extension's name is the alias name
            alias_package_source = data["alias_package_source"]
            if alias_package_source:
                logger.info(f"        {name} {alias_package_source} v{version}{extra}")
            else:
                logger.info(f"        {name} v{version}{extra}")
            if errors:
                error_accumulator[name] = (version, errors)

        # Write all errors at end:
        _log_multiple_compat_errors(logger, error_accumulator, self.verbose)

        # Write a blank line separator
        logger.info("")

    def _list_federated_extensions(self):
        self._ensure_disabled_info()
        info = self.info
        logger = self.logger

        error_accumulator = {}

        ext_dirs = dict.fromkeys(self.labextensions_path, False)
        for value in info["federated_extensions"].values():
            ext_dirs[value["ext_dir"]] = True

        for ext_dir, has_exts in ext_dirs.items():
            if not has_exts:
                continue
            logger.info(ext_dir)
            for name in info["federated_extensions"]:
                data = info["federated_extensions"][name]
                if data["ext_dir"] != ext_dir:
                    continue
                version = data["version"]
                errors = info["compat_errors"][name]
                extra = self._compose_extra_status(name, info, data, errors)

                install = data.get("install")
                if install:
                    extra += " ({}, {})".format(install["packageManager"], install["packageName"])
                logger.info(f"        {name} v{version}{extra}")
                if errors:
                    error_accumulator[name] = (version, errors)
            # Add a spacer line after
            logger.info("")

        # Write all errors at end:
        _log_multiple_compat_errors(logger, error_accumulator, self.verbose)

    def _compose_extra_status(self, name: str, info: dict, data: dict, errors) -> str:
        extra = ""
        if _is_disabled(name, info["disabled"]):
            extra += f" {RED_DISABLED}"
        else:
            extra += f" {GREEN_ENABLED}"
        if errors:
            extra += f" {RED_X}"
        else:
            extra += f" {GREEN_OK}"
        if data["is_local"]:
            extra += "*"
        lock_status = _is_locked(name, info["locked"])
        if lock_status.entire_extension_locked:
            extra += " 🔒 (all plugins locked)"
        elif lock_status.locked_plugins:
            plugin_list = ", ".join(sorted(lock_status.locked_plugins))
            extra += f" 🔒 (plugins: {plugin_list} locked)"
        return extra

    def _read_build_config(self):
        """Get the build config data for the app dir."""
        target = pjoin(self.app_dir, "settings", "build_config.json")
        if not osp.exists(target):
            return {}
        else:
            with open(target) as fid:
                return json.load(fid)

    def _write_build_config(self, config):
        """Write the build config to the app dir."""
        self._ensure_app_dirs()
        target = pjoin(self.app_dir, "settings", "build_config.json")
        with open(target, "w") as fid:
            json.dump(config, fid, indent=4)

    def _get_local_data(self, source):
        """Get the local data for extensions or linked packages."""
        config = self._read_build_config()

        data = config.setdefault(source, {})
        dead = []
        for name, source_ in data.items():
            if not osp.exists(source_):
                dead.append(name)

        for name in dead:
            link_type = source.replace("_", " ")
            msg = f'**Note: Removing dead {link_type} "{name}"'
            self.logger.warning(msg)
            del data[name]

        if dead:
            self._write_build_config(config)

        return data

    def _install_extension(self, extension, tempdir, pin=None):
        """Install an extension with validation and return the name and path."""
        info = self._extract_package(extension, tempdir, pin=pin)
        data = info["data"]

        # Check for compatible version unless:
        # - A specific version was requested (@ in name,
        #   but after first char to allow for scope marker).
        # - Package is locally installed.
        allow_fallback = "@" not in extension[1:] and not info["is_dir"]
        name = info["name"]

        # Verify that the package is an extension.
        messages = _validate_extension(data)
        if messages:
            all_messages = "\n".join(messages)
            msg = f'"{extension}" is not a valid extension:\n{all_messages}'
            if allow_fallback:
                try:
                    version = self._latest_compatible_package_version(name)
                except URLError:
                    raise ValueError(msg) from None
            else:
                raise ValueError(msg)

        # Verify package compatibility.
        deps = data.get("dependencies", {})
        errors = _validate_compatibility(extension, deps, self.core_data)
        if errors:
            msg = _format_compatibility_errors(data["name"], data["version"], errors)
            if allow_fallback:
                try:
                    version = self._latest_compatible_package_version(name)
                except URLError:
                    # We cannot add any additional information to error message
                    raise ValueError(msg) from None

                if version and name:
                    self.logger.debug(f"Incompatible extension:\n{name}")
                    self.logger.debug(f"Found compatible version: {version}")
                    with TemporaryDirectory() as tempdir2:
                        return self._install_extension(f"{name}@{version}", tempdir2)

                # Extend message to better guide the user what to do:
                conflicts = "\n".join(msg.splitlines()[2:])
                msg = "".join((self._format_no_compatible_package_version(name), "\n\n", conflicts))

            raise ValueError(msg)

        # Move the file to the app directory.
        target = pjoin(self.app_dir, "extensions", info["filename"])
        if osp.exists(target):
            os.remove(target)

        shutil.move(info["path"], target)

        info["path"] = target
        return info

    def _extract_package(self, source, tempdir, pin=None):
        """Call `npm pack` for an extension.

        The pack command will download the package tar if `source` is
        a package name, or run `npm pack` locally if `source` is a
        directory.
        """
        is_dir = osp.exists(source) and osp.isdir(source)
        if is_dir and not osp.exists(pjoin(source, "node_modules")):
            self._run(["node", YARN_PATH, "install"], cwd=source)

        info = {"source": source, "is_dir": is_dir}

        ret = self._run([which("npm"), "pack", source], cwd=tempdir)
        if ret != 0:
            msg = f'"{source}" is not a valid npm package'
            raise ValueError(msg)

        path = glob(pjoin(tempdir, "*.tgz"))[0]
        info["data"] = read_package(path)
        if is_dir:
            info["sha"] = sha = _tarsum(path)
            target = path.replace(".tgz", f"-{sha}.tgz")
            shutil.move(path, target)
            info["path"] = target
        else:
            info["path"] = path
        if pin:
            old_path = info["path"]
            new_path = pjoin(osp.dirname(old_path), f"{PIN_PREFIX}{pin}.tgz")
            shutil.move(old_path, new_path)
            info["path"] = new_path

        info["filename"] = osp.basename(info["path"])
        info["name"] = info["data"]["name"]
        info["version"] = info["data"]["version"]

        return info

    def _latest_compatible_package_version(self, name):
        """Get the latest compatible version of a package"""
        core_data = self.info["core_data"]
        try:
            metadata = _fetch_package_metadata(self.registry, name, self.logger)
        except URLError:
            return
        versions = metadata.get("versions", {})

        # Sort pre-release first, as we will reverse the sort:
        def sort_key(key_value):
            return _semver_key(key_value[0], prerelease_first=True)

        for version, data in sorted(versions.items(), key=sort_key, reverse=True):
            deps = data.get("dependencies", {})
            errors = _validate_compatibility(name, deps, core_data)
            if not errors:
                # Found a compatible version
                # skip deprecated versions
                if "deprecated" in data:
                    self.logger.debug(
                        f"Disregarding compatible version of package as it is deprecated: {name}@{version}"
                    )
                    continue
                # Verify that the version is a valid extension.
                with TemporaryDirectory() as tempdir:
                    info = self._extract_package(f"{name}@{version}", tempdir)
                if _validate_extension(info["data"]):
                    # Invalid, do not consider other versions
                    return
                # Valid
                return version

    def latest_compatible_package_versions(self, names):
        """Get the latest compatible versions of several packages

        Like _latest_compatible_package_version, but optimized for
        retrieving the latest version for several packages in one go.
        """
        core_data = self.info["core_data"]

        keys = []
        for name in names:
            try:
                metadata = _fetch_package_metadata(self.registry, name, self.logger)
            except URLError:
                continue
            versions = metadata.get("versions", {})

            # Sort pre-release first, as we will reverse the sort:
            def sort_key(key_value):
                return _semver_key(key_value[0], prerelease_first=True)

            for version, data in sorted(versions.items(), key=sort_key, reverse=True):
                # skip deprecated versions
                if "deprecated" in data:
                    continue

                deps = data.get("dependencies", {})
                errors = _validate_compatibility(name, deps, core_data)
                if not errors:
                    # Found a compatible version
                    keys.append(f"{name}@{version}")
                    break  # break inner for

        versions = {}
        if not keys:
            return versions
        with TemporaryDirectory() as tempdir:
            ret = self._run([which("npm"), "pack", *keys], cwd=tempdir)
            if ret != 0:
                msg = f'"{keys}" is not a valid npm package'
                raise ValueError(msg)

            for key in keys:
                fname = (
                    key[0].replace("@", "") + key[1:].replace("@", "-").replace("/", "-") + ".tgz"
                )
                data = read_package(osp.join(tempdir, fname))
                # Verify that the version is a valid extension.
                if not _validate_extension(data):
                    # Valid
                    versions[data["name"]] = data["version"]
        return versions

    def _format_no_compatible_package_version(self, name):
        """Get the latest compatible version of a package"""
        core_data = self.info["core_data"]
        # Whether lab version is too new:
        lab_newer_than_latest = False
        # Whether the latest version of the extension depend on a "future" version
        # of a singleton package (from the perspective of current lab version):
        latest_newer_than_lab = False
        try:
            metadata = _fetch_package_metadata(self.registry, name, self.logger)
        except URLError:
            pass
        else:
            versions = metadata.get("versions", {})

            # Sort pre-release first, as we will reverse the sort:
            def sort_key(key_value):
                return _semver_key(key_value[0], prerelease_first=True)

            store = tuple(sorted(versions.items(), key=sort_key, reverse=True))
            latest_deps = store[0][1].get("dependencies", {})
            core_deps = core_data["resolutions"]
            singletons = core_data["jupyterlab"]["singletonPackages"]

            for key, value in latest_deps.items():
                if key in singletons:
                    # Drop prereleases in comparisons to allow extension authors
                    # to not have to update their versions for each
                    # Jupyterlab prerelease version.
                    c = _compare_ranges(core_deps[key], value, drop_prerelease1=True)
                    lab_newer_than_latest = lab_newer_than_latest or c < 0
                    latest_newer_than_lab = latest_newer_than_lab or c > 0

        if lab_newer_than_latest:
            # All singleton deps in current version of lab are newer than those
            # in the latest version of the extension
            return (
                f'The extension "{name}" does not yet support the current version of JupyterLab.\n'
            )

        parts = [
            f"No version of {name} could be found that is compatible with "
            "the current version of JupyterLab."
        ]
        if latest_newer_than_lab:
            parts.extend(
                (
                    "However, it seems to support a new version of JupyterLab.",
                    "Consider upgrading JupyterLab.",
                )
            )

        return " ".join(parts)

    def _run(self, cmd, **kwargs):
        """Run the command using our logger and abort callback.

        Returns the exit code.
        """
        if self.kill_event.is_set():
            msg = "Command was killed"
            raise ValueError(msg)

        kwargs["logger"] = self.logger
        kwargs["kill_event"] = self.kill_event
        proc = ProgressProcess(cmd, **kwargs)
        return proc.wait()


def _node_check(logger):
    """Check for the existence of nodejs with the correct version."""
    node = which("node")
    try:
        output = subprocess.check_output([node, "node-version-check.js"], cwd=HERE)  # noqa S603
        logger.debug(output.decode("utf-8"))
    except Exception:
        data = CoreConfig()._data
        ver = data["engines"]["node"]
        msg = f"Please install nodejs {ver} before continuing. nodejs may be installed using conda or directly from the nodejs website."
        raise ValueError(msg) from None


def _yarn_config(logger):
    """Get the yarn configuration.

    Returns
    -------
    {"yarn config": dict, "npm config": dict}
    if unsuccessful, the subdictionaries are empty
    """
    configuration = {"yarn config": {}, "npm config": {}}
    try:
        node = which("node")
    except ValueError:  # Node not found == user with no need for building jupyterlab
        logger.debug("NodeJS was not found. Yarn user configuration is ignored.")
        return configuration

    try:
        output_binary = subprocess.check_output(  # noqa S603
            [node, YARN_PATH, "config", "--json"],
            stderr=subprocess.PIPE,
            cwd=HERE,
        )
        output = output_binary.decode("utf-8")
        lines = iter(output.splitlines())
        try:
            for line in lines:
                info = json.loads(line)
                if info["type"] == "info":
                    key = info["data"]
                    inspect = json.loads(next(lines))
                    if inspect["type"] == "inspect":
                        configuration[key] = inspect["data"]
        except StopIteration:
            pass
        logger.debug("Yarn configuration loaded.")
    except subprocess.CalledProcessError as e:
        logger.error(
            "Fail to get yarn configuration. {!s}{!s}".format(
                e.stderr.decode("utf-8"), e.output.decode("utf-8")
            )
        )
    except Exception as e:
        logger.error(f"Fail to get yarn configuration. {e!s}")

    return configuration


def _ensure_logger(logger=None):
    """Ensure that we have a logger"""
    return logger or logging.getLogger("jupyterlab")


def _normalize_path(extension):
    """Normalize a given extension if it is a path."""
    extension = osp.expanduser(extension)
    if osp.exists(extension):
        extension = osp.abspath(extension)
    return extension


def _rmtree(path, logger):
    """Remove a tree, logging errors"""

    def onerror(*exc_info):
        logger.debug("Error in shutil.rmtree", exc_info=exc_info)

    shutil.rmtree(path, onerror=onerror)


def _unlink(path, logger):
    """Remove a file, logging errors"""
    try:
        os.unlink(path)
    except Exception:
        logger.debug("Error in os.unlink", exc_info=sys.exc_info())


def _rmtree_star(path, logger):
    """Remove all files/trees within a dir, logging errors"""
    for filename in os.listdir(path):
        file_path = osp.join(path, filename)
        if osp.isfile(file_path) or osp.islink(file_path):
            _unlink(file_path, logger)
        elif osp.isdir(file_path):
            _rmtree(file_path, logger)


def _validate_extension(data):  # noqa
    """Detect if a package is an extension using its metadata.

    Returns any problems it finds.
    """
    jlab = data.get("jupyterlab", None)
    if jlab is None:
        return ["No `jupyterlab` key"]
    if not isinstance(jlab, dict):
        return ["The `jupyterlab` key must be a JSON object"]
    extension = jlab.get("extension", False)
    mime_extension = jlab.get("mimeExtension", False)
    theme_path = jlab.get("themePath", "")
    schema_dir = jlab.get("schemaDir", "")

    messages = []
    if not extension and not mime_extension:
        messages.append("No `extension` or `mimeExtension` key present")

    if extension == mime_extension:
        msg = "`mimeExtension` and `extension` must point to different modules"
        messages.append(msg)

    files = data["jupyterlab_extracted_files"]
    main = data.get("main", "index.js")
    if not main.endswith(".js"):
        main += ".js"

    if extension is True:
        extension = main
    elif extension and not extension.endswith(".js"):
        extension += ".js"

    if mime_extension is True:
        mime_extension = main
    elif mime_extension and not mime_extension.endswith(".js"):
        mime_extension += ".js"

    if extension and extension not in files:
        messages.append(f'Missing extension module "{extension}"')

    if mime_extension and mime_extension not in files:
        messages.append(f'Missing mimeExtension module "{mime_extension}"')

    if theme_path and not any(f.startswith(str(Path(theme_path))) for f in files):
        messages.append(f'themePath is empty: "{theme_path}"')

    if schema_dir and not any(f.startswith(str(Path(schema_dir))) for f in files):
        messages.append(f'schemaDir is empty: "{schema_dir}"')

    return messages


def _tarsum(input_file):
    """
    Compute the recursive sha sum of a tar file.
    """
    chunk_size = 100 * 1024
    h = hashlib.new("sha1")  # noqa: S324

    with tarfile.open(input_file, "r") as tar:
        for member in tar:
            if not member.isfile():
                continue
            with tar.extractfile(member) as f:
                if f:  # Check if f is not None (safety check)
                    data = f.read(chunk_size)
                    while data:
                        h.update(data)
                        data = f.read(chunk_size)

    return h.hexdigest()


def _get_static_data(app_dir):
    """Get the data for the app static dir."""
    target = pjoin(app_dir, "static", "package.json")
    if osp.exists(target):
        with open(target) as fid:
            return json.load(fid)
    else:
        return None


def _validate_compatibility(extension, deps, core_data):
    """Validate the compatibility of an extension."""
    core_deps = core_data["resolutions"]
    singletons = core_data["jupyterlab"]["singletonPackages"]

    errors = []

    for key, value in deps.items():
        if key in singletons:
            # Drop prereleases in comparisons to allow extension authors
            # to not have to update their versions for each
            # Jupyterlab prerelease version.
            overlap = _test_overlap(core_deps[key], value, drop_prerelease1=True)
            if overlap is False:
                errors.append((key, core_deps[key], value))

    return errors


def _test_overlap(spec1, spec2, drop_prerelease1=False, drop_prerelease2=False):
    """Test whether two version specs overlap.

    Returns `None` if we cannot determine compatibility,
    otherwise whether there is an overlap
    """
    cmp = _compare_ranges(
        spec1, spec2, drop_prerelease1=drop_prerelease1, drop_prerelease2=drop_prerelease2
    )
    if cmp is None:
        return
    return cmp == 0


def _compare_ranges(spec1, spec2, drop_prerelease1=False, drop_prerelease2=False):  # noqa
    """Test whether two version specs overlap.

    Returns `None` if we cannot determine compatibility,
    otherwise return 0 if there is an overlap, 1 if
    spec1 is lower/older than spec2, and -1 if spec1
    is higher/newer than spec2.
    """
    # Test for overlapping semver ranges.
    r1 = Range(spec1, True)
    r2 = Range(spec2, True)

    # If either range is empty, we cannot verify.
    if not r1.range or not r2.range:
        return

    # Set return_value to a sentinel value
    return_value = False

    # r1.set may be a list of ranges if the range involved an ||, so we need to test for overlaps between each pair.
    for r1set, r2set in itertools.product(r1.set, r2.set):
        x1 = r1set[0].semver
        x2 = r1set[-1].semver
        y1 = r2set[0].semver
        y2 = r2set[-1].semver

        if x1.prerelease and drop_prerelease1:
            x1 = x1.inc("patch")

        if y1.prerelease and drop_prerelease2:
            y1 = y1.inc("patch")

        o1 = r1set[0].operator
        o2 = r2set[0].operator

        # We do not handle (<) specifiers.
        if o1.startswith("<") or o2.startswith("<"):
            continue

        # Handle single value specifiers.
        lx = lte if x1 == x2 else lt
        ly = lte if y1 == y2 else lt
        gx = gte if x1 == x2 else gt
        gy = gte if x1 == x2 else gt

        # Handle unbounded (>) specifiers.
        def noop(x, y, z):
            return True

        if x1 == x2 and o1.startswith(">"):
            lx = noop
        if y1 == y2 and o2.startswith(">"):
            ly = noop

        # Check for overlap.
        if (
            (gte(x1, y1, True) and ly(x1, y2, True))
            or (gy(x2, y1, True) and ly(x2, y2, True))
            or (gte(y1, x1, True) and lx(y1, x2, True))
            or (gx(y2, x1, True) and lx(y2, x2, True))
        ):
            # if we ever find an overlap, we can return immediately
            return 0

        if gte(y1, x2, True):
            if return_value is False:
                # We can possibly return 1
                return_value = 1
            elif return_value == -1:
                # conflicting information, so we must return None
                return_value = None
            continue

        if gte(x1, y2, True):
            if return_value is False:
                return_value = -1
            elif return_value == 1:
                # conflicting information, so we must return None
                return_value = None
            continue

        msg = "Unexpected case comparing version ranges"
        raise AssertionError(msg)

    if return_value is False:
        return_value = None
    return return_value


def _is_disabled(name, disabled=None):
    """Test whether the package is disabled."""
    disabled = disabled or {}
    for pattern, value in disabled.items():
        # skip packages explicitly marked as not disabled
        if value is False:
            continue
        if name == pattern:
            return True
        if re.compile(pattern).match(name) is not None:
            return True
    return False


@dataclass(frozen=True)
class LockStatus:
    entire_extension_locked: bool
    # locked plugins are only given if extension is not locked as a whole
    locked_plugins: Optional[frozenset[str]] = None


def _is_locked(name, locked=None) -> LockStatus:
    """Test whether the package is locked.

    If only a subset of extension plugins is locked return them.
    """
    locked = locked or {}
    locked_plugins = set()
    for lock, value in locked.items():
        # skip packages explicitly marked as not locked
        if value is False:
            continue
        if name == lock:
            return LockStatus(entire_extension_locked=True)
        extension_part = lock.partition(":")[0]
        if name == extension_part:
            locked_plugins.add(lock)

    return LockStatus(entire_extension_locked=False, locked_plugins=locked_plugins)


def _format_compatibility_errors(name, version, errors):
    """Format a message for compatibility errors."""
    msgs = []
    l0 = 10
    l1 = 10
    for error in errors:
        pkg, jlab, ext = error
        jlab = str(Range(jlab, True))
        ext = str(Range(ext, True))
        msgs.append((pkg, jlab, ext))
        l0 = max(l0, len(pkg) + 1)
        l1 = max(l1, len(jlab) + 1)

    msg = f'\n"{name}@{version}" is not compatible with the current JupyterLab'
    msg += "\nConflicting Dependencies:\n"
    msg += "JupyterLab".ljust(l0)
    msg += "Extension".ljust(l1)
    msg += "Package\n"

    for pkg, jlab, ext in msgs:
        msg += jlab.ljust(l0) + ext.ljust(l1) + pkg + "\n"

    return msg


def _log_multiple_compat_errors(logger, errors_map, verbose: bool):
    """Log compatibility errors for multiple extensions at once"""

    outdated = []

    for name, (_, errors) in errors_map.items():
        age = _compat_error_age(errors)
        if age > 0:
            outdated.append(name)

    if outdated:
        logger.warning(
            "\n        ".join(
                [
                    "\n   The following extensions may be outdated or specify dependencies that are incompatible with the current version of jupyterlab:",
                    *outdated,
                    "\n   If you are a user, check if an update is available for these packages.\n"
                    + (
                        "   If you are a developer, re-run with `--verbose` flag for more details.\n"
                        if not verbose
                        else "   See below for the details.\n"
                    ),
                ]
            )
        )

    # Print out compatibility errors for all extensions, even the ones inferred
    # to be possibly outdated, to guide developers upgrading their extensions.
    for name, (version, errors) in errors_map.items():
        if name in outdated and not verbose:
            continue
        msg = _format_compatibility_errors(name, version, errors)
        logger.warning(f"{msg}\n")


def _log_single_compat_errors(logger, name, version, errors):
    """Log compatibility errors for a single extension"""

    age = _compat_error_age(errors)
    if age > 0:
        logger.warning(f'The extension "{name}" is outdated.\n')
    else:
        msg = _format_compatibility_errors(name, version, errors)
        logger.warning(f"{msg}\n")


def _compat_error_age(errors):
    """Compare all incompatibilities for an extension.

    Returns a number > 0 if all extensions are older than that supported by lab.
    Returns a number < 0 if all extensions are newer than that supported by lab.
    Returns 0 otherwise (i.e. a mix).
    """
    # Do any extensions depend on too old lab packages?
    any_older = False
    # Do any extensions depend on too new lab packages?
    any_newer = False

    for _, jlab, ext in errors:
        # Drop prereleases in comparisons to allow extension authors
        # to not have to update their versions for each
        # Jupyterlab prerelease version.
        c = _compare_ranges(ext, jlab, drop_prerelease1=True)
        any_newer = any_newer or c < 0
        any_older = any_older or c > 0
    if any_older and not any_newer:
        return 1
    elif any_newer and not any_older:
        return -1
    return 0


def _get_core_extensions(core_data):
    """Get the core extensions."""
    data = core_data["jupyterlab"]
    return list(data["extensions"]) + list(data["mimeExtensions"])


def _semver_prerelease_key(prerelease):
    """Sort key for prereleases.

    Precedence for two pre-release versions with the same
    major, minor, and patch version MUST be determined by
    comparing each dot separated identifier from left to
    right until a difference is found as follows:
    identifiers consisting of only digits are compare
    numerically and identifiers with letters or hyphens
    are compared lexically in ASCII sort order. Numeric
    identifiers always have lower precedence than non-
    numeric identifiers. A larger set of pre-release
    fields has a higher precedence than a smaller set,
    if all of the preceding identifiers are equal.
    """
    for entry in prerelease:
        if isinstance(entry, int):
            # Assure numerics always sort before string
            yield ("", entry)
        else:
            # Use ASCII compare:
            yield (entry,)


def _semver_key(version, prerelease_first=False):
    """A sort key-function for sorting semver version string.

    The default sorting order is ascending (0.x -> 1.x -> 2.x).

    If `prerelease_first`, pre-releases will come before
    ALL other semver keys (not just those with same version).
    I.e (1.0-pre, 2.0-pre -> 0.x -> 1.x -> 2.x).

    Otherwise it will sort in the standard way that it simply
    comes before any release with shared version string
    (0.x -> 1.0-pre -> 1.x -> 2.0-pre -> 2.x).
    """
    v = make_semver(version, True)
    key = ((0,) if v.prerelease else (1,)) if prerelease_first else ()
    key = (*key, v.major, v.minor, v.patch)
    if not prerelease_first:
        #  NOT having a prerelease is > having one
        key = (*key, 0) if v.prerelease else (1,)
    if v.prerelease:
        key = key + tuple(_semver_prerelease_key(v.prerelease))

    return key


def _fetch_package_metadata(registry, name, logger):
    """Fetch the metadata for a package from the npm registry"""
    req = Request(  # noqa S310
        urljoin(registry, quote(name, safe="@")),
        headers={
            "Accept": ("application/vnd.npm.install-v1+json; q=1.0, application/json; q=0.8, */*")
        },
    )
    try:
        logger.debug(f"Fetching URL: {req.full_url}")
    except AttributeError:
        logger.debug(f"Fetching URL: {req.get_full_url()}")
    try:
        with contextlib.closing(urlopen(req)) as response:  # noqa S310
            return json.loads(response.read().decode("utf-8"))
    except URLError as exc:
        logger.warning(f"Failed to fetch package metadata for {name!r}: {exc!r}")
        raise


if __name__ == "__main__":
    watch_dev(HERE)
