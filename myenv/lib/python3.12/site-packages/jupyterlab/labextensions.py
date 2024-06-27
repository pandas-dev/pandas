"""Jupyter LabExtension Entry Points."""

# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.

import os
import sys
from copy import copy

from jupyter_core.application import JupyterApp, base_aliases, base_flags
from traitlets import Bool, Instance, List, Unicode, default

from jupyterlab.coreconfig import CoreConfig
from jupyterlab.debuglog import DebugLogFileMixin

from .commands import (
    HERE,
    AppOptions,
    build,
    check_extension,
    disable_extension,
    enable_extension,
    get_app_version,
    install_extension,
    link_package,
    list_extensions,
    lock_extension,
    uninstall_extension,
    unlink_package,
    unlock_extension,
    update_extension,
)
from .federated_labextensions import build_labextension, develop_labextension_py, watch_labextension
from .labapp import LabApp

flags = dict(base_flags)
flags["no-build"] = (
    {"BaseExtensionApp": {"should_build": False}},
    "Defer building the app after the action.",
)
flags["dev-build"] = (
    {"BaseExtensionApp": {"dev_build": True}},
    "Build in development mode.",
)
flags["no-minimize"] = (
    {"BaseExtensionApp": {"minimize": False}},
    "Do not minimize a production build.",
)
flags["clean"] = (
    {"BaseExtensionApp": {"should_clean": True}},
    "Cleanup intermediate files after the action.",
)
flags["splice-source"] = (
    {"BaseExtensionApp": {"splice_source": True}},
    "Splice source packages into app directory.",
)

check_flags = copy(flags)
check_flags["installed"] = (
    {"CheckLabExtensionsApp": {"should_check_installed_only": True}},
    "Check only if the extension is installed.",
)

develop_flags = copy(flags)
develop_flags["overwrite"] = (
    {"DevelopLabExtensionApp": {"overwrite": True}},
    "Overwrite files",
)

update_flags = copy(flags)
update_flags["all"] = (
    {"UpdateLabExtensionApp": {"all": True}},
    "Update all extensions",
)

uninstall_flags = copy(flags)
uninstall_flags["all"] = (
    {"UninstallLabExtensionApp": {"all": True}},
    "Uninstall all extensions",
)

list_flags = copy(flags)
list_flags["verbose"] = (
    {"ListLabExtensionsApp": {"verbose": True}},
    "Increase verbosity level",
)

aliases = dict(base_aliases)
aliases["app-dir"] = "BaseExtensionApp.app_dir"
aliases["dev-build"] = "BaseExtensionApp.dev_build"
aliases["minimize"] = "BaseExtensionApp.minimize"
aliases["debug-log-path"] = "DebugLogFileMixin.debug_log_path"

install_aliases = copy(aliases)
install_aliases["pin-version-as"] = "InstallLabExtensionApp.pin"

enable_aliases = copy(aliases)
enable_aliases["level"] = "EnableLabExtensionsApp.level"

disable_aliases = copy(aliases)
disable_aliases["level"] = "DisableLabExtensionsApp.level"

lock_aliases = copy(aliases)
lock_aliases["level"] = "LockLabExtensionsApp.level"

unlock_aliases = copy(aliases)
unlock_aliases["level"] = "UnlockLabExtensionsApp.level"

VERSION = get_app_version()

LABEXTENSION_COMMAND_WARNING = "Users should manage prebuilt extensions with package managers like pip and conda, and extension authors are encouraged to distribute their extensions as prebuilt packages"


class BaseExtensionApp(JupyterApp, DebugLogFileMixin):
    version = VERSION
    flags = flags
    aliases = aliases
    name = "lab"

    # Not configurable!
    core_config = Instance(CoreConfig, allow_none=True)

    app_dir = Unicode("", config=True, help="The app directory to target")

    should_build = Bool(True, config=True, help="Whether to build the app after the action")

    dev_build = Bool(
        None,
        allow_none=True,
        config=True,
        help="Whether to build in dev mode. Defaults to True (dev mode) if there are any locally linked extensions, else defaults to False (production mode).",
    )

    minimize = Bool(
        True,
        config=True,
        help="Whether to minimize a production build (defaults to True).",
    )

    should_clean = Bool(
        False,
        config=True,
        help="Whether temporary files should be cleaned up after building jupyterlab",
    )

    splice_source = Bool(False, config=True, help="Splice source packages into app directory.")

    labextensions_path = List(
        Unicode(),
        help="The standard paths to look in for prebuilt JupyterLab extensions",
    )

    @default("labextensions_path")
    def _default_labextensions_path(self):
        lab = LabApp()
        lab.load_config_file()
        return lab.extra_labextensions_path + lab.labextensions_path

    @default("splice_source")
    def _default_splice_source(self):
        version = get_app_version(AppOptions(app_dir=self.app_dir))
        return version.endswith("-spliced")

    def start(self):
        if self.app_dir and self.app_dir.startswith(HERE):
            msg = "Cannot run lab extension commands in core app"
            raise ValueError(msg)
        with self.debug_logging():
            ans = self.run_task()
            if ans and self.should_build:
                production = None if self.dev_build is None else not self.dev_build
                app_options = AppOptions(
                    app_dir=self.app_dir,
                    logger=self.log,
                    core_config=self.core_config,
                    splice_source=self.splice_source,
                )
                build(
                    clean_staging=self.should_clean,
                    production=production,
                    minimize=self.minimize,
                    app_options=app_options,
                )

    def run_task(self):
        pass

    def deprecation_warning(self, msg):
        return self.log.warning(
            "\033[33m(Deprecated) %s\n\n%s \033[0m", msg, LABEXTENSION_COMMAND_WARNING
        )

    def _log_format_default(self):
        """A default format for messages"""
        return "%(message)s"


class InstallLabExtensionApp(BaseExtensionApp):
    description = """Install labextension(s)

     Usage

        jupyter labextension install [--pin-version-as <alias,...>] <package...>

    This installs JupyterLab extensions similar to yarn add or npm install.

    Pass a list of comma separate names to the --pin-version-as flag
    to use as aliases for the packages providers. This is useful to
    install multiple versions of the same extension.
    These can be uninstalled with the alias you provided
    to the flag, similar to the "alias" feature of yarn add.
    """
    aliases = install_aliases

    pin = Unicode("", config=True, help="Pin this version with a certain alias")

    def run_task(self):
        self.deprecation_warning(
            "Installing extensions with the jupyter labextension install command is now deprecated and will be removed in a future major version of JupyterLab."
        )
        pinned_versions = self.pin.split(",")
        self.extra_args = self.extra_args or [os.getcwd()]
        return any(
            install_extension(
                arg,
                # Pass in pinned alias if we have it
                pin=pinned_versions[i] if i < len(pinned_versions) else None,
                app_options=AppOptions(
                    app_dir=self.app_dir,
                    logger=self.log,
                    core_config=self.core_config,
                    labextensions_path=self.labextensions_path,
                ),
            )
            for i, arg in enumerate(self.extra_args)
        )


class DevelopLabExtensionApp(BaseExtensionApp):
    description = "(developer) Develop labextension"
    flags = develop_flags

    user = Bool(False, config=True, help="Whether to do a user install")
    sys_prefix = Bool(True, config=True, help="Use the sys.prefix as the prefix")
    overwrite = Bool(False, config=True, help="Whether to overwrite files")
    symlink = Bool(True, config=False, help="Whether to use a symlink")

    labextensions_dir = Unicode(
        "",
        config=True,
        help="Full path to labextensions dir (probably use prefix or user)",
    )

    def run_task(self):
        """Add config for this labextension"""
        self.extra_args = self.extra_args or [os.getcwd()]
        for arg in self.extra_args:
            develop_labextension_py(
                arg,
                user=self.user,
                sys_prefix=self.sys_prefix,
                labextensions_dir=self.labextensions_dir,
                logger=self.log,
                overwrite=self.overwrite,
                symlink=self.symlink,
            )


class BuildLabExtensionApp(BaseExtensionApp):
    description = "(developer) Build labextension"

    static_url = Unicode("", config=True, help="Sets the url for static assets when building")

    development = Bool(False, config=True, help="Build in development mode")

    source_map = Bool(False, config=True, help="Generate source maps")

    core_path = Unicode(
        os.path.join(HERE, "staging"),
        config=True,
        help="Directory containing core application package.json file",
    )

    aliases = {
        "static-url": "BuildLabExtensionApp.static_url",
        "development": "BuildLabExtensionApp.development",
        "source-map": "BuildLabExtensionApp.source_map",
        "core-path": "BuildLabExtensionApp.core_path",
    }

    def run_task(self):
        self.extra_args = self.extra_args or [os.getcwd()]
        build_labextension(
            self.extra_args[0],
            logger=self.log,
            development=self.development,
            static_url=self.static_url or None,
            source_map=self.source_map,
            core_path=self.core_path or None,
        )


class WatchLabExtensionApp(BaseExtensionApp):
    description = "(developer) Watch labextension"

    development = Bool(True, config=True, help="Build in development mode")

    source_map = Bool(False, config=True, help="Generate source maps")

    core_path = Unicode(
        os.path.join(HERE, "staging"),
        config=True,
        help="Directory containing core application package.json file",
    )

    aliases = {
        "core-path": "WatchLabExtensionApp.core_path",
        "development": "WatchLabExtensionApp.development",
        "source-map": "WatchLabExtensionApp.source_map",
    }

    def run_task(self):
        self.extra_args = self.extra_args or [os.getcwd()]
        labextensions_path = self.labextensions_path
        watch_labextension(
            self.extra_args[0],
            labextensions_path,
            logger=self.log,
            development=self.development,
            source_map=self.source_map,
            core_path=self.core_path or None,
        )


class UpdateLabExtensionApp(BaseExtensionApp):
    description = "Update labextension(s)"
    flags = update_flags

    all = Bool(False, config=True, help="Whether to update all extensions")

    def run_task(self):
        self.deprecation_warning(
            "Updating extensions with the jupyter labextension update command is now deprecated and will be removed in a future major version of JupyterLab."
        )
        if not self.all and not self.extra_args:
            self.log.warning(
                "Specify an extension to update, or use --all to update all extensions"
            )
            return False
        app_options = AppOptions(
            app_dir=self.app_dir,
            logger=self.log,
            core_config=self.core_config,
            labextensions_path=self.labextensions_path,
        )
        if self.all:
            return update_extension(all_=True, app_options=app_options)
        return any(update_extension(name=arg, app_options=app_options) for arg in self.extra_args)


class LinkLabExtensionApp(BaseExtensionApp):
    description = """
    Link local npm packages that are not lab extensions.

    Links a package to the JupyterLab build process. A linked
    package is manually re-installed from its source location when
    `jupyter lab build` is run.
    """
    should_build = Bool(True, config=True, help="Whether to build the app after the action")

    def run_task(self):
        self.extra_args = self.extra_args or [os.getcwd()]
        options = AppOptions(
            app_dir=self.app_dir,
            logger=self.log,
            labextensions_path=self.labextensions_path,
            core_config=self.core_config,
        )
        return any(link_package(arg, app_options=options) for arg in self.extra_args)


class UnlinkLabExtensionApp(BaseExtensionApp):
    description = "Unlink packages by name or path"

    def run_task(self):
        self.extra_args = self.extra_args or [os.getcwd()]
        options = AppOptions(
            app_dir=self.app_dir,
            logger=self.log,
            labextensions_path=self.labextensions_path,
            core_config=self.core_config,
        )
        return any(unlink_package(arg, app_options=options) for arg in self.extra_args)


class UninstallLabExtensionApp(BaseExtensionApp):
    description = "Uninstall labextension(s) by name"
    flags = uninstall_flags

    all = Bool(False, config=True, help="Whether to uninstall all extensions")

    def run_task(self):
        self.deprecation_warning(
            "Uninstalling extensions with the jupyter labextension uninstall command is now deprecated and will be removed in a future major version of JupyterLab."
        )
        self.extra_args = self.extra_args or [os.getcwd()]

        options = AppOptions(
            app_dir=self.app_dir,
            logger=self.log,
            labextensions_path=self.labextensions_path,
            core_config=self.core_config,
        )
        return any(
            uninstall_extension(arg, all_=self.all, app_options=options) for arg in self.extra_args
        )


class ListLabExtensionsApp(BaseExtensionApp):
    description = "List the installed labextensions"
    verbose = Bool(False, help="Increase verbosity level.").tag(config=True)
    flags = list_flags

    def run_task(self):
        list_extensions(
            app_options=AppOptions(
                app_dir=self.app_dir,
                logger=self.log,
                core_config=self.core_config,
                labextensions_path=self.labextensions_path,
                verbose=self.verbose,
            )
        )


class EnableLabExtensionsApp(BaseExtensionApp):
    description = "Enable labextension(s) by name"
    aliases = enable_aliases

    level = Unicode("sys_prefix", help="Level at which to enable: sys_prefix, user, system").tag(
        config=True
    )

    def run_task(self):
        app_options = AppOptions(
            app_dir=self.app_dir,
            logger=self.log,
            core_config=self.core_config,
            labextensions_path=self.labextensions_path,
        )
        [
            enable_extension(arg, app_options=app_options, level=self.level)
            for arg in self.extra_args
        ]


class DisableLabExtensionsApp(BaseExtensionApp):
    description = "Disable labextension(s) by name"
    aliases = disable_aliases

    level = Unicode("sys_prefix", help="Level at which to disable: sys_prefix, user, system").tag(
        config=True
    )

    def run_task(self):
        app_options = AppOptions(
            app_dir=self.app_dir,
            logger=self.log,
            core_config=self.core_config,
            labextensions_path=self.labextensions_path,
        )
        [
            disable_extension(arg, app_options=app_options, level=self.level)
            for arg in self.extra_args
        ]
        self.log.info(
            "Starting with JupyterLab 4.1 individual plugins can be re-enabled"
            " in the user interface. While all plugins which were previously"
            " disabled have been locked, you need to explicitly lock any newly"
            " disabled plugins by using `jupyter labextension lock` command."
        )


class LockLabExtensionsApp(BaseExtensionApp):
    description = "Lock labextension(s) by name"
    aliases = lock_aliases

    level = Unicode("sys_prefix", help="Level at which to lock: sys_prefix, user, system").tag(
        config=True
    )

    def run_task(self):
        app_options = AppOptions(
            app_dir=self.app_dir,
            logger=self.log,
            core_config=self.core_config,
            labextensions_path=self.labextensions_path,
        )
        [lock_extension(arg, app_options=app_options, level=self.level) for arg in self.extra_args]


class UnlockLabExtensionsApp(BaseExtensionApp):
    description = "Unlock labextension(s) by name"
    aliases = unlock_aliases

    level = Unicode("sys_prefix", help="Level at which to unlock: sys_prefix, user, system").tag(
        config=True
    )

    def run_task(self):
        app_options = AppOptions(
            app_dir=self.app_dir,
            logger=self.log,
            core_config=self.core_config,
            labextensions_path=self.labextensions_path,
        )
        [
            unlock_extension(arg, app_options=app_options, level=self.level)
            for arg in self.extra_args
        ]


class CheckLabExtensionsApp(BaseExtensionApp):
    description = "Check labextension(s) by name"
    flags = check_flags

    should_check_installed_only = Bool(
        False,
        config=True,
        help="Whether it should check only if the extensions is installed",
    )

    def run_task(self):
        app_options = AppOptions(
            app_dir=self.app_dir,
            logger=self.log,
            core_config=self.core_config,
            labextensions_path=self.labextensions_path,
        )
        all_enabled = all(
            check_extension(
                arg, installed=self.should_check_installed_only, app_options=app_options
            )
            for arg in self.extra_args
        )
        if not all_enabled:
            self.exit(1)


_EXAMPLES = """
jupyter labextension list                        # list all configured labextensions
jupyter labextension install <extension name>    # install a labextension
jupyter labextension uninstall <extension name>  # uninstall a labextension
jupyter labextension develop                     # (developer) develop a prebuilt labextension
jupyter labextension build                       # (developer) build a prebuilt labextension
jupyter labextension watch                       # (developer) watch a prebuilt labextension
"""


class LabExtensionApp(JupyterApp):
    """Base jupyter labextension command entry point"""

    name = "jupyter labextension"
    version = VERSION
    description = "Work with JupyterLab extensions"
    examples = _EXAMPLES

    subcommands = {
        "install": (InstallLabExtensionApp, "Install labextension(s)"),
        "update": (UpdateLabExtensionApp, "Update labextension(s)"),
        "uninstall": (UninstallLabExtensionApp, "Uninstall labextension(s)"),
        "list": (ListLabExtensionsApp, "List labextensions"),
        "link": (LinkLabExtensionApp, "Link labextension(s)"),
        "unlink": (UnlinkLabExtensionApp, "Unlink labextension(s)"),
        "enable": (EnableLabExtensionsApp, "Enable labextension(s)"),
        "disable": (DisableLabExtensionsApp, "Disable labextension(s)"),
        "lock": (LockLabExtensionsApp, "Lock labextension(s)"),
        "unlock": (UnlockLabExtensionsApp, "Unlock labextension(s)"),
        "check": (CheckLabExtensionsApp, "Check labextension(s)"),
        "develop": (DevelopLabExtensionApp, "(developer) Develop labextension(s)"),
        "build": (BuildLabExtensionApp, "(developer) Build labextension"),
        "watch": (WatchLabExtensionApp, "(developer) Watch labextension"),
    }

    def start(self):
        """Perform the App's functions as configured"""
        super().start()

        # The above should have called a subcommand and raised NoStart; if we
        # get here, it didn't, so we should self.log.info a message.
        subcmds = ", ".join(sorted(self.subcommands))
        self.exit("Please supply at least one subcommand: %s" % subcmds)


main = LabExtensionApp.launch_instance

if __name__ == "__main__":
    sys.exit(main())
