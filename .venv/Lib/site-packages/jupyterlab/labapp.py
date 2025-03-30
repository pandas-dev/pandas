"""A tornado based Jupyter lab server."""

# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.

import dataclasses
import json
import os
import sys

from jupyter_core.application import JupyterApp, NoStart, base_aliases, base_flags
from jupyter_server._version import version_info as jpserver_version_info
from jupyter_server.serverapp import flags
from jupyter_server.utils import url_path_join as ujoin
from jupyterlab_server import (
    LabServerApp,
    LicensesApp,
    WorkspaceExportApp,
    WorkspaceImportApp,
    WorkspaceListApp,
)
from notebook_shim.shim import NotebookConfigShimMixin
from traitlets import Bool, Instance, Type, Unicode, default

from ._version import __version__
from .commands import (
    DEV_DIR,
    HERE,
    AppOptions,
    build,
    clean,
    ensure_app,
    ensure_core,
    ensure_dev,
    get_app_dir,
    get_app_version,
    get_user_settings_dir,
    get_workspaces_dir,
    pjoin,
    watch,
    watch_dev,
)
from .coreconfig import CoreConfig
from .debuglog import DebugLogFileMixin
from .extensions import MANAGERS as EXT_MANAGERS
from .extensions.manager import PluginManager
from .extensions.readonly import ReadOnlyExtensionManager
from .handlers.announcements import (
    CheckForUpdate,
    CheckForUpdateABC,
    CheckForUpdateHandler,
    NewsHandler,
    check_update_handler_path,
    news_handler_path,
)
from .handlers.build_handler import Builder, BuildHandler, build_path
from .handlers.error_handler import ErrorHandler
from .handlers.extension_manager_handler import ExtensionHandler, extensions_handler_path
from .handlers.plugin_manager_handler import PluginHandler, plugins_handler_path

DEV_NOTE = """You're running JupyterLab from source.
If you're working on the TypeScript sources of JupyterLab, try running

    jupyter lab --dev-mode --watch


to have the system incrementally watch and build JupyterLab for you, as you
make changes.
"""


CORE_NOTE = """
Running the core application with no additional extensions or settings
"""

build_aliases = dict(base_aliases)
build_aliases["app-dir"] = "LabBuildApp.app_dir"
build_aliases["name"] = "LabBuildApp.name"
build_aliases["version"] = "LabBuildApp.version"
build_aliases["dev-build"] = "LabBuildApp.dev_build"
build_aliases["minimize"] = "LabBuildApp.minimize"
build_aliases["debug-log-path"] = "DebugLogFileMixin.debug_log_path"

build_flags = dict(base_flags)

build_flags["dev-build"] = (
    {"LabBuildApp": {"dev_build": True}},
    "Build in development mode.",
)
build_flags["no-minimize"] = (
    {"LabBuildApp": {"minimize": False}},
    "Do not minimize a production build.",
)
build_flags["splice-source"] = (
    {"LabBuildApp": {"splice_source": True}},
    "Splice source packages into app directory.",
)


version = __version__
app_version = get_app_version()
if version != app_version:
    version = f"{__version__} (dev), {app_version} (app)"

build_failure_msg = """Build failed.
Troubleshooting: If the build failed due to an out-of-memory error, you
may be able to fix it by disabling the `dev_build` and/or `minimize` options.

If you are building via the `jupyter lab build` command, you can disable
these options like so:

jupyter lab build --dev-build=False --minimize=False

You can also disable these options for all JupyterLab builds by adding these
lines to a Jupyter config file named `jupyter_config.py`:

c.LabBuildApp.minimize = False
c.LabBuildApp.dev_build = False

If you don't already have a `jupyter_config.py` file, you can create one by
adding a blank file of that name to any of the Jupyter config directories.
The config directories can be listed by running:

jupyter --paths

Explanation:

- `dev-build`: This option controls whether a `dev` or a more streamlined
`production` build is used. This option will default to `False` (i.e., the
`production` build) for most users. However, if you have any labextensions
installed from local files, this option will instead default to `True`.
Explicitly setting `dev-build` to `False` will ensure that the `production`
build is used in all circumstances.

- `minimize`: This option controls whether your JS bundle is minified
during the Webpack build, which helps to improve JupyterLab's overall
performance. However, the minifier plugin used by Webpack is very memory
intensive, so turning it off may help the build finish successfully in
low-memory environments.
"""


class LabBuildApp(JupyterApp, DebugLogFileMixin):
    version = version
    description = """
    Build the JupyterLab application

    The application is built in the JupyterLab app directory in `/staging`.
    When the build is complete it is put in the JupyterLab app `/static`
    directory, where it is used to serve the application.
    """
    aliases = build_aliases
    flags = build_flags

    # Not configurable!
    core_config = Instance(CoreConfig, allow_none=True)

    app_dir = Unicode("", config=True, help="The app directory to build in")

    name = Unicode("JupyterLab", config=True, help="The name of the built application")

    version = Unicode("", config=True, help="The version of the built application")

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

    pre_clean = Bool(
        False, config=True, help="Whether to clean before building (defaults to False)"
    )

    splice_source = Bool(False, config=True, help="Splice source packages into app directory.")

    def start(self):
        app_dir = self.app_dir or get_app_dir()
        app_options = AppOptions(
            app_dir=app_dir,
            logger=self.log,
            core_config=self.core_config,
            splice_source=self.splice_source,
        )
        self.log.info(f"JupyterLab {version}")
        with self.debug_logging():
            if self.pre_clean:
                self.log.info(f"Cleaning {app_dir}")
                clean(app_options=app_options)
            self.log.info(f"Building in {app_dir}")
            try:
                production = None if self.dev_build is None else not self.dev_build
                build(
                    name=self.name,
                    version=self.version,
                    app_options=app_options,
                    production=production,
                    minimize=self.minimize,
                )
            except Exception as e:
                self.log.error(build_failure_msg)
                raise e


clean_aliases = dict(base_aliases)
clean_aliases["app-dir"] = "LabCleanApp.app_dir"

ext_warn_msg = "WARNING: this will delete all of your extensions, which will need to be reinstalled"

clean_flags = dict(base_flags)
clean_flags["extensions"] = (
    {"LabCleanApp": {"extensions": True}},
    f"Also delete <app-dir>/extensions.\n{ext_warn_msg}",
)
clean_flags["settings"] = (
    {"LabCleanApp": {"settings": True}},
    "Also delete <app-dir>/settings",
)
clean_flags["static"] = (
    {"LabCleanApp": {"static": True}},
    "Also delete <app-dir>/static",
)
clean_flags["all"] = (
    {"LabCleanApp": {"all": True}},
    f"Delete the entire contents of the app directory.\n{ext_warn_msg}",
)


class LabCleanAppOptions(AppOptions):
    extensions = Bool(False)
    settings = Bool(False)
    staging = Bool(True)
    static = Bool(False)
    all = Bool(False)


class LabCleanApp(JupyterApp):
    version = version
    description = """
    Clean the JupyterLab application

    This will clean the app directory by removing the `staging` directories.
    Optionally, the `extensions`, `settings`, and/or `static` directories,
    or the entire contents of the app directory, can also be removed.
    """
    aliases = clean_aliases
    flags = clean_flags

    # Not configurable!
    core_config = Instance(CoreConfig, allow_none=True)

    app_dir = Unicode("", config=True, help="The app directory to clean")

    extensions = Bool(False, config=True, help=f"Also delete <app-dir>/extensions.\n{ext_warn_msg}")

    settings = Bool(False, config=True, help="Also delete <app-dir>/settings")

    static = Bool(False, config=True, help="Also delete <app-dir>/static")

    all = Bool(
        False,
        config=True,
        help=f"Delete the entire contents of the app directory.\n{ext_warn_msg}",
    )

    def start(self):
        app_options = LabCleanAppOptions(
            logger=self.log,
            core_config=self.core_config,
            app_dir=self.app_dir,
            extensions=self.extensions,
            settings=self.settings,
            static=self.static,
            all=self.all,
        )
        clean(app_options=app_options)


class LabPathApp(JupyterApp):
    version = version
    description = """
    Print the configured paths for the JupyterLab application

    The application path can be configured using the JUPYTERLAB_DIR
        environment variable.
    The user settings path can be configured using the JUPYTERLAB_SETTINGS_DIR
        environment variable or it will fall back to
        `/lab/user-settings` in the default Jupyter configuration directory.
    The workspaces path can be configured using the JUPYTERLAB_WORKSPACES_DIR
        environment variable or it will fall back to
        '/lab/workspaces' in the default Jupyter configuration directory.
    """

    def start(self):
        print(f"Application directory:   {get_app_dir()}")
        print(f"User Settings directory: {get_user_settings_dir()}")
        print(f"Workspaces directory: {get_workspaces_dir()}")


class LabWorkspaceExportApp(WorkspaceExportApp):
    version = version

    @default("workspaces_dir")
    def _default_workspaces_dir(self):
        return get_workspaces_dir()


class LabWorkspaceImportApp(WorkspaceImportApp):
    version = version

    @default("workspaces_dir")
    def _default_workspaces_dir(self):
        return get_workspaces_dir()


class LabWorkspaceListApp(WorkspaceListApp):
    version = version

    @default("workspaces_dir")
    def _default_workspaces_dir(self):
        return get_workspaces_dir()


class LabWorkspaceApp(JupyterApp):
    version = version
    description = """
    Import or export a JupyterLab workspace or list all the JupyterLab workspaces

    There are three sub-commands for export, import or listing of workspaces. This app
        should not otherwise do any work.
    """
    subcommands = {}
    subcommands["export"] = (
        LabWorkspaceExportApp,
        LabWorkspaceExportApp.description.splitlines()[0],
    )
    subcommands["import"] = (
        LabWorkspaceImportApp,
        LabWorkspaceImportApp.description.splitlines()[0],
    )
    subcommands["list"] = (
        LabWorkspaceListApp,
        LabWorkspaceListApp.description.splitlines()[0],
    )

    def start(self):
        try:
            super().start()
            self.log.error("One of `export`, `import` or `list` must be specified.")
            self.exit(1)
        except NoStart:
            pass
        self.exit(0)


class LabLicensesApp(LicensesApp):
    version = version

    dev_mode = Bool(
        False,
        config=True,
        help="""Whether to start the app in dev mode. Uses the unpublished local
        JavaScript packages in the `dev_mode` folder.  In this case JupyterLab will
        show a red stripe at the top of the page.  It can only be used if JupyterLab
        is installed as `pip install -e .`.
        """,
    )

    app_dir = Unicode("", config=True, help="The app directory for which to show licenses")

    aliases = {
        **LicensesApp.aliases,
        "app-dir": "LabLicensesApp.app_dir",
    }

    flags = {
        **LicensesApp.flags,
        "dev-mode": (
            {"LabLicensesApp": {"dev_mode": True}},
            "Start the app in dev mode for running from source.",
        ),
    }

    @default("app_dir")
    def _default_app_dir(self):
        return get_app_dir()

    @default("static_dir")
    def _default_static_dir(self):
        return pjoin(self.app_dir, "static")


aliases = dict(base_aliases)
aliases.update(
    {
        "ip": "ServerApp.ip",
        "port": "ServerApp.port",
        "port-retries": "ServerApp.port_retries",
        "keyfile": "ServerApp.keyfile",
        "certfile": "ServerApp.certfile",
        "client-ca": "ServerApp.client_ca",
        "notebook-dir": "ServerApp.root_dir",
        "browser": "ServerApp.browser",
        "pylab": "ServerApp.pylab",
    }
)


class LabApp(NotebookConfigShimMixin, LabServerApp):
    version = version

    name = "lab"
    app_name = "JupyterLab"

    # Should your extension expose other server extensions when launched directly?
    load_other_extensions = True

    description = """
    JupyterLab - An extensible computational environment for Jupyter.

    This launches a Tornado based HTML Server that serves up an
    HTML5/Javascript JupyterLab client.

    JupyterLab has three different modes of running:

    * Core mode (`--core-mode`): in this mode JupyterLab will run using the JavaScript
      assets contained in the installed `jupyterlab` Python package. In core mode, no
      extensions are enabled. This is the default in a stable JupyterLab release if you
      have no extensions installed.
    * Dev mode (`--dev-mode`): uses the unpublished local JavaScript packages in the
      `dev_mode` folder.  In this case JupyterLab will show a red stripe at the top of
      the page.  It can only be used if JupyterLab is installed as `pip install -e .`.
    * App mode: JupyterLab allows multiple JupyterLab "applications" to be
      created by the user with different combinations of extensions. The `--app-dir` can
      be used to set a directory for different applications. The default application
      path can be found using `jupyter lab path`.
    """

    examples = """
        jupyter lab                       # start JupyterLab
        jupyter lab --dev-mode            # start JupyterLab in development mode, with no extensions
        jupyter lab --core-mode           # start JupyterLab in core mode, with no extensions
        jupyter lab --app-dir=~/myjupyterlabapp # start JupyterLab with a particular set of extensions
        jupyter lab --certfile=mycert.pem # use SSL/TLS certificate
    """

    aliases = aliases
    aliases.update(
        {
            "watch": "LabApp.watch",
        }
    )
    aliases["app-dir"] = "LabApp.app_dir"

    flags = flags
    flags["core-mode"] = (
        {"LabApp": {"core_mode": True}},
        "Start the app in core mode.",
    )
    flags["dev-mode"] = (
        {"LabApp": {"dev_mode": True}},
        "Start the app in dev mode for running from source.",
    )
    flags["skip-dev-build"] = (
        {"LabApp": {"skip_dev_build": True}},
        "Skip the initial install and JS build of the app in dev mode.",
    )
    flags["watch"] = ({"LabApp": {"watch": True}}, "Start the app in watch mode.")
    flags["splice-source"] = (
        {"LabApp": {"splice_source": True}},
        "Splice source packages into app directory.",
    )
    flags["expose-app-in-browser"] = (
        {"LabApp": {"expose_app_in_browser": True}},
        "Expose the global app instance to browser via window.jupyterapp.",
    )
    flags["extensions-in-dev-mode"] = (
        {"LabApp": {"extensions_in_dev_mode": True}},
        "Load prebuilt extensions in dev-mode.",
    )
    flags["collaborative"] = (
        {"LabApp": {"collaborative": True}},
        """To enable real-time collaboration, you must install the extension `jupyter_collaboration`.
        You can install it using pip for example:

            python -m pip install jupyter_collaboration

        This flag is now deprecated and will be removed in JupyterLab v5.""",
    )
    flags["custom-css"] = (
        {"LabApp": {"custom_css": True}},
        "Load custom CSS in template html files. Default is False",
    )

    subcommands = {
        "build": (LabBuildApp, LabBuildApp.description.splitlines()[0]),
        "clean": (LabCleanApp, LabCleanApp.description.splitlines()[0]),
        "path": (LabPathApp, LabPathApp.description.splitlines()[0]),
        "paths": (LabPathApp, LabPathApp.description.splitlines()[0]),
        "workspace": (LabWorkspaceApp, LabWorkspaceApp.description.splitlines()[0]),
        "workspaces": (LabWorkspaceApp, LabWorkspaceApp.description.splitlines()[0]),
        "licenses": (LabLicensesApp, LabLicensesApp.description.splitlines()[0]),
    }

    default_url = Unicode("/lab", config=True, help="The default URL to redirect to from `/`")

    override_static_url = Unicode(
        config=True, help=("The override url for static lab assets, typically a CDN.")
    )

    override_theme_url = Unicode(
        config=True,
        help=("The override url for static lab theme assets, typically a CDN."),
    )

    app_dir = Unicode(None, config=True, help="The app directory to launch JupyterLab from.")

    user_settings_dir = Unicode(
        get_user_settings_dir(), config=True, help="The directory for user settings."
    )

    workspaces_dir = Unicode(get_workspaces_dir(), config=True, help="The directory for workspaces")

    core_mode = Bool(
        False,
        config=True,
        help="""Whether to start the app in core mode. In this mode, JupyterLab
        will run using the JavaScript assets that are within the installed
        JupyterLab Python package. In core mode, third party extensions are disabled.
        The `--dev-mode` flag is an alias to this to be used when the Python package
        itself is installed in development mode (`pip install -e .`).
        """,
    )

    dev_mode = Bool(
        False,
        config=True,
        help="""Whether to start the app in dev mode. Uses the unpublished local
        JavaScript packages in the `dev_mode` folder.  In this case JupyterLab will
        show a red stripe at the top of the page.  It can only be used if JupyterLab
        is installed as `pip install -e .`.
        """,
    )

    extensions_in_dev_mode = Bool(
        False,
        config=True,
        help="""Whether to load prebuilt extensions in dev mode. This may be
        useful to run and test prebuilt extensions in development installs of
        JupyterLab. APIs in a JupyterLab development install may be
        incompatible with published packages, so prebuilt extensions compiled
        against published packages may not work correctly.""",
    )

    extension_manager = Unicode(
        "pypi",
        config=True,
        help="""The extension manager factory to use. The default options are:
        "readonly" for a manager without installation capability or "pypi" for
        a manager using PyPi.org and pip to install extensions.""",
    )

    watch = Bool(False, config=True, help="Whether to serve the app in watch mode")

    skip_dev_build = Bool(
        False,
        config=True,
        help="Whether to skip the initial install and JS build of the app in dev mode",
    )

    splice_source = Bool(False, config=True, help="Splice source packages into app directory.")

    expose_app_in_browser = Bool(
        False,
        config=True,
        help="Whether to expose the global app instance to browser via window.jupyterapp",
    )

    custom_css = Bool(
        False,
        config=True,
        help="""Whether custom CSS is loaded on the page.
    Defaults to False.
    """,
    )

    collaborative = Bool(
        False,
        config=True,
        help="""To enable real-time collaboration, you must install the extension `jupyter_collaboration`.
        You can install it using pip for example:

            python -m pip install jupyter_collaboration

        This flag is now deprecated and will be removed in JupyterLab v5.""",
    )

    news_url = Unicode(
        "https://jupyterlab.github.io/assets/feed.xml",
        allow_none=True,
        help="""URL that serves news Atom feed; by default the JupyterLab organization announcements will be fetched. Set to None to turn off fetching announcements.""",
        config=True,
    )

    lock_all_plugins = Bool(
        False,
        config=True,
        help="Whether all plugins are locked (cannot be enabled/disabled from the UI)",
    )

    check_for_updates_class = Type(
        default_value=CheckForUpdate,
        klass=CheckForUpdateABC,
        config=True,
        help="""A callable class that receives the current version at instantiation and calling it must return asynchronously a string indicating which version is available and how to install or None if no update is available. The string supports Markdown format.""",
    )

    @default("app_dir")
    def _default_app_dir(self):
        app_dir = get_app_dir()
        if self.core_mode:
            app_dir = HERE
        elif self.dev_mode:
            app_dir = DEV_DIR
        return app_dir

    @default("app_settings_dir")
    def _default_app_settings_dir(self):
        return pjoin(self.app_dir, "settings")

    @default("app_version")
    def _default_app_version(self):
        return app_version

    @default("cache_files")
    def _default_cache_files(self):
        return False

    @default("schemas_dir")
    def _default_schemas_dir(self):
        return pjoin(self.app_dir, "schemas")

    @default("templates_dir")
    def _default_templates_dir(self):
        return pjoin(self.app_dir, "static")

    @default("themes_dir")
    def _default_themes_dir(self):
        if self.override_theme_url:
            return ""
        return pjoin(self.app_dir, "themes")

    @default("static_dir")
    def _default_static_dir(self):
        return pjoin(self.app_dir, "static")

    @default("static_url_prefix")
    def _default_static_url_prefix(self):
        if self.override_static_url:
            return self.override_static_url
        else:
            static_url = f"/static/{self.name}/"
            return ujoin(self.serverapp.base_url, static_url)

    @default("theme_url")
    def _default_theme_url(self):
        if self.override_theme_url:
            return self.override_theme_url
        return ""

    def initialize_templates(self):
        # Determine which model to run JupyterLab
        if self.core_mode or self.app_dir.startswith(HERE + os.sep):
            self.core_mode = True
            self.log.info("Running JupyterLab in core mode")

        if self.dev_mode or self.app_dir.startswith(DEV_DIR + os.sep):
            self.dev_mode = True
            self.log.info("Running JupyterLab in dev mode")

        if self.watch and self.core_mode:
            self.log.warning("Cannot watch in core mode, did you mean --dev-mode?")
            self.watch = False

        if self.core_mode and self.dev_mode:
            self.log.warning("Conflicting modes, choosing dev_mode over core_mode")
            self.core_mode = False

        # Set the paths based on JupyterLab's mode.
        if self.dev_mode:
            dev_static_dir = ujoin(DEV_DIR, "static")
            self.static_paths = [dev_static_dir]
            self.template_paths = [dev_static_dir]
            if not self.extensions_in_dev_mode:
                # Add an exception for @jupyterlab/galata-extension
                galata_extension = pjoin(HERE, "galata")
                self.labextensions_path = (
                    [galata_extension]
                    if galata_extension in map(os.path.abspath, self.labextensions_path)
                    else []
                )
                self.extra_labextensions_path = (
                    [galata_extension]
                    if galata_extension in map(os.path.abspath, self.extra_labextensions_path)
                    else []
                )
        elif self.core_mode:
            dev_static_dir = ujoin(HERE, "static")
            self.static_paths = [dev_static_dir]
            self.template_paths = [dev_static_dir]
            self.labextensions_path = []
            self.extra_labextensions_path = []
        else:
            self.static_paths = [self.static_dir]
            self.template_paths = [self.templates_dir]

    def _prepare_templates(self):
        super()._prepare_templates()
        self.jinja2_env.globals.update(custom_css=self.custom_css)

    def initialize_handlers(self):  # noqa
        handlers = []

        # Set config for Jupyterlab
        page_config = self.serverapp.web_app.settings.setdefault("page_config_data", {})
        page_config.setdefault("buildAvailable", not self.core_mode and not self.dev_mode)
        page_config.setdefault("buildCheck", not self.core_mode and not self.dev_mode)
        page_config["devMode"] = self.dev_mode
        page_config["token"] = self.serverapp.identity_provider.token
        page_config["exposeAppInBrowser"] = self.expose_app_in_browser
        page_config["quitButton"] = self.serverapp.quit_button
        page_config["allow_hidden_files"] = self.serverapp.contents_manager.allow_hidden

        # Client-side code assumes notebookVersion is a JSON-encoded string
        page_config["notebookVersion"] = json.dumps(jpserver_version_info)

        self.log.info(f"JupyterLab extension loaded from {HERE!s}")
        self.log.info(f"JupyterLab application directory is {self.app_dir!s}")

        if self.custom_css:
            handlers.append(
                (
                    r"/custom/(.*)(?<!\.js)$",
                    self.serverapp.web_app.settings["static_handler_class"],
                    {
                        "path": self.serverapp.web_app.settings["static_custom_path"],
                        "no_cache_paths": ["/"],  # don't cache anything in custom
                    },
                )
            )

        app_options = AppOptions(
            logger=self.log,
            app_dir=self.app_dir,
            labextensions_path=self.extra_labextensions_path + self.labextensions_path,
            splice_source=self.splice_source,
        )
        builder = Builder(self.core_mode, app_options=app_options)
        build_handler = (build_path, BuildHandler, {"builder": builder})
        handlers.append(build_handler)

        errored = False

        if self.core_mode:
            self.log.info(CORE_NOTE.strip())
            ensure_core(self.log)
        elif self.dev_mode:
            if not (self.watch or self.skip_dev_build):
                ensure_dev(self.log)
                self.log.info(DEV_NOTE)
        else:
            if self.splice_source:
                ensure_dev(self.log)
            msgs = ensure_app(self.app_dir)
            if msgs:
                [self.log.error(msg) for msg in msgs]
                handler = (self.app_url, ErrorHandler, {"messages": msgs})
                handlers.append(handler)
                errored = True

        if self.watch:
            self.log.info("Starting JupyterLab watch mode...")
            if self.dev_mode:
                watch_dev(self.log)
            else:
                watch(app_options=app_options)
                page_config["buildAvailable"] = False
            self.cache_files = False

        if not self.core_mode and not errored:
            # Add extension management handlers
            provider = self.extension_manager
            entry_point = EXT_MANAGERS.get(provider)
            if entry_point is None:
                self.log.error(f"Extension Manager: No manager defined for provider '{provider}'.")
                raise NotImplementedError()
            else:
                self.log.info(f"Extension Manager is '{provider}'.")
            manager_factory = entry_point.load()
            config = self.settings.get("config", {}).get("LabServerApp", {})

            blocked_extensions_uris = config.get("blocked_extensions_uris", "")
            allowed_extensions_uris = config.get("allowed_extensions_uris", "")

            if (blocked_extensions_uris) and (allowed_extensions_uris):
                self.log.error(
                    "Simultaneous LabServerApp.blocked_extensions_uris and LabServerApp.allowed_extensions_uris is not supported. Please define only one of those."
                )
                import sys

                sys.exit(-1)

            listings_config = {
                "blocked_extensions_uris": set(
                    filter(lambda uri: len(uri) > 0, blocked_extensions_uris.split(","))
                ),
                "allowed_extensions_uris": set(
                    filter(lambda uri: len(uri) > 0, allowed_extensions_uris.split(","))
                ),
                "listings_refresh_seconds": config.get("listings_refresh_seconds", 60 * 60),
                "listings_tornado_options": config.get("listings_tornado_options", {}),
            }
            if len(listings_config["blocked_extensions_uris"]) or len(
                listings_config["allowed_extensions_uris"]
            ):
                self.log.debug(f"Extension manager will be constrained by {listings_config}")

            try:
                ext_manager = manager_factory(app_options, listings_config, self)
                metadata = dataclasses.asdict(ext_manager.metadata)
            except Exception as err:
                self.log.warning(
                    f"Failed to instantiate the extension manager {provider}. Falling back to read-only manager.",
                    exc_info=err,
                )
                ext_manager = ReadOnlyExtensionManager(app_options, listings_config, self)
                metadata = dataclasses.asdict(ext_manager.metadata)

            page_config["extensionManager"] = metadata
            ext_handler = (
                extensions_handler_path,
                ExtensionHandler,
                {"manager": ext_manager},
            )
            handlers.append(ext_handler)

            # Add plugin manager handlers
            lock_rules = frozenset(
                {rule for rule, value in page_config.get("lockedExtensions", {}).items() if value}
            )
            handlers.append(
                (
                    plugins_handler_path,
                    PluginHandler,
                    {
                        "manager": PluginManager(
                            app_options=app_options,
                            ext_options={
                                "lock_rules": lock_rules,
                                "all_locked": self.lock_all_plugins,
                            },
                            parent=self,
                        )
                    },
                )
            )

            # Add announcement handlers
            page_config["news"] = {"disabled": self.news_url is None}
            handlers.extend(
                [
                    (
                        check_update_handler_path,
                        CheckForUpdateHandler,
                        {
                            "update_checker": self.check_for_updates_class(__version__),
                        },
                    ),
                    (
                        news_handler_path,
                        NewsHandler,
                        {
                            "news_url": self.news_url,
                        },
                    ),
                ]
            )

        # If running under JupyterHub, add more metadata.
        if "hub_prefix" in self.serverapp.tornado_settings:
            tornado_settings = self.serverapp.tornado_settings
            hub_prefix = tornado_settings["hub_prefix"]
            page_config["hubPrefix"] = hub_prefix
            page_config["hubHost"] = tornado_settings["hub_host"]
            page_config["hubUser"] = tornado_settings["user"]
            page_config["shareUrl"] = ujoin(hub_prefix, "user-redirect")
            # Assume the server_name property indicates running JupyterHub 1.0.
            if hasattr(self.serverapp, "server_name"):
                page_config["hubServerName"] = self.serverapp.server_name
            # avoid setting API token in page config
            # $JUPYTERHUB_API_TOKEN identifies the server, not the client
            # but at least make sure we don't use the token
            # if the serverapp set one
            page_config["token"] = ""

        # Update Jupyter Server's webapp settings with jupyterlab settings.
        self.serverapp.web_app.settings["page_config_data"] = page_config

        # Extend Server handlers with jupyterlab handlers.
        self.handlers.extend(handlers)
        super().initialize_handlers()

    def initialize(self, argv=None):
        """Subclass because the ExtensionApp.initialize() method does not take arguments"""
        super().initialize()
        if self.collaborative:
            try:
                import jupyter_collaboration  # noqa
            except ImportError:
                self.log.critical(
                    """\
Jupyter Lab cannot start, because `jupyter_collaboration` was configured but cannot be `import`ed.

To fix this, either:

1) install the extension `jupyter-collaboration`; for example: `python -m pip install jupyter-collaboration`

2) disable collaboration; for example, remove the `--collaborative` flag from the commandline.  To see more ways to adjust the collaborative behavior, see https://jupyterlab-realtime-collaboration.readthedocs.io/en/latest/configuration.html .
"""
                )
                sys.exit(1)


# -----------------------------------------------------------------------------
# Main entry point
# -----------------------------------------------------------------------------

main = launch_new_instance = LabApp.launch_instance

if __name__ == "__main__":
    main()
