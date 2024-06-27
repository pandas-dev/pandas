"""An extension application."""

from __future__ import annotations

import logging
import re
import sys
import typing as t

from jinja2 import Environment, FileSystemLoader
from jupyter_core.application import JupyterApp, NoStart
from tornado.log import LogFormatter
from tornado.web import RedirectHandler
from traitlets import Any, Bool, Dict, HasTraits, List, Unicode, default
from traitlets.config import Config

from jupyter_server.serverapp import ServerApp
from jupyter_server.transutils import _i18n
from jupyter_server.utils import is_namespace_package, url_path_join

from .handler import ExtensionHandlerMixin

# -----------------------------------------------------------------------------
# Util functions and classes.
# -----------------------------------------------------------------------------


def _preparse_for_subcommand(application_klass, argv):
    """Preparse command line to look for subcommands."""
    # Read in arguments from command line.
    if len(argv) == 0:
        return None

    # Find any subcommands.
    if application_klass.subcommands and len(argv) > 0:
        # we have subcommands, and one may have been specified
        subc, subargv = argv[0], argv[1:]
        if re.match(r"^\w(\-?\w)*$", subc) and subc in application_klass.subcommands:
            # it's a subcommand, and *not* a flag or class parameter
            app = application_klass()
            app.initialize_subcommand(subc, subargv)
            return app.subapp


def _preparse_for_stopping_flags(application_klass, argv):
    """Looks for 'help', 'version', and 'generate-config; commands
    in command line. If found, raises the help and version of
    current Application.

    This is useful for traitlets applications that have to parse
    the command line multiple times, but want to control when
    when 'help' and 'version' is raised.
    """
    # Arguments after a '--' argument are for the script IPython may be
    # about to run, not IPython iteslf. For arguments parsed here (help and
    # version), we want to only search the arguments up to the first
    # occurrence of '--', which we're calling interpreted_argv.
    try:
        interpreted_argv = argv[: argv.index("--")]
    except ValueError:
        interpreted_argv = argv

    # Catch any help calls.
    if any(x in interpreted_argv for x in ("-h", "--help-all", "--help")):
        app = application_klass()
        app.print_help("--help-all" in interpreted_argv)
        app.exit(0)

    # Catch version commands
    if "--version" in interpreted_argv or "-V" in interpreted_argv:
        app = application_klass()
        app.print_version()
        app.exit(0)

    # Catch generate-config commands.
    if "--generate-config" in interpreted_argv:
        app = application_klass()
        app.write_default_config()
        app.exit(0)


class ExtensionAppJinjaMixin(HasTraits):
    """Use Jinja templates for HTML templates on top of an ExtensionApp."""

    jinja2_options = Dict(
        help=_i18n(
            """Options to pass to the jinja2 environment for this
        """
        )
    ).tag(config=True)

    @t.no_type_check
    def _prepare_templates(self):
        """Get templates defined in a subclass."""
        self.initialize_templates()
        # Add templates to web app settings if extension has templates.
        if len(self.template_paths) > 0:
            self.settings.update({f"{self.name}_template_paths": self.template_paths})

        # Create a jinja environment for logging html templates.
        self.jinja2_env = Environment(
            loader=FileSystemLoader(self.template_paths),
            extensions=["jinja2.ext.i18n"],
            autoescape=True,
            **self.jinja2_options,
        )

        # Add the jinja2 environment for this extension to the tornado settings.
        self.settings.update({f"{self.name}_jinja2_env": self.jinja2_env})


# -----------------------------------------------------------------------------
# ExtensionApp
# -----------------------------------------------------------------------------


class JupyterServerExtensionException(Exception):
    """Exception class for raising for Server extensions errors."""


# -----------------------------------------------------------------------------
# ExtensionApp
# -----------------------------------------------------------------------------


class ExtensionApp(JupyterApp):
    """Base class for configurable Jupyter Server Extension Applications.

    ExtensionApp subclasses can be initialized two ways:

    - Extension is listed as a jpserver_extension, and ServerApp calls
      its load_jupyter_server_extension classmethod. This is the
      classic way of loading a server extension.

    - Extension is launched directly by calling its `launch_instance`
      class method. This method can be set as a entry_point in
      the extensions setup.py.
    """

    # Subclasses should override this trait. Tells the server if
    # this extension allows other other extensions to be loaded
    # side-by-side when launched directly.
    load_other_extensions = True

    # A useful class property that subclasses can override to
    # configure the underlying Jupyter Server when this extension
    # is launched directly (using its `launch_instance` method).
    serverapp_config: dict[str, t.Any] = {}

    # Some subclasses will likely override this trait to flip
    # the default value to False if they don't offer a browser
    # based frontend.
    open_browser = Bool(
        help="""Whether to open in a browser after starting.
        The specific browser used is platform dependent and
        determined by the python standard library `webbrowser`
        module, unless it is overridden using the --browser
        (ServerApp.browser) configuration option.
        """
    ).tag(config=True)

    @default("open_browser")
    def _default_open_browser(self):
        assert self.serverapp is not None
        return self.serverapp.config["ServerApp"].get("open_browser", True)

    @property
    def config_file_paths(self):
        """Look on the same path as our parent for config files"""
        # rely on parent serverapp, which should control all config loading
        assert self.serverapp is not None
        return self.serverapp.config_file_paths

    # The extension name used to name the jupyter config
    # file, jupyter_{name}_config.
    # This should also match the jupyter subcommand used to launch
    # this extension from the CLI, e.g. `jupyter {name}`.
    name: str | Unicode[str, str] = "ExtensionApp"  # type:ignore[assignment]

    @classmethod
    def get_extension_package(cls):
        """Get an extension package."""
        parts = cls.__module__.split(".")
        if is_namespace_package(parts[0]):
            # in this case the package name is `<namespace>.<package>`.
            return ".".join(parts[0:2])
        return parts[0]

    @classmethod
    def get_extension_point(cls):
        """Get an extension point."""
        return cls.__module__

    # Extension URL sets the default landing page for this extension.
    extension_url = "/"

    default_url = Unicode().tag(config=True)

    @default("default_url")
    def _default_url(self):
        return self.extension_url

    file_url_prefix = Unicode("notebooks")

    # Is this linked to a serverapp yet?
    _linked = Bool(False)

    # Extension can configure the ServerApp from the command-line
    classes = [
        ServerApp,
    ]

    # A ServerApp is not defined yet, but will be initialized below.
    serverapp: ServerApp | None = Any()  # type:ignore[assignment]

    @default("serverapp")
    def _default_serverapp(self):
        # load the current global instance, if any
        if ServerApp.initialized():
            try:
                return ServerApp.instance()
            except Exception:
                # error retrieving instance, e.g. MultipleInstanceError
                pass

        # serverapp accessed before it was defined,
        # declare an empty one
        return ServerApp()

    _log_formatter_cls = LogFormatter  # type:ignore[assignment]

    @default("log_level")
    def _default_log_level(self):
        return logging.INFO

    @default("log_format")
    def _default_log_format(self):
        """override default log format to include date & time"""
        return (
            "%(color)s[%(levelname)1.1s %(asctime)s.%(msecs).03d %(name)s]%(end_color)s %(message)s"
        )

    static_url_prefix = Unicode(
        help="""Url where the static assets for the extension are served."""
    ).tag(config=True)

    @default("static_url_prefix")
    def _default_static_url_prefix(self):
        static_url = f"static/{self.name}/"
        assert self.serverapp is not None
        return url_path_join(self.serverapp.base_url, static_url)

    static_paths = List(
        Unicode(),
        help="""paths to search for serving static files.

        This allows adding javascript/css to be available from the notebook server machine,
        or overriding individual files in the IPython
        """,
    ).tag(config=True)

    template_paths = List(
        Unicode(),
        help=_i18n(
            """Paths to search for serving jinja templates.

        Can be used to override templates from notebook.templates."""
        ),
    ).tag(config=True)

    settings = Dict(help=_i18n("""Settings that will passed to the server.""")).tag(config=True)

    handlers: List[tuple[t.Any, ...]] = List(
        help=_i18n("""Handlers appended to the server.""")
    ).tag(config=True)

    def _config_file_name_default(self):
        """The default config file name."""
        if not self.name:
            return ""
        return "jupyter_{}_config".format(self.name.replace("-", "_"))

    def initialize_settings(self):
        """Override this method to add handling of settings."""

    def initialize_handlers(self):
        """Override this method to append handlers to a Jupyter Server."""

    def initialize_templates(self):
        """Override this method to add handling of template files."""

    def _prepare_config(self):
        """Builds a Config object from the extension's traits and passes
        the object to the webapp's settings as `<name>_config`.
        """
        traits = self.class_own_traits().keys()
        self.extension_config = Config({t: getattr(self, t) for t in traits})
        self.settings[f"{self.name}_config"] = self.extension_config

    def _prepare_settings(self):
        """Prepare the settings."""
        # Make webapp settings accessible to initialize_settings method
        assert self.serverapp is not None
        webapp = self.serverapp.web_app
        self.settings.update(**webapp.settings)

        # Add static and template paths to settings.
        self.settings.update(
            {
                f"{self.name}_static_paths": self.static_paths,
                f"{self.name}": self,
            }
        )

        # Get setting defined by subclass using initialize_settings method.
        self.initialize_settings()

        # Update server settings with extension settings.
        webapp.settings.update(**self.settings)

    def _prepare_handlers(self):
        """Prepare the handlers."""
        assert self.serverapp is not None
        webapp = self.serverapp.web_app

        # Get handlers defined by extension subclass.
        self.initialize_handlers()

        # prepend base_url onto the patterns that we match
        new_handlers = []
        for handler_items in self.handlers:
            # Build url pattern including base_url
            pattern = url_path_join(webapp.settings["base_url"], handler_items[0])
            handler = handler_items[1]

            # Get handler kwargs, if given
            kwargs: dict[str, t.Any] = {}
            if issubclass(handler, ExtensionHandlerMixin):
                kwargs["name"] = self.name

            try:
                kwargs.update(handler_items[2])
            except IndexError:
                pass

            new_handler = (pattern, handler, kwargs)
            new_handlers.append(new_handler)

        # Add static endpoint for this extension, if static paths are given.
        if len(self.static_paths) > 0:
            # Append the extension's static directory to server handlers.
            static_url = url_path_join(self.static_url_prefix, "(.*)")

            # Construct handler.
            handler = (
                static_url,
                webapp.settings["static_handler_class"],
                {"path": self.static_paths},
            )
            new_handlers.append(handler)

        webapp.add_handlers(".*$", new_handlers)

    def _prepare_templates(self):
        """Add templates to web app settings if extension has templates."""
        if len(self.template_paths) > 0:
            self.settings.update({f"{self.name}_template_paths": self.template_paths})
        self.initialize_templates()

    def _jupyter_server_config(self):
        """The jupyter server config."""
        base_config = {
            "ServerApp": {
                "default_url": self.default_url,
                "open_browser": self.open_browser,
                "file_url_prefix": self.file_url_prefix,
            }
        }
        base_config["ServerApp"].update(self.serverapp_config)
        return base_config

    def _link_jupyter_server_extension(self, serverapp: ServerApp) -> None:
        """Link the ExtensionApp to an initialized ServerApp.

        The ServerApp is stored as an attribute and config
        is exchanged between ServerApp and `self` in case
        the command line contains traits for the ExtensionApp
        or the ExtensionApp's config files have server
        settings.

        Note, the ServerApp has not initialized the Tornado
        Web Application yet, so do not try to affect the
        `web_app` attribute.
        """
        self.serverapp = serverapp
        # Load config from an ExtensionApp's config files.
        self.load_config_file()
        # ServerApp's config might have picked up
        # config for the ExtensionApp. We call
        # update_config to update ExtensionApp's
        # traits with these values found in ServerApp's
        # config.
        # ServerApp config ---> ExtensionApp traits
        self.update_config(self.serverapp.config)
        # Use ExtensionApp's CLI parser to find any extra
        # args that passed through ServerApp and
        # now belong to ExtensionApp.
        self.parse_command_line(self.serverapp.extra_args)
        # If any config should be passed upstream to the
        # ServerApp, do it here.
        # i.e. ServerApp traits <--- ExtensionApp config
        self.serverapp.update_config(self.config)
        # Acknowledge that this extension has been linked.
        self._linked = True

    def initialize(self):
        """Initialize the extension app. The
        corresponding server app and webapp should already
        be initialized by this step.

        - Appends Handlers to the ServerApp,
        - Passes config and settings from ExtensionApp
          to the Tornado web application
        - Points Tornado Webapp to templates and static assets.
        """
        if not self.serverapp:
            msg = (
                "This extension has no attribute `serverapp`. "
                "Try calling `.link_to_serverapp()` before calling "
                "`.initialize()`."
            )
            raise JupyterServerExtensionException(msg)

        self._prepare_config()
        self._prepare_templates()
        self._prepare_settings()
        self._prepare_handlers()

    def start(self):
        """Start the underlying Jupyter server.

        Server should be started after extension is initialized.
        """
        super().start()
        # Start the server.
        assert self.serverapp is not None
        self.serverapp.start()

    def current_activity(self):
        """Return a list of activity happening in this extension."""
        return

    async def stop_extension(self):
        """Cleanup any resources managed by this extension."""

    def stop(self):
        """Stop the underlying Jupyter server."""
        assert self.serverapp is not None
        self.serverapp.stop()
        self.serverapp.clear_instance()

    @classmethod
    def _load_jupyter_server_extension(cls, serverapp):
        """Initialize and configure this extension, then add the extension's
        settings and handlers to the server's web application.
        """
        extension_manager = serverapp.extension_manager
        try:
            # Get loaded extension from serverapp.
            point = extension_manager.extension_points[cls.name]
            extension = point.app
        except KeyError:
            extension = cls()
            extension._link_jupyter_server_extension(serverapp)
        extension.initialize()
        return extension

    @classmethod
    def load_classic_server_extension(cls, serverapp):
        """Enables extension to be loaded as classic Notebook (jupyter/notebook) extension."""
        extension = cls()
        extension.serverapp = serverapp
        extension.load_config_file()
        extension.update_config(serverapp.config)
        extension.parse_command_line(serverapp.extra_args)
        # Add redirects to get favicons from old locations in the classic notebook server
        extension.handlers.extend(
            [
                (
                    r"/static/favicons/favicon.ico",
                    RedirectHandler,
                    {"url": url_path_join(serverapp.base_url, "static/base/images/favicon.ico")},
                ),
                (
                    r"/static/favicons/favicon-busy-1.ico",
                    RedirectHandler,
                    {
                        "url": url_path_join(
                            serverapp.base_url, "static/base/images/favicon-busy-1.ico"
                        )
                    },
                ),
                (
                    r"/static/favicons/favicon-busy-2.ico",
                    RedirectHandler,
                    {
                        "url": url_path_join(
                            serverapp.base_url, "static/base/images/favicon-busy-2.ico"
                        )
                    },
                ),
                (
                    r"/static/favicons/favicon-busy-3.ico",
                    RedirectHandler,
                    {
                        "url": url_path_join(
                            serverapp.base_url, "static/base/images/favicon-busy-3.ico"
                        )
                    },
                ),
                (
                    r"/static/favicons/favicon-file.ico",
                    RedirectHandler,
                    {
                        "url": url_path_join(
                            serverapp.base_url, "static/base/images/favicon-file.ico"
                        )
                    },
                ),
                (
                    r"/static/favicons/favicon-notebook.ico",
                    RedirectHandler,
                    {
                        "url": url_path_join(
                            serverapp.base_url,
                            "static/base/images/favicon-notebook.ico",
                        )
                    },
                ),
                (
                    r"/static/favicons/favicon-terminal.ico",
                    RedirectHandler,
                    {
                        "url": url_path_join(
                            serverapp.base_url,
                            "static/base/images/favicon-terminal.ico",
                        )
                    },
                ),
                (
                    r"/static/logo/logo.png",
                    RedirectHandler,
                    {"url": url_path_join(serverapp.base_url, "static/base/images/logo.png")},
                ),
            ]
        )
        extension.initialize()

    serverapp_class = ServerApp

    @classmethod
    def make_serverapp(cls, **kwargs: t.Any) -> ServerApp:
        """Instantiate the ServerApp

        Override to customize the ServerApp before it loads any configuration
        """
        return cls.serverapp_class.instance(**kwargs)

    @classmethod
    def initialize_server(cls, argv=None, load_other_extensions=True, **kwargs):
        """Creates an instance of ServerApp and explicitly sets
        this extension to enabled=True (i.e. superseding disabling
        found in other config from files).

        The `launch_instance` method uses this method to initialize
        and start a server.
        """
        jpserver_extensions = {cls.get_extension_package(): True}
        find_extensions = cls.load_other_extensions
        if "jpserver_extensions" in cls.serverapp_config:
            jpserver_extensions.update(cls.serverapp_config["jpserver_extensions"])
            cls.serverapp_config["jpserver_extensions"] = jpserver_extensions
            find_extensions = False
        serverapp = cls.make_serverapp(jpserver_extensions=jpserver_extensions, **kwargs)
        serverapp.aliases.update(cls.aliases)  # type:ignore[has-type]
        serverapp.initialize(
            argv=argv or [],
            starter_extension=cls.name,
            find_extensions=find_extensions,
        )
        return serverapp

    @classmethod
    def launch_instance(cls, argv=None, **kwargs):
        """Launch the extension like an application. Initializes+configs a stock server
        and appends the extension to the server. Then starts the server and routes to
        extension's landing page.
        """
        # Handle arguments.
        if argv is None:  # noqa: SIM108
            args = sys.argv[1:]  # slice out extension config.
        else:
            args = argv

        # Handle all "stops" that could happen before
        # continuing to launch a server+extension.
        subapp = _preparse_for_subcommand(cls, args)
        if subapp:
            subapp.start()
            return

        # Check for help, version, and generate-config arguments
        # before initializing server to make sure these
        # arguments trigger actions from the extension not the server.
        _preparse_for_stopping_flags(cls, args)
        serverapp = cls.initialize_server(argv=args)

        # Log if extension is blocking other extensions from loading.
        if not cls.load_other_extensions:
            serverapp.log.info(f"{cls.name} is running without loading other extensions.")
        # Start the server.
        try:
            serverapp.start()
        except NoStart:
            pass
