"""An extension handler."""

from __future__ import annotations

from logging import Logger
from typing import TYPE_CHECKING, Any, cast

from jinja2 import Template
from jinja2.exceptions import TemplateNotFound

from jupyter_server.base.handlers import FileFindHandler

if TYPE_CHECKING:
    from traitlets.config import Config

    from jupyter_server.extension.application import ExtensionApp
    from jupyter_server.serverapp import ServerApp


class ExtensionHandlerJinjaMixin:
    """Mixin class for ExtensionApp handlers that use jinja templating for
    template rendering.
    """

    def get_template(self, name: str) -> Template:
        """Return the jinja template object for a given name"""
        try:
            env = f"{self.name}_jinja2_env"  # type:ignore[attr-defined]
            template = cast(Template, self.settings[env].get_template(name))  # type:ignore[attr-defined]
            return template
        except TemplateNotFound:
            return cast(Template, super().get_template(name))  # type:ignore[misc]


class ExtensionHandlerMixin:
    """Base class for Jupyter server extension handlers.

    Subclasses can serve static files behind a namespaced
    endpoint: "<base_url>/static/<name>/"

    This allows multiple extensions to serve static files under
    their own namespace and avoid intercepting requests for
    other extensions.
    """

    settings: dict[str, Any]

    def initialize(self, name: str, *args: Any, **kwargs: Any) -> None:
        self.name = name
        try:
            super().initialize(*args, **kwargs)  # type:ignore[misc]
        except TypeError:
            pass

    @property
    def extensionapp(self) -> ExtensionApp:
        return cast("ExtensionApp", self.settings[self.name])

    @property
    def serverapp(self) -> ServerApp:
        key = "serverapp"
        return cast("ServerApp", self.settings[key])

    @property
    def log(self) -> Logger:
        if not hasattr(self, "name"):
            return cast(Logger, super().log)  # type:ignore[misc]
        # Attempt to pull the ExtensionApp's log, otherwise fall back to ServerApp.
        try:
            return cast(Logger, self.extensionapp.log)
        except AttributeError:
            return cast(Logger, self.serverapp.log)

    @property
    def config(self) -> Config:
        return cast("Config", self.settings[f"{self.name}_config"])

    @property
    def server_config(self) -> Config:
        return cast("Config", self.settings["config"])

    @property
    def base_url(self) -> str:
        return cast(str, self.settings.get("base_url", "/"))

    def render_template(self, name: str, **ns) -> str:
        """Override render template to handle static_paths

        If render_template is called with a template from the base environment
        (e.g. default error pages)
        make sure our extension-specific static_url is _not_ used.
        """
        template = cast(Template, self.get_template(name))  # type:ignore[attr-defined]
        ns.update(self.template_namespace)  # type:ignore[attr-defined]
        if template.environment is self.settings["jinja2_env"]:
            # default template environment, use default static_url
            ns["static_url"] = super().static_url  # type:ignore[misc]
        return cast(str, template.render(**ns))

    @property
    def static_url_prefix(self) -> str:
        return self.extensionapp.static_url_prefix

    @property
    def static_path(self) -> str:
        return cast(str, self.settings[f"{self.name}_static_paths"])

    def static_url(self, path: str, include_host: bool | None = None, **kwargs: Any) -> str:
        """Returns a static URL for the given relative static file path.
        This method requires you set the ``{name}_static_path``
        setting in your extension (which specifies the root directory
        of your static files).
        This method returns a versioned url (by default appending
        ``?v=<signature>``), which allows the static files to be
        cached indefinitely.  This can be disabled by passing
        ``include_version=False`` (in the default implementation;
        other static file implementations are not required to support
        this, but they may support other options).
        By default this method returns URLs relative to the current
        host, but if ``include_host`` is true the URL returned will be
        absolute.  If this handler has an ``include_host`` attribute,
        that value will be used as the default for all `static_url`
        calls that do not pass ``include_host`` as a keyword argument.
        """
        key = f"{self.name}_static_paths"
        try:
            self.require_setting(key, "static_url")  # type:ignore[attr-defined]
        except Exception as e:
            if key in self.settings:
                msg = (
                    "This extension doesn't have any static paths listed. Check that the "
                    "extension's `static_paths` trait is set."
                )
                raise Exception(msg) from None
            else:
                raise e

        get_url = self.settings.get("static_handler_class", FileFindHandler).make_static_url

        if include_host is None:
            include_host = getattr(self, "include_host", False)

        base = ""
        if include_host:
            base = self.request.protocol + "://" + self.request.host  # type:ignore[attr-defined]

        # Hijack settings dict to send extension templates to extension
        # static directory.
        settings = {
            "static_path": self.static_path,
            "static_url_prefix": self.static_url_prefix,
        }

        return base + cast(str, get_url(settings, path, **kwargs))
