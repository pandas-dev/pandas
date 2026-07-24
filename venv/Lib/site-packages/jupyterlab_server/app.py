"""JupyterLab Server Application"""

# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.

from glob import glob
from os.path import relpath
from typing import Any

from jupyter_server.extension.application import ExtensionApp, ExtensionAppJinjaMixin
from jupyter_server.utils import url_path_join as ujoin
from traitlets import Dict, Integer, Unicode, observe

from ._version import __version__
from .handlers import LabConfig, add_handlers


class LabServerApp(ExtensionAppJinjaMixin, LabConfig, ExtensionApp):
    """A Lab Server Application that runs out-of-the-box"""

    name = "jupyterlab_server"
    extension_url = "/lab"
    app_name = "JupyterLab Server Application"  # type:ignore[assignment]
    file_url_prefix = "/lab/tree"  # type:ignore[assignment]

    @property
    def app_namespace(self) -> str:  # type:ignore[override]
        return self.name

    default_url = Unicode("/lab", help="The default URL to redirect to from `/`")

    # Should your extension expose other server extensions when launched directly?
    load_other_extensions = True

    app_version = Unicode("", help="The version of the application.").tag(default=__version__)

    blacklist_uris = Unicode(
        "", config=True, help="Deprecated, use `LabServerApp.blocked_extensions_uris`"
    )

    blocked_extensions_uris = Unicode(
        "",
        config=True,
        help="""
        A list of comma-separated URIs to get the blocked extensions list

        .. versionchanged:: 2.0.0
            `LabServerApp.blacklist_uris` renamed to `blocked_extensions_uris`
        """,
    )

    whitelist_uris = Unicode(
        "", config=True, help="Deprecated, use `LabServerApp.allowed_extensions_uris`"
    )

    allowed_extensions_uris = Unicode(
        "",
        config=True,
        help="""
        "A list of comma-separated URIs to get the allowed extensions list

        .. versionchanged:: 2.0.0
            `LabServerApp.whitetlist_uris` renamed to `allowed_extensions_uris`
        """,
    )

    listings_refresh_seconds = Integer(
        60 * 60, config=True, help="The interval delay in seconds to refresh the lists"
    )

    listings_request_options = Dict(
        {},
        config=True,
        help="The optional kwargs to use for the listings HTTP requests \
            as described on https://2.python-requests.org/en/v2.7.0/api/#requests.request",
    )

    _deprecated_aliases = {
        "blacklist_uris": ("blocked_extensions_uris", "1.2"),
        "whitelist_uris": ("allowed_extensions_uris", "1.2"),
    }

    # Method copied from
    # https://github.com/jupyterhub/jupyterhub/blob/d1a85e53dccfc7b1dd81b0c1985d158cc6b61820/jupyterhub/auth.py#L143-L161
    @observe(*list(_deprecated_aliases))
    def _deprecated_trait(self, change: Any) -> None:
        """observer for deprecated traits"""
        old_attr = change.name
        new_attr, version = self._deprecated_aliases.get(old_attr)  # type:ignore[misc]
        new_value = getattr(self, new_attr)
        if new_value != change.new:
            # only warn if different
            # protects backward-compatible config from warnings
            # if they set the same value under both names
            self.log.warning(
                "%s.%s is deprecated in JupyterLab %s, use %s.%s instead",
                self.__class__.__name__,
                old_attr,
                version,
                self.__class__.__name__,
                new_attr,
            )

            setattr(self, new_attr, change.new)

    def initialize_settings(self) -> None:
        """Initialize the settings:

        set the static files as immutable, since they should have all hashed name.
        """
        immutable_cache = set(self.settings.get("static_immutable_cache", []))

        # Set lab static files as immutables
        immutable_cache.add(self.static_url_prefix)

        # Set extensions static files as immutables
        for extension_path in self.labextensions_path + self.extra_labextensions_path:
            extensions_url = [
                ujoin(self.labextensions_url, relpath(path, extension_path))
                for path in glob(f"{extension_path}/**/static", recursive=True)
            ]

            immutable_cache.update(extensions_url)

        self.settings.update({"static_immutable_cache": list(immutable_cache)})
        if self.serverapp:
            untracked_message_types = getattr(
                self.serverapp.kernel_manager, "untracked_message_types", None
            )
            if untracked_message_types:
                web_app = self.serverapp.web_app
                page_config_data = web_app.settings.setdefault("page_config_data", {})
                page_config_data["untracked_message_types"] = list(untracked_message_types)

    def initialize_templates(self) -> None:
        """Initialize templates."""
        self.static_paths = [self.static_dir]
        self.template_paths = [self.templates_dir]

    def initialize_handlers(self) -> None:
        """Initialize handlers."""
        add_handlers(self.handlers, self)


main = launch_new_instance = LabServerApp.launch_instance
