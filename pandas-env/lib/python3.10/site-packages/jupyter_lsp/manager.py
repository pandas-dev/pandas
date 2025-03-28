""" A configurable frontend for stdio-based Language Servers
"""

import asyncio
import os
import sys
import traceback
from typing import Dict, Text, Tuple, cast

# See compatibility note on `group` keyword in
# https://docs.python.org/3/library/importlib.metadata.html#entry-points
if sys.version_info < (3, 10):  # pragma: no cover
    from importlib_metadata import entry_points
else:  # pragma: no cover
    from importlib.metadata import entry_points

from jupyter_core.paths import jupyter_config_path
from jupyter_server.services.config import ConfigManager

try:
    from jupyter_server.transutils import _i18n as _
except ImportError:  # pragma: no cover
    from jupyter_server.transutils import _

from traitlets import Bool
from traitlets import Dict as Dict_
from traitlets import Instance
from traitlets import List as List_
from traitlets import Unicode, default

from .constants import (
    APP_CONFIG_D_SECTIONS,
    EP_LISTENER_ALL_V1,
    EP_LISTENER_CLIENT_V1,
    EP_LISTENER_SERVER_V1,
    EP_SPEC_V1,
)
from .schema import LANGUAGE_SERVER_SPEC_MAP
from .session import LanguageServerSession
from .trait_types import LoadableCallable, Schema
from .types import (
    KeyedLanguageServerSpecs,
    LanguageServerManagerAPI,
    MessageScope,
    SpecBase,
    SpecMaker,
)


class LanguageServerManager(LanguageServerManagerAPI):
    """Manage language servers"""

    conf_d_language_servers = Schema(  # type:ignore[assignment]
        validator=LANGUAGE_SERVER_SPEC_MAP,
        help=_("extra language server specs, keyed by implementation, from conf.d"),
    )  # type: KeyedLanguageServerSpecs

    language_servers = Schema(  # type:ignore[assignment]
        validator=LANGUAGE_SERVER_SPEC_MAP,
        help=_("a dict of language server specs, keyed by implementation"),
    ).tag(
        config=True
    )  # type: KeyedLanguageServerSpecs

    autodetect: bool = Bool(  # type:ignore[assignment]
        True, help=_("try to find known language servers in sys.prefix (and elsewhere)")
    ).tag(config=True)

    sessions: Dict[Tuple[Text], LanguageServerSession] = (
        Dict_(  # type:ignore[assignment]
            trait=Instance(LanguageServerSession),
            default_value={},
            help="sessions keyed by language server name",
        )
    )

    virtual_documents_dir = Unicode(
        help="""Path to virtual documents relative to the content manager root
        directory.

        Its default value can be set with JP_LSP_VIRTUAL_DIR and fallback to
        '.virtual_documents'.
        """
    ).tag(config=True)

    _ready = Bool(
        help="""Whether the manager has been initialized""", default_value=False
    )

    all_listeners = List_(  # type:ignore[var-annotated]
        trait=LoadableCallable  # type:ignore[arg-type]
    ).tag(config=True)
    server_listeners = List_(  # type:ignore[var-annotated]
        trait=LoadableCallable  # type:ignore[arg-type]
    ).tag(config=True)
    client_listeners = List_(  # type:ignore[var-annotated]
        trait=LoadableCallable  # type:ignore[arg-type]
    ).tag(config=True)

    @default("language_servers")
    def _default_language_servers(self):
        return {}

    @default("virtual_documents_dir")
    def _default_virtual_documents_dir(self):
        return os.getenv("JP_LSP_VIRTUAL_DIR", ".virtual_documents")

    @default("conf_d_language_servers")
    def _default_conf_d_language_servers(self) -> KeyedLanguageServerSpecs:
        language_servers: KeyedLanguageServerSpecs = {}

        manager = ConfigManager(read_config_path=jupyter_config_path())

        for app in APP_CONFIG_D_SECTIONS:
            language_servers.update(
                **manager.get(f"jupyter{app}config")
                .get(self.__class__.__name__, {})
                .get("language_servers", {})
            )

        return language_servers

    def __init__(self, **kwargs: Dict):
        """Before starting, perform all necessary configuration"""
        self.all_language_servers: KeyedLanguageServerSpecs = {}
        self._language_servers_from_config: KeyedLanguageServerSpecs = {}
        super().__init__(**kwargs)

    def initialize(self, *args, **kwargs):
        self.init_language_servers()
        self.init_listeners()
        self.init_sessions()
        self._ready = True

    async def ready(self):
        while not self._ready:  # pragma: no cover
            await asyncio.sleep(0.1)
        return True

    def init_language_servers(self) -> None:
        """determine the final language server configuration."""
        # copy the language servers before anybody monkeys with them
        self._language_servers_from_config = dict(self.language_servers)
        self.language_servers = self._collect_language_servers(only_installed=True)
        self.all_language_servers = self._collect_language_servers(only_installed=False)

    def _collect_language_servers(
        self, only_installed: bool
    ) -> KeyedLanguageServerSpecs:
        language_servers: KeyedLanguageServerSpecs = {}

        language_servers_from_config = dict(self._language_servers_from_config)
        language_servers_from_config.update(self.conf_d_language_servers)

        if self.autodetect:
            language_servers.update(
                self._autodetect_language_servers(only_installed=only_installed)
            )

        # restore config
        language_servers.update(language_servers_from_config)

        # coalesce the servers, allowing a user to opt-out by specifying `[]`
        return {key: spec for key, spec in language_servers.items() if spec.get("argv")}

    def init_sessions(self):
        """create, but do not initialize all sessions"""
        sessions = {}
        for language_server, spec in self.language_servers.items():
            sessions[language_server] = LanguageServerSession(
                language_server=language_server, spec=spec, parent=self
            )
        self.sessions = sessions

    def init_listeners(self):
        """register traitlets-configured listeners"""

        scopes = {
            MessageScope.ALL: [self.all_listeners, EP_LISTENER_ALL_V1],
            MessageScope.CLIENT: [self.client_listeners, EP_LISTENER_CLIENT_V1],
            MessageScope.SERVER: [self.server_listeners, EP_LISTENER_SERVER_V1],
        }
        for scope, trt_ep in scopes.items():
            listeners, entry_point = trt_ep

            for ept in entry_points(group=entry_point):  # pragma: no cover
                try:
                    listeners.append(ept.load())
                except Exception as err:
                    self.log.warning("Failed to load entry point %s: %s", ept.name, err)

            for listener in listeners:
                self.__class__.register_message_listener(scope=scope.value)(listener)

    def subscribe(self, handler):
        """subscribe a handler to session, or sta"""
        session = self.sessions.get(handler.language_server)

        if session is None:
            self.log.error(
                "[{}] no session: handler subscription failed".format(
                    handler.language_server
                )
            )
            return

        session.handlers = set([handler]) | session.handlers

    async def on_client_message(self, message, handler):
        await self.wait_for_listeners(
            MessageScope.CLIENT, message, handler.language_server
        )
        session = self.sessions.get(handler.language_server)

        if session is None:
            self.log.error(
                "[{}] no session: client message dropped".format(
                    handler.language_server
                )
            )
            return

        session.write(message)

    async def on_server_message(self, message, session):
        language_servers = [
            ls_key for ls_key, sess in self.sessions.items() if sess == session
        ]

        for language_servers in language_servers:
            await self.wait_for_listeners(
                MessageScope.SERVER, message, language_servers
            )

        for handler in session.handlers:
            handler.write_message(message)

    def unsubscribe(self, handler):
        session = self.sessions.get(handler.language_server)

        if session is None:
            self.log.error(
                "[{}] no session: handler unsubscription failed".format(
                    handler.language_server
                )
            )
            return

        session.handlers = [h for h in session.handlers if h != handler]

    def _autodetect_language_servers(self, only_installed: bool):
        _entry_points = None

        try:
            _entry_points = entry_points(group=EP_SPEC_V1)
        except Exception:  # pragma: no cover
            self.log.exception("Failed to load entry_points")

        skipped_servers = []

        for ep in _entry_points or []:
            try:
                spec_finder: SpecMaker = ep.load()
            except Exception as err:  # pragma: no cover
                self.log.warning(
                    _("Failed to load language server spec finder `{}`: \n{}").format(
                        ep.name, err
                    )
                )
                continue

            try:
                if only_installed:
                    if hasattr(spec_finder, "is_installed"):
                        spec_finder_from_base = cast(SpecBase, spec_finder)
                        if not spec_finder_from_base.is_installed(self):
                            skipped_servers.append(ep.name)
                            continue
                specs = spec_finder(self) or {}
            except Exception as err:  # pragma: no cover
                self.log.warning(
                    _(
                        "Failed to fetch commands from language server spec finder"
                        " `{}`:\n{}"
                    ).format(ep.name, err)
                )
                traceback.print_exc()

                continue

            errors = list(LANGUAGE_SERVER_SPEC_MAP.iter_errors(specs))

            if errors:  # pragma: no cover
                self.log.warning(
                    _(
                        "Failed to validate commands from language server spec finder"
                        " `{}`:\n{}"
                    ).format(ep.name, errors)
                )
                continue

            for key, spec in specs.items():
                yield key, spec

        if skipped_servers:
            self.log.info(
                _("Skipped non-installed server(s): {}").format(
                    ", ".join(skipped_servers)
                )
            )


# the listener decorator
lsp_message_listener = LanguageServerManager.register_message_listener  # noqa
