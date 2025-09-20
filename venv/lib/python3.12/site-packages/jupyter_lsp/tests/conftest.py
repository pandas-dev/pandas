import json
import os
import pathlib
import shutil
from pathlib import Path
from typing import Text

from jupyter_server.serverapp import ServerApp
from pytest import fixture
from tornado.httpserver import HTTPRequest
from tornado.httputil import HTTPServerRequest
from tornado.queues import Queue
from tornado.web import Application

# local imports
from jupyter_lsp import LanguageServerManager
from jupyter_lsp.constants import APP_CONFIG_D_SECTIONS
from jupyter_lsp.handlers import LanguageServersHandler, LanguageServerWebSocketHandler

# these should always be available in a test environment
KNOWN_SERVERS = [
    "bash-language-server",
    "dockerfile-language-server-nodejs",
    "typescript-language-server",
    "pylsp",
    "unified-language-server",
    "sql-language-server",
    "vscode-css-languageserver-bin",
    "vscode-html-languageserver-bin",
    "vscode-json-languageserver-bin",
    "yaml-language-server",
]

CMD_BASED_SERVERS = {
    "Rscript": ["r-languageserver"],
    "texlab": ["texlab"],
    "jedi-language-server": ["jedi-language-server"],
    "julia": ["julia-language-server"],
}

KNOWN_SERVERS += sum(
    [langs for cmd, langs in CMD_BASED_SERVERS.items() if shutil.which(cmd)], []
)

KNOWN_UNKNOWN_SERVERS = ["foo-language-server"]


def extra_node_roots():
    root = Path(os.environ.get("JLSP_TEST_ROOT") or Path.cwd())
    return dict(extra_node_roots=[str(root)] if root else [])


@fixture
def manager() -> LanguageServerManager:
    return LanguageServerManager(**extra_node_roots())


@fixture
def echo_spec():
    return {"argv": ["echo", "no server here"], "languages": ["klingon"], "version": 2}


@fixture
def echo_conf_json(echo_spec) -> str:
    return json.dumps(
        {"LanguageServerManager": {"language_servers": {"_echo_": echo_spec}}},
        indent=2,
        sort_keys=True,
    )


@fixture(params=sorted(APP_CONFIG_D_SECTIONS))
def app_config_d(request, tmp_path, monkeypatch) -> pathlib.Path:
    conf_d = tmp_path / f"jupyter{request.param}config.d"
    conf_d.mkdir()
    monkeypatch.setenv("JUPYTER_CONFIG_PATH", f"{tmp_path}")
    return conf_d


@fixture(params=sorted(KNOWN_SERVERS))
def known_server(request):
    return request.param


@fixture(params=sorted(KNOWN_UNKNOWN_SERVERS))
def known_unknown_server(request):
    return request.param


@fixture
def handlers(manager):
    ws_handler = MockWebsocketHandler()
    ws_handler.initialize(manager)
    handler = MockHandler()
    handler.initialize(manager)
    return handler, ws_handler


@fixture
def jsonrpc_init_msg():
    return json.dumps(
        {
            "id": 0,
            "jsonrpc": "2.0",
            "method": "initialize",
            "params": {
                "capabilities": {
                    # see: https://github.com/julia-vscode/LanguageServer.jl/issues/1008
                    # LanguageServer.jl assumes that it is not missing
                    "workspace": {"didChangeConfiguration": {}},
                    # LanguageServer.jl assumes that it is not missing
                    "textDocument": {},
                },
                "initializationOptions": None,
                "processId": None,
                "rootUri": pathlib.Path(__file__).parent.as_uri(),
                "workspaceFolders": None,
            },
        }
    )


@fixture
def app():
    return MockServerApp()


# mocks
class MockWebsocketHandler(LanguageServerWebSocketHandler):
    _messages_wrote = None  # type: Queue
    _ping_sent = None  # type: bool

    def __init__(self):
        self.request = HTTPServerRequest()
        self.application = Application()

    def initialize(self, manager):
        super().initialize(manager)
        self._messages_wrote = Queue()
        self._ping_sent = False

    def write_message(self, message: Text) -> None:  # type: ignore
        self.log.warning("write_message %s", message)
        self._messages_wrote.put_nowait(message)

    def send_ping(self):
        self._ping_sent = True


class MockHandler(LanguageServersHandler):
    _payload = None
    _jupyter_current_user = "foo"  # type:ignore[assignment]

    def __init__(self):
        self.request = HTTPRequest("GET")
        self.application = Application()

    def finish(self, payload):
        self._payload = payload


class MockServerApp(ServerApp):
    pass
