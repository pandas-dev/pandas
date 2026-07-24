# flake8: noqa: F401
from ._version import __version__
from .manager import LanguageServerManager, lsp_message_listener
from .serverextension import load_jupyter_server_extension
from .specs.utils import NodeModuleSpec, ShellSpec
from .types import (
    KeyedLanguageServerSpecs,
    LanguageServerManagerAPI,
    LanguageServerSpec,
)


def _jupyter_server_extension_paths():
    return [{"module": "jupyter_lsp"}]


_jupyter_server_extension_points = _jupyter_server_extension_paths
