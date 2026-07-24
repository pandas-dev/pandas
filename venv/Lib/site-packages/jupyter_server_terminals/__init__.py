from typing import Any, Dict, List

from ._version import __version__  # noqa:F401

try:
    from jupyter_server._version import version_info
except ModuleNotFoundError:
    msg = "Jupyter Server must be installed to use this extension."
    raise ModuleNotFoundError(msg) from None

if int(version_info[0]) < 2:  # type:ignore[call-overload]
    msg = "Jupyter Server Terminals requires Jupyter Server 2.0+"
    raise RuntimeError(msg)

from .app import TerminalsExtensionApp


def _jupyter_server_extension_points() -> List[Dict[str, Any]]:  # pragma: no cover
    return [
        {
            "module": "jupyter_server_terminals.app",
            "app": TerminalsExtensionApp,
        },
    ]
