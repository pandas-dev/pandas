from __future__ import annotations

from typing import Any

from ._version import __version__, version_info  # noqa: F401


def _jupyter_server_extension_paths() -> list[dict[str, str]]:
    return [{"module": "notebook"}]


def _jupyter_server_extension_points() -> list[dict[str, Any]]:
    from .app import JupyterNotebookApp

    return [{"module": "notebook", "app": JupyterNotebookApp}]


def _jupyter_labextension_paths() -> list[dict[str, str]]:
    return [{"src": "labextension", "dest": "@jupyter-notebook/lab-extension"}]
