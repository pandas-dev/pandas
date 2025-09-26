"""A JupyterLite kernel powered by Pyodide."""

from ._version import __version__

__all__ = ["__version__", "_jupyter_labextension_paths"]


def _jupyter_labextension_paths():
    from .constants import PYODIDE_KERNEL_NPM_NAME

    return [{"src": "labextension", "dest": PYODIDE_KERNEL_NPM_NAME}]
