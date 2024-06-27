try:
    from ._version import __version__  # noqa
except ImportError:
    # Fallback when using the package in dev mode without installing
    # in editable mode with pip. Here this is particularly important
    # to be able to run the generate_css.py script.
    __version__ = "dev"
from .style import JupyterStyle  # noqa


def _jupyter_labextension_paths():
    return [{
        "src": "labextension",
        "dest": "jupyterlab_pygments"
    }]
