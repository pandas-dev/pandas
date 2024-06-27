"""\
The ``ansicolor`` object provides ``bg`` and ``fg`` to access colors,
and attributes like bold and
underlined text. It also provides ``reset`` to recover the normal font.
"""

import sys

from .factories import StyleFactory
from .styles import ANSIStyle, ColorNotFound, HTMLStyle, Style

__all__ = (
    "ANSIStyle",
    "ColorNotFound",
    "HTMLStyle",
    "Style",
    "StyleFactory",
    "ansicolors",
    "htmlcolors",
    "load_ipython_extension",
    "main",
)

ansicolors = StyleFactory(ANSIStyle)
htmlcolors = StyleFactory(HTMLStyle)


def load_ipython_extension(ipython):  # pragma: no cover
    try:
        from ._ipython_ext import OutputMagics  # pylint:disable=import-outside-toplevel
    except ImportError:
        print("IPython required for the IPython extension to be loaded.")  # noqa: T201
        raise

    ipython.push({"colors": htmlcolors})
    ipython.register_magics(OutputMagics)


def main():  # pragma: no cover
    """Color changing script entry. Call using
    python3 -m plumbum.colors, will reset if no arguments given."""
    color = " ".join(sys.argv[1:]) if len(sys.argv) > 1 else ""
    ansicolors.use_color = True
    ansicolors.get_colors_from_string(color).now()
