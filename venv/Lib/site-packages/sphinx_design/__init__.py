"""A sphinx extension for designing beautiful, view size responsive web components."""

from typing import TYPE_CHECKING

__version__ = "0.7.0"

if TYPE_CHECKING:
    from sphinx.application import Sphinx


def setup(app: "Sphinx") -> dict:
    from .extension import setup_extension  # noqa: PLC0415

    setup_extension(app)
    return {
        "version": __version__,
        "parallel_read_safe": True,
        "parallel_write_safe": True,
    }
