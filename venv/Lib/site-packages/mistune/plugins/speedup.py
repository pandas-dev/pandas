from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from ..markdown import Markdown

__all__ = ["speedup"]


def speedup(md: "Markdown") -> None:
    """Compatibility plugin for the former parser speedups.

    The paragraph and inline text fast paths are now part of the core parsers,
    so installing this plugin intentionally leaves the Markdown instance
    unchanged.
    """
    return None
