__all__ = ("commonmark", "default", "gfm_like", "gfm_like2", "js_default", "zero")

from ..utils import PresetType
from . import commonmark, default, zero

js_default = default


class gfm_like:  # noqa: N801
    """GitHub Flavoured Markdown (GFM) like.

    This adds the linkify, table and strikethrough components to CommmonMark.

    Note, it lacks task-list items and raw HTML filtering,
    to meet the the full GFM specification
    (see https://github.github.com/gfm/#autolinks-extension-).
    """

    @staticmethod
    def make() -> PresetType:
        config = commonmark.make()
        config["components"]["core"]["rules"].append("linkify")
        config["components"]["block"]["rules"].append("table")
        config["components"]["inline"]["rules"].extend(["strikethrough", "linkify"])
        config["components"]["inline"]["rules2"].append("strikethrough")
        config["options"]["linkify"] = True
        config["options"]["html"] = True
        return config


class gfm_like2:  # noqa: N801
    """GitHub Flavoured Markdown (GFM) like, extended.

    Builds on ``gfm-like`` and additionally enables:

    - Task lists (``- [x] done``)
    - Alerts (``> [!NOTE]``)
    - Single-tilde strikethrough (``~text~`` in addition to ``~~text~~``)
    """

    @staticmethod
    def make() -> PresetType:
        config = gfm_like.make()
        config["options"]["tasklists"] = True
        config["options"]["tasklists_editable"] = False
        config["options"]["alerts"] = True
        config["options"]["strikethrough_single_tilde"] = True
        return config
