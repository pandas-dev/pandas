"""Support for the `bgcl` (Apple Color Emoji background) table.

This table stores a JSON blob. We decode it to a Python object on
decompile and emit human readable JSON in the TTX via a <json> element.
On compile we re-encode the JSON to UTF-8 bytes which is what Apple code
reads via CTFontCopyTable then JSONDecoder.

On iOS 16 and later the `bgcl` payload is used as the wallpaper
background when the user selects an emoji wallpaper.

Fields:

- ``colors``: list of palette entries. Each entry is an array of four
    integers ``[R, G, B, A]``. R/G/B are 0-255. A is 0-1.

- ``emojicolors``: list of per-emoji palettes. Each item is an array of
    three sublists: primary/dominant, accent, contextual
    (names inferred). Each sublist contains integer indexes referencing
    entries in ``colors``; the runtime uses these to assemble layered
    background tints for an emoji.

- ``indexmap``: mapping (glyph identifier → palette index). The map
    maps a glyph identity to an integer index selecting an entry in
    ``emojicolors``. The font/UI uses this to pick the correct palette
    for a glyph when rendering wallpaper backgrounds.

- ``version``: integer table version used for parsing/compatibility.

Runtime usage summary:

- The system fetches the table bytes with ``CTFontCopyTable('bgcl')``,
    decodes the bytes as UTF‑8 JSON and runs the JSON through the app's
    decoder into an internal ``BgclTable`` structure. The app looks up a
    glyph's entry in ``indexmap``, retrieves the corresponding
    ``emojicolors`` palette, then converts the referenced ``colors``
    entries into color objects (normalizing channels/alpha as needed).
    This resulting color(s) drive the wallpaper background appearance for
    emoji wallpapers on supported iOS versions.

This implementation preserves the JSON payload and exposes convenience
attributes ``colors``, ``emojicolors``, ``indexmap`` and ``version``
when available.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from fontTools.misc.textTools import tostr, strjoin
from . import DefaultTable
import json

if TYPE_CHECKING:
    from fontTools.misc.xmlWriter import XMLWriter
    from fontTools.ttLib import TTFont


class table__b_g_c_l(DefaultTable.DefaultTable):
    """bgcl table: stores a JSON blob describing background palettes.

    The JSON structure typically contains the top-level keys:
      - colors: [[R,G,B,A], ...]
      - emojicolors: [[[dominant...],[accent...],[contextual...]], ...]
      - indexmap: {"glyphIndex": emojicolors_index, ...}
      - version: int
    """

    def decompile(self, data: bytes, ttFont: TTFont) -> None:
        """Store raw bytes and attempt to parse JSON.

        The JSON commonly includes palette/lookup data used at runtime to
        construct wallpaper/background colors for emoji wallpapers on
        recent iOS releases.
        """
        self.data = data
        try:
            text = tostr(data, "utf_8")
            self.json = json.loads(text)
        except Exception as e:  # keep table decompilation robust
            self.json = None
            self.ERROR = f"bgcl JSON parse error: {e!r}"
            return
        # convenient attributes
        self.colors = self.json.get("colors")
        self.emojicolors = self.json.get("emojicolors")
        self.indexmap = self.json.get("indexmap")
        self.version = self.json.get("version")

    def compile(self, ttFont: TTFont) -> bytes:
        """Encode the JSON object to UTF-8 bytes for font binary storage."""
        if getattr(self, "json", None) is None:
            # fallback to raw bytes if parsing failed earlier
            return getattr(self, "data", b"")
        # use compact representation for binary table
        return json.dumps(self.json, separators=(",", ":"), ensure_ascii=False).encode(
            "utf_8"
        )

    def toXML(self, writer: XMLWriter, ttFont: TTFont) -> None:
        """Emit pretty-printed JSON inside a <json> element for human inspection."""
        if getattr(self, "json", None) is None:
            # fallback to default hex output
            DefaultTable.DefaultTable.toXML(self, writer, ttFont)
            return
        writer.begintag("json")
        writer.newline()
        pretty = json.dumps(self.json, indent=2, ensure_ascii=False)
        writer.writecdata(pretty)
        writer.newline()
        writer.endtag("json")
        writer.newline()

    def fromXML(
        self, name: str, attrs: dict[str, str], content, ttFont: TTFont
    ) -> None:
        """Read JSON from the <json> element. `content` may be a list.

        This mirrors SVG/other table `fromXML` handlers which accept a
        list of content chunks.
        """
        if name != "json":
            # fall back to DefaultTable behavior for unknown elements
            return DefaultTable.DefaultTable.fromXML(self, name, attrs, content, ttFont)
        text = strjoin(content).strip()
        try:
            self.json = json.loads(text)
            # keep raw bytes in sync
            self.data = text.encode("utf_8")
            self.colors = self.json.get("colors")
            self.emojicolors = self.json.get("emojicolors")
            self.indexmap = self.json.get("indexmap")
            self.version = self.json.get("version")
        except Exception as e:
            # store error and fall back to raw
            self.json = None
            self.ERROR = f"bgcl JSON parse error in fromXML: {e!r}"
