import re
import types
from typing import TYPE_CHECKING, List, Match, Tuple

from ..helpers import PREVENT_BACKSLASH
from ..util import escape

if TYPE_CHECKING:
    from ..block_parser import BlockParser
    from ..core import BaseRenderer, BlockState, InlineState
    from ..inline_parser import InlineParser
    from ..markdown import Markdown

__all__ = ["abbr"]
TextSegment = Tuple[str, bool]

# https://michelf.ca/projects/php-markdown/extra/#abbr
REF_ABBR = (
    r"^ {0,3}\*\[(?P<abbr_key>[^\]]+)" + PREVENT_BACKSLASH + r"\]:"
    r"(?P<abbr_text>(?:[ \t]*\n(?: {3,}|\t)[^\n]+)|(?:[^\n]*))$"
)


def parse_ref_abbr(block: "BlockParser", m: Match[str], state: "BlockState") -> int:
    ref = state.env.get("ref_abbrs")
    if not ref:
        ref = {}
    key = m.group("abbr_key")
    text = m.group("abbr_text")
    ref[key] = text.strip()
    state.env["ref_abbrs"] = ref
    # abbr definition can split paragraph
    state.append_token({"type": "blank_line"})
    return m.end() + 1


def _append_text(
    inline: "InlineParser",
    text: str,
    state: "InlineState",
    parse_emphasis: bool,
) -> None:
    type(inline).process_text(inline, text, state, parse_emphasis=parse_emphasis)


def _append_text_segments(
    inline: "InlineParser",
    segments: List[TextSegment],
    state: "InlineState",
    start: int,
    end: int,
) -> None:
    offset = 0
    for value, parse_emphasis in segments:
        next_offset = offset + len(value)
        if start < next_offset and end > offset:
            part_start = max(start - offset, 0)
            part_end = min(end - offset, len(value))
            _append_text(inline, value[part_start:part_end], state, parse_emphasis)
        offset = next_offset


def process_text(
    inline: "InlineParser",
    text: str,
    state: "InlineState",
    parse_emphasis: bool = True,
) -> None:
    ref = state.env.get("ref_abbrs")
    if not ref:
        return _append_text(inline, text, state, parse_emphasis)

    segments: List[TextSegment] = []
    if state.tokens:
        last = state.tokens[-1]
        if last["type"] == "text":
            state.tokens.pop()
            segments.append((last["raw"], last.get("_emphasis", True)))
    segments.append((text, parse_emphasis))
    text = "".join(value for value, _ in segments)

    abbrs_re = state.env.get("abbrs_re")
    if not abbrs_re:
        abbrs_re = re.compile(r"|".join(re.escape(k) for k in ref.keys()))
        state.env["abbrs_re"] = abbrs_re

    pos = 0
    while pos < len(text):
        m = abbrs_re.search(text, pos)
        if not m:
            break

        end_pos = m.start()
        if end_pos > pos:
            _append_text_segments(inline, segments, state, pos, end_pos)

        label = m.group(0)
        state.append_token(
            {"type": "abbr", "children": [{"type": "text", "raw": label}], "attrs": {"title": ref[label]}}
        )
        pos = m.end()

    if pos == 0:
        # special case, just pure text
        _append_text_segments(inline, segments, state, 0, len(text))
    elif pos < len(text):
        _append_text_segments(inline, segments, state, pos, len(text))


def render_abbr(renderer: "BaseRenderer", text: str, title: str) -> str:
    if not title:
        return "<abbr>" + text + "</abbr>"
    return '<abbr title="' + escape(title) + '">' + text + "</abbr>"


def abbr(md: "Markdown") -> None:
    """A mistune plugin to support abbreviations, spec defined at
    https://michelf.ca/projects/php-markdown/extra/#abbr

    Here is an example:

    .. code-block:: text

        The HTML specification
        is maintained by the W3C.

        *[HTML]: Hyper Text Markup Language
        *[W3C]:  World Wide Web Consortium

    It will be converted into HTML:

    .. code-block:: html

        The <abbr title="Hyper Text Markup Language">HTML</abbr> specification
        is maintained by the <abbr title="World Wide Web Consortium">W3C</abbr>.

    :param md: Markdown instance
    """
    md.block.register("ref_abbr", REF_ABBR, parse_ref_abbr, before="paragraph")
    # replace process_text
    md.inline.process_text = types.MethodType(process_text, md.inline)  # type: ignore[method-assign]
    if md.renderer and md.renderer.NAME == "html":
        md.renderer.register("abbr", render_abbr)
