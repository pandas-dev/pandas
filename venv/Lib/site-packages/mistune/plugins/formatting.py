from typing import TYPE_CHECKING, Match, Optional

from ..helpers import PREVENT_BACKSLASH

if TYPE_CHECKING:
    from ..core import BaseRenderer, InlineState
    from ..inline_parser import InlineParser
    from ..markdown import Markdown

__all__ = ["strikethrough", "mark", "insert", "superscript", "subscript"]

SUPERSCRIPT_PATTERN = r"\^(?:" + PREVENT_BACKSLASH + r"\\\^|\S|\\ )+?\^"
SUBSCRIPT_PATTERN = r"~(?:" + PREVENT_BACKSLASH + r"\\~|\S|\\ )+?~"


def parse_strikethrough(inline: "InlineParser", m: Match[str], state: "InlineState") -> Optional[int]:
    return _parse_to_end(inline, m, state, "strikethrough", "~~")


def render_strikethrough(renderer: "BaseRenderer", text: str) -> str:
    return "<del>" + text + "</del>"


def parse_mark(inline: "InlineParser", m: Match[str], state: "InlineState") -> Optional[int]:
    return _parse_to_end(inline, m, state, "mark", "==")


def render_mark(renderer: "BaseRenderer", text: str) -> str:
    return "<mark>" + text + "</mark>"


def parse_insert(inline: "InlineParser", m: Match[str], state: "InlineState") -> Optional[int]:
    return _parse_to_end(inline, m, state, "insert", "^^")


def render_insert(renderer: "BaseRenderer", text: str) -> str:
    return "<ins>" + text + "</ins>"


def parse_superscript(inline: "InlineParser", m: Match[str], state: "InlineState") -> int:
    return _parse_script(inline, m, state, "superscript")


def render_superscript(renderer: "BaseRenderer", text: str) -> str:
    return "<sup>" + text + "</sup>"


def parse_subscript(inline: "InlineParser", m: Match[str], state: "InlineState") -> int:
    return _parse_script(inline, m, state, "subscript")


def render_subscript(renderer: "BaseRenderer", text: str) -> str:
    return "<sub>" + text + "</sub>"


def _parse_to_end(
    inline: "InlineParser",
    m: Match[str],
    state: "InlineState",
    tok_type: str,
    marker: str,
) -> Optional[int]:
    pos = m.end()
    cache_key = (id(state.src), marker)
    cache = state.formatting_no_end.get(cache_key)
    if cache is not None and cache[0] is state.src and pos <= cache[1]:
        return None

    end_pos = _find_end_marker(state.src, pos, marker)
    if end_pos is None:
        state.formatting_no_end[cache_key] = (state.src, len(state.src))
        return None
    text = state.src[pos : end_pos - 2]
    new_state = state.copy()
    new_state.src = text
    children = inline.render(new_state)
    state.append_token({"type": tok_type, "children": children})
    return end_pos


def _find_end_marker(src: str, pos: int, marker: str) -> Optional[int]:
    c = marker[0]
    marker_len = len(marker)
    end = src.find(marker, pos)
    while end != -1:
        marker_end = end + marker_len
        if marker_end < len(src) and src[marker_end] == c and not src.startswith(marker, marker_end):
            end = src.find(marker, end + 1)
            continue

        if end > pos:
            prev = src[end - 1]
            escaped_marker_before = (
                prev == c and end >= pos + 2 and src[end - 2] == "\\" and not _has_odd_backslashes(src, end - 2)
            )
            if (not prev.isspace() and prev != c) or escaped_marker_before:
                return marker_end

        end = src.find(marker, end + 1)
    return None


def _has_odd_backslashes(src: str, pos: int) -> bool:
    count = 0
    pos -= 1
    while pos >= 0 and src[pos] == "\\":
        count += 1
        pos -= 1
    return count % 2 == 1


def _parse_script(inline: "InlineParser", m: Match[str], state: "InlineState", tok_type: str) -> int:
    text = m.group(0)
    new_state = state.copy()
    new_state.src = text[1:-1].replace("\\ ", " ")
    children = inline.render(new_state)
    state.append_token({"type": tok_type, "children": children})
    return m.end()


def strikethrough(md: "Markdown") -> None:
    """A mistune plugin to support strikethrough. Spec defined by
    GitHub flavored Markdown and commonly used by many parsers:

    .. code-block:: text

        ~~This was mistaken text~~

    It will be converted into HTML:

    .. code-block:: html

        <del>This was mistaken text</del>

    :param md: Markdown instance
    """
    md.inline.register(
        "strikethrough",
        r"~~(?=[^\s~])",
        parse_strikethrough,
        before="link",
    )
    if md.renderer and md.renderer.NAME == "html":
        md.renderer.register("strikethrough", render_strikethrough)


def mark(md: "Markdown") -> None:
    """A mistune plugin to add ``<mark>`` tag. Spec defined at
    https://facelessuser.github.io/pymdown-extensions/extensions/mark/:

    .. code-block:: text

        ==mark me== ==mark \\=\\= equal==

    :param md: Markdown instance
    """
    md.inline.register(
        "mark",
        r"==(?=[^\s=])",
        parse_mark,
        before="link",
    )
    if md.renderer and md.renderer.NAME == "html":
        md.renderer.register("mark", render_mark)


def insert(md: "Markdown") -> None:
    """A mistune plugin to add ``<ins>`` tag. Spec defined at
    https://facelessuser.github.io/pymdown-extensions/extensions/caret/#insert:

    .. code-block:: text

        ^^insert me^^

    :param md: Markdown instance
    """
    md.inline.register(
        "insert",
        r"\^\^(?=[^\s\^])",
        parse_insert,
        before="link",
    )
    if md.renderer and md.renderer.NAME == "html":
        md.renderer.register("insert", render_insert)


def superscript(md: "Markdown") -> None:
    """A mistune plugin to add ``<sup>`` tag. Spec defined at
    https://pandoc.org/MANUAL.html#superscripts-and-subscripts:

    .. code-block:: text

        2^10^ is 1024.

    :param md: Markdown instance
    """
    md.inline.register("superscript", SUPERSCRIPT_PATTERN, parse_superscript, before="linebreak")
    if md.renderer and md.renderer.NAME == "html":
        md.renderer.register("superscript", render_superscript)


def subscript(md: "Markdown") -> None:
    """A mistune plugin to add ``<sub>`` tag. Spec defined at
    https://pandoc.org/MANUAL.html#superscripts-and-subscripts:

    .. code-block:: text

        H~2~O is a liquid.

    :param md: Markdown instance
    """
    md.inline.register("subscript", SUBSCRIPT_PATTERN, parse_subscript, before="linebreak")
    if md.renderer and md.renderer.NAME == "html":
        md.renderer.register("subscript", render_subscript)
