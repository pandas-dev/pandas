import re
import string
from typing import Any, Dict, Tuple, Union

from .util import escape_url

PREVENT_BACKSLASH = r"(?<!\\)(?:\\\\)*"
PUNCTUATION = r"[" + re.escape(string.punctuation) + r"]"

LINK_LABEL = r"(?:[^\\\[\]]|\\.){0,500}"

ASCII_WHITESPACE = " \t\n\r\f"

HTML_TAGNAME = r"[A-Za-z][A-Za-z0-9-]*"
HTML_ATTRIBUTES = (
    r"(?:\s+[A-Za-z_:][A-Za-z0-9_.:-]*"
    r'(?:\s*=\s*(?:[^ !"\'=<>`]+|\'[^\']*?\'|"[^\"]*?"))?)*'
)

BLOCK_TAGS = (
    "address",
    "article",
    "aside",
    "base",
    "basefont",
    "blockquote",
    "body",
    "caption",
    "center",
    "col",
    "colgroup",
    "dd",
    "details",
    "dialog",
    "dir",
    "div",
    "dl",
    "dt",
    "fieldset",
    "figcaption",
    "figure",
    "footer",
    "form",
    "frame",
    "frameset",
    "h1",
    "h2",
    "h3",
    "h4",
    "h5",
    "h6",
    "head",
    "header",
    "hr",
    "html",
    "iframe",
    "legend",
    "li",
    "link",
    "main",
    "menu",
    "menuitem",
    "meta",
    "nav",
    "noframes",
    "ol",
    "optgroup",
    "option",
    "p",
    "param",
    "section",
    "source",
    "summary",
    "table",
    "tbody",
    "td",
    "tfoot",
    "th",
    "thead",
    "title",
    "tr",
    "track",
    "ul",
)
PRE_TAGS = ("pre", "script", "style", "textarea")

_INLINE_LINK_LABEL_RE = re.compile(LINK_LABEL + r"\]")
_INLINE_SQUARE_BRACKET_RE = re.compile(PREVENT_BACKSLASH + r"[\[\]]")
_ESCAPE_CHAR_RE = re.compile(r"\\(" + PUNCTUATION + r")")


def unescape_char(text: str) -> str:
    return _ESCAPE_CHAR_RE.sub(r"\1", text)


def parse_link_text(src: str, pos: int) -> Union[Tuple[str, int], Tuple[None, int]]:
    level = 1
    found = False
    start_pos = pos

    while pos < len(src):
        m = _INLINE_SQUARE_BRACKET_RE.search(src, pos)
        if not m:
            pos = len(src)  # FIX: record we scanned to end
            break

        pos = m.end()
        marker = m.group(0)
        if marker == "]":
            level -= 1
            if level == 0:
                found = True
                break
        else:
            level += 1

    if found:
        text = src[start_pos : pos - 1]
        return text, pos
    return None, pos  # FIX: return pos instead of None


def parse_link_label(src: str, start_pos: int) -> Union[Tuple[str, int], Tuple[None, None]]:
    m = _INLINE_LINK_LABEL_RE.match(src, start_pos)
    if m:
        label = m.group(0)[:-1]
        return label, m.end()
    return None, None


def parse_link_href(src: str, start_pos: int, block: bool = False) -> Union[Tuple[str, int], Tuple[None, None]]:
    href, href_pos, _end_pos = _parse_link_href(src, start_pos, block=block)
    if href is None:
        return None, None
    assert href_pos is not None
    return href, href_pos


def _parse_link_href(
    src: str, start_pos: int, block: bool = False
) -> Tuple[Union[str, None], Union[int, None], int]:
    pos = _skip_link_start_whitespace(src, start_pos)
    if pos >= len(src):
        return None, None, pos

    if src[pos] == "<":
        href, href_pos = _parse_angle_link_href(src, pos)
        if href is None:
            return None, None, pos
        assert href_pos is not None
        return href, href_pos, href_pos
    if block and src[pos] in ASCII_WHITESPACE:
        return None, None, pos

    start = pos
    level = 0
    while pos < len(src):
        c = src[pos]
        if c in ASCII_WHITESPACE:
            break
        if c == "\x00":
            return None, None, pos
        if c == "\\" and pos + 1 < len(src) and src[pos + 1] in string.punctuation:
            pos = min(pos + 2, len(src))
            continue
        if not block:
            if c == "(":
                level += 1
            elif c == ")":
                if level == 0:
                    break
                level -= 1
        pos += 1

    if not block and level != 0:
        return None, None, pos
    return src[start:pos], pos, pos


def parse_link_title(src: str, start_pos: int, max_pos: int) -> Union[Tuple[str, int], Tuple[None, None]]:
    pos = start_pos
    if pos >= max_pos or src[pos] not in ASCII_WHITESPACE:
        return None, None

    pos = _skip_ascii_whitespace(src, pos, max_pos)
    if pos >= max_pos:
        return None, None

    opener = src[pos]
    closer = {"'": "'", '"': '"', "(": ")"}.get(opener)
    if closer is None:
        return None, None

    pos += 1
    title = []
    while pos < max_pos:
        c = src[pos]
        if c == "\x00":
            return None, None
        if c == "\\":
            if pos + 1 < max_pos:
                title.append(src[pos : pos + 2])
                pos += 2
                continue
            return None, None
        if c == closer:
            return unescape_char("".join(title)), pos + 1
        title.append(src[pos])
        pos += 1
    return None, None


def parse_link(src: str, pos: int) -> Union[Tuple[Dict[str, Any], int], Tuple[None, None]]:
    attrs, next_pos, _end_pos = parse_link_with_end(src, pos)
    if attrs is None:
        return None, None
    assert next_pos is not None
    return attrs, next_pos


def parse_link_with_end(
    src: str, pos: int
) -> Tuple[Union[Dict[str, Any], None], Union[int, None], int]:
    href, href_pos, scan_end = _parse_link_href(src, pos)
    if href is None:
        return None, None, scan_end
    assert href_pos is not None
    title, title_pos = parse_link_title(src, href_pos, len(src))
    next_pos = title_pos or href_pos
    next_pos = _skip_ascii_whitespace(src, next_pos)
    if next_pos >= len(src) or src[next_pos] != ")":
        return None, None, next_pos

    href = unescape_char(href)
    attrs = {"url": escape_url(href)}
    if title:
        attrs["title"] = title
    return attrs, next_pos + 1, next_pos + 1


def _skip_ascii_whitespace(src: str, pos: int, max_pos: Union[int, None] = None) -> int:
    if max_pos is None:
        max_pos = len(src)
    while pos < max_pos and src[pos] in ASCII_WHITESPACE:
        pos += 1
    return pos


def _skip_link_start_whitespace(src: str, pos: int) -> int:
    while pos < len(src) and src[pos] in " \t":
        pos += 1
    if pos < len(src) and src[pos] in "\n\r":
        if src[pos] == "\r" and pos + 1 < len(src) and src[pos + 1] == "\n":
            pos += 2
        else:
            pos += 1
        while pos < len(src) and src[pos] in " \t":
            pos += 1
    return pos


def _parse_angle_link_href(src: str, pos: int) -> Union[Tuple[str, int], Tuple[None, None]]:
    start = pos + 1
    pos = start
    while pos < len(src):
        c = src[pos]
        if c == ">":
            return src[start:pos], pos + 1
        if c in "<\\\n\r\x00":
            return None, None
        pos += 1
    return None, None
