import re
from typing import TYPE_CHECKING, Any, Dict, Iterable, List, Match, Optional, Tuple

from ..util import strip_end

if TYPE_CHECKING:
    from ..block_parser import BlockParser
    from ..core import BaseRenderer, BlockState
    from ..markdown import Markdown

__all__ = ["def_list"]

# https://michelf.ca/projects/php-markdown/extra/#def-list

DEF_PATTERN = r"^:[ \t]+.*(?:\n|$)"
DD_START_RE = re.compile(r"^:[ \t]+", re.M)
DD_MARKER_RE = re.compile(r"^:[ \t]")
TRIM_RE = re.compile(r"^ {0,4}", re.M)
HAS_BLANK_LINE_RE = re.compile(r"\n[ \t]*\n$")


def parse_def_list(block: "BlockParser", m: Match[str], state: "BlockState") -> Optional[int]:
    head = _get_previous_paragraph(state)
    if head is None:
        return None

    definitions, end_pos = _collect_definitions(state.src, state.cursor, state.cursor_max, head[1])
    if not definitions:
        return None

    children = list(_parse_def_item(block, head[0], definitions))
    _replace_previous_paragraph(state, children)
    return end_pos


def _parse_def_item(
    block: "BlockParser",
    head: str,
    definitions: List[Tuple[str, bool]],
) -> Iterable[Dict[str, Any]]:
    for line in head.splitlines():
        yield {
            "type": "def_list_head",
            "text": line,
        }

    for text, loose in definitions:
        children = _process_text(block, text, loose)
        yield {
            "type": "def_list_item",
            "children": children,
        }


def _process_text(block: "BlockParser", text: str, loose: bool) -> List[Any]:
    text = TRIM_RE.sub("", text)
    state = block.state_cls()
    state.process(strip_end(text))
    # use default list rules
    block.parse(state, block.list_rules)
    tokens = state.tokens
    if not loose and len(tokens) == 1 and tokens[0]["type"] == "paragraph":
        tokens[0]["type"] = "block_text"
    return tokens


def _get_previous_paragraph(state: "BlockState") -> Optional[Tuple[str, bool]]:
    if not state.tokens:
        return None

    last_token = state.tokens[-1]
    if last_token["type"] == "paragraph":
        return last_token["text"], False

    if last_token["type"] == "blank_line" and len(state.tokens) > 1:
        prev_token = state.tokens[-2]
        if prev_token["type"] == "paragraph":
            return prev_token["text"], True

    return None


def _replace_previous_paragraph(state: "BlockState", children: List[Dict[str, Any]]) -> None:
    if state.tokens[-1]["type"] == "blank_line":
        state.tokens.pop()
    state.tokens.pop()

    if state.tokens and state.tokens[-1]["type"] == "def_list":
        state.tokens[-1]["children"].extend(children)
    else:
        state.append_token(
            {
                "type": "def_list",
                "children": children,
            }
        )


def _collect_definitions(src: str, pos: int, max_pos: int, loose: bool) -> Tuple[List[Tuple[str, bool]], int]:
    definitions = []
    while pos < max_pos:
        line = _get_line(src, pos, max_pos)
        if not DD_START_RE.match(line):
            break

        start = pos
        pos += len(line)
        pos = _scan_definition_tail(src, pos, max_pos)
        # drop the ':' marker plus its trailing space or tab; a tab here would
        # otherwise survive TRIM_RE (spaces only) and make the line a code block
        text = DD_MARKER_RE.sub("  ", src[start:pos], count=1)
        definitions.append((text, loose))
        loose = bool(HAS_BLANK_LINE_RE.search(text))

    return definitions, pos


def _scan_definition_tail(src: str, pos: int, max_pos: int) -> int:
    while pos < max_pos:
        line = _get_line(src, pos, max_pos)
        if DD_START_RE.match(line):
            break

        if line.strip():
            pos += len(line)
            continue

        while pos < max_pos:
            line = _get_line(src, pos, max_pos)
            if line.strip():
                break
            pos += len(line)

        if pos >= max_pos:
            break

        line = _get_line(src, pos, max_pos)
        if DD_START_RE.match(line):
            break
        if line.startswith((" ", "\t")):
            continue
        return pos

    return pos


def _get_line(src: str, pos: int, max_pos: int) -> str:
    end = src.find("\n", pos, max_pos)
    if end == -1:
        return src[pos:max_pos]
    return src[pos : end + 1]


def render_def_list(renderer: "BaseRenderer", text: str) -> str:
    return "<dl>\n" + text + "</dl>\n"


def render_def_list_head(renderer: "BaseRenderer", text: str) -> str:
    return "<dt>" + text + "</dt>\n"


def render_def_list_item(renderer: "BaseRenderer", text: str) -> str:
    return "<dd>" + text + "</dd>\n"


def def_list(md: "Markdown") -> None:
    """A mistune plugin to support def list, spec defined at
    https://michelf.ca/projects/php-markdown/extra/#def-list

    Here is an example:

    .. code-block:: text

        Apple
        :   Pomaceous fruit of plants of the genus Malus in
            the family Rosaceae.

        Orange
        :   The fruit of an evergreen tree of the genus Citrus.

    It will be converted into HTML:

    .. code-block:: html

        <dl>
        <dt>Apple</dt>
        <dd>Pomaceous fruit of plants of the genus Malus in
        the family Rosaceae.</dd>

        <dt>Orange</dt>
        <dd>The fruit of an evergreen tree of the genus Citrus.</dd>
        </dl>

    :param md: Markdown instance
    """
    md.block.register("def_list", DEF_PATTERN, parse_def_list, before="paragraph")
    if md.renderer and md.renderer.NAME == "html":
        md.renderer.register("def_list", render_def_list)
        md.renderer.register("def_list_head", render_def_list_head)
        md.renderer.register("def_list_item", render_def_list_item)
