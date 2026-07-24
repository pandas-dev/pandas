import re
from typing import (
    TYPE_CHECKING,
    Any,
    Dict,
    List,
    Match,
    Optional,
    Tuple,
    Union,
)

if TYPE_CHECKING:
    from ..block_parser import BlockParser
    from ..core import BaseRenderer, BlockState
    from ..markdown import Markdown

# https://michelf.ca/projects/php-markdown/extra/#table

__all__ = ["table", "table_in_quote", "table_in_list"]


TABLE_PATTERN = r"^ {0,3}\|[^\n]*\|[ \t]*(?:\n|$)"
NP_TABLE_PATTERN = r"^ {0,3}\S[^\n]*\|[^\n]*(?:\n|$)"

ALIGN_CENTER = re.compile(r"^ *:-+: *$")
ALIGN_LEFT = re.compile(r"^ *:-+ *$")
ALIGN_RIGHT = re.compile(r"^ *-+: *$")
ALIGN_NONE = re.compile(r"^ *-+ *$")


def parse_table(block: "BlockParser", m: Match[str], state: "BlockState") -> Optional[int]:
    pos = m.end()
    header = _strip_pipe_table_row(m.group(0))
    if header is None:
        return None

    align_line = state.get_line(pos)
    align = _strip_pipe_table_row(align_line)
    if align is None:
        return None

    thead, aligns = _process_thead(header, align)
    if not thead:
        return _parse_invalid_pipe_table(state, pos + len(align_line))
    assert aligns is not None
    pos += len(align_line)

    rows = []
    while pos < state.cursor_max:
        line = state.get_line(pos)
        text = _strip_pipe_table_row(line)
        if text is None:
            break

        row = _process_row(text, aligns)
        if not row:
            return _parse_invalid_pipe_table(state, pos + len(line))
        rows.append(row)
        pos += len(line)

    children = [thead, {"type": "table_body", "children": rows}]
    state.append_token({"type": "table", "children": children})
    return pos


def parse_nptable(block: "BlockParser", m: Match[str], state: "BlockState") -> Optional[int]:
    pos = m.end()
    header = _strip_table_line(m.group(0))
    if header is None:
        return None

    align_line = state.get_line(pos)
    align = _strip_table_line(align_line)
    if align is None:
        return None

    thead, aligns = _process_thead(header, align)
    if not thead:
        return None
    assert aligns is not None
    pos += len(align_line)

    rows = []
    while pos < state.cursor_max:
        line = state.get_line(pos)
        text = _strip_table_line(line)
        if text is None:
            break

        row = _process_row(text, aligns)
        if not row:
            return None
        rows.append(row)
        pos += len(line)

    children = [thead, {"type": "table_body", "children": rows}]
    state.append_token({"type": "table", "children": children})
    return pos


def _process_thead(header: str, align: str) -> Union[Tuple[None, None], Tuple[Dict[str, Any], List[Optional[str]]]]:
    headers = _split_table_cells(header)
    raw_aligns = _split_table_cells(align)
    if len(headers) != len(raw_aligns):
        return None, None

    aligns: List[Optional[str]] = []
    for v in raw_aligns:
        if ALIGN_CENTER.match(v):
            aligns.append("center")
        elif ALIGN_LEFT.match(v):
            aligns.append("left")
        elif ALIGN_RIGHT.match(v):
            aligns.append("right")
        elif ALIGN_NONE.match(v) or not v.strip():
            aligns.append(None)
        else:
            # a delimiter cell must be dashes (optionally colon-flanked) or empty;
            # anything else means this is not a delimiter row, so not a table
            return None, None

    children: List[Dict[str, Any]] = [
        {"type": "table_cell", "text": text.strip(), "attrs": {"align": aligns[i], "head": True}}
        for i, text in enumerate(headers)
    ]
    thead: Dict[str, Any] = {"type": "table_head", "children": children}
    return thead, aligns


def _process_row(text: str, aligns: List[Optional[str]]) -> Optional[Dict[str, Any]]:
    cells = _split_table_cells(text)
    if len(cells) != len(aligns):
        return None

    children: List[Dict[str, Any]] = [
        {"type": "table_cell", "text": text.strip(), "attrs": {"align": aligns[i], "head": False}}
        for i, text in enumerate(cells)
    ]
    return {"type": "table_row", "children": children}


def _strip_pipe_table_row(line: str) -> Optional[str]:
    text = line.rstrip("\n").rstrip(" \t")
    if not text.startswith("|") and text.startswith((" ", "\t")):
        text = text.lstrip(" ")
    if not text.startswith("|") or not text.endswith("|"):
        return None
    return text[1:-1]


def _parse_invalid_pipe_table(state: "BlockState", pos: int) -> int:
    while pos < state.cursor_max:
        line = state.get_line(pos)
        if _strip_pipe_table_row(line) is None:
            break
        pos += len(line)
    state.add_paragraph(state.src[state.cursor : pos])
    return pos


def _strip_table_line(line: str) -> Optional[str]:
    text = line.rstrip("\n").rstrip(" \t")
    if not text or "|" not in text:
        return None
    return text


def _split_table_cells(text: str) -> List[str]:
    cells = []
    start = 0
    pos = 0
    while pos < len(text):
        if text[pos] == "|" and not _is_escaped_pipe(text, pos):
            cells.append(text[start:pos].strip())
            start = pos + 1
        pos += 1
    cells.append(text[start:].strip())
    return cells


def _is_escaped_pipe(text: str, pos: int) -> bool:
    backslashes = 0
    pos -= 1
    while pos >= 0 and text[pos] == "\\":
        backslashes += 1
        pos -= 1
    return backslashes % 2 == 1


def render_table(renderer: "BaseRenderer", text: str) -> str:
    return "<table>\n" + text + "</table>\n"


def render_table_head(renderer: "BaseRenderer", text: str) -> str:
    return "<thead>\n<tr>\n" + text + "</tr>\n</thead>\n"


def render_table_body(renderer: "BaseRenderer", text: str) -> str:
    return "<tbody>\n" + text + "</tbody>\n"


def render_table_row(renderer: "BaseRenderer", text: str) -> str:
    return "<tr>\n" + text + "</tr>\n"


def render_table_cell(renderer: "BaseRenderer", text: str, align: Optional[str] = None, head: bool = False) -> str:
    if head:
        tag = "th"
    else:
        tag = "td"

    html = "  <" + tag
    if align:
        html += ' style="text-align:' + align + '"'

    return html + ">" + text + "</" + tag + ">\n"


def table(md: "Markdown") -> None:
    """A mistune plugin to support table, spec defined at
    https://michelf.ca/projects/php-markdown/extra/#table

    Here is an example:

    .. code-block:: text

        First Header  | Second Header
        ------------- | -------------
        Content Cell  | Content Cell
        Content Cell  | Content Cell

    :param md: Markdown instance
    """
    md.block.register("table", TABLE_PATTERN, parse_table, before="paragraph")
    md.block.register("nptable", NP_TABLE_PATTERN, parse_nptable, before="paragraph")

    if md.renderer and md.renderer.NAME == "html":
        md.renderer.register("table", render_table)
        md.renderer.register("table_head", render_table_head)
        md.renderer.register("table_body", render_table_body)
        md.renderer.register("table_row", render_table_row)
        md.renderer.register("table_cell", render_table_cell)


def table_in_quote(md: "Markdown") -> None:
    """Enable table plugin in block quotes."""
    md.block.insert_rule(md.block.block_quote_rules, "table", before="paragraph")
    md.block.insert_rule(md.block.block_quote_rules, "nptable", before="paragraph")


def table_in_list(md: "Markdown") -> None:
    """Enable table plugin in list."""
    md.block.insert_rule(md.block.list_rules, "table", before="paragraph")
    md.block.insert_rule(md.block.list_rules, "nptable", before="paragraph")
