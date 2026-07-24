"""because list is complex, split list parser in a new file"""

from __future__ import annotations

import re
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any, Iterable, Optional, Match, Pattern, cast
from .util import strip_end

if TYPE_CHECKING:
    from .block_parser import BlockParser
    from .core import BlockState

LIST_PATTERN = (
    r"^(?P<list_1> {0,3})"
    r"(?P<list_2>[\*\+-]|\d{1,9}[.)])"
    r"(?P<list_3>[ \t]*|[ \t].+)$"
)

_LINE_HAS_TEXT = re.compile(r"(\s*)\S")


@dataclass
class _ListMarker:
    spaces: str
    marker: str
    text: str

    @property
    def leading_width(self) -> int:
        return len(self.spaces) + len(self.marker)

    @property
    def bullet(self) -> str:
        return self.marker[-1]

    @property
    def ordered(self) -> bool:
        return len(self.marker) > 1


@dataclass
class _ListItemLines:
    src: str
    next_item: Optional[_ListMarker] = None
    loose: bool = False
    end_pos: Optional[int] = None
    token_index: Optional[int] = None


def parse_list(block: "BlockParser", m: Match[str], state: "BlockState") -> int:
    """Parse tokens for ordered and unordered list."""
    item = _create_list_marker(m, "list")
    text = item.text
    if not text.strip():
        # Example 285
        # an empty list item cannot interrupt a paragraph
        end_pos = state.append_paragraph()
        if end_pos:
            return end_pos

    marker = item.marker
    depth = state.depth()
    token: dict[str, Any] = {
        "type": "list",
        "children": [],
        "tight": True,
        "bullet": item.bullet,
        "attrs": {
            "depth": depth,
            "ordered": item.ordered,
        },
    }
    if item.ordered:
        start = int(marker[:-1])
        if start != 1:
            # Example 304
            # we allow only lists starting with 1 to interrupt paragraphs
            end_pos = state.append_paragraph()
            if end_pos:
                return end_pos
            token["attrs"]["start"] = start

    state.cursor = m.end() + 1
    item_or_none: Optional[_ListMarker] = item

    if depth >= block.max_nested_level - 1:
        # At the nesting limit, stop descending into any further container
        # blocks. Trimming only "list" still allowed lists and block quotes to
        # recurse into each other without bound (RecursionError).
        rules = [rule for rule in block.list_rules if rule not in ("list", "block_quote")]
    else:
        rules = block.list_rules

    bullet = _get_list_bullet(item.bullet)
    while item_or_none:
        item_or_none = _parse_list_item(block, bullet, item_or_none, token, state, rules)

    end_pos = cast(Optional[int], token.pop("_end_pos", None))
    _transform_tight_list(token)
    if end_pos:
        index = cast(int, token.pop("_tok_index"))
        state.tokens.insert(index, token)
        return end_pos

    state.append_token(token)
    return state.cursor


def _transform_tight_list(token: dict[str, Any]) -> None:
    if token["tight"]:
        # reset tight list item
        for list_item in token["children"]:
            for tok in list_item["children"]:
                if tok["type"] == "paragraph":
                    tok["type"] = "block_text"
                elif tok["type"] == "list":
                    _transform_tight_list(tok)


def _parse_list_item(
    block: "BlockParser",
    bullet: str,
    item: _ListMarker,
    token: dict[str, Any],
    state: "BlockState",
    rules: list[str],
) -> _ListMarker | None:
    text = item.text
    leading_width = item.leading_width
    text, continue_width = _compile_continue_width(text, leading_width)
    list_item_re = re.compile(_compile_list_item_pattern(bullet, leading_width))
    break_sc = _compile_list_break_sc(block, leading_width)

    lines = _collect_list_item_lines(block, list_item_re, break_sc, state, text, continue_width)
    if lines.loose:
        token["tight"] = False
    if lines.end_pos is not None:
        token["_tok_index"] = lines.token_index
        token["_end_pos"] = lines.end_pos

    child = state.child_state(_build_list_item_source(text, lines.src, continue_width))

    block.parse(child, rules)

    if token["tight"] and _is_loose_list(child.tokens):
        token["tight"] = False

    token["children"].append(
        {
            "type": "list_item",
            "children": child.tokens,
        }
    )
    if lines.next_item:
        return lines.next_item

    return None


def _collect_list_item_lines(
    block: "BlockParser",
    list_item_re: Pattern[str],
    break_sc: Pattern[str],
    state: "BlockState",
    text: str,
    continue_width: int,
) -> _ListItemLines:
    src = ""
    next_item = None
    prev_blank_line = False
    while state.cursor < state.cursor_max:
        raw_line = state.get_line(state.cursor)
        next_pos = state.cursor + len(raw_line)
        if block.BLANK_LINE.match(raw_line):
            src += "\n"
            prev_blank_line = True
            state.cursor = next_pos
            continue

        has_continuation = _has_continuation_indent(raw_line, continue_width)
        if has_continuation:
            if prev_blank_line and not text and not src.strip():
                # Example 280
                # A list item can begin with at most one blank line
                break

            src += raw_line
            prev_blank_line = False
            state.cursor = next_pos
            continue

        line = _expand_leading_tabs(raw_line)
        line_break = _match_list_item_break(list_item_re, break_sc, state, line)
        if line_break:
            tok_type, m = line_break
            if tok_type == "list_item":
                next_item = _create_list_marker(m, "listitem")
                state.cursor = next_pos
                return _ListItemLines(src, next_item=next_item, loose=prev_blank_line)

            if tok_type == "list":
                break

            tok_index = len(state.tokens)
            end_pos = block.parse_method(m, state)
            if end_pos:
                return _ListItemLines(src, end_pos=end_pos, token_index=tok_index)

        if prev_blank_line and not has_continuation:
            # not a continue line, and previous line is blank
            break

        src += raw_line
        state.cursor = next_pos

    return _ListItemLines(src)


def _create_list_marker(m: Match[str], prefix: str) -> _ListMarker:
    return _ListMarker(
        spaces=m.group(prefix + "_1"),
        marker=m.group(prefix + "_2"),
        text=m.group(prefix + "_3"),
    )


def _build_list_item_source(text: str, src: str, continue_width: int) -> str:
    text += _clean_list_item_text(src, continue_width)
    return strip_end(text)


def _compile_list_break_sc(block: "BlockParser", leading_width: int) -> Pattern[str]:
    pairs = [(name, block.specification[name]) for name in _get_list_break_rules(block)]
    if leading_width < 3:
        # Relax the leading indent bound only. Matching on a bare "3" would
        # rewrite the first quantifier of any rule that has no indent prefix --
        # e.g. a fenced directive's "{3,}" marker run.
        _repl = " {0,%d}" % leading_width
        pairs = [(n, p.replace(" {0,3}", _repl, 1)) for n, p in pairs]

    regex = "|".join(r"(?P<%s>(?<=\n)%s)" % pair for pair in pairs)
    return re.compile(regex, re.M)


def _get_list_break_rules(block: "BlockParser") -> list[str]:
    rules = [
        "thematic_break",
        "fenced_code",
        "atx_heading",
        "block_quote",
        "block_html",
        "list",
    ]
    if "fenced_directive" in block.specification:
        rules.insert(1, "fenced_directive")
    return rules


def _match_list_item_break(
    list_item_re: Pattern[str],
    break_sc: Pattern[str],
    state: "BlockState",
    line: str,
) -> tuple[str, Match[str]] | None:
    m = break_sc.match(state.src, state.cursor)
    if m and m.lastgroup == "thematic_break":
        return "thematic_break", m

    m2 = list_item_re.match(line)
    if m2:
        return "list_item", m2

    if m:
        tok_type = m.lastgroup
        assert tok_type is not None
        return tok_type, m
    return None


def _get_list_bullet(c: str) -> str:
    if c == ".":
        bullet = r"\d{0,9}\."
    elif c == ")":
        bullet = r"\d{0,9}\)"
    elif c == "*":
        bullet = r"\*"
    elif c == "+":
        bullet = r"\+"
    else:
        bullet = "-"
    return bullet


def _compile_list_item_pattern(bullet: str, leading_width: int) -> str:
    if leading_width > 3:
        leading_width = 3
    return (
        r"^(?P<listitem_1> {0," + str(leading_width) + "})"
        r"(?P<listitem_2>" + bullet + ")"
        r"(?P<listitem_3>[ \t]*|[ \t][^\n]+)$"
    )


def _compile_continue_width(text: str, leading_width: int) -> tuple[str, int]:
    text = _expand_leading_tabs(text, leading_width)

    m2 = _LINE_HAS_TEXT.match(text)
    if m2:
        # indent code, startswith 5 spaces
        indent = _count_indent(text)
        if indent >= 5:
            space_width = 1
        else:
            space_width = indent

        text = text[space_width:] + "\n"
    else:
        space_width = 1
        text = ""

    continue_width = leading_width + space_width
    return text, continue_width


def _clean_list_item_text(src: str, continue_width: int) -> str:
    rv = []
    lines = src.split("\n")
    for line in lines:
        if _has_continuation_indent(line, continue_width):
            rv.append(_strip_continuation_indent(line, continue_width))
        else:
            rv.append(_expand_leading_tabs(line))

    return "\n".join(rv)


def _has_continuation_indent(line: str, columns: int) -> bool:
    return _count_indent(line) >= columns


def _strip_continuation_indent(line: str, columns: int) -> str:
    expanded = _expand_leading_tabs(line)
    if len(expanded) >= columns:
        return expanded[columns:]
    return ""


def _expand_leading_tabs(line: str, start_column: int = 0) -> str:
    column = start_column
    parts = []
    index = 0
    while index < len(line):
        c = line[index]
        if c == " ":
            parts.append(" ")
            column += 1
        elif c == "\t":
            size = 4 - column % 4
            parts.append(" " * size)
            column += size
        else:
            break
        index += 1
    return "".join(parts) + line[index:]


def _count_indent(text: str) -> int:
    column = 0
    for c in text:
        if c == " ":
            column += 1
        elif c == "\t":
            column += 4 - column % 4
        else:
            break
    return column


def _is_loose_list(tokens: Iterable[dict[str, Any]]) -> bool:
    paragraph_count = 0
    for tok in tokens:
        if tok["type"] == "blank_line":
            return True
        if tok["type"] == "paragraph":
            paragraph_count += 1
            if paragraph_count > 1:
                return True
    return False
