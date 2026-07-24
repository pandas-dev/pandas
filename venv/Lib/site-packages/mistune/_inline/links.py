from __future__ import annotations

from bisect import bisect_left
from typing import TYPE_CHECKING, Dict, List, Match, Optional, Tuple

from ..core import InlineState
from ..helpers import (
    parse_link as parse_link_destination,
    parse_link_label,
    parse_link_with_end,
)
from ..util import unikey

if TYPE_CHECKING:
    from ..inline_parser import InlineParser


def parse_link(inline: "InlineParser", m: Match[str], state: InlineState) -> Optional[int]:
    pos = m.end()

    marker = m.group(0)
    is_image = marker[0] == "!"
    if is_image and inline.max_image_depth > 0 and state.image_depth >= inline.max_image_depth:
        state.append_token({"type": "text", "raw": marker + state.src[pos:]})
        return len(state.src)
    if not is_image and state.in_link:
        state.append_token({"type": "text", "raw": marker})
        return pos
    if not is_image and pos <= state.no_link_before:
        state.append_token({"type": "text", "raw": marker})
        return pos
    if is_image and pos <= state.no_image_before:
        state.append_token({"type": "text", "raw": marker})
        return pos

    text = None
    text_start = pos
    text_end = pos
    label, end_pos = parse_link_label(state.src, pos)
    if label is None:
        if pos <= state.no_close_bracket_before:
            state.append_token({"type": "text", "raw": marker})
            return pos
        close_pos = find_closing_bracket(state, pos)
        if close_pos is None:
            if len(state.src) > state.no_close_bracket_before:
                state.no_close_bracket_before = len(state.src)
            return None
        text_start = pos
        text_end = close_pos
        end_pos = close_pos + 1

    assert end_pos is not None

    if label is not None:
        text = label
        text_start = pos
        text_end = end_pos - 1

    body_end_pos = end_pos

    has_nested_link = not is_image and label_contains_link(state, text_start, text_end)
    if has_nested_link:
        return None
    if end_pos >= len(state.src) and label is None:
        mark_no_link_before(state, body_end_pos)
        return None

    if not is_image:
        rules = ["codespan", "prec_auto_link", "prec_inline_html"]
        prec_pos = inline.precedence_scan(m, state, end_pos, rules)
        if prec_pos:
            return prec_pos

    if end_pos < len(state.src):
        char = state.src[end_pos]
        if char == "(":
            attrs, pos2, scan_end = parse_link_with_end(state.src, end_pos + 1)
            if pos2:
                if text is None:
                    text = state.src[text_start:text_end]
                token = build_link_token(inline, is_image, text, attrs, state)
                state.append_token(token)
                return pos2
            if scan_end > body_end_pos:
                if is_image:
                    mark_no_image_before(state, scan_end)
                else:
                    mark_no_link_before(state, scan_end)
        elif char == "[":
            label2, pos2 = parse_link_label(state.src, end_pos + 1)
            if pos2:
                end_pos = pos2
                if label2:
                    label = label2

    if label is None:
        ref_links = state.env.get("ref_links")
        if not ref_links:
            mark_no_link_before(state, body_end_pos)
            return None
        if text is None:
            text = state.src[text_start:text_end]
        label = text
    ref_links = state.env.get("ref_links")
    if not ref_links:
        mark_no_link_before(state, body_end_pos)
        return None

    key = unikey(label)
    env = ref_links.get(key)
    if env:
        if text is None:
            text = state.src[text_start:text_end]
        attrs = {"url": env["url"], "title": env.get("title")}
        token = build_link_token(inline, is_image, text, attrs, state)
        token["ref"] = key
        token["label"] = label
        state.append_token(token)
        return end_pos
    mark_no_link_before(state, body_end_pos)
    return None


def build_link_token(
    inline: "InlineParser",
    is_image: bool,
    text: str,
    attrs: Optional[Dict[str, object]],
    state: InlineState,
) -> Dict[str, object]:
    new_state = state.copy()
    new_state.src = text
    if is_image:
        new_state.in_image = True
        new_state.image_depth += 1
        return {
            "type": "image",
            "children": inline.render(new_state),
            "attrs": attrs,
        }
    new_state.in_link = True
    return {
        "type": "link",
        "children": inline.render(new_state),
        "attrs": attrs,
    }


def mark_no_link_before(state: InlineState, end_pos: int) -> None:
    if end_pos > state.no_link_before:
        state.no_link_before = end_pos


def mark_no_image_before(state: InlineState, end_pos: int) -> None:
    if end_pos > state.no_image_before:
        state.no_image_before = end_pos


def find_closing_bracket(state: InlineState, pos: int) -> Optional[int]:
    return get_closing_bracket_map(state).get(pos)


def label_contains_link(state: InlineState, start: int, end: int) -> bool:
    if start >= end:
        return False

    starts, suffix_min_ends = get_link_range_index(state)
    index = bisect_left(starts, start)
    return index < len(starts) and starts[index] < end and suffix_min_ends[index] <= end


def get_link_range_index(state: InlineState) -> Tuple[List[int], List[int]]:
    cache = state.link_ranges.get(id(state.src))
    if cache is not None and cache[0] is state.src:
        return cache[1], cache[2]

    pairs = get_closing_bracket_map(state)
    ranges: List[Tuple[int, int]] = []
    for label_start, close_pos in pairs.items():
        opener = label_start - 1
        if opener > 0 and state.src[opener - 1] == "!":
            continue
        link_end = find_link_range_end(state.src, label_start, close_pos, state)
        if link_end is not None:
            ranges.append((opener, link_end))

    ranges.sort()
    starts = [start for start, _end in ranges]
    suffix_min_ends = [0] * len(ranges)
    min_end = len(state.src) + 1
    for index in range(len(ranges) - 1, -1, -1):
        end = ranges[index][1]
        if end < min_end:
            min_end = end
        suffix_min_ends[index] = min_end

    state.link_ranges[id(state.src)] = (state.src, starts, suffix_min_ends)
    return starts, suffix_min_ends


def get_closing_bracket_map(state: InlineState) -> Dict[int, int]:
    cache = state.link_brackets.get(id(state.src))
    if cache is not None and cache[0] is state.src:
        return cache[1]

    pairs = build_closing_bracket_map(state.src)
    state.link_brackets[id(state.src)] = (state.src, pairs)
    return pairs


def find_link_range_end(src: str, label_start: int, close_pos: int, state: InlineState) -> Optional[int]:
    end_pos = close_pos + 1
    if end_pos < len(src):
        marker = src[end_pos]
        if marker == "(":
            _attrs, new_pos = parse_link_destination(src, end_pos + 1)
            return new_pos

        if marker == "[":
            label, new_pos = parse_link_label(src, end_pos + 1)
            if not new_pos:
                return None
            ref_label = label or src[label_start:close_pos]
            ref_links = state.env.get("ref_links")
            if ref_links and unikey(ref_label) in ref_links:
                return new_pos
            return None

    ref_links = state.env.get("ref_links")
    if ref_links and unikey(src[label_start:close_pos]) in ref_links:
        return end_pos
    return None


def build_closing_bracket_map(src: str) -> Dict[int, int]:
    pairs: Dict[int, int] = {}
    stack: List[int] = []
    pos = 0
    while pos < len(src):
        char = src[pos]
        if char == "\\":
            pos += 2
            continue
        if char == "[":
            stack.append(pos + 1)
        elif char == "]" and stack:
            pairs[stack.pop()] = pos
        pos += 1
    return pairs
