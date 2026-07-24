from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Any, Dict, List, Tuple, cast

_CHARREF_PREFIX = re.compile(r"(#[0-9]{1,7};|#[xX][0-9a-fA-F]+;|[^\t\n\f <&#;]{1,32};)")


def is_entity_boundary(left: str, right: str) -> bool:
    return left.endswith("&") and _CHARREF_PREFIX.match(right) is not None


@dataclass
class _Delimiter:
    index: int
    marker: str
    length: int
    can_open: bool
    can_close: bool
    #: original run length, used by the CommonMark multiple-of-3 rule even
    #: after the run has been partially consumed
    orig_length: int = 0
    order: int = 0


class _DelimiterIndex:
    """Track part positions without rewriting every delimiter after a splice."""

    def __init__(self, delimiters: List[_Delimiter], part_count: int) -> None:
        self._tree = [0] * (part_count + 1)
        self._next = list(range(len(delimiters) + 1))

    def current(self, delimiter: _Delimiter) -> int:
        index = delimiter.index
        total = 0
        cursor = index + 1
        while cursor:
            total += self._tree[cursor]
            cursor -= cursor & -cursor
        return index - total

    def collapse(self, closer: _Delimiter, removed: int) -> None:
        if not removed:
            return
        cursor = closer.index + 1
        while cursor < len(self._tree):
            self._tree[cursor] += removed
            cursor += cursor & -cursor

    def deactivate(self, order: int) -> None:
        self._next[order] = self._find(order + 1)

    def deactivate_range(self, delimiters: List[_Delimiter], start: int, end: int) -> None:
        order = self._find(start)
        while order < end:
            delimiters[order].length = 0
            self._next[order] = self._find(order + 1)
            order = self._find(order)

    def _find(self, order: int) -> int:
        root = order
        while self._next[root] != root:
            root = self._next[root]
        while self._next[order] != order:
            parent = self._next[order]
            self._next[order] = root
            order = parent
        return root


def finalize_emphasis_tokens(
    tokens: List[Dict[str, Any]],
    enabled: bool,
    max_depth: int,
) -> List[Dict[str, Any]]:
    if not enabled:
        return _clean_emphasis_tokens(tokens)
    if not _contains_emphasis_marker(tokens):
        return _clean_emphasis_tokens(tokens)

    parts: List[Dict[str, Any]] = []
    delimiters: List[_Delimiter] = []
    source = _emphasis_source_text(tokens)
    source_pos = 0
    for token in tokens:
        if token["type"] == "text" and token.get("_emphasis", True):
            _split_text_token(token, source, source_pos, parts, delimiters)
        else:
            parts.append(_clean_emphasis_token(token))
        source_pos += _emphasis_source_length(token)

    if _process_dense_emphasis(parts, delimiters):
        return _merge_text_tokens(parts)

    _process_emphasis_delimiters(parts, delimiters, max_depth)
    return _merge_text_tokens(parts)


def _process_dense_emphasis(parts: List[Dict[str, Any]], delimiters: List[_Delimiter]) -> bool:
    """Process a flat run such as ``*a*a*a`` without repeated list scans."""
    if len(delimiters) < 4:
        return False

    marker = delimiters[0].marker
    if marker not in ("*", "_") or len(parts) != len(delimiters) * 2:
        return False

    for index, delimiter in enumerate(delimiters):
        if (
            delimiter.marker != marker
            or delimiter.length != 1
            or delimiter.index != index * 2
            or parts[delimiter.index]["type"] != "text"
            or parts[delimiter.index + 1]["type"] != "text"
        ):
            return False

    pair_count = len(delimiters) // 2
    processed: List[Dict[str, Any]] = []
    for pair_index in range(pair_count):
        opener = delimiters[pair_index * 2]
        closer = delimiters[pair_index * 2 + 1]
        if (
            not opener.can_open
            or not closer.can_close
            or not _can_match_emphasis_delimiters(opener, closer)
            or not _has_emphasis_content(parts, opener.index + 1, closer.index)
        ):
            return False

        if pair_index:
            processed.append(parts[opener.index - 1])
        processed.append(
            {
                "type": "emphasis",
                "children": [parts[opener.index + 1]],
            }
        )

    if pair_count:
        last_close = delimiters[pair_count * 2 - 1].index
        parts[:] = processed + parts[last_close + 1 :]
    return True


def _contains_emphasis_marker(tokens: List[Dict[str, Any]]) -> bool:
    for token in tokens:
        if token["type"] == "text" and token.get("_emphasis", True):
            raw = token["raw"]
            if "*" in raw or "_" in raw:
                return True
    return False


def _clean_emphasis_tokens(tokens: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    return [_clean_emphasis_token(token) for token in tokens]


def _clean_emphasis_token(token: Dict[str, Any]) -> Dict[str, Any]:
    if "_emphasis" not in token:
        return token
    token = token.copy()
    token.pop("_emphasis", None)
    return token


def _emphasis_source_text(tokens: List[Dict[str, Any]]) -> str:
    values = []
    for token in tokens:
        if token["type"] == "text":
            values.append(token["raw"])
        elif token["type"] in ("softbreak", "linebreak"):
            values.append("\n")
        else:
            values.append("\ufffc")
    return "".join(values)


def _emphasis_source_length(token: Dict[str, Any]) -> int:
    if token["type"] == "text":
        return len(token["raw"])
    return 1


def _split_text_token(
    token: Dict[str, Any],
    source: str,
    source_start: int,
    parts: List[Dict[str, Any]],
    delimiters: List[_Delimiter],
) -> None:
    text = token["raw"]
    pos = 0
    while pos < len(text):
        if text[pos] not in "*_":
            end = _next_delimiter_run(text, pos)
            parts.append({"type": "text", "raw": text[pos:end]})
            pos = end
            continue

        marker = text[pos]
        end = pos
        while end < len(text) and text[end] == marker:
            end += 1
        length = end - pos
        absolute = source_start + pos
        can_open = _can_open_emphasis(source, absolute, length, marker)
        can_close = _can_close_emphasis(source, absolute, length, marker)
        index = len(parts)
        parts.append({"type": "text", "raw": text[pos:end]})
        if can_open or can_close:
            delimiters.append(_Delimiter(index, marker, length, can_open, can_close, length, len(delimiters)))
        pos = end


def _next_delimiter_run(text: str, pos: int) -> int:
    while pos < len(text) and text[pos] not in "*_":
        pos += 1
    return pos


def _process_emphasis_delimiters(
    parts: List[Dict[str, Any]],
    delimiters: List[_Delimiter],
    max_depth: int,
) -> None:
    index_map = _DelimiterIndex(delimiters, len(parts))
    closer_pos = 0
    openers_bottom: Dict[Tuple[str, int, bool], int] = {}
    while closer_pos < len(delimiters):
        closer = delimiters[closer_pos]
        if not closer.can_close or closer.length == 0:
            closer_pos += 1
            continue

        opener_key = (closer.marker, closer.length % 3, closer.can_open)
        opener_pos = closer_pos - 1
        opener_bottom = openers_bottom.get(opener_key, 0)
        opener = None
        while opener_pos >= opener_bottom:
            candidate = delimiters[opener_pos]
            if (
                candidate.marker == closer.marker
                and candidate.can_open
                and candidate.length > 0
                and _can_match_emphasis_delimiters(candidate, closer)
            ):
                opener = candidate
                break
            opener_pos -= 1

        if opener is None:
            openers_bottom[opener_key] = closer_pos
            closer_pos += 1
            continue

        opener_index = index_map.current(opener)
        closer_index = index_map.current(closer)
        if opener.length >= 2 and closer.length >= 2:
            use_length = 2
        else:
            use_length = 1
        if use_length == 2 and not _has_strong_enabled(parts, opener_index, closer_index):
            use_length = 1
        if use_length == 1 and not _has_emphasis_enabled(parts, opener_index, closer_index):
            closer_pos += 1
            continue
        if not _has_emphasis_content(parts, opener_index + 1, closer_index):
            closer_pos += 1
            continue

        opener_text = parts[opener_index]
        closer_text = parts[closer_index]
        if opener_text["type"] != "text" or closer_text["type"] != "text":
            closer_pos += 1
            continue

        children = parts[opener_index + 1 : closer_index]
        if max_depth > 0 and _emphasis_depth(children) >= max_depth:
            closer_pos += 1
            continue
        opener_text["raw"] = opener_text["raw"][:-use_length]
        closer_text["raw"] = closer_text["raw"][use_length:]
        if use_length == 2:
            node = {"type": "strong", "children": children}
        else:
            node = {"type": "emphasis", "children": children}

        old_closer_index = closer_index
        parts[opener_index + 1 : old_closer_index] = [node]

        removed = old_closer_index - opener_index - 2
        if removed:
            index_map.deactivate_range(delimiters, opener_pos + 1, closer_pos)
            index_map.collapse(closer, removed)

        opener.length -= use_length
        closer.length -= use_length
        if opener.length == 0:
            opener.can_open = False
            index_map.deactivate(opener.order)
        if closer.length == 0:
            closer.can_close = False
            index_map.deactivate(closer.order)

        if opener.can_open or closer.can_close:
            closer_pos = max(opener_pos, openers_bottom.get(opener_key, 0))
        else:
            closer_pos += 1


def _emphasis_depth(tokens: List[Dict[str, Any]]) -> int:
    max_depth = 0
    stack = [(token, 0) for token in tokens]
    while stack:
        token, depth = stack.pop()
        token_type = token["type"]
        if token_type in ("emphasis", "strong"):
            depth += 1
            if depth > max_depth:
                max_depth = depth
        for child in token.get("children", ()):
            stack.append((child, depth))
    return max_depth


def _has_strong_enabled(parts: List[Dict[str, Any]], opener_index: int, closer_index: int) -> bool:
    return len(_text_raw(parts[opener_index])) >= 2 and len(_text_raw(parts[closer_index])) >= 2


def _has_emphasis_enabled(parts: List[Dict[str, Any]], opener_index: int, closer_index: int) -> bool:
    return bool(_text_raw(parts[opener_index]) and _text_raw(parts[closer_index]))


def _text_raw(token: Dict[str, Any]) -> str:
    if token["type"] == "text":
        return cast(str, token["raw"])
    return ""


def _has_emphasis_content(parts: List[Dict[str, Any]], start: int, end: int) -> bool:
    for part in parts[start:end]:
        if part["type"] != "text" or part["raw"] != "":
            return True
    return False


def _can_match_emphasis_delimiters(opener: _Delimiter, closer: _Delimiter) -> bool:
    if opener.can_close or closer.can_open:
        open_len = opener.orig_length
        close_len = closer.orig_length
        return (open_len + close_len) % 3 != 0 or open_len % 3 == 0 and close_len % 3 == 0
    return True


def _can_open_emphasis(text: str, start: int, size: int, marker: str) -> bool:
    previous = text[start - 1] if start > 0 else "\n"
    next_char = text[start + size] if start + size < len(text) else "\n"
    if marker == "_" and previous.isalnum() and next_char.isalnum():
        return False
    if next_char.isspace():
        return False
    if _is_punctuation(next_char) and not previous.isspace() and not _is_punctuation(previous):
        return False
    return True


def _can_close_emphasis(text: str, start: int, size: int, marker: str) -> bool:
    previous = text[start - 1] if start > 0 else "\n"
    next_char = text[start + size] if start + size < len(text) else "\n"
    if marker == "_" and previous.isalnum() and next_char.isalnum():
        return False
    if previous.isspace():
        return False
    if _is_punctuation(previous) and not next_char.isspace() and not _is_punctuation(next_char):
        return False
    return True


def _is_punctuation(char: str) -> bool:
    return not char.isspace() and not char.isalnum()


def _merge_text_tokens(tokens: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    result: List[Dict[str, Any]] = []
    for token in tokens:
        if token["type"] == "text" and token["raw"] == "":
            continue
        if token["type"] == "text" and result and result[-1]["type"] == "text":
            if not is_entity_boundary(result[-1]["raw"], token["raw"]):
                result[-1]["raw"] += token["raw"]
                continue
        result.append(_clean_emphasis_token(token))
    return result
