"""Token-related utilities"""

# Copyright (c) IPython Development Team.
# Distributed under the terms of the Modified BSD License.
from __future__ import annotations

import itertools
import tokenize
from io import StringIO
from keyword import iskeyword
from tokenize import TokenInfo
from typing import NamedTuple
from collections.abc import Generator


class Token(NamedTuple):
    token: int
    text: str
    start: int
    end: int
    line: str


def generate_tokens(readline) -> Generator[TokenInfo, None, None]:
    """wrap generate_tkens to catch EOF errors"""
    try:
        yield from tokenize.generate_tokens(readline)
    except tokenize.TokenError:
        # catch EOF error
        return


def generate_tokens_catch_errors(
    readline, extra_errors_to_catch: list[str] | None = None
):
    default_errors_to_catch = [
        "unterminated string literal",
        "invalid non-printable character",
        "after line continuation character",
    ]
    assert extra_errors_to_catch is None or isinstance(extra_errors_to_catch, list)
    errors_to_catch = default_errors_to_catch + (extra_errors_to_catch or [])

    tokens: list[TokenInfo] = []
    try:
        for token in tokenize.generate_tokens(readline):
            tokens.append(token)
            yield token
    except tokenize.TokenError as exc:
        if any(error in exc.args[0] for error in errors_to_catch):
            if tokens:
                start = tokens[-1].start[0], tokens[-1].end[0]
                end = start
                line = tokens[-1].line
            else:
                start = end = (1, 0)
                line = ""
            yield TokenInfo(tokenize.ERRORTOKEN, "", start, end, line)
        else:
            # Catch EOF
            raise


def line_at_cursor(cell: str, cursor_pos: int = 0) -> tuple[str, int]:
    """Return the line in a cell at a given cursor position

    Used for calling line-based APIs that don't support multi-line input, yet.

    Parameters
    ----------
    cell : str
        multiline block of text
    cursor_pos : integer
        the cursor position

    Returns
    -------
    (line, offset): (string, integer)
        The line with the current cursor, and the character offset of the start of the line.
    """
    offset = 0
    lines = cell.splitlines(True)
    for line in lines:
        next_offset = offset + len(line)
        if not line.endswith("\n"):
            # If the last line doesn't have a trailing newline, treat it as if
            # it does so that the cursor at the end of the line still counts
            # as being on that line.
            next_offset += 1
        if next_offset > cursor_pos:
            break
        offset = next_offset
    else:
        line = ""
    return line, offset


def token_at_cursor(cell: str, cursor_pos: int = 0) -> str:
    """Get the token at a given cursor

    Used for introspection.

    Function calls are prioritized, so the token for the callable will be returned
    if the cursor is anywhere inside the call.

    Parameters
    ----------
    cell : str
        A block of Python code
    cursor_pos : int
        The location of the cursor in the block where the token should be found
    """
    names: list[str] = []
    call_names: list[str] = []
    closing_call_name: str | None = None
    most_recent_outer_name: str | None = None

    offsets = {1: 0}  # lines start at 1
    intersects_with_cursor = False
    cur_token_is_name = False
    tokens: list[Token | None] = [
        Token(*tup) for tup in generate_tokens(StringIO(cell).readline)
    ]
    if not tokens:
        return ""
    for prev_tok, (tok, next_tok) in zip(
        [None] + tokens, itertools.pairwise(tokens + [None])
    ):
        # token, text, start, end, line = tup
        start_line, start_col = tok.start
        end_line, end_col = tok.end
        if end_line + 1 not in offsets:
            # keep track of offsets for each line
            lines = tok.line.splitlines(True)
            for lineno, line in enumerate(lines, start_line + 1):
                if lineno not in offsets:
                    offsets[lineno] = offsets[lineno - 1] + len(line)

        closing_call_name = None

        offset = offsets[start_line]
        if offset + start_col > cursor_pos:
            # current token starts after the cursor,
            # don't consume it
            break

        if cur_token_is_name := tok.token == tokenize.NAME and not iskeyword(tok.text):
            if (
                names
                and prev_tok
                and prev_tok.token == tokenize.OP
                and prev_tok.text == "."
            ):
                names[-1] = "%s.%s" % (names[-1], tok.text)
            else:
                names.append(tok.text)
            if (
                next_tok is not None
                and next_tok.token == tokenize.OP
                and next_tok.text == "="
            ):
                # don't inspect the lhs of an assignment
                names.pop(-1)
                cur_token_is_name = False
            if not call_names:
                most_recent_outer_name = names[-1] if names else None
        elif tok.token == tokenize.OP:
            if tok.text == "(" and names:
                # if we are inside a function call, inspect the function
                call_names.append(names[-1])
            elif tok.text == ")" and call_names:
                # keep track of the most recently popped call_name from the stack
                closing_call_name = call_names.pop(-1)

        if offsets[end_line] + end_col > cursor_pos:
            # we found the cursor, stop reading
            # if the current token intersects directly, use it instead of the call token
            intersects_with_cursor = offsets[start_line] + start_col <= cursor_pos
            break

    if cur_token_is_name and intersects_with_cursor:
        return names[-1]
    # if the cursor isn't directly over a name token, use the most recent
    # call name if we can find one
    elif closing_call_name:
        # if we're on a ")", use the most recently popped call name
        return closing_call_name
    elif call_names:
        # otherwise, look for the most recent call name in the stack
        return call_names[-1]
    elif most_recent_outer_name:
        # if we've popped all the call names, use the most recently-seen
        # outer name
        return most_recent_outer_name
    elif names:
        # failing that, use the most recently seen name
        return names[-1]
    else:
        # give up
        return ""
