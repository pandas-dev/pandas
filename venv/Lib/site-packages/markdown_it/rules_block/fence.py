# fences (``` lang, ~~~ lang)
from __future__ import annotations

from collections.abc import Callable
import logging

from .state_block import StateBlock

LOGGER = logging.getLogger(__name__)


def make_fence_rule(
    *,
    markers: tuple[str, ...] = ("~", "`"),
    token_type: str = "fence",
    exact_match: bool = False,
    disallow_marker_in_info: tuple[str, ...] = ("`",),
    min_markers: int = 3,
) -> Callable[[StateBlock, int, int, bool], bool]:
    """Create a fence parsing rule with configurable options.

    :param markers: Tuple of single characters that can be used as fence markers.
    :param token_type: The token type name to emit (e.g. "fence", "colon_fence").
    :param exact_match: If True, the closing fence must have exactly the same
        number of marker characters as the opening fence (not "at least as many").
        This enables nesting of fences with different marker counts.
    :param disallow_marker_in_info: Tuple of marker characters that are not allowed
        to appear in the info string. The check only applies when the actual opening
        marker is in this tuple (e.g. a tilde fence is unaffected by ``"`"`` being
        listed). Per CommonMark, backtick fences cannot have backticks in the info
        string. Use ``()`` to disable this restriction.
    :param min_markers: Minimum number of marker characters to form a fence.
    :return: A block rule function with signature
        ``(state, startLine, endLine, silent) -> bool``.
    """

    closing_matcher: Callable[[int, int], bool]
    if exact_match:
        # closing code fence must have exactly the same number of markers as the opening one
        closing_matcher = lambda opening_len, closing_len: closing_len == opening_len  # noqa: E731
    else:
        # closing code fence must be at least as long as the opening one
        closing_matcher = lambda opening_len, closing_len: closing_len >= opening_len  # noqa: E731

    def _fence_rule(
        state: StateBlock, startLine: int, endLine: int, silent: bool
    ) -> bool:
        LOGGER.debug(
            "entering fence: %s, %s, %s, %s", state, startLine, endLine, silent
        )

        haveEndMarker = False
        pos = state.bMarks[startLine] + state.tShift[startLine]
        maximum = state.eMarks[startLine]

        if state.is_code_block(startLine):
            return False

        if pos + min_markers > maximum:
            return False

        marker = state.src[pos]

        if marker not in markers:
            return False

        # scan marker length
        mem = pos
        pos = state.skipCharsStr(pos, marker)

        length = pos - mem

        if length < min_markers:
            return False

        markup = state.src[mem:pos]
        params = state.src[pos:maximum]

        if marker in disallow_marker_in_info and marker in params:
            return False

        # Since start is found, we can report success here in validation mode
        if silent:
            return True

        # search end of block
        nextLine = startLine

        while True:
            nextLine += 1
            if nextLine >= endLine:
                # unclosed block should be autoclosed by end of document.
                # also block seems to be autoclosed by end of parent
                break

            pos = mem = state.bMarks[nextLine] + state.tShift[nextLine]
            maximum = state.eMarks[nextLine]

            if pos < maximum and state.sCount[nextLine] < state.blkIndent:
                # non-empty line with negative indent should stop the list:
                # - ```
                #  test
                break

            try:
                if state.src[pos] != marker:
                    continue
            except IndexError:
                break

            if state.is_code_block(nextLine):
                continue

            pos = state.skipCharsStr(pos, marker)

            if not closing_matcher(length, pos - mem):
                continue

            # make sure tail has spaces only
            pos = state.skipSpaces(pos)

            if pos < maximum:
                continue

            haveEndMarker = True
            # found!
            break

        # If a fence has heading spaces, they should be removed from its inner block
        length = state.sCount[startLine]

        state.line = nextLine + (1 if haveEndMarker else 0)

        token = state.push(token_type, "code", 0)
        token.info = params
        token.content = state.getLines(startLine + 1, nextLine, length, True)
        token.markup = markup
        token.map = [startLine, state.line]

        return True

    return _fence_rule


#: The default fence rule (backtick and tilde markers, CommonMark compliant).
fence = make_fence_rule()
