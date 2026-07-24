# Skip text characters for text token, place those to pending buffer
# and increment current pos
from .state_inline import StateInline

# Rule to skip pure text


def text(state: StateInline, silent: bool) -> bool:
    pos = state.pos
    posMax = state.posMax

    terminator_char = state.md.inline.terminator_re.search(state.src, pos)
    pos = terminator_char.start() if terminator_char else posMax

    if pos == state.pos:
        return False

    if not silent:
        state.pending += state.src[state.pos : pos]

    state.pos = pos

    return True
