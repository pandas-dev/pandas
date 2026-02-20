r"""
Terminal escape sequence patterns.

This module provides regex patterns for matching terminal escape sequences. All patterns match
sequences that begin with ESC (``\x1b``). Before calling re.match with these patterns, callers
should first check that the character at the current position is ESC for optimal performance.
"""
# std imports
import re

# Zero-width escape sequences (SGR, OSC, CSI, etc.). This table, like INDETERMINATE_EFFECT_SEQUENCE,
# originated from the 'blessed' library.
ZERO_WIDTH_PATTERN = re.compile(
    # CSI sequences
    r'\x1b\[[\x30-\x3f]*[\x20-\x2f]*[\x40-\x7e]|'
    # OSC sequences
    r'\x1b\][^\x07\x1b]*(?:\x07|\x1b\\)|'
    # APC sequences
    r'\x1b_[^\x1b\x07]*(?:\x07|\x1b\\)|'
    # DCS sequences
    r'\x1bP[^\x1b\x07]*(?:\x07|\x1b\\)|'
    # PM sequences
    r'\x1b\^[^\x1b\x07]*(?:\x07|\x1b\\)|'
    # Character set designation
    r'\x1b[()].|'
    # Fe sequences
    r'\x1b[\x40-\x5f]|'
    # Fp sequences
    r'\x1b[78=>g]'
)

# Cursor right movement: CSI [n] C, parameter may be parsed by width()
CURSOR_RIGHT_SEQUENCE = re.compile(r'\x1b\[(\d*)C')

# Cursor left movement: CSI [n] D, parameter may be parsed by width()
CURSOR_LEFT_SEQUENCE = re.compile(r'\x1b\[(\d*)D')

# Indeterminate effect sequences - raise ValueError in 'strict' mode. The effects of these sequences
# are likely to be undesirable, moving the cursor vertically or to any unknown position, and
# otherwise not managed by the 'width' method of this library.
#
# This table was created initially with code generation by extraction of termcap library with
# techniques used at 'blessed' library runtime for 'xterm', 'alacritty', 'kitty', ghostty',
# 'screen', 'tmux', and others. Then, these common capabilities were merged into the list below.
INDETERMINATE_EFFECT_SEQUENCE = re.compile(
    '|'.join(f'(?:{_pattern})' for _pattern in (
        r'\x1b\[\d+;\d+r',           # change_scroll_region
        r'\x1b\[\d*K',               # erase_in_line (clr_eol, clr_bol)
        r'\x1b\[\d*J',               # erase_in_display (clr_eos, erase_display)
        r'\x1b\[\d*G',               # column_address
        r'\x1b\[\d+;\d+H',           # cursor_address
        r'\x1b\[\d*H',               # cursor_home
        r'\x1b\[\d*A',               # cursor_up
        r'\x1b\[\d*B',               # cursor_down
        r'\x1b\[\d*P',               # delete_character
        r'\x1b\[\d*M',               # delete_line
        r'\x1b\[\d*L',               # insert_line
        r'\x1b\[\d*@',               # insert_character
        r'\x1b\[\d+X',               # erase_chars
        r'\x1b\[\d*S',               # scroll_up (parm_index)
        r'\x1b\[\d*T',               # scroll_down (parm_rindex)
        r'\x1b\[\d*d',               # row_address
        r'\x1b\[\?1049[hl]',         # alternate screen buffer
        r'\x1b\[\?47[hl]',           # alternate screen (legacy)
        r'\x1b8',                    # restore_cursor
        r'\x1bD',                    # scroll_forward (index)
        r'\x1bM',                    # scroll_reverse (reverse index)
    ))
)
