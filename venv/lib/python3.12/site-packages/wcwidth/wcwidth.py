"""
This is a python implementation of wcwidth() and wcswidth().

https://github.com/jquast/wcwidth

from Markus Kuhn's C code, retrieved from:

    http://www.cl.cam.ac.uk/~mgk25/ucs/wcwidth.c

This is an implementation of wcwidth() and wcswidth() (defined in
IEEE Std 1002.1-2001) for Unicode.

http://www.opengroup.org/onlinepubs/007904975/functions/wcwidth.html
http://www.opengroup.org/onlinepubs/007904975/functions/wcswidth.html

In fixed-width output devices, Latin characters all occupy a single
"cell" position of equal width, whereas ideographic CJK characters
occupy two such cells. Interoperability between terminal-line
applications and (teletype-style) character terminals using the
UTF-8 encoding requires agreement on which character should advance
the cursor by how many cell positions. No established formal
standards exist at present on which Unicode character shall occupy
how many cell positions on character terminals. These routines are
a first attempt of defining such behavior based on simple rules
applied to data provided by the Unicode Consortium.

For some graphical characters, the Unicode standard explicitly
defines a character-cell width via the definition of the East Asian
FullWidth (F), Wide (W), Half-width (H), and Narrow (Na) classes.
In all these cases, there is no ambiguity about which width a
terminal shall use. For characters in the East Asian Ambiguous (A)
class, the width choice depends purely on a preference of backward
compatibility with either historic CJK or Western practice.
Choosing single-width for these characters is easy to justify as
the appropriate long-term solution, as the CJK practice of
displaying these characters as double-width comes from historic
implementation simplicity (8-bit encoded characters were displayed
single-width and 16-bit ones double-width, even for Greek,
Cyrillic, etc.) and not any typographic considerations.

Much less clear is the choice of width for the Not East Asian
(Neutral) class. Existing practice does not dictate a width for any
of these characters. It would nevertheless make sense
typographically to allocate two character cells to characters such
as for instance EM SPACE or VOLUME INTEGRAL, which cannot be
represented adequately with a single-width glyph. The following
routines at present merely assign a single-cell width to all
neutral characters, in the interest of simplicity. This is not
entirely satisfactory and should be reconsidered before
establishing a formal standard in this area. At the moment, the
decision which Not East Asian (Neutral) characters should be
represented by double-width glyphs cannot yet be answered by
applying a simple rule from the Unicode database content. Setting
up a proper standard for the behavior of UTF-8 character terminals
will require a careful analysis not only of each Unicode character,
but also of each presentation form, something the author of these
routines has avoided to do so far.

http://www.unicode.org/unicode/reports/tr11/

Latest version: http://www.cl.cam.ac.uk/~mgk25/ucs/wcwidth.c
"""

from __future__ import annotations

# std imports
from functools import lru_cache

from typing import TYPE_CHECKING

# local
from .bisearch import bisearch as _bisearch
from .grapheme import iter_graphemes
from .table_mc import CATEGORY_MC
from .sgr_state import (_SGR_PATTERN,
                        _SGR_STATE_DEFAULT,
                        _sgr_state_update,
                        _sgr_state_is_active,
                        _sgr_state_to_sequence)
from .table_vs16 import VS16_NARROW_TO_WIDE
from .table_wide import WIDE_EASTASIAN
from .table_zero import ZERO_WIDTH
from .control_codes import ILLEGAL_CTRL, VERTICAL_CTRL, HORIZONTAL_CTRL, ZERO_WIDTH_CTRL
from .table_grapheme import ISC_CONSONANT, EXTENDED_PICTOGRAPHIC, GRAPHEME_REGIONAL_INDICATOR
from .table_ambiguous import AMBIGUOUS_EASTASIAN
from .escape_sequences import (ZERO_WIDTH_PATTERN,
                               CURSOR_LEFT_SEQUENCE,
                               CURSOR_RIGHT_SEQUENCE,
                               INDETERMINATE_EFFECT_SEQUENCE)
from .unicode_versions import list_versions

if TYPE_CHECKING:  # pragma: no cover
    # std imports
    from collections.abc import Iterator

    from typing import Literal

# Pre-compute table references for the latest (and only) Unicode version.
_LATEST_VERSION = list_versions()[-1]
_ZERO_WIDTH_TABLE = ZERO_WIDTH[_LATEST_VERSION]
_WIDE_EASTASIAN_TABLE = WIDE_EASTASIAN[_LATEST_VERSION]
_AMBIGUOUS_TABLE = AMBIGUOUS_EASTASIAN[next(iter(AMBIGUOUS_EASTASIAN))]
_CATEGORY_MC_TABLE = CATEGORY_MC[_LATEST_VERSION]
_REGIONAL_INDICATOR_SET = frozenset(
    range(GRAPHEME_REGIONAL_INDICATOR[0][0], GRAPHEME_REGIONAL_INDICATOR[0][1] + 1)
)
_EMOJI_ZWJ_SET = frozenset(
    cp for lo, hi in EXTENDED_PICTOGRAPHIC for cp in range(lo, hi + 1)
) | _REGIONAL_INDICATOR_SET
_FITZPATRICK_RANGE = (0x1F3FB, 0x1F3FF)
# Indic_Syllabic_Category=Virama codepoints, from IndicSyllabicCategory.txt.
# These are structurally tied to their scripts and not expected to change.
# https://www.unicode.org/Public/UCD/latest/ucd/IndicSyllabicCategory.txt
_ISC_VIRAMA_SET = frozenset((
    0x094D,   # DEVANAGARI SIGN VIRAMA
    0x09CD,   # BENGALI SIGN VIRAMA
    0x0A4D,   # GURMUKHI SIGN VIRAMA
    0x0ACD,   # GUJARATI SIGN VIRAMA
    0x0B4D,   # ORIYA SIGN VIRAMA
    0x0BCD,   # TAMIL SIGN VIRAMA
    0x0C4D,   # TELUGU SIGN VIRAMA
    0x0CCD,   # KANNADA SIGN VIRAMA
    0x0D4D,   # MALAYALAM SIGN VIRAMA
    0x0DCA,   # SINHALA SIGN AL-LAKUNA
    0x1B44,   # BALINESE ADEG ADEG
    0xA806,   # SYLOTI NAGRI SIGN HASANTA
    0xA8C4,   # SAURASHTRA SIGN VIRAMA
    0xA9C0,   # JAVANESE PANGKON
    0x11046,  # BRAHMI VIRAMA
    0x110B9,  # KAITHI SIGN VIRAMA
    0x111C0,  # SHARADA SIGN VIRAMA
    0x11235,  # KHOJKI SIGN VIRAMA
    0x1134D,  # GRANTHA SIGN VIRAMA
    0x11442,  # NEWA SIGN VIRAMA
    0x114C2,  # TIRHUTA SIGN VIRAMA
    0x115BF,  # SIDDHAM SIGN VIRAMA
    0x1163F,  # MODI SIGN VIRAMA
    0x116B6,  # TAKRI SIGN VIRAMA
    0x11839,  # DOGRA SIGN VIRAMA
    0x119E0,  # NANDINAGARI SIGN VIRAMA
    0x11C3F,  # BHAIKSUKI SIGN VIRAMA
))
_ISC_CONSONANT_TABLE = ISC_CONSONANT

# In 'parse' mode, strings longer than this are checked for cursor-movement
# controls (BS, TAB, CR, cursor sequences); when absent, mode downgrades to
# 'ignore' to skip character-by-character parsing. The detection scan cost is
# negligible for long strings but wasted on short ones like labels or headings.
_WIDTH_FAST_PATH_MIN_LEN = 20

# Translation table to strip C0/C1 control characters for fast 'ignore' mode.
_CONTROL_CHAR_TABLE = str.maketrans('', '', (
    ''.join(chr(c) for c in range(0x00, 0x20)) +   # C0: NUL through US (including tab)
    '\x7f' +                                       # DEL
    ''.join(chr(c) for c in range(0x80, 0xa0))     # C1: U+0080-U+009F
))

# Unlike wcwidth.__all__, wcwidth.wcwidth.__all__ is NOT for the purpose of defining a public API,
# or what we prefer to be imported with statement, "from wcwidth.wcwidth import *".  Explicitly
# re-export imports here for no other reason than to satisfy the type checkers (mypy). Yak shavings.
__all__ = (
    'ZERO_WIDTH',
    'WIDE_EASTASIAN',
    'AMBIGUOUS_EASTASIAN',
    'VS16_NARROW_TO_WIDE',
    'list_versions',
    'wcwidth',
    'wcswidth',
    'width',
    'iter_sequences',
    'ljust',
    'rjust',
    'center',
    'clip',
    'strip_sequences',
    '_wcmatch_version',
    '_wcversion_value',
)


# maxsize=1024: western scripts need ~64 unique codepoints per session, but
# CJK sessions may use ~2000 of ~3500 common hanzi/kanji. 1024 accommodates
# heavy CJK use. Performance floor at 32; bisearch is ~100ns per miss.

@lru_cache(maxsize=1024)
def wcwidth(wc: str, unicode_version: str = 'auto', ambiguous_width: int = 1) -> int:  # pylint: disable=unused-argument
    r"""
    Given one Unicode codepoint, return its printable length on a terminal.

    :param wc: A single Unicode character.
    :param unicode_version: Ignored. Retained for backwards compatibility.

        .. deprecated:: 0.3.0
           Only the latest Unicode version is now shipped.

    :param ambiguous_width: Width to use for East Asian Ambiguous (A)
        characters. Default is ``1`` (narrow). Set to ``2`` for CJK contexts
        where ambiguous characters display as double-width. See
        :ref:`ambiguous_width` for details.
    :returns: The width, in cells, necessary to display the character of
        Unicode string character, ``wc``.  Returns 0 if the ``wc`` argument has
        no printable effect on a terminal (such as NUL '\0'), -1 if ``wc`` is
        not printable, or has an indeterminate effect on the terminal, such as
        a control character.  Otherwise, the number of column positions the
        character occupies on a graphic terminal (1 or 2) is returned.

    See :ref:`Specification` for details of cell measurement.
    """
    ucs = ord(wc) if wc else 0

    # small optimization: early return of 1 for printable ASCII, this provides
    # approximately 40% performance improvement for mostly-ascii documents, with
    # less than 1% impact to others.
    if 32 <= ucs < 0x7f:
        return 1

    # C0/C1 control characters are -1 for compatibility with POSIX-like calls
    if ucs and ucs < 32 or 0x07F <= ucs < 0x0A0:
        return -1

    # Zero width
    if _bisearch(ucs, _ZERO_WIDTH_TABLE):
        return 0

    # Wide (F/W categories)
    if _bisearch(ucs, _WIDE_EASTASIAN_TABLE):
        return 2

    # Ambiguous width (A category) - only when ambiguous_width=2
    if ambiguous_width == 2 and _bisearch(ucs, _AMBIGUOUS_TABLE):
        return 2

    return 1


def wcswidth(
    pwcs: str,
    n: int | None = None,
    unicode_version: str = 'auto',
    ambiguous_width: int = 1,
) -> int:
    """
    Given a unicode string, return its printable length on a terminal.

    :param pwcs: Measure width of given unicode string.
    :param n: When ``n`` is None (default), return the length of the entire
        string, otherwise only the first ``n`` characters are measured.

        Better to use string slicing capability, ``wcswidth(pwcs[:n])``, instead,
        for performance.  This argument is a holdover from the POSIX function for
        matching signatures. Be careful that ``n`` is at grapheme boundaries.

    :param unicode_version: Ignored. Retained for backwards compatibility.

        .. deprecated:: 0.3.0
           Only the latest Unicode version is now shipped.

    :param ambiguous_width: Width to use for East Asian Ambiguous (A)
        characters. Default is ``1`` (narrow). Set to ``2`` for CJK contexts.
    :returns: The width, in cells, needed to display the first ``n`` characters
        of the unicode string ``pwcs``.  Returns ``-1`` for C0 and C1 control
        characters!

    See :ref:`Specification` for details of cell measurement.
    """
    # pylint: disable=unused-argument,too-many-locals,too-many-statements
    # pylint: disable=too-complex,too-many-branches
    # This function intentionally kept long without delegating functions to reduce function calls in
    # "hot path", the overhead per-character adds up.

    # Fast path: pure ASCII printable strings are always width == length
    if n is None and pwcs.isascii() and pwcs.isprintable():
        return len(pwcs)

    # Select wcwidth call pattern for best lru_cache performance:
    # - ambiguous_width=1 (default): single-arg calls share cache with direct wcwidth() calls
    # - ambiguous_width=2: full positional args needed (results differ, separate cache is correct)
    _wcwidth = wcwidth if ambiguous_width == 1 else lambda c: wcwidth(c, 'auto', ambiguous_width)

    end = len(pwcs) if n is None else n
    total_width = 0
    idx = 0
    last_measured_idx = -2  # Track index of last measured char for VS16
    last_measured_ucs = -1  # Codepoint of last measured char (for deferred emoji check)
    last_was_virama = False  # Virama conjunct formation state
    conjunct_pending = False  # Deferred +1 for bare conjuncts (no trailing Mc)
    while idx < end:
        char = pwcs[idx]
        ucs = ord(char)
        if ucs == 0x200D:
            if last_was_virama:
                # ZWJ after virama requests explicit half-form rendering but
                # does not change cell count — consume ZWJ only, let the next
                # consonant be handled by the virama conjunct rule.
                idx += 1
            elif idx + 1 < end:
                # Emoji ZWJ: skip next character unconditionally.
                idx += 2
                last_was_virama = False
            else:
                idx += 1
                last_was_virama = False
            continue
        if ucs == 0xFE0F and last_measured_idx >= 0:
            # VS16 following a measured character: add 1 if that character is
            # known to be converted from narrow to wide by VS16.
            total_width += _bisearch(ord(pwcs[last_measured_idx]),
                                     VS16_NARROW_TO_WIDE["9.0.0"])
            last_measured_idx = -2  # Prevent double application
            # VS16 preserves emoji context: last_measured_ucs stays as the base
            idx += 1
            continue
        # Regional Indicator & Fitzpatrick: both above BMP (U+1F1E6+)
        if ucs > 0xFFFF:
            if ucs in _REGIONAL_INDICATOR_SET:
                # Lazy RI pairing: count preceding consecutive RIs only when the last one is
                # received, because RI's are received so rarely its better than per-loop tracking of
                # 'last char was an RI'.
                ri_before = 0
                j = idx - 1
                while j >= 0 and ord(pwcs[j]) in _REGIONAL_INDICATOR_SET:
                    ri_before += 1
                    j -= 1
                if ri_before % 2 == 1:
                    # Second RI in pair: contributes 0 (pair = one 2-cell flag) using an even-or-odd
                    # check to determine, 'CAUS' would be two flags, but 'CAU' would be 1 flag
                    # and wide 'U'.
                    idx += 1
                    last_measured_ucs = ucs
                    continue
                # First or unpaired RI: measured normally (width 2 from table)
            # Fitzpatrick modifier: zero-width when following emoji base
            elif (_FITZPATRICK_RANGE[0] <= ucs <= _FITZPATRICK_RANGE[1]
                  and last_measured_ucs in _EMOJI_ZWJ_SET):
                idx += 1
                continue
        # Virama conjunct formation: consonant following virama contributes 0 width.
        # See https://www.unicode.org/reports/tr44/#Indic_Syllabic_Category
        if last_was_virama and _bisearch(ucs, _ISC_CONSONANT_TABLE):
            last_measured_idx = idx
            last_measured_ucs = ucs
            last_was_virama = False
            conjunct_pending = True
            idx += 1
            continue
        wcw = _wcwidth(char)
        if wcw < 0:
            # early return -1 on C0 and C1 control characters
            return wcw
        if wcw > 0:
            if conjunct_pending:
                total_width += 1
                conjunct_pending = False
            last_measured_idx = idx
            last_measured_ucs = ucs
            last_was_virama = False
        elif last_measured_idx >= 0 and _bisearch(ucs, _CATEGORY_MC_TABLE):
            # Spacing Combining Mark (Mc) following a base character adds 1
            wcw = 1
            last_measured_idx = -2
            last_was_virama = False
            conjunct_pending = False
        else:
            last_was_virama = ucs in _ISC_VIRAMA_SET
        total_width += wcw
        idx += 1
    if conjunct_pending:
        total_width += 1
    return total_width


# NOTE: _wcversion_value and _wcmatch_version are no longer used internally
# by wcwidth since version 0.5.0 (only the latest Unicode version is shipped).
#
# They are retained for API compatibility with external tools like ucs-detect
# that may use these private functions.


@lru_cache(maxsize=128)
def _wcversion_value(ver_string: str) -> tuple[int, ...]:  # pragma: no cover
    """
    Integer-mapped value of given dotted version string.

    .. deprecated:: 0.3.0

        This function is no longer used internally by wcwidth but is retained
        for API compatibility with external tools.

    :param ver_string: Unicode version string, of form ``n.n.n``.
    :returns: tuple of digit tuples, ``tuple(int, [...])``.
    """
    retval = tuple(map(int, (ver_string.split('.'))))
    return retval


@lru_cache(maxsize=8)
def _wcmatch_version(given_version: str) -> str:  # pylint: disable=unused-argument
    """
    Return the supported Unicode version level.

    .. deprecated:: 0.3.0
        This function now always returns the latest version.

        This function is no longer used internally by wcwidth but is retained
        for API compatibility with external tools.

    :param given_version: Ignored. Any value is accepted for compatibility.
    :returns: The latest unicode version string.
    """
    return _LATEST_VERSION


def iter_sequences(text: str) -> Iterator[tuple[str, bool]]:
    r"""
    Iterate through text, yielding segments with sequence identification.

    This generator yields tuples of ``(segment, is_sequence)`` for each part
    of the input text, where ``is_sequence`` is ``True`` if the segment is
    a recognized terminal escape sequence.

    :param text: String to iterate through.
    :returns: Iterator of (segment, is_sequence) tuples.

    .. versionadded:: 0.3.0

    Example::

        >>> list(iter_sequences('hello'))
        [('hello', False)]
        >>> list(iter_sequences('\x1b[31mred'))
        [('\x1b[31m', True), ('red', False)]
        >>> list(iter_sequences('\x1b[1m\x1b[31m'))
        [('\x1b[1m', True), ('\x1b[31m', True)]
    """
    idx = 0
    text_len = len(text)
    segment_start = 0

    while idx < text_len:
        char = text[idx]

        if char == '\x1b':
            # Yield any accumulated non-sequence text
            if idx > segment_start:
                yield (text[segment_start:idx], False)

            # Try to match an escape sequence
            match = ZERO_WIDTH_PATTERN.match(text, idx)
            if match:
                yield (match.group(), True)
                idx = match.end()
            else:
                # Lone ESC or unrecognized - yield as sequence anyway
                yield (char, True)
                idx += 1
            segment_start = idx
        else:
            idx += 1

    # Yield any remaining text
    if segment_start < text_len:
        yield (text[segment_start:], False)


def _width_ignored_codes(text: str, ambiguous_width: int = 1) -> int:
    """
    Fast path for width() with control_codes='ignore'.

    Strips escape sequences and control characters, then measures remaining text.
    """
    return wcswidth(
        strip_sequences(text).translate(_CONTROL_CHAR_TABLE),
        ambiguous_width=ambiguous_width
    )


def width(
    text: str,
    *,
    control_codes: Literal['parse', 'strict', 'ignore'] = 'parse',
    tabsize: int = 8,
    ambiguous_width: int = 1,
) -> int:
    r"""
    Return printable width of text containing many kinds of control codes and sequences.

    Unlike :func:`wcswidth`, this function handles most control characters and many popular terminal
    output sequences.  Never returns -1.

    :param text: String to measure.
    :param control_codes: How to handle control characters and sequences:

        - ``'parse'`` (default): Track horizontal cursor movement from BS ``\b``, CR ``\r``, TAB
          ``\t``, and cursor left and right movement sequences.  Vertical movement (LF, VT, FF) and
          indeterminate sequences are zero-width. Never raises.
        - ``'strict'``: Like parse, but raises :exc:`ValueError` on control characters with
          indeterminate results of the screen or cursor, like clear or vertical movement. Generally,
          these should be handled with a virtual terminal emulator (like 'pyte').
        - ``'ignore'``: All C0 and C1 control characters and escape sequences are measured as
          width 0. This is the fastest measurement for text already filtered or known not to contain
          any kinds of control codes or sequences. TAB ``\t`` is zero-width; for tab expansion,
          pre-process: ``text.replace('\t', ' ' * 8)``.

    :param tabsize: Tab stop width for ``'parse'`` and ``'strict'`` modes. Default is 8.
        Must be positive. Has no effect when ``control_codes='ignore'``.
    :param ambiguous_width: Width to use for East Asian Ambiguous (A)
        characters. Default is ``1`` (narrow). Set to ``2`` for CJK contexts.
    :returns: Maximum cursor position reached, "extent", accounting for cursor movement sequences
        present in ``text`` according to given parameters.  This represents the rightmost column the
        cursor reaches.  Always a non-negative integer.

    :raises ValueError: If ``control_codes='strict'`` and control characters with indeterminate
        effects, such as vertical movement or clear sequences are encountered, or on unexpected
        C0 or C1 control code. Also raised when ``control_codes`` is not one of the valid values.

    .. versionadded:: 0.3.0

    Examples::

        >>> width('hello')
        5
        >>> width('コンニチハ')
        10
        >>> width('\x1b[31mred\x1b[0m')
        3
        >>> width('\x1b[31mred\x1b[0m', control_codes='ignore')  # same result (ignored)
        3
        >>> width('123\b4')     # backspace overwrites previous cell (outputs '124')
        3
        >>> width('abc\t')      # tab caused cursor to move to column 8
        8
        >>> width('1\x1b[10C')  # '1' + cursor right 10, cursor ends on column 11
        11
        >>> width('1\x1b[10C', control_codes='ignore')   # faster but wrong in this case
        1
    """
    # pylint: disable=too-complex,too-many-branches,too-many-statements,too-many-locals
    # This could be broken into sub-functions (#1, #3, and 6 especially), but for reduced overhead
    # considering this function is a likely "hot path", they are inlined, breaking many of our
    # complexity rules.

    # Fast path for ASCII printable (no tabs, escapes, or control chars)
    if text.isascii() and text.isprintable():
        return len(text)

    # Fast parse: if no horizontal cursor movements are possible, switch to 'ignore' mode.
    # Only check for longer strings - the detection overhead hurts short string performance.
    if control_codes == 'parse' and len(text) > _WIDTH_FAST_PATH_MIN_LEN:
        # Check for cursor-affecting control characters
        if '\b' not in text and '\t' not in text and '\r' not in text:
            # Check for escape sequences - if none, or only non-cursor-movement sequences
            if '\x1b' not in text or (
                not CURSOR_RIGHT_SEQUENCE.search(text) and
                not CURSOR_LEFT_SEQUENCE.search(text)
            ):
                control_codes = 'ignore'

    # Fast path for ignore mode -- this is useful if you know the text is already "clean"
    if control_codes == 'ignore':
        return _width_ignored_codes(text, ambiguous_width)

    strict = control_codes == 'strict'
    # Track absolute positions: tab stops need modulo on absolute column, CR resets to 0.
    # Initialize max_extent to 0 so backward movement (CR, BS) won't yield negative width.
    current_col = 0
    max_extent = 0
    idx = 0
    last_measured_idx = -2  # Track index of last measured char for VS16; -2 can never match idx-1
    last_measured_ucs = -1  # Codepoint of last measured char (for deferred emoji check)
    last_was_virama = False  # Virama conjunct formation state
    conjunct_pending = False  # Deferred +1 for bare conjuncts (no trailing Mc)
    text_len = len(text)

    # Select wcwidth call pattern for best lru_cache performance:
    # - ambiguous_width=1 (default): single-arg calls share cache with direct wcwidth() calls
    # - ambiguous_width=2: full positional args needed (results differ, separate cache is correct)
    _wcwidth = wcwidth if ambiguous_width == 1 else lambda c: wcwidth(c, 'auto', ambiguous_width)

    while idx < text_len:
        char = text[idx]

        # 1. Handle ESC sequences
        if char == '\x1b':
            match = ZERO_WIDTH_PATTERN.match(text, idx)
            if match:
                seq = match.group()
                if strict and INDETERMINATE_EFFECT_SEQUENCE.match(seq):
                    raise ValueError(f"Indeterminate cursor sequence at position {idx}")
                # Apply cursor movement
                right = CURSOR_RIGHT_SEQUENCE.match(seq)
                if right:
                    current_col += int(right.group(1) or 1)
                else:
                    left = CURSOR_LEFT_SEQUENCE.match(seq)
                    if left:
                        current_col = max(0, current_col - int(left.group(1) or 1))
                idx = match.end()
            else:
                idx += 1
            max_extent = max(max_extent, current_col)
            continue

        # 2. Handle illegal and vertical control characters (zero width, error in strict)
        if char in ILLEGAL_CTRL:
            if strict:
                raise ValueError(f"Illegal control character {ord(char):#x} at position {idx}")
            idx += 1
            continue

        if char in VERTICAL_CTRL:
            if strict:
                raise ValueError(f"Vertical movement character {ord(char):#x} at position {idx}")
            idx += 1
            continue

        # 3. Handle horizontal movement characters
        if char in HORIZONTAL_CTRL:
            if char == '\x09' and tabsize > 0:  # Tab
                current_col += tabsize - (current_col % tabsize)
            elif char == '\x08':  # Backspace
                if current_col > 0:
                    current_col -= 1
            elif char == '\x0d':  # Carriage return
                current_col = 0
            max_extent = max(max_extent, current_col)
            idx += 1
            continue

        # 4. Handle ZWJ
        if char == '\u200D':
            if last_was_virama:
                # ZWJ after virama requests explicit half-form rendering but
                # does not change cell count — consume ZWJ only, let the next
                # consonant be handled by the virama conjunct rule.
                idx += 1
            elif idx + 1 < text_len:
                # Emoji ZWJ: skip next character unconditionally.
                idx += 2
                last_was_virama = False
            else:
                idx += 1
                last_was_virama = False
            continue

        # 5. Handle other zero-width characters (control chars)
        if char in ZERO_WIDTH_CTRL:
            idx += 1
            continue

        ucs = ord(char)

        # 6. Handle VS16: converts preceding narrow character to wide
        if ucs == 0xFE0F:
            if last_measured_idx == idx - 1:
                if _bisearch(ord(text[last_measured_idx]), VS16_NARROW_TO_WIDE["9.0.0"]):
                    current_col += 1
                    max_extent = max(max_extent, current_col)
            # VS16 preserves emoji context: last_measured_ucs stays as the base
            idx += 1
            continue

        # 6b. Regional Indicator & Fitzpatrick: both above BMP (U+1F1E6+)
        if ucs > 0xFFFF:
            if ucs in _REGIONAL_INDICATOR_SET:
                # Lazy RI pairing: count preceding consecutive RIs
                ri_before = 0
                j = idx - 1
                while j >= 0 and ord(text[j]) in _REGIONAL_INDICATOR_SET:
                    ri_before += 1
                    j -= 1
                if ri_before % 2 == 1:
                    last_measured_ucs = ucs
                    idx += 1
                    continue
            # 6c. Fitzpatrick modifier: zero-width when following emoji base
            elif (_FITZPATRICK_RANGE[0] <= ucs <= _FITZPATRICK_RANGE[1]
                  and last_measured_ucs in _EMOJI_ZWJ_SET):
                idx += 1
                continue

        # 7. Virama conjunct formation: consonant following virama contributes 0 width.
        # See https://www.unicode.org/reports/tr44/#Indic_Syllabic_Category
        if last_was_virama and _bisearch(ucs, _ISC_CONSONANT_TABLE):
            last_measured_idx = idx
            last_measured_ucs = ucs
            last_was_virama = False
            conjunct_pending = True
            idx += 1
            continue

        # 8. Normal characters: measure with wcwidth
        w = _wcwidth(char)
        if w > 0:
            if conjunct_pending:
                current_col += 1
                conjunct_pending = False
            current_col += w
            max_extent = max(max_extent, current_col)
            last_measured_idx = idx
            last_measured_ucs = ucs
            last_was_virama = False
        elif last_measured_idx >= 0 and _bisearch(ucs, _CATEGORY_MC_TABLE):
            # Spacing Combining Mark (Mc) following a base character adds 1
            current_col += 1
            max_extent = max(max_extent, current_col)
            last_measured_idx = -2
            last_was_virama = False
            conjunct_pending = False
        else:
            last_was_virama = ucs in _ISC_VIRAMA_SET
        idx += 1

    if conjunct_pending:
        current_col += 1
        max_extent = max(max_extent, current_col)
    return max_extent


def ljust(
    text: str,
    dest_width: int,
    fillchar: str = ' ',
    *,
    control_codes: Literal['parse', 'strict', 'ignore'] = 'parse',
    ambiguous_width: int = 1,
) -> str:
    r"""
    Return text left-justified in a string of given display width.

    :param text: String to justify, may contain terminal sequences.
    :param dest_width: Total display width of result in terminal cells.
    :param fillchar: Single character for padding (default space). Must have
        display width of 1 (not wide, not zero-width, not combining). Unicode
        characters like ``'·'`` are acceptable. The width is not validated.
    :param control_codes: How to handle control sequences when measuring.
        Passed to :func:`width` for measurement.
    :param ambiguous_width: Width to use for East Asian Ambiguous (A)
        characters. Default is ``1`` (narrow). Set to ``2`` for CJK contexts.
    :returns: Text padded on the right to reach ``dest_width``.

    .. versionadded:: 0.3.0

    Example::

        >>> wcwidth.ljust('hi', 5)
        'hi   '
        >>> wcwidth.ljust('\x1b[31mhi\x1b[0m', 5)
        '\x1b[31mhi\x1b[0m   '
        >>> wcwidth.ljust('\U0001F468\u200D\U0001F469\u200D\U0001F467', 6)
        '👨‍👩‍👧    '
    """
    if text.isascii() and text.isprintable():
        text_width = len(text)
    else:
        text_width = width(text, control_codes=control_codes, ambiguous_width=ambiguous_width)
    padding_cells = max(0, dest_width - text_width)
    return text + fillchar * padding_cells


def rjust(
    text: str,
    dest_width: int,
    fillchar: str = ' ',
    *,
    control_codes: Literal['parse', 'strict', 'ignore'] = 'parse',
    ambiguous_width: int = 1,
) -> str:
    r"""
    Return text right-justified in a string of given display width.

    :param text: String to justify, may contain terminal sequences.
    :param dest_width: Total display width of result in terminal cells.
    :param fillchar: Single character for padding (default space). Must have
        display width of 1 (not wide, not zero-width, not combining). Unicode
        characters like ``'·'`` are acceptable. The width is not validated.
    :param control_codes: How to handle control sequences when measuring.
        Passed to :func:`width` for measurement.
    :param ambiguous_width: Width to use for East Asian Ambiguous (A)
        characters. Default is ``1`` (narrow). Set to ``2`` for CJK contexts.
    :returns: Text padded on the left to reach ``dest_width``.

    .. versionadded:: 0.3.0

    Example::

        >>> wcwidth.rjust('hi', 5)
        '   hi'
        >>> wcwidth.rjust('\x1b[31mhi\x1b[0m', 5)
        '   \x1b[31mhi\x1b[0m'
        >>> wcwidth.rjust('\U0001F468\u200D\U0001F469\u200D\U0001F467', 6)
        '    👨‍👩‍👧'
    """
    if text.isascii() and text.isprintable():
        text_width = len(text)
    else:
        text_width = width(text, control_codes=control_codes, ambiguous_width=ambiguous_width)
    padding_cells = max(0, dest_width - text_width)
    return fillchar * padding_cells + text


def center(
    text: str,
    dest_width: int,
    fillchar: str = ' ',
    *,
    control_codes: Literal['parse', 'strict', 'ignore'] = 'parse',
    ambiguous_width: int = 1,
) -> str:
    r"""
    Return text centered in a string of given display width.

    :param text: String to center, may contain terminal sequences.
    :param dest_width: Total display width of result in terminal cells.
    :param fillchar: Single character for padding (default space). Must have
        display width of 1 (not wide, not zero-width, not combining). Unicode
        characters like ``'·'`` are acceptable. The width is not validated.
    :param control_codes: How to handle control sequences when measuring.
        Passed to :func:`width` for measurement.
    :param ambiguous_width: Width to use for East Asian Ambiguous (A)
        characters. Default is ``1`` (narrow). Set to ``2`` for CJK contexts.
    :returns: Text padded on both sides to reach ``dest_width``.

    For odd-width padding, the extra cell goes on the right (matching
    Python's :meth:`str.center` behavior).

    .. versionadded:: 0.3.0

    Example::

        >>> wcwidth.center('hi', 6)
        '  hi  '
        >>> wcwidth.center('\x1b[31mhi\x1b[0m', 6)
        '  \x1b[31mhi\x1b[0m  '
        >>> wcwidth.center('\U0001F468\u200D\U0001F469\u200D\U0001F467', 6)
        '  👨‍👩‍👧  '
    """
    if text.isascii() and text.isprintable():
        text_width = len(text)
    else:
        text_width = width(text, control_codes=control_codes, ambiguous_width=ambiguous_width)
    total_padding = max(0, dest_width - text_width)
    # matching https://jazcap53.github.io/pythons-eccentric-strcenter.html
    left_pad = total_padding // 2 + (total_padding & dest_width & 1)
    right_pad = total_padding - left_pad
    return fillchar * left_pad + text + fillchar * right_pad


def strip_sequences(text: str) -> str:
    r"""
    Return text with all terminal escape sequences removed.

    Unknown or incomplete ESC sequences are preserved.

    :param text: String that may contain terminal escape sequences.
    :returns: The input text with all escape sequences stripped.

    .. versionadded:: 0.3.0

    Example::

        >>> strip_sequences('\x1b[31mred\x1b[0m')
        'red'
        >>> strip_sequences('hello')
        'hello'
        >>> strip_sequences('\x1b[1m\x1b[31mbold red\x1b[0m text')
        'bold red text'
    """
    return ZERO_WIDTH_PATTERN.sub('', text)


def clip(
    text: str,
    start: int,
    end: int,
    *,
    fillchar: str = ' ',
    tabsize: int = 8,
    ambiguous_width: int = 1,
    propagate_sgr: bool = True,
) -> str:
    r"""
    Clip text to display columns ``(start, end)`` while preserving all terminal sequences.

    This function extracts a substring based on visible column positions rather than
    character indices. Terminal escape sequences are preserved in the output since
    they have zero display width. If a wide character (width 2) would be split at
    either boundary, it is replaced with ``fillchar``.

    TAB characters (``\t``) are expanded to spaces up to the next tab stop,
    controlled by the ``tabsize`` parameter.

    Other cursor movement characters (backspace, carriage return) and cursor
    movement sequences are passed through unchanged as zero-width.

    :param text: String to clip, may contain terminal escape sequences.
    :param start: Absolute starting column (inclusive, 0-indexed).
    :param end: Absolute ending column (exclusive).
    :param fillchar: Character to use when a wide character must be split at
        a boundary (default space). Must have display width of 1.
    :param tabsize: Tab stop width (default 8). Set to 0 to pass tabs through
        as zero-width (preserved in output but don't advance column position).
    :param ambiguous_width: Width to use for East Asian Ambiguous (A)
        characters. Default is ``1`` (narrow). Set to ``2`` for CJK contexts.
    :param propagate_sgr: If True (default), SGR (terminal styling) sequences
        are propagated. The result begins with any active style at the start
        position and ends with a reset sequence if styles are active.
    :returns: Substring of ``text`` spanning display columns ``(start, end)``,
        with all terminal sequences preserved and wide characters at boundaries
        replaced with ``fillchar``.

    SGR (terminal styling) sequences are propagated by default. The result
    begins with any active style and ends with a reset::

        >>> clip('\x1b[1;34mHello world\x1b[0m', 6, 11)
        '\x1b[1;34mworld\x1b[0m'

    Set ``propagate_sgr=False`` to disable this behavior.

    .. versionadded:: 0.3.0

    .. versionchanged:: 0.5.0
       Added ``propagate_sgr`` parameter (default True).

    Example::

        >>> clip('hello world', 0, 5)
        'hello'
        >>> clip('中文字', 0, 3)  # Wide char split at column 3
        '中 '
        >>> clip('a\tb', 0, 10)  # Tab expanded to spaces
        'a       b'
    """
    # pylint: disable=too-complex,too-many-locals,too-many-branches,too-many-statements,too-many-nested-blocks
    # Again, for 'hot path', we avoid additional delegate functions and accept the cost
    # of complexity for improved python performance.
    start = max(start, 0)
    if end <= start:
        return ''

    # Fast path: printable ASCII only (no tabs, escape sequences, or wide or zero-width chars)
    if text.isascii() and text.isprintable():
        return text[start:end]

    # Fast path: no escape sequences means no SGR tracking needed
    if propagate_sgr and '\x1b' not in text:
        propagate_sgr = False

    # SGR tracking state (only when propagate_sgr=True)
    sgr_at_clip_start = None  # state when first visible char emitted (None = not yet)
    if propagate_sgr:
        sgr = _SGR_STATE_DEFAULT  # current SGR state, updated by all sequences

    output: list[str] = []
    col = 0
    idx = 0

    while idx < len(text):
        char = text[idx]

        # Early exit: past visible region, SGR captured, no escape ahead
        if col >= end and sgr_at_clip_start is not None and char != '\x1b':
            break

        # Handle escape sequences
        if char == '\x1b' and (match := ZERO_WIDTH_PATTERN.match(text, idx)):
            seq = match.group()
            if propagate_sgr and _SGR_PATTERN.match(seq):
                # Update SGR state; will be applied as prefix when visible content starts
                sgr = _sgr_state_update(sgr, seq)
            else:
                # Non-SGR sequences always preserved
                output.append(seq)
            idx = match.end()
            continue

        # Handle bare ESC (not a valid sequence)
        if char == '\x1b':
            output.append(char)
            idx += 1
            continue

        # TAB expansion
        if char == '\t':
            if tabsize > 0:
                next_tab = col + (tabsize - (col % tabsize))
                while col < next_tab:
                    if start <= col < end:
                        output.append(' ')
                        if propagate_sgr and sgr_at_clip_start is None:
                            sgr_at_clip_start = sgr
                    col += 1
            else:
                output.append(char)
            idx += 1
            continue

        # Grapheme clustering for everything else
        grapheme = next(iter_graphemes(text, start=idx))
        w = width(grapheme, ambiguous_width=ambiguous_width)

        if w == 0:
            if start <= col < end:
                output.append(grapheme)
        elif col >= start and col + w <= end:
            # Fully visible
            output.append(grapheme)
            if propagate_sgr and sgr_at_clip_start is None:
                sgr_at_clip_start = sgr
            col += w
        elif col < end and col + w > start:
            # Partially visible (wide char at boundary)
            output.append(fillchar * (min(end, col + w) - max(start, col)))
            if propagate_sgr and sgr_at_clip_start is None:
                sgr_at_clip_start = sgr
            col += w
        else:
            col += w

        idx += len(grapheme)

    result = ''.join(output)

    # Apply SGR prefix/suffix
    if sgr_at_clip_start is not None:
        if prefix := _sgr_state_to_sequence(sgr_at_clip_start):
            result = prefix + result
        if _sgr_state_is_active(sgr_at_clip_start):
            result += '\x1b[0m'

    return result
