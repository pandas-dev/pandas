"""
Translation of Python ``re`` patterns into RE2 patterns with the same meaning.

The ``.str`` methods take Python regular expressions, but pyarrow evaluates
regular expressions with RE2, which differs from ``re`` both in syntax and in the
meaning of shared syntax: ``\\w``, ``\\d``, ``\\s`` and ``\\b`` are ASCII-only,
``$`` does not match before a trailing newline, case folding differs, and so on.
A pattern is therefore never handed to pyarrow as written. It is parsed with
``re``'s own parser, and every node is written out in a form that means the same
in RE2 as it does in ``re``:

- character classes, categories and case-insensitive literals become explicit
  codepoint ranges, computed by running ``re`` itself over the codepoints;
- anchors become the RE2 assertion with the same meaning;
- quantifiers and groups are written out with RE2's syntax.

Constructs that have no such form (backreferences, lookaround, conditionals,
atomic groups, possessive quantifiers, repetition counts above RE2's limit) make
the translation fail and callers fall back to ``re``. Some constructs translate
only for data that cannot tell the two engines apart (e.g. a Unicode ``\\b`` on
ASCII-only strings); :class:`RE2Pattern` records these conditions, and the
callers check them against the data.
"""

from __future__ import annotations

from dataclasses import dataclass
import functools
import re
from re import _parser  # type: ignore[attr-defined]
import sys

import numpy as np

from pandas.compat import HAS_PYARROW

if HAS_PYARROW:
    import pyarrow as pa
    import pyarrow.compute as pc

_NUM_CODEPOINTS = 0x110000
_SURROGATES = slice(0xD800, 0xE000)
# RE2 rejects counted repetitions above this, including nested products
_RE2_MAX_REPEAT = 1000

_CATEGORY_SOURCES = {
    _parser.CATEGORY_DIGIT: r"\d",
    _parser.CATEGORY_NOT_DIGIT: r"\D",
    _parser.CATEGORY_SPACE: r"\s",
    _parser.CATEGORY_NOT_SPACE: r"\S",
    _parser.CATEGORY_WORD: r"\w",
    _parser.CATEGORY_NOT_WORD: r"\W",
}


@dataclass(frozen=True)
class RE2Pattern:
    """
    An RE2 pattern equivalent to a Python pattern, with the caveats that apply.

    Attributes
    ----------
    pattern : str
        The RE2 pattern. It carries no flags; all flags of the Python pattern are
        already applied.
    num_groups : int
        Number of capture groups; RE2 numbers them as ``re`` does.
    nullable : bool
        Whether the pattern can match the empty string. pyarrow's count, replace
        and split kernels do not step over empty matches as ``re`` does.
    has_assertions : bool
        Whether the pattern contains any zero-width assertion.
    start_assertions : bool
        Whether the pattern contains an assertion that looks at the text before
        the match (``^``, ``\\A``, ``\\b``, ``\\B``). pyarrow's count and split
        kernels restart the search on the remainder of the string, losing it.
    unicode_boundary : bool
        Whether the pattern contains a Unicode ``\\b`` or ``\\B``. RE2 has only
        the ASCII ones, which agree with ``re`` on ASCII-only strings.
    end_before_newline : bool
        Whether the pattern contains a ``$`` without MULTILINE. It is translated
        to ``\\z``, which agrees with ``re`` on strings that do not end in
        ``"\\n"``.
    non_boundary : bool
        Whether the pattern contains a ``\\B``. Before Python 3.14 it does not
        match in an empty string, while RE2's does.
    nullable_loop : bool
        Whether the pattern quantifies an expression that can match the empty
        string. ``re`` and RE2 handle the empty iterations differently, so they
        can find different matches and submatches, though for the same strings.
    optional_groups : bool
        Whether some group may not take part in a match. ``re`` reports such a
        group as None, pyarrow's extract kernel as an empty string.
    """

    pattern: str
    num_groups: int
    nullable: bool
    has_assertions: bool
    start_assertions: bool
    unicode_boundary: bool
    end_before_newline: bool
    non_boundary: bool
    nullable_loop: bool
    optional_groups: bool

    def matches_python_on(self, values: pa.ChunkedArray | pa.Array) -> bool:
        """
        Whether this pattern behaves as the Python pattern on ``values``.
        """
        # min_count=0 so that an empty or all-NA array reports True instead of NA
        if (
            self.unicode_boundary
            and not pc.all(pc.string_is_ascii(values), min_count=0).as_py()
        ):
            return False
        if (
            self.end_before_newline
            and pc.any(pc.ends_with(values, "\n"), min_count=0).as_py()
        ):
            return False
        if (
            self.non_boundary
            and sys.version_info < (3, 14)
            and pc.any(pc.equal(pc.binary_length(values), 0), min_count=0).as_py()
        ):
            return False
        return True


def translate(pattern: str, flags: int = 0, named_groups: bool = False):
    """
    Translate a Python pattern into an RE2 pattern with the same meaning.

    Parameters
    ----------
    pattern : str
        Python regular expression.
    flags : int
        ``re`` flags.
    named_groups : bool
        Write capture groups as named groups ``g1``, ``g2``, ... as pyarrow's
        extract kernel requires.

    Returns
    -------
    RE2Pattern or None
        None if the pattern cannot be expressed in RE2.

    Raises
    ------
    re.error
        If ``pattern`` is not a valid Python pattern.
    """
    # validate exactly as the Python implementation would, so that both raise
    #  the same error for an invalid pattern
    compiled = re.compile(pattern, flags)
    if not isinstance(pattern, str):
        return None
    return _translate(pattern, compiled.flags, compiled.groups, named_groups)


@functools.lru_cache(maxsize=256)
def _translate(
    pattern: str, flags: int, num_groups: int, named_groups: bool
) -> RE2Pattern | None:
    parsed = _parser.parse(pattern, flags)
    translator = _Translator(named_groups)
    text, nullable = translator.sequence(list(parsed), parsed.state.flags, 1)
    if translator.untranslatable or not _compiles_in_re2(text):
        # e.g. the expanded character classes exceed RE2's memory budget
        return None
    return RE2Pattern(
        pattern=text,
        num_groups=num_groups,
        nullable=nullable,
        has_assertions=translator.has_assertions,
        start_assertions=translator.start_assertions,
        unicode_boundary=translator.unicode_boundary,
        end_before_newline=translator.end_before_newline,
        non_boundary=translator.non_boundary,
        nullable_loop=translator.nullable_loop,
        optional_groups=translator.optional_groups,
    )


def _compiles_in_re2(pattern: str) -> bool:
    try:
        pc.match_substring_regex(pa.array([], type=pa.string()), pattern)
    except pa.ArrowInvalid:
        return False
    return True


def translate_template(compiled: re.Pattern[str], repl: str) -> tuple[str, bool] | None:
    """
    Translate a Python replacement template into an RE2 rewrite string.

    Returns
    -------
    tuple of (str, bool) or None
        The rewrite string and whether it refers to any group, or None if it
        cannot be expressed in RE2 (RE2 can refer to groups 0-9 only).

    Raises
    ------
    re.error
        If ``repl`` is not a valid template for ``compiled``.
    """
    if sys.version_info >= (3, 12):
        items = _parser.parse_template(repl, compiled)
    else:
        group_positions, items = _parser.parse_template(repl, compiled)
        for position, group in group_positions:
            items[position] = group

    parts = []
    uses_groups = False
    for item in items:
        if isinstance(item, int):
            if item > 9:
                return None
            parts.append(f"\\{item}")
            uses_groups = True
        elif item:
            parts.append(item.replace("\\", "\\\\"))
    return "".join(parts), uses_groups


class _Translator:
    def __init__(self, named_groups: bool) -> None:
        self.named_groups = named_groups
        # set when the pattern uses a construct RE2 has no equivalent for
        self.untranslatable = False
        self.has_assertions = False
        self.start_assertions = False
        self.unicode_boundary = False
        self.end_before_newline = False
        self.non_boundary = False
        self.nullable_loop = False
        self.optional_groups = False
        # number of enclosing alternatives / quantifiers with a minimum of zero
        self._optional_depth = 0

    def sequence(self, items, flags: int, repeat_product: int) -> tuple[str, bool]:
        parts = []
        nullable = True
        for op, av in items:
            text, node_nullable = self.node(op, av, flags, repeat_product)
            parts.append(text)
            nullable = nullable and node_nullable
        return "".join(parts), nullable

    def node(self, op, av, flags: int, repeat_product: int) -> tuple[str, bool]:
        if op is _parser.LITERAL:
            if flags & _parser.SRE_FLAG_IGNORECASE and av in _cased_codepoints_set():
                return self.char_class(((_parser.LITERAL, av),), False, flags), False
            return _re2_char(av), False
        if op is _parser.NOT_LITERAL:
            return self.char_class(((_parser.LITERAL, av),), True, flags), False
        if op is _parser.IN:
            negate = bool(av) and av[0][0] is _parser.NEGATE
            items = tuple(av[1:] if negate else av)
            return self.char_class(items, negate, flags), False
        if op is _parser.ANY:
            return ("(?s:.)" if flags & _parser.SRE_FLAG_DOTALL else r"[^\n]"), False
        if op is _parser.AT:
            return self.assertion(av, flags), True
        if op is _parser.SUBPATTERN:
            group, add_flags, del_flags, items = av
            inner_flags = (flags | add_flags) & ~del_flags
            text, nullable = self.sequence(items, inner_flags, repeat_product)
            if group is None:
                return f"(?:{text})", nullable
            if self._optional_depth:
                self.optional_groups = True
            if self.named_groups:
                return f"(?P<g{group}>{text})", nullable
            return f"({text})", nullable
        if op is _parser.BRANCH:
            _, alternatives = av
            if len(alternatives) > 1:
                self._optional_depth += 1
            translated = [
                self.sequence(items, flags, repeat_product) for items in alternatives
            ]
            if len(alternatives) > 1:
                self._optional_depth -= 1
            text = "|".join(text for text, _ in translated)
            return f"(?:{text})", any(nullable for _, nullable in translated)
        if op is _parser.MAX_REPEAT or op is _parser.MIN_REPEAT:
            return self.repeat(op, av, flags, repeat_product)
        # GROUPREF, GROUPREF_EXISTS, ASSERT, ASSERT_NOT, ATOMIC_GROUP,
        #  POSSESSIVE_REPEAT: RE2 has no equivalent
        self.untranslatable = True
        return "", False

    def char_class(self, items: tuple, negate: bool, flags: int) -> str:
        text = _char_class(items, negate, flags)
        if text is None:
            self.untranslatable = True
            return ""
        return text

    def repeat(self, op, av, flags: int, repeat_product: int) -> tuple[str, bool]:
        low, high, items = av
        unbounded = high is _parser.MAXREPEAT
        count = low if unbounded else high
        if count > _RE2_MAX_REPEAT or repeat_product * max(count, 1) > (
            _RE2_MAX_REPEAT
        ):
            self.untranslatable = True
            return "", False

        nullable_operand = _is_nullable(items)
        if nullable_operand:
            self.nullable_loop = True
        self._optional_depth += low == 0
        text, _ = self.sequence(items, flags, repeat_product * max(count, 1))
        self._optional_depth -= low == 0

        if unbounded:
            quantifier = {0: "*", 1: "+"}.get(low, f"{{{low},}}")
        elif low == high:
            quantifier = f"{{{low}}}"
        elif (low, high) == (0, 1):
            quantifier = "?"
        else:
            quantifier = f"{{{low},{high}}}"
        if op is _parser.MIN_REPEAT:
            quantifier += "?"
        return f"(?:{text}){quantifier}", low == 0 or nullable_operand

    def assertion(self, code, flags: int) -> str:
        self.has_assertions = True
        multiline = flags & _parser.SRE_FLAG_MULTILINE
        if code is _parser.AT_BEGINNING:
            self.start_assertions = True
            return "(?m:^)" if multiline else r"\A"
        if code is _parser.AT_BEGINNING_STRING:
            self.start_assertions = True
            return r"\A"
        if code is _parser.AT_END:
            if multiline:
                return "(?m:$)"
            # re's "$" also matches before a newline that ends the string; RE2
            #  cannot say this without a lookahead
            self.end_before_newline = True
            return r"\z"
        if code is _parser.AT_END_STRING:
            return r"\z"
        if code is _parser.AT_BOUNDARY or code is _parser.AT_NON_BOUNDARY:
            self.start_assertions = True
            if not flags & _parser.SRE_FLAG_ASCII:
                self.unicode_boundary = True
            if code is _parser.AT_BOUNDARY:
                return r"\b"
            self.non_boundary = True
            return r"\B"
        self.untranslatable = True
        return ""


def _is_nullable(items) -> bool:
    for op, av in items:
        if op is _parser.AT:
            continue
        if op is _parser.SUBPATTERN:
            if not _is_nullable(av[3]):
                return False
        elif op is _parser.BRANCH:
            if not any(_is_nullable(alternative) for alternative in av[1]):
                return False
        elif op is _parser.MAX_REPEAT or op is _parser.MIN_REPEAT:
            if av[0] != 0 and not _is_nullable(av[2]):
                return False
        else:
            return False
    return True


def _re2_char(codepoint: int) -> str:
    char = chr(codepoint)
    if char.isascii() and (char.isalnum() or char == "_"):
        return char
    return f"\\x{{{codepoint:X}}}"


def _python_char(codepoint: int) -> str:
    return f"\\U{codepoint:08x}"


@functools.cache
def _all_chars() -> str:
    """
    Every codepoint except the surrogates, which cannot occur in Arrow strings.
    """
    codepoints = np.arange(_NUM_CODEPOINTS, dtype="<u4")
    codepoints = codepoints[(codepoints < 0xD800) | (codepoints >= 0xE000)]
    return codepoints.tobytes().decode("utf-32-le")


def _codepoint_at(index: int) -> int:
    # inverse of the surrogate gap in _all_chars
    return index if index < 0xD800 else index + 0x800


@functools.cache
def _python_class_mask(source: str, flags: int) -> np.ndarray:
    """
    Boolean mask over all codepoints of the characters ``source`` matches in ``re``.
    """
    mask = np.zeros(_NUM_CODEPOINTS, dtype=bool)
    for match in re.finditer(f"(?:{source})+", _all_chars(), flags):
        start, stop = match.span()
        mask[_codepoint_at(start) : _codepoint_at(stop - 1) + 1] = True
    mask[_SURROGATES] = False
    mask.flags.writeable = False
    return mask


@functools.cache
def _cased_codepoints() -> np.ndarray:
    """
    Codepoints whose meaning in a pattern can change under IGNORECASE.

    Every other character matches, under IGNORECASE, exactly what it matches
    without it.
    """
    return np.array(
        [
            ord(char)
            for char in _all_chars()
            if char.lower() != char or char.upper() != char
        ],
        dtype=np.int64,
    )


@functools.cache
def _cased_codepoints_set() -> frozenset[int]:
    return frozenset(_cased_codepoints().tolist())


@functools.lru_cache(maxsize=1024)
def _char_class(items: tuple, negate: bool, flags: int) -> str | None:
    """
    RE2 character class matching what the ``re`` set ``[items]`` matches.

    Returns None if ``items`` contains something other than characters, ranges
    and the categories of ``_CATEGORY_SOURCES``.
    """
    ascii_flag = flags & _parser.SRE_FLAG_ASCII
    mask = np.zeros(_NUM_CODEPOINTS, dtype=bool)
    sources = []
    for op, av in items:
        if op is _parser.LITERAL:
            mask[av] = True
            sources.append(_python_char(av))
        elif op is _parser.RANGE:
            low, high = av
            mask[low : high + 1] = True
            sources.append(f"{_python_char(low)}-{_python_char(high)}")
        elif op is _parser.CATEGORY and av in _CATEGORY_SOURCES:
            source = _CATEGORY_SOURCES[av]
            mask |= _python_class_mask(source, ascii_flag)
            sources.append(source)
        else:
            return None

    if flags & _parser.SRE_FLAG_IGNORECASE:
        # ask re which of the characters with case mappings it matches
        compiled = re.compile(
            f"[{''.join(sources)}]", ascii_flag | _parser.SRE_FLAG_IGNORECASE
        )
        cased = _cased_codepoints()
        mask[cased] = [compiled.match(chr(cp)) is not None for cp in cased.tolist()]

    if negate:
        mask = ~mask
    mask[_SURROGATES] = False
    return _re2_class(mask)


def _re2_class(mask: np.ndarray) -> str:
    padded = np.concatenate(([False], mask, [False]))
    edges = np.flatnonzero(padded[1:] != padded[:-1])
    starts, stops = edges[::2], edges[1::2] - 1
    if len(starts) == 0:
        # RE2 has no empty class; this one matches no character
        return r"[^\x00-\x{10FFFF}]"
    if (
        len(starts) == 2
        and starts[0] == stops[0]
        and starts[1] == stops[1]
        and ord("A") <= starts[0] <= ord("Z")
        and starts[1] == starts[0] + 32
    ):
        # RE2 turns a class of an ASCII letter in both cases into a
        #  case-insensitive literal, which in some contexts it then folds with
        #  Unicode rules: [Kk] can match the Kelvin sign and [Ss] a long s
        return f"(?:{chr(starts[0])}|{chr(starts[1])})"
    parts = [
        _re2_char(start) if start == stop else f"{_re2_char(start)}-{_re2_char(stop)}"
        for start, stop in zip(starts.tolist(), stops.tolist(), strict=True)
    ]
    return f"[{''.join(parts)}]"
