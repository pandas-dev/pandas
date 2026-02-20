"""
Type aliases for the `strings` standard library.

See Also:
    - https://docs.python.org/3/library/string.html
"""

from typing import Final, Literal as L, TypeAlias  # noqa: N817

__all__ = (
    "DIGITS",
    "DIGITS_BIN",
    "DIGITS_HEX",
    "DIGITS_OCT",
    "LETTERS",
    "LETTERS_LOWER",
    "LETTERS_UPPER",
    "PRINTABLE",
    "PUNCTUATION",
    "WHITESPACE",
    "Digit",
    "DigitBin",
    "DigitHex",
    "DigitOct",
    "Letter",
    "LetterLower",
    "LetterUpper",
    "Printable",
    "Punctuation",
    "Whitespace",
)


DigitBin: TypeAlias = L["0", "1"]
DIGITS_BIN: Final = "0", "1"

# compatible with `string.octdigits`
DigitOct: TypeAlias = L["0", "1", "2", "3", "4", "5", "6", "7"]
DIGITS_OCT: Final = "0", "1", "2", "3", "4", "5", "6", "7"

# compatible with `string.hexdigits`
DigitHex: TypeAlias = L[
    "0", "1", "2", "3", "4", "5", "6", "7", "8", "9",
    "a", "b", "c", "d", "e", "f",
    "A", "B", "C", "D", "E", "F",
]  # fmt: skip
DIGITS_HEX: Final = (
    "0", "1", "2", "3", "4", "5", "6", "7", "8", "9",
    "a", "b", "c", "d", "e", "f", "A", "B", "C", "D", "E", "F",
)  # fmt: skip

# compatible with `string.digits`
Digit: TypeAlias = L["0", "1", "2", "3", "4", "5", "6", "7", "8", "9"]
DIGITS: Final = "0", "1", "2", "3", "4", "5", "6", "7", "8", "9"


# compatible with `string.ascii_letters`
LetterLower: TypeAlias = L[
    "a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m",
    "n", "o", "p", "q", "r", "s", "t", "u", "v", "w", "x", "y", "z",
]  # fmt: skip
LETTERS_LOWER: Final = (
    "a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m",
    "n", "o", "p", "q", "r", "s", "t", "u", "v", "w", "x", "y", "z",
)  # fmt: skip

# compatible with `string.ascii_lowercase`
LetterUpper: TypeAlias = L[
    "A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M",
    "N", "O", "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y", "Z",
]  # fmt: skip
LETTERS_UPPER: Final = (
    "A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M",
    "N", "O", "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y", "Z",
)  # fmt: skip

# compatible with `string.ascii_letters`
Letter: TypeAlias = LetterLower | LetterUpper
LETTERS: Final = LETTERS_LOWER + LETTERS_UPPER


# compatible with `string.punctuation`
Punctuation: TypeAlias = L[
    "!", '"', "#", "$", "%", "&", "'", "(",
    ")", "*", "+", ",", "-", ".", "/", ":",
    ";", "<", "=", ">", "?", "@", "[", "\\",
    "]", "^", "_", "`", "{", "|", "}", "~",
]  # fmt: skip
PUNCTUATION: Final = (
    "!", '"', "#", "$", "%", "&", "'", "(",
    ")", "*", "+", ",", "-", ".", "/", ":",
    ";", "<", "=", ">", "?", "@", "[", "\\",
    "]", "^", "_", "`", "{", "|", "}", "~",
)  # fmt: skip

# compatible with `string.whitespace`
Whitespace: TypeAlias = L[" ", "\t", "\n", "\r", "\v", "\f"]
WHITESPACE: Final = " ", "\t", "\n", "\r", "\v", "\f"

# compatible with `string.printable`
Printable: TypeAlias = Digit | Letter | Punctuation | Whitespace
PRINTABLE: Final = DIGITS + LETTERS + PUNCTUATION + WHITESPACE
