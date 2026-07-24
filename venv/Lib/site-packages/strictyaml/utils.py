from strictyaml.ruamel.comments import CommentedSeq, CommentedMap
from strictyaml import exceptions
from re import compile
import decimal
import sys

if sys.version_info[:2] > (3, 4):
    from collections.abc import Iterable
else:
    from collections import Iterable

if sys.version_info[0] == 3:
    unicode = str


def flatten(items):
    """
    Yield items from any nested iterable.

    >>> list(flatten([[1, 2, 3], [[4, 5], 6, 7]]))
    [1, 2, 3, 4, 5, 6, 7]
    """
    for x in items:
        if isinstance(x, Iterable) and not isinstance(x, (str, bytes)):
            for sub_x in flatten(x):
                yield sub_x
        else:
            yield x


def has_number_type(value):
    """
    Is a value a number or a non-number?

    >>> has_number_type(3.5)
    True

    >>> has_number_type(3)
    True

    >>> has_number_type(decimal.Decimal("3.5"))
    True

    >>> has_number_type("3.5")
    False

    >>> has_number_type(True)
    False
    """
    return isinstance(value, (int, float, decimal.Decimal)) and not isinstance(
        value, bool
    )


def is_string(value):
    """
    Python 2/3 compatible way of checking if a value is a string.
    """
    return isinstance(value, unicode) or str(type(value)) in (
        "<type 'unicode'>",
        "<type 'str'>",
        "<class 'str'>",
    )


def is_integer(value):
    """
    Is a string a string of an integer?

    >>> is_integer("4")
    True

    >>> is_integer("4_000")
    True

    >>> is_integer("3.4")
    False
    """
    return compile(r"^[-+]?[0-9_]+$").match(value) is not None


def is_hexadecimal(value):
    """
    Is a string a string of a hexademcial integer?

    >>> is_hexadecimal("0xa1")
    True

    >>> is_hexadecimal("0XA1")
    True

    >>> is_hexadecimal("0xa1x")
    False

    >>> is_hexadecimal("xa1")
    False

    >>> is_hexadecimal("a1")
    False

    >>> is_hexadecimal("1")
    False
    """
    return compile(r"^0[xX]+[a-fA-F0-9]+$").match(value) is not None


def is_decimal(value):
    """
    Is a string a decimal?

    >>> is_decimal("4")
    True

    >>> is_decimal("4_000")
    True

    >>> is_decimal("3.5")
    True

    >>> is_decimal("4.")
    True

    >>> is_decimal("4.000_001")
    True

    >>> is_decimal("blah")
    False
    """
    return (
        compile(r"^[-+]?[0-9_]*(\.[0-9_]*)?([eE][-+]?[0-9_]+)?$").match(value)
        is not None
    )


def is_infinity(value):
    """
    Is string a valid representation for positive or negative infinity?

    Valid formats are:
    [+/-]inf, [+/-]INF, [+/-]Inf, [+/-].inf, [+/-].INF and [+/-].Inf

    >>> is_infinity(".inf")
    True

    >>> is_infinity("+.INF")
    True

    >>> is_infinity("-.Inf")
    True

    >>> is_infinity("Inf")
    True

    >>> is_infinity("INF")
    True

    >>> is_infinity("-INF")
    True

    >>> is_infinity("infinitesimal")
    False
    """
    return compile(r"^[-+]?\.?(?:inf|Inf|INF)$").match(value) is not None


def is_not_a_number(value):
    """
    Is string a valid representation for 'not a number'?

    Valid formats are: nan, NaN, NAN, .nan, .NaN, .NAN.

    >>> is_not_a_number(".nan")
    True

    >>> is_not_a_number(".NaN")
    True

    >>> is_not_a_number("NAN")
    True

    >>> is_not_a_number("nan")
    True

    >>> is_not_a_number("nanan")
    False

    >>> is_not_a_number("1e5")
    False
    """
    return compile(r"^\.?(?:nan|NaN|NAN)$").match(value) is not None


def comma_separated_positions(text):
    """
    Start and end positions of comma separated text items.

    Commas and trailing spaces should not be included.

    >>> comma_separated_positions("ABC, 2,3")
    [(0, 3), (5, 6), (7, 8)]
    """
    chunks = []
    start = 0
    end = 0
    for item in text.split(","):
        space_increment = 1 if item[0] == " " else 0
        start += space_increment  # Is there a space after the comma to ignore? ", "
        end += len(item.lstrip()) + space_increment
        chunks.append((start, end))
        start += len(item.lstrip()) + 1  # Plus comma
        end = start
    return chunks


def ruamel_structure(data, validator=None):
    """
    Take dicts and lists and return a strictyaml.ruamel style
    structure of CommentedMaps, CommentedSeqs and
    data.

    If a validator is presented and the type is unknown,
    it is checked against the validator to see if it will
    turn it back in to YAML.
    """
    if isinstance(data, dict):
        if len(data) == 0:
            raise exceptions.CannotBuildDocumentsFromEmptyDictOrList(
                "Document must be built with non-empty dicts and lists"
            )
        return CommentedMap(
            [
                (ruamel_structure(key), ruamel_structure(value))
                for key, value in data.items()
            ]
        )
    elif isinstance(data, list):
        if len(data) == 0:
            raise exceptions.CannotBuildDocumentsFromEmptyDictOrList(
                "Document must be built with non-empty dicts and lists"
            )
        return CommentedSeq([ruamel_structure(item) for item in data])
    elif isinstance(data, bool):
        return "yes" if data else "no"
    elif isinstance(data, (int, float)):
        return str(data)
    else:
        if not is_string(data):
            raise exceptions.CannotBuildDocumentFromInvalidData(
                (
                    "Document must be built from a combination of:\n"
                    "string, int, float, bool or nonempty list/dict\n\n"
                    "Instead, found variable with type '{}': '{}'"
                ).format(type(data).__name__, data)
            )
        return data
