"""
Normalization utilities for color values.

"""

# SPDX-License-Identifier: BSD-3-Clause

from ._definitions import _HEX_COLOR_RE
from ._types import IntegerRGB, IntTuple, PercentRGB, PercentTuple


def normalize_hex(hex_value: str) -> str:
    """
    Normalize a hexadecimal color value to a string consisting of the character `#`
    followed by six lowercase hexadecimal digits (what HTML5 terms a "valid lowercase
    simple color").

    If the supplied value cannot be interpreted as a hexadecimal color value,
    :exc:`ValueError` is raised. See :ref:`the conventions used by this module
    <conventions>` for information on acceptable formats for hexadecimal values.

    Examples:

    .. doctest::

        >>> normalize_hex("#0099cc")
        '#0099cc'
        >>> normalize_hex("#0099CC")
        '#0099cc'
        >>> normalize_hex("#09c")
        '#0099cc'
        >>> normalize_hex("#09C")
        '#0099cc'
        >>> normalize_hex("#0099gg")
        Traceback (most recent call last):
            ...
        ValueError: '#0099gg' is not a valid hexadecimal color value.
        >>> normalize_hex("0099cc")
        Traceback (most recent call last):
            ...
        ValueError: '0099cc' is not a valid hexadecimal color value.

    :param hex_value: The hexadecimal color value to normalize.
    :raises ValueError: when the input is not a valid hexadecimal color value.

    """
    if (match := _HEX_COLOR_RE.match(hex_value)) is None:
        raise ValueError(f'"{hex_value}" is not a valid hexadecimal color value.')
    hex_digits = match.group(1)
    if len(hex_digits) == 3:
        hex_digits = "".join(2 * s for s in hex_digits)
    return f"#{hex_digits.lower()}"


def _normalize_integer_rgb(value: int) -> int:
    """
    Internal normalization function for clipping integer values into the permitted
    range (0-255, inclusive).

    """
    return 0 if value < 0 else 255 if value > 255 else value


def normalize_integer_triplet(rgb_triplet: IntTuple) -> IntegerRGB:
    """
    Normalize an integer ``rgb()`` triplet so that all values are within the range
    0..255.

    Examples:

    .. doctest::

        >>> normalize_integer_triplet((128, 128, 128))
        IntegerRGB(red=128, green=128, blue=128)
        >>> normalize_integer_triplet((0, 0, 0))
        IntegerRGB(red=0, green=0, blue=0)
        >>> normalize_integer_triplet((255, 255, 255))
        IntegerRGB(red=255, green=255, blue=255)
        >>> normalize_integer_triplet((270, -20, -0))
        IntegerRGB(red=255, green=0, blue=0)

    :param rgb_triplet: The percentage `rgb()` triplet to normalize.

    """
    return IntegerRGB._make(_normalize_integer_rgb(value) for value in rgb_triplet)


def _normalize_percent_rgb(value: str) -> str:
    """
    Internal normalization function for clipping percent values into the permitted
    range (0%-100%, inclusive).

    """
    value = value.split("%")[0]
    percent = float(value) if "." in value else int(value)

    return "0%" if percent < 0 else "100%" if percent > 100 else f"{percent}%"


def normalize_percent_triplet(rgb_triplet: PercentTuple) -> PercentRGB:
    """
    Normalize a percentage ``rgb()`` triplet so that all values are within the range
    0%..100%.

    Examples:

    .. doctest::

       >>> normalize_percent_triplet(("50%", "50%", "50%"))
       PercentRGB(red='50%', green='50%', blue='50%')
       >>> normalize_percent_triplet(("0%", "100%", "0%"))
       PercentRGB(red='0%', green='100%', blue='0%')
       >>> normalize_percent_triplet(("-10%", "-0%", "500%"))
       PercentRGB(red='0%', green='0%', blue='100%')

    :param rgb_triplet: The percentage `rgb()` triplet to normalize.

    """
    return PercentRGB._make(_normalize_percent_rgb(value) for value in rgb_triplet)


def _percent_to_integer(percent: str) -> int:
    """
    Internal helper for converting a percentage value to an integer between 0 and
    255 inclusive.

    """
    return int(round(float(percent.split("%")[0]) / 100 * 255))
