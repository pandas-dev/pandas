"""
HTML5 color algorithms.

Note that these functions are written in a way that may seem strange to developers
familiar with Python, because they do not use the most efficient or idiomatic way of
accomplishing their tasks. This is because, for compliance, these functions are written
as literal translations into Python of the algorithms in HTML5:

https://html.spec.whatwg.org/multipage/common-microsyntaxes.html#colours

For ease of understanding, the relevant steps of the algorithm from the standard are
included as comments interspersed in the implementation.

"""

# SPDX-License-Identifier: BSD-3-Clause

import string

from ._definitions import _CSS3_NAMES_TO_HEX
from ._types import HTML5SimpleColor, IntTuple


def html5_parse_simple_color(value: str) -> HTML5SimpleColor:
    """
    Apply the HTML5 simple color parsing algorithm.

    Examples:

    .. doctest::

        >>> html5_parse_simple_color("#ffffff")
        HTML5SimpleColor(red=255, green=255, blue=255)
        >>> html5_parse_simple_color("#fff")
        Traceback (most recent call last):
            ...
        ValueError: An HTML5 simple color must be a string seven characters long.

    :param value: The color to parse.
    :type value: :class:`str`, which must consist of exactly the character ``"#"``
        followed by six hexadecimal digits.
    :raises ValueError: when the given value is not a Unicode string of length 7,
       consisting of exactly the character ``#`` followed by six hexadecimal digits.

    """
    # 1. Let input be the string being parsed.
    #
    # 2. If input is not exactly seven characters long, then return an error.
    if not isinstance(value, str) or len(value) != 7:
        raise ValueError(
            "An HTML5 simple color must be a Unicode string seven characters long."
        )

    # 3. If the first character in input is not a U+0023 NUMBER SIGN character (#), then
    #    return an error.
    if not value.startswith("#"):
        raise ValueError(
            "An HTML5 simple color must begin with the character '#' (U+0023)."
        )

    # 4. If the last six characters of input are not all ASCII hex digits, then return
    #    an error.
    if not all(c in string.hexdigits for c in value[1:]):
        raise ValueError(
            "An HTML5 simple color must contain exactly six ASCII hex digits."
        )

    # 5. Let result be a simple color.
    #
    # 6. Interpret the second and third characters as a hexadecimal number and let the
    #    result be the red component of result.
    #
    # 7. Interpret the fourth and fifth characters as a hexadecimal number and let the
    #    result be the green component of result.
    #
    # 8. Interpret the sixth and seventh characters as a hexadecimal number and let the
    #    result be the blue component of result.
    #
    # 9. Return result.
    return HTML5SimpleColor(
        int(value[1:3], 16), int(value[3:5], 16), int(value[5:7], 16)
    )


def html5_serialize_simple_color(simple_color: IntTuple) -> str:
    """
    Apply the HTML5 simple color serialization algorithm.

    Examples:

    .. doctest::

        >>> html5_serialize_simple_color((0, 0, 0))
        '#000000'
        >>> html5_serialize_simple_color((255, 255, 255))
        '#ffffff'

    :param simple_color: The color to serialize.

    """
    red, green, blue = simple_color

    # 1. Let result be a string consisting of a single "#" (U+0023) character.
    #
    # 2. Convert the red, green, and blue components in turn to two-digit hexadecimal
    #    numbers using lowercase ASCII hex digits, zero-padding if necessary, and append
    #    these numbers to result, in the order red, green, blue.
    #
    # 3. Return result, which will be a valid lowercase simple color.
    return f"#{red:02x}{green:02x}{blue:02x}"


def html5_parse_legacy_color(value: str) -> HTML5SimpleColor:
    """
    Apply the HTML5 legacy color parsing algorithm.

    Note that, since this algorithm is intended to handle many types of
    malformed color values present in real-world Web documents, it is
    *extremely* forgiving of input, but the results of parsing inputs
    with high levels of "junk" (i.e., text other than a color value)
    may be surprising.

    Examples:

    .. doctest::

        >>> html5_parse_legacy_color("black")
        HTML5SimpleColor(red=0, green=0, blue=0)
        >>> html5_parse_legacy_color("chucknorris")
        HTML5SimpleColor(red=192, green=0, blue=0)
        >>> html5_parse_legacy_color("Window")
        HTML5SimpleColor(red=0, green=13, blue=0)

    :param value: The color to parse.

    :raises ValueError: when the given value is not a Unicode string, when it is the
       empty string, or when it is precisely the string ``"transparent"``.

    """
    # 1. Let input be the string being parsed.
    if not isinstance(value, str):
        raise ValueError(
            "HTML5 legacy color parsing requires a Unicode string as input."
        )

    # 2. If input is the empty string, then return an error.
    if value == "":
        raise ValueError("HTML5 legacy color parsing forbids empty string as a value.")

    # 3. Strip leading and trailing ASCII whitespace from input.
    value = value.strip()

    # 4. If input is an ASCII case-insensitive match for the string "transparent", then
    #    return an error.
    if value.lower() == "transparent":
        raise ValueError('HTML5 legacy color parsing forbids "transparent" as a value.')

    # 5. If input is an ASCII case-insensitive match for one of the named colors, then
    #    return the simple color corresponding to that keyword.
    #
    #    Note: CSS2 System Colors are not recognized.
    if keyword_hex := _CSS3_NAMES_TO_HEX.get(value.lower()):
        return html5_parse_simple_color(keyword_hex)

    # 6. If input's code point length is four, and the first character in input is
    #    U+0023 (#), and the last three characters of input are all ASCII hex digits,
    #    then:
    if (
        len(value) == 4
        and value.startswith("#")
        and all(c in string.hexdigits for c in value[1:])
    ):
        # 1. Let result be a simple color.
        #
        # 2. Interpret the second character of input as a hexadecimal digit; let the red
        #    component of result be the resulting number multiplied by 17.
        #
        # 3. Interpret the third character of input as a hexadecimal digit; let the
        #    green component of result be the resulting number multiplied by 17.
        #
        # 4. Interpret the fourth character of input as a hexadecimal digit; let the
        #    blue component of result be the resulting number multiplied by 17.
        result = HTML5SimpleColor(
            int(value[1], 16) * 17, int(value[2], 16) * 17, int(value[3], 16) * 17
        )

        # 5. Return result.
        return result

    # 7. Replace any code points greater than U+FFFF in input (i.e., any characters that
    #    are not in the basic multilingual plane) with the two-character string "00".
    value = "".join("00" if ord(c) > 0xFFFF else c for c in value)

    # 8. If input's code point length is greater than 128, truncate input, leaving only
    #    the first 128 characters.
    if len(value) > 128:
        value = value[:128]

    # 9. If the first character in input is a U+0023 NUMBER SIGN character (#), remove
    #    it.
    if value.startswith("#"):
        value = value[1:]

    # 10. Replace any character in input that is not an ASCII hex digit with the
    # character U+0030 DIGIT ZERO (0).
    value = "".join(c if c in string.hexdigits else "0" for c in value)

    # 11. While input's code point length is zero or not a multiple of three, append a
    #     U+0030 DIGIT ZERO (0) character to input.
    while (len(value) == 0) or (len(value) % 3 != 0):
        value += "0"

    # 12. Split input into three strings of equal code point length, to obtain three
    #     components. Let length be the code point length that all of those components
    #     have (one third the code point length of input).
    length = int(len(value) / 3)
    red = value[:length]
    green = value[length : length * 2]
    blue = value[length * 2 :]

    # 13. If length is greater than 8, then remove the leading length-8 characters in
    #     each component, and let length be 8.
    if length > 8:
        red, green, blue = (red[length - 8 :], green[length - 8 :], blue[length - 8 :])
        length = 8

    # 14. While length is greater than two and the first character in each component is
    #     a U+0030 DIGIT ZERO (0) character, remove that character and reduce length by
    #     one.
    while (length > 2) and (red[0] == "0" and green[0] == "0" and blue[0] == "0"):
        red, green, blue = (red[1:], green[1:], blue[1:])
        length -= 1

    # 15. If length is still greater than two, truncate each component, leaving only the
    #     first two characters in each.
    if length > 2:
        red, green, blue = (red[:2], green[:2], blue[:2])

    # 16. Let result be a simple color.
    #
    # 17. Interpret the first component as a hexadecimal number; let the red component
    #     of result be the resulting number.
    #
    # 18. Interpret the second component as a hexadecimal number; let the green
    #     component of result be the resulting number.
    #
    # 19. Interpret the third component as a hexadecimal number; let the blue component
    #     of result be the resulting number.
    #
    # 20. Return result.
    return HTML5SimpleColor(int(red, 16), int(green, 16), int(blue, 16))
