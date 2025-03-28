"""
Types and type aliases used to represent colors in various formats.

"""

# SPDX-License-Identifier: BSD-3-Clause

import typing


class IntegerRGB(typing.NamedTuple):
    """
    :class:`~typing.NamedTuple` representing an integer RGB triplet.

    Has three fields, each of type :class:`int` and in the range 0-255 inclusive:

    .. attribute:: red

       The red portion of the color value.

    .. attribute:: green

       The green portion of the color value.

    .. attribute:: blue

       The blue portion of the color value.

    """

    red: int
    green: int
    blue: int


class PercentRGB(typing.NamedTuple):
    """
    :class:`~typing.NamedTuple` representing a percentage RGB triplet.

    Has three fields, each of type :class:`str` and representing a percentage value in
    the range 0%-100% inclusive:

    .. attribute:: red

       The red portion of the color value.

    .. attribute:: green

       The green portion of the color value.

    .. attribute:: blue

       The blue portion of the color value.

    """

    red: str
    green: str
    blue: str


class HTML5SimpleColor(typing.NamedTuple):
    """
    :class:`~typing.NamedTuple` representing an HTML5 simple color.

    Has three fields, each of type :class:`int` and in the range 0-255 inclusive:

    .. attribute:: red

       The red portion of the color value.

    .. attribute:: green

       The green portion of the color value.

    .. attribute:: blue

       The blue portion of the color value.

    """

    red: int
    green: int
    blue: int


# Union type representing the possible types of an integer RGB tuple.
IntTuple = typing.Union[IntegerRGB, HTML5SimpleColor, typing.Tuple[int, int, int]]

# Union type representing the possible types of a percentage RGB tuple.
PercentTuple = typing.Union[PercentRGB, typing.Tuple[str, str, str]]
