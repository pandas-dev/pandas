"""
Functions which convert between various types of color values.

"""

# SPDX-License-Identifier: BSD-3-Clause

from ._definitions import CSS3, _get_hex_to_name_map, _get_name_to_hex_map
from ._normalization import (
    _percent_to_integer,
    normalize_hex,
    normalize_integer_triplet,
    normalize_percent_triplet,
)
from ._types import IntegerRGB, IntTuple, PercentRGB, PercentTuple

# Conversions from color names to other formats.
# --------------------------------------------------------------------------------


def name_to_hex(name: str, spec: str = CSS3) -> str:
    """
    Convert a color name to a normalized hexadecimal color value.

    The color name will be normalized to lower-case before being looked up.

    Examples:

    .. doctest::

        >>> name_to_hex("white")
        '#ffffff'
        >>> name_to_hex("navy")
        '#000080'
        >>> name_to_hex("goldenrod")
        '#daa520'
        >>> name_to_hex("goldenrod", spec=HTML4)
        Traceback (most recent call last):
            ...
        ValueError: "goldenrod" is not defined as a named color in html4.

    :param name: The color name to convert.
    :param spec: The specification from which to draw the list of color names. Default
       is :data:`CSS3`.
    :raises ValueError: when the given name has no definition in the given spec.

    """
    color_map = _get_name_to_hex_map(spec)
    if hex_value := color_map.get(name.lower()):
        return hex_value
    raise ValueError(f'"{name}" is not defined as a named color in {spec}')


def name_to_rgb(name: str, spec: str = CSS3) -> IntegerRGB:
    """
    Convert a color name to a 3-:class:`tuple` of :class:`int` suitable for use in
    an ``rgb()`` triplet specifying that color.

    The color name will be normalized to lower-case before being looked up.

    Examples:

    .. doctest::

        >>> name_to_rgb("white")
        IntegerRGB(red=255, green=255, blue=255)
        >>> name_to_rgb("navy")
        IntegerRGB(red=0, green=0, blue=128)
        >>> name_to_rgb("goldenrod")
        IntegerRGB(red=218, green=165, blue=32)

    :param name: The color name to convert.
    :param spec: The specification from which to draw the list of color names. Default
       is :data:`CSS3.`
    :raises ValueError: when the given name has no definition in the given spec.

    """
    return hex_to_rgb(name_to_hex(name, spec=spec))


def name_to_rgb_percent(name: str, spec: str = CSS3) -> PercentRGB:
    """
    Convert a color name to a 3-:class:`tuple` of percentages suitable for use in an
    ``rgb()`` triplet specifying that color.

    The color name will be normalized to lower-case before being looked up.

    Examples:

    .. doctest::

        >>> name_to_rgb_percent("white")
        PercentRGB(red='100%', green='100%', blue='100%')
        >>> name_to_rgb_percent("navy")
        PercentRGB(red='0%', green='0%', blue='50%')
        >>> name_to_rgb_percent("goldenrod")
        PercentRGB(red='85.49%', green='64.71%', blue='12.5%')

    :param name: The color name to convert.
    :param spec: The specification from which to draw the list of color names. Default
       is :data:`CSS3`.
    :raises ValueError: when the given name has no definition in the given spec.

    """
    return rgb_to_rgb_percent(name_to_rgb(name, spec=spec))


# Conversions from hexadecimal color values to other formats.
# --------------------------------------------------------------------------------


def hex_to_name(hex_value: str, spec: str = CSS3) -> str:
    """
    Convert a hexadecimal color value to its corresponding normalized color name, if
    any such name exists.

    The hexadecimal value will be normalized before being looked up.

    .. note:: **Spelling variants**

       Some values representing named gray colors can map to either of two names in
       CSS3, because it supports both ``"gray"`` and ``"grey"`` spelling variants for
       those colors. This function will always return the variant spelled ``"gray"``
       (such as ``"lightgray"`` instead of ``"lightgrey"``). See :ref:`the documentation
       on name conventions <color-name-conventions>` for details.

    Examples:

    .. doctest::

        >>> hex_to_name("#ffffff")
        'white'
        >>> hex_to_name("#fff")
        'white'
        >>> hex_to_name("#000080")
        'navy'
        >>> hex_to_name("#daa520")
        'goldenrod'
        >>> hex_to_name("#daa520", spec=HTML4)
        Traceback (most recent call last):
            ...
        ValueError: "#daa520" has no defined color name in html4.

    :param hex_value: The hexadecimal color value to convert.
    :param spec: The specification from which to draw the list of color names. Default
       is :data:`CSS3`.
    :raises ValueError: when the given color has no name in the given spec, or when the
       supplied hex value is invalid.

    """
    color_map = _get_hex_to_name_map(spec)
    if name := color_map.get(normalize_hex(hex_value)):
        return name
    raise ValueError(f'"{hex_value}" has no defined color name in {spec}.')


def hex_to_rgb(hex_value: str) -> IntegerRGB:
    """
    Convert a hexadecimal color value to a 3-:class:`tuple` of :class:`int` suitable
    for use in an ``rgb()`` triplet specifying that color.

    The hexadecimal value will be normalized before being converted.

    Examples:

    .. doctest::

        >>> hex_to_rgb("#fff")
        IntegerRGB(red=255, green=255, blue=255)
        >>> hex_to_rgb("#000080")
        IntegerRGB(red=0, green=0, blue=128)

    :param hex_value: The hexadecimal color value to convert.
    :raises ValueError: when the supplied hex value is invalid.

    """
    int_value = int(normalize_hex(hex_value)[1:], 16)
    return IntegerRGB(int_value >> 16, int_value >> 8 & 0xFF, int_value & 0xFF)


def hex_to_rgb_percent(hex_value: str) -> PercentRGB:
    """
    Convert a hexadecimal color value to a 3-:class:`tuple` of percentages suitable
    for use in an ``rgb()`` triplet representing that color.

    The hexadecimal value will be normalized before being converted.

    Examples:

    .. doctest::

        >>> hex_to_rgb_percent("#ffffff")
        PercentRGB(red='100%', green='100%', blue='100%')
        >>> hex_to_rgb_percent("#000080")
        PercentRGB(red='0%', green='0%', blue='50%')

    :param hex_value: The hexadecimal color value to convert.
    :raises ValueError: when the supplied hex value is invalid.

    """
    return rgb_to_rgb_percent(hex_to_rgb(hex_value))


# Conversions from  integer rgb() triplets to other formats.
# --------------------------------------------------------------------------------


def rgb_to_name(rgb_triplet: IntTuple, spec: str = CSS3) -> str:
    """
    Convert a 3-:class:`tuple` of :class:`int`, suitable for use in an ``rgb()``
    color triplet, to its corresponding normalized color name, if any such name exists.

    To determine the name, the triplet will be converted to a normalized hexadecimal
    value.

    .. note:: **Spelling variants**

       Some values representing named gray colors can map to either of two names in
       CSS3, because it supports both ``"gray"`` and ``"grey"`` spelling variants for
       those colors. This function will always return the variant spelled ``"gray"``
       (such as ``"lightgray"`` instead of ``"lightgrey"``). See :ref:`the documentation
       on name conventions <color-name-conventions>` for details.

    Examples:

    .. doctest::

        >>> rgb_to_name((255, 255, 255))
        'white'
        >>> rgb_to_name((0, 0, 128))
        'navy'

    :param rgb_triplet: The ``rgb()`` triplet.
    :param spec: The specification from which to draw the list of color names. Default
       is :data:`CSS3`.
    :raises ValueError: when the given color has no name in the given spec.

    """
    return hex_to_name(rgb_to_hex(normalize_integer_triplet(rgb_triplet)), spec=spec)


def rgb_to_hex(rgb_triplet: IntTuple) -> str:
    """
    Convert a 3-:class:`tuple` of :class:`int`, suitable for use in an ``rgb()``
    color triplet, to a normalized hexadecimal value for that color.

    Examples:

    .. doctest::

        >>> rgb_to_hex((255, 255, 255))
        '#ffffff'
        >>> rgb_to_hex((0, 0, 128))
        '#000080'

    :param rgb_triplet: The ``rgb()`` triplet.

    """
    red, green, blue = normalize_integer_triplet(rgb_triplet)
    return f"#{red:02x}{green:02x}{blue:02x}"


def rgb_to_rgb_percent(rgb_triplet: IntTuple) -> PercentRGB:
    """
    Convert a 3-:class:`tuple` of :class:`int`, suitable for use in an ``rgb()``
    color triplet, to a 3-:class:`tuple` of percentages suitable for use in representing
    that color.

    .. note:: **Floating-point precision**

       This function makes some trade-offs in terms of the accuracy of the final
       representation. For some common integer values, special-case logic is used to
       ensure a precise result (e.g., integer 128 will always convert to ``"50%"``,
       integer 32 will always convert to ``"12.5%"``), but for all other values a
       standard Python :class:`float` is used and rounded to two decimal places, which
       may result in a loss of precision for some values due to the inherent imprecision
       of `IEEE floating-point numbers <https://en.wikipedia.org/wiki/IEEE_754>`_.

    Examples:

    .. doctest::

        >>> rgb_to_rgb_percent((255, 255, 255))
        PercentRGB(red='100%', green='100%', blue='100%')
        >>> rgb_to_rgb_percent((0, 0, 128))
        PercentRGB(red='0%', green='0%', blue='50%')
        >>> rgb_to_rgb_percent((218, 165, 32))
        PercentRGB(red='85.49%', green='64.71%', blue='12.5%')

    :param rgb_triplet: The ``rgb()`` triplet.

    """
    # In order to maintain precision for common values,
    # special-case them.
    specials = {
        255: "100%",
        128: "50%",
        64: "25%",
        32: "12.5%",
        16: "6.25%",
        0: "0%",
    }
    return PercentRGB._make(
        specials.get(d, f"{d / 255.0 * 100:.02f}%")
        for d in normalize_integer_triplet(rgb_triplet)
    )


# Conversions from percentage rgb() triplets to other formats.
# --------------------------------------------------------------------------------


def rgb_percent_to_name(rgb_percent_triplet: PercentTuple, spec: str = CSS3) -> str:
    """
    Convert a 3-:class:`tuple` of percentages, suitable for use in an ``rgb()``
    color triplet, to its corresponding normalized color name, if any such name exists.

    To determine the name, the triplet will be converted to a normalized hexadecimal
    value.

    .. note:: **Spelling variants**

       Some values representing named gray colors can map to either of two names in
       CSS3, because it supports both ``"gray"`` and ``"grey"`` spelling variants for
       those colors. This function will always return the variant spelled ``"gray"``
       (such as ``"lightgray"`` instead of ``"lightgrey"``). See :ref:`the documentation
       on name conventions <color-name-conventions>` for details.

    Examples:

    .. doctest::

        >>> rgb_percent_to_name(("100%", "100%", "100%"))
        'white'
        >>> rgb_percent_to_name(("0%", "0%", "50%"))
        'navy'
        >>> rgb_percent_to_name(("85.49%", "64.71%", "12.5%"))
        'goldenrod'

    :param rgb_percent_triplet: The ``rgb()`` triplet.
    :param spec: The specification from which to draw the list of color names. Default
        is :data:`CSS3`.
    :raises ValueError: when the given color has no name in the given spec.

    """
    return rgb_to_name(
        rgb_percent_to_rgb(normalize_percent_triplet(rgb_percent_triplet)),
        spec=spec,
    )


def rgb_percent_to_hex(rgb_percent_triplet: PercentTuple) -> str:
    """
    Convert a 3-:class:`tuple` of percentages, suitable for use in an ``rgb()``
    color triplet, to a normalized hexadecimal color value for that color.

    Examples:

    .. doctest::

        >>> rgb_percent_to_hex(("100%", "100%", "0%"))
        '#ffff00'
        >>> rgb_percent_to_hex(("0%", "0%", "50%"))
        '#000080'
        >>> rgb_percent_to_hex(("85.49%", "64.71%", "12.5%"))
        '#daa520'

    :param rgb_percent_triplet: The ``rgb()`` triplet.

    """
    return rgb_to_hex(
        rgb_percent_to_rgb(normalize_percent_triplet(rgb_percent_triplet))
    )


def rgb_percent_to_rgb(rgb_percent_triplet: PercentTuple) -> IntegerRGB:
    """
    Convert a 3-:class:`tuple` of percentages, suitable for use in an ``rgb()``
    color triplet, to a 3-:class:`tuple` of :class:`int` suitable for use in
    representing that color.

    Some precision may be lost in this conversion. See the note regarding precision for
    :func:`~webcolors.rgb_to_rgb_percent` for details.

    Examples:

    .. doctest::

        >>> rgb_percent_to_rgb(("100%", "100%", "100%"))
        IntegerRGB(red=255, green=255, blue=255)
        >>> rgb_percent_to_rgb(("0%", "0%", "50%"))
        IntegerRGB(red=0, green=0, blue=128)
        >>> rgb_percent_to_rgb(("85.49%", "64.71%", "12.5%"))
        IntegerRGB(red=218, green=165, blue=32)

    :param rgb_percent_triplet: The ``rgb()`` triplet.

    """
    return IntegerRGB._make(
        map(
            _percent_to_integer,  # pylint: disable=protected-access
            normalize_percent_triplet(rgb_percent_triplet),
        )
    )
