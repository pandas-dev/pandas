"""
Definitions of valid formats and values for colors.

"""

# SPDX-License-Identifier: BSD-3-Clause

import re
from typing import List


def _reversedict(dict_to_reverse: dict) -> dict:
    """
    Internal helper for generating reverse mappings; given a dictionary, returns a
    new dictionary with keys and values swapped.

    """
    return {value: key for key, value in dict_to_reverse.items()}


_HEX_COLOR_RE = re.compile(r"^#([a-fA-F0-9]{3}|[a-fA-F0-9]{6})$")

HTML4 = "html4"
CSS2 = "css2"
CSS21 = "css21"
CSS3 = "css3"

_SUPPORTED_SPECIFICATIONS = (HTML4, CSS2, CSS21, CSS3)

_SPECIFICATION_ERROR_TEMPLATE = (
    f"{{spec}} is not a supported specification for color name lookups; "
    f"supported specifications are: {_SUPPORTED_SPECIFICATIONS}."
)

# Mappings of color names to normalized hexadecimal color values.
# --------------------------------------------------------------------------------

# The HTML 4 named colors.
#
# The canonical source for these color definitions is the HTML 4 specification:
#
# http://www.w3.org/TR/html401/types.html#h-6.5
#
# The file tests/definitions.py in the source distribution of this module downloads a
# copy of the HTML 4 standard and parses out the color names to ensure the values below
# are correct.
_HTML4_NAMES_TO_HEX = {
    "aqua": "#00ffff",
    "black": "#000000",
    "blue": "#0000ff",
    "fuchsia": "#ff00ff",
    "green": "#008000",
    "gray": "#808080",
    "lime": "#00ff00",
    "maroon": "#800000",
    "navy": "#000080",
    "olive": "#808000",
    "purple": "#800080",
    "red": "#ff0000",
    "silver": "#c0c0c0",
    "teal": "#008080",
    "white": "#ffffff",
    "yellow": "#ffff00",
}

# CSS2 used the same list as HTML 4.
_CSS2_NAMES_TO_HEX = _HTML4_NAMES_TO_HEX

# CSS2.1 added orange.
_CSS21_NAMES_TO_HEX = {"orange": "#ffa500", **_HTML4_NAMES_TO_HEX}

# The CSS3/SVG named colors.
#
# The canonical source for these color definitions is the SVG specification's color list
# (which was adopted as CSS 3's color definition):
#
# http://www.w3.org/TR/SVG11/types.html#ColorKeywords
#
# CSS3 also provides definitions of these colors:
#
# http://www.w3.org/TR/css3-color/#svg-color
#
# SVG provides the definitions as RGB triplets. CSS3 provides them both as RGB triplets
# and as hexadecimal. Since hex values are more common in real-world HTML and CSS, the
# mapping below is to hex values instead. The file tests/definitions.py in the source
# distribution of this module downloads a copy of the CSS3 color module and parses out
# the color names to ensure the values below are correct.
_CSS3_NAMES_TO_HEX = {
    "aliceblue": "#f0f8ff",
    "antiquewhite": "#faebd7",
    "aqua": "#00ffff",
    "aquamarine": "#7fffd4",
    "azure": "#f0ffff",
    "beige": "#f5f5dc",
    "bisque": "#ffe4c4",
    "black": "#000000",
    "blanchedalmond": "#ffebcd",
    "blue": "#0000ff",
    "blueviolet": "#8a2be2",
    "brown": "#a52a2a",
    "burlywood": "#deb887",
    "cadetblue": "#5f9ea0",
    "chartreuse": "#7fff00",
    "chocolate": "#d2691e",
    "coral": "#ff7f50",
    "cornflowerblue": "#6495ed",
    "cornsilk": "#fff8dc",
    "crimson": "#dc143c",
    "cyan": "#00ffff",
    "darkblue": "#00008b",
    "darkcyan": "#008b8b",
    "darkgoldenrod": "#b8860b",
    "darkgray": "#a9a9a9",
    "darkgrey": "#a9a9a9",
    "darkgreen": "#006400",
    "darkkhaki": "#bdb76b",
    "darkmagenta": "#8b008b",
    "darkolivegreen": "#556b2f",
    "darkorange": "#ff8c00",
    "darkorchid": "#9932cc",
    "darkred": "#8b0000",
    "darksalmon": "#e9967a",
    "darkseagreen": "#8fbc8f",
    "darkslateblue": "#483d8b",
    "darkslategray": "#2f4f4f",
    "darkslategrey": "#2f4f4f",
    "darkturquoise": "#00ced1",
    "darkviolet": "#9400d3",
    "deeppink": "#ff1493",
    "deepskyblue": "#00bfff",
    "dimgray": "#696969",
    "dimgrey": "#696969",
    "dodgerblue": "#1e90ff",
    "firebrick": "#b22222",
    "floralwhite": "#fffaf0",
    "forestgreen": "#228b22",
    "fuchsia": "#ff00ff",
    "gainsboro": "#dcdcdc",
    "ghostwhite": "#f8f8ff",
    "gold": "#ffd700",
    "goldenrod": "#daa520",
    "gray": "#808080",
    "grey": "#808080",
    "green": "#008000",
    "greenyellow": "#adff2f",
    "honeydew": "#f0fff0",
    "hotpink": "#ff69b4",
    "indianred": "#cd5c5c",
    "indigo": "#4b0082",
    "ivory": "#fffff0",
    "khaki": "#f0e68c",
    "lavender": "#e6e6fa",
    "lavenderblush": "#fff0f5",
    "lawngreen": "#7cfc00",
    "lemonchiffon": "#fffacd",
    "lightblue": "#add8e6",
    "lightcoral": "#f08080",
    "lightcyan": "#e0ffff",
    "lightgoldenrodyellow": "#fafad2",
    "lightgray": "#d3d3d3",
    "lightgrey": "#d3d3d3",
    "lightgreen": "#90ee90",
    "lightpink": "#ffb6c1",
    "lightsalmon": "#ffa07a",
    "lightseagreen": "#20b2aa",
    "lightskyblue": "#87cefa",
    "lightslategray": "#778899",
    "lightslategrey": "#778899",
    "lightsteelblue": "#b0c4de",
    "lightyellow": "#ffffe0",
    "lime": "#00ff00",
    "limegreen": "#32cd32",
    "linen": "#faf0e6",
    "magenta": "#ff00ff",
    "maroon": "#800000",
    "mediumaquamarine": "#66cdaa",
    "mediumblue": "#0000cd",
    "mediumorchid": "#ba55d3",
    "mediumpurple": "#9370db",
    "mediumseagreen": "#3cb371",
    "mediumslateblue": "#7b68ee",
    "mediumspringgreen": "#00fa9a",
    "mediumturquoise": "#48d1cc",
    "mediumvioletred": "#c71585",
    "midnightblue": "#191970",
    "mintcream": "#f5fffa",
    "mistyrose": "#ffe4e1",
    "moccasin": "#ffe4b5",
    "navajowhite": "#ffdead",
    "navy": "#000080",
    "oldlace": "#fdf5e6",
    "olive": "#808000",
    "olivedrab": "#6b8e23",
    "orange": "#ffa500",
    "orangered": "#ff4500",
    "orchid": "#da70d6",
    "palegoldenrod": "#eee8aa",
    "palegreen": "#98fb98",
    "paleturquoise": "#afeeee",
    "palevioletred": "#db7093",
    "papayawhip": "#ffefd5",
    "peachpuff": "#ffdab9",
    "peru": "#cd853f",
    "pink": "#ffc0cb",
    "plum": "#dda0dd",
    "powderblue": "#b0e0e6",
    "purple": "#800080",
    "red": "#ff0000",
    "rosybrown": "#bc8f8f",
    "royalblue": "#4169e1",
    "saddlebrown": "#8b4513",
    "salmon": "#fa8072",
    "sandybrown": "#f4a460",
    "seagreen": "#2e8b57",
    "seashell": "#fff5ee",
    "sienna": "#a0522d",
    "silver": "#c0c0c0",
    "skyblue": "#87ceeb",
    "slateblue": "#6a5acd",
    "slategray": "#708090",
    "slategrey": "#708090",
    "snow": "#fffafa",
    "springgreen": "#00ff7f",
    "steelblue": "#4682b4",
    "tan": "#d2b48c",
    "teal": "#008080",
    "thistle": "#d8bfd8",
    "tomato": "#ff6347",
    "turquoise": "#40e0d0",
    "violet": "#ee82ee",
    "wheat": "#f5deb3",
    "white": "#ffffff",
    "whitesmoke": "#f5f5f5",
    "yellow": "#ffff00",
    "yellowgreen": "#9acd32",
}


# Mappings of normalized hexadecimal color values to color names.
# --------------------------------------------------------------------------------

_HTML4_HEX_TO_NAMES = _reversedict(_HTML4_NAMES_TO_HEX)

_CSS2_HEX_TO_NAMES = _HTML4_HEX_TO_NAMES

_CSS21_HEX_TO_NAMES = _reversedict(_CSS21_NAMES_TO_HEX)

_CSS3_HEX_TO_NAMES = _reversedict(_CSS3_NAMES_TO_HEX)

# CSS3 defines both "gray" and "grey", as well as defining either spelling variant for
# other related colors like "darkgray"/"darkgrey", etc. For a "forward" lookup from
# name to hex, this is straightforward, but a "reverse" lookup from hex to name requires
# picking one spelling and being consistent about it.
#
# Since "gray" was the only spelling supported in HTML 4, CSS1, and CSS2, "gray" and its
# variants are chosen here.
_CSS3_HEX_TO_NAMES["#a9a9a9"] = "darkgray"
_CSS3_HEX_TO_NAMES["#2f4f4f"] = "darkslategray"
_CSS3_HEX_TO_NAMES["#696969"] = "dimgray"
_CSS3_HEX_TO_NAMES["#808080"] = "gray"
_CSS3_HEX_TO_NAMES["#d3d3d3"] = "lightgray"
_CSS3_HEX_TO_NAMES["#778899"] = "lightslategray"
_CSS3_HEX_TO_NAMES["#708090"] = "slategray"


_names_to_hex = {
    HTML4: _HTML4_NAMES_TO_HEX,
    CSS2: _CSS2_NAMES_TO_HEX,
    CSS21: _CSS21_NAMES_TO_HEX,
    CSS3: _CSS3_NAMES_TO_HEX,
}

_hex_to_names = {
    HTML4: _HTML4_HEX_TO_NAMES,
    CSS2: _CSS2_HEX_TO_NAMES,
    CSS21: _CSS21_HEX_TO_NAMES,
    CSS3: _CSS3_HEX_TO_NAMES,
}


def _get_name_to_hex_map(spec: str) -> dict:
    """
    Return the name-to-hex mapping for the given specification.

    :raises ValueError: when the given spec is not supported.

    """
    if spec not in _SUPPORTED_SPECIFICATIONS:
        raise ValueError(_SPECIFICATION_ERROR_TEMPLATE.format(spec=spec))
    return _names_to_hex[spec]


def _get_hex_to_name_map(spec: str) -> dict:
    """
    Return the hex-to-name mapping for the given specification.

    :raises ValueError: when the given spec is not supported.

    """
    if spec not in _SUPPORTED_SPECIFICATIONS:
        raise ValueError(_SPECIFICATION_ERROR_TEMPLATE.format(spec=spec))
    return _hex_to_names[spec]


def names(spec: str = CSS3) -> List[str]:
    """
    Return the list of valid color names for the given specification.

    The color names will be normalized to all-lowercase, and will be returned in
    alphabetical order.

    .. note:: **Spelling variants**

       Some values representing named gray colors can map to either of two names in
       CSS3, because it supports both ``"gray"`` and ``"grey"`` spelling variants for
       those colors. Functions which produce a name from a color value in other formats
       all normalize to the ``"gray"`` spelling for consistency with earlier CSS and
       HTML specifications which only supported ``"gray"``. Here, however, *all* valid
       names are returned, including -- for CSS3 -- both variant spellings for each of
       the affected ``"gray"``/``"grey"`` colors.

    Examples:

    .. doctest::

        >>> names(spec=HTML4)
        ['aqua', 'black', 'blue', 'fuchsia', 'gray', 'green',
         'lime', 'maroon', 'navy', 'olive', 'purple', 'red',
         'silver', 'teal', 'white', 'yellow']
        >>> names(spec="CSS1")
        Traceback (most recent call last):
            ...
        ValueError: "CSS1" is not a supported specification ...


    :raises ValueError: when the given spec is not supported.

    """
    if spec not in _SUPPORTED_SPECIFICATIONS:
        raise ValueError(_SPECIFICATION_ERROR_TEMPLATE.format(spec=spec))
    mapping = _names_to_hex[spec]
    return list(sorted(mapping.keys()))
