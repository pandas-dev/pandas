"""
Functions for working with the color names and color value formats defined by the
HTML and CSS specifications for use in documents on the web.

See documentation (in docs/ directory of source distribution) for details of the
supported formats, conventions and conversions.

"""

# SPDX-License-Identifier: BSD-3-Clause

from ._conversion import (
    hex_to_name,
    hex_to_rgb,
    hex_to_rgb_percent,
    name_to_hex,
    name_to_rgb,
    name_to_rgb_percent,
    rgb_percent_to_hex,
    rgb_percent_to_name,
    rgb_percent_to_rgb,
    rgb_to_hex,
    rgb_to_name,
    rgb_to_rgb_percent,
)
from ._definitions import CSS2, CSS3, CSS21, HTML4
from ._html5 import (
    html5_parse_legacy_color,
    html5_parse_simple_color,
    html5_serialize_simple_color,
)
from ._normalization import (
    normalize_hex,
    normalize_integer_triplet,
    normalize_percent_triplet,
)
from ._types import HTML5SimpleColor, IntegerRGB, IntTuple, PercentRGB, PercentTuple

__version__ = "24.6.0"

__all__ = [
    "HTML4",
    "CSS2",
    "CSS21",
    "CSS3",
    "name_to_hex",
    "name_to_rgb",
    "name_to_rgb_percent",
    "hex_to_name",
    "hex_to_rgb",
    "hex_to_rgb_percent",
    "rgb_to_hex",
    "rgb_to_name",
    "rgb_to_rgb_percent",
    "rgb_percent_to_hex",
    "rgb_percent_to_name",
    "rgb_percent_to_rgb",
    "html5_parse_simple_color",
    "html5_parse_legacy_color",
    "html5_serialize_simple_color",
    "normalize_hex",
    "normalize_integer_triplet",
    "normalize_percent_triplet",
    "IntegerRGB",
    "PercentRGB",
    "HTML5SimpleColor",
    "IntTuple",
    "PercentTuple",
]
