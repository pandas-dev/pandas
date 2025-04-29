"""
Names for the standard and extended color set.
Extended set is similar to `vim wiki <http://vim.wikia.com/wiki/Xterm256_color_names_for_console_Vim>`_, `colored <https://pypi.python.org/pypi/colored>`_, etc. Colors based on `wikipedia <https://en.wikipedia.org/wiki/ANSI_escape_code#Colors>`_.

You can access the index of the colors with names.index(name). You can access the
rgb values with ``r=int(html[n][1:3],16)``, etc.
"""

from __future__ import annotations

color_names = """\
black
red
green
yellow
blue
magenta
cyan
light_gray
dark_gray
light_red
light_green
light_yellow
light_blue
light_magenta
light_cyan
white
grey_0
navy_blue
dark_blue
blue_3
blue_3a
blue_1
dark_green
deep_sky_blue_4
deep_sky_blue_4a
deep_sky_blue_4b
dodger_blue_3
dodger_blue_2
green_4
spring_green_4
turquoise_4
deep_sky_blue_3
deep_sky_blue_3a
dodger_blue_1
green_3
spring_green_3
dark_cyan
light_sea_green
deep_sky_blue_2
deep_sky_blue_1
green_3a
spring_green_3a
spring_green_2
cyan_3
dark_turquoise
turquoise_2
green_1
spring_green_2a
spring_green_1
medium_spring_green
cyan_2
cyan_1
dark_red
deep_pink_4
purple_4
purple_4a
purple_3
blue_violet
orange_4
grey_37
medium_purple_4
slate_blue_3
slate_blue_3a
royal_blue_1
chartreuse_4
dark_sea_green_4
pale_turquoise_4
steel_blue
steel_blue_3
cornflower_blue
chartreuse_3
dark_sea_green_4a
cadet_blue
cadet_blue_a
sky_blue_3
steel_blue_1
chartreuse_3a
pale_green_3
sea_green_3
aquamarine_3
medium_turquoise
steel_blue_1a
chartreuse_2a
sea_green_2
sea_green_1
sea_green_1a
aquamarine_1
dark_slate_gray_2
dark_red_a
deep_pink_4a
dark_magenta
dark_magenta_a
dark_violet
purple
orange_4a
light_pink_4
plum_4
medium_purple_3
medium_purple_3a
slate_blue_1
yellow_4
wheat_4
grey_53
light_slate_grey
medium_purple
light_slate_blue
yellow_4_a
dark_olive_green_3
dark_sea_green
light_sky_blue_3
light_sky_blue_3a
sky_blue_2
chartreuse_2
dark_olive_green_3a
pale_green_3a
dark_sea_green_3
dark_slate_gray_3
sky_blue_1
chartreuse_1
light_green_a
light_green_b
pale_green_1
aquamarine_1a
dark_slate_gray_1
red_3
deep_pink_4b
medium_violet_red
magenta_3
dark_violet_a
purple_a
dark_orange_3
indian_red
hot_pink_3
medium_orchid_3
medium_orchid
medium_purple_2
dark_goldenrod
light_salmon_3
rosy_brown
grey_63
medium_purple_2a
medium_purple_1
gold_3
dark_khaki
navajo_white_3
grey_69
light_steel_blue_3
light_steel_blue
yellow_3
dark_olive_green_3b
dark_sea_green_3a
dark_sea_green_2
light_cyan_3
light_sky_blue_1
green_yellow
dark_olive_green_2
pale_green_1a
dark_sea_green_2a
dark_sea_green_1
pale_turquoise_1
red_3a
deep_pink_3
deep_pink_3a
magenta_3a
magenta_3b
magenta_2
dark_orange_3a
indian_red_a
hot_pink_3a
hot_pink_2
orchid
medium_orchid_1
orange_3
light_salmon_3a
light_pink_3
pink_3
plum_3
violet
gold_3a
light_goldenrod_3
tan
misty_rose_3
thistle_3
plum_2
yellow_3a
khaki_3
light_goldenrod_2
light_yellow_3
grey_84
light_steel_blue_1
yellow_2
dark_olive_green_1
dark_olive_green_1a
dark_sea_green_1a
honeydew_2
light_cyan_1
red_1
deep_pink_2
deep_pink_1
deep_pink_1a
magenta_2a
magenta_1
orange_red_1
indian_red_1
indian_red_1a
hot_pink
hot_pink_a
medium_orchid_1a
dark_orange
salmon_1
light_coral
pale_violet_red_1
orchid_2
orchid_1
orange_1
sandy_brown
light_salmon_1
light_pink_1
pink_1
plum_1
gold_1
light_goldenrod_2a
light_goldenrod_2b
navajo_white_1
misty_rose_1
thistle_1
yellow_1
light_goldenrod_1
khaki_1
wheat_1
cornsilk_1
grey_10_0
grey_3
grey_7
grey_11
grey_15
grey_19
grey_23
grey_27
grey_30
grey_35
grey_39
grey_42
grey_46
grey_50
grey_54
grey_58
grey_62
grey_66
grey_70
grey_74
grey_78
grey_82
grey_85
grey_89
grey_93""".split()

EMPTY_SLICE = slice(None, None, None)

_greys = (
    3.4,
    7.4,
    11,
    15,
    19,
    23,
    26.7,
    30.49,
    34.6,
    38.6,
    42.4,
    46.4,
    50,
    54,
    58,
    62,
    66,
    69.8,
    73.8,
    77.7,
    81.6,
    85.3,
    89.3,
    93,
)
_grey_vals = [int(x / 100.0 * 16 * 16) for x in _greys]

_grey_html = ["#" + format(x, "02x") * 3 for x in _grey_vals]

_normals = [int(x, 16) for x in "0 5f 87 af d7 ff".split()]
_normal_html = [
    "#"
    + format(_normals[n // 36], "02x")
    + format(_normals[n // 6 % 6], "02x")
    + format(_normals[n % 6], "02x")
    for n in range(16 - 16, 232 - 16)
]

_base_pattern = [(n // 4, n // 2 % 2, n % 2) for n in range(8)]
_base_html = (
    [f"#{x[2] * 192:02x}{x[1] * 192:02x}{x[0] * 192:02x}" for x in _base_pattern]
    + ["#808080"]
    + [f"#{x[2] * 255:02x}{x[1] * 255:02x}{x[0] * 255:02x}" for x in _base_pattern][1:]
)
color_html = _base_html + _normal_html + _grey_html

color_codes_simple = list(range(8)) + list(range(60, 68))
"""Simple colors, remember that reset is #9, second half is non as common."""

# Attributes
attributes_ansi = {
    "bold": 1,
    "dim": 2,
    "italics": 3,
    "underline": 4,
    "reverse": 7,
    "hidden": 8,
    "strikeout": 9,
}

# Stylesheet
default_styles = {
    "warn": "fg red",
    "title": "fg cyan underline bold",
    "fatal": "fg red bold",
    "highlight": "bg yellow",
    "info": "fg blue",
    "success": "fg green",
}

# Functions to be used for color name operations


class FindNearest:
    """This is a class for finding the nearest color given rgb values.
    Different find methods are available."""

    def __init__(self, r: int, g: int, b: int) -> None:
        self.r = r
        self.b = b
        self.g = g

    def only_basic(self):
        """This will only return the first 8 colors!
        Breaks the colorspace into cubes, returns color"""
        midlevel = 0x40  # Since bright is not included

        # The colors are organised so that it is a
        # 3D cube, black at 0,0,0, white at 1,1,1
        # Compressed to linear_integers r,g,b
        # [[[0,1],[2,3]],[[4,5],[6,7]]]
        # r*1 + g*2 + b*4
        return (
            (self.r >= midlevel) * 1
            + (self.g >= midlevel) * 2
            + (self.b >= midlevel) * 4
        )

    def all_slow(self, color_slice: slice = EMPTY_SLICE) -> int:
        """This is a slow way to find the nearest color."""
        distances = [
            self._distance_to_color(color) for color in color_html[color_slice]
        ]
        return min(range(len(distances)), key=distances.__getitem__)

    def _distance_to_color(self, color: str) -> int:
        """This computes the distance to a color, should be minimized."""
        rgb = (int(color[1:3], 16), int(color[3:5], 16), int(color[5:7], 16))
        return (self.r - rgb[0]) ** 2 + (self.g - rgb[1]) ** 2 + (self.b - rgb[2]) ** 2

    def _distance_to_color_number(self, n: int) -> int:
        color = color_html[n]
        return self._distance_to_color(color)

    def only_colorblock(self) -> int:
        """This finds the nearest color based on block system, only works
        for 17-232 color values."""
        rint = min(
            range(len(_normals)), key=[abs(x - self.r) for x in _normals].__getitem__
        )
        bint = min(
            range(len(_normals)), key=[abs(x - self.b) for x in _normals].__getitem__
        )
        gint = min(
            range(len(_normals)), key=[abs(x - self.g) for x in _normals].__getitem__
        )
        return 16 + 36 * rint + 6 * gint + bint

    def only_simple(self) -> int:
        """Finds the simple color-block color."""
        return self.all_slow(slice(0, 16, None))

    def only_grey(self) -> int:
        """Finds the greyscale color."""
        rawval = (self.r + self.b + self.g) / 3
        n = min(
            range(len(_grey_vals)),
            key=[abs(x - rawval) for x in _grey_vals].__getitem__,
        )
        return n + 232

    def all_fast(self) -> int:
        """Runs roughly 8 times faster than the slow version."""
        colors = [self.only_simple(), self.only_colorblock(), self.only_grey()]
        distances = [self._distance_to_color_number(n) for n in colors]
        return colors[min(range(len(distances)), key=distances.__getitem__)]


def from_html(color: str) -> tuple[int, int, int]:
    """Convert html hex code to rgb."""
    if len(color) != 7 or color[0] != "#":
        raise ValueError("Invalid length of html code")
    return (int(color[1:3], 16), int(color[3:5], 16), int(color[5:7], 16))


def to_html(r, g, b):
    """Convert rgb to html hex code."""
    return f"#{r:02x}{g:02x}{b:02x}"
