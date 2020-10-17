# being a bit too dynamic
import warnings

import matplotlib.cm as cm
import matplotlib.colors
import matplotlib.pyplot as plt
import numpy as np

from pandas.core.dtypes.common import is_list_like

import pandas.core.common as com


def get_standard_colors(
    num_colors: int,
    colormap=None,
    color_type: str = "default",
    color=None,
):
    colors = _get_colors(
        color=color,
        colormap=colormap,
        color_type=color_type,
        num_colors=num_colors,
    )

    return _cycle_colors(colors, num_colors=num_colors)


def _get_colors(*, color, colormap, color_type, num_colors):
    if color is None and colormap is not None:
        return _get_colors_from_colormap(colormap, num_colors=num_colors)
    elif color is not None:
        if colormap is not None:
            warnings.warn(
                "'color' and 'colormap' cannot be used simultaneously. Using 'color'"
            )
        return _get_colors_from_color(color)
    else:
        return _get_colors_from_color_type(color_type, num_colors=num_colors)


def _cycle_colors(colors, num_colors):
    # Append more colors by cycling if there is not enough color.
    # Extra colors will be ignored by matplotlib if there are more colors
    # than needed and nothing needs to be done here.
    if len(colors) < num_colors:
        try:
            multiple = num_colors // len(colors) - 1
        except ZeroDivisionError:
            raise ValueError("Invalid color argument: ''")
        mod = num_colors % len(colors)

        colors += multiple * colors
        colors += colors[:mod]

    return colors


def _get_colors_from_colormap(colormap, num_colors):
    colormap = _get_cmap_instance(colormap)
    return [colormap(num) for num in np.linspace(0, 1, num=num_colors)]


def _get_cmap_instance(colormap):
    if isinstance(colormap, str):
        cmap = colormap
        colormap = cm.get_cmap(colormap)
        if colormap is None:
            raise ValueError(f"Colormap {cmap} is not recognized")
    return colormap


def _get_colors_from_color(color):
    if is_list_like(color) and not isinstance(color, dict):
        return list(color)

    if _is_single_color(color):
        # GH #36972
        return [color]

    return color


def _get_colors_from_color_type(color_type, num_colors):
    if color_type == "default":
        return _get_default_colors(num_colors)
    elif color_type == "random":
        return _get_random_colors(num_colors)
    else:
        raise ValueError("color_type must be either 'default' or 'random'")


def _get_default_colors(num_colors):
    # need to call list() on the result to copy so we don't
    # modify the global rcParams below
    try:
        colors = [c["color"] for c in list(plt.rcParams["axes.prop_cycle"])]
    except KeyError:
        colors = list(plt.rcParams.get("axes.color_cycle", list("bgrcmyk")))

    if isinstance(colors, str):
        colors = list(colors)

    return colors[0:num_colors]


def _get_random_colors(num_colors):
    return [_random_color(num) for num in range(num_colors)]


def _random_color(column):
    """ Returns a random color represented as a list of length 3"""
    # GH17525 use common._random_state to avoid resetting the seed
    rs = com.random_state(column)
    return rs.rand(3).tolist()


def _is_single_color(color: str) -> bool:
    """Check if ``color`` is a single color.

    Examples of single colors:
        - 'r'
        - 'g'
        - 'red'
        - 'green'
        - 'C3'

    Parameters
    ----------
    color : string
        Color string.

    Returns
    -------
    bool
        True if ``color`` looks like a valid color.
        False otherwise.
    """
    conv = matplotlib.colors.ColorConverter()
    try:
        conv.to_rgba(color)
    except ValueError:
        return False
    else:
        return True
