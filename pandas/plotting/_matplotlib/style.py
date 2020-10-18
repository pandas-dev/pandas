from typing import (
    TYPE_CHECKING,
    Collection,
    List,
    Optional,
    Sequence,
    Union,
    cast,
)
import warnings

import matplotlib.cm as cm
import matplotlib.colors
import numpy as np

import pandas.core.common as com
from pandas.core.dtypes.common import is_list_like

if TYPE_CHECKING:
    from matplotlib.colors import Colormap


Color = Union[str, Sequence[float]]


def get_standard_colors(
    num_colors: int,
    colormap=None,
    color_type: str = "default",
    color=None,
):
    if isinstance(color, dict):
        return color

    colors = _get_colors(
        color=color,
        colormap=colormap,
        color_type=color_type,
        num_colors=num_colors,
    )

    return _cycle_colors(colors, num_colors=num_colors)


def _get_colors(
    *,
    color: Optional[Union[Color, Collection[Color]]],
    colormap: Optional[Union[str, "Colormap"]],
    color_type: str,
    num_colors: int,
) -> List[Color]:
    """Get colors from user input."""
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


def _cycle_colors(colors: List[Color], num_colors: int) -> List[Color]:
    """Append more colors by cycling if there is not enough color.

    Extra colors will be ignored by matplotlib if there are more colors
    than needed and nothing needs to be done here.
    """
    if len(colors) < num_colors:
        multiple = num_colors // len(colors) - 1
        mod = num_colors % len(colors)
        colors += multiple * colors
        colors += colors[:mod]

    return colors


def _get_colors_from_colormap(
    colormap: Union[str, "Colormap"],
    num_colors: int,
) -> List[Color]:
    """Get colors from colormap."""
    colormap = _get_cmap_instance(colormap)
    return [colormap(num) for num in np.linspace(0, 1, num=num_colors)]


def _get_cmap_instance(colormap: Union[str, "Colormap"]) -> "Colormap":
    """Get instance of matplotlib colormap."""
    if isinstance(colormap, str):
        cmap = colormap
        colormap = cm.get_cmap(colormap)
        if colormap is None:
            raise ValueError(f"Colormap {cmap} is not recognized")
    return colormap


def _get_colors_from_color(
    color: Union[Color, Collection[Color]],
) -> List[Color]:
    """Get colors from user input color."""
    if len(color) == 0:
        raise ValueError(f"Invalid color argument: {color}")

    if isinstance(color, str):
        if _is_single_color(color):
            # GH #36972
            return [color]
        else:
            return list(color)

    if _is_floats_color(color):
        color = cast(Sequence[float], color)
        return [color]

    color = cast(Collection[Color], color)
    colors = []
    for x in color:
        if _is_single_color(x):
            colors.append(x)
        else:
            raise ValueError(f"Invalid color {x}")
    return colors


def _is_floats_color(color: Union[Color, Collection[Color]]) -> bool:
    """Check if color comprises a sequence of floats representing color."""
    return bool(
        is_list_like(color)
        and (len(color) == 3 or len(color) == 4)
        and all([isinstance(x, float) for x in color])
    )


def _get_colors_from_color_type(color_type: str, num_colors: int) -> List[Color]:
    """Get colors from user input color type."""
    if color_type == "default":
        return _get_default_colors(num_colors)
    elif color_type == "random":
        return _get_random_colors(num_colors)
    else:
        raise ValueError("color_type must be either 'default' or 'random'")


def _get_default_colors(num_colors: int) -> List[Color]:
    """Get ``num_colors`` of default colors from matplotlib rc params."""
    import matplotlib.pyplot as plt

    # need to call list() on the result to copy so we don't
    # modify the global rcParams below
    try:
        colors = [c["color"] for c in list(plt.rcParams["axes.prop_cycle"])]
    except KeyError:
        colors = list(plt.rcParams.get("axes.color_cycle", list("bgrcmyk")))

    return colors[0:num_colors]


def _get_random_colors(num_colors: int) -> List[Color]:
    """Get ``num_colors`` of random colors."""
    return [_random_color(num) for num in range(num_colors)]


def _random_color(column: int) -> List[float]:
    """Get a random color represented as a list of length 3"""
    # GH17525 use common._random_state to avoid resetting the seed
    rs = com.random_state(column)
    return rs.rand(3).tolist()


def _is_single_color(color: Color) -> bool:
    """Check if ``color`` is a single color.

    Examples of single colors:
        - 'r'
        - 'g'
        - 'red'
        - 'green'
        - 'C3'

    Parameters
    ----------
    color : Color
        Color string or sequence of floats.

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
