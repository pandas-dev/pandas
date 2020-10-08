# being a bit too dynamic
from typing import Sequence, Union
import warnings

import matplotlib.cm as cm
import matplotlib.colors
import numpy as np

from pandas.core.dtypes.common import is_list_like

import pandas.core.common as com


def get_standard_colors(
    num_colors: int, colormap=None, color_type: str = "default", color=None
):
    import matplotlib.pyplot as plt

    if color is None and colormap is not None:
        if isinstance(colormap, str):
            cmap = colormap
            colormap = cm.get_cmap(colormap)
            if colormap is None:
                raise ValueError(f"Colormap {cmap} is not recognized")
        colors = [colormap(num) for num in np.linspace(0, 1, num=num_colors)]
    elif color is not None:
        if colormap is not None:
            warnings.warn(
                "'color' and 'colormap' cannot be used simultaneously. Using 'color'"
            )
        colors = (
            list(color)
            if is_list_like(color) and not isinstance(color, dict)
            else color
        )
    else:
        if color_type == "default":
            # need to call list() on the result to copy so we don't
            # modify the global rcParams below
            try:
                colors = [c["color"] for c in list(plt.rcParams["axes.prop_cycle"])]
            except KeyError:
                colors = list(plt.rcParams.get("axes.color_cycle", list("bgrcmyk")))
            if isinstance(colors, str):
                colors = list(colors)

            colors = colors[0:num_colors]
        elif color_type == "random":

            def random_color(column):
                """ Returns a random color represented as a list of length 3"""
                # GH17525 use common._random_state to avoid resetting the seed
                rs = com.random_state(column)
                return rs.rand(3).tolist()

            colors = [random_color(num) for num in range(num_colors)]
        else:
            raise ValueError("color_type must be either 'default' or 'random'")

    if isinstance(colors, str):
        conv = matplotlib.colors.ColorConverter()

        def _is_cn_color(color: str) -> bool:
            return bool(color in ["C" + str(x) for x in range(10)])

        def _is_single_color(colors: Union[str, Sequence[str]]) -> bool:
            if _is_cn_color(colors):
                return True
            try:
                conv.to_rgba(colors)
            except ValueError:
                return False
            else:
                return True

        def _is_color_cycle(colors: Union[str, Sequence[str]]) -> bool:
            if _is_cn_color(colors):
                return False
            return all([_is_single_color(c) for c in colors])

        if _is_single_color(colors) and _is_color_cycle(colors) and len(colors) > 1:
            hex_color = [c["color"] for c in list(plt.rcParams["axes.prop_cycle"])]
            colors = [hex_color[int(colors[1])]]
        elif _is_single_color(colors):
            colors = [colors]
        else:
            # ``colors`` is regarded as color cycle.
            # mpl will raise error any of them is invalid
            pass

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
