from typing import Literal, overload

from matplotlib.colors import LinearSegmentedColormap

__all__ = [
    "choose_colorbrewer_palette",
    "choose_cubehelix_palette",
    "choose_dark_palette",
    "choose_light_palette",
    "choose_diverging_palette",
]

@overload
def choose_colorbrewer_palette(
    data_type: Literal["sequential", "diverging", "qualitative"], as_cmap: Literal[True]
) -> LinearSegmentedColormap: ...
@overload
def choose_colorbrewer_palette(
    data_type: Literal["sequential", "diverging", "qualitative"], as_cmap: Literal[False] = False
) -> list[tuple[float, float, float]]: ...
@overload
def choose_dark_palette(input: str = "husl", *, as_cmap: Literal[True]) -> LinearSegmentedColormap: ...
@overload
def choose_dark_palette(input: str = "husl", as_cmap: Literal[False] = False) -> list[tuple[float, float, float]]: ...
@overload
def choose_light_palette(input: str = "husl", *, as_cmap: Literal[True]) -> LinearSegmentedColormap: ...
@overload
def choose_light_palette(input: str = "husl", as_cmap: Literal[False] = False) -> list[tuple[float, float, float]]: ...
@overload
def choose_diverging_palette(as_cmap: Literal[True]) -> LinearSegmentedColormap: ...
@overload
def choose_diverging_palette(as_cmap: Literal[False] = False) -> list[tuple[float, float, float]]: ...
@overload
def choose_cubehelix_palette(as_cmap: Literal[True]) -> LinearSegmentedColormap: ...
@overload
def choose_cubehelix_palette(as_cmap: Literal[False] = False) -> list[tuple[float, float, float]]: ...
