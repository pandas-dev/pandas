import datetime as dt
from _typeshed import Incomplete, SupportsGetItem
from collections.abc import Callable, Iterable, Mapping, Sequence
from typing import Any, Literal, SupportsIndex, TypeVar, overload
from typing_extensions import TypeAlias, deprecated

import numpy as np
import pandas as pd
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.legend import Legend
from matplotlib.text import Text
from matplotlib.ticker import Locator
from matplotlib.typing import ColorType
from numpy.typing import ArrayLike, NDArray
from pandas import DataFrame

from .axisgrid import Grid

__all__ = [
    "desaturate",
    "saturate",
    "set_hls_values",
    "move_legend",
    "despine",
    "get_dataset_names",
    "get_data_home",
    "load_dataset",
]

_VectorT = TypeVar("_VectorT", bound=SupportsGetItem[Any, Any])

# Type aliases used heavily throughout seaborn
_ErrorBar: TypeAlias = str | tuple[str, float] | Callable[[Iterable[float]], tuple[float, float]]  # noqa: Y047
_Estimator: TypeAlias = str | Callable[..., Incomplete]  # noqa: Y047
_Legend: TypeAlias = Literal["auto", "brief", "full"] | bool  # noqa: Y047
_LogScale: TypeAlias = bool | float | tuple[bool | float, bool | float]  # noqa: Y047
# `palette` requires dict but we use mapping to avoid a very long union because dict is invariant in its value type
_Palette: TypeAlias = str | Sequence[ColorType] | Mapping[Any, ColorType]  # noqa: Y047
_Seed: TypeAlias = int | np.random.Generator | np.random.RandomState  # noqa: Y047
_Scalar: TypeAlias = (
    # numeric
    float
    | complex
    | np.number[Any]
    # categorical
    | bool
    | str
    | bytes
    | None
    # dates
    | dt.date
    | dt.datetime
    | dt.timedelta
    | pd.Timestamp
    | pd.Timedelta
)
_Vector: TypeAlias = Iterable[_Scalar]
_DataSourceWideForm: TypeAlias = (  # noqa: Y047
    # Mapping of keys to "convertible to pd.Series" vectors
    Mapping[Any, _Vector]
    # Sequence of "convertible to pd.Series" vectors
    | Sequence[_Vector]
    # A "convertible to pd.DataFrame" table
    | Mapping[Any, Mapping[Any, _Scalar]]
    | NDArray[Any]
    # Flat "convertible to pd.Series" vector of scalars
    | Sequence[_Scalar]
)

DATASET_SOURCE: str
DATASET_NAMES_URL: str

def ci_to_errsize(cis: ArrayLike, heights: ArrayLike) -> NDArray[np.float64]: ...
def desaturate(color: ColorType, prop: float) -> tuple[float, float, float]: ...
def saturate(color: ColorType) -> tuple[float, float, float]: ...
def set_hls_values(
    color: ColorType, h: float | None = None, l: float | None = None, s: float | None = None
) -> tuple[float, float, float]: ...
@deprecated("Function `axlabel` is deprecated and will be removed in a future version")
def axlabel(xlabel: str, ylabel: str, **kwargs: Any) -> None: ...
def remove_na(vector: _VectorT) -> _VectorT: ...
def get_color_cycle() -> list[str]: ...

# `despine` should be kept roughly in line with `seaborn.axisgrid.FacetGrid.despine`
def despine(
    fig: Figure | None = None,
    ax: Axes | None = None,
    top: bool = True,
    right: bool = True,
    left: bool = False,
    bottom: bool = False,
    offset: int | Mapping[str, int] | None = None,
    trim: bool = False,
) -> None: ...
def move_legend(obj: Grid | Axes | Figure, loc: str | int, **kwargs: Any) -> None: ...
def ci(
    a: float | ArrayLike, which: float | ArrayLike = 95, axis: SupportsIndex | Sequence[SupportsIndex] | None = None
) -> NDArray[np.float64]: ...
def get_dataset_names() -> list[str]: ...
def get_data_home(data_home: str | None = None) -> str: ...
def load_dataset(name: str, cache: bool = True, data_home: str | None = None, **kws: Any) -> DataFrame: ...
def axis_ticklabels_overlap(labels: Iterable[Text]) -> bool: ...
def axes_ticklabels_overlap(ax: Axes) -> tuple[bool, bool]: ...
def locator_to_legend_entries(locator: Locator, limits: Iterable[float], dtype) -> tuple[list[Incomplete], list[str]]: ...
@overload
def relative_luminance(color: ColorType) -> float: ...  # type: ignore[overload-overlap]
@overload
def relative_luminance(color: Sequence[ColorType]) -> NDArray[np.float64]: ...
@overload
def relative_luminance(color: ColorType | Sequence[ColorType] | ArrayLike) -> float | NDArray[np.float64]: ...
def to_utf8(obj: object) -> str: ...
def adjust_legend_subtitles(legend: Legend) -> None: ...  # not public API
