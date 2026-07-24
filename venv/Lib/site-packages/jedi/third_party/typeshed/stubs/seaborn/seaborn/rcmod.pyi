from _typeshed import Unused
from collections.abc import Callable, Sequence
from typing import Any, Literal, TypeVar
from typing_extensions import deprecated

from matplotlib.typing import ColorType

__all__ = [
    "set_theme",
    "set",
    "reset_defaults",
    "reset_orig",
    "axes_style",
    "set_style",
    "plotting_context",
    "set_context",
    "set_palette",
]

_KT = TypeVar("_KT")
_VT = TypeVar("_VT")
_F = TypeVar("_F", bound=Callable[..., Any])

def set_theme(
    context: Literal["paper", "notebook", "talk", "poster"] | dict[str, Any] = "notebook",
    style: Literal["white", "dark", "whitegrid", "darkgrid", "ticks"] | dict[str, Any] = "darkgrid",
    palette: str | Sequence[ColorType] | None = "deep",
    font: str = "sans-serif",
    font_scale: float = 1,
    color_codes: bool = True,
    rc: dict[str, Any] | None = None,
) -> None: ...
@deprecated("Function `set` is deprecated in favor of `set_theme`")
def set(
    context: Literal["paper", "notebook", "talk", "poster"] | dict[str, Any] = "notebook",
    style: Literal["white", "dark", "whitegrid", "darkgrid", "ticks"] | dict[str, Any] = "darkgrid",
    palette: str | Sequence[ColorType] | None = "deep",
    font: str = "sans-serif",
    font_scale: float = 1,
    color_codes: bool = True,
    rc: dict[str, Any] | None = None,
) -> None: ...
def reset_defaults() -> None: ...
def reset_orig() -> None: ...
def axes_style(
    style: Literal["white", "dark", "whitegrid", "darkgrid", "ticks"] | dict[str, Any] | None = None,
    rc: dict[str, Any] | None = None,
) -> _AxesStyle[str, Any]: ...
def set_style(
    style: Literal["white", "dark", "whitegrid", "darkgrid", "ticks"] | dict[str, Any] | None = None,
    rc: dict[str, Any] | None = None,
) -> None: ...
def plotting_context(
    context: Literal["paper", "notebook", "talk", "poster"] | dict[str, Any] | None = None,
    font_scale: float = 1,
    rc: dict[str, Any] | None = None,
) -> _PlottingContext[str, Any]: ...
def set_context(
    context: Literal["paper", "notebook", "talk", "poster"] | dict[str, Any] | None = None,
    font_scale: float = 1,
    rc: dict[str, Any] | None = None,
) -> None: ...

class _RCAesthetics(dict[_KT, _VT]):
    def __enter__(self) -> None: ...
    def __exit__(self, exc_type: Unused, exc_value: Unused, exc_tb: Unused) -> None: ...
    def __call__(self, func: _F) -> _F: ...

class _AxesStyle(_RCAesthetics[_KT, _VT]): ...
class _PlottingContext(_RCAesthetics[_KT, _VT]): ...

def set_palette(
    palette: str | Sequence[ColorType] | None, n_colors: int | None = None, desat: float | None = None, color_codes: bool = False
) -> None: ...
