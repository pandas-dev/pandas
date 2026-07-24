from collections.abc import Callable
from typing import Any, overload
from typing_extensions import deprecated

from numpy.typing import ArrayLike, NDArray

from .utils import _Seed

@overload
def bootstrap(
    *args: ArrayLike,
    n_boot: int = 10000,
    func: str | Callable[..., Any] = "mean",
    axis: int | None = None,
    units: ArrayLike | None = None,
    seed: _Seed | None = None,
) -> NDArray[Any]: ...
@overload
@deprecated("Parameter `random_seed` is deprecated in favor of `seed`")
def bootstrap(
    *args: ArrayLike,
    n_boot: int = 10000,
    func: str | Callable[..., Any] = "mean",
    axis: int | None = None,
    units: ArrayLike | None = None,
    seed: _Seed | None = None,
    random_seed: _Seed | None = None,
) -> NDArray[Any]: ...
