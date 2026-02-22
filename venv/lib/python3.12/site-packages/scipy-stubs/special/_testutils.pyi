__all__ = ["FuncData", "assert_func_equal", "with_special_errors"]

from collections.abc import Callable, Generator, Sequence
from types import ModuleType
from typing import Any, Final, Generic, TypeAlias, overload
from typing_extensions import TypeVar

import numpy as np
import optype.numpy as onp
import pytest

###

_FuncT = TypeVar("_FuncT", bound=Callable[..., object])
_ScalarT_co = TypeVar("_ScalarT_co", bound=np.generic, default=Any, covariant=True)

_ResultFunc: TypeAlias = Callable[..., object]
_FilterFunc: TypeAlias = Callable[..., bool]

###

class MissingModule:  # undocumented
    name: Final[str]
    def __init__(self, /, name: str) -> None: ...

#
def check_version(module: ModuleType | MissingModule, min_ver: str) -> pytest.MarkDecorator: ...  # undocumented

#
def with_special_errors(func: _FuncT) -> _FuncT: ...  # undocumented

#
def assert_func_equal(
    func: _ResultFunc,
    results: _ResultFunc | onp.ToArrayND,
    points: Generator[object] | onp.ToArrayND,
    rtol: float | None = None,
    atol: float | None = None,
    param_filter: Callable[..., bool] | tuple[Callable[..., bool] | None, ...] | None = None,
    knownfailure: str | None = None,
    vectorized: bool = True,
    dtype: onp.ToDType | None = None,  # unused
    nan_ok: bool = False,
    ignore_inf_sign: bool = False,
    distinguish_nan_and_inf: bool = True,
) -> None: ...  # undocumented

#
class FuncData(Generic[_ScalarT_co]):  # undocumented
    func: Final[_ResultFunc]
    data: onp.Array2D[_ScalarT_co]
    dataname: Final[str]
    param_columns: tuple[int, ...]
    result_columns: tuple[int, ...]
    result_func: Final[_ResultFunc | None]
    rtol: Final[float | None]
    atol: Final[float | None]
    param_filter: Final[tuple[_FilterFunc | None, ...]]
    knownfailure: Final[str | None]
    nan_ok: Final[bool]
    vectorized: Final[bool]
    ignore_inf_sign: Final[bool]
    distinguish_nan_and_inf: Final[bool]

    @overload  # result_colums, no result_func
    def __init__(
        self,
        /,
        func: _ResultFunc,
        data: onp.Array2D[_ScalarT_co],
        param_columns: int | Sequence[int],
        result_columns: int | tuple[int, ...] | None = None,
        result_func: None = None,
        rtol: float | None = None,
        atol: float | None = None,
        param_filter: _FilterFunc | tuple[_FilterFunc | None, ...] | None = None,
        knownfailure: str | None = None,
        dataname: str | None = None,
        nan_ok: bool = False,
        vectorized: bool = True,
        ignore_inf_sign: bool = False,
        distinguish_nan_and_inf: bool = True,
    ) -> None: ...
    @overload  # no result_columns, result_func
    def __init__(
        self,
        /,
        func: _ResultFunc,
        data: onp.Array2D[_ScalarT_co],
        param_columns: int | Sequence[int],
        result_columns: None = None,
        result_func: _ResultFunc | None = None,
        rtol: float | None = None,
        atol: float | None = None,
        param_filter: _FilterFunc | tuple[_FilterFunc | None, ...] | None = None,
        knownfailure: str | None = None,
        dataname: str | None = None,
        nan_ok: bool = False,
        vectorized: bool = True,
        ignore_inf_sign: bool = False,
        distinguish_nan_and_inf: bool = True,
    ) -> None: ...

    #
    def get_tolerances(self, /, dtype: onp.ToDType) -> tuple[float, float]: ...
    def check(
        self, /, data: onp.ArrayND | None = None, dtype: onp.ToDType | None = None, dtypes: Sequence[onp.ToDType] | None = None
    ) -> None: ...
