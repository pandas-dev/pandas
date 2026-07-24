from _typeshed import Incomplete
from collections.abc import Callable, Hashable, Mapping
from typing_extensions import Concatenate, ParamSpec, TypeAlias

from numpy import ufunc
from pandas import DataFrame

# pandas._typing.AggFuncTypeFrame is partially Unknown
_AggFuncTypeBase: TypeAlias = Callable[..., Incomplete] | str | ufunc
_AggFuncTypeDictFrame: TypeAlias = Mapping[Hashable, _AggFuncTypeBase | list[_AggFuncTypeBase]]
_AggFuncTypeFrame: TypeAlias = _AggFuncTypeBase | list[_AggFuncTypeBase] | _AggFuncTypeDictFrame

_P = ParamSpec("_P")

class GroupBy:
    order: dict[str, list[Incomplete] | None]
    def __init__(self, order: list[str] | dict[str, list[Incomplete] | None]) -> None: ...
    # Signature based on pandas.core.groupby.generic.DataFrameGroupBy.aggregate
    # args and kwargs possible values depend on func which itself can be
    # an attribute name, a mapping, a callable, or lead to a jitted numba function
    def agg(
        self,
        data: DataFrame,
        func: _AggFuncTypeFrame = ...,
        *args,
        engine: str | None = None,
        engine_kwargs: dict[str, bool] | None = None,
        **kwargs,
    ) -> DataFrame: ...
    def apply(
        self, data: DataFrame, func: Callable[Concatenate[DataFrame, _P], DataFrame], *args: _P.args, **kwargs: _P.kwargs
    ) -> DataFrame: ...
