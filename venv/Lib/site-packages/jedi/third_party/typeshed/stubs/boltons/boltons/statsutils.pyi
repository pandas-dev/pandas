from _typeshed import ConvertibleToFloat, Incomplete
from collections.abc import Callable, Iterable, Iterator
from typing import Any, Literal, overload
from typing_extensions import Self

class _StatsProperty:
    name: str
    func: Callable[..., Any]
    internal_name: str
    __doc__: str | None
    def __init__(self, name: str, func: Callable[..., Any]) -> None: ...
    @overload
    def __get__(self, obj: None, objtype: object = None) -> Self: ...
    @overload
    def __get__(self, obj: Stats, objtype: object = None) -> float: ...

class Stats:
    data: list[float]
    default: float
    @overload
    def __init__(self, data: list[float], default: float = 0.0, *, use_copy: Literal[False], is_sorted: bool = False) -> None: ...
    @overload
    def __init__(self, data: list[float], default: float, use_copy: Literal[False], is_sorted: bool = False) -> None: ...
    @overload
    def __init__(
        self, data: Iterable[float], default: float = 0.0, use_copy: Literal[True] = True, is_sorted: bool = False
    ) -> None: ...
    def __len__(self) -> int: ...
    def __iter__(self) -> Iterator[float]: ...
    def clear_cache(self) -> None: ...
    count: _StatsProperty
    mean: _StatsProperty
    max: _StatsProperty
    min: _StatsProperty
    median: _StatsProperty
    iqr: _StatsProperty
    trimean: _StatsProperty
    variance: _StatsProperty
    std_dev: _StatsProperty
    median_abs_dev: _StatsProperty
    mad: _StatsProperty
    rel_std_dev: _StatsProperty
    skewness: _StatsProperty
    kurtosis: _StatsProperty
    pearson_type: _StatsProperty
    def get_quantile(self, q: ConvertibleToFloat) -> float: ...
    def get_zscore(self, value: float) -> float: ...
    def trim_relative(self, amount: float = 0.15) -> None: ...
    def get_histogram_counts(self, bins: int | list[float] | None = None, **kw) -> list[tuple[float, int]]: ...
    def format_histogram(self, bins: int | list[float] | None = None, **kw) -> str: ...
    def describe(
        self, quantiles: Iterable[float] | None = None, format: str | None = None
    ) -> dict[str, float] | list[tuple[str, float]] | str: ...

def describe(
    data: Iterable[float], quantiles: Iterable[float] | None = None, format: str | None = None
) -> dict[str, float] | list[tuple[str, float]] | str: ...

mean: Incomplete
median: Incomplete
iqr: Incomplete
trimean: Incomplete
variance: Incomplete
std_dev: Incomplete
median_abs_dev: Incomplete
rel_std_dev: Incomplete
skewness: Incomplete
kurtosis: Incomplete
pearson_type: Incomplete

def format_histogram_counts(
    bin_counts: list[float], width: int | None = None, format_bin: Callable[..., Any] | None = None
) -> str: ...
