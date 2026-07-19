from collections.abc import (
    Callable,
    Hashable,
    Iterable,
    Mapping,
    Sequence,
)
from typing import (
    Any,
    Literal,
)

import numpy as np

from pandas._typing import (
    ArrayLike,
    Dtype,
    ReadCsvBuffer,
    UsecolsArgType,
    npt,
)

STR_NA_VALUES: set[str]
DEFAULT_BUFFER_HEURISTIC: int

def sanitize_objects(
    values: npt.NDArray[np.object_],
    na_values: set[Hashable],
) -> int: ...

class _PendingStringColumn:
    def __len__(self) -> int: ...
    def materialize(self) -> Any: ...  # -> pyarrow.Array

class TextReader:
    unnamed_cols: set[str]
    table_width: int  # int64_t
    leading_cols: int  # int64_t
    header: list[list[int]]  # non-negative integers
    defer_pa_wrap: bool
    def __init__(
        self,
        source: ReadCsvBuffer[str] | ReadCsvBuffer[bytes],
        delimiter: bytes | str = ...,  # single-character only
        header: int | Sequence[int] | None = ...,
        header_start: int = ...,  # int64_t
        header_end: int = ...,  # uint64_t
        index_col: Hashable | Sequence[Hashable] | Literal[False] | None = ...,
        names: Sequence[Hashable] | None = ...,
        tokenize_chunksize: int = ...,  # int64_t
        delim_whitespace: bool = ...,
        converters: Mapping[Hashable, Callable] | None = ...,
        skipinitialspace: bool = ...,
        escapechar: bytes | str | None = ...,  # single-character only
        doublequote: bool = ...,
        quotechar: str | bytes | None = ...,  # at most 1 character
        quoting: int = ...,
        lineterminator: bytes | str | None = ...,  # at most 1 character
        comment: bytes | str | None = ...,  # at most 1 character
        decimal: bytes | str = ...,  # single-character only
        thousands: bytes | str | None = ...,  # single-character only
        dtype: Dtype | dict[Hashable, Dtype] = ...,
        usecols: UsecolsArgType = ...,
        error_bad_lines: bool = ...,
        warn_bad_lines: bool = ...,
        na_filter: bool = ...,
        na_values: set[str] | Mapping[Hashable, set[str]] = ...,
        na_fvalues: set[float] | Mapping[Hashable, set[float]] = ...,
        keep_default_na: bool = ...,
        true_values: list | None = ...,
        false_values: list | None = ...,
        allow_leading_cols: bool = ...,
        skiprows: int | Iterable[int] | Callable[[Hashable], bool] | None = ...,
        skipfooter: int = ...,  # int64_t
        verbose: bool = ...,
        float_precision: Literal["round_trip", "legacy", "high"] | None = ...,
        skip_blank_lines: bool = ...,
        encoding_errors: bytes | str = ...,
    ) -> None: ...
    def set_noconvert(self, i: int) -> None: ...
    def remove_noconvert(self, i: int) -> None: ...
    def close(self) -> None: ...
    def load_buffer(self, data: bytes | memoryview, strip_bom: bool = ...) -> None: ...
    def read(self, rows: int | None = ...) -> dict[int, ArrayLike]: ...
    def read_low_memory(self, rows: int | None) -> list[dict[int, ArrayLike]]: ...

# _maybe_upcast, na_values are only exposed for testing
na_values: dict[type | np.dtype, float | int]

def _maybe_upcast(
    arr: np.ndarray, use_dtype_backend: bool = ..., dtype_backend: str = ...
) -> np.ndarray: ...
