from collections.abc import Callable

import numpy as np

from pandas._typing import ArrayLike

def write_csv_rows(
    data: list[ArrayLike],
    data_index: np.ndarray,
    nlevels: int,
    cols: np.ndarray,
    writer: object,  # _csv.writer
) -> None: ...
def convert_json_to_lines(arr: str) -> str: ...
def build_json_rows(
    columns: list[list],
    nrows: int,
    keys: list | None = ...,
) -> list: ...
def json_zip_object(keys: bytes, values: bytes) -> bytes: ...
def json_scan(text: bytes) -> tuple[int, int, int, bool]: ...
def json_nan_to_null(text: bytes, nan_count: int) -> bytes: ...
def json_quote_elements(values: bytes) -> bytes: ...
def json_object_of_fragments(keys: bytes, values: list[bytes]) -> bytes: ...
def json_zip_records(keys: bytes, rows: bytes, labels: bytes | None = ...) -> bytes: ...
def normalize_json_objects(
    arr: np.ndarray,  # np.ndarray[object, ndim=1]
    fallback: Callable[[object], object],
    nan_to_none: bool,
    iso_dates: bool,
    date_unit: str,
) -> list: ...
def datetime64_to_json(
    values: np.ndarray,  # np.ndarray[np.int64]
    value_unit: str,
    out_unit: str,
    utc: bool,
    nat_as_string: bool = ...,
) -> bytes: ...
def timedelta64_to_json(
    values: np.ndarray,  # np.ndarray[np.int64]
    unit: str,
    nat_as_string: bool = ...,
) -> bytes: ...
def int64_to_json(
    values: np.ndarray,  # np.ndarray[np.int64]
    nat_as_string: bool = ...,
) -> bytes: ...
def json_transpose(
    columns: list[bytes],
    keys: bytes | None = ...,
    labels: bytes | None = ...,
) -> bytes: ...
def max_len_string_array(
    arr: np.ndarray,  # pandas_string[:]
) -> int: ...
def word_len(val: object) -> int: ...
def string_array_replace_from_nan_rep(
    arr: np.ndarray,  # np.ndarray[object, ndim=1]
    nan_rep: object,
) -> None: ...
