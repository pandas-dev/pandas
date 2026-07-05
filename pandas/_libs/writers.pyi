import numpy as np

from pandas._typing import ArrayLike

CSV_KIND_OBJ: int
CSV_KIND_FLOAT64: int
CSV_KIND_FLOAT32: int
CSV_KIND_INT64: int
CSV_KIND_UINT64: int
CSV_KIND_BOOL: int
CSV_KIND_DT64: int
FLOAT32_NATIVE: bool

def write_csv_rows(
    data: list[ArrayLike],
    data_index: np.ndarray,
    nlevels: int,
    cols: np.ndarray,
    writer: object,  # _csv.writer
) -> None: ...
def write_csv_chunk(
    cols: list[tuple[int, np.ndarray, int, int, bool]],
    nrows: int,
    sep: str,
    quotechar: str,
    lineterminator: str,
    na_rep: str,
    crlf_always: bool,
) -> str: ...
def convert_json_to_lines(arr: str) -> str: ...
def max_len_string_array(
    arr: np.ndarray,  # pandas_string[:]
) -> int: ...
def word_len(val: object) -> int: ...
def string_array_replace_from_nan_rep(
    arr: np.ndarray,  # np.ndarray[object, ndim=1]
    nan_rep: object,
) -> None: ...
