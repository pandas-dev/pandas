from typing import (
    Any,
    Tuple,
    Union,
)

import numpy as np

from pandas._typing import ArrayLike

# placeholder until we can specify np.ndarray[object, ndim=2]
ndarray_obj_2d = np.ndarray

no_default: object

def item_from_zerodim(val: object) -> object: ...
def infer_dtype(value: object, skipna: bool = True) -> str: ...

def is_iterator(obj: object) -> bool: ...
def is_scalar(val: object) -> bool: ...
def is_list_like(obj: object, allow_sets: bool = True) -> bool: ...

def is_period(val: object) -> bool: ...
def is_interval(val: object) -> bool: ...
def is_decimal(val: object) -> bool: ...
def is_complex(val: object) -> bool: ...
def is_bool(val: object) -> bool: ...
def is_integer(val: object) -> bool: ...
def is_float(val: object) -> bool: ...

def is_interval_array(values: np.ndarray) -> bool: ...
def is_period_array(values: np.ndarray) -> bool: ...
def is_datetime64_array(values: np.ndarray) -> bool: ...
def is_timedelta_or_timedelta64_array(values: np.ndarray) -> bool: ...
def is_datetime_with_singletz_array(values: np.ndarray) -> bool: ...

def is_time_array(values: np.ndarray, skipna: bool = False): ...
def is_date_array(values: np.ndarray, skipna: bool = False): ...
def is_datetime_array(values: np.ndarray, skipna: bool = False): ...
def is_string_array(values: np.ndarray, skipna: bool = False): ...
def is_float_array(values: np.ndarray, skipna: bool = False): ...
def is_integer_array(values: np.ndarray, skipna: bool = False): ...
def is_bool_array(values: np.ndarray, skipna: bool = False): ...

def fast_multiget(mapping: dict, keys: np.ndarray, default=np.nan) -> ArrayLike: ...

# TODO: gen: Generator?
def fast_unique_multiple_list_gen(gen: object, sort: bool = True) -> list: ...
def fast_unique_multiple_list(lists: list, sort: bool = True) -> list: ...
def fast_unique_multiple(arrays: list, sort: bool = True) -> list: ...

# TODO: f: Callable?
def map_infer(
    arr: np.ndarray, f: object, convert: bool = True, ignore_na: bool = False
) -> ArrayLike: ...

def maybe_convert_objects(
    objects: np.ndarray[object],
    try_float: bool = False,
    safe: bool = False,
    convert_datetime: bool = False,
    convert_timedelta: bool = False,
    convert_to_nullable_integer: bool = False,
) -> ArrayLike: ...

def maybe_convert_numeric(
    values: np.ndarray[object],
    na_values: set,
    convert_empty: bool = True,
    coerce_numeric: bool = False,
) -> np.ndarray: ...

# TODO: restrict `arr`?
def ensure_string_array(
    arr,
    na_value: object = np.nan,
    convert_na_value: bool = True,
    copy: bool = True,
    skipna: bool = True,
) -> np.ndarray[object]: ...

def infer_datetimelike_array(arr: np.ndarray[object]) -> str: ...

# TODO: new_dtype -> np.dtype?
def astype_intsafe(arr: np.ndarray[object], new_dtype) -> np.ndarray: ...

def fast_zip(ndarrays: list) -> np.ndarray[object]: ...

# TODO: can we be more specific about rows?
def to_object_array_tuples(rows: object) -> ndarray_obj_2d: ...

def tuples_to_object_array(tuples: np.ndarray[object]) -> ndarray_obj_2d: ...

# TODO: can we be more specific about rows?
def to_object_array(rows: object, min_width: int = 0) -> ndarray_obj_2d: ...

def dicts_to_array(dicts: list, columns: list) -> ndarray_obj_2d: ...


def maybe_booleans_to_slice(
    mask: np.ndarray[np.uint8]
) -> Union[slice, np.ndarray[np.uint8]]: ...

def maybe_indices_to_slice(
    indices: np.ndarray[np.intp], max_len: int
) -> Union[slice, np.ndarray[np.intp]]: ...

def clean_index_list(obj: list) -> Tuple[Union[list, np.ndarray], bool]: ...


# -----------------------------------------------------------------
# Functions for which we lie, using ndarray[foo] in place of foo[:]

# actually object[:]
def memory_usage_of_objects(arr: np.ndarray[object]) -> np.int64: ...


# TODO: f: Callable?
# TODO: dtype -> DtypeObj?
def map_infer_mask(
    arr: np.ndarray,
    f: object,
    mask: np.ndarray[np.uint8],  # actually const uint8_t[:]
    convert: bool = True,
    na_value: Any = no_default,
    dtype: Any = object,
) -> ArrayLike: ...

def indices_fast(
    index: np.ndarray,
    labels: np.ndarray[np.int64],  # actually const int64_t[:]
    keys: list,
    sorted_labels: list,
) -> dict: ...

def generate_slices(
    labels: np.ndarray[np.int64],  # actually const int64_t[:]
    ngroups: int
) -> Tuple[np.ndarray[np.int64], np.ndarray[np.int64]]: ...

def count_level_2d(
    mask: np.ndarray[np.uint8],  # actually ndarray[uint8_t, ndim=2, cast=True],
    labels: np.ndarray[np.int64],  # actually const int64_t[:]
    max_bin: int,
    axis: int
) -> np.ndarray[np.int64]: ...   # actually np.ndarray[np.int64, ndim=2]

def get_level_sorter(
    label: np.ndarray[np.int64],  # actually const int64_t[:]
    starts: np.ndarray[np.int64],  # actually const int64_t[:]
) -> np.ndarray[np.int64]: ...  # actually np.ndarray[np.int64, ndim=1]


def generate_bins_dt64(
    values: np.ndarray[np.int64],
    binner: np.ndarray[np.int64],  # actually const int64_t[:]
    closed: object = "left",
    hasnans: bool = False,
) -> np.ndarray[np.int64]: ...  # actually np.ndarray[np.int64, ndim=1]


# actually object[:] for both args
def array_equivalent_object(
    left: np.ndarray[object], right: np.ndarray[object]
) -> bool: ...

# actually const float64_t[:]
def has_infs_f8(arr: np.ndarray[np.float64]) -> bool: ...

# actually const float32_t[:]
def has_infs_f4(arr: np.ndarray[np.float32]) -> bool: ...

def get_reverse_indexer(
    indexer: np.ndarray[np.int64],  # actually const int64_t[:]
    length: int,
) -> np.ndarray[np.int64]: ...
