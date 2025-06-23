"""
Factory functions for creating pandas dtypes with consistent naming conventions.
"""

from __future__ import annotations

from typing import (
    Any,
    Literal,
)

import numpy as np

from pandas._libs import missing as libmissing

from pandas.core.dtypes.dtypes import (
    CategoricalDtype,
    DatetimeTZDtype,
    IntervalDtype,
    PeriodDtype,
    SparseDtype,
)

from pandas.core.api import (
    ArrowDtype,
    BooleanDtype,
    Float32Dtype,
    Float64Dtype,
    Int8Dtype,
    Int16Dtype,
    Int32Dtype,
    Int64Dtype,
)
from pandas.core.arrays.string_ import StringDtype


def string(
    backend: Literal["python", "pyarrow"] = "python",
    large: bool = False,
    mode: Literal["string", "binary"] = "string",
    na_value: Any = libmissing.NA,
) -> StringDtype | ArrowDtype:
    """
    Create a string or binary dtype with specified backend and NA value.

    Parameters
    ----------
    backend : {"python", "pyarrow"}, default "python"
        The backend to use for string or binary data.
    large : bool, default False
        If True and backend is "pyarrow", uses pa.large_string() or pa.large_binary().
    mode : {"string", "binary"}, default "string"
        Whether to create a string or binary dtype.
    na_value : {pd.NA, np.nan}, default pd.NA
        The value to use for missing data. Must be either pd.NA or np.nan.

    Returns
    -------
    StringDtype or ArrowDtype
        A string or binary dtype with the specified configuration.

    Examples
    --------
    >>> print(string())  # Default python backend with pd.NA
    string
    >>> print(string(backend="pyarrow", mode="string"))  # PyArrow string backend
    string[pyarrow]
    >>> print(
    ...     string(backend="pyarrow", mode="string", large=True)
    ... )  # PyArrow large string
    large_string[pyarrow]
    >>> print(string(backend="pyarrow", mode="binary"))  # PyArrow binary
    binary[pyarrow]
    >>> print(
    ...     string(backend="pyarrow", mode="binary", large=True)
    ... )  # PyArrow large binary
    large_binary[pyarrow]
    """
    valid_modes = ["string", "binary"]
    if mode not in valid_modes:
        raise ValueError(f"mode must be one of {valid_modes}, got {mode}")
    if backend == "pyarrow":
        import pyarrow as pa

        if mode == "string":
            pa_type = pa.large_string() if large else pa.string()
        else:  # mode == "binary"
            pa_type = pa.large_binary() if large else pa.binary()
        return ArrowDtype(pa_type)
    if mode == "binary":
        raise ValueError("Binary mode is only supported with PyArrow backend.")
    return StringDtype(storage="python", na_value=na_value)


def datetime(
    unit: str = "ns",
    tz: Any | None = None,
    backend: Literal["numpy", "pyarrow"] = "numpy",
) -> np.dtype | DatetimeTZDtype | ArrowDtype:
    """
    Create a datetime dtype with specified unit, timezone and backend.

    Parameters
    ----------
    unit : str, default "ns"
        The datetime unit to use.
    tz : str, int, or datetime.tzinfo, optional
        The timezone to use.
    backend : {"numpy", "pyarrow"}, default "numpy"
        The backend to use for datetime storage.

    Returns
    -------
    Union[np.dtype, DatetimeTZDtype, ArrowDtype]
        A datetime dtype with the specified configuration.

    Examples
    --------
    >>> print(datetime())  # Default numpy backend with ns unit
    datetime64[ns]
    >>> print(datetime(unit="us"))  # Microsecond precision
    datetime64[us]
    >>> print(datetime(tz="UTC"))  # Timezone-aware datetime
    datetime64[ns, UTC]
    >>> print(datetime(backend="pyarrow"))  # PyArrow backend
    timestamp[ns][pyarrow]
    """
    valid_units = ["D", "h", "m", "s", "ms", "us", "ns"]
    if backend == "numpy":
        if unit not in valid_units:
            raise ValueError(f"unit must be one of {valid_units}, got {unit}")
        if tz is not None:
            return DatetimeTZDtype(unit=unit, tz=tz)
        return np.dtype(f"datetime64[{unit}]")
    else:  # pyarrow
        import pyarrow as pa

        return ArrowDtype(pa.timestamp(unit, tz=tz))


def integer(
    bits: int = 64,
    backend: Literal["numpy", "pyarrow", "pandas"] = "pandas",
) -> Int8Dtype | Int16Dtype | Int32Dtype | Int64Dtype | ArrowDtype | np.dtype[Any]:
    """
    Create an integer dtype with specified bits and backend.

    Parameters
    ----------
    bits : int, default 64
        The number of bits for the integer type. Must be one of 8, 16, 32, or 64.
    backend : {"pandas", "numpy", "pyarrow"}, default "pandas"
        The backend to use for integer storage.

    Returns
    -------
    Union[Int8Dtype, Int16Dtype, Int32Dtype, Int64Dtype, ArrowDtype]
        An integer dtype with the specified configuration.

    Examples
    --------
    >>> print(integer())  # Default: 64 bits with pandas backend
    Int64
    >>> print(integer(bits=32))  # 32-bit integer with pandas backend
    Int32
    >>> print(integer(bits=64, backend="numpy"))  # 64-bit integer with NumPy backend
    int64
    >>> print(
    ...     integer(bits=64, backend="pyarrow")
    ... )  # 64-bit integer with PyArrow backend
    int64[pyarrow]
    """
    valid_bits = [8, 16, 32, 64]
    if bits not in valid_bits:
        raise ValueError(f"bits must be one of {valid_bits}, got {bits}")

    if backend == "numpy":
        return np.dtype(f"int{bits}")
    elif backend == "pandas":
        if bits == 8:
            return Int8Dtype()
        elif bits == 16:
            return Int16Dtype()
        elif bits == 32:
            return Int32Dtype()
        else:  # bits == 64
            return Int64Dtype()
    elif backend == "pyarrow":
        import pyarrow as pa

        if bits == 8:
            return ArrowDtype(pa.int8())
        elif bits == 16:
            return ArrowDtype(pa.int16())
        elif bits == 32:
            return ArrowDtype(pa.int32())
        else:  # bits == 64
            return ArrowDtype(pa.int64())
    else:
        raise ValueError(f"Unsupported backend: {backend!r}")


def floating(
    bits: int = 64,
    backend: Literal["numpy", "pyarrow", "pandas"] = "pandas",
) -> Float32Dtype | Float64Dtype | ArrowDtype | np.dtype[Any]:
    """
    Create a floating-point dtype with specified bits and backend.

    Parameters
    ----------
    bits : int, default 64
        The number of bits for the floating-point type. Must be one of 32 or 64.
    backend : {"numpy", "pyarrow"}, default "numpy"
        The backend to use for floating-point storage.

    Returns
    -------
    Union[Float32Dtype, Float64Dtype, ArrowDtype]
        A floating-point dtype with the specified configuration.

    Examples
    --------
    >>> print(floating())  # Default: 64 bits with NumPy backend
    Float64
    >>> print(floating(bits=32))  # 32-bit float with NumPy backend
    Float32
    >>> print(floating(bits=64, backend="pyarrow"))  # 64-bit float with PyArrow backend
    double[pyarrow]
    """
    valid_bits = [32, 64]
    if bits not in valid_bits:
        raise ValueError(f"bits must be one of {valid_bits}, got {bits}")

    if backend == "numpy":
        return np.dtype(f"float{bits}")
    elif backend == "pandas":
        if bits == 32:
            return Float32Dtype()
        else:  # bits == 64
            return Float64Dtype()
    elif backend == "pyarrow":
        import pyarrow as pa

        if bits == 32:
            return ArrowDtype(pa.float32())
        else:  # bits == 64
            return ArrowDtype(pa.float64())
    else:
        raise ValueError(f"Unsupported backend: {backend!r}")


def decimal(
    precision: int,
    scale: int,
    backend: Literal["pyarrow"] = "pyarrow",
) -> ArrowDtype:
    """
    Create a decimal dtype with specified precision and scale.

    Parameters
    ----------
    precision : int
        The total number of digits in the decimal number.
    scale : int
        The number of digits to the right of the decimal point.
    backend : {"pyarrow"}, default "pyarrow"
        The backend to use for decimal storage. Only PyArrow is supported.

    Returns
    -------
    ArrowDtype
        A decimal dtype with the specified configuration.

    Examples
    --------
    >>> print(decimal(precision=10, scale=2))  # Decimal with 10 digits,
    ... # 2 after the decimal point
    decimal128(10, 2)[pyarrow]
    >>> print(decimal(precision=40, scale=5))  # Larger precision, uses decimal256
    decimal256(40, 5)[pyarrow]
    """
    if backend == "pyarrow":
        import pyarrow as pa

        if precision <= 38:
            return ArrowDtype(pa.decimal128(precision, scale))
        return ArrowDtype(pa.decimal256(precision, scale))
    raise ValueError("Decimal types are only supported with PyArrow backend.")


def boolean(
    backend: Literal["numpy", "pyarrow"] = "numpy",
) -> BooleanDtype | ArrowDtype:
    """
    Create a boolean dtype with specified backend.

    Parameters
    ----------
    backend : {"numpy", "pyarrow"}, default "numpy"
        The backend to use for boolean storage.

    Returns
    -------
    Union[BooleanDtype, ArrowDtype]
        A boolean dtype with the specified configuration.

    Examples
    --------
    >>> print(boolean())  # Default: NumPy backend
    boolean
    >>> print(boolean(backend="pyarrow"))  # PyArrow backend
    bool[pyarrow]
    """
    if backend == "numpy":
        return BooleanDtype()
    else:  # pyarrow
        import pyarrow as pa

        return ArrowDtype(pa.bool_())


def list_dtype(
    value_type: Any = None,
    large: bool = False,
    backend: Literal["numpy", "pyarrow"] = "numpy",
) -> np.dtype | ArrowDtype:
    """
    Create a list dtype with specified value type, size, and backend.

    Parameters
    ----------
    value_type : Any, optional
        The type of the list elements (e.g., pa.int64(), pa.string()). If None,
        defaults to object (NumPy) or int64 (PyArrow).
    large : bool, default False
        If True and backend is "pyarrow", uses pa.large_list() instead of pa.list_().
    backend : {"numpy", "pyarrow"}, default "numpy"
        The backend to use for list storage.

    Returns
    -------
    Union[np.dtype, ArrowDtype]
        A list dtype with the specified configuration.

    Examples
    --------
    >>> print(list_dtype())  # Default numpy backend
    object
    >>> print(list_dtype(backend="pyarrow"))  # PyArrow backend with default int64
    list<item: int64>[pyarrow]
    >>> import pyarrow as pa
    >>> print(
    ...     list_dtype(value_type=pa.string(), backend="pyarrow")
    ... )  # PyArrow with string
    list<item: string>[pyarrow]
    >>> print(
    ...     list_dtype(value_type=pa.string(), large=True, backend="pyarrow")
    ... )  # PyArrow large list
    large_list<item: string>[pyarrow]
    """
    if backend == "numpy":
        return np.dtype("object")
    else:  # pyarrow
        import pyarrow as pa

        if value_type is None:
            value_type = pa.int64()
        pa_type = pa.large_list(value_type) if large else pa.list_(value_type)
        return ArrowDtype(pa_type)


def categorical(
    categories: list[Any] | None = None,
    ordered: bool = False,
    index_type: Any = None,
    value_type: Any = None,
    backend: Literal["numpy", "pyarrow"] = "numpy",
) -> CategoricalDtype | ArrowDtype:
    """
    Create a categorical dtype with specified categories, ordering, and backend.

    Parameters
    ----------
    categories : list, optional
        The categories for the categorical dtype.
    ordered : bool, default False
        Whether the categories are ordered.
    index_type : Any, optional
        The type of the dictionary indices (PyArrow only, e.g., pa.int32()).
        Defaults to pa.int32() if None.
    value_type : Any, optional
        The type of the dictionary values (PyArrow only, e.g., pa.string()).
        Defaults to pa.string() if None.
    backend : {"numpy", "pyarrow"}, default "numpy"
        The backend to use for categorical storage.

    Returns
    -------
    Union[CategoricalDtype, ArrowDtype]
        A categorical dtype with the specified configuration.

    Examples
    --------
    >>> print(categorical())  # Default numpy backend
    category
    >>> print(categorical(categories=["a", "b", "c"]))  # With categories
    category
    >>> print(categorical(ordered=True))  # Ordered categories
    category
    >>> import pyarrow as pa
    >>> print(categorical(backend="pyarrow"))  # PyArrow backend
    dictionary<values=string, indices=int32, ordered=0>[pyarrow]
    >>> print(
    ...     categorical(index_type=pa.int64(), value_type=pa.int32(), backend="pyarrow")
    ... )
    dictionary<values=int32, indices=int64, ordered=0>[pyarrow]
    """
    if backend == "numpy":
        return CategoricalDtype(categories=categories, ordered=ordered)
    else:  # pyarrow
        import pyarrow as pa

        index_type = pa.int32() if index_type is None else index_type
        value_type = pa.string() if value_type is None else value_type
        return ArrowDtype(pa.dictionary(index_type, value_type))


def interval(
    subtype: Any = None,
    closed: Literal["left", "right", "both", "neither"] = "right",
    backend: Literal["numpy", "pyarrow"] = "numpy",
) -> IntervalDtype | ArrowDtype:
    """
    Create an interval dtype with specified subtype and closed bounds.

    Parameters
    ----------
    subtype : dtype, optional
        The dtype of the interval bounds.
    closed : {"left", "right", "both", "neither"}, default "right"
        Whether the interval is closed on the left, right, both or neither.
    backend : {"numpy", "pyarrow"}, default "numpy"
        The backend to use for interval storage.

    Returns
    -------
    Union[IntervalDtype, ArrowDtype]
        An interval dtype with the specified configuration.

    Examples
    --------
    >>> print(interval())  # Default numpy backend
    interval
    >>> print(interval(subtype="int64"))  # With specific subtype
    interval[int64, right]
    >>> print(interval(closed="both"))  # Closed on both sides
    interval
    >>> print(interval(backend="pyarrow"))  # PyArrow backend
    struct<left: double, right: double>[pyarrow]
    """
    if backend == "numpy":
        return IntervalDtype(subtype=subtype, closed=closed)
    else:  # pyarrow
        import pyarrow as pa

        if subtype is not None:
            return ArrowDtype(
                pa.struct(
                    [
                        ("left", pa.from_numpy_dtype(subtype)),
                        ("right", pa.from_numpy_dtype(subtype)),
                    ]
                )
            )
        return ArrowDtype(pa.struct([("left", pa.float64()), ("right", pa.float64())]))


def period(
    freq: str = "D",
    backend: Literal["numpy", "pyarrow"] = "numpy",
) -> PeriodDtype | ArrowDtype:
    """
    Create a period dtype with specified frequency.

    Parameters
    ----------
    freq : str, default "D"
        The frequency of the period. Common values are:
        - "D" for daily
        - "M" for monthly
        - "Y" for yearly
        - "H" for hourly
        - "T" for minute
        - "S" for second
    backend : {"numpy", "pyarrow"}, default "numpy"
        The backend to use for period storage.

    Returns
    -------
    Union[PeriodDtype, ArrowDtype]
        A period dtype with the specified configuration.

    Notes
    -----
    PyArrow backend uses `month_day_nano_interval` for periods, which represents
    intervals in terms of months, days, and nanoseconds.

    Examples
    --------
    >>> print(period())  # Default numpy backend with daily frequency
    period[D]
    >>> print(period(freq="M"))  # Monthly frequency
    period[M]
    >>> print(period(backend="pyarrow"))  # PyArrow backend
    month_day_nano_interval[pyarrow]
    """
    if backend == "numpy":
        return PeriodDtype(freq=freq)
    else:  # pyarrow
        import pyarrow as pa

        return ArrowDtype(pa.month_day_nano_interval())


def sparse(
    dtype: Any = None,
    fill_value: Any = None,
    backend: Literal["numpy"] = "numpy",
) -> np.dtype | SparseDtype:
    """
    Create a sparse dtype with specified dtype and fill value.

    Parameters
    ----------
    dtype : dtype, optional
        The dtype of the non-sparse values. If None, defaults to float64.
    fill_value : scalar, optional
        The value to use for missing values. If None, defaults to np.nan for float
        and 0 for integer dtypes.
    backend : {"numpy"}, default "numpy"
        The backend to use for sparse storage. Only NumPy is supported, as PyArrow
        does not have a native sparse type.

    Returns
    -------
    np.dtype
        A sparse dtype with the specified configuration.

    Examples
    --------
    >>> print(sparse())  # Default numpy backend
    Sparse[float64, nan]
    >>> print(sparse(dtype="int64"))  # With specific dtype
    Sparse[int64, 0]
    >>> print(sparse(fill_value=-1))  # With specific fill value
    Sparse[float64, -1]
    """
    if backend != "numpy":
        raise ValueError(
            "Sparse types are only supported with NumPy backend, as PyArrow "
            "does not have a native sparse type."
        )

    if dtype is None:
        dtype = np.float64
    if fill_value is None:
        fill_value = np.nan if np.issubdtype(dtype, np.floating) else 0
    return SparseDtype(dtype=dtype, fill_value=fill_value)


def date(
    unit: Literal["day", "ms"] = "day",
    backend: Literal["pyarrow"] = "pyarrow",
) -> ArrowDtype:
    """
    Create a date dtype with specified unit, using PyArrow backend.

    This function creates a dtype for representing dates without time components,
    suitable for calendar-based data. PyArrow provides two date types: `date32` for
    day precision (stored as days since UNIX epoch) and `date64` for millisecond
    precision (stored as milliseconds since UNIX epoch). NumPy does not natively
    support a date-only type, so only PyArrow backend is supported.

    Parameters
    ----------
    unit : {"day", "ms"}, default "day"
        The precision unit for the date:
        - "day": Uses `date32`, representing dates as days since UNIX epoch
        (1970-01-01).
        - "ms": Uses `date64`, representing dates as milliseconds since UNIX epoch.
    backend : {"pyarrow"}, default "pyarrow"
        The backend to use for date storage. Only PyArrow is supported.

    Returns
    -------
    ArrowDtype
        A date dtype with the specified configuration, wrapped in an ArrowDtype.

    Raises
    ------
    ValueError
        If a backend other than "pyarrow" is specified.

    Examples
    --------
    >>> print(date())  # Default day precision with PyArrow
    date32[day][pyarrow]
    >>> print(date(unit="ms"))  # Millisecond precision with PyArrow
    date64[ms][pyarrow]
    >>> import pandas as pd
    >>> pd.Series(
    ...     [pd.Timestamp("2023-01-01"), pd.Timestamp("2023-01-02")], dtype=date()
    ... )
    0    2023-01-01
    1    2023-01-02
    dtype: date32[day][pyarrow]
    """

    if backend != "pyarrow":
        raise ValueError("Date types are only supported with PyArrow backend.")
    import pyarrow as pa

    return ArrowDtype(pa.date32() if unit == "day" else pa.date64())


def duration(
    unit: Literal["ns", "us", "ms", "s"] = "ns",
    backend: Literal["numpy", "pyarrow"] = "pyarrow",
) -> np.dtype | ArrowDtype:
    """
    Create a duration dtype with specified unit and backend.

    Parameters
    ----------
    unit : {"ns", "us", "ms", "s"}, default "ns"
        The unit of precision for the duration:
        - "ns": Nanoseconds.
        - "us": Microseconds.
        - "ms": Milliseconds.
        - "s": Seconds.
    backend : {"numpy", "pyarrow"}, default "pyarrow"
        The backend to use for duration storage.

    Returns
    -------
    Union[np.dtype, ArrowDtype]
        A duration dtype with the specified configuration.

    Examples
    --------
    >>> print(duration())  # Default PyArrow backend
    duration[ns][pyarrow]
    >>> print(duration(unit="s", backend="numpy"))  # NumPy backend
    timedelta64[s]
    """
    valid_units = ["ns", "us", "ms", "s"]
    if unit not in valid_units:
        raise ValueError(f"Unit must be one of {valid_units}")
    if backend == "numpy":
        return np.dtype(f"timedelta64[{unit}]")
    else:  # pyarrow
        import pyarrow as pa

        return ArrowDtype(pa.duration(unit))


def map(
    index_type: Any,
    value_type: Any,
    backend: Literal["pyarrow"] = "pyarrow",
) -> ArrowDtype:
    """
    Create a map dtype with specified index and value types, using PyArrow backend.

    This function creates a dtype for representing key-value mappings (similar to a
    dictionary), where each element is a list of key-value pairs. PyArrow's `map` type
    ensures that keys are unique within each element.
    This type is not natively supported by NumPy, so only PyArrow backend is supported.

    Parameters
    ----------
    index_type : Any
        The type of the map's keys (e.g., `pa.int32()`, `pa.string()`).
    value_type : Any
        The type of the map's values (e.g., `pa.float64()`, `pa.string()`).
    backend : {"pyarrow"}, default "pyarrow"
        The backend to use for map storage. Only PyArrow is supported.

    Returns
    -------
    ArrowDtype
        A map dtype with the specified configuration, wrapped in an ArrowDtype.

    Raises
    ------
    ValueError
        If a backend other than "pyarrow" is specified.

    Examples
    --------
    >>> import pyarrow as pa
    >>> print(map(index_type=pa.int32(), value_type=pa.string()))
    map<int32, string>[pyarrow]
    >>> import pandas as pd
    >>> data = [[(1, "a"), (2, "b")], [(3, "c")]]
    >>> pd.Series(data, dtype=map(pa.int32(), pa.string()))
    0    [(1, 'a'), (2, 'b')]
    1              [(3, 'c')]
    dtype: map<int32, string>[pyarrow]
    """
    if backend != "pyarrow":
        raise ValueError("Map types are only supported with PyArrow backend.")
    import pyarrow as pa

    return ArrowDtype(pa.map_(index_type, value_type))


def struct(
    fields: list[tuple[str, Any]],
    backend: Literal["pyarrow"] = "pyarrow",
) -> ArrowDtype:
    """
    Create a struct dtype with specified fields, using PyArrow backend.

    This function creates a dtype for representing structured data, where each element
    is a record with named fields, similar to a named tuple or dictionary. Each field
    in the struct has a name and a type, defined in the `fields` parameter. PyArrow's
    `struct` type is used to store this data, allowing for nested structures. NumPy does
    not natively support structured types in the same way, so only PyArrow backend is
    supported.

    Parameters
    ----------
    fields : list of tuple[str, Any]
        A list of (name, type) tuples defining the fields of the struct, where:
        - `name` is a string representing the field name.
        - `type` is a PyArrow type (e.g., `pa.int32()`, `pa.string()`).
    backend : {"pyarrow"}, default "pyarrow"
        The backend to use for struct storage. Only PyArrow is supported.

    Returns
    -------
    ArrowDtype
        A struct dtype with the specified configuration, wrapped in an ArrowDtype.

    Raises
    ------
    ValueError
        If a backend other than "pyarrow" is specified.

    Examples
    --------
    >>> import pyarrow as pa
    >>> print(struct([("id", pa.int32()), ("name", pa.string())]))
    struct<id: int32, name: string>[pyarrow]
    >>> import pandas as pd
    >>> data = [(1, "Alice"), (2, "Bob")]
    >>> pd.Series(data, dtype=struct([("id", pa.int32()), ("name", pa.string())]))
    0    {'id': 1, 'name': 'Alice'}
    1      {'id': 2, 'name': 'Bob'}
    dtype: struct<id: int32, name: string>[pyarrow]
    """
    if backend == "pyarrow":
        import pyarrow as pa

        pa_fields = [(name, getattr(typ, "pyarrow_dtype", typ)) for name, typ in fields]
        return ArrowDtype(pa.struct(pa_fields))
    else:
        raise ValueError(f"Unsupported backend: {backend!r}")
