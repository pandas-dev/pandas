from __future__ import annotations

from typing import (
    TYPE_CHECKING,
    Literal,
)

import numpy as np

from pandas._config import using_string_dtype

from pandas._libs import lib
from pandas.compat import (
    pa_version_under18p0,
    pa_version_under19p0,
)
from pandas.compat._optional import import_optional_dependency

from pandas.core.dtypes.common import pandas_dtype

import pandas as pd

if TYPE_CHECKING:
    from collections.abc import (
        Callable,
        Hashable,
        Sequence,
    )

    import pyarrow

    from pandas._typing import (
        DtypeArg,
        DtypeBackend,
    )


def _arrow_dtype_mapping() -> dict:
    pa = import_optional_dependency("pyarrow")
    return {
        pa.int8(): pd.Int8Dtype(),
        pa.int16(): pd.Int16Dtype(),
        pa.int32(): pd.Int32Dtype(),
        pa.int64(): pd.Int64Dtype(),
        pa.uint8(): pd.UInt8Dtype(),
        pa.uint16(): pd.UInt16Dtype(),
        pa.uint32(): pd.UInt32Dtype(),
        pa.uint64(): pd.UInt64Dtype(),
        pa.bool_(): pd.BooleanDtype(),
        pa.string(): pd.StringDtype(),
        pa.float32(): pd.Float32Dtype(),
        pa.float64(): pd.Float64Dtype(),
        pa.string(): pd.StringDtype(),
        pa.large_string(): pd.StringDtype(),
    }


def _arrow_string_types_mapper() -> Callable:
    pa = import_optional_dependency("pyarrow")

    mapping = {
        pa.string(): pd.StringDtype(na_value=np.nan),
        pa.large_string(): pd.StringDtype(na_value=np.nan),
    }
    if not pa_version_under18p0:
        mapping[pa.string_view()] = pd.StringDtype(na_value=np.nan)

    return mapping.get


def arrow_table_to_pandas(
    table: pyarrow.Table,
    dtype_backend: DtypeBackend | Literal["numpy"] | lib.NoDefault = lib.no_default,
    null_to_int64: bool = False,
    to_pandas_kwargs: dict | None = None,
    dtype: DtypeArg | None = None,
    names: Sequence[Hashable] | None = None,
) -> pd.DataFrame:
    pa = import_optional_dependency("pyarrow")

    to_pandas_kwargs = {} if to_pandas_kwargs is None else to_pandas_kwargs

    types_mapper: type[pd.ArrowDtype] | None | Callable
    if dtype_backend == "numpy_nullable":
        mapping = _arrow_dtype_mapping()
        if null_to_int64:
            # Modify the default mapping to also map null to Int64
            # (to match other engines - only for CSV parser)
            mapping[pa.null()] = pd.Int64Dtype()
        types_mapper = mapping.get
    elif dtype_backend == "pyarrow":
        types_mapper = pd.ArrowDtype
    elif using_string_dtype():
        if pa_version_under19p0:
            types_mapper = _arrow_string_types_mapper()
        elif dtype is not None:
            # GH#56136 Avoid lossy conversion to float64
            # We'll convert to numpy below if
            types_mapper = {
                pa.int8(): pd.Int8Dtype(),
                pa.int16(): pd.Int16Dtype(),
                pa.int32(): pd.Int32Dtype(),
                pa.int64(): pd.Int64Dtype(),
            }.get
        else:
            types_mapper = None
    elif dtype_backend is lib.no_default or dtype_backend == "numpy":
        if dtype is not None:
            # GH#56136 Avoid lossy conversion to float64
            # We'll convert to numpy below if
            types_mapper = {
                pa.int8(): pd.Int8Dtype(),
                pa.int16(): pd.Int16Dtype(),
                pa.int32(): pd.Int32Dtype(),
                pa.int64(): pd.Int64Dtype(),
            }.get
        else:
            types_mapper = None
    else:
        raise NotImplementedError

    df = table.to_pandas(types_mapper=types_mapper, **to_pandas_kwargs)
    return _post_convert_dtypes(df, dtype_backend, dtype, names)


def _post_convert_dtypes(
    df: pd.DataFrame,
    dtype_backend: DtypeBackend | Literal["numpy"] | lib.NoDefault,
    dtype: DtypeArg | None,
    names: Sequence[Hashable] | None,
) -> pd.DataFrame:
    if dtype is not None and (
        dtype_backend is lib.no_default or dtype_backend == "numpy"
    ):
        # GH#56136 apply any user-provided dtype, and convert any IntegerDtype
        #  columns the user didn't explicitly ask for.
        if isinstance(dtype, dict):
            if names is not None:
                df.columns = names

            cmp_dtypes = {
                pd.Int8Dtype(),
                pd.Int16Dtype(),
                pd.Int32Dtype(),
                pd.Int64Dtype(),
            }
            for col in df.columns:
                if col not in dtype and df[col].dtype in cmp_dtypes:
                    # Any key that the user didn't explicitly specify
                    #  that got converted to IntegerDtype now gets converted
                    #  to numpy dtype.
                    dtype[col] = df[col].dtype.numpy_dtype

            # Ignore non-existent columns from dtype mapping
            # like other parsers do
            dtype = {
                key: pandas_dtype(dtype[key]) for key in dtype if key in df.columns
            }

        else:
            dtype = pandas_dtype(dtype)

        try:
            df = df.astype(dtype)
        except TypeError as err:
            # GH#44901 reraise to keep api consistent
            raise ValueError(str(err)) from err

    if (
        not using_string_dtype()
        and dtype != "str"
        and (dtype_backend is lib.no_default or dtype_backend == "numpy")
    ):
        # Convert any StringDtype columns back to object dtype (pyarrow always
        # uses string dtype even when the infer_string option is False)
        for col, dtype in zip(df.columns, df.dtypes, strict=True):
            if isinstance(dtype, pd.StringDtype) and dtype.na_value is np.nan:
                df[col] = df[col].astype("object").fillna(None)
            if isinstance(dtype, pd.CategoricalDtype):
                cat_dtype = dtype.categories.dtype
                if (
                    isinstance(cat_dtype, pd.StringDtype)
                    and cat_dtype.na_value is np.nan
                ):
                    cat_dtype = pd.CategoricalDtype(
                        categories=dtype.categories.astype("object"),
                        ordered=dtype.ordered,
                    )
                    df[col] = df[col].astype(cat_dtype)

    return df
