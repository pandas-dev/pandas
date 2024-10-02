from __future__ import annotations

from typing import (
    TYPE_CHECKING,
    Literal,
)

import numpy as np

from pandas._libs import (
    lib,
    missing as libmissing,
)
from pandas.util._validators import check_dtype_backend

from pandas.core.dtypes.cast import maybe_downcast_numeric
from pandas.core.dtypes.common import (
    ensure_object,
    is_bool_dtype,
    is_decimal,
    is_integer_dtype,
    is_number,
    is_numeric_dtype,
    is_scalar,
    is_string_dtype,
    needs_i8_conversion,
)
from pandas.core.dtypes.dtypes import ArrowDtype
from pandas.core.dtypes.generic import (
    ABCIndex,
    ABCSeries,
)

from pandas.core.arrays import BaseMaskedArray
from pandas.core.arrays.string_ import StringDtype

if TYPE_CHECKING:
    from pandas._typing import (
        DateTimeErrorChoices,
        DtypeBackend,
        npt,
    )


def parse_numeric(value):
    if isinstance(value, str):
        try:
            return int(value, 0)  # Automatically detect radix
        except ValueError:
            try:
                return float(value)
            except ValueError:
                return libmissing.NA
    return value


def to_numeric(
        arg,
        errors: DateTimeErrorChoices = "raise",
        downcast: Literal["integer", "signed", "unsigned", "float"] | None = None,
        dtype_backend: DtypeBackend | lib.NoDefault = lib.no_default,
):
    """
    Convert argument to a numeric type.
    ...
    """
    if downcast not in (None, "integer", "signed", "unsigned", "float"):
        raise ValueError("invalid downcasting method provided")

    if errors not in ("raise", "coerce"):
        raise ValueError("invalid error value specified")

    check_dtype_backend(dtype_backend)

    is_series = False
    is_index = False
    is_scalars = False

    if isinstance(arg, ABCSeries):
        is_series = True
        values = arg.values
    elif isinstance(arg, ABCIndex):
        is_index = True
        if needs_i8_conversion(arg.dtype):
            values = arg.view("i8")
        else:
            values = arg.values
    elif isinstance(arg, (list, tuple)):
        values = np.array(arg, dtype="O")
    elif is_scalar(arg):
        if is_decimal(arg):
            return float(arg)
        if is_number(arg):
            return arg
        is_scalars = True
        values = np.array([arg], dtype="O")
    elif getattr(arg, "ndim", 1) > 1:
        raise TypeError("arg must be a list, tuple, 1-d array, or Series")
    else:
        values = arg

    mask: npt.NDArray[np.bool_] | None = None
    if isinstance(values, BaseMaskedArray):
        mask = values._mask
        values = values._data[~mask]

    values_dtype = getattr(values, "dtype", None)
    if isinstance(values_dtype, ArrowDtype):
        mask = values.isna()
        values = values.dropna().to_numpy()
    new_mask: np.ndarray | None = None

    if is_numeric_dtype(values_dtype):
        pass
    elif lib.is_np_dtype(values_dtype, "mM"):
        values = values.view(np.int64)
    else:
        values = ensure_object(values)
        parsed_values = []
        new_mask = []
        for idx, x in enumerate(values):
            parsed_value = parse_numeric(x)
            if libmissing.checknull(parsed_values):
                if errors == "raise":
                    raise ValueError(f"Unable to parse string '{x}' at position {idx}")
                elif errors == "coerce":
                    parsed_values.append(libmissing.NA)
                    new_mask.append(True)
                    continue
            else:
                parsed_values.append(parsed_value)
                new_mask.append(False)

        values = np.array(parsed_values, dtype=object)
        new_mask = np.array(new_mask, dtype=bool)

    if new_mask is not None:
        values = values[~new_mask]
    elif (
            dtype_backend is not lib.no_default
            and new_mask is None
            or isinstance(values_dtype, StringDtype)
            and values_dtype.na_value is libmissing.NA
    ):
        new_mask = np.zeros(values.shape, dtype=np.bool_)

    if downcast is not None and is_numeric_dtype(values.dtype):
        typecodes: str | None = None

        if downcast in ("integer", "signed"):
            typecodes = np.typecodes["Integer"]
        elif downcast == "unsigned" and (not len(values) or np.min(values) >= 0):
            typecodes = np.typecodes["UnsignedInteger"]
        elif downcast == "float":
            typecodes = np.typecodes["Float"]
            float_32_char = np.dtype(np.float32).char
            float_32_ind = typecodes.index(float_32_char)
            typecodes = typecodes[float_32_ind:]

        if typecodes is not None:
            for typecode in typecodes:
                dtype = np.dtype(typecode)
                if dtype.itemsize <= values.dtype.itemsize:
                    # Only downcast if values are all integers
                    if downcast in ("integer", "signed", "unsigned") and not np.isin(np.mod(values, 1), 0).all():
                        continue  # Skip downcasting if there are any float values
                    values = maybe_downcast_numeric(values, dtype)
                    if values.dtype == dtype:
                        break

    if (mask is not None or new_mask is not None) and not is_string_dtype(values.dtype):
        if mask is None or (new_mask is not None and new_mask.shape == mask.shape):
            mask = new_mask
        else:
            mask = mask.copy()
        assert isinstance(mask, np.ndarray)
        data = np.zeros(mask.shape, dtype=values.dtype)
        data[~mask] = values

        from pandas.core.arrays import (
            ArrowExtensionArray,
            BooleanArray,
            FloatingArray,
            IntegerArray,
        )

        klass: type[IntegerArray | BooleanArray | FloatingArray]
        if is_integer_dtype(data.dtype):
            klass = IntegerArray
        elif is_bool_dtype(data.dtype):
            klass = BooleanArray
        else:
            klass = FloatingArray
        values = klass(data, mask)

        if dtype_backend == "pyarrow" or isinstance(values_dtype, ArrowDtype):
            values = ArrowExtensionArray(values.__arrow_array__())

    if is_series:
        return arg._constructor(values, index=arg.index, name=arg.name)
    elif is_index:
        from pandas import Index
        return Index(values, name=arg.name)
    elif is_scalars:
        return values[0]
    else:
        return values


if __name__ == "__main__":
    import numpy as np

    test_data = ["0x1A", "0b1010", "0o17", "25", "3.14", "invalid"]
    result = to_numeric(test_data, errors="coerce")
    print("Inputs:", test_data)
    print("ParseResult:", result)
