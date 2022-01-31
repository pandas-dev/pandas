import pytest

from genomictypes import Variant, VariantArray
import numpy as np
import pandas as pd

from pandas.core.arrays import StructArray
from pandas.core.arrays import StructDtype

# from pandas.core.dtypes.dtypes import StructDtype


def test_variantdtype():
    var_dtype = StructDtype({
        "chrom": "string",
        "start": "int32",
        "end": "int32",
        "ref": "string",
        "alt": "string"
    }, nullable=True)
    print(var_dtype)

    var_scalar_type = var_dtype.type

    var1 = var_scalar_type("chr1", 10, 11, "A", "G")
    var2 = var_scalar_type("chr1", 12, 15, "AAT", "C")

    assert tuple(var1) == var1
    assert var2 == (
        var2.chrom,
        var2.start,
        var2.end,
        var2.ref,
        var2.alt
    )

    # test array
    var_array = StructArray([
        var1,
        ("chr1", 10, 11, "A", "G"),
        var2,
        None
    ], dtype=var_dtype)
    list(var_array)
    print(var_array)

    assert len(var_array) == 4
    assert all(var_array[:2] == StructArray([var1, var1]))

    scalar = var_array[0]
    assert isinstance(scalar, var_scalar_type)
    assert scalar == var1

    assert all(var_array.isna() == [False, False, False, True])
    assert all(
        (np.asarray(var_array) == np.array([var1, var1, var2, None], dtype=object)) | pd.isna(var_array)
    )

    # TODO: proper test values
    assert all(
        StructArray.from_df(
            var_array.as_frame()
        ) == var_array
    )

    assert len(var_array.unique()) == 3
    assert set(var_array.unique()) == {var1, var2, pd.NA}

    # don't explicitly check the repr, only test if it works
    str(var_array)

    assert all(pd.array(var_array) == var_array)

    # # test pyarrow
    # arrow_array = var_array.__arrow_array__()
    # arrow_array.to_pandas()


def test_structarray():
    dtype = StructDtype({"x": "int", "y": "float"}, nullable=True)
    x = StructArray([(0, 3), (1, 2)], dtype=dtype)

    assert type(x[0]) == dtype.type
    assert x[0] == (0, 3)

    assert all(StructArray.from_df(x.as_frame()) == x)

    # x == y
    assert all((x == (0, 1)) == [False, False])
    assert all((x == (0, 3)) == [True, False])
    # x != y
    assert all((x != (0, 1)) == [True, True])
    assert all((x != (0, 3)) == [False, True])
    # x < y
    assert all((x < (0, 3)) == [False, False])
    assert all((x < (0, 5)) == [True, False])
    assert all((x < (3, 5)) == [True, True])
    # x <= y
    assert all((x <= (0, 3)) == [True, False])
    assert all((x <= (0, 5)) == [True, False])
    assert all((x <= (3, 5)) == [True, True])
    # x >= y
    assert all((x >= (0, 3)) == [True, True])
    assert all((x >= (0, 5)) == [False, True])
    assert all((x >= (3, 5)) == [False, False])
    # x > y
    assert all((x > (0, 3)) == [False, True])
    assert all((x > (0, 5)) == [False, True])
    assert all((x > (3, 5)) == [False, False])

    # test type inference
    assert pd.api.types.is_dtype_equal(
        StructArray([(0, 3), (1, 2.), None]).dtype,
        StructDtype({"f_0": "int64", "f_1": "float64"}, nullable=True)
    )
    assert not pd.api.types.is_dtype_equal(
        StructArray([(0, 3), (1, 2.), None]).dtype,
        StructDtype({"f_0": "int64", "f_1": "float64"}, nullable=False)
    )
