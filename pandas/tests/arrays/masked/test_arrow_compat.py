import pytest

import pandas.util._test_decorators as td

import pandas as pd
import pandas._testing as tm

arrays = [pd.array([1, 2, 3, None], dtype=dtype) for dtype in tm.ALL_EA_INT_DTYPES]
arrays += [pd.array([0.1, 0.2, 0.3, None], dtype=dtype) for dtype in tm.FLOAT_EA_DTYPES]
arrays += [pd.array([True, False, True, None], dtype="boolean")]


@pytest.fixture(params=arrays, ids=[a.dtype.name for a in arrays])
def data(request):
    return request.param


@td.skip_if_no("pyarrow", min_version="0.15.0")
def test_arrow_array(data):
    # protocol added in 0.15.0
    import pyarrow as pa

    arr = pa.array(data)
    expected = pa.array(
        data.to_numpy(object, na_value=None),
        type=pa.from_numpy_dtype(data.dtype.numpy_dtype),
    )
    assert arr.equals(expected)


@td.skip_if_no("pyarrow", min_version="0.16.0")
def test_arrow_roundtrip(data):
    # roundtrip possible from arrow 0.16.0
    import pyarrow as pa

    df = pd.DataFrame({"a": data})
    table = pa.table(df)
    assert table.field("a").type == str(data.dtype.numpy_dtype)
    result = table.to_pandas()
    assert result["a"].dtype == data.dtype
    tm.assert_frame_equal(result, df)


@td.skip_if_no("pyarrow", min_version="0.15.1.dev")
def test_arrow_load_from_zero_chunks(data):
    # GH-41040
    import pyarrow as pa

    df = pd.DataFrame({"a": data[0:0]})
    table = pa.table(df)
    assert table.field("a").type == str(data.dtype.numpy_dtype)
    table = pa.table(
        [pa.chunked_array([], type=table.field("a").type)], schema=table.schema
    )
    result = table.to_pandas()
    assert result["a"].dtype == data.dtype
    tm.assert_frame_equal(result, df)


@td.skip_if_no("pyarrow", min_version="0.16.0")
def test_arrow_from_arrow_uint():
    # https://github.com/pandas-dev/pandas/issues/31896
    # possible mismatch in types
    import pyarrow as pa

    dtype = pd.UInt32Dtype()
    result = dtype.__from_arrow__(pa.array([1, 2, 3, 4, None], type="int64"))
    expected = pd.array([1, 2, 3, 4, None], dtype="UInt32")

    tm.assert_extension_array_equal(result, expected)


@td.skip_if_no("pyarrow", min_version="0.16.0")
def test_arrow_sliced(data):
    # https://github.com/pandas-dev/pandas/issues/38525
    import pyarrow as pa

    df = pd.DataFrame({"a": data})
    table = pa.table(df)
    result = table.slice(2, None).to_pandas()
    expected = df.iloc[2:].reset_index(drop=True)
    tm.assert_frame_equal(result, expected)

    # no missing values
    df2 = df.fillna(data[0])
    table = pa.table(df2)
    result = table.slice(2, None).to_pandas()
    expected = df2.iloc[2:].reset_index(drop=True)
    tm.assert_frame_equal(result, expected)


@td.skip_if_no("pyarrow", min_version="0.16.0")
def test_from_arrow_type_error(request, data):
    # ensure that __from_arrow__ returns a TypeError when getting a wrong
    # array type
    import pyarrow as pa

    if data.dtype != "boolean":
        # TODO numeric dtypes cast any incoming array to the correct dtype
        # instead of erroring
        request.node.add_marker(
            pytest.mark.xfail(reason="numeric dtypes don't error but cast")
        )

    arr = pa.array(data).cast("string")
    with pytest.raises(TypeError, match=None):
        # we don't test the exact error message, only the fact that it raises
        # a TypeError is relevant
        data.dtype.__from_arrow__(arr)
