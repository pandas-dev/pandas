import pytest
import math
import ctypes


@pytest.mark.parametrize(
    "test_data",
    [
        {"a": ["foo", "bar"], "b": ["baz", "qux"]},
        {"a": [1.5, 2.5, 3.5], "b": [9.2, 10.5, 11.8]},
        {"A": [1, 2, 3, 4], "B": [1, 2, 3, 4]},
    ],
    ids=["str_data", "float_data", "int_data"],
)
def test_only_one_dtype(test_data, df_from_dict):
    columns = list(test_data.keys())
    df = df_from_dict(test_data)
    dfX = df.__dataframe__()

    column_size = len(test_data[columns[0]])
    for column in columns:
        assert dfX.get_column_by_name(column).null_count == 0
        assert dfX.get_column_by_name(column).size == column_size
        assert dfX.get_column_by_name(column).offset == 0


def test_float_int(df_from_dict):
    df = df_from_dict(
        {
            "a": [1, 2, 3],
            "b": [3, 4, 5],
            "c": [1.5, 2.5, 3.5],
            "d": [9, 10, 11],
            "e": [True, False, True],
            "f": ["a", "", "c"],
        }
    )
    dfX = df.__dataframe__()
    columns = {"a": 0, "b": 0, "c": 2, "d": 0, "e": 20, "f": 21}

    for column, kind in columns.items():
        colX = dfX.get_column_by_name(column)
        assert colX.null_count == 0
        assert colX.size == 3
        assert colX.offset == 0

        assert colX.dtype[0] == kind

    assert dfX.get_column_by_name("c").dtype[1] == 64


def test_na_float(df_from_dict):
    df = df_from_dict({"a": [1.0, math.nan, 2.0]})
    dfX = df.__dataframe__()
    colX = dfX.get_column_by_name("a")
    assert colX.null_count == 1


def test_noncategorical(df_from_dict):
    df = df_from_dict({"a": [1, 2, 3]})
    dfX = df.__dataframe__()
    colX = dfX.get_column_by_name("a")
    with pytest.raises(TypeError):
        colX.describe_categorical


def test_categorical(df_from_dict):
    df = df_from_dict(
        {"weekday": ["Mon", "Tue", "Mon", "Wed", "Mon", "Thu", "Fri", "Sat", "Sun"]},
        is_categorical=True,
    )

    colX = df.__dataframe__().get_column_by_name("weekday")
    is_ordered, is_dictionary, _ = colX.describe_categorical
    assert isinstance(is_ordered, bool)
    assert isinstance(is_dictionary, bool)


def test_dataframe(df_from_dict):
    df = df_from_dict(
        {"x": [True, True, False], "y": [1, 2, 0], "z": [9.2, 10.5, 11.8]}
    )
    dfX = df.__dataframe__()

    assert dfX.num_columns() == 3
    assert dfX.num_rows() == 3
    assert dfX.num_chunks() == 1
    assert dfX.column_names() == ["x", "y", "z"]
    assert (
        dfX.select_columns((0, 2)).column_names()
        == dfX.select_columns_by_name(("x", "z")).column_names()
    )


@pytest.mark.parametrize(["size", "n_chunks"], [(10, 3), (12, 3), (12, 5)])
def test_df_get_chunks(size, n_chunks, df_from_dict):
    df = df_from_dict({"x": list(range(size))})
    dfX = df.__dataframe__()
    chunks = list(dfX.get_chunks(n_chunks))
    assert len(chunks) == n_chunks
    assert sum(chunk.num_rows() for chunk in chunks) == size


@pytest.mark.parametrize(["size", "n_chunks"], [(10, 3), (12, 3), (12, 5)])
def test_column_get_chunks(size, n_chunks, df_from_dict):
    df = df_from_dict({"x": list(range(size))})
    dfX = df.__dataframe__()
    chunks = list(dfX.get_column(0).get_chunks(n_chunks))
    assert len(chunks) == n_chunks
    assert sum(chunk.size for chunk in chunks) == size


def test_get_columns(df_from_dict):
    df = df_from_dict({"a": [0, 1], "b": [2.5, 3.5]})
    dfX = df.__dataframe__()
    for colX in dfX.get_columns():
        assert colX.size == 2
        assert colX.num_chunks() == 1
    assert dfX.get_column(0).dtype[0] == 0
    assert dfX.get_column(1).dtype[0] == 2


def test_buffer(df_from_dict):
    arr = [0, 1, -1]
    df = df_from_dict({"a": arr})
    dfX = df.__dataframe__()
    colX = dfX.get_column(0)
    bufX = colX.get_buffers()

    dataBuf, dataDtype = bufX["data"]

    assert dataBuf.bufsize > 0
    assert dataBuf.ptr != 0
    device, _ = dataBuf.__dlpack_device__

    assert dataDtype[0] == 0

    if device == 1:  # CPU-only as we're going to directly read memory here
        bitwidth = dataDtype[1]
        ctype = {
            8: ctypes.c_int8,
            16: ctypes.c_int16,
            32: ctypes.c_int32,
            64: ctypes.c_int64,
        }[bitwidth]

        for idx, truth in enumerate(arr):
            val = ctype.from_address(dataBuf.ptr + idx * (bitwidth // 8)).value
            assert val == truth, f"Buffer at index {idx} mismatch"
