import pytest
import numpy as np

@pytest.mark.parametrize("test_data",
                         [
                             ({'a': [np.array([1, 2, 3]), np.array([4, 5, 6])],
                               'b': [np.array([1.5, 2.0, 3.2]), np.array([4.1, 5.7, 6.9])]},
                              np.object_, None),
                             ({'a': [1.5, 2.5, 3.5], 'b': [9.2, 10.5, 11.8]}, np.float64, None),
                             ({'A': [1, 2, 3, 4], 'B': [1, 2, 3, 4]}, np.int64, np.float64)
                         ],
                         ids=["array_data", "float_data", "int_data"])
def test_only_one_data(test_data, create_df_from_dict):
    data, dtype, new_dtype = test_data
    columns = list(data.keys())
    df = create_df_from_dict(data)
    df2 = df.__dataframe__()
    new_dtype = dtype if new_dtype is None else new_dtype
    assert df.columns.values.tolist() == columns
    val = len(df[columns[0]])-1
    column_size = df.size
    for column in columns:
        assert df[column].tolist() == df[column].tolist()
        assert df[column].dtype.type is dtype
        assert df2.get_column_by_name(column).null_count == 0
        assert df2.get_column_by_name(column).size == column_size
        assert df2.get_column_by_name(column).offset == 0
        assert not df2["x"].is_masked
        n = np.random.randint(0, val)
        (df[column])[n] = None
        assert df[column].dtype.type is new_dtype
        assert df2.get_column_by_name(column).null_count == 1


def test_float_int(create_df_from_dict):
    df = create_df_from_dict({'a': [1, 2, 3], 'b': [3, 4, 5],
                       'c': [1.5, 2.5, 3.5], 'd': [9, 10, 11]})
    df2 = df.__dataframe__()
    columns = ['a', 'b', 'c', 'd']
    assert df.columns.values.tolist() == columns
    for column in columns:
        assert df[column].tolist() == df[column].tolist()
        if column is 'c':
            assert df[column].dtype.type is np.float64
        else:
            assert df[column].dtype.type is np.int64

        assert df2.get_column_by_name(column).null_count == 0
        assert df2.get_column_by_name(column).size == 3
        assert df2.get_column_by_name(column).offset == 0

        n = np.random.randint(0, 2)
        (df[column])[n] = None
        assert df[column].dtype.type is np.float64
        assert df2.get_column_by_name(column).null_count == 1


def test_mixed_intfloatbool(create_df_from_dict):
    df = create_df_from_dict({"x": np.array([True, True, False]),
                       "y": np.array([1, 2, 0]),
                       "z": np.array([9.2, 10.5, 11.8])})
    df2 = df.__dataframe__()
    columns = ['x', 'y', 'z']
    assert df.columns.values.tolist() == columns
    for column in columns:
        assert df[column].tolist() == df[column].tolist()
        assert df2.get_column_by_name(column).null_count == 0
        assert df2.get_column_by_name(column).size == 3
        assert df2.get_column_by_name(column).offset == 0

    assert df["x"].dtype.type is np.bool_
    assert df["y"].dtype.type is np.int32
    assert df["z"].dtype.type is np.float64

    assert df2.get_column_by_name("x")._allow_copy == True

    for column in columns:
        n = np.random.randint(0, 2)
        (df[column])[n] = None
        if column is "x":
            assert df[column].dtype.type is np.object_
        else:
            assert df[column].dtype.type is np.float64
        assert df2.get_column_by_name(column).null_count == 1


def test_string_dtype(create_df_from_dict):
    df = create_df_from_dict({"A": ["a", "b", "cdef", "", "g"]})
    df2 = df.__dataframe__()
    columns = ['A']
    assert df.columns.values.tolist() == columns
    for column in columns:
        assert df[column].tolist() == df[column].tolist()
        assert df[column].dtype.type is np.object_
        assert df2.get_column_by_name(column).null_count == 0


def test_categorical(create_df_from_dict):
    df = create_df_from_dict({"year": [2012, 2013, 2015, 2019], "weekday": [0, 1, 4, 6]})
    df = df.categorize("year", min_value=2012, max_value=2019)
    df = df.categorize("weekday", labels=["Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun"])
    # Some detailed testing for correctness of dtype and null handling:
    col = df.__dataframe__().get_column_by_name("year")
    assert col.describe_categorical == (False, True, {0: 2012, 1: 2013, 2: 2014, 3: 2015, 4: 2016, 5: 2017, 6: 2018, 7: 2019})
    assert col.describe_null == (0, None)
    col2 = df.__dataframe__().get_column_by_name("weekday")
    assert col2.describe_categorical == (False, True, {0: "Mon", 1: "Tue", 2: "Wed", 3: "Thu", 4: "Fri", 5: "Sat", 6: "Sun"})
    assert col2.describe_null == (0, None)


def test_dataframe(create_df_from_dict):
    df = create_df_from_dict({"x": [True, True, False], "y": [1, 2, 0], "z": [9.2, 10.5, 11.8]})
    df2 = df.__dataframe__()
    assert df2._allow_copy == True
    assert df2.num_columns() == 3
    assert df2.num_rows() == 3
    assert df2.num_chunks() == 1
    assert df2.column_names() == ["x", "y", "z"]
    assert df2.select_columns((0, 2))._df[:, 0].tolist() == df2.select_columns_by_name(("x", "z"))._df[:, 0].tolist()
    assert df2.select_columns((0, 2))._df[:, 1].tolist() == df2.select_columns_by_name(("x", "z"))._df[:, 1].tolist()


def test_chunks(create_df_from_dict):
    df = create_df_from_dict({"x": [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]})
    df2 = df.__dataframe__()
    chunk_iter = iter(df2.get_chunks(3))
    chunk = next(chunk_iter)
    assert chunk.num_rows() == 4
    chunk = next(chunk_iter)
    assert chunk.num_rows() == 4
    chunk = next(chunk_iter)
    assert chunk.num_rows() == 2
    with pytest.raises(StopIteration):
        chunk = next(chunk_iter)


def test_get_chunks(create_df_from_dict):
    df = create_df_from_dict({"x": [1]})
    df2 = df.__dataframe__()
    assert df2.get_chunks() == 1
