import random

import numpy as np
import pytest

import pandas as pd
from pandas.core.exchange.dataframe_protocol import (
    ColumnNullType,
    DtypeKind,
)
from pandas.core.exchange.from_dataframe import from_dataframe
from pandas.testing import assert_frame_equal

test_data_categorical = {
    "ordered": pd.Categorical(list("testdata") * 30, ordered=True),
    "unordered": pd.Categorical(list("testdata") * 30, ordered=False),
}

NCOLS, NROWS = 100, 200

int_data = {
    f"col{int((i - NCOLS / 2) % NCOLS + 1)}": [
        random.randint(0, 100) for _ in range(NROWS)
    ]
    for i in range(NCOLS)
}

bool_data = {
    f"col{int((i - NCOLS / 2) % NCOLS + 1)}": [
        random.choice([True, False]) for _ in range(NROWS)
    ]
    for i in range(NCOLS)
}

float_data = {
    f"col{int((i - NCOLS / 2) % NCOLS + 1)}": [
        random.random() for _ in range(NROWS)
    ]
    for i in range(NCOLS)
}

string_data = {
    "separator data": [
        "abC|DeF,Hik",
        "234,3245.67",
        "gSaf,qWer|Gre",
        "asd3,4sad|",
        np.NaN,
    ]
}


@pytest.mark.parametrize("data", [("ordered", True), ("unordered", False)])
def test_categorical_dtype(data):
    df = pd.DataFrame({"A": (test_data_categorical[data[0]])})

    col = df.__dataframe__().get_column_by_name("A")
    assert col.dtype[0] == DtypeKind.CATEGORICAL
    assert col.null_count == 0
    assert col.describe_null == (ColumnNullType.USE_SENTINEL, -1)
    assert col.num_chunks() == 1
    assert col.describe_categorical == {
        "is_ordered": data[1],
        "is_dictionary": True,
        "mapping": {0: "a", 1: "d", 2: "e", 3: "s", 4: "t"},
    }

    assert_frame_equal(df, from_dataframe(df.__dataframe__()))


@pytest.mark.parametrize("data", [int_data, float_data, bool_data])
def test_dataframe(data):
    df = pd.DataFrame(data)

    df2 = df.__dataframe__()

    assert df2._allow_copy is True
    assert df2.num_columns() == NCOLS
    assert df2.num_rows() == NROWS

    assert list(df2.column_names()) == list(data.keys())

    indices = (0, 2)
    names = tuple(list(data.keys())[idx] for idx in indices)

    assert_frame_equal(
        from_dataframe(df2.select_columns(indices)),
        from_dataframe(df2.select_columns_by_name(names)),
    )


def test_missing_from_masked():
    df = pd.DataFrame(
        {
            "x": np.array([1, 2, 3, 4, 0]),
            "y": np.array([1.5, 2.5, 3.5, 4.5, 0]),
            "z": np.array([True, False, True, True, True]),
        }
    )

    df2 = df.__dataframe__()

    rng = np.random.RandomState(42)
    dict_null = {col: rng.randint(low=0, high=len(df)) for col in df.columns}
    for col, num_nulls in dict_null.items():
        null_idx = df.index[
            rng.choice(np.arange(len(df)), size=num_nulls, replace=False)
        ]
        df.loc[null_idx, col] = None

    df2 = df.__dataframe__()

    assert df2.get_column_by_name("x").null_count == dict_null["x"]
    assert df2.get_column_by_name("y").null_count == dict_null["y"]
    assert df2.get_column_by_name("z").null_count == dict_null["z"]


@pytest.mark.parametrize(
    "data",
    [
        {"x": [1.5, 2.5, 3.5], "y": [9.2, 10.5, 11.8]},
        {"x": [1, 2, 0], "y": [9.2, 10.5, 11.8]},
        {
            "x": np.array([True, True, False]),
            "y": np.array([1, 2, 0]),
            "z": np.array([9.2, 10.5, 11.8]),
        },
    ],
)
def test_mixed_data(data):
    df = pd.DataFrame(data)
    df2 = df.__dataframe__()

    for col_name in df.columns:
        assert df2.get_column_by_name(col_name).null_count == 0


def test_mixed_missing():
    df = pd.DataFrame(
        {
            "x": np.array([True, None, False, None, True]),
            "y": np.array([None, 2, None, 1, 2]),
            "z": np.array([9.2, 10.5, None, 11.8, None]),
        }
    )

    df2 = df.__dataframe__()

    for col_name in df.columns:
        assert df2.get_column_by_name(col_name).null_count == 2


def test_select_columns_error():
    df = pd.DataFrame(int_data)

    df2 = df.__dataframe__()

    with pytest.raises(ValueError):
        assert from_dataframe(df2.select_columns(np.array([0, 2]))) == from_dataframe(
            df2.select_columns_by_name(("col33", "col35"))
        )


def test_select_columns_by_name_error():
    df = pd.DataFrame(int_data)

    df2 = df.__dataframe__()

    with pytest.raises(ValueError):
        assert from_dataframe(
            df2.select_columns_by_name(np.array(["col33", "col35"]))
        ) == from_dataframe(df2.select_columns((0, 2)))


def test_string():
    test_str_data = string_data["separator data"] + [""]
    df = pd.DataFrame({"A": test_str_data})
    col = df.__dataframe__().get_column_by_name("A")

    assert col.size == 6
    assert col.null_count == 1
    assert col.dtype[0] == DtypeKind.STRING
    assert col.describe_null == (ColumnNullType.USE_BYTEMASK, 0)

    df_sliced = df[1:]
    col = df_sliced.__dataframe__().get_column_by_name("A")
    assert col.size == 5
    assert col.null_count == 1
    assert col.dtype[0] == DtypeKind.STRING
    assert col.describe_null == (ColumnNullType.USE_BYTEMASK, 0)
