"""
Tests dtype specification during parsing
for all of the parsers defined in parsers.py
"""
from io import StringIO
import os

import numpy as np
import pytest

from pandas.errors import ParserWarning

from pandas.core.dtypes.dtypes import CategoricalDtype

import pandas as pd
from pandas import Categorical, DataFrame, Index, MultiIndex, Series, Timestamp, concat
import pandas._testing as tm


@pytest.mark.parametrize("dtype", [str, object])
@pytest.mark.parametrize("check_orig", [True, False])
def test_dtype_all_columns(all_parsers, dtype, check_orig):
    # see gh-3795, gh-6607
    parser = all_parsers

    df = DataFrame(
        np.random.rand(5, 2).round(4),
        columns=list("AB"),
        index=["1A", "1B", "1C", "1D", "1E"],
    )

    with tm.ensure_clean("__passing_str_as_dtype__.csv") as path:
        df.to_csv(path)

        result = parser.read_csv(path, dtype=dtype, index_col=0)

        if check_orig:
            expected = df.copy()
            result = result.astype(float)
        else:
            expected = df.astype(str)

        tm.assert_frame_equal(result, expected)


def test_dtype_all_columns_empty(all_parsers):
    # see gh-12048
    parser = all_parsers
    result = parser.read_csv(StringIO("A,B"), dtype=str)

    expected = DataFrame({"A": [], "B": []}, index=[], dtype=str)
    tm.assert_frame_equal(result, expected)


def test_dtype_per_column(all_parsers):
    parser = all_parsers
    data = """\
one,two
1,2.5
2,3.5
3,4.5
4,5.5"""
    expected = DataFrame(
        [[1, "2.5"], [2, "3.5"], [3, "4.5"], [4, "5.5"]], columns=["one", "two"]
    )
    expected["one"] = expected["one"].astype(np.float64)
    expected["two"] = expected["two"].astype(object)

    result = parser.read_csv(StringIO(data), dtype={"one": np.float64, 1: str})
    tm.assert_frame_equal(result, expected)


def test_invalid_dtype_per_column(all_parsers):
    parser = all_parsers
    data = """\
one,two
1,2.5
2,3.5
3,4.5
4,5.5"""

    with pytest.raises(TypeError, match="data type [\"']foo[\"'] not understood"):
        parser.read_csv(StringIO(data), dtype={"one": "foo", 1: "int"})


@pytest.mark.parametrize(
    "dtype",
    [
        "category",
        CategoricalDtype(),
        {"a": "category", "b": "category", "c": CategoricalDtype()},
    ],
)
def test_categorical_dtype(all_parsers, dtype):
    # see gh-10153
    parser = all_parsers
    data = """a,b,c
1,a,3.4
1,a,3.4
2,b,4.5"""
    expected = DataFrame(
        {
            "a": Categorical(["1", "1", "2"]),
            "b": Categorical(["a", "a", "b"]),
            "c": Categorical(["3.4", "3.4", "4.5"]),
        }
    )
    actual = parser.read_csv(StringIO(data), dtype=dtype)
    tm.assert_frame_equal(actual, expected)


@pytest.mark.parametrize("dtype", [{"b": "category"}, {1: "category"}])
def test_categorical_dtype_single(all_parsers, dtype):
    # see gh-10153
    parser = all_parsers
    data = """a,b,c
1,a,3.4
1,a,3.4
2,b,4.5"""
    expected = DataFrame(
        {"a": [1, 1, 2], "b": Categorical(["a", "a", "b"]), "c": [3.4, 3.4, 4.5]}
    )
    actual = parser.read_csv(StringIO(data), dtype=dtype)
    tm.assert_frame_equal(actual, expected)


def test_categorical_dtype_unsorted(all_parsers):
    # see gh-10153
    parser = all_parsers
    data = """a,b,c
1,b,3.4
1,b,3.4
2,a,4.5"""
    expected = DataFrame(
        {
            "a": Categorical(["1", "1", "2"]),
            "b": Categorical(["b", "b", "a"]),
            "c": Categorical(["3.4", "3.4", "4.5"]),
        }
    )
    actual = parser.read_csv(StringIO(data), dtype="category")
    tm.assert_frame_equal(actual, expected)


def test_categorical_dtype_missing(all_parsers):
    # see gh-10153
    parser = all_parsers
    data = """a,b,c
1,b,3.4
1,nan,3.4
2,a,4.5"""
    expected = DataFrame(
        {
            "a": Categorical(["1", "1", "2"]),
            "b": Categorical(["b", np.nan, "a"]),
            "c": Categorical(["3.4", "3.4", "4.5"]),
        }
    )
    actual = parser.read_csv(StringIO(data), dtype="category")
    tm.assert_frame_equal(actual, expected)


@pytest.mark.slow
def test_categorical_dtype_high_cardinality_numeric(all_parsers):
    # see gh-18186
    parser = all_parsers
    data = np.sort([str(i) for i in range(524289)])
    expected = DataFrame({"a": Categorical(data, ordered=True)})

    actual = parser.read_csv(StringIO("a\n" + "\n".join(data)), dtype="category")
    actual["a"] = actual["a"].cat.reorder_categories(
        np.sort(actual.a.cat.categories), ordered=True
    )
    tm.assert_frame_equal(actual, expected)


def test_categorical_dtype_latin1(all_parsers, csv_dir_path):
    # see gh-10153
    pth = os.path.join(csv_dir_path, "unicode_series.csv")
    parser = all_parsers
    encoding = "latin-1"

    expected = parser.read_csv(pth, header=None, encoding=encoding)
    expected[1] = Categorical(expected[1])

    actual = parser.read_csv(pth, header=None, encoding=encoding, dtype={1: "category"})
    tm.assert_frame_equal(actual, expected)


def test_categorical_dtype_utf16(all_parsers, csv_dir_path):
    # see gh-10153
    pth = os.path.join(csv_dir_path, "utf16_ex.txt")
    parser = all_parsers
    encoding = "utf-16"
    sep = ","

    expected = parser.read_csv(pth, sep=sep, encoding=encoding)
    expected = expected.apply(Categorical)

    actual = parser.read_csv(pth, sep=sep, encoding=encoding, dtype="category")
    tm.assert_frame_equal(actual, expected)


def test_categorical_dtype_chunksize_infer_categories(all_parsers):
    # see gh-10153
    parser = all_parsers
    data = """a,b
1,a
1,b
1,b
2,c"""
    expecteds = [
        DataFrame({"a": [1, 1], "b": Categorical(["a", "b"])}),
        DataFrame({"a": [1, 2], "b": Categorical(["b", "c"])}, index=[2, 3]),
    ]
    actuals = parser.read_csv(StringIO(data), dtype={"b": "category"}, chunksize=2)

    for actual, expected in zip(actuals, expecteds):
        tm.assert_frame_equal(actual, expected)


def test_categorical_dtype_chunksize_explicit_categories(all_parsers):
    # see gh-10153
    parser = all_parsers
    data = """a,b
1,a
1,b
1,b
2,c"""
    cats = ["a", "b", "c"]
    expecteds = [
        DataFrame({"a": [1, 1], "b": Categorical(["a", "b"], categories=cats)}),
        DataFrame(
            {"a": [1, 2], "b": Categorical(["b", "c"], categories=cats)}, index=[2, 3]
        ),
    ]
    dtype = CategoricalDtype(cats)
    actuals = parser.read_csv(StringIO(data), dtype={"b": dtype}, chunksize=2)

    for actual, expected in zip(actuals, expecteds):
        tm.assert_frame_equal(actual, expected)


@pytest.mark.parametrize("ordered", [False, True])
@pytest.mark.parametrize(
    "categories",
    [["a", "b", "c"], ["a", "c", "b"], ["a", "b", "c", "d"], ["c", "b", "a"]],
)
def test_categorical_category_dtype(all_parsers, categories, ordered):
    parser = all_parsers
    data = """a,b
1,a
1,b
1,b
2,c"""
    expected = DataFrame(
        {
            "a": [1, 1, 1, 2],
            "b": Categorical(
                ["a", "b", "b", "c"], categories=categories, ordered=ordered
            ),
        }
    )

    dtype = {"b": CategoricalDtype(categories=categories, ordered=ordered)}
    result = parser.read_csv(StringIO(data), dtype=dtype)
    tm.assert_frame_equal(result, expected)


def test_categorical_category_dtype_unsorted(all_parsers):
    parser = all_parsers
    data = """a,b
1,a
1,b
1,b
2,c"""
    dtype = CategoricalDtype(["c", "b", "a"])
    expected = DataFrame(
        {
            "a": [1, 1, 1, 2],
            "b": Categorical(["a", "b", "b", "c"], categories=["c", "b", "a"]),
        }
    )

    result = parser.read_csv(StringIO(data), dtype={"b": dtype})
    tm.assert_frame_equal(result, expected)


def test_categorical_coerces_numeric(all_parsers):
    parser = all_parsers
    dtype = {"b": CategoricalDtype([1, 2, 3])}

    data = "b\n1\n1\n2\n3"
    expected = DataFrame({"b": Categorical([1, 1, 2, 3])})

    result = parser.read_csv(StringIO(data), dtype=dtype)
    tm.assert_frame_equal(result, expected)


def test_categorical_coerces_datetime(all_parsers):
    parser = all_parsers
    dtype = {"b": CategoricalDtype(pd.date_range("2017", "2019", freq="AS"))}

    data = "b\n2017-01-01\n2018-01-01\n2019-01-01"
    expected = DataFrame({"b": Categorical(dtype["b"].categories)})

    result = parser.read_csv(StringIO(data), dtype=dtype)
    tm.assert_frame_equal(result, expected)


def test_categorical_coerces_timestamp(all_parsers):
    parser = all_parsers
    dtype = {"b": CategoricalDtype([Timestamp("2014")])}

    data = "b\n2014-01-01\n2014-01-01T00:00:00"
    expected = DataFrame({"b": Categorical([Timestamp("2014")] * 2)})

    result = parser.read_csv(StringIO(data), dtype=dtype)
    tm.assert_frame_equal(result, expected)


def test_categorical_coerces_timedelta(all_parsers):
    parser = all_parsers
    dtype = {"b": CategoricalDtype(pd.to_timedelta(["1H", "2H", "3H"]))}

    data = "b\n1H\n2H\n3H"
    expected = DataFrame({"b": Categorical(dtype["b"].categories)})

    result = parser.read_csv(StringIO(data), dtype=dtype)
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize(
    "data",
    [
        "b\nTrue\nFalse\nNA\nFalse",
        "b\ntrue\nfalse\nNA\nfalse",
        "b\nTRUE\nFALSE\nNA\nFALSE",
        "b\nTrue\nFalse\nNA\nFALSE",
    ],
)
def test_categorical_dtype_coerces_boolean(all_parsers, data):
    # see gh-20498
    parser = all_parsers
    dtype = {"b": CategoricalDtype([False, True])}
    expected = DataFrame({"b": Categorical([True, False, None, False])})

    result = parser.read_csv(StringIO(data), dtype=dtype)
    tm.assert_frame_equal(result, expected)


def test_categorical_unexpected_categories(all_parsers):
    parser = all_parsers
    dtype = {"b": CategoricalDtype(["a", "b", "d", "e"])}

    data = "b\nd\na\nc\nd"  # Unexpected c
    expected = DataFrame({"b": Categorical(list("dacd"), dtype=dtype["b"])})

    result = parser.read_csv(StringIO(data), dtype=dtype)
    tm.assert_frame_equal(result, expected)


def test_empty_pass_dtype(all_parsers):
    parser = all_parsers

    data = "one,two"
    result = parser.read_csv(StringIO(data), dtype={"one": "u1"})

    expected = DataFrame(
        {"one": np.empty(0, dtype="u1"), "two": np.empty(0, dtype=np.object)},
        index=Index([], dtype=object),
    )
    tm.assert_frame_equal(result, expected)


def test_empty_with_index_pass_dtype(all_parsers):
    parser = all_parsers

    data = "one,two"
    result = parser.read_csv(
        StringIO(data), index_col=["one"], dtype={"one": "u1", 1: "f"}
    )

    expected = DataFrame(
        {"two": np.empty(0, dtype="f")}, index=Index([], dtype="u1", name="one")
    )
    tm.assert_frame_equal(result, expected)


def test_empty_with_multi_index_pass_dtype(all_parsers):
    parser = all_parsers

    data = "one,two,three"
    result = parser.read_csv(
        StringIO(data), index_col=["one", "two"], dtype={"one": "u1", 1: "f8"}
    )

    exp_idx = MultiIndex.from_arrays(
        [np.empty(0, dtype="u1"), np.empty(0, dtype=np.float64)], names=["one", "two"]
    )
    expected = DataFrame({"three": np.empty(0, dtype=np.object)}, index=exp_idx)
    tm.assert_frame_equal(result, expected)


def test_empty_with_mangled_column_pass_dtype_by_names(all_parsers):
    parser = all_parsers

    data = "one,one"
    result = parser.read_csv(StringIO(data), dtype={"one": "u1", "one.1": "f"})

    expected = DataFrame(
        {"one": np.empty(0, dtype="u1"), "one.1": np.empty(0, dtype="f")},
        index=Index([], dtype=object),
    )
    tm.assert_frame_equal(result, expected)


def test_empty_with_mangled_column_pass_dtype_by_indexes(all_parsers):
    parser = all_parsers

    data = "one,one"
    result = parser.read_csv(StringIO(data), dtype={0: "u1", 1: "f"})

    expected = DataFrame(
        {"one": np.empty(0, dtype="u1"), "one.1": np.empty(0, dtype="f")},
        index=Index([], dtype=object),
    )
    tm.assert_frame_equal(result, expected)


def test_empty_with_dup_column_pass_dtype_by_indexes(all_parsers):
    # see gh-9424
    parser = all_parsers
    expected = concat(
        [Series([], name="one", dtype="u1"), Series([], name="one.1", dtype="f")],
        axis=1,
    )
    expected.index = expected.index.astype(object)

    data = "one,one"
    result = parser.read_csv(StringIO(data), dtype={0: "u1", 1: "f"})
    tm.assert_frame_equal(result, expected)


def test_empty_with_dup_column_pass_dtype_by_indexes_raises(all_parsers):
    # see gh-9424
    parser = all_parsers
    expected = concat(
        [Series([], name="one", dtype="u1"), Series([], name="one.1", dtype="f")],
        axis=1,
    )
    expected.index = expected.index.astype(object)

    with pytest.raises(ValueError, match="Duplicate names"):
        data = ""
        parser.read_csv(StringIO(data), names=["one", "one"], dtype={0: "u1", 1: "f"})


def test_raise_on_passed_int_dtype_with_nas(all_parsers):
    # see gh-2631
    parser = all_parsers
    data = """YEAR, DOY, a
2001,106380451,10
2001,,11
2001,106380451,67"""

    msg = (
        "Integer column has NA values"
        if parser.engine == "c"
        else "Unable to convert column DOY"
    )
    with pytest.raises(ValueError, match=msg):
        parser.read_csv(StringIO(data), dtype={"DOY": np.int64}, skipinitialspace=True)


def test_dtype_with_converters(all_parsers):
    parser = all_parsers
    data = """a,b
1.1,2.2
1.2,2.3"""

    # Dtype spec ignored if converted specified.
    with tm.assert_produces_warning(ParserWarning):
        result = parser.read_csv(
            StringIO(data), dtype={"a": "i8"}, converters={"a": lambda x: str(x)}
        )
    expected = DataFrame({"a": ["1.1", "1.2"], "b": [2.2, 2.3]})
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize(
    "dtype,expected",
    [
        (np.float64, DataFrame(columns=["a", "b"], dtype=np.float64)),
        ("category", DataFrame({"a": Categorical([]), "b": Categorical([])}, index=[])),
        (
            dict(a="category", b="category"),
            DataFrame({"a": Categorical([]), "b": Categorical([])}, index=[]),
        ),
        ("datetime64[ns]", DataFrame(columns=["a", "b"], dtype="datetime64[ns]")),
        (
            "timedelta64[ns]",
            DataFrame(
                {
                    "a": Series([], dtype="timedelta64[ns]"),
                    "b": Series([], dtype="timedelta64[ns]"),
                },
                index=[],
            ),
        ),
        (
            dict(a=np.int64, b=np.int32),
            DataFrame(
                {"a": Series([], dtype=np.int64), "b": Series([], dtype=np.int32)},
                index=[],
            ),
        ),
        (
            {0: np.int64, 1: np.int32},
            DataFrame(
                {"a": Series([], dtype=np.int64), "b": Series([], dtype=np.int32)},
                index=[],
            ),
        ),
        (
            {"a": np.int64, 1: np.int32},
            DataFrame(
                {"a": Series([], dtype=np.int64), "b": Series([], dtype=np.int32)},
                index=[],
            ),
        ),
    ],
)
def test_empty_dtype(all_parsers, dtype, expected):
    # see gh-14712
    parser = all_parsers
    data = "a,b"

    result = parser.read_csv(StringIO(data), header=0, dtype=dtype)
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize(
    "dtype", list(np.typecodes["AllInteger"] + np.typecodes["Float"])
)
def test_numeric_dtype(all_parsers, dtype):
    data = "0\n1"
    parser = all_parsers
    expected = DataFrame([0, 1], dtype=dtype)

    result = parser.read_csv(StringIO(data), header=None, dtype=dtype)
    tm.assert_frame_equal(expected, result)


def test_boolean_dtype(all_parsers):
    parser = all_parsers
    data = "\n".join(
        [
            "a",
            "True",
            "TRUE",
            "true",
            "False",
            "FALSE",
            "false",
            "NaN",
            "nan",
            "NA",
            "null",
            "NULL",
        ]
    )

    result = parser.read_csv(StringIO(data), dtype="boolean")
    expected = pd.DataFrame(
        {
            "a": pd.array(
                [True, True, True, False, False, False, None, None, None, None, None],
                dtype="boolean",
            )
        }
    )

    tm.assert_frame_equal(result, expected)
