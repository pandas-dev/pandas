"""
Tests dtype specification during parsing
for all of the parsers defined in parsers.py
"""

import csv
from io import StringIO

import numpy as np
import pytest

from pandas._libs import parsers as libparsers
from pandas.errors import Pandas4Warning

from pandas.core.dtypes.dtypes import CategoricalDtype

import pandas as pd
from pandas import (
    Categorical,
    DataFrame,
    Timestamp,
)
import pandas._testing as tm

pytestmark = pytest.mark.filterwarnings(
    "ignore:Passing a BlockManager to DataFrame:DeprecationWarning"
)

xfail_pyarrow = pytest.mark.usefixtures("pyarrow_xfail")


@pytest.mark.parametrize(
    "dtype",
    [
        "category",
        CategoricalDtype(),
        {"a": "category", "b": "category", "c": CategoricalDtype()},
    ],
)
def test_categorical_dtype(all_parsers, dtype):
    # see gh-10153, gh-56044
    # categories should have inferred types (not strings) across all engines
    parser = all_parsers
    data = """a,b,c
1,a,3.4
1,a,3.4
2,b,4.5"""
    expected = DataFrame(
        {
            "a": Categorical([1, 1, 2]),
            "b": Categorical(["a", "a", "b"]),
            "c": Categorical([3.4, 3.4, 4.5]),
        }
    )
    actual = parser.read_csv(StringIO(data), dtype=dtype)
    tm.assert_frame_equal(actual, expected)


@pytest.mark.parametrize("dtype", [{"b": "category"}, {1: "category"}])
def test_categorical_dtype_single(all_parsers, dtype, request):
    # see gh-10153
    parser = all_parsers
    data = """a,b,c
1,a,3.4
1,a,3.4
2,b,4.5"""
    expected = DataFrame(
        {"a": [1, 1, 2], "b": Categorical(["a", "a", "b"]), "c": [3.4, 3.4, 4.5]}
    )
    if parser.engine == "pyarrow" and any(isinstance(key, int) for key in dtype):
        mark = pytest.mark.xfail(
            reason="pyarrow doesn't support specifying dtype by column index",
        )
        request.applymarker(mark)

    actual = parser.read_csv(StringIO(data), dtype=dtype)
    tm.assert_frame_equal(actual, expected)


def test_categorical_dtype_unsorted(all_parsers):
    # see gh-10153, gh-56044
    parser = all_parsers
    data = """a,b,c
1,b,3.4
1,b,3.4
2,a,4.5"""
    expected = DataFrame(
        {
            "a": Categorical([1, 1, 2]),
            "b": Categorical(["b", "b", "a"]),
            "c": Categorical([3.4, 3.4, 4.5]),
        }
    )
    actual = parser.read_csv(StringIO(data), dtype="category")
    tm.assert_frame_equal(actual, expected)


def test_categorical_dtype_missing(all_parsers):
    # see gh-10153, gh-56044
    parser = all_parsers
    data = """a,b,c
1,b,3.4
1,nan,3.4
2,a,4.5"""
    expected = DataFrame(
        {
            "a": Categorical([1, 1, 2]),
            "b": Categorical(["b", np.nan, "a"]),
            "c": Categorical([3.4, 3.4, 4.5]),
        }
    )
    actual = parser.read_csv(StringIO(data), dtype="category")
    tm.assert_frame_equal(actual, expected)


def test_categorical_dtype_numeric_duplicates(all_parsers):
    # GH#56044 distinct strings that convert to the same number should be
    #  merged into a single category
    parser = all_parsers
    data = "a\n1\n1.0\n2"
    expected = DataFrame({"a": Categorical([1.0, 1.0, 2.0])})
    actual = parser.read_csv(StringIO(data), dtype="category")
    tm.assert_frame_equal(actual, expected)


@pytest.mark.parametrize(
    "data",
    [
        "a\nTrue\nFalse\nTrue",
        "a\ntrue\nfalse\ntrue",
        "a\nTRUE\nFALSE\nTRUE",
    ],
)
def test_categorical_dtype_infers_boolean(all_parsers, data):
    # GH#56044 boolean-looking columns infer bool categories across all
    #  engines, matching non-categorical parsing
    parser = all_parsers
    expected = DataFrame({"a": Categorical([True, False, True])})
    actual = parser.read_csv(StringIO(data), dtype="category")
    tm.assert_frame_equal(actual, expected)


def test_categorical_dtype_boolean_duplicates(all_parsers):
    # GH#56044 distinct strings that convert to the same bool should be
    #  merged into a single category
    parser = all_parsers
    data = "a\nTrue\nTRUE\ntrue\nFalse"
    expected = DataFrame({"a": Categorical([True, True, True, False])})
    actual = parser.read_csv(StringIO(data), dtype="category")
    tm.assert_frame_equal(actual, expected)


def test_categorical_dtype_boolean_custom_values(all_parsers):
    # GH#56044 true_values/false_values are honored when inferring bool
    #  categories
    parser = all_parsers
    data = "a\nyes\nno\nyes"
    expected = DataFrame({"a": Categorical([True, False, True])})
    actual = parser.read_csv(
        StringIO(data), dtype="category", true_values=["yes"], false_values=["no"]
    )
    tm.assert_frame_equal(actual, expected)


@xfail_pyarrow  # ValueError: The 'thousands' option is not supported
def test_categorical_dtype_thousands(all_parsers):
    # GH#56044 numeric inference is unaware of the thousands option,
    #  so categories stay strings rather than mis-parsing "1.000" as 1.0
    parser = all_parsers
    data = "a\n1.000\n2.000"
    expected = DataFrame({"a": Categorical(["1.000", "2.000"])})
    actual = parser.read_csv(StringIO(data), dtype="category", thousands=".")
    tm.assert_frame_equal(actual, expected)


@xfail_pyarrow  # pyarrow parses with decimal_point and infers numeric
def test_categorical_dtype_decimal(all_parsers):
    # GH#56044 numeric inference is unaware of the decimal option,
    #  so categories stay strings rather than mis-parsing "1,5"
    parser = all_parsers
    data = "a;b\n1;1,5\n2;2,5"
    expected = DataFrame(
        {"a": Categorical(["1", "2"]), "b": Categorical(["1,5", "2,5"])}
    )
    actual = parser.read_csv(StringIO(data), dtype="category", sep=";", decimal=",")
    tm.assert_frame_equal(actual, expected)


@xfail_pyarrow  # pyarrow casts to float64, losing precision
def test_categorical_dtype_large_integers(all_parsers):
    # GH#56044 integers too large for int64/uint64 give object categories
    parser = all_parsers
    data = "a\n99999999999999999999999999\n1"
    expected = DataFrame({"a": Categorical([99999999999999999999999999, 1])})
    actual = parser.read_csv(StringIO(data), dtype="category")
    tm.assert_frame_equal(actual, expected)


def test_categorical_dtype_explicit_integer_ea_categories(all_parsers):
    # GH#56136 explicitly-requested IntegerDtype categories are preserved
    parser = all_parsers
    cat_dtype = CategoricalDtype(pd.array([1, 2], dtype="Int64"))
    data = "a\n1\n2\n1"
    expected = DataFrame({"a": Categorical.from_codes([0, 1, 0], dtype=cat_dtype)})
    actual = parser.read_csv(StringIO(data), dtype={"a": cat_dtype})
    tm.assert_frame_equal(actual, expected)


def test_categorical_dtype_non_default_dtype_backend(all_parsers, dtype_backend):
    # GH#56044 categories are inferred, but the c and python engines do not yet
    #  honor dtype_backend and give numpy categories where a non-categorical
    #  read would give Int64/int64[pyarrow]; pin the current behavior, see
    #  GH#66382
    parser = all_parsers
    data = "a\n1\n2"
    result = parser.read_csv(
        StringIO(data), dtype="category", dtype_backend=dtype_backend
    )
    cat_dtype = result["a"].cat.categories.dtype
    if parser.engine == "pyarrow" and dtype_backend == "numpy_nullable":
        assert cat_dtype == pd.Int64Dtype()
    elif dtype_backend == "pyarrow" and parser.engine in ("pyarrow", "python"):
        # the python parser's arrow-backed strings convert straight to arrow
        pyarrow = pytest.importorskip("pyarrow")
        assert cat_dtype == pd.ArrowDtype(pyarrow.int64())
    else:
        assert cat_dtype == np.dtype("int64")


@xfail_pyarrow  # ValueError: The 'quoting' option is not supported
def test_categorical_dtype_quote_nonnumeric(all_parsers):
    # GH#56044 with QUOTE_NONNUMERIC, non-categorical columns parse as
    #  float64; the c engine currently infers int64 categories, while the
    #  python engine (whose tokenizer produces floats) gives float64
    parser = all_parsers
    data = '"a"\n1\n2'
    result = parser.read_csv(
        StringIO(data), dtype="category", quoting=csv.QUOTE_NONNUMERIC
    )
    if parser.engine == "python":
        expected = DataFrame({"a": Categorical([1.0, 2.0])})
    else:
        expected = DataFrame({"a": Categorical([1, 2])})
    tm.assert_frame_equal(result, expected)


def test_categorical_dtype_low_memory_mixed_numeric_chunks(all_parsers, monkeypatch):
    # GH#56044 in low-memory mode chunks can disagree on the inferred category
    #  dtype (int64 until "1.5" appears, then float64), which used to break
    #  union_categoricals; inference now happens after concatenation
    parser = all_parsers
    heuristic = 2**5
    ints = [str(i) for i in range(40)]
    rows = [*ints, "1.5", *ints]
    expected = DataFrame({"a": Categorical([float(x) for x in rows])})
    with monkeypatch.context() as m:
        m.setattr(libparsers, "DEFAULT_BUFFER_HEURISTIC", heuristic)
        actual = parser.read_csv(StringIO("a\n" + "\n".join(rows)), dtype="category")
    tm.assert_frame_equal(actual, expected)


@pytest.mark.parametrize(
    "rows",
    [
        [*(str(i) for i in range(40)), "apple", *(str(i) for i in range(40))],
        [*["True", "False"] * 20, "maybe", *["True", "False"] * 20],
        [*["True", "False"] * 20, "3", *["True", "False"] * 20],
    ],
)
def test_categorical_dtype_low_memory_mixed_type_chunks(all_parsers, monkeypatch, rows):
    # GH#56044 a chunk whose values do not all parse as numeric or boolean
    #  keeps every chunk's categories strings
    parser = all_parsers
    heuristic = 2**5
    expected = DataFrame({"a": Categorical(rows)})
    with monkeypatch.context() as m:
        m.setattr(libparsers, "DEFAULT_BUFFER_HEURISTIC", heuristic)
        actual = parser.read_csv(StringIO("a\n" + "\n".join(rows)), dtype="category")
    # category order for unconverted string chunks follows chunk-union
    #  order; normalize to sorted before comparing
    actual["a"] = actual["a"].cat.reorder_categories(
        actual["a"].cat.categories.sort_values()
    )
    tm.assert_frame_equal(actual, expected)


@xfail_pyarrow  # pyarrow treats "" as null regardless of na_filter
def test_categorical_dtype_empty_string_na_filter_false(all_parsers):
    # GH#56044 to_numeric converts "" to NaN, which cannot be a category;
    #  keep string categories
    parser = all_parsers
    data = "a,b\n,1\n2,3"
    expected = DataFrame({"a": Categorical(["", "2"]), "b": Categorical([1, 3])})
    actual = parser.read_csv(StringIO(data), dtype="category", na_filter=False)
    tm.assert_frame_equal(actual, expected)


def test_categorical_dtype_empty_string_keep_default_na(all_parsers):
    # GH#56044 as above, with "" surviving via keep_default_na=False
    parser = all_parsers
    data = "a,b\n,1\n2,3"
    expected = DataFrame({"a": Categorical(["", "2"]), "b": Categorical([1, 3])})
    actual = parser.read_csv(StringIO(data), dtype="category", keep_default_na=False)
    tm.assert_frame_equal(actual, expected)


@xfail_pyarrow  # pyarrow gives float64 categories for all-NA columns
def test_categorical_dtype_all_na(all_parsers):
    # GH#56044 empty inferred categories keep object dtype
    parser = all_parsers
    data = "a,b\n,1\n,2"
    expected = DataFrame(
        {
            "a": Categorical.from_codes([-1, -1], pd.Index([], dtype=object)),
            "b": Categorical([1, 2]),
        }
    )
    actual = parser.read_csv(StringIO(data), dtype="category")
    tm.assert_frame_equal(actual, expected)


@pytest.mark.slow
def test_categorical_dtype_high_cardinality_numeric(all_parsers, monkeypatch):
    # see gh-18186, gh-56044
    # was an issue with C parser, due to DEFAULT_BUFFER_HEURISTIC
    parser = all_parsers
    heuristic = 2**5
    data = np.sort([str(i) for i in range(heuristic + 1)])
    csv_data = "a\n" + "\n".join(data)
    int_data = np.array([int(x) for x in data])
    expected = DataFrame({"a": Categorical(int_data, ordered=True)})
    with monkeypatch.context() as m:
        m.setattr(libparsers, "DEFAULT_BUFFER_HEURISTIC", heuristic)
        actual = parser.read_csv(StringIO(csv_data), dtype="category")
    actual["a"] = actual["a"].cat.reorder_categories(
        np.sort(actual["a"].cat.categories), ordered=True
    )
    tm.assert_frame_equal(actual, expected)


def test_categorical_dtype_utf16(all_parsers, datapath):
    # see gh-10153
    pth = datapath("io", "parser", "data", "utf16_ex.txt")
    parser = all_parsers
    encoding = "utf-16"
    sep = "\t"

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

    if parser.engine == "pyarrow":
        msg = "The 'chunksize' option is not supported with the 'pyarrow' engine"
        with pytest.raises(ValueError, match=msg):
            parser.read_csv(StringIO(data), dtype={"b": "category"}, chunksize=2)
        return

    with parser.read_csv(
        StringIO(data), dtype={"b": "category"}, chunksize=2
    ) as actuals:
        for actual, expected in zip(actuals, expecteds, strict=True):
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
            {"a": [1, 2], "b": Categorical(["b", "c"], categories=cats)},
            index=[2, 3],
        ),
    ]
    dtype = CategoricalDtype(cats)

    if parser.engine == "pyarrow":
        msg = "The 'chunksize' option is not supported with the 'pyarrow' engine"
        with pytest.raises(ValueError, match=msg):
            parser.read_csv(StringIO(data), dtype={"b": dtype}, chunksize=2)
        return

    with parser.read_csv(StringIO(data), dtype={"b": dtype}, chunksize=2) as actuals:
        for actual, expected in zip(actuals, expecteds, strict=True):
            tm.assert_frame_equal(actual, expected)


def test_categorical_dtype_latin1(all_parsers, datapath):
    # see gh-10153
    pth = datapath("io", "parser", "data", "unicode_series.csv")
    parser = all_parsers
    encoding = "latin-1"

    expected = parser.read_csv(pth, header=None, encoding=encoding)
    expected[1] = Categorical(expected[1])

    actual = parser.read_csv(pth, header=None, encoding=encoding, dtype={1: "category"})
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
    dti = pd.DatetimeIndex(["2017-01-01", "2018-01-01", "2019-01-01"], freq=None)
    dtype = {"b": CategoricalDtype(dti)}

    data = "b\n2017-01-01\n2018-01-01\n2019-01-01"
    expected = DataFrame({"b": Categorical(dtype["b"].categories)})

    result = parser.read_csv(StringIO(data), dtype=dtype)
    tm.assert_frame_equal(result, expected)


def test_categorical_coerces_timestamp(all_parsers):
    parser = all_parsers
    dtype = {"b": CategoricalDtype([Timestamp("2014")])}

    data = "b\n2014-01-01\n2014-01-01"
    expected = DataFrame({"b": Categorical([Timestamp("2014")] * 2)})

    result = parser.read_csv(StringIO(data), dtype=dtype)
    tm.assert_frame_equal(result, expected)


def test_categorical_coerces_timedelta(all_parsers):
    parser = all_parsers
    dtype = {"b": CategoricalDtype(pd.to_timedelta(["1h", "2h", "3h"]))}

    data = "b\n1h\n2h\n3h"
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
    expected = DataFrame({"b": Categorical(["d", "a", None, "d"], dtype=dtype["b"])})

    msg = "Constructing a Categorical with a dtype and values containing"
    with tm.assert_produces_warning(Pandas4Warning, match=msg, check_stacklevel=False):
        result = parser.read_csv(StringIO(data), dtype=dtype)
    tm.assert_frame_equal(result, expected)
