"""
Tests dtype specification during parsing
for all of the parsers defined in parsers.py
"""
from io import StringIO

import numpy as np
import pytest

from pandas.errors import ParserWarning

import pandas as pd
from pandas import DataFrame
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
            "1",
            "1.0",
            "False",
            "FALSE",
            "false",
            "0",
            "0.0",
            "NaN",
            "nan",
            "NA",
            "null",
            "NULL",
        ]
    )

    result = parser.read_csv(StringIO(data), dtype="boolean")
    expected = DataFrame(
        {
            "a": pd.array(
                [
                    True,
                    True,
                    True,
                    True,
                    True,
                    False,
                    False,
                    False,
                    False,
                    False,
                    None,
                    None,
                    None,
                    None,
                    None,
                ],
                dtype="boolean",
            )
        }
    )

    tm.assert_frame_equal(result, expected)
