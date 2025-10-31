"""
Tests that work on both the Python and C engines but do not have a
specific classification into the other test modules.
"""

from io import StringIO

import numpy as np
import pytest

from pandas import DataFrame
import pandas._testing as tm

pytestmark = pytest.mark.filterwarnings(
    "ignore:Passing a BlockManager to DataFrame:DeprecationWarning"
)
skip_pyarrow = pytest.mark.usefixtures("pyarrow_skip")


@skip_pyarrow  # ParserError: CSV parse error: Empty CSV file or block
def test_float_parser(all_parsers):
    # see gh-9565
    parser = all_parsers
    data = "45e-1,4.5,45.,inf,-inf"
    result = parser.read_csv(StringIO(data), header=None)

    expected = DataFrame([[float(s) for s in data.split(",")]])
    tm.assert_frame_equal(result, expected)


def test_scientific_no_exponent(all_parsers_all_precisions):
    # see gh-12215
    df = DataFrame.from_dict({"w": ["2e"], "x": ["3E"], "y": ["42e"], "z": ["632E"]})
    data = df.to_csv(index=False)
    parser, precision = all_parsers_all_precisions

    df_roundtrip = parser.read_csv(StringIO(data), float_precision=precision)
    tm.assert_frame_equal(df_roundtrip, df)


@pytest.mark.parametrize(
    "value, expected_value",
    [
        ("0E-617", 0.0),
        ("0E99999999", 0.0),
        ("-0E99999999", 0.0),
        ("-0E-99999999", 0.0),
        ("10E-617", 0.0),
        ("10E-100000", 0.0),
        ("-10E-100000", 0.0),
        ("10e-99999999999", 0.0),
        ("10e-999999999999", 0.0),
        ("10e-9999999999999", 0.0),
        ("10E999", np.inf),
        ("-10e99999999999", -np.inf),
        ("10e99999999999", np.inf),
        ("10e999999999999", np.inf),
        ("10e9999999999999", np.inf),
        ("50060e8007123400", np.inf),
        ("-50060e8007123400", -np.inf),
    ],
)
def test_large_exponent(all_parsers_all_precisions, value, expected_value):
    # GH#38753; GH#38794; GH#62740
    parser, precision = all_parsers_all_precisions

    data = f"data\n{value}"
    result = parser.read_csv(StringIO(data), float_precision=precision)
    expected = DataFrame({"data": [expected_value]})
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize(
    "value, expected_value",
    [
        ("32.0", 32.0),
        ("32e0", 32.0),
        ("3.2e1", 32.0),
        ("3.2e80", 3.2e80),
        ("3.2e-80", 3.2e-80),
        ("18446744073709551616.0", float(1 << 64)),  # loses precision
        ("18446744073709551616.5", float(1 << 64)),  # loses precision
        ("36893488147419103232.3", float(1 << 65)),  # loses precision
    ],
)
def test_small_int_followed_by_float(
    all_parsers_all_precisions, value, expected_value, request
):
    # GH#51295
    parser, precision = all_parsers_all_precisions
    data = f"""data
    42
    {value}"""
    result = parser.read_csv(StringIO(data), float_precision=precision)
    expected = DataFrame({"data": [42.0, expected_value]})

    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize(
    "value", ["81e31d04049863b72", "d81e31d04049863b72", "81e3104049863b72"]
)
def test_invalid_float_number(all_parsers_all_precisions, value):
    # GH#62617
    parser, precision = all_parsers_all_precisions
    data = f"h1,h2,h3\ndata1,{value},data3"

    result = parser.read_csv(StringIO(data), float_precision=precision)
    expected = DataFrame({"h1": ["data1"], "h2": [value], "h3": "data3"})
    tm.assert_frame_equal(result, expected)
