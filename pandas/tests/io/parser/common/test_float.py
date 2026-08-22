"""
Tests that work on both the Python and C engines but do not have a
specific classification into the other test modules.
"""

from io import StringIO

import numpy as np
import pytest

from pandas.errors import Pandas4Warning

from pandas import DataFrame
import pandas._testing as tm

pytestmark = pytest.mark.filterwarnings(
    "ignore:Passing a BlockManager to DataFrame:DeprecationWarning"
)
skip_pyarrow = pytest.mark.usefixtures("pyarrow_skip")

depr_msg = "float_precision"


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

    warn = Pandas4Warning if precision is not None else None
    with tm.assert_produces_warning(warn, match=depr_msg, check_stacklevel=False):
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
    warn = Pandas4Warning if precision is not None else None
    with tm.assert_produces_warning(warn, match=depr_msg, check_stacklevel=False):
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
    warn = Pandas4Warning if precision is not None else None
    with tm.assert_produces_warning(warn, match=depr_msg, check_stacklevel=False):
        result = parser.read_csv(StringIO(data), float_precision=precision)
    expected = DataFrame({"data": [42.0, expected_value]})

    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize(
    "value",
    [
        "90071992547409930.0",
        "90071992547409931.0",
        "90071992547409935.0",
        "90071992547409938.0",
    ],
)
def test_precise_xstrtod_large_mantissa(c_parser_only, value):
    # GH#64357
    # When a 17-digit mantissa's 16-digit prefix crosses 2^53
    # (= 9007199254740992), the old per-digit FP accumulation
    #   number = number * 10. + digit
    # introduced a rounding error at step 16 that shifted the final result
    # by one ULP.  Specifically, the 16-digit prefix 9007199254740993
    # (= 2^53 + 1) was not representable as a double and rounded down to
    # 9007199254740992, causing subsequent arithmetic to produce a value
    # ~10 units below the true mantissa.  With ULP = 16 near 9e16, this
    # 10-unit error was enough to round to the wrong double.
    #
    # The fix accumulates mantissa digits in uint64_t and converts to
    # double once, so there is at most one rounding instead of up to
    # max_digits=17.
    parser = c_parser_only
    data = f"val\n{value}"
    with tm.assert_produces_warning(
        Pandas4Warning, match=depr_msg, check_stacklevel=False
    ):
        result = parser.read_csv(StringIO(data), float_precision="high")["val"][0]
    assert result == float(value)


@pytest.mark.parametrize(
    "value",
    [
        "000000000010084566.0",
        # leading zeros before the decimal point must not eat into the budget
        # either, or the low digits of the fractional part are lost
        "00000000000.00564728024629623",
        "000.0469911460454665528",
        "000000000000000000000000.03704296608561985",
        "0.12345678901234567890",
    ],
)
def test_precise_xstrtod_leading_zeros(c_parser_only, value):
    # GH#64184
    parser = c_parser_only
    data = f"val\n{value}\n"
    result = parser.read_csv(StringIO(data), thousands=",")
    expected = DataFrame({"val": [float(value)]})
    # check_exact: the digits lost to the bug are well inside the default rtol
    tm.assert_frame_equal(result, expected, check_exact=True)


@pytest.mark.parametrize("value", ["0", "000", "0.", "0.0", "-0", "0e5"])
def test_precise_xstrtod_all_zero_mantissa(c_parser_only, value):
    # GH#64184: an all-zero mantissa is a valid zero, not an unparsable string.
    # Guards the branch that skipping the leading zeros makes necessary.
    parser = c_parser_only
    data = f"val\n{value}\n"
    result = parser.read_csv(StringIO(data), thousands=",")
    expected = DataFrame({"val": [float(value)]})
    tm.assert_frame_equal(result, expected, check_dtype=False, check_exact=True)


def test_precise_xstrtod_leading_zero_matches_bare_point(c_parser_only):
    # GH#64184: "0.x" must get the same significant-digit budget as ".x"
    parser = c_parser_only
    digits = "12345678901234567890"
    data = f"val\n0.{digits}\n.{digits}\n"
    result = parser.read_csv(StringIO(data), thousands=",")["val"]
    assert result[0] == result[1] == float(f"0.{digits}")


@pytest.mark.parametrize(
    "value",
    [
        # up to 19 significant digits, parsed straight from the digits
        "1.7976931348623157e308",  # largest finite double
        "2.2250738585072014e-308",  # smallest normal double
        "9007199254740993.0",  # 2**53 + 1, not representable; rounds to 2**53
        "12345.678901234567",
        # more than 19 significant digits, so the mantissa has to be truncated
        # before it can be rounded
        "2.22507385850720113605740979670913197593481954635164565e-308",
        "1.00000000000000000000000000000000000000000000000000001",
        "123456789012345678901234567890.5",
        # subnormals, where picking the last bit requires comparing against the
        # full decimal expansion rather than a truncated mantissa
        "4.9406564584124654e-324",  # smallest positive subnormal
        "2.4703282292062327e-324",  # just under the halfway point -> 0.0
        "2.4703282292062328e-324",  # just over it -> smallest subnormal
    ],
)
def test_float_correctly_rounded(all_parsers, value):
    # GH#66457
    parser = all_parsers
    result = parser.read_csv(StringIO(f"data\n{value}"))

    expected = DataFrame({"data": [float(value)]})
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize(
    "value", ["81e31d04049863b72", "d81e31d04049863b72", "81e3104049863b72"]
)
def test_invalid_float_number(all_parsers_all_precisions, value):
    # GH#62617
    parser, precision = all_parsers_all_precisions
    data = f"h1,h2,h3\ndata1,{value},data3"

    warn = Pandas4Warning if precision is not None else None
    with tm.assert_produces_warning(warn, match=depr_msg, check_stacklevel=False):
        result = parser.read_csv(StringIO(data), float_precision=precision)
    expected = DataFrame({"h1": ["data1"], "h2": [value], "h3": "data3"})
    tm.assert_frame_equal(result, expected)
