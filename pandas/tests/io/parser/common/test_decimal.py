"""
Tests that work on both the Python and C engines but do not have a
specific classification into the other test modules.
"""
from io import StringIO

import pytest

from pandas import DataFrame
import pandas._testing as tm


@pytest.mark.parametrize(
    "data,thousands,decimal",
    [
        (
            """A|B|C
1|2,334.01|5
10|13|10.
""",
            ",",
            ".",
        ),
        (
            """A|B|C
1|2.334,01|5
10|13|10,
""",
            ".",
            ",",
        ),
    ],
)
def test_1000_sep_with_decimal(all_parsers, data, thousands, decimal):
    parser = all_parsers
    expected = DataFrame({"A": [1, 10], "B": [2334.01, 13], "C": [5, 10.0]})

    result = parser.read_csv(
        StringIO(data), sep="|", thousands=thousands, decimal=decimal
    )
    tm.assert_frame_equal(result, expected)


def test_euro_decimal_format(all_parsers):
    parser = all_parsers
    data = """Id;Number1;Number2;Text1;Text2;Number3
1;1521,1541;187101,9543;ABC;poi;4,738797819
2;121,12;14897,76;DEF;uyt;0,377320872
3;878,158;108013,434;GHI;rez;2,735694704"""

    result = parser.read_csv(StringIO(data), sep=";", decimal=",")
    expected = DataFrame(
        [
            [1, 1521.1541, 187101.9543, "ABC", "poi", 4.738797819],
            [2, 121.12, 14897.76, "DEF", "uyt", 0.377320872],
            [3, 878.158, 108013.434, "GHI", "rez", 2.735694704],
        ],
        columns=["Id", "Number1", "Number2", "Text1", "Text2", "Number3"],
    )
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize("thousands", [None, "."])
@pytest.mark.parametrize(
    "value",
    ["e11,2", "1e11,2", "1,2,2", "1,2.1", "1,2e-10e1", "--1,2", "1a.2,1", "1..2,3"],
)
def test_decimal_and_exponential_erroneous(all_parsers, thousands, value):
    # GH#31920
    data = StringIO(
        f"""a	b
    1,1	{value}
    """
    )
    parser = all_parsers
    result = parser.read_csv(data, "\t", decimal=",", thousands=thousands)
    expected = DataFrame({"a": [1.1], "b": [value]})
    tm.assert_frame_equal(result, expected)
