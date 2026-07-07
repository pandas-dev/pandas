"""
Tests that work on both the Python and C engines but do not have a
specific classification into the other test modules.
"""

from io import StringIO

import pytest

from pandas import DataFrame
import pandas._testing as tm

pytestmark = pytest.mark.filterwarnings(
    "ignore:Passing a BlockManager to DataFrame:DeprecationWarning"
)


@pytest.mark.parametrize("na_filter", [True, False])
def test_inf_parsing(all_parsers, na_filter):
    parser = all_parsers
    data = """\
,A
a,inf
b,-inf
c,+Inf
d,-Inf
e,INF
f,-INF
g,+INf
h,-INf
i,inF
j,-inF"""
    expected = DataFrame(
        {"A": [float("inf"), float("-inf")] * 5},
        index=["a", "b", "c", "d", "e", "f", "g", "h", "i", "j"],
    )

    if parser.engine == "pyarrow" and not na_filter:
        # the pyarrow engine does not support na_filter=False
        msg = "The 'na_filter' option is not supported with the 'pyarrow' engine"
        with pytest.raises(ValueError, match=msg):
            parser.read_csv(StringIO(data), index_col=0, na_filter=na_filter)
        return

    result = parser.read_csv(StringIO(data), index_col=0, na_filter=na_filter)
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize("na_filter", [True, False])
def test_infinity_parsing(all_parsers, na_filter):
    parser = all_parsers
    data = """\
,A
a,Infinity
b,-Infinity
c,+Infinity
"""
    expected = DataFrame(
        {"A": [float("infinity"), float("-infinity"), float("+infinity")]},
        index=["a", "b", "c"],
    )

    if parser.engine == "pyarrow" and not na_filter:
        # the pyarrow engine does not support na_filter=False
        msg = "The 'na_filter' option is not supported with the 'pyarrow' engine"
        with pytest.raises(ValueError, match=msg):
            parser.read_csv(StringIO(data), index_col=0, na_filter=na_filter)
        return

    result = parser.read_csv(StringIO(data), index_col=0, na_filter=na_filter)
    tm.assert_frame_equal(result, expected)
