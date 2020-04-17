import pytest

import pandas as pd
import pandas._testing as tm


@pytest.mark.parametrize(
    "to_concat_dtypes, result_dtype",
    [
        (["Int64", "Int64"], "Int64"),
        (["UInt64", "UInt64"], "UInt64"),
        (["Int8", "Int8"], "Int8"),
        (["Int8", "Int16"], "Int16"),
        (["UInt8", "Int8"], "Int16"),
        (["Int32", "UInt32"], "Int64"),
        # this still gives object (awaiting float extension dtype)
        (["Int64", "UInt64"], "object"),
    ],
)
def test_concat_series(to_concat_dtypes, result_dtype):

    result = pd.concat([pd.Series([1, 2, pd.NA], dtype=t) for t in to_concat_dtypes])
    expected = pd.concat([pd.Series([1, 2, pd.NA], dtype=object)] * 2).astype(
        result_dtype
    )
    tm.assert_series_equal(result, expected)
