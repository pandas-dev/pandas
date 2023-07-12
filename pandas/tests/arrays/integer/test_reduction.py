import numpy as np
import pytest

import pandas as pd
import pandas._testing as tm


@pytest.mark.parametrize("op", ["sum", "min", "max", "prod"])
def test_preserve_dtypes(op):
    # TODO(#22346): preserve Int64 dtype
    # for ops that enable (mean would actually work here
    # but generally it is a float return value)
    df = pd.DataFrame(
        {
            "A": ["a", "b", "b"],
            "B": [1, None, 3],
            "C": pd.array([1, None, 3], dtype="Int64"),
        }
    )

    # op
    result = getattr(df.C, op)()
    assert isinstance(result, np.int64)

    # groupby
    result = getattr(df.groupby("A"), op)()

    expected = pd.DataFrame(
        {"B": np.array([1.0, 3.0]), "C": pd.array([1, 3], dtype="Int64")},
        index=pd.Index(["a", "b"], name="A"),
    )
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize("op", ["sum", "prod", "min", "max"])
def test_dataframe_reductions(op):
    # https://github.com/pandas-dev/pandas/pull/32867
    # ensure the integers are not cast to float during reductions
    df = pd.DataFrame({"a": pd.array([1, 2], dtype="Int64")})
    result = getattr(df, op)()
    assert isinstance(result["a"], np.int64)
