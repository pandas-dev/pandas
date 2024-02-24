import pytest

import pandas as pd
import pandas._testing as tm


@pytest.mark.parametrize(
    "meth, kwargs",
    [
        ("truncate", {}),
        ("tz_convert", {"tz": "UTC"}),
        ("tz_localize", {"tz": "UTC"}),
        ("infer_objects", {}),
        ("astype", {"dtype": "float64"}),
        ("reindex", {"index": [2, 0, 1]}),
    ],
)
def test_copy_deprecation(meth, kwargs):
    df = pd.DataFrame({"a": [1, 2, 3], "b": [4, 5, 6]})

    if meth in ("tz_convert", "tz_localize"):
        tz = None if meth == "tz_localize" else "US/Eastern"
        df.index = pd.date_range("2020-01-01", freq="D", periods=len(df), tz=tz)

    with tm.assert_produces_warning(DeprecationWarning, match="copy"):
        getattr(df, meth)(copy=False, **kwargs)

    with tm.assert_produces_warning(DeprecationWarning, match="copy"):
        getattr(df.a, meth)(copy=False, **kwargs)


def test_copy_deprecation_reindex_like_align():
    df = pd.DataFrame({"a": [1, 2, 3], "b": [4, 5, 6]})
    # Somehow the stack level check is incorrect here
    with tm.assert_produces_warning(
        DeprecationWarning, match="copy", check_stacklevel=False
    ):
        df.reindex_like(df, copy=False)

    with tm.assert_produces_warning(
        DeprecationWarning, match="copy", check_stacklevel=False
    ):
        df.a.reindex_like(df.a, copy=False)

    with tm.assert_produces_warning(
        DeprecationWarning, match="copy", check_stacklevel=False
    ):
        df.align(df, copy=False)

    with tm.assert_produces_warning(
        DeprecationWarning, match="copy", check_stacklevel=False
    ):
        df.a.align(df.a, copy=False)
