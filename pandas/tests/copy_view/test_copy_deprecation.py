import pytest

from pandas.errors import Pandas4Warning

import pandas as pd
from pandas import (
    concat,
    merge,
)
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
        ("transpose", {}),
        ("set_axis", {"labels": [1, 2, 3]}),
        ("rename", {"index": {1: 2}}),
        ("set_flags", {}),
        ("to_period", {}),
        ("to_timestamp", {}),
        ("swaplevel", {"i": 0, "j": 1}),
    ],
)
def test_copy_deprecation(meth, kwargs):
    df = pd.DataFrame({"a": [1, 2, 3], "b": [4, 5, 6], "c": 1})

    if meth in ("tz_convert", "tz_localize", "to_period"):
        tz = None if meth in ("tz_localize", "to_period") else "US/Eastern"
        df.index = pd.date_range("2020-01-01", freq="D", periods=len(df), tz=tz)
    elif meth == "to_timestamp":
        df.index = pd.period_range("2020-01-01", freq="D", periods=len(df))
    elif meth == "swaplevel":
        df = df.set_index(["b", "c"])

    if meth != "swaplevel":
        with tm.assert_produces_warning(Pandas4Warning, match="copy"):
            getattr(df, meth)(copy=False, **kwargs)

    if meth != "transpose":
        with tm.assert_produces_warning(Pandas4Warning, match="copy"):
            getattr(df.a, meth)(copy=False, **kwargs)


def test_copy_deprecation_reindex_like_align():
    df = pd.DataFrame({"a": [1, 2, 3], "b": [4, 5, 6]})
    # Somehow the stack level check is incorrect here
    with tm.assert_produces_warning(
        Pandas4Warning, match="copy", check_stacklevel=False
    ):
        df.reindex_like(df, copy=False)

    with tm.assert_produces_warning(
        Pandas4Warning, match="copy", check_stacklevel=False
    ):
        df.a.reindex_like(df.a, copy=False)

    with tm.assert_produces_warning(
        Pandas4Warning, match="copy", check_stacklevel=False
    ):
        df.align(df, copy=False)

    with tm.assert_produces_warning(
        Pandas4Warning, match="copy", check_stacklevel=False
    ):
        df.a.align(df.a, copy=False)


def test_copy_deprecation_merge_concat():
    df = pd.DataFrame({"a": [1, 2, 3], "b": [4, 5, 6]})

    with tm.assert_produces_warning(
        Pandas4Warning, match="copy", check_stacklevel=False
    ):
        df.merge(df, copy=False)

    with tm.assert_produces_warning(
        Pandas4Warning, match="copy", check_stacklevel=False
    ):
        merge(df, df, copy=False)

    with tm.assert_produces_warning(
        Pandas4Warning, match="copy", check_stacklevel=False
    ):
        concat([df, df], copy=False)


@pytest.mark.parametrize("value", [False, True, "warn"])
def test_copy_on_write_deprecation_option(value):
    msg = "Copy-on-Write can no longer be disabled"
    # stacklevel points to contextlib due to use of context manager.
    with tm.assert_produces_warning(Pandas4Warning, match=msg, check_stacklevel=False):
        with pd.option_context("mode.copy_on_write", value):
            pass
