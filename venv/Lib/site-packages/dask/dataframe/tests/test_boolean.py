from __future__ import annotations

import pandas as pd

import dask.dataframe as dd


def test_meta():
    values = pd.array([True, False, None], dtype="boolean")
    ds = dd.from_pandas(pd.Series(values), 2)
    assert ds.dtype == pd.BooleanDtype()

    dd.utils.assert_eq(ds._meta_nonempty, pd.Series([True, pd.NA], dtype="boolean"))

    ddf = dd.from_pandas(pd.DataFrame({"A": values}), 2)
    assert ddf.dtypes["A"] == pd.BooleanDtype()

    dd.utils.assert_eq(
        ddf._meta_nonempty,
        pd.DataFrame({"A": pd.array([True, pd.NA], dtype="boolean")}),
    )


def test_ops():
    s1 = pd.Series(pd.array([True, False, None] * 3, dtype="boolean"))
    s2 = pd.Series(pd.array([True] * 3 + [False] * 3 + [None] * 3, dtype="boolean"))

    ds1 = dd.from_pandas(s1, 2)
    ds2 = dd.from_pandas(s2, 2)

    dd.utils.assert_eq(ds1 | ds2, s1 | s2)
    dd.utils.assert_eq(ds1 & ds2, s1 & s2)
    dd.utils.assert_eq(ds1 ^ ds2, s1 ^ s2)
