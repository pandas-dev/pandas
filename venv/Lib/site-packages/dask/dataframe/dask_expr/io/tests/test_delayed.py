from __future__ import annotations

import numpy as np
import pytest

import dask
from dask import delayed
from dask.dataframe.dask_expr import from_delayed, from_dict
from dask.dataframe.dask_expr.tests._util import _backend_library, assert_eq

pd = _backend_library()


def test_from_delayed_optimizing():
    parts = from_dict({"a": np.arange(300)}, npartitions=30).to_delayed()
    result = from_delayed(parts[0], meta=pd.DataFrame({"a": pd.Series(dtype=np.int64)}))
    assert len(result.optimize().dask) == 2
    assert_eq(result, pd.DataFrame({"a": pd.Series(np.arange(10))}))


@pytest.mark.parametrize("prefix", [None, "foo"])
def test_from_delayed(prefix):
    pdf = pd.DataFrame(
        data=np.random.normal(size=(10, 4)), columns=["a", "b", "c", "d"]
    )
    parts = [pdf.iloc[:1], pdf.iloc[1:3], pdf.iloc[3:6], pdf.iloc[6:10]]
    dfs = [delayed(parts.__getitem__)(i) for i in range(4)]

    df = from_delayed(dfs, meta=pdf.head(0), divisions=None, prefix=prefix)
    assert_eq(df, pdf)
    assert len({k[0] for k in df.optimize().dask}) == 2
    if prefix:
        assert df._name.startswith(prefix)

    divisions = tuple([p.index[0] for p in parts] + [parts[-1].index[-1]])
    df = from_delayed(dfs, meta=pdf.head(0), divisions=divisions, prefix=prefix)
    assert_eq(df, pdf)
    if prefix:
        assert df._name.startswith(prefix)


def test_from_delayed_dask():
    df = pd.DataFrame(data=np.random.normal(size=(10, 4)), columns=list("abcd"))
    parts = [df.iloc[:1], df.iloc[1:3], df.iloc[3:6], df.iloc[6:10]]
    dfs = [delayed(parts.__getitem__)(i) for i in range(4)]
    meta = dfs[0].compute()

    my_len = lambda x: pd.Series([len(x)])

    for divisions in [None, [0, 1, 3, 6, 10]]:
        ddf = from_delayed(dfs, meta=meta, divisions=divisions)
        assert_eq(ddf, df)
        assert list(ddf.map_partitions(my_len).compute()) == [1, 2, 3, 4]
        assert ddf.known_divisions == (divisions is not None)

        s = from_delayed([d.a for d in dfs], meta=meta.a, divisions=divisions)
        assert_eq(s, df.a)
        assert list(s.map_partitions(my_len).compute()) == [1, 2, 3, 4]
        assert ddf.known_divisions == (divisions is not None)


@dask.delayed
def _load(x):
    return pd.DataFrame({"x": x, "y": [1, 2, 3]})


def test_from_delayed_fusion():
    func = lambda x: None
    df = from_delayed([_load(x) for x in range(10)], meta={"x": "int64", "y": "int64"})
    result = df.map_partitions(func, meta={}).optimize().dask
    expected = df.map_partitions(func, meta={}).optimize(fuse=False).dask
    assert result.keys() == expected.keys()

    expected = df.map_partitions(func, meta={}).lower_completely().dask
    assert result.keys() == expected.keys()
