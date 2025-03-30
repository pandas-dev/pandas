from __future__ import annotations

import os

import pytest

from dask.dataframe._compat import PYARROW_GE_1500
from dask.dataframe.dask_expr import read_parquet
from dask.dataframe.dask_expr.tests._util import _backend_library, assert_eq

distributed = pytest.importorskip("distributed")

from distributed.utils_test import client as c  # noqa F401
from distributed.utils_test import gen_cluster

import dask.dataframe as dd

pd = _backend_library()


@pytest.fixture(params=["arrow"])
def filesystem(request):
    return request.param


def _make_file(dir, df=None, filename="myfile.parquet", **kwargs):
    fn = os.path.join(str(dir), filename)
    if df is None:
        df = pd.DataFrame({c: range(10) for c in "abcde"})
    df.to_parquet(fn, **kwargs)
    return fn


@gen_cluster(client=True)
async def test_io_fusion_merge(c, s, a, b, tmpdir):
    pdf = pd.DataFrame({c: range(100) for c in "abcdefghij"})
    p = dd.from_pandas(pdf, 2).to_parquet(tmpdir, compute=False)
    _ = await c.gather(c.compute(p))

    df = dd.read_parquet(tmpdir).merge(
        dd.read_parquet(tmpdir).add_suffix("_x"), left_on="a", right_on="a_x"
    )[["a_x", "b_x", "b"]]
    out = await c.gather(c.compute(df.optimize()))
    pd.testing.assert_frame_equal(
        out.sort_values(by="a_x", ignore_index=True),
        pdf.merge(pdf.add_suffix("_x"), left_on="a", right_on="a_x")[
            ["a_x", "b_x", "b"]
        ],
    )


@pytest.mark.skipif(not PYARROW_GE_1500, reason="requires 15.0.0")
@pytest.mark.filterwarnings("error")
@gen_cluster(client=True)
async def test_parquet_distriuted(c, s, a, b, tmpdir, filesystem):
    pdf = pd.DataFrame({"x": [1, 4, 3, 2, 0, 5]})
    df = read_parquet(_make_file(tmpdir, df=pdf), filesystem=filesystem)
    assert_eq(await c.gather(c.compute(df.optimize())), pdf)


@pytest.mark.skipif(not PYARROW_GE_1500, reason="requires 15.0.0")
def test_pickle_size(tmpdir, filesystem):
    pdf = pd.DataFrame({"x": [1, 4, 3, 2, 0, 5]})
    [_make_file(tmpdir, df=pdf, filename=f"{x}.parquet") for x in range(10)]
    df = read_parquet(tmpdir, filesystem=filesystem)
    from distributed.protocol import dumps

    assert (low_level := len(b"".join(dumps(df.optimize().dask)))) <= 10_000
    assert len(b"".join(dumps(df.optimize()))) <= low_level
