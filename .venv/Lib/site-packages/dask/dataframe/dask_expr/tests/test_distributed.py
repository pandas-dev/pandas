from __future__ import annotations

import numpy as np
import pytest

import dask
from dask.dataframe.dask_expr import Repartition, from_pandas, map_partitions, merge
from dask.dataframe.dask_expr._merge import BroadcastJoin
from dask.dataframe.dask_expr._shuffle import P2PShuffle, TaskShuffle
from dask.dataframe.dask_expr.tests._util import _backend_library

distributed = pytest.importorskip("distributed")

from distributed import Client, LocalCluster, SchedulerPlugin
from distributed.shuffle._core import id_from_key
from distributed.utils_test import cleanup  # noqa F401
from distributed.utils_test import client as c  # noqa F401
from distributed.utils_test import gen_cluster, loop, loop_in_thread  # noqa F401

import dask.dataframe as dd

# Set DataFrame backend for this module
pd = _backend_library()


@pytest.fixture
def pdf():
    return pd.DataFrame({"x": list(range(20)) * 5, "y": range(100)})


@pytest.fixture
def df(pdf):
    return from_pandas(pdf, npartitions=10)


@pytest.mark.parametrize("npartitions", [None, 1, 20])
@gen_cluster(client=True)
async def test_p2p_shuffle(c, s, a, b, npartitions):
    df = dd.dask_expr.datasets.timeseries(
        start="2000-01-01",
        end="2000-01-10",
        dtypes={"x": float, "y": float},
        freq="10 s",
    )
    out = df.shuffle("x", shuffle_method="p2p", npartitions=npartitions)
    if npartitions is None:
        assert out.npartitions == df.npartitions
    else:
        assert out.npartitions == npartitions
    x, y, z = c.compute([df.x.size, out.x.size, out.partitions[-1].x.size])
    x = await x
    y = await y
    z = await z
    assert x == y
    if npartitions != 1:
        assert x > z


@pytest.mark.parametrize("npartitions_left", [5, 6])
@gen_cluster(client=True)
async def test_merge_p2p_shuffle(c, s, a, b, npartitions_left):
    df_left = pd.DataFrame({"a": [1, 2, 3] * 100, "b": 2})
    df_right = pd.DataFrame({"a": [4, 2, 3] * 100, "c": 2})
    left = from_pandas(df_left, npartitions=npartitions_left)
    right = from_pandas(df_right, npartitions=5)

    out = left.merge(right, shuffle_method="p2p")
    assert out.npartitions == npartitions_left
    x = c.compute(out)
    x = await x
    pd.testing.assert_frame_equal(x.reset_index(drop=True), df_left.merge(df_right))


@gen_cluster(client=True)
async def test_self_merge_p2p_shuffle(c, s, a, b):
    pdf = pd.DataFrame({"a": range(100), "b": range(0, 200, 2)})
    ddf = from_pandas(pdf, npartitions=5)

    out = ddf.merge(ddf, left_on="a", right_on="b", shuffle_method="p2p")
    # Generate unique shuffle IDs if the input frame is the same but parameters differ
    assert sum(id_from_key(k) is not None for k in out.dask) == 2
    x = await c.compute(out)
    expected = pdf.merge(pdf, left_on="a", right_on="b")
    pd.testing.assert_frame_equal(
        x.sort_values("a_x", ignore_index=True),
        expected.sort_values("a_x", ignore_index=True),
    )


@gen_cluster(client=True)
@pytest.mark.parametrize("name", ["a", None])
@pytest.mark.parametrize("shuffle", ["tasks", "disk", "p2p"])
async def test_merge_index_precedence(c, s, a, b, shuffle, name):
    pdf = pd.DataFrame(
        {"a": [1, 2, 3, 4, 5, 6]}, index=pd.Index([6, 5, 4, 3, 2, 1], name=name)
    )
    pdf2 = pd.DataFrame(
        {"b": [1, 2, 3, 4, 5, 6]}, index=pd.Index([1, 2, 7, 4, 5, 6], name=name)
    )
    df = from_pandas(pdf, npartitions=2, sort=False)
    df2 = from_pandas(pdf2, npartitions=3, sort=False)

    result = df.join(df2, shuffle_method=shuffle)
    x = await c.compute(result)
    assert result.npartitions == 3
    pd.testing.assert_frame_equal(x.sort_index(ascending=False), pdf.join(pdf2))


@pytest.mark.filterwarnings("ignore:DataFrameGroupBy.apply operated on")
@gen_cluster(client=True)
@pytest.mark.parametrize("shuffle", ["disk", "p2p", "tasks"])
async def test_shuffle_nulls_introduced(c, s, a, b, shuffle):
    df1 = pd.DataFrame([[True], [False]] * 50, columns=["A"])
    df1["B"] = list(range(100))

    df2 = pd.DataFrame(
        [[2, 3], [109, 2], [345, 3], [50, 7], [95, 1]], columns=["B", "C"]
    )

    ddf1 = from_pandas(df1, npartitions=10)
    ddf2 = from_pandas(df2, npartitions=1)
    meta = pd.Series(dtype=int, index=pd.Index([], dtype=bool, name="A"), name="A")
    result = (
        merge(ddf1, ddf2, how="outer", on="B", shuffle_method=shuffle)
        .groupby("A")
        .apply(lambda df: len(df), meta=meta, shuffle_method=shuffle)
    )
    expected = (
        pd.merge(df1, df2, how="outer", on="B").groupby("A").apply(lambda df: len(df))
    )
    x = await c.compute(result)

    pd.testing.assert_series_equal(
        x.sort_index(), expected.sort_index(), check_names=False
    )


@gen_cluster(client=True)
@pytest.mark.parametrize("shuffle", ["tasks", "p2p"])
@pytest.mark.parametrize("broadcast", [True, 0.6])
@pytest.mark.parametrize("how", ["left", "inner"])
async def test_merge_broadcast(c, s, a, b, shuffle, broadcast, how):
    pdf = pd.DataFrame({"a": [1, 2, 3, 4, 5, 6] * 5, "c": 1})
    pdf2 = pd.DataFrame({"b": [1, 2, 3, 4, 5, 6]})
    df = from_pandas(pdf, npartitions=15)
    df2 = from_pandas(pdf2, npartitions=2)

    result = df.merge(
        df2,
        left_on="a",
        right_on="b",
        shuffle_method=shuffle,
        broadcast=broadcast,
        how=how,
    )
    q = result.optimize()
    assert len(list(q.find_operations(BroadcastJoin))) > 0
    x = await c.compute(result)
    assert result.npartitions == 15
    pd.testing.assert_frame_equal(
        x.sort_values(by="a", ignore_index=True),
        pdf.merge(pdf2, left_on="a", right_on="b", how=how).sort_values(
            by="a", ignore_index=True
        ),
    )


@gen_cluster(client=True)
async def test_merge_p2p_shuffle_reused_dataframe_with_different_parameters(c, s, a, b):
    pdf1 = pd.DataFrame({"a": range(100), "b": range(0, 200, 2)})
    pdf2 = pd.DataFrame({"x": range(200), "y": [1, 2, 3, 4] * 50})
    ddf1 = from_pandas(pdf1, npartitions=5)
    ddf2 = from_pandas(pdf2, npartitions=10)

    out = (
        ddf1.merge(ddf2, left_on="a", right_on="x", shuffle_method="p2p")
        # Vary the number of output partitions for the shuffles of dd2
        .repartition(npartitions=20).merge(
            ddf2, left_on="b", right_on="x", shuffle_method="p2p"
        )
    )
    # Generate unique shuffle IDs if the input frame is the same but
    # parameters differ. Reusing shuffles in merges is dangerous because of the
    # required coordination and complexity introduced through dynamic clusters.
    assert sum(id_from_key(k) is not None for k in out.dask) == 4
    x = await c.compute(out)
    expected = pdf1.merge(pdf2, left_on="a", right_on="x").merge(
        pdf2, left_on="b", right_on="x"
    )
    pd.testing.assert_frame_equal(
        x.sort_values("a", ignore_index=True),
        expected.sort_values("a", ignore_index=True),
    )


@gen_cluster(client=True)
async def test_merge_p2p_shuffle_reused_dataframe_with_same_parameters(c, s, a, b):
    pdf1 = pd.DataFrame({"a": range(100), "b": range(0, 200, 2)})
    pdf2 = pd.DataFrame({"x": range(200), "y": [1, 2, 3, 4] * 50})
    ddf1 = from_pandas(pdf1, npartitions=5)
    ddf2 = from_pandas(pdf2, npartitions=10)

    # This performs two shuffles:
    #   * ddf1 is shuffled on `a`
    #   * ddf2 is shuffled on `x`
    ddf3 = ddf1.merge(ddf2, left_on="a", right_on="x", shuffle_method="p2p")

    # This performs one shuffle:
    #   * ddf3 is shuffled on `b`
    # We can reuse the shuffle of dd2 on `x` from the previous merge.
    out = ddf2.merge(
        ddf3,
        left_on="x",
        right_on="b",
        shuffle_method="p2p",
    )
    # Generate unique shuffle IDs if the input frame is the same and all its
    # parameters match. Reusing shuffles in merges is dangerous because of the
    # required coordination and complexity introduced through dynamic clusters.
    assert sum(id_from_key(k) is not None for k in out.dask) == 4
    x = await c.compute(out)
    expected = pdf2.merge(
        pdf1.merge(pdf2, left_on="a", right_on="x"), left_on="x", right_on="b"
    )
    pd.testing.assert_frame_equal(
        x.sort_values("a", ignore_index=True),
        expected.sort_values("a", ignore_index=True),
    )


@pytest.mark.parametrize("npartitions_left", [5, 6])
@gen_cluster(client=True)
async def test_index_merge_p2p_shuffle(c, s, a, b, npartitions_left):
    df_left = pd.DataFrame({"a": [1, 2, 3] * 100, "b": 2}).set_index("a")
    df_right = pd.DataFrame({"a": [4, 2, 3] * 100, "c": 2})
    left = from_pandas(df_left, npartitions=npartitions_left, sort=False)
    right = from_pandas(df_right, npartitions=5)

    out = left.merge(right, left_index=True, right_on="a", shuffle_method="p2p")
    assert out.npartitions == npartitions_left
    x = c.compute(out)
    x = await x
    pd.testing.assert_frame_equal(
        x.sort_index(),
        df_left.merge(df_right, left_index=True, right_on="a").sort_index(),
    )


@gen_cluster(client=True)
async def test_merge_p2p_shuffle(c, s, a, b):
    df_left = pd.DataFrame({"a": [1, 2, 3] * 100, "b": 2, "e": 2})
    df_right = pd.DataFrame({"a": [4, 2, 3] * 100, "c": 2})
    left = from_pandas(df_left, npartitions=6)
    right = from_pandas(df_right, npartitions=5)

    out = left.merge(right, shuffle_method="p2p")[["b", "c"]]
    assert out.npartitions == 6
    x = c.compute(out)
    x = await x
    pd.testing.assert_frame_equal(
        x.reset_index(drop=True), df_left.merge(df_right)[["b", "c"]]
    )


@gen_cluster(client=True)
async def test_merge_p2p_shuffle_projection_error(c, s, a, b):
    pdf1 = pd.DataFrame({"a": [1, 2, 3], "b": 1})
    pdf2 = pd.DataFrame({"x": [1, 2, 3, 4, 5, 6], "y": 1})
    df1 = from_pandas(pdf1, npartitions=2)
    df2 = from_pandas(pdf2, npartitions=3)
    df = df1.merge(df2, left_on="a", right_on="x")
    min_val = df.groupby("x")["y"].sum().reset_index()
    result = df.merge(min_val)
    expected = pd.DataFrame(
        {"a": [2, 3, 1], "b": 1, "x": [2, 3, 1], "y": 1}, index=[0, 0, 1]
    )
    x = c.compute(result)
    x = await x
    pd.testing.assert_frame_equal(
        x.sort_values("a", ignore_index=True),
        expected.sort_values("a", ignore_index=True),
    )


def test_sort_values():
    with LocalCluster(processes=False, n_workers=2, dashboard_address=":0") as cluster:
        with Client(cluster) as client:  # noqa: F841
            pdf = pd.DataFrame({"a": [5] + list(range(100)), "b": 2})
            df = from_pandas(pdf, npartitions=10)

            out = df.sort_values(by="a").compute()
    pd.testing.assert_frame_equal(
        out.reset_index(drop=True),
        pdf.sort_values(by="a", ignore_index=True),
    )


@pytest.mark.parametrize("add_repartition", [True, False])
def test_merge_combine_similar_squash_merges(add_repartition):
    with LocalCluster(processes=False, n_workers=2, dashboard_address=":0") as cluster:
        with Client(cluster) as client:  # noqa: F841
            pdf = pd.DataFrame(
                {
                    "a": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10] * 2,
                    "b": 1,
                    "c": 1,
                }
            )
            pdf2 = pd.DataFrame(
                {"m": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10] * 2, "n": 1, "o": 2, "p": 3}
            )
            pdf3 = pd.DataFrame(
                {"x": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10] * 2, "y": 1, "z": 1, "zz": 2}
            )

            df = from_pandas(pdf, npartitions=5)
            df2 = from_pandas(pdf2, npartitions=10)
            df3 = from_pandas(pdf3, npartitions=10)

            df = df[df.a > 1]
            df2 = df2[df2.m > 1]
            df3 = df3[df3.x > 1]
            if add_repartition:
                df = df.repartition(npartitions=df.npartitions // 2)
                df2 = df2.repartition(npartitions=df2.npartitions // 2)
            q = df.merge(df2, left_on="a", right_on="m")
            if add_repartition:
                df3 = df3.repartition(npartitions=df3.npartitions // 2)
                q = q.repartition(npartitions=q.npartitions // 2)
            q = q.merge(df3, left_on="n", right_on="x")
            q["revenue"] = q.y * (1 - q.z)
            result = q[["x", "n", "o", "revenue"]]
            result_q = result.optimize(fuse=False)

            assert (
                result_q.expr.frame.frame.frame._name
                == result_q.expr.frame.operands[2].left.frame._name
            )
            out = result.compute()

    pdf = pdf[pdf.a > 1]
    pdf2 = pdf2[pdf2.m > 1]
    pdf3 = pdf3[pdf3.x > 1]
    q = pdf.merge(pdf2, left_on="a", right_on="m")
    q = q.merge(pdf3, left_on="n", right_on="x")
    q["revenue"] = q.y * (1 - q.z)
    expected = q[["x", "n", "o", "revenue"]]

    pd.testing.assert_frame_equal(
        out.reset_index(drop=True),
        expected,
    )


@gen_cluster(client=True)
async def test_future_in_map_partitions(c, s, a, b):
    # xgboost uses this pattern

    def test_func(n):
        import pandas as pd

        return pd.DataFrame({"a": list(range(n))})

    df = from_pandas(pd.DataFrame({"a": [1, 2, 3, 4]}), npartitions=2)

    f = c.submit(test_func, 100)
    q = map_partitions(lambda x, y: y + x.sum(), f, df, meta=df._meta)
    result = c.compute(q)
    result = await result
    expected = pd.DataFrame({"a": [4951, 4952, 4953, 4954]})
    pd.testing.assert_frame_equal(result, expected)


@gen_cluster(client=True)
async def test_shuffle_consistency_checks(c, s, a, b):
    pdf = pd.DataFrame({"x": [1, 2, 3]})
    df = from_pandas(pdf, npartitions=2)
    df2 = df.astype(np.csingle)
    with pytest.raises(TypeError, match="p2p does not support data of type"):
        df2.shuffle("x")

    df.columns = [1]
    with pytest.raises(TypeError, match="p2p requires all column names to be str"):
        df.shuffle(1)


@gen_cluster(client=True)
async def test_merge_indicator(c, s, a, b):
    data = {
        "id": ["101-a", "102-a", "103-a"],
        "test": ["val101a", "val102a", "val103a"],
    }
    pdf = pd.DataFrame(data, dtype="string[pyarrow]")
    df = from_pandas(pdf, npartitions=2)
    result = df.merge(df, on="id", how="outer", indicator=True)
    x = c.compute(result)
    x = await x
    expected = pdf.merge(pdf, on="id", how="outer", indicator=True)

    pd.testing.assert_frame_equal(
        x.sort_values("id", ignore_index=True),
        expected.sort_values("id", ignore_index=True),
    )


@gen_cluster(client=True)
async def test_shuffle_partition_reduction(c, s, a, b):
    pdf = pd.DataFrame({"a": [1, 2, 3, 4] * 100, "b": 1})
    df = from_pandas(pdf, npartitions=10)
    result = df.shuffle(on="a", npartitions=4)
    q = result.optimize(fuse=False)
    assert not any(isinstance(op, Repartition) for op in q.walk())
    x = c.compute(result)
    x = await x

    pd.testing.assert_frame_equal(
        x.sort_values("a", ignore_index=True),
        pdf.sort_values("a", ignore_index=True),
    )


@gen_cluster(client=True)
async def test_p2p_and_merge_shuffle(c, s, a, b):
    pdf = pd.DataFrame({"a": np.random.randint(1, 100, (100,)), "b": 1})

    df = from_pandas(pdf, npartitions=15)

    pdf2 = pd.DataFrame({"a": np.random.randint(1, 100, (100,)), "c": 1})
    df2 = from_pandas(pdf2, npartitions=10)

    pdf3 = pd.DataFrame({"a": np.random.randint(1, 100, (100,)), "d": 2})
    df3 = from_pandas(pdf3, npartitions=5)

    result = df.merge(df2)
    result = result.merge(df3)
    x = c.compute(result)
    x = await x
    pd.testing.assert_frame_equal(
        x.sort_values("a", ignore_index=True),
        pdf.merge(pdf2).merge(pdf3).sort_values("a", ignore_index=True),
    )


@pytest.mark.parametrize("func", ["set_index", "sort_values", "shuffle"])
def test_respect_context_shuffle(df, pdf, func):
    with dask.config.set({"dataframe.shuffle.method": "tasks"}):
        q = getattr(df, func)("x")
    result = q.optimize(fuse=False)
    assert len([x for x in result.walk() if isinstance(x, TaskShuffle)]) > 0

    q = getattr(df, func)("x")
    with dask.config.set({"dataframe.shuffle.method": "tasks"}):
        result = q.optimize(fuse=False)
    assert len([x for x in result.walk() if isinstance(x, TaskShuffle)]) > 0

    with dask.config.set({"dataframe.shuffle.method": "p2p"}):
        q = getattr(df, func)("x")
    with dask.config.set({"dataframe.shuffle.method": "tasks"}):
        result = q.optimize(fuse=False)
    assert len([x for x in result.walk() if isinstance(x, P2PShuffle)]) > 0


@pytest.mark.parametrize("concatenate", [True, False])
def test_compute_concatenates(loop, concatenate):
    pdf = pd.DataFrame({"a": np.random.randint(1, 100, (100,)), "b": 1})
    df = from_pandas(pdf, npartitions=10)

    class Plugin(SchedulerPlugin):
        def start(self, *args, **kwargs):
            self.repartition_in_tasks = False

        def update_graph(
            self,
            scheduler,
            *,
            client,
            keys,
            tasks,
            annotations,
            priority,
            dependencies,
            **kwargs,
        ):
            for key in dependencies:
                if not isinstance(key, tuple):
                    continue
                group = key[0]
                if not isinstance(group, str):
                    continue
                self.repartition_in_tasks |= group.startswith("repartitiontofewer")

    with Client(loop=loop) as c:
        c.register_plugin(Plugin(), name="tracker")
        df.compute(fuse=False, concatenate=concatenate)
        plugin = c.cluster.scheduler.plugins["tracker"]
        assert plugin.repartition_in_tasks is concatenate
