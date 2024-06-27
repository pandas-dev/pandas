from __future__ import annotations

import os

import pytest

pytest.importorskip("distributed")

import sys
from operator import getitem

from distributed import Client, SchedulerPlugin
from distributed.utils_test import cluster, loop  # noqa F401

from dask.highlevelgraph import HighLevelGraph
from dask.layers import ArrayChunkShapeDep, ArraySliceDep, fractional_slice


class SchedulerImportCheck(SchedulerPlugin):
    """Plugin to help record which modules are imported on the scheduler"""

    name = "import-check"

    def __init__(self, pattern):
        self.pattern = pattern

    async def start(self, scheduler):
        # Record the modules that have been imported when the scheduler starts
        self.start_modules = set()
        for mod in set(sys.modules):
            if not mod.startswith(self.pattern):
                self.start_modules.add(mod)
            else:
                # Manually remove the target library
                sys.modules.pop(mod)


def test_array_chunk_shape_dep():
    pytest.importorskip("numpy")
    dac = pytest.importorskip("dask.array.core")
    d = 2  # number of chunks in x,y
    chunk = (2, 3)  # chunk shape
    shape = tuple(d * n for n in chunk)  # array shape
    chunks = dac.normalize_chunks(chunk, shape)
    array_deps = ArrayChunkShapeDep(chunks)

    def check(i, j):
        chunk_shape = array_deps[(i, j)]
        assert chunk_shape == chunk

    for i in range(d):
        for j in range(d):
            check(i, j)


def test_array_slice_deps():
    pytest.importorskip("numpy")
    dac = pytest.importorskip("dask.array.core")
    d = 2  # number of chunks in x,y
    chunk = (2, 3)  # chunk shape
    shape = tuple(d * n for n in chunk)  # array shape
    chunks = dac.normalize_chunks(chunk, shape)
    array_deps = ArraySliceDep(chunks)

    def check(i, j):
        slices = array_deps[(i, j)]
        assert slices == (
            slice(chunk[0] * i, chunk[0] * (i + 1), None),
            slice(chunk[1] * j, chunk[1] * (j + 1), None),
        )

    for i in range(d):
        for j in range(d):
            check(i, j)


def _dataframe_shuffle(tmpdir):
    pd = pytest.importorskip("pandas")
    dd = pytest.importorskip("dask.dataframe")

    # Perform a computation using an HLG-based shuffle
    df = pd.DataFrame({"a": range(10), "b": range(10, 20)})
    return dd.from_pandas(df, npartitions=2).shuffle("a", shuffle_method="tasks")


def _dataframe_tree_reduction(tmpdir):
    pd = pytest.importorskip("pandas")
    dd = pytest.importorskip("dask.dataframe")

    # Perform a computation using an HLG-based tree reduction
    df = pd.DataFrame({"a": range(10), "b": range(10, 20)})
    return dd.from_pandas(df, npartitions=2).mean()


def _dataframe_broadcast_join(tmpdir):
    pd = pytest.importorskip("pandas")
    dd = pytest.importorskip("dask.dataframe")

    # Perform a computation using an HLG-based broadcast join
    df = pd.DataFrame({"a": range(10), "b": range(10, 20)})
    ddf1 = dd.from_pandas(df, npartitions=4)
    ddf2 = dd.from_pandas(df, npartitions=1)
    return ddf1.merge(ddf2, how="left", broadcast=True, shuffle_method="tasks")


def _array_creation(tmpdir):
    da = pytest.importorskip("dask.array")

    # Perform a computation using HLG-based array creation
    return da.ones((100,)) + da.zeros((100,))


def _array_map_overlap(tmpdir):
    da = pytest.importorskip("dask.array")
    array = da.ones((100,))
    return array.map_overlap(lambda x: x, depth=1, boundary="none")


def test_fractional_slice():
    assert fractional_slice(("x", 4.9), {0: 2}) == (getitem, ("x", 5), (slice(0, 2),))

    assert fractional_slice(("x", 3, 5.1), {0: 2, 1: 3}) == (
        getitem,
        ("x", 3, 5),
        (slice(None, None, None), slice(-3, None)),
    )

    assert fractional_slice(("x", 2.9, 5.1), {0: 2, 1: 3}) == (
        getitem,
        ("x", 3, 5),
        (slice(0, 2), slice(-3, None)),
    )

    fs = fractional_slice(("x", 4.9), {0: 2})
    assert isinstance(fs[1][1], int)


def _pq_pyarrow(tmpdir):
    pytest.importorskip("pyarrow")
    pd = pytest.importorskip("pandas")
    dd = pytest.importorskip("dask.dataframe")

    dd.from_pandas(pd.DataFrame({"a": range(10)}), npartitions=2).to_parquet(
        str(tmpdir)
    )
    filters = [(("a", "<=", 2))]

    ddf1 = dd.read_parquet(str(tmpdir), filters=filters)
    return ddf1


def _pq_fastparquet(tmpdir):
    pytest.importorskip("fastparquet")
    pd = pytest.importorskip("pandas")
    dd = pytest.importorskip("dask.dataframe")

    df = dd.from_pandas(pd.DataFrame({"a": range(10)}), npartitions=2)
    with pytest.warns(FutureWarning):
        df.to_parquet(str(tmpdir), engine="fastparquet")
        return dd.read_parquet(str(tmpdir), engine="fastparquet")


def _read_csv(tmpdir):
    pd = pytest.importorskip("pandas")
    dd = pytest.importorskip("dask.dataframe")

    dd.from_pandas(
        pd.DataFrame({"a": range(10)}),
        npartitions=2,
    ).to_csv(str(tmpdir))
    return dd.read_csv(os.path.join(str(tmpdir), "*"))


@pytest.mark.xfail(reason="#8480")
@pytest.mark.parametrize(
    "op,lib",
    [
        (_dataframe_shuffle, "pandas."),
        (_dataframe_tree_reduction, "pandas."),
        (_dataframe_broadcast_join, "pandas."),
        (_pq_pyarrow, "pandas."),
        (_pq_fastparquet, "pandas."),
        (_read_csv, "pandas."),
        (_array_creation, "numpy."),
        (_array_map_overlap, "numpy."),
    ],
)
@pytest.mark.parametrize("optimize_graph", [True, False])
def test_scheduler_highlevel_graph_unpack_import(op, lib, optimize_graph, loop, tmpdir):
    # Test that array/dataframe-specific modules are not imported
    # on the scheduler when an HLG layers are unpacked/materialized.

    with cluster(scheduler_kwargs={"plugins": [SchedulerImportCheck(lib)]}) as (
        scheduler,
        workers,
    ):
        with Client(scheduler["address"], loop=loop) as c:
            # Perform a computation using a HighLevelGraph Layer
            c.compute(op(tmpdir), optimize_graph=optimize_graph)

            # Get the new modules which were imported on the scheduler during the computation
            end_modules = c.run_on_scheduler(lambda: set(sys.modules))
            start_modules = c.run_on_scheduler(
                lambda dask_scheduler: dask_scheduler.plugins[
                    SchedulerImportCheck.name
                ].start_modules
            )
            new_modules = end_modules - start_modules

            # Check that the scheduler didn't start with `lib`
            # (otherwise we aren't testing anything)
            assert not any(module.startswith(lib) for module in start_modules)

            # Check whether we imported `lib` on the scheduler
            assert not any(module.startswith(lib) for module in new_modules)


def _shuffle_op(ddf):
    return ddf.shuffle("x", shuffle_method="tasks")


def _groupby_op(ddf):
    return ddf.groupby("name").agg({"x": "mean"})


@pytest.mark.parametrize("op", [_shuffle_op, _groupby_op])
def test_dataframe_cull_key_dependencies(op):
    # Test that HighLevelGraph.cull does not populate the
    # output graph with incorrect key_dependencies for
    # "complex" DataFrame Layers
    # See: https://github.com/dask/dask/pull/9267
    pytest.importorskip("pandas")
    dd = pytest.importorskip("dask.dataframe")
    if dd._dask_expr_enabled():
        pytest.skip("not supported")
    datasets = pytest.importorskip("dask.datasets")

    result = op(datasets.timeseries(end="2000-01-15")).count()
    graph = result.dask
    culled_graph = graph.cull(result.__dask_keys__())

    assert graph.get_all_dependencies() == culled_graph.get_all_dependencies()


def test_dataframe_cull_key_dependencies_materialized():
    # Test that caching of MaterializedLayer
    # dependencies during culling doesn't break
    # the result of ``get_all_dependencies``
    pytest.importorskip("pandas")
    datasets = pytest.importorskip("dask.datasets")
    dd = pytest.importorskip("dask.dataframe")
    if dd._dask_expr_enabled():
        pytest.skip("not supported")

    ddf = datasets.timeseries(end="2000-01-15")

    # Build a custom layer to ensure
    # MaterializedLayer is used
    name = "custom_graph_test"
    name_0 = "custom_graph_test_0"
    dsk = {}
    for i in range(ddf.npartitions):
        dsk[(name_0, i)] = (lambda x: x, (ddf._name, i))
        dsk[(name, i)] = (lambda x: x, (name_0, i))
    dsk = HighLevelGraph.from_collections(name, dsk, dependencies=[ddf])
    result = dd.core.new_dd_object(dsk, name, ddf._meta, ddf.divisions)
    graph = result.dask

    # HLG cull
    culled_keys = [k for k in result.__dask_keys__() if k != (name, 0)]
    culled_graph = graph.cull(culled_keys)

    # Check that culled_deps are cached
    # See: https://github.com/dask/dask/issues/9389
    cached_deps = culled_graph.key_dependencies.copy()
    deps = culled_graph.get_all_dependencies()
    assert cached_deps == deps

    # Manual cull
    deps0 = graph.get_all_dependencies()
    for name, i in list(deps0.keys()):
        if i == 0:
            deps0.pop((name, i))

    # Check that get_all_dependencies results match
    assert deps0 == deps
