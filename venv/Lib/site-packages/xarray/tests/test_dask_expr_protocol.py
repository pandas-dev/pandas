from __future__ import annotations

from importlib import import_module
from typing import Any

import numpy as np
import pytest

import xarray as xr
from xarray import DataArray, Dataset
from xarray.testing import assert_identical
from xarray.tests import requires_scipy_or_netCDF4

dask = pytest.importorskip("dask")
da = pytest.importorskip("dask.array")
dask_array = pytest.importorskip("dask_array")

dask_expr = pytest.importorskip("dask._expr")
CompositeExpr: Any = getattr(dask_expr, "CompositeExpr", None)
HLGExpr: Any = getattr(dask_expr, "HLGExpr", None)
if CompositeExpr is None or HLGExpr is None:
    pytest.skip("requires Dask composite expressions", allow_module_level=True)


def test_dataset_composite_expr_protocol_simple():
    x = dask_array.arange(3, chunks=(1,))
    ds = Dataset({"x": ("i", x + 1)})

    expr = dask.base.collections_to_expr(ds)

    assert isinstance(expr, CompositeExpr)
    assert len(expr.exprs) == 1
    assert_identical(
        dask.compute(ds, scheduler="single-threaded")[0],
        Dataset({"x": ("i", np.arange(3) + 1)}),
    )


def test_variable_expr_protocol_avoids_graph_materialization(monkeypatch):
    x = dask_array.arange(6, chunks=(3,))
    var = xr.Variable(("i",), x + 1)

    def raise_if_materialized(self):
        raise AssertionError("__dask_graph__ should not be materialized")

    monkeypatch.setattr(dask_array.Array, "__dask_graph__", raise_if_materialized)

    assert dask.is_dask_collection(var)
    exprs = var.__dask_exprs__()
    assert exprs is not None
    assert len(exprs) == 1
    assert isinstance(dask.base.collections_to_expr(Dataset({"x": var})), CompositeExpr)


def test_dataarray_composite_expr_protocol_includes_chunked_coord():
    x = dask_array.arange(6, chunks=(3,))
    arr = DataArray(x + 1, dims=("i",), coords={"coord": ("i", x + 10)}, name="z")

    expr = dask.base.collections_to_expr(arr)

    assert isinstance(expr, CompositeExpr)
    assert len(expr.exprs) == 2
    assert_identical(
        dask.compute(arr, scheduler="single-threaded")[0],
        DataArray(
            np.arange(6) + 1,
            dims=("i",),
            coords={"coord": ("i", np.arange(6) + 10)},
            name="z",
        ),
    )


def test_dataset_compute_persist_optimize_end_to_end():
    x = dask_array.arange(6, chunks=(3,))
    ds = Dataset(
        {"foo": ("i", x + 1)},
        coords={"coord": ("i", x + 10)},
        attrs={"source": "test"},
    )
    ds["foo"].encoding["example"] = "kept"

    expected = Dataset(
        {"foo": ("i", np.arange(6) + 1)},
        coords={"coord": ("i", np.arange(6) + 10)},
        attrs={"source": "test"},
    )

    computed_ds, computed_x = dask.compute(ds, x, scheduler="single-threaded")
    assert_identical(computed_ds, expected)
    np.testing.assert_array_equal(computed_x, np.arange(6))
    assert computed_ds["foo"].encoding["example"] == "kept"

    persisted = dask.persist(ds, scheduler="single-threaded")[0]
    assert isinstance(persisted["foo"].data, dask_array.Array)
    assert_identical(persisted.compute(), expected)
    assert persisted["foo"].encoding["example"] == "kept"

    optimized = dask.optimize(ds)[0]
    assert isinstance(optimized["foo"].data, dask_array.Array)
    assert_identical(optimized.compute(), expected)


@requires_scipy_or_netCDF4
def test_open_mfdataset_end_to_end(tmp_path):
    paths = []
    for i in range(2):
        path = tmp_path / f"part-{i}.nc"
        Dataset(
            {"x": ("t", np.arange(i * 3, i * 3 + 3))},
            coords={"t": np.arange(i * 3, i * 3 + 3)},
        ).to_netcdf(path)
        paths.append(path)

    with xr.open_mfdataset(paths, chunks={"t": 2}, combine="by_coords") as ds:
        assert isinstance(ds["x"].data, dask_array.Array)
        assert isinstance(dask.base.collections_to_expr(ds), CompositeExpr)

        expected = Dataset({"x": ("t", np.arange(6))}, coords={"t": np.arange(6)})
        assert_identical(ds.compute(scheduler="single-threaded"), expected)
        assert_identical(dask.compute(ds, scheduler="single-threaded")[0], expected)
        assert_identical(
            ds.persist(scheduler="single-threaded").compute(
                scheduler="single-threaded"
            ),
            expected,
        )
        assert_identical(
            dask.optimize(ds)[0].compute(scheduler="single-threaded"), expected
        )

        mapped = xr.map_blocks(lambda block: block + 1, ds)
        assert isinstance(mapped["x"].data, dask_array.Array)
        assert isinstance(dask.base.collections_to_expr(mapped), CompositeExpr)
        assert_identical(
            mapped.compute(scheduler="single-threaded"),
            Dataset({"x": ("t", np.arange(6) + 1)}, coords={"t": np.arange(6)}),
        )


@requires_scipy_or_netCDF4
def test_open_dataset_rechunk_optimization_crosses_composite_expr(tmp_path):
    Rechunk = import_module("dask_array._rechunk").Rechunk

    path = tmp_path / "data.nc"
    Dataset({"x": ("t", np.arange(12))}, coords={"t": np.arange(12)}).to_netcdf(path)

    with xr.open_dataset(path, chunks={"t": 3}) as ds:
        out = ds.chunk({"t": 4})
        expr = dask.base.collections_to_expr(out)

        assert isinstance(expr, CompositeExpr)
        assert len(expr.exprs) == 1
        assert list(expr.exprs[0].find_operations(Rechunk))

        optimized_expr = expr.optimize()

        assert not list(optimized_expr.exprs[0].find_operations(Rechunk))

        optimized = dask.optimize(out)[0]
        assert isinstance(optimized["x"].data, dask_array.Array)
        assert not list(optimized["x"].data.expr.find_operations(Rechunk))
        assert_identical(
            optimized.compute(scheduler="single-threaded"),
            Dataset({"x": ("t", np.arange(12))}, coords={"t": np.arange(12)}),
        )


def test_map_blocks_dataarray_end_to_end():
    x = dask_array.arange(6, chunks=(3,))
    arr = DataArray(x, dims="t", name="x")
    out = xr.map_blocks(lambda block: block + 1, arr)

    assert isinstance(out.data, dask_array.Array)
    assert isinstance(dask.base.collections_to_expr(out), CompositeExpr)

    expected = DataArray(np.arange(6) + 1, dims="t", name="x")
    assert_identical(out.compute(scheduler="single-threaded"), expected)
    assert_identical(
        out.persist(scheduler="single-threaded").compute(scheduler="single-threaded"),
        expected,
    )
    assert_identical(
        dask.optimize(out)[0].compute(scheduler="single-threaded"), expected
    )


def test_map_blocks_dataset_outputs_share_block_calls():
    calls = []
    x = dask_array.arange(6, chunks=(3,))
    ds = Dataset({"x": ("t", x)}, coords={"qc": ("t", x + 10)})
    template = Dataset(
        {"a": ("t", x), "b": ("t", x)},
        coords={"qc": ("t", x + 10)},
        attrs={"kind": "mapped"},
    )

    def func(block):
        calls.append(block.sizes["t"])
        return Dataset(
            {"a": block["x"] + 1, "b": block["x"] + 2},
            coords={"qc": block["qc"]},
            attrs={"kind": "mapped"},
        )

    out = xr.map_blocks(func, ds, template=template)

    assert isinstance(dask.base.collections_to_expr(out), CompositeExpr)
    assert all(isinstance(out[name].data, dask_array.Array) for name in out)
    assert sorted(calls) == []

    expected = Dataset(
        {"a": ("t", np.arange(6) + 1), "b": ("t", np.arange(6) + 2)},
        coords={"qc": ("t", np.arange(6) + 10)},
        attrs={"kind": "mapped"},
    )
    assert_identical(out.compute(scheduler="single-threaded"), expected)
    assert sorted(calls) == [3, 3]


def test_map_blocks_preserves_scalar_coords():
    x = dask_array.arange(6, chunks=(3,)).reshape((3, 2))
    arr = DataArray(
        x,
        dims=("x", "y"),
        coords={"label": ("x", ["a", "b", "c"]), "scale": 2},
        name="z",
    )

    out = xr.map_blocks(lambda block: block + block.scale, arr)

    assert isinstance(dask.base.collections_to_expr(out), CompositeExpr)
    assert_identical(
        out.compute(scheduler="single-threaded"),
        DataArray(
            np.arange(6).reshape((3, 2)) + 2,
            dims=("x", "y"),
            coords={"label": ("x", ["a", "b", "c"]), "scale": 2},
        ),
    )


def test_map_blocks_reduces_single_chunk_dimension():
    x = dask_array.arange(12, chunks=(12,)).reshape((3, 4)).rechunk((3, 2))
    arr = DataArray(x, dims=("x", "y"), name="z")

    out = xr.map_blocks(lambda block: block.sum("x"), arr)

    assert isinstance(dask.base.collections_to_expr(out), CompositeExpr)
    assert_identical(
        out.compute(scheduler="single-threaded"),
        DataArray(np.arange(12).reshape(3, 4).sum(axis=0), dims="y", name="z"),
    )


def test_map_blocks_slice_pushdown_equivalent_to_preselection():
    x = dask_array.arange(12, chunks=(3,))
    ds = Dataset({"x": ("t", x)}, coords={"t": np.arange(12)})
    post_calls: list[tuple[int, ...]] = []
    pre_calls: list[tuple[int, ...]] = []

    def add_one(calls):
        def func(block):
            if block.sizes["t"]:
                calls.append(tuple(block["t"].values.tolist()))
            return block + 1

        return func

    post_selected = xr.map_blocks(add_one(post_calls), ds).sel(t=slice(3, 8))
    pre_selected = xr.map_blocks(add_one(pre_calls), ds.sel(t=slice(3, 8)))

    assert isinstance(dask.base.collections_to_expr(post_selected), CompositeExpr)
    assert_identical(
        post_selected.compute(scheduler="single-threaded"),
        pre_selected.compute(scheduler="single-threaded"),
    )
    assert sorted(post_calls) == [(3, 4, 5), (6, 7, 8)]
    assert sorted(pre_calls) == [(3, 4, 5), (6, 7, 8)]


def test_mixed_legacy_inputs_do_not_use_composite_path():
    ds = Dataset(
        {
            "x": ("i", dask_array.arange(3, chunks=(1,))),
            "legacy": ("i", da.arange(3, chunks=(1,))),
        }
    )
    expected = Dataset(
        {"x": ("i", np.arange(3)), "legacy": ("i", np.arange(3))},
    )

    assert ds.__dask_exprs__() is None
    assert isinstance(dask.base.collections_to_expr(ds), HLGExpr)
    assert_identical(dask.compute(ds, scheduler="single-threaded")[0], expected)

    persisted = dask.persist(ds, scheduler="single-threaded")[0]
    assert isinstance(persisted["x"].data, dask_array.Array)
    assert isinstance(persisted["legacy"].data, da.Array)
    assert_identical(
        dask.compute(persisted, scheduler="single-threaded")[0],
        expected,
    )

    optimized = dask.optimize(ds)[0]
    assert isinstance(optimized["x"].data, dask_array.Array)
    assert isinstance(optimized["legacy"].data, da.Array)
    assert_identical(
        dask.compute(optimized, scheduler="single-threaded")[0],
        expected,
    )

    with pytest.raises(TypeError, match=r"cannot mix dask_array\.Array"):
        xr.map_blocks(lambda block: block + 1, ds)

    arr = DataArray(dask_array.arange(6, chunks=(3,)), dims="i")
    other = DataArray(da.arange(6, chunks=(3,)), dims="i")
    with pytest.raises(TypeError, match=r"cannot mix dask_array\.Array"):
        xr.map_blocks(lambda a, b: a + b, arr, args=[other])


def test_shared_subexpressions_optimize_without_cross_contamination():
    from dask.core import flatten

    x = dask_array.arange(6, chunks=(3,))
    ds = Dataset({"foo": ("i", x + 1), "bar": ("i", x + 2)})

    optimized = dask.optimize(ds)[0]

    foo_graph_keys = set(optimized["foo"].data.__dask_graph__())
    bar_output_keys = set(flatten(optimized["bar"].data.__dask_keys__()))
    assert not foo_graph_keys & bar_output_keys
    assert_identical(
        optimized.compute(),
        Dataset({"foo": ("i", np.arange(6) + 1), "bar": ("i", np.arange(6) + 2)}),
    )


def test_rechunk_reduction_chain_uses_composite_expr():
    x = dask_array.arange(12, chunks=(4,)).reshape((3, 4))
    out = Dataset({"x": (("a", "b"), x)}).chunk({"a": 1, "b": 2}).sum("b") + 1

    assert isinstance(dask.base.collections_to_expr(out), CompositeExpr)
    expected = Dataset({"x": ("a", np.arange(12).reshape(3, 4).sum(axis=1) + 1)})
    assert_identical(out.compute(scheduler="single-threaded"), expected)
    assert_identical(
        dask.optimize(out)[0].compute(scheduler="single-threaded"), expected
    )


def test_apply_ufunc_parallelized_uses_composite_expr():
    x = dask_array.arange(6, chunks=(3,))
    arr = DataArray(x, dims="t", name="x")
    out = xr.apply_ufunc(lambda z: z + 2, arr, dask="parallelized", output_dtypes=[int])

    assert isinstance(dask.base.collections_to_expr(out), CompositeExpr)
    assert_identical(
        out.compute(scheduler="single-threaded"),
        DataArray(np.arange(6) + 2, dims="t", name="x"),
    )


@pytest.mark.xfail_with_dask_array(
    reason="flox groupby currently builds legacy dask arrays"
)
def test_groupby_sum_uses_composite_expr():
    x = dask_array.arange(6, chunks=(3,))
    arr = DataArray(
        x,
        dims="t",
        coords={"label": ("t", np.array(["a", "a", "b", "b", "a", "b"]))},
        name="x",
    )

    out = arr.groupby("label").sum()

    assert isinstance(out.data, dask_array.Array)
    assert isinstance(dask.base.collections_to_expr(out), CompositeExpr)
    assert_identical(
        out.compute(scheduler="single-threaded"),
        DataArray([5, 10], dims="label", coords={"label": ["a", "b"]}, name="x"),
    )
