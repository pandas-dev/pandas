import numpy as np
import pytest

import xarray as xr
from xarray.indexes import NDPointIndex
from xarray.tests import assert_identical

pytest.importorskip("scipy")


def test_tree_index_init() -> None:
    from xarray.indexes.nd_point_index import ScipyKDTreeAdapter

    xx, yy = np.meshgrid([1.0, 2.0], [3.0, 4.0])
    ds = xr.Dataset(coords={"xx": (("y", "x"), xx), "yy": (("y", "x"), yy)})

    ds_indexed1 = ds.set_xindex(("xx", "yy"), NDPointIndex)
    assert "xx" in ds_indexed1.xindexes
    assert "yy" in ds_indexed1.xindexes
    assert isinstance(ds_indexed1.xindexes["xx"], NDPointIndex)
    assert ds_indexed1.xindexes["xx"] is ds_indexed1.xindexes["yy"]

    ds_indexed2 = ds.set_xindex(
        ("xx", "yy"), NDPointIndex, tree_adapter_cls=ScipyKDTreeAdapter
    )
    assert ds_indexed1.xindexes["xx"].equals(ds_indexed2.xindexes["yy"])


def test_tree_index_init_errors() -> None:
    xx, yy = np.meshgrid([1.0, 2.0], [3.0, 4.0])
    ds = xr.Dataset(coords={"xx": (("y", "x"), xx), "yy": (("y", "x"), yy)})

    with pytest.raises(ValueError, match="number of variables"):
        ds.set_xindex("xx", NDPointIndex)

    ds2 = ds.assign_coords(yy=(("u", "v"), [[3.0, 3.0], [4.0, 4.0]]))

    with pytest.raises(ValueError, match="same dimensions"):
        ds2.set_xindex(("xx", "yy"), NDPointIndex)


def test_tree_index_sel() -> None:
    xx, yy = np.meshgrid([1.0, 2.0], [3.0, 4.0])
    ds = xr.Dataset(coords={"xx": (("y", "x"), xx), "yy": (("y", "x"), yy)}).set_xindex(
        ("xx", "yy"), NDPointIndex
    )

    # 1-dimensional labels
    actual = ds.sel(
        xx=xr.Variable("u", [1.1, 1.1, 1.1]),
        yy=xr.Variable("u", [3.1, 3.1, 3.1]),
        method="nearest",
    )
    expected = xr.Dataset(
        coords={"xx": ("u", [1.0, 1.0, 1.0]), "yy": ("u", [3.0, 3.0, 3.0])}
    )
    assert_identical(actual, expected)

    # 2-dimensional labels
    actual = ds.sel(
        xx=xr.Variable(("u", "v"), [[1.1, 1.1, 1.1], [1.9, 1.9, 1.9]]),
        yy=xr.Variable(("u", "v"), [[3.1, 3.1, 3.1], [3.9, 3.9, 3.9]]),
        method="nearest",
    )
    expected = xr.Dataset(
        coords={
            "xx": (("u", "v"), [[1.0, 1.0, 1.0], [2.0, 2.0, 2.0]]),
            "yy": (("u", "v"), [[3.0, 3.0, 3.0], [4.0, 4.0, 4.0]]),
        },
    )
    assert_identical(actual, expected)

    # all scalar labels
    actual = ds.sel(xx=1.1, yy=3.1, method="nearest")
    expected = xr.Dataset(coords={"xx": 1.0, "yy": 3.0})
    assert_identical(actual, expected)

    # broadcast scalar to label shape and dimensions
    actual = ds.sel(xx=1.1, yy=xr.Variable("u", [3.1, 3.1, 3.1]), method="nearest")
    expected = ds.sel(
        xx=xr.Variable("u", [1.1, 1.1, 1.1]),
        yy=xr.Variable("u", [3.1, 3.1, 3.1]),
        method="nearest",
    )
    assert_identical(actual, expected)

    # broadcast orthogonal 1-dimensional labels
    actual = ds.sel(
        xx=xr.Variable("u", [1.1, 1.1]),
        yy=xr.Variable("v", [3.1, 3.1]),
        method="nearest",
    )
    expected = xr.Dataset(
        coords={
            "xx": (("u", "v"), [[1.0, 1.0], [1.0, 1.0]]),
            "yy": (("u", "v"), [[3.0, 3.0], [3.0, 3.0]]),
        },
    )
    assert_identical(actual, expected)

    # implicit dimension array-like labels
    actual = ds.sel(
        xx=[[1.1, 1.1, 1.1], [1.9, 1.9, 1.9]],
        yy=[[3.1, 3.1, 3.1], [3.9, 3.9, 3.9]],
        method="nearest",
    )
    expected = ds.sel(
        xx=xr.Variable(ds.xx.dims, [[1.1, 1.1, 1.1], [1.9, 1.9, 1.9]]),
        yy=xr.Variable(ds.yy.dims, [[3.1, 3.1, 3.1], [3.9, 3.9, 3.9]]),
        method="nearest",
    )
    assert_identical(actual, expected)


def test_tree_index_sel_errors() -> None:
    xx, yy = np.meshgrid([1.0, 2.0], [3.0, 4.0])
    ds = xr.Dataset(coords={"xx": (("y", "x"), xx), "yy": (("y", "x"), yy)}).set_xindex(
        ("xx", "yy"), NDPointIndex
    )

    with pytest.raises(ValueError, match="method='nearest'"):
        ds.sel(xx=1.1, yy=3.1)

    with pytest.raises(ValueError, match="missing labels"):
        ds.sel(xx=1.1, method="nearest")

    with pytest.raises(ValueError, match="invalid label value"):
        # invalid array-like dimensions
        ds.sel(xx=[1.1, 1.9], yy=[3.1, 3.9], method="nearest")

    # error while trying to broadcast labels
    with pytest.raises(xr.AlignmentError, match=".*conflicting dimension sizes"):
        ds.sel(
            xx=xr.Variable("u", [1.1, 1.1, 1.1]),
            yy=xr.Variable("u", [3.1, 3.1]),
            method="nearest",
        )


def test_tree_index_equals() -> None:
    xx1, yy1 = np.meshgrid([1.0, 2.0], [3.0, 4.0])
    ds1 = xr.Dataset(
        coords={"xx": (("y", "x"), xx1), "yy": (("y", "x"), yy1)}
    ).set_xindex(("xx", "yy"), NDPointIndex)

    xx2, yy2 = np.meshgrid([1.0, 2.0], [3.0, 4.0])
    ds2 = xr.Dataset(
        coords={"xx": (("y", "x"), xx2), "yy": (("y", "x"), yy2)}
    ).set_xindex(("xx", "yy"), NDPointIndex)

    xx3, yy3 = np.meshgrid([10.0, 20.0], [30.0, 40.0])
    ds3 = xr.Dataset(
        coords={"xx": (("y", "x"), xx3), "yy": (("y", "x"), yy3)}
    ).set_xindex(("xx", "yy"), NDPointIndex)

    assert ds1.xindexes["xx"].equals(ds2.xindexes["xx"])
    assert not ds1.xindexes["xx"].equals(ds3.xindexes["xx"])


def test_tree_index_rename() -> None:
    xx, yy = np.meshgrid([1.0, 2.0], [3.0, 4.0])
    ds = xr.Dataset(coords={"xx": (("y", "x"), xx), "yy": (("y", "x"), yy)}).set_xindex(
        ("xx", "yy"), NDPointIndex
    )

    ds_renamed = ds.rename_dims(y="u").rename_vars(yy="uu")
    assert "uu" in ds_renamed.xindexes
    assert isinstance(ds_renamed.xindexes["uu"], NDPointIndex)
    assert ds_renamed.xindexes["xx"] is ds_renamed.xindexes["uu"]

    # test via sel() with implicit dimension array-like labels, which relies on
    # NDPointIndex._coord_names and NDPointIndex._dims internal attrs
    actual = ds_renamed.sel(
        xx=[[1.1, 1.1, 1.1], [1.9, 1.9, 1.9]],
        uu=[[3.1, 3.1, 3.1], [3.9, 3.9, 3.9]],
        method="nearest",
    )
    expected = ds_renamed.sel(
        xx=xr.Variable(ds_renamed.xx.dims, [[1.1, 1.1, 1.1], [1.9, 1.9, 1.9]]),
        uu=xr.Variable(ds_renamed.uu.dims, [[3.1, 3.1, 3.1], [3.9, 3.9, 3.9]]),
        method="nearest",
    )
    assert_identical(actual, expected)
