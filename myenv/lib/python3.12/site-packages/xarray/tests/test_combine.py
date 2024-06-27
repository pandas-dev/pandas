from __future__ import annotations

from itertools import product

import numpy as np
import pytest

from xarray import (
    DataArray,
    Dataset,
    MergeError,
    combine_by_coords,
    combine_nested,
    concat,
    merge,
)
from xarray.core import dtypes
from xarray.core.combine import (
    _check_shape_tile_ids,
    _combine_all_along_first_dim,
    _combine_nd,
    _infer_concat_order_from_coords,
    _infer_concat_order_from_positions,
    _new_tile_id,
)
from xarray.tests import assert_equal, assert_identical, requires_cftime
from xarray.tests.test_dataset import create_test_data


def assert_combined_tile_ids_equal(dict1, dict2):
    assert len(dict1) == len(dict2)
    for k, v in dict1.items():
        assert k in dict2.keys()
        assert_equal(dict1[k], dict2[k])


class TestTileIDsFromNestedList:
    def test_1d(self):
        ds = create_test_data
        input = [ds(0), ds(1)]

        expected = {(0,): ds(0), (1,): ds(1)}
        actual = _infer_concat_order_from_positions(input)
        assert_combined_tile_ids_equal(expected, actual)

    def test_2d(self):
        ds = create_test_data
        input = [[ds(0), ds(1)], [ds(2), ds(3)], [ds(4), ds(5)]]

        expected = {
            (0, 0): ds(0),
            (0, 1): ds(1),
            (1, 0): ds(2),
            (1, 1): ds(3),
            (2, 0): ds(4),
            (2, 1): ds(5),
        }
        actual = _infer_concat_order_from_positions(input)
        assert_combined_tile_ids_equal(expected, actual)

    def test_3d(self):
        ds = create_test_data
        input = [
            [[ds(0), ds(1)], [ds(2), ds(3)], [ds(4), ds(5)]],
            [[ds(6), ds(7)], [ds(8), ds(9)], [ds(10), ds(11)]],
        ]

        expected = {
            (0, 0, 0): ds(0),
            (0, 0, 1): ds(1),
            (0, 1, 0): ds(2),
            (0, 1, 1): ds(3),
            (0, 2, 0): ds(4),
            (0, 2, 1): ds(5),
            (1, 0, 0): ds(6),
            (1, 0, 1): ds(7),
            (1, 1, 0): ds(8),
            (1, 1, 1): ds(9),
            (1, 2, 0): ds(10),
            (1, 2, 1): ds(11),
        }
        actual = _infer_concat_order_from_positions(input)
        assert_combined_tile_ids_equal(expected, actual)

    def test_single_dataset(self):
        ds = create_test_data(0)
        input = [ds]

        expected = {(0,): ds}
        actual = _infer_concat_order_from_positions(input)
        assert_combined_tile_ids_equal(expected, actual)

    def test_redundant_nesting(self):
        ds = create_test_data
        input = [[ds(0)], [ds(1)]]

        expected = {(0, 0): ds(0), (1, 0): ds(1)}
        actual = _infer_concat_order_from_positions(input)
        assert_combined_tile_ids_equal(expected, actual)

    def test_ignore_empty_list(self):
        ds = create_test_data(0)
        input = [ds, []]
        expected = {(0,): ds}
        actual = _infer_concat_order_from_positions(input)
        assert_combined_tile_ids_equal(expected, actual)

    def test_uneven_depth_input(self):
        # Auto_combine won't work on ragged input
        # but this is just to increase test coverage
        ds = create_test_data
        input = [ds(0), [ds(1), ds(2)]]

        expected = {(0,): ds(0), (1, 0): ds(1), (1, 1): ds(2)}
        actual = _infer_concat_order_from_positions(input)
        assert_combined_tile_ids_equal(expected, actual)

    def test_uneven_length_input(self):
        # Auto_combine won't work on ragged input
        # but this is just to increase test coverage
        ds = create_test_data
        input = [[ds(0)], [ds(1), ds(2)]]

        expected = {(0, 0): ds(0), (1, 0): ds(1), (1, 1): ds(2)}
        actual = _infer_concat_order_from_positions(input)
        assert_combined_tile_ids_equal(expected, actual)

    def test_infer_from_datasets(self):
        ds = create_test_data
        input = [ds(0), ds(1)]

        expected = {(0,): ds(0), (1,): ds(1)}
        actual = _infer_concat_order_from_positions(input)
        assert_combined_tile_ids_equal(expected, actual)


class TestTileIDsFromCoords:
    def test_1d(self):
        ds0 = Dataset({"x": [0, 1]})
        ds1 = Dataset({"x": [2, 3]})

        expected = {(0,): ds0, (1,): ds1}
        actual, concat_dims = _infer_concat_order_from_coords([ds1, ds0])
        assert_combined_tile_ids_equal(expected, actual)
        assert concat_dims == ["x"]

    def test_2d(self):
        ds0 = Dataset({"x": [0, 1], "y": [10, 20, 30]})
        ds1 = Dataset({"x": [2, 3], "y": [10, 20, 30]})
        ds2 = Dataset({"x": [0, 1], "y": [40, 50, 60]})
        ds3 = Dataset({"x": [2, 3], "y": [40, 50, 60]})
        ds4 = Dataset({"x": [0, 1], "y": [70, 80, 90]})
        ds5 = Dataset({"x": [2, 3], "y": [70, 80, 90]})

        expected = {
            (0, 0): ds0,
            (1, 0): ds1,
            (0, 1): ds2,
            (1, 1): ds3,
            (0, 2): ds4,
            (1, 2): ds5,
        }
        actual, concat_dims = _infer_concat_order_from_coords(
            [ds1, ds0, ds3, ds5, ds2, ds4]
        )
        assert_combined_tile_ids_equal(expected, actual)
        assert concat_dims == ["x", "y"]

    def test_no_dimension_coords(self):
        ds0 = Dataset({"foo": ("x", [0, 1])})
        ds1 = Dataset({"foo": ("x", [2, 3])})
        with pytest.raises(ValueError, match=r"Could not find any dimension"):
            _infer_concat_order_from_coords([ds1, ds0])

    def test_coord_not_monotonic(self):
        ds0 = Dataset({"x": [0, 1]})
        ds1 = Dataset({"x": [3, 2]})
        with pytest.raises(
            ValueError,
            match=r"Coordinate variable x is neither monotonically increasing nor",
        ):
            _infer_concat_order_from_coords([ds1, ds0])

    def test_coord_monotonically_decreasing(self):
        ds0 = Dataset({"x": [3, 2]})
        ds1 = Dataset({"x": [1, 0]})

        expected = {(0,): ds0, (1,): ds1}
        actual, concat_dims = _infer_concat_order_from_coords([ds1, ds0])
        assert_combined_tile_ids_equal(expected, actual)
        assert concat_dims == ["x"]

    def test_no_concatenation_needed(self):
        ds = Dataset({"foo": ("x", [0, 1])})
        expected = {(): ds}
        actual, concat_dims = _infer_concat_order_from_coords([ds])
        assert_combined_tile_ids_equal(expected, actual)
        assert concat_dims == []

    def test_2d_plus_bystander_dim(self):
        ds0 = Dataset({"x": [0, 1], "y": [10, 20, 30], "t": [0.1, 0.2]})
        ds1 = Dataset({"x": [2, 3], "y": [10, 20, 30], "t": [0.1, 0.2]})
        ds2 = Dataset({"x": [0, 1], "y": [40, 50, 60], "t": [0.1, 0.2]})
        ds3 = Dataset({"x": [2, 3], "y": [40, 50, 60], "t": [0.1, 0.2]})

        expected = {(0, 0): ds0, (1, 0): ds1, (0, 1): ds2, (1, 1): ds3}
        actual, concat_dims = _infer_concat_order_from_coords([ds1, ds0, ds3, ds2])
        assert_combined_tile_ids_equal(expected, actual)
        assert concat_dims == ["x", "y"]

    def test_string_coords(self):
        ds0 = Dataset({"person": ["Alice", "Bob"]})
        ds1 = Dataset({"person": ["Caroline", "Daniel"]})

        expected = {(0,): ds0, (1,): ds1}
        actual, concat_dims = _infer_concat_order_from_coords([ds1, ds0])
        assert_combined_tile_ids_equal(expected, actual)
        assert concat_dims == ["person"]

    # Decided against natural sorting of string coords GH #2616
    def test_lexicographic_sort_string_coords(self):
        ds0 = Dataset({"simulation": ["run8", "run9"]})
        ds1 = Dataset({"simulation": ["run10", "run11"]})

        expected = {(0,): ds1, (1,): ds0}
        actual, concat_dims = _infer_concat_order_from_coords([ds1, ds0])
        assert_combined_tile_ids_equal(expected, actual)
        assert concat_dims == ["simulation"]

    def test_datetime_coords(self):
        ds0 = Dataset(
            {"time": np.array(["2000-03-06", "2000-03-07"], dtype="datetime64[ns]")}
        )
        ds1 = Dataset(
            {"time": np.array(["1999-01-01", "1999-02-04"], dtype="datetime64[ns]")}
        )

        expected = {(0,): ds1, (1,): ds0}
        actual, concat_dims = _infer_concat_order_from_coords([ds0, ds1])
        assert_combined_tile_ids_equal(expected, actual)
        assert concat_dims == ["time"]


@pytest.fixture(scope="module")
def create_combined_ids():
    return _create_combined_ids


def _create_combined_ids(shape):
    tile_ids = _create_tile_ids(shape)
    nums = range(len(tile_ids))
    return {tile_id: create_test_data(num) for tile_id, num in zip(tile_ids, nums)}


def _create_tile_ids(shape):
    tile_ids = product(*(range(i) for i in shape))
    return list(tile_ids)


class TestNewTileIDs:
    @pytest.mark.parametrize(
        "old_id, new_id",
        [((3, 0, 1), (0, 1)), ((0, 0), (0,)), ((1,), ()), ((0,), ()), ((1, 0), (0,))],
    )
    def test_new_tile_id(self, old_id, new_id):
        ds = create_test_data
        assert _new_tile_id((old_id, ds)) == new_id

    def test_get_new_tile_ids(self, create_combined_ids):
        shape = (1, 2, 3)
        combined_ids = create_combined_ids(shape)

        expected_tile_ids = sorted(combined_ids.keys())
        actual_tile_ids = _create_tile_ids(shape)
        assert expected_tile_ids == actual_tile_ids


class TestCombineND:
    @pytest.mark.parametrize("concat_dim", ["dim1", "new_dim"])
    def test_concat_once(self, create_combined_ids, concat_dim):
        shape = (2,)
        combined_ids = create_combined_ids(shape)
        ds = create_test_data
        result = _combine_all_along_first_dim(
            combined_ids,
            dim=concat_dim,
            data_vars="all",
            coords="different",
            compat="no_conflicts",
        )

        expected_ds = concat([ds(0), ds(1)], dim=concat_dim)
        assert_combined_tile_ids_equal(result, {(): expected_ds})

    def test_concat_only_first_dim(self, create_combined_ids):
        shape = (2, 3)
        combined_ids = create_combined_ids(shape)
        result = _combine_all_along_first_dim(
            combined_ids,
            dim="dim1",
            data_vars="all",
            coords="different",
            compat="no_conflicts",
        )

        ds = create_test_data
        partway1 = concat([ds(0), ds(3)], dim="dim1")
        partway2 = concat([ds(1), ds(4)], dim="dim1")
        partway3 = concat([ds(2), ds(5)], dim="dim1")
        expected_datasets = [partway1, partway2, partway3]
        expected = {(i,): ds for i, ds in enumerate(expected_datasets)}

        assert_combined_tile_ids_equal(result, expected)

    @pytest.mark.parametrize("concat_dim", ["dim1", "new_dim"])
    def test_concat_twice(self, create_combined_ids, concat_dim):
        shape = (2, 3)
        combined_ids = create_combined_ids(shape)
        result = _combine_nd(combined_ids, concat_dims=["dim1", concat_dim])

        ds = create_test_data
        partway1 = concat([ds(0), ds(3)], dim="dim1")
        partway2 = concat([ds(1), ds(4)], dim="dim1")
        partway3 = concat([ds(2), ds(5)], dim="dim1")
        expected = concat([partway1, partway2, partway3], dim=concat_dim)

        assert_equal(result, expected)


class TestCheckShapeTileIDs:
    def test_check_depths(self):
        ds = create_test_data(0)
        combined_tile_ids = {(0,): ds, (0, 1): ds}
        with pytest.raises(
            ValueError, match=r"sub-lists do not have consistent depths"
        ):
            _check_shape_tile_ids(combined_tile_ids)

    def test_check_lengths(self):
        ds = create_test_data(0)
        combined_tile_ids = {(0, 0): ds, (0, 1): ds, (0, 2): ds, (1, 0): ds, (1, 1): ds}
        with pytest.raises(
            ValueError, match=r"sub-lists do not have consistent lengths"
        ):
            _check_shape_tile_ids(combined_tile_ids)


class TestNestedCombine:
    def test_nested_concat(self):
        objs = [Dataset({"x": [0]}), Dataset({"x": [1]})]
        expected = Dataset({"x": [0, 1]})
        actual = combine_nested(objs, concat_dim="x")
        assert_identical(expected, actual)
        actual = combine_nested(objs, concat_dim=["x"])
        assert_identical(expected, actual)

        actual = combine_nested([actual], concat_dim=None)
        assert_identical(expected, actual)

        actual = combine_nested([actual], concat_dim="x")
        assert_identical(expected, actual)

        objs = [Dataset({"x": [0, 1]}), Dataset({"x": [2]})]
        actual = combine_nested(objs, concat_dim="x")
        expected = Dataset({"x": [0, 1, 2]})
        assert_identical(expected, actual)

        # ensure combine_nested handles non-sorted variables
        objs = [
            Dataset({"x": ("a", [0]), "y": ("a", [0])}),
            Dataset({"y": ("a", [1]), "x": ("a", [1])}),
        ]
        actual = combine_nested(objs, concat_dim="a")
        expected = Dataset({"x": ("a", [0, 1]), "y": ("a", [0, 1])})
        assert_identical(expected, actual)

        objs = [Dataset({"x": [0], "y": [0]}), Dataset({"x": [1]})]
        actual = combine_nested(objs, concat_dim="x")
        expected = Dataset({"x": [0, 1], "y": [0]})
        assert_identical(expected, actual)

    @pytest.mark.parametrize(
        "join, expected",
        [
            ("outer", Dataset({"x": [0, 1], "y": [0, 1]})),
            ("inner", Dataset({"x": [0, 1], "y": []})),
            ("left", Dataset({"x": [0, 1], "y": [0]})),
            ("right", Dataset({"x": [0, 1], "y": [1]})),
        ],
    )
    def test_combine_nested_join(self, join, expected):
        objs = [Dataset({"x": [0], "y": [0]}), Dataset({"x": [1], "y": [1]})]
        actual = combine_nested(objs, concat_dim="x", join=join)
        assert_identical(expected, actual)

    def test_combine_nested_join_exact(self):
        objs = [Dataset({"x": [0], "y": [0]}), Dataset({"x": [1], "y": [1]})]
        with pytest.raises(ValueError, match=r"cannot align.*join.*exact"):
            combine_nested(objs, concat_dim="x", join="exact")

    def test_empty_input(self):
        assert_identical(Dataset(), combine_nested([], concat_dim="x"))

    # Fails because of concat's weird treatment of dimension coords, see #2975
    @pytest.mark.xfail
    def test_nested_concat_too_many_dims_at_once(self):
        objs = [Dataset({"x": [0], "y": [1]}), Dataset({"y": [0], "x": [1]})]
        with pytest.raises(ValueError, match="not equal across datasets"):
            combine_nested(objs, concat_dim="x", coords="minimal")

    def test_nested_concat_along_new_dim(self):
        objs = [
            Dataset({"a": ("x", [10]), "x": [0]}),
            Dataset({"a": ("x", [20]), "x": [0]}),
        ]
        expected = Dataset({"a": (("t", "x"), [[10], [20]]), "x": [0]})
        actual = combine_nested(objs, concat_dim="t")
        assert_identical(expected, actual)

        # Same but with a DataArray as new dim, see GH #1988 and #2647
        dim = DataArray([100, 150], name="baz", dims="baz")
        expected = Dataset(
            {"a": (("baz", "x"), [[10], [20]]), "x": [0], "baz": [100, 150]}
        )
        actual = combine_nested(objs, concat_dim=dim)
        assert_identical(expected, actual)

    def test_nested_merge(self):
        data = Dataset({"x": 0})
        actual = combine_nested([data, data, data], concat_dim=None)
        assert_identical(data, actual)

        ds1 = Dataset({"a": ("x", [1, 2]), "x": [0, 1]})
        ds2 = Dataset({"a": ("x", [2, 3]), "x": [1, 2]})
        expected = Dataset({"a": ("x", [1, 2, 3]), "x": [0, 1, 2]})
        actual = combine_nested([ds1, ds2], concat_dim=None)
        assert_identical(expected, actual)
        actual = combine_nested([ds1, ds2], concat_dim=[None])
        assert_identical(expected, actual)

        tmp1 = Dataset({"x": 0})
        tmp2 = Dataset({"x": np.nan})
        actual = combine_nested([tmp1, tmp2], concat_dim=None)
        assert_identical(tmp1, actual)
        actual = combine_nested([tmp1, tmp2], concat_dim=[None])
        assert_identical(tmp1, actual)

        # Single object, with a concat_dim explicitly provided
        # Test the issue reported in GH #1988
        objs = [Dataset({"x": 0, "y": 1})]
        dim = DataArray([100], name="baz", dims="baz")
        actual = combine_nested(objs, concat_dim=[dim])
        expected = Dataset({"x": ("baz", [0]), "y": ("baz", [1])}, {"baz": [100]})
        assert_identical(expected, actual)

        # Just making sure that auto_combine is doing what is
        # expected for non-scalar values, too.
        objs = [Dataset({"x": ("z", [0, 1]), "y": ("z", [1, 2])})]
        dim = DataArray([100], name="baz", dims="baz")
        actual = combine_nested(objs, concat_dim=[dim])
        expected = Dataset(
            {"x": (("baz", "z"), [[0, 1]]), "y": (("baz", "z"), [[1, 2]])},
            {"baz": [100]},
        )
        assert_identical(expected, actual)

    def test_concat_multiple_dims(self):
        objs = [
            [Dataset({"a": (("x", "y"), [[0]])}), Dataset({"a": (("x", "y"), [[1]])})],
            [Dataset({"a": (("x", "y"), [[2]])}), Dataset({"a": (("x", "y"), [[3]])})],
        ]
        actual = combine_nested(objs, concat_dim=["x", "y"])
        expected = Dataset({"a": (("x", "y"), [[0, 1], [2, 3]])})
        assert_identical(expected, actual)

    def test_concat_name_symmetry(self):
        """Inspired by the discussion on GH issue #2777"""

        da1 = DataArray(name="a", data=[[0]], dims=["x", "y"])
        da2 = DataArray(name="b", data=[[1]], dims=["x", "y"])
        da3 = DataArray(name="a", data=[[2]], dims=["x", "y"])
        da4 = DataArray(name="b", data=[[3]], dims=["x", "y"])

        x_first = combine_nested([[da1, da2], [da3, da4]], concat_dim=["x", "y"])
        y_first = combine_nested([[da1, da3], [da2, da4]], concat_dim=["y", "x"])

        assert_identical(x_first, y_first)

    def test_concat_one_dim_merge_another(self):
        data = create_test_data(add_attrs=False)

        data1 = data.copy(deep=True)
        data2 = data.copy(deep=True)

        objs = [
            [data1.var1.isel(dim2=slice(4)), data2.var1.isel(dim2=slice(4, 9))],
            [data1.var2.isel(dim2=slice(4)), data2.var2.isel(dim2=slice(4, 9))],
        ]

        expected = data[["var1", "var2"]]
        actual = combine_nested(objs, concat_dim=[None, "dim2"])
        assert_identical(expected, actual)

    def test_auto_combine_2d(self):
        ds = create_test_data

        partway1 = concat([ds(0), ds(3)], dim="dim1")
        partway2 = concat([ds(1), ds(4)], dim="dim1")
        partway3 = concat([ds(2), ds(5)], dim="dim1")
        expected = concat([partway1, partway2, partway3], dim="dim2")

        datasets = [[ds(0), ds(1), ds(2)], [ds(3), ds(4), ds(5)]]
        result = combine_nested(datasets, concat_dim=["dim1", "dim2"])
        assert_equal(result, expected)

    def test_auto_combine_2d_combine_attrs_kwarg(self):
        ds = lambda x: create_test_data(x, add_attrs=False)

        partway1 = concat([ds(0), ds(3)], dim="dim1")
        partway2 = concat([ds(1), ds(4)], dim="dim1")
        partway3 = concat([ds(2), ds(5)], dim="dim1")
        expected = concat([partway1, partway2, partway3], dim="dim2")

        expected_dict = {}
        expected_dict["drop"] = expected.copy(deep=True)
        expected_dict["drop"].attrs = {}
        expected_dict["no_conflicts"] = expected.copy(deep=True)
        expected_dict["no_conflicts"].attrs = {
            "a": 1,
            "b": 2,
            "c": 3,
            "d": 4,
            "e": 5,
            "f": 6,
        }
        expected_dict["override"] = expected.copy(deep=True)
        expected_dict["override"].attrs = {"a": 1}
        f = lambda attrs, context: attrs[0]
        expected_dict[f] = expected.copy(deep=True)
        expected_dict[f].attrs = f([{"a": 1}], None)

        datasets = [[ds(0), ds(1), ds(2)], [ds(3), ds(4), ds(5)]]

        datasets[0][0].attrs = {"a": 1}
        datasets[0][1].attrs = {"a": 1, "b": 2}
        datasets[0][2].attrs = {"a": 1, "c": 3}
        datasets[1][0].attrs = {"a": 1, "d": 4}
        datasets[1][1].attrs = {"a": 1, "e": 5}
        datasets[1][2].attrs = {"a": 1, "f": 6}

        with pytest.raises(ValueError, match=r"combine_attrs='identical'"):
            result = combine_nested(
                datasets, concat_dim=["dim1", "dim2"], combine_attrs="identical"
            )

        for combine_attrs in expected_dict:
            result = combine_nested(
                datasets, concat_dim=["dim1", "dim2"], combine_attrs=combine_attrs
            )
            assert_identical(result, expected_dict[combine_attrs])

    def test_combine_nested_missing_data_new_dim(self):
        # Your data includes "time" and "station" dimensions, and each year's
        # data has a different set of stations.
        datasets = [
            Dataset({"a": ("x", [2, 3]), "x": [1, 2]}),
            Dataset({"a": ("x", [1, 2]), "x": [0, 1]}),
        ]
        expected = Dataset(
            {"a": (("t", "x"), [[np.nan, 2, 3], [1, 2, np.nan]])}, {"x": [0, 1, 2]}
        )
        actual = combine_nested(datasets, concat_dim="t")
        assert_identical(expected, actual)

    def test_invalid_hypercube_input(self):
        ds = create_test_data

        datasets = [[ds(0), ds(1), ds(2)], [ds(3), ds(4)]]
        with pytest.raises(
            ValueError, match=r"sub-lists do not have consistent lengths"
        ):
            combine_nested(datasets, concat_dim=["dim1", "dim2"])

        datasets = [[ds(0), ds(1)], [[ds(3), ds(4)]]]
        with pytest.raises(
            ValueError, match=r"sub-lists do not have consistent depths"
        ):
            combine_nested(datasets, concat_dim=["dim1", "dim2"])

        datasets = [[ds(0), ds(1)], [ds(3), ds(4)]]
        with pytest.raises(ValueError, match=r"concat_dims has length"):
            combine_nested(datasets, concat_dim=["dim1"])

    def test_merge_one_dim_concat_another(self):
        objs = [
            [Dataset({"foo": ("x", [0, 1])}), Dataset({"bar": ("x", [10, 20])})],
            [Dataset({"foo": ("x", [2, 3])}), Dataset({"bar": ("x", [30, 40])})],
        ]
        expected = Dataset({"foo": ("x", [0, 1, 2, 3]), "bar": ("x", [10, 20, 30, 40])})

        actual = combine_nested(objs, concat_dim=["x", None], compat="equals")
        assert_identical(expected, actual)

        # Proving it works symmetrically
        objs = [
            [Dataset({"foo": ("x", [0, 1])}), Dataset({"foo": ("x", [2, 3])})],
            [Dataset({"bar": ("x", [10, 20])}), Dataset({"bar": ("x", [30, 40])})],
        ]
        actual = combine_nested(objs, concat_dim=[None, "x"], compat="equals")
        assert_identical(expected, actual)

    def test_combine_concat_over_redundant_nesting(self):
        objs = [[Dataset({"x": [0]}), Dataset({"x": [1]})]]
        actual = combine_nested(objs, concat_dim=[None, "x"])
        expected = Dataset({"x": [0, 1]})
        assert_identical(expected, actual)

        objs = [[Dataset({"x": [0]})], [Dataset({"x": [1]})]]
        actual = combine_nested(objs, concat_dim=["x", None])
        expected = Dataset({"x": [0, 1]})
        assert_identical(expected, actual)

        objs = [[Dataset({"x": [0]})]]
        actual = combine_nested(objs, concat_dim=[None, None])
        expected = Dataset({"x": [0]})
        assert_identical(expected, actual)

    @pytest.mark.parametrize("fill_value", [dtypes.NA, 2, 2.0, {"a": 2, "b": 1}])
    def test_combine_nested_fill_value(self, fill_value):
        datasets = [
            Dataset({"a": ("x", [2, 3]), "b": ("x", [-2, 1]), "x": [1, 2]}),
            Dataset({"a": ("x", [1, 2]), "b": ("x", [3, -1]), "x": [0, 1]}),
        ]
        if fill_value == dtypes.NA:
            # if we supply the default, we expect the missing value for a
            # float array
            fill_value_a = fill_value_b = np.nan
        elif isinstance(fill_value, dict):
            fill_value_a = fill_value["a"]
            fill_value_b = fill_value["b"]
        else:
            fill_value_a = fill_value_b = fill_value
        expected = Dataset(
            {
                "a": (("t", "x"), [[fill_value_a, 2, 3], [1, 2, fill_value_a]]),
                "b": (("t", "x"), [[fill_value_b, -2, 1], [3, -1, fill_value_b]]),
            },
            {"x": [0, 1, 2]},
        )
        actual = combine_nested(datasets, concat_dim="t", fill_value=fill_value)
        assert_identical(expected, actual)

    def test_combine_nested_unnamed_data_arrays(self):
        unnamed_array = DataArray(data=[1.0, 2.0], coords={"x": [0, 1]}, dims="x")

        actual = combine_nested([unnamed_array], concat_dim="x")
        expected = unnamed_array
        assert_identical(expected, actual)

        unnamed_array1 = DataArray(data=[1.0, 2.0], coords={"x": [0, 1]}, dims="x")
        unnamed_array2 = DataArray(data=[3.0, 4.0], coords={"x": [2, 3]}, dims="x")

        actual = combine_nested([unnamed_array1, unnamed_array2], concat_dim="x")
        expected = DataArray(
            data=[1.0, 2.0, 3.0, 4.0], coords={"x": [0, 1, 2, 3]}, dims="x"
        )
        assert_identical(expected, actual)

        da1 = DataArray(data=[[0.0]], coords={"x": [0], "y": [0]}, dims=["x", "y"])
        da2 = DataArray(data=[[1.0]], coords={"x": [0], "y": [1]}, dims=["x", "y"])
        da3 = DataArray(data=[[2.0]], coords={"x": [1], "y": [0]}, dims=["x", "y"])
        da4 = DataArray(data=[[3.0]], coords={"x": [1], "y": [1]}, dims=["x", "y"])
        objs = [[da1, da2], [da3, da4]]

        expected = DataArray(
            data=[[0.0, 1.0], [2.0, 3.0]],
            coords={"x": [0, 1], "y": [0, 1]},
            dims=["x", "y"],
        )
        actual = combine_nested(objs, concat_dim=["x", "y"])
        assert_identical(expected, actual)

    # TODO aijams - Determine if this test is appropriate.
    def test_nested_combine_mixed_datasets_arrays(self):
        objs = [
            DataArray([0, 1], dims=("x"), coords=({"x": [0, 1]})),
            Dataset({"x": [2, 3]}),
        ]
        with pytest.raises(
            ValueError, match=r"Can't combine datasets with unnamed arrays."
        ):
            combine_nested(objs, "x")


class TestCombineDatasetsbyCoords:
    def test_combine_by_coords(self):
        objs = [Dataset({"x": [0]}), Dataset({"x": [1]})]
        actual = combine_by_coords(objs)
        expected = Dataset({"x": [0, 1]})
        assert_identical(expected, actual)

        actual = combine_by_coords([actual])
        assert_identical(expected, actual)

        objs = [Dataset({"x": [0, 1]}), Dataset({"x": [2]})]
        actual = combine_by_coords(objs)
        expected = Dataset({"x": [0, 1, 2]})
        assert_identical(expected, actual)

        # ensure auto_combine handles non-sorted variables
        objs = [
            Dataset({"x": ("a", [0]), "y": ("a", [0]), "a": [0]}),
            Dataset({"x": ("a", [1]), "y": ("a", [1]), "a": [1]}),
        ]
        actual = combine_by_coords(objs)
        expected = Dataset({"x": ("a", [0, 1]), "y": ("a", [0, 1]), "a": [0, 1]})
        assert_identical(expected, actual)

        objs = [Dataset({"x": [0], "y": [0]}), Dataset({"y": [1], "x": [1]})]
        actual = combine_by_coords(objs)
        expected = Dataset({"x": [0, 1], "y": [0, 1]})
        assert_equal(actual, expected)

        objs = [Dataset({"x": 0}), Dataset({"x": 1})]
        with pytest.raises(
            ValueError, match=r"Could not find any dimension coordinates"
        ):
            combine_by_coords(objs)

        objs = [Dataset({"x": [0], "y": [0]}), Dataset({"x": [0]})]
        with pytest.raises(ValueError, match=r"Every dimension needs a coordinate"):
            combine_by_coords(objs)

    def test_empty_input(self):
        assert_identical(Dataset(), combine_by_coords([]))

    @pytest.mark.parametrize(
        "join, expected",
        [
            ("outer", Dataset({"x": [0, 1], "y": [0, 1]})),
            ("inner", Dataset({"x": [0, 1], "y": []})),
            ("left", Dataset({"x": [0, 1], "y": [0]})),
            ("right", Dataset({"x": [0, 1], "y": [1]})),
        ],
    )
    def test_combine_coords_join(self, join, expected):
        objs = [Dataset({"x": [0], "y": [0]}), Dataset({"x": [1], "y": [1]})]
        actual = combine_nested(objs, concat_dim="x", join=join)
        assert_identical(expected, actual)

    def test_combine_coords_join_exact(self):
        objs = [Dataset({"x": [0], "y": [0]}), Dataset({"x": [1], "y": [1]})]
        with pytest.raises(ValueError, match=r"cannot align.*join.*exact.*"):
            combine_nested(objs, concat_dim="x", join="exact")

    @pytest.mark.parametrize(
        "combine_attrs, expected",
        [
            ("drop", Dataset({"x": [0, 1], "y": [0, 1]}, attrs={})),
            (
                "no_conflicts",
                Dataset({"x": [0, 1], "y": [0, 1]}, attrs={"a": 1, "b": 2}),
            ),
            ("override", Dataset({"x": [0, 1], "y": [0, 1]}, attrs={"a": 1})),
            (
                lambda attrs, context: attrs[1],
                Dataset({"x": [0, 1], "y": [0, 1]}, attrs={"a": 1, "b": 2}),
            ),
        ],
    )
    def test_combine_coords_combine_attrs(self, combine_attrs, expected):
        objs = [
            Dataset({"x": [0], "y": [0]}, attrs={"a": 1}),
            Dataset({"x": [1], "y": [1]}, attrs={"a": 1, "b": 2}),
        ]
        actual = combine_nested(
            objs, concat_dim="x", join="outer", combine_attrs=combine_attrs
        )
        assert_identical(expected, actual)

        if combine_attrs == "no_conflicts":
            objs[1].attrs["a"] = 2
            with pytest.raises(ValueError, match=r"combine_attrs='no_conflicts'"):
                actual = combine_nested(
                    objs, concat_dim="x", join="outer", combine_attrs=combine_attrs
                )

    def test_combine_coords_combine_attrs_identical(self):
        objs = [
            Dataset({"x": [0], "y": [0]}, attrs={"a": 1}),
            Dataset({"x": [1], "y": [1]}, attrs={"a": 1}),
        ]
        expected = Dataset({"x": [0, 1], "y": [0, 1]}, attrs={"a": 1})
        actual = combine_nested(
            objs, concat_dim="x", join="outer", combine_attrs="identical"
        )
        assert_identical(expected, actual)

        objs[1].attrs["b"] = 2

        with pytest.raises(ValueError, match=r"combine_attrs='identical'"):
            actual = combine_nested(
                objs, concat_dim="x", join="outer", combine_attrs="identical"
            )

    def test_combine_nested_combine_attrs_drop_conflicts(self):
        objs = [
            Dataset({"x": [0], "y": [0]}, attrs={"a": 1, "b": 2, "c": 3}),
            Dataset({"x": [1], "y": [1]}, attrs={"a": 1, "b": 0, "d": 3}),
        ]
        expected = Dataset({"x": [0, 1], "y": [0, 1]}, attrs={"a": 1, "c": 3, "d": 3})
        actual = combine_nested(
            objs, concat_dim="x", join="outer", combine_attrs="drop_conflicts"
        )
        assert_identical(expected, actual)

    @pytest.mark.parametrize(
        "combine_attrs, attrs1, attrs2, expected_attrs, expect_exception",
        [
            (
                "no_conflicts",
                {"a": 1, "b": 2},
                {"a": 1, "c": 3},
                {"a": 1, "b": 2, "c": 3},
                False,
            ),
            ("no_conflicts", {"a": 1, "b": 2}, {}, {"a": 1, "b": 2}, False),
            ("no_conflicts", {}, {"a": 1, "c": 3}, {"a": 1, "c": 3}, False),
            (
                "no_conflicts",
                {"a": 1, "b": 2},
                {"a": 4, "c": 3},
                {"a": 1, "b": 2, "c": 3},
                True,
            ),
            ("drop", {"a": 1, "b": 2}, {"a": 1, "c": 3}, {}, False),
            ("identical", {"a": 1, "b": 2}, {"a": 1, "b": 2}, {"a": 1, "b": 2}, False),
            ("identical", {"a": 1, "b": 2}, {"a": 1, "c": 3}, {"a": 1, "b": 2}, True),
            (
                "override",
                {"a": 1, "b": 2},
                {"a": 4, "b": 5, "c": 3},
                {"a": 1, "b": 2},
                False,
            ),
            (
                "drop_conflicts",
                {"a": 1, "b": 2, "c": 3},
                {"b": 1, "c": 3, "d": 4},
                {"a": 1, "c": 3, "d": 4},
                False,
            ),
        ],
    )
    def test_combine_nested_combine_attrs_variables(
        self, combine_attrs, attrs1, attrs2, expected_attrs, expect_exception
    ):
        """check that combine_attrs is used on data variables and coords"""
        data1 = Dataset(
            {
                "a": ("x", [1, 2], attrs1),
                "b": ("x", [3, -1], attrs1),
                "x": ("x", [0, 1], attrs1),
            }
        )
        data2 = Dataset(
            {
                "a": ("x", [2, 3], attrs2),
                "b": ("x", [-2, 1], attrs2),
                "x": ("x", [2, 3], attrs2),
            }
        )

        if expect_exception:
            with pytest.raises(MergeError, match="combine_attrs"):
                combine_by_coords([data1, data2], combine_attrs=combine_attrs)
        else:
            actual = combine_by_coords([data1, data2], combine_attrs=combine_attrs)
            expected = Dataset(
                {
                    "a": ("x", [1, 2, 2, 3], expected_attrs),
                    "b": ("x", [3, -1, -2, 1], expected_attrs),
                },
                {"x": ("x", [0, 1, 2, 3], expected_attrs)},
            )

            assert_identical(actual, expected)

    @pytest.mark.parametrize(
        "combine_attrs, attrs1, attrs2, expected_attrs, expect_exception",
        [
            (
                "no_conflicts",
                {"a": 1, "b": 2},
                {"a": 1, "c": 3},
                {"a": 1, "b": 2, "c": 3},
                False,
            ),
            ("no_conflicts", {"a": 1, "b": 2}, {}, {"a": 1, "b": 2}, False),
            ("no_conflicts", {}, {"a": 1, "c": 3}, {"a": 1, "c": 3}, False),
            (
                "no_conflicts",
                {"a": 1, "b": 2},
                {"a": 4, "c": 3},
                {"a": 1, "b": 2, "c": 3},
                True,
            ),
            ("drop", {"a": 1, "b": 2}, {"a": 1, "c": 3}, {}, False),
            ("identical", {"a": 1, "b": 2}, {"a": 1, "b": 2}, {"a": 1, "b": 2}, False),
            ("identical", {"a": 1, "b": 2}, {"a": 1, "c": 3}, {"a": 1, "b": 2}, True),
            (
                "override",
                {"a": 1, "b": 2},
                {"a": 4, "b": 5, "c": 3},
                {"a": 1, "b": 2},
                False,
            ),
            (
                "drop_conflicts",
                {"a": 1, "b": 2, "c": 3},
                {"b": 1, "c": 3, "d": 4},
                {"a": 1, "c": 3, "d": 4},
                False,
            ),
        ],
    )
    def test_combine_by_coords_combine_attrs_variables(
        self, combine_attrs, attrs1, attrs2, expected_attrs, expect_exception
    ):
        """check that combine_attrs is used on data variables and coords"""
        data1 = Dataset(
            {"x": ("a", [0], attrs1), "y": ("a", [0], attrs1), "a": ("a", [0], attrs1)}
        )
        data2 = Dataset(
            {"x": ("a", [1], attrs2), "y": ("a", [1], attrs2), "a": ("a", [1], attrs2)}
        )

        if expect_exception:
            with pytest.raises(MergeError, match="combine_attrs"):
                combine_by_coords([data1, data2], combine_attrs=combine_attrs)
        else:
            actual = combine_by_coords([data1, data2], combine_attrs=combine_attrs)
            expected = Dataset(
                {
                    "x": ("a", [0, 1], expected_attrs),
                    "y": ("a", [0, 1], expected_attrs),
                    "a": ("a", [0, 1], expected_attrs),
                }
            )

            assert_identical(actual, expected)

    def test_infer_order_from_coords(self):
        data = create_test_data()
        objs = [data.isel(dim2=slice(4, 9)), data.isel(dim2=slice(4))]
        actual = combine_by_coords(objs)
        expected = data
        assert expected.broadcast_equals(actual)

    def test_combine_leaving_bystander_dimensions(self):
        # Check non-monotonic bystander dimension coord doesn't raise
        # ValueError on combine (https://github.com/pydata/xarray/issues/3150)
        ycoord = ["a", "c", "b"]

        data = np.random.rand(7, 3)

        ds1 = Dataset(
            data_vars=dict(data=(["x", "y"], data[:3, :])),
            coords=dict(x=[1, 2, 3], y=ycoord),
        )

        ds2 = Dataset(
            data_vars=dict(data=(["x", "y"], data[3:, :])),
            coords=dict(x=[4, 5, 6, 7], y=ycoord),
        )

        expected = Dataset(
            data_vars=dict(data=(["x", "y"], data)),
            coords=dict(x=[1, 2, 3, 4, 5, 6, 7], y=ycoord),
        )

        actual = combine_by_coords((ds1, ds2))
        assert_identical(expected, actual)

    def test_combine_by_coords_previously_failed(self):
        # In the above scenario, one file is missing, containing the data for
        # one year's data for one variable.
        datasets = [
            Dataset({"a": ("x", [0]), "x": [0]}),
            Dataset({"b": ("x", [0]), "x": [0]}),
            Dataset({"a": ("x", [1]), "x": [1]}),
        ]
        expected = Dataset({"a": ("x", [0, 1]), "b": ("x", [0, np.nan])}, {"x": [0, 1]})
        actual = combine_by_coords(datasets)
        assert_identical(expected, actual)

    def test_combine_by_coords_still_fails(self):
        # concat can't handle new variables (yet):
        # https://github.com/pydata/xarray/issues/508
        datasets = [Dataset({"x": 0}, {"y": 0}), Dataset({"x": 1}, {"y": 1, "z": 1})]
        with pytest.raises(ValueError):
            combine_by_coords(datasets, "y")

    def test_combine_by_coords_no_concat(self):
        objs = [Dataset({"x": 0}), Dataset({"y": 1})]
        actual = combine_by_coords(objs)
        expected = Dataset({"x": 0, "y": 1})
        assert_identical(expected, actual)

        objs = [Dataset({"x": 0, "y": 1}), Dataset({"y": np.nan, "z": 2})]
        actual = combine_by_coords(objs)
        expected = Dataset({"x": 0, "y": 1, "z": 2})
        assert_identical(expected, actual)

    def test_check_for_impossible_ordering(self):
        ds0 = Dataset({"x": [0, 1, 5]})
        ds1 = Dataset({"x": [2, 3]})
        with pytest.raises(
            ValueError,
            match=r"does not have monotonic global indexes along dimension x",
        ):
            combine_by_coords([ds1, ds0])

    def test_combine_by_coords_incomplete_hypercube(self):
        # test that this succeeds with default fill_value
        x1 = Dataset({"a": (("y", "x"), [[1]])}, coords={"y": [0], "x": [0]})
        x2 = Dataset({"a": (("y", "x"), [[1]])}, coords={"y": [1], "x": [0]})
        x3 = Dataset({"a": (("y", "x"), [[1]])}, coords={"y": [0], "x": [1]})
        actual = combine_by_coords([x1, x2, x3])
        expected = Dataset(
            {"a": (("y", "x"), [[1, 1], [1, np.nan]])},
            coords={"y": [0, 1], "x": [0, 1]},
        )
        assert_identical(expected, actual)

        # test that this fails if fill_value is None
        with pytest.raises(ValueError):
            combine_by_coords([x1, x2, x3], fill_value=None)


class TestCombineMixedObjectsbyCoords:
    def test_combine_by_coords_mixed_unnamed_dataarrays(self):
        named_da = DataArray(name="a", data=[1.0, 2.0], coords={"x": [0, 1]}, dims="x")
        unnamed_da = DataArray(data=[3.0, 4.0], coords={"x": [2, 3]}, dims="x")

        with pytest.raises(
            ValueError, match="Can't automatically combine unnamed DataArrays with"
        ):
            combine_by_coords([named_da, unnamed_da])

        da = DataArray([0, 1], dims="x", coords=({"x": [0, 1]}))
        ds = Dataset({"x": [2, 3]})
        with pytest.raises(
            ValueError,
            match="Can't automatically combine unnamed DataArrays with",
        ):
            combine_by_coords([da, ds])

    def test_combine_coords_mixed_datasets_named_dataarrays(self):
        da = DataArray(name="a", data=[4, 5], dims="x", coords=({"x": [0, 1]}))
        ds = Dataset({"b": ("x", [2, 3])})
        actual = combine_by_coords([da, ds])
        expected = Dataset(
            {"a": ("x", [4, 5]), "b": ("x", [2, 3])}, coords={"x": ("x", [0, 1])}
        )
        assert_identical(expected, actual)

    def test_combine_by_coords_all_unnamed_dataarrays(self):
        unnamed_array = DataArray(data=[1.0, 2.0], coords={"x": [0, 1]}, dims="x")

        actual = combine_by_coords([unnamed_array])
        expected = unnamed_array
        assert_identical(expected, actual)

        unnamed_array1 = DataArray(data=[1.0, 2.0], coords={"x": [0, 1]}, dims="x")
        unnamed_array2 = DataArray(data=[3.0, 4.0], coords={"x": [2, 3]}, dims="x")

        actual = combine_by_coords([unnamed_array1, unnamed_array2])
        expected = DataArray(
            data=[1.0, 2.0, 3.0, 4.0], coords={"x": [0, 1, 2, 3]}, dims="x"
        )
        assert_identical(expected, actual)

    def test_combine_by_coords_all_named_dataarrays(self):
        named_da = DataArray(name="a", data=[1.0, 2.0], coords={"x": [0, 1]}, dims="x")

        actual = combine_by_coords([named_da])
        expected = named_da.to_dataset()
        assert_identical(expected, actual)

        named_da1 = DataArray(name="a", data=[1.0, 2.0], coords={"x": [0, 1]}, dims="x")
        named_da2 = DataArray(name="b", data=[3.0, 4.0], coords={"x": [2, 3]}, dims="x")

        actual = combine_by_coords([named_da1, named_da2])
        expected = Dataset(
            {
                "a": DataArray(data=[1.0, 2.0], coords={"x": [0, 1]}, dims="x"),
                "b": DataArray(data=[3.0, 4.0], coords={"x": [2, 3]}, dims="x"),
            }
        )
        assert_identical(expected, actual)

    def test_combine_by_coords_all_dataarrays_with_the_same_name(self):
        named_da1 = DataArray(name="a", data=[1.0, 2.0], coords={"x": [0, 1]}, dims="x")
        named_da2 = DataArray(name="a", data=[3.0, 4.0], coords={"x": [2, 3]}, dims="x")

        actual = combine_by_coords([named_da1, named_da2])
        expected = merge([named_da1, named_da2])
        assert_identical(expected, actual)


@requires_cftime
def test_combine_by_coords_distant_cftime_dates():
    # Regression test for https://github.com/pydata/xarray/issues/3535
    import cftime

    time_1 = [cftime.DatetimeGregorian(4500, 12, 31)]
    time_2 = [cftime.DatetimeGregorian(4600, 12, 31)]
    time_3 = [cftime.DatetimeGregorian(5100, 12, 31)]

    da_1 = DataArray([0], dims=["time"], coords=[time_1], name="a").to_dataset()
    da_2 = DataArray([1], dims=["time"], coords=[time_2], name="a").to_dataset()
    da_3 = DataArray([2], dims=["time"], coords=[time_3], name="a").to_dataset()

    result = combine_by_coords([da_1, da_2, da_3])

    expected_time = np.concatenate([time_1, time_2, time_3])
    expected = DataArray(
        [0, 1, 2], dims=["time"], coords=[expected_time], name="a"
    ).to_dataset()
    assert_identical(result, expected)


@requires_cftime
def test_combine_by_coords_raises_for_differing_calendars():
    # previously failed with uninformative StopIteration instead of TypeError
    # https://github.com/pydata/xarray/issues/4495

    import cftime

    time_1 = [cftime.DatetimeGregorian(2000, 1, 1)]
    time_2 = [cftime.DatetimeProlepticGregorian(2001, 1, 1)]

    da_1 = DataArray([0], dims=["time"], coords=[time_1], name="a").to_dataset()
    da_2 = DataArray([1], dims=["time"], coords=[time_2], name="a").to_dataset()

    error_msg = (
        "Cannot combine along dimension 'time' with mixed types."
        " Found:.*"
        " If importing data directly from a file then setting"
        " `use_cftime=True` may fix this issue."
    )
    with pytest.raises(TypeError, match=error_msg):
        combine_by_coords([da_1, da_2])


def test_combine_by_coords_raises_for_differing_types():
    # str and byte cannot be compared
    da_1 = DataArray([0], dims=["time"], coords=[["a"]], name="a").to_dataset()
    da_2 = DataArray([1], dims=["time"], coords=[[b"b"]], name="a").to_dataset()

    with pytest.raises(
        TypeError, match=r"Cannot combine along dimension 'time' with mixed types."
    ):
        combine_by_coords([da_1, da_2])
