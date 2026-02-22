import re

import numpy as np
import pytest

import xarray as xr
from xarray.core.datatree_mapping import map_over_datasets
from xarray.core.treenode import TreeIsomorphismError
from xarray.testing import assert_equal, assert_identical

empty = xr.Dataset()


class TestMapOverSubTree:
    def test_no_trees_passed(self):
        with pytest.raises(TypeError, match="must pass at least one tree object"):
            map_over_datasets(lambda x: x, "dt")

    def test_not_isomorphic(self, create_test_datatree):
        dt1 = create_test_datatree()
        dt2 = create_test_datatree()
        dt2["set1/set2/extra"] = xr.DataTree(name="extra")

        with pytest.raises(
            TreeIsomorphismError,
            match=re.escape(
                r"children at node 'set1/set2' do not match: [] vs ['extra']"
            ),
        ):
            map_over_datasets(lambda x, y: None, dt1, dt2)

    def test_no_trees_returned(self, create_test_datatree):
        dt1 = create_test_datatree()
        dt2 = create_test_datatree()
        expected = xr.DataTree.from_dict(dict.fromkeys(dt1.to_dict()))
        actual = map_over_datasets(lambda x, y: None, dt1, dt2)
        assert_equal(expected, actual)

    def test_single_tree_arg(self, create_test_datatree):
        dt = create_test_datatree()
        expected = create_test_datatree(modify=lambda x: 10.0 * x)
        result_tree = map_over_datasets(lambda x: 10 * x, dt)
        assert_equal(result_tree, expected)

    def test_single_tree_arg_plus_arg(self, create_test_datatree):
        dt = create_test_datatree()
        expected = create_test_datatree(modify=lambda ds: (10.0 * ds))
        result_tree = map_over_datasets(lambda x, y: x * y, dt, 10.0)
        assert_equal(result_tree, expected)

        result_tree = map_over_datasets(lambda x, y: x * y, 10.0, dt)
        assert_equal(result_tree, expected)

    def test_single_tree_arg_plus_kwarg(self, create_test_datatree):
        dt = create_test_datatree()
        expected = create_test_datatree(modify=lambda ds: (10.0 * ds))

        def multiply_by_kwarg(ds, **kwargs):
            ds = ds * kwargs.pop("multiplier")
            return ds

        result_tree = map_over_datasets(
            multiply_by_kwarg, dt, kwargs=dict(multiplier=10.0)
        )
        assert_equal(result_tree, expected)

    def test_multiple_tree_args(self, create_test_datatree):
        dt1 = create_test_datatree()
        dt2 = create_test_datatree()
        expected = create_test_datatree(modify=lambda ds: 2.0 * ds)
        result = map_over_datasets(lambda x, y: x + y, dt1, dt2)
        assert_equal(result, expected)

    def test_return_multiple_trees(self, create_test_datatree):
        dt = create_test_datatree()
        dt_min, dt_max = map_over_datasets(lambda x: (x.min(), x.max()), dt)
        expected_min = create_test_datatree(modify=lambda ds: ds.min())
        assert_equal(dt_min, expected_min)
        expected_max = create_test_datatree(modify=lambda ds: ds.max())
        assert_equal(dt_max, expected_max)

    def test_return_wrong_type(self, simple_datatree):
        dt1 = simple_datatree

        with pytest.raises(
            TypeError,
            match=re.escape(
                "the result of calling func on the node at position '.' is not a "
                "Dataset or None or a tuple of such types"
            ),
        ):
            map_over_datasets(lambda x: "string", dt1)  # type: ignore[arg-type,return-value]

    def test_return_tuple_of_wrong_types(self, simple_datatree):
        dt1 = simple_datatree

        with pytest.raises(
            TypeError,
            match=re.escape(
                "the result of calling func on the node at position '.' is not a "
                "Dataset or None or a tuple of such types"
            ),
        ):
            map_over_datasets(lambda x: (x, "string"), dt1)  # type: ignore[arg-type,return-value]

    def test_return_inconsistent_number_of_results(self, simple_datatree):
        dt1 = simple_datatree

        with pytest.raises(
            TypeError,
            match=re.escape(
                r"Calling func on the nodes at position set1 returns a tuple "
                "of 0 datasets, whereas calling func on the nodes at position "
                ". instead returns a tuple of 2 datasets."
            ),
        ):
            # Datasets in simple_datatree have different numbers of dims
            map_over_datasets(lambda ds: tuple((None,) * len(ds.dims)), dt1)

    def test_wrong_number_of_arguments_for_func(self, simple_datatree):
        dt = simple_datatree

        with pytest.raises(
            TypeError, match="takes 1 positional argument but 2 were given"
        ):
            map_over_datasets(lambda x: 10 * x, dt, dt)

    def test_map_single_dataset_against_whole_tree(self, create_test_datatree):
        dt = create_test_datatree()

        def nodewise_merge(node_ds, fixed_ds):
            return xr.merge([node_ds, fixed_ds])

        other_ds = xr.Dataset({"z": ("z", [0])})
        expected = create_test_datatree(modify=lambda ds: xr.merge([ds, other_ds]))
        result_tree = map_over_datasets(nodewise_merge, dt, other_ds)
        assert_equal(result_tree, expected)

    @pytest.mark.xfail
    def test_trees_with_different_node_names(self):
        # TODO test this after I've got good tests for renaming nodes
        raise NotImplementedError

    def test_tree_method(self, create_test_datatree):
        dt = create_test_datatree()

        def multiply(ds, times):
            return times * ds

        expected = create_test_datatree(modify=lambda ds: 10.0 * ds)
        result_tree = dt.map_over_datasets(multiply, 10.0)
        assert_equal(result_tree, expected)

    def test_tree_method_with_kwarg(self, create_test_datatree):
        dt = create_test_datatree()

        def multiply(ds, **kwargs):
            return kwargs.pop("times") * ds

        expected = create_test_datatree(modify=lambda ds: 10.0 * ds)
        result_tree = dt.map_over_datasets(multiply, kwargs=dict(times=10.0))
        assert_equal(result_tree, expected)

    def test_discard_ancestry(self, create_test_datatree):
        # Check for datatree GH issue https://github.com/xarray-contrib/datatree/issues/48
        dt = create_test_datatree()
        subtree = dt["set1"]
        expected = create_test_datatree(modify=lambda ds: 10.0 * ds)["set1"]
        result_tree = map_over_datasets(lambda x: 10.0 * x, subtree)
        assert_equal(result_tree, expected)

    def test_keep_attrs_on_empty_nodes(self, create_test_datatree):
        # GH278
        dt = create_test_datatree()
        dt["set1/set2"].attrs["foo"] = "bar"

        def empty_func(ds):
            return ds

        result = dt.map_over_datasets(empty_func)
        assert result["set1/set2"].attrs == dt["set1/set2"].attrs

    def test_error_contains_path_of_offending_node(self, create_test_datatree):
        dt = create_test_datatree()
        dt["set1"]["bad_var"] = 0
        print(dt)

        def fail_on_specific_node(ds):
            if "bad_var" in ds:
                raise ValueError("Failed because 'bar_var' present in dataset")

        with pytest.raises(
            ValueError,
            match=re.escape(
                r"Raised whilst mapping function over node(s) with path 'set1'"
            ),
        ):
            dt.map_over_datasets(fail_on_specific_node)

    def test_inherited_coordinates_with_index(self):
        root = xr.Dataset(coords={"x": [1, 2]})
        child = xr.Dataset({"foo": ("x", [0, 1])})  # no coordinates
        tree = xr.DataTree.from_dict({"/": root, "/child": child})
        actual = tree.map_over_datasets(lambda ds: ds)  # identity
        assert isinstance(actual, xr.DataTree)
        assert_identical(tree, actual)

        actual_child = actual.children["child"].to_dataset(inherit=False)
        assert_identical(actual_child, child)


class TestMutableOperations:
    def test_construct_using_type(self):
        # from datatree GH issue https://github.com/xarray-contrib/datatree/issues/188
        # xarray's .weighted is unusual because it uses type() to create a Dataset/DataArray

        a = xr.DataArray(
            np.random.rand(3, 4, 10),
            dims=["x", "y", "time"],
            coords={"area": (["x", "y"], np.random.rand(3, 4))},
        ).to_dataset(name="data")
        b = xr.DataArray(
            np.random.rand(2, 6, 14),
            dims=["x", "y", "time"],
            coords={"area": (["x", "y"], np.random.rand(2, 6))},
        ).to_dataset(name="data")
        dt = xr.DataTree.from_dict({"a": a, "b": b})

        def weighted_mean(ds):
            if "area" not in ds.coords:
                return None
            return ds.weighted(ds.area).mean(["x", "y"])

        dt.map_over_datasets(weighted_mean)

    def test_alter_inplace_forbidden(self):
        simpsons = xr.DataTree.from_dict(
            {
                "/": xr.Dataset({"age": 83}),
                "/Herbert": xr.Dataset({"age": 40}),
                "/Homer": xr.Dataset({"age": 39}),
                "/Homer/Bart": xr.Dataset({"age": 10}),
                "/Homer/Lisa": xr.Dataset({"age": 8}),
                "/Homer/Maggie": xr.Dataset({"age": 1}),
            },
            name="Abe",
        )

        def fast_forward(ds: xr.Dataset, years: float) -> xr.Dataset:
            """Add some years to the age, but by altering the given dataset"""
            ds["age"] = ds["age"] + years
            return ds

        with pytest.raises(AttributeError):
            simpsons.map_over_datasets(fast_forward, 10)
