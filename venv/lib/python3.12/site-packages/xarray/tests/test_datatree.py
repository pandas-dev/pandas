import re
import sys
import typing
from collections.abc import Callable, Mapping
from copy import copy, deepcopy
from textwrap import dedent

import numpy as np
import pytest

import xarray as xr
from xarray import DataArray, Dataset
from xarray.core.coordinates import DataTreeCoordinates
from xarray.core.datatree import DataTree
from xarray.core.treenode import NotFoundInTreeError
from xarray.testing import assert_equal, assert_identical
from xarray.tests import (
    assert_array_equal,
    create_test_data,
    requires_dask,
    source_ndarray,
)

ON_WINDOWS = sys.platform == "win32"


class TestTreeCreation:
    def test_empty(self) -> None:
        dt = DataTree(name="root")
        assert dt.name == "root"
        assert dt.parent is None
        assert dt.children == {}
        assert_identical(dt.to_dataset(), xr.Dataset())

    def test_name(self) -> None:
        dt = DataTree()
        assert dt.name is None

        dt = DataTree(name="foo")
        assert dt.name == "foo"

        dt.name = "bar"
        assert dt.name == "bar"

        dt = DataTree(children={"foo": DataTree()})
        assert dt["/foo"].name == "foo"
        with pytest.raises(
            ValueError, match="cannot set the name of a node which already has a parent"
        ):
            dt["/foo"].name = "bar"

        detached = dt["/foo"].copy()
        assert detached.name == "foo"
        detached.name = "bar"
        assert detached.name == "bar"

    def test_bad_names(self) -> None:
        with pytest.raises(TypeError):
            DataTree(name=5)  # type: ignore[arg-type]

        with pytest.raises(ValueError):
            DataTree(name="folder/data")

    def test_data_arg(self) -> None:
        ds = xr.Dataset({"foo": 42})
        tree: DataTree = DataTree(dataset=ds)
        assert_identical(tree.to_dataset(), ds)

        with pytest.raises(TypeError):
            DataTree(dataset=xr.DataArray(42, name="foo"))  # type: ignore[arg-type]

    def test_child_data_not_copied(self) -> None:
        # regression test for https://github.com/pydata/xarray/issues/9683
        class NoDeepCopy:
            def __deepcopy__(self, memo):
                raise TypeError("class can't be deepcopied")

        da = xr.DataArray(NoDeepCopy())
        ds = xr.Dataset({"var": da})
        dt1 = xr.DataTree(ds)
        dt2 = xr.DataTree(ds, children={"child": dt1})
        dt3 = xr.DataTree.from_dict({"/": ds, "child": ds})
        assert_identical(dt2, dt3)


class TestFamilyTree:
    def test_dont_modify_children_inplace(self) -> None:
        # GH issue 9196
        child = DataTree()
        DataTree(children={"child": child})
        assert child.parent is None

    def test_create_two_children(self) -> None:
        root_data = xr.Dataset({"a": ("y", [6, 7, 8]), "set0": ("x", [9, 10])})
        set1_data = xr.Dataset({"a": 0, "b": 1})
        root = DataTree.from_dict(
            {"/": root_data, "/set1": set1_data, "/set1/set2": None}
        )
        assert root["/set1"].name == "set1"
        assert root["/set1/set2"].name == "set2"

    def test_create_full_tree(self, simple_datatree) -> None:
        d = simple_datatree.to_dict()
        d_keys = list(d.keys())

        expected_keys = [
            "/",
            "/set1",
            "/set2",
            "/set3",
            "/set1/set1",
            "/set1/set2",
            "/set2/set1",
        ]

        assert d_keys == expected_keys


class TestNames:
    def test_child_gets_named_on_attach(self) -> None:
        sue = DataTree()
        mary = DataTree(children={"Sue": sue})
        assert mary.children["Sue"].name == "Sue"

    def test_dataset_containing_slashes(self) -> None:
        xda: xr.DataArray = xr.DataArray(
            [[1, 2]],
            coords={"label": ["a"], "R30m/y": [30, 60]},
        )
        xds: xr.Dataset = xr.Dataset({"group/subgroup/my_variable": xda})
        with pytest.raises(
            ValueError,
            match=re.escape(
                "Given variables have names containing the '/' character: "
                "['R30m/y', 'group/subgroup/my_variable']. "
                "Variables stored in DataTree objects cannot have names containing '/' characters, "
                "as this would make path-like access to variables ambiguous."
            ),
        ):
            DataTree(xds)


class TestPaths:
    def test_path_property(self) -> None:
        john = DataTree.from_dict(
            {
                "/Mary/Sue": DataTree(),
            }
        )
        assert john["/Mary/Sue"].path == "/Mary/Sue"
        assert john.path == "/"

    def test_path_roundtrip(self) -> None:
        john = DataTree.from_dict(
            {
                "/Mary/Sue": DataTree(),
            }
        )
        assert john["/Mary/Sue"].name == "Sue"

    def test_same_tree(self) -> None:
        john = DataTree.from_dict(
            {
                "/Mary": DataTree(),
                "/Kate": DataTree(),
            }
        )
        mary = john.children["Mary"]
        kate = john.children["Kate"]
        assert mary.same_tree(kate)

    def test_relative_paths(self) -> None:
        john = DataTree.from_dict(
            {
                "/Mary/Sue": DataTree(),
                "/Annie": DataTree(),
            }
        )
        sue = john.children["Mary"].children["Sue"]
        annie = john.children["Annie"]
        assert sue.relative_to(john) == "Mary/Sue"
        assert john.relative_to(sue) == "../.."
        assert annie.relative_to(sue) == "../../Annie"
        assert sue.relative_to(annie) == "../Mary/Sue"
        assert sue.relative_to(sue) == "."

        evil_kate = DataTree()
        with pytest.raises(
            NotFoundInTreeError, match="nodes do not lie within the same tree"
        ):
            sue.relative_to(evil_kate)


class TestStoreDatasets:
    def test_create_with_data(self) -> None:
        dat = xr.Dataset({"a": 0})
        john = DataTree(name="john", dataset=dat)

        assert_identical(john.to_dataset(), dat)

        with pytest.raises(TypeError):
            DataTree(name="mary", dataset="junk")  # type: ignore[arg-type]

    def test_set_data(self) -> None:
        john = DataTree(name="john")
        dat = xr.Dataset({"a": 0})
        john.dataset = dat  # type: ignore[assignment,unused-ignore]

        assert_identical(john.to_dataset(), dat)

        with pytest.raises(TypeError):
            john.dataset = "junk"  # type: ignore[assignment]

    def test_has_data(self) -> None:
        john = DataTree(name="john", dataset=xr.Dataset({"a": 0}))
        assert john.has_data

        john_no_data = DataTree(name="john", dataset=None)
        assert not john_no_data.has_data

    def test_is_hollow(self) -> None:
        john = DataTree(dataset=xr.Dataset({"a": 0}))
        assert john.is_hollow

        eve = DataTree(children={"john": john})
        assert eve.is_hollow

        eve.dataset = xr.Dataset({"a": 1})  # type: ignore[assignment,unused-ignore]
        assert not eve.is_hollow


class TestToDataset:
    def test_to_dataset_inherited(self) -> None:
        base = xr.Dataset(coords={"a": [1], "b": 2})
        sub = xr.Dataset(coords={"c": [3]})
        tree = DataTree.from_dict({"/": base, "/sub": sub})
        subtree = typing.cast(DataTree, tree["sub"])

        assert_identical(tree.to_dataset(inherit=False), base)
        assert_identical(subtree.to_dataset(inherit=False), sub)

        sub_and_base = xr.Dataset(coords={"a": [1], "c": [3]})  # no "b"
        assert_identical(tree.to_dataset(inherit=True), base)
        assert_identical(subtree.to_dataset(inherit=True), sub_and_base)


class TestVariablesChildrenNameCollisions:
    def test_parent_already_has_variable_with_childs_name(self) -> None:
        with pytest.raises(KeyError, match="already contains a variable named a"):
            DataTree.from_dict({"/": xr.Dataset({"a": [0], "b": 1}), "/a": None})

    def test_parent_already_has_variable_with_childs_name_update(self) -> None:
        dt = DataTree(dataset=xr.Dataset({"a": [0], "b": 1}))
        with pytest.raises(ValueError, match="already contains a variable named a"):
            dt.update({"a": DataTree()})

    def test_assign_when_already_child_with_variables_name(self) -> None:
        dt = DataTree.from_dict(
            {
                "/a": DataTree(),
            }
        )

        with pytest.raises(ValueError, match="node already contains a variable"):
            dt.dataset = xr.Dataset({"a": 0})  # type: ignore[assignment,unused-ignore]

        dt.dataset = xr.Dataset()  # type: ignore[assignment,unused-ignore]

        new_ds = dt.to_dataset().assign(a=xr.DataArray(0))
        with pytest.raises(ValueError, match="node already contains a variable"):
            dt.dataset = new_ds  # type: ignore[assignment,unused-ignore]


class TestGet: ...


class TestGetItem:
    def test_getitem_node(self) -> None:
        folder1 = DataTree.from_dict(
            {
                "/results/highres": DataTree(),
            }
        )

        assert folder1["results"].name == "results"
        assert folder1["results/highres"].name == "highres"

    def test_getitem_self(self) -> None:
        dt = DataTree()
        assert dt["."] is dt

    def test_getitem_single_data_variable(self) -> None:
        data = xr.Dataset({"temp": [0, 50]})
        results = DataTree(name="results", dataset=data)
        assert_identical(results["temp"], data["temp"])

    def test_getitem_single_data_variable_from_node(self) -> None:
        data = xr.Dataset({"temp": [0, 50]})
        folder1 = DataTree.from_dict(
            {
                "/results/highres": data,
            }
        )
        assert_identical(folder1["results/highres/temp"], data["temp"])

    def test_getitem_nonexistent_node(self) -> None:
        folder1 = DataTree.from_dict({"/results": DataTree()}, name="folder1")
        with pytest.raises(KeyError):
            folder1["results/highres"]

    def test_getitem_nonexistent_variable(self) -> None:
        data = xr.Dataset({"temp": [0, 50]})
        results = DataTree(name="results", dataset=data)
        with pytest.raises(KeyError):
            results["pressure"]

    @pytest.mark.xfail(reason="Should be deprecated in favour of .subset")
    def test_getitem_multiple_data_variables(self) -> None:
        data = xr.Dataset({"temp": [0, 50], "p": [5, 8, 7]})
        results = DataTree(name="results", dataset=data)
        assert_identical(results[["temp", "p"]], data[["temp", "p"]])  # type: ignore[index]

    @pytest.mark.xfail(
        reason="Indexing needs to return whole tree (GH https://github.com/xarray-contrib/datatree/issues/77)"
    )
    def test_getitem_dict_like_selection_access_to_dataset(self) -> None:
        data = xr.Dataset({"temp": [0, 50]})
        results = DataTree(name="results", dataset=data)
        assert_identical(results[{"temp": 1}], data[{"temp": 1}])  # type: ignore[index]


class TestUpdate:
    def test_update(self) -> None:
        dt = DataTree()
        dt.update({"foo": xr.DataArray(0), "a": DataTree()})
        expected = DataTree.from_dict({"/": xr.Dataset({"foo": 0}), "a": None})
        assert_equal(dt, expected)
        assert dt.groups == ("/", "/a")

    def test_update_new_named_dataarray(self) -> None:
        da = xr.DataArray(name="temp", data=[0, 50])
        folder1 = DataTree(name="folder1")
        folder1.update({"results": da})
        expected = da.rename("results")
        assert_equal(folder1["results"], expected)

    def test_update_doesnt_alter_child_name(self) -> None:
        dt = DataTree()
        dt.update({"foo": xr.DataArray(0), "a": DataTree(name="b")})
        assert "a" in dt.children
        child = dt["a"]
        assert child.name == "a"

    def test_update_overwrite(self) -> None:
        actual = DataTree.from_dict({"a": DataTree(xr.Dataset({"x": 1}))})
        actual.update({"a": DataTree(xr.Dataset({"x": 2}))})
        expected = DataTree.from_dict({"a": DataTree(xr.Dataset({"x": 2}))})
        assert_equal(actual, expected)

    def test_update_coordinates(self) -> None:
        expected = DataTree.from_dict({"/": xr.Dataset(coords={"a": 1})})
        actual = DataTree.from_dict({"/": xr.Dataset()})
        actual.update(xr.Dataset(coords={"a": 1}))
        assert_equal(actual, expected)

    def test_update_inherited_coords(self) -> None:
        expected = DataTree.from_dict(
            {
                "/": xr.Dataset(coords={"a": 1}),
                "/b": xr.Dataset(coords={"c": 1}),
            }
        )
        actual = DataTree.from_dict(
            {
                "/": xr.Dataset(coords={"a": 1}),
                "/b": xr.Dataset(),
            }
        )
        actual["/b"].update(xr.Dataset(coords={"c": 1}))
        assert_identical(actual, expected)

        # DataTree.identical() currently does not require that non-inherited
        # coordinates are defined identically, so we need to check this
        # explicitly
        actual_node = actual.children["b"].to_dataset(inherit=False)
        expected_node = expected.children["b"].to_dataset(inherit=False)
        assert_identical(actual_node, expected_node)


class TestCopy:
    def test_copy(self, create_test_datatree) -> None:
        dt = create_test_datatree()

        for node in dt.root.subtree:
            node.attrs["Test"] = [1, 2, 3]

        for copied in [dt.copy(deep=False), copy(dt)]:
            assert_identical(dt, copied)

            for node, copied_node in zip(
                dt.root.subtree, copied.root.subtree, strict=True
            ):
                assert node.encoding == copied_node.encoding
                # Note: IndexVariable objects with string dtype are always
                # copied because of xarray.core.util.safe_cast_to_index.
                # Limiting the test to data variables.
                for k in node.data_vars:
                    v0 = node.variables[k]
                    v1 = copied_node.variables[k]
                    assert source_ndarray(v0.data) is source_ndarray(v1.data)
                copied_node["foo"] = xr.DataArray(data=np.arange(5), dims="z")
                assert "foo" not in node

                copied_node.attrs["foo"] = "bar"
                assert "foo" not in node.attrs
                assert node.attrs["Test"] is copied_node.attrs["Test"]

    def test_copy_subtree(self) -> None:
        dt = DataTree.from_dict({"/level1/level2/level3": xr.Dataset()})

        actual = dt["/level1/level2"].copy()
        expected = DataTree.from_dict({"/level3": xr.Dataset()}, name="level2")

        assert_identical(actual, expected)

    def test_copy_coord_inheritance(self) -> None:
        tree = DataTree.from_dict(
            {"/": xr.Dataset(coords={"x": [0, 1]}), "/c": DataTree()}
        )
        actual = tree.copy()
        node_ds = actual.children["c"].to_dataset(inherit=False)
        assert_identical(node_ds, xr.Dataset())

        actual = tree.children["c"].copy()
        expected = DataTree(Dataset(coords={"x": [0, 1]}), name="c")
        assert_identical(expected, actual)

        actual = tree.children["c"].copy(inherit=False)
        expected = DataTree(name="c")
        assert_identical(expected, actual)

    def test_deepcopy(self, create_test_datatree) -> None:
        dt = create_test_datatree()

        for node in dt.root.subtree:
            node.attrs["Test"] = [1, 2, 3]

        for copied in [dt.copy(deep=True), deepcopy(dt)]:
            assert_identical(dt, copied)

            for node, copied_node in zip(
                dt.root.subtree, copied.root.subtree, strict=True
            ):
                assert node.encoding == copied_node.encoding
                # Note: IndexVariable objects with string dtype are always
                # copied because of xarray.core.util.safe_cast_to_index.
                # Limiting the test to data variables.
                for k in node.data_vars:
                    v0 = node.variables[k]
                    v1 = copied_node.variables[k]
                    assert source_ndarray(v0.data) is not source_ndarray(v1.data)
                copied_node["foo"] = xr.DataArray(data=np.arange(5), dims="z")
                assert "foo" not in node

                copied_node.attrs["foo"] = "bar"
                assert "foo" not in node.attrs
                assert node.attrs["Test"] is not copied_node.attrs["Test"]

    @pytest.mark.xfail(reason="data argument not yet implemented")
    def test_copy_with_data(self, create_test_datatree) -> None:
        orig = create_test_datatree()
        # TODO use .data_vars once that property is available
        data_vars = {
            k: v for k, v in orig.variables.items() if k not in orig._coord_names
        }
        new_data = {k: np.random.randn(*v.shape) for k, v in data_vars.items()}
        actual = orig.copy(data=new_data)

        expected = orig.copy()
        for k, v in new_data.items():
            expected[k].data = v
        assert_identical(expected, actual)

        # TODO test parents and children?


class TestSetItem:
    def test_setitem_new_child_node(self) -> None:
        john = DataTree(name="john")
        mary = DataTree(name="mary")
        john["mary"] = mary

        grafted_mary = john["mary"]
        assert grafted_mary.parent is john
        assert grafted_mary.name == "mary"

    def test_setitem_unnamed_child_node_becomes_named(self) -> None:
        john2 = DataTree(name="john2")
        john2["sonny"] = DataTree()
        assert john2["sonny"].name == "sonny"

    def test_setitem_new_grandchild_node(self) -> None:
        john = DataTree.from_dict({"/Mary/Rose": DataTree()})
        new_rose = DataTree(dataset=xr.Dataset({"x": 0}))
        john["Mary/Rose"] = new_rose

        grafted_rose = john["Mary/Rose"]
        assert grafted_rose.parent is john["/Mary"]
        assert grafted_rose.name == "Rose"

    def test_grafted_subtree_retains_name(self) -> None:
        subtree = DataTree(name="original_subtree_name")
        root = DataTree(name="root")
        root["new_subtree_name"] = subtree
        assert subtree.name == "original_subtree_name"

    def test_setitem_new_empty_node(self) -> None:
        john = DataTree(name="john")
        john["mary"] = DataTree()
        mary = john["mary"]
        assert isinstance(mary, DataTree)
        assert_identical(mary.to_dataset(), xr.Dataset())

    def test_setitem_overwrite_data_in_node_with_none(self) -> None:
        john = DataTree.from_dict({"/mary": xr.Dataset()}, name="john")

        john["mary"] = DataTree()
        assert_identical(john["mary"].to_dataset(), xr.Dataset())

        john.dataset = xr.Dataset()  # type: ignore[assignment,unused-ignore]
        with pytest.raises(ValueError, match="has no name"):
            john["."] = DataTree()

    @pytest.mark.xfail(reason="assigning Datasets doesn't yet create new nodes")
    def test_setitem_dataset_on_this_node(self) -> None:
        data = xr.Dataset({"temp": [0, 50]})
        results = DataTree(name="results")
        results["."] = data
        assert_identical(results.to_dataset(), data)

    def test_setitem_dataset_as_new_node(self) -> None:
        data = xr.Dataset({"temp": [0, 50]})
        folder1 = DataTree(name="folder1")
        folder1["results"] = data
        assert_identical(folder1["results"].to_dataset(), data)

    def test_setitem_dataset_as_new_node_requiring_intermediate_nodes(self) -> None:
        data = xr.Dataset({"temp": [0, 50]})
        folder1 = DataTree(name="folder1")
        folder1["results/highres"] = data
        assert_identical(folder1["results/highres"].to_dataset(), data)

    def test_setitem_named_dataarray(self) -> None:
        da = xr.DataArray(name="temp", data=[0, 50])
        folder1 = DataTree(name="folder1")
        folder1["results"] = da
        expected = da.rename("results")
        assert_equal(folder1["results"], expected)

    def test_setitem_unnamed_dataarray(self) -> None:
        data = xr.DataArray([0, 50])
        folder1 = DataTree(name="folder1")
        folder1["results"] = data
        assert_equal(folder1["results"], data)

    def test_setitem_variable(self) -> None:
        var = xr.Variable(data=[0, 50], dims="x")
        folder1 = DataTree(name="folder1")
        folder1["results"] = var
        assert_equal(folder1["results"], xr.DataArray(var))

    def test_setitem_coerce_to_dataarray(self) -> None:
        folder1 = DataTree(name="folder1")
        folder1["results"] = 0
        assert_equal(folder1["results"], xr.DataArray(0))

    def test_setitem_add_new_variable_to_empty_node(self) -> None:
        results = DataTree(name="results")
        results["pressure"] = xr.DataArray(data=[2, 3])
        assert "pressure" in results.dataset
        results["temp"] = xr.Variable(data=[10, 11], dims=["x"])
        assert "temp" in results.dataset

        # What if there is a path to traverse first?
        results_with_path = DataTree(name="results")
        results_with_path["highres/pressure"] = xr.DataArray(data=[2, 3])
        assert "pressure" in results_with_path["highres"].dataset
        results_with_path["highres/temp"] = xr.Variable(data=[10, 11], dims=["x"])
        assert "temp" in results_with_path["highres"].dataset

    def test_setitem_dataarray_replace_existing_node(self) -> None:
        t = xr.Dataset({"temp": [0, 50]})
        results = DataTree(name="results", dataset=t)
        p = xr.DataArray(data=[2, 3])
        results["pressure"] = p
        expected = t.assign(pressure=p)
        assert_identical(results.to_dataset(), expected)


class TestCoords:
    def test_properties(self) -> None:
        # use int64 for repr consistency on windows
        ds = Dataset(
            data_vars={
                "foo": (["x", "y"], np.random.randn(2, 3)),
            },
            coords={
                "x": ("x", np.array([-1, -2], "int64")),
                "y": ("y", np.array([0, 1, 2], "int64")),
                "a": ("x", np.array([4, 5], "int64")),
                "b": np.int64(-10),
            },
        )
        dt = DataTree(dataset=ds)
        dt["child"] = DataTree()

        coords = dt.coords
        assert isinstance(coords, DataTreeCoordinates)

        # len
        assert len(coords) == 4

        # iter
        assert list(coords) == ["x", "y", "a", "b"]

        assert_identical(coords["x"].variable, dt["x"].variable)
        assert_identical(coords["y"].variable, dt["y"].variable)

        assert "x" in coords
        assert "a" in coords
        assert 0 not in coords
        assert "foo" not in coords
        assert "child" not in coords

        with pytest.raises(KeyError):
            coords["foo"]

        # TODO this currently raises a ValueError instead of a KeyError
        # with pytest.raises(KeyError):
        #     coords[0]

        # repr
        expected = dedent(
            """\
        Coordinates:
          * x        (x) int64 16B -1 -2
          * y        (y) int64 24B 0 1 2
            a        (x) int64 16B 4 5
            b        int64 8B -10"""
        )
        actual = repr(coords)
        assert expected == actual

        # dims
        assert coords.sizes == {"x": 2, "y": 3}

        # dtypes
        assert coords.dtypes == {
            "x": np.dtype("int64"),
            "y": np.dtype("int64"),
            "a": np.dtype("int64"),
            "b": np.dtype("int64"),
        }

    def test_modify(self) -> None:
        ds = Dataset(
            data_vars={
                "foo": (["x", "y"], np.random.randn(2, 3)),
            },
            coords={
                "x": ("x", np.array([-1, -2], "int64")),
                "y": ("y", np.array([0, 1, 2], "int64")),
                "a": ("x", np.array([4, 5], "int64")),
                "b": np.int64(-10),
            },
        )
        dt = DataTree(dataset=ds)
        dt["child"] = DataTree()

        actual = dt.copy(deep=True)
        actual.coords["x"] = ("x", ["a", "b"])
        assert_array_equal(actual["x"], ["a", "b"])

        actual = dt.copy(deep=True)
        actual.coords["z"] = ("z", ["a", "b"])
        assert_array_equal(actual["z"], ["a", "b"])

        actual = dt.copy(deep=True)
        with pytest.raises(ValueError, match=r"conflicting dimension sizes"):
            actual.coords["x"] = ("x", [-1])
        assert_identical(actual, dt)  # should not be modified

        # TODO: re-enable after implementing reset_coords()
        # actual = dt.copy()
        # del actual.coords["b"]
        # expected = dt.reset_coords("b", drop=True)
        # assert_identical(expected, actual)

        with pytest.raises(KeyError):
            del dt.coords["not_found"]

        with pytest.raises(KeyError):
            del dt.coords["foo"]

        # TODO: re-enable after implementing assign_coords()
        # actual = dt.copy(deep=True)
        # actual.coords.update({"c": 11})
        # expected = dt.assign_coords({"c": 11})
        # assert_identical(expected, actual)

        # # regression test for GH3746
        # del actual.coords["x"]
        # assert "x" not in actual.xindexes

        # test that constructors can also handle the `DataTreeCoordinates` object
        ds2 = Dataset(coords=dt.coords)
        assert_identical(ds2.coords, dt.coords)
        da = DataArray(coords=dt.coords)
        assert_identical(da.coords, dt.coords)

        # DataTree constructor doesn't accept coords= but should still be able to handle DatasetCoordinates
        dt2 = DataTree(dataset=dt.coords)
        assert_identical(dt2.coords, dt.coords)

    def test_inherited(self) -> None:
        ds = Dataset(
            data_vars={
                "foo": (["x", "y"], np.random.randn(2, 3)),
            },
            coords={
                "x": ("x", np.array([-1, -2], "int64")),
                "y": ("y", np.array([0, 1, 2], "int64")),
                "a": ("x", np.array([4, 5], "int64")),
                "b": np.int64(-10),
            },
        )
        dt = DataTree(dataset=ds)
        dt["child"] = DataTree()
        child = dt["child"]

        assert set(dt.coords) == {"x", "y", "a", "b"}
        assert set(child.coords) == {"x", "y"}

        actual = child.copy(deep=True)
        actual.coords["x"] = ("x", ["a", "b"])
        assert_array_equal(actual["x"], ["a", "b"])

        actual = child.copy(deep=True)
        actual.coords.update({"c": 11})
        expected = child.copy(deep=True)
        expected.coords["c"] = 11
        # check we have only altered the child node
        assert_identical(expected.root, actual.root)

        with pytest.raises(KeyError):
            # cannot delete inherited coordinate from child node
            del child["x"]

        # TODO requires a fix for #9472
        # actual = child.copy(deep=True)
        # actual.coords.update({"c": 11})
        # expected = child.assign_coords({"c": 11})
        # assert_identical(expected, actual)


def test_delitem() -> None:
    ds = Dataset({"a": 0}, coords={"x": ("x", [1, 2]), "z": "a"})
    dt = DataTree(ds, children={"c": DataTree()})

    with pytest.raises(KeyError):
        del dt["foo"]

    # test delete children
    del dt["c"]
    assert dt.children == {}
    assert set(dt.variables) == {"x", "z", "a"}
    with pytest.raises(KeyError):
        del dt["c"]

    # test delete variables
    del dt["a"]
    assert set(dt.coords) == {"x", "z"}
    with pytest.raises(KeyError):
        del dt["a"]

    # test delete coordinates
    del dt["z"]
    assert set(dt.coords) == {"x"}
    with pytest.raises(KeyError):
        del dt["z"]

    # test delete indexed coordinates
    del dt["x"]
    assert dt.variables == {}
    assert dt.coords == {}
    assert dt.indexes == {}
    with pytest.raises(KeyError):
        del dt["x"]


class TestTreeFromDict:
    def test_data_in_root(self) -> None:
        dat = xr.Dataset()
        dt = DataTree.from_dict({"/": dat})
        assert dt.name is None
        assert dt.parent is None
        assert dt.children == {}
        assert_identical(dt.to_dataset(), dat)

    def test_one_layer(self) -> None:
        dat1, dat2 = xr.Dataset({"a": 1}), xr.Dataset({"b": 2})
        dt = DataTree.from_dict({"run1": dat1, "run2": dat2})
        assert_identical(dt.to_dataset(), xr.Dataset())
        assert dt.name is None
        assert_identical(dt["run1"].to_dataset(), dat1)
        assert dt["run1"].children == {}
        assert_identical(dt["run2"].to_dataset(), dat2)
        assert dt["run2"].children == {}

    def test_two_layers(self) -> None:
        dat1, dat2 = xr.Dataset({"a": 1}), xr.Dataset({"a": [1, 2]})
        dt = DataTree.from_dict({"highres/run": dat1, "lowres/run": dat2})
        assert "highres" in dt.children
        assert "lowres" in dt.children
        highres_run = dt["highres/run"]
        assert_identical(highres_run.to_dataset(), dat1)

    def test_nones(self) -> None:
        dt = DataTree.from_dict({"d": None, "d/e": None})
        assert [node.name for node in dt.subtree] == [None, "d", "e"]
        assert [node.path for node in dt.subtree] == ["/", "/d", "/d/e"]
        assert_identical(dt["d/e"].to_dataset(), xr.Dataset())

    def test_full(self, simple_datatree) -> None:
        dt = simple_datatree
        paths = [node.path for node in dt.subtree]
        assert paths == [
            "/",
            "/set1",
            "/set2",
            "/set3",
            "/set1/set1",
            "/set1/set2",
            "/set2/set1",
        ]

    def test_datatree_values(self) -> None:
        dat1 = DataTree(dataset=xr.Dataset({"a": 1}))
        expected = DataTree()
        expected["a"] = dat1

        actual = DataTree.from_dict({"a": dat1})

        assert_identical(actual, expected)

    def test_roundtrip_to_dict(self, simple_datatree) -> None:
        tree = simple_datatree
        roundtrip = DataTree.from_dict(tree.to_dict())
        assert_identical(tree, roundtrip)

    def test_to_dict(self):
        tree = DataTree.from_dict({"/a/b/c": None})
        roundtrip = DataTree.from_dict(tree.to_dict())
        assert_identical(tree, roundtrip)

        roundtrip = DataTree.from_dict(tree.to_dict(relative=True))
        assert_identical(tree, roundtrip)

        roundtrip = DataTree.from_dict(tree.children["a"].to_dict(relative=False))
        assert_identical(tree, roundtrip)

        expected = DataTree.from_dict({"b/c": None})
        actual = DataTree.from_dict(tree.children["a"].to_dict(relative=True))
        assert_identical(expected, actual)

    def test_roundtrip_unnamed_root(self, simple_datatree) -> None:
        # See GH81

        dt = simple_datatree
        dt.name = "root"
        roundtrip = DataTree.from_dict(dt.to_dict())
        assert roundtrip.equals(dt)

    def test_insertion_order(self) -> None:
        # regression test for GH issue #9276
        reversed = DataTree.from_dict(
            {
                "/Homer/Lisa": xr.Dataset({"age": 8}),
                "/Homer/Bart": xr.Dataset({"age": 10}),
                "/Homer": xr.Dataset({"age": 39}),
                "/": xr.Dataset({"age": 83}),
            }
        )
        expected = DataTree.from_dict(
            {
                "/": xr.Dataset({"age": 83}),
                "/Homer": xr.Dataset({"age": 39}),
                "/Homer/Lisa": xr.Dataset({"age": 8}),
                "/Homer/Bart": xr.Dataset({"age": 10}),
            }
        )
        assert reversed.equals(expected)

        # Check that Bart and Lisa's order is still preserved within the group,
        # despite 'Bart' coming before 'Lisa' when sorted alphabetically
        assert list(reversed["Homer"].children.keys()) == ["Lisa", "Bart"]

    def test_array_values(self) -> None:
        data = {"foo": xr.DataArray(1, name="bar")}
        with pytest.raises(TypeError):
            DataTree.from_dict(data)  # type: ignore[arg-type]

    def test_relative_paths(self) -> None:
        tree = DataTree.from_dict({".": None, "foo": None, "./bar": None, "x/y": None})
        paths = [node.path for node in tree.subtree]
        assert paths == [
            "/",
            "/foo",
            "/bar",
            "/x",
            "/x/y",
        ]

    def test_root_keys(self):
        ds = Dataset({"x": 1})
        expected = DataTree(dataset=ds)

        actual = DataTree.from_dict({"": ds})
        assert_identical(actual, expected)

        actual = DataTree.from_dict({".": ds})
        assert_identical(actual, expected)

        actual = DataTree.from_dict({"/": ds})
        assert_identical(actual, expected)

        actual = DataTree.from_dict({"./": ds})
        assert_identical(actual, expected)

        with pytest.raises(
            ValueError, match="multiple entries found corresponding to the root node"
        ):
            DataTree.from_dict({"": ds, "/": ds})

    def test_name(self):
        tree = DataTree.from_dict({"/": None}, name="foo")
        assert tree.name == "foo"

        tree = DataTree.from_dict({"/": DataTree()}, name="foo")
        assert tree.name == "foo"

        tree = DataTree.from_dict({"/": DataTree(name="bar")}, name="foo")
        assert tree.name == "foo"


class TestDatasetView:
    def test_view_contents(self) -> None:
        ds = create_test_data()
        dt = DataTree(dataset=ds)
        assert ds.identical(
            dt.dataset
        )  # this only works because Dataset.identical doesn't check types
        assert isinstance(dt.dataset, xr.Dataset)

    def test_immutability(self) -> None:
        # See issue https://github.com/xarray-contrib/datatree/issues/38
        dt = DataTree.from_dict(
            {
                "/": None,
                "/a": None,
            },
            name="root",
        )

        with pytest.raises(
            AttributeError, match="Mutation of the DatasetView is not allowed"
        ):
            dt.dataset["a"] = xr.DataArray(0)

        with pytest.raises(
            AttributeError, match="Mutation of the DatasetView is not allowed"
        ):
            dt.dataset.update({"a": 0})

        # TODO are there any other ways you can normally modify state (in-place)?
        # (not attribute-like assignment because that doesn't work on Dataset anyway)

    def test_methods(self) -> None:
        ds = create_test_data()
        dt = DataTree(dataset=ds)
        assert ds.mean().identical(dt.dataset.mean())
        assert isinstance(dt.dataset.mean(), xr.Dataset)

    def test_arithmetic(self, create_test_datatree) -> None:
        dt = create_test_datatree()
        expected = create_test_datatree(modify=lambda ds: 10.0 * ds)[
            "set1"
        ].to_dataset()
        result = 10.0 * dt["set1"].dataset
        assert result.identical(expected)

    def test_init_via_type(self) -> None:
        # from datatree GH issue https://github.com/xarray-contrib/datatree/issues/188
        # xarray's .weighted is unusual because it uses type() to create a Dataset/DataArray

        a = xr.DataArray(
            np.random.rand(3, 4, 10),
            dims=["x", "y", "time"],
            coords={"area": (["x", "y"], np.random.rand(3, 4))},
        ).to_dataset(name="data")
        dt = DataTree(dataset=a)

        def weighted_mean(ds):
            return ds.weighted(ds.area).mean(["x", "y"])

        weighted_mean(dt.dataset)

    def test_map_keep_attrs(self) -> None:
        # test DatasetView.map(..., keep_attrs=...)
        data = xr.DataArray([1, 2, 3], dims="x", attrs={"da": "attrs"})
        ds = xr.Dataset({"data": data}, attrs={"ds": "attrs"})
        dt = DataTree(ds)

        def func_keep(ds):
            # x.mean() removes the attrs of the data_vars
            return ds.map(lambda x: x.mean(), keep_attrs=True)

        result = xr.map_over_datasets(func_keep, dt)
        expected = dt.mean(keep_attrs=True)
        xr.testing.assert_identical(result, expected)

        # per default DatasetView.map does not keep attrs
        def func(ds):
            # x.mean() removes the attrs of the data_vars
            return ds.map(lambda x: x.mean())

        result = xr.map_over_datasets(func, dt)
        expected = dt.mean()
        xr.testing.assert_identical(result, expected.mean())


class TestAccess:
    def test_attribute_access(self, create_test_datatree) -> None:
        dt = create_test_datatree()

        # vars / coords
        for key in ["a", "set0"]:
            assert_equal(dt[key], getattr(dt, key))
            assert key in dir(dt)

        # dims
        assert_equal(dt["a"]["y"], dt.a.y)
        assert "y" in dir(dt["a"])

        # children
        for key in ["set1", "set2", "set3"]:
            assert_equal(dt[key], getattr(dt, key))
            assert key in dir(dt)

        # attrs
        dt.attrs["meta"] = "NASA"
        assert dt.attrs["meta"] == "NASA"
        assert "meta" in dir(dt)

    def test_ipython_key_completions_complex(self, create_test_datatree) -> None:
        dt = create_test_datatree()
        key_completions = dt._ipython_key_completions_()

        node_keys = [node.path[1:] for node in dt.descendants]
        assert all(node_key in key_completions for node_key in node_keys)

        var_keys = list(dt.variables.keys())
        assert all(var_key in key_completions for var_key in var_keys)

    def test_ipython_key_completitions_subnode(self) -> None:
        tree = xr.DataTree.from_dict({"/": None, "/a": None, "/a/b/": None})
        expected = ["b"]
        actual = tree["a"]._ipython_key_completions_()
        assert expected == actual

    def test_operation_with_attrs_but_no_data(self) -> None:
        # tests bug from xarray-datatree GH262
        xs = xr.Dataset({"testvar": xr.DataArray(np.ones((2, 3)))})
        dt = DataTree.from_dict({"node1": xs, "node2": xs})
        dt.attrs["test_key"] = 1  # sel works fine without this line
        dt.sel(dim_0=0)


class TestRepr:
    def test_repr_four_nodes(self) -> None:
        dt = DataTree.from_dict(
            {
                "/": xr.Dataset(
                    {"e": (("x",), [1.0, 2.0])},
                    coords={"x": [2.0, 3.0]},
                ),
                "/b": xr.Dataset({"f": (("y",), [3.0])}),
                "/b/c": xr.Dataset(),
                "/b/d": xr.Dataset({"g": 4.0}),
            }
        )

        result = repr(dt)
        expected = dedent(
            """
            <xarray.DataTree>
            Group: /
            │   Dimensions:  (x: 2)
            │   Coordinates:
            │     * x        (x) float64 16B 2.0 3.0
            │   Data variables:
            │       e        (x) float64 16B 1.0 2.0
            └── Group: /b
                │   Dimensions:  (y: 1)
                │   Dimensions without coordinates: y
                │   Data variables:
                │       f        (y) float64 8B 3.0
                ├── Group: /b/c
                └── Group: /b/d
                        Dimensions:  ()
                        Data variables:
                            g        float64 8B 4.0
            """
        ).strip()
        assert result == expected

        result = repr(dt.b)
        expected = dedent(
            """
            <xarray.DataTree 'b'>
            Group: /b
            │   Dimensions:  (x: 2, y: 1)
            │   Inherited coordinates:
            │     * x        (x) float64 16B 2.0 3.0
            │   Dimensions without coordinates: y
            │   Data variables:
            │       f        (y) float64 8B 3.0
            ├── Group: /b/c
            └── Group: /b/d
                    Dimensions:  ()
                    Data variables:
                        g        float64 8B 4.0
            """
        ).strip()
        assert result == expected

        result = repr(dt.b.d)
        expected = dedent(
            """
            <xarray.DataTree 'd'>
            Group: /b/d
                Dimensions:  (x: 2, y: 1)
                Inherited coordinates:
                  * x        (x) float64 16B 2.0 3.0
                Dimensions without coordinates: y
                Data variables:
                    g        float64 8B 4.0
            """
        ).strip()
        assert result == expected

    def test_repr_two_children(self) -> None:
        tree = DataTree.from_dict(
            {
                "/": Dataset(coords={"x": [1.0]}),
                "/first_child": None,
                "/second_child": Dataset({"foo": ("x", [0.0])}, coords={"z": 1.0}),
            }
        )

        result = repr(tree)
        expected = dedent(
            """
            <xarray.DataTree>
            Group: /
            │   Dimensions:  (x: 1)
            │   Coordinates:
            │     * x        (x) float64 8B 1.0
            ├── Group: /first_child
            └── Group: /second_child
                    Dimensions:  (x: 1)
                    Coordinates:
                        z        float64 8B 1.0
                    Data variables:
                        foo      (x) float64 8B 0.0
            """
        ).strip()
        assert result == expected

        result = repr(tree["first_child"])
        expected = dedent(
            """
            <xarray.DataTree 'first_child'>
            Group: /first_child
                Dimensions:  (x: 1)
                Inherited coordinates:
                  * x        (x) float64 8B 1.0
            """
        ).strip()
        assert result == expected

        result = repr(tree["second_child"])
        expected = dedent(
            """
            <xarray.DataTree 'second_child'>
            Group: /second_child
                Dimensions:  (x: 1)
                Coordinates:
                    z        float64 8B 1.0
                Inherited coordinates:
                  * x        (x) float64 8B 1.0
                Data variables:
                    foo      (x) float64 8B 0.0
            """
        ).strip()
        assert result == expected

    def test_repr_truncates_nodes(self) -> None:
        # construct a datatree with 50 nodes
        number_of_files = 10
        number_of_groups = 5
        tree_dict = {}
        for f in range(number_of_files):
            for g in range(number_of_groups):
                tree_dict[f"file_{f}/group_{g}"] = Dataset({"g": f * g})

        tree = DataTree.from_dict(tree_dict)
        with xr.set_options(display_max_children=3):
            result = repr(tree)

        expected = dedent(
            """
            <xarray.DataTree>
            Group: /
            ├── Group: /file_0
            │   ├── Group: /file_0/group_0
            │   │       Dimensions:  ()
            │   │       Data variables:
            │   │           g        int64 8B 0
            │   ├── Group: /file_0/group_1
            │   │       Dimensions:  ()
            │   │       Data variables:
            │   │           g        int64 8B 0
            │   ...
            │   └── Group: /file_0/group_4
            │           Dimensions:  ()
            │           Data variables:
            │               g        int64 8B 0
            ├── Group: /file_1
            │   ├── Group: /file_1/group_0
            │   │       Dimensions:  ()
            │   │       Data variables:
            │   │           g        int64 8B 0
            │   ├── Group: /file_1/group_1
            │   │       Dimensions:  ()
            │   │       Data variables:
            │   │           g        int64 8B 1
            │   ...
            │   └── Group: /file_1/group_4
            │           Dimensions:  ()
            │           Data variables:
            │               g        int64 8B 4
            ...
            └── Group: /file_9
                ├── Group: /file_9/group_0
                │       Dimensions:  ()
                │       Data variables:
                │           g        int64 8B 0
                ├── Group: /file_9/group_1
                │       Dimensions:  ()
                │       Data variables:
                │           g        int64 8B 9
                ...
                └── Group: /file_9/group_4
                        Dimensions:  ()
                        Data variables:
                            g        int64 8B 36
            """
        ).strip()
        assert expected == result

        with xr.set_options(display_max_children=10):
            result = repr(tree)

        for key in tree_dict:
            assert key in result

    def test_repr_inherited_dims(self) -> None:
        tree = DataTree.from_dict(
            {
                "/": Dataset({"foo": ("x", [1.0])}),
                "/child": Dataset({"bar": ("y", [2.0])}),
            }
        )

        result = repr(tree)
        expected = dedent(
            """
            <xarray.DataTree>
            Group: /
            │   Dimensions:  (x: 1)
            │   Dimensions without coordinates: x
            │   Data variables:
            │       foo      (x) float64 8B 1.0
            └── Group: /child
                    Dimensions:  (y: 1)
                    Dimensions without coordinates: y
                    Data variables:
                        bar      (y) float64 8B 2.0
            """
        ).strip()
        assert result == expected

        result = repr(tree["child"])
        expected = dedent(
            """
            <xarray.DataTree 'child'>
            Group: /child
                Dimensions:  (x: 1, y: 1)
                Dimensions without coordinates: x, y
                Data variables:
                    bar      (y) float64 8B 2.0
            """
        ).strip()
        assert result == expected

    @pytest.mark.skipif(
        ON_WINDOWS, reason="windows (pre NumPy2) uses int32 instead of int64"
    )
    def test_doc_example(self) -> None:
        # regression test for https://github.com/pydata/xarray/issues/9499
        time = xr.DataArray(
            data=np.array(["2022-01", "2023-01"], dtype="<U7"), dims="time"
        )
        stations = xr.DataArray(
            data=np.array(list("abcdef"), dtype="<U1"), dims="station"
        )
        lon = [-100, -80, -60]
        lat = [10, 20, 30]
        # Set up fake data
        wind_speed = xr.DataArray(np.ones((2, 6)) * 2, dims=("time", "station"))
        pressure = xr.DataArray(np.ones((2, 6)) * 3, dims=("time", "station"))
        air_temperature = xr.DataArray(np.ones((2, 6)) * 4, dims=("time", "station"))
        dewpoint = xr.DataArray(np.ones((2, 6)) * 5, dims=("time", "station"))
        infrared = xr.DataArray(np.ones((2, 3, 3)) * 6, dims=("time", "lon", "lat"))
        true_color = xr.DataArray(np.ones((2, 3, 3)) * 7, dims=("time", "lon", "lat"))
        tree = xr.DataTree.from_dict(
            {
                "/": xr.Dataset(
                    coords={"time": time},
                ),
                "/weather": xr.Dataset(
                    coords={"station": stations},
                    data_vars={
                        "wind_speed": wind_speed,
                        "pressure": pressure,
                    },
                ),
                "/weather/temperature": xr.Dataset(
                    data_vars={
                        "air_temperature": air_temperature,
                        "dewpoint": dewpoint,
                    },
                ),
                "/satellite": xr.Dataset(
                    coords={"lat": lat, "lon": lon},
                    data_vars={
                        "infrared": infrared,
                        "true_color": true_color,
                    },
                ),
            },
        )

        result = repr(tree)
        expected = dedent(
            """
            <xarray.DataTree>
            Group: /
            │   Dimensions:  (time: 2)
            │   Coordinates:
            │     * time     (time) <U7 56B '2022-01' '2023-01'
            ├── Group: /weather
            │   │   Dimensions:     (station: 6, time: 2)
            │   │   Coordinates:
            │   │     * station     (station) <U1 24B 'a' 'b' 'c' 'd' 'e' 'f'
            │   │   Data variables:
            │   │       wind_speed  (time, station) float64 96B 2.0 2.0 2.0 2.0 ... 2.0 2.0 2.0 2.0
            │   │       pressure    (time, station) float64 96B 3.0 3.0 3.0 3.0 ... 3.0 3.0 3.0 3.0
            │   └── Group: /weather/temperature
            │           Dimensions:          (time: 2, station: 6)
            │           Data variables:
            │               air_temperature  (time, station) float64 96B 4.0 4.0 4.0 4.0 ... 4.0 4.0 4.0
            │               dewpoint         (time, station) float64 96B 5.0 5.0 5.0 5.0 ... 5.0 5.0 5.0
            └── Group: /satellite
                    Dimensions:     (lat: 3, lon: 3, time: 2)
                    Coordinates:
                      * lat         (lat) int64 24B 10 20 30
                      * lon         (lon) int64 24B -100 -80 -60
                    Data variables:
                        infrared    (time, lon, lat) float64 144B 6.0 6.0 6.0 6.0 ... 6.0 6.0 6.0
                        true_color  (time, lon, lat) float64 144B 7.0 7.0 7.0 7.0 ... 7.0 7.0 7.0
            """
        ).strip()
        assert result == expected

        result = repr(tree["weather"])
        expected = dedent(
            """
            <xarray.DataTree 'weather'>
            Group: /weather
            │   Dimensions:     (time: 2, station: 6)
            │   Coordinates:
            │     * station     (station) <U1 24B 'a' 'b' 'c' 'd' 'e' 'f'
            │   Inherited coordinates:
            │     * time        (time) <U7 56B '2022-01' '2023-01'
            │   Data variables:
            │       wind_speed  (time, station) float64 96B 2.0 2.0 2.0 2.0 ... 2.0 2.0 2.0 2.0
            │       pressure    (time, station) float64 96B 3.0 3.0 3.0 3.0 ... 3.0 3.0 3.0 3.0
            └── Group: /weather/temperature
                    Dimensions:          (time: 2, station: 6)
                    Data variables:
                        air_temperature  (time, station) float64 96B 4.0 4.0 4.0 4.0 ... 4.0 4.0 4.0
                        dewpoint         (time, station) float64 96B 5.0 5.0 5.0 5.0 ... 5.0 5.0 5.0
            """
        ).strip()
        assert result == expected


def _exact_match(message: str) -> str:
    return re.escape(dedent(message).strip())


class TestInheritance:
    def test_inherited_dims(self) -> None:
        dt = DataTree.from_dict(
            {
                "/": xr.Dataset({"d": (("x",), [1, 2])}),
                "/b": xr.Dataset({"e": (("y",), [3])}),
                "/c": xr.Dataset({"f": (("y",), [3, 4, 5])}),
            }
        )
        assert dt.sizes == {"x": 2}
        # nodes should include inherited dimensions
        assert dt.b.sizes == {"x": 2, "y": 1}
        assert dt.c.sizes == {"x": 2, "y": 3}
        # dataset objects created from nodes should not
        assert dt.b.dataset.sizes == {"y": 1}
        assert dt.b.to_dataset(inherit=True).sizes == {"y": 1}
        assert dt.b.to_dataset(inherit=False).sizes == {"y": 1}

    def test_inherited_coords_index(self) -> None:
        dt = DataTree.from_dict(
            {
                "/": xr.Dataset({"d": (("x",), [1, 2])}, coords={"x": [2, 3]}),
                "/b": xr.Dataset({"e": (("y",), [3])}),
            }
        )
        assert "x" in dt["/b"].indexes
        assert "x" in dt["/b"].coords
        xr.testing.assert_identical(dt["/x"], dt["/b/x"])

    def test_inherit_only_index_coords(self) -> None:
        dt = DataTree.from_dict(
            {
                "/": xr.Dataset(coords={"x": [1], "y": 2}),
                "/b": xr.Dataset(coords={"z": 3}),
            }
        )
        assert dt.coords.keys() == {"x", "y"}
        xr.testing.assert_equal(
            dt["/x"], xr.DataArray([1], dims=["x"], coords={"x": [1], "y": 2})
        )
        xr.testing.assert_equal(dt["/y"], xr.DataArray(2, coords={"y": 2}))
        assert dt["/b"].coords.keys() == {"x", "z"}
        xr.testing.assert_equal(
            dt["/b/x"], xr.DataArray([1], dims=["x"], coords={"x": [1], "z": 3})
        )
        xr.testing.assert_equal(dt["/b/z"], xr.DataArray(3, coords={"z": 3}))

    def test_inherited_coords_with_index_are_deduplicated(self) -> None:
        dt = DataTree.from_dict(
            {
                "/": xr.Dataset(coords={"x": [1, 2]}),
                "/b": xr.Dataset(coords={"x": [1, 2]}),
            }
        )
        child_dataset = dt.children["b"].to_dataset(inherit=False)
        expected = xr.Dataset()
        assert_identical(child_dataset, expected)

        dt["/c"] = xr.Dataset({"foo": ("x", [4, 5])}, coords={"x": [1, 2]})
        child_dataset = dt.children["c"].to_dataset(inherit=False)
        expected = xr.Dataset({"foo": ("x", [4, 5])})
        assert_identical(child_dataset, expected)

    def test_deduplicated_after_setitem(self) -> None:
        # regression test for GH #9601
        dt = DataTree.from_dict(
            {
                "/": xr.Dataset(coords={"x": [1, 2]}),
                "/b": None,
            }
        )
        dt["b/x"] = dt["x"]
        child_dataset = dt.children["b"].to_dataset(inherit=False)
        expected = xr.Dataset()
        assert_identical(child_dataset, expected)

    def test_inconsistent_dims(self) -> None:
        expected_msg = _exact_match(
            """
            group '/b' is not aligned with its parents:
            Group:
                Dimensions:  (x: 1)
                Dimensions without coordinates: x
                Data variables:
                    c        (x) float64 8B 3.0
            From parents:
                Dimensions:  (x: 2)
                Dimensions without coordinates: x
            """
        )

        with pytest.raises(ValueError, match=expected_msg):
            DataTree.from_dict(
                {
                    "/": xr.Dataset({"a": (("x",), [1.0, 2.0])}),
                    "/b": xr.Dataset({"c": (("x",), [3.0])}),
                }
            )

        dt = DataTree()
        dt["/a"] = xr.DataArray([1.0, 2.0], dims=["x"])
        with pytest.raises(ValueError, match=expected_msg):
            dt["/b/c"] = xr.DataArray([3.0], dims=["x"])

        b = DataTree(dataset=xr.Dataset({"c": (("x",), [3.0])}))
        with pytest.raises(ValueError, match=expected_msg):
            DataTree(
                dataset=xr.Dataset({"a": (("x",), [1.0, 2.0])}),
                children={"b": b},
            )

    def test_inconsistent_child_indexes(self) -> None:
        expected_msg = _exact_match(
            """
            group '/b' is not aligned with its parents:
            Group:
                Dimensions:  (x: 1)
                Coordinates:
                  * x        (x) float64 8B 2.0
                Data variables:
                    *empty*
            From parents:
                Dimensions:  (x: 1)
                Coordinates:
                  * x        (x) float64 8B 1.0
            """
        )

        with pytest.raises(ValueError, match=expected_msg):
            DataTree.from_dict(
                {
                    "/": xr.Dataset(coords={"x": [1.0]}),
                    "/b": xr.Dataset(coords={"x": [2.0]}),
                }
            )

        dt = DataTree()
        dt.dataset = xr.Dataset(coords={"x": [1.0]})  # type: ignore[assignment,unused-ignore]
        dt["/b"] = DataTree()
        with pytest.raises(ValueError, match=expected_msg):
            dt["/b"].dataset = xr.Dataset(coords={"x": [2.0]})

        b = DataTree(xr.Dataset(coords={"x": [2.0]}))
        with pytest.raises(ValueError, match=expected_msg):
            DataTree(dataset=xr.Dataset(coords={"x": [1.0]}), children={"b": b})

    def test_inconsistent_grandchild_indexes(self) -> None:
        expected_msg = _exact_match(
            """
            group '/b/c' is not aligned with its parents:
            Group:
                Dimensions:  (x: 1)
                Coordinates:
                  * x        (x) float64 8B 2.0
                Data variables:
                    *empty*
            From parents:
                Dimensions:  (x: 1)
                Coordinates:
                  * x        (x) float64 8B 1.0
            """
        )

        with pytest.raises(ValueError, match=expected_msg):
            DataTree.from_dict(
                {
                    "/": xr.Dataset(coords={"x": [1.0]}),
                    "/b/c": xr.Dataset(coords={"x": [2.0]}),
                }
            )

        dt = DataTree()
        dt.dataset = xr.Dataset(coords={"x": [1.0]})  # type: ignore[assignment,unused-ignore]
        dt["/b/c"] = DataTree()
        with pytest.raises(ValueError, match=expected_msg):
            dt["/b/c"].dataset = xr.Dataset(coords={"x": [2.0]})

        c = DataTree(xr.Dataset(coords={"x": [2.0]}))
        b = DataTree(children={"c": c})
        with pytest.raises(ValueError, match=expected_msg):
            DataTree(dataset=xr.Dataset(coords={"x": [1.0]}), children={"b": b})

    def test_inconsistent_grandchild_dims(self) -> None:
        expected_msg = _exact_match(
            """
            group '/b/c' is not aligned with its parents:
            Group:
                Dimensions:  (x: 1)
                Dimensions without coordinates: x
                Data variables:
                    d        (x) float64 8B 3.0
            From parents:
                Dimensions:  (x: 2)
                Dimensions without coordinates: x
            """
        )

        with pytest.raises(ValueError, match=expected_msg):
            DataTree.from_dict(
                {
                    "/": xr.Dataset({"a": (("x",), [1.0, 2.0])}),
                    "/b/c": xr.Dataset({"d": (("x",), [3.0])}),
                }
            )

        dt = DataTree()
        dt["/a"] = xr.DataArray([1.0, 2.0], dims=["x"])
        with pytest.raises(ValueError, match=expected_msg):
            dt["/b/c/d"] = xr.DataArray([3.0], dims=["x"])


class TestRestructuring:
    def test_drop_nodes(self) -> None:
        sue = DataTree.from_dict({"Mary": None, "Kate": None, "Ashley": None})

        # test drop just one node
        dropped_one = sue.drop_nodes(names="Mary")
        assert "Mary" not in dropped_one.children

        # test drop multiple nodes
        dropped = sue.drop_nodes(names=["Mary", "Kate"])
        assert not {"Mary", "Kate"}.intersection(set(dropped.children))
        assert "Ashley" in dropped.children

        # test raise
        with pytest.raises(KeyError, match="nodes {'Mary'} not present"):
            dropped.drop_nodes(names=["Mary", "Ashley"])

        # test ignore
        childless = dropped.drop_nodes(names=["Mary", "Ashley"], errors="ignore")
        assert childless.children == {}

    def test_assign(self) -> None:
        dt = DataTree()
        expected = DataTree.from_dict({"/": xr.Dataset({"foo": 0}), "/a": None})

        # kwargs form
        result = dt.assign(foo=xr.DataArray(0), a=DataTree())
        assert_equal(result, expected)

        # dict form
        result = dt.assign({"foo": xr.DataArray(0), "a": DataTree()})
        assert_equal(result, expected)

    def test_filter_like(self) -> None:
        flower_tree = DataTree.from_dict(
            {"root": None, "trunk": None, "leaves": None, "flowers": None}
        )
        fruit_tree = DataTree.from_dict(
            {"root": None, "trunk": None, "leaves": None, "fruit": None}
        )
        barren_tree = DataTree.from_dict({"root": None, "trunk": None})

        # test filter_like tree
        filtered_tree = flower_tree.filter_like(barren_tree)

        assert filtered_tree.equals(barren_tree)
        assert "flowers" not in filtered_tree.children

        # test symmetrical pruning results in isomorphic trees
        assert flower_tree.filter_like(fruit_tree).isomorphic(
            fruit_tree.filter_like(flower_tree)
        )

        # test "deep" pruning
        dt = DataTree.from_dict(
            {"/a/A": None, "/a/B": None, "/b/A": None, "/b/B": None}
        )
        other = DataTree.from_dict({"/a/A": None, "/b/A": None})

        filtered = dt.filter_like(other)
        assert filtered.equals(other)


class TestPipe:
    def test_noop(self, create_test_datatree: Callable[[], DataTree]) -> None:
        dt = create_test_datatree()

        actual = dt.pipe(lambda tree: tree)
        assert actual.identical(dt)

    def test_args(self, create_test_datatree: Callable[[], DataTree]) -> None:
        dt = create_test_datatree()

        def f(tree: DataTree, x: int, y: int) -> DataTree:
            return tree.assign(
                arr_with_attrs=xr.Variable("dim0", [], attrs=dict(x=x, y=y))
            )

        actual = dt.pipe(f, 1, 2)
        assert actual["arr_with_attrs"].attrs == dict(x=1, y=2)

    def test_kwargs(self, create_test_datatree: Callable[[], DataTree]) -> None:
        dt = create_test_datatree()

        def f(tree: DataTree, *, x: int, y: int, z: int) -> DataTree:
            return tree.assign(
                arr_with_attrs=xr.Variable("dim0", [], attrs=dict(x=x, y=y, z=z))
            )

        attrs = {"x": 1, "y": 2, "z": 3}

        actual = dt.pipe(f, **attrs)
        assert actual["arr_with_attrs"].attrs == attrs

    def test_args_kwargs(self, create_test_datatree: Callable[[], DataTree]) -> None:
        dt = create_test_datatree()

        def f(tree: DataTree, x: int, *, y: int, z: int) -> DataTree:
            return tree.assign(
                arr_with_attrs=xr.Variable("dim0", [], attrs=dict(x=x, y=y, z=z))
            )

        attrs = {"x": 1, "y": 2, "z": 3}

        actual = dt.pipe(f, attrs["x"], y=attrs["y"], z=attrs["z"])
        assert actual["arr_with_attrs"].attrs == attrs

    def test_named_self(self, create_test_datatree: Callable[[], DataTree]) -> None:
        dt = create_test_datatree()

        def f(x: int, tree: DataTree, y: int):
            tree.attrs.update({"x": x, "y": y})
            return tree

        attrs = {"x": 1, "y": 2}

        actual = dt.pipe((f, "tree"), **attrs)

        assert actual is dt and actual.attrs == attrs


class TestIsomorphicEqualsAndIdentical:
    def test_isomorphic(self):
        tree = DataTree.from_dict({"/a": None, "/a/b": None, "/c": None})

        diff_data = DataTree.from_dict(
            {"/a": None, "/a/b": None, "/c": xr.Dataset({"foo": 1})}
        )
        assert tree.isomorphic(diff_data)

        diff_order = DataTree.from_dict({"/c": None, "/a": None, "/a/b": None})
        assert tree.isomorphic(diff_order)

        diff_nodes = DataTree.from_dict({"/a": None, "/a/b": None, "/d": None})
        assert not tree.isomorphic(diff_nodes)

        more_nodes = DataTree.from_dict(
            {"/a": None, "/a/b": None, "/c": None, "/d": None}
        )
        assert not tree.isomorphic(more_nodes)

    def test_minimal_variations(self):
        tree = DataTree.from_dict(
            {
                "/": Dataset({"x": 1}),
                "/child": Dataset({"x": 2}),
            }
        )
        assert tree.equals(tree)
        assert tree.identical(tree)

        child = tree.children["child"]
        assert child.equals(child)
        assert child.identical(child)

        new_child = DataTree(dataset=Dataset({"x": 2}), name="child")
        assert child.equals(new_child)
        assert child.identical(new_child)

        anonymous_child = DataTree(dataset=Dataset({"x": 2}))
        # TODO: re-enable this after fixing .equals() not to require matching
        # names on the root node (i.e., after switching to use zip_subtrees)
        # assert child.equals(anonymous_child)
        assert not child.identical(anonymous_child)

        different_variables = DataTree.from_dict(
            {
                "/": Dataset(),
                "/other": Dataset({"x": 2}),
            }
        )
        assert not tree.equals(different_variables)
        assert not tree.identical(different_variables)

        different_root_data = DataTree.from_dict(
            {
                "/": Dataset({"x": 4}),
                "/child": Dataset({"x": 2}),
            }
        )
        assert not tree.equals(different_root_data)
        assert not tree.identical(different_root_data)

        different_child_data = DataTree.from_dict(
            {
                "/": Dataset({"x": 1}),
                "/child": Dataset({"x": 3}),
            }
        )
        assert not tree.equals(different_child_data)
        assert not tree.identical(different_child_data)

        different_child_node_attrs = DataTree.from_dict(
            {
                "/": Dataset({"x": 1}),
                "/child": Dataset({"x": 2}, attrs={"foo": "bar"}),
            }
        )
        assert tree.equals(different_child_node_attrs)
        assert not tree.identical(different_child_node_attrs)

        different_child_variable_attrs = DataTree.from_dict(
            {
                "/": Dataset({"x": 1}),
                "/child": Dataset({"x": ((), 2, {"foo": "bar"})}),
            }
        )
        assert tree.equals(different_child_variable_attrs)
        assert not tree.identical(different_child_variable_attrs)

        different_name = DataTree.from_dict(
            {
                "/": Dataset({"x": 1}),
                "/child": Dataset({"x": 2}),
            },
            name="different",
        )
        # TODO: re-enable this after fixing .equals() not to require matching
        # names on the root node (i.e., after switching to use zip_subtrees)
        # assert tree.equals(different_name)
        assert not tree.identical(different_name)

    def test_differently_inherited_coordinates(self):
        root = DataTree.from_dict(
            {
                "/": Dataset(coords={"x": [1, 2]}),
                "/child": Dataset(),
            }
        )
        child = root.children["child"]
        assert child.equals(child)
        assert child.identical(child)

        new_child = DataTree(dataset=Dataset(coords={"x": [1, 2]}), name="child")
        assert child.equals(new_child)
        assert not child.identical(new_child)

        deeper_root = DataTree(children={"root": root})
        grandchild = deeper_root.children["root"].children["child"]
        assert child.equals(grandchild)
        assert child.identical(grandchild)


class TestSubset:
    def test_match(self) -> None:
        # TODO is this example going to cause problems with case sensitivity?
        dt = DataTree.from_dict(
            {
                "/a/A": None,
                "/a/B": None,
                "/b/A": None,
                "/b/B": None,
            }
        )
        result = dt.match("*/B")
        expected = DataTree.from_dict(
            {
                "/a/B": None,
                "/b/B": None,
            }
        )
        assert_identical(result, expected)

        result = dt.children["a"].match("B")
        expected = DataTree.from_dict({"/B": None}, name="a")
        assert_identical(result, expected)

    def test_filter(self) -> None:
        simpsons = DataTree.from_dict(
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
        expected = DataTree.from_dict(
            {
                "/": xr.Dataset({"age": 83}),
                "/Herbert": xr.Dataset({"age": 40}),
                "/Homer": xr.Dataset({"age": 39}),
            },
            name="Abe",
        )
        elders = simpsons.filter(lambda node: node["age"].item() > 18)
        assert_identical(elders, expected)

        expected = DataTree.from_dict({"/Bart": xr.Dataset({"age": 10})}, name="Homer")
        actual = simpsons.children["Homer"].filter(
            lambda node: node["age"].item() == 10
        )
        assert_identical(actual, expected)

    def test_prune_basic(self) -> None:
        tree = DataTree.from_dict(
            {"/a": xr.Dataset({"foo": ("x", [1, 2])}), "/b": xr.Dataset()}
        )

        pruned = tree.prune()

        assert "a" in pruned.children
        assert "b" not in pruned.children
        assert_identical(
            pruned.children["a"].to_dataset(), tree.children["a"].to_dataset()
        )

    def test_prune_with_zero_size_vars(self) -> None:
        tree = DataTree.from_dict(
            {
                "/a": xr.Dataset({"foo": ("x", [1, 2])}),
                "/b": xr.Dataset({"empty": ("dim", [])}),
                "/c": xr.Dataset(),
            }
        )

        pruned_default = tree.prune()
        expected_default = DataTree.from_dict(
            {
                "/a": xr.Dataset({"foo": ("x", [1, 2])}),
                "/b": xr.Dataset({"empty": ("dim", [])}),
            }
        )
        assert_identical(pruned_default, expected_default)

        pruned_strict = tree.prune(drop_size_zero_vars=True)
        expected_strict = DataTree.from_dict(
            {
                "/a": xr.Dataset({"foo": ("x", [1, 2])}),
            }
        )
        assert_identical(pruned_strict, expected_strict)

    def test_prune_with_intermediate_nodes(self) -> None:
        tree = DataTree.from_dict(
            {
                "/": xr.Dataset(),
                "/group1": xr.Dataset(),
                "/group1/subA": xr.Dataset({"temp": ("x", [1, 2])}),
                "/group1/subB": xr.Dataset(),
                "/group2": xr.Dataset({"empty": ("dim", [])}),
            }
        )
        pruned = tree.prune()
        expected_tree = DataTree.from_dict(
            {
                "/group1/subA": xr.Dataset({"temp": ("x", [1, 2])}),
                "/group2": xr.Dataset({"empty": ("dim", [])}),
            }
        )
        assert_identical(pruned, expected_tree)

    def test_prune_after_filtering(self) -> None:
        from pandas import date_range

        ds1 = xr.Dataset(
            {"foo": ("time", [1, 2, 3, 4, 5])},
            coords={"time": date_range("2023-01-01", periods=5, freq="D")},
        )
        ds2 = xr.Dataset(
            {"var": ("time", [1, 2, 3, 4, 5])},
            coords={"time": date_range("2023-01-04", periods=5, freq="D")},
        )

        tree = DataTree.from_dict({"a": ds1, "b": ds2})
        filtered = tree.sel(time=slice("2023-01-01", "2023-01-03"))

        pruned = filtered.prune(drop_size_zero_vars=True)
        expected_tree = DataTree.from_dict(
            {"a": ds1.sel(time=slice("2023-01-01", "2023-01-03"))}
        )
        assert_identical(pruned, expected_tree)


class TestIndexing:
    def test_isel_siblings(self) -> None:
        tree = DataTree.from_dict(
            {
                "/first": xr.Dataset({"a": ("x", [1, 2])}),
                "/second": xr.Dataset({"b": ("x", [1, 2, 3])}),
            }
        )

        expected = DataTree.from_dict(
            {
                "/first": xr.Dataset({"a": 2}),
                "/second": xr.Dataset({"b": 3}),
            }
        )
        actual = tree.isel(x=-1)
        assert_identical(actual, expected)

        expected = DataTree.from_dict(
            {
                "/first": xr.Dataset({"a": ("x", [1])}),
                "/second": xr.Dataset({"b": ("x", [1])}),
            }
        )
        actual = tree.isel(x=slice(1))
        assert_identical(actual, expected)

        actual = tree.isel(x=[0])
        assert_identical(actual, expected)

        actual = tree.isel(x=slice(None))
        assert_identical(actual, tree)

    def test_isel_inherited(self) -> None:
        tree = DataTree.from_dict(
            {
                "/": xr.Dataset(coords={"x": [1, 2]}),
                "/child": xr.Dataset({"foo": ("x", [3, 4])}),
            }
        )

        expected = DataTree.from_dict(
            {
                "/": xr.Dataset(coords={"x": 2}),
                "/child": xr.Dataset({"foo": 4}),
            }
        )
        actual = tree.isel(x=-1)
        assert_identical(actual, expected)

        expected = DataTree.from_dict(
            {
                "/child": xr.Dataset({"foo": 4}),
            }
        )
        actual = tree.isel(x=-1, drop=True)
        assert_identical(actual, expected)

        expected = DataTree.from_dict(
            {
                "/": xr.Dataset(coords={"x": [1]}),
                "/child": xr.Dataset({"foo": ("x", [3])}),
            }
        )
        actual = tree.isel(x=[0])
        assert_identical(actual, expected)

        actual = tree.isel(x=slice(None))

        # TODO: re-enable after the fix to copy() from #9628 is submitted
        # actual = tree.children["child"].isel(x=slice(None))
        # expected = tree.children["child"].copy()
        # assert_identical(actual, expected)

        actual = tree.children["child"].isel(x=0)
        expected = DataTree(
            dataset=xr.Dataset({"foo": 3}, coords={"x": 1}),
            name="child",
        )
        assert_identical(actual, expected)

    def test_sel(self) -> None:
        tree = DataTree.from_dict(
            {
                "/first": xr.Dataset({"a": ("x", [1, 2, 3])}, coords={"x": [1, 2, 3]}),
                "/second": xr.Dataset({"b": ("x", [4, 5])}, coords={"x": [2, 3]}),
            }
        )
        expected = DataTree.from_dict(
            {
                "/first": xr.Dataset({"a": 2}, coords={"x": 2}),
                "/second": xr.Dataset({"b": 4}, coords={"x": 2}),
            }
        )
        actual = tree.sel(x=2)
        assert_identical(actual, expected)

        actual = tree.children["first"].sel(x=2)
        expected = DataTree(
            dataset=xr.Dataset({"a": 2}, coords={"x": 2}),
            name="first",
        )
        assert_identical(actual, expected)

    def test_sel_isel_error_has_node_info(self) -> None:
        tree = DataTree.from_dict(
            {
                "/first": xr.Dataset({"a": ("x", [1, 2, 3])}, coords={"x": [1, 2, 3]}),
                "/second": xr.Dataset({"b": ("x", [4, 5])}, coords={"x": [2, 3]}),
            }
        )

        with pytest.raises(
            KeyError,
            match="Raised whilst mapping function over node with path 'second'",
        ):
            tree.sel(x=1)

        with pytest.raises(
            IndexError,
            match="Raised whilst mapping function over node with path 'first'",
        ):
            tree.isel(x=4)


class TestAggregations:
    def test_reduce_method(self) -> None:
        ds = xr.Dataset({"a": ("x", [False, True, False])})
        dt = DataTree.from_dict({"/": ds, "/results": ds})

        expected = DataTree.from_dict({"/": ds.any(), "/results": ds.any()})

        result = dt.any()
        assert_equal(result, expected)

    def test_nan_reduce_method(self) -> None:
        ds = xr.Dataset({"a": ("x", [1, 2, 3])})
        dt = DataTree.from_dict({"/": ds, "/results": ds})

        expected = DataTree.from_dict({"/": ds.mean(), "/results": ds.mean()})

        result = dt.mean()
        assert_equal(result, expected)

    def test_cum_method(self) -> None:
        ds = xr.Dataset({"a": ("x", [1, 2, 3])})
        dt = DataTree.from_dict({"/": ds, "/results": ds})

        expected = DataTree.from_dict(
            {
                "/": ds.cumsum(),
                "/results": ds.cumsum(),
            }
        )

        result = dt.cumsum()
        assert_equal(result, expected)

    def test_dim_argument(self) -> None:
        dt = DataTree.from_dict(
            {
                "/a": xr.Dataset({"A": ("x", [1, 2])}),
                "/b": xr.Dataset({"B": ("y", [1, 2])}),
            }
        )

        expected = DataTree.from_dict(
            {
                "/a": xr.Dataset({"A": 1.5}),
                "/b": xr.Dataset({"B": 1.5}),
            }
        )
        actual = dt.mean()
        assert_equal(expected, actual)

        actual = dt.mean(dim=...)
        assert_equal(expected, actual)

        expected = DataTree.from_dict(
            {
                "/a": xr.Dataset({"A": 1.5}),
                "/b": xr.Dataset({"B": ("y", [1.0, 2.0])}),
            }
        )
        actual = dt.mean("x")
        assert_equal(expected, actual)

        with pytest.raises(
            ValueError,
            match=re.escape("Dimension(s) 'invalid' do not exist."),
        ):
            dt.mean("invalid")

    def test_subtree(self) -> None:
        tree = DataTree.from_dict(
            {
                "/child": Dataset({"a": ("x", [1, 2])}),
            }
        )
        expected = DataTree(dataset=Dataset({"a": 1.5}), name="child")
        actual = tree.children["child"].mean()
        assert_identical(expected, actual)


class TestOps:
    def test_unary_op(self) -> None:
        ds1 = xr.Dataset({"a": [5], "b": [3]})
        ds2 = xr.Dataset({"x": [0.1, 0.2], "y": [10, 20]})
        dt = DataTree.from_dict({"/": ds1, "/subnode": ds2})

        expected = DataTree.from_dict({"/": (-ds1), "/subnode": (-ds2)})

        result = -dt
        assert_equal(result, expected)

    def test_unary_op_inherited_coords(self) -> None:
        tree = DataTree(xr.Dataset(coords={"x": [1, 2, 3]}))
        tree["/foo"] = DataTree(xr.Dataset({"bar": ("x", [4, 5, 6])}))
        actual = -tree

        actual_dataset = actual.children["foo"].to_dataset(inherit=False)
        assert "x" not in actual_dataset.coords

        expected = tree.copy()
        # unary ops are not applied to coordinate variables, only data variables
        expected["/foo/bar"].data = np.array([-4, -5, -6])
        assert_identical(actual, expected)

    def test_binary_op_on_int(self) -> None:
        ds1 = xr.Dataset({"a": [5], "b": [3]})
        ds2 = xr.Dataset({"x": [0.1, 0.2], "y": [10, 20]})
        dt = DataTree.from_dict({"/": ds1, "/subnode": ds2})

        expected = DataTree.from_dict({"/": ds1 * 5, "/subnode": ds2 * 5})

        result = dt * 5
        assert_equal(result, expected)

    def test_binary_op_on_dataarray(self) -> None:
        ds1 = xr.Dataset({"a": [5], "b": [3]})
        ds2 = xr.Dataset({"x": [0.1, 0.2], "y": [10, 20]})
        dt = DataTree.from_dict(
            {
                "/": ds1,
                "/subnode": ds2,
            }
        )

        other_da = xr.DataArray(name="z", data=[0.1, 0.2], dims="z")

        expected = DataTree.from_dict(
            {
                "/": ds1 * other_da,
                "/subnode": ds2 * other_da,
            }
        )

        result = dt * other_da
        assert_equal(result, expected)

    def test_binary_op_on_dataset(self) -> None:
        ds1 = xr.Dataset({"a": [5], "b": [3]})
        ds2 = xr.Dataset({"x": [0.1, 0.2], "y": [10, 20]})
        dt = DataTree.from_dict(
            {
                "/": ds1,
                "/subnode": ds2,
            }
        )

        other_ds = xr.Dataset({"z": ("z", [0.1, 0.2])})

        expected = DataTree.from_dict(
            {
                "/": ds1 * other_ds,
                "/subnode": ds2 * other_ds,
            }
        )

        result = dt * other_ds
        assert_equal(result, expected)

    def test_binary_op_on_datatree(self) -> None:
        ds1 = xr.Dataset({"a": [5], "b": [3]})
        ds2 = xr.Dataset({"x": [0.1, 0.2], "y": [10, 20]})

        dt = DataTree.from_dict({"/": ds1, "/subnode": ds2})

        expected = DataTree.from_dict({"/": ds1 * ds1, "/subnode": ds2 * ds2})

        result = dt * dt
        assert_equal(result, expected)

    def test_binary_op_order_invariant(self) -> None:
        tree_ab = DataTree.from_dict({"/a": Dataset({"a": 1}), "/b": Dataset({"b": 2})})
        tree_ba = DataTree.from_dict({"/b": Dataset({"b": 2}), "/a": Dataset({"a": 1})})
        expected = DataTree.from_dict(
            {"/a": Dataset({"a": 2}), "/b": Dataset({"b": 4})}
        )
        actual = tree_ab + tree_ba
        assert_identical(expected, actual)

    def test_arithmetic_inherited_coords(self) -> None:
        tree = DataTree(xr.Dataset(coords={"x": [1, 2, 3]}))
        tree["/foo"] = DataTree(xr.Dataset({"bar": ("x", [4, 5, 6])}))
        actual = 2 * tree

        actual_dataset = actual.children["foo"].to_dataset(inherit=False)
        assert "x" not in actual_dataset.coords

        expected = tree.copy()
        expected["/foo/bar"].data = np.array([8, 10, 12])
        assert_identical(actual, expected)

    def test_binary_op_commutativity_with_dataset(self) -> None:
        # regression test for #9365

        ds1 = xr.Dataset({"a": [5], "b": [3]})
        ds2 = xr.Dataset({"x": [0.1, 0.2], "y": [10, 20]})
        dt = DataTree.from_dict(
            {
                "/": ds1,
                "/subnode": ds2,
            }
        )

        other_ds = xr.Dataset({"z": ("z", [0.1, 0.2])})

        expected = DataTree.from_dict(
            {
                "/": ds1 * other_ds,
                "/subnode": ds2 * other_ds,
            }
        )

        result = other_ds * dt
        assert_equal(result, expected)

    def test_inplace_binary_op(self) -> None:
        ds1 = xr.Dataset({"a": [5], "b": [3]})
        ds2 = xr.Dataset({"x": [0.1, 0.2], "y": [10, 20]})
        dt = DataTree.from_dict({"/": ds1, "/subnode": ds2})

        expected = DataTree.from_dict({"/": ds1 + 1, "/subnode": ds2 + 1})

        dt += 1
        assert_equal(dt, expected)

    def test_dont_broadcast_single_node_tree(self) -> None:
        # regression test for https://github.com/pydata/xarray/issues/9365#issuecomment-2291622577
        ds1 = xr.Dataset({"a": [5], "b": [3]})
        ds2 = xr.Dataset({"x": [0.1, 0.2], "y": [10, 20]})
        dt = DataTree.from_dict({"/": ds1, "/subnode": ds2})
        node = dt["/subnode"]

        with pytest.raises(
            xr.TreeIsomorphismError,
            match=re.escape(r"children at root node do not match: ['subnode'] vs []"),
        ):
            dt * node


class TestUFuncs:
    @pytest.mark.xfail(reason="__array_ufunc__ not implemented yet")
    def test_tree(self, create_test_datatree):
        dt = create_test_datatree()
        expected = create_test_datatree(modify=np.sin)
        result_tree = np.sin(dt)
        assert_equal(result_tree, expected)


class Closer:
    def __init__(self):
        self.closed = False

    def close(self):
        if self.closed:
            raise RuntimeError("already closed")
        self.closed = True


@pytest.fixture
def tree_and_closers():
    tree = DataTree.from_dict({"/child/grandchild": None})
    closers = {
        "/": Closer(),
        "/child": Closer(),
        "/child/grandchild": Closer(),
    }
    for path, closer in closers.items():
        tree[path].set_close(closer.close)
    return tree, closers


class TestClose:
    def test_close(self, tree_and_closers):
        tree, closers = tree_and_closers
        assert not any(closer.closed for closer in closers.values())
        tree.close()
        assert all(closer.closed for closer in closers.values())
        tree.close()  # should not error

    def test_context_manager(self, tree_and_closers):
        tree, closers = tree_and_closers
        assert not any(closer.closed for closer in closers.values())
        with tree:
            pass
        assert all(closer.closed for closer in closers.values())

    def test_close_child(self, tree_and_closers):
        tree, closers = tree_and_closers
        assert not any(closer.closed for closer in closers.values())
        tree["child"].close()  # should only close descendants
        assert not closers["/"].closed
        assert closers["/child"].closed
        assert closers["/child/grandchild"].closed

    def test_close_datasetview(self, tree_and_closers):
        tree, _ = tree_and_closers

        with pytest.raises(
            AttributeError,
            match=re.escape(
                r"cannot close a DatasetView(). Close the associated DataTree node instead"
            ),
        ):
            tree.dataset.close()

        with pytest.raises(
            AttributeError, match=re.escape(r"cannot modify a DatasetView()")
        ):
            tree.dataset.set_close(None)

    def test_close_dataset(self, tree_and_closers):
        tree, closers = tree_and_closers
        ds = tree.to_dataset()  # should discard closers
        ds.close()
        assert not closers["/"].closed

    # with tree:
    #     pass


@requires_dask
class TestDask:
    def test_chunksizes(self):
        ds1 = xr.Dataset({"a": ("x", np.arange(10))})
        ds2 = xr.Dataset({"b": ("y", np.arange(5))})
        ds3 = xr.Dataset({"c": ("z", np.arange(4))})
        ds4 = xr.Dataset({"d": ("x", np.arange(-5, 5))})

        groups = {
            "/": ds1.chunk({"x": 5}),
            "/group1": ds2.chunk({"y": 3}),
            "/group2": ds3.chunk({"z": 2}),
            "/group1/subgroup1": ds4.chunk({"x": 5}),
        }

        tree = xr.DataTree.from_dict(groups)

        expected_chunksizes = {path: node.chunksizes for path, node in groups.items()}

        assert tree.chunksizes == expected_chunksizes

    def test_load(self):
        ds1 = xr.Dataset({"a": ("x", np.arange(10))})
        ds2 = xr.Dataset({"b": ("y", np.arange(5))})
        ds3 = xr.Dataset({"c": ("z", np.arange(4))})
        ds4 = xr.Dataset({"d": ("x", np.arange(-5, 5))})

        groups = {"/": ds1, "/group1": ds2, "/group2": ds3, "/group1/subgroup1": ds4}

        expected = xr.DataTree.from_dict(groups)
        tree = xr.DataTree.from_dict(
            {
                "/": ds1.chunk({"x": 5}),
                "/group1": ds2.chunk({"y": 3}),
                "/group2": ds3.chunk({"z": 2}),
                "/group1/subgroup1": ds4.chunk({"x": 5}),
            }
        )
        expected_chunksizes: Mapping[str, Mapping]
        expected_chunksizes = {node.path: {} for node in tree.subtree}
        actual = tree.load()

        assert_identical(actual, expected)
        assert tree.chunksizes == expected_chunksizes
        assert actual.chunksizes == expected_chunksizes

        tree = xr.DataTree.from_dict(groups)
        actual = tree.load()
        assert_identical(actual, expected)
        assert actual.chunksizes == expected_chunksizes

    def test_compute(self):
        ds1 = xr.Dataset({"a": ("x", np.arange(10))})
        ds2 = xr.Dataset({"b": ("y", np.arange(5))})
        ds3 = xr.Dataset({"c": ("z", np.arange(4))})
        ds4 = xr.Dataset({"d": ("x", np.arange(-5, 5))})

        expected = xr.DataTree.from_dict(
            {"/": ds1, "/group1": ds2, "/group2": ds3, "/group1/subgroup1": ds4}
        )
        tree = xr.DataTree.from_dict(
            {
                "/": ds1.chunk({"x": 5}),
                "/group1": ds2.chunk({"y": 3}),
                "/group2": ds3.chunk({"z": 2}),
                "/group1/subgroup1": ds4.chunk({"x": 5}),
            }
        )
        original_chunksizes = tree.chunksizes
        expected_chunksizes: Mapping[str, Mapping]
        expected_chunksizes = {node.path: {} for node in tree.subtree}
        actual = tree.compute()

        assert_identical(actual, expected)

        assert actual.chunksizes == expected_chunksizes, "mismatching chunksizes"
        assert tree.chunksizes == original_chunksizes, "original tree was modified"

    def test_persist(self):
        ds1 = xr.Dataset({"a": ("x", np.arange(10))})
        ds2 = xr.Dataset({"b": ("y", np.arange(5))})
        ds3 = xr.Dataset({"c": ("z", np.arange(4))})
        ds4 = xr.Dataset({"d": ("x", np.arange(-5, 5))})

        def fn(x):
            return 2 * x

        expected = xr.DataTree.from_dict(
            {
                "/": fn(ds1).chunk({"x": 5}),
                "/group1": fn(ds2).chunk({"y": 3}),
                "/group2": fn(ds3).chunk({"z": 2}),
                "/group1/subgroup1": fn(ds4).chunk({"x": 5}),
            }
        )
        # Add trivial second layer to the task graph, persist should reduce to one
        tree = xr.DataTree.from_dict(
            {
                "/": fn(ds1.chunk({"x": 5})),
                "/group1": fn(ds2.chunk({"y": 3})),
                "/group2": fn(ds3.chunk({"z": 2})),
                "/group1/subgroup1": fn(ds4.chunk({"x": 5})),
            }
        )
        original_chunksizes = tree.chunksizes
        original_hlg_depths = {
            node.path: len(node.dataset.__dask_graph__().layers)
            for node in tree.subtree
        }

        actual = tree.persist()
        actual_hlg_depths = {
            node.path: len(node.dataset.__dask_graph__().layers)
            for node in actual.subtree
        }

        assert_identical(actual, expected)

        assert actual.chunksizes == original_chunksizes, "chunksizes were modified"
        assert tree.chunksizes == original_chunksizes, (
            "original chunksizes were modified"
        )
        assert all(d == 1 for d in actual_hlg_depths.values()), (
            "unexpected dask graph depth"
        )
        assert all(d == 2 for d in original_hlg_depths.values()), (
            "original dask graph was modified"
        )

    def test_chunk(self):
        ds1 = xr.Dataset({"a": ("x", np.arange(10))})
        ds2 = xr.Dataset({"b": ("y", np.arange(5))})
        ds3 = xr.Dataset({"c": ("z", np.arange(4))})
        ds4 = xr.Dataset({"d": ("x", np.arange(-5, 5))})

        expected = xr.DataTree.from_dict(
            {
                "/": ds1.chunk({"x": 5}),
                "/group1": ds2.chunk({"y": 3}),
                "/group2": ds3.chunk({"z": 2}),
                "/group1/subgroup1": ds4.chunk({"x": 5}),
            }
        )

        tree = xr.DataTree.from_dict(
            {"/": ds1, "/group1": ds2, "/group2": ds3, "/group1/subgroup1": ds4}
        )
        actual = tree.chunk({"x": 5, "y": 3, "z": 2})

        assert_identical(actual, expected)
        assert actual.chunksizes == expected.chunksizes

        with pytest.raises(TypeError, match="invalid type"):
            tree.chunk(None)

        with pytest.raises(TypeError, match="invalid type"):
            tree.chunk((1, 2))

        with pytest.raises(ValueError, match="not found in data dimensions"):
            tree.chunk({"u": 2})
