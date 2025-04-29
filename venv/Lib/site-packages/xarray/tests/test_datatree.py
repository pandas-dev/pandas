import re
import typing
from copy import copy, deepcopy
from textwrap import dedent

import numpy as np
import pytest

import xarray as xr
from xarray import DataArray, Dataset
from xarray.core.coordinates import DataTreeCoordinates
from xarray.core.datatree import DataTree
from xarray.core.datatree_ops import _MAPPED_DOCSTRING_ADDENDUM, insert_doc_addendum
from xarray.core.treenode import NotFoundInTreeError
from xarray.testing import assert_equal, assert_identical
from xarray.tests import assert_array_equal, create_test_data, source_ndarray


class TestTreeCreation:
    def test_empty(self):
        dt = DataTree(name="root")
        assert dt.name == "root"
        assert dt.parent is None
        assert dt.children == {}
        assert_identical(dt.to_dataset(), xr.Dataset())

    def test_unnamed(self):
        dt = DataTree()
        assert dt.name is None

    def test_bad_names(self):
        with pytest.raises(TypeError):
            DataTree(name=5)  # type: ignore[arg-type]

        with pytest.raises(ValueError):
            DataTree(name="folder/data")

    def test_data_arg(self):
        ds = xr.Dataset({"foo": 42})
        tree: DataTree = DataTree(dataset=ds)
        assert_identical(tree.to_dataset(), ds)

        with pytest.raises(TypeError):
            DataTree(dataset=xr.DataArray(42, name="foo"))  # type: ignore[arg-type]


class TestFamilyTree:
    def test_dont_modify_children_inplace(self):
        # GH issue 9196
        child = DataTree()
        DataTree(children={"child": child})
        assert child.parent is None

    def test_create_two_children(self):
        root_data = xr.Dataset({"a": ("y", [6, 7, 8]), "set0": ("x", [9, 10])})
        set1_data = xr.Dataset({"a": 0, "b": 1})
        root = DataTree.from_dict(
            {"/": root_data, "/set1": set1_data, "/set1/set2": None}
        )
        assert root["/set1"].name == "set1"
        assert root["/set1/set2"].name == "set2"

    def test_create_full_tree(self, simple_datatree):
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
    def test_child_gets_named_on_attach(self):
        sue = DataTree()
        mary = DataTree(children={"Sue": sue})  # noqa
        assert mary.children["Sue"].name == "Sue"


class TestPaths:
    def test_path_property(self):
        john = DataTree.from_dict(
            {
                "/Mary/Sue": DataTree(),
            }
        )
        assert john["/Mary/Sue"].path == "/Mary/Sue"
        assert john.path == "/"

    def test_path_roundtrip(self):
        john = DataTree.from_dict(
            {
                "/Mary/Sue": DataTree(),
            }
        )
        assert john["/Mary/Sue"].name == "Sue"

    def test_same_tree(self):
        john = DataTree.from_dict(
            {
                "/Mary": DataTree(),
                "/Kate": DataTree(),
            }
        )
        assert john["/Mary"].same_tree(john["/Kate"])

    def test_relative_paths(self):
        john = DataTree.from_dict(
            {
                "/Mary/Sue": DataTree(),
                "/Annie": DataTree(),
            }
        )
        sue_result = john["Mary/Sue"]
        if isinstance(sue_result, DataTree):
            sue: DataTree = sue_result

        annie_result = john["Annie"]
        if isinstance(annie_result, DataTree):
            annie: DataTree = annie_result

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
    def test_create_with_data(self):
        dat = xr.Dataset({"a": 0})
        john = DataTree(name="john", dataset=dat)

        assert_identical(john.to_dataset(), dat)

        with pytest.raises(TypeError):
            DataTree(name="mary", dataset="junk")  # type: ignore[arg-type]

    def test_set_data(self):
        john = DataTree(name="john")
        dat = xr.Dataset({"a": 0})
        john.dataset = dat  # type: ignore[assignment]

        assert_identical(john.to_dataset(), dat)

        with pytest.raises(TypeError):
            john.dataset = "junk"  # type: ignore[assignment]

    def test_has_data(self):
        john = DataTree(name="john", dataset=xr.Dataset({"a": 0}))
        assert john.has_data

        john_no_data = DataTree(name="john", dataset=None)
        assert not john_no_data.has_data

    def test_is_hollow(self):
        john = DataTree(dataset=xr.Dataset({"a": 0}))
        assert john.is_hollow

        eve = DataTree(children={"john": john})
        assert eve.is_hollow

        eve.dataset = xr.Dataset({"a": 1})  # type: ignore[assignment]
        assert not eve.is_hollow


class TestToDataset:
    def test_to_dataset(self):
        base = xr.Dataset(coords={"a": 1})
        sub = xr.Dataset(coords={"b": 2})
        tree = DataTree.from_dict({"/": base, "/sub": sub})
        subtree = typing.cast(DataTree, tree["sub"])

        assert_identical(tree.to_dataset(inherited=False), base)
        assert_identical(subtree.to_dataset(inherited=False), sub)

        sub_and_base = xr.Dataset(coords={"a": 1, "b": 2})
        assert_identical(tree.to_dataset(inherited=True), base)
        assert_identical(subtree.to_dataset(inherited=True), sub_and_base)


class TestVariablesChildrenNameCollisions:
    def test_parent_already_has_variable_with_childs_name(self):
        with pytest.raises(KeyError, match="already contains a variable named a"):
            DataTree.from_dict({"/": xr.Dataset({"a": [0], "b": 1}), "/a": None})

    def test_parent_already_has_variable_with_childs_name_update(self):
        dt = DataTree(dataset=xr.Dataset({"a": [0], "b": 1}))
        with pytest.raises(ValueError, match="already contains a variable named a"):
            dt.update({"a": DataTree()})

    def test_assign_when_already_child_with_variables_name(self):
        dt = DataTree.from_dict(
            {
                "/a": DataTree(),
            }
        )

        with pytest.raises(ValueError, match="node already contains a variable"):
            dt.dataset = xr.Dataset({"a": 0})  # type: ignore[assignment]

        dt.dataset = xr.Dataset()  # type: ignore[assignment]

        new_ds = dt.to_dataset().assign(a=xr.DataArray(0))
        with pytest.raises(ValueError, match="node already contains a variable"):
            dt.dataset = new_ds  # type: ignore[assignment]


class TestGet: ...


class TestGetItem:
    def test_getitem_node(self):
        folder1 = DataTree.from_dict(
            {
                "/results/highres": DataTree(),
            }
        )

        assert folder1["results"].name == "results"
        assert folder1["results/highres"].name == "highres"

    def test_getitem_self(self):
        dt = DataTree()
        assert dt["."] is dt

    def test_getitem_single_data_variable(self):
        data = xr.Dataset({"temp": [0, 50]})
        results = DataTree(name="results", dataset=data)
        assert_identical(results["temp"], data["temp"])

    def test_getitem_single_data_variable_from_node(self):
        data = xr.Dataset({"temp": [0, 50]})
        folder1 = DataTree.from_dict(
            {
                "/results/highres": data,
            }
        )
        assert_identical(folder1["results/highres/temp"], data["temp"])

    def test_getitem_nonexistent_node(self):
        folder1 = DataTree.from_dict({"/results": DataTree()}, name="folder1")
        with pytest.raises(KeyError):
            folder1["results/highres"]

    def test_getitem_nonexistent_variable(self):
        data = xr.Dataset({"temp": [0, 50]})
        results = DataTree(name="results", dataset=data)
        with pytest.raises(KeyError):
            results["pressure"]

    @pytest.mark.xfail(reason="Should be deprecated in favour of .subset")
    def test_getitem_multiple_data_variables(self):
        data = xr.Dataset({"temp": [0, 50], "p": [5, 8, 7]})
        results = DataTree(name="results", dataset=data)
        assert_identical(results[["temp", "p"]], data[["temp", "p"]])  # type: ignore[index]

    @pytest.mark.xfail(
        reason="Indexing needs to return whole tree (GH https://github.com/xarray-contrib/datatree/issues/77)"
    )
    def test_getitem_dict_like_selection_access_to_dataset(self):
        data = xr.Dataset({"temp": [0, 50]})
        results = DataTree(name="results", dataset=data)
        assert_identical(results[{"temp": 1}], data[{"temp": 1}])  # type: ignore[index]


class TestUpdate:
    def test_update(self):
        dt = DataTree()
        dt.update({"foo": xr.DataArray(0), "a": DataTree()})
        expected = DataTree.from_dict({"/": xr.Dataset({"foo": 0}), "a": None})
        assert_equal(dt, expected)
        assert dt.groups == ("/", "/a")

    def test_update_new_named_dataarray(self):
        da = xr.DataArray(name="temp", data=[0, 50])
        folder1 = DataTree(name="folder1")
        folder1.update({"results": da})
        expected = da.rename("results")
        assert_equal(folder1["results"], expected)

    def test_update_doesnt_alter_child_name(self):
        dt = DataTree()
        dt.update({"foo": xr.DataArray(0), "a": DataTree(name="b")})
        assert "a" in dt.children
        child = dt["a"]
        assert child.name == "a"

    def test_update_overwrite(self):
        actual = DataTree.from_dict({"a": DataTree(xr.Dataset({"x": 1}))})
        actual.update({"a": DataTree(xr.Dataset({"x": 2}))})
        expected = DataTree.from_dict({"a": DataTree(xr.Dataset({"x": 2}))})
        assert_equal(actual, expected)

    def test_update_coordinates(self):
        expected = DataTree.from_dict({"/": xr.Dataset(coords={"a": 1})})
        actual = DataTree.from_dict({"/": xr.Dataset()})
        actual.update(xr.Dataset(coords={"a": 1}))
        assert_equal(actual, expected)

    def test_update_inherited_coords(self):
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
        actual_node = actual.children["b"].to_dataset(inherited=False)
        expected_node = expected.children["b"].to_dataset(inherited=False)
        assert_identical(actual_node, expected_node)


class TestCopy:
    def test_copy(self, create_test_datatree):
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

    def test_copy_subtree(self):
        dt = DataTree.from_dict({"/level1/level2/level3": xr.Dataset()})

        actual = dt["/level1/level2"].copy()
        expected = DataTree.from_dict({"/level3": xr.Dataset()}, name="level2")

        assert_identical(actual, expected)

    def test_copy_coord_inheritance(self) -> None:
        tree = DataTree.from_dict(
            {"/": xr.Dataset(coords={"x": [0, 1]}), "/c": DataTree()}
        )
        tree2 = tree.copy()
        node_ds = tree2.children["c"].to_dataset(inherited=False)
        assert_identical(node_ds, xr.Dataset())

    def test_deepcopy(self, create_test_datatree):
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
    def test_copy_with_data(self, create_test_datatree):
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
    def test_setitem_new_child_node(self):
        john = DataTree(name="john")
        mary = DataTree(name="mary")
        john["mary"] = mary

        grafted_mary = john["mary"]
        assert grafted_mary.parent is john
        assert grafted_mary.name == "mary"

    def test_setitem_unnamed_child_node_becomes_named(self):
        john2 = DataTree(name="john2")
        john2["sonny"] = DataTree()
        assert john2["sonny"].name == "sonny"

    def test_setitem_new_grandchild_node(self):
        john = DataTree.from_dict({"/Mary/Rose": DataTree()})
        new_rose = DataTree(dataset=xr.Dataset({"x": 0}))
        john["Mary/Rose"] = new_rose

        grafted_rose = john["Mary/Rose"]
        assert grafted_rose.parent is john["/Mary"]
        assert grafted_rose.name == "Rose"

    def test_grafted_subtree_retains_name(self):
        subtree = DataTree(name="original_subtree_name")
        root = DataTree(name="root")
        root["new_subtree_name"] = subtree  # noqa
        assert subtree.name == "original_subtree_name"

    def test_setitem_new_empty_node(self):
        john = DataTree(name="john")
        john["mary"] = DataTree()
        mary = john["mary"]
        assert isinstance(mary, DataTree)
        assert_identical(mary.to_dataset(), xr.Dataset())

    def test_setitem_overwrite_data_in_node_with_none(self):
        john = DataTree.from_dict({"/mary": xr.Dataset()}, name="john")

        john["mary"] = DataTree()
        assert_identical(john["mary"].to_dataset(), xr.Dataset())

        john.dataset = xr.Dataset()  # type: ignore[assignment]
        with pytest.raises(ValueError, match="has no name"):
            john["."] = DataTree()

    @pytest.mark.xfail(reason="assigning Datasets doesn't yet create new nodes")
    def test_setitem_dataset_on_this_node(self):
        data = xr.Dataset({"temp": [0, 50]})
        results = DataTree(name="results")
        results["."] = data
        assert_identical(results.to_dataset(), data)

    @pytest.mark.xfail(reason="assigning Datasets doesn't yet create new nodes")
    def test_setitem_dataset_as_new_node(self):
        data = xr.Dataset({"temp": [0, 50]})
        folder1 = DataTree(name="folder1")
        folder1["results"] = data
        assert_identical(folder1["results"].to_dataset(), data)

    @pytest.mark.xfail(reason="assigning Datasets doesn't yet create new nodes")
    def test_setitem_dataset_as_new_node_requiring_intermediate_nodes(self):
        data = xr.Dataset({"temp": [0, 50]})
        folder1 = DataTree(name="folder1")
        folder1["results/highres"] = data
        assert_identical(folder1["results/highres"].to_dataset(), data)

    def test_setitem_named_dataarray(self):
        da = xr.DataArray(name="temp", data=[0, 50])
        folder1 = DataTree(name="folder1")
        folder1["results"] = da
        expected = da.rename("results")
        assert_equal(folder1["results"], expected)

    def test_setitem_unnamed_dataarray(self):
        data = xr.DataArray([0, 50])
        folder1 = DataTree(name="folder1")
        folder1["results"] = data
        assert_equal(folder1["results"], data)

    def test_setitem_variable(self):
        var = xr.Variable(data=[0, 50], dims="x")
        folder1 = DataTree(name="folder1")
        folder1["results"] = var
        assert_equal(folder1["results"], xr.DataArray(var))

    def test_setitem_coerce_to_dataarray(self):
        folder1 = DataTree(name="folder1")
        folder1["results"] = 0
        assert_equal(folder1["results"], xr.DataArray(0))

    def test_setitem_add_new_variable_to_empty_node(self):
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

    def test_setitem_dataarray_replace_existing_node(self):
        t = xr.Dataset({"temp": [0, 50]})
        results = DataTree(name="results", dataset=t)
        p = xr.DataArray(data=[2, 3])
        results["pressure"] = p
        expected = t.assign(pressure=p)
        assert_identical(results.to_dataset(), expected)


class TestCoords:
    def test_properties(self):
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

    def test_modify(self):
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

        actual = dt.copy()
        del actual.coords["b"]
        expected = dt.reset_coords("b", drop=True)
        assert_identical(expected, actual)

        with pytest.raises(KeyError):
            del dt.coords["not_found"]

        with pytest.raises(KeyError):
            del dt.coords["foo"]

        actual = dt.copy(deep=True)
        actual.coords.update({"c": 11})
        expected = dt.assign_coords({"c": 11})
        assert_identical(expected, actual)

        # regression test for GH3746
        del actual.coords["x"]
        assert "x" not in actual.xindexes

        # test that constructors can also handle the `DataTreeCoordinates` object
        ds2 = Dataset(coords=dt.coords)
        assert_identical(ds2.coords, dt.coords)
        da = DataArray(coords=dt.coords)
        assert_identical(da.coords, dt.coords)

        # DataTree constructor doesn't accept coords= but should still be able to handle DatasetCoordinates
        dt2 = DataTree(dataset=dt.coords)
        assert_identical(dt2.coords, dt.coords)

    def test_inherited(self):
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

        assert set(child.coords) == {"x", "y", "a", "b"}

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
            del child["b"]

        # TODO requires a fix for #9472
        # actual = child.copy(deep=True)
        # actual.coords.update({"c": 11})
        # expected = child.assign_coords({"c": 11})
        # assert_identical(expected, actual)


def test_delitem():
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
    def test_data_in_root(self):
        dat = xr.Dataset()
        dt = DataTree.from_dict({"/": dat})
        assert dt.name is None
        assert dt.parent is None
        assert dt.children == {}
        assert_identical(dt.to_dataset(), dat)

    def test_one_layer(self):
        dat1, dat2 = xr.Dataset({"a": 1}), xr.Dataset({"b": 2})
        dt = DataTree.from_dict({"run1": dat1, "run2": dat2})
        assert_identical(dt.to_dataset(), xr.Dataset())
        assert dt.name is None
        assert_identical(dt["run1"].to_dataset(), dat1)
        assert dt["run1"].children == {}
        assert_identical(dt["run2"].to_dataset(), dat2)
        assert dt["run2"].children == {}

    def test_two_layers(self):
        dat1, dat2 = xr.Dataset({"a": 1}), xr.Dataset({"a": [1, 2]})
        dt = DataTree.from_dict({"highres/run": dat1, "lowres/run": dat2})
        assert "highres" in dt.children
        assert "lowres" in dt.children
        highres_run = dt["highres/run"]
        assert_identical(highres_run.to_dataset(), dat1)

    def test_nones(self):
        dt = DataTree.from_dict({"d": None, "d/e": None})
        assert [node.name for node in dt.subtree] == [None, "d", "e"]
        assert [node.path for node in dt.subtree] == ["/", "/d", "/d/e"]
        assert_identical(dt["d/e"].to_dataset(), xr.Dataset())

    def test_full(self, simple_datatree):
        dt = simple_datatree
        paths = list(node.path for node in dt.subtree)
        assert paths == [
            "/",
            "/set1",
            "/set2",
            "/set3",
            "/set1/set1",
            "/set1/set2",
            "/set2/set1",
        ]

    def test_datatree_values(self):
        dat1 = DataTree(dataset=xr.Dataset({"a": 1}))
        expected = DataTree()
        expected["a"] = dat1

        actual = DataTree.from_dict({"a": dat1})

        assert_identical(actual, expected)

    def test_roundtrip(self, simple_datatree):
        dt = simple_datatree
        roundtrip = DataTree.from_dict(dt.to_dict())
        assert roundtrip.equals(dt)

    @pytest.mark.xfail
    def test_roundtrip_unnamed_root(self, simple_datatree):
        # See GH81

        dt = simple_datatree
        dt.name = "root"
        roundtrip = DataTree.from_dict(dt.to_dict())
        assert roundtrip.equals(dt)

    def test_insertion_order(self):
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

    def test_array_values(self):
        data = {"foo": xr.DataArray(1, name="bar")}
        with pytest.raises(TypeError):
            DataTree.from_dict(data)  # type: ignore[arg-type]


class TestDatasetView:
    def test_view_contents(self):
        ds = create_test_data()
        dt = DataTree(dataset=ds)
        assert ds.identical(
            dt.dataset
        )  # this only works because Dataset.identical doesn't check types
        assert isinstance(dt.dataset, xr.Dataset)

    def test_immutability(self):
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

    def test_methods(self):
        ds = create_test_data()
        dt = DataTree(dataset=ds)
        assert ds.mean().identical(dt.dataset.mean())
        assert isinstance(dt.dataset.mean(), xr.Dataset)

    def test_arithmetic(self, create_test_datatree):
        dt = create_test_datatree()
        expected = create_test_datatree(modify=lambda ds: 10.0 * ds)[
            "set1"
        ].to_dataset()
        result = 10.0 * dt["set1"].dataset
        assert result.identical(expected)

    def test_init_via_type(self):
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


class TestAccess:
    def test_attribute_access(self, create_test_datatree):
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

    def test_ipython_key_completions(self, create_test_datatree):
        dt = create_test_datatree()
        key_completions = dt._ipython_key_completions_()

        node_keys = [node.path[1:] for node in dt.subtree]
        assert all(node_key in key_completions for node_key in node_keys)

        var_keys = list(dt.variables.keys())
        assert all(var_key in key_completions for var_key in var_keys)

    def test_operation_with_attrs_but_no_data(self):
        # tests bug from xarray-datatree GH262
        xs = xr.Dataset({"testvar": xr.DataArray(np.ones((2, 3)))})
        dt = DataTree.from_dict({"node1": xs, "node2": xs})
        dt.attrs["test_key"] = 1  # sel works fine without this line
        dt.sel(dim_0=0)


class TestRepr:

    def test_repr_four_nodes(self):
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

    def test_repr_two_children(self):
        tree = DataTree.from_dict(
            {
                "/": Dataset(coords={"x": [1.0]}),
                "/first_child": None,
                "/second_child": Dataset({"foo": ("x", [0.0])}),
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
                Inherited coordinates:
                  * x        (x) float64 8B 1.0
                Data variables:
                    foo      (x) float64 8B 0.0
            """
        ).strip()
        assert result == expected

    def test_repr_inherited_dims(self):
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


def _exact_match(message: str) -> str:
    return re.escape(dedent(message).strip())
    return "^" + re.escape(dedent(message.rstrip())) + "$"


class TestInheritance:
    def test_inherited_dims(self):
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
        assert dt.b.to_dataset(inherited=True).sizes == {"y": 1}
        assert dt.b.to_dataset(inherited=False).sizes == {"y": 1}

    def test_inherited_coords_index(self):
        dt = DataTree.from_dict(
            {
                "/": xr.Dataset({"d": (("x",), [1, 2])}, coords={"x": [2, 3]}),
                "/b": xr.Dataset({"e": (("y",), [3])}),
            }
        )
        assert "x" in dt["/b"].indexes
        assert "x" in dt["/b"].coords
        xr.testing.assert_identical(dt["/x"], dt["/b/x"])

    def test_inherited_coords_override(self):
        dt = DataTree.from_dict(
            {
                "/": xr.Dataset(coords={"x": 1, "y": 2}),
                "/b": xr.Dataset(coords={"x": 4, "z": 3}),
            }
        )
        assert dt.coords.keys() == {"x", "y"}
        root_coords = {"x": 1, "y": 2}
        sub_coords = {"x": 4, "y": 2, "z": 3}
        xr.testing.assert_equal(dt["/x"], xr.DataArray(1, coords=root_coords))
        xr.testing.assert_equal(dt["/y"], xr.DataArray(2, coords=root_coords))
        assert dt["/b"].coords.keys() == {"x", "y", "z"}
        xr.testing.assert_equal(dt["/b/x"], xr.DataArray(4, coords=sub_coords))
        xr.testing.assert_equal(dt["/b/y"], xr.DataArray(2, coords=sub_coords))
        xr.testing.assert_equal(dt["/b/z"], xr.DataArray(3, coords=sub_coords))

    def test_inconsistent_dims(self):
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

    def test_inconsistent_child_indexes(self):
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
        dt.dataset = xr.Dataset(coords={"x": [1.0]})  # type: ignore[assignment]
        dt["/b"] = DataTree()
        with pytest.raises(ValueError, match=expected_msg):
            dt["/b"].dataset = xr.Dataset(coords={"x": [2.0]})

        b = DataTree(xr.Dataset(coords={"x": [2.0]}))
        with pytest.raises(ValueError, match=expected_msg):
            DataTree(dataset=xr.Dataset(coords={"x": [1.0]}), children={"b": b})

    def test_inconsistent_grandchild_indexes(self):
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
        dt.dataset = xr.Dataset(coords={"x": [1.0]})  # type: ignore[assignment]
        dt["/b/c"] = DataTree()
        with pytest.raises(ValueError, match=expected_msg):
            dt["/b/c"].dataset = xr.Dataset(coords={"x": [2.0]})

        c = DataTree(xr.Dataset(coords={"x": [2.0]}))
        b = DataTree(children={"c": c})
        with pytest.raises(ValueError, match=expected_msg):
            DataTree(dataset=xr.Dataset(coords={"x": [1.0]}), children={"b": b})

    def test_inconsistent_grandchild_dims(self):
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
    def test_drop_nodes(self):
        sue = DataTree.from_dict({"Mary": None, "Kate": None, "Ashley": None})

        # test drop just one node
        dropped_one = sue.drop_nodes(names="Mary")
        assert "Mary" not in dropped_one.children

        # test drop multiple nodes
        dropped = sue.drop_nodes(names=["Mary", "Kate"])
        assert not set(["Mary", "Kate"]).intersection(set(dropped.children))
        assert "Ashley" in dropped.children

        # test raise
        with pytest.raises(KeyError, match="nodes {'Mary'} not present"):
            dropped.drop_nodes(names=["Mary", "Ashley"])

        # test ignore
        childless = dropped.drop_nodes(names=["Mary", "Ashley"], errors="ignore")
        assert childless.children == {}

    def test_assign(self):
        dt = DataTree()
        expected = DataTree.from_dict({"/": xr.Dataset({"foo": 0}), "/a": None})

        # kwargs form
        result = dt.assign(foo=xr.DataArray(0), a=DataTree())
        assert_equal(result, expected)

        # dict form
        result = dt.assign({"foo": xr.DataArray(0), "a": DataTree()})
        assert_equal(result, expected)


class TestPipe:
    def test_noop(self, create_test_datatree):
        dt = create_test_datatree()

        actual = dt.pipe(lambda tree: tree)
        assert actual.identical(dt)

    def test_params(self, create_test_datatree):
        dt = create_test_datatree()

        def f(tree, **attrs):
            return tree.assign(arr_with_attrs=xr.Variable("dim0", [], attrs=attrs))

        attrs = {"x": 1, "y": 2, "z": 3}

        actual = dt.pipe(f, **attrs)
        assert actual["arr_with_attrs"].attrs == attrs

    def test_named_self(self, create_test_datatree):
        dt = create_test_datatree()

        def f(x, tree, y):
            tree.attrs.update({"x": x, "y": y})
            return tree

        attrs = {"x": 1, "y": 2}

        actual = dt.pipe((f, "tree"), **attrs)

        assert actual is dt and actual.attrs == attrs


class TestSubset:
    def test_match(self):
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

    def test_filter(self):
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


class TestDSMethodInheritance:
    def test_dataset_method(self):
        ds = xr.Dataset({"a": ("x", [1, 2, 3])})
        dt = DataTree.from_dict(
            {
                "/": ds,
                "/results": ds,
            }
        )

        expected = DataTree.from_dict(
            {
                "/": ds.isel(x=1),
                "/results": ds.isel(x=1),
            }
        )

        result = dt.isel(x=1)
        assert_equal(result, expected)

    def test_reduce_method(self):
        ds = xr.Dataset({"a": ("x", [False, True, False])})
        dt = DataTree.from_dict({"/": ds, "/results": ds})

        expected = DataTree.from_dict({"/": ds.any(), "/results": ds.any()})

        result = dt.any()
        assert_equal(result, expected)

    def test_nan_reduce_method(self):
        ds = xr.Dataset({"a": ("x", [1, 2, 3])})
        dt = DataTree.from_dict({"/": ds, "/results": ds})

        expected = DataTree.from_dict({"/": ds.mean(), "/results": ds.mean()})

        result = dt.mean()
        assert_equal(result, expected)

    def test_cum_method(self):
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


class TestOps:
    def test_binary_op_on_int(self):
        ds1 = xr.Dataset({"a": [5], "b": [3]})
        ds2 = xr.Dataset({"x": [0.1, 0.2], "y": [10, 20]})
        dt = DataTree.from_dict({"/": ds1, "/subnode": ds2})

        expected = DataTree.from_dict({"/": ds1 * 5, "/subnode": ds2 * 5})

        # TODO: Remove ignore when ops.py is migrated?
        result: DataTree = dt * 5  # type: ignore[assignment,operator]
        assert_equal(result, expected)

    def test_binary_op_on_dataset(self):
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

    def test_binary_op_on_datatree(self):
        ds1 = xr.Dataset({"a": [5], "b": [3]})
        ds2 = xr.Dataset({"x": [0.1, 0.2], "y": [10, 20]})

        dt = DataTree.from_dict({"/": ds1, "/subnode": ds2})

        expected = DataTree.from_dict({"/": ds1 * ds1, "/subnode": ds2 * ds2})

        # TODO: Remove ignore when ops.py is migrated?
        result = dt * dt  # type: ignore[operator]
        assert_equal(result, expected)


class TestUFuncs:
    def test_tree(self, create_test_datatree):
        dt = create_test_datatree()
        expected = create_test_datatree(modify=lambda ds: np.sin(ds))
        result_tree = np.sin(dt)
        assert_equal(result_tree, expected)


class TestDocInsertion:
    """Tests map_over_subtree docstring injection."""

    def test_standard_doc(self):

        dataset_doc = dedent(
            """\
            Manually trigger loading and/or computation of this dataset's data
                    from disk or a remote source into memory and return this dataset.
                    Unlike compute, the original dataset is modified and returned.

                    Normally, it should not be necessary to call this method in user code,
                    because all xarray functions should either work on deferred data or
                    load data automatically. However, this method can be necessary when
                    working with many file objects on disk.

                    Parameters
                    ----------
                    **kwargs : dict
                        Additional keyword arguments passed on to ``dask.compute``.

                    See Also
                    --------
                    dask.compute"""
        )

        expected_doc = dedent(
            """\
            Manually trigger loading and/or computation of this dataset's data
                    from disk or a remote source into memory and return this dataset.
                    Unlike compute, the original dataset is modified and returned.

                    .. note::
                        This method was copied from xarray.Dataset, but has been altered to
                        call the method on the Datasets stored in every node of the
                        subtree. See the `map_over_subtree` function for more details.

                    Normally, it should not be necessary to call this method in user code,
                    because all xarray functions should either work on deferred data or
                    load data automatically. However, this method can be necessary when
                    working with many file objects on disk.

                    Parameters
                    ----------
                    **kwargs : dict
                        Additional keyword arguments passed on to ``dask.compute``.

                    See Also
                    --------
                    dask.compute"""
        )

        wrapped_doc = insert_doc_addendum(dataset_doc, _MAPPED_DOCSTRING_ADDENDUM)

        assert expected_doc == wrapped_doc

    def test_one_liner(self):
        mixin_doc = "Same as abs(a)."

        expected_doc = dedent(
            """\
            Same as abs(a).

            This method was copied from xarray.Dataset, but has been altered to call the
                method on the Datasets stored in every node of the subtree. See the
                `map_over_subtree` function for more details."""
        )

        actual_doc = insert_doc_addendum(mixin_doc, _MAPPED_DOCSTRING_ADDENDUM)
        assert expected_doc == actual_doc

    def test_none(self):
        actual_doc = insert_doc_addendum(None, _MAPPED_DOCSTRING_ADDENDUM)
        assert actual_doc is None
