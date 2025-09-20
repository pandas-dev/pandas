from __future__ import annotations

import sys
from textwrap import dedent

import numpy as np
import pandas as pd
import pytest

import xarray as xr
from xarray.core import formatting
from xarray.core.indexes import Index
from xarray.tests import requires_cftime, requires_dask, requires_netCDF4


class CustomIndex(Index):
    names: tuple[str, ...]

    def __init__(self, names: tuple[str, ...]):
        self.names = names

    def __repr__(self):
        return f"CustomIndex(coords={self.names})"


class TestFormatting:
    def test_get_indexer_at_least_n_items(self) -> None:
        cases = [
            ((20,), (slice(10),), (slice(-10, None),)),
            ((3, 20), (0, slice(10)), (-1, slice(-10, None))),
            ((2, 10), (0, slice(10)), (-1, slice(-10, None))),
            ((2, 5), (slice(2), slice(None)), (slice(-2, None), slice(None))),
            ((1, 2, 5), (0, slice(2), slice(None)), (-1, slice(-2, None), slice(None))),
            ((2, 3, 5), (0, slice(2), slice(None)), (-1, slice(-2, None), slice(None))),
            (
                (1, 10, 1),
                (0, slice(10), slice(None)),
                (-1, slice(-10, None), slice(None)),
            ),
            (
                (2, 5, 1),
                (slice(2), slice(None), slice(None)),
                (slice(-2, None), slice(None), slice(None)),
            ),
            ((2, 5, 3), (0, slice(4), slice(None)), (-1, slice(-4, None), slice(None))),
            (
                (2, 3, 3),
                (slice(2), slice(None), slice(None)),
                (slice(-2, None), slice(None), slice(None)),
            ),
        ]
        for shape, start_expected, end_expected in cases:
            actual = formatting._get_indexer_at_least_n_items(shape, 10, from_end=False)
            assert start_expected == actual
            actual = formatting._get_indexer_at_least_n_items(shape, 10, from_end=True)
            assert end_expected == actual

    def test_first_n_items(self) -> None:
        array = np.arange(100).reshape(10, 5, 2)
        for n in [3, 10, 13, 100, 200]:
            actual = formatting.first_n_items(array, n)
            expected = array.flat[:n]
            assert (expected == actual).all()

        with pytest.raises(ValueError, match=r"at least one item"):
            formatting.first_n_items(array, 0)

    def test_last_n_items(self) -> None:
        array = np.arange(100).reshape(10, 5, 2)
        for n in [3, 10, 13, 100, 200]:
            actual = formatting.last_n_items(array, n)
            expected = array.flat[-n:]
            assert (expected == actual).all()

        with pytest.raises(ValueError, match=r"at least one item"):
            formatting.first_n_items(array, 0)

    def test_last_item(self) -> None:
        array = np.arange(100)

        reshape = ((10, 10), (1, 100), (2, 2, 5, 5))
        expected = np.array([99])

        for r in reshape:
            result = formatting.last_item(array.reshape(r))
            assert result == expected

    def test_format_item(self) -> None:
        cases = [
            (pd.Timestamp("2000-01-01T12"), "2000-01-01T12:00:00"),
            (pd.Timestamp("2000-01-01"), "2000-01-01"),
            (pd.Timestamp("NaT"), "NaT"),
            (pd.Timedelta("10 days 1 hour"), "10 days 01:00:00"),
            (pd.Timedelta("-3 days"), "-3 days +00:00:00"),
            (pd.Timedelta("3 hours"), "0 days 03:00:00"),
            (pd.Timedelta("NaT"), "NaT"),
            ("foo", "'foo'"),
            (b"foo", "b'foo'"),
            (1, "1"),
            (1.0, "1.0"),
            (np.float16(1.1234), "1.123"),
            (np.float32(1.0111111), "1.011"),
            (np.float64(22.222222), "22.22"),
            (np.zeros((1, 1)), "[[0.]]"),
            (np.zeros(2), "[0. 0.]"),
            (np.zeros((2, 2)), "[[0. 0.]\n [0. 0.]]"),
        ]
        for item, expected in cases:
            actual = formatting.format_item(item)
            assert expected == actual

    def test_format_items(self) -> None:
        cases = [
            (np.arange(4) * np.timedelta64(1, "D"), "0 days 1 days 2 days 3 days"),
            (
                np.arange(4) * np.timedelta64(3, "h"),
                "00:00:00 03:00:00 06:00:00 09:00:00",
            ),
            (
                np.arange(4) * np.timedelta64(500, "ms"),
                "00:00:00 00:00:00.500000 00:00:01 00:00:01.500000",
            ),
            (pd.to_timedelta(["NaT", "0s", "1s", "NaT"]), "NaT 00:00:00 00:00:01 NaT"),  # type: ignore[arg-type, unused-ignore]
            (
                pd.to_timedelta(["1 day 1 hour", "1 day", "0 hours"]),  # type: ignore[arg-type, unused-ignore]
                "1 days 01:00:00 1 days 00:00:00 0 days 00:00:00",
            ),
            ([1, 2, 3], "1 2 3"),
        ]
        for item, expected in cases:
            actual = " ".join(formatting.format_items(item))
            assert expected == actual

    def test_format_array_flat(self) -> None:
        actual = formatting.format_array_flat(np.arange(100), 2)
        expected = "..."
        assert expected == actual

        actual = formatting.format_array_flat(np.arange(100), 9)
        expected = "0 ... 99"
        assert expected == actual

        actual = formatting.format_array_flat(np.arange(100), 10)
        expected = "0 1 ... 99"
        assert expected == actual

        actual = formatting.format_array_flat(np.arange(100), 13)
        expected = "0 1 ... 98 99"
        assert expected == actual

        actual = formatting.format_array_flat(np.arange(100), 15)
        expected = "0 1 2 ... 98 99"
        assert expected == actual

        # NB: Probably not ideal; an alternative would be cutting after the
        # first ellipsis
        actual = formatting.format_array_flat(np.arange(100.0), 11)
        expected = "0.0 ... ..."
        assert expected == actual

        actual = formatting.format_array_flat(np.arange(100.0), 12)
        expected = "0.0 ... 99.0"
        assert expected == actual

        actual = formatting.format_array_flat(np.arange(3), 5)
        expected = "0 1 2"
        assert expected == actual

        actual = formatting.format_array_flat(np.arange(4.0), 11)
        expected = "0.0 ... 3.0"
        assert expected == actual

        actual = formatting.format_array_flat(np.arange(0), 0)
        expected = ""
        assert expected == actual

        actual = formatting.format_array_flat(np.arange(1), 1)
        expected = "0"
        assert expected == actual

        actual = formatting.format_array_flat(np.arange(2), 3)
        expected = "0 1"
        assert expected == actual

        actual = formatting.format_array_flat(np.arange(4), 7)
        expected = "0 1 2 3"
        assert expected == actual

        actual = formatting.format_array_flat(np.arange(5), 7)
        expected = "0 ... 4"
        assert expected == actual

        long_str = [" ".join(["hello world" for _ in range(100)])]
        actual = formatting.format_array_flat(np.asarray([long_str]), 21)
        expected = "'hello world hello..."
        assert expected == actual

    def test_pretty_print(self) -> None:
        assert formatting.pretty_print("abcdefghij", 8) == "abcde..."
        assert formatting.pretty_print("ß", 1) == "ß"

    def test_maybe_truncate(self) -> None:
        assert formatting.maybe_truncate("ß", 10) == "ß"

    def test_format_timestamp_invalid_pandas_format(self) -> None:
        expected = "2021-12-06 17:00:00 00"
        with pytest.raises(ValueError):
            formatting.format_timestamp(expected)

    def test_format_timestamp_out_of_bounds(self) -> None:
        from datetime import datetime

        date = datetime(1300, 12, 1)
        expected = "1300-12-01"
        result = formatting.format_timestamp(date)
        assert result == expected

        date = datetime(2300, 12, 1)
        expected = "2300-12-01"
        result = formatting.format_timestamp(date)
        assert result == expected

    def test_attribute_repr(self) -> None:
        short = formatting.summarize_attr("key", "Short string")
        long = formatting.summarize_attr("key", 100 * "Very long string ")
        newlines = formatting.summarize_attr("key", "\n\n\n")
        tabs = formatting.summarize_attr("key", "\t\t\t")
        assert short == "    key: Short string"
        assert len(long) <= 80
        assert long.endswith("...")
        assert "\n" not in newlines
        assert "\t" not in tabs

    def test_index_repr(self) -> None:
        coord_names = ("x", "y")
        index = CustomIndex(coord_names)
        names = ("x",)

        normal = formatting.summarize_index(names, index, col_width=20)
        assert names[0] in normal
        assert len(normal.splitlines()) == len(names)
        assert "CustomIndex" in normal

        class IndexWithInlineRepr(CustomIndex):
            def _repr_inline_(self, max_width: int):
                return f"CustomIndex[{', '.join(self.names)}]"

        index = IndexWithInlineRepr(coord_names)
        inline = formatting.summarize_index(names, index, col_width=20)
        assert names[0] in inline
        assert index._repr_inline_(max_width=40) in inline

    @pytest.mark.parametrize(
        "names",
        (
            ("x",),
            ("x", "y"),
            ("x", "y", "z"),
            ("x", "y", "z", "a"),
        ),
    )
    def test_index_repr_grouping(self, names) -> None:
        index = CustomIndex(names)

        normal = formatting.summarize_index(names, index, col_width=20)
        assert all(name in normal for name in names)
        assert len(normal.splitlines()) == len(names)
        assert "CustomIndex" in normal

        hint_chars = [line[2] for line in normal.splitlines()]

        if len(names) <= 1:
            assert hint_chars == [" "]
        else:
            assert hint_chars[0] == "┌" and hint_chars[-1] == "└"
            assert len(names) == 2 or hint_chars[1:-1] == ["│"] * (len(names) - 2)

    def test_diff_array_repr(self) -> None:
        da_a = xr.DataArray(
            np.array([[1, 2, 3], [4, 5, 6]], dtype="int64"),
            dims=("x", "y"),
            coords={
                "x": np.array(["a", "b"], dtype="U1"),
                "y": np.array([1, 2, 3], dtype="int64"),
            },
            attrs={"units": "m", "description": "desc"},
        )

        da_b = xr.DataArray(
            np.array([1, 2], dtype="int64"),
            dims="x",
            coords={
                "x": np.array(["a", "c"], dtype="U1"),
                "label": ("x", np.array([1, 2], dtype="int64")),
            },
            attrs={"units": "kg"},
        )

        byteorder = "<" if sys.byteorder == "little" else ">"
        expected = dedent(
            f"""\
        Left and right DataArray objects are not identical
        Differing dimensions:
            (x: 2, y: 3) != (x: 2)
        Differing values:
        L
            array([[1, 2, 3],
                   [4, 5, 6]], dtype=int64)
        R
            array([1, 2], dtype=int64)
        Differing coordinates:
        L * x        (x) {byteorder}U1 8B 'a' 'b'
        R * x        (x) {byteorder}U1 8B 'a' 'c'
        Coordinates only on the left object:
          * y        (y) int64 24B 1 2 3
        Coordinates only on the right object:
            label    (x) int64 16B 1 2
        Differing attributes:
        L   units: m
        R   units: kg
        Attributes only on the left object:
            description: desc"""
        )

        actual = formatting.diff_array_repr(da_a, da_b, "identical")
        try:
            assert actual == expected
        except AssertionError:
            # depending on platform, dtype may not be shown in numpy array repr
            assert actual == expected.replace(", dtype=int64", "")

        da_a = xr.DataArray(
            np.array([[1, 2, 3], [4, 5, 6]], dtype="int8"),
            dims=("x", "y"),
            coords=xr.Coordinates(
                {
                    "x": np.array([True, False], dtype="bool"),
                    "y": np.array([1, 2, 3], dtype="int16"),
                },
                indexes={"y": CustomIndex(("y",))},
            ),
        )

        da_b = xr.DataArray(
            np.array([1, 2], dtype="int8"),
            dims="x",
            coords=xr.Coordinates(
                {
                    "x": np.array([True, False], dtype="bool"),
                    "label": ("x", np.array([1, 2], dtype="int16")),
                },
                indexes={"label": CustomIndex(("label",))},
            ),
        )

        expected = dedent(
            """\
            Left and right DataArray objects are not equal
            Differing dimensions:
                (x: 2, y: 3) != (x: 2)
            Differing values:
            L
                array([[1, 2, 3],
                       [4, 5, 6]], dtype=int8)
            R
                array([1, 2], dtype=int8)
            Coordinates only on the left object:
              * y        (y) int16 6B 1 2 3
            Coordinates only on the right object:
              * label    (x) int16 4B 1 2
            """.rstrip()
        )

        actual = formatting.diff_array_repr(da_a, da_b, "equals")
        assert actual == expected

        va = xr.Variable(
            "x", np.array([1, 2, 3], dtype="int64"), {"title": "test Variable"}
        )
        vb = xr.Variable(("x", "y"), np.array([[1, 2, 3], [4, 5, 6]], dtype="int64"))

        expected = dedent(
            """\
        Left and right Variable objects are not equal
        Differing dimensions:
            (x: 3) != (x: 2, y: 3)
        Differing values:
        L
            array([1, 2, 3], dtype=int64)
        R
            array([[1, 2, 3],
                   [4, 5, 6]], dtype=int64)"""
        )

        actual = formatting.diff_array_repr(va, vb, "equals")
        try:
            assert actual == expected
        except AssertionError:
            assert actual == expected.replace(", dtype=int64", "")

    @pytest.mark.filterwarnings("error")
    def test_diff_attrs_repr_with_array(self) -> None:
        attrs_a = {"attr": np.array([0, 1])}

        attrs_b = {"attr": 1}
        expected = dedent(
            """\
            Differing attributes:
            L   attr: [0 1]
            R   attr: 1
            """
        ).strip()
        actual = formatting.diff_attrs_repr(attrs_a, attrs_b, "equals")
        assert expected == actual

        attrs_c = {"attr": np.array([-3, 5])}
        expected = dedent(
            """\
            Differing attributes:
            L   attr: [0 1]
            R   attr: [-3  5]
            """
        ).strip()
        actual = formatting.diff_attrs_repr(attrs_a, attrs_c, "equals")
        assert expected == actual

        # should not raise a warning
        attrs_c = {"attr": np.array([0, 1, 2])}
        expected = dedent(
            """\
            Differing attributes:
            L   attr: [0 1]
            R   attr: [0 1 2]
            """
        ).strip()
        actual = formatting.diff_attrs_repr(attrs_a, attrs_c, "equals")
        assert expected == actual

    def test__diff_mapping_repr_array_attrs_on_variables(self) -> None:
        a = {
            "a": xr.DataArray(
                dims="x",
                data=np.array([1], dtype="int16"),
                attrs={"b": np.array([1, 2], dtype="int8")},
            )
        }
        b = {
            "a": xr.DataArray(
                dims="x",
                data=np.array([1], dtype="int16"),
                attrs={"b": np.array([2, 3], dtype="int8")},
            )
        }
        actual = formatting.diff_data_vars_repr(a, b, compat="identical", col_width=8)
        expected = dedent(
            """\
            Differing data variables:
            L   a   (x) int16 2B 1
                Differing variable attributes:
                    b: [1 2]
            R   a   (x) int16 2B 1
                Differing variable attributes:
                    b: [2 3]
            """.rstrip()
        )

        assert actual == expected

    def test_diff_dataset_repr(self) -> None:
        ds_a = xr.Dataset(
            data_vars={
                "var1": (("x", "y"), np.array([[1, 2, 3], [4, 5, 6]], dtype="int64")),
                "var2": ("x", np.array([3, 4], dtype="int64")),
            },
            coords={
                "x": (
                    "x",
                    np.array(["a", "b"], dtype="U1"),
                    {"foo": "bar", "same": "same"},
                ),
                "y": np.array([1, 2, 3], dtype="int64"),
            },
            attrs={"title": "mytitle", "description": "desc"},
        )

        ds_b = xr.Dataset(
            data_vars={"var1": ("x", np.array([1, 2], dtype="int64"))},
            coords={
                "x": (
                    "x",
                    np.array(["a", "c"], dtype="U1"),
                    {"source": 0, "foo": "baz", "same": "same"},
                ),
                "label": ("x", np.array([1, 2], dtype="int64")),
            },
            attrs={"title": "newtitle"},
        )

        byteorder = "<" if sys.byteorder == "little" else ">"
        expected = dedent(
            f"""\
        Left and right Dataset objects are not identical
        Differing dimensions:
            (x: 2, y: 3) != (x: 2)
        Differing coordinates:
        L * x        (x) {byteorder}U1 8B 'a' 'b'
            Differing variable attributes:
                foo: bar
        R * x        (x) {byteorder}U1 8B 'a' 'c'
            Differing variable attributes:
                source: 0
                foo: baz
        Coordinates only on the left object:
          * y        (y) int64 24B 1 2 3
        Coordinates only on the right object:
            label    (x) int64 16B 1 2
        Differing data variables:
        L   var1     (x, y) int64 48B 1 2 3 4 5 6
        R   var1     (x) int64 16B 1 2
        Data variables only on the left object:
            var2     (x) int64 16B 3 4
        Differing attributes:
        L   title: mytitle
        R   title: newtitle
        Attributes only on the left object:
            description: desc"""
        )

        actual = formatting.diff_dataset_repr(ds_a, ds_b, "identical")
        assert actual == expected

    def test_array_repr(self) -> None:
        ds = xr.Dataset(
            coords={
                "foo": np.array([1, 2, 3], dtype=np.uint64),
                "bar": np.array([1, 2, 3], dtype=np.uint64),
            }
        )
        ds[(1, 2)] = xr.DataArray(np.array([0], dtype=np.uint64), dims="test")
        ds_12 = ds[(1, 2)]

        # Test repr function behaves correctly:
        actual = formatting.array_repr(ds_12)

        expected = dedent(
            """\
        <xarray.DataArray (1, 2) (test: 1)> Size: 8B
        array([0], dtype=uint64)
        Dimensions without coordinates: test"""
        )

        assert actual == expected

        # Test repr, str prints returns correctly as well:
        assert repr(ds_12) == expected
        assert str(ds_12) == expected

        # f-strings (aka format(...)) by default should use the repr:
        actual = f"{ds_12}"
        assert actual == expected

        with xr.set_options(display_expand_data=False):
            actual = formatting.array_repr(ds[(1, 2)])
            expected = dedent(
                """\
            <xarray.DataArray (1, 2) (test: 1)> Size: 8B
            0
            Dimensions without coordinates: test"""
            )

            assert actual == expected

    def test_array_repr_variable(self) -> None:
        var = xr.Variable("x", [0, 1])

        formatting.array_repr(var)

        with xr.set_options(display_expand_data=False):
            formatting.array_repr(var)

    def test_array_repr_recursive(self) -> None:
        # GH:issue:7111

        # direct recursion
        var = xr.Variable("x", [0, 1])
        var.attrs["x"] = var
        formatting.array_repr(var)

        da = xr.DataArray([0, 1], dims=["x"])
        da.attrs["x"] = da
        formatting.array_repr(da)

        # indirect recursion
        var.attrs["x"] = da
        da.attrs["x"] = var
        formatting.array_repr(var)
        formatting.array_repr(da)

    @requires_dask
    def test_array_scalar_format(self) -> None:
        # Test numpy scalars:
        var = xr.DataArray(np.array(0))
        assert format(var, "") == repr(var)
        assert format(var, "d") == "0"
        assert format(var, ".2f") == "0.00"

        # Test dask scalars, not supported however:
        import dask.array as da

        var = xr.DataArray(da.array(0))
        assert format(var, "") == repr(var)
        with pytest.raises(TypeError) as excinfo:
            format(var, ".2f")
        assert "unsupported format string passed to" in str(excinfo.value)

        # Test numpy arrays raises:
        var = xr.DataArray([0.1, 0.2])
        with pytest.raises(NotImplementedError) as excinfo:  # type: ignore[assignment]
            format(var, ".2f")
        assert "Using format_spec is only supported" in str(excinfo.value)

    def test_datatree_print_empty_node(self):
        dt: xr.DataTree = xr.DataTree(name="root")
        printout = str(dt)
        assert printout == "<xarray.DataTree 'root'>\nGroup: /"

    def test_datatree_print_empty_node_with_attrs(self):
        dat = xr.Dataset(attrs={"note": "has attrs"})
        dt: xr.DataTree = xr.DataTree(name="root", dataset=dat)
        printout = str(dt)
        assert printout == dedent(
            """\
            <xarray.DataTree 'root'>
            Group: /
                Attributes:
                    note:     has attrs"""
        )

    def test_datatree_print_node_with_data(self):
        dat = xr.Dataset({"a": [0, 2]})
        dt: xr.DataTree = xr.DataTree(name="root", dataset=dat)
        printout = str(dt)
        expected = [
            "<xarray.DataTree 'root'>",
            "Group: /",
            "Dimensions",
            "Coordinates",
            "a",
        ]
        for expected_line, printed_line in zip(
            expected, printout.splitlines(), strict=True
        ):
            assert expected_line in printed_line

    def test_datatree_printout_nested_node(self):
        dat = xr.Dataset({"a": [0, 2]})
        root = xr.DataTree.from_dict(
            {
                "/results": dat,
            }
        )
        printout = str(root)
        assert printout.splitlines()[3].startswith("    ")

    def test_datatree_repr_of_node_with_data(self):
        dat = xr.Dataset({"a": [0, 2]})
        dt: xr.DataTree = xr.DataTree(name="root", dataset=dat)
        assert "Coordinates" in repr(dt)

    def test_diff_datatree_repr_different_groups(self):
        dt_1: xr.DataTree = xr.DataTree.from_dict({"a": None})
        dt_2: xr.DataTree = xr.DataTree.from_dict({"b": None})

        expected = dedent(
            """\
            Left and right DataTree objects are not identical

            Children at root node do not match: ['a'] vs ['b']"""
        )
        actual = formatting.diff_datatree_repr(dt_1, dt_2, "identical")
        assert actual == expected

    def test_diff_datatree_repr_different_subgroups(self):
        dt_1: xr.DataTree = xr.DataTree.from_dict({"a": None, "a/b": None, "a/c": None})
        dt_2: xr.DataTree = xr.DataTree.from_dict({"a": None, "a/b": None})

        expected = dedent(
            """\
            Left and right DataTree objects are not isomorphic

            Children at node 'a' do not match: ['b', 'c'] vs ['b']"""
        )
        actual = formatting.diff_datatree_repr(dt_1, dt_2, "isomorphic")
        assert actual == expected

    def test_diff_datatree_repr_node_data(self):
        # casting to int64 explicitly ensures that int64s are created on all architectures
        ds1 = xr.Dataset({"u": np.int64(0), "v": np.int64(1)})
        ds3 = xr.Dataset({"w": np.int64(5)})
        dt_1: xr.DataTree = xr.DataTree.from_dict({"a": ds1, "a/b": ds3})
        ds2 = xr.Dataset({"u": np.int64(0)})
        ds4 = xr.Dataset({"w": np.int64(6)})
        dt_2: xr.DataTree = xr.DataTree.from_dict({"a": ds2, "a/b": ds4}, name="foo")

        expected = dedent(
            """\
            Left and right DataTree objects are not identical

            Differing names:
                None != 'foo'

            Data at node 'a' does not match:
                Data variables only on the left object:
                    v        int64 8B 1

            Data at node 'a/b' does not match:
                Differing data variables:
                L   w        int64 8B 5
                R   w        int64 8B 6"""
        )
        actual = formatting.diff_datatree_repr(dt_1, dt_2, "identical")
        assert actual == expected

    def test_diff_datatree_repr_equals(self) -> None:
        ds1 = xr.Dataset(data_vars={"data": ("y", [5, 2])})
        ds2 = xr.Dataset(data_vars={"data": (("x", "y"), [[5, 2]])})
        dt1 = xr.DataTree.from_dict({"node": ds1})
        dt2 = xr.DataTree.from_dict({"node": ds2})

        expected = dedent(
            """\
            Left and right DataTree objects are not equal

            Data at node 'node' does not match:
                Differing dimensions:
                    (y: 2) != (x: 1, y: 2)
                Differing data variables:
                L   data     (y) int64 16B 5 2
                R   data     (x, y) int64 16B 5 2"""
        )

        actual = formatting.diff_datatree_repr(dt1, dt2, "equals")
        assert actual == expected


def test_inline_variable_array_repr_custom_repr() -> None:
    class CustomArray:
        def __init__(self, value, attr):
            self.value = value
            self.attr = attr

        def _repr_inline_(self, width):
            formatted = f"({self.attr}) {self.value}"
            if len(formatted) > width:
                formatted = f"({self.attr}) ..."

            return formatted

        def __array_namespace__(self, *args, **kwargs):
            return NotImplemented

        @property
        def shape(self) -> tuple[int, ...]:
            return self.value.shape

        @property
        def dtype(self):
            return self.value.dtype

        @property
        def ndim(self):
            return self.value.ndim

    value = CustomArray(np.array([20, 40]), "m")
    variable = xr.Variable("x", value)

    max_width = 10
    actual = formatting.inline_variable_array_repr(variable, max_width=10)

    assert actual == value._repr_inline_(max_width)


def test_set_numpy_options() -> None:
    original_options = np.get_printoptions()
    with formatting.set_numpy_options(threshold=10):
        assert len(repr(np.arange(500))) < 200
    # original options are restored
    assert np.get_printoptions() == original_options


def test_short_array_repr() -> None:
    cases = [
        np.random.randn(500),
        np.random.randn(20, 20),
        np.random.randn(5, 10, 15),
        np.random.randn(5, 10, 15, 3),
        np.random.randn(100, 5, 1),
    ]
    # number of lines:
    # for default numpy repr: 167, 140, 254, 248, 599
    # for short_array_repr: 1, 7, 24, 19, 25
    for array in cases:
        num_lines = formatting.short_array_repr(array).count("\n") + 1
        assert num_lines < 30

    # threshold option (default: 200)
    array2 = np.arange(100)
    assert "..." not in formatting.short_array_repr(array2)
    with xr.set_options(display_values_threshold=10):
        assert "..." in formatting.short_array_repr(array2)


def test_large_array_repr_length() -> None:
    da = xr.DataArray(np.random.randn(100, 5, 1))

    result = repr(da).splitlines()
    assert len(result) < 50


@requires_netCDF4
def test_repr_file_collapsed(tmp_path) -> None:
    arr_to_store = xr.DataArray(np.arange(300, dtype=np.int64), dims="test")
    arr_to_store.to_netcdf(tmp_path / "test.nc", engine="netcdf4")

    with (
        xr.open_dataarray(tmp_path / "test.nc") as arr,
        xr.set_options(display_expand_data=False),
    ):
        actual = repr(arr)
        expected = dedent(
            """\
        <xarray.DataArray (test: 300)> Size: 2kB
        [300 values with dtype=int64]
        Dimensions without coordinates: test"""
        )

        assert actual == expected

        arr_loaded = arr.compute()
        actual = arr_loaded.__repr__()
        expected = dedent(
            """\
        <xarray.DataArray (test: 300)> Size: 2kB
        0 1 2 3 4 5 6 7 8 9 10 11 12 ... 288 289 290 291 292 293 294 295 296 297 298 299
        Dimensions without coordinates: test"""
        )

        assert actual == expected


@pytest.mark.parametrize(
    "display_max_rows, n_vars, n_attr",
    [(50, 40, 30), (35, 40, 30), (11, 40, 30), (1, 40, 30)],
)
def test__mapping_repr(display_max_rows, n_vars, n_attr) -> None:
    long_name = "long_name"
    a = np.char.add(long_name, np.arange(0, n_vars).astype(str))
    b = np.char.add("attr_", np.arange(0, n_attr).astype(str))
    c = np.char.add("coord", np.arange(0, n_vars).astype(str))
    attrs = dict.fromkeys(b, 2)
    coords = {_c: np.array([0, 1], dtype=np.uint64) for _c in c}
    data_vars = dict()
    for v, _c in zip(a, coords.items(), strict=True):
        data_vars[v] = xr.DataArray(
            name=v,
            data=np.array([3, 4], dtype=np.uint64),
            dims=[_c[0]],
            coords=dict([_c]),
        )

    ds = xr.Dataset(data_vars)
    ds.attrs = attrs

    with xr.set_options(display_max_rows=display_max_rows):
        # Parse the data_vars print and show only data_vars rows:
        summary = formatting.dataset_repr(ds).split("\n")
        summary = [v for v in summary if long_name in v]
        # The length should be less than or equal to display_max_rows:
        len_summary = len(summary)
        data_vars_print_size = min(display_max_rows, len_summary)
        assert len_summary == data_vars_print_size

        summary = formatting.data_vars_repr(ds.data_vars).split("\n")
        summary = [v for v in summary if long_name in v]
        # The length should be equal to the number of data variables
        len_summary = len(summary)
        assert len_summary == n_vars

        summary = formatting.coords_repr(ds.coords).split("\n")
        summary = [v for v in summary if "coord" in v]
        # The length should be equal to the number of data variables
        len_summary = len(summary)
        assert len_summary == n_vars

    with xr.set_options(
        display_max_rows=display_max_rows,
        display_expand_coords=False,
        display_expand_data_vars=False,
        display_expand_attrs=False,
    ):
        actual = formatting.dataset_repr(ds)
        col_width = formatting._calculate_col_width(ds.variables)
        dims_start = formatting.pretty_print("Dimensions:", col_width)
        dims_values = formatting.dim_summary_limited(
            ds.sizes, col_width=col_width + 1, max_rows=display_max_rows
        )
        expected_size = "1kB"
        expected = f"""\
<xarray.Dataset> Size: {expected_size}
{dims_start}({dims_values})
Coordinates: ({n_vars})
Data variables: ({n_vars})
Attributes: ({n_attr})"""
        expected = dedent(expected)
        assert actual == expected


def test__mapping_repr_recursive() -> None:
    # GH:issue:7111

    # direct recursion
    ds = xr.Dataset({"a": ("x", [1, 2, 3])})
    ds.attrs["ds"] = ds
    formatting.dataset_repr(ds)

    # indirect recursion
    ds2 = xr.Dataset({"b": ("y", [1, 2, 3])})
    ds.attrs["ds"] = ds2
    ds2.attrs["ds"] = ds
    formatting.dataset_repr(ds2)


def test__element_formatter(n_elements: int = 100) -> None:
    expected = """\
    Dimensions without coordinates: dim_0: 3, dim_1: 3, dim_2: 3, dim_3: 3,
                                    dim_4: 3, dim_5: 3, dim_6: 3, dim_7: 3,
                                    dim_8: 3, dim_9: 3, dim_10: 3, dim_11: 3,
                                    dim_12: 3, dim_13: 3, dim_14: 3, dim_15: 3,
                                    dim_16: 3, dim_17: 3, dim_18: 3, dim_19: 3,
                                    dim_20: 3, dim_21: 3, dim_22: 3, dim_23: 3,
                                    ...
                                    dim_76: 3, dim_77: 3, dim_78: 3, dim_79: 3,
                                    dim_80: 3, dim_81: 3, dim_82: 3, dim_83: 3,
                                    dim_84: 3, dim_85: 3, dim_86: 3, dim_87: 3,
                                    dim_88: 3, dim_89: 3, dim_90: 3, dim_91: 3,
                                    dim_92: 3, dim_93: 3, dim_94: 3, dim_95: 3,
                                    dim_96: 3, dim_97: 3, dim_98: 3, dim_99: 3"""
    expected = dedent(expected)

    intro = "Dimensions without coordinates: "
    elements = [
        f"{k}: {v}" for k, v in {f"dim_{k}": 3 for k in np.arange(n_elements)}.items()
    ]
    values = xr.core.formatting._element_formatter(
        elements, col_width=len(intro), max_rows=12
    )
    actual = intro + values
    assert expected == actual


def test_lazy_array_wont_compute() -> None:
    from xarray.core.indexing import LazilyIndexedArray

    class LazilyIndexedArrayNotComputable(LazilyIndexedArray):
        def __array__(
            self, dtype: np.typing.DTypeLike = None, /, *, copy: bool | None = None
        ) -> np.ndarray:
            raise NotImplementedError("Computing this array is not possible.")

    arr = LazilyIndexedArrayNotComputable(np.array([1, 2]))
    var = xr.DataArray(arr)

    # These will crash if var.data are converted to numpy arrays:
    var.__repr__()
    var._repr_html_()


@pytest.mark.parametrize("as_dataset", (False, True))
def test_format_xindexes_none(as_dataset: bool) -> None:
    # ensure repr for empty xindexes can be displayed #8367

    expected = """\
    Indexes:
        *empty*"""
    expected = dedent(expected)

    obj: xr.DataArray | xr.Dataset = xr.DataArray()
    obj = obj._to_temp_dataset() if as_dataset else obj

    actual = repr(obj.xindexes)
    assert actual == expected


@pytest.mark.parametrize("as_dataset", (False, True))
def test_format_xindexes(as_dataset: bool) -> None:
    expected = """\
    Indexes:
        x        PandasIndex"""
    expected = dedent(expected)

    obj: xr.DataArray | xr.Dataset = xr.DataArray([1], coords={"x": [1]})
    obj = obj._to_temp_dataset() if as_dataset else obj

    actual = repr(obj.xindexes)
    assert actual == expected


@requires_cftime
def test_empty_cftimeindex_repr() -> None:
    index = xr.coding.cftimeindex.CFTimeIndex([])

    expected = """\
    Indexes:
        time     CFTimeIndex([], dtype='object', length=0, calendar=None, freq=None)"""
    expected = dedent(expected)

    da = xr.DataArray([], coords={"time": index})

    actual = repr(da.indexes)
    assert actual == expected


def test_display_nbytes() -> None:
    xds = xr.Dataset(
        {
            "foo": np.arange(1200, dtype=np.int16),
            "bar": np.arange(111, dtype=np.int16),
        }
    )

    # Note: int16 is used to ensure that dtype is shown in the
    # numpy array representation for all OSes included Windows

    actual = repr(xds)
    expected = """
<xarray.Dataset> Size: 3kB
Dimensions:  (foo: 1200, bar: 111)
Coordinates:
  * foo      (foo) int16 2kB 0 1 2 3 4 5 6 ... 1194 1195 1196 1197 1198 1199
  * bar      (bar) int16 222B 0 1 2 3 4 5 6 7 ... 104 105 106 107 108 109 110
Data variables:
    *empty*
    """.strip()
    assert actual == expected

    actual = repr(xds["foo"])
    array_repr = repr(xds.foo.data).replace("\n     ", "")
    expected = f"""
<xarray.DataArray 'foo' (foo: 1200)> Size: 2kB
{array_repr}
Coordinates:
  * foo      (foo) int16 2kB 0 1 2 3 4 5 6 ... 1194 1195 1196 1197 1198 1199
""".strip()
    assert actual == expected


def test_array_repr_dtypes():
    # These dtypes are expected to be represented similarly
    # on Ubuntu, macOS and Windows environments of the CI.
    # Unsigned integer could be used as easy replacements
    # for tests where the data-type does not matter,
    # but the repr does, including the size
    # (size of a int == size of an uint)

    # Signed integer dtypes

    ds = xr.DataArray(np.array([0], dtype="int8"), dims="x")
    actual = repr(ds)
    expected = """
<xarray.DataArray (x: 1)> Size: 1B
array([0], dtype=int8)
Dimensions without coordinates: x
        """.strip()
    assert actual == expected

    ds = xr.DataArray(np.array([0], dtype="int16"), dims="x")
    actual = repr(ds)
    expected = """
<xarray.DataArray (x: 1)> Size: 2B
array([0], dtype=int16)
Dimensions without coordinates: x
        """.strip()
    assert actual == expected

    # Unsigned integer dtypes

    ds = xr.DataArray(np.array([0], dtype="uint8"), dims="x")
    actual = repr(ds)
    expected = """
<xarray.DataArray (x: 1)> Size: 1B
array([0], dtype=uint8)
Dimensions without coordinates: x
        """.strip()
    assert actual == expected

    ds = xr.DataArray(np.array([0], dtype="uint16"), dims="x")
    actual = repr(ds)
    expected = """
<xarray.DataArray (x: 1)> Size: 2B
array([0], dtype=uint16)
Dimensions without coordinates: x
        """.strip()
    assert actual == expected

    ds = xr.DataArray(np.array([0], dtype="uint32"), dims="x")
    actual = repr(ds)
    expected = """
<xarray.DataArray (x: 1)> Size: 4B
array([0], dtype=uint32)
Dimensions without coordinates: x
        """.strip()
    assert actual == expected

    ds = xr.DataArray(np.array([0], dtype="uint64"), dims="x")
    actual = repr(ds)
    expected = """
<xarray.DataArray (x: 1)> Size: 8B
array([0], dtype=uint64)
Dimensions without coordinates: x
        """.strip()
    assert actual == expected

    # Float dtypes

    ds = xr.DataArray(np.array([0.0]), dims="x")
    actual = repr(ds)
    expected = """
<xarray.DataArray (x: 1)> Size: 8B
array([0.])
Dimensions without coordinates: x
        """.strip()
    assert actual == expected

    ds = xr.DataArray(np.array([0], dtype="float16"), dims="x")
    actual = repr(ds)
    expected = """
<xarray.DataArray (x: 1)> Size: 2B
array([0.], dtype=float16)
Dimensions without coordinates: x
        """.strip()
    assert actual == expected

    ds = xr.DataArray(np.array([0], dtype="float32"), dims="x")
    actual = repr(ds)
    expected = """
<xarray.DataArray (x: 1)> Size: 4B
array([0.], dtype=float32)
Dimensions without coordinates: x
        """.strip()
    assert actual == expected

    ds = xr.DataArray(np.array([0], dtype="float64"), dims="x")
    actual = repr(ds)
    expected = """
<xarray.DataArray (x: 1)> Size: 8B
array([0.])
Dimensions without coordinates: x
        """.strip()
    assert actual == expected

    # Signed integer dtypes

    array = np.array([0])
    ds = xr.DataArray(array, dims="x")
    actual = repr(ds)
    expected = f"""
<xarray.DataArray (x: 1)> Size: {array.dtype.itemsize}B
{array!r}
Dimensions without coordinates: x
        """.strip()
    assert actual == expected

    array = np.array([0], dtype="int32")
    ds = xr.DataArray(array, dims="x")
    actual = repr(ds)
    expected = f"""
<xarray.DataArray (x: 1)> Size: 4B
{array!r}
Dimensions without coordinates: x
        """.strip()
    assert actual == expected

    array = np.array([0], dtype="int64")
    ds = xr.DataArray(array, dims="x")
    actual = repr(ds)
    expected = f"""
<xarray.DataArray (x: 1)> Size: 8B
{array!r}
Dimensions without coordinates: x
        """.strip()
    assert actual == expected


def test_repr_pandas_range_index() -> None:
    # lazy data repr but values shown in inline repr
    xidx = xr.indexes.PandasIndex(pd.RangeIndex(10), "x")
    ds = xr.Dataset(coords=xr.Coordinates.from_xindex(xidx))
    actual = repr(ds.x)
    expected = """
<xarray.DataArray 'x' (x: 10)> Size: 80B
[10 values with dtype=int64]
Coordinates:
  * x        (x) int64 80B 0 1 2 3 4 5 6 7 8 9
    """.strip()
    assert actual == expected


def test_repr_pandas_multi_index() -> None:
    # lazy data repr but values shown in inline repr
    midx = pd.MultiIndex.from_product([["a", "b"], [1, 2]], names=["foo", "bar"])
    coords = xr.Coordinates.from_pandas_multiindex(midx, "x")
    ds = xr.Dataset(coords=coords)

    actual = repr(ds.x)
    expected = """
<xarray.DataArray 'x' (x: 4)> Size: 32B
[4 values with dtype=object]
Coordinates:
  * x        (x) object 32B MultiIndex
  * foo      (x) object 32B 'a' 'a' 'b' 'b'
  * bar      (x) int64 32B 1 2 1 2
    """.strip()
    assert actual == expected

    actual = repr(ds.foo)
    expected = """
<xarray.DataArray 'foo' (x: 4)> Size: 32B
[4 values with dtype=object]
Coordinates:
  * x        (x) object 32B MultiIndex
  * foo      (x) object 32B 'a' 'a' 'b' 'b'
  * bar      (x) int64 32B 1 2 1 2
    """.strip()
    assert actual == expected
