from __future__ import annotations

import re
import warnings

import numpy as np
import pandas as pd
import pytest

import xarray as xr
from xarray.core import dtypes
from xarray.core.options import set_options
from xarray.structure import merge
from xarray.structure.merge import MergeError
from xarray.testing import assert_equal, assert_identical
from xarray.tests.test_dataset import create_test_data


class TestMergeInternals:
    def test_broadcast_dimension_size(self):
        actual = merge.broadcast_dimension_size(
            [xr.Variable("x", [1]), xr.Variable("y", [2, 1])]
        )
        assert actual == {"x": 1, "y": 2}

        actual = merge.broadcast_dimension_size(
            [xr.Variable(("x", "y"), [[1, 2]]), xr.Variable("y", [2, 1])]
        )
        assert actual == {"x": 1, "y": 2}

        with pytest.raises(ValueError):
            merge.broadcast_dimension_size(
                [xr.Variable(("x", "y"), [[1, 2]]), xr.Variable("y", [2])]
            )


class TestMergeFunction:
    def test_merge_arrays(self):
        data = create_test_data(add_attrs=False)

        actual = xr.merge([data.var1, data.var2])
        expected = data[["var1", "var2"]]
        assert_identical(actual, expected)

    @pytest.mark.parametrize("use_new_combine_kwarg_defaults", [True, False])
    def test_merge_datasets(self, use_new_combine_kwarg_defaults):
        with set_options(use_new_combine_kwarg_defaults=use_new_combine_kwarg_defaults):
            data = create_test_data(add_attrs=False, use_extension_array=True)

            actual = xr.merge([data[["var1"]], data[["var2"]]])
            expected = data[["var1", "var2"]]
            assert_identical(actual, expected)

            actual = xr.merge([data, data])
            assert_identical(actual, data)

    def test_merge_dataarray_unnamed(self):
        data = xr.DataArray([1, 2], dims="x")
        with pytest.raises(ValueError, match=r"without providing an explicit name"):
            xr.merge([data])

    def test_merge_arrays_attrs_default(self):
        var1_attrs = {"a": 1, "b": 2}
        var2_attrs = {"a": 1, "c": 3}
        expected_attrs = {"a": 1, "b": 2}

        data = create_test_data(add_attrs=False)
        expected = data[["var1", "var2"]].copy()
        expected.var1.attrs = var1_attrs
        expected.var2.attrs = var2_attrs
        expected.attrs = expected_attrs

        data.var1.attrs = var1_attrs
        data.var2.attrs = var2_attrs
        actual = xr.merge([data.var1, data.var2])
        assert_identical(actual, expected)

    @pytest.mark.parametrize(
        "combine_attrs, var1_attrs, var2_attrs, expected_attrs, expect_exception",
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
            (
                "drop_conflicts",
                {"a": 1, "b": np.array([2]), "c": np.array([3])},
                {"b": 1, "c": np.array([3]), "d": 4},
                {"a": 1, "c": np.array([3]), "d": 4},
                False,
            ),
            (
                lambda attrs, context: attrs[1],
                {"a": 1, "b": 2, "c": 3},
                {"a": 4, "b": 3, "c": 1},
                {"a": 4, "b": 3, "c": 1},
                False,
            ),
        ],
    )
    def test_merge_arrays_attrs(
        self, combine_attrs, var1_attrs, var2_attrs, expected_attrs, expect_exception
    ):
        data1 = xr.Dataset(attrs=var1_attrs)
        data2 = xr.Dataset(attrs=var2_attrs)
        if expect_exception:
            with pytest.raises(MergeError, match="combine_attrs"):
                actual = xr.merge([data1, data2], combine_attrs=combine_attrs)
        else:
            actual = xr.merge([data1, data2], combine_attrs=combine_attrs)
            expected = xr.Dataset(attrs=expected_attrs)

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
            (
                lambda attrs, context: attrs[1],
                {"a": 1, "b": 2, "c": 3},
                {"a": 4, "b": 3, "c": 1},
                {"a": 4, "b": 3, "c": 1},
                False,
            ),
        ],
    )
    def test_merge_arrays_attrs_variables(
        self, combine_attrs, attrs1, attrs2, expected_attrs, expect_exception
    ):
        """check that combine_attrs is used on data variables and coords"""
        input_attrs1 = attrs1.copy()
        data1 = xr.Dataset(
            {"var1": ("dim1", [], attrs1)}, coords={"dim1": ("dim1", [], attrs1)}
        )
        input_attrs2 = attrs2.copy()
        data2 = xr.Dataset(
            {"var1": ("dim1", [], attrs2)}, coords={"dim1": ("dim1", [], attrs2)}
        )

        if expect_exception:
            with pytest.raises(MergeError, match="combine_attrs"):
                with pytest.warns(
                    FutureWarning,
                    match="will change from compat='no_conflicts' to compat='override'",
                ):
                    actual = xr.merge([data1, data2], combine_attrs=combine_attrs)
        else:
            actual = xr.merge(
                [data1, data2], compat="no_conflicts", combine_attrs=combine_attrs
            )
            expected = xr.Dataset(
                {"var1": ("dim1", [], expected_attrs)},
                coords={"dim1": ("dim1", [], expected_attrs)},
            )

            assert_identical(actual, expected)

            # Check also that input attributes weren't modified
            assert data1["var1"].attrs == input_attrs1
            assert data1.coords["dim1"].attrs == input_attrs1
            assert data2["var1"].attrs == input_attrs2
            assert data2.coords["dim1"].attrs == input_attrs2

    def test_merge_attrs_override_copy(self):
        ds1 = xr.Dataset(attrs={"x": 0})
        ds2 = xr.Dataset(attrs={"x": 1})
        ds3 = xr.merge([ds1, ds2], combine_attrs="override")
        ds3.attrs["x"] = 2
        assert ds1.x == 0

    def test_merge_attrs_drop_conflicts(self):
        ds1 = xr.Dataset(attrs={"a": 0, "b": 0, "c": 0})
        ds2 = xr.Dataset(attrs={"b": 0, "c": 1, "d": 0})
        ds3 = xr.Dataset(attrs={"a": 0, "b": 1, "c": 0, "e": 0})

        actual = xr.merge([ds1, ds2, ds3], combine_attrs="drop_conflicts")
        expected = xr.Dataset(attrs={"a": 0, "d": 0, "e": 0})
        assert_identical(actual, expected)

    def test_merge_attrs_drop_conflicts_numpy_arrays(self):
        """Test drop_conflicts with numpy arrays."""
        # Test with numpy arrays (which return arrays from ==)
        arr1 = np.array([1, 2, 3])
        arr2 = np.array([1, 2, 3])
        arr3 = np.array([4, 5, 6])

        ds1 = xr.Dataset(attrs={"arr": arr1, "scalar": 1})
        ds2 = xr.Dataset(attrs={"arr": arr2, "scalar": 1})  # Same array values
        ds3 = xr.Dataset(attrs={"arr": arr3, "other": 2})  # Different array values

        # Arrays are considered equivalent if they have the same values
        actual = xr.merge([ds1, ds2], combine_attrs="drop_conflicts")
        assert "arr" in actual.attrs  # Should keep the array since they're equivalent
        assert actual.attrs["scalar"] == 1

        # Different arrays cause the attribute to be dropped
        actual = xr.merge([ds1, ds3], combine_attrs="drop_conflicts")
        assert "arr" not in actual.attrs  # Should drop due to conflict
        assert "other" in actual.attrs

    def test_merge_attrs_drop_conflicts_custom_eq_returns_array(self):
        """Test drop_conflicts with custom objects that return arrays from __eq__."""

        # Test with custom objects that return non-bool from __eq__
        class CustomEq:
            """Object whose __eq__ returns a non-bool value."""

            def __init__(self, value):
                self.value = value

            def __eq__(self, other):
                if not isinstance(other, CustomEq):
                    return False
                # Return a numpy array (truthy if all elements are non-zero)
                return np.array([self.value == other.value])

            def __repr__(self):
                return f"CustomEq({self.value})"

        obj1 = CustomEq(42)
        obj2 = CustomEq(42)  # Same value
        obj3 = CustomEq(99)  # Different value

        ds4 = xr.Dataset(attrs={"custom": obj1, "x": 1})
        ds5 = xr.Dataset(attrs={"custom": obj2, "x": 1})
        ds6 = xr.Dataset(attrs={"custom": obj3, "y": 2})

        # Suppress DeprecationWarning from numpy < 2.0 about ambiguous truth values
        # when our custom __eq__ returns arrays that are evaluated in boolean context
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=DeprecationWarning)

            # Objects returning arrays are dropped (non-boolean return)
            actual = xr.merge([ds4, ds5], combine_attrs="drop_conflicts")
            assert "custom" not in actual.attrs  # Dropped - returns array, not bool
            assert actual.attrs["x"] == 1

            # Different values also dropped (returns array, not bool)
            actual = xr.merge([ds4, ds6], combine_attrs="drop_conflicts")
            assert "custom" not in actual.attrs  # Dropped - returns non-boolean
            assert actual.attrs["x"] == 1
            assert actual.attrs["y"] == 2

    def test_merge_attrs_drop_conflicts_ambiguous_array_returns(self):
        """Test drop_conflicts with objects returning ambiguous arrays from __eq__."""

        # Test edge case: object whose __eq__ returns empty array (ambiguous truth value)
        class EmptyArrayEq:
            def __eq__(self, other):
                if not isinstance(other, EmptyArrayEq):
                    return False
                return np.array([])  # Empty array has ambiguous truth value

            def __repr__(self):
                return "EmptyArrayEq()"

        empty_obj1 = EmptyArrayEq()
        empty_obj2 = EmptyArrayEq()

        ds7 = xr.Dataset(attrs={"empty": empty_obj1})
        ds8 = xr.Dataset(attrs={"empty": empty_obj2})

        # With new behavior: ambiguous truth values are treated as non-equivalent
        # So the attribute is dropped instead of raising an error
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=DeprecationWarning)
            actual = xr.merge([ds7, ds8], combine_attrs="drop_conflicts")
            assert "empty" not in actual.attrs  # Dropped due to ambiguous comparison

        # Test with object that returns multi-element array (also ambiguous)
        class MultiArrayEq:
            def __eq__(self, other):
                if not isinstance(other, MultiArrayEq):
                    return False
                return np.array([True, False])  # Multi-element array is ambiguous

            def __repr__(self):
                return "MultiArrayEq()"

        multi_obj1 = MultiArrayEq()
        multi_obj2 = MultiArrayEq()

        ds9 = xr.Dataset(attrs={"multi": multi_obj1})
        ds10 = xr.Dataset(attrs={"multi": multi_obj2})

        # With new behavior: ambiguous arrays are treated as non-equivalent
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=DeprecationWarning)
            actual = xr.merge([ds9, ds10], combine_attrs="drop_conflicts")
            assert "multi" not in actual.attrs  # Dropped due to ambiguous comparison

    def test_merge_attrs_drop_conflicts_all_true_array(self):
        """Test drop_conflicts with all-True multi-element array from __eq__."""

        # Test with all-True multi-element array (unambiguous truthy)
        class AllTrueArrayEq:
            def __eq__(self, other):
                if not isinstance(other, AllTrueArrayEq):
                    return False
                return np.array([True, True, True])  # All True, but still multi-element

            def __repr__(self):
                return "AllTrueArrayEq()"

        alltrue1 = AllTrueArrayEq()
        alltrue2 = AllTrueArrayEq()

        ds11 = xr.Dataset(attrs={"alltrue": alltrue1})
        ds12 = xr.Dataset(attrs={"alltrue": alltrue2})

        # Multi-element arrays are ambiguous even if all True
        actual = xr.merge([ds11, ds12], combine_attrs="drop_conflicts")
        assert "alltrue" not in actual.attrs  # Dropped due to ambiguous comparison

    def test_merge_attrs_drop_conflicts_nested_arrays(self):
        """Test drop_conflicts with NumPy object arrays containing nested arrays."""
        # Test 1: NumPy object arrays with nested arrays
        # These can have complex comparison behavior
        x = np.array([None], dtype=object)
        x[0] = np.arange(3)
        y = np.array([None], dtype=object)
        y[0] = np.arange(10, 13)

        ds1 = xr.Dataset(attrs={"nested_array": x, "common": 1})
        ds2 = xr.Dataset(attrs={"nested_array": y, "common": 1})

        # Different nested arrays should cause attribute to be dropped
        actual = xr.merge([ds1, ds2], combine_attrs="drop_conflicts")
        assert (
            "nested_array" not in actual.attrs
        )  # Dropped due to different nested arrays
        assert actual.attrs["common"] == 1

        # Test with identical nested arrays
        # Note: Even identical nested arrays will be dropped because comparison
        # raises ValueError due to ambiguous truth value
        z = np.array([None], dtype=object)
        z[0] = np.arange(3)  # Same as x
        ds3 = xr.Dataset(attrs={"nested_array": z, "other": 2})

        actual = xr.merge([ds1, ds3], combine_attrs="drop_conflicts")
        assert (
            "nested_array" not in actual.attrs
        )  # Dropped due to ValueError in comparison
        assert actual.attrs["other"] == 2

    def test_merge_attrs_drop_conflicts_dataset_attrs(self):
        """Test drop_conflicts with xarray.Dataset objects as attributes."""
        # xarray.Dataset objects as attributes (raises TypeError in equivalent)
        attr_ds1 = xr.Dataset({"foo": 1})
        attr_ds2 = xr.Dataset({"bar": 1})  # Different dataset
        attr_ds3 = xr.Dataset({"foo": 1})  # Same as attr_ds1

        ds4 = xr.Dataset(attrs={"dataset_attr": attr_ds1, "scalar": 42})
        ds5 = xr.Dataset(attrs={"dataset_attr": attr_ds2, "scalar": 42})
        ds6 = xr.Dataset(attrs={"dataset_attr": attr_ds3, "other": 99})

        # Different datasets raise TypeError and should be dropped
        actual = xr.merge([ds4, ds5], combine_attrs="drop_conflicts")
        assert "dataset_attr" not in actual.attrs  # Dropped due to TypeError
        assert actual.attrs["scalar"] == 42

        # Identical datasets are also dropped (comparison returns Dataset, not bool)
        actual = xr.merge([ds4, ds6], combine_attrs="drop_conflicts")
        assert "dataset_attr" not in actual.attrs  # Dropped - returns Dataset, not bool
        assert actual.attrs["other"] == 99

    def test_merge_attrs_drop_conflicts_pandas_series(self):
        """Test drop_conflicts with Pandas Series as attributes."""
        # Pandas Series (raises ValueError due to ambiguous truth value)
        series1 = pd.Series([1, 2])
        series2 = pd.Series([3, 4])  # Different values
        series3 = pd.Series([1, 2])  # Same as series1

        ds7 = xr.Dataset(attrs={"series": series1, "value": "a"})
        ds8 = xr.Dataset(attrs={"series": series2, "value": "a"})
        ds9 = xr.Dataset(attrs={"series": series3, "value": "a"})

        # Suppress potential warnings from pandas comparisons
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=DeprecationWarning)
            warnings.filterwarnings("ignore", category=FutureWarning)

            # Different series raise ValueError and get dropped
            actual = xr.merge([ds7, ds8], combine_attrs="drop_conflicts")
            assert "series" not in actual.attrs  # Dropped due to ValueError
            assert actual.attrs["value"] == "a"

            # Even identical series raise ValueError in equivalent() and get dropped
            # because Series comparison returns another Series with ambiguous truth value
            actual = xr.merge([ds7, ds9], combine_attrs="drop_conflicts")
            assert "series" not in actual.attrs  # Dropped due to ValueError
            assert actual.attrs["value"] == "a"

    def test_merge_attrs_drop_conflicts_eq_returns_string(self):
        """Test objects whose __eq__ returns strings are dropped."""

        # Case 1: Objects whose __eq__ returns non-boolean strings
        class ReturnsString:
            def __init__(self, value):
                self.value = value

            def __eq__(self, other):
                # Always returns a string (non-boolean)
                return "comparison result"

        obj1 = ReturnsString("A")
        obj2 = ReturnsString("B")  # Different object

        ds1 = xr.Dataset(attrs={"obj": obj1})
        ds2 = xr.Dataset(attrs={"obj": obj2})

        actual = xr.merge([ds1, ds2], combine_attrs="drop_conflicts")

        # Strict behavior: drops attribute because __eq__ returns non-boolean
        assert "obj" not in actual.attrs

    def test_merge_attrs_drop_conflicts_eq_returns_number(self):
        """Test objects whose __eq__ returns numbers are dropped."""

        # Case 2: Objects whose __eq__ returns numbers
        class ReturnsZero:
            def __init__(self, value):
                self.value = value

            def __eq__(self, other):
                # Always returns 0 (non-boolean)
                return 0

        obj3 = ReturnsZero("same")
        obj4 = ReturnsZero("same")  # Different object, same value

        ds3 = xr.Dataset(attrs={"zero": obj3})
        ds4 = xr.Dataset(attrs={"zero": obj4})

        actual = xr.merge([ds3, ds4], combine_attrs="drop_conflicts")

        # Strict behavior: drops attribute because __eq__ returns non-boolean
        assert "zero" not in actual.attrs

    def test_merge_attrs_no_conflicts_compat_minimal(self):
        """make sure compat="minimal" does not silence errors"""
        ds1 = xr.Dataset({"a": ("x", [], {"a": 0})})
        ds2 = xr.Dataset({"a": ("x", [], {"a": 1})})

        with pytest.raises(xr.MergeError, match="combine_attrs"):
            xr.merge([ds1, ds2], combine_attrs="no_conflicts", compat="minimal")

    def test_merge_dicts_simple(self):
        actual = xr.merge([{"foo": 0}, {"bar": "one"}, {"baz": 3.5}])
        expected = xr.Dataset({"foo": 0, "bar": "one", "baz": 3.5})
        assert_identical(actual, expected)

    def test_merge_dicts_dims(self):
        actual = xr.merge([{"y": ("x", [13])}, {"x": [12]}])
        expected = xr.Dataset({"x": [12], "y": ("x", [13])})
        assert_identical(actual, expected)

    def test_merge_coordinates(self):
        coords1 = xr.Coordinates({"x": ("x", [0, 1, 2])})
        coords2 = xr.Coordinates({"y": ("y", [3, 4, 5])})
        expected = xr.Dataset(coords={"x": [0, 1, 2], "y": [3, 4, 5]})
        actual = xr.merge([coords1, coords2])
        assert_identical(actual, expected)

    def test_merge_error(self):
        ds = xr.Dataset({"x": 0})
        with pytest.raises(xr.MergeError):
            xr.merge([ds, ds + 1], compat="no_conflicts")

    def test_merge_alignment_error(self):
        ds = xr.Dataset(coords={"x": [1, 2]})
        other = xr.Dataset(coords={"x": [2, 3]})
        with pytest.raises(ValueError, match=r"cannot align.*join.*exact.*not equal.*"):
            xr.merge([ds, other], join="exact")

    def test_merge_wrong_input_error(self):
        with pytest.raises(TypeError, match=r"objects must be an iterable"):
            xr.merge([1])  # type: ignore[list-item]
        ds = xr.Dataset(coords={"x": [1, 2]})
        with pytest.raises(TypeError, match=r"objects must be an iterable"):
            xr.merge({"a": ds})  # type: ignore[dict-item]
        with pytest.raises(TypeError, match=r"objects must be an iterable"):
            xr.merge([ds, 1])  # type: ignore[list-item]

    def test_merge_no_conflicts_single_var(self):
        ds1 = xr.Dataset({"a": ("x", [1, 2]), "x": [0, 1]})
        ds2 = xr.Dataset({"a": ("x", [2, 3]), "x": [1, 2]})
        expected = xr.Dataset({"a": ("x", [1, 2, 3]), "x": [0, 1, 2]})
        assert expected.identical(
            xr.merge([ds1, ds2], compat="no_conflicts", join="outer")
        )
        assert expected.identical(
            xr.merge([ds2, ds1], compat="no_conflicts", join="outer")
        )
        assert ds1.identical(xr.merge([ds1, ds2], compat="no_conflicts", join="left"))
        assert ds2.identical(xr.merge([ds1, ds2], compat="no_conflicts", join="right"))
        expected = xr.Dataset({"a": ("x", [2]), "x": [1]})
        assert expected.identical(
            xr.merge([ds1, ds2], compat="no_conflicts", join="inner")
        )

        with pytest.raises(xr.MergeError):
            ds3 = xr.Dataset({"a": ("x", [99, 3]), "x": [1, 2]})
            xr.merge([ds1, ds3], compat="no_conflicts", join="outer")

        with pytest.raises(xr.MergeError):
            ds3 = xr.Dataset({"a": ("y", [2, 3]), "y": [1, 2]})
            xr.merge([ds1, ds3], compat="no_conflicts", join="outer")

    def test_merge_no_conflicts_multi_var(self):
        data = create_test_data(add_attrs=False)
        data1 = data.copy(deep=True)
        data2 = data.copy(deep=True)

        expected = data[["var1", "var2"]]
        actual = xr.merge([data1.var1, data2.var2], compat="no_conflicts")
        assert_identical(expected, actual)

        data1["var1"][:, :5] = np.nan
        data2["var1"][:, 5:] = np.nan
        data1["var2"][:4, :] = np.nan
        data2["var2"][4:, :] = np.nan
        del data2["var3"]

        actual = xr.merge([data1, data2], compat="no_conflicts")
        assert_equal(data, actual)

    def test_merge_no_conflicts_preserve_attrs(self):
        data = xr.Dataset({"x": ([], 0, {"foo": "bar"})})
        actual = xr.merge([data, data], combine_attrs="no_conflicts")
        assert_identical(data, actual)

    def test_merge_no_conflicts_broadcast(self):
        datasets = [xr.Dataset({"x": ("y", [0])}), xr.Dataset({"x": np.nan})]
        actual = xr.merge(datasets, compat="no_conflicts")
        expected = xr.Dataset({"x": ("y", [0])})
        assert_identical(expected, actual)

        datasets = [xr.Dataset({"x": ("y", [np.nan])}), xr.Dataset({"x": 0})]
        actual = xr.merge(datasets, compat="no_conflicts")
        assert_identical(expected, actual)


class TestMergeMethod:
    def test_merge(self):
        data = create_test_data()
        ds1 = data[["var1"]]
        ds2 = data[["var3"]]
        expected = data[["var1", "var3"]]
        actual = ds1.merge(ds2)
        assert_identical(expected, actual)

        actual = ds2.merge(ds1)
        assert_identical(expected, actual)

        actual = data.merge(data)
        assert_identical(data, actual)
        actual = data.reset_coords(drop=True).merge(data)
        assert_identical(data, actual)
        actual = data.merge(data.reset_coords(drop=True))
        assert_identical(data, actual)

        with pytest.raises(ValueError, match="conflicting values for variable"):
            ds1.merge(ds2.rename({"var3": "var1"}), compat="no_conflicts")
        with pytest.raises(ValueError, match=r"should be coordinates or not"):
            data.reset_coords().merge(data)
        with pytest.raises(ValueError, match=r"should be coordinates or not"):
            data.merge(data.reset_coords())

    @pytest.mark.parametrize(
        "join", ["outer", "inner", "left", "right", "exact", "override"]
    )
    def test_merge_drop_attrs(self, join):
        data = create_test_data()
        ds1 = data[["var1"]]
        ds2 = data[["var3"]]
        ds1.coords["dim2"].attrs["keep me"] = "example"
        ds2.coords["numbers"].attrs["foo"] = "bar"
        actual = ds1.merge(ds2, combine_attrs="drop", join=join)
        assert actual.coords["dim2"].attrs == {}
        assert actual.coords["numbers"].attrs == {}
        assert ds1.coords["dim2"].attrs["keep me"] == "example"
        assert ds2.coords["numbers"].attrs["foo"] == "bar"

    def test_merge_compat_broadcast_equals(self):
        ds1 = xr.Dataset({"x": 0})
        ds2 = xr.Dataset({"x": ("y", [0, 0])})
        actual = ds1.merge(ds2, compat="broadcast_equals")
        assert_identical(ds2, actual)

        actual = ds2.merge(ds1, compat="broadcast_equals")
        assert_identical(ds2, actual)

        actual = ds1.copy()
        actual.update(ds2)
        assert_identical(ds2, actual)

        ds1 = xr.Dataset({"x": np.nan})
        ds2 = xr.Dataset({"x": ("y", [np.nan, np.nan])})
        actual = ds1.merge(ds2, compat="broadcast_equals")
        assert_identical(ds2, actual)

    def test_merge_compat(self):
        ds1 = xr.Dataset({"x": 0})
        ds2 = xr.Dataset({"x": 1})
        for compat in ["broadcast_equals", "equals", "identical", "no_conflicts"]:
            with pytest.raises(xr.MergeError):
                ds1.merge(ds2, compat=compat)  # type: ignore[arg-type]

        ds2 = xr.Dataset({"x": [0, 0]})
        for compat in ["equals", "identical"]:
            with pytest.raises(ValueError, match=r"should be coordinates or not"):
                ds1.merge(ds2, compat=compat)  # type: ignore[arg-type]

        ds2 = xr.Dataset({"x": ((), 0, {"foo": "bar"})})
        with pytest.raises(xr.MergeError):
            ds1.merge(ds2, compat="identical")

        with pytest.raises(ValueError, match=r"compat=.* invalid"):
            ds1.merge(ds2, compat="foobar")  # type: ignore[arg-type]

        assert ds1.identical(ds1.merge(ds2, compat="override"))

    def test_merge_compat_minimal(self) -> None:
        """Test that we drop the conflicting bar coordinate."""
        # https://github.com/pydata/xarray/issues/7405
        # https://github.com/pydata/xarray/issues/7588
        ds1 = xr.Dataset(coords={"foo": [1, 2, 3], "bar": 4})
        ds2 = xr.Dataset(coords={"foo": [1, 2, 3], "bar": 5})

        actual = xr.merge([ds1, ds2], compat="minimal")
        expected = xr.Dataset(coords={"foo": [1, 2, 3]})
        assert_identical(actual, expected)

    def test_merge_join_outer(self):
        ds1 = xr.Dataset({"a": ("x", [1, 2]), "x": [0, 1]})
        ds2 = xr.Dataset({"b": ("x", [3, 4]), "x": [1, 2]})
        expected = xr.Dataset(
            {"a": ("x", [1, 2, np.nan]), "b": ("x", [np.nan, 3, 4])}, {"x": [0, 1, 2]}
        )
        assert expected.identical(ds1.merge(ds2, join="outer"))
        assert expected.identical(ds2.merge(ds1, join="outer"))

        expected = expected.isel(x=slice(2))
        assert expected.identical(ds1.merge(ds2, join="left"))
        assert expected.identical(ds2.merge(ds1, join="right"))

        expected = expected.isel(x=slice(1, 2))
        assert expected.identical(ds1.merge(ds2, join="inner"))
        assert expected.identical(ds2.merge(ds1, join="inner"))

    @pytest.mark.parametrize("fill_value", [dtypes.NA, 2, 2.0, {"a": 2, "b": 1}])
    def test_merge_fill_value(self, fill_value):
        ds1 = xr.Dataset({"a": ("x", [1, 2]), "x": [0, 1]})
        ds2 = xr.Dataset({"b": ("x", [3, 4]), "x": [1, 2]})
        if fill_value == dtypes.NA:
            # if we supply the default, we expect the missing value for a
            # float array
            fill_value_a = fill_value_b = np.nan
        elif isinstance(fill_value, dict):
            fill_value_a = fill_value["a"]
            fill_value_b = fill_value["b"]
        else:
            fill_value_a = fill_value_b = fill_value

        expected = xr.Dataset(
            {"a": ("x", [1, 2, fill_value_a]), "b": ("x", [fill_value_b, 3, 4])},
            {"x": [0, 1, 2]},
        )
        assert expected.identical(ds1.merge(ds2, join="outer", fill_value=fill_value))
        assert expected.identical(ds2.merge(ds1, join="outer", fill_value=fill_value))
        assert expected.identical(
            xr.merge([ds1, ds2], join="outer", fill_value=fill_value)
        )

    def test_merge_no_conflicts(self):
        ds1 = xr.Dataset({"a": ("x", [1, 2]), "x": [0, 1]})
        ds2 = xr.Dataset({"a": ("x", [2, 3]), "x": [1, 2]})
        expected = xr.Dataset({"a": ("x", [1, 2, 3]), "x": [0, 1, 2]})

        assert expected.identical(ds1.merge(ds2, compat="no_conflicts", join="outer"))
        assert expected.identical(ds2.merge(ds1, compat="no_conflicts", join="outer"))

        assert ds1.identical(ds1.merge(ds2, compat="no_conflicts", join="left"))

        assert ds2.identical(ds1.merge(ds2, compat="no_conflicts", join="right"))

        expected2 = xr.Dataset({"a": ("x", [2]), "x": [1]})
        assert expected2.identical(ds1.merge(ds2, compat="no_conflicts", join="inner"))

        with pytest.raises(xr.MergeError):
            ds3 = xr.Dataset({"a": ("x", [99, 3]), "x": [1, 2]})
            ds1.merge(ds3, compat="no_conflicts", join="outer")

        with pytest.raises(xr.MergeError):
            ds3 = xr.Dataset({"a": ("y", [2, 3]), "y": [1, 2]})
            ds1.merge(ds3, compat="no_conflicts", join="outer")

    def test_merge_dataarray(self):
        ds = xr.Dataset({"a": 0})
        da = xr.DataArray(data=1, name="b")

        assert_identical(ds.merge(da), xr.merge([ds, da]))

    @pytest.mark.parametrize(
        ["combine_attrs", "attrs1", "attrs2", "expected_attrs", "expect_error"],
        # don't need to test thoroughly
        (
            ("drop", {"a": 0, "b": 1, "c": 2}, {"a": 1, "b": 2, "c": 3}, {}, False),
            (
                "drop_conflicts",
                {"a": 0, "b": 1, "c": 2},
                {"b": 2, "c": 2, "d": 3},
                {"a": 0, "c": 2, "d": 3},
                False,
            ),
            ("override", {"a": 0, "b": 1}, {"a": 1, "b": 2}, {"a": 0, "b": 1}, False),
            ("no_conflicts", {"a": 0, "b": 1}, {"a": 0, "b": 2}, None, True),
            ("identical", {"a": 0, "b": 1}, {"a": 0, "b": 2}, None, True),
        ),
    )
    def test_merge_combine_attrs(
        self, combine_attrs, attrs1, attrs2, expected_attrs, expect_error
    ):
        ds1 = xr.Dataset(attrs=attrs1)
        ds2 = xr.Dataset(attrs=attrs2)

        if expect_error:
            with pytest.raises(xr.MergeError):
                ds1.merge(ds2, combine_attrs=combine_attrs)
        else:
            actual = ds1.merge(ds2, combine_attrs=combine_attrs)
            expected = xr.Dataset(attrs=expected_attrs)
            assert_identical(actual, expected)


class TestNewDefaults:
    def test_merge_datasets_false_warning(self):
        data = create_test_data(add_attrs=False, use_extension_array=True)

        with set_options(use_new_combine_kwarg_defaults=False):
            old = xr.merge([data, data])

        with set_options(use_new_combine_kwarg_defaults=True):
            new = xr.merge([data, data])

        assert_identical(old, new)

    def test_merge(self):
        data = create_test_data()
        ds1 = data[["var1"]]
        ds2 = data[["var3"]]
        expected = data[["var1", "var3"]]
        with set_options(use_new_combine_kwarg_defaults=True):
            actual = ds1.merge(ds2)
            assert_identical(expected, actual)

            actual = ds2.merge(ds1)
            assert_identical(expected, actual)

            actual = data.merge(data)
            assert_identical(data, actual)

            ds1.merge(ds2.rename({"var3": "var1"}))

            with pytest.raises(ValueError, match=r"should be coordinates or not"):
                data.reset_coords().merge(data)
            with pytest.raises(ValueError, match=r"should be coordinates or not"):
                data.merge(data.reset_coords())

    def test_merge_broadcast_equals(self):
        ds1 = xr.Dataset({"x": 0})
        ds2 = xr.Dataset({"x": ("y", [0, 0])})

        with set_options(use_new_combine_kwarg_defaults=False):
            with pytest.warns(
                FutureWarning,
                match="will change from compat='no_conflicts' to compat='override'",
            ):
                old = ds1.merge(ds2)

        with set_options(use_new_combine_kwarg_defaults=True):
            new = ds1.merge(ds2)

        assert_identical(ds2, old)
        with pytest.raises(AssertionError):
            assert_identical(old, new)

    def test_merge_auto_align(self):
        ds1 = xr.Dataset({"a": ("x", [1, 2]), "x": [0, 1]})
        ds2 = xr.Dataset({"b": ("x", [3, 4]), "x": [1, 2]})
        expected = xr.Dataset(
            {"a": ("x", [1, 2, np.nan]), "b": ("x", [np.nan, 3, 4])}, {"x": [0, 1, 2]}
        )
        with set_options(use_new_combine_kwarg_defaults=False):
            with pytest.warns(
                FutureWarning, match="will change from join='outer' to join='exact'"
            ):
                assert expected.identical(ds1.merge(ds2))
            with pytest.warns(
                FutureWarning, match="will change from join='outer' to join='exact'"
            ):
                assert expected.identical(ds2.merge(ds1))

        with set_options(use_new_combine_kwarg_defaults=True):
            with pytest.raises(ValueError, match="might be related to new default"):
                expected.identical(ds2.merge(ds1))


class TestMergeDataTree:
    def test_mixed(self) -> None:
        tree = xr.DataTree()
        ds = xr.Dataset()
        with pytest.raises(
            TypeError,
            match="merge does not support mixed type arguments when one argument is a DataTree",
        ):
            xr.merge([tree, ds])  # type: ignore[list-item]

    def test_distinct(self) -> None:
        tree1 = xr.DataTree.from_dict({"/a/b/c": 1})
        tree2 = xr.DataTree.from_dict({"/a/d/e": 2})
        expected = xr.DataTree.from_dict({"/a/b/c": 1, "/a/d/e": 2})
        merged = xr.merge([tree1, tree2])
        assert_equal(merged, expected)

    def test_overlap(self) -> None:
        tree1 = xr.DataTree.from_dict({"/a/b": 1})
        tree2 = xr.DataTree.from_dict({"/a/c": 2})
        tree3 = xr.DataTree.from_dict({"/a/d": 3})
        expected = xr.DataTree.from_dict({"/a/b": 1, "/a/c": 2, "/a/d": 3})
        merged = xr.merge([tree1, tree2, tree3])
        assert_equal(merged, expected)

    def test_inherited(self) -> None:
        tree1 = xr.DataTree.from_dict({"/a/b": ("x", [1])}, coords={"x": [0]})
        tree2 = xr.DataTree.from_dict({"/a/c": ("x", [2])})
        expected = xr.DataTree.from_dict(
            {"/a/b": ("x", [1]), "a/c": ("x", [2])}, coords={"x": [0]}
        )
        merged = xr.merge([tree1, tree2])
        assert_equal(merged, expected)

    def test_inherited_join(self) -> None:
        tree1 = xr.DataTree.from_dict({"/a/b": ("x", [0, 1])}, coords={"x": [0, 1]})
        tree2 = xr.DataTree.from_dict({"/a/c": ("x", [1, 2])}, coords={"x": [1, 2]})

        expected = xr.DataTree.from_dict(
            {"/a/b": ("x", [0, 1]), "a/c": ("x", [np.nan, 1])}, coords={"x": [0, 1]}
        )
        merged = xr.merge([tree1, tree2], join="left")
        assert_equal(merged, expected)

        expected = xr.DataTree.from_dict(
            {"/a/b": ("x", [1, np.nan]), "a/c": ("x", [1, 2])}, coords={"x": [1, 2]}
        )
        merged = xr.merge([tree1, tree2], join="right")
        assert_equal(merged, expected)

        expected = xr.DataTree.from_dict(
            {"/a/b": ("x", [1]), "a/c": ("x", [1])}, coords={"x": [1]}
        )
        merged = xr.merge([tree1, tree2], join="inner")
        assert_equal(merged, expected)

        expected = xr.DataTree.from_dict(
            {"/a/b": ("x", [0, 1, np.nan]), "a/c": ("x", [np.nan, 1, 2])},
            coords={"x": [0, 1, 2]},
        )
        merged = xr.merge([tree1, tree2], join="outer")
        assert_equal(merged, expected)

        with pytest.raises(
            xr.AlignmentError,
            match=re.escape("cannot align objects with join='exact'"),
        ):
            xr.merge([tree1, tree2], join="exact")

    def test_merge_error_includes_path(self) -> None:
        tree1 = xr.DataTree.from_dict({"/a/b": ("x", [0, 1])})
        tree2 = xr.DataTree.from_dict({"/a/b": ("x", [1, 2])})
        with pytest.raises(
            xr.MergeError,
            match=re.escape(
                "Raised whilst mapping function over node(s) with path 'a'"
            ),
        ):
            xr.merge([tree1, tree2], join="exact", compat="no_conflicts")

    def test_fill_value_errors(self) -> None:
        trees = [xr.DataTree(), xr.DataTree()]

        with pytest.raises(
            NotImplementedError,
            match=re.escape(
                "fill_value is not yet supported for DataTree objects in merge"
            ),
        ):
            xr.merge(trees, fill_value=None)
