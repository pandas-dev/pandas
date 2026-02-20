"""Tests for Dataset.eval() functionality."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

import xarray as xr
from xarray import DataArray, Dataset
from xarray.tests import (
    assert_equal,
    assert_identical,
    raise_if_dask_computes,
    requires_dask,
)


def test_eval(ds) -> None:
    """Test basic eval functionality."""
    actual = ds.eval("z1 + 5")
    expect = ds["z1"] + 5
    assert_identical(expect, actual)

    # Use bitwise operators for element-wise operations on arrays
    actual = ds.eval("(z1 > 5) & (z2 > 0)")
    expect = (ds["z1"] > 5) & (ds["z2"] > 0)
    assert_identical(expect, actual)


def test_eval_parser_deprecated(ds) -> None:
    """Test that passing parser= raises a FutureWarning."""
    with pytest.warns(FutureWarning, match="parser.*deprecated"):
        ds.eval("z1 + 5", parser="pandas")


def test_eval_logical_operators(ds) -> None:
    """Test that 'and'/'or'/'not' are transformed for query() consistency.

    These operators are transformed to '&'/'|'/'~' to match pd.eval() behavior,
    which query() uses. This ensures syntax that works in query() also works in
    eval().
    """
    # 'and' transformed to '&'
    actual = ds.eval("(z1 > 5) and (z2 > 0)")
    expect = (ds["z1"] > 5) & (ds["z2"] > 0)
    assert_identical(expect, actual)

    # 'or' transformed to '|'
    actual = ds.eval("(z1 > 5) or (z2 > 0)")
    expect = (ds["z1"] > 5) | (ds["z2"] > 0)
    assert_identical(expect, actual)

    # 'not' transformed to '~'
    actual = ds.eval("not (z1 > 5)")
    expect = ~(ds["z1"] > 5)
    assert_identical(expect, actual)


def test_eval_ndimensional() -> None:
    """Test that eval works with N-dimensional data where N > 2."""
    # Create a 3D dataset - this previously failed with pd.eval
    rng = np.random.default_rng(42)
    ds = Dataset(
        {
            "x": (["time", "lat", "lon"], rng.random((3, 4, 5))),
            "y": (["time", "lat", "lon"], rng.random((3, 4, 5))),
        }
    )

    # Basic arithmetic
    actual = ds.eval("x + y")
    expect = ds["x"] + ds["y"]
    assert_identical(expect, actual)

    # Assignment
    actual = ds.eval("z = x + y")
    assert "z" in actual.data_vars
    assert_equal(ds["x"] + ds["y"], actual["z"])

    # Complex expression
    actual = ds.eval("x * 2 + y ** 2")
    expect = ds["x"] * 2 + ds["y"] ** 2
    assert_identical(expect, actual)

    # Comparison
    actual = ds.eval("x > y")
    expect = ds["x"] > ds["y"]
    assert_identical(expect, actual)

    # Use bitwise operators for element-wise boolean operations
    actual = ds.eval("(x > 0.5) & (y < 0.5)")
    expect = (ds["x"] > 0.5) & (ds["y"] < 0.5)
    assert_identical(expect, actual)


def test_eval_chained_comparisons() -> None:
    """Test that chained comparisons are transformed for query() consistency.

    Chained comparisons like 'a < b < c' are transformed to '(a < b) & (b < c)'
    to match pd.eval() behavior, which query() uses.
    """
    ds = Dataset({"x": ("dim", np.arange(10))})

    # Basic chained comparison: 2 < x < 7
    actual = ds.eval("2 < x < 7")
    expect = (ds["x"] > 2) & (ds["x"] < 7)
    assert_identical(expect, actual)

    # Mixed operators: 0 <= x < 5
    actual = ds.eval("0 <= x < 5")
    expect = (ds["x"] >= 0) & (ds["x"] < 5)
    assert_identical(expect, actual)

    # Explicit bitwise operators also work
    actual = ds.eval("(x > 2) & (x < 7)")
    expect = (ds["x"] > 2) & (ds["x"] < 7)
    assert_identical(expect, actual)


def test_eval_restricted_syntax() -> None:
    """Test that eval blocks certain syntax to emulate pd.eval() behavior."""
    ds = Dataset({"a": ("x", [1, 2, 3])})

    # Private attribute access is not allowed (consistent with pd.eval)
    with pytest.raises(ValueError, match="Access to private attributes is not allowed"):
        ds.eval("a.__class__")

    with pytest.raises(ValueError, match="Access to private attributes is not allowed"):
        ds.eval("a._private")

    # Lambda expressions are not allowed (pd.eval: "Only named functions are supported")
    with pytest.raises(ValueError, match="Lambda expressions are not allowed"):
        ds.eval("(lambda x: x + 1)(a)")

    # These builtins are not in the namespace
    with pytest.raises(NameError):
        ds.eval("__import__('os')")

    with pytest.raises(NameError):
        ds.eval("open('file.txt')")


def test_eval_unsupported_statements() -> None:
    """Test that unsupported statement types produce clear errors."""
    ds = Dataset({"a": ("x", [1, 2, 3])})

    # Augmented assignment is not supported
    with pytest.raises(ValueError, match="Unsupported statement type"):
        ds.eval("a += 1")


def test_eval_functions() -> None:
    """Test that numpy and other functions work in eval."""
    ds = Dataset({"a": ("x", [0.0, 1.0, 4.0])})

    # numpy functions via np namespace should work
    result = ds.eval("np.sqrt(a)")
    assert_equal(result, np.sqrt(ds["a"]))

    result = ds.eval("np.sin(a) + np.cos(a)")
    assert_equal(result, np.sin(ds["a"]) + np.cos(ds["a"]))

    # pandas namespace should work
    result = ds.eval("pd.isna(a)")
    # pd.isna returns ndarray, not DataArray
    np.testing.assert_array_equal(result, pd.isna(ds["a"].values))

    # xarray namespace should work
    result = ds.eval("xr.where(a > 1, a, 0)")

    assert_equal(result, xr.where(ds["a"] > 1, ds["a"], 0))

    # Common builtins should work
    result = ds.eval("abs(a - 2)")
    assert_equal(result, abs(ds["a"] - 2))

    result = ds.eval("round(float(a.mean()))")
    assert result == round(float(ds["a"].mean()))

    result = ds.eval("len(a)")
    assert result == 3

    result = ds.eval("pow(a, 2)")
    assert_equal(result, ds["a"] ** 2)

    # Attribute access on DataArrays should work
    result = ds.eval("a.values")
    assert isinstance(result, np.ndarray)

    # Method calls on DataArrays should work
    result = ds.eval("a.mean()")
    assert float(result) == np.mean([0.0, 1.0, 4.0])


def test_eval_extended_builtins() -> None:
    """Test extended builtins available in eval namespace.

    These builtins are safe (no I/O, no code execution) and commonly needed
    for typical xarray operations like slicing, type conversion, and iteration.
    """
    ds = Dataset(
        {"a": ("x", [1.0, 2.0, 3.0, 4.0, 5.0])},
        coords={"time": pd.date_range("2019-01-01", periods=5)},
    )

    # slice - essential for .sel() with ranges
    result = ds.eval("a.sel(x=slice(1, 3))")
    expected = ds["a"].sel(x=slice(1, 3))
    assert_equal(result, expected)

    # str - type constructor
    result = ds.eval("str(int(a.mean()))")
    assert result == "3"

    # list, tuple - type constructors
    result = ds.eval("list(range(3))")
    assert result == [0, 1, 2]

    result = ds.eval("tuple(range(3))")
    assert result == (0, 1, 2)

    # dict, set - type constructors
    result = ds.eval("dict(x=1, y=2)")
    assert result == {"x": 1, "y": 2}

    result = ds.eval("set([1, 2, 2, 3])")
    assert result == {1, 2, 3}

    # range - iteration
    result = ds.eval("list(range(3))")
    assert result == [0, 1, 2]

    # zip, enumerate - iteration helpers
    result = ds.eval("list(zip([1, 2], [3, 4]))")
    assert result == [(1, 3), (2, 4)]

    result = ds.eval("list(enumerate(['a', 'b']))")
    assert result == [(0, "a"), (1, "b")]

    # map, filter - functional helpers
    result = ds.eval("list(map(abs, [-1, -2, 3]))")
    assert result == [1, 2, 3]

    result = ds.eval("list(filter(bool, [0, 1, 0, 2]))")
    assert result == [1, 2]

    # any, all - aggregation
    result = ds.eval("any([False, True, False])")
    assert result is True

    result = ds.eval("all([True, True, True])")
    assert result is True

    result = ds.eval("all([True, False, True])")
    assert result is False


def test_eval_data_variable_priority() -> None:
    """Test that data variables take priority over builtin functions.

    Users may have data variables named 'sum', 'abs', 'min', etc. When they
    reference these in eval(), they should get their data, not the Python builtins.
    The builtins should still be accessible via the np namespace (np.sum, np.abs).
    """
    # Create dataset with data variables that shadow builtins
    ds = Dataset(
        {
            "sum": ("x", [10.0, 20.0, 30.0]),  # shadows builtin sum
            "abs": ("x", [1.0, 2.0, 3.0]),  # shadows builtin abs
            "min": ("x", [100.0, 200.0, 300.0]),  # shadows builtin min
            "other": ("x", [5.0, 10.0, 15.0]),
        }
    )

    # Data variables should take priority - user data wins
    result = ds.eval("sum + other")
    expected = ds["sum"] + ds["other"]
    assert_equal(result, expected)

    # Should get the data variable, not builtin sum applied to something
    result = ds.eval("sum * 2")
    expected = ds["sum"] * 2
    assert_equal(result, expected)

    # abs as data variable should work
    result = ds.eval("abs + 1")
    expected = ds["abs"] + 1
    assert_equal(result, expected)

    # min as data variable should work
    result = ds.eval("min - 50")
    expected = ds["min"] - 50
    assert_equal(result, expected)

    # np namespace should still provide access to actual functions
    result = ds.eval("np.abs(other - 10)")
    expected = abs(ds["other"] - 10)
    assert_equal(result, expected)

    # np.sum should work even when 'sum' is a data variable
    result = ds.eval("np.sum(other)")
    expected = np.sum(ds["other"])
    assert result == expected


def test_eval_coordinate_priority() -> None:
    """Test that coordinates also take priority over builtins."""
    ds = Dataset(
        {"data": ("x", [1.0, 2.0, 3.0])},
        coords={"sum": ("x", [10.0, 20.0, 30.0])},  # coordinate named 'sum'
    )

    # Coordinate should be accessible and take priority over builtin
    result = ds.eval("data + sum")
    expected = ds["data"] + ds.coords["sum"]
    assert_equal(result, expected)


# Error message tests


def test_eval_error_undefined_variable() -> None:
    """Test error message when referencing an undefined variable."""
    ds = Dataset({"a": ("x", [1, 2, 3])})
    with pytest.raises(NameError, match="undefined_var"):
        ds.eval("undefined_var + a")


def test_eval_error_syntax() -> None:
    """Test error message for malformed expressions."""
    ds = Dataset({"a": ("x", [1, 2, 3])})
    with pytest.raises(ValueError, match="Invalid"):
        ds.eval("a +")


def test_eval_error_invalid_assignment() -> None:
    """Test error message when assignment target is invalid."""
    ds = Dataset({"a": ("x", [1, 2, 3])})
    # "1 = a" should fail during parsing - can't assign to a literal
    with pytest.raises(ValueError, match="Invalid"):
        ds.eval("1 = a")


def test_eval_error_dunder_access() -> None:
    """Test error message when trying to access dunder attributes."""
    ds = Dataset({"a": ("x", [1, 2, 3])})
    with pytest.raises(ValueError, match="private attributes"):
        ds.eval("a.__class__")


def test_eval_error_missing_method() -> None:
    """Test error message when calling a nonexistent method."""
    ds = Dataset({"a": ("x", [1, 2, 3])})
    # This should raise AttributeError from the DataArray
    with pytest.raises(AttributeError, match="nonexistent_method"):
        ds.eval("a.nonexistent_method()")


def test_eval_error_type_mismatch() -> None:
    """Test error message when types are incompatible."""
    ds = Dataset({"a": ("x", [1, 2, 3])})
    # Adding string to numeric array should raise TypeError or similar
    with pytest.raises((TypeError, np.exceptions.DTypePromotionError)):
        ds.eval("a + 'string'")


# Edge case tests


def test_eval_empty_expression() -> None:
    """Test handling of empty expression string."""
    ds = Dataset({"a": ("x", [1, 2, 3])})
    with pytest.raises(ValueError):
        ds.eval("")


def test_eval_whitespace_only_expression() -> None:
    """Test handling of whitespace-only expression."""
    ds = Dataset({"a": ("x", [1, 2, 3])})
    with pytest.raises(ValueError):
        ds.eval("   ")


def test_eval_just_variable_name() -> None:
    """Test that just a variable name returns the variable."""
    ds = Dataset({"a": ("x", [1, 2, 3])})
    result = ds.eval("a")
    expected = ds["a"]
    assert_equal(result, expected)


def test_eval_unicode_variable_names() -> None:
    """Test that unicode variable names work in expressions."""
    # Greek letters are valid Python identifiers
    ds = Dataset({"α": ("x", [1.0, 2.0, 3.0]), "β": ("x", [4.0, 5.0, 6.0])})
    result = ds.eval("α + β")
    expected = ds["α"] + ds["β"]
    assert_equal(result, expected)


def test_eval_long_expression() -> None:
    """Test that very long expressions work correctly."""
    ds = Dataset({"a": ("x", [1.0, 2.0, 3.0])})
    # Build a long expression: a + a + a + ... (50 times)
    long_expr = " + ".join(["a"] * 50)
    result = ds.eval(long_expr)
    expected = ds["a"] * 50
    assert_equal(result, expected)


# Dask tests


@requires_dask
def test_eval_dask_basic_arithmetic() -> None:
    """Test that basic arithmetic with dask arrays returns dask-backed result."""
    from xarray.core.utils import is_duck_dask_array

    ds = Dataset(
        {"a": ("x", np.arange(10.0)), "b": ("x", np.linspace(0, 1, 10))}
    ).chunk({"x": 5})

    with raise_if_dask_computes():
        result = ds.eval("a + b")

    assert isinstance(result, DataArray)
    assert is_duck_dask_array(result.data)

    # Verify correctness when computed
    expected = ds["a"] + ds["b"]
    assert_equal(result, expected)


@requires_dask
def test_eval_dask_assignment() -> None:
    """Test that assignments with dask arrays preserve lazy evaluation."""
    from xarray.core.utils import is_duck_dask_array

    ds = Dataset(
        {"a": ("x", np.arange(10.0)), "b": ("x", np.linspace(0, 1, 10))}
    ).chunk({"x": 5})

    with raise_if_dask_computes():
        result = ds.eval("z = a + b")

    assert isinstance(result, Dataset)
    assert "z" in result.data_vars
    assert is_duck_dask_array(result["z"].data)

    # Verify correctness when computed
    expected = ds["a"] + ds["b"]
    assert_equal(result["z"], expected)


@requires_dask
def test_eval_dask_method_chaining() -> None:
    """Test that method chaining works with dask arrays."""
    ds = Dataset({"a": (("x", "y"), np.arange(20.0).reshape(4, 5))}).chunk(
        {"x": 2, "y": 5}
    )

    # Calling .mean() should still be lazy
    result = ds.eval("a.mean(dim='x')")
    # Calling .compute() should return numpy-backed result
    computed = result.compute()

    expected = ds["a"].mean(dim="x").compute()
    assert_equal(computed, expected)


@requires_dask
def test_eval_dask_xr_where() -> None:
    """Test that xr.where() with dask arrays preserves lazy evaluation."""
    from xarray.core.utils import is_duck_dask_array

    ds = Dataset({"a": ("x", np.arange(-5, 5, dtype=float))}).chunk({"x": 5})

    with raise_if_dask_computes():
        result = ds.eval("xr.where(a > 0, a, 0)")

    assert isinstance(result, DataArray)
    assert is_duck_dask_array(result.data)

    # Verify correctness when computed
    expected = xr.where(ds["a"] > 0, ds["a"], 0)
    assert_equal(result, expected)


@requires_dask
def test_eval_dask_complex_expression() -> None:
    """Test that complex expressions preserve dask backing."""
    from xarray.core.utils import is_duck_dask_array

    rng = np.random.default_rng(42)
    ds = Dataset(
        {
            "x": (["time", "lat", "lon"], rng.random((3, 4, 5))),
            "y": (["time", "lat", "lon"], rng.random((3, 4, 5))),
        }
    ).chunk({"time": 1, "lat": 2, "lon": 5})

    with raise_if_dask_computes():
        result = ds.eval("x * 2 + y ** 2")

    assert is_duck_dask_array(result.data)

    # Verify correctness when computed
    expected = ds["x"] * 2 + ds["y"] ** 2
    assert_equal(result, expected)


@requires_dask
def test_eval_dask_mixed_backends() -> None:
    """Test expressions with mixed dask and numpy arrays."""
    from xarray.core.utils import is_duck_dask_array

    ds = Dataset(
        {
            "dask_var": ("x", np.arange(10.0)),
            "numpy_var": ("x", np.linspace(0, 1, 10)),
        }
    )
    # Only chunk one variable
    ds["dask_var"] = ds["dask_var"].chunk({"x": 5})

    with raise_if_dask_computes():
        result = ds.eval("dask_var + numpy_var")

    # Result should be dask-backed when any input is dask
    assert is_duck_dask_array(result.data)

    # Verify correctness
    expected = ds["dask_var"] + ds["numpy_var"]
    assert_equal(result, expected)


@requires_dask
def test_eval_dask_np_functions() -> None:
    """Test that numpy functions via np namespace preserve dask."""
    from xarray.core.utils import is_duck_dask_array

    ds = Dataset({"a": ("x", np.arange(1.0, 11.0))}).chunk({"x": 5})

    with raise_if_dask_computes():
        result = ds.eval("np.sqrt(a)")

    assert is_duck_dask_array(result.data)

    # Verify correctness
    expected = np.sqrt(ds["a"])
    assert_equal(result, expected)


@requires_dask
def test_eval_dask_comparison() -> None:
    """Test that comparison operations preserve dask backing."""
    from xarray.core.utils import is_duck_dask_array

    ds = Dataset(
        {"a": ("x", np.arange(10.0)), "b": ("x", np.arange(10.0)[::-1])}
    ).chunk({"x": 5})

    with raise_if_dask_computes():
        result = ds.eval("a > b")

    assert is_duck_dask_array(result.data)

    # Verify correctness
    expected = ds["a"] > ds["b"]
    assert_equal(result, expected)


@requires_dask
def test_eval_dask_boolean_operators() -> None:
    """Test that bitwise boolean operators preserve dask."""
    from xarray.core.utils import is_duck_dask_array

    ds = Dataset(
        {"a": ("x", np.arange(10.0)), "b": ("x", np.arange(10.0)[::-1])}
    ).chunk({"x": 5})

    with raise_if_dask_computes():
        result = ds.eval("(a > 3) & (b < 7)")

    assert is_duck_dask_array(result.data)

    # Verify correctness
    expected = (ds["a"] > 3) & (ds["b"] < 7)
    assert_equal(result, expected)


@requires_dask
def test_eval_dask_chained_comparisons() -> None:
    """Test that chained comparisons preserve dask backing."""
    from xarray.core.utils import is_duck_dask_array

    ds = Dataset({"x": ("dim", np.arange(10.0))}).chunk({"dim": 5})

    with raise_if_dask_computes():
        result = ds.eval("2 < x < 7")

    assert is_duck_dask_array(result.data)

    # Verify correctness
    expected = (ds["x"] > 2) & (ds["x"] < 7)
    assert_equal(result, expected)
