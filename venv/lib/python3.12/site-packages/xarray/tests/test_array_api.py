from __future__ import annotations

import pytest

import xarray as xr
from xarray.testing import assert_equal

np = pytest.importorskip("numpy", minversion="1.22")
xp = pytest.importorskip("array_api_strict")

from array_api_strict._array_object import Array  # isort:skip # type: ignore[no-redef]


@pytest.fixture
def arrays() -> tuple[xr.DataArray, xr.DataArray]:
    np_arr = xr.DataArray(
        np.array([[1.0, 2.0, 3.0], [4.0, 5.0, np.nan]]),
        dims=("x", "y"),
        coords={"x": [10, 20]},
    )
    xp_arr = xr.DataArray(
        xp.asarray([[1.0, 2.0, 3.0], [4.0, 5.0, np.nan]]),
        dims=("x", "y"),
        coords={"x": [10, 20]},
    )
    assert isinstance(xp_arr.data, Array)
    return np_arr, xp_arr


def test_arithmetic(arrays: tuple[xr.DataArray, xr.DataArray]) -> None:
    np_arr, xp_arr = arrays
    expected = np_arr + 7
    actual = xp_arr + 7
    assert isinstance(actual.data, Array)
    assert_equal(actual, expected)


def test_aggregation(arrays: tuple[xr.DataArray, xr.DataArray]) -> None:
    np_arr, xp_arr = arrays
    expected = np_arr.sum()
    actual = xp_arr.sum()
    assert isinstance(actual.data, Array)
    assert_equal(actual, expected)


def test_aggregation_skipna(arrays) -> None:
    np_arr, xp_arr = arrays
    expected = np_arr.sum(skipna=False)
    actual = xp_arr.sum(skipna=False)
    assert isinstance(actual.data, Array)
    assert_equal(actual, expected)


# casting nan warns
@pytest.mark.filterwarnings("ignore:invalid value encountered in cast")
def test_astype(arrays) -> None:
    np_arr, xp_arr = arrays
    expected = np_arr.astype(np.int64)
    actual = xp_arr.astype(xp.int64)
    assert actual.dtype == xp.int64
    assert isinstance(actual.data, Array)
    assert_equal(actual, expected)


def test_broadcast(arrays: tuple[xr.DataArray, xr.DataArray]) -> None:
    np_arr, xp_arr = arrays
    np_arr2 = xr.DataArray(np.array([1.0, 2.0]), dims="x")
    xp_arr2 = xr.DataArray(xp.asarray([1.0, 2.0]), dims="x")

    expected = xr.broadcast(np_arr, np_arr2)
    actual = xr.broadcast(xp_arr, xp_arr2)
    assert len(actual) == len(expected)
    for a, e in zip(actual, expected, strict=True):
        assert isinstance(a.data, Array)
        assert_equal(a, e)


def test_broadcast_during_arithmetic(arrays: tuple[xr.DataArray, xr.DataArray]) -> None:
    np_arr, xp_arr = arrays
    np_arr2 = xr.DataArray(np.array([1.0, 2.0]), dims="x")
    xp_arr2 = xr.DataArray(xp.asarray([1.0, 2.0]), dims="x")

    expected = np_arr * np_arr2
    actual = xp_arr * xp_arr2
    assert isinstance(actual.data, Array)
    assert_equal(actual, expected)

    expected = np_arr2 * np_arr
    actual = xp_arr2 * xp_arr
    assert isinstance(actual.data, Array)
    assert_equal(actual, expected)


def test_concat(arrays: tuple[xr.DataArray, xr.DataArray]) -> None:
    np_arr, xp_arr = arrays
    expected = xr.concat((np_arr, np_arr), dim="x")
    actual = xr.concat((xp_arr, xp_arr), dim="x")
    assert isinstance(actual.data, Array)
    assert_equal(actual, expected)


def test_indexing(arrays: tuple[xr.DataArray, xr.DataArray]) -> None:
    np_arr, xp_arr = arrays
    expected = np_arr[:, 0]
    actual = xp_arr[:, 0]
    assert isinstance(actual.data, Array)
    assert_equal(actual, expected)


def test_properties(arrays: tuple[xr.DataArray, xr.DataArray]) -> None:
    np_arr, xp_arr = arrays

    expected = np_arr.data.nbytes
    assert np_arr.nbytes == expected
    assert xp_arr.nbytes == expected


def test_reorganizing_operation(arrays: tuple[xr.DataArray, xr.DataArray]) -> None:
    np_arr, xp_arr = arrays
    expected = np_arr.transpose()
    actual = xp_arr.transpose()
    assert isinstance(actual.data, Array)
    assert_equal(actual, expected)


def test_stack(arrays: tuple[xr.DataArray, xr.DataArray]) -> None:
    np_arr, xp_arr = arrays
    expected = np_arr.stack(z=("x", "y"))
    actual = xp_arr.stack(z=("x", "y"))
    assert isinstance(actual.data, Array)
    assert_equal(actual, expected)


def test_unstack(arrays: tuple[xr.DataArray, xr.DataArray]) -> None:
    np_arr, xp_arr = arrays
    expected = np_arr.stack(z=("x", "y")).unstack()
    actual = xp_arr.stack(z=("x", "y")).unstack()
    assert isinstance(actual.data, Array)
    assert_equal(actual, expected)


def test_where() -> None:
    np_arr = xr.DataArray(np.array([1, 0]), dims="x")
    xp_arr = xr.DataArray(xp.asarray([1, 0]), dims="x")
    expected = xr.where(np_arr, 1, 0)
    actual = xr.where(xp_arr, 1, 0)
    assert isinstance(actual.data, Array)
    assert_equal(actual, expected)
