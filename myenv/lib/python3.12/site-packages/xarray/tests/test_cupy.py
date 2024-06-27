from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

import xarray as xr

cp = pytest.importorskip("cupy")


@pytest.fixture
def toy_weather_data():
    """Construct the example DataSet from the Toy weather data example.

    https://docs.xarray.dev/en/stable/examples/weather-data.html

    Here we construct the DataSet exactly as shown in the example and then
    convert the numpy arrays to cupy.

    """
    np.random.seed(123)
    times = pd.date_range("2000-01-01", "2001-12-31", name="time")
    annual_cycle = np.sin(2 * np.pi * (times.dayofyear.values / 365.25 - 0.28))

    base = 10 + 15 * annual_cycle.reshape(-1, 1)
    tmin_values = base + 3 * np.random.randn(annual_cycle.size, 3)
    tmax_values = base + 10 + 3 * np.random.randn(annual_cycle.size, 3)

    ds = xr.Dataset(
        {
            "tmin": (("time", "location"), tmin_values),
            "tmax": (("time", "location"), tmax_values),
        },
        {"time": times, "location": ["IA", "IN", "IL"]},
    )

    ds.tmax.data = cp.asarray(ds.tmax.data)
    ds.tmin.data = cp.asarray(ds.tmin.data)

    return ds


def test_cupy_import() -> None:
    """Check the import worked."""
    assert cp


def test_check_data_stays_on_gpu(toy_weather_data) -> None:
    """Perform some operations and check the data stays on the GPU."""
    freeze = (toy_weather_data["tmin"] <= 0).groupby("time.month").mean("time")
    assert isinstance(freeze.data, cp.ndarray)


def test_where() -> None:
    from xarray.core.duck_array_ops import where

    data = cp.zeros(10)

    output = where(data < 1, 1, data).all()
    assert output
    assert isinstance(output, cp.ndarray)
