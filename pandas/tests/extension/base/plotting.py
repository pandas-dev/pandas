"""
Tests for plotting compatibility.
"""

from typing import Literal

import pytest

mpl = pytest.importorskip("matplotlib", reason="test requires matplotlib")
from matplotlib import (
    axes as maxes,
    axis as maxis,
    lines as mlines,
    units as munits,
)
import numpy as np

import pandas as pd
import pandas._testing as tm
from pandas.api.types import (
    PeriodDtype,
    is_any_real_numeric_dtype,
    is_bool_dtype,
    is_datetime64_any_dtype,
)
from pandas.core.arrays import ExtensionArray
from pandas.tests.plotting.common import _check_plot_works

from pandas.plotting._matplotlib.converter import pandas_converters


def _get_plot_df(data: ExtensionArray) -> pd.DataFrame:
    """Helper function to create a DataFrame for plotting tests.

    Parameters
    ----------
    data : ExtensionArray
        The ExtensionArray to be used as the 'Data' column in the DataFrame.

    Returns
    -------
    pd.DataFrame
        A DataFrame with two columns: 'Data' containing the provided ExtensionArray and
        'Numeric' containing a range of integers.
    """
    return pd.DataFrame({"Data": data, "Numeric": np.arange(len(data))}).dropna()


def _check_plot_data(
    ax: maxes.Axes,
    ser: pd.Series,
    axis: Literal["x", "y"],
) -> None:
    """Check that the data plotted matches the expected converted data from the Series.

    Parameters
    ----------
    ax : maxes.Axes
        The Axes object containing the plot.
    ser : pd.Series
        The Series containing the expected data.
    axis : Literal["x", "y"]
        The axis to check ("x" or "y").
    """
    arr = ser.to_numpy()
    if munits._is_natively_supported(arr) or is_bool_dtype(arr):
        # Convert natively or boolean just to float
        converted_data = arr.astype(np.float64)
    else:
        # Need to get registered converter for non-natively
        type_ = ser.dtype.type
        with pandas_converters():
            converter = munits.registry[type_]
        axis_: maxis.Axis = getattr(ax, f"{axis}axis")
        unit = converter.default_units(arr, axis_)
        converted_data = np.array(converter.convert(arr, unit, axis_)).astype(
            np.float64
        )

    # Get the data plotted on the specified axis from the first line in the Axes
    line: mlines.Line2D = ax.lines[0]
    line_data = getattr(line, f"get_{axis}data")(orig=False)

    # Assert that the plotted data matches the expected converted data
    tm.assert_almost_equal(line_data, converted_data)


def _plot(
    data: ExtensionArray,
    x: str,
    y: str,
    **kwargs,
) -> None:
    """Helper function to plot a DataFrame with the given x and y columns.

    Parameters
    ----------
    data : ExtensionArray
        The ExtensionArray to plot.
    x : str
        The name of the column to use for the x-axis.
    y : str
        The name of the column to use for the y-axis.
    kind : str, optional
        The kind of plot to create. If None, the default plot type is used.
    **kwargs
        Additional keyword arguments to pass to the plot function.
    """
    # Create a DataFrame for plotting using the provided ExtensionArray data
    plot_df = _get_plot_df(data)

    # Check that the plot works with the specified x and y columns and any additional
    # keyword arguments
    ax = _check_plot_works(plot_df.plot, x=x, y=y, **kwargs)

    # Check that the data plotted on the specified axes matches the expected data from
    # the DataFrame
    x_dtype = plot_df[x].dtype
    for axis, col in zip(["x", "y"], [x, y], strict=True):
        # Skip checking the x-axis if the data type is not numeric, period, or datetime
        # as pandas uses just numeric values for all other types and use the actual
        # values only as labels.
        if axis == "x" and not (
            is_any_real_numeric_dtype(x_dtype)
            or isinstance(x_dtype, PeriodDtype)
            or is_datetime64_any_dtype(x_dtype)
        ):
            continue
        _check_plot_data(ax, plot_df[col], axis)  # type: ignore[arg-type]


class BasePlottingTests:
    # Note: these are ONLY for ExtensionArray subclasses that support plotting.

    def test_plot_on_x_axis(self, data):
        """Test that EA data can be plotted on the x-axis."""
        _plot(data, x="Data", y="Numeric")

    def test_plot_on_y_axis(self, data, **kwargs):
        """Test that EA data can be plotted on the y-axis."""
        _plot(data, x="Numeric", y="Data", **kwargs)
