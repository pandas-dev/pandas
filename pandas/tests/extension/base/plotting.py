"""
Tests for plotting compatibility.
"""

from __future__ import annotations

from typing import (
    TYPE_CHECKING,
    Literal,
)
import warnings

if TYPE_CHECKING:
    from matplotlib import (
        axes as maxes,
        axis as maxis,
        lines as mlines,
    )
    from pandas.core.arrays import ExtensionArray

import numpy as np
import pytest

import pandas as pd
import pandas._testing as tm
from pandas.api.types import (
    PeriodDtype,
    is_any_real_numeric_dtype,
    is_bool_dtype,
    is_datetime64_any_dtype,
)
from pandas.tests.plotting.common import _check_plot_works
from pandas.util.version import Version


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
    from matplotlib import units as munits

    from pandas.plotting._matplotlib.converter import pandas_converters

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
    plot_data: pd.DataFrame,
    x: str,
    y: str,
    **kwargs,
) -> None:
    """Helper function to plot a DataFrame with the given x and y columns.

    Parameters
    ----------
    data : DataFrame
        A DataFrame with an ExtensionArray in a "Data" column.
    x : str
        The name of the column to use for the x-axis.
    y : str
        The name of the column to use for the y-axis.
    kind : str, optional
        The kind of plot to create. If None, the default plot type is used.
    **kwargs
        Additional keyword arguments to pass to the plot function.
    """
    # Set include_bool flag for boolean dtypes
    if is_bool_dtype(plot_data["Data"].dtype) and "include_bool" not in kwargs:
        kwargs["include_bool"] = True

    # Check that the plot works with the specified x and y columns and any additional
    # keyword arguments
    ax = _check_plot_works(plot_data.plot, x=x, y=y, **kwargs)

    # Check that the data plotted on the specified axes matches the expected data from
    # the DataFrame
    x_dtype = plot_data[x].dtype
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
        _check_plot_data(ax, plot_data[col], axis)  # type: ignore[arg-type]


class BasePlottingTests:
    # Note: these are ONLY for ExtensionArray subclasses that support plotting.

    @pytest.fixture(params=[True, False], ids=["skipna", "no_skipna"])
    def plot_data(self, data: ExtensionArray, request) -> pd.DataFrame:
        df = pd.DataFrame({"Data": data, "Numeric": np.arange(len(data))})
        if request.param:
            df = df.dropna()
        return df

    def skip_if_no_matplotlib(self):
        """Skips a test if matplotlib dependency not fulfilled.

        Also adds a filter for warnings raised by pyparsing in minimum-version CI for
        pyparsing>=3.3.0 and matplotlib<3.11,
        see https://github.com/matplotlib/matplotlib/pull/29745
        """
        pyparsing = pytest.importorskip(
            "pyparsing", reason="matplotlib requires pyparsing"
        )

        if Version(pyparsing.__version__) >= Version("3.3.0"):
            from pyparsing.warnings import PyparsingDeprecationWarning

            warnings.filterwarnings(
                "ignore",
                message=(
                    "'enablePackrat|oneOf|parseString|resetCache' deprecated - "
                    "use 'enable_packrat|one_of|parse_string|reset_cache'"
                ),
                category=PyparsingDeprecationWarning,
                module=r"matplotlib.*",
            )
        pytest.importorskip("matplotlib", reason="test requires matplotlib")

    def test_plot_on_x_axis(self, plot_data):
        """Test that EA data can be plotted on the x-axis."""
        self.skip_if_no_matplotlib()
        _plot(plot_data, x="Data", y="Numeric")

    def test_plot_on_y_axis(self, plot_data, **kwargs):
        """Test that EA data can be plotted on the y-axis."""
        self.skip_if_no_matplotlib()
        _plot(plot_data, x="Numeric", y="Data", **kwargs)
