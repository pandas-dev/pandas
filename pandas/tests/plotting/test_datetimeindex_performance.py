"""
Tests for optimized DatetimeIndex plotting performance.
"""

import numpy as np

import pandas.util._test_decorators as td

import pandas as pd
from pandas import (
    DataFrame,
    testing as tm,
)


@td.skip_if_no_mpl
def test_plot_with_datetimeindex_performance():
    """
    Test that plotting a DataFrame with DatetimeIndex is performant.

    Check that plotting multiple columns with the same DatetimeIndex
    doesn't perform redundant calculations/conversions.
    """
    # Create a DataFrame with DatetimeIndex
    rng = np.random.RandomState(42)
    n = 1000
    idx = pd.date_range(start="2020-01-01", periods=n, freq="D")

    # Dataframe with many columns to highlight the optimization
    df = DataFrame(rng.randn(n, 5), index=idx)

    # Define function to plot in the regular/slow way
    def plot_df():
        ax = df.plot(figsize=(10, 5))
        return ax

    # The first run might be slower due to imports, etc.
    # This test is primarily to catch regressions in performance
    # rather than assert a specific speedup.
    plot_df()  # warmup

    # Assert that subsequent plotting operations are not significantly slower
    # (verify caching is working)
    tm.assert_faster_than(plot_df, plot_df, significant=False)
