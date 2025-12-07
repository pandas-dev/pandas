"""
Ensure wide DataFrame.line plots use a single LineCollection
instead of one Line2D per column (GH #61764).
"""

import numpy as np
import pytest

import pandas as pd

# Skip this entire module if matplotlib is not installed
mpl = pytest.importorskip("matplotlib")
plt = pytest.importorskip("matplotlib.pyplot")
from matplotlib.collections import LineCollection


def test_linecollection_used_for_wide_dataframe():
    rng = np.random.default_rng(0)
    df = pd.DataFrame(rng.standard_normal((10, 201)).cumsum(axis=0))

    ax = df.plot(legend=False)

    # exactly one LineCollection, and no Line2D artists
    assert sum(isinstance(c, LineCollection) for c in ax.collections) == 1
    assert len(ax.lines) == 0

    plt.close(ax.figure)
