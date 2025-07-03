"""
Ensure wide DataFrame.line plots use a single LineCollection
instead of one Line2D per column (PR #61764).
"""

from matplotlib.collections import LineCollection
import matplotlib.pyplot as plt
import numpy as np

import pandas as pd


def test_linecollection_used_for_wide_dataframe():
    rng = np.random.default_rng(0)
    df = pd.DataFrame(rng.standard_normal((10, 201)).cumsum(axis=0))

    ax = df.plot(legend=False)

    # one LineCollection, zero Line2D objects
    assert sum(isinstance(c, LineCollection) for c in ax.collections) == 1
    assert len(ax.lines) == 0

    plt.close(ax.figure)
