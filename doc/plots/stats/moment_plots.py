import numpy as np

import matplotlib.pyplot as plt
import pandas.util.testing as t
import pandas.stats.moments as m

def test_series(n=1000):
    t.N = n
    s = t.makeTimeSeries()
    return s

def plot_timeseries(*args, **kwds):
    n = len(args)

    fig, axes = plt.subplots(n, 1, figsize=kwds.get('size', (10, 5)),
                             sharex=True)
    titles = kwds.get('titles', None)

    for k in range(1, n + 1):
        ax = axes[k-1]
        ts = args[k-1]
        ax.plot(ts.index, ts.values)

        if titles:
            ax.set_title(titles[k-1])

    fig.autofmt_xdate()
    fig.subplots_adjust(bottom=0.10, top=0.95)

