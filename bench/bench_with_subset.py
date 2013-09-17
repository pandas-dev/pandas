#!/usr/bin/env python

"""
Microbenchmarks for comparison with R's "with" and "subset" functions
"""

from __future__ import print_function
import numpy as np
from numpy import array
from timeit import repeat as timeit
from pandas.compat import range, zip
from pandas import DataFrame


setup_common = """from pandas import DataFrame
from numpy.random import randn
df = DataFrame(randn(%d, 3), columns=list('abc'))
%s"""


setup_with = "s = 'a + b * (c ** 2 + b ** 2 - a) / (a * c) ** 3'"


def bench_with(n, times=10, repeat=3, engine='numexpr'):
    return np.array(timeit('df.eval(s, engine=%r)' % engine,
                           setup=setup_common % (n, setup_with),
                           repeat=repeat, number=times)) / times


setup_subset = "s = 'a <= b <= c ** 2 + b ** 2 - a and b > c'"


def bench_subset(n, times=10, repeat=3, engine='numexpr'):
    return np.array(timeit('df.query(s, engine=%r)' % engine,
                           setup=setup_common % (n, setup_subset),
                           repeat=repeat, number=times)) / times


def bench(mn=1, mx=7, num=100, engines=('python', 'numexpr'), verbose=False):
    r = np.logspace(mn, mx, num=num).round().astype(int)

    ev = DataFrame(np.empty((num, len(engines))), columns=engines)
    qu = ev.copy(deep=True)

    ev['size'] = qu['size'] = r

    for engine in engines:
        for i, n in enumerate(r):
            if verbose:
                print('engine: %r, i == %d' % (engine, i))
            ev.loc[i, engine] = bench_with(n, times=1, repeat=1, engine=engine)
            qu.loc[i, engine] = bench_subset(n, times=1, repeat=1,
                                             engine=engine)

    return ev, qu


def plot_perf(df, engines, title, filename=None):
    from matplotlib.pyplot import figure, rc

    try:
        from mpltools import style
    except ImportError:
        pass
    else:
        style.use('ggplot')

    rc('text', usetex=True)

    fig = figure(figsize=(4, 3), dpi=100)
    ax = fig.add_subplot(111)

    for engine in engines:
        ax.plot(df.size, df[engine], label=engine, lw=2)

    ax.set_xlabel('Number of Rows')
    ax.set_ylabel('Time (s)')
    ax.set_title(title)
    ax.legend(loc='best')
    ax.tick_params(top=False, right=False)

    fig.tight_layout()

    if filename is not None:
        fig.savefig(filename)


if __name__ == '__main__':
    import os
    import pandas as pd

    pandas_dir = os.path.dirname(os.path.abspath(os.path.dirname(__file__)))
    static_path = os.path.join(pandas_dir, 'doc', 'source', '_static')

    join = lambda p: os.path.join(static_path, p)

    fn = join('eval-query-perf-data.h5')

    engines = 'python', 'numexpr'

    if not os.path.exists(fn):
        ev, qu = bench(verbose=True)
        ev.to_hdf(fn, 'eval')
        qu.to_hdf(fn, 'query')
    else:
        ev = pd.read_hdf(fn, 'eval')
        qu = pd.read_hdf(fn, 'query')

    plot_perf(ev, engines, 'DataFrame.eval()', filename=join('eval-perf.png'))
    plot_perf(qu, engines, 'DataFrame.query()',
              filename=join('query-perf.png'))

    plot_perf(ev[ev.size <= 50000], engines, 'DataFrame.eval()',
              filename=join('eval-perf-small.png'))
    plot_perf(qu[qu.size <= 500000], engines, 'DataFrame.query()',
              filename=join('query-perf-small.png'))
