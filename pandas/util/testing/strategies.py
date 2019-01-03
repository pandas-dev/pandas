"""
strategies for creating test data.

`strategies` name derived from `hypothesis`
"""
from datetime import datetime
import string
import warnings

import numpy as np

import pandas.compat as compat
from pandas.compat import lrange, lmap, Counter, lzip, u

from pandas.core.dtypes.common import is_sequence

import pandas as pd


N = 30
K = 4


RANDS_CHARS = np.array(list(string.ascii_letters + string.digits),
                       dtype=(np.str_, 1))
RANDU_CHARS = np.array(list(u("").join(map(unichr, lrange(1488, 1488 + 26))) +
                            string.digits), dtype=(np.unicode_, 1))


def rands_array(nchars, size, dtype='O'):
    """Generate an array of byte strings."""
    retval = (np.random.choice(RANDS_CHARS, size=nchars * np.prod(size))
              .view((np.str_, nchars)).reshape(size))
    if dtype is None:
        return retval
    else:
        return retval.astype(dtype)


def randu_array(nchars, size, dtype='O'):
    """Generate an array of unicode strings."""
    retval = (np.random.choice(RANDU_CHARS, size=nchars * np.prod(size))
              .view((np.unicode_, nchars)).reshape(size))
    if dtype is None:
        return retval
    else:
        return retval.astype(dtype)


def getCols(k):
    return string.ascii_uppercase[:k]


def makeStringIndex(k=10, name=None):
    return pd.Index(rands_array(nchars=10, size=k), name=name)


def makeUnicodeIndex(k=10, name=None):
    return pd.Index(randu_array(nchars=10, size=k), name=name)


def makeCategoricalIndex(k=10, n=3, name=None, **kwargs):
    """ make a length k index or n categories """
    x = rands_array(nchars=4, size=n)
    return pd.CategoricalIndex(np.random.choice(x, k), name=name, **kwargs)


def makeIntervalIndex(k=10, name=None, **kwargs):
    """ make a length k IntervalIndex """
    x = np.linspace(0, 100, num=(k + 1))
    return pd.IntervalIndex.from_breaks(x, name=name, **kwargs)


def makeBoolIndex(k=10, name=None):
    if k == 1:
        return pd.Index([True], name=name)
    elif k == 2:
        return pd.Index([False, True], name=name)
    return pd.Index([False, True] + [False] * (k - 2), name=name)


def makeIntIndex(k=10, name=None):
    return pd.Index(lrange(k), name=name)


def makeUIntIndex(k=10, name=None):
    return pd.Index([2**63 + i for i in lrange(k)], name=name)


def makeRangeIndex(k=10, name=None, **kwargs):
    return pd.RangeIndex(0, k, 1, name=name, **kwargs)


def makeFloatIndex(k=10, name=None):
    values = sorted(np.random.random_sample(k)) - np.random.random_sample(1)
    return pd.Index(values * (10 ** np.random.randint(0, 9)), name=name)


def makeDateIndex(k=10, freq='B', name=None, **kwargs):
    dt = datetime(2000, 1, 1)
    dr = pd.bdate_range(dt, periods=k, freq=freq, name=name)
    return pd.DatetimeIndex(dr, name=name, **kwargs)


def makeTimedeltaIndex(k=10, freq='D', name=None, **kwargs):
    return pd.timedelta_range(start='1 day', periods=k, freq=freq,
                              name=name, **kwargs)


def makePeriodIndex(k=10, name=None, **kwargs):
    dt = datetime(2000, 1, 1)
    dr = pd.period_range(start=dt, periods=k, freq='B', name=name, **kwargs)
    return dr


def makeMultiIndex(k=10, names=None, **kwargs):
    return pd.MultiIndex.from_product(
        (('foo', 'bar'), (1, 2)), names=names, **kwargs)


def makeFloatSeries(name=None):
    index = makeStringIndex(N)
    return pd.Series(np.random.randn(N), index=index, name=name)


def makeStringSeries(name=None):
    index = makeStringIndex(N)
    return pd.Series(np.random.randn(N), index=index, name=name)


def makeObjectSeries(name=None):
    dateIndex = makeDateIndex(N)
    dateIndex = pd.Index(dateIndex, dtype=object)
    index = makeStringIndex(N)
    return pd.Series(dateIndex, index=index, name=name)


def getSeriesData():
    index = makeStringIndex(N)
    return {c: pd.Series(np.random.randn(N), index=index) for c in getCols(K)}


def makeTimeSeries(nper=None, freq='B', name=None):
    if nper is None:
        nper = N
    return pd.Series(np.random.randn(nper),
                     index=makeDateIndex(nper, freq=freq),
                     name=name)


def makePeriodSeries(nper=None, name=None):
    if nper is None:
        nper = N
    return pd.Series(np.random.randn(nper),
                     index=makePeriodIndex(nper),
                     name=name)


def getTimeSeriesData(nper=None, freq='B'):
    return {c: makeTimeSeries(nper, freq) for c in getCols(K)}


def getPeriodData(nper=None):
    return {c: makePeriodSeries(nper) for c in getCols(K)}


def makeTimeDataFrame(nper=None, freq='B'):
    data = getTimeSeriesData(nper, freq)
    return pd.DataFrame(data)


def makeDataFrame():
    data = getSeriesData()
    return pd.DataFrame(data)


def getMixedTypeDict():
    index = pd.Index(['a', 'b', 'c', 'd', 'e'])

    data = {
        'A': [0., 1., 2., 3., 4.],
        'B': [0., 1., 0., 1., 0.],
        'C': ['foo1', 'foo2', 'foo3', 'foo4', 'foo5'],
        'D': pd.bdate_range('1/1/2009', periods=5)
    }

    return index, data


def makeMixedDataFrame():
    return pd.DataFrame(getMixedTypeDict()[1])


def makePeriodFrame(nper=None):
    data = getPeriodData(nper)
    return pd.DataFrame(data)


def makePanel(nper=None):
    with warnings.catch_warnings(record=True):
        warnings.filterwarnings("ignore", "\\nPanel", FutureWarning)
        cols = ['Item' + c for c in string.ascii_uppercase[:K - 1]]
        data = {c: makeTimeDataFrame(nper) for c in cols}
        return pd.Panel.fromDict(data)


def makePeriodPanel(nper=None):
    with warnings.catch_warnings(record=True):
        warnings.filterwarnings("ignore", "\\nPanel", FutureWarning)
        cols = ['Item' + c for c in string.ascii_uppercase[:K - 1]]
        data = {c: makePeriodFrame(nper) for c in cols}
        return pd.Panel.fromDict(data)


def makeCustomIndex(nentries, nlevels, prefix='#', names=False, ndupe_l=None,
                    idx_type=None):
    """Create an index/multindex with given dimensions, levels, names, etc'

    nentries - number of entries in index
    nlevels - number of levels (> 1 produces multindex)
    prefix - a string prefix for labels
    names - (Optional), bool or list of strings. if True will use default
       names, if false will use no names, if a list is given, the name of
       each level in the index will be taken from the list.
    ndupe_l - (Optional), list of ints, the number of rows for which the
       label will repeated at the corresponding level, you can specify just
       the first few, the rest will use the default ndupe_l of 1.
       len(ndupe_l) <= nlevels.
    idx_type - "i"/"f"/"s"/"u"/"dt"/"p"/"td".
       If idx_type is not None, `idx_nlevels` must be 1.
       "i"/"f" creates an integer/float index,
       "s"/"u" creates a string/unicode index
       "dt" create a datetime index.
       "td" create a datetime index.

        if unspecified, string labels will be generated.
    """

    if ndupe_l is None:
        ndupe_l = [1] * nlevels
    assert (is_sequence(ndupe_l) and len(ndupe_l) <= nlevels)
    assert (names is None or names is False or
            names is True or len(names) is nlevels)
    assert idx_type is None or (idx_type in ('i', 'f', 's', 'u',
                                             'dt', 'p', 'td')
                                and nlevels == 1)

    if names is True:
        # build default names
        names = [prefix + str(i) for i in range(nlevels)]
    if names is False:
        # pass None to index constructor for no name
        names = None

    # make singelton case uniform
    if isinstance(names, compat.string_types) and nlevels == 1:
        names = [names]

    # specific 1D index type requested?
    idx_func = dict(i=makeIntIndex, f=makeFloatIndex,
                    s=makeStringIndex, u=makeUnicodeIndex,
                    dt=makeDateIndex, td=makeTimedeltaIndex,
                    p=makePeriodIndex).get(idx_type)
    if idx_func:
        idx = idx_func(nentries)
        # but we need to fill in the name
        if names:
            idx.name = names[0]
        return idx
    elif idx_type is not None:
        raise ValueError('"{idx_type}" is not a legal value for `idx_type`, '
                         'use  "i"/"f"/"s"/"u"/"dt/"p"/"td".'
                         .format(idx_type=idx_type))

    if len(ndupe_l) < nlevels:
        ndupe_l.extend([1] * (nlevels - len(ndupe_l)))
    assert len(ndupe_l) == nlevels

    assert all(x > 0 for x in ndupe_l)

    tuples = []
    for i in range(nlevels):
        def keyfunc(x):
            import re
            numeric_tuple = re.sub(r"[^\d_]_?", "", x).split("_")
            return lmap(int, numeric_tuple)

        # build a list of lists to create the index from
        div_factor = nentries // ndupe_l[i] + 1
        cnt = Counter()
        for j in range(div_factor):
            label = '{prefix}_l{i}_g{j}'.format(prefix=prefix, i=i, j=j)
            cnt[label] = ndupe_l[i]
        # cute Counter trick
        result = list(sorted(cnt.elements(), key=keyfunc))[:nentries]
        tuples.append(result)

    tuples = lzip(*tuples)

    # convert tuples to index
    if nentries == 1:
        # we have a single level of tuples, i.e. a regular Index
        index = pd.Index(tuples[0], name=names[0])
    elif nlevels == 1:
        name = None if names is None else names[0]
        index = pd.Index((x[0] for x in tuples), name=name)
    else:
        index = pd.MultiIndex.from_tuples(tuples, names=names)
    return index


def makeCustomDataframe(nrows, ncols, c_idx_names=True, r_idx_names=True,
                        c_idx_nlevels=1, r_idx_nlevels=1, data_gen_f=None,
                        c_ndupe_l=None, r_ndupe_l=None, dtype=None,
                        c_idx_type=None, r_idx_type=None):
    """
   nrows,  ncols - number of data rows/cols
   c_idx_names, idx_names  - False/True/list of strings,  yields No names ,
        default names or uses the provided names for the levels of the
        corresponding index. You can provide a single string when
        c_idx_nlevels ==1.
   c_idx_nlevels - number of levels in columns index. > 1 will yield MultiIndex
   r_idx_nlevels - number of levels in rows index. > 1 will yield MultiIndex
   data_gen_f - a function f(row,col) which return the data value
        at that position, the default generator used yields values of the form
        "RxCy" based on position.
   c_ndupe_l, r_ndupe_l - list of integers, determines the number
        of duplicates for each label at a given level of the corresponding
        index. The default `None` value produces a multiplicity of 1 across
        all levels, i.e. a unique index. Will accept a partial list of length
        N < idx_nlevels, for just the first N levels. If ndupe doesn't divide
        nrows/ncol, the last label might have lower multiplicity.
   dtype - passed to the DataFrame constructor as is, in case you wish to
        have more control in conjuncion with a custom `data_gen_f`
   r_idx_type, c_idx_type -  "i"/"f"/"s"/"u"/"dt"/"td".
       If idx_type is not None, `idx_nlevels` must be 1.
       "i"/"f" creates an integer/float index,
       "s"/"u" creates a string/unicode index
       "dt" create a datetime index.
       "td" create a timedelta index.

        if unspecified, string labels will be generated.

    Examples:

    # 5 row, 3 columns, default names on both, single index on both axis
    >> makeCustomDataframe(5,3)

    # make the data a random int between 1 and 100
    >> mkdf(5,3,data_gen_f=lambda r,c:randint(1,100))

    # 2-level multiindex on rows with each label duplicated
    # twice on first level, default names on both axis, single
    # index on both axis
    >> a=makeCustomDataframe(5,3,r_idx_nlevels=2,r_ndupe_l=[2])

    # DatetimeIndex on row, index with unicode labels on columns
    # no names on either axis
    >> a=makeCustomDataframe(5,3,c_idx_names=False,r_idx_names=False,
                             r_idx_type="dt",c_idx_type="u")

    # 4-level multindex on rows with names provided, 2-level multindex
    # on columns with default labels and default names.
    >> a=makeCustomDataframe(5,3,r_idx_nlevels=4,
                             r_idx_names=["FEE","FI","FO","FAM"],
                             c_idx_nlevels=2)

    >> a=mkdf(5,3,r_idx_nlevels=2,c_idx_nlevels=4)
    """

    assert c_idx_nlevels > 0
    assert r_idx_nlevels > 0
    assert r_idx_type is None or (r_idx_type in ('i', 'f', 's',
                                                 'u', 'dt', 'p', 'td')
                                  and r_idx_nlevels == 1)
    assert c_idx_type is None or (c_idx_type in ('i', 'f', 's',
                                                 'u', 'dt', 'p', 'td')
                                  and c_idx_nlevels == 1)

    columns = makeCustomIndex(ncols, nlevels=c_idx_nlevels, prefix='C',
                              names=c_idx_names, ndupe_l=c_ndupe_l,
                              idx_type=c_idx_type)
    index = makeCustomIndex(nrows, nlevels=r_idx_nlevels, prefix='R',
                            names=r_idx_names, ndupe_l=r_ndupe_l,
                            idx_type=r_idx_type)

    # by default, generate data based on location
    if data_gen_f is None:
        data_gen_f = lambda r, c: "R{rows}C{cols}".format(rows=r, cols=c)

    data = [[data_gen_f(r, c) for c in range(ncols)] for r in range(nrows)]

    return pd.DataFrame(data, index, columns, dtype=dtype)


def _create_missing_idx(nrows, ncols, density, random_state=None):
    if random_state is None:
        random_state = np.random
    else:
        random_state = np.random.RandomState(random_state)

    # below is cribbed from scipy.sparse
    size = int(np.round((1 - density) * nrows * ncols))
    # generate a few more to ensure unique values
    min_rows = 5
    fac = 1.02
    extra_size = min(size + min_rows, fac * size)

    def _gen_unique_rand(rng, _extra_size):
        ind = rng.rand(int(_extra_size))
        return np.unique(np.floor(ind * nrows * ncols))[:size]

    ind = _gen_unique_rand(random_state, extra_size)
    while ind.size < size:
        extra_size *= 1.05
        ind = _gen_unique_rand(random_state, extra_size)

    j = np.floor(ind * 1. / nrows).astype(int)
    i = (ind - j * nrows).astype(int)
    return i.tolist(), j.tolist()


def makeMissingCustomDataframe(nrows, ncols, density=.9, random_state=None,
                               c_idx_names=True, r_idx_names=True,
                               c_idx_nlevels=1, r_idx_nlevels=1,
                               data_gen_f=None,
                               c_ndupe_l=None, r_ndupe_l=None, dtype=None,
                               c_idx_type=None, r_idx_type=None):
    """
    Parameters
    ----------
    Density : float, optional
        Float in (0, 1) that gives the percentage of non-missing numbers in
        the DataFrame.
    random_state : {np.random.RandomState, int}, optional
        Random number generator or random seed.

    See makeCustomDataframe for descriptions of the rest of the parameters.
    """
    df = makeCustomDataframe(nrows, ncols, c_idx_names=c_idx_names,
                             r_idx_names=r_idx_names,
                             c_idx_nlevels=c_idx_nlevels,
                             r_idx_nlevels=r_idx_nlevels,
                             data_gen_f=data_gen_f,
                             c_ndupe_l=c_ndupe_l, r_ndupe_l=r_ndupe_l,
                             dtype=dtype, c_idx_type=c_idx_type,
                             r_idx_type=r_idx_type)

    i, j = _create_missing_idx(nrows, ncols, density, random_state)
    df.values[i, j] = np.nan
    return df


def makeMissingDataframe(density=.9, random_state=None):
    df = makeDataFrame()
    i, j = _create_missing_idx(*df.shape, density=density,
                               random_state=random_state)
    df.values[i, j] = np.nan
    return df
