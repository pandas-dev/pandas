"""
Quantilization functions and related stuff
"""

import numpy as np


def cut(data, bins, closed='right', labels=None,
        precision=3):
    """

    Parameters
    ----------
    data : Series or ndarray
    bins : array or integer
        Integer number of bins (equal-size) or array of break points
    closed : {'right', 'left'}, default 'right'
        Side of half-open interval which is closed
    precision : integer, default 3


    """
    bins = np.asarray(bins)
    labels = bins.searchsorted(data.values, side=closed)

    if labels is None or labels is not False:
        if closed == 'left':
            strings = ['[%s, %s)' % tup for tup in zip(bins, bins[1:])]
        else:
            strings = ['(%s, %s]' % tup for tup in zip(bins, bins[1:])]
        strings = np.asarray(strings, dtype=object)
        return strings.take(labels - 1)

    return labels, data.values


def bucket(series, k, by=None):
    """
    Produce DataFrame representing quantiles of a Series

    Parameters
    ----------
    series : Series
    k : int
        number of quantiles
    by : Series or same-length array
        bucket by value

    Returns
    -------
    DataFrame
    """
    if by is None:
        by = series
    else:
        by = by.reindex(series.index)

    split = _split_quantile(by, k)
    mat = np.empty((len(series), k), dtype=float) * np.NaN

    for i, v in enumerate(split):
        mat[:, i][v] = series.take(v)

    return DataFrame(mat, index=series.index, columns=np.arange(k) + 1)

def _split_quantile(arr, k):
    arr = np.asarray(arr)
    mask = np.isfinite(arr)
    order = arr[mask].argsort()
    n = len(arr)

    return np.array_split(np.arange(n)[mask].take(order), k)

def bucketcat(series, cats):
    """
    Produce DataFrame representing quantiles of a Series

    Parameters
    ----------
    series : Series
    cat : Series or same-length array
        bucket by category; mutually exxlusive with 'by'

    Returns
    -------
    DataFrame
    """
    if not isinstance(series, Series):
        series = Series(series, index=np.arange(len(series)))

    cats = np.asarray(cats)

    unique_labels = np.unique(cats)
    unique_labels = unique_labels[notnull(unique_labels)]

    # group by
    data = {}

    for i, label in enumerate(unique_labels):
        data[label] = series[cats == label]

    return DataFrame(data, columns=unique_labels)

def bucketpanel(series, bins=None, by=None, cat=None):
    """
    Bucket data by two Series to create summary panel

    Parameters
    ----------
    series : Series
    bins : tuple (length-2)
        e.g. (2, 2)
    by : tuple of Series
        bucket by value
    cat : tuple of Series
        bucket by category; mutually exxlusive with 'by'

    Returns
    -------
    DataFrame
    """
    use_by = by is not None
    use_cat = cat is not None

    if use_by and use_cat:
        raise Exception('must specify by or cat, but not both')
    elif use_by:
        if len(by) != 2:
            raise Exception('must provide two bucketing series')

        xby, yby = by
        xbins, ybins = bins

        return _bucketpanel_by(series, xby, yby, xbins, ybins)

    elif use_cat:
        xcat, ycat = cat
        return _bucketpanel_cat(series, xcat, ycat)
    else:
        raise Exception('must specify either values or categories to bucket by')

def _bucketpanel_by(series, xby, yby, xbins, ybins):
    xby = xby.reindex(series.index)
    yby = yby.reindex(series.index)

    n = len(series)
    # indices = np.arange(n)

    xlabels = _bucket_labels(xby.reindex(series.index), xbins)
    ylabels = _bucket_labels(yby.reindex(series.index), ybins)

    labels = _uniquify(xlabels, ylabels, xbins, ybins)

    mask = isnull(labels)
    labels[mask] = -1

    unique_labels = np.unique(labels)
    bucketed = bucketcat(series, labels)

    _ulist = list(labels)
    index_map = dict((x, _ulist.index(x)) for x in unique_labels)

    def relabel(key):
        pos = index_map[key]

        xlab = xlabels[pos]
        ylab = ylabels[pos]

        return '%sx%s' % (int(xlab) if notnull(xlab) else 'NULL',
                          int(ylab) if notnull(ylab) else 'NULL')

    return bucketed.rename(columns=relabel)

def _bucketpanel_cat(series, xcat, ycat):
    xlabels, xmapping = _intern(xcat)
    ylabels, ymapping = _intern(ycat)

    shift = 10 ** (np.ceil(np.log10(ylabels.max())))
    labels = xlabels * shift + ylabels

    sorter = labels.argsort()
    sorted_labels = labels.take(sorter)
    sorted_xlabels = xlabels.take(sorter)
    sorted_ylabels = ylabels.take(sorter)

    unique_labels = np.unique(labels)
    unique_labels = unique_labels[notnull(unique_labels)]

    locs = sorted_labels.searchsorted(unique_labels)
    xkeys = sorted_xlabels.take(locs)
    ykeys = sorted_ylabels.take(locs)

    stringified = ['(%s, %s)' % arg
                   for arg in zip(xmapping.take(xkeys), ymapping.take(ykeys))]

    result = bucketcat(series, labels)
    result.columns = stringified

    return result

def _intern(values):
    # assumed no NaN values
    values = np.asarray(values)

    uniqued = np.unique(values)
    labels = uniqued.searchsorted(values)
    return labels, uniqued

def _intern_fast(values):
    pass

def _uniquify(xlabels, ylabels, xbins, ybins):
    # encode the stuff, create unique label
    shifter = 10 ** max(xbins, ybins)
    _xpiece = xlabels * shifter
    _ypiece = ylabels

    return _xpiece + _ypiece

def _cat_labels(labels):
    # group by
    data = {}

    unique_labels = np.unique(labels)
    unique_labels = unique_labels[notnull(unique_labels)]

    for label in unique_labels:
        mask = labels == label
        data[stringified] = series[mask]

    return DataFrame(data, index=series.index)

def _bucket_labels(series, k):
    arr = np.asarray(series)
    mask = np.isfinite(arr)
    order = arr[mask].argsort()
    n = len(series)

    split = np.array_split(np.arange(n)[mask].take(order), k)

    bucketsize = n / k

    mat = np.empty(n, dtype=float) * np.NaN
    for i, v in enumerate(split):
        mat[v] = i

    return mat + 1

def makeQuantiles(series, n):
    """
    Compute quantiles of input series.

    Parameters
    ----------
    series: Series
        Must have 'order' method and index
    n: int
        Number of quantile buckets

    Returns
    -------
    (edges, quantiles)
       edges: ith bucket --> (left edge, right edge)
       quantiles: ith bucket --> set of values
    """
    series = remove_na(series).copy()
    series = series.order()
    quantiles = {}
    edges = {}
    T = float(len(series))
    inc = T / n
    for i in range(n):
        theSlice = series[inc*i:(i+1)*inc]
        quantiles[i+1] = theSlice
        edges[i+1] = theSlice[0], theSlice[-1]
    return edges, quantiles

def quantileTS(frame, percentile):
    """
    Return score at percentile for each point in time (cross-section)

    Parameters
    ----------
    frame: DataFrame
    percentile: int
       nth percentile

    Returns
    -------
    Series (or TimeSeries)
    """
    def func(x):
        x = np.asarray(x.valid())
        if x.any():
            return scoreatpercentile(x, percentile)
        else:
            return NaN
    return frame.apply(func, axis=1)
