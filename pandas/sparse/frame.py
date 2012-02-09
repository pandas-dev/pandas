"""
Data structures for sparse float data. Life is made simpler by dealing only with
float64 data
"""

# pylint: disable=E1101,E1103,W0231,E0202

from numpy import nan
import numpy as np

from pandas.core.common import _pickle_array, _unpickle_array, _try_sort
from pandas.core.index import Index, MultiIndex, NULL_INDEX, _ensure_index
from pandas.core.series import Series
from pandas.core.frame import (DataFrame, extract_index, _prep_ndarray,
                               _default_index)
from pandas.util.decorators import cache_readonly
import pandas.core.common as com
import pandas.core.datetools as datetools

from pandas.sparse.series import SparseSeries
from pandas.util.decorators import Appender


class _SparseMockBlockManager(object):

    def __init__(self, sp_frame):
        self.sp_frame = sp_frame

    def get(self, item):
        return self.sp_frame[item].values

    @property
    def shape(self):
        x, y = self.sp_frame.shape
        return y, x

    @property
    def axes(self):
        return [self.sp_frame.columns, self.sp_frame.index]

class SparseDataFrame(DataFrame):
    """
    DataFrame containing sparse floating point data in the form of SparseSeries
    objects

    Parameters
    ----------
    data : same types as can be passed to DataFrame
    index : array-like, optional
    column : array-like, optional
    default_kind : {'block', 'integer'}, default 'block'
        Default sparse kind for converting Series to SparseSeries. Will not
        override SparseSeries passed into constructor
    default_fill_value : float
        Default fill_value for converting Series to SparseSeries. Will not
        override SparseSeries passed in
    """
    _verbose_info = False
    _columns = None
    _series = None
    _is_mixed_type = False
    ndim = 2

    def __init__(self, data=None, index=None, columns=None,
                 default_kind='block', default_fill_value=None):
        if default_fill_value is None:
            default_fill_value = np.nan

        self.default_kind = default_kind
        self.default_fill_value = default_fill_value

        if isinstance(data, dict):
            sdict, columns, index = self._init_dict(data, index, columns)
        elif isinstance(data, (np.ndarray, list)):
            sdict, columns, index = self._init_matrix(data, index, columns)
        elif isinstance(data, DataFrame):
            sdict, columns, index = self._init_dict(data, data.index,
                                                    data.columns)
        elif data is None:
            sdict = {}

            if index is None:
                index = NULL_INDEX

            if columns is None:
                columns = NULL_INDEX
            else:
                for c in columns:
                    sdict[c] = Series(np.nan, index=index)

        self._series = sdict
        self.columns = columns
        self.index = index

    def _from_axes(self, data, axes):
        columns, index = axes
        return self._constructor(data, index=index, columns=columns)

    @cache_readonly
    def _data(self):
        return _SparseMockBlockManager(self)

    def _consolidate_inplace(self):
        # do nothing when DataFrame calls this method
        pass

    def convert_objects(self):
        # XXX
        return self

    @property
    def _constructor(self):
        def wrapper(data, index=None, columns=None):
            return SparseDataFrame(data, index=index, columns=columns,
                                   default_fill_value=self.default_fill_value,
                                   default_kind=self.default_kind)
        return wrapper

    def _init_dict(self, data, index, columns, dtype=None):
        # pre-filter out columns if we passed it
        if columns is not None:
            columns = _ensure_index(columns)
            data = dict((k, v) for k, v in data.iteritems() if k in columns)
        else:
            columns = Index(_try_sort(data.keys()))

        if index is None:
            index = extract_index(data)

        sp_maker = lambda x: SparseSeries(x, index=index,
                                          kind=self.default_kind,
                                          fill_value=self.default_fill_value,
                                          copy=True)

        sdict = {}
        for k, v in data.iteritems():
            if isinstance(v, Series):
                # Force alignment, no copy necessary
                if not v.index.equals(index):
                    v = v.reindex(index)

                if not isinstance(v, SparseSeries):
                    v = sp_maker(v)
            else:
                if isinstance(v, dict):
                    v = [v.get(i, nan) for i in index]

                v = sp_maker(v)
            sdict[k] = v

        # TODO: figure out how to handle this case, all nan's?
        # add in any other columns we want to have (completeness)
        nan_vec = np.empty(len(index))
        nan_vec.fill(nan)
        for c in columns:
            if c not in sdict:
                sdict[c] = sp_maker(nan_vec)

        return sdict, columns, index

    def _init_matrix(self, data, index, columns, dtype=None):
        data = _prep_ndarray(data, copy=False)
        N, K = data.shape
        if index is None:
            index = _default_index(N)
        if columns is None:
            columns = _default_index(K)

        if len(columns) != K:
            raise Exception('Column length mismatch: %d vs. %d' %
                            (len(columns), K))
        if len(index) != N:
            raise Exception('Index length mismatch: %d vs. %d' %
                            (len(index), N))

        data = dict([(idx, data[:, i]) for i, idx in enumerate(columns)])
        return self._init_dict(data, index, columns, dtype)

    def __array_wrap__(self, result):
        return SparseDataFrame(result, index=self.index, columns=self.columns,
                               default_kind=self.default_kind,
                               default_fill_value=self.default_fill_value)

    def __getstate__(self):
        series = dict((k, (v.sp_index, v.sp_values))
                      for k, v in self.iteritems())
        columns = _pickle_array(self.columns)
        index = _pickle_array(self.index)

        return (series, columns, index, self.default_fill_value,
                self.default_kind)

    def __setstate__(self, state):
        series, cols, idx, fv, kind = state
        columns = _unpickle_array(cols)
        index = _unpickle_array(idx)

        series_dict = {}
        for col, (sp_index, sp_values) in series.iteritems():
            series_dict[col] = SparseSeries(sp_values, sparse_index=sp_index,
                                            fill_value=fv)

        self._series = series_dict
        self.index = index
        self.columns = columns
        self.default_fill_value = fv
        self.default_kind = kind

    def to_dense(self):
        """
        Convert to dense DataFrame

        Returns
        -------
        df : DataFrame
        """
        data = dict((k, v.to_dense()) for k, v in self.iteritems())
        return DataFrame(data, index=self.index)

    def astype(self, dtype):
        raise NotImplementedError

    def copy(self, deep=True):
        """
        Make a copy of this SparseDataFrame
        """
        series = dict((k, v.copy()) for k, v in self.iteritems())
        return SparseDataFrame(series, index=self.index, columns=self.columns,
                               default_fill_value=self.default_fill_value,
                               default_kind=self.default_kind)

    @property
    def density(self):
        """
        Ratio of non-sparse points to total (dense) data points
        represented in the frame
        """
        tot_nonsparse = sum([ser.sp_index.npoints
                             for _, ser in self.iteritems()])
        tot = len(self.index) * len(self.columns)
        return tot_nonsparse / float(tot)

    #----------------------------------------------------------------------
    # Support different internal rep'n of SparseDataFrame

    def _set_item(self, key, value):
        sp_maker = lambda x: SparseSeries(x, index=self.index,
                                          fill_value=self.default_fill_value,
                                          kind=self.default_kind)
        if hasattr(value, '__iter__'):
            if isinstance(value, Series):
                clean_series = value.reindex(self.index)
                if not isinstance(value, SparseSeries):
                    clean_series = sp_maker(clean_series)
            else:
                clean_series = sp_maker(value)

            self._series[key] = clean_series
        # Scalar
        else:
            self._series[key] = sp_maker(value)

        if key not in self.columns:
            self._insert_column(key)

    def _insert_column(self, key):
        self.columns = Index(np.concatenate((self.columns, [key])))

    def __delitem__(self, key):
        """
        Delete column from DataFrame
        """
        loc = self.columns.get_loc(key)
        del self._series[key]
        self._delete_column_index(loc)

    def _delete_column_index(self, loc):
        if loc == len(self.columns) - 1:
            new_columns = self.columns[:loc]
        else:
            new_columns = Index(np.concatenate((self.columns[:loc],
                                               self.columns[loc+1:])))
        self.columns = new_columns

    _index = None
    def _set_index(self, index):
        self._index = _ensure_index(index)
        for v in self._series.values():
            v.index = self._index

    def _get_index(self):
        return self._index

    def _get_columns(self):
        return self._columns

    def _set_columns(self, cols):
        if len(cols) != len(self._series):
            raise Exception('Columns length %d did not match data %d!' %
                            (len(cols), len(self._series)))
        self._columns = _ensure_index(cols)

    index = property(fget=_get_index, fset=_set_index)
    columns = property(fget=_get_columns, fset=_set_columns)

    def __getitem__(self, item):
        """
        Retrieve column or slice from DataFrame
        """
        try:
            # unsure about how kludgy this is
            s = self._series[item]
            s.name = item
            return s
        except (TypeError, KeyError):
            if isinstance(item, slice):
                dateRange = self.index[item]
                return self.reindex(dateRange)

            elif isinstance(item, np.ndarray):
                if len(item) != len(self.index):
                    raise Exception('Item wrong length %d instead of %d!' %
                                    (len(item), len(self.index)))
                newIndex = self.index[item]
                return self.reindex(newIndex)
            else: # pragma: no cover
                raise

    @Appender(DataFrame.get_value.__doc__, indents=0)
    def get_value(self, index, col):
        s = self._series[col]
        return s.get_value(index)

    def set_value(self, index, col, value):
        """
        Put single value at passed column and index

        Parameters
        ----------
        index : row label
        col : column label
        value : scalar value

        Notes
        -----
        This method *always* returns a new object. It is currently not
        particularly efficient (and potentially very expensive) but is provided
        for API compatibility with DataFrame

        Returns
        -------
        frame : DataFrame
        """
        dense = self.to_dense().set_value(index, col, value)
        return dense.to_sparse(kind=self.default_kind,
                               fill_value=self.default_fill_value)

    def _slice(self, slobj, axis=0):
        if axis == 0:
            new_index = self.index[slobj]
            new_columns = self.columns
        else:
            new_index = self.index
            new_columns = self.columns[slobj]

        return self.reindex(index=new_index, columns=new_columns)

    def as_matrix(self, columns=None):
        """
        Convert the frame to its Numpy-array matrix representation

        Columns are presented in sorted order unless a specific list
        of columns is provided.
        """
        if columns is None:
            columns = self.columns

        if len(columns) == 0:
            return np.zeros((len(self.index), 0), dtype=float)

        return np.array([self[col].values for col in columns]).T

    values = property(as_matrix)

    def xs(self, key, axis=0, copy=False):
        """
        Returns a row (cross-section) from the SparseDataFrame as a Series
        object.

        Parameters
        ----------
        key : some index contained in the index

        Returns
        -------
        xs : Series
        """
        if axis == 1:
            data = self[key]
            return data

        i = self.index.get_loc(key)
        series = self._series
        values = [series[k][i] for k in self.columns]
        return Series(values, index=self.columns)

    #----------------------------------------------------------------------
    # Arithmetic-related methods

    def _combine_frame(self, other, func, fill_value=None, level=None):
        this, other = self.align(other, join='outer', level=level,
                                 copy=False)
        new_index, new_columns = this.index, this.columns

        if fill_value is not None or level is not None:
            raise NotImplementedError

        if not self and not other:
            return SparseDataFrame(index=new_index)

        new_data = {}
        for col in new_columns:
            if col in this and col in other:
                new_data[col] = func(this[col], other[col])

        return self._constructor(data=new_data, index=new_index,
                                 columns=new_columns)

    def _combine_match_index(self, other, func, fill_value=None):
        new_data = {}

        if fill_value is not None:
            raise NotImplementedError

        new_index = self.index.union(other.index)
        this = self
        if self.index is not new_index:
            this = self.reindex(new_index)

        if other.index is not new_index:
            other = other.reindex(new_index)

        for col, series in this.iteritems():
            new_data[col] = func(series.values, other.values)

        return self._constructor(new_data, index=new_index,
                                 columns=self.columns)

    def _combine_match_columns(self, other, func, fill_value):
        # patched version of DataFrame._combine_match_columns to account for
        # NumPy circumventing __rsub__ with float64 types, e.g.: 3.0 - series,
        # where 3.0 is numpy.float64 and series is a SparseSeries. Still
        # possible for this to happen, which is bothersome

        if fill_value is not None:
            raise NotImplementedError

        new_data = {}

        union = intersection = self.columns

        if not union.equals(other.index):
            union = other.index.union(self.columns)
            intersection = other.index.intersection(self.columns)

        for col in intersection:
            new_data[col] = func(self[col], float(other[col]))

        return self._constructor(new_data, index=self.index,
                                 columns=union)

    def _combine_const(self, other, func):
        new_data = {}
        for col, series in self.iteritems():
            new_data[col] = func(series, other)

        return self._constructor(data=new_data, index=self.index,
                                 columns=self.columns)

    def _reindex_index(self, index, method, copy, level):
        if level is not None:
            raise Exception('Reindex by level not supported for sparse')

        if self.index.equals(index):
            if copy:
                return self.copy()
            else:
                return self

        if len(self.index) == 0:
            return SparseDataFrame(index=index, columns=self.columns)

        indexer = self.index.get_indexer(index, method)
        mask = indexer == -1
        need_mask = mask.any()

        new_series = {}
        for col, series in self.iteritems():
            values = series.values
            new = values.take(indexer)

            if need_mask:
                np.putmask(new, mask, nan)

            new_series[col] = new

        return SparseDataFrame(new_series, index=index, columns=self.columns,
                               default_fill_value=self.default_fill_value)

    def _reindex_columns(self, columns, copy, level):
        if level is not None:
            raise Exception('Reindex by level not supported for sparse')

        # TODO: fill value handling
        sdict = dict((k, v) for k, v in self.iteritems() if k in columns)
        return SparseDataFrame(sdict, index=self.index, columns=columns,
                               default_fill_value=self.default_fill_value)

    def _reindex_with_indexers(self, index, row_indexer, columns, col_indexer,
                               copy):
        if columns is None:
            columns = self.columns

        new_arrays = {}
        for col in columns:
            if col not in self:
                continue
            if row_indexer is not None:
                new_arrays[col] = com.take_1d(self[col].values, row_indexer)
            else:
                new_arrays[col] = self[col]

        return self._constructor(new_arrays, index=index, columns=columns)

    def _rename_index_inplace(self, mapper):
        self.index = [mapper(x) for x in self.index]

    def _rename_columns_inplace(self, mapper):
        new_series = {}
        new_columns = []

        for col in self.columns:
            new_col = mapper(col)
            if new_col in new_series: # pragma: no cover
                raise Exception('Non-unique mapping!')
            new_series[new_col] = self[col]
            new_columns.append(new_col)

        self.columns = new_columns
        self._series = new_series

    def take(self, indices, axis=0):
        """
        Analogous to ndarray.take, return SparseDataFrame corresponding to
        requested indices along an axis

        Parameters
        ----------
        indices : list / array of ints
        axis : {0, 1}

        Returns
        -------
        taken : SparseDataFrame
        """
        new_values = self.values.take(indices, axis=axis)
        if axis == 0:
            new_columns = self.columns
            new_index = self.index.take(indices)
        else:
            new_columns = self.columns.take(indices)
            new_index = self.index
        return self._constructor(new_values, index=new_index,
                                 columns=new_columns)

    def add_prefix(self, prefix):
        f = (('%s' % prefix) + '%s').__mod__
        return self.rename(columns=f)

    def add_suffix(self, suffix):
        f = ('%s' + ('%s' % suffix)).__mod__
        return self.rename(columns=f)

    def _join_compat(self, other, on=None, how='left', lsuffix='', rsuffix='',
                     sort=False):
        if on is not None:
            raise NotImplementedError
        else:
            return self._join_index(other, how, lsuffix, rsuffix)

    def _join_index(self, other, how, lsuffix, rsuffix):
        if isinstance(other, Series):
            assert(other.name is not None)
            other = SparseDataFrame({other.name : other},
                                    default_fill_value=self.default_fill_value)

        join_index = self.index.join(other.index, how=how)

        this = self.reindex(join_index)
        other = other.reindex(join_index)

        this, other = this._maybe_rename_join(other, lsuffix, rsuffix)

        result_series = this._series
        other_series = other._series
        result_series.update(other_series)

        return self._constructor(result_series, index=join_index)

    def _maybe_rename_join(self, other, lsuffix, rsuffix):
        intersection = self.columns.intersection(other.columns)

        if len(intersection) > 0:
            if not lsuffix and not rsuffix:
                raise Exception('columns overlap: %s' % intersection)

            def lrenamer(x):
                if x in intersection:
                    return '%s%s' % (x, lsuffix)
                return x

            def rrenamer(x):
                if x in intersection:
                    return '%s%s' % (x, rsuffix)
                return x

            this = self.rename(columns=lrenamer)
            other = other.rename(columns=rrenamer)
        else:
            this = self

        return this, other

    def transpose(self):
        """
        Returns a DataFrame with the rows/columns switched.
        """
        return SparseDataFrame(self.values.T, index=self.columns,
                               columns=self.index,
                               default_fill_value=self.default_fill_value,
                               default_kind=self.default_kind)
    T = property(transpose)

    @Appender(DataFrame.count.__doc__)
    def count(self, axis=0, **kwds):
        return self.apply(lambda x: x.count(), axis=axis)

    def cumsum(self, axis=0):
        """
        Return SparseDataFrame of cumulative sums over requested axis.

        Parameters
        ----------
        axis : {0, 1}
            0 for row-wise, 1 for column-wise

        Returns
        -------
        y : SparseDataFrame
        """
        return self.apply(lambda x: x.cumsum(), axis=axis)

    def shift(self, periods, offset=None, timeRule=None):
        """
        Analogous to DataFrame.shift
        """
        if timeRule is not None and offset is None:
            offset = datetools.getOffset(timeRule)

        new_series = {}
        if offset is None:
            new_index = self.index
            for col, s in self.iteritems():
                new_series[col] = s.shift(periods)
        else:
            new_index = self.index.shift(periods, offset)
            for col, s in self.iteritems():
                new_series[col] = SparseSeries(s.sp_values, index=new_index,
                                               sparse_index=s.sp_index,
                                               fill_value=s.fill_value)

        return SparseDataFrame(new_series, index=new_index,
                               columns=self.columns,
                               default_fill_value=self.default_fill_value,
                               default_kind=self.default_kind)

    def apply(self, func, axis=0, broadcast=False):
        """
        Analogous to DataFrame.apply, for SparseDataFrame

        Parameters
        ----------
        func : function
            Function to apply to each column
        axis : {0, 1}
        broadcast : bool, default False
            For aggregation functions, return object of same size with values
            propagated

        Returns
        -------
        applied : Series or SparseDataFrame
        """
        if not len(self.columns):
            return self

        if isinstance(func, np.ufunc):
            new_series = {}
            for k, v in self.iteritems():
                applied = func(v)
                applied.fill_value = func(applied.fill_value)
                new_series[k] = applied
            return SparseDataFrame(new_series, index=self.index,
                                   columns=self.columns,
                                   default_fill_value=self.default_fill_value,
                                   default_kind=self.default_kind)
        else:
            if not broadcast:
                return self._apply_standard(func, axis)
            else:
                return self._apply_broadcast(func, axis)

    def fillna(self, *args, **kwargs):
        raise NotImplementedError

def stack_sparse_frame(frame):
    """
    Only makes sense when fill_value is NaN
    """
    lengths = [s.sp_index.npoints for _, s in frame.iteritems()]
    nobs = sum(lengths)

    # this is pretty fast
    minor_labels = np.repeat(np.arange(len(frame.columns)), lengths)

    inds_to_concat = []
    vals_to_concat = []
    for _, series in frame.iteritems():
        if not np.isnan(series.fill_value):
            raise Exception('This routine assumes NaN fill value')

        int_index = series.sp_index.to_int_index()
        inds_to_concat.append(int_index.indices)
        vals_to_concat.append(series.sp_values)

    major_labels = np.concatenate(inds_to_concat)
    stacked_values = np.concatenate(vals_to_concat)
    index = MultiIndex(levels=[frame.index, frame.columns],
                       labels=[major_labels, minor_labels])

    lp = DataFrame(stacked_values.reshape((nobs, 1)), index=index,
                   columns=['foo'])
    return lp.sortlevel(level=0)


def homogenize(series_dict):
    """
    Conform a set of SparseSeries (with NaN fill_value) to a common SparseIndex
    corresponding to the locations where they all have data

    Parameters
    ----------
    series_dict : dict or DataFrame

    Notes
    -----
    Using the dumbest algorithm I could think of. Should put some more thought
    into this

    Returns
    -------
    homogenized : dict of SparseSeries
    """
    index = None

    need_reindex = False

    for _, series in series_dict.iteritems():
        if not np.isnan(series.fill_value):
            raise Exception('this method is only valid with NaN fill values')

        if index is None:
            index = series.sp_index
        elif not series.sp_index.equals(index):
            need_reindex = True
            index = index.intersect(series.sp_index)

    if need_reindex:
        output = {}
        for name, series in series_dict.iteritems():
            if not series.sp_index.equals(index):
                series = series.sparse_reindex(index)

            output[name] = series
    else:
        output = series_dict

    return output
