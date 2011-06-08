# pylint: disable=E1101,E1103
# pylint: disable=W0212,W0703,W0231,W0622

from cStringIO import StringIO
import sys

from numpy import NaN
import numpy as np

from pandas.core.common import (_pickle_array, _unpickle_array, _try_sort)
from pandas.core.frame import (DataFrame, extract_index, _homogenize_series,
                               _default_index, _ensure_index, _prep_ndarray)
from pandas.core.index import Index, NULL_INDEX
from pandas.core.internals import BlockManager, Block
from pandas.core.series import Series
import pandas.core.common as common
import pandas.core.datetools as datetools
import pandas.lib.tseries as tseries

#-------------------------------------------------------------------------------
# DataMatrix class

class DataMatrix(DataFrame):
    """
    Matrix version of DataFrame, optimized for cross-section operations,
    numerical computation, and other operations that do not require the frame to
    change size.

    Parameters
    ----------
    data : numpy ndarray or dict of sequence-like objects
        Dict can contain Series, arrays, or list-like objects
        Constructor can understand various kinds of inputs
    index : Index or array-like
        Index to use for resulting frame (optional if provided dict of Series)
    columns : Index or array-like
        Required if data is ndarray
    dtype : dtype, default None (infer)
        Data type to force

    Notes
    -----
    Most operations are faster with DataMatrix. You should use it primarily
    unless you are doing a lot of column insertion / deletion (which causes the
    underlying ndarray to have to be reallocated!).
    """
    def __init__(self, data=None, index=None, columns=None):
        if data is None:
            data = {}

        if isinstance(data, BlockManager):
            mgr = data
            index = data.index
            columns = data.columns
        elif isinstance(data, dict):
            index, columns, mgr = _init_dict(data, index, columns)
        elif isinstance(data, (np.ndarray, list)):
            index, columns, mgr = _init_matrix(data, index, columns)
        else:
            raise Exception('DataMatrix constructor not properly called!')

        self._data = mgr

    def _set_columns(self, cols):
        if len(cols) != self.values.shape[1]:
            raise Exception('Columns length %d did not match values %d!' %
                            (len(cols), self.values.shape[1]))

        self._data.columns = _ensure_index(cols)

    def _set_index(self, index):
        if len(index) > 0:
            if len(index) != self.values.shape[0]:
                raise Exception('Index length %d did not match values %d!' %
                                (len(index), self.values.shape[0]))

        self._data.index = _ensure_index(index)

    def _get_index(self):
        return self._data.index

    def _get_columns(self):
        return self._data.columns

    def _get_values(self):
        return self._data.as_matrix()

    def _set_values(self, values):
        raise Exception('Values cannot be assigned to')

    values = property(fget=_get_values)

    def _consolidate_inplace(self):
        self._data = self._data.consolidate()

    def consolidate(self):
        #TODO
        raise NotImplementedError

    @property
    def _constructor(self):
        return DataMatrix

    def __array__(self):
        return self.values

    def __array_wrap__(self, result):
        return DataMatrix(result, index=self.index, columns=self.columns)

#-------------------------------------------------------------------------------
# DataMatrix-specific implementation of private API

    # TODO!
    def _join_on(self, other, on):
        if len(other.index) == 0:
            return self

        if on not in self:
            raise Exception('%s column not contained in this frame!' % on)

        fillVec, mask = tseries.getMergeVec(self[on],
                                            other.index.indexMap)
        notmask = -mask

        tmpMatrix = other.values.take(fillVec, axis=0)
        tmpMatrix[notmask] = NaN

        seriesDict = dict((col, tmpMatrix[:, j])
                           for j, col in enumerate(other.columns))

        if getattr(other, 'objects'):
            objects = other.objects

            tmpMat = objects.values.take(fillVec, axis=0)
            tmpMat[notmask] = NaN
            objDict = dict((col, tmpMat[:, j])
                           for j, col in enumerate(objects.columns))

            seriesDict.update(objDict)

        filledFrame = DataFrame(data=seriesDict, index=self.index)

        return self.join(filledFrame, how='left')

    def _reindex_index(self, new_index, method):
        if new_index is self.index:
            return self.copy()

        # TODO: want to preserve dtypes though...
        new_data = self._data.reindex_index(new_index, method)
        return DataMatrix(new_data)

    def _reindex_columns(self, new_columns):
        if len(new_columns) == 0:
            return DataMatrix(index=self.index)

        new_data = self._data.reindex_columns(new_columns)
        return DataMatrix(new_data)

        # indexer, mask = self.columns.get_indexer(columns)
        # mat = self.values.take(indexer, axis=1)

        # notmask = -mask
        # if len(mask) > 0:
        #     if notmask.any():
        #         if issubclass(mat.dtype.type, np.int_):
        #             mat = mat.astype(float)
        #         elif issubclass(mat.dtype.type, np.bool_):
        #             mat = mat.astype(float)

        #         common.null_out_axis(mat, notmask, 1)

        # return DataMatrix(mat, index=self.index, columns=columns,
        #                   objects=objects)

    def _rename_columns_inplace(self, mapper):
        self.columns = [mapper(x) for x in self.columns]

        if self.objects is not None:
            self.objects._rename_columns_inplace(mapper)

    def _combine_frame(self, other, func):
        """
        Methodology, briefly
        - Really concerned here about speed, space

        - Get new index
        - Reindex to new index
        - Determine new_columns and commonColumns
        - Add common columns over all (new) indices
        - Fill to new set of columns

        Could probably deal with some Cython action in here at some point
        """
        new_index = self._union_index(other)

        if not self and not other:
            return DataMatrix(index=new_index)
        elif not self:
            return other * NaN
        elif not other:
            return self * NaN

        need_reindex = False
        new_columns = self._union_columns(other)
        need_reindex = (need_reindex or new_index is not self.index
                        or new_index is not other.index)
        need_reindex = (need_reindex or new_columns is not self.columns
                        or new_columns is not other.columns)

        this = self
        if need_reindex:
            this = self.reindex(index=new_index, columns=new_columns)
            other = other.reindex(index=new_index, columns=new_columns)

        return DataMatrix(func(this.values, other.values),
                          index=new_index, columns=new_columns)

    def _combine_match_index(self, other, func):
        new_index = self._union_index(other)
        values = self.values
        other_vals = other.values

        # Operate row-wise
        if not other.index.equals(new_index):
            other_vals = other.reindex(new_index).values

        if not self.index.equals(new_index):
            values = self.reindex(new_index).values

        return DataMatrix(func(values.T, other_vals).T,
                          index=new_index, columns=self.columns)

    def _combine_match_columns(self, other, func):
        newCols = self.columns.union(other.index)

        # Operate column-wise
        this = self.reindex(columns=newCols)
        other = other.reindex(newCols).values

        return DataMatrix(func(this.values, other),
                          index=self.index, columns=newCols)

    def _combine_const(self, other, func):
        if not self:
            return self

        # TODO: deal with objects
        return DataMatrix(func(self.values, other), index=self.index,
                          columns=self.columns)

#-------------------------------------------------------------------------------
# "Magic methods"

    def __getstate__(self):
        if self.objects is not None:
            objects = self.objects._matrix_state(pickle_index=False)
        else:
            objects = None

        state = self._matrix_state()

        return (state, objects)

    def _matrix_state(self, pickle_index=True):
        columns = _pickle_array(self.columns)

        if pickle_index:
            index = _pickle_array(self.index)
        else:
            index = None

        return self.values, index, columns

    def __setstate__(self, state):
        (vals, idx, cols), object_state = state

        self.values = vals
        self.index = _unpickle_array(idx)
        self.columns = _unpickle_array(cols)

        if object_state:
            ovals, _, ocols = object_state
            self.objects = DataMatrix(ovals,
                                      index=self.index,
                                      columns=_unpickle_array(ocols))
        else:
            self.objects = None

    def __getitem__(self, item):
        """
        Retrieve column, slice, or subset from DataMatrix.

        Possible inputs
        ---------------
        single value : retrieve a column as a Series
        slice : reindex to indices specified by slice
        boolean vector : like slice but more general, reindex to indices
          where the input vector is True

        Examples
        --------
        column = dm['A']

        dmSlice = dm[:20] # First 20 rows

        dmSelect = dm[dm.count(axis=1) > 10]

        Notes
        -----
        This is a magic method. Do NOT call explicity.
        """
        if isinstance(item, slice):
            new_data = self._data.get_slice(item)
            return DataMatrix(new_data)
        elif isinstance(item, np.ndarray):
            if len(item) != len(self.index):
                raise Exception('Item wrong length %d instead of %d!' %
                                (len(item), len(self.index)))
            new_index = self.index[item]
            return self.reindex(new_index)
        else:
            values = self._data.get(item)
            return Series(values, index=self.index)

    # __setitem__ logic

    def _boolean_set(self, key, value):
        mask = key.values
        if mask.dtype != np.bool_:
            raise Exception('Must pass DataFrame with boolean values only')

        self.values[mask] = value

    def _insert_item(self, key, value):
        """
        Add series to DataMatrix in specified column.

        If series is a numpy-array (not a Series/TimeSeries), it must be the
        same length as the DataMatrix's index or an error will be thrown.

        Series/TimeSeries will be conformed to the DataMatrix's index to
        ensure homogeneity.
        """
        if hasattr(value, '__iter__'):
            if isinstance(value, Series):
                if value.index.equals(self.index):
                    # no need to copy
                    value = value.values
                else:
                    value = value.reindex(self.index).values
            else:
                assert(len(value) == len(self.index))

                if not isinstance(value, np.ndarray):
                    value = np.array(value)
                    if value.dtype.type == np.str_:
                        value = np.array(value, dtype=object)
        else:
            value = np.repeat(value, len(self.index))

        self._data.set(key, value)

    def __delitem__(self, key):
        """
        Delete column from DataMatrix
        """
        self._data.delete(key)

    # to support old APIs
    @property
    def _series(self):
        return self._data.get_series_dict(self.index)

#-------------------------------------------------------------------------------
# Public methods

    def apply(self, func, axis=0, broadcast=False):
        """
        Applies func to columns (Series) of this DataMatrix and returns either
        a DataMatrix (if the function produces another series) or a Series
        indexed on the column names of the DataFrame if the function produces
        a value.

        Parameters
        ----------
        func : function
            Function to apply to each column
        broadcast : bool, default False
            For aggregation functions, return object of same size with values
            propagated

        Examples
        --------
        >>> df.apply(numpy.sqrt) --> DataMatrix
        >>> df.apply(numpy.sum) --> Series

        N.B.: Do NOT use functions that might toy with the index.
        """
        if not len(self.columns):
            return self

        if isinstance(func, np.ufunc):
            results = func(self.values)
            return DataMatrix(data=results, index=self.index,
                              columns=self.columns)
        else:
            return DataFrame.apply(self, func, axis=axis,
                                   broadcast=broadcast)

    def applymap(self, func):
        """
        Apply a function to a DataMatrix that is intended to operate
        elementwise, i.e. like doing
            map(func, series) for each series in the DataMatrix

        Parameters
        ----------
        func : function
            Python function, returns a single value from a single value

        Note : try to avoid using this function if you can, very slow.
        """
        npfunc = np.frompyfunc(func, 1, 1)
        results = npfunc(self.values)
        try:
            results = results.astype(self.values.dtype)
        except Exception:
            pass

        return DataMatrix(results, index=self.index, columns=self.columns)

    def append(self, other):
        """
        Glue together DataFrame objects having non-overlapping indices

        Parameters
        ----------
        other : DataFrame
        """
        if not other:
            return self.copy()

        if not self:
            return other.copy()

        # TODO: with blocks
        # idx = Index(np.concatenate([self.index, other.index]))
        # mat = np.vstack((self.values, other.values))
        return DataFrame.append(self, other)

    def asMatrix(self, columns=None):
        """
        Convert the DataMatrix to its Numpy-array matrix representation

        Columns are presented in sorted order unless a specific list
        of columns is provided.

        Parameters
        ----------
        columns : list-like
            columns to use in producing matrix, must all be contained

        Returns
        -------
        ndarray
        """
        return self._data.as_matrix(columns=columns)

    def copy(self):
        """
        Make a copy of this DataMatrix
        """
        return DataMatrix(self._data.copy())

    def cumsum(self, axis=0):
        """
        Return DataMatrix of cumulative sums over requested axis.

        Parameters
        ----------
        axis : {0, 1}
            0 for row-wise, 1 for column-wise

        Returns
        -------
        y : DataMatrix
        """
        y = np.array(self.values, subok=True)
        if not issubclass(y.dtype.type, np.int_):
            mask = np.isnan(self.values)
            y[mask] = 0
            result = y.cumsum(axis)
            has_obs = (-mask).astype(int).cumsum(axis) > 0
            result[-has_obs] = np.NaN
        else:
            result = y.cumsum(axis)

        return DataMatrix(result, index=self.index,
                          columns=self.columns, objects=self.objects)

    def min(self, axis=0):
        """
        Return array or Series of minimums over requested axis.

        Parameters
        ----------
        axis : {0, 1}
            0 for row-wise, 1 for column-wise

        Returns
        -------
        Series or TimeSeries
        """
        values = self.values.copy()
        np.putmask(values, -np.isfinite(values), np.inf)
        return Series(values.min(axis), index=self._get_agg_axis(axis))

    def max(self, axis=0):
        """
        Return array or Series of maximums over requested axis.

        Parameters
        ----------
        axis : {0, 1}
            0 for row-wise, 1 for column-wise

        Returns
        -------
        Series or TimeSeries
        """
        values = self.values.copy()
        np.putmask(values, -np.isfinite(values), -np.inf)
        return Series(values.max(axis), index=self._get_agg_axis(axis))

    def fillna(self, value=None, method='pad'):
        """
        Fill NaN values using the specified method.

        Member Series / TimeSeries are filled separately.

        Parameters
        ----------
        value : any kind (should be same type as array)
            Value to use to fill holes (e.g. 0)

        method : {'backfill', 'pad', None}
            Method to use for filling holes in new inde

        Returns
        -------
        y : DataMatrix

        See also
        --------
        DataMatrix.reindex, DataMatrix.asfreq
        """
        if value is None:
            result = {}
            series = self._series
            for col, s in series.iteritems():
                result[col] = s.fillna(method=method, value=value)

            return DataMatrix(result, index=self.index, objects=self.objects)
        else:
            # Float type values
            if len(self.columns) == 0:
                return self

            vals = self.values.copy()
            vals.flat[common.isnull(vals.ravel())] = value

            objects = None

            if self.objects is not None:
                objects = self.objects.copy()

            return DataMatrix(vals, index=self.index, columns=self.columns,
                              objects=objects)

    def xs(self, key, copy=True):
        """
        Returns a row from the DataMatrix as a Series object.

        Parameters
        ----------
        key : some index contained in the index

        Returns
        -------
        Series
        """
        if key not in self.index:
            raise Exception('No cross-section for %s' % key)

        self._consolidate_inplace()

        loc = self.index.get_loc(key)
        xs = self._data.xs(loc, copy=copy)
        return Series(xs, index=self.columns)

    @property
    def T(self):
        """
        Returns a DataMatrix with the rows/columns switched.
        """
        return DataMatrix(data=self.values.T, index=self.columns,
                          columns=self.index)

    def shift(self, periods, offset=None, timeRule=None):
        """
        Shift the underlying series of the DataMatrix and Series objects within
        by given number (positive or negative) of periods.

        Parameters
        ----------
        periods : int (+ or -)
            Number of periods to move
        offset : DateOffset, optional
            Increment to use from datetools module
        timeRule : string
            Time rule to use by name

        Returns
        -------
        DataMatrix
        """
        if periods == 0:
            return self

        if timeRule is not None and offset is None:
            offset = datetools.getOffset(timeRule)

        if offset is None:
            indexer = self._shift_indexer(periods)
            new_values = self.values.take(indexer, axis=0)
            new_index = self.index

            new_values = common.ensure_float(new_values)

            if periods > 0:
                new_values[:periods] = NaN
            else:
                new_values[periods:] = NaN
        else:
            new_index = self.index.shift(periods, offset)
            new_values = self.values.copy()

        if self.objects is not None:
            shifted_objects = self.objects.shift(periods, offset=offset,
                                                 timeRule=timeRule)

            shifted_objects.index = new_index
        else:
            shifted_objects = None

        return DataMatrix(data=new_values, index=new_index,
                          columns=self.columns, objects=shifted_objects)

_data_types = [np.float_, np.int_]

def _filter_out(data, columns):
    if columns is not None:
        colset = set(columns)
        data = dict((k, v) for k, v in data.iteritems() if k in colset)

    return data


def _group_dtypes(data, columns):
    import itertools

    chunk_cols = []
    chunks = []
    for dtype, gp_cols in itertools.groupby(columns, lambda x: data[x].dtype):
        chunk = np.vstack([data[k] for k in gp_cols]).T

        chunks.append(chunk)
        chunk_cols.append(gp_cols)

    return chunks, chunk_cols

def _init_dict(data, index, columns):
    """
    Segregate Series based on type and coerce into matrices.

    Needs to handle a lot of exceptional cases.

    Somehow this got outrageously complicated
    """
    # figure out the index, if necessary
    if index is None:
        index = extract_index(data)

    # TODO: deal with emptiness!
    # TODO: dtype casting?

    # don't force copy because getting jammed in an ndarray anyway
    homogenized = _homogenize_series(data, index, force_copy=False)
    # segregates dtypes and forms blocks matching to columns
    blocks, columns = _form_blocks(homogenized, columns)
    mgr = BlockManager(blocks, index, columns)
    return index, columns, mgr

def _form_blocks(data, columns):
    # pre-filter out columns if we passed it
    if columns is None:
        columns = Index(_try_sort(data.keys()))
        extra_columns = NULL_INDEX
    else:
        columns = _ensure_index(columns)
        extra_columns = columns - Index(data.keys())

    # prefilter
    data = dict((k, v) for k, v in data.iteritems() if k in columns)

    # put "leftover" columns in float bucket, where else?
    # generalize?
    float_dict = {}
    object_dict = {}
    for k, v in data.iteritems():
        if issubclass(v.dtype.type, (np.floating, np.integer)):
            float_dict[k] = v
        else:
            object_dict[k] = v

    blocks = []

    # oof
    if len(float_dict) > 0 or len(extra_columns) > 0:
        if len(extra_columns):
            bcolumns = extra_columns.union(float_dict.keys())
            float_block = _blockify(float_dict, np.float64,
                                    columns=bcolumns)
        else:
            float_block = _blockify(float_dict, np.float64)
        blocks.append(float_block)

    if len(object_dict) > 0:
        object_block = _blockify(object_dict, np.object_)
        blocks.append(object_block)

    return blocks, columns

def _blockify(dct, dtype, columns=None):
    dict_columns = Index(_try_sort(dct))
    stacked_series = np.vstack([dct[k] for k in dict_columns]).T
    if columns is None:
        columns = dict_columns
        values = stacked_series
        # CHECK DTYPE
    else:
        n = len(dct.values()[0])
        k = len(columns)
        values = np.empty((n, k), dtype=dtype)

        indexer, mask = columns.get_indexer(dict_columns)
        assert(mask.all())
        values[:, indexer] = stacked_series

    # do something with dtype?
    return Block(values, columns)

def _init_matrix(values, index, columns):
    values = _prep_ndarray(values)

    if values.ndim == 1:
        N = values.shape[0]
        if N == 0:
            values = values.reshape((values.shape[0], 0))
        else:
            values = values.reshape((values.shape[0], 1))

    # address this can of worms later
    # if dtype is not None:
    #     try:
    #         values = values.astype(dtype)
    #     except Exception:
    #         pass

    N, K = values.shape

    if index is None:
        index = _default_index(N)

    if columns is None:
        columns = _default_index(K)

    columns = _ensure_index(columns)
    block = Block(values, columns)
    mgr = BlockManager([block], index, columns)

    return index, columns, mgr

def _reorder_columns(mat, current, desired):
    indexer, mask = common.get_indexer(current, desired, None)
    return mat.take(indexer[mask], axis=1)

if __name__ == '__main__':
    pass
