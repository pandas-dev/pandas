"""
DataFrame
---------
An efficient 2D container for potentially mixed-type time series or other
labeled data series.

Similar to its R counterpart, data.frame, except providing automatic data
alignment and a host of useful data manipulation methods having to do with the
labeling information
"""

# pylint: disable=E1101,E1103
# pylint: disable=W0212,W0231,W0703,W0622

from cStringIO import StringIO
from datetime import datetime
import operator
import sys
import warnings

from numpy import nan
import numpy as np

from pandas.core.common import (isnull, notnull, PandasError, _ensure_index,
                                _try_sort, _pfixed, _default_index,
                                _infer_dtype)
from pandas.core.daterange import DateRange
from pandas.core.generic import AxisProperty, NDFrame
from pandas.core.index import Index, MultiIndex, NULL_INDEX
from pandas.core.internals import BlockManager, make_block, form_blocks
from pandas.core.series import Series, _is_bool_indexer
import pandas.core.common as common
import pandas.core.datetools as datetools
import pandas._tseries as _tseries

#-------------------------------------------------------------------------------
# Factory helper methods

_arith_doc ="""
Arithmetic method: %s

Parameters
----------
other : Series, DataFrame, or constant
axis : {0, 1, 'index', 'columns'}
    For Series input, axis to match Series index on
fill_value : None or float value, default None
    Fill missing (NaN) values with this value. If both DataFrame locations are
    missing, the result will be missing

Notes
-----
Mismatched indices will be unioned together

Returns
-------
result : DataFrame
"""

def _arith_method(func, name, default_axis='columns'):
    def f(self, other, axis=default_axis, fill_value=None):
        if isinstance(other, DataFrame):    # Another DataFrame
            return self._combine_frame(other, func, fill_value)
        elif isinstance(other, Series):
            if axis is not None:
                axis = self._get_axis_name(axis)
                if axis == 'index':
                    return self._combine_match_index(other, func, fill_value)
                else:
                    return self._combine_match_columns(other, func, fill_value)
            return self._combine_series_infer(other, func, fill_value)
        else:
            return self._combine_const(other, func)

    f.__name__ = name
    f.__doc__ = _arith_doc % name

    return f

def comp_method(func, name):
    def f(self, other):
        if isinstance(other, DataFrame):    # Another DataFrame
            return self._compare_frame(other, func)
        elif isinstance(other, Series):
            return self._combine_series_infer(other, func)
        else:
            return self._combine_const(other, func)

    f.__name__ = name
    f.__doc__ = 'Wrapper for comparison method %s' % name

    return f

#-------------------------------------------------------------------------------
# DataFrame class

class DataFrame(NDFrame):
    """
    Homogenously indexed table with named columns, with intelligent arithmetic
    operations, slicing, reindexing, aggregation, etc. Can function
    interchangeably as a dictionary.

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
    copy : boolean, default False
        Copy data from inputs. Only affects DataFrame / 2d ndarray input

    Examples
    --------
        >>> d = {'col1' : ts1, 'col2' : ts2}
        >>> df = DataFrame(data=d, index=someIndex)
    """
    _auto_consolidate = True

    _AXIS_NUMBERS = {
        'index' : 0,
        'columns' : 1
    }

    _AXIS_NAMES = dict((v, k) for k, v in _AXIS_NUMBERS.iteritems())

    def __init__(self, data=None, index=None, columns=None, dtype=None,
                 copy=False):

        if data is None:
            data = {}

        if isinstance(data, DataFrame):
            data = data._data

        if isinstance(data, BlockManager):
            # do not copy BlockManager unless explicitly done
            mgr = data
            if copy and dtype is None:
                mgr = mgr.copy()
            elif dtype is not None:
                # no choice but to copy
                mgr = mgr.cast(dtype)
        elif isinstance(data, dict):
            mgr = self._init_dict(data, index, columns, dtype=dtype)
        elif isinstance(data, np.ndarray):
            if data.dtype.names:
                data_columns, data = _rec_to_dict(data)
                if columns is None:
                    columns = data_columns
                mgr = self._init_dict(data, index, columns, dtype=dtype)
            else:
                mgr = self._init_ndarray(data, index, columns, dtype=dtype,
                                         copy=copy)
        elif isinstance(data, list):
            mgr = self._init_ndarray(data, index, columns, dtype=dtype,
                                     copy=copy)
        else:
            raise PandasError('DataFrame constructor not properly called!')

        self._data = mgr

    def _init_dict(self, data, index, columns, dtype=None):
        """
        Segregate Series based on type and coerce into matrices.
        Needs to handle a lot of exceptional cases.
        """
        # prefilter if columns passed
        if columns is not None:
            columns = _ensure_index(columns)
            data = dict((k, v) for k, v in data.iteritems() if k in columns)
        else:
            columns = Index(_try_sort(data.keys()))

        # figure out the index, if necessary
        if index is None:
            index = extract_index(data)
        else:
            index = _ensure_index(index)

        # don't force copy because getting jammed in an ndarray anyway
        homogenized = _homogenize(data, index, columns, dtype)

        # from BlockManager perspective
        axes = [columns, index]

        # segregates dtypes and forms blocks matching to columns
        blocks = form_blocks(homogenized, axes)

        # consolidate for now
        mgr = BlockManager(blocks, axes)
        return mgr.consolidate()

    def _init_ndarray(self, values, index, columns, dtype=None,
                      copy=False):
        values = _prep_ndarray(values, copy=copy)

        if dtype is not None:
            try:
                values = values.astype(dtype)
            except Exception:
                raise ValueError('failed to cast to %s' % dtype)

        N, K = values.shape

        if index is None:
            index = _default_index(N)

        if columns is None:
            columns = _default_index(K)

        columns = _ensure_index(columns)
        block = make_block(values.T, columns, columns)
        return BlockManager([block], [columns, index])

    def astype(self, dtype):
        """
        Cast DataFrame to input numpy.dtype

        Parameters
        ----------
        dtype : numpy.dtype or Python type

        Returns
        -------
        casted : DataFrame
        """
        return self._constructor(self._data, dtype=dtype)

    def _wrap_array(self, arr, axes, copy=False):
        index, columns = axes
        return self._constructor(arr, index=index, columns=columns, copy=copy)

    @property
    def axes(self):
        return [self.index, self.columns]

    @property
    def _constructor(self):
        return DataFrame

    #----------------------------------------------------------------------
    # Class behavior

    def __nonzero__(self):
        # e.g. "if frame: ..."
        return len(self.columns) > 0 and len(self.index) > 0

    def __repr__(self):
        """
        Return a string representation for a particular DataFrame
        """
        buf = StringIO()
        if len(self.index) < 500 and len(self.columns) < 10:
            self.toString(buf=buf)
        else:
            self.info(buf=buf)

        return buf.getvalue()

    def __iter__(self):
        """
        Iterate over columns of the frame.
        """
        return iter(self.columns)

    def iteritems(self):
        """Iterator over (column, series) pairs"""
        series = self._series
        return ((k, series[k]) for k in self.columns)

    def __len__(self):
        """
        Returns number of columns/Series inside
        """
        return len(self.index)

    def __contains__(self, key):
        """
        True if DataFrame has this column
        """
        return key in self.columns

    def copy(self):
        """
        Make a copy of this DataFrame
        """
        return self._constructor(self._data.copy())

    #----------------------------------------------------------------------
    # Arithmetic methods

    add = _arith_method(operator.add, 'add')
    mul = _arith_method(operator.mul, 'multiply')
    sub = _arith_method(operator.sub, 'subtract')
    div = _arith_method(operator.div, 'divide')

    radd = _arith_method(operator.add, 'add')
    rmul = _arith_method(operator.mul, 'multiply')
    rsub = _arith_method(lambda x, y: y - x, 'subtract')
    rdiv = _arith_method(lambda x, y: y / x, 'divide')

    __add__ = _arith_method(operator.add, '__add__', default_axis=None)
    __sub__ = _arith_method(operator.sub, '__sub__', default_axis=None)
    __mul__ = _arith_method(operator.mul, '__mul__', default_axis=None)
    __div__ = _arith_method(operator.div, '__div__', default_axis=None)
    __truediv__ = _arith_method(operator.truediv, '__truediv__',
                               default_axis=None)
    __pow__ = _arith_method(operator.pow, '__pow__', default_axis=None)

    __radd__ = _arith_method(operator.add, '__radd__', default_axis=None)
    __rmul__ = _arith_method(operator.mul, '__rmul__', default_axis=None)
    __rsub__ = _arith_method(lambda x, y: y - x, '__rsub__', default_axis=None)
    __rdiv__ = _arith_method(lambda x, y: y / x, '__rdiv__', default_axis=None)
    __rtruediv__ = _arith_method(lambda x, y: y / x, '__rtruediv__',
                                default_axis=None)
    __rpow__ = _arith_method(lambda x, y: y ** x, '__rpow__', default_axis=None)

    def __neg__(self):
        return self * -1

    #----------------------------------------------------------------------
    # Comparison methods

    __eq__ = comp_method(operator.eq, '__eq__')
    __ne__ = comp_method(operator.ne, '__ne__')
    __lt__ = comp_method(operator.lt, '__lt__')
    __gt__ = comp_method(operator.gt, '__gt__')
    __le__ = comp_method(operator.le, '__le__')
    __ge__ = comp_method(operator.ge, '__ge__')

    #----------------------------------------------------------------------
    # IO methods (to / from other formats)

    def toDict(self):
        """
        Convert DataFrame to nested dictionary (non-pandas)

        Return
        ------
        nested dict mapping: {column -> index -> value}
        """
        return dict((k, v.toDict()) for k, v in self.iteritems())

    @classmethod
    def from_records(cls, data, indexField=None):
        """
        Convert structured or record ndarray to DataFrame

        Parameters
        ----------
        input : NumPy structured array

        Returns
        -------
        DataFrame
        """
        if not data.dtype.names:
            raise Exception('Input was not a structured array!')

        columns, sdict = _rec_to_dict(data)
        if indexField is not None:
            index = sdict.pop(indexField)
            columns.remove(indexField)
        else:
            index = np.arange(len(data))

        return cls(sdict, index=index, columns=columns)

    def to_records(self, index=True):
        """
        Convert DataFrame to record array. Index will be put in the
        'index' field of the record array.

        Returns
        -------
        y : recarray
        """
        if index:
            arrays = [self.index] + [self[c] for c in self.columns]
            names = ['index'] + list(self.columns)
        else:
            arrays = [self[c] for c in self.columns]
            names = list(self.columns)

        return np.rec.fromarrays(arrays, names=names)

    @classmethod
    def fromcsv(cls, path, header=0, delimiter=',', index_col=0):
        """
        Read delimited file into DataFrame

        Parameters
        ----------
        path : string
        header : int, default 0
            Row to use at header (skip prior rows)
        delimiter : string, default ','
        index_col : int
            Column to use for index

        Notes
        -----
        Will attempt to convert index to datetimes for time series
        data. Uses numpy.genfromtxt to do the actual parsing into
        ndarray

        Returns
        -------
        y : DataFrame or DataFrame
        """
        from pandas.io.parsers import read_table
        df = read_table(path, header=header, sep=delimiter,
                        index_col=index_col)
        return df

    def to_sparse(self, fill_value=None, kind='block'):
        """
        Convert to SparseDataFrame

        Parametpers
        ----------
        fill_value : float, default NaN
        kind : {'block', 'integer'}

        Returns
        -------
        y : SparseDataFrame
        """
        from pandas.core.sparse import SparseDataFrame
        return SparseDataFrame(self._series, index=self.index,
                               default_kind=kind, default_fill_value=fill_value)

    def toCSV(self, path, nanRep='', cols=None, header=True,
              index=True, mode='wb'):
        """
        Write the DataFrame to a CSV file

        Parameters
        ----------
        path : string
            File path
        nanRep : string, default ''
            Missing data rep'n
        cols : sequence, optional
        header : boolean, default True
            Write out column names
        index : boolean, default True
            Write row names (index)
        """
        f = open(path, mode)

        if cols is None:
            cols = self.columns

        series = self._series
        if header:
            joined_cols = ','.join([str(c) for c in cols])
            if index:
                # this could be dangerous
                f.write('index,%s' % joined_cols)
            else:
                f.write(joined_cols)
            f.write('\n')

        for idx in self.index:
            if index:
                f.write(str(idx))
            for i, col in enumerate(cols):
                val = series[col].get(idx)
                if isnull(val):
                    val = nanRep
                else:
                    val = str(val)

                if i > 0 or index:
                    f.write(',%s' % val)
                else:
                    f.write('%s' % val)

            f.write('\n')

        f.close()

    def toString(self, buf=sys.stdout, columns=None, colSpace=None,
                 nanRep='NaN', formatters=None, float_format=None):
        from pandas.core.common import _format, adjoin

        if colSpace is None:
            def _myformat(v):
                return _format(v, nanRep=nanRep,
                               float_format=float_format)
        else:
            def _myformat(v):
                return _pfixed(v, colSpace, nanRep=nanRep,
                               float_format=float_format)

        if formatters is None:
            formatters = {}

        def _stringify(col):
            formatter = formatters.get(col, _myformat)
            return [formatter(x) for x in self[col]]

        if columns is None:
            columns = self.columns
        else:
            columns = [c for c in columns if c in self]

        if len(columns) == 0 or len(self.index) == 0:
            print >> buf, 'Empty %s' % type(self).__name__
            print >> buf, repr(self.index)
        else:
            fmt_index = self.index.format().split('\n')
            fmt_columns = self.columns
            str_index = [''] + fmt_index
            stringified = [[' %s' % str(c)] + _stringify(c) for c in columns]
            print >> buf, adjoin(2, str_index, *stringified)

    def info(self, verbose=True, buf=sys.stdout):
        """
        Concise summary of a DataFrame, used in __repr__ when very large.
        """
        print >> buf, str(type(self))
        print >> buf, self.index.summary()

        if len(self.columns) == 0:
            print >> buf, 'Empty %s' % type(self).__name__
            return

        cols = self.columns

        if verbose:
            print >> buf, 'Data columns:'
            space = max([len(str(k)) for k in self.columns]) + 4
            col_counts = []
            counts = self.count()
            assert(len(cols) == len(counts))
            for col, count in counts.iteritems():
                col_counts.append('%s%d  non-null values' %
                                  (_put_str(str(col), space), count))

            print >> buf, '\n'.join(col_counts)
        else:
            if len(cols) <= 2:
                print >> buf, 'Columns: %s' % repr(cols)
            else:
                print >> buf, 'Columns: %s to %s' % (cols[0], cols[-1])

        counts = self._get_dtype_counts()
        dtypes = ['%s(%d)' % k for k in sorted(counts.iteritems())]
        buf.write('dtypes: %s' % ', '.join(dtypes))

    def _get_dtype_counts(self):
        counts = {}
        for _, series in self.iteritems():
            if series.dtype in counts:
                counts[series.dtype] += 1
            else:
                counts[series.dtype] = 1

        return counts

    #----------------------------------------------------------------------
    # properties for index and columns

    # reference underlying BlockManager
    columns = AxisProperty(0)
    index = AxisProperty(1)

    def as_matrix(self, columns=None):
        """
        Convert the frame to its Numpy-array matrix representation

        Columns are presented in sorted order unless a specific list
        of columns is provided.
        """
        self._consolidate_inplace()
        return self._data.as_matrix(columns).T

    values = property(fget=as_matrix)

    def transpose(self):
        """
        Returns a DataFrame with the rows/columns switched. Copy of data is not
        made by default
        """
        return self._constructor(data=self.values.T, index=self.columns,
                                 columns=self.index, copy=False)
    T = property(transpose)

    #----------------------------------------------------------------------
    # Picklability

    def __getstate__(self):
        return self._data

    def __setstate__(self, state):
        # old DataFrame pickle
        if isinstance(state, BlockManager):
            self._data = state
        elif isinstance(state[0], dict): # pragma: no cover
            self._unpickle_frame_compat(state)
        else: # pragma: no cover
            # old pickling format, for compatibility
            self._unpickle_matrix_compat(state)

    def _unpickle_frame_compat(self, state): # pragma: no cover
        from pandas.core.common import _unpickle_array
        if len(state) == 2: # pragma: no cover
            series, idx = state
            columns = sorted(series)
        else:
            series, cols, idx = state
            columns = _unpickle_array(cols)

        index = _unpickle_array(idx)
        self._data = self._init_dict(series, index, columns, None)

    def _unpickle_matrix_compat(self, state): # pragma: no cover
        from pandas.core.common import _unpickle_array
        # old unpickling
        (vals, idx, cols), object_state = state

        index = _unpickle_array(idx)
        dm = DataFrame(vals, index=index, columns=_unpickle_array(cols),
                       copy=False)

        if object_state is not None:
            ovals, _, ocols = object_state
            objects = DataFrame(ovals, index=index,
                                columns=_unpickle_array(ocols),
                                copy=False)

            dm = dm.join(objects)

        self._data = dm._data

    #----------------------------------------------------------------------
    # Private helper methods

    def _intersect_index(self, other):
        common_index = self.index

        if not common_index.equals(other.index):
            common_index = common_index.intersection(other.index)

        return common_index

    def _intersect_columns(self, other):
        common_cols = self.columns

        if not common_cols.equals(other.columns):
            common_cols = common_cols.intersection(other.columns)

        return common_cols

    #----------------------------------------------------------------------
    # Array interface

    def __array__(self, dtype=None):
        return self.values

    def __array_wrap__(self, result):
        return self._constructor(result, index=self.index, columns=self.columns,
                                 copy=False)

    #----------------------------------------------------------------------
    # getitem/setitem related

    def __getitem__(self, item):
        """
        Retrieve column, slice, or subset from DataFrame.

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
            new_data = self._data.get_slice(item, axis=1)
            return self._constructor(new_data)
        elif isinstance(item, np.ndarray):
            if len(item) != len(self.index):
                raise ValueError('Item wrong length %d instead of %d!' %
                                 (len(item), len(self.index)))

            # also raises Exception if object array with NA values
            if _is_bool_indexer(item):
                item = np.asarray(item, dtype=bool)

            new_index = self.index[item]
            return self.reindex(new_index)
        else:
            values = self._data.get(item)
            return Series(values, index=self.index)

    def __setitem__(self, key, value):
        """
        Add series to DataFrame in specified column.

        If series is a numpy-array (not a Series/TimeSeries), it must be the
        same length as the DataFrame's index or an error will be thrown.

        Series/TimeSeries will be conformed to the DataFrame's index to
        ensure homogeneity.
        """
        # Array
        if isinstance(key, DataFrame):
            if not (key.index.equals(self.index) and
                    key.columns.equals(self.columns)):
                raise PandasError('Can only index with like-indexed '
                                  'DataFrame objects')

            self._boolean_set(key, value)
        else:
            self._set_item(key, value)

    def _boolean_set(self, key, value):
        mask = key.values
        if mask.dtype != np.bool_:
            raise ValueError('Must pass DataFrame with boolean values only')

        if self._data.is_mixed_dtype():
            raise ValueError('Boolean setting not possible on mixed-type frame')

        self.values[mask] = value

    def insert(self, loc, column, value):
        """
        Insert column into DataFrame at specified location. Raises Exception if
        column is already contained in the DataFrame

        Parameters
        ----------
        loc : int
            Must have 0 <= loc <= len(columns)
        column : object
        value : int, Series, or array-like
        """
        value = self._sanitize_column(value)
        value = np.atleast_2d(value) # is this a hack?
        self._data.insert(loc, column, value)

    def _set_item(self, key, value):
        """
        Add series to DataFrame in specified column.

        If series is a numpy-array (not a Series/TimeSeries), it must be the
        same length as the DataFrame's index or an error will be thrown.

        Series/TimeSeries will be conformed to the DataFrame's index to
        ensure homogeneity.
        """
        value = self._sanitize_column(value)
        value = np.atleast_2d(value) # is this a hack?
        self._data.set(key, value)

    def _sanitize_column(self, value):
        # Need to make sure new columns (which go into the BlockManager as new
        # blocks) are always copied
        if hasattr(value, '__iter__'):
            if isinstance(value, Series):
                if value.index.equals(self.index):
                    # copy the values
                    value = value.values.copy()
                else:
                    value = value.reindex(self.index).values
            else:
                assert(len(value) == len(self.index))

                if not isinstance(value, np.ndarray):
                    value = np.array(value)
                    if value.dtype.type == np.str_:
                        value = np.array(value, dtype=object)
                else:
                    value = value.copy()
        else:
            value = np.repeat(value, len(self.index))

        return value

    def __delitem__(self, key):
        """
        Delete column from DataFrame
        """
        self._data.delete(key)

    def pop(self, item):
        """
        Return column and drop from frame. Raise KeyError if not
        found.

        Returns
        -------
        Series
        """
        result = self[item]
        del self[item]
        return result

    # to support old APIs
    @property
    def _series(self):
        return self._data.get_series_dict()

    def xs(self, key, copy=True):
        """
        Returns a row from the DataFrame as a Series object.

        Parameters
        ----------
        key : some index contained in the index

        Returns
        -------
        xs : Series
        """
        if key not in self.index:
            raise Exception('No cross-section for %s' % key)

        self._consolidate_inplace()
        values = self._data.xs(key, axis=1, copy=copy)
        return Series(values.as_matrix(), index=self.columns)

    #----------------------------------------------------------------------
    # Reindexing

    def reindex(self, index=None, columns=None, method=None):
        """
        Reindex data inside, optionally filling according to some rule.

        Parameters
        ----------
        index : array-like, optional
            preferably an Index object (to avoid duplicating data)
        columns : array-like, optional
        method : {'backfill', 'bfill', 'pad', 'ffill', None}, default None
            Method to use for filling holes in reindexed Series

            pad / ffill: propagate last valid observation forward to next valid
            backfill / bfill: use NEXT valid observation to fill gap

        Returns
        -------
        y : same type as calling instance
        """
        self._consolidate_inplace()
        frame = self

        if index is not None:
            index = _ensure_index(index)
            frame = frame._reindex_index(index, method)

        if columns is not None:
            columns = _ensure_index(columns)
            frame = frame._reindex_columns(columns)

        return frame

    def _reindex_index(self, new_index, method):
        if new_index is self.index:
            return self.copy()
        new_data = self._data.reindex_axis(new_index, method, axis=1)
        return self._constructor(new_data)

    def _reindex_columns(self, new_columns):
        new_data = self._data.reindex_axis(new_columns, axis=0)
        return self._constructor(new_data)

    def reindex_like(self, other, method=None):
        """
        Reindex DataFrame to match indices of another DataFrame

        Parameters
        ----------
        other : DataFrame
        method : string or None

        Notes
        -----
        Like calling s.reindex(index=other.index, columns=other.columns)

        Returns
        -------
        reindexed : DataFrame
        """
        # todo: object columns
        return self.reindex(index=other.index, columns=other.columns,
                            method=method)

    #----------------------------------------------------------------------
    # Reindex-based selection methods

    def filter(self, items=None, like=None, regex=None):
        """
        Restrict frame's columns to set of items or wildcard

        Parameters
        ----------
        items : list-like
            List of columns to restrict to (must not all be present)
        like : string
            Keep columns where "arg in col == True"
        regex : string (regular expression)
            Keep columns with re.search(regex, col) == True

        Notes
        -----
        Arguments are mutually exclusive!

        Returns
        -------
        DataFrame with filtered columns
        """
        import re
        if items is not None:
            return self.reindex(columns=[r for r in items if r in self])
        elif like:
            return self.select(lambda x: like in x, axis=1)
        elif regex:
            matcher = re.compile(regex)
            return self.select(lambda x: matcher.match(x) is not None, axis=1)
        else:
            raise ValueError('items was None!')

    def dropna(self, axis=0, how='any', thresh=None, subset=None):
        """
        Return object with labels on given axis omitted where alternately any or
        all of the data are missing

        Parameters
        ----------
        axis : int
        how : {'any', 'all'}
            any : if any NA values are present, drop that label
            all : if all values are NA, drop that label
        thresh : int, default None
            int value : require that many non-NA values
        subset : array-like

        Returns
        -------
        dropped : DataFrame
        """
        axis_name = self._get_axis_name(axis)

        if axis == 0:
            agg_axis = 1
        elif axis == 1:
            agg_axis = 0
        else: # pragma: no cover
            raise ValueError('axis must be 0 or 1')

        agg_obj = self
        if subset is not None:
            agg_axis_name = self._get_axis_name(agg_axis)
            agg_obj = self.reindex(**{agg_axis_name : subset})

        count = agg_obj.count(axis=agg_axis)

        if thresh is not None:
            mask = count >= thresh
        elif how == 'any':
            mask = count == len(agg_obj._get_axis(agg_axis))
        elif how == 'all':
            mask = count > 0
        else:
            if how is not None:
                raise ValueError('do not recognize %s' % how)
            else:
                raise ValueError('must specify how or thresh')

        labels = self._get_axis(axis)
        new_labels = labels[mask]
        return self.reindex(**{axis_name : new_labels})

    #----------------------------------------------------------------------
    # Sorting

    def sort(self, column=None, axis=0, ascending=True):
        """
        Sort DataFrame either by index (default) by the values in a column

        Parameters
        ----------
        columns : object
            Column name in frame
        ascending : boolean, default True
            Sort ascending vs. descending

        Returns
        -------
        sorted : DataFrame
        """
        if column:
            series = self[column].order(na_last=False)
            sort_index = series.index
        else:
            index = np.asarray(self.index)
            argsorted = np.argsort(index)
            sort_index = index[argsorted.astype(int)]

        if not ascending:
            sort_index = sort_index[::-1]

        return self.reindex(sort_index)

    def sortlevel(self, level=0, axis=0, ascending=True):
        """
        Sort multilevel index by chosen axis and primary level. Data will be
        lexicographically sorted by the chosen level followed by the other
        levels (in order)

        Parameters
        ----------
        level : int
        axis : int
        ascending : bool, default True

        Returns
        -------
        sorted : DataFrame
        """
        the_axis = self._get_axis(axis)
        if not isinstance(the_axis, MultiIndex):
            raise Exception('can only sort by level with a hierarchical index')

        new_axis, indexer = the_axis.sortlevel(level, ascending=ascending)
        new_values = self.values.take(indexer, axis=0)

        if axis == 0:
            index = new_axis
            columns = self.columns
        else:
            index = self.index
            columns = new_axis

        return self._constructor(new_values, index=index, columns=columns)

    #----------------------------------------------------------------------
    # Filling NA's

    def fillna(self, value=None, method='pad'):
        """
        Fill nan values using the specified method.

        Member Series / TimeSeries are filled separately.

        Parameters
        ----------
        method : {'backfill', 'bfill', 'pad', 'ffill', None}, default 'pad'
            Method to use for filling holes in reindexed Series

            pad / ffill: propagate last valid observation forward to next valid
            backfill / bfill: use NEXT valid observation to fill gap

        value : any kind (should be same type as array)
            Value to use to fill holes (e.g. 0)

        Returns
        -------
        y : DataFrame

        See also
        --------
        DataFrame.reindex, DataFrame.asfreq
        """
        if value is None:
            result = {}
            series = self._series
            for col, s in series.iteritems():
                result[col] = s.fillna(method=method, value=value)
            return self._constructor(result, index=self.index,
                                     columns=self.columns)
        else:
            # Float type values
            if len(self.columns) == 0:
                return self

            new_data = self._data.fillna(value)
            return self._constructor(new_data, index=self.index,
                                     columns=self.columns)

    #----------------------------------------------------------------------
    # Rename

    def rename(self, index=None, columns=None):
        """
        Alter index and / or columns using input function or functions

        Parameters
        ----------
        index : dict-like or function, optional
            Transformation to apply to index values
        columns : dict-like or function, optional
            Transformation to apply to column values

        See also
        --------
        Series.rename

        Notes
        -----
        Function / dict values must be unique (1-to-1)

        Returns
        -------
        y : DataFrame (new object)
        """
        if isinstance(index, (dict, Series)):
            index = index.__getitem__

        if isinstance(columns, (dict, Series)):
            columns = columns.__getitem__

        if index is None and columns is None:
            raise Exception('must pass either index or columns')

        self._consolidate_inplace()

        result = self.copy()

        if index is not None:
            result._rename_index_inplace(index)

        if columns is not None:
            result._rename_columns_inplace(columns)

        return result

    def _rename_index_inplace(self, mapper):
        self._data = self._data.rename_axis(mapper, axis=1)

    def _rename_columns_inplace(self, mapper):
        self._data = self._data.rename_items(mapper)

    #----------------------------------------------------------------------
    # Arithmetic / combination related

    def _combine_frame(self, other, func, fill_value=None):
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
        new_index = _union_indices(self.index, other.index)

        # some shortcuts
        if fill_value is None:
            if not self and not other:
                return self._constructor(index=new_index)
            elif not self:
                return other * nan
            elif not other:
                return self * nan

        need_reindex = False
        new_columns = _union_indices(self.columns, other.columns)
        need_reindex = (need_reindex or not new_index.equals(self.index)
                        or not new_index.equals(other.index))
        need_reindex = (need_reindex or not new_columns.equals(self.columns)
                        or not new_columns.equals(other.columns))

        this = self
        if need_reindex:
            this = self.reindex(index=new_index, columns=new_columns)
            other = other.reindex(index=new_index, columns=new_columns)

        this_vals = this.values
        other_vals = other.values

        if fill_value is not None:
            this_mask = isnull(this_vals)
            other_mask = isnull(other_vals)
            this_vals = this_vals.copy()
            other_vals = other_vals.copy()

            # one but not both
            mask = this_mask ^ other_mask
            this_vals[this_mask & mask] = fill_value
            other_vals[other_mask & mask] = fill_value

        result = func(this_vals, other_vals)
        return self._constructor(result, index=new_index, columns=new_columns,
                                 copy=False)

    def _indexed_same(self, other):
        same_index = self.index.equals(other.index)
        same_columns = self.columns.equals(other.columns)
        return same_index and same_columns

    def _combine_series_infer(self, other, func, fill_value=None):
        if len(other) == 0:
            return self * nan

        if len(self) == 0:
            # Ambiguous case, use _series so works with DataFrame
            return self._constructor(data=self._series, index=self.index,
                                     columns=self.columns)

        # teeny hack because one does DataFrame + TimeSeries all the time
        if self.index.is_all_dates() and other.index.is_all_dates():
            return self._combine_match_index(other, func, fill_value)
        else:
            return self._combine_match_columns(other, func, fill_value)

    def _combine_match_index(self, other, func, fill_value=None):
        new_index = _union_indices(self.index, other.index)
        values = self.values
        other_vals = other.values

        # Operate row-wise
        if not other.index.equals(new_index):
            other_vals = other.reindex(new_index).values

        if not self.index.equals(new_index):
            values = self.reindex(new_index).values

        if fill_value is not None:
            raise NotImplementedError

        return self._constructor(func(values.T, other_vals).T, index=new_index,
                                 columns=self.columns, copy=False)

    def _combine_match_columns(self, other, func, fill_value=None):
        newCols = self.columns.union(other.index)

        # Operate column-wise
        this = self.reindex(columns=newCols)
        other = other.reindex(newCols).values

        if fill_value is not None:
            raise NotImplementedError

        return self._constructor(func(this.values, other), index=self.index,
                                 columns=newCols, copy=False)

    def _combine_const(self, other, func):
        if not self:
            return self

        return self._constructor(func(self.values, other), index=self.index,
                                 columns=self.columns, copy=False)

    def _compare_frame(self, other, func):
        if not self._indexed_same(other):
            raise Exception('Can only compare identically-labeled '
                            'DataFrame objects')

        new_data = {}
        for col in self.columns:
            new_data[col] = func(self[col], other[col])

        return self._constructor(data=new_data, index=self.index,
                                 columns=self.columns, copy=False)

    def combine(self, other, func, fill_value=None):
        """
        Add two DataFrame objects and do not propagate NaN values, so if for a
        (column, time) one frame is missing a value, it will default to the
        other frame's value (which might be NaN as well)

        Parameters
        ----------
        other : DataFrame

        Returns
        -------
        DataFrame
        """
        if not other:
            return self.copy()

        if not self:
            return other.copy()

        new_index = self.index
        this = self

        if not self.index.equals(other.index):
            new_index = self.index + other.index
            this = self.reindex(new_index)
            other = other.reindex(new_index)

        new_columns = _try_sort(set(this.columns + other.columns))
        do_fill = fill_value is not None

        result = {}
        for col in new_columns:
            if col in this and col in other:
                series = this[col].values
                otherSeries = other[col].values

                if do_fill:
                    this_mask = isnull(series)
                    other_mask = isnull(otherSeries)
                    series = series.copy()
                    otherSeries = otherSeries.copy()
                    series[this_mask] = fill_value
                    otherSeries[other_mask] = fill_value

                arr = func(series, otherSeries)

                if do_fill:
                    arr = common.ensure_float(arr)
                    arr[this_mask & other_mask] = nan

                result[col] = arr

            elif col in this:
                result[col] = this[col]
            elif col in other:
                result[col] = other[col]

        return self._constructor(result, index=new_index, columns=new_columns)

    def combineFirst(self, other):
        """
        Combine two DataFrame objects and default to value
        in frame calling the method.

        Parameters
        ----------
        otherFrame : DataFrame

        Examples
        --------
        a.combineFirst(b)
            a's values prioritized, use values from b to fill holes

        Returns
        -------
        DataFrame
        """
        combiner = lambda x, y: np.where(isnull(x), y, x)
        return self.combine(other, combiner)

    #----------------------------------------------------------------------
    # Misc methods

    def first_valid_index(self):
        return self.index[self.count(1) > 0][0]

    def last_valid_index(self):
        return self.index[self.count(1) > 0][-1]

    def head(self):
        return self[:5]

    def tail(self):
        return self[-5:]

    #----------------------------------------------------------------------
    # Data reshaping

    def pivot(self, index=None, columns=None, values=None):
        """
        Produce 'pivot' table this DataFrame. Uses unique values from index /
        columns to form axes and return either DataFrame or WidePanel, depending
        on whether you request a single value column (DataFrame) or all columns
        (WidePanel)

        Parameters
        ----------
        index : string or object
            Column name to use to make new frame's index
        columns : string or object
            Column name to use to make new frame's columns
        values : string or object, optional
            Column name to use for populating new frame's values

        Examples
        --------
        >>> df
            foo   bar  baz
        0   one   A    1.
        1   one   B    2.
        2   one   C    3.
        3   two   A    4.
        4   two   B    5.
        5   two   C    6.

        >>> df.pivot('foo', 'bar', 'baz')
             A   B   C
        one  1   2   3
        two  4   5   6

        >>> df.pivot('foo', 'bar')['baz']
             A   B   C
        one  1   2   3
        two  4   5   6

        Returns
        -------
        pivoted : DataFrame (value column specified) or WidePanel (no value
        column specified)
        """
        from pandas.core.panel import _make_long_index, LongPanel

        index_vals = self[index]
        column_vals = self[columns]
        long_index = _make_long_index(index_vals, column_vals)

        if values is None:
            items = self.columns - [index, columns]
            mat = self.reindex(columns=items).values
        else:
            items = [values]
            mat = np.atleast_2d(self[values].values).T

        lp = LongPanel(mat, index=long_index, columns=items)
        lp = lp.sortlevel(level=0)

        wp = lp.to_wide()
        if values is not None:
            return wp[values]
        else:
            return wp

    def stack(self):
        """
        Convert DataFrame to Series with multi-level Index

        Returns
        -------
        stacked : Series
        """
        N, K = len(self.index), len(self.columns)
        ilabels = np.arange(N).repeat(K)
        clabels = np.tile(np.arange(K), N).ravel()
        index = MultiIndex(levels=[self.index, self.columns],
                           labels=[ilabels, clabels])
        return Series(self.values.ravel(), index=index)

    def delevel(self):
        """
        For DataFrame with multi-level index, return new DataFrame with labeling
        information in the columns under names 'level_0', 'level_1', etc.

        Note: experimental, subject ot API change

        Returns
        -------
        deleveled : DataFrame
        """
        if not isinstance(self.index, MultiIndex):
            raise Exception('this DataFrame does not have a multi-level index')

        new_obj = self.copy()

        zipped = zip(self.index.levels, self.index.labels)
        for i, (lev, lab) in reversed(list(enumerate(zipped))):
            new_obj.insert(0, 'label_%d' % i, np.asarray(lev).take(lab))

        new_obj.index = np.arange(len(new_obj))

        return new_obj

    #----------------------------------------------------------------------
    # Time series-related

    def asfreq(self, freq, method=None):
        """
        Convert all TimeSeries inside to specified frequency using
        DateOffset objects. Optionally provide fill method to pad or
        backfill missing values.

        Parameters
        ----------
        offset : DateOffset object, or string in {'WEEKDAY', 'EOM'}
            DateOffset object or subclass (e.g. monthEnd)

        method : {'backfill', 'bfill', 'pad', 'ffill', None}
            Method to use for filling holes in reindexed Series

            pad / ffill: propagate last valid observation forward to next valid
            backfill / bfill: use NEXT valid observation to fill methdo
        """
        if len(self.index) == 0:
            return self.copy()

        if isinstance(freq, datetools.DateOffset):
            dateRange = DateRange(self.index[0], self.index[-1], offset=freq)
        else:
            dateRange = DateRange(self.index[0], self.index[-1], timeRule=freq)

        return self.reindex(dateRange, method=method)

    def diff(self, periods=1):
        return self - self.shift(periods)

    def shift(self, periods, offset=None, timeRule=None):
        """
        Shift the underlying series of the DataFrame and Series objects within
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
        DataFrame
        """
        if periods == 0:
            return self

        if timeRule is not None and offset is None:
            offset = datetools.getOffset(timeRule)

        def _shift_block(blk, indexer):
            new_values = blk.values.take(indexer, axis=1)
            # convert integer to float if necessary. need to do a lot more than
            # that, handle boolean etc also
            new_values = common.ensure_float(new_values)
            if periods > 0:
                new_values[:, :periods] = nan
            else:
                new_values[:, periods:] = nan
            return make_block(new_values, blk.items, blk.ref_items)

        if offset is None:
            indexer = self._shift_indexer(periods)
            new_blocks = [_shift_block(b, indexer) for b in self._data.blocks]
            new_data = BlockManager(new_blocks, [self.columns, self.index])
        else:
            new_data = self._data.copy()
            new_data.axes[1] = self.index.shift(periods, offset)

        return self._constructor(new_data)

    def _shift_indexer(self, periods):
        # small reusable utility
        N = len(self)
        indexer = np.zeros(N, dtype=int)

        if periods > 0:
            indexer[periods:] = np.arange(N - periods)
        else:
            indexer[:periods] = np.arange(-periods, N)

        return indexer

    #----------------------------------------------------------------------
    # Function application

    def apply(self, func, axis=0, broadcast=False):
        """
        Applies func to columns (Series) of this DataFrame and returns either
        a DataFrame (if the function produces another series) or a Series
        indexed on the column names of the DataFrame if the function produces
        a value.

        Parameters
        ----------
        func : function
            Function to apply to each column
        axis : {0, 1}
        broadcast : bool, default False
            For aggregation functions, return object of same size with values
            propagated

        Examples
        --------
        >>> df.apply(numpy.sqrt) --> DataFrame
        >>> df.apply(numpy.sum) --> Series

        Notes
        -----
        Functions altering the index are not supported (yet)
        """
        if not len(self.columns):
            return self

        if isinstance(func, np.ufunc):
            results = func(self.values)
            return self._constructor(data=results, index=self.index,
                                     columns=self.columns, copy=False)
        else:
            if not broadcast:
                return self._apply_standard(func, axis)
            else:
                return self._apply_broadcast(func, axis)

    def _apply_standard(self, func, axis):
        if axis == 0:
            target = self
            agg_index = self.columns
        elif axis == 1:
            target = self.T
            agg_index = self.index

        results = {}
        for k in target.columns:
            results[k] = func(target[k])

        if hasattr(results.values()[0], '__iter__'):
            result = self._constructor(data=results, index=target.index,
                                       columns=target.columns)

            if axis == 1:
                result = result.T

            return result
        else:
            return Series(results, index=agg_index)

    def _apply_broadcast(self, func, axis):
        if axis == 0:
            target = self
        elif axis == 1:
            target = self.T

        result_values = np.empty_like(target.values)
        columns = target.columns
        for i, col in enumerate(columns):
            result_values[:, i] = func(target[col])

        result = self._constructor(result_values, index=target.index,
                                   columns=target.columns)

        if axis == 1:
            result = result.T

        return result

    def _apply_level(self, f, axis='major', broadcast=False):
        from pandas.core.panel import LongPanel

        if axis == 'major':
            panel = self.swapaxes()
            result = panel._apply_level(f, axis='minor', broadcast=broadcast)
            if broadcast:
                result = result.swapaxes()

            return result

        bounds = self.index._bounds
        values = self.values
        N, _ = values.shape
        result = group_agg(values, bounds, f)

        if broadcast:
            repeater = np.concatenate((np.diff(bounds), [N - bounds[-1]]))
            panel = LongPanel(result.repeat(repeater, axis=0),
                              columns=self.items, index=self.index)
        else:
            panel = DataFrame(result, index=self.major_axis,
                              columns=self.items)

        return panel

    def tapply(self, func):
        """
        Apply func to the transposed DataFrame, results as per apply
        """
        return self.apply(func, axis=1)

    def applymap(self, func):
        """
        Apply a function to a DataFrame that is intended to operate
        elementwise, i.e. like doing
            map(func, series) for each series in the DataFrame

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

        return self._constructor(results, index=self.index,
                                 columns=self.columns)

    #----------------------------------------------------------------------
    # Merging / joining methods

    def append(self, other):
        """
        Append columns of other to end of this frame's columns and index.

        Columns not in this frame are added as new columns.
        """
        if not other:
            return self.copy()
        if not self:
            return other.copy()

        new_index = np.concatenate((self.index, other.index))
        new_data = {}

        new_columns = self.columns

        if not new_columns.equals(other.columns):
            new_columns = self.columns + other.columns

        for column, series in self.iteritems():
            values = series.values
            if column in other:
                other_values = other[column].values
                new_data[column] = np.concatenate((values, other_values))
            else:
                new_data[column] = series

        for column, series in other.iteritems():
            if column not in self:
                new_data[column] = series

        return self._constructor(data=new_data, index=new_index,
                                 columns=new_columns)

    def join(self, other, on=None, how=None):
        """
        Join columns with other DataFrame either on index or on a key
        column

        Parameters
        ----------
        other : DataFrame
            Index should be similar to one of the columns in this one
        on : string, default None
            Column name to use, otherwise join on index
        how : {'left', 'right', 'outer', 'inner'}
            How to handle indexes of the two objects. Default: 'left'
            for joining on index, None otherwise
            * left: use calling frame's index
            * right: use input frame's index
            * outer: form union of indexes
            * inner: use intersection of indexes

        Returns
        -------
        joined : DataFrame
        """
        if on is not None:
            if how is not None:
                raise Exception('how parameter is not valid when '
                                '*on* specified')
            return self._join_on(other, on)
        else:
            if how is None:
                how = 'left'
            return self._join_index(other, how)

    def _join_on(self, other, on):
        if len(other.index) == 0:
            return self

        if on not in self:
            raise Exception('%s column not contained in this frame!' % on)

        new_data = self._data.join_on(other._data, self[on], axis=1)
        return self._constructor(new_data)

    def _join_index(self, other, how):
        join_index = self._get_join_index(other, how)
        this_data = self.reindex(join_index)._data
        other_data = other.reindex(join_index)._data

        # merge blocks
        merged_data = this_data.merge(other_data)
        assert(merged_data.axes[1] is join_index) # maybe unnecessary
        return self._constructor(merged_data)

    def _get_join_index(self, other, how):
        if how == 'left':
            join_index = self.index
        elif how == 'right':
            join_index = other.index
        elif how == 'inner':
            join_index = self.index.intersection(other.index)
        elif how == 'outer':
            join_index = self.index.union(other.index)
        else:
            raise Exception('do not recognize join method %s' % how)

        return join_index

    #----------------------------------------------------------------------
    # groupby

    def groupby(self, by, axis=0, column=None):
        """
        Goup series using mapper (dict or key function, apply given function to
        group, return result as series) or by a series of columns

        Parameters
        ----------
        by : mapping function, dict, Series, or tuple / list of column names
            Called on each element of the object index to determine
            the groups.  If a dict or Series is passed, the Series or
            dict VALUES will be used to determine the groups
        axis : int, default 0

        Examples
        --------
        # DataFrame result
        >>> data.groupby(func, axis=0).mean()

        # DataFrame result
        >>> data.groupby(['col1', 'col2'])['col3'].mean()

        # WidePanel result
        >>> data.groupby(['col1', 'col2']).mean()

        Returns
        -------
        GroupBy object
        """
        from pandas.core.groupby import groupby
        return groupby(self, by, axis=axis)

    def tgroupby(self, keyfunc, applyfunc):
        """
        Aggregate columns based on passed function

        Parameters
        ----------
        keyfunc : function
        applyfunc : function

        Returns
        -------
        y : DataFrame
        """
        return self.T.groupby(keyfunc).aggregate(applyfunc).T

    #----------------------------------------------------------------------
    # Statistical methods, etc.

    def corr(self):
        """
        Compute pairwise correlation of columns, excluding NaN values

        Returns
        -------
        y : DataFrame
        """
        cols = self.columns
        mat = self.as_matrix(cols).T
        baseCov = np.cov(mat)

        sigma = np.sqrt(np.diag(baseCov))
        correl = baseCov / np.outer(sigma, sigma)

        # Get the covariance with items that have NaN values
        for i, A in enumerate(mat):
            aok = np.isfinite(A)
            if not aok.all():
                for j, B in enumerate(mat):
                    commonVec = aok & np.isfinite(B)
                    if commonVec.any():
                        ac, bc = A[commonVec], B[commonVec]
                        c = np.corrcoef(ac, bc)[0, 1]
                        correl[i, j] = c
                        correl[j, i] = c

        return self._constructor(correl, index=cols, columns=cols)

    def corrwith(self, other, axis=0, drop=False):
        """
        Compute pairwise correlation between rows or columns of two DataFrame
        objects.

        Parameters
        ----------
        other : DataFrame
        axis : int
        drop : boolean, default False
            Drop missing indices from result, default returns union of all

        Returns
        -------
        correls : Series
        """
        com_index = self._intersect_index(other)
        com_cols = self._intersect_columns(other)

        # feels hackish
        if axis == 0:
            result_index = com_index
            if not drop:
                result_index = self.columns.union(other.columns)
        else:
            result_index = com_cols
            if not drop:
                result_index = self.index.union(other.index)

        left = self.reindex(index=com_index, columns=com_cols)
        right = other.reindex(index=com_index, columns=com_cols)

        # mask missing values
        left = left + right * 0
        right = right + left * 0

        # demeaned data
        ldem = left - left.mean(axis)
        rdem = right - right.mean(axis)

        num = (ldem * rdem).sum(axis)
        dom = (left.count(axis) - 1) * left.std(axis) * right.std(axis)

        correl = num / dom

        if not drop:
            correl = correl.reindex(result_index)

        return correl

    def describe(self):
        """
        Generate various summary statistics of columns, excluding NaN values

        Returns
        -------
        DataFrame
        """
        cols = self._get_numeric_columns()
        tmp = self.reindex(columns=cols)

        cols_destat = ['count', 'mean', 'std', 'min',
                       '10%', '50%', '90%', 'max']

        data = [tmp.count(), tmp.mean(), tmp.std(), tmp.min(),
                tmp.quantile(.1), tmp.median(),
                tmp.quantile(.9), tmp.max()]

        return self._constructor(data, index=cols_destat, columns=cols)

    #----------------------------------------------------------------------
    # ndarray-like stats methods

    def count(self, axis=0, numeric_only=False):
        """
        Return array or Series of # observations over requested axis.

        Parameters
        ----------
        axis : {0, 1}
            0 for row-wise, 1 for column-wise
        numeric_only : boolean, default False
            Include only float, int, boolean data

        Notes
        -----
        Also examines non-float data and checks for None and NaN in such data

        Returns
        -------
        Series or TimeSeries
        """
        try:
            y, axis_labels = self._get_agg_data(axis, numeric_only=numeric_only)
            mask = np.empty(y.shape, dtype=bool)
            mask.flat = notnull(y.ravel())
            return Series(mask.sum(axis), index=axis_labels)
        except Exception:
            f = lambda s: notnull(s).sum()
            return self.apply(f, axis=axis)

    def sum(self, axis=0, numeric_only=False):
        """
        Return array or Series of sums over requested axis.

        Parameters
        ----------
        axis : {0, 1}
            0 for row-wise, 1 for column-wise
        numeric_only : boolean, default False
            Include only float, int, boolean data

        Returns
        -------
        Series or TimeSeries

        Examples
        --------
        >>> df
            c1  c2
        a   1   0
        b   0   2
        c   3   0
        d   0   4

        >>> df.sum(axis=0)
        c1    4
        c2    6
        """
        y, axis_labels = self._get_agg_data(axis, numeric_only=numeric_only)

        if len(axis_labels) == 0:
            return Series([], index=[])

        if y.dtype == np.object_:
            the_sum = y.sum(axis)
        else:
            mask = np.isfinite(y)
            if not issubclass(y.dtype.type, np.int_):
                y[-mask] = 0
            the_sum = y.sum(axis)
            the_count = mask.sum(axis)

            ct_mask = the_count == 0
            if ct_mask.any():
                the_sum[ct_mask] = nan

        return Series(the_sum, index=axis_labels)

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

    def product(self, axis=0):
        """
        Return array or Series of products over requested axis.

        Parameters
        ----------
        axis : {0, 1}
            0 for row-wise, 1 for column-wise

        Returns
        -------
        Series or TimeSeries
        """
        y = np.array(self.values, subok=True)
        try:

            if not issubclass(y.dtype.type, np.int_):
                y[np.isnan(y)] = 1
            theProd = y.prod(axis)
            theCount = self.count(axis)
            theProd[theCount == 0] = nan
        except Exception:
            theProd = self.apply(np.prod, axis=axis)

        return Series(theProd, index=self._get_agg_axis(axis))

    def mean(self, axis=0):
        """
        Return array or Series of means over requested axis.

        Parameters
        ----------
        axis : {0, 1}
            0 for row-wise, 1 for column-wise

        Notes
        -----

        Returns
        -------
        Series or TimeSeries
        """
        summed = self.sum(axis, numeric_only=True)
        count = self.count(axis, numeric_only=True).astype(float)
        return summed / count

    def quantile(self, q=0.5, axis=0):
        """
        Return array or Series of values at the given quantile over requested
        axis.

        Parameters
        ----------
        q : quantile
            0 <= q <= 1

        axis : {0, 1}
            0 for row-wise, 1 for column-wise

        Returns
        -------
        quantiles : Series
        """
        from scipy.stats import scoreatpercentile
        per = q * 100
        def f(arr):
            arr = arr.values
            if arr.dtype != np.float_:
                arr = arr.astype(float)
            arr = arr[notnull(arr)]
            if len(arr) == 0:
                return nan
            else:
                return scoreatpercentile(arr, per)

        return self.apply(f, axis=axis)

    def median(self, axis=0):
        """
        Return array or Series of medians over requested axis.

        Parameters
        ----------
        axis : {0, 1}
            0 for row-wise, 1 for column-wise

        Returns
        -------
        Series or TimeSeries
        """
        def f(arr):
            if arr.dtype != np.float_:
                arr = arr.astype(float)
            return _tseries.median(arr[notnull(arr)])

        if axis == 0:
            med = [f(self[col].values) for col in self.columns]
            return Series(med, index=self.columns)
        elif axis == 1:
            med = [f(self.xs(k).values) for k in self.index]
            return Series(med, index=self.index)
        else:
            raise Exception('Must have 0<= axis <= 1')

    def mad(self, axis=0):
        """
        Return array or Series of mean absolute deviation over
        requested axis.

        Parameters
        ----------
        axis : {0, 1}
            0 for row-wise, 1 for column-wise

        Returns
        -------
        Series or TimeSeries
        """
        if axis == 0:
            demeaned = self-self.mean(axis=axis)
        else:
            demeaned = (self.T-self.mean(axis=axis)).T

        y = np.array(demeaned.values, subok=True)

        if not issubclass(y.dtype.type, np.int_):
            y[np.isnan(y)] = 0

        result = np.abs(y).mean(axis=axis)

        if axis == 0:
            return Series(result, demeaned.columns)
        else:
            return Series(result, demeaned.index)

    def var(self, axis=0):
        """
        Return array or Series of unbiased variance over requested axis.

        Parameters
        ----------
        axis : {0, 1}
            0 for row-wise, 1 for column-wise

        Returns
        -------
        Series or TimeSeries
        """
        y, axis_labels = self._get_agg_data(axis, numeric_only=True)

        mask = np.isnan(y)
        count = (y.shape[axis] - mask.sum(axis)).astype(float)
        y[mask] = 0

        X = y.sum(axis)
        XX = (y**2).sum(axis)

        theVar = (XX - X**2 / count) / (count - 1)

        return Series(theVar, index=axis_labels)

    def std(self, axis=0):
        """
        Return array or Series of unbiased std deviation over requested axis.

        Parameters
        ----------
        axis : {0, 1}
            0 for row-wise, 1 for column-wise

        Returns
        -------
        Series or TimeSeries
        """
        return np.sqrt(self.var(axis=axis))

    def skew(self, axis=0):
        """
        Return array or Series of unbiased skewness over requested axis.

        Parameters
        ----------
        axis : {0, 1}
            0 for row-wise, 1 for column-wise

        Returns
        -------
        Series or TimeSeries
        """
        y, axis_labels = self._get_agg_data(axis, numeric_only=True)

        mask = np.isnan(y)
        count = (y.shape[axis] - mask.sum(axis)).astype(float)
        y[mask] = 0

        A = y.sum(axis) / count
        B = (y**2).sum(axis) / count  - A**2
        C = (y**3).sum(axis) / count - A**3 - 3*A*B

        theSkew = (np.sqrt((count**2-count))*C) / ((count-2)*np.sqrt(B)**3)

        return Series(theSkew, index=axis_labels)

    def _get_agg_data(self, axis, numeric_only=True):
        num_cols = self._get_numeric_columns()

        if len(num_cols) < len(self.columns) and numeric_only:
            y = self.as_matrix(num_cols)
            if axis == 0:
                axis_labels = num_cols
            else:
                axis_labels = self.index
        else:
            y = self.values.copy()
            axis_labels = self._get_agg_axis(axis)

        return y, axis_labels

    def _get_agg_axis(self, axis_num):
        if axis_num == 0:
            return self.columns
        elif axis_num == 1:
            return self.index
        else:
            raise Exception('Must have 0<= axis <= 1')

    def _get_numeric_columns(self):
        from pandas.core.internals import ObjectBlock

        cols = []
        for col, blk in zip(self.columns, self._data.block_id_vector):
            if not isinstance(self._data.blocks[blk], ObjectBlock):
                cols.append(col)

        return cols

    def clip(self, upper=None, lower=None):
        """
        Trim values at input threshold(s)

        Parameters
        ----------
        lower : float, default None
        upper : float, default None

        Returns
        -------
        y : DataFrame
        """
        return self.apply(lambda x: x.clip(lower=lower, upper=upper))

    def clip_upper(self, threshold):
        """
        Trim values above threshold

        Returns
        -------
        y : DataFrame
        """
        return self.apply(lambda x: x.clip_upper(threshold))

    def clip_lower(self, threshold):
        """
        Trim values below threshold

        Returns
        -------
        y : DataFrame
        """
        return self.apply(lambda x: x.clip_lower(threshold))

    #----------------------------------------------------------------------
    # Plotting

    def plot(self, subplots=False, sharex=True,
             sharey=False, use_index=True, **kwds): # pragma: no cover
        """
        Make line plot of DataFrame's series with the index on the x-axis using
        matplotlib / pylab.

        Parameters
        ----------
        subplots : boolean, default False
            Make separate subplots for each time series
        sharex : boolean, default True
            In case subplots=True, share x axis
        sharey : boolean, default False
            In case subplots=True, share y axis
        use_index : boolean, default True
            Use index as ticks for x axis
        kwds : keywords
            Options to pass to Axis.plot

        Notes
        -----
        This method doesn't make much sense for cross-sections,
        and will error.
        """
        import matplotlib.pyplot as plt

        if subplots:
            _, axes = plt.subplots(nrows=len(self.columns),
                                   sharex=sharex, sharey=sharey)
        else:
            fig = plt.figure()
            ax = fig.add_subplot(111)

        if use_index:
            x = self.index
        else:
            x = range(len(self))

        for i, col in enumerate(_try_sort(self.columns)):
            if subplots:
                ax = axes[i]
                ax.plot(x, self[col].values, 'k', label=col, **kwds)
                ax.legend(loc='best')
            else:
                ax.plot(x, self[col].values, label=col, **kwds)

    def hist(self): # pragma: no cover
        """
        Draw Histogram the DataFrame's series using matplotlib / pylab.

        Parameters
        ----------
        kwds : other plotting keyword arguments

        """
        import matplotlib.pyplot as plt

        n = len(self.columns)
        k = 1
        while k**2 < n:
            k += 1
        _, axes = plt.subplots(nrows=k, ncols=k)

        for i, col in enumerate(_try_sort(self.columns)):
            ax = axes[i / k][i % k]
            ax.hist(self[col].values)
            ax.set_title(col)

    #----------------------------------------------------------------------
    # Deprecated stuff

    def combineAdd(self, other):
        """
        Add two DataFrame objects and do not propagate
        NaN values, so if for a (column, time) one frame is missing a
        value, it will default to the other frame's value (which might
        be NaN as well)

        Parameters
        ----------
        other : DataFrame

        Returns
        -------
        DataFrame
        """
        return self.add(other, fill_value=0.)

    def combineMult(self, other):
        """
        Multiply two DataFrame objects and do not propagate NaN values, so if
        for a (column, time) one frame is missing a value, it will default to
        the other frame's value (which might be NaN as well)

        Parameters
        ----------
        other : DataFrame

        Returns
        -------
        DataFrame
        """
        return self.mul(other, fill_value=1.)

    def toDataMatrix(self): # pragma: no cover
        warnings.warn("toDataMatrix will disappear in next release "
                      "as there is no longer a DataMatrix class", FutureWarning)
        return self.copy()

    def rows(self): # pragma: no cover
        """Alias for the frame's index"""
        warnings.warn("Replace usage of .rows() with .index, will be removed "
                      "in next release", FutureWarning)
        return self.index

    def cols(self): # pragma: no cover
        """Return sorted list of frame's columns"""
        warnings.warn("Replace usage of .cols() with .columns, will be removed "
                      "in next release", FutureWarning)
        return list(self.columns)

    def getXS(self, key): # pragma: no cover
        warnings.warn("'getXS' is deprecated. Use 'xs' instead",
                      FutureWarning)
        return self.xs(key)

    def merge(self, *args, **kwargs): # pragma: no cover
        warnings.warn("merge is deprecated. Use 'join' instead",
                      FutureWarning)
        return self.join(*args, **kwargs)

    def asMatrix(self, *args, **kwargs): # pragma: no cover
        warnings.warn("asMatrix is deprecated. Use 'as_matrix' or .values "
                      "instead", FutureWarning)
        return self.as_matrix(*args, **kwargs)

    def toRecords(self, *args, **kwargs): # pragma: no cover
        warnings.warn("toRecords is deprecated. Use 'to_records' "
                      "instead", FutureWarning)
        return self.to_records(*args, **kwargs)

    @classmethod
    def fromRecords(cls, *args, **kwargs): # pragma: no cover
        warnings.warn("fromRecords is deprecated. Use 'from_records' "
                      "instead", FutureWarning)
        return cls.from_records(*args, **kwargs)

    def _firstTimeWithValue(self): # pragma: no cover
        warnings.warn("_firstTimeWithValue is deprecated. Use "
                      "first_valid_index instead", FutureWarning)
        return self.first_valid_index()

    def _lastTimeWithValue(self): # pragma: no cover
        warnings.warn("_firstTimeWithValue is deprecated. Use "
                      "last_valid_index instead", FutureWarning)
        return self.last_valid_index()

    def dropEmptyRows(self, specificColumns=None): # pragma: no cover
        """
        Return DataFrame with rows omitted containing ALL NaN values
        for optionally specified set of columns.

        Parameters
        ----------
        specificColumns : list-like, optional keyword
            Columns to consider in removing NaN values. As a typical
            application, you might provide the list of the columns involved in
            a regression to exlude all the missing data in one shot.

        Returns
        -------
        This DataFrame with rows containing any NaN values deleted
        """
        warnings.warn("dropEmptyRows is deprecated. Use dropna with how='all'",
                      FutureWarning)
        return self.dropna(axis=0, subset=specificColumns, how='all')

    def dropIncompleteRows(self, specificColumns=None,
                           minObs=None): # pragma: no cover
        """
        Return DataFrame with rows omitted containing ANY NaN values for
        optionally specified set of columns.

        Parameters
        ----------
        minObs : int or None (default)
           Instead of requiring all the columns to have observations, require
           only minObs observations
        specificColumns : list-like, optional keyword
            Columns to consider in removing NaN values. As a typical
            application, you might provide the list of the columns involved in
            a regression to exlude all the missing data in one shot.

        Returns
        -------
        This DataFrame with rows containing any NaN values deleted

        """
        warnings.warn("dropEmptyRows is deprecated. Use dropna",
                      FutureWarning)
        if minObs is None:
            return self.dropna(axis=0, subset=specificColumns, how='any')
        else:
            return self.dropna(axis=0, subset=specificColumns, thresh=minObs)

    #----------------------------------------------------------------------
    # Fancy indexing

    _ix = None
    @property
    def ix(self):
        from pandas.core.indexing import _DataFrameIndexer
        if self._ix is None:
            self._ix = _DataFrameIndexer(self)

        return self._ix


def group_agg(values, bounds, f):
    """
    R-style aggregator

    Parameters
    ----------
    values : N-length or N x K ndarray
    bounds : B-length ndarray
    f : ndarray aggregation function

    Returns
    -------
    ndarray with same length as bounds array
    """
    if values.ndim == 1:
        N = len(values)
        result = np.empty(len(bounds), dtype=float)
    elif values.ndim == 2:
        N, K = values.shape
        result = np.empty((len(bounds), K), dtype=float)

    testagg = f(values[:min(1, len(values))])
    if isinstance(testagg, np.ndarray) and testagg.ndim == 2:
        raise Exception('Passed function does not aggregate!')

    for i, left_bound in enumerate(bounds):
        if i == len(bounds) - 1:
            right_bound = N
        else:
            right_bound = bounds[i + 1]

        result[i] = f(values[left_bound : right_bound])

    return result


def factor_agg(factor, vec, func):
    """
    Aggregate array based on Factor

    Parameters
    ----------
    factor : Factor
        length n
    vec : sequence
        length n
    func : function
        1D array aggregation function

    Returns
    -------
    ndarray corresponding to Factor levels
    """
    indexer = np.argsort(factor.labels)
    unique_labels = np.arange(len(factor.levels))

    ordered_labels = factor.labels.take(indexer)
    ordered_vec = np.asarray(vec).take(indexer)
    bounds = ordered_labels.searchsorted(unique_labels)

    return group_agg(ordered_vec, bounds, func)


def _union_indices(a, b):
    if len(a) == 0:
        return b
    elif len(b) == 0:
        return a
    if not a.equals(b):
        return a.union(b)
    return a

def extract_index(data):
    def _union_if(index, new_index):
        if index is None:
            index = new_index
        else:
            index = index.union(new_index)
        return index

    def _get_index(obj):
        if isinstance(v, Series):
            return v.index
        elif isinstance(v, dict):
            return Index(_try_sort(v))

    index = None
    if len(data) == 0:
        index = NULL_INDEX
    elif len(data) > 0 and index is None:
        have_raw_arrays = _check_data_types(data)

        # this is still kludgier than I'd like
        if have_raw_arrays:
            lengths = list(set(len(x) for x in data.values()))
            if len(lengths) > 1:
                raise ValueError('arrays must all be same length')
            index = Index(np.arange(lengths[0]))
        else:
            for v in data.values():
                index = _union_if(index, _get_index(v))

    if len(index) == 0:
        index = NULL_INDEX

    return _ensure_index(index)

def _check_data_types(data):
    have_raw_arrays = False
    have_series = False
    for v in data.values():
        if not isinstance(v, (dict, Series)):
            have_raw_arrays = True
        else:
            have_series = True

    if have_series and have_raw_arrays:
        raise Exception('Cannot mix Series / dict objects'
                        ' with ndarray / sequence input')

    return have_raw_arrays

def _prep_ndarray(values, copy=True):
    if not isinstance(values, np.ndarray):
        values = np.asarray(values)
        # NumPy strings are a pain, convert to object
        if issubclass(values.dtype.type, basestring):
            values = np.array(values, dtype=object, copy=True)
    else:
        if copy:
            values = values.copy()

    if values.ndim == 1:
        N = values.shape[0]
        if N == 0:
            values = values.reshape((values.shape[0], 0))
        else:
            values = values.reshape((values.shape[0], 1))
    elif values.ndim != 2:
        raise Exception('Must pass 2-d input')

    return values


def _rec_to_dict(arr):
    columns = list(arr.dtype.names)
    sdict = dict((k, arr[k]) for k in columns)
    return columns, sdict

def _homogenize(data, index, columns, dtype=None):
    homogenized = {}

    if dtype is not None:
        dtype = np.dtype(dtype)

    for k in columns:
        if k not in data:
            # no obvious "empty" int column
            if dtype is not None and issubclass(dtype.type, np.integer):
                continue

            v = np.empty(len(index), dtype=dtype)
            v.fill(nan)
        else:
            v = data[k]

        if isinstance(v, Series):
            if dtype is not None:
                v = v.astype(dtype)
            if v.index is not index:
                # Forces alignment. No need to copy data since we
                # are putting it into an ndarray later
                v = v.reindex(index)
        else:
            if isinstance(v, dict):
                v = [v.get(i, nan) for i in index]
            elif np.isscalar(v):
                _v = np.empty(len(index), dtype=_infer_dtype(v))
                _v.fill(v)
                v = _v
            else:
                assert(len(v) == len(index))

            # only *attempt* to cast to dtype
            try:
                arr = np.asarray(v, dtype=dtype)
                if issubclass(arr.dtype.type, basestring):
                    arr = np.array(v, dtype=object, copy=False)
                v = arr
            except Exception:
                v = np.asarray(v)

        homogenized[k] = v

    return homogenized

def _put_str(s, space):
    return ('%s' % s)[:space].ljust(space)

if __name__ == '__main__':
    import nose
    nose.runmodule(argv=[__file__,'-vvs','-x','--pdb', '--pdb-failure'],
                   exit=False)
