# pylint: disable=E1101,E1103
# pylint: disable=W0212,W0231,W0703,W0622

from cStringIO import StringIO
from datetime import datetime
import operator
import sys
import warnings

from numpy import nan
import numpy as np

from pandas.core.common import (isnull, notnull, _check_step, _is_list_like,
                                _need_slice, _is_label_slice, _ensure_index,
                                _try_sort, _pfixed)
from pandas.core.daterange import DateRange
from pandas.core.generic import PandasGeneric
from pandas.core.index import Index, NULL_INDEX
from pandas.core.internals import BlockManager, make_block
from pandas.core.series import Series
import pandas.core.common as common
import pandas.core.datetools as datetools
import pandas.lib.tseries as tseries

#-------------------------------------------------------------------------------
# Factory helper methods

_arith_doc ="""
Arithmetic method: %s

Parameters
----------
other : Series, DataFrame, or constant
axis : {0, 1, 'index', 'columns'}
    For Series input, axis to match Series index on

Notes
-----
Mismatched indices will be unioned together

Returns
-------
result : DataFrame
"""

def _arith_method(func, name, default_axis='columns'):
    def f(self, other, axis=default_axis):
        if isinstance(other, DataFrame):    # Another DataFrame
            return self._combine_frame(other, func)
        elif isinstance(other, Series):
            if axis is not None:
                axis = self._get_axis_name(axis)
                if axis == 'index':
                    return self._combine_match_index(other, func)
                else:
                    return self._combine_match_columns(other, func)
            return self._combine_series_infer(other, func)
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

class DataFrame(PandasGeneric):
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
    copy : boolean, default True
        Copy data from inputs

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
                 copy=True):

        if data is None:
            data = {}

        if isinstance(data, BlockManager):
            mgr = data
            if copy and dtype is None:
                mgr = mgr.copy()
            if dtype is not None:
                # no choice but to copy
                mgr = mgr.cast(dtype)
        elif isinstance(data, DataFrame):
            mgr = data._data
            if copy and dtype is None:
                mgr = mgr.copy()
            elif dtype is not None:
                # no choice but to copy
                mgr = mgr.cast(dtype)
        elif isinstance(data, dict):
            mgr = self._init_dict(data, index, columns, dtype=dtype)
        elif isinstance(data, (np.ndarray, list)):
            mgr = self._init_matrix(data, index, columns, dtype=dtype,
                                    copy=copy)
        else:
            raise Exception('DataFrame constructor not properly called!')

        self._data = mgr

    def _init_dict(self, data, index, columns, dtype=None, copy=True):
        """
        Segregate Series based on type and coerce into matrices.

        Needs to handle a lot of exceptional cases.

        Somehow this got outrageously complicated
        """
        from pandas.core.internals import form_blocks

        # TODO: deal with emptiness!
        # prefilter if columns passed
        if columns is not None:
            columns = _ensure_index(columns)
            data = dict((k, v) for k, v in data.iteritems() if k in columns)

        # figure out the index, if necessary
        if index is None:
            index = extract_index(data)

        # don't force copy because getting jammed in an ndarray anyway
        homogenized = _homogenize_series(data, index, dtype, force_copy=False)
        # segregates dtypes and forms blocks matching to columns
        blocks, columns = form_blocks(homogenized, index, columns)

        # TODO: need consolidate here?
        return BlockManager(blocks, index, columns).consolidate()

    def _init_matrix(self, values, index, columns, dtype=None,
                     copy=True):
        from pandas.core.internals import make_block

        values = _prep_ndarray(values)

        if values.ndim == 1:
            N = values.shape[0]
            if N == 0:
                values = values.reshape((values.shape[0], 0))
            else:
                values = values.reshape((values.shape[0], 1))

        if dtype is not None:
            try:
                values = values.astype(dtype)
            except Exception:
                pass

        N, K = values.shape

        if index is None:
            index = _default_index(N)

        if columns is None:
            columns = _default_index(K)

        columns = _ensure_index(columns)
        block = make_block(values, columns, columns)
        return BlockManager([block], index, columns)

    def astype(self, dtype):
        new_data = self._data.cast(dtype)
        return DataFrame(new_data, copy=False)

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
        return DataFrame(self._data.copy(), copy=False)

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
    def fromRecords(cls, data, indexField=None):
        """
        Convert structured or record ndarray to DataFrame

        Parameters
        ----------
        input : NumPy structured array

        Returns
        -------
        DataFrame
        """
        # Dtype when you have records
        if not issubclass(data.dtype.type, np.void):
            raise Exception('Input was not a structured array!')

        columns = data.dtype.names
        sdict = dict((k, data[k]) for k in columns)

        if indexField is not None:
            index = sdict.pop(indexField)
            columns = [c for c in columns if c != indexField]
        else:
            index = np.arange(len(data))

        return cls(sdict, index=index, columns=columns)

    def toRecords(self, index=True):
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
        data = np.genfromtxt(path, delimiter=delimiter, dtype=None,
                             skip_header=header, names=True)

        if index_col is not None:
            field = data.dtype.names[index_col]
            df = cls.fromRecords(data, indexField=field)

            # have dates?
            test_val = datetools.to_datetime(df.index[0])
            if isinstance(test_val, datetime):
                df = df.rename(index=datetools.to_datetime)
        else:
            df = cls.fromRecords(data, indexField=None)

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
        """Output a tab-separated version of this DataFrame"""
        series = self._series
        if columns is None:
            columns = self.columns
        else:
            columns = [c for c in columns if c in self]

        formatters = formatters or {}
        ident = lambda x: x

        if colSpace is None:
            colSpace = {}

            for c in columns:
                if np.issctype(series[c].dtype):
                    colSpace[c] = max(len(str(c)) + 4, 12)
                else:
                    # HACK
                    colSpace[c] = 15
        else:
            colSpace = dict((k, colSpace) for k in columns)

        if len(columns) == 0 or len(self.index) == 0:
            print >> buf, 'Empty %s' % type(self).__name__
            print >> buf, repr(self.index)
        else:
            idxSpace = max([len(str(idx)) for idx in self.index]) + 4
            head = ' ' * idxSpace

            for h in columns:
                head += _put_str(h, colSpace[h])

            print >> buf, head

            for idx in self.index:
                ot = _put_str(idx, idxSpace - 1)
                for k in columns:
                    formatter = formatters.get(k, ident)
                    ot += _pfixed(formatter(series[k][idx]),
                                  colSpace[k], nanRep=nanRep,
                                  float_format=float_format)
                print >> buf, ot

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
                                  (_put_str(col, space), count))

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

    def _set_columns(self, cols):
        if len(cols) != len(self.columns):
            raise Exception('Length mismatch (%d vs %d)'
                            % (len(cols), len(self.columns)))

        self._data.columns = _ensure_index(cols)

    def _set_index(self, index):
        if len(index) > 0:
            if len(index) != len(self.index):
                raise Exception('Length mismatch (%d vs %d)'
                                % (len(index), len(self.index)))

        self._data.index = _ensure_index(index)

    def _get_index(self):
        return self._data.index

    def _get_columns(self):
        return self._data.columns

    def _get_values(self):
        self._consolidate_inplace()
        return self._data.as_matrix()

    index = property(fget=lambda self: self._get_index(),
                     fset=lambda self, x: self._set_index(x))
    columns = property(fget=lambda self: self._get_columns(),
                     fset=lambda self, x: self._set_columns(x))
    values = property(fget=_get_values)

    def as_matrix(self, columns=None):
        """
        Convert the frame to its Numpy-array matrix representation

        Columns are presented in sorted order unless a specific list
        of columns is provided.
        """
        return self._data.as_matrix(columns)

    asMatrix = as_matrix
    # For DataFrame compatibility

    def transpose(self):
        """
        Returns a DataFrame with the rows/columns switched.
        """
        return DataFrame(data=self.values.T, index=self.columns,
                         columns=self.index)
    T = property(transpose)

    #----------------------------------------------------------------------
    # Picklability

    def __getstate__(self):
        return self._data

    def __setstate__(self, state):
        # old DataFrame pickle
        if len(state) == 3:
            self._unpickle_frame_compat(state)
        # old DataFrame pickle
        elif len(state) == 2: # pragma: no cover
            # old pickling format, for compatibility
            self._unpickle_matrix_compat(state)
        else:
            assert(isinstance(state, BlockManager))
            self._data = state

    def _unpickle_frame_compat(self, state): # pragma: no cover
        from pandas.core.common import _unpickle_array
        series, cols, idx = state
        columns = _unpickle_array(cols)
        index = _unpickle_array(idx)
        self._data = self._init_dict(series, index, columns, None)

    def _unpickle_matrix_compat(self, state): # pragma: no cover
        from pandas.core.common import _unpickle_array
        # old unpickling
        (vals, idx, cols), object_state = state

        index = _unpickle_array(idx)
        dm = DataFrame(vals, index=index,
                        columns=_unpickle_array(cols))

        if object_state is not None:
            ovals, _, ocols = object_state
            objects = DataFrame(ovals, index=index,
                                 columns=_unpickle_array(ocols))

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

    def _union_columns(self, other):
        union_cols = self.columns

        if not union_cols.equals(other.columns):
            union_cols = union_cols.union(other.columns)

        return union_cols

    def _union_index(self, other):
        union_index = self.index

        if not union_index.equals(other.index):
            union_index = union_index.union(other.index)

        return union_index

    #----------------------------------------------------------------------
    # Consolidation of internals

    def _consolidate_inplace(self):
        self._data = self._data.consolidate()

    def consolidate(self):
        """
        Compute DataFrame with "consolidated" internals (data of each dtype
        grouped together in a single ndarray). Mainly an internal API function,
        but available here to the savvy user

        Returns
        -------
        consolidated : DataFrame
        """
        cons_data = self._data.consolidate()
        if cons_data is self._data:
            cons_data = cons_data.copy()
        return DataFrame(cons_data, copy=False)

    #----------------------------------------------------------------------
    # Array interface

    def __array__(self):
        return self.values

    def __array_wrap__(self, result):
        return DataFrame(result, index=self.index, columns=self.columns)

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
            new_data = self._data.get_slice(item)
            return DataFrame(new_data)
        elif isinstance(item, np.ndarray):
            if len(item) != len(self.index):
                raise Exception('Item wrong length %d instead of %d!' %
                                (len(item), len(self.index)))
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
                raise Exception('Can only index with like-indexed '
                                'DataFrame objects')

            self._boolean_set(key, value)
        else:
            self._insert_item(key, value)

    def _boolean_set(self, key, value):
        mask = key.values
        if mask.dtype != np.bool_:
            raise Exception('Must pass DataFrame with boolean values only')

        self.values[mask] = value

    def _insert_item(self, key, value):
        """
        Add series to DataFrame in specified column.

        If series is a numpy-array (not a Series/TimeSeries), it must be the
        same length as the DataFrame's index or an error will be thrown.

        Series/TimeSeries will be conformed to the DataFrame's index to
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
        return self._data.get_series_dict(self.index)

    def xs(self, key, copy=True):
        """
        Returns a row from the DataFrame as a Series object.

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
        return self._data.xs(loc, copy=copy)

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

        # TODO: want to preserve dtypes though...
        new_data = self._data.reindex_index(new_index, method)
        return DataFrame(new_data)

    def _reindex_columns(self, new_columns):
        new_data = self._data.reindex_columns(new_columns)
        return DataFrame(new_data)

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
            raise Exception('items was None!')

    def select(self, crit, axis=0):
        """
        Return data corresponding to axis labels matching criteria

        Parameters
        ----------
        crit : function
            To be called on each index (label). Should return True or False
        axis : {0, 1}

        Returns
        -------
        selection : DataFrame
        """
        return self._select_generic(crit, axis=axis)

    def dropEmptyRows(self, specificColumns=None):
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
        if specificColumns:
            theCount = self.filter(items=specificColumns).count(axis=1)
        else:
            theCount = self.count(axis=1)

        return self.reindex(self.index[theCount != 0])

    def dropIncompleteRows(self, specificColumns=None, minObs=None):
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
        N = len(self.columns)

        if specificColumns:
            colSet = set(specificColumns)
            intersection = set(self.columns) & colSet

            N = len(intersection)

            filtered = self.filter(items=intersection)
            theCount = filtered.count(axis=1)
        else:
            theCount = self.count(axis=1)

        if minObs is None:
            minObs = N

        return self.reindex(self.index[theCount >= minObs])

    #----------------------------------------------------------------------
    # Sorting

    def sort(self, column=None, ascending=True):
        if column:
            series = self[column].order(missingAtEnd=False)
            sort_index = series.index
        else:
            index = np.asarray(self.index)
            argsorted = np.argsort(index)
            sort_index = index[argsorted.astype(int)]

        if not ascending:
            sort_index = sort_index[::-1]

        return self.reindex(sort_index)

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
            return DataFrame(result, index=self.index)
        else:
            # Float type values
            if len(self.columns) == 0:
                return self

            new_data = self._data.fillna(value)
            return DataFrame(new_data)

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

        result = self.copy()

        if index is not None:
            result._rename_index_inplace(index)

        if columns is not None:
            result._rename_columns_inplace(columns)

        return result

    def _rename_index_inplace(self, mapper):
        self._data = self._data.rename_index(mapper)

    def _rename_columns_inplace(self, mapper):
        self._data = self._data.rename_columns(mapper)

    #----------------------------------------------------------------------
    # Arithmetic / combination related

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
            return DataFrame(index=new_index)
        elif not self:
            return other * nan
        elif not other:
            return self * nan

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

        return DataFrame(func(this.values, other.values),
                          index=new_index, columns=new_columns)

    def _indexed_same(self, other):
        same_index = self.index.equals(other.index)
        same_columns = self.columns.equals(other.columns)
        return same_index and same_columns

    def _combine_series_infer(self, other, func):
        if len(other) == 0:
            return self * nan

        if len(self) == 0:
            # Ambiguous case, use _series so works with DataFrame
            return self._constructor(data=self._series, index=self.index,
                                     columns=self.columns)

        # teeny hack because one does DataFrame + TimeSeries all the time
        if self.index.is_all_dates() and other.index.is_all_dates():
            return self._combine_match_index(other, func)
        else:
            return self._combine_match_columns(other, func)

    def _combine_match_index(self, other, func):
        new_index = self._union_index(other)
        values = self.values
        other_vals = other.values

        # Operate row-wise
        if not other.index.equals(new_index):
            other_vals = other.reindex(new_index).values

        if not self.index.equals(new_index):
            values = self.reindex(new_index).values

        return DataFrame(func(values.T, other_vals).T,
                          index=new_index, columns=self.columns)

    def _combine_match_columns(self, other, func):
        newCols = self.columns.union(other.index)

        # Operate column-wise
        this = self.reindex(columns=newCols)
        other = other.reindex(newCols).values

        return DataFrame(func(this.values, other),
                          index=self.index, columns=newCols)

    def _combine_const(self, other, func):
        if not self:
            return self

        # TODO: deal with objects
        return DataFrame(func(self.values, other), index=self.index,
                          columns=self.columns)

    def _compare_frame(self, other, func):
        if not self._indexed_same(other):
            raise Exception('Can only compare identically-labeled '
                            'DataFrame objects')

        new_data = {}
        for col in self.columns:
            new_data[col] = func(self[col], other[col])

        return self._constructor(data=new_data, index=self.index,
                                 columns=self.columns)

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
        return self.combine(other, np.add, fill_value=0.)

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
        return self.combine(other, np.multiply, fill_value=1.)

    #----------------------------------------------------------------------
    # Misc methods

    def first_valid_index(self):
        return self.index[self.count(1) > 0][0]

    def last_valid_index(self):
        return self.index[self.count(1) > 0][-1]

    # to avoid API breakage
    _firstTimeWithValue = first_valid_index
    _lastTimeWithValue = last_valid_index

    def head(self):
        return self[:5]

    def tail(self):
        return self[-5:]

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
        from pandas.core.internals import make_block

        if periods == 0:
            return self

        if timeRule is not None and offset is None:
            offset = datetools.getOffset(timeRule)

        def _shift_block(blk, indexer):
            new_values = blk.values.take(indexer, axis=0)
            # convert integer to float if necessary. need to do a lot more than
            # that, handle boolean etc also
            new_values = common.ensure_float(new_values)
            if periods > 0:
                new_values[:periods] = nan
            else:
                new_values[periods:] = nan
            return make_block(new_values, blk.columns, blk.ref_columns)

        if offset is None:
            indexer = self._shift_indexer(periods)
            new_blocks = [_shift_block(b, indexer) for b in self._data.blocks]
            new_data = BlockManager(new_blocks, self.index, self.columns)
        else:
            new_data = self._data.copy()
            new_data.index = self.index.shift(periods, offset)

        return DataFrame(new_data)

    def _shift_indexer(self, periods):
        # small reusable utility
        N = len(self)
        indexer = np.zeros(N, dtype=int)

        if periods > 0:
            indexer[periods:] = np.arange(N - periods)
        else:
            indexer[:periods] = np.arange(-periods, N)

        return indexer

    def truncate(self, before=None, after=None):
        """Function truncate a sorted DataFrame before and/or after
        some particular dates.

        Parameters
        ----------
        before : date
            Truncate before date
        after : date
            Truncate after date

        Returns
        -------
        DataFrame
        """
        before = datetools.to_datetime(before)
        after = datetools.to_datetime(after)
        return self.ix[before:after]

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
            return DataFrame(data=results, index=self.index,
                             columns=self.columns)
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

        return DataFrame(results, index=self.index, columns=self.columns)

    #----------------------------------------------------------------------
    # Merging / joining methods

    def append(self, other):
        """
        Append columns of other to end of this frame's columns and index.

        Columns not in this frame are added as new columns.
        """
        # TODO: with blocks

        if not other:
            return self.copy()

        if not self:
            return other.copy()

        new_index = np.concatenate((self.index, other.index))
        new_columns = self.columns
        new_data = {}

        if not new_columns.equals(other.columns):
            new_columns = self.columns + other.columns

        for column, series in self.iteritems():
            if column in other:
                new_data[column] = series.append(other[column])
            else:
                new_data[column] = series

        for column, series in other.iteritems():
            if column not in self:
                new_data[column] = series

        # TODO: column ordering issues?
        return DataFrame(data=new_data, index=new_index,
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

        new_data = self._data.join_on(other._data, self[on])
        return DataFrame(new_data)

    # def _join_on(self, other, on):
    #     # Check for column overlap
    #     overlap = set(self.columns) & set(other.columns)

    #     if overlap:
    #         raise Exception('Columns overlap: %s' % _try_sort(overlap))

    #     if len(other.index) == 0:
    #         result = self.copy()

    #         for col in other:
    #             result[col] = nan

    #         return result

    #     indexer, mask = other.index.get_indexer(self[on])
    #     notmask = -mask
    #     need_mask = notmask.any()

    #     new_data = {}

    #     for col, series in other.iteritems():
    #         arr = series.view(np.ndarray).take(indexer)

    #         if need_mask:
    #             arr = common.ensure_float(arr)
    #             arr[notmask] = nan

    #         new_data[col] = arr

    #     new_data.update(self._series)

    #     return self._constructor(new_data, index=self.index)

    def _join_index(self, other, how):
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

        result_series = self.reindex(join_index)._series
        other_series = other.reindex(join_index)._series

        for col in other_series:
            if col in result_series:
                raise Exception('Overlapping columns!')

        result_series.update(other_series)

        return self._constructor(result_series, index=join_index)

    #----------------------------------------------------------------------
    # Data reshaping

    def pivot(self, index=None, columns=None, values=None):
        """
        Produce 'pivot' table based on 3 columns of this DataFrame.
        Uses unique values from index / columns and fills with values.

        Parameters
        ----------
        index : string or object
            Column name to use to make new frame's index
        columns : string or object
            Column name to use to make new frame's columns
        values : string or object
            Column name to use for populating new frame's values
        """
        from pandas.core.panel import pivot
        return pivot(self[index], self[columns], self[values])

    #----------------------------------------------------------------------
    # groupby

    def groupby(self, mapper, axis=0):
        """
        Goup series using mapper (dict or key function, apply given
        function to group, return result as series).

        Parameters
        ----------
        mapper : function, dict or Series
            Called on each element of the object index to determine
            the groups.  If a dict or Series is passed, the Series or
            dict VALUES will be used to determine the groups

        Returns
        -------
        GroupBy object
        """
        from pandas.core.groupby import groupby
        return groupby(self, mapper, axis=axis)

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
        mat = self.asMatrix(cols).T
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
        left += right * 0
        right += left * 0

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

    def _get_agg_axis(self, axis_num):
        if axis_num == 0:
            return self.columns
        elif axis_num == 1:
            return self.index
        else:
            raise Exception('Must have 0<= axis <= 1')

    def count(self, axis=0):
        """
        Return array or Series of # observations over requested axis.

        Parameters
        ----------
        axis : {0, 1}
            0 for row-wise, 1 for column-wise

        Notes
        -----
        Also examines non-float data and checks for None and NaN in such data

        Returns
        -------
        Series or TimeSeries
        """
        try:
            cols = self.columns
            values = self.asMatrix(cols)

            if axis == 0:
                axis_labels = cols
            else:
                axis_labels = self.index

            mask = np.empty(values.shape, dtype=bool)
            mask.flat = notnull(values.ravel())
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
        num_cols = self._get_numeric_columns()

        if len(num_cols) < len(self.columns) and numeric_only:
            y = self.as_matrix(num_cols)
            axis_labels = num_cols
        else:
            y = self.values.copy()
            axis_labels = self._get_agg_axis(axis)

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
            the_sum[the_count == 0] = nan

        return Series(the_sum, index=axis_labels)

    def _get_numeric_columns(self):
        return [col for col in self.columns
                if issubclass(self[col].dtype.type, np.number)]

    def _get_object_columns(self):
        return [col for col in self.columns if self[col].dtype == np.object_]


    def cumsum(self, axis=0):
        """
        Return DataFrame of cumulative sums over requested axis.

        Parameters
        ----------
        axis : {0, 1}
            0 for row-wise, 1 for column-wise

        Returns
        -------
        y : DataFrame
        """
        y = np.array(self.values, subok=True)
        if not issubclass(y.dtype.type, np.int_):
            mask = np.isnan(self.values)
            y[mask] = 0
            result = y.cumsum(axis)
            has_obs = (-mask).astype(int).cumsum(axis) > 0
            result[-has_obs] = np.nan
        else:
            result = y.cumsum(axis)
        return DataFrame(result, index=self.index,
                          columns=self.columns)

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

    def cumprod(self, axis=0):
        """
        Return cumulative product over requested axis as DataFrame

        Parameters
        ----------
        axis : {0, 1}
            0 for row-wise, 1 for column-wise

        Returns
        -------
        y : DataFrame
        """
        def get_cumprod(y):
            y = np.array(y)
            mask = isnull(y)
            if not issubclass(y.dtype.type, np.int_):
                y[mask] = 1
            result = y.cumprod()

            return result

        return self.apply(get_cumprod, axis=axis)

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

        Returns
        -------
        Series or TimeSeries
        """
        summed = self.sum(axis, numeric_only=True)
        count = self.count(axis).astype(float)

        if not count.index.equals(summed.index):
            count = count.reindex(summed.index)

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
            if arr.dtype != np.float_:
                arr = arr.astype(float)
            return scoreatpercentile(arr[notnull(arr)], per)

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
            return tseries.median(arr[notnull(arr)])

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
        y = np.asarray(self.values)
        mask = np.isnan(y)
        count = (y.shape[axis] - mask.sum(axis)).astype(float)
        y[mask] = 0

        X = y.sum(axis)
        XX = (y**2).sum(axis)

        theVar = (XX - X**2 / count) / (count - 1)

        return Series(theVar, index=self._get_agg_axis(axis))

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
        y = np.asarray(self.values)
        mask = np.isnan(y)
        count = (y.shape[axis] - mask.sum(axis)).astype(float)
        y[mask] = 0

        A = y.sum(axis) / count
        B = (y**2).sum(axis) / count  - A**2
        C = (y**3).sum(axis) / count - A**3 - 3*A*B

        theSkew = (np.sqrt((count**2-count))*C) / ((count-2)*np.sqrt(B)**3)

        return Series(theSkew, index=self._get_agg_axis(axis))

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

    def plot(self, kind='line', subplots=False, sharex=True,
             sharey=False, use_index=True, **kwds): # pragma: no cover
        """
        Plot the DataFrame's series with the index on the x-axis using
        matplotlib / pylab.

        Parameters
        ----------
        kind : {'line', 'bar', 'hist'}
            Default: line for TimeSeries, hist for Series

        kwds : other plotting keyword arguments

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

    def toDataMatrix(self):
        warnings.warn("toDataMatrix will disappear in next release "
                      "as there is no longer a DataMatrix class", FutureWarning)
        return self.copy()

    def rows(self):
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

    def merge(self, *args, **kwargs):
        warnings.warn("merge is deprecated. Use 'join' instead",
                      FutureWarning)
        return self.join(*args, **kwargs)

    #----------------------------------------------------------------------
    # Fancy indexing

    _ix = None
    @property
    def ix(self):
        if self._ix is None:
            self._ix = _DataFrameIndexer(self)

        return self._ix

    def _fancy_getitem(self, key, axis=0):
        labels = self._get_axis(axis)
        axis_name = self._get_axis_name(axis)

        # asarray can be unsafe, NumPy strings are weird
        isbool = np.asarray(key).dtype == np.bool_
        if isbool:
            if isinstance(key, Series):
                if not key.index.equals(labels):
                    raise Exception('Cannot use boolean index with misaligned '
                                    'or unequal labels')
            return self.reindex(**{axis_name : labels[key]})
        else:
            return self.reindex(**{axis_name : key})

    def _fancy_getitem_tuple(self, rowkey, colkey):
        result = self._fancy_getitem_axis(colkey, axis=1)

        if isinstance(result, Series):
            result = result[rowkey]
        else:
            result = result._fancy_getitem_axis(rowkey, axis=0)

        return result

    def _fancy_getitem_axis(self, key, axis=0):
        if isinstance(key, slice):
            return self._slice_axis(key, axis=axis)
        elif _is_list_like(key):
            return self._fancy_getitem(key, axis=axis)
        elif axis == 0:
            idx = key
            if isinstance(key, int):
                idx = self.index[key]

            return self.xs(idx)
        else:
            col = key
            if isinstance(key, int):
                col = self.columns[key]

            return self[col]

    def _slice_axis(self, slice_obj, axis=0):
        _check_step(slice_obj)

        if not _need_slice(slice_obj):
            return self

        axis_name = self._get_axis_name(axis)

        labels = getattr(self, axis_name)
        if _is_label_slice(labels, slice_obj):
            i, j = labels.slice_locs(slice_obj.start, slice_obj.stop)
            new_labels = labels[i:j]
        else:
            new_labels = labels[slice_obj]

        return self.reindex(**{axis_name : new_labels})

class _DataFrameIndexer(object):
    """
    Class to support fancy indexing, potentially using labels of DataFrame

    Notes
    -----
    Indexing based on labels is INCLUSIVE
    Slicing uses PYTHON SEMANTICS (endpoint is excluded)

    If Index contains int labels, these will be used rather than the locations,
    so be very careful (ambiguous).

    Examples
    --------
    >>> frame.ix[5:10, ['A', 'B']]
    >>> frame.ix[date1:date2, 'A']
    """

    def __init__(self, frame):
        self.frame = frame

    def __getitem__(self, key):
        frame = self.frame
        if isinstance(key, slice):
            return frame._fancy_getitem_axis(key, axis=0)
        elif isinstance(key, tuple):
            if len(key) != 2:
                raise Exception('only length 2 tuple supported')
            return frame._fancy_getitem_tuple(*key)
        elif _is_list_like(key):
            return frame._fancy_getitem(key, axis=0)
        else:
            return frame._fancy_getitem_axis(key, axis=0)

    def __setitem__(self, key, value):
        raise NotImplementedError

def extract_index(data):
    def _union_if(index, new_index):
        if index is None:
            index = new_index
        else:
            index = index.union(new_index)
        return index

    index = None
    if len(data) == 0:
        index = NULL_INDEX
    elif len(data) > 0 and index is None:
        _check_data_types(data)

        # this is still kludgier than I'd like
        for v in data.values():
            if isinstance(v, Series):
                index = _union_if(index, v.index)
            elif isinstance(v, dict):
                index = _union_if(index, Index(_try_sort(v)))
            else: # not dict-like, assign integer labels
                index = Index(np.arange(len(v)))

    if len(index) == 0:
        index = NULL_INDEX

    return _ensure_index(index)

def _check_data_types(data):
    have_series = False
    have_raw_arrays = False
    for v in data.values():
        if isinstance(v, (dict, Series)):
            have_series = True
        else:
            have_raw_arrays = True

    if have_series and have_raw_arrays:
        raise Exception('Cannot mix Series / dict objects'
                        ' with ndarray / sequence input')

def _prep_ndarray(values):
    if not isinstance(values, np.ndarray):
        values = np.asarray(values)
        # NumPy strings are a pain, convert to object
        if issubclass(values.dtype.type, basestring):
            values = np.array(values, dtype=object, copy=True)
    return values

def _homogenize_series(data, index, dtype=None, force_copy=True):
    homogenized = {}

    for k, v in data.iteritems():
        if isinstance(v, Series):
            if dtype is not None:
                v = v.astype(dtype)
            if v.index is not index:
                # Forces alignment. No need to copy data since we
                # are putting it into an ndarray later
                v = v.reindex(index)
            elif force_copy:
                # same index, but want to copy
                v = v.copy()
        else:
            if isinstance(v, dict):
                v = [v.get(i, nan) for i in index]
            else:
                assert(len(v) == len(index))
            try:
                v = Series(v, dtype=dtype, index=index)
            except Exception:
                v = Series(v, index=index)

            if force_copy:
                v = v.copy()

        # # OK, I will relent for now.
        # if not issubclass(v.dtype.type, (float, int)):
        # #     v = v.astype(np.float64)
        # # else:
        #     v = v.astype(object)

        homogenized[k] = v

    return homogenized

def _default_index(n):
    if n == 0:
        return NULL_INDEX
    else:
        return np.arange(n)

def _put_str(s, space):
    return ('%s' % s)[:space].ljust(space)

def _reorder_columns(mat, current, desired):
    indexer, mask = common.get_indexer(current, desired, None)
    return mat.take(indexer[mask], axis=1)

if __name__ == '__main__':
    import nose
    nose.runmodule(argv=[__file__,'-vvs','-x','--pdb', '--pdb-failure'],
                   exit=False)
