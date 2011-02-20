# pylint: disable=E1101,E1103
# pylint: disable=W0212,W0231,W0703,W0622

from cStringIO import StringIO
from datetime import datetime
import operator
import sys
import warnings

from numpy import NaN
import numpy as np

from pandas.core.common import (_pickle_array, _unpickle_array, _pfixed,
                                isnull, notnull)
from pandas.core.daterange import DateRange
from pandas.core.index import Index, NULL_INDEX
from pandas.core.mixins import Picklable, Groupable
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

def arith_method(func, name, default_axis='columns'):
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

_AXIS_NUMBERS = {
    'index' : 0,
    'columns' : 1
}

_AXIS_NAMES = dict((v, k) for k, v in _AXIS_NUMBERS.iteritems())

#-------------------------------------------------------------------------------
# DataFrame class

class DataFrame(Picklable, Groupable):
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

    Notes
    -----
    Data contained within is COPIED from input arrays, this is to prevent silly
    behavior like altering the original arrays and having those changes
    reflected in the frame.

    See also
    --------
    DataMatrix: more efficient version of DataFrame for most operations

    Examples
    --------
        >>> d = {'col1' : ts1, 'col2' : ts2}
        >>> df = DataFrame(data=d, index=someIndex)
    """
    _columns = None

    def __init__(self, data=None, index=None, columns=None, dtype=None):
        if isinstance(data, dict):
            sdict, columns, index = self._init_dict(data, index, columns, dtype)
        elif isinstance(data, (np.ndarray, list)):
            sdict, columns, index = self._init_matrix(data, index, columns,
                                                      dtype)
        elif isinstance(data, DataFrame):
            sdict = data._series.copy()

            if dtype is not None:
                sdict = dict((k, v.astype(dtype)) for k, v in data.iteritems())
            index = data.index
            columns = data.columns
        elif data is None:
            sdict = {}

            if index is None:
                index = NULL_INDEX

            if columns is None:
                columns = NULL_INDEX
            else:
                for c in columns:
                    sdict[c] = Series(np.NaN, index=index)

        self._series = sdict
        self.columns = columns
        self.index = index

    def _init_dict(self, data, index, columns, dtype):
        # pre-filter out columns if we passed it
        if columns is not None:
            if not isinstance(columns, Index):
                columns = Index(columns)

            data = dict((k, v) for k, v in data.iteritems() if k in columns)
        else:
            columns = Index(_try_sort(data.keys()))

        index = _extract_index(data, index)

        sdict = {}
        for k, v in data.iteritems():
            if isinstance(v, Series):
                # Forces alignment and copies data
                sdict[k] = v.reindex(index)
            else:
                if isinstance(v, dict):
                    v = [v.get(i, NaN) for i in index]

                try:
                    v = Series(v, dtype=dtype, index=index)
                except Exception:
                    v = Series(v, index=index)

                sdict[k] = v.copy()

        # add in any other columns we want to have (completeness)
        for c in columns:
            if c not in sdict:
                sdict[c] = Series(np.NaN, index=index)

        return sdict, columns, index

    def _init_matrix(self, data, index, columns, dtype):
        if not isinstance(data, np.ndarray):
            arr = np.array(data)
            if issubclass(arr.dtype.type, basestring):
                arr = np.array(data, dtype=object, copy=True)

            data = arr

        if data.ndim == 1:
            data = data.reshape((len(data), 1))
        elif data.ndim != 2:
            raise Exception('Must pass 2-d input!')

        N, K = data.shape

        if index is None:
            index = _default_index(N)

        if columns is None:
            columns = _default_index(K)

        if len(columns) != K:
            raise Exception('Column length mismatch: %d vs. %d' %
                            (len(columns), K))

        data = dict([(idx, data[:, i]) for i, idx in enumerate(columns)])
        return self._init_dict(data, index, columns, dtype)

    @classmethod
    def _get_axis_number(cls, axis):
        if axis in (0, 1):
            return axis
        else:
            return _AXIS_NUMBERS[axis]

    @classmethod
    def _get_axis_name(cls, axis):
        if axis in _AXIS_NUMBERS:
            return axis
        else:
            return _AXIS_NAMES[axis]

    def _get_axis(self, axis):
        results = {
            0 : self.index,
            1 : self.columns,
        }

        return results[self._get_axis_number(axis)]

    @property
    def _constructor(self):
        return DataFrame

    def __getstate__(self):
        series = dict((k, v.values) for k, v in self.iteritems())
        columns = _pickle_array(self.columns)
        index = _pickle_array(self.index)

        return series, columns, index

    def __setstate__(self, state):
        # for compatibility with old pickle files
        if len(state) == 2: # pragma: no cover
            series, idx = state
            columns = sorted(series)
        else:
            series, cols, idx = state
            columns = _unpickle_array(cols)

        index = _unpickle_array(idx)
        self._series = dict((k, Series(v, index=index))
                            for k, v in series.iteritems())
        self.index = index
        self.columns = columns

    _index = None
    def _set_index(self, index):
        if isinstance(index, Index):
            self._index = index
        else:
            self._index = Index(index)

        for v in self._series.values():
            v.index = self._index

    def _get_index(self):
        return self._index

    index = property(fget=lambda self: self._get_index(),
                     fset=lambda self, x: self._set_index(x))

    def _get_columns(self):
        return self._columns

    def _set_columns(self, cols):
        if len(cols) != len(self._series):
            raise Exception('Columns length %d did not match data %d!' %
                            (len(cols), len(self._series)))

        if not isinstance(cols, Index):
            cols = Index(cols)

        self._columns = cols

    def _insert_column_index(self, key, loc):
        if loc == len(self.columns):
            columns = Index(np.concatenate((self.columns, [key])))
        elif loc == 0:
            columns = Index(np.concatenate(([key], self.columns)))
        else:
            columns = Index(np.concatenate((self.columns[:loc], [key],
                                            self.columns[loc:])))

        self.columns = columns

    def _get_insert_loc(self, key):
        try:
            loc = self.columns.searchsorted(key)
        except TypeError:
            loc = len(self.columns)

        return loc

    def _delete_column_index(self, loc):
        if loc == len(self.columns) - 1:
            new_columns = self.columns[:loc]
        else:
            new_columns = Index(np.concatenate((self.columns[:loc],
                                               self.columns[loc+1:])))
        self.columns = new_columns

    columns = property(fget=lambda self: self._get_columns(),
                     fset=lambda self, x: self._set_columns(x))

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
        y : DataFrame or DataMatrix
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

    def toRecords(self, index=True):
        """
        Convert DataFrame to record array. Index will be put in the
        'index' field of the record array.

        Returns
        -------
        y : recarray
        """
        if index:
            arrays = [self.index] + [self[c] for c in self.cols()]
            names = ['index'] + list(self.cols())
        else:
            arrays = [self[c] for c in self.cols()]
            names = list(self.cols())

        return np.rec.fromarrays(arrays, names=names)

#-------------------------------------------------------------------------------
# Magic methods

    def __array__(self):
        return self.values

    def __array_wrap__(self, result):
        return DataFrame(result, index=self.index, columns=self.columns)

    def __nonzero__(self):
        return len(self.columns) > 0 and len(self.index) > 0

    def __repr__(self):
        """
        Return a string representation for a particular DataFrame
        """
        buf = StringIO()
        if len(self.index) < 500 and len(self.columns) < 10:
            self.toString(buffer=buf)
        else:
            buf.write(str(self.__class__) + '\n')
            self.info(buffer=buf)

        return buf.getvalue()

    def __getitem__(self, item):
        """
        Retrieve column or slice from DataFrame
        """
        try:
            return self._series[item]
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
            else:
                raise

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
        columns = self.columns
        values = self.values

        if mask.dtype != np.bool_:
            raise Exception('Must pass DataFrame with boolean values only')

        values[mask] = value
        values = values.T
        self._series = dict((c, Series(values[i], index=self.index))
                            for i, c in enumerate(columns))

    def _insert_item(self, key, value):
        if hasattr(value, '__iter__'):
            if isinstance(value, Series):
                cleanSeries = value.reindex(self.index)
            else:
                cleanSeries = Series(value, index=self.index)

            self._series[key] = cleanSeries
        # Scalar
        else:
            self._series[key] = Series(value, index=self.index)

        if key not in self.columns:
            loc = self._get_insert_loc(key)
            self._insert_column_index(key, loc)

    def __delitem__(self, key):
        """
        Delete column from DataFrame
        """
        loc = self.columns.indexMap[key]
        del self._series[key]
        self._delete_column_index(loc)

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

    def __iter__(self):
        """
        Iterate over columns of the frame.
        """
        return iter(self.columns)

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

#-------------------------------------------------------------------------------
# Arithmetic methods

    add = arith_method(operator.add, 'add')
    mul = arith_method(operator.mul, 'multiply')
    sub = arith_method(operator.sub, 'subtract')
    div = arith_method(operator.div, 'divide')

    radd = arith_method(operator.add, 'add')
    rmul = arith_method(operator.mul, 'multiply')
    rsub = arith_method(lambda x, y: y - x, 'subtract')
    rdiv = arith_method(lambda x, y: y / x, 'divide')

    __add__ = arith_method(operator.add, '__add__', default_axis=None)
    __sub__ = arith_method(operator.sub, '__sub__', default_axis=None)
    __mul__ = arith_method(operator.mul, '__mul__', default_axis=None)
    __div__ = arith_method(operator.div, '__div__', default_axis=None)
    __pow__ = arith_method(operator.pow, '__pow__', default_axis=None)

    __radd__ = arith_method(operator.add, '__radd__', default_axis=None)
    __rmul__ = arith_method(operator.mul, '__rmul__', default_axis=None)
    __rsub__ = arith_method(lambda x, y: y - x, '__rsub__', default_axis=None)
    __rdiv__ = arith_method(lambda x, y: y / x, '__rdiv__', default_axis=None)
    __rpow__ = arith_method(lambda x, y: y ** x, '__rpow__', default_axis=None)

    def __neg__(self):
        return self * -1

#-------------------------------------------------------------------------------
# Comparison methods

    __eq__ = comp_method(operator.eq, '__eq__')
    __lt__ = comp_method(operator.lt, '__lt__')
    __gt__ = comp_method(operator.gt, '__gt__')
    __le__ = comp_method(operator.le, '__le__')
    __ge__ = comp_method(operator.ge, '__ge__')

#-------------------------------------------------------------------------------
# Private / helper methods

    def first_valid_index(self):
        return self.index[self.count(1) > 0][0]

    def last_valid_index(self):
        return self.index[self.count(1) > 0][-1]

    # to avoid API breakage
    _firstTimeWithValue = first_valid_index
    _lastTimeWithValue = last_valid_index

    def _combine_frame(self, other, func):
        new_index = self._union_index(other)
        new_columns = self._union_columns(other)

        this = self
        if self.index is not new_index:
            this = self.reindex(new_index)
            other = other.reindex(new_index)

        if not self and not other:
            return DataFrame(index=new_index)

        if not other:
            return self * NaN

        if not self:
            return other * NaN

        new_data = {}
        for col in new_columns:
            if col in this and col in other:
                new_data[col] = func(this[col], other[col])

        return DataFrame(data=new_data, index=new_index, columns=new_columns)

    def _compare_frame(self, other, func):
        if not self._indexed_same(other):
            raise Exception('Can only compare identically-labeled '
                            'DataFrame objects')

        new_data = {}
        for col in self.columns:
            new_data[col] = func(self[col], other[col])

        return DataFrame(data=new_data, index=self.index,
                         columns=self.columns)

    def _indexed_same(self, other):
        same_index = self.index.equals(other.index)

        # for DataMatrix compat
        same_columns = Index(self.cols()).equals(Index(other.cols()))

        return same_index and same_columns

    def _combine_series_infer(self, other, func):
        if len(other) == 0:
            return self * NaN

        if len(self) == 0:
            # Ambiguous case, use _series so works with DataMatrix
            return self._constructor(data=self._series, index=self.index,
                                     columns=self.columns)

        # teeny hack because one does DataFrame + TimeSeries all the time
        if self.index._allDates and other.index._allDates:
            return self._combine_match_index(other, func)
        else:
            return self._combine_match_columns(other, func)

    def _combine_match_index(self, other, func):
        new_data = {}

        new_index = self._union_index(other)
        this = self
        if self.index is not new_index:
            this = self.reindex(new_index)
            other = other.reindex(new_index)

        for col, series in this.iteritems():
            new_data[col] = func(series, other)

        return DataFrame(new_data, index=new_index, columns=self.columns)

    def _combine_match_columns(self, other, func):
        new_data = {}

        union = intersection = self.columns

        if not union.equals(other.index):
            union = other.index.union(self.columns)
            intersection = other.index.intersection(self.columns)

        for col in intersection:
            new_data[col] = func(self[col], other[col])

        return DataFrame(new_data, index=self.index, columns=union)

    def _combine_const(self, other, func):
        new_data = {}
        for col, series in self.iteritems():
            new_data[col] = func(series, other)

        return DataFrame(data=new_data, index=self.index,
                         columns=self.columns)

#-------------------------------------------------------------------------------
# Public methods

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
            cols = self.cols()

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
            for col in cols:
                val = series[col].get(idx)
                if isnull(val):
                    val = nanRep
                else:
                    val = str(val)
                f.write(',%s' % val)
            f.write('\n')

        f.close()

    def toDataMatrix(self):
        from pandas.core.matrix import DataMatrix
        return DataMatrix(self._series, index=self.index)

    def toString(self, buffer=sys.stdout, columns=None, colSpace=15,
                 nanRep='NaN', formatters=None, float_format=None):
        """Output a tab-separated version of this DataFrame"""
        series = self._series
        if columns is None:
            columns = self.columns
        else:
            columns = [c for c in columns if c in self]

        formatters = formatters or {}
        ident = lambda x: x

        if len(columns) == 0 or len(self.index) == 0:
            print >> buffer, 'Empty DataFrame'
            print >> buffer, repr(self.index)
        else:
            idxSpace = max([len(str(idx)) for idx in self.index]) + 4
            head = _pfixed('', idxSpace)

            for h in columns:
                head += _pfixed(h, colSpace)

            print >> buffer, head

            for idx in self.index:
                ot = _pfixed(idx, idxSpace)
                for k in columns:
                    formatter = formatters.get(k, ident)
                    ot += _pfixed(formatter(series[k][idx]),
                                  colSpace, nanRep=nanRep,
                                  float_format=float_format)
                print >> buffer, ot

    def head(self, buffer=sys.stdout):
        chunk = self[:5]
        if len(self.cols()) > 6:
            print 'Probably too wide to display, transposing'
            chunk = chunk.T

        chunk.toString(buffer=buffer)

    def tail(self, buffer=sys.stdout):
        chunk = self[-5:]
        if len(self.cols()) > 6:
            print 'Probably too wide to display, transposing'
            chunk = chunk.T

        chunk.toString(buffer=buffer)

    def info(self, buffer=sys.stdout):
        """Concise summary of a DataFrame, used in __repr__ when very large."""
        print >> buffer, 'Index: %s entries' % len(self.index),
        if len(self.index) > 0:
            print >> buffer, ', %s to %s' % (self.index[0], self.index[-1])
        else:
            print >> buffer, ''

        if len(self._series) == 0:
            print >> buffer, 'DataFrame is empty!'
            return

        series = self._series
        columns = self.columns
        space = max([len(str(k)) for k in columns]) + 4
        for k in columns:
            out = _pfixed(k, space)
            out += '%d  non-null values' % series[k].count()
            print >> buffer, out

    def rows(self):
        """Alias for the frame's index"""
        return self.index

    def cols(self):
        """Return sorted list of frame's columns"""
        return list(self.columns)

    def iteritems(self):
        """Iterator over (column, series) pairs"""
        return ((k, self._series[k]) for k in self.columns)

    def append(self, other):
        """
        Append columns of other to end of this frame's columns and index.

        Columns not in this frame are added as new columns.
        """
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

    def asfreq(self, freq, method=None, fillMethod=None):
        """
        Convert all TimeSeries inside to specified frequency using
        DateOffset objects. Optionally provide fill method to pad or
        backfill missing values.

        Parameters
        ----------
        offset : DateOffset object, or string in {'WEEKDAY', 'EOM'}
            DateOffset object or subclass (e.g. monthEnd)

        fillMethod : {'backfill', 'pad', None}
                    Method to use for filling holes in new inde
        """
        if len(self.index) == 0:
            return self.copy()

        if isinstance(freq, datetools.DateOffset):
            dateRange = DateRange(self.index[0], self.index[-1], offset=freq)
        else:
            dateRange = DateRange(self.index[0], self.index[-1], timeRule=freq)

        return self.reindex(dateRange, method=method, fillMethod=fillMethod)

    def as_matrix(self, columns=None):
        """
        Convert the frame to its Numpy-array matrix representation

        Columns are presented in sorted order unless a specific list
        of columns is provided.
        """
        if columns is None:
            columns = self.cols()

        if len(columns) == 0:
            return np.zeros((0, 0))

        return np.array([self[col] for col in columns]).T

    asMatrix = as_matrix
    # For DataMatrix compatibility
    values = property(as_matrix)

    def copy(self):
        """
        Make a deep copy of this frame
        """
        return DataFrame(self._series, index=self.index,
                         columns=self.columns)

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
        N = len(self.cols())

        if specificColumns:
            colSet = set(specificColumns)
            intersection = set(self.cols()) & colSet

            N = len(intersection)

            filtered = self.filter(items=intersection)
            theCount = filtered.count(axis=1)
        else:
            theCount = self.count(axis=1)

        if minObs is None:
            minObs = N

        return self.reindex(self.index[theCount >= minObs])

    def fill(self, value=None, method='pad'): # pragma: no cover
        warnings.warn("fill is being replaced by fillna, and the fill function "
                      "behavior will disappear in the next release: please "
                      "modify your code accordingly",
                      FutureWarning)
        return self.fillna(value=value, method=method)

    def fillna(self, value=None, method='pad'):
        """
        Fill NaN values using the specified method.

        Member Series / TimeSeries are filled separately.

        Parameters
        ----------
        method : {'backfill', 'pad', None}
            Method to use for filling holes in new inde

        value : any kind (should be same type as array)
            Value to use to fill holes (e.g. 0)

        Returns
        -------
        DataFrame with NaN's filled

        See also
        --------
        reindex, asfreq
        """
        mycopy = self.copy()
        for col in mycopy._series.keys():
            series = mycopy._series[col]
            filledSeries = series.fillna(method=method, value=value)
            mycopy._series[col] = filledSeries

        return mycopy

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
        beg_slice, end_slice = self._getIndices(before, after)

        return self[beg_slice:end_slice]

    def _getIndices(self, before, after):
        before = datetools.to_datetime(before)
        after = datetools.to_datetime(after)

        if before is None:
            beg_slice = 0
        else:
            beg_slice = self.index.searchsorted(before, side='left')

        if after is None:
            end_slice = len(self.index)
        else:
            end_slice = self.index.searchsorted(after, side='right')

        return beg_slice, end_slice

    def xs(self, key):
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

        subset = self.cols()
        rowValues = [self._series[k][key] for k in subset]

        if len(set((type(x) for x in rowValues))) > 1:
            return Series(np.array(rowValues, dtype=np.object_), index=subset)
        else:
            return Series(np.array(rowValues), index=subset)

    def getXS(self, key): # pragma: no cover
        warnings.warn("'getXS' is deprecated. Use 'xs' instead",
                      FutureWarning)

        return self.xs(key)

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

    def reindex(self, index=None, columns=None, method=None, fillMethod=None):
        """
        Reindex data inside, optionally filling according to some rule.

        Parameters
        ----------
        index : array-like, optional
            preferably an Index object (to avoid duplicating data)
        columns : array-like, optional
        method : {'backfill', 'pad', None}
            Method to use for filling data holes using the index. See
            Series.reindex for more information

        Returns
        -------
        y : same type as calling instance
        """
        # TODO: remove this on next release
        if fillMethod is not None: # pragma: no cover
            warnings.warn("'fillMethod' is deprecated. Use 'method' instead",
                          FutureWarning)

            method = fillMethod

        frame = self

        if index is not None:
            frame = frame._reindex_index(index, method)

        if columns is not None:
            frame = frame._reindex_columns(columns)

        return frame

    def _reindex_index(self, index, method):
        if self.index.equals(index):
            return self.copy()

        if not isinstance(index, Index):
            index = Index(index)

        if len(self.index) == 0:
            return DataFrame(index=index, columns=self.columns)

        indexer, mask = common.get_indexer(self.index, index, method)

        # Maybe this is a bit much? Wish I had more unit tests...
        typeHierarchy = [
            (float, float),
            (int, float),
            (bool, float),
            (np.bool_, float),
            (basestring, object),
            (object, object)
        ]

        missingValue = {
            float  : NaN,
            object : NaN,
            np.bool_ : False
        }

        notmask = -mask
        need_cast = notmask.any()

        newSeries = {}
        for col, series in self.iteritems():
            series = series.view(np.ndarray)
            for klass, dest in typeHierarchy:
                if issubclass(series.dtype.type, klass):
                    new = series.take(indexer)

                    if need_cast:
                        new = new.astype(dest)
                        new[notmask] = missingValue[dest]

                    newSeries[col] = new
                    break

        return DataFrame(newSeries, index=index, columns=self.columns)

    def _reindex_columns(self, columns):
        if not isinstance(columns, Index):
            columns = Index(columns)

        sdict = dict((k, v) for k, v in self.iteritems() if k in columns)
        return DataFrame(sdict, index=self.index, columns=columns)

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
        self.index = [mapper(x) for x in self.index]

    def _rename_columns_inplace(self, mapper):
        new_series = {}
        new_columns = []

        for col in self.columns:
            new_col = mapper(col)
            if new_col in new_series:
                raise Exception('Non-unique mapping!')
            new_series[new_col] = self[col]
            new_columns.append(new_col)

        self.columns = new_columns
        self._series = new_series

    @property
    def T(self):
        """
        Returns a DataFrame with the rows/columns switched.
        """
        # Need to do some 'type inference' to avoid casting
        # float to string in weird cases
        dtypes = list(set([x.dtype for x in self._series.values()]))
        if len(dtypes) > 1:
            theDtype = np.object_
        else:
            theDtype = dtypes[0]

        valuesT = np.array([self[col] for col in self.columns],
                           dtype=theDtype).T

        return DataFrame(data=dict(zip(self.index, valuesT)),
                         index=self.columns, columns=self.index)

    def diff(self, periods=1):
        return self - self.shift(periods)

    def shift(self, periods, offset=None, timeRule=None):
        """
        Shift the underlying series of the DataFrame and Series objects within
        by given number (positive or negative) of business/weekdays.

        Note, nan values inserted at beginning of series.
        """
        if periods == 0:
            return self

        if timeRule and not offset:
            offset = datetools.getOffset(timeRule)

        if offset is None:
            new_index = self.index
            indexer = self._shift_indexer(periods)

            if periods > 0:
                def do_shift(series):
                    values = np.asarray(series).take(indexer)
                    values = common.ensure_float(values)
                    values[:periods] = NaN
                    return values
            else:
                def do_shift(series):
                    values = np.asarray(series).take(indexer)
                    values = common.ensure_float(values)
                    values[periods:] = NaN
                    return values

            new_data = dict([(col, do_shift(series))
                              for col, series in self.iteritems()])
        else:
            new_index = self.index.shift(periods, offset)
            new_data = dict([(col, np.asarray(series))
                               for col, series in self.iteritems()])

        return DataFrame(data=new_data, index=new_index,
                         columns=self.columns)

    def _shift_indexer(self, periods):
        # small reusable utility
        N = len(self)
        indexer = np.zeros(N, dtype=int)

        if periods > 0:
            indexer[periods:] = np.arange(N - periods)
        else:
            indexer[:periods] = np.arange(-periods, N)

        return indexer

    def apply(self, func, axis=0):
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

        Examples
        --------
            >>> df.apply(numpy.sqrt) --> DataFrame
            >>> df.apply(numpy.sum) --> Series

        Notes
        -----
        Functions altering the index are not supported (yet)
        """
        if not len(self.cols()):
            return self

        if axis == 0:
            target = self
            agg_index = self.cols()
        elif axis == 1:
            target = self.T
            agg_index = self.index

        results = {}
        for k in target.cols():
            results[k] = func(target[k])

        if hasattr(results.values()[0], '__iter__'):
            result = self._constructor(data=results, index=target.index,
                                       columns=target.columns)

            if axis == 1:
                result = result.T

            return result
        else:
            return Series(results, index=agg_index)

    def tapply(self, func):
        """
        Apply func to the transposed DataFrame, results as per apply
        """
        return self.apply(func, axis=1)

    def applymap(self, func):
        """
        Apply a function to a DataFrame that is intended to operate elementwise

        Please try to use apply if possible

        Parameters
        ----------
        func : function
            Python function to apply to each element
        """
        results = {}
        for col, series in self.iteritems():
            results[col] = [func(v) for v in series]
        return DataFrame(data=results, index=self.index, columns=self.columns)

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
            columns = [r for r in items if r in self]
        elif like:
            columns = [c for c in self.cols() if like in c]
        elif regex:
            matcher = re.compile(regex)
            columns = [c for c in self.cols() if matcher.match(c)]
        else:
            raise Exception('items was None!')

        return self.reindex(columns=columns)

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

    def combine(self, other, func, fill_value=None):
        """
        Add two DataFrame / DataMatrix objects and do not propagate NaN values,
        so if for a (column, time) one frame is missing a value, it will
        default to the other frame's value (which might be NaN as well)

        Parameters
        ----------
        other : DataFrame / Matrix

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

        new_columns = _try_sort(set(this.cols() + other.cols()))
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
                    arr[this_mask & other_mask] = np.NaN

                result[col] = arr

            elif col in this:
                result[col] = this[col]
            elif col in other:
                result[col] = other[col]

        return self._constructor(result, index=new_index, columns=new_columns)

    def combineFirst(self, other):
        """
        Combine two DataFrame / DataMatrix objects and default to value
        in frame calling the method.

        Parameters
        ----------
        otherFrame : DataFrame / Matrix

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
        Add two DataFrame / DataMatrix objects and do not propagate
        NaN values, so if for a (column, time) one frame is missing a
        value, it will default to the other frame's value (which might
        be NaN as well)

        Parameters
        ----------
        other : DataFrame / Matrix

        Returns
        -------
        DataFrame
        """
        return self.combine(other, np.add, fill_value=0.)

    def combineMult(self, other):
        """
        Multiply two DataFrame / DataMatrix objects and do not
        propagate NaN values, so if for a (column, time) one frame is
        missing a value, it will default to the other frame's value
        (which might be NaN as well)

        Parameters
        ----------
        other : DataFrame / Matrix

        Returns
        -------
        DataFrame
        """
        return self.combine(other, np.multiply, fill_value=1.)

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

    merge = join

    def _join_on(self, other, on):
        # Check for column overlap
        overlap = set(self.cols()) & set(other.cols())

        if overlap:
            raise Exception('Columns overlap: %s' % _try_sort(overlap))

        if len(other.index) == 0:
            result = self.copy()

            for col in other:
                result[col] = np.NaN

            return result

        indexer, mask = tseries.getMergeVec(self[on], other.index.indexMap)
        notmask = -mask
        need_mask = notmask.any()

        new_data = {}

        for col, series in other.iteritems():
            arr = series.view(np.ndarray).take(indexer)

            if need_mask:
                arr = common.ensure_float(arr)
                arr[notmask] = NaN

            new_data[col] = arr

        new_data.update(self._series)

        return self._constructor(new_data, index=self.index)

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

    def plot(self, kind='line', subplots=False, sharex=True, sharey=False,
             **kwds): # pragma: no cover
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

        for i, col in enumerate(_try_sort(self.columns)):
            if subplots:
                ax = axes[i]
                ax.plot(self.index, self[col].values, 'k', label=col,
                        **kwds)
                ax.legend(loc='best')
            else:
                ax.plot(self.index, self[col].values, label=col,
                        **kwds)

    def _get_agg_axis(self, axis_num):
        if axis_num == 0:
            return self.columns
        elif axis_num == 1:
            return self.index
        else:
            raise Exception('Must have 0<= axis <= 1')

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

    # ndarray-like stats methods
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
            cols = self.cols()
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

        if len(num_cols) < len(self.cols()) and numeric_only:
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
            the_sum[the_count == 0] = NaN

        return Series(the_sum, index=axis_labels)

    def _get_numeric_columns(self):
        return [col for col in self.cols()
                if issubclass(self[col].dtype.type, np.number)]

    def _get_object_columns(self):
        return [col for col in self.cols() if self[col].dtype == np.object_]

    def cumsum(self, axis=0):
        """
        Return cumulative sum over requested axis as DataFrame

        Parameters
        ----------
        axis : {0, 1}
            0 for row-wise, 1 for column-wise

        Returns
        -------
        y : DataFrame
        """
        def get_cumsum(y):
            y = np.array(y)

            if not issubclass(y.dtype.type, np.int_):
                mask = isnull(y)
                y[mask] = 0
                result = y.cumsum()

                has_obs = (-mask).astype(int).cumsum() > 0
                result[-has_obs] = np.NaN
            else:
                result = y.cumsum()

            return result

        return self.apply(get_cumsum, axis=axis)

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
            theProd[theCount == 0] = NaN
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
        return self.apply(Series.min, axis=axis)

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
        return self.apply(Series.max, axis=axis)

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
            return Series(result, demeaned.cols())
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
        y = np.array(self.values, subok=True)
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
        y = np.array(self.values, subok=True)
        mask = np.isnan(y)
        count = (y.shape[axis] - mask.sum(axis)).astype(float)
        y[mask] = 0

        A = y.sum(axis) / count
        B = (y**2).sum(axis) / count  - A**2
        C = (y**3).sum(axis) / count - A**3 - 3*A*B

        theSkew = (np.sqrt((count**2-count))*C) / ((count-2)*np.sqrt(B)**3)

        return Series(theSkew, index=self._get_agg_axis(axis))

def _try_sort(iterable):
    listed = list(iterable)
    try:
        return sorted(listed)
    except Exception:
        return listed

def _extract_index(data, index):
    if len(data) == 0:
        if index is None:
            index = NULL_INDEX
    elif len(data) > 0 and index is None:
        # aggregate union of indices
        need_labels = False

        # this is pretty kludgy, better way?
        for v in data.values():
            if isinstance(v, Series):
                if index is None:
                    index = v.index
                elif need_labels:
                    raise Exception('Cannot mix Series / dict objects'
                                    ' with ndarray / sequence input')
                elif not index.equals(v.index):
                    index = index + v.index

            elif isinstance(v, dict):
                if index is None:
                    index = Index(_try_sort(v))
                elif need_labels:
                    raise Exception('Cannot mix Series / dict objects'
                                    ' with ndarray / sequence input')
                else:
                    index = index + Index(v.keys())

            else: # not dict-like, assign integer labels
                if index is not None and not need_labels:
                    raise Exception('Cannot mix Series / dict objects'
                                    ' with ndarray / sequence input')

                need_labels = True
                index = Index(np.arange(len(v)))

    if len(index) == 0 or index is None:
        index = NULL_INDEX

    if not isinstance(index, Index):
        index = Index(index)

    return index

def _default_index(n):
    if n == 0:
        return NULL_INDEX
    else:
        return np.arange(n)
