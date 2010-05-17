# pylint: disable-msg=E1101,E1103
# pylint: disable-msg=W0212,W0231,W0703,W0622

from cStringIO import StringIO
import operator
import sys

from numpy import NaN
import numpy as np

from pandas.core.common import _pickle_array, _unpickle_array, _pfixed
from pandas.core.daterange import DateRange
from pandas.core.index import Index, NULL_INDEX
from pandas.core.mixins import Picklable, Groupable
from pandas.core.series import Series
from pandas.lib.tseries import isnull, notnull
import pandas.core.datetools as datetools
import pandas.lib.tseries as tseries

#-------------------------------------------------------------------------------
# Factory helper methods

def arith_method(func, name):
    def f(self, other):
        return self._combineFunc(other, func)

    f.__name__ = name
    f.__doc__ = 'Wrapper for arithmetic method %s' % name

    return f

#-------------------------------------------------------------------------------
# DataFrame class

class DataFrame(Picklable, Groupable):
    """
    Homogenously indexed table with named columns, with intelligent
    arithmetic operations, slicing, reindexing, aggregation, etc. Can
    function interchangeably as a dictionary.

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
    Data contained within is COPIED from input arrays, this is to
    prevent silly behavior like altering the original arrays and
    having those changes reflected in the frame.

    See also
    --------
    DataMatrix: more efficient version of DataFrame for most operations

    Examples
    --------
        >>> d = {'col1' : ts1, 'col2' : ts2}
        >>> df = DataFrame(data=d, index=someIndex)
    """

    def __init__(self, data=None, index=None, columns=None, dtype=None):
        if isinstance(data, dict):
            self._series, self.index = self._initDict(data, index,
                                                      columns, dtype)

        elif isinstance(data, (np.ndarray, list)):
            self._series, self.index = self._initMatrix(data, index,
                                                        columns, dtype)
        elif isinstance(data, DataFrame):
            self._series = data._series.copy()
            self.index = data.index
        elif data is None:
            if index is None:
                index = NULL_INDEX

            self._series, self.index = {}, index

    def _initDict(self, data, index, columns, dtype):
        # pre-filter out columns if we passed it
        if columns is not None:
            colset = set(columns)
            data = dict((k, v) for k, v in data.iteritems() if k in colset)

        index = _extract_index(data, index)

        series = {}
        for k, v in data.iteritems():
            if isinstance(v, Series):
                # Forces alignment and copies data
                series[k] = v.reindex(index)
            else:
                if isinstance(v, dict):
                    v = [v.get(i, NaN) for i in index]

                try:
                    v = Series(v, dtype=dtype, index=index)
                except Exception:
                    v = Series(v, index=index)

                series[k] = v.copy()

        # add in any other columns we want to have (completeness)
        if columns is not None:
            for c in columns:
                if c not in series:
                    series[c] = Series.fromValue(np.NaN, index=index)

        return series, index

    def _initMatrix(self, data, index, columns, dtype):
        if not isinstance(data, np.ndarray):
            arr = np.array(data)
            if issubclass(arr.dtype.type, basestring):
                arr = np.array(data, dtype=object, copy=True)

            data = arr

        if data.ndim == 1:
            data = data.reshape((len(data), 1))
        elif data.ndim != 2:
            raise Exception('Must pass 2-d input!')

        if columns is None:
            raise Exception('Must pass column names')
        if index is None:
            raise Exception('Must pass index')

        N, K = data.shape

        if len(index) != N:
            raise Exception('Index length mismatch: %d vs. %d' %
                            (len(index), N))
        if len(columns) != K:
            raise Exception('Index length mismatch: %d vs. %d' %
                            (len(columns), K))

        data = dict([(idx, data[:, i]) for i, idx in enumerate(columns)])
        return self._initDict(data, index, columns, dtype)

    @property
    def _constructor(self):
        return DataFrame

    def __getstate__(self):
        series = dict((k, v.values()) for k, v in self.iteritems())
        index = _pickle_array(self.index)

        return series, index

    def __setstate__(self, state):
        series, idx = state

        index = _unpickle_array(idx)
        self._series = dict((k, Series(v, index=index))
                            for k, v in series.iteritems())
        self.index = index

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

    index = property(fget=_get_index, fset=_set_index)

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

        dataDict = dict((k, data[k]) for k in data.dtype.names)

        if indexField is not None:
            index = dataDict.pop(indexField)
        else:
            index = np.arange(len(data))

        return cls(dataDict, index=index)

    def toRecords(self):
        """
        Convert DataFrame to record array. Index will be put in the
        'index' field of the record array.

        Returns
        -------
        recarray
        """
        arrays = [self.index] + [self[c] for c in self.cols()]
        names = ['index'] + list(self.cols())

        return np.rec.fromarrays(arrays, names=names)

#-------------------------------------------------------------------------------
# Magic methods

    def __array__(self):
        return self.values

    def __array_wrap__(self, result):
        return DataFrame(result, index=self.index, columns=self.columns)

    def __nonzero__(self):
        return len(self._series) > 0 and len(self.index) > 0

    def __repr__(self):
        """
        Return a string representation for a particular DataFrame
        """
        buf = StringIO()
        if len(self.index) < 500 and len(self._series) < 10:
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
        if hasattr(value, '__iter__'):
            if isinstance(value, Series):
                cleanSeries = value.reindex(self.index)
            else:
                cleanSeries = Series(value, index=self.index)

            self._series[key] = cleanSeries
        # Scalar
        else:
            self._series[key] = Series.fromValue(value, index=self.index)

    def __delitem__(self, key):
        """
        Delete column from DataFrame
        """
        del self._series[key]

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
        return iter(self._series)

    def __len__(self):
        """
        Returns number of columns/Series inside
        """
        return len(self.index)

    def __contains__(self, key):
        """
        True if DataFrame has this column
        """
        return key in self._series

    __add__ = arith_method(operator.add, '__add__')
    __sub__ = arith_method(operator.sub, '__sub__')
    __mul__ = arith_method(operator.mul, '__mul__')
    __div__ = arith_method(operator.div, '__div__')
    __pow__ = arith_method(operator.pow, '__pow__')

    __radd__ = arith_method(operator.add, '__radd__')
    __rmul__ = arith_method(operator.mul, '__rmul__')
    __rsub__ = arith_method(lambda x, y: y - x, '__rsub__')
    __rdiv__ = arith_method(lambda x, y: y / x, '__rdiv__')
    __rpow__ = arith_method(lambda x, y: y ** x, '__rpow__')

    def __neg__(self):
        return self * -1

#-------------------------------------------------------------------------------
# Private / helper methods

    def _firstTimeWithValue(self):
        return self.index[self.count(1) > 0][0]

    def _lastTimeWithValue(self):
        return self.index[self.count(1) > 0][-1]

    def _combineFrame(self, other, func):
        newColumns = {}
        newIndex = self.index

        if self.index.equals(other.index):
            newIndex = self.index

            this = self
        else:
            newIndex = self.index + other.index

            this = self.reindex(newIndex)
            other = other.reindex(newIndex)

        if not self and not other:
            return DataFrame(index=newIndex)

        if not other:
            return self * NaN

        if not self:
            return other * NaN

        for col, series in this.iteritems():
            if col in other:
                newColumns[col] = func(series, other[col])
            else:
                newColumns[col] = series.fromValue(np.NaN, index=newIndex)

        for col, series in other.iteritems():
            if col not in self:
                newColumns[col] = series.fromValue(np.NaN, index=newIndex)

        return DataFrame(data=newColumns, index=newIndex)

    def _combineSeries(self, other, func):
        newColumns = {}
        newIndex = self.index

        if len(other) == 0:
            return self * NaN

        if len(self) == 0:
            # Ambiguous case
            return DataFrame(index=self.index, columns=self.cols())

        if self.index._allDates and other.index._allDates:
            if self.index.equals(other.index):
                newIndex = self.index
                this = self
            else:
                newIndex = self.index + other.index
                this = self.reindex(newIndex)
                other = other.reindex(newIndex)

            for col, series in this.iteritems():
                newColumns[col] = func(series, other)

            result = DataFrame(newColumns, index=newIndex)
        else:
            union = other.index.union(self.cols())
            intersection = other.index.intersection(self.cols())

            for col in intersection:
                newColumns[col] = func(self[col], other[col])

            result = DataFrame(newColumns, index=self.index)

            for col in (x for x in union if x not in intersection):
                result[col] = NaN

        return result

    def _combineFunc(self, other, func):
        """
        Combine DataFrame objects or a single DataFrame and a constant
        or other object using the supplied function.

        This is the core method used for all the 'magic' DataFrame methods.

        Parameters
        ----------
        other : constant, array, or DataFrame/Matrix
        func : function taking two arguments

        Examples
        --------
        frame._combineFunc(otherFrame, lambda x, y: x + y)
        """
        if isinstance(other, DataFrame):    # Another DataFrame
            return self._combineFrame(other, func)
        elif isinstance(other, Series):
            return self._combineSeries(other, func)
        else:
            newColumns = {}
            newIndex = self.index

            for col, series in self.iteritems():
                newColumns[col] = func(series, other)

            return DataFrame(data=newColumns, index=newIndex)

#-------------------------------------------------------------------------------
# Public methods

    def toCSV(self, path, nanRep='', cols=None, header=True,
              index=True, verbose=False):
        """
        Write the DataFrame to a CSV file
        """
        f = open(path, 'w')

        if cols is None:
            cols = self.cols()

        if header:
            if index:
                f.write(',' + ','.join([str(c) for c in cols]))
            else:
                f.write(','.join([str(c) for c in cols]))
            f.write('\n')

        for idx in self.index:
            if index:
                f.write(str(idx) + ',')
            for col in cols:
                val = self._series[col].get(idx)
                if isnull(val):
                    val = nanRep
                else:
                    val = str(val)
                f.write(val + ',')
            f.write('\n')

        f.close()

        if verbose: # pragma: no cover
            print 'CSV file written successfully: %s' % path

    def toDataMatrix(self):
        from pandas.core.matrix import DataMatrix

        return DataMatrix(self._series, index=self.index)

    def toString(self, buffer=sys.stdout, columns=None, colSpace=15,
                 nanRep='NaN', formatters=None, float_format=None):
        """Output a tab-separated version of this DataFrame"""
        series = self._series

        if columns is None:
            columns = _try_sort(series)
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
        columns = _try_sort(self.cols())
        space = max([len(str(k)) for k in columns]) + 4
        for k in columns:
            out = _pfixed(k, space)
            N = notnull(series[k]).sum()
            out += '%d  non-null values' % N
            print >> buffer, out

    def rows(self):
        """Alias for the frame's index"""
        return self.index

    def cols(self):
        """Return sorted list of frame's columns"""
        return _try_sort(self._series)

    # For DataMatrix compatibility
    columns = property(lambda self: Index(self.cols()))

    def iteritems(self):
        """Iterator over (column, series) pairs"""
        return self._series.iteritems()

    def append(self, otherFrame):
        """
        Append columns of otherFrame to end of this frame's columns and index.

        Columns not in this frame are added as new columns.
        """
        newIndex = np.concatenate((self.index, otherFrame.index))
        newValues = {}

        for column, series in self.iteritems():
            if column in otherFrame:
                newValues[column] = series.append(otherFrame[column])
            else:
                newValues[column] = series

        for column, series in otherFrame.iteritems():
            if column not in self:
                newValues[column] = series

        return DataFrame(data=newValues, index=newIndex)

    def asfreq(self, freq, fillMethod=None):
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
        if isinstance(freq, datetools.DateOffset):
            dateRange = DateRange(self.index[0], self.index[-1], offset=freq)
        else:
            dateRange = DateRange(self.index[0], self.index[-1], timeRule=freq)

        return self.reindex(dateRange, fillMethod=fillMethod)

    def asMatrix(self, columns=None):
        """
        Convert the frame to its Numpy-array matrix representation

        Columns are presented in sorted order unless a specific list
        of columns is provided.
        """
        if columns is None:
            return np.array([self[col] for col in self.cols()]).T
        else:
            return np.array([self[col] for col in columns]).T

    # For DataMatrix compatibility
    values = property(asMatrix)

    def copy(self):
        """
        Make a deep copy of this frame
        """
        newFrame = DataFrame(index=self.index)
        newFrame._series = dict((k, v.copy()) for k, v in self.iteritems())
        return newFrame

    def corr(self):
        """
        Compute correlation of columns. Returns DataFrame of result matrix.

        In the presence of NaN values in any column, tries to compute the
        pairwise correlation between the column and the other columns.
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

        newIndex = self.index[theCount != 0]
        return self.reindex(newIndex)

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

        newIndex = self.index[theCount >= minObs]
        return self.reindex(newIndex)

    def fill(self, value=None, method='pad'):
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
            filledSeries = series.fill(method=method, value=value)
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

    def getXS(self, key):
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
        from pandas.core.panel import _slow_pivot

        return _slow_pivot(self[index], self[columns], self[values])

    def reindex(self, index=None, columns=None, fillMethod=None):
        """
        Reindex data inside, optionally filling according to some rule.

        Parameters
        ----------
        index : array-like, optional
            preferably an Index object (to avoid duplicating data)
        columns : array-like, optional
        fillMethod : {'backfill', 'pad', None}
            Method to use for filling data holes using the index

        Returns
        -------
        y : same type as calling instance
        """
        if fillMethod:
            fillMethod = fillMethod.upper()

        frame = self

        if index is not None:
            frame = frame._reindex_index(index, fillMethod)

        if columns is not None:
            frame = frame._reindex_columns(columns)

        return frame

    def _reindex_index(self, index, method):
        if self.index.equals(index):
            return self.copy()

        if not isinstance(index, Index):
            index = Index(index)

        if len(self.index) == 0:
            return DataFrame(index=index)

        fillVec, mask = tseries.getFillVec(self.index, index,
                                           self.index.indexMap,
                                           index.indexMap, method)

        # Maybe this is a bit much? Wish I had unit tests...
        typeHierarchy = [
            (float, float),
            (int, float),
            (bool, np.bool_),
            (np.bool_, np.bool_),
            (basestring, object),
            (object, object)
        ]

        missingValue = {
            float  : NaN,
            object : NaN,
            np.bool_ : False
        }

        newSeries = {}
        for col, series in self.iteritems():
            series = series.view(np.ndarray)
            for klass, dest in typeHierarchy:
                if issubclass(series.dtype.type, klass):
                    new = series.take(fillVec).astype(dest)
                    new[-mask] = missingValue[dest]
                    newSeries[col] = new
                    break

        return DataFrame(newSeries, index=index)

    def _reindex_columns(self, columns):
        result = DataFrame(index=self.index)

        if len(columns) == 0:
            return result

        for col in columns:
            if col in self:
                result[col] = self[col]
            else:
                result[col] = NaN

        return result

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
        self._series = dict((mapper(k), v) for k, v in self._series.iteritems())

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

        selfM = np.array([self[col] for col in self.cols()], dtype=theDtype)
        idxMap = self.index.indexMap
        return DataFrame(data=dict([(idx, selfM[:, idxMap[idx]])
                                        for idx in self.index]),
                                        index=Index(self.cols()))

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

        N = len(self)

        if offset is None:
            newIndex = self.index

            indexer = np.zeros(N, dtype=int)
            if periods > 0:
                indexer[periods:] = np.arange(N - periods)
                def do_shift(series):
                    values = np.asarray(series).take(indexer)
                    values[:periods] = NaN
                    return values

            else:
                indexer[:periods] = np.arange(-periods, N)
                def do_shift(series):
                    values = np.asarray(series).take(indexer)
                    values[periods:] = NaN
                    return values

            newValues = dict([(col, do_shift(series))
                              for col, series in self.iteritems()])
        else:
            offset = periods * offset
            newIndex = Index([idx + offset for idx in self.index])
            newValues = dict([(col, np.asarray(series))
                               for col, series in self.iteritems()])

        return DataFrame(data=newValues, index=newIndex)

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

        Note
        ----
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
            return self._constructor(data=results, index=target.index)
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
        return DataFrame(data=results, index=self.index)

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

    def sortUp(self, column=None):
        """
        Sort DataFrame in ascending order according to specified column,
        otherwise by the index.
        """
        if column:
            series = self[column].order(missingAtEnd=True)
            return self.reindex(series.index)
        else:
            idx = np.array(np.argsort(self.index))
            newIndex = self.index[idx.astype(int)]
            return self.reindex(newIndex)

    def sortDown(self, column=None):
        """
        Sort DataFrame in ascending order according to specified column,
        otherwise by the index.
        """
        if column:
            series = self[column].order(missingAtEnd=False)
            return self.reindex(series.index[::-1])
        else:
            idx = np.array(np.argsort(self.index))
            idx = idx[::-1]  # Reverses array
            newIndex = self.index[idx.astype(int)]
            return self.reindex(newIndex)

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
            return self

        if not self:
            return other

        if self.index is not other.index:
            unionIndex = self.index + other.index
            frame = self.reindex(unionIndex)
            other = other.reindex(unionIndex)
        else:
            unionIndex = self.index
            frame = self

        do_fill = fill_value is not None
        unionCols = _try_sort(set(frame.cols() + other.cols()))

        result = {}
        for col in unionCols:
            if col in frame and col in other:
                series = frame[col].values()
                otherSeries = other[col].values()

                if do_fill:
                    this_mask = isnull(series)
                    other_mask = isnull(otherSeries)
                    series = series.copy()
                    otherSeries = otherSeries.copy()
                    series[this_mask] = fill_value
                    otherSeries[other_mask] = fill_value

                result[col] = func(series, otherSeries)

                if do_fill:
                    result[col][this_mask & other_mask] = np.NaN

            elif col in frame:
                result[col] = frame[col]
            elif col in other:
                result[col] = other[col]

        return self._constructor(result, index=unionIndex)

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
            default: 'left' for joining on index, None otherwise
            How to handle indexes of the two objects.
              * left: use calling frame's index
              * right: use input frame's index
              * outer: form union of indexes
              * inner: use intersection of indexes

        Examples
        --------
        This frame         Other frame
            c1                 q1
        a   1              0   v1
        b   0              1   v2
        c   1
        d   0
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

        fillVec, mask = tseries.getMergeVec(self[on], other.index.indexMap)

        newSeries = {}

        for col, series in other.iteritems():
            arr = series.view(np.ndarray).take(fillVec)
            arr[-mask] = NaN

            newSeries[col] = arr

        newSeries.update(self._series)

        return self._constructor(newSeries, index=self.index)

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

    def plot(self, kind='line', **kwds): # pragma: no cover
        """
        Plot the DataFrame's series with the index on the x-axis using
        matplotlib / pylab.

        Parameters
        ----------
        kind : {'line', 'bar', 'hist'}
            Default: line for TimeSeries, hist for Series

        kwds : other plotting keyword arguments

        Note
        ----
        This method doesn't make much sense for cross-sections,
        and will error.
        """
        from pylab import plot

        for col in _try_sort(self.columns):
            plot(self.index, self[col].values(), label=col)

    def _get_agg_axis(self, axis_num):
        if axis_num == 0:
            return self.columns
        elif axis_num == 1:
            return self.index
        else:
            raise Exception('Must have 0<= axis <= 1')

    def cap(self, threshold):
        """
        Trim values at threshold

        Returns
        -------
        DataFrame
        """
        return self.apply(lambda x: x.cap(threshold))

    def floor(self, threshold):
        """
        Trim values below threshold

        Returns
        -------
        DataFrame
        """
        return self.apply(lambda x: x.floor(threshold))

    # ndarray-like stats methods
    def count(self, axis=0):
        """
        Return array or Series of # observations over requested axis.

        Parameters
        ----------
        axis : {0, 1}
            0 for row-wise, 1 for column-wise

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

    def sum(self, axis=0):
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
        self:
            c1  c2
        a   1   0
        b   0   2
        c   3   0
        d   0   4

        >>>self.sum(axis=0)
        Series:
        c1: 4
        c2: 6
        """
        try:
            y = self.values.copy()
            mask = np.isfinite(y)
            if not issubclass(y.dtype.type, np.int_):
                y[-mask] = 0
            theSum = y.sum(axis)
            theCount = mask.sum(axis)
            theSum[theCount == 0] = NaN

        except Exception:
            theSum = self.apply(np.sum, axis=axis)

        axis_labels = self._get_agg_axis(axis)
        return Series(theSum, index=axis_labels)

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
        summed = self.sum(axis)
        count = self.count(axis).astype(float)

        return summed / count.reindex(summed.index)

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
            return tseries.median(arr[notnull(arr)])

        if axis == 0:
            med = [f(self[col].values()) for col in self.columns]
            return Series(med, index=self.columns)
        elif axis == 1:
            med = [f(self.getXS(k).values()) for k in self.index]
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
        if axis == 0:
            med = [np.min(self[col].valid()) for col in self.columns]
            return Series(med, index=self.columns)
        elif axis == 1:
            med = [np.min(self.getXS(k).valid()) for k in self.index]
            return Series(med, index=self.index)
        else:
            raise Exception('Must have 0<= axis <= 1')

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
        if axis == 0:
            med = [np.max(self[col].valid()) for col in self.columns]
            return Series(med, index=self.columns)
        elif axis == 1:
            med = [np.max(self.getXS(k).valid()) for k in self.index]
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

