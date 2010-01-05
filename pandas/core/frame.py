# pylint: disable-msg=E1101
# pylint: disable-msg=E1103
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
    data : dict
        Mapping of column name --> array or Series/TimeSeries objects
    index : array-like, optional
        Specific index to use for the Frame, Series will be conformed
        to this if you provide it. If not input, index will be
        inferred from input Series
    columns : array-like, optional

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
    def __init__(self, data=None, index=None, columns=None):
        self._series = {}
        if data is not None and len(data) > 0:
            if index is None:
                s = data.values()[0]
                if isinstance(s, Series):
                    self.index = s.index
                else:
                    self.index = np.arange(len(s))
            else:
                self.index = index

            for k, v in data.iteritems():
                if isinstance(v, Series):
                    # Forces homogeneity and copies data
                    self._series[k] = v.reindex(self.index)
                else:
                    # Copies data and checks length
                    self._series[k] = Series(v, index=self.index)

        elif index is not None:
            self.index = index
        else:
            self.index = NULL_INDEX

    def __getstate__(self):
        series = dict((k, v.values()) for k, v in self.iteritems())
        index = _pickle_array(self.index)

        return series, index

    def __setstate__(self, state):
        series, idx = state

        self.index = index = _unpickle_array(idx)
        self._series = dict((k, Series(v, index=index))
                            for k, v in series.iteritems())

    _index = None
    def _set_index(self, index):
        if isinstance(index, Index):
            self._index = index
        else:
            self._index = Index(index)

    def _get_index(self):
        return self._index

    index = property(fget=_get_index, fset=_set_index)

    # Alternate constructors
    @classmethod
    def fromDict(cls, inputDict=None, castFloat=True, **kwds):
        """
        Convert a two-level tree representation of a series or time series
        to a DataFrame.

        tree is structured as:
            {'col1' : {
                idx1 : value1,
                ...
                idxn : valueN
                    },
            ...}
        e.g. tree['returns'][curDate] --> return value for curDate

        Parameters
        ----------
        input : dict object
            Keys become column names of returned frame
        kwds : optionally provide arguments as keywords

        Returns
        -------
        DataFrame

        Examples
        --------
        df1 = DataFrame.fromDict(myDict)
        df2 = DataFrame.fromDict(A=seriesA, B=seriesB)
        """
        if inputDict is None:
            inputDict = {}
        else:
            if not hasattr(inputDict, 'iteritems'):
                raise Exception('Input must be a dict or dict-like!')
            inputDict = inputDict.copy()

        inputDict.update(kwds)

        if len(inputDict) == 0:
            return DataFrame(index=NULL_INDEX)
        elif len(inputDict) == 1:
            index = inputDict.values()[0].keys()
            if not isinstance(index, Index):
                index = Index(sorted(index))
        else:
            # GET set of indices
            indices = set([])
            for key, branch in inputDict.iteritems():
                indices = indices | set(branch.keys())
            index = Index(sorted(indices))

        columns = {}
        for key, branch in inputDict.iteritems():
            if isinstance(branch, Series):
                tmp = branch.reindex(index)
            else:
                tmp = [branch.get(i, NaN) for i in index]

            try:
                if castFloat:
                    columns[key] = Series(tmp, dtype=float, index=index)
                else:
                    columns[key] = Series(tmp, index=index)
            except Exception:
                columns[key] = Series(tmp, index=index)

        return DataFrame(data=columns, index=index)

    def toDict(self):
        """
        Simpler pseudo-inverse operation of DataFrame.fromDict, NaN
        values will be included in the resulting dict-tree.

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
        if data.dtype.type != np.void:
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

    @classmethod
    def fromMatrix(cls, mat, colNames, rowNames):
        """
        Convert input matrix to DataFrame given column and row names (index)

        Parameters
        ----------
        mat : ndarray
            Dimension T x N
        colNames : iterable
            Dimension N
        rowNames : iterable
            Dimension T

        Returns
        -------
        DataFrame
        """
        rows, cols = mat.shape
        try:
            assert(rows == len(rowNames))
            assert(cols == len(colNames))
        except AssertionError:
            raise Exception('Dimensions do not match: %s, %s, %s' %
                            (mat.shape, len(rowNames), len(colNames)))

        index = Index(rowNames)
        colIndex = Index(colNames)

        idxMap = colIndex.indexMap

        data = dict([(idx, mat[:, idxMap[idx]]) for idx in colIndex])
        return DataFrame(data=data, index=index)

#-------------------------------------------------------------------------------
# Magic methods

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

    def _firstTimeWithNValues(self):
        # Need to test this!
        N = len(self._series)
        theCount = np.isfinite(self.asMatrix()).sum(1)
        selector = (theCount == N)
        if not selector.any():
            raise Exception('No time has %d values!' % N)

        return self.index[selector][0]

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

        if self.index._allDates and other.index._allDates:
            if not self:
                return DataFrame(index=other.index)

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
        newColumns = {}
        newIndex = self.index

        if isinstance(other, DataFrame):    # Another DataFrame
            return self._combineFrame(other, func)
        elif isinstance(other, Series):
            return self._combineSeries(other, func)
        else:
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

        cols = self.cols() if cols is None else cols

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

    def toString(self, buffer=sys.stdout, verbose=False,
                 columns=None, colSpace=15, nanRep='NaN',
                 formatters=None, float_format=None):
        """Output a tab-separated version of this DataFrame"""
        series = self._series

        if columns is None:
            columns = sorted(series.keys())
        else:
            columns = [c for c in columns if c in self]

        formatters = formatters or {}

        # TODO
        ident = lambda x: x

        if len(columns) == 0 or len(self.index) == 0:
            print >> buffer, 'Empty DataFrame'
            print >> buffer, repr(self.index)
        else:
            idxSpace = max([len(str(idx)) for idx in self.index]) + 4
            head = _pfixed('', idxSpace)
            if verbose:
                colSpace = max([len(c) for c in columns]) + 4

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
        print >> buffer, 'Index: %s entries, %s to %s' % (len(self.index),
                                                          min(self.index),
                                                          max(self.index))
        print >> buffer, 'Data columns:'

        if len(self._series) == 0:
            print >> buffer, 'DataFrame is empty!'
            return

        series = self._series
        columns = sorted(self.cols())
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
        return sorted(self._series.keys())

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

        return self.fromMatrix(correl, cols, cols)

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
        newIndex = self.index[self.count(1) != 0]
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
        N = len(self._series)

        if specificColumns:
            colSet = set(specificColumns)
            intersection = set(self.cols()) & colSet

            N = len(intersection)

            filtered = self.filterItems(intersection)
            theCount = filtered.count(axis=1)
        else:
            theCount = self.count(axis=1)

        if minObs is None:
            minObs = N

        newIndex = self.index[theCount >= N]
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

    def getTS(self, colName=None, fromDate=None, toDate=None, nPeriods=None):
        """
        Return a DataFrame / TimeSeries corresponding to given arguments

        Parameters
        ----------
        colName : particular column name requested
        fromDate : datetime
        toDate : datetime
        nPeriods : int/float

        Note
        ----
        Error thrown if all of fromDate, toDate, nPeriods specified.
        """
        beg_slice, end_slice = self._getIndices(fromDate, toDate)

        if nPeriods:
            if fromDate and toDate:
                raise Exception('fromDate/toDate, toDate/nPeriods, ' + \
                                ' fromDate/nPeriods are mutually exclusive')
            elif fromDate:
                end_slice = min(len(self), beg_slice + nPeriods)
            elif toDate:
                beg_slice = max(0, end_slice - nPeriods)
            else:
                raise Exception('Not enough arguments provided to getTS')

        dateRange = self.index[beg_slice:end_slice]

        if colName:
            return self[colName].reindex(dateRange)
        else:
            return self.reindex(dateRange)

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
        before = arg_before = datetools.to_datetime(before)
        after = arg_after = datetools.to_datetime(after)

        if before is None:
            before = self.index[0]
        elif before not in self.index:
            loc = self.index.searchsorted(before, side='left')
            before = self.index[loc]

        if after is None:
            after = self.index[-1]
        elif after not in self.index:
            loc = self.index.searchsorted(after, side='right') - 1
            loc = loc if loc < len(self.index) else -1
            after = self.index[loc]

        beg_slice = self.index.indexMap[before]
        end_slice = self.index.indexMap[after] + 1

        return beg_slice, end_slice

    def getXS(self, key, subset=None):
        """
        Returns a row from the DataFrame as a Series object.

        Parameters
        ----------
        key : some index contained in the index
        subset : iterable (list, array, set, etc.), optional
            columns to be included

        Note
        ----
        Will try to return a TimeSeries if the columns are dates.
        """
        subset = list(set(subset) if subset else set(self._series.keys()))
        subset.sort()

        rowValues = [self._series[k].get(key) for k in subset]

        if len(set(map(type, rowValues))) > 1:
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
        fillMethod = fillMethod.upper() if fillMethod else ''

        if fillMethod not in ['BACKFILL', 'PAD', '']:
            raise Exception("Don't recognize fillMethod: %s" % fillMethod)

        frame = self

        if index is not None:
            frame = frame._reindex_index(index, fillMethod)

        if columns is not None:
            frame = frame._reindex_columns(columns)

        return frame

    def _reindex_index(self, index, method):
        if self.index.equals(index):
            return self.copy()

        if len(index) == 0:
            return DataFrame(index=NULL_INDEX)

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
            object : None,
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

        if timeRule is not None and offset is None:
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
        elif axis == 1:
            target = self.T

        results = dict([(k, func(target[k])) for k in target.columns])

        if hasattr(results.values()[0], '__iter__'):
            return DataFrame(data=results, index=target.index)
        else:
            return Series.fromDict(results)

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
        TODO
        """
        import re

        if items is not None:
            data = dict([(r, self[r]) for r in items if r in self])
            return DataFrame(data=data, index=self.index)
        elif like:
            columns = [c for c in self.cols() if like in c]
            return self.reindex(columns=columns)
        elif regex:
            matcher = re.compile(regex)
            columns = [c for c in self.cols() if matcher.match(c)]
            return self.reindex(columns=columns)

    def filterItems(self, items):
        """
        Restrict frame's columns to input set of items.

        Parameters
        ----------
        items : list-like
            List of columns to restrict to (must not all be present)

        Returns
        -------
        DataFrame with filtered columns
        """
        return self.filter(items=items)

    def filterLike(self, arg):
        """
        Filter to columns partially matching the import argument.

        Keep columns where "arg in col == True"

        Parameter
        ---------
        arg : string

        Return
        ------
        DataFrame with matching columns
        """
        return self.filter(like=arg)

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
        unionCols = sorted(set(frame.cols() + other.cols()))

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

        return DataFrame(result, index = unionIndex)

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

    def outerJoin(self, *frames):
        unionIndex = self.index
        for frame in frames:
            unionIndex  = unionIndex + frame.index

        joined = self.reindex(unionIndex)
        for frame in frames:
            frame = frame.reindex(unionIndex)
            for col, series in frame.iteritems():
                if col in joined:
                    raise Exception('Overlapping columns!')
                joined[col] = series

        return joined

    def leftJoin(self, *frames):
        """
        Insert columns of input DataFrames / dicts into this one.

        Columns must not overlap. Returns a copy.

        Parameters
        ----------
        *frames : list-like
            List of frames (DataMatrix or DataFrame) as function arguments

        Returns
        -------
        DataFrame
        """
        joined = self.copy()

        for frame in frames:
            for col, series in frame.iteritems():
                if col in joined:
                    raise Exception('Overlapping columns!')
                joined[col] = series.copy()

        return joined

    def merge(self, other, on=None):
        """
        Merge DataFrame or DataMatrix with this one on some many-to-one index

        Parameters
        ----------
        other : DataFrame
            Index should be similar to one of the columns in this one
        on : string
            Column name to use

        Examples
        --------
        This frame         Other frame
            c1                 q1
        a   1              0   v1
        b   0              1   v2
        c   1
        d   0
        """
        if len(other.index) == 0:
            return self

        if on not in self:
            raise Exception('%s column not contained in this frame!' % on)

        # Check for column overlap
        overlap = set(self.cols()) & set(other.cols())

        if any(overlap):
            raise Exception('Columns overlap: %s' % sorted(overlap))

        fillVec, mask = tseries.getMergeVec(self[on], other.index.indexMap)

        newSeries = {}

        for col, series in other.iteritems():
            arr = series.view(np.ndarray).take(fillVec)
            arr[-mask] = NaN

            newSeries[col] = arr

        newSeries.update(self._series)

        return DataFrame(newSeries, index=self.index)

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

        for col in sorted(self.columns):
            s = self[col]
            plot(s.index, s, label=col)

    def _get_axis(self, axis_num):
        if axis_num == 0:
            return self.columns
        elif axis_num == 1:
            return self.index
        else:
            raise Exception('Must have 0<= axis <= 1')

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
            theCount = np.isfinite(self.values).sum(axis)
        except Exception:
            f = lambda s: notnull(s).sum()
            theCount = self.apply(f, axis=axis)

        return Series(theCount, index=self._get_axis(axis))

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
        y = np.array(self.values, subok=True)

        try:
            if not issubclass(y.dtype.type, np.int_):
                y[np.isnan(y)] = 0
            theSum = y.sum(axis)
            theCount = self.count(axis)
            theSum[theCount == 0] = NaN
        except Exception:
            theSum = self.apply(np.sum, axis=axis)

        return Series(theSum, index=self._get_axis(axis))

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
                y[np.isnan(y)] = 0
            return y.cumsum()

        return self.apply(get_cumsum, axis=axis)

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

        return Series(theProd, index=self._get_axis(axis))

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
        return self.sum(axis) / self.count(axis).values().astype(float)

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
        if axis == 0:
            med = [np.median(self[col].valid()) for col in self.columns]
            return Series(med, index=self.columns)
        elif axis == 1:
            med = [np.median(self.getXS(k).valid()) for k in self.index]
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

        return Series(theVar, index=self._get_axis(axis))

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

        return Series(theSkew, index=self._get_axis(axis))

