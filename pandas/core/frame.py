# pylint: disable-msg=E1101
# pylint: disable-msg=E1103
# pylint: disable-msg=W0212

import operator

from numpy import NaN
import numpy as np

from pandas.core.daterange import DateRange
from pandas.core.datetools import DateOffset, to_datetime
from pandas.core.index import Index, NULL_INDEX
from pandas.core.mixins import Picklable, Groupable
from pandas.core.series import Series, remove_na
from pandas.lib.tseries import isnull, notnull
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
    data: dict
        Mapping of column name --> array or Series/TimeSeries objects
    index: array-like
        Specific index to use for the Frame, Series will be conformed to this
        if you provide it.

    Notes
    -----
    Data contained within is COPIED from input arrays, this is to
    prevent silly behavior like altering the original arrays and
    having those changes reflected in the frame.

    See also
    --------
    DataMatrix: more efficient version of DataFrame for most operations

    Example usage
    -------------
        >>> d = {'col1' : ts1, 'col2' : ts2}
        >>> df = DataFrame(data=d, index=someIndex)
    """
    def __init__(self, data = None, index = None):
        series = data
        self._series = {}
        if series is not None and len(series) > 0:
            if index is None:
                for s in series.values(): break
                if hasattr(s, 'index'):
                    self.index = s.index
                else:
                    self.index = Index(np.arange(len(s)))
            else:
                if isinstance(index, Index):
                    self.index = index
                else:
                    self.index = Index(index)
            for k, v in series.iteritems():
                if isinstance(v, Series):
                    self._series[k] = v.reindex(self.index) # Forces homogoneity
                else:
                    assert(len(v) == len(self.index))
                    s = Series(v, index=self.index)
                    self._series[k] = s
        elif index is not None:
            if isinstance(index, Index):
                self.index = index
            else:
                self.index = Index(index)
        else:
            raise Exception('DataFrame constructor not properly called!')

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
        input: dict object
            Keys become column names of returned frame
        kwds: optionally provide arguments as keywords

        Returns
        -------
        DataFrame

        Example
        -------
        df1 = DataFrame.fromDict(myDict)
        df2 = DataFrame.fromDict(A=seriesA, B=seriesB)
        """
        if inputDict is None:
            inputDict = {}
        else:
            if not hasattr(inputDict, 'iteritems'):
                raise Exception('Input must be dict or dict-like!')
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
            except Exception, e:
                columns[key] = Series(tmp, index=index)

        return DataFrame(data=columns, index=index)

    @classmethod
    def fromRecords(cls, data, indexField=None):
        """
        Convert structured or record ndarray to DataFrame

        Parameters
        ----------
        input: NumPy structured array

        Returns
        -------
        DataFrame
        """
        # Dtype when you have records
        if data.dtype.type != np.void:
            raise Exception('Input was not a structured array!')

        columns = data.dtype.names
        dataDict = dict((k, data[k]) for k in columns)

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
        mat: ndarray
            Dimension T x N
        colNames: iterable
            Dimension N
        rowNames: iterable
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

        return DataFrame(data = dict([(idx, mat[:, idxMap[idx]])
                                      for idx in colIndex]),
                         index = index)

    @classmethod
    def load(cls, baseFile):
        """
        Load DataFrame from file.

        Parameters
        ----------
        baseFile: string
            Filename base where index/values are stored.
            e.g. baseFile='myfile' --> 'myfile_index.npy', 'myfile_values.npy'

        Returns
        -------
        DataFrame
        """
        cacheLoad = np.load(baseFile + '.npz')

        values = cacheLoad['v']
        index = Index(cacheLoad['i'])
        cols = cacheLoad['c']

        return cls.fromMatrix(values, cols, index)

    def save(self, baseFile):
        """
        Write DataFrame efficiently to file using NumPy serialization,
        which is easily 100x faster than cPickle.

        Note
        ----
        Saves data to 3 files, one for index, columns, and values matrix.
        """
        np.savez(baseFile, i=self.index, v=self.values, c=self.columns)

#-------------------------------------------------------------------------------
# Magic methods

    def __nonzero__(self):
        return len(self._series) > 0 and len(self.index) > 0

    def __repr__(self):
        """
        Return a string representation for a particular DataFrame
        """
        if len(self.index) < 1000 and len(self._series) < 10:
            return self.toString(to_stdout=False)
        else:
            output = str(self.__class__) + '\n'
            return output + self.info(to_stdout=False)

    def __getitem__(self, item):
        """
        Retrieve column or slice from DataFrame
        """
        try:
            return self._series[item]
        except (TypeError, KeyError):
            if isinstance(item, slice):
                start, stop = item.start, item.stop
                start = 0 if start is None else start
                stop = len(self) if stop is None else stop
                if start < 0:
                    start += len(self)
                if stop < 0:
                    stop += len(self)

                dateRange = self.index[start:stop]
                newColumns = {}
                for col, series in self.iteritems():
                    newColumns[col] = series[start:stop]
                return DataFrame(data=newColumns, index=dateRange)
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
        try:
            if hasattr(value, '__iter__'):
                if hasattr(value, 'reindex'):
                    cleanSeries = value.reindex(self.index)
                else:
                    cleanSeries = Series(value, index=self.index)
                assert(len(cleanSeries) == len(self.index))
                self._series[key] = cleanSeries

            # Scalar
            else:
                self._series[key] = Series.fromValue(value, index=self.index)
        except AssertionError:
            raise
        except Exception, e:
            raise Exception('Could not put key, value pair in Frame!')

    def __delitem__(self, key):
        """
        Delete column from DataFrame (only deletes the reference)
        """
        r = self._series.pop(key, None)

    def __iter__(self):
        """
        Iterate over columns of the frame.
        """
        return self._series.__iter__()

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
        mycopy = self.copy()
        myseries = mycopy._series
        for col in myseries:
            mycopy[col] = -myseries[col]
        return mycopy

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
        else:
            newIndex = self.index + other.index

        if not self and not other:
            return DataFrame(index=newIndex)

        if not other:
            return self * NaN

        if not self:
            return other * NaN

        for col, series in self.iteritems():
            if col in other:
                newSeries = func(series, other[col])
                newColumns[col] = newSeries.reindex(newIndex)
            else:
                cls = series.__class__
                newColumns[col] = cls(np.repeat(NaN, len(newIndex)),
                                          index=newIndex)
        for col, series in other.iteritems():
            if col not in self:
                cls = series.__class__
                newColumns[col] = cls(np.repeat(NaN, len(newIndex)),
                                      index=newIndex)

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
            else:
                newIndex = self.index + other.index

            other = other.reindex(newIndex)
            for col, series in self.iteritems():
                newColumns[col] = func(series.reindex(newIndex), other)

        else:
            for col, series in self.iteritems():
                if col in other.index:
                    newColumns[col] = func(series, other[col])
                else:
                    cls = series.__class__
                    newColumns[col] = cls(np.repeat(NaN, len(self.index)),
                                          index=self.index)

        return DataFrame(data=newColumns, index=newIndex)

    def _combineFunc(self, other, func):
        """
        Combine DataFrame objects or a single DataFrame and a constant
        or other object using the supplied function.

        This is the core method used for all the 'magic' DataFrame methods.

        Parameters
        ----------
        other: constant, array, or DataFrame/Matrix
        func: function taking two arguments

        Example
        -------
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

    def toCSV(self, path=None, nanRep='', cols=None, inclHeader=True,
              inclIndex=True, verbose=False):
        """
        Write the DataFrame to a CSV file
        """
        if path is None:
            import sys
            f = sys.stdout
        else:
            f = open(path, 'w')
        cols = self.cols() if cols is None else cols
        if inclHeader:
            if inclIndex:
                f.write(',' + ','.join([str(c) for c in cols]))
            else:
                f.write(','.join([str(c) for c in cols]))
            f.write('\n')
        for idx in self.index:
            if inclIndex:
                f.write(str(idx) + ',')
            for col in cols:
                val = self._series[col].get(idx)
                if isnull(val):
                    val = nanRep
                else:
                    val = str(val)
                f.write(val + ',')
            f.write('\n')
        if path is not None:
            f.close()

        if verbose:
            print 'CSV file written successfully: %s' % path

    def toDict(self):
        """
        Simpler pseudo-inverse operation of dictToDataFrame, NaN values will be
        included in the resulting dict-tree.

        Return
        ------
        nested dict mapping: {column -> index -> value}
        """
        tree = {}
        for col, series in self.iteritems():
            tree[col] = branch = {}
            for i in self.index:
                branch[i] = series[i]
        return tree

    def toDataMatrix(self):
        from pandas.core.matrix import DataMatrix

        return DataMatrix(self._series, index=self.index)

    def toString(self, to_stdout=True, verbose=False, colSpace=15, nanRep=None):
        """Output a tab-separated version of this DataFrame"""
        series = self._series
        skeys = sorted(series.keys())
        if len(skeys) == 0 or len(self.index) == 0:
            output = 'Empty DataFrame\n'
            output += self.index.__repr__()
        else:
            idxSpace = max([len(str(idx)) for idx in self.index]) + 4
            head = _pfixed('', idxSpace)
            if verbose:
                colSpace = max([len(c) for c in self.columns]) + 4
            for h in skeys:
                head += _pfixed(h, colSpace)
            output = head + '\n'
            for idx in self.index:
                ot = _pfixed(idx, idxSpace)
                for k in skeys:
                    ot += _pfixed(series[k][idx], colSpace, nanRep=nanRep)
                output += ot + '\n'
        if to_stdout:
            print output
        else:
            return output

    def info(self, to_stdout=True):
        """Concise summary of a DataFrame, used in __repr__ when very large."""
        if len(self._series) == 0:
            output = 'DataFrame is empty!\n'
            output += self.index.__repr__()
            return output

        output = 'Index: %s entries, %s to %s\n' % (len(self.index),
                                                    min(self.index),
                                                    max(self.index))
        output += 'Columns:\n'
        series = self._series
        skeys = sorted(self.cols())
        space = max([len(str(k)) for k in skeys]) + 4
        for k in skeys:
            out = _pfixed(k, space)
            N = notnull(series[k]).sum()
            out += '%d  non-null values\n' % N
            output += out
        if to_stdout:
            print output
        else:
            return output

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

    def asfreq(self, freq, fillMethod = None):
        """
        Convert all TimeSeries inside to specified frequency using DateOffset
        objects. Optionally provide fill method to pad/backfill/interpolate
        missing values.

        Parameters
        ----------
        offset: DateOffset object, or string in {'WEEKDAY', 'EOM'}
            DateOffset object or subclass (e.g. monthEnd)

        fillMethod: {'backfill', 'pad', 'interpolate', None}
                    Method to use for filling holes in new inde
        """
        if not isinstance(freq, datetools.DateOffset):
            raise Exception('Must pass DateOffset!')

        dateRange = DateRange(self.index[0], self.index[-1], offset=freq)

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
            bada = np.isnan(A)
            if np.any(bada):
                for j, B in enumerate(mat):
                    commonVec = (- bada) & (- np.isnan(B))
                    if any(commonVec):
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
        specificColumns: list-like, optional keyword
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
        minObs: int or None (default)
           Instead of requiring all the columns to have observations, require
           only minObs observations
        specificColumns: list-like, optional keyword
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

            #if len(cols) < N:
                #diff = str(set(specificColumns) - set(cols))
                #raise Exception('Missing columns: %s' % diff)

            filtered = self.filterItems(intersection)
            theCount = filtered.count(axis=1, asarray=True)
        else:
            theCount = self.count(axis=1, asarray=True)

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
        method: {'backfill', 'pad', None}
            Method to use for filling holes in new inde

        value: any kind (should be same type as array)
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
        colName: particular column name requested
        fromDate: datetime
        toDate: datetime
        nPeriods: int/float

        NOTE: Error thrown if all of fromDate, toDate, nPeriods specified.
        """
        if toDate:
            if toDate not in self.index:
                if toDate > self.index[0]:
                    toDate = self.index.asOfDate(toDate)
                else:
                    raise Exception('End date after last date in this index!')
        if fromDate:
            if fromDate not in self.index:
                if fromDate < self.index[-1]:
                    fromDate = self.index.asOfDate(fromDate)
                else:
                    raise Exception('Begin date after last date in this index!')
        if fromDate and toDate:
            if nPeriods:
                raise Exception('fromDate/toDate, toDate/nPeriods, ' + \
                                ' fromDate/nPeriods are mutually exclusive')
            beg_slice = self.index.indexMap[fromDate]
            end_slice = self.index.indexMap[toDate] + 1
        elif fromDate and nPeriods:
            beg_slice = self.index.indexMap[fromDate]
            end_slice = self.index.indexMap[fromDate] + nPeriods
        elif toDate and nPeriods:
            beg_slice = self.index.indexMap[toDate] - nPeriods + 1
            end_slice = self.index.indexMap[toDate] + 1
        else:
            raise Exception('Not enough arguments provided to getTS')

        # Fix indices in case they fall out of the boundaries
        beg_slice = max(0, beg_slice)
        end_slice = min(len(self.index), end_slice)
        dateRange = self.index[beg_slice:end_slice]

        if colName:
            return self[colName][beg_slice:end_slice]
        else:
            newColumns = {}
            for col, series in self.iteritems():
                newColumns[col] = series[beg_slice:end_slice]
            return DataFrame(data=newColumns, index=dateRange)

    def truncate(self, before = None, after = None):
        """
        Placeholder for documentation
        """
        import bisect

        before = to_datetime(before)
        after = to_datetime(after)

        if before is not None:
            binsearch = bisect.bisect_left(self.index, before)
            cur = self.index[binsearch]
            next = self.index[min(binsearch + 1, len(self.index) - 1)]
            leftDate = cur if cur >= before else next
        else:
            leftDate = self.index[0]

        if after is not None:
            if after < self.index[-1]:
                binsearch = bisect.bisect_right(self.index, after)
                cur = self.index[binsearch]
                prior = self.index[max(binsearch - 1, 0)]
                rightDate = cur if cur <= after else prior
            else:
                rightDate = self.index[-1]
        else:
            rightDate = self.index[-1]

        return self.getTS(fromDate=leftDate, toDate=rightDate)

    def getXS(self, key, subset=None, asOf=False):
        """
        Returns a row from the DataFrame as a Series object.

        Parameters
        ----------
        key: some index contained in the index
        subset: iterable (list, array, set, etc.), optional
            columns to be included
        asOf: boolean, optional
            Whether to use asOf values for TimeSeries objects
            Won't do anything for Series objects.

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
        index: string or object
            Column name to use to make new frame's index
        columns: string or object
            Column name to use to make new frame's columns
        values: string or object
            Column name to use for populating new frame's values
        """
        from pandas.core.panel import pivot, _slow_pivot

        return _slow_pivot(self[index], self[columns], self[values])

    def reindex(self, newIndex, fillMethod = None):
        """
        Reindex data inside, optionally filling according to some rule.

        Parameters
        ----------
        newIndex:   array-like
            preferably an Index object (to avoid duplicating data)
        fillMethod: {'backfill', 'pad', 'interpolate', None}
            Method to use for filling holes in reindexed DataFrame
        """
        if newIndex is self.index:
            return self.copy()

        if len(newIndex) == 0:
            return DataFrame(index=NULL_INDEX)

        if not isinstance(newIndex, Index):
            newIndex = Index(newIndex)

        if len(self.index) == 0:
            return DataFrame(index=newIndex)

        oldMap = self.index.indexMap
        newMap = newIndex.indexMap

        fillMethod = fillMethod.upper() if fillMethod else ''
        if fillMethod not in ['BACKFILL', 'PAD', '']:
            raise Exception("Don't recognize fillMethod: %s" % fillMethod)

        fillVec, mask = tseries.getFillVec(self.index, newIndex, oldMap,
                                           newMap, fillMethod)

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
            for type, dest in typeHierarchy:
                if issubclass(series.dtype.type, type):
                    new = series.take(fillVec).astype(dest)
                    new[-mask] = missingValue[dest]
                    newSeries[col] = new
                    break

        return DataFrame(newSeries, index=newIndex)

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

    def diff(self, periods = 1):
        temp = self.values
        temp = temp[periods:] - temp[:-periods]
        return self.fromMatrix(temp, self.cols(), self.index[periods:])

    def shift(self, periods, offset=None):
        """
        Shift the underlying series of the DataFrame and Series objects within
        by given number (positive or negative) of business/weekdays.

        Note, nan values inserted at beginning of series.
        """
        if periods == 0:
            return self

        if offset is None:
            if periods > 0:
                newIndex = self.index[periods:]
                newValues = dict([(col, np.asarray(series)[:-periods])
                                   for col, series in self.iteritems()])
            else:
                newIndex = self.index[:periods]
                newValues = dict([(col, np.asarray(series)[-periods:])
                                   for col, series in self.iteritems()])
        else:
            offset = periods * offset
            newIndex = Index([idx + offset for idx in self.index])
            newValues = dict([(col, np.asarray(series))
                               for col, series in self.iteritems()])
        return DataFrame(data = newValues, index= newIndex)

    def apply(self, func):
        """
        Applies func to columns (Series) of this DataFrame and returns either
        a DataFrame (if the function produces another series) or a Series
        indexed on the column names of the DataFrame if the function produces
        a value.

        Parameters
        ----------
        func: function
            Function to apply to each column

        Example
        -------

            >>> df.apply(numpy.sqrt) --> DataFrame
            >>> df.apply(numpy.sum) --> Series

        N.B.: Do NOT use functions that might toy with the index.
        """
        if not len(self.cols()):
            return self

        results = {}
        for col, series in self.iteritems():
            result = func(series)

            # WORKAROUND FOR NUMPY/SCIPY FUNCTIONS RETURNING UNSIZED NDARRAY
            if isinstance(result, np.ndarray) and result.ndim == 0:
                result = result.item()

            results[col] = result

        if hasattr(results.values()[0], '__iter__'):
            return DataFrame(data=results, index=self.index)
        else:
            keyArray = np.asarray(sorted(set(results.keys())), dtype=object)
            newIndex = Index(keyArray)

            arr = np.array([results[idx] for idx in newIndex])
            return Series(arr, index=newIndex)

    def tapply(self, func):
        """
        Apply func to the transposed DataFrame, results as per apply
        """
        return self.T.apply(func)

    def applymap(self, func):
        """
        Apply a function to a DataFrame that is intended to operate elementwise

        Please try to use apply if possible

        Parameters
        ----------
        func: function
            Python function to apply to each element
        """
        results = {}
        for col, series in self.iteritems():
            results[col] = map(func, series)
        return DataFrame(data=results, index=self.index)

    def tgroupby(self, keyfunc, applyfunc):
        """
        Call groupby on transposed frame
        """
        return self.T.groupby(keyfunc).aggregate(applyfunc).T

    def filterItems(self, items):
        """
        Restrict frame's columns to input set of items.

        Parameters
        ----------
        items: list-like
            List of columns to restrict to (must not all be present)

        Returns
        -------
        DataFrame with filtered columns
        """
        data = dict([(r, self[r]) for r in items if r in self])
        return DataFrame(data=data, index=self.index)

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

    def filterLike(self, arg):
        """
        Filter to columns partially matching the import argument.

        Keep columns where "arg in col == True"

        Parameter
        ---------
        arg: string

        Return
        ------
        DataFrame with matching columns
        """
        mycopy = self.copy()
        for col in mycopy._series.keys():
            series = mycopy._series.pop(col)
            if arg in col:
                mycopy._series[col] = series
        return mycopy

    def combineFirst(self, otherFrame):
        """
        Combine two DataFrame / DataMatrix objects and default to value
        in frame calling the method.

        Example: a.combineFirst(b)
            a's values prioritized, use values from b to fill holes

        Parameters
        ----------
        otherFrame: DataFrame / Matrix

        Returns
        -------
        DataFrame
        """
        if not otherFrame:
            return self

        if not self:
            return otherFrame

        if self.index is not otherFrame.index:
            unionIndex = self.index + otherFrame.index
            frame = self.reindex(unionIndex)
            otherFrame = otherFrame.reindex(unionIndex)
        else:
            unionIndex = self.index
            frame = self

        result = {}
        for col, series in frame.iteritems():
            otherSeries = otherFrame[col] if col in otherFrame else None
            if otherSeries is not None:
                result[col] = series.__class__(np.where(isnull(series),
                                                        otherSeries, series),
                                               index=unionIndex)
            else:
                result[col] = series

        for col, series in otherFrame.iteritems():
            if col not in self:
                result[col] = series

        return DataFrame(result, index = unionIndex)

    def combineAdd(self, otherFrame):
        """
        Add two DataFrame / DataMatrix objects and do not propagate NaN values,
        so if for a (column, time) one frame is missing a value, it will
        default to the other frame's value (which might be NaN as well)

        Parameters
        ----------
        otherFrame: DataFrame / Matrix

        Returns
        -------
        DataFrame
        """
        if not otherFrame:
            return self

        if not self:
            return otherFrame

        if self.index is not otherFrame.index:
            unionIndex = self.index + otherFrame.index
            frame = self.reindex(unionIndex)
            otherFrame = otherFrame.reindex(unionIndex)
        else:
            unionIndex = self.index
            frame = self

        unionCols = sorted(set(frame.cols() + otherFrame.cols()))

        result = {}
        for col in unionCols:
            if col in frame and col in otherFrame:
                series = frame[col].view(np.ndarray)
                otherSeries = otherFrame[col].view(np.ndarray)
                sok = np.isfinite(series)
                ook = np.isfinite(otherSeries)

                result[col] = np.where(sok & ook, series + otherSeries,
                                       np.where(sok, series, otherSeries))

            elif col in frame:
                result[col] = frame[col]
            elif col in otherFrame:
                result[col] = otherFrame[col]
            else:
                raise Exception('Phantom column, be very afraid')

        return DataFrame(result, index = unionIndex)

    def combineMult(self, otherFrame):
        return (self * otherFrame).combineFirst(self)

    def outerJoin(self, *frames):
        mergedSeries = self._series.copy()

        unionIndex = self.index
        for frame in frames:
            unionIndex  = unionIndex + frame.index

        for frame in frames:
            for col, series in frame.iteritems():
                if col in mergedSeries:
                    raise Exception('Overlapping columns!')
                mergedSeries[col] = series

        return DataFrame.fromDict(mergedSeries)

    def leftJoin(self, *frames):
        """
        Insert columns of input DataFrames / dicts into this one.

        Columns must not overlap. Returns a copy.

        Parameters
        ----------
        *frames: list-like
            List of frames (DataMatrix or DataFrame) as function arguments

        Returns
        -------
        DataFrame
        """
        mergedSeries = DataFrame(index=self.index)
        mergedSeries._series = self._series.copy()

        for frame in frames:
            for col, series in frame.iteritems():
                if col in mergedSeries:
                    raise Exception('Overlapping columns!')
                mergedSeries[col] = series

        return mergedSeries

    def merge(self, other, on=None):
        """
        Merge DataFrame or DataMatrix with this one on some many-to-one index

        Parameters
        ----------
        other: DataFrame
            Index should be similar to one of the columns in this one
        on: string
            Column name to use

        Example
        -------
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
            arr = series.view(ndarray).take(fillVec)
            arr[-mask] = NaN

            newSeries[col] = arr

        newSeries.update(self._series)

        return DataFrame(newSeries, index=self.index)

    def plot(self, kind='line', **kwds):
        """
        Plot the DataFrame's series with the index on the x-axis using
        matplotlib / pylab.

        Params
        ------
        kind: {'line', 'bar', 'hist'}
            Default: line for TimeSeries, hist for Series

        kwds: other plotting keyword arguments

        NOTE: This method doesn't make much sense for cross-sections,
        and will error.
        """
        try:
            plot
        except Exception, e:
            from pylab import plot

        for col in sorted(self.columns):
            s = self[col]
            plot(s.index, s, label=col)

    # ndarray-like stats methods
    def count(self, axis=0, asarray=False):
        """
        Return array or Series of # observations over requested axis.

        Parameters
        ----------
        axis: {0, 1}
            0 for row-wise, 1 for column-wise
        asarray: boolean, default False
            Choose to return as ndarray or have index attached

        Returns
        -------
        Series or TimeSeries
        """
        try:
            theCount = np.isfinite(self.values).sum(axis)
        except Exception, e:
            f = lambda s: notnull(s).sum()
            if axis == 0:
                theCount = self.apply(f)
            elif axis == 1:
                theCount = self.tapply(f)

        if asarray:
            return theCount
        else:
            if axis == 0:
                return Series(theCount, index=self.columns)
            elif axis == 1:
                return Series(theCount, index=self.index)
            else:
                raise Exception('Must have 0<= axis <= 1')

    def sum(self, axis=0, asarray=False):
        """
        Return array or Series of sums over requested axis.

        Parameters
        ----------
        axis: {0, 1}
            0 for row-wise, 1 for column-wise
        asarray: boolean, default False
            Choose to return as ndarray or have index attached

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
        except Exception, e:
            if axis == 0:
                theSum = self.apply(np.sum)
            else:
                theSum = self.tapply(np.sum)

        if asarray:
            return theSum
        if axis == 0:
            return Series(theSum, index=self.columns)
        elif axis == 1:
            return Series(theSum, index=self.index)
        else:
            raise Exception('Must have 0<= axis <= 1')

    def product(self, axis=0, asarray=False):
        """
        Return array or Series of products over requested axis.

        Parameters
        ----------
        axis: {0, 1}
            0 for row-wise, 1 for column-wise
        asarray: boolean, default False
            Choose to return as ndarray or have index attached

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
        except Exception, e:
            if axis == 0:
                theProd = self.apply(np.prod)
            else:
                theProd = self.tapply(np.prod)

        if asarray:
            return theProd
        if axis == 0:
            return Series(theProd, index=self.columns)
        elif axis == 1:
            return Series(theProd, index=self.index)
        else:
            raise Exception('Must have 0<= axis <= 1')

    def mean(self, axis=0):
        """
        Return array or Series of means over requested axis.

        Parameters
        ----------
        axis: {0, 1}
            0 for row-wise, 1 for column-wise

        Returns
        -------
        Series or TimeSeries
        """
        return self.sum(axis) / self.count(axis, asarray=True).astype(float)

    def median(self, axis=0):
        """
        Return array or Series of medians over requested axis.

        Parameters
        ----------
        axis: {0, 1}
            0 for row-wise, 1 for column-wise

        Returns
        -------
        Series or TimeSeries
        """
        if axis == 0:
            med = [np.median(remove_na(self[col])) for col in self.columns]
            return Series(med, index=self.columns)
        elif axis == 1:
            med = [np.median(remove_na(self.getXS(k))) for k in self.index]
            return Series(med, index=self.index)
        else:
            raise Exception('Must have 0<= axis <= 1')

    def min(self, axis=0):
        """
        Return array or Series of minimums over requested axis.

        Parameters
        ----------
        axis: {0, 1}
            0 for row-wise, 1 for column-wise

        Returns
        -------
        Series or TimeSeries
        """
        if axis == 0:
            med = [np.min(remove_na(self[col])) for col in self.columns]
            return Series(med, index=self.columns)
        elif axis == 1:
            med = [np.min(remove_na(self.getXS(k))) for k in self.index]
            return Series(med, index=self.index)
        else:
            raise Exception('Must have 0<= axis <= 1')

    def max(self, axis=0):
        """
        Return array or Series of maximums over requested axis.

        Parameters
        ----------
        axis: {0, 1}
            0 for row-wise, 1 for column-wise

        Returns
        -------
        Series or TimeSeries
        """
        if axis == 0:
            med = [np.max(remove_na(self[col])) for col in self.columns]
            return Series(med, index=self.columns)
        elif axis == 1:
            med = [np.max(remove_na(self.getXS(k))) for k in self.index]
            return Series(med, index=self.index)
        else:
            raise Exception('Must have 0<= axis <= 1')

    def mad(self, axis=0, asarray=False):
        """
        Return array or Series of mean absolute deviation over
        requested axis.

        Parameters
        ----------
        axis: {0, 1}
            0 for row-wise, 1 for column-wise
        asarray: boolean, default False
            Choose to return as ndarray or have index attached

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

        if asarray:
            return np.sum(np.abs(y), axis)

        if axis == 0:
            return Series(np.sum(np.abs(y), axis), self.cols())
        else:
            return Series(np.sum(np.abs(y), axis), self.index)

    def var(self, axis=0, asarray=False):
        """
        Return array or Series of unbiased variance over requested axis.

        Parameters
        ----------
        axis: {0, 1}
            0 for row-wise, 1 for column-wise
        asarray: boolean, default False
            Choose to return as ndarray or have index attached

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
        if asarray:
            return theVar
        if axis == 0:
            return Series(theVar, index=self.columns)
        elif axis == 1:
            return Series(theVar, index=self.index)
        else:
            raise Exception('Must have 0<= axis <= 1')

    def std(self, axis=0, asarray=False):
        """
        Return array or Series of unbiased std deviation over requested axis.

        Parameters
        ----------
        axis: {0, 1}
            0 for row-wise, 1 for column-wise
        asarray: boolean, default False
            Choose to return as ndarray or have index attached

        Returns
        -------
        Series or TimeSeries
        """
        return np.sqrt(self.var(axis=axis, asarray=asarray))

    def skew(self, axis=0, asarray=False):
        """
        Return array or Series of unbiased skewness over requested axis.

        Parameters
        ----------
        axis: {0, 1}
            0 for row-wise, 1 for column-wise
        asarray: boolean, default False
            Choose to return as ndarray or have index attached

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

        if asarray:
            return theSkew
        if axis == 0:
            return Series(theSkew, index=self.columns)
        elif axis == 1:
            return Series(theSkew, index=self.index)
        else:
            raise Exception('Must have 0<= axis <= 1')

    # TODO
    def kurtosis(self, axis=0):
        """
        Return array or Series of unbiased kurtosis over requested axis.

        Parameters
        ----------
        axis: {0, 1}
            0 for row-wise, 1 for column-wise
        asarray: boolean, default False
            Choose to return as ndarray or have index attached

        Returns
        -------
        Series or TimeSeries
        """
        raise Exception('Not implemented yet!')

    def _withColumns(self, newCols):
        """
        Utility method, force values matrix to have particular columns
        Can make this as cute as we like
        """
        if len(newCols) == 0:
            return DataFrame(index=self.index)

        newFrame = self.filterItems(newCols)

        for col in newCols:
            if col not in newFrame:
                newFrame[col] = NaN

        return newFrame
def _pfixed(s, space, nanRep=None):
    if isinstance(s, float):
        fstring = '%-' + str(space-4) + 'g'
        if nanRep is not None and isnull(s):
            return nanRep.ljust(space)
        return (fstring % s).ljust(space)
    else:
        return str(s)[:space-4].ljust(space)

