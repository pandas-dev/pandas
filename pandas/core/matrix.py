# pylint: disable-msg=E1101,E1103
# pylint: disable-msg=W0212,W0703,W0231,W0622

from cStringIO import StringIO
import sys

from numpy import NaN
import numpy as np

from pandas.core.common import _pfixed, _pickle_array, _unpickle_array
from pandas.core.frame import DataFrame, _try_sort, _extract_index
from pandas.core.index import Index, NULL_INDEX
from pandas.core.series import Series
from pandas.lib.tseries import isnull
import pandas.core.datetools as datetools
import pandas.lib.tseries as tseries

#-------------------------------------------------------------------------------
# DataMatrix class

class DataMatrix(DataFrame):
    """
    Matrix version of DataFrame, optimized for cross-section operations,
    numerical computation, and other operations that do not require the
    frame to change size.

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
    Transposing is much faster in this regime, as is calling getXS, so please
    take note of this.
    """
    objects = None
    def __init__(self, data=None, index=None, columns=None, dtype=None,
                 objects=None):

        if isinstance(data, dict) and len(data) > 0:
            (index, columns,
             values, objects) = self._initDict(data, index, columns, objects,
                                               dtype)
        elif isinstance(data, (np.ndarray, list)):
            (index, columns, values) = self._initMatrix(data, index,
                                                        columns, dtype)

            if objects is not None:
                if isinstance(objects, DataMatrix):
                    if not objects.index.equals(index):
                        objects = objects.reindex(index)
                else:
                    objects = DataMatrix(objects, index=index)

        elif data is None or len(data) == 0:
            # this is a touch convoluted...
            if objects is not None:
                if isinstance(objects, DataMatrix):
                    if index is not None and objects.index is not index:
                        objects = objects.reindex(index)
                else:
                    objects = DataMatrix(objects, index=index)

                index = objects.index

            if index is None:
                N = 0
                index = NULL_INDEX
            else:
                N = len(index)

            if columns is None:
                K = 0
                columns = NULL_INDEX
            else:
                K = len(columns)

            values = np.empty((N, K), dtype=dtype)
            values[:] = NaN
        else:
            raise Exception('DataMatrix constructor not properly called!')

        self.values = values
        self.index = index
        self.columns = columns
        self.objects = objects

    def _initDict(self, data, index, columns, objects, dtype):
        """
        Segregate Series based on type and coerce into matrices.

        Needs to handle a lot of exceptional cases.

        Somehow this got outrageously complicated
        """
        # pre-filter out columns if we passed it
        if columns is not None:
            colset = set(columns)
            data = dict((k, v) for k, v in data.iteritems() if k in colset)

        index = _extract_index(data, index)

        objectDict = {}
        if objects is not None and isinstance(objects, dict):
            objectDict.update(objects)

        valueDict = {}
        for k, v in data.iteritems():
            if isinstance(v, Series):
                # Forces alignment, copies data
                v = v.reindex(index)
            else:
                if isinstance(v, dict):
                    v = [v.get(i, NaN) for i in index]
                else:
                    assert(len(v) == len(index))

                try:
                    v = Series(v, dtype=dtype, index=index)
                except Exception:
                    v = Series(v, index=index)

                # copy data
                v = v.copy()

            if issubclass(v.dtype.type, (np.bool_, float, int)):
                valueDict[k] = v
            else:
                objectDict[k] = v

        if columns is None:
            columns = Index(_try_sort(valueDict))
            objectColumns = Index(_try_sort(objectDict))
        else:
            objectColumns = Index([c for c in columns if c in objectDict])
            columns = Index([c for c in columns if c not in objectDict])

        if len(valueDict) == 0:
            dtype = np.object_
            valueDict = objectDict
            columns = objectColumns
        else:
            dtype = np.float_
            if len(objectDict) > 0:
                new_objects = DataMatrix(objectDict,
                                         dtype=np.object_,
                                         index=index,
                                         columns=objectColumns)
                if isinstance(objects, DataMatrix):
                    objects = objects.join(new_objects, how='left')
                else:
                    objects = new_objects

        values = np.empty((len(index), len(columns)), dtype=dtype)

        for i, col in enumerate(columns):
            if col in valueDict:
                values[:, i] = valueDict[col]
            else:
                values[:, i] = np.NaN

        return index, columns, values, objects

    def _initMatrix(self, values, index, columns, dtype):
        if not isinstance(values, np.ndarray):
            arr = np.array(values)
            if issubclass(arr.dtype.type, basestring):
                arr = np.array(values, dtype=object, copy=True)

            values = arr

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

        if index is None:
            if values.shape[0] == 0:
                index = NULL_INDEX
            else:
                raise Exception('Must pass index!')

        if columns is None:
            if values.shape[1] == 0:
                columns = NULL_INDEX
            else:
                raise Exception('Must pass columns!')

        return index, columns, values

    @property
    def _constructor(self):
        return DataMatrix

    # Because of DataFrame property
    values = None

    def __array__(self):
        return self.values

    def __array_wrap__(self, result):
        return DataMatrix(result, index=self.index, columns=self.columns)

#-------------------------------------------------------------------------------
# DataMatrix-specific implementation of private API

    def _join_on(self, other, on):
        if len(other.index) == 0:
            return self

        if on not in self:
            raise Exception('%s column not contained in this frame!' % on)

        fillVec, mask = tseries.getMergeVec(self[on],
                                            other.index.indexMap)

        tmpMatrix = other.values.take(fillVec, axis=0)
        tmpMatrix[-mask] = NaN

        seriesDict = dict((col, tmpMatrix[:, j])
                           for j, col in enumerate(other.columns))

        if getattr(other, 'objects'):
            objects = other.objects

            tmpMat = objects.values.take(fillVec, axis=0)
            tmpMat[-mask] = NaN
            objDict = dict((col, tmpMat[:, j])
                           for j, col in enumerate(objects.columns))

            seriesDict.update(objDict)

        filledFrame = DataFrame(data=seriesDict, index=self.index)

        return self.join(filledFrame, how='left')

    def _reindex_index(self, index, method):
        if index is self.index:
            return self.copy()

        if not isinstance(index, Index):
            index = Index(index)

        if len(self.index) == 0:
            return DataMatrix(index=index, columns=self.columns)

        fillVec, mask = tseries.getFillVec(self.index, index,
                                           self.index.indexMap,
                                           index.indexMap, method)

        mat = self.values.take(fillVec, axis=0)

        notmask = -mask
        if len(index) > 0:
            if notmask.any():
                if issubclass(mat.dtype.type, np.int_):
                    mat = mat.astype(float)
                elif issubclass(mat.dtype.type, np.bool_):
                    mat = mat.astype(object)

                mat[-mask] = NaN

        if self.objects is not None and len(self.objects.columns) > 0:
            newObjects = self.objects.reindex(index)
        else:
            newObjects = None

        return DataMatrix(mat, index=index, columns=self.columns,
                          objects=newObjects)

    def _reindex_columns(self, columns):
        if len(columns) == 0:
            return DataMatrix(index=self.index)

        if not isinstance(columns, Index):
            columns = Index(columns)

        if self.objects is not None:
            object_columns = columns.intersection(self.objects.columns)
            columns = columns - object_columns

            objects = self.objects._reindex_columns(object_columns)
        else:
            objects = None

        if len(columns) > 0 and len(self.columns) == 0:
            return DataMatrix(index=self.index, columns=columns,
                              objects=objects)

        indexer, mask = tseries.getFillVec(self.columns, columns,
                                           self.columns.indexMap,
                                           columns.indexMap, None)

        mat = self.values.take(indexer, axis=1)

        notmask = -mask
        if len(mask) > 0:
            if notmask.any():
                if issubclass(mat.dtype.type, np.int_):
                    mat = mat.astype(float)
                elif issubclass(mat.dtype.type, np.bool_):
                    mat = mat.astype(object)

                mat[:, -mask] = NaN

        return DataMatrix(mat, index=self.index, columns=columns,
                          objects=objects)

    def _rename_columns_inplace(self, mapper):
        self.columns = [mapper(x) for x in self.columns]

        if self.objects is not None:
            self.objects._rename_columns_inplace(mapper)

    def _combineFrame(self, other, func):
        """
        Methodology, briefly
        - Really concerned here about speed, space

        - Get new index
        - Reindex to new index
        - Determine newColumns and commonColumns
        - Add common columns over all (new) indices
        - Fill to new set of columns

        Could probably deal with some Cython action in here at some point
        """
        if self.index.equals(other.index):
            newIndex = self.index
            myReindex = self
            hisReindex = other
        else:
            newIndex = self.index.union(other.index)
            myReindex = self.reindex(newIndex)
            hisReindex = other.reindex(newIndex)

        if not self and not other:
            return DataMatrix(index=newIndex)
        elif not self:
            return other * NaN
        elif not other:
            return self * NaN

        myValues = myReindex.values
        safe = self.columns.equals(other.columns)

        if safe:
            newCols = self.columns
            commonCols = self.columns
        else:
            newCols = self.columns.union(other.columns)
            commonCols = self.columns.intersection(other.columns)

        hisValues = hisReindex.values
        hisCols = hisReindex.columns

        if safe:
            resultMatrix = func(myValues, hisValues)
        else:
            T, N = len(newIndex), len(newCols)
            resultMatrix = np.empty((T, N), dtype=self.values.dtype)
            resultMatrix.fill(NaN)

            myIndexer = [self.columns.indexMap[idx] for idx in commonCols]
            hisIndexer =  [hisCols.indexMap[idx] for idx in commonCols]
            resultIndexer = [newCols.indexMap[idx] for idx in commonCols]

            resultMatrix[:, resultIndexer] = func(myValues[:, myIndexer],
                                                  hisValues[:, hisIndexer])

        # TODO: deal with objects
        return DataMatrix(resultMatrix, index=newIndex, columns=newCols)

    def _combineSeries(self, other, func):
        newIndex = self.index
        newCols = self.columns

        if len(self) == 0:
            # Ambiguous case
            return DataMatrix(index=self.index, columns=self.columns,
                              objects=self.objects)

        if self.index._allDates and other.index._allDates:
            # Operate row-wise
            if self.index.equals(other.index):
                newIndex = self.index
            else:
                newIndex = self.index + other.index

            other = other.reindex(newIndex).view(np.ndarray)
            myReindex = self.reindex(newIndex)
            resultMatrix = func(myReindex.values.T, other).T
        else:
            if len(other) == 0:
                return self * NaN

            newCols = self.columns.union(other.index)

            # Operate column-wise
            this = self.reindex(columns=newCols)
            other = other.reindex(newCols).values()

            resultMatrix = func(this.values, other)

        # TODO: deal with objects
        return DataMatrix(resultMatrix, index=newIndex, columns=newCols)

    def _combineFunc(self, other, func):
        """
        Combine DataMatrix objects with other Series- or DataFrame-like objects

        This is the core method used for the overloaded arithmetic methods

        Result hierarchy
        ----------------
        DataMatrix + DataFrame --> DataMatrix
        DataMatrix + DataMatrix --> DataMatrix
        DataMatrix + Series --> DataMatrix
        DataMatrix + constant --> DataMatrix

        The reason for 'upcasting' the result is that if addition succeed,
        we can assume that the input DataFrame was homogeneous.
        """
        newIndex = self.index
        if isinstance(other, DataFrame):
            return self._combineFrame(other, func)

        elif isinstance(other, Series):
            return self._combineSeries(other, func)

        else:
            if not self:
                return self

            # Constant of some kind
            newCols = self.columns
            resultMatrix = func(self.values, other)

        # TODO: deal with objects
        return DataMatrix(resultMatrix, index=newIndex, columns=newCols)

#-------------------------------------------------------------------------------
# Properties for index and columns

    _columns = None
    def _get_columns(self):
        return self._columns

    def _set_columns(self, cols):
        if len(cols) != self.values.shape[1]:
            raise Exception('Columns length %d did not match values %d!' %
                            (len(cols), self.values.shape[1]))

        if not isinstance(cols, Index):
            cols = Index(cols)

        self._columns = cols

    columns = property(fget=_get_columns, fset=_set_columns)

    def _set_index(self, index):
        if len(index) > 0:
            if len(index) != self.values.shape[0]:
                raise Exception('Index length %d did not match values %d!' %
                                (len(index), self.values.shape[0]))

        if not isinstance(index, Index):
            index = Index(index)

        self._index = index

        if self.objects is not None:
            self.objects._index = index

    def _get_index(self):
        return self._index

    index = property(fget=_get_index, fset=_set_index)

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

    def __nonzero__(self):
        N, K = self.values.shape
        if N == 0 or K == 0:
            if self.objects is None:
                return False
            else:
                return self.objects.__nonzero__()
        else:
            return True

    def __neg__(self):
        mycopy = self.copy()
        mycopy.values = -mycopy.values
        return mycopy

    def __repr__(self):
        """Return a string representation for a particular DataMatrix"""
        buffer = StringIO()

        if len(self.cols()) == 0:
            buffer.write('Empty DataMatrix\nIndex: %s' % repr(self.index))
        elif 0 < len(self.index) < 500 and self.values.shape[1] < 10:
            self.toString(buffer=buffer)
        else:
            print >> buffer, str(self.__class__)
            self.info(buffer=buffer)

        return buffer.getvalue()

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

        Note
        ----
        This is a magic method. Do NOT call explicity.
        """
        if isinstance(item, slice):
            indexRange = self.index[item]
            return self.reindex(indexRange)

        elif isinstance(item, np.ndarray):
            if len(item) != len(self.index):
                raise Exception('Item wrong length %d instead of %d!' %
                                (len(item), len(self.index)))
            newIndex = self.index[item]
            return self.reindex(newIndex)
        else:
            if self.objects is not None and item in self.objects:
                return self.objects[item]
            else:
                return self._getSeries(item)

    _dataTypes = [np.float_, np.bool_, np.int_]
    def __setitem__(self, key, value):
        """
        Add series to DataMatrix in specified column.

        If series is a numpy-array (not a Series/TimeSeries), it must be the
        same length as the DataMatrix's index or an error will be thrown.

        Series/TimeSeries will be conformed to the DataMatrix's index to
        ensure homogeneity.
        """
        import bisect

        isObject = False
        if hasattr(value, '__iter__'):
            if isinstance(value, Series):
                value = np.asarray(value.reindex(self.index))

            else:
                assert(len(value) == len(self.index))

                if not isinstance(value, np.ndarray):
                    value = np.array(value)
                    if value.dtype.type == np.str_:
                        value = np.array(value, dtype=object)
        else:
            value = np.repeat(value, len(self.index))

        if value.dtype not in self._dataTypes:
            isObject = True

        if self.values.dtype == np.object_:
            if key in self.columns:
                loc = self.columns.indexMap[key]
                self.values[:, loc] = value
            elif len(self.columns) == 0:
                self.values = value.reshape((len(value), 1)).copy()
                self.columns = Index([key])
            else:
                try:
                    loc = bisect.bisect_right(self.columns, key)
                except TypeError:
                    loc = len(self.columns)

                if loc == self.values.shape[1]:
                    newValues = np.c_[self.values, value]
                    newColumns = Index(np.concatenate((self.columns, [key])))
                elif loc == 0:
                    newValues = np.c_[value, self.values]
                    newColumns = Index(np.concatenate(([key], self.columns)))
                else:
                    newValues = np.c_[self.values[:, :loc], value,
                                      self.values[:, loc:]]
                    toConcat = (self.columns[:loc], [key], self.columns[loc:])
                    newColumns = Index(np.concatenate(toConcat))
                self.values = newValues
                self.columns = newColumns
        else:
            if key in self.columns:
                loc = self.columns.indexMap[key]
                self.values[:, loc] = value
            elif isObject:
                if self.objects is None:
                    self.objects = DataMatrix({key : value},
                                              index=self.index)
                else:
                    self.objects[key] = value
            elif len(self.columns) == 0:
                self.values = value.reshape((len(value), 1)).astype(np.float)
                self.columns = Index([key])
            else:
                try:
                    loc = bisect.bisect_right(self.columns, key)
                except TypeError:
                    loc = len(self.columns)

                if loc == self.values.shape[1]:
                    newValues = np.c_[self.values, value]
                    newColumns = Index(np.concatenate((self.columns, [key])))
                elif loc == 0:
                    newValues = np.c_[value, self.values]
                    newColumns = Index(np.concatenate(([key], self.columns)))
                else:
                    newValues = np.c_[self.values[:, :loc], value,
                                      self.values[:, loc:]]
                    toConcat = (self.columns[:loc], [key], self.columns[loc:])
                    newColumns = Index(np.concatenate(toConcat))
                self.values = newValues
                self.columns = newColumns

    def __delitem__(self, key):
        """
        Delete column from DataMatrix
        """
        if key in self.columns:
            loc = self.columns.indexMap[key]
            if loc == self.values.shape[1] - 1:
                newValues = self.values[:, :loc]
                newColumns = self.columns[:loc]
            else:
                newValues = np.c_[self.values[:, :loc], self.values[:, loc+1:]]
                newColumns = Index(np.concatenate((self.columns[:loc],
                                                   self.columns[loc+1:])))
            self.values = newValues
            self.columns = newColumns
        else:
            if self.objects is not None and key in self.objects:
                del self.objects[key]
            else:
                raise KeyError('%s' % key)

    def __iter__(self):
        """Iterate over columns of the frame."""
        return iter(self.columns)

    def __contains__(self, key):
        """True if DataMatrix has this column"""
        hasCol = key in self.columns
        if hasCol:
            return True
        else:
            if self.objects is not None and key in self.objects:
                return True
            return False

    def iteritems(self):
        return self._series.iteritems()

#-------------------------------------------------------------------------------
# Helper methods

    # For DataFrame compatibility
    def _getSeries(self, item=None, loc=None):
        if loc is None:
            try:
                loc = self.columns.indexMap[item]
            except KeyError:
                raise Exception('%s not here!' % item)
        return Series(self.values[:, loc], index=self.index)

    def _getSeriesDict(self):
        series = {}
        for i, col in enumerate(self.columns):
            series[col] = self._getSeries(loc=i)
        if self.objects is not None:
            for i, col in enumerate(self.objects.columns):
                series[col] = self.objects._getSeries(loc=i)

        return series
    _series = property(_getSeriesDict)

#-------------------------------------------------------------------------------
# Outputting

    def toCSV(self, path, nanRep='', writeMode='wb', index=True,
              header=True, cols=None, verbose=False):
        """
        Write the DataMatrix to a CSV file

        Parameters
        ----------
        path : string
            Output file path
        nanRep : string, default=''
            Appearance of NaN values in output
        index : boolean, default=True
            Prints index if True
        header : boolean, default=True
            Prints header if True
        cols : list of strings
            Prints the values in order specified by cols.
            By default, prints all columns in lexicographical order.
        """
        f = open(path, writeMode)

        if cols is None:
            cols = self.cols()
        series = self._series

        if header:
            if index:
                f.write(',')
            f.write(','.join([str(c) for c in cols]))
            f.write('\n')

        for idx in self.index:
            if index:
                f.write(str(idx) + ',')

            for col in cols:
                val = series[col][idx]
                if isnull(val):
                    val = nanRep
                else:
                    val = str(val)
                f.write(val + ',')
            f.write('\n')

        f.close()

        if verbose: # pragma: no cover
            print 'CSV file written successfully: %s' % path

    def toString(self, buffer=sys.stdout, columns=None, colSpace=15,
                 nanRep='NaN', formatters=None, float_format=None):
        """
        Output a string version of this DataMatrix
        """
        formatters = formatters or {}

        if columns is None:
            columns = self.columns
            values = self.values
            if self.objects:
                columns = list(columns) + list(self.objects.columns)
                values = np.column_stack((values.astype(object),
                                          self.objects.values))
        else:
            columns = [c for c in columns if c in self]
            values = self.asMatrix(columns)

        ident = lambda x: x

        idxSpace = max([len(str(idx)) for idx in self.index]) + 4

        if len(self.cols()) == 0:
            buffer.write('DataMatrix is empty!\n')
            buffer.write(repr(self.index))
        else:
            buffer.write(_pfixed('', idxSpace))
            for h in columns:
                buffer.write(_pfixed(h, colSpace))
            buffer.write('\n')

            for i, idx in enumerate(self.index):
                buffer.write(_pfixed(idx, idxSpace))
                for j, col in enumerate(columns):
                    formatter = formatters.get(col, ident)
                    buffer.write(_pfixed(formatter(values[i, j]), colSpace,
                                         float_format=float_format,
                                         nanRep=nanRep))
                buffer.write('\n')

    def info(self, buffer=sys.stdout):
        """
        Concise summary of a DataMatrix, used in __repr__ when very large.
        """
        print >> buffer, 'Index: %s entries' % len(self.index),
        if len(self.index) > 0:
            print >> buffer, ', %s to %s' % (self.index[0], self.index[-1])
        else:
            print >> buffer, ''

        if len(self.columns) == 0:
            print >> buffer, 'DataMatrix is empty!'
            print >> buffer, repr(self.index)
            return

        print >> buffer, 'Data columns:'
        space = max([len(str(k)) for k in self.cols()]) + 4

        counts = self.count()

        columns = []
        for j, col in enumerate(self.columns):
            columns.append((col, '%s%d  non-null values' %
                           (_pfixed(col, space), counts[j])))

        if self.objects is not None and len(self.objects.columns) > 0:
            n = len(self.objects.index)
            for col in self.objects:
                line = '%s%d  non-null values' % (_pfixed(col, space), n)
                columns.append((col, line))

        try:
            columns = [c[1] for c in sorted(columns)]
        except TypeError:
            columns = sorted([c[1] for c in columns])

        dtypeLine = ''

        nf = len(self.columns)
        df = self.values.dtype
        if self.objects is not None:
            no = len(self.objects.columns)
            do = self.objects.values.dtype
            dtypeLine = '\ndtypes: %s(%d), %s(%d)' % (df, nf, do, no)
        else:
            dtypeLine = '\ndtype: %s(%d)' % (df, nf)

        buffer.write('\n'.join(columns) + dtypeLine)


#-------------------------------------------------------------------------------
# Public methods

    def apply(self, func, axis=0):
        """
        Applies func to columns (Series) of this DataMatrix and returns either
        a DataMatrix (if the function produces another series) or a Series
        indexed on the column names of the DataFrame if the function produces
        a value.

        Parameters
        ----------
        func : function
            Function to apply to each column

        Examples
        --------

            >>> df.apply(numpy.sqrt) --> DataMatrix
            >>> df.apply(numpy.sum) --> Series

        N.B.: Do NOT use functions that might toy with the index.
        """
        if not len(self.cols()):
            return self

        if isinstance(func, np.ufunc):
            results = func(self.values)
            return DataMatrix(data=results, index=self.index,
                              columns=self.columns, objects=self.objects)
        else:
            return DataFrame.apply(self, func, axis=axis)

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

        if (isinstance(other, DataMatrix) and
            self.columns.equals(other.columns)):

            idx = Index(np.concatenate([self.index, other.index]))
            mat = np.vstack((self.values, other.values))

            if other.objects is None:
                objects = self.objects
            elif self.objects is None:
                objects = other.objects
            else:
                objects = self.objects.append(other.objects)

            if objects:
                objects = objects.reindex(idx)

            dm = DataMatrix(mat, idx, self.columns, objects=objects)
            return dm
        else:
            return super(DataMatrix, self).append(other)

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
        if columns is None:
            values = self.values.copy()

            if self.objects:
                values = np.column_stack((values, self.objects.values))

            return values
        else:
            if not isinstance(columns, Index):
                columns = Index(columns)

            values = self.values
            order = self.columns

            if self.objects:
                idxMap = self.objects.columns.indexMap
                indexer = [idxMap[col] for col in columns if col in idxMap]

                obj_values = self.objects.values.take(indexer, axis=1)

                values = np.column_stack((values, obj_values))
                order = Index(np.concatenate((order, self.objects.columns)))

                # now put in the right order

            values = _reorder_columns(values, order, columns)

            return values

    def cols(self):
        """Return sorted list of frame's columns"""
        if self.objects is not None and len(self.objects.columns) > 0:
            return list(self.columns.union(self.objects.columns))
        else:
            return list(self.columns)

    def copy(self):
        """
        Make a copy of this DataMatrix
        """
        if self.objects:
            objects = self.objects.copy()
        else:
            objects = None

        return DataMatrix(self.values.copy(), index=self.index,
                          columns=self.columns, objects=objects)

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

    def fill(self, value=None, method='pad'):
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
                result[col] = s.fill(method=method, value=value)

            return DataMatrix(result, index=self.index, objects=self.objects)
        else:
            if (isinstance(value, (int, float))
                and self.values.dtype == np.float64):
                # Float type values
                if len(self.columns) == 0:
                    return self

                vals = self.values.copy()
                vals.flat[isnull(vals.ravel())] = value

                objects = None

                if self.objects is not None:
                    objects = self.objects.copy()

                return DataMatrix(vals, index=self.index, columns=self.columns,
                                  objects=objects)

            else:
                # Object type values
                if len(self.columns) == 0:
                    return self

                # XXX

                myCopy = self.copy()
                vals = myCopy.values
                vals = self.values.copy()
                vals.flat[isnull(vals.ravel())] = value

                return myCopy

    def getXS(self, key):
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

        loc = self.index.indexMap[key]
        theSlice = self.values[loc, :].copy()
        xsIndex = self.columns

        result = Series(theSlice, index=xsIndex)

        if self.objects is not None and len(self.objects.columns) > 0:
            result = result.append(self.objects.getXS(key))

        return result

    @property
    def T(self):
        """
        Returns a DataMatrix with the rows/columns switched.
        """
        if self.objects is not None:
            objectsT = self.objects.values.T
            valuesT = self.values.T
            newValues = np.concatenate((valuesT, objectsT), axis=0)
            newIndex = Index(np.concatenate((self.columns,
                                             self.objects.columns)))

            return DataMatrix(newValues, index=newIndex, columns=self.index)
        else:
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

        N = len(self)

        if offset is None:
            newIndex = self.index
            indexer = np.zeros(N, dtype=int)
            if periods > 0:
                indexer[periods:] = np.arange(N - periods)
                newValues = self.values.take(indexer, axis=0)
                newValues[:periods] = NaN
            else:
                indexer[:periods] = np.arange(-periods, N)
                newValues = self.values.take(indexer, axis=0)
                newValues[periods:] = NaN
        else:
            offset = periods * offset
            newIndex = Index([x + offset for x in self.index])
            newValues = self.values.copy()

        if self.objects is not None:
            shifted_objects = self.objects.shift(periods, offset=offset,
                                                 timeRule=timeRule)

            shifted_objects.index = newIndex
        else:
            shifted_objects = None

        return DataMatrix(data=newValues, index=newIndex, columns=self.columns,
                          objects=shifted_objects)

    def cap(self, threshold):
        """
        Trim values at threshold

        Returns
        -------
        DataMatrix
        """
        return DataMatrix(np.where(self.values > threshold,
                                   threshold, self.values),
                          index=self.index, columns=self.columns,
                          objects=self.objects)

    def floor(self, threshold):
        """
        Trim values below threshold

        Returns
        -------
        DataMatrix
        """
        return DataMatrix(np.where(self.values < threshold,
                                   threshold, self.values),
                          index=self.index, columns=self.columns,
                          objects=self.objects)

def _reorder_columns(mat, current, desired):
    fillVec, mask = tseries.getFillVec(current, desired, current.indexMap,
                                       desired.indexMap, None)

    fillVec = fillVec[mask]

    return mat.take(fillVec, axis=1)
