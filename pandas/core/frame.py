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

from StringIO import StringIO
import operator
import warnings

from numpy import nan
import numpy as np

from pandas.core.common import (isnull, notnull, PandasError,
                                _try_sort, _pfixed, _default_index,
                                _infer_dtype, _stringify, _maybe_upcast)
from pandas.core.daterange import DateRange
from pandas.core.generic import AxisProperty, NDFrame
from pandas.core.index import Index, MultiIndex, NULL_INDEX, _ensure_index
from pandas.core.indexing import _NDFrameIndexer, _maybe_droplevels
from pandas.core.internals import BlockManager, make_block, form_blocks
from pandas.core.series import Series, _is_bool_indexer
from pandas.util.decorators import deprecate
import pandas.core.common as common
import pandas.core.datetools as datetools
import pandas._tseries as _tseries

#----------------------------------------------------------------------
# Factory helper methods

_arith_doc = """
Binary operator %s with support to substitute a fill_value for missing data in
one of the inputs

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
            return self._combine_series(other, func, fill_value, axis)
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

#----------------------------------------------------------------------
# DataFrame class


class DataFrame(NDFrame):
    _auto_consolidate = True
    _verbose_info = True
    _het_axis = 1

    _AXIS_NUMBERS = {
        'index' : 0,
        'columns' : 1
    }

    _AXIS_NAMES = dict((v, k) for k, v in _AXIS_NUMBERS.iteritems())

    def __init__(self, data=None, index=None, columns=None, dtype=None,
                 copy=False):
        """Two-dimensional size-mutable, potentially heterogeneous tabular data
        structure with labeled axes (rows and columns). Arithmetic operations
        align on both row and column labels. Can be thought of as a dict-like
        container for Series objects. The primary pandas data structure

        Parameters
        ----------
        data : numpy ndarray (structured or homogeneous), dict, or DataFrame
            Dict can contain Series, arrays, constants, or list-like objects
        index : Index or array-like
            Index to use for resulting frame. Will default to np.arange(n) if no
            indexing information part of input data and no index provided
        columns : Index or array-like
            Will default to np.arange(n) if not column labels provided
        dtype : dtype, default None
            Data type to force, otherwise infer
        copy : boolean, default False
            Copy data from inputs. Only affects DataFrame / 2d ndarray input

        Examples
        --------
        >>> d = {'col1' : ts1, 'col2' : ts2}
        >>> df = DataFrame(data=d, index=index)
        >>> df2 = DataFrame(np.random.randn(10, 5))
        >>> df3 = DataFrame(np.random.randn(10, 5),
                            columns=['a', 'b', 'c', 'd', 'e'])
        """

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
                mgr = mgr.astype(dtype)
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
        self._series_cache = {}

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

    def _wrap_array(self, arr, axes, copy=False):
        index, columns = axes
        return self._constructor(arr, index=index, columns=columns, copy=copy)

    @property
    def axes(self):
        return [self.index, self.columns]

    @property
    def _constructor(self):
        return DataFrame

    # Fancy indexing
    _ix = None

    @property
    def ix(self):
        if self._ix is None:
            self._ix = _NDFrameIndexer(self)

        return self._ix

    @property
    def shape(self):
        return (len(self.index), len(self.columns))

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
        if len(self.index) < 500 and len(self.columns) <= 10:
            self.to_string(buf=buf)
        else:
            self.info(buf=buf, verbose=self._verbose_info)

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
        """Returns length of index"""
        return len(self.index)

    def __contains__(self, key):
        """True if DataFrame has this column"""
        return key in self.columns

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
    __rpow__ = _arith_method(lambda x, y: y ** x, '__rpow__',
                             default_axis=None)

    def __neg__(self):
        return self * -1

    # Comparison methods
    __eq__ = comp_method(operator.eq, '__eq__')
    __ne__ = comp_method(operator.ne, '__ne__')
    __lt__ = comp_method(operator.lt, '__lt__')
    __gt__ = comp_method(operator.gt, '__gt__')
    __le__ = comp_method(operator.le, '__le__')
    __ge__ = comp_method(operator.ge, '__ge__')

    #----------------------------------------------------------------------
    # IO methods (to / from other formats)

    def to_dict(self):
        """
        Convert DataFrame to nested dictionary

        Returns
        -------
        result : dict like {column -> {index -> value}}
        """
        return dict((k, v.to_dict()) for k, v in self.iteritems())

    @classmethod
    def from_records(cls, data, index=None, indexField=None,
                     exclude=None):
        """
        Convert structured or record ndarray to DataFrame

        Parameters
        ----------
        data : NumPy structured array
        index : string, list of fields, array-like
            Field of array to use as the index, alternately a specific set of
            input labels to use

        Returns
        -------
        df : DataFrame
        """
        if indexField is not None:  # pragma: no cover
            warnings.warn("indexField argument is deprecated. Use index "
                          "instead", FutureWarning)
            index = indexField

        columns, sdict = _rec_to_dict(data)

        if exclude is None:
            exclude = set()
        else:
            exclude = set(exclude)

        for col in exclude:
            del sdict[col]
            columns.remove(col)

        if index is not None:
            if isinstance(index, basestring):
                result_index = sdict.pop(index)
                columns.remove(index)
            else:
                try:
                    arrays = []
                    for field in index:
                        arrays.append(sdict[field])
                    for field in index:
                        del sdict[field]
                        columns.remove(field)
                    result_index = MultiIndex.from_arrays(arrays)
                except Exception:
                    result_index = index
        else:
            result_index = np.arange(len(data))

        return cls(sdict, index=result_index, columns=columns)

    def to_records(self, index=True):
        """
        Convert DataFrame to record array. Index will be put in the
        'index' field of the record array if requested

        Parameters
        ----------
        index : boolean, default True
            Include index in resulting record array, stored in 'index' field

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
    def from_csv(cls, path, header=0, delimiter=',', index_col=0):
        """
        Read delimited file into DataFrame

        Parameters
        ----------
        path : string
        header : int, default 0
            Row to use at header (skip prior rows)
        delimiter : string, default ','
        index_col : int or sequence, default 0
            Column to use for index. If a sequence is given, a MultiIndex
            is used.

        Notes
        -----
        Will attempt to convert index to datetimes for time series
        data. Use read_csv for more options

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

        Parameters
        ----------
        fill_value : float, default NaN
        kind : {'block', 'integer'}

        Returns
        -------
        y : SparseDataFrame
        """
        from pandas.core.sparse import SparseDataFrame
        return SparseDataFrame(self._series, index=self.index,
                               default_kind=kind,
                               default_fill_value=fill_value)

    def to_csv(self, path, nanRep='', cols=None, header=True,
              index=True, index_label=None, mode='wb'):
        """
        Write DataFrame to a comma-separated values (csv) file

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
        index_label : string or sequence, default None
            Column label for index column(s) if desired. If None is given, and
            `header` and `index` are True, then the index names are used. A
            sequence should be given if the DataFrame uses MultiIndex.
        mode : Python write mode, default 'wb'
        """
        f = open(path, mode)

        if cols is None:
            cols = self.columns

        series = self._series
        if header:
            joined_cols = ','.join([str(c) for c in cols])
            if index:
                # should write something for index label
                if index_label is None:
                    index_label = getattr(self.index, 'names', ['index'])
                elif not isinstance(index_label, (list, tuple, np.ndarray)):
                    # given a string for a DF with Index
                    index_label = [index_label]
                f.write('%s,%s' % (",".join(index_label), joined_cols))
            else:
                f.write(joined_cols)
            f.write('\n')

        nlevels = getattr(self.index, 'nlevels', 1)
        for idx in self.index:
            if index:
                if nlevels == 1:
                    f.write(str(idx))
                else: # handle MultiIndex
                    f.write(",".join([str(i) for i in idx]))
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

    def to_string(self, buf=None, columns=None, colSpace=None,
                  nanRep='NaN', formatters=None, float_format=None,
                  sparsify=True):
        from pandas.core.common import _format, adjoin
        import sys

        if buf is None:  # pragma: no cover
            buf = sys.stdout

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

        def _format_col(col):
            formatter = formatters.get(col, _myformat)
            return [formatter(x) for x in self[col]]

        if columns is None:
            columns = self.columns
        else:
            columns = [c for c in columns if c in self]

        to_write = []

        if len(columns) == 0 or len(self.index) == 0:
            to_write.append('Empty %s' % type(self).__name__)
            to_write.append(repr(self.index))
        else:
            (str_index,
             str_columns) = self._get_formatted_labels(sparsify=sparsify)
            stringified = [str_columns[i] + _format_col(c)
                           for i, c in enumerate(columns)]
            to_write.append(adjoin(1, str_index, *stringified))

        for s in to_write:
            if isinstance(s, unicode):
                to_write = [unicode(s) for s in to_write]
                break

        for s in to_write:
            print >> buf, s

    def _get_formatted_labels(self, sparsify=True):
        from pandas.core.index import _sparsify

        if isinstance(self.index, MultiIndex):
            fmt_index = self.index.format(sparsify=sparsify)
        else:
            fmt_index = self.index.format()

        if isinstance(self.columns, MultiIndex):
            fmt_columns = self.columns.format(sparsify=False, adjoin=False)
            str_columns = zip(*[[' %s' % y for y in x]
                                for x in zip(*fmt_columns)])
            if sparsify:
                str_columns = _sparsify(str_columns)

            str_columns = [list(x) for x in zip(*str_columns)]
            str_index = [''] * self.columns.nlevels + fmt_index
        else:
            str_columns = [[' %s' % x] for x in self.columns.format()]
            str_index = [''] + fmt_index

        return str_index, str_columns

    def info(self, verbose=True, buf=None):
        """
        Concise summary of a DataFrame, used in __repr__ when very large.

        Parameters
        ----------
        verbose : boolean, default True
            If False, don't print column count summary
        buf : writable buffer, defaults to sys.stdout
        """
        import sys
        if buf is None:  # pragma: no cover
            buf = sys.stdout

        print >> buf, str(type(self))
        print >> buf, self.index.summary()

        if len(self.columns) == 0:
            print >> buf, 'Empty %s' % type(self).__name__
            return

        cols = self.columns

        if verbose:
            print >> buf, unicode('Data columns:')
            space = max([len(_stringify(k)) for k in self.columns]) + 4
            col_counts = []
            counts = self.count()
            assert(len(cols) == len(counts))
            for col, count in counts.iteritems():
                colstr = _stringify(col)
                col_counts.append('%s%d  non-null values' %
                                  (_put_str(colstr, space), count))
            print >> buf, unicode('\n'.join(col_counts))
        else:
            if len(cols) <= 2:
                print >> buf, unicode('Columns: %s' % repr(cols))
            else:
                print >> buf, unicode('Columns: %s to %s'
                                      % (_stringify(cols[0]),
                                         _stringify(cols[-1])))

        counts = self.get_dtype_counts()
        dtypes = ['%s(%d)' % k for k in sorted(counts.iteritems())]
        buf.write(u'dtypes: %s' % ', '.join(dtypes))

    @property
    def dtypes(self):
        return self.apply(lambda x: x.dtype)

    def get_dtype_counts(self):
        counts = {}
        for _, series in self.iteritems():
            if series.dtype in counts:
                counts[series.dtype] += 1
            else:
                counts[series.dtype] = 1

        return Series(counts)

    #----------------------------------------------------------------------
    # properties for index and columns

    def _get_columns(self):
        return self._data.axes[0]

    def _set_columns(self, value):
        self._data.set_axis(0, value)
        self._series_cache.clear()
    columns = property(fset=_set_columns, fget=_get_columns)

    # reference underlying BlockManager
    index = AxisProperty(1)

    def as_matrix(self, columns=None):
        """
        Convert the frame to its Numpy-array matrix representation. Columns
        are presented in sorted order unless a specific list of columns is
        provided.

        Parameters
        ----------
        columns : array-like
            Specific column order

        Returns
        -------
        values : ndarray
            If the DataFrame is heterogeneous and contains booleans or objects,
            the result will be of dtype=object
        """
        self._consolidate_inplace()
        return self._data.as_matrix(columns).T

    values = property(fget=as_matrix)

    def transpose(self):
        """
        Returns a DataFrame with the rows/columns switched. If the DataFrame is
        homogeneously-typed, the data is not copied
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
        elif isinstance(state[0], dict):  # pragma: no cover
            self._unpickle_frame_compat(state)
        else:  # pragma: no cover
            # old pickling format, for compatibility
            self._unpickle_matrix_compat(state)
        self._series_cache = {}

    # legacy pickle formats
    def _unpickle_frame_compat(self, state):  # pragma: no cover
        from pandas.core.common import _unpickle_array
        if len(state) == 2:  # pragma: no cover
            series, idx = state
            columns = sorted(series)
        else:
            series, cols, idx = state
            columns = _unpickle_array(cols)

        index = _unpickle_array(idx)
        self._data = self._init_dict(series, index, columns, None)

    def _unpickle_matrix_compat(self, state):  # pragma: no cover
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
        return self._constructor(result, index=self.index,
                                 columns=self.columns, copy=False)

    #----------------------------------------------------------------------
    # getitem/setitem related

    def __getitem__(self, key):
        # slice rows
        if isinstance(key, slice):
            new_data = self._data.get_slice(key, axis=1)
            return self._constructor(new_data)
        # either boolean or fancy integer index
        elif isinstance(key, np.ndarray):
            if len(key) != len(self.index):
                raise ValueError('Item wrong length %d instead of %d!' %
                                 (len(key), len(self.index)))

            # also raises Exception if object array with NA values
            if _is_bool_indexer(key):
                key = np.asarray(key, dtype=bool)

            new_index = self.index[key]
            return self.reindex(new_index)
        elif isinstance(self.columns, MultiIndex):
            return self._getitem_multilevel(key)
        else:
            return self._getitem_single(key)

    def _slice(self, slobj, axis=0):
        if axis == 0:
            mgr_axis = 1
        else:
            mgr_axis = 0

        new_data = self._data.get_slice(slobj, axis=mgr_axis)
        return self._constructor(new_data)

    def _getitem_multilevel(self, key):
        loc = self.columns.get_loc(key)
        if isinstance(loc, (slice, np.ndarray)):
            new_columns = self.columns[loc]
            result_columns = _maybe_droplevels(new_columns, key)
            if self._is_mixed_type:
                result = self.reindex(columns=new_columns)
                result.columns = result_columns
            else:
                new_values = self.values[:, loc]
                result = DataFrame(new_values, index=self.index,
                                   columns=result_columns)
            return result
        else:
            return self._getitem_single(key)

    def _getitem_single(self, key):
        res = self._series_cache.get(key)
        if res is not None:
            return res

        values = self._data.get(key)
        res = Series(values, index=self.index)
        self._series_cache[key] = res
        return res

    def __setitem__(self, key, value):
        # support boolean setting with DataFrame input, e.g.
        # df[df > df2] = 0
        if isinstance(key, DataFrame):
            if not (key.index.equals(self.index) and
                    key.columns.equals(self.columns)):
                raise PandasError('Can only index with like-indexed '
                                  'DataFrame objects')

            self._boolean_set(key, value)
        else:
            # set column
            self._set_item(key, value)

    def _boolean_set(self, key, value):
        mask = key.values
        if mask.dtype != np.bool_:
            raise ValueError('Must pass DataFrame with boolean values only')

        if self._data.is_mixed_dtype():
            raise ValueError('Cannot do boolean setting on mixed-type frame')

        if isinstance(value, DataFrame):
            assert(value._indexed_same(self))
            np.putmask(self.values, mask, value.values)
        else:
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
        value = np.atleast_2d(value)
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
        value = np.atleast_2d(value)
        self._data.set(key, value)

        try:
            del self._series_cache[key]
        except KeyError:
            pass

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

        try:
            del self._series_cache[key]
        except KeyError:
            pass

    def pop(self, item):
        """
        Return column and drop from frame. Raise KeyError if not found.

        Returns
        -------
        column : Series
        """
        result = self[item]
        del self[item]
        return result

    # to support old APIs
    @property
    def _series(self):
        return self._data.get_series_dict()

    def xs(self, key, axis=0, copy=True):
        """
        Returns a cross-section (row or column) from the DataFrame as a Series
        object. Defaults to returning a row (axis 0)

        Parameters
        ----------
        key : object
            Some label contained in the index, or partially in a MultiIndex
        axis : int, default 0
            Axis to retrieve cross-section on
        copy : boolean, default True
            Whether to make a copy of the data

        Returns
        -------
        xs : Series
        """
        if axis == 1:
            data = self[key]
            if copy:
                data = data.copy()
            return data

        self._consolidate_inplace()
        new_data = self._data.xs(key, axis=1, copy=copy)
        if new_data.ndim == 1:
            return Series(new_data.as_matrix(), index=self.columns)
        else:
            result = DataFrame(new_data)
            result.index = _maybe_droplevels(result.index, key)
            return result

    #----------------------------------------------------------------------
    # Reindexing

    def reindex(self, index=None, columns=None, method=None, copy=True):
        """Conform Series to new index with optional filling logic, placing
        NA/NaN in locations having no value in the previous index. A new object
        is produced unless the new index is equivalent to the current one and
        copy=False

        Parameters
        ----------
        index : array-like, optional
            New labels / index to conform to. Preferably an Index object to
            avoid duplicating data
        columns : array-like, optional
            Same usage as index argument
        method : {'backfill', 'bfill', 'pad', 'ffill', None}, default None
            Method to use for filling holes in reindexed DataFrame
            pad / ffill: propagate last valid observation forward to next valid
            backfill / bfill: use NEXT valid observation to fill gap
        copy : boolean, default True
            Return a new object, even if the passed indexes are the same

        Examples
        --------
        >>> df.reindex(index=[date1, date2, date3], columns=['A', 'B', 'C'])

        Returns
        -------
        reindexed : same type as calling instance
        """
        self._consolidate_inplace()
        frame = self

        if index is not None:
            index = _ensure_index(index)
            frame = frame._reindex_index(index, method, copy)

        if columns is not None:
            columns = _ensure_index(columns)
            frame = frame._reindex_columns(columns, copy)

        return frame

    def _reindex_index(self, new_index, method, copy):
        if new_index.equals(self.index):
            if copy:
                return self.copy()
            else:
                return self
        new_data = self._data.reindex_axis(new_index, method, axis=1)
        return self._constructor(new_data)

    def _reindex_columns(self, new_columns, copy):
        if new_columns.equals(self.columns):
            if copy:
                return self.copy()
            else:
                return self
        new_data = self._data.reindex_axis(new_columns, axis=0)
        return self._constructor(new_data)

    def reindex_like(self, other, method=None, copy=True):
        """
        Reindex DataFrame to match indices of another DataFrame, optionally
        with filling logic

        Parameters
        ----------
        other : DataFrame
        method : string or None
        copy : boolean, default True

        Notes
        -----
        Like calling s.reindex(index=other.index, columns=other.columns,
                               method=...)

        Returns
        -------
        reindexed : DataFrame
        """
        return self.reindex(index=other.index, columns=other.columns,
                            method=method, copy=copy)

    def take(self, indices, axis=0):
        """
        Analogous to ndarray.take, return DataFrame corresponding to requested
        indices along an axis

        Parameters
        ----------
        indices : list / array of ints
        axis : {0, 1}

        Returns
        -------
        taken : DataFrame
        """
        if self._data.is_mixed_dtype():
            if axis == 0:
                new_data = self._data.take(indices, axis=1)
                return DataFrame(new_data)
            else:
                new_columns = self.columns.take(indices)
                return self.reindex(columns=new_columns)
        else:
            new_values = self.values.take(indices, axis=axis)
            if axis == 0:
                new_columns = self.columns
                new_index = self.index.take(indices)
            else:
                new_columns = self.columns.take(indices)
                new_index = self.index
            return DataFrame(new_values, index=new_index,
                             columns=new_columns)

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
        Arguments are mutually exclusive, but this is not checked for

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
        Return object with labels on given axis omitted where alternately any
        or all of the data are missing

        Parameters
        ----------
        axis : {0, 1}
        how : {'any', 'all'}
            any : if any NA values are present, drop that label
            all : if all values are NA, drop that label
        thresh : int, default None
            int value : require that many non-NA values
        subset : array-like
            Labels along other axis to consider, e.g. if you are dropping rows
            these would be a list of columns to include

        Returns
        -------
        dropped : DataFrame
        """
        axis_name = self._get_axis_name(axis)

        if axis == 0:
            agg_axis = 1
        elif axis == 1:
            agg_axis = 0
        else:  # pragma: no cover
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
        Sort DataFrame either by labels (along either axis) or by the values in
        a column

        Parameters
        ----------
        columns : object
            Column name in frame
        ascending : boolean, default True
            Sort ascending vs. descending
        axis : {0, 1}
            Sort index/rows versus columns

        Returns
        -------
        sorted : DataFrame
        """
        by = None
        if column:
            assert(axis == 0)
            by = self[column].values
        return self.sort_index(by=by, axis=axis, ascending=ascending)

    def sort_index(self, axis=0, by=None, ascending=True):
        """
        Sort DataFrame either by labels (along either axis) or by the values in
        a column

        Parameters
        ----------
        axis : {0, 1}
            Sort index/rows versus columns
        by : object
            Column name in frame
        ascending : boolean, default True
            Sort ascending vs. descending

        Returns
        -------
        sorted : DataFrame
        """
        labels = self._get_axis(axis)

        if by is not None:
            try:
                if by in self.columns:
                    assert(axis == 0)
                by = self[by].values
            except Exception:
                pass

            assert(len(by) == len(labels))
            sort_index = Series(by, index=labels).order().index
        else:
            sort_index = labels.take(labels.argsort())

        if not ascending:
            sort_index = sort_index[::-1]

        if axis == 0:
            return self.reindex(sort_index)
        else:
            return self.reindex(columns=sort_index)

    def sortlevel(self, level=0, axis=0, ascending=True):
        """
        Sort multilevel index by chosen axis and primary level. Data will be
        lexicographically sorted by the chosen level followed by the other
        levels (in order)

        Parameters
        ----------
        level : int
        axis : {0, 1}
        ascending : bool, default True

        Returns
        -------
        sorted : DataFrame
        """
        the_axis = self._get_axis(axis)
        if not isinstance(the_axis, MultiIndex):
            raise Exception('can only sort by level with a hierarchical index')

        new_axis, indexer = the_axis.sortlevel(level, ascending=ascending)

        if self._data.is_mixed_dtype():
            if axis == 0:
                return self.reindex(index=new_axis)
            else:
                return self.reindex(columns=new_axis)

        if axis == 0:
            index = new_axis
            columns = self.columns
        else:
            index = self.index
            columns = new_axis
        new_values = self.values.take(indexer, axis=axis)
        return self._constructor(new_values, index=index, columns=columns)

    def swaplevel(self, i, j, axis=0):
        """
        Swap levels i and j in a MultiIndex on a particular axis

        Returns
        -------
        swapped : type of caller (new object)
        """
        result = self.copy()

        if axis == 0:
            result.index = result.index.swaplevel(i, j)
        else:
            result.columns = result.columns.swaplevel(i, j)
        return result

    #----------------------------------------------------------------------
    # Filling NA's

    def fillna(self, value=None, method='pad'):
        """
        Fill NA/NaN values using the specified method. Member Series /
        TimeSeries are filled separately

        Parameters
        ----------
        method : {'backfill', 'bfill', 'pad', 'ffill', None}, default 'pad'
            Method to use for filling holes in reindexed Series
            pad / ffill: propagate last valid observation forward to next valid
            backfill / bfill: use NEXT valid observation to fill gap
        value : any kind (should be same type as array)
            Value to use to fill holes (e.g. 0)

        See also
        --------
        reindex, asfreq

        Returns
        -------
        filled : DataFrame
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

    def rename(self, index=None, columns=None, copy=True):
        """
        Alter index and / or columns using input function or
        functions. Function / dict values must be unique (1-to-1). Labels not
        contained in a dict / Series will be left as-is.

        Parameters
        ----------
        index : dict-like or function, optional
            Transformation to apply to index values
        columns : dict-like or function, optional
            Transformation to apply to column values
        copy : boolean, default True
            Also copy underlying data

        See also
        --------
        Series.rename

        Returns
        -------
        renamed : DataFrame (new object)
        """
        if isinstance(index, (dict, Series)):
            def index_f(x):
                if x in index:
                    return index[x]
                else:
                    return x
        else:
            index_f = index

        if isinstance(columns, (dict, Series)):
            def columns_f(x):
                if x in columns:
                    return columns[x]
                else:
                    return x
        else:
            columns_f = columns

        if index is None and columns is None:
            raise Exception('must pass either index or columns')

        self._consolidate_inplace()

        result = self.copy(deep=copy)

        if index is not None:
            result._rename_index_inplace(index_f)

        if columns is not None:
            result._rename_columns_inplace(columns_f)

        return result

    def _rename_index_inplace(self, mapper):
        self._data = self._data.rename_axis(mapper, axis=1)
        self._series_cache.clear()

    def _rename_columns_inplace(self, mapper):
        self._data = self._data.rename_items(mapper, copydata=False)
        self._series_cache.clear()

    #----------------------------------------------------------------------
    # Arithmetic / combination related

    def _combine_frame(self, other, func, fill_value=None):
        new_index = self.index.union(other.index)

        # some shortcuts
        if fill_value is None:
            if not self and not other:
                return self._constructor(index=new_index)
            elif not self:
                return other * nan
            elif not other:
                return self * nan

        need_reindex = False
        new_columns = self.columns.union(other.columns)
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

    def _combine_series(self, other, func, fill_value=None, axis=None):
        if axis is not None:
            axis = self._get_axis_name(axis)
            if axis == 'index':
                return self._combine_match_index(other, func, fill_value)
            else:
                return self._combine_match_columns(other, func, fill_value)
        return self._combine_series_infer(other, func, fill_value)

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
        new_index = self.index.union(other.index)
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
        func : function
        fill_value : scalar value

        Returns
        -------
        result : DataFrame
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

        # sorts if possible
        new_columns = this.columns.union(other.columns)
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

    def combine_first(self, other):
        """
        Combine two DataFrame objects and default to non-null values in frame
        calling the method. Result index will be the union of the two indexes

        Parameters
        ----------
        other : DataFrame

        Examples
        --------
        >>> a.combine_first(b)
            a's values prioritized, use values from b to fill holes

        Returns
        -------
        combined : DataFrame
        """
        combiner = lambda x, y: np.where(isnull(x), y, x)
        return self.combine(other, combiner)

    #----------------------------------------------------------------------
    # Misc methods

    def first_valid_index(self):
        """
        Return label for first non-NA/null value
        """
        return self.index[self.count(1) > 0][0]

    def last_valid_index(self):
        """
        Return label for last non-NA/null value
        """
        return self.index[self.count(1) > 0][-1]

    def head(self, n=5):
        """Returns first n rows of DataFrame
        """
        return self[:n]

    def tail(self, n=5):
        """Returns last n rows of DataFrame
        """
        return self[-n:]

    #----------------------------------------------------------------------
    # Data reshaping

    def pivot(self, index=None, columns=None, values=None):
        """
        Reshape data (produce a "pivot" table) based on column values. Uses
        unique values from index / columns to form axes and return either
        DataFrame or Panel, depending on whether you request a single value
        column (DataFrame) or all columns (Panel)

        Parameters
        ----------
        index : string or object
            Column name to use to make new frame's index
        columns : string or object
            Column name to use to make new frame's columns
        values : string or object, optional
            Column name to use for populating new frame's values

        Notes
        -----
        For finer-tuned control, see hierarchical indexing documentation along
        with the related stack/unstack methods

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
        pivoted : DataFrame
            If no values column specified, will have hierarchically indexed
            columns
        """
        from pandas.core.reshape import pivot
        return pivot(self, index=index, columns=columns, values=values)

    def stack(self, level=-1, dropna=True):
        """
        Convert DataFrame to Series with multi-level Index. Columns become the
        second level of the resulting hierarchical index

        Returns
        -------
        stacked : Series
        """
        from pandas.core.reshape import stack
        return stack(self, level=level, dropna=dropna)

    def unstack(self, level=-1):
        """
        "Unstack" level from MultiLevel index to produce reshaped DataFrame

        Parameters
        ----------
        level : int, default last level
            Level to unstack

        Examples
        --------
        >>> s
        one  a   1.
        one  b   2.
        two  a   3.
        two  b   4.

        >>> s.unstack(level=-1)
             a   b
        one  1.  2.
        two  3.  4.

        >>> s.unstack(level=0)
           one  two
        a  1.   2.
        b  3.   4.

        Returns
        -------
        unstacked : DataFrame
        """
        from pandas.core.reshape import _Unstacker
        unstacker = _Unstacker(self.values, self.index, level=level,
                               value_columns=self.columns)
        return unstacker.get_result()

    def delevel(self):
        """
        For DataFrame with multi-level index, return new DataFrame with
        labeling information in the columns under names 'level_0', 'level_1',
        etc.

        Notes
        -----
        Experimental, subject to API change

        Returns
        -------
        deleveled : DataFrame
        """
        if not isinstance(self.index, MultiIndex):
            raise Exception('this DataFrame does not have a multi-level index')

        new_obj = self.copy()
        names = self.index.names

        zipped = zip(self.index.levels, self.index.labels)
        for i, (lev, lab) in reversed(list(enumerate(zipped))):
            new_obj.insert(0, names[i], np.asarray(lev).take(lab))

        new_obj.index = np.arange(len(new_obj))

        return new_obj

    #----------------------------------------------------------------------
    # Time series-related

    def asfreq(self, freq, method=None):
        """
        Convert all TimeSeries inside to specified frequency using DateOffset
        objects. Optionally provide fill method to pad/backfill missing values.

        Parameters
        ----------
        offset : DateOffset object, or string in {'WEEKDAY', 'EOM'}
            DateOffset object or subclass (e.g. monthEnd)
        method : {'backfill', 'bfill', 'pad', 'ffill', None}
            Method to use for filling holes in reindexed Series
            pad / ffill: propagate last valid observation forward to next valid
            backfill / bfill: use NEXT valid observation to fill methdo

        Returns
        -------
        converted : DataFrame
        """
        if len(self.index) == 0:
            return self.copy()

        if isinstance(freq, datetools.DateOffset):
            dateRange = DateRange(self.index[0], self.index[-1], offset=freq)
        else:
            dateRange = DateRange(self.index[0], self.index[-1], time_rule=freq)

        return self.reindex(dateRange, method=method)

    def diff(self, periods=1):
        """
        1st discrete difference of object

        Parameters
        ----------
        periods : int, default 1
            Periods to shift for forming difference

        Returns
        -------
        diffed : DataFrame
        """
        return self - self.shift(periods)

    def shift(self, periods, offset=None, **kwds):
        """
        Shift the index of the DataFrame by desired number of periods with an
        optional time offset

        Parameters
        ----------
        periods : int
            Number of periods to move, can be positive or negative
        offset : DateOffset, timedelta, or time rule string, optional
            Increment to use from datetools module or time rule (e.g. 'EOM')

        Returns
        -------
        shifted : DataFrame
        """
        if periods == 0:
            return self

        offset = kwds.get('timeRule', offset)
        if isinstance(offset, basestring):
            offset = datetools.getOffset(offset)

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
        Applies function along input axis of DataFrame. Objects passed to
        functions are Series objects having index either the DataFrame's index
        (axis=0) or the columns (axis=1). Returns either a DataFrame (if the
        function produces another Series) or a Series indexed on either the
        index or columns if the function produces an aggregated value.

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
        >>> df.apply(numpy.sum, axis=0) # equiv to df.sum(0)
        >>> df.apply(numpy.sum, axis=1) # equiv to df.sum(1)

        Notes
        -----
        Functions should not alter the index of the Series passed to them

        Returns
        -------
        applied : Series or DataFrame
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

    def applymap(self, func):
        """
        Apply a function to a DataFrame that is intended to operate
        elementwise, i.e. like doing map(func, series) for each series in the
        DataFrame

        Parameters
        ----------
        func : function
            Python function, returns a single value from a single value

        Returns
        -------
        applied : DataFrame
        """
        npfunc = np.frompyfunc(func, 1, 1)

        def f(x):
            result = npfunc(x)
            try:
                result = result.astype(x.dtype)
            except Exception:
                pass
            return result

        return self.apply(f)

    #----------------------------------------------------------------------
    # Merging / joining methods

    def append(self, other, ignore_index=False):
        """
        Append columns of other to end of this frame's columns and index.
        Columns not in this frame are added as new columns.

        Parameters
        ----------
        other : DataFrame
        ignore_index : boolean, default False
            If True do not use the index labels. Useful for gluing together
            record arrays

        Returns
        -------
        appended : DataFrame
        """
        if not other:
            return self.copy()
        if not self:
            return other.copy()

        if ignore_index:
            new_index = None
        else:
            new_index = np.concatenate((self.index, other.index))

        if self.columns.equals(other.columns):
            return self._append_same_columns(other, new_index)
        else:
            return self._append_different_columns(other, new_index)

    def _append_different_columns(self, other, new_index):
        new_columns = self.columns + other.columns
        new_data = self._append_column_by_column(other)
        return self._constructor(data=new_data, index=new_index,
                                 columns=new_columns)

    def _append_same_columns(self, other, new_index):
        if self._is_mixed_type:
            new_data = self._append_column_by_column(other)
        else:
            new_data= np.concatenate((self.values, other.values), axis=0)
        return self._constructor(new_data, index=new_index,
                                 columns=self.columns)

    def _append_column_by_column(self, other):
        def _concat_missing(values, n):
            values = _maybe_upcast(values)
            missing_values = np.empty(n, dtype=values.dtype)
            missing_values.fill(np.nan)
            return values, missing_values

        new_data = {}
        for col in self:
            values = self._get_raw_column(col)
            if col in other:
                other_values = other._get_raw_column(col)
            else:
                values, other_values = _concat_missing(values, len(other))
            new_data[col] = np.concatenate((values, other_values))

        for col in other:
            values = other._get_raw_column(col)
            if col not in self:
                values, missing_values = _concat_missing(values, len(self))
                new_data[col] = np.concatenate((missing_values, values))

        return new_data

    def _get_raw_column(self, col):
        return self._data.get(col)

    def join(self, other, on=None, how=None, lsuffix='', rsuffix=''):
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
        lsuffix : string
            Suffix to use from left frame's overlapping columns
        rsuffix : string
            Suffix to use from right frame's overlapping columns

        Returns
        -------
        joined : DataFrame
        """
        if on is not None:
            if how is not None:
                raise Exception('how parameter is not valid when '
                                '*on* specified')
            return self._join_on(other, on, lsuffix, rsuffix)
        else:
            if how is None:
                how = 'left'
            return self._join_index(other, how, lsuffix, rsuffix)

    def _join_on(self, other, on, lsuffix, rsuffix):
        if len(other.index) == 0:
            return self

        new_data = self._data.join_on(other._data, self[on], axis=1,
                                      lsuffix=lsuffix, rsuffix=rsuffix)
        return self._constructor(new_data)

    def _join_index(self, other, how, lsuffix, rsuffix):
        from pandas.core.internals import join_managers

        thisdata, otherdata = self._data._maybe_rename_join(
            other._data, lsuffix, rsuffix, copydata=False)

        # this will always ensure copied data
        merged_data = join_managers(thisdata, otherdata, axis=1, how=how)
        return self._constructor(merged_data)

    #----------------------------------------------------------------------
    # Statistical methods, etc.

    def corr(self):
        """
        Compute pairwise correlation of columns, excluding NA/null values

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
        mask = np.isfinite(mat)
        for i, A in enumerate(mat):
            if not mask[i].all():
                for j, B in enumerate(mat):
                    in_common = mask[i] & mask[j]
                    if in_common.any():
                        ac, bc = A[in_common], B[in_common]
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
        axis : {0, 1}
            0 to compute column-wise, 1 for row-wise
        drop : boolean, default False
            Drop missing indices from result, default returns union of all

        Returns
        -------
        correls : Series
        """
        this = self._get_numeric_data()
        other = other._get_numeric_data()

        com_index = this._intersect_index(other)
        com_cols = this._intersect_columns(other)

        # feels hackish
        if axis == 0:
            result_index = com_index
            if not drop:
                result_index = this.columns.union(other.columns)
        else:
            result_index = com_cols
            if not drop:
                result_index = this.index.union(other.index)

        left = this.reindex(index=com_index, columns=com_cols)
        right = other.reindex(index=com_index, columns=com_cols)

        # mask missing values
        left = left + right * 0
        right = right + left * 0

        if axis == 1:
            left = left.T
            right = right.T

        # demeaned data
        ldem = left - left.mean()
        rdem = right - right.mean()

        num = (ldem * rdem).sum()
        dom = (left.count() - 1) * left.std() * right.std()

        correl = num / dom

        if not drop:
            correl = correl.reindex(result_index)

        return correl

    def describe(self):
        """
        Generate various summary statistics of each column, excluding NaN
        values. These include: count, mean, std, min, max, and 10%/50%/90%
        quantiles

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

    def count(self, axis=0, level=None, numeric_only=False):
        """
        Return Series with number of non-NA/null observations over requested
        axis. Works with non-floating point data as well (detects NaN and None)

        Parameters
        ----------
        axis : {0, 1}
            0 for row-wise, 1 for column-wise
        level : int, default None
            If the axis is a MultiIndex (hierarchical), count along a
            particular level, collapsing into a DataFrame
        numeric_only : boolean, default False
            Include only float, int, boolean data

        Returns
        -------
        count : Series (or DataFrame if level specified)
        """
        if level is not None:
            return self._count_level(level, axis=axis,
                                     numeric_only=numeric_only)

        y, axis_labels = self._get_agg_data(axis, numeric_only=numeric_only,
                                            copy=False)
        mask = notnull(y)
        return Series(mask.sum(axis), index=axis_labels)

    def _count_level(self, level, axis=0, numeric_only=False):
        # TODO: deal with sortedness??
        obj = self.sortlevel(level, axis=axis)
        axis_index = obj._get_axis(axis)
        y, _ = self._get_agg_data(axis, numeric_only=numeric_only)
        mask = notnull(y)

        level_index = axis_index.levels[level]

        if len(self) == 0:
            return DataFrame(np.zeros((len(level_index),
                                       len(self.columns)), dtype=int),
                             index=level_index, columns=self.columns)

        n = len(level_index)
        locs = axis_index.labels[level].searchsorted(np.arange(n))

        # WORKAROUND: reduceat fusses about the endpoints. should file ticket?
        start = locs.searchsorted(0, side='right') - 1
        end = locs.searchsorted(len(mask), side='left')

        if axis == 0:
            index = level_index
            columns = self.columns
            result = np.zeros((n, len(self.columns)), dtype=int)
            out = result[start:end]
            np.add.reduceat(mask, locs[start:end], axis=axis, out=out)
        else:
            index = self.index
            columns = level_index
            result = np.zeros((len(self.index), n), dtype=int)
            out = result[:, start:end]
            np.add.reduceat(mask, locs[start:end], axis=axis, out=out)

        # WORKAROUND: to see why, try this
        # arr = np.ones((10, 4), dtype=bool)
        # np.add.reduceat(arr, [0, 3, 3, 7, 9], axis=0)

        # this stinks
        if len(locs) > 1:
            workaround_mask = locs[:-1] == locs[1:]
            if axis == 0:
                result[:-1][workaround_mask] = 0
            else:
                result[:, :-1][:, workaround_mask] = 0

        return DataFrame(result, index=index, columns=columns)

    def sum(self, axis=0, numeric_only=False, skipna=True):
        """
        Return sum over requested axis

        Parameters
        ----------
        axis : {0, 1}
            0 for row-wise, 1 for column-wise
        numeric_only : boolean, default False
            Include only float, int, boolean data
        skipna : boolean, default True
            Exclude NA/null values. If an entire row/column is NA, the result
            will be NA

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

        Returns
        -------
        sum : Series
        """
        y, axis_labels = self._get_agg_data(axis, numeric_only=numeric_only)

        if len(axis_labels) == 0:
            return Series([], index=[])

        if y.dtype == np.object_:
            the_sum = y.sum(axis)
        else:
            mask = np.isfinite(y)

            if skipna:
                if not issubclass(y.dtype.type, np.int_):
                    np.putmask(y, -mask, 0)

            the_sum = y.sum(axis)
            the_count = mask.sum(axis)

            ct_mask = the_count == 0
            if ct_mask.any():
                the_sum[ct_mask] = nan

        return Series(the_sum, index=axis_labels)

    def min(self, axis=0, skipna=True):
        """
        Return minimum over requested axis. NA/null values are excluded

        Parameters
        ----------
        axis : {0, 1}
            0 for row-wise, 1 for column-wise
        skipna : boolean, default True
            Exclude NA/null values. If an entire row/column is NA, the result
            will be NA

        Returns
        -------
        min : Series
        """
        values = self.values.copy()
        if skipna:
            np.putmask(values, -np.isfinite(values), np.inf)
        return Series(values.min(axis), index=self._get_agg_axis(axis))

    def max(self, axis=0, skipna=True):
        """
        Return maximum over requested axis. NA/null values are excluded

        Parameters
        ----------
        axis : {0, 1}
            0 for row-wise, 1 for column-wise
        skipna : boolean, default True
            Exclude NA/null values. If an entire row/column is NA, the result
            will be NA

        Returns
        -------
        max : Series
        """
        values = self.values.copy()
        if skipna:
            np.putmask(values, -np.isfinite(values), -np.inf)
        return Series(values.max(axis), index=self._get_agg_axis(axis))

    def prod(self, axis=0, skipna=True):
        """
        Return product over requested axis. NA/null values are treated as 1

        Parameters
        ----------
        axis : {0, 1}
            0 for row-wise, 1 for column-wise
        skipna : boolean, default True
            Exclude NA/null values. If an entire row/column is NA, the result
            will be NA

        Returns
        -------
        product : Series
        """
        y = np.array(self.values, subok=True)
        if skipna:
            if not issubclass(y.dtype.type, np.int_):
                y[np.isnan(y)] = 1
        result = y.prod(axis)
        count = self.count(axis)
        result[count == 0] = nan
        return Series(result, index=self._get_agg_axis(axis))

    product = prod

    def mean(self, axis=0, skipna=True):
        """
        Return mean over requested axis. NA/null values are excluded

        Parameters
        ----------
        axis : {0, 1}
            0 for row-wise, 1 for column-wise
        skipna : boolean, default True
            Exclude NA/null values. If an entire row/column is NA, the result
            will be NA

        Returns
        -------
        mean : Series
        """
        summed = self.sum(axis, numeric_only=True, skipna=skipna)
        count = self.count(axis, numeric_only=True).astype(float)
        return summed / count

    def quantile(self, q=0.5, axis=0):
        """
        Return values at the given quantile over requested axis, a la
        scoreatpercentile in scipy.stats

        Parameters
        ----------
        q : quantile, default 0.5 (50% quantile)
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

    def median(self, axis=0, skipna=True):
        """
        Return median over requested axis, NA/null are exluded

        Parameters
        ----------
        axis : {0, 1}
            0 for row-wise, 1 for column-wise
        skipna : boolean, default True
            Exclude NA/null values. If an entire row/column is NA, the result
            will be NA

        Returns
        -------
        Series or TimeSeries
        """
        if axis == 0:
            med = [self[col].median(skipna=skipna) for col in self.columns]
            return Series(med, index=self.columns)
        elif axis == 1:
            med = [self.xs(k).median(skipna=skipna) for k in self.index]
            return Series(med, index=self.index)
        else:
            raise Exception('Must have 0<= axis <= 1')

    def mad(self, axis=0, skipna=True):
        """
        Return mean absolute deviation over requested axis

        Parameters
        ----------
        axis : {0, 1}
            0 for row-wise, 1 for column-wise
        skipna : boolean, default True
            Exclude NA/null values. If an entire row/column is NA, the result
            will be NA

        Returns
        -------
        mad : Series
        """
        if axis == 0:
            demeaned = self - self.mean(axis=0)
        else:
            demeaned = self.sub(self.mean(axis=1), axis=0)
        return np.abs(demeaned).mean(axis=axis, skipna=skipna)

    def var(self, axis=0, skipna=True):
        """
        Return unbiased variance over requested axis

        Parameters
        ----------
        axis : {0, 1}
            0 for row-wise, 1 for column-wise
        skipna : boolean, default True
            Exclude NA/null values. If an entire row/column is NA, the result
            will be NA

        Returns
        -------
        var : Series
        """
        y, axis_labels = self._get_agg_data(axis, numeric_only=True)

        mask = np.isnan(y)
        count = (y.shape[axis] - mask.sum(axis)).astype(float)

        if skipna:
            np.putmask(y, mask, 0)

        X = y.sum(axis)
        XX = (y ** 2).sum(axis)

        theVar = (XX - X ** 2 / count) / (count - 1)

        return Series(theVar, index=axis_labels)

    def std(self, axis=0, skipna=True):
        """
        Return unbiased std deviation over requested axis

        Parameters
        ----------
        axis : {0, 1}
            0 for row-wise, 1 for column-wise
        skipna : boolean, default True
            Exclude NA/null values. If an entire row/column is NA, the result
            will be NA

        Returns
        -------
        std : Series
        """
        return np.sqrt(self.var(axis=axis, skipna=skipna))

    def skew(self, axis=0, skipna=True):
        """
        Return unbiased skewness over requested axis

        Parameters
        ----------
        axis : {0, 1}
            0 for row-wise, 1 for column-wise
        skipna : boolean, default True
            Exclude NA/null values. If an entire row/column is NA, the result
            will be NA

        Returns
        -------
        skew : Series
        """
        y, axis_labels = self._get_agg_data(axis, numeric_only=True)

        mask = np.isnan(y)
        count = (y.shape[axis] - mask.sum(axis)).astype(float)

        if skipna:
            np.putmask(y, mask, 0)

        A = y.sum(axis) / count
        B = (y ** 2).sum(axis) / count - A ** 2
        C = (y ** 3).sum(axis) / count - A ** 3 - 3 * A * B

        # floating point error
        B = np.where(np.abs(B) < 1e-14, 0, B)
        C = np.where(np.abs(C) < 1e-14, 0, C)

        result = ((np.sqrt((count ** 2 - count)) * C) /
                  ((count - 2) * np.sqrt(B) ** 3))

        result = np.where(B == 0, 0, result)

        return Series(result, index=axis_labels)

    def _get_agg_data(self, axis, numeric_only=True, copy=True):
        num_cols = self._get_numeric_columns()

        if len(num_cols) < len(self.columns) and numeric_only:
            y = self.as_matrix(num_cols)
            if axis == 0:
                axis_labels = num_cols
            else:
                axis_labels = self.index
        else:
            y = self.values
            if copy:
                y = y.copy()
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

    def _get_numeric_data(self):
        if self._is_mixed_type:
            return self.ix[:, self._get_numeric_columns()]
        else:
            if self.values.dtype != np.object_:
                return self
            else:
                return self.ix[:, []]

    def clip(self, upper=None, lower=None):
        """
        Trim values at input threshold(s)

        Parameters
        ----------
        lower : float, default None
        upper : float, default None

        Returns
        -------
        clipped : DataFrame
        """
        return self.apply(lambda x: x.clip(lower=lower, upper=upper))

    def clip_upper(self, threshold):
        """
        Trim values above threshold

        Returns
        -------
        clipped : DataFrame
        """
        return self.apply(lambda x: x.clip_upper(threshold))

    def clip_lower(self, threshold):
        """
        Trim values below threshold

        Returns
        -------
        clipped : DataFrame
        """
        return self.apply(lambda x: x.clip_lower(threshold))

    #----------------------------------------------------------------------
    # Plotting

    def plot(self, subplots=False, sharex=True, sharey=False, use_index=True,
             figsize=None, grid=True, **kwds):  # pragma: no cover
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
                                   sharex=sharex, sharey=sharey,
                                   figsize=figsize)
        else:
            fig = plt.figure(figsize=figsize)
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

            ax.grid(grid)

        # try to make things prettier
        try:
            fig = plt.gcf()
            fig.autofmt_xdate()
        except Exception:
            pass

        plt.draw_if_interactive()

    def hist(self, grid=True, **kwds):  # pragma: no cover
        """
        Draw Histogram the DataFrame's series using matplotlib / pylab.

        Parameters
        ----------
        kwds : other plotting keyword arguments
            To be passed to hist function
        """
        import matplotlib.pyplot as plt

        n = len(self.columns)
        k = 1
        while k ** 2 < n:
            k += 1
        _, axes = plt.subplots(nrows=k, ncols=k)

        for i, col in enumerate(_try_sort(self.columns)):
            ax = axes[i / k][i % k]
            ax.hist(self[col].dropna().values, **kwds)
            ax.set_title(col)
            ax.grid(grid)

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

    def toDataMatrix(self):  # pragma: no cover
        warnings.warn("toDataMatrix will disappear in next release "
                      "as there is no longer a DataMatrix class",
                      FutureWarning)
        return self.copy()

    def rows(self):  # pragma: no cover
        """Alias for the frame's index"""
        warnings.warn("Replace usage of .rows() with .index, will be removed "
                      "in next release", FutureWarning)
        return self.index

    def cols(self):  # pragma: no cover
        """Return sorted list of frame's columns"""
        warnings.warn("Replace usage of .cols() with .columns, will be "
                      "removed in next release", FutureWarning)
        return list(self.columns)

    def asMatrix(self, *args, **kwargs):  # pragma: no cover
        warnings.warn("asMatrix is deprecated. Use 'as_matrix' or .values "
                      "instead", FutureWarning)
        return self.as_matrix(*args, **kwargs)

    @classmethod
    def fromRecords(cls, *args, **kwargs):  # pragma: no cover
        warnings.warn("fromRecords is deprecated. Use 'from_records' "
                      "instead", FutureWarning)
        return cls.from_records(*args, **kwargs)

    @classmethod
    def fromcsv(cls, *args, **kwargs):  # pragma: no cover
        warnings.warn("fromcsv is deprecated. Use 'from_csv' "
                      "instead", FutureWarning)
        return cls.from_csv(*args, **kwargs)

    combineFirst = deprecate('combineFirst', combine_first)
    getXS = deprecate('getXS', xs)
    merge = deprecate('merge', join)
    toRecords = deprecate('toRecords', to_records)
    toDict = deprecate('toDict', to_dict)
    toString = deprecate('toString', to_string)
    _firstTimeWithValue = deprecate('_firstTimeWithValue', first_valid_index)
    _lastTimeWithValue = deprecate('_lastTimeWithValue', last_valid_index)
    toCSV = deprecate('toCSV', to_csv)

    def dropEmptyRows(self, specificColumns=None):  # pragma: no cover
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
        warnings.warn("dropEmptyRows is deprecated. Use dropna(how='all')",
                      FutureWarning)
        return self.dropna(axis=0, subset=specificColumns, how='all')

    def dropIncompleteRows(self, specificColumns=None,
                           minObs=None):  # pragma: no cover
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
        warnings.warn("dropEmptyRows is deprecated. Use dropna()",
                      FutureWarning)
        if minObs is None:
            return self.dropna(axis=0, subset=specificColumns, how='any')
        else:
            return self.dropna(axis=0, subset=specificColumns, thresh=minObs)

    def tapply(self, func):  # pragma: no cover
        """
        Apply func to the transposed DataFrame, results as per apply
        """
        warnings.warn("tapply is deprecated. Use apply(f, axis=1)",
                      FutureWarning)
        return self.apply(func, axis=1)

    def tgroupby(self, keyfunc, applyfunc):  # pragma: no cover
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
        warnings.warn("tgroupby is deprecated. Use groupby with axis=1",
                      FutureWarning)
        return self.T.groupby(keyfunc).aggregate(applyfunc).T

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


def extract_index(data):
    def _get_index(v):
        if isinstance(v, Series):
            return v.index
        elif isinstance(v, dict):
            return Index(_try_sort(v))

    index = None
    if len(data) == 0:
        index = NULL_INDEX
    elif len(data) > 0 and index is None:
        have_raw_arrays = _check_data_types(data)

        indexes = []

        # this is still kludgier than I'd like
        if have_raw_arrays:
            lengths = list(set(len(x) for x in data.values()))
            if len(lengths) > 1:
                raise ValueError('arrays must all be same length')
            indexes.append(Index(np.arange(lengths[0])))
        else:
            for v in data.values():
                indexes.append(_get_index(v))

        index = _union_indexes(indexes)

    if len(index) == 0:
        index = NULL_INDEX

    return _ensure_index(index)


def _union_indexes(indexes):
    if len(indexes) == 1:
        index = indexes[0]
    if _any_special_indexes(indexes):
        result = indexes[0]
        for other in indexes[1:]:
            result = result.union(other)
        return result
    else:
        index = indexes[0]
        for other in indexes[1:]:
            if not index.equals(other):
                return Index(_tseries.fast_unique_multiple(indexes))

        return index


def _any_special_indexes(indexes):
    for index in indexes:
        if type(index) != Index:
            return True
    return False


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
        arr = np.asarray(values)
        # NumPy strings are a pain, convert to object
        if issubclass(arr.dtype.type, basestring):
            arr = np.array(values, dtype=object, copy=True)
        values = arr
    else:
        # drop subclass info, do not copy data
        values = np.asarray(values)
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
    if isinstance(arr, np.ndarray):
        columns = list(arr.dtype.names)
        sdict = dict((k, arr[k]) for k in columns)
    elif isinstance(arr, DataFrame):
        columns = list(arr.columns)
        sdict = arr._series
    elif isinstance(arr, dict):
        columns = sorted(arr)
        sdict = arr.copy()
    else:  # pragma: no cover
        raise TypeError('%s' % type(arr))

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
                v = v.reindex(index, copy=False)
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

                # prevent NumPy from casting things to string when it shouldn't
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
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
