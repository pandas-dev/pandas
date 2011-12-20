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

from itertools import izip
from StringIO import StringIO
import csv
import operator
import sys

from numpy import nan
import numpy as np
import numpy.ma as ma

from pandas.core.common import (isnull, notnull, PandasError, _try_sort,
                                _default_index, _stringify, _maybe_upcast,
                                _max_rows, _max_columns)
from pandas.core.daterange import DateRange
from pandas.core.generic import NDFrame, AxisProperty
from pandas.core.index import Index, MultiIndex, NULL_INDEX, _ensure_index
from pandas.core.indexing import _NDFrameIndexer, _maybe_droplevels
from pandas.core.internals import BlockManager, make_block, form_blocks
from pandas.core.series import Series, _is_bool_indexer
from pandas.util import py3compat
from pandas.util.terminal import get_terminal_size
import pandas.core.nanops as nanops
import pandas.core.common as com
import pandas.core.datetools as datetools
import pandas._tseries as lib

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

_stat_doc = """
Return %(name)s over requested axis.
%(na_action)s

Parameters
----------
axis : {0, 1}
    0 for row-wise, 1 for column-wise
skipna : boolean, default True
    Exclude NA/null values. If an entire row/column is NA, the result
    will be NA
level : int, default None
    If the axis is a MultiIndex (hierarchical), count along a
    particular level, collapsing into a DataFrame
%(extras)s
Returns
-------
%(shortname)s : Series (or DataFrame if level specified)
"""

_doc_exclude_na = "NA/null values are excluded"

_numeric_only_doc = """numeric_only : boolean, default None
    Include only float, int, boolean data. If None, will attempt to use
    everything, then use only numeric data
"""

def _add_stat_doc(f, name, shortname, na_action=_doc_exclude_na,
                  extras=''):
    doc = _stat_doc % {'name' : name,
                       'shortname' : shortname,
                       'na_action' : na_action,
                       'extras' : extras}
    f.__doc__ = doc

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
        ...                 columns=['a', 'b', 'c', 'd', 'e'])

        See also
        --------
        DataFrame.from_records: constructor from tuples, also record arrays
        DataFrame.from_dict: from dicts of Series, arrays, or dicts
        DataFrame.from_csv: from CSV files
        DataFrame.from_items: from sequence of (key, value) pairs
        read_csv / read_table / read_clipboard
        """
        if data is None:
            data = {}

        if isinstance(data, DataFrame):
            data = data._data

        if isinstance(data, BlockManager):
            mgr = self._init_mgr(data, index, columns, dtype=dtype, copy=copy)
        elif isinstance(data, dict):
            mgr = self._init_dict(data, index, columns, dtype=dtype)
        elif isinstance(data, ma.MaskedArray):
            mask = ma.getmaskarray(data)
            datacopy = ma.copy(data)
            datacopy[mask] = np.nan
            mgr = self._init_ndarray(datacopy, index, columns, dtype=dtype,
                                     copy=copy)
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
            if len(data) > 0 and isinstance(data[0], (list, tuple)):
                data, columns = _list_to_sdict(data, columns)
                mgr = self._init_dict(data, index, columns, dtype=dtype)
            else:
                mgr = self._init_ndarray(data, index, columns, dtype=dtype,
                                         copy=copy)
        else:
            raise PandasError('DataFrame constructor not properly called!')

        NDFrame.__init__(self, mgr)

    def _init_mgr(self, mgr, index, columns, dtype=None, copy=False):
        if columns is not None:
            mgr = mgr.reindex_axis(columns, axis=0)
        if index is not None:
            mgr = mgr.reindex_axis(index, axis=1)
        # do not copy BlockManager unless explicitly done
        if copy and dtype is None:
            mgr = mgr.copy()
        elif dtype is not None:
            # no choice but to copy
            mgr = mgr.astype(dtype)
        return mgr

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
        if isinstance(values, Series):
            if columns is None and values.name is not None:
                columns = [values.name]
            if index is None:
                index = values.index
            else:
                values = values.reindex(index)

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
        terminal_width, terminal_height = get_terminal_size()
        max_rows = terminal_height if com._max_rows == 0 else com._max_rows
        max_columns = com._max_columns

        if max_columns > 0:
            buf = StringIO()
            if len(self.index) < max_rows and \
                    len(self.columns) <= max_columns:
                self.to_string(buf=buf)
            else:
                self.info(buf=buf, verbose=self._verbose_info)
            return buf.getvalue()
        else:
            if len(self.index) > max_rows:
                buf = StringIO()
                self.info(buf=buf, verbose=self._verbose_info)
                return buf.getvalue()
            else:
                buf = StringIO()
                self.to_string(buf=buf)
                value = buf.getvalue()
                if max([len(l) for l in value.split('\n')]) <= terminal_width:
                    return value
                else:
                    buf = StringIO()
                    self.info(buf=buf, verbose=self._verbose_info)
                    return buf.getvalue()

    def __iter__(self):
        """
        Iterate over columns of the frame.
        """
        return iter(self.columns)

    def iteritems(self):
        """Iterator over (column, series) pairs"""
        return ((k, self[k]) for k in self.columns)

    iterkv = iteritems
    if py3compat.PY3:  # pragma: no cover
        items = iteritems

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
    div = _arith_method(lambda x, y: x / y, 'divide')

    radd = _arith_method(lambda x, y: y + x, 'radd')
    rmul = _arith_method(operator.mul, 'rmultiply')
    rsub = _arith_method(lambda x, y: y - x, 'rsubtract')
    rdiv = _arith_method(lambda x, y: y / x, 'rdivide')

    __add__ = _arith_method(operator.add, '__add__', default_axis=None)
    __sub__ = _arith_method(operator.sub, '__sub__', default_axis=None)
    __mul__ = _arith_method(operator.mul, '__mul__', default_axis=None)
    __truediv__ = _arith_method(operator.truediv, '__truediv__',
                               default_axis=None)
    __floordiv__ = _arith_method(operator.floordiv, '__floordiv__',
                               default_axis=None)
    __pow__ = _arith_method(operator.pow, '__pow__', default_axis=None)

    __radd__ = _arith_method(lambda x, y: y + x, '__radd__', default_axis=None)
    __rmul__ = _arith_method(operator.mul, '__rmul__', default_axis=None)
    __rsub__ = _arith_method(lambda x, y: y - x, '__rsub__', default_axis=None)
    __rtruediv__ = _arith_method(lambda x, y: y / x, '__rtruediv__',
                                default_axis=None)
    __rfloordiv__ = _arith_method(lambda x, y: y // x, '__rfloordiv__',
                               default_axis=None)
    __rpow__ = _arith_method(lambda x, y: y ** x, '__rpow__',
                             default_axis=None)

    # boolean operators
    __and__ = _arith_method(operator.and_, '__and__')
    __or__ = _arith_method(operator.or_, '__or__')
    __xor__ = _arith_method(operator.xor, '__xor__')

    # Python 2 division methods
    if not py3compat.PY3:
        __div__ = _arith_method(operator.div, '__div__', default_axis=None)
        __rdiv__ = _arith_method(lambda x, y: y / x, '__rdiv__', default_axis=None)

    def __neg__(self):
        return self * -1

    # Comparison methods
    __eq__ = comp_method(operator.eq, '__eq__')
    __ne__ = comp_method(operator.ne, '__ne__')
    __lt__ = comp_method(operator.lt, '__lt__')
    __gt__ = comp_method(operator.gt, '__gt__')
    __le__ = comp_method(operator.le, '__le__')
    __ge__ = comp_method(operator.ge, '__ge__')

    def dot(self, other):
        """
        Matrix multiplication with DataFrame objects. Does no data alignment

        Parameters
        ----------
        other : DataFrame

        Returns
        -------
        dot_product : DataFrame
        """
        lvals = self.values
        rvals = other.values
        result = np.dot(lvals, rvals)
        return DataFrame(result, index=self.index, columns=other.columns)

    #----------------------------------------------------------------------
    # IO methods (to / from other formats)

    @classmethod
    def from_dict(cls, data, orient='columns', dtype=None):
        """
        Construct Panel from dict of DataFrame objects

        Parameters
        ----------
        data : dict
            {field : DataFrame}
        intersect : boolean
            Intersect indexes of input DataFrames
        orient : {'columns', 'index'}, default 'items'
            The "orientation" of the data. If the keys of the passed dict
            should be the items of the result panel, pass 'items'
            (default). Otherwise if the columns of the values of the passed
            DataFrame objects should be the items (which in the case of
            mixed-dtype data you should do), instead pass 'minor'


        Returns
        -------
        Panel
        """
        from collections import defaultdict

        orient = orient.lower()
        if orient == 'index':
            # TODO: this should be seriously cythonized
            new_data = defaultdict(dict)
            for index, s in data.iteritems():
                for col, v in s.iteritems():
                    new_data[col][index] = v
            data = new_data
        elif orient != 'columns':  # pragma: no cover
            raise ValueError('only recognize index or columns for orient')

        return DataFrame(data, dtype=dtype)

    def to_dict(self):
        """
        Convert DataFrame to nested dictionary

        Returns
        -------
        result : dict like {column -> {index -> value}}
        """
        return dict((k, v.to_dict()) for k, v in self.iteritems())

    @classmethod
    def from_records(cls, data, index=None, exclude=None, columns=None,
                     names=None):
        """
        Convert structured or record ndarray to DataFrame

        Parameters
        ----------
        data : ndarray (structured dtype) or list of tuples
        index : string, list of fields, array-like
            Field of array to use as the index, alternately a specific set of
            input labels to use
        exclude: sequence, default None
            Columns or fields to exclude
        columns : sequence, default None
            Column names to use, replacing any found in passed data

        Returns
        -------
        df : DataFrame
        """
        import warnings

        if names is not None:  # pragma: no cover
            columns = names
            warnings.warn("'names' parameter to DataFrame.from_records is "
                          "being renamed to 'columns', 'names' will be "
                          "removed in 0.8.0", FutureWarning)

        if isinstance(data, (np.ndarray, DataFrame, dict)):
            columns, sdict = _rec_to_dict(data)
        else:
            sdict, columns = _list_to_sdict(data, columns)

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
    def from_items(cls, items, columns=None, orient='columns'):
        """
        Convert (key, value) pairs to DataFrame. The keys will be the axis
        index (usually the columns, but depends on the specified
        orientation). The values should be arrays or Series

        Parameters
        ----------
        items : sequence of (key, value) pairs
            Values should be arrays or Series
        columns : sequence, optional
            Must be passed in the
        orient : {'columns', 'index'}, default 'items'
            The "orientation" of the data. If the keys of the passed dict
            should be the items of the result panel, pass 'items'
            (default). Otherwise if the columns of the values of the passed
            DataFrame objects should be the items (which in the case of
            mixed-dtype data you should do), instead pass 'minor'

        Returns
        -------
        frame : DataFrame
        """
        keys, values = zip(*items)

        if orient == 'columns':
            cols_to_use = columns if columns is not None else keys
            # iterable may have been consumed
            return DataFrame(dict(zip(keys, values)), columns=cols_to_use)
        elif orient == 'index':
            if columns is None:
                raise ValueError("Must pass columns with orient='index'")

            arr = np.array(values, dtype=object).T
            new_data = dict((k, lib.maybe_convert_objects(v))
                            for k, v in zip(columns, arr))
            return DataFrame(new_data, index=keys, columns=columns)
        elif orient != 'columns':  # pragma: no cover
            raise ValueError('only recognize index or columns for orient')

    @classmethod
    def from_csv(cls, path, header=0, sep=',', index_col=0,
                 parse_dates=True):
        """
        Read delimited file into DataFrame

        Parameters
        ----------
        path : string
        header : int, default 0
            Row to use at header (skip prior rows)
        sep : string, default ','
            Field delimiter
        index_col : int or sequence, default 0
            Column to use for index. If a sequence is given, a MultiIndex
            is used. Different default from read_table
        parse_dates : boolean, default True
            Parse dates. Different default from read_table

        Notes
        -----
        Preferable to use read_table for most general purposes but from_csv
        makes for an easy roundtrip to and from file, especially with a
        DataFrame of time series data

        Returns
        -------
        y : DataFrame
        """
        from pandas.io.parsers import read_table
        return read_table(path, header=header, sep=sep,
                          parse_dates=parse_dates, index_col=index_col)

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

    def to_csv(self, path, sep=",", na_rep='', cols=None, header=True,
              index=True, index_label=None, mode='w', nanRep=None):
        """
        Write DataFrame to a comma-separated values (csv) file

        Parameters
        ----------
        path : string
            File path
        nanRep : string, default ''
            Missing data rep'n
        cols : sequence, optional
            Columns to write
        header : boolean, default True
            Write out column names
        index : boolean, default True
            Write row names (index)
        index_label : string or sequence, default None
            Column label for index column(s) if desired. If None is given, and
            `header` and `index` are True, then the index names are used. A
            sequence should be given if the DataFrame uses MultiIndex.
        mode : Python write mode, default 'w'
        sep : character, default ","
            Field delimiter for the output file.
        """
        f = open(path, mode)
        csvout = csv.writer(f, lineterminator='\n', delimiter=sep)

        if nanRep is not None:  # pragma: no cover
            import warnings
            warnings.warn("nanRep is deprecated, use na_rep",
                          FutureWarning)
            na_rep = nanRep

        if cols is None:
            cols = self.columns

        series = self._series
        if header:
            if index:
                # should write something for index label
                if index_label is None:
                    if isinstance(self.index, MultiIndex):
                        index_label = []
                        for i, name in enumerate(self.index.names):
                            if name is None:
                                name = 'level_%d' % i
                            index_label.append(name)
                    else:
                        index_label = self.index.name
                        if index_label is None:
                            index_label = ['index']
                        else:
                            index_label = [index_label]
                elif not isinstance(index_label, (list, tuple, np.ndarray)):
                    # given a string for a DF with Index
                    index_label = [index_label]
                csvout.writerow(list(index_label) + list(cols))
            else:
                csvout.writerow(cols)

        nlevels = getattr(self.index, 'nlevels', 1)
        for idx in self.index:
            row_fields = []
            if index:
                if nlevels == 1:
                    row_fields = [idx]
                else: # handle MultiIndex
                    row_fields = list(idx)
            for i, col in enumerate(cols):
                val = series[col].get(idx)
                if isnull(val):
                    val = na_rep

                row_fields.append(val)

            csvout.writerow(row_fields)

        f.close()

    def to_string(self, buf=None, columns=None, colSpace=None,
                  na_rep='NaN', formatters=None, float_format=None,
                  sparsify=True, nanRep=None, index_names=True):
        """
        Render a DataFrame to a console-friendly tabular output.
        """
        from pandas.core.format import DataFrameFormatter

        if nanRep is not None:  # pragma: no cover
            import warnings
            warnings.warn("nanRep is deprecated, use na_rep",
                          FutureWarning)
            na_rep = nanRep


        formatter = DataFrameFormatter(self, buf=buf, columns=columns,
                                       col_space=colSpace, na_rep=na_rep,
                                       formatters=formatters,
                                       float_format=float_format,
                                       sparsify=sparsify,
                                       index_names=index_names)
        formatter.to_string()

        if buf is None:
            return formatter.buf.getvalue()

    def to_html(self, buf=None, columns=None, colSpace=None,
                na_rep='NaN', formatters=None, float_format=None,
                sparsify=True, index_names=True):
        """
        Render a DataFrame to a html table.
        """
        from pandas.core.format import DataFrameFormatter

        formatter = DataFrameFormatter(self, buf=buf, columns=columns,
                                       col_space=colSpace, na_rep=na_rep,
                                       formatters=formatters,
                                       float_format=float_format,
                                       sparsify=sparsify,
                                       index_names=index_names)
        formatter.to_html()

        if buf is None:
            return formatter.buf.getvalue()

    def info(self, verbose=True, buf=None):
        """
        Concise summary of a DataFrame, used in __repr__ when very large.

        Parameters
        ----------
        verbose : boolean, default True
            If False, don't print column count summary
        buf : writable buffer, defaults to sys.stdout
        """
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

    def convert_objects(self):
        """
        Attempt to infer better dtype for object columns

        Returns
        -------
        converted : DataFrame
        """
        new_data = {}

        # TODO: could be more efficient taking advantage of the block
        for col, s in self.iteritems():
            if s.dtype == np.object_:
                new_data[col] = lib.maybe_convert_objects(s)
            else:
                new_data[col] = s

        return DataFrame(new_data, index=self.index, columns=self.columns)

    def get_dtype_counts(self):
        counts = {}
        for _, series in self.iterkv():
            if series.dtype in counts:
                counts[series.dtype] += 1
            else:
                counts[series.dtype] = 1

        return Series(counts)

    #----------------------------------------------------------------------
    # properties for index and columns


    columns = AxisProperty(0)
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

        # ordinarily created in NDFrame
        self._item_cache = {}

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
    # Array interface

    def __array__(self, dtype=None):
        return self.values

    def __array_wrap__(self, result):
        return self._constructor(result, index=self.index,
                                 columns=self.columns, copy=False)

    #----------------------------------------------------------------------
    # Getting and setting elements

    def get_value(self, index, col):
        """
        Quickly retrieve single value at passed column and index

        Parameters
        ----------
        index : row label
        col : column label

        Returns
        -------
        value : scalar value
        """
        series = self._get_item_cache(col)
        engine = self.index._engine
        return engine.get_value(series, index)

    def set_value(self, index, col, value):
        """
        Put single value at passed column and index

        Parameters
        ----------
        index : row label
        col : column label
        value : scalar value

        Returns
        -------
        frame : DataFrame
            If label pair is contained, will be reference to calling DataFrame,
            otherwise a new object
        """
        try:
            series = self._get_item_cache(col)
            engine = self.index._engine
            engine.set_value(series, index, value)
            return self
        except KeyError:
            new_index, new_columns = self._expand_axes((index, col))
            result = self.reindex(index=new_index, columns=new_columns,
                                  copy=False)
            likely_dtype = com._infer_dtype(value)
            if result[col].dtype != likely_dtype:
                result[col] = result[col].astype(likely_dtype)
            result.set_value(index, col, value)
            return result

    def __getitem__(self, key):
        # slice rows
        if isinstance(key, slice):
            new_data = self._data.get_slice(key, axis=1)
            return self._constructor(new_data)
        # either boolean or fancy integer index
        elif isinstance(key, (np.ndarray, list)):
            if isinstance(key, list):
                key = np.array(key, dtype=object)

            # also raises Exception if object array with NA values
            if _is_bool_indexer(key):
                key = np.asarray(key, dtype=bool)
            return self._getitem_array(key)
        elif isinstance(self.columns, MultiIndex):
            return self._getitem_multilevel(key)
        else:
            return self._get_item_cache(key)

    def _getitem_array(self, key):
        if key.dtype == np.bool_:
            if len(key) != len(self.index):
                raise ValueError('Item wrong length %d instead of %d!' %
                                 (len(key), len(self.index)))

            new_index = self.index[key]
            return self.reindex(new_index)
        else:
            indexer = self.columns.get_indexer(key)
            mask = indexer == -1
            if mask.any():
                raise Exception("No column(s) named: %s" % str(key[mask]))
            return self.reindex(columns=key)

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
            return self._get_item_cache(key)

    def _box_item_values(self, key, values):
        return Series(values, index=self.index, name=key)

    def __getattr__(self, name):
        """After regular attribute access, try looking up the name of a column.
        This allows simpler access to columns for interactive use."""
        if name in self.columns:
            return self[name]
        raise AttributeError("'%s' object has no attribute '%s'" % \
                                                    (type(self).__name__, name))

    def __setitem__(self, key, value):
        # support boolean setting with DataFrame input, e.g.
        # df[df > df2] = 0
        if isinstance(key, DataFrame):
            if not (key.index.equals(self.index) and
                    key.columns.equals(self.columns)):
                raise PandasError('Can only index with like-indexed '
                                  'DataFrame objects')

            self._boolean_set(key, value)
        elif isinstance(key, (np.ndarray, list)):
            return self._set_item_multiple(key, value)
        else:
            # set column
            self._set_item(key, value)

    def _boolean_set(self, key, value):
        mask = key.values
        if mask.dtype != np.bool_:
            raise ValueError('Must pass DataFrame with boolean values only')

        if self._is_mixed_type:
            raise ValueError('Cannot do boolean setting on mixed-type frame')

        if isinstance(value, DataFrame):
            assert(value._indexed_same(self))
            np.putmask(self.values, mask, value.values)
        else:
            self.values[mask] = value

    def _set_item_multiple(self, keys, value):
        if isinstance(value, DataFrame):
            assert(len(value.columns) == len(keys))
            for k1, k2 in zip(keys, value.columns):
                self[k1] = value[k2]
        else:
            self.ix[:, keys] = value

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
        NDFrame._set_item(self, key, value)

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

    def _sanitize_column(self, value):
        # Need to make sure new columns (which go into the BlockManager as new
        # blocks) are always copied
        if _is_sequence(value):
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

    def pop(self, item):
        """
        Return column and drop from frame. Raise KeyError if not found.

        Returns
        -------
        column : Series
        """
        return NDFrame.pop(self, item)

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
        loc = self.index.get_loc(key)
        if np.isscalar(loc):
            new_values = self._data.fast_2d_xs(loc, copy=copy)
            return Series(new_values, index=self.columns, name=key)
        else:
            new_data = self._data.xs(key, axis=1, copy=copy)
            result = DataFrame(new_data)
            result.index = _maybe_droplevels(result.index, key)
            return result

    #----------------------------------------------------------------------
    # Reindexing and alignment

    def align(self, other, join='outer', axis=None, copy=True):
        """
        Align two DataFrame object on their index and columns with the specified
        join method for each axis Index

        Parameters
        ----------
        other : DataFrame or Series
        join : {'outer', 'inner', 'left', 'right'}, default 'outer'
        axis : {0, 1, None}, default None
            Align on index (0), columns (1), or both (None)

        Returns
        -------
        (left, right) : (Series, Series)
            Aligned Series
        """
        if isinstance(other, DataFrame):
            return self._align_frame(other, join=join, axis=axis, copy=copy)
        elif isinstance(other, Series):
            return self._align_series(other, join=join, axis=axis, copy=copy)
        else:  # pragma: no cover
            raise TypeError('unsupported type: %s' % type(other))

    def _align_frame(self, other, join='outer', axis=None, copy=True):
        # defaults
        join_index = self.index
        join_columns = self.columns
        ilidx, iridx = None, None
        clidx, cridx = None, None

        if axis is None or axis == 0:
            if not self.index.equals(other.index):
                join_index, ilidx, iridx = self.index.join(other.index, how=join,
                                                           return_indexers=True)

        if axis is None or axis == 1:
            if not self.columns.equals(other.columns):
                join_columns, clidx, cridx = self.columns.join(other.columns,
                                                               how=join,
                                                               return_indexers=True)

        def _align(frame, row_idx, col_idx):
            new_data = frame._data
            if row_idx is not None:
                new_data = new_data.reindex_indexer(join_index, row_idx, axis=1)

            if col_idx is not None:
                # TODO: speed up on homogeneous DataFrame objects
                new_data = new_data.reindex_items(join_columns)

            if copy and new_data is frame._data:
                new_data = new_data.copy()

            return DataFrame(new_data)

        left = _align(self, ilidx, clidx)
        right = _align(other, iridx, cridx)
        return left, right

    def _align_series(self, other, join='outer', axis=None, copy=True):
        fdata = self._data
        if axis == 0:
            join_index = self.index
            lidx, ridx = None, None
            if not self.index.equals(other.index):
                join_index, lidx, ridx = self.index.join(other.index, how=join,
                                                         return_indexers=True)

            if lidx is not None:
                fdata = fdata.reindex_indexer(join_index, lidx, axis=1)
        elif axis == 1:
            join_index = self.columns
            lidx, ridx = None, None
            if not self.columns.equals(other.index):
                join_index, lidx, ridx = self.columns.join(other.index, how=join,
                                                           return_indexers=True)

            if lidx is not None:
                fdata = fdata.reindex_items(join_index)
        else:
            raise ValueError('Must specify axis=0 or 1')

        if copy and fdata is self._data:
            fdata = fdata.copy()

        left_result = DataFrame(fdata)
        right_result = other if ridx is None else other.reindex(join_index)
        return left_result, right_result

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
                result = self.copy()
                result.index = new_index
                return result
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

    def set_index(self, col_or_cols, drop=True, inplace=False):
        """
        Set the DataFrame index (row labels) using one or more existing
        columns. By default yields a new object.

        Parameters
        ----------
        col_or_cols : column label or list of column labels
        drop : boolean, default True
            Delete columns to be used as the new index
        inplace : boolean, default False
            Modify the DataFrame in place (do not create a new object)

        Returns
        -------
        dataframe : DataFrame
        """
        cols = col_or_cols
        if not isinstance(col_or_cols, (list, tuple)):
            cols = [col_or_cols]

        if inplace:
            frame = self

        else:
            frame = self.copy()

        arrays = []
        for col in cols:
            level = frame[col]
            if drop:
                del frame[col]
            arrays.append(level)

        index = MultiIndex.from_arrays(arrays, names=cols)

        if not index._verify_integrity():
            duplicates = index._get_duplicates()
            raise Exception('Index has duplicate keys: %s' % duplicates)

        # clear up memory usage
        index._cleanup()

        frame.index = index
        return frame

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

    def drop_duplicates(self, col_or_columns=None, take_last=False):
        """
        Return DataFrame with duplicate rows removed, optionally only
        considering certain columns

        Parameters
        ----------
        col_or_columns : column label or sequence of labels, optional
            Only consider certain columns for identifying duplicates, by
            default use all of the columns
        take_last : boolean, default False
            Take the last observed row in a row. Defaults to the first row

        Returns
        -------
        deduplicated : DataFrame
        """
        duplicated = self.duplicated(col_or_columns, take_last=take_last)
        return self[-duplicated]

    def duplicated(self, col_or_columns=None, take_last=False):
        """
        Return boolean Series denoting duplicate rows, optionally only
        considering certain columns

        Parameters
        ----------
        col_or_columns : column label or sequence of labels, optional
            Only consider certain columns for identifying duplicates, by
            default use all of the columns
        take_last : boolean, default False
            Take the last observed row in a row. Defaults to the first row

        Returns
        -------
        duplicated : Series
        """
        if col_or_columns is not None:
            if isinstance(col_or_columns, list):
                keys = zip(*[self[x] for x in col_or_columns])
            else:
                keys = list(self[col_or_columns])
        else:
            keys = zip(*self.values.T)

        duplicated = lib.duplicated(keys, take_last=take_last)
        return Series(duplicated, index=self.index)

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
        return self.sort_index(by=column, axis=axis, ascending=ascending)

    def sort_index(self, axis=0, by=None, ascending=True):
        """
        Sort DataFrame either by labels (along either axis) or by the values in
        a column

        Parameters
        ----------
        axis : {0, 1}
            Sort index/rows versus columns
        by : object
            Column names in frame
        ascending : boolean, default True
            Sort ascending vs. descending

        Returns
        -------
        sorted : DataFrame
        """
        labels = self._get_axis(axis)

        if by is not None:
            assert(axis == 0)
            if isinstance(by, (tuple, list)):
                to_sort = lib.fast_zip([self[x] for x in by])
            else:
                to_sort = self[by].values

            indexer = to_sort.argsort()
        else:
            indexer = labels.argsort()

        if not ascending:
            indexer = indexer[::-1]

        return self.take(indexer, axis=axis)

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
            return self._constructor(new_data)

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
        from pandas.core.series import _get_rename_function

        if index is None and columns is None:
            raise Exception('must pass either index or columns')

        index_f = _get_rename_function(index)
        columns_f = _get_rename_function(columns)

        self._consolidate_inplace()

        result = self.copy(deep=copy)

        if index is not None:
            result._rename_index_inplace(index_f)

        if columns is not None:
            result._rename_columns_inplace(columns_f)

        return result

    def _rename_index_inplace(self, mapper):
        self._data = self._data.rename_axis(mapper, axis=1)
        self._clear_item_cache()

    def _rename_columns_inplace(self, mapper):
        self._data = self._data.rename_items(mapper, copydata=False)
        self._clear_item_cache()

    #----------------------------------------------------------------------
    # Arithmetic / combination related

    def _combine_frame(self, other, func, fill_value=None):
        this, other = self.align(other, join='outer', copy=False)
        new_index, new_columns = this.index, this.columns

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
        left, right = self.align(other, join='outer', axis=0, copy=False)
        if fill_value is not None:
            raise NotImplementedError
        return self._constructor(func(left.values.T, right.values).T,
                                 index=left.index,
                                 columns=self.columns, copy=False)

    def _combine_match_columns(self, other, func, fill_value=None):
        left, right = self.align(other, join='outer', axis=1, copy=False)
        if fill_value is not None:
            raise NotImplementedError

        return self._constructor(func(left.values, right.values),
                                 index=self.index,
                                 columns=left.columns, copy=False)

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

        this, other = self.align(other, axis=0, copy=False)
        new_index = this.index

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
                    arr = com.ensure_float(arr)
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

        Parameters
        ----------
        level : int, string, or list of these, default last level
            Level(s) to stack, can pass level name

        Returns
        -------
        stacked : Series
        """
        from pandas.core.reshape import stack

        if isinstance(level, (tuple, list)):
            result = self
            for lev in level:
                result = stack(result, lev, dropna=dropna)
            return result
        else:
            return stack(self, level, dropna=dropna)

    def unstack(self, level=-1):
        """
        "Unstack" level from MultiLevel index to produce reshaped DataFrame

        Parameters
        ----------
        level : int, string, or list of these, default last level
            Level(s) to unstack, can pass level name

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
        from pandas.core.reshape import unstack
        if isinstance(level, (tuple, list)):
            result = self
            to_unstack = level
            while to_unstack:
                lev = to_unstack[0]
                result = unstack(result, lev)
                to_unstack = [other - 1 if other > lev else other
                              for other in to_unstack[1:]]
            return result
        else:
            return unstack(self, level)

    def delevel(self):
        """
        For DataFrame with multi-level index, return new DataFrame with labeling
        information in the columns under names 'level_0', 'level_1', etc.

        Notes
        -----
        Experimental, subject to API change

        Returns
        -------
        deleveled : DataFrame
        """
        new_obj = self.copy()
        if isinstance(self.index, MultiIndex):
            names = self.index.names
            zipped = zip(self.index.levels, self.index.labels)
            for i, (lev, lab) in reversed(list(enumerate(zipped))):
                col_name = names[i]
                if col_name is None:
                    col_name = 'level_%d' % i

                # to ndarray and maybe infer different dtype
                level_values = lev.values
                level_values = lib.maybe_convert_objects(level_values)
                new_obj.insert(0, col_name, level_values.take(lab))
        else:
            if self.index.name is None:
                raise Exception('Must have name set')
            new_obj.insert(0, self.index.name, self.index.values)
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
            new_values = com.ensure_float(new_values)
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

    def apply(self, func, axis=0, broadcast=False, raw=False,
              args=(), **kwds):
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
        raw : boolean, default False
            If False, convert each row or column into a Series. If raw=True the
            passed function will receive ndarray objects instead. If you are
            just applying a NumPy reduction function this will achieve much
            better performance
        args : tuple
            Positional arguments to pass to function in addition to the
            array/series
        Additional keyword arguments will be passed as keywords to the function

        Examples
        --------
        >>> df.apply(numpy.sqrt) # returns DataFrame
        >>> df.apply(numpy.sum, axis=0) # equiv to df.sum(0)
        >>> df.apply(numpy.sum, axis=1) # equiv to df.sum(1)

        Notes
        -----
        Function passed should not have side effects. If the result is a Series,
        it should have the same index

        Returns
        -------
        applied : Series or DataFrame
        """
        if len(self.columns) == 0 and len(self.index) == 0:
            return self

        if kwds or args and not isinstance(func, np.ufunc):
            f = lambda x: func(x, *args, **kwds)
        else:
            f = func

        if isinstance(f, np.ufunc):
            results = f(self.values)
            return self._constructor(data=results, index=self.index,
                                     columns=self.columns, copy=False)
        else:
            if not broadcast:
                if not all(self.shape):
                    is_reduction = not isinstance(f(_EMPTY_SERIES),
                                                  np.ndarray)
                    if is_reduction:
                        return Series(np.nan, index=self._get_agg_axis(axis))
                    else:
                        return self.copy()

                if raw and not self._is_mixed_type:
                    return self._apply_raw(f, axis)
                else:
                    return self._apply_standard(f, axis)
            else:
                return self._apply_broadcast(f, axis)

    def _apply_raw(self, func, axis):
        try:
            result = lib.reduce(self.values, func, axis=axis)
        except Exception:
            result = np.apply_along_axis(func, axis, self.values)

        # TODO: mixed type case
        if result.ndim == 2:
            return DataFrame(result, index=self.index,
                             columns=self.columns)
        else:
            return Series(result, index=self._get_agg_axis(axis))

    def _apply_standard(self, func, axis, ignore_failures=False):
        try:

            assert(not self._is_mixed_type)  # maybe a hack for now
            values = self.values
            dummy = Series(np.nan, index=self._get_axis(axis),
                           dtype=values.dtype)
            result = lib.reduce(values, func, axis=axis, dummy=dummy)
            return Series(result, index=self._get_agg_axis(axis))
        except Exception:
            pass

        if axis == 0:
            series_gen = ((c, self[c]) for c in self.columns)
            res_index = self.columns
            res_columns = self.index
        elif axis == 1:
            res_index = self.index
            res_columns = self.columns
            series_gen = ((i, Series(v, self.columns))
                          for i, v in izip(self.index, self.values))

        results = {}
        if ignore_failures:
            successes = []
            for i, (k, v) in enumerate(series_gen):
                try:
                    results[k] = func(v)
                    successes.append(i)
                except Exception:
                    pass
            # so will work with MultiIndex, need test
            if len(successes) < len(res_index):
                res_index = res_index.take(successes)
        else:
            for k, v in series_gen:
                results[k] = func(v)

        if len(results) > 0 and _is_sequence(results.values()[0]):
            if not isinstance(results.values()[0], Series):
                index = res_columns
            else:
                index = None

            result = self._constructor(data=results, index=index)

            if axis == 1:
                result = result.T

            return result
        else:
            return Series(results, index=res_index)

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
        return self.apply(lambda x: lib.map_infer(x, func))

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
        if isinstance(other, Series):
            other = other.reindex(self.columns, copy=False)
            other = DataFrame(other.values.reshape((1, len(other))),
                              columns=self.columns)
            if not ignore_index:
                raise Exception('Can only append a Series if ignore_index=True')

        if not other:
            return self.copy()
        if not self:
            return other.copy()

        if ignore_index:
            new_index = None
        else:
            new_index = self.index.append(other.index)
            assert(new_index._verify_integrity())

        if self.columns.equals(other.columns):
            return self._append_same_columns(other, new_index)
        else:
            return self._append_different_columns(other, new_index)

    def _append_different_columns(self, other, new_index):
        indexer = self.columns.get_indexer(other.columns)

        if not (indexer == -1).any():
            new_columns = self.columns
        else:
            new_columns = self.columns.union(other.columns)

        new_data = self._append_column_by_column(other)
        return self._constructor(data=new_data, index=new_index,
                                 columns=new_columns)

    def _append_same_columns(self, other, new_index):
        if self._is_mixed_type:
            new_data = self._append_column_by_column(other)
        else:
            new_data = np.concatenate((self.values, other.values), axis=0)
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
        column.

        Parameters
        ----------
        other : DataFrame, or Series with name field set
            Index should be similar to one of the columns in this one. If a
            Series is passed, its name attribute must be set, and that will be
            used as the column name in the resulting joined DataFrame
        on : column name, tuple/list of column names, or array-like
            Column(s) to use for joining, otherwise join on index. If multiples
            columns given, the passed DataFrame must have a MultiIndex. Can
            pass an array as the join key if not already contained in the
            calling DataFrame. Like an Excel VLOOKUP operation
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
        if how is None:
            how = 'left'
        if on is not None:
            return self._join_on(other, on, how, lsuffix, rsuffix)
        else:
            return self._join_index(other, how, lsuffix, rsuffix)

    def _join_on(self, other, on, how, lsuffix, rsuffix):
        if how not in ('left', 'inner'):  # pragma: no cover
            raise Exception('Only inner / left joins currently supported')

        if isinstance(other, Series):
            assert(other.name is not None)
            other = DataFrame({other.name : other})

        if isinstance(on, (list, tuple)):
            if len(on) == 1:
                join_key = self[on[0]].values
            else:
                join_key = lib.fast_zip([self[k] for k in on])
        elif isinstance(on, np.ndarray) and len(on) == len(self):
            join_key = on
        else:
            join_key = self[on].values

        new_data = self._data.join_on(other._data, join_key, how=how, axis=1,
                                      lsuffix=lsuffix, rsuffix=rsuffix)
        return self._constructor(new_data)

    def _join_index(self, other, how, lsuffix, rsuffix):
        from pandas.core.internals import join_managers

        if isinstance(other, Series):
            assert(other.name is not None)
            other = DataFrame({other.name : other})

        thisdata, otherdata = self._data._maybe_rename_join(
            other._data, lsuffix, rsuffix, copydata=False)

        # this will always ensure copied data
        merged_data = join_managers(thisdata, otherdata, axis=1, how=how)
        return self._constructor(merged_data)

    #----------------------------------------------------------------------
    # Statistical methods, etc.

    def corr(self, method='pearson'):
        """
        Compute pairwise correlation of columns, excluding NA/null values

        Parameters
        ----------
        method : {'pearson', 'kendall', 'spearman'}
            pearson : standard correlation coefficient
            kendall : Kendall Tau correlation coefficient
            spearman : Spearman rank correlation

        Returns
        -------
        y : DataFrame
        """
        cols = self._get_numeric_columns()
        mat = self.as_matrix(cols).T
        corrf = nanops.get_corr_func(method)
        K = len(cols)
        correl = np.empty((K, K), dtype=float)
        mask = np.isfinite(mat)
        for i, ac in enumerate(mat):
            for j, bc  in enumerate(mat):
                valid = mask[i] & mask[j]
                if not valid.all():
                    c = corrf(ac[valid], bc[valid])
                else:
                    c = corrf(ac, bc)
                correl[i, j] = c
                correl[j, i] = c

        return self._constructor(correl, index=cols, columns=cols)

    def cov(self):
        """
        Compute pairwise covariance of columns, excluding NA/null values

        Returns
        -------
        y : DataFrame
        """
        cols = self._get_numeric_columns()
        mat = self.as_matrix(cols).T
        baseCov = np.cov(mat)

        for i, j, ac, bc in self._cov_helper(mat):
            c = np.cov(ac, bc)[0, 1]
            baseCov[i, j] = c
            baseCov[j, i] = c

        return self._constructor(baseCov, index=cols, columns=cols)

    def _cov_helper(self, mat):
        # Get the covariance with items that have NaN values
        mask = np.isfinite(mat)
        for i, A in enumerate(mat):
            if not mask[i].all():
                for j, B in enumerate(mat):
                    in_common = mask[i] & mask[j]
                    if in_common.any():
                        yield i, j, A[in_common], B[in_common]

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

        left, right = this.align(other, join='inner', copy=False)

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
            raxis = 1 if axis == 0 else 0
            result_index = this._get_axis(raxis).union(other._get_axis(raxis))
            correl = correl.reindex(result_index)

        return correl

    def describe(self):
        """
        Generate various summary statistics of each column, excluding NaN
        values. These include: count, mean, std, min, max, and 10%/50%/90%
        quantiles

        Returns
        -------
        DataFrame of summary statistics
        """
        numeric_columns = self._get_numeric_columns()

        if len(numeric_columns) == 0:
            return DataFrame(dict((k, v.describe())
                                  for k, v in self.iteritems()),
                                  columns=self.columns)

        destat_columns = ['count', 'mean', 'std', 'min',
                          '25%', '50%', '75%', 'max']

        destat = []

        for column in numeric_columns:
            series = self[column]
            destat.append([series.count(), series.mean(), series.std(),
                           series.min(), series.quantile(.25), series.median(),
                           series.quantile(.75), series.max()])

        return self._constructor(map(list, zip(*destat)), index=destat_columns, columns=numeric_columns)

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

        if numeric_only:
            frame = self.ix[:, self._get_numeric_columns()]
        else:
            frame = self

        # GH #423
        if len(frame._get_axis(axis)) == 0:
            result = Series(0, index=frame._get_agg_axis(axis))
        else:
            if axis == 1:
                counts = notnull(frame.values).sum(1)
                result = Series(counts, index=frame._get_agg_axis(axis))
            else:
                result = DataFrame.apply(frame, Series.count, axis=axis)

        # what happens with empty DataFrame
        if isinstance(result, DataFrame):
            result = Series({})

        return result

    def _count_level(self, level, axis=0, numeric_only=False):
        if numeric_only:
            frame = self._get_numeric_data()
        else:
            frame = self

        if axis == 1:
            frame = frame.T

        mask = notnull(frame.values).view(np.uint8) # python 2.5

        level_index = frame.index.levels[level]
        counts = lib.count_level_2d(mask, frame.index.labels[level],
                                    len(level_index))

        result = DataFrame(counts, index=level_index,
                           columns=frame.columns)

        if axis == 1:
            return result.T
        else:
            return result

    def sum(self, axis=0, numeric_only=None, skipna=True, level=None):
        if level is not None:
            return self._agg_by_level('sum', axis=axis, level=level,
                                      skipna=skipna)
        return self._reduce(nanops.nansum, axis=axis, skipna=skipna,
                            numeric_only=numeric_only)
    _add_stat_doc(sum, 'sum', 'sum', extras=_numeric_only_doc)

    def mean(self, axis=0, skipna=True, level=None):
        if level is not None:
            return self._agg_by_level('mean', axis=axis, level=level,
                                      skipna=skipna)
        return self._reduce(nanops.nanmean, axis=axis, skipna=skipna,
                            numeric_only=None)
    _add_stat_doc(mean, 'mean', 'mean')

    def min(self, axis=0, skipna=True, level=None):
        if level is not None:
            return self._agg_by_level('min', axis=axis, level=level,
                                      skipna=skipna)
        return self._reduce(nanops.nanmin, axis=axis, skipna=skipna,
                            numeric_only=None)
    _add_stat_doc(min, 'minimum', 'min')

    def max(self, axis=0, skipna=True, level=None):
        if level is not None:
            return self._agg_by_level('max', axis=axis, level=level,
                                      skipna=skipna)
        return self._reduce(nanops.nanmax, axis=axis, skipna=skipna,
                            numeric_only=None)
    _add_stat_doc(max, 'maximum', 'max')

    def prod(self, axis=0, skipna=True, level=None):
        if level is not None:
            return self._agg_by_level('prod', axis=axis, level=level,
                                      skipna=skipna)
        return self._reduce(nanops.nanprod, axis=axis, skipna=skipna,
                            numeric_only=None)
    _add_stat_doc(prod, 'product', 'product',
                  na_action='NA/null values are treated as 1')
    product = prod

    def median(self, axis=0, skipna=True, level=None):
        if level is not None:
            return self._agg_by_level('median', axis=axis, level=level,
                                      skipna=skipna)
        return self._reduce(nanops.nanmedian, axis=axis, skipna=skipna,
                            numeric_only=None)
    _add_stat_doc(median, 'median', 'median')

    def mad(self, axis=0, skipna=True, level=None):
        if level is not None:
            return self._agg_by_level('mad', axis=axis, level=level,
                                      skipna=skipna)

        frame = self._get_numeric_data()

        if axis == 0:
            demeaned = frame - frame.mean(axis=0)
        else:
            demeaned = frame.sub(frame.mean(axis=1), axis=0)
        return np.abs(demeaned).mean(axis=axis, skipna=skipna)
    _add_stat_doc(mad, 'mean absolute deviation', 'mad')

    def var(self, axis=0, skipna=True, level=None):
        if level is not None:
            return self._agg_by_level('var', axis=axis, level=level,
                                      skipna=skipna)
        return self._reduce(nanops.nanvar, axis=axis, skipna=skipna,
                            numeric_only=None)
    _add_stat_doc(var, 'unbiased variance', 'var')

    def std(self, axis=0, skipna=True, level=None):
        if level is not None:
            return self._agg_by_level('std', axis=axis, level=level,
                                      skipna=skipna)
        return np.sqrt(self.var(axis=axis, skipna=skipna))
    _add_stat_doc(std, 'unbiased standard deviation', 'std')

    def skew(self, axis=0, skipna=True, level=None):
        if level is not None:
            return self._agg_by_level('skew', axis=axis, level=level,
                                      skipna=skipna)
        return self._reduce(nanops.nanskew, axis=axis, skipna=skipna,
                            numeric_only=None)
    _add_stat_doc(skew, 'unbiased skewness', 'skew')

    def _reduce(self, op, axis=0, skipna=True, numeric_only=None):
        f = lambda x: op(x, axis=axis, skipna=skipna, copy=True)
        labels = self._get_agg_axis(axis)
        if numeric_only is None:
            try:
                values = self.values
                if not self._is_mixed_type:
                    values = values.copy()
                result = f(values)
            except Exception:
                data = self._get_numeric_data()
                result = f(data.values)
                labels = data._get_agg_axis(axis)
        else:
            if numeric_only:
                data = self._get_numeric_data()
                values = data.values
                labels = data._get_agg_axis(axis)
            else:
                values = self.values
            result = f(values)

        if result.dtype == np.object_:
            try:
                result = result.astype('f8')
            except (ValueError, TypeError):
                pass

        return Series(result, index=labels)

    def idxmin(self, axis=0, skipna=True):
        """
        Return index of first occurence of minimum over requested axis.
        NA/null values are excluded.

        Parameters
        ----------
        axis : {0, 1}
            0 for row-wise, 1 for column-wise
        skipna : boolean, default True
            Exclude NA/null values. If an entire row/column is NA, the result
            will be NA

        Returns
        -------
        idxmin : Series
        """
        indices = nanops.nanargmin(self.values, axis=axis, skipna=skipna)
        index = self._get_axis(axis)
        result = [index[i] if i >= 0 else np.nan for i in indices]
        return Series(result, index=self._get_agg_axis(axis))

    def idxmax(self, axis=0, skipna=True):
        """
        Return index of first occurence of maximum over requested axis.
        NA/null values are excluded.

        Parameters
        ----------
        axis : {0, 1}
            0 for row-wise, 1 for column-wise
        skipna : boolean, default True
            Exclude NA/null values. If an entire row/column is NA, the result
            will be first index.

        Returns
        -------
        idxmax : Series
        """
        indices = nanops.nanargmax(self.values, axis=axis, skipna=skipna)
        index = self._get_axis(axis)
        result = [index[i] if i >= 0 else np.nan for i in indices]
        return Series(result, index=self._get_agg_axis(axis))

    def _agg_by_level(self, name, axis=0, level=0, skipna=True):
        method = getattr(type(self), name)
        applyf = lambda x: method(x, axis=axis, skipna=skipna)
        return self.groupby(level=level, axis=axis).aggregate(applyf)

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

    def rank(self, axis=0):
        """
        Compute numericaldata ranks (1 through n) along axis. Equal values are
        assigned a rank that is the average of the ranks of those values

        Parameters
        ----------
        axis : {0, 1}, default 0
            Ranks over columns (0) or rows (1)

        Returns
        -------
        ranks : DataFrame
        """
        data = self._get_numeric_data()
        ranks = lib.rank_2d_float64(data.values.astype('f8'), axis=axis)
        return DataFrame(ranks, index=data.index, columns=data.columns)

    #----------------------------------------------------------------------
    # Plotting

    def boxplot(self, column=None, by=None, ax=None, fontsize=None,
                rot=0, grid=True, **kwds):
        """
        Make a box plot from DataFrame column/columns optionally grouped
        (stratified) by one or more columns

        Parameters
        ----------
        data : DataFrame
        column : column names or list of names, or vector
            Can be any valid input to groupby
        by : string or sequence
            Column in the DataFrame to group by
        fontsize : int or string

        Returns
        -------
        ax : matplotlib.axes.AxesSubplot
        """
        import pandas.tools.plotting as plots
        import matplotlib.pyplot as plt
        ax = plots.boxplot(self, column=column, by=by, ax=ax, fontsize=fontsize,
                           grid=grid, rot=rot)
        plt.draw_if_interactive()
        return ax

    def plot(self, subplots=False, sharex=True, sharey=False, use_index=True,
             figsize=None, grid=True, legend=True, rot=30, ax=None,
             kind='line', **kwds):
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
        kind : {'line', 'bar'}
        kwds : keywords
            Options to pass to Axis.plot

        Notes
        -----
        This method doesn't make much sense for cross-sections,
        and will error.
        """
        import matplotlib.pyplot as plt

        if subplots:
            fig, axes = plt.subplots(nrows=len(self.columns),
                                   sharex=sharex, sharey=sharey,
                                   figsize=figsize)
        else:
            if ax is None:
                fig = plt.figure(figsize=figsize)
                ax = fig.add_subplot(111)
                axes = [ax]
            else:
                fig = ax.get_figure()

        if kind == 'line':
            if use_index:
                x = self.index
            else:
                x = range(len(self))

            for i, col in enumerate(_try_sort(self.columns)):
                empty = self[col].count() == 0
                y = self[col].values if not empty else np.zeros(x.shape)

                try:
                    if subplots:
                        ax = axes[i]
                        ax.plot(x, y, 'k', label=str(col), **kwds)
                        ax.legend(loc='best')
                    else:
                        ax.plot(x, y, label=str(col), **kwds)
                except Exception, e:
                    msg = ('Unable to plot data %s vs index %s,\n'
                           'error was: %s' % (str(y), str(x), str(e)))
                    raise Exception(msg)

                ax.grid(grid)

            if legend and not subplots:
                ax.legend(loc='best')
        elif kind == 'bar':
            self._bar_plot(axes, subplots=subplots, grid=grid, rot=rot,
                           legend=legend)

        # try to make things prettier
        try:
            fig.autofmt_xdate()
        except Exception:  # pragma: no cover
            pass

        plt.draw_if_interactive()
        if subplots:
            return axes
        else:
            return ax

    def _bar_plot(self, axes, subplots=False, use_index=True, grid=True,
                  rot=30, legend=True, **kwds):
        N, K = self.shape
        xinds = np.arange(N) + 0.25
        colors = 'rgbyk'
        rects = []
        labels = []

        if not subplots:
            ax = axes[0]

        for i, col in enumerate(_try_sort(self.columns)):
            empty = self[col].count() == 0
            y = self[col].values if not empty else np.zeros(len(self))
            if subplots:
                ax = axes[i]
                ax.bar(xinds, y, 0.5,
                       bottom=np.zeros(N), linewidth=1, **kwds)
                ax.set_title(col)
            else:
                rects.append(ax.bar(xinds + i * 0.5/K, y, 0.5/K,
                                    bottom=np.zeros(N),
                                    color=colors[i % len(colors)], **kwds))
                labels.append(col)

        if N < 10:
            fontsize = 12
        else:
            fontsize = 10

        ax.set_xticks(xinds + 0.25)
        ax.set_xticklabels(self.index, rotation=rot, fontsize=fontsize)

        if legend and not subplots:
            fig = ax.get_figure()
            fig.legend([r[0] for r in rects], labels, loc='upper center',
                       fancybox=True, ncol=6, mode='expand')

        import matplotlib.pyplot as plt
        plt.subplots_adjust(top=0.8)

    def hist(self, grid=True, **kwds):
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

        return axes
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


_EMPTY_SERIES = Series([])


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
            return v.keys()

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
    if len(indexes) == 0:
        return Index([])

    if len(indexes) == 1:
        result = indexes[0]
        if isinstance(result, list):
            result = Index(sorted(result))
        return result

    indexes, kind = _sanitize_and_check(indexes)

    if kind == 'special':
        result = indexes[0]
        for other in indexes[1:]:
            result = result.union(other)
        return result
    elif kind == 'array':
        index = indexes[0]
        for other in indexes[1:]:
            if not index.equals(other):
                return Index(lib.fast_unique_multiple(indexes))

        return index
    else:
        return Index(lib.fast_unique_multiple_list(indexes))


def _sanitize_and_check(indexes):
    kinds = list(set([type(index) for index in indexes]))

    if list in kinds:
        if len(kinds) > 1:
            indexes = [Index(_try_sort(x)) if not isinstance(x, Index) else x
                       for x in indexes]
            kinds.remove(list)
        else:
            return indexes, 'list'


    if len(kinds) > 1 or Index not in kinds:
        return indexes, 'special'
    else:
        return indexes, 'array'


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

def _list_to_sdict(data, columns):
    if len(data) > 0 and isinstance(data[0], tuple):
        content = list(lib.to_object_array_tuples(data).T)
    elif len(data) > 0:
        # list of lists
        content = list(lib.to_object_array(data).T)
    else:
        if columns is None:
            columns = []
        return {}, columns

    if columns is None:
        columns = range(len(content))
    else:
        if len(columns) != len(content):
            raise AssertionError('%d columns passed, passed data had %s '
                                 'columns' % (len(columns), len(content)))

    sdict = dict((c, lib.maybe_convert_objects(vals))
                 for c, vals in zip(columns, content))
    return sdict, columns

def _homogenize(data, index, columns, dtype=None):
    from pandas.core.series import _sanitize_array

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

            v = _sanitize_array(v, index, dtype=dtype, copy=False,
                                raise_cast_failure=False)

        homogenized[k] = v

    return homogenized

def _put_str(s, space):
    return ('%s' % s)[:space].ljust(space)

def _is_sequence(x):
    try:
        iter(x)
        assert(not isinstance(x, basestring))
        return True
    except Exception:
        return False

def install_ipython_completers():  # pragma: no cover
    """Register the DataFrame type with IPython's tab completion machinery, so
    that it knows about accessing column names as attributes."""
    from IPython.utils.generics import complete_object

    @complete_object.when_type(DataFrame)
    def complete_dataframe(obj, prev_completions):
        return prev_completions + [c for c in obj.columns \
                    if isinstance(c, basestring) and py3compat.isidentifier(c)]

# Importing IPython brings in about 200 modules, so we want to avoid it unless
# we're in IPython (when those modules are loaded anyway).
if "IPython" in sys.modules:  # pragma: no cover
    try:
        install_ipython_completers()
    except Exception:
        pass


if __name__ == '__main__':
    import nose
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
