from __future__ import with_statement

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
                                _default_index, _stringify)
from pandas.core.generic import NDFrame
from pandas.core.index import Index, MultiIndex, _ensure_index
from pandas.core.indexing import _NDFrameIndexer, _maybe_droplevels
from pandas.core.internals import BlockManager, make_block, form_blocks
from pandas.core.series import Series, _radd_compat
from pandas.compat.scipy import scoreatpercentile as _quantile
from pandas.util import py3compat
from pandas.util.terminal import get_terminal_size
from pandas.util.decorators import deprecate, Appender, Substitution

from pandas.tseries.period import PeriodIndex

import pandas.core.algorithms as algos
import pandas.core.datetools as datetools
import pandas.core.common as com
import pandas.core.format as fmt
import pandas.core.generic as generic
import pandas.core.nanops as nanops
import pandas.lib as lib

#----------------------------------------------------------------------
# Docstring templates

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
level : int or name
    Broadcast across a level, matching Index values on the
    passed MultiIndex level

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

_merge_doc = """
Merge DataFrame objects by performing a database-style join operation by
columns or indexes.

If joining columns on columns, the DataFrame indexes *will be
ignored*. Otherwise if joining indexes on indexes or indexes on a column or
columns, the index will be passed on.

Parameters
----------%s
right : DataFrame
how : {'left', 'right', 'outer', 'inner'}, default 'inner'
    * left: use only keys from left frame (SQL: left outer join)
    * right: use only keys from right frame (SQL: right outer join)
    * outer: use union of keys from both frames (SQL: full outer join)
    * inner: use intersection of keys from both frames (SQL: inner join)
on : label or list
    Field names to join on. Must be found in both DataFrames. If on is
    None and not merging on indexes, then it merges on the intersection of
    the columns by default.
left_on : label or list, or array-like
    Field names to join on in left DataFrame. Can be a vector or list of
    vectors of the length of the DataFrame to use a particular vector as
    the join key instead of columns
right_on : label or list, or array-like
    Field names to join on in right DataFrame or vector/list of vectors per
    left_on docs
left_index : boolean, default False
    Use the index from the left DataFrame as the join key(s). If it is a
    MultiIndex, the number of keys in the other DataFrame (either the index
    or a number of columns) must match the number of levels
right_index : boolean, default False
    Use the index from the right DataFrame as the join key. Same caveats as
    left_index
sort : boolean, default True
    Sort the join keys lexicographically in the result DataFrame
suffixes : 2-length sequence (tuple, list, ...)
    Suffix to apply to overlapping column names in the left and right
    side, respectively
copy : boolean, default True
    If False, do not copy data unnecessarily

Examples
--------

>>> A              >>> B
    lkey value         rkey value
0   foo  1         0   foo  5
1   bar  2         1   bar  6
2   baz  3         2   qux  7
3   foo  4         3   bar  8

>>> merge(A, B, left_on='lkey', right_on='rkey', how='outer')
   lkey  value_x  rkey  value_y
0  bar   2        bar   6
1  bar   2        bar   8
2  baz   3        NaN   NaN
3  foo   1        foo   5
4  foo   4        foo   5
5  NaN   NaN      qux   7

Returns
-------
merged : DataFrame
"""

# Custom error class for update

class DataConflictError(Exception): pass

#----------------------------------------------------------------------
# Factory helper methods

def _arith_method(op, name, default_axis='columns'):
    def na_op(x, y):
        try:
            result = op(x, y)
        except TypeError:
            xrav = x.ravel()
            result = np.empty(x.size, dtype=x.dtype)
            if isinstance(y, np.ndarray):
                yrav = y.ravel()
                mask = notnull(xrav) & notnull(yrav)
                result[mask] = op(xrav[mask], yrav[mask])
            else:
                mask = notnull(xrav)
                result[mask] = op(xrav[mask], y)

            np.putmask(result, -mask, np.nan)
            result = result.reshape(x.shape)

        return result

    @Appender(_arith_doc % name)
    def f(self, other, axis=default_axis, level=None, fill_value=None):
        if isinstance(other, DataFrame):    # Another DataFrame
            return self._combine_frame(other, na_op, fill_value, level)
        elif isinstance(other, Series):
            return self._combine_series(other, na_op, fill_value, axis, level)
        elif isinstance(other, (list, tuple)):
            if axis is not None and self._get_axis_name(axis) == 'index':
                casted = Series(other, index=self.index)
            else:
                casted = Series(other, index=self.columns)
            return self._combine_series(casted, na_op, fill_value, axis, level)
        elif isinstance(other, np.ndarray):
            if other.ndim == 1:
                if axis is not None and self._get_axis_name(axis) == 'index':
                    casted = Series(other, index=self.index)
                else:
                    casted = Series(other, index=self.columns)
                return self._combine_series(casted, na_op, fill_value,
                                            axis, level)
            elif other.ndim == 2:
                casted = DataFrame(other, index=self.index,
                                   columns=self.columns)
                return self._combine_frame(casted, na_op, fill_value, level)
            else:  # pragma: no cover
                raise ValueError("Bad argument shape")
        else:
            return self._combine_const(other, na_op)

    f.__name__ = name

    return f

def _flex_comp_method(op, name, default_axis='columns'):

    def na_op(x, y):
        try:
            result = op(x, y)
        except TypeError:
            xrav = x.ravel()
            result = np.empty(x.size, dtype=x.dtype)
            if isinstance(y, np.ndarray):
                yrav = y.ravel()
                mask = notnull(xrav) & notnull(yrav)
                result[mask] = op(np.array(list(xrav[mask])),
                                  np.array(list(yrav[mask])))
            else:
                mask = notnull(xrav)
                result[mask] = op(np.array(list(xrav[mask])), y)

            if op == operator.ne:  # pragma: no cover
                np.putmask(result, -mask, True)
            else:
                np.putmask(result, -mask, False)
            result = result.reshape(x.shape)

        return result

    @Appender('Wrapper for flexible comparison methods %s' % name)
    def f(self, other, axis=default_axis, level=None):
        if isinstance(other, DataFrame):    # Another DataFrame
            return self._flex_compare_frame(other, na_op, level)

        elif isinstance(other, Series):
            return self._combine_series(other, na_op, None, axis, level)

        elif isinstance(other, (list, tuple)):
            if axis is not None and self._get_axis_name(axis) == 'index':
                casted = Series(other, index=self.index)
            else:
                casted = Series(other, index=self.columns)

            return self._combine_series(casted, na_op, None, axis, level)

        elif isinstance(other, np.ndarray):
            if other.ndim == 1:
                if axis is not None and self._get_axis_name(axis) == 'index':
                    casted = Series(other, index=self.index)
                else:
                    casted = Series(other, index=self.columns)

                return self._combine_series(casted, na_op, None, axis, level)

            elif other.ndim == 2:
                casted = DataFrame(other, index=self.index,
                                   columns=self.columns)

                return self._flex_compare_frame(casted, na_op, level)

            else:  # pragma: no cover
                raise ValueError("Bad argument shape")

        else:
            return self._combine_const(other, na_op)

    f.__name__ = name

    return f


def _comp_method(func, name):
    @Appender('Wrapper for comparison method %s' % name)
    def f(self, other):
        if isinstance(other, DataFrame):    # Another DataFrame
            return self._compare_frame(other, func)
        elif isinstance(other, Series):
            return self._combine_series_infer(other, func)
        else:
            return self._combine_const(other, func)

    f.__name__ = name

    return f


#----------------------------------------------------------------------
# DataFrame class


class DataFrame(NDFrame):
    _auto_consolidate = True
    _verbose_info = True
    _het_axis = 1

    _AXIS_NUMBERS = {
        'index': 0,
        'columns': 1
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
            Index to use for resulting frame. Will default to np.arange(n) if
            no indexing information part of input data and no index provided
        columns : Index or array-like
            Will default to np.arange(n) if not column labels provided
        dtype : dtype, default None
            Data type to force, otherwise infer
        copy : boolean, default False
            Copy data from inputs. Only affects DataFrame / 2d ndarray input

        Examples
        --------
        >>> d = {'col1': ts1, 'col2': ts2}
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
            if issubclass(data.dtype.type, np.datetime64):
                datacopy[mask] = lib.iNaT
            else:
                datacopy = com._maybe_upcast(datacopy)
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
            if len(data) > 0:
                if index is None and isinstance(data[0], Series):
                    index = _get_names_from_index(data)

                if isinstance(data[0], (list, tuple, dict, Series)):
                    conv_data, columns = _to_sdict(data, columns)
                    if isinstance(conv_data, dict):
                        if len(conv_data) == 0 and index is None:
                            index = np.arange(len(data))
                        mgr = self._init_dict(conv_data, index, columns,
                                              dtype=dtype)
                    else:
                        mgr = self._init_ndarray(conv_data, index, columns,
                                                 dtype=dtype, copy=copy)
                else:
                    mgr = self._init_ndarray(data, index, columns, dtype=dtype,
                                             copy=copy)
            else:
                mgr = self._init_ndarray(data, index, columns, dtype=dtype,
                                         copy=copy)
        else:
            raise PandasError('DataFrame constructor not properly called!')

        NDFrame.__init__(self, mgr)

    @classmethod
    def _from_axes(cls, data, axes):
        # for construction from BlockManager
        if isinstance(data, BlockManager):
            return cls(data)
        else:
            columns, index = axes
            return cls(data, index=index, columns=columns, copy=False)

    def _init_mgr(self, mgr, index, columns, dtype=None, copy=False):
        if columns is not None:
            mgr = mgr.reindex_axis(columns, axis=0, copy=False)
        if index is not None:
            mgr = mgr.reindex_axis(index, axis=1, copy=False)
        # do not copy BlockManager unless explicitly done
        if copy and dtype is None:
            mgr = mgr.copy()
        elif dtype is not None:
            # avoid copy if we can
            if len(mgr.blocks) > 1 or mgr.blocks[0].values.dtype != dtype:
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
            if values.dtype != dtype:
                try:
                    values = values.astype(dtype)
                except Exception:
                    raise ValueError('failed to cast to %s' % dtype)

        N, K = values.shape

        if index is None:
            index = _default_index(N)
        else:
            index = _ensure_index(index)

        if columns is None:
            columns = _default_index(K)
        else:
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

    @property
    def empty(self):
        return not (len(self.columns) > 0 and len(self.index) > 0)

    def __nonzero__(self):
        raise ValueError("Cannot call bool() on DataFrame.")

    def _need_info_repr_(self):
        """
        Check if it is needed to use info/summary view to represent a
        particular DataFrame.
        """
        config = fmt.print_config

        terminal_width, terminal_height = get_terminal_size()
        max_rows = (terminal_height if config.max_rows == 0
                    else config.max_rows)
        max_columns = config.max_columns

        if max_columns > 0:
            if len(self.index) <= max_rows and \
                    len(self.columns) <= max_columns:
                return False
            else:
                return True
        else:
            # save us
            if (len(self.index) > max_rows or
                len(self.columns) > terminal_width // 2):
                return True
            else:
                buf = StringIO()
                self.to_string(buf=buf)
                value = buf.getvalue()
                if max([len(l) for l in value.split('\n')]) > terminal_width:
                    return True
                else:
                    return False

    def __repr__(self):
        """
        Return a string representation for a particular DataFrame
        """
        buf = StringIO()
        if self._need_info_repr_():
            self.info(buf=buf, verbose=self._verbose_info)
        else:
            self.to_string(buf=buf)
        value = buf.getvalue()
        return com.console_encode(value)

    def _repr_html_(self):
        """
        Return a html representation for a particular DataFrame.
        Mainly for IPython notebook.
        """
        if fmt.print_config.notebook_repr_html:
            if self._need_info_repr_():
                return None
            else:
                return ('<div style="max-height:1000px;'
                        'max-width:1500px;overflow:auto;">\n' +
                        self.to_html() + '\n</div>')
        else:
            return None

    def __iter__(self):
        """
        Iterate over columns of the frame.
        """
        return iter(self.columns)

    def keys(self):
        return self.columns

    def iteritems(self):
        """Iterator over (column, series) pairs"""
        return ((k, self[k]) for k in self.columns)

    def iterrows(self):
        """
        Iterate over rows of DataFrame as (index, Series) pairs
        """
        columns = self.columns
        for k, v in izip(self.index, self.values):
            s = v.view(Series)
            s.index = columns
            s.name = k
            yield k, s

    def itertuples(self, index=True):
        """
        Iterate over rows of DataFrame as tuples, with index value
        as first element of the tuple
        """
        arrays = []
        if index:
            arrays.append(self.index)
        arrays.extend(self[k] for k in self.columns)
        return izip(*arrays)

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
    div = divide = _arith_method(lambda x, y: x / y, 'divide')

    radd = _arith_method(_radd_compat, 'radd')
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

    __radd__ = _arith_method(_radd_compat, '__radd__', default_axis=None)
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
        __rdiv__ = _arith_method(lambda x, y: y / x, '__rdiv__',
                                 default_axis=None)

    def __neg__(self):
        arr = operator.neg(self.values)
        return self._wrap_array(arr, self.axes, copy=False)

    # Comparison methods
    __eq__ = _comp_method(operator.eq, '__eq__')
    __ne__ = _comp_method(operator.ne, '__ne__')
    __lt__ = _comp_method(operator.lt, '__lt__')
    __gt__ = _comp_method(operator.gt, '__gt__')
    __le__ = _comp_method(operator.le, '__le__')
    __ge__ = _comp_method(operator.ge, '__ge__')

    eq = _flex_comp_method(operator.eq, 'eq')
    ne = _flex_comp_method(operator.ne, 'ne')
    gt = _flex_comp_method(operator.gt, 'gt')
    lt = _flex_comp_method(operator.lt, 'lt')
    ge = _flex_comp_method(operator.ge, 'ge')
    le = _flex_comp_method(operator.le, 'le')

    def dot(self, other):
        """
        Matrix multiplication with DataFrame or Series objects

        Parameters
        ----------
        other : DataFrame or Series

        Returns
        -------
        dot_product : DataFrame or Series
        """
        common = self.columns.union(other.index)
        if len(common) > len(self.columns) or len(common) > len(other.index):
            raise ValueError('matrices are not aligned')
        left = self.reindex(columns=common, copy=False)
        right = other.reindex(index=common, copy=False)
        lvals = left.values
        rvals = right.values
        if isinstance(other, DataFrame):
            return DataFrame(np.dot(lvals, rvals), index=self.index, columns=other.columns)
        elif isinstance(other, Series):
            return Series(np.dot(lvals, rvals), index=left.index)
        else:
            raise TypeError('unsupported type: %s' % type(other))

    #----------------------------------------------------------------------
    # IO methods (to / from other formats)

    @classmethod
    def from_dict(cls, data, orient='columns', dtype=None):
        """
        Construct DataFrame from dict of array-like or dicts

        Parameters
        ----------
        data : dict
            {field : array-like} or {field : dict}
        orient : {'columns', 'index'}, default 'columns'
            The "orientation" of the data. If the keys of the passed dict
            should be the columns of the resulting DataFrame, pass 'columns'
            (default). Otherwise if the keys should be rows, pass 'index'.

        Returns
        -------
        DataFrame
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

    def to_dict(self, outtype='dict'):
        """
        Convert DataFrame to dictionary.

        Parameters
        ----------
        outtype : str {'dict', 'list', 'series'}
            Determines the type of the values of the dictionary. The
            default `dict` is a nested dictionary {column -> {index -> value}}.
            `list` returns {column -> list(values)}. `series` returns
            {column -> Series(values)}.
            Abbreviations are allowed.


        Returns
        -------
        result : dict like {column -> {index -> value}}
        """
        if outtype.lower().startswith('d'):
            return dict((k, v.to_dict()) for k, v in self.iteritems())
        elif outtype.lower().startswith('l'):
            return dict((k, v.tolist()) for k, v in self.iteritems())
        elif outtype.lower().startswith('s'):
            return dict((k, v) for k,v in self.iteritems())
        else: # pragma: no cover
            raise ValueError("outtype %s not understood" % outtype)

    @classmethod
    def from_records(cls, data, index=None, exclude=None, columns=None,
                     names=None, coerce_float=False):
        """
        Convert structured or record ndarray to DataFrame

        Parameters
        ----------
        data : ndarray (structured dtype), list of tuples, or DataFrame
        index : string, list of fields, array-like
            Field of array to use as the index, alternately a specific set of
            input labels to use
        exclude: sequence, default None
            Columns or fields to exclude
        columns : sequence, default None
            Column names to use, replacing any found in passed data
        coerce_float : boolean, default False
            Attempt to convert values to non-string, non-numeric objects (like
            decimal.Decimal) to floating point, useful for SQL result sets

        Returns
        -------
        df : DataFrame
        """
        import warnings

        # Make a copy of the input columns so we can modify it
        if columns is not None:
            columns = list(columns)

            if len(algos.unique(columns)) < len(columns):
                raise ValueError('Non-unique columns not yet supported in from_records')

        if names is not None:  # pragma: no cover
            columns = names
            warnings.warn("'names' parameter to DataFrame.from_records is "
                          "being renamed to 'columns', 'names' will be "
                          "removed in 0.8.0",
                          FutureWarning)

        if isinstance(data, (np.ndarray, DataFrame, dict)):
            columns, sdict = _rec_to_dict(data)
        else:
            sdict, columns = _to_sdict(data, columns,
                                       coerce_float=coerce_float)

        if exclude is None:
            exclude = set()
        else:
            exclude = set(exclude)

        for col in exclude:
            del sdict[col]
            columns.remove(col)

        if index is not None:
            if (isinstance(index, basestring) or
                not hasattr(index, "__iter__")):
                result_index = sdict.pop(index)
                result_index = Index(result_index, name=index)
                columns.remove(index)
            else:
                try:
                    arrays = []
                    for field in index:
                        arrays.append(sdict[field])
                    for field in index:
                        del sdict[field]
                        columns.remove(field)
                    result_index = MultiIndex.from_arrays(arrays, names=index)
                except Exception:
                    result_index = index
        elif isinstance(data, dict) and len(data) > 0:
            # utilize first element of sdict to get length
            result_index = np.arange(len(data.values()[0]))
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
            arrays = [self.index.values] + [self[c].values
                                            for c in self.columns]
            names = ['index'] + list(map(str, self.columns))
        else:
            arrays = [self[c].values for c in self.columns]
            names = list(map(str, self.columns))

        dtype = np.dtype([(x, v.dtype) for x, v in zip(names, arrays)])
        return np.rec.fromarrays(arrays, dtype=dtype, names=names)

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
                 parse_dates=True, encoding=None):
        """
        Read delimited file into DataFrame

        Parameters
        ----------
        path : string file path or file handle / StringIO
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
                          parse_dates=parse_dates, index_col=index_col,
                          encoding=encoding)

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

    def to_panel(self):
        """
        Transform long (stacked) format (DataFrame) into wide (3D, Panel)
        format.

        Currently the index of the DataFrame must be a 2-level MultiIndex. This
        may be generalized later

        Returns
        -------
        panel : Panel
        """
        from pandas.core.panel import Panel
        from pandas.core.reshape import block2d_to_block3d

        # only support this kind for now
        assert(isinstance(self.index, MultiIndex) and
               len(self.index.levels) == 2)

        self._consolidate_inplace()

        # minor axis must be sorted
        if self.index.lexsort_depth < 2:
            selfsorted = self.sortlevel(0)
        else:
            selfsorted = self

        major_axis, minor_axis = selfsorted.index.levels
        major_labels, minor_labels = selfsorted.index.labels

        shape = len(major_axis), len(minor_axis)

        new_blocks = []
        for block in selfsorted._data.blocks:
            newb = block2d_to_block3d(block.values.T, block.items, shape,
                                      major_labels, minor_labels,
                                      ref_items=selfsorted.columns)
            new_blocks.append(newb)

        new_axes = [selfsorted.columns, major_axis, minor_axis]
        new_mgr = BlockManager(new_blocks, new_axes)

        return Panel(new_mgr)

    to_wide = deprecate('to_wide', to_panel)

    def _helper_csvexcel(self, writer, na_rep=None, cols=None,
                         header=True, index=True,
                         index_label=None, float_format=None):
        if cols is None:
            cols = self.columns

        series = {}
        for k, v in self._series.iteritems():
            series[k] = v.values

        has_aliases = isinstance(header, (tuple, list, np.ndarray))
        if has_aliases or header:
            if index:
                # should write something for index label
                if index_label is not False:
                    if index_label is None:
                        if isinstance(self.index, MultiIndex):
                            index_label = []
                            for i, name in enumerate(self.index.names):
                                if name is None:
                                    name = ''
                                index_label.append(name)
                        else:
                            index_label = self.index.name
                            if index_label is None:
                                index_label = ['']
                            else:
                                index_label = [index_label]
                    elif not isinstance(index_label, (list, tuple, np.ndarray)):
                        # given a string for a DF with Index
                        index_label = [index_label]

                    encoded_labels = list(index_label)
                else:
                    encoded_labels = []

                if has_aliases:
                    if len(header) != len(cols):
                        raise ValueError(('Writing %d cols but got %d aliases'
                                          % (len(cols), len(header))))
                    else:
                        write_cols = header
                else:
                    write_cols = cols
                encoded_cols = list(write_cols)

                writer.writerow(encoded_labels + encoded_cols)
            else:
                encoded_cols = list(cols)
                writer.writerow(encoded_cols)

        nlevels = getattr(self.index, 'nlevels', 1)
        for j, idx in enumerate(self.index):
            row_fields = []
            if index:
                if nlevels == 1:
                    row_fields = [idx]
                else:  # handle MultiIndex
                    row_fields = list(idx)
            for i, col in enumerate(cols):
                val = series[col][j]
                if lib.checknull(val):
                    val = na_rep

                if float_format is not None and com.is_float(val):
                    val = float_format % val
                elif isinstance(val, np.datetime64):
                    val = lib.Timestamp(val)._repr_base

                row_fields.append(val)

            writer.writerow(row_fields)

    def to_csv(self, path_or_buf, sep=",", na_rep='', float_format=None,
               cols=None, header=True, index=True, index_label=None,
               mode='w', nanRep=None, encoding=None, quoting=None):
        """
        Write DataFrame to a comma-separated values (csv) file

        Parameters
        ----------
        path_or_buf : string or file handle / StringIO
            File path
        na_rep : string, default ''
            Missing data representation
        float_format : string, default None
            Format string for floating point numbers
        cols : sequence, optional
            Columns to write
        header : boolean or list of string, default True
            Write out column names. If a list of string is given it is
            assumed to be aliases for the column names
        index : boolean, default True
            Write row names (index)
        index_label : string or sequence, or False, default None
            Column label for index column(s) if desired. If None is given, and
            `header` and `index` are True, then the index names are used. A
            sequence should be given if the DataFrame uses MultiIndex.  If
            False do not print fields for index names. Use index_label=False
            for easier importing in R
        mode : Python write mode, default 'w'
        sep : character, default ","
            Field delimiter for the output file.
        encoding : string, optional
            a string representing the encoding to use if the contents are
            non-ascii, for python versions prior to 3
        """
        if nanRep is not None:  # pragma: no cover
            import warnings
            warnings.warn("nanRep is deprecated, use na_rep",
                          FutureWarning)
            na_rep = nanRep

        if hasattr(path_or_buf, 'read'):
            f = path_or_buf
            close = False
        else:
            f = com._get_handle(path_or_buf, mode, encoding=encoding)
            close = True


        if quoting is None:
            quoting = csv.QUOTE_MINIMAL

        try:
            if encoding is not None:
                csvout = com.UnicodeWriter(f, lineterminator='\n',
                                           delimiter=sep, encoding=encoding,
                                           quoting=quoting)
            else:
                csvout = csv.writer(f, lineterminator='\n', delimiter=sep,
                                    quoting=quoting)
            self._helper_csvexcel(csvout, na_rep=na_rep,
                                  float_format=float_format, cols=cols,
                                  header=header, index=index,
                                  index_label=index_label)

        finally:
            if close:
                f.close()

    def to_excel(self, excel_writer, sheet_name='sheet1', na_rep='',
                 float_format=None, cols=None, header=True, index=True,
                 index_label=None):
        """
        Write DataFrame to a excel sheet

        Parameters
        ----------
        excel_writer : string or ExcelWriter object
            File path or existing ExcelWriter
        sheet_name : string, default 'sheet1'
            Name of sheet which will contain DataFrame
        na_rep : string, default ''
            Missing data representation
        float_format : string, default None
            Format string for floating point numbers
        cols : sequence, optional
            Columns to write
        header : boolean or list of string, default True
            Write out column names. If a list of string is given it is
            assumed to be aliases for the column names
        index : boolean, default True
            Write row names (index)
        index_label : string or sequence, default None
            Column label for index column(s) if desired. If None is given, and
            `header` and `index` are True, then the index names are used. A
            sequence should be given if the DataFrame uses MultiIndex.

        Notes
        -----
        If passing an existing ExcelWriter object, then the sheet will be added
        to the existing workbook.  This can be used to save different
        DataFrames to one workbook
        >>> writer = ExcelWriter('output.xlsx')
        >>> df1.to_excel(writer,'sheet1')
        >>> df2.to_excel(writer,'sheet2')
        >>> writer.save()
        """
        from pandas.io.parsers import ExcelWriter
        need_save = False
        if isinstance(excel_writer, basestring):
            excel_writer = ExcelWriter(excel_writer)
            need_save = True
        excel_writer.cur_sheet = sheet_name
        self._helper_csvexcel(excel_writer, na_rep=na_rep,
                              float_format=float_format, cols=cols,
                              header=header, index=index,
                              index_label=index_label)
        if need_save:
            excel_writer.save()

    @Appender(fmt.docstring_to_string, indents=1)
    def to_string(self, buf=None, columns=None, col_space=None, colSpace=None,
                  header=True, index=True, na_rep='NaN', formatters=None,
                  float_format=None, sparsify=None, nanRep=None,
                  index_names=True, justify=None, force_unicode=False):
        """
        Render a DataFrame to a console-friendly tabular output.
        """

        if nanRep is not None:  # pragma: no cover
            import warnings
            warnings.warn("nanRep is deprecated, use na_rep",
                          FutureWarning)
            na_rep = nanRep

        if colSpace is not None:  # pragma: no cover
            import warnings
            warnings.warn("colSpace is deprecated, use col_space",
                          FutureWarning)
            col_space = colSpace

        formatter = fmt.DataFrameFormatter(self, buf=buf, columns=columns,
                                           col_space=col_space, na_rep=na_rep,
                                           formatters=formatters,
                                           float_format=float_format,
                                           sparsify=sparsify,
                                           justify=justify,
                                           index_names=index_names,
                                           header=header, index=index)
        formatter.to_string(force_unicode=force_unicode)

        if buf is None:
            result = formatter.buf.getvalue()
            if not force_unicode:
                try:
                    result = str(result)
                except ValueError:
                    pass
            return result

    @Appender(fmt.docstring_to_string, indents=1)
    def to_html(self, buf=None, columns=None, col_space=None, colSpace=None,
                header=True, index=True, na_rep='NaN', formatters=None,
                float_format=None, sparsify=None, index_names=True,
                justify=None, force_unicode=False, bold_rows=True,
                classes=None):
        """
        to_html-specific options
        bold_rows : boolean, default True
            Make the row labels bold in the output

        Render a DataFrame to an html table.
        """

        if colSpace is not None:  # pragma: no cover
            import warnings
            warnings.warn("colSpace is deprecated, use col_space",
                          FutureWarning)
            col_space = colSpace

        formatter = fmt.DataFrameFormatter(self, buf=buf, columns=columns,
                                           col_space=col_space, na_rep=na_rep,
                                           formatters=formatters,
                                           float_format=float_format,
                                           sparsify=sparsify,
                                           justify=justify,
                                           index_names=index_names,
                                           header=header, index=index,
                                           bold_rows=bold_rows)
        formatter.to_html(classes=classes)

        if buf is None:
            return formatter.buf.getvalue()

    @Appender(fmt.docstring_to_string, indents=1)
    def to_latex(self, buf=None, columns=None, col_space=None, colSpace=None,
                 header=True, index=True, na_rep='NaN', formatters=None,
                 float_format=None, sparsify=None, index_names=True,
                 bold_rows=True):
        """
        to_latex-specific options
        bold_rows : boolean, default True
            Make the row labels bold in the output

        Render a DataFrame to a tabular environment table.
        You can splice this into a LaTeX document.
        """
        formatter = fmt.DataFrameFormatter(self, buf=buf, columns=columns,
                                           col_space=col_space, na_rep=na_rep,
                                           header=header, index=index,
                                           formatters=formatters,
                                           float_format=float_format,
                                           bold_rows=bold_rows,
                                           sparsify=sparsify,
                                           index_names=index_names)
        formatter.to_latex()

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
        from pandas.core.format import _put_lines

        if buf is None:  # pragma: no cover
            buf = sys.stdout

        lines = []

        lines.append(str(type(self)))
        lines.append(self.index.summary())

        if len(self.columns) == 0:
            lines.append('Empty %s' % type(self).__name__)
            _put_lines(buf, lines)
            return

        cols = self.columns

        # hack
        if verbose and len(self.columns) < 100:
            lines.append('Data columns:')
            space = max([len(_stringify(k)) for k in self.columns]) + 4
            counts = self.count()
            assert(len(cols) == len(counts))
            for col, count in counts.iteritems():
                if not isinstance(col, basestring):
                    col = str(col)
                lines.append(_put_str(col, space) +
                             '%d  non-null values' % count)
        else:
            lines.append(self.columns.summary(name='Columns'))

        counts = self.get_dtype_counts()
        dtypes = ['%s(%d)' % k for k in sorted(counts.iteritems())]
        lines.append('dtypes: %s' % ', '.join(dtypes))
        _put_lines(buf, lines)

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
        for i in range(len(self.columns)):
            series = self.icol(i)
            # endianness can cause dtypes to look different
            dtype_str = str(series.dtype)
            if dtype_str in counts:
                counts[dtype_str] += 1
            else:
                counts[dtype_str] = 1
        return Series(counts)

    #----------------------------------------------------------------------
    # properties for index and columns

    columns = lib.AxisProperty(0)
    index = lib.AxisProperty(1)

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

    def swapaxes(self, i, j):
        """
        Like ndarray.swapaxes, equivalent to transpose

        Returns
        -------
        swapped : DataFrame
            View on original data (no copy)
        """
        if i in (0, 1) and j in (0, 1):
            if i == j:
                return self
            return self._constructor(data=self.values.T, index=self.columns,
                                     columns=self.index, copy=False)
        else:
            raise ValueError('Axis numbers must be in (0, 1)')

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

            made_bigger = not np.array_equal(new_columns, self.columns)

            # how to make this logic simpler?
            if made_bigger:
                com._possibly_cast_item(result, col, likely_dtype)

            return result.set_value(index, col, value)

    def irow(self, i, copy=False):
        """
        Retrieve the i-th row or rows of the DataFrame by location

        Parameters
        ----------
        i : int, slice, or sequence of integers

        Notes
        -----
        If slice passed, the resulting data will be a view

        Returns
        -------
        row : Series (int) or DataFrame (slice, sequence)
        """
        if isinstance(i, slice):
            return self[i]
        else:
            label = self.index[i]
            if isinstance(label, Index):
                return self.reindex(label)
            else:
                try:
                    new_values = self._data.fast_2d_xs(i, copy=copy)
                except:
                    new_values = self._data.fast_2d_xs(i, copy=True)
                return Series(new_values, index=self.columns,
                              name=self.index[i])

    def icol(self, i):
        """
        Retrieve the i-th column or columns of the DataFrame by location

        Parameters
        ----------
        i : int, slice, or sequence of integers

        Notes
        -----
        If slice passed, the resulting data will be a view

        Returns
        -------
        column : Series (int) or DataFrame (slice, sequence)
        """
        label = self.columns[i]
        if isinstance(i, slice):
            # need to return view
            lab_slice = slice(label[0], label[-1])
            return self.ix[:, lab_slice]
        else:
            label = self.columns[i]
            if isinstance(label, Index):
                return self.reindex(columns=label)

            values = self._data.iget(i)
            return Series(values, index=self.index, name=label)

    def _ixs(self, i, axis=0):
        if axis == 0:
            return self.irow(i)
        else:
            return self.icol(i)

    def iget_value(self, i, j):
        """
        Return scalar value stored at row i and column j, where i and j are
        integers

        Parameters
        ----------
        i : int
        j : int

        Returns
        -------
        value : scalar value
        """
        row = self.index[i]
        col = self.columns[j]
        return self.get_value(row, col)

    def __getitem__(self, key):
        # slice rows
        if isinstance(key, slice):
            from pandas.core.indexing import _is_index_slice
            idx_type = self.index.inferred_type
            if idx_type == 'floating':
                indexer = self.ix._convert_to_indexer(key, axis=0)
            elif idx_type == 'integer' or _is_index_slice(key):
                indexer = key
            else:
                indexer = self.ix._convert_to_indexer(key, axis=0)
            new_data = self._data.get_slice(indexer, axis=1)
            return self._constructor(new_data)
        # either boolean or fancy integer index
        elif isinstance(key, (np.ndarray, list)):
            if isinstance(key, list):
                key = lib.list_to_object_array(key)

            # also raises Exception if object array with NA values
            if com._is_bool_indexer(key):
                key = np.asarray(key, dtype=bool)
            return self._getitem_array(key)
        elif isinstance(self.columns, MultiIndex):
            return self._getitem_multilevel(key)
        elif isinstance(key, DataFrame):
            values = key.values
            if values.dtype == bool:
                return self.values[values]
            else:
                raise ValueError('Cannot index using non-boolean DataFrame')
        else:
            return self._get_item_cache(key)

    def _getitem_array(self, key):
        if key.dtype == np.bool_:
            if len(key) != len(self.index):
                raise ValueError('Item wrong length %d instead of %d!' %
                                 (len(key), len(self.index)))

            inds, = key.nonzero()
            return self.take(inds)
        else:
            if self.columns.is_unique:
                indexer = self.columns.get_indexer(key)
                mask = indexer == -1
                if mask.any():
                    raise KeyError("No column(s) named: %s" % str(key[mask]))
                result = self.reindex(columns=key)
                if result.columns.name is None:
                    result.columns.name = self.columns.name
                return result
            else:
                mask = self.columns.isin(key)
                return self.take(mask.nonzero()[0], axis=1)

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
            if len(result.columns) == 1:
                top = result.columns[0]
                if (type(top) == str and top == '' or
                        type(top) == tuple and top[0] == ''):
                    result = Series(result[''], index=self.index, name=key)
            return result
        else:
            return self._get_item_cache(key)

    def _box_item_values(self, key, values):
        items = self.columns[self.columns.get_loc(key)]
        if values.ndim == 2:
            return DataFrame(values.T, columns=items, index=self.index)
        else:
            return Series(values, index=self.index, name=items)

    def __getattr__(self, name):
        """After regular attribute access, try looking up the name of a column.
        This allows simpler access to columns for interactive use."""
        if name in self.columns:
            return self[name]
        raise AttributeError("'%s' object has no attribute '%s'" %
                             (type(self).__name__, name))

    def __setattr__(self, name, value):
        """After regular attribute access, try looking up the name of a column.
        This allows simpler access to columns for interactive use."""
        if name == '_data':
            super(DataFrame, self).__setattr__(name, value)
        else:
            try:
                existing = getattr(self, name)
                if isinstance(existing, Index):
                    super(DataFrame, self).__setattr__(name, value)
                elif name in self.columns:
                    self[name] = value
                else:
                    object.__setattr__(self, name, value)
            except (AttributeError, TypeError):
                object.__setattr__(self, name, value)

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
            if isinstance(keys, np.ndarray) and keys.dtype == np.bool_:
                # boolean slicing should happen on rows, consistent with
                # behavior of getitem
                self.ix[keys, :] = value
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
        value = self._sanitize_column(key, value)
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
        value = self._sanitize_column(column, value)
        self._data.insert(loc, column, value)

    def _sanitize_column(self, key, value):
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
                    value = com._asarray_tuplesafe(value)
                else:
                    value = value.copy()
        else:
            value = np.repeat(value, len(self.index))
            if key in self.columns:
                existing_column = self[key]
                # special case for now
                if (com.is_float_dtype(existing_column) and
                    com.is_integer_dtype(value)):
                    value = value.astype(np.float64)

        return np.atleast_2d(np.asarray(value))

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

    def xs(self, key, axis=0, level=None, copy=True):
        """
        Returns a cross-section (row(s) or column(s)) from the DataFrame.
        Defaults to cross-section on the rows (axis=0).

        Parameters
        ----------
        key : object
            Some label contained in the index, or partially in a MultiIndex
        axis : int, default 0
            Axis to retrieve cross-section on
        level : object, defaults to first n levels (n=1 or len(key))
            In case of a key partially contained in a MultiIndex, indicate
            which levels are used. Levels can be referred by label or position.
        copy : boolean, default True
            Whether to make a copy of the data

        Examples
        --------
        >>> df
           A  B  C
        a  4  5  2
        b  4  0  9
        c  9  7  3
        >>> df.xs('a')
        A    4
        B    5
        C    2
        Name: a
        >>> df.xs('C', axis=1)
        a    2
        b    9
        c    3
        Name: C
        >>> s = df.xs('a', copy=False)
        >>> s['A'] = 100
        >>> df
             A  B  C
        a  100  5  2
        b    4  0  9
        c    9  7  3


        >>> df
                            A  B  C  D
        first second third
        bar   one    1      4  1  8  9
              two    1      7  5  5  0
        baz   one    1      6  6  8  0
              three  2      5  3  5  3
        >>> df.xs(('baz', 'three'))
               A  B  C  D
        third
        2      5  3  5  3
        >>> df.xs('one', level=1)
                     A  B  C  D
        first third
        bar   1      4  1  8  9
        baz   1      6  6  8  0
        >>> df.xs(('baz', 2), level=[0, 'third'])
                A  B  C  D
        second
        three   5  3  5  3

        Returns
        -------
        xs : Series or DataFrame
        """
        labels = self._get_axis(axis)
        if level is not None:
            loc, new_ax = labels.get_loc_level(key, level=level)

            # level = 0
            if not isinstance(loc, slice):
                indexer = [slice(None, None)] * 2
                indexer[axis] = loc
                indexer = tuple(indexer)
            else:
                indexer = loc

            result = self.ix[indexer]
            setattr(result, result._get_axis_name(axis), new_ax)
            return result

        if axis == 1:
            data = self[key]
            if copy:
                data = data.copy()
            return data

        self._consolidate_inplace()

        index = self.index
        if isinstance(index, MultiIndex):
            loc, new_index = self.index.get_loc_level(key)
        else:
            loc = self.index.get_loc(key)

            if isinstance(loc, np.ndarray):
                if loc.dtype == np.bool_:
                    inds, = loc.nonzero()
                    if len(inds) == 1:
                        loc = inds[0]
                    else:
                        return self.take(inds, axis=axis)
                else:
                    return self.take(loc, axis=axis)

            if not np.isscalar(loc):
                new_index = self.index[loc]

        if np.isscalar(loc):
            new_values = self._data.fast_2d_xs(loc, copy=copy)
            return Series(new_values, index=self.columns,
                          name=self.index[loc])
        else: # isinstance(loc, slice) or loc.dtype == np.bool_:
            result = self[loc]
            result.index = new_index
            return result
        # else:
        #     return self.take(loc)

    def lookup(self, row_labels, col_labels):
        """
        Label-based "fancy indexing" function for DataFrame. Given equal-length
        arrays of row and column labels, return an array of the values
        corresponding to each (row, col)  pair.

        Parameters
        ----------
        row_labels : sequence
        col_labels : sequence

        Notes
        -----
        Akin to

        result = []
        for row, col in zip(row_labels, col_labels):
            result.append(df.get_value(row, col))

        Example
        -------
        values : ndarray
        """
        from itertools import izip

        n = len(row_labels)
        assert(n == len(col_labels))

        thresh = 1000
        if not self._is_mixed_type or n > thresh:
            values = self.values
            ridx = self.index.get_indexer(row_labels)
            cidx = self.columns.get_indexer(col_labels)
            if (ridx == -1).any():
                raise ValueError('One or more row labels was not found')
            if (cidx == -1).any():
                raise ValueError('One or more column labels was not found')
            flat_index = ridx * len(self.columns) + cidx
            result = values.flat[flat_index]
        else:
            result = np.empty(n, dtype='O')
            for i, (r, c) in enumerate(izip(row_labels, col_labels)):
                result[i] = self.get_value(r, c)

        if result.dtype == 'O':
            result = lib.maybe_convert_objects(result)

        return result

    #----------------------------------------------------------------------
    # Reindexing and alignment

    def align(self, other, join='outer', axis=None, level=None, copy=True,
              fill_value=np.nan, method=None, limit=None, fill_axis=0):
        """
        Align two DataFrame object on their index and columns with the
        specified join method for each axis Index

        Parameters
        ----------
        other : DataFrame or Series
        join : {'outer', 'inner', 'left', 'right'}, default 'outer'
        axis : {0, 1, None}, default None
            Align on index (0), columns (1), or both (None)
        level : int or name
            Broadcast across a level, matching Index values on the
            passed MultiIndex level
        copy : boolean, default True
            Always returns new objects. If copy=False and no reindexing is
            required then original objects are returned.
        fill_value : scalar, default np.NaN
            Value to use for missing values. Defaults to NaN, but can be any
            "compatible" value
        method : str, default None
        limit : int, default None
        fill_axis : {0, 1}, default 0
            Filling axis, method and limit

        Returns
        -------
        (left, right) : (DataFrame, type of other)
            Aligned objects
        """
        if isinstance(other, DataFrame):
            return self._align_frame(other, join=join, axis=axis, level=level,
                                     copy=copy, fill_value=fill_value,
                                     method=method, limit=limit,
                                     fill_axis=fill_axis)
        elif isinstance(other, Series):
            return self._align_series(other, join=join, axis=axis, level=level,
                                      copy=copy, fill_value=fill_value,
                                      method=method, limit=limit,
                                      fill_axis=fill_axis)
        else:  # pragma: no cover
            raise TypeError('unsupported type: %s' % type(other))

    def _align_frame(self, other, join='outer', axis=None, level=None,
                     copy=True, fill_value=np.nan, method=None, limit=None,
                     fill_axis=0):
        # defaults
        join_index, join_columns = None, None
        ilidx, iridx = None, None
        clidx, cridx = None, None

        if axis is None or axis == 0:
            if not self.index.equals(other.index):
                join_index, ilidx, iridx = \
                    self.index.join(other.index, how=join, level=level,
                                    return_indexers=True)

        if axis is None or axis == 1:
            if not self.columns.equals(other.columns):
                join_columns, clidx, cridx = \
                    self.columns.join(other.columns, how=join, level=level,
                                      return_indexers=True)

        left = self._reindex_with_indexers(join_index, ilidx,
                                           join_columns, clidx, copy,
                                           fill_value=fill_value)
        right = other._reindex_with_indexers(join_index, iridx,
                                             join_columns, cridx, copy,
                                             fill_value=fill_value)

        if method is not None:
            left = left.fillna(axis=fill_axis, method=method, limit=limit)
            right = right.fillna(axis=fill_axis, method=method, limit=limit)

        return left, right

    def _align_series(self, other, join='outer', axis=None, level=None,
                      copy=True, fill_value=None, method=None, limit=None,
                      fill_axis=0):
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
                join_index, lidx, ridx = \
                    self.columns.join(other.index, how=join,
                                      return_indexers=True)

            if lidx is not None:
                fdata = fdata.reindex_indexer(join_index, lidx, axis=0)
        else:
            raise ValueError('Must specify axis=0 or 1')

        if copy and fdata is self._data:
            fdata = fdata.copy()

        left_result = DataFrame(fdata)
        right_result = other if ridx is None else other.reindex(join_index)

        fill_na = notnull(fill_value) or (method is not None)
        if fill_na:
            return (left_result.fillna(fill_value, method=method, limit=limit,
                                       axis=fill_axis),
                    right_result.fillna(fill_value, method=method, limit=limit))
        else:
            return left_result, right_result

    def reindex(self, index=None, columns=None, method=None, level=None,
                fill_value=np.nan, limit=None, copy=True):
        """Conform DataFrame to new index with optional filling logic, placing
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
        level : int or name
            Broadcast across a level, matching Index values on the
            passed MultiIndex level
        fill_value : scalar, default np.NaN
            Value to use for missing values. Defaults to NaN, but can be any
            "compatible" value
        limit : int, default None
            Maximum size gap to forward or backward fill

        Examples
        --------
        >>> df.reindex(index=[date1, date2, date3], columns=['A', 'B', 'C'])

        Returns
        -------
        reindexed : same type as calling instance
        """
        self._consolidate_inplace()
        frame = self

        if (index is not None and columns is not None
            and method is None and level is None
            and not self._is_mixed_type):
            return self._reindex_multi(index, columns, copy, fill_value)

        if columns is not None:
            frame = frame._reindex_columns(columns, copy, level,
                                           fill_value, limit)

        if index is not None:
            frame = frame._reindex_index(index, method, copy, level,
                                         fill_value, limit)

        return frame

    def reindex_axis(self, labels, axis=0, method=None, level=None, copy=True,
                     limit=None, fill_value=np.nan):
        """Conform DataFrame to new index with optional filling logic, placing
        NA/NaN in locations having no value in the previous index. A new object
        is produced unless the new index is equivalent to the current one and
        copy=False

        Parameters
        ----------
        index : array-like, optional
            New labels / index to conform to. Preferably an Index object to
            avoid duplicating data
        axis : {0, 1}
            0 -> index (rows)
            1 -> columns
        method : {'backfill', 'bfill', 'pad', 'ffill', None}, default None
            Method to use for filling holes in reindexed DataFrame
            pad / ffill: propagate last valid observation forward to next valid
            backfill / bfill: use NEXT valid observation to fill gap
        copy : boolean, default True
            Return a new object, even if the passed indexes are the same
        level : int or name
            Broadcast across a level, matching Index values on the
            passed MultiIndex level
        limit : int, default None
            Maximum size gap to forward or backward fill

        Examples
        --------
        >>> df.reindex_axis(['A', 'B', 'C'], axis=1)

        See also
        --------
        DataFrame.reindex, DataFrame.reindex_like

        Returns
        -------
        reindexed : same type as calling instance
        """
        self._consolidate_inplace()
        if axis == 0:
            return self._reindex_index(labels, method, copy, level,
                                       fill_value=fill_value,
                                       limit=limit)
        elif axis == 1:
            return self._reindex_columns(labels, copy, level,
                                         fill_value=fill_value,
                                         limit=limit)
        else:  # pragma: no cover
            raise ValueError('Must specify axis=0 or 1')

    def _reindex_multi(self, new_index, new_columns, copy, fill_value):
        new_index, row_indexer = self.index.reindex(new_index)
        new_columns, col_indexer = self.columns.reindex(new_columns)

        if row_indexer is not None and col_indexer is not None:
            new_values = com.take_2d_multi(self.values, row_indexer,
                                           col_indexer, fill_value=fill_value)
            return DataFrame(new_values, index=new_index, columns=new_columns)
        elif row_indexer is not None:
            return self._reindex_with_indexers(new_index, row_indexer,
                                               None, None, copy, fill_value)
        elif col_indexer is not None:
            return self._reindex_with_indexers(None, None,
                                               new_columns, col_indexer,
                                               copy, fill_value)
        else:
            return self.copy() if copy else self

    def _reindex_index(self, new_index, method, copy, level, fill_value=np.nan,
                       limit=None):
        new_index, indexer = self.index.reindex(new_index, method, level,
                                                limit=limit)
        return self._reindex_with_indexers(new_index, indexer, None, None,
                                           copy, fill_value)

    def _reindex_columns(self, new_columns, copy, level, fill_value=np.nan,
                         limit=None):
        new_columns, indexer = self.columns.reindex(new_columns, level=level,
                                                    limit=limit)
        return self._reindex_with_indexers(None, None, new_columns, indexer,
                                           copy, fill_value)

    def _reindex_with_indexers(self, index, row_indexer, columns, col_indexer,
                               copy, fill_value):
        new_data = self._data
        if row_indexer is not None:
            row_indexer = com._ensure_int64(row_indexer)
            new_data = new_data.reindex_indexer(index, row_indexer, axis=1,
                                                fill_value=fill_value)
        elif index is not None and index is not new_data.axes[1]:
            new_data = new_data.copy(deep=copy)
            new_data.axes[1] = index

        if col_indexer is not None:
            # TODO: speed up on homogeneous DataFrame objects
            col_indexer = com._ensure_int64(col_indexer)
            new_data = new_data.reindex_indexer(columns, col_indexer, axis=0,
                                                fill_value=fill_value)
        elif columns is not None and columns is not new_data.axes[0]:
            new_data = new_data.reindex_items(columns, copy=copy,
                                              fill_value=fill_value)

        if copy and new_data is self._data:
            new_data = new_data.copy()

        return DataFrame(new_data)

    def reindex_like(self, other, method=None, copy=True, limit=None):
        """
        Reindex DataFrame to match indices of another DataFrame, optionally
        with filling logic

        Parameters
        ----------
        other : DataFrame
        method : string or None
        copy : boolean, default True
        limit : int, default None
            Maximum size gap to forward or backward fill

        Notes
        -----
        Like calling s.reindex(index=other.index, columns=other.columns,
                               method=...)

        Returns
        -------
        reindexed : DataFrame
        """
        return self.reindex(index=other.index, columns=other.columns,
                            method=method, copy=copy, limit=limit)

    truncate = generic.truncate

    def set_index(self, keys, drop=True, append=False, inplace=False,
                  verify_integrity=False):
        """
        Set the DataFrame index (row labels) using one or more existing
        columns. By default yields a new object.

        Parameters
        ----------
        keys : column label or list of column labels / arrays
        drop : boolean, default True
            Delete columns to be used as the new index
        append : boolean, default False
            Whether to append columns to existing index
        inplace : boolean, default False
            Modify the DataFrame in place (do not create a new object)
        verify_integrity : boolean, default False
            Check the new index for duplicates. Otherwise defer the check until
            necessary. Setting to False will improve the performance of this
            method

        Examples
        --------
        indexed_df = df.set_index(['A', 'B'])
        indexed_df2 = df.set_index(['A', [0, 1, 2, 0, 1, 2]])
        indexed_df3 = df.set_index([[0, 1, 2, 0, 1, 2]])

        Returns
        -------
        dataframe : DataFrame
        """
        if not isinstance(keys, list):
            keys = [keys]

        if inplace:
            frame = self

        else:
            frame = self.copy()

        arrays = []
        names = []
        if append:
            names = [x for x in self.index.names]
            if isinstance(self.index, MultiIndex):
                for i in range(self.index.nlevels):
                    arrays.append(self.index.get_level_values(i))
            else:
                arrays.append(np.asarray(self.index))

        to_remove = []
        for col in keys:
            if isinstance(col, Series):
                level = col.values
                names.append(col.name)
            elif isinstance(col, (list, np.ndarray)):
                level = col
                names.append(None)
            else:
                level = frame[col].values
                names.append(col)
                if drop:
                    to_remove.append(col)
            arrays.append(level)

        index = MultiIndex.from_arrays(arrays, names=names)

        if verify_integrity and not index.is_unique:
            duplicates = index.get_duplicates()
            raise Exception('Index has duplicate keys: %s' % duplicates)

        for c in to_remove:
            del frame[c]

        # clear up memory usage
        index._cleanup()

        frame.index = index
        return frame

    def reset_index(self, level=None, drop=False, inplace=False, col_level=0,
                    col_fill=''):
        """
        For DataFrame with multi-level index, return new DataFrame with
        labeling information in the columns under the index names, defaulting
        to 'level_0', 'level_1', etc. if any are None. For a standard index,
        the index name will be used (if set), otherwise a default 'index' or
        'level_0' (if 'index' is already taken) will be used.

        Parameters
        ----------
        level : int, str, tuple, or list, default None
            Only remove the given levels from the index. Removes all levels by
            default
        drop : boolean, default False
            Do not try to insert index into dataframe columns. This resets
            the index to the default integer index.
        inplace : boolean, default False
            Modify the DataFrame in place (do not create a new object)
        col_level : int or str, default 0
            If the columns have multiple levels, determines which level the
            labels are inserted into. By default it is inserted into the first
            level.
        col_fill : object, default ''
            If the columns have multiple levels, determines how the other levels
            are named. If None then the index name is repeated.

        Returns
        -------
        resetted : DataFrame
        """
        if inplace:
            new_obj  = self
        else:
            new_obj = self.copy()

        def _maybe_cast(values):
            if values.dtype == np.object_:
                values = lib.maybe_convert_objects(values)
            return values

        new_index = np.arange(len(new_obj))
        if isinstance(self.index, MultiIndex):
            if level is not None:
                if not isinstance(level, (tuple, list)):
                    level = [level]
                level = [self.index._get_level_number(lev) for lev in level]
                if len(level) < len(self.index.levels):
                    new_index = self.index.droplevel(level)

            if not drop:
                names = self.index.names
                zipped = zip(self.index.levels, self.index.labels)

                multi_col = isinstance(self.columns, MultiIndex)
                for i, (lev, lab) in reversed(list(enumerate(zipped))):
                    col_name = names[i]
                    if col_name is None:
                        col_name = 'level_%d' % i

                    if multi_col:
                        if col_fill is None:
                            col_name = tuple([col_name] *
                                             self.columns.nlevels)
                        else:
                            name_lst = [col_fill] * self.columns.nlevels
                            lev_num = self.columns._get_level_number(col_level)
                            name_lst[lev_num] = col_name
                            col_name = tuple(name_lst)

                    # to ndarray and maybe infer different dtype
                    level_values = _maybe_cast(lev.values)
                    if level is None or i in level:
                        new_obj.insert(0, col_name, level_values.take(lab))

        elif not drop:
            name = self.index.name
            if name is None or name == 'index':
                name = 'index' if 'index' not in self else 'level_0'
            if isinstance(self.columns, MultiIndex):
                if col_fill is None:
                    name = tuple([name] * self.columns.nlevels)
                else:
                    name_lst = [col_fill] * self.columns.nlevels
                    lev_num = self.columns._get_level_number(col_level)
                    name_lst[lev_num] = name
                    name = tuple(name_lst)
            new_obj.insert(0, name, _maybe_cast(self.index.values))

        new_obj.index = new_index
        return new_obj

    delevel = deprecate('delevel', reset_index)

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
        if isinstance(indices, list):
            indices = np.array(indices)
        if self._data.is_mixed_dtype():
            if axis == 0:
                new_data = self._data.take(indices, axis=1)
                return DataFrame(new_data)
            else:
                new_columns = self.columns.take(indices)
                return self.reindex(columns=new_columns)
        else:
            new_values = com.take_2d(self.values,
                                     com._ensure_int64(indices),
                                     axis=axis)
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
        axis : {0, 1}, or tuple/list thereof
            Pass tuple or list to drop on multiple axes
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
        if isinstance(axis, (tuple, list)):
            result = self
            for ax in axis:
                result = result.dropna(how=how, thresh=thresh,
                                       subset=subset, axis=ax)
            return result

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
            agg_obj = self.reindex(**{agg_axis_name: subset})

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
        return self.reindex(**{axis_name: new_labels})

    def drop_duplicates(self, cols=None, take_last=False, inplace=False):
        """
        Return DataFrame with duplicate rows removed, optionally only
        considering certain columns

        Parameters
        ----------
        cols : column label or sequence of labels, optional
            Only consider certain columns for identifying duplicates, by
            default use all of the columns
        take_last : boolean, default False
            Take the last observed row in a row. Defaults to the first row
        inplace : boolean, default False
            Whether to drop duplicates in place or to return a copy

        Returns
        -------
        deduplicated : DataFrame
        """

        duplicated = self.duplicated(cols, take_last=take_last)

        if inplace:
            inds, = (-duplicated).nonzero()
            self._data = self._data.take(inds)
            self._clear_item_cache()
            return self
        else:
            return self[-duplicated]

    def duplicated(self, cols=None, take_last=False):
        """
        Return boolean Series denoting duplicate rows, optionally only
        considering certain columns

        Parameters
        ----------
        cols : column label or sequence of labels, optional
            Only consider certain columns for identifying duplicates, by
            default use all of the columns
        take_last : boolean, default False
            Take the last observed row in a row. Defaults to the first row

        Returns
        -------
        duplicated : Series
        """
        # kludge for #1833
        def _m8_to_i8(x):
            if issubclass(x.dtype.type, np.datetime64):
                return x.view(np.int64)
            return x

        if cols is None:
            values = list(_m8_to_i8(self.values.T))
        else:
            if np.iterable(cols) and not isinstance(cols, basestring):
                if isinstance(cols, tuple):
                    if cols in self.columns:
                        values = [self[cols]]
                    else:
                        values = [_m8_to_i8(self[x].values) for x in cols]
                else:
                    values = [_m8_to_i8(self[x].values) for x in cols]
            else:
                values = [self[cols]]

        keys = lib.fast_zip_fillna(values)
        duplicated = lib.duplicated(keys, take_last=take_last)
        return Series(duplicated, index=self.index)

    #----------------------------------------------------------------------
    # Sorting

    def sort(self, columns=None, column=None, axis=0, ascending=True,
             inplace=False):
        """
        Sort DataFrame either by labels (along either axis) or by the values in
        column(s)

        Parameters
        ----------
        columns : object
            Column name(s) in frame. Accepts a column name or a list or tuple
            for a nested sort.
        ascending : boolean, default True
            Sort ascending vs. descending
        axis : {0, 1}
            Sort index/rows versus columns
        inplace : boolean, default False
            Sort the DataFrame without creating a new instance

        Returns
        -------
        sorted : DataFrame
        """
        if column is not None: # pragma: no cover
            import warnings
            warnings.warn("column is deprecated, use columns", FutureWarning)
            columns = column
        return self.sort_index(by=columns, axis=axis, ascending=ascending,
                               inplace=inplace)

    def sort_index(self, axis=0, by=None, ascending=True, inplace=False):
        """
        Sort DataFrame either by labels (along either axis) or by the values in
        a column

        Parameters
        ----------
        axis : {0, 1}
            Sort index/rows versus columns
        by : object
            Column name(s) in frame. Accepts a column name or a list or tuple
            for a nested sort.
        ascending : boolean, default True
            Sort ascending vs. descending
        inplace : boolean, default False
            Sort the DataFrame without creating a new instance

        Returns
        -------
        sorted : DataFrame
        """
        from pandas.core.groupby import _lexsort_indexer

        if axis not in [0, 1]:
            raise ValueError('Axis must be 0 or 1, got %s' % str(axis))

        labels = self._get_axis(axis)

        if by is not None:
            assert(axis == 0)
            if isinstance(by, (tuple, list)):
                keys = [self[x].values for x in by]
                indexer = _lexsort_indexer(keys)
            else:
                indexer = self[by].values.argsort()
        else:
            indexer = labels.argsort()

        if not ascending:
            indexer = indexer[::-1]

        if inplace:
            if axis == 1:
                self._data = self._data.reindex_items(self._data.items[indexer],
                                                      copy=False)
            elif axis == 0:
                self._data = self._data.take(indexer)

            self._clear_item_cache()
            return self
        else:
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

    def reorder_levels(self, order, axis=0):
        """
        Rearrange index levels using input order.
        May not drop or duplicate levels

        Parameters
        ----------
        order: list of int representing new level order.
               (reference level by number not by key)
        axis: where to reorder levels

        Returns
        -------
        type of caller (new object)
        """
        if not isinstance(self._get_axis(axis),
                          MultiIndex):  # pragma: no cover
            raise Exception('Can only reorder levels on a hierarchical axis.')

        result = self.copy()

        if axis == 0:
            result.index = result.index.reorder_levels(order)
        else:
            result.columns = result.columns.reorder_levels(order)
        return result

    #----------------------------------------------------------------------
    # Filling NA's

    def fillna(self, value=None, method='pad', axis=0, inplace=False,
               limit=None):
        """
        Fill NA/NaN values using the specified method

        Parameters
        ----------
        method : {'backfill', 'bfill', 'pad', 'ffill', None}, default 'pad'
            Method to use for filling holes in reindexed Series
            pad / ffill: propagate last valid observation forward to next valid
            backfill / bfill: use NEXT valid observation to fill gap
        value : scalar or dict
            Value to use to fill holes (e.g. 0), alternately a dict of values
            specifying which value to use for each column (columns not in the
            dict will not be filled)
        axis : {0, 1}, default 0
            0: fill column-by-column
            1: fill row-by-row
        inplace : boolean, default False
            If True, fill the DataFrame in place. Note: this will modify any
            other views on this DataFrame, like if you took a no-copy slice of
            an existing DataFrame, for example a column in a DataFrame. Returns
            a reference to the filled object, which is self if inplace=True
        limit : int, default None
            Maximum size gap to forward or backward fill

        See also
        --------
        reindex, asfreq

        Returns
        -------
        filled : DataFrame
        """
        self._consolidate_inplace()

        if value is None:
            if self._is_mixed_type and axis == 1:
                return self.T.fillna(method=method, limit=limit).T

            new_blocks = []
            method = com._clean_fill_method(method)
            for block in self._data.blocks:
                if block._can_hold_na:
                    newb = block.interpolate(method, axis=axis,
                                             limit=limit, inplace=inplace)
                else:
                    newb = block if inplace else block.copy()
                new_blocks.append(newb)

            new_data = BlockManager(new_blocks, self._data.axes)
        else:
            # Float type values
            if len(self.columns) == 0:
                return self
            if isinstance(value, (dict, Series)):
                if axis == 1:
                    raise NotImplementedError('Currently only can fill '
                                              'with dict/Series column '
                                              'by column')

                result = self if inplace else self.copy()
                for k, v in value.iteritems():
                    if k not in result:
                        continue
                    result[k].fillna(v, inplace=True)
                return result
            else:
                new_data = self._data.fillna(value, inplace=inplace)

        if inplace:
            self._data = new_data
            return self
        else:
            return self._constructor(new_data)

    def replace(self, to_replace, value=None, method='pad', axis=0,
                inplace=False, limit=None):
        """
        Replace values given in 'to_replace' with 'value' or using 'method'

        Parameters
        ----------
        value : scalar or dict, default None
            Value to use to fill holes (e.g. 0), alternately a dict of values
            specifying which value to use for each column (columns not in the
            dict will not be filled)
        method : {'backfill', 'bfill', 'pad', 'ffill', None}, default 'pad'
            Method to use for filling holes in reindexed Series
            pad / ffill: propagate last valid observation forward to next valid
            backfill / bfill: use NEXT valid observation to fill gap
        axis : {0, 1}, default 0
            0: fill column-by-column
            1: fill row-by-row
        inplace : boolean, default False
            If True, fill the DataFrame in place. Note: this will modify any
            other views on this DataFrame, like if you took a no-copy slice of
            an existing DataFrame, for example a column in a DataFrame. Returns
            a reference to the filled object, which is self if inplace=True
        limit : int, default None
            Maximum size gap to forward or backward fill

        See also
        --------
        reindex, asfreq

        Returns
        -------
        filled : DataFrame
        """
        self._consolidate_inplace()

        if value is None:
            return self._interpolate(to_replace, method, axis, inplace, limit)
        else:
            if len(self.columns) == 0:
                return self

            if isinstance(to_replace, dict):
                if isinstance(value, dict): # {'A' : np.nan} -> {'A' : 0}
                    return self._replace_both_dict(to_replace, value, inplace)

                elif not isinstance(value, (list, np.ndarray)):
                    return self._replace_src_dict(to_replace, value, inplace)

                raise ValueError('Fill value must be scalar or dict')

            elif isinstance(to_replace, (list, np.ndarray)):
                # [np.nan, ''] -> [0, 'missing']
                if isinstance(value, (list, np.ndarray)):
                    if len(to_replace) != len(value):
                        raise ValueError('Replacement lists must match '
                                         'in length. Expecting %d got %d ' %
                                         (len(to_replace), len(value)))

                    new_data = self._data if inplace else self.copy()._data
                    new_data._replace_list(to_replace, value)

                else: # [np.nan, ''] -> 0
                    new_data = self._data.replace(to_replace, value,
                                                  inplace=inplace)

                if inplace:
                    self._data = new_data
                    return self
                else:
                    return self._constructor(new_data)
            else:
                if isinstance(value, dict): # np.nan -> {'A' : 0, 'B' : -1}
                    return self._replace_dest_dict(to_replace, value, inplace)
                elif not isinstance(value, (list, np.ndarray)): # np.nan -> 0
                    new_data = self._data.replace(to_replace, value,
                                                  inplace=inplace)
                    if inplace:
                        self._data = new_data
                        return self
                    else:
                        return self._constructor(new_data)

            raise ValueError('Invalid to_replace type: %s' %
                             type(to_replace)) # pragma: no cover

    def _interpolate(self, to_replace, method, axis, inplace, limit):
        if self._is_mixed_type and axis == 1:
            return self.T.replace(to_replace, method=method, limit=limit).T

        method = com._clean_fill_method(method)

        if isinstance(to_replace, dict):
            if axis == 1:
                return self.T.replace(to_replace, method=method,
                                      limit=limit).T

            rs = self if inplace else self.copy()
            for k, v in to_replace.iteritems():
                if k in rs:
                    rs[k].replace(v, method=method, limit=limit,
                                  inplace=True)
            return rs

        else:
            new_blocks = []
            for block in self._data.blocks:
                newb = block.interpolate(method, axis=axis,
                                         limit=limit, inplace=inplace,
                                         missing=to_replace)
                new_blocks.append(newb)
            new_data = BlockManager(new_blocks, self._data.axes)

            if inplace:
                self._data = new_data
                return self
            else:
                return self._constructor(new_data)

    def _replace_dest_dict(self, to_replace, value, inplace):
        rs = self if inplace else self.copy()
        for k, v in value.iteritems():
            if k in rs:
                rs[k].replace(to_replace, v, inplace=True)
        return rs

    def _replace_src_dict(self, to_replace, value, inplace):
        rs = self if inplace else self.copy()
        for k, src in to_replace.iteritems():
            if k in rs:
                rs[k].replace(src, value, inplace=True)
        return rs

    def _replace_both_dict(self, to_replace, value, inplace):
        rs = self if inplace else self.copy()
        for c, src in to_replace.iteritems():
            if c in value and c in rs:
                rs[c].replace(src, value[c], inplace=True)
        return rs

    #----------------------------------------------------------------------
    # Rename

    def rename(self, index=None, columns=None, copy=True, inplace=False):
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
        inplace : boolean, default False
            Whether to return a new DataFrame. If True then value of copy is
            ignored.

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

        result = self if inplace else self.copy(deep=copy)

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

    def _combine_frame(self, other, func, fill_value=None, level=None):
        this, other = self.align(other, join='outer', level=level, copy=False)
        new_index, new_columns = this.index, this.columns

        def _arith_op(left, right):
            if fill_value is not None:
                left_mask = isnull(left)
                right_mask = isnull(right)
                left = left.copy()
                right = right.copy()

                # one but not both
                mask = left_mask ^ right_mask
                left[left_mask & mask] = fill_value
                right[right_mask & mask] = fill_value

            return func(left, right)

        if this._is_mixed_type or other._is_mixed_type:
            # XXX no good for duplicate columns
            result = {}
            for col in this:
                result[col] = func(this[col].values, other[col].values)
        else:
            result = _arith_op(this.values, other.values)

        return self._constructor(result, index=new_index,
                                 columns=new_columns, copy=False)

    def _indexed_same(self, other):
        same_index = self.index.equals(other.index)
        same_columns = self.columns.equals(other.columns)
        return same_index and same_columns

    def _combine_series(self, other, func, fill_value=None, axis=None,
                        level=None):
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
        if self.index.is_all_dates and other.index.is_all_dates:
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
        if self.empty:
            return self

        result_values = func(self.values, other)

        if not isinstance(result_values, np.ndarray):
            raise TypeError('Could not compare %s with DataFrame values'
                            % repr(other))

        return self._constructor(result_values, index=self.index,
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

    def _flex_compare_frame(self, other, func, level):
        if not self._indexed_same(other):
            self, other = self.align(other, 'outer', level=level)

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
        if other.empty:
            return self.copy()

        if self.empty:
            return other.copy()

        this, other = self.align(other, copy=False)
        new_index = this.index

        # sorts if possible
        new_columns = this.columns.union(other.columns)
        do_fill = fill_value is not None

        result = {}
        for col in new_columns:
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

    def update(self, other, join='left', overwrite=True, filter_func=None,
                     raise_conflict=False):
        """
        Modify DataFrame in place using non-NA values from passed
        DataFrame. Aligns on indices

        Parameters
        ----------
        other : DataFrame, or object coercible into a DataFrame
        join : {'left', 'right', 'outer', 'inner'}, default 'left'
        overwrite : boolean, default True
            If True then overwrite values for common keys in the calling frame
        filter_func : callable(1d-array) -> 1d-array<boolean>, default None
            Can choose to replace values other than NA. Return True for values
            that should be updated
        raise_conflict : bool
            If True, will raise an error if the DataFrame and other both
            contain data in the same place.
        """
        if join != 'left':
            raise NotImplementedError

        if not isinstance(other, DataFrame):
            other = DataFrame(other)

        other = other.reindex_like(self)

        for col in self.columns:
            this = self[col].values
            that = other[col].values
            if filter_func is not None:
                mask = -filter_func(this) | isnull(that)
            else:
                if raise_conflict:
                    mask_this = notnull(that)
                    mask_that = notnull(this)
                    if any(mask_this & mask_that):
                        raise DataConflictError("Data overlaps.")

                if overwrite:
                    mask = isnull(that)
                else:
                    mask = notnull(this)
            self[col] = np.where(mask, this, that)

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
        Pivot a level of the (possibly hierarchical) column labels, returning a
        DataFrame (or Series in the case of an object with a single level of
        column labels) having a hierarchical index with a new inner-most level
        of row labels.

        Parameters
        ----------
        level : int, string, or list of these, default last level
            Level(s) to stack, can pass level name
        dropna : boolean, default True
            Whether to drop rows in the resulting Frame/Series with no valid
            values

        Examples
        ----------
        >>> s
             a   b
        one  1.  2.
        two  3.  4.

        >>> s.stack()
        one a    1
            b    2
        two a    3
            b    4

        Returns
        -------
        stacked : DataFrame or Series
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
        Pivot a level of the (necessarily hierarchical) index labels, returning
        a DataFrame having a new level of column labels whose inner-most level
        consists of the pivoted index labels. If the index is not a MultiIndex,
        the output will be a Series (the analogue of stack when the columns are
        not a MultiIndex)

        Parameters
        ----------
        level : int, string, or list of these, default last level
            Level(s) of index to unstack, can pass level name

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

        >>> df = s.unstack(level=0)
        >>> df
           one  two
        a  1.   2.
        b  3.   4.

        >>> df.unstack()
        one  a  1.
             b  3.
        two  a  2.
             b  4.

        Returns
        -------
        unstacked : DataFrame or Series
        """
        from pandas.core.reshape import unstack
        return unstack(self, level)

    #----------------------------------------------------------------------
    # Time series-related

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

    def shift(self, periods=1, freq=None, **kwds):
        """
        Shift the index of the DataFrame by desired number of periods with an
        optional time freq

        Parameters
        ----------
        periods : int
            Number of periods to move, can be positive or negative
        freq : DateOffset, timedelta, or time rule string, optional
            Increment to use from datetools module or time rule (e.g. 'EOM')

        Notes
        -----
        If freq is specified then the index values are shifted but the data
        if not realigned

        Returns
        -------
        shifted : DataFrame
        """
        from pandas.core.series import _resolve_offset

        if periods == 0:
            return self

        offset = _resolve_offset(freq, kwds)

        if isinstance(offset, basestring):
            offset = datetools.to_offset(offset)

        def _shift_block(blk, indexer):
            new_values = blk.values.take(indexer, axis=1)
            # convert integer to float if necessary. need to do a lot more than
            # that, handle boolean etc also
            new_values = com._maybe_upcast(new_values)
            if periods > 0:
                new_values[:, :periods] = nan
            else:
                new_values[:, periods:] = nan
            return make_block(new_values, blk.items, blk.ref_items)

        if offset is None:
            indexer = self._shift_indexer(periods)
            new_blocks = [_shift_block(b, indexer) for b in self._data.blocks]
            new_data = BlockManager(new_blocks, [self.columns, self.index])
        elif isinstance(self.index, PeriodIndex):
            orig_offset = datetools.to_offset(self.index.freq)
            if offset == orig_offset:
                new_data = self._data.copy()
                new_data.axes[1] = self.index.shift(periods)
            else:
                msg = ('Given freq %s does not match PeriodIndex freq %s' %
                       (offset.rule_code, orig_offset.rule_code))
                raise ValueError(msg)
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
        (axis=0) or the columns (axis=1). Return type depends on whether passed
        function aggregates

        Parameters
        ----------
        func : function
            Function to apply to each column
        axis : {0, 1}
            0 : apply function to each column
            1 : apply function to each row
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
        To apply a function elementwise, use applymap

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
                        return Series(np.nan,
                                      index=self._get_agg_axis(axis))
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

            labels = self._get_agg_axis(axis)
            result = lib.reduce(values, func, axis=axis, dummy=dummy,
                                labels=labels)
            return Series(result, index=self._get_agg_axis(axis))
        except Exception:
            pass

        if axis == 0:
            series_gen = (self.icol(i) for i in range(len(self.columns)))
            res_index = self.columns
            res_columns = self.index
        elif axis == 1:
            res_index = self.index
            res_columns = self.columns
            values = self.values
            series_gen = (Series.from_array(arr, index=res_columns, name=name)
                          for i, (arr, name) in
                          enumerate(izip(values, res_index)))
        else:
            raise ValueError('Axis must be 0 or 1, got %s' % str(axis))

        keys = []
        results = {}
        if ignore_failures:
            successes = []
            for i, v in enumerate(series_gen):
                try:
                    results[i] = func(v)
                    keys.append(v.name)
                    successes.append(i)
                except Exception:
                    pass
            # so will work with MultiIndex, need test
            if len(successes) < len(res_index):
                res_index = res_index.take(successes)
        else:
            try:
                for i, v in enumerate(series_gen):
                    results[i] = func(v)
                    keys.append(v.name)
            except Exception, e:
                try:
                    if hasattr(e, 'args'):
                        k = res_index[i]
                        e.args = e.args + ('occurred at index %s' % str(k),)
                except NameError: # pragma: no cover
                    # no k defined yet
                    pass
                raise

        if len(results) > 0 and _is_sequence(results[0]):
            if not isinstance(results[0], Series):
                index = res_columns
            else:
                index = None

            result = self._constructor(data=results, index=index)
            result.rename(columns=dict(zip(range(len(res_index)), res_index)),
                                       inplace=True)

            if axis == 1:
                result = result.T
            result = result.convert_objects()

            return result
        else:
            s = Series(results)
            s.index = res_index
            return s

    def _apply_broadcast(self, func, axis):
        if axis == 0:
            target = self
        elif axis == 1:
            target = self.T
        else:
            raise ValueError('Axis must be 0 or 1, got %s' % str(axis))

        result_values = np.empty_like(target.values)
        columns = target.columns
        for i, col in enumerate(columns):
            result_values[:, i] = func(target[col])

        result = self._constructor(result_values, index=target.index,
                                   columns=target.columns)

        if axis == 1:
            result = result.T

        return result

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

    def append(self, other, ignore_index=False, verify_integrity=False):
        """
        Append columns of other to end of this frame's columns and index,
        returning a new object.  Columns not in this frame are added as new
        columns.

        Parameters
        ----------
        other : DataFrame or list of Series/dict-like objects
        ignore_index : boolean, default False
            If True do not use the index labels. Useful for gluing together
            record arrays
        verify_integrity : boolean, default False
            If True, raise Exception on creating index with duplicates

        Notes
        -----
        If a list of dict is passed and the keys are all contained in the
        DataFrame's index, the order of the columns in the resulting DataFrame
        will be unchanged

        Returns
        -------
        appended : DataFrame
        """
        if isinstance(other, (Series, dict)):
            if isinstance(other, dict):
                other = Series(other)
            if other.name is None and not ignore_index:
                raise Exception('Can only append a Series if '
                                'ignore_index=True')

            index = None if other.name is None else [other.name]
            other = other.reindex(self.columns, copy=False)
            other = DataFrame(other.values.reshape((1, len(other))),
                              index=index, columns=self.columns)
        elif isinstance(other, list) and not isinstance(other[0], DataFrame):
            other = DataFrame(other)
            if (self.columns.get_indexer(other.columns) >= 0).all():
                other = other.ix[:, self.columns]

        from pandas.tools.merge import concat
        if isinstance(other, (list, tuple)):
            to_concat = [self] + other
        else:
            to_concat = [self, other]
        return concat(to_concat, ignore_index=ignore_index,
                      verify_integrity=verify_integrity)

    def join(self, other, on=None, how='left', lsuffix='', rsuffix='',
             sort=False):
        """
        Join columns with other DataFrame either on index or on a key
        column. Efficiently Join multiple DataFrame objects by index at once by
        passing a list.

        Parameters
        ----------
        other : DataFrame, Series with name field set, or list of DataFrame
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
        sort : boolean, default False
            Order result DataFrame lexicographically by the join key. If False,
            preserves the index order of the calling (left) DataFrame

        Notes
        -----
        on, lsuffix, and rsuffix options are not supported when passing a list
        of DataFrame objects

        Returns
        -------
        joined : DataFrame
        """
        # For SparseDataFrame's benefit
        return self._join_compat(other, on=on, how=how, lsuffix=lsuffix,
                                 rsuffix=rsuffix, sort=sort)

    def _join_compat(self, other, on=None, how='left', lsuffix='', rsuffix='',
                     sort=False):
        from pandas.tools.merge import merge, concat

        if isinstance(other, Series):
            assert(other.name is not None)
            other = DataFrame({other.name: other})

        if isinstance(other, DataFrame):
            return merge(self, other, left_on=on, how=how,
                         left_index=on is None, right_index=True,
                         suffixes=(lsuffix, rsuffix), sort=sort)
        else:
            if on is not None:
                raise ValueError('Joining multiple DataFrames only supported'
                                 ' for joining on index')

            # join indexes only using concat
            if how == 'left':
                how = 'outer'
                join_axes = [self.index]
            else:
                join_axes = None

            frames = [self] + list(other)

            can_concat = all(df.index.is_unique for df in frames)

            if can_concat:
                return concat(frames, axis=1, join=how, join_axes=join_axes,
                              verify_integrity=True)

            joined = frames[0]

            for frame in frames[1:]:
                joined = merge(joined, frame, how=how,
                               left_index=True, right_index=True)

            return joined

    @Substitution('')
    @Appender(_merge_doc, indents=2)
    def merge(self, right, how='inner', on=None, left_on=None, right_on=None,
              left_index=False, right_index=False, sort=True,
              suffixes=('_x', '_y'), copy=True):
        from pandas.tools.merge import merge
        return merge(self, right, how=how, on=on,
                     left_on=left_on, right_on=right_on,
                     left_index=left_index, right_index=right_index, sort=sort,
                     suffixes=suffixes, copy=copy)

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
        numeric_df = self._get_numeric_data()
        cols = numeric_df.columns
        mat = numeric_df.values

        if method == 'pearson':
            correl = lib.nancorr(com._ensure_float64(mat))
        else:
            mat = mat.T
            corrf = nanops.get_corr_func(method)
            K = len(cols)
            correl = np.empty((K, K), dtype=float)
            mask = np.isfinite(mat)
            for i, ac in enumerate(mat):
                for j, bc  in enumerate(mat):
                    valid = mask[i] & mask[j]
                    if not valid.any():
                        c = np.nan
                    elif not valid.all():
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

        y contains the covariance matrix of the DataFrame's time series.
        The covariance is normalized by N-1 (unbiased estimator).
        """
        numeric_df = self._get_numeric_data()
        cols = numeric_df.columns
        mat = numeric_df.values

        if notnull(mat).all():
            baseCov = np.cov(mat.T)
        else:
            baseCov = lib.nancorr(com._ensure_float64(mat), cov=True)

        return self._constructor(baseCov, index=cols, columns=cols)

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
        if isinstance(other, Series):
            return self.apply(other.corr, axis=axis)

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

    def describe(self, percentile_width=50):
        """
        Generate various summary statistics of each column, excluding
        NaN values. These include: count, mean, std, min, max, and
        lower%/50%/upper% percentiles

        Parameters
        ----------
        percentile_width : float, optional
            width of the desired uncertainty interval, default is 50,
            which corresponds to lower=25, upper=75

        Returns
        -------
        DataFrame of summary statistics
        """
        numdata = self._get_numeric_data()

        if len(numdata.columns) == 0:
            return DataFrame(dict((k, v.describe())
                                  for k, v in self.iteritems()),
                                  columns=self.columns)

        lb = .5 * (1. - percentile_width/100.)
        ub = 1. - lb

        def pretty_name(x):
            x *= 100
            if x == int(x):
                return '%.0f%%' % x
            else:
                return '%.1f%%' % x

        destat_columns = ['count', 'mean', 'std', 'min',
                          pretty_name(lb), '50%', pretty_name(ub),
                          'max']

        destat = []

        for column in numdata.columns:
            series = self[column]
            ser_desc = series.describe()
            destat.append([series.count(), series.mean(), series.std(),
                           series.min(), series.quantile(lb), series.median(),
                           series.quantile(ub), series.max()])

        return self._constructor(map(list, zip(*destat)), index=destat_columns,
                                 columns=numdata.columns)

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
            frame = self._get_numeric_data()
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

        return result

    def _count_level(self, level, axis=0, numeric_only=False):
        if numeric_only:
            frame = self._get_numeric_data()
        else:
            frame = self

        if axis == 1:
            frame = frame.T

        # python 2.5
        mask = notnull(frame.values).view(np.uint8)

        if isinstance(level, basestring):
            level = self.index._get_level_number(level)

        level_index = frame.index.levels[level]
        labels = com._ensure_int64(frame.index.labels[level])
        counts = lib.count_level_2d(mask, labels, len(level_index))

        result = DataFrame(counts, index=level_index,
                           columns=frame.columns)

        if axis == 1:
            return result.T
        else:
            return result

    def any(self, axis=0, bool_only=None, skipna=True, level=None):
        """
        Return whether any element is True over requested axis.
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
        bool_only : boolean, default None
            Only include boolean data.

        Returns
        -------
        any : Series (or DataFrame if level specified)
        """
        if level is not None:
            return self._agg_by_level('any', axis=axis, level=level,
                                      skipna=skipna)
        return self._reduce(nanops.nanany, axis=axis, skipna=skipna,
                            numeric_only=bool_only, filter_type='bool')

    def all(self, axis=0, bool_only=None, skipna=True, level=None):
        """
        Return whether any element is True over requested axis.
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
        bool_only : boolean, default None
            Only include boolean data.

        Returns
        -------
        any : Series (or DataFrame if level specified)
        """
        if level is not None:
            return self._agg_by_level('all', axis=axis, level=level,
                                      skipna=skipna)
        return self._reduce(nanops.nanall, axis=axis, skipna=skipna,
                            numeric_only=bool_only, filter_type='bool')

    @Substitution(name='sum', shortname='sum', na_action=_doc_exclude_na,
                  extras=_numeric_only_doc)
    @Appender(_stat_doc)
    def sum(self, axis=0, numeric_only=None, skipna=True, level=None):
        if level is not None:
            return self._agg_by_level('sum', axis=axis, level=level,
                                      skipna=skipna)
        return self._reduce(nanops.nansum, axis=axis, skipna=skipna,
                            numeric_only=numeric_only)

    @Substitution(name='mean', shortname='mean', na_action=_doc_exclude_na,
                  extras='')
    @Appender(_stat_doc)
    def mean(self, axis=0, skipna=True, level=None):
        if level is not None:
            return self._agg_by_level('mean', axis=axis, level=level,
                                      skipna=skipna)
        return self._reduce(nanops.nanmean, axis=axis, skipna=skipna,
                            numeric_only=None)

    @Substitution(name='minimum', shortname='min', na_action=_doc_exclude_na,
                  extras='')
    @Appender(_stat_doc)
    def min(self, axis=0, skipna=True, level=None):
        if level is not None:
            return self._agg_by_level('min', axis=axis, level=level,
                                      skipna=skipna)
        return self._reduce(nanops.nanmin, axis=axis, skipna=skipna,
                            numeric_only=None)

    @Substitution(name='maximum', shortname='max', na_action=_doc_exclude_na,
                  extras='')
    @Appender(_stat_doc)
    def max(self, axis=0, skipna=True, level=None):
        if level is not None:
            return self._agg_by_level('max', axis=axis, level=level,
                                      skipna=skipna)
        return self._reduce(nanops.nanmax, axis=axis, skipna=skipna,
                            numeric_only=None)

    @Substitution(name='product', shortname='product',
                  na_action='NA/null values are treated as 1', extras='')
    @Appender(_stat_doc)
    def prod(self, axis=0, skipna=True, level=None):
        if level is not None:
            return self._agg_by_level('prod', axis=axis, level=level,
                                      skipna=skipna)
        return self._reduce(nanops.nanprod, axis=axis, skipna=skipna,
                            numeric_only=None)

    product = prod

    @Substitution(name='median', shortname='median', na_action=_doc_exclude_na,
                  extras='')
    @Appender(_stat_doc)
    def median(self, axis=0, skipna=True, level=None):
        if level is not None:
            return self._agg_by_level('median', axis=axis, level=level,
                                      skipna=skipna)
        return self._reduce(nanops.nanmedian, axis=axis, skipna=skipna,
                            numeric_only=None)

    @Substitution(name='mean absolute deviation', shortname='mad',
                  na_action=_doc_exclude_na, extras='')
    @Appender(_stat_doc)
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

    @Substitution(name='variance', shortname='var',
                  na_action=_doc_exclude_na, extras='')
    @Appender(_stat_doc +
        """
        Normalized by N-1 (unbiased estimator).
        """)
    def var(self, axis=0, skipna=True, level=None, ddof=1):
        if level is not None:
            return self._agg_by_level('var', axis=axis, level=level,
                                      skipna=skipna, ddof=ddof)
        return self._reduce(nanops.nanvar, axis=axis, skipna=skipna,
                            numeric_only=None, ddof=ddof)

    @Substitution(name='standard deviation', shortname='std',
                  na_action=_doc_exclude_na, extras='')
    @Appender(_stat_doc +
        """
        Normalized by N-1 (unbiased estimator).
        """)
    def std(self, axis=0, skipna=True, level=None, ddof=1):
        if level is not None:
            return self._agg_by_level('std', axis=axis, level=level,
                                      skipna=skipna, ddof=ddof)
        return np.sqrt(self.var(axis=axis, skipna=skipna, ddof=ddof))

    @Substitution(name='unbiased skewness', shortname='skew',
                  na_action=_doc_exclude_na, extras='')
    @Appender(_stat_doc)
    def skew(self, axis=0, skipna=True, level=None):
        if level is not None:
            return self._agg_by_level('skew', axis=axis, level=level,
                                      skipna=skipna)
        return self._reduce(nanops.nanskew, axis=axis, skipna=skipna,
                            numeric_only=None)


    @Substitution(name='unbiased kurtosis', shortname='kurt',
                  na_action=_doc_exclude_na, extras='')
    @Appender(_stat_doc)
    def kurt(self, axis=0, skipna=True, level=None):
        if level is not None:
            return self._agg_by_level('kurt', axis=axis, level=level,
                                      skipna=skipna)
        return self._reduce(nanops.nankurt, axis=axis, skipna=skipna,
                            numeric_only=None)

    def _agg_by_level(self, name, axis=0, level=0, skipna=True, **kwds):
        grouped = self.groupby(level=level, axis=axis)
        if hasattr(grouped, name) and skipna:
            return getattr(grouped, name)(**kwds)
        method = getattr(type(self), name)
        applyf = lambda x: method(x, axis=axis, skipna=skipna, **kwds)
        return grouped.aggregate(applyf)

    def _reduce(self, op, axis=0, skipna=True, numeric_only=None,
                filter_type=None, **kwds):
        f = lambda x: op(x, axis=axis, skipna=skipna, **kwds)
        labels = self._get_agg_axis(axis)
        if numeric_only is None:
            try:
                values = self.values
                result = f(values)
            except Exception:
                if filter_type is None or filter_type == 'numeric':
                    data = self._get_numeric_data()
                elif filter_type == 'bool':
                    data = self._get_bool_data()
                else:
                    raise ValueError('Invalid filter_type %s ' %
                                     str(filter_type))
                result = f(data.values)
                labels = data._get_agg_axis(axis)
        else:
            if numeric_only:
                if filter_type is None or filter_type == 'numeric':
                    data = self._get_numeric_data()
                elif filter_type == 'bool':
                    data = self._get_bool_data()
                else:
                    raise ValueError('Invalid filter_type %s ' %
                                     str(filter_type))
                values = data.values
                labels = data._get_agg_axis(axis)
            else:
                values = self.values
            result = f(values)

        if result.dtype == np.object_:
            try:
                if filter_type is None or filter_type == 'numeric':
                    result = result.astype(np.float64)
                elif filter_type == 'bool' and notnull(result).all():
                    result = result.astype(np.bool_)
                else:
                    raise ValueError('Invalid dtype %s ' % str(filter_type))

            except (ValueError, TypeError):
                pass

        return Series(result, index=labels)

    def idxmin(self, axis=0, skipna=True):
        """
        Return index of first occurrence of minimum over requested axis.
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
        Return index of first occurrence of maximum over requested axis.
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

    def _get_agg_axis(self, axis_num):
        if axis_num == 0:
            return self.columns
        elif axis_num == 1:
            return self.index
        else:
            raise Exception('Must have 0<= axis <= 1')

    def _get_numeric_data(self):
        if self._is_mixed_type:
            num_data = self._data.get_numeric_data()
            return DataFrame(num_data, copy=False)
        else:
            if (self.values.dtype != np.object_ and
                not issubclass(self.values.dtype.type, np.datetime64)):
                return self
            else:
                return self.ix[:, []]

    def _get_bool_data(self):
        if self._is_mixed_type:
            bool_data = self._data.get_bool_data()
            return DataFrame(bool_data, copy=False)
        else:
            if self.values.dtype == np.bool_:
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
        per = q * 100

        def f(arr):
            arr = arr.values
            if arr.dtype != np.float_:
                arr = arr.astype(float)
            arr = arr[notnull(arr)]
            if len(arr) == 0:
                return nan
            else:
                return _quantile(arr, per)

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

    def rank(self, axis=0, numeric_only=None, method='average',
             na_option='keep', ascending=True):
        """
        Compute numerical data ranks (1 through n) along axis. Equal values are
        assigned a rank that is the average of the ranks of those values

        Parameters
        ----------
        axis : {0, 1}, default 0
            Ranks over columns (0) or rows (1)
        numeric_only : boolean, default None
            Include only float, int, boolean data
        method : {'average', 'min', 'max', 'first'}
            average: average rank of group
            min: lowest rank in group
            max: highest rank in group
            first: ranks assigned in order they appear in the array
        na_option : {'keep'}
            keep: leave NA values where they are
        ascending : boolean, default True
            False for ranks by high (1) to low (N)

        Returns
        -------
        ranks : DataFrame
        """
        if numeric_only is None:
            try:
                ranks = algos.rank(self.values, axis=axis, method=method,
                                   ascending=ascending)
                return DataFrame(ranks, index=self.index, columns=self.columns)
            except TypeError:
                numeric_only = True

        if numeric_only:
            data = self._get_numeric_data()
        else:
            data = self
        ranks = algos.rank(data.values, axis=axis, method=method,
                           ascending=ascending)
        return DataFrame(ranks, index=data.index, columns=data.columns)

    def to_timestamp(self, freq=None, how='start', axis=0, copy=True):
        """
        Cast to DatetimeIndex of timestamps, at *beginning* of period

        Parameters
        ----------
        freq : string, default frequency of PeriodIndex
            Desired frequency
        how : {'s', 'e', 'start', 'end'}
            Convention for converting period to timestamp; start of period
            vs. end
        axis : {0, 1} default 0
            The axis to convert (the index by default)
        copy : boolean, default True
            If false then underlying input data is not copied

        Returns
        -------
        df : DataFrame with DatetimeIndex
        """
        new_data = self._data
        if copy:
            new_data = new_data.copy()

        if axis == 0:
            new_data.set_axis(1, self.index.to_timestamp(freq=freq, how=how))
        elif axis == 1:
            new_data.set_axis(0, self.columns.to_timestamp(freq=freq, how=how))
        else:
            raise ValueError('Axis must be 0 or 1. Got %s' % str(axis))

        return DataFrame(new_data)

    def to_period(self, freq=None, axis=0, copy=True):
        """
        Convert DataFrame from DatetimeIndex to PeriodIndex with desired
        frequency (inferred from index if not passed)

        Parameters
        ----------
        freq : string, default
        axis : {0, 1}, default 0
            The axis to convert (the index by default)
        copy : boolean, default True
            If False then underlying input data is not copied

        Returns
        -------
        ts : TimeSeries with PeriodIndex
        """
        new_data = self._data
        if copy:
            new_data = new_data.copy()

        if axis == 0:
            if freq is None:
                freq = self.index.freqstr or self.index.inferred_freq
            new_data.set_axis(1, self.index.to_period(freq=freq))
        elif axis == 1:
            if freq is None:
                freq = self.columns.freqstr or self.columns.inferred_freq
            new_data.set_axis(0, self.columns.to_period(freq=freq))
        else:
            raise ValueError('Axis must be 0 or 1. Got %s' % str(axis))

        return DataFrame(new_data)

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

        result[i] = f(values[left_bound:right_bound])

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
    from pandas.core.index import _union_indexes

    index = None
    if len(data) == 0:
        index = Index([])
    elif len(data) > 0 and index is None:
        raw_lengths = []
        indexes = []

        have_raw_arrays = False
        have_series = False
        have_dicts = False

        for v in data.values():
            if isinstance(v, Series):
                have_series = True
                indexes.append(v.index)
            elif isinstance(v, dict):
                have_dicts = True
                indexes.append(v.keys())
            elif isinstance(v, (list, tuple, np.ndarray)):
                have_raw_arrays = True
                raw_lengths.append(len(v))

        if not indexes and not raw_lengths:
            raise ValueError('If use all scalar values, must pass index')

        if have_series or have_dicts:
            index = _union_indexes(indexes)

        if have_raw_arrays:
            lengths = list(set(raw_lengths))
            if len(lengths) > 1:
                raise ValueError('arrays must all be same length')

            if have_dicts:
                raise ValueError('Mixing dicts with non-Series may lead to '
                                 'ambiguous ordering.')

            if have_series:
                if lengths[0] != len(index):
                    msg = ('array length %d does not match index length %d'
                          % (lengths[0], len(index)))
                    raise ValueError(msg)
            else:
                index = Index(np.arange(lengths[0]))

    return _ensure_index(index)


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
        sdict = dict((k, v.values) for k, v in arr.iteritems())
    elif isinstance(arr, dict):
        columns = sorted(arr)
        sdict = arr.copy()
    else:  # pragma: no cover
        raise TypeError('%s' % type(arr))

    return columns, sdict


def _to_sdict(data, columns, coerce_float=False):
    if len(data) == 0:
        return {}, columns
    if isinstance(data[0], (list, tuple)):
        return _list_to_sdict(data, columns, coerce_float=coerce_float)
    elif isinstance(data[0], dict):
        return _list_of_dict_to_sdict(data, columns, coerce_float=coerce_float)
    elif isinstance(data[0], Series):
        return _list_of_series_to_sdict(data, columns,
                                        coerce_float=coerce_float)
    else:
        # last ditch effort
        data = map(tuple, data)
        return _list_to_sdict(data, columns, coerce_float=coerce_float)

def _list_to_sdict(data, columns, coerce_float=False):
    if len(data) > 0 and isinstance(data[0], tuple):
        content = list(lib.to_object_array_tuples(data).T)
    elif len(data) > 0:
        # list of lists
        content = list(lib.to_object_array(data).T)
    else:
        if columns is None:
            columns = []
        return {}, columns
    return _convert_object_array(content, columns,
                                 coerce_float=coerce_float)

def _list_of_series_to_sdict(data, columns, coerce_float=False):
    from pandas.core.index import _get_combined_index

    if columns is None:
        columns = _get_combined_index([s.index for s in data])

    indexer_cache = {}

    aligned_values = []
    for s in data:
        index = s.index
        if id(index) in indexer_cache:
            indexer = indexer_cache[id(index)]
        else:
            indexer = indexer_cache[id(index)] = index.get_indexer(columns)
        aligned_values.append(com.take_1d(s.values, indexer))

    values = np.vstack(aligned_values)

    if values.dtype == np.object_:
        content = list(values.T)
        return _convert_object_array(content, columns,
                                     coerce_float=coerce_float)
    else:
        return values, columns


def _list_of_dict_to_sdict(data, columns, coerce_float=False):
    if columns is None:
        gen = (x.keys() for x in data)
        columns = lib.fast_unique_multiple_list_gen(gen)

    # assure that they are of the base dict class and not of derived
    # classes
    data = [(type(d) is dict) and d or dict(d)
            for d in data]

    content = list(lib.dicts_to_array(data, list(columns)).T)
    return _convert_object_array(content, columns,
                                 coerce_float=coerce_float)


def _convert_object_array(content, columns, coerce_float=False):
    if columns is None:
        columns = range(len(content))
    else:
        if len(columns) != len(content):
            raise AssertionError('%d columns passed, passed data had %s '
                                 'columns' % (len(columns), len(content)))

    sdict = dict((c, lib.maybe_convert_objects(vals, try_float=coerce_float))
                 for c, vals in zip(columns, content))
    return sdict, columns

def _get_names_from_index(data):
    index = range(len(data))
    has_some_name = any([s.name is not None for s in data])
    if not has_some_name:
        return index

    count = 0
    for i, s in enumerate(data):
        n = s.name
        if n is not None:
            index[i] = n
        else:
            index[i] = 'Unnamed %d' % count
            count += 1

    return index

def _homogenize(data, index, columns, dtype=None):
    from pandas.core.series import _sanitize_array

    homogenized = {}

    if dtype is not None:
        dtype = np.dtype(dtype)

    oindex = None

    for k in columns:
        if k not in data:
            # no obvious "empty" int column
            if dtype is not None and issubclass(dtype.type, np.integer):
                continue

            if dtype is None:
                # #1783
                v = np.empty(len(index), dtype=object)
            else:
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
                if oindex is None:
                    oindex = index.astype('O')
                if type(v) == dict:
                    # fast cython method
                    v = lib.fast_multiget(v, oindex, default=np.nan)
                else:
                    v = lib.map_infer(oindex, v.get)

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

#----------------------------------------------------------------------
# Add plotting methods to DataFrame

import pandas.tools.plotting as gfx

DataFrame.plot = gfx.plot_frame
DataFrame.hist = gfx.hist_frame

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
    ax = plots.boxplot(self, column=column, by=by, ax=ax,
                       fontsize=fontsize, grid=grid, rot=rot, **kwds)
    plt.draw_if_interactive()
    return ax
DataFrame.boxplot = boxplot


if __name__ == '__main__':
    import nose
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
