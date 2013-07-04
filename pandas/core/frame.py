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
import operator
import sys
import collections
import itertools

from numpy import nan as NA
import numpy as np
import numpy.ma as ma

from pandas.core.common import (isnull, notnull, PandasError, _try_sort,
                                _default_index, _maybe_upcast, _is_sequence,
                                _infer_dtype_from_scalar)
from pandas.core.generic import NDFrame
from pandas.core.index import Index, MultiIndex, _ensure_index
from pandas.core.indexing import (_NDFrameIndexer, _maybe_droplevels,
                                  _convert_to_index_sliceable, _check_bool_indexer,
                                  _maybe_convert_indices)
from pandas.core.internals import (BlockManager,
                                   create_block_manager_from_arrays,
                                   create_block_manager_from_blocks)
from pandas.core.series import Series, _radd_compat
import pandas.core.expressions as expressions
from pandas.compat.scipy import scoreatpercentile as _quantile
from pandas.util.compat import OrderedDict
from pandas.util import py3compat
from pandas.util.terminal import get_terminal_size
from pandas.util.decorators import deprecate, Appender, Substitution

from pandas.tseries.period import PeriodIndex
from pandas.tseries.index import DatetimeIndex

import pandas.core.algorithms as algos
import pandas.core.datetools as datetools
import pandas.core.common as com
import pandas.core.format as fmt
import pandas.core.generic as generic
import pandas.core.nanops as nanops

import pandas.lib as lib
import pandas.tslib as tslib
import pandas.algos as _algos

from pandas.core.config import get_option, set_option

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
sort : boolean, default False
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


class DataConflictError(Exception):
    pass

#----------------------------------------------------------------------
# Factory helper methods


def _arith_method(op, name, str_rep = None, default_axis='columns', fill_zeros=None, **eval_kwargs):
    def na_op(x, y):
        try:
            result = expressions.evaluate(op, str_rep, x, y, raise_on_error=True, **eval_kwargs)
            result = com._fill_zeros(result,y,fill_zeros)

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

            result, changed = com._maybe_upcast_putmask(result,-mask,np.nan)
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


def _flex_comp_method(op, name, str_rep = None, default_axis='columns'):

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
            return self._flex_compare_frame(other, na_op, str_rep, level)

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

                return self._flex_compare_frame(casted, na_op, str_rep, level)

            else:  # pragma: no cover
                raise ValueError("Bad argument shape")

        else:
            return self._combine_const(other, na_op)

    f.__name__ = name

    return f


def _comp_method(func, name, str_rep):
    @Appender('Wrapper for comparison method %s' % name)
    def f(self, other):
        if isinstance(other, DataFrame):    # Another DataFrame
            return self._compare_frame(other, func, str_rep)
        elif isinstance(other, Series):
            return self._combine_series_infer(other, func)
        else:

            # straight boolean comparisions we want to allow all columns
            # (regardless of dtype to pass thru)
            return self._combine_const(other, func, raise_on_error = False).fillna(True).astype(bool)

    f.__name__ = name

    return f


#----------------------------------------------------------------------
# DataFrame class


class DataFrame(NDFrame):
    """ Two-dimensional size-mutable, potentially heterogeneous tabular data
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
    _auto_consolidate = True
    _het_axis = 1
    _info_axis = 'columns'
    _col_klass = Series

    _AXIS_NUMBERS = {
        'index': 0,
        'columns': 1
    }

    _AXIS_NAMES = dict((v, k) for k, v in _AXIS_NUMBERS.iteritems())

    def __init__(self, data=None, index=None, columns=None, dtype=None,
                 copy=False):
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
            if mask.any():
                data, fill_value = _maybe_upcast(data, copy=True)
                data[mask] = fill_value
            else:
                data = data.copy()
            mgr = self._init_ndarray(data, index, columns, dtype=dtype,
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

                if isinstance(data[0], (list, tuple, collections.Mapping, Series)):
                    arrays, columns = _to_arrays(data, columns, dtype=dtype)
                    columns = _ensure_index(columns)

                    if index is None:
                        index = _default_index(len(data))
                    mgr = _arrays_to_mgr(arrays, columns, index, columns,
                                         dtype=dtype)
                else:
                    mgr = self._init_ndarray(data, index, columns, dtype=dtype,
                                             copy=copy)
            else:
                mgr = self._init_ndarray(data, index, columns, dtype=dtype,
                                         copy=copy)
        else:
            try:
                arr = np.array(data, dtype=dtype, copy=copy)
            except (ValueError, TypeError):
                raise PandasError('DataFrame constructor called with '
                                  'incompatible data and dtype')

            if arr.ndim == 0 and index is not None and columns is not None:
                if isinstance(data, basestring) and dtype is None:
                    dtype = np.object_
                if dtype is None:
                    dtype, data = _infer_dtype_from_scalar(data)

                values = np.empty((len(index), len(columns)), dtype=dtype)
                values.fill(data)
                mgr = self._init_ndarray(values, index, columns, dtype=dtype,
                                         copy=False)
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
        if dtype is not None:
            dtype = np.dtype(dtype)

        if columns is not None:
            columns = _ensure_index(columns)

            # prefilter if columns passed

            data = dict((k, v) for k, v in data.iteritems() if k in columns)

            if index is None:
                index = extract_index(data.values())
            else:
                index = _ensure_index(index)

            arrays = []
            data_names = []
            for k in columns:
                if k not in data:
                    # no obvious "empty" int column
                    if dtype is not None and issubclass(dtype.type,
                                                        np.integer):
                        continue

                    if dtype is None:
                        # #1783
                        v = np.empty(len(index), dtype=object)
                    else:
                        v = np.empty(len(index), dtype=dtype)

                    v.fill(NA)
                else:
                    v = data[k]
                data_names.append(k)
                arrays.append(v)
        else:
            keys = data.keys()
            if not isinstance(data, OrderedDict):
                keys = _try_sort(data.keys())
            columns = data_names = Index(keys)
            arrays = [data[k] for k in columns]

        return _arrays_to_mgr(arrays, data_names, index, columns,
                              dtype=dtype)

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

        return create_block_manager_from_blocks([ values.T ], [ columns, index ])

    def _wrap_array(self, arr, axes, copy=False):
        index, columns = axes
        return self._constructor(arr, index=index, columns=columns, copy=copy)

    @property
    def _verbose_info(self):
        import warnings
        warnings.warn('The _verbose_info property will be removed in version '
                      '0.13. please use "max_info_rows"', FutureWarning)
        return get_option('display.max_info_rows') is None

    @_verbose_info.setter
    def _verbose_info(self, value):
        import warnings
        warnings.warn('The _verbose_info property will be removed in version '
                      '0.13. please use "max_info_rows"', FutureWarning)

        value = None if value else 1000000
        set_option('display.max_info_rows', value)

    @property
    def axes(self):
        return [self.index, self.columns]

    @property
    def shape(self):
        return (len(self.index), len(self.columns))

    #----------------------------------------------------------------------
    # Class behavior
    def _repr_fits_vertical_(self):
        """
        Check length against max_rows.
        """
        max_rows = get_option("display.max_rows")
        return len(self) <= max_rows

    def _repr_fits_horizontal_(self,ignore_width=False):
        """
        Check if full repr fits in horizontal boundaries imposed by the display
        options width and max_columns. In case off non-interactive session, no
        boundaries apply.

        ignore_width is here so ipnb+HTML output can behave the way
        users expect. display.max_columns remains in effect.
        GH3541, GH3573
        """

        width, height = fmt.get_console_size()
        max_columns = get_option("display.max_columns")
        nb_columns = len(self.columns)

        # exceed max columns
        if ((max_columns and nb_columns > max_columns) or
            ((not ignore_width) and width and nb_columns > (width // 2))):
            return False

        if (ignore_width  # used by repr_html under IPython notebook
            or not com.in_interactive_session()): # scripts ignore terminal dims
            return True

        if (get_option('display.width') is not None or
            com.in_ipython_frontend()):
            # check at least the column row for excessive width
            max_rows = 1
        else:
            max_rows = get_option("display.max_rows")

        # when auto-detecting, so width=None and not in ipython front end
        # check whether repr fits horizontal by actualy checking
        # the width of the rendered repr
        buf = StringIO()

        # only care about the stuff we'll actually print out
        # and to_string on entire frame may be expensive
        d = self

        if not (max_rows is None): # unlimited rows
            # min of two, where one may be None
            d=d.iloc[:min(max_rows,len(d))]
        else:
            return True

        d.to_string(buf=buf)
        value = buf.getvalue()
        repr_width = max([len(l) for l in value.split('\n')])

        return repr_width < width

    def __unicode__(self):
        """
        Return a string representation for a particular DataFrame

        Invoked by unicode(df) in py2 only. Yields a Unicode String in both
        py2/py3.
        """
        buf = StringIO(u"")
        fits_vertical = self._repr_fits_vertical_()
        fits_horizontal = False
        if fits_vertical:
            # This needs to compute the entire repr
            # so don't do it unless rownum is bounded
            fits_horizontal = self._repr_fits_horizontal_()

        if fits_vertical and fits_horizontal:
            self.to_string(buf=buf)
        else:
            width, _ = fmt.get_console_size()
            max_columns = get_option("display.max_columns")
            expand_repr = get_option("display.expand_frame_repr")
            # within max_cols and max_rows, but cols exceed width
            # of terminal, then use expand_repr
            if (fits_vertical and
                expand_repr and
                len(self.columns) <= max_columns):
                self.to_string(buf=buf, line_width=width)
            else:
                max_info_rows = get_option('display.max_info_rows')
                verbose = (max_info_rows is None or
                           self.shape[0] <= max_info_rows)
                self.info(buf=buf, verbose=verbose)

        value = buf.getvalue()
        if not  type(value) == unicode:
            raise AssertionError()

        return value

    def _repr_html_(self):
        """
        Return a html representation for a particular DataFrame.
        Mainly for IPython notebook.
        """
        # ipnb in html repr mode allows scrolling
        # users strongly prefer to h-scroll a wide HTML table in the browser
        # then to get a summary view. GH3541, GH3573
        ipnbh = com.in_ipnb() and get_option('display.notebook_repr_html')

        # qtconsole doesn't report it's line width, and also
        # behaves badly when outputting an HTML table
        # that doesn't fit the window, so disable it.
        if com.in_qtconsole():
            raise ValueError('Disable HTML output in QtConsole')

        if get_option("display.notebook_repr_html"):
            fits_vertical = self._repr_fits_vertical_()
            fits_horizontal = False
            if fits_vertical:
                fits_horizontal = self._repr_fits_horizontal_(ignore_width=ipnbh)

            if fits_horizontal and fits_vertical:
                return ('<div style="max-height:1000px;'
                        'max-width:1500px;overflow:auto;">\n' +
                        self.to_html() + '\n</div>')
            else:
                buf = StringIO(u"")
                max_info_rows = get_option('display.max_info_rows')
                verbose = (max_info_rows is None or
                           self.shape[0] <= max_info_rows)
                self.info(buf=buf, verbose=verbose)
                info = buf.getvalue()
                info = info.replace('&', r'&amp;')
                info = info.replace('<', r'&lt;')
                info = info.replace('>', r'&gt;')
                return ('<pre>\n' + info + '\n</pre>')
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
        if self.columns.is_unique and hasattr(self, '_item_cache'):
            for k in self.columns:
                yield k, self._get_item_cache(k)
        else:
            for i, k in enumerate(self.columns):
                yield k, self.icol(i)

    def iterrows(self):
        """
        Iterate over rows of DataFrame as (index, Series) pairs.

        Notes
        -----

        * ``iterrows`` does **not** preserve dtypes across the rows (dtypes
          are preserved across columns for DataFrames). For example,

            >>> df = DataFrame([[1, 1.0]], columns=['x', 'y'])
            >>> row = next(df.iterrows())[1]
            >>> print row['x'].dtype
            float64
            >>> print df['x'].dtype
            int64

        Returns
        -------
        it : generator
            A generator that iterates over the rows of the frame.
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

        # use integer indexing because of possible duplicate column names
        arrays.extend(self.iloc[:, k] for k in xrange(len(self.columns)))
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

    add = _arith_method(operator.add, 'add', '+')
    mul = _arith_method(operator.mul, 'multiply', '*')
    sub = _arith_method(operator.sub, 'subtract', '-')
    div = divide = _arith_method(lambda x, y: x / y, 'divide', '/')
    pow = _arith_method(operator.pow, 'pow', '**')
    mod = _arith_method(lambda x, y: x % y, 'mod')

    radd = _arith_method(_radd_compat, 'radd')
    rmul = _arith_method(operator.mul, 'rmultiply')
    rsub = _arith_method(lambda x, y: y - x, 'rsubtract')
    rdiv = _arith_method(lambda x, y: y / x, 'rdivide')
    rpow = _arith_method(lambda x, y: y ** x, 'rpow')
    rmod = _arith_method(lambda x, y: y % x, 'rmod')

    __add__ = _arith_method(operator.add, '__add__', '+', default_axis=None)
    __sub__ = _arith_method(operator.sub, '__sub__', '-', default_axis=None)
    __mul__ = _arith_method(operator.mul, '__mul__', '*', default_axis=None)
    __truediv__ = _arith_method(operator.truediv, '__truediv__', '/',
                                default_axis=None, fill_zeros=np.inf, truediv=True)
    # numexpr produces a different value (python/numpy: 0.000, numexpr: inf)
    # when dividing by zero, so can't use floordiv speed up (yet)
    # __floordiv__ = _arith_method(operator.floordiv, '__floordiv__', '//',
    __floordiv__ = _arith_method(operator.floordiv, '__floordiv__',
                                 default_axis=None, fill_zeros=np.inf)
    __pow__ = _arith_method(operator.pow, '__pow__', '**', default_axis=None)

    # currently causes a floating point exception to occur - so sticking with unaccelerated for now
    # __mod__ = _arith_method(operator.mod, '__mod__', '%', default_axis=None, fill_zeros=np.nan)
    __mod__ = _arith_method(operator.mod, '__mod__', default_axis=None, fill_zeros=np.nan)

    __radd__ = _arith_method(_radd_compat, '__radd__', default_axis=None)
    __rmul__ = _arith_method(operator.mul, '__rmul__', default_axis=None)
    __rsub__ = _arith_method(lambda x, y: y - x, '__rsub__', default_axis=None)
    __rtruediv__ = _arith_method(lambda x, y: y / x, '__rtruediv__',
                                 default_axis=None, fill_zeros=np.inf)
    __rfloordiv__ = _arith_method(lambda x, y: y // x, '__rfloordiv__',
                                  default_axis=None, fill_zeros=np.inf)
    __rpow__ = _arith_method(lambda x, y: y ** x, '__rpow__',
                             default_axis=None)
    __rmod__ = _arith_method(lambda x, y: y % x, '__rmod__', default_axis=None,
                             fill_zeros=np.nan)

    # boolean operators
    __and__ = _arith_method(operator.and_, '__and__', '&')
    __or__ = _arith_method(operator.or_, '__or__', '|')
    __xor__ = _arith_method(operator.xor, '__xor__')

    # Python 2 division methods
    if not py3compat.PY3:
        __div__ = _arith_method(operator.div, '__div__', '/',
                                default_axis=None, fill_zeros=np.inf, truediv=False)
        __rdiv__ = _arith_method(lambda x, y: y / x, '__rdiv__',
                                 default_axis=None, fill_zeros=np.inf)

    def __neg__(self):
        arr = operator.neg(self.values)
        return self._wrap_array(arr, self.axes, copy=False)

    def __invert__(self):
        arr = operator.inv(self.values)
        return self._wrap_array(arr, self.axes, copy=False)

    # Comparison methods
    __eq__ = _comp_method(operator.eq, '__eq__', '==')
    __ne__ = _comp_method(operator.ne, '__ne__', '!=')
    __lt__ = _comp_method(operator.lt, '__lt__', '<' )
    __gt__ = _comp_method(operator.gt, '__gt__', '>' )
    __le__ = _comp_method(operator.le, '__le__', '<=')
    __ge__ = _comp_method(operator.ge, '__ge__', '>=')

    eq = _flex_comp_method(operator.eq, 'eq', '==')
    ne = _flex_comp_method(operator.ne, 'ne', '!=')
    lt = _flex_comp_method(operator.lt, 'lt', '<')
    gt = _flex_comp_method(operator.gt, 'gt', '>')
    le = _flex_comp_method(operator.le, 'le', '<=')
    ge = _flex_comp_method(operator.ge, 'ge', '>=')

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
        if isinstance(other, (Series, DataFrame)):
            common = self.columns.union(other.index)
            if (len(common) > len(self.columns) or
                    len(common) > len(other.index)):
                raise ValueError('matrices are not aligned')

            left = self.reindex(columns=common, copy=False)
            right = other.reindex(index=common, copy=False)
            lvals = left.values
            rvals = right.values
        else:
            left = self
            lvals = self.values
            rvals = np.asarray(other)
            if lvals.shape[1] != rvals.shape[0]:
                raise Exception('Dot product shape mismatch, %s vs %s' %
                                (lvals.shape, rvals.shape))

        if isinstance(other, DataFrame):
            return self._constructor(np.dot(lvals, rvals),
                                     index=left.index,
                                     columns=other.columns)
        elif isinstance(other, Series):
            return Series(np.dot(lvals, rvals), index=left.index)
        elif isinstance(rvals, np.ndarray):
            result = np.dot(lvals, rvals)
            if result.ndim == 2:
                return self._constructor(result, index=left.index)
            else:
                return Series(result, index=left.index)
        else:  # pragma: no cover
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
        index, columns = None, None
        orient = orient.lower()
        if orient == 'index':
            if len(data) > 0:
                # TODO speed up Series case
                if isinstance(data.values()[0], (Series, dict)):
                    data = _from_nested_dict(data)
                else:
                    data, index = data.values(), data.keys()
        elif orient != 'columns':  # pragma: no cover
            raise ValueError('only recognize index or columns for orient')

        return cls(data, index=index, columns=columns, dtype=dtype)

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
        import warnings
        if not self.columns.is_unique:
            warnings.warn("DataFrame columns are not unique, some "
                          "columns will be omitted.", UserWarning)
        if outtype.lower().startswith('d'):
            return dict((k, v.to_dict()) for k, v in self.iteritems())
        elif outtype.lower().startswith('l'):
            return dict((k, v.tolist()) for k, v in self.iteritems())
        elif outtype.lower().startswith('s'):
            return dict((k, v) for k, v in self.iteritems())
        else:  # pragma: no cover
            raise ValueError("outtype %s not understood" % outtype)

    @classmethod
    def from_records(cls, data, index=None, exclude=None, columns=None,
                     coerce_float=False, nrows=None):
        """
        Convert structured or record ndarray to DataFrame

        Parameters
        ----------
        data : ndarray (structured dtype), list of tuples, dict, or DataFrame
        index : string, list of fields, array-like
            Field of array to use as the index, alternately a specific set of
            input labels to use
        exclude: sequence, default None
            Columns or fields to exclude
        columns : sequence, default None
            Column names to use. If the passed data do not have named
            associated with them, this argument provides names for the
            columns. Otherwise this argument indicates the order of the columns
            in the result (any names not found in the data will become all-NA
            columns)
        coerce_float : boolean, default False
            Attempt to convert values to non-string, non-numeric objects (like
            decimal.Decimal) to floating point, useful for SQL result sets

        Returns
        -------
        df : DataFrame
        """
        # Make a copy of the input columns so we can modify it
        if columns is not None:
            columns = _ensure_index(columns)

        if com.is_iterator(data):
            if nrows == 0:
                return cls()

            try:
                if py3compat.PY3:
                    first_row = next(data)
                else:
                    first_row = data.next()
            except StopIteration:
                return cls(index=index, columns=columns)

            dtype = None
            if hasattr(first_row, 'dtype') and first_row.dtype.names:
                dtype = first_row.dtype

            values = [first_row]

            i = 1
            for row in data:
                values.append(row)
                i += 1
                if i >= nrows:
                    break

            if dtype is not None:
                data = np.array(values, dtype=dtype)
            else:
                data = values

        if isinstance(data, dict):
            if columns is None:
                columns = arr_columns = _ensure_index(sorted(data))
                arrays = [data[k] for k in columns]
            else:
                arrays = []
                arr_columns = []
                for k, v in data.iteritems():
                    if k in columns:
                        arr_columns.append(k)
                        arrays.append(v)

                # reorder according to the columns
                if len(columns) and len(arr_columns):
                    indexer     = _ensure_index(arr_columns).get_indexer(columns)
                    arr_columns = _ensure_index([ arr_columns[i] for i in indexer ])
                    arrays      = [ arrays[i] for i in indexer ]

        elif isinstance(data, (np.ndarray, DataFrame)):
            arrays, columns = _to_arrays(data, columns)
            if columns is not None:
                columns = _ensure_index(columns)
            arr_columns = columns
        else:
            arrays, arr_columns = _to_arrays(data, columns,
                                             coerce_float=coerce_float)

            arr_columns = _ensure_index(arr_columns)
            if columns is not None:
                columns = _ensure_index(columns)
            else:
                columns = arr_columns

        if exclude is None:
            exclude = set()
        else:
            exclude = set(exclude)

        result_index = None
        if index is not None:
            if (isinstance(index, basestring) or
                    not hasattr(index, "__iter__")):
                i = columns.get_loc(index)
                exclude.add(index)
                if len(arrays) > 0:
                    result_index = Index(arrays[i], name=index)
                else:
                    result_index = Index([], name=index)
            else:
                try:
                    to_remove = [arr_columns.get_loc(field) for field in index]

                    result_index = MultiIndex.from_arrays(
                        [arrays[i] for i in to_remove], names=index)

                    exclude.update(index)
                except Exception:
                    result_index = index

        if any(exclude):
            arr_exclude = [x for x in exclude if x in arr_columns]
            to_remove = [arr_columns.get_loc(col) for col in arr_exclude]
            arrays = [v for i, v in enumerate(arrays) if i not in to_remove]

            arr_columns = arr_columns.drop(arr_exclude)
            columns = columns.drop(exclude)

        mgr = _arrays_to_mgr(arrays, arr_columns, result_index,
                             columns)

        return cls(mgr)

    def to_records(self, index=True, convert_datetime64=True):
        """
        Convert DataFrame to record array. Index will be put in the
        'index' field of the record array if requested

        Parameters
        ----------
        index : boolean, default True
            Include index in resulting record array, stored in 'index' field
        convert_datetime64 : boolean, default True
            Whether to convert the index to datetime.datetime if it is a
            DatetimeIndex

        Returns
        -------
        y : recarray
        """
        if index:
            if com.is_datetime64_dtype(self.index) and convert_datetime64:
                ix_vals = [self.index.to_pydatetime()]
            else:
                if isinstance(self.index, MultiIndex):
                    # array of tuples to numpy cols. copy copy copy
                    ix_vals = map(np.array,zip(*self.index.values))
                else:
                    ix_vals = [self.index.values]

            arrays = ix_vals+ [self[c].values for c in self.columns]

            count = 0
            index_names = self.index.names
            if isinstance(self.index, MultiIndex):
                for i, n in enumerate(index_names):
                    if n is None:
                        index_names[i] = 'level_%d' % count
                        count += 1
            elif index_names[0] is None:
                index_names = ['index']
            names = index_names + list(map(str, self.columns))
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
        orientation). The values should be arrays or Series.

        Parameters
        ----------
        items : sequence of (key, value) pairs
            Values should be arrays or Series.
        columns : sequence of column labels, optional
            Must be passed if orient='index'.
        orient : {'columns', 'index'}, default 'columns'
            The "orientation" of the data. If the keys of the
            input correspond to column labels, pass 'columns'
            (default). Otherwise if the keys correspond to the index,
            pass 'index'.

        Returns
        -------
        frame : DataFrame
        """
        keys, values = zip(*items)

        if orient == 'columns':
            if columns is not None:
                columns = _ensure_index(columns)

                idict = dict(items)
                if len(idict) < len(items):
                    if not columns.equals(_ensure_index(keys)):
                        raise ValueError('With non-unique item names, passed '
                                         'columns must be identical')
                    arrays = values
                else:
                    arrays = [idict[k] for k in columns if k in idict]
            else:
                columns = _ensure_index(keys)
                arrays = values

            return cls._from_arrays(arrays, columns, None)
        elif orient == 'index':
            if columns is None:
                raise ValueError("Must pass columns with orient='index'")

            keys = _ensure_index(keys)

            arr = np.array(values, dtype=object).T
            data = [lib.maybe_convert_objects(v) for v in arr]
            return cls._from_arrays(data, columns, keys)
        else:  # pragma: no cover
            raise ValueError("'orient' must be either 'columns' or 'index'")

    @classmethod
    def _from_arrays(cls, arrays, columns, index, dtype=None):
        mgr = _arrays_to_mgr(arrays, columns, index, columns, dtype=dtype)
        return cls(mgr)

    @classmethod
    def from_csv(cls, path, header=0, sep=',', index_col=0,
                 parse_dates=True, encoding=None, tupleize_cols=False):
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
        tupleize_cols : boolean, default True
            write multi_index columns as a list of tuples (if True)
            or new (expanded format) if False)

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
                          encoding=encoding,tupleize_cols=False)

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
        from pandas.core.reshape import block2d_to_blocknd

        # only support this kind for now
        if (not isinstance(self.index, MultiIndex) or
                len(self.index.levels) != 2):
            raise AssertionError('Must have 2-level MultiIndex')

        if not self.index.is_unique:
            raise Exception("Can't convert non-uniquely indexed "
                            "DataFrame to Panel")

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
            newb = block2d_to_blocknd(block.values.T, block.items, shape,
                                      [ major_labels, minor_labels ],
                                      ref_items=selfsorted.columns)
            new_blocks.append(newb)

        # preserve names, if any
        major_axis = major_axis.copy()
        major_axis.name = self.index.names[0]

        minor_axis = minor_axis.copy()
        minor_axis.name = self.index.names[1]

        new_axes = [selfsorted.columns, major_axis, minor_axis]
        new_mgr = create_block_manager_from_blocks(new_blocks, new_axes)

        return Panel(new_mgr)

    to_wide = deprecate('to_wide', to_panel)

    def to_csv(self, path_or_buf, sep=",", na_rep='', float_format=None,
               cols=None, header=True, index=True, index_label=None,
               mode='w', nanRep=None, encoding=None, quoting=None,
               line_terminator='\n', chunksize=None,
               tupleize_cols=True, **kwds):
        r"""Write DataFrame to a comma-separated values (csv) file

        Parameters
        ----------
        path_or_buf : string or file handle / StringIO
            File path
        sep : character, default ","
            Field delimiter for the output file.
        na_rep : string, default ''
            Missing data representation
        float_format : string, default None
            Format string for floating point numbers
        cols : sequence, optional
            Columns to write
        header : boolean or list of string, default True
            Write out column names. If a list of string is given it is assumed
            to be aliases for the column names
        index : boolean, default True
            Write row names (index)
        index_label : string or sequence, or False, default None
            Column label for index column(s) if desired. If None is given, and
            `header` and `index` are True, then the index names are used. A
            sequence should be given if the DataFrame uses MultiIndex.  If
            False do not print fields for index names. Use index_label=False
            for easier importing in R
        nanRep : None
            deprecated, use na_rep
        mode : str
            Python write mode, default 'w'
        encoding : string, optional
            a string representing the encoding to use if the contents are
            non-ascii, for python versions prior to 3
        line_terminator : string, default '\\n'
            The newline character or character sequence to use in the output
            file
        quoting : optional constant from csv module
            defaults to csv.QUOTE_MINIMAL
        chunksize : int or None
            rows to write at a time
        tupleize_cols : boolean, default True
            write multi_index columns as a list of tuples (if True)
            or new (expanded format) if False)
        """
        if nanRep is not None:  # pragma: no cover
            import warnings
            warnings.warn("nanRep is deprecated, use na_rep",
                          FutureWarning)
            na_rep = nanRep

        formatter = fmt.CSVFormatter(self, path_or_buf,
                                     line_terminator=line_terminator,
                                     sep=sep, encoding=encoding,
                                     quoting=quoting,na_rep=na_rep,
                                     float_format=float_format, cols=cols,
                                     header=header, index=index,
                                     index_label=index_label,mode=mode,
                                     chunksize=chunksize,engine=kwds.get("engine"),
                                     tupleize_cols=tupleize_cols)
        formatter.save()

    def to_excel(self, excel_writer, sheet_name='sheet1', na_rep='',
                 float_format=None, cols=None, header=True, index=True,
                 index_label=None, startrow=0, startcol=0):
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
        startow : upper left cell row to dump data frame
        startcol : upper left cell column to dump data frame


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
        from pandas.io.excel import ExcelWriter
        need_save = False
        if isinstance(excel_writer, basestring):
            excel_writer = ExcelWriter(excel_writer)
            need_save = True

        formatter = fmt.ExcelFormatter(self,
                                       na_rep=na_rep,
                                       cols=cols,
                                       header=header,
                                       float_format=float_format,
                                       index=index,
                                       index_label=index_label)
        formatted_cells = formatter.get_formatted_cells()
        excel_writer.write_cells(formatted_cells, sheet_name,
                                 startrow=startrow, startcol=startcol)
        if need_save:
            excel_writer.save()

    def to_stata(self, fname, convert_dates=None, write_index=True, encoding="latin-1",
                 byteorder=None):
        """
        A class for writing Stata binary dta files from array-like objects

        Parameters
        ----------
        fname : file path or buffer
            Where to save the dta file.
        convert_dates : dict
            Dictionary mapping column of datetime types to the stata internal
            format that you want to use for the dates. Options are
            'tc', 'td', 'tm', 'tw', 'th', 'tq', 'ty'. Column can be either a
            number or a name.
        encoding : str
            Default is latin-1. Note that Stata does not support unicode.
        byteorder : str
            Can be ">", "<", "little", or "big". The default is None which uses
            `sys.byteorder`

        Examples
        --------
        >>> writer = StataWriter('./data_file.dta', data)
        >>> writer.write_file()

        Or with dates

        >>> writer = StataWriter('./date_data_file.dta', data, {2 : 'tw'})
        >>> writer.write_file()
        """
        from pandas.io.stata import StataWriter
        writer = StataWriter(fname,self,convert_dates=convert_dates, encoding=encoding, byteorder=byteorder)
        writer.write_file()

    def to_sql(self, name, con, flavor='sqlite', if_exists='fail', **kwargs):
        """
        Write records stored in a DataFrame to a SQL database.

        Parameters
        ----------
        name: name of SQL table
        conn: an open SQL database connection object
        flavor: {'sqlite', 'mysql', 'oracle'}, default 'sqlite'
        if_exists: {'fail', 'replace', 'append'}, default 'fail'
            fail: If table exists, do nothing.
            replace: If table exists, drop it, recreate it, and insert data.
            append: If table exists, insert data. Create if does not exist.
        """
        from pandas.io.sql import write_frame
        write_frame(self, name, con, flavor=flavor, if_exists=if_exists, **kwargs)

    @Appender(fmt.docstring_to_string, indents=1)
    def to_string(self, buf=None, columns=None, col_space=None, colSpace=None,
                  header=True, index=True, na_rep='NaN', formatters=None,
                  float_format=None, sparsify=None, nanRep=None,
                  index_names=True, justify=None, force_unicode=None,
                  line_width=None):
        """
        Render a DataFrame to a console-friendly tabular output.
        """
        import warnings
        if force_unicode is not None:  # pragma: no cover
            warnings.warn("force_unicode is deprecated, it will have no "
                          "effect", FutureWarning)

        if nanRep is not None:  # pragma: no cover
            warnings.warn("nanRep is deprecated, use na_rep",
                          FutureWarning)
            na_rep = nanRep

        if colSpace is not None:  # pragma: no cover
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
                                           line_width=line_width)
        formatter.to_string()

        if buf is None:
            result = formatter.buf.getvalue()
            return result

    @Appender(fmt.docstring_to_string, indents=1)
    def to_html(self, buf=None, columns=None, col_space=None, colSpace=None,
                header=True, index=True, na_rep='NaN', formatters=None,
                float_format=None, sparsify=None, index_names=True,
                justify=None, force_unicode=None, bold_rows=True,
                classes=None, escape=True):
        """
        to_html-specific options

        bold_rows : boolean, default True
            Make the row labels bold in the output
        classes : str or list or tuple, default None
            CSS class(es) to apply to the resulting html table
        escape : boolean, default True
            Convert the characters <, >, and & to HTML-safe sequences.

        Render a DataFrame as an HTML table.
        """

        import warnings
        if force_unicode is not None:  # pragma: no cover
            warnings.warn("force_unicode is deprecated, it will have no "
                          "effect", FutureWarning)

        if colSpace is not None:  # pragma: no cover
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
                                           bold_rows=bold_rows,
                                           escape=escape)
        formatter.to_html(classes=classes)

        if buf is None:
            return formatter.buf.getvalue()

    @Appender(fmt.docstring_to_string, indents=1)
    def to_latex(self, buf=None, columns=None, col_space=None, colSpace=None,
                 header=True, index=True, na_rep='NaN', formatters=None,
                 float_format=None, sparsify=None, index_names=True,
                 bold_rows=True, force_unicode=None):
        """
        to_latex-specific options
        bold_rows : boolean, default True
            Make the row labels bold in the output

        Render a DataFrame to a tabular environment table.
        You can splice this into a LaTeX document.
        """

        import warnings
        if force_unicode is not None:  # pragma: no cover
            warnings.warn("force_unicode is deprecated, it will have no "
                          "effect", FutureWarning)

        if colSpace is not None:  # pragma: no cover
            warnings.warn("colSpace is deprecated, use col_space",
                          FutureWarning)
            col_space = colSpace

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

    def info(self, verbose=True, buf=None, max_cols=None):
        """
        Concise summary of a DataFrame, used in __repr__ when very large.

        Parameters
        ----------
        verbose : boolean, default True
            If False, don't print column count summary
        buf : writable buffer, defaults to sys.stdout
        max_cols : int, default None
            Determines whether full summary or short summary is printed
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
        if max_cols is None:
            max_cols = get_option('display.max_info_columns',len(self.columns)+1)

        if verbose and len(self.columns) <= max_cols:
            lines.append('Data columns (total %d columns):' % len(self.columns))
            space = max([len(com.pprint_thing(k)) for k in self.columns]) + 4
            counts = self.count()
            if len(cols) != len(counts):
                raise AssertionError('Columns must equal counts')
            for col, count in counts.iteritems():
                col = com.pprint_thing(col)
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

    def convert_objects(self, convert_dates=True, convert_numeric=False, copy=True):
        """
        Attempt to infer better dtype for object columns

        Parameters
        ----------
        convert_dates : if True, attempt to soft convert_dates, if 'coerce', force conversion (and non-convertibles get NaT)
        convert_numeric : if True attempt to coerce to numerbers (including strings), non-convertibles get NaN
        copy : boolean, return a copy if True (True by default)

        Returns
        -------
        converted : DataFrame
        """
        return self._constructor(self._data.convert(convert_dates=convert_dates,
                                                    convert_numeric=convert_numeric,
                                                    copy=copy))

    #----------------------------------------------------------------------
    # properties for index and columns

    columns = lib.AxisProperty(0)
    index = lib.AxisProperty(1)

    def as_matrix(self, columns=None):
        """
        Convert the frame to its Numpy-array matrix representation. Columns
        are presented in sorted order unless a specific list of columns is
        provided.

        NOTE: the dtype will be a lower-common-denominator dtype (implicit upcasting)
              that is to say if the dtypes (even of numeric types) are mixed, the one that accomodates all will be chosen
              use this with care if you are not dealing with the blocks

              e.g. if the dtypes are float16,float32         -> float32
                                     float16,float32,float64 -> float64
                                     int32,uint8             -> int32

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

    def as_blocks(self, columns=None):
        """
        Convert the frame to a dict of dtype -> DataFrames that each has a homogeneous dtype.
        are presented in sorted order unless a specific list of columns is
        provided.

        NOTE: the dtypes of the blocks WILL BE PRESERVED HERE (unlike in as_matrix)

        Parameters
        ----------
        columns : array-like
            Specific column order

        Returns
        -------
        values : a list of DataFrames
        """
        self._consolidate_inplace()

        bd = dict()
        for b in self._data.blocks:
            b = b.reindex_items_from(columns or b.items)
            bd[str(b.dtype)] = DataFrame(BlockManager([ b ], [ b.items, self.index ]))
        return bd

    blocks = property(fget=as_blocks)

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
            likely_dtype, value = _infer_dtype_from_scalar(value)

            made_bigger = not np.array_equal(new_columns, self.columns)

            # how to make this logic simpler?
            if made_bigger:
                com._possibly_cast_item(result, col, likely_dtype)

            return result.set_value(index, col, value)

    def irow(self, i, copy=False):
        return self._ixs(i,axis=0)

    def icol(self, i):
        return self._ixs(i,axis=1)

    def _ixs(self, i, axis=0, copy=False):
        """
        i : int, slice, or sequence of integers
        axis : int
        """

        # irow
        if axis == 0:

            """
            Notes
            -----
            If slice passed, the resulting data will be a view
            """

            if isinstance(i, slice):
                return self[i]
            else:
                label = self.index[i]
                if isinstance(label, Index):

                    # a location index by definition
                    i = _maybe_convert_indices(i, len(self._get_axis(axis)))
                    return self.reindex(i, takeable=True)
                else:
                    try:
                        new_values = self._data.fast_2d_xs(i, copy=copy)
                    except:
                        new_values = self._data.fast_2d_xs(i, copy=True)
                    return Series(new_values, index=self.columns,
                                  name=self.index[i])

        # icol
        else:

            """
            Notes
            -----
            If slice passed, the resulting data will be a view
            """

            label = self.columns[i]
            if isinstance(i, slice):
                # need to return view
                lab_slice = slice(label[0], label[-1])
                return self.ix[:, lab_slice]
            else:
                label = self.columns[i]
                if isinstance(label, Index):
                    return self.take(i, axis=1, convert=True)

                values = self._data.iget(i)
                return self._col_klass.from_array(values, index=self.index,
                                                  name=label)

    def iget_value(self, i, j):
        return self.iat[i,j]

    def __getitem__(self, key):

        # see if we can slice the rows
        indexer = _convert_to_index_sliceable(self, key)
        if indexer is not None:
            return self._getitem_slice(indexer)

        if isinstance(key, (np.ndarray, list)):
            # either boolean or fancy integer index
            return self._getitem_array(key)
        elif isinstance(key, DataFrame):
            return self._getitem_frame(key)
        elif isinstance(self.columns, MultiIndex):
            return self._getitem_multilevel(key)
        else:
            # get column
            if self.columns.is_unique:
                return self._get_item_cache(key)

            # duplicate columns
            return self._constructor(self._data.get(key))

    def _getitem_slice(self, key):
        return self._slice(key, axis=0)

    def _getitem_array(self, key):
        # also raises Exception if object array with NA values
        if com._is_bool_indexer(key):
            # warning here just in case -- previously __setitem__ was
            # reindexing but __getitem__ was not; it seems more reasonable to
            # go with the __setitem__ behavior since that is more consistent
            # with all other indexing behavior
            if isinstance(key, Series) and not key.index.equals(self.index):
                import warnings
                warnings.warn("Boolean Series key will be reindexed to match "
                              "DataFrame index.", UserWarning)
            elif len(key) != len(self.index):
                raise ValueError('Item wrong length %d instead of %d!' %
                                 (len(key), len(self.index)))
            # _check_bool_indexer will throw exception if Series key cannot
            # be reindexed to match DataFrame rows
            key = _check_bool_indexer(self.index, key)
            indexer = key.nonzero()[0]
            return self.take(indexer, axis=0, convert=False)
        else:
            indexer = self.ix._convert_to_indexer(key, axis=1)
            return self.take(indexer, axis=1, convert=True)

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
                if ((type(top) == str and top == '') or
                        (type(top) == tuple and top[0] == '')):
                    result = result['']
                    if isinstance(result, Series):
                        result = Series(result, index=self.index, name=key)

            return result
        else:
            return self._get_item_cache(key)

    def _getitem_frame(self, key):
        if key.values.dtype != np.bool_:
            raise ValueError('Must pass DataFrame with boolean values only')
        return self.where(key)

    def _slice(self, slobj, axis=0, raise_on_error=False):
        if axis == 0:
            mgr_axis = 1
        else:
            mgr_axis = 0

        self._consolidate_inplace()
        new_data = self._data.get_slice(slobj, axis=mgr_axis,
                                        raise_on_error=raise_on_error)

        return self._constructor(new_data)

    def _box_item_values(self, key, values):
        items = self.columns[self.columns.get_loc(key)]
        if values.ndim == 2:
            return self._constructor(values.T, columns=items, index=self.index)
        else:
            return Series.from_array(values, index=self.index, name=items)

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
        # see if we can slice the rows
        indexer = _convert_to_index_sliceable(self, key)
        if indexer is not None:
            return self._setitem_slice(indexer, value)

        if isinstance(key, (np.ndarray, list)):
            self._setitem_array(key, value)
        elif isinstance(key, DataFrame):
            self._setitem_frame(key, value)
        else:
            # set column
            self._set_item(key, value)

    def _setitem_slice(self, key, value):
        self.ix._setitem_with_indexer(key, value)

    def _setitem_array(self, key, value):
        # also raises Exception if object array with NA values
        if com._is_bool_indexer(key):
            if len(key) != len(self.index):
                raise ValueError('Item wrong length %d instead of %d!' %
                                 (len(key), len(self.index)))
            key = _check_bool_indexer(self.index, key)
            indexer = key.nonzero()[0]
            self.ix._setitem_with_indexer(indexer, value)
        else:
            if isinstance(value, DataFrame):
                if len(value.columns) != len(key):
                    raise AssertionError('Columns must be same length as key')
                for k1, k2 in zip(key, value.columns):
                    self[k1] = value[k2]
            else:
                indexer = self.ix._convert_to_indexer(key, axis=1)
                self.ix._setitem_with_indexer((slice(None), indexer), value)

    def _setitem_frame(self, key, value):
        # support boolean setting with DataFrame input, e.g.
        # df[df > df2] = 0
        if key.values.dtype != np.bool_:
            raise ValueError('Must pass DataFrame with boolean values only')

        if self._is_mixed_type:
            if not self._is_numeric_mixed_type:
                raise ValueError('Cannot do boolean setting on mixed-type frame')

        self.where(-key, value, inplace=True)

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

    def insert(self, loc, column, value, allow_duplicates=False):
        """
        Insert column into DataFrame at specified location.
        if allow_duplicates is False, Raises Exception if column is already contained in the DataFrame

        Parameters
        ----------
        loc : int
            Must have 0 <= loc <= len(columns)
        column : object
        value : int, Series, or array-like
        """
        value = self._sanitize_column(column, value)
        self._data.insert(loc, column, value, allow_duplicates=allow_duplicates)

    def _sanitize_column(self, key, value):
        # Need to make sure new columns (which go into the BlockManager as new
        # blocks) are always copied
        if _is_sequence(value):
            is_frame = isinstance(value, DataFrame)
            if isinstance(value, Series) or is_frame:
                if value.index.equals(self.index):
                    # copy the values
                    value = value.values.copy()
                else:

                    # GH 4107
                    try:
                        value = value.reindex(self.index).values
                    except:
                        raise TypeError('incompatible index of inserted column '
                                        'with frame index')

                if is_frame:
                    value = value.T
            else:
                if len(value) != len(self.index):
                    raise AssertionError('Length of values does not match '
                                         'length of index')

                if not isinstance(value, np.ndarray):
                    if isinstance(value, list) and len(value) > 0:
                        value = com._possibly_convert_platform(value)
                    else:
                        value = com._asarray_tuplesafe(value)
                elif isinstance(value, PeriodIndex):
                    value = value.asobject
                elif value.ndim == 2:
                    value = value.copy().T
                else:
                    value = value.copy()

            # Broadcasting funtimes
            if key in self.columns and value.ndim == 1:
                existing_piece = self[key]
                if isinstance(existing_piece, DataFrame):
                    value = np.tile(value, (len(existing_piece.columns), 1))
        else:
            if key in self.columns:
                existing_piece = self[key]

                # upcast the scalar
                dtype, value = _infer_dtype_from_scalar(value)

                # transpose hack
                if isinstance(existing_piece, DataFrame):
                    shape = (len(existing_piece.columns), len(self.index))
                    value = np.repeat(value, np.prod(shape)).reshape(shape)
                else:
                    value = np.repeat(value, len(self.index))

                value = value.astype(dtype)

            else:
                # upcast the scalar
                dtype, value = _infer_dtype_from_scalar(value)
                value = np.array(np.repeat(value, len(self.index)), dtype=dtype)

            value = com._possibly_cast_to_datetime(value, dtype)
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
        axis = self._get_axis_number(axis)
        labels = self._get_axis(axis)
        if level is not None:
            loc, new_ax = labels.get_loc_level(key, level=level)

            if not copy and not isinstance(loc, slice):
                raise ValueError('Cannot retrieve view (copy=False)')

            # level = 0
            loc_is_slice = isinstance(loc, slice)
            if not loc_is_slice:
                indexer = [slice(None)] * 2
                indexer[axis] = loc
                indexer = tuple(indexer)
            else:
                indexer = loc
                lev_num = labels._get_level_number(level)
                if labels.levels[lev_num].inferred_type == 'integer':
                    indexer = self.index[loc]

            # select on the correct axis
            if axis == 1 and loc_is_slice:
                indexer = slice(None), indexer
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
                    return self.take(inds, axis=axis, convert=False)
                else:
                    return self.take(loc, axis=axis, convert=True)

            if not np.isscalar(loc):
                new_index = self.index[loc]

        if np.isscalar(loc):
            new_values = self._data.fast_2d_xs(loc, copy=copy)
            return Series(new_values, index=self.columns,
                          name=self.index[loc])
        else:
            result = self[loc]
            result.index = new_index
            return result

    _xs = xs

    def lookup(self, row_labels, col_labels):
        """Label-based "fancy indexing" function for DataFrame. Given
        equal-length arrays of row and column labels, return an array of the
        values corresponding to each (row, col)  pair.

        Parameters
        ----------
        row_labels : sequence
            The row labels to use for lookup
        col_labels : sequence
            The column labels to use for lookup

        Notes
        -----
        Akin to

            .. code-block:: python

                result = []
                for row, col in zip(row_labels, col_labels):
                    result.append(df.get_value(row, col))

        Examples
        --------
        values : ndarray
            The found values

        """
        from itertools import izip

        n = len(row_labels)
        if n != len(col_labels):
            raise AssertionError('Row labels must have same size as '
                                 'column labels')

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
              fill_value=NA, method=None, limit=None, fill_axis=0):
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
        if axis is not None:
            axis = self._get_axis_number(axis)
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
                     copy=True, fill_value=NA, method=None, limit=None,
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
                    right_result.fillna(fill_value, method=method,
                                        limit=limit))
        else:
            return left_result, right_result

    def reindex(self, index=None, columns=None, method=None, level=None,
                fill_value=NA, limit=None, copy=True, takeable=False):
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
        takeable : the labels are locations (and not labels)

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
                                           fill_value, limit, takeable)

        if index is not None:
            frame = frame._reindex_index(index, method, copy, level,
                                         fill_value, limit, takeable)

        return frame

    def reindex_axis(self, labels, axis=0, method=None, level=None, copy=True,
                     limit=None, fill_value=NA):
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
        axis = self._get_axis_number(axis)
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
            indexer = row_indexer, col_indexer
            new_values = com.take_2d_multi(self.values, indexer,
                                           fill_value=fill_value)
            return self._constructor(new_values, index=new_index,
                                     columns=new_columns)
        elif row_indexer is not None:
            return self._reindex_with_indexers(new_index, row_indexer,
                                               None, None, copy, fill_value)
        elif col_indexer is not None:
            return self._reindex_with_indexers(None, None,
                                               new_columns, col_indexer,
                                               copy, fill_value)
        else:
            return self.copy() if copy else self

    def _reindex_index(self, new_index, method, copy, level, fill_value=NA,
                       limit=None, takeable=False):
        new_index, indexer = self.index.reindex(new_index, method, level,
                                                limit=limit, copy_if_needed=True,
                                                takeable=takeable)
        return self._reindex_with_indexers(new_index, indexer, None, None,
                                           copy, fill_value)

    def _reindex_columns(self, new_columns, copy, level, fill_value=NA,
                         limit=None, takeable=False):
        new_columns, indexer = self.columns.reindex(new_columns, level=level,
                                                    limit=limit, copy_if_needed=True,
                                                    takeable=takeable)
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

    def reindex_like(self, other, method=None, copy=True, limit=None,
                     fill_value=NA):
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
                            method=method, copy=copy, limit=limit,
                            fill_value=fill_value)

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
        >>> indexed_df = df.set_index(['A', 'B'])
        >>> indexed_df2 = df.set_index(['A', [0, 1, 2, 0, 1, 2]])
        >>> indexed_df3 = df.set_index([[0, 1, 2, 0, 1, 2]])

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

        if not inplace:
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
            If the columns have multiple levels, determines how the other
            levels are named. If None then the index name is repeated.

        Returns
        -------
        resetted : DataFrame
        """
        if inplace:
            new_obj = self
        else:
            new_obj = self.copy()

        def _maybe_cast(values, labels=None):

            if values.dtype == np.object_:
                values = lib.maybe_convert_objects(values)

            # if we have the labels, extract the values with a mask
            if labels is not None:
                mask = labels == -1
                values = values.take(labels)
                if mask.any():
                    values, changed = com._maybe_upcast_putmask(values,mask,np.nan)

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
                    level_values = _maybe_cast(lev.values, lab)
                    if level is None or i in level:
                        new_obj.insert(0, col_name, level_values)

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
            if isinstance(self.index, PeriodIndex):
                values = self.index.asobject
            elif (isinstance(self.index, DatetimeIndex) and
                  self.index.tz is not None):
                values = self.index.asobject
            else:
                values = _maybe_cast(self.index.values)
            new_obj.insert(0, name, values)

        new_obj.index = new_index
        if not inplace:
            return new_obj

    delevel = deprecate('delevel', reset_index)

    def take(self, indices, axis=0, convert=True):
        """
        Analogous to ndarray.take, return DataFrame corresponding to requested
        indices along an axis

        Parameters
        ----------
        indices : list / array of ints
        axis : {0, 1}
        convert : convert indices for negative values, check bounds, default True
                  mainly useful for an user routine calling

        Returns
        -------
        taken : DataFrame
        """

        # check/convert indicies here
        if convert:
            axis = self._get_axis_number(axis)
            indices = _maybe_convert_indices(indices, len(self._get_axis(axis)))

        if self._is_mixed_type:
            if axis == 0:
                new_data = self._data.take(indices, axis=1, verify=False)
                return DataFrame(new_data)
            else:
                new_columns = self.columns.take(indices)
                return self.reindex(columns=new_columns)
        else:
            new_values = com.take_nd(self.values,
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
            matchf = lambda x: (like in x if isinstance(x, basestring)
                                else like in str(x))
            return self.select(matchf, axis=1)
        elif regex:
            matcher = re.compile(regex)
            return self.select(lambda x: matcher.search(x) is not None, axis=1)
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

        axis = self._get_axis_number(axis)
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

        return self.take(mask.nonzero()[0], axis=axis, convert=False)

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
        ascending : boolean or list, default True
            Sort ascending vs. descending. Specify list for multiple sort
            orders
        axis : {0, 1}
            Sort index/rows versus columns
        inplace : boolean, default False
            Sort the DataFrame without creating a new instance

        Examples
        --------
        >>> result = df.sort(['A', 'B'], ascending=[1, 0])

        Returns
        -------
        sorted : DataFrame
        """
        if column is not None:  # pragma: no cover
            import warnings
            warnings.warn("column is deprecated, use columns", FutureWarning)
            columns = column
        return self.sort_index(by=columns, axis=axis, ascending=ascending,
                               inplace=inplace)

    def sort_index(self, axis=0, by=None, ascending=True, inplace=False,
                   kind='quicksort'):
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
        ascending : boolean or list, default True
            Sort ascending vs. descending. Specify list for multiple sort
            orders
        inplace : boolean, default False
            Sort the DataFrame without creating a new instance

        Examples
        --------
        >>> result = df.sort_index(by=['A', 'B'], ascending=[1, 0])

        Returns
        -------
        sorted : DataFrame
        """
        from pandas.core.groupby import _lexsort_indexer

        axis = self._get_axis_number(axis)
        if axis not in [0, 1]:
            raise ValueError('Axis must be 0 or 1, got %s' % str(axis))

        labels = self._get_axis(axis)

        if by is not None:
            if axis != 0:
                raise AssertionError('Axis must be 0')
            if not isinstance(by, (tuple, list)):
                by = [by]

            if len(by) > 1:
                keys = []
                for x in by:
                    k = self[x].values
                    if k.ndim == 2:
                        raise ValueError('Cannot sort by duplicate column %s'
                                         % str(x))
                    keys.append(k)

                def trans(v):
                    if com.needs_i8_conversion(v):
                        return v.view('i8')
                    return v

                keys = [trans(self[x].values) for x in by]
                indexer = _lexsort_indexer(keys, orders=ascending)
                indexer = com._ensure_platform_int(indexer)
            else:
                by = by[0]
                k = self[by].values
                if k.ndim == 2:
                    raise ValueError('Cannot sort by duplicate column %s'
                                     % str(by))
                indexer = k.argsort(kind=kind)
                if not ascending:
                    indexer = indexer[::-1]
        elif isinstance(labels, MultiIndex):
            indexer = _lexsort_indexer(labels.labels, orders=ascending)
            indexer = com._ensure_platform_int(indexer)
        else:
            indexer = labels.argsort(kind=kind)
            if not ascending:
                indexer = indexer[::-1]

        if inplace:
            if axis == 1:
                self._data = self._data.reindex_items(
                    self._data.items[indexer],
                    copy=False)
            elif axis == 0:
                self._data = self._data.take(indexer)

            self._clear_item_cache()
        else:
            return self.take(indexer, axis=axis, convert=False)

    def sortlevel(self, level=0, axis=0, ascending=True, inplace=False):
        """
        Sort multilevel index by chosen axis and primary level. Data will be
        lexicographically sorted by the chosen level followed by the other
        levels (in order)

        Parameters
        ----------
        level : int
        axis : {0, 1}
        ascending : bool, default True
        inplace : boolean, default False
            Sort the DataFrame without creating a new instance

        Returns
        -------
        sorted : DataFrame
        """
        axis = self._get_axis_number(axis)
        the_axis = self._get_axis(axis)
        if not isinstance(the_axis, MultiIndex):
            raise Exception('can only sort by level with a hierarchical index')

        new_axis, indexer = the_axis.sortlevel(level, ascending=ascending)

        if self._is_mixed_type and not inplace:
            if axis == 0:
                return self.reindex(index=new_axis)
            else:
                return self.reindex(columns=new_axis)

        if inplace:
            if axis == 1:
                self._data = self._data.reindex_items(
                    self._data.items[indexer],
                    copy=False)
            elif axis == 0:
                self._data = self._data.take(indexer)

            self._clear_item_cache()
        else:
            return self.take(indexer, axis=axis, convert=False)

    def swaplevel(self, i, j, axis=0):
        """
        Swap levels i and j in a MultiIndex on a particular axis

        Parameters
        ----------
        i, j : int, string (can be mixed)
            Level of index to be swapped. Can pass level name as string.

        Returns
        -------
        swapped : type of caller (new object)
        """
        result = self.copy()

        axis = self._get_axis_number(axis)
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
        axis = self._get_axis_number(axis)
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

    def fillna(self, value=None, method=None, axis=0, inplace=False,
               limit=None, downcast=None):
        """
        Fill NA/NaN values using the specified method

        Parameters
        ----------
        method : {'backfill', 'bfill', 'pad', 'ffill', None}, default None
            Method to use for filling holes in reindexed Series
            pad / ffill: propagate last valid observation forward to next valid
            backfill / bfill: use NEXT valid observation to fill gap
        value : scalar or dict
            Value to use to fill holes (e.g. 0), alternately a dict of values
            specifying which value to use for each column (columns not in the
            dict will not be filled). This value cannot be a list.
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
        downcast : dict, default is None, a dict of item->dtype of what to
            downcast if possible

        See also
        --------
        reindex, asfreq

        Returns
        -------
        filled : DataFrame
        """
        if isinstance(value, (list, tuple)):
            raise TypeError('"value" parameter must be a scalar or dict, but '
                            'you passed a "{0}"'.format(type(value).__name__))
        self._consolidate_inplace()

        axis = self._get_axis_number(axis)
        if value is None:
            if method is None:
                raise ValueError('must specify a fill method or value')
            if self._is_mixed_type and axis == 1:
                if inplace:
                    raise NotImplementedError()
                return self.T.fillna(method=method, limit=limit).T

            method = com._clean_fill_method(method)
            new_data = self._data.interpolate(method  = method,
                                              axis    = axis,
                                              limit   = limit,
                                              inplace = inplace,
                                              coerce  = True)
        else:
            if method is not None:
                raise ValueError('cannot specify both a fill method and value')
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
                new_data = self._data.fillna(value, inplace=inplace,
                                             downcast=downcast)

        if inplace:
            self._data = new_data
        else:
            return self._constructor(new_data)

    def ffill(self, axis=0, inplace=False, limit=None):
        return self.fillna(method='ffill', axis=axis, inplace=inplace,
                           limit=limit)

    def bfill(self, axis=0, inplace=False, limit=None):
        return self.fillna(method='bfill', axis=axis, inplace=inplace,
                           limit=limit)

    def replace(self, to_replace=None, value=None, inplace=False, limit=None,
                regex=False, method=None, axis=None):
        """
        Replace values given in 'to_replace' with 'value'.

        Parameters
        ----------
        to_replace : str, regex, list, dict, Series, numeric, or None

            * str or regex:

                - str: string exactly matching `to_replace` will be replaced
                  with `value`
                - regex: regexs matching `to_replace` will be replaced with
                  `value`

            * list of str, regex, or numeric:

                - First, if `to_replace` and `value` are both lists, they
                  **must** be the same length.
                - Second, if ``regex=True`` then all of the strings in **both**
                  lists will be interpreted as regexs otherwise they will match
                  directly. This doesn't matter much for `value` since there
                  are only a few possible substitution regexes you can use.
                - str and regex rules apply as above.

            * dict:

                - Nested dictionaries, e.g., {'a': {'b': nan}}, are read as
                  follows: look in column 'a' for the value 'b' and replace it
                  with nan. You can nest regular expressions as well. Note that
                  column names (the top-level dictionary keys in a nested
                  dictionary) **cannot** be regular expressions.
                - Keys map to column names and values map to substitution
                  values. You can treat this as a special case of passing two
                  lists except that you are specifying the column to search in.

            * None:

                - This means that the ``regex`` argument must be a string,
                  compiled regular expression, or list, dict, ndarray or Series
                  of such elements. If `value` is also ``None`` then this
                  **must** be a nested dictionary or ``Series``.

            See the examples section for examples of each of these.
        value : scalar, dict, list, str, regex, default None
            Value to use to fill holes (e.g. 0), alternately a dict of values
            specifying which value to use for each column (columns not in the
            dict will not be filled). Regular expressions, strings and lists or
            dicts of such objects are also allowed.
        inplace : boolean, default False
            If True, fill the DataFrame in place. Note: this will modify any
            other views on this DataFrame, like if you took a no-copy slice of
            an existing DataFrame, for example a column in a DataFrame. Returns
            a reference to the filled object, which is self if inplace=True
        limit : int, default None
            Maximum size gap to forward or backward fill
        regex : bool or same types as `to_replace`, default False
            Whether to interpret `to_replace` and/or `value` as regular
            expressions. If this is ``True`` then `to_replace` *must* be a
            string. Otherwise, `to_replace` must be ``None`` because this
            parameter will be interpreted as a regular expression or a list,
            dict, or array of regular expressions.

        See also
        --------
        reindex, asfreq, fillna

        Returns
        -------
        filled : DataFrame

        Raises
        ------
        AssertionError
            * If `regex` is not a ``bool`` and `to_replace` is not ``None``.
        TypeError
            * If `to_replace` is a ``dict`` and `value` is not a ``list``,
              ``dict``, ``ndarray``, or ``Series``
            * If `to_replace` is ``None`` and `regex` is not compilable into a
              regular expression or is a list, dict, ndarray, or Series.
        ValueError
            * If `to_replace` and `value` are ``list`` s or ``ndarray`` s, but
              they are not the same length.

        Notes
        -----
        * Regex substitution is performed under the hood with ``re.sub``. The
          rules for substitution for ``re.sub`` are the same.
        * Regular expressions will only substitute on strings, meaning you
          cannot provide, for example, a regular expression matching floating
          point numbers and expect the columns in your frame that have a
          numeric dtype to be matched. However, if those floating point numbers
          *are* strings, then you can do this.
        * This method has *a lot* of options. You are encouraged to experiment
          and play with this method to gain intuition about how it works.

        """
        if not com.is_bool(regex) and to_replace is not None:
            raise AssertionError("'to_replace' must be 'None' if 'regex' is "
                                 "not a bool")
        if method is not None:
            from warnings import warn
            warn('the "method" argument is deprecated and will be removed in'
                 'v0.13; this argument has no effect')

        if axis is not None:
            from warnings import warn
            warn('the "axis" argument is deprecated and will be removed in'
                 'v0.13; this argument has no effect')

        self._consolidate_inplace()

        if value is None:
            if not isinstance(to_replace, (dict, Series)):
                if not isinstance(regex, (dict, Series)):
                    raise TypeError('If "to_replace" and "value" are both None'
                                    ' then regex must be a mapping')
                to_replace = regex
                regex = True

            items = to_replace.items()
            keys, values = itertools.izip(*items)

            are_mappings = [isinstance(v, (dict, Series)) for v in values]

            if any(are_mappings):
                if not all(are_mappings):
                    raise TypeError("If a nested mapping is passed, all values"
                                    " of the top level mapping must be "
                                    "mappings")
                # passed a nested dict/Series
                to_rep_dict = {}
                value_dict = {}

                for k, v in items:
                    to_rep_dict[k] = v.keys()
                    value_dict[k] = v.values()

                to_replace, value = to_rep_dict, value_dict
            else:
                to_replace, value = keys, values

            return self.replace(to_replace, value, inplace=inplace,
                                limit=limit, regex=regex)
        else:
            if not len(self.columns):
                return self

            new_data = self._data
            if isinstance(to_replace, (dict, Series)):
                if isinstance(value, (dict, Series)):  # {'A' : NA} -> {'A' : 0}
                    new_data = self._data
                    for c, src in to_replace.iteritems():
                        if c in value and c in self:
                            new_data = new_data.replace(src, value[c],
                                                        filter=[c],
                                                        inplace=inplace,
                                                        regex=regex)

                elif not isinstance(value, (list, np.ndarray)):  # {'A': NA} -> 0
                    new_data = self._data
                    for k, src in to_replace.iteritems():
                        if k in self:
                            new_data = new_data.replace(src, value,
                                                        filter=[k],
                                                        inplace=inplace,
                                                        regex=regex)
                else:
                    raise TypeError('Fill value must be scalar, dict, or '
                                    'Series')

            elif isinstance(to_replace, (list, np.ndarray)):
                # [NA, ''] -> [0, 'missing']
                if isinstance(value, (list, np.ndarray)):
                    if len(to_replace) != len(value):
                        raise ValueError('Replacement lists must match '
                                         'in length. Expecting %d got %d ' %
                                         (len(to_replace), len(value)))

                    new_data = self._data.replace_list(to_replace, value,
                                                       inplace=inplace,
                                                       regex=regex)

                else:  # [NA, ''] -> 0
                    new_data = self._data.replace(to_replace, value,
                                                  inplace=inplace, regex=regex)
            elif to_replace is None:
                if not (com.is_re_compilable(regex) or
                        isinstance(regex, (list, dict, np.ndarray, Series))):
                    raise TypeError("'regex' must be a string or a compiled "
                                    "regular expression or a list or dict of "
                                    "strings or regular expressions, you "
                                    "passed a {0}".format(type(regex)))
                return self.replace(regex, value, inplace=inplace, limit=limit,
                                    regex=True)
            else:

                # dest iterable dict-like
                if isinstance(value, (dict, Series)):  # NA -> {'A' : 0, 'B' : -1}
                    new_data = self._data

                    for k, v in value.iteritems():
                        if k in self:
                            new_data = new_data.replace(to_replace, v,
                                                        filter=[k],
                                                        inplace=inplace,
                                                        regex=regex)

                elif not isinstance(value, (list, np.ndarray)):  # NA -> 0
                    new_data = self._data.replace(to_replace, value,
                                                  inplace=inplace, regex=regex)
                else:
                    raise TypeError('Invalid "to_replace" type: '
                                    '{0}'.format(type(to_replace)))  # pragma: no cover

        new_data = new_data.convert(copy=not inplace, convert_numeric=False)

        if inplace:
            self._data = new_data
        else:
            return self._constructor(new_data)

    def interpolate(self, to_replace, method='pad', axis=0, inplace=False,
                    limit=None):
        """Interpolate values according to different methods.

        Parameters
        ----------
        to_replace : dict, Series
        method : str
        axis : int
        inplace : bool
        limit : int, default None

        Returns
        -------
        frame : interpolated

        See Also
        --------
        reindex, replace, fillna
        """
        from warnings import warn
        warn('DataFrame.interpolate will be removed in v0.13, please use '
             'either DataFrame.fillna or DataFrame.replace instead',
             FutureWarning)
        if self._is_mixed_type and axis == 1:
            return self.T.replace(to_replace, method=method, limit=limit).T

        method = com._clean_fill_method(method)

        if isinstance(to_replace, (dict, Series)):
            if axis == 0:
                return self.replace(to_replace, method=method, inplace=inplace,
                                    limit=limit, axis=axis)
            elif axis == 1:
                obj = self.T
                if inplace:
                    obj.replace(to_replace, method=method, limit=limit,
                                inplace=inplace, axis=0)
                    return obj.T
                return obj.replace(to_replace, method=method, limit=limit,
                                   inplace=inplace, axis=0).T
            else:
                raise ValueError('Invalid value for axis')
        else:
            new_data = self._data.interpolate(method=method, axis=axis,
                                              limit=limit, inplace=inplace,
                                              missing=to_replace, coerce=False)

            if inplace:
                self._data = new_data
            else:
                return self._constructor(new_data)

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

        if not inplace:
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
            # but cannot outer join in align if dups anyways?
            result = {}
            for col in this:
                result[col] = _arith_op(this[col].values, other[col].values)
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
            return self * NA

        if len(self) == 0:
            # Ambiguous case, use _series so works with DataFrame
            return self._constructor(data=self._series, index=self.index,
                                     columns=self.columns)

        # teeny hack because one does DataFrame + TimeSeries all the time
        if self.index.is_all_dates and other.index.is_all_dates:
            import warnings
            warnings.warn(("TimeSeries broadcasting along DataFrame index "
                           "by default is deprecated. Please use "
                           "DataFrame.<op> to explicitly broadcast arithmetic "
                           "operations along the index"),
                          FutureWarning)
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

        new_data = left._data.eval(func, right, axes = [left.columns, self.index])
        return self._constructor(new_data)

    def _combine_const(self, other, func, raise_on_error = True):
        if self.empty:
            return self

        new_data = self._data.eval(func, other, raise_on_error=raise_on_error)
        return self._constructor(new_data)

    def _compare_frame(self, other, func, str_rep):
        if not self._indexed_same(other):
            raise Exception('Can only compare identically-labeled '
                            'DataFrame objects')

        def _compare(a, b):
            return dict([ (col,func(a[col], b[col])) for col in a.columns ])
        new_data = expressions.evaluate(_compare, str_rep, self, other)

        return self._constructor(data=new_data, index=self.index,
                                 columns=self.columns, copy=False)

    def _flex_compare_frame(self, other, func, str_rep, level):
        if not self._indexed_same(other):
            self, other = self.align(other, 'outer', level=level)

        def _compare(a, b):
            return dict([ (col,func(a[col], b[col])) for col in a.columns ])
        new_data = expressions.evaluate(_compare, str_rep, self, other)

        return self._constructor(data=new_data, index=self.index,
                                 columns=self.columns, copy=False)

    def combine(self, other, func, fill_value=None, overwrite=True):
        """
        Add two DataFrame objects and do not propagate NaN values, so if for a
        (column, time) one frame is missing a value, it will default to the
        other frame's value (which might be NaN as well)

        Parameters
        ----------
        other : DataFrame
        func : function
        fill_value : scalar value
        overwrite : boolean, default True
            If True then overwrite values for common keys in the calling frame

        Returns
        -------
        result : DataFrame
        """

        other_idxlen = len(other.index)  # save for compare

        this, other = self.align(other, copy=False)
        new_index = this.index

        if other.empty and len(new_index) == len(self.index):
            return self.copy()

        if self.empty and len(other) == other_idxlen:
            return other.copy()

        # sorts if possible
        new_columns = this.columns.union(other.columns)
        do_fill = fill_value is not None

        result = {}
        for col in new_columns:
            series = this[col]
            otherSeries = other[col]

            this_dtype = series.dtype
            other_dtype = otherSeries.dtype

            this_mask = isnull(series)
            other_mask = isnull(otherSeries)

            # don't overwrite columns unecessarily
            # DO propogate if this column is not in the intersection
            if not overwrite and other_mask.all():
                result[col] = this[col].copy()
                continue

            if do_fill:
                series = series.copy()
                otherSeries = otherSeries.copy()
                series[this_mask] = fill_value
                otherSeries[other_mask] = fill_value

            # if we have different dtypes, possibily promote
            new_dtype = this_dtype
            if this_dtype != other_dtype:
                new_dtype = com._lcd_dtypes(this_dtype,other_dtype)
                series = series.astype(new_dtype)
                otherSeries = otherSeries.astype(new_dtype)

            # see if we need to be represented as i8 (datetimelike)
            # try to keep us at this dtype
            needs_i8_conversion = com.needs_i8_conversion(new_dtype)
            if needs_i8_conversion:
                this_dtype = new_dtype
                arr = func(series, otherSeries, True)
            else:
                arr = func(series, otherSeries)

            if do_fill:
                arr = com.ensure_float(arr)
                arr[this_mask & other_mask] = NA

            # try to downcast back to the original dtype
            if needs_i8_conversion:
                arr = com._possibly_cast_to_datetime(arr, this_dtype)
            else:
                arr = com._possibly_downcast_to_dtype(arr, this_dtype)

            result[col] = arr

        # convert_objects just in case
        return self._constructor(result,
                                 index=new_index,
                                 columns=new_columns).convert_objects(
            convert_dates=True,
            copy=False)

    def combine_first(self, other):
        """
        Combine two DataFrame objects and default to non-null values in frame
        calling the method. Result index columns will be the union of the
        respective indexes and columns

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
        def combiner(x, y, needs_i8_conversion=False):
            x_values = x.values if hasattr(x,'values') else x
            y_values = y.values if hasattr(y,'values') else y
            if needs_i8_conversion:
                mask = isnull(x)
                x_values = x_values.view('i8')
                y_values = y_values.view('i8')
            else:
                mask = isnull(x_values)

            return expressions.where(mask, y_values, x_values, raise_on_error=True)

        return self.combine(other, combiner, overwrite=False)

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

                    # don't overwrite columns unecessarily
                    if mask.all():
                        continue
                else:
                    mask = notnull(this)

            self[col] = expressions.where(mask, this, that, raise_on_error=True)

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
        new_data = self._data.diff(periods)
        return self._constructor(new_data)

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

        if offset is None:
            indexer = com._shift_indexer(len(self), periods)
            new_data = self._data.shift(indexer, periods)
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

        See also
        --------
        DataFrame.applymap: For elementwise operations

        Returns
        -------
        applied : Series or DataFrame
        """
        if len(self.columns) == 0 and len(self.index) == 0:
            return self

        axis = self._get_axis_number(axis)
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
                    # How to determine this better?
                    is_reduction = False
                    try:
                        is_reduction = not isinstance(f(_EMPTY_SERIES),
                                                      np.ndarray)
                    except Exception:
                        pass

                    if is_reduction:
                        return Series(NA, index=self._get_agg_axis(axis))
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

            if self._is_mixed_type:  # maybe a hack for now
                raise AssertionError('Must be mixed type DataFrame')
            values = self.values
            dummy = Series(NA, index=self._get_axis(axis),
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
            # so will work with MultiIndex
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
                        e.args = e.args + ('occurred at index %s' %
                                           com.pprint_thing(k),)
                except (NameError, UnboundLocalError):  # pragma: no cover
                    # no k defined yet
                    pass
                raise e


        if len(results) > 0 and _is_sequence(results[0]):
            if not isinstance(results[0], Series):
                index = res_columns
            else:
                index = None

            result = self._constructor(data=results, index=index)
            result.columns = res_index

            if axis == 1:
                result = result.T
            result = result.convert_objects(copy=False)

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
        else:  # pragma: no cover
            raise ValueError('Axis must be 0 or 1, got %s' % axis)

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

        # if we have a dtype == 'M8[ns]', provide boxed values
        def infer(x):
            if com.is_datetime64_dtype(x):
                x = lib.map_infer(x, lib.Timestamp)
            return lib.map_infer(x, func)
        return self.apply(infer)

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
            if other.name is None:
                raise AssertionError('Other Series must have a name')
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
              left_index=False, right_index=False, sort=False,
              suffixes=('_x', '_y'), copy=True):
        from pandas.tools.merge import merge
        return merge(self, right, how=how, on=on,
                     left_on=left_on, right_on=right_on,
                     left_index=left_index, right_index=right_index, sort=sort,
                     suffixes=suffixes, copy=copy)

    #----------------------------------------------------------------------
    # Statistical methods, etc.

    def corr(self, method='pearson', min_periods=1):
        """
        Compute pairwise correlation of columns, excluding NA/null values

        Parameters
        ----------
        method : {'pearson', 'kendall', 'spearman'}
            pearson : standard correlation coefficient
            kendall : Kendall Tau correlation coefficient
            spearman : Spearman rank correlation
        min_periods : int, optional
            Minimum number of observations required per pair of columns
            to have a valid result. Currently only available for pearson
            and spearman correlation

        Returns
        -------
        y : DataFrame
        """
        numeric_df = self._get_numeric_data()
        cols = numeric_df.columns
        mat = numeric_df.values

        if method == 'pearson':
            correl = _algos.nancorr(com._ensure_float64(mat),
                                    minp=min_periods)
        elif method == 'spearman':
            correl = _algos.nancorr_spearman(com._ensure_float64(mat),
                                             minp=min_periods)
        else:
            if min_periods is None:
                min_periods = 1
            mat = mat.T
            corrf = nanops.get_corr_func(method)
            K = len(cols)
            correl = np.empty((K, K), dtype=float)
            mask = np.isfinite(mat)
            for i, ac in enumerate(mat):
                for j, bc in enumerate(mat):
                    valid = mask[i] & mask[j]
                    if valid.sum() < min_periods:
                        c = NA
                    elif not valid.all():
                        c = corrf(ac[valid], bc[valid])
                    else:
                        c = corrf(ac, bc)
                    correl[i, j] = c
                    correl[j, i] = c

        return self._constructor(correl, index=cols, columns=cols)

    def cov(self, min_periods=None):
        """
        Compute pairwise covariance of columns, excluding NA/null values

        Parameters
        ----------
        min_periods : int, optional
            Minimum number of observations required per pair of columns
            to have a valid result.

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
            if min_periods is not None and min_periods > len(mat):
                baseCov = np.empty((mat.shape[1], mat.shape[1]))
                baseCov.fill(np.nan)
            else:
                baseCov = np.cov(mat.T)
            baseCov = baseCov.reshape((len(cols),len(cols)))
        else:
            baseCov = _algos.nancorr(com._ensure_float64(mat), cov=True,
                                     minp=min_periods)

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
        axis = self._get_axis_number(axis)
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

        lb = .5 * (1. - percentile_width / 100.)
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
        axis = self._get_axis_number(axis)
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
        Return whether all elements are True over requested axis.
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
        """
        Notes
        -----
        This method returns the minimum of the values in the DataFrame. If you
        want the *index* of the minimum, use ``DataFrame.idxmin``. This is the
        equivalent of the ``numpy.ndarray`` method ``argmin``.

        See Also
        --------
        DataFrame.idxmin
        Series.idxmin
        """
        if level is not None:
            return self._agg_by_level('min', axis=axis, level=level,
                                      skipna=skipna)
        return self._reduce(nanops.nanmin, axis=axis, skipna=skipna,
                            numeric_only=None)

    @Substitution(name='maximum', shortname='max', na_action=_doc_exclude_na,
                  extras='')
    @Appender(_stat_doc)
    def max(self, axis=0, skipna=True, level=None):
        """
        Notes
        -----
        This method returns the maximum of the values in the DataFrame. If you
        want the *index* of the maximum, use ``DataFrame.idxmax``. This is the
        equivalent of the ``numpy.ndarray`` method ``argmax``.

        See Also
        --------
        DataFrame.idxmax
        Series.idxmax
        """
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

        axis = self._get_axis_number(axis)
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
        axis = self._get_axis_number(axis)
        method = getattr(type(self), name)
        applyf = lambda x: method(x, axis=axis, skipna=skipna, **kwds)
        return grouped.aggregate(applyf)

    def _reduce(self, op, axis=0, skipna=True, numeric_only=None,
                filter_type=None, **kwds):
        axis = self._get_axis_number(axis)
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
                    raise NotImplementedError
                result = f(data.values)
                labels = data._get_agg_axis(axis)
        else:
            if numeric_only:
                if filter_type is None or filter_type == 'numeric':
                    data = self._get_numeric_data()
                elif filter_type == 'bool':
                    data = self._get_bool_data()
                else:
                    raise NotImplementedError
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
                # otherwise, accept it
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

        Notes
        -----
        This method is the DataFrame version of ``ndarray.argmin``.

        See Also
        --------
        Series.idxmin
        """
        axis = self._get_axis_number(axis)
        indices = nanops.nanargmin(self.values, axis=axis, skipna=skipna)
        index = self._get_axis(axis)
        result = [index[i] if i >= 0 else NA for i in indices]
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

        Notes
        -----
        This method is the DataFrame version of ``ndarray.argmax``.

        See Also
        --------
        Series.idxmax
        """
        axis = self._get_axis_number(axis)
        indices = nanops.nanargmax(self.values, axis=axis, skipna=skipna)
        index = self._get_axis(axis)
        result = [index[i] if i >= 0 else NA for i in indices]
        return Series(result, index=self._get_agg_axis(axis))

    def _get_agg_axis(self, axis_num):
        if axis_num == 0:
            return self.columns
        elif axis_num == 1:
            return self.index
        else:
            raise Exception('Must have 0<= axis <= 1')

    def _get_numeric_data(self):
        return self._constructor(self._data.get_numeric_data(), index=self.index, copy=False)

    def _get_bool_data(self):
        return self._constructor(self._data.get_bool_data(), index=self.index, copy=False)

    def quantile(self, q=0.5, axis=0, numeric_only=True):
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
                return NA
            else:
                return _quantile(arr, per)

        data = self._get_numeric_data() if numeric_only else self
        return data.apply(f, axis=axis)

    def clip(self, lower=None, upper=None):
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

        # GH 2747 (arguments were reversed)
        if lower is not None and upper is not None:
            lower, upper = min(lower,upper), max(lower,upper)

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
        na_option : {'keep', 'top', 'bottom'}
            keep: leave NA values where they are
            top: smallest rank if ascending
            bottom: smallest rank if descending
        ascending : boolean, default True
            False for ranks by high (1) to low (N)

        Returns
        -------
        ranks : DataFrame
        """
        axis = self._get_axis_number(axis)
        if numeric_only is None:
            try:
                ranks = algos.rank(self.values, axis=axis, method=method,
                                   ascending=ascending, na_option=na_option)
                return self._constructor(ranks, index=self.index,
                                         columns=self.columns)
            except TypeError:
                numeric_only = True

        if numeric_only:
            data = self._get_numeric_data()
        else:
            data = self
        ranks = algos.rank(data.values, axis=axis, method=method,
                           ascending=ascending, na_option=na_option)
        return self._constructor(ranks, index=data.index, columns=data.columns)

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

        axis = self._get_axis_number(axis)
        if axis == 0:
            new_data.set_axis(1, self.index.to_timestamp(freq=freq, how=how))
        elif axis == 1:
            new_data.set_axis(0, self.columns.to_timestamp(freq=freq, how=how))
        else:
            raise ValueError('Axis must be 0 or 1. Got %s' % str(axis))

        return self._constructor(new_data)

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

        axis = self._get_axis_number(axis)
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

        return self._constructor(new_data)

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

    def where(self, cond, other=NA, inplace=False, try_cast=False, raise_on_error=True):
        """
        Return a DataFrame with the same shape as self and whose corresponding
        entries are from self where cond is True and otherwise are from other.

        Parameters
        ----------
        cond : boolean DataFrame or array
        other : scalar or DataFrame
        inplace : boolean, default False
            Whether to perform the operation in place on the data
        try_cast : boolean, default False
            try to cast the result back to the input type (if possible),
        raise_on_error : boolean, default True
            Whether to raise on invalid data types (e.g. trying to where on
            strings)

        Returns
        -------
        wh : DataFrame
        """
        if isinstance(cond, DataFrame):
            # this already checks for index/column equality
            cond = cond.reindex(self.index, columns=self.columns)
        else:
            if not hasattr(cond, 'shape'):
                raise ValueError('where requires an ndarray like object for its '
                                 'condition')
            if cond.shape != self.shape:
                raise ValueError('Array conditional must be same shape as self')
            cond = self._constructor(cond, index=self.index,
                                     columns=self.columns)

        if inplace:
            cond = -(cond.fillna(True).astype(bool))
        else:
            cond = cond.fillna(False).astype(bool)

        if isinstance(other, DataFrame):
            _, other = self.align(other, join='left', fill_value=NA)
        elif isinstance(other,np.ndarray):
            if other.shape != self.shape:
                raise ValueError('other must be the same shape as self '
                                 'when an ndarray')
            other = self._constructor(other, self.index, self.columns)

        if inplace:
            # we may have different type blocks come out of putmask, so
            # reconstruct the block manager
            self._data = self._data.putmask(cond,other,inplace=True)

        else:
            new_data = self._data.where(other, cond,
                                        raise_on_error=raise_on_error,
                                        try_cast=try_cast)

            return self._constructor(new_data)

    def mask(self, cond):
        """
        Returns copy of self whose values are replaced with nan if the
        inverted condition is True

        Parameters
        ----------
        cond: boolean DataFrame or array

        Returns
        -------
        wh: DataFrame
        """
        return self.where(~cond, NA)

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
        raise AssertionError('Function must reduce')

    for i, left_bound in enumerate(bounds):
        if i == len(bounds) - 1:
            right_bound = N
        else:
            right_bound = bounds[i + 1]

        result[i] = f(values[left_bound:right_bound])

    return result


def factor_agg(factor, vec, func):
    """
    Aggregate array based on Categorical

    Parameters
    ----------
    factor : Categorical
        length n
    vec : sequence
        length n
    func : function
        1D array aggregation function

    Returns
    -------
    ndarray corresponding to factor levels

    See Also
    --------
    pandas.Categorical
    """
    indexer = np.argsort(factor.labels)
    unique_labels = np.arange(len(factor.levels))

    ordered_labels = factor.labels.take(indexer)
    ordered_vec = np.asarray(vec).take(indexer)
    bounds = ordered_labels.searchsorted(unique_labels)

    return group_agg(ordered_vec, bounds, func)


def _arrays_to_mgr(arrays, arr_names, index, columns, dtype=None):
    """
    Segregate Series based on type and coerce into matrices.
    Needs to handle a lot of exceptional cases.
    """
    # figure out the index, if necessary
    if index is None:
        index = extract_index(arrays)
    else:
        index = _ensure_index(index)

    # don't force copy because getting jammed in an ndarray anyway
    arrays = _homogenize(arrays, index, dtype)

    # from BlockManager perspective
    axes = [_ensure_index(columns), _ensure_index(index)]

    return create_block_manager_from_arrays(arrays, arr_names, axes)

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

        for v in data:
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
            raise ValueError('If using all scalar values, you must must pass'
                             ' an index')

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
        if len(values) == 0:
            return np.empty((0, 0), dtype=object)

        def convert(v):
            return com._possibly_convert_platform(v)

        # we could have a 1-dim or 2-dim list here
        # this is equiv of np.asarray, but does object conversion
        # and platform dtype preservation
        if com.is_list_like(values[0]) or hasattr(values[0],'len'):
            values = np.array([ convert(v) for v in values])
        else:
            values = convert(values)

    else:
        # drop subclass info, do not copy data
        values = np.asarray(values)
        if copy:
            values = values.copy()

    if values.ndim == 1:
        values = values.reshape((values.shape[0], 1))
    elif values.ndim != 2:
        raise ValueError('Must pass 2-d input')

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


def _to_arrays(data, columns, coerce_float=False, dtype=None):
    """
    Return list of arrays, columns
    """
    if isinstance(data, DataFrame):
        if columns is not None:
            arrays = [data.icol(i).values for i, col in enumerate(data.columns)
                      if col in columns]
        else:
            columns = data.columns
            arrays = [data.icol(i).values for i in range(len(columns))]

        return arrays, columns

    if not len(data):
        if isinstance(data, np.ndarray):
            columns = data.dtype.names
            if columns is not None:
                return [[]] * len(columns), columns
        return [], []  # columns if columns is not None else []
    if isinstance(data[0], (list, tuple)):
        return _list_to_arrays(data, columns, coerce_float=coerce_float,
                               dtype=dtype)
    elif isinstance(data[0], collections.Mapping):
        return _list_of_dict_to_arrays(data, columns,
                                       coerce_float=coerce_float,
                                       dtype=dtype)
    elif isinstance(data[0], Series):
        return _list_of_series_to_arrays(data, columns,
                                         coerce_float=coerce_float,
                                         dtype=dtype)
    elif isinstance(data, np.ndarray):
        columns = list(data.dtype.names)
        arrays = [data[k] for k in columns]
        return arrays, columns
    else:
        # last ditch effort
        data = map(tuple, data)
        return _list_to_arrays(data, columns,
                               coerce_float=coerce_float,
                               dtype=dtype)


def _list_to_arrays(data, columns, coerce_float=False, dtype=None):
    if len(data) > 0 and isinstance(data[0], tuple):
        content = list(lib.to_object_array_tuples(data).T)
    else:
        # list of lists
        content = list(lib.to_object_array(data).T)
    return _convert_object_array(content, columns, dtype=dtype,
                                 coerce_float=coerce_float)


def _list_of_series_to_arrays(data, columns, coerce_float=False, dtype=None):
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
        return _convert_object_array(content, columns, dtype=dtype,
                                     coerce_float=coerce_float)
    else:
        return values.T, columns


def _list_of_dict_to_arrays(data, columns, coerce_float=False, dtype=None):
    if columns is None:
        gen = (x.keys() for x in data)
        columns = lib.fast_unique_multiple_list_gen(gen)

    # assure that they are of the base dict class and not of derived
    # classes
    data = [(type(d) is dict) and d or dict(d)
            for d in data]

    content = list(lib.dicts_to_array(data, list(columns)).T)
    return _convert_object_array(content, columns, dtype=dtype,
                                 coerce_float=coerce_float)


def _convert_object_array(content, columns, coerce_float=False, dtype=None):
    if columns is None:
        columns = _default_index(len(content))
    else:
        if len(columns) != len(content):
            raise AssertionError('%d columns passed, passed data had %s '
                                 'columns' % (len(columns), len(content)))

    arrays = [lib.maybe_convert_objects(arr, try_float=coerce_float)
              if dtype != object and dtype != np.object else arr
              for arr in content]

    return arrays, columns


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


def _homogenize(data, index, dtype=None):
    from pandas.core.series import _sanitize_array

    if dtype is not None:
        dtype = np.dtype(dtype)

    oindex = None
    homogenized = []

    for v in data:
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
                    v = lib.fast_multiget(v, oindex, default=NA)
                else:
                    v = lib.map_infer(oindex, v.get)

            v = _sanitize_array(v, index, dtype=dtype, copy=False,
                                raise_cast_failure=False)

        homogenized.append(v)

    return homogenized

def _from_nested_dict(data):
    # TODO: this should be seriously cythonized
    new_data = OrderedDict()
    for index, s in data.iteritems():
        for col, v in s.iteritems():
            new_data[col] = new_data.get(col, OrderedDict())
            new_data[col][index] = v
    return new_data


def _put_str(s, space):
    return ('%s' % s)[:space].ljust(space)


def install_ipython_completers():  # pragma: no cover
    """Register the DataFrame type with IPython's tab completion machinery, so
    that it knows about accessing column names as attributes."""
    from IPython.utils.generics import complete_object

    @complete_object.when_type(DataFrame)
    def complete_dataframe(obj, prev_completions):
        return prev_completions + [c for c in obj.columns
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
	ax : matplotlib axis object, default None
    fontsize : int or string
	rot : int, default None
        Rotation for ticks
	grid : boolean, default None (matlab style default)
        Axis grid lines

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
