"""
DataFrame
---------
An efficient 2D container for potentially mixed-type time series or other
labeled data series.

Similar to its R counterpart, data.frame, except providing automatic data
alignment and a host of useful data manipulation methods having to do with the
labeling information
"""
from __future__ import division
# pylint: disable=E1101,E1103
# pylint: disable=W0212,W0231,W0703,W0622

import functools
import collections
import itertools
import sys
import types
import warnings

from numpy import nan as NA
import numpy as np
import numpy.ma as ma

from pandas.types.cast import (_maybe_upcast, _infer_dtype_from_scalar,
                               _possibly_cast_to_datetime,
                               _possibly_infer_to_datetimelike,
                               _possibly_convert_platform,
                               _possibly_downcast_to_dtype,
                               _invalidate_string_dtypes,
                               _coerce_to_dtypes,
                               _maybe_upcast_putmask,
                               _find_common_type)
from pandas.types.common import (is_categorical_dtype,
                                 is_object_dtype,
                                 is_extension_type,
                                 is_datetimetz,
                                 is_datetime64_any_dtype,
                                 is_datetime64tz_dtype,
                                 is_bool_dtype,
                                 is_integer_dtype,
                                 is_float_dtype,
                                 is_integer,
                                 is_scalar,
                                 is_dtype_equal,
                                 needs_i8_conversion,
                                 _get_dtype_from_object,
                                 _ensure_float,
                                 _ensure_float64,
                                 _ensure_int64,
                                 _ensure_platform_int,
                                 is_list_like,
                                 is_iterator,
                                 is_sequence,
                                 is_named_tuple)
from pandas.types.missing import isnull, notnull

from pandas.core.common import (PandasError, _try_sort,
                                _default_index,
                                _values_from_object,
                                _maybe_box_datetimelike,
                                _dict_compat)
from pandas.core.generic import NDFrame, _shared_docs
from pandas.core.index import Index, MultiIndex, _ensure_index
from pandas.core.indexing import (maybe_droplevels, convert_to_index_sliceable,
                                  check_bool_indexer)
from pandas.core.internals import (BlockManager,
                                   create_block_manager_from_arrays,
                                   create_block_manager_from_blocks)
from pandas.core.series import Series
from pandas.core.categorical import Categorical
import pandas.computation.expressions as expressions
import pandas.core.algorithms as algorithms
from pandas.computation.eval import eval as _eval
from pandas.compat import (range, map, zip, lrange, lmap, lzip, StringIO, u,
                           OrderedDict, raise_with_traceback)
from pandas import compat
from pandas.compat.numpy import function as nv
from pandas.util.decorators import (deprecate_kwarg, Appender,
                                    Substitution)
from pandas.util.validators import validate_bool_kwarg

from pandas.tseries.period import PeriodIndex
from pandas.tseries.index import DatetimeIndex
from pandas.tseries.tdi import TimedeltaIndex

import pandas.core.base as base
import pandas.core.common as com
import pandas.core.nanops as nanops
import pandas.core.ops as ops
import pandas.formats.format as fmt
from pandas.formats.printing import pprint_thing
import pandas.tools.plotting as gfx

from pandas._libs import lib, algos as libalgos

from pandas.core.config import get_option

# ---------------------------------------------------------------------
# Docstring templates

_shared_doc_kwargs = dict(
    axes='index, columns', klass='DataFrame',
    axes_single_arg="{0 or 'index', 1 or 'columns'}",
    optional_by="""
        by : str or list of str
            Name or list of names which refer to the axis items.""",
    versionadded_to_excel='')

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
indicator : boolean or string, default False
    If True, adds a column to output DataFrame called "_merge" with
    information on the source of each row.
    If string, column with information on source of each row will be added to
    output DataFrame, and column will be named value of string.
    Information column is Categorical-type and takes on a value of "left_only"
    for observations whose merge key only appears in 'left' DataFrame,
    "right_only" for observations whose merge key only appears in 'right'
    DataFrame, and "both" if the observation's merge key is found in both.

    .. versionadded:: 0.17.0

Examples
--------

>>> A              >>> B
    lkey value         rkey value
0   foo  1         0   foo  5
1   bar  2         1   bar  6
2   baz  3         2   qux  7
3   foo  4         3   bar  8

>>> A.merge(B, left_on='lkey', right_on='rkey', how='outer')
   lkey  value_x  rkey  value_y
0  foo   1        foo   5
1  foo   4        foo   5
2  bar   2        bar   6
3  bar   2        bar   8
4  baz   3        NaN   NaN
5  NaN   NaN      qux   7

Returns
-------
merged : DataFrame
    The output type will the be same as 'left', if it is a subclass
    of DataFrame.

See also
--------
merge_ordered
merge_asof

"""

# -----------------------------------------------------------------------
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
        Column labels to use for resulting frame. Will default to
        np.arange(n) if no column labels are provided
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
    DataFrame.from_records : constructor from tuples, also record arrays
    DataFrame.from_dict : from dicts of Series, arrays, or dicts
    DataFrame.from_items : from sequence of (key, value) pairs
    pandas.read_csv, pandas.read_table, pandas.read_clipboard
    """

    @property
    def _constructor(self):
        return DataFrame

    _constructor_sliced = Series

    @property
    def _constructor_expanddim(self):
        from pandas.core.panel import Panel
        return Panel

    def __init__(self, data=None, index=None, columns=None, dtype=None,
                 copy=False):
        if data is None:
            data = {}
        if dtype is not None:
            dtype = self._validate_dtype(dtype)

        if isinstance(data, DataFrame):
            data = data._data

        if isinstance(data, BlockManager):
            mgr = self._init_mgr(data, axes=dict(index=index, columns=columns),
                                 dtype=dtype, copy=copy)
        elif isinstance(data, dict):
            mgr = self._init_dict(data, index, columns, dtype=dtype)
        elif isinstance(data, ma.MaskedArray):
            import numpy.ma.mrecords as mrecords
            # masked recarray
            if isinstance(data, mrecords.MaskedRecords):
                mgr = _masked_rec_array_to_mgr(data, index, columns, dtype,
                                               copy)

            # a masked array
            else:
                mask = ma.getmaskarray(data)
                if mask.any():
                    data, fill_value = _maybe_upcast(data, copy=True)
                    data[mask] = fill_value
                else:
                    data = data.copy()
                mgr = self._init_ndarray(data, index, columns, dtype=dtype,
                                         copy=copy)

        elif isinstance(data, (np.ndarray, Series, Index)):
            if data.dtype.names:
                data_columns = list(data.dtype.names)
                data = dict((k, data[k]) for k in data_columns)
                if columns is None:
                    columns = data_columns
                mgr = self._init_dict(data, index, columns, dtype=dtype)
            elif getattr(data, 'name', None):
                mgr = self._init_dict({data.name: data}, index, columns,
                                      dtype=dtype)
            else:
                mgr = self._init_ndarray(data, index, columns, dtype=dtype,
                                         copy=copy)
        elif isinstance(data, (list, types.GeneratorType)):
            if isinstance(data, types.GeneratorType):
                data = list(data)
            if len(data) > 0:
                if is_list_like(data[0]) and getattr(data[0], 'ndim', 1) == 1:
                    if is_named_tuple(data[0]) and columns is None:
                        columns = data[0]._fields
                    arrays, columns = _to_arrays(data, columns, dtype=dtype)
                    columns = _ensure_index(columns)

                    # set the index
                    if index is None:
                        if isinstance(data[0], Series):
                            index = _get_names_from_index(data)
                        elif isinstance(data[0], Categorical):
                            index = _default_index(len(data[0]))
                        else:
                            index = _default_index(len(data))

                    mgr = _arrays_to_mgr(arrays, columns, index, columns,
                                         dtype=dtype)
                else:
                    mgr = self._init_ndarray(data, index, columns, dtype=dtype,
                                             copy=copy)
            else:
                mgr = self._init_dict({}, index, columns, dtype=dtype)
        elif isinstance(data, collections.Iterator):
            raise TypeError("data argument can't be an iterator")
        else:
            try:
                arr = np.array(data, dtype=dtype, copy=copy)
            except (ValueError, TypeError) as e:
                exc = TypeError('DataFrame constructor called with '
                                'incompatible data and dtype: %s' % e)
                raise_with_traceback(exc)

            if arr.ndim == 0 and index is not None and columns is not None:
                if isinstance(data, compat.string_types) and dtype is None:
                    dtype = np.object_
                if dtype is None:
                    dtype, data = _infer_dtype_from_scalar(data)

                values = np.empty((len(index), len(columns)), dtype=dtype)
                values.fill(data)
                mgr = self._init_ndarray(values, index, columns, dtype=dtype,
                                         copy=False)
            else:
                raise PandasError('DataFrame constructor not properly called!')

        NDFrame.__init__(self, mgr, fastpath=True)

    def _init_dict(self, data, index, columns, dtype=None):
        """
        Segregate Series based on type and coerce into matrices.
        Needs to handle a lot of exceptional cases.
        """
        if columns is not None:
            columns = _ensure_index(columns)

            # GH10856
            # raise ValueError if only scalars in dict
            if index is None:
                extract_index(list(data.values()))

            # prefilter if columns passed
            data = dict((k, v) for k, v in compat.iteritems(data)
                        if k in columns)

            if index is None:
                index = extract_index(list(data.values()))

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
                        # 1783
                        v = np.empty(len(index), dtype=object)
                    elif np.issubdtype(dtype, np.flexible):
                        v = np.empty(len(index), dtype=object)
                    else:
                        v = np.empty(len(index), dtype=dtype)

                    v.fill(NA)
                else:
                    v = data[k]
                data_names.append(k)
                arrays.append(v)

        else:
            keys = list(data.keys())
            if not isinstance(data, OrderedDict):
                keys = _try_sort(keys)
            columns = data_names = Index(keys)
            arrays = [data[k] for k in keys]

        return _arrays_to_mgr(arrays, data_names, index, columns, dtype=dtype)

    def _init_ndarray(self, values, index, columns, dtype=None, copy=False):
        # input must be a ndarray, list, Series, index

        if isinstance(values, Series):
            if columns is None:
                if values.name is not None:
                    columns = [values.name]
            if index is None:
                index = values.index
            else:
                values = values.reindex(index)

            # zero len case (GH #2234)
            if not len(values) and columns is not None and len(columns):
                values = np.empty((0, 1), dtype=object)

        # helper to create the axes as indexes
        def _get_axes(N, K, index=index, columns=columns):
            # return axes or defaults

            if index is None:
                index = _default_index(N)
            else:
                index = _ensure_index(index)

            if columns is None:
                columns = _default_index(K)
            else:
                columns = _ensure_index(columns)
            return index, columns

        # we could have a categorical type passed or coerced to 'category'
        # recast this to an _arrays_to_mgr
        if (is_categorical_dtype(getattr(values, 'dtype', None)) or
                is_categorical_dtype(dtype)):

            if not hasattr(values, 'dtype'):
                values = _prep_ndarray(values, copy=copy)
                values = values.ravel()
            elif copy:
                values = values.copy()

            index, columns = _get_axes(len(values), 1)
            return _arrays_to_mgr([values], columns, index, columns,
                                  dtype=dtype)
        elif is_datetimetz(values):
            return self._init_dict({0: values}, index, columns, dtype=dtype)

        # by definition an array here
        # the dtypes will be coerced to a single dtype
        values = _prep_ndarray(values, copy=copy)

        if dtype is not None:
            if values.dtype != dtype:
                try:
                    values = values.astype(dtype)
                except Exception as orig:
                    e = ValueError("failed to cast to '%s' (Exception was: %s)"
                                   % (dtype, orig))
                    raise_with_traceback(e)

        index, columns = _get_axes(*values.shape)
        values = values.T

        # if we don't have a dtype specified, then try to convert objects
        # on the entire block; this is to convert if we have datetimelike's
        # embedded in an object type
        if dtype is None and is_object_dtype(values):
            values = _possibly_infer_to_datetimelike(values)

        return create_block_manager_from_blocks([values], [columns, index])

    @property
    def axes(self):
        """
        Return a list with the row axis labels and column axis labels as the
        only members. They are returned in that order.
        """
        return [self.index, self.columns]

    @property
    def shape(self):
        """
        Return a tuple representing the dimensionality of the DataFrame.
        """
        return len(self.index), len(self.columns)

    def _repr_fits_vertical_(self):
        """
        Check length against max_rows.
        """
        max_rows = get_option("display.max_rows")
        return len(self) <= max_rows

    def _repr_fits_horizontal_(self, ignore_width=False):
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

        # used by repr_html under IPython notebook or scripts ignore terminal
        # dims
        if ignore_width or not com.in_interactive_session():
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

        if not (max_rows is None):  # unlimited rows
            # min of two, where one may be None
            d = d.iloc[:min(max_rows, len(d))]
        else:
            return True

        d.to_string(buf=buf)
        value = buf.getvalue()
        repr_width = max([len(l) for l in value.split('\n')])

        return repr_width < width

    def _info_repr(self):
        """True if the repr should show the info view."""
        info_repr_option = (get_option("display.large_repr") == "info")
        return info_repr_option and not (self._repr_fits_horizontal_() and
                                         self._repr_fits_vertical_())

    def __unicode__(self):
        """
        Return a string representation for a particular DataFrame

        Invoked by unicode(df) in py2 only. Yields a Unicode String in both
        py2/py3.
        """
        buf = StringIO(u(""))
        if self._info_repr():
            self.info(buf=buf)
            return buf.getvalue()

        max_rows = get_option("display.max_rows")
        max_cols = get_option("display.max_columns")
        show_dimensions = get_option("display.show_dimensions")
        if get_option("display.expand_frame_repr"):
            width, _ = fmt.get_console_size()
        else:
            width = None
        self.to_string(buf=buf, max_rows=max_rows, max_cols=max_cols,
                       line_width=width, show_dimensions=show_dimensions)

        return buf.getvalue()

    def _repr_html_(self):
        """
        Return a html representation for a particular DataFrame.
        Mainly for IPython notebook.
        """
        # qtconsole doesn't report its line width, and also
        # behaves badly when outputting an HTML table
        # that doesn't fit the window, so disable it.
        # XXX: In IPython 3.x and above, the Qt console will not attempt to
        # display HTML, so this check can be removed when support for
        # IPython 2.x is no longer needed.
        if com.in_qtconsole():
            # 'HTML output is disabled in QtConsole'
            return None

        if self._info_repr():
            buf = StringIO(u(""))
            self.info(buf=buf)
            # need to escape the <class>, should be the first line.
            val = buf.getvalue().replace('<', r'&lt;', 1)
            val = val.replace('>', r'&gt;', 1)
            return '<pre>' + val + '</pre>'

        if get_option("display.notebook_repr_html"):
            max_rows = get_option("display.max_rows")
            max_cols = get_option("display.max_columns")
            show_dimensions = get_option("display.show_dimensions")

            return self.to_html(max_rows=max_rows, max_cols=max_cols,
                                show_dimensions=show_dimensions, notebook=True)
        else:
            return None

    def _repr_latex_(self):
        """
        Returns a LaTeX representation for a particular Dataframe.
        Mainly for use with nbconvert (jupyter notebook conversion to pdf).
        """
        if get_option('display.latex.repr'):
            return self.to_latex()
        else:
            return None

    @property
    def style(self):
        """
        Property returning a Styler object containing methods for
        building a styled HTML representation fo the DataFrame.

        See Also
        --------
        pandas.formats.style.Styler
        """
        from pandas.formats.style import Styler
        return Styler(self)

    def iteritems(self):
        """
        Iterator over (column name, Series) pairs.

        See also
        --------
        iterrows : Iterate over DataFrame rows as (index, Series) pairs.
        itertuples : Iterate over DataFrame rows as namedtuples of the values.

        """
        if self.columns.is_unique and hasattr(self, '_item_cache'):
            for k in self.columns:
                yield k, self._get_item_cache(k)
        else:
            for i, k in enumerate(self.columns):
                yield k, self._ixs(i, axis=1)

    def iterrows(self):
        """
        Iterate over DataFrame rows as (index, Series) pairs.

        Notes
        -----

        1. Because ``iterrows`` returns a Series for each row,
           it does **not** preserve dtypes across the rows (dtypes are
           preserved across columns for DataFrames). For example,

           >>> df = pd.DataFrame([[1, 1.5]], columns=['int', 'float'])
           >>> row = next(df.iterrows())[1]
           >>> row
           int      1.0
           float    1.5
           Name: 0, dtype: float64
           >>> print(row['int'].dtype)
           float64
           >>> print(df['int'].dtype)
           int64

           To preserve dtypes while iterating over the rows, it is better
           to use :meth:`itertuples` which returns namedtuples of the values
           and which is generally faster than ``iterrows``.

        2. You should **never modify** something you are iterating over.
           This is not guaranteed to work in all cases. Depending on the
           data types, the iterator returns a copy and not a view, and writing
           to it will have no effect.

        Returns
        -------
        it : generator
            A generator that iterates over the rows of the frame.

        See also
        --------
        itertuples : Iterate over DataFrame rows as namedtuples of the values.
        iteritems : Iterate over (column name, Series) pairs.

        """
        columns = self.columns
        klass = self._constructor_sliced
        for k, v in zip(self.index, self.values):
            s = klass(v, index=columns, name=k)
            yield k, s

    def itertuples(self, index=True, name="Pandas"):
        """
        Iterate over DataFrame rows as namedtuples, with index value as first
        element of the tuple.

        Parameters
        ----------
        index : boolean, default True
            If True, return the index as the first element of the tuple.
        name : string, default "Pandas"
            The name of the returned namedtuples or None to return regular
            tuples.

        Notes
        -----
        The column names will be renamed to positional names if they are
        invalid Python identifiers, repeated, or start with an underscore.
        With a large number of columns (>255), regular tuples are returned.

        See also
        --------
        iterrows : Iterate over DataFrame rows as (index, Series) pairs.
        iteritems : Iterate over (column name, Series) pairs.

        Examples
        --------

        >>> df = pd.DataFrame({'col1': [1, 2], 'col2': [0.1, 0.2]},
                              index=['a', 'b'])
        >>> df
           col1  col2
        a     1   0.1
        b     2   0.2
        >>> for row in df.itertuples():
        ...     print(row)
        ...
        Pandas(Index='a', col1=1, col2=0.10000000000000001)
        Pandas(Index='b', col1=2, col2=0.20000000000000001)

        """
        arrays = []
        fields = []
        if index:
            arrays.append(self.index)
            fields.append("Index")

        # use integer indexing because of possible duplicate column names
        arrays.extend(self.iloc[:, k] for k in range(len(self.columns)))

        # Python 3 supports at most 255 arguments to constructor, and
        # things get slow with this many fields in Python 2
        if name is not None and len(self.columns) + index < 256:
            # `rename` is unsupported in Python 2.6
            try:
                itertuple = collections.namedtuple(name,
                                                   fields + list(self.columns),
                                                   rename=True)
                return map(itertuple._make, zip(*arrays))
            except Exception:
                pass

        # fallback to regular tuples
        return zip(*arrays)

    if compat.PY3:  # pragma: no cover
        items = iteritems

    def __len__(self):
        """Returns length of info axis, but here we use the index """
        return len(self.index)

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
                raise ValueError('Dot product shape mismatch, %s vs %s' %
                                 (lvals.shape, rvals.shape))

        if isinstance(other, DataFrame):
            return self._constructor(np.dot(lvals, rvals), index=left.index,
                                     columns=other.columns)
        elif isinstance(other, Series):
            return Series(np.dot(lvals, rvals), index=left.index)
        elif isinstance(rvals, (np.ndarray, Index)):
            result = np.dot(lvals, rvals)
            if result.ndim == 2:
                return self._constructor(result, index=left.index)
            else:
                return Series(result, index=left.index)
        else:  # pragma: no cover
            raise TypeError('unsupported type: %s' % type(other))

    # ----------------------------------------------------------------------
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
        dtype : dtype, default None
            Data type to force, otherwise infer

        Returns
        -------
        DataFrame
        """
        index, columns = None, None
        orient = orient.lower()
        if orient == 'index':
            if len(data) > 0:
                # TODO speed up Series case
                if isinstance(list(data.values())[0], (Series, dict)):
                    data = _from_nested_dict(data)
                else:
                    data, index = list(data.values()), list(data.keys())
        elif orient != 'columns':  # pragma: no cover
            raise ValueError('only recognize index or columns for orient')

        return cls(data, index=index, columns=columns, dtype=dtype)

    def to_dict(self, orient='dict'):
        """Convert DataFrame to dictionary.

        Parameters
        ----------
        orient : str {'dict', 'list', 'series', 'split', 'records', 'index'}
            Determines the type of the values of the dictionary.

            - dict (default) : dict like {column -> {index -> value}}
            - list : dict like {column -> [values]}
            - series : dict like {column -> Series(values)}
            - split : dict like
              {index -> [index], columns -> [columns], data -> [values]}
            - records : list like
              [{column -> value}, ... , {column -> value}]
            - index : dict like {index -> {column -> value}}

              .. versionadded:: 0.17.0

            Abbreviations are allowed. `s` indicates `series` and `sp`
            indicates `split`.

        Returns
        -------
        result : dict like {column -> {index -> value}}
        """
        if not self.columns.is_unique:
            warnings.warn("DataFrame columns are not unique, some "
                          "columns will be omitted.", UserWarning)
        if orient.lower().startswith('d'):
            return dict((k, v.to_dict()) for k, v in compat.iteritems(self))
        elif orient.lower().startswith('l'):
            return dict((k, v.tolist()) for k, v in compat.iteritems(self))
        elif orient.lower().startswith('sp'):
            return {'index': self.index.tolist(),
                    'columns': self.columns.tolist(),
                    'data': lib.map_infer(self.values.ravel(),
                                          _maybe_box_datetimelike)
                    .reshape(self.values.shape).tolist()}
        elif orient.lower().startswith('s'):
            return dict((k, _maybe_box_datetimelike(v))
                        for k, v in compat.iteritems(self))
        elif orient.lower().startswith('r'):
            return [dict((k, _maybe_box_datetimelike(v))
                         for k, v in zip(self.columns, row))
                    for row in self.values]
        elif orient.lower().startswith('i'):
            return dict((k, v.to_dict()) for k, v in self.iterrows())
        else:
            raise ValueError("orient '%s' not understood" % orient)

    def to_gbq(self, destination_table, project_id, chunksize=10000,
               verbose=True, reauth=False, if_exists='fail', private_key=None):
        """Write a DataFrame to a Google BigQuery table.

        The main method a user calls to export pandas DataFrame contents to
        Google BigQuery table.

        Google BigQuery API Client Library v2 for Python is used.
        Documentation is available `here
        <https://developers.google.com/api-client-library/python/apis/bigquery/v2>`__

        Authentication to the Google BigQuery service is via OAuth 2.0.

        - If "private_key" is not provided:

          By default "application default credentials" are used.

          If default application credentials are not found or are restrictive,
          user account credentials are used. In this case, you will be asked to
          grant permissions for product name 'pandas GBQ'.

        - If "private_key" is provided:

          Service account credentials will be used to authenticate.

        Parameters
        ----------
        dataframe : DataFrame
            DataFrame to be written
        destination_table : string
            Name of table to be written, in the form 'dataset.tablename'
        project_id : str
            Google BigQuery Account project ID.
        chunksize : int (default 10000)
            Number of rows to be inserted in each chunk from the dataframe.
        verbose : boolean (default True)
            Show percentage complete
        reauth : boolean (default False)
            Force Google BigQuery to reauthenticate the user. This is useful
            if multiple accounts are used.
        if_exists : {'fail', 'replace', 'append'}, default 'fail'
            'fail': If table exists, do nothing.
            'replace': If table exists, drop it, recreate it, and insert data.
            'append': If table exists, insert data. Create if does not exist.
        private_key : str (optional)
            Service account private key in JSON format. Can be file path
            or string contents. This is useful for remote server
            authentication (eg. jupyter iPython notebook on remote host)
        """

        from pandas.io import gbq
        return gbq.to_gbq(self, destination_table, project_id=project_id,
                          chunksize=chunksize, verbose=verbose, reauth=reauth,
                          if_exists=if_exists, private_key=private_key)

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
        exclude : sequence, default None
            Columns or fields to exclude
        columns : sequence, default None
            Column names to use. If the passed data do not have names
            associated with them, this argument provides names for the
            columns. Otherwise this argument indicates the order of the columns
            in the result (any names not found in the data will become all-NA
            columns)
        coerce_float : boolean, default False
            Attempt to convert values of non-string, non-numeric objects (like
            decimal.Decimal) to floating point, useful for SQL result sets

        Returns
        -------
        df : DataFrame
        """

        # Make a copy of the input columns so we can modify it
        if columns is not None:
            columns = _ensure_index(columns)

        if is_iterator(data):
            if nrows == 0:
                return cls()

            try:
                first_row = next(data)
            except StopIteration:
                return cls(index=index, columns=columns)

            dtype = None
            if hasattr(first_row, 'dtype') and first_row.dtype.names:
                dtype = first_row.dtype

            values = [first_row]

            if nrows is None:
                values += data
            else:
                values.extend(itertools.islice(data, nrows - 1))

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
                for k, v in compat.iteritems(data):
                    if k in columns:
                        arr_columns.append(k)
                        arrays.append(v)

                arrays, arr_columns = _reorder_arrays(arrays, arr_columns,
                                                      columns)

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
            if (isinstance(index, compat.string_types) or
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

        mgr = _arrays_to_mgr(arrays, arr_columns, result_index, columns)

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
            if is_datetime64_any_dtype(self.index) and convert_datetime64:
                ix_vals = [self.index.to_pydatetime()]
            else:
                if isinstance(self.index, MultiIndex):
                    # array of tuples to numpy cols. copy copy copy
                    ix_vals = lmap(np.array, zip(*self.index.values))
                else:
                    ix_vals = [self.index.values]

            arrays = ix_vals + [self[c].get_values() for c in self.columns]

            count = 0
            index_names = list(self.index.names)
            if isinstance(self.index, MultiIndex):
                for i, n in enumerate(index_names):
                    if n is None:
                        index_names[i] = 'level_%d' % count
                        count += 1
            elif index_names[0] is None:
                index_names = ['index']
            names = (lmap(compat.text_type, index_names) +
                     lmap(compat.text_type, self.columns))
        else:
            arrays = [self[c].get_values() for c in self.columns]
            names = lmap(compat.text_type, self.columns)

        formats = [v.dtype for v in arrays]
        return np.rec.fromarrays(
            arrays,
            dtype={'names': names, 'formats': formats}
        )

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
        keys, values = lzip(*items)

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
                raise TypeError("Must pass columns with orient='index'")

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
    def from_csv(cls, path, header=0, sep=',', index_col=0, parse_dates=True,
                 encoding=None, tupleize_cols=False,
                 infer_datetime_format=False):
        """
        Read CSV file (DISCOURAGED, please use :func:`pandas.read_csv`
        instead).

        It is preferable to use the more powerful :func:`pandas.read_csv`
        for most general purposes, but ``from_csv`` makes for an easy
        roundtrip to and from a file (the exact counterpart of
        ``to_csv``), especially with a DataFrame of time series data.

        This method only differs from the preferred :func:`pandas.read_csv`
        in some defaults:

        - `index_col` is ``0`` instead of ``None`` (take first column as index
          by default)
        - `parse_dates` is ``True`` instead of ``False`` (try parsing the index
          as datetime by default)

        So a ``pd.DataFrame.from_csv(path)`` can be replaced by
        ``pd.read_csv(path, index_col=0, parse_dates=True)``.

        Parameters
        ----------
        path : string file path or file handle / StringIO
        header : int, default 0
            Row to use as header (skip prior rows)
        sep : string, default ','
            Field delimiter
        index_col : int or sequence, default 0
            Column to use for index. If a sequence is given, a MultiIndex
            is used. Different default from read_table
        parse_dates : boolean, default True
            Parse dates. Different default from read_table
        tupleize_cols : boolean, default False
            write multi_index columns as a list of tuples (if True)
            or new (expanded format) if False)
        infer_datetime_format: boolean, default False
            If True and `parse_dates` is True for a column, try to infer the
            datetime format based on the first datetime string. If the format
            can be inferred, there often will be a large parsing speed-up.

        See also
        --------
        pandas.read_csv

        Returns
        -------
        y : DataFrame

        """
        from pandas.io.parsers import read_table
        return read_table(path, header=header, sep=sep,
                          parse_dates=parse_dates, index_col=index_col,
                          encoding=encoding, tupleize_cols=tupleize_cols,
                          infer_datetime_format=infer_datetime_format)

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
                               columns=self.columns, default_kind=kind,
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
        # only support this kind for now
        if (not isinstance(self.index, MultiIndex) or  # pragma: no cover
                len(self.index.levels) != 2):
            raise NotImplementedError('Only 2-level MultiIndex are supported.')

        if not self.index.is_unique:
            raise ValueError("Can't convert non-uniquely indexed "
                             "DataFrame to Panel")

        self._consolidate_inplace()

        # minor axis must be sorted
        if self.index.lexsort_depth < 2:
            selfsorted = self.sort_index(level=0)
        else:
            selfsorted = self

        major_axis, minor_axis = selfsorted.index.levels
        major_labels, minor_labels = selfsorted.index.labels
        shape = len(major_axis), len(minor_axis)

        # preserve names, if any
        major_axis = major_axis.copy()
        major_axis.name = self.index.names[0]

        minor_axis = minor_axis.copy()
        minor_axis.name = self.index.names[1]

        # create new axes
        new_axes = [selfsorted.columns, major_axis, minor_axis]

        # create new manager
        new_mgr = selfsorted._data.reshape_nd(axes=new_axes,
                                              labels=[major_labels,
                                                      minor_labels],
                                              shape=shape,
                                              ref_items=selfsorted.columns)

        return self._constructor_expanddim(new_mgr)

    def to_csv(self, path_or_buf=None, sep=",", na_rep='', float_format=None,
               columns=None, header=True, index=True, index_label=None,
               mode='w', encoding=None, compression=None, quoting=None,
               quotechar='"', line_terminator='\n', chunksize=None,
               tupleize_cols=False, date_format=None, doublequote=True,
               escapechar=None, decimal='.'):
        r"""Write DataFrame to a comma-separated values (csv) file

        Parameters
        ----------
        path_or_buf : string or file handle, default None
            File path or object, if None is provided the result is returned as
            a string.
        sep : character, default ','
            Field delimiter for the output file.
        na_rep : string, default ''
            Missing data representation
        float_format : string, default None
            Format string for floating point numbers
        columns : sequence, optional
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
        mode : str
            Python write mode, default 'w'
        encoding : string, optional
            A string representing the encoding to use in the output file,
            defaults to 'ascii' on Python 2 and 'utf-8' on Python 3.
        compression : string, optional
            a string representing the compression to use in the output file,
            allowed values are 'gzip', 'bz2', 'xz',
            only used when the first argument is a filename
        line_terminator : string, default ``'\n'``
            The newline character or character sequence to use in the output
            file
        quoting : optional constant from csv module
            defaults to csv.QUOTE_MINIMAL. If you have set a `float_format`
            then floats are converted to strings and thus csv.QUOTE_NONNUMERIC
            will treat them as non-numeric
        quotechar : string (length 1), default '\"'
            character used to quote fields
        doublequote : boolean, default True
            Control quoting of `quotechar` inside a field
        escapechar : string (length 1), default None
            character used to escape `sep` and `quotechar` when appropriate
        chunksize : int or None
            rows to write at a time
        tupleize_cols : boolean, default False
            write multi_index columns as a list of tuples (if True)
            or new (expanded format) if False)
        date_format : string, default None
            Format string for datetime objects
        decimal: string, default '.'
            Character recognized as decimal separator. E.g. use ',' for
            European data

            .. versionadded:: 0.16.0

        """
        formatter = fmt.CSVFormatter(self, path_or_buf,
                                     line_terminator=line_terminator, sep=sep,
                                     encoding=encoding,
                                     compression=compression, quoting=quoting,
                                     na_rep=na_rep, float_format=float_format,
                                     cols=columns, header=header, index=index,
                                     index_label=index_label, mode=mode,
                                     chunksize=chunksize, quotechar=quotechar,
                                     tupleize_cols=tupleize_cols,
                                     date_format=date_format,
                                     doublequote=doublequote,
                                     escapechar=escapechar, decimal=decimal)
        formatter.save()

        if path_or_buf is None:
            return formatter.path_or_buf.getvalue()

    @Appender(_shared_docs['to_excel'] % _shared_doc_kwargs)
    def to_excel(self, excel_writer, sheet_name='Sheet1', na_rep='',
                 float_format=None, columns=None, header=True, index=True,
                 index_label=None, startrow=0, startcol=0, engine=None,
                 merge_cells=True, encoding=None, inf_rep='inf', verbose=True,
                 freeze_panes=None):
        from pandas.io.excel import ExcelWriter
        need_save = False
        if encoding is None:
            encoding = 'ascii'

        if isinstance(excel_writer, compat.string_types):
            excel_writer = ExcelWriter(excel_writer, engine=engine)
            need_save = True

        formatter = fmt.ExcelFormatter(self, na_rep=na_rep, cols=columns,
                                       header=header,
                                       float_format=float_format, index=index,
                                       index_label=index_label,
                                       merge_cells=merge_cells,
                                       inf_rep=inf_rep)

        formatted_cells = formatter.get_formatted_cells()
        excel_writer.write_cells(formatted_cells, sheet_name,
                                 startrow=startrow, startcol=startcol,
                                 freeze_panes=freeze_panes)
        if need_save:
            excel_writer.save()

    def to_stata(self, fname, convert_dates=None, write_index=True,
                 encoding="latin-1", byteorder=None, time_stamp=None,
                 data_label=None, variable_labels=None):
        """
        A class for writing Stata binary dta files from array-like objects

        Parameters
        ----------
        fname : str or buffer
            String path of file-like object
        convert_dates : dict
            Dictionary mapping columns containing datetime types to stata
            internal format to use when wirting the dates. Options are 'tc',
            'td', 'tm', 'tw', 'th', 'tq', 'ty'. Column can be either an integer
            or a name. Datetime columns that do not have a conversion type
            specified will be converted to 'tc'. Raises NotImplementedError if
            a datetime column has timezone information
        write_index : bool
            Write the index to Stata dataset.
        encoding : str
            Default is latin-1. Unicode is not supported
        byteorder : str
            Can be ">", "<", "little", or "big". default is `sys.byteorder`
        time_stamp : datetime
            A datetime to use as file creation date.  Default is the current
            time.
        dataset_label : str
            A label for the data set.  Must be 80 characters or smaller.
        variable_labels : dict
            Dictionary containing columns as keys and variable labels as
            values. Each label must be 80 characters or smaller.

            .. versionadded:: 0.19.0

        Raises
        ------
        NotImplementedError
            * If datetimes contain timezone information
            * Column dtype is not representable in Stata
        ValueError
            * Columns listed in convert_dates are noth either datetime64[ns]
              or datetime.datetime
            * Column listed in convert_dates is not in DataFrame
            * Categorical label contains more than 32,000 characters

            .. versionadded:: 0.19.0

        Examples
        --------
        >>> writer = StataWriter('./data_file.dta', data)
        >>> writer.write_file()

        Or with dates

        >>> writer = StataWriter('./date_data_file.dta', data, {2 : 'tw'})
        >>> writer.write_file()
        """
        from pandas.io.stata import StataWriter
        writer = StataWriter(fname, self, convert_dates=convert_dates,
                             encoding=encoding, byteorder=byteorder,
                             time_stamp=time_stamp, data_label=data_label,
                             write_index=write_index,
                             variable_labels=variable_labels)
        writer.write_file()

    def to_feather(self, fname):
        """
        write out the binary feather-format for DataFrames

        .. versionadded:: 0.20.0

        Parameters
        ----------
        fname : str
            string file path

        """
        from pandas.io.feather_format import to_feather
        to_feather(self, fname)

    @Appender(fmt.docstring_to_string, indents=1)
    def to_string(self, buf=None, columns=None, col_space=None, header=True,
                  index=True, na_rep='NaN', formatters=None, float_format=None,
                  sparsify=None, index_names=True, justify=None,
                  line_width=None, max_rows=None, max_cols=None,
                  show_dimensions=False):
        """
        Render a DataFrame to a console-friendly tabular output.
        """

        formatter = fmt.DataFrameFormatter(self, buf=buf, columns=columns,
                                           col_space=col_space, na_rep=na_rep,
                                           formatters=formatters,
                                           float_format=float_format,
                                           sparsify=sparsify, justify=justify,
                                           index_names=index_names,
                                           header=header, index=index,
                                           line_width=line_width,
                                           max_rows=max_rows,
                                           max_cols=max_cols,
                                           show_dimensions=show_dimensions)
        formatter.to_string()

        if buf is None:
            result = formatter.buf.getvalue()
            return result

    @Appender(fmt.docstring_to_string, indents=1)
    def to_html(self, buf=None, columns=None, col_space=None, header=True,
                index=True, na_rep='NaN', formatters=None, float_format=None,
                sparsify=None, index_names=True, justify=None, bold_rows=True,
                classes=None, escape=True, max_rows=None, max_cols=None,
                show_dimensions=False, notebook=False, decimal='.',
                border=None):
        """
        Render a DataFrame as an HTML table.

        `to_html`-specific options:

        bold_rows : boolean, default True
            Make the row labels bold in the output
        classes : str or list or tuple, default None
            CSS class(es) to apply to the resulting html table
        escape : boolean, default True
            Convert the characters <, >, and & to HTML-safe sequences.=
        max_rows : int, optional
            Maximum number of rows to show before truncating. If None, show
            all.
        max_cols : int, optional
            Maximum number of columns to show before truncating. If None, show
            all.
        decimal : string, default '.'
            Character recognized as decimal separator, e.g. ',' in Europe

            .. versionadded:: 0.18.0
        border : int
            A ``border=border`` attribute is included in the opening
            `<table>` tag. Default ``pd.options.html.border``.

            .. versionadded:: 0.19.0
        """

        formatter = fmt.DataFrameFormatter(self, buf=buf, columns=columns,
                                           col_space=col_space, na_rep=na_rep,
                                           formatters=formatters,
                                           float_format=float_format,
                                           sparsify=sparsify, justify=justify,
                                           index_names=index_names,
                                           header=header, index=index,
                                           bold_rows=bold_rows, escape=escape,
                                           max_rows=max_rows,
                                           max_cols=max_cols,
                                           show_dimensions=show_dimensions,
                                           decimal=decimal)
        # TODO: a generic formatter wld b in DataFrameFormatter
        formatter.to_html(classes=classes, notebook=notebook, border=border)

        if buf is None:
            return formatter.buf.getvalue()

    @Appender(fmt.common_docstring + fmt.return_docstring, indents=1)
    def to_latex(self, buf=None, columns=None, col_space=None, header=True,
                 index=True, na_rep='NaN', formatters=None, float_format=None,
                 sparsify=None, index_names=True, bold_rows=True,
                 column_format=None, longtable=None, escape=None,
                 encoding=None, decimal='.', multicolumn=None,
                 multicolumn_format=None, multirow=None):
        r"""
        Render a DataFrame to a tabular environment table. You can splice
        this into a LaTeX document. Requires \usepackage{booktabs}.

        `to_latex`-specific options:

        bold_rows : boolean, default True
            Make the row labels bold in the output
        column_format : str, default None
            The columns format as specified in `LaTeX table format
            <https://en.wikibooks.org/wiki/LaTeX/Tables>`__ e.g 'rcl' for 3
            columns
        longtable : boolean, default will be read from the pandas config module
            Default: False.
            Use a longtable environment instead of tabular. Requires adding
            a \usepackage{longtable} to your LaTeX preamble.
        escape : boolean, default will be read from the pandas config module
            Default: True.
            When set to False prevents from escaping latex special
            characters in column names.
        encoding : str, default None
            A string representing the encoding to use in the output file,
            defaults to 'ascii' on Python 2 and 'utf-8' on Python 3.
        decimal : string, default '.'
            Character recognized as decimal separator, e.g. ',' in Europe.

            .. versionadded:: 0.18.0

        multicolumn : boolean, default True
            Use \multicolumn to enhance MultiIndex columns.
            The default will be read from the config module.

            .. versionadded:: 0.20.0

        multicolumn_format : str, default 'l'
            The alignment for multicolumns, similar to `column_format`
            The default will be read from the config module.

            .. versionadded:: 0.20.0

        multirow : boolean, default False
            Use \multirow to enhance MultiIndex rows.
            Requires adding a \usepackage{multirow} to your LaTeX preamble.
            Will print centered labels (instead of top-aligned)
            across the contained rows, separating groups via clines.
            The default will be read from the pandas config module.

            .. versionadded:: 0.20.0

        """
        # Get defaults from the pandas config
        if longtable is None:
            longtable = get_option("display.latex.longtable")
        if escape is None:
            escape = get_option("display.latex.escape")
        if multicolumn is None:
            multicolumn = get_option("display.latex.multicolumn")
        if multicolumn_format is None:
            multicolumn_format = get_option("display.latex.multicolumn_format")
        if multirow is None:
            multirow = get_option("display.latex.multirow")

        formatter = fmt.DataFrameFormatter(self, buf=buf, columns=columns,
                                           col_space=col_space, na_rep=na_rep,
                                           header=header, index=index,
                                           formatters=formatters,
                                           float_format=float_format,
                                           bold_rows=bold_rows,
                                           sparsify=sparsify,
                                           index_names=index_names,
                                           escape=escape, decimal=decimal)
        formatter.to_latex(column_format=column_format, longtable=longtable,
                           encoding=encoding, multicolumn=multicolumn,
                           multicolumn_format=multicolumn_format,
                           multirow=multirow)

        if buf is None:
            return formatter.buf.getvalue()

    def info(self, verbose=None, buf=None, max_cols=None, memory_usage=None,
             null_counts=None):
        """
        Concise summary of a DataFrame.

        Parameters
        ----------
        verbose : {None, True, False}, optional
            Whether to print the full summary.
            None follows the `display.max_info_columns` setting.
            True or False overrides the `display.max_info_columns` setting.
        buf : writable buffer, defaults to sys.stdout
        max_cols : int, default None
            Determines whether full summary or short summary is printed.
            None follows the `display.max_info_columns` setting.
        memory_usage : boolean/string, default None
            Specifies whether total memory usage of the DataFrame
            elements (including index) should be displayed. None follows
            the `display.memory_usage` setting. True or False overrides
            the `display.memory_usage` setting. A value of 'deep' is equivalent
            of True, with deep introspection. Memory usage is shown in
            human-readable units (base-2 representation).
        null_counts : boolean, default None
            Whether to show the non-null counts

            - If None, then only show if the frame is smaller than
              max_info_rows and max_info_columns.
            - If True, always show counts.
            - If False, never show counts.

        """
        from pandas.formats.format import _put_lines

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
            max_cols = get_option('display.max_info_columns',
                                  len(self.columns) + 1)

        max_rows = get_option('display.max_info_rows', len(self) + 1)

        if null_counts is None:
            show_counts = ((len(self.columns) <= max_cols) and
                           (len(self) < max_rows))
        else:
            show_counts = null_counts
        exceeds_info_cols = len(self.columns) > max_cols

        def _verbose_repr():
            lines.append('Data columns (total %d columns):' %
                         len(self.columns))
            space = max([len(pprint_thing(k)) for k in self.columns]) + 4
            counts = None

            tmpl = "%s%s"
            if show_counts:
                counts = self.count()
                if len(cols) != len(counts):  # pragma: no cover
                    raise AssertionError('Columns must equal counts (%d != %d)'
                                         % (len(cols), len(counts)))
                tmpl = "%s non-null %s"

            dtypes = self.dtypes
            for i, col in enumerate(self.columns):
                dtype = dtypes.iloc[i]
                col = pprint_thing(col)

                count = ""
                if show_counts:
                    count = counts.iloc[i]

                lines.append(_put_str(col, space) + tmpl % (count, dtype))

        def _non_verbose_repr():
            lines.append(self.columns.summary(name='Columns'))

        def _sizeof_fmt(num, size_qualifier):
            # returns size in human readable format
            for x in ['bytes', 'KB', 'MB', 'GB', 'TB']:
                if num < 1024.0:
                    return "%3.1f%s %s" % (num, size_qualifier, x)
                num /= 1024.0
            return "%3.1f%s %s" % (num, size_qualifier, 'PB')

        if verbose:
            _verbose_repr()
        elif verbose is False:  # specifically set to False, not nesc None
            _non_verbose_repr()
        else:
            if exceeds_info_cols:
                _non_verbose_repr()
            else:
                _verbose_repr()

        counts = self.get_dtype_counts()
        dtypes = ['%s(%d)' % k for k in sorted(compat.iteritems(counts))]
        lines.append('dtypes: %s' % ', '.join(dtypes))

        if memory_usage is None:
            memory_usage = get_option('display.memory_usage')
        if memory_usage:
            # append memory usage of df to display
            size_qualifier = ''
            if memory_usage == 'deep':
                deep = True
            else:
                # size_qualifier is just a best effort; not guaranteed to catch
                # all cases (e.g., it misses categorical data even with object
                # categories)
                deep = False
                if ('object' in counts or
                        self.index._is_memory_usage_qualified()):
                    size_qualifier = '+'
            mem_usage = self.memory_usage(index=True, deep=deep).sum()
            lines.append("memory usage: %s\n" %
                         _sizeof_fmt(mem_usage, size_qualifier))
        _put_lines(buf, lines)

    def memory_usage(self, index=True, deep=False):
        """Memory usage of DataFrame columns.

        Parameters
        ----------
        index : bool
            Specifies whether to include memory usage of DataFrame's
            index in returned Series. If `index=True` (default is False)
            the first index of the Series is `Index`.
        deep : bool
            Introspect the data deeply, interrogate
            `object` dtypes for system-level memory consumption

        Returns
        -------
        sizes : Series
            A series with column names as index and memory usage of
            columns with units of bytes.

        Notes
        -----
        Memory usage does not include memory consumed by elements that
        are not components of the array if deep=False

        See Also
        --------
        numpy.ndarray.nbytes
        """
        result = Series([c.memory_usage(index=False, deep=deep)
                         for col, c in self.iteritems()], index=self.columns)
        if index:
            result = Series(self.index.memory_usage(deep=deep),
                            index=['Index']).append(result)
        return result

    def transpose(self, *args, **kwargs):
        """Transpose index and columns"""
        nv.validate_transpose(args, dict())
        return super(DataFrame, self).transpose(1, 0, **kwargs)

    T = property(transpose)

    # ----------------------------------------------------------------------
    # Picklability

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
                                columns=_unpickle_array(ocols), copy=False)

            dm = dm.join(objects)

        self._data = dm._data

    # ----------------------------------------------------------------------
    # Getting and setting elements

    def get_value(self, index, col, takeable=False):
        """
        Quickly retrieve single value at passed column and index

        Parameters
        ----------
        index : row label
        col : column label
        takeable : interpret the index/col as indexers, default False

        Returns
        -------
        value : scalar value
        """

        if takeable:
            series = self._iget_item_cache(col)
            return _maybe_box_datetimelike(series._values[index])

        series = self._get_item_cache(col)
        engine = self.index._engine
        return engine.get_value(series.get_values(), index)

    def set_value(self, index, col, value, takeable=False):
        """
        Put single value at passed column and index

        Parameters
        ----------
        index : row label
        col : column label
        value : scalar value
        takeable : interpret the index/col as indexers, default False

        Returns
        -------
        frame : DataFrame
            If label pair is contained, will be reference to calling DataFrame,
            otherwise a new object
        """
        try:
            if takeable is True:
                series = self._iget_item_cache(col)
                return series.set_value(index, value, takeable=True)

            series = self._get_item_cache(col)
            engine = self.index._engine
            engine.set_value(series._values, index, value)
            return self
        except (KeyError, TypeError):

            # set using a non-recursive method & reset the cache
            self.loc[index, col] = value
            self._item_cache.pop(col, None)

            return self

    def _ixs(self, i, axis=0):
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
                    result = self.take(i, axis=axis)
                    copy = True
                else:
                    new_values = self._data.fast_xs(i)
                    if is_scalar(new_values):
                        return new_values

                    # if we are a copy, mark as such
                    copy = (isinstance(new_values, np.ndarray) and
                            new_values.base is None)
                    result = self._constructor_sliced(new_values,
                                                      index=self.columns,
                                                      name=self.index[i],
                                                      dtype=new_values.dtype)
                result._set_is_copy(self, copy=copy)
                return result

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
                return self.loc[:, lab_slice]
            else:
                if isinstance(label, Index):
                    return self.take(i, axis=1, convert=True)

                index_len = len(self.index)

                # if the values returned are not the same length
                # as the index (iow a not found value), iget returns
                # a 0-len ndarray. This is effectively catching
                # a numpy error (as numpy should really raise)
                values = self._data.iget(i)

                if index_len and not len(values):
                    values = np.array([np.nan] * index_len, dtype=object)
                result = self._constructor_sliced.from_array(values,
                                                             index=self.index,
                                                             name=label,
                                                             fastpath=True)

                # this is a cached value, mark it so
                result._set_as_cached(label, self)

                return result

    def __getitem__(self, key):
        key = com._apply_if_callable(key, self)

        # shortcut if we are an actual column
        is_mi_columns = isinstance(self.columns, MultiIndex)
        try:
            if key in self.columns and not is_mi_columns:
                return self._getitem_column(key)
        except:
            pass

        # see if we can slice the rows
        indexer = convert_to_index_sliceable(self, key)
        if indexer is not None:
            return self._getitem_slice(indexer)

        if isinstance(key, (Series, np.ndarray, Index, list)):
            # either boolean or fancy integer index
            return self._getitem_array(key)
        elif isinstance(key, DataFrame):
            return self._getitem_frame(key)
        elif is_mi_columns:
            return self._getitem_multilevel(key)
        else:
            return self._getitem_column(key)

    def _getitem_column(self, key):
        """ return the actual column """

        # get column
        if self.columns.is_unique:
            return self._get_item_cache(key)

        # duplicate columns & possible reduce dimensionality
        result = self._constructor(self._data.get(key))
        if result.columns.is_unique:
            result = result[key]

        return result

    def _getitem_slice(self, key):
        return self._slice(key, axis=0)

    def _getitem_array(self, key):
        # also raises Exception if object array with NA values
        if com.is_bool_indexer(key):
            # warning here just in case -- previously __setitem__ was
            # reindexing but __getitem__ was not; it seems more reasonable to
            # go with the __setitem__ behavior since that is more consistent
            # with all other indexing behavior
            if isinstance(key, Series) and not key.index.equals(self.index):
                warnings.warn("Boolean Series key will be reindexed to match "
                              "DataFrame index.", UserWarning, stacklevel=3)
            elif len(key) != len(self.index):
                raise ValueError('Item wrong length %d instead of %d.' %
                                 (len(key), len(self.index)))
            # check_bool_indexer will throw exception if Series key cannot
            # be reindexed to match DataFrame rows
            key = check_bool_indexer(self.index, key)
            indexer = key.nonzero()[0]
            return self.take(indexer, axis=0, convert=False)
        else:
            indexer = self.loc._convert_to_indexer(key, axis=1)
            return self.take(indexer, axis=1, convert=True)

    def _getitem_multilevel(self, key):
        loc = self.columns.get_loc(key)
        if isinstance(loc, (slice, Series, np.ndarray, Index)):
            new_columns = self.columns[loc]
            result_columns = maybe_droplevels(new_columns, key)
            if self._is_mixed_type:
                result = self.reindex(columns=new_columns)
                result.columns = result_columns
            else:
                new_values = self.values[:, loc]
                result = self._constructor(new_values, index=self.index,
                                           columns=result_columns)
                result = result.__finalize__(self)
            if len(result.columns) == 1:
                top = result.columns[0]
                if ((type(top) == str and top == '') or
                        (type(top) == tuple and top[0] == '')):
                    result = result['']
                    if isinstance(result, Series):
                        result = self._constructor_sliced(result,
                                                          index=self.index,
                                                          name=key)

            result._set_is_copy(self)
            return result
        else:
            return self._get_item_cache(key)

    def _getitem_frame(self, key):
        if key.values.size and not is_bool_dtype(key.values):
            raise ValueError('Must pass DataFrame with boolean values only')
        return self.where(key)

    def query(self, expr, inplace=False, **kwargs):
        """Query the columns of a frame with a boolean expression.

        .. versionadded:: 0.13

        Parameters
        ----------
        expr : string
            The query string to evaluate.  You can refer to variables
            in the environment by prefixing them with an '@' character like
            ``@a + b``.
        inplace : bool
            Whether the query should modify the data in place or return
            a modified copy

            .. versionadded:: 0.18.0

        kwargs : dict
            See the documentation for :func:`pandas.eval` for complete details
            on the keyword arguments accepted by :meth:`DataFrame.query`.

        Returns
        -------
        q : DataFrame

        Notes
        -----
        The result of the evaluation of this expression is first passed to
        :attr:`DataFrame.loc` and if that fails because of a
        multidimensional key (e.g., a DataFrame) then the result will be passed
        to :meth:`DataFrame.__getitem__`.

        This method uses the top-level :func:`pandas.eval` function to
        evaluate the passed query.

        The :meth:`~pandas.DataFrame.query` method uses a slightly
        modified Python syntax by default. For example, the ``&`` and ``|``
        (bitwise) operators have the precedence of their boolean cousins,
        :keyword:`and` and :keyword:`or`. This *is* syntactically valid Python,
        however the semantics are different.

        You can change the semantics of the expression by passing the keyword
        argument ``parser='python'``. This enforces the same semantics as
        evaluation in Python space. Likewise, you can pass ``engine='python'``
        to evaluate an expression using Python itself as a backend. This is not
        recommended as it is inefficient compared to using ``numexpr`` as the
        engine.

        The :attr:`DataFrame.index` and
        :attr:`DataFrame.columns` attributes of the
        :class:`~pandas.DataFrame` instance are placed in the query namespace
        by default, which allows you to treat both the index and columns of the
        frame as a column in the frame.
        The identifier ``index`` is used for the frame index; you can also
        use the name of the index to identify it in a query.

        For further details and examples see the ``query`` documentation in
        :ref:`indexing <indexing.query>`.

        See Also
        --------
        pandas.eval
        DataFrame.eval

        Examples
        --------
        >>> from numpy.random import randn
        >>> from pandas import DataFrame
        >>> df = DataFrame(randn(10, 2), columns=list('ab'))
        >>> df.query('a > b')
        >>> df[df.a > df.b]  # same result as the previous expression
        """
        inplace = validate_bool_kwarg(inplace, 'inplace')
        if not isinstance(expr, compat.string_types):
            msg = "expr must be a string to be evaluated, {0} given"
            raise ValueError(msg.format(type(expr)))
        kwargs['level'] = kwargs.pop('level', 0) + 1
        kwargs['target'] = None
        res = self.eval(expr, **kwargs)

        try:
            new_data = self.loc[res]
        except ValueError:
            # when res is multi-dimensional loc raises, but this is sometimes a
            # valid query
            new_data = self[res]

        if inplace:
            self._update_inplace(new_data)
        else:
            return new_data

    def eval(self, expr, inplace=None, **kwargs):
        """Evaluate an expression in the context of the calling DataFrame
        instance.

        Parameters
        ----------
        expr : string
            The expression string to evaluate.
        inplace : bool
            If the expression contains an assignment, whether to return a new
            DataFrame or mutate the existing.

            WARNING: inplace=None currently falls back to to True, but
            in a future version, will default to False.  Use inplace=True
            explicitly rather than relying on the default.

            .. versionadded:: 0.18.0

        kwargs : dict
            See the documentation for :func:`~pandas.eval` for complete details
            on the keyword arguments accepted by
            :meth:`~pandas.DataFrame.query`.

        Returns
        -------
        ret : ndarray, scalar, or pandas object

        See Also
        --------
        pandas.DataFrame.query
        pandas.DataFrame.assign
        pandas.eval

        Notes
        -----
        For more details see the API documentation for :func:`~pandas.eval`.
        For detailed examples see :ref:`enhancing performance with eval
        <enhancingperf.eval>`.

        Examples
        --------
        >>> from numpy.random import randn
        >>> from pandas import DataFrame
        >>> df = DataFrame(randn(10, 2), columns=list('ab'))
        >>> df.eval('a + b')
        >>> df.eval('c = a + b')
        """
        inplace = validate_bool_kwarg(inplace, 'inplace')
        resolvers = kwargs.pop('resolvers', None)
        kwargs['level'] = kwargs.pop('level', 0) + 1
        if resolvers is None:
            index_resolvers = self._get_index_resolvers()
            resolvers = dict(self.iteritems()), index_resolvers
        if 'target' not in kwargs:
            kwargs['target'] = self
        kwargs['resolvers'] = kwargs.get('resolvers', ()) + tuple(resolvers)
        return _eval(expr, inplace=inplace, **kwargs)

    def select_dtypes(self, include=None, exclude=None):
        """Return a subset of a DataFrame including/excluding columns based on
        their ``dtype``.

        Parameters
        ----------
        include, exclude : list-like
            A list of dtypes or strings to be included/excluded. You must pass
            in a non-empty sequence for at least one of these.

        Raises
        ------
        ValueError
            * If both of ``include`` and ``exclude`` are empty
            * If ``include`` and ``exclude`` have overlapping elements
            * If any kind of string dtype is passed in.
        TypeError
            * If either of ``include`` or ``exclude`` is not a sequence

        Returns
        -------
        subset : DataFrame
            The subset of the frame including the dtypes in ``include`` and
            excluding the dtypes in ``exclude``.

        Notes
        -----
        * To select all *numeric* types use the numpy dtype ``numpy.number``
        * To select strings you must use the ``object`` dtype, but note that
          this will return *all* object dtype columns
        * See the `numpy dtype hierarchy
          <http://docs.scipy.org/doc/numpy/reference/arrays.scalars.html>`__
        * To select datetimes, use np.datetime64, 'datetime' or 'datetime64'
        * To select timedeltas, use np.timedelta64, 'timedelta' or
          'timedelta64'
        * To select Pandas categorical dtypes, use 'category'
        * To select Pandas datetimetz dtypes, use 'datetimetz' (new in 0.20.0),
          or a 'datetime64[ns, tz]' string

        Examples
        --------
        >>> df = pd.DataFrame({'a': np.random.randn(6).astype('f4'),
        ...                    'b': [True, False] * 3,
        ...                    'c': [1.0, 2.0] * 3})
        >>> df
                a      b  c
        0  0.3962   True  1
        1  0.1459  False  2
        2  0.2623   True  1
        3  0.0764  False  2
        4 -0.9703   True  1
        5 -1.2094  False  2
        >>> df.select_dtypes(include=['float64'])
           c
        0  1
        1  2
        2  1
        3  2
        4  1
        5  2
        >>> df.select_dtypes(exclude=['floating'])
               b
        0   True
        1  False
        2   True
        3  False
        4   True
        5  False
        """
        include, exclude = include or (), exclude or ()
        if not (is_list_like(include) and is_list_like(exclude)):
            raise TypeError('include and exclude must both be non-string'
                            ' sequences')
        selection = tuple(map(frozenset, (include, exclude)))

        if not any(selection):
            raise ValueError('at least one of include or exclude must be '
                             'nonempty')

        # convert the myriad valid dtypes object to a single representation
        include, exclude = map(
            lambda x: frozenset(map(_get_dtype_from_object, x)), selection)
        for dtypes in (include, exclude):
            _invalidate_string_dtypes(dtypes)

        # can't both include AND exclude!
        if not include.isdisjoint(exclude):
            raise ValueError('include and exclude overlap on %s' %
                             (include & exclude))

        # empty include/exclude -> defaults to True
        # three cases (we've already raised if both are empty)
        # case 1: empty include, nonempty exclude
        # we have True, True, ... True for include, same for exclude
        # in the loop below we get the excluded
        # and when we call '&' below we get only the excluded
        # case 2: nonempty include, empty exclude
        # same as case 1, but with include
        # case 3: both nonempty
        # the "union" of the logic of case 1 and case 2:
        # we get the included and excluded, and return their logical and
        include_these = Series(not bool(include), index=self.columns)
        exclude_these = Series(not bool(exclude), index=self.columns)

        def is_dtype_instance_mapper(column, dtype):
            return column, functools.partial(issubclass, dtype.type)

        for column, f in itertools.starmap(is_dtype_instance_mapper,
                                           self.dtypes.iteritems()):
            if include:  # checks for the case of empty include or exclude
                include_these[column] = any(map(f, include))
            if exclude:
                exclude_these[column] = not any(map(f, exclude))

        dtype_indexer = include_these & exclude_these
        return self.loc[com._get_info_slice(self, dtype_indexer)]

    def _box_item_values(self, key, values):
        items = self.columns[self.columns.get_loc(key)]
        if values.ndim == 2:
            return self._constructor(values.T, columns=items, index=self.index)
        else:
            return self._box_col_values(values, items)

    def _box_col_values(self, values, items):
        """ provide boxed values for a column """
        return self._constructor_sliced.from_array(values, index=self.index,
                                                   name=items, fastpath=True)

    def __setitem__(self, key, value):
        key = com._apply_if_callable(key, self)

        # see if we can slice the rows
        indexer = convert_to_index_sliceable(self, key)
        if indexer is not None:
            return self._setitem_slice(indexer, value)

        if isinstance(key, (Series, np.ndarray, list, Index)):
            self._setitem_array(key, value)
        elif isinstance(key, DataFrame):
            self._setitem_frame(key, value)
        else:
            # set column
            self._set_item(key, value)

    def _setitem_slice(self, key, value):
        self._check_setitem_copy()
        self.loc._setitem_with_indexer(key, value)

    def _setitem_array(self, key, value):
        # also raises Exception if object array with NA values
        if com.is_bool_indexer(key):
            if len(key) != len(self.index):
                raise ValueError('Item wrong length %d instead of %d!' %
                                 (len(key), len(self.index)))
            key = check_bool_indexer(self.index, key)
            indexer = key.nonzero()[0]
            self._check_setitem_copy()
            self.loc._setitem_with_indexer(indexer, value)
        else:
            if isinstance(value, DataFrame):
                if len(value.columns) != len(key):
                    raise ValueError('Columns must be same length as key')
                for k1, k2 in zip(key, value.columns):
                    self[k1] = value[k2]
            else:
                indexer = self.loc._convert_to_indexer(key, axis=1)
                self._check_setitem_copy()
                self.loc._setitem_with_indexer((slice(None), indexer), value)

    def _setitem_frame(self, key, value):
        # support boolean setting with DataFrame input, e.g.
        # df[df > df2] = 0
        if key.values.size and not is_bool_dtype(key.values):
            raise TypeError('Must pass DataFrame with boolean values only')

        self._check_inplace_setting(value)
        self._check_setitem_copy()
        self._where(-key, value, inplace=True)

    def _ensure_valid_index(self, value):
        """
        ensure that if we don't have an index, that we can create one from the
        passed value
        """
        # GH5632, make sure that we are a Series convertible
        if not len(self.index) and is_list_like(value):
            try:
                value = Series(value)
            except:
                raise ValueError('Cannot set a frame with no defined index '
                                 'and a value that cannot be converted to a '
                                 'Series')

            self._data = self._data.reindex_axis(value.index.copy(), axis=1,
                                                 fill_value=np.nan)

    def _set_item(self, key, value):
        """
        Add series to DataFrame in specified column.

        If series is a numpy-array (not a Series/TimeSeries), it must be the
        same length as the DataFrames index or an error will be thrown.

        Series/TimeSeries will be conformed to the DataFrames index to
        ensure homogeneity.
        """

        self._ensure_valid_index(value)
        value = self._sanitize_column(key, value)
        NDFrame._set_item(self, key, value)

        # check if we are modifying a copy
        # try to set first as we want an invalid
        # value exception to occur first
        if len(self):
            self._check_setitem_copy()

    def insert(self, loc, column, value, allow_duplicates=False):
        """
        Insert column into DataFrame at specified location.

        If `allow_duplicates` is False, raises Exception if column
        is already contained in the DataFrame.

        Parameters
        ----------
        loc : int
            Must have 0 <= loc <= len(columns)
        column : object
        value : scalar, Series, or array-like
        """
        self._ensure_valid_index(value)
        value = self._sanitize_column(column, value, broadcast=False)
        self._data.insert(loc, column, value,
                          allow_duplicates=allow_duplicates)

    def assign(self, **kwargs):
        """
        Assign new columns to a DataFrame, returning a new object
        (a copy) with all the original columns in addition to the new ones.

        .. versionadded:: 0.16.0

        Parameters
        ----------
        kwargs : keyword, value pairs
            keywords are the column names. If the values are
            callable, they are computed on the DataFrame and
            assigned to the new columns. The callable must not
            change input DataFrame (though pandas doesn't check it).
            If the values are not callable, (e.g. a Series, scalar, or array),
            they are simply assigned.

        Returns
        -------
        df : DataFrame
            A new DataFrame with the new columns in addition to
            all the existing columns.

        Notes
        -----
        Since ``kwargs`` is a dictionary, the order of your
        arguments may not be preserved. To make things predicatable,
        the columns are inserted in alphabetical order, at the end of
        your DataFrame. Assigning multiple columns within the same
        ``assign`` is possible, but you cannot reference other columns
        created within the same ``assign`` call.

        Examples
        --------
        >>> df = DataFrame({'A': range(1, 11), 'B': np.random.randn(10)})

        Where the value is a callable, evaluated on `df`:

        >>> df.assign(ln_A = lambda x: np.log(x.A))
            A         B      ln_A
        0   1  0.426905  0.000000
        1   2 -0.780949  0.693147
        2   3 -0.418711  1.098612
        3   4 -0.269708  1.386294
        4   5 -0.274002  1.609438
        5   6 -0.500792  1.791759
        6   7  1.649697  1.945910
        7   8 -1.495604  2.079442
        8   9  0.549296  2.197225
        9  10 -0.758542  2.302585

        Where the value already exists and is inserted:

        >>> newcol = np.log(df['A'])
        >>> df.assign(ln_A=newcol)
            A         B      ln_A
        0   1  0.426905  0.000000
        1   2 -0.780949  0.693147
        2   3 -0.418711  1.098612
        3   4 -0.269708  1.386294
        4   5 -0.274002  1.609438
        5   6 -0.500792  1.791759
        6   7  1.649697  1.945910
        7   8 -1.495604  2.079442
        8   9  0.549296  2.197225
        9  10 -0.758542  2.302585
        """
        data = self.copy()

        # do all calculations first...
        results = {}
        for k, v in kwargs.items():
            results[k] = com._apply_if_callable(v, data)

        # ... and then assign
        for k, v in sorted(results.items()):
            data[k] = v

        return data

    def _sanitize_column(self, key, value, broadcast=True):
        """
        Ensures new columns (which go into the BlockManager as new blocks) are
        always copied and converted into an array.

        Parameters
        ----------
        key : object
        value : scalar, Series, or array-like
        broadcast : bool, default True
            If ``key`` matches multiple duplicate column names in the
            DataFrame, this parameter indicates whether ``value`` should be
            tiled so that the returned array contains a (duplicated) column for
            each occurrence of the key. If False, ``value`` will not be tiled.

        Returns
        -------
        sanitized_column : numpy-array
        """

        def reindexer(value):
            # reindex if necessary

            if value.index.equals(self.index) or not len(self.index):
                value = value._values.copy()
            else:

                # GH 4107
                try:
                    value = value.reindex(self.index)._values
                except Exception as e:

                    # duplicate axis
                    if not value.index.is_unique:
                        raise e

                    # other
                    raise TypeError('incompatible index of inserted column '
                                    'with frame index')
            return value

        if isinstance(value, Series):
            value = reindexer(value)

        elif isinstance(value, DataFrame):
            # align right-hand-side columns if self.columns
            # is multi-index and self[key] is a sub-frame
            if isinstance(self.columns, MultiIndex) and key in self.columns:
                loc = self.columns.get_loc(key)
                if isinstance(loc, (slice, Series, np.ndarray, Index)):
                    cols = maybe_droplevels(self.columns[loc], key)
                    if len(cols) and not cols.equals(value.columns):
                        value = value.reindex_axis(cols, axis=1)
            # now align rows
            value = reindexer(value).T

        elif isinstance(value, Categorical):
            value = value.copy()

        elif isinstance(value, Index) or is_sequence(value):
            from pandas.core.series import _sanitize_index

            # turn me into an ndarray
            value = _sanitize_index(value, self.index, copy=False)
            if not isinstance(value, (np.ndarray, Index)):
                if isinstance(value, list) and len(value) > 0:
                    value = _possibly_convert_platform(value)
                else:
                    value = com._asarray_tuplesafe(value)
            elif value.ndim == 2:
                value = value.copy().T
            elif isinstance(value, Index):
                value = value.copy(deep=True)
            else:
                value = value.copy()

            # possibly infer to datetimelike
            if is_object_dtype(value.dtype):
                value = _possibly_infer_to_datetimelike(value)

        else:
            # upcast the scalar
            dtype, value = _infer_dtype_from_scalar(value)
            value = np.repeat(value, len(self.index)).astype(dtype)
            value = _possibly_cast_to_datetime(value, dtype)

        # return internal types directly
        if is_extension_type(value):
            return value

        # broadcast across multiple columns if necessary
        if broadcast and key in self.columns and value.ndim == 1:
            if (not self.columns.is_unique or
                    isinstance(self.columns, MultiIndex)):
                existing_piece = self[key]
                if isinstance(existing_piece, DataFrame):
                    value = np.tile(value, (len(existing_piece.columns), 1))

        return np.atleast_2d(np.asarray(value))

    @property
    def _series(self):
        result = {}
        for idx, item in enumerate(self.columns):
            result[item] = Series(self._data.iget(idx), index=self.index,
                                  name=item)
        return result

    def lookup(self, row_labels, col_labels):
        """Label-based "fancy indexing" function for DataFrame.
        Given equal-length arrays of row and column labels, return an
        array of the values corresponding to each (row, col) pair.

        Parameters
        ----------
        row_labels : sequence
            The row labels to use for lookup
        col_labels : sequence
            The column labels to use for lookup

        Notes
        -----
        Akin to::

            result = []
            for row, col in zip(row_labels, col_labels):
                result.append(df.get_value(row, col))

        Examples
        --------
        values : ndarray
            The found values

        """
        n = len(row_labels)
        if n != len(col_labels):
            raise ValueError('Row labels must have same size as column labels')

        thresh = 1000
        if not self._is_mixed_type or n > thresh:
            values = self.values
            ridx = self.index.get_indexer(row_labels)
            cidx = self.columns.get_indexer(col_labels)
            if (ridx == -1).any():
                raise KeyError('One or more row labels was not found')
            if (cidx == -1).any():
                raise KeyError('One or more column labels was not found')
            flat_index = ridx * len(self.columns) + cidx
            result = values.flat[flat_index]
        else:
            result = np.empty(n, dtype='O')
            for i, (r, c) in enumerate(zip(row_labels, col_labels)):
                result[i] = self.get_value(r, c)

        if is_object_dtype(result):
            result = lib.maybe_convert_objects(result)

        return result

    # ----------------------------------------------------------------------
    # Reindexing and alignment

    def _reindex_axes(self, axes, level, limit, tolerance, method, fill_value,
                      copy):
        frame = self

        columns = axes['columns']
        if columns is not None:
            frame = frame._reindex_columns(columns, method, copy, level,
                                           fill_value, limit, tolerance)

        index = axes['index']
        if index is not None:
            frame = frame._reindex_index(index, method, copy, level,
                                         fill_value, limit, tolerance)

        return frame

    def _reindex_index(self, new_index, method, copy, level, fill_value=NA,
                       limit=None, tolerance=None):
        new_index, indexer = self.index.reindex(new_index, method=method,
                                                level=level, limit=limit,
                                                tolerance=tolerance)
        return self._reindex_with_indexers({0: [new_index, indexer]},
                                           copy=copy, fill_value=fill_value,
                                           allow_dups=False)

    def _reindex_columns(self, new_columns, method, copy, level, fill_value=NA,
                         limit=None, tolerance=None):
        new_columns, indexer = self.columns.reindex(new_columns, method=method,
                                                    level=level, limit=limit,
                                                    tolerance=tolerance)
        return self._reindex_with_indexers({1: [new_columns, indexer]},
                                           copy=copy, fill_value=fill_value,
                                           allow_dups=False)

    def _reindex_multi(self, axes, copy, fill_value):
        """ we are guaranteed non-Nones in the axes! """

        new_index, row_indexer = self.index.reindex(axes['index'])
        new_columns, col_indexer = self.columns.reindex(axes['columns'])

        if row_indexer is not None and col_indexer is not None:
            indexer = row_indexer, col_indexer
            new_values = algorithms.take_2d_multi(self.values, indexer,
                                                  fill_value=fill_value)
            return self._constructor(new_values, index=new_index,
                                     columns=new_columns)
        else:
            return self._reindex_with_indexers({0: [new_index, row_indexer],
                                                1: [new_columns, col_indexer]},
                                               copy=copy,
                                               fill_value=fill_value)

    @Appender(_shared_docs['align'] % _shared_doc_kwargs)
    def align(self, other, join='outer', axis=None, level=None, copy=True,
              fill_value=None, method=None, limit=None, fill_axis=0,
              broadcast_axis=None):
        return super(DataFrame, self).align(other, join=join, axis=axis,
                                            level=level, copy=copy,
                                            fill_value=fill_value,
                                            method=method, limit=limit,
                                            fill_axis=fill_axis,
                                            broadcast_axis=broadcast_axis)

    @Appender(_shared_docs['reindex'] % _shared_doc_kwargs)
    def reindex(self, index=None, columns=None, **kwargs):
        return super(DataFrame, self).reindex(index=index, columns=columns,
                                              **kwargs)

    @Appender(_shared_docs['reindex_axis'] % _shared_doc_kwargs)
    def reindex_axis(self, labels, axis=0, method=None, level=None, copy=True,
                     limit=None, fill_value=np.nan):
        return super(DataFrame,
                     self).reindex_axis(labels=labels, axis=axis,
                                        method=method, level=level, copy=copy,
                                        limit=limit, fill_value=fill_value)

    @Appender(_shared_docs['rename'] % _shared_doc_kwargs)
    def rename(self, index=None, columns=None, **kwargs):
        return super(DataFrame, self).rename(index=index, columns=columns,
                                             **kwargs)

    @Appender(_shared_docs['fillna'] % _shared_doc_kwargs)
    def fillna(self, value=None, method=None, axis=None, inplace=False,
               limit=None, downcast=None, **kwargs):
        return super(DataFrame,
                     self).fillna(value=value, method=method, axis=axis,
                                  inplace=inplace, limit=limit,
                                  downcast=downcast, **kwargs)

    @Appender(_shared_docs['shift'] % _shared_doc_kwargs)
    def shift(self, periods=1, freq=None, axis=0):
        return super(DataFrame, self).shift(periods=periods, freq=freq,
                                            axis=axis)

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
        inplace = validate_bool_kwarg(inplace, 'inplace')
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
                    arrays.append(self.index._get_level_values(i))
            else:
                arrays.append(self.index)

        to_remove = []
        for col in keys:
            if isinstance(col, MultiIndex):
                # append all but the last column so we don't have to modify
                # the end of this loop
                for n in range(col.nlevels - 1):
                    arrays.append(col._get_level_values(n))

                level = col._get_level_values(col.nlevels - 1)
                names.extend(col.names)
            elif isinstance(col, Series):
                level = col._values
                names.append(col.name)
            elif isinstance(col, Index):
                level = col
                names.append(col.name)
            elif isinstance(col, (list, np.ndarray, Index)):
                level = col
                names.append(None)
            else:
                level = frame[col]._values
                names.append(col)
                if drop:
                    to_remove.append(col)
            arrays.append(level)

        index = MultiIndex.from_arrays(arrays, names=names)

        if verify_integrity and not index.is_unique:
            duplicates = index.get_duplicates()
            raise ValueError('Index has duplicate keys: %s' % duplicates)

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
        inplace = validate_bool_kwarg(inplace, 'inplace')
        if inplace:
            new_obj = self
        else:
            new_obj = self.copy()

        def _maybe_casted_values(index, labels=None):
            if isinstance(index, PeriodIndex):
                values = index.asobject.values
            elif isinstance(index, DatetimeIndex) and index.tz is not None:
                values = index
            else:
                values = index.values
                if values.dtype == np.object_:
                    values = lib.maybe_convert_objects(values)

            # if we have the labels, extract the values with a mask
            if labels is not None:
                mask = labels == -1

                # we can have situations where the whole mask is -1,
                # meaning there is nothing found in labels, so make all nan's
                if mask.all():
                    values = np.empty(len(mask))
                    values.fill(np.nan)
                else:
                    values = values.take(labels)
                    if mask.any():
                        values, changed = _maybe_upcast_putmask(values, mask,
                                                                np.nan)
            return values

        new_index = _default_index(len(new_obj))
        if isinstance(self.index, MultiIndex):
            if level is not None:
                if not isinstance(level, (tuple, list)):
                    level = [level]
                level = [self.index._get_level_number(lev) for lev in level]
                if len(level) < len(self.index.levels):
                    new_index = self.index.droplevel(level)

            if not drop:
                names = self.index.names
                zipped = lzip(self.index.levels, self.index.labels)

                multi_col = isinstance(self.columns, MultiIndex)
                for i, (lev, lab) in reversed(list(enumerate(zipped))):
                    col_name = names[i]
                    if col_name is None:
                        col_name = 'level_%d' % i

                    if multi_col:
                        if col_fill is None:
                            col_name = tuple([col_name] * self.columns.nlevels)
                        else:
                            name_lst = [col_fill] * self.columns.nlevels
                            lev_num = self.columns._get_level_number(col_level)
                            name_lst[lev_num] = col_name
                            col_name = tuple(name_lst)

                    # to ndarray and maybe infer different dtype
                    level_values = _maybe_casted_values(lev, lab)
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
            values = _maybe_casted_values(self.index)
            new_obj.insert(0, name, values)

        new_obj.index = new_index
        if not inplace:
            return new_obj

    # ----------------------------------------------------------------------
    # Reindex-based selection methods

    def dropna(self, axis=0, how='any', thresh=None, subset=None,
               inplace=False):
        """
        Return object with labels on given axis omitted where alternately any
        or all of the data are missing

        Parameters
        ----------
        axis : {0 or 'index', 1 or 'columns'}, or tuple/list thereof
            Pass tuple or list to drop on multiple axes
        how : {'any', 'all'}
            * any : if any NA values are present, drop that label
            * all : if all values are NA, drop that label
        thresh : int, default None
            int value : require that many non-NA values
        subset : array-like
            Labels along other axis to consider, e.g. if you are dropping rows
            these would be a list of columns to include
        inplace : boolean, default False
            If True, do operation inplace and return None.

        Returns
        -------
        dropped : DataFrame

        Examples
        --------
        >>> df = pd.DataFrame([[np.nan, 2, np.nan, 0], [3, 4, np.nan, 1],
                               [np.nan, np.nan, np.nan, 5]],
                              columns=list('ABCD'))
        >>> df
        >>>     A    B   C    D
            0   NaN  2.0 NaN  0
            1   3.0  4.0 NaN  1
            2   NaN  NaN NaN  5

        Drop the columns where all elements are nan

        >>> df.dropna(axis=1, how='all')
        >>>     A   B    D
            0  NaN  2.0  0
            1  3.0  4.0  1
            2  NaN  NaN  5
        >>> df.dropna(axis=1, how='any')

        Drop the columns where any of the elements is nan

        >>>    D
            0  0
            1  1
            2  5
        Drop the rows where all of the elements is nan
        (there is no row to drop, so df stays the same)

        >>> df.dropna(axis=0, how='all')
        >>>     A    B   C    D
            0   NaN  2.0 NaN  0
            1   3.0  4.0 NaN  1
            2   NaN  NaN NaN  5

        Keep only the rows with at least 2 non-na values
        
        >>> df.dropna(thresh=2)
        >>>    A    B     C    D
            0  NaN  2.0   NaN  0
            1  3.0  4.0   NaN  1


DataFrame._setup_axes(['index', 'columns'], info_axis=1, stat_axis=0,
                      axes_are_reversed=True, aliases={'rows': 0})
DataFrame._add_numeric_operations()
DataFrame._add_series_or_dataframe_operations()

_EMPTY_SERIES = Series([])


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
    elif len(data) > 0:
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
                indexes.append(list(v.keys()))
            elif is_list_like(v) and getattr(v, 'ndim', 1) == 1:
                have_raw_arrays = True
                raw_lengths.append(len(v))

        if not indexes and not raw_lengths:
            raise ValueError('If using all scalar values, you must pass'
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
                    msg = ('array length %d does not match index length %d' %
                           (lengths[0], len(index)))
                    raise ValueError(msg)
            else:
                index = _default_index(lengths[0])

    return _ensure_index(index)


def _prep_ndarray(values, copy=True):
    if not isinstance(values, (np.ndarray, Series, Index)):
        if len(values) == 0:
            return np.empty((0, 0), dtype=object)

        def convert(v):
            return _possibly_convert_platform(v)

        # we could have a 1-dim or 2-dim list here
        # this is equiv of np.asarray, but does object conversion
        # and platform dtype preservation
        try:
            if is_list_like(values[0]) or hasattr(values[0], 'len'):
                values = np.array([convert(v) for v in values])
            else:
                values = convert(values)
        except:
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


def _to_arrays(data, columns, coerce_float=False, dtype=None):
    """
    Return list of arrays, columns
    """
    if isinstance(data, DataFrame):
        if columns is not None:
            arrays = [data._ixs(i, axis=1).values
                      for i, col in enumerate(data.columns) if col in columns]
        else:
            columns = data.columns
            arrays = [data._ixs(i, axis=1).values for i in range(len(columns))]

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
                                       coerce_float=coerce_float, dtype=dtype)
    elif isinstance(data[0], Series):
        return _list_of_series_to_arrays(data, columns,
                                         coerce_float=coerce_float,
                                         dtype=dtype)
    elif isinstance(data[0], Categorical):
        if columns is None:
            columns = _default_index(len(data))
        return data, columns
    elif (isinstance(data, (np.ndarray, Series, Index)) and
          data.dtype.names is not None):

        columns = list(data.dtype.names)
        arrays = [data[k] for k in columns]
        return arrays, columns
    else:
        # last ditch effort
        data = lmap(tuple, data)
        return _list_to_arrays(data, columns, coerce_float=coerce_float,
                               dtype=dtype)


def _masked_rec_array_to_mgr(data, index, columns, dtype, copy):
    """ extract from a masked rec array and create the manager """

    # essentially process a record array then fill it
    fill_value = data.fill_value
    fdata = ma.getdata(data)
    if index is None:
        index = _get_names_from_index(fdata)
        if index is None:
            index = _default_index(len(data))
    index = _ensure_index(index)

    if columns is not None:
        columns = _ensure_index(columns)
    arrays, arr_columns = _to_arrays(fdata, columns)

    # fill if needed
    new_arrays = []
    for fv, arr, col in zip(fill_value, arrays, arr_columns):
        mask = ma.getmaskarray(data[col])
        if mask.any():
            arr, fv = _maybe_upcast(arr, fill_value=fv, copy=True)
            arr[mask] = fv
        new_arrays.append(arr)

    # create the manager
    arrays, arr_columns = _reorder_arrays(new_arrays, arr_columns, columns)
    if columns is None:
        columns = arr_columns

    mgr = _arrays_to_mgr(arrays, arr_columns, index, columns)

    if copy:
        mgr = mgr.copy()
    return mgr


def _reorder_arrays(arrays, arr_columns, columns):
    # reorder according to the columns
    if (columns is not None and len(columns) and arr_columns is not None and
            len(arr_columns)):
        indexer = _ensure_index(arr_columns).get_indexer(columns)
        arr_columns = _ensure_index([arr_columns[i] for i in indexer])
        arrays = [arrays[i] for i in indexer]
    return arrays, arr_columns


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
        columns = _get_combined_index([
            s.index for s in data if getattr(s, 'index', None) is not None
        ])

    indexer_cache = {}

    aligned_values = []
    for s in data:
        index = getattr(s, 'index', None)
        if index is None:
            index = _default_index(len(s))

        if id(index) in indexer_cache:
            indexer = indexer_cache[id(index)]
        else:
            indexer = indexer_cache[id(index)] = index.get_indexer(columns)

        values = _values_from_object(s)
        aligned_values.append(algorithms.take_1d(values, indexer))

    values = np.vstack(aligned_values)

    if values.dtype == np.object_:
        content = list(values.T)
        return _convert_object_array(content, columns, dtype=dtype,
                                     coerce_float=coerce_float)
    else:
        return values.T, columns


def _list_of_dict_to_arrays(data, columns, coerce_float=False, dtype=None):
    if columns is None:
        gen = (list(x.keys()) for x in data)
        sort = not any(isinstance(d, OrderedDict) for d in data)
        columns = lib.fast_unique_multiple_list_gen(gen, sort=sort)

    # assure that they are of the base dict class and not of derived
    # classes
    data = [(type(d) is dict) and d or dict(d) for d in data]

    content = list(lib.dicts_to_array(data, list(columns)).T)
    return _convert_object_array(content, columns, dtype=dtype,
                                 coerce_float=coerce_float)


def _convert_object_array(content, columns, coerce_float=False, dtype=None):
    if columns is None:
        columns = _default_index(len(content))
    else:
        if len(columns) != len(content):  # pragma: no cover
            # caller's responsibility to check for this...
            raise AssertionError('%d columns passed, passed data had %s '
                                 'columns' % (len(columns), len(content)))

    # provide soft conversion of object dtypes
    def convert(arr):
        if dtype != object and dtype != np.object:
            arr = lib.maybe_convert_objects(arr, try_float=coerce_float)
            arr = _possibly_cast_to_datetime(arr, dtype)
        return arr

    arrays = [convert(arr) for arr in content]

    return arrays, columns


def _get_names_from_index(data):
    has_some_name = any([getattr(s, 'name', None) is not None for s in data])
    if not has_some_name:
        return _default_index(len(data))

    index = lrange(len(data))
    count = 0
    for i, s in enumerate(data):
        n = getattr(s, 'name', None)
        if n is not None:
            index[i] = n
        else:
            index[i] = 'Unnamed %d' % count
            count += 1

    return index


def _homogenize(data, index, dtype=None):
    from pandas.core.series import _sanitize_array

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

                if isinstance(index, (DatetimeIndex, TimedeltaIndex)):
                    v = _dict_compat(v)
                else:
                    v = dict(v)
                v = lib.fast_multiget(v, oindex.values, default=NA)
            v = _sanitize_array(v, index, dtype=dtype, copy=False,
                                raise_cast_failure=False)

        homogenized.append(v)

    return homogenized


def _from_nested_dict(data):
    # TODO: this should be seriously cythonized
    new_data = OrderedDict()
    for index, s in compat.iteritems(data):
        for col, v in compat.iteritems(s):
            new_data[col] = new_data.get(col, OrderedDict())
            new_data[col][index] = v
    return new_data


def _put_str(s, space):
    return ('%s' % s)[:space].ljust(space)


# ----------------------------------------------------------------------
# Add plotting methods to DataFrame
DataFrame.plot = base.AccessorProperty(gfx.FramePlotMethods,
                                       gfx.FramePlotMethods)
DataFrame.hist = gfx.hist_frame


@Appender(_shared_docs['boxplot'] % _shared_doc_kwargs)
def boxplot(self, column=None, by=None, ax=None, fontsize=None, rot=0,
            grid=True, figsize=None, layout=None, return_type=None, **kwds):
    import pandas.tools.plotting as plots
    import matplotlib.pyplot as plt
    ax = plots.boxplot(self, column=column, by=by, ax=ax, fontsize=fontsize,
                       grid=grid, rot=rot, figsize=figsize, layout=layout,
                       return_type=return_type, **kwds)
    plt.draw_if_interactive()
    return ax


DataFrame.boxplot = boxplot

ops.add_flex_arithmetic_methods(DataFrame, **ops.frame_flex_funcs)
ops.add_special_arithmetic_methods(DataFrame, **ops.frame_special_funcs)
