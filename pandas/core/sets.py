import numpy as np

from functools import reduce
from operator import __or__, __xor__, __and__, __sub__

from pandas.core.dtypes.generic import ABCSeries
from pandas.core.dtypes.missing import isna
from pandas.core.dtypes.common import is_list_like

from pandas.core.base import NoNewAttributesMixin
from pandas.util._decorators import Appender
import pandas.compat as compat

_shared_docs = dict()


class SetMethods(NoNewAttributesMixin):
    """
    Vectorized set functions for Series. NAs get turned to empty sets by
    default - this behavior can be changed by using the 'fill_value'-parameter.
    All methods have an 'errors'-parameter that determines how values are
    converted to sets.

    Examples
    --------
    >>> s.set.union()
    >>> s.set.intersect()
    """

    def __init__(self, data):
        self._data = data
        self._freeze()

    @staticmethod
    def _validate(data, errors='raise', fill_value=None):
        """
        TODO
        """

        # signature following GH13877
        if not isinstance(data, ABCSeries):
            raise ValueError("Must pass Series for validating inputs of set "
                             "accessor operations")

        if fill_value is not None and fill_value is not np.nan:
            err_str = ("The parameter 'fill_value' must be list-like!")
            if not is_list_like(fill_value):
                raise ValueError(err_str)
            fill_value = set(fill_value)

        data = data.copy()  # avoid changing original input
        na_mask = data.isna()

        if errors == 'raise':
            forbidden = ~data.loc[~na_mask].map(is_list_like)
            if forbidden.any():
                raise ValueError("By default, can only use .set accessor with "
                                 "values that can be mapped to sets. For more "
                                 "permissive error-handling, set 'errors'="
                                 "'ignore'|'coerce'|'wrap'.")
        elif errors == 'ignore':
            ignore = ~data.loc[~na_mask].map(is_list_like)
            # everything that's not list-like gets set to na
            na_mask.loc[~na_mask] = ignore
        elif errors == 'coerce':
            permitted = lambda x: (isinstance(x, compat.string_types)
                                   or is_list_like(x))
            # everything that's not a string or container gets set to na
            na_mask.loc[~na_mask] = ~data.loc[~na_mask].map(permitted)
        elif errors == 'wrap':
            singletons = ~na_mask & ~data.map(is_list_like)
            data.loc[singletons] = data.loc[singletons].map(lambda x: [x])
        elif errors == 'skip':
            pass
        else:
            raise ValueError("Received illegal value for parameter 'errors'; "
                             "allowed values are {'raise'|'ignore'|"
                             "'coerce'|'wrap'|'skip'}")

        if errors != 'skip':
            data.loc[na_mask] = np.nan
            # everything else gets mapped to sets
            data.loc[~na_mask] = data.loc[~na_mask].map(set)

        # cannot use fillna due to GH21329
        if fill_value is not None and na_mask.any():
            data.loc[na_mask] = [fill_value] * na_mask.sum()

        return data

    def _wrap_result(self, result, name=None, index=None):
        """
        TODO
        """

        from pandas import Series

        if name is None:
            name = getattr(result, 'name', None)
        if name is None:
            name = self._data.name

        if not hasattr(result, 'ndim') or not hasattr(result, 'dtype'):
            return result
        assert result.ndim < 3

        index = self._data.index if index is None else index
        return Series(result, name=name, index=index)

    def _get_series_list(self, others):
        """
        Auxiliary function for set-accessor functions. Turn potentially mixed
        input into a list of Series.

        Parameters
        ----------
        others : Series, DataFrame, np.ndarray, or list-like of objects that
            are either Series or np.ndarray (1-dim). If it is a list-like that
            *only* contains scalar values, this list-like object will be
            broadcast to every element of a Series with the same index as the
            calling Series.

        Returns
        -------
        list : others transformed into list of Series
        """

        from pandas import Series, DataFrame

        idx = self._data.index

        err_msg = ('others must be Series, DataFrame, np.ndarrary or '
                   'list-like (containing either only scalar values, or only '
                   'objects of type Series/np.ndarray)!')

        # np.ndarray inherits the index `idx` of the calling Series - i.e. must
        # have matching length. Series/DataFrame keep their own index.
        # List-likes must contain only Series or 1-dim np.ndarray
        if isinstance(others, Series):
            return [others]
        elif isinstance(others, DataFrame):
            return [others[x] for x in others]
        elif isinstance(others, np.ndarray) and others.ndim == 1:
            return [Series(others, index=idx)]
        elif isinstance(others, np.ndarray) and others.ndim == 2:
            others = DataFrame(others, index=idx)
            return [others[x] for x in others]
        elif is_list_like(others):
            others = list(others)  # ensure iterators do not get read twice etc

            # in case of list-like `others`, all elements must be either be
            # scalar or Series/np.ndarray
            if all(not is_list_like(x) for x in others):  # True if empty
                # in this case, we broadcast others to every element of a new
                # Series with the same index as the caller
                return [Series([others] * len(idx), index=idx)]

            check = lambda x: (isinstance(x, Series)
                               or (isinstance(x, np.ndarray) and x.ndim == 1))
            if all(check(x) for x in others):
                los = []
                # iterate through list and append list of series for each
                # element (which we check to be one-dimensional)
                while others:
                    nxt = others.pop(0)  # Series or np.ndarray by the above
                    los = los + self._get_series_list(nxt)
                return los
        raise TypeError(err_msg)

    def _apply_op(self, others, operator, errors, fill_value, join):
        """
        TODO
        """

        from pandas import concat

        data = self._validate(self._data, errors, fill_value)

        # concatenate Series/Index with itself if no "others"
        if others is None:
            return reduce(operator, data.dropna())

        try:
            # turn others into list of series -- necessary for concat/align
            others = self._get_series_list(others)
        except ValueError:  # do not catch TypeError raised by _get_series_list
            raise ValueError('If `others` contains arrays, these must all be '
                             'of the same length as the calling Series.')
        # check if all series are legal for set ops; raise/convert otherwise
        others = [self._validate(x, errors, fill_value) for x in others]

        # Need to add keys for uniqueness in case of duplicate columns
        others = concat(others, axis=1,
                        join=(join if join == 'inner' else 'outer'),
                        keys=range(len(others)))
        data, others = data.align(others, join=join)
        allcols = [data] + [others[x] for x in others]  # again list of Series

        # if alignment introduced NaNs anywhere, need to re-apply fill_value
        if fill_value is not None and any(x.isna().any() for x in allcols):
            allcols = [self._validate(x, 'skip', fill_value) for x in allcols]

        result = self._apply_op_core(allcols, operator)
        index = others.index if join == 'right' else data.index
        return self._wrap_result(result, index=index)

    def _apply_op_core(self, list_of_series, operator):
        """
        TODO
        """
        # list_of_series: must be aligned already!
        masks = np.array([isna(x).values for x in list_of_series])
        na_mask = np.logical_or.reduce(masks, axis=0)
        result = np.empty(len(na_mask), dtype=object)
        np.putmask(result, na_mask, np.nan)

        # apply operator over columns left-to-right; everything aligned already
        result[~na_mask] = reduce(operator,
                                  [x.values[~na_mask] for x in list_of_series])
        return result

    _shared_docs['set_ops'] = ("""
    Calculate %(op)s for Series.

    If `others` is specified, this method applies the %(op)s per element. If
    `others` is not passed, then the %(op)s is applied to the elements in the
    Series%(add)s.

    Parameters
    ----------
    others : Series, DataFrame, np.ndarray, or list-like, default None
        np.ndarray (one- or two-dimensional) must have the same length as the
        calling Series; Series and DataFrame get matched on index and therefore
        do not have to match in length.

        If `others` is a list-like, it may contain either:

        - Only Series or np.ndarray. The latter must have the same length as
          the calling Series
        - Only scalars. In this case, this list-like object will be used as the
          right-hand side of the %(op)s for all elements of the calling Series.

        If `others` is None, the method applies the %(op)s%(add)s to the
        elements of the calling Series.
    join : {'left', 'right', 'outer', 'inner'}, default 'left'
        Determines the join-style between the calling Series and any Series or
        DataFrame in `others` (np.ndarrays need to match the length of the
        calling Series). To disable alignment, use `.values` on any Series or
        DataFrame in `others`.
    errors : {'raise', 'ignore', 'coerce', 'wrap', 'skip'}, default 'raise'
        Determines how values that are not sets are treated, both in the
        calling Series, as well as any column in `others`. All options ignore
        missing values, and all options except 'raise' and 'skip' will set
        elements that they cannot map to `np.nan` - these values can be further
        processed using the `fill_value`-parameter.

        - 'raise': Raise error for any element that cannot be unambiguously
          mapped to a set (including strings).
        - 'ignore': Ignore all elements that cannot be unambiguously mapped to
          a set (including strings).
        - 'coerce': Forcefully map everything possible to a set. In particular,
          strings will be mapped to the set of their characters.
        - 'wrap': Maps list-likes to sets, and wraps all other elements
          (including strings; except missing values) into singleton sets.
        - 'skip': Do not run any checks or conversions to `set`, if
          performance is critical (`fill_values` will work as usual). In this
          case, it is up to the user that all non-null elements are compatible
          with the respective `numpy` set-methods.
    fill_value : list-like, default None
        Value to use for missing values in the calling Series and any column in
        `others`.

    Returns
    -------
    result : set or Series/Index of objects
        If `others` is None, `set` is returned, otherwise a `Series` of objects
        is returned.

    See Also
    --------
    %(also)s
    Examples
    --------
    If `others` is not specified, the operation will be applied%(add)s to all
    elements of the Series.

    >>> s = pd.Series([{1, 2}, {2, 4}, {3, 1}])
    >>> s
    0    {1, 2}
    1    {2, 4}
    2    {1, 3}
    dtype: object
    %(ex_no_others)s
    If `others` is a Series (or np.ndarray of the correct length), the
    operation will be applied element-wise.

    >>> t = pd.Series([{2, 3}, {1, 2}, np.nan])
    >>> t
    0    {2, 3}
    1    {1, 2}
    2       NaN
    dtype: object
    %(ex_with_others)s
    By default, missing values in any of the input columns will remain missing
    in the result. To change this, use the `fill_value` parameter, which fills
    the columns *before* applying the operation.
    %(ex_with_fill)s
    Finally, if `others` is a list-like containing only scalar values, this
    list-like object will be used as the right-hand side of the %(op)s for all
    elements of the calling Series. (as in the corresponding `numpy` methods).

    >>> s
    0    {1, 2}
    1    {2, 4}
    2    {1, 3}
    dtype: object
    %(ex_scalar)s
    For more examples, see :ref:`here <set.accessor>`.
    """)

    also = '''
    intersect : Calculate intersection
    diff : Calculate set difference
    xor : Calculate symmetric set difference
    '''
    ex_no_others = '''>>>
    >>> s.set.union()
    {1, 2, 3, 4}
    '''
    ex_with_others = '''>>>
    >>> s.set.union(t)
    0    {1, 2, 3}
    1    {1, 2, 4}
    2          NaN
    dtype: object
    '''
    ex_with_fill = '''
    >>> s.set.union(t, fill_value=set())  # equivalent fill values: [], {}, ...
    0    {1, 2, 3}
    1    {1, 2, 4}
    2       {1, 3}
    dtype: object
    >>>
    >>> s.set.union(t, fill_value={1, 3, 5})
    0    {1, 2, 3}
    1    {1, 2, 4}
    2    {1, 3, 5}
    dtype: object
    '''
    ex_scalar = '''>>>
    >>> s.set.union({1})
    0       {1, 2}
    1    {1, 2, 4}
    2       {1, 3}
    dtype: object
    '''

    @Appender(_shared_docs['set_ops'] % {
        'op': 'union',
        'add': '',
        'also': also,
        'ex_no_others': ex_no_others,
        'ex_with_others': ex_with_others,
        'ex_with_fill': ex_with_fill,
        'ex_scalar': ex_scalar
    })
    def union(self, others=None, join='left', errors='raise', fill_value=None):
        return self._apply_op(others, __or__, errors, fill_value, join)

    also = '''
    union : Calculate union
    diff : Calculate set difference
    xor : Calculate symmetric set difference
    '''
    ex_no_others = '''>>>
    >>> s.set.intersect()
    set()
    '''
    ex_with_others = '''>>>
    >>> s.set.intersect(t)
    0    {2}
    1    {2}
    2    NaN
    dtype: object
    '''
    ex_with_fill = '''
    >>> s.set.intersect(t, fill_value=set())  # equiv. fill values: [], {}, ...
    0    {2}
    1    {2}
    2     {}
    dtype: object
    >>>
    >>> s.set.intersect(t, fill_value={1, 3, 5})
    0       {2}
    1       {2}
    2    {1, 3}
    dtype: object
    '''
    ex_scalar = '''>>>
    >>> s.set.intersect({1})
    0    {1}
    1     {}
    2    {1}
    dtype: object
    '''

    @Appender(_shared_docs['set_ops'] % {
        'op': 'intersection',
        'add': '',
        'also': also,
        'ex_no_others': ex_no_others,
        'ex_with_others': ex_with_others,
        'ex_with_fill': ex_with_fill,
        'ex_scalar': ex_scalar
    })
    def intersect(self, others=None, join='left',
                  errors='raise', fill_value=None):
        return self._apply_op(others, __and__, errors, fill_value, join)

    also = '''
    intersect : Calculate intersection
    union : Calculate union
    xor : Calculate symmetric set difference
    '''
    ex_no_others = '''>>>
    >>> s.set.diff()
    set()
    '''
    ex_with_others = '''>>>
    >>> s.set.diff(t)
    0    {1}
    1    {4}
    2    NaN
    dtype: object
    '''
    ex_with_fill = '''
    >>> s.set.diff(t, fill_value=set())  # equivalent fill values: [], {}, ...
    0       {1}
    1       {4}
    2    {1, 3}
    dtype: object
    >>>
    >>> s.set.diff(t, fill_value={1, 3, 5})
    0    {1}
    1    {4}
    2     {}
    dtype: object
    '''
    ex_scalar = '''>>>
    >>> s.set.diff({1})
    0       {2}
    1    {2, 4}
    2       {3}
    dtype: object
    '''

    @Appender(_shared_docs['set_ops'] % {
        'op': 'set difference',
        'add': ' sequentially',
        'also': also,
        'ex_no_others': ex_no_others,
        'ex_with_others': ex_with_others,
        'ex_with_fill': ex_with_fill,
        'ex_scalar': ex_scalar
    })
    def diff(self, others=None, join='left', errors='raise', fill_value=None):
        return self._apply_op(others, __sub__, errors, fill_value, join)

    also = '''
    intersect : Calculate intersection
    union : Calculate union
    diff : Calculate symmetric set difference
    '''
    ex_no_others = '''>>>
    >>> s.set.xor()
    {3, 4}
    '''
    ex_with_others = '''>>>
    >>> s.set.xor(t)
    0    {1, 3}
    1    {1, 4}
    2       NaN
    dtype: object
    '''
    ex_with_fill = '''
    >>> s.set.xor(t, fill_value=set())  # equivalent fill values: [], {}, ...
    0    {1, 3}
    1    {1, 4}
    2    {1, 3}
    dtype: object
    >>>
    >>> s.set.xor(t, fill_value={1, 3, 5})
    0    {1, 3}
    1    {1, 4}
    2       {5}
    dtype: object
    '''
    ex_scalar = '''>>>
    >>> s.set.xor({1})
    0          {2}
    1    {1, 2, 4}
    2          {3}
    dtype: object
    '''

    @Appender(_shared_docs['set_ops'] % {
        'op': 'symmetric set difference',
        'add': ' sequentially',
        'also': also,
        'ex_no_others': ex_no_others,
        'ex_with_others': ex_with_others,
        'ex_with_fill': ex_with_fill,
        'ex_scalar': ex_scalar
    })
    def xor(self, others=None, join='left', errors='raise', fill_value=None):
        return self._apply_op(others, __xor__, errors, fill_value, join)

    def len(self, errors='raise', fill_value=None):
        """
        TODO
        """

        from pandas import Series

        data = self._validate(self._data, errors, fill_value)
        na_mask = data.isna()
        result = Series(index=data.index)
        result.loc[~na_mask] = data.dropna().map(len)
        return result

    @classmethod
    def _make_accessor(cls, data):
        cls._validate(data)
        return cls(data)
