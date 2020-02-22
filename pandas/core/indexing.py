from typing import Hashable, List, Tuple, Union

import numpy as np

from pandas._libs.indexing import _NDFrameIndexerBase
from pandas._libs.lib import item_from_zerodim
from pandas.errors import AbstractMethodError
from pandas.util._decorators import Appender

from pandas.core.dtypes.common import (
    is_float,
    is_integer,
    is_iterator,
    is_list_like,
    is_numeric_dtype,
    is_object_dtype,
    is_scalar,
    is_sequence,
)
from pandas.core.dtypes.concat import concat_compat
from pandas.core.dtypes.generic import ABCDataFrame, ABCMultiIndex, ABCSeries
from pandas.core.dtypes.missing import _infer_fill_value, isna

import pandas.core.common as com
from pandas.core.indexers import (
    check_array_indexer,
    is_list_like_indexer,
    length_of_indexer,
)
from pandas.core.indexes.api import Index, InvalidIndexError

# "null slice"
_NS = slice(None, None)


# the public IndexSlicerMaker
class _IndexSlice:
    """
    Create an object to more easily perform multi-index slicing.

    See Also
    --------
    MultiIndex.remove_unused_levels : New MultiIndex with no unused levels.

    Notes
    -----
    See :ref:`Defined Levels <advanced.shown_levels>`
    for further info on slicing a MultiIndex.

    Examples
    --------
    >>> midx = pd.MultiIndex.from_product([['A0','A1'], ['B0','B1','B2','B3']])
    >>> columns = ['foo', 'bar']
    >>> dfmi = pd.DataFrame(np.arange(16).reshape((len(midx), len(columns))),
                            index=midx, columns=columns)

    Using the default slice command:

    >>> dfmi.loc[(slice(None), slice('B0', 'B1')), :]
               foo  bar
        A0 B0    0    1
           B1    2    3
        A1 B0    8    9
           B1   10   11

    Using the IndexSlice class for a more intuitive command:

    >>> idx = pd.IndexSlice
    >>> dfmi.loc[idx[:, 'B0':'B1'], :]
               foo  bar
        A0 B0    0    1
           B1    2    3
        A1 B0    8    9
           B1   10   11
    """

    def __getitem__(self, arg):
        return arg


IndexSlice = _IndexSlice()


class IndexingError(Exception):
    pass


class IndexingMixin:
    """
    Mixin for adding .loc/.iloc/.at/.iat to Datafames and Series.
    """

    @property
    def iloc(self) -> "_iLocIndexer":
        """
        Purely integer-location based indexing for selection by position.

        ``.iloc[]`` is primarily integer position based (from ``0`` to
        ``length-1`` of the axis), but may also be used with a boolean
        array.

        Allowed inputs are:

        - An integer, e.g. ``5``.
        - A list or array of integers, e.g. ``[4, 3, 0]``.
        - A slice object with ints, e.g. ``1:7``.
        - A boolean array.
        - A ``callable`` function with one argument (the calling Series or
          DataFrame) and that returns valid output for indexing (one of the above).
          This is useful in method chains, when you don't have a reference to the
          calling object, but would like to base your selection on some value.

        ``.iloc`` will raise ``IndexError`` if a requested indexer is
        out-of-bounds, except *slice* indexers which allow out-of-bounds
        indexing (this conforms with python/numpy *slice* semantics).

        See more at :ref:`Selection by Position <indexing.integer>`.

        See Also
        --------
        DataFrame.iat : Fast integer location scalar accessor.
        DataFrame.loc : Purely label-location based indexer for selection by label.
        Series.iloc : Purely integer-location based indexing for
                       selection by position.

        Examples
        --------
        >>> mydict = [{'a': 1, 'b': 2, 'c': 3, 'd': 4},
        ...           {'a': 100, 'b': 200, 'c': 300, 'd': 400},
        ...           {'a': 1000, 'b': 2000, 'c': 3000, 'd': 4000 }]
        >>> df = pd.DataFrame(mydict)
        >>> df
              a     b     c     d
        0     1     2     3     4
        1   100   200   300   400
        2  1000  2000  3000  4000

        **Indexing just the rows**

        With a scalar integer.

        >>> type(df.iloc[0])
        <class 'pandas.core.series.Series'>
        >>> df.iloc[0]
        a    1
        b    2
        c    3
        d    4
        Name: 0, dtype: int64

        With a list of integers.

        >>> df.iloc[[0]]
           a  b  c  d
        0  1  2  3  4
        >>> type(df.iloc[[0]])
        <class 'pandas.core.frame.DataFrame'>

        >>> df.iloc[[0, 1]]
             a    b    c    d
        0    1    2    3    4
        1  100  200  300  400

        With a `slice` object.

        >>> df.iloc[:3]
              a     b     c     d
        0     1     2     3     4
        1   100   200   300   400
        2  1000  2000  3000  4000

        With a boolean mask the same length as the index.

        >>> df.iloc[[True, False, True]]
              a     b     c     d
        0     1     2     3     4
        2  1000  2000  3000  4000

        With a callable, useful in method chains. The `x` passed
        to the ``lambda`` is the DataFrame being sliced. This selects
        the rows whose index label even.

        >>> df.iloc[lambda x: x.index % 2 == 0]
              a     b     c     d
        0     1     2     3     4
        2  1000  2000  3000  4000

        **Indexing both axes**

        You can mix the indexer types for the index and columns. Use ``:`` to
        select the entire axis.

        With scalar integers.

        >>> df.iloc[0, 1]
        2

        With lists of integers.

        >>> df.iloc[[0, 2], [1, 3]]
              b     d
        0     2     4
        2  2000  4000

        With `slice` objects.

        >>> df.iloc[1:3, 0:3]
              a     b     c
        1   100   200   300
        2  1000  2000  3000

        With a boolean array whose length matches the columns.

        >>> df.iloc[:, [True, False, True, False]]
              a     c
        0     1     3
        1   100   300
        2  1000  3000

        With a callable function that expects the Series or DataFrame.

        >>> df.iloc[:, lambda df: [0, 2]]
              a     c
        0     1     3
        1   100   300
        2  1000  3000
        """
        return _iLocIndexer("iloc", self)

    @property
    def loc(self) -> "_LocIndexer":
        """
        Access a group of rows and columns by label(s) or a boolean array.

        ``.loc[]`` is primarily label based, but may also be used with a
        boolean array.

        Allowed inputs are:

        - A single label, e.g. ``5`` or ``'a'``, (note that ``5`` is
          interpreted as a *label* of the index, and **never** as an
          integer position along the index).
        - A list or array of labels, e.g. ``['a', 'b', 'c']``.
        - A slice object with labels, e.g. ``'a':'f'``.

          .. warning:: Note that contrary to usual python slices, **both** the
              start and the stop are included

        - A boolean array of the same length as the axis being sliced,
          e.g. ``[True, False, True]``.
        - A ``callable`` function with one argument (the calling Series or
          DataFrame) and that returns valid output for indexing (one of the above)

        See more at :ref:`Selection by Label <indexing.label>`

        Raises
        ------
        KeyError
            If any items are not found.

        See Also
        --------
        DataFrame.at : Access a single value for a row/column label pair.
        DataFrame.iloc : Access group of rows and columns by integer position(s).
        DataFrame.xs : Returns a cross-section (row(s) or column(s)) from the
            Series/DataFrame.
        Series.loc : Access group of values using labels.

        Examples
        --------
        **Getting values**

        >>> df = pd.DataFrame([[1, 2], [4, 5], [7, 8]],
        ...      index=['cobra', 'viper', 'sidewinder'],
        ...      columns=['max_speed', 'shield'])
        >>> df
                    max_speed  shield
        cobra               1       2
        viper               4       5
        sidewinder          7       8

        Single label. Note this returns the row as a Series.

        >>> df.loc['viper']
        max_speed    4
        shield       5
        Name: viper, dtype: int64

        List of labels. Note using ``[[]]`` returns a DataFrame.

        >>> df.loc[['viper', 'sidewinder']]
                    max_speed  shield
        viper               4       5
        sidewinder          7       8

        Single label for row and column

        >>> df.loc['cobra', 'shield']
        2

        Slice with labels for row and single label for column. As mentioned
        above, note that both the start and stop of the slice are included.

        >>> df.loc['cobra':'viper', 'max_speed']
        cobra    1
        viper    4
        Name: max_speed, dtype: int64

        Boolean list with the same length as the row axis

        >>> df.loc[[False, False, True]]
                    max_speed  shield
        sidewinder          7       8

        Conditional that returns a boolean Series

        >>> df.loc[df['shield'] > 6]
                    max_speed  shield
        sidewinder          7       8

        Conditional that returns a boolean Series with column labels specified

        >>> df.loc[df['shield'] > 6, ['max_speed']]
                    max_speed
        sidewinder          7

        Callable that returns a boolean Series

        >>> df.loc[lambda df: df['shield'] == 8]
                    max_speed  shield
        sidewinder          7       8

        **Setting values**

        Set value for all items matching the list of labels

        >>> df.loc[['viper', 'sidewinder'], ['shield']] = 50
        >>> df
                    max_speed  shield
        cobra               1       2
        viper               4      50
        sidewinder          7      50

        Set value for an entire row

        >>> df.loc['cobra'] = 10
        >>> df
                    max_speed  shield
        cobra              10      10
        viper               4      50
        sidewinder          7      50

        Set value for an entire column

        >>> df.loc[:, 'max_speed'] = 30
        >>> df
                    max_speed  shield
        cobra              30      10
        viper              30      50
        sidewinder         30      50

        Set value for rows matching callable condition

        >>> df.loc[df['shield'] > 35] = 0
        >>> df
                    max_speed  shield
        cobra              30      10
        viper               0       0
        sidewinder          0       0

        **Getting values on a DataFrame with an index that has integer labels**

        Another example using integers for the index

        >>> df = pd.DataFrame([[1, 2], [4, 5], [7, 8]],
        ...      index=[7, 8, 9], columns=['max_speed', 'shield'])
        >>> df
           max_speed  shield
        7          1       2
        8          4       5
        9          7       8

        Slice with integer labels for rows. As mentioned above, note that both
        the start and stop of the slice are included.

        >>> df.loc[7:9]
           max_speed  shield
        7          1       2
        8          4       5
        9          7       8

        **Getting values with a MultiIndex**

        A number of examples using a DataFrame with a MultiIndex

        >>> tuples = [
        ...    ('cobra', 'mark i'), ('cobra', 'mark ii'),
        ...    ('sidewinder', 'mark i'), ('sidewinder', 'mark ii'),
        ...    ('viper', 'mark ii'), ('viper', 'mark iii')
        ... ]
        >>> index = pd.MultiIndex.from_tuples(tuples)
        >>> values = [[12, 2], [0, 4], [10, 20],
        ...         [1, 4], [7, 1], [16, 36]]
        >>> df = pd.DataFrame(values, columns=['max_speed', 'shield'], index=index)
        >>> df
                             max_speed  shield
        cobra      mark i           12       2
                   mark ii           0       4
        sidewinder mark i           10      20
                   mark ii           1       4
        viper      mark ii           7       1
                   mark iii         16      36

        Single label. Note this returns a DataFrame with a single index.

        >>> df.loc['cobra']
                 max_speed  shield
        mark i          12       2
        mark ii          0       4

        Single index tuple. Note this returns a Series.

        >>> df.loc[('cobra', 'mark ii')]
        max_speed    0
        shield       4
        Name: (cobra, mark ii), dtype: int64

        Single label for row and column. Similar to passing in a tuple, this
        returns a Series.

        >>> df.loc['cobra', 'mark i']
        max_speed    12
        shield        2
        Name: (cobra, mark i), dtype: int64

        Single tuple. Note using ``[[]]`` returns a DataFrame.

        >>> df.loc[[('cobra', 'mark ii')]]
                       max_speed  shield
        cobra mark ii          0       4

        Single tuple for the index with a single label for the column

        >>> df.loc[('cobra', 'mark i'), 'shield']
        2

        Slice from index tuple to single label

        >>> df.loc[('cobra', 'mark i'):'viper']
                             max_speed  shield
        cobra      mark i           12       2
                   mark ii           0       4
        sidewinder mark i           10      20
                   mark ii           1       4
        viper      mark ii           7       1
                   mark iii         16      36

        Slice from index tuple to index tuple

        >>> df.loc[('cobra', 'mark i'):('viper', 'mark ii')]
                            max_speed  shield
        cobra      mark i          12       2
                   mark ii          0       4
        sidewinder mark i          10      20
                   mark ii          1       4
        viper      mark ii          7       1
        """
        return _LocIndexer("loc", self)

    @property
    def at(self) -> "_AtIndexer":
        """
        Access a single value for a row/column label pair.

        Similar to ``loc``, in that both provide label-based lookups. Use
        ``at`` if you only need to get or set a single value in a DataFrame
        or Series.

        Raises
        ------
        KeyError
            If 'label' does not exist in DataFrame.

        See Also
        --------
        DataFrame.iat : Access a single value for a row/column pair by integer
            position.
        DataFrame.loc : Access a group of rows and columns by label(s).
        Series.at : Access a single value using a label.

        Examples
        --------
        >>> df = pd.DataFrame([[0, 2, 3], [0, 4, 1], [10, 20, 30]],
        ...                   index=[4, 5, 6], columns=['A', 'B', 'C'])
        >>> df
            A   B   C
        4   0   2   3
        5   0   4   1
        6  10  20  30

        Get value at specified row/column pair

        >>> df.at[4, 'B']
        2

        Set value at specified row/column pair

        >>> df.at[4, 'B'] = 10
        >>> df.at[4, 'B']
        10

        Get value within a Series

        >>> df.loc[5].at['B']
        4
        """
        return _AtIndexer("at", self)

    @property
    def iat(self) -> "_iAtIndexer":
        """
        Access a single value for a row/column pair by integer position.

        Similar to ``iloc``, in that both provide integer-based lookups. Use
        ``iat`` if you only need to get or set a single value in a DataFrame
        or Series.

        Raises
        ------
        IndexError
            When integer position is out of bounds.

        See Also
        --------
        DataFrame.at : Access a single value for a row/column label pair.
        DataFrame.loc : Access a group of rows and columns by label(s).
        DataFrame.iloc : Access a group of rows and columns by integer position(s).

        Examples
        --------
        >>> df = pd.DataFrame([[0, 2, 3], [0, 4, 1], [10, 20, 30]],
        ...                   columns=['A', 'B', 'C'])
        >>> df
            A   B   C
        0   0   2   3
        1   0   4   1
        2  10  20  30

        Get value at specified row/column pair

        >>> df.iat[1, 2]
        1

        Set value at specified row/column pair

        >>> df.iat[1, 2] = 10
        >>> df.iat[1, 2]
        10

        Get value within a series

        >>> df.loc[0].iat[1]
        2
        """
        return _iAtIndexer("iat", self)


class _LocationIndexer(_NDFrameIndexerBase):
    _valid_types: str
    axis = None

    def __call__(self, axis=None):
        # we need to return a copy of ourselves
        new_self = type(self)(self.name, self.obj)

        if axis is not None:
            axis = self.obj._get_axis_number(axis)
        new_self.axis = axis
        return new_self

    def _get_setitem_indexer(self, key):
        """
        Convert a potentially-label-based key into a positional indexer.
        """
        if self.axis is not None:
            return self._convert_tuple(key, is_setter=True)

        ax = self.obj._get_axis(0)

        if isinstance(ax, ABCMultiIndex) and self.name != "iloc":
            try:
                return ax.get_loc(key)
            except (TypeError, KeyError, InvalidIndexError):
                # TypeError e.g. passed a bool
                pass

        if isinstance(key, tuple):
            try:
                return self._convert_tuple(key, is_setter=True)
            except IndexingError:
                pass

        if isinstance(key, range):
            return list(key)

        try:
            return self._convert_to_indexer(key, axis=0, is_setter=True)
        except TypeError as e:

            # invalid indexer type vs 'other' indexing errors
            if "cannot do" in str(e):
                raise
            raise IndexingError(key)

    def __setitem__(self, key, value):
        if isinstance(key, tuple):
            key = tuple(com.apply_if_callable(x, self.obj) for x in key)
        else:
            key = com.apply_if_callable(key, self.obj)
        indexer = self._get_setitem_indexer(key)
        self._has_valid_setitem_indexer(key)

        iloc = self if self.name == "iloc" else self.obj.iloc
        iloc._setitem_with_indexer(indexer, value)

    def _validate_key(self, key, axis: int):
        """
        Ensure that key is valid for current indexer.

        Parameters
        ----------
        key : scalar, slice or list-like
            Key requested.
        axis : int
            Dimension on which the indexing is being made.

        Raises
        ------
        TypeError
            If the key (or some element of it) has wrong type.
        IndexError
            If the key (or some element of it) is out of bounds.
        KeyError
            If the key was not found.
        """
        raise AbstractMethodError(self)

    def _has_valid_tuple(self, key: Tuple):
        """
        Check the key for valid keys across my indexer.
        """
        for i, k in enumerate(key):
            if i >= self.ndim:
                raise IndexingError("Too many indexers")
            try:
                self._validate_key(k, i)
            except ValueError:
                raise ValueError(
                    "Location based indexing can only have "
                    f"[{self._valid_types}] types"
                )

    def _is_nested_tuple_indexer(self, tup: Tuple) -> bool:
        """
        Returns
        -------
        bool
        """
        if any(isinstance(ax, ABCMultiIndex) for ax in self.obj.axes):
            return any(is_nested_tuple(tup, ax) for ax in self.obj.axes)
        return False

    def _convert_tuple(self, key, is_setter: bool = False):
        keyidx = []
        if self.axis is not None:
            axis = self.obj._get_axis_number(self.axis)
            for i in range(self.ndim):
                if i == axis:
                    keyidx.append(
                        self._convert_to_indexer(key, axis=axis, is_setter=is_setter)
                    )
                else:
                    keyidx.append(slice(None))
        else:
            for i, k in enumerate(key):
                if i >= self.ndim:
                    raise IndexingError("Too many indexers")
                idx = self._convert_to_indexer(k, axis=i, is_setter=is_setter)
                keyidx.append(idx)
        return tuple(keyidx)

    def _getitem_tuple_same_dim(self, tup: Tuple):
        """
        Index with indexers that should return an object of the same dimension
        as self.obj.

        This is only called after a failed call to _getitem_lowerdim.
        """
        retval = self.obj
        for i, key in enumerate(tup):
            if com.is_null_slice(key):
                continue

            retval = getattr(retval, self.name)._getitem_axis(key, axis=i)
            # We should never have retval.ndim < self.ndim, as that should
            #  be handled by the _getitem_lowerdim call above.
            assert retval.ndim == self.ndim

        return retval

    def _getitem_lowerdim(self, tup: Tuple):

        # we can directly get the axis result since the axis is specified
        if self.axis is not None:
            axis = self.obj._get_axis_number(self.axis)
            return self._getitem_axis(tup, axis=axis)

        # we may have a nested tuples indexer here
        if self._is_nested_tuple_indexer(tup):
            return self._getitem_nested_tuple(tup)

        # we maybe be using a tuple to represent multiple dimensions here
        ax0 = self.obj._get_axis(0)
        # ...but iloc should handle the tuple as simple integer-location
        # instead of checking it as multiindex representation (GH 13797)
        if isinstance(ax0, ABCMultiIndex) and self.name != "iloc":
            result = self._handle_lowerdim_multi_index_axis0(tup)
            if result is not None:
                return result

        if len(tup) > self.ndim:
            raise IndexingError("Too many indexers. handle elsewhere")

        for i, key in enumerate(tup):
            if is_label_like(key) or isinstance(key, tuple):
                section = self._getitem_axis(key, axis=i)

                # we have yielded a scalar ?
                if not is_list_like_indexer(section):
                    return section

                elif section.ndim == self.ndim:
                    # we're in the middle of slicing through a MultiIndex
                    # revise the key wrt to `section` by inserting an _NS
                    new_key = tup[:i] + (_NS,) + tup[i + 1 :]

                else:
                    # Note: the section.ndim == self.ndim check above
                    #  rules out having DataFrame here, so we dont need to worry
                    #  about transposing.
                    new_key = tup[:i] + tup[i + 1 :]

                    if len(new_key) == 1:
                        new_key = new_key[0]

                # Slices should return views, but calling iloc/loc with a null
                # slice returns a new object.
                if com.is_null_slice(new_key):
                    return section
                # This is an elided recursive call to iloc/loc/etc'
                return getattr(section, self.name)[new_key]

        raise IndexingError("not applicable")

    def _getitem_nested_tuple(self, tup: Tuple):
        # we have a nested tuple so have at least 1 multi-index level
        # we should be able to match up the dimensionality here

        # we have too many indexers for our dim, but have at least 1
        # multi-index dimension, try to see if we have something like
        # a tuple passed to a series with a multi-index
        if len(tup) > self.ndim:
            if self.name != "loc":
                # This should never be reached, but lets be explicit about it
                raise ValueError("Too many indices")
            result = self._handle_lowerdim_multi_index_axis0(tup)
            if result is not None:
                return result

            # this is a series with a multi-index specified a tuple of
            # selectors
            axis = self.axis or 0
            return self._getitem_axis(tup, axis=axis)

        # handle the multi-axis by taking sections and reducing
        # this is iterative
        obj = self.obj
        axis = 0
        for i, key in enumerate(tup):

            if com.is_null_slice(key):
                axis += 1
                continue

            current_ndim = obj.ndim
            obj = getattr(obj, self.name)._getitem_axis(key, axis=axis)
            axis += 1

            # if we have a scalar, we are done
            if is_scalar(obj) or not hasattr(obj, "ndim"):
                break

            # has the dim of the obj changed?
            # GH 7199
            if obj.ndim < current_ndim:
                axis -= 1

        return obj

    def _convert_to_indexer(self, key, axis: int, is_setter: bool = False):
        raise AbstractMethodError(self)

    def __getitem__(self, key):
        if type(key) is tuple:
            key = tuple(com.apply_if_callable(x, self.obj) for x in key)
            if self._is_scalar_access(key):
                try:
                    return self.obj._get_value(*key, takeable=self._takeable)
                except (KeyError, IndexError, AttributeError):
                    # AttributeError for IntervalTree get_value
                    pass
            return self._getitem_tuple(key)
        else:
            # we by definition only have the 0th axis
            axis = self.axis or 0

            maybe_callable = com.apply_if_callable(key, self.obj)
            return self._getitem_axis(maybe_callable, axis=axis)

    def _is_scalar_access(self, key: Tuple):
        raise NotImplementedError()

    def _getitem_tuple(self, tup: Tuple):
        raise AbstractMethodError(self)

    def _getitem_axis(self, key, axis: int):
        raise NotImplementedError()

    def _has_valid_setitem_indexer(self, indexer) -> bool:
        raise AbstractMethodError(self)

    def _getbool_axis(self, key, axis: int):
        # caller is responsible for ensuring non-None axis
        labels = self.obj._get_axis(axis)
        key = check_bool_indexer(labels, key)
        inds = key.nonzero()[0]
        return self.obj._take_with_is_copy(inds, axis=axis)


@Appender(IndexingMixin.loc.__doc__)
class _LocIndexer(_LocationIndexer):
    _takeable: bool = False
    _valid_types = (
        "labels (MUST BE IN THE INDEX), slices of labels (BOTH "
        "endpoints included! Can be slices of integers if the "
        "index is integers), listlike of labels, boolean"
    )

    # -------------------------------------------------------------------
    # Key Checks

    @Appender(_LocationIndexer._validate_key.__doc__)
    def _validate_key(self, key, axis: int):

        # valid for a collection of labels (we check their presence later)
        # slice of labels (where start-end in labels)
        # slice of integers (only if in the labels)
        # boolean

        if isinstance(key, slice):
            return

        if com.is_bool_indexer(key):
            return

        if not is_list_like_indexer(key):
            labels = self.obj._get_axis(axis)
            labels._convert_scalar_indexer(key, kind="loc")

    def _has_valid_setitem_indexer(self, indexer) -> bool:
        return True

    def _is_scalar_access(self, key: Tuple) -> bool:
        """
        Returns
        -------
        bool
        """
        # this is a shortcut accessor to both .loc and .iloc
        # that provide the equivalent access of .at and .iat
        # a) avoid getting things via sections and (to minimize dtype changes)
        # b) provide a performant path
        if len(key) != self.ndim:
            return False

        for i, k in enumerate(key):
            if not is_scalar(k):
                return False

            ax = self.obj.axes[i]
            if isinstance(ax, ABCMultiIndex):
                return False

            if isinstance(k, str) and ax._supports_partial_string_indexing:
                # partial string indexing, df.loc['2000', 'A']
                # should not be considered scalar
                return False

            if not ax.is_unique:
                return False

        return True

    # -------------------------------------------------------------------
    # MultiIndex Handling

    def _multi_take_opportunity(self, tup: Tuple) -> bool:
        """
        Check whether there is the possibility to use ``_multi_take``.

        Currently the limit is that all axes being indexed, must be indexed with
        list-likes.

        Parameters
        ----------
        tup : tuple
            Tuple of indexers, one per axis.

        Returns
        -------
        bool
            Whether the current indexing,
            can be passed through `_multi_take`.
        """
        if not all(is_list_like_indexer(x) for x in tup):
            return False

        # just too complicated
        if any(com.is_bool_indexer(x) for x in tup):
            return False

        return True

    def _multi_take(self, tup: Tuple):
        """
        Create the indexers for the passed tuple of keys, and
        executes the take operation. This allows the take operation to be
        executed all at once, rather than once for each dimension.
        Improving efficiency.

        Parameters
        ----------
        tup : tuple
            Tuple of indexers, one per axis.

        Returns
        -------
        values: same type as the object being indexed
        """
        # GH 836
        d = {
            axis: self._get_listlike_indexer(key, axis)
            for (key, axis) in zip(tup, self.obj._AXIS_ORDERS)
        }
        return self.obj._reindex_with_indexers(d, copy=True, allow_dups=True)

    # -------------------------------------------------------------------

    def _get_partial_string_timestamp_match_key(self, key, labels):
        """
        Translate any partial string timestamp matches in key, returning the
        new key.

        (GH 10331)
        """
        if isinstance(labels, ABCMultiIndex):
            if (
                isinstance(key, str)
                and labels.levels[0]._supports_partial_string_indexing
            ):
                # Convert key '2016-01-01' to
                # ('2016-01-01'[, slice(None, None, None)]+)
                key = tuple([key] + [slice(None)] * (len(labels.levels) - 1))

            if isinstance(key, tuple):
                # Convert (..., '2016-01-01', ...) in tuple to
                # (..., slice('2016-01-01', '2016-01-01', None), ...)
                new_key = []
                for i, component in enumerate(key):
                    if (
                        isinstance(component, str)
                        and labels.levels[i]._supports_partial_string_indexing
                    ):
                        new_key.append(slice(component, component, None))
                    else:
                        new_key.append(component)
                key = tuple(new_key)

        return key

    def _getitem_iterable(self, key, axis: int):
        """
        Index current object with an an iterable collection of keys.

        Parameters
        ----------
        key : iterable
            Targeted labels.
        axis: int
            Dimension on which the indexing is being made.

        Raises
        ------
        KeyError
            If no key was found. Will change in the future to raise if not all
            keys were found.

        Returns
        -------
        scalar, DataFrame, or Series: indexed value(s).
        """
        # we assume that not com.is_bool_indexer(key), as that is
        #  handled before we get here.
        self._validate_key(key, axis)

        # A collection of keys
        keyarr, indexer = self._get_listlike_indexer(key, axis, raise_missing=False)
        return self.obj._reindex_with_indexers(
            {axis: [keyarr, indexer]}, copy=True, allow_dups=True
        )

    def _getitem_tuple(self, tup: Tuple):
        try:
            return self._getitem_lowerdim(tup)
        except IndexingError:
            pass

        # no multi-index, so validate all of the indexers
        self._has_valid_tuple(tup)

        # ugly hack for GH #836
        if self._multi_take_opportunity(tup):
            return self._multi_take(tup)

        return self._getitem_tuple_same_dim(tup)

    def _get_label(self, label, axis: int):
        if self.ndim == 1:
            # for perf reasons we want to try _xs first
            # as its basically direct indexing
            # but will fail when the index is not present
            # see GH5667
            return self.obj._xs(label, axis=axis)
        elif isinstance(label, tuple) and isinstance(label[axis], slice):
            raise IndexingError("no slices here, handle elsewhere")

        return self.obj._xs(label, axis=axis)

    def _handle_lowerdim_multi_index_axis0(self, tup: Tuple):
        # we have an axis0 multi-index, handle or raise
        axis = self.axis or 0
        try:
            # fast path for series or for tup devoid of slices
            return self._get_label(tup, axis=axis)
        except TypeError:
            # slices are unhashable
            pass
        except KeyError as ek:
            # raise KeyError if number of indexers match
            # else IndexingError will be raised
            if len(tup) <= self.obj.index.nlevels and len(tup) > self.ndim:
                raise ek

        return None

    def _getitem_axis(self, key, axis: int):
        key = item_from_zerodim(key)
        if is_iterator(key):
            key = list(key)

        labels = self.obj._get_axis(axis)
        key = self._get_partial_string_timestamp_match_key(key, labels)

        if isinstance(key, slice):
            self._validate_key(key, axis)
            return self._get_slice_axis(key, axis=axis)
        elif com.is_bool_indexer(key):
            return self._getbool_axis(key, axis=axis)
        elif is_list_like_indexer(key):

            # convert various list-like indexers
            # to a list of keys
            # we will use the *values* of the object
            # and NOT the index if its a PandasObject
            if isinstance(labels, ABCMultiIndex):

                if isinstance(key, (ABCSeries, np.ndarray)) and key.ndim <= 1:
                    # Series, or 0,1 ndim ndarray
                    # GH 14730
                    key = list(key)
                elif isinstance(key, ABCDataFrame):
                    # GH 15438
                    raise NotImplementedError(
                        "Indexing a MultiIndex with a "
                        "DataFrame key is not "
                        "implemented"
                    )
                elif hasattr(key, "ndim") and key.ndim > 1:
                    raise NotImplementedError(
                        "Indexing a MultiIndex with a "
                        "multidimensional key is not "
                        "implemented"
                    )

                if (
                    not isinstance(key, tuple)
                    and len(key)
                    and not isinstance(key[0], tuple)
                ):
                    key = tuple([key])

            # an iterable multi-selection
            if not (isinstance(key, tuple) and isinstance(labels, ABCMultiIndex)):

                if hasattr(key, "ndim") and key.ndim > 1:
                    raise ValueError("Cannot index with multidimensional key")

                return self._getitem_iterable(key, axis=axis)

            # nested tuple slicing
            if is_nested_tuple(key, labels):
                locs = labels.get_locs(key)
                indexer = [slice(None)] * self.ndim
                indexer[axis] = locs
                return self.obj.iloc[tuple(indexer)]

        # fall thru to straight lookup
        self._validate_key(key, axis)
        return self._get_label(key, axis=axis)

    def _get_slice_axis(self, slice_obj: slice, axis: int):
        """
        This is pretty simple as we just have to deal with labels.
        """
        # caller is responsible for ensuring non-None axis
        obj = self.obj
        if not need_slice(slice_obj):
            return obj.copy(deep=False)

        labels = obj._get_axis(axis)
        indexer = labels.slice_indexer(
            slice_obj.start, slice_obj.stop, slice_obj.step, kind="loc"
        )

        if isinstance(indexer, slice):
            return self.obj._slice(indexer, axis=axis)
        else:
            # DatetimeIndex overrides Index.slice_indexer and may
            #  return a DatetimeIndex instead of a slice object.
            return self.obj.take(indexer, axis=axis)

    def _convert_to_indexer(self, key, axis: int, is_setter: bool = False):
        """
        Convert indexing key into something we can use to do actual fancy
        indexing on a ndarray.

        Examples
        ix[:5] -> slice(0, 5)
        ix[[1,2,3]] -> [1,2,3]
        ix[['foo', 'bar', 'baz']] -> [i, j, k] (indices of foo, bar, baz)

        Going by Zen of Python?
        'In the face of ambiguity, refuse the temptation to guess.'
        raise AmbiguousIndexError with integer labels?
        - No, prefer label-based indexing
        """
        labels = self.obj._get_axis(axis)

        if isinstance(key, slice):
            return labels._convert_slice_indexer(key, kind="loc")

        if is_scalar(key):
            # try to find out correct indexer, if not type correct raise
            try:
                key = labels._convert_scalar_indexer(key, kind="loc")
            except TypeError:
                # but we will allow setting
                if not is_setter:
                    raise

        # see if we are positional in nature
        is_int_index = labels.is_integer()
        is_int_positional = is_integer(key) and not is_int_index

        if is_scalar(key) or isinstance(labels, ABCMultiIndex):
            # Otherwise get_loc will raise InvalidIndexError

            # if we are a label return me
            try:
                return labels.get_loc(key)
            except LookupError:
                if isinstance(key, tuple) and isinstance(labels, ABCMultiIndex):
                    if len(key) == labels.nlevels:
                        return {"key": key}
                    raise
            except TypeError:
                pass
            except ValueError:
                if not is_int_positional:
                    raise

        # a positional
        if is_int_positional:

            # if we are setting and its not a valid location
            # its an insert which fails by definition

            # always valid
            return {"key": key}

        if is_nested_tuple(key, labels):
            return labels.get_locs(key)

        elif is_list_like_indexer(key):

            if com.is_bool_indexer(key):
                key = check_bool_indexer(labels, key)
                (inds,) = key.nonzero()
                return inds
            else:
                # When setting, missing keys are not allowed, even with .loc:
                return self._get_listlike_indexer(key, axis, raise_missing=True)[1]
        else:
            try:
                return labels.get_loc(key)
            except LookupError:
                # allow a not found key only if we are a setter
                if not is_list_like_indexer(key):
                    return {"key": key}
                raise

    def _get_listlike_indexer(self, key, axis: int, raise_missing: bool = False):
        """
        Transform a list-like of keys into a new index and an indexer.

        Parameters
        ----------
        key : list-like
            Targeted labels.
        axis: int
            Dimension on which the indexing is being made.
        raise_missing: bool, default False
            Whether to raise a KeyError if some labels were not found.
            Will be removed in the future, and then this method will always behave as
            if ``raise_missing=True``.

        Raises
        ------
        KeyError
            If at least one key was requested but none was found, and
            raise_missing=True.

        Returns
        -------
        keyarr: Index
            New index (coinciding with 'key' if the axis is unique).
        values : array-like
            Indexer for the return object, -1 denotes keys not found.
        """
        ax = self.obj._get_axis(axis)

        # Have the index compute an indexer or return None
        # if it cannot handle:
        indexer, keyarr = ax._convert_listlike_indexer(key)
        # We only act on all found values:
        if indexer is not None and (indexer != -1).all():
            self._validate_read_indexer(key, indexer, axis, raise_missing=raise_missing)
            return ax[indexer], indexer

        if ax.is_unique and not getattr(ax, "is_overlapping", False):
            indexer = ax.get_indexer_for(key)
            keyarr = ax.reindex(keyarr)[0]
        else:
            keyarr, indexer, new_indexer = ax._reindex_non_unique(keyarr)

        self._validate_read_indexer(keyarr, indexer, axis, raise_missing=raise_missing)
        return keyarr, indexer

    def _validate_read_indexer(
        self, key, indexer, axis: int, raise_missing: bool = False
    ):
        """
        Check that indexer can be used to return a result.

        e.g. at least one element was found,
        unless the list of keys was actually empty.

        Parameters
        ----------
        key : list-like
            Targeted labels (only used to show correct error message).
        indexer: array-like of booleans
            Indices corresponding to the key,
            (with -1 indicating not found).
        axis: int
            Dimension on which the indexing is being made.
        raise_missing: bool
            Whether to raise a KeyError if some labels are not found. Will be
            removed in the future, and then this method will always behave as
            if raise_missing=True.

        Raises
        ------
        KeyError
            If at least one key was requested but none was found, and
            raise_missing=True.
        """
        ax = self.obj._get_axis(axis)

        if len(key) == 0:
            return

        # Count missing values:
        missing = (indexer < 0).sum()

        if missing:
            if missing == len(indexer):
                axis_name = self.obj._get_axis_name(axis)
                raise KeyError(f"None of [{key}] are in the [{axis_name}]")

            # We (temporarily) allow for some missing keys with .loc, except in
            # some cases (e.g. setting) in which "raise_missing" will be False
            if not (self.name == "loc" and not raise_missing):
                not_found = list(set(key) - set(ax))
                raise KeyError(f"{not_found} not in index")

            # we skip the warning on Categorical/Interval
            # as this check is actually done (check for
            # non-missing values), but a bit later in the
            # code, so we want to avoid warning & then
            # just raising
            if not (ax.is_categorical() or ax.is_interval()):
                raise KeyError(
                    "Passing list-likes to .loc or [] with any missing labels "
                    "is no longer supported, see "
                    "https://pandas.pydata.org/pandas-docs/stable/user_guide/indexing.html#deprecate-loc-reindex-listlike"  # noqa:E501
                )


@Appender(IndexingMixin.iloc.__doc__)
class _iLocIndexer(_LocationIndexer):
    _valid_types = (
        "integer, integer slice (START point is INCLUDED, END "
        "point is EXCLUDED), listlike of integers, boolean array"
    )
    _takeable = True

    # -------------------------------------------------------------------
    # Key Checks

    def _validate_key(self, key, axis: int):
        if com.is_bool_indexer(key):
            if hasattr(key, "index") and isinstance(key.index, Index):
                if key.index.inferred_type == "integer":
                    raise NotImplementedError(
                        "iLocation based boolean "
                        "indexing on an integer type "
                        "is not available"
                    )
                raise ValueError(
                    "iLocation based boolean indexing cannot use "
                    "an indexable as a mask"
                )
            return

        if isinstance(key, slice):
            return
        elif is_integer(key):
            self._validate_integer(key, axis)
        elif isinstance(key, tuple):
            # a tuple should already have been caught by this point
            # so don't treat a tuple as a valid indexer
            raise IndexingError("Too many indexers")
        elif is_list_like_indexer(key):
            arr = np.array(key)
            len_axis = len(self.obj._get_axis(axis))

            # check that the key has a numeric dtype
            if not is_numeric_dtype(arr.dtype):
                raise IndexError(f".iloc requires numeric indexers, got {arr}")

            # check that the key does not exceed the maximum size of the index
            if len(arr) and (arr.max() >= len_axis or arr.min() < -len_axis):
                raise IndexError("positional indexers are out-of-bounds")
        else:
            raise ValueError(f"Can only index by location with a [{self._valid_types}]")

    def _has_valid_setitem_indexer(self, indexer):
        self._has_valid_positional_setitem_indexer(indexer)

    def _has_valid_positional_setitem_indexer(self, indexer) -> bool:
        """
        Validate that a positional indexer cannot enlarge its target
        will raise if needed, does not modify the indexer externally.

        Returns
        -------
        bool
        """
        if isinstance(indexer, dict):
            raise IndexError(f"{self.name} cannot enlarge its target object")
        else:
            if not isinstance(indexer, tuple):
                indexer = _tuplify(self.ndim, indexer)
            for ax, i in zip(self.obj.axes, indexer):
                if isinstance(i, slice):
                    # should check the stop slice?
                    pass
                elif is_list_like_indexer(i):
                    # should check the elements?
                    pass
                elif is_integer(i):
                    if i >= len(ax):
                        raise IndexError(
                            f"{self.name} cannot enlarge its target object"
                        )
                elif isinstance(i, dict):
                    raise IndexError(f"{self.name} cannot enlarge its target object")

        return True

    def _is_scalar_access(self, key: Tuple) -> bool:
        """
        Returns
        -------
        bool
        """
        # this is a shortcut accessor to both .loc and .iloc
        # that provide the equivalent access of .at and .iat
        # a) avoid getting things via sections and (to minimize dtype changes)
        # b) provide a performant path
        if len(key) != self.ndim:
            return False

        for i, k in enumerate(key):
            if not is_integer(k):
                return False

            ax = self.obj.axes[i]
            if not ax.is_unique:
                return False

        return True

    def _validate_integer(self, key: int, axis: int) -> None:
        """
        Check that 'key' is a valid position in the desired axis.

        Parameters
        ----------
        key : int
            Requested position.
        axis : int
            Desired axis.

        Raises
        ------
        IndexError
            If 'key' is not a valid position in axis 'axis'.
        """
        len_axis = len(self.obj._get_axis(axis))
        if key >= len_axis or key < -len_axis:
            raise IndexError("single positional indexer is out-of-bounds")

    # -------------------------------------------------------------------

    def _getitem_tuple(self, tup: Tuple):

        self._has_valid_tuple(tup)
        try:
            return self._getitem_lowerdim(tup)
        except IndexingError:
            pass

        return self._getitem_tuple_same_dim(tup)

    def _get_list_axis(self, key, axis: int):
        """
        Return Series values by list or array of integers.

        Parameters
        ----------
        key : list-like positional indexer
        axis : int

        Returns
        -------
        Series object

        Notes
        -----
        `axis` can only be zero.
        """
        try:
            return self.obj._take_with_is_copy(key, axis=axis)
        except IndexError:
            # re-raise with different error message
            raise IndexError("positional indexers are out-of-bounds")

    def _getitem_axis(self, key, axis: int):
        if isinstance(key, slice):
            return self._get_slice_axis(key, axis=axis)

        if isinstance(key, list):
            key = np.asarray(key)

        if com.is_bool_indexer(key):
            self._validate_key(key, axis)
            return self._getbool_axis(key, axis=axis)

        # a list of integers
        elif is_list_like_indexer(key):
            return self._get_list_axis(key, axis=axis)

        # a single integer
        else:
            key = item_from_zerodim(key)
            if not is_integer(key):
                raise TypeError("Cannot index by location index with a non-integer key")

            # validate the location
            self._validate_integer(key, axis)

            return self.obj._ixs(key, axis=axis)

    def _get_slice_axis(self, slice_obj: slice, axis: int):
        # caller is responsible for ensuring non-None axis
        obj = self.obj

        if not need_slice(slice_obj):
            return obj.copy(deep=False)

        labels = obj._get_axis(axis)
        labels._validate_positional_slice(slice_obj)
        return self.obj._slice(slice_obj, axis=axis)

    def _convert_to_indexer(self, key, axis: int, is_setter: bool = False):
        """
        Much simpler as we only have to deal with our valid types.
        """
        labels = self.obj._get_axis(axis)

        # make need to convert a float key
        if isinstance(key, slice):
            labels._validate_positional_slice(key)
            return key

        elif is_float(key):
            labels._validate_indexer("positional", key, "iloc")
            return key

        self._validate_key(key, axis)
        return key

    # -------------------------------------------------------------------

    def _setitem_with_indexer(self, indexer, value):
        """
        _setitem_with_indexer is for setting values on a Series/DataFrame
        using positional indexers.

        If the relevant keys are not present, the Series/DataFrame may be
        expanded.

        This method is currently broken when dealing with non-unique Indexes,
        since it goes from positional indexers back to labels when calling
        BlockManager methods, see GH#12991, GH#22046, GH#15686.
        """

        # also has the side effect of consolidating in-place
        from pandas import Series

        info_axis = self.obj._info_axis_number

        # maybe partial set
        take_split_path = self.obj._is_mixed_type

        # if there is only one block/type, still have to take split path
        # unless the block is one-dimensional or it can hold the value
        if not take_split_path and self.obj._data.blocks:
            (blk,) = self.obj._data.blocks
            if 1 < blk.ndim:  # in case of dict, keys are indices
                val = list(value.values()) if isinstance(value, dict) else value
                take_split_path = not blk._can_hold_element(val)

        # if we have any multi-indexes that have non-trivial slices
        # (not null slices) then we must take the split path, xref
        # GH 10360, GH 27841
        if isinstance(indexer, tuple) and len(indexer) == len(self.obj.axes):
            for i, ax in zip(indexer, self.obj.axes):
                if isinstance(ax, ABCMultiIndex) and not (
                    is_integer(i) or com.is_null_slice(i)
                ):
                    take_split_path = True
                    break

        if isinstance(indexer, tuple):
            nindexer = []
            for i, idx in enumerate(indexer):
                if isinstance(idx, dict):

                    # reindex the axis to the new value
                    # and set inplace
                    key, _ = convert_missing_indexer(idx)

                    # if this is the items axes, then take the main missing
                    # path first
                    # this correctly sets the dtype and avoids cache issues
                    # essentially this separates out the block that is needed
                    # to possibly be modified
                    if self.ndim > 1 and i == self.obj._info_axis_number:

                        # add the new item, and set the value
                        # must have all defined axes if we have a scalar
                        # or a list-like on the non-info axes if we have a
                        # list-like
                        len_non_info_axes = (
                            len(_ax) for _i, _ax in enumerate(self.obj.axes) if _i != i
                        )
                        if any(not l for l in len_non_info_axes):
                            if not is_list_like_indexer(value):
                                raise ValueError(
                                    "cannot set a frame with no "
                                    "defined index and a scalar"
                                )
                            self.obj[key] = value
                            return

                        # add a new item with the dtype setup
                        self.obj[key] = _infer_fill_value(value)

                        new_indexer = convert_from_missing_indexer_tuple(
                            indexer, self.obj.axes
                        )
                        self._setitem_with_indexer(new_indexer, value)

                        return

                    # reindex the axis
                    # make sure to clear the cache because we are
                    # just replacing the block manager here
                    # so the object is the same
                    index = self.obj._get_axis(i)
                    labels = index.insert(len(index), key)
                    self.obj._data = self.obj.reindex(labels, axis=i)._data
                    self.obj._maybe_update_cacher(clear=True)
                    self.obj._is_copy = None

                    nindexer.append(labels.get_loc(key))

                else:
                    nindexer.append(idx)

            indexer = tuple(nindexer)
        else:

            indexer, missing = convert_missing_indexer(indexer)

            if missing:
                self._setitem_with_indexer_missing(indexer, value)
                return

        # set
        item_labels = self.obj._get_axis(info_axis)

        # align and set the values
        if take_split_path:
            # Above we only set take_split_path to True for 2D cases
            assert self.ndim == 2
            assert info_axis == 1

            if not isinstance(indexer, tuple):
                indexer = _tuplify(self.ndim, indexer)

            if isinstance(value, ABCSeries):
                value = self._align_series(indexer, value)

            info_idx = indexer[info_axis]
            if is_integer(info_idx):
                info_idx = [info_idx]
            labels = item_labels[info_idx]

            if len(labels) == 1:
                # We can operate on a single column
                item = labels[0]
                idx = indexer[0]

                plane_indexer = tuple([idx])
                lplane_indexer = length_of_indexer(plane_indexer[0], self.obj.index)
                # lplane_indexer gives the expected length of obj[idx]

                # require that we are setting the right number of values that
                # we are indexing
                if is_list_like_indexer(value) and 0 != lplane_indexer != len(value):
                    # Exclude zero-len for e.g. boolean masking that is all-false
                    raise ValueError(
                        "cannot set using a multi-index "
                        "selection indexer with a different "
                        "length than the value"
                    )

            # non-mi
            else:
                plane_indexer = indexer[:1]
                lplane_indexer = length_of_indexer(plane_indexer[0], self.obj.index)

            def setter(item, v):
                ser = self.obj[item]
                pi = plane_indexer[0] if lplane_indexer == 1 else plane_indexer

                # perform the equivalent of a setitem on the info axis
                # as we have a null slice or a slice with full bounds
                # which means essentially reassign to the columns of a
                # multi-dim object
                # GH6149 (null slice), GH10408 (full bounds)
                if isinstance(pi, tuple) and all(
                    com.is_null_slice(idx) or com.is_full_slice(idx, len(self.obj))
                    for idx in pi
                ):
                    ser = v
                else:
                    # set the item, possibly having a dtype change
                    ser._consolidate_inplace()
                    ser = ser.copy()
                    ser._data = ser._data.setitem(indexer=pi, value=v)
                    ser._maybe_update_cacher(clear=True)

                # reset the sliced object if unique
                self.obj[item] = ser

            # we need an iterable, with a ndim of at least 1
            # eg. don't pass through np.array(0)
            if is_list_like_indexer(value) and getattr(value, "ndim", 1) > 0:

                # we have an equal len Frame
                if isinstance(value, ABCDataFrame):
                    sub_indexer = list(indexer)
                    multiindex_indexer = isinstance(labels, ABCMultiIndex)

                    for item in labels:
                        if item in value:
                            sub_indexer[info_axis] = item
                            v = self._align_series(
                                tuple(sub_indexer), value[item], multiindex_indexer
                            )
                        else:
                            v = np.nan

                        setter(item, v)

                # we have an equal len ndarray/convertible to our labels
                # hasattr first, to avoid coercing to ndarray without reason.
                # But we may be relying on the ndarray coercion to check ndim.
                # Why not just convert to an ndarray earlier on if needed?
                elif np.ndim(value) == 2:

                    # note that this coerces the dtype if we are mixed
                    # GH 7551
                    value = np.array(value, dtype=object)
                    if len(labels) != value.shape[1]:
                        raise ValueError(
                            "Must have equal len keys and value "
                            "when setting with an ndarray"
                        )

                    for i, item in enumerate(labels):

                        # setting with a list, recoerces
                        setter(item, value[:, i].tolist())

                # we have an equal len list/ndarray
                elif _can_do_equal_len(
                    labels, value, plane_indexer, lplane_indexer, self.obj
                ):
                    setter(labels[0], value)

                # per label values
                else:

                    if len(labels) != len(value):
                        raise ValueError(
                            "Must have equal len keys and value "
                            "when setting with an iterable"
                        )

                    for item, v in zip(labels, value):
                        setter(item, v)
            else:

                # scalar
                for item in labels:
                    setter(item, value)

        else:
            if isinstance(indexer, tuple):
                indexer = maybe_convert_ix(*indexer)

                # if we are setting on the info axis ONLY
                # set using those methods to avoid block-splitting
                # logic here
                if (
                    len(indexer) > info_axis
                    and is_integer(indexer[info_axis])
                    and all(
                        com.is_null_slice(idx)
                        for i, idx in enumerate(indexer)
                        if i != info_axis
                    )
                    and item_labels.is_unique
                ):
                    self.obj[item_labels[indexer[info_axis]]] = value
                    return

            if isinstance(value, (ABCSeries, dict)):
                # TODO(EA): ExtensionBlock.setitem this causes issues with
                # setting for extensionarrays that store dicts. Need to decide
                # if it's worth supporting that.
                value = self._align_series(indexer, Series(value))

            elif isinstance(value, ABCDataFrame):
                value = self._align_frame(indexer, value)

            # check for chained assignment
            self.obj._check_is_chained_assignment_possible()

            # actually do the set
            self.obj._consolidate_inplace()
            self.obj._data = self.obj._data.setitem(indexer=indexer, value=value)
            self.obj._maybe_update_cacher(clear=True)

    def _setitem_with_indexer_missing(self, indexer, value):
        """
        Insert new row(s) or column(s) into the Series or DataFrame.
        """
        from pandas import Series

        # reindex the axis to the new value
        # and set inplace
        if self.ndim == 1:
            index = self.obj.index
            new_index = index.insert(len(index), indexer)

            # we have a coerced indexer, e.g. a float
            # that matches in an Int64Index, so
            # we will not create a duplicate index, rather
            # index to that element
            # e.g. 0.0 -> 0
            # GH#12246
            if index.is_unique:
                new_indexer = index.get_indexer([new_index[-1]])
                if (new_indexer != -1).any():
                    return self._setitem_with_indexer(new_indexer, value)

            # this preserves dtype of the value
            new_values = Series([value])._values
            if len(self.obj._values):
                # GH#22717 handle casting compatibility that np.concatenate
                #  does incorrectly
                new_values = concat_compat([self.obj._values, new_values])
            self.obj._data = self.obj._constructor(
                new_values, index=new_index, name=self.obj.name
            )._data
            self.obj._maybe_update_cacher(clear=True)

        elif self.ndim == 2:

            if not len(self.obj.columns):
                # no columns and scalar
                raise ValueError("cannot set a frame with no defined columns")

            if isinstance(value, ABCSeries):
                # append a Series
                value = value.reindex(index=self.obj.columns, copy=True)
                value.name = indexer

            else:
                # a list-list
                if is_list_like_indexer(value):
                    # must have conforming columns
                    if len(value) != len(self.obj.columns):
                        raise ValueError("cannot set a row with mismatched columns")

                value = Series(value, index=self.obj.columns, name=indexer)

            self.obj._data = self.obj.append(value)._data
            self.obj._maybe_update_cacher(clear=True)

    def _align_series(self, indexer, ser: ABCSeries, multiindex_indexer: bool = False):
        """
        Parameters
        ----------
        indexer : tuple, slice, scalar
            Indexer used to get the locations that will be set to `ser`.
        ser : pd.Series
            Values to assign to the locations specified by `indexer`.
        multiindex_indexer : boolean, optional
            Defaults to False. Should be set to True if `indexer` was from
            a `pd.MultiIndex`, to avoid unnecessary broadcasting.

        Returns
        -------
        `np.array` of `ser` broadcast to the appropriate shape for assignment
        to the locations selected by `indexer`
        """
        if isinstance(indexer, (slice, np.ndarray, list, Index)):
            indexer = tuple([indexer])

        if isinstance(indexer, tuple):

            # flatten np.ndarray indexers
            def ravel(i):
                return i.ravel() if isinstance(i, np.ndarray) else i

            indexer = tuple(map(ravel, indexer))

            aligners = [not com.is_null_slice(idx) for idx in indexer]
            sum_aligners = sum(aligners)
            single_aligner = sum_aligners == 1
            is_frame = self.ndim == 2
            obj = self.obj

            # are we a single alignable value on a non-primary
            # dim (e.g. panel: 1,2, or frame: 0) ?
            # hence need to align to a single axis dimension
            # rather that find all valid dims

            # frame
            if is_frame:
                single_aligner = single_aligner and aligners[0]

            # we have a frame, with multiple indexers on both axes; and a
            # series, so need to broadcast (see GH5206)
            if sum_aligners == self.ndim and all(is_sequence(_) for _ in indexer):
                ser = ser.reindex(obj.axes[0][indexer[0]], copy=True)._values

                # single indexer
                if len(indexer) > 1 and not multiindex_indexer:
                    len_indexer = len(indexer[1])
                    ser = np.tile(ser, len_indexer).reshape(len_indexer, -1).T

                return ser

            for i, idx in enumerate(indexer):
                ax = obj.axes[i]

                # multiple aligners (or null slices)
                if is_sequence(idx) or isinstance(idx, slice):
                    if single_aligner and com.is_null_slice(idx):
                        continue
                    new_ix = ax[idx]
                    if not is_list_like_indexer(new_ix):
                        new_ix = Index([new_ix])
                    else:
                        new_ix = Index(new_ix)
                    if ser.index.equals(new_ix) or not len(new_ix):
                        return ser._values.copy()

                    return ser.reindex(new_ix)._values

                # 2 dims
                elif single_aligner:

                    # reindex along index
                    ax = self.obj.axes[1]
                    if ser.index.equals(ax) or not len(ax):
                        return ser._values.copy()
                    return ser.reindex(ax)._values

        elif is_scalar(indexer):
            ax = self.obj._get_axis(1)

            if ser.index.equals(ax):
                return ser._values.copy()

            return ser.reindex(ax)._values

        raise ValueError("Incompatible indexer with Series")

    def _align_frame(self, indexer, df: ABCDataFrame):
        is_frame = self.ndim == 2

        if isinstance(indexer, tuple):

            idx, cols = None, None
            sindexers = []
            for i, ix in enumerate(indexer):
                ax = self.obj.axes[i]
                if is_sequence(ix) or isinstance(ix, slice):
                    if isinstance(ix, np.ndarray):
                        ix = ix.ravel()
                    if idx is None:
                        idx = ax[ix]
                    elif cols is None:
                        cols = ax[ix]
                    else:
                        break
                else:
                    sindexers.append(i)

            if idx is not None and cols is not None:

                if df.index.equals(idx) and df.columns.equals(cols):
                    val = df.copy()._values
                else:
                    val = df.reindex(idx, columns=cols)._values
                return val

        elif (isinstance(indexer, slice) or is_list_like_indexer(indexer)) and is_frame:
            ax = self.obj.index[indexer]
            if df.index.equals(ax):
                val = df.copy()._values
            else:

                # we have a multi-index and are trying to align
                # with a particular, level GH3738
                if (
                    isinstance(ax, ABCMultiIndex)
                    and isinstance(df.index, ABCMultiIndex)
                    and ax.nlevels != df.index.nlevels
                ):
                    raise TypeError(
                        "cannot align on a multi-index with out "
                        "specifying the join levels"
                    )

                val = df.reindex(index=ax)._values
            return val

        raise ValueError("Incompatible indexer with DataFrame")


class _ScalarAccessIndexer(_NDFrameIndexerBase):
    """
    Access scalars quickly.
    """

    def _convert_key(self, key, is_setter: bool = False):
        raise AbstractMethodError(self)

    def __getitem__(self, key):
        if not isinstance(key, tuple):

            # we could have a convertible item here (e.g. Timestamp)
            if not is_list_like_indexer(key):
                key = tuple([key])
            else:
                raise ValueError("Invalid call for scalar access (getting)!")

        key = self._convert_key(key)
        return self.obj._get_value(*key, takeable=self._takeable)

    def __setitem__(self, key, value):
        if isinstance(key, tuple):
            key = tuple(com.apply_if_callable(x, self.obj) for x in key)
        else:
            # scalar callable may return tuple
            key = com.apply_if_callable(key, self.obj)

        if not isinstance(key, tuple):
            key = _tuplify(self.ndim, key)
        if len(key) != self.ndim:
            raise ValueError("Not enough indexers for scalar access (setting)!")
        key = list(self._convert_key(key, is_setter=True))
        self.obj._set_value(*key, value=value, takeable=self._takeable)


@Appender(IndexingMixin.at.__doc__)
class _AtIndexer(_ScalarAccessIndexer):
    _takeable = False

    def _convert_key(self, key, is_setter: bool = False):
        """
        Require they keys to be the same type as the index. (so we don't
        fallback)
        """
        # allow arbitrary setting
        if is_setter:
            return list(key)

        lkey = list(key)
        for n, (ax, i) in enumerate(zip(self.obj.axes, key)):
            lkey[n] = ax._convert_scalar_indexer(i, kind="loc")

        return tuple(lkey)


@Appender(IndexingMixin.iat.__doc__)
class _iAtIndexer(_ScalarAccessIndexer):
    _takeable = True

    def _convert_key(self, key, is_setter: bool = False):
        """
        Require integer args. (and convert to label arguments)
        """
        for a, i in zip(self.obj.axes, key):
            if not is_integer(i):
                raise ValueError("iAt based indexing can only have integer indexers")
        return key


def _tuplify(ndim: int, loc: Hashable) -> Tuple[Union[Hashable, slice], ...]:
    """
    Given an indexer for the first dimension, create an equivalent tuple
    for indexing over all dimensions.

    Parameters
    ----------
    ndim : int
    loc : object

    Returns
    -------
    tuple
    """
    _tup: List[Union[Hashable, slice]]
    _tup = [slice(None, None) for _ in range(ndim)]
    _tup[0] = loc
    return tuple(_tup)


def convert_to_index_sliceable(obj, key):
    """
    If we are index sliceable, then return my slicer, otherwise return None.
    """
    idx = obj.index
    if isinstance(key, slice):
        return idx._convert_slice_indexer(key, kind="getitem")

    elif isinstance(key, str):

        # we are an actual column
        if key in obj._data.items:
            return None

        # We might have a datetimelike string that we can translate to a
        # slice here via partial string indexing
        if idx._supports_partial_string_indexing:
            try:
                return idx._get_string_slice(key)
            except (KeyError, ValueError, NotImplementedError):
                return None

    return None


def check_bool_indexer(index: Index, key) -> np.ndarray:
    """
    Check if key is a valid boolean indexer for an object with such index and
    perform reindexing or conversion if needed.

    This function assumes that is_bool_indexer(key) == True.

    Parameters
    ----------
    index : Index
        Index of the object on which the indexing is done.
    key : list-like
        Boolean indexer to check.

    Returns
    -------
    np.array
        Resulting key.

    Raises
    ------
    IndexError
        If the key does not have the same length as index.
    IndexingError
        If the index of the key is unalignable to index.
    """
    result = key
    if isinstance(key, ABCSeries) and not key.index.equals(index):
        result = result.reindex(index)
        mask = isna(result._values)
        if mask.any():
            raise IndexingError(
                "Unalignable boolean Series provided as "
                "indexer (index of the boolean Series and of "
                "the indexed object do not match)."
            )
        result = result.astype(bool)._values
    elif is_object_dtype(key):
        # key might be object-dtype bool, check_array_indexer needs bool array
        result = np.asarray(result, dtype=bool)
        result = check_array_indexer(index, result)
    else:
        result = check_array_indexer(index, result)

    return result


def convert_missing_indexer(indexer):
    """
    Reverse convert a missing indexer, which is a dict
    return the scalar indexer and a boolean indicating if we converted
    """
    if isinstance(indexer, dict):

        # a missing key (but not a tuple indexer)
        indexer = indexer["key"]

        if isinstance(indexer, bool):
            raise KeyError("cannot use a single bool to index into setitem")
        return indexer, True

    return indexer, False


def convert_from_missing_indexer_tuple(indexer, axes):
    """
    Create a filtered indexer that doesn't have any missing indexers.
    """

    def get_indexer(_i, _idx):
        return axes[_i].get_loc(_idx["key"]) if isinstance(_idx, dict) else _idx

    return tuple(get_indexer(_i, _idx) for _i, _idx in enumerate(indexer))


def maybe_convert_ix(*args):
    """
    We likely want to take the cross-product.
    """
    ixify = True
    for arg in args:
        if not isinstance(arg, (np.ndarray, list, ABCSeries, Index)):
            ixify = False

    if ixify:
        return np.ix_(*args)
    else:
        return args


def is_nested_tuple(tup, labels) -> bool:
    """
    Returns
    -------
    bool
    """
    # check for a compatible nested tuple and multiindexes among the axes
    if not isinstance(tup, tuple):
        return False

    for i, k in enumerate(tup):

        if is_list_like(k) or isinstance(k, slice):
            return isinstance(labels, ABCMultiIndex)

    return False


def is_label_like(key) -> bool:
    """
    Returns
    -------
    bool
    """
    # select a label or row
    return not isinstance(key, slice) and not is_list_like_indexer(key)


def need_slice(obj) -> bool:
    """
    Returns
    -------
    bool
    """
    return (
        obj.start is not None
        or obj.stop is not None
        or (obj.step is not None and obj.step != 1)
    )


def _non_reducing_slice(slice_):
    """
    Ensurse that a slice doesn't reduce to a Series or Scalar.

    Any user-paseed `subset` should have this called on it
    to make sure we're always working with DataFrames.
    """
    # default to column slice, like DataFrame
    # ['A', 'B'] -> IndexSlices[:, ['A', 'B']]
    kinds = (ABCSeries, np.ndarray, Index, list, str)
    if isinstance(slice_, kinds):
        slice_ = IndexSlice[:, slice_]

    def pred(part) -> bool:
        """
        Returns
        -------
        bool
            True if slice does *not* reduce,
            False if `part` is a tuple.
        """
        # true when slice does *not* reduce, False when part is a tuple,
        # i.e. MultiIndex slice
        return (isinstance(part, slice) or is_list_like(part)) and not isinstance(
            part, tuple
        )

    if not is_list_like(slice_):
        if not isinstance(slice_, slice):
            # a 1-d slice, like df.loc[1]
            slice_ = [[slice_]]
        else:
            # slice(a, b, c)
            slice_ = [slice_]  # to tuplize later
    else:
        slice_ = [part if pred(part) else [part] for part in slice_]
    return tuple(slice_)


def _maybe_numeric_slice(df, slice_, include_bool=False):
    """
    Want nice defaults for background_gradient that don't break
    with non-numeric data. But if slice_ is passed go with that.
    """
    if slice_ is None:
        dtypes = [np.number]
        if include_bool:
            dtypes.append(bool)
        slice_ = IndexSlice[:, df.select_dtypes(include=dtypes).columns]
    return slice_


def _can_do_equal_len(labels, value, plane_indexer, lplane_indexer, obj) -> bool:
    """
    Returns
    -------
    bool
        True if we have an equal len settable.
    """
    if not len(labels) == 1 or not np.iterable(value) or is_scalar(plane_indexer[0]):
        return False

    item = labels[0]
    index = obj[item].index

    values_len = len(value)
    # equal len list/ndarray
    if len(index) == values_len:
        return True
    elif lplane_indexer == values_len:
        return True

    return False
