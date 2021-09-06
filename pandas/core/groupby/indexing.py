from __future__ import annotations

import numpy as np

from pandas.util._decorators import doc


class GroupByIndexingMixin:
    """
    Mixin for adding .rows to GroupBy.
    """

    @property
    def rows(self) -> _rowsGroupByIndexer:
        """
        Purely integer-location based indexing for selection by position per group.

        ``.rows[]`` is primarily integer position based (from ``0`` to
        ``length-1`` of the axis),

        Allowed inputs for the index are:

        - An integer valued iterable, e.g. ``range(2, 4)``.
        - A comma separated list of integers and slices, e.g. ``5``, ``2, 4``, ``2:4``.

        Note: the slice step cannot be negative.

        The output format is the same as GroupBy.head and GroupBy.tail, namely a subset
        of the original DataFrame or Series with the index and order preserved.

        The effect of ``grouped.rows[i:j]`` is similar to

            ``grouped.apply(lambda x: x.iloc[i:j])``

        but much faster and preserving the original index and order.

        The behaviour is different from GroupBy.nth:

        - Input to rows can include one or more slices whereas nth just handles
          a list of indexes.
        - Output from rows is in the same order as the original grouped DataFrame
          or Series.
        - Output from rows has the same index columns as the original grouped DataFrame
          or Series. (nth behaves like an aggregator and removes the non-grouped
          indexes.)
        - GroupBy.rows can  define a slice relative to the last row of each group.
        - GroupBy.rows is faster than nth.
        - GroupBy.rows does not handle dropna.

        An important use case for GroupBy.rows is a multi-indexed DataFrame with a
        large primary index (Date, say) and a secondary index sorted to a different
        order for each Date.
        To reduce the DataFrame to a middle slice of each Date:

            ``df.groupby("Date").rows[5:-5]``

        This returns a subset of df containing just the middle rows for each Date
        and with its original order and indexing preserved.

        To reduce the DataFrame to the remaining rows:

            ``df.groupby("Date").rows[:5, -5:]``

        Returns
        -------
        Series
            The filtered subset of the original grouped Series.
        DataFrame
            The filtered subset of the original grouped DataFrame.

        See Also
        --------
        DataFrame.iloc : Purely integer-location based indexing for selection by
            position.
        GroupBy.head : Return first n rows of each group.
        GroupBy.tail : Return last n rows of each group.
        GroupBy.nth : Take the nth row from each group if n is an int, or a
            subset of rows, if n is a list of ints.

        Examples
        --------
            >>> df = pd.DataFrame([["a", 1], ["a", 2], ["a", 3], ["b", 4], ["b", 5]],
            ...                   columns=["A", "B"])
            >>> df.groupby("A").rows[1:2]
               A  B
            1  a  2
            4  b  5

            >>> df.groupby("A").rows[1, -1]
               A  B
            1  a  2
            2  a  3
            4  b  5
        """
        return _rowsGroupByIndexer(self)


@doc(GroupByIndexingMixin.rows)
class _rowsGroupByIndexer:
    def __init__(self, grouped):
        self.grouped = grouped

    def __getitem__(self, arg):
        self.grouped._reset_group_selection()
        self._cached_ascending_count = None
        self._cached_descending_count = None

        if isinstance(arg, tuple):
            mask = self._handle_tuple(arg)

        elif isinstance(arg, slice):
            mask = self._handle_slice(arg)

        elif isinstance(arg, int):
            mask = self._handle_int(arg)

        elif isinstance(arg, list):
            mask = self._handle_list(arg)

        else:
            try:
                list_arg = list(arg)

            except TypeError:
                raise ValueError(
                    f"Invalid index {type(arg)}. Must be iterable or a list of "
                    "integers and slices"
                )

            mask = self._handle_list(list_arg)

        if mask is None or mask is True:
            mask = slice(None)

        return self.grouped._selected_obj.iloc[mask]

    def _handle_int(self, arg):
        if arg >= 0:
            return self._ascending_count == arg

        else:
            return self._descending_count == (-arg - 1)

    def _handle_list(self, args):
        positive = [arg for arg in args if arg >= 0]
        negative = [-arg - 1 for arg in args if arg < 0]

        if positive:
            mask = np.isin(self._ascending_count, positive)

        else:
            mask = False

        if negative:
            mask |= np.isin(self._descending_count, negative)

        return mask

    def _handle_tuple(self, args):
        mask = False

        for arg in args:
            if isinstance(arg, int):
                mask |= self._handle_int(arg)

            elif isinstance(arg, slice):
                mask |= self._handle_slice(arg)

            else:
                raise ValueError(
                    f"Invalid argument {type(arg)}. Should be int or slice."
                )

        return mask

    def _handle_slice(self, arg):
        start = arg.start
        stop = arg.stop
        step = arg.step

        if step is not None and step < 0:
            raise ValueError(f"Invalid step {step}. Must be non-negative")

        mask = True
        if step is None:
            step = 1

        if start is None:
            if step > 1:
                mask = self._ascending_count % step == 0

        else:
            if start >= 0:
                mask = self._ascending_count >= start

                if step > 1:
                    mask &= (self._ascending_count - start) % step == 0

            else:
                mask = self._descending_count < -start

                offset_array = self._descending_count + start + 1
                limit_array = (
                    self._ascending_count + self._descending_count + (start + 1)
                ) < 0
                offset_array = np.where(
                    limit_array, self._ascending_count, offset_array
                )

                mask &= offset_array % step == 0

        if stop is not None:
            if stop >= 0:
                mask &= self._ascending_count < stop

            else:
                mask &= self._descending_count >= -stop

        return mask

    @property
    def _ascending_count(self):
        if self._cached_ascending_count is None:
            self._cached_ascending_count = self.grouped._cumcount_array()

        return self._cached_ascending_count

    @property
    def _descending_count(self):
        if self._cached_descending_count is None:
            self._cached_descending_count = self.grouped._cumcount_array(
                ascending=False
            )

        return self._cached_descending_count
