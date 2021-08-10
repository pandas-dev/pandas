from __future__ import annotations
from unittest.mock import PropertyMock

from pandas.util._decorators import doc
import numpy as np


class GroupByIndexingMixin:
    """
    Mixin for adding .iloc to GroupBy.
    """

    @property
    def iloc(self) -> _ilocGroupByIndexer:
        """
        Integer location-based indexing for selection by position per group.

        Similar to ``.apply(lambda x: x.iloc[i:j, k:l])``, but much faster and returns
        a subset of rows from the original DataFrame with the original index and order
        preserved.

        The output is compatible with head() and tail()
        The output is different from take() and nth() which do not preserve the index or order

        Inputs
        ------
        Allowed inputs for the first index are:

        - An integer, e.g. ``5``.
        - A slice object with ints and positive step, e.g. ``1:``, ``4:-3:2``.

        Allowed inputs for the second index are as for DataFrame.iloc, namely:

        - An integer, e.g. ``5``.
        - A list or array of integers, e.g. ``[4, 3, 0]``.
        - A slice object with ints, e.g. ``1:7``.
        - A boolean array.
        - A ``callable`` function with one argument (the calling Series or
          DataFrame) and that returns valid output for indexing (one of the above).

        Returns
        -------
        Series or DataFrame

        Note
        ----
        Neither GroupBy.nth() nor GroupBy.take() take a slice argument and
        neither of them preserve the original DataFrame order and index.
        They are both slow for large integer lists and take() is very slow for large group counts.

        Use Case
        --------
        Supose that we have a multi-indexed DataFrame with a large primary index and a secondary sorted
        to a different order for each primary.
        To reduce the DataFrame to a middle slice of each secondary, group by the primary and then
        use iloc.
        This preserves the original DataFrame's order and indexing.
        (See tests/groupby/test_groupby_iloc)

        Examples
        --------
        >>> df = pd.DataFrame([['a', 1], ['a', 2], ['a', 3], ['b', 4], ['b', 5]],
        ...                   columns=['A', 'B'])
        >>> df.groupby('A').iloc[1:2]
           A  B
        1  a  2
        4  b  5
        >>> df.groupby('A').iloc[:-1, -1:]
           B
        0  1
        1  2
        3  4
        """
        return _ilocGroupByIndexer(self)


@doc(GroupByIndexingMixin.iloc)
class _ilocGroupByIndexer:
    def __init__(self, grouped):
        self.grouped = grouped
        self.reversed = False
        self._cached_ascending_count = None
        self._cached_descending_count = None

    def __getitem__(self, arg):
        self.reversed = False

        if type(arg) == tuple:
            return self._handle_item(arg[0], arg[1])

        else:
            return self._handle_item(arg, None)

    def _handle_item(self, arg0, arg1):
        typeof_arg = type(arg0)

        if typeof_arg == slice:
            start = arg0.start
            stop = arg0.stop
            step = arg0.step

            if step is not None and step < 0:
                raise ValueError(
                    f'GroupBy.iloc row slice step must be positive. Slice was {start}:{stop}:{step}'
                )
                # self.reversed = True
                # start = None if start is None else -start - 1
                # stop = None if stop is None else -stop - 1
                # step = -step

            return self._handle_slice(start, stop, step, arg1)

        elif typeof_arg == int:
            return self._handle_slice(arg0, arg0 + 1, 1, arg1)

        else:
            raise ValueError(
                f'GroupBy.iloc row must be an integer or a slice, not a {typeof_arg}'
            )

    def _handle_slice(self, start, stop, step, arg1):
        mask = None
        if step is None:
            step = 1

        self.grouped._reset_group_selection()

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

                if step > 1:
                    #
                    # if start is -ve and -start excedes the length of a group
                    # then step must count from the
                    # first row of that group rather than the calculated offset
                    #
                    # count_array + reverse_array gives the length of the
                    # current group enabling to switch between
                    # the offset_array and the count_array depening on whether
                    #  -start excedes the group size
                    #
                    offset_array = self._descending_count + start + 1
                    limit_array = (self._ascending_count + self._descending_count + (start + 1)) < 0
                    offset_array = np.where(
                        limit_array, self._ascending_count, offset_array
                    )

                    mask &= offset_array % step == 0

        if stop is not None:
            if stop >= 0:
                if mask is None:
                    mask = self._ascending_count < stop

                else:
                    mask &= self._ascending_count < stop
            else:
                if mask is None:
                    mask = self._descending_count >= -stop

                else:
                    mask &= self._descending_count >= -stop

        if mask is None:
            arg0 = slice(None)

        else:
            arg0 = mask

        if arg1 is None:
            return self._selected_obj.iloc[arg0]

        else:
            return self._selected_obj.iloc[arg0, arg1]

    @property
    def _ascending_count(self):
        if self._cached_ascending_count is None:
            self._cached_ascending_count = self.grouped._cumcount_array()
            if self.reversed:
                self._cached_ascending_count = self._cached_ascending_count[::-1]

        return self._cached_ascending_count

    @property
    def _descending_count(self):
        if self._cached_descending_count is None:
            self._cached_descending_count = self.grouped._cumcount_array(
                ascending=False
            )
            if self.reversed:
                self._cached_descending_count = self._cached_descending_count[::-1]

        return self._cached_descending_count

    @property
    def _selected_obj(self):
        if self.reversed:
            return self.grouped._selected_obj.iloc[::-1]

        else:
            return self.grouped._selected_obj
