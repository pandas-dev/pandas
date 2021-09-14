from __future__ import annotations

from typing import (
    Iterable,
    cast,
)

import numpy as np

from pandas._typing import (
    FrameOrSeries,
    PositionalIndexer,
)
from pandas.util._decorators import (
    cache_readonly,
    doc,
)

from pandas.core.dtypes.common import (
    is_integer,
    is_list_like,
)

from pandas.core.groupby import groupby
from pandas.core.indexes.api import CategoricalIndex


class GroupByIndexingMixin:
    """
    Mixin for adding ._rows to GroupBy.
    """

    @property
    def _rows(self) -> RowsGroupByIndexer:
        return RowsGroupByIndexer(cast(groupby.GroupBy, self))


@doc(GroupByIndexingMixin._rows)
class RowsGroupByIndexer:
    def __init__(self, groupby_object: groupby.GroupBy):
        self.groupby_object = groupby_object

    def __getitem__(self, arg: PositionalIndexer | tuple) -> FrameOrSeries:
        """
        Positional index for selection by integer location per group.

        Used to implement GroupBy._rows which is used to implement GroupBy.nth
        in the case when the keyword dropna is None or absent.

        Parameters
        ----------
        arg : PositionalIndexer | tuple
            Allowed values are:
            - int
            - int valued iterable such as list or range
            - slice with step either None or positive
            - tuple of integers and slices

        Returns
        -------
        Series
            The filtered subset of the original groupby Series.
        DataFrame
            The filtered subset of the original groupby DataFrame.

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
        >>> df.groupby("A", as_index=False)._rows[1:2]
           A  B
        1  a  2
        4  b  5

        >>> df.groupby("A", as_index=False)._rows[1, -1]
           A  B
        1  a  2
        2  a  3
        4  b  5
        """
        with groupby.group_selection_context(self.groupby_object):
            if isinstance(arg, tuple):
                if all(is_integer(i) for i in arg):
                    mask = self._handle_list(arg)
                else:
                    mask = self._handle_tuple(arg)

            elif isinstance(arg, slice):
                mask = self._handle_slice(arg)
            elif is_integer(arg):
                mask = self._handle_int(cast(int, arg))
            elif is_list_like(arg):
                mask = self._handle_list(cast(Iterable[int], arg))
            else:
                raise TypeError(
                    f"Invalid index {type(arg)}. "
                    "Must be integer, list-like, slice or a tuple of "
                    "integers and slices"
                )

            ids, _, _ = self.groupby_object.grouper.group_info

            # Drop NA values in grouping
            mask &= ids != -1

            if mask is None or mask is True:
                result = self.groupby_object._selected_obj[:]
            else:
                result = self.groupby_object._selected_obj[mask]

            if self.groupby_object.as_index:
                result_index = self.groupby_object.grouper.result_index
                result.index = result_index[ids[mask]]

                if not self.groupby_object.observed and isinstance(
                    result_index, CategoricalIndex
                ):
                    result = result.reindex(result_index)

                result = self.groupby_object._reindex_output(result)
                if self.groupby_object.sort:
                    result = result.sort_index()

            return result

    def _handle_int(self, arg: int) -> bool | np.ndarray:
        if arg >= 0:
            return self._ascending_count == arg
        else:
            return self._descending_count == (-arg - 1)

    def _handle_list(self, args: Iterable[int]) -> bool | np.ndarray:
        positive = [arg for arg in args if arg >= 0]
        negative = [-arg - 1 for arg in args if arg < 0]

        mask: bool | np.ndarray = False

        if positive:
            mask |= np.isin(self._ascending_count, positive)

        if negative:
            mask |= np.isin(self._descending_count, negative)

        return mask

    def _handle_tuple(self, args: tuple) -> bool | np.ndarray:
        mask: bool | np.ndarray = False

        for arg in args:
            if is_integer(arg):
                mask |= self._handle_int(cast(int, arg))
            elif isinstance(arg, slice):
                mask |= self._handle_slice(arg)
            else:
                raise ValueError(
                    f"Invalid argument {type(arg)}. Should be int or slice."
                )

        return mask

    def _handle_slice(self, arg: slice) -> bool | np.ndarray:
        start = arg.start
        stop = arg.stop
        step = arg.step

        if step is not None and step < 0:
            raise ValueError(f"Invalid step {step}. Must be non-negative")

        mask: bool | np.ndarray = True

        if step is None:
            step = 1

        if start is None:
            if step > 1:
                mask &= self._ascending_count % step == 0

        elif start >= 0:
            mask &= self._ascending_count >= start

            if step > 1:
                mask &= (self._ascending_count - start) % step == 0

        else:
            mask &= self._descending_count < -start

            offset_array = self._descending_count + start + 1
            limit_array = (
                self._ascending_count + self._descending_count + (start + 1)
            ) < 0
            offset_array = np.where(limit_array, self._ascending_count, offset_array)

            mask &= offset_array % step == 0

        if stop is not None:
            if stop >= 0:
                mask &= self._ascending_count < stop
            else:
                mask &= self._descending_count >= -stop

        return mask

    @cache_readonly
    def _ascending_count(self) -> np.ndarray:
        return self.groupby_object._cumcount_array()

    @cache_readonly
    def _descending_count(self) -> np.ndarray:
        return self.groupby_object._cumcount_array(ascending=False)
