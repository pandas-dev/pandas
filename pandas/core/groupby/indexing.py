from __future__ import annotations

from typing import (
    TYPE_CHECKING,
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

if TYPE_CHECKING:
    from pandas.core.groupby import groupby


class GroupByIndexingMixin:
    """
    Mixin for adding ._body to GroupBy.
    """

    @cache_readonly
    def _body(self) -> BodyGroupByIndexer:
        if TYPE_CHECKING:
            self = cast(groupby.GroupBy, self)

        return BodyGroupByIndexer(self)

    def _make_mask(self, arg: PositionalIndexer | tuple) -> np.ndarray:
        if is_list_like(arg):
            if all(is_integer(i) for i in cast(Iterable, arg)):
                mask = self._handle_list(cast(Iterable[int], arg))
            else:
                mask = self._handle_tuple(cast(tuple, arg))

        elif isinstance(arg, slice):
            mask = self._handle_slice(arg)
        elif is_integer(arg):
            mask = self._handle_int(cast(int, arg))
        else:
            raise TypeError(
                f"Invalid index {type(arg)}. "
                "Must be integer, list-like, slice or a tuple of "
                "integers and slices"
            )

        if isinstance(mask, bool):
            if mask:
                mask = self._ascending_count >= 0
            else:
                mask = self._ascending_count < 0

        return cast(np.ndarray, mask)

    def _handle_int(self, arg: int) -> np.ndarray:
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

    def _apply_mask(self, mask: np.ndarray):
        if TYPE_CHECKING:
            self = cast(groupby.GroupBy, self)

        if self.axis == 0:
            return self._selected_obj[mask]
        else:
            return self._selected_obj.iloc[:, mask]

    @cache_readonly
    def _ascending_count(self) -> np.ndarray:
        if TYPE_CHECKING:
            self = cast(groupby.GroupBy, self)

        return self._cumcount_array()

    @cache_readonly
    def _descending_count(self) -> np.ndarray:
        if TYPE_CHECKING:
            self = cast(groupby.GroupBy, self)

        return self._cumcount_array(ascending=False)


@doc(GroupByIndexingMixin._body)
class BodyGroupByIndexer:
    def __init__(self, groupby_object: groupby.GroupBy):
        self.groupby_object = groupby_object

    def __getitem__(self, arg: PositionalIndexer | tuple) -> FrameOrSeries:
        """
        Positional index for selection by integer location per group.

        Used to implement GroupBy._body which is used to implement GroupBy.nth
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
        """
        self.groupby_object._reset_group_selection()
        mask = self.groupby_object._make_mask(arg)
        return self.groupby_object._apply_mask(mask)
