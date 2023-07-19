"""
Accessor namespace for Series/DataFrame methods that operate on the
index/columns and not the data itself.
"""
from __future__ import annotations

from datetime import timedelta
from typing import (
    TYPE_CHECKING,
    Generic,
)

from pandas._config import using_copy_on_write

from pandas._typing import NDFrameT

from pandas.core.indexes.api import (
    DatetimeIndex,
    MultiIndex,
    PeriodIndex,
)

if TYPE_CHECKING:
    from collections.abc import Sequence

    from pandas._typing import (
        Axis,
        Frequency,
        IndexLabel,
        Level,
        TimeAmbiguous,
        TimeNonexistent,
        ToTimestampHow,
    )


class AxisOps(Generic[NDFrameT]):
    """
    Accessor namespace for Series/DataFrame methods that operate on the
    index/columns and not the data itself.
    """

    def __init__(self, obj: NDFrameT) -> None:
        self.obj = obj

    def reorder_levels(self, order: Sequence[Level], axis: Axis = 0) -> NDFrameT:
        """
        Rearrange index levels using input order. May not drop or duplicate levels.

        Parameters
        ----------
        order : list of int or list of str
            List representing new level order. Reference level by number
            (position) or by key (label).
        axis : {0 or 'index', 1 or 'columns'}, default 0
            Where to reorder levels.

        Returns
        -------
        Series/DataFrame

        Examples
        --------
        >>> data = {
        ...     "class": ["Mammals", "Mammals", "Reptiles"],
        ...     "diet": ["Omnivore", "Carnivore", "Carnivore"],
        ...     "species": ["Humans", "Dogs", "Snakes"],
        ... }
        >>> df = pd.DataFrame(data, columns=["class", "diet", "species"])
        >>> df = df.set_index(["class", "diet"])
        >>> df
                                          species
        class      diet
        Mammals    Omnivore                Humans
                   Carnivore                 Dogs
        Reptiles   Carnivore               Snakes

        Let's reorder the levels of the index:

        >>> df.axis_ops.reorder_levels(["diet", "class"])
                                          species
        diet      class
        Omnivore  Mammals                  Humans
        Carnivore Mammals                    Dogs
                  Reptiles                 Snakes
        """
        obj = self.obj

        axis = obj._get_axis_number(axis)
        if not isinstance(obj._get_axis(axis), MultiIndex):  # pragma: no cover
            raise TypeError("Can only reorder levels on a hierarchical axis.")

        result = obj.copy(deep=None)

        if axis == 0:
            assert isinstance(result.index, MultiIndex)
            result.index = result.index.reorder_levels(order)
        else:
            assert isinstance(result.columns, MultiIndex)
            result.columns = result.columns.reorder_levels(order)
        return result

    def drop_level(self, level: IndexLabel, axis: Axis = 0) -> NDFrameT:
        """
        Return Series/Dataframe with requested index / column level(s) removed.

        Parameters
        ----------
        level : int, str, or list-like
            If a string is given, must be the name of a level
            If list-like, elements must be names or positional indexes
            of levels.

        axis : {{0 or 'index', 1 or 'columns'}}, default 0
            Axis along which the level(s) is removed:

            * 0 or 'index': remove level(s) in column.
            * 1 or 'columns': remove level(s) in row.

            For `Series` this parameter is unused and defaults to 0.

        Returns
        -------
        Series/Dataframe
            Series/Dataframe with requested index / column level(s) removed.

        Examples
        --------
        >>> df = pd.DataFrame([
        ...     [1, 2, 3, 4],
        ...     [5, 6, 7, 8],
        ...     [9, 10, 11, 12]
        ... ]).set_index([0, 1]).rename_axis(['a', 'b'])

        >>> df.columns = pd.MultiIndex.from_tuples([
        ...     ('c', 'e'), ('d', 'f')
        ... ], names=['level_1', 'level_2'])

        >>> df
        level_1   c   d
        level_2   e   f
        a b
        1 2      3   4
        5 6      7   8
        9 10    11  12

        >>> df.axis_ops.drop_level('a')
        level_1   c   d
        level_2   e   f
        b
        2        3   4
        6        7   8
        10      11  12

        >>> df.axis_ops.drop_level('level_2', axis=1)
        level_1   c   d
        a b
        1 2      3   4
        5 6      7   8
        9 10    11  12
        """
        obj = self.obj

        labels = obj._get_axis(axis)
        new_labels = labels.droplevel(level)
        return obj.set_axis(new_labels, axis=axis, copy=None)

    def swap_level(
        self, i: Level = -2, j: Level = -1, *, axis: Axis = 0, copy: bool | None = None
    ) -> NDFrameT:
        """
        Swap levels i and j in a :class:`MultiIndex`.

        Default is to swap the two innermost levels of the index.

        Parameters
        ----------
        i, j : int or str
            Levels of the indices to be swapped. Can pass level name as string.
        axis : {0 or 'index', 1 or 'columns'}, default 0
            The axis to swap levels on. 0 or 'index' for row-wise, 1 or
            'columns' for column-wise.

        Returns
        -------
        Series/DataFrame
            Series/DataFrame with levels swapped in MultiIndex.

        Examples
        --------
        >>> df = pd.DataFrame(
        ...     {"Grade": ["A", "B", "A", "C"]},
        ...     index=[
        ...         ["Final exam", "Final exam", "Coursework", "Coursework"],
        ...         ["History", "Geography", "History", "Geography"],
        ...         ["January", "February", "March", "April"],
        ...     ],
        ... )
        >>> df
                                            Grade
        Final exam  History     January      A
                    Geography   February     B
        Coursework  History     March        A
                    Geography   April        C

        In the following example, we will swap the levels of the indices.
        Here, we will swap the levels column-wise, but levels can be swapped row-wise
        in a similar manner. Note that column-wise is the default behaviour.
        By not supplying any arguments for i and j, we swap the last and second to
        last indices.

        >>> df.axis_ops.swap_level()
                                            Grade
        Final exam  January     History         A
                    February    Geography       B
        Coursework  March       History         A
                    April       Geography       C

        By supplying one argument, we can choose which index to swap the last
        index with. We can for example swap the first index with the last one as
        follows.

        >>> df.axis_ops.swap_level(0)
                                            Grade
        January     History     Final exam      A
        February    Geography   Final exam      B
        March       History     Coursework      A
        April       Geography   Coursework      C

        We can also define explicitly which indices we want to swap by supplying values
        for both i and j. Here, we for example swap the first and second indices.

        >>> df.axis_ops.swap_level(0, 1)
                                            Grade
        History     Final exam  January         A
        Geography   Final exam  February        B
        History     Coursework  March           A
        Geography   Coursework  April           C
        """
        obj = self.obj
        result = obj.copy(deep=copy and not using_copy_on_write())

        axis = obj._get_axis_number(axis)

        if not isinstance(result._get_axis(axis), MultiIndex):  # pragma: no cover
            raise TypeError("Can only swap levels on a hierarchical axis.")

        if axis == 0:
            assert isinstance(result.index, MultiIndex)
            result.index = result.index.swaplevel(i, j)
        else:
            assert isinstance(result.columns, MultiIndex)
            result.columns = result.columns.swaplevel(i, j)
        return result

    def tz_convert(
        self, tz, axis: Axis = 0, level=None, copy: bool | None = None
    ) -> NDFrameT:
        """
        Convert tz-aware axis to target time zone.

        Parameters
        ----------
        tz : str or tzinfo object or None
            Target time zone. Passing ``None`` will convert to
            UTC and remove the timezone information.
        axis : {{0 or 'index', 1 or 'columns'}}, default 0
            The axis to convert
        level : int, str, default None
            If axis is a MultiIndex, convert a specific level. Otherwise
            must be None.
        copy : bool, default True
            Also make a copy of the underlying data.

        Returns
        -------
        Series/DataFrame
            Object with time zone converted axis.

        Raises
        ------
        TypeError
            If the axis is tz-naive.

        Examples
        --------
        Change to another time zone:

        >>> ser = pd.Series(
        ...       [1],
        ...       index=pd.DatetimeIndex(['2018-09-15 01:30:00+02:00']),
        ... )
        >>> ser.axis_ops.tz_convert('Asia/Shanghai')
        2018-09-15 07:30:00+08:00    1
        dtype: int64

        Pass None to convert to UTC and get a tz-naive index:

        >>> ser = pd.Series([1],
        ...       index=pd.DatetimeIndex(['2018-09-15 01:30:00+02:00']))
        >>> ser.axis_ops.tz_convert(None)
        2018-09-14 23:30:00    1
        dtype: int64
        """
        obj = self.obj

        axis = obj._get_axis_number(axis)
        ax = obj._get_axis(axis)

        def _tz_convert(ax, tz):
            if not hasattr(ax, "tz_convert"):
                if len(ax) > 0:
                    ax_name = obj._get_axis_name(axis)
                    raise TypeError(
                        f"{ax_name} is not a valid DatetimeIndex or PeriodIndex"
                    )
                ax = DatetimeIndex([], tz=tz)
            else:
                ax = ax.tz_convert(tz)
            return ax

        # if a level is given it must be a MultiIndex level or
        # equivalent to the axis name
        if isinstance(ax, MultiIndex):
            level = ax._get_level_number(level)
            new_level = _tz_convert(ax.levels[level], tz)
            ax = ax.set_levels(new_level, level=level)
        else:
            if level not in (None, 0, ax.name):
                raise ValueError(f"The level {level} is not valid")
            ax = _tz_convert(ax, tz)

        result = obj.copy(deep=copy and not using_copy_on_write())
        result = result.set_axis(ax, axis=axis, copy=False)
        return result.__finalize__(obj, method="tz_convert")

    def tz_localize(
        self,
        tz,
        axis: Axis = 0,
        level=None,
        copy: bool | None = None,
        ambiguous: TimeAmbiguous = "raise",
        nonexistent: TimeNonexistent = "raise",
    ) -> NDFrameT:
        """
        Localize tz-naive index of a Series or DataFrame to target time zone.

        This operation localizes the Index. To localize the values in a
        timezone-naive Series, use :meth:`Series.dt.tz_localize`.

        Parameters
        ----------
        tz : str or tzinfo or None
            Time zone to localize. Passing ``None`` will remove the
            time zone information and preserve local time.
        axis : {{0 or 'index', 1 or 'columns'}}, default 0
            The axis to localize
        level : int, str, default None
            If axis ia a MultiIndex, localize a specific level. Otherwise
            must be None.
        copy : bool, default True
            Also make a copy of the underlying data.
        ambiguous : 'infer', bool-ndarray, 'NaT', default 'raise'
            When clocks moved backward due to DST, ambiguous times may arise.
            For example in Central European Time (UTC+01), when going from
            03:00 DST to 02:00 non-DST, 02:30:00 local time occurs both at
            00:30:00 UTC and at 01:30:00 UTC. In such a situation, the
            `ambiguous` parameter dictates how ambiguous times should be
            handled.

            - 'infer' will attempt to infer fall dst-transition hours based on
              order
            - bool-ndarray where True signifies a DST time, False designates
              a non-DST time (note that this flag is only applicable for
              ambiguous times)
            - 'NaT' will return NaT where there are ambiguous times
            - 'raise' will raise an AmbiguousTimeError if there are ambiguous
              times.
        nonexistent : str, default 'raise'
            A nonexistent time does not exist in a particular timezone
            where clocks moved forward due to DST. Valid values are:

            - 'shift_forward' will shift the nonexistent time forward to the
              closest existing time
            - 'shift_backward' will shift the nonexistent time backward to the
              closest existing time
            - 'NaT' will return NaT where there are nonexistent times
            - timedelta objects will shift nonexistent times by the timedelta
            - 'raise' will raise an NonExistentTimeError if there are
              nonexistent times.

        Returns
        -------
        Series/DataFrame
            Same type as the input.

        Raises
        ------
        TypeError
            If the TimeSeries is tz-aware and tz is not None.

        Examples
        --------
        Localize local times:

        >>> ser = pd.Series(
        ...       [1],
        ...       index=pd.DatetimeIndex(['2018-09-15 01:30:00']),
        ... )
        >>> ser.axis_ops.tz_localize('CET')
        2018-09-15 01:30:00+02:00    1
        dtype: int64

        Pass None to convert to tz-naive index and preserve local time:

        >>> ser = pd.Series([1],
        ...       index=pd.DatetimeIndex(['2018-09-15 01:30:00+02:00']))
        >>> ser.axis_ops.tz_localize(None)
        2018-09-15 01:30:00    1
        dtype: int64

        Be careful with DST changes. When there is sequential data, pandas
        can infer the DST time:

        >>> ser = pd.Series(range(7),
        ...                 index=pd.DatetimeIndex(['2018-10-28 01:30:00',
        ...                                         '2018-10-28 02:00:00',
        ...                                         '2018-10-28 02:30:00',
        ...                                         '2018-10-28 02:00:00',
        ...                                         '2018-10-28 02:30:00',
        ...                                         '2018-10-28 03:00:00',
        ...                                         '2018-10-28 03:30:00']))
        >>> ser.axis_ops.tz_localize('CET', ambiguous='infer')
        2018-10-28 01:30:00+02:00    0
        2018-10-28 02:00:00+02:00    1
        2018-10-28 02:30:00+02:00    2
        2018-10-28 02:00:00+01:00    3
        2018-10-28 02:30:00+01:00    4
        2018-10-28 03:00:00+01:00    5
        2018-10-28 03:30:00+01:00    6
        dtype: int64

        In some cases, inferring the DST is impossible. In such cases, you can
        pass an ndarray to the ambiguous parameter to set the DST explicitly

        >>> ser = pd.Series(range(3),
        ...                 index=pd.DatetimeIndex(['2018-10-28 01:20:00',
        ...                                         '2018-10-28 02:36:00',
        ...                                         '2018-10-28 03:46:00']))
        >>> ser.axis_ops.tz_localize('CET', ambiguous=np.array([True, True, False]))
        2018-10-28 01:20:00+02:00    0
        2018-10-28 02:36:00+02:00    1
        2018-10-28 03:46:00+01:00    2
        dtype: int64

        If the DST transition causes nonexistent times, you can shift these
        dates forward or backward with a timedelta object or `'shift_forward'`
        or `'shift_backward'`.

        >>> ser = pd.Series(range(2),
        ...                 index=pd.DatetimeIndex(['2015-03-29 02:30:00',
        ...                                         '2015-03-29 03:30:00']))
        >>> ser.axis_ops.tz_localize('Europe/Warsaw', nonexistent='shift_forward')
        2015-03-29 03:00:00+02:00    0
        2015-03-29 03:30:00+02:00    1
        dtype: int64
        >>> ser.axis_ops.tz_localize('Europe/Warsaw', nonexistent='shift_backward')
        2015-03-29 01:59:59.999999999+01:00    0
        2015-03-29 03:30:00+02:00              1
        dtype: int64
        >>> ser.axis_ops.tz_localize('Europe/Warsaw', nonexistent=pd.Timedelta('1H'))
        2015-03-29 03:30:00+02:00    0
        2015-03-29 03:30:00+02:00    1
        dtype: int64
        """
        obj = self.obj

        nonexistent_options = ("raise", "NaT", "shift_forward", "shift_backward")
        if nonexistent not in nonexistent_options and not isinstance(
            nonexistent, timedelta
        ):
            raise ValueError(
                "The nonexistent argument must be one of 'raise', "
                "'NaT', 'shift_forward', 'shift_backward' or "
                "a timedelta object"
            )

        axis = obj._get_axis_number(axis)
        ax = obj._get_axis(axis)

        def _tz_localize(ax, tz, ambiguous, nonexistent):
            if not hasattr(ax, "tz_localize"):
                if len(ax) > 0:
                    ax_name = obj._get_axis_name(axis)
                    raise TypeError(
                        f"{ax_name} is not a valid DatetimeIndex or PeriodIndex"
                    )
                ax = DatetimeIndex([], tz=tz)
            else:
                ax = ax.tz_localize(tz, ambiguous=ambiguous, nonexistent=nonexistent)
            return ax

        # if a level is given it must be a MultiIndex level or
        # equivalent to the axis name
        if isinstance(ax, MultiIndex):
            level = ax._get_level_number(level)
            new_level = _tz_localize(ax.levels[level], tz, ambiguous, nonexistent)
            ax = ax.set_levels(new_level, level=level)
        else:
            if level not in (None, 0, ax.name):
                raise ValueError(f"The level {level} is not valid")
            ax = _tz_localize(ax, tz, ambiguous, nonexistent)

        result = obj.copy(deep=copy and not using_copy_on_write())
        result = result.set_axis(ax, axis=axis, copy=False)
        return result.__finalize__(obj, method="tz_localize")

    def to_timestamp(
        self,
        freq=None,
        *,
        how: ToTimestampHow = "start",
        axis: Axis = 0,
        copy: bool | None = None,
    ) -> NDFrameT:
        """
        Cast to DatetimeIndex of Timestamps, at *beginning* of period.

        Parameters
        ----------
        freq : str, default frequency of PeriodIndex
            Desired frequency.
        how : {'s', 'e', 'start', 'end'}
            Convention for converting period to timestamp; start of period
            vs. end.
        axis : {0 or 'index', 1 or 'columns'}, default 0
        copy : bool, default True
            Whether or not to return a copy.

        Returns
        -------
        Series with DatetimeIndex

        Examples
        --------
        >>> idx = pd.PeriodIndex(['2023', '2024', '2025'], freq='Y')
        >>> s1 = pd.Series([1, 2, 3], index=idx)
        >>> s1
        2023    1
        2024    2
        2025    3
        Freq: A-DEC, dtype: int64

        The resulting frequency of the Timestamps is `YearBegin`

        >>> s1 = s1.axis_ops.to_timestamp()
        >>> s1
        2023-01-01    1
        2024-01-01    2
        2025-01-01    3
        Freq: AS-JAN, dtype: int64

        Using `freq` which is the offset that the Timestamps will have

        >>> s2 = pd.Series([1, 2, 3], index=idx)
        >>> s2 = s2.axis_ops.to_timestamp(freq='M')
        >>> s2
        2023-01-31    1
        2024-01-31    2
        2025-01-31    3
        Freq: A-JAN, dtype: int64
        """
        obj = self.obj

        axis_name = obj._get_axis_name(axis)
        old_ax = getattr(obj, axis_name)

        if not isinstance(old_ax, PeriodIndex):
            raise TypeError(f"unsupported Type {type(old_ax).__name__}")

        new_obj = obj.copy(deep=copy and not using_copy_on_write())
        new_ax = old_ax.to_timestamp(freq=freq, how=how)
        setattr(new_obj, axis_name, new_ax)
        return new_obj

    def to_period(
        self, freq: Frequency | None = None, *, axis: Axis = 0, copy: bool | None = None
    ) -> NDFrameT:
        """
        Convert Series/DataFrame from DatetimeIndex to PeriodIndex.

        Parameters
        ----------
        freq : str, default None
            Frequency associated with the PeriodIndex.
        axis : {0 or 'index', 1 or 'columns'}, default 0
        copy : bool, default True
            Whether or not to return a copy.

        Returns
        -------
        Series
            Series with index converted to PeriodIndex.

        Examples
        --------
        >>> idx = pd.DatetimeIndex(['2023', '2024', '2025'])
        >>> ser = pd.Series([1, 2, 3], index=idx)
        >>> ser = ser.axis_ops.to_period()
        >>> ser
        2023    1
        2024    2
        2025    3
        Freq: A-DEC, dtype: int64

        Viewing the index

        >>> ser.index
        PeriodIndex(['2023', '2024', '2025'], dtype='period[A-DEC]')
        """
        obj = self.obj

        axis_name = obj._get_axis_name(axis)
        old_ax = getattr(obj, axis_name)

        if not isinstance(old_ax, DatetimeIndex):
            raise TypeError(f"unsupported Type {type(old_ax).__name__}")

        new_obj = obj.copy(deep=copy and not using_copy_on_write())
        new_ax = old_ax.to_period(freq=freq)
        setattr(new_obj, axis_name, new_ax)
        return new_obj
