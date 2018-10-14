import textwrap
import numpy as np

from pandas._libs.interval import (Interval, IntervalMixin,
                                   intervals_to_interval_bounds)
from pandas.compat import add_metaclass
from pandas.compat.numpy import function as nv
import pandas.core.common as com
from pandas.core.config import get_option
from pandas.core.dtypes.cast import maybe_convert_platform
from pandas.core.dtypes.common import (is_categorical_dtype, is_float_dtype,
                                       is_integer_dtype, is_interval_dtype,
                                       is_scalar, is_string_dtype,
                                       is_datetime64_any_dtype,
                                       is_timedelta64_dtype, is_interval,
                                       pandas_dtype)
from pandas.core.dtypes.dtypes import IntervalDtype
from pandas.core.dtypes.generic import (ABCDatetimeIndex, ABCPeriodIndex,
                                        ABCSeries, ABCIntervalIndex,
                                        ABCInterval)
from pandas.core.dtypes.missing import isna, notna
from pandas.core.indexes.base import Index, ensure_index
from pandas.util._decorators import Appender
from pandas.util._doctools import _WritableDoc

from pandas.core.arrays import ExtensionArray, Categorical

_VALID_CLOSED = {'left', 'right', 'both', 'neither'}
_interval_shared_docs = {}
_shared_docs_kwargs = dict(
    klass='IntervalArray',
    name=''
)


_interval_shared_docs['class'] = """%(summary)s

.. versionadded:: %(versionadded)s

.. warning::

   The indexing behaviors are provisional and may change in
   a future version of pandas.

Parameters
----------
data : array-like (1-dimensional)
    Array-like containing Interval objects from which to build the
    %(klass)s.
closed : {'left', 'right', 'both', 'neither'}, default 'right'
    Whether the intervals are closed on the left-side, right-side, both or
    neither.
%(name)s\
copy : boolean, default False
    Copy the meta-data.
dtype : dtype or None, default None
    If None, dtype will be inferred

    .. versionadded:: 0.23.0

Attributes
----------
left
right
closed
mid
length
values
is_non_overlapping_monotonic

Methods
-------
from_arrays
from_tuples
from_breaks
set_closed
%(extra_methods)s\

%(examples)s\

Notes
------
See the `user guide
<http://pandas.pydata.org/pandas-docs/stable/advanced.html#intervalindex>`_
for more.

See Also
--------
Index : The base pandas Index type
Interval : A bounded slice-like interval; the elements of an IntervalIndex
interval_range : Function to create a fixed frequency IntervalIndex
cut, qcut : Convert arrays of continuous data into Categoricals/Series of
            Intervals
"""


@Appender(_interval_shared_docs['class'] % dict(
    klass="IntervalArray",
    summary="Pandas array for interval data that are closed on the same side",
    versionadded="0.24.0",
    name='', extra_methods='', examples='',
))
@add_metaclass(_WritableDoc)
class IntervalArray(IntervalMixin, ExtensionArray):
    dtype = IntervalDtype()
    ndim = 1
    can_hold_na = True
    _na_value = _fill_value = np.nan

    def __new__(cls, data, closed=None, dtype=None, copy=False,
                verify_integrity=True):

        if isinstance(data, ABCSeries) and is_interval_dtype(data):
            data = data.values

        if isinstance(data, (cls, ABCIntervalIndex)):
            left = data.left
            right = data.right
            closed = closed or data.closed
        else:

            # don't allow scalars
            if is_scalar(data):
                msg = ("{}(...) must be called with a collection of some kind,"
                       " {} was passed")
                raise TypeError(msg.format(cls.__name__, data))

            # might need to convert empty or purely na data
            data = maybe_convert_platform_interval(data)
            left, right, infer_closed = intervals_to_interval_bounds(
                data, validate_closed=closed is None)
            closed = closed or infer_closed

        return cls._simple_new(left, right, closed, copy=copy, dtype=dtype,
                               verify_integrity=verify_integrity)

    @classmethod
    def _simple_new(cls, left, right, closed=None,
                    copy=False, dtype=None, verify_integrity=True):
        result = IntervalMixin.__new__(cls)

        closed = closed or 'right'
        left = ensure_index(left, copy=copy)
        right = ensure_index(right, copy=copy)

        if dtype is not None:
            # GH 19262: dtype must be an IntervalDtype to override inferred
            dtype = pandas_dtype(dtype)
            if not is_interval_dtype(dtype):
                msg = 'dtype must be an IntervalDtype, got {dtype}'
                raise TypeError(msg.format(dtype=dtype))
            elif dtype.subtype is not None:
                left = left.astype(dtype.subtype)
                right = right.astype(dtype.subtype)

        # coerce dtypes to match if needed
        if is_float_dtype(left) and is_integer_dtype(right):
            right = right.astype(left.dtype)
        elif is_float_dtype(right) and is_integer_dtype(left):
            left = left.astype(right.dtype)

        if type(left) != type(right):
            msg = ('must not have differing left [{ltype}] and right '
                   '[{rtype}] types')
            raise ValueError(msg.format(ltype=type(left).__name__,
                                        rtype=type(right).__name__))
        elif is_categorical_dtype(left.dtype) or is_string_dtype(left.dtype):
            # GH 19016
            msg = ('category, object, and string subtypes are not supported '
                   'for IntervalArray')
            raise TypeError(msg)
        elif isinstance(left, ABCPeriodIndex):
            msg = 'Period dtypes are not supported, use a PeriodIndex instead'
            raise ValueError(msg)
        elif (isinstance(left, ABCDatetimeIndex) and
                str(left.tz) != str(right.tz)):
            msg = ("left and right must have the same time zone, got "
                   "'{left_tz}' and '{right_tz}'")
            raise ValueError(msg.format(left_tz=left.tz, right_tz=right.tz))

        result._left = left
        result._right = right
        result._closed = closed
        if verify_integrity:
            result._validate()
        return result

    @classmethod
    def _from_sequence(cls, scalars, dtype=None, copy=False):
        return cls(scalars, dtype=dtype, copy=copy)

    @classmethod
    def _from_factorized(cls, values, original):
        return cls(values, closed=original.closed)

    _interval_shared_docs['from_breaks'] = """
    Construct an %(klass)s from an array of splits.

    Parameters
    ----------
    breaks : array-like (1-dimensional)
        Left and right bounds for each interval.
    closed : {'left', 'right', 'both', 'neither'}, default 'right'
        Whether the intervals are closed on the left-side, right-side, both
        or neither.
    copy : boolean, default False
        copy the data
    dtype : dtype or None, default None
        If None, dtype will be inferred

        .. versionadded:: 0.23.0

    Examples
    --------
    >>> pd.%(klass)s.from_breaks([0, 1, 2, 3])
    %(klass)s([(0, 1], (1, 2], (2, 3]]
                  closed='right',
                  dtype='interval[int64]')

    See Also
    --------
    interval_range : Function to create a fixed frequency IntervalIndex
    %(klass)s.from_arrays : Construct from a left and right array
    %(klass)s.from_tuples : Construct from a sequence of tuples
    """

    @classmethod
    @Appender(_interval_shared_docs['from_breaks'] % _shared_docs_kwargs)
    def from_breaks(cls, breaks, closed='right', copy=False, dtype=None):
        breaks = maybe_convert_platform_interval(breaks)

        return cls.from_arrays(breaks[:-1], breaks[1:], closed, copy=copy,
                               dtype=dtype)

    _interval_shared_docs['from_arrays'] = """
        Construct from two arrays defining the left and right bounds.

        Parameters
        ----------
        left : array-like (1-dimensional)
            Left bounds for each interval.
        right : array-like (1-dimensional)
            Right bounds for each interval.
        closed : {'left', 'right', 'both', 'neither'}, default 'right'
            Whether the intervals are closed on the left-side, right-side, both
            or neither.
        copy : boolean, default False
            Copy the data.
        dtype : dtype, optional
            If None, dtype will be inferred.

            .. versionadded:: 0.23.0

        Returns
        -------
        %(klass)s

        Notes
        -----
        Each element of `left` must be less than or equal to the `right`
        element at the same position. If an element is missing, it must be
        missing in both `left` and `right`. A TypeError is raised when
        using an unsupported type for `left` or `right`. At the moment,
        'category', 'object', and 'string' subtypes are not supported.

        Raises
        ------
        ValueError
            When a value is missing in only one of `left` or `right`.
            When a value in `left` is greater than the corresponding value
            in `right`.

        See Also
        --------
        interval_range : Function to create a fixed frequency IntervalIndex.
        %(klass)s.from_breaks : Construct an %(klass)s from an array of
            splits.
        %(klass)s.from_tuples : Construct an %(klass)s from an
            array-like of tuples.


        Examples
        --------
        >>> %(klass)s.from_arrays([0, 1, 2], [1, 2, 3])
        %(klass)s([(0, 1], (1, 2], (2, 3]]
                     closed='right',
                     dtype='interval[int64]')
        """

    @classmethod
    @Appender(_interval_shared_docs['from_arrays'] % _shared_docs_kwargs)
    def from_arrays(cls, left, right, closed='right', copy=False, dtype=None):
        left = maybe_convert_platform_interval(left)
        right = maybe_convert_platform_interval(right)

        return cls._simple_new(left, right, closed, copy=copy,
                               dtype=dtype, verify_integrity=True)

    _interval_shared_docs['from_intervals'] = """
    Construct an %(klass)s from a 1d array of Interval objects

    .. deprecated:: 0.23.0

    Parameters
    ----------
    data : array-like (1-dimensional)
        Array of Interval objects. All intervals must be closed on the same
        sides.
    copy : boolean, default False
        by-default copy the data, this is compat only and ignored
    dtype : dtype or None, default None
        If None, dtype will be inferred

        ..versionadded:: 0.23.0

    Examples
    --------
    >>> pd.%(klass)s.from_intervals([pd.Interval(0, 1),
    ...                                  pd.Interval(1, 2)])
    %(klass)s([(0, 1], (1, 2]]
                  closed='right', dtype='interval[int64]')

    The generic Index constructor work identically when it infers an array
    of all intervals:

    >>> pd.Index([pd.Interval(0, 1), pd.Interval(1, 2)])
    %(klass)s([(0, 1], (1, 2]]
                  closed='right', dtype='interval[int64]')

    See Also
    --------
    interval_range : Function to create a fixed frequency IntervalIndex
    %(klass)s.from_arrays : Construct an %(klass)s from a left and
                                right array
    %(klass)s.from_breaks : Construct an %(klass)s from an array of
                                splits
    %(klass)s.from_tuples : Construct an %(klass)s from an
                                array-like of tuples
    """

    _interval_shared_docs['from_tuples'] = """
    Construct an %(klass)s from an array-like of tuples

    Parameters
    ----------
    data : array-like (1-dimensional)
        Array of tuples
    closed : {'left', 'right', 'both', 'neither'}, default 'right'
        Whether the intervals are closed on the left-side, right-side, both
        or neither.
    copy : boolean, default False
        by-default copy the data, this is compat only and ignored
    dtype : dtype or None, default None
        If None, dtype will be inferred

        ..versionadded:: 0.23.0


    Examples
    --------
    >>>  pd.%(klass)s.from_tuples([(0, 1), (1, 2)])
    %(klass)s([(0, 1], (1, 2]],
                closed='right', dtype='interval[int64]')

    See Also
    --------
    interval_range : Function to create a fixed frequency IntervalIndex
    %(klass)s.from_arrays : Construct an %(klass)s from a left and
                                right array
    %(klass)s.from_breaks : Construct an %(klass)s from an array of
                                splits
    """

    @classmethod
    @Appender(_interval_shared_docs['from_tuples'] % _shared_docs_kwargs)
    def from_tuples(cls, data, closed='right', copy=False, dtype=None):
        if len(data):
            left, right = [], []
        else:
            # ensure that empty data keeps input dtype
            left = right = data

        for d in data:
            if isna(d):
                lhs = rhs = np.nan
            else:
                name = cls.__name__
                try:
                    # need list of length 2 tuples, e.g. [(0, 1), (1, 2), ...]
                    lhs, rhs = d
                except ValueError:
                    msg = ('{name}.from_tuples requires tuples of '
                           'length 2, got {tpl}').format(name=name, tpl=d)
                    raise ValueError(msg)
                except TypeError:
                    msg = ('{name}.from_tuples received an invalid '
                           'item, {tpl}').format(name=name, tpl=d)
                    raise TypeError(msg)
            left.append(lhs)
            right.append(rhs)

        return cls.from_arrays(left, right, closed, copy=False,
                               dtype=dtype)

    def _validate(self):
        """Verify that the IntervalArray is valid.

        Checks that

        * closed is valid
        * left and right match lengths
        * left and right have the same missing values
        * left is always below right
        """
        if self.closed not in _VALID_CLOSED:
            raise ValueError("invalid option for 'closed': {closed}"
                             .format(closed=self.closed))
        if len(self.left) != len(self.right):
            raise ValueError('left and right must have the same length')
        left_mask = notna(self.left)
        right_mask = notna(self.right)
        if not (left_mask == right_mask).all():
            raise ValueError('missing values must be missing in the same '
                             'location both left and right sides')
        if not (self.left[left_mask] <= self.right[left_mask]).all():
            raise ValueError('left side of interval must be <= right side')

    # ---------
    # Interface
    # ---------
    def __iter__(self):
        return iter(np.asarray(self))

    def __len__(self):
        return len(self.left)

    def __getitem__(self, value):
        left = self.left[value]
        right = self.right[value]

        # scalar
        if not isinstance(left, Index):
            if isna(left):
                return self._fill_value
            return Interval(left, right, self.closed)

        return self._shallow_copy(left, right)

    def __setitem__(self, key, value):
        # na value: need special casing to set directly on numpy arrays
        needs_float_conversion = False
        if is_scalar(value) and isna(value):
            if is_integer_dtype(self.dtype.subtype):
                # can't set NaN on a numpy integer array
                needs_float_conversion = True
            elif is_datetime64_any_dtype(self.dtype.subtype):
                # need proper NaT to set directly on the numpy array
                value = np.datetime64('NaT')
            elif is_timedelta64_dtype(self.dtype.subtype):
                # need proper NaT to set directly on the numpy array
                value = np.timedelta64('NaT')
            value_left, value_right = value, value

        # scalar interval
        elif is_interval_dtype(value) or isinstance(value, ABCInterval):
            self._check_closed_matches(value, name="value")
            value_left, value_right = value.left, value.right

        else:
            # list-like of intervals
            try:
                array = IntervalArray(value)
                value_left, value_right = array.left, array.right
            except TypeError:
                # wrong type: not interval or NA
                msg = "'value' should be an interval type, got {} instead."
                raise TypeError(msg.format(type(value)))

        # Need to ensure that left and right are updated atomically, so we're
        # forced to copy, update the copy, and swap in the new values.
        left = self.left.copy(deep=True)
        if needs_float_conversion:
            left = left.astype('float')
        left.values[key] = value_left
        self._left = left

        right = self.right.copy(deep=True)
        if needs_float_conversion:
            right = right.astype('float')
        right.values[key] = value_right
        self._right = right

    def fillna(self, value=None, method=None, limit=None):
        """
        Fill NA/NaN values using the specified method.

        Parameters
        ----------
        value : scalar, dict, Series
            If a scalar value is passed it is used to fill all missing values.
            Alternatively, a Series or dict can be used to fill in different
            values for each index. The value should not be a list. The
            value(s) passed should be either Interval objects or NA/NaN.
        method : {'backfill', 'bfill', 'pad', 'ffill', None}, default None
            (Not implemented yet for IntervalArray)
            Method to use for filling holes in reindexed Series
        limit : int, default None
            (Not implemented yet for IntervalArray)
            If method is specified, this is the maximum number of consecutive
            NaN values to forward/backward fill. In other words, if there is
            a gap with more than this number of consecutive NaNs, it will only
            be partially filled. If method is not specified, this is the
            maximum number of entries along the entire axis where NaNs will be
            filled.

        Returns
        -------
        filled : IntervalArray with NA/NaN filled
        """
        if method is not None:
            raise TypeError('Filling by method is not supported for '
                            'IntervalArray.')
        if limit is not None:
            raise TypeError('limit is not supported for IntervalArray.')

        if not isinstance(value, ABCInterval):
            msg = ("'IntervalArray.fillna' only supports filling with a "
                   "scalar 'pandas.Interval'. Got a '{}' instead."
                   .format(type(value).__name__))
            raise TypeError(msg)

        value = getattr(value, '_values', value)
        self._check_closed_matches(value, name="value")

        left = self.left.fillna(value=value.left)
        right = self.right.fillna(value=value.right)
        return self._shallow_copy(left, right)

    @property
    def dtype(self):
        return IntervalDtype(self.left.dtype)

    def astype(self, dtype, copy=True):
        """
        Cast to an ExtensionArray or NumPy array with dtype 'dtype'.

        Parameters
        ----------
        dtype : str or dtype
            Typecode or data-type to which the array is cast.

        copy : bool, default True
            Whether to copy the data, even if not necessary. If False,
            a copy is made only if the old dtype does not match the
            new dtype.

        Returns
        -------
        array : ExtensionArray or ndarray
            ExtensionArray or NumPy ndarray with 'dtype' for its dtype.
        """
        dtype = pandas_dtype(dtype)
        if is_interval_dtype(dtype):
            if dtype == self.dtype:
                return self.copy() if copy else self

            # need to cast to different subtype
            try:
                new_left = self.left.astype(dtype.subtype)
                new_right = self.right.astype(dtype.subtype)
            except TypeError:
                msg = ('Cannot convert {dtype} to {new_dtype}; subtypes are '
                       'incompatible')
                raise TypeError(msg.format(dtype=self.dtype, new_dtype=dtype))
            return self._shallow_copy(new_left, new_right)
        elif is_categorical_dtype(dtype):
            return Categorical(np.asarray(self))
        # TODO: This try/except will be repeated.
        try:
            return np.asarray(self).astype(dtype, copy=copy)
        except (TypeError, ValueError):
            msg = 'Cannot cast {name} to dtype {dtype}'
            raise TypeError(msg.format(name=type(self).__name__, dtype=dtype))

    @classmethod
    def _concat_same_type(cls, to_concat):
        """
        Concatenate multiple IntervalArray

        Parameters
        ----------
        to_concat : sequence of IntervalArray

        Returns
        -------
        IntervalArray
        """
        closed = {interval.closed for interval in to_concat}
        if len(closed) != 1:
            raise ValueError("Intervals must all be closed on the same side.")
        closed = closed.pop()

        left = np.concatenate([interval.left for interval in to_concat])
        right = np.concatenate([interval.right for interval in to_concat])
        return cls._simple_new(left, right, closed=closed, copy=False)

    def _shallow_copy(self, left=None, right=None, closed=None):
        """
        Return a new IntervalArray with the replacement attributes

        Parameters
        ----------
        left : array-like
            Values to be used for the left-side of the the intervals.
            If None, the existing left and right values will be used.

        right : array-like
            Values to be used for the right-side of the the intervals.
            If None and left is IntervalArray-like, the left and right
            of the IntervalArray-like will be used.

        closed : {'left', 'right', 'both', 'neither'}, optional
            Whether the intervals are closed on the left-side, right-side, both
            or neither.  If None, the existing closed will be used.
        """
        if left is None:

            # no values passed
            left, right = self.left, self.right

        elif right is None:

            # only single value passed, could be an IntervalArray
            # or array of Intervals
            if not isinstance(left, (type(self), ABCIntervalIndex)):
                left = type(self)(left)

            left, right = left.left, left.right
        else:

            # both left and right are values
            pass

        closed = closed or self.closed
        return self._simple_new(
            left, right, closed=closed, verify_integrity=False)

    def copy(self, deep=False):
        """
        Return a copy of the array.

        Parameters
        ----------
        deep : bool, default False
            Also copy the underlying data backing this array.

        Returns
        -------
        IntervalArray
        """
        left = self.left.copy(deep=True) if deep else self.left
        right = self.right.copy(deep=True) if deep else self.right
        closed = self.closed
        # TODO: Could skip verify_integrity here.
        return type(self).from_arrays(left, right, closed=closed)

    def _formatting_values(self):
        return np.asarray(self)

    def isna(self):
        return isna(self.left)

    @property
    def nbytes(self):
        return self.left.nbytes + self.right.nbytes

    @property
    def size(self):
        # Avoid materializing self.values
        return self.left.size

    @property
    def shape(self):
        return self.left.shape

    def take(self, indices, allow_fill=False, fill_value=None, axis=None,
             **kwargs):
        """
        Take elements from the IntervalArray.

        Parameters
        ----------
        indices : sequence of integers
            Indices to be taken.

        allow_fill : bool, default False
            How to handle negative values in `indices`.

            * False: negative values in `indices` indicate positional indices
              from the right (the default). This is similar to
              :func:`numpy.take`.

            * True: negative values in `indices` indicate
              missing values. These values are set to `fill_value`. Any other
              other negative values raise a ``ValueError``.

        fill_value : Interval or NA, optional
            Fill value to use for NA-indices when `allow_fill` is True.
            This may be ``None``, in which case the default NA value for
            the type, ``self.dtype.na_value``, is used.

            For many ExtensionArrays, there will be two representations of
            `fill_value`: a user-facing "boxed" scalar, and a low-level
            physical NA value. `fill_value` should be the user-facing version,
            and the implementation should handle translating that to the
            physical version for processing the take if necessary.

        axis : any, default None
            Present for compat with IntervalIndex; does nothing.

        Returns
        -------
        IntervalArray

        Raises
        ------
        IndexError
            When the indices are out of bounds for the array.
        ValueError
            When `indices` contains negative values other than ``-1``
            and `allow_fill` is True.
        """
        from pandas.core.algorithms import take

        nv.validate_take(tuple(), kwargs)

        fill_left = fill_right = fill_value
        if allow_fill:
            if fill_value is None:
                fill_left = fill_right = self.left._na_value
            elif is_interval(fill_value):
                self._check_closed_matches(fill_value, name='fill_value')
                fill_left, fill_right = fill_value.left, fill_value.right
            elif not is_scalar(fill_value) and notna(fill_value):
                msg = ("'IntervalArray.fillna' only supports filling with a "
                       "'scalar pandas.Interval or NA'. Got a '{}' instead."
                       .format(type(fill_value).__name__))
                raise ValueError(msg)

        left_take = take(self.left, indices,
                         allow_fill=allow_fill, fill_value=fill_left)
        right_take = take(self.right, indices,
                          allow_fill=allow_fill, fill_value=fill_right)

        return self._shallow_copy(left_take, right_take)

    def value_counts(self, dropna=True):
        """
        Returns a Series containing counts of each interval.

        Parameters
        ----------
        dropna : boolean, default True
            Don't include counts of NaN.

        Returns
        -------
        counts : Series

        See Also
        --------
        Series.value_counts
        """
        # TODO: implement this is a non-naive way!
        from pandas.core.algorithms import value_counts
        return value_counts(np.asarray(self), dropna=dropna)

    # Formatting

    def _format_data(self):

        # TODO: integrate with categorical and make generic
        # name argument is unused here; just for compat with base / categorical
        n = len(self)
        max_seq_items = min((get_option(
            'display.max_seq_items') or n) // 10, 10)

        formatter = str

        if n == 0:
            summary = '[]'
        elif n == 1:
            first = formatter(self[0])
            summary = '[{first}]'.format(first=first)
        elif n == 2:
            first = formatter(self[0])
            last = formatter(self[-1])
            summary = '[{first}, {last}]'.format(first=first, last=last)
        else:

            if n > max_seq_items:
                n = min(max_seq_items // 2, 10)
                head = [formatter(x) for x in self[:n]]
                tail = [formatter(x) for x in self[-n:]]
                summary = '[{head} ... {tail}]'.format(
                    head=', '.join(head), tail=', '.join(tail))
            else:
                tail = [formatter(x) for x in self]
                summary = '[{tail}]'.format(tail=', '.join(tail))

        return summary

    def __repr__(self):
        tpl = textwrap.dedent("""\
        {cls}({data},
        {lead}closed='{closed}',
        {lead}dtype='{dtype}')""")
        return tpl.format(cls=self.__class__.__name__,
                          data=self._format_data(),
                          lead=' ' * len(self.__class__.__name__) + ' ',
                          closed=self.closed, dtype=self.dtype)

    def _format_space(self):
        space = ' ' * (len(self.__class__.__name__) + 1)
        return "\n{space}".format(space=space)

    @property
    def left(self):
        """
        Return the left endpoints of each Interval in the IntervalArray as
        an Index
        """
        return self._left

    @property
    def right(self):
        """
        Return the right endpoints of each Interval in the IntervalArray as
        an Index
        """
        return self._right

    @property
    def closed(self):
        """
        Whether the intervals are closed on the left-side, right-side, both or
        neither
        """
        return self._closed

    _interval_shared_docs['set_closed'] = """
        Return an %(klass)s identical to the current one, but closed on the
        specified side

        .. versionadded:: 0.24.0

        Parameters
        ----------
        closed : {'left', 'right', 'both', 'neither'}
            Whether the intervals are closed on the left-side, right-side, both
            or neither.

        Returns
        -------
        new_index : %(klass)s

        Examples
        --------
        >>>  index = pd.interval_range(0, 3)
        >>>  index
        %(klass)s([(0, 1], (1, 2], (2, 3]]
              closed='right',
              dtype='interval[int64]')
        >>>  index.set_closed('both')
        %(klass)s([[0, 1], [1, 2], [2, 3]]
              closed='both',
              dtype='interval[int64]')
        """

    @Appender(_interval_shared_docs['set_closed'] % _shared_docs_kwargs)
    def set_closed(self, closed):
        if closed not in _VALID_CLOSED:
            msg = "invalid option for 'closed': {closed}"
            raise ValueError(msg.format(closed=closed))

        return self._shallow_copy(closed=closed)

    @property
    def length(self):
        """
        Return an Index with entries denoting the length of each Interval in
        the IntervalArray
        """
        try:
            return self.right - self.left
        except TypeError:
            # length not defined for some types, e.g. string
            msg = ('IntervalArray contains Intervals without defined length, '
                   'e.g. Intervals with string endpoints')
            raise TypeError(msg)

    @property
    def mid(self):
        """
        Return the midpoint of each Interval in the IntervalArray as an Index
        """
        try:
            return 0.5 * (self.left + self.right)
        except TypeError:
            # datetime safe version
            return self.left + 0.5 * self.length

    @property
    def is_non_overlapping_monotonic(self):
        """
        Return True if the IntervalArray is non-overlapping (no Intervals share
        points) and is either monotonic increasing or monotonic decreasing,
        else False
        """
        # must be increasing  (e.g., [0, 1), [1, 2), [2, 3), ... )
        # or decreasing (e.g., [-1, 0), [-2, -1), [-3, -2), ...)
        # we already require left <= right

        # strict inequality for closed == 'both'; equality implies overlapping
        # at a point when both sides of intervals are included
        if self.closed == 'both':
            return bool((self.right[:-1] < self.left[1:]).all() or
                        (self.left[:-1] > self.right[1:]).all())

        # non-strict inequality when closed != 'both'; at least one side is
        # not included in the intervals, so equality does not imply overlapping
        return bool((self.right[:-1] <= self.left[1:]).all() or
                    (self.left[:-1] >= self.right[1:]).all())

    # Conversion
    def __array__(self, dtype=None):
        """
        Return the IntervalArray's data as a numpy array of Interval
        objects (with dtype='object')
        """
        left = self.left
        right = self.right
        mask = self.isna()
        closed = self._closed

        result = np.empty(len(left), dtype=object)
        for i in range(len(left)):
            if mask[i]:
                result[i] = np.nan
            else:
                result[i] = Interval(left[i], right[i], closed)
        return result

    _interval_shared_docs['to_tuples'] = """\
        Return an %(return_type)s of tuples of the form (left, right)

        Parameters
        ----------
        na_tuple : boolean, default True
            Returns NA as a tuple if True, ``(nan, nan)``, or just as the NA
            value itself if False, ``nan``.

            ..versionadded:: 0.23.0

        Returns
        -------
        tuples: %(return_type)s
        %(examples)s\
    """

    @Appender(_interval_shared_docs['to_tuples'] % dict(
        return_type='ndarray',
        examples='',
    ))
    def to_tuples(self, na_tuple=True):
        tuples = com.asarray_tuplesafe(zip(self.left, self.right))
        if not na_tuple:
            # GH 18756
            tuples = np.where(~self.isna(), tuples, np.nan)
        return tuples

    def repeat(self, repeats, **kwargs):
        """
        Repeat elements of an IntervalArray.

        Returns a new IntervalArray where each element of the current
        IntervalArray is repeated consecutively a given number of times.

        Parameters
        ----------
        repeats : int
            The number of repetitions for each element.

        **kwargs
            Additional keywords have no effect but might be accepted for
            compatibility with numpy.

        Returns
        -------
        IntervalArray
            Newly created IntervalArray with repeated elements.

        See Also
        --------
        Index.repeat : Equivalent function for Index
        Series.repeat : Equivalent function for Series
        numpy.repeat : Underlying implementation
        """
        left_repeat = self.left.repeat(repeats, **kwargs)
        right_repeat = self.right.repeat(repeats, **kwargs)
        return self._shallow_copy(left=left_repeat, right=right_repeat)


def maybe_convert_platform_interval(values):
    """
    Try to do platform conversion, with special casing for IntervalArray.
    Wrapper around maybe_convert_platform that alters the default return
    dtype in certain cases to be compatible with IntervalArray.  For example,
    empty lists return with integer dtype instead of object dtype, which is
    prohibited for IntervalArray.

    Parameters
    ----------
    values : array-like

    Returns
    -------
    array
    """
    if isinstance(values, (list, tuple)) and len(values) == 0:
        # GH 19016
        # empty lists/tuples get object dtype by default, but this is not
        # prohibited for IntervalArray, so coerce to integer instead
        return np.array([], dtype=np.int64)
    elif is_categorical_dtype(values):
        values = np.asarray(values)

    return maybe_convert_platform(values)
