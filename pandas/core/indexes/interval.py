""" define the IntervalIndex """
from operator import le, lt
import textwrap
from typing import Any, Optional, Tuple, Union

import numpy as np

from pandas._config import get_option

from pandas._libs import Timedelta, Timestamp, lib
from pandas._libs.interval import Interval, IntervalMixin, IntervalTree
from pandas._typing import AnyArrayLike
from pandas.util._decorators import Appender, Substitution, cache_readonly
from pandas.util._exceptions import rewrite_exception

from pandas.core.dtypes.cast import (
    find_common_type,
    infer_dtype_from_scalar,
    maybe_downcast_to_dtype,
)
from pandas.core.dtypes.common import (
    ensure_platform_int,
    is_categorical,
    is_datetime64tz_dtype,
    is_datetime_or_timedelta_dtype,
    is_dtype_equal,
    is_float,
    is_float_dtype,
    is_integer,
    is_integer_dtype,
    is_interval_dtype,
    is_list_like,
    is_number,
    is_object_dtype,
    is_scalar,
)
from pandas.core.dtypes.missing import isna

from pandas.core.algorithms import take_1d
from pandas.core.arrays.interval import IntervalArray, _interval_shared_docs
import pandas.core.common as com
import pandas.core.indexes.base as ibase
from pandas.core.indexes.base import (
    Index,
    InvalidIndexError,
    _index_shared_docs,
    default_pprint,
    ensure_index,
    maybe_extract_name,
)
from pandas.core.indexes.datetimes import DatetimeIndex, date_range
from pandas.core.indexes.extension import ExtensionIndex, inherit_names
from pandas.core.indexes.multi import MultiIndex
from pandas.core.indexes.timedeltas import TimedeltaIndex, timedelta_range
from pandas.core.ops import get_op_result_name

from pandas.tseries.frequencies import to_offset
from pandas.tseries.offsets import DateOffset

_VALID_CLOSED = {"left", "right", "both", "neither"}
_index_doc_kwargs = dict(ibase._index_doc_kwargs)

_index_doc_kwargs.update(
    dict(
        klass="IntervalIndex",
        qualname="IntervalIndex",
        target_klass="IntervalIndex or list of Intervals",
        name=textwrap.dedent(
            """\
         name : object, optional
              Name to be stored in the index.
         """
        ),
    )
)


def _get_next_label(label):
    dtype = getattr(label, "dtype", type(label))
    if isinstance(label, (Timestamp, Timedelta)):
        dtype = "datetime64"
    if is_datetime_or_timedelta_dtype(dtype) or is_datetime64tz_dtype(dtype):
        return label + np.timedelta64(1, "ns")
    elif is_integer_dtype(dtype):
        return label + 1
    elif is_float_dtype(dtype):
        return np.nextafter(label, np.infty)
    else:
        raise TypeError(f"cannot determine next label for type {repr(type(label))}")


def _get_prev_label(label):
    dtype = getattr(label, "dtype", type(label))
    if isinstance(label, (Timestamp, Timedelta)):
        dtype = "datetime64"
    if is_datetime_or_timedelta_dtype(dtype) or is_datetime64tz_dtype(dtype):
        return label - np.timedelta64(1, "ns")
    elif is_integer_dtype(dtype):
        return label - 1
    elif is_float_dtype(dtype):
        return np.nextafter(label, -np.infty)
    else:
        raise TypeError(f"cannot determine next label for type {repr(type(label))}")


def _new_IntervalIndex(cls, d):
    """
    This is called upon unpickling, rather than the default which doesn't have
    arguments and breaks __new__.
    """
    return cls.from_arrays(**d)


class SetopCheck:
    """
    This is called to decorate the set operations of IntervalIndex
    to perform the type check in advance.
    """

    def __init__(self, op_name):
        self.op_name = op_name

    def __call__(self, setop):
        def func(intvidx_self, other, sort=False):
            intvidx_self._assert_can_do_setop(other)
            other = ensure_index(other)

            if not isinstance(other, IntervalIndex):
                result = getattr(intvidx_self.astype(object), self.op_name)(other)
                if self.op_name in ("difference",):
                    result = result.astype(intvidx_self.dtype)
                return result
            elif intvidx_self.closed != other.closed:
                raise ValueError(
                    "can only do set operations between two IntervalIndex "
                    "objects that are closed on the same side"
                )

            # GH 19016: ensure set op will not return a prohibited dtype
            subtypes = [intvidx_self.dtype.subtype, other.dtype.subtype]
            common_subtype = find_common_type(subtypes)
            if is_object_dtype(common_subtype):
                raise TypeError(
                    f"can only do {self.op_name} between two IntervalIndex "
                    "objects that have compatible dtypes"
                )

            return setop(intvidx_self, other, sort)

        return func


@Appender(
    _interval_shared_docs["class"]
    % dict(
        klass="IntervalIndex",
        summary="Immutable index of intervals that are closed on the same side.",
        name=_index_doc_kwargs["name"],
        versionadded="0.20.0",
        extra_attributes="is_overlapping\nvalues\n",
        extra_methods="",
        examples=textwrap.dedent(
            """\
    Examples
    --------
    A new ``IntervalIndex`` is typically constructed using
    :func:`interval_range`:

    >>> pd.interval_range(start=0, end=5)
    IntervalIndex([(0, 1], (1, 2], (2, 3], (3, 4], (4, 5]],
                  closed='right',
                  dtype='interval[int64]')

    It may also be constructed using one of the constructor
    methods: :meth:`IntervalIndex.from_arrays`,
    :meth:`IntervalIndex.from_breaks`, and :meth:`IntervalIndex.from_tuples`.

    See further examples in the doc strings of ``interval_range`` and the
    mentioned constructor methods.
    """
        ),
    )
)
@inherit_names(["set_closed", "to_tuples"], IntervalArray, wrap=True)
@inherit_names(
    [
        "__len__",
        "__array__",
        "overlaps",
        "contains",
        "size",
        "dtype",
        "left",
        "right",
        "length",
    ],
    IntervalArray,
)
@inherit_names(
    ["is_non_overlapping_monotonic", "mid", "_ndarray_values", "closed"],
    IntervalArray,
    cache=True,
)
class IntervalIndex(IntervalMixin, ExtensionIndex):
    _typ = "intervalindex"
    _comparables = ["name"]
    _attributes = ["name", "closed"]

    # we would like our indexing holder to defer to us
    _defer_to_indexing = True

    # Immutable, so we are able to cache computations like isna in '_mask'
    _mask = None

    _data: IntervalArray
    # --------------------------------------------------------------------
    # Constructors

    def __new__(
        cls,
        data,
        closed=None,
        dtype=None,
        copy: bool = False,
        name=None,
        verify_integrity: bool = True,
    ):

        name = maybe_extract_name(name, data, cls)

        with rewrite_exception("IntervalArray", cls.__name__):
            array = IntervalArray(
                data,
                closed=closed,
                copy=copy,
                dtype=dtype,
                verify_integrity=verify_integrity,
            )

        return cls._simple_new(array, name)

    @classmethod
    def _simple_new(cls, array, name, closed=None):
        """
        Construct from an IntervalArray

        Parameters
        ----------
        array : IntervalArray
        name : str
            Attached as result.name
        closed : Any
            Ignored.
        """
        assert isinstance(array, IntervalArray), type(array)

        result = IntervalMixin.__new__(cls)
        result._data = array
        result.name = name
        result._no_setting_name = False
        result._reset_identity()
        return result

    @classmethod
    @Appender(
        _interval_shared_docs["from_breaks"]
        % dict(
            klass="IntervalIndex",
            examples=textwrap.dedent(
                """\
        Examples
        --------
        >>> pd.IntervalIndex.from_breaks([0, 1, 2, 3])
        IntervalIndex([(0, 1], (1, 2], (2, 3]],
                      closed='right',
                      dtype='interval[int64]')
        """
            ),
        )
    )
    def from_breaks(
        cls, breaks, closed: str = "right", name=None, copy: bool = False, dtype=None
    ):
        with rewrite_exception("IntervalArray", cls.__name__):
            array = IntervalArray.from_breaks(
                breaks, closed=closed, copy=copy, dtype=dtype
            )
        return cls._simple_new(array, name=name)

    @classmethod
    @Appender(
        _interval_shared_docs["from_arrays"]
        % dict(
            klass="IntervalIndex",
            examples=textwrap.dedent(
                """\
        Examples
        --------
        >>> pd.IntervalIndex.from_arrays([0, 1, 2], [1, 2, 3])
        IntervalIndex([(0, 1], (1, 2], (2, 3]],
                      closed='right',
                      dtype='interval[int64]')
        """
            ),
        )
    )
    def from_arrays(
        cls,
        left,
        right,
        closed: str = "right",
        name=None,
        copy: bool = False,
        dtype=None,
    ):
        with rewrite_exception("IntervalArray", cls.__name__):
            array = IntervalArray.from_arrays(
                left, right, closed, copy=copy, dtype=dtype
            )
        return cls._simple_new(array, name=name)

    @classmethod
    @Appender(
        _interval_shared_docs["from_tuples"]
        % dict(
            klass="IntervalIndex",
            examples=textwrap.dedent(
                """\
        Examples
        --------
        >>> pd.IntervalIndex.from_tuples([(0, 1), (1, 2)])
        IntervalIndex([(0, 1], (1, 2]],
                       closed='right',
                       dtype='interval[int64]')
        """
            ),
        )
    )
    def from_tuples(
        cls, data, closed: str = "right", name=None, copy: bool = False, dtype=None
    ):
        with rewrite_exception("IntervalArray", cls.__name__):
            arr = IntervalArray.from_tuples(data, closed=closed, copy=copy, dtype=dtype)
        return cls._simple_new(arr, name=name)

    # --------------------------------------------------------------------

    @Appender(Index._shallow_copy.__doc__)
    def _shallow_copy(self, left=None, right=None, **kwargs):
        result = self._data._shallow_copy(left=left, right=right)
        attributes = self._get_attributes_dict()
        attributes.update(kwargs)
        return self._simple_new(result, **attributes)

    @cache_readonly
    def _isnan(self):
        """
        Return a mask indicating if each value is NA.
        """
        if self._mask is None:
            self._mask = isna(self.left)
        return self._mask

    @cache_readonly
    def _engine(self):
        left = self._maybe_convert_i8(self.left)
        right = self._maybe_convert_i8(self.right)
        return IntervalTree(left, right, closed=self.closed)

    def __contains__(self, key: Any) -> bool:
        """
        return a boolean if this key is IN the index
        We *only* accept an Interval

        Parameters
        ----------
        key : Interval

        Returns
        -------
        bool
        """
        hash(key)
        if not isinstance(key, Interval):
            return False

        try:
            self.get_loc(key)
            return True
        except KeyError:
            return False

    @cache_readonly
    def _multiindex(self) -> MultiIndex:
        return MultiIndex.from_arrays([self.left, self.right], names=["left", "right"])

    @cache_readonly
    def values(self) -> IntervalArray:
        """
        Return the IntervalIndex's data as an IntervalArray.
        """
        return self._data

    @property
    def _has_complex_internals(self) -> bool:
        # used to avoid libreduction code paths, which raise or require conversion
        return True

    def __array_wrap__(self, result, context=None):
        # we don't want the superclass implementation
        return result

    def __reduce__(self):
        d = dict(left=self.left, right=self.right)
        d.update(self._get_attributes_dict())
        return _new_IntervalIndex, (type(self), d), None

    @Appender(Index.astype.__doc__)
    def astype(self, dtype, copy=True):
        with rewrite_exception("IntervalArray", type(self).__name__):
            new_values = self.values.astype(dtype, copy=copy)
        if is_interval_dtype(new_values):
            return self._shallow_copy(new_values.left, new_values.right)
        return Index.astype(self, dtype, copy=copy)

    @property
    def inferred_type(self) -> str:
        """Return a string of the type inferred from the values"""
        return "interval"

    @Appender(Index.memory_usage.__doc__)
    def memory_usage(self, deep: bool = False) -> int:
        # we don't use an explicit engine
        # so return the bytes here
        return self.left.memory_usage(deep=deep) + self.right.memory_usage(deep=deep)

    # IntervalTree doesn't have a is_monotonic_decreasing, so have to override
    #  the Index implemenation
    @cache_readonly
    def is_monotonic_decreasing(self) -> bool:
        """
        Return True if the IntervalIndex is monotonic decreasing (only equal or
        decreasing values), else False
        """
        return self[::-1].is_monotonic_increasing

    @cache_readonly
    def is_unique(self):
        """
        Return True if the IntervalIndex contains unique elements, else False.
        """
        left = self.left
        right = self.right

        if self.isna().sum() > 1:
            return False

        if left.is_unique or right.is_unique:
            return True

        seen_pairs = set()
        check_idx = np.where(left.duplicated(keep=False))[0]
        for idx in check_idx:
            pair = (left[idx], right[idx])
            if pair in seen_pairs:
                return False
            seen_pairs.add(pair)

        return True

    @property
    def is_overlapping(self) -> bool:
        """
        Return True if the IntervalIndex has overlapping intervals, else False.

        Two intervals overlap if they share a common point, including closed
        endpoints. Intervals that only have an open endpoint in common do not
        overlap.

        .. versionadded:: 0.24.0

        Returns
        -------
        bool
            Boolean indicating if the IntervalIndex has overlapping intervals.

        See Also
        --------
        Interval.overlaps : Check whether two Interval objects overlap.
        IntervalIndex.overlaps : Check an IntervalIndex elementwise for
            overlaps.

        Examples
        --------
        >>> index = pd.IntervalIndex.from_tuples([(0, 2), (1, 3), (4, 5)])
        >>> index
        IntervalIndex([(0, 2], (1, 3], (4, 5]],
              closed='right',
              dtype='interval[int64]')
        >>> index.is_overlapping
        True

        Intervals that share closed endpoints overlap:

        >>> index = pd.interval_range(0, 3, closed='both')
        >>> index
        IntervalIndex([[0, 1], [1, 2], [2, 3]],
              closed='both',
              dtype='interval[int64]')
        >>> index.is_overlapping
        True

        Intervals that only have an open endpoint in common do not overlap:

        >>> index = pd.interval_range(0, 3, closed='left')
        >>> index
        IntervalIndex([[0, 1), [1, 2), [2, 3)],
              closed='left',
              dtype='interval[int64]')
        >>> index.is_overlapping
        False
        """
        # GH 23309
        return self._engine.is_overlapping

    def holds_integer(self):
        return self.dtype.subtype.kind not in ["m", "M"]
        # TODO: There must already exist something for this?

    @Appender(Index._convert_scalar_indexer.__doc__)
    def _convert_scalar_indexer(self, key, kind=None):
        if kind == "iloc":
            return super()._convert_scalar_indexer(key, kind=kind)
        return key

    def _maybe_cast_slice_bound(self, label, side, kind):
        return getattr(self, side)._maybe_cast_slice_bound(label, side, kind)

    @Appender(Index._convert_list_indexer.__doc__)
    def _convert_list_indexer(self, keyarr, kind=None):
        """
        we are passed a list-like indexer. Return the
        indexer for matching intervals.
        """
        locs = self.get_indexer_for(keyarr)

        # we have missing values
        if (locs == -1).any():
            raise KeyError

        return locs

    def _can_reindex(self, indexer: np.ndarray) -> None:
        """
        Check if we are allowing reindexing with this particular indexer.

        Parameters
        ----------
        indexer : an integer indexer

        Raises
        ------
        ValueError if its a duplicate axis
        """

        # trying to reindex on an axis with duplicates
        if self.is_overlapping and len(indexer):
            raise ValueError("cannot reindex from an overlapping axis")

    def _needs_i8_conversion(self, key) -> bool:
        """
        Check if a given key needs i8 conversion. Conversion is necessary for
        Timestamp, Timedelta, DatetimeIndex, and TimedeltaIndex keys. An
        Interval-like requires conversion if it's endpoints are one of the
        aforementioned types.

        Assumes that any list-like data has already been cast to an Index.

        Parameters
        ----------
        key : scalar or Index-like
            The key that should be checked for i8 conversion

        Returns
        -------
        bool
        """
        if is_interval_dtype(key) or isinstance(key, Interval):
            return self._needs_i8_conversion(key.left)

        i8_types = (Timestamp, Timedelta, DatetimeIndex, TimedeltaIndex)
        return isinstance(key, i8_types)

    def _maybe_convert_i8(self, key):
        """
        Maybe convert a given key to it's equivalent i8 value(s). Used as a
        preprocessing step prior to IntervalTree queries (self._engine), which
        expects numeric data.

        Parameters
        ----------
        key : scalar or list-like
            The key that should maybe be converted to i8.

        Returns
        -------
        scalar or list-like
            The original key if no conversion occurred, int if converted scalar,
            Int64Index if converted list-like.
        """
        original = key
        if is_list_like(key):
            key = ensure_index(key)

        if not self._needs_i8_conversion(key):
            return original

        scalar = is_scalar(key)
        if is_interval_dtype(key) or isinstance(key, Interval):
            # convert left/right and reconstruct
            left = self._maybe_convert_i8(key.left)
            right = self._maybe_convert_i8(key.right)
            constructor = Interval if scalar else IntervalIndex.from_arrays
            return constructor(left, right, closed=self.closed)

        if scalar:
            # Timestamp/Timedelta
            key_dtype, key_i8 = infer_dtype_from_scalar(key, pandas_dtype=True)
        else:
            # DatetimeIndex/TimedeltaIndex
            key_dtype, key_i8 = key.dtype, Index(key.asi8)
            if key.hasnans:
                # convert NaT from it's i8 value to np.nan so it's not viewed
                # as a valid value, maybe causing errors (e.g. is_overlapping)
                key_i8 = key_i8.where(~key._isnan)

        # ensure consistency with IntervalIndex subtype
        subtype = self.dtype.subtype

        if not is_dtype_equal(subtype, key_dtype):
            raise ValueError(
                f"Cannot index an IntervalIndex of subtype {subtype} with "
                f"values of dtype {key_dtype}"
            )

        return key_i8

    def _check_method(self, method):
        if method is None:
            return

        if method in ["bfill", "backfill", "pad", "ffill", "nearest"]:
            raise NotImplementedError(
                f"method {method} not yet implemented for IntervalIndex"
            )

        raise ValueError("Invalid fill method")

    def _searchsorted_monotonic(self, label, side, exclude_label=False):
        if not self.is_non_overlapping_monotonic:
            raise KeyError(
                "can only get slices from an IntervalIndex if bounds are "
                "non-overlapping and all monotonic increasing or decreasing"
            )

        if isinstance(label, IntervalMixin):
            raise NotImplementedError("Interval objects are not currently supported")

        # GH 20921: "not is_monotonic_increasing" for the second condition
        # instead of "is_monotonic_decreasing" to account for single element
        # indexes being both increasing and decreasing
        if (side == "left" and self.left.is_monotonic_increasing) or (
            side == "right" and not self.left.is_monotonic_increasing
        ):
            sub_idx = self.right
            if self.open_right or exclude_label:
                label = _get_next_label(label)
        else:
            sub_idx = self.left
            if self.open_left or exclude_label:
                label = _get_prev_label(label)

        return sub_idx._searchsorted_monotonic(label, side)

    def get_loc(
        self, key, method: Optional[str] = None, tolerance=None
    ) -> Union[int, slice, np.ndarray]:
        """
        Get integer location, slice or boolean mask for requested label.

        Parameters
        ----------
        key : label
        method : {None}, optional
            * default: matches where the label is within an interval only.

        Returns
        -------
        int if unique index, slice if monotonic index, else mask

        Examples
        --------
        >>> i1, i2 = pd.Interval(0, 1), pd.Interval(1, 2)
        >>> index = pd.IntervalIndex([i1, i2])
        >>> index.get_loc(1)
        0

        You can also supply a point inside an interval.

        >>> index.get_loc(1.5)
        1

        If a label is in several intervals, you get the locations of all the
        relevant intervals.

        >>> i3 = pd.Interval(0, 2)
        >>> overlapping_index = pd.IntervalIndex([i1, i2, i3])
        >>> overlapping_index.get_loc(0.5)
        array([ True, False,  True])

        Only exact matches will be returned if an interval is provided.

        >>> index.get_loc(pd.Interval(0, 1))
        0
        """
        self._check_method(method)

        if not is_scalar(key):
            raise InvalidIndexError(key)

        if isinstance(key, Interval):
            if self.closed != key.closed:
                raise KeyError(key)
            mask = (self.left == key.left) & (self.right == key.right)
        else:
            # assume scalar
            op_left = le if self.closed_left else lt
            op_right = le if self.closed_right else lt
            try:
                mask = op_left(self.left, key) & op_right(key, self.right)
            except TypeError:
                # scalar is not comparable to II subtype --> invalid label
                raise KeyError(key)

        matches = mask.sum()
        if matches == 0:
            raise KeyError(key)
        elif matches == 1:
            return mask.argmax()
        return lib.maybe_booleans_to_slice(mask.view("u1"))

    @Substitution(
        **dict(
            _index_doc_kwargs,
            **{
                "raises_section": textwrap.dedent(
                    """
        Raises
        ------
        NotImplementedError
            If any method argument other than the default of
            None is specified as these are not yet implemented.
        """
                )
            },
        )
    )
    @Appender(_index_shared_docs["get_indexer"])
    def get_indexer(
        self,
        target: AnyArrayLike,
        method: Optional[str] = None,
        limit: Optional[int] = None,
        tolerance: Optional[Any] = None,
    ) -> np.ndarray:

        self._check_method(method)

        if self.is_overlapping:
            raise InvalidIndexError(
                "cannot handle overlapping indices; "
                "use IntervalIndex.get_indexer_non_unique"
            )

        target_as_index = ensure_index(target)

        if isinstance(target_as_index, IntervalIndex):
            # equal indexes -> 1:1 positional match
            if self.equals(target_as_index):
                return np.arange(len(self), dtype="intp")

            # different closed or incompatible subtype -> no matches
            common_subtype = find_common_type(
                [self.dtype.subtype, target_as_index.dtype.subtype]
            )
            if self.closed != target_as_index.closed or is_object_dtype(common_subtype):
                return np.repeat(np.intp(-1), len(target_as_index))

            # non-overlapping -> at most one match per interval in target_as_index
            # want exact matches -> need both left/right to match, so defer to
            # left/right get_indexer, compare elementwise, equality -> match
            left_indexer = self.left.get_indexer(target_as_index.left)
            right_indexer = self.right.get_indexer(target_as_index.right)
            indexer = np.where(left_indexer == right_indexer, left_indexer, -1)
        elif is_categorical(target_as_index):
            # get an indexer for unique categories then propagate to codes via take_1d
            categories_indexer = self.get_indexer(target_as_index.categories)
            indexer = take_1d(categories_indexer, target_as_index.codes, fill_value=-1)
        elif not is_object_dtype(target_as_index):
            # homogeneous scalar index: use IntervalTree
            target_as_index = self._maybe_convert_i8(target_as_index)
            indexer = self._engine.get_indexer(target_as_index.values)
        else:
            # heterogeneous scalar index: defer elementwise to get_loc
            # (non-overlapping so get_loc guarantees scalar of KeyError)
            indexer = []
            for key in target_as_index:
                try:
                    loc = self.get_loc(key)
                except KeyError:
                    loc = -1
                except InvalidIndexError:
                    # i.e. non-scalar key
                    raise TypeError(key)
                indexer.append(loc)

        return ensure_platform_int(indexer)

    @Appender(_index_shared_docs["get_indexer_non_unique"] % _index_doc_kwargs)
    def get_indexer_non_unique(
        self, target: AnyArrayLike
    ) -> Tuple[np.ndarray, np.ndarray]:
        target_as_index = ensure_index(target)

        # check that target_as_index IntervalIndex is compatible
        if isinstance(target_as_index, IntervalIndex):
            common_subtype = find_common_type(
                [self.dtype.subtype, target_as_index.dtype.subtype]
            )
            if self.closed != target_as_index.closed or is_object_dtype(common_subtype):
                # different closed or incompatible subtype -> no matches
                return (
                    np.repeat(-1, len(target_as_index)),
                    np.arange(len(target_as_index)),
                )

        if is_object_dtype(target_as_index) or isinstance(
            target_as_index, IntervalIndex
        ):
            # target_as_index might contain intervals: defer elementwise to get_loc
            indexer, missing = [], []
            for i, key in enumerate(target_as_index):
                try:
                    locs = self.get_loc(key)
                    if isinstance(locs, slice):
                        locs = np.arange(locs.start, locs.stop, locs.step, dtype="intp")
                    locs = np.array(locs, ndmin=1)
                except KeyError:
                    missing.append(i)
                    locs = np.array([-1])
                indexer.append(locs)
            indexer = np.concatenate(indexer)
        else:
            target_as_index = self._maybe_convert_i8(target_as_index)
            indexer, missing = self._engine.get_indexer_non_unique(
                target_as_index.values
            )

        return ensure_platform_int(indexer), ensure_platform_int(missing)

    def get_indexer_for(self, target: AnyArrayLike, **kwargs) -> np.ndarray:
        """
        Guaranteed return of an indexer even when overlapping.

        This dispatches to get_indexer or get_indexer_non_unique
        as appropriate.

        Returns
        -------
        numpy.ndarray
            List of indices.
        """
        if self.is_overlapping:
            return self.get_indexer_non_unique(target)[0]
        return self.get_indexer(target, **kwargs)

    def _convert_slice_indexer(self, key: slice, kind=None):
        if not (key.step is None or key.step == 1):
            raise ValueError("cannot support not-default step in a slice")
        return super()._convert_slice_indexer(key, kind)

    @Appender(Index.where.__doc__)
    def where(self, cond, other=None):
        if other is None:
            other = self._na_value
        values = np.where(cond, self.values, other)
        return self._shallow_copy(values)

    def delete(self, loc):
        """
        Return a new IntervalIndex with passed location(-s) deleted

        Returns
        -------
        IntervalIndex
        """
        new_left = self.left.delete(loc)
        new_right = self.right.delete(loc)
        return self._shallow_copy(new_left, new_right)

    def insert(self, loc, item):
        """
        Return a new IntervalIndex inserting new item at location. Follows
        Python list.append semantics for negative values.  Only Interval
        objects and NA can be inserted into an IntervalIndex

        Parameters
        ----------
        loc : int
        item : object

        Returns
        -------
        IntervalIndex
        """
        if isinstance(item, Interval):
            if item.closed != self.closed:
                raise ValueError(
                    "inserted item must be closed on the same side as the index"
                )
            left_insert = item.left
            right_insert = item.right
        elif is_scalar(item) and isna(item):
            # GH 18295
            left_insert = right_insert = item
        else:
            raise ValueError(
                "can only insert Interval objects and NA into an IntervalIndex"
            )

        new_left = self.left.insert(loc, left_insert)
        new_right = self.right.insert(loc, right_insert)
        return self._shallow_copy(new_left, new_right)

    def _concat_same_dtype(self, to_concat, name):
        """
        assert that we all have the same .closed
        we allow a 0-len index here as well
        """
        if not len({i.closed for i in to_concat if len(i)}) == 1:
            raise ValueError(
                "can only append two IntervalIndex objects "
                "that are closed on the same side"
            )
        return super()._concat_same_dtype(to_concat, name)

    @Appender(_index_shared_docs["take"] % _index_doc_kwargs)
    def take(self, indices, axis=0, allow_fill=True, fill_value=None, **kwargs):
        result = self._data.take(
            indices, axis=axis, allow_fill=allow_fill, fill_value=fill_value, **kwargs
        )
        return self._shallow_copy(result)

    def __getitem__(self, value):
        result = self._data[value]
        if isinstance(result, IntervalArray):
            return self._shallow_copy(result)
        else:
            # scalar
            return result

    # --------------------------------------------------------------------
    # Rendering Methods
    # __repr__ associated methods are based on MultiIndex

    def _format_with_header(self, header, **kwargs):
        return header + list(self._format_native_types(**kwargs))

    def _format_native_types(self, na_rep="NaN", quoting=None, **kwargs):
        # GH 28210: use base method but with different default na_rep
        return super()._format_native_types(na_rep=na_rep, quoting=quoting, **kwargs)

    def _format_data(self, name=None):

        # TODO: integrate with categorical and make generic
        # name argument is unused here; just for compat with base / categorical
        n = len(self)
        max_seq_items = min((get_option("display.max_seq_items") or n) // 10, 10)

        formatter = str

        if n == 0:
            summary = "[]"
        elif n == 1:
            first = formatter(self[0])
            summary = f"[{first}]"
        elif n == 2:
            first = formatter(self[0])
            last = formatter(self[-1])
            summary = f"[{first}, {last}]"
        else:

            if n > max_seq_items:
                n = min(max_seq_items // 2, 10)
                head = [formatter(x) for x in self[:n]]
                tail = [formatter(x) for x in self[-n:]]
                head_joined = ", ".join(head)
                tail_joined = ", ".join(tail)
                summary = f"[{head_joined} ... {tail_joined}]"
            else:
                tail = [formatter(x) for x in self]
                joined = ", ".join(tail)
                summary = f"[{joined}]"

        return summary + "," + self._format_space()

    def _format_attrs(self):
        attrs = [("closed", repr(self.closed))]
        if self.name is not None:
            attrs.append(("name", default_pprint(self.name)))
        attrs.append(("dtype", f"'{self.dtype}'"))
        return attrs

    def _format_space(self) -> str:
        space = " " * (len(type(self).__name__) + 1)
        return f"\n{space}"

    # --------------------------------------------------------------------

    def argsort(self, *args, **kwargs) -> np.ndarray:
        return np.lexsort((self.right, self.left))

    def equals(self, other) -> bool:
        """
        Determines if two IntervalIndex objects contain the same elements.
        """
        if self.is_(other):
            return True

        # if we can coerce to an II
        # then we can compare
        if not isinstance(other, IntervalIndex):
            if not is_interval_dtype(other):
                return False
            other = Index(other)

        return (
            self.left.equals(other.left)
            and self.right.equals(other.right)
            and self.closed == other.closed
        )

    @Appender(Index.intersection.__doc__)
    @SetopCheck(op_name="intersection")
    def intersection(
        self, other: "IntervalIndex", sort: bool = False
    ) -> "IntervalIndex":
        if self.left.is_unique and self.right.is_unique:
            taken = self._intersection_unique(other)
        elif other.left.is_unique and other.right.is_unique and self.isna().sum() <= 1:
            # Swap other/self if other is unique and self does not have
            # multiple NaNs
            taken = other._intersection_unique(self)
        else:
            # duplicates
            taken = self._intersection_non_unique(other)

        if sort is None:
            taken = taken.sort_values()

        return taken

    def _intersection_unique(self, other: "IntervalIndex") -> "IntervalIndex":
        """
        Used when the IntervalIndex does not have any common endpoint,
        no mater left or right.
        Return the intersection with another IntervalIndex.

        Parameters
        ----------
        other : IntervalIndex

        Returns
        -------
        IntervalIndex
        """
        lindexer = self.left.get_indexer(other.left)
        rindexer = self.right.get_indexer(other.right)

        match = (lindexer == rindexer) & (lindexer != -1)
        indexer = lindexer.take(match.nonzero()[0])

        return self.take(indexer)

    def _intersection_non_unique(self, other: "IntervalIndex") -> "IntervalIndex":
        """
        Used when the IntervalIndex does have some common endpoints,
        on either sides.
        Return the intersection with another IntervalIndex.

        Parameters
        ----------
        other : IntervalIndex

        Returns
        -------
        IntervalIndex
        """
        mask = np.zeros(len(self), dtype=bool)

        if self.hasnans and other.hasnans:
            first_nan_loc = np.arange(len(self))[self.isna()][0]
            mask[first_nan_loc] = True

        other_tups = set(zip(other.left, other.right))
        for i, tup in enumerate(zip(self.left, self.right)):
            if tup in other_tups:
                mask[i] = True

        return self[mask]

    def _setop(op_name: str, sort=None):
        @SetopCheck(op_name=op_name)
        def func(self, other, sort=sort):
            result = getattr(self._multiindex, op_name)(other._multiindex, sort=sort)
            result_name = get_op_result_name(self, other)

            # GH 19101: ensure empty results have correct dtype
            if result.empty:
                result = result.values.astype(self.dtype.subtype)
            else:
                result = result.values

            return type(self).from_tuples(result, closed=self.closed, name=result_name)

        return func

    @property
    def is_all_dates(self) -> bool:
        """
        This is False even when left/right contain datetime-like objects,
        as the check is done on the Interval itself
        """
        return False

    union = _setop("union")
    difference = _setop("difference")
    symmetric_difference = _setop("symmetric_difference")

    # TODO: arithmetic operations

    # GH#30817 until IntervalArray implements inequalities, get them from Index
    def __lt__(self, other):
        return Index.__lt__(self, other)

    def __le__(self, other):
        return Index.__le__(self, other)

    def __gt__(self, other):
        return Index.__gt__(self, other)

    def __ge__(self, other):
        return Index.__ge__(self, other)


IntervalIndex._add_logical_methods_disabled()


def _is_valid_endpoint(endpoint) -> bool:
    """
    Helper for interval_range to check if start/end are valid types.
    """
    return any(
        [
            is_number(endpoint),
            isinstance(endpoint, Timestamp),
            isinstance(endpoint, Timedelta),
            endpoint is None,
        ]
    )


def _is_type_compatible(a, b) -> bool:
    """
    Helper for interval_range to check type compat of start/end/freq.
    """
    is_ts_compat = lambda x: isinstance(x, (Timestamp, DateOffset))
    is_td_compat = lambda x: isinstance(x, (Timedelta, DateOffset))
    return (
        (is_number(a) and is_number(b))
        or (is_ts_compat(a) and is_ts_compat(b))
        or (is_td_compat(a) and is_td_compat(b))
        or com.any_none(a, b)
    )


def interval_range(
    start=None, end=None, periods=None, freq=None, name=None, closed="right"
):
    """
    Return a fixed frequency IntervalIndex.

    Parameters
    ----------
    start : numeric or datetime-like, default None
        Left bound for generating intervals.
    end : numeric or datetime-like, default None
        Right bound for generating intervals.
    periods : int, default None
        Number of periods to generate.
    freq : numeric, str, or DateOffset, default None
        The length of each interval. Must be consistent with the type of start
        and end, e.g. 2 for numeric, or '5H' for datetime-like.  Default is 1
        for numeric and 'D' for datetime-like.
    name : str, default None
        Name of the resulting IntervalIndex.
    closed : {'left', 'right', 'both', 'neither'}, default 'right'
        Whether the intervals are closed on the left-side, right-side, both
        or neither.

    Returns
    -------
    IntervalIndex

    See Also
    --------
    IntervalIndex : An Index of intervals that are all closed on the same side.

    Notes
    -----
    Of the four parameters ``start``, ``end``, ``periods``, and ``freq``,
    exactly three must be specified. If ``freq`` is omitted, the resulting
    ``IntervalIndex`` will have ``periods`` linearly spaced elements between
    ``start`` and ``end``, inclusively.

    To learn more about datetime-like frequency strings, please see `this link
    <https://pandas.pydata.org/pandas-docs/stable/user_guide/timeseries.html#offset-aliases>`__.

    Examples
    --------
    Numeric ``start`` and  ``end`` is supported.

    >>> pd.interval_range(start=0, end=5)
    IntervalIndex([(0, 1], (1, 2], (2, 3], (3, 4], (4, 5]],
                  closed='right', dtype='interval[int64]')

    Additionally, datetime-like input is also supported.

    >>> pd.interval_range(start=pd.Timestamp('2017-01-01'),
    ...                   end=pd.Timestamp('2017-01-04'))
    IntervalIndex([(2017-01-01, 2017-01-02], (2017-01-02, 2017-01-03],
                   (2017-01-03, 2017-01-04]],
                  closed='right', dtype='interval[datetime64[ns]]')

    The ``freq`` parameter specifies the frequency between the left and right.
    endpoints of the individual intervals within the ``IntervalIndex``.  For
    numeric ``start`` and ``end``, the frequency must also be numeric.

    >>> pd.interval_range(start=0, periods=4, freq=1.5)
    IntervalIndex([(0.0, 1.5], (1.5, 3.0], (3.0, 4.5], (4.5, 6.0]],
                  closed='right', dtype='interval[float64]')

    Similarly, for datetime-like ``start`` and ``end``, the frequency must be
    convertible to a DateOffset.

    >>> pd.interval_range(start=pd.Timestamp('2017-01-01'),
    ...                   periods=3, freq='MS')
    IntervalIndex([(2017-01-01, 2017-02-01], (2017-02-01, 2017-03-01],
                   (2017-03-01, 2017-04-01]],
                  closed='right', dtype='interval[datetime64[ns]]')

    Specify ``start``, ``end``, and ``periods``; the frequency is generated
    automatically (linearly spaced).

    >>> pd.interval_range(start=0, end=6, periods=4)
    IntervalIndex([(0.0, 1.5], (1.5, 3.0], (3.0, 4.5], (4.5, 6.0]],
              closed='right',
              dtype='interval[float64]')

    The ``closed`` parameter specifies which endpoints of the individual
    intervals within the ``IntervalIndex`` are closed.

    >>> pd.interval_range(end=5, periods=4, closed='both')
    IntervalIndex([[1, 2], [2, 3], [3, 4], [4, 5]],
                  closed='both', dtype='interval[int64]')
    """
    start = com.maybe_box_datetimelike(start)
    end = com.maybe_box_datetimelike(end)
    endpoint = start if start is not None else end

    if freq is None and com.any_none(periods, start, end):
        freq = 1 if is_number(endpoint) else "D"

    if com.count_not_none(start, end, periods, freq) != 3:
        raise ValueError(
            "Of the four parameters: start, end, periods, and "
            "freq, exactly three must be specified"
        )

    if not _is_valid_endpoint(start):
        raise ValueError(f"start must be numeric or datetime-like, got {start}")
    elif not _is_valid_endpoint(end):
        raise ValueError(f"end must be numeric or datetime-like, got {end}")

    if is_float(periods):
        periods = int(periods)
    elif not is_integer(periods) and periods is not None:
        raise TypeError(f"periods must be a number, got {periods}")

    if freq is not None and not is_number(freq):
        try:
            freq = to_offset(freq)
        except ValueError:
            raise ValueError(
                f"freq must be numeric or convertible to DateOffset, got {freq}"
            )

    # verify type compatibility
    if not all(
        [
            _is_type_compatible(start, end),
            _is_type_compatible(start, freq),
            _is_type_compatible(end, freq),
        ]
    ):
        raise TypeError("start, end, freq need to be type compatible")

    # +1 to convert interval count to breaks count (n breaks = n-1 intervals)
    if periods is not None:
        periods += 1

    if is_number(endpoint):
        # force consistency between start/end/freq (lower end if freq skips it)
        if com.all_not_none(start, end, freq):
            end -= (end - start) % freq

        # compute the period/start/end if unspecified (at most one)
        if periods is None:
            periods = int((end - start) // freq) + 1
        elif start is None:
            start = end - (periods - 1) * freq
        elif end is None:
            end = start + (periods - 1) * freq

        breaks = np.linspace(start, end, periods)
        if all(is_integer(x) for x in com.not_none(start, end, freq)):
            # np.linspace always produces float output
            breaks = maybe_downcast_to_dtype(breaks, "int64")
    else:
        # delegate to the appropriate range function
        if isinstance(endpoint, Timestamp):
            range_func = date_range
        else:
            range_func = timedelta_range

        breaks = range_func(start=start, end=end, periods=periods, freq=freq)

    return IntervalIndex.from_breaks(breaks, name=name, closed=closed)
