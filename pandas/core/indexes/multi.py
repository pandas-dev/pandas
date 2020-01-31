import datetime
from sys import getsizeof
from typing import Any, Hashable, List, Optional, Sequence, Union
import warnings

import numpy as np

from pandas._config import get_option

from pandas._libs import Timestamp, algos as libalgos, index as libindex, lib, tslibs
from pandas._libs.hashtable import duplicated_int64
from pandas.compat.numpy import function as nv
from pandas.errors import PerformanceWarning, UnsortedIndexError
from pandas.util._decorators import Appender, cache_readonly

from pandas.core.dtypes.cast import coerce_indexer_dtype
from pandas.core.dtypes.common import (
    ensure_int64,
    ensure_platform_int,
    is_categorical_dtype,
    is_hashable,
    is_integer,
    is_iterator,
    is_list_like,
    is_object_dtype,
    is_scalar,
    pandas_dtype,
)
from pandas.core.dtypes.dtypes import ExtensionDtype
from pandas.core.dtypes.generic import ABCDataFrame
from pandas.core.dtypes.missing import array_equivalent, isna

import pandas.core.algorithms as algos
from pandas.core.arrays import Categorical
from pandas.core.arrays.categorical import factorize_from_iterables
import pandas.core.common as com
import pandas.core.indexes.base as ibase
from pandas.core.indexes.base import (
    Index,
    InvalidIndexError,
    _index_shared_docs,
    ensure_index,
)
from pandas.core.indexes.frozen import FrozenList
import pandas.core.missing as missing
from pandas.core.sorting import (
    get_group_index,
    indexer_from_factorized,
    lexsort_indexer,
)
from pandas.core.util.hashing import hash_tuple, hash_tuples

from pandas.io.formats.printing import (
    format_object_attrs,
    format_object_summary,
    pprint_thing,
)

_index_doc_kwargs = dict(ibase._index_doc_kwargs)
_index_doc_kwargs.update(
    dict(klass="MultiIndex", target_klass="MultiIndex or list of tuples")
)


class MultiIndexUIntEngine(libindex.BaseMultiIndexCodesEngine, libindex.UInt64Engine):
    """
    This class manages a MultiIndex by mapping label combinations to positive
    integers.
    """

    _base = libindex.UInt64Engine

    def _codes_to_ints(self, codes):
        """
        Transform combination(s) of uint64 in one uint64 (each), in a strictly
        monotonic way (i.e. respecting the lexicographic order of integer
        combinations): see BaseMultiIndexCodesEngine documentation.

        Parameters
        ----------
        codes : 1- or 2-dimensional array of dtype uint64
            Combinations of integers (one per row)

        Returns
        -------
        scalar or 1-dimensional array, of dtype uint64
            Integer(s) representing one combination (each).
        """
        # Shift the representation of each level by the pre-calculated number
        # of bits:
        codes <<= self.offsets

        # Now sum and OR are in fact interchangeable. This is a simple
        # composition of the (disjunct) significant bits of each level (i.e.
        # each column in "codes") in a single positive integer:
        if codes.ndim == 1:
            # Single key
            return np.bitwise_or.reduce(codes)

        # Multiple keys
        return np.bitwise_or.reduce(codes, axis=1)


class MultiIndexPyIntEngine(libindex.BaseMultiIndexCodesEngine, libindex.ObjectEngine):
    """
    This class manages those (extreme) cases in which the number of possible
    label combinations overflows the 64 bits integers, and uses an ObjectEngine
    containing Python integers.
    """

    _base = libindex.ObjectEngine

    def _codes_to_ints(self, codes):
        """
        Transform combination(s) of uint64 in one Python integer (each), in a
        strictly monotonic way (i.e. respecting the lexicographic order of
        integer combinations): see BaseMultiIndexCodesEngine documentation.

        Parameters
        ----------
        codes : 1- or 2-dimensional array of dtype uint64
            Combinations of integers (one per row)

        Returns
        -------
        int, or 1-dimensional array of dtype object
            Integer(s) representing one combination (each).
        """

        # Shift the representation of each level by the pre-calculated number
        # of bits. Since this can overflow uint64, first make sure we are
        # working with Python integers:
        codes = codes.astype("object") << self.offsets

        # Now sum and OR are in fact interchangeable. This is a simple
        # composition of the (disjunct) significant bits of each level (i.e.
        # each column in "codes") in a single positive integer (per row):
        if codes.ndim == 1:
            # Single key
            return np.bitwise_or.reduce(codes)

        # Multiple keys
        return np.bitwise_or.reduce(codes, axis=1)


class MultiIndex(Index):
    """
    A multi-level, or hierarchical, index object for pandas objects.

    Parameters
    ----------
    levels : sequence of arrays
        The unique labels for each level.
    codes : sequence of arrays
        Integers for each level designating which label at each location.

        .. versionadded:: 0.24.0
    sortorder : optional int
        Level of sortedness (must be lexicographically sorted by that
        level).
    names : optional sequence of objects
        Names for each of the index levels. (name is accepted for compat).
    copy : bool, default False
        Copy the meta-data.
    verify_integrity : bool, default True
        Check that the levels/codes are consistent and valid.

    Attributes
    ----------
    names
    levels
    codes
    nlevels
    levshape

    Methods
    -------
    from_arrays
    from_tuples
    from_product
    from_frame
    set_levels
    set_codes
    to_frame
    to_flat_index
    is_lexsorted
    sortlevel
    droplevel
    swaplevel
    reorder_levels
    remove_unused_levels
    get_locs

    See Also
    --------
    MultiIndex.from_arrays  : Convert list of arrays to MultiIndex.
    MultiIndex.from_product : Create a MultiIndex from the cartesian product
                              of iterables.
    MultiIndex.from_tuples  : Convert list of tuples to a MultiIndex.
    MultiIndex.from_frame   : Make a MultiIndex from a DataFrame.
    Index : The base pandas Index type.

    Notes
    -----
    See the `user guide
    <https://pandas.pydata.org/pandas-docs/stable/user_guide/advanced.html>`_
    for more.

    Examples
    --------
    A new ``MultiIndex`` is typically constructed using one of the helper
    methods :meth:`MultiIndex.from_arrays`, :meth:`MultiIndex.from_product`
    and :meth:`MultiIndex.from_tuples`. For example (using ``.from_arrays``):

    >>> arrays = [[1, 1, 2, 2], ['red', 'blue', 'red', 'blue']]
    >>> pd.MultiIndex.from_arrays(arrays, names=('number', 'color'))
    MultiIndex([(1,  'red'),
                (1, 'blue'),
                (2,  'red'),
                (2, 'blue')],
               names=['number', 'color'])

    See further examples for how to construct a MultiIndex in the doc strings
    of the mentioned helper methods.
    """

    _deprecations = Index._deprecations | frozenset()

    # initialize to zero-length tuples to make everything work
    _typ = "multiindex"
    _names = FrozenList()
    _levels = FrozenList()
    _codes = FrozenList()
    _comparables = ["names"]
    rename = Index.set_names

    # --------------------------------------------------------------------
    # Constructors

    def __new__(
        cls,
        levels=None,
        codes=None,
        sortorder=None,
        names=None,
        dtype=None,
        copy=False,
        name=None,
        verify_integrity: bool = True,
        _set_identity: bool = True,
    ):

        # compat with Index
        if name is not None:
            names = name
        if levels is None or codes is None:
            raise TypeError("Must pass both levels and codes")
        if len(levels) != len(codes):
            raise ValueError("Length of levels and codes must be the same.")
        if len(levels) == 0:
            raise ValueError("Must pass non-zero number of levels/codes")

        result = object.__new__(MultiIndex)

        # we've already validated levels and codes, so shortcut here
        result._set_levels(levels, copy=copy, validate=False)
        result._set_codes(codes, copy=copy, validate=False)

        result._names = [None] * len(levels)
        if names is not None:
            # handles name validation
            result._set_names(names)

        if sortorder is not None:
            result.sortorder = int(sortorder)
        else:
            result.sortorder = sortorder

        if verify_integrity:
            new_codes = result._verify_integrity()
            result._codes = new_codes

        if _set_identity:
            result._reset_identity()

        return result

    def _validate_codes(self, level: List, code: List):
        """
        Reassign code values as -1 if their corresponding levels are NaN.

        Parameters
        ----------
        code : list
            Code to reassign.
        level : list
            Level to check for missing values (NaN, NaT, None).

        Returns
        -------
        new code where code value = -1 if it corresponds
        to a level with missing values (NaN, NaT, None).
        """
        null_mask = isna(level)
        if np.any(null_mask):
            code = np.where(null_mask[code], -1, code)
        return code

    def _verify_integrity(
        self, codes: Optional[List] = None, levels: Optional[List] = None
    ):
        """
        Parameters
        ----------
        codes : optional list
            Codes to check for validity. Defaults to current codes.
        levels : optional list
            Levels to check for validity. Defaults to current levels.

        Raises
        ------
        ValueError
            If length of levels and codes don't match, if the codes for any
            level would exceed level bounds, or there are any duplicate levels.

        Returns
        -------
        new codes where code value = -1 if it corresponds to a
        NaN level.
        """
        # NOTE: Currently does not check, among other things, that cached
        # nlevels matches nor that sortorder matches actually sortorder.
        codes = codes or self.codes
        levels = levels or self.levels

        if len(levels) != len(codes):
            raise ValueError(
                "Length of levels and codes must match. NOTE: "
                "this index is in an inconsistent state."
            )
        codes_length = len(codes[0])
        for i, (level, level_codes) in enumerate(zip(levels, codes)):
            if len(level_codes) != codes_length:
                raise ValueError(
                    f"Unequal code lengths: {[len(code_) for code_ in codes]}"
                )
            if len(level_codes) and level_codes.max() >= len(level):
                raise ValueError(
                    f"On level {i}, code max ({level_codes.max()}) >= length of "
                    f"level ({len(level)}). NOTE: this index is in an "
                    "inconsistent state"
                )
            if len(level_codes) and level_codes.min() < -1:
                raise ValueError(f"On level {i}, code value ({level_codes.min()}) < -1")
            if not level.is_unique:
                raise ValueError(
                    f"Level values must be unique: {list(level)} on level {i}"
                )
        if self.sortorder is not None:
            if self.sortorder > self._lexsort_depth():
                raise ValueError(
                    "Value for sortorder must be inferior or equal to actual "
                    f"lexsort_depth: sortorder {self.sortorder} "
                    f"with lexsort_depth {self._lexsort_depth()}"
                )

        codes = [
            self._validate_codes(level, code) for level, code in zip(levels, codes)
        ]
        new_codes = FrozenList(codes)
        return new_codes

    @classmethod
    def from_arrays(cls, arrays, sortorder=None, names=lib.no_default):
        """
        Convert arrays to MultiIndex.

        Parameters
        ----------
        arrays : list / sequence of array-likes
            Each array-like gives one level's value for each data point.
            len(arrays) is the number of levels.
        sortorder : int or None
            Level of sortedness (must be lexicographically sorted by that
            level).
        names : list / sequence of str, optional
            Names for the levels in the index.

        Returns
        -------
        MultiIndex

        See Also
        --------
        MultiIndex.from_tuples : Convert list of tuples to MultiIndex.
        MultiIndex.from_product : Make a MultiIndex from cartesian product
                                  of iterables.
        MultiIndex.from_frame : Make a MultiIndex from a DataFrame.

        Examples
        --------
        >>> arrays = [[1, 1, 2, 2], ['red', 'blue', 'red', 'blue']]
        >>> pd.MultiIndex.from_arrays(arrays, names=('number', 'color'))
        MultiIndex([(1,  'red'),
                    (1, 'blue'),
                    (2,  'red'),
                    (2, 'blue')],
                   names=['number', 'color'])
        """
        error_msg = "Input must be a list / sequence of array-likes."
        if not is_list_like(arrays):
            raise TypeError(error_msg)
        elif is_iterator(arrays):
            arrays = list(arrays)

        # Check if elements of array are list-like
        for array in arrays:
            if not is_list_like(array):
                raise TypeError(error_msg)

        # Check if lengths of all arrays are equal or not,
        # raise ValueError, if not
        for i in range(1, len(arrays)):
            if len(arrays[i]) != len(arrays[i - 1]):
                raise ValueError("all arrays must be same length")

        codes, levels = factorize_from_iterables(arrays)
        if names is lib.no_default:
            names = [getattr(arr, "name", None) for arr in arrays]

        return MultiIndex(
            levels=levels,
            codes=codes,
            sortorder=sortorder,
            names=names,
            verify_integrity=False,
        )

    @classmethod
    def from_tuples(cls, tuples, sortorder=None, names=None):
        """
        Convert list of tuples to MultiIndex.

        Parameters
        ----------
        tuples : list / sequence of tuple-likes
            Each tuple is the index of one row/column.
        sortorder : int or None
            Level of sortedness (must be lexicographically sorted by that
            level).
        names : list / sequence of str, optional
            Names for the levels in the index.

        Returns
        -------
        MultiIndex

        See Also
        --------
        MultiIndex.from_arrays : Convert list of arrays to MultiIndex.
        MultiIndex.from_product : Make a MultiIndex from cartesian product
                                  of iterables.
        MultiIndex.from_frame : Make a MultiIndex from a DataFrame.

        Examples
        --------
        >>> tuples = [(1, 'red'), (1, 'blue'),
        ...           (2, 'red'), (2, 'blue')]
        >>> pd.MultiIndex.from_tuples(tuples, names=('number', 'color'))
        MultiIndex([(1,  'red'),
                    (1, 'blue'),
                    (2,  'red'),
                    (2, 'blue')],
                   names=['number', 'color'])
        """
        if not is_list_like(tuples):
            raise TypeError("Input must be a list / sequence of tuple-likes.")
        elif is_iterator(tuples):
            tuples = list(tuples)

        if len(tuples) == 0:
            if names is None:
                raise TypeError("Cannot infer number of levels from empty list")
            arrays = [[]] * len(names)
        elif isinstance(tuples, (np.ndarray, Index)):
            if isinstance(tuples, Index):
                tuples = tuples._values

            arrays = list(lib.tuples_to_object_array(tuples).T)
        elif isinstance(tuples, list):
            arrays = list(lib.to_object_array_tuples(tuples).T)
        else:
            arrays = zip(*tuples)

        return MultiIndex.from_arrays(arrays, sortorder=sortorder, names=names)

    @classmethod
    def from_product(cls, iterables, sortorder=None, names=lib.no_default):
        """
        Make a MultiIndex from the cartesian product of multiple iterables.

        Parameters
        ----------
        iterables : list / sequence of iterables
            Each iterable has unique labels for each level of the index.
        sortorder : int or None
            Level of sortedness (must be lexicographically sorted by that
            level).
        names : list / sequence of str, optional
            Names for the levels in the index.

            .. versionchanged:: 1.0.0

               If not explicitly provided, names will be inferred from the
               elements of iterables if an element has a name attribute

        Returns
        -------
        MultiIndex

        See Also
        --------
        MultiIndex.from_arrays : Convert list of arrays to MultiIndex.
        MultiIndex.from_tuples : Convert list of tuples to MultiIndex.
        MultiIndex.from_frame : Make a MultiIndex from a DataFrame.

        Examples
        --------
        >>> numbers = [0, 1, 2]
        >>> colors = ['green', 'purple']
        >>> pd.MultiIndex.from_product([numbers, colors],
        ...                            names=['number', 'color'])
        MultiIndex([(0,  'green'),
                    (0, 'purple'),
                    (1,  'green'),
                    (1, 'purple'),
                    (2,  'green'),
                    (2, 'purple')],
                   names=['number', 'color'])
        """
        from pandas.core.reshape.util import cartesian_product

        if not is_list_like(iterables):
            raise TypeError("Input must be a list / sequence of iterables.")
        elif is_iterator(iterables):
            iterables = list(iterables)

        codes, levels = factorize_from_iterables(iterables)
        if names is lib.no_default:
            names = [getattr(it, "name", None) for it in iterables]

        codes = cartesian_product(codes)
        return MultiIndex(levels, codes, sortorder=sortorder, names=names)

    @classmethod
    def from_frame(cls, df, sortorder=None, names=None):
        """
        Make a MultiIndex from a DataFrame.

        .. versionadded:: 0.24.0

        Parameters
        ----------
        df : DataFrame
            DataFrame to be converted to MultiIndex.
        sortorder : int, optional
            Level of sortedness (must be lexicographically sorted by that
            level).
        names : list-like, optional
            If no names are provided, use the column names, or tuple of column
            names if the columns is a MultiIndex. If a sequence, overwrite
            names with the given sequence.

        Returns
        -------
        MultiIndex
            The MultiIndex representation of the given DataFrame.

        See Also
        --------
        MultiIndex.from_arrays : Convert list of arrays to MultiIndex.
        MultiIndex.from_tuples : Convert list of tuples to MultiIndex.
        MultiIndex.from_product : Make a MultiIndex from cartesian product
                                  of iterables.

        Examples
        --------
        >>> df = pd.DataFrame([['HI', 'Temp'], ['HI', 'Precip'],
        ...                    ['NJ', 'Temp'], ['NJ', 'Precip']],
        ...                   columns=['a', 'b'])
        >>> df
              a       b
        0    HI    Temp
        1    HI  Precip
        2    NJ    Temp
        3    NJ  Precip

        >>> pd.MultiIndex.from_frame(df)
        MultiIndex([('HI',   'Temp'),
                    ('HI', 'Precip'),
                    ('NJ',   'Temp'),
                    ('NJ', 'Precip')],
                   names=['a', 'b'])

        Using explicit names, instead of the column names

        >>> pd.MultiIndex.from_frame(df, names=['state', 'observation'])
        MultiIndex([('HI',   'Temp'),
                    ('HI', 'Precip'),
                    ('NJ',   'Temp'),
                    ('NJ', 'Precip')],
                   names=['state', 'observation'])
        """
        if not isinstance(df, ABCDataFrame):
            raise TypeError("Input must be a DataFrame")

        column_names, columns = zip(*df.items())
        names = column_names if names is None else names
        return cls.from_arrays(columns, sortorder=sortorder, names=names)

    # --------------------------------------------------------------------

    @property
    def levels(self):
        result = [
            x._shallow_copy(name=name) for x, name in zip(self._levels, self._names)
        ]
        for level in result:
            # disallow midx.levels[0].name = "foo"
            level._no_setting_name = True
        return FrozenList(result)

    @property
    def _values(self):
        # We override here, since our parent uses _data, which we don't use.
        return self.values

    @property
    def shape(self):
        """
        Return a tuple of the shape of the underlying data.
        """
        # overriding the base Index.shape definition to avoid materializing
        # the values (GH-27384, GH-27775)
        return (len(self),)

    @property
    def array(self):
        """
        Raises a ValueError for `MultiIndex` because there's no single
        array backing a MultiIndex.

        Raises
        ------
        ValueError
        """
        raise ValueError(
            "MultiIndex has no single backing array. Use "
            "'MultiIndex.to_numpy()' to get a NumPy array of tuples."
        )

    def _set_levels(
        self, levels, level=None, copy=False, validate=True, verify_integrity=False
    ):
        # This is NOT part of the levels property because it should be
        # externally not allowed to set levels. User beware if you change
        # _levels directly
        if validate:
            if len(levels) == 0:
                raise ValueError("Must set non-zero number of levels.")
            if level is None and len(levels) != self.nlevels:
                raise ValueError("Length of levels must match number of levels.")
            if level is not None and len(levels) != len(level):
                raise ValueError("Length of levels must match length of level.")

        if level is None:
            new_levels = FrozenList(
                ensure_index(lev, copy=copy)._shallow_copy() for lev in levels
            )
        else:
            level_numbers = [self._get_level_number(lev) for lev in level]
            new_levels = list(self._levels)
            for lev_num, lev in zip(level_numbers, levels):
                new_levels[lev_num] = ensure_index(lev, copy=copy)._shallow_copy()
            new_levels = FrozenList(new_levels)

        if verify_integrity:
            new_codes = self._verify_integrity(levels=new_levels)
            self._codes = new_codes

        names = self.names
        self._levels = new_levels
        if any(names):
            self._set_names(names)

        self._tuples = None
        self._reset_cache()

    def set_levels(self, levels, level=None, inplace=False, verify_integrity=True):
        """
        Set new levels on MultiIndex. Defaults to returning new index.

        Parameters
        ----------
        levels : sequence or list of sequence
            New level(s) to apply.
        level : int, level name, or sequence of int/level names (default None)
            Level(s) to set (None for all levels).
        inplace : bool
            If True, mutates in place.
        verify_integrity : bool, default True
            If True, checks that levels and codes are compatible.

        Returns
        -------
        new index (of same type and class...etc)

        Examples
        --------
        >>> idx = pd.MultiIndex.from_tuples([(1, 'one'), (1, 'two'),
                                            (2, 'one'), (2, 'two'),
                                            (3, 'one'), (3, 'two')],
                                            names=['foo', 'bar'])
        >>> idx.set_levels([['a', 'b', 'c'], [1, 2]])
        MultiIndex([('a', 1),
                    ('a', 2),
                    ('b', 1),
                    ('b', 2),
                    ('c', 1),
                    ('c', 2)],
                   names=['foo', 'bar'])
        >>> idx.set_levels(['a', 'b', 'c'], level=0)
        MultiIndex([('a', 'one'),
                    ('a', 'two'),
                    ('b', 'one'),
                    ('b', 'two'),
                    ('c', 'one'),
                    ('c', 'two')],
                   names=['foo', 'bar'])
        >>> idx.set_levels(['a', 'b'], level='bar')
        MultiIndex([(1, 'a'),
                    (1, 'b'),
                    (2, 'a'),
                    (2, 'b'),
                    (3, 'a'),
                    (3, 'b')],
                   names=['foo', 'bar'])

        If any of the levels passed to ``set_levels()`` exceeds the
        existing length, all of the values from that argument will
        be stored in the MultiIndex levels, though the values will
        be truncated in the MultiIndex output.

        >>> idx.set_levels([['a', 'b', 'c'], [1, 2, 3, 4]], level=[0, 1])
        MultiIndex([('a', 1),
                    ('a', 2),
                    ('b', 1),
                    ('b', 2)],
                   names=['foo', 'bar'])
        >>> idx.set_levels([['a', 'b', 'c'], [1, 2, 3, 4]], level=[0, 1]).levels
        FrozenList([['a', 'b', 'c'], [1, 2, 3, 4]])
        """
        if is_list_like(levels) and not isinstance(levels, Index):
            levels = list(levels)

        if level is not None and not is_list_like(level):
            if not is_list_like(levels):
                raise TypeError("Levels must be list-like")
            if is_list_like(levels[0]):
                raise TypeError("Levels must be list-like")
            level = [level]
            levels = [levels]
        elif level is None or is_list_like(level):
            if not is_list_like(levels) or not is_list_like(levels[0]):
                raise TypeError("Levels must be list of lists-like")

        if inplace:
            idx = self
        else:
            idx = self._shallow_copy()
        idx._reset_identity()
        idx._set_levels(
            levels, level=level, validate=True, verify_integrity=verify_integrity
        )
        if not inplace:
            return idx

    @property
    def codes(self):
        return self._codes

    def _set_codes(
        self, codes, level=None, copy=False, validate=True, verify_integrity=False
    ):
        if validate:
            if level is None and len(codes) != self.nlevels:
                raise ValueError("Length of codes must match number of levels")
            if level is not None and len(codes) != len(level):
                raise ValueError("Length of codes must match length of levels.")

        if level is None:
            new_codes = FrozenList(
                _coerce_indexer_frozen(level_codes, lev, copy=copy).view()
                for lev, level_codes in zip(self._levels, codes)
            )
        else:
            level_numbers = [self._get_level_number(lev) for lev in level]
            new_codes = list(self._codes)
            for lev_num, level_codes in zip(level_numbers, codes):
                lev = self.levels[lev_num]
                new_codes[lev_num] = _coerce_indexer_frozen(level_codes, lev, copy=copy)
            new_codes = FrozenList(new_codes)

        if verify_integrity:
            new_codes = self._verify_integrity(codes=new_codes)

        self._codes = new_codes

        self._tuples = None
        self._reset_cache()

    def set_codes(self, codes, level=None, inplace=False, verify_integrity=True):
        """
        Set new codes on MultiIndex. Defaults to returning
        new index.

        .. versionadded:: 0.24.0

           New name for deprecated method `set_labels`.

        Parameters
        ----------
        codes : sequence or list of sequence
            New codes to apply.
        level : int, level name, or sequence of int/level names (default None)
            Level(s) to set (None for all levels).
        inplace : bool
            If True, mutates in place.
        verify_integrity : bool (default True)
            If True, checks that levels and codes are compatible.

        Returns
        -------
        new index (of same type and class...etc)

        Examples
        --------
        >>> idx = pd.MultiIndex.from_tuples([(1, 'one'),
                                             (1, 'two'),
                                             (2, 'one'),
                                             (2, 'two')],
                                            names=['foo', 'bar'])
        >>> idx.set_codes([[1, 0, 1, 0], [0, 0, 1, 1]])
        MultiIndex([(2, 'one'),
                    (1, 'one'),
                    (2, 'two'),
                    (1, 'two')],
                   names=['foo', 'bar'])
        >>> idx.set_codes([1, 0, 1, 0], level=0)
        MultiIndex([(2, 'one'),
                    (1, 'two'),
                    (2, 'one'),
                    (1, 'two')],
                   names=['foo', 'bar'])
        >>> idx.set_codes([0, 0, 1, 1], level='bar')
        MultiIndex([(1, 'one'),
                    (1, 'one'),
                    (2, 'two'),
                    (2, 'two')],
                   names=['foo', 'bar'])
        >>> idx.set_codes([[1, 0, 1, 0], [0, 0, 1, 1]], level=[0, 1])
        MultiIndex([(2, 'one'),
                    (1, 'one'),
                    (2, 'two'),
                    (1, 'two')],
                   names=['foo', 'bar'])
        """
        if level is not None and not is_list_like(level):
            if not is_list_like(codes):
                raise TypeError("Codes must be list-like")
            if is_list_like(codes[0]):
                raise TypeError("Codes must be list-like")
            level = [level]
            codes = [codes]
        elif level is None or is_list_like(level):
            if not is_list_like(codes) or not is_list_like(codes[0]):
                raise TypeError("Codes must be list of lists-like")

        if inplace:
            idx = self
        else:
            idx = self._shallow_copy()
        idx._reset_identity()
        idx._set_codes(codes, level=level, verify_integrity=verify_integrity)
        if not inplace:
            return idx

    def copy(
        self,
        names=None,
        dtype=None,
        levels=None,
        codes=None,
        deep=False,
        _set_identity=False,
        **kwargs,
    ):
        """
        Make a copy of this object. Names, dtype, levels and codes can be
        passed and will be set on new copy.

        Parameters
        ----------
        names : sequence, optional
        dtype : numpy dtype or pandas type, optional
        levels : sequence, optional
        codes : sequence, optional

        Returns
        -------
        copy : MultiIndex

        Notes
        -----
        In most cases, there should be no functional difference from using
        ``deep``, but if ``deep`` is passed it will attempt to deepcopy.
        This could be potentially expensive on large MultiIndex objects.
        """
        name = kwargs.get("name")
        names = self._validate_names(name=name, names=names, deep=deep)
        if "labels" in kwargs:
            raise TypeError("'labels' argument has been removed; use 'codes' instead")
        if deep:
            from copy import deepcopy

            if levels is None:
                levels = deepcopy(self.levels)
            if codes is None:
                codes = deepcopy(self.codes)
        else:
            if levels is None:
                levels = self.levels
            if codes is None:
                codes = self.codes
        return MultiIndex(
            levels=levels,
            codes=codes,
            names=names,
            sortorder=self.sortorder,
            verify_integrity=False,
            _set_identity=_set_identity,
        )

    def __array__(self, dtype=None) -> np.ndarray:
        """ the array interface, return my values """
        return self.values

    def view(self, cls=None):
        """ this is defined as a copy with the same identity """
        result = self.copy()
        result._id = self._id
        return result

    def _shallow_copy_with_infer(self, values, **kwargs):
        # On equal MultiIndexes the difference is empty.
        # Therefore, an empty MultiIndex is returned GH13490
        if len(values) == 0:
            return MultiIndex(
                levels=[[] for _ in range(self.nlevels)],
                codes=[[] for _ in range(self.nlevels)],
                **kwargs,
            )
        return self._shallow_copy(values, **kwargs)

    @Appender(Index.__contains__.__doc__)
    def __contains__(self, key: Any) -> bool:
        hash(key)
        try:
            self.get_loc(key)
            return True
        except (LookupError, TypeError, ValueError):
            return False

    @Appender(Index._shallow_copy.__doc__)
    def _shallow_copy(self, values=None, **kwargs):
        if values is not None:
            names = kwargs.pop("names", kwargs.pop("name", self.names))
            # discards freq
            kwargs.pop("freq", None)
            return MultiIndex.from_tuples(values, names=names, **kwargs)
        return self.copy(**kwargs)

    @cache_readonly
    def dtype(self) -> np.dtype:
        return np.dtype("O")

    def _is_memory_usage_qualified(self) -> bool:
        """ return a boolean if we need a qualified .info display """

        def f(l):
            return "mixed" in l or "string" in l or "unicode" in l

        return any(f(l) for l in self._inferred_type_levels)

    @Appender(Index.memory_usage.__doc__)
    def memory_usage(self, deep: bool = False) -> int:
        # we are overwriting our base class to avoid
        # computing .values here which could materialize
        # a tuple representation unnecessarily
        return self._nbytes(deep)

    @cache_readonly
    def nbytes(self) -> int:
        """ return the number of bytes in the underlying data """
        return self._nbytes(False)

    def _nbytes(self, deep: bool = False) -> int:
        """
        return the number of bytes in the underlying data
        deeply introspect the level data if deep=True

        include the engine hashtable

        *this is in internal routine*

        """

        # for implementations with no useful getsizeof (PyPy)
        objsize = 24

        level_nbytes = sum(i.memory_usage(deep=deep) for i in self.levels)
        label_nbytes = sum(i.nbytes for i in self.codes)
        names_nbytes = sum(getsizeof(i, objsize) for i in self.names)
        result = level_nbytes + label_nbytes + names_nbytes

        # include our engine hashtable
        result += self._engine.sizeof(deep=deep)
        return result

    # --------------------------------------------------------------------
    # Rendering Methods
    def _formatter_func(self, tup):
        """
        Formats each item in tup according to its level's formatter function.
        """
        formatter_funcs = [level._formatter_func for level in self.levels]
        return tuple(func(val) for func, val in zip(formatter_funcs, tup))

    def _format_data(self, name=None):
        """
        Return the formatted data as a unicode string
        """
        return format_object_summary(
            self, self._formatter_func, name=name, line_break_each_value=True
        )

    def _format_attrs(self):
        """
        Return a list of tuples of the (attr,formatted_value).
        """
        return format_object_attrs(self, include_dtype=False)

    def _format_native_types(self, na_rep="nan", **kwargs):
        new_levels = []
        new_codes = []

        # go through the levels and format them
        for level, level_codes in zip(self.levels, self.codes):
            level = level._format_native_types(na_rep=na_rep, **kwargs)
            # add nan values, if there are any
            mask = level_codes == -1
            if mask.any():
                nan_index = len(level)
                level = np.append(level, na_rep)
                assert not level_codes.flags.writeable  # i.e. copy is needed
                level_codes = level_codes.copy()  # make writeable
                level_codes[mask] = nan_index
            new_levels.append(level)
            new_codes.append(level_codes)

        if len(new_levels) == 1:
            # a single-level multi-index
            return Index(new_levels[0].take(new_codes[0]))._format_native_types()
        else:
            # reconstruct the multi-index
            mi = MultiIndex(
                levels=new_levels,
                codes=new_codes,
                names=self.names,
                sortorder=self.sortorder,
                verify_integrity=False,
            )
            return mi.values

    def format(
        self,
        space=2,
        sparsify=None,
        adjoin=True,
        names=False,
        na_rep=None,
        formatter=None,
    ):
        if len(self) == 0:
            return []

        stringified_levels = []
        for lev, level_codes in zip(self.levels, self.codes):
            na = na_rep if na_rep is not None else _get_na_rep(lev.dtype.type)

            if len(lev) > 0:

                formatted = lev.take(level_codes).format(formatter=formatter)

                # we have some NA
                mask = level_codes == -1
                if mask.any():
                    formatted = np.array(formatted, dtype=object)
                    formatted[mask] = na
                    formatted = formatted.tolist()

            else:
                # weird all NA case
                formatted = [
                    pprint_thing(na if isna(x) else x, escape_chars=("\t", "\r", "\n"))
                    for x in algos.take_1d(lev._values, level_codes)
                ]
            stringified_levels.append(formatted)

        result_levels = []
        for lev, name in zip(stringified_levels, self.names):
            level = []

            if names:
                level.append(
                    pprint_thing(name, escape_chars=("\t", "\r", "\n"))
                    if name is not None
                    else ""
                )

            level.extend(np.array(lev, dtype=object))
            result_levels.append(level)

        if sparsify is None:
            sparsify = get_option("display.multi_sparse")

        if sparsify:
            sentinel = ""
            # GH3547
            # use value of sparsify as sentinel,  unless it's an obvious
            # "Truthy" value
            if sparsify not in [True, 1]:
                sentinel = sparsify
            # little bit of a kludge job for #1217
            result_levels = _sparsify(
                result_levels, start=int(names), sentinel=sentinel
            )

        if adjoin:
            from pandas.io.formats.format import _get_adjustment

            adj = _get_adjustment()
            return adj.adjoin(space, *result_levels).split("\n")
        else:
            return result_levels

    # --------------------------------------------------------------------

    def __len__(self) -> int:
        return len(self.codes[0])

    def _get_names(self):
        return FrozenList(self._names)

    def _set_names(self, names, level=None, validate=True):
        """
        Set new names on index. Each name has to be a hashable type.

        Parameters
        ----------
        values : str or sequence
            name(s) to set
        level : int, level name, or sequence of int/level names (default None)
            If the index is a MultiIndex (hierarchical), level(s) to set (None
            for all levels).  Otherwise level must be None
        validate : boolean, default True
            validate that the names match level lengths

        Raises
        ------
        TypeError if each name is not hashable.

        Notes
        -----
        sets names on levels. WARNING: mutates!

        Note that you generally want to set this *after* changing levels, so
        that it only acts on copies
        """
        # GH 15110
        # Don't allow a single string for names in a MultiIndex
        if names is not None and not is_list_like(names):
            raise ValueError("Names should be list-like for a MultiIndex")
        names = list(names)

        if validate:
            if level is not None and len(names) != len(level):
                raise ValueError("Length of names must match length of level.")
            if level is None and len(names) != self.nlevels:
                raise ValueError(
                    "Length of names must match number of levels in MultiIndex."
                )

        if level is None:
            level = range(self.nlevels)
        else:
            level = [self._get_level_number(lev) for lev in level]

        # set the name
        for lev, name in zip(level, names):
            if name is not None:
                # GH 20527
                # All items in 'names' need to be hashable:
                if not is_hashable(name):
                    raise TypeError(
                        f"{type(self).__name__}.name must be a hashable type"
                    )
            self._names[lev] = name

    names = property(
        fset=_set_names, fget=_get_names, doc="""\nNames of levels in MultiIndex.\n"""
    )

    @Appender(Index._get_grouper_for_level.__doc__)
    def _get_grouper_for_level(self, mapper, level):
        indexer = self.codes[level]
        level_index = self.levels[level]

        if mapper is not None:
            # Handle group mapping function and return
            level_values = self.levels[level].take(indexer)
            grouper = level_values.map(mapper)
            return grouper, None, None

        codes, uniques = algos.factorize(indexer, sort=True)

        if len(uniques) > 0 and uniques[0] == -1:
            # Handle NAs
            mask = indexer != -1
            ok_codes, uniques = algos.factorize(indexer[mask], sort=True)

            codes = np.empty(len(indexer), dtype=indexer.dtype)
            codes[mask] = ok_codes
            codes[~mask] = -1

        if len(uniques) < len(level_index):
            # Remove unobserved levels from level_index
            level_index = level_index.take(uniques)
        else:
            # break references back to us so that setting the name
            # on the output of a groupby doesn't reflect back here.
            level_index = level_index.copy()

        if level_index._can_hold_na:
            grouper = level_index.take(codes, fill_value=True)
        else:
            grouper = level_index.take(codes)

        return grouper, codes, level_index

    @property
    def _constructor(self):
        return MultiIndex.from_tuples

    @cache_readonly
    def inferred_type(self) -> str:
        return "mixed"

    def _get_level_number(self, level) -> int:
        count = self.names.count(level)
        if (count > 1) and not is_integer(level):
            raise ValueError(
                f"The name {level} occurs multiple times, use a level number"
            )
        try:
            level = self.names.index(level)
        except ValueError:
            if not is_integer(level):
                raise KeyError(f"Level {level} not found")
            elif level < 0:
                level += self.nlevels
                if level < 0:
                    orig_level = level - self.nlevels
                    raise IndexError(
                        f"Too many levels: Index has only {self.nlevels} levels, "
                        f"{orig_level} is not a valid level number"
                    )
            # Note: levels are zero-based
            elif level >= self.nlevels:
                raise IndexError(
                    f"Too many levels: Index has only {self.nlevels} levels, "
                    f"not {level + 1}"
                )
        return level

    _tuples = None

    @cache_readonly
    def _engine(self):
        # Calculate the number of bits needed to represent labels in each
        # level, as log2 of their sizes (including -1 for NaN):
        sizes = np.ceil(np.log2([len(l) + 1 for l in self.levels]))

        # Sum bit counts, starting from the _right_....
        lev_bits = np.cumsum(sizes[::-1])[::-1]

        # ... in order to obtain offsets such that sorting the combination of
        # shifted codes (one for each level, resulting in a unique integer) is
        # equivalent to sorting lexicographically the codes themselves. Notice
        # that each level needs to be shifted by the number of bits needed to
        # represent the _previous_ ones:
        offsets = np.concatenate([lev_bits[1:], [0]]).astype("uint64")

        # Check the total number of bits needed for our representation:
        if lev_bits[0] > 64:
            # The levels would overflow a 64 bit uint - use Python integers:
            return MultiIndexPyIntEngine(self.levels, self.codes, offsets)
        return MultiIndexUIntEngine(self.levels, self.codes, offsets)

    @property
    def values(self):
        if self._tuples is not None:
            return self._tuples

        values = []

        for i in range(self.nlevels):
            vals = self._get_level_values(i)
            if is_categorical_dtype(vals):
                vals = vals._internal_get_values()
            if isinstance(vals.dtype, ExtensionDtype) or hasattr(vals, "_box_values"):
                vals = vals.astype(object)
            vals = np.array(vals, copy=False)
            values.append(vals)

        self._tuples = lib.fast_zip(values)
        return self._tuples

    @property
    def _has_complex_internals(self) -> bool:
        # used to avoid libreduction code paths, which raise or require conversion
        return True

    @cache_readonly
    def is_monotonic_increasing(self) -> bool:
        """
        return if the index is monotonic increasing (only equal or
        increasing) values.
        """

        if all(x.is_monotonic for x in self.levels):
            # If each level is sorted, we can operate on the codes directly. GH27495
            return libalgos.is_lexsorted(
                [x.astype("int64", copy=False) for x in self.codes]
            )

        # reversed() because lexsort() wants the most significant key last.
        values = [
            self._get_level_values(i).values for i in reversed(range(len(self.levels)))
        ]
        try:
            sort_order = np.lexsort(values)
            return Index(sort_order).is_monotonic
        except TypeError:

            # we have mixed types and np.lexsort is not happy
            return Index(self.values).is_monotonic

    @cache_readonly
    def is_monotonic_decreasing(self) -> bool:
        """
        return if the index is monotonic decreasing (only equal or
        decreasing) values.
        """
        # monotonic decreasing if and only if reverse is monotonic increasing
        return self[::-1].is_monotonic_increasing

    @cache_readonly
    def _have_mixed_levels(self):
        """ return a boolean list indicated if we have mixed levels """
        return ["mixed" in l for l in self._inferred_type_levels]

    @cache_readonly
    def _inferred_type_levels(self):
        """ return a list of the inferred types, one for each level """
        return [i.inferred_type for i in self.levels]

    @cache_readonly
    def _hashed_values(self):
        """ return a uint64 ndarray of my hashed values """
        return hash_tuples(self)

    def _hashed_indexing_key(self, key):
        """
        validate and return the hash for the provided key

        *this is internal for use for the cython routines*

        Parameters
        ----------
        key : string or tuple

        Returns
        -------
        np.uint64

        Notes
        -----
        we need to stringify if we have mixed levels
        """

        if not isinstance(key, tuple):
            return hash_tuples(key)

        if not len(key) == self.nlevels:
            raise KeyError

        def f(k, stringify):
            if stringify and not isinstance(k, str):
                k = str(k)
            return k

        key = tuple(
            f(k, stringify) for k, stringify in zip(key, self._have_mixed_levels)
        )
        return hash_tuple(key)

    @Appender(Index.duplicated.__doc__)
    def duplicated(self, keep="first"):
        shape = map(len, self.levels)
        ids = get_group_index(self.codes, shape, sort=False, xnull=False)

        return duplicated_int64(ids, keep)

    def fillna(self, value=None, downcast=None):
        """
        fillna is not implemented for MultiIndex
        """
        raise NotImplementedError("isna is not defined for MultiIndex")

    @Appender(Index.dropna.__doc__)
    def dropna(self, how="any"):
        nans = [level_codes == -1 for level_codes in self.codes]
        if how == "any":
            indexer = np.any(nans, axis=0)
        elif how == "all":
            indexer = np.all(nans, axis=0)
        else:
            raise ValueError(f"invalid how option: {how}")

        new_codes = [level_codes[~indexer] for level_codes in self.codes]
        return self.copy(codes=new_codes, deep=True)

    def get_value(self, series, key):
        # Label-based
        s = com.values_from_object(series)
        k = com.values_from_object(key)

        def _try_mi(k):
            # TODO: what if a level contains tuples??
            loc = self.get_loc(k)
            new_values = series._values[loc]
            new_index = self[loc]
            new_index = maybe_droplevels(new_index, k)
            return series._constructor(
                new_values, index=new_index, name=series.name
            ).__finalize__(self)

        try:
            return self._engine.get_value(s, k)
        except KeyError as e1:
            try:
                return _try_mi(key)
            except KeyError:
                pass

            try:
                return libindex.get_value_at(s, k)
            except IndexError:
                raise
            except TypeError:
                # generator/iterator-like
                if is_iterator(key):
                    raise InvalidIndexError(key)
                else:
                    raise e1
            except Exception:  # pragma: no cover
                raise e1
        except TypeError:

            # a Timestamp will raise a TypeError in a multi-index
            # rather than a KeyError, try it here
            # note that a string that 'looks' like a Timestamp will raise
            # a KeyError! (GH5725)
            if isinstance(key, (datetime.datetime, np.datetime64, str)):
                try:
                    return _try_mi(key)
                except KeyError:
                    raise
                except (IndexError, ValueError, TypeError):
                    pass

                try:
                    return _try_mi(Timestamp(key))
                except (
                    KeyError,
                    TypeError,
                    IndexError,
                    ValueError,
                    tslibs.OutOfBoundsDatetime,
                ):
                    pass

            raise InvalidIndexError(key)

    def _get_level_values(self, level, unique=False):
        """
        Return vector of label values for requested level,
        equal to the length of the index

        **this is an internal method**

        Parameters
        ----------
        level : int level
        unique : bool, default False
            if True, drop duplicated values

        Returns
        -------
        values : ndarray
        """

        lev = self.levels[level]
        level_codes = self.codes[level]
        name = self._names[level]
        if unique:
            level_codes = algos.unique(level_codes)
        filled = algos.take_1d(lev._values, level_codes, fill_value=lev._na_value)
        return lev._shallow_copy(filled, name=name)

    def get_level_values(self, level):
        """
        Return vector of label values for requested level,
        equal to the length of the index.

        Parameters
        ----------
        level : int or str
            ``level`` is either the integer position of the level in the
            MultiIndex, or the name of the level.

        Returns
        -------
        values : Index
            Values is a level of this MultiIndex converted to
            a single :class:`Index` (or subclass thereof).

        Examples
        --------

        Create a MultiIndex:

        >>> mi = pd.MultiIndex.from_arrays((list('abc'), list('def')))
        >>> mi.names = ['level_1', 'level_2']

        Get level values by supplying level as either integer or name:

        >>> mi.get_level_values(0)
        Index(['a', 'b', 'c'], dtype='object', name='level_1')
        >>> mi.get_level_values('level_2')
        Index(['d', 'e', 'f'], dtype='object', name='level_2')
        """
        level = self._get_level_number(level)
        values = self._get_level_values(level)
        return values

    @Appender(Index.unique.__doc__)
    def unique(self, level=None):

        if level is None:
            return super().unique()
        else:
            level = self._get_level_number(level)
            return self._get_level_values(level=level, unique=True)

    def _to_safe_for_reshape(self):
        """ convert to object if we are a categorical """
        return self.set_levels([i._to_safe_for_reshape() for i in self.levels])

    def to_frame(self, index=True, name=None):
        """
        Create a DataFrame with the levels of the MultiIndex as columns.

        Column ordering is determined by the DataFrame constructor with data as
        a dict.

        .. versionadded:: 0.24.0

        Parameters
        ----------
        index : bool, default True
            Set the index of the returned DataFrame as the original MultiIndex.

        name : list / sequence of strings, optional
            The passed names should substitute index level names.

        Returns
        -------
        DataFrame : a DataFrame containing the original MultiIndex data.

        See Also
        --------
        DataFrame
        """

        from pandas import DataFrame

        if name is not None:
            if not is_list_like(name):
                raise TypeError("'name' must be a list / sequence of column names.")

            if len(name) != len(self.levels):
                raise ValueError(
                    "'name' should have same length as number of levels on index."
                )
            idx_names = name
        else:
            idx_names = self.names

        # Guarantee resulting column order - PY36+ dict maintains insertion order
        result = DataFrame(
            {
                (level if lvlname is None else lvlname): self._get_level_values(level)
                for lvlname, level in zip(idx_names, range(len(self.levels)))
            },
            copy=False,
        )

        if index:
            result.index = self
        return result

    def to_flat_index(self):
        """
        Convert a MultiIndex to an Index of Tuples containing the level values.

        .. versionadded:: 0.24.0

        Returns
        -------
        pd.Index
            Index with the MultiIndex data represented in Tuples.

        Notes
        -----
        This method will simply return the caller if called by anything other
        than a MultiIndex.

        Examples
        --------
        >>> index = pd.MultiIndex.from_product(
        ...     [['foo', 'bar'], ['baz', 'qux']],
        ...     names=['a', 'b'])
        >>> index.to_flat_index()
        Index([('foo', 'baz'), ('foo', 'qux'),
               ('bar', 'baz'), ('bar', 'qux')],
              dtype='object')
        """
        return Index(self.values, tupleize_cols=False)

    @property
    def is_all_dates(self) -> bool:
        return False

    def is_lexsorted(self) -> bool:
        """
        Return True if the codes are lexicographically sorted.

        Returns
        -------
        bool
        """
        return self.lexsort_depth == self.nlevels

    @cache_readonly
    def lexsort_depth(self):
        if self.sortorder is not None:
            return self.sortorder

        return self._lexsort_depth()

    def _lexsort_depth(self) -> int:
        """
        Compute and return the lexsort_depth, the number of levels of the
        MultiIndex that are sorted lexically

        Returns
        ------
        int
        """
        int64_codes = [ensure_int64(level_codes) for level_codes in self.codes]
        for k in range(self.nlevels, 0, -1):
            if libalgos.is_lexsorted(int64_codes[:k]):
                return k
        return 0

    def _sort_levels_monotonic(self):
        """
        This is an *internal* function.

        Create a new MultiIndex from the current to monotonically sorted
        items IN the levels. This does not actually make the entire MultiIndex
        monotonic, JUST the levels.

        The resulting MultiIndex will have the same outward
        appearance, meaning the same .values and ordering. It will also
        be .equals() to the original.

        Returns
        -------
        MultiIndex

        Examples
        --------

        >>> mi = pd.MultiIndex(levels=[['a', 'b'], ['bb', 'aa']],
        ...                    codes=[[0, 0, 1, 1], [0, 1, 0, 1]])
        >>> mi
        MultiIndex([('a', 'bb'),
                    ('a', 'aa'),
                    ('b', 'bb'),
                    ('b', 'aa')],
                   )

        >>> mi.sort_values()
        MultiIndex([('a', 'aa'),
                    ('a', 'bb'),
                    ('b', 'aa'),
                    ('b', 'bb')],
                   )
        """

        if self.is_lexsorted() and self.is_monotonic:
            return self

        new_levels = []
        new_codes = []

        for lev, level_codes in zip(self.levels, self.codes):

            if not lev.is_monotonic:
                try:
                    # indexer to reorder the levels
                    indexer = lev.argsort()
                except TypeError:
                    pass
                else:
                    lev = lev.take(indexer)

                    # indexer to reorder the level codes
                    indexer = ensure_int64(indexer)
                    ri = lib.get_reverse_indexer(indexer, len(indexer))
                    level_codes = algos.take_1d(ri, level_codes)

            new_levels.append(lev)
            new_codes.append(level_codes)

        return MultiIndex(
            new_levels,
            new_codes,
            names=self.names,
            sortorder=self.sortorder,
            verify_integrity=False,
        )

    def remove_unused_levels(self):
        """
        Create a new MultiIndex from the current that removes
        unused levels, meaning that they are not expressed in the labels.

        The resulting MultiIndex will have the same outward
        appearance, meaning the same .values and ordering. It will also
        be .equals() to the original.

        Returns
        -------
        MultiIndex

        Examples
        --------
        >>> mi = pd.MultiIndex.from_product([range(2), list('ab')])
        >>> mi
        MultiIndex([(0, 'a'),
                    (0, 'b'),
                    (1, 'a'),
                    (1, 'b')],
                   )

        >>> mi[2:]
        MultiIndex([(1, 'a'),
                    (1, 'b')],
                   )

        The 0 from the first level is not represented
        and can be removed

        >>> mi2 = mi[2:].remove_unused_levels()
        >>> mi2.levels
        FrozenList([[1], ['a', 'b']])
        """

        new_levels = []
        new_codes = []

        changed = False
        for lev, level_codes in zip(self.levels, self.codes):

            # Since few levels are typically unused, bincount() is more
            # efficient than unique() - however it only accepts positive values
            # (and drops order):
            uniques = np.where(np.bincount(level_codes + 1) > 0)[0] - 1
            has_na = int(len(uniques) and (uniques[0] == -1))

            if len(uniques) != len(lev) + has_na:
                # We have unused levels
                changed = True

                # Recalculate uniques, now preserving order.
                # Can easily be cythonized by exploiting the already existing
                # "uniques" and stop parsing "level_codes" when all items
                # are found:
                uniques = algos.unique(level_codes)
                if has_na:
                    na_idx = np.where(uniques == -1)[0]
                    # Just ensure that -1 is in first position:
                    uniques[[0, na_idx[0]]] = uniques[[na_idx[0], 0]]

                # codes get mapped from uniques to 0:len(uniques)
                # -1 (if present) is mapped to last position
                code_mapping = np.zeros(len(lev) + has_na)
                # ... and reassigned value -1:
                code_mapping[uniques] = np.arange(len(uniques)) - has_na

                level_codes = code_mapping[level_codes]

                # new levels are simple
                lev = lev.take(uniques[has_na:])

            new_levels.append(lev)
            new_codes.append(level_codes)

        result = self.view()

        if changed:
            result._reset_identity()
            result._set_levels(new_levels, validate=False)
            result._set_codes(new_codes, validate=False)

        return result

    @property
    def nlevels(self) -> int:
        """
        Integer number of levels in this MultiIndex.
        """
        return len(self._levels)

    @property
    def levshape(self):
        """
        A tuple with the length of each level.
        """
        return tuple(len(x) for x in self.levels)

    def __reduce__(self):
        """Necessary for making this object picklable"""
        d = dict(
            levels=list(self.levels),
            codes=list(self.codes),
            sortorder=self.sortorder,
            names=list(self.names),
        )
        return ibase._new_Index, (type(self), d), None

    def __setstate__(self, state):
        """Necessary for making this object picklable"""

        if isinstance(state, dict):
            levels = state.get("levels")
            codes = state.get("codes")
            sortorder = state.get("sortorder")
            names = state.get("names")

        elif isinstance(state, tuple):

            nd_state, own_state = state
            levels, codes, sortorder, names = own_state

        self._set_levels([Index(x) for x in levels], validate=False)
        self._set_codes(codes)
        new_codes = self._verify_integrity()
        self._set_codes(new_codes)
        self._set_names(names)
        self.sortorder = sortorder
        self._reset_identity()

    def __getitem__(self, key):
        if is_scalar(key):
            key = com.cast_scalar_indexer(key)

            retval = []
            for lev, level_codes in zip(self.levels, self.codes):
                if level_codes[key] == -1:
                    retval.append(np.nan)
                else:
                    retval.append(lev[level_codes[key]])

            return tuple(retval)
        else:
            if com.is_bool_indexer(key):
                key = np.asarray(key, dtype=bool)
                sortorder = self.sortorder
            else:
                # cannot be sure whether the result will be sorted
                sortorder = None

                if isinstance(key, Index):
                    key = np.asarray(key)

            new_codes = [level_codes[key] for level_codes in self.codes]

            return MultiIndex(
                levels=self.levels,
                codes=new_codes,
                names=self.names,
                sortorder=sortorder,
                verify_integrity=False,
            )

    @Appender(_index_shared_docs["take"] % _index_doc_kwargs)
    def take(self, indices, axis=0, allow_fill=True, fill_value=None, **kwargs):
        nv.validate_take(tuple(), kwargs)
        indices = ensure_platform_int(indices)
        taken = self._assert_take_fillable(
            self.codes,
            indices,
            allow_fill=allow_fill,
            fill_value=fill_value,
            na_value=-1,
        )
        return MultiIndex(
            levels=self.levels, codes=taken, names=self.names, verify_integrity=False
        )

    def _assert_take_fillable(
        self, values, indices, allow_fill=True, fill_value=None, na_value=None
    ):
        """ Internal method to handle NA filling of take """
        # only fill if we are passing a non-None fill_value
        if allow_fill and fill_value is not None:
            if (indices < -1).any():
                msg = (
                    "When allow_fill=True and fill_value is not None, "
                    "all indices must be >= -1"
                )
                raise ValueError(msg)
            taken = [lab.take(indices) for lab in self.codes]
            mask = indices == -1
            if mask.any():
                masked = []
                for new_label in taken:
                    label_values = new_label
                    label_values[mask] = na_value
                    masked.append(np.asarray(label_values))
                taken = masked
        else:
            taken = [lab.take(indices) for lab in self.codes]
        return taken

    def append(self, other):
        """
        Append a collection of Index options together

        Parameters
        ----------
        other : Index or list/tuple of indices

        Returns
        -------
        appended : Index
        """
        if not isinstance(other, (list, tuple)):
            other = [other]

        if all(
            (isinstance(o, MultiIndex) and o.nlevels >= self.nlevels) for o in other
        ):
            arrays = []
            for i in range(self.nlevels):
                label = self._get_level_values(i)
                appended = [o._get_level_values(i) for o in other]
                arrays.append(label.append(appended))
            return MultiIndex.from_arrays(arrays, names=self.names)

        to_concat = (self.values,) + tuple(k._values for k in other)
        new_tuples = np.concatenate(to_concat)

        # if all(isinstance(x, MultiIndex) for x in other):
        try:
            return MultiIndex.from_tuples(new_tuples, names=self.names)
        except (TypeError, IndexError):
            return Index(new_tuples)

    def argsort(self, *args, **kwargs) -> np.ndarray:
        return self.values.argsort(*args, **kwargs)

    @Appender(_index_shared_docs["repeat"] % _index_doc_kwargs)
    def repeat(self, repeats, axis=None):
        nv.validate_repeat(tuple(), dict(axis=axis))
        repeats = ensure_platform_int(repeats)
        return MultiIndex(
            levels=self.levels,
            codes=[
                level_codes.view(np.ndarray).astype(np.intp).repeat(repeats)
                for level_codes in self.codes
            ],
            names=self.names,
            sortorder=self.sortorder,
            verify_integrity=False,
        )

    def where(self, cond, other=None):
        raise NotImplementedError(".where is not supported for MultiIndex operations")

    def drop(self, codes, level=None, errors="raise"):
        """
        Make new MultiIndex with passed list of codes deleted

        Parameters
        ----------
        codes : array-like
            Must be a list of tuples
        level : int or level name, default None
        errors : str, default 'raise'

        Returns
        -------
        dropped : MultiIndex
        """
        if level is not None:
            return self._drop_from_level(codes, level, errors)

        if not isinstance(codes, (np.ndarray, Index)):
            try:
                codes = com.index_labels_to_array(codes, dtype=object)
            except ValueError:
                pass

        inds = []
        for level_codes in codes:
            try:
                loc = self.get_loc(level_codes)
                # get_loc returns either an integer, a slice, or a boolean
                # mask
                if isinstance(loc, int):
                    inds.append(loc)
                elif isinstance(loc, slice):
                    inds.extend(range(loc.start, loc.stop))
                elif com.is_bool_indexer(loc):
                    if self.lexsort_depth == 0:
                        warnings.warn(
                            "dropping on a non-lexsorted multi-index "
                            "without a level parameter may impact performance.",
                            PerformanceWarning,
                            stacklevel=3,
                        )
                    loc = loc.nonzero()[0]
                    inds.extend(loc)
                else:
                    msg = f"unsupported indexer of type {type(loc)}"
                    raise AssertionError(msg)
            except KeyError:
                if errors != "ignore":
                    raise

        return self.delete(inds)

    def _drop_from_level(self, codes, level, errors="raise"):
        codes = com.index_labels_to_array(codes)
        i = self._get_level_number(level)
        index = self.levels[i]
        values = index.get_indexer(codes)

        mask = ~algos.isin(self.codes[i], values)
        if mask.all() and errors != "ignore":
            raise KeyError(f"labels {codes} not found in level")

        return self[mask]

    def swaplevel(self, i=-2, j=-1):
        """
        Swap level i with level j.

        Calling this method does not change the ordering of the values.

        Parameters
        ----------
        i : int, str, default -2
            First level of index to be swapped. Can pass level name as string.
            Type of parameters can be mixed.
        j : int, str, default -1
            Second level of index to be swapped. Can pass level name as string.
            Type of parameters can be mixed.

        Returns
        -------
        MultiIndex
            A new MultiIndex.

        See Also
        --------
        Series.swaplevel : Swap levels i and j in a MultiIndex.
        Dataframe.swaplevel : Swap levels i and j in a MultiIndex on a
            particular axis.

        Examples
        --------
        >>> mi = pd.MultiIndex(levels=[['a', 'b'], ['bb', 'aa']],
        ...                    codes=[[0, 0, 1, 1], [0, 1, 0, 1]])
        >>> mi
        MultiIndex([('a', 'bb'),
                    ('a', 'aa'),
                    ('b', 'bb'),
                    ('b', 'aa')],
                   )
        >>> mi.swaplevel(0, 1)
        MultiIndex([('bb', 'a'),
                    ('aa', 'a'),
                    ('bb', 'b'),
                    ('aa', 'b')],
                   )
        """
        new_levels = list(self.levels)
        new_codes = list(self.codes)
        new_names = list(self.names)

        i = self._get_level_number(i)
        j = self._get_level_number(j)

        new_levels[i], new_levels[j] = new_levels[j], new_levels[i]
        new_codes[i], new_codes[j] = new_codes[j], new_codes[i]
        new_names[i], new_names[j] = new_names[j], new_names[i]

        return MultiIndex(
            levels=new_levels, codes=new_codes, names=new_names, verify_integrity=False
        )

    def reorder_levels(self, order):
        """
        Rearrange levels using input order. May not drop or duplicate levels.

        Parameters
        ----------

        Returns
        -------
        MultiIndex
        """
        order = [self._get_level_number(i) for i in order]
        if len(order) != self.nlevels:
            raise AssertionError(
                f"Length of order must be same as number of levels ({self.nlevels}), "
                f"got {len(order)}"
            )
        new_levels = [self.levels[i] for i in order]
        new_codes = [self.codes[i] for i in order]
        new_names = [self.names[i] for i in order]

        return MultiIndex(
            levels=new_levels, codes=new_codes, names=new_names, verify_integrity=False
        )

    def _get_codes_for_sorting(self):
        """
        we categorizing our codes by using the
        available categories (all, not just observed)
        excluding any missing ones (-1); this is in preparation
        for sorting, where we need to disambiguate that -1 is not
        a valid valid
        """

        def cats(level_codes):
            return np.arange(
                np.array(level_codes).max() + 1 if len(level_codes) else 0,
                dtype=level_codes.dtype,
            )

        return [
            Categorical.from_codes(level_codes, cats(level_codes), ordered=True)
            for level_codes in self.codes
        ]

    def sortlevel(self, level=0, ascending=True, sort_remaining=True):
        """
        Sort MultiIndex at the requested level. The result will respect the
        original ordering of the associated factor at that level.

        Parameters
        ----------
        level : list-like, int or str, default 0
            If a string is given, must be a name of the level.
            If list-like must be names or ints of levels.
        ascending : bool, default True
            False to sort in descending order.
            Can also be a list to specify a directed ordering.
        sort_remaining : sort by the remaining levels after level

        Returns
        -------
        sorted_index : pd.MultiIndex
            Resulting index.
        indexer : np.ndarray
            Indices of output values in original index.
        """
        if isinstance(level, (str, int)):
            level = [level]
        level = [self._get_level_number(lev) for lev in level]
        sortorder = None

        # we have a directed ordering via ascending
        if isinstance(ascending, list):
            if not len(level) == len(ascending):
                raise ValueError("level must have same length as ascending")

            indexer = lexsort_indexer(
                [self.codes[lev] for lev in level], orders=ascending
            )

        # level ordering
        else:

            codes = list(self.codes)
            shape = list(self.levshape)

            # partition codes and shape
            primary = tuple(codes[lev] for lev in level)
            primshp = tuple(shape[lev] for lev in level)

            # Reverse sorted to retain the order of
            # smaller indices that needs to be removed
            for lev in sorted(level, reverse=True):
                codes.pop(lev)
                shape.pop(lev)

            if sort_remaining:
                primary += primary + tuple(codes)
                primshp += primshp + tuple(shape)
            else:
                sortorder = level[0]

            indexer = indexer_from_factorized(primary, primshp, compress=False)

            if not ascending:
                indexer = indexer[::-1]

        indexer = ensure_platform_int(indexer)
        new_codes = [level_codes.take(indexer) for level_codes in self.codes]

        new_index = MultiIndex(
            codes=new_codes,
            levels=self.levels,
            names=self.names,
            sortorder=sortorder,
            verify_integrity=False,
        )

        return new_index, indexer

    def _convert_listlike_indexer(self, keyarr, kind=None):
        """
        Parameters
        ----------
        keyarr : list-like
            Indexer to convert.

        Returns
        -------
        tuple (indexer, keyarr)
            indexer is an ndarray or None if cannot convert
            keyarr are tuple-safe keys
        """
        indexer, keyarr = super()._convert_listlike_indexer(keyarr, kind=kind)

        # are we indexing a specific level
        if indexer is None and len(keyarr) and not isinstance(keyarr[0], tuple):
            level = 0
            _, indexer = self.reindex(keyarr, level=level)

            # take all
            if indexer is None:
                indexer = np.arange(len(self))

            check = self.levels[0].get_indexer(keyarr)
            mask = check == -1
            if mask.any():
                raise KeyError(f"{keyarr[mask]} not in index")

        return indexer, keyarr

    @Appender(_index_shared_docs["get_indexer"] % _index_doc_kwargs)
    def get_indexer(self, target, method=None, limit=None, tolerance=None):
        method = missing.clean_reindex_fill_method(method)
        target = ensure_index(target)

        # empty indexer
        if is_list_like(target) and not len(target):
            return ensure_platform_int(np.array([]))

        if not isinstance(target, MultiIndex):
            try:
                target = MultiIndex.from_tuples(target)
            except (TypeError, ValueError):

                # let's instead try with a straight Index
                if method is None:
                    return Index(self.values).get_indexer(
                        target, method=method, limit=limit, tolerance=tolerance
                    )

        if not self.is_unique:
            raise ValueError("Reindexing only valid with uniquely valued Index objects")

        if method == "pad" or method == "backfill":
            if tolerance is not None:
                raise NotImplementedError(
                    "tolerance not implemented yet for MultiIndex"
                )
            indexer = self._engine.get_indexer(target, method, limit)
        elif method == "nearest":
            raise NotImplementedError(
                "method='nearest' not implemented yet "
                "for MultiIndex; see GitHub issue 9365"
            )
        else:
            indexer = self._engine.get_indexer(target)

        return ensure_platform_int(indexer)

    @Appender(_index_shared_docs["get_indexer_non_unique"] % _index_doc_kwargs)
    def get_indexer_non_unique(self, target):
        return super().get_indexer_non_unique(target)

    def reindex(self, target, method=None, level=None, limit=None, tolerance=None):
        """
        Create index with target's values (move/add/delete values as necessary)

        Returns
        -------
        new_index : pd.MultiIndex
            Resulting index
        indexer : np.ndarray or None
            Indices of output values in original index.

        """
        # GH6552: preserve names when reindexing to non-named target
        # (i.e. neither Index nor Series).
        preserve_names = not hasattr(target, "names")

        if level is not None:
            if method is not None:
                raise TypeError("Fill method not supported if level passed")

            # GH7774: preserve dtype/tz if target is empty and not an Index.
            # target may be an iterator
            target = ibase._ensure_has_len(target)
            if len(target) == 0 and not isinstance(target, Index):
                idx = self.levels[level]
                attrs = idx._get_attributes_dict()
                attrs.pop("freq", None)  # don't preserve freq
                target = type(idx)._simple_new(np.empty(0, dtype=idx.dtype), **attrs)
            else:
                target = ensure_index(target)
            target, indexer, _ = self._join_level(
                target, level, how="right", return_indexers=True, keep_order=False
            )
        else:
            target = ensure_index(target)
            if self.equals(target):
                indexer = None
            else:
                if self.is_unique:
                    indexer = self.get_indexer(
                        target, method=method, limit=limit, tolerance=tolerance
                    )
                else:
                    raise ValueError("cannot handle a non-unique multi-index!")

        if not isinstance(target, MultiIndex):
            if indexer is None:
                target = self
            elif (indexer >= 0).all():
                target = self.take(indexer)
            else:
                # hopefully?
                target = MultiIndex.from_tuples(target)

        if (
            preserve_names
            and target.nlevels == self.nlevels
            and target.names != self.names
        ):
            target = target.copy(deep=False)
            target.names = self.names

        return target, indexer

    def get_slice_bound(
        self, label: Union[Hashable, Sequence[Hashable]], side: str, kind: str
    ) -> int:
        """
        For an ordered MultiIndex, compute slice bound
        that corresponds to given label.

        Returns leftmost (one-past-the-rightmost if `side=='right') position
        of given label.

        Parameters
        ----------
        label : object or tuple of objects
        side : {'left', 'right'}
        kind : {'loc', 'getitem'}

        Returns
        -------
        int
            Index of label.

        Notes
        -----
        This method only works if level 0 index of the MultiIndex is lexsorted.

        Examples
        --------
        >>> mi = pd.MultiIndex.from_arrays([list('abbc'), list('gefd')])

        Get the locations from the leftmost 'b' in the first level
        until the end of the multiindex:

        >>> mi.get_slice_bound('b', side="left", kind="loc")
        1

        Like above, but if you get the locations from the rightmost
        'b' in the first level and 'f' in the second level:

        >>> mi.get_slice_bound(('b','f'), side="right", kind="loc")
        3

        See Also
        --------
        MultiIndex.get_loc : Get location for a label or a tuple of labels.
        MultiIndex.get_locs : Get location for a label/slice/list/mask or a
                              sequence of such.
        """

        if not isinstance(label, tuple):
            label = (label,)
        return self._partial_tup_index(label, side=side)

    def slice_locs(self, start=None, end=None, step=None, kind=None):
        """
        For an ordered MultiIndex, compute the slice locations for input
        labels.

        The input labels can be tuples representing partial levels, e.g. for a
        MultiIndex with 3 levels, you can pass a single value (corresponding to
        the first level), or a 1-, 2-, or 3-tuple.

        Parameters
        ----------
        start : label or tuple, default None
            If None, defaults to the beginning
        end : label or tuple
            If None, defaults to the end
        step : int or None
            Slice step
        kind : string, optional, defaults None

        Returns
        -------
        (start, end) : (int, int)

        Notes
        -----
        This method only works if the MultiIndex is properly lexsorted. So,
        if only the first 2 levels of a 3-level MultiIndex are lexsorted,
        you can only pass two levels to ``.slice_locs``.

        Examples
        --------
        >>> mi = pd.MultiIndex.from_arrays([list('abbd'), list('deff')],
        ...                                names=['A', 'B'])

        Get the slice locations from the beginning of 'b' in the first level
        until the end of the multiindex:

        >>> mi.slice_locs(start='b')
        (1, 4)

        Like above, but stop at the end of 'b' in the first level and 'f' in
        the second level:

        >>> mi.slice_locs(start='b', end=('b', 'f'))
        (1, 3)

        See Also
        --------
        MultiIndex.get_loc : Get location for a label or a tuple of labels.
        MultiIndex.get_locs : Get location for a label/slice/list/mask or a
                              sequence of such.
        """
        # This function adds nothing to its parent implementation (the magic
        # happens in get_slice_bound method), but it adds meaningful doc.
        return super().slice_locs(start, end, step, kind=kind)

    def _partial_tup_index(self, tup, side="left"):
        if len(tup) > self.lexsort_depth:
            raise UnsortedIndexError(
                f"Key length ({len(tup)}) was greater than MultiIndex lexsort depth "
                f"({self.lexsort_depth})"
            )

        n = len(tup)
        start, end = 0, len(self)
        zipped = zip(tup, self.levels, self.codes)
        for k, (lab, lev, labs) in enumerate(zipped):
            section = labs[start:end]

            if lab not in lev and not isna(lab):
                if not lev.is_type_compatible(lib.infer_dtype([lab], skipna=False)):
                    raise TypeError(f"Level type mismatch: {lab}")

                # short circuit
                loc = lev.searchsorted(lab, side=side)
                if side == "right" and loc >= 0:
                    loc -= 1
                return start + section.searchsorted(loc, side=side)

            idx = self._get_loc_single_level_index(lev, lab)
            if k < n - 1:
                end = start + section.searchsorted(idx, side="right")
                start = start + section.searchsorted(idx, side="left")
            else:
                return start + section.searchsorted(idx, side=side)

    def _get_loc_single_level_index(self, level_index: Index, key: Hashable) -> int:
        """
        If key is NA value, location of index unify as -1.

        Parameters
        ----------
        level_index: Index
        key : label

        Returns
        -------
        loc : int
            If key is NA value, loc is -1
            Else, location of key in index.

        See Also
        --------
        Index.get_loc : The get_loc method for (single-level) index.
        """

        if is_scalar(key) and isna(key):
            return -1
        else:
            return level_index.get_loc(key)

    def get_loc(self, key, method=None):
        """
        Get location for a label or a tuple of labels as an integer, slice or
        boolean mask.

        Parameters
        ----------
        key : label or tuple of labels (one for each level)
        method : None

        Returns
        -------
        loc : int, slice object or boolean mask
            If the key is past the lexsort depth, the return may be a
            boolean mask array, otherwise it is always a slice or int.

        See Also
        --------
        Index.get_loc : The get_loc method for (single-level) index.
        MultiIndex.slice_locs : Get slice location given start label(s) and
                                end label(s).
        MultiIndex.get_locs : Get location for a label/slice/list/mask or a
                              sequence of such.

        Notes
        -----
        The key cannot be a slice, list of same-level labels, a boolean mask,
        or a sequence of such. If you want to use those, use
        :meth:`MultiIndex.get_locs` instead.

        Examples
        --------
        >>> mi = pd.MultiIndex.from_arrays([list('abb'), list('def')])

        >>> mi.get_loc('b')
        slice(1, 3, None)

        >>> mi.get_loc(('b', 'e'))
        1
        """
        if method is not None:
            raise NotImplementedError(
                "only the default get_loc method is "
                "currently supported for MultiIndex"
            )

        def _maybe_to_slice(loc):
            """convert integer indexer to boolean mask or slice if possible"""
            if not isinstance(loc, np.ndarray) or loc.dtype != "int64":
                return loc

            loc = lib.maybe_indices_to_slice(loc, len(self))
            if isinstance(loc, slice):
                return loc

            mask = np.empty(len(self), dtype="bool")
            mask.fill(False)
            mask[loc] = True
            return mask

        if not isinstance(key, (tuple, list)):
            # not including list here breaks some indexing, xref #30892
            loc = self._get_level_indexer(key, level=0)
            return _maybe_to_slice(loc)

        keylen = len(key)
        if self.nlevels < keylen:
            raise KeyError(
                f"Key length ({keylen}) exceeds index depth ({self.nlevels})"
            )

        if keylen == self.nlevels and self.is_unique:
            return self._engine.get_loc(key)

        # -- partial selection or non-unique index
        # break the key into 2 parts based on the lexsort_depth of the index;
        # the first part returns a continuous slice of the index; the 2nd part
        # needs linear search within the slice
        i = self.lexsort_depth
        lead_key, follow_key = key[:i], key[i:]
        start, stop = (
            self.slice_locs(lead_key, lead_key) if lead_key else (0, len(self))
        )

        if start == stop:
            raise KeyError(key)

        if not follow_key:
            return slice(start, stop)

        warnings.warn(
            "indexing past lexsort depth may impact performance.",
            PerformanceWarning,
            stacklevel=10,
        )

        loc = np.arange(start, stop, dtype="int64")

        for i, k in enumerate(follow_key, len(lead_key)):
            mask = self.codes[i][loc] == self._get_loc_single_level_index(
                self.levels[i], k
            )
            if not mask.all():
                loc = loc[mask]
            if not len(loc):
                raise KeyError(key)

        return _maybe_to_slice(loc) if len(loc) != stop - start else slice(start, stop)

    def get_loc_level(self, key, level=0, drop_level: bool = True):
        """
        Get both the location for the requested label(s) and the
        resulting sliced index.

        Parameters
        ----------
        key : label or sequence of labels
        level : int/level name or list thereof, optional
        drop_level : bool, default True
            If ``False``, the resulting index will not drop any level.

        Returns
        -------
        loc : A 2-tuple where the elements are:
              Element 0: int, slice object or boolean array
              Element 1: The resulting sliced multiindex/index. If the key
              contains all levels, this will be ``None``.

        See Also
        --------
        MultiIndex.get_loc  : Get location for a label or a tuple of labels.
        MultiIndex.get_locs : Get location for a label/slice/list/mask or a
                              sequence of such.

        Examples
        --------
        >>> mi = pd.MultiIndex.from_arrays([list('abb'), list('def')],
        ...                                names=['A', 'B'])

        >>> mi.get_loc_level('b')
        (slice(1, 3, None), Index(['e', 'f'], dtype='object', name='B'))

        >>> mi.get_loc_level('e', level='B')
        (array([False,  True, False], dtype=bool),
        Index(['b'], dtype='object', name='A'))

        >>> mi.get_loc_level(['b', 'e'])
        (1, None)
        """

        # different name to distinguish from maybe_droplevels
        def maybe_mi_droplevels(indexer, levels, drop_level: bool):
            if not drop_level:
                return self[indexer]
            # kludgearound
            orig_index = new_index = self[indexer]
            levels = [self._get_level_number(i) for i in levels]
            for i in sorted(levels, reverse=True):
                try:
                    new_index = new_index.droplevel(i)
                except ValueError:

                    # no dropping here
                    return orig_index
            return new_index

        if isinstance(level, (tuple, list)):
            if len(key) != len(level):
                raise AssertionError(
                    "Key for location must have same length as number of levels"
                )
            result = None
            for lev, k in zip(level, key):
                loc, new_index = self.get_loc_level(k, level=lev)
                if isinstance(loc, slice):
                    mask = np.zeros(len(self), dtype=bool)
                    mask[loc] = True
                    loc = mask

                result = loc if result is None else result & loc

            return result, maybe_mi_droplevels(result, level, drop_level)

        level = self._get_level_number(level)

        # kludge for #1796
        if isinstance(key, list):
            key = tuple(key)

        if isinstance(key, tuple) and level == 0:

            try:
                if key in self.levels[0]:
                    indexer = self._get_level_indexer(key, level=level)
                    new_index = maybe_mi_droplevels(indexer, [0], drop_level)
                    return indexer, new_index
            except (TypeError, InvalidIndexError):
                pass

            if not any(isinstance(k, slice) for k in key):

                # partial selection
                # optionally get indexer to avoid re-calculation
                def partial_selection(key, indexer=None):
                    if indexer is None:
                        indexer = self.get_loc(key)
                    ilevels = [
                        i for i in range(len(key)) if key[i] != slice(None, None)
                    ]
                    return indexer, maybe_mi_droplevels(indexer, ilevels, drop_level)

                if len(key) == self.nlevels and self.is_unique:
                    # Complete key in unique index -> standard get_loc
                    try:
                        return (self._engine.get_loc(key), None)
                    except KeyError as e:
                        raise KeyError(key) from e
                else:
                    return partial_selection(key)
            else:
                indexer = None
                for i, k in enumerate(key):
                    if not isinstance(k, slice):
                        k = self._get_level_indexer(k, level=i)
                        if isinstance(k, slice):
                            # everything
                            if k.start == 0 and k.stop == len(self):
                                k = slice(None, None)
                        else:
                            k_index = k

                    if isinstance(k, slice):
                        if k == slice(None, None):
                            continue
                        else:
                            raise TypeError(key)

                    if indexer is None:
                        indexer = k_index
                    else:  # pragma: no cover
                        indexer &= k_index
                if indexer is None:
                    indexer = slice(None, None)
                ilevels = [i for i in range(len(key)) if key[i] != slice(None, None)]
                return indexer, maybe_mi_droplevels(indexer, ilevels, drop_level)
        else:
            indexer = self._get_level_indexer(key, level=level)
            return indexer, maybe_mi_droplevels(indexer, [level], drop_level)

    def _get_level_indexer(self, key, level=0, indexer=None):
        # return an indexer, boolean array or a slice showing where the key is
        # in the totality of values
        # if the indexer is provided, then use this

        level_index = self.levels[level]
        level_codes = self.codes[level]

        def convert_indexer(start, stop, step, indexer=indexer, codes=level_codes):
            # given the inputs and the codes/indexer, compute an indexer set
            # if we have a provided indexer, then this need not consider
            # the entire labels set

            r = np.arange(start, stop, step)
            if indexer is not None and len(indexer) != len(codes):

                # we have an indexer which maps the locations in the labels
                # that we have already selected (and is not an indexer for the
                # entire set) otherwise this is wasteful so we only need to
                # examine locations that are in this set the only magic here is
                # that the result are the mappings to the set that we have
                # selected
                from pandas import Series

                mapper = Series(indexer)
                indexer = codes.take(ensure_platform_int(indexer))
                result = Series(Index(indexer).isin(r).nonzero()[0])
                m = result.map(mapper)._ndarray_values

            else:
                m = np.zeros(len(codes), dtype=bool)
                m[np.in1d(codes, r, assume_unique=Index(codes).is_unique)] = True

            return m

        if isinstance(key, slice):
            # handle a slice, returning a slice if we can
            # otherwise a boolean indexer

            try:
                if key.start is not None:
                    start = level_index.get_loc(key.start)
                else:
                    start = 0
                if key.stop is not None:
                    stop = level_index.get_loc(key.stop)
                else:
                    stop = len(level_index) - 1
                step = key.step
            except KeyError:

                # we have a partial slice (like looking up a partial date
                # string)
                start = stop = level_index.slice_indexer(
                    key.start, key.stop, key.step, kind="loc"
                )
                step = start.step

            if isinstance(start, slice) or isinstance(stop, slice):
                # we have a slice for start and/or stop
                # a partial date slicer on a DatetimeIndex generates a slice
                # note that the stop ALREADY includes the stopped point (if
                # it was a string sliced)
                start = getattr(start, "start", start)
                stop = getattr(stop, "stop", stop)
                return convert_indexer(start, stop, step)

            elif level > 0 or self.lexsort_depth == 0 or step is not None:
                # need to have like semantics here to right
                # searching as when we are using a slice
                # so include the stop+1 (so we include stop)
                return convert_indexer(start, stop + 1, step)
            else:
                # sorted, so can return slice object -> view
                i = level_codes.searchsorted(start, side="left")
                j = level_codes.searchsorted(stop, side="right")
                return slice(i, j, step)

        else:

            code = self._get_loc_single_level_index(level_index, key)

            if level > 0 or self.lexsort_depth == 0:
                # Desired level is not sorted
                locs = np.array(level_codes == code, dtype=bool, copy=False)
                if not locs.any():
                    # The label is present in self.levels[level] but unused:
                    raise KeyError(key)
                return locs

            i = level_codes.searchsorted(code, side="left")
            j = level_codes.searchsorted(code, side="right")
            if i == j:
                # The label is present in self.levels[level] but unused:
                raise KeyError(key)
            return slice(i, j)

    def get_locs(self, seq):
        """
        Get location for a sequence of labels.

        Parameters
        ----------
        seq : label, slice, list, mask or a sequence of such
           You should use one of the above for each level.
           If a level should not be used, set it to ``slice(None)``.

        Returns
        -------
        numpy.ndarray
            NumPy array of integers suitable for passing to iloc.

        See Also
        --------
        MultiIndex.get_loc : Get location for a label or a tuple of labels.
        MultiIndex.slice_locs : Get slice location given start label(s) and
                                end label(s).

        Examples
        --------
        >>> mi = pd.MultiIndex.from_arrays([list('abb'), list('def')])

        >>> mi.get_locs('b')  # doctest: +SKIP
        array([1, 2], dtype=int64)

        >>> mi.get_locs([slice(None), ['e', 'f']])  # doctest: +SKIP
        array([1, 2], dtype=int64)

        >>> mi.get_locs([[True, False, True], slice('e', 'f')])  # doctest: +SKIP
        array([2], dtype=int64)
        """
        from pandas.core.indexes.numeric import Int64Index

        # must be lexsorted to at least as many levels
        true_slices = [i for (i, s) in enumerate(com.is_true_slices(seq)) if s]
        if true_slices and true_slices[-1] >= self.lexsort_depth:
            raise UnsortedIndexError(
                "MultiIndex slicing requires the index to be lexsorted: slicing "
                f"on levels {true_slices}, lexsort depth {self.lexsort_depth}"
            )
        # indexer
        # this is the list of all values that we want to select
        n = len(self)
        indexer = None

        def _convert_to_indexer(r):
            # return an indexer
            if isinstance(r, slice):
                m = np.zeros(n, dtype=bool)
                m[r] = True
                r = m.nonzero()[0]
            elif com.is_bool_indexer(r):
                if len(r) != n:
                    raise ValueError(
                        "cannot index with a boolean indexer "
                        "that is not the same length as the "
                        "index"
                    )
                r = r.nonzero()[0]
            return Int64Index(r)

        def _update_indexer(idxr, indexer=indexer):
            if indexer is None:
                indexer = Index(np.arange(n))
            if idxr is None:
                return indexer
            return indexer & idxr

        for i, k in enumerate(seq):

            if com.is_bool_indexer(k):
                # a boolean indexer, must be the same length!
                k = np.asarray(k)
                indexer = _update_indexer(_convert_to_indexer(k), indexer=indexer)

            elif is_list_like(k):
                # a collection of labels to include from this level (these
                # are or'd)
                indexers = None
                for x in k:
                    try:
                        idxrs = _convert_to_indexer(
                            self._get_level_indexer(x, level=i, indexer=indexer)
                        )
                        indexers = idxrs if indexers is None else indexers | idxrs
                    except KeyError:

                        # ignore not founds
                        continue

                if indexers is not None:
                    indexer = _update_indexer(indexers, indexer=indexer)
                else:
                    # no matches we are done
                    return Int64Index([])._ndarray_values

            elif com.is_null_slice(k):
                # empty slice
                indexer = _update_indexer(None, indexer=indexer)

            elif isinstance(k, slice):

                # a slice, include BOTH of the labels
                indexer = _update_indexer(
                    _convert_to_indexer(
                        self._get_level_indexer(k, level=i, indexer=indexer)
                    ),
                    indexer=indexer,
                )
            else:
                # a single label
                indexer = _update_indexer(
                    _convert_to_indexer(
                        self.get_loc_level(k, level=i, drop_level=False)[0]
                    ),
                    indexer=indexer,
                )

        # empty indexer
        if indexer is None:
            return Int64Index([])._ndarray_values
        return indexer._ndarray_values

    def truncate(self, before=None, after=None):
        """
        Slice index between two labels / tuples, return new MultiIndex

        Parameters
        ----------
        before : label or tuple, can be partial. Default None
            None defaults to start
        after : label or tuple, can be partial. Default None
            None defaults to end

        Returns
        -------
        truncated : MultiIndex
        """
        if after and before and after < before:
            raise ValueError("after < before")

        i, j = self.levels[0].slice_locs(before, after)
        left, right = self.slice_locs(before, after)

        new_levels = list(self.levels)
        new_levels[0] = new_levels[0][i:j]

        new_codes = [level_codes[left:right] for level_codes in self.codes]
        new_codes[0] = new_codes[0] - i

        return MultiIndex(levels=new_levels, codes=new_codes, verify_integrity=False)

    def equals(self, other) -> bool:
        """
        Determines if two MultiIndex objects have the same labeling information
        (the levels themselves do not necessarily have to be the same)

        See Also
        --------
        equal_levels
        """
        if self.is_(other):
            return True

        if not isinstance(other, Index):
            return False

        if not isinstance(other, MultiIndex):
            # d-level MultiIndex can equal d-tuple Index
            if not is_object_dtype(other.dtype):
                if self.nlevels != other.nlevels:
                    return False

            other_vals = com.values_from_object(ensure_index(other))
            return array_equivalent(self._ndarray_values, other_vals)

        if self.nlevels != other.nlevels:
            return False

        if len(self) != len(other):
            return False

        for i in range(self.nlevels):
            self_codes = self.codes[i]
            self_codes = self_codes[self_codes != -1]
            self_values = algos.take_nd(
                np.asarray(self.levels[i]._values), self_codes, allow_fill=False
            )

            other_codes = other.codes[i]
            other_codes = other_codes[other_codes != -1]
            other_values = algos.take_nd(
                np.asarray(other.levels[i]._values), other_codes, allow_fill=False
            )

            # since we use NaT both datetime64 and timedelta64
            # we can have a situation where a level is typed say
            # timedelta64 in self (IOW it has other values than NaT)
            # but types datetime64 in other (where its all NaT)
            # but these are equivalent
            if len(self_values) == 0 and len(other_values) == 0:
                continue

            if not array_equivalent(self_values, other_values):
                return False

        return True

    def equal_levels(self, other) -> bool:
        """
        Return True if the levels of both MultiIndex objects are the same

        """
        if self.nlevels != other.nlevels:
            return False

        for i in range(self.nlevels):
            if not self.levels[i].equals(other.levels[i]):
                return False
        return True

    def union(self, other, sort=None):
        """
        Form the union of two MultiIndex objects

        Parameters
        ----------
        other : MultiIndex or array / Index of tuples
        sort : False or None, default None
            Whether to sort the resulting Index.

            * None : Sort the result, except when

              1. `self` and `other` are equal.
              2. `self` has length 0.
              3. Some values in `self` or `other` cannot be compared.
                 A RuntimeWarning is issued in this case.

            * False : do not sort the result.

            .. versionadded:: 0.24.0

            .. versionchanged:: 0.24.1

               Changed the default value from ``True`` to ``None``
               (without change in behaviour).

        Returns
        -------
        Index

        >>> index.union(index2)
        """
        self._validate_sort_keyword(sort)
        self._assert_can_do_setop(other)
        other, result_names = self._convert_can_do_setop(other)

        if len(other) == 0 or self.equals(other):
            return self

        # TODO: Index.union returns other when `len(self)` is 0.

        uniq_tuples = lib.fast_unique_multiple(
            [self._ndarray_values, other._ndarray_values], sort=sort
        )

        return MultiIndex.from_arrays(
            zip(*uniq_tuples), sortorder=0, names=result_names
        )

    def intersection(self, other, sort=False):
        """
        Form the intersection of two MultiIndex objects.

        Parameters
        ----------
        other : MultiIndex or array / Index of tuples
        sort : False or None, default False
            Sort the resulting MultiIndex if possible

            .. versionadded:: 0.24.0

            .. versionchanged:: 0.24.1

               Changed the default from ``True`` to ``False``, to match
               behaviour from before 0.24.0

        Returns
        -------
        Index
        """
        self._validate_sort_keyword(sort)
        self._assert_can_do_setop(other)
        other, result_names = self._convert_can_do_setop(other)

        if self.equals(other):
            return self

        self_tuples = self._ndarray_values
        other_tuples = other._ndarray_values
        uniq_tuples = set(self_tuples) & set(other_tuples)

        if sort is None:
            uniq_tuples = sorted(uniq_tuples)

        if len(uniq_tuples) == 0:
            return MultiIndex(
                levels=self.levels,
                codes=[[]] * self.nlevels,
                names=result_names,
                verify_integrity=False,
            )
        else:
            return MultiIndex.from_arrays(
                zip(*uniq_tuples), sortorder=0, names=result_names
            )

    def difference(self, other, sort=None):
        """
        Compute set difference of two MultiIndex objects

        Parameters
        ----------
        other : MultiIndex
        sort : False or None, default None
            Sort the resulting MultiIndex if possible

            .. versionadded:: 0.24.0

            .. versionchanged:: 0.24.1

               Changed the default value from ``True`` to ``None``
               (without change in behaviour).

        Returns
        -------
        diff : MultiIndex
        """
        self._validate_sort_keyword(sort)
        self._assert_can_do_setop(other)
        other, result_names = self._convert_can_do_setop(other)

        if len(other) == 0:
            return self

        if self.equals(other):
            return MultiIndex(
                levels=self.levels,
                codes=[[]] * self.nlevels,
                names=result_names,
                verify_integrity=False,
            )

        this = self._get_unique_index()

        indexer = this.get_indexer(other)
        indexer = indexer.take((indexer != -1).nonzero()[0])

        label_diff = np.setdiff1d(np.arange(this.size), indexer, assume_unique=True)
        difference = this.values.take(label_diff)
        if sort is None:
            difference = sorted(difference)

        if len(difference) == 0:
            return MultiIndex(
                levels=[[]] * self.nlevels,
                codes=[[]] * self.nlevels,
                names=result_names,
                verify_integrity=False,
            )
        else:
            return MultiIndex.from_tuples(difference, sortorder=0, names=result_names)

    @Appender(Index.astype.__doc__)
    def astype(self, dtype, copy=True):
        dtype = pandas_dtype(dtype)
        if is_categorical_dtype(dtype):
            msg = "> 1 ndim Categorical are not supported at this time"
            raise NotImplementedError(msg)
        elif not is_object_dtype(dtype):
            raise TypeError(
                f"Setting {type(self)} dtype to anything other "
                "than object is not supported"
            )
        elif copy is True:
            return self._shallow_copy()
        return self

    def _convert_can_do_setop(self, other):
        result_names = self.names

        if not hasattr(other, "names"):
            if len(other) == 0:
                other = MultiIndex(
                    levels=[[]] * self.nlevels,
                    codes=[[]] * self.nlevels,
                    verify_integrity=False,
                )
            else:
                msg = "other must be a MultiIndex or a list of tuples"
                try:
                    other = MultiIndex.from_tuples(other)
                except TypeError:
                    raise TypeError(msg)
        else:
            result_names = self.names if self.names == other.names else None
        return other, result_names

    def insert(self, loc: int, item):
        """
        Make new MultiIndex inserting new item at location

        Parameters
        ----------
        loc : int
        item : tuple
            Must be same length as number of levels in the MultiIndex

        Returns
        -------
        new_index : Index
        """
        # Pad the key with empty strings if lower levels of the key
        # aren't specified:
        if not isinstance(item, tuple):
            item = (item,) + ("",) * (self.nlevels - 1)
        elif len(item) != self.nlevels:
            raise ValueError("Item must have length equal to number of levels.")

        new_levels = []
        new_codes = []
        for k, level, level_codes in zip(item, self.levels, self.codes):
            if k not in level:
                # have to insert into level
                # must insert at end otherwise you have to recompute all the
                # other codes
                lev_loc = len(level)
                level = level.insert(lev_loc, k)
            else:
                lev_loc = level.get_loc(k)

            new_levels.append(level)
            new_codes.append(np.insert(ensure_int64(level_codes), loc, lev_loc))

        return MultiIndex(
            levels=new_levels, codes=new_codes, names=self.names, verify_integrity=False
        )

    def delete(self, loc):
        """
        Make new index with passed location deleted

        Returns
        -------
        new_index : MultiIndex
        """
        new_codes = [np.delete(level_codes, loc) for level_codes in self.codes]
        return MultiIndex(
            levels=self.levels,
            codes=new_codes,
            names=self.names,
            verify_integrity=False,
        )

    def _wrap_joined_index(self, joined, other):
        names = self.names if self.names == other.names else None
        return MultiIndex.from_tuples(joined, names=names)

    @Appender(Index.isin.__doc__)
    def isin(self, values, level=None):
        if level is None:
            values = MultiIndex.from_tuples(values, names=self.names).values
            return algos.isin(self.values, values)
        else:
            num = self._get_level_number(level)
            levs = self.get_level_values(num)

            if levs.size == 0:
                return np.zeros(len(levs), dtype=np.bool_)
            return levs.isin(values)


MultiIndex._add_numeric_methods_disabled()
MultiIndex._add_numeric_methods_add_sub_disabled()
MultiIndex._add_logical_methods_disabled()


def _sparsify(label_list, start: int = 0, sentinel=""):
    pivoted = list(zip(*label_list))
    k = len(label_list)

    result = pivoted[: start + 1]
    prev = pivoted[start]

    for cur in pivoted[start + 1 :]:
        sparse_cur = []

        for i, (p, t) in enumerate(zip(prev, cur)):
            if i == k - 1:
                sparse_cur.append(t)
                result.append(sparse_cur)
                break

            if p == t:
                sparse_cur.append(sentinel)
            else:
                sparse_cur.extend(cur[i:])
                result.append(sparse_cur)
                break

        prev = cur

    return list(zip(*result))


def _get_na_rep(dtype) -> str:
    return {np.datetime64: "NaT", np.timedelta64: "NaT"}.get(dtype, "NaN")


def maybe_droplevels(index, key):
    """
    Attempt to drop level or levels from the given index.

    Parameters
    ----------
    index: Index
    key : scalar or tuple

    Returns
    -------
    Index
    """
    # drop levels
    original_index = index
    if isinstance(key, tuple):
        for _ in key:
            try:
                index = index.droplevel(0)
            except ValueError:
                # we have dropped too much, so back out
                return original_index
    else:
        try:
            index = index.droplevel(0)
        except ValueError:
            pass

    return index


def _coerce_indexer_frozen(array_like, categories, copy: bool = False) -> np.ndarray:
    """
    Coerce the array_like indexer to the smallest integer dtype that can encode all
    of the given categories.

    Parameters
    ----------
    array_like : array-like
    categories : array-like
    copy : bool

    Returns
    -------
    np.ndarray
        Non-writeable.
    """
    array_like = coerce_indexer_dtype(array_like, categories)
    if copy:
        array_like = array_like.copy()
    array_like.flags.writeable = False
    return array_like
