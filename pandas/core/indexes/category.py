from __future__ import annotations

from typing import (
    TYPE_CHECKING,
    Any,
    Literal,
    Self,
    cast,
)

import numpy as np

from pandas._libs import index as libindex
from pandas.util._decorators import (
    cache_readonly,
    set_module,
)

from pandas.core.dtypes.common import is_scalar
from pandas.core.dtypes.dtypes import CategoricalDtype
from pandas.core.dtypes.missing import (
    is_valid_na_for_dtype,
)

from pandas.core.arrays.categorical import (
    Categorical,
    contains,
)
from pandas.core.construction import extract_array
from pandas.core.indexes.base import (
    Index,
    maybe_extract_name,
)
from pandas.core.indexes.extension import NDArrayBackedExtensionIndex

if TYPE_CHECKING:
    from collections.abc import (
        Callable,
        Hashable,
    )

    from pandas._typing import (
        ArrayLike,
        Axes,
        Dtype,
        DtypeObj,
        Level,
        ReindexMethod,
        npt,
    )

    from pandas import Series
    from pandas.core.indexes.multi import MultiIndex


@set_module("pandas")
class CategoricalIndex(NDArrayBackedExtensionIndex):
    """
    Index based on an underlying :class:`Categorical`.

    CategoricalIndex, like Categorical, can only take on a limited,
    and usually fixed, number of possible values (`categories`). Also,
    like Categorical, it might have an order, but numerical operations
    (additions, divisions, ...) are not possible.

    Parameters
    ----------
    data : array-like (1-dimensional)
        The values of the categorical. If `categories` are given, values not in
        `categories` will be replaced with NaN.
    categories : index-like, optional
        The categories for the categorical. Items need to be unique.
        If the categories are not given here (and also not in `dtype`), they
        will be inferred from the `data`.
    ordered : bool, optional
        Whether or not this categorical is treated as an ordered
        categorical. If not given here or in `dtype`, the resulting
        categorical will be unordered.
    dtype : CategoricalDtype or "category", optional
        If :class:`CategoricalDtype`, cannot be used together with
        `categories` or `ordered`.
    copy : bool, default False
        Make a copy of input ndarray.
    name : object, optional
        Name to be stored in the index.

    Attributes
    ----------
    codes
    categories
    ordered

    Methods
    -------
    rename_categories
    reorder_categories
    add_categories
    remove_categories
    remove_unused_categories
    set_categories
    as_ordered
    as_unordered
    map

    Raises
    ------
    ValueError
        If the categories do not validate.
    TypeError
        If an explicit ``ordered=True`` is given but no `categories` and the
        `values` are not sortable.

    See Also
    --------
    Index : The base pandas Index type.
    Categorical : A categorical array.
    CategoricalDtype : Type for categorical data.

    Notes
    -----
    See the `user guide
    <https://pandas.pydata.org/pandas-docs/stable/user_guide/advanced.html#categoricalindex>`__
    for more.

    Examples
    --------
    >>> pd.CategoricalIndex(["a", "b", "c", "a", "b", "c"])
    CategoricalIndex(['a', 'b', 'c', 'a', 'b', 'c'],
                     categories=['a', 'b', 'c'], ordered=False, dtype='category')

    ``CategoricalIndex`` can also be instantiated from a ``Categorical``:

    >>> c = pd.Categorical(["a", "b", "c", "a", "b", "c"])
    >>> pd.CategoricalIndex(c)
    CategoricalIndex(['a', 'b', 'c', 'a', 'b', 'c'],
                     categories=['a', 'b', 'c'], ordered=False, dtype='category')

    Ordered ``CategoricalIndex`` can have a min and max value.

    >>> ci = pd.CategoricalIndex(
    ...     ["a", "b", "c", "a", "b", "c"], ordered=True, categories=["c", "b", "a"]
    ... )
    >>> ci
    CategoricalIndex(['a', 'b', 'c', 'a', 'b', 'c'],
                     categories=['c', 'b', 'a'], ordered=True, dtype='category')
    >>> ci.min()
    'c'
    """

    _typ = "categoricalindex"
    _data_cls = Categorical

    @property
    def _can_hold_strings(self) -> bool:
        return self.categories._can_hold_strings

    @cache_readonly
    def _should_fallback_to_positional(self) -> bool:
        return self.categories._should_fallback_to_positional

    _data: Categorical
    _values: Categorical

    # --------------------------------------------------------------------
    # Categories/Codes/Ordered

    @property
    def codes(self) -> np.ndarray:
        """
        The category codes of this categorical index.

        Codes are an array of integers which are the positions of the actual
        values in the categories array. There is no setter.

        Returns
        -------
        ndarray[int]
            A non-writable view of the ``codes`` array.

        See Also
        --------
        Categorical.from_codes : Make a Categorical from codes.
        CategoricalIndex : An Index with an underlying ``Categorical``.

        Examples
        --------
        >>> ci = pd.CategoricalIndex(["a", "b", "c", "a", "b", "c"])
        >>> ci.codes
        array([0, 1, 2, 0, 1, 2], dtype=int8)

        >>> ci = pd.CategoricalIndex(["a", "c"], categories=["c", "b", "a"])
        >>> ci.codes
        array([2, 0], dtype=int8)
        """
        return self._data.codes

    @property
    def categories(self) -> Index:
        """
        The categories of this CategoricalIndex.

        These are the unique values the index may hold, in category order.
        Use the ``*_categories`` methods to return an index with changed
        categories.

        See Also
        --------
        CategoricalIndex.rename_categories : Rename categories.
        CategoricalIndex.reorder_categories : Reorder categories.
        CategoricalIndex.add_categories : Add new categories.
        CategoricalIndex.remove_categories : Remove the specified categories.
        CategoricalIndex.remove_unused_categories : Remove categories which
            are not used.
        CategoricalIndex.set_categories : Set the categories to the specified ones.

        Examples
        --------
        >>> ci = pd.CategoricalIndex(["a", "c", "b", "a", "c", "b"])
        >>> ci.categories
        Index(['a', 'b', 'c'], dtype='str')

        >>> ci = pd.CategoricalIndex(["a", "c"], categories=["c", "b", "a"])
        >>> ci.categories
        Index(['c', 'b', 'a'], dtype='str')
        """
        return self._data.categories

    @property
    def ordered(self) -> bool | None:
        """
        Whether the categories have an ordered relationship.

        This property returns True if the categories are ordered, meaning
        they have a meaningful order that allows comparison operations.

        See Also
        --------
        CategoricalIndex.as_ordered : Set the CategoricalIndex to be ordered.
        CategoricalIndex.as_unordered : Set the CategoricalIndex to be unordered.

        Examples
        --------
        >>> ci = pd.CategoricalIndex(["a", "b"], ordered=True)
        >>> ci.ordered
        True

        >>> ci = pd.CategoricalIndex(["a", "b"], ordered=False)
        >>> ci.ordered
        False
        """
        return self._data.ordered

    def _wrap_categorical_result(self, result: Categorical) -> Self:
        return type(self)._simple_new(result, name=self.name)

    def rename_categories(self, new_categories: Any) -> Self:
        """
        Rename categories.

        This method is commonly used to re-label or adjust the
        category names in categorical data without changing the
        underlying data. It is useful in situations where you want
        to modify the labels used for clarity, consistency,
        or readability.

        Parameters
        ----------
        new_categories : list-like, dict-like or callable

            New categories which will replace old categories.

            * list-like: all items must be unique and the number of items in
              the new categories must match the existing number of categories.

            * dict-like: specifies a mapping from
              old categories to new. Categories not contained in the mapping
              are passed through and extra categories in the mapping are
              ignored.

            * callable : a callable that is called on all items in the old
              categories and whose return values comprise the new categories.

        Returns
        -------
        CategoricalIndex
            CategoricalIndex with renamed categories.

        Raises
        ------
        ValueError
            If new categories are list-like and do not have the same number of
            items as the current categories or do not validate as categories

        See Also
        --------
        CategoricalIndex.reorder_categories : Reorder categories.
        CategoricalIndex.add_categories : Add new categories.
        CategoricalIndex.remove_categories : Remove the specified categories.
        CategoricalIndex.remove_unused_categories : Remove categories which
            are not used.
        CategoricalIndex.set_categories : Set the categories to the specified
            ones.

        Examples
        --------
        >>> ci = pd.CategoricalIndex(["a", "a", "b"])
        >>> ci.rename_categories([0, 1])
        CategoricalIndex([0, 0, 1], categories=[0, 1], ordered=False, dtype='category')

        For dict-like ``new_categories``, extra keys are ignored and
        categories not in the dictionary are passed through

        >>> ci.rename_categories({"a": "A", "c": "C"})
        CategoricalIndex(['A', 'A', 'b'], categories=['A', 'b'], ordered=False,
                         dtype='category')

        You may also provide a callable to create the new categories

        >>> ci.rename_categories(lambda x: x.upper())
        CategoricalIndex(['A', 'A', 'B'], categories=['A', 'B'], ordered=False,
                         dtype='category')
        """
        result = self._data.rename_categories(new_categories)
        return self._wrap_categorical_result(result)

    def reorder_categories(
        self, new_categories: Axes, ordered: bool | None = None
    ) -> Self:
        """
        Reorder categories as specified in new_categories.

        ``new_categories`` need to include all old categories and no new category
        items.

        Parameters
        ----------
        new_categories : Index-like
           The categories in new order.
        ordered : bool, optional
           Whether or not the categorical is treated as an ordered categorical.
           If not given, do not change the ordered information.

        Returns
        -------
        CategoricalIndex
            CategoricalIndex with reordered categories.

        Raises
        ------
        ValueError
            If the new categories do not contain all old category items or any
            new ones

        See Also
        --------
        CategoricalIndex.rename_categories : Rename categories.
        CategoricalIndex.add_categories : Add new categories.
        CategoricalIndex.remove_categories : Remove the specified categories.
        CategoricalIndex.remove_unused_categories : Remove categories which
            are not used.
        CategoricalIndex.set_categories : Set the categories to the specified
            ones.

        Examples
        --------
        >>> ci = pd.CategoricalIndex(["a", "b", "c", "a"])
        >>> ci
        CategoricalIndex(['a', 'b', 'c', 'a'], categories=['a', 'b', 'c'],
                         ordered=False, dtype='category')
        >>> ci.reorder_categories(["c", "b", "a"], ordered=True)
        CategoricalIndex(['a', 'b', 'c', 'a'], categories=['c', 'b', 'a'],
                         ordered=True, dtype='category')
        """
        result = self._data.reorder_categories(new_categories, ordered=ordered)
        return self._wrap_categorical_result(result)

    def add_categories(self, new_categories: Any) -> Self:
        """
        Add new categories.

        `new_categories` will be included at the last/highest place in the
        categories and will be unused directly after this call.

        Parameters
        ----------
        new_categories : category or list-like of category
            The new categories to be included.

        Returns
        -------
        CategoricalIndex
            CategoricalIndex with new categories added.

        Raises
        ------
        ValueError
            If the new categories include old categories or do not validate as
            categories

        See Also
        --------
        CategoricalIndex.rename_categories : Rename categories.
        CategoricalIndex.reorder_categories : Reorder categories.
        CategoricalIndex.remove_categories : Remove the specified categories.
        CategoricalIndex.remove_unused_categories : Remove categories which
            are not used.
        CategoricalIndex.set_categories : Set the categories to the specified
            ones.

        Examples
        --------
        >>> ci = pd.CategoricalIndex(["c", "b", "c"])
        >>> ci
        CategoricalIndex(['c', 'b', 'c'], categories=['b', 'c'], ordered=False,
                         dtype='category')

        >>> ci.add_categories(["d", "a"])
        CategoricalIndex(['c', 'b', 'c'], categories=['b', 'c', 'd', 'a'],
                         ordered=False, dtype='category')
        """
        result = self._data.add_categories(new_categories)
        return self._wrap_categorical_result(result)

    def remove_categories(self, removals: Any) -> Self:
        """
        Remove the specified categories.

        The ``removals`` argument must be a subset of the current categories.
        Any values that were part of the removed categories will be set to NaN.

        Parameters
        ----------
        removals : category or list of categories
           The categories which should be removed.

        Returns
        -------
        CategoricalIndex
            CategoricalIndex with removed categories.

        Raises
        ------
        ValueError
            If the removals are not contained in the categories

        See Also
        --------
        CategoricalIndex.rename_categories : Rename categories.
        CategoricalIndex.reorder_categories : Reorder categories.
        CategoricalIndex.add_categories : Add new categories.
        CategoricalIndex.remove_unused_categories : Remove categories which
            are not used.
        CategoricalIndex.set_categories : Set the categories to the specified
            ones.

        Examples
        --------
        >>> ci = pd.CategoricalIndex(["a", "c", "b", "c", "d"])
        >>> ci
        CategoricalIndex(['a', 'c', 'b', 'c', 'd'], categories=['a', 'b', 'c', 'd'],
                         ordered=False, dtype='category')

        >>> ci.remove_categories(["d", "a"])
        CategoricalIndex([nan, 'c', 'b', 'c', nan], categories=['b', 'c'],
                         ordered=False, dtype='category')
        """
        result = self._data.remove_categories(removals)
        return self._wrap_categorical_result(result)

    def remove_unused_categories(self) -> Self:
        """
        Remove categories which are not used.

        This method is useful when working with datasets
        that undergo dynamic changes where categories may no longer be
        relevant, allowing to maintain a clean, efficient data structure.

        Returns
        -------
        CategoricalIndex
            CategoricalIndex with unused categories dropped.

        See Also
        --------
        CategoricalIndex.rename_categories : Rename categories.
        CategoricalIndex.reorder_categories : Reorder categories.
        CategoricalIndex.add_categories : Add new categories.
        CategoricalIndex.remove_categories : Remove the specified categories.
        CategoricalIndex.set_categories : Set the categories to the specified
            ones.

        Examples
        --------
        >>> ci = pd.CategoricalIndex(["a", "c", "b", "c", "d"])
        >>> ci
        CategoricalIndex(['a', 'c', 'b', 'c', 'd'], categories=['a', 'b', 'c', 'd'],
                         ordered=False, dtype='category')

        >>> ci = ci[[0, 1, 0, 1, 1]]
        >>> ci
        CategoricalIndex(['a', 'c', 'a', 'c', 'c'], categories=['a', 'b', 'c', 'd'],
                         ordered=False, dtype='category')

        >>> ci.remove_unused_categories()
        CategoricalIndex(['a', 'c', 'a', 'c', 'c'], categories=['a', 'c'],
                         ordered=False, dtype='category')
        """
        result = self._data.remove_unused_categories()
        return self._wrap_categorical_result(result)

    def set_categories(
        self,
        new_categories: Axes,
        ordered: bool | None = None,
        rename: bool = False,
    ) -> Self:
        """
        Set the categories to the specified new categories.

        ``new_categories`` can include new categories (which will result in
        unused categories) or remove old categories (which results in values
        set to ``NaN``). If ``rename=True``, the categories will simply be renamed
        (less or more items than in old categories will result in values set to
        ``NaN`` or in unused categories respectively).

        This method can be used to perform more than one action of adding,
        removing, and reordering simultaneously and is therefore faster than
        performing the individual steps via the more specialised methods.

        On the other hand this method does not check whether the old categories
        are included in the new categories on a reorder, which can result in
        surprising changes.

        Parameters
        ----------
        new_categories : Index-like
           The categories in new order.
        ordered : bool, default None
           Whether or not the categorical is treated as an ordered categorical.
           If not given, do not change the ordered information.
        rename : bool, default False
           Whether or not the new_categories should be considered as a rename
           of the old categories or as reordered categories.

        Returns
        -------
        CategoricalIndex
            CategoricalIndex with the new categories.

        Raises
        ------
        ValueError
            If new_categories does not validate as categories

        See Also
        --------
        CategoricalIndex.rename_categories : Rename categories.
        CategoricalIndex.reorder_categories : Reorder categories.
        CategoricalIndex.add_categories : Add new categories.
        CategoricalIndex.remove_categories : Remove the specified categories.
        CategoricalIndex.remove_unused_categories : Remove categories which
            are not used.

        Examples
        --------
        >>> ci = pd.CategoricalIndex(
        ...     ["a", "b", "c", None], categories=["a", "b", "c"], ordered=True
        ... )
        >>> ci
        CategoricalIndex(['a', 'b', 'c', nan], categories=['a', 'b', 'c'],
                         ordered=True, dtype='category')

        >>> ci.set_categories(["A", "b", "c"])
        CategoricalIndex([nan, 'b', 'c', nan], categories=['A', 'b', 'c'],
                         ordered=True, dtype='category')
        >>> ci.set_categories(["A", "b", "c"], rename=True)
        CategoricalIndex(['A', 'b', 'c', nan], categories=['A', 'b', 'c'],
                         ordered=True, dtype='category')
        """
        result = self._data.set_categories(
            new_categories, ordered=ordered, rename=rename
        )
        return self._wrap_categorical_result(result)

    def as_ordered(self) -> Self:
        """
        Set the CategoricalIndex to be ordered.

        This method returns a new CategoricalIndex with the ordered attribute
        set to True, enabling comparison operations between categories.

        Returns
        -------
        CategoricalIndex
            Ordered CategoricalIndex.

        See Also
        --------
        CategoricalIndex.as_unordered : Set the CategoricalIndex to be unordered.

        Examples
        --------
        >>> ci = pd.CategoricalIndex(["a", "b", "c", "a"])
        >>> ci.ordered
        False
        >>> ci = ci.as_ordered()
        >>> ci.ordered
        True
        """
        result = self._data.as_ordered()
        return self._wrap_categorical_result(result)

    def as_unordered(self) -> Self:
        """
        Set the CategoricalIndex to be unordered.

        This method returns a new CategoricalIndex with the ordered attribute
        set to False, disabling comparison operations between categories.

        Returns
        -------
        CategoricalIndex
            Unordered CategoricalIndex.

        See Also
        --------
        CategoricalIndex.as_ordered : Set the CategoricalIndex to be ordered.

        Examples
        --------
        >>> ci = pd.CategoricalIndex(["a", "b", "c", "a"], ordered=True)
        >>> ci.ordered
        True
        >>> ci = ci.as_unordered()
        >>> ci.ordered
        False
        """
        result = self._data.as_unordered()
        return self._wrap_categorical_result(result)

    def min(self, *, skipna: bool = True, **kwargs: Any) -> Any:  # type: ignore[override]
        """
        Return the minimum value of the CategoricalIndex.

        Only an ordered CategoricalIndex has a minimum.

        Parameters
        ----------
        skipna : bool, default True
            Exclude NA/null values when showing the result.
        **kwargs
            Additional keyword arguments passed through to the reduction.

        Returns
        -------
        scalar
            The minimum of this CategoricalIndex, NA value if empty.

        Raises
        ------
        TypeError
            If the CategoricalIndex is not ordered.

        See Also
        --------
        CategoricalIndex.max : Return the maximum value of the CategoricalIndex.

        Examples
        --------
        >>> ci = pd.CategoricalIndex(["a", "b", "c", "a"], ordered=True)
        >>> ci.min()
        'a'
        """
        return self._data.min(skipna=skipna, **kwargs)

    def max(self, *, skipna: bool = True, **kwargs: Any) -> Any:  # type: ignore[override]
        """
        Return the maximum value of the CategoricalIndex.

        Only an ordered CategoricalIndex has a maximum.

        Parameters
        ----------
        skipna : bool, default True
            Exclude NA/null values when showing the result.
        **kwargs
            Additional keyword arguments passed through to the reduction.

        Returns
        -------
        scalar
            The maximum of this CategoricalIndex, NA value if empty.

        Raises
        ------
        TypeError
            If the CategoricalIndex is not ordered.

        See Also
        --------
        CategoricalIndex.min : Return the minimum value of the CategoricalIndex.

        Examples
        --------
        >>> ci = pd.CategoricalIndex(["a", "b", "c", "a"], ordered=True)
        >>> ci.max()
        'c'
        """
        return self._data.max(skipna=skipna, **kwargs)

    def _reverse_indexer(self) -> dict[Hashable, npt.NDArray[np.intp]]:
        """See Categorical._reverse_indexer."""
        return self._data._reverse_indexer()

    @property
    def _engine_type(self) -> type[libindex.IndexEngine]:
        # self.codes can have dtype int8, int16, int32 or int64, so we need
        # to return the corresponding engine type (libindex.Int8Engine, etc.).
        return {
            np.int8: libindex.Int8Engine,
            np.int16: libindex.Int16Engine,
            np.int32: libindex.Int32Engine,
            np.int64: libindex.Int64Engine,
        }[self.codes.dtype.type]

    # --------------------------------------------------------------------
    # Constructors

    def __new__(
        cls,
        data: Axes | None = None,
        categories: Axes | None = None,
        ordered: bool | None = None,
        dtype: Dtype | None = None,
        copy: bool = False,
        name: Hashable | None = None,
    ) -> Self:
        name = maybe_extract_name(name, data, cls)

        if is_scalar(data):
            # GH#38944 include None here, which pre-2.0 subbed in []
            cls._raise_scalar_data_error(data)

        data = Categorical(
            data, categories=categories, ordered=ordered, dtype=dtype, copy=copy
        )

        return cls._simple_new(data, name=name)

    # --------------------------------------------------------------------

    def _is_dtype_compat(self, other: Index) -> Categorical:
        """
        *this is an internal non-public method*

        provide a comparison between the dtype of self and other (coercing if
        needed)

        Parameters
        ----------
        other : Index

        Returns
        -------
        Categorical

        Raises
        ------
        TypeError if the dtypes are not compatible
        """
        if isinstance(other.dtype, CategoricalDtype):
            cat = extract_array(other)
            cat = cast("Categorical", cat)
            if not cat._categories_match_up_to_permutation(self._values):
                raise TypeError(
                    "categories must match existing categories when appending"
                )

        elif other._is_multi:
            # preempt raising NotImplementedError in isna call
            raise TypeError("MultiIndex is not dtype-compatible with CategoricalIndex")
        else:
            values = other

            codes = self.categories.get_indexer(values)
            if ((codes == -1) & ~values.isna()).any():
                # GH#37667 see test_equals_non_category
                raise TypeError(
                    "categories must match existing categories when appending"
                )
            cat = Categorical(other, dtype=self.dtype)
            other = CategoricalIndex(cat)
            if not other.isin(values).all():
                raise TypeError(
                    "cannot append a non-category item to a CategoricalIndex"
                )
            cat = other._values

        return cat

    def equals(self, other: object) -> bool:
        """
        Determine if two CategoricalIndex objects contain the same elements.

        The order and orderedness of elements matters. The categories matter,
        but the order of the categories matters only when ``ordered=True``.

        Parameters
        ----------
        other : object
            The CategoricalIndex object to compare with.

        Returns
        -------
        bool
            ``True`` if two :class:`pandas.CategoricalIndex` objects have equal
            elements, ``False`` otherwise.

        See Also
        --------
        Categorical.equals : Returns True if categorical arrays are equal.

        Examples
        --------
        >>> ci = pd.CategoricalIndex(["a", "b", "c", "a", "b", "c"])
        >>> ci2 = pd.CategoricalIndex(pd.Categorical(["a", "b", "c", "a", "b", "c"]))
        >>> ci.equals(ci2)
        True

        The order of elements matters.

        >>> ci3 = pd.CategoricalIndex(["c", "b", "a", "a", "b", "c"])
        >>> ci.equals(ci3)
        False

        The orderedness also matters.

        >>> ci4 = ci.as_ordered()
        >>> ci.equals(ci4)
        False

        The categories matter, but the order of the categories matters only when
        ``ordered=True``.

        >>> ci5 = ci.set_categories(["a", "b", "c", "d"])
        >>> ci.equals(ci5)
        False

        >>> ci6 = ci.set_categories(["b", "c", "a"])
        >>> ci.equals(ci6)
        True
        >>> ci_ordered = pd.CategoricalIndex(
        ...     ["a", "b", "c", "a", "b", "c"], ordered=True
        ... )
        >>> ci2_ordered = ci_ordered.set_categories(["b", "c", "a"])
        >>> ci_ordered.equals(ci2_ordered)
        False
        """
        if self.is_(other):  # type: ignore[arg-type]
            return True

        if not isinstance(other, Index):
            return False

        try:
            other = self._is_dtype_compat(other)
        except (TypeError, ValueError):
            return False

        return self._data.equals(other)

    # --------------------------------------------------------------------
    # Rendering Methods

    def _formatter_func(self, val: Hashable) -> str:
        return self.categories._formatter_func(val)

    def _format_attrs(self) -> list[tuple[str, str | int | bool | None]]:
        """
        Return a list of tuples of the (attr,formatted_value)
        """
        attrs: list[tuple[str, str | int | bool | None]]

        attrs = [
            (
                "categories",
                f"[{', '.join(self._data._repr_categories())}]",
            ),
            ("ordered", self.ordered),
        ]
        extra = super()._format_attrs()
        return attrs + extra

    # --------------------------------------------------------------------

    @property
    def inferred_type(self) -> str:
        return "categorical"

    def __contains__(self, key: Any) -> bool:
        """
        Return a boolean indicating whether the provided key is in the index.

        Parameters
        ----------
        key : label
            The key to check if it is present in the index.

        Returns
        -------
        bool
            Whether the key search is in the index.

        Raises
        ------
        TypeError
            If the key is not hashable.

        See Also
        --------
        Index.isin : Returns an ndarray of boolean dtype indicating whether the
            list-like key is in the index.

        Examples
        --------
        >>> idx = pd.Index([1, 2, 3, 4])
        >>> idx
        Index([1, 2, 3, 4], dtype='int64')

        >>> 2 in idx
        True
        >>> 6 in idx
        False
        """
        # if key is a NaN, check if any NaN is in self.
        if is_valid_na_for_dtype(key, self.categories.dtype):
            return self.hasnans
        if self.categories._typ == "rangeindex":
            container: Index | libindex.IndexEngine | libindex.ExtensionEngine = (
                self.categories
            )
        else:
            container = self._engine
        return contains(self, key, container=container)

    def reindex(
        self,
        target: Axes,
        method: ReindexMethod | None = None,
        level: Level | None = None,
        limit: int | None = None,
        tolerance: float | None = None,
    ) -> tuple[Index, npt.NDArray[np.intp] | None]:
        """
        Create index with target's values (move/add/delete values as necessary)

        Returns
        -------
        new_index : pd.Index
            Resulting index
        indexer : np.ndarray[np.intp] or None
            Indices of output values in original index

        """
        if method is not None:
            raise NotImplementedError(
                "argument method is not implemented for CategoricalIndex.reindex"
            )
        if level is not None:
            raise NotImplementedError(
                "argument level is not implemented for CategoricalIndex.reindex"
            )
        if limit is not None:
            raise NotImplementedError(
                "argument limit is not implemented for CategoricalIndex.reindex"
            )
        return super().reindex(target)

    # --------------------------------------------------------------------
    # Indexing Methods

    def _maybe_cast_indexer(self, key: Hashable) -> int:
        # GH#41933: we have to do this instead of self._data._validate_scalar
        #  because this will correctly get partial-indexing on Interval categories
        try:
            return self._data._unbox_scalar(key)
        except KeyError:
            if is_valid_na_for_dtype(key, self.categories.dtype):
                return -1
            raise

    def _maybe_cast_listlike_indexer(self, values: Axes) -> CategoricalIndex:
        if isinstance(values, CategoricalIndex):
            values = values._data
        if isinstance(values, Categorical):
            # Indexing on codes is more efficient if categories are the same,
            #  so we can apply some optimizations based on the degree of
            #  dtype-matching.
            cat = self._data._encode_with_my_categories(values)
            codes = cat._codes
        else:
            codes = self.categories.get_indexer(values)
            codes = codes.astype(self.codes.dtype, copy=False)
            cat = self._data._from_backing_data(codes)
        return type(self)._simple_new(cat)

    # --------------------------------------------------------------------

    def _is_comparable_dtype(self, dtype: DtypeObj) -> bool:
        return self.categories._is_comparable_dtype(dtype)

    def _intersection(
        self, other: Index, sort: bool = False
    ) -> Index | ArrayLike | MultiIndex:
        # Reached only via Index.intersection after dtype reconciliation, so
        # other is necessarily a CategoricalIndex with matching dtype.
        # For unordered CategoricalIndex, libjoin operates on integer codes.
        # When two indexes have the same categories in different order, matching
        # codes map to different values, giving incorrect results (GH#55335).
        # Reorder other's categories to match self so the codes align.
        other = cast("CategoricalIndex", other)
        if not self.ordered and not self.categories.equals(other.categories):
            reordered = other._data.reorder_categories(self.categories)
            other = other._shallow_copy(reordered)
        return super()._intersection(other, sort=sort)

    def _union(self, other: Index, sort: bool | None) -> Index | ArrayLike | MultiIndex:
        # See _intersection for explanation of GH#55335.
        other = cast("CategoricalIndex", other)
        if not self.ordered and not self.categories.equals(other.categories):
            reordered = other._data.reorder_categories(self.categories)
            other = other._shallow_copy(reordered)
        return super()._union(other, sort)

    def map(
        self,
        mapper: Callable[..., Any] | dict[Hashable, Any] | Series,
        na_action: Literal["ignore"] | None = None,
    ) -> Index:
        """
        Map values using input an input mapping or function.

        Maps the values (their categories, not the codes) of the index to new
        categories. If the mapping correspondence is one-to-one the result is a
        :class:`~pandas.CategoricalIndex` which has the same order property as
        the original, otherwise an :class:`~pandas.Index` is returned.

        If a `dict` or :class:`~pandas.Series` is used any unmapped category is
        mapped to `NaN`. Note that if this happens an :class:`~pandas.Index`
        will be returned.

        Parameters
        ----------
        mapper : function, dict, or Series
            Mapping correspondence.
        na_action : {None, 'ignore'}, default 'ignore'
            If 'ignore', propagate NaN values, without passing them to
            the mapping correspondence.

        Returns
        -------
        pandas.CategoricalIndex or pandas.Index
            Mapped index.

        See Also
        --------
        Index.map : Apply a mapping correspondence on an
            :class:`~pandas.Index`.
        Series.map : Apply a mapping correspondence on a
            :class:`~pandas.Series`.
        Series.apply : Apply more complex functions on a
            :class:`~pandas.Series`.

        Examples
        --------
        >>> idx = pd.CategoricalIndex(["a", "b", "c"])
        >>> idx
        CategoricalIndex(['a', 'b', 'c'], categories=['a', 'b', 'c'],
                          ordered=False, dtype='category')
        >>> idx.map(lambda x: x.upper())
        CategoricalIndex(['A', 'B', 'C'], categories=['A', 'B', 'C'],
                         ordered=False, dtype='category')
        >>> idx.map({"a": "first", "b": "second", "c": "third"})
        CategoricalIndex(['first', 'second', 'third'], categories=['first',
                         'second', 'third'], ordered=False, dtype='category')

        If the mapping is one-to-one the ordering of the categories is
        preserved:

        >>> idx = pd.CategoricalIndex(["a", "b", "c"], ordered=True)
        >>> idx
        CategoricalIndex(['a', 'b', 'c'], categories=['a', 'b', 'c'],
                         ordered=True, dtype='category')
        >>> idx.map({"a": 3, "b": 2, "c": 1})
        CategoricalIndex([3, 2, 1], categories=[3, 2, 1], ordered=True,
                         dtype='category')

        If the mapping is not one-to-one an :class:`~pandas.Index` is returned:

        >>> idx.map({"a": "first", "b": "second", "c": "first"})
        Index(['first', 'second', 'first'], dtype='str')

        If a `dict` is used, all unmapped categories are mapped to `NaN` and
        the result is an :class:`~pandas.Index`:

        >>> idx.map({"a": "first", "b": "second"})
        Index(['first', 'second', NaN], dtype='str')
        """
        mapped = self._values.map(mapper, na_action=na_action)
        return Index(mapped, name=self.name, copy=False)
