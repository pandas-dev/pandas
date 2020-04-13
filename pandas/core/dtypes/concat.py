"""
Utility functions related to concat.
"""

import numpy as np

from pandas.core.dtypes.common import (
    is_bool_dtype,
    is_categorical_dtype,
    is_datetime64_dtype,
    is_datetime64tz_dtype,
    is_dtype_equal,
    is_extension_array_dtype,
    is_object_dtype,
    is_sparse,
    is_timedelta64_dtype,
)
from pandas.core.dtypes.generic import ABCCategoricalIndex, ABCRangeIndex, ABCSeries


def get_dtype_kinds(l):
    """
    Parameters
    ----------
    l : list of arrays

    Returns
    -------
    a set of kinds that exist in this list of arrays
    """
    typs = set()
    for arr in l:

        dtype = arr.dtype
        if is_categorical_dtype(dtype):
            typ = "category"
        elif is_sparse(arr):
            typ = "sparse"
        elif isinstance(arr, ABCRangeIndex):
            typ = "range"
        elif is_datetime64tz_dtype(arr):
            # if to_concat contains different tz,
            # the result must be object dtype
            typ = str(arr.dtype)
        elif is_datetime64_dtype(dtype):
            typ = "datetime"
        elif is_timedelta64_dtype(dtype):
            typ = "timedelta"
        elif is_object_dtype(dtype):
            typ = "object"
        elif is_bool_dtype(dtype):
            typ = "bool"
        elif is_extension_array_dtype(dtype):
            typ = str(arr.dtype)
        else:
            typ = dtype.kind
        typs.add(typ)
    return typs


def concat_compat(to_concat, axis: int = 0):
    """
    provide concatenation of an array of arrays each of which is a single
    'normalized' dtypes (in that for example, if it's object, then it is a
    non-datetimelike and provide a combined dtype for the resulting array that
    preserves the overall dtype if possible)

    Parameters
    ----------
    to_concat : array of arrays
    axis : axis to provide concatenation

    Returns
    -------
    a single array, preserving the combined dtypes
    """
    # filter empty arrays
    # 1-d dtypes always are included here
    def is_nonempty(x) -> bool:
        if x.ndim <= axis:
            return True
        return x.shape[axis] > 0

    # If all arrays are empty, there's nothing to convert, just short-cut to
    # the concatenation, #3121.
    #
    # Creating an empty array directly is tempting, but the winnings would be
    # marginal given that it would still require shape & dtype calculation and
    # np.concatenate which has them both implemented is compiled.
    non_empties = [x for x in to_concat if is_nonempty(x)]
    if non_empties and axis == 0:
        to_concat = non_empties

    typs = get_dtype_kinds(to_concat)
    _contains_datetime = any(typ.startswith("datetime") for typ in typs)
    _contains_period = any(typ.startswith("period") for typ in typs)

    if "category" in typs:
        # this must be prior to concat_datetime,
        # to support Categorical + datetime-like
        from pandas import Categorical

        return Categorical._concat_arrays(to_concat, axis=axis)

    elif _contains_datetime or "timedelta" in typs or _contains_period:
        return concat_datetime(to_concat, axis=axis, typs=typs)

    # these are mandated to handle empties as well
    elif "sparse" in typs:
        from pandas.core.arrays import SparseArray

        return SparseArray._concat_arrays(to_concat, axis=axis)

    all_empty = not len(non_empties)
    single_dtype = len({x.dtype for x in to_concat}) == 1
    any_ea = any(is_extension_array_dtype(x.dtype) for x in to_concat)

    if any_ea and axis == 1:
        to_concat = [np.atleast_2d(x.astype("object")) for x in to_concat]

    elif any_ea and single_dtype and axis == 0:
        cls = type(to_concat[0])
        return cls._concat_same_type(to_concat)

    if all_empty:
        # we have all empties, but may need to coerce the result dtype to
        # object if we have non-numeric type operands (numpy would otherwise
        # cast this to float)
        typs = get_dtype_kinds(to_concat)
        if len(typs) != 1:

            if not len(typs - {"i", "u", "f"}) or not len(typs - {"bool", "i", "u"}):
                # let numpy coerce
                pass
            else:
                # coerce to object
                to_concat = [x.astype("object") for x in to_concat]

    return np.concatenate(to_concat, axis=axis)


def union_categoricals(
    to_union, sort_categories: bool = False, ignore_order: bool = False
):
    """
    Combine list-like of Categorical-like, unioning categories.

    All categories must have the same dtype.

    Parameters
    ----------
    to_union : list-like
        Categorical, CategoricalIndex, or Series with dtype='category'.
    sort_categories : bool, default False
        If true, resulting categories will be lexsorted, otherwise
        they will be ordered as they appear in the data.
    ignore_order : bool, default False
        If true, the ordered attribute of the Categoricals will be ignored.
        Results in an unordered categorical.

    Returns
    -------
    Categorical

    Raises
    ------
    TypeError
        - all inputs do not have the same dtype
        - all inputs do not have the same ordered property
        - all inputs are ordered and their categories are not identical
        - sort_categories=True and Categoricals are ordered
    ValueError
        Empty list of categoricals passed

    Notes
    -----
    To learn more about categories, see `link
    <https://pandas.pydata.org/pandas-docs/stable/user_guide/categorical.html#unioning>`__

    Examples
    --------
    >>> from pandas.api.types import union_categoricals

    If you want to combine categoricals that do not necessarily have
    the same categories, `union_categoricals` will combine a list-like
    of categoricals. The new categories will be the union of the
    categories being combined.

    >>> a = pd.Categorical(["b", "c"])
    >>> b = pd.Categorical(["a", "b"])
    >>> union_categoricals([a, b])
    [b, c, a, b]
    Categories (3, object): [b, c, a]

    By default, the resulting categories will be ordered as they appear
    in the `categories` of the data. If you want the categories to be
    lexsorted, use `sort_categories=True` argument.

    >>> union_categoricals([a, b], sort_categories=True)
    [b, c, a, b]
    Categories (3, object): [a, b, c]

    `union_categoricals` also works with the case of combining two
    categoricals of the same categories and order information (e.g. what
    you could also `append` for).

    >>> a = pd.Categorical(["a", "b"], ordered=True)
    >>> b = pd.Categorical(["a", "b", "a"], ordered=True)
    >>> union_categoricals([a, b])
    [a, b, a, b, a]
    Categories (2, object): [a < b]

    Raises `TypeError` because the categories are ordered and not identical.

    >>> a = pd.Categorical(["a", "b"], ordered=True)
    >>> b = pd.Categorical(["a", "b", "c"], ordered=True)
    >>> union_categoricals([a, b])
    Traceback (most recent call last):
        ...
    TypeError: to union ordered Categoricals, all categories must be the same

    New in version 0.20.0

    Ordered categoricals with different categories or orderings can be
    combined by using the `ignore_ordered=True` argument.

    >>> a = pd.Categorical(["a", "b", "c"], ordered=True)
    >>> b = pd.Categorical(["c", "b", "a"], ordered=True)
    >>> union_categoricals([a, b], ignore_order=True)
    [a, b, c, c, b, a]
    Categories (3, object): [a, b, c]

    `union_categoricals` also works with a `CategoricalIndex`, or `Series`
    containing categorical data, but note that the resulting array will
    always be a plain `Categorical`

    >>> a = pd.Series(["b", "c"], dtype='category')
    >>> b = pd.Series(["a", "b"], dtype='category')
    >>> union_categoricals([a, b])
    [b, c, a, b]
    Categories (3, object): [b, c, a]
    """
    from pandas import Index, Categorical
    from pandas.core.arrays.categorical import recode_for_categories

    if len(to_union) == 0:
        raise ValueError("No Categoricals to union")

    def _maybe_unwrap(x):
        if isinstance(x, (ABCCategoricalIndex, ABCSeries)):
            return x.values
        elif isinstance(x, Categorical):
            return x
        else:
            raise TypeError("all components to combine must be Categorical")

    to_union = [_maybe_unwrap(x) for x in to_union]
    first = to_union[0]

    if not all(
        is_dtype_equal(other.categories.dtype, first.categories.dtype)
        for other in to_union[1:]
    ):
        raise TypeError("dtype of categories must be the same")

    ordered = False
    if all(first.is_dtype_equal(other) for other in to_union[1:]):
        # identical categories - fastpath
        return Categorical._concat_same_dtype(
            to_union, sort_categories=sort_categories, ignore_order=ignore_order,
        )

    elif ignore_order or all(not c.ordered for c in to_union):
        # different categories - union and recode
        cats = first.categories.append([c.categories for c in to_union[1:]])
        categories = Index(cats.unique())
        if sort_categories:
            categories = categories.sort_values()

        new_codes = [
            recode_for_categories(c.codes, c.categories, categories) for c in to_union
        ]
        new_codes = np.concatenate(new_codes)
    else:
        # ordered - to show a proper error message
        if all(c.ordered for c in to_union):
            msg = "to union ordered Categoricals, all categories must be the same"
            raise TypeError(msg)
        else:
            raise TypeError("Categorical.ordered must be the same")

    if ignore_order:
        ordered = False

    return Categorical(new_codes, categories=categories, ordered=ordered, fastpath=True)


def concat_datetime(to_concat, axis=0, typs=None):
    """
    provide concatenation of an datetimelike array of arrays each of which is a
    single M8[ns], datetimet64[ns, tz] or m8[ns] dtype

    Parameters
    ----------
    to_concat : array of arrays
    axis : axis to provide concatenation
    typs : set of to_concat dtypes

    Returns
    -------
    a single array, preserving the combined dtypes
    """
    from pandas.core.arrays import datetimelike as dtl
    from pandas.core.ops.array_ops import maybe_upcast_datetimelike_array

    to_concat = [maybe_upcast_datetimelike_array(x) for x in to_concat]

    obj = [x for x in to_concat if isinstance(x, dtl.DatetimeLikeArrayMixin)][0]
    return type(obj)._concat_arrays(to_concat, axis=axis)
