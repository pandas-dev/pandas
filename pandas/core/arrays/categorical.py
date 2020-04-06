import operator
from shutil import get_terminal_size
from typing import Dict, Hashable, List, Type, Union, cast
from warnings import warn

import numpy as np

from pandas._config import get_option

from pandas._libs import algos as libalgos, hashtable as htable
from pandas._typing import ArrayLike, Dtype, Ordered, Scalar
from pandas.compat.numpy import function as nv
from pandas.util._decorators import (
    Appender,
    Substitution,
    cache_readonly,
    deprecate_kwarg,
    doc,
)
from pandas.util._validators import validate_bool_kwarg, validate_fillna_kwargs

from pandas.core.dtypes.cast import (
    coerce_indexer_dtype,
    maybe_cast_to_extension_array,
    maybe_infer_to_datetimelike,
)
from pandas.core.dtypes.common import (
    ensure_int64,
    ensure_object,
    ensure_platform_int,
    is_categorical_dtype,
    is_datetime64_dtype,
    is_dict_like,
    is_dtype_equal,
    is_extension_array_dtype,
    is_integer_dtype,
    is_iterator,
    is_list_like,
    is_object_dtype,
    is_scalar,
    is_sequence,
    is_timedelta64_dtype,
    needs_i8_conversion,
)
from pandas.core.dtypes.dtypes import CategoricalDtype
from pandas.core.dtypes.generic import ABCIndexClass, ABCSeries
from pandas.core.dtypes.inference import is_hashable
from pandas.core.dtypes.missing import isna, notna

from pandas.core import ops
from pandas.core.accessor import PandasDelegate, delegate_names
import pandas.core.algorithms as algorithms
from pandas.core.algorithms import _get_data_algo, factorize, take, take_1d, unique1d
from pandas.core.arrays.base import ExtensionArray, _extension_array_shared_docs
from pandas.core.base import NoNewAttributesMixin, PandasObject, _shared_docs
import pandas.core.common as com
from pandas.core.construction import array, extract_array, sanitize_array
from pandas.core.indexers import check_array_indexer, deprecate_ndim_indexing
from pandas.core.missing import interpolate_2d
from pandas.core.ops.common import unpack_zerodim_and_defer
from pandas.core.sorting import nargsort

from pandas.io.formats import console


def _cat_compare_op(op):
    opname = f"__{op.__name__}__"

    @unpack_zerodim_and_defer(opname)
    def func(self, other):
        if is_list_like(other) and len(other) != len(self):
            # TODO: Could this fail if the categories are listlike objects?
            raise ValueError("Lengths must match.")

        if not self.ordered:
            if opname in ["__lt__", "__gt__", "__le__", "__ge__"]:
                raise TypeError(
                    "Unordered Categoricals can only compare equality or not"
                )
        if isinstance(other, Categorical):
            # Two Categoricals can only be be compared if the categories are
            # the same (maybe up to ordering, depending on ordered)

            msg = "Categoricals can only be compared if 'categories' are the same."
            if len(self.categories) != len(other.categories):
                raise TypeError(msg + " Categories are different lengths")
            elif self.ordered and not (self.categories == other.categories).all():
                raise TypeError(msg)
            elif not set(self.categories) == set(other.categories):
                raise TypeError(msg)

            if not (self.ordered == other.ordered):
                raise TypeError(
                    "Categoricals can only be compared if 'ordered' is the same"
                )
            if not self.ordered and not self.categories.equals(other.categories):
                # both unordered and different order
                other_codes = _get_codes_for_values(other, self.categories)
            else:
                other_codes = other._codes

            f = getattr(self._codes, opname)
            ret = f(other_codes)
            mask = (self._codes == -1) | (other_codes == -1)
            if mask.any():
                # In other series, the leads to False, so do that here too
                if opname == "__ne__":
                    ret[(self._codes == -1) & (other_codes == -1)] = True
                else:
                    ret[mask] = False
            return ret

        if is_scalar(other):
            if other in self.categories:
                i = self.categories.get_loc(other)
                ret = getattr(self._codes, opname)(i)

                if opname not in {"__eq__", "__ge__", "__gt__"}:
                    # check for NaN needed if we are not equal or larger
                    mask = self._codes == -1
                    ret[mask] = False
                return ret
            else:
                if opname == "__eq__":
                    return np.zeros(len(self), dtype=bool)
                elif opname == "__ne__":
                    return np.ones(len(self), dtype=bool)
                else:
                    raise TypeError(
                        f"Cannot compare a Categorical for op {opname} with a "
                        "scalar, which is not a category."
                    )
        else:

            # allow categorical vs object dtype array comparisons for equality
            # these are only positional comparisons
            if opname in ["__eq__", "__ne__"]:
                return getattr(np.array(self), opname)(np.array(other))

            raise TypeError(
                f"Cannot compare a Categorical for op {opname} with "
                f"type {type(other)}.\nIf you want to compare values, "
                "use 'np.asarray(cat) <op> other'."
            )

    func.__name__ = opname

    return func


def contains(cat, key, container):
    """
    Helper for membership check for ``key`` in ``cat``.

    This is a helper method for :method:`__contains__`
    and :class:`CategoricalIndex.__contains__`.

    Returns True if ``key`` is in ``cat.categories`` and the
    location of ``key`` in ``categories`` is in ``container``.

    Parameters
    ----------
    cat : :class:`Categorical`or :class:`categoricalIndex`
    key : a hashable object
        The key to check membership for.
    container : Container (e.g. list-like or mapping)
        The container to check for membership in.

    Returns
    -------
    is_in : bool
        True if ``key`` is in ``self.categories`` and location of
        ``key`` in ``categories`` is in ``container``, else False.

    Notes
    -----
    This method does not check for NaN values. Do that separately
    before calling this method.
    """
    hash(key)

    # get location of key in categories.
    # If a KeyError, the key isn't in categories, so logically
    #  can't be in container either.
    try:
        loc = cat.categories.get_loc(key)
    except (KeyError, TypeError):
        return False

    # loc is the location of key in categories, but also the *value*
    # for key in container. So, `key` may be in categories,
    # but still not in `container`. Example ('b' in categories,
    # but not in values):
    # 'b' in Categorical(['a'], categories=['a', 'b'])  # False
    if is_scalar(loc):
        return loc in container
    else:
        # if categories is an IntervalIndex, loc is an array.
        return any(loc_ in container for loc_ in loc)


_codes_doc = """
The category codes of this categorical.

Level codes are an array if integer which are the positions of the real
values in the categories array.

There is not setter, use the other categorical methods and the normal item
setter to change values in the categorical.
"""


class Categorical(ExtensionArray, PandasObject):
    """
    Represent a categorical variable in classic R / S-plus fashion.

    `Categoricals` can only take on only a limited, and usually fixed, number
    of possible values (`categories`). In contrast to statistical categorical
    variables, a `Categorical` might have an order, but numerical operations
    (additions, divisions, ...) are not possible.

    All values of the `Categorical` are either in `categories` or `np.nan`.
    Assigning values outside of `categories` will raise a `ValueError`. Order
    is defined by the order of the `categories`, not lexical order of the
    values.

    Parameters
    ----------
    values : list-like
        The values of the categorical. If categories are given, values not in
        categories will be replaced with NaN.
    categories : Index-like (unique), optional
        The unique categories for this categorical. If not given, the
        categories are assumed to be the unique values of `values` (sorted, if
        possible, otherwise in the order in which they appear).
    ordered : bool, default False
        Whether or not this categorical is treated as a ordered categorical.
        If True, the resulting categorical will be ordered.
        An ordered categorical respects, when sorted, the order of its
        `categories` attribute (which in turn is the `categories` argument, if
        provided).
    dtype : CategoricalDtype
        An instance of ``CategoricalDtype`` to use for this categorical.

    Attributes
    ----------
    categories : Index
        The categories of this categorical
    codes : ndarray
        The codes (integer positions, which point to the categories) of this
        categorical, read only.
    ordered : bool
        Whether or not this Categorical is ordered.
    dtype : CategoricalDtype
        The instance of ``CategoricalDtype`` storing the ``categories``
        and ``ordered``.

    Methods
    -------
    from_codes
    __array__

    Raises
    ------
    ValueError
        If the categories do not validate.
    TypeError
        If an explicit ``ordered=True`` is given but no `categories` and the
        `values` are not sortable.

    See Also
    --------
    CategoricalDtype : Type for categorical data.
    CategoricalIndex : An Index with an underlying ``Categorical``.

    Notes
    -----
    See the `user guide
    <https://pandas.pydata.org/pandas-docs/stable/user_guide/categorical.html>`_
    for more.

    Examples
    --------
    >>> pd.Categorical([1, 2, 3, 1, 2, 3])
    [1, 2, 3, 1, 2, 3]
    Categories (3, int64): [1, 2, 3]

    >>> pd.Categorical(['a', 'b', 'c', 'a', 'b', 'c'])
    [a, b, c, a, b, c]
    Categories (3, object): [a, b, c]

    Ordered `Categoricals` can be sorted according to the custom order
    of the categories and can have a min and max value.

    >>> c = pd.Categorical(['a', 'b', 'c', 'a', 'b', 'c'], ordered=True,
    ...                    categories=['c', 'b', 'a'])
    >>> c
    [a, b, c, a, b, c]
    Categories (3, object): [c < b < a]
    >>> c.min()
    'c'
    """

    # For comparisons, so that numpy uses our implementation if the compare
    # ops, which raise
    __array_priority__ = 1000
    _dtype = CategoricalDtype(ordered=False)
    # tolist is not actually deprecated, just suppressed in the __dir__
    _deprecations = PandasObject._deprecations | frozenset(["tolist"])
    _typ = "categorical"

    def __init__(
        self, values, categories=None, ordered=None, dtype=None, fastpath=False
    ):

        dtype = CategoricalDtype._from_values_or_dtype(
            values, categories, ordered, dtype
        )
        # At this point, dtype is always a CategoricalDtype, but
        # we may have dtype.categories be None, and we need to
        # infer categories in a factorization step further below

        if fastpath:
            self._codes = coerce_indexer_dtype(values, dtype.categories)
            self._dtype = self._dtype.update_dtype(dtype)
            return

        # null_mask indicates missing values we want to exclude from inference.
        # This means: only missing values in list-likes (not arrays/ndframes).
        null_mask = np.array(False)

        # sanitize input
        if is_categorical_dtype(values):
            if dtype.categories is None:
                dtype = CategoricalDtype(values.categories, dtype.ordered)
        elif not isinstance(values, (ABCIndexClass, ABCSeries)):
            # sanitize_array coerces np.nan to a string under certain versions
            # of numpy
            values = maybe_infer_to_datetimelike(values, convert_dates=True)
            if not isinstance(values, np.ndarray):
                values = _convert_to_list_like(values)

                # By convention, empty lists result in object dtype:
                sanitize_dtype = "object" if len(values) == 0 else None
                null_mask = isna(values)
                if null_mask.any():
                    values = [values[idx] for idx in np.where(~null_mask)[0]]
                values = sanitize_array(values, None, dtype=sanitize_dtype)

        if dtype.categories is None:
            try:
                codes, categories = factorize(values, sort=True)
            except TypeError as err:
                codes, categories = factorize(values, sort=False)
                if dtype.ordered:
                    # raise, as we don't have a sortable data structure and so
                    # the user should give us one by specifying categories
                    raise TypeError(
                        "'values' is not ordered, please "
                        "explicitly specify the categories order "
                        "by passing in a categories argument."
                    ) from err
            except ValueError as err:

                # FIXME
                raise NotImplementedError(
                    "> 1 ndim Categorical are not supported at this time"
                ) from err

            # we're inferring from values
            dtype = CategoricalDtype(categories, dtype.ordered)

        elif is_categorical_dtype(values):
            old_codes = (
                values._values.codes if isinstance(values, ABCSeries) else values.codes
            )
            codes = recode_for_categories(
                old_codes, values.dtype.categories, dtype.categories
            )

        else:
            codes = _get_codes_for_values(values, dtype.categories)

        if null_mask.any():
            # Reinsert -1 placeholders for previously removed missing values
            full_codes = -np.ones(null_mask.shape, dtype=codes.dtype)
            full_codes[~null_mask] = codes
            codes = full_codes

        self._dtype = self._dtype.update_dtype(dtype)
        self._codes = coerce_indexer_dtype(codes, dtype.categories)

    @property
    def categories(self):
        """
        The categories of this categorical.

        Setting assigns new values to each category (effectively a rename of
        each individual category).

        The assigned value has to be a list-like object. All items must be
        unique and the number of items in the new categories must be the same
        as the number of items in the old categories.

        Assigning to `categories` is a inplace operation!

        Raises
        ------
        ValueError
            If the new categories do not validate as categories or if the
            number of new categories is unequal the number of old categories

        See Also
        --------
        rename_categories : Rename categories.
        reorder_categories : Reorder categories.
        add_categories : Add new categories.
        remove_categories : Remove the specified categories.
        remove_unused_categories : Remove categories which are not used.
        set_categories : Set the categories to the specified ones.
        """
        return self.dtype.categories

    @categories.setter
    def categories(self, categories):
        new_dtype = CategoricalDtype(categories, ordered=self.ordered)
        if self.dtype.categories is not None and len(self.dtype.categories) != len(
            new_dtype.categories
        ):
            raise ValueError(
                "new categories need to have the same number of "
                "items as the old categories!"
            )
        self._dtype = new_dtype

    @property
    def ordered(self) -> Ordered:
        """
        Whether the categories have an ordered relationship.
        """
        return self.dtype.ordered

    @property
    def dtype(self) -> CategoricalDtype:
        """
        The :class:`~pandas.api.types.CategoricalDtype` for this instance.
        """
        return self._dtype

    @property
    def _constructor(self) -> Type["Categorical"]:
        return Categorical

    @classmethod
    def _from_sequence(cls, scalars, dtype=None, copy=False):
        return Categorical(scalars, dtype=dtype)

    def _formatter(self, boxed=False):
        # Defer to CategoricalFormatter's formatter.
        return None

    def copy(self) -> "Categorical":
        """
        Copy constructor.
        """
        return self._constructor(
            values=self._codes.copy(), dtype=self.dtype, fastpath=True
        )

    def astype(self, dtype: Dtype, copy: bool = True) -> ArrayLike:
        """
        Coerce this type to another dtype

        Parameters
        ----------
        dtype : numpy dtype or pandas type
        copy : bool, default True
            By default, astype always returns a newly allocated object.
            If copy is set to False and dtype is categorical, the original
            object is returned.
        """
        if is_categorical_dtype(dtype):
            dtype = cast(Union[str, CategoricalDtype], dtype)

            # GH 10696/18593
            dtype = self.dtype.update_dtype(dtype)
            self = self.copy() if copy else self
            if dtype == self.dtype:
                return self
            return self._set_dtype(dtype)
        if is_extension_array_dtype(dtype):
            return array(self, dtype=dtype, copy=copy)  # type: ignore # GH 28770
        if is_integer_dtype(dtype) and self.isna().any():
            raise ValueError("Cannot convert float NaN to integer")
        return np.array(self, dtype=dtype, copy=copy)

    @cache_readonly
    def size(self) -> int:
        """
        Return the len of myself.
        """
        return self._codes.size

    @cache_readonly
    def itemsize(self) -> int:
        """
        return the size of a single category
        """
        return self.categories.itemsize

    def tolist(self) -> List[Scalar]:
        """
        Return a list of the values.

        These are each a scalar type, which is a Python scalar
        (for str, int, float) or a pandas scalar
        (for Timestamp/Timedelta/Interval/Period)
        """
        return list(self)

    to_list = tolist

    @classmethod
    def _from_inferred_categories(
        cls, inferred_categories, inferred_codes, dtype, true_values=None
    ):
        """
        Construct a Categorical from inferred values.

        For inferred categories (`dtype` is None) the categories are sorted.
        For explicit `dtype`, the `inferred_categories` are cast to the
        appropriate type.

        Parameters
        ----------
        inferred_categories : Index
        inferred_codes : Index
        dtype : CategoricalDtype or 'category'
        true_values : list, optional
            If none are provided, the default ones are
            "True", "TRUE", and "true."

        Returns
        -------
        Categorical
        """
        from pandas import Index, to_numeric, to_datetime, to_timedelta

        cats = Index(inferred_categories)
        known_categories = (
            isinstance(dtype, CategoricalDtype) and dtype.categories is not None
        )

        if known_categories:
            # Convert to a specialized type with `dtype` if specified.
            if dtype.categories.is_numeric():
                cats = to_numeric(inferred_categories, errors="coerce")
            elif is_datetime64_dtype(dtype.categories):
                cats = to_datetime(inferred_categories, errors="coerce")
            elif is_timedelta64_dtype(dtype.categories):
                cats = to_timedelta(inferred_categories, errors="coerce")
            elif dtype.categories.is_boolean():
                if true_values is None:
                    true_values = ["True", "TRUE", "true"]

                cats = cats.isin(true_values)

        if known_categories:
            # Recode from observation order to dtype.categories order.
            categories = dtype.categories
            codes = recode_for_categories(inferred_codes, cats, categories)
        elif not cats.is_monotonic_increasing:
            # Sort categories and recode for unknown categories.
            unsorted = cats.copy()
            categories = cats.sort_values()

            codes = recode_for_categories(inferred_codes, unsorted, categories)
            dtype = CategoricalDtype(categories, ordered=False)
        else:
            dtype = CategoricalDtype(cats, ordered=False)
            codes = inferred_codes

        return cls(codes, dtype=dtype, fastpath=True)

    @classmethod
    def from_codes(cls, codes, categories=None, ordered=None, dtype=None):
        """
        Make a Categorical type from codes and categories or dtype.

        This constructor is useful if you already have codes and
        categories/dtype and so do not need the (computation intensive)
        factorization step, which is usually done on the constructor.

        If your data does not follow this convention, please use the normal
        constructor.

        Parameters
        ----------
        codes : array-like of int
            An integer array, where each integer points to a category in
            categories or dtype.categories, or else is -1 for NaN.
        categories : index-like, optional
            The categories for the categorical. Items need to be unique.
            If the categories are not given here, then they must be provided
            in `dtype`.
        ordered : bool, optional
            Whether or not this categorical is treated as an ordered
            categorical. If not given here or in `dtype`, the resulting
            categorical will be unordered.
        dtype : CategoricalDtype or "category", optional
            If :class:`CategoricalDtype`, cannot be used together with
            `categories` or `ordered`.

            .. versionadded:: 0.24.0

               When `dtype` is provided, neither `categories` nor `ordered`
               should be provided.

        Returns
        -------
        Categorical

        Examples
        --------
        >>> dtype = pd.CategoricalDtype(['a', 'b'], ordered=True)
        >>> pd.Categorical.from_codes(codes=[0, 1, 0, 1], dtype=dtype)
        [a, b, a, b]
        Categories (2, object): [a < b]
        """
        dtype = CategoricalDtype._from_values_or_dtype(
            categories=categories, ordered=ordered, dtype=dtype
        )
        if dtype.categories is None:
            msg = (
                "The categories must be provided in 'categories' or "
                "'dtype'. Both were None."
            )
            raise ValueError(msg)

        if is_extension_array_dtype(codes) and is_integer_dtype(codes):
            # Avoid the implicit conversion of Int to object
            if isna(codes).any():
                raise ValueError("codes cannot contain NA values")
            codes = codes.to_numpy(dtype=np.int64)
        else:
            codes = np.asarray(codes)
        if len(codes) and not is_integer_dtype(codes):
            raise ValueError("codes need to be array-like integers")

        if len(codes) and (codes.max() >= len(dtype.categories) or codes.min() < -1):
            raise ValueError("codes need to be between -1 and len(categories)-1")

        return cls(codes, dtype=dtype, fastpath=True)

    def _get_codes(self):
        """
        Get the codes.

        Returns
        -------
        codes : integer array view
            A non writable view of the `codes` array.
        """
        v = self._codes.view()
        v.flags.writeable = False
        return v

    def _set_codes(self, codes):
        """
        Not settable by the user directly
        """
        raise ValueError("cannot set Categorical codes directly")

    codes = property(fget=_get_codes, fset=_set_codes, doc=_codes_doc)

    def _set_categories(self, categories, fastpath=False):
        """
        Sets new categories inplace

        Parameters
        ----------
        fastpath : bool, default False
           Don't perform validation of the categories for uniqueness or nulls

        Examples
        --------
        >>> c = pd.Categorical(['a', 'b'])
        >>> c
        [a, b]
        Categories (2, object): [a, b]

        >>> c._set_categories(pd.Index(['a', 'c']))
        >>> c
        [a, c]
        Categories (2, object): [a, c]
        """
        if fastpath:
            new_dtype = CategoricalDtype._from_fastpath(categories, self.ordered)
        else:
            new_dtype = CategoricalDtype(categories, ordered=self.ordered)
        if (
            not fastpath
            and self.dtype.categories is not None
            and len(new_dtype.categories) != len(self.dtype.categories)
        ):
            raise ValueError(
                "new categories need to have the same number of "
                "items than the old categories!"
            )

        self._dtype = new_dtype

    def _set_dtype(self, dtype: CategoricalDtype) -> "Categorical":
        """
        Internal method for directly updating the CategoricalDtype

        Parameters
        ----------
        dtype : CategoricalDtype

        Notes
        -----
        We don't do any validation here. It's assumed that the dtype is
        a (valid) instance of `CategoricalDtype`.
        """
        codes = recode_for_categories(self.codes, self.categories, dtype.categories)
        return type(self)(codes, dtype=dtype, fastpath=True)

    def set_ordered(self, value, inplace=False):
        """
        Set the ordered attribute to the boolean value.

        Parameters
        ----------
        value : bool
           Set whether this categorical is ordered (True) or not (False).
        inplace : bool, default False
           Whether or not to set the ordered attribute in-place or return
           a copy of this categorical with ordered set to the value.
        """
        inplace = validate_bool_kwarg(inplace, "inplace")
        new_dtype = CategoricalDtype(self.categories, ordered=value)
        cat = self if inplace else self.copy()
        cat._dtype = new_dtype
        if not inplace:
            return cat

    def as_ordered(self, inplace=False):
        """
        Set the Categorical to be ordered.

        Parameters
        ----------
        inplace : bool, default False
           Whether or not to set the ordered attribute in-place or return
           a copy of this categorical with ordered set to True.

        Returns
        -------
        Categorical
            Ordered Categorical.
        """
        inplace = validate_bool_kwarg(inplace, "inplace")
        return self.set_ordered(True, inplace=inplace)

    def as_unordered(self, inplace=False):
        """
        Set the Categorical to be unordered.

        Parameters
        ----------
        inplace : bool, default False
           Whether or not to set the ordered attribute in-place or return
           a copy of this categorical with ordered set to False.

        Returns
        -------
        Categorical
            Unordered Categorical.
        """
        inplace = validate_bool_kwarg(inplace, "inplace")
        return self.set_ordered(False, inplace=inplace)

    def set_categories(self, new_categories, ordered=None, rename=False, inplace=False):
        """
        Set the categories to the specified new_categories.

        `new_categories` can include new categories (which will result in
        unused categories) or remove old categories (which results in values
        set to NaN). If `rename==True`, the categories will simple be renamed
        (less or more items than in old categories will result in values set to
        NaN or in unused categories respectively).

        This method can be used to perform more than one action of adding,
        removing, and reordering simultaneously and is therefore faster than
        performing the individual steps via the more specialised methods.

        On the other hand this methods does not do checks (e.g., whether the
        old categories are included in the new categories on a reorder), which
        can result in surprising changes, for example when using special string
        dtypes, which does not considers a S1 string equal to a single char
        python string.

        Parameters
        ----------
        new_categories : Index-like
           The categories in new order.
        ordered : bool, default False
           Whether or not the categorical is treated as a ordered categorical.
           If not given, do not change the ordered information.
        rename : bool, default False
           Whether or not the new_categories should be considered as a rename
           of the old categories or as reordered categories.
        inplace : bool, default False
           Whether or not to reorder the categories in-place or return a copy
           of this categorical with reordered categories.

        Returns
        -------
        Categorical with reordered categories or None if inplace.

        Raises
        ------
        ValueError
            If new_categories does not validate as categories

        See Also
        --------
        rename_categories : Rename categories.
        reorder_categories : Reorder categories.
        add_categories : Add new categories.
        remove_categories : Remove the specified categories.
        remove_unused_categories : Remove categories which are not used.
        """
        inplace = validate_bool_kwarg(inplace, "inplace")
        if ordered is None:
            ordered = self.dtype.ordered
        new_dtype = CategoricalDtype(new_categories, ordered=ordered)

        cat = self if inplace else self.copy()
        if rename:
            if cat.dtype.categories is not None and len(new_dtype.categories) < len(
                cat.dtype.categories
            ):
                # remove all _codes which are larger and set to -1/NaN
                cat._codes[cat._codes >= len(new_dtype.categories)] = -1
        else:
            codes = recode_for_categories(
                cat.codes, cat.categories, new_dtype.categories
            )
            cat._codes = codes
        cat._dtype = new_dtype

        if not inplace:
            return cat

    def rename_categories(self, new_categories, inplace=False):
        """
        Rename categories.

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

            .. versionadded:: 0.23.0.

        inplace : bool, default False
            Whether or not to rename the categories inplace or return a copy of
            this categorical with renamed categories.

        Returns
        -------
        cat : Categorical or None
           With ``inplace=False``, the new categorical is returned.
           With ``inplace=True``, there is no return value.

        Raises
        ------
        ValueError
            If new categories are list-like and do not have the same number of
            items than the current categories or do not validate as categories

        See Also
        --------
        reorder_categories : Reorder categories.
        add_categories : Add new categories.
        remove_categories : Remove the specified categories.
        remove_unused_categories : Remove categories which are not used.
        set_categories : Set the categories to the specified ones.

        Examples
        --------
        >>> c = pd.Categorical(['a', 'a', 'b'])
        >>> c.rename_categories([0, 1])
        [0, 0, 1]
        Categories (2, int64): [0, 1]

        For dict-like ``new_categories``, extra keys are ignored and
        categories not in the dictionary are passed through

        >>> c.rename_categories({'a': 'A', 'c': 'C'})
        [A, A, b]
        Categories (2, object): [A, b]

        You may also provide a callable to create the new categories

        >>> c.rename_categories(lambda x: x.upper())
        [A, A, B]
        Categories (2, object): [A, B]
        """
        inplace = validate_bool_kwarg(inplace, "inplace")
        cat = self if inplace else self.copy()

        if is_dict_like(new_categories):
            cat.categories = [new_categories.get(item, item) for item in cat.categories]
        elif callable(new_categories):
            cat.categories = [new_categories(item) for item in cat.categories]
        else:
            cat.categories = new_categories
        if not inplace:
            return cat

    def reorder_categories(self, new_categories, ordered=None, inplace=False):
        """
        Reorder categories as specified in new_categories.

        `new_categories` need to include all old categories and no new category
        items.

        Parameters
        ----------
        new_categories : Index-like
           The categories in new order.
        ordered : bool, optional
           Whether or not the categorical is treated as a ordered categorical.
           If not given, do not change the ordered information.
        inplace : bool, default False
           Whether or not to reorder the categories inplace or return a copy of
           this categorical with reordered categories.

        Returns
        -------
        cat : Categorical with reordered categories or None if inplace.

        Raises
        ------
        ValueError
            If the new categories do not contain all old category items or any
            new ones

        See Also
        --------
        rename_categories : Rename categories.
        add_categories : Add new categories.
        remove_categories : Remove the specified categories.
        remove_unused_categories : Remove categories which are not used.
        set_categories : Set the categories to the specified ones.
        """
        inplace = validate_bool_kwarg(inplace, "inplace")
        if set(self.dtype.categories) != set(new_categories):
            raise ValueError(
                "items in new_categories are not the same as in old categories"
            )
        return self.set_categories(new_categories, ordered=ordered, inplace=inplace)

    def add_categories(self, new_categories, inplace=False):
        """
        Add new categories.

        `new_categories` will be included at the last/highest place in the
        categories and will be unused directly after this call.

        Parameters
        ----------
        new_categories : category or list-like of category
           The new categories to be included.
        inplace : bool, default False
           Whether or not to add the categories inplace or return a copy of
           this categorical with added categories.

        Returns
        -------
        cat : Categorical with new categories added or None if inplace.

        Raises
        ------
        ValueError
            If the new categories include old categories or do not validate as
            categories

        See Also
        --------
        rename_categories : Rename categories.
        reorder_categories : Reorder categories.
        remove_categories : Remove the specified categories.
        remove_unused_categories : Remove categories which are not used.
        set_categories : Set the categories to the specified ones.
        """
        inplace = validate_bool_kwarg(inplace, "inplace")
        if not is_list_like(new_categories):
            new_categories = [new_categories]
        already_included = set(new_categories) & set(self.dtype.categories)
        if len(already_included) != 0:
            raise ValueError(
                f"new categories must not include old categories: {already_included}"
            )
        new_categories = list(self.dtype.categories) + list(new_categories)
        new_dtype = CategoricalDtype(new_categories, self.ordered)

        cat = self if inplace else self.copy()
        cat._dtype = new_dtype
        cat._codes = coerce_indexer_dtype(cat._codes, new_dtype.categories)
        if not inplace:
            return cat

    def remove_categories(self, removals, inplace=False):
        """
        Remove the specified categories.

        `removals` must be included in the old categories. Values which were in
        the removed categories will be set to NaN

        Parameters
        ----------
        removals : category or list of categories
           The categories which should be removed.
        inplace : bool, default False
           Whether or not to remove the categories inplace or return a copy of
           this categorical with removed categories.

        Returns
        -------
        cat : Categorical with removed categories or None if inplace.

        Raises
        ------
        ValueError
            If the removals are not contained in the categories

        See Also
        --------
        rename_categories : Rename categories.
        reorder_categories : Reorder categories.
        add_categories : Add new categories.
        remove_unused_categories : Remove categories which are not used.
        set_categories : Set the categories to the specified ones.
        """
        inplace = validate_bool_kwarg(inplace, "inplace")
        if not is_list_like(removals):
            removals = [removals]

        removal_set = set(removals)
        not_included = removal_set - set(self.dtype.categories)
        new_categories = [c for c in self.dtype.categories if c not in removal_set]

        # GH 10156
        if any(isna(removals)):
            not_included = {x for x in not_included if notna(x)}
            new_categories = [x for x in new_categories if notna(x)]

        if len(not_included) != 0:
            raise ValueError(f"removals must all be in old categories: {not_included}")

        return self.set_categories(
            new_categories, ordered=self.ordered, rename=False, inplace=inplace
        )

    def remove_unused_categories(self, inplace=False):
        """
        Remove categories which are not used.

        Parameters
        ----------
        inplace : bool, default False
           Whether or not to drop unused categories inplace or return a copy of
           this categorical with unused categories dropped.

        Returns
        -------
        cat : Categorical with unused categories dropped or None if inplace.

        See Also
        --------
        rename_categories : Rename categories.
        reorder_categories : Reorder categories.
        add_categories : Add new categories.
        remove_categories : Remove the specified categories.
        set_categories : Set the categories to the specified ones.
        """
        inplace = validate_bool_kwarg(inplace, "inplace")
        cat = self if inplace else self.copy()
        idx, inv = np.unique(cat._codes, return_inverse=True)

        if idx.size != 0 and idx[0] == -1:  # na sentinel
            idx, inv = idx[1:], inv - 1

        new_categories = cat.dtype.categories.take(idx)
        new_dtype = CategoricalDtype._from_fastpath(
            new_categories, ordered=self.ordered
        )
        cat._dtype = new_dtype
        cat._codes = coerce_indexer_dtype(inv, new_dtype.categories)

        if not inplace:
            return cat

    def map(self, mapper):
        """
        Map categories using input correspondence (dict, Series, or function).

        Maps the categories to new categories. If the mapping correspondence is
        one-to-one the result is a :class:`~pandas.Categorical` which has the
        same order property as the original, otherwise a :class:`~pandas.Index`
        is returned. NaN values are unaffected.

        If a `dict` or :class:`~pandas.Series` is used any unmapped category is
        mapped to `NaN`. Note that if this happens an :class:`~pandas.Index`
        will be returned.

        Parameters
        ----------
        mapper : function, dict, or Series
            Mapping correspondence.

        Returns
        -------
        pandas.Categorical or pandas.Index
            Mapped categorical.

        See Also
        --------
        CategoricalIndex.map : Apply a mapping correspondence on a
            :class:`~pandas.CategoricalIndex`.
        Index.map : Apply a mapping correspondence on an
            :class:`~pandas.Index`.
        Series.map : Apply a mapping correspondence on a
            :class:`~pandas.Series`.
        Series.apply : Apply more complex functions on a
            :class:`~pandas.Series`.

        Examples
        --------
        >>> cat = pd.Categorical(['a', 'b', 'c'])
        >>> cat
        [a, b, c]
        Categories (3, object): [a, b, c]
        >>> cat.map(lambda x: x.upper())
        [A, B, C]
        Categories (3, object): [A, B, C]
        >>> cat.map({'a': 'first', 'b': 'second', 'c': 'third'})
        [first, second, third]
        Categories (3, object): [first, second, third]

        If the mapping is one-to-one the ordering of the categories is
        preserved:

        >>> cat = pd.Categorical(['a', 'b', 'c'], ordered=True)
        >>> cat
        [a, b, c]
        Categories (3, object): [a < b < c]
        >>> cat.map({'a': 3, 'b': 2, 'c': 1})
        [3, 2, 1]
        Categories (3, int64): [3 < 2 < 1]

        If the mapping is not one-to-one an :class:`~pandas.Index` is returned:

        >>> cat.map({'a': 'first', 'b': 'second', 'c': 'first'})
        Index(['first', 'second', 'first'], dtype='object')

        If a `dict` is used, all unmapped categories are mapped to `NaN` and
        the result is an :class:`~pandas.Index`:

        >>> cat.map({'a': 'first', 'b': 'second'})
        Index(['first', 'second', nan], dtype='object')
        """
        new_categories = self.categories.map(mapper)
        try:
            return self.from_codes(
                self._codes.copy(), categories=new_categories, ordered=self.ordered
            )
        except ValueError:
            # NA values are represented in self._codes with -1
            # np.take causes NA values to take final element in new_categories
            if np.any(self._codes == -1):
                new_categories = new_categories.insert(len(new_categories), np.nan)
            return np.take(new_categories, self._codes)

    __eq__ = _cat_compare_op(operator.eq)
    __ne__ = _cat_compare_op(operator.ne)
    __lt__ = _cat_compare_op(operator.lt)
    __gt__ = _cat_compare_op(operator.gt)
    __le__ = _cat_compare_op(operator.le)
    __ge__ = _cat_compare_op(operator.ge)

    # for Series/ndarray like compat
    @property
    def shape(self):
        """
        Shape of the Categorical.

        For internal compatibility with numpy arrays.

        Returns
        -------
        shape : tuple
        """
        return tuple([len(self._codes)])

    def shift(self, periods, fill_value=None):
        """
        Shift Categorical by desired number of periods.

        Parameters
        ----------
        periods : int
            Number of periods to move, can be positive or negative
        fill_value : object, optional
            The scalar value to use for newly introduced missing values.

            .. versionadded:: 0.24.0

        Returns
        -------
        shifted : Categorical
        """
        # since categoricals always have ndim == 1, an axis parameter
        # doesn't make any sense here.
        codes = self.codes
        if codes.ndim > 1:
            raise NotImplementedError("Categorical with ndim > 1.")
        if np.prod(codes.shape) and (periods != 0):
            codes = np.roll(codes, ensure_platform_int(periods), axis=0)
            if isna(fill_value):
                fill_value = -1
            elif fill_value in self.categories:
                fill_value = self.categories.get_loc(fill_value)
            else:
                raise ValueError(
                    f"'fill_value={fill_value}' is not present "
                    "in this Categorical's categories"
                )
            if periods > 0:
                codes[:periods] = fill_value
            else:
                codes[periods:] = fill_value

        return self.from_codes(codes, dtype=self.dtype)

    def __array__(self, dtype=None) -> np.ndarray:
        """
        The numpy array interface.

        Returns
        -------
        numpy.array
            A numpy array of either the specified dtype or,
            if dtype==None (default), the same dtype as
            categorical.categories.dtype.
        """
        ret = take_1d(self.categories.values, self._codes)
        if dtype and not is_dtype_equal(dtype, self.categories.dtype):
            return np.asarray(ret, dtype)
        if is_extension_array_dtype(ret):
            # When we're a Categorical[ExtensionArray], like Interval,
            # we need to ensure __array__ get's all the way to an
            # ndarray.
            ret = np.asarray(ret)
        return ret

    def __array_ufunc__(self, ufunc, method, *inputs, **kwargs):
        # for binary ops, use our custom dunder methods
        result = ops.maybe_dispatch_ufunc_to_dunder_op(
            self, ufunc, method, *inputs, **kwargs
        )
        if result is not NotImplemented:
            return result

        # for all other cases, raise for now (similarly as what happens in
        # Series.__array_prepare__)
        raise TypeError(
            f"Object with dtype {self.dtype} cannot perform "
            f"the numpy op {ufunc.__name__}"
        )

    def __setstate__(self, state):
        """Necessary for making this object picklable"""
        if not isinstance(state, dict):
            raise Exception("invalid pickle state")

        if "_dtype" not in state:
            state["_dtype"] = CategoricalDtype(state["_categories"], state["_ordered"])

        for k, v in state.items():
            setattr(self, k, v)

    @property
    def T(self) -> "Categorical":
        """
        Return transposed numpy array.
        """
        return self

    @property
    def nbytes(self):
        return self._codes.nbytes + self.dtype.categories.values.nbytes

    def memory_usage(self, deep=False):
        """
        Memory usage of my values

        Parameters
        ----------
        deep : bool
            Introspect the data deeply, interrogate
            `object` dtypes for system-level memory consumption

        Returns
        -------
        bytes used

        Notes
        -----
        Memory usage does not include memory consumed by elements that
        are not components of the array if deep=False

        See Also
        --------
        numpy.ndarray.nbytes
        """
        return self._codes.nbytes + self.dtype.categories.memory_usage(deep=deep)

    @doc(_shared_docs["searchsorted"], klass="Categorical")
    def searchsorted(self, value, side="left", sorter=None):
        # searchsorted is very performance sensitive. By converting codes
        # to same dtype as self.codes, we get much faster performance.
        if is_scalar(value):
            codes = self.categories.get_loc(value)
            codes = self.codes.dtype.type(codes)
        else:
            locs = [self.categories.get_loc(x) for x in value]
            codes = np.array(locs, dtype=self.codes.dtype)
        return self.codes.searchsorted(codes, side=side, sorter=sorter)

    def isna(self):
        """
        Detect missing values

        Missing values (-1 in .codes) are detected.

        Returns
        -------
        a boolean array of whether my values are null

        See Also
        --------
        isna : Top-level isna.
        isnull : Alias of isna.
        Categorical.notna : Boolean inverse of Categorical.isna.

        """
        ret = self._codes == -1
        return ret

    isnull = isna

    def notna(self):
        """
        Inverse of isna

        Both missing values (-1 in .codes) and NA as a category are detected as
        null.

        Returns
        -------
        a boolean array of whether my values are not null

        See Also
        --------
        notna : Top-level notna.
        notnull : Alias of notna.
        Categorical.isna : Boolean inverse of Categorical.notna.

        """
        return ~self.isna()

    notnull = notna

    def dropna(self):
        """
        Return the Categorical without null values.

        Missing values (-1 in .codes) are detected.

        Returns
        -------
        valid : Categorical
        """
        result = self[self.notna()]

        return result

    def value_counts(self, dropna=True):
        """
        Return a Series containing counts of each category.

        Every category will have an entry, even those with a count of 0.

        Parameters
        ----------
        dropna : bool, default True
            Don't include counts of NaN.

        Returns
        -------
        counts : Series

        See Also
        --------
        Series.value_counts
        """
        from pandas import Series, CategoricalIndex

        code, cat = self._codes, self.categories
        ncat, mask = len(cat), 0 <= code
        ix, clean = np.arange(ncat), mask.all()

        if dropna or clean:
            obs = code if clean else code[mask]
            count = np.bincount(obs, minlength=ncat or 0)
        else:
            count = np.bincount(np.where(mask, code, ncat))
            ix = np.append(ix, -1)

        ix = self._constructor(ix, dtype=self.dtype, fastpath=True)

        return Series(count, index=CategoricalIndex(ix), dtype="int64")

    def _internal_get_values(self):
        """
        Return the values.

        For internal compatibility with pandas formatting.

        Returns
        -------
        np.ndarray or Index
            A numpy array of the same dtype as categorical.categories.dtype or
            Index if datetime / periods.
        """
        # if we are a datetime and period index, return Index to keep metadata
        if needs_i8_conversion(self.categories):
            return self.categories.take(self._codes, fill_value=np.nan)
        elif is_integer_dtype(self.categories) and -1 in self._codes:
            return self.categories.astype("object").take(self._codes, fill_value=np.nan)
        return np.array(self)

    def check_for_ordered(self, op):
        """ assert that we are ordered """
        if not self.ordered:
            raise TypeError(
                f"Categorical is not ordered for operation {op}\n"
                "you can use .as_ordered() to change the "
                "Categorical to an ordered one\n"
            )

    def _values_for_argsort(self):
        return self._codes

    def argsort(self, ascending=True, kind="quicksort", **kwargs):
        """
        Return the indices that would sort the Categorical.

        .. versionchanged:: 0.25.0

           Changed to sort missing values at the end.

        Parameters
        ----------
        ascending : bool, default True
            Whether the indices should result in an ascending
            or descending sort.
        kind : {'quicksort', 'mergesort', 'heapsort'}, optional
            Sorting algorithm.
        **kwargs:
            passed through to :func:`numpy.argsort`.

        Returns
        -------
        numpy.array

        See Also
        --------
        numpy.ndarray.argsort

        Notes
        -----
        While an ordering is applied to the category values, arg-sorting
        in this context refers more to organizing and grouping together
        based on matching category values. Thus, this function can be
        called on an unordered Categorical instance unlike the functions
        'Categorical.min' and 'Categorical.max'.

        Examples
        --------
        >>> pd.Categorical(['b', 'b', 'a', 'c']).argsort()
        array([2, 0, 1, 3])

        >>> cat = pd.Categorical(['b', 'b', 'a', 'c'],
        ...                      categories=['c', 'b', 'a'],
        ...                      ordered=True)
        >>> cat.argsort()
        array([3, 0, 1, 2])

        Missing values are placed at the end

        >>> cat = pd.Categorical([2, None, 1])
        >>> cat.argsort()
        array([2, 0, 1])
        """
        return super().argsort(ascending=ascending, kind=kind, **kwargs)

    def sort_values(self, inplace=False, ascending=True, na_position="last"):
        """
        Sort the Categorical by category value returning a new
        Categorical by default.

        While an ordering is applied to the category values, sorting in this
        context refers more to organizing and grouping together based on
        matching category values. Thus, this function can be called on an
        unordered Categorical instance unlike the functions 'Categorical.min'
        and 'Categorical.max'.

        Parameters
        ----------
        inplace : bool, default False
            Do operation in place.
        ascending : bool, default True
            Order ascending. Passing False orders descending. The
            ordering parameter provides the method by which the
            category values are organized.
        na_position : {'first', 'last'} (optional, default='last')
            'first' puts NaNs at the beginning
            'last' puts NaNs at the end

        Returns
        -------
        Categorical or None

        See Also
        --------
        Categorical.sort
        Series.sort_values

        Examples
        --------
        >>> c = pd.Categorical([1, 2, 2, 1, 5])
        >>> c
        [1, 2, 2, 1, 5]
        Categories (3, int64): [1, 2, 5]
        >>> c.sort_values()
        [1, 1, 2, 2, 5]
        Categories (3, int64): [1, 2, 5]
        >>> c.sort_values(ascending=False)
        [5, 2, 2, 1, 1]
        Categories (3, int64): [1, 2, 5]

        Inplace sorting can be done as well:

        >>> c.sort_values(inplace=True)
        >>> c
        [1, 1, 2, 2, 5]
        Categories (3, int64): [1, 2, 5]
        >>>
        >>> c = pd.Categorical([1, 2, 2, 1, 5])

        'sort_values' behaviour with NaNs. Note that 'na_position'
        is independent of the 'ascending' parameter:

        >>> c = pd.Categorical([np.nan, 2, 2, np.nan, 5])
        >>> c
        [NaN, 2, 2, NaN, 5]
        Categories (2, int64): [2, 5]
        >>> c.sort_values()
        [2, 2, 5, NaN, NaN]
        Categories (2, int64): [2, 5]
        >>> c.sort_values(ascending=False)
        [5, 2, 2, NaN, NaN]
        Categories (2, int64): [2, 5]
        >>> c.sort_values(na_position='first')
        [NaN, NaN, 2, 2, 5]
        Categories (2, int64): [2, 5]
        >>> c.sort_values(ascending=False, na_position='first')
        [NaN, NaN, 5, 2, 2]
        Categories (2, int64): [2, 5]
        """
        inplace = validate_bool_kwarg(inplace, "inplace")
        if na_position not in ["last", "first"]:
            raise ValueError(f"invalid na_position: {repr(na_position)}")

        sorted_idx = nargsort(self, ascending=ascending, na_position=na_position)

        if inplace:
            self._codes = self._codes[sorted_idx]
        else:
            return self._constructor(
                values=self._codes[sorted_idx], dtype=self.dtype, fastpath=True
            )

    def _values_for_rank(self):
        """
        For correctly ranking ordered categorical data. See GH#15420

        Ordered categorical data should be ranked on the basis of
        codes with -1 translated to NaN.

        Returns
        -------
        numpy.array

        """
        from pandas import Series

        if self.ordered:
            values = self.codes
            mask = values == -1
            if mask.any():
                values = values.astype("float64")
                values[mask] = np.nan
        elif self.categories.is_numeric():
            values = np.array(self)
        else:
            #  reorder the categories (so rank can use the float codes)
            #  instead of passing an object array to rank
            values = np.array(
                self.rename_categories(Series(self.categories).rank().values)
            )
        return values

    def view(self, dtype=None):
        if dtype is not None:
            raise NotImplementedError(dtype)
        return self._constructor(values=self._codes, dtype=self.dtype, fastpath=True)

    def to_dense(self):
        """
        Return my 'dense' representation

        For internal compatibility with numpy arrays.

        Returns
        -------
        dense : array
        """
        warn(
            "Categorical.to_dense is deprecated and will be removed in "
            "a future version.  Use np.asarray(cat) instead.",
            FutureWarning,
            stacklevel=2,
        )
        return np.asarray(self)

    def fillna(self, value=None, method=None, limit=None):
        """
        Fill NA/NaN values using the specified method.

        Parameters
        ----------
        value : scalar, dict, Series
            If a scalar value is passed it is used to fill all missing values.
            Alternatively, a Series or dict can be used to fill in different
            values for each index. The value should not be a list. The
            value(s) passed should either be in the categories or should be
            NaN.
        method : {'backfill', 'bfill', 'pad', 'ffill', None}, default None
            Method to use for filling holes in reindexed Series
            pad / ffill: propagate last valid observation forward to next valid
            backfill / bfill: use NEXT valid observation to fill gap
        limit : int, default None
            (Not implemented yet for Categorical!)
            If method is specified, this is the maximum number of consecutive
            NaN values to forward/backward fill. In other words, if there is
            a gap with more than this number of consecutive NaNs, it will only
            be partially filled. If method is not specified, this is the
            maximum number of entries along the entire axis where NaNs will be
            filled.

        Returns
        -------
        filled : Categorical with NA/NaN filled
        """
        value, method = validate_fillna_kwargs(
            value, method, validate_scalar_dict_value=False
        )

        if value is None:
            value = np.nan
        if limit is not None:
            raise NotImplementedError(
                "specifying a limit for fillna has not been implemented yet"
            )

        codes = self._codes

        # pad / bfill
        if method is not None:

            # TODO: dispatch when self.categories is EA-dtype
            values = np.asarray(self).reshape(-1, len(self))
            values = interpolate_2d(values, method, 0, None, value).astype(
                self.categories.dtype
            )[0]
            codes = _get_codes_for_values(values, self.categories)

        else:

            # If value is a dict or a Series (a dict value has already
            # been converted to a Series)
            if isinstance(value, (np.ndarray, Categorical, ABCSeries)):
                # We get ndarray or Categorical if called via Series.fillna,
                #  where it will unwrap another aligned Series before getting here

                mask = ~algorithms.isin(value, self.categories)
                if not isna(value[mask]).all():
                    raise ValueError("fill value must be in categories")

                values_codes = _get_codes_for_values(value, self.categories)
                indexer = np.where(codes == -1)
                codes = codes.copy()
                codes[indexer] = values_codes[indexer]

            # If value is not a dict or Series it should be a scalar
            elif is_hashable(value):
                if not isna(value) and value not in self.categories:
                    raise ValueError("fill value must be in categories")

                mask = codes == -1
                if mask.any():
                    codes = codes.copy()
                    if isna(value):
                        codes[mask] = -1
                    else:
                        codes[mask] = self.categories.get_loc(value)

            else:
                raise TypeError(
                    f"'value' parameter must be a scalar, dict "
                    f"or Series, but you passed a {type(value).__name__}"
                )

        return self._constructor(codes, dtype=self.dtype, fastpath=True)

    def take(self, indexer, allow_fill: bool = False, fill_value=None):
        """
        Take elements from the Categorical.

        Parameters
        ----------
        indexer : sequence of int
            The indices in `self` to take. The meaning of negative values in
            `indexer` depends on the value of `allow_fill`.
        allow_fill : bool, default False
            How to handle negative values in `indexer`.

            * False: negative values in `indices` indicate positional indices
              from the right. This is similar to
              :func:`numpy.take`.

            * True: negative values in `indices` indicate missing values
              (the default). These values are set to `fill_value`. Any other
              other negative values raise a ``ValueError``.

            .. versionchanged:: 1.0.0

               Default value changed from ``True`` to ``False``.

        fill_value : object
            The value to use for `indices` that are missing (-1), when
            ``allow_fill=True``. This should be the category, i.e. a value
            in ``self.categories``, not a code.

        Returns
        -------
        Categorical
            This Categorical will have the same categories and ordered as
            `self`.

        See Also
        --------
        Series.take : Similar method for Series.
        numpy.ndarray.take : Similar method for NumPy arrays.

        Examples
        --------
        >>> cat = pd.Categorical(['a', 'a', 'b'])
        >>> cat
        [a, a, b]
        Categories (2, object): [a, b]

        Specify ``allow_fill==False`` to have negative indices mean indexing
        from the right.

        >>> cat.take([0, -1, -2], allow_fill=False)
        [a, b, a]
        Categories (2, object): [a, b]

        With ``allow_fill=True``, indices equal to ``-1`` mean "missing"
        values that should be filled with the `fill_value`, which is
        ``np.nan`` by default.

        >>> cat.take([0, -1, -1], allow_fill=True)
        [a, NaN, NaN]
        Categories (2, object): [a, b]

        The fill value can be specified.

        >>> cat.take([0, -1, -1], allow_fill=True, fill_value='a')
        [a, a, a]
        Categories (2, object): [a, b]

        Specifying a fill value that's not in ``self.categories``
        will raise a ``TypeError``.
        """
        indexer = np.asarray(indexer, dtype=np.intp)

        dtype = self.dtype

        if isna(fill_value):
            fill_value = -1
        elif allow_fill:
            # convert user-provided `fill_value` to codes
            if fill_value in self.categories:
                fill_value = self.categories.get_loc(fill_value)
            else:
                msg = (
                    f"'fill_value' ('{fill_value}') is not in this "
                    "Categorical's categories."
                )
                raise TypeError(msg)

        codes = take(self._codes, indexer, allow_fill=allow_fill, fill_value=fill_value)
        result = type(self).from_codes(codes, dtype=dtype)
        return result

    def take_nd(self, indexer, allow_fill: bool = False, fill_value=None):
        # GH#27745 deprecate alias that other EAs dont have
        warn(
            "Categorical.take_nd is deprecated, use Categorical.take instead",
            FutureWarning,
            stacklevel=2,
        )
        return self.take(indexer, allow_fill=allow_fill, fill_value=fill_value)

    def __len__(self) -> int:
        """
        The length of this Categorical.
        """
        return len(self._codes)

    def __iter__(self):
        """
        Returns an Iterator over the values of this Categorical.
        """
        return iter(self._internal_get_values().tolist())

    def __contains__(self, key) -> bool:
        """
        Returns True if `key` is in this Categorical.
        """
        # if key is a NaN, check if any NaN is in self.
        if is_scalar(key) and isna(key):
            return self.isna().any()

        return contains(self, key, container=self._codes)

    def _tidy_repr(self, max_vals=10, footer=True) -> str:
        """
        a short repr displaying only max_vals and an optional (but default
        footer)
        """
        num = max_vals // 2
        head = self[:num]._get_repr(length=False, footer=False)
        tail = self[-(max_vals - num) :]._get_repr(length=False, footer=False)

        result = f"{head[:-1]}, ..., {tail[1:]}"
        if footer:
            result = f"{result}\n{self._repr_footer()}"

        return str(result)

    def _repr_categories(self):
        """
        return the base repr for the categories
        """
        max_categories = (
            10
            if get_option("display.max_categories") == 0
            else get_option("display.max_categories")
        )
        from pandas.io.formats import format as fmt

        if len(self.categories) > max_categories:
            num = max_categories // 2
            head = fmt.format_array(self.categories[:num], None)
            tail = fmt.format_array(self.categories[-num:], None)
            category_strs = head + ["..."] + tail
        else:
            category_strs = fmt.format_array(self.categories, None)

        # Strip all leading spaces, which format_array adds for columns...
        category_strs = [x.strip() for x in category_strs]
        return category_strs

    def _repr_categories_info(self) -> str:
        """
        Returns a string representation of the footer.
        """
        category_strs = self._repr_categories()
        dtype = str(self.categories.dtype)
        levheader = f"Categories ({len(self.categories)}, {dtype}): "
        width, height = get_terminal_size()
        max_width = get_option("display.width") or width
        if console.in_ipython_frontend():
            # 0 = no breaks
            max_width = 0
        levstring = ""
        start = True
        cur_col_len = len(levheader)  # header
        sep_len, sep = (3, " < ") if self.ordered else (2, ", ")
        linesep = sep.rstrip() + "\n"  # remove whitespace
        for val in category_strs:
            if max_width != 0 and cur_col_len + sep_len + len(val) > max_width:
                levstring += linesep + (" " * (len(levheader) + 1))
                cur_col_len = len(levheader) + 1  # header + a whitespace
            elif not start:
                levstring += sep
                cur_col_len += len(val)
            levstring += val
            start = False
        # replace to simple save space by
        return levheader + "[" + levstring.replace(" < ... < ", " ... ") + "]"

    def _repr_footer(self) -> str:
        info = self._repr_categories_info()
        return f"Length: {len(self)}\n{info}"

    def _get_repr(self, length=True, na_rep="NaN", footer=True) -> str:
        from pandas.io.formats import format as fmt

        formatter = fmt.CategoricalFormatter(
            self, length=length, na_rep=na_rep, footer=footer
        )
        result = formatter.to_string()
        return str(result)

    def __repr__(self) -> str:
        """
        String representation.
        """
        _maxlen = 10
        if len(self._codes) > _maxlen:
            result = self._tidy_repr(_maxlen)
        elif len(self._codes) > 0:
            result = self._get_repr(length=len(self) > _maxlen)
        else:
            msg = self._get_repr(length=False, footer=True).replace("\n", ", ")
            result = f"[], {msg}"

        return result

    def _maybe_coerce_indexer(self, indexer):
        """
        return an indexer coerced to the codes dtype
        """
        if isinstance(indexer, np.ndarray) and indexer.dtype.kind == "i":
            indexer = indexer.astype(self._codes.dtype)
        return indexer

    def __getitem__(self, key):
        """
        Return an item.
        """
        if isinstance(key, (int, np.integer)):
            i = self._codes[key]
            if i == -1:
                return np.nan
            else:
                return self.categories[i]

        key = check_array_indexer(self, key)

        result = self._codes[key]
        if result.ndim > 1:
            deprecate_ndim_indexing(result)
            return result
        return self._constructor(result, dtype=self.dtype, fastpath=True)

    def __setitem__(self, key, value):
        """
        Item assignment.

        Raises
        ------
        ValueError
            If (one or more) Value is not in categories or if a assigned
            `Categorical` does not have the same categories
        """
        value = extract_array(value, extract_numpy=True)

        # require identical categories set
        if isinstance(value, Categorical):
            if not is_dtype_equal(self, value):
                raise ValueError(
                    "Cannot set a Categorical with another, "
                    "without identical categories"
                )
            if not self.categories.equals(value.categories):
                new_codes = recode_for_categories(
                    value.codes, value.categories, self.categories
                )
                value = Categorical.from_codes(new_codes, dtype=self.dtype)

        rvalue = value if is_list_like(value) else [value]

        from pandas import Index

        to_add = Index(rvalue).difference(self.categories)

        # no assignments of values not in categories, but it's always ok to set
        # something to np.nan
        if len(to_add) and not isna(to_add).all():
            raise ValueError(
                "Cannot setitem on a Categorical with a new "
                "category, set the categories first"
            )

        # set by position
        if isinstance(key, (int, np.integer)):
            pass

        # tuple of indexers (dataframe)
        elif isinstance(key, tuple):
            # only allow 1 dimensional slicing, but can
            # in a 2-d case be passd (slice(None),....)
            if len(key) == 2:
                if not com.is_null_slice(key[0]):
                    raise AssertionError("invalid slicing for a 1-ndim categorical")
                key = key[1]
            elif len(key) == 1:
                key = key[0]
            else:
                raise AssertionError("invalid slicing for a 1-ndim categorical")

        # slicing in Series or Categorical
        elif isinstance(key, slice):
            pass

        # else: array of True/False in Series or Categorical

        lindexer = self.categories.get_indexer(rvalue)
        lindexer = self._maybe_coerce_indexer(lindexer)

        key = check_array_indexer(self, key)
        self._codes[key] = lindexer

    def _reverse_indexer(self) -> Dict[Hashable, np.ndarray]:
        """
        Compute the inverse of a categorical, returning
        a dict of categories -> indexers.

        *This is an internal function*

        Returns
        -------
        dict of categories -> indexers

        Examples
        --------
        >>> c = pd.Categorical(list('aabca'))
        >>> c
        [a, a, b, c, a]
        Categories (3, object): [a, b, c]
        >>> c.categories
        Index(['a', 'b', 'c'], dtype='object')
        >>> c.codes
        array([0, 0, 1, 2, 0], dtype=int8)
        >>> c._reverse_indexer()
        {'a': array([0, 1, 4]), 'b': array([2]), 'c': array([3])}

        """
        categories = self.categories
        r, counts = libalgos.groupsort_indexer(
            self.codes.astype("int64"), categories.size
        )
        counts = counts.cumsum()
        _result = (r[start:end] for start, end in zip(counts, counts[1:]))
        result = dict(zip(categories, _result))
        return result

    # reduction ops #
    def _reduce(self, name, axis=0, **kwargs):
        func = getattr(self, name, None)
        if func is None:
            raise TypeError(f"Categorical cannot perform the operation {name}")
        return func(**kwargs)

    @deprecate_kwarg(old_arg_name="numeric_only", new_arg_name="skipna")
    def min(self, skipna=True):
        """
        The minimum value of the object.

        Only ordered `Categoricals` have a minimum!

        .. versionchanged:: 1.0.0

           Returns an NA value on empty arrays

        Raises
        ------
        TypeError
            If the `Categorical` is not `ordered`.

        Returns
        -------
        min : the minimum of this `Categorical`
        """
        self.check_for_ordered("min")

        if not len(self._codes):
            return self.dtype.na_value

        good = self._codes != -1
        if not good.all():
            if skipna:
                pointer = self._codes[good].min()
            else:
                return np.nan
        else:
            pointer = self._codes.min()
        return self.categories[pointer]

    @deprecate_kwarg(old_arg_name="numeric_only", new_arg_name="skipna")
    def max(self, skipna=True):
        """
        The maximum value of the object.

        Only ordered `Categoricals` have a maximum!

        .. versionchanged:: 1.0.0

           Returns an NA value on empty arrays

        Raises
        ------
        TypeError
            If the `Categorical` is not `ordered`.

        Returns
        -------
        max : the maximum of this `Categorical`
        """
        self.check_for_ordered("max")

        if not len(self._codes):
            return self.dtype.na_value

        good = self._codes != -1
        if not good.all():
            if skipna:
                pointer = self._codes[good].max()
            else:
                return np.nan
        else:
            pointer = self._codes.max()
        return self.categories[pointer]

    def mode(self, dropna=True):
        """
        Returns the mode(s) of the Categorical.

        Always returns `Categorical` even if only one value.

        Parameters
        ----------
        dropna : bool, default True
            Don't consider counts of NaN/NaT.

            .. versionadded:: 0.24.0

        Returns
        -------
        modes : `Categorical` (sorted)
        """
        codes = self._codes
        if dropna:
            good = self._codes != -1
            codes = self._codes[good]
        codes = sorted(htable.mode_int64(ensure_int64(codes), dropna))
        return self._constructor(values=codes, dtype=self.dtype, fastpath=True)

    def unique(self):
        """
        Return the ``Categorical`` which ``categories`` and ``codes`` are
        unique. Unused categories are NOT returned.

        - unordered category: values and categories are sorted by appearance
          order.
        - ordered category: values are sorted by appearance order, categories
          keeps existing order.

        Returns
        -------
        unique values : ``Categorical``

        See Also
        --------
        pandas.unique
        CategoricalIndex.unique
        Series.unique

        Examples
        --------
        An unordered Categorical will return categories in the
        order of appearance.

        >>> pd.Categorical(list("baabc")).unique()
        [b, a, c]
        Categories (3, object): [b, a, c]

        >>> pd.Categorical(list("baabc"), categories=list("abc")).unique()
        [b, a, c]
        Categories (3, object): [b, a, c]

        An ordered Categorical preserves the category ordering.

        >>> pd.Categorical(
        ...     list("baabc"), categories=list("abc"), ordered=True
        ... ).unique()
        [b, a, c]
        Categories (3, object): [a < b < c]
        """
        # unlike np.unique, unique1d does not sort
        unique_codes = unique1d(self.codes)
        cat = self.copy()

        # keep nan in codes
        cat._codes = unique_codes

        # exclude nan from indexer for categories
        take_codes = unique_codes[unique_codes != -1]
        if self.ordered:
            take_codes = np.sort(take_codes)
        return cat.set_categories(cat.categories.take(take_codes))

    def _values_for_factorize(self):
        codes = self.codes.astype("int64")
        return codes, -1

    @classmethod
    def _from_factorized(cls, uniques, original):
        return original._constructor(
            original.categories.take(uniques), dtype=original.dtype
        )

    def equals(self, other):
        """
        Returns True if categorical arrays are equal.

        Parameters
        ----------
        other : `Categorical`

        Returns
        -------
        bool
        """
        if self.is_dtype_equal(other):
            if self.categories.equals(other.categories):
                # fastpath to avoid re-coding
                other_codes = other._codes
            else:
                other_codes = recode_for_categories(
                    other.codes, other.categories, self.categories
                )
            return np.array_equal(self._codes, other_codes)
        return False

    def is_dtype_equal(self, other):
        """
        Returns True if categoricals are the same dtype
          same categories, and same ordered

        Parameters
        ----------
        other : Categorical

        Returns
        -------
        bool
        """
        try:
            return hash(self.dtype) == hash(other.dtype)
        except (AttributeError, TypeError):
            return False

    def describe(self):
        """
        Describes this Categorical

        Returns
        -------
        description: `DataFrame`
            A dataframe with frequency and counts by category.
        """
        counts = self.value_counts(dropna=False)
        freqs = counts / float(counts.sum())

        from pandas.core.reshape.concat import concat

        result = concat([counts, freqs], axis=1)
        result.columns = ["counts", "freqs"]
        result.index.name = "categories"

        return result

    @Substitution(klass="Categorical")
    @Appender(_extension_array_shared_docs["repeat"])
    def repeat(self, repeats, axis=None):
        nv.validate_repeat(tuple(), dict(axis=axis))
        codes = self._codes.repeat(repeats)
        return self._constructor(values=codes, dtype=self.dtype, fastpath=True)

    # Implement the ExtensionArray interface
    @property
    def _can_hold_na(self):
        return True

    @classmethod
    def _concat_same_type(self, to_concat):
        from pandas.core.dtypes.concat import concat_categorical

        return concat_categorical(to_concat)

    def isin(self, values):
        """
        Check whether `values` are contained in Categorical.

        Return a boolean NumPy Array showing whether each element in
        the Categorical matches an element in the passed sequence of
        `values` exactly.

        Parameters
        ----------
        values : set or list-like
            The sequence of values to test. Passing in a single string will
            raise a ``TypeError``. Instead, turn a single string into a
            list of one element.

        Returns
        -------
        isin : numpy.ndarray (bool dtype)

        Raises
        ------
        TypeError
          * If `values` is not a set or list-like

        See Also
        --------
        pandas.Series.isin : Equivalent method on Series.

        Examples
        --------
        >>> s = pd.Categorical(['lama', 'cow', 'lama', 'beetle', 'lama',
        ...                'hippo'])
        >>> s.isin(['cow', 'lama'])
        array([ True,  True,  True, False,  True, False])

        Passing a single string as ``s.isin('lama')`` will raise an error. Use
        a list of one element instead:

        >>> s.isin(['lama'])
        array([ True, False,  True, False,  True, False])
        """
        if not is_list_like(values):
            values_type = type(values).__name__
            raise TypeError(
                "only list-like objects are allowed to be passed "
                f"to isin(), you passed a [{values_type}]"
            )
        values = sanitize_array(values, None, None)
        null_mask = np.asarray(isna(values))
        code_values = self.categories.get_indexer(values)
        code_values = code_values[null_mask | (code_values >= 0)]
        return algorithms.isin(self.codes, code_values)

    def replace(self, to_replace, value, inplace: bool = False):
        """
        Replaces all instances of one value with another

        Parameters
        ----------
        to_replace: object
            The value to be replaced

        value: object
            The value to replace it with

        inplace: bool
            Whether the operation is done in-place

        Returns
        -------
        None if inplace is True, otherwise the new Categorical after replacement


        Examples
        --------
        >>> s = pd.Categorical([1, 2, 1, 3])
        >>> s.replace(1, 3)
        [3, 2, 3, 3]
        Categories (2, int64): [2, 3]
        """
        inplace = validate_bool_kwarg(inplace, "inplace")
        cat = self if inplace else self.copy()

        # build a dict of (to replace -> value) pairs
        if is_list_like(to_replace):
            # if to_replace is list-like and value is scalar
            replace_dict = {replace_value: value for replace_value in to_replace}
        else:
            # if both to_replace and value are scalar
            replace_dict = {to_replace: value}

        # other cases, like if both to_replace and value are list-like or if
        # to_replace is a dict, are handled separately in NDFrame
        for replace_value, new_value in replace_dict.items():
            if replace_value in cat.categories:
                if isna(new_value):
                    cat.remove_categories(replace_value, inplace=True)
                    continue
                categories = cat.categories.tolist()
                index = categories.index(replace_value)
                if new_value in cat.categories:
                    value_index = categories.index(new_value)
                    cat._codes[cat._codes == index] = value_index
                    cat.remove_categories(replace_value, inplace=True)
                else:
                    categories[index] = new_value
                    cat.rename_categories(categories, inplace=True)
        if not inplace:
            return cat


# The Series.cat accessor


@delegate_names(
    delegate=Categorical, accessors=["categories", "ordered"], typ="property"
)
@delegate_names(
    delegate=Categorical,
    accessors=[
        "rename_categories",
        "reorder_categories",
        "add_categories",
        "remove_categories",
        "remove_unused_categories",
        "set_categories",
        "as_ordered",
        "as_unordered",
    ],
    typ="method",
)
class CategoricalAccessor(PandasDelegate, PandasObject, NoNewAttributesMixin):
    """
    Accessor object for categorical properties of the Series values.

    Be aware that assigning to `categories` is a inplace operation, while all
    methods return new categorical data per default (but can be called with
    `inplace=True`).

    Parameters
    ----------
    data : Series or CategoricalIndex

    Examples
    --------
    >>> s = pd.Series(list("abbccc")).astype("category")
    >>> s
    0    a
    1    b
    2    b
    3    c
    4    c
    5    c
    dtype: category
    Categories (3, object): [a, b, c]

    >>> s.cat.categories
    Index(['a', 'b', 'c'], dtype='object')

    >>> s.cat.rename_categories(list("cba"))
    0    c
    1    b
    2    b
    3    a
    4    a
    5    a
    dtype: category
    Categories (3, object): [c, b, a]

    >>> s.cat.reorder_categories(list("cba"))
    0    a
    1    b
    2    b
    3    c
    4    c
    5    c
    dtype: category
    Categories (3, object): [c, b, a]

    >>> s.cat.add_categories(["d", "e"])
    0    a
    1    b
    2    b
    3    c
    4    c
    5    c
    dtype: category
    Categories (5, object): [a, b, c, d, e]

    >>> s.cat.remove_categories(["a", "c"])
    0    NaN
    1      b
    2      b
    3    NaN
    4    NaN
    5    NaN
    dtype: category
    Categories (1, object): [b]

    >>> s1 = s.cat.add_categories(["d", "e"])
    >>> s1.cat.remove_unused_categories()
    0    a
    1    b
    2    b
    3    c
    4    c
    5    c
    dtype: category
    Categories (3, object): [a, b, c]

    >>> s.cat.set_categories(list("abcde"))
    0    a
    1    b
    2    b
    3    c
    4    c
    5    c
    dtype: category
    Categories (5, object): [a, b, c, d, e]

    >>> s.cat.as_ordered()
    0    a
    1    b
    2    b
    3    c
    4    c
    5    c
    dtype: category
    Categories (3, object): [a < b < c]

    >>> s.cat.as_unordered()
    0    a
    1    b
    2    b
    3    c
    4    c
    5    c
    dtype: category
    Categories (3, object): [a, b, c]
    """

    def __init__(self, data):
        self._validate(data)
        self._parent = data.values
        self._index = data.index
        self._name = data.name
        self._freeze()

    @staticmethod
    def _validate(data):
        if not is_categorical_dtype(data.dtype):
            raise AttributeError("Can only use .cat accessor with a 'category' dtype")

    def _delegate_property_get(self, name):
        return getattr(self._parent, name)

    def _delegate_property_set(self, name, new_values):
        return setattr(self._parent, name, new_values)

    @property
    def codes(self):
        """
        Return Series of codes as well as the index.
        """
        from pandas import Series

        return Series(self._parent.codes, index=self._index)

    def _delegate_method(self, name, *args, **kwargs):
        from pandas import Series

        method = getattr(self._parent, name)
        res = method(*args, **kwargs)
        if res is not None:
            return Series(res, index=self._index, name=self._name)


# utility routines


def _get_codes_for_values(values, categories):
    """
    utility routine to turn values into codes given the specified categories
    """
    dtype_equal = is_dtype_equal(values.dtype, categories.dtype)

    if is_extension_array_dtype(categories.dtype) and is_object_dtype(values):
        # Support inferring the correct extension dtype from an array of
        # scalar objects. e.g.
        # Categorical(array[Period, Period], categories=PeriodIndex(...))
        cls = categories.dtype.construct_array_type()
        values = maybe_cast_to_extension_array(cls, values)
        if not isinstance(values, cls):
            # exception raised in _from_sequence
            values = ensure_object(values)
            categories = ensure_object(categories)
    elif not dtype_equal:
        values = ensure_object(values)
        categories = ensure_object(categories)

    hash_klass, vals = _get_data_algo(values)
    _, cats = _get_data_algo(categories)
    t = hash_klass(len(cats))
    t.map_locations(cats)
    return coerce_indexer_dtype(t.lookup(vals), cats)


def recode_for_categories(codes: np.ndarray, old_categories, new_categories):
    """
    Convert a set of codes for to a new set of categories

    Parameters
    ----------
    codes : np.ndarray
    old_categories, new_categories : Index

    Returns
    -------
    new_codes : np.ndarray[np.int64]

    Examples
    --------
    >>> old_cat = pd.Index(['b', 'a', 'c'])
    >>> new_cat = pd.Index(['a', 'b'])
    >>> codes = np.array([0, 1, 1, 2])
    >>> recode_for_categories(codes, old_cat, new_cat)
    array([ 1,  0,  0, -1], dtype=int8)
    """
    if len(old_categories) == 0:
        # All null anyway, so just retain the nulls
        return codes.copy()
    elif new_categories.equals(old_categories):
        # Same categories, so no need to actually recode
        return codes.copy()
    indexer = coerce_indexer_dtype(
        new_categories.get_indexer(old_categories), new_categories
    )
    new_codes = take_1d(indexer, codes.copy(), fill_value=-1)
    return new_codes


def _convert_to_list_like(list_like):
    if hasattr(list_like, "dtype"):
        return list_like
    if isinstance(list_like, list):
        return list_like
    if is_sequence(list_like) or isinstance(list_like, tuple) or is_iterator(list_like):
        return list(list_like)
    elif is_scalar(list_like):
        return [list_like]
    else:
        # TODO: is this reached?
        return [list_like]


def factorize_from_iterable(values):
    """
    Factorize an input `values` into `categories` and `codes`. Preserves
    categorical dtype in `categories`.

    *This is an internal function*

    Parameters
    ----------
    values : list-like

    Returns
    -------
    codes : ndarray
    categories : Index
        If `values` has a categorical dtype, then `categories` is
        a CategoricalIndex keeping the categories and order of `values`.
    """
    if not is_list_like(values):
        raise TypeError("Input must be list-like")

    if is_categorical_dtype(values):
        values = extract_array(values)
        # The Categorical we want to build has the same categories
        # as values but its codes are by def [0, ..., len(n_categories) - 1]
        cat_codes = np.arange(len(values.categories), dtype=values.codes.dtype)
        categories = Categorical.from_codes(cat_codes, dtype=values.dtype)
        codes = values.codes
    else:
        # The value of ordered is irrelevant since we don't use cat as such,
        # but only the resulting categories, the order of which is independent
        # from ordered. Set ordered to False as default. See GH #15457
        cat = Categorical(values, ordered=False)
        categories = cat.categories
        codes = cat.codes
    return codes, categories


def factorize_from_iterables(iterables):
    """
    A higher-level wrapper over `factorize_from_iterable`.

    *This is an internal function*

    Parameters
    ----------
    iterables : list-like of list-likes

    Returns
    -------
    codes_list : list of ndarrays
    categories_list : list of Indexes

    Notes
    -----
    See `factorize_from_iterable` for more info.
    """
    if len(iterables) == 0:
        # For consistency, it should return a list of 2 lists.
        return [[], []]
    return map(list, zip(*(factorize_from_iterable(it) for it in iterables)))
