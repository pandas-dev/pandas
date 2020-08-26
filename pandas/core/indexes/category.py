from typing import Any, List
import warnings

import numpy as np

from pandas._config import get_option

from pandas._libs import index as libindex
from pandas._libs.hashtable import duplicated_int64
from pandas._libs.lib import no_default
from pandas._typing import Label
from pandas.util._decorators import Appender, cache_readonly, doc

from pandas.core.dtypes.common import (
    ensure_platform_int,
    is_categorical_dtype,
    is_interval_dtype,
    is_list_like,
    is_scalar,
    pandas_dtype,
)
from pandas.core.dtypes.dtypes import CategoricalDtype
from pandas.core.dtypes.missing import is_valid_nat_for_dtype, isna, notna

from pandas.core import accessor
from pandas.core.algorithms import take_1d
from pandas.core.arrays.categorical import Categorical, contains, recode_for_categories
import pandas.core.common as com
from pandas.core.construction import extract_array
import pandas.core.indexes.base as ibase
from pandas.core.indexes.base import Index, _index_shared_docs, maybe_extract_name
from pandas.core.indexes.extension import ExtensionIndex, inherit_names
import pandas.core.missing as missing
from pandas.core.ops import get_op_result_name

_index_doc_kwargs = dict(ibase._index_doc_kwargs)
_index_doc_kwargs.update(dict(target_klass="CategoricalIndex"))


@inherit_names(
    [
        "argsort",
        "_internal_get_values",
        "tolist",
        "codes",
        "categories",
        "ordered",
        "_reverse_indexer",
        "searchsorted",
        "is_dtype_equal",
        "min",
        "max",
    ],
    Categorical,
)
@accessor.delegate_names(
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
    overwrite=True,
)
class CategoricalIndex(ExtensionIndex, accessor.PandasDelegate):
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
    <https://pandas.pydata.org/pandas-docs/stable/user_guide/advanced.html#categoricalindex>`_
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

    codes: np.ndarray
    categories: Index
    _data: Categorical

    @property
    def _engine_type(self):
        # self.codes can have dtype int8, int16, int32 or int64, so we need
        # to return the corresponding engine type (libindex.Int8Engine, etc.).
        return {
            np.int8: libindex.Int8Engine,
            np.int16: libindex.Int16Engine,
            np.int32: libindex.Int32Engine,
            np.int64: libindex.Int64Engine,
        }[self.codes.dtype.type]

    _attributes = ["name"]

    # --------------------------------------------------------------------
    # Constructors

    def __new__(
        cls, data=None, categories=None, ordered=None, dtype=None, copy=False, name=None
    ):

        dtype = CategoricalDtype._from_values_or_dtype(data, categories, ordered, dtype)

        name = maybe_extract_name(name, data, cls)

        if not is_categorical_dtype(data):
            # don't allow scalars
            # if data is None, then categories must be provided
            if is_scalar(data):
                if data is not None or categories is None:
                    raise cls._scalar_data_error(data)
                data = []

        assert isinstance(dtype, CategoricalDtype), dtype
        data = extract_array(data, extract_numpy=True)

        if not isinstance(data, Categorical):
            data = Categorical(data, dtype=dtype)
        elif isinstance(dtype, CategoricalDtype) and dtype != data.dtype:
            # we want to silently ignore dtype='category'
            data = data._set_dtype(dtype)

        data = data.copy() if copy else data

        return cls._simple_new(data, name=name)

    def _create_from_codes(self, codes, dtype=None, name=None):
        """
        *this is an internal non-public method*

        create the correct categorical from codes

        Parameters
        ----------
        codes : new codes
        dtype: CategoricalDtype, defaults to existing
        name : optional name attribute, defaults to existing

        Returns
        -------
        CategoricalIndex
        """
        if dtype is None:
            dtype = self.dtype
        if name is None:
            name = self.name
        cat = Categorical.from_codes(codes, dtype=dtype)
        return CategoricalIndex(cat, name=name)

    @classmethod
    def _simple_new(cls, values: Categorical, name: Label = None):
        assert isinstance(values, Categorical), type(values)
        result = object.__new__(cls)

        result._data = values
        result.name = name
        result._cache = {}

        result._reset_identity()
        result._no_setting_name = False
        return result

    # --------------------------------------------------------------------

    @doc(Index._shallow_copy)
    def _shallow_copy(self, values=None, name: Label = no_default):
        name = self.name if name is no_default else name

        if values is not None:
            values = Categorical(values, dtype=self.dtype)

        return super()._shallow_copy(values=values, name=name)

    def _is_dtype_compat(self, other) -> bool:
        """
        *this is an internal non-public method*

        provide a comparison between the dtype of self and other (coercing if
        needed)

        Raises
        ------
        TypeError if the dtypes are not compatible
        """
        if is_categorical_dtype(other):
            if isinstance(other, CategoricalIndex):
                other = other._values
            if not other.is_dtype_equal(self):
                raise TypeError(
                    "categories must match existing categories when appending"
                )
        else:
            values = other
            if not is_list_like(values):
                values = [values]
            cat = Categorical(other, dtype=self.dtype)
            other = CategoricalIndex(cat)
            if not other.isin(values).all():
                raise TypeError(
                    "cannot append a non-category item to a CategoricalIndex"
                )

        return other

    def equals(self, other) -> bool:
        """
        Determine if two CategoricalIndex objects contain the same elements.

        Returns
        -------
        bool
            If two CategoricalIndex objects have equal elements True,
            otherwise False.
        """
        if self.is_(other):
            return True

        if not isinstance(other, Index):
            return False

        try:
            other = self._is_dtype_compat(other)
            if isinstance(other, type(self)):
                other = other._data
            return self._data.equals(other)
        except (TypeError, ValueError):
            pass

        return False

    # --------------------------------------------------------------------
    # Rendering Methods

    @property
    def _formatter_func(self):
        return self.categories._formatter_func

    def _format_attrs(self):
        """
        Return a list of tuples of the (attr,formatted_value)
        """
        max_categories = (
            10
            if get_option("display.max_categories") == 0
            else get_option("display.max_categories")
        )
        attrs = [
            (
                "categories",
                ibase.default_pprint(self.categories, max_seq_items=max_categories),
            ),
            ("ordered", self.ordered),
        ]
        if self.name is not None:
            attrs.append(("name", ibase.default_pprint(self.name)))
        attrs.append(("dtype", f"'{self.dtype.name}'"))
        max_seq_items = get_option("display.max_seq_items") or len(self)
        if len(self) > max_seq_items:
            attrs.append(("length", len(self)))
        return attrs

    def _format_with_header(self, header: List[str], na_rep: str = "NaN") -> List[str]:
        from pandas.io.formats.printing import pprint_thing

        result = [
            pprint_thing(x, escape_chars=("\t", "\r", "\n")) if notna(x) else na_rep
            for x in self._values
        ]
        return header + result

    # --------------------------------------------------------------------

    @property
    def inferred_type(self) -> str:
        return "categorical"

    @property
    def values(self):
        """ return the underlying data, which is a Categorical """
        return self._data

    @property
    def _has_complex_internals(self) -> bool:
        # used to avoid libreduction code paths, which raise or require conversion
        return True

    @doc(Index.__contains__)
    def __contains__(self, key: Any) -> bool:
        # if key is a NaN, check if any NaN is in self.
        if is_valid_nat_for_dtype(key, self.categories.dtype):
            return self.hasnans

        return contains(self, key, container=self._engine)

    @doc(Index.astype)
    def astype(self, dtype, copy=True):
        if dtype is not None:
            dtype = pandas_dtype(dtype)

        if is_interval_dtype(dtype):
            from pandas import IntervalIndex

            return IntervalIndex(np.array(self))
        elif is_categorical_dtype(dtype):
            # GH 18630
            dtype = self.dtype.update_dtype(dtype)
            if dtype == self.dtype:
                return self.copy() if copy else self

        return Index.astype(self, dtype=dtype, copy=copy)

    @cache_readonly
    def _isnan(self):
        """ return if each value is nan"""
        return self._data.codes == -1

    @doc(Index.fillna)
    def fillna(self, value, downcast=None):
        self._assert_can_do_op(value)
        return CategoricalIndex(self._data.fillna(value), name=self.name)

    @cache_readonly
    def _engine(self):
        # we are going to look things up with the codes themselves.
        # To avoid a reference cycle, bind `codes` to a local variable, so
        # `self` is not passed into the lambda.
        codes = self.codes
        return self._engine_type(lambda: codes, len(self))

    @doc(Index.unique)
    def unique(self, level=None):
        if level is not None:
            self._validate_index_level(level)
        result = self._values.unique()
        # Use _simple_new instead of _shallow_copy to ensure we keep dtype
        #  of result, not self.
        return type(self)._simple_new(result, name=self.name)

    @doc(Index.duplicated)
    def duplicated(self, keep="first"):
        codes = self.codes.astype("i8")
        return duplicated_int64(codes, keep)

    def _to_safe_for_reshape(self):
        """ convert to object if we are a categorical """
        return self.astype("object")

    def _maybe_cast_indexer(self, key):
        code = self.categories.get_loc(key)
        code = self.codes.dtype.type(code)
        return code

    @doc(Index.where)
    def where(self, cond, other=None):
        # TODO: Investigate an alternative implementation with
        # 1. copy the underlying Categorical
        # 2. setitem with `cond` and `other`
        # 3. Rebuild CategoricalIndex.
        if other is None:
            other = self._na_value
        values = np.where(cond, self._values, other)
        cat = Categorical(values, dtype=self.dtype)
        return type(self)._simple_new(cat, name=self.name)

    def reindex(self, target, method=None, level=None, limit=None, tolerance=None):
        """
        Create index with target's values (move/add/delete values as necessary)

        Returns
        -------
        new_index : pd.Index
            Resulting index
        indexer : np.ndarray or None
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

        target = ibase.ensure_index(target)

        missing: List[int]
        if self.equals(target):
            indexer = None
            missing = []
        else:
            indexer, missing = self.get_indexer_non_unique(np.array(target))

        if len(self.codes) and indexer is not None:
            new_target = self.take(indexer)
        else:
            new_target = target

        # filling in missing if needed
        if len(missing):
            cats = self.categories.get_indexer(target)

            if (cats == -1).any():
                # coerce to a regular index here!
                result = Index(np.array(self), name=self.name)
                new_target, indexer, _ = result._reindex_non_unique(np.array(target))
            else:

                codes = new_target.codes.copy()
                codes[indexer == -1] = cats[missing]
                new_target = self._create_from_codes(codes)

        # we always want to return an Index type here
        # to be consistent with .reindex for other index types (e.g. they don't
        # coerce based on the actual values, only on the dtype)
        # unless we had an initial Categorical to begin with
        # in which case we are going to conform to the passed Categorical
        new_target = np.asarray(new_target)
        if is_categorical_dtype(target):
            new_target = target._shallow_copy(new_target, name=self.name)
        else:
            new_target = Index(new_target, name=self.name)

        return new_target, indexer

    def _reindex_non_unique(self, target):
        """
        reindex from a non-unique; which CategoricalIndex's are almost
        always
        """
        new_target, indexer = self.reindex(target)
        new_indexer = None

        check = indexer == -1
        if check.any():
            new_indexer = np.arange(len(self.take(indexer)))
            new_indexer[check] = -1

        cats = self.categories.get_indexer(target)
        if not (cats == -1).any():
            # .reindex returns normal Index. Revert to CategoricalIndex if
            # all targets are included in my categories
            new_target = self._shallow_copy(new_target)

        return new_target, indexer, new_indexer

    @Appender(_index_shared_docs["get_indexer"] % _index_doc_kwargs)
    def get_indexer(self, target, method=None, limit=None, tolerance=None):
        method = missing.clean_reindex_fill_method(method)
        target = ibase.ensure_index(target)

        if self.is_unique and self.equals(target):
            return np.arange(len(self), dtype="intp")

        if method == "pad" or method == "backfill":
            raise NotImplementedError(
                "method='pad' and method='backfill' not "
                "implemented yet for CategoricalIndex"
            )
        elif method == "nearest":
            raise NotImplementedError(
                "method='nearest' not implemented yet for CategoricalIndex"
            )

        if isinstance(target, CategoricalIndex) and self._values.is_dtype_equal(target):
            if self._values.equals(target._values):
                # we have the same codes
                codes = target.codes
            else:
                codes = recode_for_categories(
                    target.codes, target.categories, self._values.categories
                )
        else:
            if isinstance(target, CategoricalIndex):
                code_indexer = self.categories.get_indexer(target.categories)
                codes = take_1d(code_indexer, target.codes, fill_value=-1)
            else:
                codes = self.categories.get_indexer(target)

        indexer, _ = self._engine.get_indexer_non_unique(codes)
        return ensure_platform_int(indexer)

    @Appender(_index_shared_docs["get_indexer_non_unique"] % _index_doc_kwargs)
    def get_indexer_non_unique(self, target):
        target = ibase.ensure_index(target)

        if isinstance(target, CategoricalIndex):
            # Indexing on codes is more efficient if categories are the same:
            if target.categories is self.categories:
                target = target.codes
                indexer, missing = self._engine.get_indexer_non_unique(target)
                return ensure_platform_int(indexer), missing
            target = target._values

        codes = self.categories.get_indexer(target)
        indexer, missing = self._engine.get_indexer_non_unique(codes)
        return ensure_platform_int(indexer), missing

    @doc(Index._convert_list_indexer)
    def _convert_list_indexer(self, keyarr):
        # Return our indexer or raise if all of the values are not included in
        # the categories

        if self.categories._defer_to_indexing:
            indexer = self.categories._convert_list_indexer(keyarr)
            return Index(self.codes).get_indexer_for(indexer)

        indexer = self.categories.get_indexer(np.asarray(keyarr))
        if (indexer == -1).any():
            raise KeyError(
                "a list-indexer must only include values that are in the categories"
            )

        return self.get_indexer(keyarr)

    @doc(Index._convert_arr_indexer)
    def _convert_arr_indexer(self, keyarr):
        keyarr = com.asarray_tuplesafe(keyarr)

        if self.categories._defer_to_indexing:
            return keyarr

        return self._shallow_copy(keyarr)

    @doc(Index._convert_index_indexer)
    def _convert_index_indexer(self, keyarr):
        return self._shallow_copy(keyarr)

    def take_nd(self, *args, **kwargs):
        """Alias for `take`"""
        warnings.warn(
            "CategoricalIndex.take_nd is deprecated, use CategoricalIndex.take instead",
            FutureWarning,
            stacklevel=2,
        )
        return self.take(*args, **kwargs)

    @doc(Index._maybe_cast_slice_bound)
    def _maybe_cast_slice_bound(self, label, side, kind):
        if kind == "loc":
            return label

        return super()._maybe_cast_slice_bound(label, side, kind)

    def map(self, mapper):
        """
        Map values using input correspondence (a dict, Series, or function).

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
        >>> idx = pd.CategoricalIndex(['a', 'b', 'c'])
        >>> idx
        CategoricalIndex(['a', 'b', 'c'], categories=['a', 'b', 'c'],
                          ordered=False, dtype='category')
        >>> idx.map(lambda x: x.upper())
        CategoricalIndex(['A', 'B', 'C'], categories=['A', 'B', 'C'],
                         ordered=False, dtype='category')
        >>> idx.map({'a': 'first', 'b': 'second', 'c': 'third'})
        CategoricalIndex(['first', 'second', 'third'], categories=['first',
                         'second', 'third'], ordered=False, dtype='category')

        If the mapping is one-to-one the ordering of the categories is
        preserved:

        >>> idx = pd.CategoricalIndex(['a', 'b', 'c'], ordered=True)
        >>> idx
        CategoricalIndex(['a', 'b', 'c'], categories=['a', 'b', 'c'],
                         ordered=True, dtype='category')
        >>> idx.map({'a': 3, 'b': 2, 'c': 1})
        CategoricalIndex([3, 2, 1], categories=[3, 2, 1], ordered=True,
                         dtype='category')

        If the mapping is not one-to-one an :class:`~pandas.Index` is returned:

        >>> idx.map({'a': 'first', 'b': 'second', 'c': 'first'})
        Index(['first', 'second', 'first'], dtype='object')

        If a `dict` is used, all unmapped categories are mapped to `NaN` and
        the result is an :class:`~pandas.Index`:

        >>> idx.map({'a': 'first', 'b': 'second'})
        Index(['first', 'second', nan], dtype='object')
        """
        mapped = self._values.map(mapper)
        return Index(mapped, name=self.name)

    def delete(self, loc):
        """
        Make new Index with passed location(-s) deleted

        Returns
        -------
        new_index : Index
        """
        return self._create_from_codes(np.delete(self.codes, loc))

    def insert(self, loc: int, item):
        """
        Make new Index inserting new item at location. Follows
        Python list.append semantics for negative values

        Parameters
        ----------
        loc : int
        item : object

        Returns
        -------
        new_index : Index

        Raises
        ------
        ValueError if the item is not in the categories

        """
        code = self.categories.get_indexer([item])
        if (code == -1) and not (is_scalar(item) and isna(item)):
            raise TypeError(
                "cannot insert an item into a CategoricalIndex "
                "that is not already an existing category"
            )

        codes = self.codes
        codes = np.concatenate((codes[:loc], code, codes[loc:]))
        return self._create_from_codes(codes)

    def _concat(self, to_concat, name):
        # if calling index is category, don't check dtype of others
        codes = np.concatenate([self._is_dtype_compat(c).codes for c in to_concat])
        result = self._create_from_codes(codes, name=name)
        # if name is None, _create_from_codes sets self.name
        result.name = name
        return result

    def _delegate_method(self, name: str, *args, **kwargs):
        """ method delegation to the ._values """
        method = getattr(self._values, name)
        if "inplace" in kwargs:
            raise ValueError("cannot use inplace with CategoricalIndex")
        res = method(*args, **kwargs)
        if is_scalar(res):
            return res
        return CategoricalIndex(res, name=self.name)

    def _wrap_joined_index(
        self, joined: np.ndarray, other: "CategoricalIndex"
    ) -> "CategoricalIndex":
        name = get_op_result_name(self, other)
        return self._create_from_codes(joined, name=name)


CategoricalIndex._add_numeric_methods_add_sub_disabled()
CategoricalIndex._add_numeric_methods_disabled()
CategoricalIndex._add_logical_methods_disabled()
