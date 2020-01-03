import operator
from typing import Any, List

import numpy as np

from pandas._config import get_option

from pandas._libs import index as libindex
from pandas._libs.hashtable import duplicated_int64
from pandas._typing import AnyArrayLike
import pandas.compat as compat
from pandas.compat.numpy import function as nv
from pandas.util._decorators import Appender, cache_readonly

from pandas.core.dtypes.common import (
    ensure_platform_int,
    is_categorical_dtype,
    is_interval_dtype,
    is_list_like,
    is_scalar,
)
from pandas.core.dtypes.dtypes import CategoricalDtype
from pandas.core.dtypes.generic import ABCCategorical, ABCSeries
from pandas.core.dtypes.missing import isna

from pandas.core import accessor
from pandas.core.algorithms import take_1d
from pandas.core.arrays.categorical import Categorical, _recode_for_categories, contains
import pandas.core.common as com
import pandas.core.indexes.base as ibase
from pandas.core.indexes.base import Index, _index_shared_docs, maybe_extract_name
import pandas.core.missing as missing
from pandas.core.ops import get_op_result_name

_index_doc_kwargs = dict(ibase._index_doc_kwargs)
_index_doc_kwargs.update(dict(target_klass="CategoricalIndex"))


@accessor.delegate_names(
    delegate=Categorical,
    accessors=["codes", "categories", "ordered"],
    typ="property",
    overwrite=True,
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
        "min",
        "max",
        "is_dtype_equal",
        "tolist",
        "_internal_get_values",
        "_reverse_indexer",
        "searchsorted",
        "argsort",
    ],
    typ="method",
    overwrite=True,
)
class CategoricalIndex(Index, accessor.PandasDelegate):
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

        .. versionadded:: 0.21.0
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
    <http://pandas.pydata.org/pandas-docs/stable/user_guide/advanced.html#categoricalindex>`_
    for more.

    Examples
    --------
    >>> pd.CategoricalIndex(['a', 'b', 'c', 'a', 'b', 'c'])
    CategoricalIndex(['a', 'b', 'c', 'a', 'b', 'c'], categories=['a', 'b', 'c'], ordered=False, dtype='category')  # noqa

    ``CategoricalIndex`` can also be instantiated from a ``Categorical``:

    >>> c = pd.Categorical(['a', 'b', 'c', 'a', 'b', 'c'])
    >>> pd.CategoricalIndex(c)
    CategoricalIndex(['a', 'b', 'c', 'a', 'b', 'c'], categories=['a', 'b', 'c'], ordered=False, dtype='category')  # noqa

    Ordered ``CategoricalIndex`` can have a min and max value.

    >>> ci = pd.CategoricalIndex(['a','b','c','a','b','c'], ordered=True,
    ...                          categories=['c', 'b', 'a'])
    >>> ci
    CategoricalIndex(['a', 'b', 'c', 'a', 'b', 'c'], categories=['c', 'b', 'a'], ordered=True, dtype='category')  # noqa
    >>> ci.min()
    'c'
    """

    _typ = "categoricalindex"

    _raw_inherit = {
        "argsort",
        "_internal_get_values",
        "tolist",
        "codes",
        "categories",
        "ordered",
        "_reverse_indexer",
        "searchsorted",
    }

    codes: np.ndarray
    categories: Index

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
        cls,
        data=None,
        categories=None,
        ordered=None,
        dtype=None,
        copy=False,
        name=None,
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

        data = cls._create_categorical(data, dtype=dtype)

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
    def _create_categorical(cls, data, dtype=None):
        """
        *this is an internal non-public method*

        create the correct categorical from data and the properties

        Parameters
        ----------
        data : data for new Categorical
        dtype : CategoricalDtype, defaults to existing

        Returns
        -------
        Categorical
        """
        if isinstance(data, (cls, ABCSeries)) and is_categorical_dtype(data):
            data = data.values

        if not isinstance(data, ABCCategorical):
            return Categorical(data, dtype=dtype)

        if isinstance(dtype, CategoricalDtype) and dtype != data.dtype:
            # we want to silently ignore dtype='category'
            data = data._set_dtype(dtype)
        return data

    @classmethod
    def _simple_new(cls, values, name=None, dtype=None, **kwargs):
        result = object.__new__(cls)

        values = cls._create_categorical(values, dtype=dtype)
        result._data = values
        result.name = name
        for k, v in kwargs.items():
            setattr(result, k, v)

        result._reset_identity()
        result._no_setting_name = False
        return result

    # --------------------------------------------------------------------

    @Appender(_index_shared_docs["_shallow_copy"])
    def _shallow_copy(self, values=None, dtype=None, **kwargs):
        if dtype is None:
            dtype = self.dtype
        return super()._shallow_copy(values=values, dtype=dtype, **kwargs)

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
            other = CategoricalIndex(self._create_categorical(other, dtype=self.dtype))
            if not other.isin(values).all():
                raise TypeError(
                    "cannot append a non-category item to a CategoricalIndex"
                )

        return other

    def equals(self, other):
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

    # --------------------------------------------------------------------

    @property
    def inferred_type(self) -> str:
        return "categorical"

    @property
    def values(self):
        """ return the underlying data, which is a Categorical """
        return self._data

    def _wrap_setop_result(self, other, result):
        name = get_op_result_name(self, other)
        return self._shallow_copy(result, name=name)

    @Appender(_index_shared_docs["contains"] % _index_doc_kwargs)
    def __contains__(self, key) -> bool:
        # if key is a NaN, check if any NaN is in self.
        if is_scalar(key) and isna(key):
            return self.hasnans

        return contains(self, key, container=self._engine)

    def __array__(self, dtype=None):
        """ the array interface, return my values """
        return np.array(self._data, dtype=dtype)

    @Appender(_index_shared_docs["astype"])
    def astype(self, dtype, copy=True):
        if is_interval_dtype(dtype):
            from pandas import IntervalIndex

            return IntervalIndex(np.array(self))
        elif is_categorical_dtype(dtype):
            # GH 18630
            dtype = self.dtype.update_dtype(dtype)
            if dtype == self.dtype:
                return self.copy() if copy else self

        return super().astype(dtype=dtype, copy=copy)

    @cache_readonly
    def _isnan(self):
        """ return if each value is nan"""
        return self._data.codes == -1

    @Appender(ibase._index_shared_docs["fillna"])
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

    # introspection
    @cache_readonly
    def is_unique(self) -> bool:
        return self._engine.is_unique

    @property
    def is_monotonic_increasing(self):
        return self._engine.is_monotonic_increasing

    @property
    def is_monotonic_decreasing(self) -> bool:
        return self._engine.is_monotonic_decreasing

    @Appender(_index_shared_docs["index_unique"] % _index_doc_kwargs)
    def unique(self, level=None):
        if level is not None:
            self._validate_index_level(level)
        result = self.values.unique()
        # CategoricalIndex._shallow_copy keeps original dtype
        # if not otherwise specified
        return self._shallow_copy(result, dtype=result.dtype)

    @Appender(Index.duplicated.__doc__)
    def duplicated(self, keep="first"):
        codes = self.codes.astype("i8")
        return duplicated_int64(codes, keep)

    def _to_safe_for_reshape(self):
        """ convert to object if we are a categorical """
        return self.astype("object")

    def get_loc(self, key, method=None):
        """
        Get integer location, slice or boolean mask for requested label.

        Parameters
        ----------
        key : label
        method : {None}
            * default: exact matches only.

        Returns
        -------
        loc : int if unique index, slice if monotonic index, else mask

        Raises
        ------
        KeyError : if the key is not in the index

        Examples
        --------
        >>> unique_index = pd.CategoricalIndex(list('abc'))
        >>> unique_index.get_loc('b')
        1

        >>> monotonic_index = pd.CategoricalIndex(list('abbc'))
        >>> monotonic_index.get_loc('b')
        slice(1, 3, None)

        >>> non_monotonic_index = pd.CategoricalIndex(list('abcb'))
        >>> non_monotonic_index.get_loc('b')
        array([False,  True, False,  True], dtype=bool)
        """
        code = self.categories.get_loc(key)
        code = self.codes.dtype.type(code)
        try:
            return self._engine.get_loc(code)
        except KeyError:
            raise KeyError(key)

    def get_value(self, series: AnyArrayLike, key: Any):
        """
        Fast lookup of value from 1-dimensional ndarray. Only use this if you
        know what you're doing

        Parameters
        ----------
        series : Series, ExtensionArray, Index, or ndarray
            1-dimensional array to take values from
        key: : scalar
            The value of this index at the position of the desired value,
            otherwise the positional index of the desired value

        Returns
        -------
        Any
            The element of the series at the position indicated by the key
        """
        try:
            k = com.values_from_object(key)
            k = self._convert_scalar_indexer(k, kind="getitem")
            indexer = self.get_loc(k)
            return series.take([indexer])[0]
        except (KeyError, TypeError):
            pass

        # we might be a positional inexer
        return super().get_value(series, key)

    @Appender(_index_shared_docs["where"])
    def where(self, cond, other=None):
        # TODO: Investigate an alternative implementation with
        # 1. copy the underlying Categorical
        # 2. setitem with `cond` and `other`
        # 3. Rebuild CategoricalIndex.
        if other is None:
            other = self._na_value
        values = np.where(cond, self.values, other)
        cat = Categorical(values, dtype=self.dtype)
        return self._shallow_copy(cat, **self._get_attributes_dict())

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
        """ reindex from a non-unique; which CategoricalIndex's are almost
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

        if isinstance(target, CategoricalIndex) and self.values.is_dtype_equal(target):
            if self.values.equals(target.values):
                # we have the same codes
                codes = target.codes
            else:
                codes = _recode_for_categories(
                    target.codes, target.categories, self.values.categories
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
            target = target.values

        codes = self.categories.get_indexer(target)
        indexer, missing = self._engine.get_indexer_non_unique(codes)
        return ensure_platform_int(indexer), missing

    @Appender(_index_shared_docs["_convert_scalar_indexer"])
    def _convert_scalar_indexer(self, key, kind=None):
        if kind == "loc":
            try:
                return self.categories._convert_scalar_indexer(key, kind=kind)
            except TypeError:
                self._invalid_indexer("label", key)
        return super()._convert_scalar_indexer(key, kind=kind)

    @Appender(_index_shared_docs["_convert_list_indexer"])
    def _convert_list_indexer(self, keyarr, kind=None):
        # Return our indexer or raise if all of the values are not included in
        # the categories

        if self.categories._defer_to_indexing:
            indexer = self.categories._convert_list_indexer(keyarr, kind=kind)
            return Index(self.codes).get_indexer_for(indexer)

        indexer = self.categories.get_indexer(np.asarray(keyarr))
        if (indexer == -1).any():
            raise KeyError(
                "a list-indexer must only include values that are in the categories"
            )

        return self.get_indexer(keyarr)

    @Appender(_index_shared_docs["_convert_arr_indexer"])
    def _convert_arr_indexer(self, keyarr):
        keyarr = com.asarray_tuplesafe(keyarr)

        if self.categories._defer_to_indexing:
            return keyarr

        return self._shallow_copy(keyarr)

    @Appender(_index_shared_docs["_convert_index_indexer"])
    def _convert_index_indexer(self, keyarr):
        return self._shallow_copy(keyarr)

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
        return self._create_from_codes(taken)

    take_nd = take

    @Appender(_index_shared_docs["_maybe_cast_slice_bound"])
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
        return self._shallow_copy_with_infer(self.values.map(mapper))

    def delete(self, loc):
        """
        Make new Index with passed location(-s) deleted

        Returns
        -------
        new_index : Index
        """
        return self._create_from_codes(np.delete(self.codes, loc))

    def insert(self, loc, item):
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
        return CategoricalIndex._concat_same_dtype(self, to_concat, name)

    def _concat_same_dtype(self, to_concat, name):
        """
        Concatenate to_concat which has the same class
        ValueError if other is not in the categories
        """
        codes = np.concatenate([self._is_dtype_compat(c).codes for c in to_concat])
        result = self._create_from_codes(codes, name=name)
        # if name is None, _create_from_codes sets self.name
        result.name = name
        return result

    @classmethod
    def _add_comparison_methods(cls):
        """ add in comparison methods """

        def _make_compare(op):
            opname = f"__{op.__name__}__"

            def _evaluate_compare(self, other):
                with np.errstate(all="ignore"):
                    result = op(self.array, other)
                if isinstance(result, ABCSeries):
                    # Dispatch to pd.Categorical returned NotImplemented
                    # and we got a Series back; down-cast to ndarray
                    result = result._values
                return result

            return compat.set_function_name(_evaluate_compare, opname, cls)

        cls.__eq__ = _make_compare(operator.eq)
        cls.__ne__ = _make_compare(operator.ne)
        cls.__lt__ = _make_compare(operator.lt)
        cls.__gt__ = _make_compare(operator.gt)
        cls.__le__ = _make_compare(operator.le)
        cls.__ge__ = _make_compare(operator.ge)

    def _delegate_property_get(self, name, *args, **kwargs):
        """ method delegation to the ._values """
        prop = getattr(self._values, name)
        return prop  # no wrapping for now

    def _delegate_method(self, name, *args, **kwargs):
        """ method delegation to the ._values """
        method = getattr(self._values, name)
        if "inplace" in kwargs:
            raise ValueError("cannot use inplace with CategoricalIndex")
        res = method(*args, **kwargs)
        if is_scalar(res) or name in self._raw_inherit:
            return res
        return CategoricalIndex(res, name=self.name)


CategoricalIndex._add_numeric_methods_add_sub_disabled()
CategoricalIndex._add_numeric_methods_disabled()
CategoricalIndex._add_logical_methods_disabled()
CategoricalIndex._add_comparison_methods()
