import numpy as np
import pandas.lib as lib
import pandas.index as _index

from pandas import compat
from pandas.compat.numpy import function as nv
from pandas.util.decorators import (Appender, cache_readonly,
                                    deprecate_kwarg)
from pandas.core.config import get_option
from pandas.indexes.base import Index, _index_shared_docs
import pandas.core.base as base
import pandas.core.common as com
import pandas.core.missing as missing
import pandas.indexes.base as ibase


class CategoricalIndex(Index, base.PandasDelegate):
    """

    Immutable Index implementing an ordered, sliceable set. CategoricalIndex
    represents a sparsely populated Index with an underlying Categorical.

    .. versionadded:: 0.16.1

    Parameters
    ----------
    data : array-like or Categorical, (1-dimensional)
    categories : optional, array-like
        categories for the CategoricalIndex
    ordered : boolean,
        designating if the categories are ordered
    copy : bool
        Make a copy of input ndarray
    name : object
        Name to be stored in the index

    """

    _typ = 'categoricalindex'
    _engine_type = _index.Int64Engine
    _attributes = ['name']

    def __new__(cls, data=None, categories=None, ordered=None, dtype=None,
                copy=False, name=None, fastpath=False, **kwargs):

        if fastpath:
            return cls._simple_new(data, name=name)

        if isinstance(data, com.ABCCategorical):
            data = cls._create_categorical(cls, data, categories, ordered)
        elif isinstance(data, CategoricalIndex):
            data = data._data
            data = cls._create_categorical(cls, data, categories, ordered)
        else:

            # don't allow scalars
            # if data is None, then categories must be provided
            if lib.isscalar(data):
                if data is not None or categories is None:
                    cls._scalar_data_error(data)
                data = []
            data = cls._create_categorical(cls, data, categories, ordered)

        if copy:
            data = data.copy()

        return cls._simple_new(data, name=name)

    def _create_from_codes(self, codes, categories=None, ordered=None,
                           name=None):
        """
        *this is an internal non-public method*

        create the correct categorical from codes

        Parameters
        ----------
        codes : new codes
        categories : optional categories, defaults to existing
        ordered : optional ordered attribute, defaults to existing
        name : optional name attribute, defaults to existing

        Returns
        -------
        CategoricalIndex
        """

        from pandas.core.categorical import Categorical
        if categories is None:
            categories = self.categories
        if ordered is None:
            ordered = self.ordered
        if name is None:
            name = self.name
        cat = Categorical.from_codes(codes, categories=categories,
                                     ordered=self.ordered)
        return CategoricalIndex(cat, name=name)

    @staticmethod
    def _create_categorical(self, data, categories=None, ordered=None):
        """
        *this is an internal non-public method*

        create the correct categorical from data and the properties

        Parameters
        ----------
        data : data for new Categorical
        categories : optional categories, defaults to existing
        ordered : optional ordered attribute, defaults to existing

        Returns
        -------
        Categorical
        """
        if not isinstance(data, com.ABCCategorical):
            from pandas.core.categorical import Categorical
            data = Categorical(data, categories=categories, ordered=ordered)
        else:
            if categories is not None:
                data = data.set_categories(categories)
            if ordered is not None:
                data = data.set_ordered(ordered)
        return data

    @classmethod
    def _simple_new(cls, values, name=None, categories=None, ordered=None,
                    **kwargs):
        result = object.__new__(cls)

        values = cls._create_categorical(cls, values, categories, ordered)
        result._data = values
        result.name = name
        for k, v in compat.iteritems(kwargs):
            setattr(result, k, v)

        result._reset_identity()
        return result

    @Appender(_index_shared_docs['_shallow_copy'])
    def _shallow_copy(self, values=None, categories=None, ordered=None,
                      **kwargs):
        # categories and ordered can't be part of attributes,
        # as these are properties
        if categories is None:
            categories = self.categories
        if ordered is None:
            ordered = self.ordered
        return super(CategoricalIndex,
                     self)._shallow_copy(values=values, categories=categories,
                                         ordered=ordered, **kwargs)

    def _is_dtype_compat(self, other):
        """
        *this is an internal non-public method*

        provide a comparison between the dtype of self and other (coercing if
        needed)

        Raises
        ------
        TypeError if the dtypes are not compatible
        """
        if com.is_categorical_dtype(other):
            if isinstance(other, CategoricalIndex):
                other = other._values
            if not other.is_dtype_equal(self):
                raise TypeError("categories must match existing categories "
                                "when appending")
        else:
            values = other
            if not com.is_list_like(values):
                values = [values]
            other = CategoricalIndex(self._create_categorical(
                self, other, categories=self.categories, ordered=self.ordered))
            if not other.isin(values).all():
                raise TypeError("cannot append a non-category item to a "
                                "CategoricalIndex")

        return other

    def equals(self, other):
        """
        Determines if two CategorialIndex objects contain the same elements.
        """
        if self.is_(other):
            return True

        try:
            other = self._is_dtype_compat(other)
            return com.array_equivalent(self._data, other)
        except (TypeError, ValueError):
            pass

        return False

    @property
    def _formatter_func(self):
        return self.categories._formatter_func

    def _format_attrs(self):
        """
        Return a list of tuples of the (attr,formatted_value)
        """
        max_categories = (10 if get_option("display.max_categories") == 0 else
                          get_option("display.max_categories"))
        attrs = [
            ('categories',
             ibase.default_pprint(self.categories,
                                  max_seq_items=max_categories)),
            ('ordered', self.ordered)]
        if self.name is not None:
            attrs.append(('name', ibase.default_pprint(self.name)))
        attrs.append(('dtype', "'%s'" % self.dtype))
        max_seq_items = get_option('display.max_seq_items') or len(self)
        if len(self) > max_seq_items:
            attrs.append(('length', len(self)))
        return attrs

    @property
    def inferred_type(self):
        return 'categorical'

    @property
    def values(self):
        """ return the underlying data, which is a Categorical """
        return self._data

    def get_values(self):
        """ return the underlying data as an ndarray """
        return self._data.get_values()

    @property
    def codes(self):
        return self._data.codes

    @property
    def categories(self):
        return self._data.categories

    @property
    def ordered(self):
        return self._data.ordered

    def __contains__(self, key):
        hash(key)
        return key in self.values

    def __array__(self, dtype=None):
        """ the array interface, return my values """
        return np.array(self._data, dtype=dtype)

    @cache_readonly
    def _isnan(self):
        """ return if each value is nan"""
        return self._data.codes == -1

    @Appender(ibase._index_shared_docs['fillna'])
    def fillna(self, value, downcast=None):
        self._assert_can_do_op(value)
        return CategoricalIndex(self._data.fillna(value), name=self.name)

    def argsort(self, *args, **kwargs):
        return self.values.argsort(*args, **kwargs)

    @cache_readonly
    def _engine(self):

        # we are going to look things up with the codes themselves
        return self._engine_type(lambda: self.codes.astype('i8'), len(self))

    @cache_readonly
    def is_unique(self):
        return not self.duplicated().any()

    @deprecate_kwarg('take_last', 'keep', mapping={True: 'last',
                                                   False: 'first'})
    @Appender(base._shared_docs['duplicated'] % ibase._index_doc_kwargs)
    def duplicated(self, keep='first'):
        from pandas.hashtable import duplicated_int64
        return duplicated_int64(self.codes.astype('i8'), keep)

    def _to_safe_for_reshape(self):
        """ convert to object if we are a categorical """
        return self.astype('object')

    def get_loc(self, key, method=None):
        """
        Get integer location for requested label

        Parameters
        ----------
        key : label
        method : {None}
            * default: exact matches only.

        Returns
        -------
        loc : int if unique index, possibly slice or mask if not
        """
        codes = self.categories.get_loc(key)
        if (codes == -1):
            raise KeyError(key)
        return self._engine.get_loc(codes)

    def _can_reindex(self, indexer):
        """ always allow reindexing """
        pass

    def reindex(self, target, method=None, level=None, limit=None,
                tolerance=None):
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
            raise NotImplementedError("argument method is not implemented for "
                                      "CategoricalIndex.reindex")
        if level is not None:
            raise NotImplementedError("argument level is not implemented for "
                                      "CategoricalIndex.reindex")
        if limit is not None:
            raise NotImplementedError("argument limit is not implemented for "
                                      "CategoricalIndex.reindex")

        target = ibase._ensure_index(target)

        if not com.is_categorical_dtype(target) and not target.is_unique:
            raise ValueError("cannot reindex with a non-unique indexer")

        indexer, missing = self.get_indexer_non_unique(np.array(target))
        new_target = self.take(indexer)

        # filling in missing if needed
        if len(missing):
            cats = self.categories.get_indexer(target)

            if (cats == -1).any():
                # coerce to a regular index here!
                result = Index(np.array(self), name=self.name)
                new_target, indexer, _ = result._reindex_non_unique(
                    np.array(target))

            else:

                codes = new_target.codes.copy()
                codes[indexer == -1] = cats[missing]
                new_target = self._create_from_codes(codes)

        # we always want to return an Index type here
        # to be consistent with .reindex for other index types (e.g. they don't
        # coerce based on the actual values, only on the dtype)
        # unless we had an inital Categorical to begin with
        # in which case we are going to conform to the passed Categorical
        new_target = np.asarray(new_target)
        if com.is_categorical_dtype(target):
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

    def get_indexer(self, target, method=None, limit=None, tolerance=None):
        """
        Compute indexer and mask for new index given the current index. The
        indexer should be then used as an input to ndarray.take to align the
        current data to the new index. The mask determines whether labels are
        found or not in the current index

        Parameters
        ----------
        target : MultiIndex or Index (of tuples)
        method : {'pad', 'ffill', 'backfill', 'bfill'}
            pad / ffill: propagate LAST valid observation forward to next valid
            backfill / bfill: use NEXT valid observation to fill gap

        Notes
        -----
        This is a low-level method and probably should be used at your own risk

        Examples
        --------
        >>> indexer, mask = index.get_indexer(new_index)
        >>> new_values = cur_values.take(indexer)
        >>> new_values[-mask] = np.nan

        Returns
        -------
        (indexer, mask) : (ndarray, ndarray)
        """
        method = missing.clean_reindex_fill_method(method)
        target = ibase._ensure_index(target)

        if isinstance(target, CategoricalIndex):
            target = target.categories

        if method == 'pad' or method == 'backfill':
            raise NotImplementedError("method='pad' and method='backfill' not "
                                      "implemented yet for CategoricalIndex")
        elif method == 'nearest':
            raise NotImplementedError("method='nearest' not implemented yet "
                                      'for CategoricalIndex')
        else:

            codes = self.categories.get_indexer(target)
            indexer, _ = self._engine.get_indexer_non_unique(codes)

        return com._ensure_platform_int(indexer)

    def get_indexer_non_unique(self, target):
        """ this is the same for a CategoricalIndex for get_indexer; the API
        returns the missing values as well
        """
        target = ibase._ensure_index(target)

        if isinstance(target, CategoricalIndex):
            target = target.categories

        codes = self.categories.get_indexer(target)
        return self._engine.get_indexer_non_unique(codes)

    def _convert_list_indexer(self, keyarr, kind=None):
        """
        we are passed a list indexer.
        Return our indexer or raise if all of the values are not included in
        the categories
        """
        codes = self.categories.get_indexer(keyarr)
        if (codes == -1).any():
            raise KeyError("a list-indexer must only include values that are "
                           "in the categories")

        return None

    @Appender(_index_shared_docs['take'])
    def take(self, indices, axis=0, allow_fill=True,
             fill_value=None, **kwargs):
        nv.validate_take(tuple(), kwargs)
        indices = com._ensure_platform_int(indices)
        taken = self._assert_take_fillable(self.codes, indices,
                                           allow_fill=allow_fill,
                                           fill_value=fill_value,
                                           na_value=-1)
        return self._create_from_codes(taken)

    def map(self, mapper):
        """
        Apply mapper function to its categories (not codes).

        Parameters
        ----------
        mapper : callable
            Function to be applied. When all categories are mapped
            to different categories, the result will be Categorical which has
            the same order property as the original. Otherwise, the result will
            be np.ndarray.

        Returns
        -------
        applied : Categorical or np.ndarray.
        """
        return self.values.map(mapper)

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
        if (code == -1):
            raise TypeError("cannot insert an item into a CategoricalIndex "
                            "that is not already an existing category")

        codes = self.codes
        codes = np.concatenate((codes[:loc], code, codes[loc:]))
        return self._create_from_codes(codes)

    def append(self, other):
        """
        Append a collection of CategoricalIndex options together

        Parameters
        ----------
        other : Index or list/tuple of indices

        Returns
        -------
        appended : Index

        Raises
        ------
        ValueError if other is not in the categories
        """
        to_concat, name = self._ensure_compat_append(other)
        to_concat = [self._is_dtype_compat(c) for c in to_concat]
        codes = np.concatenate([c.codes for c in to_concat])
        return self._create_from_codes(codes, name=name)

    @classmethod
    def _add_comparison_methods(cls):
        """ add in comparison methods """

        def _make_compare(op):
            def _evaluate_compare(self, other):

                # if we have a Categorical type, then must have the same
                # categories
                if isinstance(other, CategoricalIndex):
                    other = other._values
                elif isinstance(other, Index):
                    other = self._create_categorical(
                        self, other._values, categories=self.categories,
                        ordered=self.ordered)

                if isinstance(other, (com.ABCCategorical, np.ndarray,
                                      com.ABCSeries)):
                    if len(self.values) != len(other):
                        raise ValueError("Lengths must match to compare")

                if isinstance(other, com.ABCCategorical):
                    if not self.values.is_dtype_equal(other):
                        raise TypeError("categorical index comparisions must "
                                        "have the same categories and ordered "
                                        "attributes")

                return getattr(self.values, op)(other)

            return _evaluate_compare

        cls.__eq__ = _make_compare('__eq__')
        cls.__ne__ = _make_compare('__ne__')
        cls.__lt__ = _make_compare('__lt__')
        cls.__gt__ = _make_compare('__gt__')
        cls.__le__ = _make_compare('__le__')
        cls.__ge__ = _make_compare('__ge__')

    def _delegate_method(self, name, *args, **kwargs):
        """ method delegation to the ._values """
        method = getattr(self._values, name)
        if 'inplace' in kwargs:
            raise ValueError("cannot use inplace with CategoricalIndex")
        res = method(*args, **kwargs)
        if lib.isscalar(res):
            return res
        return CategoricalIndex(res, name=self.name)

    @classmethod
    def _add_accessors(cls):
        """ add in Categorical accessor methods """

        from pandas.core.categorical import Categorical
        CategoricalIndex._add_delegate_accessors(
            delegate=Categorical, accessors=["rename_categories",
                                             "reorder_categories",
                                             "add_categories",
                                             "remove_categories",
                                             "remove_unused_categories",
                                             "set_categories",
                                             "as_ordered", "as_unordered",
                                             "min", "max"],
            typ='method', overwrite=True)


CategoricalIndex._add_numericlike_set_methods_disabled()
CategoricalIndex._add_numeric_methods_disabled()
CategoricalIndex._add_logical_methods_disabled()
CategoricalIndex._add_comparison_methods()
CategoricalIndex._add_accessors()
