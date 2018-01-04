from pandas.api.types import is_categorical_dtype
from pandas.core import strings
from pandas.core.accessor import (PandasDelegate, register_index_accessor,
                                  register_series_accessor)
from pandas.core.base import NoNewAttributesMixin, PandasObject
from pandas.core.categorical import Categorical
from pandas.core.indexes.accessors import CombinedDatetimelikeProperties


@register_series_accessor("cat")
class CategoricalAccessor(PandasDelegate, PandasObject, NoNewAttributesMixin):
    """
    Accessor object for categorical properties of the Series values.

    Be aware that assigning to `categories` is a inplace operation, while all
    methods return new categorical data per default (but can be called with
    `inplace=True`).

    Examples
    --------
    >>> s.cat.categories
    >>> s.cat.categories = list('abc')
    >>> s.cat.rename_categories(list('cab'))
    >>> s.cat.reorder_categories(list('cab'))
    >>> s.cat.add_categories(['d','e'])
    >>> s.cat.remove_categories(['d'])
    >>> s.cat.remove_unused_categories()
    >>> s.cat.set_categories(list('abcde'))
    >>> s.cat.as_ordered()
    >>> s.cat.as_unordered()

    """

    def __init__(self, data):
        self._validate_dtype(data)
        self.categorical = data.values
        self.index = data.index
        self.name = data.name
        self._freeze()

    @staticmethod
    def _validate_dtype(data):
        if not is_categorical_dtype(data):
            msg = "Can only use '.cat' accessor with 'category' dtype"
            raise AttributeError(msg)

    def _delegate_property_get(self, name):
        return getattr(self.categorical, name)

    def _delegate_property_set(self, name, new_values):
        return setattr(self.categorical, name, new_values)

    @property
    def codes(self):
        from pandas import Series
        return Series(self.categorical.codes, index=self.index)

    def _delegate_method(self, name, *args, **kwargs):
        from pandas import Series
        method = getattr(self.categorical, name)
        res = method(*args, **kwargs)
        if res is not None:
            return Series(res, index=self.index, name=self.name)

    @classmethod
    def _make_accessor(cls, data):
        cls._validate_dtype(data)
        return CategoricalAccessor(data.values, data.index,
                                   getattr(data, 'name', None),)


CategoricalAccessor._add_delegate_accessors(delegate=Categorical,
                                            accessors=["categories",
                                                       "ordered"],
                                            typ='property')
CategoricalAccessor._add_delegate_accessors(delegate=Categorical, accessors=[
    "rename_categories", "reorder_categories", "add_categories",
    "remove_categories", "remove_unused_categories", "set_categories",
    "as_ordered", "as_unordered"], typ='method')


# ---
# str
# ---

@register_index_accessor("str")
@register_series_accessor("str")
class StringAccessor(strings.StringMethods):
    pass


# --
# dt
# --

@register_index_accessor("dt")
@register_series_accessor("dt")
class DatetimeAccessor(CombinedDatetimelikeProperties):
    pass


__all__ = []
