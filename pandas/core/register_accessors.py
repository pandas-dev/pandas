from pandas.core import strings
from pandas.core.accessor import (register_index_accessor,
                                  register_series_accessor)
import pandas.core.categorical


@register_series_accessor("cat")
class CategoricalAccessor(pandas.core.categorical.CategoricalAccessor):
    pass


# ---
# str
# ---

@register_index_accessor("str")
@register_series_accessor("str")
class StringAccessor(strings.StringMethods):
    pass


__all__ = []
