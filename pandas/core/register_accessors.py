from pandas.core import strings
from pandas.core.accessor import (register_index_accessor,
                                  register_series_accessor)
import pandas.core.categorical
from pandas.core.indexes.accessors import CombinedDatetimelikeProperties


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


# --
# dt
# --

register_series_accessor("dt", cache=False)(CombinedDatetimelikeProperties)


__all__ = []
