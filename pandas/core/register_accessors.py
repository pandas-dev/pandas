from pandas.core import strings
from pandas.core.accessor import (register_dataframe_accessor,
                                  register_index_accessor,
                                  register_series_accessor)
import pandas.core.categorical
from pandas.core.indexes.accessors import CombinedDatetimelikeProperties
from pandas.plotting._core import SeriesPlotMethods, FramePlotMethods


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

# ----
# plot
# ----

# TODO: see if this triggers the actual mpl import...
register_series_accessor("plot")(SeriesPlotMethods)
register_dataframe_accessor("plot")(FramePlotMethods)

__all__ = []
