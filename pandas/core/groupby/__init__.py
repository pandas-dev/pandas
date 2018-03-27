import warnings

# TODO: Remove after 0.XX.x
warnings.warn("'pandas.core.groupby' has moved. Use "
              "'pandas.core.groupby.groupby' instead", FutureWarning,
              stacklevel=2)

from pandas.core.groupby.groupby import Grouper
from pandas.core.groupby.groupby import groupby
from pandas.core.groupby.groupby import BinGrouper
from pandas.core.groupby.groupby import _GroupBy
from pandas.core.groupby.groupby import GroupBy
from pandas.core.groupby.groupby import SeriesGroupBy
from pandas.core.groupby.groupby import _pipe_template
from pandas.core.groupby.groupby import PanelGroupBy
from pandas.core.groupby.groupby import Grouping
from pandas.core.groupby.groupby import SpecificationError
from pandas.core.groupby.groupby import DataError
from pandas.core.groupby.groupby import generate_bins_generic
