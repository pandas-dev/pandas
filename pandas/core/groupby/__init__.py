import warnings

# TODO: Remove after 0.XX.x
warnings.warn("'pandas.core.groupby' has moved. Use "
              "'pandas.core.groupby.groupby' instead", FutureWarning,
              stacklevel=2)

from pandas.core.groupby.groupby import Grouper  # noqa: F401
from pandas.core.groupby.groupby import groupby  # noqa: F401
from pandas.core.groupby.groupby import BinGrouper  # noqa: F401
from pandas.core.groupby.groupby import _GroupBy  # noqa: F401
from pandas.core.groupby.groupby import GroupBy  # noqa: F401
from pandas.core.groupby.groupby import SeriesGroupBy  # noqa: F401
from pandas.core.groupby.groupby import _pipe_template  # noqa: F401
from pandas.core.groupby.groupby import PanelGroupBy  # noqa: F401
from pandas.core.groupby.groupby import Grouping  # noqa: F401
from pandas.core.groupby.groupby import SpecificationError  # noqa: F401
from pandas.core.groupby.groupby import DataError  # noqa: F401
from pandas.core.groupby.groupby import generate_bins_generic  # noqa: F401
