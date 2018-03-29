import warnings

# TODO: Remove after 0.XX.x
warnings.warn("'pandas.core.groupby' has moved. Use "
              "'pandas.core.groupby.groupby' instead", FutureWarning,
              stacklevel=2)

from pandas.core.groupby.groupby import (Grouper, groups,  # noqa: F401
                                         indices, get_group)
