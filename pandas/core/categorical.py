import warnings

# TODO: Remove after 0.23.x
warnings.warn("'pandas.core' is private. Use 'pandas.Categorical'",
              FutureWarning, stacklevel=2)

from pandas.core.arrays import Categorical  # noqa
