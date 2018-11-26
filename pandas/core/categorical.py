import warnings

from pandas.core.dtypes.dtypes import CategoricalDtype  # noqa

from pandas.core.arrays import Categorical  # noqa

# TODO: Remove after 0.23.x
warnings.warn("'pandas.core' is private. Use 'pandas.Categorical'",
              FutureWarning, stacklevel=2)
