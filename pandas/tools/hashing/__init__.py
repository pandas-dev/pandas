import warnings


warnings.warn(
    "The pandas.tools.hashimg module is deprecated and will be removed in a "
    "future version. Please import from the pandas.util.hashing module instead.",
    FutureWarning, stacklevel=2
)

from pandas.util.hashing import *
