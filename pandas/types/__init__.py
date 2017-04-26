import warnings


warnings.warn(
    "The pandas.types module is deprecated and will be removed in a future "
    "version. Please import from the pandas.core.dtypes module instead.",
    FutureWarning, stacklevel=2
)
