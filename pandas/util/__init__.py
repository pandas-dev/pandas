from pandas.util._decorators import (  # noqa:F401
    Appender,
    Substitution,
    cache_readonly,
)

from pandas._core.util.hashing import (  # noqa:F401
    hash_array,
    hash_pandas_object,
)


def __getattr__(name):
    if name == "testing":
        import pandas.util.testing

        return pandas.util.testing
    else:
        raise AttributeError(f"module 'pandas.util' has no attribute '{name}'")
