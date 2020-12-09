from pandas.util._decorators import Appender, Substitution, cache_readonly  # noqa

from pandas.core.util.hashing import hash_array, hash_pandas_object  # noqa


def __getattr__(name):
    if name == "testing":
        import pandas.util.testing

        return pandas.util.testing
    else:
        raise AttributeError(f"module 'pandas.util' has no attribute '{name}'")
