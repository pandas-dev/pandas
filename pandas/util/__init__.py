from pandas.util._decorators import Appender, Substitution, cache_readonly  # noqa

from pandas import compat
from pandas.core.util.hashing import hash_array, hash_pandas_object  # noqa

# compatibility for import pandas; pandas.util.testing

if compat.PY37:

    def __getattr__(name):
        if name == "testing":
            import pandas.util.testing

            return pandas.util.testing
        else:
            raise AttributeError(f"module 'pandas.util' has no attribute '{name}'")


else:

    class _testing:
        def __getattr__(self, item):
            import pandas.util.testing

            return getattr(pandas.util.testing, item)

    testing = _testing()


del compat
