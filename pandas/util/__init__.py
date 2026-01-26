def __getattr__(key: str):
    # These imports need to be lazy to avoid circular import errors
    if key == "hash_array":
        from pandas.core.util.hashing import hash_array

        return hash_array
    if key == "hash_pandas_object":
        from pandas.core.util.hashing import hash_pandas_object

        return hash_pandas_object
    if key == "Appender":
        import warnings

        from pandas.errors import Pandas4Warning
        from pandas.util._decorators import Appender

        warnings.warn(
            "Appender is deprecated and will be removed in a future version.",
            Pandas4Warning,
            stacklevel=2,
        )
        return Appender
    if key == "Substitution":
        import warnings

        from pandas.errors import Pandas4Warning
        from pandas.util._decorators import Substitution

        warnings.warn(
            "Substitution is deprecated and will be removed in a future version.",
            Pandas4Warning,
            stacklevel=2,
        )
        return Substitution
    if key == "cache_readonly":
        from pandas.util._decorators import cache_readonly

        return cache_readonly

    raise AttributeError(f"module 'pandas.util' has no attribute '{key}'")


def __dir__() -> list[str]:
    return [*list(globals().keys()), "hash_array", "hash_pandas_object"]
