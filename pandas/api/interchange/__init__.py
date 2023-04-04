# pylint: disable=undefined-all-variable
"""
Public API for DataFrame interchange protocol.
"""


def __getattr__(key: str):
    # lazy imports to speed 'import pandas as pd'
    if key == "DataFrame":
        from pandas.core.interchange.dataframe_protocol import DataFrame

        return DataFrame
    elif key == "from_dataframe":
        from pandas.core.interchange.from_dataframe import from_dataframe

        return from_dataframe
    raise AttributeError(key)


__all__ = [  # pyright: ignore[reportUnsupportedDunderAll]
    "from_dataframe",
    "DataFrame",
]
