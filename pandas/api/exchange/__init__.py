"""
Public API for DataFrame exchange protocol.
"""

from pandas.core.exchange.dataframe_protocol import DataFrame
from pandas.core.exchange.from_dataframe import from_dataframe

__all__ = ["from_dataframe", "DataFrame"]
