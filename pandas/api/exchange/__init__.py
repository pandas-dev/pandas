"""
Public API for DataFrame exchange protocol.
"""

from pandas._core.exchange.dataframe_protocol import DataFrame
from pandas._core.exchange.from_dataframe import from_dataframe

__all__ = ["from_dataframe", "DataFrame"]
