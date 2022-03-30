"""
Public API for DataFrame exchange protocol.
"""

from pandas.core.exchange.implementation import from_dataframe
from pandas.core.exchange.dataframe_protocol import DataFrame

__all__ = ["from_dataframe", "DataFrame"]
