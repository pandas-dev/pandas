"""
Public API for DataFrame interchange protocol.
"""

from pandas.core.interchange._from_dataframe import from_dataframe
from pandas.core.interchange.dataframe_protocol import DataFrame

__all__ = ["from_dataframe", "DataFrame"]
