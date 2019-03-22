"""
Timeseries API
"""
__all__ = ["infer_freq", "offset"]

from pandas.tseries.frequencies import infer_freq
import pandas.tseries.offsets as offsets
