"""
The Panel class is completely remove from pandas. GH #27101.

This stub remains for compatibility of imports on several
downstream packages, including: xarray & statsmodels

This compatibility should be removed in 1.0
"""


class Panel:

    def __init__(self, *args, **kwargs):
        raise ImportError("The Panel class is removed from pandas")
