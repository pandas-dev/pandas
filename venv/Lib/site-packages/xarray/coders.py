"""
This module provides coder objects that encapsulate the
"encoding/decoding" process.
"""

from xarray.coding.times import CFDatetimeCoder, CFTimedeltaCoder

__all__ = ["CFDatetimeCoder", "CFTimedeltaCoder"]
