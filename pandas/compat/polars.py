"""support polars compatibility across versions"""

from __future__ import annotations

from pandas.util.version import Version

try:
    import polars as pl

    _plv = Version(Version(pl.__version__).base_version)
    HAS_POLARS = _plv >= Version("0.20.0")  # Minimum version for to_pandas compatibility
except ImportError:
    HAS_POLARS = False
