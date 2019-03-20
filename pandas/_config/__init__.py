"""
pandas._config is considered explicitly upstream of everything else in pandas,
should have no intra-pandas dependencies.
"""
__all__ = ["config"]
from pandas._config import config
