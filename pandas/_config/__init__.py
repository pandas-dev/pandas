"""
pandas._config is considered explicitly upstream of everything else in pandas,
should have no intra-pandas dependencies.
"""
__all__ = ["config", "get_option"]
from pandas._config import config
from pandas._config.config import get_option