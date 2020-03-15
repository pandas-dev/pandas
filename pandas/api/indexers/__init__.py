"""
Public API for Rolling Window Indexers.
"""

from pandas.core.indexers import check_array_indexer
from pandas.core.window.indexers import BaseIndexer

__all__ = ["check_array_indexer", "BaseIndexer"]
