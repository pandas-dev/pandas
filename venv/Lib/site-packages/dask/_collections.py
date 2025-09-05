from __future__ import annotations

from dask._dispatch import get_collection_type


def new_collection(expr):
    """Create new collection from an expr"""
    meta = expr._meta
    expr._name  # Ensure backend is imported
    return get_collection_type(meta)(expr)
