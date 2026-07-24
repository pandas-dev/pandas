"""
Various utilities
"""
from __future__ import annotations


class JupyterEventsVersionWarning(UserWarning):
    """Emitted when an event schema version is an `int` when it should be `str`."""
