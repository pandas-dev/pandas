"""Cross-process and cross-host reader/writer lock on :class:`~filelock.SoftFileLock` primitives."""

from __future__ import annotations

from ._async import AsyncAcquireSoftReadWriteReturnProxy, AsyncSoftReadWriteLock
from ._sync import SoftReadWriteLock

__all__ = [
    "AsyncAcquireSoftReadWriteReturnProxy",
    "AsyncSoftReadWriteLock",
    "SoftReadWriteLock",
]
