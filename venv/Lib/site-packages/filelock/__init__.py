"""
A platform independent file lock that supports the with-statement.

.. autodata:: filelock.__version__
    :no-value:

"""

from __future__ import annotations

import sys
import warnings
from typing import TYPE_CHECKING, Final

from ._api import AcquireReturnProxy, BaseFileLock, CloseErrorPolicy, ContextErrorPolicy, LockOptions
from ._descriptor import lock_descriptor, unlock_descriptor
from ._error import LeaseSettingsMismatch, SoftFileLockLifetimeWarning, SoftFileLockProtocolError, Timeout
from ._lease import LeaseCompromise, SoftFileLease
from ._marker import MarkerSoftFileLock, OwnerRecord

if TYPE_CHECKING:
    from ._async_read_write import (
        AsyncAcquireReadWriteReturnProxy,
        AsyncReadWriteLock,
    )
    from ._read_write import ReadWriteLock
else:
    try:
        from ._async_read_write import AsyncAcquireReadWriteReturnProxy, AsyncReadWriteLock
        from ._read_write import ReadWriteLock
    except ImportError:  # pragma: lacks sqlite3
        AsyncAcquireReadWriteReturnProxy = None
        AsyncReadWriteLock = None
        ReadWriteLock = None

from ._soft import SoftFileLock
from ._soft_rw import AsyncAcquireSoftReadWriteReturnProxy, AsyncSoftReadWriteLock, SoftReadWriteLock
from ._strict import StrictSoftFileClaim, StrictSoftFileClaimState, StrictSoftFileLock
from ._unix import UnixFileLock, has_fcntl
from ._windows import WindowsFileLock
from .asyncio import (
    AsyncAcquireReturnProxy,
    AsyncSoftFileLease,
    AsyncSoftFileLock,
    AsyncStrictSoftFileLock,
    AsyncUnixFileLock,
    AsyncWindowsFileLock,
    BaseAsyncFileLock,
)
from .version import version

#: version of the project as a string
__version__: Final[str] = version


if sys.platform == "win32":  # pragma: win32 cover
    _FileLock: type[BaseFileLock] = WindowsFileLock
    _AsyncFileLock: type[BaseAsyncFileLock] = AsyncWindowsFileLock
else:  # pragma: win32 no cover # ruff:ignore[collapsible-else-if]  # the else carries the win32 no-cover pragma
    if has_fcntl:
        _FileLock: type[BaseFileLock] = UnixFileLock
        _AsyncFileLock: type[BaseAsyncFileLock] = AsyncUnixFileLock
    else:
        _FileLock = SoftFileLock
        _AsyncFileLock = AsyncSoftFileLock
        warnings.warn("only soft file lock is available", stacklevel=2)

if TYPE_CHECKING:
    FileLock = SoftFileLock
    AsyncFileLock = AsyncSoftFileLock
else:
    #: Alias for the lock, which should be used for the current platform.
    FileLock = _FileLock
    AsyncFileLock = _AsyncFileLock


__all__ = [
    "AcquireReturnProxy",
    "AsyncAcquireReadWriteReturnProxy",
    "AsyncAcquireReturnProxy",
    "AsyncAcquireSoftReadWriteReturnProxy",
    "AsyncFileLock",
    "AsyncReadWriteLock",
    "AsyncSoftFileLease",
    "AsyncSoftFileLock",
    "AsyncSoftReadWriteLock",
    "AsyncStrictSoftFileLock",
    "AsyncUnixFileLock",
    "AsyncWindowsFileLock",
    "BaseAsyncFileLock",
    "BaseFileLock",
    "CloseErrorPolicy",
    "ContextErrorPolicy",
    "FileLock",
    "LeaseCompromise",
    "LeaseSettingsMismatch",
    "LockOptions",
    "MarkerSoftFileLock",
    "OwnerRecord",
    "ReadWriteLock",
    "SoftFileLease",
    "SoftFileLock",
    "SoftFileLockLifetimeWarning",
    "SoftFileLockProtocolError",
    "SoftReadWriteLock",
    "StrictSoftFileClaim",
    "StrictSoftFileClaimState",
    "StrictSoftFileLock",
    "Timeout",
    "UnixFileLock",
    "WindowsFileLock",
    "__version__",
    "lock_descriptor",
    "unlock_descriptor",
]
