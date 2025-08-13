from __future__ import annotations

try:
    from sqlalchemy.exc import RemovedIn20Warning
except ImportError:

    class _RemovedIn20Warning(Warning):
        pass

    RemovedIn20Warning = _RemovedIn20Warning
