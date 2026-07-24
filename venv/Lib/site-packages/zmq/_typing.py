from __future__ import annotations

import sys

if sys.version_info >= (3, 10):
    from typing import TypeAlias
else:
    try:
        from typing_extensions import TypeAlias
    except ImportError:
        TypeAlias = type  # type: ignore
