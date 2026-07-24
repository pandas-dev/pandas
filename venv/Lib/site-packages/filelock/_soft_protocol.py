from __future__ import annotations

from typing import Final

STRICT_SOFT_SENTINEL_RECORD: Final[str] = "1\nfilelock-strict-v1\x00\n0\n"

__all__ = [
    "STRICT_SOFT_SENTINEL_RECORD",
]
