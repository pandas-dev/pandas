"""Version info."""
from __future__ import annotations

import re

__version__ = "0.10.4"

# Build up version_info tuple for backwards compatibility
pattern = r"(?P<major>\d+).(?P<minor>\d+).(?P<patch>\d+)(?P<rest>.*)"
match = re.match(pattern, __version__)
if match:
    parts: list[int | str] = [int(match[part]) for part in ["major", "minor", "patch"]]
    if match["rest"]:
        parts.append(match["rest"])
else:
    parts = []
version_info = tuple(parts)
