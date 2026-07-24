"""
This module provides the base implementation for the :module:`hyperscan`
backend.

WARNING: The *pathspec._backends.hyperscan* package is not part of the public
API. Its contents and structure are likely to change.
"""
from __future__ import annotations

from typing import (
	Optional)  # Replaced by `X | None` in 3.10.

try:
	import hyperscan
	hyperscan_error = None
except ModuleNotFoundError as e:
	hyperscan = None  # type: ignore[assignment]
	hyperscan_error = e.with_traceback(None)

hyperscan_error: Optional[ModuleNotFoundError]  # type: ignore[no-redef]
"""
*hyperscan_error* (:class:`ModuleNotFoundError` or :data:`None`) is the
hyperscan import error.
"""
