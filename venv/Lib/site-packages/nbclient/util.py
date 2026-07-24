"""General utility methods"""

# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.
from __future__ import annotations

import inspect
from collections.abc import Callable
from typing import Any

from jupyter_core.utils import ensure_async, run_sync

__all__ = ["ensure_async", "run_sync", "run_hook"]


async def run_hook(hook: Callable[..., Any] | None, **kwargs: Any) -> None:
    """Run a hook callback."""
    if hook is None:
        return
    res = hook(**kwargs)
    if inspect.isawaitable(res):
        await res
