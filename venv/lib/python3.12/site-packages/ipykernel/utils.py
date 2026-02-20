"""Utilities"""

from __future__ import annotations

import asyncio
import sys
import typing as t
from collections.abc import Mapping
from contextvars import copy_context
from functools import partial, wraps

if t.TYPE_CHECKING:
    from collections.abc import Callable
    from contextvars import Context


class LazyDict(Mapping[str, t.Any]):
    """Lazy evaluated read-only dictionary.

    Initialised with a dictionary of key-value pairs where the values are either
    constants or callables. Callables are evaluated each time the respective item is
    read.
    """

    def __init__(self, dict):
        self._dict = dict

    def __getitem__(self, key):
        item = self._dict.get(key)
        return item() if callable(item) else item

    def __len__(self):
        return len(self._dict)

    def __iter__(self):
        return iter(self._dict)


T = t.TypeVar("T")
U = t.TypeVar("U")
V = t.TypeVar("V")


def _async_in_context(
    f: Callable[..., t.Coroutine[T, U, V]], context: Context | None = None
) -> Callable[..., t.Coroutine[T, U, V]]:
    """
    Wrapper to run a coroutine in a persistent ContextVar Context.

    Backports asyncio.create_task(context=...) behavior from Python 3.11
    """
    if context is None:
        context = copy_context()

    if sys.version_info >= (3, 11):

        @wraps(f)
        async def run_in_context(*args, **kwargs):
            coro = f(*args, **kwargs)
            return await asyncio.create_task(coro, context=context)

        return run_in_context

    # don't need this backport when we require 3.11
    # context_holder so we have a modifiable container for later calls
    context_holder = [context]  # type: ignore[unreachable]

    async def preserve_context(f, *args, **kwargs):
        """call a coroutine, preserving the context after it is called"""
        try:
            return await f(*args, **kwargs)
        finally:
            # persist changes to the context for future calls
            context_holder[0] = copy_context()

    @wraps(f)
    async def run_in_context_pre311(*args, **kwargs):
        ctx = context_holder[0]
        return await ctx.run(partial(asyncio.create_task, preserve_context(f, *args, **kwargs)))

    return run_in_context_pre311
