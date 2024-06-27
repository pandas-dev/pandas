# util/concurrency.py
# Copyright (C) 2005-2024 the SQLAlchemy authors and contributors
# <see AUTHORS file>
#
# This module is part of SQLAlchemy and is released under
# the MIT License: https://www.opensource.org/licenses/mit-license.php
# mypy: allow-untyped-defs, allow-untyped-calls

from __future__ import annotations

import asyncio  # noqa
import typing
from typing import Any
from typing import Callable
from typing import Coroutine
from typing import TypeVar

have_greenlet = False
greenlet_error = None
try:
    import greenlet  # type: ignore[import-untyped,unused-ignore]  # noqa: F401,E501
except ImportError as e:
    greenlet_error = str(e)
    pass
else:
    have_greenlet = True
    from ._concurrency_py3k import await_only as await_only
    from ._concurrency_py3k import await_fallback as await_fallback
    from ._concurrency_py3k import in_greenlet as in_greenlet
    from ._concurrency_py3k import greenlet_spawn as greenlet_spawn
    from ._concurrency_py3k import is_exit_exception as is_exit_exception
    from ._concurrency_py3k import AsyncAdaptedLock as AsyncAdaptedLock
    from ._concurrency_py3k import _Runner

_T = TypeVar("_T")


class _AsyncUtil:
    """Asyncio util for test suite/ util only"""

    def __init__(self) -> None:
        if have_greenlet:
            self.runner = _Runner()

    def run(
        self,
        fn: Callable[..., Coroutine[Any, Any, _T]],
        *args: Any,
        **kwargs: Any,
    ) -> _T:
        """Run coroutine on the loop"""
        return self.runner.run(fn(*args, **kwargs))

    def run_in_greenlet(
        self, fn: Callable[..., _T], *args: Any, **kwargs: Any
    ) -> _T:
        """Run sync function in greenlet. Support nested calls"""
        if have_greenlet:
            if self.runner.get_loop().is_running():
                return fn(*args, **kwargs)
            else:
                return self.runner.run(greenlet_spawn(fn, *args, **kwargs))
        else:
            return fn(*args, **kwargs)

    def close(self) -> None:
        if have_greenlet:
            self.runner.close()


if not typing.TYPE_CHECKING and not have_greenlet:

    def _not_implemented():
        # this conditional is to prevent pylance from considering
        # greenlet_spawn() etc as "no return" and dimming out code below it
        if have_greenlet:
            return None

        raise ValueError(
            "the greenlet library is required to use this function."
            " %s" % greenlet_error
            if greenlet_error
            else ""
        )

    def is_exit_exception(e):  # noqa: F811
        return not isinstance(e, Exception)

    def await_only(thing):  # type: ignore  # noqa: F811
        _not_implemented()

    def await_fallback(thing):  # type: ignore  # noqa: F811
        return thing

    def in_greenlet():  # type: ignore  # noqa: F811
        _not_implemented()

    def greenlet_spawn(fn, *args, **kw):  # type: ignore  # noqa: F811
        _not_implemented()

    def AsyncAdaptedLock(*args, **kw):  # type: ignore  # noqa: F811
        _not_implemented()

    def _util_async_run(fn, *arg, **kw):  # type: ignore  # noqa: F811
        return fn(*arg, **kw)

    def _util_async_run_coroutine_function(fn, *arg, **kw):  # type: ignore  # noqa: F811,E501
        _not_implemented()
