# Copyright 2022 Amethyst Reese
# Licensed under the MIT license

import inspect
from collections.abc import Awaitable

from typing import Protocol, Union

from .types import T


class Orderable(Protocol):  # pragma: no cover
    def __lt__(self, other): ...

    def __gt__(self, other): ...


async def maybe_await(object: Union[Awaitable[T], T]) -> T:
    if inspect.isawaitable(object):
        return await object  # type: ignore
    return object  # type: ignore
