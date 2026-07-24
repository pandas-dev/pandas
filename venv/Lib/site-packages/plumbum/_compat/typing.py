from __future__ import annotations

import sys
import typing

if sys.version_info >= (3, 11):
    from typing import Self
elif typing.TYPE_CHECKING:
    from typing_extensions import Self
else:
    Self = object()

if sys.version_info >= (3, 13):
    from typing import TypeVar
elif typing.TYPE_CHECKING:
    from typing_extensions import TypeVar
else:

    def TypeVar(
        name: str,
        *args: object,
        default: object = None,  # noqa: ARG001
        **kwargs: object,
    ) -> typing.TypeVar:
        return typing.TypeVar(name, *args, **kwargs)


__all__ = ["Self", "TypeVar"]


def __dir__() -> list[str]:
    return list(__all__)
