"""
Reversed Operations not available in the stdlib operator module.
Defining these instead of using lambdas allows us to reference them by name.
"""

from __future__ import annotations

import operator
from typing import Any


def radd(left: Any, right: Any) -> Any:
    return right + left


def rsub(left: Any, right: Any) -> Any:
    return right - left


def rmul(left: Any, right: Any) -> Any:
    return right * left


def rdiv(left: Any, right: Any) -> Any:
    return right / left


def rtruediv(left: Any, right: Any) -> Any:
    return right / left


def rfloordiv(left: Any, right: Any) -> Any:
    return right // left


def rmod(left: Any, right: Any) -> Any:
    # check if right is a string as % is the string
    # formatting operation; this is a TypeError
    # otherwise perform the op
    if isinstance(right, str):
        typ = type(left).__name__
        raise TypeError(f"{typ} cannot perform the operation mod")

    return right % left


def rdivmod(left: Any, right: Any) -> tuple[Any, Any]:
    return divmod(right, left)


def rpow(left: Any, right: Any) -> Any:
    return right**left


def rand_(left: Any, right: Any) -> Any:
    return operator.and_(right, left)


def ror_(left: Any, right: Any) -> Any:
    return operator.or_(right, left)


def rxor(left: Any, right: Any) -> Any:
    return operator.xor(right, left)
