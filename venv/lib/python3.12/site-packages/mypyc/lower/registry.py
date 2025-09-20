from __future__ import annotations

from typing import Callable, Final, Optional, TypeVar

from mypyc.ir.ops import Value
from mypyc.irbuild.ll_builder import LowLevelIRBuilder

LowerFunc = Callable[[LowLevelIRBuilder, list[Value], int], Value]
LowerFuncOpt = Callable[[LowLevelIRBuilder, list[Value], int], Optional[Value]]

lowering_registry: Final[dict[str, LowerFuncOpt]] = {}

LF = TypeVar("LF", LowerFunc, LowerFuncOpt)


def lower_primitive_op(name: str) -> Callable[[LF], LF]:
    """Register a handler that generates low-level IR for a primitive op."""

    def wrapper(f: LF) -> LF:
        assert name not in lowering_registry
        lowering_registry[name] = f
        return f

    return wrapper


# Import various modules that set up global state.
from mypyc.lower import int_ops, list_ops, misc_ops  # noqa: F401
