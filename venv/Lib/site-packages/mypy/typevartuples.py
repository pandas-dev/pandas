"""Helpers for interacting with type var tuples."""

from __future__ import annotations

from collections.abc import Sequence

from mypy.types import (
    AnyType,
    Instance,
    Type,
    TypeVarLikeType,
    TypeVarTupleType,
    UnpackType,
    split_with_prefix_and_suffix,
)


def split_with_instance(
    typ: Instance,
) -> tuple[tuple[Type, ...], tuple[Type, ...], tuple[Type, ...]]:
    assert typ.type.type_var_tuple_prefix is not None
    assert typ.type.type_var_tuple_suffix is not None
    return split_with_prefix_and_suffix(
        typ.args, typ.type.type_var_tuple_prefix, typ.type.type_var_tuple_suffix
    )


def erased_vars(type_vars: Sequence[TypeVarLikeType], type_of_any: int) -> list[Type]:
    args: list[Type] = []
    for tv in type_vars:
        # Valid erasure for *Ts is *tuple[Any, ...], not just Any.
        if isinstance(tv, TypeVarTupleType):
            args.append(UnpackType(tv.tuple_fallback.copy_modified(args=[AnyType(type_of_any)])))
        else:
            args.append(AnyType(type_of_any))
    return args
