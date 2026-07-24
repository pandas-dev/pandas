from typing import Annotated, Literal

__all__ = "AnnotatedAlias", "GenericType", "LiteralAlias", "UnionAlias"


class _C[T]: ...


# typing._GenericAlias
# NOTE: this is not the same as`types.GenericAlias`!
GenericType = type(_C[None])

# typing._AnnotatedAlias
AnnotatedAlias = type(Annotated[None, None])

# typing._LiteralGenericAlias
LiteralAlias = type(Literal[0])

# typing._UnionGenericAlias
# NOTE: this is not the same as `types.UnionType`!
UnionAlias = type(Literal[0] | None)
