from typing import Final, Literal, LiteralString, overload

__all__ = ["ConstantWarning", "find", "physical_constants", "precision", "unit", "value"]

# private

c: Final = 299792458.0
k: Final = "electric constant"

# public

physical_constants: Final[dict[str, tuple[float, str, float]]] = ...

class ConstantWarning(DeprecationWarning): ...

#
def value(key: str) -> float: ...
def unit(key: str) -> LiteralString: ...
def precision(key: str) -> float: ...

#
@overload  # disp=False (default)
def find(sub: str | None = None, disp: Literal[False, 0] = False) -> list[str]: ...
@overload  # disp=True (positional)
def find(sub: str | None, disp: Literal[True, 1]) -> None: ...
@overload  # disp=True (keyword)
def find(sub: str | None = None, *, disp: Literal[True, 1]) -> None: ...
