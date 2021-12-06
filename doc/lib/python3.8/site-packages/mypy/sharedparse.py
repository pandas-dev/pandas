from typing import Optional
from typing_extensions import Final

"""Shared logic between our three mypy parser files."""


_NON_BINARY_MAGIC_METHODS = {
    "__abs__",
    "__call__",
    "__complex__",
    "__contains__",
    "__del__",
    "__delattr__",
    "__delitem__",
    "__enter__",
    "__exit__",
    "__float__",
    "__getattr__",
    "__getattribute__",
    "__getitem__",
    "__hex__",
    "__init__",
    "__init_subclass__",
    "__int__",
    "__invert__",
    "__iter__",
    "__len__",
    "__long__",
    "__neg__",
    "__new__",
    "__nonzero__",
    "__oct__",
    "__pos__",
    "__repr__",
    "__reversed__",
    "__setattr__",
    "__setitem__",
    "__str__",
    "__unicode__",
}  # type: Final

MAGIC_METHODS_ALLOWING_KWARGS = {
    "__init__",
    "__init_subclass__",
    "__new__",
    "__call__",
}  # type: Final

BINARY_MAGIC_METHODS = {
    "__add__",
    "__and__",
    "__cmp__",
    "__divmod__",
    "__div__",
    "__eq__",
    "__floordiv__",
    "__ge__",
    "__gt__",
    "__iadd__",
    "__iand__",
    "__idiv__",
    "__ifloordiv__",
    "__ilshift__",
    "__imatmul__",
    "__imod__",
    "__imul__",
    "__ior__",
    "__ipow__",
    "__irshift__",
    "__isub__",
    "__itruediv__",
    "__ixor__",
    "__le__",
    "__lshift__",
    "__lt__",
    "__matmul__",
    "__mod__",
    "__mul__",
    "__ne__",
    "__or__",
    "__pow__",
    "__radd__",
    "__rand__",
    "__rdiv__",
    "__rfloordiv__",
    "__rlshift__",
    "__rmatmul__",
    "__rmod__",
    "__rmul__",
    "__ror__",
    "__rpow__",
    "__rrshift__",
    "__rshift__",
    "__rsub__",
    "__rtruediv__",
    "__rxor__",
    "__sub__",
    "__truediv__",
    "__xor__",
}  # type: Final

assert not (_NON_BINARY_MAGIC_METHODS & BINARY_MAGIC_METHODS)

MAGIC_METHODS = _NON_BINARY_MAGIC_METHODS | BINARY_MAGIC_METHODS  # type: Final

MAGIC_METHODS_POS_ARGS_ONLY = MAGIC_METHODS - MAGIC_METHODS_ALLOWING_KWARGS  # type: Final


def special_function_elide_names(name: str) -> bool:
    return name in MAGIC_METHODS_POS_ARGS_ONLY


def argument_elide_name(name: Optional[str]) -> bool:
    return name is not None and name.startswith("__") and not name.endswith("__")
