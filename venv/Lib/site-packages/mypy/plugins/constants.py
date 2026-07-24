"""Constant definitions for plugins kept here to help with import cycles."""

from typing import Final

from mypy.semanal_enum import ENUM_BASES

SINGLEDISPATCH_TYPE: Final = "functools._SingleDispatchCallable"
SINGLEDISPATCH_REGISTER_METHOD: Final = f"{SINGLEDISPATCH_TYPE}.register"
SINGLEDISPATCH_CALLABLE_CALL_METHOD: Final = f"{SINGLEDISPATCH_TYPE}.__call__"
SINGLEDISPATCH_REGISTER_RETURN_CLASS: Final = "_SingleDispatchRegisterCallable"
SINGLEDISPATCH_REGISTER_CALLABLE_CALL_METHOD: Final = (
    f"functools.{SINGLEDISPATCH_REGISTER_RETURN_CLASS}.__call__"
)

ENUM_NAME_ACCESS: Final = {f"{prefix}.name" for prefix in ENUM_BASES} | {
    f"{prefix}._name_" for prefix in ENUM_BASES
}
ENUM_VALUE_ACCESS: Final = {f"{prefix}.value" for prefix in ENUM_BASES} | {
    f"{prefix}._value_" for prefix in ENUM_BASES
}
