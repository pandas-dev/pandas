"""
Typing for enums.
"""
import operator
from numba.core import types
from numba.core.typing.templates import (AbstractTemplate, AttributeTemplate,
                                         signature, Registry)

registry = Registry()
infer = registry.register
infer_global = registry.register_global
infer_getattr = registry.register_attr


@infer_getattr
class EnumAttribute(AttributeTemplate):
    key = types.EnumMember

    def resolve_value(self, ty):
        return ty.dtype


@infer_getattr
class EnumClassAttribute(AttributeTemplate):
    key = types.EnumClass

    def generic_resolve(self, ty, attr):
        """
        Resolve attributes of an enum class as enum members.
        """
        if attr in ty.instance_class.__members__:
            return ty.member_type


@infer
class EnumClassStaticGetItem(AbstractTemplate):
    key = "static_getitem"

    def generic(self, args, kws):
        enum, idx = args
        if (isinstance(enum, types.EnumClass)
                and idx in enum.instance_class.__members__):
            return signature(enum.member_type, *args)


class EnumCompare(AbstractTemplate):

    def generic(self, args, kws):
        [lhs, rhs] = args
        if (isinstance(lhs, types.EnumMember)
                and isinstance(rhs, types.EnumMember)
                and lhs == rhs):
            return signature(types.boolean, lhs, rhs)


@infer_global(operator.eq)
class EnumEq(EnumCompare):
    pass



@infer_global(operator.ne)
class EnumNe(EnumCompare):
    pass
