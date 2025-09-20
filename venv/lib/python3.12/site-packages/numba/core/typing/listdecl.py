import operator
from numba.core import types
from .templates import (ConcreteTemplate, AbstractTemplate, AttributeTemplate,
                        CallableTemplate,  Registry, signature, bound_function,
                        make_callable_template)
# Ensure list is typed as a collection as well
from numba.core.typing import collections


registry = Registry()
infer = registry.register
infer_global = registry.register_global
infer_getattr = registry.register_attr


@infer_global(list)
class ListBuiltin(AbstractTemplate):

    def generic(self, args, kws):
        assert not kws
        if args:
            iterable, = args
            if isinstance(iterable, types.IterableType):
                dtype = iterable.iterator_type.yield_type
                return signature(types.List(dtype), iterable)
        else:
            return signature(types.List(types.undefined))


@infer_getattr
class ListAttribute(AttributeTemplate):
    key = types.List

    # NOTE: some of these should be Sequence / MutableSequence methods

    @bound_function("list.append")
    def resolve_append(self, list, args, kws):
        item, = args
        assert not kws
        unified = self.context.unify_pairs(list.dtype, item)
        if unified is not None:
            sig = signature(types.none, unified)
            sig = sig.replace(recvr=list.copy(dtype=unified))
            return sig

    @bound_function("list.clear")
    def resolve_clear(self, list, args, kws):
        assert not args
        assert not kws
        return signature(types.none)

    @bound_function("list.extend")
    def resolve_extend(self, list, args, kws):
        iterable, = args
        assert not kws
        if not isinstance(iterable, types.IterableType):
            return

        dtype = iterable.iterator_type.yield_type
        unified = self.context.unify_pairs(list.dtype, dtype)
        if unified is not None:
            sig = signature(types.none, iterable)
            sig = sig.replace(recvr = list.copy(dtype=unified))
            return sig

    @bound_function("list.insert")
    def resolve_insert(self, list, args, kws):
        idx, item = args
        assert not kws
        if isinstance(idx, types.Integer):
            unified = self.context.unify_pairs(list.dtype, item)
            if unified is not None:
                sig = signature(types.none, types.intp, unified)
                sig = sig.replace(recvr = list.copy(dtype=unified))
                return sig

    @bound_function("list.pop")
    def resolve_pop(self, list, args, kws):
        assert not kws
        if not args:
            return signature(list.dtype)
        else:
            idx, = args
            if isinstance(idx, types.Integer):
                return signature(list.dtype, types.intp)

@infer_global(operator.add)
class AddList(AbstractTemplate):

    def generic(self, args, kws):
        if len(args) == 2:
            a, b = args
            if isinstance(a, types.List) and isinstance(b, types.List):
                unified = self.context.unify_pairs(a, b)
                if unified is not None:
                    return signature(unified, a, b)


@infer_global(operator.iadd)
class InplaceAddList(AbstractTemplate):

    def generic(self, args, kws):
        if len(args) == 2:
            a, b = args
            if isinstance(a, types.List) and isinstance(b, types.List):
                if self.context.can_convert(b.dtype, a.dtype):
                    return signature(a, a, b)


@infer_global(operator.mul)
class MulList(AbstractTemplate):
    #key = operator.mul

    def generic(self, args, kws):
        a, b = args
        if isinstance(a, types.List) and isinstance(b, types.Integer):
            return signature(a, a, types.intp)
        elif isinstance(a, types.Integer) and isinstance(b, types.List):
            return signature(b, types.intp, b)


@infer_global(operator.imul)
class InplaceMulList(MulList): pass
    #key = operator.imul


class ListCompare(AbstractTemplate):

    def generic(self, args, kws):
        [lhs, rhs] = args
        if isinstance(lhs, types.List) and isinstance(rhs, types.List):
            # Check element-wise comparability
            res = self.context.resolve_function_type(self.key,
                                                     (lhs.dtype, rhs.dtype), {})
            if res is not None:
                return signature(types.boolean, lhs, rhs)

@infer_global(operator.eq)
class ListEq(ListCompare): pass
    #key = operator.eq
