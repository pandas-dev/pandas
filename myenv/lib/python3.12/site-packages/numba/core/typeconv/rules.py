import itertools
from .typeconv import TypeManager, TypeCastingRules
from numba.core import types


default_type_manager = TypeManager()


def dump_number_rules():
    tm = default_type_manager
    for a, b in itertools.product(types.number_domain, types.number_domain):
        print(a, '->', b, tm.check_compatible(a, b))


def _init_casting_rules(tm):
    tcr = TypeCastingRules(tm)
    tcr.safe_unsafe(types.boolean, types.int8)
    tcr.safe_unsafe(types.boolean, types.uint8)

    tcr.promote_unsafe(types.int8, types.int16)
    tcr.promote_unsafe(types.uint8, types.uint16)

    tcr.promote_unsafe(types.int16, types.int32)
    tcr.promote_unsafe(types.uint16, types.uint32)

    tcr.promote_unsafe(types.int32, types.int64)
    tcr.promote_unsafe(types.uint32, types.uint64)

    tcr.safe_unsafe(types.uint8, types.int16)
    tcr.safe_unsafe(types.uint16, types.int32)
    tcr.safe_unsafe(types.uint32, types.int64)

    tcr.safe_unsafe(types.int8, types.float16)
    tcr.safe_unsafe(types.int16, types.float32)
    tcr.safe_unsafe(types.int32, types.float64)


    tcr.unsafe_unsafe(types.int16, types.float16)
    tcr.unsafe_unsafe(types.int32, types.float32)
    # XXX this is inconsistent with the above; but we want to prefer
    # float64 over int64 when typing a heterogeneous operation,
    # e.g. `float64 + int64`.  Perhaps we need more granularity in the
    # conversion kinds.
    tcr.safe_unsafe(types.int64, types.float64)
    tcr.safe_unsafe(types.uint64, types.float64)

    tcr.promote_unsafe(types.float16, types.float32)
    tcr.promote_unsafe(types.float32, types.float64)

    tcr.safe(types.float32, types.complex64)
    tcr.safe(types.float64, types.complex128)

    tcr.promote_unsafe(types.complex64, types.complex128)

    # Allow integers to cast ot void*
    tcr.unsafe_unsafe(types.uintp, types.voidptr)

    return tcr


default_casting_rules = _init_casting_rules(default_type_manager)

