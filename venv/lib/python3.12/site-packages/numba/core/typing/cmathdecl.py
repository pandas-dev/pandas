import cmath

from numba.core import types, utils
from numba.core.typing.templates import (AbstractTemplate, ConcreteTemplate,
                                    signature, Registry)

registry = Registry()
infer_global = registry.register_global

# TODO: support non-complex arguments (floats and ints)


@infer_global(cmath.acos)
@infer_global(cmath.asin)
@infer_global(cmath.asinh)
@infer_global(cmath.atan)
@infer_global(cmath.atanh)
@infer_global(cmath.cos)
@infer_global(cmath.exp)
@infer_global(cmath.sin)
@infer_global(cmath.sqrt)
@infer_global(cmath.tan)
class CMath_unary(ConcreteTemplate):
    cases = [signature(tp, tp) for tp in sorted(types.complex_domain)]


@infer_global(cmath.isinf)
@infer_global(cmath.isnan)
class CMath_predicate(ConcreteTemplate):
    cases = [signature(types.boolean, tp) for tp in
             sorted(types.complex_domain)]


@infer_global(cmath.isfinite)
class CMath_isfinite(CMath_predicate):
    pass


@infer_global(cmath.log)
class Cmath_log(ConcreteTemplate):
    # unary cmath.log()
    cases = [signature(tp, tp) for tp in sorted(types.complex_domain)]
    # binary cmath.log()
    cases += [signature(tp, tp, tp) for tp in sorted(types.complex_domain)]
