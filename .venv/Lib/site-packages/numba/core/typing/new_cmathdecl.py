import cmath

from numba.core import types, utils
from numba.core.typing.templates import (AbstractTemplate, ConcreteTemplate,
                                    signature, Registry)

registry = Registry()
infer_global = registry.register_global

# TODO: support non-complex arguments (floats and ints)

# TODO: New Type System
# These functions are part of the Python standard library
#  and (without checking) probably accept anything which
#  is "number"-like i.e. has a __float__, __int__, or
#  __index__
# This needs fixing in the new type system


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
    cases = []


@infer_global(cmath.isinf)
@infer_global(cmath.isnan)
class CMath_predicate(ConcreteTemplate):
    cases = []


@infer_global(cmath.isfinite)
class CMath_isfinite(CMath_predicate):
    pass


@infer_global(cmath.log)
class Cmath_log(ConcreteTemplate):
    # unary cmath.log()
    cases = []
    # binary cmath.log()
    cases += []
