"""
Core Implementations for Generator/BitGenerator Models.
"""

from llvmlite import ir
from numba.core import cgutils, types, config
from numba.core.extending import (intrinsic, make_attribute_wrapper, models,
                                  overload, register_jitable,
                                  register_model)


@register_model(types.NumPyRandomBitGeneratorType)
class NumPyRngBitGeneratorModel(models.StructModel):
    def __init__(self, dmm, fe_type):
        members = [
            ('parent', types.pyobject),
            ('state_address', types.uintp),
            ('state', types.uintp),
            ('fnptr_next_uint64', types.uintp),
            ('fnptr_next_uint32', types.uintp),
            ('fnptr_next_double', types.uintp),
            ('bit_generator', types.uintp),
        ]
        super(NumPyRngBitGeneratorModel, self).__init__(dmm, fe_type, members)


_bit_gen_type = types.NumPyRandomBitGeneratorType('bit_generator')


@register_model(types.NumPyRandomGeneratorType)
class NumPyRandomGeneratorTypeModel(models.StructModel):
    def __init__(self, dmm, fe_type):
        members = [
            ('bit_generator', _bit_gen_type),
            ('meminfo', types.MemInfoPointer(types.voidptr)),
            ('parent', types.pyobject)
        ]
        super(
            NumPyRandomGeneratorTypeModel,
            self).__init__(
            dmm,
            fe_type,
            members)


# The Generator instances have a bit_generator attr
make_attribute_wrapper(
    types.NumPyRandomGeneratorType,
    'bit_generator',
    'bit_generator')


def _generate_next_binding(overloadable_function, return_type):
    """
        Generate the overloads for "next_(some type)" functions.
    """
    @intrinsic
    def intrin_NumPyRandomBitGeneratorType_next_ty(tyctx, inst):
        sig = return_type(inst)

        def codegen(cgctx, builder, sig, llargs):
            name = overloadable_function.__name__
            struct_ptr = cgutils.create_struct_proxy(inst)(cgctx, builder,
                                                           value=llargs[0])

            # Get the 'state' and 'fnptr_next_(type)' members of the struct
            state = struct_ptr.state
            next_double_addr = getattr(struct_ptr, f'fnptr_{name}')

            # LLVM IR types needed
            ll_void_ptr_t = cgctx.get_value_type(types.voidptr)
            ll_return_t = cgctx.get_value_type(return_type)
            ll_uintp_t = cgctx.get_value_type(types.uintp)

            # Convert the stored Generator function address to a pointer
            next_fn_fnptr = builder.inttoptr(
                next_double_addr, ll_void_ptr_t)
            # Add the function to the module
            fnty = ir.FunctionType(ll_return_t, (ll_uintp_t,))
            next_fn = cgutils.get_or_insert_function(
                builder.module, fnty, name)
            # Bit cast the function pointer to the function type
            fnptr_as_fntype = builder.bitcast(next_fn_fnptr, next_fn.type)
            # call it with the "state" address as the arg
            ret = builder.call(fnptr_as_fntype, (state,))
            return ret
        return sig, codegen

    @overload(overloadable_function)
    def ol_next_ty(bitgen):
        if isinstance(bitgen, types.NumPyRandomBitGeneratorType):
            def impl(bitgen):
                return intrin_NumPyRandomBitGeneratorType_next_ty(bitgen)
            return impl


# Some function stubs for "next(some type)", these will be overloaded
def next_double(bitgen):
    return bitgen.ctypes.next_double(bitgen.ctypes.state)


def next_uint32(bitgen):
    return bitgen.ctypes.next_uint32(bitgen.ctypes.state)


def next_uint64(bitgen):
    return bitgen.ctypes.next_uint64(bitgen.ctypes.state)


if config.USE_LEGACY_TYPE_SYSTEM:
    _generate_next_binding(next_double, types.double)
    _generate_next_binding(next_uint32, types.uint32)
    _generate_next_binding(next_uint64, types.uint64)

    # See: https://github.com/numpy/numpy/pull/20314
    @register_jitable
    def next_float(bitgen):
        return types.float32(types.float32(next_uint32(bitgen) >> 8)
                             * types.float32(1.0)
                             / types.float32(16777216.0))

else:
    _generate_next_binding(next_double, types.np_double)
    _generate_next_binding(next_uint32, types.np_uint32)
    _generate_next_binding(next_uint64, types.np_uint64)

    # See: https://github.com/numpy/numpy/pull/20314
    @register_jitable
    def next_float(bitgen):
        return types.np_float32(types.np_float32(next_uint32(bitgen) >> 8)
                                * types.np_float32(1.0)
                                / types.np_float32(16777216.0))
