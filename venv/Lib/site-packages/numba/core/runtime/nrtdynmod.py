"""
Dynamically generate the NRT module
"""


from numba.core import config
from numba.core import types, cgutils
from llvmlite import ir, binding


_word_type = ir.IntType(config.MACHINE_BITS)
_pointer_type = ir.PointerType(ir.IntType(8))

_meminfo_struct_type = ir.LiteralStructType([
    _word_type,     # size_t refct
    _pointer_type,  # dtor_function dtor
    _pointer_type,  # void *dtor_info
    _pointer_type,  # void *data
    _word_type,     # size_t size
    ])


incref_decref_ty = ir.FunctionType(ir.VoidType(), [_pointer_type])
meminfo_data_ty = ir.FunctionType(_pointer_type, [_pointer_type])


def _define_nrt_meminfo_data(module):
    """
    Implement NRT_MemInfo_data_fast in the module.  This allows LLVM
    to inline lookup of the data pointer.
    """
    fn = cgutils.get_or_insert_function(module, meminfo_data_ty,
                                        "NRT_MemInfo_data_fast")
    builder = ir.IRBuilder(fn.append_basic_block())
    [ptr] = fn.args
    struct_ptr = builder.bitcast(ptr, _meminfo_struct_type.as_pointer())
    data_ptr = builder.load(cgutils.gep(builder, struct_ptr, 0, 3))
    builder.ret(data_ptr)


def _define_nrt_incref(module, atomic_incr):
    """
    Implement NRT_incref in the module
    """
    fn_incref = cgutils.get_or_insert_function(module, incref_decref_ty,
                                              "NRT_incref")
    # Cannot inline this for refcount pruning to work
    fn_incref.attributes.add('noinline')
    builder = ir.IRBuilder(fn_incref.append_basic_block())
    [ptr] = fn_incref.args
    is_null = builder.icmp_unsigned("==", ptr, cgutils.get_null_value(ptr.type))
    with cgutils.if_unlikely(builder, is_null):
        builder.ret_void()

    word_ptr = builder.bitcast(ptr, atomic_incr.args[0].type)
    if config.DEBUG_NRT:
        cgutils.printf(builder, "*** NRT_Incref %zu [%p]\n", builder.load(word_ptr),
                       ptr)
    builder.call(atomic_incr, [word_ptr])
    builder.ret_void()


def _define_nrt_decref(module, atomic_decr):
    """
    Implement NRT_decref in the module
    """
    fn_decref = cgutils.get_or_insert_function(module, incref_decref_ty,
                                               "NRT_decref")
    # Cannot inline this for refcount pruning to work
    fn_decref.attributes.add('noinline')
    calldtor = ir.Function(module,
                           ir.FunctionType(ir.VoidType(), [_pointer_type]),
                           name="NRT_MemInfo_call_dtor")

    builder = ir.IRBuilder(fn_decref.append_basic_block())
    [ptr] = fn_decref.args
    is_null = builder.icmp_unsigned("==", ptr, cgutils.get_null_value(ptr.type))
    with cgutils.if_unlikely(builder, is_null):
        builder.ret_void()


    # For memory fence usage, see https://llvm.org/docs/Atomics.html

    # A release fence is used before the relevant write operation.
    # No-op on x86.  On POWER, it lowers to lwsync.
    builder.fence("release")

    word_ptr = builder.bitcast(ptr, atomic_decr.args[0].type)

    if config.DEBUG_NRT:
        cgutils.printf(builder, "*** NRT_Decref %zu [%p]\n", builder.load(word_ptr),
                       ptr)
    newrefct = builder.call(atomic_decr,
                            [word_ptr])

    refct_eq_0 = builder.icmp_unsigned("==", newrefct,
                                       ir.Constant(newrefct.type, 0))
    with cgutils.if_unlikely(builder, refct_eq_0):
        # An acquire fence is used after the relevant read operation.
        # No-op on x86.  On POWER, it lowers to lwsync.
        builder.fence("acquire")
        builder.call(calldtor, [ptr])
    builder.ret_void()


# Set this to True to measure the overhead of atomic refcounts compared
# to non-atomic.
_disable_atomicity = 0


def _define_atomic_inc_dec(module, op, ordering):
    """Define a llvm function for atomic increment/decrement to the given module
    Argument ``op`` is the operation "add"/"sub".  Argument ``ordering`` is
    the memory ordering.  The generated function returns the new value.
    """
    ftype = ir.FunctionType(_word_type, [_word_type.as_pointer()])
    fn_atomic = ir.Function(module, ftype, name="nrt_atomic_{0}".format(op))

    [ptr] = fn_atomic.args
    bb = fn_atomic.append_basic_block()
    builder = ir.IRBuilder(bb)
    ONE = ir.Constant(_word_type, 1)
    if not _disable_atomicity:
        oldval = builder.atomic_rmw(op, ptr, ONE, ordering=ordering)
        # Perform the operation on the old value so that we can pretend returning
        # the "new" value.
        res = getattr(builder, op)(oldval, ONE)
        builder.ret(res)
    else:
        oldval = builder.load(ptr)
        newval = getattr(builder, op)(oldval, ONE)
        builder.store(newval, ptr)
        builder.ret(oldval)

    return fn_atomic


def _define_atomic_cas(module, ordering):
    """Define a llvm function for atomic compare-and-swap.
    The generated function is a direct wrapper of the LLVM cmpxchg with the
    difference that the a int indicate success (1) or failure (0) is returned
    and the last argument is a output pointer for storing the old value.

    Note
    ----
    On failure, the generated function behaves like an atomic load.  The loaded
    value is stored to the last argument.
    """
    ftype = ir.FunctionType(ir.IntType(32), [_word_type.as_pointer(),
                                             _word_type, _word_type,
                                             _word_type.as_pointer()])
    fn_cas = ir.Function(module, ftype, name="nrt_atomic_cas")

    [ptr, cmp, repl, oldptr] = fn_cas.args
    bb = fn_cas.append_basic_block()
    builder = ir.IRBuilder(bb)
    outtup = builder.cmpxchg(ptr, cmp, repl, ordering=ordering)
    old, ok = cgutils.unpack_tuple(builder, outtup, 2)
    builder.store(old, oldptr)
    builder.ret(builder.zext(ok, ftype.return_type))

    return fn_cas


def _define_nrt_unresolved_abort(ctx, module):
    """
    Defines an abort function due to unresolved symbol.

    The function takes no args and will always raise an exception.
    It should be safe to call this function with incorrect number of arguments.
    """
    fnty = ctx.call_conv.get_function_type(types.none, ())
    fn = ir.Function(module, fnty, name="nrt_unresolved_abort")
    bb = fn.append_basic_block()
    builder = ir.IRBuilder(bb)
    msg = "numba jitted function aborted due to unresolved symbol"
    ctx.call_conv.return_user_exc(builder, RuntimeError, (msg,))
    return fn


def create_nrt_module(ctx):
    """
    Create an IR module defining the LLVM NRT functions.
    A (IR module, library) tuple is returned.
    """
    codegen = ctx.codegen()
    library = codegen.create_library("nrt")

    # Implement LLVM module with atomic ops
    ir_mod = library.create_ir_module("nrt_module")

    atomic_inc = _define_atomic_inc_dec(ir_mod, "add", ordering='monotonic')
    atomic_dec = _define_atomic_inc_dec(ir_mod, "sub", ordering='monotonic')
    _define_atomic_cas(ir_mod, ordering='monotonic')

    _define_nrt_meminfo_data(ir_mod)
    _define_nrt_incref(ir_mod, atomic_inc)
    _define_nrt_decref(ir_mod, atomic_dec)

    _define_nrt_unresolved_abort(ctx, ir_mod)

    return ir_mod, library


def compile_nrt_functions(ctx):
    """
    Compile all LLVM NRT functions and return a library containing them.
    The library is created using the given target context.
    """
    ir_mod, library = create_nrt_module(ctx)

    library.add_ir_module(ir_mod)
    library.finalize()

    return library
