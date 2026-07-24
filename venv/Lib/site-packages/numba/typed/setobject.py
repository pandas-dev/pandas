"""
Compiler-side implementation of the sets.
"""
import ctypes
import operator

from enum import IntEnum

from llvmlite import ir
from numba import _helperlib

from numba.core.extending import (
    overload,
    overload_method,
    intrinsic,
    register_model,
    models,
    lower_builtin,
    make_attribute_wrapper,
)
from numba.core.imputils import iternext_impl, RefType
from numba.core import types, cgutils, typing
from numba.core.types import (
    SetType,
    SetIterableType,
    SetIteratorType,
    Type
)
from numba.core.errors import TypingError
from numba.typed.typedobjectutils import (_get_incref_decref,
                                          _container_get_data,
                                          _as_bytes, _cast,
                                          _get_equal)
from numba.core.imputils import impl_ret_borrowed


ll_set_type = cgutils.voidptr_t
ll_setiter_type = cgutils.voidptr_t
ll_voidptr_type = cgutils.voidptr_t
ll_status = cgutils.int32_t
ll_ssize_t = cgutils.intp_t
ll_hash = ll_ssize_t
ll_bytes = cgutils.voidptr_t


_meminfo_setptr = types.MemInfoPointer(types.voidptr)

# The following enums must match _setobject.c


class Status(IntEnum):
    """Status code for set operations.
    """
    ENTRY_PRESENT = 1
    OK = 0
    ERR_KEY_NOT_FOUND = -1
    ERR_SET_MUTATED = -2
    ERR_ITER_EXHAUSTED = -3
    ERR_SET_EMPTY = -4
    ERR_CMP_FAILED = -5


def new_set(key):
    """Construct a new set.

    Parameters
    ----------
    key : TypeRef
        Key type of the new set.
    """
    # With JIT disabled, ignore all arguments and return a Python set.
    return set()


@register_model(SetType)
class SetModel(models.StructModel):
    def __init__(self, dmm, fe_type):
        members = [
            ('meminfo', _meminfo_setptr),
            ('data', types.voidptr),   # ptr to the C set
        ]
        super(SetModel, self).__init__(dmm, fe_type, members)


@register_model(SetIterableType)
@register_model(SetIteratorType)
class SetIterModel(models.StructModel):
    def __init__(self, dmm, fe_type):
        members = [
            ('parent', fe_type.parent),  # reference to the set
            ('state', types.voidptr),    # iterator state in C code
        ]
        super(SetIterModel, self).__init__(dmm, fe_type, members)


# Make _parent available to make len simple
make_attribute_wrapper(SetIteratorType, "parent", "_parent")
make_attribute_wrapper(SetIterableType, "parent", "_parent")


def _raise_if_error(context, builder, status, msg):
    """Raise an internal error depending on the value of *status*
    """
    ok_status = status.type(int(Status.OK))
    with builder.if_then(builder.icmp_signed('!=', status, ok_status)):
        context.call_conv.return_user_exc(builder, RuntimeError, (msg,))


@intrinsic
def _as_meminfo(typingctx, setobj):
    """Returns the MemInfoPointer of a set.
    """
    if not isinstance(setobj, types.SetType):
        raise TypingError('expected *setobj* to be a SetType')

    def codegen(context, builder, sig, args):
        [ts] = sig.args
        [s] = args
        # Incref
        context.nrt.incref(builder, ts, s)
        ctor = cgutils.create_struct_proxy(ts)
        sstruct = ctor(context, builder, value=s)
        # Returns the plain MemInfo
        return sstruct.meminfo

    sig = _meminfo_setptr(setobj)
    return sig, codegen


@intrinsic
def _set_new_sized(typingctx, n_keys, keyty):
    """Wrap numba_set_new.

    Allocate a new set object.

    Parameters
    ----------
    keyty: Type
        Type of key, respectively.

    """
    resty = types.voidptr
    sig = resty(n_keys, keyty)

    def codegen(context, builder, sig, args):
        n_keys = builder.bitcast(args[0], ll_ssize_t)

        # Determine sizeof value type
        ll_val = context.get_data_type(keyty.instance_type)
        sz_val = context.get_abi_sizeof(ll_val)
        refsetp = cgutils.alloca_once(builder, ll_set_type, zfill=True)

        fnty = ir.FunctionType(
            ll_status,
            [ll_set_type.as_pointer(), ll_ssize_t, ll_ssize_t],
        )
        fn = cgutils.get_or_insert_function(builder.module, fnty,
                                            'numba_set_new_sized')

        status = builder.call(
            fn,
            [refsetp, ll_ssize_t(sz_val), n_keys],
        )
        _raise_if_error(
            context, builder, status,
            msg="Failed to allocate set",
        )
        setp = builder.load(refsetp)
        return setp

    return sig, codegen


def _call_set_free(context, builder, ptr):
    """Call numba_set_free(ptr)
    """
    fnty = ir.FunctionType(
        ir.VoidType(),
        [ll_set_type],
    )
    free = cgutils.get_or_insert_function(builder.module, fnty,
                                          'numba_set_free')
    builder.call(free, [ptr])


def _imp_dtor(context, module):
    """Define the dtor for set
    """
    llvoidptr = context.get_value_type(types.voidptr)
    llsize = context.get_value_type(types.uintp)
    fnty = ir.FunctionType(
        ir.VoidType(),
        [llvoidptr, llsize, llvoidptr],
    )
    fname = '_numba_set_dtor'
    fn = cgutils.get_or_insert_function(module, fnty, fname)

    if fn.is_declaration:
        # Set linkage
        fn.linkage = 'linkonce_odr'
        # Define
        builder = ir.IRBuilder(fn.append_basic_block())
        setp = builder.bitcast(fn.args[0], ll_set_type.as_pointer())
        s = builder.load(setp)
        _call_set_free(context, builder, s)
        builder.ret_void()

    return fn


@intrinsic
def _make_set(typingctx, keyty, ptr):
    """Make a set struct with the given *ptr*

    Parameters
    ----------
    keyty: Type
        Type of the key, respectively.
    ptr : llvm pointer value
        Points to the set object.
    """
    set_ty = types.SetType(keyty.instance_type)

    def codegen(context, builder, signature, args):
        [_, ptr] = args
        ctor = cgutils.create_struct_proxy(set_ty)
        sstruct = ctor(context, builder)
        sstruct.data = ptr

        alloc_size = context.get_abi_sizeof(
            context.get_value_type(types.voidptr),
        )
        dtor = _imp_dtor(context, builder.module)
        meminfo = context.nrt.meminfo_alloc_dtor(
            builder,
            context.get_constant(types.uintp, alloc_size),
            dtor,
        )

        data_pointer = context.nrt.meminfo_data(builder, meminfo)
        data_pointer = builder.bitcast(data_pointer, ll_set_type.as_pointer())
        builder.store(ptr, data_pointer)

        sstruct.meminfo = meminfo

        return sstruct._getvalue()

    sig = set_ty(keyty, ptr)
    return sig, codegen


@intrinsic
def _set_set_method_table(typingctx, setp, keyty):
    """Wrap numba_set_set_method_table
    """
    resty = types.void
    sig = resty(setp, keyty)

    def codegen(context, builder, sig, args):
        tablety = ir.LiteralStructType([
            ll_voidptr_type,  # equal
            ll_voidptr_type,  # key incref
            ll_voidptr_type,  # key decref
        ])
        setmethod_fnty = ir.FunctionType(
            ir.VoidType(),
            [ll_set_type, tablety.as_pointer()]
        )
        setmethod_fn = ir.Function(
            builder.module,
            setmethod_fnty,
            name='numba_set_set_method_table',
        )
        setp = args[0]
        table = cgutils.alloca_once(builder, tablety, zfill=True)

        key_equal_ptr = cgutils.gep_inbounds(builder, table, 0, 0)
        key_incref_ptr = cgutils.gep_inbounds(builder, table, 0, 1)
        key_decref_ptr = cgutils.gep_inbounds(builder, table, 0, 2)

        dm_key = context.data_model_manager[keyty.instance_type]
        if dm_key.contains_nrt_meminfo():
            equal = _get_equal(context, builder.module, dm_key, 'set_key')
            key_incref, key_decref = _get_incref_decref(
                context, builder.module, dm_key, 'set_key'
            )
            builder.store(
                builder.bitcast(equal, key_equal_ptr.type.pointee),
                key_equal_ptr,
            )
            builder.store(
                builder.bitcast(key_incref, key_incref_ptr.type.pointee),
                key_incref_ptr,
            )
            builder.store(
                builder.bitcast(key_decref, key_decref_ptr.type.pointee),
                key_decref_ptr,
            )

        builder.call(setmethod_fn, [setp, table])

    return sig, codegen


@overload(new_set)
def impl_new_set(value, n_keys=8):
    """Creates a new set with *value* as the type
    of the set value.
    """
    if not isinstance(value, Type):
        raise TypeError("expecting *value* to be a numba Type")

    keyty = value

    def imp(value, n_keys=8):
        setp = _set_new_sized(n_keys, keyty)
        _set_set_method_table(setp, keyty)
        s = _make_set(keyty, setp)
        return s

    return imp


@intrinsic
def _set_length(typingctx, s):
    """Wrap numba_set_length

    Returns the length of the set.
    """
    resty = types.intp
    sig = resty(s)

    def codegen(context, builder, sig, args):
        fnty = ir.FunctionType(
            ll_ssize_t,
            [ll_set_type],
        )
        fn = cgutils.get_or_insert_function(builder.module, fnty,
                                            'numba_set_length')
        [s] = args
        [ts] = sig.args
        setp = _container_get_data(context, builder, ts, s)
        n = builder.call(fn, [setp])
        return n

    return sig, codegen


@overload(len)
def impl_len(setp):
    """len(set)
    """
    if not isinstance(setp, types.SetType):
        return

    def impl(setp):
        return _set_length(setp)

    return impl


@overload(len)
def impl_len_iters(s):
    if not isinstance(s, (SetIterableType)):
        return

    def impl(s):
        return _set_length(s._parent)

    return impl


@intrinsic
def _set_add(typingctx, setp, key, hash_key):
    """Wrap numba_set_add
    """
    resty = types.int32
    sig = resty(setp, setp.key_type, types.intp)

    def codegen(context, builder, sig, args):
        fnty = ir.FunctionType(
            ll_status,
            [ll_set_type, ll_bytes, ll_hash],
        )
        [s, key, hash_key] = args
        [set_ty, keyty, hash_key_ty] = sig.args
        fn = cgutils.get_or_insert_function(builder.module, fnty,
                                            'numba_set_add')

        key_datamodel = context.data_model_manager[keyty]

        key_data = key_datamodel.as_data(builder, key)
        ptr_key = cgutils.alloca_once_value(builder, key_data)
        cgutils.memset_padding(builder, ptr_key)

        set_ptr = _container_get_data(context, builder, set_ty, s)
        status = builder.call(
            fn,
            [
                set_ptr,
                _as_bytes(builder, ptr_key),
                hash_key
            ],
        )
        return status

    return sig, codegen


@intrinsic
def _set_contains(typingctx, s, key, hash_key):
    """Wrap numba_set_contains

    Returns intp
    """
    resty = types.intp
    sig = resty(s, key, hash_key)

    def codegen(context, builder, sig, args):
        fnty = ir.FunctionType(
            ll_ssize_t,
            [ll_set_type, ll_bytes, ll_hash],
        )
        [set_ty, keyty, hash_keyty] = sig.args
        [s, key, hash_key] = args
        fn = cgutils.get_or_insert_function(builder.module, fnty,
                                            'numba_set_contains')

        dm_key = context.data_model_manager[keyty]

        data_key = dm_key.as_data(builder, key)
        ptr_key = cgutils.alloca_once_value(builder, data_key)
        cgutils.memset_padding(builder, ptr_key)

        setp = _container_get_data(context, builder, set_ty, s)
        ix = builder.call(
            fn,
            [
                setp,
                _as_bytes(builder, ptr_key),
                hash_key
            ],
        )
        return ix

    return sig, codegen


@intrinsic
def _from_meminfo(typingctx, mi, settyperef):
    """Recreate a set from a MemInfoPointer
    """
    if mi != _meminfo_setptr:
        raise TypingError('expected a MemInfoPointer for set.')
    settype = settyperef.instance_type
    if not isinstance(settype, SetType):
        raise TypingError('expected a {}'.format(SetType))

    def codegen(context, builder, sig, args):
        [tmi, tdref] = sig.args
        setty = tdref.instance_type
        [mi, _] = args

        ctor = cgutils.create_struct_proxy(setty)
        sstruct = ctor(context, builder)

        data_pointer = context.nrt.meminfo_data(builder, mi)
        data_pointer = builder.bitcast(data_pointer, ll_set_type.as_pointer())

        sstruct.data = builder.load(data_pointer)
        sstruct.meminfo = mi

        return impl_ret_borrowed(
            context,
            builder,
            settype,
            sstruct._getvalue(),
        )

    sig = settype(mi, settyperef)
    return sig, codegen


@intrinsic
def _set_discard(typingctx, s, key, hash_key):
    """Wrap numba_set_discard
    """
    resty = types.int32
    sig = resty(s, s.key_type, types.intp)

    def codegen(context, builder, sig, args):
        fnty = ir.FunctionType(
            ll_status,
            [ll_set_type, ll_bytes, ll_hash],
        )
        [s, key, hash_key] = args
        [setp_ty, keyty, hash_key_ty] = sig.args
        fn = cgutils.get_or_insert_function(builder.module, fnty,
                                            'numba_set_discard')

        key_datamodel = context.data_model_manager[keyty]

        key_data = key_datamodel.as_data(builder, key)

        ptr_key = cgutils.alloca_once_value(builder, key_data)
        cgutils.memset_padding(builder, ptr_key)

        setp = _container_get_data(context, builder, setp_ty, s)
        status = builder.call(
            fn,
            [
                setp,
                _as_bytes(builder, ptr_key),
                hash_key
            ],
        )
        return status

    return sig, codegen


def _iterator_codegen(resty):
    """The common codegen for iterator intrinsics.

    Populates the iterator struct and increfs.
    """

    def codegen(context, builder, sig, args):
        [s] = args
        [set_ty] = sig.args
        iterhelper = context.make_helper(builder, resty)
        iterhelper.parent = s
        iterhelper.state = iterhelper.state.type(None)
        return impl_ret_borrowed(
            context,
            builder,
            resty,
            iterhelper._getvalue(),
        )

    return codegen


@overload_method(types.SetType, 'add')
def impl_set_add(s, key):
    if not isinstance(s, types.SetType):
        return

    keyty = s.key_type

    def impl(s, key):
        casted_key = _cast(key, keyty)
        status = _set_add(s, casted_key, hash(casted_key))
        if status == 0:
            return

    if s.is_precise():
        # Handle the precise case.
        return impl
    else:
        # Handle the imprecise case.
        s = s.refine(key)
        # Re-bind the key type to match the arguments.
        keyty = s.key_type
        # Create the signature that we wanted this impl to have.
        sig = typing.signature(types.void, s, keyty)
        return sig, impl


@overload_method(types.SetType, 'discard')
def impl_discard(s, key):
    if not isinstance(s, types.SetType):
        return

    keyty = s.key_type

    def impl(s, key):
        castedkey = _cast(key, keyty)
        return _set_discard(s, castedkey, hash(castedkey))

    return impl


@overload(operator.contains)
def impl_contains(s, key):
    if not isinstance(s, types.SetType):
        return

    keyty = s.key_type

    def impl(s, key):
        key = _cast(key, keyty)
        ix = _set_contains(s, key, hash(key))
        return ix
    return impl


@overload(operator.eq)
def impl_equal(set_a, set_b):
    if not isinstance(set_a, types.SetType):
        return
    if not isinstance(set_b, types.SetType):
        # If RHS is not a set, always returns False
        def impl_type_mismatch(set_a, set_b):
            return False
        return impl_type_mismatch

    otherkeyty = set_b.key_type

    def impl_type_matched(set_a, set_b):
        if len(set_a) != len(set_b):
            return False
        for element in set_a:
            kb = _cast(element, otherkeyty)
            ix = _set_contains(set_b, kb, hash(kb))
            if ix == 0:
                return False
        return True

    return impl_type_matched


@overload(operator.ne)
def impl_not_equal(da, db):
    if not isinstance(da, types.SetType):
        return

    def impl(da, db):
        return not (da == db)

    return impl


@overload_method(types.SetType, 'copy')
def impl_copy(s):
    if not isinstance(s, types.SetType):
        return

    key_type = s.key_type

    def impl(s):
        newd = new_set(key_type, n_keys=len(s))
        for k in s:
            newd.add(k)
        return newd

    return impl


@lower_builtin('getiter', types.SetIterableType)
def impl_iterable_getiter(context, builder, sig, args):
    iterablety = sig.args[0]
    it = context.make_helper(builder, iterablety.iterator_type, args[0])

    fnty = ir.FunctionType(
        ir.VoidType(),
        [ll_setiter_type, ll_set_type],
    )

    fn = cgutils.get_or_insert_function(builder.module, fnty,
                                        'numba_set_iter')

    proto = ctypes.CFUNCTYPE(ctypes.c_size_t)
    setiter_sizeof = proto(_helperlib.c_helpers['set_iter_sizeof'])
    state_type = ir.ArrayType(ir.IntType(8), setiter_sizeof())

    pstate = cgutils.alloca_once(builder, state_type, zfill=True)
    it.state = _as_bytes(builder, pstate)

    setp = _container_get_data(context, builder, iterablety.parent, it.parent)
    builder.call(fn, [it.state, setp])
    return impl_ret_borrowed(
        context,
        builder,
        sig.return_type,
        it._getvalue(),
    )


@lower_builtin('getiter', types.SetType)
def impl_set_getiter(context, builder, sig, args):
    [setty] = sig.args
    [s] = args
    iterablety = types.SetIterableType(setty)
    it = context.make_helper(builder, iterablety.iterator_type)

    fnty = ir.FunctionType(
        ir.VoidType(),
        [ll_setiter_type, ll_set_type],
    )

    fn = cgutils.get_or_insert_function(builder.module, fnty, 'numba_set_iter')

    proto = ctypes.CFUNCTYPE(ctypes.c_size_t)
    setiter_sizeof = proto(_helperlib.c_helpers['set_iter_sizeof'])
    state_type = ir.ArrayType(ir.IntType(8), setiter_sizeof())

    pstate = cgutils.alloca_once(builder, state_type, zfill=True)
    it.state = _as_bytes(builder, pstate)
    it.parent = s

    setp = _container_get_data(context, builder, iterablety.parent, args[0])
    builder.call(fn, [it.state, setp])
    return impl_ret_borrowed(
        context,
        builder,
        sig.return_type,
        it._getvalue(),
    )


@lower_builtin('iternext', types.SetIteratorType)
@iternext_impl(RefType.BORROWED)
def impl_iterator_iternext(context, builder, sig, args, result):
    iter_type = sig.args[0]
    it = context.make_helper(builder, iter_type, args[0])

    p2p_bytes = ll_bytes.as_pointer()

    iternext_fnty = ir.FunctionType(
        ll_status,
        [ll_bytes, p2p_bytes]
    )
    iternext = cgutils.get_or_insert_function(
        builder.module, iternext_fnty, 'numba_set_iter_next',
    )
    key_raw_ptr = cgutils.alloca_once(builder, ll_bytes)

    status = builder.call(iternext, (it.state, key_raw_ptr))

    is_valid = builder.icmp_unsigned('==', status, status.type(0))
    result.set_valid(is_valid)

    with builder.if_then(is_valid):
        keyty = iter_type.parent.key_type

        key_datamodel = context.data_model_manager[keyty]

        key_ptr = builder.bitcast(
            builder.load(key_raw_ptr),
            key_datamodel.get_data_type().as_pointer(),
        )
        key = key_datamodel.load_from_data_pointer(builder, key_ptr)

        # All set iterators use this common implementation.
        # Their differences are resolved here.
        if isinstance(iter_type.iterable, SetIterableType):
            result.yield_(key)
        else:
            # unreachable
            raise AssertionError('unknown type: {}'.format(iter_type.iterable))
