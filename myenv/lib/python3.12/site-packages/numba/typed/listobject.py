"""
Compiler-side implementation of the Numba  typed-list.
"""
import operator
from enum import IntEnum

from llvmlite import ir

from numba.core.extending import (
    overload,
    overload_method,
    overload_attribute,
    register_jitable,
    intrinsic,
    register_model,
    models,
    lower_builtin,
)
from numba.core.imputils import iternext_impl
from numba.core import types, cgutils
from numba.core.types import (
    ListType,
    ListTypeIterableType,
    ListTypeIteratorType,
    Type,
    NoneType,
)
from numba.core.imputils import impl_ret_borrowed, RefType
from numba.core.errors import TypingError
from numba.core import typing
from numba.typed.typedobjectutils import (_as_bytes, _cast, _nonoptional,
                                          _get_incref_decref,
                                          _container_get_data,
                                          _container_get_meminfo,)
from numba.cpython import listobj

ll_list_type = cgutils.voidptr_t
ll_listiter_type = cgutils.voidptr_t
ll_voidptr_type = cgutils.voidptr_t
ll_status = cgutils.int32_t
ll_ssize_t = cgutils.intp_t
ll_bytes = cgutils.voidptr_t


_meminfo_listptr = types.MemInfoPointer(types.voidptr)

INDEXTY = types.intp

index_types = types.integer_domain

DEFAULT_ALLOCATED = 0


@register_model(ListType)
class ListModel(models.StructModel):
    def __init__(self, dmm, fe_type):
        members = [
            ('meminfo', _meminfo_listptr),
            ('data', types.voidptr),   # ptr to the C list
        ]
        super(ListModel, self).__init__(dmm, fe_type, members)


@register_model(ListTypeIterableType)
@register_model(ListTypeIteratorType)
class ListIterModel(models.StructModel):
    def __init__(self, dmm, fe_type):
        members = [
            ('size', types.intp), # the size of the iteration space
            ('parent', fe_type.parent), # the parent list
            ('index', types.EphemeralPointer(types.intp)), # current index
        ]
        super(ListIterModel, self).__init__(dmm, fe_type, members)


class ListStatus(IntEnum):
    """Status code for other list operations.
    """
    LIST_OK = 0,
    LIST_ERR_INDEX = -1
    LIST_ERR_NO_MEMORY = -2
    LIST_ERR_MUTATED = -3
    LIST_ERR_ITER_EXHAUSTED = -4
    LIST_ERR_IMMUTABLE = -5


class ErrorHandler(object):
    """ErrorHandler for calling codegen functions from this file.

    Stores the state needed to raise an exception from nopython mode.
    """

    def __init__(self, context):
        self.context = context

    def __call__(self, builder, status, msg):
        ok_status = status.type(int(ListStatus.LIST_OK))
        with builder.if_then(builder.icmp_signed('!=', status, ok_status),
                             likely=True):
            self.context.call_conv.return_user_exc(
                builder, RuntimeError, (msg,))


def _check_for_none_typed(lst, method):
    if isinstance(lst.dtype, NoneType):
        raise TypingError("method support for List[None] is limited, "
                          "not supported: '{}'.".format(method))


@intrinsic
def _as_meminfo(typingctx, lstobj):
    """Returns the MemInfoPointer of a list.
    """
    if not isinstance(lstobj, types.ListType):
        raise TypingError('expected *lstobj* to be a ListType')

    def codegen(context, builder, sig, args):
        [tl] = sig.args
        [l] = args
        # Incref
        context.nrt.incref(builder, tl, l)
        ctor = cgutils.create_struct_proxy(tl)
        lstruct = ctor(context, builder, value=l)
        # Returns the plain MemInfo
        return lstruct.meminfo

    sig = _meminfo_listptr(lstobj)
    return sig, codegen


@intrinsic
def _from_meminfo(typingctx, mi, listtyperef):
    """Recreate a list from a MemInfoPointer
    """
    if mi != _meminfo_listptr:
        raise TypingError('expected a MemInfoPointer for list.')
    listtype = listtyperef.instance_type
    if not isinstance(listtype, ListType):
        raise TypingError('expected a {}'.format(ListType))

    def codegen(context, builder, sig, args):
        [tmi, tdref] = sig.args
        td = tdref.instance_type
        [mi, _] = args

        ctor = cgutils.create_struct_proxy(td)
        dstruct = ctor(context, builder)

        data_pointer = context.nrt.meminfo_data(builder, mi)
        data_pointer = builder.bitcast(data_pointer, ll_list_type.as_pointer())

        dstruct.data = builder.load(data_pointer)
        dstruct.meminfo = mi

        return impl_ret_borrowed(
            context,
            builder,
            listtype,
            dstruct._getvalue(),
        )

    sig = listtype(mi, listtyperef)
    return sig, codegen


def _list_codegen_set_method_table(context, builder, lp, itemty):
    vtablety = ir.LiteralStructType([
        ll_voidptr_type,  # item incref
        ll_voidptr_type,  # item decref
    ])
    setmethod_fnty = ir.FunctionType(
        ir.VoidType(),
        [ll_list_type, vtablety.as_pointer()]
    )

    setmethod_fn = cgutils.get_or_insert_function(
        builder.module,
        setmethod_fnty,
        'numba_list_set_method_table')
    vtable = cgutils.alloca_once(builder, vtablety, zfill=True)

    # install item incref/decref
    item_incref_ptr = cgutils.gep_inbounds(builder, vtable, 0, 0)
    item_decref_ptr = cgutils.gep_inbounds(builder, vtable, 0, 1)

    dm_item = context.data_model_manager[itemty]
    if dm_item.contains_nrt_meminfo():
        item_incref, item_decref = _get_incref_decref(
            context, builder.module, dm_item, "list"
        )
        builder.store(
            builder.bitcast(item_incref, item_incref_ptr.type.pointee),
            item_incref_ptr,
        )
        builder.store(
            builder.bitcast(item_decref, item_decref_ptr.type.pointee),
            item_decref_ptr,
        )

    builder.call(setmethod_fn, [lp, vtable])


@intrinsic
def _list_set_method_table(typingctx, lp, itemty):
    """Wrap numba_list_set_method_table
    """
    resty = types.void
    sig = resty(lp, itemty)

    def codegen(context, builder, sig, args):
        _list_codegen_set_method_table(
            context, builder, args[0], itemty.instance_type)

    return sig, codegen


@lower_builtin(operator.is_, types.ListType, types.ListType)
def list_is(context, builder, sig, args):
    a_meminfo = _container_get_meminfo(context, builder, sig.args[0], args[0])
    b_meminfo = _container_get_meminfo(context, builder, sig.args[1], args[1])
    ma = builder.ptrtoint(a_meminfo, cgutils.intp_t)
    mb = builder.ptrtoint(b_meminfo, cgutils.intp_t)
    return builder.icmp_signed('==', ma, mb)


def _call_list_free(context, builder, ptr):
    """Call numba_list_free(ptr)
    """
    fnty = ir.FunctionType(
        ir.VoidType(),
        [ll_list_type],
    )
    free = cgutils.get_or_insert_function(builder.module, fnty,
                                          'numba_list_free')
    builder.call(free, [ptr])


# FIXME: this needs a careful review
def _imp_dtor(context, module):
    """Define the dtor for list
    """
    llvoidptr = context.get_value_type(types.voidptr)
    llsize = context.get_value_type(types.uintp)
    fnty = ir.FunctionType(
        ir.VoidType(),
        [llvoidptr, llsize, llvoidptr],
    )
    fname = '_numba_list_dtor'
    fn = cgutils.get_or_insert_function(module, fnty, fname)

    if fn.is_declaration:
        # Set linkage
        fn.linkage = 'linkonce_odr'
        # Define
        builder = ir.IRBuilder(fn.append_basic_block())
        lp = builder.bitcast(fn.args[0], ll_list_type.as_pointer())
        l = builder.load(lp)
        _call_list_free(context, builder, l)
        builder.ret_void()

    return fn


def new_list(item, allocated=DEFAULT_ALLOCATED):
    """Construct a new list. (Not implemented in the interpreter yet)

    Parameters
    ----------
    item: TypeRef
        Item type of the new list.
    allocated: int
        number of items to pre-allocate

    """
    # With JIT disabled, ignore all arguments and return a Python list.
    return list()


def _add_meminfo(context, builder, lstruct):
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
    data_pointer = builder.bitcast(data_pointer, ll_list_type.as_pointer())
    builder.store(lstruct.data, data_pointer)
    lstruct.meminfo = meminfo


@intrinsic
def _make_list(typingctx, itemty, ptr):
    """Make a list struct with the given *ptr*

    Parameters
    ----------
    itemty: Type
        Type of the item.
    ptr : llvm pointer value
        Points to the list object.
    """
    list_ty = types.ListType(itemty.instance_type)

    def codegen(context, builder, signature, args):
        ptr = args[1]
        ctor = cgutils.create_struct_proxy(list_ty)
        lstruct = ctor(context, builder)
        lstruct.data = ptr
        _add_meminfo(context, builder, lstruct)
        return lstruct._getvalue()

    sig = list_ty(itemty, ptr)
    return sig, codegen


def _list_new_codegen(context, builder, itemty, new_size, error_handler):
    fnty = ir.FunctionType(
        ll_status,
        [ll_list_type.as_pointer(), ll_ssize_t, ll_ssize_t],
    )
    fn = cgutils.get_or_insert_function(builder.module, fnty, 'numba_list_new')
    # Determine sizeof item types
    ll_item = context.get_data_type(itemty)
    sz_item = context.get_abi_sizeof(ll_item)
    reflp = cgutils.alloca_once(builder, ll_list_type, zfill=True)
    status = builder.call(
        fn,
        [reflp, ll_ssize_t(sz_item), new_size],
    )
    msg = "Failed to allocate list"
    error_handler(
        builder,
        status,
        msg,
    )
    lp = builder.load(reflp)
    return lp


@intrinsic
def _list_new(typingctx, itemty, allocated):
    """Wrap numba_list_new.

    Allocate a new list object with zero capacity.

    Parameters
    ----------
    itemty: Type
        Type of the items
    allocated: int
        number of items to pre-allocate

    """
    resty = types.voidptr
    sig = resty(itemty, allocated)

    def codegen(context, builder, sig, args):
        error_handler = ErrorHandler(context)
        return _list_new_codegen(context,
                                 builder,
                                 itemty.instance_type,
                                 args[1],
                                 error_handler,
                                 )

    return sig, codegen


@overload(new_list)
def impl_new_list(item, allocated=DEFAULT_ALLOCATED):
    """Creates a new list.

    Parameters
    ----------
    item: Numba type
        type of the list item.
    allocated: int
        number of items to pre-allocate

    """
    if not isinstance(item, Type):
        raise TypeError("expecting *item* to be a numba Type")

    itemty = item

    def imp(item, allocated=DEFAULT_ALLOCATED):
        if allocated < 0:
            raise RuntimeError("expecting *allocated* to be >= 0")
        lp = _list_new(itemty, allocated)
        _list_set_method_table(lp, itemty)
        l = _make_list(itemty, lp)
        return l

    return imp


@overload(len)
def impl_len(l):
    """len(list)
    """
    if isinstance(l, types.ListType):
        def impl(l):
            return _list_length(l)

        return impl


@intrinsic
def _list_length(typingctx, l):
    """Wrap numba_list_length

    Returns the length of the list.
    """
    sig = types.intp(l)

    def codegen(context, builder, sig, args):
        [tl] = sig.args
        [l] = args
        fnty = ir.FunctionType(
            ll_ssize_t,
            [ll_list_type],
        )
        fname = 'numba_list_size_address'
        fn = cgutils.get_or_insert_function(builder.module, fnty, fname)
        fn.attributes.add('alwaysinline')
        fn.attributes.add('readonly')
        fn.attributes.add('nounwind')
        lp = _container_get_data(context, builder, tl, l)
        len_addr = builder.call(fn, [lp,],)
        ptr = builder.inttoptr(len_addr, cgutils.intp_t.as_pointer())
        return builder.load(ptr)

    return sig, codegen


@overload_method(types.ListType, "_allocated")
def impl_allocated(l):
    """list._allocated()
    """
    if isinstance(l, types.ListType):
        def impl(l):
            return _list_allocated(l)

        return impl


@intrinsic
def _list_allocated(typingctx, l):
    """Wrap numba_list_allocated

    Returns the allocation of the list.
    """
    resty = types.intp
    sig = resty(l)

    def codegen(context, builder, sig, args):
        fnty = ir.FunctionType(
            ll_ssize_t,
            [ll_list_type],
        )
        fn = cgutils.get_or_insert_function(builder.module, fnty,
                                            'numba_list_allocated')
        [l] = args
        [tl] = sig.args
        lp = _container_get_data(context, builder, tl, l)
        n = builder.call(fn, [lp])
        return n

    return sig, codegen


@overload_method(types.ListType, "_is_mutable")
def impl_is_mutable(l):
    """list._is_mutable()"""
    if isinstance(l, types.ListType):
        def impl(l):
            return bool(_list_is_mutable(l))

        return impl


@intrinsic
def _list_is_mutable(typingctx, l):
    """Wrap numba_list_is_mutable

    Returns the state of the is_mutable member
    """
    resty = types.int32
    sig = resty(l)

    def codegen(context, builder, sig, args):
        fnty = ir.FunctionType(
            ll_status,
            [ll_list_type],
        )
        fn = cgutils.get_or_insert_function(builder.module, fnty,
                                            'numba_list_is_mutable')
        [l] = args
        [tl] = sig.args
        lp = _container_get_data(context, builder, tl, l)
        n = builder.call(fn, [lp])
        return n

    return sig, codegen


@overload_method(types.ListType, "_make_mutable")
def impl_make_mutable(l):
    """list._make_mutable()"""
    if isinstance(l, types.ListType):
        def impl(l):
            _list_set_is_mutable(l, 1)

        return impl


@overload_method(types.ListType, "_make_immutable")
def impl_make_immutable(l):
    """list._make_immutable()"""
    if isinstance(l, types.ListType):
        def impl(l):
            _list_set_is_mutable(l, 0)

        return impl


@intrinsic
def _list_set_is_mutable(typingctx, l, is_mutable):
    """Wrap numba_list_set_mutable

    Sets the state of the is_mutable member.
    """
    resty = types.void
    sig = resty(l, is_mutable)

    def codegen(context, builder, sig, args):
        fnty = ir.FunctionType(
            ir.VoidType(),
            [ll_list_type, cgutils.intp_t],
        )
        fn = cgutils.get_or_insert_function(builder.module, fnty,
                                            'numba_list_set_is_mutable')
        [l, i] = args
        [tl, ti] = sig.args
        lp = _container_get_data(context, builder, tl, l)
        builder.call(fn, [lp, i])

    return sig, codegen


@intrinsic
def _list_append(typingctx, l, item):
    """Wrap numba_list_append
    """
    resty = types.int32
    sig = resty(l, l.item_type)

    def codegen(context, builder, sig, args):
        fnty = ir.FunctionType(
            ll_status,
            [ll_list_type, ll_bytes],
        )
        [l, item] = args
        [tl, titem] = sig.args
        fn = cgutils.get_or_insert_function(builder.module, fnty,
                                            'numba_list_append')

        dm_item = context.data_model_manager[titem]

        data_item = dm_item.as_data(builder, item)

        ptr_item = cgutils.alloca_once_value(builder, data_item)

        lp = _container_get_data(context, builder, tl, l)
        status = builder.call(
            fn,
            [
                lp,
                _as_bytes(builder, ptr_item),
            ],
        )
        return status

    return sig, codegen


@overload_method(types.ListType, 'append')
def impl_append(l, item):
    if not isinstance(l, types.ListType):
        return

    itemty = l.item_type

    def impl(l, item):
        casteditem = _cast(item, itemty)
        status = _list_append(l, casteditem)
        if status == ListStatus.LIST_OK:
            return
        elif status == ListStatus.LIST_ERR_IMMUTABLE:
            raise ValueError('list is immutable')
        elif status == ListStatus.LIST_ERR_NO_MEMORY:
            raise MemoryError('Unable to allocate memory to append item')
        else:
            raise RuntimeError('list.append failed unexpectedly')

    if l.is_precise():
        # Handle the precise case.
        return impl
    else:
        # Handle the imprecise case.
        l = l.refine(item)
        # Re-bind the item type to match the arguments.
        itemty = l.item_type
        # Create the signature that we wanted this impl to have.
        sig = typing.signature(types.void, l, itemty)
        return sig, impl


@intrinsic
def fix_index(tyctx, list_ty, index_ty):
    sig = types.intp(list_ty, index_ty)

    def codegen(context, builder, sig, args):
        [list_ty, index_ty] = sig.args
        [ll_list, ll_idx] = args
        is_negative = builder.icmp_signed('<', ll_idx,
                                          ir.Constant(ll_idx.type, 0))
        fast_len_sig, length_fn = _list_length._defn(context.typing_context,
                                                     list_ty)
        length = length_fn(context, builder, fast_len_sig, (ll_list,))
        # length is an intp
        # index can be any sort of int
        # indexing in general is done with a ssize_t which correlates to an
        # intp. In llvmlite sext and trunc are guarded to return the value
        # itself if the types are the same, so there's no need to handle the
        # "equal widths" case separately. This sexts/truncs the index to the
        # length type such that `add` works for the wraparound case.
        st = 'sext' if ll_idx.type.width < length.type.width else 'trunc'
        op = getattr(builder, st)
        fixedup_idx = op(ll_idx, length.type)
        wrapped_index = builder.add(fixedup_idx, length)
        return builder.select(is_negative, wrapped_index, fixedup_idx)
    return sig, codegen


@register_jitable
def handle_index(l, index):
    """Handle index.

    If the index is negative, convert it. If the index is out of range, raise
    an IndexError.
    """
    # convert negative indices to positive ones
    index = fix_index(l, index)
    # check that the index is in range
    if index < 0 or index >= len(l):
        raise IndexError("list index out of range")
    return index


@register_jitable
def handle_slice(l, s):
    """Handle slice.

    Convert a slice object for a given list into a range object that can be
    used to index the list. Many subtle caveats here, especially if the step is
    negative.
    """
    if len(l) == 0:  # ignore slice for empty list
        return range(0)
    ll, sa, so, se = len(l), s.start, s.stop, s.step
    if se > 0:
        start = max(ll + sa,  0) if s.start < 0 else min(ll, sa)
        stop = max(ll + so, 0) if so < 0 else min(ll, so)
    elif se < 0:
        start = max(ll + sa,  -1) if s.start < 0 else min(ll - 1, sa)
        stop = max(ll + so, -1) if so < 0 else min(ll, so)
    else:
        # should be caught earlier, but isn't, so we raise here
        raise ValueError("slice step cannot be zero")
    return range(start, stop, s.step)


def _gen_getitem(borrowed):

    @intrinsic
    def impl(typingctx, l_ty, index_ty):

        is_none = isinstance(l_ty.item_type, types.NoneType)
        if is_none:
            resty = types.Tuple([types.int32, l_ty.item_type])
        else:
            resty = types.Tuple([types.int32, types.Optional(l_ty.item_type)])
        sig = resty(l_ty, index_ty)

        def codegen(context, builder, sig, args):
            [tl, tindex] = sig.args
            [l, index] = args
            fnty = ir.FunctionType(
                ll_voidptr_type,
                [ll_list_type],
            )
            fname = 'numba_list_base_ptr'
            fn = cgutils.get_or_insert_function(builder.module, fnty, fname)
            fn.attributes.add('alwaysinline')
            fn.attributes.add('nounwind')
            fn.attributes.add('readonly')

            lp = _container_get_data(context, builder, tl, l)

            base_ptr = builder.call(
                fn,
                [lp,],
            )

            llty = context.get_data_type(tl.item_type)
            casted_base_ptr = builder.bitcast(base_ptr, llty.as_pointer())

            item_ptr = cgutils.gep(builder, casted_base_ptr, index)

            if is_none:
                out = builder.load(item_ptr)
            else:
                out = context.make_optional_none(builder, tl.item_type)
                pout = cgutils.alloca_once_value(builder, out)

                dm_item = context.data_model_manager[tl.item_type]
                item = dm_item.load_from_data_pointer(builder, item_ptr)
                if not borrowed:
                    context.nrt.incref(builder, tl.item_type, item)

                if is_none:
                    loaded = item
                else:
                    loaded = context.make_optional_value(builder, tl.item_type,
                                                         item)
                builder.store(loaded, pout)

                out = builder.load(pout)
            return context.make_tuple(builder, resty, [ll_status(0), out])

        return sig, codegen
    return impl


_list_getitem = _gen_getitem(False)
_list_getitem_borrowed = _gen_getitem(True)


@overload(operator.getitem)
def impl_getitem(l, index):
    if not isinstance(l, types.ListType):
        return

    indexty = INDEXTY
    itemty = l.item_type
    IS_NOT_NONE = not isinstance(l.item_type, types.NoneType)

    if index in index_types:
        if IS_NOT_NONE:
            def integer_non_none_impl(l, index):
                castedindex = _cast(index, indexty)
                handledindex = handle_index(l, castedindex)
                status, item = _list_getitem(l, handledindex)
                if status == ListStatus.LIST_OK:
                    return _nonoptional(item)
                else:
                    raise AssertionError("internal list error during getitem")
            return integer_non_none_impl
        else:
            def integer_none_impl(l, index):
                index = handle_index(l, index)
                return None
            return integer_none_impl

    elif isinstance(index, types.SliceType):
        def slice_impl(l, index):
            newl = new_list(itemty)
            for i in handle_slice(l, index):
                newl.append(l[i])
            return newl

        return slice_impl

    else:
        raise TypingError("list indices must be integers or slices")


@intrinsic
def _list_setitem(typingctx, l, index, item):
    """Wrap numba_list_setitem
    """
    resty = types.int32
    sig = resty(l, index, item)

    def codegen(context, builder, sig, args):
        fnty = ir.FunctionType(
            ll_status,
            [ll_list_type, ll_ssize_t, ll_bytes],
        )
        [l, index, item] = args
        [tl, tindex, titem] = sig.args
        fn = cgutils.get_or_insert_function(builder.module, fnty,
                                            'numba_list_setitem')

        dm_item = context.data_model_manager[titem]
        data_item = dm_item.as_data(builder, item)
        ptr_item = cgutils.alloca_once_value(builder, data_item)

        lp = _container_get_data(context, builder, tl, l)
        status = builder.call(
            fn,
            [
                lp,
                index,
                _as_bytes(builder, ptr_item),
            ],
        )
        return status

    return sig, codegen


@overload(operator.setitem)
def impl_setitem(l, index, item):
    if not isinstance(l, types.ListType):
        return

    indexty = INDEXTY
    itemty = l.item_type

    if index in index_types:
        def impl_integer(l, index, item):
            index = handle_index(l, index)
            castedindex = _cast(index, indexty)
            casteditem = _cast(item, itemty)
            status = _list_setitem(l, castedindex, casteditem)
            if status == ListStatus.LIST_OK:
                return
            elif status == ListStatus.LIST_ERR_IMMUTABLE:
                raise ValueError("list is immutable")
            else:
                raise AssertionError("internal list error during settitem")

        return impl_integer

    elif isinstance(index, types.SliceType):
        if not isinstance(item, types.IterableType):
            raise TypingError("can only assign an iterable when using a slice "
                              "with assignment/setitem")

        def impl_slice(l, index, item):
            if not l._is_mutable():
                raise ValueError("list is immutable")
            # special case "a[i:j] = a", need to copy first
            if l is item:
                item = item.copy()
            slice_range = handle_slice(l, index)
            # non-extended (simple) slices
            if slice_range.step == 1:
                # replace
                if len(item) == len(slice_range):
                    for i, j in zip(slice_range, item):
                        l[i] = j
                # replace and insert
                if len(item) > len(slice_range):
                    # do the replaces we can
                    for i, j in zip(slice_range, item[:len(slice_range)]):
                        l[i] = j
                    # insert the remaining ones
                    insert_range = range(slice_range.stop,
                                         slice_range.stop +
                                         len(item) - len(slice_range))
                    for i, k in zip(insert_range, item[len(slice_range):]):
                        # FIXME: This may be slow.  Each insert can incur a
                        # memory copy of one or more items.
                        l.insert(i, k)
                # replace and delete
                if len(item) < len(slice_range):
                    # do the replaces we can
                    replace_range = range(slice_range.start,
                                          slice_range.start + len(item))
                    for i,j in zip(replace_range, item):
                        l[i] = j
                    # delete remaining ones
                    del l[slice_range.start + len(item):slice_range.stop]
            # Extended slices
            else:
                if len(slice_range) != len(item):
                    raise ValueError("length mismatch for extended slice "
                                     "and sequence")
                # extended slice can only replace
                for i, j in zip(slice_range, item):
                    l[i] = j

        return impl_slice

    else:
        raise TypingError("list indices must be integers or slices")


@overload_method(types.ListType, 'pop')
def impl_pop(l, index=-1):
    if not isinstance(l, types.ListType):
        return

    _check_for_none_typed(l, 'pop')

    indexty = INDEXTY

    # FIXME: this type check works, but it isn't clear why and if it optimal
    if (isinstance(index, int)
            or index in index_types
            or isinstance(index, types.Omitted)):
        def impl(l, index=-1):
            if len(l) == 0:
                raise IndexError("pop from empty list")
            cindex = _cast(handle_index(l, index), indexty)
            item = l[cindex]
            del l[cindex]
            return item
        return impl

    else:
        raise TypingError("argument for pop must be an integer")


@intrinsic
def _list_delitem(typingctx, l, index):
    resty = types.int32
    sig = resty(l, index)

    def codegen(context, builder, sig, args):
        fnty = ir.FunctionType(
            ll_status,
            [ll_list_type, ll_ssize_t],
        )
        [tl, tindex] = sig.args
        [l, index] = args
        fn = cgutils.get_or_insert_function(builder.module, fnty,
                                            'numba_list_delitem')

        lp = _container_get_data(context, builder, tl, l)
        status = builder.call(fn, [lp, index])
        return status

    return sig, codegen


@intrinsic
def _list_delete_slice(typingctx, l, start, stop, step):
    """Wrap numba_list_delete_slice
    """
    resty = types.int32
    sig = resty(l, start, stop, step)

    def codegen(context, builder, sig, args):
        fnty = ir.FunctionType(
            ll_status,
            [ll_list_type, ll_ssize_t, ll_ssize_t, ll_ssize_t],
        )
        [l, start, stop, step] = args
        [tl, tstart, tstop, tstep] = sig.args
        fn = cgutils.get_or_insert_function(builder.module, fnty,
                                            'numba_list_delete_slice')

        lp = _container_get_data(context, builder, tl, l)
        status = builder.call(
            fn,
            [
                lp,
                start,
                stop,
                step,
            ],
        )
        return status

    return sig, codegen


@overload(operator.delitem)
def impl_delitem(l, index):
    if not isinstance(l, types.ListType):
        return

    _check_for_none_typed(l, 'delitem')

    if index in index_types:
        def integer_impl(l, index):
            cindex = _cast(handle_index(l, index), INDEXTY)
            status = _list_delitem(l, cindex)
            if status == ListStatus.LIST_OK:
                return
            elif status == ListStatus.LIST_ERR_IMMUTABLE:
                raise ValueError("list is immutable")
            else:
                raise AssertionError("internal list error during delitem")
        return integer_impl

    elif isinstance(index, types.SliceType):
        def slice_impl(l, index):
            slice_range = handle_slice(l, index)
            status = _list_delete_slice(
                l,
                slice_range.start,
                slice_range.stop,
                slice_range.step)
            if status == ListStatus.LIST_ERR_MUTATED:
                raise ValueError("list is immutable")
        return slice_impl

    else:
        raise TypingError("list indices must be integers or slices")


@overload(operator.contains)
def impl_contains(l, item):
    if not isinstance(l, types.ListType):
        return

    itemty = l.item_type
    _check_for_none_typed(l, "__contains__")

    def impl(l, item):
        casteditem = _cast(item, itemty)
        for i in l:
            if i == casteditem:
                return True
        else:
            return False
    return impl


@overload_method(types.ListType, 'count')
def impl_count(l, item):
    if not isinstance(l, types.ListType):
        return

    _check_for_none_typed(l, 'count')

    itemty = l.item_type

    def impl(l, item):
        casteditem = _cast(item, itemty)
        total = 0
        for i in l:
            if i == casteditem:
                total += 1
        return total

    return impl


@overload_method(types.ListType, 'extend')
def impl_extend(l, iterable):
    if not isinstance(l, types.ListType):
        return
    if not isinstance(iterable, types.IterableType):
        raise TypingError("extend argument must be iterable")

    _check_for_none_typed(l, 'extend')

    def select_impl():
        if isinstance(iterable, types.ListType):
            def impl(l, iterable):
                if not l._is_mutable():
                    raise ValueError("list is immutable")
                # guard against l.extend(l)
                if l is iterable:
                    iterable = iterable.copy()
                for i in iterable:
                    l.append(i)

            return impl
        else:
            def impl(l, iterable):
                for i in iterable:
                    l.append(i)

            return impl

    if l.is_precise():
        # Handle the precise case.
        return select_impl()
    else:
        # Handle the imprecise case, try to 'guess' the underlying type of the
        # values in the iterable.
        if hasattr(iterable, "dtype"):  # tuples and arrays
            ty = iterable.dtype
        elif hasattr(iterable, "item_type"):  # lists
            ty = iterable.item_type
        elif hasattr(iterable, "yield_type"):  # iterators and generators
            ty = iterable.yield_type
        elif isinstance(iterable, types.UnicodeType):
            ty = iterable
        else:
            raise TypingError("unable to extend list, iterable is missing "
                              "either *dtype*, *item_type* or *yield_type*.")
        l = l.refine(ty)
        # Create the signature that we wanted this impl to have
        sig = typing.signature(types.void, l, iterable)
        return sig, select_impl()


@overload_method(types.ListType, 'insert')
def impl_insert(l, index, item):
    if not isinstance(l, types.ListType):
        return

    _check_for_none_typed(l, 'insert')
    # insert can refine
    if isinstance(item, NoneType):
        raise TypingError("method support for List[None] is limited")

    if index in index_types:
        def impl(l, index, item):
            # If the index is larger than the size of the list or if the list is
            # empty, just append.
            if index >= len(l) or len(l) == 0:
                l.append(item)
            # Else, do the insert dance
            else:
                # convert negative indices
                if index < 0:
                    # if the index is still negative after conversion, use 0
                    index = max(len(l) + index, 0)
                # grow the list by one, make room for item to insert
                l.append(l[0])
                # reverse iterate over the list and shift all elements
                i = len(l) - 1
                while (i > index):
                    l[i] = l[i - 1]
                    i -= 1
                # finally, insert the item
                l[index] = item

        if l.is_precise():
            # Handle the precise case.
            return impl
        else:
            # Handle the imprecise case
            l = l.refine(item)
            # Re-bind the item type to match the arguments.
            itemty = l.item_type
            # Create the signature that we wanted this impl to have.
            sig = typing.signature(types.void, l, INDEXTY, itemty)
            return sig, impl
    else:
        raise TypingError("list insert indices must be integers")


@overload_method(types.ListType, 'remove')
def impl_remove(l, item):
    if not isinstance(l, types.ListType):
        return

    _check_for_none_typed(l, 'remove')

    itemty = l.item_type

    def impl(l, item):
        casteditem = _cast(item, itemty)
        for i, n in enumerate(l):
            if casteditem == n:
                del l[i]
                return
        else:
            raise ValueError("list.remove(x): x not in list")

    return impl


@overload_method(types.ListType, 'clear')
def impl_clear(l):
    if not isinstance(l, types.ListType):
        return

    def impl(l):
        while len(l):
            del l[-1]

    return impl


@overload_method(types.ListType, 'reverse')
def impl_reverse(l):
    if not isinstance(l, types.ListType):
        return

    _check_for_none_typed(l, 'reverse')

    def impl(l):
        if not l._is_mutable():
            raise ValueError("list is immutable")
        front = 0
        back = len(l) - 1
        while front < back:
            l[front], l[back] = l[back], l[front]
            front += 1
            back -= 1

    return impl


@overload_method(types.ListType, 'copy')
def impl_copy(l):

    _check_for_none_typed(l, 'copy')

    itemty = l.item_type

    if isinstance(l, types.ListType):
        def impl(l):
            newl = new_list(itemty, len(l))
            for i in l:
                newl.append(i)
            return newl

        return impl


@overload_method(types.ListType, 'index')
def impl_index(l, item, start=None, end=None):
    if not isinstance(l, types.ListType):
        return

    _check_for_none_typed(l, 'index')

    itemty = l.item_type

    def check_arg(arg, name):
        if not (arg is None
                or arg in index_types
                or isinstance(arg, (types.Omitted, types.NoneType))):
            raise TypingError("{} argument for index must be an integer"
                              .format(name))
    check_arg(start, "start")
    check_arg(end, "end")

    def impl(l, item, start=None, end=None):
        casteditem = _cast(item, itemty)
        for i in handle_slice(l, slice(start, end, 1)):
            if l[i] == casteditem:
                return i
        else:
            raise ValueError("item not in list")

    return impl


@overload_method(types.ListType, "sort")
def ol_list_sort(lst, key=None, reverse=False):
    # The following is mostly borrowed from listobj.ol_list_sort
    from numba.typed import List

    listobj._sort_check_key(key)
    listobj._sort_check_reverse(reverse)

    if cgutils.is_nonelike(key):
        KEY = False
        sort_f = listobj.sort_forwards
        sort_b = listobj.sort_backwards
    elif isinstance(key, types.Dispatcher):
        KEY = True
        sort_f = listobj.arg_sort_forwards
        sort_b = listobj.arg_sort_backwards

    def impl(lst, key=None, reverse=False):
        if not lst._is_mutable():
            raise ValueError("list is immutable")
        if KEY is True:
            # There's an unknown refct problem in reflected list.
            # Using an explicit loop with typedlist somehow "fixed" it.
            _lst = List()
            for x in lst:
                _lst.append(key(x))
        else:
            _lst = lst
        if reverse is False or reverse == 0:
            tmp = sort_f(_lst)
        else:
            tmp = sort_b(_lst)
        if KEY is True:
            # There's an unknown refct problem in reflected list.
            # Using an explicit loop with typedlist somehow "fixed" it.
            ordered = List()
            for i in tmp:
                ordered.append(lst[i])
            lst[:] = ordered
    return impl


@overload_method(types.ListType, "getitem_unchecked")
def ol_getitem_unchecked(lst, index):
    if not isinstance(index, types.Integer):
        return

    def impl(lst, index):
        index = fix_index(lst, index)
        castedindex = _cast(index, types.intp)
        _, item = _list_getitem(lst, castedindex)
        return _nonoptional(item)
    return impl


@overload_attribute(types.ListType, '__hash__')
def ol_list_hash(lst):
    if not isinstance(lst, types.ListType):
        return
    return lambda lst: None


@overload_attribute(types.ListType, '_dtype')
def impl_dtype(l):
    if not isinstance(l, types.ListType):
        return
    dt = l.dtype

    def impl(l):
        return dt

    return impl


def _equals_helper(this, other, OP):
    if not isinstance(this, types.ListType):
        return
    if not isinstance(other, types.ListType):
        return lambda this, other: False

    this_is_none = isinstance(this.dtype, types.NoneType)
    other_is_none = isinstance(other.dtype, types.NoneType)

    if this_is_none or other_is_none:
        def impl_some_none(this, other):
            def equals(this, other):
                # Equal if both none-typed and have equal length
                return bool(this_is_none == other_is_none
                            and len(this) == len(other))
            return OP(equals(this, other))
        return impl_some_none
    else:
        def impl_not_none(this, other):
            def equals(this, other):
                if len(this) != len(other):
                    return False
                for i in range(len(this)):
                    if this[i] != other[i]:
                        return False
                else:
                    return True
            return OP(equals(this, other))
        return impl_not_none


@overload(operator.eq)
def impl_equals(this, other):
    return _equals_helper(this, other, operator.truth)


@overload(operator.ne)
def impl_not_equals(this, other):
    return _equals_helper(this, other, operator.not_)


@register_jitable
def compare_not_none(this, other):
    """Oldschool (python 2.x) cmp.

       if this < other return -1
       if this = other return 0
       if this > other return 1
    """
    if len(this) != len(other):
        return -1 if len(this) < len(other) else 1
    for i in range(len(this)):
        this_item, other_item = this[i], other[i]
        if this_item != other_item:
            return -1 if this_item < other_item else 1
    else:
        return 0


@register_jitable
def compare_some_none(this, other, this_is_none, other_is_none):
    """Oldschool (python 2.x) cmp for None typed lists.

       if this < other return -1
       if this = other return 0
       if this > other return 1
    """
    if len(this) != len(other):
        return -1 if len(this) < len(other) else 1
    if this_is_none and other_is_none: # both none
        return 0
    # to get here there is precisely one none, and if the first is none, by
    # induction, the second cannot be
    return -1 if this_is_none else 1


def compare_helper(this, other, accepted):
    if not isinstance(this, types.ListType):
        return
    if not isinstance(other, types.ListType):
        return lambda this, other: False

    this_is_none = isinstance(this.dtype, types.NoneType)
    other_is_none = isinstance(other.dtype, types.NoneType)

    if this_is_none or other_is_none:
        def impl(this, other):
            return compare_some_none(
                this, other, this_is_none, other_is_none) in accepted
    else:
        def impl(this, other):
            return compare_not_none(this, other) in accepted
    return impl


@overload(operator.lt)
def impl_less_than(this, other):
    return compare_helper(this, other, (-1, ))


@overload(operator.le)
def impl_less_than_or_equal(this, other):
    return compare_helper(this, other, (-1, 0))


@overload(operator.gt)
def impl_greater_than(this, other):
    return compare_helper(this, other, (1,))


@overload(operator.ge)
def impl_greater_than_or_equal(this, other):
    return compare_helper(this, other, (0, 1))


class ListIterInstance(object):

    def __init__(self, context, builder, iter_type, iter_val):
        self._context = context
        self._builder = builder
        self._iter_ty = iter_type
        self._list_ty = self._iter_ty.parent
        self._iter = context.make_helper(builder, iter_type, iter_val)

    @classmethod
    def from_list(cls, context, builder, iter_type, list_val):
        self = cls(context, builder, iter_type, None)
        index = context.get_constant(types.intp, 0)
        self._iter.index = cgutils.alloca_once_value(builder, index)
        self._iter.parent = list_val
        self._iter.size = cls._size_of_list(context, builder, self._list_ty,
                                            self._iter.parent)
        return self

    @classmethod
    def _size_of_list(cls, context, builder, list_ty, ll_list):
        tyctx = context.typing_context
        fnty = tyctx.resolve_value_type(len)
        sig = fnty.get_call_type(tyctx, (list_ty,), {})
        impl = context.get_function(fnty, sig)
        return impl(builder, (ll_list,))

    @property
    def size(self):
        tyctx = self._context.typing_context
        fnty = tyctx.resolve_value_type(len)
        ty = self._list_ty
        sig = fnty.get_call_type(tyctx, (ty,), {})
        impl = self._context.get_function(fnty, sig)
        return impl(self._builder, (self._iter.parent,))

    @property
    def value(self):
        return self._iter._getvalue()

    def getitem(self, index):
        tyctx = self._context.typing_context
        ty = self._list_ty
        sig, fn = _list_getitem_borrowed._defn(tyctx, ty, types.intp)

        statnitem = fn(self._context, self._builder, sig, (self._iter.parent,
                                                           index))
        _, item = cgutils.unpack_tuple(self._builder, statnitem)
        retty = sig.return_type[1]
        if isinstance(self._list_ty.dtype, types.NoneType):
            raw_ty = self._list_ty.dtype
        else:
            raw_ty = retty.type
        raw_item = self._context.cast(self._builder, item, retty, raw_ty)
        return raw_item

    @property
    def index(self):
        return self._builder.load(self._iter.index)

    @index.setter
    def index(self, value):
        self._builder.store(value, self._iter.index)


@lower_builtin('getiter', types.ListType)
def getiter_list(context, builder, sig, args):
    inst = ListIterInstance.from_list(context, builder, sig.return_type,
                                      args[0])
    return impl_ret_borrowed(context, builder, sig.return_type, inst.value)


@lower_builtin('iternext', types.ListTypeIteratorType)
@iternext_impl(RefType.BORROWED)
def iternext_listiter(context, builder, sig, args, result):
    inst = ListIterInstance(context, builder, sig.args[0], args[0])
    index = inst.index

    nitems = inst.size # this is current size
    init_size = inst._iter.size # this is initial size

    # if the current count is different to the initial count, bail, list is
    # being mutated whilst iterated.
    is_mutated = builder.icmp_signed('!=', init_size, nitems)
    with builder.if_then(is_mutated, likely=False):
        context.call_conv.return_user_exc(
            builder, RuntimeError, ("list was mutated during iteration",))

    is_valid = builder.icmp_signed('<', index, nitems)
    result.set_valid(is_valid)
    with builder.if_then(is_valid):
        result.yield_(inst.getitem(index))
        inst.index = builder.add(index, context.get_constant(types.intp, 1))
