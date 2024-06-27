"""
Implementation of operations on Array objects and objects supporting
the buffer protocol.
"""

import functools
import math
import operator

from llvmlite import ir
from llvmlite.ir import Constant

import numpy as np

from numba import pndindex, literal_unroll
from numba.core import types, typing, errors, cgutils, extending
from numba.np.numpy_support import (as_dtype, from_dtype, carray, farray,
                                    is_contiguous, is_fortran,
                                    check_is_integer, type_is_scalar,
                                    lt_complex, lt_floats)
from numba.np.numpy_support import type_can_asarray, is_nonelike, numpy_version
from numba.core.imputils import (lower_builtin, lower_getattr,
                                 lower_getattr_generic,
                                 lower_setattr_generic,
                                 lower_cast, lower_constant,
                                 iternext_impl, impl_ret_borrowed,
                                 impl_ret_new_ref, impl_ret_untracked,
                                 RefType)
from numba.core.typing import signature
from numba.core.types import StringLiteral
from numba.core.extending import (register_jitable, overload, overload_method,
                                  intrinsic, overload_attribute)
from numba.misc import quicksort, mergesort
from numba.cpython import slicing
from numba.cpython.unsafe.tuple import tuple_setitem, build_full_slice_tuple
from numba.core.extending import overload_classmethod
from numba.core.typing.npydecl import (parse_dtype as ty_parse_dtype,
                                       parse_shape as ty_parse_shape,
                                       _parse_nested_sequence,
                                       _sequence_of_arrays,
                                       _choose_concatenation_layout)


def set_range_metadata(builder, load, lower_bound, upper_bound):
    """
    Set the "range" metadata on a load instruction.
    Note the interval is in the form [lower_bound, upper_bound).
    """
    range_operands = [Constant(load.type, lower_bound),
                      Constant(load.type, upper_bound)]
    md = builder.module.add_metadata(range_operands)
    load.set_metadata("range", md)


def mark_positive(builder, load):
    """
    Mark the result of a load instruction as positive (or zero).
    """
    upper_bound = (1 << (load.type.width - 1)) - 1
    set_range_metadata(builder, load, 0, upper_bound)


def make_array(array_type):
    """
    Return the Structure representation of the given *array_type*
    (an instance of types.ArrayCompatible).

    Note this does not call __array_wrap__ in case a new array structure
    is being created (rather than populated).
    """
    real_array_type = array_type.as_array
    base = cgutils.create_struct_proxy(real_array_type)
    ndim = real_array_type.ndim

    class ArrayStruct(base):

        def _make_refs(self, ref):
            sig = signature(real_array_type, array_type)
            try:
                array_impl = self._context.get_function('__array__', sig)
            except NotImplementedError:
                return super(ArrayStruct, self)._make_refs(ref)

            # Return a wrapped structure and its unwrapped reference
            datamodel = self._context.data_model_manager[array_type]
            be_type = self._get_be_type(datamodel)
            if ref is None:
                outer_ref = cgutils.alloca_once(self._builder, be_type,
                                                zfill=True)
            else:
                outer_ref = ref
            # NOTE: __array__ is called with a pointer and expects a pointer
            # in return!
            ref = array_impl(self._builder, (outer_ref,))
            return outer_ref, ref

        @property
        def shape(self):
            """
            Override .shape to inform LLVM that its elements are all positive.
            """
            builder = self._builder
            if ndim == 0:
                return base.__getattr__(self, "shape")

            # Unfortunately, we can't use llvm.assume as its presence can
            # seriously pessimize performance,
            # *and* the range metadata currently isn't improving anything here,
            # see https://llvm.org/bugs/show_bug.cgi?id=23848 !
            ptr = self._get_ptr_by_name("shape")
            dims = []
            for i in range(ndim):
                dimptr = cgutils.gep_inbounds(builder, ptr, 0, i)
                load = builder.load(dimptr)
                dims.append(load)
                mark_positive(builder, load)

            return cgutils.pack_array(builder, dims)

    return ArrayStruct


def get_itemsize(context, array_type):
    """
    Return the item size for the given array or buffer type.
    """
    llty = context.get_data_type(array_type.dtype)
    return context.get_abi_sizeof(llty)


def load_item(context, builder, arrayty, ptr):
    """
    Load the item at the given array pointer.
    """
    align = None if arrayty.aligned else 1
    return context.unpack_value(builder, arrayty.dtype, ptr,
                                align=align)


def store_item(context, builder, arrayty, val, ptr):
    """
    Store the item at the given array pointer.
    """
    align = None if arrayty.aligned else 1
    return context.pack_value(builder, arrayty.dtype, val, ptr, align=align)


def fix_integer_index(context, builder, idxty, idx, size):
    """
    Fix the integer index' type and value for the given dimension size.
    """
    if idxty.signed:
        ind = context.cast(builder, idx, idxty, types.intp)
        ind = slicing.fix_index(builder, ind, size)
    else:
        ind = context.cast(builder, idx, idxty, types.uintp)
    return ind


def normalize_index(context, builder, idxty, idx):
    """
    Normalize the index type and value.  0-d arrays are converted to scalars.
    """
    if isinstance(idxty, types.Array) and idxty.ndim == 0:
        assert isinstance(idxty.dtype, types.Integer)
        idxary = make_array(idxty)(context, builder, idx)
        idxval = load_item(context, builder, idxty, idxary.data)
        return idxty.dtype, idxval
    else:
        return idxty, idx


def normalize_indices(context, builder, index_types, indices):
    """
    Same as normalize_index(), but operating on sequences of
    index types and values.
    """
    if len(indices):
        index_types, indices = zip(*[normalize_index(context, builder, idxty,
                                                     idx)
                                     for idxty, idx in zip(index_types, indices)
                                     ])
    return index_types, indices


def populate_array(array, data, shape, strides, itemsize, meminfo,
                   parent=None):
    """
    Helper function for populating array structures.
    This avoids forgetting to set fields.

    *shape* and *strides* can be Python tuples or LLVM arrays.
    """
    context = array._context
    builder = array._builder
    datamodel = array._datamodel
    # doesn't matter what this array type instance is, it's just to get the
    # fields for the datamodel of the standard array type in this context
    standard_array = types.Array(types.float64, 1, 'C')
    standard_array_type_datamodel = context.data_model_manager[standard_array]
    required_fields = set(standard_array_type_datamodel._fields)
    datamodel_fields = set(datamodel._fields)
    # Make sure that the presented array object has a data model that is close
    # enough to an array for this function to proceed.
    if (required_fields & datamodel_fields) != required_fields:
        missing = required_fields - datamodel_fields
        msg = (f"The datamodel for type {array._fe_type} is missing "
               f"field{'s' if len(missing) > 1 else ''} {missing}.")
        raise ValueError(msg)

    if meminfo is None:
        meminfo = Constant(context.get_value_type(
            datamodel.get_type('meminfo')), None)

    intp_t = context.get_value_type(types.intp)
    if isinstance(shape, (tuple, list)):
        shape = cgutils.pack_array(builder, shape, intp_t)
    if isinstance(strides, (tuple, list)):
        strides = cgutils.pack_array(builder, strides, intp_t)
    if isinstance(itemsize, int):
        itemsize = intp_t(itemsize)

    attrs = dict(shape=shape,
                 strides=strides,
                 data=data,
                 itemsize=itemsize,
                 meminfo=meminfo,)

    # Set `parent` attribute
    if parent is None:
        attrs['parent'] = Constant(context.get_value_type(
            datamodel.get_type('parent')), None)
    else:
        attrs['parent'] = parent
    # Calc num of items from shape
    nitems = context.get_constant(types.intp, 1)
    unpacked_shape = cgutils.unpack_tuple(builder, shape, shape.type.count)
    # (note empty shape => 0d array therefore nitems = 1)
    for axlen in unpacked_shape:
        nitems = builder.mul(nitems, axlen, flags=['nsw'])
    attrs['nitems'] = nitems

    # Make sure that we have all the fields
    got_fields = set(attrs.keys())
    if got_fields != required_fields:
        raise ValueError("missing {0}".format(required_fields - got_fields))

    # Set field value
    for k, v in attrs.items():
        setattr(array, k, v)

    return array


def update_array_info(aryty, array):
    """
    Update some auxiliary information in *array* after some of its fields
    were changed.  `itemsize` and `nitems` are updated.
    """
    context = array._context
    builder = array._builder

    # Calc num of items from shape
    nitems = context.get_constant(types.intp, 1)
    unpacked_shape = cgutils.unpack_tuple(builder, array.shape, aryty.ndim)
    for axlen in unpacked_shape:
        nitems = builder.mul(nitems, axlen, flags=['nsw'])
    array.nitems = nitems

    array.itemsize = context.get_constant(types.intp,
                                          get_itemsize(context, aryty))


def normalize_axis(func_name, arg_name, ndim, axis):
    """Constrain axis values to valid positive values."""
    raise NotImplementedError()


@overload(normalize_axis)
def normalize_axis_overloads(func_name, arg_name, ndim, axis):
    if not isinstance(func_name, StringLiteral):
        raise errors.TypingError("func_name must be a str literal.")
    if not isinstance(arg_name, StringLiteral):
        raise errors.TypingError("arg_name must be a str literal.")

    msg = (
        f"{func_name.literal_value}: Argument {arg_name.literal_value} "
        "out of bounds for dimensions of the array"
    )

    def impl(func_name, arg_name, ndim, axis):
        if axis < 0:
            axis += ndim
        if axis < 0 or axis >= ndim:
            raise ValueError(msg)

        return axis

    return impl


@lower_builtin('getiter', types.Buffer)
def getiter_array(context, builder, sig, args):
    [arrayty] = sig.args
    [array] = args

    iterobj = context.make_helper(builder, sig.return_type)

    zero = context.get_constant(types.intp, 0)
    indexptr = cgutils.alloca_once_value(builder, zero)

    iterobj.index = indexptr
    iterobj.array = array

    # Incref array
    if context.enable_nrt:
        context.nrt.incref(builder, arrayty, array)

    res = iterobj._getvalue()

    # Note: a decref on the iterator will dereference all internal MemInfo*
    out = impl_ret_new_ref(context, builder, sig.return_type, res)
    return out


def _getitem_array_single_int(context, builder, return_type, aryty, ary, idx):
    """ Evaluate `ary[idx]`, where idx is a single int. """
    # optimized form of _getitem_array_generic
    shapes = cgutils.unpack_tuple(builder, ary.shape, count=aryty.ndim)
    strides = cgutils.unpack_tuple(builder, ary.strides, count=aryty.ndim)
    offset = builder.mul(strides[0], idx)
    dataptr = cgutils.pointer_add(builder, ary.data, offset)
    view_shapes = shapes[1:]
    view_strides = strides[1:]

    if isinstance(return_type, types.Buffer):
        # Build array view
        retary = make_view(context, builder, aryty, ary, return_type,
                           dataptr, view_shapes, view_strides)
        return retary._getvalue()
    else:
        # Load scalar from 0-d result
        assert not view_shapes
        return load_item(context, builder, aryty, dataptr)


@lower_builtin('iternext', types.ArrayIterator)
@iternext_impl(RefType.BORROWED)
def iternext_array(context, builder, sig, args, result):
    [iterty] = sig.args
    [iter] = args
    arrayty = iterty.array_type

    iterobj = context.make_helper(builder, iterty, value=iter)
    ary = make_array(arrayty)(context, builder, value=iterobj.array)

    nitems, = cgutils.unpack_tuple(builder, ary.shape, count=1)

    index = builder.load(iterobj.index)
    is_valid = builder.icmp_signed('<', index, nitems)
    result.set_valid(is_valid)

    with builder.if_then(is_valid):
        value = _getitem_array_single_int(
            context, builder, iterty.yield_type, arrayty, ary, index
        )
        result.yield_(value)
        nindex = cgutils.increment_index(builder, index)
        builder.store(nindex, iterobj.index)


# ------------------------------------------------------------------------------
# Basic indexing (with integers and slices only)

def basic_indexing(context, builder, aryty, ary, index_types, indices,
                   boundscheck=None):
    """
    Perform basic indexing on the given array.
    A (data pointer, shapes, strides) tuple is returned describing
    the corresponding view.
    """
    zero = context.get_constant(types.intp, 0)
    one = context.get_constant(types.intp, 1)

    shapes = cgutils.unpack_tuple(builder, ary.shape, aryty.ndim)
    strides = cgutils.unpack_tuple(builder, ary.strides, aryty.ndim)

    output_indices = []
    output_shapes = []
    output_strides = []

    num_newaxes = len([idx for idx in index_types if is_nonelike(idx)])
    ax = 0
    for indexval, idxty in zip(indices, index_types):
        if idxty is types.ellipsis:
            # Fill up missing dimensions at the middle
            n_missing = aryty.ndim - len(indices) + 1 + num_newaxes
            for i in range(n_missing):
                output_indices.append(zero)
                output_shapes.append(shapes[ax])
                output_strides.append(strides[ax])
                ax += 1
            continue
        # Regular index value
        if isinstance(idxty, types.SliceType):
            slice = context.make_helper(builder, idxty, value=indexval)
            slicing.guard_invalid_slice(context, builder, idxty, slice)
            slicing.fix_slice(builder, slice, shapes[ax])
            output_indices.append(slice.start)
            sh = slicing.get_slice_length(builder, slice)
            st = slicing.fix_stride(builder, slice, strides[ax])
            output_shapes.append(sh)
            output_strides.append(st)
        elif isinstance(idxty, types.Integer):
            ind = fix_integer_index(context, builder, idxty, indexval,
                                    shapes[ax])
            if boundscheck:
                cgutils.do_boundscheck(context, builder, ind, shapes[ax], ax)
            output_indices.append(ind)
        elif is_nonelike(idxty):
            output_shapes.append(one)
            output_strides.append(zero)
            ax -= 1
        else:
            raise NotImplementedError("unexpected index type: %s" % (idxty,))
        ax += 1

    # Fill up missing dimensions at the end
    assert ax <= aryty.ndim
    while ax < aryty.ndim:
        output_shapes.append(shapes[ax])
        output_strides.append(strides[ax])
        ax += 1

    # No need to check wraparound, as negative indices were already
    # fixed in the loop above.
    dataptr = cgutils.get_item_pointer(context, builder, aryty, ary,
                                       output_indices,
                                       wraparound=False, boundscheck=False)
    return (dataptr, output_shapes, output_strides)


def make_view(context, builder, aryty, ary, return_type,
              data, shapes, strides):
    """
    Build a view over the given array with the given parameters.
    """
    retary = make_array(return_type)(context, builder)
    populate_array(retary,
                   data=data,
                   shape=shapes,
                   strides=strides,
                   itemsize=ary.itemsize,
                   meminfo=ary.meminfo,
                   parent=ary.parent)
    return retary


def _getitem_array_generic(context, builder, return_type, aryty, ary,
                           index_types, indices):
    """
    Return the result of indexing *ary* with the given *indices*,
    returning either a scalar or a view.
    """
    dataptr, view_shapes, view_strides = \
        basic_indexing(context, builder, aryty, ary, index_types, indices,
                       boundscheck=context.enable_boundscheck)

    if isinstance(return_type, types.Buffer):
        # Build array view
        retary = make_view(context, builder, aryty, ary, return_type,
                           dataptr, view_shapes, view_strides)
        return retary._getvalue()
    else:
        # Load scalar from 0-d result
        assert not view_shapes
        return load_item(context, builder, aryty, dataptr)


@lower_builtin(operator.getitem, types.Buffer, types.Integer)
@lower_builtin(operator.getitem, types.Buffer, types.SliceType)
def getitem_arraynd_intp(context, builder, sig, args):
    """
    Basic indexing with an integer or a slice.
    """
    aryty, idxty = sig.args
    ary, idx = args

    assert aryty.ndim >= 1
    ary = make_array(aryty)(context, builder, ary)

    res = _getitem_array_generic(context, builder, sig.return_type,
                                 aryty, ary, (idxty,), (idx,))
    return impl_ret_borrowed(context, builder, sig.return_type, res)


@lower_builtin(operator.getitem, types.Buffer, types.BaseTuple)
def getitem_array_tuple(context, builder, sig, args):
    """
    Basic or advanced indexing with a tuple.
    """
    aryty, tupty = sig.args
    ary, tup = args
    ary = make_array(aryty)(context, builder, ary)

    index_types = tupty.types
    indices = cgutils.unpack_tuple(builder, tup, count=len(tupty))

    index_types, indices = normalize_indices(context, builder,
                                             index_types, indices)

    if any(isinstance(ty, types.Array) for ty in index_types):
        # Advanced indexing
        return fancy_getitem(context, builder, sig, args,
                             aryty, ary, index_types, indices)

    res = _getitem_array_generic(context, builder, sig.return_type,
                                 aryty, ary, index_types, indices)
    return impl_ret_borrowed(context, builder, sig.return_type, res)


@lower_builtin(operator.setitem, types.Buffer, types.Any, types.Any)
def setitem_array(context, builder, sig, args):
    """
    array[a] = scalar_or_array
    array[a,..,b] = scalar_or_array
    """
    aryty, idxty, valty = sig.args
    ary, idx, val = args

    if isinstance(idxty, types.BaseTuple):
        index_types = idxty.types
        indices = cgutils.unpack_tuple(builder, idx, count=len(idxty))
    else:
        index_types = (idxty,)
        indices = (idx,)

    ary = make_array(aryty)(context, builder, ary)

    # First try basic indexing to see if a single array location is denoted.
    index_types, indices = normalize_indices(context, builder,
                                             index_types, indices)
    try:
        dataptr, shapes, strides = \
            basic_indexing(context, builder, aryty, ary, index_types, indices,
                           boundscheck=context.enable_boundscheck)
    except NotImplementedError:
        use_fancy_indexing = True
    else:
        use_fancy_indexing = bool(shapes)

    if use_fancy_indexing:
        # Index describes a non-trivial view => use generic slice assignment
        # (NOTE: this also handles scalar broadcasting)
        return fancy_setslice(context, builder, sig, args,
                              index_types, indices)

    # Store source value the given location
    val = context.cast(builder, val, valty, aryty.dtype)
    store_item(context, builder, aryty, val, dataptr)


@lower_builtin(len, types.Buffer)
def array_len(context, builder, sig, args):
    (aryty,) = sig.args
    (ary,) = args
    arystty = make_array(aryty)
    ary = arystty(context, builder, ary)
    shapeary = ary.shape
    res = builder.extract_value(shapeary, 0)
    return impl_ret_untracked(context, builder, sig.return_type, res)


@lower_builtin("array.item", types.Array)
def array_item(context, builder, sig, args):
    aryty, = sig.args
    ary, = args
    ary = make_array(aryty)(context, builder, ary)

    nitems = ary.nitems
    with builder.if_then(builder.icmp_signed('!=', nitems, nitems.type(1)),
                         likely=False):
        msg = "item(): can only convert an array of size 1 to a Python scalar"
        context.call_conv.return_user_exc(builder, ValueError, (msg,))

    return load_item(context, builder, aryty, ary.data)


if numpy_version < (2, 0):
    @lower_builtin("array.itemset", types.Array, types.Any)
    def array_itemset(context, builder, sig, args):
        aryty, valty = sig.args
        ary, val = args
        assert valty == aryty.dtype
        ary = make_array(aryty)(context, builder, ary)

        nitems = ary.nitems
        with builder.if_then(builder.icmp_signed('!=', nitems, nitems.type(1)),
                             likely=False):
            msg = "itemset(): can only write to an array of size 1"
            context.call_conv.return_user_exc(builder, ValueError, (msg,))

        store_item(context, builder, aryty, val, ary.data)
        return context.get_dummy_value()


# ------------------------------------------------------------------------------
# Advanced / fancy indexing


class Indexer(object):
    """
    Generic indexer interface, for generating indices over a fancy indexed
    array on a single dimension.
    """

    def prepare(self):
        """
        Prepare the indexer by initializing any required variables, basic
        blocks...
        """
        raise NotImplementedError

    def get_size(self):
        """
        Return this dimension's size as an integer.
        """
        raise NotImplementedError

    def get_shape(self):
        """
        Return this dimension's shape as a tuple.
        """
        raise NotImplementedError

    def get_index_bounds(self):
        """
        Return a half-open [lower, upper) range of indices this dimension
        is guaranteed not to step out of.
        """
        raise NotImplementedError

    def loop_head(self):
        """
        Start indexation loop.  Return a (index, count) tuple.
        *index* is an integer LLVM value representing the index over this
        dimension.
        *count* is either an integer LLVM value representing the current
        iteration count, or None if this dimension should be omitted from
        the indexation result.
        """
        raise NotImplementedError

    def loop_tail(self):
        """
        Finish indexation loop.
        """
        raise NotImplementedError


class EntireIndexer(Indexer):
    """
    Compute indices along an entire array dimension.
    """

    def __init__(self, context, builder, aryty, ary, dim):
        self.context = context
        self.builder = builder
        self.aryty = aryty
        self.ary = ary
        self.dim = dim
        self.ll_intp = self.context.get_value_type(types.intp)

    def prepare(self):
        builder = self.builder
        self.size = builder.extract_value(self.ary.shape, self.dim)
        self.index = cgutils.alloca_once(builder, self.ll_intp)
        self.bb_start = builder.append_basic_block()
        self.bb_end = builder.append_basic_block()

    def get_size(self):
        return self.size

    def get_shape(self):
        return (self.size,)

    def get_index_bounds(self):
        # [0, size)
        return (self.ll_intp(0), self.size)

    def loop_head(self):
        builder = self.builder
        # Initialize loop variable
        self.builder.store(Constant(self.ll_intp, 0), self.index)
        builder.branch(self.bb_start)
        builder.position_at_end(self.bb_start)
        cur_index = builder.load(self.index)
        with builder.if_then(builder.icmp_signed('>=', cur_index, self.size),
                             likely=False):
            builder.branch(self.bb_end)
        return cur_index, cur_index

    def loop_tail(self):
        builder = self.builder
        next_index = cgutils.increment_index(builder, builder.load(self.index))
        builder.store(next_index, self.index)
        builder.branch(self.bb_start)
        builder.position_at_end(self.bb_end)


class IntegerIndexer(Indexer):
    """
    Compute indices from a single integer.
    """

    def __init__(self, context, builder, idx):
        self.context = context
        self.builder = builder
        self.idx = idx
        self.ll_intp = self.context.get_value_type(types.intp)

    def prepare(self):
        pass

    def get_size(self):
        return Constant(self.ll_intp, 1)

    def get_shape(self):
        return ()

    def get_index_bounds(self):
        # [idx, idx+1)
        return (self.idx, self.builder.add(self.idx, self.get_size()))

    def loop_head(self):
        return self.idx, None

    def loop_tail(self):
        pass


class IntegerArrayIndexer(Indexer):
    """
    Compute indices from an array of integer indices.
    """

    def __init__(self, context, builder, idxty, idxary, size):
        self.context = context
        self.builder = builder
        self.idxty = idxty
        self.idxary = idxary
        self.size = size
        assert idxty.ndim == 1
        self.ll_intp = self.context.get_value_type(types.intp)

    def prepare(self):
        builder = self.builder
        self.idx_size = cgutils.unpack_tuple(builder, self.idxary.shape)[0]
        self.idx_index = cgutils.alloca_once(builder, self.ll_intp)
        self.bb_start = builder.append_basic_block()
        self.bb_end = builder.append_basic_block()

    def get_size(self):
        return self.idx_size

    def get_shape(self):
        return (self.idx_size,)

    def get_index_bounds(self):
        # Pessimal heuristic, as we don't want to scan for the min and max
        return (self.ll_intp(0), self.size)

    def loop_head(self):
        builder = self.builder
        # Initialize loop variable
        self.builder.store(Constant(self.ll_intp, 0), self.idx_index)
        builder.branch(self.bb_start)
        builder.position_at_end(self.bb_start)
        cur_index = builder.load(self.idx_index)
        with builder.if_then(
            builder.icmp_signed('>=', cur_index, self.idx_size),
            likely=False
        ):
            builder.branch(self.bb_end)
        # Load the actual index from the array of indices
        index = _getitem_array_single_int(
            self.context, builder, self.idxty.dtype, self.idxty, self.idxary,
            cur_index
        )
        index = fix_integer_index(self.context, builder,
                                  self.idxty.dtype, index, self.size)
        return index, cur_index

    def loop_tail(self):
        builder = self.builder
        next_index = cgutils.increment_index(builder,
                                             builder.load(self.idx_index))
        builder.store(next_index, self.idx_index)
        builder.branch(self.bb_start)
        builder.position_at_end(self.bb_end)


class BooleanArrayIndexer(Indexer):
    """
    Compute indices from an array of boolean predicates.
    """

    def __init__(self, context, builder, idxty, idxary):
        self.context = context
        self.builder = builder
        self.idxty = idxty
        self.idxary = idxary
        assert idxty.ndim == 1
        self.ll_intp = self.context.get_value_type(types.intp)
        self.zero = Constant(self.ll_intp, 0)

    def prepare(self):
        builder = self.builder
        self.size = cgutils.unpack_tuple(builder, self.idxary.shape)[0]
        self.idx_index = cgutils.alloca_once(builder, self.ll_intp)
        self.count = cgutils.alloca_once(builder, self.ll_intp)
        self.bb_start = builder.append_basic_block()
        self.bb_tail = builder.append_basic_block()
        self.bb_end = builder.append_basic_block()

    def get_size(self):
        builder = self.builder
        count = cgutils.alloca_once_value(builder, self.zero)
        # Sum all true values
        with cgutils.for_range(builder, self.size) as loop:
            c = builder.load(count)
            pred = _getitem_array_single_int(
                self.context, builder, self.idxty.dtype,
                self.idxty, self.idxary, loop.index
            )
            c = builder.add(c, builder.zext(pred, c.type))
            builder.store(c, count)

        return builder.load(count)

    def get_shape(self):
        return (self.get_size(),)

    def get_index_bounds(self):
        # Pessimal heuristic, as we don't want to scan for the
        # first and last true items
        return (self.ll_intp(0), self.size)

    def loop_head(self):
        builder = self.builder
        # Initialize loop variable
        self.builder.store(self.zero, self.idx_index)
        self.builder.store(self.zero, self.count)
        builder.branch(self.bb_start)
        builder.position_at_end(self.bb_start)
        cur_index = builder.load(self.idx_index)
        cur_count = builder.load(self.count)
        with builder.if_then(builder.icmp_signed('>=', cur_index, self.size),
                             likely=False):
            builder.branch(self.bb_end)
        # Load the predicate and branch if false
        pred = _getitem_array_single_int(
            self.context, builder, self.idxty.dtype, self.idxty, self.idxary,
            cur_index
        )
        with builder.if_then(builder.not_(pred)):
            builder.branch(self.bb_tail)
        # Increment the count for next iteration
        next_count = cgutils.increment_index(builder, cur_count)
        builder.store(next_count, self.count)
        return cur_index, cur_count

    def loop_tail(self):
        builder = self.builder
        builder.branch(self.bb_tail)
        builder.position_at_end(self.bb_tail)
        next_index = cgutils.increment_index(builder,
                                             builder.load(self.idx_index))
        builder.store(next_index, self.idx_index)
        builder.branch(self.bb_start)
        builder.position_at_end(self.bb_end)


class SliceIndexer(Indexer):
    """
    Compute indices along a slice.
    """

    def __init__(self, context, builder, aryty, ary, dim, idxty, slice):
        self.context = context
        self.builder = builder
        self.aryty = aryty
        self.ary = ary
        self.dim = dim
        self.idxty = idxty
        self.slice = slice
        self.ll_intp = self.context.get_value_type(types.intp)
        self.zero = Constant(self.ll_intp, 0)

    def prepare(self):
        builder = self.builder
        # Fix slice for the dimension's size
        self.dim_size = builder.extract_value(self.ary.shape, self.dim)
        slicing.guard_invalid_slice(self.context, builder, self.idxty,
                                    self.slice)
        slicing.fix_slice(builder, self.slice, self.dim_size)
        self.is_step_negative = cgutils.is_neg_int(builder, self.slice.step)
        # Create loop entities
        self.index = cgutils.alloca_once(builder, self.ll_intp)
        self.count = cgutils.alloca_once(builder, self.ll_intp)
        self.bb_start = builder.append_basic_block()
        self.bb_end = builder.append_basic_block()

    def get_size(self):
        return slicing.get_slice_length(self.builder, self.slice)

    def get_shape(self):
        return (self.get_size(),)

    def get_index_bounds(self):
        lower, upper = slicing.get_slice_bounds(self.builder, self.slice)
        return lower, upper

    def loop_head(self):
        builder = self.builder
        # Initialize loop variable
        self.builder.store(self.slice.start, self.index)
        self.builder.store(self.zero, self.count)
        builder.branch(self.bb_start)
        builder.position_at_end(self.bb_start)
        cur_index = builder.load(self.index)
        cur_count = builder.load(self.count)
        is_finished = builder.select(self.is_step_negative,
                                     builder.icmp_signed('<=', cur_index,
                                                         self.slice.stop),
                                     builder.icmp_signed('>=', cur_index,
                                                         self.slice.stop))
        with builder.if_then(is_finished, likely=False):
            builder.branch(self.bb_end)
        return cur_index, cur_count

    def loop_tail(self):
        builder = self.builder
        next_index = builder.add(builder.load(self.index), self.slice.step,
                                 flags=['nsw'])
        builder.store(next_index, self.index)
        next_count = cgutils.increment_index(builder, builder.load(self.count))
        builder.store(next_count, self.count)
        builder.branch(self.bb_start)
        builder.position_at_end(self.bb_end)


class FancyIndexer(object):
    """
    Perform fancy indexing on the given array.
    """

    def __init__(self, context, builder, aryty, ary, index_types, indices):
        self.context = context
        self.builder = builder
        self.aryty = aryty
        self.shapes = cgutils.unpack_tuple(builder, ary.shape, aryty.ndim)
        self.strides = cgutils.unpack_tuple(builder, ary.strides, aryty.ndim)
        self.ll_intp = self.context.get_value_type(types.intp)
        self.newaxes = []

        indexers = []
        num_newaxes = len([idx for idx in index_types if is_nonelike(idx)])

        ax = 0 # keeps track of position of original axes
        new_ax = 0 # keeps track of position for inserting new axes
        for indexval, idxty in zip(indices, index_types):
            if idxty is types.ellipsis:
                # Fill up missing dimensions at the middle
                n_missing = aryty.ndim - len(indices) + 1 + num_newaxes
                for i in range(n_missing):
                    indexer = EntireIndexer(context, builder, aryty, ary, ax)
                    indexers.append(indexer)
                    ax += 1
                    new_ax += 1
                continue

            # Regular index value
            if isinstance(idxty, types.SliceType):
                slice = context.make_helper(builder, idxty, indexval)
                indexer = SliceIndexer(context, builder, aryty, ary, ax,
                                       idxty, slice)
                indexers.append(indexer)
            elif isinstance(idxty, types.Integer):
                ind = fix_integer_index(context, builder, idxty, indexval,
                                        self.shapes[ax])
                indexer = IntegerIndexer(context, builder, ind)
                indexers.append(indexer)
            elif isinstance(idxty, types.Array):
                idxary = make_array(idxty)(context, builder, indexval)
                if isinstance(idxty.dtype, types.Integer):
                    indexer = IntegerArrayIndexer(context, builder,
                                                  idxty, idxary,
                                                  self.shapes[ax])
                elif isinstance(idxty.dtype, types.Boolean):
                    indexer = BooleanArrayIndexer(context, builder,
                                                  idxty, idxary)
                else:
                    assert 0
                indexers.append(indexer)
            elif is_nonelike(idxty):
                self.newaxes.append(new_ax)
                ax -= 1
            else:
                raise AssertionError("unexpected index type: %s" % (idxty,))
            ax += 1
            new_ax += 1

        # Fill up missing dimensions at the end
        assert ax <= aryty.ndim, (ax, aryty.ndim)
        while ax < aryty.ndim:
            indexer = EntireIndexer(context, builder, aryty, ary, ax)
            indexers.append(indexer)
            ax += 1

        assert len(indexers) == aryty.ndim, (len(indexers), aryty.ndim)
        self.indexers = indexers

    def prepare(self):
        for i in self.indexers:
            i.prepare()

        one = self.context.get_constant(types.intp, 1)

        # Compute the resulting shape given by the indices
        res_shape = [i.get_shape() for i in self.indexers]

        # At every position where newaxis/None is present insert
        # one as a constant shape in the resulting list of shapes.
        for i in self.newaxes:
            res_shape.insert(i, (one,))

        # Store the shape as a tuple, we can't do a simple
        # tuple(res_shape) here since res_shape is a list
        # of tuples which may be differently sized.
        self.indexers_shape = sum(res_shape, ())

    def get_shape(self):
        """
        Get the resulting data shape as Python tuple.
        """
        return self.indexers_shape

    def get_offset_bounds(self, strides, itemsize):
        """
        Get a half-open [lower, upper) range of byte offsets spanned by
        the indexer with the given strides and itemsize.  The indexer is
        guaranteed to not go past those bounds.
        """
        assert len(strides) == self.aryty.ndim
        builder = self.builder
        is_empty = cgutils.false_bit
        zero = self.ll_intp(0)
        one = self.ll_intp(1)
        lower = zero
        upper = zero
        for indexer, shape, stride in zip(self.indexers, self.indexers_shape,
                                          strides):
            is_empty = builder.or_(is_empty,
                                   builder.icmp_unsigned('==', shape, zero))
            # Compute [lower, upper) indices on this dimension
            lower_index, upper_index = indexer.get_index_bounds()
            lower_offset = builder.mul(stride, lower_index)
            upper_offset = builder.mul(stride, builder.sub(upper_index, one))
            # Adjust total interval
            is_downwards = builder.icmp_signed('<', stride, zero)
            lower = builder.add(lower,
                                builder.select(is_downwards,
                                               upper_offset,
                                               lower_offset))
            upper = builder.add(upper,
                                builder.select(is_downwards,
                                               lower_offset,
                                               upper_offset))
        # Make interval half-open
        upper = builder.add(upper, itemsize)
        # Adjust for empty shape
        lower = builder.select(is_empty, zero, lower)
        upper = builder.select(is_empty, zero, upper)
        return lower, upper

    def begin_loops(self):
        indices, counts = zip(*(i.loop_head() for i in self.indexers))
        return indices, counts

    def end_loops(self):
        for i in reversed(self.indexers):
            i.loop_tail()


def fancy_getitem(context, builder, sig, args,
                  aryty, ary, index_types, indices):

    shapes = cgutils.unpack_tuple(builder, ary.shape)
    strides = cgutils.unpack_tuple(builder, ary.strides)
    data = ary.data

    indexer = FancyIndexer(context, builder, aryty, ary,
                           index_types, indices)
    indexer.prepare()

    # Construct output array
    out_ty = sig.return_type
    out_shapes = indexer.get_shape()

    out = _empty_nd_impl(context, builder, out_ty, out_shapes)
    out_data = out.data
    out_idx = cgutils.alloca_once_value(builder,
                                        context.get_constant(types.intp, 0))

    # Loop on source and copy to destination
    indices, _ = indexer.begin_loops()

    # No need to check for wraparound, as the indexers all ensure
    # a positive index is returned.
    ptr = cgutils.get_item_pointer2(context, builder, data, shapes, strides,
                                    aryty.layout, indices, wraparound=False,
                                    boundscheck=context.enable_boundscheck)
    val = load_item(context, builder, aryty, ptr)

    # Since the destination is C-contiguous, no need for multi-dimensional
    # indexing.
    cur = builder.load(out_idx)
    ptr = builder.gep(out_data, [cur])
    store_item(context, builder, out_ty, val, ptr)
    next_idx = cgutils.increment_index(builder, cur)
    builder.store(next_idx, out_idx)

    indexer.end_loops()

    return impl_ret_new_ref(context, builder, out_ty, out._getvalue())


@lower_builtin(operator.getitem, types.Buffer, types.Array)
def fancy_getitem_array(context, builder, sig, args):
    """
    Advanced or basic indexing with an array.
    """
    aryty, idxty = sig.args
    ary, idx = args
    ary = make_array(aryty)(context, builder, ary)
    if idxty.ndim == 0:
        # 0-d array index acts as a basic integer index
        idxty, idx = normalize_index(context, builder, idxty, idx)
        res = _getitem_array_generic(context, builder, sig.return_type,
                                     aryty, ary, (idxty,), (idx,))
        return impl_ret_borrowed(context, builder, sig.return_type, res)
    else:
        # Advanced indexing
        return fancy_getitem(context, builder, sig, args,
                             aryty, ary, (idxty,), (idx,))


def offset_bounds_from_strides(context, builder, arrty, arr, shapes, strides):
    """
    Compute a half-open range [lower, upper) of byte offsets from the
    array's data pointer, that bound the in-memory extent of the array.

    This mimics offset_bounds_from_strides() from
    numpy/core/src/private/mem_overlap.c
    """
    itemsize = arr.itemsize
    zero = itemsize.type(0)
    one = zero.type(1)
    if arrty.layout in 'CF':
        # Array is contiguous: contents are laid out sequentially
        # starting from arr.data and upwards
        lower = zero
        upper = builder.mul(itemsize, arr.nitems)
    else:
        # Non-contiguous array: need to examine strides
        lower = zero
        upper = zero
        for i in range(arrty.ndim):
            # Compute the largest byte offset on this dimension
            #   max_axis_offset = strides[i] * (shapes[i] - 1)
            # (shapes[i] == 0 is catered for by the empty array case below)
            max_axis_offset = builder.mul(strides[i],
                                          builder.sub(shapes[i], one))
            is_upwards = builder.icmp_signed('>=', max_axis_offset, zero)
            # Expand either upwards or downwards depending on stride
            upper = builder.select(is_upwards,
                                   builder.add(upper, max_axis_offset), upper)
            lower = builder.select(is_upwards,
                                   lower, builder.add(lower, max_axis_offset))
        # Return a half-open range
        upper = builder.add(upper, itemsize)
        # Adjust for empty arrays
        is_empty = builder.icmp_signed('==', arr.nitems, zero)
        upper = builder.select(is_empty, zero, upper)
        lower = builder.select(is_empty, zero, lower)

    return lower, upper


def compute_memory_extents(context, builder, lower, upper, data):
    """
    Given [lower, upper) byte offsets and a base data pointer,
    compute the memory pointer bounds as pointer-sized integers.
    """
    data_ptr_as_int = builder.ptrtoint(data, lower.type)
    start = builder.add(data_ptr_as_int, lower)
    end = builder.add(data_ptr_as_int, upper)
    return start, end


def get_array_memory_extents(context, builder, arrty, arr, shapes, strides,
                             data):
    """
    Compute a half-open range [start, end) of pointer-sized integers
    which fully contain the array data.
    """
    lower, upper = offset_bounds_from_strides(context, builder, arrty, arr,
                                              shapes, strides)
    return compute_memory_extents(context, builder, lower, upper, data)


def extents_may_overlap(context, builder, a_start, a_end, b_start, b_end):
    """
    Whether two memory extents [a_start, a_end) and [b_start, b_end)
    may overlap.
    """
    # Comparisons are unsigned, since we are really comparing pointers
    may_overlap = builder.and_(
        builder.icmp_unsigned('<', a_start, b_end),
        builder.icmp_unsigned('<', b_start, a_end),
    )
    return may_overlap


def maybe_copy_source(context, builder, use_copy,
                      srcty, src, src_shapes, src_strides, src_data):
    ptrty = src_data.type

    copy_layout = 'C'
    copy_data = cgutils.alloca_once_value(builder, src_data)
    copy_shapes = src_shapes
    copy_strides = None  # unneeded for contiguous arrays

    with builder.if_then(use_copy, likely=False):
        # Allocate temporary scratchpad
        # XXX: should we use a stack-allocated array for very small
        # data sizes?
        allocsize = builder.mul(src.itemsize, src.nitems)
        data = context.nrt.allocate(builder, allocsize)
        voidptrty = data.type
        data = builder.bitcast(data, ptrty)
        builder.store(data, copy_data)

        # Copy source data into scratchpad
        intp_t = context.get_value_type(types.intp)

        with cgutils.loop_nest(builder, src_shapes, intp_t) as indices:
            src_ptr = cgutils.get_item_pointer2(context, builder, src_data,
                                                src_shapes, src_strides,
                                                srcty.layout, indices)
            dest_ptr = cgutils.get_item_pointer2(context, builder, data,
                                                 copy_shapes, copy_strides,
                                                 copy_layout, indices)
            builder.store(builder.load(src_ptr), dest_ptr)

    def src_getitem(source_indices):
        src_ptr = cgutils.alloca_once(builder, ptrty)
        with builder.if_else(use_copy, likely=False) as (if_copy, otherwise):
            with if_copy:
                builder.store(
                    cgutils.get_item_pointer2(context, builder,
                                              builder.load(copy_data),
                                              copy_shapes, copy_strides,
                                              copy_layout, source_indices,
                                              wraparound=False),
                    src_ptr)
            with otherwise:
                builder.store(
                    cgutils.get_item_pointer2(context, builder, src_data,
                                              src_shapes, src_strides,
                                              srcty.layout, source_indices,
                                              wraparound=False),
                    src_ptr)
        return load_item(context, builder, srcty, builder.load(src_ptr))

    def src_cleanup():
        # Deallocate memory
        with builder.if_then(use_copy, likely=False):
            data = builder.load(copy_data)
            data = builder.bitcast(data, voidptrty)
            context.nrt.free(builder, data)

    return src_getitem, src_cleanup


def _bc_adjust_dimension(context, builder, shapes, strides, target_shape):
    """
    Preprocess dimension for broadcasting.
    Returns (shapes, strides) such that the ndim match *target_shape*.
    When expanding to higher ndim, the returning shapes and strides are
    prepended with ones and zeros, respectively.
    When truncating to lower ndim, the shapes are checked (in runtime).
    All extra dimension must have size of 1.
    """
    zero = context.get_constant(types.uintp, 0)
    one = context.get_constant(types.uintp, 1)

    # Adjust for broadcasting to higher dimension
    if len(target_shape) > len(shapes):
        nd_diff = len(target_shape) - len(shapes)
        # Fill missing shapes with one, strides with zeros
        shapes = [one] * nd_diff + shapes
        strides = [zero] * nd_diff + strides
    # Adjust for broadcasting to lower dimension
    elif len(target_shape) < len(shapes):
        # Accepted if all extra dims has shape 1
        nd_diff = len(shapes) - len(target_shape)
        dim_is_one = [builder.icmp_unsigned('==', sh, one)
                      for sh in shapes[:nd_diff]]
        accepted = functools.reduce(builder.and_, dim_is_one,
                                    cgutils.true_bit)
        # Check error
        with builder.if_then(builder.not_(accepted), likely=False):
            msg = "cannot broadcast source array for assignment"
            context.call_conv.return_user_exc(builder, ValueError, (msg,))
        # Truncate extra shapes, strides
        shapes = shapes[nd_diff:]
        strides = strides[nd_diff:]

    return shapes, strides


def _bc_adjust_shape_strides(context, builder, shapes, strides, target_shape):
    """
    Broadcast shapes and strides to target_shape given that their ndim already
    matches.  For each location where the shape is 1 and does not match the
    dim for target, it is set to the value at the target and the stride is
    set to zero.
    """
    bc_shapes = []
    bc_strides = []
    zero = context.get_constant(types.uintp, 0)
    one = context.get_constant(types.uintp, 1)
    # Adjust all mismatching ones in shape
    mismatch = [builder.icmp_signed('!=', tar, old)
                for tar, old in zip(target_shape, shapes)]
    src_is_one = [builder.icmp_signed('==', old, one) for old in shapes]
    preds = [builder.and_(x, y) for x, y in zip(mismatch, src_is_one)]
    bc_shapes = [builder.select(p, tar, old)
                 for p, tar, old in zip(preds, target_shape, shapes)]
    bc_strides = [builder.select(p, zero, old)
                  for p, old in zip(preds, strides)]
    return bc_shapes, bc_strides


def _broadcast_to_shape(context, builder, arrtype, arr, target_shape):
    """
    Broadcast the given array to the target_shape.
    Returns (array_type, array)
    """
    # Compute broadcasted shape and strides
    shapes = cgutils.unpack_tuple(builder, arr.shape)
    strides = cgutils.unpack_tuple(builder, arr.strides)

    shapes, strides = _bc_adjust_dimension(context, builder, shapes, strides,
                                           target_shape)
    shapes, strides = _bc_adjust_shape_strides(context, builder, shapes,
                                               strides, target_shape)
    new_arrtype = arrtype.copy(ndim=len(target_shape), layout='A')
    # Create new view
    new_arr = make_array(new_arrtype)(context, builder)
    populate_array(new_arr,
                   data=arr.data,
                   shape=cgutils.pack_array(builder, shapes),
                   strides=cgutils.pack_array(builder, strides),
                   itemsize=arr.itemsize,
                   meminfo=arr.meminfo,
                   parent=arr.parent)
    return new_arrtype, new_arr


@intrinsic
def _numpy_broadcast_to(typingctx, array, shape):
    ret = array.copy(ndim=shape.count, layout='A', readonly=True)
    sig = ret(array, shape)

    def codegen(context, builder, sig, args):
        src, shape_ = args
        srcty = sig.args[0]

        src = make_array(srcty)(context, builder, src)
        shape_ = cgutils.unpack_tuple(builder, shape_)
        _, dest = _broadcast_to_shape(context, builder, srcty, src, shape_,)

        # Hack to get np.broadcast_to to return a read-only array
        setattr(dest, 'parent', Constant(
                context.get_value_type(dest._datamodel.get_type('parent')),
                None))

        res = dest._getvalue()
        return impl_ret_borrowed(context, builder, sig.return_type, res)
    return sig, codegen


@intrinsic
def get_readonly_array(typingctx, arr):
    # returns a copy of arr which is readonly
    ret = arr.copy(readonly=True)
    sig = ret(arr)

    def codegen(context, builder, sig, args):
        [src] = args
        srcty = sig.args[0]

        dest = make_array(srcty)(context, builder, src)
        # Hack to return a read-only array
        dest.parent = cgutils.get_null_value(dest.parent.type)
        res = dest._getvalue()
        return impl_ret_borrowed(context, builder, sig.return_type, res)
    return sig, codegen


@register_jitable
def _can_broadcast(array, dest_shape):
    src_shape = array.shape
    src_ndim = len(src_shape)
    dest_ndim = len(dest_shape)
    if src_ndim > dest_ndim:
        raise ValueError('input operand has more dimensions than allowed '
                         'by the axis remapping')
    for size in dest_shape:
        if size < 0:
            raise ValueError('all elements of broadcast shape must be '
                             'non-negative')

    # based on _broadcast_onto function in numba/np/npyimpl.py
    src_index = 0
    dest_index = dest_ndim - src_ndim
    while src_index < src_ndim:
        src_dim = src_shape[src_index]
        dest_dim = dest_shape[dest_index]
        # possible cases for (src_dim, dest_dim):
        #  * (1, 1)   -> Ok
        #  * (>1, 1)  -> Error!
        #  * (>1, >1) -> src_dim == dest_dim else error!
        #  * (1, >1)  -> Ok
        if src_dim == dest_dim or src_dim == 1:
            src_index += 1
            dest_index += 1
        else:
            raise ValueError('operands could not be broadcast together '
                             'with remapped shapes')


def _default_broadcast_to_impl(array, shape):
    array = np.asarray(array)
    _can_broadcast(array, shape)
    return _numpy_broadcast_to(array, shape)


@overload(np.broadcast_to)
def numpy_broadcast_to(array, shape):
    if not type_can_asarray(array):
        raise errors.TypingError('The first argument "array" must '
                                 'be array-like')

    if isinstance(shape, types.Integer):
        def impl(array, shape):
            return np.broadcast_to(array, (shape,))
        return impl

    elif isinstance(shape, types.UniTuple):
        if not isinstance(shape.dtype, types.Integer):
            msg = 'The second argument "shape" must be a tuple of integers'
            raise errors.TypingError(msg)
        return _default_broadcast_to_impl

    elif isinstance(shape, types.Tuple) and shape.count > 0:
        # check if all types are integers
        if not all([isinstance(typ, types.IntegerLiteral) for typ in shape]):
            msg = f'"{shape}" object cannot be interpreted as an integer'
            raise errors.TypingError(msg)
        return _default_broadcast_to_impl
    elif isinstance(shape, types.Tuple) and shape.count == 0:
        is_scalar_array = isinstance(array, types.Array) and array.ndim == 0
        if type_is_scalar(array) or is_scalar_array:

            def impl(array, shape):  # broadcast_to(array, ())
                # Array type must be supported by "type_can_asarray"
                # Quick note that unicode types are not supported!
                array = np.asarray(array)
                return get_readonly_array(array)
            return impl

        else:
            msg = 'Cannot broadcast a non-scalar to a scalar array'
            raise errors.TypingError(msg)
    else:
        msg = ('The argument "shape" must be a tuple or an integer. '
               'Got %s' % shape)
        raise errors.TypingError(msg)


@register_jitable
def numpy_broadcast_shapes_list(r, m, shape):
    for i in range(len(shape)):
        k = m - len(shape) + i
        tmp = shape[i]
        if tmp < 0:
            raise ValueError("negative dimensions are not allowed")
        if tmp == 1:
            continue
        if r[k] == 1:
            r[k] = tmp
        elif r[k] != tmp:
            raise ValueError("shape mismatch: objects"
                             " cannot be broadcast"
                             " to a single shape")


@overload(np.broadcast_shapes)
def ol_numpy_broadcast_shapes(*args):
    # Based on https://github.com/numpy/numpy/blob/f702b26fff3271ba6a6ba29a021fc19051d1f007/numpy/core/src/multiarray/iterators.c#L1129-L1212  # noqa
    for idx, arg in enumerate(args):
        is_int = isinstance(arg, types.Integer)
        is_int_tuple = isinstance(arg, types.UniTuple) and \
            isinstance(arg.dtype, types.Integer)
        is_empty_tuple = isinstance(arg, types.Tuple) and len(arg.types) == 0
        if not (is_int or is_int_tuple or is_empty_tuple):
            msg = (f'Argument {idx} must be either an int or tuple[int]. '
                   f'Got {arg}')
            raise errors.TypingError(msg)

    # discover the number of dimensions
    m = 0
    for arg in args:
        if isinstance(arg, types.Integer):
            m = max(m, 1)
        elif isinstance(arg, types.BaseTuple):
            m = max(m, len(arg))

    if m == 0:
        return lambda *args: ()
    else:
        tup_init = (1,) * m

        def impl(*args):
            # propagate args
            r = [1] * m
            tup = tup_init
            for arg in literal_unroll(args):
                if isinstance(arg, tuple) and len(arg) > 0:
                    numpy_broadcast_shapes_list(r, m, arg)
                elif isinstance(arg, int):
                    numpy_broadcast_shapes_list(r, m, (arg,))
            for idx, elem in enumerate(r):
                tup = tuple_setitem(tup, idx, elem)
            return tup
        return impl


@overload(np.broadcast_arrays)
def numpy_broadcast_arrays(*args):

    for idx, arg in enumerate(args):
        if not type_can_asarray(arg):
            raise errors.TypingError(f'Argument "{idx}" must '
                                     'be array-like')

    unified_dtype = None
    dt = None
    for arg in args:
        if isinstance(arg, (types.Array, types.BaseTuple)):
            dt = arg.dtype
        else:
            dt = arg

        if unified_dtype is None:
            unified_dtype = dt
        elif unified_dtype != dt:
            raise errors.TypingError('Mismatch of argument types. Numba cannot '
                                     'broadcast arrays with different types. '
                                     f'Got {args}')

    # number of dimensions
    m = 0
    for idx, arg in enumerate(args):
        if isinstance(arg, types.ArrayCompatible):
            m = max(m, arg.ndim)
        elif isinstance(arg, (types.Number, types.Boolean, types.BaseTuple)):
            m = max(m, 1)
        else:
            raise errors.TypingError(f'Unhandled type {arg}')

    tup_init = (0,) * m

    def impl(*args):
        # find out the output shape
        # we can't call np.broadcast_shapes here since args may have arrays
        # with different shapes and it is not possible to create a list
        # with those shapes dynamically
        shape = [1] * m
        for array in literal_unroll(args):
            numpy_broadcast_shapes_list(shape, m, np.asarray(array).shape)

        tup = tup_init

        for i in range(m):
            tup = tuple_setitem(tup, i, shape[i])

        # numpy checks if the input arrays have the same shape as `shape`
        outs = []
        for array in literal_unroll(args):
            outs.append(np.broadcast_to(np.asarray(array), tup))
        return outs

    return impl


def fancy_setslice(context, builder, sig, args, index_types, indices):
    """
    Implement slice assignment for arrays.  This implementation works for
    basic as well as fancy indexing, since there's no functional difference
    between the two for indexed assignment.
    """
    aryty, _, srcty = sig.args
    ary, _, src = args

    ary = make_array(aryty)(context, builder, ary)
    dest_shapes = cgutils.unpack_tuple(builder, ary.shape)
    dest_strides = cgutils.unpack_tuple(builder, ary.strides)
    dest_data = ary.data

    indexer = FancyIndexer(context, builder, aryty, ary,
                           index_types, indices)
    indexer.prepare()

    if isinstance(srcty, types.Buffer):
        # Source is an array
        src_dtype = srcty.dtype
        index_shape = indexer.get_shape()
        src = make_array(srcty)(context, builder, src)
        # Broadcast source array to shape
        srcty, src = _broadcast_to_shape(context, builder, srcty, src,
                                         index_shape)
        src_shapes = cgutils.unpack_tuple(builder, src.shape)
        src_strides = cgutils.unpack_tuple(builder, src.strides)
        src_data = src.data

        # Check shapes are equal
        shape_error = cgutils.false_bit
        assert len(index_shape) == len(src_shapes)

        for u, v in zip(src_shapes, index_shape):
            shape_error = builder.or_(shape_error,
                                      builder.icmp_signed('!=', u, v))

        with builder.if_then(shape_error, likely=False):
            msg = "cannot assign slice from input of different size"
            context.call_conv.return_user_exc(builder, ValueError, (msg,))

        # Check for array overlap
        src_start, src_end = get_array_memory_extents(context, builder, srcty,
                                                      src, src_shapes,
                                                      src_strides, src_data)

        dest_lower, dest_upper = indexer.get_offset_bounds(dest_strides,
                                                           ary.itemsize)
        dest_start, dest_end = compute_memory_extents(context, builder,
                                                      dest_lower, dest_upper,
                                                      dest_data)

        use_copy = extents_may_overlap(context, builder, src_start, src_end,
                                       dest_start, dest_end)

        src_getitem, src_cleanup = maybe_copy_source(context, builder, use_copy,
                                                     srcty, src, src_shapes,
                                                     src_strides, src_data)

    elif isinstance(srcty, types.Sequence):
        src_dtype = srcty.dtype

        # Check shape is equal to sequence length
        index_shape = indexer.get_shape()
        assert len(index_shape) == 1
        len_impl = context.get_function(len, signature(types.intp, srcty))
        seq_len = len_impl(builder, (src,))

        shape_error = builder.icmp_signed('!=', index_shape[0], seq_len)

        with builder.if_then(shape_error, likely=False):
            msg = "cannot assign slice from input of different size"
            context.call_conv.return_user_exc(builder, ValueError, (msg,))

        def src_getitem(source_indices):
            idx, = source_indices
            getitem_impl = context.get_function(
                operator.getitem,
                signature(src_dtype, srcty, types.intp),
            )
            return getitem_impl(builder, (src, idx))

        def src_cleanup():
            pass

    else:
        # Source is a scalar (broadcast or not, depending on destination
        # shape).
        src_dtype = srcty

        def src_getitem(source_indices):
            return src

        def src_cleanup():
            pass

    zero = context.get_constant(types.uintp, 0)
    # Loop on destination and copy from source to destination
    dest_indices, counts = indexer.begin_loops()

    # Source is iterated in natural order

    # Counts represent a counter for the number of times a specified axis
    # is being accessed, during setitem they are used as source
    # indices
    counts = list(counts)

    # We need to artifically introduce the index zero wherever a
    # newaxis is present within the indexer. These always remain
    # zero.
    for i in indexer.newaxes:
        counts.insert(i, zero)

    source_indices = [c for c in counts if c is not None]

    val = src_getitem(source_indices)

    # Cast to the destination dtype (cross-dtype slice assignment is allowed)
    val = context.cast(builder, val, src_dtype, aryty.dtype)

    # No need to check for wraparound, as the indexers all ensure
    # a positive index is returned.
    dest_ptr = cgutils.get_item_pointer2(context, builder, dest_data,
                                         dest_shapes, dest_strides,
                                         aryty.layout, dest_indices,
                                         wraparound=False,
                                         boundscheck=context.enable_boundscheck)
    store_item(context, builder, aryty, val, dest_ptr)

    indexer.end_loops()

    src_cleanup()

    return context.get_dummy_value()


# ------------------------------------------------------------------------------
# Shape / layout altering

def vararg_to_tuple(context, builder, sig, args):
    aryty = sig.args[0]
    dimtys = sig.args[1:]
    # values
    ary = args[0]
    dims = args[1:]
    # coerce all types to intp
    dims = [context.cast(builder, val, ty, types.intp)
            for ty, val in zip(dimtys, dims)]
    # make a tuple
    shape = cgutils.pack_array(builder, dims, dims[0].type)

    shapety = types.UniTuple(dtype=types.intp, count=len(dims))
    new_sig = typing.signature(sig.return_type, aryty, shapety)
    new_args = ary, shape

    return new_sig, new_args


@lower_builtin('array.transpose', types.Array)
def array_transpose(context, builder, sig, args):
    return array_T(context, builder, sig.args[0], args[0])


def permute_arrays(axis, shape, strides):
    if len(axis) != len(set(axis)):
        raise ValueError("repeated axis in transpose")
    dim = len(shape)
    for x in axis:
        if x >= dim or abs(x) > dim:
            raise ValueError("axis is out of bounds for array of "
                             "given dimension")

    shape[:] = shape[axis]
    strides[:] = strides[axis]


# Transposing an array involves permuting the shape and strides of the array
# based on the given axes.
@lower_builtin('array.transpose', types.Array, types.BaseTuple)
def array_transpose_tuple(context, builder, sig, args):
    aryty = sig.args[0]
    ary = make_array(aryty)(context, builder, args[0])

    axisty, axis = sig.args[1], args[1]
    num_axis, dtype = axisty.count, axisty.dtype

    ll_intp = context.get_value_type(types.intp)
    ll_ary_size = ir.ArrayType(ll_intp, num_axis)

    # Allocate memory for axes, shapes, and strides arrays.
    arys = [axis, ary.shape, ary.strides]
    ll_arys = [cgutils.alloca_once(builder, ll_ary_size) for _ in arys]

    # Store axes, shapes, and strides arrays to the allocated memory.
    for src, dst in zip(arys, ll_arys):
        builder.store(src, dst)

    np_ary_ty = types.Array(dtype=dtype, ndim=1, layout='C')
    np_itemsize = context.get_constant(types.intp,
                                       context.get_abi_sizeof(ll_intp))

    # Form NumPy arrays for axes, shapes, and strides arrays.
    np_arys = [make_array(np_ary_ty)(context, builder) for _ in arys]

    # Roughly, `np_ary = np.array(ll_ary)` for each of axes, shapes, and strides
    for np_ary, ll_ary in zip(np_arys, ll_arys):
        populate_array(np_ary,
                       data=builder.bitcast(ll_ary, ll_intp.as_pointer()),
                       shape=[context.get_constant(types.intp, num_axis)],
                       strides=[np_itemsize],
                       itemsize=np_itemsize,
                       meminfo=None)

    # Pass NumPy arrays formed above to permute_arrays function that permutes
    # shapes and strides based on axis contents.
    context.compile_internal(builder, permute_arrays,
                             typing.signature(types.void,
                                              np_ary_ty, np_ary_ty, np_ary_ty),
                             [a._getvalue() for a in np_arys])

    # Make a new array based on permuted shape and strides and return it.
    ret = make_array(sig.return_type)(context, builder)
    populate_array(ret,
                   data=ary.data,
                   shape=builder.load(ll_arys[1]),
                   strides=builder.load(ll_arys[2]),
                   itemsize=ary.itemsize,
                   meminfo=ary.meminfo,
                   parent=ary.parent)
    res = ret._getvalue()
    return impl_ret_borrowed(context, builder, sig.return_type, res)


@lower_builtin('array.transpose', types.Array, types.VarArg(types.Any))
def array_transpose_vararg(context, builder, sig, args):
    new_sig, new_args = vararg_to_tuple(context, builder, sig, args)
    return array_transpose_tuple(context, builder, new_sig, new_args)


@overload(np.transpose)
def numpy_transpose(a, axes=None):
    if isinstance(a, types.BaseTuple):
        raise errors.UnsupportedError("np.transpose does not accept tuples")

    if axes is None:
        def np_transpose_impl(a, axes=None):
            return a.transpose()
    else:
        def np_transpose_impl(a, axes=None):
            return a.transpose(axes)

    return np_transpose_impl


@lower_getattr(types.Array, 'T')
def array_T(context, builder, typ, value):
    if typ.ndim <= 1:
        res = value
    else:
        ary = make_array(typ)(context, builder, value)
        ret = make_array(typ)(context, builder)
        shapes = cgutils.unpack_tuple(builder, ary.shape, typ.ndim)
        strides = cgutils.unpack_tuple(builder, ary.strides, typ.ndim)
        populate_array(ret,
                       data=ary.data,
                       shape=cgutils.pack_array(builder, shapes[::-1]),
                       strides=cgutils.pack_array(builder, strides[::-1]),
                       itemsize=ary.itemsize,
                       meminfo=ary.meminfo,
                       parent=ary.parent)
        res = ret._getvalue()
    return impl_ret_borrowed(context, builder, typ, res)


@overload(np.logspace)
def numpy_logspace(start, stop, num=50):
    if not isinstance(start, types.Number):
        raise errors.TypingError('The first argument "start" must be a number')
    if not isinstance(stop, types.Number):
        raise errors.TypingError('The second argument "stop" must be a number')
    if not isinstance(num, (int, types.Integer)):
        raise errors.TypingError('The third argument "num" must be an integer')

    def impl(start, stop, num=50):
        y = np.linspace(start, stop, num)
        return np.power(10.0, y)

    return impl


@overload(np.geomspace)
def numpy_geomspace(start, stop, num=50):
    if not isinstance(start, types.Number):
        msg = 'The argument "start" must be a number'
        raise errors.TypingError(msg)

    if not isinstance(stop, types.Number):
        msg = 'The argument "stop" must be a number'
        raise errors.TypingError(msg)

    if not isinstance(num, (int, types.Integer)):
        msg = 'The argument "num" must be an integer'
        raise errors.TypingError(msg)

    if any(isinstance(arg, types.Complex) for arg in [start, stop]):
        result_dtype = from_dtype(np.result_type(as_dtype(start),
                                                 as_dtype(stop), None))

        def impl(start, stop, num=50):
            if start == 0 or stop == 0:
                raise ValueError('Geometric sequence cannot include zero')
            start = result_dtype(start)
            stop = result_dtype(stop)
            both_imaginary = (start.real == 0) & (stop.real == 0)
            both_negative = (np.sign(start) == -1) & (np.sign(stop) == -1)
            out_sign = 1
            if both_imaginary:
                start = start.imag
                stop = stop.imag
                out_sign = 1j
            if both_negative:
                start = -start
                stop = -stop
                out_sign = -out_sign
            logstart = np.log10(start)
            logstop = np.log10(stop)
            result = np.logspace(logstart, logstop, num)
            # Make sure the endpoints match the start and stop arguments.
            # This is necessary because np.exp(np.log(x)) is not necessarily
            # equal to x.
            if num > 0:
                result[0] = start
                if num > 1:
                    result[-1] = stop
            return out_sign * result

    else:
        def impl(start, stop, num=50):
            if start == 0 or stop == 0:
                raise ValueError('Geometric sequence cannot include zero')
            both_negative = (np.sign(start) == -1) & (np.sign(stop) == -1)
            out_sign = 1
            if both_negative:
                start = -start
                stop = -stop
                out_sign = -out_sign
            logstart = np.log10(start)
            logstop = np.log10(stop)
            result = np.logspace(logstart, logstop, num)
            # Make sure the endpoints match the start and stop arguments.
            # This is necessary because np.exp(np.log(x)) is not necessarily
            # equal to x.
            if num > 0:
                result[0] = start
                if num > 1:
                    result[-1] = stop
            return out_sign * result

    return impl


@overload(np.rot90)
def numpy_rot90(m, k=1):
    # supporting axes argument it needs to be included in np.flip
    if not isinstance(k, (int, types.Integer)):
        raise errors.TypingError('The second argument "k" must be an integer')
    if not isinstance(m, types.Array):
        raise errors.TypingError('The first argument "m" must be an array')

    if m.ndim < 2:
        raise errors.NumbaValueError('Input must be >= 2-d.')

    def impl(m, k=1):
        k = k % 4
        if k == 0:
            return m[:]
        elif k == 1:
            return np.swapaxes(np.fliplr(m), 0, 1)
        elif k == 2:
            return np.flipud(np.fliplr(m))
        elif k == 3:
            return np.fliplr(np.swapaxes(m, 0, 1))
        else:
            raise AssertionError  # unreachable

    return impl


def _attempt_nocopy_reshape(context, builder, aryty, ary,
                            newnd, newshape, newstrides):
    """
    Call into Numba_attempt_nocopy_reshape() for the given array type
    and instance, and the specified new shape.

    Return value is non-zero if successful, and the array pointed to
    by *newstrides* will be filled up with the computed results.
    """
    ll_intp = context.get_value_type(types.intp)
    ll_intp_star = ll_intp.as_pointer()
    ll_intc = context.get_value_type(types.intc)
    fnty = ir.FunctionType(ll_intc, [
        # nd, *dims, *strides
        ll_intp, ll_intp_star, ll_intp_star,
        # newnd, *newdims, *newstrides
        ll_intp, ll_intp_star, ll_intp_star,
        # itemsize, is_f_order
        ll_intp, ll_intc])
    fn = cgutils.get_or_insert_function(builder.module, fnty,
                                        "numba_attempt_nocopy_reshape")

    nd = ll_intp(aryty.ndim)
    shape = cgutils.gep_inbounds(builder, ary._get_ptr_by_name('shape'), 0, 0)
    strides = cgutils.gep_inbounds(builder, ary._get_ptr_by_name('strides'),
                                   0, 0)
    newnd = ll_intp(newnd)
    newshape = cgutils.gep_inbounds(builder, newshape, 0, 0)
    newstrides = cgutils.gep_inbounds(builder, newstrides, 0, 0)
    is_f_order = ll_intc(0)
    res = builder.call(fn, [nd, shape, strides,
                            newnd, newshape, newstrides,
                            ary.itemsize, is_f_order])
    return res


def normalize_reshape_value(origsize, shape):
    num_neg_value = 0
    known_size = 1
    for ax, s in enumerate(shape):
        if s < 0:
            num_neg_value += 1
            neg_ax = ax
        else:
            known_size *= s

    if num_neg_value == 0:
        if origsize != known_size:
            raise ValueError("total size of new array must be unchanged")

    elif num_neg_value == 1:
        # Infer negative dimension
        if known_size == 0:
            inferred = 0
            ok = origsize == 0
        else:
            inferred = origsize // known_size
            ok = origsize % known_size == 0
        if not ok:
            raise ValueError("total size of new array must be unchanged")
        shape[neg_ax] = inferred

    else:
        raise ValueError("multiple negative shape values")


@lower_builtin('array.reshape', types.Array, types.BaseTuple)
def array_reshape(context, builder, sig, args):
    aryty = sig.args[0]
    retty = sig.return_type

    shapety = sig.args[1]
    shape = args[1]

    ll_intp = context.get_value_type(types.intp)
    ll_shape = ir.ArrayType(ll_intp, shapety.count)

    ary = make_array(aryty)(context, builder, args[0])

    # We will change the target shape in this slot
    # (see normalize_reshape_value() below)
    newshape = cgutils.alloca_once(builder, ll_shape)
    builder.store(shape, newshape)

    # Create a shape array pointing to the value of newshape.
    # (roughly, `shape_ary = np.array(ary.shape)`)
    shape_ary_ty = types.Array(dtype=shapety.dtype, ndim=1, layout='C')
    shape_ary = make_array(shape_ary_ty)(context, builder)
    shape_itemsize = context.get_constant(types.intp,
                                          context.get_abi_sizeof(ll_intp))
    populate_array(shape_ary,
                   data=builder.bitcast(newshape, ll_intp.as_pointer()),
                   shape=[context.get_constant(types.intp, shapety.count)],
                   strides=[shape_itemsize],
                   itemsize=shape_itemsize,
                   meminfo=None)

    # Compute the original array size
    size = ary.nitems

    # Call our normalizer which will fix the shape array in case of negative
    # shape value
    context.compile_internal(builder, normalize_reshape_value,
                             typing.signature(types.void,
                                              types.uintp, shape_ary_ty),
                             [size, shape_ary._getvalue()])

    # Perform reshape (nocopy)
    newnd = shapety.count
    newstrides = cgutils.alloca_once(builder, ll_shape)

    ok = _attempt_nocopy_reshape(context, builder, aryty, ary, newnd,
                                 newshape, newstrides)
    fail = builder.icmp_unsigned('==', ok, ok.type(0))

    with builder.if_then(fail):
        msg = "incompatible shape for array"
        context.call_conv.return_user_exc(builder, NotImplementedError, (msg,))

    ret = make_array(retty)(context, builder)
    populate_array(ret,
                   data=ary.data,
                   shape=builder.load(newshape),
                   strides=builder.load(newstrides),
                   itemsize=ary.itemsize,
                   meminfo=ary.meminfo,
                   parent=ary.parent)
    res = ret._getvalue()
    return impl_ret_borrowed(context, builder, sig.return_type, res)


@lower_builtin('array.reshape', types.Array, types.VarArg(types.Any))
def array_reshape_vararg(context, builder, sig, args):
    new_sig, new_args = vararg_to_tuple(context, builder, sig, args)
    return array_reshape(context, builder, new_sig, new_args)


@overload(np.reshape)
def np_reshape(a, newshape):
    def np_reshape_impl(a, newshape):
        return a.reshape(newshape)
    return np_reshape_impl


@overload(np.resize)
def numpy_resize(a, new_shape):

    if not type_can_asarray(a):
        msg = 'The argument "a" must be array-like'
        raise errors.TypingError(msg)

    if not ((isinstance(new_shape, types.UniTuple)
             and
             isinstance(new_shape.dtype, types.Integer))
            or
            isinstance(new_shape, types.Integer)):
        msg = ('The argument "new_shape" must be an integer or '
               'a tuple of integers')
        raise errors.TypingError(msg)

    def impl(a, new_shape):
        a = np.asarray(a)
        a = np.ravel(a)

        if isinstance(new_shape, tuple):
            new_size = 1
            for dim_length in np.asarray(new_shape):
                new_size *= dim_length
                if dim_length < 0:
                    msg = 'All elements of `new_shape` must be non-negative'
                    raise ValueError(msg)
        else:
            if new_shape < 0:
                msg2 = 'All elements of `new_shape` must be non-negative'
                raise ValueError(msg2)
            new_size = new_shape

        if a.size == 0:
            return np.zeros(new_shape).astype(a.dtype)

        repeats = -(-new_size // a.size)  # ceil division
        res = a
        for i in range(repeats - 1):
            res = np.concatenate((res, a))
        res = res[:new_size]

        return np.reshape(res, new_shape)

    return impl


@overload(np.append)
def np_append(arr, values, axis=None):

    if not type_can_asarray(arr):
        raise errors.TypingError('The first argument "arr" must be array-like')

    if not type_can_asarray(values):
        raise errors.TypingError('The second argument "values" must be '
                                 'array-like')

    if is_nonelike(axis):
        def impl(arr, values, axis=None):
            arr = np.ravel(np.asarray(arr))
            values = np.ravel(np.asarray(values))
            return np.concatenate((arr, values))
    else:

        if not isinstance(axis, types.Integer):
            raise errors.TypingError('The third argument "axis" must be an '
                                     'integer')

        def impl(arr, values, axis=None):
            return np.concatenate((arr, values), axis=axis)
    return impl


@lower_builtin('array.ravel', types.Array)
def array_ravel(context, builder, sig, args):
    # Only support no argument version (default order='C')
    def imp_nocopy(ary):
        """No copy version"""
        return ary.reshape(ary.size)

    def imp_copy(ary):
        """Copy version"""
        return ary.flatten()

    # If the input array is C layout already, use the nocopy version
    if sig.args[0].layout == 'C':
        imp = imp_nocopy
    # otherwise, use flatten under-the-hood
    else:
        imp = imp_copy

    res = context.compile_internal(builder, imp, sig, args)
    res = impl_ret_new_ref(context, builder, sig.return_type, res)
    return res


@lower_builtin(np.ravel, types.Array)
def np_ravel(context, builder, sig, args):
    def np_ravel_impl(a):
        return a.ravel()

    return context.compile_internal(builder, np_ravel_impl, sig, args)


@lower_builtin('array.flatten', types.Array)
def array_flatten(context, builder, sig, args):
    # Only support flattening to C layout currently.
    def imp(ary):
        return ary.copy().reshape(ary.size)

    res = context.compile_internal(builder, imp, sig, args)
    res = impl_ret_new_ref(context, builder, sig.return_type, res)
    return res


@register_jitable
def _np_clip_impl(a, a_min, a_max, out):
    # Both a_min and a_max are numpy arrays
    ret = np.empty_like(a) if out is None else out
    a_b, a_min_b, a_max_b = np.broadcast_arrays(a, a_min, a_max)
    for index in np.ndindex(a_b.shape):
        val_a = a_b[index]
        val_a_min = a_min_b[index]
        val_a_max = a_max_b[index]
        ret[index] = min(max(val_a, val_a_min), val_a_max)

    return ret


@register_jitable
def _np_clip_impl_none(a, b, use_min, out):
    for index in np.ndindex(a.shape):
        val_a = a[index]
        val_b = b[index]
        if use_min:
            out[index] = min(val_a, val_b)
        else:
            out[index] = max(val_a, val_b)
    return out


@overload(np.clip)
def np_clip(a, a_min, a_max, out=None):
    if not type_can_asarray(a):
        raise errors.TypingError('The argument "a" must be array-like')

    if (not isinstance(a_min, types.NoneType) and
            not type_can_asarray(a_min)):
        raise errors.TypingError(('The argument "a_min" must be a number '
                                 'or an array-like'))

    if (not isinstance(a_max, types.NoneType) and
            not type_can_asarray(a_max)):
        raise errors.TypingError('The argument "a_max" must be a number '
                                 'or an array-like')

    if not (isinstance(out, types.Array) or is_nonelike(out)):
        msg = 'The argument "out" must be an array if it is provided'
        raise errors.TypingError(msg)

    # TODO: support scalar a (issue #3469)
    a_min_is_none = a_min is None or isinstance(a_min, types.NoneType)
    a_max_is_none = a_max is None or isinstance(a_max, types.NoneType)

    if a_min_is_none and a_max_is_none:
        # Raises value error when both a_min and a_max are None
        def np_clip_nn(a, a_min, a_max, out=None):
            raise ValueError("array_clip: must set either max or min")

        return np_clip_nn

    a_min_is_scalar = isinstance(a_min, types.Number)
    a_max_is_scalar = isinstance(a_max, types.Number)

    if a_min_is_scalar and a_max_is_scalar:
        def np_clip_ss(a, a_min, a_max, out=None):
            # a_min and a_max are scalars
            # since their shape will be empty
            # so broadcasting is not needed at all
            ret = np.empty_like(a) if out is None else out
            for index in np.ndindex(a.shape):
                val_a = a[index]
                ret[index] = min(max(val_a, a_min), a_max)

            return ret

        return np_clip_ss
    elif a_min_is_scalar and not a_max_is_scalar:
        if a_max_is_none:
            def np_clip_sn(a, a_min, a_max, out=None):
                # a_min is a scalar
                # since its shape will be empty
                # so broadcasting is not needed at all
                ret = np.empty_like(a) if out is None else out
                for index in np.ndindex(a.shape):
                    val_a = a[index]
                    ret[index] = max(val_a, a_min)

                return ret

            return np_clip_sn
        else:
            def np_clip_sa(a, a_min, a_max, out=None):
                # a_min is a scalar
                # since its shape will be empty
                # broadcast it to shape of a
                # by using np.full_like
                a_min_full = np.full_like(a, a_min)
                return _np_clip_impl(a, a_min_full, a_max, out)

            return np_clip_sa
    elif not a_min_is_scalar and a_max_is_scalar:
        if a_min_is_none:
            def np_clip_ns(a, a_min, a_max, out=None):
                # a_max is a scalar
                # since its shape will be empty
                # so broadcasting is not needed at all
                ret = np.empty_like(a) if out is None else out
                for index in np.ndindex(a.shape):
                    val_a = a[index]
                    ret[index] = min(val_a, a_max)

                return ret

            return np_clip_ns
        else:
            def np_clip_as(a, a_min, a_max, out=None):
                # a_max is a scalar
                # since its shape will be empty
                # broadcast it to shape of a
                # by using np.full_like
                a_max_full = np.full_like(a, a_max)
                return _np_clip_impl(a, a_min, a_max_full, out)

            return np_clip_as
    else:
        # Case where exactly one of a_min or a_max is None
        if a_min_is_none:
            def np_clip_na(a, a_min, a_max, out=None):
                # a_max is a numpy array but a_min is None
                ret = np.empty_like(a) if out is None else out
                a_b, a_max_b = np.broadcast_arrays(a, a_max)
                return _np_clip_impl_none(a_b, a_max_b, True, ret)

            return np_clip_na
        elif a_max_is_none:
            def np_clip_an(a, a_min, a_max, out=None):
                # a_min is a numpy array but a_max is None
                ret = np.empty_like(a) if out is None else out
                a_b, a_min_b = np.broadcast_arrays(a, a_min)
                return _np_clip_impl_none(a_b, a_min_b, False, ret)

            return np_clip_an
        else:
            def np_clip_aa(a, a_min, a_max, out=None):
                # Both a_min and a_max are clearly arrays
                # because none of the above branches
                # returned
                return _np_clip_impl(a, a_min, a_max, out)

            return np_clip_aa


@overload_method(types.Array, 'clip')
def array_clip(a, a_min=None, a_max=None, out=None):
    def impl(a, a_min=None, a_max=None, out=None):
        return np.clip(a, a_min, a_max, out)
    return impl


def _change_dtype(context, builder, oldty, newty, ary):
    """
    Attempt to fix up *ary* for switching from *oldty* to *newty*.

    See Numpy's array_descr_set()
    (np/core/src/multiarray/getset.c).
    Attempt to fix the array's shape and strides for a new dtype.
    False is returned on failure, True on success.
    """
    assert oldty.ndim == newty.ndim
    assert oldty.layout == newty.layout

    new_layout = ord(newty.layout)
    any_layout = ord('A')
    c_layout = ord('C')
    f_layout = ord('F')

    int8 = types.int8

    def imp(nd, dims, strides, old_itemsize, new_itemsize, layout):
        # Attempt to update the layout due to limitation of the numba
        # type system.
        if layout == any_layout:
            # Test rightmost stride to be contiguous
            if strides[-1] == old_itemsize:
                # Process this as if it is C contiguous
                layout = int8(c_layout)
            # Test leftmost stride to be F contiguous
            elif strides[0] == old_itemsize:
                # Process this as if it is F contiguous
                layout = int8(f_layout)

        if old_itemsize != new_itemsize and (layout == any_layout or nd == 0):
            return False

        if layout == c_layout:
            i = nd - 1
        else:
            i = 0

        if new_itemsize < old_itemsize:
            # If it is compatible, increase the size of the dimension
            # at the end (or at the front if F-contiguous)
            if (old_itemsize % new_itemsize) != 0:
                return False

            newdim = old_itemsize // new_itemsize
            dims[i] *= newdim
            strides[i] = new_itemsize

        elif new_itemsize > old_itemsize:
            # Determine if last (or first if F-contiguous) dimension
            # is compatible
            bytelength = dims[i] * old_itemsize
            if (bytelength % new_itemsize) != 0:
                return False

            dims[i] = bytelength // new_itemsize
            strides[i] = new_itemsize

        else:
            # Same item size: nothing to do (this also works for
            # non-contiguous arrays).
            pass

        return True

    old_itemsize = context.get_constant(types.intp,
                                        get_itemsize(context, oldty))
    new_itemsize = context.get_constant(types.intp,
                                        get_itemsize(context, newty))

    nd = context.get_constant(types.intp, newty.ndim)
    shape_data = cgutils.gep_inbounds(builder, ary._get_ptr_by_name('shape'),
                                      0, 0)
    strides_data = cgutils.gep_inbounds(builder,
                                        ary._get_ptr_by_name('strides'), 0, 0)

    shape_strides_array_type = types.Array(dtype=types.intp, ndim=1, layout='C')
    arycls = context.make_array(shape_strides_array_type)

    shape_constant = cgutils.pack_array(builder,
                                        [context.get_constant(types.intp,
                                                              newty.ndim)])

    sizeof_intp = context.get_abi_sizeof(context.get_data_type(types.intp))
    sizeof_intp = context.get_constant(types.intp, sizeof_intp)
    strides_constant = cgutils.pack_array(builder, [sizeof_intp])

    shape_ary = arycls(context, builder)

    populate_array(shape_ary,
                   data=shape_data,
                   shape=shape_constant,
                   strides=strides_constant,
                   itemsize=sizeof_intp,
                   meminfo=None)

    strides_ary = arycls(context, builder)
    populate_array(strides_ary,
                   data=strides_data,
                   shape=shape_constant,
                   strides=strides_constant,
                   itemsize=sizeof_intp,
                   meminfo=None)

    shape = shape_ary._getvalue()
    strides = strides_ary._getvalue()
    args = [nd, shape, strides, old_itemsize, new_itemsize,
            context.get_constant(types.int8, new_layout)]

    sig = signature(types.boolean,
                    types.intp,  # nd
                    shape_strides_array_type,  # dims
                    shape_strides_array_type,  # strides
                    types.intp,  # old_itemsize
                    types.intp,  # new_itemsize
                    types.int8,  # layout
                    )

    res = context.compile_internal(builder, imp, sig, args)
    update_array_info(newty, ary)
    res = impl_ret_borrowed(context, builder, sig.return_type, res)
    return res


@overload(np.shape)
def np_shape(a):
    if not type_can_asarray(a):
        raise errors.TypingError("The argument to np.shape must be array-like")

    def impl(a):
        return np.asarray(a).shape
    return impl


@overload(np.size)
def np_size(a):
    if not type_can_asarray(a):
        raise errors.TypingError("The argument to np.size must be array-like")

    def impl(a):
        return np.asarray(a).size
    return impl

# ------------------------------------------------------------------------------


@overload(np.unique)
def np_unique(ar):
    def np_unique_impl(ar):
        b = np.sort(ar.ravel())
        head = list(b[:1])
        tail = [x for i, x in enumerate(b[1:]) if b[i] != x]
        return np.array(head + tail)
    return np_unique_impl


@overload(np.repeat)
def np_repeat(a, repeats):
    # Implementation for repeats being a scalar is a module global function
    # (see below) because it might be called from the implementation below.

    def np_repeat_impl_repeats_array_like(a, repeats):
        # implementation if repeats is an array like
        repeats_array = np.asarray(repeats, dtype=np.int64)
        # if it is a singleton array, invoke the scalar implementation
        if repeats_array.shape[0] == 1:
            return np_repeat_impl_repeats_scaler(a, repeats_array[0])
        if np.any(repeats_array < 0):
            raise ValueError("negative dimensions are not allowed")
        asa = np.asarray(a)
        aravel = asa.ravel()
        n = aravel.shape[0]
        if aravel.shape != repeats_array.shape:
            raise ValueError(
                "operands could not be broadcast together")
        to_return = np.empty(np.sum(repeats_array), dtype=asa.dtype)
        pos = 0
        for i in range(n):
            to_return[pos : pos + repeats_array[i]] = aravel[i]
            pos += repeats_array[i]
        return to_return

    # type checking
    if isinstance(a, (types.Array,
                      types.List,
                      types.BaseTuple,
                      types.Number,
                      types.Boolean,
                      )
                  ):
        if isinstance(repeats, types.Integer):
            return np_repeat_impl_repeats_scaler
        elif isinstance(repeats, (types.Array, types.List)):
            if isinstance(repeats.dtype, types.Integer):
                return np_repeat_impl_repeats_array_like

        raise errors.TypingError(
            "The repeats argument must be an integer "
            "or an array-like of integer dtype")


@register_jitable
def np_repeat_impl_repeats_scaler(a, repeats):
    if repeats < 0:
        raise ValueError("negative dimensions are not allowed")
    asa = np.asarray(a)
    aravel = asa.ravel()
    n = aravel.shape[0]
    if repeats == 0:
        return np.empty(0, dtype=asa.dtype)
    elif repeats == 1:
        return np.copy(aravel)
    else:
        to_return = np.empty(n * repeats, dtype=asa.dtype)
        for i in range(n):
            to_return[i * repeats : (i + 1) * repeats] = aravel[i]
        return to_return


@extending.overload_method(types.Array, 'repeat')
def array_repeat(a, repeats):
    def array_repeat_impl(a, repeats):
        return np.repeat(a, repeats)

    return array_repeat_impl


@intrinsic
def _intrin_get_itemsize(tyctx, dtype):
    """Computes the itemsize of the dtype"""
    sig = types.intp(dtype)

    def codegen(cgctx, builder, sig, llargs):
        llty = cgctx.get_data_type(sig.args[0].dtype)
        llintp = cgctx.get_data_type(sig.return_type)
        return llintp(cgctx.get_abi_sizeof(llty))
    return sig, codegen


def _compatible_view(a, dtype):
    pass


@overload(_compatible_view, target='generic')
def ol_compatible_view(a, dtype):
    """Determines if the array and dtype are compatible for forming a view."""
    # NOTE: NumPy 1.23+ uses this check.
    # Code based on:
    # https://github.com/numpy/numpy/blob/750ad21258cfc00663586d5a466e24f91b48edc7/numpy/core/src/multiarray/getset.c#L500-L555  # noqa: E501
    def impl(a, dtype):
        dtype_size = _intrin_get_itemsize(dtype)
        if dtype_size != a.itemsize:
            # catch forbidden cases
            if a.ndim == 0:
                msg1 = ("Changing the dtype of a 0d array is only supported "
                        "if the itemsize is unchanged")
                raise ValueError(msg1)
            else:
                # NumPy has a check here for subarray type conversion which
                # Numba doesn't support
                pass

            # Resize on last axis only
            axis = a.ndim - 1
            p1 = a.shape[axis] != 1
            p2 = a.size != 0
            p3 = a.strides[axis] != a.itemsize
            if (p1 and p2 and p3):
                msg2 = ("To change to a dtype of a different size, the last "
                        "axis must be contiguous")
                raise ValueError(msg2)

            if dtype_size < a.itemsize:
                if dtype_size == 0 or a.itemsize % dtype_size != 0:
                    msg3 = ("When changing to a smaller dtype, its size must "
                            "be a divisor of the size of original dtype")
                    raise ValueError(msg3)
            else:
                newdim = a.shape[axis] * a.itemsize
                if newdim % dtype_size != 0:
                    msg4 = ("When changing to a larger dtype, its size must be "
                            "a divisor of the total size in bytes of the last "
                            "axis of the array.")
                    raise ValueError(msg4)
    return impl


@lower_builtin('array.view', types.Array, types.DTypeSpec)
def array_view(context, builder, sig, args):
    aryty = sig.args[0]
    retty = sig.return_type

    ary = make_array(aryty)(context, builder, args[0])
    ret = make_array(retty)(context, builder)
    # Copy all fields, casting the "data" pointer appropriately
    fields = set(ret._datamodel._fields)
    for k in sorted(fields):
        val = getattr(ary, k)
        if k == 'data':
            ptrty = ret.data.type
            ret.data = builder.bitcast(val, ptrty)
        else:
            setattr(ret, k, val)

    if numpy_version >= (1, 23):
        # NumPy 1.23+ bans views using a dtype that is a different size to that
        # of the array when the last axis is not contiguous. For example, this
        # manifests at runtime when a dtype size altering view is requested
        # on a Fortran ordered array.

        tyctx = context.typing_context
        fnty = tyctx.resolve_value_type(_compatible_view)
        _compatible_view_sig = fnty.get_call_type(tyctx, (*sig.args,), {})
        impl = context.get_function(fnty, _compatible_view_sig)
        impl(builder, args)

    ok = _change_dtype(context, builder, aryty, retty, ret)
    fail = builder.icmp_unsigned('==', ok, Constant(ok.type, 0))

    with builder.if_then(fail):
        msg = "new type not compatible with array"
        context.call_conv.return_user_exc(builder, ValueError, (msg,))

    res = ret._getvalue()
    return impl_ret_borrowed(context, builder, sig.return_type, res)


# ------------------------------------------------------------------------------
# Array attributes

@lower_getattr(types.Array, "dtype")
def array_dtype(context, builder, typ, value):
    res = context.get_dummy_value()
    return impl_ret_untracked(context, builder, typ, res)


@lower_getattr(types.Array, "shape")
@lower_getattr(types.MemoryView, "shape")
def array_shape(context, builder, typ, value):
    arrayty = make_array(typ)
    array = arrayty(context, builder, value)
    res = array.shape
    return impl_ret_untracked(context, builder, typ, res)


@lower_getattr(types.Array, "strides")
@lower_getattr(types.MemoryView, "strides")
def array_strides(context, builder, typ, value):
    arrayty = make_array(typ)
    array = arrayty(context, builder, value)
    res = array.strides
    return impl_ret_untracked(context, builder, typ, res)


@lower_getattr(types.Array, "ndim")
@lower_getattr(types.MemoryView, "ndim")
def array_ndim(context, builder, typ, value):
    res = context.get_constant(types.intp, typ.ndim)
    return impl_ret_untracked(context, builder, typ, res)


@lower_getattr(types.Array, "size")
def array_size(context, builder, typ, value):
    arrayty = make_array(typ)
    array = arrayty(context, builder, value)
    res = array.nitems
    return impl_ret_untracked(context, builder, typ, res)


@lower_getattr(types.Array, "itemsize")
@lower_getattr(types.MemoryView, "itemsize")
def array_itemsize(context, builder, typ, value):
    arrayty = make_array(typ)
    array = arrayty(context, builder, value)
    res = array.itemsize
    return impl_ret_untracked(context, builder, typ, res)


@lower_getattr(types.Array, "nbytes")
@lower_getattr(types.MemoryView, "nbytes")
def array_nbytes(context, builder, typ, value):
    """
    nbytes = size * itemsize
    """
    arrayty = make_array(typ)
    array = arrayty(context, builder, value)
    res = builder.mul(array.nitems, array.itemsize)
    return impl_ret_untracked(context, builder, typ, res)


@lower_getattr(types.MemoryView, "contiguous")
def array_contiguous(context, builder, typ, value):
    res = context.get_constant(types.boolean, typ.is_contig)
    return impl_ret_untracked(context, builder, typ, res)


@lower_getattr(types.MemoryView, "c_contiguous")
def array_c_contiguous(context, builder, typ, value):
    res = context.get_constant(types.boolean, typ.is_c_contig)
    return impl_ret_untracked(context, builder, typ, res)


@lower_getattr(types.MemoryView, "f_contiguous")
def array_f_contiguous(context, builder, typ, value):
    res = context.get_constant(types.boolean, typ.is_f_contig)
    return impl_ret_untracked(context, builder, typ, res)


@lower_getattr(types.MemoryView, "readonly")
def array_readonly(context, builder, typ, value):
    res = context.get_constant(types.boolean, not typ.mutable)
    return impl_ret_untracked(context, builder, typ, res)


# array.ctypes

@lower_getattr(types.Array, "ctypes")
def array_ctypes(context, builder, typ, value):
    arrayty = make_array(typ)
    array = arrayty(context, builder, value)
    # Create new ArrayCType structure
    act = types.ArrayCTypes(typ)
    ctinfo = context.make_helper(builder, act)
    ctinfo.data = array.data
    ctinfo.meminfo = array.meminfo
    res = ctinfo._getvalue()
    return impl_ret_borrowed(context, builder, act, res)


@lower_getattr(types.ArrayCTypes, "data")
def array_ctypes_data(context, builder, typ, value):
    ctinfo = context.make_helper(builder, typ, value=value)
    res = ctinfo.data
    # Convert it to an integer
    res = builder.ptrtoint(res, context.get_value_type(types.intp))
    return impl_ret_untracked(context, builder, typ, res)


@lower_cast(types.ArrayCTypes, types.CPointer)
@lower_cast(types.ArrayCTypes, types.voidptr)
def array_ctypes_to_pointer(context, builder, fromty, toty, val):
    ctinfo = context.make_helper(builder, fromty, value=val)
    res = ctinfo.data
    res = builder.bitcast(res, context.get_value_type(toty))
    return impl_ret_untracked(context, builder, toty, res)


def _call_contiguous_check(checker, context, builder, aryty, ary):
    """Helper to invoke the contiguous checker function on an array

    Args
    ----
    checker :
        ``numba.numpy_supports.is_contiguous``, or
        ``numba.numpy_supports.is_fortran``.
    context : target context
    builder : llvm ir builder
    aryty : numba type
    ary : llvm value
    """
    ary = make_array(aryty)(context, builder, value=ary)
    tup_intp = types.UniTuple(types.intp, aryty.ndim)
    itemsize = context.get_abi_sizeof(context.get_value_type(aryty.dtype))
    check_sig = signature(types.bool_, tup_intp, tup_intp, types.intp)
    check_args = [ary.shape, ary.strides,
                  context.get_constant(types.intp, itemsize)]
    is_contig = context.compile_internal(builder, checker, check_sig,
                                         check_args)
    return is_contig


# array.flags

@lower_getattr(types.Array, "flags")
def array_flags(context, builder, typ, value):
    flagsobj = context.make_helper(builder, types.ArrayFlags(typ))
    flagsobj.parent = value
    res = flagsobj._getvalue()
    context.nrt.incref(builder, typ, value)
    return impl_ret_new_ref(context, builder, typ, res)


@lower_getattr(types.ArrayFlags, "contiguous")
@lower_getattr(types.ArrayFlags, "c_contiguous")
def array_flags_c_contiguous(context, builder, typ, value):
    if typ.array_type.layout != 'C':
        # any layout can still be contiguous
        flagsobj = context.make_helper(builder, typ, value=value)
        res = _call_contiguous_check(is_contiguous, context, builder,
                                     typ.array_type, flagsobj.parent)
    else:
        val = typ.array_type.layout == 'C'
        res = context.get_constant(types.boolean, val)
    return impl_ret_untracked(context, builder, typ, res)


@lower_getattr(types.ArrayFlags, "f_contiguous")
def array_flags_f_contiguous(context, builder, typ, value):
    if typ.array_type.layout != 'F':
        # any layout can still be contiguous
        flagsobj = context.make_helper(builder, typ, value=value)
        res = _call_contiguous_check(is_fortran, context, builder,
                                     typ.array_type, flagsobj.parent)
    else:
        layout = typ.array_type.layout
        val = layout == 'F' if typ.array_type.ndim > 1 else layout in 'CF'
        res = context.get_constant(types.boolean, val)
    return impl_ret_untracked(context, builder, typ, res)


# ------------------------------------------------------------------------------
# .real / .imag

@lower_getattr(types.Array, "real")
def array_real_part(context, builder, typ, value):
    if typ.dtype in types.complex_domain:
        return array_complex_attr(context, builder, typ, value, attr='real')
    elif typ.dtype in types.number_domain:
        # as an identity function
        return impl_ret_borrowed(context, builder, typ, value)
    else:
        raise NotImplementedError('unsupported .real for {}'.format(type.dtype))


@lower_getattr(types.Array, "imag")
def array_imag_part(context, builder, typ, value):
    if typ.dtype in types.complex_domain:
        return array_complex_attr(context, builder, typ, value, attr='imag')
    elif typ.dtype in types.number_domain:
        # return a readonly zero array
        sig = signature(typ.copy(readonly=True), typ)
        arrtype, shapes = _parse_empty_like_args(context, builder, sig, [value])
        ary = _empty_nd_impl(context, builder, arrtype, shapes)
        cgutils.memset(builder, ary.data, builder.mul(ary.itemsize,
                                                      ary.nitems), 0)
        return impl_ret_new_ref(context, builder, sig.return_type,
                                ary._getvalue())
    else:
        raise NotImplementedError('unsupported .imag for {}'.format(type.dtype))


def array_complex_attr(context, builder, typ, value, attr):
    """
    Given a complex array, it's memory layout is:

        R C R C R C
        ^   ^   ^

    (`R` indicates a float for the real part;
     `C` indicates a float for the imaginary part;
     the `^` indicates the start of each element)

    To get the real part, we can simply change the dtype and itemsize to that
    of the underlying float type.  The new layout is:

        R x R x R x
        ^   ^   ^

    (`x` indicates unused)

    A load operation will use the dtype to determine the number of bytes to
    load.

    To get the imaginary part, we shift the pointer by 1 float offset and
    change the dtype and itemsize.  The new layout is:

        x C x C x C
          ^   ^   ^
    """
    if attr not in ['real', 'imag'] or typ.dtype not in types.complex_domain:
        raise NotImplementedError("cannot get attribute `{}`".format(attr))

    arrayty = make_array(typ)
    array = arrayty(context, builder, value)

    # sizeof underlying float type
    flty = typ.dtype.underlying_float
    sizeof_flty = context.get_abi_sizeof(context.get_data_type(flty))
    itemsize = array.itemsize.type(sizeof_flty)

    # cast data pointer to float type
    llfltptrty = context.get_value_type(flty).as_pointer()
    dataptr = builder.bitcast(array.data, llfltptrty)

    # add offset
    if attr == 'imag':
        dataptr = builder.gep(dataptr, [ir.IntType(32)(1)])

    # make result
    resultty = typ.copy(dtype=flty, layout='A')
    result = make_array(resultty)(context, builder)
    repl = dict(data=dataptr, itemsize=itemsize)
    cgutils.copy_struct(result, array, repl)
    return impl_ret_borrowed(context, builder, resultty, result._getvalue())


@overload_method(types.Array, 'conj')
@overload_method(types.Array, 'conjugate')
def array_conj(arr):
    def impl(arr):
        return np.conj(arr)
    return impl

# ------------------------------------------------------------------------------
# DType attribute


def dtype_type(context, builder, dtypety, dtypeval):
    # Just return a dummy opaque value
    return context.get_dummy_value()


lower_getattr(types.DType, 'type')(dtype_type)
lower_getattr(types.DType, 'kind')(dtype_type)


# ------------------------------------------------------------------------------
# static_getitem on Numba numerical types to create "array" types


@lower_builtin('static_getitem', types.NumberClass, types.Any)
def static_getitem_number_clazz(context, builder, sig, args):
    """This handles the "static_getitem" when a Numba type is subscripted e.g:
    var = typed.List.empty_list(float64[::1, :])
    It only allows this on simple numerical types. Compound types, like
    records, are not supported.
    """
    retty = sig.return_type
    if isinstance(retty, types.Array):
        # This isn't used or practically accessible, but has to exist, so just
        # put in a NULL of the right type.
        res = context.get_value_type(retty)(None)
        return impl_ret_untracked(context, builder, retty, res)
    else:
        # This should be unreachable unless the implementation on the Type
        # metaclass is changed.
        msg = ("Unreachable; the definition of __getitem__ on the "
               "numba.types.abstract.Type metaclass should prevent access.")
        raise errors.LoweringError(msg)


# ------------------------------------------------------------------------------
# Structured / record lookup

@lower_getattr_generic(types.Array)
def array_record_getattr(context, builder, typ, value, attr):
    """
    Generic getattr() implementation for record arrays: fetch the given
    record member, i.e. a subarray.
    """
    arrayty = make_array(typ)
    array = arrayty(context, builder, value)

    rectype = typ.dtype
    if not isinstance(rectype, types.Record):
        raise NotImplementedError("attribute %r of %s not defined"
                                  % (attr, typ))
    dtype = rectype.typeof(attr)
    offset = rectype.offset(attr)

    if isinstance(dtype, types.NestedArray):
        resty = typ.copy(
            dtype=dtype.dtype, ndim=typ.ndim + dtype.ndim, layout='A')
    else:
        resty = typ.copy(dtype=dtype, layout='A')

    raryty = make_array(resty)

    rary = raryty(context, builder)

    constoffset = context.get_constant(types.intp, offset)

    newdataptr = cgutils.pointer_add(
        builder, array.data, constoffset,  return_type=rary.data.type,
    )
    if isinstance(dtype, types.NestedArray):
        # new shape = recarray shape + inner dimension from nestedarray
        shape = cgutils.unpack_tuple(builder, array.shape, typ.ndim)
        shape += [context.get_constant(types.intp, i) for i in dtype.shape]
        # new strides = recarray strides + strides of the inner nestedarray
        strides = cgutils.unpack_tuple(builder, array.strides, typ.ndim)
        strides += [context.get_constant(types.intp, i) for i in dtype.strides]
        # New datasize = size of elements of the nestedarray
        datasize = context.get_abi_sizeof(context.get_data_type(dtype.dtype))
    else:
        # New shape, strides, and datasize match the underlying array
        shape = array.shape
        strides = array.strides
        datasize = context.get_abi_sizeof(context.get_data_type(dtype))
    populate_array(rary,
                   data=newdataptr,
                   shape=shape,
                   strides=strides,
                   itemsize=context.get_constant(types.intp, datasize),
                   meminfo=array.meminfo,
                   parent=array.parent)
    res = rary._getvalue()
    return impl_ret_borrowed(context, builder, resty, res)


@lower_builtin('static_getitem', types.Array, types.StringLiteral)
def array_record_getitem(context, builder, sig, args):
    index = args[1]
    if not isinstance(index, str):
        # This will fallback to normal getitem
        raise NotImplementedError
    return array_record_getattr(context, builder, sig.args[0], args[0], index)


@lower_getattr_generic(types.Record)
def record_getattr(context, builder, typ, value, attr):
    """
    Generic getattr() implementation for records: get the given record member.
    """
    context.sentry_record_alignment(typ, attr)
    offset = typ.offset(attr)
    elemty = typ.typeof(attr)

    if isinstance(elemty, types.NestedArray):
        # Only a nested array's *data* is stored in a structured array,
        # so we create an array structure to point to that data.
        aryty = make_array(elemty)
        ary = aryty(context, builder)
        dtype = elemty.dtype
        newshape = [context.get_constant(types.intp, s) for s in
                    elemty.shape]
        newstrides = [context.get_constant(types.intp, s) for s in
                      elemty.strides]
        newdata = cgutils.get_record_member(builder, value, offset,
                                            context.get_data_type(dtype))
        populate_array(
            ary,
            data=newdata,
            shape=cgutils.pack_array(builder, newshape),
            strides=cgutils.pack_array(builder, newstrides),
            itemsize=context.get_constant(types.intp, elemty.size),
            meminfo=None,
            parent=None,
        )
        res = ary._getvalue()
        return impl_ret_borrowed(context, builder, typ, res)
    else:
        dptr = cgutils.get_record_member(builder, value, offset,
                                         context.get_data_type(elemty))
        align = None if typ.aligned else 1
        res = context.unpack_value(builder, elemty, dptr, align)
        return impl_ret_borrowed(context, builder, typ, res)


@lower_setattr_generic(types.Record)
def record_setattr(context, builder, sig, args, attr):
    """
    Generic setattr() implementation for records: set the given record member.
    """
    typ, valty = sig.args
    target, val = args

    context.sentry_record_alignment(typ, attr)
    offset = typ.offset(attr)
    elemty = typ.typeof(attr)

    if isinstance(elemty, types.NestedArray):
        # Copy the data from the RHS into the nested array
        val_struct = cgutils.create_struct_proxy(valty)(context, builder,
                                                        value=args[1])
        src = val_struct.data
        dest = cgutils.get_record_member(builder, target, offset,
                                         src.type.pointee)
        cgutils.memcpy(builder, dest, src,
                       context.get_constant(types.intp, elemty.nitems))
    else:
        # Set the given scalar record member
        dptr = cgutils.get_record_member(builder, target, offset,
                                         context.get_data_type(elemty))
        val = context.cast(builder, val, valty, elemty)
        align = None if typ.aligned else 1
        context.pack_value(builder, elemty, val, dptr, align=align)


@lower_builtin('static_getitem', types.Record, types.StringLiteral)
def record_static_getitem_str(context, builder, sig, args):
    """
    Record.__getitem__ redirects to getattr()
    """
    impl = context.get_getattr(sig.args[0], args[1])
    return impl(context, builder, sig.args[0], args[0], args[1])


@lower_builtin('static_getitem', types.Record, types.IntegerLiteral)
def record_static_getitem_int(context, builder, sig, args):
    """
    Record.__getitem__ redirects to getattr()
    """
    idx = sig.args[1].literal_value
    fields = list(sig.args[0].fields)
    ll_field = context.insert_const_string(builder.module, fields[idx])
    impl = context.get_getattr(sig.args[0], ll_field)
    return impl(context, builder, sig.args[0], args[0], fields[idx])


@lower_builtin('static_setitem', types.Record, types.StringLiteral, types.Any)
def record_static_setitem_str(context, builder, sig, args):
    """
    Record.__setitem__ redirects to setattr()
    """
    recty, _, valty = sig.args
    rec, idx, val = args
    getattr_sig = signature(sig.return_type, recty, valty)
    impl = context.get_setattr(idx, getattr_sig)
    assert impl is not None
    return impl(builder, (rec, val))


@lower_builtin('static_setitem', types.Record, types.IntegerLiteral, types.Any)
def record_static_setitem_int(context, builder, sig, args):
    """
    Record.__setitem__ redirects to setattr()
    """
    recty, _, valty = sig.args
    rec, idx, val = args
    getattr_sig = signature(sig.return_type, recty, valty)
    fields = list(sig.args[0].fields)
    impl = context.get_setattr(fields[idx], getattr_sig)
    assert impl is not None
    return impl(builder, (rec, val))


# ------------------------------------------------------------------------------
# Constant arrays and records


@lower_constant(types.Array)
def constant_array(context, builder, ty, pyval):
    """
    Create a constant array (mechanism is target-dependent).
    """
    return context.make_constant_array(builder, ty, pyval)


@lower_constant(types.Record)
def constant_record(context, builder, ty, pyval):
    """
    Create a record constant as a stack-allocated array of bytes.
    """
    lty = ir.ArrayType(ir.IntType(8), pyval.nbytes)
    val = lty(bytearray(pyval.tostring()))
    return cgutils.alloca_once_value(builder, val)


@lower_constant(types.Bytes)
def constant_bytes(context, builder, ty, pyval):
    """
    Create a constant array from bytes (mechanism is target-dependent).
    """
    buf = np.array(bytearray(pyval), dtype=np.uint8)
    return context.make_constant_array(builder, ty, buf)

# ------------------------------------------------------------------------------
# Comparisons


@lower_builtin(operator.is_, types.Array, types.Array)
def array_is(context, builder, sig, args):
    aty, bty = sig.args
    if aty != bty:
        return cgutils.false_bit

    def array_is_impl(a, b):
        return (a.shape == b.shape and
                a.strides == b.strides and
                a.ctypes.data == b.ctypes.data)

    return context.compile_internal(builder, array_is_impl, sig, args)

# ------------------------------------------------------------------------------
# Hash


@overload_attribute(types.Array, "__hash__")
def ol_array_hash(arr):
    return lambda arr: None


# ------------------------------------------------------------------------------
# builtin `np.flat` implementation

def make_array_flat_cls(flatiterty):
    """
    Return the Structure representation of the given *flatiterty* (an
    instance of types.NumpyFlatType).
    """
    return _make_flattening_iter_cls(flatiterty, 'flat')


def make_array_ndenumerate_cls(nditerty):
    """
    Return the Structure representation of the given *nditerty* (an
    instance of types.NumpyNdEnumerateType).
    """
    return _make_flattening_iter_cls(nditerty, 'ndenumerate')


def _increment_indices(context, builder, ndim, shape, indices, end_flag=None,
                       loop_continue=None, loop_break=None):
    zero = context.get_constant(types.intp, 0)

    bbend = builder.append_basic_block('end_increment')

    if end_flag is not None:
        builder.store(cgutils.false_byte, end_flag)

    for dim in reversed(range(ndim)):
        idxptr = cgutils.gep_inbounds(builder, indices, dim)
        idx = cgutils.increment_index(builder, builder.load(idxptr))

        count = shape[dim]
        in_bounds = builder.icmp_signed('<', idx, count)
        with cgutils.if_likely(builder, in_bounds):
            # New index is still in bounds
            builder.store(idx, idxptr)
            if loop_continue is not None:
                loop_continue(dim)
            builder.branch(bbend)
        # Index out of bounds => reset it and proceed it to outer index
        builder.store(zero, idxptr)
        if loop_break is not None:
            loop_break(dim)

    if end_flag is not None:
        builder.store(cgutils.true_byte, end_flag)
    builder.branch(bbend)

    builder.position_at_end(bbend)


def _increment_indices_array(context, builder, arrty, arr, indices,
                             end_flag=None):
    shape = cgutils.unpack_tuple(builder, arr.shape, arrty.ndim)
    _increment_indices(context, builder, arrty.ndim, shape, indices, end_flag)


def make_nditer_cls(nditerty):
    """
    Return the Structure representation of the given *nditerty* (an
    instance of types.NumpyNdIterType).
    """
    ndim = nditerty.ndim
    layout = nditerty.layout
    narrays = len(nditerty.arrays)
    nshapes = ndim if nditerty.need_shaped_indexing else 1

    class BaseSubIter(object):
        """
        Base class for sub-iterators of a nditer() instance.
        """

        def __init__(self, nditer, member_name, start_dim, end_dim):
            self.nditer = nditer
            self.member_name = member_name
            self.start_dim = start_dim
            self.end_dim = end_dim
            self.ndim = end_dim - start_dim

        def set_member_ptr(self, ptr):
            setattr(self.nditer, self.member_name, ptr)

        @functools.cached_property
        def member_ptr(self):
            return getattr(self.nditer, self.member_name)

        def init_specific(self, context, builder):
            pass

        def loop_continue(self, context, builder, logical_dim):
            pass

        def loop_break(self, context, builder, logical_dim):
            pass

    class FlatSubIter(BaseSubIter):
        """
        Sub-iterator walking a contiguous array in physical order, with
        support for broadcasting (the index is reset on the outer dimension).
        """

        def init_specific(self, context, builder):
            zero = context.get_constant(types.intp, 0)
            self.set_member_ptr(cgutils.alloca_once_value(builder, zero))

        def compute_pointer(self, context, builder, indices, arrty, arr):
            index = builder.load(self.member_ptr)
            return builder.gep(arr.data, [index])

        def loop_continue(self, context, builder, logical_dim):
            if logical_dim == self.ndim - 1:
                # Only increment index inside innermost logical dimension
                index = builder.load(self.member_ptr)
                index = cgutils.increment_index(builder, index)
                builder.store(index, self.member_ptr)

        def loop_break(self, context, builder, logical_dim):
            if logical_dim == 0:
                # At the exit of outermost logical dimension, reset index
                zero = context.get_constant(types.intp, 0)
                builder.store(zero, self.member_ptr)
            elif logical_dim == self.ndim - 1:
                # Inside innermost logical dimension, increment index
                index = builder.load(self.member_ptr)
                index = cgutils.increment_index(builder, index)
                builder.store(index, self.member_ptr)

    class TrivialFlatSubIter(BaseSubIter):
        """
        Sub-iterator walking a contiguous array in physical order,
        *without* support for broadcasting.
        """

        def init_specific(self, context, builder):
            assert not nditerty.need_shaped_indexing

        def compute_pointer(self, context, builder, indices, arrty, arr):
            assert len(indices) <= 1, len(indices)
            return builder.gep(arr.data, indices)

    class IndexedSubIter(BaseSubIter):
        """
        Sub-iterator walking an array in logical order.
        """

        def compute_pointer(self, context, builder, indices, arrty, arr):
            assert len(indices) == self.ndim
            return cgutils.get_item_pointer(context, builder, arrty, arr,
                                            indices, wraparound=False)

    class ZeroDimSubIter(BaseSubIter):
        """
        Sub-iterator "walking" a 0-d array.
        """

        def compute_pointer(self, context, builder, indices, arrty, arr):
            return arr.data

    class ScalarSubIter(BaseSubIter):
        """
        Sub-iterator "walking" a scalar value.
        """

        def compute_pointer(self, context, builder, indices, arrty, arr):
            return arr

    class NdIter(cgutils.create_struct_proxy(nditerty)):
        """
        .nditer() implementation.

        Note: 'F' layout means the shape is iterated in reverse logical order,
        so indices and shapes arrays have to be reversed as well.
        """

        @functools.cached_property
        def subiters(self):
            l = []
            factories = {'flat': FlatSubIter if nditerty.need_shaped_indexing
                         else TrivialFlatSubIter,
                         'indexed': IndexedSubIter,
                         '0d': ZeroDimSubIter,
                         'scalar': ScalarSubIter,
                         }
            for i, sub in enumerate(nditerty.indexers):
                kind, start_dim, end_dim, _ = sub
                member_name = 'index%d' % i
                factory = factories[kind]
                l.append(factory(self, member_name, start_dim, end_dim))
            return l

        def init_specific(self, context, builder, arrtys, arrays):
            """
            Initialize the nditer() instance for the specific array inputs.
            """
            zero = context.get_constant(types.intp, 0)

            # Store inputs
            self.arrays = context.make_tuple(builder, types.Tuple(arrtys),
                                             arrays)
            # Create slots for scalars
            for i, ty in enumerate(arrtys):
                if not isinstance(ty, types.Array):
                    member_name = 'scalar%d' % i
                    # XXX as_data()?
                    slot = cgutils.alloca_once_value(builder, arrays[i])
                    setattr(self, member_name, slot)

            arrays = self._arrays_or_scalars(context, builder, arrtys, arrays)

            # Extract iterator shape (the shape of the most-dimensional input)
            main_shape_ty = types.UniTuple(types.intp, ndim)
            main_shape = None
            main_nitems = None
            for i, arrty in enumerate(arrtys):
                if isinstance(arrty, types.Array) and arrty.ndim == ndim:
                    main_shape = arrays[i].shape
                    main_nitems = arrays[i].nitems
                    break
            else:
                # Only scalar inputs => synthesize a dummy shape
                assert ndim == 0
                main_shape = context.make_tuple(builder, main_shape_ty, ())
                main_nitems = context.get_constant(types.intp, 1)

            # Validate shapes of array inputs
            def check_shape(shape, main_shape):
                n = len(shape)
                for i in range(n):
                    if shape[i] != main_shape[len(main_shape) - n + i]:
                        raise ValueError("nditer(): operands could not be "
                                         "broadcast together")

            for arrty, arr in zip(arrtys, arrays):
                if isinstance(arrty, types.Array) and arrty.ndim > 0:
                    sig = signature(types.none,
                                    types.UniTuple(types.intp, arrty.ndim),
                                    main_shape_ty)
                    context.compile_internal(builder, check_shape,
                                             sig, (arr.shape, main_shape))

            # Compute shape and size
            shapes = cgutils.unpack_tuple(builder, main_shape)
            if layout == 'F':
                shapes = shapes[::-1]

            # If shape is empty, mark iterator exhausted
            shape_is_empty = builder.icmp_signed('==', main_nitems, zero)
            exhausted = builder.select(shape_is_empty, cgutils.true_byte,
                                       cgutils.false_byte)

            if not nditerty.need_shaped_indexing:
                # Flatten shape to make iteration faster on small innermost
                # dimensions (e.g. a (100000, 3) shape)
                shapes = (main_nitems,)
            assert len(shapes) == nshapes

            indices = cgutils.alloca_once(builder, zero.type, size=nshapes)
            for dim in range(nshapes):
                idxptr = cgutils.gep_inbounds(builder, indices, dim)
                builder.store(zero, idxptr)

            self.indices = indices
            self.shape = cgutils.pack_array(builder, shapes, zero.type)
            self.exhausted = cgutils.alloca_once_value(builder, exhausted)

            # Initialize subiterators
            for subiter in self.subiters:
                subiter.init_specific(context, builder)

        def iternext_specific(self, context, builder, result):
            """
            Compute next iteration of the nditer() instance.
            """
            bbend = builder.append_basic_block('end')

            # Branch early if exhausted
            exhausted = cgutils.as_bool_bit(builder,
                                            builder.load(self.exhausted))
            with cgutils.if_unlikely(builder, exhausted):
                result.set_valid(False)
                builder.branch(bbend)

            arrtys = nditerty.arrays
            arrays = cgutils.unpack_tuple(builder, self.arrays)
            arrays = self._arrays_or_scalars(context, builder, arrtys, arrays)
            indices = self.indices

            # Compute iterated results
            result.set_valid(True)
            views = self._make_views(context, builder, indices, arrtys, arrays)
            views = [v._getvalue() for v in views]
            if len(views) == 1:
                result.yield_(views[0])
            else:
                result.yield_(context.make_tuple(builder, nditerty.yield_type,
                                                 views))

            shape = cgutils.unpack_tuple(builder, self.shape)
            _increment_indices(context, builder, len(shape), shape,
                               indices, self.exhausted,
                               functools.partial(self._loop_continue,
                                                 context,
                                                 builder),
                               functools.partial(self._loop_break,
                                                 context,
                                                 builder),
                               )

            builder.branch(bbend)
            builder.position_at_end(bbend)

        def _loop_continue(self, context, builder, dim):
            for sub in self.subiters:
                if sub.start_dim <= dim < sub.end_dim:
                    sub.loop_continue(context, builder, dim - sub.start_dim)

        def _loop_break(self, context, builder, dim):
            for sub in self.subiters:
                if sub.start_dim <= dim < sub.end_dim:
                    sub.loop_break(context, builder, dim - sub.start_dim)

        def _make_views(self, context, builder, indices, arrtys, arrays):
            """
            Compute the views to be yielded.
            """
            views = [None] * narrays
            indexers = nditerty.indexers
            subiters = self.subiters
            rettys = nditerty.yield_type
            if isinstance(rettys, types.BaseTuple):
                rettys = list(rettys)
            else:
                rettys = [rettys]
            indices = [builder.load(cgutils.gep_inbounds(builder, indices, i))
                       for i in range(nshapes)]

            for sub, subiter in zip(indexers, subiters):
                _, _, _, array_indices = sub
                sub_indices = indices[subiter.start_dim:subiter.end_dim]
                if layout == 'F':
                    sub_indices = sub_indices[::-1]
                for i in array_indices:
                    assert views[i] is None
                    views[i] = self._make_view(context, builder, sub_indices,
                                               rettys[i],
                                               arrtys[i], arrays[i], subiter)
            assert all(v for v in views)
            return views

        def _make_view(self, context, builder, indices, retty, arrty, arr,
                       subiter):
            """
            Compute a 0d view for a given input array.
            """
            assert isinstance(retty, types.Array) and retty.ndim == 0

            ptr = subiter.compute_pointer(context, builder, indices, arrty, arr)
            view = context.make_array(retty)(context, builder)

            itemsize = get_itemsize(context, retty)
            shape = context.make_tuple(builder, types.UniTuple(types.intp, 0),
                                       ())
            strides = context.make_tuple(builder, types.UniTuple(types.intp, 0),
                                         ())
            # HACK: meminfo=None avoids expensive refcounting operations
            # on ephemeral views
            populate_array(view, ptr, shape, strides, itemsize, meminfo=None)
            return view

        def _arrays_or_scalars(self, context, builder, arrtys, arrays):
            # Return a list of either array structures or pointers to
            # scalar slots
            l = []
            for i, (arrty, arr) in enumerate(zip(arrtys, arrays)):
                if isinstance(arrty, types.Array):
                    l.append(context.make_array(arrty)(context,
                                                       builder,
                                                       value=arr))
                else:
                    l.append(getattr(self, "scalar%d" % i))
            return l

    return NdIter


def make_ndindex_cls(nditerty):
    """
    Return the Structure representation of the given *nditerty* (an
    instance of types.NumpyNdIndexType).
    """
    ndim = nditerty.ndim

    class NdIndexIter(cgutils.create_struct_proxy(nditerty)):
        """
        .ndindex() implementation.
        """

        def init_specific(self, context, builder, shapes):
            zero = context.get_constant(types.intp, 0)
            indices = cgutils.alloca_once(builder, zero.type,
                                          size=context.get_constant(types.intp,
                                                                    ndim))
            exhausted = cgutils.alloca_once_value(builder, cgutils.false_byte)

            for dim in range(ndim):
                idxptr = cgutils.gep_inbounds(builder, indices, dim)
                builder.store(zero, idxptr)
                # 0-sized dimensions really indicate an empty array,
                # but we have to catch that condition early to avoid
                # a bug inside the iteration logic.
                dim_size = shapes[dim]
                dim_is_empty = builder.icmp_unsigned('==', dim_size, zero)
                with cgutils.if_unlikely(builder, dim_is_empty):
                    builder.store(cgutils.true_byte, exhausted)

            self.indices = indices
            self.exhausted = exhausted
            self.shape = cgutils.pack_array(builder, shapes, zero.type)

        def iternext_specific(self, context, builder, result):
            zero = context.get_constant(types.intp, 0)

            bbend = builder.append_basic_block('end')

            exhausted = cgutils.as_bool_bit(builder,
                                            builder.load(self.exhausted))
            with cgutils.if_unlikely(builder, exhausted):
                result.set_valid(False)
                builder.branch(bbend)

            indices = [builder.load(cgutils.gep_inbounds(builder,
                                                         self.indices,
                                                         dim))
                       for dim in range(ndim)]
            for load in indices:
                mark_positive(builder, load)

            result.yield_(cgutils.pack_array(builder, indices, zero.type))
            result.set_valid(True)

            shape = cgutils.unpack_tuple(builder, self.shape, ndim)
            _increment_indices(context, builder, ndim, shape,
                               self.indices, self.exhausted)

            builder.branch(bbend)
            builder.position_at_end(bbend)

    return NdIndexIter


def _make_flattening_iter_cls(flatiterty, kind):
    assert kind in ('flat', 'ndenumerate')

    array_type = flatiterty.array_type

    if array_type.layout == 'C':
        class CContiguousFlatIter(cgutils.create_struct_proxy(flatiterty)):
            """
            .flat() / .ndenumerate() implementation for C-contiguous arrays.
            """

            def init_specific(self, context, builder, arrty, arr):
                zero = context.get_constant(types.intp, 0)
                self.index = cgutils.alloca_once_value(builder, zero)
                # We can't trust strides[-1] to always contain the right
                # step value, see
                # http://docs.scipy.org/doc/numpy-dev/release.html#npy-relaxed-strides-checking    # noqa: E501
                self.stride = arr.itemsize

                if kind == 'ndenumerate':
                    # Zero-initialize the indices array.
                    indices = cgutils.alloca_once(
                        builder, zero.type,
                        size=context.get_constant(types.intp, arrty.ndim))

                    for dim in range(arrty.ndim):
                        idxptr = cgutils.gep_inbounds(builder, indices, dim)
                        builder.store(zero, idxptr)

                    self.indices = indices

            # NOTE: Using gep() instead of explicit pointer addition helps
            # LLVM vectorize the loop (since the stride is known and
            # constant).  This is not possible in the non-contiguous case,
            # where the strides are unknown at compile-time.

            def iternext_specific(self, context, builder, arrty, arr, result):
                ndim = arrty.ndim
                nitems = arr.nitems

                index = builder.load(self.index)
                is_valid = builder.icmp_signed('<', index, nitems)
                result.set_valid(is_valid)

                with cgutils.if_likely(builder, is_valid):
                    ptr = builder.gep(arr.data, [index])
                    value = load_item(context, builder, arrty, ptr)
                    if kind == 'flat':
                        result.yield_(value)
                    else:
                        # ndenumerate(): fetch and increment indices
                        indices = self.indices
                        idxvals = [builder.load(cgutils.gep_inbounds(builder,
                                                                     indices,
                                                                     dim))
                                   for dim in range(ndim)]
                        idxtuple = cgutils.pack_array(builder, idxvals)
                        result.yield_(
                            cgutils.make_anonymous_struct(builder,
                                                          [idxtuple, value]))
                        _increment_indices_array(context, builder, arrty,
                                                 arr, indices)

                    index = cgutils.increment_index(builder, index)
                    builder.store(index, self.index)

            def getitem(self, context, builder, arrty, arr, index):
                ptr = builder.gep(arr.data, [index])
                return load_item(context, builder, arrty, ptr)

            def setitem(self, context, builder, arrty, arr, index, value):
                ptr = builder.gep(arr.data, [index])
                store_item(context, builder, arrty, value, ptr)

        return CContiguousFlatIter

    else:
        class FlatIter(cgutils.create_struct_proxy(flatiterty)):
            """
            Generic .flat() / .ndenumerate() implementation for
            non-contiguous arrays.
            It keeps track of pointers along each dimension in order to
            minimize computations.
            """

            def init_specific(self, context, builder, arrty, arr):
                zero = context.get_constant(types.intp, 0)
                data = arr.data
                ndim = arrty.ndim
                shapes = cgutils.unpack_tuple(builder, arr.shape, ndim)

                indices = cgutils.alloca_once(
                    builder, zero.type, size=context.get_constant(types.intp,
                                                                  arrty.ndim))
                pointers = cgutils.alloca_once(
                    builder, data.type, size=context.get_constant(types.intp,
                                                                  arrty.ndim))
                exhausted = cgutils.alloca_once_value(builder,
                                                      cgutils.false_byte)

                # Initialize indices and pointers with their start values.
                for dim in range(ndim):
                    idxptr = cgutils.gep_inbounds(builder, indices, dim)
                    ptrptr = cgutils.gep_inbounds(builder, pointers, dim)
                    builder.store(data, ptrptr)
                    builder.store(zero, idxptr)
                    # 0-sized dimensions really indicate an empty array,
                    # but we have to catch that condition early to avoid
                    # a bug inside the iteration logic (see issue #846).
                    dim_size = shapes[dim]
                    dim_is_empty = builder.icmp_unsigned('==', dim_size, zero)
                    with cgutils.if_unlikely(builder, dim_is_empty):
                        builder.store(cgutils.true_byte, exhausted)

                self.indices = indices
                self.pointers = pointers
                self.exhausted = exhausted

            def iternext_specific(self, context, builder, arrty, arr, result):
                ndim = arrty.ndim
                shapes = cgutils.unpack_tuple(builder, arr.shape, ndim)
                strides = cgutils.unpack_tuple(builder, arr.strides, ndim)
                indices = self.indices
                pointers = self.pointers

                zero = context.get_constant(types.intp, 0)

                bbend = builder.append_basic_block('end')

                # Catch already computed iterator exhaustion
                is_exhausted = cgutils.as_bool_bit(
                    builder, builder.load(self.exhausted))
                with cgutils.if_unlikely(builder, is_exhausted):
                    result.set_valid(False)
                    builder.branch(bbend)
                result.set_valid(True)

                # Current pointer inside last dimension
                last_ptr = cgutils.gep_inbounds(builder, pointers, ndim - 1)
                ptr = builder.load(last_ptr)
                value = load_item(context, builder, arrty, ptr)
                if kind == 'flat':
                    result.yield_(value)
                else:
                    # ndenumerate() => yield (indices, value)
                    idxvals = [builder.load(cgutils.gep_inbounds(builder,
                                                                 indices,
                                                                 dim))
                               for dim in range(ndim)]
                    idxtuple = cgutils.pack_array(builder, idxvals)
                    result.yield_(
                        cgutils.make_anonymous_struct(builder,
                                                      [idxtuple, value]))

                # Update indices and pointers by walking from inner
                # dimension to outer.
                for dim in reversed(range(ndim)):
                    idxptr = cgutils.gep_inbounds(builder, indices, dim)
                    idx = cgutils.increment_index(builder,
                                                  builder.load(idxptr))

                    count = shapes[dim]
                    stride = strides[dim]
                    in_bounds = builder.icmp_signed('<', idx, count)
                    with cgutils.if_likely(builder, in_bounds):
                        # Index is valid => pointer can simply be incremented.
                        builder.store(idx, idxptr)
                        ptrptr = cgutils.gep_inbounds(builder, pointers, dim)
                        ptr = builder.load(ptrptr)
                        ptr = cgutils.pointer_add(builder, ptr, stride)
                        builder.store(ptr, ptrptr)
                        # Reset pointers in inner dimensions
                        for inner_dim in range(dim + 1, ndim):
                            ptrptr = cgutils.gep_inbounds(builder,
                                                          pointers,
                                                          inner_dim)
                            builder.store(ptr, ptrptr)
                        builder.branch(bbend)
                    # Reset index and continue with next dimension
                    builder.store(zero, idxptr)

                # End of array
                builder.store(cgutils.true_byte, self.exhausted)
                builder.branch(bbend)

                builder.position_at_end(bbend)

            def _ptr_for_index(self, context, builder, arrty, arr, index):
                ndim = arrty.ndim
                shapes = cgutils.unpack_tuple(builder, arr.shape, count=ndim)
                strides = cgutils.unpack_tuple(builder, arr.strides, count=ndim)

                # First convert the flattened index into a regular n-dim index
                indices = []
                for dim in reversed(range(ndim)):
                    indices.append(builder.urem(index, shapes[dim]))
                    index = builder.udiv(index, shapes[dim])
                indices.reverse()

                ptr = cgutils.get_item_pointer2(context, builder, arr.data,
                                                shapes, strides, arrty.layout,
                                                indices)
                return ptr

            def getitem(self, context, builder, arrty, arr, index):
                ptr = self._ptr_for_index(context, builder, arrty, arr, index)
                return load_item(context, builder, arrty, ptr)

            def setitem(self, context, builder, arrty, arr, index, value):
                ptr = self._ptr_for_index(context, builder, arrty, arr, index)
                store_item(context, builder, arrty, value, ptr)

        return FlatIter


@lower_getattr(types.Array, "flat")
def make_array_flatiter(context, builder, arrty, arr):
    flatitercls = make_array_flat_cls(types.NumpyFlatType(arrty))
    flatiter = flatitercls(context, builder)

    flatiter.array = arr

    arrcls = context.make_array(arrty)
    arr = arrcls(context, builder, ref=flatiter._get_ptr_by_name('array'))

    flatiter.init_specific(context, builder, arrty, arr)

    res = flatiter._getvalue()
    return impl_ret_borrowed(context, builder, types.NumpyFlatType(arrty), res)


@lower_builtin('iternext', types.NumpyFlatType)
@iternext_impl(RefType.BORROWED)
def iternext_numpy_flatiter(context, builder, sig, args, result):
    [flatiterty] = sig.args
    [flatiter] = args

    flatitercls = make_array_flat_cls(flatiterty)
    flatiter = flatitercls(context, builder, value=flatiter)

    arrty = flatiterty.array_type
    arrcls = context.make_array(arrty)
    arr = arrcls(context, builder, value=flatiter.array)

    flatiter.iternext_specific(context, builder, arrty, arr, result)


@lower_builtin(operator.getitem, types.NumpyFlatType, types.Integer)
def iternext_numpy_getitem(context, builder, sig, args):
    flatiterty = sig.args[0]
    flatiter, index = args

    flatitercls = make_array_flat_cls(flatiterty)
    flatiter = flatitercls(context, builder, value=flatiter)

    arrty = flatiterty.array_type
    arrcls = context.make_array(arrty)
    arr = arrcls(context, builder, value=flatiter.array)

    res = flatiter.getitem(context, builder, arrty, arr, index)
    return impl_ret_borrowed(context, builder, sig.return_type, res)


@lower_builtin(operator.setitem, types.NumpyFlatType, types.Integer,
               types.Any)
def iternext_numpy_getitem_any(context, builder, sig, args):
    flatiterty = sig.args[0]
    flatiter, index, value = args

    flatitercls = make_array_flat_cls(flatiterty)
    flatiter = flatitercls(context, builder, value=flatiter)

    arrty = flatiterty.array_type
    arrcls = context.make_array(arrty)
    arr = arrcls(context, builder, value=flatiter.array)

    flatiter.setitem(context, builder, arrty, arr, index, value)
    return context.get_dummy_value()


@lower_builtin(len, types.NumpyFlatType)
def iternext_numpy_getitem_flat(context, builder, sig, args):
    flatiterty = sig.args[0]
    flatitercls = make_array_flat_cls(flatiterty)
    flatiter = flatitercls(context, builder, value=args[0])

    arrcls = context.make_array(flatiterty.array_type)
    arr = arrcls(context, builder, value=flatiter.array)
    return arr.nitems


@lower_builtin(np.ndenumerate, types.Array)
def make_array_ndenumerate(context, builder, sig, args):
    arrty, = sig.args
    arr, = args
    nditercls = make_array_ndenumerate_cls(types.NumpyNdEnumerateType(arrty))
    nditer = nditercls(context, builder)

    nditer.array = arr

    arrcls = context.make_array(arrty)
    arr = arrcls(context, builder, ref=nditer._get_ptr_by_name('array'))

    nditer.init_specific(context, builder, arrty, arr)

    res = nditer._getvalue()
    return impl_ret_borrowed(context, builder, sig.return_type, res)


@lower_builtin('iternext', types.NumpyNdEnumerateType)
@iternext_impl(RefType.BORROWED)
def iternext_numpy_nditer(context, builder, sig, args, result):
    [nditerty] = sig.args
    [nditer] = args

    nditercls = make_array_ndenumerate_cls(nditerty)
    nditer = nditercls(context, builder, value=nditer)

    arrty = nditerty.array_type
    arrcls = context.make_array(arrty)
    arr = arrcls(context, builder, value=nditer.array)

    nditer.iternext_specific(context, builder, arrty, arr, result)


@lower_builtin(pndindex, types.VarArg(types.Integer))
@lower_builtin(np.ndindex, types.VarArg(types.Integer))
def make_array_ndindex(context, builder, sig, args):
    """ndindex(*shape)"""
    shape = [context.cast(builder, arg, argty, types.intp)
             for argty, arg in zip(sig.args, args)]

    nditercls = make_ndindex_cls(types.NumpyNdIndexType(len(shape)))
    nditer = nditercls(context, builder)
    nditer.init_specific(context, builder, shape)

    res = nditer._getvalue()
    return impl_ret_borrowed(context, builder, sig.return_type, res)


@lower_builtin(pndindex, types.BaseTuple)
@lower_builtin(np.ndindex, types.BaseTuple)
def make_array_ndindex_tuple(context, builder, sig, args):
    """ndindex(shape)"""
    ndim = sig.return_type.ndim
    if ndim > 0:
        idxty = sig.args[0].dtype
        tup = args[0]

        shape = cgutils.unpack_tuple(builder, tup, ndim)
        shape = [context.cast(builder, idx, idxty, types.intp)
                 for idx in shape]
    else:
        shape = []

    nditercls = make_ndindex_cls(types.NumpyNdIndexType(len(shape)))
    nditer = nditercls(context, builder)
    nditer.init_specific(context, builder, shape)

    res = nditer._getvalue()
    return impl_ret_borrowed(context, builder, sig.return_type, res)


@lower_builtin('iternext', types.NumpyNdIndexType)
@iternext_impl(RefType.BORROWED)
def iternext_numpy_ndindex(context, builder, sig, args, result):
    [nditerty] = sig.args
    [nditer] = args

    nditercls = make_ndindex_cls(nditerty)
    nditer = nditercls(context, builder, value=nditer)

    nditer.iternext_specific(context, builder, result)


@lower_builtin(np.nditer, types.Any)
def make_array_nditer(context, builder, sig, args):
    """
    nditer(...)
    """
    nditerty = sig.return_type
    arrtys = nditerty.arrays

    if isinstance(sig.args[0], types.BaseTuple):
        arrays = cgutils.unpack_tuple(builder, args[0])
    else:
        arrays = [args[0]]

    nditer = make_nditer_cls(nditerty)(context, builder)
    nditer.init_specific(context, builder, arrtys, arrays)

    res = nditer._getvalue()
    return impl_ret_borrowed(context, builder, nditerty, res)


@lower_builtin('iternext', types.NumpyNdIterType)
@iternext_impl(RefType.BORROWED)
def iternext_numpy_nditer2(context, builder, sig, args, result):
    [nditerty] = sig.args
    [nditer] = args

    nditer = make_nditer_cls(nditerty)(context, builder, value=nditer)
    nditer.iternext_specific(context, builder, result)


@lower_builtin(operator.eq, types.DType, types.DType)
def dtype_eq_impl(context, builder, sig, args):
    arg1, arg2 = sig.args
    res = ir.Constant(ir.IntType(1), int(arg1 == arg2))
    return impl_ret_untracked(context, builder, sig.return_type, res)


# ------------------------------------------------------------------------------
# Numpy array constructors

def _empty_nd_impl(context, builder, arrtype, shapes):
    """Utility function used for allocating a new array during LLVM code
    generation (lowering).  Given a target context, builder, array
    type, and a tuple or list of lowered dimension sizes, returns a
    LLVM value pointing at a Numba runtime allocated array.
    """
    arycls = make_array(arrtype)
    ary = arycls(context, builder)

    datatype = context.get_data_type(arrtype.dtype)
    itemsize = context.get_constant(types.intp, get_itemsize(context, arrtype))

    # compute array length
    arrlen = context.get_constant(types.intp, 1)
    overflow = Constant(ir.IntType(1), 0)
    for s in shapes:
        arrlen_mult = builder.smul_with_overflow(arrlen, s)
        arrlen = builder.extract_value(arrlen_mult, 0)
        overflow = builder.or_(
            overflow, builder.extract_value(arrlen_mult, 1)
        )

    if arrtype.ndim == 0:
        strides = ()
    elif arrtype.layout == 'C':
        strides = [itemsize]
        for dimension_size in reversed(shapes[1:]):
            strides.append(builder.mul(strides[-1], dimension_size))
        strides = tuple(reversed(strides))
    elif arrtype.layout == 'F':
        strides = [itemsize]
        for dimension_size in shapes[:-1]:
            strides.append(builder.mul(strides[-1], dimension_size))
        strides = tuple(strides)
    else:
        raise NotImplementedError(
            "Don't know how to allocate array with layout '{0}'.".format(
                arrtype.layout))

    # Check overflow, numpy also does this after checking order
    allocsize_mult = builder.smul_with_overflow(arrlen, itemsize)
    allocsize = builder.extract_value(allocsize_mult, 0)
    overflow = builder.or_(overflow, builder.extract_value(allocsize_mult, 1))

    with builder.if_then(overflow, likely=False):
        # Raise same error as numpy, see:
        # https://github.com/numpy/numpy/blob/2a488fe76a0f732dc418d03b452caace161673da/numpy/core/src/multiarray/ctors.c#L1095-L1101    # noqa: E501
        context.call_conv.return_user_exc(
            builder, ValueError,
            ("array is too big; `arr.size * arr.dtype.itemsize` is larger than"
             " the maximum possible size.",)
        )

    dtype = arrtype.dtype
    align_val = context.get_preferred_array_alignment(dtype)
    align = context.get_constant(types.uint32, align_val)
    args = (context.get_dummy_value(), allocsize, align)

    mip = types.MemInfoPointer(types.voidptr)
    arytypeclass = types.TypeRef(type(arrtype))
    argtypes = signature(mip, arytypeclass, types.intp, types.uint32)

    meminfo = context.compile_internal(builder, _call_allocator, argtypes, args)
    data = context.nrt.meminfo_data(builder, meminfo)

    intp_t = context.get_value_type(types.intp)
    shape_array = cgutils.pack_array(builder, shapes, ty=intp_t)
    strides_array = cgutils.pack_array(builder, strides, ty=intp_t)

    populate_array(ary,
                   data=builder.bitcast(data, datatype.as_pointer()),
                   shape=shape_array,
                   strides=strides_array,
                   itemsize=itemsize,
                   meminfo=meminfo)

    return ary


@overload_classmethod(types.Array, "_allocate")
def _ol_array_allocate(cls, allocsize, align):
    """Implements a Numba-only default target (cpu) classmethod on the array
    type.
    """
    def impl(cls, allocsize, align):
        return intrin_alloc(allocsize, align)
    return impl


def _call_allocator(arrtype, size, align):
    """Trampoline to call the intrinsic used for allocation
    """
    return arrtype._allocate(size, align)


@intrinsic
def intrin_alloc(typingctx, allocsize, align):
    """Intrinsic to call into the allocator for Array
    """
    def codegen(context, builder, signature, args):
        [allocsize, align] = args
        meminfo = context.nrt.meminfo_alloc_aligned(builder, allocsize, align)
        return meminfo

    mip = types.MemInfoPointer(types.voidptr)    # return untyped pointer
    sig = signature(mip, allocsize, align)
    return sig, codegen


def _parse_shape(context, builder, ty, val):
    """
    Parse the shape argument to an array constructor.
    """
    def safecast_intp(context, builder, src_t, src):
        """Cast src to intp only if value can be maintained"""
        intp_t = context.get_value_type(types.intp)
        intp_width = intp_t.width
        intp_ir = ir.IntType(intp_width)
        maxval = Constant(intp_ir, ((1 << intp_width - 1) - 1))
        if src_t.width < intp_width:
            res = builder.sext(src, intp_ir)
        elif src_t.width >= intp_width:
            is_larger = builder.icmp_signed(">", src, maxval)
            with builder.if_then(is_larger, likely=False):
                context.call_conv.return_user_exc(
                    builder, ValueError,
                    ("Cannot safely convert value to intp",)
                )
            if src_t.width > intp_width:
                res = builder.trunc(src, intp_ir)
            else:
                res = src
        return res

    if isinstance(ty, types.Integer):
        ndim = 1
        passed_shapes = [context.cast(builder, val, ty, types.intp)]
    else:
        assert isinstance(ty, types.BaseTuple)
        ndim = ty.count
        passed_shapes = cgutils.unpack_tuple(builder, val, count=ndim)

    shapes = []
    for s in passed_shapes:
        shapes.append(safecast_intp(context, builder, s.type, s))

    zero = context.get_constant_generic(builder, types.intp, 0)
    for dim in range(ndim):
        is_neg = builder.icmp_signed('<', shapes[dim], zero)
        with cgutils.if_unlikely(builder, is_neg):
            context.call_conv.return_user_exc(
                builder, ValueError, ("negative dimensions not allowed",)
            )

    return shapes


def _parse_empty_args(context, builder, sig, args):
    """
    Parse the arguments of a np.empty(), np.zeros() or np.ones() call.
    """
    arrshapetype = sig.args[0]
    arrshape = args[0]
    arrtype = sig.return_type
    return arrtype, _parse_shape(context, builder, arrshapetype, arrshape)


def _parse_empty_like_args(context, builder, sig, args):
    """
    Parse the arguments of a np.empty_like(), np.zeros_like() or
    np.ones_like() call.
    """
    arytype = sig.args[0]
    if isinstance(arytype, types.Array):
        ary = make_array(arytype)(context, builder, value=args[0])
        shapes = cgutils.unpack_tuple(builder, ary.shape, count=arytype.ndim)
        return sig.return_type, shapes
    else:
        return sig.return_type, ()


def _check_const_str_dtype(fname, dtype):
    if isinstance(dtype, types.UnicodeType):
        msg = f"If np.{fname} dtype is a string it must be a string constant."
        raise errors.TypingError(msg)


@intrinsic
def numpy_empty_nd(tyctx, ty_shape, ty_dtype, ty_retty_ref):
    ty_retty = ty_retty_ref.instance_type
    sig = ty_retty(ty_shape, ty_dtype, ty_retty_ref)

    def codegen(cgctx, builder, sig, llargs):
        arrtype, shapes = _parse_empty_args(cgctx, builder, sig, llargs)
        ary = _empty_nd_impl(cgctx, builder, arrtype, shapes)
        return ary._getvalue()
    return sig, codegen


@overload(np.empty)
def ol_np_empty(shape, dtype=float):
    _check_const_str_dtype("empty", dtype)
    if (dtype is float or
        (isinstance(dtype, types.Function) and dtype.typing_key is float) or
            is_nonelike(dtype)): #default
        nb_dtype = types.double
    else:
        nb_dtype = ty_parse_dtype(dtype)

    ndim = ty_parse_shape(shape)
    if nb_dtype is not None and ndim is not None:
        retty = types.Array(dtype=nb_dtype, ndim=ndim, layout='C')

        def impl(shape, dtype=float):
            return numpy_empty_nd(shape, dtype, retty)
        return impl
    else:
        msg = f"Cannot parse input types to function np.empty({shape}, {dtype})"
        raise errors.TypingError(msg)


@intrinsic
def numpy_empty_like_nd(tyctx, ty_prototype, ty_dtype, ty_retty_ref):
    ty_retty = ty_retty_ref.instance_type
    sig = ty_retty(ty_prototype, ty_dtype, ty_retty_ref)

    def codegen(cgctx, builder, sig, llargs):
        arrtype, shapes = _parse_empty_like_args(cgctx, builder, sig, llargs)
        ary = _empty_nd_impl(cgctx, builder, arrtype, shapes)
        return ary._getvalue()
    return sig, codegen


@overload(np.empty_like)
def ol_np_empty_like(arr, dtype=None):
    _check_const_str_dtype("empty_like", dtype)
    if not is_nonelike(dtype):
        nb_dtype = ty_parse_dtype(dtype)
    elif isinstance(arr, types.Array):
        nb_dtype = arr.dtype
    else:
        nb_dtype = arr
    if nb_dtype is not None:
        if isinstance(arr, types.Array):
            layout = arr.layout if arr.layout != 'A' else 'C'
            retty = arr.copy(dtype=nb_dtype, layout=layout, readonly=False)
        else:
            retty = types.Array(nb_dtype, 0, 'C')
    else:
        msg = ("Cannot parse input types to function "
               f"np.empty_like({arr}, {dtype})")
        raise errors.TypingError(msg)

    def impl(arr, dtype=None):
        return numpy_empty_like_nd(arr, dtype, retty)
    return impl


@intrinsic
def _zero_fill_array_method(tyctx, self):
    sig = types.none(self)

    def codegen(cgctx, builder, sig, llargs):
        ary = make_array(sig.args[0])(cgctx, builder, llargs[0])
        cgutils.memset(builder, ary.data, builder.mul(ary.itemsize, ary.nitems),
                       0)
    return sig, codegen


@overload_method(types.Array, '_zero_fill')
def ol_array_zero_fill(self):
    """Adds a `._zero_fill` method to zero fill an array using memset."""
    def impl(self):
        _zero_fill_array_method(self)
    return impl


@overload(np.zeros)
def ol_np_zeros(shape, dtype=float):
    _check_const_str_dtype("zeros", dtype)

    def impl(shape, dtype=float):
        arr = np.empty(shape, dtype=dtype)
        arr._zero_fill()
        return arr
    return impl


@overload(np.zeros_like)
def ol_np_zeros_like(a, dtype=None):
    _check_const_str_dtype("zeros_like", dtype)

    # NumPy uses 'a' as the arg name for the array-like
    def impl(a, dtype=None):
        arr = np.empty_like(a, dtype=dtype)
        arr._zero_fill()
        return arr
    return impl


@overload(np.ones_like)
def ol_np_ones_like(a, dtype=None):
    _check_const_str_dtype("ones_like", dtype)

    # NumPy uses 'a' as the arg name for the array-like
    def impl(a, dtype=None):
        arr = np.empty_like(a, dtype=dtype)
        arr_flat = arr.flat
        for idx in range(len(arr_flat)):
            arr_flat[idx] = 1
        return arr
    return impl


@overload(np.full)
def impl_np_full(shape, fill_value, dtype=None):
    _check_const_str_dtype("full", dtype)
    if not is_nonelike(dtype):
        nb_dtype = ty_parse_dtype(dtype)
    else:
        nb_dtype = fill_value

    def full(shape, fill_value, dtype=None):
        arr = np.empty(shape, nb_dtype)
        arr_flat = arr.flat
        for idx in range(len(arr_flat)):
            arr_flat[idx] = fill_value
        return arr
    return full


@overload(np.full_like)
def impl_np_full_like(a, fill_value, dtype=None):
    _check_const_str_dtype("full_like", dtype)

    def full_like(a, fill_value, dtype=None):
        arr = np.empty_like(a, dtype)
        arr_flat = arr.flat
        for idx in range(len(arr_flat)):
            arr_flat[idx] = fill_value
        return arr

    return full_like


@overload(np.ones)
def ol_np_ones(shape, dtype=None):
    # for some reason the NumPy default for dtype is None in the source but
    # ends up as np.float64 by definition.
    _check_const_str_dtype("ones", dtype)

    def impl(shape, dtype=None):
        arr = np.empty(shape, dtype=dtype)
        arr_flat = arr.flat
        for idx in range(len(arr_flat)):
            arr_flat[idx] = 1
        return arr
    return impl


@overload(np.identity)
def impl_np_identity(n, dtype=None):
    _check_const_str_dtype("identity", dtype)
    if not is_nonelike(dtype):
        nb_dtype = ty_parse_dtype(dtype)
    else:
        nb_dtype = types.double

    def identity(n, dtype=None):
        arr = np.zeros((n, n), nb_dtype)
        for i in range(n):
            arr[i, i] = 1
        return arr
    return identity


def _eye_none_handler(N, M):
    pass


@extending.overload(_eye_none_handler)
def _eye_none_handler_impl(N, M):
    if isinstance(M, types.NoneType):
        def impl(N, M):
            return N
    else:
        def impl(N, M):
            return M
    return impl


@extending.overload(np.eye)
def numpy_eye(N, M=None, k=0, dtype=float):

    if dtype is None or isinstance(dtype, types.NoneType):
        dt = np.dtype(float)
    elif isinstance(dtype, (types.DTypeSpec, types.Number)):
        # dtype or instance of dtype
        dt = as_dtype(getattr(dtype, 'dtype', dtype))
    else:
        dt = np.dtype(dtype)

    def impl(N, M=None, k=0, dtype=float):
        _M = _eye_none_handler(N, M)
        arr = np.zeros((N, _M), dt)
        if k >= 0:
            d = min(N, _M - k)
            for i in range(d):
                arr[i, i + k] = 1
        else:
            d = min(N + k, _M)
            for i in range(d):
                arr[i - k, i] = 1
        return arr
    return impl


@overload(np.diag)
def impl_np_diag(v, k=0):
    if not type_can_asarray(v):
        raise errors.TypingError('The argument "v" must be array-like')

    if isinstance(v, types.Array):
        if v.ndim not in (1, 2):
            raise errors.NumbaTypeError("Input must be 1- or 2-d.")

        def diag_impl(v, k=0):
            if v.ndim == 1:
                s = v.shape
                n = s[0] + abs(k)
                ret = np.zeros((n, n), v.dtype)
                if k >= 0:
                    for i in range(n - k):
                        ret[i, k + i] = v[i]
                else:
                    for i in range(n + k):
                        ret[i - k, i] = v[i]
                return ret
            else:  # 2-d
                rows, cols = v.shape
                if k < 0:
                    rows = rows + k
                if k > 0:
                    cols = cols - k
                n = max(min(rows, cols), 0)
                ret = np.empty(n, v.dtype)
                if k >= 0:
                    for i in range(n):
                        ret[i] = v[i, k + i]
                else:
                    for i in range(n):
                        ret[i] = v[i - k, i]
                return ret
        return diag_impl


@overload(np.indices)
def numpy_indices(dimensions):
    if not isinstance(dimensions, types.UniTuple):
        msg = 'The argument "dimensions" must be a tuple of integers'
        raise errors.TypingError(msg)

    if not isinstance(dimensions.dtype, types.Integer):
        msg = 'The argument "dimensions" must be a tuple of integers'
        raise errors.TypingError(msg)

    N = len(dimensions)
    shape = (1,) * N

    def impl(dimensions):
        res = np.empty((N,) + dimensions, dtype=np.int64)
        i = 0
        for dim in dimensions:
            idx = np.arange(dim, dtype=np.int64).reshape(
                tuple_setitem(shape, i, dim)
            )
            res[i] = idx
            i += 1

        return res

    return impl


@overload(np.diagflat)
def numpy_diagflat(v, k=0):
    if not type_can_asarray(v):
        msg = 'The argument "v" must be array-like'
        raise errors.TypingError(msg)

    if not isinstance(k, (int, types.Integer)):
        msg = 'The argument "k" must be an integer'
        raise errors.TypingError(msg)

    def impl(v, k=0):
        v = np.asarray(v)
        v = v.ravel()
        s = len(v)
        abs_k = abs(k)
        n = s + abs_k
        res = np.zeros((n, n), v.dtype)
        i = np.maximum(0, -k)
        j = np.maximum(0, k)
        for t in range(s):
            res[i + t, j + t] = v[t]

        return res

    return impl


@overload(np.take)
@overload_method(types.Array, 'take')
def numpy_take(a, indices):

    if isinstance(a, types.Array) and isinstance(indices, types.Integer):
        def take_impl(a, indices):
            if indices > (a.size - 1) or indices < -a.size:
                raise IndexError("Index out of bounds")
            return a.ravel()[indices]
        return take_impl

    if all(isinstance(arg, types.Array) for arg in [a, indices]):
        F_order = indices.layout == 'F'

        def take_impl(a, indices):
            ret = np.empty(indices.size, dtype=a.dtype)
            if F_order:
                walker = indices.copy()  # get C order
            else:
                walker = indices
            it = np.nditer(walker)
            i = 0
            flat = a.ravel()
            for x in it:
                if x > (a.size - 1) or x < -a.size:
                    raise IndexError("Index out of bounds")
                ret[i] = flat[x]
                i = i + 1
            return ret.reshape(indices.shape)
        return take_impl

    if isinstance(a, types.Array) and \
            isinstance(indices, (types.List, types.BaseTuple)):
        def take_impl(a, indices):
            convert = np.array(indices)
            ret = np.empty(convert.size, dtype=a.dtype)
            it = np.nditer(convert)
            i = 0
            flat = a.ravel()
            for x in it:
                if x > (a.size - 1) or x < -a.size:
                    raise IndexError("Index out of bounds")
                ret[i] = flat[x]
                i = i + 1
            return ret.reshape(convert.shape)
        return take_impl


def _arange_dtype(*args):
    bounds = [a for a in args if not isinstance(a, types.NoneType)]

    if any(isinstance(a, types.Complex) for a in bounds):
        dtype = types.complex128
    elif any(isinstance(a, types.Float) for a in bounds):
        dtype = types.float64
    else:
        # `np.arange(10).dtype` is always `np.dtype(int)`, aka `np.int_`, which
        # in all released versions of numpy corresponds to the C `long` type.
        # Windows 64 is broken by default here because Numba (as of 0.47) does
        # not differentiate between Python and NumPy integers, so a `typeof(1)`
        # on w64 is `int64`, i.e. `intp`. This means an arange(<some int>) will
        # be typed as arange(int64) and the following will yield int64 opposed
        # to int32. Example: without a load of analysis to work out of the args
        # were wrapped in NumPy int*() calls it's not possible to detect the
        # difference between `np.arange(10)` and `np.arange(np.int64(10)`.
        NPY_TY = getattr(types, "int%s" % (8 * np.dtype(int).itemsize))

        # unliteral these types such that `max` works.
        unliteral_bounds = [types.unliteral(x) for x in bounds]
        dtype = max(unliteral_bounds + [NPY_TY,])

    return dtype


@overload(np.arange)
def np_arange(start, / ,stop=None, step=None, dtype=None):
    if isinstance(stop, types.Optional):
        stop = stop.type
    if isinstance(step, types.Optional):
        step = step.type
    if isinstance(dtype, types.Optional):
        dtype = dtype.type

    if stop is None:
        stop = types.none
    if step is None:
        step = types.none
    if dtype is None:
        dtype = types.none

    if (not isinstance(start, types.Number) or
        not isinstance(stop, (types.NoneType, types.Number)) or
        not isinstance(step, (types.NoneType, types.Number)) or
            not isinstance(dtype, (types.NoneType, types.DTypeSpec))):

        return

    if isinstance(dtype, types.NoneType):
        true_dtype = _arange_dtype(start, stop, step)
    else:
        true_dtype = dtype.dtype

    use_complex = any([isinstance(x, types.Complex)
                       for x in (start, stop, step)])

    start_value = getattr(start, "literal_value", None)
    stop_value = getattr(stop, "literal_value", None)
    step_value = getattr(step, "literal_value", None)

    def impl(start, /, stop=None, step=None, dtype=None):
        # Allow for improved performance if given literal arguments.
        lit_start = start_value if start_value is not None else start
        lit_stop = stop_value if stop_value is not None else stop
        lit_step = step_value if step_value is not None else step

        _step = lit_step if lit_step is not None else 1
        if lit_stop is None:
            _start, _stop = 0, lit_start
        else:
            _start, _stop = lit_start, lit_stop

        if _step == 0:
            raise ValueError("Maximum allowed size exceeded")

        nitems_c = (_stop - _start) / _step
        nitems_r = int(math.ceil(nitems_c.real))

        # Binary operator needed for compiler branch pruning.
        if use_complex is True:
            nitems_i = int(math.ceil(nitems_c.imag))
            nitems = max(min(nitems_i, nitems_r), 0)
        else:
            nitems = max(nitems_r, 0)
        arr = np.empty(nitems, true_dtype)
        val = _start
        for i in range(nitems):
            arr[i] = val + (i * _step)
        return arr

    return impl


@overload(np.linspace)
def numpy_linspace(start, stop, num=50):
    if not all(isinstance(arg, types.Number) for arg in [start, stop]):
        return

    if not isinstance(num, (int, types.Integer)):
        msg = 'The argument "num" must be an integer'
        raise errors.TypingError(msg)

    if any(isinstance(arg, types.Complex) for arg in [start, stop]):
        dtype = types.complex128
    else:
        dtype = types.float64

    # Implementation based on https://github.com/numpy/numpy/blob/v1.20.0/numpy/core/function_base.py#L24 # noqa: E501
    def linspace(start, stop, num=50):
        arr = np.empty(num, dtype)
        # The multiply by 1.0 mirrors
        # https://github.com/numpy/numpy/blob/v1.20.0/numpy/core/function_base.py#L125-L128  # noqa: E501
        # the side effect of this is important... start and stop become the same
        # type as `dtype` i.e. 64/128 bits wide (float/complex). This is
        # important later when used in the `np.divide`.
        start = start * 1.0
        stop = stop * 1.0
        if num == 0:
            return arr
        div = num - 1
        if div > 0:
            delta = stop - start
            step = np.divide(delta, div)
            for i in range(0, num):
                arr[i] = start + (i * step)
        else:
            arr[0] = start
        if num > 1:
            arr[-1] = stop
        return arr
    return linspace


def _array_copy(context, builder, sig, args):
    """
    Array copy.
    """
    arytype = sig.args[0]
    ary = make_array(arytype)(context, builder, value=args[0])
    shapes = cgutils.unpack_tuple(builder, ary.shape)

    rettype = sig.return_type
    ret = _empty_nd_impl(context, builder, rettype, shapes)

    src_data = ary.data
    dest_data = ret.data

    assert rettype.layout in "CF"
    if arytype.layout == rettype.layout:
        # Fast path: memcpy
        cgutils.raw_memcpy(builder, dest_data, src_data, ary.nitems,
                           ary.itemsize, align=1)

    else:
        src_strides = cgutils.unpack_tuple(builder, ary.strides)
        dest_strides = cgutils.unpack_tuple(builder, ret.strides)
        intp_t = context.get_value_type(types.intp)

        with cgutils.loop_nest(builder, shapes, intp_t) as indices:
            src_ptr = cgutils.get_item_pointer2(context, builder, src_data,
                                                shapes, src_strides,
                                                arytype.layout, indices)
            dest_ptr = cgutils.get_item_pointer2(context, builder, dest_data,
                                                 shapes, dest_strides,
                                                 rettype.layout, indices)
            builder.store(builder.load(src_ptr), dest_ptr)

    return impl_ret_new_ref(context, builder, sig.return_type, ret._getvalue())


@intrinsic
def _array_copy_intrinsic(typingctx, a):
    assert isinstance(a, types.Array)
    layout = 'F' if a.layout == 'F' else 'C'
    ret = a.copy(layout=layout, readonly=False)
    sig = ret(a)
    return sig, _array_copy


@lower_builtin("array.copy", types.Array)
def array_copy(context, builder, sig, args):
    return _array_copy(context, builder, sig, args)


@overload(np.copy)
def impl_numpy_copy(a):
    if isinstance(a, types.Array):
        def numpy_copy(a):
            return _array_copy_intrinsic(a)
    return numpy_copy


def _as_layout_array(context, builder, sig, args, output_layout):
    """
    Common logic for layout conversion function;
    e.g. ascontiguousarray and asfortranarray
    """
    retty = sig.return_type
    aryty = sig.args[0]
    assert retty.layout == output_layout, 'return-type has incorrect layout'

    if aryty.ndim == 0:
        # 0-dim input => asfortranarray() returns a 1-dim array
        assert retty.ndim == 1
        ary = make_array(aryty)(context, builder, value=args[0])
        ret = make_array(retty)(context, builder)

        shape = context.get_constant_generic(
            builder, types.UniTuple(types.intp, 1), (1,),
        )
        strides = context.make_tuple(builder,
                                     types.UniTuple(types.intp, 1),
                                     (ary.itemsize,))
        populate_array(ret, ary.data, shape, strides, ary.itemsize,
                       ary.meminfo, ary.parent)
        return impl_ret_borrowed(context, builder, retty, ret._getvalue())

    elif (retty.layout == aryty.layout
            or (aryty.ndim == 1 and aryty.layout in 'CF')):
        # 1-dim contiguous input => return the same array
        return impl_ret_borrowed(context, builder, retty, args[0])

    else:
        if aryty.layout == 'A':
            # There's still chance the array is in contiguous layout,
            # just that we don't know at compile time.
            # We can do a runtime check.

            # Prepare and call is_contiguous or is_fortran
            assert output_layout in 'CF'
            check_func = is_contiguous if output_layout == 'C' else is_fortran
            is_contig = _call_contiguous_check(check_func,
                                               context,
                                               builder,
                                               aryty,
                                               args[0])
            with builder.if_else(is_contig) as (then, orelse):
                # If the array is already contiguous, just return it
                with then:
                    out_then = impl_ret_borrowed(context, builder, retty,
                                                 args[0])
                    then_blk = builder.block
                # Otherwise, copy to a new contiguous region
                with orelse:
                    out_orelse = _array_copy(context, builder, sig, args)
                    orelse_blk = builder.block
            # Phi node for the return value
            ret_phi = builder.phi(out_then.type)
            ret_phi.add_incoming(out_then, then_blk)
            ret_phi.add_incoming(out_orelse, orelse_blk)
            return ret_phi

        else:
            # Return a copy with the right layout
            return _array_copy(context, builder, sig, args)


@intrinsic
def _as_layout_array_intrinsic(typingctx, a, output_layout):
    if not isinstance(output_layout, types.StringLiteral):
        raise errors.RequireLiteralValue(output_layout)

    ret = a.copy(layout=output_layout.literal_value, ndim=max(a.ndim, 1))
    sig = ret(a, output_layout)

    return sig, lambda c, b, s, a: _as_layout_array(
        c, b, s, a, output_layout=output_layout.literal_value)


@overload(np.ascontiguousarray)
def array_ascontiguousarray(a):
    if not type_can_asarray(a):
        raise errors.TypingError('The argument "a" must be array-like')

    if isinstance(a, (types.Number, types.Boolean,)):
        def impl(a):
            return np.ascontiguousarray(np.array(a))
    elif isinstance(a, types.Array):
        def impl(a):
            return _as_layout_array_intrinsic(a, 'C')
    return impl


@overload(np.asfortranarray)
def array_asfortranarray(a):
    if not type_can_asarray(a):
        raise errors.TypingError('The argument "a" must be array-like')

    if isinstance(a, (types.Number, types.Boolean,)):
        def impl(a):
            return np.asfortranarray(np.array(a))
        return impl
    elif isinstance(a, types.Array):
        def impl(a):
            return _as_layout_array_intrinsic(a, 'F')
        return impl


@lower_builtin("array.astype", types.Array, types.DTypeSpec)
@lower_builtin("array.astype", types.Array, types.StringLiteral)
def array_astype(context, builder, sig, args):
    arytype = sig.args[0]
    ary = make_array(arytype)(context, builder, value=args[0])
    shapes = cgutils.unpack_tuple(builder, ary.shape)

    rettype = sig.return_type
    ret = _empty_nd_impl(context, builder, rettype, shapes)

    src_data = ary.data
    dest_data = ret.data

    src_strides = cgutils.unpack_tuple(builder, ary.strides)
    dest_strides = cgutils.unpack_tuple(builder, ret.strides)
    intp_t = context.get_value_type(types.intp)

    with cgutils.loop_nest(builder, shapes, intp_t) as indices:
        src_ptr = cgutils.get_item_pointer2(context, builder, src_data,
                                            shapes, src_strides,
                                            arytype.layout, indices)
        dest_ptr = cgutils.get_item_pointer2(context, builder, dest_data,
                                             shapes, dest_strides,
                                             rettype.layout, indices)
        item = load_item(context, builder, arytype, src_ptr)
        item = context.cast(builder, item, arytype.dtype, rettype.dtype)
        store_item(context, builder, rettype, item, dest_ptr)

    return impl_ret_new_ref(context, builder, sig.return_type, ret._getvalue())


@intrinsic
def np_frombuffer(typingctx, buffer, dtype, retty):
    ty = retty.instance_type
    sig = ty(buffer, dtype, retty)

    def codegen(context, builder, sig, args):
        bufty = sig.args[0]
        aryty = sig.return_type

        buf = make_array(bufty)(context, builder, value=args[0])
        out_ary_ty = make_array(aryty)
        out_ary = out_ary_ty(context, builder)
        out_datamodel = out_ary._datamodel

        itemsize = get_itemsize(context, aryty)
        ll_itemsize = Constant(buf.itemsize.type, itemsize)
        nbytes = builder.mul(buf.nitems, buf.itemsize)

        # Check that the buffer size is compatible
        rem = builder.srem(nbytes, ll_itemsize)
        is_incompatible = cgutils.is_not_null(builder, rem)
        with builder.if_then(is_incompatible, likely=False):
            msg = "buffer size must be a multiple of element size"
            context.call_conv.return_user_exc(builder, ValueError, (msg,))

        shape = cgutils.pack_array(builder, [builder.sdiv(nbytes, ll_itemsize)])
        strides = cgutils.pack_array(builder, [ll_itemsize])
        data = builder.bitcast(
            buf.data, context.get_value_type(out_datamodel.get_type('data'))
        )

        populate_array(out_ary,
                       data=data,
                       shape=shape,
                       strides=strides,
                       itemsize=ll_itemsize,
                       meminfo=buf.meminfo,
                       parent=buf.parent,)

        res = out_ary._getvalue()
        return impl_ret_borrowed(context, builder, sig.return_type, res)
    return sig, codegen


@overload(np.frombuffer)
def impl_np_frombuffer(buffer, dtype=float):
    _check_const_str_dtype("frombuffer", dtype)

    if not isinstance(buffer, types.Buffer) or buffer.layout != 'C':
        msg = f'Argument "buffer" must be buffer-like. Got {buffer}'
        raise errors.TypingError(msg)

    if (dtype is float or
        (isinstance(dtype, types.Function) and dtype.typing_key is float) or
            is_nonelike(dtype)): #default
        nb_dtype = types.double
    else:
        nb_dtype = ty_parse_dtype(dtype)

    if nb_dtype is not None:
        retty = types.Array(dtype=nb_dtype, ndim=1, layout='C',
                            readonly=not buffer.mutable)
    else:
        msg = ("Cannot parse input types to function "
               f"np.frombuffer({buffer}, {dtype})")
        raise errors.TypingError(msg)

    def impl(buffer, dtype=float):
        return np_frombuffer(buffer, dtype, retty)
    return impl


@overload(carray)
def impl_carray(ptr, shape, dtype=None):
    if is_nonelike(dtype):
        intrinsic_cfarray = get_cfarray_intrinsic('C', None)

        def impl(ptr, shape, dtype=None):
            return intrinsic_cfarray(ptr, shape)
        return impl
    elif isinstance(dtype, types.DTypeSpec):
        intrinsic_cfarray = get_cfarray_intrinsic('C', dtype)

        def impl(ptr, shape, dtype=None):
            return intrinsic_cfarray(ptr, shape)
        return impl


@overload(farray)
def impl_farray(ptr, shape, dtype=None):
    if is_nonelike(dtype):
        intrinsic_cfarray = get_cfarray_intrinsic('F', None)

        def impl(ptr, shape, dtype=None):
            return intrinsic_cfarray(ptr, shape)
        return impl
    elif isinstance(dtype, types.DTypeSpec):
        intrinsic_cfarray = get_cfarray_intrinsic('F', dtype)

        def impl(ptr, shape, dtype=None):
            return intrinsic_cfarray(ptr, shape)
        return impl


def get_cfarray_intrinsic(layout, dtype_):
    @intrinsic
    def intrinsic_cfarray(typingctx, ptr, shape):
        if ptr is types.voidptr:
            ptr_dtype = None
        elif isinstance(ptr, types.CPointer):
            ptr_dtype = ptr.dtype
        else:
            msg = f"pointer argument expected, got '{ptr}'"
            raise errors.NumbaTypeError(msg)

        if dtype_ is None:
            if ptr_dtype is None:
                msg = "explicit dtype required for void* argument"
                raise errors.NumbaTypeError(msg)
            dtype = ptr_dtype
        elif isinstance(dtype_, types.DTypeSpec):
            dtype = dtype_.dtype
            if ptr_dtype is not None and dtype != ptr_dtype:
                msg = f"mismatching dtype '{dtype}' for pointer type '{ptr}'"
                raise errors.NumbaTypeError(msg)
        else:
            msg = f"invalid dtype spec '{dtype_}'"
            raise errors.NumbaTypeError(msg)

        ndim = ty_parse_shape(shape)
        if ndim is None:
            msg = f"invalid shape '{shape}'"
            raise errors.NumbaTypeError(msg)

        retty = types.Array(dtype, ndim, layout)
        sig = signature(retty, ptr, shape)
        return sig, np_cfarray
    return intrinsic_cfarray


def np_cfarray(context, builder, sig, args):
    """
    numba.numpy_support.carray(...) and
    numba.numpy_support.farray(...).
    """
    ptrty, shapety = sig.args[:2]
    ptr, shape = args[:2]

    aryty = sig.return_type
    assert aryty.layout in 'CF'

    out_ary = make_array(aryty)(context, builder)

    itemsize = get_itemsize(context, aryty)
    ll_itemsize = cgutils.intp_t(itemsize)

    if isinstance(shapety, types.BaseTuple):
        shapes = cgutils.unpack_tuple(builder, shape)
    else:
        shapety = (shapety,)
        shapes = (shape,)
    shapes = [context.cast(builder, value, fromty, types.intp)
              for fromty, value in zip(shapety, shapes)]

    off = ll_itemsize
    strides = []
    if aryty.layout == 'F':
        for s in shapes:
            strides.append(off)
            off = builder.mul(off, s)
    else:
        for s in reversed(shapes):
            strides.append(off)
            off = builder.mul(off, s)
        strides.reverse()

    data = builder.bitcast(ptr,
                           context.get_data_type(aryty.dtype).as_pointer())

    populate_array(out_ary,
                   data=data,
                   shape=shapes,
                   strides=strides,
                   itemsize=ll_itemsize,
                   # Array is not memory-managed
                   meminfo=None,
                   )

    res = out_ary._getvalue()
    return impl_ret_new_ref(context, builder, sig.return_type, res)


def _get_seq_size(context, builder, seqty, seq):
    if isinstance(seqty, types.BaseTuple):
        return context.get_constant(types.intp, len(seqty))
    elif isinstance(seqty, types.Sequence):
        len_impl = context.get_function(len, signature(types.intp, seqty,))
        return len_impl(builder, (seq,))
    else:
        assert 0


def _get_borrowing_getitem(context, seqty):
    """
    Return a getitem() implementation that doesn't incref its result.
    """
    retty = seqty.dtype
    getitem_impl = context.get_function(operator.getitem,
                                        signature(retty, seqty, types.intp))

    def wrap(builder, args):
        ret = getitem_impl(builder, args)
        if context.enable_nrt:
            context.nrt.decref(builder, retty, ret)
        return ret

    return wrap


def compute_sequence_shape(context, builder, ndim, seqty, seq):
    """
    Compute the likely shape of a nested sequence (possibly 0d).
    """
    intp_t = context.get_value_type(types.intp)
    zero = Constant(intp_t, 0)

    def get_first_item(seqty, seq):
        if isinstance(seqty, types.BaseTuple):
            if len(seqty) == 0:
                return None, None
            else:
                return seqty[0], builder.extract_value(seq, 0)
        else:
            getitem_impl = _get_borrowing_getitem(context, seqty)
            return seqty.dtype, getitem_impl(builder, (seq, zero))

    # Compute shape by traversing the first element of each nested
    # sequence
    shapes = []
    innerty, inner = seqty, seq

    for i in range(ndim):
        if i > 0:
            innerty, inner = get_first_item(innerty, inner)
        shapes.append(_get_seq_size(context, builder, innerty, inner))

    return tuple(shapes)


def check_sequence_shape(context, builder, seqty, seq, shapes):
    """
    Check the nested sequence matches the given *shapes*.
    """

    def _fail():
        context.call_conv.return_user_exc(builder, ValueError,
                                          ("incompatible sequence shape",))

    def check_seq_size(seqty, seq, shapes):
        if len(shapes) == 0:
            return

        size = _get_seq_size(context, builder, seqty, seq)
        expected = shapes[0]
        mismatch = builder.icmp_signed('!=', size, expected)
        with builder.if_then(mismatch, likely=False):
            _fail()

        if len(shapes) == 1:
            return

        if isinstance(seqty, types.Sequence):
            getitem_impl = _get_borrowing_getitem(context, seqty)
            with cgutils.for_range(builder, size) as loop:
                innerty = seqty.dtype
                inner = getitem_impl(builder, (seq, loop.index))
                check_seq_size(innerty, inner, shapes[1:])

        elif isinstance(seqty, types.BaseTuple):
            for i in range(len(seqty)):
                innerty = seqty[i]
                inner = builder.extract_value(seq, i)
                check_seq_size(innerty, inner, shapes[1:])

        else:
            assert 0, seqty

    check_seq_size(seqty, seq, shapes)


def assign_sequence_to_array(context, builder, data, shapes, strides,
                             arrty, seqty, seq):
    """
    Assign a nested sequence contents to an array.  The shape must match
    the sequence's structure.
    """

    def assign_item(indices, valty, val):
        ptr = cgutils.get_item_pointer2(context, builder, data, shapes, strides,
                                        arrty.layout, indices, wraparound=False)
        val = context.cast(builder, val, valty, arrty.dtype)
        store_item(context, builder, arrty, val, ptr)

    def assign(seqty, seq, shapes, indices):
        if len(shapes) == 0:
            assert not isinstance(seqty, (types.Sequence, types.BaseTuple))
            assign_item(indices, seqty, seq)
            return

        size = shapes[0]

        if isinstance(seqty, types.Sequence):
            getitem_impl = _get_borrowing_getitem(context, seqty)
            with cgutils.for_range(builder, size) as loop:
                innerty = seqty.dtype
                inner = getitem_impl(builder, (seq, loop.index))
                assign(innerty, inner, shapes[1:], indices + (loop.index,))

        elif isinstance(seqty, types.BaseTuple):
            for i in range(len(seqty)):
                innerty = seqty[i]
                inner = builder.extract_value(seq, i)
                index = context.get_constant(types.intp, i)
                assign(innerty, inner, shapes[1:], indices + (index,))

        else:
            assert 0, seqty

    assign(seqty, seq, shapes, ())


def np_array_typer(typingctx, object, dtype):
    ndim, seq_dtype = _parse_nested_sequence(typingctx, object)
    if is_nonelike(dtype):
        dtype = seq_dtype
    else:
        dtype = ty_parse_dtype(dtype)
        if dtype is None:
            return
    return types.Array(dtype, ndim, 'C')


@intrinsic
def np_array(typingctx, obj, dtype):
    _check_const_str_dtype("array", dtype)
    ret = np_array_typer(typingctx, obj, dtype)
    sig = ret(obj, dtype)

    def codegen(context, builder, sig, args):
        arrty = sig.return_type
        ndim = arrty.ndim
        seqty = sig.args[0]
        seq = args[0]

        shapes = compute_sequence_shape(context, builder, ndim, seqty, seq)
        assert len(shapes) == ndim

        check_sequence_shape(context, builder, seqty, seq, shapes)
        arr = _empty_nd_impl(context, builder, arrty, shapes)
        assign_sequence_to_array(context, builder, arr.data, shapes,
                                 arr.strides, arrty, seqty, seq)

        return impl_ret_new_ref(context, builder, sig.return_type,
                                arr._getvalue())

    return sig, codegen


@overload(np.array)
def impl_np_array(object, dtype=None):
    _check_const_str_dtype("array", dtype)
    if not type_can_asarray(object):
        raise errors.TypingError('The argument "object" must '
                                 'be array-like')
    if not is_nonelike(dtype) and ty_parse_dtype(dtype) is None:
        msg = 'The argument "dtype" must be a data-type if it is provided'
        raise errors.TypingError(msg)

    def impl(object, dtype=None):
        return np_array(object, dtype)
    return impl


def _normalize_axis(context, builder, func_name, ndim, axis):
    zero = axis.type(0)
    ll_ndim = axis.type(ndim)

    # Normalize negative axis
    is_neg_axis = builder.icmp_signed('<', axis, zero)
    axis = builder.select(is_neg_axis, builder.add(axis, ll_ndim), axis)

    # Check axis for bounds
    axis_out_of_bounds = builder.or_(
        builder.icmp_signed('<', axis, zero),
        builder.icmp_signed('>=', axis, ll_ndim))
    with builder.if_then(axis_out_of_bounds, likely=False):
        msg = "%s(): axis out of bounds" % func_name
        context.call_conv.return_user_exc(builder, IndexError, (msg,))

    return axis


def _insert_axis_in_shape(context, builder, orig_shape, ndim, axis):
    """
    Compute shape with the new axis inserted
    e.g. given original shape (2, 3, 4) and axis=2,
    the returned new shape is (2, 3, 1, 4).
    """
    assert len(orig_shape) == ndim - 1

    ll_shty = ir.ArrayType(cgutils.intp_t, ndim)
    shapes = cgutils.alloca_once(builder, ll_shty)

    one = cgutils.intp_t(1)

    # 1. copy original sizes at appropriate places
    for dim in range(ndim - 1):
        ll_dim = cgutils.intp_t(dim)
        after_axis = builder.icmp_signed('>=', ll_dim, axis)
        sh = orig_shape[dim]
        idx = builder.select(after_axis,
                             builder.add(ll_dim, one),
                             ll_dim)
        builder.store(sh, cgutils.gep_inbounds(builder, shapes, 0, idx))

    # 2. insert new size (1) at axis dimension
    builder.store(one, cgutils.gep_inbounds(builder, shapes, 0, axis))

    return cgutils.unpack_tuple(builder, builder.load(shapes))


def _insert_axis_in_strides(context, builder, orig_strides, ndim, axis):
    """
    Same as _insert_axis_in_shape(), but with a strides array.
    """
    assert len(orig_strides) == ndim - 1

    ll_shty = ir.ArrayType(cgutils.intp_t, ndim)
    strides = cgutils.alloca_once(builder, ll_shty)

    one = cgutils.intp_t(1)
    zero = cgutils.intp_t(0)

    # 1. copy original strides at appropriate places
    for dim in range(ndim - 1):
        ll_dim = cgutils.intp_t(dim)
        after_axis = builder.icmp_signed('>=', ll_dim, axis)
        idx = builder.select(after_axis,
                             builder.add(ll_dim, one),
                             ll_dim)
        builder.store(orig_strides[dim],
                      cgutils.gep_inbounds(builder, strides, 0, idx))

    # 2. insert new stride at axis dimension
    # (the value is indifferent for a 1-sized dimension, we use 0)
    builder.store(zero, cgutils.gep_inbounds(builder, strides, 0, axis))

    return cgutils.unpack_tuple(builder, builder.load(strides))


def expand_dims(context, builder, sig, args, axis):
    """
    np.expand_dims() with the given axis.
    """
    retty = sig.return_type
    ndim = retty.ndim
    arrty = sig.args[0]

    arr = make_array(arrty)(context, builder, value=args[0])
    ret = make_array(retty)(context, builder)

    shapes = cgutils.unpack_tuple(builder, arr.shape)
    strides = cgutils.unpack_tuple(builder, arr.strides)

    new_shapes = _insert_axis_in_shape(context, builder, shapes, ndim, axis)
    new_strides = _insert_axis_in_strides(context, builder, strides, ndim, axis)

    populate_array(ret,
                   data=arr.data,
                   shape=new_shapes,
                   strides=new_strides,
                   itemsize=arr.itemsize,
                   meminfo=arr.meminfo,
                   parent=arr.parent)

    return ret._getvalue()


@intrinsic
def np_expand_dims(typingctx, a, axis):
    layout = a.layout if a.ndim <= 1 else 'A'
    ret = a.copy(ndim=a.ndim + 1, layout=layout)
    sig = ret(a, axis)

    def codegen(context, builder, sig, args):
        axis = context.cast(builder, args[1], sig.args[1], types.intp)
        axis = _normalize_axis(context, builder, "np.expand_dims",
                               sig.return_type.ndim, axis)

        ret = expand_dims(context, builder, sig, args, axis)
        return impl_ret_borrowed(context, builder, sig.return_type, ret)

    return sig, codegen


@overload(np.expand_dims)
def impl_np_expand_dims(a, axis):
    if not isinstance(a, types.Array):
        msg = f'First argument "a" must be an array. Got {a}'
        raise errors.TypingError(msg)

    if not isinstance(axis, types.Integer):
        msg = f'Argument "axis" must be an integer. Got {axis}'
        raise errors.TypingError(msg)

    def impl(a, axis):
        return np_expand_dims(a, axis)
    return impl


def _atleast_nd(minimum, axes):
    @intrinsic
    def impl(typingcontext, *args):
        arrtys = args
        rettys = [arg.copy(ndim=max(arg.ndim, minimum)) for arg in args]

        def codegen(context, builder, sig, args):
            transform = _atleast_nd_transform(minimum, axes)
            arrs = cgutils.unpack_tuple(builder, args[0])

            rets = [transform(context, builder, arr, arrty, retty)
                    for arr, arrty, retty in zip(arrs, arrtys, rettys)]

            if len(rets) > 1:
                ret = context.make_tuple(builder, sig.return_type, rets)
            else:
                ret = rets[0]
            return impl_ret_borrowed(context, builder, sig.return_type, ret)

        return signature(types.Tuple(rettys) if len(rettys) > 1 else rettys[0],
                         types.StarArgTuple.from_types(args)), codegen

    return lambda *args: impl(*args)


def _atleast_nd_transform(min_ndim, axes):
    """
    Return a callback successively inserting 1-sized dimensions at the
    following axes.
    """
    assert min_ndim == len(axes)

    def transform(context, builder, arr, arrty, retty):
        for i in range(min_ndim):
            ndim = i + 1
            if arrty.ndim < ndim:
                axis = cgutils.intp_t(axes[i])
                newarrty = arrty.copy(ndim=arrty.ndim + 1)
                arr = expand_dims(context, builder,
                                  typing.signature(newarrty, arrty), (arr,),
                                  axis)
                arrty = newarrty

        return arr

    return transform


@overload(np.atleast_1d)
def np_atleast_1d(*args):
    if all(isinstance(arg, types.Array) for arg in args):
        return _atleast_nd(1, [0])


@overload(np.atleast_2d)
def np_atleast_2d(*args):
    if all(isinstance(arg, types.Array) for arg in args):
        return _atleast_nd(2, [0, 0])


@overload(np.atleast_3d)
def np_atleast_3d(*args):
    if all(isinstance(arg, types.Array) for arg in args):
        return _atleast_nd(3, [0, 0, 2])


def _do_concatenate(context, builder, axis,
                    arrtys, arrs, arr_shapes, arr_strides,
                    retty, ret_shapes):
    """
    Concatenate arrays along the given axis.
    """
    assert len(arrtys) == len(arrs) == len(arr_shapes) == len(arr_strides)

    zero = cgutils.intp_t(0)

    # Allocate return array
    ret = _empty_nd_impl(context, builder, retty, ret_shapes)
    ret_strides = cgutils.unpack_tuple(builder, ret.strides)

    # Compute the offset by which to bump the destination pointer
    # after copying each input array.
    # Morally, we need to copy each input array at different start indices
    # into the destination array; bumping the destination pointer
    # is simply easier than offsetting all destination indices.
    copy_offsets = []

    for arr_sh in arr_shapes:
        # offset = ret_strides[axis] * input_shape[axis]
        offset = zero
        for dim, (size, stride) in enumerate(zip(arr_sh, ret_strides)):
            is_axis = builder.icmp_signed('==', axis.type(dim), axis)
            addend = builder.mul(size, stride)
            offset = builder.select(is_axis,
                                    builder.add(offset, addend),
                                    offset)
        copy_offsets.append(offset)

    # Copy input arrays into the return array
    ret_data = ret.data

    for arrty, arr, arr_sh, arr_st, offset in zip(arrtys, arrs, arr_shapes,
                                                  arr_strides, copy_offsets):
        arr_data = arr.data

        # Do the copy loop
        # Note the loop nesting is optimized for the destination layout
        loop_nest = cgutils.loop_nest(builder, arr_sh, cgutils.intp_t,
                                      order=retty.layout)

        with loop_nest as indices:
            src_ptr = cgutils.get_item_pointer2(context, builder, arr_data,
                                                arr_sh, arr_st,
                                                arrty.layout, indices)
            val = load_item(context, builder, arrty, src_ptr)
            val = context.cast(builder, val, arrty.dtype, retty.dtype)
            dest_ptr = cgutils.get_item_pointer2(context, builder, ret_data,
                                                 ret_shapes, ret_strides,
                                                 retty.layout, indices)
            store_item(context, builder, retty, val, dest_ptr)

        # Bump destination pointer
        ret_data = cgutils.pointer_add(builder, ret_data, offset)

    return ret


def _np_concatenate(context, builder, arrtys, arrs, retty, axis):
    ndim = retty.ndim

    arrs = [make_array(aty)(context, builder, value=a)
            for aty, a in zip(arrtys, arrs)]

    axis = _normalize_axis(context, builder, "np.concatenate", ndim, axis)

    # Get input shapes
    arr_shapes = [cgutils.unpack_tuple(builder, arr.shape) for arr in arrs]
    arr_strides = [cgutils.unpack_tuple(builder, arr.strides) for arr in arrs]

    # Compute return shape:
    # - the dimension for the concatenation axis is summed over all inputs
    # - other dimensions must match exactly for each input
    ret_shapes = [cgutils.alloca_once_value(builder, sh)
                  for sh in arr_shapes[0]]

    for dim in range(ndim):
        is_axis = builder.icmp_signed('==', axis.type(dim), axis)
        ret_shape_ptr = ret_shapes[dim]
        ret_sh = builder.load(ret_shape_ptr)
        other_shapes = [sh[dim] for sh in arr_shapes[1:]]

        with builder.if_else(is_axis) as (on_axis, on_other_dim):
            with on_axis:
                sh = functools.reduce(
                    builder.add,
                    other_shapes + [ret_sh])
                builder.store(sh, ret_shape_ptr)

            with on_other_dim:
                is_ok = cgutils.true_bit
                for sh in other_shapes:
                    is_ok = builder.and_(is_ok,
                                         builder.icmp_signed('==', sh, ret_sh))
                with builder.if_then(builder.not_(is_ok), likely=False):
                    context.call_conv.return_user_exc(
                        builder, ValueError,
                        ("np.concatenate(): input sizes over "
                         "dimension %d do not match" % dim,))

    ret_shapes = [builder.load(sh) for sh in ret_shapes]

    ret = _do_concatenate(context, builder, axis,
                          arrtys, arrs, arr_shapes, arr_strides,
                          retty, ret_shapes)
    return impl_ret_new_ref(context, builder, retty, ret._getvalue())


def _np_stack(context, builder, arrtys, arrs, retty, axis):
    ndim = retty.ndim

    zero = cgutils.intp_t(0)
    one = cgutils.intp_t(1)
    ll_narrays = cgutils.intp_t(len(arrs))

    arrs = [make_array(aty)(context, builder, value=a)
            for aty, a in zip(arrtys, arrs)]

    axis = _normalize_axis(context, builder, "np.stack", ndim, axis)

    # Check input arrays have the same shape
    orig_shape = cgutils.unpack_tuple(builder, arrs[0].shape)

    for arr in arrs[1:]:
        is_ok = cgutils.true_bit
        for sh, orig_sh in zip(cgutils.unpack_tuple(builder, arr.shape),
                               orig_shape):
            is_ok = builder.and_(is_ok, builder.icmp_signed('==', sh, orig_sh))
            with builder.if_then(builder.not_(is_ok), likely=False):
                context.call_conv.return_user_exc(
                    builder, ValueError,
                    ("np.stack(): all input arrays must have the same shape",))

    orig_strides = [cgutils.unpack_tuple(builder, arr.strides) for arr in arrs]

    # Compute input shapes and return shape with the new axis inserted
    # e.g. given 5 input arrays of shape (2, 3, 4) and axis=1,
    # corrected input shape is (2, 1, 3, 4) and return shape is (2, 5, 3, 4).
    ll_shty = ir.ArrayType(cgutils.intp_t, ndim)

    input_shapes = cgutils.alloca_once(builder, ll_shty)
    ret_shapes = cgutils.alloca_once(builder, ll_shty)

    # 1. copy original sizes at appropriate places
    for dim in range(ndim - 1):
        ll_dim = cgutils.intp_t(dim)
        after_axis = builder.icmp_signed('>=', ll_dim, axis)
        sh = orig_shape[dim]
        idx = builder.select(after_axis,
                             builder.add(ll_dim, one),
                             ll_dim)
        builder.store(sh, cgutils.gep_inbounds(builder, input_shapes, 0, idx))
        builder.store(sh, cgutils.gep_inbounds(builder, ret_shapes, 0, idx))

    # 2. insert new size at axis dimension
    builder.store(one, cgutils.gep_inbounds(builder, input_shapes, 0, axis))
    builder.store(ll_narrays, cgutils.gep_inbounds(builder,
                                                   ret_shapes,
                                                   0,
                                                   axis))

    input_shapes = cgutils.unpack_tuple(builder, builder.load(input_shapes))
    input_shapes = [input_shapes] * len(arrs)
    ret_shapes = cgutils.unpack_tuple(builder, builder.load(ret_shapes))

    # Compute input strides for each array with the new axis inserted
    input_strides = [cgutils.alloca_once(builder, ll_shty)
                     for i in range(len(arrs))]

    # 1. copy original strides at appropriate places
    for dim in range(ndim - 1):
        ll_dim = cgutils.intp_t(dim)
        after_axis = builder.icmp_signed('>=', ll_dim, axis)
        idx = builder.select(after_axis,
                             builder.add(ll_dim, one),
                             ll_dim)
        for i in range(len(arrs)):
            builder.store(orig_strides[i][dim],
                          cgutils.gep_inbounds(builder, input_strides[i], 0,
                                               idx))

    # 2. insert new stride at axis dimension
    # (the value is indifferent for a 1-sized dimension, we put 0)
    for i in range(len(arrs)):
        builder.store(zero, cgutils.gep_inbounds(builder, input_strides[i], 0,
                                                 axis))

    input_strides = [cgutils.unpack_tuple(builder, builder.load(st))
                     for st in input_strides]

    # Create concatenated array
    ret = _do_concatenate(context, builder, axis,
                          arrtys, arrs, input_shapes, input_strides,
                          retty, ret_shapes)
    return impl_ret_new_ref(context, builder, retty, ret._getvalue())


def np_concatenate_typer(typingctx, arrays, axis):
    if axis is not None and not isinstance(axis, types.Integer):
        # Note Numpy allows axis=None, but it isn't documented:
        # https://github.com/numpy/numpy/issues/7968
        return

    # does type checking
    dtype, ndim = _sequence_of_arrays(typingctx,
                                      "np.concatenate", arrays)
    if ndim == 0:
        raise TypeError("zero-dimensional arrays cannot be concatenated")

    layout = _choose_concatenation_layout(arrays)

    return types.Array(dtype, ndim, layout)


@intrinsic
def np_concatenate(typingctx, arrays, axis):
    ret = np_concatenate_typer(typingctx, arrays, axis)
    assert isinstance(ret, types.Array)
    sig = ret(arrays, axis)

    def codegen(context, builder, sig, args):
        axis = context.cast(builder, args[1], sig.args[1], types.intp)
        return _np_concatenate(context, builder,
                               list(sig.args[0]),
                               cgutils.unpack_tuple(builder, args[0]),
                               sig.return_type,
                               axis)

    return sig, codegen


@overload(np.concatenate)
def impl_np_concatenate(arrays, axis=0):
    if isinstance(arrays, types.BaseTuple):
        def impl(arrays, axis=0):
            return np_concatenate(arrays, axis)
        return impl


def _column_stack_dims(context, func_name, arrays):
    # column_stack() allows stacking 1-d and 2-d arrays together
    for a in arrays:
        if a.ndim < 1 or a.ndim > 2:
            raise TypeError("np.column_stack() is only defined on "
                            "1-d and 2-d arrays")
    return 2


@intrinsic
def np_column_stack(typingctx, tup):
    dtype, ndim = _sequence_of_arrays(typingctx,
                                      "np.column_stack", tup,
                                      dim_chooser=_column_stack_dims)
    layout = _choose_concatenation_layout(tup)
    ret = types.Array(dtype, ndim, layout)
    sig = ret(tup)

    def codegen(context, builder, sig, args):
        orig_arrtys = list(sig.args[0])
        orig_arrs = cgutils.unpack_tuple(builder, args[0])

        arrtys = []
        arrs = []

        axis = context.get_constant(types.intp, 1)

        for arrty, arr in zip(orig_arrtys, orig_arrs):
            if arrty.ndim == 2:
                arrtys.append(arrty)
                arrs.append(arr)
            else:
                # Convert 1d array to 2d column array: np.expand_dims(a, 1)
                assert arrty.ndim == 1
                newty = arrty.copy(ndim=2)
                expand_sig = typing.signature(newty, arrty)
                newarr = expand_dims(context, builder, expand_sig, (arr,), axis)

                arrtys.append(newty)
                arrs.append(newarr)

        return _np_concatenate(context, builder, arrtys, arrs,
                               sig.return_type, axis)

    return sig, codegen


@overload(np.column_stack)
def impl_column_stack(tup):
    if isinstance(tup, types.BaseTuple):
        def impl(tup):
            return np_column_stack(tup)
        return impl


def _np_stack_common(context, builder, sig, args, axis):
    """
    np.stack() with the given axis value.
    """
    return _np_stack(context, builder,
                     list(sig.args[0]),
                     cgutils.unpack_tuple(builder, args[0]),
                     sig.return_type,
                     axis)


@intrinsic
def np_stack_common(typingctx, arrays, axis):
    # does type checking
    dtype, ndim = _sequence_of_arrays(typingctx,
                                      "np.stack", arrays)
    layout = 'F' if all(a.layout == 'F' for a in arrays) else 'C'
    ret = types.Array(dtype, ndim + 1, layout)
    sig = ret(arrays, axis)

    def codegen(context, builder, sig, args):
        axis = context.cast(builder, args[1], sig.args[1], types.intp)
        return _np_stack_common(context, builder, sig, args, axis)

    return sig, codegen


@overload(np.stack)
def impl_np_stack(arrays, axis=0):
    if isinstance(arrays, types.BaseTuple):
        def impl(arrays, axis=0):
            return np_stack_common(arrays, axis)
        return impl


def NdStack_typer(typingctx, func_name, arrays, ndim_min):
    # does type checking
    dtype, ndim = _sequence_of_arrays(typingctx, func_name, arrays)
    ndim = max(ndim, ndim_min)
    layout = _choose_concatenation_layout(arrays)
    ret = types.Array(dtype, ndim, layout)
    return ret


@intrinsic
def _np_hstack(typingctx, tup):
    ret = NdStack_typer(typingctx, "np.hstack", tup, 1)
    sig = ret(tup)

    def codegen(context, builder, sig, args):
        tupty = sig.args[0]
        ndim = tupty[0].ndim

        if ndim == 0:
            # hstack() on 0-d arrays returns a 1-d array
            axis = context.get_constant(types.intp, 0)
            return _np_stack_common(context, builder, sig, args, axis)

        else:
            # As a special case, dimension 0 of 1-dimensional arrays
            # is "horizontal"
            axis = 0 if ndim == 1 else 1

            def np_hstack_impl(arrays):
                return np.concatenate(arrays, axis=axis)

            return context.compile_internal(builder, np_hstack_impl, sig, args)

    return sig, codegen


@overload(np.hstack)
def impl_np_hstack(tup):
    if isinstance(tup, types.BaseTuple):
        def impl(tup):
            return _np_hstack(tup)
        return impl


@intrinsic
def _np_vstack(typingctx, tup):
    ret = NdStack_typer(typingctx, "np.vstack", tup, 2)
    sig = ret(tup)

    def codegen(context, builder, sig, args):
        tupty = sig.args[0]
        ndim = tupty[0].ndim

        if ndim == 0:
            def np_vstack_impl(arrays):
                return np.expand_dims(np.hstack(arrays), 1)

        elif ndim == 1:
            # np.stack(arrays, axis=0)
            axis = context.get_constant(types.intp, 0)
            return _np_stack_common(context, builder, sig, args, axis)

        else:
            def np_vstack_impl(arrays):
                return np.concatenate(arrays, axis=0)

        return context.compile_internal(builder, np_vstack_impl, sig, args)

    return sig, codegen


@overload(np.vstack)
def impl_np_vstack(tup):
    if isinstance(tup, types.BaseTuple):
        def impl(tup):
            return _np_vstack(tup)
        return impl


if numpy_version >= (2, 0):
    overload(np.row_stack)(impl_np_vstack)


@intrinsic
def _np_dstack(typingctx, tup):
    ret = NdStack_typer(typingctx, "np.dstack", tup, 3)
    sig = ret(tup)

    def codegen(context, builder, sig, args):
        tupty = sig.args[0]
        retty = sig.return_type
        ndim = tupty[0].ndim

        if ndim == 0:
            def np_vstack_impl(arrays):
                return np.hstack(arrays).reshape(1, 1, -1)

            return context.compile_internal(builder, np_vstack_impl, sig, args)

        elif ndim == 1:
            # np.expand_dims(np.stack(arrays, axis=1), axis=0)
            axis = context.get_constant(types.intp, 1)
            stack_retty = retty.copy(ndim=retty.ndim - 1)
            stack_sig = typing.signature(stack_retty, *sig.args)
            stack_ret = _np_stack_common(context, builder, stack_sig, args,
                                         axis)

            axis = context.get_constant(types.intp, 0)
            expand_sig = typing.signature(retty, stack_retty)
            return expand_dims(context, builder, expand_sig, (stack_ret,), axis)

        elif ndim == 2:
            # np.stack(arrays, axis=2)
            axis = context.get_constant(types.intp, 2)
            return _np_stack_common(context, builder, sig, args, axis)

        else:
            def np_vstack_impl(arrays):
                return np.concatenate(arrays, axis=2)

            return context.compile_internal(builder, np_vstack_impl, sig, args)

    return sig, codegen


@overload(np.dstack)
def impl_np_dstack(tup):
    if isinstance(tup, types.BaseTuple):
        def impl(tup):
            return _np_dstack(tup)
        return impl


@extending.overload_method(types.Array, 'fill')
def arr_fill(arr, val):

    def fill_impl(arr, val):
        arr[:] = val
        return None

    return fill_impl


@extending.overload_method(types.Array, 'dot')
def array_dot(arr, other):
    def dot_impl(arr, other):
        return np.dot(arr, other)

    return dot_impl


@overload(np.fliplr)
def np_flip_lr(m):

    if not type_can_asarray(m):
        raise errors.TypingError("Cannot np.fliplr on %s type" % m)

    def impl(m):
        A = np.asarray(m)
        # this handling is superfluous/dead as < 2d array cannot be indexed as
        # present below and so typing fails. If the typing doesn't fail due to
        # some future change, this will catch it.
        if A.ndim < 2:
            raise ValueError('Input must be >= 2-d.')
        return A[::, ::-1, ...]
    return impl


@overload(np.flipud)
def np_flip_ud(m):

    if not type_can_asarray(m):
        raise errors.TypingError("Cannot np.flipud on %s type" % m)

    def impl(m):
        A = np.asarray(m)
        # this handling is superfluous/dead as a 0d array cannot be indexed as
        # present below and so typing fails. If the typing doesn't fail due to
        # some future change, this will catch it.
        if A.ndim < 1:
            raise ValueError('Input must be >= 1-d.')
        return A[::-1, ...]
    return impl


@intrinsic
def _build_flip_slice_tuple(tyctx, sz):
    """ Creates a tuple of slices for np.flip indexing like
    `(slice(None, None, -1),) * sz` """
    if not isinstance(sz, types.IntegerLiteral):
        raise errors.RequireLiteralValue(sz)
    size = int(sz.literal_value)
    tuple_type = types.UniTuple(dtype=types.slice3_type, count=size)
    sig = tuple_type(sz)

    def codegen(context, builder, signature, args):
        def impl(length, empty_tuple):
            out = empty_tuple
            for i in range(length):
                out = tuple_setitem(out, i, slice(None, None, -1))
            return out

        inner_argtypes = [types.intp, tuple_type]
        inner_sig = typing.signature(tuple_type, *inner_argtypes)
        ll_idx_type = context.get_value_type(types.intp)
        # Allocate an empty tuple
        empty_tuple = context.get_constant_undef(tuple_type)
        inner_args = [ll_idx_type(size), empty_tuple]

        res = context.compile_internal(builder, impl, inner_sig, inner_args)
        return res

    return sig, codegen


@overload(np.flip)
def np_flip(m):
    # a constant value is needed for the tuple slice, types.Array.ndim can
    # provide this and so at presnet only type.Array is support
    if not isinstance(m, types.Array):
        raise errors.TypingError("Cannot np.flip on %s type" % m)

    def impl(m):
        sl = _build_flip_slice_tuple(m.ndim)
        return m[sl]

    return impl


@overload(np.array_split)
def np_array_split(ary, indices_or_sections, axis=0):
    if isinstance(ary, (types.UniTuple, types.ListType, types.List)):
        def impl(ary, indices_or_sections, axis=0):
            return np.array_split(
                np.asarray(ary),
                indices_or_sections,
                axis=axis
            )

        return impl

    if isinstance(indices_or_sections, types.Integer):
        def impl(ary, indices_or_sections, axis=0):
            l, rem = divmod(ary.shape[axis], indices_or_sections)
            indices = np.cumsum(np.array(
                [l + 1] * rem +
                [l] * (indices_or_sections - rem - 1)
            ))
            return np.array_split(ary, indices, axis=axis)

        return impl

    elif (
        isinstance(indices_or_sections, types.IterableType)
        and isinstance(
            indices_or_sections.iterator_type.yield_type,
            types.Integer
        )
    ):
        def impl(ary, indices_or_sections, axis=0):
            slice_tup = build_full_slice_tuple(ary.ndim)
            axis = normalize_axis("np.split", "axis", ary.ndim, axis)
            out = []
            prev = 0
            for cur in indices_or_sections:
                idx = tuple_setitem(slice_tup, axis, slice(prev, cur))
                out.append(ary[idx])
                prev = cur
            out.append(ary[tuple_setitem(slice_tup, axis, slice(cur, None))])
            return out

        return impl

    elif (
        isinstance(indices_or_sections, types.Tuple)
        and all(isinstance(t, types.Integer) for t in indices_or_sections.types)
    ):
        def impl(ary, indices_or_sections, axis=0):
            slice_tup = build_full_slice_tuple(ary.ndim)
            axis = normalize_axis("np.split", "axis", ary.ndim, axis)
            out = []
            prev = 0
            for cur in literal_unroll(indices_or_sections):
                idx = tuple_setitem(slice_tup, axis, slice(prev, cur))
                out.append(ary[idx])
                prev = cur
            out.append(ary[tuple_setitem(slice_tup, axis, slice(cur, None))])
            return out

        return impl


@overload(np.split)
def np_split(ary, indices_or_sections, axis=0):
    # This is just a wrapper of array_split, but with an extra error if
    # indices is an int.
    if isinstance(ary, (types.UniTuple, types.ListType, types.List)):
        def impl(ary, indices_or_sections, axis=0):
            return np.split(np.asarray(ary), indices_or_sections, axis=axis)

        return impl

    if isinstance(indices_or_sections, types.Integer):
        def impl(ary, indices_or_sections, axis=0):
            _, rem = divmod(ary.shape[axis], indices_or_sections)
            if rem != 0:
                raise ValueError(
                    "array split does not result in an equal division"
                )
            return np.array_split(
                ary, indices_or_sections, axis=axis
            )

        return impl

    else:
        return np_array_split(ary, indices_or_sections, axis=axis)


@overload(np.vsplit)
def numpy_vsplit(ary, indices_or_sections):
    if not isinstance(ary, types.Array):
        msg = 'The argument "ary" must be an array'
        raise errors.TypingError(msg)

    if not isinstance(indices_or_sections, (types.Integer, types.Array,
                                            types.List, types.UniTuple)):
        msg = ('The argument "indices_or_sections" must be int or 1d-array')
        raise errors.TypingError(msg)

    def impl(ary, indices_or_sections):
        if ary.ndim < 2:
            raise ValueError(('vsplit only works on '
                              'arrays of 2 or more dimensions'))
        return np.split(ary, indices_or_sections, axis=0)

    return impl


@overload(np.hsplit)
def numpy_hsplit(ary, indices_or_sections):
    if not isinstance(ary, types.Array):
        msg = 'The argument "ary" must be an array'
        raise errors.TypingError(msg)

    if not isinstance(indices_or_sections, (types.Integer, types.Array,
                                            types.List, types.UniTuple)):
        msg = ('The argument "indices_or_sections" must be int or 1d-array')
        raise errors.TypingError(msg)

    def impl(ary, indices_or_sections):
        if ary.ndim == 0:
            raise ValueError(('hsplit only works on '
                              'arrays of 1 or more dimensions'))
        if ary.ndim > 1:
            return np.split(ary, indices_or_sections, axis=1)
        return np.split(ary, indices_or_sections, axis=0)

    return impl


@overload(np.dsplit)
def numpy_dsplit(ary, indices_or_sections):
    if not isinstance(ary, types.Array):
        msg = 'The argument "ary" must be an array'
        raise errors.TypingError(msg)

    if not isinstance(indices_or_sections, (types.Integer, types.Array,
                                            types.List, types.UniTuple)):
        msg = ('The argument "indices_or_sections" must be int or 1d-array')
        raise errors.TypingError(msg)

    def impl(ary, indices_or_sections):
        if ary.ndim < 3:
            raise ValueError('dsplit only works on arrays of 3 or more '
                             'dimensions')
        return np.split(ary, indices_or_sections, axis=2)

    return impl


# -----------------------------------------------------------------------------
# Sorting

_sorts = {}


def default_lt(a, b):
    """
    Trivial comparison function between two keys.
    """
    return a < b


def get_sort_func(kind, lt_impl, is_argsort=False):
    """
    Get a sort implementation of the given kind.
    """
    key = kind, lt_impl.__name__, is_argsort

    try:
        return _sorts[key]
    except KeyError:
        if kind == 'quicksort':
            sort = quicksort.make_jit_quicksort(
                lt=lt_impl,
                is_argsort=is_argsort,
                is_np_array=True)
            func = sort.run_quicksort
        elif kind == 'mergesort':
            sort = mergesort.make_jit_mergesort(
                lt=lt_impl,
                is_argsort=is_argsort)
            func = sort.run_mergesort
        _sorts[key] = func
        return func


def lt_implementation(dtype):
    if isinstance(dtype, types.Float):
        return lt_floats
    elif isinstance(dtype, types.Complex):
        return lt_complex
    else:
        return default_lt


@lower_builtin("array.sort", types.Array)
def array_sort(context, builder, sig, args):
    arytype = sig.args[0]

    sort_func = get_sort_func(kind='quicksort',
                              lt_impl=lt_implementation(arytype.dtype))

    def array_sort_impl(arr):
        # Note we clobber the return value
        sort_func(arr)

    return context.compile_internal(builder, array_sort_impl, sig, args)


@overload(np.sort)
def impl_np_sort(a):
    if not type_can_asarray(a):
        raise errors.TypingError('Argument "a" must '
                                 'be array-like')

    def np_sort_impl(a):
        res = a.copy()
        res.sort()
        return res
    return np_sort_impl


@lower_builtin("array.argsort", types.Array, types.StringLiteral)
@lower_builtin(np.argsort, types.Array, types.StringLiteral)
def array_argsort(context, builder, sig, args):
    arytype, kind = sig.args

    sort_func = get_sort_func(kind=kind.literal_value,
                              lt_impl=lt_implementation(arytype.dtype),
                              is_argsort=True)

    def array_argsort_impl(arr):
        return sort_func(arr)

    innersig = sig.replace(args=sig.args[:1])
    innerargs = args[:1]
    return context.compile_internal(builder, array_argsort_impl,
                                    innersig, innerargs)


# ------------------------------------------------------------------------------
# Implicit cast

@lower_cast(types.Array, types.Array)
def array_to_array(context, builder, fromty, toty, val):
    # Type inference should have prevented illegal array casting.
    assert fromty.mutable != toty.mutable or toty.layout == 'A'
    return val


@lower_cast(types.Array, types.UnicodeCharSeq)
@lower_cast(types.Array, types.Float)
@lower_cast(types.Array, types.Integer)
@lower_cast(types.Array, types.Complex)
@lower_cast(types.Array, types.Boolean)
@lower_cast(types.Array, types.NPTimedelta)
@lower_cast(types.Array, types.NPDatetime)
def array0d_to_scalar(context, builder, fromty, toty, val):
    def impl(a):
        # a is an array(T, 0d, O), T is type, O is order
        return a.take(0)

    sig = signature(toty, fromty)
    res = context.compile_internal(builder, impl, sig, [val])
    return impl_ret_untracked(context, builder, sig.return_type, res)


@lower_cast(types.Array, types.UnicodeCharSeq)
def array_to_unichrseq(context, builder, fromty, toty, val):
    def impl(a):
        return str(a[()])

    sig = signature(toty, fromty)
    res = context.compile_internal(builder, impl, sig, [val])
    return impl_ret_borrowed(context, builder, sig.return_type, res)


# ------------------------------------------------------------------------------
# Stride tricks

def reshape_unchecked(a, shape, strides):
    """
    An intrinsic returning a derived array with the given shape and strides.
    """
    raise NotImplementedError


@extending.type_callable(reshape_unchecked)
def type_reshape_unchecked(context):
    def check_shape(shape):
        return (isinstance(shape, types.BaseTuple) and
                all(isinstance(v, types.Integer) for v in shape))

    def typer(a, shape, strides):
        if not isinstance(a, types.Array):
            return
        if not check_shape(shape) or not check_shape(strides):
            return
        if len(shape) != len(strides):
            return
        return a.copy(ndim=len(shape), layout='A')

    return typer


@lower_builtin(reshape_unchecked, types.Array, types.BaseTuple, types.BaseTuple)
def impl_shape_unchecked(context, builder, sig, args):
    aryty = sig.args[0]
    retty = sig.return_type

    ary = make_array(aryty)(context, builder, args[0])
    out = make_array(retty)(context, builder)
    shape = cgutils.unpack_tuple(builder, args[1])
    strides = cgutils.unpack_tuple(builder, args[2])

    populate_array(out,
                   data=ary.data,
                   shape=shape,
                   strides=strides,
                   itemsize=ary.itemsize,
                   meminfo=ary.meminfo,
                   )

    res = out._getvalue()
    return impl_ret_borrowed(context, builder, retty, res)


@extending.overload(np.lib.stride_tricks.as_strided)
def as_strided(x, shape=None, strides=None):
    if shape in (None, types.none):
        @register_jitable
        def get_shape(x, shape):
            return x.shape
    else:
        @register_jitable
        def get_shape(x, shape):
            return shape

    if strides in (None, types.none):
        # When *strides* is not passed, as_strided() does a non-size-checking
        # reshape(), possibly changing the original strides.  This is too
        # cumbersome to support right now, and a Web search shows all example
        # use cases of as_strided() pass explicit *strides*.
        raise NotImplementedError("as_strided() strides argument is mandatory")
    else:
        @register_jitable
        def get_strides(x, strides):
            return strides

    def as_strided_impl(x, shape=None, strides=None):
        x = reshape_unchecked(x, get_shape(x, shape), get_strides(x, strides))
        return x

    return as_strided_impl


@extending.overload(np.lib.stride_tricks.sliding_window_view)
def sliding_window_view(x, window_shape, axis=None):

    # Window shape must be given as either an integer or tuple of integers.
    # We also need to generate buffer tuples we can modify to contain the
    # final shape and strides (reshape_unchecked does not accept lists).
    if isinstance(window_shape, types.Integer):
        shape_buffer = tuple(range(x.ndim + 1))
        stride_buffer = tuple(range(x.ndim + 1))

        @register_jitable
        def get_window_shape(window_shape):
            return (window_shape,)

    elif (isinstance(window_shape, types.UniTuple) and
            isinstance(window_shape.dtype, types.Integer)):
        shape_buffer = tuple(range(x.ndim + len(window_shape)))
        stride_buffer = tuple(range(x.ndim + len(window_shape)))

        @register_jitable
        def get_window_shape(window_shape):
            return window_shape

    else:
        raise errors.TypingError(
            "window_shape must be an integer or tuple of integers"
        )

    # Axis must be integer, tuple of integers, or None for all axes.
    if is_nonelike(axis):
        @register_jitable
        def get_axis(window_shape, axis, ndim):
            return list(range(ndim))

    elif isinstance(axis, types.Integer):
        @register_jitable
        def get_axis(window_shape, axis, ndim):
            return [
                normalize_axis("sliding_window_view", "axis", ndim, axis)
            ]

    elif (isinstance(axis, types.UniTuple) and
            isinstance(axis.dtype, types.Integer)):
        @register_jitable
        def get_axis(window_shape, axis, ndim):
            return [normalize_axis("sliding_window_view", "axis", ndim, a)
                    for a in axis]

    else:
        raise errors.TypingError(
            "axis must be None, an integer or tuple of integers"
        )

    def sliding_window_view_impl(x, window_shape, axis=None):
        window_shape = get_window_shape(window_shape)
        axis = get_axis(window_shape, axis, x.ndim)
        if len(window_shape) != len(axis):
            raise ValueError(
                "Must provide matching length window_shape and axis"
            )

        # Initialise view details with shape and strides of x.
        out_shape = shape_buffer
        out_strides = stride_buffer
        for i in range(x.ndim):
            out_shape = tuple_setitem(out_shape, i, x.shape[i])
            out_strides = tuple_setitem(out_strides, i, x.strides[i])

        # Trim the dimensions being windowed and set the window shape and
        # strides. Note: the same axis can be windowed repeatedly.
        i = x.ndim
        for ax, dim in zip(axis, window_shape):
            if dim < 0:
                raise ValueError(
                    "`window_shape` cannot contain negative values"
                )
            if out_shape[ax] < dim:
                raise ValueError(
                    "window_shape cannot be larger than input array shape"
                )

            trimmed = out_shape[ax] - dim + 1
            out_shape = tuple_setitem(out_shape, ax, trimmed)
            out_shape = tuple_setitem(out_shape, i, dim)
            out_strides = tuple_setitem(out_strides, i, x.strides[ax])
            i += 1

        # The NumPy version calls as_strided, but our implementation of
        # as_strided is effectively a wrapper for reshape_unchecked.
        view = reshape_unchecked(x, out_shape, out_strides)
        return view

    return sliding_window_view_impl


@overload(bool)
def ol_bool(arr):
    if isinstance(arr, types.Array):
        def impl(arr):
            if arr.size == 0:
                return False # this is deprecated
            elif arr.size == 1:
                return bool(arr.take(0))
            else:
                msg = ("The truth value of an array with more than one element "
                       "is ambiguous. Use a.any() or a.all()")
                raise ValueError(msg)
        return impl


@overload(np.swapaxes)
def numpy_swapaxes(a, axis1, axis2):
    if not isinstance(axis1, (int, types.Integer)):
        raise errors.TypingError('The second argument "axis1" must be an '
                                 'integer')
    if not isinstance(axis2, (int, types.Integer)):
        raise errors.TypingError('The third argument "axis2" must be an '
                                 'integer')
    if not isinstance(a, types.Array):
        raise errors.TypingError('The first argument "a" must be an array')

    # create tuple list for transpose
    ndim = a.ndim
    axes_list = tuple(range(ndim))

    def impl(a, axis1, axis2):
        axis1 = normalize_axis("np.swapaxes", "axis1", ndim, axis1)
        axis2 = normalize_axis("np.swapaxes", "axis2", ndim, axis2)

        # to ensure tuple_setitem support of negative values
        if axis1 < 0:
            axis1 += ndim
        if axis2 < 0:
            axis2 += ndim

        axes_tuple = tuple_setitem(axes_list, axis1, axis2)
        axes_tuple = tuple_setitem(axes_tuple, axis2, axis1)
        return np.transpose(a, axes_tuple)

    return impl


@register_jitable
def _take_along_axis_impl(
        arr, indices, axis, Ni_orig, Nk_orig, indices_broadcast_shape
):
    # Based on example code in
    # https://github.com/numpy/numpy/blob/623bc1fae1d47df24e7f1e29321d0c0ba2771ce0/numpy/lib/shape_base.py#L90-L103
    # With addition of pre-broadcasting:
    # https://github.com/numpy/numpy/issues/19704

    # Wrap axis, it's used in tuple_setitem so must be (axis >= 0) to ensure
    # the GEP is in bounds.
    axis = normalize_axis("np.take_along_axis", "axis", arr.ndim, axis)

    # Broadcast the two arrays to matching shapes:
    arr_shape = list(arr.shape)
    arr_shape[axis] = 1
    for i, (d1, d2) in enumerate(zip(arr_shape, indices.shape)):
        if d1 == 1:
            new_val = d2
        elif d2 == 1:
            new_val = d1
        else:
            if d1 != d2:
                raise ValueError(
                    "`arr` and `indices` dimensions don't match"
                )
            new_val = d1
        indices_broadcast_shape = tuple_setitem(
            indices_broadcast_shape, i, new_val
        )
    arr_broadcast_shape = tuple_setitem(
        indices_broadcast_shape, axis, arr.shape[axis]
    )
    arr = np.broadcast_to(arr, arr_broadcast_shape)
    indices = np.broadcast_to(indices, indices_broadcast_shape)

    Ni = Ni_orig
    if len(Ni_orig) > 0:
        for i in range(len(Ni)):
            Ni = tuple_setitem(Ni, i, arr.shape[i])
    Nk = Nk_orig
    if len(Nk_orig) > 0:
        for i in range(len(Nk)):
            Nk = tuple_setitem(Nk, i, arr.shape[axis + 1 + i])

    J = indices.shape[axis]  # Need not equal M
    out = np.empty(Ni + (J,) + Nk, arr.dtype)

    np_s_ = (slice(None, None, None),)

    for ii in np.ndindex(Ni):
        for kk in np.ndindex(Nk):
            a_1d = arr[ii + np_s_ + kk]
            indices_1d = indices[ii + np_s_ + kk]
            out_1d = out[ii + np_s_ + kk]
            for j in range(J):
                out_1d[j] = a_1d[indices_1d[j]]
    return out


@overload(np.take_along_axis)
def arr_take_along_axis(arr, indices, axis):
    if not isinstance(arr, types.Array):
        raise errors.TypingError('The first argument "arr" must be an array')
    if not isinstance(indices, types.Array):
        raise errors.TypingError(
            'The second argument "indices" must be an array')
    if not isinstance(indices.dtype, types.Integer):
        raise errors.TypingError('The indices array must contain integers')
    if is_nonelike(axis):
        arr_ndim = 1
    else:
        arr_ndim = arr.ndim
    if arr_ndim != indices.ndim:
        # Matches NumPy error:
        raise errors.TypingError(
            "`indices` and `arr` must have the same number of dimensions"
        )

    indices_broadcast_shape = tuple(range(indices.ndim))
    if is_nonelike(axis):
        def take_along_axis_impl(arr, indices, axis):
            return _take_along_axis_impl(arr.flatten(), indices, 0, (), (),
                                         indices_broadcast_shape)
    else:
        check_is_integer(axis, "axis")
        if not isinstance(axis, types.IntegerLiteral):
            raise errors.NumbaValueError("axis must be a literal value")
        axis = axis.literal_value
        if axis < 0:
            axis = arr.ndim + axis

        if axis < 0 or axis >= arr.ndim:
            raise errors.NumbaValueError("axis is out of bounds")

        Ni = tuple(range(axis))
        Nk = tuple(range(axis + 1, arr.ndim))

        def take_along_axis_impl(arr, indices, axis):
            return _take_along_axis_impl(arr, indices, axis, Ni, Nk,
                                         indices_broadcast_shape)
    return take_along_axis_impl


@overload(np.nan_to_num)
def nan_to_num_impl(x, copy=True, nan=0.0):
    if isinstance(x, types.Number):
        if isinstance(x, types.Integer):
            # Integers do not have nans or infs
            def impl(x, copy=True, nan=0.0):
                return x

        elif isinstance(x, types.Float):
            def impl(x, copy=True, nan=0.0):
                if np.isnan(x):
                    return nan
                elif np.isneginf(x):
                    return np.finfo(type(x)).min
                elif np.isposinf(x):
                    return np.finfo(type(x)).max
                return x
        elif isinstance(x, types.Complex):
            def impl(x, copy=True, nan=0.0):
                r = np.nan_to_num(x.real, nan=nan)
                c = np.nan_to_num(x.imag, nan=nan)
                return complex(r, c)
        else:
            raise errors.TypingError(
                "Only Integer, Float, and Complex values are accepted"
            )

    elif type_can_asarray(x):
        if isinstance(x.dtype, types.Integer):
            # Integers do not have nans or infs
            def impl(x, copy=True, nan=0.0):
                return x
        elif isinstance(x.dtype, types.Float):
            def impl(x, copy=True, nan=0.0):
                min_inf = np.finfo(x.dtype).min
                max_inf = np.finfo(x.dtype).max

                x_ = np.asarray(x)
                output = np.copy(x_) if copy else x_

                output_flat = output.flat
                for i in range(output.size):
                    if np.isnan(output_flat[i]):
                        output_flat[i] = nan
                    elif np.isneginf(output_flat[i]):
                        output_flat[i] = min_inf
                    elif np.isposinf(output_flat[i]):
                        output_flat[i] = max_inf
                return output
        elif isinstance(x.dtype, types.Complex):
            def impl(x, copy=True, nan=0.0):
                x_ = np.asarray(x)
                output = np.copy(x_) if copy else x_

                np.nan_to_num(output.real, copy=False, nan=nan)
                np.nan_to_num(output.imag, copy=False, nan=nan)
                return output
        else:
            raise errors.TypingError(
                "Only Integer, Float, and Complex values are accepted"
            )
    else:
        raise errors.TypingError("The first argument must be a scalar or an "
                                 "array-like")
    return impl
