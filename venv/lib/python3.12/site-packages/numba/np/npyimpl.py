"""
Implementation of functions in the Numpy package.
"""


import math
import sys
import itertools
from collections import namedtuple

import llvmlite.ir as ir

import numpy as np
import operator

from numba.np import arrayobj, ufunc_db, numpy_support
from numba.np.ufunc.sigparse import parse_signature
from numba.core.imputils import (Registry, impl_ret_new_ref, force_error_model, impl_ret_borrowed)
from numba.core import typing, types, utils, cgutils, callconv, config
from numba.np.numpy_support import (
    ufunc_find_matching_loop, select_array_wrapper, from_dtype, _ufunc_loop_sig
)
from numba.np.arrayobj import _getitem_array_generic
from numba.core.typing import npydecl
from numba.core.extending import overload, intrinsic

from numba.core import errors

registry = Registry('npyimpl')


########################################################################

# In the way we generate code, ufuncs work with scalar as well as
# with array arguments. The following helper classes help dealing
# with scalar and array arguments in a regular way.
#
# In short, the classes provide a uniform interface. The interface
# handles the indexing of as many dimensions as the array may have.
# For scalars, all indexing is ignored and when the value is read,
# the scalar is returned. For arrays code for actual indexing is
# generated and reading performs the appropriate indirection.

class _ScalarIndexingHelper(object):
    def update_indices(self, loop_indices, name):
        pass

    def as_values(self):
        pass


class _ScalarHelper(object):
    """Helper class to handle scalar arguments (and result).
    Note that store_data is only used when generating code for
    a scalar ufunc and to write the output value.

    For loading, the value is directly used without having any
    kind of indexing nor memory backing it up. This is the use
    for input arguments.

    For storing, a variable is created in the stack where the
    value will be written.

    Note that it is not supported (as it is unneeded for our
    current use-cases) reading back a stored value. This class
    will always "load" the original value it got at its creation.
    """
    def __init__(self, ctxt, bld, val, ty):
        self.context = ctxt
        self.builder = bld
        self.val = val
        self.base_type = ty
        intpty = ctxt.get_value_type(types.intp)
        self.shape = [ir.Constant(intpty, 1)]

        lty = ctxt.get_data_type(ty) if ty != types.boolean else ir.IntType(1)
        self._ptr = cgutils.alloca_once(bld, lty)

    def create_iter_indices(self):
        return _ScalarIndexingHelper()

    def load_data(self, indices):
        return self.val

    def store_data(self, indices, val):
        self.builder.store(val, self._ptr)

    @property
    def return_val(self):
        return self.builder.load(self._ptr)


class _ArrayIndexingHelper(namedtuple('_ArrayIndexingHelper',
                                      ('array', 'indices'))):
    def update_indices(self, loop_indices, name):
        bld = self.array.builder
        intpty = self.array.context.get_value_type(types.intp)
        ONE = ir.Constant(ir.IntType(intpty.width), 1)

        # we are only interested in as many inner dimensions as dimensions
        # the indexed array has (the outer dimensions are broadcast, so
        # ignoring the outer indices produces the desired result.
        indices = loop_indices[len(loop_indices) - len(self.indices):]
        for src, dst, dim in zip(indices, self.indices, self.array.shape):
            cond = bld.icmp_unsigned('>', dim, ONE)
            with bld.if_then(cond):
                bld.store(src, dst)

    def as_values(self):
        """
        The indexing helper is built using alloca for each value, so it
        actually contains pointers to the actual indices to load. Note
        that update_indices assumes the same. This method returns the
        indices as values
        """
        bld = self.array.builder
        return [bld.load(index) for index in self.indices]


class _ArrayHelper(namedtuple('_ArrayHelper', ('context', 'builder',
                                               'shape', 'strides', 'data',
                                               'layout', 'base_type', 'ndim',
                                               'return_val'))):
    """Helper class to handle array arguments/result.
    It provides methods to generate code loading/storing specific
    items as well as support code for handling indices.
    """
    def create_iter_indices(self):
        intpty = self.context.get_value_type(types.intp)
        ZERO = ir.Constant(ir.IntType(intpty.width), 0)

        indices = []
        for i in range(self.ndim):
            x = cgutils.alloca_once(self.builder, ir.IntType(intpty.width))
            self.builder.store(ZERO, x)
            indices.append(x)
        return _ArrayIndexingHelper(self, indices)

    def _load_effective_address(self, indices):
        return cgutils.get_item_pointer2(self.context,
                                         self.builder,
                                         data=self.data,
                                         shape=self.shape,
                                         strides=self.strides,
                                         layout=self.layout,
                                         inds=indices)

    def load_data(self, indices):
        model = self.context.data_model_manager[self.base_type]
        ptr = self._load_effective_address(indices)
        return model.load_from_data_pointer(self.builder, ptr)

    def store_data(self, indices, value):
        ctx = self.context
        bld = self.builder
        store_value = ctx.get_value_as_data(bld, self.base_type, value)
        assert ctx.get_data_type(self.base_type) == store_value.type
        bld.store(store_value, self._load_effective_address(indices))


class _ArrayGUHelper(namedtuple('_ArrayHelper', ('context', 'builder',
                                                 'shape', 'strides', 'data',
                                                 'layout', 'base_type', 'ndim',
                                                 'inner_arr_ty', 'is_input_arg'))):
    """Helper class to handle array arguments/result.
    It provides methods to generate code loading/storing specific
    items as well as support code for handling indices.

    Contrary to _ArrayHelper, this class can create a view to a subarray
    """
    def create_iter_indices(self):
        intpty = self.context.get_value_type(types.intp)
        ZERO = ir.Constant(ir.IntType(intpty.width), 0)

        indices = []
        for i in range(self.ndim - self.inner_arr_ty.ndim):
            x = cgutils.alloca_once(self.builder, ir.IntType(intpty.width))
            self.builder.store(ZERO, x)
            indices.append(x)
        return _ArrayIndexingHelper(self, indices)

    def _load_effective_address(self, indices):
        context = self.context
        builder = self.builder
        arr_ty = types.Array(self.base_type, self.ndim, self.layout)
        arr = context.make_array(arr_ty)(context, builder, self.data)

        return cgutils.get_item_pointer2(context,
                                         builder,
                                         data=arr.data,
                                         shape=self.shape,
                                         strides=self.strides,
                                         layout=self.layout,
                                         inds=indices)

    def load_data(self, indices):
        context, builder = self.context, self.builder

        if self.inner_arr_ty.ndim == 0 and self.is_input_arg:
            # scalar case for input arguments
            model = context.data_model_manager[self.base_type]
            ptr = self._load_effective_address(indices)
            return model.load_from_data_pointer(builder, ptr)
        elif self.inner_arr_ty.ndim == 0 and not self.is_input_arg:
            # Output arrays are handled as 1d with shape=(1,) when its
            # signature represents a scalar. For instance: "(n),(m) -> ()"
            intpty = context.get_value_type(types.intp)
            one = intpty(1)

            fromty = types.Array(self.base_type, self.ndim, self.layout)
            toty = types.Array(self.base_type, 1, self.layout)
            itemsize = intpty(arrayobj.get_itemsize(context, fromty))

            # create a view from the original ndarray to a 1d array
            arr_from = self.context.make_array(fromty)(context,
                                                       builder,
                                                       self.data)
            arr_to = self.context.make_array(toty)(context, builder)
            arrayobj.populate_array(
                arr_to,
                data=self._load_effective_address(indices),
                shape=cgutils.pack_array(builder, [one]),
                strides=cgutils.pack_array(builder, [itemsize]),
                itemsize=arr_from.itemsize,
                meminfo=arr_from.meminfo,
                parent=arr_from.parent)
            return arr_to._getvalue()
        else:
            # generic case
            # getitem n-dim array -> m-dim array, where N > M
            index_types = (types.int64,) * (self.ndim - self.inner_arr_ty.ndim)
            arrty = types.Array(self.base_type, self.ndim, self.layout)
            arr = self.context.make_array(arrty)(context, builder, self.data)
            res = _getitem_array_generic(context, builder,
                                         self.inner_arr_ty, arrty, arr,
                                         index_types, indices)
            # NOTE: don't call impl_ret_borrowed since the caller doesn't handle
            #       references; but this is a borrow.
            return res

    def guard_shape(self, loopshape):
        inner_ndim = self.inner_arr_ty.ndim
        def raise_impl(loop_shape, array_shape):
            # This would in fact be a test for broadcasting.
            # Broadcast would fail if, ignoring the core dimensions, the
            # remaining ones are different than indices given by loop shape.

            remaining = len(array_shape) - inner_ndim
            _raise = (remaining > len(loop_shape))
            if not _raise:
                for i in range(remaining):
                    _raise |= (array_shape[i] != loop_shape[i])
            if _raise:
                # Ideally we should call `np.broadcast_shapes` with loop and
                # array shapes. But since broadcasting is not supported here,
                # we just raise an error
                # TODO: check why raising a dynamic exception here fails
                raise ValueError('Loop and array shapes are incompatible')

        context, builder = self.context, self.builder
        sig = types.none(
            types.UniTuple(types.intp, len(loopshape)),
            types.UniTuple(types.intp, len(self.shape)),
        )
        tup = (context.make_tuple(builder, sig.args[0], loopshape),
               context.make_tuple(builder, sig.args[1], self.shape))
        context.compile_internal(builder, raise_impl, sig, tup)

    def guard_match_core_dims(self, other: '_ArrayGUHelper', ndims: int):
        # arguments with the same signature should match their core dimensions
        #
        # @guvectorize('(n,m), (n,m) -> (n)')
        # def foo(x, y, res):
        #     ...
        #
        # x and y should have the same core (2D) dimensions
        def raise_impl(self_shape, other_shape):
            same = True
            a, b = len(self_shape) - ndims, len(other_shape) - ndims
            for i in range(ndims):
                same &= self_shape[a + i] == other_shape[b + i]
            if not same:
                # NumPy raises the following:
                # ValueError: gufunc: Input operand 1 has a mismatch in its
                # core dimension 0, with gufunc signature (n),(n) -> ()
                # (size 3 is different from 2)
                # But since we cannot raise a dynamic exception here, we just
                # (try) something meaninful
                msg = ('Operand has a mismatch in one of its core dimensions. '
                       'Please, check if all arguments to a @guvectorize '
                       'function have the same core dimensions.')
                raise ValueError(msg)

        context, builder = self.context, self.builder
        sig = types.none(
            types.UniTuple(types.intp, len(self.shape)),
            types.UniTuple(types.intp, len(other.shape)),
        )
        tup = (context.make_tuple(builder, sig.args[0], self.shape),
               context.make_tuple(builder, sig.args[1], other.shape),)
        context.compile_internal(builder, raise_impl, sig, tup)


def _prepare_argument(ctxt, bld, inp, tyinp, where='input operand'):
    """returns an instance of the appropriate Helper (either
    _ScalarHelper or _ArrayHelper) class to handle the argument.
    using the polymorphic interface of the Helper classes, scalar
    and array cases can be handled with the same code"""

    # first un-Optional Optionals
    if isinstance(tyinp, types.Optional):
        oty = tyinp
        tyinp = tyinp.type
        inp = ctxt.cast(bld, inp, oty, tyinp)

    # then prepare the arg for a concrete instance
    if isinstance(tyinp, types.ArrayCompatible):
        ary     = ctxt.make_array(tyinp)(ctxt, bld, inp)
        shape   = cgutils.unpack_tuple(bld, ary.shape, tyinp.ndim)
        strides = cgutils.unpack_tuple(bld, ary.strides, tyinp.ndim)
        return _ArrayHelper(ctxt, bld, shape, strides, ary.data,
                            tyinp.layout, tyinp.dtype, tyinp.ndim, inp)
    elif (types.unliteral(tyinp) in types.number_domain | {types.boolean}
          or isinstance(tyinp, types.scalars._NPDatetimeBase)):
        return _ScalarHelper(ctxt, bld, inp, tyinp)
    else:
        raise NotImplementedError('unsupported type for {0}: {1}'.format(where,
                                  str(tyinp)))


if config.USE_LEGACY_TYPE_SYSTEM:
    _broadcast_onto_sig = types.intp(types.intp, types.CPointer(types.intp),
                                    types.intp, types.CPointer(types.intp))
else:
    _broadcast_onto_sig = types.np_intp(types.np_intp, types.CPointer(types.np_intp),
                                    types.np_intp, types.CPointer(types.np_intp))

def _broadcast_onto(src_ndim, src_shape, dest_ndim, dest_shape):
    '''Low-level utility function used in calculating a shape for
    an implicit output array.  This function assumes that the
    destination shape is an LLVM pointer to a C-style array that was
    already initialized to a size of one along all axes.

    Returns an integer value:
    >= 1  :  Succeeded.  Return value should equal the number of dimensions in
             the destination shape.
    0     :  Failed to broadcast because source shape is larger than the
             destination shape (this case should be weeded out at type
             checking).
    < 0   :  Failed to broadcast onto destination axis, at axis number ==
             -(return_value + 1).
    '''
    if src_ndim > dest_ndim:
        # This check should have been done during type checking, but
        # let's be defensive anyway...
        return 0
    else:
        src_index = 0
        dest_index = dest_ndim - src_ndim
        while src_index < src_ndim:
            src_dim_size = src_shape[src_index]
            dest_dim_size = dest_shape[dest_index]
            # Check to see if we've already mutated the destination
            # shape along this axis.
            if dest_dim_size != 1:
                # If we have mutated the destination shape already,
                # then the source axis size must either be one,
                # or the destination axis size.
                if src_dim_size != dest_dim_size and src_dim_size != 1:
                    return -(dest_index + 1)
            elif src_dim_size != 1:
                # If the destination size is still its initial
                dest_shape[dest_index] = src_dim_size
            src_index += 1
            dest_index += 1
    return dest_index

def _build_array(context, builder, array_ty, input_types, inputs):
    """Utility function to handle allocation of an implicit output array
    given the target context, builder, output array type, and a list of
    _ArrayHelper instances.
    """
    # First, strip optional types, ufunc loops are typed on concrete types
    input_types = [x.type if isinstance(x, types.Optional) else x
                   for x in input_types]

    intp_ty = context.get_value_type(types.intp)
    def make_intp_const(val):
        return context.get_constant(types.intp, val)

    ZERO = make_intp_const(0)
    ONE = make_intp_const(1)

    src_shape = cgutils.alloca_once(builder, intp_ty, array_ty.ndim,
                                    "src_shape")
    dest_ndim = make_intp_const(array_ty.ndim)
    dest_shape = cgutils.alloca_once(builder, intp_ty, array_ty.ndim,
                                     "dest_shape")
    dest_shape_addrs = tuple(cgutils.gep_inbounds(builder, dest_shape, index)
                             for index in range(array_ty.ndim))

    # Initialize the destination shape with all ones.
    for dest_shape_addr in dest_shape_addrs:
        builder.store(ONE, dest_shape_addr)

    # For each argument, try to broadcast onto the destination shape,
    # mutating along any axis where the argument shape is not one and
    # the destination shape is one.
    for arg_number, arg in enumerate(inputs):
        if not hasattr(arg, "ndim"): # Skip scalar arguments
            continue
        arg_ndim = make_intp_const(arg.ndim)
        for index in range(arg.ndim):
            builder.store(arg.shape[index],
                          cgutils.gep_inbounds(builder, src_shape, index))
        arg_result = context.compile_internal(
            builder, _broadcast_onto, _broadcast_onto_sig,
            [arg_ndim, src_shape, dest_ndim, dest_shape])
        with cgutils.if_unlikely(builder,
                                 builder.icmp_signed('<', arg_result, ONE)):
            msg = "unable to broadcast argument %d to output array" % (
                arg_number,)

            loc = errors.loc_info.get('loc', None)
            if loc is not None:
                msg += '\nFile "%s", line %d, ' % (loc.filename, loc.line)

            context.call_conv.return_user_exc(builder, ValueError, (msg,))

    real_array_ty = array_ty.as_array

    dest_shape_tup = tuple(builder.load(dest_shape_addr)
                           for dest_shape_addr in dest_shape_addrs)
    array_val = arrayobj._empty_nd_impl(context, builder, real_array_ty,
                                        dest_shape_tup)

    # Get the best argument to call __array_wrap__ on
    array_wrapper_index = select_array_wrapper(input_types)
    array_wrapper_ty = input_types[array_wrapper_index]
    try:
        # __array_wrap__(source wrapped array, out array) -> out wrapped array
        array_wrap = context.get_function('__array_wrap__',
                                          array_ty(array_wrapper_ty, real_array_ty))
    except NotImplementedError:
        # If it's the same priority as a regular array, assume we
        # should use the allocated array unchanged.
        if array_wrapper_ty.array_priority != types.Array.array_priority:
            raise
        out_val = array_val._getvalue()
    else:
        wrap_args = (inputs[array_wrapper_index].return_val, array_val._getvalue())
        out_val = array_wrap(builder, wrap_args)

    ndim = array_ty.ndim
    shape   = cgutils.unpack_tuple(builder, array_val.shape, ndim)
    strides = cgutils.unpack_tuple(builder, array_val.strides, ndim)
    return _ArrayHelper(context, builder, shape, strides, array_val.data,
                        array_ty.layout, array_ty.dtype, ndim,
                        out_val)

# ufuncs either return a single result when nout == 1, else a tuple of results

def _unpack_output_types(ufunc, sig):
    if ufunc.nout == 1:
        return [sig.return_type]
    else:
        return list(sig.return_type)


def _unpack_output_values(ufunc, builder, values):
    if ufunc.nout == 1:
        return [values]
    else:
        return cgutils.unpack_tuple(builder, values)


def _pack_output_values(ufunc, context, builder, typ, values):
    if ufunc.nout == 1:
        return values[0]
    else:
        return context.make_tuple(builder, typ, values)


def numpy_ufunc_kernel(context, builder, sig, args, ufunc, kernel_class):
    # This is the code generator that builds all the looping needed
    # to execute a numpy functions over several dimensions (including
    # scalar cases).
    #
    # context - the code generation context
    # builder - the code emitter
    # sig - signature of the ufunc
    # args - the args to the ufunc
    # ufunc - the ufunc itself
    # kernel_class -  a code generating subclass of _Kernel that provides

    arguments = [_prepare_argument(context, builder, arg, tyarg)
                 for arg, tyarg in zip(args, sig.args)]

    if len(arguments) < ufunc.nin:
        raise RuntimeError(
            "Not enough inputs to {}, expected {} got {}"
            .format(ufunc.__name__, ufunc.nin, len(arguments)))

    for out_i, ret_ty in enumerate(_unpack_output_types(ufunc, sig)):
        if ufunc.nin + out_i >= len(arguments):
            # this out argument is not provided
            if isinstance(ret_ty, types.ArrayCompatible):
                output = _build_array(context, builder, ret_ty, sig.args, arguments)
            else:
                output = _prepare_argument(
                    context, builder,
                    ir.Constant(context.get_value_type(ret_ty), None), ret_ty)
            arguments.append(output)
        elif context.enable_nrt:
            # Incref the output
            context.nrt.incref(builder, ret_ty, args[ufunc.nin + out_i])

    inputs = arguments[:ufunc.nin]
    outputs = arguments[ufunc.nin:]
    assert len(outputs) == ufunc.nout

    outer_sig = _ufunc_loop_sig(
        [a.base_type for a in outputs],
        [a.base_type for a in inputs]
    )
    kernel = kernel_class(context, builder, outer_sig)
    intpty = context.get_value_type(types.intp)

    indices = [inp.create_iter_indices() for inp in inputs]

    # assume outputs are all the same size, which numpy requires

    loopshape = outputs[0].shape

    # count the number of C and F layout arrays, respectively
    input_layouts = [inp.layout for inp in inputs
                     if isinstance(inp, _ArrayHelper)]
    num_c_layout = len([x for x in input_layouts if x == 'C'])
    num_f_layout = len([x for x in input_layouts if x == 'F'])

    # Only choose F iteration order if more arrays are in F layout.
    # Default to C order otherwise.
    # This is a best effort for performance. NumPy has more fancy logic that
    # uses array iterators in non-trivial cases.
    if num_f_layout > num_c_layout:
        order = 'F'
    else:
        order = 'C'

    with cgutils.loop_nest(builder, loopshape, intp=intpty, order=order) as loop_indices:
        vals_in = []
        for i, (index, arg) in enumerate(zip(indices, inputs)):
            index.update_indices(loop_indices, i)
            vals_in.append(arg.load_data(index.as_values()))

        vals_out = _unpack_output_values(ufunc, builder, kernel.generate(*vals_in))
        for val_out, output in zip(vals_out, outputs):
            output.store_data(loop_indices, val_out)

    out = _pack_output_values(ufunc, context, builder, sig.return_type, [o.return_val for o in outputs])
    return impl_ret_new_ref(context, builder, sig.return_type, out)


def numpy_gufunc_kernel(context, builder, sig, args, ufunc, kernel_class):
    arguments = []
    expected_ndims = kernel_class.dufunc.expected_ndims()
    expected_ndims = expected_ndims[0] + expected_ndims[1]
    is_input = [True] * ufunc.nin + [False] * ufunc.nout
    for arg, ty, exp_ndim, is_inp in zip(args, sig.args, expected_ndims, is_input):  # noqa: E501
        if isinstance(ty, types.ArrayCompatible):
            # Create an array helper that iteration returns a subarray
            # with ndim specified by "exp_ndim"
            arr = context.make_array(ty)(context, builder, arg)
            shape = cgutils.unpack_tuple(builder, arr.shape, ty.ndim)
            strides = cgutils.unpack_tuple(builder, arr.strides, ty.ndim)
            inner_arr_ty = ty.copy(ndim=exp_ndim)
            ndim = ty.ndim
            layout = ty.layout
            base_type = ty.dtype
            array_helper = _ArrayGUHelper(context, builder,
                                          shape, strides, arg,
                                          layout, base_type, ndim,
                                          inner_arr_ty, is_inp)
            arguments.append(array_helper)
        else:
            scalar_helper = _ScalarHelper(context, builder, arg, ty)
            arguments.append(scalar_helper)
    kernel = kernel_class(context, builder, sig)

    layouts = [arg.layout for arg in arguments
               if isinstance(arg, _ArrayGUHelper)]
    num_c_layout = len([x for x in layouts if x == 'C'])
    num_f_layout = len([x for x in layouts if x == 'F'])

    # Only choose F iteration order if more arrays are in F layout.
    # Default to C order otherwise.
    # This is a best effort for performance. NumPy has more fancy logic that
    # uses array iterators in non-trivial cases.
    if num_f_layout > num_c_layout:
        order = 'F'
    else:
        order = 'C'

    outputs = arguments[ufunc.nin:]
    intpty = context.get_value_type(types.intp)
    indices = [inp.create_iter_indices() for inp in arguments]
    loopshape_ndim = outputs[0].ndim - outputs[0].inner_arr_ty.ndim
    loopshape = outputs[0].shape[ : loopshape_ndim]

    _sig = parse_signature(ufunc.gufunc_builder.signature)
    for (idx_a, sig_a), (idx_b, sig_b) in itertools.combinations(
            zip(range(len(arguments)),
            _sig[0] + _sig[1]),
            r = 2
    ):
        # For each pair of arguments, both inputs and outputs, must match their
        # inner dimensions if their signatures are the same.
        arg_a, arg_b = arguments[idx_a], arguments[idx_b]
        if sig_a == sig_b and \
                all(isinstance(x, _ArrayGUHelper) for x in (arg_a, arg_b)):
            arg_a, arg_b = arguments[idx_a], arguments[idx_b]
            arg_a.guard_match_core_dims(arg_b, len(sig_a))

    for arg in arguments[:ufunc.nin]:
        if isinstance(arg, _ArrayGUHelper):
            arg.guard_shape(loopshape)

    with cgutils.loop_nest(builder,
                           loopshape,
                           intp=intpty,
                           order=order) as loop_indices:
        vals_in = []
        for i, (index, arg) in enumerate(zip(indices, arguments)):
            index.update_indices(loop_indices, i)
            vals_in.append(arg.load_data(index.as_values()))

        kernel.generate(*vals_in)


# Kernels are the code to be executed inside the multidimensional loop.
class _Kernel(object):
    def __init__(self, context, builder, outer_sig):
        self.context = context
        self.builder = builder
        self.outer_sig = outer_sig

    def cast(self, val, fromty, toty):
        """Numpy uses cast semantics that are different from standard Python
        (for example, it does allow casting from complex to float).

        This method acts as a patch to context.cast so that it allows
        complex to real/int casts.

        """
        if (isinstance(fromty, types.Complex) and
            not isinstance(toty, types.Complex)):
            # attempt conversion of the real part to the specified type.
            # note that NumPy issues a warning in this kind of conversions
            newty = fromty.underlying_float
            attr = self.context.get_getattr(fromty, 'real')
            val = attr(self.context, self.builder, fromty, val, 'real')
            fromty = newty
            # let the regular cast do the rest...

        return self.context.cast(self.builder, val, fromty, toty)

    def generate(self, *args):
        isig = self.inner_sig
        osig = self.outer_sig
        cast_args = [self.cast(val, inty, outty)
                     for val, inty, outty in
                     zip(args, osig.args, isig.args)]
        if self.cres.objectmode:
            func_type = self.context.call_conv.get_function_type(
                types.pyobject, [types.pyobject] * len(isig.args))
        else:
            func_type = self.context.call_conv.get_function_type(
                isig.return_type, isig.args)
        module = self.builder.block.function.module
        entry_point = cgutils.get_or_insert_function(
            module, func_type,
            self.cres.fndesc.llvm_func_name)
        entry_point.attributes.add("alwaysinline")

        _, res = self.context.call_conv.call_function(
            self.builder, entry_point, isig.return_type, isig.args,
            cast_args)
        return self.cast(res, isig.return_type, osig.return_type)


def _ufunc_db_function(ufunc):
    """Use the ufunc loop type information to select the code generation
    function from the table provided by the dict_of_kernels. The dict
    of kernels maps the loop identifier to a function with the
    following signature: (context, builder, signature, args).

    The loop type information has the form 'AB->C'. The letters to the
    left of '->' are the input types (specified as NumPy letter
    types).  The letters to the right of '->' are the output
    types. There must be 'ufunc.nin' letters to the left of '->', and
    'ufunc.nout' letters to the right.

    For example, a binary float loop resulting in a float, will have
    the following signature: 'ff->f'.

    A given ufunc implements many loops. The list of loops implemented
    for a given ufunc can be accessed using the 'types' attribute in
    the ufunc object. The NumPy machinery selects the first loop that
    fits a given calling signature (in our case, what we call the
    outer_sig). This logic is mimicked by 'ufunc_find_matching_loop'.
    """

    class _KernelImpl(_Kernel):
        def __init__(self, context, builder, outer_sig):
            super(_KernelImpl, self).__init__(context, builder, outer_sig)
            loop = ufunc_find_matching_loop(
                ufunc, outer_sig.args + tuple(_unpack_output_types(ufunc, outer_sig)))
            self.fn = context.get_ufunc_info(ufunc).get(loop.ufunc_sig)
            self.inner_sig = _ufunc_loop_sig(loop.outputs, loop.inputs)

            if self.fn is None:
                msg = "Don't know how to lower ufunc '{0}' for loop '{1}'"
                raise NotImplementedError(msg.format(ufunc.__name__, loop))

        def generate(self, *args):
            isig = self.inner_sig
            osig = self.outer_sig

            cast_args = [self.cast(val, inty, outty)
                         for val, inty, outty in zip(args, osig.args,
                                                     isig.args)]
            with force_error_model(self.context, 'numpy'):
                res = self.fn(self.context, self.builder, isig, cast_args)
            dmm = self.context.data_model_manager
            res = dmm[isig.return_type].from_return(self.builder, res)
            return self.cast(res, isig.return_type, osig.return_type)

    return _KernelImpl


################################################################################
# Helper functions that register the ufuncs

def register_ufunc_kernel(ufunc, kernel, lower):
    def do_ufunc(context, builder, sig, args):
        return numpy_ufunc_kernel(context, builder, sig, args, ufunc, kernel)

    _any = types.Any
    in_args = (_any,) * ufunc.nin

    # Add a lowering for each out argument that is missing.
    for n_explicit_out in range(ufunc.nout + 1):
        out_args = (types.Array,) * n_explicit_out
        lower(ufunc, *in_args, *out_args)(do_ufunc)

    return kernel


def register_unary_operator_kernel(operator, ufunc, kernel, lower,
                                   inplace=False):
    assert not inplace  # are there any inplace unary operators?
    def lower_unary_operator(context, builder, sig, args):
        return numpy_ufunc_kernel(context, builder, sig, args, ufunc, kernel)
    _arr_kind = types.Array
    lower(operator, _arr_kind)(lower_unary_operator)


def register_binary_operator_kernel(op, ufunc, kernel, lower, inplace=False):
    def lower_binary_operator(context, builder, sig, args):
        return numpy_ufunc_kernel(context, builder, sig, args, ufunc, kernel)

    def lower_inplace_operator(context, builder, sig, args):
        # The visible signature is (A, B) -> A
        # The implementation's signature (with explicit output)
        # is (A, B, A) -> A
        args = tuple(args) + (args[0],)
        sig = typing.signature(sig.return_type, *sig.args + (sig.args[0],))
        return numpy_ufunc_kernel(context, builder, sig, args, ufunc, kernel)

    _any = types.Any
    _arr_kind = types.Array
    formal_sigs = [(_arr_kind, _arr_kind), (_any, _arr_kind), (_arr_kind, _any)]
    for sig in formal_sigs:
        if not inplace:
            lower(op, *sig)(lower_binary_operator)
        else:
            lower(op, *sig)(lower_inplace_operator)


################################################################################
# Use the contents of ufunc_db to initialize the supported ufuncs

@registry.lower(operator.pos, types.Array)
def array_positive_impl(context, builder, sig, args):
    '''Lowering function for +(array) expressions.  Defined here
    (numba.targets.npyimpl) since the remaining array-operator
    lowering functions are also registered in this module.
    '''
    class _UnaryPositiveKernel(_Kernel):
        def generate(self, *args):
            [val] = args
            return val

    return numpy_ufunc_kernel(context, builder, sig, args, np.positive,
                              _UnaryPositiveKernel)


def register_ufuncs(ufuncs, lower):
    kernels = {}
    for ufunc in ufuncs:
        db_func = _ufunc_db_function(ufunc)
        kernels[ufunc] = register_ufunc_kernel(ufunc, db_func, lower)

    for _op_map in (npydecl.NumpyRulesUnaryArrayOperator._op_map,
                    npydecl.NumpyRulesArrayOperator._op_map,
                    ):
        for operator, ufunc_name in _op_map.items():
            ufunc = getattr(np, ufunc_name)
            kernel = kernels[ufunc]
            if ufunc.nin == 1:
                register_unary_operator_kernel(operator, ufunc, kernel, lower)
            elif ufunc.nin == 2:
                register_binary_operator_kernel(operator, ufunc, kernel, lower)
            else:
                raise RuntimeError("There shouldn't be any non-unary or binary operators")

    for _op_map in (npydecl.NumpyRulesInplaceArrayOperator._op_map,
                    ):
        for operator, ufunc_name in _op_map.items():
            ufunc = getattr(np, ufunc_name)
            kernel = kernels[ufunc]
            if ufunc.nin == 1:
                register_unary_operator_kernel(operator, ufunc, kernel, lower,
                                               inplace=True)
            elif ufunc.nin == 2:
                register_binary_operator_kernel(operator, ufunc, kernel, lower,
                                                inplace=True)
            else:
                raise RuntimeError("There shouldn't be any non-unary or binary operators")


register_ufuncs(ufunc_db.get_ufuncs(), registry.lower)


@intrinsic
def _make_dtype_object(typingctx, desc):
    """Given a string or NumberClass description *desc*, returns the dtype object.
    """
    def from_nb_type(nb_type):
        return_type = types.DType(nb_type)
        sig = return_type(desc)

        def codegen(context, builder, signature, args):
            # All dtype objects are dummy values in LLVM.
            # They only exist in the type level.
            return context.get_dummy_value()

        return sig, codegen

    if isinstance(desc, types.Literal):
        # Convert the str description into np.dtype then to numba type.
        nb_type = from_dtype(np.dtype(desc.literal_value))
        return from_nb_type(nb_type)
    elif isinstance(desc, types.functions.NumberClass):
        thestr = str(desc.dtype)
        # Convert the str description into np.dtype then to numba type.
        nb_type = from_dtype(np.dtype(thestr))
        return from_nb_type(nb_type)

@overload(np.dtype)
def numpy_dtype(desc):
    """Provide an implementation so that numpy.dtype function can be lowered.
    """
    if isinstance(desc, (types.Literal, types.functions.NumberClass)):
        def imp(desc):
            return _make_dtype_object(desc)
        return imp
    else:
        raise errors.NumbaTypeError('unknown dtype descriptor: {}'.format(desc))
