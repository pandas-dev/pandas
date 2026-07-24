from collections import namedtuple

import numpy as np

from llvmlite.ir import Constant, IRBuilder
from llvmlite import ir

from numba.core import types, cgutils
from numba.core.compiler_lock import global_compiler_lock
from numba.core.caching import make_library_cache, NullCache


_wrapper_info = namedtuple('_wrapper_info', ['library', 'env', 'name'])


def _build_ufunc_loop_body(load, store, context, func, builder, arrays, out,
                           offsets, store_offset, signature, pyapi, env):
    elems = load()

    # Compute
    status, retval = context.call_conv.call_function(builder, func,
                                                     signature.return_type,
                                                     signature.args, elems)

    # Store
    with builder.if_else(status.is_ok, likely=True) as (if_ok, if_error):
        with if_ok:
            store(retval)
        with if_error:
            gil = pyapi.gil_ensure()
            context.call_conv.raise_error(builder, pyapi, status)
            pyapi.gil_release(gil)

    # increment indices
    for off, ary in zip(offsets, arrays):
        builder.store(builder.add(builder.load(off), ary.step), off)

    builder.store(builder.add(builder.load(store_offset), out.step),
                  store_offset)

    return status.code


def _build_ufunc_loop_body_objmode(load, store, context, func, builder,
                                   arrays, out, offsets, store_offset,
                                   signature, env, pyapi):
    elems = load()

    # Compute
    _objargs = [types.pyobject] * len(signature.args)
    # We need to push the error indicator to avoid it messing with
    # the ufunc's execution.  We restore it unless the ufunc raised
    # a new error.
    with pyapi.err_push(keep_new=True):
        status, retval = context.call_conv.call_function(builder, func,
                                                         types.pyobject,
                                                         _objargs, elems)
        # Release owned reference to arguments
        for elem in elems:
            pyapi.decref(elem)
    # NOTE: if an error occurred, it will be caught by the Numpy machinery

    # Store
    store(retval)

    # increment indices
    for off, ary in zip(offsets, arrays):
        builder.store(builder.add(builder.load(off), ary.step), off)

    builder.store(builder.add(builder.load(store_offset), out.step),
                  store_offset)

    return status.code


def build_slow_loop_body(context, func, builder, arrays, out, offsets,
                         store_offset, signature, pyapi, env):
    def load():
        elems = [ary.load_direct(builder.load(off))
                 for off, ary in zip(offsets, arrays)]
        return elems

    def store(retval):
        out.store_direct(retval, builder.load(store_offset))

    return _build_ufunc_loop_body(load, store, context, func, builder, arrays,
                                  out, offsets, store_offset, signature, pyapi,
                                  env=env)


def build_obj_loop_body(context, func, builder, arrays, out, offsets,
                        store_offset, signature, pyapi, envptr, env):
    env_body = context.get_env_body(builder, envptr)
    env_manager = pyapi.get_env_manager(env, env_body, envptr)

    def load():
        # Load
        elems = [ary.load_direct(builder.load(off))
                 for off, ary in zip(offsets, arrays)]
        # Box
        elems = [pyapi.from_native_value(t, v, env_manager)
                 for v, t in zip(elems, signature.args)]
        return elems

    def store(retval):
        is_ok = cgutils.is_not_null(builder, retval)
        # If an error is raised by the object mode ufunc, it will
        # simply get caught by the Numpy ufunc machinery.
        with builder.if_then(is_ok, likely=True):
            # Unbox
            native = pyapi.to_native_value(signature.return_type, retval)
            assert native.cleanup is None
            # Store
            out.store_direct(native.value, builder.load(store_offset))
            # Release owned reference
            pyapi.decref(retval)

    return _build_ufunc_loop_body_objmode(load, store, context, func, builder,
                                          arrays, out, offsets, store_offset,
                                          signature, envptr, pyapi)


def build_fast_loop_body(context, func, builder, arrays, out, offsets,
                         store_offset, signature, ind, pyapi, env):
    def load():
        elems = [ary.load_aligned(ind)
                 for ary in arrays]
        return elems

    def store(retval):
        out.store_aligned(retval, ind)

    return _build_ufunc_loop_body(load, store, context, func, builder, arrays,
                                  out, offsets, store_offset, signature, pyapi,
                                  env=env)


def build_ufunc_wrapper(library, context, fname, signature, objmode, cres):
    """
    Wrap the scalar function with a loop that iterates over the arguments

    Returns
    -------
    (library, env, name)
    """
    assert isinstance(fname, str)
    byte_t = ir.IntType(8)
    byte_ptr_t = ir.PointerType(byte_t)
    byte_ptr_ptr_t = ir.PointerType(byte_ptr_t)
    intp_t = context.get_value_type(types.intp)
    intp_ptr_t = ir.PointerType(intp_t)

    fnty = ir.FunctionType(ir.VoidType(), [byte_ptr_ptr_t, intp_ptr_t,
                                           intp_ptr_t, byte_ptr_t])

    wrapperlib = context.codegen().create_library('ufunc_wrapper')
    wrapper_module = wrapperlib.create_ir_module('')
    if objmode:
        func_type = context.call_conv.get_function_type(
            types.pyobject, [types.pyobject] * len(signature.args))
    else:
        func_type = context.call_conv.get_function_type(
            signature.return_type, signature.args)

    func = ir.Function(wrapper_module, func_type, name=fname)
    func.attributes.add("alwaysinline")

    wrapper = ir.Function(wrapper_module, fnty, "__ufunc__." + func.name)
    arg_args, arg_dims, arg_steps, arg_data = wrapper.args
    arg_args.name = "args"
    arg_dims.name = "dims"
    arg_steps.name = "steps"
    arg_data.name = "data"

    builder = IRBuilder(wrapper.append_basic_block("entry"))

    # Prepare Environment
    envname = context.get_env_name(cres.fndesc)
    env = cres.environment
    envptr = builder.load(context.declare_env_global(builder.module, envname))

    # Emit loop
    loopcount = builder.load(arg_dims, name="loopcount")

    # Prepare inputs
    arrays = []
    for i, typ in enumerate(signature.args):
        arrays.append(UArrayArg(context, builder, arg_args, arg_steps, i, typ))

    # Prepare output
    out = UArrayArg(context, builder, arg_args, arg_steps, len(arrays),
                    signature.return_type)

    # Setup indices
    offsets = []
    zero = context.get_constant(types.intp, 0)
    for _ in arrays:
        p = cgutils.alloca_once(builder, intp_t)
        offsets.append(p)
        builder.store(zero, p)

    store_offset = cgutils.alloca_once(builder, intp_t)
    builder.store(zero, store_offset)

    unit_strided = cgutils.true_bit
    for ary in arrays:
        unit_strided = builder.and_(unit_strided, ary.is_unit_strided)

    pyapi = context.get_python_api(builder)
    if objmode:
        # General loop
        gil = pyapi.gil_ensure()
        with cgutils.for_range(builder, loopcount, intp=intp_t):
            build_obj_loop_body(
                context, func, builder, arrays, out, offsets,
                store_offset, signature, pyapi, envptr, env,
            )
        pyapi.gil_release(gil)
        builder.ret_void()

    else:
        with builder.if_else(unit_strided) as (is_unit_strided, is_strided):
            with is_unit_strided:
                with cgutils.for_range(builder, loopcount, intp=intp_t) as loop:
                    build_fast_loop_body(
                        context, func, builder, arrays, out, offsets,
                        store_offset, signature, loop.index, pyapi,
                        env=envptr,
                    )

            with is_strided:
                # General loop
                with cgutils.for_range(builder, loopcount, intp=intp_t):
                    build_slow_loop_body(
                        context, func, builder, arrays, out, offsets,
                        store_offset, signature, pyapi,
                        env=envptr,
                    )

        builder.ret_void()
    del builder

    # Link and finalize
    wrapperlib.add_ir_module(wrapper_module)
    wrapperlib.add_linking_library(library)
    return _wrapper_info(library=wrapperlib, env=env, name=wrapper.name)


class UArrayArg(object):
    def __init__(self, context, builder, args, steps, i, fe_type):
        self.context = context
        self.builder = builder
        self.fe_type = fe_type
        offset = self.context.get_constant(types.intp, i)
        offseted_args = self.builder.load(builder.gep(args, [offset]))
        data_type = context.get_data_type(fe_type)
        self.dataptr = self.builder.bitcast(offseted_args,
                                            data_type.as_pointer())
        sizeof = self.context.get_abi_sizeof(data_type)
        self.abisize = self.context.get_constant(types.intp, sizeof)
        offseted_step = self.builder.gep(steps, [offset])
        self.step = self.builder.load(offseted_step)
        self.is_unit_strided = builder.icmp_unsigned('==',
                                                     self.abisize, self.step)
        self.builder = builder

    def load_direct(self, byteoffset):
        """
        Generic load from the given *byteoffset*.  load_aligned() is
        preferred if possible.
        """
        ptr = cgutils.pointer_add(self.builder, self.dataptr, byteoffset)
        return self.context.unpack_value(self.builder, self.fe_type, ptr)

    def load_aligned(self, ind):
        # Using gep() instead of explicit pointer addition helps LLVM
        # vectorize the loop.
        ptr = self.builder.gep(self.dataptr, [ind])
        return self.context.unpack_value(self.builder, self.fe_type, ptr)

    def store_direct(self, value, byteoffset):
        ptr = cgutils.pointer_add(self.builder, self.dataptr, byteoffset)
        self.context.pack_value(self.builder, self.fe_type, value, ptr)

    def store_aligned(self, value, ind):
        ptr = self.builder.gep(self.dataptr, [ind])
        self.context.pack_value(self.builder, self.fe_type, value, ptr)


GufWrapperCache = make_library_cache('guf')


class _GufuncWrapper(object):
    def __init__(self, py_func, cres, sin, sout, cache, is_parfors):
        """
        The *is_parfors* argument is a boolean that indicates if the GUfunc
        being built is to be used as a ParFors kernel. If True, it disables
        the caching on the wrapper as a separate unit because it will be linked
        into the caller function and cached along with it.
        """
        self.py_func = py_func
        self.cres = cres
        self.sin = sin
        self.sout = sout
        self.is_objectmode = self.signature.return_type == types.pyobject
        self.cache = (GufWrapperCache(py_func=self.py_func)
                      if cache else NullCache())
        self.is_parfors = bool(is_parfors)

    @property
    def library(self):
        return self.cres.library

    @property
    def context(self):
        return self.cres.target_context

    @property
    def call_conv(self):
        return self.context.call_conv

    @property
    def signature(self):
        return self.cres.signature

    @property
    def fndesc(self):
        return self.cres.fndesc

    @property
    def env(self):
        return self.cres.environment

    def _wrapper_function_type(self):
        byte_t = ir.IntType(8)
        byte_ptr_t = ir.PointerType(byte_t)
        byte_ptr_ptr_t = ir.PointerType(byte_ptr_t)
        intp_t = self.context.get_value_type(types.intp)
        intp_ptr_t = ir.PointerType(intp_t)

        fnty = ir.FunctionType(ir.VoidType(), [byte_ptr_ptr_t, intp_ptr_t,
                                               intp_ptr_t, byte_ptr_t])
        return fnty

    def _build_wrapper(self, library, name):
        """
        The LLVM IRBuilder code to create the gufunc wrapper.
        The *library* arg is the CodeLibrary to which the wrapper should
        be added.  The *name* arg is the name of the wrapper function being
        created.
        """
        intp_t = self.context.get_value_type(types.intp)
        fnty = self._wrapper_function_type()

        wrapper_module = library.create_ir_module('_gufunc_wrapper')
        func_type = self.call_conv.get_function_type(self.fndesc.restype,
                                                     self.fndesc.argtypes)
        fname = self.fndesc.llvm_func_name
        func = ir.Function(wrapper_module, func_type, name=fname)

        func.attributes.add("alwaysinline")
        wrapper = ir.Function(wrapper_module, fnty, name)
        # The use of weak_odr linkage avoids the function being dropped due
        # to the order in which the wrappers and the user function are linked.
        wrapper.linkage = 'weak_odr'
        arg_args, arg_dims, arg_steps, arg_data = wrapper.args
        arg_args.name = "args"
        arg_dims.name = "dims"
        arg_steps.name = "steps"
        arg_data.name = "data"

        builder = IRBuilder(wrapper.append_basic_block("entry"))
        loopcount = builder.load(arg_dims, name="loopcount")
        pyapi = self.context.get_python_api(builder)

        # Unpack shapes
        unique_syms = set()
        for grp in (self.sin, self.sout):
            for syms in grp:
                unique_syms |= set(syms)

        sym_map = {}
        for syms in self.sin:
            for s in syms:
                if s not in sym_map:
                    sym_map[s] = len(sym_map)

        sym_dim = {}
        for s, i in sym_map.items():
            sym_dim[s] = builder.load(builder.gep(arg_dims,
                                                  [self.context.get_constant(
                                                      types.intp,
                                                      i + 1)]))

        # Prepare inputs
        arrays = []
        step_offset = len(self.sin) + len(self.sout)
        for i, (typ, sym) in enumerate(zip(self.signature.args,
                                           self.sin + self.sout)):
            ary = GUArrayArg(self.context, builder, arg_args,
                             arg_steps, i, step_offset, typ, sym, sym_dim)
            step_offset += len(sym)
            arrays.append(ary)

        bbreturn = builder.append_basic_block('.return')

        # Prologue
        self.gen_prologue(builder, pyapi)

        # Loop
        with cgutils.for_range(builder, loopcount, intp=intp_t) as loop:
            args = [a.get_array_at_offset(loop.index) for a in arrays]
            innercall, error = self.gen_loop_body(builder, pyapi, func, args)
            # If error, escape
            cgutils.cbranch_or_continue(builder, error, bbreturn)

        builder.branch(bbreturn)
        builder.position_at_end(bbreturn)

        # Epilogue
        self.gen_epilogue(builder, pyapi)

        builder.ret_void()

        # Link
        library.add_ir_module(wrapper_module)
        library.add_linking_library(self.library)

    def _compile_wrapper(self, wrapper_name):
        # Gufunc created by Parfors?
        if self.is_parfors:
            # No wrapper caching for parfors
            wrapperlib = self.context.codegen().create_library(str(self))
            # Build wrapper
            self._build_wrapper(wrapperlib, wrapper_name)
        # Non-parfors?
        else:
            # Use cache and compiler in a critical section
            wrapperlib = self.cache.load_overload(
                self.cres.signature, self.cres.target_context,
            )
            if wrapperlib is None:
                # Create library and enable caching
                wrapperlib = self.context.codegen().create_library(str(self))
                wrapperlib.enable_object_caching()
                # Build wrapper
                self._build_wrapper(wrapperlib, wrapper_name)
                # Cache
                self.cache.save_overload(self.cres.signature, wrapperlib)

        return wrapperlib

    @global_compiler_lock
    def build(self):
        wrapper_name = "__gufunc__." + self.fndesc.mangled_name
        wrapperlib = self._compile_wrapper(wrapper_name)
        return _wrapper_info(
            library=wrapperlib, env=self.env, name=wrapper_name,
        )

    def gen_loop_body(self, builder, pyapi, func, args):
        status, retval = self.call_conv.call_function(
            builder, func, self.signature.return_type, self.signature.args,
            args)

        with builder.if_then(status.is_error, likely=False):
            gil = pyapi.gil_ensure()
            self.context.call_conv.raise_error(builder, pyapi, status)
            pyapi.gil_release(gil)

        return status.code, status.is_error

    def gen_prologue(self, builder, pyapi):
        pass        # Do nothing

    def gen_epilogue(self, builder, pyapi):
        pass        # Do nothing


class _GufuncObjectWrapper(_GufuncWrapper):
    def gen_loop_body(self, builder, pyapi, func, args):
        innercall, error = _prepare_call_to_object_mode(self.context,
                                                        builder, pyapi, func,
                                                        self.signature,
                                                        args)
        return innercall, error

    def gen_prologue(self, builder, pyapi):
        # Acquire the GIL
        self.gil = pyapi.gil_ensure()

    def gen_epilogue(self, builder, pyapi):
        # Release GIL
        pyapi.gil_release(self.gil)


def build_gufunc_wrapper(py_func, cres, sin, sout, cache, is_parfors):
    signature = cres.signature
    wrapcls = (_GufuncObjectWrapper
               if signature.return_type == types.pyobject
               else _GufuncWrapper)
    return wrapcls(
        py_func, cres, sin, sout, cache, is_parfors=is_parfors,
    ).build()


def _prepare_call_to_object_mode(context, builder, pyapi, func,
                                 signature, args):
    mod = builder.module

    bb_core_return = builder.append_basic_block('ufunc.core.return')

    # Call to
    # PyObject* ndarray_new(int nd,
    #       npy_intp *dims,   /* shape */
    #       npy_intp *strides,
    #       void* data,
    #       int type_num,
    #       int itemsize)

    ll_int = context.get_value_type(types.int32)
    ll_intp = context.get_value_type(types.intp)
    ll_intp_ptr = ir.PointerType(ll_intp)
    ll_voidptr = context.get_value_type(types.voidptr)
    ll_pyobj = context.get_value_type(types.pyobject)
    fnty = ir.FunctionType(ll_pyobj, [ll_int, ll_intp_ptr,
                                      ll_intp_ptr, ll_voidptr,
                                      ll_int, ll_int])

    fn_array_new = cgutils.get_or_insert_function(mod, fnty,
                                                  "numba_ndarray_new")

    # Convert each llarray into pyobject
    error_pointer = cgutils.alloca_once(builder, ir.IntType(1), name='error')
    builder.store(cgutils.true_bit, error_pointer)

    # The PyObject* arguments to the kernel function
    object_args = []
    object_pointers = []

    for i, (arg, argty) in enumerate(zip(args, signature.args)):
        # Allocate NULL-initialized slot for this argument
        objptr = cgutils.alloca_once(builder, ll_pyobj, zfill=True)
        object_pointers.append(objptr)

        if isinstance(argty, types.Array):
            # Special case arrays: we don't need full-blown NRT reflection
            # since the argument will be gone at the end of the kernel
            arycls = context.make_array(argty)
            array = arycls(context, builder, value=arg)

            zero = Constant(ll_int, 0)

            # Extract members of the llarray
            nd = Constant(ll_int, argty.ndim)
            dims = builder.gep(array._get_ptr_by_name('shape'), [zero, zero])
            strides = builder.gep(array._get_ptr_by_name('strides'),
                                  [zero, zero])
            data = builder.bitcast(array.data, ll_voidptr)
            dtype = np.dtype(str(argty.dtype))

            # Prepare other info for reconstruction of the PyArray
            type_num = Constant(ll_int, dtype.num)
            itemsize = Constant(ll_int, dtype.itemsize)

            # Call helper to reconstruct PyArray objects
            obj = builder.call(fn_array_new, [nd, dims, strides, data,
                                              type_num, itemsize])
        else:
            # Other argument types => use generic boxing
            obj = pyapi.from_native_value(argty, arg)

        builder.store(obj, objptr)
        object_args.append(obj)

        obj_is_null = cgutils.is_null(builder, obj)
        builder.store(obj_is_null, error_pointer)
        cgutils.cbranch_or_continue(builder, obj_is_null, bb_core_return)

    # Call ufunc core function
    object_sig = [types.pyobject] * len(object_args)

    status, retval = context.call_conv.call_function(
        builder, func, types.pyobject, object_sig,
        object_args)
    builder.store(status.is_error, error_pointer)

    # Release returned object
    pyapi.decref(retval)

    builder.branch(bb_core_return)
    # At return block
    builder.position_at_end(bb_core_return)

    # Release argument objects
    for objptr in object_pointers:
        pyapi.decref(builder.load(objptr))

    innercall = status.code
    return innercall, builder.load(error_pointer)


class GUArrayArg(object):
    def __init__(self, context, builder, args, steps, i, step_offset,
                 typ, syms, sym_dim):

        self.context = context
        self.builder = builder

        offset = context.get_constant(types.intp, i)

        data = builder.load(builder.gep(args, [offset], name="data.ptr"),
                            name="data")
        self.data = data

        core_step_ptr = builder.gep(steps, [offset], name="core.step.ptr")
        core_step = builder.load(core_step_ptr)

        if isinstance(typ, types.Array):
            as_scalar = not syms

            # number of symbol in the shape spec should match the dimension
            # of the array type.
            if len(syms) != typ.ndim:
                if len(syms) == 0 and typ.ndim == 1:
                    # This is an exception for handling scalar argument.
                    # The type can be 1D array for scalar.
                    # In the future, we may deprecate this exception.
                    pass
                else:
                    raise TypeError("type and shape signature mismatch for arg "
                                    "#{0}".format(i + 1))

            ndim = typ.ndim
            shape = [sym_dim[s] for s in syms]
            strides = []

            for j in range(ndim):
                stepptr = builder.gep(steps,
                                      [context.get_constant(types.intp,
                                                            step_offset + j)],
                                      name="step.ptr")
                step = builder.load(stepptr)
                strides.append(step)

            ldcls = (_ArrayAsScalarArgLoader
                     if as_scalar
                     else _ArrayArgLoader)

            self._loader = ldcls(dtype=typ.dtype,
                                 ndim=ndim,
                                 core_step=core_step,
                                 as_scalar=as_scalar,
                                 shape=shape,
                                 strides=strides)
        else:
            # If typ is not an array
            if syms:
                raise TypeError("scalar type {0} given for non scalar "
                                "argument #{1}".format(typ, i + 1))
            self._loader = _ScalarArgLoader(dtype=typ, stride=core_step)

    def get_array_at_offset(self, ind):
        return self._loader.load(context=self.context, builder=self.builder,
                                 data=self.data, ind=ind)


class _ScalarArgLoader(object):
    """
    Handle GFunc argument loading where a scalar type is used in the core
    function.
    Note: It still has a stride because the input to the gufunc can be an array
          for this argument.
    """

    def __init__(self, dtype, stride):
        self.dtype = dtype
        self.stride = stride

    def load(self, context, builder, data, ind):
        # Load at base + ind * stride
        data = builder.gep(data, [builder.mul(ind, self.stride)])
        dptr = builder.bitcast(data,
                               context.get_data_type(self.dtype).as_pointer())
        return builder.load(dptr)


class _ArrayArgLoader(object):
    """
    Handle GUFunc argument loading where an array is expected.
    """

    def __init__(self, dtype, ndim, core_step, as_scalar, shape, strides):
        self.dtype = dtype
        self.ndim = ndim
        self.core_step = core_step
        self.as_scalar = as_scalar
        self.shape = shape
        self.strides = strides

    def load(self, context, builder, data, ind):
        arytyp = types.Array(dtype=self.dtype, ndim=self.ndim, layout="A")
        arycls = context.make_array(arytyp)

        array = arycls(context, builder)
        offseted_data = cgutils.pointer_add(builder,
                                            data,
                                            builder.mul(self.core_step,
                                                        ind))

        shape, strides = self._shape_and_strides(context, builder)

        itemsize = context.get_abi_sizeof(context.get_data_type(self.dtype))
        context.populate_array(array,
                               data=builder.bitcast(offseted_data,
                                                    array.data.type),
                               shape=shape,
                               strides=strides,
                               itemsize=context.get_constant(types.intp,
                                                             itemsize),
                               meminfo=None)

        return array._getvalue()

    def _shape_and_strides(self, context, builder):
        shape = cgutils.pack_array(builder, self.shape)
        strides = cgutils.pack_array(builder, self.strides)
        return shape, strides


class _ArrayAsScalarArgLoader(_ArrayArgLoader):
    """
    Handle GUFunc argument loading where the shape signature specifies
    a scalar "()" but a 1D array is used for the type of the core function.
    """

    def _shape_and_strides(self, context, builder):
        # Set shape and strides for a 1D size 1 array
        one = context.get_constant(types.intp, 1)
        zero = context.get_constant(types.intp, 0)
        shape = cgutils.pack_array(builder, [one])
        strides = cgutils.pack_array(builder, [zero])
        return shape, strides
