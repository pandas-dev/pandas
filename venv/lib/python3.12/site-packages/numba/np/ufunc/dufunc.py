import functools
import operator
import warnings

import numpy as np

from numba import jit, typeof
from numba.core import cgutils, types, serialize, sigutils, errors
from numba.core.extending import (is_jitted, overload_attribute,
                                  overload_method, register_jitable,
                                  intrinsic)
from numba.core.typing import npydecl
from numba.core.typing.templates import AbstractTemplate, signature
from numba.cpython.unsafe.tuple import tuple_setitem
from numba.np.ufunc import _internal
from numba.np.ufunc.ufunc_base import UfuncBase, UfuncLowererBase
from numba.parfors import array_analysis
from numba.np.ufunc import ufuncbuilder
from numba.np import numpy_support
from typing import Callable
from llvmlite import ir
from numba.core.compiler_lock import global_compiler_lock


class UfuncAtIterator:

    def __init__(self, ufunc, a, a_ty, indices, indices_ty, b=None, b_ty=None):
        self.ufunc = ufunc
        self.a = a
        self.a_ty = a_ty
        self.indices = indices
        self.indices_ty = indices_ty
        self.b = b
        self.b_ty = b_ty

    def run(self, context, builder):
        self._prepare(context, builder)
        loop_indices, _ = self.indexer.begin_loops()
        self._call_ufunc(context, builder, loop_indices)
        self.indexer.end_loops()

    def need_advanced_indexing(self):
        return isinstance(self.indices_ty, types.BaseTuple)

    def _prepare(self, context, builder):
        from numba.np.arrayobj import normalize_indices, FancyIndexer

        a, indices = self.a, self.indices
        a_ty, indices_ty = self.a_ty, self.indices_ty

        zero = context.get_value_type(types.intp)(0)

        if self.b is not None:
            self.b_indice = cgutils.alloca_once_value(builder, zero)

        if self.need_advanced_indexing():
            indices = cgutils.unpack_tuple(builder, indices,
                                           count=len(indices_ty))
            index_types = indices_ty.types
            index_types, indices = normalize_indices(context, builder,
                                                     index_types, indices)
        else:
            indices = (indices,)
            index_types = (indices_ty,)
            index_types, indices = normalize_indices(context, builder,
                                                     index_types, indices)

        self.indexer = FancyIndexer(context, builder, a_ty, a,
                                    index_types, indices)
        self.indexer.prepare()
        self.cres = self._compile_ufunc(context, builder)

    def _load_val(self, context, builder, loop_indices, array, array_ty):
        from numba.np.arrayobj import load_item
        shapes = cgutils.unpack_tuple(builder, array.shape)
        strides = cgutils.unpack_tuple(builder, array.strides)
        data = array.data

        ptr = cgutils.get_item_pointer2(context, builder, data, shapes, strides,
                                        array_ty.layout, loop_indices)
        val = load_item(context, builder, array_ty, ptr)
        return ptr, val

    def _load_flat(self, context, builder, indices, array, array_ty):
        idx = builder.load(indices)
        sig = array_ty.dtype(array_ty, types.intp)
        impl = context.get_function(operator.getitem, sig)
        val = impl(builder, (array, idx))

        # increment indices
        one = context.get_value_type(types.intp)(1)
        idx = builder.add(idx, one)
        builder.store(idx, indices)

        return None, val

    def _store_val(self, context, builder, array, array_ty, ptr, val):
        from numba.np.arrayobj import store_item
        fromty = self.cres.signature.return_type
        toty = array_ty.dtype
        val = context.cast(builder, val, fromty, toty)
        store_item(context, builder, array_ty, val, ptr)

    def _compile_ufunc(self, context, builder):
        ufunc = self.ufunc.key[0]

        if self.b is None:
            sig = (self.a_ty.dtype,)
        else:
            sig = (self.a_ty.dtype, self.b_ty.dtype)

        cres = ufunc.add(sig)
        context.add_linking_libs((cres.library,))
        return cres

    def _call_ufunc(self, context, builder, loop_indices):
        cres = self.cres
        a, a_ty = self.a, self.a_ty

        ptr, val = self._load_val(context, builder, loop_indices, a, a_ty)

        if self.b is None:
            args = (val,)
        else:
            b, b_ty, b_idx = self.b, self.b_ty, self.b_indice
            _, val_b = self._load_flat(context, builder, b_idx, b, b_ty)
            args = (val, val_b)

        res = context.call_internal(builder, cres.fndesc, cres.signature,
                                    args)
        self._store_val(context, builder, a, a_ty, ptr, res)


def make_dufunc_kernel(_dufunc):
    from numba.np import npyimpl

    class DUFuncKernel(npyimpl._Kernel):
        """
        npyimpl._Kernel subclass responsible for lowering a DUFunc kernel
        (element-wise function) inside a broadcast loop (which is
        generated by npyimpl.numpy_ufunc_kernel()).
        """
        dufunc = _dufunc

        def __init__(self, context, builder, outer_sig):
            super().__init__(context, builder, outer_sig)
            self.inner_sig, self.cres = self.dufunc.find_ewise_function(
                outer_sig.args)

    DUFuncKernel.__name__ += _dufunc.ufunc.__name__
    return DUFuncKernel


class DUFuncLowerer(UfuncLowererBase):
    '''Callable class responsible for lowering calls to a specific DUFunc.
    '''
    def __init__(self, dufunc):
        from numba.np import npyimpl
        super().__init__(dufunc,
                         make_dufunc_kernel,
                         npyimpl.numpy_ufunc_kernel)


class DUFunc(serialize.ReduceMixin, _internal._DUFunc, UfuncBase):
    """
    Dynamic universal function (DUFunc) intended to act like a normal
    Numpy ufunc, but capable of call-time (just-in-time) compilation
    of fast loops specialized to inputs.
    """
    # NOTE: __base_kwargs must be kept in synch with the kwlist in
    # _internal.c:dufunc_init()
    __base_kwargs = set(('identity', '_keepalive', 'nin', 'nout'))

    def __init__(self, py_func, identity=None, cache=False, targetoptions=None):
        if targetoptions is None:
            targetoptions = {}
        if is_jitted(py_func):
            py_func = py_func.py_func
        with ufuncbuilder._suppress_deprecation_warning_nopython_not_supplied():
            dispatcher = jit(_target='npyufunc',
                             cache=cache,
                             **targetoptions)(py_func)
        self._initialize(dispatcher, identity)
        functools.update_wrapper(self, py_func)

    def _initialize(self, dispatcher, identity):
        identity = ufuncbuilder.parse_identity(identity)
        super(DUFunc, self).__init__(dispatcher, identity=identity)
        # Loop over a copy of the keys instead of the keys themselves,
        # since we're changing the dictionary while looping.
        self.reorderable = (identity != _internal.PyUFunc_None)
        self.__name__ = dispatcher.py_func.__name__
        self.__doc__ = dispatcher.py_func.__doc__
        self._lower_me = DUFuncLowerer(self)
        self._install_cg()
        self._install_type()

    def _reduce_states(self):
        """
        NOTE: part of ReduceMixin protocol
        """
        siglist = list(self._dispatcher.overloads.keys())
        return dict(
            dispatcher=self._dispatcher,
            identity=self.identity,
            frozen=self._frozen,
            siglist=siglist,
        )

    @classmethod
    def _rebuild(cls, dispatcher, identity, frozen, siglist):
        """
        NOTE: part of ReduceMixin protocol
        """
        self = _internal._DUFunc.__new__(cls)
        self._initialize(dispatcher, identity)
        # Re-add signatures
        for sig in siglist:
            self.add(sig)
        if frozen:
            self.disable_compile()
        return self

    def build_ufunc(self):
        """
        For compatibility with the various *UFuncBuilder classes.
        """
        return self

    @property
    def targetoptions(self):
        return self._dispatcher.targetoptions

    @property
    def nin(self):
        return self.ufunc.nin

    @property
    def nout(self):
        return self.ufunc.nout

    @property
    def nargs(self):
        return self.ufunc.nargs

    @property
    def ntypes(self):
        return self.ufunc.ntypes

    @property
    def types(self):
        return self.ufunc.types

    @property
    def identity(self):
        return self.ufunc.identity

    @property
    def signature(self):
        return self.ufunc.signature

    def disable_compile(self):
        """
        Disable the compilation of new signatures at call time.
        """
        # If disabling compilation then there must be at least one signature
        assert len(self._dispatcher.overloads) > 0
        self._frozen = True

    def add(self, sig):
        """
        Compile the DUFunc for the given signature.
        """
        args, return_type = sigutils.normalize_signature(sig)
        return self._compile_for_argtys(args, return_type)

    def __call__(self, *args, **kws):
        """
        Allow any argument that has overridden __array_ufunc__ (NEP-18)
        to take control of DUFunc.__call__.
        """
        default = numpy_support.np.ndarray.__array_ufunc__

        for arg in args + tuple(kws.values()):
            if getattr(type(arg), "__array_ufunc__", default) is not default:
                output = arg.__array_ufunc__(self, "__call__", *args, **kws)
                if output is not NotImplemented:
                    return output
        else:
            return super().__call__(*args, **kws)

    def _compile_for_args(self, *args, **kws):
        nin = self.ufunc.nin
        if kws:
            if 'out' in kws:
                out = kws.pop('out')
                args += (out,)
            if kws:
                raise TypeError("unexpected keyword arguments to ufunc: %s"
                                % ", ".join(repr(k) for k in sorted(kws)))

        args_len = len(args)
        assert (args_len == nin) or (args_len == nin + self.ufunc.nout)
        assert not kws
        argtys = []
        for arg in args[:nin]:
            argty = typeof(arg)
            if isinstance(argty, types.Array):
                argty = argty.dtype
            else:
                # To avoid a mismatch in how Numba types scalar values as
                # opposed to Numpy, we need special logic for scalars.
                # For example, on 64-bit systems, numba.typeof(3) => int32, but
                # np.array(3).dtype => int64.

                # Note: this will not handle numpy "duckarrays" correctly,
                # including but not limited to those implementing `__array__`
                # and `__array_ufunc__`.
                argty = numpy_support.map_arrayscalar_type(arg)
            argtys.append(argty)
        return self._compile_for_argtys(tuple(argtys))

    @global_compiler_lock
    def _compile_for_argtys(self, argtys, return_type=None):
        """
        Given a tuple of argument types (these should be the array
        dtypes, and not the array types themselves), compile the
        element-wise function for those inputs, generate a UFunc loop
        wrapper, and register the loop with the Numpy ufunc object for
        this DUFunc.
        """
        if self._frozen:
            raise RuntimeError("compilation disabled for %s" % (self,))
        assert isinstance(argtys, tuple)
        if return_type is None:
            sig = argtys
        else:
            sig = return_type(*argtys)

        for k, cres in self._dispatcher.overloads.items():
            if argtys == k.args:
                msg = ("Compilation requested for previously compiled argument"
                       f" types ({argtys}). This has no effect and perhaps "
                       "indicates a bug in the calling code (compiling a "
                       "ufunc more than once for the same signature")
                warnings.warn(msg, errors.NumbaWarning)
                return cres

        cres, argtys, return_type = ufuncbuilder._compile_element_wise_function(
            self._dispatcher, self.targetoptions, sig)
        actual_sig = ufuncbuilder._finalize_ufunc_signature(
            cres, argtys, return_type)
        dtypenums, ptr, env = ufuncbuilder._build_element_wise_ufunc_wrapper(
            cres, actual_sig)
        self._add_loop(int(ptr), dtypenums)
        self._keepalive.append((ptr, cres.library, env))
        self._lower_me.libs.append(cres.library)
        return cres

    def match_signature(self, ewise_types, sig):
        return sig.args == ewise_types

    def _install_ufunc_attributes(self, template) -> None:

        def get_attr_fn(attr: str) -> Callable:

            def impl(ufunc):
                val = getattr(ufunc.key[0], attr)
                return lambda ufunc: val
            return impl

        # ntypes/types needs "at" to be a BoundFunction rather than a Function
        # But this fails as it cannot a weak reference to an ufunc due to NumPy
        # not setting the "tp_weaklistoffset" field. See:
        # https://github.com/numpy/numpy/blob/7fc72776b972bfbfdb909e4b15feb0308cf8adba/numpy/core/src/umath/ufunc_object.c#L6968-L6983  # noqa: E501

        at = types.Function(template)
        attributes = ('nin', 'nout', 'nargs', # 'ntypes', # 'types',
                      'identity', 'signature')
        for attr in attributes:
            attr_fn = get_attr_fn(attr)
            overload_attribute(at, attr)(attr_fn)

    def _install_ufunc_methods(self, template) -> None:
        self._install_ufunc_reduce(template)
        self._install_ufunc_reduceat(template)
        self._install_ufunc_at(template)

    def _install_ufunc_at(self, template) -> None:
        at = types.Function(template)

        @overload_method(at, 'at')
        def ol_at(ufunc, a, indices, b=None):
            warnings.warn("ufunc.at feature is experimental",
                          category=errors.NumbaExperimentalFeatureWarning)

            if not isinstance(a, types.Array):
                msg = 'The first argument "a" must be array-like'
                raise errors.NumbaTypeError(msg)

            indices_arr = isinstance(indices, types.Array)
            indices_list = isinstance(indices, types.List)
            indices_tuple = isinstance(indices, types.Tuple)
            indices_slice = isinstance(indices, types.SliceType)
            indices_scalar = not (indices_arr or indices_slice or indices_tuple)
            indices_empty_tuple = indices_tuple and len(indices) == 0
            b_array = isinstance(b, (types.Array, types.Sequence, types.List,
                                     types.Tuple))
            b_none = cgutils.is_nonelike(b)
            b_scalar = not (b_array or b_none)
            need_cast = any([indices_list])

            nin = self.ufunc.nin

            # missing second argument?
            if nin == 2 and cgutils.is_nonelike(b):
                raise errors.TypingError('second operand needed for ufunc')

            # extra second argument
            if nin == 1 and not cgutils.is_nonelike(b):
                msg = 'second operand provided when ufunc is unary'
                raise errors.TypingError(msg)

            if cgutils.is_nonelike(b):
                self.add((a.dtype,))
            elif b_scalar:
                self.add((a.dtype, b))
            else:
                self.add((a.dtype, b.dtype))

            def apply_ufunc_codegen(context, builder, sig, args):
                from numba.np.arrayobj import make_array

                if len(args) == 4:
                    _, aty, idxty, bty = sig.args
                    _, a, indices, b = args
                else:
                    _, aty, idxty, bty = sig.args + (None,)
                    _, a, indices, b = args + (None,)

                a = make_array(aty)(context, builder, a)
                at_iter = UfuncAtIterator(ufunc, a, aty, indices, idxty, b, bty)
                at_iter.run(context, builder)

            @intrinsic
            def apply_a_b_ufunc(typingctx, ufunc, a, indices, b):
                sig = types.none(ufunc, a, indices, b)
                return sig, apply_ufunc_codegen

            @intrinsic
            def apply_a_ufunc(typingctx, ufunc, a, indices):
                sig = types.none(ufunc, a, indices)
                return sig, apply_ufunc_codegen

            def impl_cast(ufunc, a, indices, b=None):
                if b_none:
                    return ufunc.at(a, np.asarray(indices))
                else:
                    return ufunc.at(a,
                                    np.asarray(indices),
                                    np.asarray(b))

            def impl_generic(ufunc, a, indices, b=None):
                if b_none:
                    apply_a_ufunc(ufunc, a, indices,)
                else:
                    b_ = np.asarray(b)
                    a_ = a[indices]
                    b_ = np.broadcast_to(b_, a_.shape)
                    apply_a_b_ufunc(ufunc, a, indices, b_.flat)

            def impl_indices_empty_b_scalar(ufunc, a, indices, b=None):
                a[()] = ufunc(a[()], b)

            def impl_scalar_scalar(ufunc, a, indices, b=None):
                if b_none:
                    a[indices] = ufunc(a[indices])
                else:
                    a[indices] = ufunc(a[indices], b)

            if need_cast:
                return impl_cast
            elif indices_empty_tuple and b_scalar:
                return impl_indices_empty_b_scalar
            elif indices_scalar and b_scalar:
                return impl_scalar_scalar
            else:
                return impl_generic

    def _install_ufunc_reduce(self, template) -> None:
        at = types.Function(template)

        @overload_method(at, 'reduce')
        def ol_reduce(ufunc, array, axis=0, dtype=None, initial=None):

            warnings.warn("ufunc.reduce feature is experimental",
                          category=errors.NumbaExperimentalFeatureWarning)

            if not isinstance(array, types.Array):
                msg = 'The first argument "array" must be array-like'
                raise errors.NumbaTypeError(msg)

            axis_int_tuple = isinstance(axis, types.UniTuple) and \
                isinstance(axis.dtype, types.Integer)
            axis_empty_tuple = isinstance(axis, types.Tuple) and len(axis) == 0
            axis_none = cgutils.is_nonelike(axis)

            identity_none = self.ufunc.identity is None
            ufunc_name = self.ufunc.__name__

            # In NumPy, a ufunc is reorderable if its identity type is **not**
            # PyUfunc_None.
            if not self.reorderable and axis_int_tuple and len(axis) > 1:
                msg = (f"reduction operation '{ufunc_name}' is not "
                       "reorderable, so at most one axis may be specified")
                raise errors.NumbaTypeError(msg)

            tup_init = (0,) * (array.ndim)
            tup_init_m1 = (0,) * (array.ndim - 1)
            nb_dtype = array.dtype if cgutils.is_nonelike(dtype) else dtype
            identity = self.identity

            id_none = cgutils.is_nonelike(identity)
            init_none = cgutils.is_nonelike(initial)

            @register_jitable
            def tuple_slice(tup, pos):
                # Same as
                # tup = tup[0 : pos] + tup[pos + 1:]
                s = tup_init_m1
                i = 0
                for j, e in enumerate(tup):
                    if j == pos:
                        continue
                    s = tuple_setitem(s, i, e)
                    i += 1
                return s

            @register_jitable
            def tuple_slice_append(tup, pos, val):
                # Same as
                # tup = tup[0 : pos] + val + tup[pos + 1:]
                s = tup_init
                i, j, sz = 0, 0, len(s)
                while j < sz:
                    if j == pos:
                        s = tuple_setitem(s, j, val)
                    else:
                        e = tup[i]
                        s = tuple_setitem(s, j, e)
                        i += 1
                    j += 1
                return s

            @intrinsic
            def compute_flat_idx(typingctx, strides, itemsize, idx, axis):
                sig = types.intp(strides, itemsize, idx, axis)
                len_idx = len(idx)

                def gen_block(builder, block_pos, block_name, bb_end, args):
                    strides, _, idx, _ = args
                    bb = builder.append_basic_block(name=block_name)

                    with builder.goto_block(bb):
                        zero = ir.IntType(64)(0)
                        flat_idx = zero

                        if block_pos == 0:
                            for i in range(1, len_idx):
                                stride = builder.extract_value(strides, i - 1)
                                idx_i = builder.extract_value(idx, i)
                                m = builder.mul(stride, idx_i)
                                flat_idx = builder.add(flat_idx, m)
                        elif 0 < block_pos < len_idx - 1:
                            for i in range(0, block_pos):
                                stride = builder.extract_value(strides, i)
                                idx_i = builder.extract_value(idx, i)
                                m = builder.mul(stride, idx_i)
                                flat_idx = builder.add(flat_idx, m)

                            for i in range(block_pos + 1, len_idx):
                                stride = builder.extract_value(strides, i - 1)
                                idx_i = builder.extract_value(idx, i)
                                m = builder.mul(stride, idx_i)
                                flat_idx = builder.add(flat_idx, m)
                        else:
                            for i in range(0, len_idx - 1):
                                stride = builder.extract_value(strides, i)
                                idx_i = builder.extract_value(idx, i)
                                m = builder.mul(stride, idx_i)
                                flat_idx = builder.add(flat_idx, m)

                        builder.branch(bb_end)

                    return bb, flat_idx

                def codegen(context, builder, sig, args):
                    strides, itemsize, idx, axis = args

                    bb = builder.basic_block
                    switch_end = builder.append_basic_block(name='axis_end')
                    l = []
                    for i in range(len_idx):
                        block, flat_idx = gen_block(builder, i, f"axis_{i}",
                                                    switch_end, args)
                        l.append((block, flat_idx))

                    with builder.goto_block(bb):
                        switch = builder.switch(axis, l[-1][0])
                        for i in range(len_idx):
                            switch.add_case(i, l[i][0])

                    builder.position_at_end(switch_end)
                    phi = builder.phi(l[0][1].type)
                    for block, value in l:
                        phi.add_incoming(value, block)
                    return builder.sdiv(phi, itemsize)

                return sig, codegen

            @register_jitable
            def fixup_axis(axis, ndim):
                ax = axis
                for i in range(len(axis)):
                    val = axis[i] + ndim if axis[i] < 0 else axis[i]
                    ax = tuple_setitem(ax, i, val)
                return ax

            @register_jitable
            def find_min(tup):
                idx, e = 0, tup[0]
                for i in range(len(tup)):
                    if tup[i] < e:
                        idx, e = i, tup[i]
                return idx, e

            def impl_1d(ufunc, array, axis=0, dtype=None, initial=None):
                if identity_none and initial is None and len(array) == 0:
                    msg = ('zero-size array to reduction operation '
                           f'{ufunc_name} which has no identity')
                    raise ValueError(msg)

                start = 0
                if init_none and id_none:
                    start = 1
                    r = array[0]
                elif init_none:
                    r = identity
                else:
                    r = initial

                sz = array.shape[0]
                for i in range(start, sz):
                    r = ufunc(r, array[i])
                return r

            def impl_nd_axis_int(ufunc,
                                 array,
                                 axis=0,
                                 dtype=None,
                                 initial=None):
                if axis is None:
                    raise ValueError("'axis' must be specified")

                if axis < 0:
                    axis += array.ndim

                if axis < 0 or axis >= array.ndim:
                    raise ValueError("Invalid axis")

                if identity_none and initial is None and array.shape[axis] == 0:
                    msg = ('zero-size array to reduction operation '
                           f'{ufunc_name} which has no identity')
                    raise ValueError(msg)

                # create result array
                shape = tuple_slice(array.shape, axis)

                if initial is None and identity is None:
                    r = np.empty(shape, dtype=nb_dtype)
                    for idx, _ in np.ndenumerate(r):
                        # shape[0:axis] + 0 + shape[axis:]
                        result_idx = tuple_slice_append(idx, axis, 0)
                        r[idx] = array[result_idx]
                elif initial is None and identity is not None:
                    # Checking if identity is not none is redundant but required
                    # compile this block
                    r = np.full(shape, fill_value=identity, dtype=nb_dtype)
                else:
                    r = np.full(shape, fill_value=initial, dtype=nb_dtype)

                # One approach to implement reduce is to remove the axis index
                # from the indexing tuple returned by "np.ndenumerate". For
                # instance, if idx = (X, Y, Z) and axis=1, the result index
                # is (X, Y).
                # Another way is to compute the result index using strides,
                # which is faster than manipulating tuples.
                view = r.ravel()
                if initial is None and identity is None:
                    for idx, val in np.ndenumerate(array):
                        if idx[axis] == 0:
                            continue
                        else:
                            flat_pos = compute_flat_idx(r.strides, r.itemsize,
                                                        idx, axis)
                            lhs, rhs = view[flat_pos], val
                            view[flat_pos] = ufunc(lhs, rhs)
                else:
                    for idx, val in np.ndenumerate(array):
                        if initial is None and identity is None and \
                                idx[axis] == 0:
                            continue
                        flat_pos = compute_flat_idx(r.strides, r.itemsize,
                                                    idx, axis)
                        lhs, rhs = view[flat_pos], val
                        view[flat_pos] = ufunc(lhs, rhs)
                return r

            def impl_nd_axis_tuple(ufunc,
                                   array,
                                   axis=0,
                                   dtype=None,
                                   initial=None):
                axis_ = fixup_axis(axis, array.ndim)
                for i in range(0, len(axis_)):
                    if axis_[i] < 0 or axis_[i] >= array.ndim:
                        raise ValueError("Invalid axis")

                    for j in range(i + 1, len(axis_)):
                        if axis_[i] == axis_[j]:
                            raise ValueError("duplicate value in 'axis'")

                min_idx, min_elem = find_min(axis_)
                r = ufunc.reduce(array,
                                 axis=min_elem,
                                 dtype=dtype,
                                 initial=initial)
                if len(axis) == 1:
                    return r
                elif len(axis) == 2:
                    return ufunc.reduce(r, axis=axis_[(min_idx + 1) % 2] - 1)
                else:
                    ax = axis_tup
                    for i in range(len(ax)):
                        if i != min_idx:
                            ax = tuple_setitem(ax, i, axis_[i])
                    return ufunc.reduce(r, axis=ax)

            def impl_axis_empty_tuple(ufunc,
                                      array,
                                      axis=0,
                                      dtype=None,
                                      initial=None):
                return array

            def impl_axis_none(ufunc,
                               array,
                               axis=0,
                               dtype=None,
                               initial=None):
                return ufunc.reduce(array, axis_tup, dtype, initial)

            if array.ndim == 1 and not axis_empty_tuple:
                return impl_1d
            elif axis_empty_tuple:
                # ufunc(array, axis=())
                return impl_axis_empty_tuple
            elif axis_none:
                # ufunc(array, axis=None)
                axis_tup = tuple(range(array.ndim))
                return impl_axis_none
            elif axis_int_tuple:
                # axis is tuple of integers
                # ufunc(array, axis=(1, 2, ...))
                axis_tup = (0,) * (len(axis) - 1)
                return impl_nd_axis_tuple
            elif axis == 0 or isinstance(axis, (types.Integer,
                                                types.Omitted,
                                                types.IntegerLiteral)):
                # axis is default value (0) or an integer
                # ufunc(array, axis=0)
                return impl_nd_axis_int

    def _install_ufunc_reduceat(self, template) -> None:
        at = types.Function(template)

        @overload_method(at, 'reduceat')
        def ol_reduceat(ufunc, array, indices, axis=0, dtype=None, out=None):

            warnings.warn("ufunc.reduceat feature is experimental",
                          category=errors.NumbaExperimentalFeatureWarning)

            if self.ufunc.nin != 2:
                msg = 'reduceat only supported for binary functions'
                raise errors.NumbaTypeError(msg)

            if not numpy_support.type_can_asarray(array):
                msg = 'The first argument "array" must be array-like'
                raise errors.NumbaTypeError(msg)

            if not numpy_support.type_can_asarray(indices):
                msg = 'The second argument "indices" must be array-like'
                raise errors.NumbaTypeError(msg)

            if not isinstance(axis, (types.Integer, types.Omitted,
                                     types.IntegerLiteral, int)):
                msg = '"axis" must be an integer'
                raise errors.NumbaTypeError(msg)

            out_none = cgutils.is_nonelike(out)
            if not out_none and not isinstance(out, types.Array):
                raise errors.NumbaTypeError('output must be an array')

            array_arr = isinstance(array, types.Array)
            array_ndim = array.ndim if array_arr else 0
            tup_m1 = (0,) * (array_ndim - 1)
            indices_arr = isinstance(indices, types.Array)
            dt = array.dtype if cgutils.is_nonelike(dtype) else dtype

            need_cast = not (array_arr and indices_arr)

            if indices_arr and indices.ndim != 1:
                msg = ('Expect "indices" array to have at most 1 dimension. '
                       f'Got {indices.ndim}')
                raise errors.NumbaValueError(msg)

            if need_cast:
                def impl_cast(ufunc, array, indices, axis=0, dtype=None, out=None):  # noqa: E501
                    return ufunc.reduceat(np.asarray(array),
                                          np.asarray(indices),
                                          axis,
                                          dtype=dtype,
                                          out=out)
                return impl_cast

            @register_jitable
            def tuple_slice(tup, pos):
                s = tup_m1
                i, j = 0, 0
                while i < len(tup):
                    if i != pos:
                        s = tuple_setitem(s, j, tup[i])
                        j += 1
                    i += 1
                return s

            # importing here as an import at the top scope brings unwanted
            # stuff. See numba/tests/test_import.py::test_no_impl_import
            from numba.np.arrayobj import generate_getitem_setitem_with_axis
            setitem = generate_getitem_setitem_with_axis(array.ndim, 'setitem')

            def impl(ufunc, array, indices, axis=0, dtype=None, out=None):
                sz = indices.shape[0]

                ax = axis
                if axis < 0:
                    axis += array_ndim

                if axis < 0 or axis >= array_ndim:
                    raise ValueError(f"axis {ax} is out of bounds for array "
                                     f"of dimension {array_ndim}")

                shape = tuple_setitem(array.shape, axis, sz)

                if not out_none and out.shape != shape:
                    # TODO: improve error message once #9524 gets merged
                    msg = ('operands could not be broadcast together with '
                           'remmaped shapes')
                    raise ValueError(msg)

                if out_none:
                    out = np.zeros(shape, dtype=dt)

                # short-circuit to avoid overflow on Windows
                if len(indices) == 0:
                    return out

                j = 0
                for i in range(len(indices) - 1):
                    if indices[i] < indices[i + 1]:
                        idx = np.arange(indices[i], indices[i + 1])
                        if array_ndim > 1:
                            arr_slice = np.take(array, idx, axis)
                        else:
                            arr_slice = array[idx]
                        arr_reduce = ufunc.reduce(arr_slice, axis=axis)
                        setitem(out, j, axis, arr_reduce)
                    elif indices[i] >= indices[i + 1]:
                        idx = indices[i]
                        arr_slice = np.take(array, idx, axis)
                        if array_ndim > 1:
                            _slice = tuple_slice(array.shape, axis)
                            arr_slice = arr_slice.reshape(_slice)
                        setitem(out, j, axis, arr_slice)
                    elif indices[i] >= sz or indices[i] < 0:
                        raise ValueError('Invalid value for indices')
                    j += 1

                # last index
                idx = np.arange(indices[-1], array.shape[axis])
                if array_ndim > 1:
                    arr_slice = np.take(array, idx, axis)
                else:
                    arr_slice = array[idx]
                arr_reduce = ufunc.reduce(arr_slice, axis)
                setitem(out, j, axis, arr_reduce)

                return out

            return impl

    def at(self, a, indices, b=None):
        # dynamic compile ufunc.at
        args = (a,) if cgutils.is_nonelike(b) else (a, b)
        argtys = (typeof(arg) for arg in args)
        ewise_types = tuple(arg.dtype if isinstance(arg, types.Array) else arg
                            for arg in argtys)

        if self.find_ewise_function(ewise_types) == (None, None):
            # cannot find a matching function and compilation is disabled
            if self._frozen:
                msg = "compilation disabled for %s.at(...)" % (self,)
                raise RuntimeError(msg)

            self._compile_for_args(*args)

        # all good, just dispatch to the function
        if cgutils.is_nonelike(b):
            return super().at(a, indices)
        else:
            return super().at(*(a, indices, b))

    def _install_type(self, typingctx=None):
        """Constructs and installs a typing class for a DUFunc object in the
        input typing context.  If no typing context is given, then
        _install_type() installs into the typing context of the
        dispatcher object (should be same default context used by
        jit() and njit()).
        """
        if typingctx is None:
            typingctx = self._dispatcher.targetdescr.typing_context
        _ty_cls = type('DUFuncTyping_' + self.ufunc.__name__,
                       (AbstractTemplate,),
                       dict(key=self, generic=self._type_me))
        typingctx.insert_user_function(self, _ty_cls)
        self._install_ufunc_attributes(_ty_cls)
        self._install_ufunc_methods(_ty_cls)

    def find_ewise_function(self, ewise_types):
        """
        Given a tuple of element-wise argument types, find a matching
        signature in the dispatcher.

        Return a 2-tuple containing the matching signature, and
        compilation result.  Will return two None's if no matching
        signature was found.
        """
        if self._frozen:
            # If we cannot compile, coerce to the best matching loop
            loop = numpy_support.ufunc_find_matching_loop(self, ewise_types)
            if loop is None:
                return None, None
            ewise_types = tuple(loop.inputs + loop.outputs)[:len(ewise_types)]
        for sig, cres in self._dispatcher.overloads.items():
            if sig.args == ewise_types:
                return sig, cres
        return None, None

    def _type_me(self, argtys, kwtys):
        """
        Implement AbstractTemplate.generic() for the typing class
        built by DUFunc._install_type().

        Return the call-site signature after either validating the
        element-wise signature or compiling for it.
        """
        assert not kwtys
        ufunc = self.ufunc
        _handle_inputs_result = npydecl.Numpy_rules_ufunc._handle_inputs(
            ufunc, argtys, kwtys)
        base_types, explicit_outputs, ndims, layout = _handle_inputs_result
        explicit_output_count = len(explicit_outputs)
        if explicit_output_count > 0:
            ewise_types = tuple(base_types[:-len(explicit_outputs)])
        else:
            ewise_types = tuple(base_types)
        sig, cres = self.find_ewise_function(ewise_types)
        if sig is None:
            # Matching element-wise signature was not found; must
            # compile.
            if self._frozen:
                raise errors.NumbaTypeError("cannot call %s with types %s"
                                            % (self, argtys))
            self._compile_for_argtys(ewise_types)
            sig, cres = self.find_ewise_function(ewise_types)
            assert sig is not None
        if explicit_output_count > 0:
            outtys = list(explicit_outputs)
        elif ufunc.nout == 1:
            if ndims > 0:
                outtys = [types.Array(sig.return_type, ndims, layout)]
            else:
                outtys = [sig.return_type]
        else:
            raise errors.NumbaNotImplementedError("typing gufuncs (nout > 1)")
        outtys.extend(argtys)
        return signature(*outtys)


array_analysis.MAP_TYPES.append(DUFunc)
