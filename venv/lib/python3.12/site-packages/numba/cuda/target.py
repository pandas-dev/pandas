import re
from functools import cached_property
import llvmlite.binding as ll
from llvmlite import ir

from numba.core import (cgutils, config, debuginfo, itanium_mangler, types,
                        typing, utils)
from numba.core.dispatcher import Dispatcher
from numba.core.base import BaseContext
from numba.core.callconv import BaseCallConv, MinimalCallConv
from numba.core.typing import cmathdecl
from numba.core import datamodel

from .cudadrv import nvvm
from numba.cuda import codegen, nvvmutils, ufuncs
from numba.cuda.models import cuda_data_manager

# -----------------------------------------------------------------------------
# Typing


class CUDATypingContext(typing.BaseContext):
    def load_additional_registries(self):
        from . import cudadecl, cudamath, libdevicedecl, vector_types
        from numba.core.typing import enumdecl, cffi_utils

        self.install_registry(cudadecl.registry)
        self.install_registry(cffi_utils.registry)
        self.install_registry(cudamath.registry)
        self.install_registry(cmathdecl.registry)
        self.install_registry(libdevicedecl.registry)
        self.install_registry(enumdecl.registry)
        self.install_registry(vector_types.typing_registry)

    def resolve_value_type(self, val):
        # treat other dispatcher object as another device function
        from numba.cuda.dispatcher import CUDADispatcher
        if (isinstance(val, Dispatcher) and not
                isinstance(val, CUDADispatcher)):
            try:
                # use cached device function
                val = val.__dispatcher
            except AttributeError:
                if not val._can_compile:
                    raise ValueError('using cpu function on device '
                                     'but its compilation is disabled')
                targetoptions = val.targetoptions.copy()
                targetoptions['device'] = True
                targetoptions['debug'] = targetoptions.get('debug', False)
                targetoptions['opt'] = targetoptions.get('opt', True)
                disp = CUDADispatcher(val.py_func, targetoptions)
                # cache the device function for future use and to avoid
                # duplicated copy of the same function.
                val.__dispatcher = disp
                val = disp

        # continue with parent logic
        return super(CUDATypingContext, self).resolve_value_type(val)

# -----------------------------------------------------------------------------
# Implementation


VALID_CHARS = re.compile(r'[^a-z0-9]', re.I)


class CUDATargetContext(BaseContext):
    implement_powi_as_math_call = True
    strict_alignment = True

    def __init__(self, typingctx, target='cuda'):
        super().__init__(typingctx, target)
        self.data_model_manager = cuda_data_manager.chain(
            datamodel.default_manager
        )

    @property
    def DIBuilder(self):
        return debuginfo.DIBuilder

    @property
    def enable_boundscheck(self):
        # Unconditionally disabled
        return False

    # Overrides
    def create_module(self, name):
        return self._internal_codegen._create_empty_module(name)

    def init(self):
        self._internal_codegen = codegen.JITCUDACodegen("numba.cuda.jit")
        self._target_data = None

    def load_additional_registries(self):
        # side effect of import needed for numba.cpython.*, the builtins
        # registry is updated at import time.
        from numba.cpython import numbers, tupleobj, slicing # noqa: F401
        from numba.cpython import rangeobj, iterators, enumimpl # noqa: F401
        from numba.cpython import unicode, charseq # noqa: F401
        from numba.cpython import cmathimpl
        from numba.misc import cffiimpl
        from numba.np import arrayobj # noqa: F401
        from numba.np import npdatetime # noqa: F401
        from . import (
            cudaimpl, printimpl, libdeviceimpl, mathimpl, vector_types
        )
        # fix for #8940
        from numba.np.unsafe import ndarray # noqa F401

        self.install_registry(cudaimpl.registry)
        self.install_registry(cffiimpl.registry)
        self.install_registry(printimpl.registry)
        self.install_registry(libdeviceimpl.registry)
        self.install_registry(cmathimpl.registry)
        self.install_registry(mathimpl.registry)
        self.install_registry(vector_types.impl_registry)

    def codegen(self):
        return self._internal_codegen

    @property
    def target_data(self):
        if self._target_data is None:
            self._target_data = ll.create_target_data(nvvm.NVVM().data_layout)
        return self._target_data

    @cached_property
    def nonconst_module_attrs(self):
        """
        Some CUDA intrinsics are at the module level, but cannot be treated as
        constants, because they are loaded from a special register in the PTX.
        These include threadIdx, blockDim, etc.
        """
        from numba import cuda
        nonconsts = ('threadIdx', 'blockDim', 'blockIdx', 'gridDim', 'laneid',
                     'warpsize')
        nonconsts_with_mod = tuple([(types.Module(cuda), nc)
                                    for nc in nonconsts])
        return nonconsts_with_mod

    @cached_property
    def call_conv(self):
        return CUDACallConv(self)

    def mangler(self, name, argtypes, *, abi_tags=(), uid=None):
        return itanium_mangler.mangle(name, argtypes, abi_tags=abi_tags,
                                      uid=uid)

    def prepare_cuda_kernel(self, codelib, fndesc, debug, lineinfo,
                            nvvm_options, filename, linenum,
                            max_registers=None):
        """
        Adapt a code library ``codelib`` with the numba compiled CUDA kernel
        with name ``fname`` and arguments ``argtypes`` for NVVM.
        A new library is created with a wrapper function that can be used as
        the kernel entry point for the given kernel.

        Returns the new code library and the wrapper function.

        Parameters:

        codelib:       The CodeLibrary containing the device function to wrap
                       in a kernel call.
        fndesc:        The FunctionDescriptor of the source function.
        debug:         Whether to compile with debug.
        lineinfo:      Whether to emit line info.
        nvvm_options:  Dict of NVVM options used when compiling the new library.
        filename:      The source filename that the function is contained in.
        linenum:       The source line that the function is on.
        max_registers: The max_registers argument for the code library.
        """
        kernel_name = itanium_mangler.prepend_namespace(
            fndesc.llvm_func_name, ns='cudapy',
        )
        library = self.codegen().create_library(f'{codelib.name}_kernel_',
                                                entry_name=kernel_name,
                                                nvvm_options=nvvm_options,
                                                max_registers=max_registers)
        library.add_linking_library(codelib)
        wrapper = self.generate_kernel_wrapper(library, fndesc, kernel_name,
                                               debug, lineinfo, filename,
                                               linenum)
        return library, wrapper

    def generate_kernel_wrapper(self, library, fndesc, kernel_name, debug,
                                lineinfo, filename, linenum):
        """
        Generate the kernel wrapper in the given ``library``.
        The function being wrapped is described by ``fndesc``.
        The wrapper function is returned.
        """

        argtypes = fndesc.argtypes
        arginfo = self.get_arg_packer(argtypes)
        argtys = list(arginfo.argument_types)
        wrapfnty = ir.FunctionType(ir.VoidType(), argtys)
        wrapper_module = self.create_module("cuda.kernel.wrapper")
        fnty = ir.FunctionType(ir.IntType(32),
                               [self.call_conv.get_return_type(types.pyobject)]
                               + argtys)
        func = ir.Function(wrapper_module, fnty, fndesc.llvm_func_name)

        prefixed = itanium_mangler.prepend_namespace(func.name, ns='cudapy')
        wrapfn = ir.Function(wrapper_module, wrapfnty, prefixed)
        builder = ir.IRBuilder(wrapfn.append_basic_block(''))

        if debug or lineinfo:
            directives_only = lineinfo and not debug
            debuginfo = self.DIBuilder(module=wrapper_module,
                                       filepath=filename,
                                       cgctx=self,
                                       directives_only=directives_only)
            debuginfo.mark_subprogram(
                wrapfn, kernel_name, fndesc.args, argtypes, linenum,
            )
            debuginfo.mark_location(builder, linenum)

        # Define error handling variable
        def define_error_gv(postfix):
            name = wrapfn.name + postfix
            gv = cgutils.add_global_variable(wrapper_module, ir.IntType(32),
                                             name)
            gv.initializer = ir.Constant(gv.type.pointee, None)
            return gv

        gv_exc = define_error_gv("__errcode__")
        gv_tid = []
        gv_ctaid = []
        for i in 'xyz':
            gv_tid.append(define_error_gv("__tid%s__" % i))
            gv_ctaid.append(define_error_gv("__ctaid%s__" % i))

        callargs = arginfo.from_arguments(builder, wrapfn.args)
        status, _ = self.call_conv.call_function(
            builder, func, types.void, argtypes, callargs)

        if debug:
            # Check error status
            with cgutils.if_likely(builder, status.is_ok):
                builder.ret_void()

            with builder.if_then(builder.not_(status.is_python_exc)):
                # User exception raised
                old = ir.Constant(gv_exc.type.pointee, None)

                # Use atomic cmpxchg to prevent rewriting the error status
                # Only the first error is recorded

                xchg = builder.cmpxchg(gv_exc, old, status.code,
                                       'monotonic', 'monotonic')
                changed = builder.extract_value(xchg, 1)

                # If the xchange is successful, save the thread ID.
                sreg = nvvmutils.SRegBuilder(builder)
                with builder.if_then(changed):
                    for dim, ptr, in zip("xyz", gv_tid):
                        val = sreg.tid(dim)
                        builder.store(val, ptr)

                    for dim, ptr, in zip("xyz", gv_ctaid):
                        val = sreg.ctaid(dim)
                        builder.store(val, ptr)

        builder.ret_void()

        nvvm.set_cuda_kernel(wrapfn)
        library.add_ir_module(wrapper_module)
        if debug or lineinfo:
            debuginfo.finalize()
        library.finalize()

        if config.DUMP_LLVM:
            utils.dump_llvm(fndesc, wrapper_module)

        return library.get_function(wrapfn.name)

    def make_constant_array(self, builder, aryty, arr):
        """
        Unlike the parent version.  This returns a a pointer in the constant
        addrspace.
        """

        lmod = builder.module

        constvals = [
            self.get_constant(types.byte, i)
            for i in iter(arr.tobytes(order='A'))
        ]
        constaryty = ir.ArrayType(ir.IntType(8), len(constvals))
        constary = ir.Constant(constaryty, constvals)

        addrspace = nvvm.ADDRSPACE_CONSTANT
        gv = cgutils.add_global_variable(lmod, constary.type, "_cudapy_cmem",
                                         addrspace=addrspace)
        gv.linkage = 'internal'
        gv.global_constant = True
        gv.initializer = constary

        # Preserve the underlying alignment
        lldtype = self.get_data_type(aryty.dtype)
        align = self.get_abi_sizeof(lldtype)
        gv.align = 2 ** (align - 1).bit_length()

        # Convert to generic address-space
        ptrty = ir.PointerType(ir.IntType(8))
        genptr = builder.addrspacecast(gv, ptrty, 'generic')

        # Create array object
        ary = self.make_array(aryty)(self, builder)
        kshape = [self.get_constant(types.intp, s) for s in arr.shape]
        kstrides = [self.get_constant(types.intp, s) for s in arr.strides]
        self.populate_array(ary, data=builder.bitcast(genptr, ary.data.type),
                            shape=kshape,
                            strides=kstrides,
                            itemsize=ary.itemsize, parent=ary.parent,
                            meminfo=None)

        return ary._getvalue()

    def insert_const_string(self, mod, string):
        """
        Unlike the parent version.  This returns a a pointer in the constant
        addrspace.
        """
        text = cgutils.make_bytearray(string.encode("utf-8") + b"\x00")
        name = '$'.join(["__conststring__",
                         itanium_mangler.mangle_identifier(string)])
        # Try to reuse existing global
        gv = mod.globals.get(name)
        if gv is None:
            # Not defined yet
            gv = cgutils.add_global_variable(mod, text.type, name,
                                             addrspace=nvvm.ADDRSPACE_CONSTANT)
            gv.linkage = 'internal'
            gv.global_constant = True
            gv.initializer = text

        # Cast to a i8* pointer
        charty = gv.type.pointee.element
        return gv.bitcast(charty.as_pointer(nvvm.ADDRSPACE_CONSTANT))

    def insert_string_const_addrspace(self, builder, string):
        """
        Insert a constant string in the constant addresspace and return a
        generic i8 pointer to the data.

        This function attempts to deduplicate.
        """
        lmod = builder.module
        gv = self.insert_const_string(lmod, string)
        charptrty = ir.PointerType(ir.IntType(8))
        return builder.addrspacecast(gv, charptrty, 'generic')

    def optimize_function(self, func):
        """Run O1 function passes
        """
        pass
        ## XXX skipped for now
        # fpm = lp.FunctionPassManager.new(func.module)
        #
        # lp.PassManagerBuilder.new().populate(fpm)
        #
        # fpm.initialize()
        # fpm.run(func)
        # fpm.finalize()

    def get_ufunc_info(self, ufunc_key):
        return ufuncs.get_ufunc_info(ufunc_key)


class CUDACallConv(MinimalCallConv):
    pass


class CUDACABICallConv(BaseCallConv):
    """
    Calling convention aimed at matching the CUDA C/C++ ABI. The implemented
    function signature is:

        <Python return type> (<Python arguments>)

    Exceptions are unsupported in this convention.
    """

    def _make_call_helper(self, builder):
        # Call helpers are used to help report exceptions back to Python, so
        # none is required here.
        return None

    def return_value(self, builder, retval):
        return builder.ret(retval)

    def return_user_exc(self, builder, exc, exc_args=None, loc=None,
                        func_name=None):
        msg = "Python exceptions are unsupported in the CUDA C/C++ ABI"
        raise NotImplementedError(msg)

    def return_status_propagate(self, builder, status):
        msg = "Return status is unsupported in the CUDA C/C++ ABI"
        raise NotImplementedError(msg)

    def get_function_type(self, restype, argtypes):
        """
        Get the LLVM IR Function type for *restype* and *argtypes*.
        """
        arginfo = self._get_arg_packer(argtypes)
        argtypes = list(arginfo.argument_types)
        fnty = ir.FunctionType(self.get_return_type(restype), argtypes)
        return fnty

    def decorate_function(self, fn, args, fe_argtypes, noalias=False):
        """
        Set names and attributes of function arguments.
        """
        assert not noalias
        arginfo = self._get_arg_packer(fe_argtypes)
        arginfo.assign_names(self.get_arguments(fn),
                             ['arg.' + a for a in args])

    def get_arguments(self, func):
        """
        Get the Python-level arguments of LLVM *func*.
        """
        return func.args

    def call_function(self, builder, callee, resty, argtys, args):
        """
        Call the Numba-compiled *callee*.
        """
        arginfo = self._get_arg_packer(argtys)
        realargs = arginfo.as_arguments(builder, args)
        code = builder.call(callee, realargs)
        # No status required as we don't support exceptions or a distinct None
        # value in a C ABI.
        status = None
        out = self.context.get_returned_value(builder, resty, code)
        return status, out

    def get_return_type(self, ty):
        return self.context.data_model_manager[ty].get_return_type()
