import platform

import llvmlite.binding as ll
from llvmlite import ir

from numba import _dynfunc
from numba.core.callwrapper import PyCallWrapper
from numba.core.base import BaseContext
from numba.core import (utils, types, config, cgutils, callconv, codegen,
                        externals, fastmathpass, intrinsics)
from numba.core.options import TargetOptions, include_default_options
from numba.core.runtime import rtsys
from numba.core.compiler_lock import global_compiler_lock
import numba.core.entrypoints
# Re-export these options, they are used from the cpu module throughout the code
# base.
from numba.core.cpu_options import (ParallelOptions, # noqa F401
                                    FastMathOptions, InlineOptions) # noqa F401
from numba.np import ufunc_db

# Keep those structures in sync with _dynfunc.c.


class ClosureBody(cgutils.Structure):
    _fields = [('env', types.pyobject)]


class EnvBody(cgutils.Structure):
    _fields = [
        ('globals', types.pyobject),
        ('consts', types.pyobject),
    ]


class CPUContext(BaseContext):
    """
    Changes BaseContext calling convention
    """
    allow_dynamic_globals = True

    def __init__(self, typingctx, target='cpu'):
        super().__init__(typingctx, target)

    # Overrides
    def create_module(self, name):
        return self._internal_codegen._create_empty_module(name)

    @global_compiler_lock
    def init(self):
        self.is32bit = (utils.MACHINE_BITS == 32)
        self._internal_codegen = codegen.JITCPUCodegen("numba.exec")

        # Add ARM ABI functions from libgcc_s
        if platform.machine() == 'armv7l':
            ll.load_library_permanently('libgcc_s.so.1')

        # Map external C functions.
        externals.c_math_functions.install(self)

    def load_additional_registries(self):
        # Only initialize the NRT once something is about to be compiled. The
        # "initialized" state doesn't need to be threadsafe, there's a lock
        # around the internal compilation and the rtsys.initialize call can be
        # made multiple times, worse case init just gets called a bit more often
        # than optimal.
        rtsys.initialize(self)

        # Add implementations that work via import
        from numba.cpython import (builtins, charseq, enumimpl, # noqa F401
                                   hashing, heapq, iterators, # noqa F401
                                   listobj, numbers, rangeobj, # noqa F401
                                   setobj, slicing, tupleobj, # noqa F401
                                   unicode,) # noqa F401
        from numba.core import optional, inline_closurecall # noqa F401
        from numba.misc import gdb_hook, literal # noqa F401
        from numba.np import linalg, arraymath, arrayobj # noqa F401
        from numba.np.random import generator_core, generator_methods # noqa F401
        from numba.np.polynomial import polynomial_core, polynomial_functions # noqa F401
        from numba.typed import typeddict, dictimpl # noqa F401
        from numba.typed import typedlist, listobject # noqa F401
        from numba.experimental import jitclass, function_type # noqa F401
        from numba.np import npdatetime # noqa F401

        # Add target specific implementations
        from numba.np import npyimpl
        from numba.cpython import cmathimpl, mathimpl, printimpl, randomimpl
        from numba.misc import cffiimpl
        from numba.experimental.jitclass.base import ClassBuilder as \
            jitclassimpl
        self.install_registry(cmathimpl.registry)
        self.install_registry(cffiimpl.registry)
        self.install_registry(mathimpl.registry)
        self.install_registry(npyimpl.registry)
        self.install_registry(printimpl.registry)
        self.install_registry(randomimpl.registry)
        self.install_registry(jitclassimpl.class_impl_registry)

        # load 3rd party extensions
        numba.core.entrypoints.init_all()

        # fix for #8940
        from numba.np.unsafe import ndarray # noqa F401

    @property
    def target_data(self):
        return self._internal_codegen.target_data

    def with_aot_codegen(self, name, **aot_options):
        aot_codegen = codegen.AOTCPUCodegen(name, **aot_options)
        return self.subtarget(_internal_codegen=aot_codegen,
                              aot_mode=True)

    def codegen(self):
        return self._internal_codegen

    @property
    def call_conv(self):
        return callconv.CPUCallConv(self)

    def get_env_body(self, builder, envptr):
        """
        From the given *envptr* (a pointer to a _dynfunc.Environment object),
        get a EnvBody allowing structured access to environment fields.
        """
        body_ptr = cgutils.pointer_add(
            builder, envptr, _dynfunc._impl_info['offsetof_env_body'])
        return EnvBody(self, builder, ref=body_ptr, cast_ref=True)

    def get_env_manager(self, builder, return_pyobject=False):
        envgv = self.declare_env_global(builder.module,
                                        self.get_env_name(self.fndesc))
        envarg = builder.load(envgv)
        pyapi = self.get_python_api(builder)
        pyapi.emit_environment_sentry(
            envarg,
            return_pyobject=return_pyobject,
            debug_msg=self.fndesc.env_name,
        )
        env_body = self.get_env_body(builder, envarg)
        return pyapi.get_env_manager(self.environment, env_body, envarg)

    def get_generator_state(self, builder, genptr, return_type):
        """
        From the given *genptr* (a pointer to a _dynfunc.Generator object),
        get a pointer to its state area.
        """
        return cgutils.pointer_add(
            builder, genptr, _dynfunc._impl_info['offsetof_generator_state'],
            return_type=return_type)

    def build_list(self, builder, list_type, items):
        """
        Build a list from the Numba *list_type* and its initial *items*.
        """
        from numba.cpython import listobj
        return listobj.build_list(self, builder, list_type, items)

    def build_set(self, builder, set_type, items):
        """
        Build a set from the Numba *set_type* and its initial *items*.
        """
        from numba.cpython import setobj
        return setobj.build_set(self, builder, set_type, items)

    def build_map(self, builder, dict_type, item_types, items):
        from numba.typed import dictobject

        return dictobject.build_map(self, builder, dict_type, item_types, items)

    def post_lowering(self, mod, library):
        if self.fastmath:
            fastmathpass.rewrite_module(mod, self.fastmath)

        if self.is32bit:
            # 32-bit machine needs to replace all 64-bit div/rem to avoid
            # calls to compiler-rt
            intrinsics.fix_divmod(mod)

        library.add_linking_library(rtsys.library)

    def create_cpython_wrapper(self, library, fndesc, env, call_helper,
                               release_gil=False):
        wrapper_module = self.create_module("wrapper")
        fnty = self.call_conv.get_function_type(fndesc.restype, fndesc.argtypes)
        wrapper_callee = ir.Function(wrapper_module, fnty,
                                     fndesc.llvm_func_name)
        builder = PyCallWrapper(self, wrapper_module, wrapper_callee,
                                fndesc, env, call_helper=call_helper,
                                release_gil=release_gil)
        builder.build()
        library.add_ir_module(wrapper_module)

    def create_cfunc_wrapper(self, library, fndesc, env, call_helper):
        wrapper_module = self.create_module("cfunc_wrapper")
        fnty = self.call_conv.get_function_type(fndesc.restype, fndesc.argtypes)
        wrapper_callee = ir.Function(wrapper_module, fnty,
                                     fndesc.llvm_func_name)

        ll_argtypes = [self.get_value_type(ty) for ty in fndesc.argtypes]
        ll_return_type = self.get_value_type(fndesc.restype)
        wrapty = ir.FunctionType(ll_return_type, ll_argtypes)
        wrapfn = ir.Function(wrapper_module, wrapty,
                             fndesc.llvm_cfunc_wrapper_name)
        builder = ir.IRBuilder(wrapfn.append_basic_block('entry'))

        status, out = self.call_conv.call_function(
            builder, wrapper_callee, fndesc.restype, fndesc.argtypes,
            wrapfn.args, attrs=('noinline',))

        with builder.if_then(status.is_error, likely=False):
            # If (and only if) an error occurred, acquire the GIL
            # and use the interpreter to write out the exception.
            pyapi = self.get_python_api(builder)
            gil_state = pyapi.gil_ensure()
            self.call_conv.raise_error(builder, pyapi, status)
            cstr = self.insert_const_string(builder.module, repr(self))
            strobj = pyapi.string_from_string(cstr)
            pyapi.err_write_unraisable(strobj)
            pyapi.decref(strobj)
            pyapi.gil_release(gil_state)

        builder.ret(out)
        library.add_ir_module(wrapper_module)

    def get_executable(self, library, fndesc, env):
        """
        Returns
        -------
        (cfunc, fnptr)

        - cfunc
            callable function (Can be None)
        - fnptr
            callable function address
        - env
            an execution environment (from _dynfunc)
        """
        # Code generation
        fnptr = library.get_pointer_to_function(
            fndesc.llvm_cpython_wrapper_name)

        # Note: we avoid reusing the original docstring to avoid encoding
        # issues on Python 2, see issue #1908
        doc = "compiled wrapper for %r" % (fndesc.qualname,)
        cfunc = _dynfunc.make_function(fndesc.lookup_module(),
                                       fndesc.qualname.split('.')[-1],
                                       doc, fnptr, env,
                                       # objects to keepalive with the function
                                       (library,)
                                       )
        library.codegen.set_env(self.get_env_name(fndesc), env)
        return cfunc

    def calc_array_sizeof(self, ndim):
        '''
        Calculate the size of an array struct on the CPU target
        '''
        aryty = types.Array(types.int32, ndim, 'A')
        return self.get_abi_sizeof(self.get_value_type(aryty))

    # Overrides
    def get_ufunc_info(self, ufunc_key):
        return ufunc_db.get_ufunc_info(ufunc_key)


# ----------------------------------------------------------------------------
# TargetOptions

_options_mixin = include_default_options(
    "nopython",
    "forceobj",
    "looplift",
    "_nrt",
    "debug",
    "boundscheck",
    "nogil",
    "no_rewrites",
    "no_cpython_wrapper",
    "no_cfunc_wrapper",
    "parallel",
    "fastmath",
    "error_model",
    "inline",
    "forceinline",
    "_dbg_extend_lifetimes",
    "_dbg_optnone",
)


class CPUTargetOptions(_options_mixin, TargetOptions):
    def finalize(self, flags, options):
        if not flags.is_set("enable_pyobject"):
            flags.enable_pyobject = True

        if not flags.is_set("enable_looplift"):
            flags.enable_looplift = True

        flags.inherit_if_not_set("nrt", default=True)

        if not flags.is_set("debuginfo"):
            flags.debuginfo = config.DEBUGINFO_DEFAULT

        if not flags.is_set("dbg_extend_lifetimes"):
            if flags.debuginfo:
                # auto turn on extend-lifetimes if debuginfo is on and
                # dbg_extend_lifetimes is not set
                flags.dbg_extend_lifetimes = True
            else:
                # set flag using env-var config
                flags.dbg_extend_lifetimes = config.EXTEND_VARIABLE_LIFETIMES

        if not flags.is_set("boundscheck"):
            flags.boundscheck = flags.debuginfo

        flags.enable_pyobject_looplift = True

        flags.inherit_if_not_set("fastmath")

        flags.inherit_if_not_set("error_model", default="python")

        flags.inherit_if_not_set("forceinline")

        if flags.forceinline:
            # forceinline turns off optnone, just like clang.
            flags.dbg_optnone = False
