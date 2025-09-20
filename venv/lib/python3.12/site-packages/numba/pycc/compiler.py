# -*- coding: utf-8 -*-

import logging
import os
import sys

from llvmlite import ir
from llvmlite.binding import Linkage

from numba.pycc import llvm_types as lt
from numba.core.cgutils import create_constant_array
from numba.core.compiler import compile_extra, Flags
from numba.core.compiler_lock import global_compiler_lock

from numba.core.registry import cpu_target
from numba.core.runtime import nrtdynmod
from numba.core import cgutils
from numba.core.environment import lookup_environment


logger = logging.getLogger(__name__)

__all__ = ['Compiler']

NULL = ir.Constant(lt._void_star, None)
ZERO = ir.Constant(lt._int32, 0)
ONE = ir.Constant(lt._int32, 1)
METH_VARARGS_AND_KEYWORDS = ir.Constant(lt._int32, 1|2)


def get_header():
    import numpy
    import textwrap

    return textwrap.dedent("""\
    #include <stdint.h>

    #ifndef HAVE_LONGDOUBLE
        #define HAVE_LONGDOUBLE %d
    #endif

    typedef struct {
        float real;
        float imag;
    } complex64;

    typedef struct {
        double real;
        double imag;
    } complex128;

    #if HAVE_LONGDOUBLE
    typedef struct {
        long double real;
        long double imag;
    } complex256;
    #endif

    typedef float float32;
    typedef double float64;
    #if HAVE_LONGDOUBLE
    typedef long double float128;
    #endif
    """ % hasattr(numpy, 'complex256'))


class ExportEntry(object):
    """
    A simple record for exporting symbols.
    """

    def __init__(self, symbol, signature, function):
        self.symbol = symbol
        self.signature = signature
        self.function = function

    def __repr__(self):
        return "ExportEntry(%r, %r)" % (self.symbol, self.signature)


class _ModuleCompiler(object):
    """A base class to compile Python modules to a single shared library or
    extension module.

    :param export_entries: a list of ExportEntry instances.
    :param module_name: the name of the exported module.
    """

    #: Structure used to describe a method of an extension type.
    #: struct PyMethodDef {
    #:     const char  *ml_name;       /* The name of the built-in function/method */
    #:     PyCFunction  ml_meth;       /* The C function that implements it */
    #:     int          ml_flags;      /* Combination of METH_xxx flags, which mostly
    #:                                    describe the args expected by the C func */
    #:     const char  *ml_doc;        /* The __doc__ attribute, or NULL */
    #: };
    method_def_ty = ir.LiteralStructType((lt._int8_star,
                                          lt._void_star,
                                          lt._int32,
                                          lt._int8_star))

    method_def_ptr = ir.PointerType(method_def_ty)
    # The structure type constructed by PythonAPI.serialize_uncached()
    # when updating this, also make sure to update `env_def_t` struct in
    # numba/pycc/modulemixin.c
    env_def_ty = ir.LiteralStructType((lt._void_star,
                                       lt._int32,
                                       lt._void_star,
                                       lt._void_star,
                                       lt._int32))
    env_def_ptr = ir.PointerType(env_def_ty)

    def __init__(self, export_entries, module_name, use_nrt=False,
                 **aot_options):
        self.module_name = module_name
        self.export_python_wrap = False
        self.dll_exports = []
        self.export_entries = export_entries
        # Used by the CC API but not the legacy API
        self.external_init_function = None
        self.use_nrt = use_nrt

        self.typing_context = cpu_target.typing_context
        self.context = cpu_target.target_context.with_aot_codegen(
            self.module_name, **aot_options)

    def _mangle_method_symbol(self, func_name):
        return "._pycc_method_%s" % (func_name,)

    def _emit_python_wrapper(self, llvm_module):
        """Emit generated Python wrapper and extension module code.
        """
        raise NotImplementedError

    @global_compiler_lock
    def _cull_exports(self):
        """Read all the exported functions/modules in the translator
        environment, and join them into a single LLVM module.
        """
        self.exported_function_types = {}
        self.function_environments = {}
        self.environment_gvs = {}
        self.extra_environments = {}

        codegen = self.context.codegen()
        library = codegen.create_library(self.module_name)

        # Generate IR for all exported functions
        flags = Flags()
        flags.no_compile = True
        if not self.export_python_wrap:
            flags.no_cpython_wrapper = True
            flags.no_cfunc_wrapper = True
        if self.use_nrt:
            flags.nrt = True
            # Compile NRT helpers
            nrt_module, _ = nrtdynmod.create_nrt_module(self.context)
            library.add_ir_module(nrt_module)

        for entry in self.export_entries:
            cres = compile_extra(self.typing_context, self.context,
                                entry.function,
                                entry.signature.args,
                                entry.signature.return_type, flags,
                                locals={}, library=library)

            # Fix up dynamic exc globals
            module = library._final_module
            for gv in module.functions:
                if gv.name.startswith("__excinfo_unwrap_args"):
                    gv.linkage = "linkonce_odr"

            func_name = cres.fndesc.llvm_func_name
            llvm_func = cres.library.get_function(func_name)

            if self.export_python_wrap:
                llvm_func.linkage = 'internal'
                wrappername = cres.fndesc.llvm_cpython_wrapper_name
                wrapper = cres.library.get_function(wrappername)
                wrapper.name = self._mangle_method_symbol(entry.symbol)
                wrapper.linkage = 'external'
                fnty = cres.target_context.call_conv.get_function_type(
                    cres.fndesc.restype, cres.fndesc.argtypes)
                self.exported_function_types[entry] = fnty
                self.function_environments[entry] = cres.environment
                self.environment_gvs[entry] = cres.fndesc.env_name

                # Search for extra environments from linked libraries
                for linkedlib in library._linking_libraries:
                    linkedmod = linkedlib._final_module
                    # Find environments
                    for gv in linkedmod.global_variables:
                        gvn = gv.name
                        if gvn.startswith("_ZN08NumbaEnv"):
                            env = lookup_environment(gvn)
                            if env is not None:
                                if env.can_cache():
                                    self.extra_environments[gvn] = env
            else:
                llvm_func.name = entry.symbol
                self.dll_exports.append(entry.symbol)

        if self.export_python_wrap:
            wrapper_module = library.create_ir_module("wrapper")
            self._emit_python_wrapper(wrapper_module)
            library.add_ir_module(wrapper_module)

        # Hide all functions in the DLL except those explicitly exported
        library.finalize()
        for fn in library.get_defined_functions():
            if fn.name not in self.dll_exports:
                if fn.linkage in {Linkage.private, Linkage.internal}:
                    # Private/Internal linkage must have "default" visibility
                    fn.visibility = "default"
                else:
                    fn.visibility = 'hidden'
        return library

    def write_llvm_bitcode(self, output, wrap=False, **kws):
        self.export_python_wrap = wrap
        library = self._cull_exports()
        with open(output, 'wb') as fout:
            fout.write(library.emit_bitcode())

    def write_native_object(self, output, wrap=False, **kws):
        self.export_python_wrap = wrap
        library = self._cull_exports()
        with open(output, 'wb') as fout:
            fout.write(library.emit_native_object())

    def emit_type(self, tyobj):
        ret_val = str(tyobj)
        if 'int' in ret_val:
            if ret_val.endswith(('8', '16', '32', '64')):
                ret_val += "_t"
        return ret_val

    def emit_header(self, output):
        fname, ext = os.path.splitext(output)
        with open(fname + '.h', 'w') as fout:
            fout.write(get_header())
            fout.write("\n/* Prototypes */\n")
            for export_entry in self.export_entries:
                name = export_entry.symbol
                restype = self.emit_type(export_entry.signature.return_type)
                args = ", ".join(self.emit_type(argtype)
                                 for argtype in export_entry.signature.args)
                fout.write("extern %s %s(%s);\n" % (restype, name, args))

    def _emit_method_array(self, llvm_module):
        """
        Collect exported methods and emit a PyMethodDef array.

        :returns: a pointer to the PyMethodDef array.
        """
        method_defs = []
        for entry in self.export_entries:
            name = entry.symbol
            llvm_func_name = self._mangle_method_symbol(name)
            fnty = self.exported_function_types[entry]
            lfunc = ir.Function(llvm_module, fnty, llvm_func_name)

            method_name = self.context.insert_const_string(llvm_module, name)
            method_def_const = ir.Constant.literal_struct(
                (method_name,
                 ir.Constant.bitcast(lfunc, lt._void_star),
                 METH_VARARGS_AND_KEYWORDS,
                 NULL))
            method_defs.append(method_def_const)

        sentinel = ir.Constant.literal_struct([NULL, NULL, ZERO, NULL])
        method_defs.append(sentinel)
        method_array_init = create_constant_array(self.method_def_ty, method_defs)
        method_array = cgutils.add_global_variable(llvm_module,
                                                   method_array_init.type,
                                                   '.module_methods')
        method_array.initializer = method_array_init
        method_array.linkage = 'internal'
        method_array_ptr = ir.Constant.gep(method_array, [ZERO, ZERO])
        return method_array_ptr

    def _emit_environment_array(self, llvm_module, builder, pyapi):
        """
        Emit an array of env_def_t structures (see modulemixin.c)
        storing the pickled environment constants for each of the
        exported functions.
        """
        env_defs = []
        for entry in self.export_entries:
            env = self.function_environments[entry]
            # Constants may be unhashable so avoid trying to cache them
            env_def = pyapi.serialize_uncached(env.consts)
            env_defs.append(env_def)
        # Append extra environments
        env_defs.extend([
            pyapi.serialize_uncached(env.consts)
            for env in self.extra_environments.values()
        ])
        env_defs_init = create_constant_array(self.env_def_ty, env_defs)
        gv = self.context.insert_unique_const(llvm_module,
                                              '.module_environments',
                                              env_defs_init)
        return gv.gep([ZERO, ZERO])

    def _emit_envgvs_array(self, llvm_module, builder, pyapi):
        """
        Emit an array of Environment pointers that needs to be filled at
        initialization.
        """
        env_setters = []
        for entry in self.export_entries:
            envgv_name = self.environment_gvs[entry]
            gv = self.context.declare_env_global(llvm_module, envgv_name)
            envgv = gv.bitcast(lt._void_star)
            env_setters.append(envgv)
        # Append extra environments
        env_setters.extend([
            self.context.declare_env_global(
                llvm_module,
                envgv_name
            ).bitcast(lt._void_star)
            for envgv_name in self.extra_environments
        ])
        # The array ends with NULL
        env_setters.append(lt._void_star(None))

        env_setters_init = create_constant_array(lt._void_star, env_setters)
        gv = self.context.insert_unique_const(llvm_module,
                                              '.module_envgvs',
                                              env_setters_init)
        return gv.gep([ZERO, ZERO])

    def _emit_module_init_code(self, llvm_module, builder, modobj,
                               method_array, env_array, envgv_array):
        """
        Emit call to "external" init function, if any.
        """
        if self.external_init_function:
            fnty = ir.FunctionType(lt._int32,
                                   [modobj.type, self.method_def_ptr,
                                    self.env_def_ptr, envgv_array.type])
            fn = ir.Function(llvm_module, fnty, self.external_init_function)
            return builder.call(fn, [modobj, method_array, env_array,
                                     envgv_array])
        else:
            return None


class ModuleCompiler(_ModuleCompiler):

    _ptr_fun = lambda ret, *args: ir.PointerType(ir.FunctionType(ret, args))

    #: typedef int (*visitproc)(PyObject *, void *);
    visitproc_ty = _ptr_fun(lt._int8,
                            lt._pyobject_head_p)

    #: typedef int (*inquiry)(PyObject *);
    inquiry_ty = _ptr_fun(lt._int8,
                          lt._pyobject_head_p)

    #: typedef int (*traverseproc)(PyObject *, visitproc, void *);
    traverseproc_ty = _ptr_fun(lt._int8,
                               lt._pyobject_head_p,
                               visitproc_ty,
                               lt._void_star)

    #  typedef void (*freefunc)(void *)
    freefunc_ty = _ptr_fun(lt._int8,
                           lt._void_star)

    # PyObject* (*m_init)(void);
    m_init_ty = _ptr_fun(lt._int8)

    _char_star = lt._int8_star

    #: typedef struct PyModuleDef_Base {
    #:   PyObject_HEAD
    #:   PyObject* (*m_init)(void);
    #:   Py_ssize_t m_index;
    #:   PyObject* m_copy;
    #: } PyModuleDef_Base;
    module_def_base_ty = ir.LiteralStructType(
        (
            lt._pyobject_head,
            m_init_ty,
            lt._llvm_py_ssize_t,
            lt._pyobject_head_p
        ))

    #: This struct holds all information that is needed to create a module object.
    #: typedef struct PyModuleDef{
    #:   PyModuleDef_Base m_base;
    #:   const char* m_name;
    #:   const char* m_doc;
    #:   Py_ssize_t m_size;
    #:   PyMethodDef *m_methods;
    #:   inquiry m_reload;
    #:   traverseproc m_traverse;
    #:   inquiry m_clear;
    #:   freefunc m_free;
    #: }PyModuleDef;
    module_def_ty = ir.LiteralStructType(
        (
            module_def_base_ty,
            _char_star,
            _char_star,
            lt._llvm_py_ssize_t,
            _ModuleCompiler.method_def_ptr,
            inquiry_ty,
            traverseproc_ty,
            inquiry_ty,
            freefunc_ty
        ))

    @property
    def module_create_definition(self):
        """
        Return the signature and name of the Python C API function to
        initialize the module.
        """
        signature = ir.FunctionType(lt._pyobject_head_p,
                                    (ir.PointerType(self.module_def_ty),
                                     lt._int32))

        name = "PyModule_Create2"
        if lt._trace_refs_:
            name += "TraceRefs"

        return signature, name

    @property
    def module_init_definition(self):
        """
        Return the name and signature of the module's initialization function.
        """
        signature = ir.FunctionType(lt._pyobject_head_p, ())

        return signature, "PyInit_" + self.module_name

    def _emit_python_wrapper(self, llvm_module):
        # Figure out the Python C API module creation function, and
        # get a LLVM function for it.
        create_module_fn = ir.Function(llvm_module, *self.module_create_definition)
        create_module_fn.linkage = 'external'

        # Define a constant string for the module name.
        mod_name_const = self.context.insert_const_string(llvm_module,
                                                          self.module_name)

        mod_def_base_init = ir.Constant.literal_struct(
            (lt._pyobject_head_init,                        # PyObject_HEAD
             ir.Constant(self.m_init_ty, None),             # m_init
             ir.Constant(lt._llvm_py_ssize_t, None),        # m_index
             ir.Constant(lt._pyobject_head_p, None),        # m_copy
            )
        )
        mod_def_base = cgutils.add_global_variable(llvm_module,
                                                   mod_def_base_init.type,
                                                   '.module_def_base')
        mod_def_base.initializer = mod_def_base_init
        mod_def_base.linkage = 'internal'

        method_array = self._emit_method_array(llvm_module)

        mod_def_init = ir.Constant.literal_struct(
            (mod_def_base_init,                              # m_base
             mod_name_const,                                 # m_name
             ir.Constant(self._char_star, None),             # m_doc
             ir.Constant(lt._llvm_py_ssize_t, -1),           # m_size
             method_array,                                   # m_methods
             ir.Constant(self.inquiry_ty, None),             # m_reload
             ir.Constant(self.traverseproc_ty, None),        # m_traverse
             ir.Constant(self.inquiry_ty, None),             # m_clear
             ir.Constant(self.freefunc_ty, None)             # m_free
            )
        )

        # Define a constant string for the module name.
        mod_def = cgutils.add_global_variable(llvm_module, mod_def_init.type,
                                              '.module_def')
        mod_def.initializer = mod_def_init
        mod_def.linkage = 'internal'

        # Define the module initialization function.
        mod_init_fn = ir.Function(llvm_module, *self.module_init_definition)
        entry = mod_init_fn.append_basic_block('Entry')
        builder = ir.IRBuilder(entry)
        pyapi = self.context.get_python_api(builder)

        mod = builder.call(create_module_fn,
                           (mod_def,
                            ir.Constant(lt._int32, sys.api_version)))

        # Test if module has been created correctly.
        # (XXX for some reason comparing with the NULL constant fails llvm
        #  with an assertion in pydebug mode)
        with builder.if_then(cgutils.is_null(builder, mod)):
            builder.ret(NULL.bitcast(mod_init_fn.type.pointee.return_type))

        env_array = self._emit_environment_array(llvm_module, builder, pyapi)
        envgv_array = self._emit_envgvs_array(llvm_module, builder, pyapi)
        ret = self._emit_module_init_code(llvm_module, builder, mod,
                                          method_array, env_array, envgv_array)
        if ret is not None:
            with builder.if_then(cgutils.is_not_null(builder, ret)):
                # Init function errored out
                builder.ret(ir.Constant(mod.type, None))

        builder.ret(mod)

        self.dll_exports.append(mod_init_fn.name)

