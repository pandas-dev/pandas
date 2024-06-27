from collections import defaultdict
import copy
import sys
from itertools import permutations, takewhile
from contextlib import contextmanager
from functools import cached_property

from llvmlite import ir as llvmir
from llvmlite.ir import Constant
import llvmlite.binding as ll

from numba.core import types, utils, datamodel, debuginfo, funcdesc, config, cgutils, imputils
from numba.core import event, errors, targetconfig
from numba import _dynfunc, _helperlib
from numba.core.compiler_lock import global_compiler_lock
from numba.core.pythonapi import PythonAPI
from numba.core.imputils import (user_function, user_generator,
                       builtin_registry, impl_ret_borrowed,
                       RegistryLoader)
from numba.cpython import builtins

GENERIC_POINTER = llvmir.PointerType(llvmir.IntType(8))
PYOBJECT = GENERIC_POINTER
void_ptr = GENERIC_POINTER


class OverloadSelector(object):
    """
    An object matching an actual signature against a registry of formal
    signatures and choosing the best candidate, if any.

    In the current implementation:
    - a "signature" is a tuple of type classes or type instances
    - the "best candidate" is the most specific match
    """

    def __init__(self):
        # A list of (formal args tuple, value)
        self.versions = []
        self._cache = {}

    def find(self, sig):
        out = self._cache.get(sig)
        if out is None:
            out = self._find(sig)
            self._cache[sig] = out
        return out

    def _find(self, sig):
        candidates = self._select_compatible(sig)
        if candidates:
            return candidates[self._best_signature(candidates)]
        else:
            raise errors.NumbaNotImplementedError(f'{self}, {sig}')

    def _select_compatible(self, sig):
        """
        Select all compatible signatures and their implementation.
        """
        out = {}
        for ver_sig, impl in self.versions:
            if self._match_arglist(ver_sig, sig):
                out[ver_sig] = impl
        return out

    def _best_signature(self, candidates):
        """
        Returns the best signature out of the candidates
        """
        ordered, genericity = self._sort_signatures(candidates)
        # check for ambiguous signatures
        if len(ordered) > 1:
            firstscore = genericity[ordered[0]]
            same = list(takewhile(lambda x: genericity[x] == firstscore,
                                  ordered))
            if len(same) > 1:
                msg = ["{n} ambiguous signatures".format(n=len(same))]
                for sig in same:
                    msg += ["{0} => {1}".format(sig, candidates[sig])]
                raise errors.NumbaTypeError('\n'.join(msg))
        return ordered[0]

    def _sort_signatures(self, candidates):
        """
        Sort signatures in ascending level of genericity.

        Returns a 2-tuple:

            * ordered list of signatures
            * dictionary containing genericity scores
        """
        # score by genericity
        genericity = defaultdict(int)
        for this, other in permutations(candidates.keys(), r=2):
            matched = self._match_arglist(formal_args=this, actual_args=other)
            if matched:
                # genericity score +1 for every another compatible signature
                genericity[this] += 1
        # order candidates in ascending level of genericity
        ordered = sorted(candidates.keys(), key=lambda x: genericity[x])
        return ordered, genericity

    def _match_arglist(self, formal_args, actual_args):
        """
        Returns True if the signature is "matching".
        A formal signature is "matching" if the actual signature matches exactly
        or if the formal signature is a compatible generic signature.
        """
        # normalize VarArg
        if formal_args and isinstance(formal_args[-1], types.VarArg):
            ndiff = len(actual_args) - len(formal_args) + 1
            formal_args = formal_args[:-1] + (formal_args[-1].dtype,) * ndiff

        if len(formal_args) != len(actual_args):
            return False

        for formal, actual in zip(formal_args, actual_args):
            if not self._match(formal, actual):
                return False

        return True

    def _match(self, formal, actual):
        if formal == actual:
            # formal argument matches actual arguments
            return True
        elif types.Any == formal:
            # formal argument is any
            return True
        elif isinstance(formal, type) and issubclass(formal, types.Type):
            if isinstance(actual, type) and issubclass(actual, formal):
                # formal arg is a type class and actual arg is a subclass
                return True
            elif isinstance(actual, formal):
                # formal arg is a type class of which actual arg is an instance
                return True

    def append(self, value, sig):
        """
        Add a formal signature and its associated value.
        """
        assert isinstance(sig, tuple), (value, sig)
        self.versions.append((sig, value))
        self._cache.clear()


@utils.runonce
def _load_global_helpers():
    """
    Execute once to install special symbols into the LLVM symbol table.
    """
    # This is Py_None's real C name
    ll.add_symbol("_Py_NoneStruct", id(None))

    # Add Numba C helper functions
    for c_helpers in (_helperlib.c_helpers, _dynfunc.c_helpers):
        for py_name, c_address in c_helpers.items():
            c_name = "numba_" + py_name
            ll.add_symbol(c_name, c_address)

    # Add all built-in exception classes
    for obj in utils.builtins.__dict__.values():
        if isinstance(obj, type) and issubclass(obj, BaseException):
            ll.add_symbol("PyExc_%s" % (obj.__name__), id(obj))


class BaseContext(object):
    """

    Notes on Structure
    ------------------

    Most objects are lowered as plain-old-data structure in the generated
    llvm.  They are passed around by reference (a pointer to the structure).
    Only POD structure can live across function boundaries by copying the
    data.
    """
    # True if the target requires strict alignment
    # Causes exception to be raised if the record members are not aligned.
    strict_alignment = False

    # Force powi implementation as math.pow call
    implement_powi_as_math_call = False
    implement_pow_as_math_call = False

    # Emit Debug info
    enable_debuginfo = False
    DIBuilder = debuginfo.DIBuilder

    # Bound checking
    @property
    def enable_boundscheck(self):
        if config.BOUNDSCHECK is not None:
            return config.BOUNDSCHECK
        return self._boundscheck

    @enable_boundscheck.setter
    def enable_boundscheck(self, value):
        self._boundscheck = value

    # NRT
    enable_nrt = False

    # Auto parallelization
    auto_parallel = False

    # PYCC
    aot_mode = False

    # Error model for various operations (only FP exceptions currently)
    error_model = None

    # Whether dynamic globals (CPU runtime addresses) is allowed
    allow_dynamic_globals = False

    # Fast math flags
    fastmath = False

    # python execution environment
    environment = None

    # the function descriptor
    fndesc = None

    def __init__(self, typing_context, target):
        _load_global_helpers()

        self.address_size = utils.MACHINE_BITS
        self.typing_context = typing_context
        from numba.core.target_extension import target_registry
        self.target_name = target
        self.target = target_registry[target]

        # A mapping of installed registries to their loaders
        self._registries = {}
        # Declarations loaded from registries and other sources
        self._defns = defaultdict(OverloadSelector)
        self._getattrs = defaultdict(OverloadSelector)
        self._setattrs = defaultdict(OverloadSelector)
        self._casts = OverloadSelector()
        self._get_constants = OverloadSelector()
        # Other declarations
        self._generators = {}
        self.special_ops = {}
        self.cached_internal_func = {}
        self._pid = None
        self._codelib_stack = []

        self._boundscheck = False

        self.data_model_manager = datamodel.default_manager

        # Initialize
        self.init()

    def init(self):
        """
        For subclasses to add initializer
        """

    def refresh(self):
        """
        Refresh context with new declarations from known registries.
        Useful for third-party extensions.
        """
        # load target specific registries
        self.load_additional_registries()

        # Populate the builtin registry, this has to happen after loading
        # additional registries as some of the "additional" registries write
        # their implementations into the builtin_registry and would be missed if
        # this ran first.
        self.install_registry(builtin_registry)

        # Also refresh typing context, since @overload declarations can
        # affect it.
        self.typing_context.refresh()

    def load_additional_registries(self):
        """
        Load target-specific registries.  Can be overridden by subclasses.
        """

    def mangler(self, name, types, *, abi_tags=(), uid=None):
        """
        Perform name mangling.
        """
        return funcdesc.default_mangler(name, types, abi_tags=abi_tags, uid=uid)

    def get_env_name(self, fndesc):
        """Get the environment name given a FunctionDescriptor.

        Use this instead of the ``fndesc.env_name`` so that the target-context
        can provide necessary mangling of the symbol to meet ABI requirements.
        """
        return fndesc.env_name

    def declare_env_global(self, module, envname):
        """Declare the Environment pointer as a global of the module.

        The pointer is initialized to NULL.  It must be filled by the runtime
        with the actual address of the Env before the associated function
        can be executed.

        Parameters
        ----------
        module :
            The LLVM Module
        envname : str
            The name of the global variable.
        """
        if envname not in module.globals:
            gv = llvmir.GlobalVariable(module, cgutils.voidptr_t, name=envname)
            gv.linkage = 'common'
            gv.initializer = cgutils.get_null_value(gv.type.pointee)

        return module.globals[envname]

    def get_arg_packer(self, fe_args):
        return datamodel.ArgPacker(self.data_model_manager, fe_args)

    def get_data_packer(self, fe_types):
        return datamodel.DataPacker(self.data_model_manager, fe_types)

    @property
    def target_data(self):
        raise NotImplementedError

    @cached_property
    def nonconst_module_attrs(self):
        """
        All module attrs are constant for targets using BaseContext.
        """
        return tuple()

    @cached_property
    def nrt(self):
        from numba.core.runtime.context import NRTContext
        return NRTContext(self, self.enable_nrt)

    def subtarget(self, **kws):
        obj = copy.copy(self)  # shallow copy
        for k, v in kws.items():
            if not hasattr(obj, k):
                raise NameError("unknown option {0!r}".format(k))
            setattr(obj, k, v)
        if obj.codegen() is not self.codegen():
            # We can't share functions across different codegens
            obj.cached_internal_func = {}
        return obj

    def install_registry(self, registry):
        """
        Install a *registry* (a imputils.Registry instance) of function
        and attribute implementations.
        """
        try:
            loader = self._registries[registry]
        except KeyError:
            loader = RegistryLoader(registry)
            self._registries[registry] = loader
        self.insert_func_defn(loader.new_registrations('functions'))
        self._insert_getattr_defn(loader.new_registrations('getattrs'))
        self._insert_setattr_defn(loader.new_registrations('setattrs'))
        self._insert_cast_defn(loader.new_registrations('casts'))
        self._insert_get_constant_defn(loader.new_registrations('constants'))

    def insert_func_defn(self, defns):
        for impl, func, sig in defns:
            self._defns[func].append(impl, sig)

    def _insert_getattr_defn(self, defns):
        for impl, attr, sig in defns:
            self._getattrs[attr].append(impl, sig)

    def _insert_setattr_defn(self, defns):
        for impl, attr, sig in defns:
            self._setattrs[attr].append(impl, sig)

    def _insert_cast_defn(self, defns):
        for impl, sig in defns:
            self._casts.append(impl, sig)

    def _insert_get_constant_defn(self, defns):
        for impl, sig in defns:
            self._get_constants.append(impl, sig)

    def insert_user_function(self, func, fndesc, libs=()):
        impl = user_function(fndesc, libs)
        self._defns[func].append(impl, impl.signature)

    def insert_generator(self, genty, gendesc, libs=()):
        assert isinstance(genty, types.Generator)
        impl = user_generator(gendesc, libs)
        self._generators[genty] = gendesc, impl

    def remove_user_function(self, func):
        """
        Remove user function *func*.
        KeyError is raised if the function isn't known to us.
        """
        del self._defns[func]

    def get_external_function_type(self, fndesc):
        argtypes = [self.get_argument_type(aty)
                    for aty in fndesc.argtypes]
        # don't wrap in pointer
        restype = self.get_argument_type(fndesc.restype)
        fnty = llvmir.FunctionType(restype, argtypes)
        return fnty

    def declare_function(self, module, fndesc):
        fnty = self.call_conv.get_function_type(fndesc.restype, fndesc.argtypes)
        fn = cgutils.get_or_insert_function(module, fnty, fndesc.mangled_name)
        self.call_conv.decorate_function(fn, fndesc.args, fndesc.argtypes, noalias=fndesc.noalias)
        if fndesc.inline:
            fn.attributes.add('alwaysinline')
            # alwaysinline overrides optnone
            fn.attributes.discard('noinline')
            fn.attributes.discard('optnone')
        return fn

    def declare_external_function(self, module, fndesc):
        fnty = self.get_external_function_type(fndesc)
        fn = cgutils.get_or_insert_function(module, fnty, fndesc.mangled_name)
        assert fn.is_declaration
        for ak, av in zip(fndesc.args, fn.args):
            av.name = "arg.%s" % ak
        return fn

    def insert_const_string(self, mod, string):
        """
        Insert constant *string* (a str object) into module *mod*.
        """
        stringtype = GENERIC_POINTER
        name = ".const.%s" % string
        text = cgutils.make_bytearray(string.encode("utf-8") + b"\x00")
        gv = self.insert_unique_const(mod, name, text)
        return Constant.bitcast(gv, stringtype)

    def insert_const_bytes(self, mod, bytes, name=None):
        """
        Insert constant *byte* (a `bytes` object) into module *mod*.
        """
        stringtype = GENERIC_POINTER
        name = ".bytes.%s" % (name or hash(bytes))
        text = cgutils.make_bytearray(bytes)
        gv = self.insert_unique_const(mod, name, text)
        return Constant.bitcast(gv, stringtype)

    def insert_unique_const(self, mod, name, val):
        """
        Insert a unique internal constant named *name*, with LLVM value
        *val*, into module *mod*.
        """
        try:
            gv = mod.get_global(name)
        except KeyError:
            return cgutils.global_constant(mod, name, val)
        else:
            return gv

    def get_argument_type(self, ty):
        return self.data_model_manager[ty].get_argument_type()

    def get_return_type(self, ty):
        return self.data_model_manager[ty].get_return_type()

    def get_data_type(self, ty):
        """
        Get a LLVM data representation of the Numba type *ty* that is safe
        for storage.  Record data are stored as byte array.

        The return value is a llvmlite.ir.Type object, or None if the type
        is an opaque pointer (???).
        """
        return self.data_model_manager[ty].get_data_type()

    def get_value_type(self, ty):
        return self.data_model_manager[ty].get_value_type()

    def pack_value(self, builder, ty, value, ptr, align=None):
        """
        Pack value into the array storage at *ptr*.
        If *align* is given, it is the guaranteed alignment for *ptr*
        (by default, the standard ABI alignment).
        """
        dataval = self.data_model_manager[ty].as_data(builder, value)
        builder.store(dataval, ptr, align=align)

    def unpack_value(self, builder, ty, ptr, align=None):
        """
        Unpack value from the array storage at *ptr*.
        If *align* is given, it is the guaranteed alignment for *ptr*
        (by default, the standard ABI alignment).
        """
        dm = self.data_model_manager[ty]
        return dm.load_from_data_pointer(builder, ptr, align)

    def get_constant_generic(self, builder, ty, val):
        """
        Return a LLVM constant representing value *val* of Numba type *ty*.
        """
        try:
            impl = self._get_constants.find((ty,))
            return impl(self, builder, ty, val)
        except NotImplementedError:
            raise NotImplementedError("Cannot lower constant of type '%s'" % (ty,))

    def get_constant(self, ty, val):
        """
        Same as get_constant_generic(), but without specifying *builder*.
        Works only for simple types.
        """
        # HACK: pass builder=None to preserve get_constant() API
        return self.get_constant_generic(None, ty, val)

    def get_constant_undef(self, ty):
        lty = self.get_value_type(ty)
        return Constant(lty, llvmir.Undefined)

    def get_constant_null(self, ty):
        lty = self.get_value_type(ty)
        return Constant(lty, None)

    def get_function(self, fn, sig, _firstcall=True):
        """
        Return the implementation of function *fn* for signature *sig*.
        The return value is a callable with the signature (builder, args).
        """
        assert sig is not None
        sig = sig.as_function()
        if isinstance(fn, types.Callable):
            key = fn.get_impl_key(sig)
            overloads = self._defns[key]
        else:
            key = fn
            overloads = self._defns[key]

        try:
            return _wrap_impl(overloads.find(sig.args), self, sig)
        except errors.NumbaNotImplementedError:
            pass
        if isinstance(fn, types.Type):
            # It's a type instance => try to find a definition for the type class
            try:
                return self.get_function(type(fn), sig)
            except NotImplementedError:
                # Raise exception for the type instance, for a better error message
                pass

        # Automatically refresh the context to load new registries if we are
        # calling the first time.
        if _firstcall:
            self.refresh()
            return self.get_function(fn, sig, _firstcall=False)

        raise NotImplementedError("No definition for lowering %s%s" % (key, sig))

    def get_generator_desc(self, genty):
        """
        """
        return self._generators[genty][0]

    def get_generator_impl(self, genty):
        """
        """
        res = self._generators[genty][1]
        self.add_linking_libs(getattr(res, 'libs', ()))
        return res

    def get_bound_function(self, builder, obj, ty):
        assert self.get_value_type(ty) == obj.type
        return obj

    def get_getattr(self, typ, attr):
        """
        Get the getattr() implementation for the given type and attribute name.
        The return value is a callable with the signature
        (context, builder, typ, val, attr).
        """
        const_attr = (typ, attr) not in self.nonconst_module_attrs
        is_module = isinstance(typ, types.Module)
        if is_module and const_attr:
            # Implement getattr for module-level globals that we treat as
            # constants.
            # XXX We shouldn't have to retype this
            attrty = self.typing_context.resolve_module_constants(typ, attr)
            if attrty is None or isinstance(attrty, types.Dummy):
                # No implementation required for dummies (functions, modules...),
                # which are dealt with later
                return None
            else:
                pyval = getattr(typ.pymod, attr)
                def imp(context, builder, typ, val, attr):
                    llval = self.get_constant_generic(builder, attrty, pyval)
                    return impl_ret_borrowed(context, builder, attrty, llval)
                return imp

        # Lookup specific getattr implementation for this type and attribute
        overloads = self._getattrs[attr]
        try:
            return overloads.find((typ,))
        except errors.NumbaNotImplementedError:
            pass
        # Lookup generic getattr implementation for this type
        overloads = self._getattrs[None]
        try:
            return overloads.find((typ,))
        except errors.NumbaNotImplementedError:
            pass

        raise NotImplementedError("No definition for lowering %s.%s" % (typ, attr))

    def get_setattr(self, attr, sig):
        """
        Get the setattr() implementation for the given attribute name
        and signature.
        The return value is a callable with the signature (builder, args).
        """
        assert len(sig.args) == 2
        typ = sig.args[0]
        valty = sig.args[1]

        def wrap_setattr(impl):
            def wrapped(builder, args):
                return impl(self, builder, sig, args, attr)
            return wrapped

        # Lookup specific setattr implementation for this type and attribute
        overloads = self._setattrs[attr]
        try:
            return wrap_setattr(overloads.find((typ, valty)))
        except errors.NumbaNotImplementedError:
            pass
        # Lookup generic setattr implementation for this type
        overloads = self._setattrs[None]
        try:
            return wrap_setattr(overloads.find((typ, valty)))
        except errors.NumbaNotImplementedError:
            pass

        raise NotImplementedError("No definition for lowering %s.%s = %s"
                                  % (typ, attr, valty))

    def get_argument_value(self, builder, ty, val):
        """
        Argument representation to local value representation
        """
        return self.data_model_manager[ty].from_argument(builder, val)

    def get_returned_value(self, builder, ty, val):
        """
        Return value representation to local value representation
        """
        return self.data_model_manager[ty].from_return(builder, val)

    def get_return_value(self, builder, ty, val):
        """
        Local value representation to return type representation
        """
        return self.data_model_manager[ty].as_return(builder, val)

    def get_value_as_argument(self, builder, ty, val):
        """Prepare local value representation as argument type representation
        """
        return self.data_model_manager[ty].as_argument(builder, val)

    def get_value_as_data(self, builder, ty, val):
        return self.data_model_manager[ty].as_data(builder, val)

    def get_data_as_value(self, builder, ty, val):
        return self.data_model_manager[ty].from_data(builder, val)

    def pair_first(self, builder, val, ty):
        """
        Extract the first element of a heterogeneous pair.
        """
        pair = self.make_helper(builder, ty, val)
        return pair.first

    def pair_second(self, builder, val, ty):
        """
        Extract the second element of a heterogeneous pair.
        """
        pair = self.make_helper(builder, ty, val)
        return pair.second

    def cast(self, builder, val, fromty, toty):
        """
        Cast a value of type *fromty* to type *toty*.
        This implements implicit conversions as can happen due to the
        granularity of the Numba type system, or lax Python semantics.
        """
        if fromty is types._undef_var:
            # Special case for undefined variable
            return self.get_constant_null(toty)
        elif fromty == toty or toty == types.Any:
            return val
        try:
            impl = self._casts.find((fromty, toty))
            return impl(self, builder, fromty, toty, val)
        except errors.NumbaNotImplementedError:
            raise errors.NumbaNotImplementedError(
                "Cannot cast %s to %s: %s" % (fromty, toty, val))

    def generic_compare(self, builder, key, argtypes, args):
        """
        Compare the given LLVM values of the given Numba types using
        the comparison *key* (e.g. '==').  The values are first cast to
        a common safe conversion type.
        """
        at, bt = argtypes
        av, bv = args
        ty = self.typing_context.unify_types(at, bt)
        assert ty is not None
        cav = self.cast(builder, av, at, ty)
        cbv = self.cast(builder, bv, bt, ty)
        fnty = self.typing_context.resolve_value_type(key)
        # the sig is homogeneous in the unified casted type
        cmpsig = fnty.get_call_type(self.typing_context, (ty, ty), {})
        cmpfunc = self.get_function(fnty, cmpsig)
        self.add_linking_libs(getattr(cmpfunc, 'libs', ()))
        return cmpfunc(builder, (cav, cbv))

    def make_optional_none(self, builder, valtype):
        optval = self.make_helper(builder, types.Optional(valtype))
        optval.valid = cgutils.false_bit
        return optval._getvalue()

    def make_optional_value(self, builder, valtype, value):
        optval = self.make_helper(builder, types.Optional(valtype))
        optval.valid = cgutils.true_bit
        optval.data = value
        return optval._getvalue()

    def is_true(self, builder, typ, val):
        """
        Return the truth value of a value of the given Numba type.
        """
        fnty = self.typing_context.resolve_value_type(bool)
        sig = fnty.get_call_type(self.typing_context, (typ,), {})
        impl = self.get_function(fnty, sig)
        return impl(builder, (val,))

    def get_c_value(self, builder, typ, name, dllimport=False):
        """
        Get a global value through its C-accessible *name*, with the given
        LLVM type.
        If *dllimport* is true, the symbol will be marked as imported
        from a DLL (necessary for AOT compilation under Windows).
        """
        module = builder.function.module
        try:
            gv = module.globals[name]
        except KeyError:
            gv = cgutils.add_global_variable(module, typ, name)
            if dllimport and self.aot_mode and sys.platform == 'win32':
                gv.storage_class = "dllimport"
        return gv

    def call_external_function(self, builder, callee, argtys, args):
        args = [self.get_value_as_argument(builder, ty, arg)
                for ty, arg in zip(argtys, args)]
        retval = builder.call(callee, args)
        return retval

    def get_function_pointer_type(self, typ):
        return self.data_model_manager[typ].get_data_type()

    def call_function_pointer(self, builder, funcptr, args, cconv=None):
        return builder.call(funcptr, args, cconv=cconv)

    def print_string(self, builder, text):
        mod = builder.module
        cstring = GENERIC_POINTER
        fnty = llvmir.FunctionType(llvmir.IntType(32), [cstring])
        puts = cgutils.get_or_insert_function(mod, fnty, "puts")
        return builder.call(puts, [text])

    def debug_print(self, builder, text):
        mod = builder.module
        cstr = self.insert_const_string(mod, str(text))
        self.print_string(builder, cstr)

    def printf(self, builder, format_string, *args):
        mod = builder.module
        if isinstance(format_string, str):
            cstr = self.insert_const_string(mod, format_string)
        else:
            cstr = format_string
        fnty = llvmir.FunctionType(llvmir.IntType(32), (GENERIC_POINTER,), var_arg=True)
        fn = cgutils.get_or_insert_function(mod, fnty, "printf")
        return builder.call(fn, (cstr,) + tuple(args))

    def get_struct_type(self, struct):
        """
        Get the LLVM struct type for the given Structure class *struct*.
        """
        fields = [self.get_value_type(v) for _, v in struct._fields]
        return llvmir.LiteralStructType(fields)

    def get_dummy_value(self):
        return Constant(self.get_dummy_type(), None)

    def get_dummy_type(self):
        return GENERIC_POINTER

    def _compile_subroutine_no_cache(self, builder, impl, sig, locals={},
                                     flags=None):
        """
        Invoke the compiler to compile a function to be used inside a
        nopython function, but without generating code to call that
        function.

        Note this context's flags are not inherited.
        """
        # Compile
        from numba.core import compiler

        with global_compiler_lock:
            codegen = self.codegen()
            library = codegen.create_library(impl.__name__)
            if flags is None:

                cstk = targetconfig.ConfigStack()
                flags = compiler.Flags()
                if cstk:
                    tls_flags = cstk.top()
                    if tls_flags.is_set("nrt") and tls_flags.nrt:
                        flags.nrt = True

            flags.no_compile = True
            flags.no_cpython_wrapper = True
            flags.no_cfunc_wrapper = True

            cres = compiler.compile_internal(self.typing_context, self,
                                             library,
                                             impl, sig.args,
                                             sig.return_type, flags,
                                             locals=locals)

            # Allow inlining the function inside callers.
            self.active_code_library.add_linking_library(cres.library)
            return cres

    def compile_subroutine(self, builder, impl, sig, locals={}, flags=None,
                           caching=True):
        """
        Compile the function *impl* for the given *sig* (in nopython mode).
        Return an instance of CompileResult.

        If *caching* evaluates True, the function keeps the compiled function
        for reuse in *.cached_internal_func*.
        """
        cache_key = (impl.__code__, sig, type(self.error_model))
        if not caching:
            cached = None
        else:
            if impl.__closure__:
                # XXX This obviously won't work if a cell's value is
                # unhashable.
                cache_key += tuple(c.cell_contents for c in impl.__closure__)
            cached = self.cached_internal_func.get(cache_key)
        if cached is None:
            cres = self._compile_subroutine_no_cache(builder, impl, sig,
                                                     locals=locals,
                                                     flags=flags)
            self.cached_internal_func[cache_key] = cres

        cres = self.cached_internal_func[cache_key]
        # Allow inlining the function inside callers.
        self.active_code_library.add_linking_library(cres.library)
        return cres

    def compile_internal(self, builder, impl, sig, args, locals={}):
        """
        Like compile_subroutine(), but also call the function with the given
        *args*.
        """
        cres = self.compile_subroutine(builder, impl, sig, locals)
        return self.call_internal(builder, cres.fndesc, sig, args)

    def call_internal(self, builder, fndesc, sig, args):
        """
        Given the function descriptor of an internally compiled function,
        emit a call to that function with the given arguments.
        """
        status, res = self.call_internal_no_propagate(builder, fndesc, sig, args)
        with cgutils.if_unlikely(builder, status.is_error):
            self.call_conv.return_status_propagate(builder, status)

        res = imputils.fix_returning_optional(self, builder, sig, status, res)
        return res

    def call_internal_no_propagate(self, builder, fndesc, sig, args):
        """Similar to `.call_internal()` but does not handle or propagate
        the return status automatically.
        """
        # Add call to the generated function
        llvm_mod = builder.module
        fn = self.declare_function(llvm_mod, fndesc)
        status, res = self.call_conv.call_function(builder, fn, sig.return_type,
                                                   sig.args, args)
        return status, res

    def call_unresolved(self, builder, name, sig, args):
        """
        Insert a function call to an unresolved symbol with the given *name*.

        Note: this is used for recursive call.

        In the mutual recursion case::

            @njit
            def foo():
                ...  # calls bar()

            @njit
            def bar():
                ... # calls foo()

            foo()

        When foo() is called, the compilation of bar() is fully completed
        (codegen'ed and loaded) before foo() is. Since MCJIT's eager compilation
        doesn't allow loading modules with declare-only functions (which is
        needed for foo() in bar()), the call_unresolved injects a global
        variable that the "linker" can update even after the module is loaded by
        MCJIT. The linker would allocate space for the global variable before
        the bar() module is loaded. When later foo() module is defined, it will
        update bar()'s reference to foo().

        The legacy lazy JIT and the new ORC JIT would allow a declare-only
        function be used in a module as long as it is defined by the time of its
        first use.
        """
        # Insert an unresolved reference to the function being called.
        codegen = self.codegen()
        fnty = self.call_conv.get_function_type(sig.return_type, sig.args)
        fn = codegen.insert_unresolved_ref(builder, fnty, name)
        # Normal call sequence
        status, res = self.call_conv.call_function(builder, fn, sig.return_type,
                                                   sig.args, args)
        with cgutils.if_unlikely(builder, status.is_error):
            self.call_conv.return_status_propagate(builder, status)

        res = imputils.fix_returning_optional(self, builder, sig, status, res)
        return res

    def get_executable(self, func, fndesc, env):
        raise NotImplementedError

    def get_python_api(self, builder):
        return PythonAPI(self, builder)

    def sentry_record_alignment(self, rectyp, attr):
        """
        Assumes offset starts from a properly aligned location
        """
        if self.strict_alignment:
            offset = rectyp.offset(attr)
            elemty = rectyp.typeof(attr)
            if isinstance(elemty, types.NestedArray):
                # For a NestedArray we need to consider the data type of
                # elements of the array for alignment, not the array structure
                # itself
                elemty = elemty.dtype
            align = self.get_abi_alignment(self.get_data_type(elemty))
            if offset % align:
                msg = "{rec}.{attr} of type {type} is not aligned".format(
                    rec=rectyp, attr=attr, type=elemty)
                raise TypeError(msg)

    def get_helper_class(self, typ, kind='value'):
        """
        Get a helper class for the given *typ*.
        """
        # XXX handle all types: complex, array, etc.
        # XXX should it be a method on the model instead? this would allow a default kind...
        return cgutils.create_struct_proxy(typ, kind)

    def _make_helper(self, builder, typ, value=None, ref=None, kind='value'):
        cls = self.get_helper_class(typ, kind)
        return cls(self, builder, value=value, ref=ref)

    def make_helper(self, builder, typ, value=None, ref=None):
        """
        Get a helper object to access the *typ*'s members,
        for the given value or reference.
        """
        return self._make_helper(builder, typ, value, ref, kind='value')

    def make_data_helper(self, builder, typ, ref=None):
        """
        As make_helper(), but considers the value as stored in memory,
        rather than a live value.
        """
        return self._make_helper(builder, typ, ref=ref, kind='data')

    def make_array(self, typ):
        from numba.np import arrayobj
        return arrayobj.make_array(typ)

    def populate_array(self, arr, **kwargs):
        """
        Populate array structure.
        """
        from numba.np import arrayobj
        return arrayobj.populate_array(arr, **kwargs)

    def make_complex(self, builder, typ, value=None):
        """
        Get a helper object to access the given complex numbers' members.
        """
        assert isinstance(typ, types.Complex), typ
        return self.make_helper(builder, typ, value)

    def make_tuple(self, builder, typ, values):
        """
        Create a tuple of the given *typ* containing the *values*.
        """
        tup = self.get_constant_undef(typ)
        for i, val in enumerate(values):
            tup = builder.insert_value(tup, val, i)
        return tup

    def make_constant_array(self, builder, typ, ary):
        """
        Create an array structure reifying the given constant array.
        A low-level contiguous array constant is created in the LLVM IR.
        """
        datatype = self.get_data_type(typ.dtype)
        # don't freeze ary of non-contig or bigger than 1MB
        size_limit = 10**6

        if (self.allow_dynamic_globals and
                (typ.layout not in 'FC' or ary.nbytes > size_limit)):
            # get pointer from the ary
            dataptr = ary.ctypes.data
            data = self.add_dynamic_addr(builder, dataptr, info=str(type(dataptr)))
            rt_addr = self.add_dynamic_addr(builder, id(ary), info=str(type(ary)))
        else:
            # Handle data: reify the flattened array in "C" or "F" order as a
            # global array of bytes.
            flat = ary.flatten(order=typ.layout)
            # Note: we use `bytearray(flat.data)` instead of `bytearray(flat)` to
            #       workaround issue #1850 which is due to numpy issue #3147
            consts = cgutils.create_constant_array(llvmir.IntType(8), bytearray(flat.data))
            data = cgutils.global_constant(builder, ".const.array.data", consts)
            # Ensure correct data alignment (issue #1933)
            data.align = self.get_abi_alignment(datatype)
            # No reference to parent ndarray
            rt_addr = None

        # Handle shape
        llintp = self.get_value_type(types.intp)
        shapevals = [self.get_constant(types.intp, s) for s in ary.shape]
        cshape = cgutils.create_constant_array(llintp, shapevals)

        # Handle strides
        stridevals = [self.get_constant(types.intp, s) for s in ary.strides]
        cstrides = cgutils.create_constant_array(llintp, stridevals)

        # Create array structure
        cary = self.make_array(typ)(self, builder)

        intp_itemsize = self.get_constant(types.intp, ary.dtype.itemsize)
        self.populate_array(cary,
                            data=builder.bitcast(data, cary.data.type),
                            shape=cshape,
                            strides=cstrides,
                            itemsize=intp_itemsize,
                            parent=rt_addr,
                            meminfo=None)

        return cary._getvalue()

    def add_dynamic_addr(self, builder, intaddr, info):
        """
        Returns dynamic address as a void pointer `i8*`.

        Internally, a global variable is added to inform the lowerer about
        the usage of dynamic addresses.  Caching will be disabled.
        """
        assert self.allow_dynamic_globals, "dyn globals disabled in this target"
        assert isinstance(intaddr, int), 'dyn addr not of int type'
        mod = builder.module
        llvoidptr = self.get_value_type(types.voidptr)
        addr = self.get_constant(types.uintp, intaddr).inttoptr(llvoidptr)
        # Use a unique name by embedding the address value
        symname = 'numba.dynamic.globals.{:x}'.format(intaddr)
        gv = cgutils.add_global_variable(mod, llvoidptr, symname)
        # Use linkonce linkage to allow merging with other GV of the same name.
        # And, avoid optimization from assuming its value.
        gv.linkage = 'linkonce'
        gv.initializer = addr
        return builder.load(gv)

    def get_abi_sizeof(self, ty):
        """
        Get the ABI size of LLVM type *ty*.
        """
        assert isinstance(ty, llvmir.Type), "Expected LLVM type"
        return ty.get_abi_size(self.target_data)

    def get_abi_alignment(self, ty):
        """
        Get the ABI alignment of LLVM type *ty*.
        """
        assert isinstance(ty, llvmir.Type), "Expected LLVM type"
        return ty.get_abi_alignment(self.target_data)

    def get_preferred_array_alignment(context, ty):
        """
        Get preferred array alignment for Numba type *ty*.
        """
        # AVX prefers 32-byte alignment
        return 32

    def post_lowering(self, mod, library):
        """Run target specific post-lowering transformation here.
        """

    def create_module(self, name):
        """Create a LLVM module

        The default implementation in BaseContext always raises a
        ``NotImplementedError`` exception. Subclasses should implement
        this method.
        """
        raise NotImplementedError

    @property
    def active_code_library(self):
        """Get the active code library
        """
        return self._codelib_stack[-1]

    @contextmanager
    def push_code_library(self, lib):
        """Push the active code library for the context
        """
        self._codelib_stack.append(lib)
        try:
            yield
        finally:
            self._codelib_stack.pop()

    def add_linking_libs(self, libs):
        """Add iterable of linking libraries to the *active_code_library*.
        """
        colib = self.active_code_library
        for lib in libs:
            colib.add_linking_library(lib)

    def get_ufunc_info(self, ufunc_key):
        """Get the ufunc implementation for a given ufunc object.

        The default implementation in BaseContext always raises a
        ``NotImplementedError`` exception. Subclasses may raise ``KeyError``
        to signal that the given ``ufunc_key`` is not available.

        Parameters
        ----------
        ufunc_key : NumPy ufunc

        Returns
        -------
        res : dict[str, callable]
            A mapping of a NumPy ufunc type signature to a lower-level
            implementation.
        """
        raise NotImplementedError(f"{self} does not support ufunc")

class _wrap_impl(object):
    """
    A wrapper object to call an implementation function with some predefined
    (context, signature) arguments.
    The wrapper also forwards attribute queries, which is important.
    """

    def __init__(self, imp, context, sig):
        self._callable = _wrap_missing_loc(imp)
        self._imp = self._callable()
        self._context = context
        self._sig = sig

    def __call__(self, builder, args, loc=None):
        res = self._imp(self._context, builder, self._sig, args, loc=loc)
        self._context.add_linking_libs(getattr(self, 'libs', ()))
        return res

    def __getattr__(self, item):
        return getattr(self._imp, item)

    def __repr__(self):
        return "<wrapped %s>" % repr(self._callable)

def _has_loc(fn):
    """Does function *fn* take ``loc`` argument?
    """
    sig = utils.pysignature(fn)
    return 'loc' in sig.parameters


class _wrap_missing_loc(object):

    def __init__(self, fn):
        self.func = fn # store this to help with debug

    def __call__(self):
        """Wrap function for missing ``loc`` keyword argument.
        Otherwise, return the original *fn*.
        """
        fn = self.func
        if not _has_loc(fn):
            def wrapper(*args, **kwargs):
                kwargs.pop('loc')     # drop unused loc
                return fn(*args, **kwargs)

            # Copy the following attributes from the wrapped.
            # Following similar implementation as functools.wraps but
            # ignore attributes if not available (i.e fix py2.7)
            attrs = '__name__', 'libs'
            for attr in attrs:
                try:
                    val = getattr(fn, attr)
                except AttributeError:
                    pass
                else:
                    setattr(wrapper, attr, val)

            return wrapper
        else:
            return fn

    def __repr__(self):
        return "<wrapped %s>" % self.func


@utils.runonce
def _initialize_llvm_lock_event():
    """Initial event triggers for LLVM lock
    """
    def enter_fn():
        event.start_event("numba:llvm_lock")

    def exit_fn():
        event.end_event("numba:llvm_lock")

    ll.ffi.register_lock_callback(enter_fn, exit_fn)


_initialize_llvm_lock_event()
