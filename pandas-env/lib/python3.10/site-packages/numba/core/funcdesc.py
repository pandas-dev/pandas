"""
Function descriptors.
"""

from collections import defaultdict
import importlib

from numba.core import types, itanium_mangler
from numba.core.utils import _dynamic_modname, _dynamic_module


def default_mangler(name, argtypes, *, abi_tags=(), uid=None):
    return itanium_mangler.mangle(name, argtypes, abi_tags=abi_tags, uid=uid)


def qualifying_prefix(modname, qualname):
    """
    Returns a new string that is used for the first half of the mangled name.
    """
    # XXX choose a different convention for object mode
    return '{}.{}'.format(modname, qualname) if modname else qualname


class FunctionDescriptor(object):
    """
    Base class for function descriptors: an object used to carry
    useful metadata about a natively callable function.

    Note that while `FunctionIdentity` denotes a Python function
    which is being concretely compiled by Numba, `FunctionDescriptor`
    may be more "abstract".
    """
    __slots__ = ('native', 'modname', 'qualname', 'doc', 'typemap',
                 'calltypes', 'args', 'kws', 'restype', 'argtypes',
                 'mangled_name', 'unique_name', 'env_name', 'global_dict',
                 'inline', 'noalias', 'abi_tags', 'uid')

    def __init__(self, native, modname, qualname, unique_name, doc,
                 typemap, restype, calltypes, args, kws, mangler=None,
                 argtypes=None, inline=False, noalias=False, env_name=None,
                 global_dict=None, abi_tags=(), uid=None):
        self.native = native
        self.modname = modname
        self.global_dict = global_dict
        self.qualname = qualname
        self.unique_name = unique_name
        self.doc = doc
        # XXX typemap and calltypes should be on the compile result,
        # not the FunctionDescriptor
        self.typemap = typemap
        self.calltypes = calltypes
        self.args = args
        self.kws = kws
        self.restype = restype
        # Argument types
        if argtypes is not None:
            assert isinstance(argtypes, tuple), argtypes
            self.argtypes = argtypes
        else:
            # Get argument types from the type inference result
            # (note the "arg.FOO" convention as used in typeinfer
            self.argtypes = tuple(self.typemap['arg.' + a] for a in args)
        mangler = default_mangler if mangler is None else mangler
        # The mangled name *must* be unique, else the wrong function can
        # be chosen at link time.
        qualprefix = qualifying_prefix(self.modname, self.qualname)
        self.uid = uid
        self.mangled_name = mangler(
            qualprefix, self.argtypes, abi_tags=abi_tags, uid=uid,
        )
        if env_name is None:
            env_name = mangler(".NumbaEnv.{}".format(qualprefix),
                               self.argtypes, abi_tags=abi_tags, uid=uid)
        self.env_name = env_name
        self.inline = inline
        self.noalias = noalias
        self.abi_tags = abi_tags

    def lookup_globals(self):
        """
        Return the global dictionary of the function.
        It may not match the Module's globals if the function is created
        dynamically (i.e. exec)
        """
        return self.global_dict or self.lookup_module().__dict__

    def lookup_module(self):
        """
        Return the module in which this function is supposed to exist.
        This may be a dummy module if the function was dynamically
        generated or the module can't be found.
        """
        if self.modname == _dynamic_modname:
            return _dynamic_module
        else:
            try:
                # ensure module exist
                return importlib.import_module(self.modname)
            except ImportError:
                return _dynamic_module

    def lookup_function(self):
        """
        Return the original function object described by this object.
        """
        return getattr(self.lookup_module(), self.qualname)

    @property
    def llvm_func_name(self):
        """
        The LLVM-registered name for the raw function.
        """
        return self.mangled_name

    # XXX refactor this

    @property
    def llvm_cpython_wrapper_name(self):
        """
        The LLVM-registered name for a CPython-compatible wrapper of the
        raw function (i.e. a PyCFunctionWithKeywords).
        """
        return itanium_mangler.prepend_namespace(self.mangled_name,
                                                 ns='cpython')

    @property
    def llvm_cfunc_wrapper_name(self):
        """
        The LLVM-registered name for a C-compatible wrapper of the
        raw function.
        """
        return 'cfunc.' + self.mangled_name

    def __repr__(self):
        return "<function descriptor %r>" % (self.unique_name)

    @classmethod
    def _get_function_info(cls, func_ir):
        """
        Returns
        -------
        qualname, unique_name, modname, doc, args, kws, globals

        ``unique_name`` must be a unique name.
        """
        func = func_ir.func_id.func
        qualname = func_ir.func_id.func_qualname
        # XXX to func_id
        modname = func.__module__
        doc = func.__doc__ or ''
        args = tuple(func_ir.arg_names)
        kws = ()        # TODO
        global_dict = None

        if modname is None:
            # Dynamically generated function.
            modname = _dynamic_modname
            # Retain a reference to the dictionary of the function.
            # This disables caching, serialization and pickling.
            global_dict = func_ir.func_id.func.__globals__

        unique_name = func_ir.func_id.unique_name

        return qualname, unique_name, modname, doc, args, kws, global_dict

    @classmethod
    def _from_python_function(cls, func_ir, typemap, restype,
                              calltypes, native, mangler=None,
                              inline=False, noalias=False, abi_tags=()):
        (qualname, unique_name, modname, doc, args, kws, global_dict,
         ) = cls._get_function_info(func_ir)

        self = cls(native, modname, qualname, unique_name, doc,
                   typemap, restype, calltypes,
                   args, kws, mangler=mangler, inline=inline, noalias=noalias,
                   global_dict=global_dict, abi_tags=abi_tags,
                   uid=func_ir.func_id.unique_id)
        return self


class PythonFunctionDescriptor(FunctionDescriptor):
    """
    A FunctionDescriptor subclass for Numba-compiled functions.
    """
    __slots__ = ()

    @classmethod
    def from_specialized_function(cls, func_ir, typemap, restype, calltypes,
                                  mangler, inline, noalias, abi_tags):
        """
        Build a FunctionDescriptor for a given specialization of a Python
        function (in nopython mode).
        """
        return cls._from_python_function(func_ir, typemap, restype, calltypes,
                                         native=True, mangler=mangler,
                                         inline=inline, noalias=noalias,
                                         abi_tags=abi_tags)

    @classmethod
    def from_object_mode_function(cls, func_ir):
        """
        Build a FunctionDescriptor for an object mode variant of a Python
        function.
        """
        typemap = defaultdict(lambda: types.pyobject)
        calltypes = typemap.copy()
        restype = types.pyobject
        return cls._from_python_function(func_ir, typemap, restype, calltypes,
                                         native=False)


class ExternalFunctionDescriptor(FunctionDescriptor):
    """
    A FunctionDescriptor subclass for opaque external functions
    (e.g. raw C functions).
    """
    __slots__ = ()

    def __init__(self, name, restype, argtypes):
        args = ["arg%d" % i for i in range(len(argtypes))]

        def mangler(a, x, abi_tags, uid=None):
            return a
        super(ExternalFunctionDescriptor, self
              ).__init__(native=True, modname=None, qualname=name,
                         unique_name=name, doc='', typemap=None,
                         restype=restype, calltypes=None, args=args,
                         kws=None,
                         mangler=mangler,
                         argtypes=argtypes)
