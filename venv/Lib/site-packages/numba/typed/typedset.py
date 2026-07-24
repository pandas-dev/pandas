"""
Python wrapper that connects CPython interpreter to the numba setobject.
"""
from collections.abc import MutableSet
from numba.core.types import SetType
from numba.core.imputils import numba_typeref_ctor
from numba import njit, typeof
from numba.core import types, errors, config, cgutils
from numba.core.extending import (
    overload,
    box,
    unbox,
    NativeValue,
    type_callable,
    overload_classmethod,
)
from numba.typed import setobject
from numba.core.typing import signature


@njit
def _make_set(keyty):
    return setobject._as_meminfo(setobject.new_set(keyty))


@njit
def _length(s):
    return len(s)


@njit
def _copy(s):
    return s.copy()


@njit
def _additem(s, key):
    s.add(key)


@njit
def _set_contains(s, key):
    return (key in s)


@njit
def _discarditem(s, key):
    return s.discard(key)


@njit
def _iter(s):
    return list(s)


def _from_meminfo_ptr(ptr, settype):
    return Set(meminfo=ptr, settype=settype)


class Set(MutableSet):
    """A typed-set usable in Numba compiled functions.

    Implements the MutableSet interface.
    """

    def __new__(cls, settype=None, meminfo=None):
        if config.DISABLE_JIT:
            return set.__new__(set)
        else:
            return object.__new__(cls)

    @classmethod
    def empty(cls, key_type):
        """Create a new empty Set with *key_type*
        as the types for the values of the set.
        """
        if config.DISABLE_JIT:
            return set()
        else:
            return cls(settype=SetType(key_type))

    def __init__(self, **kwargs):
        """
        For users, the constructor does not take any parameters.
        The keyword arguments are for internal use only.

        Parameters
        ----------
        settype : numba.core.types.SetType; keyword-only
            Used internally for the set type.
        meminfo : MemInfo; keyword-only
            Used internally to pass the MemInfo object when boxing.
        """
        if kwargs:
            self._set_type, self._opaque = self._parse_arg(**kwargs)
        else:
            self._set_type = None

    def _parse_arg(self, settype, meminfo=None):
        if not isinstance(settype, SetType):
            raise TypeError('*settype* must be a SetType')

        if meminfo is not None:
            opaque = meminfo
        else:
            opaque = _make_set(settype.key_type)
        return settype, opaque

    @property
    def _numba_type_(self):
        if self._set_type is None:
            raise TypeError("invalid operation on untyped set")
        return self._set_type

    @property
    def _typed(self):
        """Returns True if the set is typed.
        """
        return self._set_type is not None

    def _initialise_set(self, key):
        settype = types.SetType(typeof(key))
        self._set_type, self._opaque = self._parse_arg(settype)

    def __len__(self):
        if not self._typed:
            return 0
        else:
            return _length(self)

    def copy(self):
        if not self._typed:
            raise TypeError("Cannot copy an untyped set")

        return _copy(self)

    def __contains__(self, key):
        if len(self) == 0:
            return False

        return bool(_set_contains(self, key))

    def __iter__(self):
        if not self._typed:
            return iter(())

        return iter(_iter(self))

    def add(self, key):
        if not self._typed:
            self._initialise_set(key)

        return _additem(self, key)

    def discard(self, key):
        if not self._typed:
            raise TypeError("Cannot discard from an untyped set")

        _discarditem(self, key)

    def __str__(self):
        buf = []
        for k in self:
            buf.append("{}".format(k))
        return '{{{0}}}'.format(', '.join(buf))

    def __repr__(self):
        body = str(self)
        prefix = str(self._set_type) if self._set_type else "UntypedSet"
        return "{prefix}({body})".format(prefix=prefix, body=body)


@overload_classmethod(types.SetType, 'empty')
def typedset_empty(cls, key_type):
    if cls.instance_type is not SetType:
        return

    def impl(cls, key_type):
        return setobject.new_set(key_type)

    return impl


@box(types.SetType)
def box_settype(typ, val, c):
    context = c.context
    builder = c.builder

    ctor = cgutils.create_struct_proxy(typ)
    dstruct = ctor(context, builder, value=val)
    # Returns the plain MemInfo
    boxed_meminfo = c.box(
        types.MemInfoPointer(types.voidptr),
        dstruct.meminfo,
    )

    modname = c.context.insert_const_string(
        c.builder.module, 'numba.typed.typedset',
    )
    typedset_mod = c.pyapi.import_module(modname)
    fmp_fn = c.pyapi.object_getattr_string(typedset_mod, '_from_meminfo_ptr')

    settype_obj = c.pyapi.unserialize(c.pyapi.serialize_object(typ))

    result_var = builder.alloca(c.pyapi.pyobj)
    builder.store(cgutils.get_null_value(c.pyapi.pyobj), result_var)
    with builder.if_then(cgutils.is_not_null(builder, settype_obj)):
        res = c.pyapi.call_function_objargs(
            fmp_fn, (boxed_meminfo, settype_obj),
        )
        c.pyapi.decref(fmp_fn)
        c.pyapi.decref(typedset_mod)
        c.pyapi.decref(boxed_meminfo)
        builder.store(res, result_var)
    return builder.load(result_var)


@unbox(types.SetType)
def unbox_settype(typ, val, c):
    context = c.context

    set_type = c.pyapi.unserialize(c.pyapi.serialize_object(Set))
    key_type = c.pyapi.object_type(val)
    same_type = c.builder.icmp_unsigned("==", key_type, set_type)

    with c.builder.if_else(same_type) as (then, orelse):
        with then:
            miptr = c.pyapi.object_getattr_string(val, '_opaque')

            mip_type = types.MemInfoPointer(types.voidptr)
            native = c.unbox(mip_type, miptr)

            mi = native.value

            argtypes = mip_type, typeof(typ)

            def convert(mi, typ):
                return setobject._from_meminfo(mi, typ)

            sig = signature(typ, *argtypes)
            nil_typeref = context.get_constant_null(argtypes[1])
            args = (mi, nil_typeref)
            is_error, setobj = c.pyapi.call_jit_code(convert , sig, args)
            # decref here because we are stealing a reference.
            c.context.nrt.decref(c.builder, typ, setobj)

            c.pyapi.decref(miptr)
            bb_unboxed = c.builder.basic_block

        with orelse:
            # Raise error on incorrect type
            c.pyapi.err_format(
                "PyExc_TypeError",
                "can't unbox a %S as a %S",
                key_type, set_type,
            )
            bb_else = c.builder.basic_block

    # Phi nodes to gather the output
    setobj_res = c.builder.phi(setobj.type)
    is_error_res = c.builder.phi(is_error.type)

    setobj_res.add_incoming(setobj, bb_unboxed)
    setobj_res.add_incoming(setobj.type(None), bb_else)

    is_error_res.add_incoming(is_error, bb_unboxed)
    is_error_res.add_incoming(cgutils.true_bit, bb_else)

    # cleanup
    c.pyapi.decref(set_type)
    c.pyapi.decref(key_type)

    return NativeValue(setobj_res, is_error=is_error_res)


@type_callable(SetType)
def typedset_call(context):
    """
    Defines typing logic for ``Set()``.
    Produces set[undefined]
    """
    def typer():
        return types.SetType(types.undefined)
    return typer


@overload(numba_typeref_ctor)
def impl_numba_typeref_ctor(cls):
    """
    Defines ``Set()``, the type-inferred version of the set ctor.

    Parameters
    ----------
    cls : TypeRef
        Expecting a TypeRef of a precise SetType.

    See also: `redirect_type_ctor` in numba/cpython/bulitins.py
    """
    set_ty = cls.instance_type
    if not isinstance(set_ty, types.SetType):
        msg = "expecting a SetType but got {}".format(set_ty)
        return  # reject
    # Ensure the set is precisely typed.
    if not set_ty.is_precise():
        msg = "expecting a precise SetType but got {}".format(set_ty)
        raise errors.LoweringError(msg)

    key_type = types.TypeRef(set_ty.key_type)

    def impl(cls):
        # Simply call .empty() with the value types from *cls*
        return Set.empty(key_type)

    return impl
