"""
Python wrapper that connects CPython interpreter to the Numba typed-list.

This is the code that is used when creating typed lists outside of a `@jit`
context and when returning a typed-list from a `@jit` decorated function. It
basically a Python class that has a Numba allocated typed-list under the hood
and uses `@jit` functions to access it. Since it inherits from MutableSequence
it should really quack like the CPython `list`.

"""
from collections.abc import MutableSequence

from numba.core.types import ListType
from numba.core.imputils import numba_typeref_ctor
from numba.core.dispatcher import Dispatcher
from numba.core import types, config, cgutils
from numba import njit, typeof
from numba.core.extending import (
    overload,
    box,
    unbox,
    NativeValue,
    type_callable,
    overload_classmethod,
)
from numba.typed import listobject
from numba.core.errors import TypingError, LoweringError
from numba.core.typing.templates import Signature
import typing as pt


Int_or_Slice = pt.Union["pt.SupportsIndex", slice]


T_co = pt.TypeVar('T_co', covariant=True)


class _Sequence(pt.Protocol[T_co]):
    def __getitem__(self, i: int) -> T_co:
        ...

    def __len__(self) -> int:
        ...


DEFAULT_ALLOCATED = listobject.DEFAULT_ALLOCATED


@njit
def _make_list(itemty, allocated=DEFAULT_ALLOCATED):
    return listobject._as_meminfo(listobject.new_list(itemty,
                                                      allocated=allocated))


@njit
def _length(l):
    return len(l)


@njit
def _allocated(l):
    return l._allocated()


@njit
def _is_mutable(l):
    return l._is_mutable()


@njit
def _make_mutable(l):
    return l._make_mutable()


@njit
def _make_immutable(l):
    return l._make_immutable()


@njit
def _append(l, item):
    l.append(item)


@njit
def _setitem(l, i, item):
    l[i] = item


@njit
def _getitem(l, i):
    return l[i]


@njit
def _contains(l, item):
    return item in l


@njit
def _count(l, item):
    return l.count(item)


@njit
def _pop(l, i):
    return l.pop(i)


@njit
def _delitem(l, i):
    del l[i]


@njit
def _extend(l, iterable):
    return l.extend(iterable)


@njit
def _insert(l, i, item):
    l.insert(i, item)


@njit
def _remove(l, item):
    l.remove(item)


@njit
def _clear(l):
    l.clear()


@njit
def _reverse(l):
    l.reverse()


@njit
def _copy(l):
    return l.copy()


@njit
def _eq(t, o):
    return t == o


@njit
def _ne(t, o):
    return t != o


@njit
def _lt(t, o):
    return t < o


@njit
def _le(t, o):
    return t <= o


@njit
def _gt(t, o):
    return t > o


@njit
def _ge(t, o):
    return t >= o


@njit
def _index(l, item, start, end):
    return l.index(item, start, end)


@njit
def _sort(l, key, reverse):
    return l.sort(key, reverse)


def _from_meminfo_ptr(ptr, listtype):
    return List(meminfo=ptr, lsttype=listtype)


T = pt.TypeVar('T')
T_or_ListT = pt.Union[T, 'List[T]']


class List(MutableSequence, pt.Generic[T]):
    """A typed-list usable in Numba compiled functions.

    Implements the MutableSequence interface.
    """

    _legal_kwargs = ["lsttype", "meminfo", "allocated"]

    def __new__(cls,
                *args,
                lsttype=None,
                meminfo=None,
                allocated=DEFAULT_ALLOCATED,
                **kwargs):
        if config.DISABLE_JIT:
            return list(*args, **kwargs)
        else:
            return object.__new__(cls)

    @classmethod
    def empty_list(cls, item_type, allocated=DEFAULT_ALLOCATED):
        """Create a new empty List.

        Parameters
        ----------
        item_type: Numba type
            type of the list item.
        allocated: int
            number of items to pre-allocate
        """
        if config.DISABLE_JIT:
            return list()
        else:
            return cls(lsttype=ListType(item_type), allocated=allocated)

    def __init__(self, *args, **kwargs):
        """
        For users, the constructor does not take any parameters.
        The keyword arguments are for internal use only.

        Parameters
        ----------
        args: iterable
            The iterable to initialize the list from
        lsttype : numba.core.types.ListType; keyword-only
            Used internally for the list type.
        meminfo : MemInfo; keyword-only
            Used internally to pass the MemInfo object when boxing.
        allocated: int; keyword-only
            Used internally to pre-allocate space for items
        """
        illegal_kwargs = any((kw not in self._legal_kwargs for kw in kwargs))
        if illegal_kwargs or args and kwargs:
            raise TypeError("List() takes no keyword arguments")
        if kwargs:
            self._list_type, self._opaque = self._parse_arg(**kwargs)
        else:
            self._list_type = None
            if args:
                if not 0 <= len(args) <= 1:
                    raise TypeError(
                        "List() expected at most 1 argument, got {}"
                        .format(len(args))
                    )
                iterable = args[0]
                # Special case Numpy scalars or anything that quacks like a
                # NumPy Array.
                if hasattr(iterable, "ndim") and iterable.ndim == 0:
                    self.append(iterable.item())
                else:
                    try:
                        iter(iterable)
                    except TypeError:
                        raise TypeError("List() argument must be iterable")
                    for i in args[0]:
                        self.append(i)

    def _parse_arg(self, lsttype, meminfo=None, allocated=DEFAULT_ALLOCATED):
        if not isinstance(lsttype, ListType):
            raise TypeError('*lsttype* must be a ListType')

        if meminfo is not None:
            opaque = meminfo
        else:
            opaque = _make_list(lsttype.item_type, allocated=allocated)
        return lsttype, opaque

    @property
    def _numba_type_(self):
        if self._list_type is None:
            raise TypeError("invalid operation on untyped list")
        return self._list_type

    @property
    def _typed(self):
        """Returns True if the list is typed.
        """
        return self._list_type is not None

    @property
    def _dtype(self):
        if not self._typed:
            raise RuntimeError("invalid operation on untyped list")
        return self._list_type.dtype

    def _initialise_list(self, item):
        lsttype = types.ListType(typeof(item))
        self._list_type, self._opaque = self._parse_arg(lsttype)

    def __len__(self) -> int:
        if not self._typed:
            return 0
        else:
            return _length(self)

    def _allocated(self):
        if not self._typed:
            return DEFAULT_ALLOCATED
        else:
            return _allocated(self)

    def _is_mutable(self):
        return _is_mutable(self)

    def _make_mutable(self):
        return _make_mutable(self)

    def _make_immutable(self):
        return _make_immutable(self)

    def __eq__(self, other):
        return _eq(self, other)

    def __ne__(self, other):
        return _ne(self, other)

    def __lt__(self, other):
        return _lt(self, other)

    def __le__(self, other):
        return _le(self, other)

    def __gt__(self, other):
        return _gt(self, other)

    def __ge__(self, other):
        return _ge(self, other)

    def append(self, item: T) -> None:
        if not self._typed:
            self._initialise_list(item)
        _append(self, item)

    # noqa F811 comments required due to github.com/PyCQA/pyflakes/issues/592
    # noqa E704 required to follow overload style of using ... in the same line
    @pt.overload  # type: ignore[override]
    def __setitem__(self, i: int, o: T) -> None: ...  # noqa: F811, E704
    @pt.overload
    def __setitem__(self, s: slice, o: 'List[T]') -> None: ...  # noqa: F811, E704, E501

    def __setitem__(self, i: Int_or_Slice, item: T_or_ListT) -> None:  # noqa: F811, E501
        if not self._typed:
            self._initialise_list(item)
        _setitem(self, i, item)

    # noqa F811 comments required due to github.com/PyCQA/pyflakes/issues/592
    # noqa E704 required to follow overload style of using ... in the same line
    @pt.overload
    def __getitem__(self, i: int) -> T: ...  # noqa: F811, E704
    @pt.overload
    def __getitem__(self, i: slice) -> 'List[T]': ...  # noqa: F811, E704

    def __getitem__(self, i: Int_or_Slice) -> T_or_ListT:  # noqa: F811
        if not self._typed:
            raise IndexError
        else:
            return _getitem(self, i)

    def __iter__(self) -> pt.Iterator[T]:
        for i in range(len(self)):
            yield self[i]

    def __contains__(self, item: T) -> bool:  # type: ignore[override]
        return _contains(self, item)

    def __delitem__(self, i: Int_or_Slice) -> None:
        _delitem(self, i)

    def insert(self, i: int, item: T) -> None:
        if not self._typed:
            self._initialise_list(item)
        _insert(self, i, item)

    def count(self, item: T) -> int:
        return _count(self, item)

    def pop(self, i: "pt.SupportsIndex" = -1) -> T:
        return _pop(self, i)

    def extend(self, iterable: "_Sequence[T]") -> None: #type: ignore[override]
        # Empty iterable, do nothing
        if len(iterable) == 0:
            return None
        if not self._typed:
            # Need to get the first element of the iterable to initialise the
            # type of the list. FIXME: this may be a problem if the iterable
            # can not be sliced.
            self._initialise_list(iterable[0])
        return _extend(self, iterable)

    def remove(self, item: T) -> None:
        return _remove(self, item)

    def clear(self):
        return _clear(self)

    def reverse(self):
        return _reverse(self)

    def copy(self):
        return _copy(self)

    def index(self, item: T, start: pt.Optional[int] = None,
              stop: pt.Optional[int] = None) -> int:
        return _index(self, item, start, stop)

    def sort(self, key=None, reverse=False):
        """Sort the list inplace.

        See also ``list.sort()``
        """
        # If key is not already a dispatcher object, make it so
        if callable(key) and not isinstance(key, Dispatcher):
            key = njit(key)
        return _sort(self, key, reverse)

    def __str__(self):
        buf = []
        for x in self:
            buf.append("{}".format(x))
        # Check whether the code was invoked from IPython shell
        try:
            get_ipython
            return '[{0}, ...]'.format(', '.join(buf[:1000]))
        except (NameError, IndexError):
            return '[{0}]'.format(', '.join(buf))

    def __repr__(self):
        body = str(self)
        prefix = str(self._list_type) if self._typed else "ListType[Undefined]"
        return "{prefix}({body})".format(prefix=prefix, body=body)


@overload_classmethod(ListType, 'empty_list')
def typedlist_empty(cls, item_type, allocated=DEFAULT_ALLOCATED):
    if cls.instance_type is not ListType:
        return

    def impl(cls, item_type, allocated=DEFAULT_ALLOCATED):
        return listobject.new_list(item_type, allocated=allocated)

    return impl


@box(types.ListType)
def box_lsttype(typ, val, c):
    context = c.context
    builder = c.builder

    # XXX deduplicate
    ctor = cgutils.create_struct_proxy(typ)
    lstruct = ctor(context, builder, value=val)
    # Returns the plain MemInfo
    boxed_meminfo = c.box(
        types.MemInfoPointer(types.voidptr),
        lstruct.meminfo,
    )

    modname = c.context.insert_const_string(
        c.builder.module, 'numba.typed.typedlist',
    )
    typedlist_mod = c.pyapi.import_module_noblock(modname)
    fmp_fn = c.pyapi.object_getattr_string(typedlist_mod, '_from_meminfo_ptr')

    lsttype_obj = c.pyapi.unserialize(c.pyapi.serialize_object(typ))

    result_var = builder.alloca(c.pyapi.pyobj)
    builder.store(cgutils.get_null_value(c.pyapi.pyobj), result_var)

    with builder.if_then(cgutils.is_not_null(builder, lsttype_obj)):
        res = c.pyapi.call_function_objargs(
            fmp_fn, (boxed_meminfo, lsttype_obj),
        )
        c.pyapi.decref(fmp_fn)
        c.pyapi.decref(typedlist_mod)
        c.pyapi.decref(boxed_meminfo)
        builder.store(res, result_var)
    return builder.load(result_var)


@unbox(types.ListType)
def unbox_listtype(typ, val, c):
    context = c.context
    builder = c.builder

    # Check that `type(val) is Dict`
    list_type = c.pyapi.unserialize(c.pyapi.serialize_object(List))
    valtype = c.pyapi.object_type(val)
    same_type = builder.icmp_unsigned("==", valtype, list_type)

    with c.builder.if_else(same_type) as (then, orelse):
        with then:
            miptr = c.pyapi.object_getattr_string(val, '_opaque')

            native = c.unbox(types.MemInfoPointer(types.voidptr), miptr)

            mi = native.value
            ctor = cgutils.create_struct_proxy(typ)
            lstruct = ctor(context, builder)

            data_pointer = context.nrt.meminfo_data(builder, mi)
            data_pointer = builder.bitcast(
                data_pointer,
                listobject.ll_list_type.as_pointer(),
            )

            lstruct.data = builder.load(data_pointer)
            lstruct.meminfo = mi

            lstobj = lstruct._getvalue()
            c.pyapi.decref(miptr)
            bb_unboxed = c.builder.basic_block

        with orelse:
            # Raise error on incorrect type
            c.pyapi.err_format(
                "PyExc_TypeError",
                "can't unbox a %S as a %S",
                valtype, list_type,
            )
            bb_else = c.builder.basic_block

    # Phi nodes to gather the output
    lstobj_res = c.builder.phi(lstobj.type)
    is_error_res = c.builder.phi(cgutils.bool_t)

    lstobj_res.add_incoming(lstobj, bb_unboxed)
    lstobj_res.add_incoming(lstobj.type(None), bb_else)

    is_error_res.add_incoming(cgutils.false_bit, bb_unboxed)
    is_error_res.add_incoming(cgutils.true_bit, bb_else)

    # cleanup
    c.pyapi.decref(list_type)
    c.pyapi.decref(valtype)

    return NativeValue(lstobj_res, is_error=is_error_res)


#
# The following contains the logic for the type-inferred constructor
#

def _guess_dtype(iterable):
    """Guess the correct dtype of the iterable type. """
    if not isinstance(iterable, types.IterableType):
        raise TypingError(
            "List() argument must be iterable")
    # Special case for nested NumPy arrays.
    elif isinstance(iterable, types.Array) and iterable.ndim > 1:
        return iterable.copy(ndim=iterable.ndim - 1, layout='A')
    elif hasattr(iterable, "dtype"):
        return iterable.dtype
    elif hasattr(iterable, "yield_type"):
        return iterable.yield_type
    elif isinstance(iterable, types.UnicodeType):
        return iterable
    elif isinstance(iterable, types.DictType):
        return iterable.key_type
    else:
        # This should never happen, since the 'dtype' of any iterable
        # should have determined above.
        raise TypingError(
            "List() argument does not have a suitable dtype")


@type_callable(ListType)
def typedlist_call(context):
    """Defines typing logic for ``List()`` and ``List(iterable)``.

    If no argument is given, the returned typer types a new typed-list with an
    undefined item type. If a single argument is given it must be iterable with
    a guessable 'dtype'. In this case, the typer types a new typed-list with
    the type set to the 'dtype' of the iterable arg.

    Parameters
    ----------
    arg : single iterable (optional)
        The single optional argument.

    Returns
    -------
    typer : function
        A typer suitable to type constructor calls.

    Raises
    ------
    The returned typer raises a TypingError in case of unsuitable arguments.

    """

    class Typer(object):

        def attach_sig(self):
            from inspect import signature as mypysig

            def mytyper(iterable):
                pass
            self.pysig = mypysig(mytyper)

        def __call__(self, *args, **kwargs):
            if kwargs:
                raise TypingError(
                    "List() takes no keyword arguments"
                )
            elif args:
                if not 0 <= len(args) <= 1:
                    raise TypingError(
                        "List() expected at most 1 argument, got {}"
                        .format(len(args))
                    )
                rt = types.ListType(_guess_dtype(args[0]))
                self.attach_sig()
                return Signature(rt, args, None, pysig=self.pysig)
            else:
                item_type = types.undefined
                return types.ListType(item_type)

    return Typer()


@overload(numba_typeref_ctor)
def impl_numba_typeref_ctor(cls, *args):
    """Defines lowering for ``List()`` and ``List(iterable)``.

    This defines the lowering logic to instantiate either an empty typed-list
    or a typed-list initialised with values from a single iterable argument.

    Parameters
    ----------
    cls : TypeRef
        Expecting a TypeRef of a precise ListType.
    args: tuple
        A tuple that contains a single iterable (optional)

    Returns
    -------
    impl : function
        An implementation suitable for lowering the constructor call.

    See also: `redirect_type_ctor` in numba/cpython/bulitins.py
    """
    list_ty = cls.instance_type
    if not isinstance(list_ty, types.ListType):
        return  # reject
    # Ensure the list is precisely typed.
    if not list_ty.is_precise():
        msg = "expecting a precise ListType but got {}".format(list_ty)
        raise LoweringError(msg)

    item_type = types.TypeRef(list_ty.item_type)
    if args:
        # special case 0d Numpy arrays
        if isinstance(args[0], types.Array) and args[0].ndim == 0:
            def impl(cls, *args):
                # Instantiate an empty list and populate it with the single
                # value from the array.
                r = List.empty_list(item_type)
                r.append(args[0].item())
                return r
        else:
            def impl(cls, *args):
                # Instantiate an empty list and populate it with values from
                # the iterable.
                r = List.empty_list(item_type)
                for i in args[0]:
                    r.append(i)
                return r
    else:
        def impl(cls, *args):
            # Simply call .empty_list with the item type from *cls*
            return List.empty_list(item_type)

    return impl
