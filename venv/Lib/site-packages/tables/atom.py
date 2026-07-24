"""Atom classes for describing dataset contents."""

from __future__ import annotations

import re
import pickle
import inspect
import warnings
from typing import Any, NoReturn, dataclass_transform
from collections.abc import Callable

import numpy as np
from numpy.typing import DTypeLike

from .utils import SizeType
from .misc.enum import Enum
from .exceptions import FlavorWarning

Shape = tuple[SizeType, ...]

__docformat__ = "reStructuredText"
"""The format of documentation strings in this module."""

all_types = set()  # filled as atom classes are created
"""Set of all PyTables types."""

atom_map: dict[str, Any] = {}  # filled as atom classes are created
"""Maps atom kinds to item sizes and atom classes.

If there is a fixed set of possible item sizes for a given kind, the
kind maps to another mapping from item size in bytes to atom class.
Otherwise, the kind maps directly to the atom class.
"""

deftype_from_kind = {}  # filled as atom classes are created
"""Maps atom kinds to their default atom type (if any)."""


_type_re = re.compile(r"^([a-z]+)([0-9]*)$")


def split_type(type_: str) -> tuple[str, int | None]:
    """Split a PyTables type into a PyTables kind and an item size.

    Returns a tuple of (kind, itemsize). If no item size is present in the type
    (in the form of a precision), the returned item size is None::

        >>> split_type('int32')
        ('int', 4)
        >>> split_type('string')
        ('string', None)
        >>> split_type('int20')
        Traceback (most recent call last):
        ...
        ValueError: precision must be a multiple of 8: 20
        >>> split_type('foo bar')
        Traceback (most recent call last):
        ...
        ValueError: malformed type: 'foo bar'

    """
    match = _type_re.match(type_)
    if not match:
        raise ValueError("malformed type: %r" % type_)
    kind, precision = match.groups()
    itemsize = None
    if precision:
        precision = int(precision)
        itemsize, remainder = divmod(precision, 8)
        if remainder:  # 0 could be a valid item size
            raise ValueError(
                "precision must be a multiple of 8: %d" % precision
            )
    return (kind, itemsize)


def _invalid_itemsize_error(
    kind: str, itemsize: int, itemsizes: list[int]
) -> ValueError:
    isizes = sorted(itemsizes)
    return ValueError(
        "invalid item size for kind ``%s``: %r; "
        "it must be one of ``%r``" % (kind, itemsize, isizes)
    )


def _normalize_shape(shape: Shape | np.integer | int) -> Shape:
    """Check that the `shape` is safe to be used and return it as a tuple."""
    if isinstance(shape, (np.integer, int)):
        if shape < 1:
            raise ValueError("shape value must be greater than 0: %d" % shape)
        shape = (shape,)  # N is a shorthand for (N,)
    try:
        shape = tuple(shape)
    except TypeError:
        raise TypeError(f"shape must be an integer or sequence: {shape!r}")

    # XXX Get from HDF5 library if possible.
    # HDF5 does not support ranks greater than 32
    if len(shape) > 32:
        raise ValueError(f"shapes with rank > 32 are not supported: {shape!r}")

    return tuple(SizeType(s) for s in shape)


def _normalize_default(value: Any, dtype: DTypeLike) -> np.ndarray:
    """Return `value` as a valid default of NumPy type `dtype`."""
    # Create NumPy objects as defaults
    # This is better in order to serialize them as attributes
    if value is None:
        value = 0
    basedtype = dtype.base
    try:
        default = np.array(value, dtype=basedtype)
    except ValueError:
        array = np.array(value)
        if array.shape != basedtype.shape:
            raise
        # Maybe nested dtype with "scalar" value.
        default = np.array(value, dtype=basedtype.base)
    # 0-dim arrays will be representented as NumPy scalars
    # (PyTables attribute convention)
    if default.shape == ():
        default = default[()]
    return default


def _cmp_dispatcher(other_method_name: str) -> Callable[[Any, Any], bool]:
    """Dispatch comparisons to a method of the *other* object.

    Returns a new *rich comparison* method which dispatches calls to
    the method `other_method_name` of the *other* object.  If there is
    no such method in the object, ``False`` is returned.

    This is part of the implementation of a double dispatch pattern.
    """

    def dispatched_cmp(self, other) -> bool:
        try:
            other_method: Callable[[Any], bool] = getattr(
                other, other_method_name
            )
        except AttributeError:
            return False
        return other_method(self)

    return dispatched_cmp


@dataclass_transform()
class MetaAtom(type):
    """Atom metaclass.

    This metaclass ensures that data about atom classes gets inserted
    into the suitable registries.

    """

    kind: str

    def __init__(cls, name: str, bases: tuple, dict_: dict[str, Any]) -> None:
        super().__init__(name, bases, dict_)

        kind = dict_.get("kind")
        itemsize = dict_.get("itemsize")
        type_ = dict_.get("type")
        deftype = dict_.get("_deftype")

        if kind and deftype:
            deftype_from_kind[kind] = deftype

        if type_:
            all_types.add(type_)

        if kind and itemsize and not hasattr(itemsize, "__int__"):
            # Atom classes with a non-fixed item size do have an
            # ``itemsize``, but it's not a number (e.g. property).
            atom_map[kind] = cls
            return

        if kind:  # first definition of kind, make new entry
            atom_map[kind] = {}

        if itemsize and hasattr(itemsize, "__int__"):  # fixed
            kind = cls.kind  # maybe from superclasses
            atom_map[kind][int(itemsize)] = cls


class Atom(metaclass=MetaAtom):
    """Defines the type of atomic cells stored in a dataset.

    The meaning of *atomic* is that individual elements of a cell can
    not be extracted directly by indexing (i.e.  __getitem__()) the
    dataset; e.g. if a dataset has shape (2, 2) and its atoms have
    shape (3,), to get the third element of the cell at (1, 0) one
    should use dataset[1,0][2] instead of dataset[1,0,2].

    The Atom class is meant to declare the different properties of the
    *base element* (also known as *atom*) of CArray, EArray and
    VLArray datasets, although they are also used to describe the base
    elements of Array datasets. Atoms have the property that their
    length is always the same.  However, you can grow datasets along
    the extensible dimension in the case of EArray or put a variable
    number of them on a VLArray row. Moreover, they are not restricted
    to scalar values, and they can be *fully multidimensional
    objects*.

    Parameters
    ----------
    nptype : str or np.dtype
        Sets the Numpy data type of the atom.
    shape : tuple
        Sets the shape of the atom. An integer shape of
        N is equivalent to the tuple (N,).
    dflt : Any
        Sets the default value for the atom.

    The following are the public methods and attributes of the Atom class.

    Notes
    -----
    A series of descendant classes are offered in order to make the
    use of these element descriptions easier. You should use a
    particular Atom descendant class whenever you know the exact type
    you will need when writing your code. Otherwise, you may use one
    of the Atom.from_*() factory Methods.

    .. rubric:: Atom attributes

    .. attribute:: dflt

        The default value of the atom.

        If the user does not supply a value for an element while
        filling a dataset, this default value will be written to disk.
        If the user supplies a scalar value for a multidimensional
        atom, this value is automatically *broadcast* to all the items
        in the atom cell. If dflt is not supplied, an appropriate zero
        value (or *null* string) will be chosen by default.  Please
        note that default values are kept internally as NumPy objects.

    .. attribute:: dtype

        The NumPy dtype that most closely matches this atom.

    .. attribute:: itemsize

        Size in bytes of a single item in the atom.
        Specially useful for atoms of the string kind.

    .. attribute:: kind

        The PyTables kind of the atom (a string).

    .. attribute:: shape

        The shape of the atom (a tuple for scalar atoms).

    .. attribute:: type

        The PyTables type of the atom (a string).

        Atoms can be compared with atoms and other objects for
        strict (in)equality without having to compare individual
        attributes::

            >>> atom1 = StringAtom(itemsize=10)  # same as ``atom2``
            >>> atom2 = Atom.from_kind('string', 10)  # same as ``atom1``
            >>> atom3 = IntAtom()
            >>> bool(atom1 == 'foo')
            False
            >>> bool(atom1 == atom2)
            True
            >>> bool(atom2 != atom1)
            False
            >>> bool(atom1 == atom3)
            False
            >>> bool(atom3 != atom2)
            True

    """

    dflt: Any

    dtype: np.dtype

    itemsize: int

    kind: str

    shape: Shape

    type: str  # noqa: A003

    @classmethod
    def prefix(cls) -> str:
        """Return the atom class prefix."""
        cname = cls.__name__
        return cname[: cname.rfind("Atom")]

    @classmethod
    def from_sctype(
        cls, sctype: str | np.dtype, shape: Shape = (), dflt: Any = None
    ) -> Atom:
        """Create an Atom from a NumPy scalar type sctype.

        Optional shape and default value may be specified as the
        shape and dflt
        arguments, respectively. Information in the
        sctype not represented in an Atom is ignored::

            >>> import numpy as np
            >>> Atom.from_sctype(np.int16, shape=(2, 2))
            Int16Atom(shape=(2, 2), dflt=0)
            >>> Atom.from_sctype('S5', dflt='hello')
            Traceback (most recent call last):
            ...
            ValueError: unknown NumPy scalar type: 'S5'
            >>> Atom.from_sctype('float64')
            Float64Atom(shape=(), dflt=0.0)

        """
        if not isinstance(sctype, type) or not issubclass(sctype, np.generic):
            assert isinstance(sctype, str)
            if "," in sctype:
                raise ValueError(f"unknown NumPy scalar type: {sctype!r}")
            try:
                dtype = np.dtype(sctype)
            except TypeError:
                raise ValueError(
                    f"unknown NumPy scalar type: {sctype!r}"
                ) from None
            if issubclass(dtype.type, np.flexible) and dtype.itemsize > 0:
                raise ValueError(f"unknown NumPy scalar type: {sctype!r}")

            sctype_resolved = dtype.type
        else:
            sctype_resolved = sctype
        return cls.from_dtype(np.dtype((sctype_resolved, shape)), dflt)

    @classmethod
    def from_dtype(cls, dtype: np.dtype, dflt: Any = None) -> Atom:
        """Create an Atom from a NumPy dtype.

        An optional default value may be specified as the dflt
        argument. Information in the dtype not represented in an Atom is
        ignored::

            >>> import numpy as np
            >>> Atom.from_dtype(np.dtype((np.int16, (2, 2))))
            Int16Atom(shape=(2, 2), dflt=0)
            >>> Atom.from_dtype(np.dtype('float64'))
            Float64Atom(shape=(), dflt=0.0)

        Note: for easier use in Python 3, where all strings lead to the
        Unicode dtype, this dtype will also generate a StringAtom. Since
        this is only viable for strings that are castable as ascii, a
        warning is issued.

            >>> Atom.from_dtype(np.dtype('U20')) # doctest: +SKIP
            Atom.py:392: FlavorWarning: support for unicode type is very
                limited, and only works for strings that can be cast as ascii
            StringAtom(itemsize=20, shape=(), dflt=b'')

        """
        basedtype = dtype.base
        shape = tuple(SizeType(i) for i in dtype.shape)
        if basedtype.names:
            raise ValueError(
                "compound data types are not supported: %r" % dtype
            )
        if basedtype.shape != ():
            raise ValueError("nested data types are not supported: %r" % dtype)
        if basedtype.kind == "S":  # can not reuse something like 'string80'
            itemsize = basedtype.itemsize
            return cls.from_kind("string", itemsize, shape, dflt)
        elif basedtype.kind == "U":
            # workaround for unicode type (standard string type in Python 3)
            warnings.warn(
                "support for unicode type is very limited, and "
                "only works for strings that can be cast as ascii",
                FlavorWarning,
            )
            itemsize = basedtype.itemsize // 4
            assert (
                str(itemsize) in basedtype.str
            ), "something went wrong in handling unicode."
            return cls.from_kind("string", itemsize, shape, dflt)
        # Most NumPy types have direct correspondence with PyTables types.
        return cls.from_type(basedtype.name, shape, dflt)

    @classmethod
    def from_type(
        cls, type_: str, shape: Shape = (), dflt: Any = None
    ) -> Atom:
        """Create an Atom from a PyTables type.

        Optional shape and default value may be specified as the
        shape and dflt arguments, respectively::

            >>> Atom.from_type('bool')
            BoolAtom(shape=(), dflt=False)
            >>> Atom.from_type('int16', shape=(2, 2))
            Int16Atom(shape=(2, 2), dflt=0)
            >>> Atom.from_type('string40', dflt='hello')
            Traceback (most recent call last):
            ...
            ValueError: unknown type: 'string40'
            >>> Atom.from_type('Float64')
            Traceback (most recent call last):
            ...
            ValueError: unknown type: 'Float64'

        """
        if type_ not in all_types:
            raise ValueError(f"unknown type: {type_!r}")
        kind, itemsize = split_type(type_)
        return cls.from_kind(kind, itemsize, shape, dflt)

    @classmethod
    def from_kind(
        cls,
        kind: str,
        itemsize: int | None = None,
        shape: Shape = (),
        dflt: Any = None,
    ) -> Atom:
        """Create an Atom from a PyTables kind.

        Optional item size, shape and default value may be
        specified as the itemsize, shape and dflt
        arguments, respectively. Bear in mind that not all atoms support
        a default item size::

            >>> Atom.from_kind('int', itemsize=2, shape=(2, 2))
            Int16Atom(shape=(2, 2), dflt=0)
            >>> Atom.from_kind('int', shape=(2, 2))
            Int32Atom(shape=(2, 2), dflt=0)
            >>> Atom.from_kind('int', shape=1)
            Int32Atom(shape=(1,), dflt=0)
            >>> Atom.from_kind('string', dflt=b'hello')
            Traceback (most recent call last):
            ...
            ValueError: no default item size for kind ``string``
            >>> Atom.from_kind('Float')
            Traceback (most recent call last):
            ...
            ValueError: unknown kind: 'Float'

        Moreover, some kinds with atypical constructor signatures
        are not supported; you need to use the proper
        constructor::

            >>> Atom.from_kind('enum') #doctest: +ELLIPSIS
            Traceback (most recent call last):
            ...
            ValueError: the ``enum`` kind is not supported...

        """
        kwargs: dict[str, Any] = {"shape": shape}
        if kind not in atom_map:
            raise ValueError(f"unknown kind: {kind!r}")
        # This incompatibility detection may get out-of-date and is
        # too hard-wired, but I couldn't come up with something
        # smarter.  -- Ivan (2007-02-08)
        if kind in ["enum"]:
            raise ValueError(
                "the ``%s`` kind is not supported; "
                "please use the appropriate constructor" % kind
            )
        # If no `itemsize` is given, try to get the default type of the
        # kind (which has a fixed item size).
        if itemsize is None:
            if kind not in deftype_from_kind:
                raise ValueError("no default item size for kind ``%s``" % kind)
            type_ = deftype_from_kind[kind]
            kind, itemsize = split_type(type_)
        kdata = atom_map[kind]
        # Look up the class and set a possible item size.
        if hasattr(kdata, "kind"):  # atom class: non-fixed item size
            atomclass = kdata
            kwargs["itemsize"] = itemsize
        else:  # dictionary: fixed item size
            if itemsize not in kdata:
                raise _invalid_itemsize_error(kind, itemsize, kdata)
            atomclass = kdata[itemsize]
        # Only set a `dflt` argument if given (`None` may not be understood).
        if dflt is not None:
            kwargs["dflt"] = dflt

        return atomclass(**kwargs)

    @property
    def size(self) -> int:
        """Total size in bytes of the atom."""
        return self.dtype.itemsize

    @property
    def recarrtype(self) -> str:
        """Return the string type to be used in `numpy.rec.array()`."""
        return str(self.dtype.shape) + self.dtype.base.str[1:]

    @property
    def ndim(self) -> int:
        """Return the number of dimensions of the atom.

        .. versionadded:: 2.4
        """
        return len(self.shape)

    def __init__(
        self, nptype: str | np.dtype, shape: Shape, dflt: Any
    ) -> None:
        if not hasattr(self, "type"):
            raise NotImplementedError(
                f"``{self.__class__.__name__}`` is an abstract class; "
                f"please use one of its subclasses"
            )
        self.shape = shape = _normalize_shape(shape)
        """The shape of the atom (a tuple for scalar atoms)."""
        # Curiously enough, NumPy isn't generally able to accept NumPy
        # integers in a shape. ;(
        npshape = tuple(int(s) for s in shape)
        self.dtype = dtype = np.dtype((nptype, npshape))
        """The NumPy dtype that most closely matches this atom."""
        self.dflt = _normalize_default(dflt, dtype)
        """The default value of the atom.

        If the user does not supply a value for an element while
        filling a dataset, this default value will be written to
        disk. If the user supplies a scalar value for a
        multidimensional atom, this value is automatically *broadcast*
        to all the items in the atom cell. If dflt is not supplied, an
        appropriate zero value (or *null* string) will be chosen by
        default.  Please note that default values are kept internally
        as NumPy objects."""

    def __repr__(self) -> str:
        args = f"shape={self.shape}, dflt={self.dflt!r}"
        if not hasattr(self.__class__.itemsize, "__int__"):  # non-fixed
            args = f"itemsize={self.itemsize}, {args}"
        return f"{self.__class__.__name__}({args})"

    __eq__ = _cmp_dispatcher("_is_equal_to_atom")

    def __ne__(self, other: Atom) -> bool:
        return not self.__eq__(other)

    # XXX: API incompatible change for PyTables 3 line
    # Overriding __eq__ blocks inheritance of __hash__ in 3.x
    # def __hash__(self):
    #    return hash((self.__class__, self.type, self.shape, self.itemsize,
    #                 self.dflt))

    def copy(self, **override) -> Atom:
        """Get a copy of the atom, possibly overriding some arguments.

        Constructor arguments to be overridden must be passed as
        keyword arguments::

            >>> atom1 = Int32Atom(shape=12)
            >>> atom2 = atom1.copy()
            >>> print(atom1)
            Int32Atom(shape=(12,), dflt=0)
            >>> print(atom2)
            Int32Atom(shape=(12,), dflt=0)
            >>> atom1 is atom2
            False
            >>> atom3 = atom1.copy(shape=(2, 2))
            >>> print(atom3)
            Int32Atom(shape=(2, 2), dflt=0)
            >>> atom1.copy(foobar=42) #doctest: +ELLIPSIS
            Traceback (most recent call last):
            ...
            TypeError: ...__init__() got an unexpected keyword argument ...

        """
        newargs = self._get_init_args()
        newargs.update(override)
        return self.__class__(**newargs)

    def _get_init_args(self) -> dict[str, Any]:
        """Get a dictionary of instance constructor arguments.

        This implementation works on classes which use the same names
        for both constructor arguments and instance attributes.

        """
        signature = inspect.signature(self.__init__)
        parameters = signature.parameters
        args = [
            arg
            for arg, p in parameters.items()
            if p.kind is p.POSITIONAL_OR_KEYWORD
        ]

        return {arg: getattr(self, arg) for arg in args if arg != "self"}

    def _is_equal_to_atom(self, atom: Atom) -> bool:
        """Return True if the object is equal to the given `atom`."""
        return (
            self.type == atom.type
            and self.shape == atom.shape
            and self.itemsize == atom.itemsize
            and np.all(self.dflt == atom.dflt)
        )


def _abstract_atom_init(
    deftype: str, defvalue: Any
) -> Callable[[Atom, int | None, Shape, Any], None]:
    """Return a constructor for an abstract `Atom` class."""
    defitemsize = split_type(deftype)[1]

    def __init__(  # noqa: N807
        self: Atom,
        itemsize: int | None = defitemsize,
        shape: Shape = (),
        dflt: Any = defvalue,
    ) -> None:
        assert self.kind in atom_map
        try:
            atomclass = atom_map[self.kind][itemsize]
        except KeyError:
            raise _invalid_itemsize_error(
                self.kind, itemsize, atom_map[self.kind]
            )
        self.__class__ = atomclass
        atomclass.__init__(self, shape, dflt)

    return __init__


class StringAtom(Atom):  # type: ignore[misc]
    """Defines an atom of type string.

    The item size is the *maximum* length in characters of strings.

    """

    kind: str = "string"
    type: str = "string"  # noqa: A003
    _defvalue: bytes = b""

    @property  # type: ignore[misc]
    def itemsize(self) -> int:  # type: ignore[override]
        """Size in bytes of a sigle item in the atom."""
        return self.dtype.base.itemsize

    def __init__(
        self, itemsize: int, shape: Shape = (), dflt: str | bytes = _defvalue
    ) -> None:
        if not hasattr(itemsize, "__int__") or int(itemsize) < 0:
            raise ValueError(
                f"invalid item size for kind ``string``: {itemsize!r}; "
                f"it must be a positive integer"
            )
        Atom.__init__(self, f"S{itemsize}", shape, dflt)


class BoolAtom(Atom):  # type: ignore[misc]
    """Defines an atom of type bool."""

    kind: str = "bool"
    itemsize: int = 1
    type: str = "bool"  # noqa: A003
    _deftype = "bool8"
    _defvalue = False

    def __init__(self, shape: Shape = (), dflt: bool = _defvalue) -> None:
        Atom.__init__(self, self.type, shape, dflt)


class IntAtom(Atom):  # type: ignore[misc]
    """Defines an atom of a signed integral type (int kind)."""

    kind: str = "int"
    signed: bool = True
    _deftype = "int32"
    _defvalue = 0
    __init__ = _abstract_atom_init(
        _deftype, _defvalue
    )  # type: ignore[assignment]


class UIntAtom(Atom):  # type: ignore[misc]
    """Defines an atom of an unsigned integral type (uint kind)."""

    kind: str = "uint"
    signed: bool = False
    _deftype = "uint32"
    _defvalue = 0
    __init__ = _abstract_atom_init(
        _deftype, _defvalue
    )  # type: ignore[assignment]


class FloatAtom(Atom):  # type: ignore[misc]
    """Defines an atom of a floating point type (float kind)."""

    kind: str = "float"
    _deftype = "float64"
    _defvalue = 0.0
    __init__ = _abstract_atom_init(
        _deftype, _defvalue
    )  # type: ignore[assignment]


class Int8Atom(IntAtom):  # type: ignore[misc]
    """Atom for 8 bit integers."""

    itemsize: int = 1
    type: str = "int8"  # noqa: A003

    def __init__(self, shape: Shape = (), dflt: int = 0) -> None:
        Atom.__init__(self, "int8", shape, dflt)


class Int16Atom(IntAtom):  # type: ignore[misc]
    """Atom for 12 bit integers."""

    itemsize: int = 2
    type: str = "int16"  # noqa: A003

    def __init__(self, shape: Shape = (), dflt: int = 0) -> None:
        Atom.__init__(self, "int16", shape, dflt)


class Int32Atom(IntAtom):  # type: ignore[misc]
    """Atom for 32 bit integers."""

    itemsize: int = 4
    type: str = "int32"  # noqa: A003

    def __init__(self, shape: Shape = (), dflt: int = 0) -> None:
        Atom.__init__(self, "int32", shape, dflt)


class Int64Atom(IntAtom):  # type: ignore[misc]
    """Atom for 64 bit integers."""

    itemsize: int = 8
    type: str = "int64"  # noqa: A003

    def __init__(self, shape: Shape = (), dflt: int = 0) -> None:
        Atom.__init__(self, "int64", shape, dflt)


class UInt8Atom(UIntAtom):  # type: ignore[misc]
    """Atom for 8 bit unsoged integers."""

    itemsize: int = 1
    type: str = "uint8"  # noqa: A003

    def __init__(self, shape: Shape = (), dflt: int = 0) -> None:
        Atom.__init__(self, "uint8", shape, dflt)


class UInt16Atom(UIntAtom):  # type: ignore[misc]
    """Atom for 16 bit unsigned integers."""

    itemsize: int = 2
    type: str = "uint16"  # noqa: A003

    def __init__(self, shape: Shape = (), dflt: int = 0) -> None:
        Atom.__init__(self, "uint16", shape, dflt)


class UInt32Atom(UIntAtom):  # type: ignore[misc]
    """Atom for 32 bit unsigned integers."""

    itemsize: int = 4
    type: str = "uint32"  # noqa: A003

    def __init__(self, shape: Shape = (), dflt: int = 0) -> None:
        Atom.__init__(self, "uint32", shape, dflt)


class UInt64Atom(UIntAtom):  # type: ignore[misc]
    """Atom for 16 bit unsigned integers."""

    itemsize: int = 8
    type: str = "uint64"  # noqa: A003

    def __init__(self, shape: Shape = (), dflt: int = 0) -> None:
        Atom.__init__(self, "uint64", shape, dflt)


if hasattr(np, "float16"):

    class Float16Atom(FloatAtom):  # type: ignore[misc]
        """FLoat 16 atom."""

        itemsize: int = 2
        type: str = "float16"  # noqa: A003

        def __init__(self, shape: Shape = (), dflt: float = 0.0) -> None:
            Atom.__init__(self, "float16", shape, dflt)


class Float32Atom(FloatAtom):  # type: ignore[misc]
    """Float 32 atom."""

    itemsize: int = 4
    type: str = "float32"  # noqa: A003

    def __init__(self, shape: Shape = (), dflt: float = 0.0) -> None:
        Atom.__init__(self, "float32", shape, dflt)


class Float64Atom(FloatAtom):  # type: ignore[misc]
    """Float 64 atom."""

    itemsize: int = 8
    type: str = "float64"  # noqa: A003

    def __init__(self, shape: Shape = (), dflt: float = 0.0) -> None:
        Atom.__init__(self, "float64", shape, dflt)


if hasattr(np, "float96"):

    class Float96Atom(FloatAtom):  # type: ignore[misc]
        """Float 96 atom."""

        itemsize: int = 12
        type: str = "float96"  # noqa: A003

        def __init__(self, shape: Shape = (), dflt: float = 0.0) -> None:
            Atom.__init__(self, "float96", shape, dflt)


if hasattr(np, "float128"):

    class Float128Atom(FloatAtom):  # type: ignore[misc]
        """Float 128 atom."""

        itemsize: int = 16
        type: str = "float128"  # noqa: A003

        def __init__(self, shape: Shape = (), dflt: float = 0.0) -> None:
            Atom.__init__(self, "float128", shape, dflt)


class ComplexAtom(Atom):
    """Defines an atom of kind complex.

    Allowed item sizes are 8 (single precision) and 16 (double precision). This
    class must be used instead of more concrete ones to avoid confusions with
    numarray-like precision specifications used in PyTables 1.X.

    """

    # This definition is a little more complex (no pun intended)
    # because, although the complex kind is a normal numerical one,
    # the usage of bottom-level classes is artificially forbidden.
    # Everything will be back to normality when people has stopped
    # using the old bottom-level complex classes.

    kind = "complex"
    _deftype = "complex128"
    _defvalue = 0j
    _isizes = [8, 16]

    @property  # type: ignore[misc]
    def itemsize(self) -> int:  # type: ignore[override]
        """Size in bytes of a sigle item in the atom."""
        return self.dtype.base.itemsize

    # Only instances have a `type` attribute, so complex types must be
    # registered by hand.
    all_types.add("complex64")
    all_types.add("complex128")
    if hasattr(np, "complex192"):
        all_types.add("complex192")
        _isizes.append(24)
    if hasattr(np, "complex256"):
        all_types.add("complex256")
        _isizes.append(32)

    def __init__(
        self, itemsize: int, shape: Shape = (), dflt: Any = _defvalue
    ) -> None:
        if itemsize not in self._isizes:
            raise _invalid_itemsize_error("complex", itemsize, self._isizes)
        self.type = "%s%d" % (self.kind, itemsize * 8)
        Atom.__init__(self, self.type, shape, dflt)


class _ComplexErrorAtom(ComplexAtom, metaclass=type):
    """Reminds the user to stop using the old complex atom names."""

    def __init__(
        self, shape: Shape = (), dflt=ComplexAtom._defvalue
    ) -> NoReturn:
        raise TypeError(
            "to avoid confusions with PyTables 1.X complex atom names, "
            "please use ``ComplexAtom(itemsize=N)``, "
            "where N=8 for single precision complex atoms, "
            "and N=16 for double precision complex atoms"
        )


Complex32Atom = Complex64Atom = Complex128Atom = _ComplexErrorAtom
if hasattr(np, "complex192"):
    Complex192Atom = _ComplexErrorAtom
if hasattr(np, "complex256"):
    Complex256Atom = _ComplexErrorAtom


class TimeAtom(Atom):  # type: ignore[misc]
    """Defines an atom of time type (time kind).

    There are two distinct supported types of time: a 32 bit integer value and
    a 64 bit floating point value. Both of them reflect the number of seconds
    since the Unix epoch. This atom has the property of being stored using the
    HDF5 time datatypes.

    """

    kind: str = "time"
    _deftype = "time32"
    _defvalue: int | float = 0
    __init__ = _abstract_atom_init(
        _deftype, _defvalue
    )  # type: ignore[assignment]


class Time32Atom(TimeAtom):  # type: ignore[misc]
    """Defines an atom of type time32."""

    itemsize: int = 4
    type: str = "time32"  # noqa: A003
    _defvalue = 0

    def __init__(self, shape: Shape = (), dflt=_defvalue) -> None:
        Atom.__init__(self, "int32", shape, dflt)


class Time64Atom(TimeAtom):  # type: ignore[misc]
    """Defines an atom of type time64."""

    itemsize: int = 8
    type: str = "time64"  # noqa: A003
    _defvalue: float = 0.0

    def __init__(self, shape: Shape = (), dflt: float = _defvalue) -> None:
        Atom.__init__(self, "float64", shape, dflt)


class EnumAtom(Atom):
    """Description of an atom of an enumerated type.

    Instances of this class describe the atom type used to store enumerated
    values. Those values belong to an enumerated type, defined by the first
    argument (enum) in the constructor of the atom, which accepts the same
    kinds of arguments as the Enum class (see :ref:`EnumClassDescr`).  The
    enumerated type is stored in the enum attribute of the atom.

    A default value must be specified as the second argument (dflt) in the
    constructor; it must be the *name* (a string) of one of the enumerated
    values in the enumerated type. When the atom is created, the corresponding
    concrete value is broadcast and stored in the dflt attribute (setting
    different default values for items in a multidimensional atom is not
    supported yet). If the name does not match any value in the enumerated
    type, a KeyError is raised.

    Another atom must be specified as the base argument in order to determine
    the base type used for storing the values of enumerated values in memory
    and disk. This *storage atom* is kept in the base attribute of the created
    atom. As a shorthand, you may specify a PyTables type instead of the
    storage atom, implying that this has a scalar shape.

    The storage atom should be able to represent each and every concrete value
    in the enumeration. If it is not, a TypeError is raised. The default value
    of the storage atom is ignored.

    The type attribute of enumerated atoms is always enum.

    Enumerated atoms also support comparisons with other objects::

        >>> enum = ['T0', 'T1', 'T2']
        >>> atom1 = EnumAtom(enum, 'T0', 'int8')  # same as ``atom2``
        >>> atom2 = EnumAtom(enum, 'T0', Int8Atom())  # same as ``atom1``
        >>> atom3 = EnumAtom(enum, 'T0', 'int16')
        >>> atom4 = Int8Atom()
        >>> atom1 == enum
        False
        >>> atom1 == atom2
        True
        >>> atom2 != atom1
        False
        >>> atom1 == atom3
        False
        >>> atom1 == atom4
        False
        >>> atom4 != atom1
        True

    Examples
    --------
    The next C enum construction::

        enum myEnum {
            T0,
            T1,
            T2
        };

    would correspond to the following PyTables
    declaration::

        >>> my_enum_atom = EnumAtom(['T0', 'T1', 'T2'], 'T0', 'int32')

    Please note the dflt argument with a value of 'T0'. Since the concrete
    value matching T0 is unknown right now (we have not used explicit concrete
    values), using the name is the only option left for defining a default
    value for the atom.

    The chosen representation of values for this enumerated atom uses unsigned
    32-bit integers, which surely wastes quite a lot of memory. Another size
    could be selected by using the base argument (this time with a full-blown
    storage atom)::

        >>> my_enum_atom = EnumAtom(['T0', 'T1', 'T2'], 'T0', UInt8Atom())

    You can also define multidimensional arrays for data elements::

        >>> my_enum_atom = EnumAtom(
        ...    ['T0', 'T1', 'T2'], 'T0', base='uint32', shape=(3,2))

    for 3x2 arrays of uint32.

    """

    # Registering this class in the class map may be a little wrong,
    # since the ``Atom.from_kind()`` method fails miserably with
    # enumerations, as they don't support an ``itemsize`` argument.
    # However, resetting ``__metaclass__`` to ``type`` doesn't seem to
    # work and I don't feel like creating a subclass of ``MetaAtom``.

    kind = "enum"
    type = "enum"  # noqa: A003

    @property  # type: ignore[misc]
    def itemsize(self) -> int:  # type: ignore[override]
        """Size in bytes of a single item in the atom."""
        return self.dtype.base.itemsize

    def _checkbase(self, base: Atom) -> None:
        """Check the `base` storage atom."""
        if base.kind == "enum":
            raise TypeError(
                "can not use an enumerated atom "
                "as a storage atom: %r" % base
            )

        # Check whether the storage atom can represent concrete values
        # in the enumeration...
        basedtype = base.dtype
        pyvalues = [value for (name, value) in self.enum]
        try:
            npgenvalues = np.array(pyvalues)
        except ValueError:
            raise TypeError("concrete values are not uniformly-shaped")
        try:
            npvalues = np.array(npgenvalues, dtype=basedtype.base)
        except ValueError:
            raise TypeError(
                "storage atom type is incompatible with "
                "concrete values in the enumeration"
            )
        if npvalues.shape[1:] != basedtype.shape:
            raise TypeError(
                "storage atom shape does not match that of "
                "concrete values in the enumeration"
            )
        if npvalues.tolist() != npgenvalues.tolist():
            raise TypeError(
                "storage atom type lacks precision for "
                "concrete values in the enumeration"
            )

        # ...with some implementation limitations.
        if npvalues.dtype.kind not in ["i", "u"]:
            raise NotImplementedError(
                "only integer concrete values "
                "are supported for the moment, sorry"
            )
        if len(npvalues.shape) > 1:
            raise NotImplementedError(
                "only scalar concrete values "
                "are supported for the moment, sorry"
            )

    def _get_init_args(self) -> dict[str, Any]:
        """Get a dictionary of instance constructor arguments."""
        return {
            "enum": self.enum,
            "dflt": self._defname,
            "base": self.base,
            "shape": self.shape,
        }

    def _is_equal_to_atom(self, atom) -> bool:
        """Return True if the object is equal to the given `atom`."""
        return False

    def _is_equal_to_enumatom(self, enumatom: EnumAtom) -> bool:
        """Return True if the object is equal to the given `enumatom`."""
        return (
            self.enum == enumatom.enum
            and self.shape == enumatom.shape
            and np.all(self.dflt == enumatom.dflt)
            and self.base == enumatom.base
        )

    def __init__(
        self, enum: Enum | Any, dflt: Any, base: Atom | str, shape: Shape = ()
    ) -> None:
        if not isinstance(enum, Enum):
            enum = Enum(enum)
        self.enum = enum

        if isinstance(base, str):
            base = Atom.from_type(base)

        self._checkbase(base)
        self.base = base
        assert isinstance(self.base, Atom)

        default = enum[dflt]  # check default value
        self._defname = dflt  # kept for representation purposes

        # These are kept to ease dumping this particular
        # representation of the enumeration to storage.
        names, values = [], []
        for name, value in enum:
            names.append(name)
            values.append(value)
        basedtype = self.base.dtype

        self._names = names
        self._values = np.array(values, dtype=basedtype.base)

        Atom.__init__(self, basedtype, shape, default)

    def __repr__(self) -> str:
        return "EnumAtom(enum={!r}, dflt={!r}, base={!r}, shape={!r})".format(
            self.enum,
            self._defname,
            self.base,
            self.shape,
        )

    __eq__ = _cmp_dispatcher("_is_equal_to_enumatom")

    # XXX: API incompatible change for PyTables 3 line
    # Overriding __eq__ blocks inheritance of __hash__ in 3.x
    # def __hash__(self):
    #    return hash((self.__class__, self.enum, self.shape, self.dflt,
    #                 self.base))


class ReferenceAtom(Atom):
    """Defines an atom of type object to read references.

    This atom is read-only.
    """

    kind = "reference"
    type = "object"  # noqa: A003
    _deftype = "NoneType"
    _defvalue = None

    @property  # type: ignore[misc]
    def itemsize(self) -> int:  # type: ignore[override]
        """Size in bytes of a single item in the atom."""
        return self.dtype.base.itemsize

    def __init__(self, shape: Shape = ()) -> None:
        Atom.__init__(self, self.type, shape, self._defvalue)

    def __repr__(self) -> str:
        return f"ReferenceAtom(shape={self.shape})"


# Pseudo-atom classes
# ===================
#
# Now, there come three special classes, `ObjectAtom`, `VLStringAtom`
# and `VLUnicodeAtom`, that actually do not descend from `Atom`, but
# which goal is so similar that they should be described here.
# Pseudo-atoms can only be used with `VLArray` datasets, and they do
# not support multidimensional values, nor multiple values per row.
#
# They can be recognised because they also have ``kind``, ``type`` and
# ``shape`` attributes, but no ``size``, ``itemsize`` or ``dflt``
# ones.  Instead, they have a ``base`` atom which defines the elements
# used for storage.
#
# See ``examples/vlarray1.py`` and ``examples/vlarray2.py`` for
# further examples on `VLArray` datasets, including object
# serialization and string management.


class PseudoAtom:
    """Pseudo-atoms can only be used in ``VLArray`` nodes.

    They can be recognised because they also have `kind`, `type` and
    `shape` attributes, but no `size`, `itemsize` or `dflt` ones.
    Instead, they have a `base` atom which defines the elements used
    for storage.
    """

    base: Atom

    def __repr__(self) -> str:
        return "%s()" % self.__class__.__name__

    def toarray(self, object_: Any) -> NoReturn:
        """Convert an `object_` into an array of base atoms."""
        raise NotImplementedError

    def fromarray(self, array: Any) -> NoReturn:
        """Convert an `array` of base atoms into an object."""
        raise NotImplementedError


class _BufferedAtom(PseudoAtom):
    """Pseudo-atom which stores data as a buffer (flat array of uints)."""

    shape = ()

    def toarray(self, object_: Any) -> np.ndarray:
        buffer_ = self._tobuffer(object_)
        array = np.ndarray(
            buffer=buffer_, dtype=self.base.dtype, shape=len(buffer_)
        )
        return array

    def _tobuffer(self, object_: Any) -> NoReturn:
        """Convert an `object_` into a buffer."""
        raise NotImplementedError


class VLStringAtom(_BufferedAtom):
    """Defines an atom of type ``vlstring``.

    This class describes a *row* of the VLArray class, rather than an atom. It
    differs from the StringAtom class in that you can only add *one instance of
    it to one specific row*, i.e. the :meth:`VLArray.append` method only
    accepts one object when the base atom is of this type.

    This class stores bytestrings. It does not make assumptions on the
    encoding of the string, and raw bytes are stored as is. To store a string
    you will need to *explicitly* convert it to a bytestring before you can
    save them::

        >>> s = 'A unicode string: hbar = \u210f'
        >>> bytestring = s.encode('utf-8')
        >>> VLArray.append(bytestring) # noqa: F821  # doctest: +SKIP

    For full Unicode support, using VLUnicodeAtom (see :ref:`VLUnicodeAtom`) is
    recommended.

    Variable-length string atoms do not accept parameters and they cause the
    reads of rows to always return Python bytestrings.  You can regard vlstring
    atoms as an easy way to save generic variable length strings.

    """

    kind = "vlstring"
    type = "vlstring"  # noqa: A003
    base = UInt8Atom()

    def _tobuffer(self, object_: bytes) -> np.bytes_:
        if not isinstance(object_, bytes):
            raise TypeError(f"object is not bytes: {object_!r}")
        return np.bytes_(object_)

    def fromarray(self, array: np.ndarray) -> bytes:
        """Convert array data into bytes."""
        return array.tobytes()


class VLUnicodeAtom(_BufferedAtom):
    """Defines an atom of type vlunicode.

    This class describes a *row* of the VLArray class, rather than an atom.  It
    is very similar to VLStringAtom (see :ref:`VLStringAtom`), but it stores
    Unicode strings (using 32-bit characters a la UCS-4, so all strings of the
    same length also take up the same space).

    This class does not make assumptions on the encoding of plain input
    strings.  Plain strings are supported as long as no character is out of the
    ASCII set; otherwise, you will need to *explicitly* convert them to Unicode
    before you can save them.

    Variable-length Unicode atoms do not accept parameters and they cause the
    reads of rows to always return Python Unicode strings.  You can regard
    vlunicode atoms as an easy way to save variable length Unicode strings.

    """

    kind = "vlunicode"
    type = "vlunicode"  # noqa: A003
    base = UInt32Atom()

    # numpy.unicode_ no more implements the buffer interface in Python 3
    #
    # When the Python build is UCS-2, we need to promote the
    # Unicode string to UCS-4.  We *must* use a 0-d array since
    # NumPy scalars inherit the UCS-2 encoding from Python (see
    # NumPy ticket #525).  Since ``_tobuffer()`` can't return an
    # array, we must override ``toarray()`` itself.
    def toarray(self, object_: str) -> np.ndarray:
        """Convert a string into a numpy array."""
        if not isinstance(object_, str):
            raise TypeError(f"object is not a string: {object_!r}")
        ustr = str(object_)
        uarr = np.array(ustr, dtype="U")
        return np.ndarray(buffer=uarr, dtype=self.base.dtype, shape=len(ustr))

    def _tobuffer(self, object_: str) -> np.str_:
        # This works (and is used) only with UCS-4 builds of Python,
        # where the width of the internal representation of a
        # character matches that of the base atoms.
        if not isinstance(object_, str):
            raise TypeError(f"object is not a string: {object_!r}")
        return np.str_(object_)

    def fromarray(self, array: np.ndarray) -> str:
        """Convert array data into a string."""
        length = len(array)
        if length == 0:
            return ""  # ``array.view('U0')`` raises a `TypeError`
        return array.view("U%d" % length).item()


class ObjectAtom(_BufferedAtom):
    """Defines an atom of type object.

    This class is meant to fit *any* kind of Python object in a row of a
    VLArray dataset by using pickle behind the scenes. Due to the fact that
    you can not foresee how long will be the output of the pickle
    serialization (i.e. the atom already has a *variable* length), you can only
    fit *one object per row*. However, you can still group several objects in a
    single tuple or list and pass it to the :meth:`VLArray.append` method.

    Object atoms do not accept parameters and they cause the reads of rows to
    always return Python objects. You can regard object atoms as an easy way to
    save an arbitrary number of generic Python objects in a VLArray dataset.

    """

    kind = "object"
    type = "object"  # noqa: A003
    base = UInt8Atom()

    def _tobuffer(self, object_: object) -> bytes:
        return pickle.dumps(object_, pickle.HIGHEST_PROTOCOL)

    def fromarray(self, array: np.ndarray) -> Any | None:
        """Deserialize data contained in the input array.

        A Python object is returned.
        """
        # We have to check for an empty array because of a possible
        # bug in HDF5 which makes it claim that a dataset has one
        # record when in fact it is empty.
        if array.size == 0:
            return None
        return pickle.loads(array.tobytes())
