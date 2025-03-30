"""Classes for describing columns for ``Table`` objects."""

from __future__ import annotations

import copy
import warnings
from typing import Any, Literal
from collections.abc import Callable, Generator, Sequence

import numpy as np
import numpy.typing as npt

from . import atom
from .path import check_name_validity

__docformat__ = "reStructuredText"
"""The format of documentation strings in this module."""


def same_position(
    oldmethod: Callable[[Col, Col], bool]
) -> Callable[[Col, Col], bool]:
    """Decorate `oldmethod` to also compare the `_v_pos` attribute."""

    def newmethod(self: Col, other: Col) -> bool:
        try:
            other._v_pos
        except AttributeError:
            return False  # not a column definition
        return self._v_pos == other._v_pos and oldmethod(self, other)

    newmethod.__name__ = oldmethod.__name__
    newmethod.__doc__ = oldmethod.__doc__
    return newmethod


class Col(atom.Atom, metaclass=type):
    """Defines a non-nested column.

    Col instances are used as a means to declare the different properties of a
    non-nested column in a table or nested column.  Col classes are descendants
    of their equivalent Atom classes (see :ref:`AtomClassDescr`), but their
    instances have an additional _v_pos attribute that is used to decide the
    position of the column inside its parent table or nested column (see the
    IsDescription class in :ref:`IsDescriptionClassDescr` for more information
    on column positions).

    In the same fashion as Atom, you should use a particular Col descendant
    class whenever you know the exact type you will need when writing your
    code. Otherwise, you may use one of the Col.from_*() factory methods.

    Each factory method inherited from the Atom class is available with the
    same signature, plus an additional pos parameter (placed in last position)
    which defaults to None and that may take an integer value.  This parameter
    might be used to specify the position of the column in the table.

    Besides, there are the next additional factory methods, available only for
    Col objects.

    The following parameters are available for most Col-derived constructors.

    Parameters
    ----------
    itemsize : int
        For types with a non-fixed size, this sets the size in bytes of
        individual items in the column.
    shape : tuple
        Sets the shape of the column. An integer shape of N is equivalent to
        the tuple (N,).
    dflt
        Sets the default value for the column.
    pos : int
        Sets the position of column in table.  If unspecified, the position
        will be randomly selected.
    attrs : dict
        Attribute metadata stored in the column (see
        :ref:`AttributeSetClassDescr`).

    """

    # filled as column classes are created
    _class_from_prefix: dict[str, type[Col]] = {}
    """Maps column prefixes to column classes."""

    @classmethod
    def prefix(cls) -> str:
        """Return the column class prefix."""
        cname = cls.__name__
        return cname[: cname.rfind("Col")]

    @classmethod
    def from_atom(
        cls,
        atom: atom.Atom,
        pos: int | None = None,
        _offset: int | None = None,
    ) -> Col:
        """Create a Col definition from a PyTables atom.

        An optional position may be specified as the pos argument.

        """
        prefix = atom.prefix()
        kwargs = atom._get_init_args()
        colclass = cls._class_from_prefix[prefix]
        return colclass(pos=pos, _offset=_offset, **kwargs)

    @classmethod
    def from_sctype(
        cls,
        sctype: str | np.dtype,
        shape: tuple[int, ...] = (),
        dflt: Any | None = None,
        pos: int | None = None,
    ) -> Col:
        """Create a `Col` definition from a NumPy scalar type `sctype`.

        Optional shape, default value and position may be specified as
        the `shape`, `dflt` and `pos` arguments, respectively.
        Information in the `sctype` not represented in a `Col` is
        ignored.

        """
        newatom = atom.Atom.from_sctype(sctype, shape, dflt)
        return cls.from_atom(newatom, pos=pos)

    @classmethod
    def from_dtype(
        cls,
        dtype,
        dflt: Any | None = None,
        pos: int | None = None,
        _offset: int | None = None,
    ) -> Col:
        """Create a `Col` definition from a NumPy `dtype`.

        Optional default value and position may be specified as the
        `dflt` and `pos` arguments, respectively.  The `dtype` must have
        a byte order which is irrelevant or compatible with that of the
        system.  Information in the `dtype` not represented in a `Col`
        is ignored.

        """
        newatom = atom.Atom.from_dtype(dtype, dflt)
        return cls.from_atom(newatom, pos=pos, _offset=_offset)

    @classmethod
    def from_type(
        cls,
        type_,
        shape: tuple[int, ...] = (),
        dflt: Any | None = None,
        pos: int | None = None,
    ) -> Col:
        """Create a `Col` definition from a PyTables `type`.

        Optional shape, default value and position may be specified as
        the `shape`, `dflt` and `pos` arguments, respectively.

        """
        newatom = atom.Atom.from_type(type_, shape, dflt)
        return cls.from_atom(newatom, pos=pos)

    @classmethod
    def from_kind(
        cls,
        kind,
        itemsize=None,
        shape: tuple[int, ...] = (),
        dflt: Any | None = None,
        pos: int | None = None,
    ) -> Col:
        """Create a `Col` definition from a PyTables `kind`.

        Optional item size, shape, default value and position may be
        specified as the `itemsize`, `shape`, `dflt` and `pos`
        arguments, respectively.  Bear in mind that not all columns
        support a default item size.

        """
        newatom = atom.Atom.from_kind(kind, itemsize, shape, dflt)
        return cls.from_atom(newatom, pos=pos)

    @classmethod
    def _subclass_from_prefix(cls, prefix: str) -> type[Col]:
        """Get a column subclass for the given `prefix`."""
        cname = "%sCol" % prefix
        class_from_prefix = cls._class_from_prefix
        if cname in class_from_prefix:
            return class_from_prefix[cname]
        atombase = getattr(atom, "%sAtom" % prefix)

        class NewCol(cls, atombase):
            """Defines a non-nested column of a particular type.

            The constructor accepts the same arguments as the equivalent
            `Atom` class, plus an additional ``pos`` argument for
            position information, which is assigned to the `_v_pos`
            attribute and an ``attrs`` argument for storing additional metadata
            similar to `table.attrs`, which is assigned to the `_v_col_attrs`
            attribute.

            """

            def __init__(self, *args, **kwargs) -> None:
                pos = kwargs.pop("pos", None)
                col_attrs = kwargs.pop("attrs", {})
                offset = kwargs.pop("_offset", None)
                class_from_prefix = self._class_from_prefix
                atombase.__init__(self, *args, **kwargs)
                # The constructor of an abstract atom may have changed
                # the class of `self` to something different of `NewCol`
                # and `atombase` (that's why the prefix map is saved).
                if self.__class__ is not NewCol:
                    colclass = class_from_prefix[self.prefix()]
                    self.__class__ = colclass
                self._v_pos = pos
                self._v_offset = offset
                self._v_col_attrs = col_attrs

            __eq__ = same_position(atombase.__eq__)
            _is_equal_to_atom = same_position(atombase._is_equal_to_atom)

            # XXX: API incompatible change for PyTables 3 line
            # Overriding __eq__ blocks inheritance of __hash__ in 3.x
            # def __hash__(self):
            #    return hash((self._v_pos, self.atombase))

            if prefix == "Enum":
                _is_equal_to_enumatom = same_position(
                    atombase._is_equal_to_enumatom
                )

        NewCol.__name__ = cname

        class_from_prefix[prefix] = NewCol
        return NewCol

    def __repr__(self) -> str:
        # Reuse the atom representation.
        atomrepr = super().__repr__()
        lpar = atomrepr.index("(")
        rpar = atomrepr.rindex(")")
        atomargs = atomrepr[lpar + 1 : rpar]
        classname = self.__class__.__name__
        if self._v_col_attrs:
            return (
                f"{classname}({atomargs}, pos={self._v_pos}"
                f", attrs={self._v_col_attrs})"
            )
        return f"{classname}({atomargs}, pos={self._v_pos})"

    def _get_init_args(self) -> dict[str, Any]:
        """Get a dictionary of instance constructor arguments."""
        kwargs = {arg: getattr(self, arg) for arg in ("shape", "dflt")}
        kwargs["pos"] = getattr(self, "_v_pos", None)
        return kwargs


def _generate_col_classes() -> Generator[type[Col]]:
    """Generate all column classes."""
    # Abstract classes are not in the class map.
    cprefixes = ["Int", "UInt", "Float", "Time"]
    for kind, kdata in atom.atom_map.items():
        if hasattr(kdata, "kind"):  # atom class: non-fixed item size
            atomclass = kdata
            cprefixes.append(atomclass.prefix())
        else:  # dictionary: fixed item size
            for atomclass in kdata.values():
                cprefixes.append(atomclass.prefix())

    # Bottom-level complex classes are not in the type map, of course.
    # We still want the user to get the compatibility warning, though.
    cprefixes.extend(["Complex32", "Complex64", "Complex128"])
    if hasattr(atom, "Complex192Atom"):
        cprefixes.append("Complex192")
    if hasattr(atom, "Complex256Atom"):
        cprefixes.append("Complex256")

    for cprefix in cprefixes:
        newclass = Col._subclass_from_prefix(cprefix)
        yield newclass


# Create all column classes.
# for _newclass in _generate_col_classes():
#     exec('%s = _newclass' % _newclass.__name__)
# del _newclass

StringCol = Col._subclass_from_prefix("String")
BoolCol = Col._subclass_from_prefix("Bool")
EnumCol = Col._subclass_from_prefix("Enum")
IntCol = Col._subclass_from_prefix("Int")
Int8Col = Col._subclass_from_prefix("Int8")
Int16Col = Col._subclass_from_prefix("Int16")
Int32Col = Col._subclass_from_prefix("Int32")
Int64Col = Col._subclass_from_prefix("Int64")
UIntCol = Col._subclass_from_prefix("UInt")
UInt8Col = Col._subclass_from_prefix("UInt8")
UInt16Col = Col._subclass_from_prefix("UInt16")
UInt32Col = Col._subclass_from_prefix("UInt32")
UInt64Col = Col._subclass_from_prefix("UInt64")

FloatCol = Col._subclass_from_prefix("Float")
if hasattr(atom, "Float16Atom"):
    Float16Col = Col._subclass_from_prefix("Float16")
Float32Col = Col._subclass_from_prefix("Float32")
Float64Col = Col._subclass_from_prefix("Float64")
if hasattr(atom, "Float96Atom"):
    Float96Col = Col._subclass_from_prefix("Float96")
if hasattr(atom, "Float128Atom"):
    Float128Col = Col._subclass_from_prefix("Float128")

ComplexCol = Col._subclass_from_prefix("Complex")
Complex32Col = Col._subclass_from_prefix("Complex32")
Complex64Col = Col._subclass_from_prefix("Complex64")
Complex128Col = Col._subclass_from_prefix("Complex128")
if hasattr(atom, "Complex192Atom"):
    Complex192Col = Col._subclass_from_prefix("Complex192")
if hasattr(atom, "Complex256Atom"):
    Complex256Col = Col._subclass_from_prefix("Complex256")

TimeCol = Col._subclass_from_prefix("Time")
Time32Col = Col._subclass_from_prefix("Time32")
Time64Col = Col._subclass_from_prefix("Time64")


# Table description classes
# =========================
class Description:
    """This class represents descriptions of the structure of tables.

    An instance of this class is automatically bound to Table (see
    :ref:`TableClassDescr`) objects when they are created.  It provides a
    browseable representation of the structure of the table, made of non-nested
    (Col - see :ref:`ColClassDescr`) and nested (Description) columns.

    Column definitions under a description can be accessed as attributes of it
    (*natural naming*). For instance, if table.description is a Description
    instance with a column named col1 under it, the later can be accessed as
    table.description.col1. If col1 is nested and contains a col2 column, this
    can be accessed as table.description.col1.col2. Because of natural naming,
    the names of members start with special prefixes, like in the Group class
    (see :ref:`GroupClassDescr`).

    .. rubric:: Description attributes

    .. attribute:: _v_colobjects

        A dictionary mapping the names of the columns hanging
        directly from the associated table or nested column to their
        respective descriptions (Col - see :ref:`ColClassDescr` or
        Description - see :ref:`DescriptionClassDescr` instances).

        .. versionchanged:: 3.0
           The *_v_colObjects* attribute has been renamed into
           *_v_colobjects*.

    .. attribute:: _v_dflts

        A dictionary mapping the names of non-nested columns
        hanging directly from the associated table or nested column
        to their respective default values.

    .. attribute:: _v_dtype

        The NumPy type which reflects the structure of this
        table or nested column.  You can use this as the
        dtype argument of NumPy array factories.

    .. attribute:: _v_dtypes

        A dictionary mapping the names of non-nested columns
        hanging directly from the associated table or nested column
        to their respective NumPy types.

    .. attribute:: _v_is_nested

        Whether the associated table or nested column contains
        further nested columns or not.

    .. attribute:: _v_itemsize

        The size in bytes of an item in this table or nested column.

    .. attribute:: _v_name

        The name of this description group. The name of the
        root group is '/'.

    .. attribute:: _v_names

        A list of the names of the columns hanging directly
        from the associated table or nested column. The order of the
        names matches the order of their respective columns in the
        containing table.

    .. attribute:: _v_nested_descr

        A nested list of pairs of (name, format) tuples for all the columns
        under this table or nested column. You can use this as the dtype and
        descr arguments of NumPy array factories.

        .. versionchanged:: 3.0
           The *_v_nestedDescr* attribute has been renamed into
           *_v_nested_descr*.

    .. attribute:: _v_nested_formats

        A nested list of the NumPy string formats (and shapes) of all the
        columns under this table or nested column. You can use this as the
        formats argument of NumPy array factories.

        .. versionchanged:: 3.0
           The *_v_nestedFormats* attribute has been renamed into
           *_v_nested_formats*.

    .. attribute:: _v_nestedlvl

        The level of the associated table or nested column in the nested
        datatype.

    .. attribute:: _v_nested_names

        A nested list of the names of all the columns under this table or
        nested column. You can use this as the names argument of NumPy array
        factories.

        .. versionchanged:: 3.0
           The *_v_nestedNames* attribute has been renamed into
           *_v_nested_names*.

    .. attribute:: _v_pathname

        Pathname of the table or nested column.

    .. attribute:: _v_pathnames

        A list of the pathnames of all the columns under this table or nested
        column (in preorder).  If it does not contain nested columns, this is
        exactly the same as the :attr:`Description._v_names` attribute.

    .. attribute:: _v_types

        A dictionary mapping the names of non-nested columns hanging directly
        from the associated table or nested column to their respective PyTables
        types.

    .. attribute:: _v_offsets

        A list of offsets for all the columns.  If the list is empty, means
        that there are no padding in the data structure.  However, the support
        for offsets is currently limited to flat tables; for nested tables, the
        potential padding is always removed (exactly the same as in pre-3.5
        versions), and this variable is set to empty.

        .. versionadded:: 3.5
           Previous to this version all the compound types were converted
           internally to 'packed' types, i.e. with no padding between the
           component types.  Starting with 3.5, the holes in native HDF5
           types (non-nested) are honored and replicated during dataset
           and attribute copies.
    """

    def __init__(
        self,
        classdict: dict[str, Any],
        nestedlvl: int = -1,
        validate: bool = True,
        ptparams: dict[str, Any] | None = None,
    ) -> None:

        if not classdict:
            raise ValueError("cannot create an empty data type")

        # Do a shallow copy of classdict just in case this is going to
        # be shared by other instances
        newdict = self.__dict__
        newdict["_v_name"] = "/"  # The name for root descriptor
        newdict["_v_names"] = []
        newdict["_v_dtypes"] = {}
        newdict["_v_types"] = {}
        newdict["_v_dflts"] = {}
        newdict["_v_colobjects"] = {}
        newdict["_v_is_nested"] = False
        nested_formats = []
        nested_dtype = []

        if not hasattr(newdict, "_v_nestedlvl"):
            newdict["_v_nestedlvl"] = nestedlvl + 1

        cols_with_pos = []  # colum (position, name) pairs
        cols_no_pos = []  # just column names
        cols_offsets = []  # the offsets of the columns
        valid_offsets = False  # by default there a no valid offsets

        # Check for special variables and convert column descriptions
        for name, descr in classdict.items():
            if name.startswith("_v_"):
                if name in newdict:
                    # print("Warning!")
                    # special methods &c: copy to newdict, warn about conflicts
                    warnings.warn(
                        f"Can't set attr {name!r} in description "
                        f"class {self!r}"
                    )
                else:
                    # print("Special variable!-->", name, classdict[name])
                    newdict[name] = descr
                continue  # This variable is not needed anymore

            columns = None
            if type(descr) is type(IsDescription) and issubclass(
                descr, IsDescription
            ):
                # print("Nested object (type I)-->", name)
                columns = descr().columns
            elif type(descr.__class__) is type(IsDescription) and issubclass(
                descr.__class__, IsDescription
            ):
                # print("Nested object (type II)-->", name)
                columns = descr.columns
            elif isinstance(descr, dict):
                # print("Nested object (type III)-->", name)
                columns = descr
            else:
                # print("Nested object (type IV)-->", name)
                descr = copy.copy(descr)
            # The copies above and below ensure that the structures
            # provided by the user will remain unchanged even if we
            # tamper with the values of ``_v_pos`` here.
            if columns is not None:
                descr = Description(
                    copy.copy(columns), self._v_nestedlvl, ptparams=ptparams
                )
            classdict[name] = descr

            pos = getattr(descr, "_v_pos", None)
            if pos is None:
                cols_no_pos.append(name)
            else:
                cols_with_pos.append((pos, name))
            offset = getattr(descr, "_v_offset", None)
            if offset is not None:
                cols_offsets.append(offset)

        # Sort field names:
        #
        # 1. Fields with explicit positions, according to their
        #    positions (and their names if coincident).
        # 2. Fields with no position, in alphabetical order.
        cols_with_pos.sort()
        cols_no_pos.sort()
        keys = [name for (pos, name) in cols_with_pos] + cols_no_pos

        pos = 0
        nested = False
        # Get properties for compound types
        for k in keys:
            if validate:
                # Check for key name validity
                check_name_validity(k)
            # Class variables
            obj = classdict[k]
            newdict[k] = obj  # To allow natural naming
            if not isinstance(obj, (Col, Description)):
                raise TypeError(
                    f"Passing an incorrect value to a table column."
                    f" Expected a Col (or subclass) instance and "
                    f'got: "{obj}". Please make use of the Col(), or '
                    f"descendant, constructor to properly "
                    f"initialize columns."
                )
            obj._v_pos = pos  # Set the position of this object
            obj._v_parent = self  # The parent description
            pos += 1
            newdict["_v_colobjects"][k] = obj
            newdict["_v_names"].append(k)
            obj.__dict__["_v_name"] = k

            if not isinstance(k, str):
                # numpy only accepts "str" for field names
                # Python 3.x: bytes --> str (unicode)
                kk = k.decode()
            else:
                kk = k

            if isinstance(obj, Col):
                dtype = obj.dtype
                newdict["_v_dtypes"][k] = dtype
                newdict["_v_types"][k] = obj.type
                newdict["_v_dflts"][k] = obj.dflt
                nested_formats.append(obj.recarrtype)
                baserecarrtype = dtype.base.str[1:]
                nested_dtype.append((kk, baserecarrtype, dtype.shape))
            else:  # A description
                nested_formats.append(obj._v_nested_formats)
                nested_dtype.append((kk, obj._v_dtype))
                nested = True

        # Useful for debugging purposes
        # import traceback
        # if ptparams is None:
        #     print("*** print_stack:")
        #     traceback.print_stack()

        # Check whether we are gonna use padding or not.  Two possibilities:
        #   1) Make padding True by default (except if ALLOW_PADDING is set
        #      to False)
        #   2) Make padding False by default (except if ALLOW_PADDING is set
        #       to True)
        # Currently we choose 1) because it favours honoring padding even on
        # unhandled situations (should be very few).
        # However, for development, option 2) is recommended as it catches
        # most of the unhandled situations.
        allow_padding = ptparams is None or ptparams["ALLOW_PADDING"]
        # allow_padding = ptparams is not None and ptparams['ALLOW_PADDING']
        if (
            allow_padding
            and len(cols_offsets) > 1
            and len(keys) == len(cols_with_pos)
            and len(keys) == len(cols_offsets)
            and not nested
        ):  # TODO: support offsets with nested types
            # We have to sort the offsets too, as they must follow the column
            # order. As the offsets and the pos should be place in the same
            # order, a single sort is enough here.
            cols_offsets.sort()
            valid_offsets = True
        else:
            newdict["_v_offsets"] = []

        # Assign the format list to _v_nested_formats
        newdict["_v_nested_formats"] = nested_formats

        if self._v_nestedlvl == 0:
            # Get recursively nested _v_nested_names and _v_nested_descr attrs
            self._g_set_nested_names_descr()
            # Get pathnames for nested groups
            self._g_set_path_names()
            # Check the _v_byteorder has been used an issue an Error
            if hasattr(self, "_v_byteorder"):
                raise ValueError(
                    "Using a ``_v_byteorder`` in the description is obsolete. "
                    "Use the byteorder parameter in the constructor instead."
                )

        # Compute the dtype with offsets or without
        # print("offsets ->", cols_offsets, nestedDType, nested, valid_offsets)
        if valid_offsets:
            # TODO: support offsets within nested types
            dtype_fields = {
                "names": newdict["_v_names"],
                "formats": nested_formats,
                "offsets": cols_offsets,
            }
            itemsize = newdict.get("_v_itemsize", None)
            if itemsize is not None:
                dtype_fields["itemsize"] = itemsize
            dtype = np.dtype(dtype_fields)
        else:
            dtype = np.dtype(nested_dtype)
        newdict["_v_dtype"] = dtype
        newdict["_v_itemsize"] = dtype.itemsize
        newdict["_v_offsets"] = [dtype.fields[name][1] for name in dtype.names]

    def _g_set_nested_names_descr(self) -> None:
        """Compute the nested names and descriptions for nested datatypes."""
        names = self._v_names
        fmts = self._v_nested_formats
        self._v_nested_names = names[:]  # Important to do a copy!
        self._v_nested_descr = list(zip(names, fmts))
        for i, name in enumerate(names):
            new_object = self._v_colobjects[name]
            if isinstance(new_object, Description):
                new_object._g_set_nested_names_descr()
                # replace the column nested name by a correct tuple
                self._v_nested_names[i] = (name, new_object._v_nested_names)
                self._v_nested_descr[i] = (name, new_object._v_nested_descr)
                # set the _v_is_nested flag
                self._v_is_nested = True

    def _g_set_path_names(self) -> None:
        """Compute the pathnames for arbitrary nested descriptions.

        This method sets the ``_v_pathname`` and ``_v_pathnames``
        attributes of all the elements (both descriptions and columns)
        in this nested description.

        """

        def get_cols_in_order(description: Description) -> list[Col]:
            return [
                description._v_colobjects[colname]
                for colname in description._v_names
            ]

        def join_paths(path1: str, path2: str) -> str:
            if not path1:
                return path2
            return f"{path1}/{path2}"

        # The top of the stack always has a nested description
        # and a list of its child columns
        # (be they nested ``Description`` or non-nested ``Col`` objects).
        # In the end, the list contains only a list of column paths
        # under this one.
        #
        # For instance, given this top of the stack::
        #
        #   (<Description X>, [<Column A>, <Column B>])
        #
        # After computing the rest of the stack, the top is::
        #
        #   (<Description X>, ['a', 'a/m', 'a/n', ... , 'b', ...])

        stack: list[tuple[Description, list[Col]]] = []

        # We start by pushing the top-level description
        # and its child columns.
        self._v_pathname = ""
        stack.append((self, get_cols_in_order(self)))

        while stack:
            desc, cols = stack.pop()
            head = cols[0]

            # What's the first child in the list?
            if isinstance(head, Description):
                # A nested description.  We remove it from the list and
                # push it with its child columns.  This will be the next
                # handled description.
                head._v_pathname = join_paths(desc._v_pathname, head._v_name)
                stack.append((desc, cols[1:]))  # alter the top
                stack.append((head, get_cols_in_order(head)))  # new top
            elif isinstance(head, Col):
                # A non-nested column.  We simply remove it from the
                # list and append its name to it.
                head._v_pathname = join_paths(desc._v_pathname, head._v_name)
                cols.append(head._v_name)  # alter the top
                stack.append((desc, cols[1:]))  # alter the top
            else:
                # Since paths and names are appended *to the end* of
                # children lists, a string signals that no more children
                # remain to be processed, so we are done with the
                # description at the top of the stack.
                assert isinstance(head, str)
                # Assign the computed set of descendent column paths.
                desc._v_pathnames = cols
                if len(stack) > 0:
                    # Compute the paths with respect to the parent node
                    # (including the path of the current description)
                    # and append them to its list.
                    desc_name = desc._v_name
                    col_paths = [join_paths(desc_name, path) for path in cols]
                    col_paths.insert(0, desc_name)
                    parent_cols = stack[-1][1]
                    parent_cols.extend(col_paths)
                # (Nothing is pushed, we are done with this description.)

    def _f_walk(
        self,
        type: Literal["All", "Col", "Description"] = "All",  # noqa: A002
    ) -> Generator[Col | Description]:
        """Iterate over nested columns.

        If type is 'All' (the default), all column description objects (Col and
        Description instances) are yielded in top-to-bottom order (preorder).

        If type is 'Col' or 'Description', only column descriptions of that
        type are yielded.

        """
        if type not in ["All", "Col", "Description"]:
            raise ValueError(
                "type can only take the parameters 'All', 'Col' or "
                "'Description'."
            )

        stack: list[Description] = [self]
        while stack:
            obj = stack.pop(0)  # pop at the front so as to ensure the order
            if type in ["All", "Description"]:
                yield obj  # yield description
            for name in obj._v_names:
                new_object = obj._v_colobjects[name]
                if isinstance(new_object, Description):
                    stack.append(new_object)
                else:
                    if type in ["All", "Col"]:
                        yield new_object  # yield column

    def __repr__(self) -> str:
        """Give a detailed Description column representation."""
        rep = [
            f'{"  " * self._v_nestedlvl}"{k}": {self._v_colobjects[k]!r}'
            for k in self._v_names
        ]
        return "{\n  %s}" % (",\n  ".join(rep))

    def __str__(self) -> str:
        """Give a brief Description representation."""
        return f"Description({self._v_nested_descr})"


class MetaIsDescription(type):
    """Helper metaclass to return the class variables as a dictionary."""

    def __new__(
        cls, classname: str, bases: Sequence, classdict: dict[str, Any]
    ) -> MetaIsDescription:
        """Return a new class with a "columns" attribute filled."""
        newdict = {
            "columns": {},
        }
        if "__doc__" in classdict:
            newdict["__doc__"] = classdict["__doc__"]
        for b in bases:
            if "columns" in b.__dict__:
                newdict["columns"].update(b.__dict__["columns"])
        for k in classdict:
            # if not (k.startswith('__') or k.startswith('_v_')):
            # We let pass _v_ variables to configure class behaviour
            if not (k.startswith("__")):
                newdict["columns"][k] = classdict[k]

        # Return a new class with the "columns" attribute filled
        return type.__new__(cls, classname, bases, newdict)


class IsDescription(metaclass=MetaIsDescription):
    """Description of the structure of a table or nested column.

    This class is designed to be used as an easy, yet meaningful way to
    describe the structure of new Table (see :ref:`TableClassDescr`) datasets
    or nested columns through the definition of *derived classes*. In order to
    define such a class, you must declare it as descendant of IsDescription,
    with as many attributes as columns you want in your table. The name of each
    attribute will become the name of a column, and its value will hold a
    description of it.

    Ordinary columns can be described using instances of the Col class (see
    :ref:`ColClassDescr`). Nested columns can be described by using classes
    derived from IsDescription, instances of it, or name-description
    dictionaries. Derived classes can be declared in place (in which case the
    column takes the name of the class) or referenced by name.

    Nested columns can have a _v_pos special attribute which sets the
    *relative* position of the column among sibling columns *also having
    explicit positions*.  The pos constructor argument of Col instances is used
    for the same purpose.  Columns with no explicit position will be placed
    afterwards in alphanumeric order.

    Once you have created a description object, you can pass it to the Table
    constructor, where all the information it contains will be used to define
    the table structure.

    .. rubric:: IsDescription attributes

    .. attribute:: _v_pos

        Sets the position of a possible nested column description among its
        sibling columns.  This attribute can be specified *when declaring*
        an IsDescription subclass to complement its *metadata*.

    .. attribute:: columns

        Maps the name of each column in the description to its own descriptive
        object. This attribute is *automatically created* when an IsDescription
        subclass is declared.  Please note that declared columns can no longer
        be accessed as normal class variables after its creation.

    """


def descr_from_dtype(
    dtype_: npt.DTypeLike, ptparams: dict[str, Any] | None = None
) -> tuple[Description, str]:
    """Get a description instance and byteorder from a (nested) NumPy dtype."""
    fields = {}
    fbyteorder = "|"
    for name in dtype_.names:
        dtype, offset = dtype_.fields[name][:2]
        kind = dtype.base.kind
        byteorder = dtype.base.byteorder
        if byteorder in "><=":
            if fbyteorder not in ["|", byteorder]:
                raise NotImplementedError(
                    "structured arrays with mixed byteorders "
                    "are not supported yet, sorry"
                )
            fbyteorder = byteorder
        # Non-nested column
        if kind in "biufSUc":
            col = Col.from_dtype(dtype, pos=offset, _offset=offset)
        # Nested column
        elif kind == "V" and dtype.shape in [(), (1,)]:
            if dtype.shape != ():
                warnings.warn(
                    "nested descriptions will be converted to scalar"
                )
            col, _ = descr_from_dtype(dtype.base, ptparams=ptparams)
            col._v_pos = offset
            col._v_offset = offset
        else:
            raise NotImplementedError(
                "structured arrays with columns with type description ``%s`` "
                "are not supported yet, sorry" % dtype
            )
        fields[name] = col

    return Description(fields, ptparams=ptparams), fbyteorder


def dtype_from_descr(
    descr: dict | type[IsDescription] | IsDescription,
    byteorder: str | None = None,
    ptparams: dict[str, Any] | None = None,
) -> np.dtype:
    """Get a (nested) NumPy dtype from a description instance and byteorder.

    The descr parameter can be a Description or IsDescription
    instance, sub-class of IsDescription or a dictionary.

    """
    if isinstance(descr, dict):
        descr = Description(descr, ptparams=ptparams)
    elif type(descr) is type(IsDescription) and issubclass(
        descr, IsDescription
    ):
        descr = Description(descr().columns, ptparams=ptparams)
    elif isinstance(descr, IsDescription):
        descr = Description(descr.columns, ptparams=ptparams)
    elif not isinstance(descr, Description):
        raise ValueError(f"invalid description: {descr!r}")

    dtype_ = descr._v_dtype

    if byteorder and byteorder != "|":
        dtype_ = dtype_.newbyteorder(byteorder)

    return dtype_


if __name__ == "__main__":
    """Test code."""

    class Info(IsDescription):  # noqa: D101
        _v_pos = 2
        Name = UInt32Col()
        Value = Float64Col()

    class Test(IsDescription):
        """A description that has several columns."""

        x = Col.from_type("int32", 2, 0, pos=0)
        y = Col.from_kind("float", dflt=1, shape=(2, 3))
        z = UInt8Col(dflt=1)
        color = StringCol(2, dflt=" ")
        # color = UInt32Col(2)
        Info = Info()

        class LInfo(IsDescription):  # noqa: D106
            _v_pos = 1
            name = UInt32Col()
            value = Float64Col(pos=0)
            y2 = Col.from_kind("float", dflt=1, shape=(2, 3), pos=1)
            z2 = UInt8Col(dflt=1)

            class LInfo2(IsDescription):  # noqa: D106
                y3 = Col.from_kind("float", dflt=1, shape=(2, 3))
                z3 = UInt8Col(dflt=1)
                name = UInt32Col()
                value = Float64Col()

                class LInfo3(IsDescription):  # noqa: D106
                    name = UInt32Col()
                    value = Float64Col()
                    y4 = Col.from_kind("float", dflt=1, shape=(2, 3))
                    z4 = UInt8Col(dflt=1)

    #     class Info(IsDescription):
    #         _v_pos = 2
    #         Name = StringCol(itemsize=2)
    #         Value = ComplexCol(itemsize=16)

    #     class Test(IsDescription):
    #         """A description that has several columns"""
    #         x = Col.from_type("int32", 2, 0, pos=0)
    #         y = Col.from_kind('float', dflt=1, shape=(2,3))
    #         z = UInt8Col(dflt=1)
    #         color = StringCol(2, dflt=" ")
    #         Info = Info()
    #         class info(IsDescription):
    #             _v_pos = 1
    #             name = StringCol(itemsize=2)
    #             value = ComplexCol(itemsize=16, pos=0)
    #             y2 = Col.from_kind('float', dflt=1, shape=(2,3), pos=1)
    #             z2 = UInt8Col(dflt=1)
    #             class info2(IsDescription):
    #                 y3 = Col.from_kind('float', dflt=1, shape=(2,3))
    #                 z3 = UInt8Col(dflt=1)
    #                 name = StringCol(itemsize=2)
    #                 value = ComplexCol(itemsize=16)
    #                 class info3(IsDescription):
    #                     name = StringCol(itemsize=2)
    #                     value = ComplexCol(itemsize=16)
    #                     y4 = Col.from_kind('float', dflt=1, shape=(2,3))
    #                     z4 = UInt8Col(dflt=1)

    # example cases of class Test
    klass = Test()
    # klass = Info()
    desc = Description(klass.columns)
    print("Description representation (short) ==>", desc)
    print("Description representation (long) ==>", repr(desc))
    print("Column names ==>", desc._v_names)
    print("Column x ==>", desc.x)
    print("Column Info ==>", desc.Info)
    print("Column Info.value ==>", desc.Info.Value)
    print("Nested column names  ==>", desc._v_nested_names)
    print("Defaults ==>", desc._v_dflts)
    print("Nested Formats ==>", desc._v_nested_formats)
    print("Nested Descriptions ==>", desc._v_nested_descr)
    print("Nested Descriptions (info) ==>", desc.info._v_nested_descr)
    print("Total size ==>", desc._v_dtype.itemsize)

    # check _f_walk
    for obj in desc._f_walk():
        if isinstance(obj, Description):
            print("******begin object*************", end=" ")
            print("name -->", obj._v_name)
            # print("name -->", object._v_dtype.name)
            # print("object childs-->", object._v_names)
            # print("object nested childs-->", object._v_nested_names)
            print("totalsize-->", obj._v_dtype.itemsize)
        else:
            # pass
            print("leaf -->", obj._v_name, obj.dtype)

    class TestDescParent(IsDescription):  # noqa: D101
        c = Int32Col()

    class TestDesc(TestDescParent):  # noqa: D101
        pass

    assert "c" in TestDesc.columns
