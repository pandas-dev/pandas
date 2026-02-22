from collections.abc import Iterable
from collections.abc import Sequence as pySequence
from types import MappingProxyType

from .abstract import (
    ConstSized,
    Container,
    Hashable,
    MutableSequence,
    Sequence,
    Type,
    TypeRef,
    Literal,
    InitialValue,
    Poison,
)
from .common import (
    Buffer,
    IterableType,
    SimpleIterableType,
    SimpleIteratorType,
)
from .misc import Undefined, unliteral, Optional, NoneType
from ..typeconv import Conversion
from ..errors import TypingError
from .. import utils


class Pair(Type):
    """
    A heterogeneous pair.
    """

    def __init__(self, first_type, second_type):
        self.first_type = first_type
        self.second_type = second_type
        name = "pair<%s, %s>" % (first_type, second_type)
        super(Pair, self).__init__(name=name)

    @property
    def key(self):
        return self.first_type, self.second_type

    def unify(self, typingctx, other):
        if isinstance(other, Pair):
            first = typingctx.unify_pairs(self.first_type, other.first_type)
            second = typingctx.unify_pairs(self.second_type, other.second_type)
            if first is not None and second is not None:
                return Pair(first, second)


class BaseContainerIterator(SimpleIteratorType):
    """
    Convenience base class for some container iterators.

    Derived classes must implement the *container_class* attribute.
    """

    def __init__(self, container):
        assert isinstance(container, self.container_class), container
        self.container = container
        yield_type = container.dtype
        name = "iter(%s)" % container
        super(BaseContainerIterator, self).__init__(name, yield_type)

    def unify(self, typingctx, other):
        cls = type(self)
        if isinstance(other, cls):
            container = typingctx.unify_pairs(self.container, other.container)
            if container is not None:
                return cls(container)

    @property
    def key(self):
        return self.container


class BaseContainerPayload(Type):
    """
    Convenience base class for some container payloads.

    Derived classes must implement the *container_class* attribute.
    """

    def __init__(self, container):
        assert isinstance(container, self.container_class)
        self.container = container
        name = "payload(%s)" % container
        super(BaseContainerPayload, self).__init__(name)

    @property
    def key(self):
        return self.container


class Bytes(Buffer):
    """
    Type class for Python 3.x bytes objects.
    """

    mutable = False
    # Actually true but doesn't matter since bytes is immutable
    slice_is_copy = False


class ByteArray(Buffer):
    """
    Type class for bytearray objects.
    """

    slice_is_copy = True


class PyArray(Buffer):
    """
    Type class for array.array objects.
    """

    slice_is_copy = True


class MemoryView(Buffer):
    """
    Type class for memoryview objects.
    """


def is_homogeneous(*tys):
    """Are the types homogeneous?
    """
    if tys:
        first, tys = tys[0], tys[1:]
        return not any(t != first for t in tys)
    else:
        # *tys* is empty.
        return False


class BaseTuple(ConstSized, Hashable):
    """
    The base class for all tuple types (with a known size).
    """

    @classmethod
    def from_types(cls, tys, pyclass=None):
        """
        Instantiate the right tuple type for the given element types.
        """
        if pyclass is not None and pyclass is not tuple:
            # A subclass => is it a namedtuple?
            assert issubclass(pyclass, tuple)
            if hasattr(pyclass, "_asdict"):
                tys = tuple(map(unliteral, tys))
                homogeneous = is_homogeneous(*tys)
                if homogeneous:
                    return NamedUniTuple(tys[0], len(tys), pyclass)
                else:
                    return NamedTuple(tys, pyclass)
        else:
            dtype = utils.unified_function_type(tys)
            if dtype is not None:
                return UniTuple(dtype, len(tys))
            # non-named tuple
            homogeneous = is_homogeneous(*tys)
            if homogeneous:
                return cls._make_homogeneous_tuple(tys[0], len(tys))
            else:
                return cls._make_heterogeneous_tuple(tys)

    @classmethod
    def _make_homogeneous_tuple(cls, dtype, count):
        return UniTuple(dtype, count)

    @classmethod
    def _make_heterogeneous_tuple(cls, tys):
        return Tuple(tys)


class BaseAnonymousTuple(BaseTuple):
    """
    Mixin for non-named tuples.
    """

    def can_convert_to(self, typingctx, other):
        """
        Convert this tuple to another one.  Note named tuples are rejected.
        """
        if not isinstance(other, BaseAnonymousTuple):
            return
        if len(self) != len(other):
            return
        if len(self) == 0:
            return Conversion.safe
        if isinstance(other, BaseTuple):
            kinds = [
                typingctx.can_convert(ta, tb) for ta, tb in zip(self, other)
            ]
            if any(kind is None for kind in kinds):
                return
            return max(kinds)

    def __unliteral__(self):
        return type(self).from_types([unliteral(t) for t in self])


class _HomogeneousTuple(Sequence, BaseTuple):
    @property
    def iterator_type(self):
        return UniTupleIter(self)

    def __getitem__(self, i):
        """
        Return element at position i
        """
        return self.dtype

    def __iter__(self):
        return iter([self.dtype] * self.count)

    def __len__(self):
        return self.count

    @property
    def types(self):
        return (self.dtype,) * self.count


class UniTuple(BaseAnonymousTuple, _HomogeneousTuple, Sequence):
    """
    Type class for homogeneous tuples.
    """

    def __init__(self, dtype, count):
        self.dtype = dtype
        self.count = count
        name = "%s(%s x %d)" % (self.__class__.__name__, dtype, count,)
        super(UniTuple, self).__init__(name)

    @property
    def mangling_args(self):
        return self.__class__.__name__, (self.dtype, self.count)

    @property
    def key(self):
        return self.dtype, self.count

    def unify(self, typingctx, other):
        """
        Unify UniTuples with their dtype
        """
        if isinstance(other, UniTuple) and len(self) == len(other):
            dtype = typingctx.unify_pairs(self.dtype, other.dtype)
            if dtype is not None:
                return UniTuple(dtype=dtype, count=self.count)

    def __unliteral__(self):
        return type(self)(dtype=unliteral(self.dtype), count=self.count)

    def __repr__(self):
        return f"UniTuple({repr(self.dtype)}, {self.count})"


class UniTupleIter(BaseContainerIterator):
    """
    Type class for homogeneous tuple iterators.
    """

    container_class = _HomogeneousTuple


class _HeterogeneousTuple(BaseTuple):
    def __getitem__(self, i):
        """
        Return element at position i
        """
        return self.types[i]

    def __len__(self):
        # Beware: this makes Tuple(()) false-ish
        return len(self.types)

    def __iter__(self):
        return iter(self.types)

    @staticmethod
    def is_types_iterable(types):
        # issue 4463 - check if argument 'types' is iterable
        if not isinstance(types, Iterable):
            raise TypingError("Argument 'types' is not iterable")


class UnionType(Type):
    def __init__(self, types):
        self.types = tuple(sorted(set(types), key=lambda x: x.name))
        name = "Union[{}]".format(",".join(map(str, self.types)))
        super(UnionType, self).__init__(name=name)

    def get_type_tag(self, typ):
        return self.types.index(typ)


class Tuple(BaseAnonymousTuple, _HeterogeneousTuple):
    def __new__(cls, types):

        t = utils.unified_function_type(types, require_precise=True)
        if t is not None:
            return UniTuple(dtype=t, count=len(types))

        _HeterogeneousTuple.is_types_iterable(types)

        if types and all(t == types[0] for t in types[1:]):
            return UniTuple(dtype=types[0], count=len(types))
        else:
            return object.__new__(Tuple)

    def __init__(self, types):
        self.types = tuple(types)
        self.count = len(self.types)
        self.dtype = UnionType(types)
        name = "%s(%s)" % (
            self.__class__.__name__,
            ", ".join(str(i) for i in self.types),
        )
        super(Tuple, self).__init__(name)

    @property
    def mangling_args(self):
        return self.__class__.__name__, tuple(t for t in self.types)

    @property
    def key(self):
        return self.types

    def unify(self, typingctx, other):
        """
        Unify elements of Tuples/UniTuples
        """
        # Other is UniTuple or Tuple
        if isinstance(other, BaseTuple) and len(self) == len(other):
            unified = [
                typingctx.unify_pairs(ta, tb) for ta, tb in zip(self, other)
            ]

            if all(t is not None for t in unified):
                return Tuple(unified)

    def __repr__(self):
        return f"Tuple({tuple(ty for ty in self.types)})"


class _StarArgTupleMixin:
    @classmethod
    def _make_homogeneous_tuple(cls, dtype, count):
        return StarArgUniTuple(dtype, count)

    @classmethod
    def _make_heterogeneous_tuple(cls, tys):
        return StarArgTuple(tys)


class StarArgTuple(_StarArgTupleMixin, Tuple):
    """To distinguish from Tuple() used as argument to a `*args`.
    """

    def __new__(cls, types):
        _HeterogeneousTuple.is_types_iterable(types)

        if types and all(t == types[0] for t in types[1:]):
            return StarArgUniTuple(dtype=types[0], count=len(types))
        else:
            return object.__new__(StarArgTuple)


class StarArgUniTuple(_StarArgTupleMixin, UniTuple):
    """To distinguish from UniTuple() used as argument to a `*args`.
    """


class BaseNamedTuple(BaseTuple):
    pass


class NamedUniTuple(_HomogeneousTuple, BaseNamedTuple):
    def __init__(self, dtype, count, cls):
        self.dtype = dtype
        self.count = count
        self.fields = tuple(cls._fields)
        self.instance_class = cls
        name = "%s(%s x %d)" % (cls.__name__, dtype, count)
        super(NamedUniTuple, self).__init__(name)

    @property
    def iterator_type(self):
        return UniTupleIter(self)

    @property
    def key(self):
        return self.instance_class, self.dtype, self.count


class NamedTuple(_HeterogeneousTuple, BaseNamedTuple):
    def __init__(self, types, cls):
        _HeterogeneousTuple.is_types_iterable(types)

        self.types = tuple(types)
        self.count = len(self.types)
        self.fields = tuple(cls._fields)
        self.instance_class = cls
        name = "%s(%s)" % (cls.__name__, ", ".join(str(i) for i in self.types))
        super(NamedTuple, self).__init__(name)

    @property
    def key(self):
        return self.instance_class, self.types


class List(MutableSequence, InitialValue):
    """
    Type class for (arbitrary-sized) homogeneous lists.
    """

    def __init__(self, dtype, reflected=False, initial_value=None):
        dtype = unliteral(dtype)
        self.dtype = dtype
        self.reflected = reflected
        cls_name = "reflected list" if reflected else "list"
        name = "%s(%s)<iv=%s>" % (cls_name, self.dtype, initial_value)
        super(List, self).__init__(name=name)
        InitialValue.__init__(self, initial_value)

    def copy(self, dtype=None, reflected=None):
        if dtype is None:
            dtype = self.dtype
        if reflected is None:
            reflected = self.reflected
        return List(dtype, reflected, self.initial_value)

    def unify(self, typingctx, other):
        if isinstance(other, List):
            dtype = typingctx.unify_pairs(self.dtype, other.dtype)
            reflected = self.reflected or other.reflected
            if dtype is not None:
                siv = self.initial_value
                oiv = other.initial_value
                if siv is not None and oiv is not None:
                    use = siv
                    if siv is None:
                        use = oiv
                    return List(dtype, reflected, use)
                else:
                    return List(dtype, reflected)

    @property
    def key(self):
        return self.dtype, self.reflected, str(self.initial_value)

    @property
    def iterator_type(self):
        return ListIter(self)

    def is_precise(self):
        return self.dtype.is_precise()

    def __getitem__(self, args):
        """
        Overrides the default __getitem__ from Type.
        """
        return self.dtype

    def __unliteral__(self):
        return List(self.dtype, reflected=self.reflected,
                    initial_value=None)

    def __repr__(self):
        return f"List({self.dtype}, {self.reflected})"


class LiteralList(Literal, ConstSized, Hashable):
    """A heterogeneous immutable list (basically a tuple with list semantics).
    """

    mutable = False

    def __init__(self, literal_value):
        self.is_types_iterable(literal_value)
        self._literal_init(list(literal_value))
        self.types = tuple(literal_value)
        self.count = len(self.types)
        self.name = "LiteralList({})".format(literal_value)

    def __getitem__(self, i):
        """
        Return element at position i
        """
        return self.types[i]

    def __len__(self):
        return len(self.types)

    def __iter__(self):
        return iter(self.types)

    @classmethod
    def from_types(cls, tys):
        return LiteralList(tys)

    @staticmethod
    def is_types_iterable(types):
        if not isinstance(types, Iterable):
            raise TypingError("Argument 'types' is not iterable")

    @property
    def iterator_type(self):
        return ListIter(self)

    def __unliteral__(self):
        return Poison(self)

    def unify(self, typingctx, other):
        """
        Unify this with the *other* one.
        """
        if isinstance(other, LiteralList) and self.count == other.count:
            tys = []
            for i1, i2 in zip(self.types, other.types):
                tys.append(typingctx.unify_pairs(i1, i2))
            if all(tys):
                return LiteralList(tys)


class ListIter(BaseContainerIterator):
    """
    Type class for list iterators.
    """

    container_class = List


class ListPayload(BaseContainerPayload):
    """
    Internal type class for the dynamically-allocated payload of a list.
    """

    container_class = List


class Set(Container):
    """
    Type class for homogeneous sets.
    """

    mutable = True

    def __init__(self, dtype, reflected=False):
        assert isinstance(dtype, (Hashable, Undefined))
        self.dtype = dtype
        self.reflected = reflected
        cls_name = "reflected set" if reflected else "set"
        name = "%s(%s)" % (cls_name, self.dtype)
        super(Set, self).__init__(name=name)

    @property
    def key(self):
        return self.dtype, self.reflected

    @property
    def iterator_type(self):
        return SetIter(self)

    def is_precise(self):
        return self.dtype.is_precise()

    def copy(self, dtype=None, reflected=None):
        if dtype is None:
            dtype = self.dtype
        if reflected is None:
            reflected = self.reflected
        return Set(dtype, reflected)

    def unify(self, typingctx, other):
        if isinstance(other, Set):
            dtype = typingctx.unify_pairs(self.dtype, other.dtype)
            reflected = self.reflected or other.reflected
            if dtype is not None:
                return Set(dtype, reflected)

    def __repr__(self):
        return f"Set({self.dtype}, {self.reflected})"


class SetIter(BaseContainerIterator):
    """
    Type class for set iterators.
    """

    container_class = Set


class SetPayload(BaseContainerPayload):
    """
    Internal type class for the dynamically-allocated payload of a set.
    """

    container_class = Set


class SetEntry(Type):
    """
    Internal type class for the entries of a Set's hash table.
    """

    def __init__(self, set_type):
        self.set_type = set_type
        name = "entry(%s)" % set_type
        super(SetEntry, self).__init__(name)

    @property
    def key(self):
        return self.set_type


class ListType(IterableType):
    """
    List type
    """

    mutable = True

    def __init__(self, itemty):
        assert not isinstance(itemty, TypeRef)
        itemty = unliteral(itemty)
        if isinstance(itemty, Optional):
            fmt = "List.item_type cannot be of type {}"
            raise TypingError(fmt.format(itemty))
        # FIXME: _sentry_forbidden_types(itemty)
        self.item_type = itemty
        self.dtype = itemty
        name = "{}[{}]".format(self.__class__.__name__, itemty,)
        super(ListType, self).__init__(name)

    @property
    def key(self):
        return self.item_type

    def is_precise(self):
        return not isinstance(self.item_type, Undefined)

    @property
    def iterator_type(self):
        return ListTypeIterableType(self).iterator_type

    @classmethod
    def refine(cls, itemty):
        """Refine to a precise list type
        """
        res = cls(itemty)
        assert res.is_precise()
        return res

    def unify(self, typingctx, other):
        """
        Unify this with the *other* list.
        """
        # If other is list
        if isinstance(other, ListType):
            if not other.is_precise():
                return self

    def __repr__(self):
        return f"ListType({repr(self.item_type)})"


class ListTypeIterableType(SimpleIterableType):
    """
    List iterable type
    """

    def __init__(self, parent):
        assert isinstance(parent, ListType)
        self.parent = parent
        self.yield_type = self.parent.item_type
        name = "list[{}]".format(self.parent.name)
        iterator_type = ListTypeIteratorType(self)
        super(ListTypeIterableType, self).__init__(name, iterator_type)


class ListTypeIteratorType(SimpleIteratorType):
    def __init__(self, iterable):
        self.parent = iterable.parent
        self.iterable = iterable
        yield_type = iterable.yield_type
        name = "iter[{}->{}]".format(iterable.parent, yield_type)
        super(ListTypeIteratorType, self).__init__(name, yield_type)


def _sentry_forbidden_types(key, value):
    # Forbids List and Set for now
    if isinstance(key, (Set, List)):
        raise TypingError("{} as key is forbidden".format(key))
    if isinstance(value, (Set, List)):
        raise TypingError("{} as value is forbidden".format(value))


class DictType(IterableType, InitialValue):
    """Dictionary type
    """

    def __init__(self, keyty, valty, initial_value=None):
        assert not isinstance(keyty, TypeRef)
        assert not isinstance(valty, TypeRef)
        keyty = unliteral(keyty)
        valty = unliteral(valty)
        if isinstance(keyty, (Optional, NoneType)):
            fmt = "Dict.key_type cannot be of type {}"
            raise TypingError(fmt.format(keyty))
        if isinstance(valty, (Optional, NoneType)):
            fmt = "Dict.value_type cannot be of type {}"
            raise TypingError(fmt.format(valty))
        _sentry_forbidden_types(keyty, valty)
        self.key_type = keyty
        self.value_type = valty
        self.keyvalue_type = Tuple([keyty, valty])
        name = "{}[{},{}]<iv={}>".format(
            self.__class__.__name__, keyty, valty, initial_value
        )
        super(DictType, self).__init__(name)
        InitialValue.__init__(self, initial_value)

    def is_precise(self):
        return not any(
            (
                isinstance(self.key_type, Undefined),
                isinstance(self.value_type, Undefined),
            )
        )

    @property
    def iterator_type(self):
        return DictKeysIterableType(self).iterator_type

    @classmethod
    def refine(cls, keyty, valty):
        """
        Refine to a precise dictionary type
        """
        res = cls(keyty, valty)
        assert res.is_precise()
        return res

    def unify(self, typingctx, other):
        """
        Unify this with the *other* dictionary.
        """
        # If other is dict
        if isinstance(other, DictType):
            if not other.is_precise():
                return self
            else:
                ukey_type = self.key_type == other.key_type
                uvalue_type = self.value_type == other.value_type
                if ukey_type and uvalue_type:
                    siv = self.initial_value
                    oiv = other.initial_value
                    siv_none = siv is None
                    oiv_none = oiv is None
                    if not siv_none and not oiv_none:
                        if siv == oiv:
                            return DictType(self.key_type, other.value_type,
                                            siv)
                    return DictType(self.key_type, other.value_type)

    @property
    def key(self):
        return self.key_type, self.value_type, str(self.initial_value)

    def __unliteral__(self):
        return DictType(self.key_type, self.value_type)

    def __repr__(self):
        return f"DictType({self.key_type}, {self.value_type})"


class LiteralStrKeyDict(Literal, ConstSized, Hashable):
    """A Dictionary of string keys to heterogeneous values (basically a
    namedtuple with dict semantics).
    """

    class FakeNamedTuple(pySequence):
        # This is namedtuple-like and is a workaround for #6518 and #7416.
        # This has the couple of namedtuple properties that are used by Numba's
        # internals but avoids use of an actual namedtuple as it cannot have
        # numeric field names, i.e. `namedtuple('foo', '0 1')` is invalid.
        def __init__(self, name, keys):
            self.__name__ = name
            self._fields = tuple(keys)
            super(LiteralStrKeyDict.FakeNamedTuple, self).__init__()

        def __len__(self):
            return len(self._fields)

        def __getitem__(self, key):
            return self._fields[key]

    mutable = False

    def __init__(self, literal_value, value_index=None):
        self._literal_init(literal_value)
        self.value_index = value_index
        strkeys = [x.literal_value for x in literal_value.keys()]
        self.tuple_ty = self.FakeNamedTuple("_ntclazz", strkeys)
        tys = [x for x in literal_value.values()]
        self.types = tuple(tys)
        self.count = len(self.types)
        self.fields = tuple(self.tuple_ty._fields)
        self.instance_class = self.tuple_ty
        self.name = "LiteralStrKey[Dict]({})".format(literal_value)

    def __unliteral__(self):
        return Poison(self)

    def unify(self, typingctx, other):
        """
        Unify this with the *other* one.
        """
        if isinstance(other, LiteralStrKeyDict):
            tys = []
            for (k1, v1), (k2, v2) in zip(
                self.literal_value.items(), other.literal_value.items()
            ):
                if k1 != k2:  # keys must be same
                    break
                tys.append(typingctx.unify_pairs(v1, v2))
            else:
                if all(tys):
                    d = {k: v for k, v in zip(self.literal_value.keys(), tys)}
                    return LiteralStrKeyDict(d)

    def __len__(self):
        return len(self.types)

    def __iter__(self):
        return iter(self.types)

    @property
    def key(self):
        # use the namedtuple fields not the namedtuple itself as it's created
        # locally in the ctor and comparison would always be False.
        return self.tuple_ty._fields, self.types, str(self.literal_value)


class DictItemsIterableType(SimpleIterableType):
    """Dictionary iterable type for .items()
    """

    def __init__(self, parent):
        assert isinstance(parent, DictType)
        self.parent = parent
        self.yield_type = self.parent.keyvalue_type
        name = "items[{}]".format(self.parent.name)
        self.name = name
        iterator_type = DictIteratorType(self)
        super(DictItemsIterableType, self).__init__(name, iterator_type)


class DictKeysIterableType(SimpleIterableType):
    """Dictionary iterable type for .keys()
    """

    def __init__(self, parent):
        assert isinstance(parent, DictType)
        self.parent = parent
        self.yield_type = self.parent.key_type
        name = "keys[{}]".format(self.parent.name)
        self.name = name
        iterator_type = DictIteratorType(self)
        super(DictKeysIterableType, self).__init__(name, iterator_type)


class DictValuesIterableType(SimpleIterableType):
    """Dictionary iterable type for .values()
    """

    def __init__(self, parent):
        assert isinstance(parent, DictType)
        self.parent = parent
        self.yield_type = self.parent.value_type
        name = "values[{}]".format(self.parent.name)
        self.name = name
        iterator_type = DictIteratorType(self)
        super(DictValuesIterableType, self).__init__(name, iterator_type)


class DictIteratorType(SimpleIteratorType):
    def __init__(self, iterable):
        self.parent = iterable.parent
        self.iterable = iterable
        yield_type = iterable.yield_type
        name = "iter[{}->{}],{}".format(
            iterable.parent, yield_type, iterable.name
        )
        super(DictIteratorType, self).__init__(name, yield_type)


class StructRef(Type):
    """A mutable struct.
    """

    def __init__(self, fields):
        """
        Parameters
        ----------
        fields : Sequence
            A sequence of field descriptions, which is a 2-tuple-like object
            containing `(name, type)`, where `name` is a `str` for the field
            name, and `type` is a numba type for the field type.
        """

        def check_field_pair(fieldpair):
            name, typ = fieldpair
            if not isinstance(name, str):
                msg = "expecting a str for field name"
                raise ValueError(msg)
            if not isinstance(typ, Type):
                msg = "expecting a Numba Type for field type"
                raise ValueError(msg)
            return name, typ

        fields = tuple(map(check_field_pair, fields))
        self._fields = tuple(map(check_field_pair,
                                 self.preprocess_fields(fields)))
        self._typename = self.__class__.__qualname__
        name = f"numba.{self._typename}{self._fields}"
        super().__init__(name=name)

    def preprocess_fields(self, fields):
        """Subclasses can override this to do additional clean up on fields.

        The default is an identity function.

        Parameters:
        -----------
        fields : Sequence[Tuple[str, Type]]
        """
        return fields

    @property
    def field_dict(self):
        """Return an immutable mapping for the field names and their
        corresponding types.
        """
        return MappingProxyType(dict(self._fields))

    def get_data_type(self):
        """Get the payload type for the actual underlying structure referred
        to by this struct reference.

        See also: `ClassInstanceType.get_data_type`
        """
        return StructRefPayload(
            typename=self.__class__.__name__, fields=self._fields,
        )


class StructRefPayload(Type):
    """The type of the payload of a mutable struct.
    """

    mutable = True

    def __init__(self, typename, fields):
        self._typename = typename
        self._fields = tuple(fields)
        super().__init__(name=f"numba.{typename}{self._fields}.payload")

    @property
    def field_dict(self):
        return MappingProxyType(dict(self._fields))
