import collections
import warnings
from functools import cached_property

from llvmlite import ir

from .abstract import DTypeSpec, IteratorType, MutableSequence, Number, Type
from .common import Buffer, Opaque, SimpleIteratorType
from numba.core.typeconv import Conversion
from numba.core import utils
from .misc import UnicodeType
from .containers import Bytes
import numpy as np

class CharSeq(Type):
    """
    A fixed-length 8-bit character sequence.
    """
    mutable = True

    def __init__(self, count):
        self.count = count
        name = "[char x %d]" % count
        super(CharSeq, self).__init__(name)

    @property
    def key(self):
        return self.count

    def can_convert_from(self, typingctx, other):
        if isinstance(other, Bytes):
            return Conversion.safe


class UnicodeCharSeq(Type):
    """
    A fixed-length unicode character sequence.
    """
    mutable = True

    def __init__(self, count):
        self.count = count
        name = "[unichr x %d]" % count
        super(UnicodeCharSeq, self).__init__(name)

    @property
    def key(self):
        return self.count

    def can_convert_to(self, typingctx, other):
        if isinstance(other, UnicodeCharSeq):
            return Conversion.safe

    def can_convert_from(self, typingctx, other):
        if isinstance(other, UnicodeType):
            # Assuming that unicode_type itemsize is not greater than
            # numpy.dtype('U1').itemsize that UnicodeCharSeq is based
            # on.
            return Conversion.safe

    def __repr__(self):
        return f"UnicodeCharSeq({self.count})"


_RecordField = collections.namedtuple(
    '_RecordField',
    'type,offset,alignment,title',
)


class Record(Type):
    """
    A Record datatype can be mapped to a NumPy structured dtype.
    A record is very flexible since it is laid out as a list of bytes.
    Fields can be mapped to arbitrary points inside it, even if they overlap.

    *fields* is a list of `(name:str, data:dict)`.
        Where `data` is `{ type: Type, offset: int }`
    *size* is an int; the record size
    *aligned* is a boolean; whether the record is ABI aligned.
    """
    mutable = True

    @classmethod
    def make_c_struct(cls, name_types):
        """Construct a Record type from a list of (name:str, type:Types).
        The layout of the structure will follow C.

        Note: only scalar types are supported currently.
        """
        from numba.core.registry import cpu_target

        ctx = cpu_target.target_context
        offset = 0
        fields = []
        lltypes = []
        for k, ty in name_types:
            if not isinstance(ty, (Number, NestedArray)):
                msg = "Only Number and NestedArray types are supported, found: {}. "
                raise TypeError(msg.format(ty))
            if isinstance(ty, NestedArray):
                datatype = ctx.data_model_manager[ty].as_storage_type()
            else:
                datatype = ctx.get_data_type(ty)
            lltypes.append(datatype)
            size = ctx.get_abi_sizeof(datatype)
            align = ctx.get_abi_alignment(datatype)
            # align
            misaligned = offset % align
            if misaligned:
                offset += align - misaligned
            fields.append((k, {
                'type': ty, 'offset': offset, 'alignment': align,
            }))
            offset += size
        # Adjust sizeof structure
        abi_size = ctx.get_abi_sizeof(ir.LiteralStructType(lltypes))
        return Record(fields, size=abi_size, aligned=True)

    def __init__(self, fields, size, aligned):
        fields = self._normalize_fields(fields)
        self.fields = dict(fields)
        self.size = size
        self.aligned = aligned

        # Create description
        descbuf = []
        fmt = "{}[type={};offset={}{}]"
        for k, infos in fields:
            extra = ""
            if infos.alignment is not None:
                extra += ';alignment={}'.format(infos.alignment)
            elif infos.title is not None:
                extra += ';title={}'.format(infos.title)
            descbuf.append(fmt.format(k, infos.type, infos.offset, extra))

        desc = ','.join(descbuf)
        name = 'Record({};{};{})'.format(desc, self.size, self.aligned)
        super(Record, self).__init__(name)

        self.bitwidth = self.dtype.itemsize * 8

    @classmethod
    def _normalize_fields(cls, fields):
        """
        fields:
            [name: str,
             value: {
                 type: Type,
                 offset: int,
                 [ alignment: int ],
                 [ title : str],
             }]
        """
        res = []
        for name, infos in sorted(fields, key=lambda x: (x[1]['offset'], x[0])):
            fd = _RecordField(
                type=infos['type'],
                offset=infos['offset'],
                alignment=infos.get('alignment'),
                title=infos.get('title'),
            )
            res.append((name, fd))
        return res

    @property
    def key(self):
        # Numpy dtype equality doesn't always succeed, use the name instead
        # (https://github.com/numpy/numpy/issues/5715)
        return self.name

    @property
    def mangling_args(self):
        return self.__class__.__name__, (self._code,)

    def __len__(self):
        """Returns the number of fields
        """
        return len(self.fields)

    def offset(self, key):
        """Get the byte offset of a field from the start of the structure.
        """
        return self.fields[key].offset

    def typeof(self, key):
        """Get the type of a field.
        """
        return self.fields[key].type

    def alignof(self, key):
        """Get the specified alignment of the field.

        Since field alignment is optional, this may return None.
        """
        return self.fields[key].alignment

    def has_titles(self):
        """Returns True the record uses titles.
        """
        return any(fd.title is not None for fd in self.fields.values())

    def is_title(self, key):
        """Returns True if the field named *key* is a title.
        """
        return self.fields[key].title == key

    @property
    def members(self):
        """An ordered list of (name, type) for the fields.
        """
        ordered = sorted(self.fields.items(), key=lambda x: x[1].offset)
        return [(k, v.type) for k, v in ordered]

    @property
    def dtype(self):
        from numba.np.numpy_support import as_struct_dtype

        return as_struct_dtype(self)

    def can_convert_to(self, typingctx, other):
        """
        Convert this Record to the *other*.

        This method only implements width subtyping for records.
        """
        from numba.core.errors import NumbaExperimentalFeatureWarning

        if isinstance(other, Record):
            if len(other.fields) > len(self.fields):
                return
            for other_fd, self_fd in zip(other.fields.items(),
                                         self.fields.items()):
                if not other_fd == self_fd:
                    return
            warnings.warn(f"{self} has been considered a subtype of {other} "
                          f" This is an experimental feature.",
                          category=NumbaExperimentalFeatureWarning)
            return Conversion.safe

    def __repr__(self):
        fields = [f"('{f_name}', " +
                  f"{{'type': {repr(f_info.type)}, " +
                  f"'offset': {f_info.offset}, " +
                  f"'alignment': {f_info.alignment}, " +
                  f"'title': {f_info.title}, " +
                  f"}}" +
                  ")"
                  for f_name, f_info in self.fields.items()
                  ]
        fields = "[" + ", ".join(fields) + "]"
        return f"Record({fields}, {self.size}, {self.aligned})"

class DType(DTypeSpec, Opaque):
    """
    Type class associated with the `np.dtype`.

    i.e. :code:`assert type(np.dtype('int32')) == np.dtype`

    np.dtype('int32')
    """

    def __init__(self, dtype):
        assert isinstance(dtype, Type)
        self._dtype = dtype
        name = "dtype(%s)" % (dtype,)
        super(DTypeSpec, self).__init__(name)

    @property
    def key(self):
        return self.dtype

    @property
    def dtype(self):
        return self._dtype

    def __getitem__(self, arg):
        res = super(DType, self).__getitem__(arg)
        return res.copy(dtype=self.dtype)


class NumpyFlatType(SimpleIteratorType, MutableSequence):
    """
    Type class for `ndarray.flat()` objects.
    """

    def __init__(self, arrty):
        self.array_type = arrty
        yield_type = arrty.dtype
        self.dtype = yield_type
        name = "array.flat({arrayty})".format(arrayty=arrty)
        super(NumpyFlatType, self).__init__(name, yield_type)

    @property
    def key(self):
        return self.array_type


class NumpyNdEnumerateType(SimpleIteratorType):
    """
    Type class for `np.ndenumerate()` objects.
    """

    def __init__(self, arrty):
        from . import Tuple, UniTuple, intp
        self.array_type = arrty
        yield_type = Tuple((UniTuple(intp, arrty.ndim), arrty.dtype))
        name = "ndenumerate({arrayty})".format(arrayty=arrty)
        super(NumpyNdEnumerateType, self).__init__(name, yield_type)

    @property
    def key(self):
        return self.array_type


class NumpyNdIterType(IteratorType):
    """
    Type class for `np.nditer()` objects.

    The layout denotes in which order the logical shape is iterated on.
    "C" means logical order (corresponding to in-memory order in C arrays),
    "F" means reverse logical order (corresponding to in-memory order in
    F arrays).
    """

    def __init__(self, arrays):
        # Note inputs arrays can also be scalars, in which case they are
        # broadcast.
        self.arrays = tuple(arrays)
        self.layout = self._compute_layout(self.arrays)
        self.dtypes = tuple(getattr(a, 'dtype', a) for a in self.arrays)
        self.ndim = max(getattr(a, 'ndim', 0) for a in self.arrays)
        name = "nditer(ndim={ndim}, layout={layout}, inputs={arrays})".format(
            ndim=self.ndim, layout=self.layout, arrays=self.arrays)
        super(NumpyNdIterType, self).__init__(name)

    @classmethod
    def _compute_layout(cls, arrays):
        c = collections.Counter()
        for a in arrays:
            if not isinstance(a, Array):
                continue
            if a.layout in 'CF' and a.ndim == 1:
                c['C'] += 1
                c['F'] += 1
            elif a.ndim >= 1:
                c[a.layout] += 1
        return 'F' if c['F'] > c['C'] else 'C'

    @property
    def key(self):
        return self.arrays

    @property
    def views(self):
        """
        The views yielded by the iterator.
        """
        return [Array(dtype, 0, 'C') for dtype in self.dtypes]

    @property
    def yield_type(self):
        from . import BaseTuple
        views = self.views
        if len(views) > 1:
            return BaseTuple.from_types(views)
        else:
            return views[0]

    @cached_property
    def indexers(self):
        """
        A list of (kind, start_dim, end_dim, indices) where:
        - `kind` is either "flat", "indexed", "0d" or "scalar"
        - `start_dim` and `end_dim` are the dimension numbers at which
          this indexing takes place
        - `indices` is the indices of the indexed arrays in self.arrays
        """
        d = collections.OrderedDict()
        layout = self.layout
        ndim = self.ndim
        assert layout in 'CF'
        for i, a in enumerate(self.arrays):
            if not isinstance(a, Array):
                indexer = ('scalar', 0, 0)
            elif a.ndim == 0:
                indexer = ('0d', 0, 0)
            else:
                if a.layout == layout or (a.ndim == 1 and a.layout in 'CF'):
                    kind = 'flat'
                else:
                    kind = 'indexed'
                if layout == 'C':
                    # If iterating in C order, broadcasting is done on the outer indices
                    indexer = (kind, ndim - a.ndim, ndim)
                else:
                    indexer = (kind, 0, a.ndim)
            d.setdefault(indexer, []).append(i)
        return list(k + (v,) for k, v in d.items())

    @cached_property
    def need_shaped_indexing(self):
        """
        Whether iterating on this iterator requires keeping track of
        individual indices inside the shape.  If False, only a single index
        over the equivalent flat shape is required, which can make the
        iterator more efficient.
        """
        for kind, start_dim, end_dim, _ in self.indexers:
            if kind in ('0d', 'scalar'):
                pass
            elif kind == 'flat':
                if (start_dim, end_dim) != (0, self.ndim):
                    # Broadcast flat iteration needs shaped indexing
                    # to know when to restart iteration.
                    return True
            else:
                return True
        return False


class NumpyNdIndexType(SimpleIteratorType):
    """
    Type class for `np.ndindex()` objects.
    """

    def __init__(self, ndim):
        from . import UniTuple, intp
        self.ndim = ndim
        yield_type = UniTuple(intp, self.ndim)
        name = "ndindex(ndim={ndim})".format(ndim=ndim)
        super(NumpyNdIndexType, self).__init__(name, yield_type)

    @property
    def key(self):
        return self.ndim


class Array(Buffer):
    """
    Type class for Numpy arrays.
    """

    def __init__(self, dtype, ndim, layout, readonly=False, name=None,
                 aligned=True):
        if readonly:
            self.mutable = False
        if (not aligned or
            (isinstance(dtype, Record) and not dtype.aligned)):
            self.aligned = False
        if isinstance(dtype, NestedArray):
            ndim += dtype.ndim
            dtype = dtype.dtype
        if name is None:
            type_name = "array"
            if not self.mutable:
                type_name = "readonly " + type_name
            if not self.aligned:
                type_name = "unaligned " + type_name
            name = "%s(%s, %sd, %s)" % (type_name, dtype, ndim, layout)
        super(Array, self).__init__(dtype, ndim, layout, name=name)

    @property
    def mangling_args(self):
        args = [self.dtype, self.ndim, self.layout,
                'mutable' if self.mutable else 'readonly',
                'aligned' if self.aligned else 'unaligned']
        return self.__class__.__name__, args

    def copy(self, dtype=None, ndim=None, layout=None, readonly=None):
        if dtype is None:
            dtype = self.dtype
        if ndim is None:
            ndim = self.ndim
        if layout is None:
            layout = self.layout
        if readonly is None:
            readonly = not self.mutable
        return Array(dtype=dtype, ndim=ndim, layout=layout, readonly=readonly,
                     aligned=self.aligned)

    @property
    def key(self):
        return self.dtype, self.ndim, self.layout, self.mutable, self.aligned

    def unify(self, typingctx, other):
        """
        Unify this with the *other* Array.
        """
        # If other is array and the ndim matches
        if isinstance(other, Array) and other.ndim == self.ndim:
            # If dtype matches or other.dtype is undefined (inferred)
            if other.dtype == self.dtype or not other.dtype.is_precise():
                if self.layout == other.layout:
                    layout = self.layout
                else:
                    layout = 'A'
                readonly = not (self.mutable and other.mutable)
                aligned = self.aligned and other.aligned
                return Array(dtype=self.dtype, ndim=self.ndim, layout=layout,
                             readonly=readonly, aligned=aligned)

    def can_convert_to(self, typingctx, other):
        """
        Convert this Array to the *other*.
        """
        if (isinstance(other, Array) and other.ndim == self.ndim
            and other.dtype == self.dtype):
            if (other.layout in ('A', self.layout)
                and (self.mutable or not other.mutable)
                and (self.aligned or not other.aligned)):
                return Conversion.safe

    def is_precise(self):
        return self.dtype.is_precise()

    @property
    def box_type(self):
        """Returns the Python type to box to.
        """
        return np.ndarray

    def __repr__(self):
        return (
            f"Array({repr(self.dtype)}, {self.ndim}, '{self.layout}', "
            f"{not self.mutable}, aligned={self.aligned})"
                )

class ArrayCTypes(Type):
    """
    This is the type for `np.ndarray.ctypes`.
    """
    def __init__(self, arytype):
        # This depends on the ndim for the shape and strides attributes,
        # even though they are not implemented, yet.
        self.dtype = arytype.dtype
        self.ndim = arytype.ndim
        name = "ArrayCTypes(dtype={0}, ndim={1})".format(self.dtype, self.ndim)
        super(ArrayCTypes, self).__init__(name)

    @property
    def key(self):
        return self.dtype, self.ndim

    def can_convert_to(self, typingctx, other):
        """
        Convert this type to the corresponding pointer type.
        This allows passing a array.ctypes object to a C function taking
        a raw pointer.

        Note that in pure Python, the array.ctypes object can only be
        passed to a ctypes function accepting a c_void_p, not a typed
        pointer.
        """
        from . import CPointer, voidptr
        # XXX what about readonly
        if isinstance(other, CPointer) and other.dtype == self.dtype:
            return Conversion.safe
        elif other == voidptr:
            return Conversion.safe


class ArrayFlags(Type):
    """
    This is the type for `np.ndarray.flags`.
    """
    def __init__(self, arytype):
        self.array_type = arytype
        name = "ArrayFlags({0})".format(self.array_type)
        super(ArrayFlags, self).__init__(name)

    @property
    def key(self):
        return self.array_type


class NestedArray(Array):
    """
    A NestedArray is an array nested within a structured type (which are "void"
    type in NumPy parlance). Unlike an Array, the shape, and not just the number
    of dimensions is part of the type of a NestedArray.
    """

    def __init__(self, dtype, shape):
        if isinstance(dtype, NestedArray):
            tmp = Array(dtype.dtype, dtype.ndim, 'C')
            shape += dtype.shape
            dtype = tmp.dtype
        assert dtype.bitwidth % 8 == 0, \
            "Dtype bitwidth must be a multiple of bytes"
        self._shape = shape
        name = "nestedarray(%s, %s)" % (dtype, shape)
        ndim = len(shape)
        super(NestedArray, self).__init__(dtype, ndim, 'C', name=name)

    @property
    def shape(self):
        return self._shape

    @property
    def nitems(self):
        l = 1
        for s in self.shape:
            l = l * s
        return l

    @property
    def size(self):
        return self.dtype.bitwidth // 8

    @property
    def strides(self):
        stride = self.size
        strides = []
        for i in reversed(self._shape):
             strides.append(stride)
             stride *= i
        return tuple(reversed(strides))

    @property
    def key(self):
        return self.dtype, self.shape

    def __repr__(self):
        return f"NestedArray({repr(self.dtype)}, {self.shape})"


class NumPyRandomBitGeneratorType(Type):
    def __init__(self, *args, **kwargs):
        super(NumPyRandomBitGeneratorType, self).__init__(*args, **kwargs)
        self.name = 'NumPyRandomBitGeneratorType'


class NumPyRandomGeneratorType(Type):
    def __init__(self, *args, **kwargs):
        super(NumPyRandomGeneratorType, self).__init__(*args, **kwargs)
        self.name = 'NumPyRandomGeneratorType'


class PolynomialType(Type):
    def __init__(self, coef, domain=None, window=None, n_args=1):
        super(PolynomialType, self).__init__(name=f'PolynomialType({coef}, {domain}, {domain}, {n_args})')
        self.coef = coef
        self.domain = domain
        self.window = window
        # We use n_args to keep track of the number of arguments in the
        # constructor, since the types of domain and window arguments depend on
        # that and we need that information when boxing
        self.n_args = n_args
