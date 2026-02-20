import sys
from _typeshed import Incomplete
from builtins import memoryview as py_memoryview
from types import EllipsisType
from typing import Any, Final, Literal, Never, Self, SupportsIndex, TypeAlias, TypedDict, type_check_only
from typing_extensions import CapsuleType, ReadOnly

@type_check_only
class CApiDict(TypedDict):
    __pyx_collections_abc_Sequence: ReadOnly[CapsuleType]
    generic: ReadOnly[CapsuleType]
    strided: ReadOnly[CapsuleType]
    indirect: ReadOnly[CapsuleType]
    contiguous: ReadOnly[CapsuleType]
    indirect_contiguous: ReadOnly[CapsuleType]
    __pyx_memoryview_thread_locks_used: ReadOnly[CapsuleType]
    __pyx_memoryview_thread_locks: ReadOnly[CapsuleType]
    _allocate_buffer: ReadOnly[CapsuleType]
    array_cwrapper: ReadOnly[CapsuleType]
    memoryview_cwrapper: ReadOnly[CapsuleType]
    memoryview_check: ReadOnly[CapsuleType]
    _unellipsify: ReadOnly[CapsuleType]
    assert_direct_dimensions: ReadOnly[CapsuleType]
    memview_slice: ReadOnly[CapsuleType]
    slice_memviewslice: ReadOnly[CapsuleType]
    pybuffer_index: ReadOnly[CapsuleType]
    transpose_memslice: ReadOnly[CapsuleType]
    memoryview_fromslice: ReadOnly[CapsuleType]
    get_slice_from_memview: ReadOnly[CapsuleType]
    slice_copy: ReadOnly[CapsuleType]
    memoryview_copy: ReadOnly[CapsuleType]
    memoryview_copy_from_slice: ReadOnly[CapsuleType]
    abs_py_ssize_t: ReadOnly[CapsuleType]
    get_best_order: ReadOnly[CapsuleType]
    _copy_strided_to_strided: ReadOnly[CapsuleType]
    copy_strided_to_strided: ReadOnly[CapsuleType]
    slice_get_size: ReadOnly[CapsuleType]
    fill_contig_strides_array: ReadOnly[CapsuleType]
    copy_data_to_temp: ReadOnly[CapsuleType]
    _err_extents: ReadOnly[CapsuleType]
    _err_dim: ReadOnly[CapsuleType]
    _err: ReadOnly[CapsuleType]
    _err_no_memory: ReadOnly[CapsuleType]
    memoryview_copy_contents: ReadOnly[CapsuleType]
    broadcast_leading: ReadOnly[CapsuleType]
    refcount_copying: ReadOnly[CapsuleType]
    refcount_objects_in_slice_with_gil: ReadOnly[CapsuleType]
    refcount_objects_in_slice: ReadOnly[CapsuleType]
    slice_assign_scalar: ReadOnly[CapsuleType]
    _slice_assign_scalar: ReadOnly[CapsuleType]
    __pyx_unpickle_Enum__set_state: ReadOnly[CapsuleType]
    format_from_typeinfo: ReadOnly[CapsuleType]

_ToIndex: TypeAlias = SupportsIndex | tuple[SupportsIndex, ...] | EllipsisType

###

__pyx_capi__: Final[CApiDict] = ...  # undocumented
__test__: Final[dict[Any, Any]] = ...  # undocumented

class Enum:  # undocumented
    def __init__(self, /, name: str) -> None: ...
    def __setstate__(self, state: tuple[str], /) -> None: ...

class memoryview:  # undocumented
    __pyx_vtable__: Final[CapsuleType]
    base: Final[array]
    itemsize: Final[int]
    nbytes: Final[int]
    ndim: Final[int]
    size: int
    shape: Final[tuple[int, ...]]
    strides: Final[tuple[int, ...]]
    suboffsets: Final[tuple[int, ...]]
    @property
    def T(self, /) -> Self: ...
    def copy(self, /) -> Self: ...
    def copy_fortran(self, /) -> Self: ...
    def is_c_contig(self, /) -> bool: ...
    def is_f_contig(self, /) -> bool: ...

class _memoryviewslize(memoryview):  # undocumented
    def count(self, /, value: Incomplete) -> int: ...
    def index(self, /, value: Incomplete, start: SupportsIndex = 0, stop: SupportsIndex | None = None) -> int: ...

class array:  # undocumented
    memview: Final[memoryview]

    def __new__(
        array,  # pyright: ignore[reportSelfClsParameterName]
        shape: tuple[int, *tuple[int, ...]],
        itemsize: int,
        format: bytes | str,
        mode: Literal["c", "fortran"],
        allocate_buffer: bool = False,
    ) -> Self: ...
    def __len__(self) -> int: ...
    def __getitem__(self, index: _ToIndex, /) -> memoryview | Incomplete: ...
    def __setitem__(self, index: _ToIndex, value: Incomplete, /) -> None: ...
    def __delitem__(self, index: _ToIndex, /) -> None: ...
    def __setstate__(self, state: Never, /) -> None: ...  # will raise `TypeError`
    if sys.version_info >= (3, 12):
        def __buffer__(self, flags: int, /) -> py_memoryview: ...
    def count(self, /, value: Incomplete) -> int: ...
    def index(self, /, value: Incomplete, start: SupportsIndex = 0, stop: SupportsIndex | None = None) -> int: ...
