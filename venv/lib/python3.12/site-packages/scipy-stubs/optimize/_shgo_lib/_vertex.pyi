import abc
from collections import OrderedDict
from collections.abc import Callable, Iterable, Iterator, Sequence
from typing import Concatenate, Final, Generic, TypeAlias, type_check_only
from typing_extensions import TypeVar

import numpy as np
import optype.numpy as onp

_VT_co = TypeVar("_VT_co", bound=VertexBase, default=VertexBase, covariant=True)

_Float1D: TypeAlias = onp.Array1D[np.float64]
_Fun0D: TypeAlias = Callable[Concatenate[_Float1D, ...], onp.ToFloat]
_Fun1D: TypeAlias = Callable[Concatenate[_Float1D, ...], onp.ToFloat1D]

@type_check_only
class _VertexMixin:
    def connect(self, /, v: VertexBase) -> None: ...
    def disconnect(self, /, v: VertexBase) -> None: ...

###

class VertexBase(abc.ABC):
    x: onp.ToFloat1D
    x_a: _Float1D  # lazy
    hash: Final[int]
    index: Final[int | None]
    nn: set[VertexBase]
    st: set[VertexBase]  # might not be set
    feasible: bool  # might not be set

    def __init__(self, /, x: onp.ToFloat1D, nn: Iterable[VertexBase] | None = None, index: int | None = None) -> None: ...
    @abc.abstractmethod
    def connect(self, /, v: VertexBase) -> None: ...
    @abc.abstractmethod
    def disconnect(self, /, v: VertexBase) -> None: ...
    def star(self, /) -> set[VertexBase]: ...

class VertexScalarField(_VertexMixin, VertexBase):
    check_min: bool
    check_max: bool

    def __init__(
        self,
        /,
        x: onp.ToFloat1D,
        field: _Fun0D | None = None,
        nn: Iterable[VertexBase] | None = None,
        index: int | None = None,
        field_args: tuple[object, ...] = (),
        g_cons: _Fun1D | None = None,
        g_cons_args: tuple[object, ...] = (),
    ) -> None: ...
    #
    def minimiser(self, /) -> bool: ...
    def maximiser(self, /) -> bool: ...

class VertexVectorField(_VertexMixin, VertexBase):
    # NOTE: The implementaiton is a WIP
    def __init__(
        self,
        /,
        x: onp.ToFloat1D,
        sfield: _Fun0D | None = None,
        vfield: _Fun1D | None = None,
        field_args: tuple[object, ...] = (),
        vfield_args: tuple[object, ...] = (),
        g_cons: _Fun1D | None = None,
        g_cons_args: tuple[object, ...] = (),
        nn: Iterable[VertexBase] | None = None,
        index: int | None = None,
    ) -> None: ...

class VertexCube(_VertexMixin, VertexBase):
    def __init__(self, /, x: onp.ToFloat1D, nn: Iterable[VertexBase] | None = None, index: int | None = None) -> None: ...

class VertexCacheBase(Generic[_VT_co]):
    cache: OrderedDict[onp.ToFloat1D, _VT_co]
    nfev: int
    index: int

    def __init__(self, /) -> None: ...
    def __iter__(self, /) -> Iterator[_VT_co]: ...
    def size(self, /) -> int: ...
    def print_out(self, /) -> None: ...

class VertexCacheIndex(VertexCacheBase[VertexCube]):
    Vertex: Final[type[VertexCube]]

    def __getitem__(self, x: onp.ToFloat1D, /, nn: None = None) -> VertexCube: ...

class VertexCacheField(VertexCacheBase[VertexScalarField]):
    Vertex: Final[type[VertexScalarField]]

    field: Final[_Fun0D]
    field_args: Final[tuple[object, ...]]
    wfield: Final[FieldWrapper]
    fpool: Final[set[VertexScalarField]]
    process_fpool: Final[Callable[[], None]]

    g_cons: Final[Sequence[_Fun1D]]
    g_cons_args: Final[tuple[object, ...]]
    wgcons: Final[ConstraintWrapper]
    gpool: Final[set[VertexScalarField]]
    process_gpool: Final[Callable[[], None]]

    workers: Final[int]
    index: int
    sfc_lock: bool

    def __init__(
        self,
        /,
        field: _Fun0D | _Fun1D | None = None,
        field_args: tuple[object, ...] = (),
        g_cons: Sequence[_Fun1D] | None = None,
        g_cons_args: tuple[object, ...] = (),
        workers: int = 1,
    ) -> None: ...
    def __getitem__(self, x: onp.ToFloat1D, /, nn: Iterable[VertexBase] | None = None) -> VertexScalarField: ...
    def process_pools(self, /) -> None: ...
    def feasibility_check(self, /, v: VertexBase) -> bool: ...
    def compute_sfield(self, /, v: VertexBase) -> None: ...
    def proc_gpool(self, /) -> None: ...
    def pproc_gpool(self, /) -> None: ...
    def proc_fpool_g(self, /) -> None: ...
    def proc_fpool_nog(self, /) -> None: ...
    def pproc_fpool_g(self, /) -> None: ...
    def pproc_fpool_nog(self, /) -> None: ...
    def proc_minimisers(self, /) -> None: ...

class ConstraintWrapper:
    g_cons: Sequence[_Fun1D]
    g_cons_args: Sequence[tuple[object, ...]]

    def __init__(self, /, g_cons: Sequence[_Fun1D], g_cons_args: Sequence[tuple[object, ...]]) -> None: ...
    def gcons(self, /, v_x_a: _Float1D) -> bool: ...

class FieldWrapper:
    field: _Fun0D | _Fun1D
    field_args: tuple[object, ...]

    def __init__(self, /, field: _Fun0D | _Fun1D, field_args: tuple[object, ...]) -> None: ...
    def func(self, /, v_x_a: _Float1D) -> onp.ToFloat | onp.ToFloat1D: ...
