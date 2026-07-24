from collections.abc import Callable, Iterator, Sequence
from types import GenericAlias
from typing import Any, Generic, Never, Self, SupportsIndex, overload, override
from typing_extensions import TypeVar

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc

from ._rotation import Rotation

__all__ = ["RigidTransform"]

###

_ShapeT = TypeVar("_ShapeT", bound=tuple[int, ...])
_ShapeT_co = TypeVar("_ShapeT_co", bound=tuple[int, ...], default=tuple[Any, ...], covariant=True)
_RigidTransformT = TypeVar("_RigidTransformT", bound=RigidTransform)

type _JustAnyShape = tuple[Never, Never, Never]
type _ToFloatStrictND = onp.ArrayND[npc.floating | npc.integer, _JustAnyShape]

type _Transformable[ShapeT: tuple[int, ...]] = RigidTransform[ShapeT] | Rotation[ShapeT]

###

def normalize_dual_quaternion(dual_quat: onp.ArrayND[np.float64]) -> onp.ArrayND[np.float64]: ...  # undocumented
def _create_transformation_matrix(
    translations: onp.ToFloatND, rotation_matrices: onp.ToFloatND, single: bool
) -> onp.ArrayND[np.float64]: ...  # undocumented

class RigidTransform(Generic[_ShapeT_co]):
    @classmethod
    def __class_getitem__(cls, t: object | type, /) -> GenericAlias: ...

    #
    @property
    def single(self, /) -> bool: ...
    @property
    def shape(self, /) -> _ShapeT_co: ...
    @property
    def rotation(self, /) -> Rotation[_ShapeT_co]: ...
    @property
    def translation(self, /) -> onp.ArrayND[np.float64]: ...

    #
    @overload
    def __init__(self, /, matrix: _ToFloatStrictND, normalize: bool = True, copy: bool = True) -> None: ...
    @overload
    def __init__(
        self: RigidTransform[tuple[()]], /, matrix: onp.ToFloatStrict2D, normalize: bool = True, copy: bool = True
    ) -> None: ...
    @overload
    def __init__(
        self: RigidTransform[tuple[int]], /, matrix: onp.ToFloatStrict3D, normalize: bool = True, copy: bool = True
    ) -> None: ...
    @overload
    def __init__(self, /, matrix: onp.ToFloatND, normalize: bool = True, copy: bool = True) -> None: ...

    #
    def __len__(self: RigidTransform[tuple[int, *tuple[int, ...]]], /) -> int: ...

    #
    @overload
    def __iter__(self: RigidTransform[_JustAnyShape], /) -> Iterator[RigidTransform]: ...
    @overload
    def __iter__(self: RigidTransform[tuple[int]], /) -> Iterator[RigidTransform[tuple[()]]]: ...
    @overload
    def __iter__(self: RigidTransform[tuple[int, int]], /) -> Iterator[RigidTransform[tuple[int]]]: ...
    @overload
    def __iter__(self: RigidTransform[tuple[int, *tuple[int, ...]]], /) -> Iterator[RigidTransform]: ...

    #
    @overload
    def __getitem__(self, indexer: onp.ToInt1D | slice, /) -> RigidTransform[_ShapeT_co]: ...
    @overload
    def __getitem__(self: RigidTransform[_JustAnyShape], indexer: SupportsIndex, /) -> RigidTransform: ...
    @overload
    def __getitem__(self: RigidTransform[tuple[int]], indexer: SupportsIndex, /) -> RigidTransform[tuple[()]]: ...
    @overload
    def __getitem__(self: RigidTransform[tuple[int, int]], indexer: SupportsIndex, /) -> RigidTransform[tuple[int]]: ...
    @overload
    def __getitem__(self: RigidTransform[tuple[int, *tuple[int, ...]]], indexer: SupportsIndex, /) -> RigidTransform: ...

    #
    def __setitem__(
        self: RigidTransform[tuple[int, *tuple[int, ...]]], indexer: SupportsIndex | onp.ToInt1D | slice, value: RigidTransform, /
    ) -> None: ...

    #
    @overload
    def __mul__(self, other: _Transformable[_JustAnyShape], /) -> RigidTransform: ...
    @overload
    def __mul__(self, other: _Transformable[tuple[()]], /) -> RigidTransform[_ShapeT_co]: ...
    @overload
    def __mul__(self: RigidTransform[_JustAnyShape], other: _Transformable[tuple[int, ...]], /) -> RigidTransform: ...
    @overload
    def __mul__(self: RigidTransform[tuple[()]], other: _Transformable[_ShapeT], /) -> RigidTransform[_ShapeT]: ...
    @overload
    def __mul__(self: RigidTransform[tuple[int]], other: _Transformable[tuple[int]], /) -> RigidTransform[tuple[int]]: ...
    @overload
    def __mul__(self, other: _Transformable[tuple[int, ...]], /) -> RigidTransform: ...

    #
    @overload
    def __rmul__(self, other: Rotation[_JustAnyShape], /) -> RigidTransform: ...
    @overload
    def __rmul__(self, other: Rotation[tuple[()]], /) -> RigidTransform[_ShapeT_co]: ...
    @overload
    def __rmul__(self: RigidTransform[_JustAnyShape], other: Rotation[tuple[int, ...]], /) -> RigidTransform: ...
    @overload
    def __rmul__(self: RigidTransform[tuple[()]], other: Rotation[_ShapeT], /) -> RigidTransform[_ShapeT]: ...
    @overload
    def __rmul__(self: RigidTransform[tuple[int]], other: Rotation[tuple[int]], /) -> RigidTransform[tuple[int]]: ...
    @overload
    def __rmul__(self, other: Rotation[tuple[int, ...]], /) -> RigidTransform: ...

    #
    def __pow__(self, n: onp.ToFloat, /) -> Self: ...

    #
    @override
    def __reduce__(self, /) -> tuple[Callable[[Self, onp.ToFloatND], RigidTransform], tuple[onp.ArrayND[np.float64]]]: ...

    #
    @overload
    @staticmethod
    def from_matrix(matrix: _ToFloatStrictND) -> RigidTransform: ...
    @overload
    @staticmethod
    def from_matrix(matrix: onp.ToFloatStrict2D) -> RigidTransform[tuple[()]]: ...
    @overload
    @staticmethod
    def from_matrix(matrix: onp.ToFloatStrict3D) -> RigidTransform[tuple[int]]: ...
    @overload
    @staticmethod
    def from_matrix(matrix: onp.ToFloatND) -> RigidTransform: ...

    #
    @staticmethod
    def from_rotation(rotation: Rotation[_ShapeT]) -> RigidTransform[_ShapeT]: ...

    #
    @overload
    @staticmethod
    def from_translation(translation: _ToFloatStrictND) -> RigidTransform: ...
    @overload
    @staticmethod
    def from_translation(translation: onp.ToFloatStrict1D) -> RigidTransform[tuple[()]]: ...
    @overload
    @staticmethod
    def from_translation(translation: onp.ToFloatStrict2D) -> RigidTransform[tuple[int]]: ...
    @overload
    @staticmethod
    def from_translation(translation: onp.ToFloatND) -> RigidTransform: ...

    #
    @overload
    @staticmethod
    def from_components(translation: _ToFloatStrictND, rotation: Rotation) -> RigidTransform: ...
    @overload
    @staticmethod
    def from_components(translation: onp.ToFloatStrict1D, rotation: Rotation[_ShapeT]) -> RigidTransform[_ShapeT]: ...
    @overload
    @staticmethod
    def from_components(translation: onp.ToFloatStrict2D, rotation: Rotation) -> RigidTransform[tuple[int]]: ...
    @overload
    @staticmethod
    def from_components(translation: onp.ToFloatND, rotation: Rotation) -> RigidTransform: ...

    #
    @overload
    @staticmethod
    def from_exp_coords(exp_coords: _ToFloatStrictND) -> RigidTransform: ...
    @overload
    @staticmethod
    def from_exp_coords(exp_coords: onp.ToFloatStrict1D) -> RigidTransform[tuple[()]]: ...
    @overload
    @staticmethod
    def from_exp_coords(exp_coords: onp.ToFloatStrict2D) -> RigidTransform[tuple[int]]: ...
    @overload
    @staticmethod
    def from_exp_coords(exp_coords: onp.ToFloatND) -> RigidTransform: ...

    #
    @overload
    @staticmethod
    def from_dual_quat(dual_quat: _ToFloatStrictND, *, scalar_first: bool = False) -> RigidTransform: ...
    @overload
    @staticmethod
    def from_dual_quat(dual_quat: onp.ToFloatStrict1D, *, scalar_first: bool = False) -> RigidTransform[tuple[()]]: ...
    @overload
    @staticmethod
    def from_dual_quat(dual_quat: onp.ToFloatStrict2D, *, scalar_first: bool = False) -> RigidTransform[tuple[int]]: ...
    @overload
    @staticmethod
    def from_dual_quat(dual_quat: onp.ToFloatND, *, scalar_first: bool = False) -> RigidTransform: ...

    #
    @overload
    @staticmethod
    def identity(num: None = None, *, shape: None = None) -> RigidTransform[tuple[()]]: ...
    @overload
    @staticmethod
    def identity(num: int, *, shape: None = None) -> RigidTransform[tuple[int]]: ...
    @overload
    @staticmethod
    def identity(num: None = None, *, shape: int) -> RigidTransform[tuple[int]]: ...
    @overload
    @staticmethod
    def identity(num: None = None, *, shape: _ShapeT) -> RigidTransform[_ShapeT]: ...

    #
    @overload
    @staticmethod
    def concatenate(transforms: _RigidTransformT) -> _RigidTransformT: ...
    @overload
    @staticmethod
    def concatenate(transforms: Sequence[RigidTransform[_JustAnyShape]]) -> RigidTransform: ...
    @overload
    @staticmethod
    def concatenate(transforms: Sequence[RigidTransform[tuple[()]]]) -> RigidTransform[tuple[int]]: ...
    @overload
    @staticmethod
    def concatenate(transforms: Sequence[RigidTransform[tuple[int]]]) -> RigidTransform[tuple[int, int]]: ...
    @overload
    @staticmethod
    def concatenate(transforms: Sequence[RigidTransform]) -> RigidTransform: ...  # pyright: ignore[reportOverlappingOverload]

    #
    @overload
    def as_matrix(self: RigidTransform[_JustAnyShape], /) -> onp.ArrayND[np.float64]: ...
    @overload
    def as_matrix(self: RigidTransform[tuple[()]], /) -> onp.Array2D[np.float64]: ...
    @overload
    def as_matrix(self: RigidTransform[tuple[int]], /) -> onp.Array3D[np.float64]: ...
    @overload
    def as_matrix(self, /) -> onp.ArrayND[np.float64]: ...

    #
    @overload
    def as_components(self: RigidTransform[_JustAnyShape], /) -> tuple[onp.ArrayND[np.float64], Rotation]: ...
    @overload
    def as_components(self: RigidTransform[tuple[()]], /) -> tuple[onp.Array1D[np.float64], Rotation[tuple[()]]]: ...
    @overload
    def as_components(self: RigidTransform[tuple[int]], /) -> tuple[onp.Array2D[np.float64], Rotation[tuple[int]]]: ...
    @overload
    def as_components(self, /) -> tuple[onp.ArrayND[np.float64], Rotation[_ShapeT_co]]: ...

    #
    @overload
    def as_exp_coords(self: RigidTransform[_JustAnyShape], /) -> onp.ArrayND[np.float64]: ...
    @overload
    def as_exp_coords(self: RigidTransform[tuple[()]], /) -> onp.Array1D[np.float64]: ...
    @overload
    def as_exp_coords(self: RigidTransform[tuple[int]], /) -> onp.Array2D[np.float64]: ...
    @overload
    def as_exp_coords(self, /) -> onp.ArrayND[np.float64]: ...

    #
    @overload
    def as_dual_quat(self: RigidTransform[_JustAnyShape], /, *, scalar_first: bool = False) -> onp.ArrayND[np.float64]: ...
    @overload
    def as_dual_quat(self: RigidTransform[tuple[()]], /, *, scalar_first: bool = False) -> onp.Array1D[np.float64]: ...
    @overload
    def as_dual_quat(self: RigidTransform[tuple[int]], /, *, scalar_first: bool = False) -> onp.Array2D[np.float64]: ...
    @overload
    def as_dual_quat(self, /, *, scalar_first: bool = False) -> onp.ArrayND[np.float64]: ...

    #
    def inv(self, /) -> Self: ...

    #
    @overload
    def mean(self, /, weights: onp.ToFloatND | None = None, axis: None = None) -> RigidTransform[tuple[()]]: ...
    @overload
    def mean(self, /, weights: onp.ToFloatND | None, axis: int | tuple[int, ...]) -> RigidTransform: ...
    @overload
    def mean(self, /, weights: onp.ToFloatND | None = None, *, axis: int | tuple[int, ...]) -> RigidTransform: ...

    #
    @overload
    def apply(
        self: RigidTransform[_JustAnyShape], /, vector: onp.ToFloatND, inverse: bool = False
    ) -> onp.ArrayND[np.float64]: ...
    @overload
    def apply(
        self: RigidTransform[tuple[()]], /, vector: onp.ToFloatStrict1D, inverse: bool = False
    ) -> onp.Array1D[np.float64]: ...
    @overload
    def apply(
        self: RigidTransform[tuple[int]], /, vector: onp.ToFloatStrict2D, inverse: bool = False
    ) -> onp.Array2D[np.float64]: ...
    @overload
    def apply(self, /, vector: onp.ToFloatND, inverse: bool = False) -> onp.ArrayND[np.float64]: ...
