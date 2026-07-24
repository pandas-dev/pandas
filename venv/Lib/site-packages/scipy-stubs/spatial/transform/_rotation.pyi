from collections.abc import Iterator, Sequence
from types import EllipsisType, GenericAlias, ModuleType
from typing import Any, Final, Generic, Literal as L, Never, Self, SupportsIndex, overload, override
from typing_extensions import TypeVar

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc

###

type _RotOrder = L["e", "extrinsic", "i", "intrinsic"]
type _RotGroup = L["I", "O", "T", "D", "Dn", "C", "Cn"]
type _RotAxisSeq = L[
    "xyz", "xzy", "yxz", "yzx", "zxy", "zyx",
    "xyx", "xzx", "yxy", "yzy", "zxz", "zyz",
    "XYZ", "XZY", "YXZ", "YZX", "ZXY", "ZYX",
    "XYX", "XZX", "YXY", "YZY", "ZXZ", "ZYZ",
]  # fmt: skip
type _RotAxis = L["X", "Y", "Z"]

type _JustAnyShape = tuple[Never, Never, Never, Never]
type _ToFloatStrictND = onp.ArrayND[npc.floating | npc.integer, _JustAnyShape]

_ShapeT_co = TypeVar("_ShapeT_co", bound=tuple[int, ...], default=tuple[Any, ...], covariant=True)

###

backend_registry: Final[dict[ModuleType, ModuleType]] = ...  # undocumented
rotation_extra_note: Final[str] = ...  # undocumented

# either `._rotation_cy` or `._rotation_xp`
def select_backend(xp: ModuleType, cython_compatible: bool) -> ModuleType: ...  # undocumented

class Rotation(Generic[_ShapeT_co]):
    @classmethod
    def __class_getitem__(cls, t: object | type, /) -> GenericAlias: ...

    #
    @property
    def single(self, /) -> bool: ...
    @property
    def shape(self, /) -> _ShapeT_co: ...

    #
    @overload
    def __init__(
        self, /, quat: _ToFloatStrictND, normalize: bool = True, copy: bool = True, scalar_first: bool = False
    ) -> None: ...
    @overload
    def __init__(
        self: Rotation[tuple[()]],
        /,
        quat: onp.ToFloatStrict1D,
        normalize: bool = True,
        copy: bool = True,
        scalar_first: bool = False,
    ) -> None: ...
    @overload
    def __init__(
        self: Rotation[tuple[int]],
        /,
        quat: onp.ToFloatStrict2D,
        normalize: bool = True,
        copy: bool = True,
        scalar_first: bool = False,
    ) -> None: ...
    @overload
    def __init__(self, /, quat: onp.ToFloatND, normalize: bool = True, copy: bool = True, scalar_first: bool = False) -> None: ...

    #
    @override
    def __getstate__(self) -> tuple[onp.ArrayND[np.float64], bool]: ...
    def __setstate__(self, state: tuple[onp.ArrayND[np.float64], bool]) -> None: ...

    #
    def __bool__(self, /) -> L[True]: ...
    def __len__(self, /) -> int: ...
    def __getitem__(self, key: SupportsIndex | onp.SequenceND[SupportsIndex] | onp.ToIntND | slice, /) -> Rotation: ...
    def __setitem__(self, indexer: int | slice | EllipsisType | None, value: Rotation, /) -> None: ...
    def __iter__(self) -> Iterator[Rotation[tuple[()]]]: ...
    def __mul__(self, other: Rotation, /) -> Rotation: ...
    def __pow__(self, n: onp.ToFloat, /, modulus: None = None) -> Rotation[_ShapeT_co]: ...

    #
    def inv(self, /) -> Self: ...

    #
    @overload
    def mean(self, /, weights: onp.ToFloatND | None = None, axis: None = None) -> Rotation[tuple[()]]: ...
    @overload
    def mean(self, /, weights: onp.ToFloatND | None, axis: int | tuple[int, ...]) -> Rotation: ...
    @overload
    def mean(self, /, weights: onp.ToFloatND | None = None, *, axis: int | tuple[int, ...]) -> Rotation: ...

    #
    @overload
    def as_matrix(self: Rotation[_JustAnyShape], /) -> onp.ArrayND[np.float64]: ...
    @overload
    def as_matrix(self: Rotation[tuple[()]], /) -> onp.Array2D[np.float64]: ...
    @overload
    def as_matrix(self: Rotation[tuple[int]], /) -> onp.Array3D[np.float64]: ...
    @overload
    def as_matrix(self, /) -> onp.ArrayND[np.float64]: ...

    #
    @overload
    def as_mrp(self: Rotation[_JustAnyShape], /) -> onp.ArrayND[np.float64]: ...
    @overload
    def as_mrp(self: Rotation[tuple[()]], /) -> onp.Array1D[np.float64]: ...
    @overload
    def as_mrp(self: Rotation[tuple[int]], /) -> onp.Array2D[np.float64]: ...
    @overload
    def as_mrp(self, /) -> onp.ArrayND[np.float64]: ...

    #
    @overload
    def as_rotvec(self: Rotation[_JustAnyShape], /, degrees: bool = False) -> onp.ArrayND[np.float64]: ...
    @overload
    def as_rotvec(self: Rotation[tuple[()]], /, degrees: bool = False) -> onp.Array1D[np.float64]: ...
    @overload
    def as_rotvec(self: Rotation[tuple[int]], /, degrees: bool = False) -> onp.Array2D[np.float64]: ...
    @overload
    def as_rotvec(self, /, degrees: bool = False) -> onp.ArrayND[np.float64]: ...

    #
    @overload
    def as_euler(
        self: Rotation[_JustAnyShape], /, seq: _RotAxisSeq, degrees: bool = False, *, suppress_warnings: bool = False
    ) -> onp.ArrayND[np.float64]: ...
    @overload
    def as_euler(
        self: Rotation[tuple[()]], /, seq: _RotAxisSeq, degrees: bool = False, *, suppress_warnings: bool = False
    ) -> onp.Array1D[np.float64]: ...
    @overload
    def as_euler(
        self: Rotation[tuple[int]], /, seq: _RotAxisSeq, degrees: bool = False, *, suppress_warnings: bool = False
    ) -> onp.Array2D[np.float64]: ...
    @overload
    def as_euler(
        self, /, seq: _RotAxisSeq, degrees: bool = False, *, suppress_warnings: bool = False
    ) -> onp.ArrayND[np.float64]: ...

    #
    @overload
    def as_quat(
        self: Rotation[_JustAnyShape], /, canonical: bool = False, *, scalar_first: bool = False
    ) -> onp.ArrayND[np.float64]: ...
    @overload
    def as_quat(
        self: Rotation[tuple[()]], /, canonical: bool = False, *, scalar_first: bool = False
    ) -> onp.Array1D[np.float64]: ...
    @overload
    def as_quat(
        self: Rotation[tuple[int]], /, canonical: bool = False, *, scalar_first: bool = False
    ) -> onp.Array2D[np.float64]: ...
    @overload
    def as_quat(self, /, canonical: bool = False, *, scalar_first: bool = False) -> onp.ArrayND[np.float64]: ...

    #
    @overload
    def as_davenport(
        self: Rotation[_JustAnyShape],
        /,
        axes: onp.ToFloatND,
        order: _RotOrder,
        degrees: bool = False,
        *,
        suppress_warnings: bool = False,
    ) -> onp.ArrayND[np.float64]: ...
    @overload
    def as_davenport(
        self: Rotation[tuple[()]],
        /,
        axes: onp.ToFloatND,
        order: _RotOrder,
        degrees: bool = False,
        *,
        suppress_warnings: bool = False,
    ) -> onp.Array1D[np.float64]: ...
    @overload
    def as_davenport(
        self: Rotation[tuple[int]],
        /,
        axes: onp.ToFloatND,
        order: _RotOrder,
        degrees: bool = False,
        *,
        suppress_warnings: bool = False,
    ) -> onp.Array2D[np.float64]: ...
    @overload
    def as_davenport(
        self, /, axes: onp.ToFloatND, order: _RotOrder, degrees: bool = False, *, suppress_warnings: bool = False
    ) -> onp.ArrayND[np.float64]: ...

    #
    @overload
    def magnitude(self: Rotation[_JustAnyShape], /) -> np.float64 | onp.ArrayND[np.float64]: ...
    @overload
    def magnitude(self: Rotation[tuple[()]], /) -> np.float64: ...
    @overload
    def magnitude(self: Rotation[tuple[int]], /) -> onp.Array1D[np.float64]: ...
    @overload
    def magnitude(self: Rotation[tuple[int, int]], /) -> onp.Array2D[np.float64]: ...
    @overload
    def magnitude(self, /) -> np.float64 | onp.ArrayND[np.float64]: ...

    #
    @overload
    def apply(self: Rotation[_JustAnyShape], /, vectors: onp.ToFloatND, inverse: bool = False) -> onp.ArrayND[np.float64]: ...
    @overload
    def apply(self: Rotation[tuple[()]], /, vectors: onp.ToFloatStrict1D, inverse: bool = False) -> onp.Array1D[np.float64]: ...
    @overload
    def apply(self: Rotation[tuple[int]], /, vectors: onp.ToFloatStrict2D, inverse: bool = False) -> onp.Array2D[np.float64]: ...
    @overload
    def apply(self, /, vectors: onp.ToFloatND, inverse: bool = False) -> onp.ArrayND[np.float64]: ...

    # returns `np.bool` if both are single since 1.18.0
    @overload
    def approx_equal(
        self: Rotation[_JustAnyShape], /, other: Rotation[_JustAnyShape], atol: float | None = None, degrees: bool = False
    ) -> onp.ArrayND[np.bool] | np.bool: ...
    @overload
    def approx_equal(
        self: Rotation[tuple[()]], /, other: Rotation[tuple[()]], atol: float | None = None, degrees: bool = False
    ) -> np.bool: ...
    @overload
    def approx_equal(
        self: Rotation[tuple[int] | tuple[int, int] | tuple[int, int, int]],
        /,
        other: Rotation,
        atol: float | None = None,
        degrees: bool = False,
    ) -> onp.ArrayND[np.bool]: ...
    @overload
    def approx_equal(
        self: Rotation,
        /,
        other: Rotation[tuple[int] | tuple[int, int] | tuple[int, int, int]],
        atol: float | None = None,
        degrees: bool = False,
    ) -> onp.ArrayND[np.bool]: ...

    #
    @overload
    def reduce(
        self, /, left: Rotation | None = None, right: Rotation | None = None, return_indices: L[False] = False
    ) -> Self: ...
    @overload
    def reduce(
        self, /, left: Rotation | None, right: Rotation | None, return_indices: L[True]
    ) -> tuple[Self, onp.ArrayND[np.int32 | np.int64], onp.ArrayND[np.int32 | np.int64]]: ...
    @overload
    def reduce(
        self, /, left: Rotation | None = None, right: Rotation | None = None, *, return_indices: L[True]
    ) -> tuple[Self, onp.ArrayND[np.int32 | np.int64], onp.ArrayND[np.int32 | np.int64]]: ...

    #
    @overload
    @staticmethod
    def from_matrix(matrix: _ToFloatStrictND, *, assume_valid: bool = False) -> Rotation: ...
    @overload
    @staticmethod
    def from_matrix(matrix: onp.ToFloatStrict2D, *, assume_valid: bool = False) -> Rotation[tuple[()]]: ...
    @overload
    @staticmethod
    def from_matrix(matrix: onp.ToFloatStrict3D, *, assume_valid: bool = False) -> Rotation[tuple[int]]: ...
    @overload
    @staticmethod
    def from_matrix(matrix: onp.ToFloatND, *, assume_valid: bool = False) -> Rotation: ...

    #
    @overload
    @staticmethod
    def from_mrp(mrp: _ToFloatStrictND) -> Rotation: ...
    @overload
    @staticmethod
    def from_mrp(mrp: onp.ToFloatStrict1D) -> Rotation[tuple[()]]: ...
    @overload
    @staticmethod
    def from_mrp(mrp: onp.ToFloatStrict2D) -> Rotation[tuple[int]]: ...
    @overload
    @staticmethod
    def from_mrp(mrp: onp.ToFloatND) -> Rotation: ...

    #
    @overload
    @staticmethod
    def from_quat(quat: _ToFloatStrictND, *, scalar_first: bool = False) -> Rotation: ...
    @overload
    @staticmethod
    def from_quat(quat: onp.ToFloatStrict1D, *, scalar_first: bool = False) -> Rotation[tuple[()]]: ...
    @overload
    @staticmethod
    def from_quat(quat: onp.ToFloatStrict2D, *, scalar_first: bool = False) -> Rotation[tuple[int]]: ...
    @overload
    @staticmethod
    def from_quat(quat: onp.ToFloatND, *, scalar_first: bool = False) -> Rotation: ...

    #
    @overload
    @staticmethod
    def from_rotvec(rotvec: _ToFloatStrictND, degrees: bool = False) -> Rotation: ...
    @overload
    @staticmethod
    def from_rotvec(rotvec: onp.ToFloatStrict1D, degrees: bool = False) -> Rotation[tuple[()]]: ...
    @overload
    @staticmethod
    def from_rotvec(rotvec: onp.ToFloatStrict2D, degrees: bool = False) -> Rotation[tuple[int]]: ...
    @overload
    @staticmethod
    def from_rotvec(rotvec: onp.ToFloatND, degrees: bool = False) -> Rotation: ...

    #
    @overload
    @staticmethod
    def from_euler(seq: str, angles: onp.ToFloat, degrees: bool = False) -> Rotation[tuple[()]]: ...
    @overload
    @staticmethod
    def from_euler[ShapeT: tuple[int, ...]](
        seq: str, angles: onp.ArrayND[npc.floating | npc.integer, ShapeT], degrees: bool = False
    ) -> Rotation[ShapeT]: ...
    @overload
    @staticmethod
    def from_euler(seq: str, angles: onp.ToFloatND, degrees: bool = False) -> Rotation: ...

    #
    @overload
    @staticmethod
    def from_davenport(
        axes: onp.ToFloatND, order: _RotOrder, angles: onp.ToFloat, degrees: bool = False
    ) -> Rotation[tuple[()]]: ...
    @overload
    @staticmethod
    def from_davenport[ShapeT: tuple[int, ...]](
        axes: onp.ToFloatND, order: _RotOrder, angles: onp.ArrayND[npc.floating | npc.integer, ShapeT], degrees: bool = False
    ) -> Rotation[ShapeT]: ...
    @overload
    @staticmethod
    def from_davenport(axes: onp.ToFloatND, order: _RotOrder, angles: onp.ToFloatND, degrees: bool = False) -> Rotation: ...

    #
    @overload
    @staticmethod
    def concatenate[RotationT: Rotation](rotations: RotationT) -> RotationT: ...
    @overload
    @staticmethod
    def concatenate(rotations: Sequence[Rotation[_JustAnyShape]]) -> Rotation: ...
    @overload
    @staticmethod
    def concatenate(rotations: Sequence[Rotation[tuple[()]]]) -> Rotation[tuple[int]]: ...
    @overload
    @staticmethod
    def concatenate(rotations: Sequence[Rotation[tuple[int]]]) -> Rotation[tuple[int, int]]: ...
    @overload
    @staticmethod
    def concatenate(rotations: Sequence[Rotation]) -> Rotation: ...  # pyright: ignore[reportOverlappingOverload]

    #
    @classmethod
    def create_group(cls, group: _RotGroup, axis: _RotAxis = "Z") -> Rotation[tuple[()]]: ...

    #
    @overload
    @staticmethod
    def identity(num: None = None, *, shape: None = None) -> Rotation[tuple[()]]: ...
    @overload
    @staticmethod
    def identity(num: int, *, shape: None = None) -> Rotation[tuple[int]]: ...
    @overload
    @staticmethod
    def identity(num: None = None, *, shape: int) -> Rotation[tuple[int]]: ...
    @overload
    @staticmethod
    def identity[ShapeT: tuple[int, ...]](num: None = None, *, shape: ShapeT) -> Rotation[ShapeT]: ...

    #
    @overload
    @staticmethod
    def random(
        num: None = None, rng: onp.random.ToRNG | None = None, *, shape: None = None, random_state: onp.random.ToRNG | None = None
    ) -> Rotation[tuple[()]]: ...
    @overload
    @staticmethod
    def random(
        num: int, rng: onp.random.ToRNG | None = None, *, shape: None = None, random_state: onp.random.ToRNG | None = None
    ) -> Rotation[tuple[int]]: ...
    @overload
    @staticmethod
    def random[ShapeT: tuple[int, ...]](
        num: None = None, rng: onp.random.ToRNG | None = None, *, shape: ShapeT, random_state: onp.random.ToRNG | None = None
    ) -> Rotation[ShapeT]: ...

    #
    @overload
    @staticmethod
    def align_vectors(
        a: onp.ToFloat1D | onp.ToFloat2D,
        b: onp.ToFloat1D | onp.ToFloat2D,
        weights: onp.ToFloat1D | None = None,
        return_sensitivity: L[False] = False,
    ) -> tuple[Rotation[tuple[()]], float]: ...
    @overload
    @staticmethod
    def align_vectors(
        a: onp.ToFloat1D | onp.ToFloat2D,
        b: onp.ToFloat1D | onp.ToFloat2D,
        weights: onp.ToFloat1D | None = None,
        *,
        return_sensitivity: L[True],
    ) -> tuple[Rotation[tuple[()]], float, onp.Array2D[np.float64]]: ...
    @overload
    @staticmethod
    def align_vectors(
        a: onp.ToFloat1D | onp.ToFloat2D,
        b: onp.ToFloat1D | onp.ToFloat2D,
        weights: onp.ToFloat1D | None,
        return_sensitivity: L[True],
    ) -> tuple[Rotation[tuple[()]], float, onp.Array2D[np.float64]]: ...

class Slerp:
    times: Final[onp.Array1D[np.float64]]
    timedelta: Final[onp.Array1D[np.float64]]
    rotations: Final[Rotation[tuple[int]]]
    rotvecs: Final[onp.Array2D[np.float64]]

    def __init__(self, /, times: onp.ToFloat1D, rotations: Rotation[tuple[int]]) -> None: ...

    #
    @overload
    def __call__(self, /, times: onp.ToFloat) -> Rotation[tuple[()]]: ...
    @overload
    def __call__(self, /, times: onp.ToFloat1D) -> Rotation[tuple[int]]: ...
