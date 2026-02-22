from collections.abc import Iterator, Sequence
from types import EllipsisType, ModuleType
from typing import Final, Literal as L, Self, TypeAlias, overload
from typing_extensions import TypeVar

import numpy as np
import optype as op
import optype.numpy as onp

_Float1D: TypeAlias = onp.Array1D[np.float64]
_Float2D: TypeAlias = onp.Array2D[np.float64]
_Float3D: TypeAlias = onp.Array3D[np.float64]

_Float1D2D: TypeAlias = _Float1D | _Float2D
_Float2D3D: TypeAlias = _Float2D | _Float3D

_RotOrder: TypeAlias = L["e", "extrinsic", "i", "intrinsic"]
_RotGroup: TypeAlias = L["I", "O", "T", "D", "Dn", "C", "Cn"]
_RotAxisSeq: TypeAlias = L[
    "xyz", "xzy", "yxz", "yzx", "zxy", "zyx",
    "xyx", "xzx", "yxy", "yzy", "zxz", "zyz",
    "XYZ", "XZY", "YXZ", "YZX", "ZXY", "ZYX",
    "XYX", "XZX", "YXY", "YZY", "ZXZ", "ZYZ",
]  # fmt: skip
_RotAxis: TypeAlias = L["X", "Y", "Z"]

_RotationT = TypeVar("_RotationT", bound=Rotation)

###

backend_registry: Final[dict[ModuleType, ModuleType]] = ...  # undocumented

# either `._rotation_cy` or `._rotation_xp`
def select_backend(xp: ModuleType, cython_compatible: bool) -> ModuleType: ...  # undocumented

class Rotation:
    @property
    def single(self, /) -> bool: ...
    @property
    def shape(self, /) -> tuple[int, ...]: ...

    #
    def __init__(
        self, /, quat: onp.ToFloat1D | onp.ToFloat2D, normalize: bool = True, copy: bool = True, scalar_first: bool = False
    ) -> None: ...

    #
    def __bool__(self, /) -> L[True]: ...
    def __len__(self, /) -> int: ...
    def __getitem__(self, key: op.CanIndex | onp.SequenceND[op.CanIndex] | onp.ToIntND | slice, /) -> Rotation: ...
    def __setitem__(self, indexer: int | slice | EllipsisType | None, value: Rotation, /) -> None: ...
    def __iter__(self) -> Iterator[Rotation]: ...
    def __mul__(self, other: Rotation, /) -> Rotation: ...
    def __pow__(self, n: onp.ToFloat, /, modulus: None = None) -> Rotation: ...
    def inv(self, /) -> Self: ...
    def mean(self, /, weights: onp.ToFloat1D | None = None, axis: int | tuple[int, ...] | None = None) -> Rotation: ...

    #
    def as_matrix(self, /) -> _Float2D3D: ...
    def as_mrp(self, /) -> _Float1D2D: ...
    def as_rotvec(self, /, degrees: bool = False) -> _Float1D2D: ...
    def as_euler(self, /, seq: _RotAxisSeq, degrees: bool = False, *, suppress_warnings: bool = False) -> _Float1D2D: ...
    def as_quat(self, /, canonical: bool = False, *, scalar_first: bool = False) -> _Float1D2D: ...
    def as_davenport(
        self, /, axes: onp.ToFloat1D | onp.ToFloat2D, order: _RotOrder, degrees: bool = False, *, suppress_warnings: bool = False
    ) -> _Float1D2D: ...

    #
    def magnitude(self, /) -> float | np.float64 | _Float1D: ...
    def apply(self, /, vectors: onp.ToFloat1D | onp.ToFloat2D, inverse: bool = False) -> _Float1D2D: ...
    def approx_equal(
        self, /, other: Rotation, atol: onp.ToFloat | None = None, degrees: bool = False
    ) -> bool | np.bool_ | onp.Array1D[np.bool_]: ...

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
    @staticmethod
    def from_matrix(matrix: onp.ToFloat2D | onp.ToFloat3D, *, assume_valid: bool = False) -> Rotation: ...
    @staticmethod
    def from_mrp(mrp: onp.ToFloat1D | onp.ToFloat2D) -> Rotation: ...
    @staticmethod
    def from_quat(quat: onp.ToFloat1D | onp.ToFloat2D, *, scalar_first: bool = False) -> Rotation: ...
    @staticmethod
    def from_rotvec(rotvec: onp.ToFloat1D | onp.ToFloat2D, degrees: bool = False) -> Rotation: ...
    @staticmethod
    def from_euler(seq: str, angles: onp.ToFloat | onp.ToFloat1D | onp.ToFloat2D, degrees: bool = False) -> Rotation: ...
    @staticmethod
    def from_davenport(
        axes: onp.ToFloat1D | onp.ToFloat2D,
        order: _RotOrder,
        angles: onp.ToFloat | onp.ToFloat1D | onp.ToFloat2D,
        degrees: bool = False,
    ) -> Rotation: ...

    #
    @overload
    @staticmethod
    def concatenate(rotations: _RotationT) -> _RotationT: ...
    @overload
    @staticmethod
    def concatenate(rotations: Sequence[Rotation]) -> Rotation: ...

    #
    @classmethod
    def create_group(cls, group: _RotGroup, axis: _RotAxis = "Z") -> Self: ...
    @staticmethod
    def identity(num: int | None = None, *, shape: int | tuple[int, ...] | None = None) -> Rotation: ...
    @staticmethod
    def random(
        num: int | None = None,
        rng: onp.random.ToRNG | None = None,
        *,
        shape: tuple[int, ...] | None = None,
        random_state: onp.random.ToRNG | None = None,
    ) -> Rotation: ...

    #
    @overload
    @staticmethod
    def align_vectors(
        a: onp.ToFloat1D | onp.ToFloat2D,
        b: onp.ToFloat1D | onp.ToFloat2D,
        weights: onp.ToFloat1D | None = None,
        return_sensitivity: L[False] = False,
    ) -> tuple[Rotation, float]: ...
    @overload
    @staticmethod
    def align_vectors(
        a: onp.ToFloat1D | onp.ToFloat2D,
        b: onp.ToFloat1D | onp.ToFloat2D,
        weights: onp.ToFloat1D | None = None,
        *,
        return_sensitivity: L[True],
    ) -> tuple[Rotation, float, _Float2D]: ...
    @overload
    @staticmethod
    def align_vectors(
        a: onp.ToFloat1D | onp.ToFloat2D,
        b: onp.ToFloat1D | onp.ToFloat2D,
        weights: onp.ToFloat1D | None,
        return_sensitivity: L[True],
    ) -> tuple[Rotation, float, _Float2D]: ...

class Slerp:
    times: Final[onp.Array1D[np.int32 | np.int64 | np.float32 | np.float64]]
    timedelta: Final[onp.Array1D[np.int32 | np.int64 | np.float32 | np.float64]]
    rotations: Final[Rotation]
    rotvecs: Final[_Float2D]

    def __init__(self, /, times: onp.ToFloat1D, rotations: Rotation) -> None: ...
    def __call__(self, /, times: onp.ToFloat1D) -> Rotation: ...
