# https://github.com/scipy/scipy/blob/maintenance/1.16.x/scipy/spatial/transform/_rotation.pyx

from collections.abc import Sequence
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

class Rotation:
    @property
    def single(self, /) -> bool: ...

    #
    def __init__(
        self, /, quat: onp.ToFloat1D | onp.ToFloat2D, normalize: bool = True, copy: bool = False
    ) -> None: ...  # undocumented
    def __setstate_cython__(self, pyx_state: object, /) -> None: ...  # undocumented
    def __reduce_cython__(self, /) -> None: ...  # undocumented

    #
    def __len__(self, /) -> int: ...
    def __getitem__(self, key: op.CanIndex | onp.SequenceND[op.CanIndex] | onp.ToIntND | slice, /) -> Rotation: ...
    def __mul__(self, other: Rotation, /) -> Rotation: ...
    def __pow__(self, n: onp.ToFloat, /, modulus: None) -> Rotation: ...
    def inv(self, /) -> Self: ...
    def mean(self, /, weights: onp.ToFloat1D | None = None) -> Self: ...

    #
    def as_matrix(self, /) -> _Float2D3D: ...
    def as_mrp(self, /) -> _Float1D2D: ...
    def as_rotvec(self, /, degrees: bool = False) -> _Float1D2D: ...
    def as_euler(self, /, seq: _RotAxisSeq, degrees: bool = False) -> _Float1D2D: ...
    def as_quat(self, /, canonical: bool = False, *, scalar_first: bool = False) -> _Float1D2D: ...
    def as_davenport(self, /, axes: onp.ToFloat1D | onp.ToFloat2D, order: _RotOrder, degrees: bool = False) -> _Float1D2D: ...

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
    @classmethod
    def from_matrix(cls, matrix: onp.ToFloat2D | onp.ToFloat3D) -> Self: ...
    @classmethod
    def from_mrp(cls, mrp: onp.ToFloat1D | onp.ToFloat2D) -> Self: ...
    @classmethod
    def from_quat(cls, quat: onp.ToFloat1D | onp.ToFloat2D, *, scalar_first: bool = False) -> Self: ...
    @classmethod
    def from_rotvec(cls, rotvec: onp.ToFloat1D | onp.ToFloat2D, degrees: bool = False) -> Self: ...
    @classmethod
    def from_euler(cls, seq: str, angles: onp.ToFloat | onp.ToFloat1D | onp.ToFloat2D, degrees: bool = False) -> Self: ...
    @classmethod
    def from_davenport(
        cls,
        axes: onp.ToFloat1D | onp.ToFloat2D,
        order: _RotOrder,
        angles: onp.ToFloat | onp.ToFloat1D | onp.ToFloat2D,
        degrees: bool = False,
    ) -> Self: ...

    #
    @overload
    @classmethod
    def concatenate(cls, rotations: _RotationT) -> _RotationT: ...
    @overload
    @classmethod
    def concatenate(cls, rotations: Sequence[Rotation]) -> Self: ...

    #
    @classmethod
    def create_group(cls, group: _RotGroup, axis: _RotAxis = "Z") -> Self: ...
    @classmethod
    def identity(cls, num: int | None = None) -> Self: ...
    @classmethod
    def random(
        cls, num: int | None = None, rng: onp.random.ToRNG | None = None, *, random_state: onp.random.ToRNG | None = None
    ) -> Self: ...
    #
    @overload
    @classmethod
    def align_vectors(
        cls,
        a: onp.ToFloat1D | onp.ToFloat2D,
        b: onp.ToFloat1D | onp.ToFloat2D,
        weights: onp.ToFloat1D | None = None,
        return_sensitivity: L[False] = False,
    ) -> tuple[Self, float]: ...
    @overload
    @classmethod
    def align_vectors(
        cls,
        a: onp.ToFloat1D | onp.ToFloat2D,
        b: onp.ToFloat1D | onp.ToFloat2D,
        weights: onp.ToFloat1D | None = None,
        *,
        return_sensitivity: L[True],
    ) -> tuple[Self, float, _Float2D]: ...
    @overload
    @classmethod
    def align_vectors(
        cls,
        a: onp.ToFloat1D | onp.ToFloat2D,
        b: onp.ToFloat1D | onp.ToFloat2D,
        weights: onp.ToFloat1D | None,
        return_sensitivity: L[True],
    ) -> tuple[Self, float, _Float2D]: ...

class Slerp:
    times: Final[onp.Array1D[np.int32 | np.int64 | np.float32 | np.float64]]
    timedelta: Final[onp.Array1D[np.int32 | np.int64 | np.float32 | np.float64]]
    rotations: Final[Rotation]
    rotvecs: Final[_Float2D]

    def __init__(self, /, times: onp.ToFloat1D, rotations: Rotation) -> None: ...
    def __call__(self, /, times: onp.ToFloat1D) -> Rotation: ...

# (double[:, :], double[:, :]) -> noexcept double[:, :]
def compose_quat(p: onp.ArrayND[np.float64], q: onp.ArrayND[np.float64]) -> onp.Array2D[np.float64]: ...  # undocumented
