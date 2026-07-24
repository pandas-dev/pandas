from collections.abc import Iterable
from typing import Any, Literal, SupportsIndex, overload
from typing_extensions import TypeVar

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc

__all__ = [
    "binary_closing",
    "binary_dilation",
    "binary_erosion",
    "binary_fill_holes",
    "binary_hit_or_miss",
    "binary_opening",
    "binary_propagation",
    "black_tophat",
    "distance_transform_bf",
    "distance_transform_cdt",
    "distance_transform_edt",
    "generate_binary_structure",
    "grey_closing",
    "grey_dilation",
    "grey_erosion",
    "grey_opening",
    "iterate_structure",
    "morphological_gradient",
    "morphological_laplace",
    "white_tophat",
]

###

type _Mode = Literal["reflect", "constant", "nearest", "mirror", "wrap"]
type _Metric1 = Literal["euclidean"]
type _Metric2 = Literal["taxicab", "cityblock", "manhattan"]
type _Metric3 = Literal["chessboard"]

type _Origin = int | tuple[int, ...]

_ShapeT = TypeVar("_ShapeT", bound=tuple[int, ...], default=tuple[Any, ...])
_OriginScalarT = TypeVar("_OriginScalarT", bound=int | npc.integer, default=int)

###

# NOTE: On numpy<2.1, pyright reports 2 false positive incompatible overload errors here.
# pyright: reportOverlappingOverload=false

@overload  # known shape
def iterate_structure(
    structure: onp.ArrayND[np.bool | npc.integer, _ShapeT], iterations: int, origin: None = None
) -> onp.ArrayND[np.bool, _ShapeT]: ...
@overload  # known shape, origin=<given>
def iterate_structure(
    structure: onp.ArrayND[np.bool | npc.integer, _ShapeT], iterations: int, origin: _OriginScalarT | Iterable[_OriginScalarT]
) -> tuple[onp.ArrayND[np.bool, _ShapeT], list[_OriginScalarT]]: ...
@overload  # unknown shape
def iterate_structure(structure: onp.ToIntND, iterations: int, origin: None = None) -> onp.ArrayND[np.bool]: ...
@overload  # unknown shape, origin=<given>
def iterate_structure(
    structure: onp.ToIntND, iterations: int, origin: _OriginScalarT | Iterable[_OriginScalarT]
) -> tuple[onp.ArrayND[np.bool], list[_OriginScalarT]]: ...

#
@overload
def generate_binary_structure(rank: Literal[0, -1, -2, -3], connectivity: int) -> onp.Array0D[np.bool]: ...
@overload
def generate_binary_structure(rank: Literal[1], connectivity: int) -> onp.Array1D[np.bool]: ...
@overload
def generate_binary_structure(rank: Literal[2], connectivity: int) -> onp.Array2D[np.bool]: ...
@overload
def generate_binary_structure(rank: Literal[3], connectivity: int) -> onp.Array3D[np.bool]: ...
@overload
def generate_binary_structure(rank: int, connectivity: int) -> onp.ArrayND[np.bool]: ...

#
@overload
def binary_erosion(
    input: onp.ToFloatND,
    structure: onp.ToIntND | None = None,
    iterations: SupportsIndex = 1,
    mask: onp.ToIntND | None = None,
    output: None = None,
    border_value: int = 0,
    origin: _Origin = 0,
    brute_force: bool = False,
    *,
    axes: tuple[int, ...] | None = None,
) -> onp.ArrayND[np.bool]: ...
@overload
def binary_erosion[OutputArrayT: onp.ArrayND[np.bool | npc.integer | npc.floating]](
    input: onp.ToFloatND,
    structure: onp.ToIntND | None = None,
    iterations: SupportsIndex = 1,
    mask: onp.ToIntND | None = None,
    *,
    output: OutputArrayT,
    border_value: int = 0,
    origin: _Origin = 0,
    brute_force: bool = False,
    axes: tuple[int, ...] | None = None,
) -> OutputArrayT: ...

# keep in sync with `binary_erosion`
@overload
def binary_dilation(
    input: onp.ToFloatND,
    structure: onp.ToIntND | None = None,
    iterations: SupportsIndex = 1,
    mask: onp.ToIntND | None = None,
    output: None = None,
    border_value: int = 0,
    origin: _Origin = 0,
    brute_force: bool = False,
    *,
    axes: tuple[int, ...] | None = None,
) -> onp.ArrayND[np.bool]: ...
@overload
def binary_dilation[OutputArrayT: onp.ArrayND[np.bool | npc.integer | npc.floating]](
    input: onp.ToFloatND,
    structure: onp.ToIntND | None = None,
    iterations: SupportsIndex = 1,
    mask: onp.ToIntND | None = None,
    *,
    output: OutputArrayT,
    border_value: int = 0,
    origin: _Origin = 0,
    brute_force: bool = False,
    axes: tuple[int, ...] | None = None,
) -> OutputArrayT: ...

# keep in sync with `binary_erosion` (but with shuffled `mask`, `output` and `origin`)
@overload
def binary_opening(
    input: onp.ToFloatND,
    structure: onp.ToIntND | None = None,
    iterations: SupportsIndex = 1,
    output: None = None,
    origin: _Origin = 0,
    mask: onp.ToIntND | None = None,
    border_value: int = 0,
    brute_force: bool = False,
    *,
    axes: tuple[int, ...] | None = None,
) -> onp.ArrayND[np.bool]: ...
@overload
def binary_opening[OutputArrayT: onp.ArrayND[np.bool | npc.integer | npc.floating]](
    input: onp.ToFloatND,
    structure: onp.ToIntND | None = None,
    iterations: SupportsIndex = 1,
    *,
    output: OutputArrayT,
    origin: _Origin = 0,
    mask: onp.ToIntND | None = None,
    border_value: int = 0,
    brute_force: bool = False,
    axes: tuple[int, ...] | None = None,
) -> OutputArrayT: ...

# keep in sync with `binary_erosion`
@overload
def binary_closing(
    input: onp.ToFloatND,
    structure: onp.ToIntND | None = None,
    iterations: SupportsIndex = 1,
    output: None = None,
    origin: _Origin = 0,
    mask: onp.ToIntND | None = None,
    border_value: int = 0,
    brute_force: bool = False,
    *,
    axes: tuple[int, ...] | None = None,
) -> onp.ArrayND[np.bool]: ...
@overload
def binary_closing[OutputArrayT: onp.ArrayND[np.bool | npc.integer | npc.floating]](
    input: onp.ToFloatND,
    structure: onp.ToIntND | None = None,
    iterations: SupportsIndex = 1,
    *,
    output: OutputArrayT,
    origin: _Origin = 0,
    mask: onp.ToIntND | None = None,
    border_value: int = 0,
    brute_force: bool = False,
    axes: tuple[int, ...] | None = None,
) -> OutputArrayT: ...

# keep roughly in sync with `binary_erosion`
@overload
def binary_hit_or_miss(
    input: onp.ToFloatND,
    structure1: onp.ToIntND | None = None,
    structure2: onp.ToIntND | None = None,
    output: None = None,
    origin1: _Origin = 0,
    origin2: _Origin | None = None,
    *,
    axes: tuple[int, ...] | None = None,
) -> onp.ArrayND[np.bool]: ...
@overload
def binary_hit_or_miss[OutputArrayT: onp.ArrayND[np.bool | npc.integer | npc.floating]](
    input: onp.ToFloatND,
    structure1: onp.ToIntND | None = None,
    structure2: onp.ToIntND | None = None,
    *,
    output: OutputArrayT,
    origin1: _Origin = 0,
    origin2: _Origin | None = None,
    axes: tuple[int, ...] | None = None,
) -> OutputArrayT: ...

# keep roughly in sync with `binary_erosion`
@overload
def binary_propagation(
    input: onp.ToFloatND,
    structure: onp.ToIntND | None = None,
    mask: onp.ToIntND | None = None,
    output: None = None,
    border_value: int = 0,
    origin: _Origin = 0,
    *,
    axes: tuple[int, ...] | None = None,
) -> onp.ArrayND[np.bool]: ...
@overload
def binary_propagation[OutputArrayT: onp.ArrayND[np.bool | npc.integer | npc.floating]](
    input: onp.ToFloatND,
    structure: onp.ToIntND | None = None,
    mask: onp.ToIntND | None = None,
    *,
    output: OutputArrayT,
    border_value: int = 0,
    origin: _Origin = 0,
    axes: tuple[int, ...] | None = None,
) -> OutputArrayT: ...

# keep roughly in sync with `binary_erosion`
@overload
def binary_fill_holes(
    input: onp.ToFloatND,
    structure: onp.ToIntND | None = None,
    output: None = None,
    origin: _Origin = 0,
    *,
    axes: tuple[int, ...] | None = None,
) -> onp.ArrayND[np.bool]: ...
@overload
def binary_fill_holes[OutputArrayT: onp.ArrayND[np.bool | npc.integer | npc.floating]](
    input: onp.ToFloatND,
    structure: onp.ToIntND | None = None,
    *,
    output: OutputArrayT,
    origin: _Origin = 0,
    axes: tuple[int, ...] | None = None,
) -> OutputArrayT: ...

#
@overload
def grey_erosion[OutputArrayT: onp.ArrayND[np.bool | npc.integer | npc.floating]](
    input: OutputArrayT,
    size: tuple[int, ...] | None = None,
    footprint: onp.ToIntND | None = None,
    structure: onp.ToIntND | None = None,
    output: None = None,
    mode: _Mode = "reflect",
    cval: float = 0.0,
    origin: _Origin = 0,
    *,
    axes: tuple[int, ...] | None = None,
) -> OutputArrayT: ...
@overload
def grey_erosion[OutputArrayT: onp.ArrayND[np.bool | npc.integer | npc.floating]](
    input: onp.ToFloatND,
    size: tuple[int, ...] | None = None,
    footprint: onp.ToIntND | None = None,
    structure: onp.ToIntND | None = None,
    *,
    output: OutputArrayT,
    mode: _Mode = "reflect",
    cval: float = 0.0,
    origin: _Origin = 0,
    axes: tuple[int, ...] | None = None,
) -> OutputArrayT: ...
@overload
def grey_erosion(
    input: onp.SequenceND[bool],
    size: tuple[int, ...] | None = None,
    footprint: onp.ToIntND | None = None,
    structure: onp.ToIntND | None = None,
    output: None = None,
    mode: _Mode = "reflect",
    cval: float = 0.0,
    origin: _Origin = 0,
    *,
    axes: tuple[int, ...] | None = None,
) -> onp.ArrayND[np.bool]: ...
@overload
def grey_erosion(
    input: onp.SequenceND[list[int]] | list[int],
    size: tuple[int, ...] | None = None,
    footprint: onp.ToIntND | None = None,
    structure: onp.ToIntND | None = None,
    output: None = None,
    mode: _Mode = "reflect",
    cval: float = 0.0,
    origin: _Origin = 0,
    *,
    axes: tuple[int, ...] | None = None,
) -> onp.ArrayND[np.int_]: ...
@overload
def grey_erosion(
    input: onp.SequenceND[list[float]] | list[float],
    size: tuple[int, ...] | None = None,
    footprint: onp.ToIntND | None = None,
    structure: onp.ToIntND | None = None,
    output: None = None,
    mode: _Mode = "reflect",
    cval: float = 0.0,
    origin: _Origin = 0,
    *,
    axes: tuple[int, ...] | None = None,
) -> onp.ArrayND[np.float64]: ...

# keep in sync with `grey_erosion`
@overload
def grey_dilation[OutputArrayT: onp.ArrayND[np.bool | npc.integer | npc.floating]](
    input: OutputArrayT,
    size: tuple[int, ...] | None = None,
    footprint: onp.ToIntND | None = None,
    structure: onp.ToIntND | None = None,
    output: None = None,
    mode: _Mode = "reflect",
    cval: float = 0.0,
    origin: _Origin = 0,
    *,
    axes: tuple[int, ...] | None = None,
) -> OutputArrayT: ...
@overload
def grey_dilation[OutputArrayT: onp.ArrayND[np.bool | npc.integer | npc.floating]](
    input: onp.ToFloatND,
    size: tuple[int, ...] | None = None,
    footprint: onp.ToIntND | None = None,
    structure: onp.ToIntND | None = None,
    *,
    output: OutputArrayT,
    mode: _Mode = "reflect",
    cval: float = 0.0,
    origin: _Origin = 0,
    axes: tuple[int, ...] | None = None,
) -> OutputArrayT: ...
@overload
def grey_dilation(
    input: onp.SequenceND[bool],
    size: tuple[int, ...] | None = None,
    footprint: onp.ToIntND | None = None,
    structure: onp.ToIntND | None = None,
    output: None = None,
    mode: _Mode = "reflect",
    cval: float = 0.0,
    origin: _Origin = 0,
    *,
    axes: tuple[int, ...] | None = None,
) -> onp.ArrayND[np.bool]: ...
@overload
def grey_dilation(
    input: onp.SequenceND[list[int]] | list[int],
    size: tuple[int, ...] | None = None,
    footprint: onp.ToIntND | None = None,
    structure: onp.ToIntND | None = None,
    output: None = None,
    mode: _Mode = "reflect",
    cval: float = 0.0,
    origin: _Origin = 0,
    *,
    axes: tuple[int, ...] | None = None,
) -> onp.ArrayND[np.int_]: ...
@overload
def grey_dilation(
    input: onp.SequenceND[list[float]] | list[float],
    size: tuple[int, ...] | None = None,
    footprint: onp.ToIntND | None = None,
    structure: onp.ToIntND | None = None,
    output: None = None,
    mode: _Mode = "reflect",
    cval: float = 0.0,
    origin: _Origin = 0,
    *,
    axes: tuple[int, ...] | None = None,
) -> onp.ArrayND[np.float64]: ...

# keep in sync with `grey_erosion`
@overload
def grey_opening[OutputArrayT: onp.ArrayND[np.bool | npc.integer | npc.floating]](
    input: OutputArrayT,
    size: tuple[int, ...] | None = None,
    footprint: onp.ToIntND | None = None,
    structure: onp.ToIntND | None = None,
    output: None = None,
    mode: _Mode = "reflect",
    cval: float = 0.0,
    origin: _Origin = 0,
    *,
    axes: tuple[int, ...] | None = None,
) -> OutputArrayT: ...
@overload
def grey_opening[OutputArrayT: onp.ArrayND[np.bool | npc.integer | npc.floating]](
    input: onp.ToFloatND,
    size: tuple[int, ...] | None = None,
    footprint: onp.ToIntND | None = None,
    structure: onp.ToIntND | None = None,
    *,
    output: OutputArrayT,
    mode: _Mode = "reflect",
    cval: float = 0.0,
    origin: _Origin = 0,
    axes: tuple[int, ...] | None = None,
) -> OutputArrayT: ...
@overload
def grey_opening(
    input: onp.SequenceND[bool],
    size: tuple[int, ...] | None = None,
    footprint: onp.ToIntND | None = None,
    structure: onp.ToIntND | None = None,
    output: None = None,
    mode: _Mode = "reflect",
    cval: float = 0.0,
    origin: _Origin = 0,
    *,
    axes: tuple[int, ...] | None = None,
) -> onp.ArrayND[np.bool]: ...
@overload
def grey_opening(
    input: onp.SequenceND[list[int]] | list[int],
    size: tuple[int, ...] | None = None,
    footprint: onp.ToIntND | None = None,
    structure: onp.ToIntND | None = None,
    output: None = None,
    mode: _Mode = "reflect",
    cval: float = 0.0,
    origin: _Origin = 0,
    *,
    axes: tuple[int, ...] | None = None,
) -> onp.ArrayND[np.int_]: ...
@overload
def grey_opening(
    input: onp.SequenceND[list[float]] | list[float],
    size: tuple[int, ...] | None = None,
    footprint: onp.ToIntND | None = None,
    structure: onp.ToIntND | None = None,
    output: None = None,
    mode: _Mode = "reflect",
    cval: float = 0.0,
    origin: _Origin = 0,
    *,
    axes: tuple[int, ...] | None = None,
) -> onp.ArrayND[np.float64]: ...

# keep in sync with `grey_erosion`
@overload
def grey_closing[OutputArrayT: onp.ArrayND[np.bool | npc.integer | npc.floating]](
    input: OutputArrayT,
    size: tuple[int, ...] | None = None,
    footprint: onp.ToIntND | None = None,
    structure: onp.ToIntND | None = None,
    output: None = None,
    mode: _Mode = "reflect",
    cval: float = 0.0,
    origin: _Origin = 0,
    *,
    axes: tuple[int, ...] | None = None,
) -> OutputArrayT: ...
@overload
def grey_closing[OutputArrayT: onp.ArrayND[np.bool | npc.integer | npc.floating]](
    input: onp.ToFloatND,
    size: tuple[int, ...] | None = None,
    footprint: onp.ToIntND | None = None,
    structure: onp.ToIntND | None = None,
    *,
    output: OutputArrayT,
    mode: _Mode = "reflect",
    cval: float = 0.0,
    origin: _Origin = 0,
    axes: tuple[int, ...] | None = None,
) -> OutputArrayT: ...
@overload
def grey_closing(
    input: onp.SequenceND[bool],
    size: tuple[int, ...] | None = None,
    footprint: onp.ToIntND | None = None,
    structure: onp.ToIntND | None = None,
    output: None = None,
    mode: _Mode = "reflect",
    cval: float = 0.0,
    origin: _Origin = 0,
    *,
    axes: tuple[int, ...] | None = None,
) -> onp.ArrayND[np.bool]: ...
@overload
def grey_closing(
    input: onp.SequenceND[list[int]] | list[int],
    size: tuple[int, ...] | None = None,
    footprint: onp.ToIntND | None = None,
    structure: onp.ToIntND | None = None,
    output: None = None,
    mode: _Mode = "reflect",
    cval: float = 0.0,
    origin: _Origin = 0,
    *,
    axes: tuple[int, ...] | None = None,
) -> onp.ArrayND[np.int_]: ...
@overload
def grey_closing(
    input: onp.SequenceND[list[float]] | list[float],
    size: tuple[int, ...] | None = None,
    footprint: onp.ToIntND | None = None,
    structure: onp.ToIntND | None = None,
    output: None = None,
    mode: _Mode = "reflect",
    cval: float = 0.0,
    origin: _Origin = 0,
    *,
    axes: tuple[int, ...] | None = None,
) -> onp.ArrayND[np.float64]: ...

# keep in sync with `grey_erosion`
@overload
def morphological_gradient[OutputArrayT: onp.ArrayND[np.bool | npc.integer | npc.floating]](
    input: OutputArrayT,
    size: tuple[int, ...] | None = None,
    footprint: onp.ToIntND | None = None,
    structure: onp.ToIntND | None = None,
    output: None = None,
    mode: _Mode = "reflect",
    cval: float = 0.0,
    origin: _Origin = 0,
    *,
    axes: tuple[int, ...] | None = None,
) -> OutputArrayT: ...
@overload
def morphological_gradient[OutputArrayT: onp.ArrayND[np.bool | npc.integer | npc.floating]](
    input: onp.ToFloatND,
    size: tuple[int, ...] | None = None,
    footprint: onp.ToIntND | None = None,
    structure: onp.ToIntND | None = None,
    *,
    output: OutputArrayT,
    mode: _Mode = "reflect",
    cval: float = 0.0,
    origin: _Origin = 0,
    axes: tuple[int, ...] | None = None,
) -> OutputArrayT: ...
@overload
def morphological_gradient(
    input: onp.SequenceND[bool],
    size: tuple[int, ...] | None = None,
    footprint: onp.ToIntND | None = None,
    structure: onp.ToIntND | None = None,
    output: None = None,
    mode: _Mode = "reflect",
    cval: float = 0.0,
    origin: _Origin = 0,
    *,
    axes: tuple[int, ...] | None = None,
) -> onp.ArrayND[np.bool]: ...
@overload
def morphological_gradient(
    input: onp.SequenceND[list[int]] | list[int],
    size: tuple[int, ...] | None = None,
    footprint: onp.ToIntND | None = None,
    structure: onp.ToIntND | None = None,
    output: None = None,
    mode: _Mode = "reflect",
    cval: float = 0.0,
    origin: _Origin = 0,
    *,
    axes: tuple[int, ...] | None = None,
) -> onp.ArrayND[np.int_]: ...
@overload
def morphological_gradient(
    input: onp.SequenceND[list[float]] | list[float],
    size: tuple[int, ...] | None = None,
    footprint: onp.ToIntND | None = None,
    structure: onp.ToIntND | None = None,
    output: None = None,
    mode: _Mode = "reflect",
    cval: float = 0.0,
    origin: _Origin = 0,
    *,
    axes: tuple[int, ...] | None = None,
) -> onp.ArrayND[np.float64]: ...

# keep in sync with `grey_erosion`
@overload
def morphological_laplace[OutputArrayT: onp.ArrayND[np.bool | npc.integer | npc.floating]](
    input: OutputArrayT,
    size: tuple[int, ...] | None = None,
    footprint: onp.ToIntND | None = None,
    structure: onp.ToIntND | None = None,
    output: None = None,
    mode: _Mode = "reflect",
    cval: float = 0.0,
    origin: _Origin = 0,
    *,
    axes: tuple[int, ...] | None = None,
) -> OutputArrayT: ...
@overload
def morphological_laplace[OutputArrayT: onp.ArrayND[np.bool | npc.integer | npc.floating]](
    input: onp.ToFloatND,
    size: tuple[int, ...] | None = None,
    footprint: onp.ToIntND | None = None,
    structure: onp.ToIntND | None = None,
    *,
    output: OutputArrayT,
    mode: _Mode = "reflect",
    cval: float = 0.0,
    origin: _Origin = 0,
    axes: tuple[int, ...] | None = None,
) -> OutputArrayT: ...
@overload
def morphological_laplace(
    input: onp.SequenceND[bool],
    size: tuple[int, ...] | None = None,
    footprint: onp.ToIntND | None = None,
    structure: onp.ToIntND | None = None,
    output: None = None,
    mode: _Mode = "reflect",
    cval: float = 0.0,
    origin: _Origin = 0,
    *,
    axes: tuple[int, ...] | None = None,
) -> onp.ArrayND[np.bool]: ...
@overload
def morphological_laplace(
    input: onp.SequenceND[list[int]] | list[int],
    size: tuple[int, ...] | None = None,
    footprint: onp.ToIntND | None = None,
    structure: onp.ToIntND | None = None,
    output: None = None,
    mode: _Mode = "reflect",
    cval: float = 0.0,
    origin: _Origin = 0,
    *,
    axes: tuple[int, ...] | None = None,
) -> onp.ArrayND[np.int_]: ...
@overload
def morphological_laplace(
    input: onp.SequenceND[list[float]] | list[float],
    size: tuple[int, ...] | None = None,
    footprint: onp.ToIntND | None = None,
    structure: onp.ToIntND | None = None,
    output: None = None,
    mode: _Mode = "reflect",
    cval: float = 0.0,
    origin: _Origin = 0,
    *,
    axes: tuple[int, ...] | None = None,
) -> onp.ArrayND[np.float64]: ...

# keep in sync with `grey_erosion`
@overload
def white_tophat[OutputArrayT: onp.ArrayND[np.bool | npc.integer | npc.floating]](
    input: OutputArrayT,
    size: tuple[int, ...] | None = None,
    footprint: onp.ToIntND | None = None,
    structure: onp.ToIntND | None = None,
    output: None = None,
    mode: _Mode = "reflect",
    cval: float = 0.0,
    origin: _Origin = 0,
    *,
    axes: tuple[int, ...] | None = None,
) -> OutputArrayT: ...
@overload
def white_tophat[OutputArrayT: onp.ArrayND[np.bool | npc.integer | npc.floating]](
    input: onp.ToFloatND,
    size: tuple[int, ...] | None = None,
    footprint: onp.ToIntND | None = None,
    structure: onp.ToIntND | None = None,
    *,
    output: OutputArrayT,
    mode: _Mode = "reflect",
    cval: float = 0.0,
    origin: _Origin = 0,
    axes: tuple[int, ...] | None = None,
) -> OutputArrayT: ...
@overload
def white_tophat(
    input: onp.SequenceND[bool],
    size: tuple[int, ...] | None = None,
    footprint: onp.ToIntND | None = None,
    structure: onp.ToIntND | None = None,
    output: None = None,
    mode: _Mode = "reflect",
    cval: float = 0.0,
    origin: _Origin = 0,
    *,
    axes: tuple[int, ...] | None = None,
) -> onp.ArrayND[np.bool]: ...
@overload
def white_tophat(
    input: onp.SequenceND[list[int]] | list[int],
    size: tuple[int, ...] | None = None,
    footprint: onp.ToIntND | None = None,
    structure: onp.ToIntND | None = None,
    output: None = None,
    mode: _Mode = "reflect",
    cval: float = 0.0,
    origin: _Origin = 0,
    *,
    axes: tuple[int, ...] | None = None,
) -> onp.ArrayND[np.int_]: ...
@overload
def white_tophat(
    input: onp.SequenceND[list[float]] | list[float],
    size: tuple[int, ...] | None = None,
    footprint: onp.ToIntND | None = None,
    structure: onp.ToIntND | None = None,
    output: None = None,
    mode: _Mode = "reflect",
    cval: float = 0.0,
    origin: _Origin = 0,
    *,
    axes: tuple[int, ...] | None = None,
) -> onp.ArrayND[np.float64]: ...

# keep in sync with `grey_erosion`
@overload
def black_tophat[OutputArrayT: onp.ArrayND[np.bool | npc.integer | npc.floating]](
    input: OutputArrayT,
    size: tuple[int, ...] | None = None,
    footprint: onp.ToIntND | None = None,
    structure: onp.ToIntND | None = None,
    output: None = None,
    mode: _Mode = "reflect",
    cval: float = 0.0,
    origin: _Origin = 0,
    *,
    axes: tuple[int, ...] | None = None,
) -> OutputArrayT: ...
@overload
def black_tophat[OutputArrayT: onp.ArrayND[np.bool | npc.integer | npc.floating]](
    input: onp.ToFloatND,
    size: tuple[int, ...] | None = None,
    footprint: onp.ToIntND | None = None,
    structure: onp.ToIntND | None = None,
    *,
    output: OutputArrayT,
    mode: _Mode = "reflect",
    cval: float = 0.0,
    origin: _Origin = 0,
    axes: tuple[int, ...] | None = None,
) -> OutputArrayT: ...
@overload
def black_tophat(
    input: onp.SequenceND[bool],
    size: tuple[int, ...] | None = None,
    footprint: onp.ToIntND | None = None,
    structure: onp.ToIntND | None = None,
    output: None = None,
    mode: _Mode = "reflect",
    cval: float = 0.0,
    origin: _Origin = 0,
    *,
    axes: tuple[int, ...] | None = None,
) -> onp.ArrayND[np.bool]: ...
@overload
def black_tophat(
    input: onp.SequenceND[list[int]] | list[int],
    size: tuple[int, ...] | None = None,
    footprint: onp.ToIntND | None = None,
    structure: onp.ToIntND | None = None,
    output: None = None,
    mode: _Mode = "reflect",
    cval: float = 0.0,
    origin: _Origin = 0,
    *,
    axes: tuple[int, ...] | None = None,
) -> onp.ArrayND[np.int_]: ...
@overload
def black_tophat(
    input: onp.SequenceND[list[float]] | list[float],
    size: tuple[int, ...] | None = None,
    footprint: onp.ToIntND | None = None,
    structure: onp.ToIntND | None = None,
    output: None = None,
    mode: _Mode = "reflect",
    cval: float = 0.0,
    origin: _Origin = 0,
    *,
    axes: tuple[int, ...] | None = None,
) -> onp.ArrayND[np.float64]: ...

#
@overload  # return_distances=False, return_indices=False (default)
def distance_transform_bf(
    input: onp.ToFloatND,
    metric: _Metric1 | _Metric2 | _Metric3 = "euclidean",
    sampling: onp.ToFloat | onp.ToFloat1D | None = None,
    *,
    return_distances: Literal[False],
    return_indices: Literal[False] = False,
    distances: None = None,
    indices: onp.ArrayND[np.int32] | None = None,
) -> None: ...
@overload  # return_distances=False, return_indices=True
def distance_transform_bf(
    input: onp.ToFloatND,
    metric: _Metric1 | _Metric2 | _Metric3 = "euclidean",
    sampling: onp.ToFloat | onp.ToFloat1D | None = None,
    *,
    return_distances: Literal[False],
    return_indices: Literal[True],
    distances: None = None,
    indices: onp.ArrayND[np.int32] | None = None,
) -> onp.ArrayND[np.int32]: ...
@overload  # metric == "euclidean" (default), distances=<given>, return_indices=False (default)
def distance_transform_bf(
    input: onp.ToFloatND,
    metric: _Metric1 = "euclidean",
    sampling: onp.ToFloat | onp.ToFloat1D | None = None,
    return_distances: bool = True,
    return_indices: Literal[False] = False,
    *,
    distances: onp.ArrayND[np.float64],
    indices: onp.ArrayND[np.int32] | None = None,
) -> None: ...
@overload  # metric == "euclidean" (default), distances=<given>, return_indices=True
def distance_transform_bf(
    input: onp.ToFloatND,
    metric: _Metric1 = "euclidean",
    sampling: onp.ToFloat | onp.ToFloat1D | None = None,
    return_distances: bool = True,
    *,
    return_indices: Literal[True],
    distances: onp.ArrayND[np.float64, _ShapeT],
    indices: None = None,
) -> onp.ArrayND[np.int32, _ShapeT]: ...
@overload  # metric == "euclidean" (default), return_distances=True (default), return_indices=False (default)
def distance_transform_bf(
    input: onp.ToFloatND,
    metric: _Metric1 = "euclidean",
    sampling: onp.ToFloat | onp.ToFloat1D | None = None,
    return_distances: Literal[True] = True,
    return_indices: Literal[False] = False,
    distances: None = None,
    indices: onp.ArrayND[np.int32, _ShapeT] | None = None,
) -> onp.ArrayND[np.float64, _ShapeT]: ...
@overload  # metric == "euclidean" (default), return_distances=True (default), return_indices=True
def distance_transform_bf(
    input: onp.ToFloatND,
    metric: _Metric1 = "euclidean",
    sampling: onp.ToFloat | onp.ToFloat1D | None = None,
    return_distances: Literal[True] = True,
    *,
    return_indices: Literal[True],
    distances: None = None,
    indices: None = None,
) -> tuple[onp.ArrayND[np.float64], onp.ArrayND[np.int32]]: ...
@overload  # metric != "euclidean", distances=<given>, return_indices=False (default)
def distance_transform_bf(
    input: onp.ToComplex | onp.ToComplexND,
    metric: _Metric2 | _Metric3,
    sampling: None = None,
    return_distances: bool = True,
    return_indices: Literal[False] = False,
    *,
    distances: onp.ArrayND[np.uint32],
    indices: onp.ArrayND[np.int32] | None = None,
) -> None: ...
@overload  # metric != "euclidean", distances=<given>, return_indices=True
def distance_transform_bf(
    input: onp.ToComplex | onp.ToComplexND,
    metric: _Metric2 | _Metric3,
    sampling: None = None,
    return_distances: bool = True,
    *,
    return_indices: Literal[True],
    distances: onp.ArrayND[np.uint32, _ShapeT],
    indices: None = None,
) -> onp.ArrayND[np.int32, _ShapeT]: ...
@overload  # metric != "euclidean", return_distances=True (default), return_indices=False (default)
def distance_transform_bf(
    input: onp.ToComplex | onp.ToComplexND,
    metric: _Metric2 | _Metric3,
    sampling: None = None,
    return_distances: Literal[True] = True,
    return_indices: Literal[False] = False,
    distances: None = None,
    indices: onp.ArrayND[np.int32, _ShapeT] | None = None,
) -> onp.ArrayND[np.uint32, _ShapeT]: ...
@overload  # metric != "euclidean", return_distances=True (default), return_indices=True
def distance_transform_bf(
    input: onp.ToComplex | onp.ToComplexND,
    metric: _Metric2 | _Metric3,
    sampling: None = None,
    return_distances: Literal[True] = True,
    *,
    return_indices: Literal[True],
    distances: None = None,
    indices: None = None,
) -> tuple[onp.ArrayND[np.uint32], onp.ArrayND[np.int32]]: ...

#
@overload  # return_distances=False, return_indices=False (default)
def distance_transform_cdt(
    input: onp.ToFloatND,
    metric: _Metric2 | _Metric3 | onp.ToFloatND = "chessboard",
    *,
    return_distances: Literal[False],
    return_indices: Literal[False] = False,
    distances: None = None,
    indices: onp.ArrayND[np.int32] | None = None,
) -> None: ...
@overload  # return_distances=False, return_indices=True
def distance_transform_cdt(
    input: onp.ToFloatND,
    metric: _Metric2 | _Metric3 | onp.ToFloatND = "chessboard",
    *,
    return_distances: Literal[False],
    return_indices: Literal[True],
    distances: None = None,
    indices: None = None,
) -> onp.ArrayND[np.int32]: ...
@overload  # distances=<given>, return_indices=False (default)
def distance_transform_cdt(
    input: onp.ToFloatND,
    metric: _Metric2 | _Metric3 | onp.ToFloatND = "chessboard",
    return_distances: bool = True,
    return_indices: Literal[False] = False,
    *,
    distances: onp.ArrayND[np.int32],
    indices: onp.ArrayND[np.int32] | None = None,
) -> None: ...
@overload  # distances=<given>, return_indices=True
def distance_transform_cdt(
    input: onp.ToFloatND,
    metric: _Metric2 | _Metric3 | onp.ToFloatND = "chessboard",
    return_distances: bool = True,
    *,
    return_indices: Literal[True],
    distances: onp.ArrayND[np.int32, _ShapeT],
    indices: None = None,
) -> onp.ArrayND[np.int32, _ShapeT]: ...
@overload  # return_distances=True (default), return_indices=False (default)
def distance_transform_cdt(
    input: onp.ToFloatND,
    metric: _Metric2 | _Metric3 | onp.ToFloatND = "chessboard",
    return_distances: Literal[True] = True,
    return_indices: Literal[False] = False,
    distances: None = None,
    indices: onp.ArrayND[np.int32, _ShapeT] | None = None,
) -> onp.ArrayND[np.int32, _ShapeT]: ...
@overload  # return_distances=True (default), return_indices=True
def distance_transform_cdt(
    input: onp.ToFloatND,
    metric: _Metric2 | _Metric3 | onp.ToFloatND = "chessboard",
    return_distances: Literal[True] = True,
    *,
    return_indices: Literal[True],
    distances: None = None,
    indices: None = None,
) -> tuple[onp.ArrayND[np.int32], onp.ArrayND[np.int32]]: ...

#
@overload  # return_distances=False, return_indices=False (default)
def distance_transform_edt(
    input: onp.ToFloatND,
    sampling: onp.ToFloat | onp.ToFloat1D | None = None,
    *,
    return_distances: Literal[False],
    return_indices: Literal[False] = False,
    distances: None = None,
    indices: onp.ArrayND[np.int32] | None = None,
) -> None: ...
@overload  # return_distances=False, return_indices=True
def distance_transform_edt(
    input: onp.ToFloatND,
    sampling: onp.ToFloat | onp.ToFloat1D | None = None,
    *,
    return_distances: Literal[False],
    return_indices: Literal[True],
    distances: None = None,
    indices: None = None,
) -> onp.ArrayND[np.int32]: ...
@overload  # distances=<given>, return_indices=False (default)
def distance_transform_edt(
    input: onp.ToFloatND,
    sampling: onp.ToFloat | onp.ToFloat1D | None = None,
    return_distances: bool = True,
    return_indices: Literal[False] = False,
    *,
    distances: onp.ArrayND[np.float64],
    indices: onp.ArrayND[np.int32] | None = None,
) -> None: ...
@overload  # distances=<given>, return_indices=True
def distance_transform_edt(
    input: onp.ToFloatND,
    sampling: onp.ToFloat | onp.ToFloat1D | None = None,
    return_distances: bool = True,
    *,
    return_indices: Literal[True],
    distances: onp.ArrayND[np.float64, _ShapeT],
    indices: None = None,
) -> onp.ArrayND[np.int32, _ShapeT]: ...
@overload  # return_distances=True (default), return_indices=False (default)
def distance_transform_edt(
    input: onp.ToFloatND,
    sampling: onp.ToFloat | onp.ToFloat1D | None = None,
    return_distances: Literal[True] = True,
    return_indices: Literal[False] = False,
    distances: None = None,
    indices: onp.ArrayND[np.int32, _ShapeT] | None = None,
) -> onp.ArrayND[np.float64, _ShapeT]: ...
@overload  # return_distances=True (default), return_indices=True
def distance_transform_edt(
    input: onp.ToFloatND,
    sampling: onp.ToFloat | onp.ToFloat1D | None = None,
    return_distances: Literal[True] = True,
    *,
    return_indices: Literal[True],
    distances: None = None,
    indices: None = None,
) -> tuple[onp.ArrayND[np.float64], onp.ArrayND[np.int32]]: ...
