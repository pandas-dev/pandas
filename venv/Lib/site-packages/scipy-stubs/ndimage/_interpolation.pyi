from collections.abc import Callable
from typing import Any, Concatenate, Literal, overload

import numpy as np
import optype as op
import optype.numpy as onp
import optype.numpy.compat as npc

__all__ = [
    "affine_transform",
    "geometric_transform",
    "map_coordinates",
    "rotate",
    "shift",
    "spline_filter",
    "spline_filter1d",
    "zoom",
]

###

type _Order = Literal[0, 1, 2, 3, 4, 5]
type _Mode = Literal["reflect", "grid-mirror", "constant", "grid-constant", "nearest", "mirror", "wrap", "grid-wrap"]
type _MappingFunc = Callable[Concatenate[tuple[int, ...], ...], tuple[onp.ToFloat, ...]]
type _ArrayOrDType[ScalarT: np.generic] = onp.ArrayND[ScalarT] | type[ScalarT] | np.dtype[ScalarT]

#
@overload
def spline_filter1d(
    input: onp.ToFloatND,
    order: _Order = 3,
    axis: int = -1,
    output: type[float | np.float64] = np.float64,  # noqa: PYI011
    mode: _Mode = "mirror",
) -> onp.ArrayND[np.float64]: ...
@overload
def spline_filter1d(
    input: onp.ToComplexND, order: _Order = 3, axis: int = -1, *, output: type[op.JustComplex], mode: _Mode = "mirror"
) -> onp.ArrayND[np.complex128]: ...
@overload
def spline_filter1d[ScalarT: np.generic](
    input: onp.ToComplexND, order: _Order, axis: int, output: _ArrayOrDType[ScalarT], mode: _Mode = "mirror"
) -> onp.ArrayND[ScalarT]: ...
@overload
def spline_filter1d[ScalarT: np.generic](
    input: onp.ToComplexND, order: _Order = 3, axis: int = -1, *, output: _ArrayOrDType[ScalarT], mode: _Mode = "mirror"
) -> onp.ArrayND[ScalarT]: ...

#
@overload
def spline_filter(
    input: onp.ToFloatND,
    order: _Order = 3,
    output: type[float | np.float64] = np.float64,  # noqa: PYI011
    mode: _Mode = "mirror",
) -> onp.ArrayND[np.float64]: ...
@overload
def spline_filter(
    input: onp.ToComplexND, order: _Order = 3, *, output: type[op.JustComplex], mode: _Mode = "mirror"
) -> onp.ArrayND[np.complex128]: ...
@overload
def spline_filter[ScalarT: np.generic](
    input: onp.ToComplexND, order: _Order, output: _ArrayOrDType[ScalarT], mode: _Mode = "mirror"
) -> onp.ArrayND[ScalarT]: ...
@overload
def spline_filter[ScalarT: np.generic](
    input: onp.ToComplexND, order: _Order = 3, *, output: _ArrayOrDType[ScalarT], mode: _Mode = "mirror"
) -> onp.ArrayND[ScalarT]: ...

#
@overload
def geometric_transform[ArrayT: onp.ArrayND[np.bool | npc.integer | npc.inexact32 | npc.inexact64]](
    input: ArrayT,
    mapping: _MappingFunc,
    output_shape: tuple[int, ...] | None = None,
    output: None = None,
    order: _Order = 3,
    mode: _Mode = "constant",
    cval: onp.ToFloat = 0.0,
    prefilter: bool = True,
    extra_arguments: tuple[object, ...] = (),
    extra_keywords: dict[str, Any] | None = None,
) -> ArrayT: ...
@overload
def geometric_transform(
    input: onp.SequenceND[int],
    mapping: _MappingFunc,
    output_shape: tuple[int, ...] | None = None,
    output: None = None,
    order: _Order = 3,
    mode: _Mode = "constant",
    cval: onp.ToFloat = 0.0,
    prefilter: bool = True,
    extra_arguments: tuple[object, ...] = (),
    extra_keywords: dict[str, Any] | None = None,
) -> onp.ArrayND[np.int_]: ...
@overload
def geometric_transform(
    input: onp.SequenceND[list[float]] | list[float],
    mapping: _MappingFunc,
    output_shape: tuple[int, ...] | None = None,
    output: None = None,
    order: _Order = 3,
    mode: _Mode = "constant",
    cval: onp.ToFloat = 0.0,
    prefilter: bool = True,
    extra_arguments: tuple[object, ...] = (),
    extra_keywords: dict[str, Any] | None = None,
) -> onp.ArrayND[np.float64]: ...
@overload
def geometric_transform(
    input: onp.SequenceND[list[complex]] | list[complex],
    mapping: _MappingFunc,
    output_shape: tuple[int, ...] | None = None,
    output: None = None,
    order: _Order = 3,
    mode: _Mode = "constant",
    cval: onp.ToComplex = 0.0,
    prefilter: bool = True,
    extra_arguments: tuple[object, ...] = (),
    extra_keywords: dict[str, Any] | None = None,
) -> onp.ArrayND[np.complex128]: ...
@overload
def geometric_transform(
    input: onp.ToComplexND,
    mapping: _MappingFunc,
    output_shape: tuple[int, ...] | None = None,
    output: None = None,
    order: _Order = 3,
    mode: _Mode = "constant",
    cval: onp.ToComplex = 0.0,
    prefilter: bool = True,
    extra_arguments: tuple[object, ...] = (),
    extra_keywords: dict[str, Any] | None = None,
) -> onp.ArrayND[Any]: ...
@overload
def geometric_transform[ScalarT: np.generic](
    input: onp.ToComplexND,
    mapping: _MappingFunc,
    output_shape: tuple[int, ...] | None,
    output: _ArrayOrDType[ScalarT],
    order: _Order = 3,
    mode: _Mode = "constant",
    cval: onp.ToComplex = 0.0,
    prefilter: bool = True,
    extra_arguments: tuple[object, ...] = (),
    extra_keywords: dict[str, Any] | None = None,
) -> onp.ArrayND[ScalarT]: ...
@overload
def geometric_transform[ScalarT: np.generic](
    input: onp.ToComplexND,
    mapping: _MappingFunc,
    output_shape: tuple[int, ...] | None = None,
    *,
    output: _ArrayOrDType[ScalarT],
    order: _Order = 3,
    mode: _Mode = "constant",
    cval: onp.ToComplex = 0.0,
    prefilter: bool = True,
    extra_arguments: tuple[object, ...] = (),
    extra_keywords: dict[str, Any] | None = None,
) -> onp.ArrayND[ScalarT]: ...
@overload
def geometric_transform(
    input: onp.ToComplexND,
    mapping: _MappingFunc,
    output_shape: tuple[int, ...] | None = None,
    *,
    output: type[int],
    order: _Order = 3,
    mode: _Mode = "constant",
    cval: onp.ToComplex = 0.0,
    prefilter: bool = True,
    extra_arguments: tuple[object, ...] = (),
    extra_keywords: dict[str, Any] | None = None,
) -> onp.ArrayND[np.int_]: ...
@overload
def geometric_transform(
    input: onp.ToComplexND,
    mapping: _MappingFunc,
    output_shape: tuple[int, ...] | None = None,
    *,
    output: type[op.JustFloat],
    order: _Order = 3,
    mode: _Mode = "constant",
    cval: onp.ToComplex = 0.0,
    prefilter: bool = True,
    extra_arguments: tuple[object, ...] = (),
    extra_keywords: dict[str, Any] | None = None,
) -> onp.ArrayND[np.float64]: ...
@overload
def geometric_transform(
    input: onp.ToComplexND,
    mapping: _MappingFunc,
    output_shape: tuple[int, ...] | None = None,
    *,
    output: type[op.JustComplex],
    order: _Order = 3,
    mode: _Mode = "constant",
    cval: onp.ToComplex = 0.0,
    prefilter: bool = True,
    extra_arguments: tuple[object, ...] = (),
    extra_keywords: dict[str, Any] | None = None,
) -> onp.ArrayND[np.complex128]: ...

#
@overload
def map_coordinates[NumberT: npc.number](
    input: onp.ArrayND[NumberT],
    coordinates: onp.ToFloatND,
    output: None = None,
    order: _Order = 3,
    mode: _Mode = "constant",
    cval: onp.ToFloat = 0.0,
    prefilter: bool = True,
) -> onp.ArrayND[NumberT]: ...
@overload
def map_coordinates(
    input: onp.SequenceND[int],
    coordinates: onp.ToFloatND,
    output: None = None,
    order: _Order = 3,
    mode: _Mode = "constant",
    cval: onp.ToFloat = 0.0,
    prefilter: bool = True,
) -> onp.ArrayND[np.int_]: ...
@overload
def map_coordinates(
    input: onp.SequenceND[list[float]] | list[float],
    coordinates: onp.ToFloatND,
    output: None = None,
    order: _Order = 3,
    mode: _Mode = "constant",
    cval: onp.ToFloat = 0.0,
    prefilter: bool = True,
) -> onp.ArrayND[np.float64]: ...
@overload
def map_coordinates(
    input: onp.SequenceND[list[complex]] | list[complex],
    coordinates: onp.ToFloatND,
    output: None = None,
    order: _Order = 3,
    mode: _Mode = "constant",
    cval: onp.ToFloat = 0.0,
    prefilter: bool = True,
) -> onp.ArrayND[np.complex128]: ...
@overload
def map_coordinates(
    input: onp.ToComplexND,
    coordinates: onp.ToFloatND,
    output: None = None,
    order: _Order = 3,
    mode: _Mode = "constant",
    cval: onp.ToComplex = 0.0,
    prefilter: bool = True,
) -> onp.ArrayND[Any]: ...
@overload
def map_coordinates[ScalarT: np.generic](
    input: onp.ToComplexND,
    coordinates: onp.ToFloatND,
    output: _ArrayOrDType[ScalarT],
    order: _Order = 3,
    mode: _Mode = "constant",
    cval: onp.ToComplex = 0.0,
    prefilter: bool = True,
) -> onp.ArrayND[ScalarT]: ...
@overload
def map_coordinates(
    input: onp.ToComplexND,
    coordinates: onp.ToFloatND,
    output: type[bool],
    order: _Order = 3,
    mode: _Mode = "constant",
    cval: onp.ToFloat = 0.0,
    prefilter: bool = True,
) -> onp.ArrayND[np.bool]: ...
@overload
def map_coordinates(
    input: onp.ToComplexND,
    coordinates: onp.ToFloatND,
    output: type[op.JustInt],
    order: _Order = 3,
    mode: _Mode = "constant",
    cval: onp.ToFloat = 0.0,
    prefilter: bool = True,
) -> onp.ArrayND[np.int_]: ...
@overload
def map_coordinates(
    input: onp.ToComplexND,
    coordinates: onp.ToFloatND,
    output: type[op.JustFloat],
    order: _Order = 3,
    mode: _Mode = "constant",
    cval: onp.ToFloat = 0.0,
    prefilter: bool = True,
) -> onp.ArrayND[np.float64]: ...
@overload
def map_coordinates(
    input: onp.ToComplexND,
    coordinates: onp.ToFloatND,
    output: type[op.JustComplex],
    order: _Order = 3,
    mode: _Mode = "constant",
    cval: onp.ToComplex = 0.0,
    prefilter: bool = True,
) -> onp.ArrayND[np.complex128]: ...

#
@overload
def affine_transform[NumberT: npc.number](
    input: onp.ArrayND[NumberT],
    matrix: onp.ToFloat1D | onp.ToFloat2D,
    offset: onp.ToFloat | onp.ToFloat1D = 0.0,
    output_shape: tuple[int, ...] | None = None,
    output: None = None,
    order: _Order = 3,
    mode: _Mode = "constant",
    cval: onp.ToFloat = 0.0,
    prefilter: bool = True,
) -> onp.ArrayND[NumberT]: ...
@overload
def affine_transform(
    input: onp.SequenceND[int],
    matrix: onp.ToFloat1D | onp.ToFloat2D,
    offset: onp.ToFloat | onp.ToFloat1D = 0.0,
    output_shape: tuple[int, ...] | None = None,
    output: None = None,
    order: _Order = 3,
    mode: _Mode = "constant",
    cval: onp.ToFloat = 0.0,
    prefilter: bool = True,
) -> onp.ArrayND[np.int_]: ...
@overload
def affine_transform(
    input: onp.SequenceND[list[float]] | list[float],
    matrix: onp.ToFloat1D | onp.ToFloat2D,
    offset: onp.ToFloat | onp.ToFloat1D = 0.0,
    output_shape: tuple[int, ...] | None = None,
    output: None = None,
    order: _Order = 3,
    mode: _Mode = "constant",
    cval: onp.ToFloat = 0.0,
    prefilter: bool = True,
) -> onp.ArrayND[np.float64]: ...
@overload
def affine_transform(
    input: onp.SequenceND[list[complex]] | list[complex],
    matrix: onp.ToFloat1D | onp.ToFloat2D,
    offset: onp.ToFloat | onp.ToFloat1D = 0.0,
    output_shape: tuple[int, ...] | None = None,
    output: None = None,
    order: _Order = 3,
    mode: _Mode = "constant",
    cval: onp.ToFloat = 0.0,
    prefilter: bool = True,
) -> onp.ArrayND[np.complex128]: ...
@overload
def affine_transform(
    input: onp.ToComplexND,
    matrix: onp.ToFloat1D | onp.ToFloat2D,
    offset: onp.ToFloat | onp.ToFloat1D = 0.0,
    output_shape: tuple[int, ...] | None = None,
    output: None = None,
    order: _Order = 3,
    mode: _Mode = "constant",
    cval: onp.ToComplex = 0.0,
    prefilter: bool = True,
) -> onp.ArrayND[Any]: ...
@overload
def affine_transform[ScalarT: np.generic](
    input: onp.ToComplexND,
    matrix: onp.ToFloat1D | onp.ToFloat2D,
    offset: onp.ToFloat | onp.ToFloat1D = 0.0,
    output_shape: tuple[int, ...] | None = None,
    *,
    output: _ArrayOrDType[ScalarT],
    order: _Order = 3,
    mode: _Mode = "constant",
    cval: onp.ToComplex = 0.0,
    prefilter: bool = True,
) -> onp.ArrayND[ScalarT]: ...
@overload
def affine_transform(
    input: onp.ToComplexND,
    matrix: onp.ToFloat1D | onp.ToFloat2D,
    offset: onp.ToFloat | onp.ToFloat1D = 0.0,
    output_shape: tuple[int, ...] | None = None,
    *,
    output: type[int],
    order: _Order = 3,
    mode: _Mode = "constant",
    cval: onp.ToFloat = 0.0,
    prefilter: bool = True,
) -> onp.ArrayND[np.int_]: ...
@overload
def affine_transform(
    input: onp.ToComplexND,
    matrix: onp.ToFloat1D | onp.ToFloat2D,
    offset: onp.ToFloat | onp.ToFloat1D = 0.0,
    output_shape: tuple[int, ...] | None = None,
    *,
    output: type[op.JustFloat],
    order: _Order = 3,
    mode: _Mode = "constant",
    cval: onp.ToFloat = 0.0,
    prefilter: bool = True,
) -> onp.ArrayND[np.float64]: ...
@overload
def affine_transform(
    input: onp.ToComplexND,
    matrix: onp.ToFloat1D | onp.ToFloat2D,
    offset: onp.ToFloat | onp.ToFloat1D = 0.0,
    output_shape: tuple[int, ...] | None = None,
    *,
    output: type[op.JustComplex],
    order: _Order = 3,
    mode: _Mode = "constant",
    cval: onp.ToComplex = 0.0,
    prefilter: bool = True,
) -> onp.ArrayND[np.complex128]: ...

#
@overload
def shift[ArrayT: onp.ArrayND[np.bool | npc.integer | npc.inexact32 | npc.inexact64]](
    input: ArrayT,
    shift: onp.ToFloat | onp.ToFloatND,
    output: None = None,
    order: _Order = 3,
    mode: _Mode = "constant",
    cval: onp.ToFloat = 0.0,
    prefilter: bool = True,
) -> ArrayT: ...
@overload
def shift(
    input: onp.SequenceND[int],
    shift: onp.ToFloat | onp.ToFloatND,
    output: None = None,
    order: _Order = 3,
    mode: _Mode = "constant",
    cval: onp.ToFloat = 0.0,
    prefilter: bool = True,
) -> onp.ArrayND[np.int_]: ...
@overload
def shift(
    input: onp.SequenceND[list[float]] | list[float],
    shift: onp.ToFloat | onp.ToFloatND,
    output: None = None,
    order: _Order = 3,
    mode: _Mode = "constant",
    cval: onp.ToFloat = 0.0,
    prefilter: bool = True,
) -> onp.ArrayND[np.float64]: ...
@overload
def shift(
    input: onp.SequenceND[list[complex]] | list[complex],
    shift: onp.ToFloat | onp.ToFloatND,
    output: None = None,
    order: _Order = 3,
    mode: _Mode = "constant",
    cval: onp.ToFloat = 0.0,
    prefilter: bool = True,
) -> onp.ArrayND[np.complex128]: ...
@overload
def shift(
    input: onp.ToComplexND,
    shift: onp.ToFloat | onp.ToFloatND,
    output: None = None,
    order: _Order = 3,
    mode: _Mode = "constant",
    cval: onp.ToComplex = 0.0,
    prefilter: bool = True,
) -> onp.ArrayND[Any]: ...
@overload
def shift[ScalarT: np.generic](
    input: onp.ToComplexND,
    shift: onp.ToFloat | onp.ToFloatND,
    output: _ArrayOrDType[ScalarT],
    order: _Order = 3,
    mode: _Mode = "constant",
    cval: onp.ToComplex = 0.0,
    prefilter: bool = True,
) -> onp.ArrayND[ScalarT]: ...
@overload
def shift(
    input: onp.ToComplexND,
    shift: onp.ToFloat | onp.ToFloatND,
    output: type[int],
    order: _Order = 3,
    mode: _Mode = "constant",
    cval: onp.ToFloat = 0.0,
    prefilter: bool = True,
) -> onp.ArrayND[np.int_]: ...
@overload
def shift(
    input: onp.ToComplexND,
    shift: onp.ToFloat | onp.ToFloatND,
    output: type[op.JustFloat],
    order: _Order = 3,
    mode: _Mode = "constant",
    cval: onp.ToFloat = 0.0,
    prefilter: bool = True,
) -> onp.ArrayND[np.float64]: ...
@overload
def shift(
    input: onp.ToComplexND,
    shift: onp.ToFloat | onp.ToFloatND,
    output: type[op.JustComplex],
    order: _Order = 3,
    mode: _Mode = "constant",
    cval: onp.ToComplex = 0.0,
    prefilter: bool = True,
) -> onp.ArrayND[np.complex128]: ...

#
@overload
def zoom[ArrayT: onp.ArrayND[np.bool | npc.integer | npc.inexact32 | npc.inexact64]](
    input: ArrayT,
    zoom: onp.ToFloat | onp.ToFloatND,
    output: None = None,
    order: _Order = 3,
    mode: _Mode = "constant",
    cval: onp.ToFloat = 0.0,
    prefilter: bool = True,
    *,
    grid_mode: bool = False,
) -> ArrayT: ...
@overload
def zoom(
    input: onp.SequenceND[int],
    zoom: onp.ToFloat | onp.ToFloatND,
    output: None = None,
    order: _Order = 3,
    mode: _Mode = "constant",
    cval: onp.ToFloat = 0.0,
    prefilter: bool = True,
    *,
    grid_mode: bool = False,
) -> onp.ArrayND[np.int_]: ...
@overload
def zoom(
    input: onp.SequenceND[list[float]] | list[float],
    zoom: onp.ToFloat | onp.ToFloatND,
    output: None = None,
    order: _Order = 3,
    mode: _Mode = "constant",
    cval: onp.ToFloat = 0.0,
    prefilter: bool = True,
    *,
    grid_mode: bool = False,
) -> onp.ArrayND[np.float64]: ...
@overload
def zoom(
    input: onp.SequenceND[list[complex]] | list[complex],
    zoom: onp.ToFloat | onp.ToFloatND,
    output: None = None,
    order: _Order = 3,
    mode: _Mode = "constant",
    cval: onp.ToFloat = 0.0,
    prefilter: bool = True,
    *,
    grid_mode: bool = False,
) -> onp.ArrayND[np.complex128]: ...
@overload
def zoom(
    input: onp.ToComplexND,
    zoom: onp.ToFloat | onp.ToFloatND,
    output: None = None,
    order: _Order = 3,
    mode: _Mode = "constant",
    cval: onp.ToComplex = 0.0,
    prefilter: bool = True,
    *,
    grid_mode: bool = False,
) -> onp.ArrayND[Any]: ...
@overload
def zoom[ScalarT: np.generic](
    input: onp.ToComplexND,
    zoom: onp.ToFloat | onp.ToFloatND,
    output: _ArrayOrDType[ScalarT],
    order: _Order = 3,
    mode: _Mode = "constant",
    cval: onp.ToComplex = 0.0,
    prefilter: bool = True,
    *,
    grid_mode: bool = False,
) -> onp.ArrayND[ScalarT]: ...
@overload
def zoom(
    input: onp.ToComplexND,
    zoom: onp.ToFloat | onp.ToFloatND,
    output: type[int],
    order: _Order = 3,
    mode: _Mode = "constant",
    cval: onp.ToFloat = 0.0,
    prefilter: bool = True,
    *,
    grid_mode: bool = False,
) -> onp.ArrayND[np.int_]: ...
@overload
def zoom(
    input: onp.ToComplexND,
    zoom: onp.ToFloat | onp.ToFloatND,
    output: type[op.JustFloat],
    order: _Order = 3,
    mode: _Mode = "constant",
    cval: onp.ToFloat = 0.0,
    prefilter: bool = True,
    *,
    grid_mode: bool = False,
) -> onp.ArrayND[np.float64]: ...
@overload
def zoom(
    input: onp.ToComplexND,
    zoom: onp.ToFloat | onp.ToFloatND,
    output: type[op.JustComplex],
    order: _Order = 3,
    mode: _Mode = "constant",
    cval: onp.ToComplex = 0.0,
    prefilter: bool = True,
    *,
    grid_mode: bool = False,
) -> onp.ArrayND[np.complex128]: ...

#
@overload
def rotate[ArrayT: onp.ArrayND[np.bool | npc.integer | npc.inexact32 | npc.inexact64]](
    input: ArrayT,
    angle: onp.ToFloat,
    axes: tuple[int, int] = (1, 0),
    reshape: bool = True,
    output: None = None,
    order: _Order = 3,
    mode: _Mode = "constant",
    cval: onp.ToFloat = 0.0,
    prefilter: bool = True,
) -> ArrayT: ...
@overload
def rotate(
    input: onp.SequenceND[int],
    angle: onp.ToFloat,
    axes: tuple[int, int] = (1, 0),
    reshape: bool = True,
    output: None = None,
    order: _Order = 3,
    mode: _Mode = "constant",
    cval: onp.ToFloat = 0.0,
    prefilter: bool = True,
) -> onp.ArrayND[np.int_]: ...
@overload
def rotate(
    input: onp.SequenceND[list[float]] | list[float],
    angle: onp.ToFloat,
    axes: tuple[int, int] = (1, 0),
    reshape: bool = True,
    output: None = None,
    order: _Order = 3,
    mode: _Mode = "constant",
    cval: onp.ToFloat = 0.0,
    prefilter: bool = True,
) -> onp.ArrayND[np.float64]: ...
@overload
def rotate(
    input: onp.SequenceND[list[complex]] | list[complex],
    angle: onp.ToFloat,
    axes: tuple[int, int] = (1, 0),
    reshape: bool = True,
    output: None = None,
    order: _Order = 3,
    mode: _Mode = "constant",
    cval: onp.ToFloat = 0.0,
    prefilter: bool = True,
) -> onp.ArrayND[np.complex128]: ...
@overload
def rotate(
    input: onp.ToComplexND,
    angle: onp.ToFloat,
    axes: tuple[int, int] = (1, 0),
    reshape: bool = True,
    output: None = None,
    order: _Order = 3,
    mode: _Mode = "constant",
    cval: onp.ToComplex = 0.0,
    prefilter: bool = True,
) -> onp.ArrayND[Any]: ...
@overload
def rotate[ScalarT: np.generic](
    input: onp.ToComplexND,
    angle: onp.ToFloat,
    axes: tuple[int, int] = (1, 0),
    reshape: bool = True,
    *,
    output: _ArrayOrDType[ScalarT],
    order: _Order = 3,
    mode: _Mode = "constant",
    cval: onp.ToComplex = 0.0,
    prefilter: bool = True,
) -> onp.ArrayND[ScalarT]: ...
@overload
def rotate(
    input: onp.ToComplexND,
    angle: onp.ToFloat,
    axes: tuple[int, int] = (1, 0),
    reshape: bool = True,
    *,
    output: type[int],
    order: _Order = 3,
    mode: _Mode = "constant",
    cval: onp.ToComplex = 0.0,
    prefilter: bool = True,
) -> onp.ArrayND[np.int_]: ...
@overload
def rotate(
    input: onp.ToComplexND,
    angle: onp.ToFloat,
    axes: tuple[int, int] = (1, 0),
    reshape: bool = True,
    *,
    output: type[op.JustFloat],
    order: _Order = 3,
    mode: _Mode = "constant",
    cval: onp.ToComplex = 0.0,
    prefilter: bool = True,
) -> onp.ArrayND[np.float64]: ...
@overload
def rotate(
    input: onp.ToComplexND,
    angle: onp.ToFloat,
    axes: tuple[int, int] = (1, 0),
    reshape: bool = True,
    *,
    output: type[op.JustComplex],
    order: _Order = 3,
    mode: _Mode = "constant",
    cval: onp.ToComplex = 0.0,
    prefilter: bool = True,
) -> onp.ArrayND[np.complex128]: ...
