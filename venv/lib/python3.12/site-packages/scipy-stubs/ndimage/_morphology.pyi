from typing import Literal, TypeAlias

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

_Mode: TypeAlias = Literal["reflect", "constant", "nearest", "mirror", "wrap"]
_MetricCDT: TypeAlias = Literal["chessboard", "taxicab"]
_MetricBF: TypeAlias = Literal["euclidean", _MetricCDT]
_BorderValue: TypeAlias = onp.ToInt | np.bool_

_BoolArrayOut: TypeAlias = onp.ArrayND[np.bool_]
_Origin: TypeAlias = int | tuple[int, ...]
_ScalarArrayOut: TypeAlias = onp.ArrayND[npc.number | np.bool_]

def iterate_structure(
    structure: onp.ToInt | onp.ToIntND, iterations: onp.ToInt, origin: _Origin | None = None
) -> _BoolArrayOut: ...
def generate_binary_structure(rank: int, connectivity: int) -> _BoolArrayOut: ...

#
def binary_erosion(
    input: onp.ToComplex | onp.ToComplexND,
    structure: onp.ToInt | onp.ToIntND | None = None,
    iterations: onp.ToInt = 1,
    mask: onp.ToInt | onp.ToIntND | None = None,
    output: _BoolArrayOut | type[bool | np.bool_] | None = None,
    border_value: _BorderValue = 0,
    origin: _Origin = 0,
    brute_force: onp.ToBool = False,
    *,
    axes: tuple[int, ...] | None = None,
) -> _BoolArrayOut: ...
def binary_dilation(
    input: onp.ToComplex | onp.ToComplexND,
    structure: onp.ToInt | onp.ToIntND | None = None,
    iterations: onp.ToInt = 1,
    mask: onp.ToInt | onp.ToIntND | None = None,
    output: _BoolArrayOut | type[bool | np.bool_] | None = None,
    border_value: _BorderValue = 0,
    origin: _Origin = 0,
    brute_force: onp.ToBool = False,
    *,
    axes: tuple[int, ...] | None = None,
) -> _BoolArrayOut: ...
def binary_opening(
    input: onp.ToComplex | onp.ToComplexND,
    structure: onp.ToInt | onp.ToIntND | None = None,
    iterations: onp.ToInt = 1,
    output: _BoolArrayOut | type[bool | np.bool_] | None = None,
    origin: _Origin = 0,
    mask: onp.ToInt | onp.ToIntND | None = None,
    border_value: _BorderValue = 0,
    brute_force: onp.ToBool = False,
    *,
    axes: tuple[int, ...] | None = None,
) -> _BoolArrayOut: ...
def binary_closing(
    input: onp.ToComplex | onp.ToComplexND,
    structure: onp.ToInt | onp.ToIntND | None = None,
    iterations: onp.ToInt = 1,
    output: _BoolArrayOut | type[bool | np.bool_] | None = None,
    origin: _Origin = 0,
    mask: onp.ToInt | onp.ToIntND | None = None,
    border_value: _BorderValue = 0,
    brute_force: onp.ToBool = False,
    *,
    axes: tuple[int, ...] | None = None,
) -> _BoolArrayOut: ...

#
def binary_hit_or_miss(
    input: onp.ToComplex | onp.ToComplexND,
    structure1: onp.ToInt | onp.ToIntND | None = None,
    structure2: onp.ToInt | onp.ToIntND | None = None,
    output: _BoolArrayOut | type[bool | np.bool_] | None = None,
    origin1: _Origin = 0,
    origin2: _Origin | None = None,
    *,
    axes: tuple[int, ...] | None = None,
) -> _BoolArrayOut: ...
def binary_propagation(
    input: onp.ToComplex | onp.ToComplexND,
    structure: onp.ToInt | onp.ToIntND | None = None,
    mask: onp.ToInt | onp.ToIntND | None = None,
    output: _BoolArrayOut | type[bool | np.bool_] | None = None,
    border_value: _BorderValue = 0,
    origin: _Origin = 0,
    *,
    axes: tuple[int, ...] | None = None,
) -> _BoolArrayOut: ...
def binary_fill_holes(
    input: onp.ToComplex | onp.ToComplexND,
    structure: onp.ToInt | onp.ToIntND | None = None,
    output: _BoolArrayOut | None = None,
    origin: _Origin = 0,
    *,
    axes: tuple[int, ...] | None = None,
) -> _BoolArrayOut: ...

#
def grey_erosion(
    input: onp.ToComplex | onp.ToComplexND,
    size: tuple[int, ...] | None = None,
    footprint: onp.ToScalar | onp.ToArrayND | None = None,
    structure: onp.ToInt | onp.ToIntND | None = None,
    output: _ScalarArrayOut | None = None,
    mode: _Mode = "reflect",
    cval: onp.ToComplex = 0.0,
    origin: _Origin = 0,
    *,
    axes: tuple[int, ...] | None = None,
) -> _ScalarArrayOut: ...
def grey_dilation(
    input: onp.ToComplex | onp.ToComplexND,
    size: tuple[int, ...] | None = None,
    footprint: onp.ToScalar | onp.ToArrayND | None = None,
    structure: onp.ToInt | onp.ToIntND | None = None,
    output: _ScalarArrayOut | None = None,
    mode: _Mode = "reflect",
    cval: onp.ToComplex = 0.0,
    origin: _Origin = 0,
    *,
    axes: tuple[int, ...] | None = None,
) -> _ScalarArrayOut: ...
def grey_opening(
    input: onp.ToComplex | onp.ToComplexND,
    size: tuple[int, ...] | None = None,
    footprint: onp.ToScalar | onp.ToArrayND | None = None,
    structure: onp.ToInt | onp.ToIntND | None = None,
    output: _ScalarArrayOut | None = None,
    mode: _Mode = "reflect",
    cval: onp.ToComplex = 0.0,
    origin: _Origin = 0,
    *,
    axes: tuple[int, ...] | None = None,
) -> _ScalarArrayOut: ...
def grey_closing(
    input: onp.ToComplex | onp.ToComplexND,
    size: tuple[int, ...] | None = None,
    footprint: onp.ToScalar | onp.ToArrayND | None = None,
    structure: onp.ToInt | onp.ToIntND | None = None,
    output: _ScalarArrayOut | None = None,
    mode: _Mode = "reflect",
    cval: onp.ToComplex = 0.0,
    origin: _Origin = 0,
    *,
    axes: tuple[int, ...] | None = None,
) -> _ScalarArrayOut: ...

#
def morphological_gradient(
    input: onp.ToComplex | onp.ToComplexND,
    size: tuple[int, ...] | None = None,
    footprint: onp.ToScalar | onp.ToArrayND | None = None,
    structure: onp.ToInt | onp.ToIntND | None = None,
    output: _ScalarArrayOut | None = None,
    mode: _Mode = "reflect",
    cval: onp.ToComplex = 0.0,
    origin: _Origin = 0,
    *,
    axes: tuple[int, ...] | None = None,
) -> _ScalarArrayOut: ...
def morphological_laplace(
    input: onp.ToComplex | onp.ToComplexND,
    size: tuple[int, ...] | None = None,
    footprint: onp.ToScalar | onp.ToArrayND | None = None,
    structure: onp.ToInt | onp.ToIntND | None = None,
    output: _ScalarArrayOut | None = None,
    mode: _Mode = "reflect",
    cval: onp.ToComplex = 0.0,
    origin: _Origin = 0,
    *,
    axes: tuple[int, ...] | None = None,
) -> _ScalarArrayOut: ...

#
def white_tophat(
    input: onp.ToComplex | onp.ToComplexND,
    size: tuple[int, ...] | None = None,
    footprint: onp.ToScalar | onp.ToArrayND | None = None,
    structure: onp.ToInt | onp.ToIntND | None = None,
    output: _ScalarArrayOut | None = None,
    mode: _Mode = "reflect",
    cval: onp.ToComplex = 0.0,
    origin: _Origin = 0,
    *,
    axes: tuple[int, ...] | None = None,
) -> _ScalarArrayOut: ...
def black_tophat(
    input: onp.ToComplex | onp.ToComplexND,
    size: tuple[int, ...] | None = None,
    footprint: onp.ToScalar | onp.ToArrayND | None = None,
    structure: onp.ToInt | onp.ToIntND | None = None,
    output: _ScalarArrayOut | None = None,
    mode: _Mode = "reflect",
    cval: onp.ToComplex = 0.0,
    origin: onp.ToComplex = 0,
    *,
    axes: tuple[int, ...] | None = None,
) -> _ScalarArrayOut: ...

#
def distance_transform_bf(
    input: onp.ToComplex | onp.ToComplexND,
    metric: _MetricBF = "euclidean",
    sampling: onp.ToFloat | onp.ToFloatND | None = None,
    return_distances: onp.ToBool = True,
    return_indices: onp.ToBool = False,
    distances: onp.ArrayND[np.float64 | np.uint32] | None = None,
    indices: onp.ArrayND[np.int32] | None = None,
) -> _ScalarArrayOut | onp.ArrayND[np.int32] | tuple[_ScalarArrayOut, onp.ArrayND[np.int32]]: ...
def distance_transform_cdt(
    input: onp.ToComplex | onp.ToComplexND,
    metric: _MetricCDT | onp.ToScalar | onp.ToArrayND = "chessboard",
    return_distances: onp.ToBool = True,
    return_indices: onp.ToBool = False,
    distances: onp.ArrayND[np.int32] | None = None,
    indices: onp.ArrayND[np.int32] | None = None,
) -> onp.ArrayND[np.int32] | tuple[onp.ArrayND[np.int32], onp.ArrayND[np.int32]]: ...
def distance_transform_edt(
    input: onp.ToComplex | onp.ToComplexND,
    sampling: onp.ToScalar | onp.ToArrayND | None = None,
    return_distances: onp.ToBool = True,
    return_indices: onp.ToBool = False,
    distances: onp.ArrayND[np.float64] | None = None,
    indices: onp.ArrayND[np.int32] | None = None,
) -> onp.ArrayND[np.float64] | onp.ArrayND[np.int32] | tuple[onp.ArrayND[np.float64], onp.ArrayND[np.int32]]: ...
