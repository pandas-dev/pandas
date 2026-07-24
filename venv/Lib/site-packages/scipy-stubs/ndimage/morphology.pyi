# This module is not meant for public use and will be removed in SciPy v2.0.0.
from typing_extensions import deprecated

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

#
@deprecated("will be removed in SciPy v2.0.0")
def binary_closing(
    input: object,
    structure: object = None,
    iterations: object = 1,
    output: object = None,
    origin: object = 0,
    mask: object = None,
    border_value: object = 0,
    brute_force: object = False,
    *,
    axes: tuple[int, ...] | None = None,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def binary_dilation(
    input: object,
    structure: object = None,
    iterations: object = 1,
    mask: object = None,
    output: object = None,
    border_value: object = 0,
    origin: object = 0,
    brute_force: object = False,
    *,
    axes: tuple[int, ...] | None = None,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def binary_erosion(
    input: object,
    structure: object = None,
    iterations: object = 1,
    mask: object = None,
    output: object = None,
    border_value: object = 0,
    origin: object = 0,
    brute_force: object = False,
    *,
    axes: tuple[int, ...] | None = None,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def binary_opening(
    input: object,
    structure: object = None,
    iterations: object = 1,
    output: object = None,
    origin: object = 0,
    mask: object = None,
    border_value: object = 0,
    brute_force: object = False,
    *,
    axes: tuple[int, ...] | None = None,
) -> object: ...

#
@deprecated("will be removed in SciPy v2.0.0")
def binary_fill_holes(
    input: object, structure: object = None, output: object = None, origin: object = 0, *, axes: tuple[int, ...] | None = None
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def binary_hit_or_miss(
    input: object,
    structure1: object = None,
    structure2: object = None,
    output: object = None,
    origin1: object = 0,
    origin2: object = None,
    *,
    axes: tuple[int, ...] | None = None,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def binary_propagation(
    input: object,
    structure: object = None,
    mask: object = None,
    output: object = None,
    border_value: object = 0,
    origin: object = 0,
    *,
    axes: tuple[int, ...] | None = None,
) -> object: ...

#
@deprecated("will be removed in SciPy v2.0.0")
def distance_transform_bf(
    input: object,
    metric: object = "euclidean",
    sampling: object = None,
    return_distances: object = True,
    return_indices: object = False,
    distances: object = None,
    indices: object = None,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def distance_transform_cdt(
    input: object,
    metric: object = "chessboard",
    return_distances: object = True,
    return_indices: object = False,
    distances: object = None,
    indices: object = None,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def distance_transform_edt(
    input: object,
    sampling: object = None,
    return_distances: object = True,
    return_indices: object = False,
    distances: object = None,
    indices: object = None,
) -> object: ...

#
@deprecated("will be removed in SciPy v2.0.0")
def grey_closing(
    input: object,
    size: object = None,
    footprint: object = None,
    structure: object = None,
    output: object = None,
    mode: object = "reflect",
    cval: object = 0.0,
    origin: object = 0,
    *,
    axes: tuple[int, ...] | None = None,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def grey_dilation(
    input: object,
    size: object = None,
    footprint: object = None,
    structure: object = None,
    output: object = None,
    mode: object = "reflect",
    cval: object = 0.0,
    origin: object = 0,
    *,
    axes: tuple[int, ...] | None = None,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def grey_erosion(
    input: object,
    size: object = None,
    footprint: object = None,
    structure: object = None,
    output: object = None,
    mode: object = "reflect",
    cval: object = 0.0,
    origin: object = 0,
    *,
    axes: tuple[int, ...] | None = None,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def grey_opening(
    input: object,
    size: object = None,
    footprint: object = None,
    structure: object = None,
    output: object = None,
    mode: object = "reflect",
    cval: object = 0.0,
    origin: object = 0,
    *,
    axes: tuple[int, ...] | None = None,
) -> object: ...

#
@deprecated("will be removed in SciPy v2.0.0")
def morphological_gradient(
    input: object,
    size: object = None,
    footprint: object = None,
    structure: object = None,
    output: object = None,
    mode: object = "reflect",
    cval: object = 0.0,
    origin: object = 0,
    *,
    axes: tuple[int, ...] | None = None,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def morphological_laplace(
    input: object,
    size: object = None,
    footprint: object = None,
    structure: object = None,
    output: object = None,
    mode: object = "reflect",
    cval: object = 0.0,
    origin: object = 0,
    *,
    axes: tuple[int, ...] | None = None,
) -> object: ...

#
@deprecated("will be removed in SciPy v2.0.0")
def white_tophat(
    input: object,
    size: object = None,
    footprint: object = None,
    structure: object = None,
    output: object = None,
    mode: object = "reflect",
    cval: object = 0.0,
    origin: object = 0,
    *,
    axes: tuple[int, ...] | None = None,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def black_tophat(
    input: object,
    size: object = None,
    footprint: object = None,
    structure: object = None,
    output: object = None,
    mode: object = "reflect",
    cval: object = 0.0,
    origin: object = 0,
    *,
    axes: tuple[int, ...] | None = None,
) -> object: ...

#
@deprecated("will be removed in SciPy v2.0.0")
def iterate_structure(structure: object, iterations: object, origin: object = None) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def generate_binary_structure(rank: object, connectivity: object) -> object: ...
