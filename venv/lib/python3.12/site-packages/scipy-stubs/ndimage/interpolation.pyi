# This module is not meant for public use and will be removed in SciPy v2.0.0.
from typing_extensions import deprecated

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

@deprecated("will be removed in SciPy v2.0.0")
def affine_transform(
    input: object,
    matrix: object,
    offset: object = ...,
    output_shape: object = ...,
    output: object = ...,
    order: object = ...,
    mode: object = ...,
    cval: object = ...,
    prefilter: object = ...,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def geometric_transform(
    input: object,
    mapping: object,
    output_shape: object = ...,
    output: object = ...,
    order: object = ...,
    mode: object = ...,
    cval: object = ...,
    prefilter: object = ...,
    extra_arguments: object = ...,
    extra_keywords: object = ...,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def map_coordinates(
    input: object,
    coordinates: object,
    output: object = ...,
    order: object = ...,
    mode: object = ...,
    cval: object = ...,
    prefilter: object = ...,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def rotate(
    input: object,
    angle: object,
    axes: object = ...,
    reshape: object = ...,
    output: object = ...,
    order: object = ...,
    mode: object = ...,
    cval: object = ...,
    prefilter: object = ...,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def shift(
    input: object,
    shift: object,
    output: object = ...,
    order: object = ...,
    mode: object = ...,
    cval: object = ...,
    prefilter: object = ...,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def spline_filter(input: object, order: object = ..., output: object = ..., mode: object = ...) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def spline_filter1d(
    input: object, order: object = ..., axis: object = ..., output: object = ..., mode: object = ...
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def zoom(
    input: object,
    zoom: object,
    output: object = ...,
    order: object = ...,
    mode: object = ...,
    cval: object = ...,
    prefilter: object = ...,
    *,
    grid_mode: object = ...,
) -> object: ...
