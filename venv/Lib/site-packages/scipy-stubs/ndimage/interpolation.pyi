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
    offset: object = 0.0,
    output_shape: object = None,
    output: object = None,
    order: object = 3,
    mode: object = "constant",
    cval: object = 0.0,
    prefilter: object = True,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def geometric_transform(
    input: object,
    mapping: object,
    output_shape: object = None,
    output: object = None,
    order: object = 3,
    mode: object = "constant",
    cval: object = 0.0,
    prefilter: object = True,
    extra_arguments: object = (),
    extra_keywords: object = None,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def map_coordinates(
    input: object,
    coordinates: object,
    output: object = None,
    order: object = 3,
    mode: object = "constant",
    cval: object = 0.0,
    prefilter: object = True,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def rotate(
    input: object,
    angle: object,
    axes: object = (1, 0),
    reshape: object = True,
    output: object = None,
    order: object = 3,
    mode: object = "constant",
    cval: object = 0.0,
    prefilter: object = True,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def shift(
    input: object,
    shift: object,
    output: object = None,
    order: object = 3,
    mode: object = "constant",
    cval: object = 0.0,
    prefilter: object = True,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def spline_filter(input: object, order: object = 3, output: object = ..., mode: object = "mirror") -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def spline_filter1d(
    input: object, order: object = 3, axis: object = -1, output: object = ..., mode: object = "mirror"
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def zoom(
    input: object,
    zoom: object,
    output: object = None,
    order: object = 3,
    mode: object = "constant",
    cval: object = 0.0,
    prefilter: object = True,
    *,
    grid_mode: object = False,
) -> object: ...
