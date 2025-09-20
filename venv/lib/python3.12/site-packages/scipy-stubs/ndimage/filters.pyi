# This module is not meant for public use and will be removed in SciPy v2.0.0.
from typing_extensions import deprecated

__all__ = [
    "convolve",
    "convolve1d",
    "correlate",
    "correlate1d",
    "gaussian_filter",
    "gaussian_filter1d",
    "gaussian_gradient_magnitude",
    "gaussian_laplace",
    "generic_filter",
    "generic_filter1d",
    "generic_gradient_magnitude",
    "generic_laplace",
    "laplace",
    "maximum_filter",
    "maximum_filter1d",
    "median_filter",
    "minimum_filter",
    "minimum_filter1d",
    "percentile_filter",
    "prewitt",
    "rank_filter",
    "sobel",
    "uniform_filter",
    "uniform_filter1d",
]

@deprecated("will be removed in SciPy v2.0.0")
def convolve(
    input: object,
    weights: object,
    output: object = ...,
    mode: object = ...,
    cval: object = ...,
    origin: object = ...,
    *,
    axes: tuple[int, ...] | None = None,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def convolve1d(
    input: object,
    weights: object,
    axis: object = ...,
    output: object = ...,
    mode: object = ...,
    cval: object = ...,
    origin: object = ...,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def correlate(
    input: object,
    weights: object,
    output: object = ...,
    mode: object = ...,
    cval: object = ...,
    origin: object = ...,
    *,
    axes: tuple[int, ...] | None = None,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def correlate1d(
    input: object,
    weights: object,
    axis: object = ...,
    output: object = ...,
    mode: object = ...,
    cval: object = ...,
    origin: object = ...,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def gaussian_filter(
    input: object,
    sigma: object,
    order: object = ...,
    output: object = ...,
    mode: object = ...,
    cval: object = ...,
    truncate: object = ...,
    *,
    radius: object = ...,
    axes: object = ...,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def gaussian_filter1d(
    input: object,
    sigma: object,
    axis: object = ...,
    order: object = ...,
    output: object = ...,
    mode: object = ...,
    cval: object = ...,
    truncate: object = ...,
    *,
    radius: object = ...,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def gaussian_gradient_magnitude(
    input: object,
    sigma: object,
    output: object = ...,
    mode: object = ...,
    cval: object = ...,
    *,
    axes: tuple[int, ...] | None = None,
    **kwargs: object,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def gaussian_laplace(
    input: object,
    sigma: object,
    output: object = ...,
    mode: object = ...,
    cval: object = ...,
    *,
    axes: tuple[int, ...] | None = None,
    **kwargs: object,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def generic_filter(
    input: object,
    function: object,
    size: object = ...,
    footprint: object = ...,
    output: object = ...,
    mode: object = ...,
    cval: object = ...,
    origin: object = ...,
    extra_arguments: object = ...,
    extra_keywords: object = ...,
    *,
    axes: tuple[int, ...] | None = None,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def generic_filter1d(
    input: object,
    function: object,
    filter_size: object,
    axis: object = ...,
    output: object = ...,
    mode: object = ...,
    cval: object = ...,
    origin: object = ...,
    extra_arguments: object = ...,
    extra_keywords: object = ...,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def generic_gradient_magnitude(
    input: object,
    derivative: object,
    output: object = ...,
    mode: object = ...,
    cval: object = ...,
    extra_arguments: object = ...,
    extra_keywords: object = ...,
    *,
    axes: tuple[int, ...] | None = None,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def generic_laplace(
    input: object,
    derivative2: object,
    output: object = ...,
    mode: object = ...,
    cval: object = ...,
    extra_arguments: object = ...,
    extra_keywords: object = ...,
    *,
    axes: tuple[int, ...] | None = None,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def laplace(
    input: object, output: object = ..., mode: object = ..., cval: object = ..., *, axes: tuple[int, ...] | None = None
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def maximum_filter(
    input: object,
    size: object = ...,
    footprint: object = ...,
    output: object = ...,
    mode: object = ...,
    cval: object = ...,
    origin: object = ...,
    *,
    axes: object = ...,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def maximum_filter1d(
    input: object,
    size: object,
    axis: object = ...,
    output: object = ...,
    mode: object = ...,
    cval: object = ...,
    origin: object = ...,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def median_filter(
    input: object,
    size: object = ...,
    footprint: object = ...,
    output: object = ...,
    mode: object = ...,
    cval: object = ...,
    origin: object = ...,
    *,
    axes: object = ...,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def minimum_filter(
    input: object,
    size: object = ...,
    footprint: object = ...,
    output: object = ...,
    mode: object = ...,
    cval: object = ...,
    origin: object = ...,
    *,
    axes: object = ...,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def minimum_filter1d(
    input: object,
    size: object,
    axis: object = ...,
    output: object = ...,
    mode: object = ...,
    cval: object = ...,
    origin: object = ...,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def percentile_filter(
    input: object,
    percentile: object,
    size: object = ...,
    footprint: object = ...,
    output: object = ...,
    mode: object = ...,
    cval: object = ...,
    origin: object = ...,
    *,
    axes: object = ...,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def rank_filter(
    input: object,
    rank: object,
    size: object = ...,
    footprint: object = ...,
    output: object = ...,
    mode: object = ...,
    cval: object = ...,
    origin: object = ...,
    *,
    axes: object = ...,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def uniform_filter(
    input: object,
    size: object = ...,
    output: object = ...,
    mode: object = ...,
    cval: object = ...,
    origin: object = ...,
    *,
    axes: object = ...,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def uniform_filter1d(
    input: object,
    size: object,
    axis: object = ...,
    output: object = ...,
    mode: object = ...,
    cval: object = ...,
    origin: object = ...,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def prewitt(input: object, axis: object = ..., output: object = ..., mode: object = ..., cval: object = ...) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def sobel(input: object, axis: object = ..., output: object = ..., mode: object = ..., cval: object = ...) -> object: ...
