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
    output: object = None,
    mode: object = "reflect",
    cval: object = 0.0,
    origin: object = 0,
    *,
    axes: tuple[int, ...] | None = None,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def convolve1d(
    input: object,
    weights: object,
    axis: object = -1,
    output: object = None,
    mode: object = "reflect",
    cval: object = 0.0,
    origin: object = 0,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def correlate(
    input: object,
    weights: object,
    output: object = None,
    mode: object = "reflect",
    cval: object = 0.0,
    origin: object = 0,
    *,
    axes: tuple[int, ...] | None = None,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def correlate1d(
    input: object,
    weights: object,
    axis: object = -1,
    output: object = None,
    mode: object = "reflect",
    cval: object = 0.0,
    origin: object = 0,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def gaussian_filter(
    input: object,
    sigma: object,
    order: object = 0,
    output: object = None,
    mode: object = "reflect",
    cval: object = 0.0,
    truncate: object = 4.0,
    *,
    radius: object = None,
    axes: object = None,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def gaussian_filter1d(
    input: object,
    sigma: object,
    axis: object = -1,
    order: object = 0,
    output: object = None,
    mode: object = "reflect",
    cval: object = 0.0,
    truncate: object = 4.0,
    *,
    radius: object = None,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def gaussian_gradient_magnitude(
    input: object,
    sigma: object,
    output: object = None,
    mode: object = "reflect",
    cval: object = 0.0,
    *,
    axes: tuple[int, ...] | None = None,
    **kwargs: object,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def gaussian_laplace(
    input: object,
    sigma: object,
    output: object = None,
    mode: object = "reflect",
    cval: object = 0.0,
    *,
    axes: tuple[int, ...] | None = None,
    **kwargs: object,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def generic_filter(
    input: object,
    function: object,
    size: object = None,
    footprint: object = None,
    output: object = None,
    mode: object = "reflect",
    cval: object = 0.0,
    origin: object = 0,
    extra_arguments: object = (),
    extra_keywords: object = None,
    *,
    axes: tuple[int, ...] | None = None,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def generic_filter1d(
    input: object,
    function: object,
    filter_size: object,
    axis: object = -1,
    output: object = None,
    mode: object = "reflect",
    cval: object = 0.0,
    origin: object = 0,
    extra_arguments: object = (),
    extra_keywords: object = None,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def generic_gradient_magnitude(
    input: object,
    derivative: object,
    output: object = None,
    mode: object = "reflect",
    cval: object = 0.0,
    extra_arguments: object = (),
    extra_keywords: object = None,
    *,
    axes: tuple[int, ...] | None = None,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def generic_laplace(
    input: object,
    derivative2: object,
    output: object = None,
    mode: object = "reflect",
    cval: object = 0.0,
    extra_arguments: object = (),
    extra_keywords: object = None,
    *,
    axes: tuple[int, ...] | None = None,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def laplace(
    input: object, output: object = None, mode: object = "reflect", cval: object = 0.0, *, axes: tuple[int, ...] | None = None
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def maximum_filter(
    input: object,
    size: object = None,
    footprint: object = None,
    output: object = None,
    mode: object = "reflect",
    cval: object = 0.0,
    origin: object = 0,
    *,
    axes: object = None,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def maximum_filter1d(
    input: object,
    size: object,
    axis: object = -1,
    output: object = None,
    mode: object = "reflect",
    cval: object = 0.0,
    origin: object = 0,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def median_filter(
    input: object,
    size: object = None,
    footprint: object = None,
    output: object = None,
    mode: object = "reflect",
    cval: object = 0.0,
    origin: object = 0,
    *,
    axes: object = None,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def minimum_filter(
    input: object,
    size: object = None,
    footprint: object = None,
    output: object = None,
    mode: object = "reflect",
    cval: object = 0.0,
    origin: object = 0,
    *,
    axes: object = None,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def minimum_filter1d(
    input: object,
    size: object,
    axis: object = -1,
    output: object = None,
    mode: object = "reflect",
    cval: object = 0.0,
    origin: object = 0,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def percentile_filter(
    input: object,
    percentile: object,
    size: object = None,
    footprint: object = None,
    output: object = None,
    mode: object = "reflect",
    cval: object = 0.0,
    origin: object = 0,
    *,
    axes: object = None,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def rank_filter(
    input: object,
    rank: object,
    size: object = None,
    footprint: object = None,
    output: object = None,
    mode: object = "reflect",
    cval: object = 0.0,
    origin: object = 0,
    *,
    axes: object = None,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def uniform_filter(
    input: object,
    size: object = 3,
    output: object = None,
    mode: object = "reflect",
    cval: object = 0.0,
    origin: object = 0,
    *,
    axes: object = None,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def uniform_filter1d(
    input: object,
    size: object,
    axis: object = -1,
    output: object = None,
    mode: object = "reflect",
    cval: object = 0.0,
    origin: object = 0,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def prewitt(input: object, axis: object = -1, output: object = None, mode: object = "reflect", cval: object = 0.0) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def sobel(input: object, axis: object = -1, output: object = None, mode: object = "reflect", cval: object = 0.0) -> object: ...
