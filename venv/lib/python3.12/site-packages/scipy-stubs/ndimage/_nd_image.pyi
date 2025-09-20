# defined in scipy/ndimage/src/nd_image.c

from _typeshed import Incomplete
from collections.abc import Callable, Mapping
from typing import overload
from typing_extensions import CapsuleType

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc

from scipy._lib._ccallback import LowLevelCallable

def correlate1d(
    input: onp.ArrayND[npc.number],
    weights: onp.ArrayND[np.float64],
    axis: int,
    output: onp.ArrayND[npc.number],
    mode: int,
    cval: float,
    origin: int,
    /,
) -> None: ...
def correlate(
    input: onp.ArrayND[npc.number],
    weights: onp.ArrayND[np.float64],
    output: onp.ArrayND[npc.number],
    mode: int,
    cval: float,
    origin: int | onp.ArrayND[np.intp],
    /,
) -> None: ...
def uniform_filter1d(
    input: onp.ArrayND[npc.number],
    filter_size: int,
    axis: int,
    output: onp.ArrayND[npc.number],
    mode: int,
    cval: float,
    origin: int,
    /,
) -> None: ...
def min_or_max_filter1d(
    input: onp.ArrayND[npc.number],
    filter_size: int,
    axis: int,
    output: onp.ArrayND[npc.number],
    mode: int,
    cval: float,
    origin: int,
    minimum: int,
    /,
) -> None: ...
def min_or_max_filter(
    input: onp.ArrayND[npc.number],
    footprint: onp.ArrayND[npc.number],
    structure: onp.ArrayND[npc.number] | None,
    output: onp.ArrayND[npc.number],
    mode: int,
    cval: float,
    origin: int | onp.ArrayND[np.intp],
    minimum: int,
    /,
) -> None: ...
def rank_filter(
    input: onp.ArrayND[npc.number],
    rank: int,
    footprint: onp.ArrayND[npc.number],
    output: onp.ArrayND[npc.number],
    mode: int,
    cval: float,
    origin: onp.ArrayND[np.intp],
    /,
) -> None: ...
def generic_filter1d(
    input: onp.ArrayND[npc.number],
    fnc: Callable[..., Incomplete] | LowLevelCallable,
    filter_size: int,
    axis: int,
    output: onp.ArrayND[npc.number],
    mode: int,
    cval: float,
    origin: int,
    extra_arguments: tuple[object, ...],
    extra_keywords: Mapping[str, object],
    /,
) -> None: ...
def generic_filter(
    input: onp.ArrayND[npc.number],
    fnc: Callable[..., Incomplete] | LowLevelCallable,
    footprint: onp.ArrayND[npc.number],
    output: onp.ArrayND[npc.number],
    mode: int,
    cval: float,
    origin: int | onp.ArrayND[np.intp],
    extra_arguments: tuple[object, ...],
    extra_keywords: Mapping[str, object],
    /,
) -> None: ...
def fourier_filter(
    input: onp.ArrayND[npc.number],
    parameters: onp.ArrayND[npc.number],
    n: int,
    axis: int,
    output: onp.ArrayND[npc.number],
    filter_type: int,
    /,
) -> None: ...
def fourier_shift(
    input: onp.ArrayND[npc.number], order: int, axis: int, output: onp.ArrayND[npc.number], mode: int, /
) -> None: ...
def spline_filter1d(
    input: onp.ArrayND[npc.number], order: int, axis: int, output: onp.ArrayND[npc.number], mode: int, /
) -> None: ...
def geometric_transform(
    input: onp.ArrayND[npc.number],
    fnc: Callable[..., Incomplete] | LowLevelCallable,
    coordinates: onp.ArrayND[npc.number] | None,
    matrix: onp.ArrayND[npc.number] | None,
    shift: onp.ArrayND[npc.number] | None,
    output: onp.ArrayND[npc.number],
    order: int,
    mode: int,
    cval: float,
    nprepad: int,
    extra_arguments: tuple[object, ...],
    extra_keywords: Mapping[str, object],
    /,
) -> None: ...
def zoom_shift(
    input: onp.ArrayND[npc.number],
    zoom: onp.ArrayND[npc.number] | None,
    shift: onp.ArrayND[npc.number] | None,
    output: onp.ArrayND[npc.number],
    order: int,
    mode: int,
    cval: float,
    nprepad: int,
    grid_mode: int,
    /,
) -> None: ...

#
def find_objects(input: onp.ArrayND[npc.number], max_label: int) -> None: ...
def value_indices(arr: onp.ArrayND[npc.integer], ignoreValIsNone: bool, ignorevalArr: onp.ArrayND, /) -> None: ...

#
def watershed_ift(
    input: onp.ArrayND[npc.number],
    markers: onp.ArrayND[npc.number],
    strct: onp.ArrayND[npc.number],
    output: onp.ArrayND[npc.number],
) -> None: ...
def distance_transform_bf(
    input: onp.ArrayND[npc.number],
    metric: int,
    sampling: onp.ArrayND[npc.number] | None,
    output: onp.ArrayND[npc.number] | None,
    features: onp.ArrayND[npc.number] | None,
) -> None: ...
def distance_transform_op(
    input: onp.ArrayND[npc.number],
    strct: onp.ArrayND[npc.number],
    distances: onp.ArrayND[npc.number],
    features: onp.ArrayND[npc.number] | None,
) -> None: ...
def euclidean_feature_transform(
    input: onp.ArrayND[npc.number], sampling: onp.ArrayND[npc.number] | None, features: onp.ArrayND[npc.number]
) -> None: ...

#
@overload
def binary_erosion(
    input: onp.ArrayND[npc.number],
    strct: onp.ArrayND[npc.number],
    mask: onp.ArrayND[npc.number] | None,
    output: onp.ArrayND[npc.number],
    border_value: int,
    origin: onp.ArrayND[np.intp],
    invert: int,
    center_is_true: int,
    return_coordinates: onp.ToFalse,
) -> int: ...
@overload
def binary_erosion(
    input: onp.ArrayND[npc.number],
    strct: onp.ArrayND[npc.number],
    mask: onp.ArrayND[npc.number] | None,
    output: onp.ArrayND[npc.number],
    border_value: int,
    origin: int | onp.ArrayND[np.intp],
    invert: int,
    center_is_true: int,
    return_coordinates: onp.ToTrue,
) -> tuple[int, CapsuleType]: ...

#
def binary_erosion2(
    array: onp.ArrayND[npc.number],
    strct: onp.ArrayND[npc.number],
    mask: onp.ArrayND[npc.number] | None,
    niter: int,
    origin: int | onp.ArrayND[np.intp],
    invert: int,
    cobj: CapsuleType | None,
) -> None: ...
