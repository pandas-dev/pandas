from collections.abc import Iterable
from typing import Literal, overload

import optype.numpy as onp
import optype.numpy.compat as npc
import optype.typing as opt

__all__ = ["cossin"]

type _Tuple2[T] = tuple[T, T]
type _Tuple3[T] = tuple[T, T, T]

type _Float1D = onp.Array1D[npc.floating]
type _Float2D = onp.Array2D[npc.floating]
type _FloatND = onp.ArrayND[npc.floating]
type _Inexact2D = onp.Array2D[npc.inexact]
type _InexactND = onp.ArrayND[npc.inexact]

@overload  # (float[:, :], separate=False) -> float[:, :]**3
def cossin(
    X: onp.ToFloatStrict2D | Iterable[onp.ToFloatStrict2D],
    p: opt.AnyInt | None = None,
    q: opt.AnyInt | None = None,
    separate: Literal[False] = False,
    swap_sign: bool = False,
    compute_u: bool = True,
    compute_vh: bool = True,
) -> _Tuple3[_Float2D]: ...
@overload  # (float[:, :, ...], separate=False) -> float[:, :, ...]**3
def cossin(
    X: onp.ToFloatND | Iterable[onp.ToFloatND],
    p: opt.AnyInt | None = None,
    q: opt.AnyInt | None = None,
    separate: Literal[False] = False,
    swap_sign: bool = False,
    compute_u: bool = True,
    compute_vh: bool = True,
) -> _Tuple3[_FloatND]: ...
@overload  # (float[:, :], *, separate=True) -> (float[:, :]**2, float[:], float[:, :]**2)
def cossin(
    X: onp.ToFloatStrict2D | Iterable[onp.ToFloatStrict2D],
    p: opt.AnyInt | None = None,
    q: opt.AnyInt | None = None,
    *,
    separate: Literal[True],
    swap_sign: bool = False,
    compute_u: bool = True,
    compute_vh: bool = True,
) -> tuple[_Tuple2[_Float2D], _Float1D, _Tuple2[_Float2D]]: ...
@overload  # (float[:, :, ...], *, separate=True) -> (float[:, :, ...]**2, float[:, ...], float[:, :, ...]**2)
def cossin(
    X: onp.ToFloatND | Iterable[onp.ToFloatND],
    p: opt.AnyInt | None = None,
    q: opt.AnyInt | None = None,
    *,
    separate: Literal[True],
    swap_sign: bool = False,
    compute_u: bool = True,
    compute_vh: bool = True,
) -> tuple[_Tuple2[_FloatND], _FloatND, _Tuple2[_FloatND]]: ...
@overload  # (complex[:, :], separate=False) -> complex[:, :]**3
def cossin(
    X: onp.ToComplexStrict2D | Iterable[onp.ToComplexStrict2D],
    p: opt.AnyInt | None = None,
    q: opt.AnyInt | None = None,
    separate: Literal[False] = False,
    swap_sign: bool = False,
    compute_u: bool = True,
    compute_vh: bool = True,
) -> _Tuple3[_Inexact2D]: ...
@overload  # (complex[:, :, ...], separate=False) -> complex[:, :, ...]**3
def cossin(
    X: onp.ToComplexND | Iterable[onp.ToComplexND],
    p: opt.AnyInt | None = None,
    q: opt.AnyInt | None = None,
    separate: Literal[False] = False,
    swap_sign: bool = False,
    compute_u: bool = True,
    compute_vh: bool = True,
) -> _Tuple3[_InexactND]: ...
@overload  # (complex[:, :], separate=True) -> (complex[:, :]**2, float[:], complex[:, :]**2)
def cossin(
    X: onp.ToComplexStrict2D | Iterable[onp.ToComplexStrict2D],
    p: opt.AnyInt | None = None,
    q: opt.AnyInt | None = None,
    *,
    separate: Literal[True],
    swap_sign: bool = False,
    compute_u: bool = True,
    compute_vh: bool = True,
) -> tuple[_Tuple2[_Inexact2D], _Float1D, _Tuple2[_Inexact2D]]: ...
@overload  # (complex[:, :, ...], separate=True) -> (complex[:, :, ...]**2, float[:, ...], complex[:, :, ...]**2)
def cossin(
    X: onp.ToComplexND | Iterable[onp.ToComplexND],
    p: opt.AnyInt | None = None,
    q: opt.AnyInt | None = None,
    *,
    separate: Literal[True],
    swap_sign: bool = False,
    compute_u: bool = True,
    compute_vh: bool = True,
) -> tuple[_Tuple2[_InexactND], _FloatND, _Tuple2[_InexactND]]: ...
