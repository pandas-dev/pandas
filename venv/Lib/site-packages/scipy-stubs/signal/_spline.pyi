from typing import overload

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc

###

@overload  # f64
def sepfir2d(
    input: onp.Array2D[npc.floating64 | npc.integer],
    hrow: onp.Array1D[npc.floating | npc.integer],
    hcol: onp.Array1D[npc.floating | npc.integer],
) -> onp.Array2D[np.float64]: ...
@overload  # f32
def sepfir2d(
    input: onp.Array2D[npc.floating32 | npc.floating16],
    hrow: onp.Array1D[npc.floating32 | npc.floating16],
    hcol: onp.Array1D[npc.floating32 | npc.floating16],
) -> onp.Array2D[np.float32]: ...
@overload  # c128
def sepfir2d(
    input: onp.Array2D[npc.complexfloating128], hrow: onp.Array1D[npc.number], hcol: onp.Array1D[npc.number]
) -> onp.Array2D[np.complex128]: ...
@overload  # c64
def sepfir2d(
    input: onp.Array2D[npc.complexfloating64],
    hrow: onp.Array1D[npc.inexact32 | np.float16],
    hcol: onp.Array1D[npc.inexact32 | np.float16],
) -> onp.Array2D[np.complex64]: ...
@overload  # fallback
def sepfir2d(input: onp.Array2D[npc.number], hrow: onp.Array1D[npc.number], hcol: onp.Array1D[npc.number]) -> onp.Array2D: ...

# undocumented
@overload
def symiirorder1_ic(
    input: onp.ArrayND[npc.floating64 | npc.integer, tuple[int] | tuple[int, int]], z1: float, precision: float = -1.0, /
) -> onp.Array2D[np.float64]: ...
@overload
def symiirorder1_ic(
    input: onp.ArrayND[npc.floating32 | npc.floating16, tuple[int] | tuple[int, int]], z1: float, precision: float = -1.0, /
) -> onp.Array2D[np.float32]: ...
@overload
def symiirorder1_ic(
    input: onp.ArrayND[npc.complexfloating128, tuple[int] | tuple[int, int]], z1: float, precision: float = -1.0, /
) -> onp.Array2D[np.complex128]: ...
@overload
def symiirorder1_ic(
    input: onp.ArrayND[npc.complexfloating64, tuple[int] | tuple[int, int]], z1: float, precision: float = -1.0, /
) -> onp.Array2D[np.complex64]: ...

# undocumented
@overload
def symiirorder2_ic_fwd(
    input: onp.ArrayND[npc.floating64 | npc.integer, tuple[int] | tuple[int, int]],
    r: float,
    omega: float,
    precision: float = -1.0,
    /,
) -> onp.Array2D[np.float64]: ...
@overload
def symiirorder2_ic_fwd(
    input: onp.ArrayND[npc.floating32 | npc.floating16, tuple[int] | tuple[int, int]],
    r: float,
    omega: float,
    precision: float = -1.0,
    /,
) -> onp.Array2D[np.float32]: ...
@overload
def symiirorder2_ic_fwd(
    input: onp.ArrayND[npc.complexfloating128, tuple[int] | tuple[int, int]], r: float, omega: float, precision: float = -1.0, /
) -> onp.Array2D[np.complex128]: ...
@overload
def symiirorder2_ic_fwd(
    input: onp.ArrayND[npc.complexfloating64, tuple[int] | tuple[int, int]], r: float, omega: float, precision: float = -1.0, /
) -> onp.Array2D[np.complex64]: ...

# undocumented
@overload
def symiirorder2_ic_bwd(
    input: onp.ArrayND[npc.floating64 | npc.integer, tuple[int] | tuple[int, int]],
    r: float,
    omega: float,
    precision: float = -1.0,
    /,
) -> onp.Array2D[np.float64]: ...
@overload
def symiirorder2_ic_bwd(
    input: onp.ArrayND[npc.floating32 | npc.floating16, tuple[int] | tuple[int, int]],
    r: float,
    omega: float,
    precision: float = -1.0,
    /,
) -> onp.Array2D[np.float32]: ...
@overload
def symiirorder2_ic_bwd(
    input: onp.ArrayND[npc.complexfloating128, tuple[int] | tuple[int, int]], r: float, omega: float, precision: float = -1.0, /
) -> onp.Array2D[np.complex128]: ...
@overload
def symiirorder2_ic_bwd(
    input: onp.ArrayND[npc.complexfloating64, tuple[int] | tuple[int, int]], r: float, omega: float, precision: float = -1.0, /
) -> onp.Array2D[np.complex64]: ...
