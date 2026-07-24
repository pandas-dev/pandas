# This file is not meant for public use and will be removed in SciPy v2.0.0.

from typing_extensions import deprecated

import numpy as np

__all__ = [
    "airy",
    "c_roots",
    "cg_roots",
    "chebyc",
    "chebys",
    "chebyt",
    "chebyu",
    "gegenbauer",
    "genlaguerre",
    "h_roots",
    "he_roots",
    "hermite",
    "hermitenorm",
    "j_roots",
    "jacobi",
    "js_roots",
    "l_roots",
    "la_roots",
    "laguerre",
    "legendre",
    "p_roots",
    "ps_roots",
    "roots_chebyc",
    "roots_chebys",
    "roots_chebyt",
    "roots_chebyu",
    "roots_gegenbauer",
    "roots_genlaguerre",
    "roots_hermite",
    "roots_hermitenorm",
    "roots_jacobi",
    "roots_laguerre",
    "roots_legendre",
    "roots_sh_chebyt",
    "roots_sh_chebyu",
    "roots_sh_jacobi",
    "roots_sh_legendre",
    "s_roots",
    "sh_chebyt",
    "sh_chebyu",
    "sh_jacobi",
    "sh_legendre",
    "t_roots",
    "ts_roots",
    "u_roots",
    "us_roots",
]

airy: np.ufunc

@deprecated("will be removed in SciPy v2.0.0")
def roots_jacobi(n: object, alpha: object, beta: object, mu: object = False) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def roots_sh_jacobi(n: object, p1: object, q1: object, mu: object = False) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def roots_genlaguerre(n: object, alpha: object, mu: object = False) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def roots_laguerre(n: object, mu: object = False) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def roots_hermite(n: object, mu: object = False) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def roots_hermitenorm(n: object, mu: object = False) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def roots_gegenbauer(n: object, alpha: object, mu: object = False) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def roots_chebyt(n: object, mu: object = False) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def roots_chebyu(n: object, mu: object = False) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def roots_chebyc(n: object, mu: object = False) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def roots_chebys(n: object, mu: object = False) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def roots_sh_chebyt(n: object, mu: object = False) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def roots_sh_chebyu(n: object, mu: object = False) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def roots_legendre(n: object, mu: object = False) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def roots_sh_legendre(n: object, mu: object = False) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def legendre(n: object, monic: object = False) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def chebyt(n: object, monic: object = False) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def chebyu(n: object, monic: object = False) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def chebyc(n: object, monic: object = False) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def chebys(n: object, monic: object = False) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def jacobi(n: object, alpha: object, beta: object, monic: object = False) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def laguerre(n: object, monic: object = False) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def genlaguerre(n: object, alpha: object, monic: object = False) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def hermite(n: object, monic: object = False) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def hermitenorm(n: object, monic: object = False) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def gegenbauer(n: object, alpha: object, monic: object = False) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def sh_legendre(n: object, monic: object = False) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def sh_chebyt(n: object, monic: object = False) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def sh_chebyu(n: object, monic: object = False) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def sh_jacobi(n: object, p: object, q: object, monic: object = False) -> object: ...

# pyright: reportDeprecated=false
p_roots = roots_legendre
t_roots = roots_chebyt
u_roots = roots_chebyu
c_roots = roots_chebyc
s_roots = roots_chebys
j_roots = roots_jacobi
l_roots = roots_laguerre
la_roots = roots_genlaguerre
h_roots = roots_hermite
he_roots = roots_hermitenorm
cg_roots = roots_gegenbauer
ps_roots = roots_sh_legendre
ts_roots = roots_sh_chebyt
us_roots = roots_sh_chebyu
js_roots = roots_sh_jacobi
