# This file is not meant for public use and will be removed in SciPy v2.0.0.

from typing_extensions import deprecated

import numpy as np

__all__ = [
    "ai_zeros",
    "assoc_laguerre",
    "bei_zeros",
    "beip_zeros",
    "ber_zeros",
    "bernoulli",
    "berp_zeros",
    "bi_zeros",
    "comb",
    "digamma",
    "diric",
    "erf_zeros",
    "euler",
    "factorial",
    "factorial2",
    "factorialk",
    "fresnel_zeros",
    "fresnelc_zeros",
    "fresnels_zeros",
    "gamma",
    "h1vp",
    "h2vp",
    "hankel1",
    "hankel2",
    "iv",
    "ivp",
    "jn_zeros",
    "jnjnp_zeros",
    "jnp_zeros",
    "jnyn_zeros",
    "jv",
    "jvp",
    "kei_zeros",
    "keip_zeros",
    "kelvin_zeros",
    "ker_zeros",
    "kerp_zeros",
    "kv",
    "kvp",
    "lmbda",
    "lqmn",
    "lqn",
    "mathieu_a",
    "mathieu_b",
    "mathieu_even_coef",
    "mathieu_odd_coef",
    "obl_cv_seq",
    "pbdn_seq",
    "pbdv_seq",
    "pbvv_seq",
    "perm",
    "polygamma",
    "pro_cv_seq",
    "psi",
    "riccati_jn",
    "riccati_yn",
    "sinc",
    "y0_zeros",
    "y1_zeros",
    "y1p_zeros",
    "yn_zeros",
    "ynp_zeros",
    "yv",
    "yvp",
    "zeta",
]

@deprecated("will be removed in SciPy v2.0.0")
def diric(x: object, n: object) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def jnjnp_zeros(nt: object) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def jnyn_zeros(n: object, nt: object) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def jn_zeros(n: object, nt: object) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def jnp_zeros(n: object, nt: object) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def yn_zeros(n: object, nt: object) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def ynp_zeros(n: object, nt: object) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def y0_zeros(nt: object, complex: object = False) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def y1_zeros(nt: object, complex: object = False) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def y1p_zeros(nt: object, complex: object = False) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def jvp(v: object, z: object, n: object = 1) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def yvp(v: object, z: object, n: object = 1) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def kvp(v: object, z: object, n: object = 1) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def ivp(v: object, z: object, n: object = 1) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def h1vp(v: object, z: object, n: object = 1) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def h2vp(v: object, z: object, n: object = 1) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def riccati_jn(n: object, x: object) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def riccati_yn(n: object, x: object) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def erf_zeros(nt: object) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def fresnelc_zeros(nt: object) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def fresnels_zeros(nt: object) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def fresnel_zeros(nt: object) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def assoc_laguerre(x: object, n: object, k: object = 0.0) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def polygamma(n: object, x: object) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def mathieu_even_coef(m: object, q: object) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def mathieu_odd_coef(m: object, q: object) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def lqmn(m: object, n: object, z: object) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def bernoulli(n: object) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def euler(n: object) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def lqn(n: object, z: object) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def ai_zeros(nt: object) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def bi_zeros(nt: object) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def lmbda(v: object, x: object) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def pbdv_seq(v: object, x: object) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def pbvv_seq(v: object, x: object) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def pbdn_seq(n: object, z: object) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def ber_zeros(nt: object) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def bei_zeros(nt: object) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def ker_zeros(nt: object) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def kei_zeros(nt: object) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def berp_zeros(nt: object) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def beip_zeros(nt: object) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def kerp_zeros(nt: object) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def keip_zeros(nt: object) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def kelvin_zeros(nt: object) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def pro_cv_seq(m: object, n: object, c: object) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def obl_cv_seq(m: object, n: object, c: object) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def comb(N: object, k: object, *, exact: object = False, repetition: object = False) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def perm(N: object, k: object, exact: object = False) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def factorial(n: object, exact: object = False, extend: object = "zero") -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def factorial2(n: object, exact: object = False, extend: object = "zero") -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def factorialk(n: object, k: object, exact: object = False, extend: object = "zero") -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def zeta(x: object, q: object = None, out: object = None) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def sinc(x: object) -> object: ...

# Deprecated
gamma: np.ufunc
hankel1: np.ufunc
hankel2: np.ufunc
iv: np.ufunc
jv: np.ufunc
kv: np.ufunc
mathieu_a: np.ufunc
mathieu_b: np.ufunc
psi: np.ufunc
yv: np.ufunc
digamma: np.ufunc
