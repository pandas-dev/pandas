from typing import Final, Literal

import numpy as np

from ._ufuncs import (
    airy as airy,
    airye as airye,
    bei as bei,
    beip as beip,
    ber as ber,
    berp as berp,
    besselpoly as besselpoly,
    beta as beta,
    betaln as betaln,
    binom as binom,
    cbrt as cbrt,
    cosdg as cosdg,
    cosm1 as cosm1,
    cotdg as cotdg,
    dawsn as dawsn,
    ellipe as ellipe,
    ellipeinc as ellipeinc,
    ellipj as ellipj,
    ellipk as ellipk,
    ellipkinc as ellipkinc,
    ellipkm1 as ellipkm1,
    erf as erf,
    erfc as erfc,
    erfcx as erfcx,
    erfi as erfi,
    exp1 as exp1,
    exp2 as exp2,
    exp10 as exp10,
    expi as expi,
    expit as expit,
    expm1 as expm1,
    exprel as exprel,
    fresnel as fresnel,
    gamma as gamma,
    gammainc as gammainc,
    gammaincc as gammaincc,
    gammainccinv as gammainccinv,
    gammaincinv as gammaincinv,
    gammaln as gammaln,
    gammasgn as gammasgn,
    hankel1 as hankel1,
    hankel1e as hankel1e,
    hankel2 as hankel2,
    hankel2e as hankel2e,
    hyp2f1 as hyp2f1,
    i0 as i0,
    i0e as i0e,
    i1 as i1,
    i1e as i1e,
    it2i0k0 as it2i0k0,
    it2j0y0 as it2j0y0,
    it2struve0 as it2struve0,
    itairy as itairy,
    iti0k0 as iti0k0,
    itj0y0 as itj0y0,
    itmodstruve0 as itmodstruve0,
    itstruve0 as itstruve0,
    iv as iv,
    ive as ive,
    j0 as j0,
    j1 as j1,
    jv as jv,
    jve as jve,
    k0 as k0,
    k0e as k0e,
    k1 as k1,
    k1e as k1e,
    kei as kei,
    keip as keip,
    kelvin as kelvin,
    ker as ker,
    kerp as kerp,
    kv as kv,
    kve as kve,
    log1p as log1p,
    log_expit as log_expit,
    log_ndtr as log_ndtr,
    log_wright_bessel as log_wright_bessel,
    loggamma as loggamma,
    logit as logit,
    mathieu_a as mathieu_a,
    mathieu_b as mathieu_b,
    mathieu_cem as mathieu_cem,
    mathieu_modcem1 as mathieu_modcem1,
    mathieu_modcem2 as mathieu_modcem2,
    mathieu_modsem1 as mathieu_modsem1,
    mathieu_modsem2 as mathieu_modsem2,
    mathieu_sem as mathieu_sem,
    modfresnelm as modfresnelm,
    modfresnelp as modfresnelp,
    modstruve as modstruve,
    ndtr as ndtr,
    obl_ang1 as obl_ang1,
    obl_ang1_cv as obl_ang1_cv,
    obl_cv as obl_cv,
    obl_rad1 as obl_rad1,
    obl_rad1_cv as obl_rad1_cv,
    obl_rad2 as obl_rad2,
    obl_rad2_cv as obl_rad2_cv,
    pbdv as pbdv,
    pbvv as pbvv,
    pbwa as pbwa,
    pro_ang1 as pro_ang1,
    pro_ang1_cv as pro_ang1_cv,
    pro_cv as pro_cv,
    pro_rad1 as pro_rad1,
    pro_rad1_cv as pro_rad1_cv,
    pro_rad2 as pro_rad2,
    pro_rad2_cv as pro_rad2_cv,
    psi as psi,
    radian as radian,
    rgamma as rgamma,
    sindg as sindg,
    struve as struve,
    tandg as tandg,
    voigt_profile as voigt_profile,
    wofz as wofz,
    wright_bessel as wright_bessel,
    xlog1py as xlog1py,
    xlogy as xlogy,
    y0 as y0,
    y1 as y1,
    yv as yv,
    yve as yve,
    zetac as zetac,
)

_cospi: Final[np.ufunc] = ...  # undocumented
_iv_ratio: Final[np.ufunc] = ...  # undocumented
_iv_ratio_c: Final[np.ufunc] = ...  # undocumented
_lambertw: Final[np.ufunc] = ...  # undocumented
_log1mexp: Final[np.ufunc] = ...  # undocumented
_log1pmx: Final[np.ufunc] = ...  # undocumented
_riemann_zeta: Final[np.ufunc] = ...  # undocumented
_scaled_exp1: Final[np.ufunc] = ...  # undocumented
_sinpi: Final[np.ufunc] = ...  # undocumented
_spherical_in: Final[np.ufunc] = ...  # undocumented
_spherical_in_d: Final[np.ufunc] = ...  # undocumented
_spherical_jn: Final[np.ufunc] = ...  # undocumented
_spherical_jn_d: Final[np.ufunc] = ...  # undocumented
_spherical_kn: Final[np.ufunc] = ...  # undocumented
_spherical_kn_d: Final[np.ufunc] = ...  # undocumented
_spherical_yn: Final[np.ufunc] = ...  # undocumented
_spherical_yn_d: Final[np.ufunc] = ...  # undocumented
_zeta: Final[np.ufunc] = ...  # undocumented

assoc_legendre_p: Final[dict[tuple[bool, Literal[0, 1, 2]], np.ufunc]] = ...
legendre_p: Final[tuple[np.ufunc, np.ufunc, np.ufunc]] = ...
sph_harm_y: Final[tuple[np.ufunc, np.ufunc, np.ufunc]] = ...
sph_legendre_p: Final[tuple[np.ufunc, np.ufunc, np.ufunc]] = ...

def _set_action(code: int, action: int, /) -> None: ...  # undocumented
