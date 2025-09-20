from . import dop, lsoda, odepack, quadpack, vode  # deprecated namespaces
from ._bvp import solve_bvp
from ._cubature import cubature
from ._ivp import BDF, DOP853, LSODA, RK23, RK45, DenseOutput, OdeSolution, OdeSolver, Radau, solve_ivp
from ._lebedev import lebedev_rule
from ._ode import complex_ode, ode
from ._odepack_py import ODEintWarning, odeint
from ._quad_vec import quad_vec
from ._quadpack_py import IntegrationWarning, dblquad, nquad, quad, tplquad
from ._quadrature import cumulative_simpson, cumulative_trapezoid, fixed_quad, newton_cotes, qmc_quad, romb, simpson, trapezoid
from ._tanhsinh import nsum, tanhsinh

__all__ = [
    "BDF",
    "DOP853",
    "LSODA",
    "RK23",
    "RK45",
    "DenseOutput",
    "IntegrationWarning",
    "ODEintWarning",
    "OdeSolution",
    "OdeSolver",
    "Radau",
    "complex_ode",
    "cubature",
    "cumulative_simpson",
    "cumulative_trapezoid",
    "dblquad",
    "dop",
    "fixed_quad",
    "lebedev_rule",
    "lsoda",
    "newton_cotes",
    "nquad",
    "nsum",
    "ode",
    "odeint",
    "odepack",
    "qmc_quad",
    "quad",
    "quad_vec",
    "quadpack",
    "romb",
    "simpson",
    "solve_bvp",
    "solve_ivp",
    "tanhsinh",
    "tplquad",
    "trapezoid",
    "vode",
]
