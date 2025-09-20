class error(Exception): ...

# _qagie(fun, bound, inf,                   | args, full_output, epsabs, epsrel, limit) -> {out}
# _qagse(fun, a, b,                         | args, full_output, epsabs, epsrel, limit) -> {out}
# _qagpe(fun, a, b, points,                 | args, full_output, epsabs, epsrel, limit) -> {out}
# _qawce(fun, a, b, c,                      | args, full_output, epsabs, epsrel, limit) -> {out}
# _qawse(fun, a, b, (alpha, beta), integr,  | args, full_output, epsabs, epsrel, limit) -> {out}
# _qawfe(fun, a, omega, integr,             | args, full_output, epsabs, limlst, limit, maxp1) -> {out}
# _qawoe(fun, a, b, omega, integr,          | args, full_output, epsabs, epsrel, limit, maxp1, icall, momcom, chebmo) -> {out}

# {out} := [result, abserr, infodict, ier]
