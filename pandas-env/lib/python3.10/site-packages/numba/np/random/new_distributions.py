"""
Algorithmic implementations for generating different types
of random distributions.
"""

import numpy as np

from numba.core.extending import register_jitable
from numba.np.random._constants import (wi_double, ki_double,
                                        ziggurat_nor_r, fi_double,
                                        wi_float, ki_float,
                                        ziggurat_nor_inv_r_f,
                                        ziggurat_nor_r_f, fi_float,
                                        we_double, ke_double,
                                        ziggurat_exp_r, fe_double,
                                        we_float, ke_float,
                                        ziggurat_exp_r_f, fe_float,
                                        INT64_MAX, ziggurat_nor_inv_r)
from numba.np.random.generator_core import (next_double, next_float,
                                            next_uint32, next_uint64)
# All of the following implementations are direct translations from:
# https://github.com/numpy/numpy/blob/7cfef93c77599bd387ecc6a15d186c5a46024dac/numpy/random/src/distributions/distributions.c


@register_jitable
def np_log1p(x):
    return np.log1p(x)


@register_jitable
def np_log1pf(x):
    return np.log1p(np.float32(x))


@register_jitable
def random_rayleigh(bitgen, mode):
    return mode * np.sqrt(2.0 * random_standard_exponential(bitgen))


@register_jitable
def np_expm1(x):
    return np.expm1(x)


@register_jitable
def random_standard_normal(bitgen):
    while 1:
        r = next_uint64(bitgen)
        idx = r & 0xff
        r >>= 8
        sign = r & 0x1
        rabs = (r >> 1) & 0x000fffffffffffff
        x = rabs * wi_double[idx]
        if (sign & 0x1):
            x = -x
        if rabs < ki_double[idx]:
            return x
        if idx == 0:
            while 1:
                xx = -ziggurat_nor_inv_r * np.log1p(-next_double(bitgen))
                yy = -np.log1p(-next_double(bitgen))
                if (yy + yy > xx * xx):
                    if ((rabs >> 8) & 0x1):
                        return -(ziggurat_nor_r + xx)
                    else:
                        return ziggurat_nor_r + xx
        else:
            if (((fi_double[idx - 1] - fi_double[idx]) *
                    next_double(bitgen) + fi_double[idx]) <
                    np.exp(-0.5 * x * x)):
                return x


@register_jitable
def random_standard_normal_f(bitgen):
    while 1:
        r = next_uint32(bitgen)
        idx = r & 0xff
        sign = (r >> 8) & 0x1
        rabs = (r >> 9) & 0x0007fffff
        x = np.float32(np.float32(rabs) * wi_float[idx])
        if (sign & 0x1):
            x = -x
        if (rabs < ki_float[idx]):
            return x
        if (idx == 0):
            while 1:
                xx = np.float32(-ziggurat_nor_inv_r_f *
                             np_log1pf(-next_float(bitgen)))
                yy = np.float32(-np_log1pf(-next_float(bitgen)))
                if (np.float32(yy + yy) > np.float32(xx * xx)):
                    if ((rabs >> 8) & 0x1):
                        return -np.float32(ziggurat_nor_r_f + xx)
                    else:
                        return np.float32(ziggurat_nor_r_f + xx)
        else:
            if (((fi_float[idx - 1] - fi_float[idx]) * next_float(bitgen) +
                 fi_float[idx]) < np.float32(np.exp(-np.float32(0.5) * x * x))):
                return x


@register_jitable
def random_standard_exponential(bitgen):
    while 1:
        ri = next_uint64(bitgen)
        ri >>= 3
        idx = ri & 0xFF
        ri >>= 8
        x = ri * we_double[idx]
        if (ri < ke_double[idx]):
            return x
        else:
            if idx == 0:
                return ziggurat_exp_r - np_log1p(-next_double(bitgen))
            elif ((fe_double[idx - 1] - fe_double[idx]) * next_double(bitgen) +
                  fe_double[idx] < np.exp(-x)):
                return x


@register_jitable
def random_standard_exponential_f(bitgen):
    while 1:
        ri = next_uint32(bitgen)
        ri >>= 1
        idx = ri & 0xFF
        ri >>= 8
        x = np.float32(np.float32(ri) * we_float[idx])
        if (ri < ke_float[idx]):
            return x
        else:
            if (idx == 0):
                return np.float32(ziggurat_exp_r_f -
                               np.float32(np_log1pf(-next_float(bitgen))))
            elif ((fe_float[idx - 1] - fe_float[idx]) * next_float(bitgen) +
                  fe_float[idx] < np.float32(np.exp(np.float32(-x)))):
                return x


@register_jitable
def random_standard_exponential_inv(bitgen):
    return -np_log1p(-next_double(bitgen))


@register_jitable
def random_standard_exponential_inv_f(bitgen):
    return -np.log(np.float32(1.0) - next_float(bitgen))


@register_jitable
def random_standard_gamma(bitgen, shape):
    if (shape == 1.0):
        return random_standard_exponential(bitgen)
    elif (shape == 0.0):
        return 0.0
    elif (shape < 1.0):
        while 1:
            U = next_double(bitgen)
            V = random_standard_exponential(bitgen)
            if (U <= 1.0 - shape):
                X = pow(U, 1. / shape)
                if (X <= V):
                    return X
            else:
                Y = -np.log((1 - U) / shape)
                X = pow(1.0 - shape + shape * Y, 1. / shape)
                if (X <= (V + Y)):
                    return X
    else:
        b = shape - 1. / 3.
        c = 1. / np.sqrt(9 * b)
        while 1:
            while 1:
                X = random_standard_normal(bitgen)
                V = 1.0 + c * X
                if (V > 0.0):
                    break

            V = V * V * V
            U = next_double(bitgen)
            if (U < 1.0 - 0.0331 * (X * X) * (X * X)):
                return (b * V)

            if (np.log(U) < 0.5 * X * X + b * (1. - V + np.log(V))):
                return (b * V)


@register_jitable
def random_standard_gamma_f(bitgen, shape):
    f32_one = np.float32(1.0)
    shape = np.float32(shape)
    if (shape == f32_one):
        return random_standard_exponential_f(bitgen)
    elif (shape == np.float32(0.0)):
        return np.float32(0.0)
    elif (shape < f32_one):
        while 1:
            U = next_float(bitgen)
            V = random_standard_exponential_f(bitgen)
            if (U <= f32_one - shape):
                X = np.float32(pow(U, np.float32(f32_one / shape)))
                if (X <= V):
                    return X
            else:
                Y = np.float32(-np.log(np.float32((f32_one - U) / shape)))
                X = np.float32(pow(f32_one - shape + np.float32(shape * Y),
                            np.float32(f32_one / shape)))
                if (X <= (V + Y)):
                    return X
    else:
        b = shape - f32_one / np.float32(3.0)
        c = np.float32(f32_one / np.float32(np.sqrt(np.float32(9.0) * b)))
        while 1:
            while 1:
                X = np.float32(random_standard_normal_f(bitgen))
                V = np.float32(f32_one + c * X)
                if (V > np.float32(0.0)):
                    break

            V = np.float32(V * V * V)
            U = next_float(bitgen)
            if (U < f32_one - np.float32(0.0331) * (X * X) * (X * X)):
                return np.float32(b * V)

            if (np.log(U) < np.float32(0.5) * X * X + b *
                    (f32_one - V + np.log(V))):
                return np.float32(b * V)


@register_jitable
def random_normal(bitgen, loc, scale):
    scaled_normal = scale * random_standard_normal(bitgen)
    return loc + scaled_normal


@register_jitable
def random_normal_f(bitgen, loc, scale):
    scaled_normal = np.float32(scale * random_standard_normal_f(bitgen))
    return np.float32(loc + scaled_normal)


@register_jitable
def random_exponential(bitgen, scale):
    return scale * random_standard_exponential(bitgen)


@register_jitable
def random_uniform(bitgen, lower, range):
    scaled_uniform = range * next_double(bitgen)
    return lower + scaled_uniform


@register_jitable
def random_gamma(bitgen, shape, scale):
    return scale * random_standard_gamma(bitgen, shape)


@register_jitable
def random_gamma_f(bitgen, shape, scale):
    return np.float32(scale * random_standard_gamma_f(bitgen, shape))


@register_jitable
def random_beta(bitgen, a, b):
    if a <= 1.0 and b <= 1.0:
        while 1:
            U = next_double(bitgen)
            V = next_double(bitgen)
            X = pow(U, 1.0 / a)
            Y = pow(V, 1.0 / b)
            XpY = X + Y
            if XpY <= 1.0 and XpY > 0.0:
                if (X + Y > 0):
                    return X / XpY
                else:
                    logX = np.log(U) / a
                    logY = np.log(V) / b
                    logM = min(logX, logY)
                    logX -= logM
                    logY -= logM

                    return np.exp(logX - np.log(np.exp(logX) + np.exp(logY)))
    else:
        Ga = random_standard_gamma(bitgen, a)
        Gb = random_standard_gamma(bitgen, b)
        return Ga / (Ga + Gb)


@register_jitable
def random_chisquare(bitgen, df):
    return 2.0 * random_standard_gamma(bitgen, df / 2.0)


@register_jitable
def random_f(bitgen, dfnum, dfden):
    return ((random_chisquare(bitgen, dfnum) * dfden) /
            (random_chisquare(bitgen, dfden) * dfnum))


@register_jitable
def random_standard_cauchy(bitgen):
    return random_standard_normal(bitgen) / random_standard_normal(bitgen)


@register_jitable
def random_pareto(bitgen, a):
    return np_expm1(random_standard_exponential(bitgen) / a)


@register_jitable
def random_weibull(bitgen, a):
    if (a == 0.0):
        return 0.0
    return pow(random_standard_exponential(bitgen), 1. / a)


@register_jitable
def random_power(bitgen, a):
    return pow(-np_expm1(-random_standard_exponential(bitgen)), 1. / a)


@register_jitable
def random_laplace(bitgen, loc, scale):
    U = next_double(bitgen)
    while U <= 0:
        U = next_double(bitgen)
    if (U >= 0.5):
        U = loc - scale * np.log(2.0 - U - U)
    elif (U > 0.0):
        U = loc + scale * np.log(U + U)
    return U


@register_jitable
def random_logistic(bitgen, loc, scale):
    U = next_double(bitgen)
    while U <= 0.0:
        U = next_double(bitgen)
    return loc + scale * np.log(U / (1.0 - U))


@register_jitable
def random_lognormal(bitgen, mean, sigma):
    return np.exp(random_normal(bitgen, mean, sigma))


@register_jitable
def random_standard_t(bitgen, df):
    num = random_standard_normal(bitgen)
    denom = random_standard_gamma(bitgen, df / 2)
    return np.sqrt(df / 2) * num / np.sqrt(denom)


@register_jitable
def random_wald(bitgen, mean, scale):
    mu_2l = mean / (2 * scale)
    Y = random_standard_normal(bitgen)
    Y = mean * Y * Y
    X = mean + mu_2l * (Y - np.sqrt(4 * scale * Y + Y * Y))
    U = next_double(bitgen)
    if (U <= mean / (mean + X)):
        return X
    else:
        return mean * mean / X


@register_jitable
def random_geometric_search(bitgen, p):
    X = 1
    sum = prod = p
    q = 1.0 - p
    U = next_double(bitgen)
    while (U > sum):
        prod *= q
        sum += prod
        X = X + 1
    return X


@register_jitable
def random_geometric_inversion(bitgen, p):
    return np.ceil(-random_standard_exponential(bitgen) / np.log1p(-p))


@register_jitable
def random_geometric(bitgen, p):
    if (p >= 0.333333333333333333333333):
        return random_geometric_search(bitgen, p)
    else:
        return random_geometric_inversion(bitgen, p)


@register_jitable
def random_zipf(bitgen, a):
    am1 = a - 1.0
    b = pow(2.0, am1)
    while 1:
        U = 1.0 - next_double(bitgen)
        V = next_double(bitgen)
        X = np.floor(pow(U, -1.0 / am1))
        if (X > INT64_MAX or X < 1.0):
            continue

        T = pow(1.0 + 1.0 / X, am1)
        if (V * X * (T - 1.0) / (b - 1.0) <= T / b):
            return X


@register_jitable
def random_triangular(bitgen, left, mode,
                      right):
    base = right - left
    leftbase = mode - left
    ratio = leftbase / base
    leftprod = leftbase * base
    rightprod = (right - mode) * base

    U = next_double(bitgen)
    if (U <= ratio):
        return left + np.sqrt(U * leftprod)
    else:
        return right - np.sqrt((1.0 - U) * rightprod)


@register_jitable
def random_loggam(x):
    a = [8.333333333333333e-02, -2.777777777777778e-03,
         7.936507936507937e-04, -5.952380952380952e-04,
         8.417508417508418e-04, -1.917526917526918e-03,
         6.410256410256410e-03, -2.955065359477124e-02,
         1.796443723688307e-01, -1.39243221690590e+00]

    if ((x == 1.0) or (x == 2.0)):
        return 0.0
    elif (x < 7.0):
        n = int(7 - x)
    else:
        n = 0

    x0 = x + n
    x2 = (1.0 / x0) * (1.0 / x0)
    # /* log(2 * M_PI) */
    lg2pi = 1.8378770664093453e+00
    gl0 = a[9]

    for k in range(0, 9):
        gl0 *= x2
        gl0 += a[8 - k]

    gl = gl0 / x0 + 0.5 * lg2pi + (x0 - 0.5) * np.log(x0) - x0
    if (x < 7.0):
        for k in range(1, n + 1):
            gl = gl - np.log(x0 - 1.0)
            x0 = x0 - 1.0

    return gl


@register_jitable
def random_poisson_mult(bitgen, lam):
    enlam = np.exp(-lam)
    X = 0
    prod = 1.0
    while (1):
        U = next_double(bitgen)
        prod *= U
        if (prod > enlam):
            X += 1
        else:
            return X


@register_jitable
def random_poisson_ptrs(bitgen, lam):

    slam = np.sqrt(lam)
    loglam = np.log(lam)
    b = 0.931 + 2.53 * slam
    a = -0.059 + 0.02483 * b
    invalpha = 1.1239 + 1.1328 / (b - 3.4)
    vr = 0.9277 - 3.6224 / (b - 2)

    while (1):
        U = next_double(bitgen) - 0.5
        V = next_double(bitgen)
        us = 0.5 - np.fabs(U)
        k = int((2 * a / us + b) * U + lam + 0.43)
        if ((us >= 0.07) and (V <= vr)):
            return k

        if ((k < 0) or ((us < 0.013) and (V > us))):
            continue

        # /* log(V) == log(0.0) ok here */
        # /* if U==0.0 so that us==0.0, log is ok since always returns */
        if ((np.log(V) + np.log(invalpha) - np.log(a / (us * us) + b)) <=
           (-lam + k * loglam - random_loggam(k + 1))):
            return k


@register_jitable
def random_poisson(bitgen, lam):
    if (lam >= 10):
        return random_poisson_ptrs(bitgen, lam)
    elif (lam == 0):
        return 0
    else:
        return random_poisson_mult(bitgen, lam)


@register_jitable
def random_negative_binomial(bitgen, n, p):
    Y = random_gamma(bitgen, n, (1 - p) / p)
    return random_poisson(bitgen, Y)


@register_jitable
def random_noncentral_chisquare(bitgen, df, nonc):
    if np.isnan(nonc):
        return np.nan

    if nonc == 0:
        return random_chisquare(bitgen, df)

    if 1 < df:
        Chi2 = random_chisquare(bitgen, df - 1)
        n = random_standard_normal(bitgen) + np.sqrt(nonc)
        return Chi2 + n * n
    else:
        i = random_poisson(bitgen, nonc / 2.0)
        return random_chisquare(bitgen, df + 2 * i)


@register_jitable
def random_noncentral_f(bitgen, dfnum, dfden, nonc):
    t = random_noncentral_chisquare(bitgen, dfnum, nonc) * dfden
    return t / (random_chisquare(bitgen, dfden) * dfnum)


@register_jitable
def random_logseries(bitgen, p):
    r = np_log1p(-p)

    while 1:
        V = next_double(bitgen)
        if (V >= p):
            return 1
        U = next_double(bitgen)
        q = -np.expm1(r * U)
        if (V <= q * q):
            result = np.int64(np.floor(1 + np.log(V) / np.log(q)))
            if result < 1 or V == 0.0:
                continue
            else:
                return result
        if (V >= q):
            return 1
        else:
            return 2


@register_jitable
def random_binomial_btpe(bitgen, n, p):
    r = min(p, 1.0 - p)
    q = 1.0 - r
    fm = n * r + r
    m = int(np.floor(fm))
    p1 = int(np.floor(2.195 * np.sqrt(n * r * q) - 4.6 * q) + 0.5)
    xm = m + 0.5
    xl = xm - p1
    xr = xm + p1
    c = 0.134 + 20.5 / (15.3 + m)
    a = (fm - xl) / (fm - xl * r)
    laml = a * (1.0 + a / 2.0)
    a = (xr - fm) / (xr * q)
    lamr = a * (1.0 + a / 2.0)
    p2 = p1 * (1.0 + 2.0 * c)
    p3 = p2 + c / laml
    p4 = p3 + c / lamr

    case = 10
    y = k = 0
    while 1:
        if case == 10:
            nrq = n * r * q
            u = next_double(bitgen) * p4
            v = next_double(bitgen)
            if (u > p1):
                case = 20
                continue
            y = int(np.floor(xm - p1 * v + u))
            case = 60
            continue
        elif case == 20:
            if (u > p2):
                case = 30
                continue
            x = xl + (u - p1) / c
            v = v * c + 1.0 - np.fabs(m - x + 0.5) / p1
            if (v > 1.0):
                case = 10
                continue
            y = int(np.floor(x))
            case = 50
            continue
        elif case == 30:
            if (u > p3):
                case = 40
                continue
            y = int(np.floor(xl + np.log(v) / laml))
            if ((y < 0) or (v == 0.0)):
                case = 10
                continue
            v = v * (u - p2) * laml
            case = 50
            continue
        elif case == 40:
            y = int(np.floor(xr - np.log(v) / lamr))
            if ((y > n) or (v == 0.0)):
                case = 10
                continue
            v = v * (u - p3) * lamr
            case = 50
            continue
        elif case == 50:
            k = abs(y - m)
            if ((k > 20) and (k < ((nrq) / 2.0 - 1))):
                case = 52
                continue
            s = r / q
            a = s * (n + 1)
            F = 1.0
            if (m < y):
                for i in range(m + 1, y + 1):
                    F = F * (a / i - s)
            elif (m > y):
                for i in range(y + 1, m + 1):
                    F = F / (a / i - s)
            if (v > F):
                case = 10
                continue
            case = 60
            continue
        elif case == 52:
            rho = (k / (nrq)) * \
                  ((k * (k / 3.0 + 0.625) + 0.16666666666666666) /
                   nrq + 0.5)
            t = -k * k / (2 * nrq)
            A = np.log(v)
            if (A < (t - rho)):
                case = 60
                continue
            if (A > (t + rho)):
                case = 10
                continue
            x1 = y + 1
            f1 = m + 1
            z = n + 1 - m
            w = n - y + 1
            x2 = x1 * x1
            f2 = f1 * f1
            z2 = z * z
            w2 = w * w
            if (A > (xm * np.log(f1 / x1) + (n - m + 0.5) * np.log(z / w) +
                     (y - m) * np.log(w * r / (x1 * q)) +
                     (13680. - (462. - (132. - (99. - 140. / f2) / f2) / f2)
                      / f2) / f1 / 166320. +
                     (13680. - (462. - (132. - (99. - 140. / z2) / z2) / z2)
                      / z2) / z / 166320. +
                     (13680. - (462. - (132. - (99. - 140. / x2) / x2) / x2)
                      / x2) / x1 / 166320. +
                     (13680. - (462. - (132. - (99. - 140. / w2) / w2) / w2)
                      / w2) / w / 66320.)):
                case = 10
                continue
        elif case == 60:
            if (p > 0.5):
                y = n - y
            return y


@register_jitable
def random_binomial_inversion(bitgen, n, p):
    q = 1.0 - p
    qn = np.exp(n * np.log(q))
    _np = n * p
    bound = min(n, _np + 10.0 * np.sqrt(_np * q + 1))

    X = 0
    px = qn
    U = next_double(bitgen)
    while (U > px):
        X = X + 1
        if (X > bound):
            X = 0
            px = qn
            U = next_double(bitgen)
        else:
            U -= px
            px = ((n - X + 1) * p * px) / (X * q)

    return X


@register_jitable
def random_binomial(bitgen, n, p):
    if ((n == 0) or (p == 0.0)):
        return 0

    if (p <= 0.5):
        if (p * n <= 30.0):
            return random_binomial_inversion(bitgen, n, p)
        else:
            return random_binomial_btpe(bitgen, n, p)
    else:
        q = 1.0 - p
        if (q * n <= 30.0):
            return n - random_binomial_inversion(bitgen, n, q)
        else:
            return n - random_binomial_btpe(bitgen, n, q)
