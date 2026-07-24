# This file is part of Hypothesis, which may be found at
# https://github.com/HypothesisWorks/hypothesis/
#
# Copyright the Hypothesis Authors.
# Individual contributors are listed in AUTHORS.rst and the git log.
#
# This Source Code Form is subject to the terms of the Mozilla Public License,
# v. 2.0. If a copy of the MPL was not distributed with this file, You can
# obtain one at https://mozilla.org/MPL/2.0/.

import abc
import math

from hypothesis.internal.floats import clamp

# Liam: please be aware that an LLM wrote the stdtr and stdtrit functions, with some
# guidance from me on comparing to scipy's existing implementation.
#
# This code is strongly in-distribution for AI, and is tested against a strong
# oracle (scipy's implementation), which are the only reasons I am willing to
# accept this code without deeply understanding it.


def stdtr(df: int, t: float) -> float:  # pragma: no cover  # covered by tests/scipy
    """Student's t CDF for integer df >= 1, evaluated at t.

    Closed-form finite sum from Abramowitz & Stegun 26.7.7-8.
    """
    assert isinstance(df, int)
    if df < 1:
        raise ValueError(f"stdtr requires integer df >= 1, got {df}")
    if t == 0.0:
        return 0.5
    abs_t = abs(t)
    z = 1.0 + abs_t * abs_t / df
    if df % 2 == 1:
        # odd: includes an arctan term
        u = abs_t / math.sqrt(df)
        p = math.atan(u)
        if df > 1:
            f = 1.0
            tz = 1.0
            j = 3
            while j <= df - 2:
                tz *= (j - 1) / (z * j)
                f += tz
                j += 2
            p += f * u / z
        p *= 2.0 / math.pi
    else:
        # even: simple finite sum, no arctan
        f = 1.0
        tz = 1.0
        j = 2
        while j <= df - 2:
            tz *= (j - 1) / (z * j)
            f += tz
            j += 2
        p = f * abs_t / math.sqrt(z * df)
    return 0.5 - 0.5 * p if t < 0 else 0.5 + 0.5 * p


def stdtrit(
    df: int, p: float, *, eps: float = 1e-10, max_iter: int = 50
) -> float:  # pragma: no cover  # covered by tests/scipy
    """Inverse Student's t CDF (quantile) for integer df >= 1.

    df ∈ {1, 2}: closed-form analytic quantile (Shaw 2006, eq 35-36).

    df >= 3: bracketed Newton iteration on stdtr — doubling search from t=1
    to bracket the target, then Newton steps with bisection fallback. Always
    converges; ~5-10 iterations typical.
    """
    if df < 1:
        raise ValueError(f"stdtrit requires integer df >= 1, got {df}")
    if p == 0.5:
        return 0.0
    if not 0.0 < p < 1.0:
        raise ValueError(f"stdtrit requires 0 < p < 1, got {p}")
    if df == 1:
        # Cauchy: F^{-1}(p) = -cot(pi p). Reflect via 1-p when p > 0.5 because
        # sin(pi p) near pi suffers cancellation; sin(pi (1-p)) near 0 is exact.
        if p > 0.5:
            return math.cos(math.pi * (1 - p)) / math.sin(math.pi * (1 - p))
        return -math.cos(math.pi * p) / math.sin(math.pi * p)
    if df == 2:
        return (2 * p - 1) / math.sqrt(2 * p * (1 - p))

    sign = 1.0 if p > 0.5 else -1.0
    q = p if p > 0.5 else 1 - p

    lo, hi = 0.0, 1.0
    while stdtr(df, hi) < q:
        hi *= 2

    log_norm = (
        math.lgamma(0.5 * (df + 1))
        - 0.5 * math.log(df * math.pi)
        - math.lgamma(0.5 * df)
    )
    t = 0.5 * (lo + hi)
    for _ in range(max_iter):
        F = stdtr(df, t)
        if F < q:
            lo = t
        else:
            hi = t
        log_f = log_norm - 0.5 * (df + 1) * math.log1p(t * t / df)
        f = math.exp(log_f)
        if f == 0.0:
            t = 0.5 * (lo + hi)
        else:
            t_newton = t - (F - q) / f
            t = t_newton if lo <= t_newton <= hi else 0.5 * (lo + hi)
        if hi - lo < eps * (1 + abs(t)):
            break
    return sign * t


# Liam: I'm not convinced this abstraction pays for itself. If it's a hindrance in the
# future, feel free to refactor. This class is mainly here to concretize the minimal set
# of abstractions our integer sampling code expects a distribution to define, and to
# provide stronger guardrails around composition.
#
# Note that our code makes the explicit assumption that any _Distribution
# implementation is symmetric around 0. This implies inverse_cdf(0.5) = 0, for example.
# Do not break this assumption without verifying carefully against callers.
class _Distribution(abc.ABC):
    @abc.abstractmethod
    def cdf(self, x: float) -> float:
        raise NotImplementedError

    @abc.abstractmethod
    def inverse_cdf(self, u: float) -> float:
        raise NotImplementedError

    @abc.abstractmethod
    def pdf(self, x: float) -> float:
        raise NotImplementedError


class UniformDistribution(_Distribution):
    """Uniform distribution on [-half_width, half_width]."""

    def __init__(self, *, half_width: float) -> None:
        self.half_width = half_width

    def cdf(self, x: float) -> float:
        if x < -self.half_width:
            return 0.0  # pragma: no cover
        if x > self.half_width:
            return 1.0  # pragma: no cover
        return (x + self.half_width) / (2 * self.half_width)

    def inverse_cdf(self, u: float) -> float:
        return -self.half_width + 2 * self.half_width * u

    def pdf(self, x: float) -> float:
        if -self.half_width <= x <= self.half_width:
            return 1 / (2 * self.half_width)
        return 0.0  # pragma: no cover


class LogStudentTDistribution(_Distribution):
    """Student's t distribution, in the transformed domain of log_2(x).

    Y = sign(x) * log_2(1 + |x|) ~ scale_bits * t(df).

    Note that we only support integer df. This is to reduce surface area in the vendored
    mathematical functions, which are simpler when df is an integer. There is no
    fundamental difficulty to supporting float df, other than the necessary due diligence
    for the vendored code.
    """

    LN2 = math.log(2)

    def __init__(self, *, scale_bits: float, df: int) -> None:
        self.scale_bits = scale_bits
        self.df = df
        self._t_coef = math.gamma((df + 1) / 2) / (
            math.sqrt(df * math.pi) * math.gamma(df / 2)
        )

    def cdf(self, x: float) -> float:
        y = math.copysign(math.log2(1 + abs(x)), x) / self.scale_bits
        return float(stdtr(self.df, y))

    def inverse_cdf(self, u: float) -> float:
        y = self.scale_bits * float(stdtrit(self.df, u))
        # 2^1023 is the largest power of 2 below sys.float_info.max. This makes y=1023
        # the largest value we can turn into a float to return here without overflowing.
        #
        # In the rare case that we draw such an extreme value, clamp to the bounds.
        y = clamp(-1023, y, 1023)
        return math.copysign(math.expm1(abs(y) * self.LN2), y)

    def pdf(self, x: float) -> float:
        y = math.copysign(math.log2(1 + abs(x)), x) / self.scale_bits
        f_t = self._t_coef * (1 + y * y / self.df) ** (-(self.df + 1) / 2)
        return f_t / (self.scale_bits * (1 + abs(x)) * self.LN2)


class PiecewiseDistribution:
    """Two-region splice: `inner` on (-switchover, switchover),
    `outer` on |x| >= switchover.

    Each region is normalized so that the resulting density is continuous at
    ±switchover and integrates to 1. Both inner and outer must be symmetric around 0.
    """

    def __init__(
        self, *, inner: _Distribution, outer: _Distribution, switchover: float
    ) -> None:
        # Note that this code could be made simpler by taking advantage of our assumption that
        # `inner` and `outer` are symmetric around 0. The code is currently generic enough
        # to support asymmetric distributions, but doesn't need to be,
        self.inner = inner
        self.outer = outer
        self.switchover = switchover
        self._inner_g_neg = inner.cdf(-switchover)
        self._inner_g_pos = inner.cdf(switchover)
        self._outer_g_neg = outer.cdf(-switchover)
        self._outer_g_pos = outer.cdf(switchover)
        outer_outer_mass = 1 - (self._outer_g_pos - self._outer_g_neg)
        inner_inner_mass = self._inner_g_pos - self._inner_g_neg
        # density continuity: beta * inner.pdf(c) = alpha * outer.pdf(c)
        # plus total mass 1; solve for (alpha, beta).
        inner_pdf = inner.pdf(switchover)
        outer_pdf = outer.pdf(switchover)
        assert inner_pdf != 0, (
            f"inner.pdf(switchover={switchover}) == 0; cannot density-match. inner "
            "must have support at the boundary."
        )
        self._alpha = 1 / (outer_pdf * inner_inner_mass / inner_pdf + outer_outer_mass)
        self._beta = self._alpha * outer_pdf / inner_pdf
        self._inner_mass = self._beta * inner_inner_mass
        self._left_mass = self._alpha * self._outer_g_neg

    def cdf(self, x: float) -> float:
        if x <= -self.switchover:
            return self._alpha * self.outer.cdf(x)
        if x < self.switchover:
            return self._left_mass + self._beta * (
                self.inner.cdf(x) - self._inner_g_neg
            )
        return (
            self._left_mass
            + self._inner_mass
            + self._alpha * (self.outer.cdf(x) - self._outer_g_pos)
        )

    def inverse_cdf(self, u: float) -> float:
        if u <= self._left_mass:
            return self.outer.inverse_cdf(u / self._alpha)
        if u < self._left_mass + self._inner_mass:
            target = self._inner_g_neg + (u - self._left_mass) / self._beta
            return self.inner.inverse_cdf(target)
        return self.outer.inverse_cdf(
            (u - self._left_mass - self._inner_mass) / self._alpha + self._outer_g_pos
        )
