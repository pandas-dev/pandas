import math

import numpy as np

import unittest
from numba import njit
from numba.extending import register_jitable
from numba.tests.support import TestCase


RISKFREE = 0.02
VOLATILITY = 0.30


A1 = 0.31938153
A2 = -0.356563782
A3 = 1.781477937
A4 = -1.821255978
A5 = 1.330274429
RSQRT2PI = 0.39894228040143267793994605993438


@register_jitable
def cnd_array(d):
    K = 1.0 / (1.0 + 0.2316419 * np.abs(d))
    ret_val = (RSQRT2PI * np.exp(-0.5 * d * d) *
               (K * (A1 + K * (A2 + K * (A3 + K * (A4 + K * A5))))))
    return np.where(d > 0, 1.0 - ret_val, ret_val)


@register_jitable
def cnd(d):
    K = 1.0 / (1.0 + 0.2316419 * math.fabs(d))
    ret_val = (RSQRT2PI * math.exp(-0.5 * d * d) *
               (K * (A1 + K * (A2 + K * (A3 + K * (A4 + K * A5))))))
    if d > 0:
        ret_val = 1.0 - ret_val
    return ret_val


@njit
def blackscholes_arrayexpr(stockPrice, optionStrike, optionYears, Riskfree,
                           Volatility):
    S = stockPrice
    X = optionStrike
    T = optionYears
    R = Riskfree
    V = Volatility
    sqrtT = np.sqrt(T)
    d1 = (np.log(S / X) + (R + 0.5 * V * V) * T) / (V * sqrtT)
    d2 = d1 - V * sqrtT
    cndd1 = cnd_array(d1)
    cndd2 = cnd_array(d2)

    expRT = np.exp(- R * T)

    callResult = (S * cndd1 - X * expRT * cndd2)
    putResult = (X * expRT * (1.0 - cndd2) - S * (1.0 - cndd1))
    return callResult, putResult


@njit
def blackscholes_scalar(callResult, putResult, stockPrice, optionStrike,
                        optionYears, Riskfree, Volatility):
    S = stockPrice
    X = optionStrike
    T = optionYears
    R = Riskfree
    V = Volatility
    for i in range(len(S)):
        sqrtT = math.sqrt(T[i])
        d1 = (math.log(S[i] / X[i]) + (R + 0.5 * V * V) * T[i]) / (V * sqrtT)
        d2 = d1 - V * sqrtT
        cndd1 = cnd(d1)
        cndd2 = cnd(d2)

        expRT = math.exp((-1. * R) * T[i])
        callResult[i] = (S[i] * cndd1 - X[i] * expRT * cndd2)
        putResult[i] = (X[i] * expRT * (1.0 - cndd2) - S[i] * (1.0 - cndd1))


def randfloat(rand_var, low, high):
    return (1.0 - rand_var) * low + rand_var * high


class TestBlackScholes(TestCase):
    def test_array_expr(self):
        OPT_N = 400

        stockPrice = randfloat(self.random.random_sample(OPT_N), 5.0, 30.0)
        optionStrike = randfloat(self.random.random_sample(OPT_N), 1.0, 100.0)
        optionYears = randfloat(self.random.random_sample(OPT_N), 0.25, 10.0)

        args = stockPrice, optionStrike, optionYears, RISKFREE, VOLATILITY

        callResultGold, putResultGold = blackscholes_arrayexpr.py_func(*args)
        callResultNumba, putResultNumba = blackscholes_arrayexpr(*args)

        delta = np.abs(callResultGold - callResultNumba)
        self.assertAlmostEqual(delta.max(), 0)

    def test_scalar(self):
        OPT_N = 400

        callResultGold = np.zeros(OPT_N)
        putResultGold = np.zeros(OPT_N)

        callResultNumba = np.zeros(OPT_N)
        putResultNumba = np.zeros(OPT_N)

        stockPrice = randfloat(self.random.random_sample(OPT_N), 5.0, 30.0)
        optionStrike = randfloat(self.random.random_sample(OPT_N), 1.0, 100.0)
        optionYears = randfloat(self.random.random_sample(OPT_N), 0.25, 10.0)

        args = stockPrice, optionStrike, optionYears, RISKFREE, VOLATILITY

        blackscholes_scalar.py_func(callResultGold, putResultGold, *args)
        blackscholes_scalar(callResultNumba, putResultNumba, *args)

        delta = np.abs(callResultGold - callResultNumba)
        self.assertAlmostEqual(delta.max(), 0)


if __name__ == "__main__":
    unittest.main()
