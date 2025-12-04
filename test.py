# mypy: ignore-errors
import numpy as np

from backtrace import trace_calls
import pandas as pd


def multiply_elements(a, b):
    return a * b


@trace_calls()
def safe_multiply_series(
    s1: pd.Series, s2: pd.Series, threshold: float = 1e15
) -> pd.Series:
    result = []
    for a, b in zip(s1, s2, strict=False):
        product = multiply_elements(a, b)
        if abs(product) > threshold:
            res = float(product)
        else:
            res = product
        result.append(res)
    return pd.Series(result)


def test_precision_difference():
    val1 = 116.8000030517578151564387
    val2 = 342.4357239847609283475098

    # Float64 (pandas default)
    s1 = pd.Series([val1], dtype=np.float64)
    s2 = pd.Series([val2], dtype=np.float64)
    result_f64 = ((s1 * s2) + s1).values[0]

    # Float128 computation
    val1_f128 = np.float128(val1)
    val2_f128 = np.float128(val2)
    expected_f128 = (val1_f128 * val2_f128) + val1_f128

    print(f"PD Float64 result:  {result_f64:.20f}")
    print(f"NP Float128 result: {expected_f128:.20f}")
    print(f"Difference:      {abs((expected_f128) - result_f64):.32e}")

    # # Float64 computation for comparison
    # expected_f64 = np.float64(val1) * np.float64(val2) + np.float64(val1)
    # print(f"NP Float64 result:  {expected_f64:.20f}")
    # print(f"\nTypes: {type(result_f64)}, {type(expected_f128)}")


if __name__ == "__main__":
    test_precision_difference()
