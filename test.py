# mypy: ignore-errors
import timeit

from pandas._libs.double_double import DoubleDouble


def test_summation_precision():
    large = 1.0
    small = 1e-16
    n = 1000
    expected = large + n * small

    def float64_sum():
        result = large
        for _ in range(n):
            result += small
        return result

    def dd_sum():
        result = DoubleDouble(large, 0.0)
        small_dd = DoubleDouble(small, 0.0)
        for _ in range(n):
            result = result + small_dd
        return float(result)

    methods = {"Float64": float64_sum, "DoubleDouble": dd_sum}

    print(f"Summation test (1.0 + {n} * 1e-16):")
    print(f"  Expected: {expected:.20f}\n")

    for name, func in methods.items():
        result = func()
        elapsed = timeit.timeit(func, number=1000) / 1000
        error = abs(result - expected)
        print(f"  {name:12s}: {result:.20f}")
        print(f"    Error: {error:.6e} | Time: {elapsed:.6f}s\n")


if __name__ == "__main__":
    test_summation_precision()
