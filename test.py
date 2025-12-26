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
        result = DoubleDouble(large)
        small_dd = DoubleDouble(small)
        for _ in range(n):
            result = result + small_dd
        return float(result)

    methods = {"Float64": float64_sum, "DoubleDouble": dd_sum}

    print(f"Expected: {expected:.20f}")

    for name, func in methods.items():
        result = func()
        elapsed = timeit.timeit(func, number=1000) / 1000
        error = abs(result - expected)
        print(f"{name}: {result:.20f}")
        print(f"    Error: {error:.6e} | Time: {elapsed * 1e6:.3f} µs")


def test_catastrophic_collision():
    a = 1e17
    b = 1

    def f64_op():
        return (a + b) - a

    def dd_op():
        dd_a = DoubleDouble(a)
        dd_b = DoubleDouble(b)
        return (dd_a + dd_b) - dd_a

    n = 1000
    time_f64 = timeit.timeit(f64_op, number=n) / n
    time_dd = timeit.timeit(dd_op, number=n) / n

    result_f64 = f64_op()
    result_dd = dd_op()
    DoubleDouble(a)
    print("Catastrophic Collision ((1e17 + 1) - 1e17)")
    print(f"Float64: {result_f64} ({time_f64 * 1e6:.3f} µs)")
    print(f"DoubleDouble: {result_dd.collapse()} ({time_dd * 1e6:.3f} µs)")


if __name__ == "__main__":
    tests = {
        "Sum Precision": test_summation_precision,
        "Catastriphic Collision": test_catastrophic_collision,
    }

    for name, func in tests.items():
        print("=" * 40)
        print(f"{name}:")
        print("=" * 40)
        func()

        print()
