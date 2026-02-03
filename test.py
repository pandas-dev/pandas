# mypy: ignore-errors
import timeit

import numpy as np

from pandas._libs.double_double import DoubleDouble


def test_basic_ops():
    test_cases = {
        "large_large": {"a": 10, "b": 15},
        "large_small": {"a": 1.0, "b": 1e-16},
    }

    ops = {
        "add": lambda a, b: a + b,
        "sub": lambda a, b: a - b,
        "mul": lambda a, b: a * b,
    }

    def run_float64(op, a, b):
        return op(a, b)

    def run_dd(op, a, b):
        return float(op(DoubleDouble(a), DoubleDouble(b)))

    for case_name, data in test_cases.items():
        a, b = data["a"], data["b"]
        print(f"\n  ======= {case_name} (a={a}, b={b}) =======")
        print(
            f"    {'Op':<6} | {'Float64 (µs)':<14} | "
            f"{'DoubleDouble (µs)':<18} | {'DD/F64':>6}"
        )
        print("  " + "-" * 55)

        for op_name, op in ops.items():
            f64_val = run_float64(op, a, b)
            dd_val = run_dd(op, a, b)
            # Ensures DD value matches F64
            assert f64_val == dd_val, f"{op_name}: F64={f64_val} != DD={dd_val}"

            f64_time = timeit.timeit(lambda: run_float64(op, a, b), number=1000) / 1000
            dd_time = timeit.timeit(lambda: run_dd(op, a, b), number=1000) / 1000

            print(
                f"    {op_name:<6} | {f64_time * 1e6:<14.3f} | "
                f"{dd_time * 1e6:<18.3f} | {dd_time / f64_time:>5.1f}x"
            )


def test_precision():
    large = 1.0
    small = 1e-15
    n = 1000
    expected = np.float128(large) + n * np.float128(small)

    def float64_sum():
        result = large
        for _ in range(n):
            result += small
        return result

    def dd_sum():
        result = DoubleDouble(large)
        small_dd = DoubleDouble(small)
        for _ in range(n):
            result += small_dd
        return float(result)

    methods = {"Float64": float64_sum, "DoubleDouble": dd_sum}

    print(f"\n  ======= Precision Test (large={large}, small={small}, n={n}) =======")
    print(f"  Expected: {expected:.20f}")
    print(f"{'  Method':<12} | {'Result':<24} | {'Error':<12} | {'Time (µs)':<10}")
    print("  " + "-" * 65)

    for name, func in methods.items():
        result = func()
        elapsed = timeit.timeit(func, number=1000) / 1000
        error = abs(result - expected)
        print(
            f"  {name:<12} | {result:<24.20f} |{error:<12.6e} | {elapsed * 1e6:<10.3f}"
        )


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

    print(f"  {'Method':<12} | {'Result':<12} |{'Time (µs)':<10}")
    print("  " + "-" * 40)
    print(f"  {'Float64':<12} | {result_f64:<12f} | {time_f64 * 1e6:<10.3f}")
    print(
        f"  {'DoubleDouble':<12} | {result_dd.collapse():<12.0f} | "
        f"{time_dd * 1e6:<10.3f}"
    )


def test_welford():
    """Test Welford's algorithm for online variance and correlation."""
    test_cases = {
        "Original Issue": {
            "x": [30, 30.100000381469727],
            "y": [116.80000305175781, 116.8000030517578],
        },
        "Extended": {
            "x": [30, 30.100000381469727, 30.300000381469727],
            "y": [116.80000305175781, 116.8000030517578, 116.80000305175783],
        },
    }

    def welford(xs, ys, use_dd=False):
        wrap = DoubleDouble if use_dd else lambda x: x
        s = {k: wrap(0.0) for k in ["mean_x", "mean_y", "m2_x", "m2_y", "co_moment"]}
        count = 0
        for x, y in zip(xs, ys, strict=False):
            count += 1
            x_val, y_val = wrap(x), wrap(y)
            inv = wrap(1.0 / count) if use_dd else 1.0 / count

            dx = x_val - s["mean_x"]
            s["mean_x"] = s["mean_x"] + dx * inv
            dx2 = x_val - s["mean_x"]
            s["m2_x"] = s["m2_x"] + dx * dx2

            dy = y_val - s["mean_y"]
            s["mean_y"] = s["mean_y"] + dy * inv
            dy2 = y_val - s["mean_y"]
            s["m2_y"] = s["m2_y"] + dy * dy2

            s["co_moment"] = s["co_moment"] + dx * dy2

        f = float if use_dd else lambda x: x
        var_x, var_y = f(s["m2_x"]) / count, f(s["m2_y"]) / count
        denom = f(s["m2_x"]) * f(s["m2_y"])
        corr = f(s["co_moment"]) / denom**0.5 if denom > 0 else 0.0
        return {"corr": corr, "var_x": var_x, "var_y": var_y}

    methods = {
        "Float64": lambda xs, ys: welford(xs, ys, use_dd=False),
        "DoubleDouble": lambda xs, ys: welford(xs, ys, use_dd=True),
    }

    for case_name, data in test_cases.items():
        data_x, data_y = data["x"], data["y"]
        x_f128 = np.array(data_x, dtype=np.float128)
        y_f128 = np.array(data_y, dtype=np.float128)
        true = {
            "var_x": np.var(x_f128),
            "var_y": np.var(y_f128),
            "corr": np.corrcoef(x_f128, y_f128)[0, 1],
        }

        print(f"\n  === {case_name} ===")
        for k, v in true.items():
            print(f"  True {k}: {v:.20f}")

        print(
            f"\n  {'Method':<12} | {'Metric':<8} | {'Value':<24} | "
            f"{'Error %':<12} | {'Time (µs)':<10}"
        )

        for name, func in methods.items():
            results = func(data_x, data_y)
            elapsed = timeit.timeit(lambda f=func: f(data_x, data_y), number=100) / 100
            for i, (k, v) in enumerate(results.items()):
                pct_err = (
                    (np.float128(v) - true[k]) * 100 / true[k] if true[k] != 0 else 0
                )
                if i == 0:
                    print("  " + "-" * 77)
                    print(
                        f"  {name:<12} | {k:<8} | {v:<24.20f} | "
                        f"{pct_err:<12.6f} | {elapsed * 1e6:<10.3f}"
                    )
                else:
                    print(f"  {'':<12} | {k:<8} | {v:<24.20f} | {pct_err:<12.6f} |")


if __name__ == "__main__":
    # Put print statements into file -> 1
    STDOUT2FILE = 0

    tests = {
        "Basic Ops (exclu div)": test_basic_ops,
        "Simple Precision Calculation": test_precision,
        "Catastriphic Collision": test_catastrophic_collision,
        "Welford Calculation": test_welford,
    }

    if STDOUT2FILE:
        import sys

        orig_stdout = sys.stdout
        f = open("test_out.md", "w+")
        sys.stdout = f

    for name, func in tests.items():
        print("=" * 60)
        print(f"{name}:")
        print("=" * 60)
        func()

        print()

    if STDOUT2FILE:
        sys.stdout = orig_stdout
        f.close()
        print('Output written to "test_out.md"')
