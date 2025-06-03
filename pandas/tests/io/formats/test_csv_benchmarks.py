import numpy as np
import pytest

from pandas import DataFrame

pytestmark = pytest.mark.usefixtures("benchmark")

# Create a single generator instance for all tests
rng = np.random.default_rng(seed=42)


def test_benchmark_old_style_format(benchmark):
    df = DataFrame(rng.random((1000, 1000)))
    benchmark(lambda: df.to_csv(float_format="%.6f"))


def test_benchmark_new_style_format(benchmark):
    df = DataFrame(rng.random((1000, 1000)))
    benchmark(lambda: df.to_csv(float_format="{:.6f}"))


def test_benchmark_new_style_thousands(benchmark):
    df = DataFrame(rng.random((1000, 1000)))
    benchmark(lambda: df.to_csv(float_format="{:,.2f}"))


def test_benchmark_callable_format(benchmark):
    df = DataFrame(rng.random((1000, 1000)))
    benchmark(lambda: df.to_csv(float_format=lambda x: f"{x:.6f}"))
