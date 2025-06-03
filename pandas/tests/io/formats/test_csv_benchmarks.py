import pytest
import numpy as np
import pandas as pd

pytestmark = pytest.mark.usefixtures("benchmark")

def test_benchmark_old_style_format(benchmark):
    df = pd.DataFrame(np.random.rand(1000, 1000))
    benchmark(lambda: df.to_csv(float_format="%.6f"))

def test_benchmark_new_style_format(benchmark):
    df = pd.DataFrame(np.random.rand(1000, 1000))
    benchmark(lambda: df.to_csv(float_format="{:.6f}"))

def test_benchmark_new_style_thousands(benchmark):
    df = pd.DataFrame(np.random.rand(1000, 1000))
    benchmark(lambda: df.to_csv(float_format="{:,.2f}"))

def test_benchmark_callable_format(benchmark):
    df = pd.DataFrame(np.random.rand(1000, 1000))
    benchmark(lambda: df.to_csv(float_format=lambda x: f"{x:.6f}"))