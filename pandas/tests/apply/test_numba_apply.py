import pytest
import numpy as np
import pandas as pd
import math

# NumPy 2.5 compatibility patch for Numba in test runtime
np.row_stack = np.vstack

pytestmark = [
    pytest.mark.filterwarnings("ignore::DeprecationWarning"),
    pytest.mark.filterwarnings("ignore::UserWarning"),
]


def test_numba_apply_correctness():
    # Large numeric series (> 50k elements) to trigger JIT
    s = pd.Series(np.random.randn(60_000))
    
    def udf(x):
        if x > 0:
            return math.sin(x) * 2.0
        else:
            return math.cos(x) - 1.0

    res_jit = s.apply(udf)
    
    # Force interpreter by temporarily mocking maybe_run_numba_apply
    import pandas.core.util.numba_ as nb_util
    orig = nb_util.maybe_run_numba_apply
    try:
        nb_util.maybe_run_numba_apply = lambda series, func: None
        res_interpreter = s.apply(udf)
    finally:
        nb_util.maybe_run_numba_apply = orig
    
    pd.testing.assert_series_equal(res_jit, res_interpreter)


def test_numba_apply_small_series():
    # Small numeric series (< 50k elements) should bypass JIT
    s = pd.Series(np.random.randn(100))
    
    def udf(x):
        return x * 2.0 + 1.0

    # Ensure it doesn't fail and returns correct results
    res = s.apply(udf)
    expected = s * 2.0 + 1.0
    pd.testing.assert_series_equal(res, expected)


def test_numba_apply_unsupported_ops_fallback():
    # Large numeric series to trigger JIT check
    s = pd.Series(np.random.randn(60_000))
    
    # A function that uses list operations unsupported by Numba in nopython mode
    # will fail compilation and should fallback gracefully
    def udf_fallback(x):
        l = []
        l.append(x)
        return len(l) * 2.5 + x

    res = s.apply(udf_fallback)
    expected = s + 2.5
    pd.testing.assert_series_equal(res, expected)


def test_numba_apply_non_numeric_dtype():
    # Object series should bypass JIT check
    s = pd.Series(["a", "b", "c"] * 20_000, dtype=object)
    
    def udf_str(x):
        return x + "_suffix"

    res = s.apply(udf_str)
    expected = s + "_suffix"
    pd.testing.assert_series_equal(res, expected, check_dtype=False)
