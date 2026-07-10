import pytest
import numpy as np
import pandas as pd
import math
from pandas.errors import NumbaUtilError

# NumPy 2.5 compatibility patch for Numba in test runtime
np.row_stack = np.vstack

pytest.importorskip("numba")

pytestmark = [
    pytest.mark.filterwarnings("ignore::DeprecationWarning"),
    pytest.mark.filterwarnings("ignore::UserWarning"),
]


def test_numba_apply_opt_in():
    # Explicitly opting in to JIT should work regardless of Series size
    s = pd.Series(np.random.randn(100))
    
    def udf(x):
        return math.sin(x) * 2.0
        
    res_jit = s.apply(udf, engine="numba")
    res_py = s.apply(udf, engine="python")
    
    pd.testing.assert_series_equal(res_jit, res_py)


def test_numba_apply_global_option():
    # Test transparent auto-JIT using the global configuration option
    s = pd.Series(np.random.randn(60_000))
    
    def udf(x):
        return x * 2.0 + 1.0

    # Turn global option ON
    with pd.option_context("compute.use_numba", True):
        res_auto_jit = s.apply(udf)
        
    # Turn global option OFF
    with pd.option_context("compute.use_numba", False):
        res_standard = s.apply(udf)
        
    pd.testing.assert_series_equal(res_auto_jit, res_standard)


def test_numba_apply_unsupported_ops_raises_error():
    # If the user explicitly opts in, compilation errors must raise rather than falling back silently
    s = pd.Series(np.random.randn(100))
    
    def udf_invalid(x):
        # eval is unsupported by Numba's nopython compiler
        eval("x + 1")
        return x

    with pytest.raises(NumbaUtilError, match="Failed to execute Numba JIT loop"):
        s.apply(udf_invalid, engine="numba")


def test_numba_apply_non_numeric_dtype_raises_error():
    # Explicit JIT on non-numeric series should raise ValueError
    s = pd.Series(["a", "b", "c"] * 100, dtype=object)
    
    def udf_str(x):
        return x + "_suffix"

    with pytest.raises(ValueError, match="Numba engine only supports numeric/datetime dtypes"):
        s.apply(udf_str, engine="numba")
