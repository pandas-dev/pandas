"""
Tests for DataFrame reductions that are DataFrame-specific, i.e. cannot
be shared in tests.reductions.
"""
import pandas as pd
import numpy as np


if True:#def test_blockwise_reduction():
	arr = np.arange(10)
	tdarr = arr.astype("m8[s]")
	df = pd.DataFrame({"A": arr, "B": tdarr})
