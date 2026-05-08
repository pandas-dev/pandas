"""
Benchmarks for _Unstacker on complete Cartesian-product MultiIndexes.

Run with:
    asv run --bench unstack_cartesian
"""

import numpy as np

import pandas as pd


class UnstackCartesian:
    """Series.unstack on a fully-populated Cartesian-product MultiIndex."""

    params = [
        [(5, 4), (100, 90), (200, 92), (10, 10_000), (10, 10, 10), (30, 30, 30)],
    ]
    param_names = ["shape"]

    def setup(self, shape):
        if len(shape) == 2:
            n_outer, n_inner = shape
            x0 = np.tile(np.arange(n_inner, dtype=float), n_outer)
            x1 = np.repeat(np.arange(n_outer, dtype=float), n_inner)
            idx = pd.MultiIndex.from_arrays([x0, x1], names=["x0", "x1"])
            self._last_level = "x0"
        else:
            idx = pd.MultiIndex.from_product(
                [np.arange(n) for n in shape],
                names=[f"x{i}" for i in range(len(shape))],
            )
            self._last_level = f"x{len(shape) - 1}"
        self.s = pd.Series(np.random.rand(len(idx)), index=idx)

    def time_unstack(self, shape):
        self.s.unstack(self._last_level)
