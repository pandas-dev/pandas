from __future__ import annotations

import numpy as np
import pytest

pytestmark = pytest.mark.gpu

import dask.array as da
from dask.array.gufunc import apply_gufunc
from dask.array.utils import assert_eq

cupy = pytest.importorskip("cupy")


def test_apply_gufunc_axis():
    def mydiff(x):
        return np.diff(x)

    a = cupy.random.default_rng().standard_normal((3, 6, 4))
    da_ = da.from_array(a, chunks=2, asarray=False)

    m = np.diff(a, axis=1)
    dm = apply_gufunc(
        mydiff, "(i)->(i)", da_, axis=1, output_sizes={"i": 5}, allow_rechunk=True
    )
    assert_eq(m, dm)
