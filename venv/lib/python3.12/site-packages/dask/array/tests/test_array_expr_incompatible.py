from __future__ import annotations

from unittest import mock

import numpy as np
import pytest

import dask.array as da
from dask.array.core import concatenate3
from dask.array.tests.test_dispatch import EncapsulateNDArray

if da._array_expr_enabled():
    pytest.skip("parametrize using unsupported functions", allow_module_level=True)


@pytest.mark.parametrize("one_d", [True, False])
@mock.patch.object(da.core, "_concatenate2", wraps=da.core._concatenate2)
def test_concatenate3_nep18_dispatching(mock_concatenate2, one_d):
    x = EncapsulateNDArray(np.arange(10))
    concat = [x, x] if one_d else [[x[None]], [x[None]]]
    result = concatenate3(concat)
    assert type(result) is type(x)
    mock_concatenate2.assert_called()
    mock_concatenate2.reset_mock()

    # When all the inputs are supported by plain `np.concatenate`, we should take the concatenate3
    # fastpath of allocating the full array up front and writing blocks into it.
    concat = [x.arr, x.arr] if one_d else [[x.arr[None]], [x.arr[None]]]
    plain_np_result = concatenate3(concat)
    mock_concatenate2.assert_not_called()
    assert type(plain_np_result) is np.ndarray
