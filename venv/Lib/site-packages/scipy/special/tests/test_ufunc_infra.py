import numpy as np
import pytest

from scipy.special import mathieu_sem
# from scipy.special._ufuncs import _mathieu_sem


@pytest.mark.parametrize("m_shape,q_shape,x_shape,where_shape", [
    ((5, 1), (1, 1), (1, 10), ()),
    ((2, 3), (2, 3), (2, 3), ()),
    ((1, 5), (1, 5), (10, 5), ()),
    ((5, 1), (1, 1), (1, 10), (1, 10)),
    ((2, 3), (2, 3), (2, 3), (2, 3)),
    ((1, 5), (1, 5), (10, 5), (10, 5)),
    ((5, 1), (1, 1), (1, 10), (5, 1)),
])
def test_out(m_shape, q_shape, x_shape, where_shape):
    rng = np.random.default_rng(1234)

    m = rng.integers(1, 20, m_shape)
    q = rng.uniform(0, 10, q_shape)
    x = rng.uniform(0, 90, x_shape)

    if where_shape == ():
        batch_shape = np.broadcast_shapes(m_shape, q_shape, x_shape)
        where = True
    else:
        where = rng.choice([True, False], size=where_shape)
        batch_shape = np.broadcast_shapes(m_shape, q_shape, x_shape, where_shape)

    out0 = np.empty(batch_shape)
    out1 = np.empty(batch_shape)
    res0, res1 = mathieu_sem(m, q, x, out=(out0, out1), where=where)
    assert res0 is out0 and res1 is out1
