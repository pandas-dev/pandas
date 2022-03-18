import numpy as np

from pandas._libs.tslibs.np_datetime import py_get_unit_from_dtype


def test_get_unit_from_dtype():
    # datetime64
    assert py_get_unit_from_dtype(np.dtype("M8[Y]")) == 0
    assert py_get_unit_from_dtype(np.dtype("M8[M]")) == 1
    assert py_get_unit_from_dtype(np.dtype("M8[W]")) == 2
    # B has been deprecated and removed -> no 3
    assert py_get_unit_from_dtype(np.dtype("M8[D]")) == 4
    assert py_get_unit_from_dtype(np.dtype("M8[h]")) == 5
    assert py_get_unit_from_dtype(np.dtype("M8[m]")) == 6
    assert py_get_unit_from_dtype(np.dtype("M8[s]")) == 7
    assert py_get_unit_from_dtype(np.dtype("M8[ms]")) == 8
    assert py_get_unit_from_dtype(np.dtype("M8[us]")) == 9
    assert py_get_unit_from_dtype(np.dtype("M8[ns]")) == 10
    assert py_get_unit_from_dtype(np.dtype("M8[ps]")) == 11
    assert py_get_unit_from_dtype(np.dtype("M8[fs]")) == 12
    assert py_get_unit_from_dtype(np.dtype("M8[as]")) == 13

    # timedelta64
    assert py_get_unit_from_dtype(np.dtype("m8[Y]")) == 0
    assert py_get_unit_from_dtype(np.dtype("m8[M]")) == 1
    assert py_get_unit_from_dtype(np.dtype("m8[W]")) == 2
    # B has been deprecated and removed -> no 3
    assert py_get_unit_from_dtype(np.dtype("m8[D]")) == 4
    assert py_get_unit_from_dtype(np.dtype("m8[h]")) == 5
    assert py_get_unit_from_dtype(np.dtype("m8[m]")) == 6
    assert py_get_unit_from_dtype(np.dtype("m8[s]")) == 7
    assert py_get_unit_from_dtype(np.dtype("m8[ms]")) == 8
    assert py_get_unit_from_dtype(np.dtype("m8[us]")) == 9
    assert py_get_unit_from_dtype(np.dtype("m8[ns]")) == 10
    assert py_get_unit_from_dtype(np.dtype("m8[ps]")) == 11
    assert py_get_unit_from_dtype(np.dtype("m8[fs]")) == 12
    assert py_get_unit_from_dtype(np.dtype("m8[as]")) == 13
