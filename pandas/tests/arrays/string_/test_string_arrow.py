import re

import numpy as np
import pytest

import pandas as pd
import pandas.testing as tm

pa = pytest.importorskip("pyarrow", minversion="1.0.0")

from pandas.core.arrays.string_arrow import ArrowStringArray


def test_eq_all_na():
    a = pd.array([pd.NA, pd.NA], dtype=pd.StringDtype("pyarrow"))
    result = a == a
    expected = pd.array([pd.NA, pd.NA], dtype="boolean")
    tm.assert_extension_array_equal(result, expected)


def test_config():
    # python by default
    assert pd.StringDtype().storage == "python"
    arr = pd.array(["a", "b"])
    assert arr.dtype.storage == "python"

    with pd.option_context("mode.string_storage", "pyarrow"):
        assert pd.StringDtype().storage == "pyarrow"
        arr = pd.array(["a", "b"])
        assert arr.dtype.storage == "pyarrow"

    msg = re.escape("Value must be one of python|pyarrow")
    with pytest.raises(ValueError, match=msg):
        pd.options.mode.string_storage = "foo"


@pytest.mark.parametrize("chunked", [True, False])
@pytest.mark.parametrize("np_or_pa", [np, pa])
def test_constructor_not_string_type_raises(np_or_pa, chunked):
    arr = np_or_pa.array([1, 2, 3])
    if chunked:
        if np_or_pa is np:
            pytest.skip("chunked not applicable to numpy array")
        arr = pa.chunked_array(arr)
    if np_or_pa is np:
        msg = "Unsupported type '<class 'numpy.ndarray'>' for ArrowStringArray"
    else:
        msg = re.escape(
            "ArrowStringArray requires a PyArrow (chunked) array of string type"
        )
    with pytest.raises(ValueError, match=msg):
        ArrowStringArray(arr)
