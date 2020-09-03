import pytest

import pandas as pd
import pandas.testing as tm


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

    with pytest.raises(ValueError):
        pd.options.mode.string_storage = "foo"
