import pytest

import pandas.util._test_decorators as td

import pandas as pd


@td.skip_if_installed("numba")
def test_numba_not_installed_option_context():
    with pytest.raises(ImportError, match="`Import numba` failed"):
        with pd.option_context("compute.use_numba", True):
            pass
