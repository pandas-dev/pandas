import importlib

import pytest

import pandas as pd
import pandas._testing as tm


def test_configs_deprecated():
    # GH#57872
    with tm.assert_produces_warning(FutureWarning, match="is deprecated"):
        pd.options.mode.copy_on_write = True

    with tm.assert_produces_warning(FutureWarning, match="is deprecated"):
        pd.options.future.no_silent_downcasting = True


def test_env_variable():
    # GH#57872
    pytest.MonkeyPatch().setenv("PANDAS_COPY_ON_WRITE", "1")
    with tm.assert_produces_warning(
        FutureWarning, match="deprecated", check_stacklevel=False
    ):
        importlib.reload(pd)
