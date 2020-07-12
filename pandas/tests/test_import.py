import pytest


@pytest.mark.parametrize("action", ["error", "always"])
def test_import_error(action):
    # GH#26367 Verify whether Pandas import succeeds when setting filterwarnings
    import warnings

    warnings.filterwarnings(action)
    import pandas as pd  # noqa: F401
