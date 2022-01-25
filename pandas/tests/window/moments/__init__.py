import os

import pytest

from pandas.compat import is_platform_windows

pytestmark = pytest.mark.skipif(
    os.environ.get("PANDAS_CI", "0") == "1" and is_platform_windows(),
    reason="Causes flaky timeouts possibly due to test teardown in the CI",
)
