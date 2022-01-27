import os

import pytest

from pandas.compat import is_platform_windows

pytestmark = pytest.mark.skipif(
    os.environ.get("PANDAS_CI", "0") == "1" and is_platform_windows(),
    reason="Any test in this directory can hang on the multi-process "
    "CI Windows environment",
)
