import pytest

from pandas.compat import (
    is_ci_environment,
    is_platform_windows,
)

pytestmark = pytest.mark.skipif(
    is_platform_windows() and is_ci_environment(), reason="Causes pytest INTERNALERROR"
)
