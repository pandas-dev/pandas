"""
Package for testing excel-related functionality.

All tests are marked as slow via ``pytestmark = pytest.mark.slow``
at the top of each test module.
"""

import pytest

pytestmark = [
    pytest.mark.filterwarnings(
        # Looks like tree.getiterator is deprecated in favor of tree.iter
        "ignore:This method will be removed in future versions:"
        "PendingDeprecationWarning"
    ),
    pytest.mark.filterwarnings(
        "ignore:This method will be removed in future versions:DeprecationWarning"
    ),
]
