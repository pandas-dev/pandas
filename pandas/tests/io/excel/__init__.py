import pytest

pytestmark = pytest.mark.filterwarnings(
    # Looks like tree.getiterator is deprecated in favor of tree.iter
    "ignore:This method will be removed in future versions:PendingDeprecationWarning"
)
