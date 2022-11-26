import pytest

pytestmark = [
    # fastparquet
    pytest.mark.filterwarnings(
        "ignore:PY_SSIZE_T_CLEAN will be required.*:DeprecationWarning"
    ),
    pytest.mark.filterwarnings(
        r"ignore:`np\.bool` is a deprecated alias:DeprecationWarning"
    ),
    # xlrd
    pytest.mark.filterwarnings(
        "ignore:This method will be removed in future versions:DeprecationWarning"
    ),
    pytest.mark.filterwarnings(
        "ignore:This method will be removed in future versions.  "
        r"Use 'tree.iter\(\)' or 'list\(tree.iter\(\)\)' instead."
        ":PendingDeprecationWarning"
    ),
]
