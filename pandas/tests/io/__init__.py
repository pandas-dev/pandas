import pytest

pytestmark = [
    # fastparquet
    pytest.mark.filterwarnings(
        "ignore:PY_SSIZE_T_CLEAN will be required.*:DeprecationWarning"
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
    # GH 26552
    pytest.mark.filterwarnings(
        "ignore:As the xlwt package is no longer maintained:FutureWarning"
    ),
]
