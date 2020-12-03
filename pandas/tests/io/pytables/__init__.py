import pytest

pytestmark = [
    # pytables https://github.com/PyTables/PyTables/issues/822
    pytest.mark.filterwarnings(
        "ignore:a closed node found in the registry:UserWarning"
    ),
    pytest.mark.filterwarnings(r"ignore:tostring\(\) is deprecated:DeprecationWarning"),
]
