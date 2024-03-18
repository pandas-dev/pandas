"""
Just a conftest file for doctest stuff

The main conftest file is in pandas/tests/conftest.py
"""

import pytest


# https://github.com/pytest-dev/pytest/issues/11873
# Would like to avoid autouse=True, but cannot as of pytest 8.0.0
@pytest.fixture(autouse=True)
def add_doctest_imports(doctest_namespace) -> None:
    """
    Make `np` and `pd` names available for doctests.
    """
    import numpy as np

    import pandas as pd

    doctest_namespace["np"] = np
    doctest_namespace["pd"] = pd
