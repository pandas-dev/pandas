import os

import pytest


@pytest.fixture
def datapath(request):
    """Get the path to a data file.

    Parameters
    ----------
    path : str
        Path to the file, relative to ``pandas/tests/``

    Returns
    -------
    path : path including ``pandas/tests``.

    Raises
    ------
    ValueError
        If the path doesn't exist and the --strict-data-files option is set.
    """
    def deco(*args):
        path = os.path.join('pandas', 'tests', *args)
        if not os.path.exists(path):
            if request.config.getoption("--strict-data-files"):
                raise ValueError("Failed.")
            else:
                pytest.skip("{} not included in pandas distribution."
                            .format(path))
        return path
    return deco
