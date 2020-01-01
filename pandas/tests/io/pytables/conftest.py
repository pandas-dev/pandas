import pytest

import pandas.util.testing as tm


@pytest.fixture
def setup_path():
    """Fixture for setup path"""
    return "tmp.__{}__.h5".format(tm.rands(10))


@pytest.fixture(scope="module", autouse=True)
def setup_mode():
    """ Reset testing mode fixture"""
    tm.reset_testing_mode()
    yield
    tm.set_testing_mode()
