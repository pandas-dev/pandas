import pytest


def pytest_addoption(parser):
    parser.addoption("--skip-slow", action="store_true",
                     help="skip slow tests")
    parser.addoption("--skip-network", action="store_true",
                     help="run network tests")
    parser.addoption("--run-disabled", action="store_false",
                     help="run disabled tests")


def pytest_runtest_setup(item):
    if 'slow' in item.keywords and item.config.getoption("--skip-slow"):
        pytest.skip("skipping due to --skip-slow")

    if 'skip' in item.keywords and item.config.getoption("--skip-network"):
        pytest.skip("skipping due to --skip-network")

    if 'disabled' in item.keywords and item.config.getoption("--run-disabled"):
        pytest.skip("need --run-disabled option to run")
