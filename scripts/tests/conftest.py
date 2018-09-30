def pytest_addoption(parser):
    parser.addoption("--strict-data-files", action="store_true",
                     help="Unused. For compat with setup.cfg.")
