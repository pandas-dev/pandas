# pyproject.toml defines addopts: --strict-data-files
# strict-data-files is defined & used in pandas/conftest.py
def pytest_addoption(parser) -> None:
    parser.addoption(
        "--strict-data-files",
        action="store_true",
        help="Unused",
    )
