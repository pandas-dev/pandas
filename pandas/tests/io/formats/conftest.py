import pytest


@pytest.fixture
def read_file(datapath):
    """fixture factory to read text files from tests/io/formats/data"""
    def _read_file(filename):
        filepath = datapath('io', 'formats', 'data', filename)
        with open(filepath) as f:
            contents = f.read()
        return contents
    return _read_file
