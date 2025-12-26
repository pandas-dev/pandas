import uuid

import pytest

from pandas.io.pytables import HDFStore


@pytest.fixture
def temp_h5_path(tmp_path):
    """Fixture for HDF5 path"""
    file_path = tmp_path / f"{uuid.uuid4()}.h5"
    file_path.touch()
    return file_path


@pytest.fixture
def temp_hdfstore(temp_h5_path):
    with HDFStore(temp_h5_path, mode="a") as store:
        yield store
