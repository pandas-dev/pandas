import uuid

import pytest


@pytest.fixture(name="setup_path")
def fixture_setup_path():
    """Fixture for setup path"""
    return f"tmp.__{uuid.uuid4()}__.h5"
