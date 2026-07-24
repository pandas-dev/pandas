import toolz


def test_has_version():
    # If this test fails, then toolz probably isn't installed properly.
    # For local development, try `pip install -e .` from the project directory.
    version = toolz.__version__
    assert isinstance(version, str)
    assert version.startswith("1.")
