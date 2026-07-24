import sys



def is_pytest_compatible() -> bool:
    """ returns true if executing can be used for expressions inside assert statements which are rewritten by pytest
    """
    if sys.version_info < (3, 11):
        return False

    try:
        import pytest
    except ImportError:
        return False

    return pytest.version_tuple >= (8, 3, 4)
