"""basic smoke tests of jupyterlite_core infrastructure"""

import jupyterlite_core


def test_is_documented():
    """TODO: improve the definition of documented"""
    assert jupyterlite_core.__doc__


def test_is_versioned():
    """TODO: test the version agrees with the version mangling from npm"""
    assert jupyterlite_core.__version__
