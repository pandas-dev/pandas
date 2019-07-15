import pytest

import pandas


def test_matplotlib_backend_error():
    msg = (
        "matplotlib is required for plotting when the default backend "
        '"matplotlib" is selected.'
    )
    try:
        import matplotlib  # noqa
    except ImportError:
        with pytest.raises(ImportError, match=msg):
            pandas.set_option("plotting.backend", "matplotlib")


def test_backend_is_not_module():
    msg = (
        '"not_an_existing_module" does not seem to be an installed module. '
        "A pandas plotting backend must be a module that can be imported"
    )
    with pytest.raises(ValueError, match=msg):
        pandas.set_option("plotting.backend", "not_an_existing_module")


def test_backend_is_correct(monkeypatch):
    monkeypatch.setattr(
        "pandas.core.config_init.importlib.import_module", lambda name: None
    )
    pandas.set_option("plotting.backend", "correct_backend")
    assert pandas.get_option("plotting.backend") == "correct_backend"

    # Restore backend for other tests (matplotlib can be not installed)
    try:
        pandas.set_option("plotting.backend", "matplotlib")
    except ImportError:
        pass
