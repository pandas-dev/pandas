import sys
import types

import pkg_resources
import pytest

import pandas.util._test_decorators as td

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


@td.skip_if_no_mpl
def test_register_entrypoint():
    mod = types.ModuleType("my_backend")
    mod.plot = lambda *args, **kwargs: 1

    backends = pkg_resources.get_entry_map("pandas")
    my_entrypoint = pkg_resources.EntryPoint(
        "pandas_plotting_backend",
        mod.__name__,
        dist=pkg_resources.get_distribution("pandas"),
    )
    backends["pandas_plotting_backends"]["my_backend"] = my_entrypoint
    # TODO: the docs recommend importlib.util.module_from_spec. But this works for now.
    sys.modules["my_backend"] = mod

    result = pandas.plotting._core._get_plot_backend("my_backend")
    assert result is mod

    # TODO: https://github.com/pandas-dev/pandas/issues/27517
    # Remove the td.skip_if_no_mpl
    with pandas.option_context("plotting.backend", "my_backend"):
        result = pandas.plotting._core._get_plot_backend()

    assert result is mod


def test_register_import():
    mod = types.ModuleType("my_backend2")
    mod.plot = lambda *args, **kwargs: 1
    sys.modules["my_backend2"] = mod

    result = pandas.plotting._core._get_plot_backend("my_backend2")
    assert result is mod


@td.skip_if_mpl
def test_no_matplotlib_ok():
    with pytest.raises(ImportError):
        pandas.plotting._core._get_plot_backend("matplotlib")
