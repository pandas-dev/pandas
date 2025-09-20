
import io
import logging
import pytest

from traitlets import default
from .mockextension import MockExtensionApp
from notebook_shim import shim


@pytest.fixture
def read_app_logs(capsys):
    """Fixture that returns a callable to read
    the current output from the application's logs
    that was printed to sys.stderr.
    """
    def _inner():
        captured = capsys.readouterr()
        return captured.err
    return _inner


@pytest.fixture
def jp_server_config(capsys):
    return {
        "ServerApp": {
            "jpserver_extensions": {
                "notebook_shim": True,
                "notebook_shim.tests.mockextension": True
            }
        }
    }


@pytest.fixture
def extensionapp(jp_serverapp):
    return jp_serverapp.extension_manager.extension_points["mockextension"].app


def list_test_params(param_input):
    """"""
    params = []
    for test in param_input:
        name, value = test[0], test[1]
        option = (
            '--MockExtensionApp.'
            '{name}={value}'
            .format(name=name, value=value)
        )
        params.append([[option], name, value])
    return params


@pytest.mark.parametrize(
    'jp_argv,trait_name,trait_value',
    list_test_params([
        ('enable_mathjax', False)
    ])
)
def test_EXTAPP_AND_NBAPP_SHIM_MSG(
    read_app_logs,
    extensionapp,
    jp_argv,
    trait_name,
    trait_value
):
    log = read_app_logs()
    # Verify a shim warning appeared.
    log_msg = shim.EXTAPP_AND_NBAPP_SHIM_MSG(trait_name, 'MockExtensionApp')
    assert log_msg in log
    # Verify the trait was updated.
    assert getattr(extensionapp, trait_name) == trait_value


@pytest.mark.parametrize(
    'jp_argv,trait_name,trait_value',
    list_test_params([
        ('allow_origin', ''),
        ('allow_origin_pat', ''),
    ])
)
def test_EXTAPP_AND_SVAPP_SHIM_MSG(
    read_app_logs,
    extensionapp,
    jp_argv,
    trait_name,
    trait_value
):
    log = read_app_logs()
    # Verify a shim warning appeared.
    log_msg = shim.EXTAPP_AND_SVAPP_SHIM_MSG(trait_name, 'MockExtensionApp')
    assert log_msg in log
    # Verify the trait was updated.
    assert getattr(extensionapp, trait_name) == trait_value


@pytest.mark.parametrize(
    'jp_argv,trait_name,trait_value',
    list_test_params([
        ('jinja_environment_options', {}),
        ('jinja_template_vars', {}),
        ('extra_template_paths', []),
        ('quit_button', True),
    ])
)
def test_NOT_EXTAPP_NBAPP_AND_SVAPP_SHIM_MSG(
    read_app_logs,
    extensionapp,
    jp_argv,
    trait_name,
    trait_value
):
    log = read_app_logs()
    # Verify a shim warning appeared.
    log_msg = shim.NOT_EXTAPP_NBAPP_AND_SVAPP_SHIM_MSG(trait_name, 'MockExtensionApp')
    assert log_msg in log
    # Verify the trait was updated.
    assert getattr(extensionapp.serverapp, trait_name) == trait_value


@pytest.mark.parametrize(
    'jp_argv,trait_name,trait_value',
    list_test_params([
        ('allow_credentials', False),
    ])
)
def test_EXTAPP_TO_SVAPP_SHIM_MSG(
    read_app_logs,
    extensionapp,
    jp_argv,
    trait_name,
    trait_value
):
    log = read_app_logs()
    # Verify a shim warning appeared.
    log_msg = shim.EXTAPP_TO_SVAPP_SHIM_MSG(trait_name, 'MockExtensionApp')
    assert log_msg in log
    # Verify the trait was updated.
    assert getattr(extensionapp.serverapp, trait_name) == trait_value


@pytest.mark.parametrize(
    'jp_argv,trait_name,trait_value',
    list_test_params([
        ('mathjax_config', 'TEST'),
        ('mathjax_url', 'TEST')
    ])
)
def test_EXTAPP_TO_NBAPP_SHIM_MSG(
    read_app_logs,
    extensionapp,
    jp_argv,
    trait_name,
    trait_value
):
    log = read_app_logs()
    # Verify a shim warning appeared.
    log_msg = shim.EXTAPP_TO_NBAPP_SHIM_MSG(trait_name, 'MockExtensionApp')
    assert log_msg in log
