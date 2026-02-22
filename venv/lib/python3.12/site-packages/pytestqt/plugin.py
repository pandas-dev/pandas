import warnings

import pytest

from pytestqt.exceptions import (
    _is_exception_capture_enabled,
    _QtExceptionCaptureManager,
)
from pytestqt.logging import QtLoggingPlugin, _QtMessageCapture
from pytestqt.qt_compat import qt_api
from pytestqt.qtbot import QtBot, _close_widgets


@pytest.fixture(scope="session")
def qapp_args(pytestconfig):
    """
    Fixture that provides QApplication arguments to use.

    You can override this fixture to pass different arguments to
    ``QApplication``:

    .. code-block:: python

       @pytest.fixture(scope="session")
       def qapp_args():
           return ["prog_name", "--arg=foo"]


    Note that it can only be overridden once at session scope.
    It is not possible to override this per unit test since a QApplication
    cannot be destroyed and recreated within the same app.

    The default value is a list with one element which is determined the same
    way as for ``QApplication.applicationName()``,
    see :ref:`qapp fixture<setting-qapp-name>` for more information.

    """
    return [pytestconfig.getini("qt_qapp_name")]


@pytest.fixture(scope="session")
def qapp_cls():
    """
    Fixture that provides the QApplication subclass to use.

    You can override this fixture to use a custom QApplication subclass from
    your application for tests:

    .. code-block:: python

       @pytest.fixture(scope="session")
       def qapp_cls():
           return myapp.Application

    Or use a ``QCoreApplication`` if you want to test a non-gui Qt application:

       @pytest.fixture(scope="session")
       def qapp_cls():
           return qt_api.QtCore.QCoreApplication
    """
    return qt_api.QtWidgets.QApplication


@pytest.fixture(scope="session")
def qapp(qapp_args, qapp_cls, pytestconfig):
    """
    Fixture that instantiates the QApplication instance that will be used by
    the tests.

    You can use the ``qapp`` fixture in tests which require a ``QApplication``
    to run, but where you don't need full ``qtbot`` functionality.
    """
    app = qt_api.QtWidgets.QApplication.instance()
    if app is None:
        global _qapp_instance
        _qapp_instance = qapp_cls(qapp_args)
        name = pytestconfig.getini("qt_qapp_name")
        _qapp_instance.setApplicationName(name)
        return _qapp_instance
    else:
        if not isinstance(app, qapp_cls):
            warnings.warn(
                f"Existing QApplication {app} is not an instance of qapp_cls: "
                f"{qapp_cls}"
            )
        return app


# holds a global QApplication instance created in the qapp fixture; keeping
# this reference alive avoids it being garbage collected too early
_qapp_instance = None


@pytest.fixture
def qtbot(qapp, request):
    """
    Fixture used to create a QtBot instance for using during testing.

    Make sure to call addWidget for each top-level widget you create to ensure
    that they are properly closed after the test ends.
    """
    result = QtBot(request)
    return result


@pytest.fixture
def qtlog(request):
    """Fixture that can access messages captured during testing"""
    if hasattr(request._pyfuncitem, "qt_log_capture"):
        return request._pyfuncitem.qt_log_capture
    else:
        return _QtMessageCapture([])  # pragma: no cover


@pytest.fixture
def qtmodeltester(request):
    """
    Fixture used to create a ModelTester instance to test models.
    """
    from pytestqt.modeltest import ModelTester

    tester = ModelTester(request.config)
    yield tester
    tester._cleanup()


def pytest_addoption(parser):
    parser.addini("qt_api", 'Qt api version to use: "pyside6" , "pyqt6", "pyqt5"')
    parser.addini("qt_no_exception_capture", "disable automatic exception capture")
    parser.addini(
        "qt_default_raising",
        "Default value for the raising parameter of qtbot.waitSignal/waitCallback",
    )
    parser.addini(
        "qt_qapp_name", "The Qt application name to use", default="pytest-qt-qapp"
    )

    default_log_fail = QtLoggingPlugin.LOG_FAIL_OPTIONS[0]
    parser.addini(
        "qt_log_level_fail",
        'log level in which tests can fail: {} (default: "{}")'.format(
            QtLoggingPlugin.LOG_FAIL_OPTIONS, default_log_fail
        ),
        default=default_log_fail,
    )
    parser.addini(
        "qt_log_ignore",
        "list of regexes for messages that should not cause a tests " "to fails",
        type="linelist",
    )

    group = parser.getgroup("qt", "qt testing")
    group.addoption(
        "--no-qt-log",
        dest="qt_log",
        action="store_false",
        default=True,
        help="disable pytest-qt logging capture",
    )
    group.addoption(
        "--qt-log-format",
        dest="qt_log_format",
        default=None,
        help="defines how qt log messages are displayed.",
    )


@pytest.hookimpl(wrapper=True, tryfirst=True)
def pytest_runtest_setup(item):
    """
    Hook called after before test setup starts, to start capturing exceptions
    as early as possible.
    """
    capture_enabled = _is_exception_capture_enabled(item)
    if capture_enabled:
        item.qt_exception_capture_manager = _QtExceptionCaptureManager()
        item.qt_exception_capture_manager.start()
    result = yield
    _process_events()
    if capture_enabled:
        item.qt_exception_capture_manager.fail_if_exceptions_occurred("SETUP")
    return result


@pytest.hookimpl(wrapper=True, tryfirst=True)
def pytest_runtest_call(item):
    result = yield
    _process_events()
    capture_enabled = _is_exception_capture_enabled(item)
    if capture_enabled:
        item.qt_exception_capture_manager.fail_if_exceptions_occurred("CALL")
    return result


@pytest.hookimpl(wrapper=True, trylast=True)
def pytest_runtest_teardown(item):
    """
    Hook called after each test tear down, to process any pending events and
    avoiding leaking events to the next test. Also, if exceptions have
    been captured during fixtures teardown, fail the test.
    """
    _process_events()
    _close_widgets(item)
    _process_events()
    result = yield
    _process_events()
    capture_enabled = _is_exception_capture_enabled(item)
    if capture_enabled:
        item.qt_exception_capture_manager.fail_if_exceptions_occurred("TEARDOWN")
        item.qt_exception_capture_manager.finish()
    return result


def _process_events():
    """Calls app.processEvents() while taking care of capturing exceptions
    or not based on the given item's configuration.
    """
    app = qt_api.QtWidgets.QApplication.instance()
    if app is not None:
        app.processEvents()


def pytest_configure(config):
    config.addinivalue_line(
        "markers",
        "qt_no_exception_capture: Disables pytest-qt's automatic exception "
        "capture for just one test item.",
    )

    config.addinivalue_line(
        "markers", "qt_log_level_fail: overrides qt_log_level_fail ini option."
    )
    config.addinivalue_line(
        "markers", "qt_log_ignore: overrides qt_log_ignore ini option."
    )
    config.addinivalue_line("markers", "no_qt_log: Turn off Qt logging capture.")

    if config.getoption("qt_log") and config.getoption("capture") != "no":
        config.pluginmanager.register(QtLoggingPlugin(config), "_qt_logging")

    qt_api.set_qt_api(config.getini("qt_api"))


def pytest_report_header():
    from pytestqt.qt_compat import qt_api

    v = qt_api.get_versions()
    fields = [
        f"{v.qt_api} {v.qt_api_version}",
        "Qt runtime %s" % v.runtime,
        "Qt compiled %s" % v.compiled,
    ]
    version_line = " -- ".join(fields)
    return [version_line]
