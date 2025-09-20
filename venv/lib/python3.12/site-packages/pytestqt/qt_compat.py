"""
Provide a common way to import Qt classes used by pytest-qt in a unique manner,
abstracting API differences between PyQt5/6 and PySide6.

.. note:: This module is not part of pytest-qt public API, hence its interface
may change between releases and users should not rely on it.

Based on from https://github.com/epage/PythonUtils.
"""

from collections import namedtuple, OrderedDict
import os
import sys

import pytest


VersionTuple = namedtuple("VersionTuple", "qt_api, qt_api_version, runtime, compiled")

QT_APIS = OrderedDict()
QT_APIS["pyside6"] = "PySide6"
QT_APIS["pyqt6"] = "PyQt6"
QT_APIS["pyqt5"] = "PyQt5"


def _import(name):
    """Thin call so we can mock it during testing"""
    return __import__(name)


def _is_library_loaded(name):
    return name in sys.modules


class _QtApi:
    """
    Interface to the underlying Qt API currently configured for pytest-qt.

    This object lazily loads all class references and other objects when the ``set_qt_api`` method
    gets called, providing a uniform way to access the Qt classes.
    """

    def __init__(self):
        self._import_errors = {}

    def _get_qt_api_from_env(self):
        api = os.environ.get("PYTEST_QT_API")
        supported_apis = QT_APIS.keys()

        if api is not None:
            api = api.lower()
            if api not in supported_apis:  # pragma: no cover
                msg = f"Invalid value for $PYTEST_QT_API: {api}, expected one of {supported_apis}"
                raise pytest.UsageError(msg)
        return api

    def _get_already_loaded_backend(self):
        for api, backend in QT_APIS.items():
            if _is_library_loaded(backend):
                return api
        return None

    def _guess_qt_api(self):  # pragma: no cover
        def _can_import(name):
            try:
                _import(name)
                return True
            except ModuleNotFoundError as e:
                self._import_errors[name] = str(e)
                return False

        # Note, not importing only the root namespace because when uninstalling from conda,
        # the namespace can still be there.
        for api, backend in QT_APIS.items():
            if _can_import(f"{backend}.QtCore"):
                return api
        return None

    def set_qt_api(self, api):
        self.pytest_qt_api = (
            self._get_qt_api_from_env()
            or api
            or self._get_already_loaded_backend()
            or self._guess_qt_api()
        )

        self.is_pyside = self.pytest_qt_api in ["pyside6"]
        self.is_pyqt = self.pytest_qt_api in ["pyqt5", "pyqt6"]

        if not self.pytest_qt_api:  # pragma: no cover
            errors = "\n".join(
                f"  {module}: {reason}"
                for module, reason in sorted(self._import_errors.items())
            )
            msg = (
                "pytest-qt requires either PySide6, PyQt5 or PyQt6 installed.\n"
                + errors
            )
            raise pytest.UsageError(msg)

        _root_module = QT_APIS[self.pytest_qt_api]

        def _import_module(module_name):
            m = __import__(_root_module, globals(), locals(), [module_name], 0)
            return getattr(m, module_name)

        self.QtCore = QtCore = _import_module("QtCore")
        self.QtGui = _import_module("QtGui")
        self.QtTest = _import_module("QtTest")
        self.QtWidgets = _import_module("QtWidgets")

        self._check_qt_api_version()

        # qInfo is not exposed in PySide6 < 6.8.2 (#232)
        if hasattr(QtCore, "qInfo"):
            self.qInfo = QtCore.qInfo
        elif hasattr(QtCore, "QMessageLogger"):
            self.qInfo = lambda msg: QtCore.QMessageLogger().info(msg)
        else:
            self.qInfo = None

        self.qDebug = QtCore.qDebug
        self.qWarning = QtCore.qWarning
        self.qCritical = QtCore.qCritical
        self.qFatal = QtCore.qFatal

        if self.is_pyside:
            self.Signal = QtCore.Signal
            self.Slot = QtCore.Slot
            self.Property = QtCore.Property
        elif self.is_pyqt:
            self.Signal = QtCore.pyqtSignal
            self.Slot = QtCore.pyqtSlot
            self.Property = QtCore.pyqtProperty
        else:
            assert False, "Expected either is_pyqt or is_pyside"

    def _check_qt_api_version(self):
        if not self.is_pyqt:
            # We support all PySide versions
            return

        if self.QtCore.PYQT_VERSION == 0x060000:  # 6.0.0
            raise pytest.UsageError(
                "PyQt 6.0 is not supported by pytest-qt, use 6.1+ instead."
            )
        elif self.QtCore.PYQT_VERSION < 0x050B00:  # 5.11.0
            raise pytest.UsageError(
                "PyQt < 5.11 is not supported by pytest-qt, use 5.11+ instead."
            )

    def exec(self, obj, *args, **kwargs):
        # exec was a keyword in Python 2, so PySide6 6.0
        # names the corresponding method "exec_" instead.
        #
        # The old _exec() alias is removed in PyQt6 and also deprecated as of
        # PySide 6.1:
        # https://codereview.qt-project.org/c/pyside/pyside-setup/+/342095
        if hasattr(obj, "exec"):
            return obj.exec(*args, **kwargs)
        return obj.exec_(*args, **kwargs)

    def get_versions(self):
        if self.pytest_qt_api == "pyside6":
            import PySide6  # type: ignore[import-not-found,unused-ignore]

            version = PySide6.__version__

            return VersionTuple(
                "PySide6", version, self.QtCore.qVersion(), self.QtCore.__version__
            )
        elif self.pytest_qt_api == "pyqt6":
            return VersionTuple(
                "PyQt6",
                self.QtCore.PYQT_VERSION_STR,
                self.QtCore.qVersion(),
                self.QtCore.QT_VERSION_STR,
            )
        elif self.pytest_qt_api == "pyqt5":
            return VersionTuple(
                "PyQt5",
                self.QtCore.PYQT_VERSION_STR,
                self.QtCore.qVersion(),
                self.QtCore.QT_VERSION_STR,
            )

        assert False, f"Internal error, unknown pytest_qt_api: {self.pytest_qt_api}"


qt_api = _QtApi()
