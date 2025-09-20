import functools
import sys
import traceback
from contextlib import contextmanager
from types import TracebackType

import pytest
from pytestqt.utils import get_marker

CapturedException = tuple[type[BaseException], BaseException, TracebackType]
CapturedExceptions = list[CapturedException]


@contextmanager
def capture_exceptions():
    """
    Context manager that captures exceptions that happen insides its context,
    and returns them as a list of (type, value, traceback) after the
    context ends.
    """
    manager = _QtExceptionCaptureManager()
    manager.start()
    try:
        yield manager.exceptions
    finally:
        manager.finish()


def _except_hook(type_, value, tback, exceptions=None):
    """Hook functions installed by _QtExceptionCaptureManager"""
    exceptions.append((type_, value, tback))
    sys.stderr.write(format_captured_exceptions([(type_, value, tback)]))


class _QtExceptionCaptureManager:
    """
    Manages exception capture context.
    """

    def __init__(self):
        self.old_hook = None
        self.exceptions = []

    def start(self):
        """Start exception capturing by installing a hook into sys.excepthook
        that records exceptions received into ``self.exceptions``.
        """
        self.old_hook = sys.excepthook
        sys.excepthook = functools.partial(_except_hook, exceptions=self.exceptions)

    def finish(self):
        """Stop exception capturing, restoring the original hook.

        Can be called multiple times.
        """
        if self.old_hook is not None:
            sys.excepthook = self.old_hook
            self.old_hook = None

    def fail_if_exceptions_occurred(self, when):
        """calls pytest.fail() with an informative message if exceptions
        have been captured so far. Before pytest.fail() is called, also
        finish capturing.
        """
        if self.exceptions:
            self.finish()
            exceptions = self.exceptions
            self.exceptions = []
            prefix = "%s ERROR: " % when
            msg = prefix + format_captured_exceptions(exceptions)
            del exceptions[:]  # Don't keep exceptions alive longer.
            if hasattr(sys, "last_exc"):
                sys.last_exc = None
            pytest.fail(msg, pytrace=False)


def format_captured_exceptions(exceptions):
    """
    Formats exceptions given as (type, value, traceback) into a string
    suitable to display as a test failure.
    """
    from io import StringIO

    stream = StringIO()
    stream.write("Exceptions caught in Qt event loop:\n")
    sep = "_" * 80 + "\n"
    stream.write(sep)
    for exc_type, value, tback in exceptions:
        traceback.print_exception(exc_type, value, tback, file=stream)
        stream.write(sep)
    return stream.getvalue()


def _is_exception_capture_enabled(item):
    """returns if exception capture is disabled for the given test item."""
    disabled = get_marker(item, "qt_no_exception_capture") or item.config.getini(
        "qt_no_exception_capture"
    )
    return not disabled


class TimeoutError(Exception):
    """
    .. versionadded:: 2.1

    Exception thrown by :class:`pytestqt.qtbot.QtBot` methods.

    Access via ``qtbot.TimeoutError``.
    """


class ScreenshotError(Exception):
    """
    .. versionadded:: 4.1

    Exception thrown by :meth:`pytestqt.qtbot.QtBot.screenshot` if taking the
    screenshot failed.

    .. versionchanged:: 4.2

        Access via ``qtbot.ScreenshotError``.
    """
