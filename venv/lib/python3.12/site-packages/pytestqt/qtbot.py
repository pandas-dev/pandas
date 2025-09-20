from collections.abc import Callable
import contextlib
from types import TracebackType
import weakref
import warnings
from typing import (
    TYPE_CHECKING,
    Generator,
    Iterator,
    Literal,
    Optional,
    Any,
    cast,
)
from pathlib import Path
import pytest
from typing_extensions import Self, TypeAlias

from pytestqt.exceptions import TimeoutError, ScreenshotError
from pytestqt.qt_compat import qt_api
from pytestqt.wait_signal import (
    SignalBlocker,
    MultiSignalBlocker,
    SignalEmittedSpy,
    SignalEmittedError,
    CallbackBlocker,
    CallbackCalledTwiceError,
    CheckParamsCb,
)

from pytest import FixtureRequest

# Type hint objects until figuring out how to import across qt
# versions possibly using 'qtpy' library.
QWidget: TypeAlias = Any
SignalInstance: TypeAlias = Any
QRect: TypeAlias = Any
QKeySequence: TypeAlias = Any

if TYPE_CHECKING:
    # Keep local import behavior the same.
    from pytestqt.exceptions import CapturedExceptions

BeforeCloseFunc = Callable[[QWidget], None]
WaitSignalsOrder = Literal["none", "simple", "strict"]


def _parse_ini_boolean(value: Any) -> bool:
    if value in (True, False):
        return cast("bool", value)
    try:
        return {"true": True, "false": False}[str(value).lower()]
    except KeyError:
        raise ValueError("unknown string for bool: %r" % value)


class QtBot:
    """
    Instances of this class are responsible for sending events to `Qt` objects (usually widgets),
    simulating user input.

    .. important:: Instances of this class should be accessed only by using a ``qtbot`` fixture,
                    never instantiated directly.

    **Widgets**

    .. automethod:: addWidget
    .. automethod:: captureExceptions
    .. automethod:: waitActive
    .. automethod:: waitExposed
    .. automethod:: waitForWindowShown
    .. automethod:: stop
    .. automethod:: screenshot
    .. automethod:: wait

    **Signals and Events**

    .. automethod:: waitSignal
    .. automethod:: waitSignals
    .. automethod:: assertNotEmitted
    .. automethod:: waitUntil

    **Raw QTest API**

    Methods below provide very low level functions, as sending a single mouse click or a key event.
    Those methods are just forwarded directly to the `QTest API`_. Consult the documentation for more
    information.

    .. note::
        These methods should be rarely be used, in general prefer to interact with widgets
        using their own methods such as ``QComboBox.setCurrentText``, ``QLineEdit.setText``, etc.
        Doing so will have the same effect as users interacting with the widget, but are more reliable.

        See :ref:`this note in the tutorial <note-about-qtbot-methods>` for more information.

    ---

    Below are methods used to simulate sending key events to widgets:

    .. staticmethod:: keyClick (widget, key[, modifier=Qt.KeyboardModifier.NoModifier[, delay=-1]])
    .. staticmethod:: keyClicks (widget, key_sequence[, modifier=Qt.KeyboardModifier.NoModifier[, delay=-1]])
    .. staticmethod:: keyEvent (action, widget, key[, modifier=Qt.KeyboardModifier.NoModifier[, delay=-1]])
    .. staticmethod:: keyPress (widget, key[, modifier=Qt.KeyboardModifier.NoModifier[, delay=-1]])
    .. staticmethod:: keyRelease (widget, key[, modifier=Qt.KeyboardModifier.NoModifier[, delay=-1]])

        Sends one or more keyboard events to a widget.

        :param QWidget widget: the widget that will receive the event

        :param str|int key: key to send, it can be either a Qt.Key.Key_* constant or a single character string.

        .. _keyboard modifiers:

        :param Qt.KeyboardModifier modifier: flags OR'ed together representing other modifier keys
            also pressed. Possible flags are:

            * ``Qt.KeyboardModifier.NoModifier``: No modifier key is pressed.
            * ``Qt.KeyboardModifier.ShiftModifier``: A Shift key on the keyboard is pressed.
            * ``Qt.KeyboardModifier.ControlModifier``: A Ctrl key on the keyboard is pressed.
            * ``Qt.KeyboardModifier.AltModifier``: An Alt key on the keyboard is pressed.
            * ``Qt.KeyboardModifier.MetaModifier``: A Meta key on the keyboard is pressed.
            * ``Qt.KeyboardModifier.KeypadModifier``: A keypad button is pressed.
            * ``Qt.KeyboardModifier.GroupSwitchModifier``: X11 only. A Mode_switch key on the keyboard is pressed.

        :param int delay: after the event, delay the test for this milliseconds (if > 0).


    .. staticmethod:: keyToAscii (key)

        Auxiliary method that converts the given constant to its equivalent ascii.

        :param Qt.Key.Key_* key: one of the constants for keys in the Qt namespace.

        :return type: str
        :returns: the equivalent character string.

        .. note:: This method is not available in PyQt.

    ---

    Below are methods used to simulate sending mouse events to widgets.

    .. staticmethod:: mouseClick (widget, button[, modifier=0[, pos=QPoint()[, delay=-1]]])
    .. staticmethod:: mouseDClick (widget, button[, modifier=0[, pos=QPoint()[, delay=-1]]])
    .. staticmethod:: mouseMove (widget[, pos=QPoint()[, delay=-1]])
    .. staticmethod:: mousePress (widget, button[, modifier=0[, pos=QPoint()[, delay=-1]]])
    .. staticmethod:: mouseRelease (widget, button[, modifier=0[, pos=QPoint()[, delay=-1]]])

        Sends a mouse moves and clicks to a widget.

        :param QWidget widget: the widget that will receive the event

        :param Qt.MouseButton button: flags OR'ed together representing the button pressed.
            Possible flags are:

            * ``Qt.MouseButton.NoButton``: The button state does not refer to any button
              (see QMouseEvent.button()).
            * ``Qt.MouseButton.LeftButton``: The left button is pressed, or an event refers to the left button.
              (The left button may be the right button on left-handed mice.)
            * ``Qt.MouseButton.RightButton``: The right button.
            * ``Qt.MouseButton.MidButton``: The middle button.
            * ``Qt.MouseButton.MiddleButton``: The middle button.
            * ``Qt.MouseButton.XButton1``: The first X button.
            * ``Qt.MouseButton.XButton2``: The second X button.

        :param Qt.KeyboardModifier modifier: flags OR'ed together representing other modifier keys
            also pressed. See `keyboard modifiers`_.

        :param QPoint pos: position of the mouse pointer.

        :param int delay: after the event, delay the test for this milliseconds (if > 0).

        .. note:: In the PySide bindings, the *modifier* argument is called *stateKey*.


    .. _QTest API: http://doc.qt.io/qt-5/qtest.html

    """

    def __init__(self, request: FixtureRequest) -> None:
        self._request = request
        # pep8 aliases. Set here to automatically use implementations defined in sub-classes for alias creation
        self.add_widget = self.addWidget
        self.capture_exceptions = self.captureExceptions
        self.wait_active = self.waitActive
        self.wait_exposed = self.waitExposed
        self.wait_for_window_shown = self.waitForWindowShown
        self.wait_signal = self.waitSignal
        self.wait_signals = self.waitSignals
        self.assert_not_emitted = self.assertNotEmitted
        self.wait_until = self.waitUntil
        self.wait_callback = self.waitCallback

    def _should_raise(self, raising_arg: Optional[bool]) -> bool:
        ini_val = self._request.config.getini("qt_default_raising")

        if raising_arg is not None:
            return raising_arg
        elif ini_val:
            return _parse_ini_boolean(ini_val)
        else:
            return True

    def addWidget(
        self, widget: QWidget, *, before_close_func: Optional[BeforeCloseFunc] = None
    ) -> None:
        """
        Adds a widget to be tracked by this bot. This is not required, but will ensure that the
        widget gets closed by the end of the test, so it is highly recommended.

        :param QWidget widget:
            Widget to keep track of.

        :kwparam before_close_func:
            A function that receives the widget as single parameter, which is called just before
            the ``.close()`` method gets called.

        .. note:: This method is also available as ``add_widget`` (pep-8 alias)
        """
        if not isinstance(widget, qt_api.QtWidgets.QWidget):
            raise TypeError(f"Need to pass a QWidget to addWidget: {widget!r}")
        _add_widget(self._request.node, widget, before_close_func=before_close_func)

    def waitActive(
        self, widget: QWidget, *, timeout: int = 5000
    ) -> "_WaitWidgetContextManager":
        """
        Context manager that waits for ``timeout`` milliseconds or until the window is active.
        If window is not exposed within ``timeout`` milliseconds, raise
        :class:`qtbot.TimeoutError <pytestqt.exceptions.TimeoutError>`

        This is mainly useful for asynchronous systems like X11, where a window will be mapped to screen
        some time after  being asked to show itself on the screen.

        .. code-block:: python

            with qtbot.waitActive(widget, timeout=500):
                show_action()

        :param QWidget widget:
            Widget to wait for.

        :param int|None timeout:
            How many milliseconds to wait for.

        .. note:: This method is also available as ``wait_active`` (pep-8 alias)
        """
        __tracebackhide__ = True
        return _WaitWidgetContextManager(
            "qWaitForWindowActive", "activated", widget, timeout
        )

    def waitExposed(
        self, widget: QWidget, *, timeout: int = 5000
    ) -> "_WaitWidgetContextManager":
        """
        Context manager that waits for ``timeout`` milliseconds or until the window is exposed.
        If the window is not exposed within ``timeout`` milliseconds, raise
        :class:`qtbot.TimeoutError <pytestqt.exceptions.TimeoutError>`

        This is mainly useful for asynchronous systems like X11, where a window will be mapped to screen
        some time after  being asked to show itself on the screen.

        .. code-block:: python

            with qtbot.waitExposed(splash, timeout=500):
                startup()

        :param QWidget widget:
            Widget to wait for.

        :param int|None timeout:
            How many milliseconds to wait for.

        .. note:: This method is also available as ``wait_exposed`` (pep-8 alias)
        """
        __tracebackhide__ = True
        return _WaitWidgetContextManager(
            "qWaitForWindowExposed", "exposed", widget, timeout
        )

    def waitForWindowShown(self, widget: QWidget) -> bool:
        """
        Waits until the window is shown in the screen. This is mainly useful for asynchronous
        systems like X11, where a window will be mapped to screen some time after being asked to
        show itself on the screen.

        .. warning::
            This method does **not** raise :class:`qtbot.TimeoutError <pytestqt.exceptions.TimeoutError>` if
            the window wasn't shown.

        .. deprecated:: 4.0
            Use the ``qtbot.waitExposed`` context manager instead.

        :param QWidget widget:
            Widget to wait on.

        :returns:
            ``True`` if the window was shown, ``False`` if ``.show()`` was never
            called or a timeout occurred.

        .. note:: This method is also available as ``wait_for_window_shown`` (pep-8 alias)
        """
        warnings.warn(
            "waitForWindowShown is deprecated, as the underlying Qt method was "
            "obsoleted in Qt 5.0 and removed in Qt 6.0. Its name is imprecise and "
            "the pytest-qt wrapper does not raise qtbot.TimeoutError if the window "
            "wasn't shown. Please use the qtbot.waitExposed context manager "
            "instead.",
            DeprecationWarning,
        )
        return qt_api.QtTest.QTest.qWaitForWindowExposed(widget)

    def stop(self) -> None:
        """
        Stops the current test flow, letting the user interact with any visible widget.

        This is mainly useful so that you can verify the current state of the program while writing
        tests.

        Closing the windows should resume the test run, with ``qtbot`` attempting to restore visibility
        of the widgets as they were before this call.
        """
        widget_and_visibility = []
        for weak_widget in _iter_widgets(self._request.node):
            widget = weak_widget()
            if widget is not None:
                widget_and_visibility.append((widget, widget.isVisible()))

        qt_api.exec(qt_api.QtWidgets.QApplication.instance())

        for widget, visible in widget_and_visibility:
            widget.setVisible(visible)

    def waitSignal(
        self,
        signal: SignalInstance,
        *,
        timeout: int = 5000,
        raising: Optional[bool] = None,
        check_params_cb: Optional[CheckParamsCb] = None,
    ) -> "SignalBlocker":
        """
        .. versionadded:: 1.2

        Stops current test until a signal is triggered.

        Used to stop the control flow of a test until a signal is emitted, or
        a number of milliseconds, specified by ``timeout``, has elapsed.

        Best used as a context manager::

           with qtbot.waitSignal(signal, timeout=1000):
               long_function_that_calls_signal()

        Also, you can use the :class:`SignalBlocker` directly if the context
        manager form is not convenient::

           blocker = qtbot.waitSignal(signal, timeout=1000)
           blocker.connect(another_signal)
           long_function_that_calls_signal()
           blocker.wait()

        Any additional signal, when triggered, will make :meth:`wait` return.

        .. versionadded:: 1.4
           The *raising* parameter.

        .. versionadded:: 2.0
           The *check_params_cb* parameter.

        :param Signal signal:
            A signal to wait for, or a tuple ``(signal, signal_name_as_str)`` to improve the error message that is part
            of :class:`qtbot.TimeoutError <pytestqt.exceptions.TimeoutError>`.
        :param int timeout:
            How many milliseconds to wait before resuming control flow.
        :param bool raising:
            If :class:`qtbot.TimeoutError <pytestqt.exceptions.TimeoutError>`
            should be raised if a timeout occurred.
            This defaults to ``True`` unless ``qt_default_raising = false``
            is set in the config.
        :param Callable check_params_cb:
            Optional ``callable`` that compares the provided signal parameters to some expected parameters.
            It has to match the signature of ``signal`` (just like a slot function would) and return ``True`` if
            parameters match, ``False`` otherwise.
        :returns:
            ``SignalBlocker`` object. Call ``SignalBlocker.wait()`` to wait.

        .. note::
            This method is also available as ``wait_signal`` (pep-8 alias)
        """
        if signal is None:
            raise ValueError(
                f"Passing None as signal isn't supported anymore, use qtbot.wait({timeout}) instead."
            )
        raising = self._should_raise(raising)
        blocker = SignalBlocker(
            timeout=timeout, raising=raising, check_params_cb=check_params_cb
        )
        blocker.connect(signal)
        return blocker

    def waitSignals(
        self,
        signals: list[SignalInstance],
        *,
        timeout: int = 5000,
        raising: Optional[bool] = None,
        check_params_cbs: Optional[list[CheckParamsCb]] = None,
        order: WaitSignalsOrder = "none",
    ) -> "MultiSignalBlocker":
        """
        .. versionadded:: 1.4

        Stops current test until all given signals are triggered.

        Used to stop the control flow of a test until all (and only
        all) signals are emitted or the number of milliseconds specified by
        ``timeout`` has elapsed.

        Best used as a context manager::

           with qtbot.waitSignals([signal1, signal2], timeout=1000):
               long_function_that_calls_signals()

        Also, you can use the :class:`MultiSignalBlocker` directly if the
        context manager form is not convenient::

           blocker = qtbot.waitSignals(signals, timeout=1000)
           long_function_that_calls_signal()
           blocker.wait()

        :param list signals:
            A list of :class:`Signal` objects to wait for. Alternatively: a list of (``Signal, str``) tuples of the form
            ``(signal, signal_name_as_str)`` to improve the error message that is part of ``qtbot.TimeoutError``.
        :param int timeout:
            How many milliseconds to wait before resuming control flow.
        :param bool raising:
            If :class:`qtbot.TimeoutError <pytestqt.exceptions.TimeoutError>`
            should be raised if a timeout occurred.
            This defaults to ``True`` unless ``qt_default_raising = false``
            is set in the config.
        :param list check_params_cbs:
            optional list of callables that compare the provided signal parameters to some expected parameters.
            Each callable has to match the signature of the corresponding signal in ``signals`` (just like a slot
            function would) and return ``True`` if parameters match, ``False`` otherwise.
            Instead of a specific callable, ``None`` can be provided, to disable parameter checking for the
            corresponding signal.
            If the number of callbacks doesn't match the number of signals ``ValueError`` will be raised.
        :param str order:
            Determines the order in which to expect signals:

            - ``"none"``: no order is enforced
            - ``"strict"``: signals have to be emitted strictly in the provided order
              (e.g. fails when expecting signals [a, b] and [a, a, b] is emitted)
            - ``"simple"``: like "strict", but signals may be emitted in-between the provided ones, e.g. expected
              ``signals == [a, b, c]`` and actually emitted ``signals = [a, a, b, a, c]`` works
              (would fail with ``"strict"``).

        :returns:
            ``MultiSignalBlocker`` object. Call ``MultiSignalBlocker.wait()``
            to wait.

        .. note:: This method is also available as ``wait_signals`` (pep-8 alias)
        """
        if order not in ["none", "simple", "strict"]:
            raise ValueError("order has to be set to 'none', 'simple' or 'strict'")

        if not signals:
            raise ValueError(
                f"Passing {signals} as signals isn't supported anymore, consider using qtbot.wait({timeout}) instead."
            )

        raising = self._should_raise(raising)

        if check_params_cbs:
            if len(check_params_cbs) != len(signals):
                raise ValueError(
                    "Number of callbacks ({}) does not "
                    "match number of signals ({})!".format(
                        len(check_params_cbs), len(signals)
                    )
                )
        blocker = MultiSignalBlocker(
            timeout=timeout,
            raising=raising,
            order=order,
            check_params_cbs=check_params_cbs,
        )
        blocker.add_signals(signals)
        return blocker

    def wait(self, ms: int) -> None:
        """
        .. versionadded:: 1.9

        Waits for ``ms`` milliseconds.

        While waiting, events will be processed and your test will stay
        responsive to user interface events or network communication.
        """
        blocker = MultiSignalBlocker(timeout=ms, raising=False)
        blocker.wait()

    @contextlib.contextmanager
    def assertNotEmitted(
        self, signal: SignalInstance, *, wait: int = 0
    ) -> Generator[None, None, None]:
        """
        .. versionadded:: 1.11

        Make sure the given ``signal`` doesn't get emitted.

        :param int wait:
            How many milliseconds to wait to make sure the signal isn't emitted
            asynchronously. By default, this method returns immediately and only
            catches signals emitted inside the ``with``-block.

        This is intended to be used as a context manager.

        .. note:: This method is also available as ``assert_not_emitted``
                  (pep-8 alias)
        """
        spy = SignalEmittedSpy(signal)
        with spy, self.waitSignal(signal, timeout=wait, raising=False):
            yield
        spy.assert_not_emitted()

    def waitUntil(
        self, callback: Callable[[], Optional[bool]], *, timeout: int = 5000
    ) -> None:
        """
        .. versionadded:: 2.0

        Wait in a busy loop, calling the given callback periodically until timeout is reached.

        ``callback()`` should raise ``AssertionError`` to indicate that the desired condition
        has not yet been reached, or just return ``None`` when it does. Useful to ``assert`` until
        some condition is satisfied:

        .. code-block:: python

            def view_updated():
                assert view_model.count() > 10


            qtbot.waitUntil(view_updated)

        Another possibility is for ``callback()`` to return ``True`` when the desired condition
        is met, ``False`` otherwise. Useful specially with ``lambda`` for terser code, but keep
        in mind that the error message in those cases is usually not very useful because it is
        not using an ``assert`` expression.

        .. code-block:: python

            qtbot.waitUntil(lambda: view_model.count() > 10)

        Note that this usage only accepts returning actual ``True`` and ``False`` values,
        so returning an empty list to express "falseness" raises a ``ValueError``.

        :param callback: callable that will be called periodically.
        :param timeout: timeout value in ms.
        :raises ValueError: if the return value from the callback is anything other than ``None``,
            ``True`` or ``False``.

        .. note:: This method is also available as ``wait_until`` (pep-8 alias)
        """
        __tracebackhide__ = True
        import time

        start = time.time()

        def timed_out():
            elapsed = time.time() - start
            elapsed_ms = elapsed * 1000
            return elapsed_ms > timeout

        timeout_msg = f"waitUntil timed out in {timeout} milliseconds"

        while True:
            try:
                result = callback()
            except AssertionError as e:
                if timed_out():
                    raise TimeoutError(timeout_msg) from e
            else:
                if result not in (None, True, False):
                    msg = "waitUntil() callback must return None, True or False, returned %r"
                    raise ValueError(msg % result)

                # 'assert' form
                if result is None:
                    return

                # 'True/False' form
                if result:
                    return
                if timed_out():
                    raise TimeoutError(timeout_msg)
            self.wait(10)

    def waitCallback(
        self, *, timeout: int = 5000, raising: Optional[bool] = None
    ) -> "CallbackBlocker":
        """
        .. versionadded:: 3.1

        Stops current test until a callback is called.

        Used to stop the control flow of a test until the returned callback is
        called, or a number of milliseconds, specified by ``timeout``, has
        elapsed.

        Best used as a context manager::

           with qtbot.waitCallback() as callback:
               function_taking_a_callback(callback)
           assert callback.args == [True]

        Also, you can use the :class:`CallbackBlocker` directly if the
        context manager form is not convenient::

           blocker = qtbot.waitCallback(timeout=1000)
           function_calling_a_callback(blocker)
           blocker.wait()


        :param int timeout:
            How many milliseconds to wait before resuming control flow.
        :param bool raising:
            If :class:`qtbot.TimeoutError <pytestqt.exceptions.TimeoutError>`
            should be raised if a timeout occurred.
            This defaults to ``True`` unless ``qt_default_raising = false``
            is set in the config.
        :returns:
            A ``CallbackBlocker`` object which can be used directly as a
            callback as it implements ``__call__``.

        .. note:: This method is also available as ``wait_callback`` (pep-8 alias)
        """
        raising = self._should_raise(raising)
        blocker = CallbackBlocker(timeout=timeout, raising=raising)
        return blocker

    @contextlib.contextmanager
    def captureExceptions(self) -> Iterator["CapturedExceptions"]:
        """
        .. versionadded:: 2.1

        Context manager that captures Qt virtual method exceptions that happen in block inside
        context.

        .. code-block:: python

            with qtbot.capture_exceptions() as exceptions:
                qtbot.mouseClick(
                    button, QtCore.MouseButton.LeftButton  # for PyQt6
                )  # ... or QtCore.LeftButton in PyQt5

            # exception is a list of sys.exc_info tuples
            assert len(exceptions) == 1

        .. note:: This method is also available as ``capture_exceptions`` (pep-8 alias)
        """
        from pytestqt.exceptions import capture_exceptions

        with capture_exceptions() as exceptions:
            yield exceptions

    def screenshot(
        self, widget: QWidget, suffix: str = "", region: Optional[QRect] = None
    ) -> Path:
        """
        .. versionadded:: 4.1

        Take a screenshot of the given widget and save it.

        The file is saved in a test-specific directory using pytest's ``tmp_path``
        fixture. The filename is ensured to be unique using a counter, and contains the
        ``objectName()`` of the widget if set, as well as its class name. A custom
        ``suffix`` can be given to add to the generated name.

        Raises :class:`qtbot.ScreenshotError <pytestqt.exceptions.ScreenshotError>`
        if taking the screenshot or saving the file failed.

        :param QWidget widget:
            The widget to take a screenshot of.
        :param str suffix:
            An optional suffix to add to the filename.
        :param QRect region:
            The region of the widget to screeshot. By default, the entire widget is
            contained.
        :returns:
            A ``pathlib.Path`` object with the taken screenshot.
        """
        pixmap = widget.grab() if region is None else widget.grab(region)
        if pixmap.isNull():
            raise ScreenshotError("Got null pixmap from Qt")

        tmp_path = self._request.getfixturevalue("tmp_path")

        parts = ["screenshot", widget.__class__.__name__]
        name = widget.objectName()
        if name:
            parts.append(name)
        if suffix:
            parts.append(suffix)

        for i in range(1, 500):
            counter = [] if i == 1 else [str(i)]

            path = tmp_path / ("_".join(parts + counter) + ".png")
            if path.exists():
                continue

            ok = pixmap.save(str(path))
            if not ok:
                raise ScreenshotError(f"Saving to {path} failed")

            return path

        raise ScreenshotError(f"Failed to find unique filename, last try: {path}")

    @staticmethod
    def keyClick(*args, **kwargs):
        qt_api.QtTest.QTest.keyClick(*args, **kwargs)

    @staticmethod
    def keyClicks(*args, **kwargs):
        qt_api.QtTest.QTest.keyClicks(*args, **kwargs)

    @staticmethod
    def keyEvent(*args, **kwargs):
        qt_api.QtTest.QTest.keyEvent(*args, **kwargs)

    @staticmethod
    def keyPress(*args, **kwargs):
        qt_api.QtTest.QTest.keyPress(*args, **kwargs)

    @staticmethod
    def keyRelease(*args, **kwargs):
        qt_api.QtTest.QTest.keyRelease(*args, **kwargs)

    @staticmethod
    def keySequence(widget: QWidget, key_sequence: QKeySequence) -> None:
        if not hasattr(qt_api.QtTest.QTest, "keySequence"):
            raise NotImplementedError("This method is available from Qt 5.10 upwards.")
        qt_api.QtTest.QTest.keySequence(widget, key_sequence)

    @staticmethod
    def keyToAscii(key: Any) -> None:
        if not hasattr(qt_api.QtTest.QTest, "keyToAscii"):
            raise NotImplementedError("This method isn't available on PyQt5.")
        qt_api.QtTest.QTest.keyToAscii(key)

    @staticmethod
    def mouseClick(*args, **kwargs):
        qt_api.QtTest.QTest.mouseClick(*args, **kwargs)

    @staticmethod
    def mouseDClick(*args, **kwargs):
        qt_api.QtTest.QTest.mouseDClick(*args, **kwargs)

    @staticmethod
    def mouseMove(*args, **kwargs):
        qt_api.QtTest.QTest.mouseMove(*args, **kwargs)

    @staticmethod
    def mousePress(*args, **kwargs):
        qt_api.QtTest.QTest.mousePress(*args, **kwargs)

    @staticmethod
    def mouseRelease(*args, **kwargs):
        qt_api.QtTest.QTest.mouseRelease(*args, **kwargs)

    # provide easy access to exceptions to qtbot fixtures
    SignalEmittedError = SignalEmittedError
    TimeoutError = TimeoutError
    ScreenshotError = ScreenshotError
    CallbackCalledTwiceError = CallbackCalledTwiceError


def _add_widget(
    item: pytest.Item,
    widget: QWidget,
    *,
    before_close_func: Optional[BeforeCloseFunc] = None,
) -> None:
    """
    Register a widget into the given pytest item for later closing.
    """
    qt_widgets = getattr(item, "qt_widgets", [])
    qt_widgets.append((weakref.ref(widget), before_close_func))
    item.qt_widgets = qt_widgets  # type: ignore[attr-defined]


def _close_widgets(item: pytest.Item) -> None:
    """
    Close all widgets registered in the pytest item.
    """
    widgets = getattr(item, "qt_widgets", None)
    if widgets:
        for w, before_close_func in item.qt_widgets:  # type: ignore[attr-defined]
            w = w()
            if w is not None:
                if before_close_func is not None:
                    before_close_func(w)
                w.close()
                w.deleteLater()
        del item.qt_widgets  # type: ignore[attr-defined]


def _iter_widgets(item: pytest.Item) -> Iterator[weakref.ReferenceType[QWidget]]:
    """
    Iterates over widgets registered in the given pytest item.
    """
    qt_widgets = getattr(item, "qt_widgets", [])
    return (w for (w, _) in qt_widgets)


WaitAdjectiveName = Literal["activated", "exposed"]


class _WaitWidgetContextManager:
    """
    Context manager implementation used by ``waitActive`` and ``waitExposed`` methods.
    """

    def __init__(
        self,
        method_name: str,
        adjective_name: WaitAdjectiveName,
        widget: QWidget,
        timeout: int,
    ) -> None:
        """
        :param str method_name: name to the ``QtTest`` method to call to check if widget is active/exposed.
        :param str adjective_name: "activated" or "exposed".
        :param widget:
        :param timeout:
        """
        self._method_name = method_name
        self._adjective_name = adjective_name
        self._widget = widget
        self._timeout = timeout

    def __enter__(self) -> Self:
        __tracebackhide__ = True
        return self

    def __exit__(
        self,
        exc_type: Optional[type[BaseException]],
        exc_val: Optional[BaseException],
        exc_tb: Optional[TracebackType],
    ) -> None:
        __tracebackhide__ = True
        try:
            if exc_type is None:
                method = getattr(qt_api.QtTest.QTest, self._method_name)
                r = method(self._widget, self._timeout)
                if not r:
                    msg = "widget {} not {} in {} ms.".format(
                        self._widget, self._adjective_name, self._timeout
                    )
                    raise TimeoutError(msg)
        finally:
            self._widget = None
