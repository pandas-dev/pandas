"""
========================
Widget testing utilities
========================

See also :mod:`matplotlib.tests.test_widgets`.
"""

from unittest import mock

from matplotlib import _api
from matplotlib.backend_bases import MouseEvent, KeyEvent
import matplotlib.pyplot as plt


def get_ax():
    """Create a plot and return its Axes."""
    fig, ax = plt.subplots(1, 1)
    ax.plot([0, 200], [0, 200])
    ax.set_aspect(1.0)
    fig.canvas.draw()
    return ax


def noop(*args, **kwargs):
    pass


@_api.deprecated("3.11", alternative="MouseEvent or KeyEvent")
def mock_event(ax, button=1, xdata=0, ydata=0, key=None, step=1):
    r"""
    Create a mock event that can stand in for `.Event` and its subclasses.

    This event is intended to be used in tests where it can be passed into
    event handling functions.

    Parameters
    ----------
    ax : `~matplotlib.axes.Axes`
        The Axes the event will be in.
    xdata : float
        x coord of mouse in data coords.
    ydata : float
        y coord of mouse in data coords.
    button : None or `MouseButton` or {'up', 'down'}
        The mouse button pressed in this event (see also `.MouseEvent`).
    key : None or str
        The key pressed when the mouse event triggered (see also `.KeyEvent`).
    step : int
        Number of scroll steps (positive for 'up', negative for 'down').

    Returns
    -------
    event
        A `.Event`\-like Mock instance.
    """
    event = mock.Mock()
    event.button = button
    event.x, event.y = ax.transData.transform([(xdata, ydata),
                                               (xdata, ydata)])[0]
    event.xdata, event.ydata = xdata, ydata
    event.inaxes = ax
    event.canvas = ax.get_figure(root=True).canvas
    event.key = key
    event.step = step
    event.guiEvent = None
    event.name = 'Custom'
    return event


@_api.deprecated("3.11", alternative="callbacks.process(event)")
def do_event(tool, etype, button=1, xdata=0, ydata=0, key=None, step=1):
    """
    Trigger an event on the given tool.

    Parameters
    ----------
    tool : matplotlib.widgets.AxesWidget
    etype : str
        The event to trigger.
    xdata : float
        x coord of mouse in data coords.
    ydata : float
        y coord of mouse in data coords.
    button : None or `MouseButton` or {'up', 'down'}
        The mouse button pressed in this event (see also `.MouseEvent`).
    key : None or str
        The key pressed when the mouse event triggered (see also `.KeyEvent`).
    step : int
        Number of scroll steps (positive for 'up', negative for 'down').
    """
    event = mock_event(tool.ax, button, xdata, ydata, key, step)
    func = getattr(tool, etype)
    func(event)


def click_and_drag(tool, start, end, key=None):
    """
    Helper to simulate a mouse drag operation.

    Parameters
    ----------
    tool : `~matplotlib.widgets.Widget`
    start : [float, float]
        Starting point in data coordinates.
    end : [float, float]
        End point in data coordinates.
    key : None or str
         An optional key that is pressed during the whole operation
         (see also `.KeyEvent`).
    """
    ax = tool.ax
    if key is not None:  # Press key
        KeyEvent._from_ax_coords("key_press_event", ax, start, key)._process()
    # Click, move, and release mouse
    MouseEvent._from_ax_coords("button_press_event", ax, start, 1)._process()
    MouseEvent._from_ax_coords("motion_notify_event", ax, end, 1)._process()
    MouseEvent._from_ax_coords("button_release_event", ax, end, 1)._process()
    if key is not None:  # Release key
        KeyEvent._from_ax_coords("key_release_event", ax, end, key)._process()
