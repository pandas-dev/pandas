# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.

"""Test Widget."""

import inspect

import pytest
from IPython.core.interactiveshell import InteractiveShell
from IPython.display import display
from IPython.utils.capture import capture_output

from .. import widget
from ..widget import Widget
from ..widget_button import Button
import copy


def test_no_widget_view():
    # ensure IPython shell is instantiated
    # otherwise display() just calls print
    shell = InteractiveShell.instance()

    with capture_output() as cap:
        w = Widget()
        display(w)

    assert len(cap.outputs) == 1, "expect 1 output"
    mime_bundle = cap.outputs[0].data
    assert mime_bundle["text/plain"] == repr(w), "expected plain text output"
    assert (
        "application/vnd.jupyter.widget-view+json" not in mime_bundle
    ), "widget has no view"
    assert cap.stdout == "", repr(cap.stdout)
    assert cap.stderr == "", repr(cap.stderr)


def test_widget_view():
    # ensure IPython shell is instantiated
    # otherwise display() just calls print
    shell = InteractiveShell.instance()

    with capture_output() as cap:
        w = Button()
        display(w)

    assert len(cap.outputs) == 1, "expect 1 output"
    mime_bundle = cap.outputs[0].data
    assert mime_bundle["text/plain"] == repr(w), "expected plain text output"
    assert (
        "application/vnd.jupyter.widget-view+json" in mime_bundle
    ), "widget should have have a view"
    assert cap.stdout == "", repr(cap.stdout)
    assert cap.stderr == "", repr(cap.stderr)


def test_close_all():
    # create a couple of widgets
    widgets = [Button() for i in range(10)]

    assert len(widget._instances) > 0, "expect active widgets"
    assert widget._instances[widgets[0].model_id] is widgets[0]
    # close all the widgets
    Widget.close_all()

    assert len(widget._instances) == 0, "active widgets should be cleared"


def test_compatibility():
    button = Button()
    assert widget._instances[button.model_id] is button
    with pytest.deprecated_call() as record:
        assert widget._instances is widget.Widget.widgets
        assert widget._instances is widget.Widget._active_widgets
        assert widget._registry is widget.Widget.widget_types
        assert widget._registry is widget.Widget._widget_types

        Widget.close_all()
        assert not widget.Widget.widgets
        assert not widget.Widget._active_widgets
    caller_path = inspect.stack(context=0)[1].filename
    assert all(x.filename == caller_path for x in record)
    assert len(record) == 6


def test_widget_copy():
    button = Button()
    with pytest.raises(NotImplementedError):
        copy.copy(button)
    with pytest.raises(NotImplementedError):
        copy.deepcopy(button)