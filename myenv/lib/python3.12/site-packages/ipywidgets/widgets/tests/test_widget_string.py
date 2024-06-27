# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.

import inspect
import pytest

from ..widget_string import Combobox, Text

def test_combobox_creation_blank():
    w = Combobox()
    assert w.value == ''
    assert w.options == ()
    assert w.ensure_option == False


def test_combobox_creation_kwargs():
    w = Combobox(
        value='Chocolate',
        options=[
            "Chocolate",
            "Coconut",
            "Mint",
            "Strawberry",
            "Vanilla",
        ],
        ensure_option=True
    )
    assert w.value == 'Chocolate'
    assert w.options == (
            "Chocolate",
            "Coconut",
            "Mint",
            "Strawberry",
            "Vanilla",
        )
    assert w.ensure_option == True

def test_tooltip_deprecation():
    caller_path = inspect.stack(context=0)[1].filename
    with pytest.deprecated_call() as record:
        w = Text(description_tooltip="testing")
    assert len(record) == 1
    assert record[0].filename == caller_path

    with pytest.deprecated_call() as record:
        w.description_tooltip
    assert len(record) == 1
    assert record[0].filename == caller_path

    with pytest.deprecated_call() as record:
        w.description_tooltip == "testing"
    assert len(record) == 1
    assert record[0].filename == caller_path

    with pytest.deprecated_call() as record:
        w.description_tooltip = "second value"
    assert len(record) == 1
    assert record[0].filename == caller_path
    assert w.tooltip == "second value"

def test_on_submit_deprecation():
    with pytest.deprecated_call() as record:
        Text().on_submit(lambda *args: ...)
    assert len(record) == 1
    assert record[0].filename == inspect.stack(context=0)[1].filename
