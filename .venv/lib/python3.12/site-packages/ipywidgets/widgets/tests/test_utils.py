# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.

import inspect
import pytest

from ..utils import deprecation
from .utils import call_method

CALL_PATH = inspect.getfile(call_method)

def test_deprecation():
    caller_path = inspect.stack(context=0)[1].filename
    with pytest.deprecated_call() as record:
        deprecation('Deprecated call')
    # Make sure the deprecation pointed to the external function calling this test function
    assert len(record) == 1
    assert record[0].filename == caller_path

    with pytest.deprecated_call() as record:
        deprecation('Deprecated call', ['ipywidgets/widgets/tests'])
    # Make sure the deprecation pointed to the external function calling this test function
    assert len(record) == 1
    assert record[0].filename == caller_path

    with pytest.deprecated_call() as record:
        deprecation('Deprecated call', 'ipywidgets/widgets/tests')
    # Make sure the deprecation pointed to the external function calling this test function
    assert len(record) == 1
    assert record[0].filename == caller_path

    with pytest.deprecated_call() as record:
        deprecation('Deprecated call', [])
    # Make sure the deprecation pointed to *this* file
    assert len(record) == 1
    assert record[0].filename == __file__

def test_deprecation_indirect():
    # If the line that calls "deprecation" is not internal, it is considered the source:
    with pytest.warns(DeprecationWarning) as record:
        call_method(deprecation, "test message", [])
    assert len(record) == 1
    assert record[0].filename == CALL_PATH

def test_deprecation_indirect_internal():
    # If the line that calls "deprecation" is internal, it is not considered the source:
    with pytest.warns(DeprecationWarning) as record:
        call_method(deprecation, "test message", [CALL_PATH])
    assert len(record) == 1
    assert record[0].filename == __file__

def test_deprecation_nested1():
    def level1():
        deprecation("test message", [])

    with pytest.warns(DeprecationWarning) as record:
        call_method(level1)

    assert len(record) == 1
    assert record[0].filename == __file__

def test_deprecation_nested2():
    def level2():
        deprecation("test message", [])
    def level1():
        level2()

    with pytest.warns(DeprecationWarning) as record:
        call_method(level1)

    assert len(record) == 1
    assert record[0].filename == __file__
