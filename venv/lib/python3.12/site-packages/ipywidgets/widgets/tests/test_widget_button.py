# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.

import inspect
import pytest
from ipywidgets import Button

def test_deprecation_fa_icons():
    with pytest.deprecated_call() as record:
        Button(icon='fa-home')
    assert len(record) == 1
    assert record[0].filename == inspect.stack(context=0)[1].filename
