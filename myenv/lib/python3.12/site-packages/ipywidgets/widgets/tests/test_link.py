# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.

import pytest

from .. import jslink, jsdlink, ToggleButton
from .utils import setup, teardown

def test_jslink_args():
    with pytest.raises(TypeError):
        jslink()
    w1 = ToggleButton()
    with pytest.raises(TypeError):
        jslink((w1, 'value'))

    w2 = ToggleButton()
    jslink((w1, 'value'), (w2, 'value'))

    with pytest.raises(TypeError):
        jslink((w1, 'value'), (w2, 'nosuchtrait'))

    with pytest.raises(TypeError):
        jslink((w1, 'value'), (w2, 'traits'))

def test_jsdlink_args():
    with pytest.raises(TypeError):
        jsdlink()
    w1 = ToggleButton()
    with pytest.raises(TypeError):
        jsdlink((w1, 'value'))

    w2 = ToggleButton()
    jsdlink((w1, 'value'), (w2, 'value'))

    with pytest.raises(TypeError):
        jsdlink((w1, 'value'), (w2, 'nosuchtrait'))

    with pytest.raises(TypeError):
        jsdlink((w1, 'value'), (w2, 'traits'))
