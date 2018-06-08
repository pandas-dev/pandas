# -*- coding: utf-8 -*-

import numpy as np
import pytest
from pandas import Index, MultiIndex


@pytest.fixture
def _index():
    major_axis = Index(['foo', 'bar', 'baz', 'qux'])
    minor_axis = Index(['one', 'two'])

    major_labels = np.array([0, 0, 1, 2, 3, 3])
    minor_labels = np.array([0, 1, 0, 1, 0, 1])
    index_names = ['first', 'second']
    idx = MultiIndex(
        levels=[major_axis, minor_axis],
        labels=[major_labels, minor_labels],
        names=index_names,
        verify_integrity=False
    )
    return idx


@pytest.fixture
def named_index(_index):
    return {'index': _index}


@pytest.fixture
def index_names():
    return ['first', 'second']


@pytest.fixture
def _holder():
    return MultiIndex


@pytest.fixture
def _compat_props():
    return ['shape', 'ndim', 'size']
