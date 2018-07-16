# -*- coding: utf-8 -*-

import numpy as np
import pytest
from pandas import Index, MultiIndex


@pytest.fixture
def idx():
    # a MultiIndex used to test the general functionality of the
    # general functionality of this object
    major_axis = Index(['foo', 'bar', 'baz', 'qux'])
    minor_axis = Index(['one', 'two'])

    major_labels = np.array([0, 0, 1, 2, 3, 3])
    minor_labels = np.array([0, 1, 0, 1, 0, 1])
    index_names = ['first', 'second']
    mi = MultiIndex(levels=[major_axis, minor_axis],
                    labels=[major_labels, minor_labels],
                    names=index_names, verify_integrity=False)
    return mi


@pytest.fixture
def idx_dup():
    # compare tests/indexes/multi/conftest.py
    major_axis = Index(['foo', 'bar', 'baz', 'qux'])
    minor_axis = Index(['one', 'two'])

    major_labels = np.array([0, 0, 1, 0, 1, 1])
    minor_labels = np.array([0, 1, 0, 1, 0, 1])
    index_names = ['first', 'second']
    mi = MultiIndex(levels=[major_axis, minor_axis],
                    labels=[major_labels, minor_labels],
                    names=index_names, verify_integrity=False)
    return mi


@pytest.fixture
def index_names():
    # names that match those in the idx fixture for testing equality of
    # names assigned to the idx
    return ['first', 'second']


@pytest.fixture
def holder():
    # the MultiIndex constructor used to base compatibility with pickle
    return MultiIndex


@pytest.fixture
def compat_props():
    # a MultiIndex must have these properties associated with it
    return ['shape', 'ndim', 'size']
