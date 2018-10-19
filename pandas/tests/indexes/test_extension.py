# -*- coding: utf-8 -*-

import pytest

import pandas.util.testing as tm
from pandas.core.indexes.api import Index
from .common import Base

import numpy as np

from pandas.util.testing import (
    assert_extension_array_equal, assert_index_equal)

from pandas.core.arrays import integer_array
from pandas.core.indexes.extension import ExtensionIndex


@pytest.fixture
def data():
    return integer_array([1, 2, 3, 4])


def test_constructor(data):
    result = ExtensionIndex(data, name='test')
    assert result.name == 'test'
    assert isinstance(result, ExtensionIndex)
    assert_extension_array_equal(data, result._values)

    expected = ExtensionIndex(data, name='test')
    # data and passed dtype match
    result = ExtensionIndex(data, dtype=data.dtype, name='test')
    assert_index_equal(result, expected)
    # data is converted to passed dtype
    result = ExtensionIndex(np.array(data), dtype=data.dtype, name='test')
    assert_index_equal(result, expected)
    # EA is converted to passed dtype
    expected = ExtensionIndex(integer_array(data, dtype='Int32'), name='test')
    result = ExtensionIndex(data, dtype=expected.dtype, name='test')
    assert_index_equal(result, expected)

    # no ExtensionDtype passed
    with pytest.raises(ValueError):
        ExtensionIndex(data, dtype='int64', name='test')

    with pytest.raises(ValueError):
        ExtensionIndex(data, dtype=object, name='test')

    # no ExtensionArray passed
    with pytest.raises(ValueError):
        ExtensionIndex(np.array(data), name='test')


def test_default_index_constructor(data):
    result = Index(data, name='test')
    expected = ExtensionIndex(data, name='test')
    assert_index_equal(result, expected)

    result = Index(data, dtype=data.dtype, name='test')
    assert_index_equal(result, expected)

    result = Index(np.array(data), dtype=data.dtype, name='test')
    assert_index_equal(result, expected)

    result = Index(data, dtype=object, name='test')
    expected = Index(np.array(data), dtype=object, name='test')
    assert_index_equal(result, expected)


# class TestExtensionIndex(Base):
#     _holder = ExtensionIndex

#     def setup_method(self, method):
#         self.indices = dict(
#             extIndex=ExtensionIndex(np.arange(100), dtype='Int64'))
#         self.setup_indices()

#     # def create_index(self):
#     #     if categories is None:
#     #         categories = list('cab')
#     #     return CategoricalIndex(
#     #         list('aabbca'), categories=categories, ordered=ordered)
