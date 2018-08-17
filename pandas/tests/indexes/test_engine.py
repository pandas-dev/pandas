# -*- coding: utf-8 -*-

import numpy as np
import pytest

import pandas as pd
from pandas._libs.index import (Int64Engine, UInt64Engine,
                                Float64Engine, ObjectEngine)


class TestNumericEngine(object):

    def setup_class(cls):
        cls.data = [1, 2, 3]

    def test_engine_type(self, numeric_indexing_engine):
        index = pd.Index(self.data, dtype=numeric_indexing_engine._dtype)
        if issubclass(index.dtype.type, np.signedinteger):
            assert isinstance(index._engine, Int64Engine)
        elif issubclass(index.dtype.type, np.unsignedinteger):
            assert isinstance(index._engine, UInt64Engine)
        elif issubclass(index.dtype.type, np.floating):
            assert isinstance(index._engine, Float64Engine)
        else:
            raise TypeError("unexpected dtype {}".format(index.dtype))

    def test_is_monotonic_ordered(self, numeric_indexing_engine):
        codes = np.array(self.data, dtype=numeric_indexing_engine._dtype)
        e = numeric_indexing_engine(lambda: codes, len(codes))

        assert e.is_monotonic_increasing
        assert not e.is_monotonic_decreasing

        # reverse sort order
        reversed_data = list(reversed(self.data))
        codes = np.array(reversed_data, dtype=numeric_indexing_engine._dtype)
        e = numeric_indexing_engine(lambda: codes, len(codes))

        assert not e.is_monotonic_increasing
        assert e.is_monotonic_decreasing

    def test_is_not_monotonic_ordered(self, numeric_indexing_engine):
        data = [1, 0, 2]
        codes = np.array(data, dtype=numeric_indexing_engine._dtype)
        e = numeric_indexing_engine(lambda: codes, len(codes))

        assert not e.is_monotonic_increasing
        assert not e.is_monotonic_decreasing

    @pytest.mark.parametrize('values, expected', [
        ([1, 2, 3], True),
        ([1, 1, 2], False),
    ])
    def test_is_unique(self, values, expected, numeric_indexing_engine):
        codes = np.array(values, dtype=numeric_indexing_engine._dtype)
        e = numeric_indexing_engine(lambda: codes, len(codes))

        assert e.is_unique is expected

    @pytest.mark.parametrize('values, value, expected', [
        ([1, 2, 3], 2, 1),
        ([1, 2, 2, 3], 2, slice(1, 3)),
        ([3, 2, 2, 1], 2, np.array([False, True, True, False])),
        ([1, 2, 2, 1], 2, np.array([False, True, True, False])),
        ([1, 3, 2], 2, 2),
    ])
    def test_get_loc(self, values, value, expected, numeric_indexing_engine):
        codes = np.array(values, dtype=numeric_indexing_engine._dtype)
        e = numeric_indexing_engine(lambda: codes, len(codes))
        result = e.get_loc(value)

        if isinstance(expected, np.ndarray):
            assert (result == expected).all()
        else:
            assert result == expected

    @pytest.mark.parametrize('values, value, error', [
        ([1, 2, 3], 4, KeyError),
        ([1, 2, 3], '4', KeyError),
    ])
    def test_get_loc_raises(self, values, value, error,
                            numeric_indexing_engine):
        codes = np.array(values, dtype=numeric_indexing_engine._dtype)
        e = numeric_indexing_engine(lambda: codes, len(codes))
        with pytest.raises(error):
            e.get_loc(value)


class TestObjectEngine(object):

    def setup_class(cls):
        cls.data = list('abc')
        cls.dtype = object
        cls.Engine = ObjectEngine

    def test_engine_type(self):
        index = pd.Index(self.data)
        assert isinstance(index._engine, self.Engine)

    def test_is_monotonic_ordered(self):
        codes = np.array(self.data, dtype=self.dtype)
        e = self.Engine(lambda: codes, len(codes))

        assert e.is_monotonic_increasing
        assert not e.is_monotonic_decreasing

        # reverse sort order
        reversed_data = list(reversed(self.data))
        codes = np.array(reversed_data, dtype=self.dtype)
        e = self.Engine(lambda: codes, len(codes))

        assert not e.is_monotonic_increasing
        assert e.is_monotonic_decreasing

    def test_is_not_monotonic_ordered(self):
        codes = np.array(list('cab'), dtype=self.dtype)
        e = self.Engine(lambda: codes, len(codes))

        assert not e.is_monotonic_increasing
        assert not e.is_monotonic_decreasing

    @pytest.mark.parametrize('values, expected', [
        (list('abc'), True),
        (list('aab'), False),
    ])
    def test_is_unique(self, values, expected):
        codes = np.array(values, dtype=self.dtype)
        e = self.Engine(lambda: codes, len(codes))

        assert e.is_unique is expected

    @pytest.mark.parametrize('values, value, expected', [
        (list('abc'), 'b', 1),
        (list('abbc'), 'b', slice(1, 3)),
        (list('cbba'), 'b', np.array([False, True, True, False])),
        (list('abba'), 'b', np.array([False, True, True, False])),
        (list('acb'), 'b', 2),
    ])
    def test_get_loc(self, values, value, expected):
        codes = np.array(values, dtype=self.dtype)
        e = self.Engine(lambda: codes, len(codes))
        result = e.get_loc(value)

        if isinstance(expected, np.ndarray):
            assert (result == expected).all()
        else:
            assert result == expected

    @pytest.mark.parametrize('values, value, error', [
        (list('abc'), 'd', KeyError),
        (list('abc'), 4, KeyError),
    ])
    def test_get_loc_raises(self, values, value, error):
        codes = np.array(values, dtype=self.dtype)
        e = self.Engine(lambda: codes, len(codes))
        with pytest.raises(error):
            e.get_loc(value)
