# -*- coding: utf-8 -*-


import numpy as np
import pandas.util.testing as tm
import pytest
from pandas import MultiIndex
from pandas.compat import PY3, long


def test_numeric_compat(idx):
    tm.assert_raises_regex(TypeError, "cannot perform __mul__",
                           lambda: idx * 1)
    tm.assert_raises_regex(TypeError, "cannot perform __rmul__",
                           lambda: 1 * idx)

    div_err = "cannot perform __truediv__" if PY3 \
        else "cannot perform __div__"
    tm.assert_raises_regex(TypeError, div_err, lambda: idx / 1)
    div_err = div_err.replace(' __', ' __r')
    tm.assert_raises_regex(TypeError, div_err, lambda: 1 / idx)
    tm.assert_raises_regex(TypeError, "cannot perform __floordiv__",
                           lambda: idx // 1)
    tm.assert_raises_regex(TypeError, "cannot perform __rfloordiv__",
                           lambda: 1 // idx)


def test_logical_compat(idx):
    tm.assert_raises_regex(TypeError, 'cannot perform all',
                           lambda: idx.all())
    tm.assert_raises_regex(TypeError, 'cannot perform any',
                           lambda: idx.any())


def test_boolean_context_compat(idx):

    with pytest.raises(ValueError):
        bool(idx)


def test_boolean_context_compat2():

    # boolean context compat
    # GH7897
    i1 = MultiIndex.from_tuples([('A', 1), ('A', 2)])
    i2 = MultiIndex.from_tuples([('A', 1), ('A', 3)])
    common = i1.intersection(i2)

    with pytest.raises(ValueError):
        bool(common)


def test_inplace_mutation_resets_values():
    levels = [['a', 'b', 'c'], [4]]
    levels2 = [[1, 2, 3], ['a']]
    labels = [[0, 1, 0, 2, 2, 0], [0, 0, 0, 0, 0, 0]]

    mi1 = MultiIndex(levels=levels, labels=labels)
    mi2 = MultiIndex(levels=levels2, labels=labels)
    vals = mi1.values.copy()
    vals2 = mi2.values.copy()

    assert mi1._tuples is not None

    # Make sure level setting works
    new_vals = mi1.set_levels(levels2).values
    tm.assert_almost_equal(vals2, new_vals)

    # Non-inplace doesn't kill _tuples [implementation detail]
    tm.assert_almost_equal(mi1._tuples, vals)

    # ...and values is still same too
    tm.assert_almost_equal(mi1.values, vals)

    # Inplace should kill _tuples
    mi1.set_levels(levels2, inplace=True)
    tm.assert_almost_equal(mi1.values, vals2)

    # Make sure label setting works too
    labels2 = [[0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0]]
    exp_values = np.empty((6,), dtype=object)
    exp_values[:] = [(long(1), 'a')] * 6

    # Must be 1d array of tuples
    assert exp_values.shape == (6,)
    new_values = mi2.set_labels(labels2).values

    # Not inplace shouldn't change
    tm.assert_almost_equal(mi2._tuples, vals2)

    # Should have correct values
    tm.assert_almost_equal(exp_values, new_values)

    # ...and again setting inplace should kill _tuples, etc
    mi2.set_labels(labels2, inplace=True)
    tm.assert_almost_equal(mi2.values, new_values)


def test_ndarray_compat_properties(idx, compat_props):
    assert idx.T.equals(idx)
    assert idx.transpose().equals(idx)

    values = idx.values
    for prop in compat_props:
        assert getattr(idx, prop) == getattr(values, prop)

    # test for validity
    idx.nbytes
    idx.values.nbytes


def test_compat(indices):
    assert indices.tolist() == list(indices)


def test_pickle_compat_construction(holder):
    # this is testing for pickle compat
    if holder is None:
        return

    # need an object to create with
    pytest.raises(TypeError, holder)
