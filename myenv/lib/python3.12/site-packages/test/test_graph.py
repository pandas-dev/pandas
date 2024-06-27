# Licensed under a 3-clause BSD style license - see LICENSE.rst

import os

from asv.graph import (Graph, RESAMPLED_POINTS, make_summary_graph, _fill_missing_data,
                       _combine_graph_data)
from asv import util


def test_graph_single():
    vals = [
        (1, 1, None),
        (2, 2, 1),
        (3, 3, None),
        (4, 4, 7),
        (5, 5, None),
        (6, None, None),
        (7, float('nan'), None),
    ]

    # Should give same data back, excluding missing values at edges
    g = Graph('foo', {})
    for k, v, dv in vals:
        g.add_data_point(k, v, dv)
    data = g.get_data()
    assert data == vals[:-2]

    # Should average duplicate values
    g = Graph('foo', {})
    g.add_data_point(4, 3)
    for k, v, dv in vals:
        g.add_data_point(k, v, dv)
    g.add_data_point(4, 5)
    data = g.get_data()
    assert data[3][0] == 4
    assert abs(data[3][1] - (3 + 4 + 5) / 3.) < 1e-10

    # Summary graph should be the same as the main graph
    g = Graph('foo', {})
    for k, v, dv in vals:
        g.add_data_point(k, v, dv)
    g = make_summary_graph([g])
    data = g.get_data()
    assert len(data) == len(vals) - 2
    for v, d in zip(vals, data):
        kv, xv, dv = v
        kd, xd, dd = d
        assert kv == kd
        assert dd is None
        assert abs(xv - xd) < 1e-10


def test_graph_multi():
    vals = [
        (0, [None, None, None], [None] * 3),
        (1, [1, None, float('nan')], [None] * 3),
        (2, [2, 5, 4], [1] * 3),
        (3, [3, 4, -60], [None] * 3),
        (4, [4, 3, 2], [None] * 3),
        (5, [None, 2, None], [None] * 3),
        (6, [6, 1, None], [None] * 3)
    ]

    filled_vals = [
        (1, [1, None, None]),
        (2, [2, 5, 4]),
        (3, [3, 4, -60]),
        (4, [4, 3, 2]),
        (5, [5, 2, None]),
        (6, [6, 1, None])
    ]

    # Should give same data back, with missing data at edges removed
    g = Graph('foo', {})
    for k, v, dv in vals:
        g.add_data_point(k, v, dv)
    data = g.get_data()
    assert data[0] == (1, [1, None, None], [None, None, None])
    assert data[1:] == vals[2:]

    # Should average duplicate values
    g = Graph('foo', {})
    g.add_data_point(4, [1, 2, 3])
    for k, v, dv in vals:
        g.add_data_point(k, v, dv)
    g.add_data_point(4, [3, 2, 1])
    data = g.get_data()
    assert data[3][0] == 4
    assert abs(data[3][1][0] - (1 + 4 + 3) / 3.) < 1e-10
    assert abs(data[3][1][1] - (2 + 3 + 2) / 3.) < 1e-10
    assert abs(data[3][1][2] - (3 + 2 + 1) / 3.) < 1e-10

    # The summary graph is obtained by geometric mean of filled data
    g = Graph('foo', {})
    for k, v, dv in vals:
        g.add_data_point(k, v, dv)
    g = make_summary_graph([g])
    data = g.get_data()

    for v, d in zip(filled_vals, data):
        kv, xvs = v
        kd, xd, dd = d
        assert kv == kd

        # geom mean, with some sign convention
        expected = util.geom_mean_na(xvs)
        assert (xd is None and expected is None) or abs(xd - expected) < 1e-10

    # Test summary over separate graphs -- should behave as if the
    # data was in a single graph
    g0 = Graph('foo', {})
    g1 = Graph('foo', {})
    g2 = Graph('foo', {})
    for k, v, dv in vals:
        g0.add_data_point(k, v, dv)
        g1.add_data_point(k, v[0], dv[0])
        g2.add_data_point(k, v[1:], dv[1:])

    data0 = make_summary_graph([g0]).get_data()
    data = make_summary_graph([g1, g2]).get_data()
    assert data == data0

    # Check the above is true regardless if some x-values are missing
    g0.add_data_point(7, [None, 1, None])
    g2.add_data_point(7, [1, None])

    g0.add_data_point(4.5, [9, None, None])
    g1.add_data_point(4.5, 9)

    data0 = make_summary_graph([g0]).get_data()
    data = make_summary_graph([g1, g2]).get_data()
    assert data == data0


def test_empty_graph():
    g = Graph('foo', {})
    g.add_data_point(1, None)
    g.add_data_point(2, None)
    g.add_data_point(3, None)
    data = g.get_data()
    assert data == []

    g = Graph('foo', {})
    g.add_data_point(1, None)
    g.add_data_point(1, [None, None])
    g.add_data_point(2, [None, None])
    g.add_data_point(3, None)
    g.add_data_point(4, [None, None])
    data = g.get_data()
    assert data == []


def test_nan():
    g = Graph('foo', {})
    g.add_data_point(1, 1)
    g.add_data_point(2, 2)
    g.add_data_point(2, float('nan'))
    g.add_data_point(3, 3)
    g.add_data_point(4, float('nan'))
    data = g.get_data()
    assert data == [(1, 1, None), (2, 2, None), (3, 3, None)]

    g = Graph('foo', {})
    g.add_data_point(1, None)
    g.add_data_point(1, [1, float('nan')])
    g.add_data_point(2, [2, 2])
    g.add_data_point(3, [float('nan'), float('nan')])
    g.add_data_point(4, [None, float('nan')])
    data = g.get_data()
    assert data == [(1, [1, None], [None, None]), (2, [2, 2], [None, None])]


def test_summary_graph():
    n = 2 * int(RESAMPLED_POINTS)
    g = Graph('foo', {})
    for i in range(n):
        g.add_data_point(i, 0.1)
        g.add_data_point(n + i, 0.2)
    g = make_summary_graph([g])
    data = g.get_data()
    assert len(data) == 512
    for i in range(256):
        assert abs(data[i][1] - 0.1) < 1e-7
        assert abs(data[256 + i][1] - 0.2) < 1e-7


def test_summary_graph_loop():
    n = int(RESAMPLED_POINTS)

    # Resampling shouldn't get stuck in an infinite loop
    g = Graph('foo', {})
    for j in range(n):
        g.add_data_point(j, 0.1)
    g = make_summary_graph([g])
    data = g.get_data()
    assert len(data) == 1
    assert data[0][0] == n
    assert abs(data[0][1] - 0.1) < 1e-7


def test__fill_missing_data():
    y = [None, 1, 2, None, None, 5, None, 7, None]
    filled = _fill_missing_data(y, max_gap_fraction=0.5)
    assert filled == [None, 1, 2, 3, 4, 5, 6, 7, None]

    filled = _fill_missing_data(y, max_gap_fraction=0.1)
    assert filled == [None, 1, 2, None, None, 5, 6, 7, None]


def test__combine_graph_data():
    g1 = Graph('1', {})
    g2 = Graph('2', {})
    g3 = Graph('3', {})
    for i in range(5):
        if i >= 3:
            g1.add_data_point(i, None)
        g2.add_data_point(i, i)
        if i <= 2:
            g3.add_data_point(i, [-i, -2 * i])

    x, ys = _combine_graph_data([g1, g2, g3])
    assert x == [0, 1, 2, 3, 4]
    assert ys == [[None, None, None, None, None],
                  [0, 1, 2, 3, 4],
                  [0, -1, -2, None, None],
                  [0, -2, -4, None, None]]


def test_graph_steps():
    vals = [(1, 1), (5, 1), (6, 1), (7, 1), (8, 1),
            (11, 2), (15, 2), (16, 2 + 1e-5), (17, 2), (18, 2)]

    g = Graph('foo', {})
    for x, y in vals:
        g.add_data_point(x, y)

    steps = g.get_steps()
    lastval = steps[1][4]
    assert abs(lastval - 1e-5 / 5.0) < 1e-10
    assert steps == [(1, 8 + 1, 1.0, 1.0, 0), (11, 18 + 1, 2.0, 2.0, lastval)]

    multi_g = Graph('foo', {})
    for x, y in vals:
        multi_g.add_data_point(x, [y, y, y])

    for s in multi_g.get_steps():
        assert s == steps


def test_graph_filename_sanitization():
    g = Graph('hello:world', {'a/a': 'b>b', 'c*c': 'd\0\0d'})
    assert g.path == os.path.join('graphs', 'a_a-b_b', 'c_c-d__d', 'hello_world')


def _sgn(x):
    return 1 if x >= 0 else -1
