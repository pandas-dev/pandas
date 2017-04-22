from __future__ import division

import pytest
from pandas import Interval
import pandas.util.testing as tm


class TestInterval(tm.TestCase):
    def setUp(self):
        self.interval = Interval(0, 1)

    def test_properties(self):
        self.assertEqual(self.interval.closed, 'right')
        self.assertEqual(self.interval.left, 0)
        self.assertEqual(self.interval.right, 1)
        self.assertEqual(self.interval.mid, 0.5)

    def test_repr(self):
        self.assertEqual(repr(self.interval),
                         "Interval(0, 1, closed='right')")
        self.assertEqual(str(self.interval), "(0, 1]")

        interval_left = Interval(0, 1, closed='left')
        self.assertEqual(repr(interval_left),
                         "Interval(0, 1, closed='left')")
        self.assertEqual(str(interval_left), "[0, 1)")

    def test_contains(self):
        self.assertIn(0.5, self.interval)
        self.assertIn(1, self.interval)
        self.assertNotIn(0, self.interval)
        pytest.raises(TypeError, lambda: self.interval in self.interval)

        interval = Interval(0, 1, closed='both')
        self.assertIn(0, interval)
        self.assertIn(1, interval)

        interval = Interval(0, 1, closed='neither')
        self.assertNotIn(0, interval)
        self.assertIn(0.5, interval)
        self.assertNotIn(1, interval)

    def test_equal(self):
        self.assertEqual(Interval(0, 1), Interval(0, 1, closed='right'))
        self.assertNotEqual(Interval(0, 1), Interval(0, 1, closed='left'))
        self.assertNotEqual(Interval(0, 1), 0)

    def test_comparison(self):
        with self.assertRaisesRegexp(TypeError, 'unorderable types'):
            Interval(0, 1) < 2

        self.assertTrue(Interval(0, 1) < Interval(1, 2))
        self.assertTrue(Interval(0, 1) < Interval(0, 2))
        self.assertTrue(Interval(0, 1) < Interval(0.5, 1.5))
        self.assertTrue(Interval(0, 1) <= Interval(0, 1))
        self.assertTrue(Interval(0, 1) > Interval(-1, 2))
        self.assertTrue(Interval(0, 1) >= Interval(0, 1))

    def test_hash(self):
        # should not raise
        hash(self.interval)

    def test_math_add(self):
        expected = Interval(1, 2)
        actual = self.interval + 1
        self.assertEqual(expected, actual)

        expected = Interval(1, 2)
        actual = 1 + self.interval
        self.assertEqual(expected, actual)

        actual = self.interval
        actual += 1
        self.assertEqual(expected, actual)

        with pytest.raises(TypeError):
            self.interval + Interval(1, 2)

        with pytest.raises(TypeError):
            self.interval + 'foo'

    def test_math_sub(self):
        expected = Interval(-1, 0)
        actual = self.interval - 1
        self.assertEqual(expected, actual)

        actual = self.interval
        actual -= 1
        self.assertEqual(expected, actual)

        with pytest.raises(TypeError):
            self.interval - Interval(1, 2)

        with pytest.raises(TypeError):
            self.interval - 'foo'

    def test_math_mult(self):
        expected = Interval(0, 2)
        actual = self.interval * 2
        self.assertEqual(expected, actual)

        expected = Interval(0, 2)
        actual = 2 * self.interval
        self.assertEqual(expected, actual)

        actual = self.interval
        actual *= 2
        self.assertEqual(expected, actual)

        with pytest.raises(TypeError):
            self.interval * Interval(1, 2)

        with pytest.raises(TypeError):
            self.interval * 'foo'

    def test_math_div(self):
        expected = Interval(0, 0.5)
        actual = self.interval / 2.0
        self.assertEqual(expected, actual)

        actual = self.interval
        actual /= 2.0
        self.assertEqual(expected, actual)

        with pytest.raises(TypeError):
            self.interval / Interval(1, 2)

        with pytest.raises(TypeError):
            self.interval / 'foo'
