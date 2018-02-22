import pytest
import numpy as np
from pandas import Interval, IntervalIndex
import pandas.util.testing as tm


# closed abbreviations mapping, used to compactify parameter setup
closed_abbr = {'l': 'left', 'r': 'right', 'b': 'both', 'n': 'neither'}


def make_scalar_fixture(params):
    """make a parametrized fixture that returns (index, scalar) pairs"""
    def scalar_test_id(request):
        """pytest verbose output display: closed-scalar, e.g. left-0.5"""
        closed = closed_abbr[request[0]]
        return '{c}-{s}'.format(s=request[1], c=closed)

    @pytest.fixture(ids=scalar_test_id, params=params)
    def scalar_fixture(self, request):
        """
        expects params of the form (index closed abbreviation, scalar)
        """
        idx_closed = closed_abbr[request.param[0]]
        return self.index(idx_closed), request.param[1]
    return scalar_fixture


def make_interval_fixture(params):
    """make a parametrized fixture that returns (index, interval) pairs"""
    def interval_test_id(request):
        """pytest verbose output display: closed-interval, e.g. left-[0, 1)"""
        idx_closed = closed_abbr[request[0]]
        abbr, left, right = request[1:]
        interval = Interval(left, right, closed_abbr[abbr])
        return '{idx}-{iv}'.format(idx=idx_closed, iv=interval)

    @pytest.fixture(ids=interval_test_id, params=params)
    def interval_fixture(self, request):
        """
        expects params with abbreviated closed values of the form
        (index closed, interval closed, interval left, interval right)
        """
        index = self.index(closed_abbr[request[0]])
        abbr, left, right = request[1:]
        interval = Interval(left, right, closed_abbr[abbr])
        return index, interval
    return interval_fixture


def generate_invalid_intervals(cls):
    """
    create a list of invalid intervals, with values of the form
    (index closed, interval closed, interval left, interval right)
    """
    # intervals are invalid if the bounds are the same, but different closed
    interval_tuples_invalid = [
        (c1, c2, l, r) for c1 in 'lrbn' for c2 in 'lrbn'
        for l, r in set(cls.tuples) if c1 != c2]

    # stricly contained intervals are invalid regardless of closed
    invalid_bounds = [(cls.tuples[0][0] + 0.1, cls.tuples[0][1] - 0.1)]

    # strictly covering intervals are invalid regardless of closed
    lhs = min(l for l, r in cls.tuples) - 1
    rhs = max(r for l, r in cls.tuples) + 1
    invalid_bounds.append((lhs, rhs))

    # disjoint intervals are invalid regardless of closed
    invalid_bounds.append((rhs, rhs + 1))

    # add all closed combinations to the bounds and extend
    interval_tuples_invalid += [
        (c1, c2, l, r) for c1 in 'lrbn' for c2 in 'lrbn'
        for l, r in invalid_bounds]
    return interval_tuples_invalid


def add_class_fixtures(cls):
    """
    add scalar/interval test fixutres to the test classes; scalar fixtures
    created from class defined lists, intervals generated cls.tuples
    """
    # create scalar fixtures based on the class defined valid/invalid tuples
    cls.scalars_valid = make_scalar_fixture(cls.scalar_tuples_valid)
    cls.scalars_invalid = make_scalar_fixture(cls.scalar_tuples_invalid)

    # create valid intervals fixture (only exact matches are valid)
    interval_tuples_valid = [
        (c, c, l, r) for c in 'lrbn' for l, r in set(cls.tuples)]
    cls.intervals_valid = make_interval_fixture(interval_tuples_valid)

    # create invalid intervals fixture (non-exact matches are invalid)
    interval_tuples_invalid = generate_invalid_intervals(cls)
    cls.intervals_invalid = make_interval_fixture(interval_tuples_invalid)
    return cls


class Base(object):
    """
    Common tests for IntervalIndex operations that differ depending on the
    structure of the underlying intervals, e.g. monotonic vs. non-monotonic
    """

    @pytest.fixture(params=['left', 'right', 'both', 'neither'])
    def closed(self, request):
        """Valid values of 'closed' for an InteralIndex"""
        return request.param

    @pytest.fixture
    def index(self, closed):
        """
        Create an IntervalIndex from the class specified tuples
        """
        return IntervalIndex.from_tuples(self.tuples, closed=closed)

    def test_get_loc_scalar(self, scalars_valid):
        index, scalar = scalars_valid
        result = index.get_loc(scalar)
        expected = [i for i, iv in enumerate(index) if scalar in iv]
        if len(expected) == 1:
            assert result == expected[0]
        else:
            expected = np.array(expected, dtype='int64')
            tm.assert_numpy_array_equal(result, expected)

    def test_get_loc_scalar_errors(self, scalars_invalid):
        index, scalar = scalars_invalid
        with pytest.raises(KeyError):
            index.get_loc(scalar)

    @pytest.mark.skip()
    def test_get_loc_interval(self, intervals_valid):
        index, interval = intervals_valid
        result = index.get_loc(interval)
        expected = [i for i, iv in enumerate(index) if interval == iv]
        if len(expected) == 1:
            assert result == expected[0]
        else:
            expected = np.array(expected, dtype='int64')
            tm.assert_numpy_array_equal(result, expected)

    @pytest.mark.skip()
    def test_get_loc_interval_errors(self, intervals_invalid):
        index, interval = intervals_invalid
        with pytest.raises(KeyError):
            index.get_loc(interval)


class TestNonOverlappingMonotonic(Base):
    """Tests specifics to non-overlapping monotonic IntervalIndex"""
    # tuples from which to create the underlying IntervalIndex
    tuples = [(0, 1), (2, 3)]

    # valid scalar tuples to test, of the from (closed abbr, scalar)
    # (a scalar is valid if it's contained in an interval in the index)
    scalar_tuples_valid = [
        ('l', 0), ('l', 0.5), ('l', 2), ('l', 2.5),
        ('r', 0.5), ('r', 1), ('r', 2.5), ('r', 3),
        ('b', 0), ('b', 0.5), ('b', 1), ('b', 2), ('b', 2.5), ('b', 3),
        ('n', 0.5), ('n', 2.5)]

    # invalid scalar tuples to test, of the from (closed abbr, scalar)
    scalar_tuples_invalid = [
        ('l', -0.5), ('l', 1), ('l', 1.5), ('l', 3), ('l', 3.5),
        ('r', -0.5), ('r', 0), ('r', 1.5), ('r', 2), ('r', 3.5),
        ('b', -0.5), ('b', 1.5), ('b', 3.5),
        ('n', -0.5), ('n', 0), ('n', 1), ('n', 1.5), ('n', 2), ('n', 3),
        ('n', 3.5)]


TestNonOverlappingMonotonic = add_class_fixtures(TestNonOverlappingMonotonic)


class TestOverlapping(Base):
    """Tests specifics to overlapping IntervalIndex"""

    # tuples from which to create the underlying IntervalIndex
    tuples = [(0, 2), (1, 3)]

    # valid scalar tuples to test, of the from (closed abbr, scalar)
    # (a scalar is valid if it's contained in an interval in the index)
    scalar_tuples_valid = [
        ('l', 0), ('l', 0.5), ('l', 1), ('l', 2), ('l', 2.5),
        ('r', 0.5), ('r', 1), ('r', 2), ('r', 2.5), ('r', 3),
        ('b', 0), ('b', 0.5), ('b', 1), ('b', 2), ('b', 2.5), ('b', 3),
        ('n', 0.5), ('n', 1), ('n', 1), ('n', 2.5)]

    # invalid scalar tuples to test, of the from (closed abbr, scalar)
    scalar_tuples_invalid = [
        ('l', -0.5), ('l', 3), ('l', 3.5),
        ('r', -0.5), ('r', 0), ('r', 3.5),
        ('b', -0.5), ('b', 3.5),
        ('n', -0.5), ('n', 0), ('n', 3.5)]


TestOverlapping = add_class_fixtures(TestOverlapping)


class TestNonMonotonicDupes(Base):
    """Tests specifics to non-monotonic IntervalIndex with dupes"""

    # tuples from which to create the underlying IntervalIndex
    tuples = [(0, 1), (2, 3), (0, 1)]

    # valid scalar tuples to test, of the from (closed abbr, scalar)
    # (a scalar is valid if it's contained in an interval in the index)
    scalar_tuples_valid = [
        ('l', 0), ('l', 0.5), ('l', 2), ('l', 2.5),
        ('r', 0.5), ('r', 1), ('r', 2.5), ('r', 3),
        ('b', 0), ('b', 0.5), ('b', 1), ('b', 2), ('b', 2.5), ('b', 3),
        ('n', 0.5), ('n', 2.5)]

    # invalid scalar tuples to test, of the from (closed abbr, scalar)
    scalar_tuples_invalid = [
        ('l', -0.5), ('l', 1), ('l', 1.5), ('l', 3), ('l', 3.5),
        ('r', -0.5), ('r', 0), ('r', 1.5), ('r', 2), ('r', 3.5),
        ('b', -0.5), ('b', 1.5), ('b', 3.5),
        ('n', -0.5), ('n', 0), ('n', 1), ('n', 1.5), ('n', 2), ('n', 3),
        ('n', 3.5)]


TestNonMonotonicDupes = add_class_fixtures(TestNonMonotonicDupes)
