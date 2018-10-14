"""Tests for Interval-Interval operations, such as overlaps, contains, etc."""
import pytest

from pandas import IntervalIndex
from ...arrays.interval.interval_ops import BaseOverlaps


class TestOverlaps(BaseOverlaps):

    @pytest.fixture
    def constructor(self):
        """
        Fixture for IntervalIndex class constructor (used by parent class)
        """
        return IntervalIndex
