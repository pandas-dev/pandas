"""Tests for Interval-Interval operations, such as overlaps, contains, etc."""
import pytest

from pandas.core.arrays import IntervalArray
from .interval_ops import BaseOverlaps


class TestOverlaps(BaseOverlaps):

    @pytest.fixture
    def constructor(self):
        """
        Fixture for IntervalArray class constructor (used by parent class)
        """
        return IntervalArray
