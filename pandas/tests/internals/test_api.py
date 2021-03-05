"""
Tests for the pseudo-public API implemented in internals/api.py and exposed
in core.internals
"""

from pandas.core import internals
from pandas.core.internals import api


def test_internals_api():
    assert internals.make_block is api.make_block
