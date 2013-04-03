#!/usr/bin/python
# -*- coding: utf-8 -*-
from __future__ import with_statement  # support python 2.5
import pandas as pd
import unittest
import warnings
import nose

from pandas.util.testing import assert_almost_equal

# let's get meta.

class TestUtilTesting(unittest.TestCase):
    _multiprocess_can_split_ = True

    def __init__(self, *args):
        super(TestUtilTesting, self).__init__(*args)

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_assert_almost_equal(self):
        # don't die because values are not ndarrays
        assert_almost_equal(1.1,1.1,check_less_precise=True)
