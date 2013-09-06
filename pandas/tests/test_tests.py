#!/usr/bin/python
# -*- coding: utf-8 -*-
import pandas as pd
import unittest
import warnings
import nose
import sys

from pandas.util.testing import (
    assert_almost_equal, assertRaisesRegexp, raise_with_traceback
)

# let's get meta.

class TestUtilTesting(unittest.TestCase):
    _multiprocess_can_split_ = True

    def test_assert_almost_equal(self):
        # don't die because values are not ndarrays
        assert_almost_equal(1.1,1.1,check_less_precise=True)

    def test_raise_with_traceback(self):
        with assertRaisesRegexp(LookupError, "error_text"):
            try:
                raise ValueError("THIS IS AN ERROR")
            except ValueError as e:
                e = LookupError("error_text")
                raise_with_traceback(e)
        with assertRaisesRegexp(LookupError, "error_text"):
            try:
                raise ValueError("This is another error")
            except ValueError:
                e = LookupError("error_text")
                _, _, traceback = sys.exc_info()
                raise_with_traceback(e, traceback)
