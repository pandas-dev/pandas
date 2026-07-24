from os.path import dirname
import unittest
from unittest.suite import TestSuite

from numba.testing import load_testsuite

def load_tests(loader, tests, pattern):
    suite = TestSuite()
    suite.addTests(load_testsuite(loader, dirname(__file__)))
    return suite
