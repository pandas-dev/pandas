import sys

import unittest
from unittest import TestCase

import faulthandler


try:
    # May fail in IPython Notebook with UnsupportedOperation
    faulthandler.enable()
except BaseException as e:
    msg = "Failed to enable faulthandler due to:\n{err}"
    warnings.warn(msg.format(err=e))


# Try to inject Numba's unittest customizations.
from llvmlite.tests import customize


def discover_tests(startdir):
    """Discover test under a directory
    """
    # Avoid importing unittest
    loader = unittest.TestLoader()
    suite = loader.discover(startdir)
    return suite


def run_tests(suite=None, xmloutput=None, verbosity=1):
    """
    args
    ----
    - suite [TestSuite]
        A suite of all tests to run
    - xmloutput [str or None]
        Path of XML output directory (optional)
    - verbosity [int]
        Verbosity level of tests output

    Returns the TestResult object after running the test *suite*.
    """
    if suite is None:
        suite = discover_tests("llvmlite.tests")
    if xmloutput is not None:
        import xmlrunner
        runner = xmlrunner.XMLTestRunner(output=xmloutput)
    else:
        runner = None
    prog = unittest.main(suite=suite, testRunner=runner, exit=False,
                         verbosity=verbosity)
    return prog.result


def main():
    res = run_tests()
    sys.exit(0 if res.wasSuccessful() else 1)
