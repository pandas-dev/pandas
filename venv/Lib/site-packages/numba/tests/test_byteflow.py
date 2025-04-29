"""
Test byteflow.py specific issues
"""
import unittest

from numba.tests.support import TestCase
from numba.core.compiler import run_frontend


class TestByteFlowIssues(TestCase):
    def test_issue_5087(self):
        # This is an odd issue. The exact number of print below is
        # necessary to trigger it. Too many or too few will alter the behavior.
        # Also note that the function below will not be executed. The problem
        # occurs at compilation. The definition below is invalid for execution.
        # The problem occurs in the bytecode analysis.
        def udt():
            print
            print
            print

            for i in range:
                print
                print
                print
                print
                print
                print
                print
                print
                print
                print
                print
                print
                print
                print
                print
                print
                print
                print

                for j in range:
                    print
                    print
                    print
                    print
                    print
                    print
                    print
                    for k in range:
                        for l in range:
                            print

                    print
                    print
                    print
                    print
                    print
                    print
                    print
                    print
                    print
                    if print:
                        for n in range:
                            print
                    else:
                        print

        run_frontend(udt)

    def test_issue_5097(self):
        # Inspired by https://github.com/numba/numba/issues/5097
        def udt():
            for i in range(0):
                if i > 0:
                    pass
                a = None # noqa: F841

        run_frontend(udt)

    def test_issue_5680(self):
        # From https://github.com/numba/numba/issues/5680#issuecomment-625351336
        def udt():
            for k in range(0):
                if 1 == 1:
                    ...
                if 'a' == 'a':
                    ...

        run_frontend(udt)


if __name__ == '__main__':
    unittest.main()
