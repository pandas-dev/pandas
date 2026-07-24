"""Test module for detecting uncollectable garbage in PyTables.

This test module *must* be loaded in the last place.  It just checks for
the existence of uncollectable garbage in ``gc.garbage`` after running
all the tests.

"""

import gc

from tables.tests import common


class GarbageTestCase(common.PyTablesTestCase):
    """Test for uncollectable garbage."""

    def test00(self):
        """Checking for uncollectable garbage."""

        garbageLen = len(gc.garbage)
        if garbageLen == 0:
            return  # success

        if common.verbose:
            classCount = {}
            # Count uncollected objects for each class.
            for obj in gc.garbage:
                objClass = obj.__class__.__name__
                if objClass in classCount:
                    classCount[objClass] += 1
                else:
                    classCount[objClass] = 1
            incidence = [
                "``%s``: %d" % (cls, cnt) for (cls, cnt) in classCount.items()
            ]
            print("Class incidence:", ", ".join(incidence))
        self.fail("Possible leak: %d uncollected objects." % garbageLen)


def suite():
    """Return a test suite consisting of all the test cases in the module."""

    theSuite = common.unittest.TestSuite()
    theSuite.addTest(common.make_suite(GarbageTestCase))
    return theSuite


if __name__ == "__main__":
    import sys

    common.parse_argv(sys.argv)
    common.print_versions()
    common.unittest.main(defaultTest="suite")
