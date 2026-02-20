"""Test suite consisting of all testcases."""

import sys

from tables.tests import common


def suite():
    test_modules = [
        "tables.tests.test_attributes",
        "tables.tests.test_basics",
        "tables.tests.test_create",
        "tables.tests.test_backcompat",
        "tables.tests.test_types",
        "tables.tests.test_lists",
        "tables.tests.test_tables",
        "tables.tests.test_tablesMD",
        "tables.tests.test_large_tables",
        "tables.tests.test_array",
        "tables.tests.test_earray",
        "tables.tests.test_carray",
        "tables.tests.test_vlarray",
        "tables.tests.test_tree",
        "tables.tests.test_timetype",
        "tables.tests.test_do_undo",
        "tables.tests.test_enum",
        "tables.tests.test_nestedtypes",
        "tables.tests.test_hdf5compat",
        "tables.tests.test_numpy",
        "tables.tests.test_queries",
        "tables.tests.test_expression",
        "tables.tests.test_links",
        "tables.tests.test_indexes",
        "tables.tests.test_indexvalues",
        "tables.tests.test_index_backcompat",
        "tables.tests.test_aux",
        "tables.tests.test_utils",
        "tables.tests.test_direct_chunk",
        # Sub-packages
        "tables.nodes.tests.test_filenode",
    ]

    # print('-=' * 38)

    # The test for garbage must be run *in the last place*.
    # Else, it is not as useful.
    test_modules.append("tables.tests.test_garbage")

    alltests = common.unittest.TestSuite()
    if common.show_memory:
        # Add a memory report at the beginning
        alltests.addTest(common.make_suite(common.ShowMemTime))
    for name in test_modules:
        # Unexpectedly, the following code doesn't seem to work anymore
        # in python 3
        # exec('from %s import suite as test_suite' % name)
        __import__(name)
        test_suite = sys.modules[name].suite

        alltests.addTest(test_suite())
        if common.show_memory:
            # Add a memory report after each test module
            alltests.addTest(common.make_suite(common.ShowMemTime))
    return alltests


def test(verbose=False, heavy=False):
    """Run all the tests in the test suite.

    If *verbose* is set, the test suite will emit messages with full
    verbosity (not recommended unless you are looking into a certain
    problem).

    If *heavy* is set, the test suite will be run in *heavy* mode (you
    should be careful with this because it can take a lot of time and
    resources from your computer).

    Return 0 (os.EX_OK) if all tests pass, 1 in case of failure

    """

    common.print_versions()
    common.print_heavy(heavy)

    # What a context this is!
    # oldverbose, common.verbose = common.verbose, verbose
    oldheavy, common.heavy = common.heavy, heavy
    try:
        result = common.unittest.TextTestRunner(
            verbosity=1 + int(verbose)
        ).run(suite())
        if result.wasSuccessful():
            return 0
        else:
            return 1
    finally:
        # common.verbose = oldverbose
        common.heavy = oldheavy  # there are pretty young heavies, too ;)
