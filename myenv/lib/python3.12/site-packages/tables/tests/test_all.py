"""Run all test cases."""

import sys

import numpy as np
from packaging.version import Version

import tables as tb
from tables.tests import common
from tables.tests.test_suite import suite, test


def get_tuple_version(hexversion):
    """Get a tuple from a compact version in hex."""

    h = hexversion
    return(h & 0xff0000) >> 16, (h & 0xff00) >> 8, h & 0xff


if __name__ == '__main__':

    common.parse_argv(sys.argv)

    hdf5_version = get_tuple_version(tb.which_lib_version("hdf5")[0])
    hdf5_version_str = "%s.%s.%s" % hdf5_version
    if Version(hdf5_version_str) < tb.req_versions.min_hdf5_version:
        print(f"*Warning*: HDF5 version is lower than recommended: "
              f"{hdf5_version} < {tb.req_versions.min_hdf5_version}")

    if Version(np.__version__) < tb.req_versions.min_numpy_version:
        print(f"*Warning*: NumPy version is lower than recommended: "
              f"{np.__version__} < {tb.req_versions.min_numpy_version}")

    # Handle some global flags (i.e. only useful for test_all.py)
    only_versions = 0
    args = sys.argv[:]
    for arg in args:
        # Remove 'show-versions' for PyTables 2.3 or higher
        if arg in ['--print-versions', '--show-versions']:
            only_versions = True
            sys.argv.remove(arg)
        elif arg == '--show-memory':
            common.show_memory = True
            sys.argv.remove(arg)

    common.print_versions()
    if not only_versions:
        common.print_heavy(common.heavy)
        common.unittest.main(defaultTest='tb.tests.suite')
