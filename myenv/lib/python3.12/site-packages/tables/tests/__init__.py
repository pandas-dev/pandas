"""Unit tests for PyTables.

This package contains some modules which provide a ``suite()`` function
(with no arguments) which returns a test suite for some PyTables
functionality.

"""

from tables.tests.common import print_versions
from tables.tests.test_suite import test, suite
# Necessary for the test suite
import tables.req_versions
