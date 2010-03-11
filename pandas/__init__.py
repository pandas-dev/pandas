# pylint: disable-msg=W0614,W0401,W0611,W0622

__docformat__ = 'restructuredtext'

from datetime import datetime

import numpy as np

from pandas.version import __version__
from pandas.info import __doc__

from pandas.core.api import *
from pandas.io.parsers import parseCSV, parseText, parseExcel
from pandas.stats.api import *

from numpy.testing import Tester
class NoseWrapper(Tester):
    '''
    This is simply a monkey patch for numpy.testing.Tester, so that extra_argv can
    be changed from its default None to ['--exe'] so that the tests can be run
    the same across platforms.
    '''
    def test(self, label='fast', verbose=1, extra_argv=['--exe'], doctests=False,
             coverage=False):
        ''' Run tests for module using nose

        %(test_header)s
        doctests : boolean
            If True, run doctests in module, default False
        coverage : boolean
            If True, report coverage of NumPy code, default False
            (Requires the coverage module:
             http://nedbatchelder.com/code/modules/coverage.html)
        '''

        # cap verbosity at 3 because nose becomes *very* verbose beyond that
        verbose = min(verbose, 3)

        from numpy.testing import utils
        utils.verbose = verbose

        if doctests:
            print "Running unit tests and doctests for %s" % self.package_name
        else:
            print "Running unit tests for %s" % self.package_name

        self._show_system_info()

        # reset doctest state on every run
        import doctest
        doctest.master = None

        argv, plugins = self.prepare_test_args(label, verbose, extra_argv,
                                               doctests, coverage)
        from numpy.testing.noseclasses import NumpyTestProgram
        t = NumpyTestProgram(argv=argv, exit=False, plugins=plugins)
        return t.result
test = NoseWrapper().test
