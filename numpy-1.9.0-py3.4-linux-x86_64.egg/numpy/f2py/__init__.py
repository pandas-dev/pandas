#!/usr/bin/env python
from __future__ import division, absolute_import, print_function

__all__ = ['run_main', 'compile', 'f2py_testing']

import os
import sys
import subprocess

from . import f2py2e
from . import f2py_testing
from . import diagnose

from .info import __doc__

run_main = f2py2e.run_main
main = f2py2e.main

def compile(source,
            modulename = 'untitled',
            extra_args = '',
            verbose = 1,
            source_fn = None
            ):
    ''' Build extension module from processing source with f2py.
    Read the source of this function for more information.
    '''
    from numpy.distutils.exec_command import exec_command
    import tempfile
    if source_fn is None:
        f = tempfile.NamedTemporaryFile(suffix='.f')
    else:
        f = open(source_fn, 'w')

    try:
        f.write(source)
        f.flush()

        args = ' -c -m %s %s %s'%(modulename, f.name, extra_args)
        c = '%s -c "import numpy.f2py as f2py2e;f2py2e.main()" %s' % \
                (sys.executable, args)
        s, o = exec_command(c)
    finally:
        f.close()
    return s

from numpy.testing import Tester
test = Tester().test
bench = Tester().bench
