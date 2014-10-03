from __future__ import division, absolute_import, print_function

import sys

if sys.version_info[0] < 3:
    from .__version__ import version as __version__
    # Must import local ccompiler ASAP in order to get
    # customized CCompiler.spawn effective.
    from . import ccompiler
    from . import unixccompiler

    from .info import __doc__
    from .npy_pkg_config import *

    try:
        import __config__
        _INSTALLED = True
    except ImportError:
        _INSTALLED = False
else:
    from numpy.distutils.__version__ import version as __version__
    # Must import local ccompiler ASAP in order to get
    # customized CCompiler.spawn effective.
    import numpy.distutils.ccompiler
    import numpy.distutils.unixccompiler

    from numpy.distutils.info import __doc__
    from numpy.distutils.npy_pkg_config import *

    try:
        import numpy.distutils.__config__
        _INSTALLED = True
    except ImportError:
        _INSTALLED = False

if _INSTALLED:
    from numpy.testing import Tester
    test = Tester().test
    bench = Tester().bench
