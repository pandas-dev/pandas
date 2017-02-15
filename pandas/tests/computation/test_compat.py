
# flake8: noqa

import pytest
from itertools import product
from distutils.version import LooseVersion

import pandas as pd
from pandas.util import testing as tm

from pandas.computation.engines import _engines
import pandas.computation.expr as expr
from pandas.computation import _MIN_NUMEXPR_VERSION
# Get reload for Python 3.4 and on, if not, use internal reload()
try:
    from importlib import reload
except ImportError:
    pass

ENGINES_PARSERS = list(product(_engines, expr._parsers))


def test_compat():
    # test we have compat with our version of nu

    from pandas.computation import _NUMEXPR_INSTALLED
    try:
        import numexpr as ne
        ver = ne.__version__
        if ver < LooseVersion(_MIN_NUMEXPR_VERSION):
            with tm.assert_produces_warning(UserWarning,
                                            check_stacklevel=False):
                reload(pd.computation)
            assert not _NUMEXPR_INSTALLED
        else:
            assert _NUMEXPR_INSTALLED

    except ImportError:
        pytest.skip("not testing numexpr version compat")


def test_invalid_numexpr_version():
    for engine, parser in ENGINES_PARSERS:
        yield check_invalid_numexpr_version, engine, parser


def check_invalid_numexpr_version(engine, parser):
    def testit():
        a, b = 1, 2
        res = pd.eval('a + b', engine=engine, parser=parser)
        tm.assert_equal(res, 3)

    if engine == 'numexpr':
        try:
            import numexpr as ne
        except ImportError:
            pytest.skip("no numexpr")
        else:
            if ne.__version__ < LooseVersion(_MIN_NUMEXPR_VERSION):
                with tm.assertRaises(ImportError):
                    testit()
            else:
                testit()
    else:
        testit()
