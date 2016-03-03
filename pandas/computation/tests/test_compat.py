#!/usr/bin/env python

# flake8: noqa

import nose
from itertools import product
from distutils.version import LooseVersion

import pandas as pd
from pandas.util import testing as tm

from pandas.computation.engines import _engines
import pandas.computation.expr as expr

ENGINES_PARSERS = list(product(_engines, expr._parsers))


def test_compat():
    # test we have compat with our version of nu

    from pandas.computation import _NUMEXPR_INSTALLED
    try:
        import numexpr as ne
        ver = ne.__version__
        if ver == LooseVersion('2.4.4'):
            assert not _NUMEXPR_INSTALLED
        elif ver < LooseVersion('2.1'):
            with tm.assert_produces_warning(UserWarning,
                                            check_stacklevel=False):
                assert not _NUMEXPR_INSTALLED
        else:
            assert _NUMEXPR_INSTALLED

    except ImportError:
        raise nose.SkipTest("not testing numexpr version compat")


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
            raise nose.SkipTest("no numexpr")
        else:
            if ne.__version__ < LooseVersion('2.1'):
                with tm.assertRaisesRegexp(ImportError, "'numexpr' version is "
                                           ".+, must be >= 2.1"):
                    testit()
            elif ne.__version__ == LooseVersion('2.4.4'):
                raise nose.SkipTest("numexpr version==2.4.4")
            else:
                testit()
    else:
        testit()


if __name__ == '__main__':
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
