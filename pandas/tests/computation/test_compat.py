from distutils.version import LooseVersion

import pytest

from pandas.compat._optional import VERSIONS

import pandas as pd
from pandas.core.computation.engines import ENGINES
import pandas.core.computation.expr as expr


def test_compat():
    # test we have compat with our version of nu

    from pandas.core.computation.check import NUMEXPR_INSTALLED

    try:
        import numexpr as ne

        ver = ne.__version__
        if LooseVersion(ver) < LooseVersion(VERSIONS["numexpr"]):
            assert not NUMEXPR_INSTALLED
        else:
            assert NUMEXPR_INSTALLED
    except ImportError:
        pytest.skip("not testing numexpr version compat")


@pytest.mark.parametrize("engine", ENGINES)
@pytest.mark.parametrize("parser", expr.PARSERS)
def test_invalid_numexpr_version(engine, parser):
    def testit():
        a, b = 1, 2  # noqa
        res = pd.eval("a + b", engine=engine, parser=parser)
        assert res == 3

    if engine == "numexpr":
        try:
            import numexpr as ne
        except ImportError:
            pytest.skip("no numexpr")
        else:
            if LooseVersion(ne.__version__) < LooseVersion(VERSIONS["numexpr"]):
                # TODO comment this back in once we know the exception message
                # with pytest.raises(ImportError):
                testit()
            else:
                testit()
    else:
        testit()
