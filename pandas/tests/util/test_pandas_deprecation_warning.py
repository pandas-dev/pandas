import warnings

from pandas.util._decorators import deprecate_kwarg
from pandas.util._exceptions import Pandas40DeprecationWarning

import pandas._testing as tm


def f1():
    warnings.warn("f1", Pandas40DeprecationWarning)


def test_function_warns_pandas_deprecation_warning():
    with tm.assert_produces_warning(DeprecationWarning):
        f1()


@deprecate_kwarg("old", "new")
def f2(new=0):
    return new


def test_decorator_warns_pandas_deprecation_warning():
    with tm.assert_produces_warning(DeprecationWarning):
        f2(old=1)
