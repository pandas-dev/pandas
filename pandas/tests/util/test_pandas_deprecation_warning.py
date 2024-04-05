import warnings

from pandas.errors import PandasChangeWarning
from pandas.util._decorators import deprecate_kwarg

import pandas._testing as tm


def f1():
    warnings.warn("f1", PandasChangeWarning)


def test_function_warns_pandas_deprecation_warning():
    with tm.assert_produces_warning(PandasChangeWarning):
        f1()


@deprecate_kwarg("old", klass=PandasChangeWarning, new_arg_name="new")
def f2(new=0):
    return new


def test_decorator_warns_pandas_deprecation_warning():
    with tm.assert_produces_warning(PandasChangeWarning):
        f2(old=1)
