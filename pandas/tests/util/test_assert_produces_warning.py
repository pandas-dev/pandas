import warnings

import pytest

import pandas.util.testing as tm


def f(a=FutureWarning, b=RuntimeWarning):
    warnings.warn('f1', a)
    warnings.warn('f2', b)


@pytest.mark.filterwarnings('ignore:f1:FutureWarning')
@pytest.mark.filterwarnings('ignore:f2:RuntimeWarning')
def test_assert_produces_warning_honors_filter():
    with tm.assert_produces_warning(RuntimeWarning):
        f()


@pytest.mark.filterwarnings('ignore:f1:FutureWarning')
def test_assert_produces_warning_message():
    with tm.assert_produces_warning(FutureWarning, message='f2'):
        f(FutureWarning, FutureWarning)
