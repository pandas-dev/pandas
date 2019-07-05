import warnings

import pytest

import pandas.util.testing as tm


def f():
    warnings.warn("f1", FutureWarning)
    warnings.warn("f2", RuntimeWarning)


@pytest.mark.filterwarnings("ignore:f1:FutureWarning")
def test_assert_produces_warning_honors_filter():
    # Raise by default.
    msg = r"Caused unexpected warning\(s\)"
    with pytest.raises(AssertionError, match=msg):
        with tm.assert_produces_warning(RuntimeWarning):
            f()

    with tm.assert_produces_warning(RuntimeWarning, raise_on_extra_warnings=False):
        f()
