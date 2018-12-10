from textwrap import dedent

import pytest

from pandas.util._decorators import deprecate
import pandas.util.testing as tm


def new_func():
    """
    This is the summary. The deprecate directive goes next.

    This is the extended summary. The deprecate directive goes before this.
    """
    return 1234


def new_func_no_docstring():
    pass


def new_func_with_deprecation():
    """
    This is the summary. The deprecate directive goes next.

    .. deprecated:: 1.0
        Use new_func instead.

    This is the extended summary. The deprecate directive goes before this.
    """
    pass


def test_deprecate_ok():
    depr_func = deprecate('depr_func', new_func, '1.0',
                          msg='Use new_func instead.')

    with tm.assert_produces_warning(FutureWarning):
        result = depr_func()

    assert result == 1234

    assert depr_func.__doc__ == dedent(new_func_with_deprecation.__doc__)


def test_deprecate_no_docstring():
    with pytest.raises(ValueError, match='deprecate needs a correctly '
                                         'formatted docstring'):
        deprecate('depr_func', new_func_no_docstring, '1.0',
                  msg='Use new_func instead.')
