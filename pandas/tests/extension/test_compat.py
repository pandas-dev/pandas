"""
Tests for downstream compatibility.
"""
import pandas as pd
import pandas._testing as tm


class SingleFormatter(pd.arrays.StringArray):
    """A StringArray with the old single-argument of _formatter"""

    def _formatter(self, boxed: bool = False):
        breakpoint()
        return super()._formatter(boxed)


def test_single_argument_formatter():
    s = SingleFormatter._from_sequence(["a", "b", "c"])
    with tm.assert_produces_warning(DeprecationWarning):
        str(s)
        repr(s)
        str(pd.Series(s))
        repr(pd.Series(s))
