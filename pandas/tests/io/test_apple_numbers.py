""" Test Apple Numbers read and write """
import pytest

import pandas._testing as tm

from pandas.io.apple_numbers import read_apple_numbers, to_apple_numbers  # isort:skip


@pytest.mark.single_cpu
class TestAppleNumbers:
    def check_error_on_write(self, df, exc, err_msg):
        # check that we are raising the exception on writing

        with pytest.raises(exc, match=err_msg):
            with tm.ensure_clean() as path:
                to_apple_numbers(df, path)

    def check_error_on_read(self, exc, err_msg):
        # check that we are raising the exception on reading

        with pytest.raises(exc, match=err_msg):
            with tm.ensure_clean() as path:
                _ = read_apple_numbers(path)

    def test_unimplemented(self):
        msg = "Apple Numbers support is in development"
        self.check_error_on_read(NotImplementedError, msg)
        self.check_error_on_write(object(), NotImplementedError, msg)
