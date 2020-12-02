import warnings

import pytest

import re

import pandas as pd
import pandas._testing as tm

from .base import BaseExtensionTests


def generate_error_message(dtype, op_name):
    """
    Generate an error message returned by the following

        getattr(s, op_name)(skipna=skipna)

    where `s` is an object with datatype `dtype`, `op_name` is an operation
    to be performed on `s`, and `skipna` is a boolean flag.

    Parameters
    ----------
    dtype : str
        Datatype of the object `s`
    op_name : str
        Name of operation to perform on `s`

    Returns
    -------
    str
        Error message caused by performing the operation
    """
    if dtype == "category":
        if op_name in ["min", "max"]:
            return re.escape(
                f"Categorical is not ordered for operation {op_name}\n"
                "you can use .as_ordered() to change the Categorical to an "
                "ordered one\n"
            )
        else:
            return f"'Categorical' does not implement reduction '{op_name}'"
    elif dtype in ["string", "arrow_string"]:
        return f"Cannot perform reduction '{op_name}' with string dtype"
    elif dtype == "interval":
        return rf"cannot perform {op_name} with type interval\[float64\]"
    elif dtype == "json":
        return f"cannot perform {op_name} with type json"
    elif dtype == "arrow_bool":
        return ""


class BaseReduceTests(BaseExtensionTests):
    """
    Reduction specific tests. Generally these only
    make sense for numeric/boolean operations.
    """

    def check_reduce(self, s, op_name, skipna):
        result = getattr(s, op_name)(skipna=skipna)
        expected = getattr(s.astype("float64"), op_name)(skipna=skipna)
        tm.assert_almost_equal(result, expected)


class BaseNoReduceTests(BaseReduceTests):
    """ we don't define any reductions """

    @pytest.mark.parametrize("skipna", [True, False])
    def test_reduce_series_numeric(self, data, all_numeric_reductions, skipna):
        op_name = all_numeric_reductions
        s = pd.Series(data)

        msg = generate_error_message(s.dtype.name, op_name)

        with pytest.raises(TypeError, match=msg):
            getattr(s, op_name)(skipna=skipna)

    @pytest.mark.parametrize("skipna", [True, False])
    def test_reduce_series_boolean(self, data, all_boolean_reductions, skipna):
        op_name = all_boolean_reductions
        s = pd.Series(data)

        msg = generate_error_message(s.dtype.name, op_name)

        with pytest.raises(TypeError, match=msg):
            getattr(s, op_name)(skipna=skipna)


class BaseNumericReduceTests(BaseReduceTests):
    @pytest.mark.parametrize("skipna", [True, False])
    def test_reduce_series(self, data, all_numeric_reductions, skipna):
        op_name = all_numeric_reductions
        s = pd.Series(data)

        # min/max with empty produce numpy warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            self.check_reduce(s, op_name, skipna)


class BaseBooleanReduceTests(BaseReduceTests):
    @pytest.mark.parametrize("skipna", [True, False])
    def test_reduce_series(self, data, all_boolean_reductions, skipna):
        op_name = all_boolean_reductions
        s = pd.Series(data)
        self.check_reduce(s, op_name, skipna)
