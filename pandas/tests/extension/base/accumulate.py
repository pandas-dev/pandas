import warnings

import pytest

import pandas as pd
import pandas.util.testing as tm

from .base import BaseExtensionTests


class BaseAccumulateTests(BaseExtensionTests):
    """
    Accumulation specific tests. Generally these only
    make sense for numeric/boolean operations.
    """

    def check_accumulate(self, s, op_name, skipna):
        result = getattr(s, op_name)(skipna=skipna)
        expected = getattr(s.astype("float64"), op_name)(skipna=skipna)
        tm.assert_almost_equal(result, expected, check_dtype=False)


class BaseNoAccumulateTests(BaseAccumulateTests):
    """ we don't define any accumulations """

    @pytest.mark.parametrize("skipna", [True, False])
    def test_accumulate_series_numeric(self, data, all_numeric_accumulations, skipna):
        op_name = all_numeric_accumulations
        s = pd.Series(data)

        with pytest.raises(TypeError):
            getattr(s, op_name)(skipna=skipna)


class BaseNumericAccumulateTests(BaseAccumulateTests):
    @pytest.mark.parametrize("skipna", [True, False])
    def test_accumulate_series(self, data, all_numeric_accumulations, skipna):
        op_name = all_numeric_accumulations
        s = pd.Series(data)

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            self.check_accumulate(s, op_name, skipna)
