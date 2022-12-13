import pytest

import pandas as pd
from pandas.tests.extension.base.base import BaseExtensionTests


class BaseAccumulateTests(BaseExtensionTests):
    """
    Accumulation specific tests. Generally these only
    make sense for numeric/boolean operations.
    """

    def check_accumulate(self, s, op_name, skipna):
        result = getattr(s, op_name)(skipna=skipna)

        if result.dtype == pd.Float32Dtype() and op_name == "cumprod" and skipna:
            pytest.skip(
                f"Float32 precision lead to large differences with op {op_name} "
                f"and skipna={skipna}"
            )

        expected = getattr(s.astype("float64"), op_name)(skipna=skipna)
        self.assert_series_equal(result, expected, check_dtype=False)

    @pytest.mark.parametrize("skipna", [True, False])
    def test_accumulate_series_raises(self, data, all_numeric_accumulations, skipna):
        op_name = all_numeric_accumulations
        ser = pd.Series(data)

        with pytest.raises(NotImplementedError):
            getattr(ser, op_name)(skipna=skipna)

    @pytest.mark.parametrize("skipna", [True, False])
    def test_accumulate_series(self, data, all_numeric_accumulations, skipna):
        op_name = all_numeric_accumulations
        ser = pd.Series(data)
        self.check_accumulate(ser, op_name, skipna)
