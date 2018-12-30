import numpy as np
import pytest

from pandas.compat import StringIO

import pandas as pd
import pandas.testing as tm

from .base import BaseExtensionTests


class BaseParsingTests(BaseExtensionTests):

    @pytest.mark.parametrize('engine', ['c', 'python'])
    def test_EA_types(self, engine, data):
        df = pd.DataFrame({
            'with_dtype': pd.Series(data, dtype=str(data.dtype))
        })
        csv_output = df.to_csv(index=False, na_rep=np.nan)
        result = pd.read_csv(StringIO(csv_output), dtype={
            'with_dtype': str(data.dtype)
        }, engine=engine)
        assert result is not None
        tm.assert_frame_equal(df, result)
