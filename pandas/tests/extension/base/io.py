import pandas as pd
from pandas.compat import StringIO
from pandas.core.arrays.integer import Int64Dtype
from .base import BaseExtensionTests


class ExtensionParsingTests(BaseExtensionTests):
    def test_EA_types(self):
        df = pd.DataFrame({'Int': pd.Series([1, 2, 3], dtype='Int64'),
                           'A': [1, 2, 1]})
        data = df.to_csv(index=False)
        result = pd.read_csv(StringIO(data), dtype={'Int': Int64Dtype})
        assert result is not None

        df = pd.DataFrame({'Int': pd.Series([1, 2, 3], dtype='Int8'),
                           'A': [1, 2, 1]})
        data = df.to_csv(index=False)
        result = pd.read_csv(StringIO(data), dtype={'Int': 'Int8'})
        assert result is not None
