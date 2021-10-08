"""
Tests for Indexes backed by arbitrary ExtensionArrays.
"""
import pandas as pd
from pandas.tests.extension.base.base import BaseExtensionTests


class BaseExtensionIndexTests(BaseExtensionTests):
    def test_index_from_array(self, data):
        idx = pd.Index(data)
        assert data.dtype == idx.dtype
