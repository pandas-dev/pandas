import pytest

from pandas.tests.extension.base import BaseInterfaceTests
from pandas.tests.extension.json.array import JSONArray


class JSONArray2(JSONArray):
    @property
    def values(self):
        return self.data


@pytest.fixture
def data():
    return JSONArray2([{"A": [1, 2], "B": [3, 4]}])


class TestBaseInterfaceTests(BaseInterfaceTests):
    def test_no_values_attribute(self, data):
        # GH-20735
        with pytest.raises(ValueError) as m:
            super(TestBaseInterfaceTests, self).test_no_values_attribute(data)
        assert m.match("ExtensionArray contains")
