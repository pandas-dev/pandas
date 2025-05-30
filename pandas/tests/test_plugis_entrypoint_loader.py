import pandas as pd
from pandas.core.accessor import DataFrameAccessorLoader


def test_load_dataframe_accessors(monkeypatch):
    # GH29076
    # Mocked EntryPoint to simulate a plugin
    class MockEntryPoint:
        name = "test_accessor"

        def load(self):
            class TestAccessor:
                def __init__(self, df):
                    self._df = df

                def test_method(self):
                    return "success"

            return TestAccessor

    # Mock entry_points
    def mock_entry_points(*, group):
        if group == DataFrameAccessorLoader.ENTRY_POINT_GROUP:
            return [MockEntryPoint()]
        return []

    # Patch entry_points in the correct module
    monkeypatch.setattr("pandas.core.accessor.entry_points", mock_entry_points)

    DataFrameAccessorLoader.load()

    # Create DataFrame and verify that the accessor was registered
    df = pd.DataFrame({"a": [1, 2, 3]})
    assert hasattr(df, "test_accessor")
    assert df.test_accessor.test_method() == "success"
