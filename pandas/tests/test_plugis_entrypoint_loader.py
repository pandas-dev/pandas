import pandas as pd
import pandas._testing as tm
from pandas.core.accessor import DataFrameAccessorLoader


def test_no_accessors(monkeypatch):
    # GH29076

    # Mock entry_points
    def mock_entry_points(*, group):
        return []

    # Patch entry_points in the correct module
    monkeypatch.setattr("pandas.core.accessor.entry_points", mock_entry_points)

    DataFrameAccessorLoader.load()


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


def test_duplicate_accessor_names(monkeypatch):
    # GH29076
    # Create plugin
    class MockEntryPoint1:
        name = "duplicate_accessor"

        def load(self):
            class Accessor1:
                def __init__(self, df):
                    self._df = df

                def which(self):
                    return "Accessor1"

            return Accessor1

    # Create plugin
    class MockEntryPoint2:
        name = "duplicate_accessor"

        def load(self):
            class Accessor2:
                def __init__(self, df):
                    self._df = df

                def which(self):
                    return "Accessor2"

            return Accessor2

    def mock_entry_points(*, group):
        if group == DataFrameAccessorLoader.ENTRY_POINT_GROUP:
            return [MockEntryPoint1(), MockEntryPoint2()]
        return []

    monkeypatch.setattr("pandas.core.accessor.entry_points", mock_entry_points)

    # Check that the UserWarning is raised
    with tm.assert_produces_warning(UserWarning, match="duplicate_accessor") as record:
        DataFrameAccessorLoader.load()

    messages = [str(w.message) for w in record]
    assert any("two packages with the same name" in msg for msg in messages)

    df = pd.DataFrame({"x": [1, 2, 3]})
    assert hasattr(df, "duplicate_accessor")
    assert df.duplicate_accessor.which() in {"Accessor1", "Accessor2"}


def test_unique_accessor_names(monkeypatch):
    # GH29076
    # Create plugin
    class MockEntryPoint1:
        name = "accessor1"

        def load(self):
            class Accessor1:
                def __init__(self, df):
                    self._df = df

                def which(self):
                    return "Accessor1"

            return Accessor1

    # Create plugin
    class MockEntryPoint2:
        name = "accessor2"

        def load(self):
            class Accessor2:
                def __init__(self, df):
                    self._df = df

                def which(self):
                    return "Accessor2"

            return Accessor2

    def mock_entry_points(*, group):
        if group == DataFrameAccessorLoader.ENTRY_POINT_GROUP:
            return [MockEntryPoint1(), MockEntryPoint2()]
        return []

    monkeypatch.setattr("pandas.core.accessor.entry_points", mock_entry_points)

    # Check that no UserWarning is raised
    with tm.assert_produces_warning(None, check_stacklevel=False):
        DataFrameAccessorLoader.load()

    df = pd.DataFrame({"x": [1, 2, 3]})
    assert hasattr(df, "accessor1"), "Accessor1 not registered"
    assert hasattr(df, "accessor2"), "Accessor2 not registered"
    assert df.accessor1.which() == "Accessor1", "Accessor1 method incorrect"
    assert df.accessor2.which() == "Accessor2", "Accessor2 method incorrect"


def test_duplicate_and_unique_accessor_names(monkeypatch):
    # GH29076
    # Create plugin
    class MockEntryPoint1:
        name = "duplicate_accessor"

        def load(self):
            class Accessor1:
                def __init__(self, df):
                    self._df = df

                def which(self):
                    return "Accessor1"

            return Accessor1

    # Create plugin
    class MockEntryPoint2:
        name = "duplicate_accessor"

        def load(self):
            class Accessor2:
                def __init__(self, df):
                    self._df = df

                def which(self):
                    return "Accessor2"

            return Accessor2

    # Create plugin
    class MockEntryPoint3:
        name = "unique_accessor"

        def load(self):
            class Accessor3:
                def __init__(self, df):
                    self._df = df

                def which(self):
                    return "Accessor3"

            return Accessor3

    def mock_entry_points(*, group):
        if group == DataFrameAccessorLoader.ENTRY_POINT_GROUP:
            return [MockEntryPoint1(), MockEntryPoint2(), MockEntryPoint3()]
        return []

    monkeypatch.setattr("pandas.core.accessor.entry_points", mock_entry_points)

    # Capture warnings
    with tm.assert_produces_warning(UserWarning, match="duplicate_accessor") as record:
        DataFrameAccessorLoader.load()

    messages = [str(w.message) for w in record]

    # Filter warnings for the specific message about duplicate packages
    duplicate_package_warnings = [
        msg
        for msg in messages
        if "you have two packages with the same name: 'duplicate_accessor'" in msg
    ]

    # Assert one warning about duplicate packages
    assert len(duplicate_package_warnings) == 1, (
        f"Expected exactly one warning about duplicate packages, "
        f"got {len(duplicate_package_warnings)}: {duplicate_package_warnings}"
    )

    df = pd.DataFrame({"x": [1, 2, 3]})
    assert hasattr(df, "duplicate_accessor"), "duplicate_accessor not registered"

    assert hasattr(df, "unique_accessor"), "unique_accessor not registered"

    assert df.duplicate_accessor.which() in {"Accessor1", "Accessor2"}, (
        "duplicate_accessor method incorrect"
    )
    assert df.unique_accessor.which() == "Accessor3", "unique_accessor method incorrect"
