import pandas as pd
import pandas._testing as tm
from pandas.core.accessor import accessor_entry_point_loader

PANDAS_ENTRY_POINT_GROUP: str = "pandas.accessor"


def create_mock_entry_points(entry_points):
    """
    Auxiliary function to create mock entry points for testing accessor loading.

    Parameters:
    -----------
    entry_points : list of tuple
        List of (name, accessor_class, dist_name) where:
        - name: str, the name of the accessor
        - accessor_class: class, the accessor class to be returned by load()
        - dist_name: str, the name of the distribution (package)

    Returns:
    --------
    function
        A mock_entry_points function that returns the mocked entry points.
    """

    class MockDistribution:
        def __init__(self, name):
            self.name = name

    class MockEntryPoint:
        def __init__(self, name, accessor_class, dist_name):
            self.name = name
            self._accessor_class = accessor_class
            self.dist = MockDistribution(dist_name)

        def load(self):
            return self._accessor_class

    # Create list of MockEntryPoint instances
    mock_eps = [
        MockEntryPoint(name, accessor_class, dist_name)
        for name, accessor_class, dist_name in entry_points
    ]

    def mock_entry_points(*, group):
        if group == PANDAS_ENTRY_POINT_GROUP:
            return mock_eps
        return []

    return mock_entry_points


def test_no_accessors(monkeypatch):
    # No entry points
    mock_entry_points = create_mock_entry_points([])
    monkeypatch.setattr("pandas.core.accessor.entry_points", mock_entry_points)

    accessor_entry_point_loader()


def test_load_dataframe_accessors(monkeypatch):
    class TestAccessor:
        def __init__(self, df):
            self._df = df

        def test_method(self):
            return "success"

    mock_entry_points = create_mock_entry_points(
        [("test_accessor", TestAccessor, "TestPackage")]
    )
    monkeypatch.setattr("pandas.core.accessor.entry_points", mock_entry_points)

    accessor_entry_point_loader()

    # Create DataFrame and verify that the accessor was registered
    df = pd.DataFrame({"a": [1, 2, 3]})
    assert hasattr(df, "test_accessor")
    assert df.test_accessor.test_method() == "success"


def test_duplicate_accessor_names(monkeypatch):
    class Accessor1:
        def __init__(self, df):
            self._df = df

        def which(self):
            return "Accessor1"

    class Accessor2:
        def __init__(self, df):
            self._df = df

        def which(self):
            return "Accessor2"

    mock_entry_points = create_mock_entry_points(
        [
            ("duplicate_accessor", Accessor1, "Package1"),
            ("duplicate_accessor", Accessor2, "Package2"),
        ]
    )
    monkeypatch.setattr("pandas.core.accessor.entry_points", mock_entry_points)

    # Check that the UserWarning is raised
    with tm.assert_produces_warning(UserWarning, match="duplicate_accessor") as record:
        accessor_entry_point_loader()

    messages = [str(w.message) for w in record]
    assert any("you have two accessors with the same name:" in msg for msg in messages)

    df = pd.DataFrame({"x": [1, 2, 3]})
    assert hasattr(df, "duplicate_accessor")
    assert df.duplicate_accessor.which() == "Accessor2"  # Last registered accessor


def test_unique_accessor_names(monkeypatch):
    class Accessor1:
        def __init__(self, df):
            self._df = df

        def which(self):
            return "Accessor1"

    class Accessor2:
        def __init__(self, df):
            self._df = df

        def which(self):
            return "Accessor2"

    mock_entry_points = create_mock_entry_points(
        [("accessor1", Accessor1, "Package1"), ("accessor2", Accessor2, "Package2")]
    )
    monkeypatch.setattr("pandas.core.accessor.entry_points", mock_entry_points)

    # Check that no UserWarning is raised
    with tm.assert_produces_warning(None, check_stacklevel=False):
        accessor_entry_point_loader()

    df = pd.DataFrame({"x": [1, 2, 3]})
    assert hasattr(df, "accessor1"), "Accessor1 not registered"
    assert hasattr(df, "accessor2"), "Accessor2 not registered"
    
    assert df.accessor1.which() == "Accessor1", "Accessor1 method incorrect"
    assert df.accessor2.which() == "Accessor2", "Accessor2 method incorrect"


def test_duplicate_and_unique_accessor_names(monkeypatch):
    class Accessor1:
        def __init__(self, df):
            self._df = df

        def which(self):
            return "Accessor1"

    class Accessor2:
        def __init__(self, df):
            self._df = df

        def which(self):
            return "Accessor2"

    class Accessor3:
        def __init__(self, df):
            self._df = df

        def which(self):
            return "Accessor3"

    mock_entry_points = create_mock_entry_points(
        [
            ("duplicate_accessor", Accessor1, "Package1"),
            ("duplicate_accessor", Accessor2, "Package2"),
            ("unique_accessor", Accessor3, "Package3"),
        ]
    )
    monkeypatch.setattr("pandas.core.accessor.entry_points", mock_entry_points)

    # Capture warnings
    with tm.assert_produces_warning(UserWarning, match="duplicate_accessor") as record:
        accessor_entry_point_loader()

    messages = [str(w.message) for w in record]

    # Filter warnings for the specific message about duplicate accessors
    duplicate_package_warnings = [
        msg
        for msg in messages
        if "you have two accessors with the same name: 'duplicate_accessor'" in msg
    ]

    # Assert one warning about duplicate accessors
    assert len(duplicate_package_warnings) == 1, (
        f"Expected exactly one warning about duplicate accessors, "
        f"got {len(duplicate_package_warnings)}: {duplicate_package_warnings}"
    )

    df = pd.DataFrame({"x": [1, 2, 3]})
    assert hasattr(df, "duplicate_accessor"), "duplicate_accessor not registered"
    assert hasattr(df, "unique_accessor"), "unique_accessor not registered"

    assert df.duplicate_accessor.which() == "Accessor2", (
        "duplicate_accessor should use Accessor2"
    )
    assert df.unique_accessor.which() == "Accessor3", "unique_accessor method incorrect"
