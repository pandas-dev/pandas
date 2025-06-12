from types import SimpleNamespace

import pytest

import pandas._testing as tm

from pandas.io import common


class _MockIoEngine:
    @classmethod
    def read_foo(cls, fname):
        return "third-party"


@pytest.fixture
def patch_engine(monkeypatch):
    monkeypatch.setattr(common, "_get_io_engine", lambda name: _MockIoEngine)


@pytest.fixture
def patch_entry_points(monkeypatch):
    class MockEntryPoint:
        name = "myengine"
        dist = SimpleNamespace(metadata={"Name": "mypackage"})

        @staticmethod
        def load():
            return _MockIoEngine

    class MockDuplicate1:
        name = "duplicate"
        dist = SimpleNamespace(metadata={"Name": "package1"})

        @staticmethod
        def load():
            return SimpleNamespace(read_foo=lambda fname: "dup1")

    class MockDuplicate2:
        name = "duplicate"
        dist = SimpleNamespace(metadata={"Name": "package2"})

        @staticmethod
        def load():
            return SimpleNamespace(read_foo=lambda fname: "dup1")

    monkeypatch.setattr(common, "_io_engines", None)
    monkeypatch.setattr(
        common,
        "entry_points",
        lambda: SimpleNamespace(
            select=lambda group: [MockEntryPoint, MockDuplicate1, MockDuplicate2]
        ),
    )


class TestIoEngines:
    def test_decorator_with_no_engine(self, patch_engine):
        @common.allow_third_party_engines
        def read_foo(fname, engine=None):
            return "default"

        result = read_foo("myfile.foo")
        assert result == "default"

    def test_decorator_with_skipped_engine(self, patch_engine):
        @common.allow_third_party_engines(skip_engines=["c"])
        def read_foo(fname, engine=None):
            return "default"

        result = read_foo("myfile.foo", engine="c")
        assert result == "default"

    def test_decorator_with_third_party_engine(self, patch_engine):
        @common.allow_third_party_engines
        def read_foo(fname, engine=None):
            return "default"

        result = read_foo("myfile.foo", engine="third-party")
        assert result == "third-party"

    def test_decorator_with_third_party_engine_but_no_method(self, patch_engine):
        @common.allow_third_party_engines
        def read_bar(fname, engine=None):
            return "default"

        msg = "'third-party' does not provide a 'read_bar'"
        with pytest.raises(ValueError, match=msg):
            read_bar("myfile.foo", engine="third-party")

    def test_correct_io_engine(self, patch_entry_points):
        result = common._get_io_engine("myengine")
        assert result is _MockIoEngine

    def test_unknown_io_engine(self, patch_entry_points):
        with pytest.raises(ValueError, match="'unknown' is not a known engine"):
            common._get_io_engine("unknown")

    def test_duplicate_engine(self, patch_entry_points):
        with tm.assert_produces_warning(
            RuntimeWarning,
            match="'duplicate' has been registered by the package 'package1'",
        ):
            result = common._get_io_engine("duplicate")
        assert hasattr(result, "read_foo")
