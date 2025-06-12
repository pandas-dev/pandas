import pytest

from pandas.io import common


@pytest.fixture
def patch_engine(monkeypatch):
    class MockIoEngine:
        @classmethod
        def read_foo(cls, fname):
            return "third-party"

    monkeypatch.setattr(common, "_get_io_engine", lambda name: MockIoEngine)


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
