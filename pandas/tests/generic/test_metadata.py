import pytest

import pandas as pd
from pandas.core.meta import PandasMetadata


class MyMeta(PandasMetadata):
    def finalize(self, new, other, method):
        if method == "concat":
            self.finalize_concat(new, other)
        elif method == "copy":
            self.finalize_copy(new, other)
        else:
            super().finalize(new, other, method)

    def default(self, new, other):
        new.attr = other.attr + 1

    def finalize_concat(self, new, other):
        assert isinstance(other, pd.core.reshape.concat._Concatenator)
        new.attr = sum(x.attr for x in other.objs)


mymeta = MyMeta("attr")


@pytest.fixture
def custom_meta(monkeypatch):
    original_metadata = []

    for cls in [pd.Series, pd.DataFrame]:
        original_metadata.append(cls._metadata)
        custom_metadata = cls._metadata.copy()
        custom_metadata.append("attr")

        monkeypatch.setattr(cls, "_metadata", custom_metadata)


@pytest.mark.usefixtures("custom_meta")
def test_custom_finalizer():

    df = pd.DataFrame({"A": [1, 2]})
    df.attr = 0

    result = df.copy()
    assert result.attr == 1


@pytest.mark.usefixtures("custom_meta")
def test_concat():
    df1 = pd.DataFrame({"A": [1, 2]})
    df1.attr = 2

    df2 = pd.DataFrame({"A": [1, 2]})
    df2.attr = 3

    result = pd.concat([df1, df2])
    assert result.attr == 5
