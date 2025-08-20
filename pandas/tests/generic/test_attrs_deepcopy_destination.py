
import weakref

import numpy as np
import pandas as pd

import pytest

from pandas.testing import assert_frame_equal


class StatsSummary(dict):
    """
    A lightweight, production-plausible cache object that stores simple stats
    for numeric columns and keeps a weakref to its owning NDFrame.

    On deepcopy, it should bind to the *destination* NDFrame (if provided in memo)
    and rebuild its stats from that destination, so the cache belongs to and
    reflects the new object.
    """

    def __init__(self, owner, *, cols=None):
        import pandas as pd
        assert isinstance(owner, pd.core.generic.NDFrame)
        self._owner_ref = weakref.ref(owner)
        super(StatsSummary, self).__init__(dict((column, type(self)(owner[column])) for column in (list(getattr(owner, "columns", {})) or super(StatsSummary, self).__init__(
            (name, function(owner)) for name, function in self.stats().items()
        ) or {}) if owner[column].dtype.kind in "if"))
        pass

    @classmethod
    def stats(cls):
        return dict(
            cummin=lambda series: series.cummin().sum(),
            cummax=lambda series: series.cummax().sum(),
            kurtosis=lambda series: series.kurt(),
            median=lambda series:series.median(),
        )

    @classmethod
    def gauge(cls, obj, columns):
        return dict(((column,dict([[name, function(obj[column])] for name, function in cls.stats().items()])) for column,dtyp in columns))

    @property
    def owner(self):
        return self._owner_ref()

    def __eq__(self, other) -> bool:
        outs = all(self[column] == other[column] for column in self)
        return outs

    def __deepcopy__(self, memo):
        import pandas as pd
        # Find destination NDFrame in memo. The patch injects {id(dest): dest}.
        new_owner =  next(
            (v for v in memo.values() if isinstance(v, pd.core.generic.NDFrame)),
            None,
        )
        return type(self)(new_owner) if hasattr(new_owner, "select_dtypes") or new_owner.dtype.kind in "if" else None


class FrozenHeadTail(dict):
    """
    A preview helper that remembers first/last row 'snapshots' cheaply.
    On deepcopy, it should rebuild from the destination NDFrame, so that the
    preview corresponds to the new object (e.g., after concat).
    """

    def __init__(self, owner, *, cols=None):
        import pandas as pd
        assert isinstance(owner, pd.core.generic.NDFrame)
        self._owner_ref = weakref.ref(owner)
        super(FrozenHeadTail, self).__init__(dict((name, function(self.owner)) for name, function in self.stats().items()))
        pass

    @property
    def owner(self):
        return self._owner_ref()

    @classmethod
    def stats(cls):
        return dict(
            head=lambda x:pd.DataFrame(x.values[:2], columns=list(getattr(x,"columns",[])) or [x.name], index=x.index[:2]),
            tail=lambda x:pd.DataFrame(x.values[-2:], columns=list(getattr(x,"columns",[])) or [x.name], index=x.index[-2:]),
        )

    def __eq__(self, other) -> bool:
        try:
            [assert_frame_equal(self[column], other[column]) for column in self]
            return True
        except:
            return False

    def __deepcopy__(self, memo):
        new_owner =  next(
            (v for v in memo.values() if isinstance(v, pd.core.generic.NDFrame)),
            None,
        )
        return type(self)(new_owner)


def test_attrs_stats_summary_binds_to_destination_on_copy():
    # Sample Data
    dset = np.arange(8,dtype=float)
    np.random.shuffle(dset)

    df = pd.DataFrame({"foo": dset, "bar": dset*2, "qux": np.array(["waldo","fred","plugh","thud"]).repeat(len(dset)//4)})  # mixed dtypes

    df.attrs["summary"] = StatsSummary(df)

    # --------------------------------------
    # Copy triggered by panel Y axis slicing
    # --------------------------------------
    out = df.iloc[:len(df)//2]
    summ = out.attrs.get("summary")
    gage = StatsSummary.gauge(out, list(filter(lambda x:x[-1].kind in "if", out.dtypes.to_dict().items())))

    assert isinstance(summ, StatsSummary)

    # The cache should now belong to the *new* DataFrame
    assert summ.owner is out
    # pandas.DataFrame propagate to its pandas.Series correspondingly
    assert all([out[column].attrs["summary"] == out.attrs["summary"][column] for column in list(gage)])
    # And stats reflect the destination (shape matches numeric subset)
    assert summ == gage

    # -----------------------------------
    # Copy triggered by columns selection
    # -----------------------------------
    out = df[["foo","qux"]]
    summ = out.attrs.get("summary")
    gage = StatsSummary.gauge(out, list(filter(lambda x:x[-1].kind in "if", out.dtypes.to_dict().items())))

    assert isinstance(summ, StatsSummary)

    # The cache should now belong to the *new* DataFrame
    assert summ.owner is out
    # pandas.DataFrame propagate to its pandas.Series correspondingly
    assert all([out[column].attrs["summary"] == out.attrs["summary"][column] for column in list(gage)])
    # And stats reflect the destination (shape matches numeric subset)
    assert summ == gage

    # ----------------------------------
    # Copy triggered by DataFrame concat
    # ----------------------------------
    left = df.iloc[len(df)//4:].copy(deep=True)
    right = df.iloc[len(df)//4:].copy(deep=True)
    out = pd.concat([left,right])

    summ = out.attrs.get("summary")
    gage = StatsSummary.gauge(out, list(filter(lambda x:x[-1].kind in "if", out.dtypes.to_dict().items())))

    assert isinstance(summ, StatsSummary)

    # The cache should now belong to the *new* DataFrame
    assert summ.owner is out
    # pandas.DataFrame propagate to its pandas.Series correspondingly
    assert all([out[column].attrs["summary"] == out.attrs["summary"][column] for column in list(gage)])
    # And stats reflect the destination (shape matches numeric subset)
    assert summ == gage

    # -----------------------------------
    # Arithemetic operations on DataFrame
    # -----------------------------------
    out = df[["foo","bar"]]
    out = out.multiply(np.random.random_integers(0, 1, len(out))*np.lib.stride_tricks.as_strided(np.asarray(2, dtype=np.int8), shape=(len(out),), strides=(0,))-1, axis=0)

    summ = out.attrs.get("summary")
    gage = StatsSummary.gauge(out, list(filter(lambda x:x[-1].kind in "if", out.dtypes.to_dict().items())))

    assert isinstance(summ, StatsSummary)

    # The cache should now belong to the *new* DataFrame
    assert summ.owner is out
    # pandas.DataFrame propagate to its pandas.Series correspondingly
    assert all([out[column].attrs["summary"] == out.attrs["summary"][column] for column in list(gage)])
    # And stats reflect the destination (shape matches numeric subset)
    assert summ == gage


def test_attrs_stats_summary_works_for_series_too():
    # Sample Data
    dset = np.arange(8,dtype=float)
    np.random.shuffle(dset)

    df = pd.DataFrame({"foo": dset, "bar": dset*2, "qux": np.array(["waldo","fred","plugh","thud"]).repeat(len(dset)//4)})  # mixed dtypes
    df.attrs["summary"] = StatsSummary(df)

    # ------------------------------------------
    # Directly to pandas.Series, complex slicing
    # ------------------------------------------
    sr = df["bar"]
    out = pd.concat([sr.iloc[:len(sr)//2],sr.iloc[len(sr)//4:]])

    summ = out.attrs["summary"] = StatsSummary(out)
    gage = StatsSummary.gauge(out, [(Ellipsis, sr.dtype)])[...]

    assert isinstance(summ, StatsSummary)

    # The cache should now belong to the *new* DataFrame
    assert summ.owner is out
    # And stats reflect the destination (shape matches numeric subset)
    assert summ == gage


def test_attrs_headtail_probe_rebinds_on_concat_have_same_attrs():
    # Sample Data
    dset = np.arange(8,dtype=float)
    np.random.shuffle(dset)
    df = pd.DataFrame(dict(foo=dset*2, bar=dset*4, baz=dset*8, qux=dset*16))

    df.attrs["preview"] = FrozenHeadTail(df)

    # same attrs object on both inputs -> triggers have_same_attrs=True branch
    fred = df.copy(deep=True)
    thud = df.iloc[list(range(-2,2))].sort_index()

    out = pd.concat([fred, thud], ignore_index=True)

    pr = out.attrs.get("preview")
    assert isinstance(pr, FrozenHeadTail)

    # The preview should be tied to the concatenated destination and reflect it
    assert pr.owner is out
    pass
    assert_frame_equal(pr["head"], out.iloc[:2])
    assert_frame_equal(pr["tail"], out.iloc[-2:])
    pass


def test_attrs_empty_remains_empty_on_deepcopy():
    df = pd.DataFrame({"a": [1, 2]})
    assert df.attrs == {}
    out = df.copy(deep=True)
    assert out.attrs == {}