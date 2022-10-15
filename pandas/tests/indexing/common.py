""" common utilities """
import numpy as np

from pandas import (
    DataFrame,
    MultiIndex,
    Series,
    date_range,
)
from pandas.core.api import (
    Float64Index,
    UInt64Index,
)


def _mklbl(prefix, n):
    return [f"{prefix}{i}" for i in range(n)]


def check_result(obj, method, key, axes=None, fails=None):
    if axes is None:
        axes = [0, 1]
    else:
        assert axes in [0, 1]
        axes = [axes]

    for ax in axes:
        if ax < obj.ndim:
            # create a tuple accessor
            axes = [slice(None)] * obj.ndim
            axes[ax] = key
            axified = tuple(axes)
            try:
                getattr(obj, method).__getitem__(axified)
            except (IndexError, TypeError, KeyError) as detail:
                # if we are in fails, the ok, otherwise raise it
                if fails is not None:
                    if isinstance(detail, fails):
                        return
                raise


class Base:
    """indexing comprehensive base class"""

    _kinds = {"series", "frame"}
    _typs = {
        "ints",
        "uints",
        "labels",
        "mixed",
        "ts",
        "floats",
        "empty",
        "ts_rev",
        "multi",
    }

    def setup_method(self):

        self.series_ints = Series(np.random.rand(4), index=np.arange(0, 8, 2))
        self.frame_ints = DataFrame(
            np.random.randn(4, 4), index=np.arange(0, 8, 2), columns=np.arange(0, 12, 3)
        )

        self.series_uints = Series(
            np.random.rand(4), index=UInt64Index(np.arange(0, 8, 2))
        )
        self.frame_uints = DataFrame(
            np.random.randn(4, 4),
            index=UInt64Index(range(0, 8, 2)),
            columns=UInt64Index(range(0, 12, 3)),
        )

        self.series_floats = Series(
            np.random.rand(4), index=Float64Index(range(0, 8, 2))
        )
        self.frame_floats = DataFrame(
            np.random.randn(4, 4),
            index=Float64Index(range(0, 8, 2)),
            columns=Float64Index(range(0, 12, 3)),
        )

        m_idces = [
            MultiIndex.from_product([[1, 2], [3, 4]]),
            MultiIndex.from_product([[5, 6], [7, 8]]),
            MultiIndex.from_product([[9, 10], [11, 12]]),
        ]

        self.series_multi = Series(np.random.rand(4), index=m_idces[0])
        self.frame_multi = DataFrame(
            np.random.randn(4, 4), index=m_idces[0], columns=m_idces[1]
        )

        self.series_labels = Series(np.random.randn(4), index=list("abcd"))
        self.frame_labels = DataFrame(
            np.random.randn(4, 4), index=list("abcd"), columns=list("ABCD")
        )

        self.series_mixed = Series(np.random.randn(4), index=[2, 4, "null", 8])
        self.frame_mixed = DataFrame(np.random.randn(4, 4), index=[2, 4, "null", 8])

        self.series_ts = Series(
            np.random.randn(4), index=date_range("20130101", periods=4)
        )
        self.frame_ts = DataFrame(
            np.random.randn(4, 4), index=date_range("20130101", periods=4)
        )

        dates_rev = date_range("20130101", periods=4).sort_values(ascending=False)
        self.series_ts_rev = Series(np.random.randn(4), index=dates_rev)
        self.frame_ts_rev = DataFrame(np.random.randn(4, 4), index=dates_rev)

        self.frame_empty = DataFrame()
        self.series_empty = Series(dtype=object)

        # form agglomerates
        for kind in self._kinds:
            d = {}
            for typ in self._typs:
                d[typ] = getattr(self, f"{kind}_{typ}")

            setattr(self, kind, d)
