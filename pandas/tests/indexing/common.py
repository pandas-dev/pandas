""" common utilities """
import itertools
from warnings import catch_warnings

import numpy as np

from pandas.core.dtypes.common import is_scalar

from pandas import DataFrame, Float64Index, MultiIndex, Series, UInt64Index, date_range
import pandas._testing as tm


def _mklbl(prefix, n):
    return ["{prefix}{i}".format(prefix=prefix, i=i) for i in range(n)]


def _axify(obj, key, axis):
    # create a tuple accessor
    axes = [slice(None)] * obj.ndim
    axes[axis] = key
    return tuple(axes)


class Base:
    """ indexing comprehensive base class """

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

    def setup_method(self, method):

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
            d = dict()
            for typ in self._typs:
                d[typ] = getattr(self, "{kind}_{typ}".format(kind=kind, typ=typ))

            setattr(self, kind, d)

    def generate_indices(self, f, values=False):
        """ generate the indices
        if values is True , use the axis values
        is False, use the range
        """

        axes = f.axes
        if values:
            axes = (list(range(len(ax))) for ax in axes)

        return itertools.product(*axes)

    def get_result(self, obj, method, key, axis):
        """ return the result for this obj with this key and this axis """

        if isinstance(key, dict):
            key = key[axis]

        # use an artificial conversion to map the key as integers to the labels
        # so ix can work for comparisons
        if method == "indexer":
            method = "ix"
            key = obj._get_axis(axis)[key]

        # in case we actually want 0 index slicing
        with catch_warnings(record=True):
            try:
                xp = getattr(obj, method).__getitem__(_axify(obj, key, axis))
            except AttributeError:
                xp = getattr(obj, method).__getitem__(key)

        return xp

    def get_value(self, name, f, i, values=False):
        """ return the value for the location i """

        # check against values
        if values:
            return f.values[i]

        elif name == "iat":
            return f.iloc[i]
        else:
            assert name == "at"
            return f.loc[i]

    def check_values(self, f, func, values=False):

        if f is None:
            return
        axes = f.axes
        indicies = itertools.product(*axes)

        for i in indicies:
            result = getattr(f, func)[i]

            # check against values
            if values:
                expected = f.values[i]
            else:
                expected = f
                for a in reversed(i):
                    expected = expected.__getitem__(a)

            tm.assert_almost_equal(result, expected)

    def check_result(
        self, method1, key1, method2, key2, typs=None, axes=None, fails=None,
    ):
        def _eq(axis, obj, key1, key2):
            """ compare equal for these 2 keys """
            if axis > obj.ndim - 1:
                return

            try:
                rs = getattr(obj, method1).__getitem__(_axify(obj, key1, axis))

                try:
                    xp = self.get_result(obj=obj, method=method2, key=key2, axis=axis)
                except (KeyError, IndexError):
                    # TODO: why is this allowed?
                    return

                if is_scalar(rs) and is_scalar(xp):
                    assert rs == xp
                else:
                    tm.assert_equal(rs, xp)

            except (IndexError, TypeError, KeyError) as detail:

                # if we are in fails, the ok, otherwise raise it
                if fails is not None:
                    if isinstance(detail, fails):
                        result = f"ok ({type(detail).__name__})"
                        return

                result = type(detail).__name__
                raise AssertionError(result, detail)

        if typs is None:
            typs = self._typs

        if axes is None:
            axes = [0, 1]
        elif not isinstance(axes, (tuple, list)):
            assert isinstance(axes, int)
            axes = [axes]

        # check
        for kind in self._kinds:

            d = getattr(self, kind)
            for ax in axes:
                for typ in typs:
                    if typ not in self._typs:
                        continue

                    obj = d[typ]
                    _eq(axis=ax, obj=obj, key1=key1, key2=key2)
