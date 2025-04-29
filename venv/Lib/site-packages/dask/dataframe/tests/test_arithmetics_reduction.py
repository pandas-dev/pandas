from __future__ import annotations

import contextlib
import warnings
from datetime import datetime

import numpy as np
import pandas as pd
import pytest
from pandas.api.types import is_scalar

import dask.dataframe as dd
from dask.array.numpy_compat import NUMPY_GE_125, NUMPY_GE_200
from dask.dataframe._compat import PANDAS_GE_300
from dask.dataframe.utils import assert_eq, pyarrow_strings_enabled

try:
    import scipy
except ImportError:
    scipy = None

try:
    import pyarrow as pa
    from pyarrow.lib import ArrowNotImplementedError
except ImportError:
    pa = None
    ArrowNotImplementedError = None


def test_deterministic_arithmetic_names():
    df = pd.DataFrame({"x": [1, 2, 3, 4], "y": [5, 6, 7, 8]})
    a = dd.from_pandas(df, npartitions=2)

    assert sorted((a.x + a.y**2).dask) == sorted((a.x + a.y**2).dask)
    assert sorted((a.x + a.y**2).dask) != sorted((a.x + a.y**3).dask)
    assert sorted((a.x + a.y**2).dask) != sorted((a.x - a.y**2).dask)


@pytest.mark.slow
def test_arithmetics_different_index():
    # index are different, but overwraps
    pdf1 = pd.DataFrame(
        {"a": [1, 2, 3, 4, 5], "b": [3, 5, 2, 5, 7]}, index=[1, 2, 3, 4, 5]
    )
    ddf1 = dd.from_pandas(pdf1, 2)
    pdf2 = pd.DataFrame(
        {"a": [3, 2, 6, 7, 8], "b": [9, 4, 2, 6, 2]}, index=[3, 4, 5, 6, 7]
    )
    ddf2 = dd.from_pandas(pdf2, 2)

    # index are not overwrapped
    pdf3 = pd.DataFrame(
        {"a": [1, 2, 3, 4, 5], "b": [3, 5, 2, 5, 7]}, index=[1, 2, 3, 4, 5]
    )
    ddf3 = dd.from_pandas(pdf3, 2)
    pdf4 = pd.DataFrame(
        {"a": [3, 2, 6, 7, 8], "b": [9, 4, 2, 6, 2]}, index=[10, 11, 12, 13, 14]
    )
    ddf4 = dd.from_pandas(pdf4, 2)

    # index is included in another
    pdf5 = pd.DataFrame(
        {"a": [1, 2, 3, 4, 5], "b": [3, 5, 2, 5, 7]}, index=[1, 3, 5, 7, 9]
    )
    ddf5 = dd.from_pandas(pdf5, 2)
    pdf6 = pd.DataFrame(
        {"a": [3, 2, 6, 7, 8], "b": [9, 4, 2, 6, 2]}, index=[2, 3, 4, 5, 6]
    )
    ddf6 = dd.from_pandas(pdf6, 2)

    cases = [
        (ddf1, ddf2, pdf1, pdf2),
        (ddf2, ddf1, pdf2, pdf1),
        (ddf1.repartition([1, 3, 5]), ddf2.repartition([3, 4, 7]), pdf1, pdf2),
        (ddf2.repartition([3, 4, 5, 7]), ddf1.repartition([1, 2, 4, 5]), pdf2, pdf1),
        (ddf3, ddf4, pdf3, pdf4),
        (ddf4, ddf3, pdf4, pdf3),
        (
            ddf3.repartition([1, 2, 3, 4, 5]),
            ddf4.repartition([10, 11, 12, 13, 14]),
            pdf3,
            pdf4,
        ),
        (ddf4.repartition([10, 14]), ddf3.repartition([1, 3, 4, 5]), pdf4, pdf3),
        (ddf5, ddf6, pdf5, pdf6),
        (ddf6, ddf5, pdf6, pdf5),
        (ddf5.repartition([1, 7, 8, 9]), ddf6.repartition([2, 3, 4, 6]), pdf5, pdf6),
        (ddf6.repartition([2, 6]), ddf5.repartition([1, 3, 7, 9]), pdf6, pdf5),
        # dask + pandas
        (ddf1, pdf2, pdf1, pdf2),
        (ddf2, pdf1, pdf2, pdf1),
        (ddf3, pdf4, pdf3, pdf4),
        (ddf4, pdf3, pdf4, pdf3),
        (ddf5, pdf6, pdf5, pdf6),
        (ddf6, pdf5, pdf6, pdf5),
    ]

    for l, r, el, er in cases:
        check_series_arithmetics(l.a, r.b, el.a, er.b, allow_comparison_ops=False)
        check_frame_arithmetics(l, r, el, er, allow_comparison_ops=False)

    pdf7 = pd.DataFrame(
        {"a": [1, 2, 3, 4, 5, 6, 7, 8], "b": [5, 6, 7, 8, 1, 2, 3, 4]},
        index=[0, 2, 4, 8, 9, 10, 11, 13],
    )
    pdf8 = pd.DataFrame(
        {"a": [5, 6, 7, 8, 4, 3, 2, 1], "b": [2, 4, 5, 3, 4, 2, 1, 0]},
        index=[1, 3, 4, 8, 9, 11, 12, 13],
    )
    ddf7 = dd.from_pandas(pdf7, 3)
    ddf8 = dd.from_pandas(pdf8, 2)

    pdf9 = pd.DataFrame(
        {"a": [1, 2, 3, 4, 5, 6, 7, 8], "b": [5, 6, 7, 8, 1, 2, 3, 4]},
        index=[0, 2, 4, 8, 9, 10, 11, 13],
    )
    pdf10 = pd.DataFrame(
        {"a": [5, 6, 7, 8, 4, 3, 2, 1], "b": [2, 4, 5, 3, 4, 2, 1, 0]},
        index=[0, 3, 4, 8, 9, 11, 12, 13],
    )
    ddf9 = dd.from_pandas(pdf9, 3)
    ddf10 = dd.from_pandas(pdf10, 2)

    cases = [
        (ddf7, ddf8, pdf7, pdf8),
        (ddf8, ddf7, pdf8, pdf7),
        # (ddf7.repartition([0, 13]),
        #  ddf8.repartition([0, 4, 11, 14], force=True),
        #  pdf7, pdf8),
        (
            ddf8.repartition([-5, 10, 15], force=True),
            ddf7.repartition([-1, 4, 11, 14], force=True),
            pdf8,
            pdf7,
        ),
        (
            ddf7.repartition([0, 8, 12, 13]),
            ddf8.repartition([0, 2, 8, 12, 13], force=True),
            pdf7,
            pdf8,
        ),
        (
            ddf8.repartition([-5, 0, 10, 20], force=True),
            ddf7.repartition([-1, 4, 11, 13], force=True),
            pdf8,
            pdf7,
        ),
        (ddf9, ddf10, pdf9, pdf10),
        (ddf10, ddf9, pdf10, pdf9),
        # dask + pandas
        (ddf7, pdf8, pdf7, pdf8),
        (ddf8, pdf7, pdf8, pdf7),
        (ddf9, pdf10, pdf9, pdf10),
        (ddf10, pdf9, pdf10, pdf9),
    ]

    for l, r, el, er in cases:
        check_series_arithmetics(l.a, r.b, el.a, er.b, allow_comparison_ops=False)
        check_frame_arithmetics(l, r, el, er, allow_comparison_ops=False)


def check_series_arithmetics(l, r, el, er, allow_comparison_ops=True):
    assert isinstance(l, dd.Series)
    assert isinstance(r, (dd.Series, pd.Series))
    assert isinstance(el, pd.Series)
    assert isinstance(er, pd.Series)

    # l, r may be repartitioned, test whether repartition keeps original data
    assert_eq(l, el)
    assert_eq(r, er)

    assert_eq(l + r, el + er)
    assert_eq(l * r, el * er)
    assert_eq(l - r, el - er)
    assert_eq(l / r, el / er)
    assert_eq(l // r, el // er)
    assert_eq(l**r, el**er)
    assert_eq(l % r, el % er)

    if allow_comparison_ops:
        # comparison is allowed if data have same index
        assert_eq(l & r, el & er)
        assert_eq(l | r, el | er)
        assert_eq(l ^ r, el ^ er)
        assert_eq(l > r, el > er)
        assert_eq(l < r, el < er)
        assert_eq(l >= r, el >= er)
        assert_eq(l <= r, el <= er)
        assert_eq(l == r, el == er)
        assert_eq(l != r, el != er)
        assert_eq(l.lt(r), el.lt(er))
        assert_eq(l.gt(r), el.gt(er))
        assert_eq(l.le(r), el.le(er))
        assert_eq(l.ge(r), el.ge(er))
        assert_eq(l.ne(r), el.ne(er))
        assert_eq(l.eq(r), el.eq(er))

    assert_eq(l + 2, el + 2)
    assert_eq(l * 2, el * 2)
    assert_eq(l - 2, el - 2)
    assert_eq(l / 2, el / 2)
    assert_eq(l & True, el & True)
    assert_eq(l | True, el | True)
    assert_eq(l ^ True, el ^ True)
    assert_eq(l // 2, el // 2)
    assert_eq(l**2, el**2)
    assert_eq(l % 2, el % 2)
    assert_eq(l > 2, el > 2)
    assert_eq(l < 2, el < 2)
    assert_eq(l >= 2, el >= 2)
    assert_eq(l <= 2, el <= 2)
    assert_eq(l == 2, el == 2)
    assert_eq(l != 2, el != 2)

    assert_eq(2 + r, 2 + er)
    assert_eq(2 * r, 2 * er)
    assert_eq(2 - r, 2 - er)
    assert_eq(2 / r, 2 / er)
    assert_eq(True & r, True & er)
    assert_eq(True | r, True | er)
    assert_eq(True ^ r, True ^ er)
    assert_eq(2 // r, 2 // er)
    assert_eq(2**r, 2**er)
    assert_eq(2 % r, 2 % er)
    assert_eq(2 > r, 2 > er)
    assert_eq(2 < r, 2 < er)
    assert_eq(2 >= r, 2 >= er)
    assert_eq(2 <= r, 2 <= er)
    assert_eq(2 == r, 2 == er)
    assert_eq(2 != r, 2 != er)

    assert_eq(l.lt(2), el.lt(2))
    assert_eq(l.gt(2), el.gt(2))
    assert_eq(l.le(2), el.le(2))
    assert_eq(l.ge(2), el.ge(2))
    assert_eq(l.ne(2), el.ne(2))
    assert_eq(l.eq(2), el.eq(2))

    assert_eq(-l, -el)
    assert_eq(abs(l), abs(el))

    if allow_comparison_ops:
        # comparison is allowed if data have same index
        assert_eq(~(l == r), ~(el == er))


def check_frame_arithmetics(l, r, el, er, allow_comparison_ops=True):
    assert isinstance(l, dd.DataFrame)
    assert isinstance(r, (dd.DataFrame, pd.DataFrame))
    assert isinstance(el, pd.DataFrame)
    assert isinstance(er, pd.DataFrame)
    # l, r may be repartitioned, test whether repartition keeps original data
    assert_eq(l, el)
    assert_eq(r, er)

    assert_eq(l + r, el + er)
    assert_eq(l * r, el * er)
    assert_eq(l - r, el - er)
    assert_eq(l / r, el / er)
    assert_eq(l // r, el // er)
    assert_eq(l**r, el**er)
    assert_eq(l % r, el % er)

    if allow_comparison_ops:
        # comparison is allowed if data have same index
        assert_eq(l & r, el & er)
        assert_eq(l | r, el | er)
        assert_eq(l ^ r, el ^ er)
        assert_eq(l > r, el > er)
        assert_eq(l < r, el < er)
        assert_eq(l >= r, el >= er)
        assert_eq(l <= r, el <= er)
        assert_eq(l == r, el == er)
        assert_eq(l != r, el != er)
        assert_eq(l.lt(r), el.lt(er))
        assert_eq(l.gt(r), el.gt(er))
        assert_eq(l.le(r), el.le(er))
        assert_eq(l.ge(r), el.ge(er))
        assert_eq(l.ne(r), el.ne(er))
        assert_eq(l.eq(r), el.eq(er))

    assert_eq(l + 2, el + 2)
    assert_eq(l * 2, el * 2)
    assert_eq(l - 2, el - 2)
    assert_eq(l / 2, el / 2)
    assert_eq(l & True, el & True)
    assert_eq(l | True, el | True)
    assert_eq(l ^ True, el ^ True)
    assert_eq(l // 2, el // 2)
    assert_eq(l**2, el**2)
    assert_eq(l % 2, el % 2)
    assert_eq(l > 2, el > 2)
    assert_eq(l < 2, el < 2)
    assert_eq(l >= 2, el >= 2)
    assert_eq(l <= 2, el <= 2)
    assert_eq(l == 2, el == 2)
    assert_eq(l != 2, el != 2)

    assert_eq(2 + l, 2 + el)
    assert_eq(2 * l, 2 * el)
    assert_eq(2 - l, 2 - el)
    assert_eq(2 / l, 2 / el)
    assert_eq(True & l, True & el)
    assert_eq(True | l, True | el)
    assert_eq(True ^ l, True ^ el)
    assert_eq(2 // l, 2 // el)
    assert_eq(2**l, 2**el)
    assert_eq(2 % l, 2 % el)
    assert_eq(2 > l, 2 > el)
    assert_eq(2 < l, 2 < el)
    assert_eq(2 >= l, 2 >= el)
    assert_eq(2 <= l, 2 <= el)
    assert_eq(2 == l, 2 == el)
    assert_eq(2 != l, 2 != el)

    assert_eq(l.lt(2), el.lt(2))
    assert_eq(l.gt(2), el.gt(2))
    assert_eq(l.le(2), el.le(2))
    assert_eq(l.ge(2), el.ge(2))
    assert_eq(l.ne(2), el.ne(2))
    assert_eq(l.eq(2), el.eq(2))

    assert_eq(-l, -el)
    assert_eq(abs(l), abs(el))

    if allow_comparison_ops:
        # comparison is allowed if data have same index
        assert_eq(~(l == r), ~(el == er))


def test_frame_series_arithmetic_methods():
    pdf1 = pd.DataFrame(
        {
            "A": np.arange(10),
            "B": [np.nan, 1, 2, 3, 4] * 2,
            "C": [np.nan] * 10,
            "D": np.arange(10),
        },
        index=list("abcdefghij"),
        columns=list("ABCD"),
    )
    pdf2 = pd.DataFrame(
        np.random.randn(10, 4), index=list("abcdefghjk"), columns=list("ABCX")
    )
    ps1 = pdf1.A
    ps2 = pdf2.A
    ps3 = pd.Series(np.random.randn(10), index=list("ABCDXabcde"))

    ddf1 = dd.from_pandas(pdf1, 2)
    ddf2 = dd.from_pandas(pdf2, 2)
    ds1 = ddf1.A
    ds2 = ddf2.A

    s = 4

    for l, r, el, er in [
        (ddf1, ddf2, pdf1, pdf2),
        (ds1, ds2, ps1, ps2),
        (ddf1.repartition(["a", "f", "j"]), ddf2, pdf1, pdf2),
        (ds1.repartition(["a", "b", "f", "j"]), ds2, ps1, ps2),
        (ddf1, ddf2.repartition(["a", "k"]), pdf1, pdf2),
        (ds1, ds2.repartition(["a", "b", "d", "h", "k"]), ps1, ps2),
        (ddf1, 3, pdf1, 3),
        (ds1, 3, ps1, 3),
        (ddf1, s, pdf1, 4),
        (ds1, s, ps1, 4),
    ]:
        # l, r may be repartitioned, test whether repartition keeps original data
        assert_eq(l, el)
        assert_eq(r, er)

        assert_eq(l.add(r, fill_value=0), el.add(er, fill_value=0))
        assert_eq(l.sub(r, fill_value=0), el.sub(er, fill_value=0))
        assert_eq(l.mul(r, fill_value=0), el.mul(er, fill_value=0))
        assert_eq(l.div(r, fill_value=0), el.div(er, fill_value=0))
        assert_eq(l.divide(r, fill_value=0), el.divide(er, fill_value=0))
        assert_eq(l.truediv(r, fill_value=0), el.truediv(er, fill_value=0))
        assert_eq(l.floordiv(r, fill_value=1), el.floordiv(er, fill_value=1))
        assert_eq(l.pow(r, fill_value=0), el.pow(er, fill_value=0))
        assert_eq(l.mod(r, fill_value=0), el.mod(er, fill_value=0))

        assert_eq(l.radd(r, fill_value=0), el.radd(er, fill_value=0))
        assert_eq(l.rsub(r, fill_value=0), el.rsub(er, fill_value=0))
        assert_eq(l.rmul(r, fill_value=0), el.rmul(er, fill_value=0))
        assert_eq(l.rdiv(r, fill_value=0), el.rdiv(er, fill_value=0))
        assert_eq(l.rtruediv(r, fill_value=0), el.rtruediv(er, fill_value=0))
        assert_eq(l.rpow(r, fill_value=0), el.rpow(er, fill_value=0))
        assert_eq(l.rmod(r, fill_value=0), el.rmod(er, fill_value=0))

    for l, r, el, er in [(ddf1, ds2, pdf1, ps2), (ddf1, ddf2.X, pdf1, pdf2.X)]:
        assert_eq(l, el)
        assert_eq(r, er)

        # must specify axis=0 to add Series to each column
        # axis=1 is not supported (add to each row)
        assert_eq(l.add(r, axis=0), el.add(er, axis=0))
        assert_eq(l.sub(r, axis=0), el.sub(er, axis=0))
        assert_eq(l.mul(r, axis=0), el.mul(er, axis=0))
        assert_eq(l.div(r, axis=0), el.div(er, axis=0))
        assert_eq(l.divide(r, axis=0), el.divide(er, axis=0))
        assert_eq(l.truediv(r, axis=0), el.truediv(er, axis=0))
        assert_eq(l.floordiv(r, axis=0), el.floordiv(er, axis=0))
        assert_eq(l.mod(r, axis=0), el.mod(er, axis=0))
        assert_eq(l.pow(r, axis=0), el.pow(er, axis=0))

        assert_eq(l.radd(r, axis=0), el.radd(er, axis=0))
        assert_eq(l.rsub(r, axis=0), el.rsub(er, axis=0))
        assert_eq(l.rmul(r, axis=0), el.rmul(er, axis=0))
        assert_eq(l.rdiv(r, axis=0), el.rdiv(er, axis=0))
        assert_eq(l.rtruediv(r, axis=0), el.rtruediv(er, axis=0))
        assert_eq(l.rmod(r, axis=0), el.rmod(er, axis=0))
        assert_eq(l.rpow(r, axis=0), el.rpow(er, axis=0))

        pytest.raises(ValueError, lambda l=l, r=r: l.add(r, axis=1))

    for l, r, el, er in [(ddf1, pdf2, pdf1, pdf2), (ddf1, ps3, pdf1, ps3)]:
        assert_eq(l, el)
        assert_eq(r, er)

        for axis in [0, 1, "index", "columns"]:
            assert_eq(l.add(r, axis=axis), el.add(er, axis=axis))
            assert_eq(l.sub(r, axis=axis), el.sub(er, axis=axis))
            assert_eq(l.mul(r, axis=axis), el.mul(er, axis=axis))
            assert_eq(l.div(r, axis=axis), el.div(er, axis=axis))
            assert_eq(l.divide(r, axis=axis), el.divide(er, axis=axis))
            assert_eq(l.truediv(r, axis=axis), el.truediv(er, axis=axis))
            assert_eq(l.floordiv(r, axis=axis), el.floordiv(er, axis=axis))
            assert_eq(l.mod(r, axis=axis), el.mod(er, axis=axis))
            assert_eq(l.pow(r, axis=axis), el.pow(er, axis=axis))
            assert_eq(l.rdiv(r, axis=axis), el.rdiv(er, axis=axis))
            assert_eq(l.rtruediv(r, axis=axis), el.rtruediv(er, axis=axis))
            assert_eq(l.rpow(r, axis=axis), el.rpow(er, axis=axis))
            assert_eq(l.rmod(r, axis=axis), el.rmod(er, axis=axis))
            assert_eq(l.radd(r, axis=axis), el.radd(er, axis=axis))
            assert_eq(l.rsub(r, axis=axis), el.rsub(er, axis=axis))
            assert_eq(l.rmul(r, axis=axis), el.rmul(er, axis=axis))


@pytest.mark.parametrize("split_every", [False, 2])
def test_reductions(split_every):
    dsk = {
        ("x", 0): pd.DataFrame(
            {"a": [1, 2, 3], "b": [4, 5, 6], "c": [True, True, False]}, index=[0, 1, 3]
        ),
        ("x", 1): pd.DataFrame(
            {"a": [4, 5, 6], "b": [3, 2, 1], "c": [False, False, False]},
            index=[5, 6, 8],
        ),
        ("x", 2): pd.DataFrame(
            {
                "a": [13094304034, 3489385935, 100006774],
                "b": [0, 0, 0],
                "c": [True, True, True],
            },
            index=[9, 9, 9],
        ),
    }
    ddf1 = dd.repartition(pd.concat(dsk.values()), divisions=[0, 4, 9, 9])
    pdf1 = ddf1.compute()

    nans1 = pd.Series([1] + [np.nan] * 4 + [2] + [np.nan] * 3)
    nands1 = dd.from_pandas(nans1, 2)
    nans2 = pd.Series([1] + [np.nan] * 8)
    nands2 = dd.from_pandas(nans2, 2)
    nans3 = pd.Series([np.nan] * 9)
    nands3 = dd.from_pandas(nans3, 2)

    bools = pd.Series([True, False, True, False, True], dtype=bool)
    boolds = dd.from_pandas(bools, 2)

    for dds, pds in [
        (ddf1.a, pdf1.a),
        (ddf1.b, pdf1.b),
        (ddf1.c, pdf1.c),
        (ddf1["a"], pdf1["a"]),
        (ddf1["b"], pdf1["b"]),
        (nands1, nans1),
        (nands2, nans2),
        (nands3, nans3),
        (boolds, bools),
    ]:
        assert isinstance(dds, dd.Series)
        assert isinstance(pds, pd.Series)

        assert_eq(dds.sum(split_every=split_every), pds.sum())
        assert_eq(dds.prod(split_every=split_every), pds.prod())
        assert_eq(dds.product(split_every=split_every), pds.product())
        assert_eq(dds.min(split_every=split_every), pds.min())
        assert_eq(dds.max(split_every=split_every), pds.max())
        assert_eq(dds.count(split_every=split_every), pds.count())

        if scipy:
            # pandas uses unbiased skew, need to correct for that
            n = pds.shape[0]
            bias_factor = (n * (n - 1)) ** 0.5 / (n - 2)
            assert_eq(dds.skew(), pds.skew() / bias_factor)

            # TODO: Remove this `if`-block once `axis=None` support is added.
            # https://github.com/dask/dask/issues/9915
            with pytest.raises(
                ValueError, match="`axis=None` isn't currently supported"
            ):
                dds.skew(axis=None)

        if scipy:
            # pandas uses a bias factor for kurtosis, need to correct for that
            n = pds.shape[0]
            factor = ((n - 1) * (n + 1)) / ((n - 2) * (n - 3))
            offset = (6 * (n - 1)) / ((n - 2) * (n - 3))
            assert_eq(factor * dds.kurtosis() + offset, pds.kurtosis())

            # TODO: Remove this `if`-block once `axis=None` support is added.
            # https://github.com/dask/dask/issues/9915
            with pytest.raises(
                ValueError, match="`axis=None` isn't currently supported"
            ):
                dds.kurtosis(axis=None)

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            # runtime warnings; https://github.com/dask/dask/issues/2381
            assert_eq(dds.std(split_every=split_every), pds.std())
            assert_eq(dds.var(split_every=split_every), pds.var())
            assert_eq(dds.sem(split_every=split_every), pds.sem())

        with warnings.catch_warnings():
            # dask.dataframe should probably filter this, to match pandas, but
            # it seems quite difficult.
            warnings.simplefilter("ignore", RuntimeWarning)
            assert_eq(dds.std(ddof=0, split_every=split_every), pds.std(ddof=0))
            assert_eq(dds.var(ddof=0, split_every=split_every), pds.var(ddof=0))
            assert_eq(dds.sem(ddof=0, split_every=split_every), pds.sem(ddof=0))
        assert_eq(dds.mean(split_every=split_every), pds.mean())
        assert_eq(dds.nunique(split_every=split_every), pds.nunique())

        assert_eq(dds.sum(skipna=False, split_every=split_every), pds.sum(skipna=False))
        assert_eq(
            dds.prod(skipna=False, split_every=split_every), pds.prod(skipna=False)
        )
        assert_eq(
            dds.product(skipna=False, split_every=split_every),
            pds.product(skipna=False),
        )
        assert_eq(dds.min(skipna=False, split_every=split_every), pds.min(skipna=False))
        assert_eq(dds.max(skipna=False, split_every=split_every), pds.max(skipna=False))
        assert_eq(dds.std(skipna=False, split_every=split_every), pds.std(skipna=False))
        assert_eq(dds.var(skipna=False, split_every=split_every), pds.var(skipna=False))
        assert_eq(dds.sem(skipna=False, split_every=split_every), pds.sem(skipna=False))
        assert_eq(
            dds.std(skipna=False, ddof=0, split_every=split_every),
            pds.std(skipna=False, ddof=0),
        )
        assert_eq(
            dds.var(skipna=False, ddof=0, split_every=split_every),
            pds.var(skipna=False, ddof=0),
        )
        assert_eq(
            dds.sem(skipna=False, ddof=0, split_every=split_every),
            pds.sem(skipna=False, ddof=0),
        )
        assert_eq(
            dds.mean(skipna=False, split_every=split_every), pds.mean(skipna=False)
        )

    # testing index
    assert_eq(ddf1.index.min(split_every=split_every), pdf1.index.min())
    assert_eq(ddf1.index.max(split_every=split_every), pdf1.index.max())
    assert_eq(ddf1.index.count(split_every=split_every), pd.notnull(pdf1.index).sum())


@pytest.mark.parametrize("split_every", [False, 2])
def test_reductions_timedelta(split_every):
    ds = pd.Series(pd.to_timedelta([2, 3, 4, np.nan, 5]))
    dds = dd.from_pandas(ds, 2)

    assert_eq(dds.sum(split_every=split_every), ds.sum())
    assert_eq(dds.min(split_every=split_every), ds.min())
    assert_eq(dds.max(split_every=split_every), ds.max())
    assert_eq(dds.count(split_every=split_every), ds.count())


@pytest.mark.parametrize("axis", [0, 1])
@pytest.mark.parametrize(
    "redfunc",
    ["sum", "prod", "product", "min", "max", "mean", "var", "std", "all", "any"],
)
def test_reductions_numpy_dispatch(axis, redfunc):
    pdf = pd.DataFrame({"a": [1, 2, 3], "b": [4, 5, 6]}, index=[0, 1, 3])
    df = dd.from_pandas(pdf, 3)

    if NUMPY_GE_200 and redfunc == "product":
        with pytest.raises(
            AttributeError, match="module 'numpy' has no attribute 'product'"
        ):
            np_redfunc = getattr(np, redfunc)
        pytest.skip(reason="numpy>2 removed product -- users should move to prod")

    np_redfunc = getattr(np, redfunc)
    if redfunc in ("var", "std"):
        # numpy has default ddof value 0 while
        # dask and pandas have 1, so ddof should be passed
        # explicitly when calling np.var(dask)
        expect = np_redfunc(pdf, axis=axis, ddof=1)
        actual = np_redfunc(df, axis=axis, ddof=1)
    elif NUMPY_GE_125 and redfunc == "product" and not NUMPY_GE_200:
        expect = np_redfunc(pdf, axis=axis)
        with pytest.warns(DeprecationWarning, match="`product` is deprecated"):
            actual = np_redfunc(df, axis=axis)
    else:
        expect = np_redfunc(pdf, axis=axis)
        actual = np_redfunc(df, axis=axis)

    assert_eq(expect, actual)


@pytest.mark.parametrize("split_every", [False, 2])
@pytest.mark.xfail_with_pyarrow_strings
def test_allany(split_every):
    df = pd.DataFrame(
        np.random.choice([True, False], size=(100, 4)), columns=["A", "B", "C", "D"]
    )
    df["E"] = list("abcde") * 20
    ddf = dd.from_pandas(df, 10)

    assert_eq(ddf.all(split_every=split_every), df.all())
    assert_eq(ddf.all(axis=1, split_every=split_every), df.all(axis=1))
    assert_eq(ddf.all(axis=0, split_every=split_every), df.all(axis=0))

    assert_eq(ddf.any(split_every=split_every), df.any())
    assert_eq(ddf.any(axis=1, split_every=split_every), df.any(axis=1))
    assert_eq(ddf.any(axis=0, split_every=split_every), df.any(axis=0))

    assert_eq(ddf.A.all(split_every=split_every), df.A.all())
    assert_eq(ddf.A.any(split_every=split_every), df.A.any())


@pytest.mark.parametrize("split_every", [False, 2])
def test_deterministic_reduction_names(split_every):
    df = pd.DataFrame({"x": [1, 2, 3, 4], "y": [5, 6, 7, 8]})
    ddf = dd.from_pandas(df, npartitions=2)

    for x in [ddf, ddf.x]:
        assert (
            x.sum(split_every=split_every)._name == x.sum(split_every=split_every)._name
        )
        assert (
            x.prod(split_every=split_every)._name
            == x.prod(split_every=split_every)._name
        )
        assert (
            x.product(split_every=split_every)._name
            == x.product(split_every=split_every)._name
        )
        assert (
            x.min(split_every=split_every)._name == x.min(split_every=split_every)._name
        )
        assert (
            x.max(split_every=split_every)._name == x.max(split_every=split_every)._name
        )
        assert (
            x.count(split_every=split_every)._name
            == x.count(split_every=split_every)._name
        )
        assert (
            x.std(split_every=split_every)._name == x.std(split_every=split_every)._name
        )
        assert (
            x.var(split_every=split_every)._name == x.var(split_every=split_every)._name
        )
        assert (
            x.sem(split_every=split_every)._name == x.sem(split_every=split_every)._name
        )
        assert (
            x.mean(split_every=split_every)._name
            == x.mean(split_every=split_every)._name
        )

    assert (
        ddf.x.nunique(split_every=split_every)._name
        == ddf.x.nunique(split_every=split_every)._name
    )


def test_reduction_series_invalid_axis():
    dsk = {
        ("x", 0): pd.DataFrame({"a": [1, 2, 3], "b": [4, 5, 6]}, index=[0, 1, 3]),
        ("x", 1): pd.DataFrame({"a": [4, 5, 6], "b": [3, 2, 1]}, index=[5, 6, 8]),
        ("x", 2): pd.DataFrame({"a": [7, 8, 9], "b": [0, 0, 0]}, index=[9, 9, 9]),
    }
    ddf1 = dd.repartition(pd.concat(dsk.values()), [0, 4, 9, 9])
    pdf1 = ddf1.compute()

    for axis in [1, "columns"]:
        for s in [ddf1.a, pdf1.a]:  # both must behave the same
            pytest.raises(ValueError, lambda s=s, axis=axis: s.sum(axis=axis))
            pytest.raises(ValueError, lambda s=s, axis=axis: s.prod(axis=axis))
            pytest.raises(ValueError, lambda s=s, axis=axis: s.product(axis=axis))
            pytest.raises(ValueError, lambda s=s, axis=axis: s.min(axis=axis))
            pytest.raises(ValueError, lambda s=s, axis=axis: s.max(axis=axis))
            # only count doesn't have axis keyword
            pytest.raises(
                (TypeError, ValueError), lambda s=s, axis=axis: s.count(axis=axis)
            )
            pytest.raises(ValueError, lambda s=s, axis=axis: s.std(axis=axis))
            pytest.raises(ValueError, lambda s=s, axis=axis: s.var(axis=axis))
            pytest.raises(ValueError, lambda s=s, axis=axis: s.sem(axis=axis))
            pytest.raises(ValueError, lambda s=s, axis=axis: s.mean(axis=axis))


def test_reductions_non_numeric_dtypes():
    # test non-numric blocks
    if pyarrow_strings_enabled() and not PANDAS_GE_300:
        pytest.xfail(reason="known failure with arrow strings")

    def check_raises(d, p, func):
        pytest.raises((TypeError, ValueError), lambda: getattr(d, func)().compute())
        pytest.raises((TypeError, ValueError), lambda: getattr(p, func)())

    pds = pd.Series(["a", "b", "c", "d", "e"])
    dds = dd.from_pandas(pds, 2)
    assert_eq(dds.sum(), pds.sum())
    check_raises(dds, pds, "prod")
    check_raises(dds, pds, "product")
    assert_eq(dds.min(), pds.min())
    assert_eq(dds.max(), pds.max())
    assert_eq(dds.count(), pds.count())
    check_raises(dds, pds, "std")
    check_raises(dds, pds, "var")
    check_raises(dds, pds, "sem")
    check_raises(dds, pds, "skew")
    check_raises(dds, pds, "kurtosis")
    assert_eq(dds.nunique(), pds.nunique())

    for pds in [
        pd.Series(pd.Categorical([1, 2, 3, 4, 5], ordered=True)),
        pd.Series(pd.Categorical(list("abcde"), ordered=True)),
        pd.Series(pd.date_range("2011-01-01", freq="D", periods=5)),
    ]:
        dds = dd.from_pandas(pds, 2)

        check_raises(dds, pds, "sum")
        check_raises(dds, pds, "prod")
        check_raises(dds, pds, "product")
        assert_eq(dds.min(), pds.min())
        assert_eq(dds.max(), pds.max())
        assert_eq(dds.count(), pds.count())
        if pds.dtype != "datetime64[ns]":
            # std is implemented for datetimes in pandas 1.2.0, but dask
            # implementation depends on var which isn't
            check_raises(dds, pds, "std")
        check_raises(dds, pds, "var")
        check_raises(dds, pds, "sem")
        check_raises(dds, pds, "skew")
        check_raises(dds, pds, "kurtosis")
        assert_eq(dds.nunique(), pds.nunique())

    pds = pd.Series(pd.timedelta_range("1 days", freq="D", periods=5))
    dds = dd.from_pandas(pds, 2)
    assert_eq(dds.sum(), pds.sum())
    assert_eq(dds.min(), pds.min())
    assert_eq(dds.max(), pds.max())
    assert_eq(dds.count(), pds.count())
    # both pandas and dask skew calculations do not support timedelta
    check_raises(dds, pds, "skew")
    check_raises(dds, pds, "kurtosis")

    # ToDo: pandas supports timedelta std, dask returns float64
    # assert_eq(dds.std(), pds.std())

    # ToDo: pandas supports timedelta std, otherwise dask raises:
    # TypeError: unsupported operand type(s) for *: 'float' and 'Timedelta'
    # assert_eq(dds.mean(), pds.mean())

    assert_eq(dds.nunique(), pds.nunique())


@pytest.mark.parametrize("split_every", [False, 2])
def test_reductions_frame(split_every):
    dsk = {
        ("x", 0): pd.DataFrame({"a": [1, 2, 3], "b": [4, 5, 6]}, index=[0, 1, 3]),
        ("x", 1): pd.DataFrame({"a": [4, 5, 6], "b": [3, 2, 1]}, index=[5, 6, 8]),
        ("x", 2): pd.DataFrame({"a": [7, 8, 9], "b": [0, 0, 0]}, index=[9, 9, 9]),
    }
    ddf1 = dd.repartition(pd.concat(dsk.values()), [0, 4, 9, 9])
    pdf1 = ddf1.compute()

    assert_eq(ddf1.sum(split_every=split_every), pdf1.sum())
    assert_eq(ddf1.prod(split_every=split_every), pdf1.prod())
    assert_eq(ddf1.product(split_every=split_every), pdf1.product())
    assert_eq(ddf1.min(split_every=split_every), pdf1.min())
    assert_eq(ddf1.max(split_every=split_every), pdf1.max())
    assert_eq(ddf1.count(split_every=split_every), pdf1.count())
    assert_eq(ddf1.std(split_every=split_every), pdf1.std())
    assert_eq(ddf1.var(split_every=split_every), pdf1.var())
    assert_eq(ddf1.sem(split_every=split_every), pdf1.sem())
    assert_eq(ddf1.std(ddof=0, split_every=split_every), pdf1.std(ddof=0))
    assert_eq(ddf1.var(ddof=0, split_every=split_every), pdf1.var(ddof=0))
    assert_eq(ddf1.sem(ddof=0, split_every=split_every), pdf1.sem(ddof=0))
    assert_eq(ddf1.mean(split_every=split_every), pdf1.mean())

    for axis in [0, 1, "index", "columns"]:
        assert_eq(ddf1.sum(axis=axis, split_every=split_every), pdf1.sum(axis=axis))
        assert_eq(ddf1.prod(axis=axis, split_every=split_every), pdf1.prod(axis=axis))
        assert_eq(
            ddf1.product(axis=axis, split_every=split_every), pdf1.product(axis=axis)
        )
        assert_eq(ddf1.min(axis=axis, split_every=split_every), pdf1.min(axis=axis))
        assert_eq(ddf1.max(axis=axis, split_every=split_every), pdf1.max(axis=axis))
        assert_eq(ddf1.count(axis=axis, split_every=split_every), pdf1.count(axis=axis))
        assert_eq(ddf1.std(axis=axis, split_every=split_every), pdf1.std(axis=axis))
        assert_eq(ddf1.var(axis=axis, split_every=split_every), pdf1.var(axis=axis))
        assert_eq(ddf1.sem(axis=axis, split_every=split_every), pdf1.sem(axis=axis))
        assert_eq(
            ddf1.std(axis=axis, ddof=0, split_every=split_every),
            pdf1.std(axis=axis, ddof=0),
        )
        assert_eq(
            ddf1.var(axis=axis, ddof=0, split_every=split_every),
            pdf1.var(axis=axis, ddof=0),
        )
        assert_eq(
            ddf1.sem(axis=axis, ddof=0, split_every=split_every),
            pdf1.sem(axis=axis, ddof=0),
        )
        assert_eq(ddf1.mean(axis=axis, split_every=split_every), pdf1.mean(axis=axis))

    pytest.raises(ValueError, lambda: ddf1.sum(axis="incorrect").compute())

    # min
    result = ddf1.min(axis=None, split_every=split_every)
    expected = pdf1.min(axis=None)
    assert_eq(result, expected)
    # max
    result = ddf1.max(axis=None, split_every=split_every)
    expected = pdf1.max(axis=None)
    assert_eq(result, expected)
    # mean
    result = ddf1.mean(axis=None, split_every=split_every)
    expected = pdf1.mean(axis=None)
    assert_eq(result, expected, check_dtype=False)


@pytest.mark.parametrize(
    "func, kwargs",
    [
        ("sum", None),
        ("prod", None),
        ("product", None),
        ("mean", None),
        ("std", None),
        ("std", {"ddof": 0}),
        ("std", {"skipna": False}),
        ("std", {"ddof": 0, "skipna": False}),
        ("min", None),
        ("max", None),
        ("count", None),
        ("sem", None),
        ("sem", {"ddof": 0}),
        ("sem", {"skipna": False}),
        ("sem", {"ddof": 0, "skipna": False}),
        ("var", None),
        ("var", {"ddof": 0}),
        ("var", {"skipna": False}),
        ("var", {"ddof": 0, "skipna": False}),
    ],
)
@pytest.mark.parametrize(
    "numeric_only",
    [
        None,
        True,
        pytest.param(
            False,
            marks=pytest.mark.xfail(
                True, reason="numeric_only=False not implemented", strict=False
            ),
        ),
    ],
)
def test_reductions_frame_dtypes(func, kwargs, numeric_only):
    if pyarrow_strings_enabled() and func == "sum" and numeric_only is None:
        pytest.xfail("Known failure with pyarrow strings")
    df = pd.DataFrame(
        {
            "int": [1, 2, 3, 4, 5, 6, 7, 8],
            "float": [1.0, 2.0, 3.0, 4.0, np.nan, 6.0, 7.0, 8.0],
            "dt": [pd.NaT] + [datetime(2011, i, 1) for i in range(1, 8)],
            "str": list("abcdefgh"),
            "timedelta": pd.to_timedelta([1, 2, 3, 4, 5, 6, 7, np.nan]),
            "bool": [True, False] * 4,
        }
    )

    if kwargs is None:
        kwargs = {}

    if numeric_only is False or numeric_only is None:
        if func in ("sum", "prod", "product", "mean", "median", "std", "sem", "var"):
            # datetime columns don't support some aggs
            df = df.drop(columns=["dt", "timedelta"])
        if func in ("prod", "product", "mean", "std", "sem", "var"):
            # string columns don't support some other aggs
            df = df.drop(columns=["str"])

    if numeric_only is not None:
        kwargs["numeric_only"] = numeric_only

    ddf = dd.from_pandas(df, 3)

    expected = getattr(df, func)(**kwargs)
    actual = getattr(ddf, func)(**kwargs)
    assert_eq(expected, actual)


def test_count_numeric_only_axis_one():
    df = pd.DataFrame(
        {
            "int": [1, 2, 3, 4, 5, 6, 7, 8],
            "float": [1.0, 2.0, 3.0, 4.0, np.nan, 6.0, 7.0, 8.0],
            "dt": [pd.NaT] + [datetime(2011, i, 1) for i in range(1, 8)],
            "str": list("abcdefgh"),
            "timedelta": pd.to_timedelta([1, 2, 3, 4, 5, 6, 7, np.nan]),
            "bool": [True, False] * 4,
        }
    )
    ddf = dd.from_pandas(df, npartitions=2)

    assert_eq(ddf.count(axis=1), df.count(axis=1))
    assert_eq(
        ddf.count(numeric_only=False, axis=1), df.count(numeric_only=False, axis=1)
    )
    assert_eq(ddf.count(numeric_only=True, axis=1), df.count(numeric_only=True, axis=1))


@pytest.mark.parametrize(
    "func", ["sum", "prod", "product", "min", "max", "count", "std", "var", "quantile"]
)
def test_reductions_frame_dtypes_numeric_only_supported(func):
    df = pd.DataFrame(
        {
            "int": [1, 2, 3, 4, 5, 6, 7, 8],
            "float": [1.0, 2.0, 3.0, 4.0, np.nan, 6.0, 7.0, 8.0],
            "dt": [pd.NaT] + [datetime(2011, i, 1) for i in range(1, 8)],
            "str": list("abcdefgh"),
            "timedelta": pd.to_timedelta([1, 2, 3, 4, 5, 6, 7, np.nan]),
            "bool": [True, False] * 4,
        }
    )
    npartitions = 3
    if func == "quantile":
        # bool doesn't work in pandas quantile
        df = df.drop(columns="bool")
        npartitions = 1  # https://github.com/dask/dask/issues/9227

    ddf = dd.from_pandas(df, npartitions)

    numeric_only_false_raises = ["sum", "prod", "product", "std", "var", "quantile"]

    # `numeric_only=True` is always supported
    assert_eq(
        getattr(df, func)(numeric_only=True),
        getattr(ddf, func)(numeric_only=True),
    )
    errors = (
        (ValueError, TypeError)
        if pa is None
        else (ValueError, TypeError, ArrowNotImplementedError)
    )

    # `numeric_only=False`
    if func in numeric_only_false_raises:
        with pytest.raises(
            errors,
            match="'DatetimeArray' with dtype datetime64.*|"
            "'DatetimeArray' does not implement reduction|could not convert|"
            "'ArrowStringArray' with dtype string"
            "|unsupported operand|no kernel|not supported"
            "|Cannot perform reduction",
        ):
            getattr(ddf, func)(numeric_only=False)

    else:
        assert_eq(
            getattr(df, func)(numeric_only=False),
            getattr(ddf, func)(numeric_only=False),
        )

    # `numeric_only` default value
    if func in numeric_only_false_raises:
        with pytest.raises(
            errors,
            match="'DatetimeArray' with dtype datetime64.*|"
            "'DatetimeArray' does not implement reduction|could not convert|"
            "'ArrowStringArray' with dtype string"
            "|unsupported operand|no kernel|not supported"
            "|Cannot perform reduction",
        ):
            getattr(ddf, func)()
    else:
        assert_eq(
            getattr(df, func)(),
            getattr(ddf, func)(),
        )

    num_cols = ["int", "float"]
    if func != "quantile":
        num_cols.append("bool")

    df_numerics = df[num_cols]
    ddf_numerics = ddf[num_cols]

    assert_eq(
        getattr(df_numerics, func)(),
        getattr(ddf_numerics, func)(),
    )
    assert_eq(
        getattr(df_numerics, func)(numeric_only=False),
        getattr(ddf_numerics, func)(numeric_only=False),
    )


@pytest.mark.parametrize(
    "func",
    [
        "mean",
        "sem",
    ],
)
def test_reductions_frame_dtypes_numeric_only(func):
    df = pd.DataFrame(
        {
            "int": [1, 2, 3, 4, 5, 6, 7, 8],
            "float": [1.0, 2.0, 3.0, 4.0, np.nan, 6.0, 7.0, 8.0],
            "dt": [pd.NaT] + [datetime(2011, i, 1) for i in range(1, 8)],
            "str": list("abcdefgh"),
            "timedelta": pd.to_timedelta([1, 2, 3, 4, 5, 6, 7, np.nan]),
            "bool": [True, False] * 4,
        }
    )

    ddf = dd.from_pandas(df, 3)
    kwargs = {"numeric_only": True}

    assert_eq(
        getattr(df, func)(**kwargs),
        getattr(ddf, func)(**kwargs),
    )
    assert_eq(df.sem(ddof=0, **kwargs), ddf.sem(ddof=0, **kwargs))
    assert_eq(df.std(ddof=0, **kwargs), ddf.std(ddof=0, **kwargs))
    assert_eq(df.var(ddof=0, **kwargs), ddf.var(ddof=0, **kwargs))
    assert_eq(df.var(skipna=False, **kwargs), ddf.var(skipna=False, **kwargs))
    assert_eq(
        df.var(skipna=False, ddof=0, **kwargs), ddf.var(skipna=False, ddof=0, **kwargs)
    )

    df_numerics = df[["int", "float", "bool"]]
    ddf_numerics = ddf[["int", "float", "bool"]]

    assert_eq(
        getattr(df_numerics, func)(),
        getattr(ddf_numerics, func)(),
    )


@pytest.mark.parametrize("func", ["skew", "kurtosis"])
def test_skew_kurt_numeric_only_false(func):
    pytest.importorskip("scipy.stats")
    df = pd.DataFrame(
        {
            "int": [1, 2, 3, 4, 5, 6, 7, 8],
            "float": [1.0, 2.0, 3.0, 4.0, np.nan, 6.0, 7.0, 8.0],
            "dt": [pd.NaT] + [datetime(2010, i, 1) for i in range(1, 8)],
        }
    )
    ddf = dd.from_pandas(df, npartitions=2)

    ctx = pytest.raises(TypeError, match="does not support|does not implement")

    with ctx:
        getattr(df, func)(numeric_only=False)
    with ctx:
        getattr(ddf, func)(numeric_only=False)

    with ctx:
        getattr(df, func)()
    with ctx:
        getattr(ddf, func)()


@pytest.mark.parametrize("split_every", [False, 2])
def test_reductions_frame_nan(split_every):
    df = pd.DataFrame(
        {
            "a": [1, 2, np.nan, 4, 5, 6, 7, 8],
            "b": [1, 2, np.nan, np.nan, np.nan, 5, np.nan, np.nan],
            "c": [np.nan] * 8,
        }
    )
    ddf = dd.from_pandas(df, 3)
    assert_eq(df.sum(), ddf.sum(split_every=split_every))
    assert_eq(df.prod(), ddf.prod(split_every=split_every))
    assert_eq(df.product(), ddf.product(split_every=split_every))
    assert_eq(df.min(), ddf.min(split_every=split_every))
    assert_eq(df.max(), ddf.max(split_every=split_every))
    assert_eq(df.count(), ddf.count(split_every=split_every))
    with warnings.catch_warnings():
        # dask.dataframe should probably filter this, to match pandas, but
        # it seems quite difficult.
        warnings.simplefilter("ignore", RuntimeWarning)
        assert_eq(df.std(), ddf.std(split_every=split_every))
        assert_eq(df.var(), ddf.var(split_every=split_every))
        assert_eq(df.sem(), ddf.sem(split_every=split_every))
        assert_eq(df.std(ddof=0), ddf.std(ddof=0, split_every=split_every))
        assert_eq(df.var(ddof=0), ddf.var(ddof=0, split_every=split_every))
        assert_eq(df.sem(ddof=0), ddf.sem(ddof=0, split_every=split_every))
    assert_eq(df.mean(), ddf.mean(split_every=split_every))

    with warnings.catch_warnings(record=True):
        assert_eq(df.sum(skipna=False), ddf.sum(skipna=False, split_every=split_every))
        assert_eq(
            df.prod(skipna=False), ddf.prod(skipna=False, split_every=split_every)
        )
        assert_eq(
            df.product(skipna=False), ddf.product(skipna=False, split_every=split_every)
        )
        assert_eq(df.min(skipna=False), ddf.min(skipna=False, split_every=split_every))
        assert_eq(df.max(skipna=False), ddf.max(skipna=False, split_every=split_every))
        assert_eq(df.std(skipna=False), ddf.std(skipna=False, split_every=split_every))
        assert_eq(df.var(skipna=False), ddf.var(skipna=False, split_every=split_every))
        assert_eq(df.sem(skipna=False), ddf.sem(skipna=False, split_every=split_every))
        assert_eq(
            df.std(skipna=False, ddof=0),
            ddf.std(skipna=False, ddof=0, split_every=split_every),
        )
        assert_eq(
            df.var(skipna=False, ddof=0),
            ddf.var(skipna=False, ddof=0, split_every=split_every),
        )
        assert_eq(
            df.sem(skipna=False, ddof=0),
            ddf.sem(skipna=False, ddof=0, split_every=split_every),
        )
        assert_eq(
            df.mean(skipna=False), ddf.mean(skipna=False, split_every=split_every)
        )

        assert_eq(
            df.sum(axis=1, skipna=False),
            ddf.sum(axis=1, skipna=False, split_every=split_every),
        )
        assert_eq(
            df.prod(axis=1, skipna=False),
            ddf.prod(axis=1, skipna=False, split_every=split_every),
        )
        assert_eq(
            df.product(axis=1, skipna=False),
            ddf.product(axis=1, skipna=False, split_every=split_every),
        )
        assert_eq(
            df.min(axis=1, skipna=False),
            ddf.min(axis=1, skipna=False, split_every=split_every),
        )
        assert_eq(
            df.max(axis=1, skipna=False),
            ddf.max(axis=1, skipna=False, split_every=split_every),
        )
        assert_eq(
            df.std(axis=1, skipna=False),
            ddf.std(axis=1, skipna=False, split_every=split_every),
        )
        assert_eq(
            df.var(axis=1, skipna=False),
            ddf.var(axis=1, skipna=False, split_every=split_every),
        )
        assert_eq(
            df.sem(axis=1, skipna=False),
            ddf.sem(axis=1, skipna=False, split_every=split_every),
        )
        assert_eq(
            df.std(axis=1, skipna=False, ddof=0),
            ddf.std(axis=1, skipna=False, ddof=0, split_every=split_every),
        )
        assert_eq(
            df.var(axis=1, skipna=False, ddof=0),
            ddf.var(axis=1, skipna=False, ddof=0, split_every=split_every),
        )
        assert_eq(
            df.sem(axis=1, skipna=False, ddof=0),
            ddf.sem(axis=1, skipna=False, ddof=0, split_every=split_every),
        )
        assert_eq(
            df.mean(axis=1, skipna=False),
            ddf.mean(axis=1, skipna=False, split_every=split_every),
        )


@pytest.mark.parametrize("comparison", ["lt", "gt", "le", "ge", "ne", "eq"])
def test_series_comparison_nan(comparison):
    s = pd.Series([1, 2, 3, 4, 5, 6, 7])
    s_nan = pd.Series([1, -1, 8, np.nan, 5, 6, 2.4])
    ds = dd.from_pandas(s, 3)
    ds_nan = dd.from_pandas(s_nan, 3)

    fill_value = 7
    comparison_pd = getattr(s, comparison)
    comparison_dd = getattr(ds, comparison)
    assert_eq(
        comparison_dd(ds_nan, fill_value=fill_value),
        comparison_pd(s_nan, fill_value=fill_value),
    )


def test_sum_intna():
    a = pd.Series([1, None, 2], dtype=pd.Int32Dtype())
    b = dd.from_pandas(a, 2)
    assert_eq(a.sum(), b.sum())


def test_divmod():
    df1 = pd.Series(np.random.rand(10))
    df2 = pd.Series(np.random.rand(10))

    ddf1 = dd.from_pandas(df1, npartitions=3)
    ddf2 = dd.from_pandas(df2, npartitions=3)

    result = divmod(ddf1, 2.0)
    expected = divmod(df1, 2.0)
    assert_eq(result[0], expected[0])
    assert_eq(result[1], expected[1])

    result = divmod(ddf1, ddf2)
    expected = divmod(df1, df2)
    assert_eq(result[0], expected[0])
    assert_eq(result[1], expected[1])


@pytest.mark.skipif("not scipy")
def test_moment():
    from dask.array import stats
    from dask.array.utils import assert_eq

    df = pd.Series(list(range(10)))
    ddf = dd.from_pandas(df, npartitions=2)

    assert_eq(stats.moment(ddf, 2, 0), scipy.stats.moment(df, 2, 0))


@pytest.mark.parametrize("func", ["sum", "count", "mean", "var", "sem"])
def test_empty_df_reductions(func):
    pdf = pd.DataFrame()
    ddf = dd.from_pandas(pdf, npartitions=1)

    dsk_func = getattr(ddf.__class__, func)
    pd_func = getattr(pdf.__class__, func)

    assert_eq(dsk_func(ddf), pd_func(pdf))

    idx = pd.date_range("2000", periods=4)
    pdf = pd.DataFrame(index=idx)
    ddf = dd.from_pandas(pdf, npartitions=1)

    assert_eq(dsk_func(ddf), pd_func(pdf))


@pytest.mark.parametrize("method", ["sum", "prod", "product"])
@pytest.mark.parametrize("min_count", [0, 9])
def test_series_agg_with_min_count(method, min_count):
    df = pd.DataFrame([[1]], columns=["a"])
    ddf = dd.from_pandas(df, npartitions=1)
    func = getattr(ddf["a"], method)
    result = func(min_count=min_count).compute()
    if min_count == 0:
        assert result == 1
    else:
        assert result is np.nan or pd.isna(result)


# Default absolute tolerance of 2000 nanoseconds
def assert_near_timedeltas(t1, t2, atol=2000):
    if is_scalar(t1):
        t1 = pd.Series([t1])
    if is_scalar(t2):
        t2 = pd.Series([t2])

    assert t1.dtype == t2.dtype
    assert_eq(pd.to_numeric(t1), pd.to_numeric(t2), atol=atol)


@pytest.mark.parametrize("axis", [0, 1])
@pytest.mark.parametrize("numeric_only", [True, False, None])
def test_datetime_std_creates_copy_cols(axis, numeric_only):
    pdf = pd.DataFrame(
        {
            "dt1": [
                datetime.fromtimestamp(1636426700 + (i * 250000)) for i in range(10)
            ],
            "dt2": [
                datetime.fromtimestamp(1636426700 + (i * 300000)) for i in range(10)
            ],
        }
    )

    ddf = dd.from_pandas(pdf, 3)

    kwargs = {} if numeric_only is None else {"numeric_only": numeric_only}

    # Series test (same line twice to make sure data structure wasn't mutated)
    assert_eq(ddf["dt1"].std(**kwargs), pdf["dt1"].std(**kwargs))
    assert_eq(ddf["dt1"].std(**kwargs), pdf["dt1"].std(**kwargs))

    # DataFrame test (same line twice to make sure data structure wasn't mutated)
    expected = pdf.std(axis=axis, **kwargs)
    result = ddf.std(axis=axis, **kwargs)
    assert_near_timedeltas(result.compute(), expected)

    expected = pdf.std(axis=axis, **kwargs)
    result = ddf.std(axis=axis, **kwargs)
    assert_near_timedeltas(result.compute(), expected)


@pytest.mark.filterwarnings("ignore::RuntimeWarning")
@pytest.mark.parametrize("axis", [0, 1])
@pytest.mark.parametrize("skipna", [False, True])
@pytest.mark.parametrize("numeric_only", [True, False, None])
def test_datetime_std_with_larger_dataset(axis, skipna, numeric_only):
    num_rows = 250

    dt1 = pd.concat(
        [
            pd.Series([pd.NaT] * 15, index=range(15)),
            pd.to_datetime(
                pd.Series(
                    [
                        datetime.fromtimestamp(1636426704 + (i * 250000))
                        for i in range(num_rows - 15)
                    ],
                    index=range(15, 250),
                )
            ),
        ],
        ignore_index=False,
    )

    base_numbers = [
        (1638290040706793300 + (i * 69527182702409)) for i in range(num_rows)
    ]

    pdf = pd.DataFrame(
        {"dt1": dt1, "dt2": pd.to_datetime(pd.Series(base_numbers))}, index=range(250)
    )

    for i in range(3, 8):
        pdf[f"dt{i}"] = pd.to_datetime(
            pd.Series([int(x + (0.12 * i)) for x in base_numbers])
        )

    ddf = dd.from_pandas(pdf, 8)

    kwargs = {} if numeric_only is None else {"numeric_only": numeric_only}
    kwargs["skipna"] = skipna

    # Same thing but as Series. No axis, since axis=1 raises error
    assert_near_timedeltas(ddf["dt1"].std(**kwargs).compute(), pdf["dt1"].std(**kwargs))

    # Computation on full dataset
    expected = pdf.std(axis=axis, **kwargs)
    result = ddf.std(axis=axis, **kwargs).compute()
    assert_near_timedeltas(result, expected)


@pytest.mark.parametrize("skipna", [False, True])
@pytest.mark.parametrize("numeric_only", [True, False, None])
def test_datetime_std_across_axis1_null_results(skipna, numeric_only):
    pdf = pd.DataFrame(
        {
            "dt1": [
                datetime.fromtimestamp(1636426704 + (i * 250000)) for i in range(10)
            ],
            "dt2": [
                datetime.fromtimestamp(1636426704 + (i * 217790)) for i in range(10)
            ],
            "nums": [i for i in range(10)],
        }
    )

    ddf = dd.from_pandas(pdf, 3)

    kwargs = {} if numeric_only is None else {"numeric_only": numeric_only}
    kwargs["skipna"] = skipna

    ctx = contextlib.nullcontext()
    success = True
    if numeric_only is False or numeric_only is None:
        ctx = pytest.raises(TypeError)
        success = False
    elif numeric_only is None:
        ctx = pytest.warns(FutureWarning, match="numeric_only")

    # Single column always results in NaT
    expected = pdf[["dt1"]].std(axis=1, **kwargs)
    result = ddf[["dt1"]].std(axis=1, **kwargs)
    if success:
        assert_eq(result, expected)

    # Mix of datetimes with other numeric types produces NaNs
    with ctx:
        expected = pdf.std(axis=1, **kwargs)
    with ctx:
        result = ddf.std(axis=1, **kwargs)
    if success:
        assert_eq(result, expected)

    # Test with mix of na and truthy datetimes
    pdf2 = pd.DataFrame(
        {
            "dt1": [pd.NaT]
            + [datetime.fromtimestamp(1636426704 + (i * 250000)) for i in range(10)]
            + [pd.NaT],
            "dt2": [
                datetime.fromtimestamp(1636426704 + (i * 250000)) for i in range(12)
            ],
            "dt3": [
                datetime.fromtimestamp(1636426704 + (i * 282616)) for i in range(12)
            ],
        }
    )

    ddf2 = dd.from_pandas(pdf2, 3)

    expected = pdf2.std(axis=1, **kwargs)
    result = ddf2.std(axis=1, **kwargs)
    if success:
        assert_eq(result, expected)


def test_std_raises_on_index():
    with pytest.raises(
        (NotImplementedError, AttributeError),
        match="`std` is only supported with objects that are Dataframes or Series|has no attribute",
    ):
        dd.from_pandas(pd.DataFrame({"test": [1, 2]}), npartitions=2).index.std()


def test_std_raises_with_arrow_string_ea():
    pa = pytest.importorskip("pyarrow")
    ser = pd.Series(["a", "b", "c"], dtype=pd.ArrowDtype(pa.string()))
    ds = dd.from_pandas(ser, npartitions=2)
    with pytest.raises(ValueError, match="`std` not supported with string series"):
        ds.std()


@pytest.mark.parametrize(
    "dtype",
    [
        pytest.param(
            "int64[pyarrow]",
            marks=pytest.mark.skipif(pa is None, reason="requires pyarrow installed"),
        ),
        pytest.param(
            "float64[pyarrow]",
            marks=pytest.mark.skipif(pa is None, reason="requires pyarrow installed"),
        ),
        "Int64",
        "Int32",
        "Float64",
        "UInt64",
    ],
)
@pytest.mark.parametrize("func", ["std", "var", "skew", "kurtosis"])
def test_reductions_with_pandas_and_arrow_ea(dtype, func):
    if func in ["skew", "kurtosis"]:
        pytest.importorskip("scipy")
        if "pyarrow" in dtype:
            pytest.xfail("skew/kurtosis not implemented for arrow dtypes")

    ser = pd.Series([1, 2, 3, 4], dtype=dtype)
    ds = dd.from_pandas(ser, npartitions=2)
    pd_result = getattr(ser, func)()
    dd_result = getattr(ds, func)()
    if func == "kurtosis":
        n = ser.shape[0]
        factor = ((n - 1) * (n + 1)) / ((n - 2) * (n - 3))
        offset = (6 * (n - 1)) / ((n - 2) * (n - 3))
        dd_result = factor * dd_result + offset
    # _meta is wrongly NA
    assert_eq(dd_result, pd_result, check_dtype=False)


def test_quantile_arrow_dtype():
    pdf = pd.DataFrame({"foo": [1, 1, 1, 1, 1, 1]}).convert_dtypes(
        dtype_backend="pyarrow"
    )
    df = dd.from_pandas(pdf, npartitions=2)
    assert_eq(df.quantile(0.25), pdf.quantile(0.25), check_dtype=False)
