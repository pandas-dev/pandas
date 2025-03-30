from __future__ import annotations

import contextlib
import operator
import warnings

import numpy as np
import pandas as pd
import pytest

import dask
import dask.dataframe as dd
from dask.dataframe import _compat
from dask.dataframe._compat import PANDAS_GE_210, PANDAS_GE_300, tm
from dask.dataframe._pyarrow import to_pyarrow_string
from dask.dataframe.core import _concat
from dask.dataframe.utils import assert_eq, get_string_dtype, pyarrow_strings_enabled

# Generate a list of categorical series and indices
cat_series = []
for ordered in [True, False]:
    s = pd.Series(pd.Categorical(list("bacbac"), ordered=ordered))
    ds = dd.from_pandas(s, npartitions=2)
    cat_series.append((s, ds))
s = pd.Series(range(6), index=pd.Categorical(list("bacbac")))
ds = dd.from_pandas(s, npartitions=2)
cat_series.append((ds.compute().index, ds.index))


a = pd.DataFrame(
    {
        "v": list("abcde"),
        "w": list("xxxxx"),
        "x": np.arange(5),
        "y": list("abcbc"),
        "z": np.arange(5, dtype="f8"),
    }
)

b = pd.DataFrame(
    {
        "v": list("fghij"),
        "w": list("yyyyy"),
        "x": np.arange(5, 10),
        "y": list("abbba"),
        "z": np.arange(5, 10, dtype="f8"),
    }
)

c = pd.DataFrame(
    {
        "v": list("klmno"),
        "w": list("zzzzz"),
        "x": np.arange(10, 15),
        "y": list("bcbcc"),
        "z": np.arange(10, 15, dtype="f8"),
    }
)

frames = [a, b, c]
frames2 = []
for df in frames:
    df.w = df.w.astype("category")
    df.y = df.y.astype("category")
    frames2.append(
        df.assign(
            w=df.w.cat.set_categories(list("xyz")),
            y=df.y.cat.set_categories(list("abc")),
        )
    )
frames3 = [i.set_index(i.y) for i in frames]
frames4 = [i.set_index(i.y) for i in frames2]
frames5 = [i.set_index([i.y, i.x]) for i in frames]
frames6 = [i.set_index([i.y, i.x]) for i in frames2]


def test_concat_unions_categoricals():
    # Categorical DataFrame, regular index
    tm.assert_frame_equal(_concat(frames), pd.concat(frames2))

    # Categorical Series, regular index
    tm.assert_series_equal(
        _concat([i.y for i in frames]), pd.concat([i.y for i in frames2])
    )

    # Categorical Index
    tm.assert_index_equal(
        _concat([i.index for i in frames3]), pd.concat([i for i in frames4]).index
    )

    # Categorical DataFrame, Categorical Index
    tm.assert_frame_equal(_concat(frames3), pd.concat(frames4))

    # Non-categorical DataFrame, Categorical Index
    tm.assert_frame_equal(
        _concat([i[["x", "z"]] for i in frames3]),
        pd.concat([i[["x", "z"]] for i in frames4]),
    )

    # Categorical Series, Categorical Index
    tm.assert_series_equal(
        _concat([i.z for i in frames3]), pd.concat([i.z for i in frames4])
    )

    # Non-categorical Series, Categorical Index
    tm.assert_series_equal(
        _concat([i.x for i in frames3]), pd.concat([i.x for i in frames4])
    )

    # MultiIndex with Categorical Index
    tm.assert_index_equal(
        _concat([i.index for i in frames5]), pd.concat([i for i in frames6]).index
    )

    # DataFrame, MultiIndex with CategoricalIndex
    tm.assert_frame_equal(_concat(frames5), pd.concat(frames6))


@pytest.mark.gpu
def test_unknown_categories_cudf():
    # We should always start with unknown categories
    # if `clear_known_categories` is working.
    pytest.importorskip("dask_cudf")

    with dask.config.set({"dataframe.backend": "cudf"}):
        ddf = dd.from_dict({"a": [0, 1, 0]}, npartitions=1)
    ddf["a"] = ddf["a"].astype("category")
    assert not ddf["a"].cat.known


# TODO: Remove the filterwarnings below
@pytest.mark.parametrize(
    "numeric_only",
    [
        True,
        pytest.param(
            False,
            marks=[
                pytest.mark.xfail(reason="numeric_only=False not implemented"),
            ],
        ),
        pytest.param(
            None,
            marks=pytest.mark.xfail(reason="numeric_only=False not implemented"),
        ),
    ],
)
@pytest.mark.parametrize("npartitions", [None, 10])
@pytest.mark.parametrize("split_out", [1, 4])
@pytest.mark.filterwarnings("ignore:The default value of numeric_only")
@pytest.mark.filterwarnings("ignore:Dropping")
def test_unknown_categoricals(
    shuffle_method, numeric_only, npartitions, split_out, request
):
    dsk = {("unknown", i): df for (i, df) in enumerate(frames)}
    meta = {"v": "object", "w": "category", "x": "i8", "y": "category", "z": "f8"}
    ddf = dd.from_pandas(pd.concat(dsk.values()).astype(meta), npartitions=4)
    if npartitions is not None:
        ddf = ddf.repartition(npartitions=npartitions)
    # Compute
    df = ddf.compute()

    assert_eq(ddf.w.value_counts(), df.w.value_counts())
    assert_eq(ddf.w.nunique(), df.w.nunique())

    ctx = (
        pytest.warns(FutureWarning, match="The default of observed=False")
        if PANDAS_GE_210 and not PANDAS_GE_300
        else contextlib.nullcontext()
    )
    numeric_kwargs = {} if numeric_only is None else {"numeric_only": numeric_only}
    with ctx:
        expected = df.groupby(df.w).sum(**numeric_kwargs)
    with ctx:
        result = ddf.groupby(ddf.w).sum(**numeric_kwargs)
    assert_eq(result, expected)

    with ctx:
        expected = df.groupby(df.w).y.nunique()
    with ctx:
        result = ddf.groupby(ddf.w, sort=False).y.nunique(split_out=split_out)
    assert_eq(result, expected)

    with ctx:
        expected = df.y.groupby(df.w).count()
    with ctx:
        result = ddf.y.groupby(ddf.w).count()
    assert_eq(result, expected)


def test_categorize():
    # rename y to y_ to avoid pandas future warning about ambiguous
    # levels
    pdf = frames4[0]
    if pyarrow_strings_enabled():
        # we explicitly provide meta, so it has to have pyarrow strings
        pdf = to_pyarrow_string(pdf)
    dsk = {("unknown", i): df for (i, df) in enumerate(frames3)}
    pdf = (
        pd.concat(dsk.values())
        .rename(columns={"y": "y_"})
        .astype({"w": "category", "y_": "category"})
    )
    pdf.index = pdf.index.astype("category")
    ddf = dd.from_pandas(pdf, npartitions=4, sort=False)
    ddf["w"] = ddf.w.cat.as_unknown()
    ddf["y_"] = ddf.y_.cat.as_unknown()
    ddf.index = ddf.index.cat.as_unknown()

    ddf = ddf.assign(w=ddf.w.cat.set_categories(["x", "y", "z"]))
    assert ddf.w.cat.known
    assert not ddf.y_.cat.known
    assert not ddf.index.cat.known
    df = ddf.compute()

    for index in [None, True, False]:
        known_index = index is not False
        # By default categorize object and unknown cat columns
        ddf2 = ddf.categorize(index=index)
        assert ddf2.y_.cat.known
        assert ddf2.v.cat.known
        assert ddf2.index.cat.known == known_index
        assert_eq(ddf2, df.astype({"v": "category"}), check_categorical=False)

        # Specifying split_every works
        ddf2 = ddf.categorize(index=index, split_every=2)
        assert ddf2.y_.cat.known
        assert ddf2.v.cat.known
        assert ddf2.index.cat.known == known_index
        assert_eq(ddf2, df.astype({"v": "category"}), check_categorical=False)

        # Specifying one column doesn't affect others
        ddf2 = ddf.categorize("v", index=index)
        assert not ddf2.y_.cat.known
        assert ddf2.v.cat.known
        assert ddf2.index.cat.known == known_index
        assert_eq(ddf2, df.astype({"v": "category"}), check_categorical=False)

        ddf2 = ddf.categorize("y_", index=index)
        assert ddf2.y_.cat.known
        assert ddf2.v.dtype == get_string_dtype()
        assert ddf2.index.cat.known == known_index
        assert_eq(ddf2, df)

    ddf_known_index = ddf.categorize(columns=[], index=True)
    assert ddf_known_index.index.cat.known
    assert_eq(ddf_known_index, df)

    # Specifying known categorical or no columns is a no-op:
    assert ddf.categorize(["w"], index=False) is ddf
    assert ddf.categorize([], index=False) is ddf
    assert ddf_known_index.categorize(["w"]) is ddf_known_index
    assert ddf_known_index.categorize([]) is ddf_known_index

    # Bad split_every fails
    with pytest.raises(ValueError):
        ddf.categorize(split_every=1)

    with pytest.raises(ValueError):
        ddf.categorize(split_every="foo")


def test_categorical_dtype():
    cat_dtype = dd.categorical.categorical_dtype(
        meta=a, categories=["a", "b", "c"], ordered=False
    )
    assert_eq(cat_dtype.categories, pd.Index(["a", "b", "c"]))
    assert_eq(cat_dtype.ordered, False)

    cat_dtype = dd.categorical.categorical_dtype(meta=a, categories=["a", "b", "c"])
    assert_eq(cat_dtype.categories, pd.Index(["a", "b", "c"]))
    assert_eq(cat_dtype.ordered, False)

    cat_dtype = dd.categorical.categorical_dtype(
        meta=a, categories=[1, 100, 200], ordered=True
    )
    assert_eq(cat_dtype.categories, pd.Index([1, 100, 200]))
    assert_eq(cat_dtype.ordered, True)


def test_categorize_index():
    # Object dtype
    pdf = _compat.makeDataFrame()
    ddf = dd.from_pandas(pdf, npartitions=5)
    result = ddf.compute()

    ddf2 = ddf.categorize()
    assert ddf2.index.cat.known
    assert_eq(
        ddf2,
        result.set_index(pd.CategoricalIndex(result.index)),
        check_divisions=False,
        check_categorical=False,
    )

    assert ddf.categorize(index=False) is ddf

    # Non-object dtype
    ddf = dd.from_pandas(result.set_index(result.A.rename("idx")), npartitions=5)
    result = ddf.compute()

    ddf2 = ddf.categorize(index=True)
    assert ddf2.index.cat.known
    assert_eq(
        ddf2,
        result.set_index(pd.CategoricalIndex(result.index)),
        check_divisions=False,
        check_categorical=False,
    )

    assert ddf.categorize() is ddf


def test_categorical_set_index(shuffle_method):
    df = pd.DataFrame({"x": [1, 2, 3, 4], "y": ["a", "b", "b", "c"]})
    df["y"] = pd.Categorical(df["y"], categories=["a", "b", "c"], ordered=True)
    a = dd.from_pandas(df, npartitions=2)

    with dask.config.set(scheduler="sync"):
        b = a.set_index("y", npartitions=a.npartitions)
        d1, d2 = b.get_partition(0), b.get_partition(1)
        assert list(d1.index.compute()) == ["a"]
        assert list(sorted(d2.index.compute())) == ["b", "b", "c"]

        b = a.set_index(a.y, npartitions=a.npartitions)
        d1, d2 = b.get_partition(0), b.get_partition(1)
        assert list(d1.index.compute()) == ["a"]
        assert list(sorted(d2.index.compute())) == ["b", "b", "c"]

        b = a.set_index("y", divisions=["a", "b", "c"], npartitions=a.npartitions)
        d1, d2 = b.get_partition(0), b.get_partition(1)
        assert list(d1.index.compute()) == ["a"]
        assert list(sorted(d2.index.compute())) == ["b", "b", "c"]


@pytest.mark.parametrize("ncategories", [1, 3, 6])
@pytest.mark.parametrize("npartitions", [1, 3, 6])
def test_categorical_set_index_npartitions_vs_ncategories(npartitions, ncategories):
    """https://github.com/dask/dask/issues/5343"""
    rows_per_category = 10
    n_rows = ncategories * rows_per_category

    categories = ["CAT" + str(i) for i in range(ncategories)]
    pdf = pd.DataFrame(
        {"id": categories * rows_per_category, "value": np.random.random(n_rows)}
    )
    ddf = dd.from_pandas(pdf, npartitions=npartitions)
    ddf["id"] = ddf["id"].astype("category").cat.as_ordered()
    ddf = ddf.set_index("id")
    # Test passes if this worked and didn't raise any warnings


@pytest.mark.parametrize("npartitions", [1, 4])
def test_repartition_on_categoricals(npartitions):
    df = pd.DataFrame({"x": range(10), "y": list("abababcbcb")})
    if pyarrow_strings_enabled():
        # we need this because a CategoricalDtype backed by arrow strings
        # is not the same as CategoricalDtype backed by object strings
        df = to_pyarrow_string(df)

    ddf = dd.from_pandas(df, npartitions=2)
    ddf["y"] = ddf["y"].astype("category")
    ddf2 = ddf.repartition(npartitions=npartitions)

    df = df.copy()
    df["y"] = df["y"].astype("category")
    assert_eq(df, ddf)
    assert_eq(df, ddf2)


def test_categorical_accessor_presence():
    df = pd.DataFrame({"x": list("a" * 5 + "b" * 5 + "c" * 5), "y": range(15)})
    df.x = df.x.astype("category")
    ddf = dd.from_pandas(df, npartitions=2)

    assert "cat" in dir(ddf.x)
    assert "cat" not in dir(ddf.y)
    assert hasattr(ddf.x, "cat")
    assert not hasattr(ddf.y, "cat")

    df2 = df.set_index(df.x)
    ddf2 = dd.from_pandas(df2, npartitions=2, sort=False)
    assert hasattr(ddf2.index, "categories")
    assert not hasattr(ddf.index, "categories")


def test_categorize_nan():
    df = dd.from_pandas(
        pd.DataFrame({"A": ["a", "b", "a", float("nan")]}), npartitions=2
    )
    with warnings.catch_warnings(record=True) as record:
        df.categorize().compute()
    assert not record


def get_cat(x):
    return x if isinstance(x, pd.CategoricalIndex) else x.cat


def assert_array_index_eq(left, right, check_divisions=False):
    """left and right are equal, treating index and array as equivalent"""
    assert_eq(
        left,
        pd.Index(right) if isinstance(right, np.ndarray) else right,
        check_divisions=check_divisions,
    )


def test_return_type_known_categories():
    df = pd.DataFrame({"A": ["a", "b", "c"]})
    df["A"] = df["A"].astype("category")
    dask_df = dd.from_pandas(df, 2)
    ret_type = dask_df.A.cat.as_known()
    assert isinstance(ret_type, dd.Series)


class TestCategoricalAccessor:
    @pytest.mark.parametrize("series", cat_series)
    @pytest.mark.parametrize(
        "prop, compare",
        [
            ("categories", assert_array_index_eq),
            ("ordered", assert_eq),
            ("codes", assert_array_index_eq),
        ],
    )
    def test_properties(self, series, prop, compare):
        s, ds = series
        expected = getattr(get_cat(s), prop)
        result = getattr(get_cat(ds), prop)
        compare(result, expected, check_divisions=False)

    @pytest.mark.parametrize("series", cat_series)
    @pytest.mark.parametrize(
        "method, kwargs",
        [
            ("add_categories", dict(new_categories=["d", "e"])),
            ("as_ordered", {}),
            ("as_unordered", {}),
            ("as_ordered", {}),
            ("remove_categories", dict(removals=["a"])),
            ("rename_categories", dict(new_categories=["d", "e", "f"])),
            ("reorder_categories", dict(new_categories=["a", "b", "c"])),
            ("set_categories", dict(new_categories=["a", "e", "b"])),
            ("remove_unused_categories", {}),
        ],
    )
    def test_callable(self, series, method, kwargs):
        op = operator.methodcaller(method, **kwargs)

        # Series
        s, ds = series
        expected = op(get_cat(s))
        result = op(get_cat(ds))
        assert_eq(result, expected, check_divisions=False)
        assert_eq(
            get_cat(result._meta).categories,
            get_cat(expected).categories,
            check_divisions=False,
        )
        assert_eq(
            get_cat(result._meta).ordered,
            get_cat(expected).ordered,
            check_divisions=False,
        )

    def test_categorical_empty(self):
        # GH 1705

        def make_empty():
            return pd.DataFrame({"A": pd.Categorical([np.nan, np.nan])})

        def make_full():
            return pd.DataFrame({"A": pd.Categorical(["a", "a"])})

        a = dd.from_delayed([dask.delayed(make_empty)(), dask.delayed(make_full)()])
        # Used to raise an IndexError
        a.A.cat.categories

    @pytest.mark.parametrize("series", cat_series)
    def test_unknown_categories(self, series):
        a, da = series
        assert da.cat.known
        da = da.cat.as_unknown()
        assert not da.cat.known

        with pytest.raises(NotImplementedError, match="with unknown categories"):
            da.cat.categories
        with pytest.raises(NotImplementedError, match="with unknown categories"):
            da.cat.codes
        # Also AttributeError so glob searching in IPython such as `da.cat.*?` works
        with pytest.raises(AttributeError, match="with unknown categories"):
            da.cat.categories
        with pytest.raises(AttributeError, match="with unknown categories"):
            da.cat.codes

        db = da.cat.set_categories(["a", "b", "c"])
        assert db.cat.known
        tm.assert_index_equal(db.cat.categories, get_cat(a).categories)
        assert_array_index_eq(db.cat.codes, get_cat(a).codes)

        db = da.cat.as_known()
        assert db.cat.known
        res = db.compute()
        tm.assert_index_equal(db.cat.categories, get_cat(res).categories)
        assert_array_index_eq(db.cat.codes, get_cat(res).codes)

    def test_categorical_string_ops(self):
        a = pd.Series(["a", "a", "b"], dtype="category")
        da = dd.from_pandas(a, 2)
        result = da.str.upper()
        expected = a.str.upper()
        assert_eq(result, expected)

    def test_categorical_non_string_raises(self):
        a = pd.Series([1, 2, 3], dtype="category")
        da = dd.from_pandas(a, 2)
        with pytest.raises(AttributeError):
            da.str.upper()
