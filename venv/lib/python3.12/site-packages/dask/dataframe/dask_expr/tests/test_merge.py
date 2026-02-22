from __future__ import annotations

import warnings

import numpy as np
import pytest

from dask import delayed
from dask.dataframe._compat import PANDAS_GE_220
from dask.dataframe.dask_expr import (
    Merge,
    from_delayed,
    from_dict,
    from_pandas,
    merge,
    repartition,
)
from dask.dataframe.dask_expr._expr import Filter, Projection
from dask.dataframe.dask_expr._merge import BroadcastJoin
from dask.dataframe.dask_expr._shuffle import Shuffle
from dask.dataframe.dask_expr.io import FromPandas
from dask.dataframe.dask_expr.tests._util import _backend_library, assert_eq

# Set DataFrame backend for this module
pd = _backend_library()


@pytest.mark.parametrize("how", ["left", "right", "inner", "outer"])
@pytest.mark.parametrize("shuffle_method", ["tasks", "disk"])
def test_merge(how, shuffle_method):
    # Make simple left & right dfs
    pdf1 = pd.DataFrame({"x": range(20), "y": range(20)})
    df1 = from_pandas(pdf1, 4)
    pdf2 = pd.DataFrame({"x": range(0, 20, 2), "z": range(10)})
    df2 = from_pandas(pdf2, 2)

    # Partition-wise merge with map_partitions
    df3 = df1.merge(df2, on="x", how=how, shuffle_method=shuffle_method)

    # Check result with/without fusion
    expect = pdf1.merge(pdf2, on="x", how=how)
    assert_eq(df3, expect, check_index=False)
    assert_eq(df3.optimize(), expect, check_index=False)

    df3 = merge(df1, df2, on="x", how=how, shuffle_method=shuffle_method)
    assert_eq(df3, expect, check_index=False)
    assert_eq(df3.optimize(), expect, check_index=False)


@pytest.mark.parametrize("how", ["left", "right", "inner", "outer"])
@pytest.mark.parametrize("pass_name", [True, False])
@pytest.mark.parametrize("sort", [True, False])
@pytest.mark.parametrize("shuffle_method", ["tasks", "disk"])
def test_merge_indexed(how, pass_name, sort, shuffle_method):
    # Make simple left & right dfs
    pdf1 = pd.DataFrame({"x": range(20), "y": range(20)}).set_index("x")
    df1 = from_pandas(pdf1, 4)
    pdf2 = pd.DataFrame({"x": range(0, 20, 2), "z": range(10)}).set_index("x")
    df2 = from_pandas(pdf2, 2, sort=sort)

    if pass_name:
        left_on = right_on = "x"
        left_index = right_index = False
    else:
        left_on = right_on = None
        left_index = right_index = True

    df3 = df1.merge(
        df2,
        left_index=left_index,
        left_on=left_on,
        right_index=right_index,
        right_on=right_on,
        how=how,
        shuffle_method=shuffle_method,
    )

    # Check result with/without fusion
    expect = pdf1.merge(
        pdf2,
        left_index=left_index,
        left_on=left_on,
        right_index=right_index,
        right_on=right_on,
        how=how,
    )
    assert_eq(df3, expect)
    assert_eq(df3.optimize(), expect)


@pytest.mark.parametrize("how", ["left", "right", "inner", "outer"])
@pytest.mark.parametrize("npartitions", [None, 22])
def test_broadcast_merge(how, npartitions):
    # Make simple left & right dfs
    pdf1 = pd.DataFrame({"x": range(40), "y": range(40)})
    df1 = from_pandas(pdf1, 20)
    pdf2 = pd.DataFrame({"x": range(0, 40, 2), "z": range(20)})
    df2 = from_pandas(pdf2, 2)

    df3 = df1.merge(
        df2, on="x", how=how, npartitions=npartitions, shuffle_method="tasks"
    )
    if npartitions:
        assert df3.npartitions == npartitions

    # Check that we avoid the shuffle when allowed
    if how in ("left", "inner"):
        assert all(["Shuffle" not in str(op) for op in df3.simplify().operands[:2]])

    # Check result with/without fusion
    expect = pdf1.merge(pdf2, on="x", how=how)
    assert_eq(df3, expect, check_index=False)
    assert_eq(df3.optimize(), expect, check_index=False)


def test_merge_column_projection():
    # Make simple left & right dfs
    pdf1 = pd.DataFrame({"x": range(20), "y": range(20), "z": range(20)})
    df1 = from_pandas(pdf1, 4)
    pdf2 = pd.DataFrame({"x": range(0, 20, 2), "z": range(10)})
    df2 = from_pandas(pdf2, 2)

    # Partition-wise merge with map_partitions
    df3 = df1.merge(df2, on="x")["z_x"].simplify()

    assert "y" not in df3.expr.operands[0].columns


@pytest.mark.parametrize("how", ["left", "right", "inner", "outer"])
@pytest.mark.parametrize("shuffle_method", ["tasks", "disk"])
def test_join(how, shuffle_method):
    # Make simple left & right dfs
    pdf1 = pd.DataFrame({"x": range(20), "y": range(20)})
    df1 = from_pandas(pdf1, 4)
    pdf2 = pd.DataFrame({"z": range(10)}, index=pd.Index(range(10), name="a"))
    df2 = from_pandas(pdf2, 2)

    # Partition-wise merge with map_partitions
    df3 = df1.join(df2, on="x", how=how, shuffle_method=shuffle_method)

    # Check result with/without fusion
    expect = pdf1.join(pdf2, on="x", how=how)
    assert_eq(df3.compute(), expect, check_index=False)
    assert_eq(df3.optimize(), expect, check_index=False)

    df3 = df1.join(df2.z, on="x", how=how, shuffle_method=shuffle_method)
    assert_eq(df3, expect, check_index=False)
    assert_eq(df3.optimize(), expect, check_index=False)


def test_join_recursive():
    pdf = pd.DataFrame({"x": [1, 2, 3], "y": 1}, index=pd.Index([1, 2, 3], name="a"))
    df = from_pandas(pdf, npartitions=2)

    pdf2 = pd.DataFrame(
        {"a": [1, 2, 3, 4, 5, 6], "b": 1}, index=pd.Index([1, 2, 3, 4, 5, 6], name="a")
    )
    df2 = from_pandas(pdf2, npartitions=2)

    pdf3 = pd.DataFrame({"c": [1, 2, 3], "d": 1}, index=pd.Index([1, 2, 3], name="a"))
    df3 = from_pandas(pdf3, npartitions=2)

    result = df.join([df2, df3], how="outer")
    assert_eq(result, pdf.join([pdf2, pdf3], how="outer"))

    result = df.join([df2, df3], how="left")
    # The nature of our join might cast ints to floats
    assert_eq(result, pdf.join([pdf2, pdf3], how="left"), check_dtype=False)


def test_join_recursive_raises():
    pdf = pd.DataFrame({"x": [1, 2, 3], "y": 1}, index=pd.Index([1, 2, 3], name="a"))
    df = from_pandas(pdf, npartitions=2)
    with pytest.raises(ValueError, match="other must be DataFrame"):
        df.join(["dummy"])

    with pytest.raises(ValueError, match="only supports left or outer"):
        df.join([df], how="inner")
    with pytest.raises(ValueError, match="only supports left or outer"):
        df.join([df], how="right")


def test_singleton_divisions():
    df = pd.DataFrame({"x": [1, 1, 1]}, index=[1, 2, 3])
    ddf = from_pandas(df, npartitions=2)
    ddf2 = ddf.set_index("x")

    joined = ddf2.join(ddf2, rsuffix="r")
    assert joined.divisions == (1, 1)
    joined.compute()


def test_categorical_merge_does_not_increase_npartitions():
    df1 = pd.DataFrame(data={"A": ["a", "b", "c"]}, index=["s", "v", "w"])
    df2 = pd.DataFrame(data={"B": ["t", "d", "i"]}, index=["v", "w", "r"])
    # We are npartitions=1 on both sides, so it should stay that way
    ddf1 = from_pandas(df1, npartitions=1)
    df2 = df2.astype({"B": "category"})
    assert_eq(df1.join(df2), ddf1.join(df2))


def test_merge_len():
    pdf = pd.DataFrame({"x": [1, 2, 3], "y": 1})
    df = from_pandas(pdf, npartitions=2)
    pdf2 = pd.DataFrame({"x": [1, 2, 3], "z": 1})
    df2 = from_pandas(pdf2, npartitions=2)
    assert_eq(len(df.merge(df2)), len(pdf.merge(pdf2)))
    query = df.merge(df2).index.optimize(fuse=False)
    expected = df[["x"]].merge(df2[["x"]]).index.optimize(fuse=False)
    assert query._name == expected._name


def test_merge_optimize_subset_strings():
    pdf = pd.DataFrame({"a": [1, 2], "aaa": 1})
    pdf2 = pd.DataFrame({"b": [1, 2], "aaa": 1})
    df = from_pandas(pdf)
    df2 = from_pandas(pdf2)

    query = df.merge(df2, on="aaa")[["aaa"]].optimize(fuse=False)
    exp = df[["aaa"]].merge(df2[["aaa"]], on="aaa").optimize(fuse=False)
    assert query._name == exp._name
    assert_eq(query, pdf.merge(pdf2, on="aaa")[["aaa"]])


@pytest.mark.parametrize("npartitions_left, npartitions_right", [(2, 3), (1, 1)])
def test_merge_combine_similar(npartitions_left, npartitions_right):
    pdf = pd.DataFrame(
        {
            "a": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
            "b": 1,
            "c": 1,
            "d": 1,
            "e": 1,
            "f": 1,
        }
    )
    pdf2 = pd.DataFrame({"a": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10], "x": 1})

    df = from_pandas(pdf, npartitions=npartitions_left)
    df2 = from_pandas(pdf2, npartitions=npartitions_right)

    query = df.merge(df2)
    query["new"] = query.b + query.c
    query = query.groupby(["a", "e", "x"]).new.sum()
    assert (
        len(query.optimize().__dask_graph__()) <= 30
    )  # 45 is the non-combined version

    expected = pdf.merge(pdf2)
    expected["new"] = expected.b + expected.c
    expected = expected.groupby(["a", "e", "x"]).new.sum()
    assert_eq(query, expected)


def test_merge_combine_similar_intermediate_projections():
    pdf = pd.DataFrame(
        {
            "a": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
            "b": 1,
            "c": 1,
        }
    )
    pdf2 = pd.DataFrame({"a": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10], "x": 1})
    pdf3 = pd.DataFrame({"d": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10], "e": 1, "y": 1})

    df = from_pandas(pdf, npartitions=2)
    df2 = from_pandas(pdf2, npartitions=3)
    df3 = from_pandas(pdf3, npartitions=3)

    q = df.merge(df2).merge(df3, left_on="b", right_on="d")[["b", "x", "y"]]
    q["new"] = q.b + q.x
    result = q.optimize(fuse=False)
    # Check that we have intermediate projections dropping unnecessary columns
    assert isinstance(result.expr.frame, Projection)
    assert isinstance(result.expr.frame.frame, Merge)
    assert isinstance(result.expr.frame.frame.left, Projection)
    assert isinstance(result.expr.frame.frame.left.frame, Shuffle)

    pd_result = pdf.merge(pdf2).merge(pdf3, left_on="b", right_on="d")[["b", "x", "y"]]
    pd_result["new"] = pd_result.b + pd_result.x

    assert sorted(result.expr.frame.frame.left.operand("columns")) == ["b", "x"]
    assert_eq(result, pd_result, check_index=False)


def test_categorical_merge_with_merge_column_cat_in_one_and_not_other_upcasts():
    df1 = pd.DataFrame({"A": pd.Categorical([0, 1]), "B": pd.Categorical(["a", "b"])})
    df2 = pd.DataFrame({"C": pd.Categorical(["a", "b"])})

    expected = pd.merge(df2, df1, left_index=True, right_on="A")

    ddf1 = from_pandas(df1, npartitions=2)
    ddf2 = from_pandas(df2, npartitions=2)

    actual = merge(ddf2, ddf1, left_index=True, right_on="A").compute()
    assert actual.C.dtype == "category"
    assert actual.B.dtype == "category"
    assert actual.A.dtype == "int64"
    assert actual.index.dtype == "int64"
    assert assert_eq(expected, actual)


def test_merge_combine_similar_hangs():
    var1 = 15
    var2 = "BRASS"
    var3 = "EUROPE"
    region_ds = from_pandas(
        pd.DataFrame.from_dict(
            {
                "r_regionkey": {0: 0, 1: 1},
                "r_name": {0: "AFRICA", 1: "AMERICA"},
                "r_comment": {0: "a", 1: "s "},
            }
        )
    )
    nation_filtered = from_pandas(
        pd.DataFrame.from_dict(
            {
                "n_nationkey": {0: 0, 1: 1},
                "n_name": {0: "ALGERIA", 1: "ARGENTINA"},
                "n_regionkey": {0: 0, 1: 1},
                "n_comment": {0: "fu", 1: "i"},
            }
        )
    )

    supplier_filtered = from_pandas(
        pd.DataFrame.from_dict(
            {
                "s_suppkey": {0: 1, 1: 2},
                "s_name": {0: "a#1", 1: "a#2"},
                "s_address": {0: "sdrGnX", 1: "T"},
                "s_nationkey": {0: 17, 1: 5},
                "s_phone": {0: "27-918-335-1736", 1: "15-679-861-2259"},
                "s_acctbal": {0: 5755, 1: 4032},
                "s_comment": {0: " inst", 1: " th"},
            }
        )
    )
    part_filtered = from_pandas(
        pd.DataFrame.from_dict(
            {
                "p_partkey": {0: 1, 1: 2},
                "p_name": {0: "gol", 1: "bl"},
                "p_mfgr": {0: "Manufacturer#1", 1: "Manufacturer#1"},
                "p_brand": {0: "Brand#13", 1: "Brand#13"},
                "p_type": {0: "PROM", 1: "LARG"},
                "p_size": {0: 7, 1: 1},
                "p_container": {0: "J", 1: "LG"},
                "p_retailprice": {0: 901, 1: 902},
                "p_comment": {0: "ir", 1: "ack"},
            }
        )
    )
    #
    partsupp_filtered = from_pandas(
        pd.DataFrame.from_dict(
            {
                "ps_partkey": {0: 1, 1: 1},
                "ps_suppkey": {0: 2, 1: 2502},
                "ps_availqty": {0: 3325, 1: 8076},
                "ps_supplycost": {0: 771, 1: 993},
                "ps_comment": {0: "bli", 1: "ts boo"},
            }
        )
    )

    region_filtered = region_ds[(region_ds["r_name"] == var3)]
    r_n_merged = nation_filtered.merge(
        region_filtered, left_on="n_regionkey", right_on="r_regionkey", how="inner"
    )
    s_r_n_merged = r_n_merged.merge(
        supplier_filtered,
        left_on="n_nationkey",
        right_on="s_nationkey",
        how="inner",
    )
    ps_s_r_n_merged = s_r_n_merged.merge(
        partsupp_filtered, left_on="s_suppkey", right_on="ps_suppkey", how="inner"
    )
    part_filtered = part_filtered[
        (part_filtered["p_size"] == var1)
        & (part_filtered["p_type"].astype(str).str.endswith(var2))
    ]
    merged_df = part_filtered.merge(
        ps_s_r_n_merged, left_on="p_partkey", right_on="ps_partkey", how="inner"
    )
    min_values = merged_df.groupby("p_partkey")["ps_supplycost"].min().reset_index()
    min_values.columns = ["P_PARTKEY_CPY", "MIN_SUPPLYCOST"]
    merged_df = merged_df.merge(
        min_values,
        left_on=["p_partkey", "ps_supplycost"],
        right_on=["P_PARTKEY_CPY", "MIN_SUPPLYCOST"],
        how="inner",
    )
    out = merged_df[
        [
            "s_acctbal",
            "s_name",
            "n_name",
            "p_partkey",
            "p_mfgr",
            "s_address",
            "s_phone",
            "s_comment",
        ]
    ]
    expected = pd.DataFrame(
        columns=[
            "s_acctbal",
            "s_name",
            "n_name",
            "p_partkey",
            "p_mfgr",
            "s_address",
            "s_phone",
            "s_comment",
        ]
    )
    assert_eq(out, expected, check_dtype=False)

    # Double check that these don't hang
    out.optimize(fuse=False)
    out.optimize()


def test_recursive_join():
    dfs_to_merge = []
    for i in range(10):
        df = pd.DataFrame(
            {
                f"{i}A": [5, 6, 7, 8],
                f"{i}B": [4, 3, 2, 1],
            },
            index=pd.Index([0, 1, 2, 3], name="a"),
        )
        ddf = from_pandas(df, 2)
        dfs_to_merge.append(ddf)

    ddf_loop = from_pandas(pd.DataFrame(index=pd.Index([0, 1, 3], name="a")), 3)
    for ddf in dfs_to_merge:
        ddf_loop = ddf_loop.join(ddf, how="left")

    ddf_pairwise = from_pandas(pd.DataFrame(index=pd.Index([0, 1, 3], name="a")), 3)

    ddf_pairwise = ddf_pairwise.join(dfs_to_merge, how="left")

    assert_eq(ddf_pairwise, ddf_loop)


def test_merge_repartition():
    pdf = pd.DataFrame({"a": [1, 2, 3]})
    pdf2 = pd.DataFrame({"b": [1, 2, 3]}, index=[1, 2, 3])

    df = from_pandas(pdf, npartitions=2)
    df2 = from_pandas(pdf2, npartitions=3)
    assert_eq(df.join(df2), pdf.join(pdf2))


def test_merge_reparititon_divisions():
    pdf = pd.DataFrame({"a": [1, 2, 3, 4, 5, 6]})
    pdf2 = pd.DataFrame({"b": [1, 2, 3, 4, 5, 6]}, index=[1, 2, 3, 4, 5, 6])
    pdf3 = pd.DataFrame({"c": [1, 2, 3, 4, 5, 6]}, index=[1, 2, 3, 4, 5, 6])

    df = from_pandas(pdf, npartitions=2)
    df2 = from_pandas(pdf2, npartitions=3)
    df3 = from_pandas(pdf3, npartitions=3)

    assert_eq(df.join(df2).join(df3), pdf.join(pdf2).join(pdf3))


def test_join_gives_proper_divisions():
    df = pd.DataFrame({"a": ["a", "b", "c"]}, index=[0, 1, 2])
    ddf = from_pandas(df, npartitions=1)

    right_df = pd.DataFrame({"b": [1.0, 2.0, 3.0]}, index=["a", "b", "c"])

    expected = df.join(right_df, how="inner", on="a")
    actual = ddf.join(right_df, how="inner", on="a")
    assert actual.divisions == ddf.divisions

    assert_eq(expected, actual)


@pytest.mark.parametrize("shuffle_method", ["tasks", "disk"])
@pytest.mark.parametrize("how", ["inner", "left"])
def test_merge_known_to_single(how, shuffle_method):
    partition_sizes = np.array([3, 4, 2, 5, 3, 2, 5, 9, 4, 7, 4])
    idx = [i for i, s in enumerate(partition_sizes) for _ in range(s)]
    k = [i for s in partition_sizes for i in range(s)]
    vi = range(len(k))
    pdf1 = pd.DataFrame(dict(idx=idx, k=k, v1=vi)).set_index(["idx"])

    partition_sizes = np.array([4, 2, 5, 3, 2, 5, 9, 4, 7, 4, 8])
    idx = [i for i, s in enumerate(partition_sizes) for _ in range(s)]
    k = [i for s in partition_sizes for i in range(s)]
    vi = range(len(k))
    pdf2 = pd.DataFrame(dict(idx=idx, k=k, v1=vi)).set_index(["idx"])

    df1 = repartition(pdf1, [0, 1, 3, 4, 5, 6, 7, 8, 9, 10, 11])
    df2 = from_pandas(pdf2, npartitions=1, sort=False)

    expected = pdf1.merge(pdf2, on="idx", how=how)
    result = df1.merge(df2, on="idx", how=how, shuffle_method=shuffle_method)
    assert_eq(result, expected)
    assert result.divisions == df1.divisions

    expected = pdf1.merge(pdf2, on="k", how=how)
    result = df1.merge(df2, on="k", how=how, shuffle_method=shuffle_method)
    assert_eq(result, expected, check_index=False)
    assert all(d is None for d in result.divisions)


@pytest.mark.parametrize("how", ["right", "outer"])
def test_merge_empty_left_df(how):
    left = pd.DataFrame({"a": [1, 1, 2, 2], "val": [5, 6, 7, 8]})
    right = pd.DataFrame({"a": [0, 0, 3, 3], "val": [11, 12, 13, 14]})

    dd_left = from_pandas(left, npartitions=4)
    dd_right = from_pandas(right, npartitions=4)

    merged = dd_left.merge(dd_right, on="a", how=how)
    expected = left.merge(right, on="a", how=how)
    assert_eq(merged, expected, check_index=False)

    # Check that the individual partitions have the expected shape
    merged.map_partitions(lambda x: x, meta=merged._meta).compute()


def test_merge_npartitions():
    pdf = pd.DataFrame({"a": [1, 2, 3, 4, 5, 6]})
    pdf2 = pd.DataFrame({"b": [1, 2, 3, 4, 5, 6]}, index=[1, 2, 3, 4, 5, 6])
    df = from_pandas(pdf, npartitions=1)
    df2 = from_pandas(pdf2, npartitions=3)

    result = df.join(df2, npartitions=6)
    # Ignore npartitions when broadcasting
    assert result.npartitions == 4
    assert_eq(result, pdf.join(pdf2))

    df = from_pandas(pdf, npartitions=2)
    result = df.join(df2, npartitions=6)
    # Ignore npartitions for repartition-join
    assert result.npartitions == 4
    assert_eq(result, pdf.join(pdf2))

    pdf = pd.DataFrame(
        {"a": [1, 2, 3, 4, 5, 6]}, index=pd.Index([6, 5, 4, 3, 2, 1], name="a")
    )
    pdf2 = pd.DataFrame(
        {"b": [1, 2, 3, 4, 5, 6]}, index=pd.Index([1, 2, 7, 4, 5, 6], name="a")
    )
    df = from_pandas(pdf, npartitions=2, sort=False)
    df2 = from_pandas(pdf2, npartitions=3, sort=False)

    result = df.join(df2, npartitions=6)
    assert result.npartitions == 6
    assert_eq(result, pdf.join(pdf2))


@pytest.mark.parametrize("how", ["inner", "outer", "left", "right"])
@pytest.mark.parametrize("on_index", [True, False])
def test_merge_columns_dtypes1(how, on_index):
    # tests results of merges with merge columns having different dtypes;
    # asserts that either the merge was successful or the corresponding warning is raised
    # addresses issue #4574

    df1 = pd.DataFrame(
        {"A": list(np.arange(5).astype(float)) * 2, "B": list(np.arange(5)) * 2}
    )
    df2 = pd.DataFrame({"A": np.arange(5), "B": np.arange(5)})

    a = from_pandas(df1, 2)  # merge column "A" is float
    b = from_pandas(df2, 2)  # merge column "A" is int

    on = ["A"]
    left_index = right_index = on_index

    if on_index:
        a = a.set_index("A")
        b = b.set_index("A")
        on = None

    with warnings.catch_warnings(record=True) as record:
        warnings.simplefilter("always")
        result = merge(
            a, b, on=on, how=how, left_index=left_index, right_index=right_index
        )
        warned = any("merge column data type mismatches" in str(r) for r in record)

    # result type depends on merge operation -> convert to pandas
    result = result if isinstance(result, pd.DataFrame) else result.compute()

    has_nans = result.isna().values.any()
    assert (has_nans and warned) or not has_nans


def test_merge_pandas_object():
    pdf1 = pd.DataFrame({"x": range(20), "y": range(20)})
    df1 = from_pandas(pdf1, 4)
    pdf2 = pd.DataFrame({"x": range(20), "z": range(20)})

    assert_eq(merge(df1, pdf2, on="x"), pdf1.merge(pdf2, on="x"), check_index=False)
    assert_eq(merge(pdf2, df1, on="x"), pdf2.merge(pdf1, on="x"), check_index=False)

    pdf1 = pd.DataFrame({"x": range(20), "y": range(20)}).set_index("x")
    df1 = from_pandas(pdf1, 4)
    assert_eq(
        merge(df1, pdf2, left_index=True, right_on="x"),
        pdf1.merge(pdf2, left_index=True, right_on="x"),
        check_index=False,
    )
    assert_eq(
        merge(pdf2, df1, left_on="x", right_index=True),
        pdf2.merge(pdf1, left_on="x", right_index=True),
        check_index=False,
    )


@pytest.mark.parametrize("clear_divisions", [True, False])
@pytest.mark.parametrize("how", ["left", "outer"])
@pytest.mark.parametrize("npartitions_base", [1, 2, 3])
@pytest.mark.parametrize("npartitions_other", [1, 2, 3])
def test_pairwise_merge_results_in_identical_output_df(
    how, npartitions_base, npartitions_other, clear_divisions
):
    if clear_divisions and (npartitions_other != 3 or npartitions_base != 1):
        pytest.skip(reason="Runtime still slower than I would like, so save some time")
    dfs_to_merge = []
    for i in range(10):
        df = pd.DataFrame(
            {
                f"{i}A": [5, 6, 7, 8],
                f"{i}B": [4, 3, 2, 1],
            },
            index=[0, 1, 2, 3],
        )
        ddf = from_pandas(df, npartitions_other)
        if clear_divisions:
            ddf = ddf.clear_divisions()
        dfs_to_merge.append(ddf)

    ddf_loop = from_pandas(pd.DataFrame(index=[0, 1, 3]), npartitions_base)
    if clear_divisions:
        ddf_loop = ddf_loop.clear_divisions()
    for ddf in dfs_to_merge:
        ddf_loop = ddf_loop.join(ddf, how=how)

    ddf_pairwise = from_pandas(pd.DataFrame(index=[0, 1, 3]), npartitions_base)
    if clear_divisions:
        ddf_pairwise = ddf_pairwise.clear_divisions()

    ddf_pairwise = ddf_pairwise.join(dfs_to_merge, how=how)

    # recursive join doesn't yet respect divisions in dask-expr
    assert_eq(ddf_pairwise, ddf_loop)


def test_filter_merge():
    pdf_a = pd.DataFrame(
        {
            "a": range(5),
            "b": range(5),
            "c": range(5),
            "d": [True, False, True, False, True],
        }
    )
    pdf_b = pd.DataFrame(
        {
            "c": [0, 2, 4, 6, 8],
            "x": range(5),
            "y": range(5),
            "z": [False, False, True, True, True],
        }
    )

    a = from_pandas(pdf_a)
    b = from_pandas(pdf_b)

    # Some simple cases
    df = a.merge(b)
    # A simple projection
    df = df[df.z]
    bb = b[b.z]
    expected = a.merge(bb)
    assert df.optimize()._name == expected.optimize()._name
    assert_eq(df, expected)

    # Unary op
    df = a.merge(b)
    df = df[~df.z]
    bb = b[~b.z]
    expected = a.merge(bb)
    assert df.optimize()._name == expected.optimize()._name
    assert_eq(df, expected)

    df = a.merge(b)
    df = df[df.x > 3]
    bb = b[b.x > 3]
    expected = a.merge(bb)
    assert df.optimize()._name == expected.optimize()._name
    assert_eq(df, expected)

    df = a.merge(b)
    df = df[df.b > 3]
    aa = a[a.b > 3]
    expected = aa.merge(b)
    assert df.optimize()._name == expected.optimize()._name
    assert_eq(df, expected)

    df = a.merge(b)
    df = df[3 < df.b]
    aa = a[3 < a.b]
    expected = aa.merge(b)
    assert df.optimize()._name == expected.optimize()._name
    assert_eq(df, expected)

    # Apply to both!
    df = a.merge(b)
    df = df[df.c > 3]
    aa = a[a.c > 3]
    bb = b[b.c > 3]
    expected = aa.merge(bb)
    assert df.optimize()._name == expected.optimize()._name
    assert_eq(df, expected)

    # Works with more complex expressions, and multiple columns
    df = a.merge(b)
    df = df[df.a > df.b + 1]
    aa = a[a.a > a.b + 1]
    expected = aa.merge(b)
    assert df.optimize()._name == expected.optimize()._name
    assert_eq(df, expected)

    # Only apply if all columns are in the table, not if only some are
    df = a.merge(b)
    df = df[df.c > df.x + 1]
    bb = b[b.c > b.x + 1]
    expected = a.merge(bb)
    assert df.optimize()._name == expected.optimize()._name
    assert_eq(df, expected)

    df = a.merge(b)
    df = df[df.d & df.z]
    aa = a[a.d]
    bb = b[b.z]
    expected = aa.merge(bb)
    assert df.simplify()._name == expected.simplify()._name
    assert_eq(df, expected)

    df = a.merge(b)
    df = df[(df.a > 2) & df.z]
    aa = a[a.a > 2]
    bb = b[b.z]
    expected = aa.merge(bb)
    actual = df.optimize()
    assert actual._name == expected.optimize()._name
    assert_eq(df, expected)

    # Bail if we engage non-elemwise expressions in the predicates
    df = a.merge(b)
    df = df[df.x > df.y.sum()]
    bb = b[b.x > b.y.sum()]
    not_expected = a.merge(bb)
    assert df.optimize()._name != not_expected.optimize()._name

    df = a.merge(b)
    df = df[df.d & df.z][["a"]]
    aa = a[a.d]
    bb = b[b.z]
    expected = aa.merge(bb)[["a"]]
    assert df.simplify()._name == expected.simplify()._name
    assert_eq(df, expected)

    df = a.merge(b)
    df = df[df.d & df.z & (df.b == 1)][["a"]]
    aa = a[a.d & (a.b == 1)]
    bb = b[b.z]
    expected = aa.merge(bb)[["a"]]
    assert df.simplify()._name == expected.simplify()._name
    assert_eq(df, expected)

    df = a.merge(b)
    df = df[df.d == df.z]
    assert df.simplify()._name == df._name


def test_filter_merge_suffixes():
    pdf1 = pd.DataFrame({"a": [1, 2, 3, 4], "b": 1})
    pdf2 = pd.DataFrame({"a": [1, 2, 3, 4], "b": 2})
    df1 = from_pandas(pdf1, npartitions=2)
    df2 = from_pandas(pdf2, npartitions=2)
    q = df1.merge(df2, on="a", suffixes=("", "_right"))
    result = q[q.b < 2]
    expected = pdf1.merge(pdf2, on="a", suffixes=("", "_right"))
    assert_eq(result, expected[expected.b < 2], check_index=False)
    result = q[q.b > 1]
    expected = df1.merge(df2, on="a", suffixes=("", "_right"))
    assert_eq(result, expected[expected.b > 1], check_index=False)

    q = df1.merge(df2, on="a", suffixes=("_left", "_right"))
    result = q[q.b_left < 2]
    # Don't do anything for now
    assert result._name == result.simplify()._name
    expected = df1.merge(df2, on="a", suffixes=("_left", "_right"))
    assert_eq(result, expected[expected.b_left < 2], check_index=False)


def test_merge_filter_renaming_columns():
    pdf1 = pd.DataFrame({"a": [1, 2, 3, 4], "left": 4})
    df1 = from_pandas(pdf1, npartitions=2)

    pdf2 = pd.DataFrame({"a": [1, 2, 3, 4], "right": [6, 6, 5, 6]})
    df2 = from_pandas(pdf2, npartitions=2)

    result = df1.merge(df2)
    result = result[result[["right"]].rename(columns={"right": "left"})["left"] == 6]
    expected = pdf1.merge(pdf2)
    expected = expected[
        expected[["right"]].rename(columns={"right": "left"})["left"] == 6
    ]
    assert_eq(result, expected, check_index=False)

    q = result.simplify()
    assert isinstance(q.expr.right, Filter)
    assert isinstance(q.expr.left, FromPandas)


def test_merge_avoid_overeager_filter_pushdown():
    df = pd.DataFrame({"a": [1, 2, 3], "b": 1})
    ddf = from_pandas(df, npartitions=2)
    df2 = pd.DataFrame({"a": [2, 3, 4], "c": 1})
    ddf2 = from_pandas(df2, npartitions=2)
    merged = ddf.merge(ddf2, on="a", how="left")
    rhs = merged.c.sum()
    q = merged[merged.a > 1].assign(c=rhs)
    result = q.simplify()
    assert q._name == result._name
    assert isinstance(result.expr.frame, Filter)
    assert isinstance(result.expr.frame.frame, Merge)

    merged = ddf.merge(ddf2, on="a", how="left")
    rhs = merged.a.sum()
    q = merged[merged.a > 1].assign(c=rhs)
    result = q.simplify()
    assert q._name == result._name
    assert isinstance(result.expr.frame, Filter)
    assert isinstance(result.expr.frame.frame, Merge)


@pytest.mark.skipif(not PANDAS_GE_220, reason="pandas bug")
def test_join_consistent_index_names():
    pdf1 = pd.DataFrame(index=["a", "b", "c"], data=dict(a=[1, 2, 3]))
    pdf1.index.name = "test"
    df1 = from_pandas(pdf1, 2)

    pdf2 = pd.DataFrame(index=["a", "b", "d"], data=dict(b=[1, 2, 3]))
    df2 = from_pandas(pdf2, 2)
    result = df1.join(df2, how="outer")
    expected = pdf1.join(pdf2, how="outer")
    assert_eq(result, expected, check_index=False)
    assert result.index.name == result._meta.index.name
    assert result.partitions[0].compute().index.name == result._meta.index.name
    assert expected.index.name == result.index.name


@pytest.mark.parametrize("how", ["left", "inner", "right", "outer"])
def test_isin_filter_pushdown(how):
    pdf1 = pd.DataFrame(
        {
            "o_orderkey": [1, 2, 3, 4, 5],
        }
    )
    pdf2 = pd.DataFrame(
        {
            "l_orderkey": [1, 2, 3, 4, 5],
            "l_shipmode": ["MAIL", "SHIP", "RAIL", "bla", "MAIL"],
            "l_commitdate": [1, 2, 3, 4, 5],
            "l_receiptdate": [2, 2, 3, 5, 6],
        }
    )

    df1 = from_pandas(pdf1, npartitions=2)
    df2 = from_pandas(pdf2, npartitions=2)

    table = df1.merge(df2, left_on="o_orderkey", right_on="l_orderkey", how=how)
    result = table[
        (table.l_shipmode.isin(("MAIL", "SHIP")))
        & (table.l_commitdate < table.l_receiptdate)
    ]

    table = pdf1.merge(pdf2, left_on="o_orderkey", right_on="l_orderkey", how=how)
    expected = table[
        (table.l_shipmode.isin(("MAIL", "SHIP")))
        & (table.l_commitdate < table.l_receiptdate)
    ].sort_values(by="o_orderkey", ascending=False)
    assert_eq(result, expected, check_index=False)

    table = df2.merge(df1, left_on="l_orderkey", right_on="o_orderkey", how=how)
    result = table[
        (table.l_shipmode.isin(("MAIL", "SHIP")))
        & (table.l_commitdate < table.l_receiptdate)
    ]

    table = pdf2.merge(pdf1, left_on="l_orderkey", right_on="o_orderkey", how=how)
    expected = table[
        (table.l_shipmode.isin(("MAIL", "SHIP")))
        & (table.l_commitdate < table.l_receiptdate)
    ].sort_values(by="o_orderkey", ascending=False)
    assert_eq(result, expected, check_index=False)


def test_merge_filter_pushdown_broadcast():
    pdf = pd.DataFrame({"a": [1, 2, 3], "b": 1})
    pdf2 = pd.DataFrame({"c": [1, 2, 3], "b": 1})
    df1 = from_pandas(pdf, npartitions=2)
    df2 = from_pandas(pdf2, npartitions=2)
    result = df1.merge(df2, broadcast=True, shuffle_method="tasks")
    result = result[result.a > 1]
    result.optimize().pprint()
    assert len(list(result.optimize().find_operations(BroadcastJoin))) > 0


def test_merge_filter_stuck_between_merges():
    pdf = pd.DataFrame({"a": [1, 2, 3, 4, 5, 6, 7] * 100, "b": 1, "c": 2})
    pdf2 = pd.DataFrame({"a": [1, 2, 3, 4, 5, 6, 7] * 100, "d": 1, "e": 1, "f": 1})
    df1 = from_pandas(pdf, npartitions=2)
    df2 = from_pandas(pdf2, npartitions=2)
    result = df1.merge(df2)
    q = result.groupby("a").b.sum(split_out=True).reset_index()
    result = result.merge(q)
    result = result[(result.d == 1) & (result.e == 1)]
    e = df1.merge(df2)
    expected = e[(e.d == 1) & (e.e == 1)]
    expected_q = e.groupby("a").b.sum(split_out=True).reset_index()
    expected = expected.merge(expected_q)
    assert result.optimize()._name == expected.optimize()._name


def test_merge_or_two_branches():
    pdf = pd.DataFrame({"a": [1, 2, 3, 4, 5, 6, 7] * 100, "b": 1, "c": 2})
    pdf2 = pd.DataFrame({"a": [1, 2, 3, 4, 5, 6, 7] * 100, "d": 1, "e": 1, "f": 1})
    df1 = from_pandas(pdf, npartitions=2)
    df2 = from_pandas(pdf2, npartitions=2)

    q = df1.merge(df2)
    q = q[((q.b == 1) & (q.d == 1)) | ((q.b == 2) & (q.d == 1))]
    result = (q.c * (q.f)).sum()

    left = df1[(df1.b == 1) | (df1.b == 2)]
    right = df2[df2.d == 1]
    expected = left.merge(right)
    expected = (expected.c * (expected.f)).sum()
    assert result.optimize()._name == expected.optimize()._name


def test_merge_or_two_branches_no_pushdown():
    pdf = pd.DataFrame({"a": [1, 2, 3, 4, 5, 6, 7] * 100, "b": 1, "c": 2})
    pdf2 = pd.DataFrame({"a": [1, 2, 3, 4, 5, 6, 7] * 100, "d": 1, "e": 1, "f": 1})
    df1 = from_pandas(pdf, npartitions=2)
    df2 = from_pandas(pdf2, npartitions=2)

    q = df1.merge(df2)
    result = q[((q.b == 1) & (q.d == 1)) | ((q.b == 2) & (q.d == 1))]
    result = (result.c * (result.f)).sum() + q.a.sum()

    assert isinstance(result.expr.right.frame.frame, Merge)
    # Filter wasn't pushed through the merges
    assert isinstance(result.expr.left.frame.left.frame, Filter)
    assert isinstance(result.expr.left.frame.left.frame.frame, Merge)


def test_merge_filter_column_used_multiple_times():
    pdf = pd.DataFrame({"a": [1, 2, 3, 4, 5, 6, 7] * 100, "b": 1, "c": 2})
    pdf2 = pd.DataFrame({"a": [1, 2, 3, 4, 5, 6, 7] * 100, "d": 1, "e": 1, "f": 1})

    df1 = from_pandas(pdf, npartitions=2)
    df2 = from_pandas(pdf2, npartitions=2)
    result = df1.merge(df2)
    result = result[(result.b <= result.c) & (result.b == 1)]
    expected = df1[(df1.b <= df1.c) & (df1.b == 1)]
    expected = expected.merge(df2)
    assert result.optimize()._name == expected.optimize()._name


def test_merge_scalar_comparison():
    pdf = pd.DataFrame({"a": [1, 2, 3], "b": 1})
    pdf2 = pd.DataFrame({"c": [1, 2, 3], "b": 1})
    df = from_pandas(pdf, npartitions=2)
    df2 = from_pandas(pdf2, npartitions=2)
    result = df.merge(df2)
    result = result[result.a > df.a.mean()]
    expected = pdf.merge(pdf2)
    expected = expected[expected.a > pdf.a.mean()]
    assert_eq(result, expected, check_index=False)


def test_merge_leftsemi():
    pdf1 = pd.DataFrame({"aa": [1, 2, 3, 4, 5, 6, 1, 2, 3], "bb": 1})
    pdf2 = pd.DataFrame({"aa": [1, 2, 2, 4, 4, 10], "cc": 1})

    df1 = from_pandas(pdf1, npartitions=2)
    df2 = from_pandas(pdf2, npartitions=2)
    assert_eq(
        df1.merge(df2, how="leftsemi"),
        pdf1[pdf1.aa.isin(pdf2.aa)],
        check_index=False,
    )
    df2 = df2.rename(columns={"aa": "dd"})
    assert_eq(
        df1.merge(df2, how="leftsemi", left_on="aa", right_on="dd"),
        pdf1[pdf1.aa.isin(pdf2.aa)],
        check_index=False,
    )
    with pytest.raises(NotImplementedError, match="right_index=True"):
        df1.merge(df2, how="leftsemi")

    pdf2 = pdf2.set_index("aa")
    df2 = from_pandas(pdf2, npartitions=2)
    with pytest.raises(NotImplementedError, match="on columns from the index"):
        df1.merge(df2, how="leftsemi", on="aa")


def test_merge_suffix_projections():
    df = pd.DataFrame(
        {
            "a": [1, 2, 3] * 5,
            "b": [1, 2, 3] * 5,
            "c": ["A"] * 15,
        },
    )
    ddf = from_pandas(df, npartitions=1)

    merged_ddf = merge(ddf, ddf, on="a")
    result = merged_ddf[merged_ddf["c_x"] == "A"]["c_y"]
    expected = df.merge(df, on="a")
    expected = expected[expected["c_x"] == "A"]["c_y"]
    assert_eq(result, expected)


@pytest.mark.parametrize("index", [True, False])
def test_merge_after_rename(index):
    pleft = pd.Series(range(10))
    pleft = pleft.index if index else pleft
    pleft = pleft.drop_duplicates().to_frame()
    pleft.columns = ["a"]

    left = from_pandas(pd.Series(range(10)), npartitions=2)
    left = left.index if index else left
    left = left.drop_duplicates().to_frame()
    left.columns = ["a"]

    right = pd.DataFrame({"a": [1, 2] * 5})
    expected = pleft.merge(right, how="inner")
    result = left.merge(right, how="inner")
    assert_eq(result, expected, check_index=False)


def test_mismatching_meta():
    df1 = from_pandas(pd.DataFrame({"value": [1, 2, 3]}))

    def creator(i):
        return pd.DataFrame(
            {"value2": [1, 2, 3]}, index=pd.Index([1, 2, 3], name="index")
        )

    df2 = from_delayed(
        [delayed(creator)(i) for i in range(3)],
        meta=pd.DataFrame({"value2": [1, 2, 3]}),
    )
    result = df1.join(df2)
    assert_eq(result, df1.compute().join(df2.compute()), check_index=False)


def test_merge_tuple_left_on():
    df = pd.DataFrame(
        {
            "a": [1, 2, 3] * 5,
            "b": [1, 2, 3] * 5,
            "c": ["A"] * 15,
        },
    )
    ddf = from_pandas(df, npartitions=2)
    assert_eq(
        ddf.merge(ddf, left_on=("a",), right_on=("a",)),
        df.merge(df, left_on=("a",), right_on=("a",)),
        check_index=False,
    )
    assert_eq(
        ddf.merge(ddf, on=("a",)),
        df.merge(df, on=("a",)),
        check_index=False,
    )


def test_merged_partitions_filtered():
    a = from_dict(
        {"x": range(1000), "y": [1, 2, 3, 4] * 250}, npartitions=10
    ).partitions[:5]
    b = from_dict({"xx": range(100), "yy": [1, 2] * 50}, npartitions=3)
    result = a.merge(b, left_on=["y"], right_on=["yy"], how="inner", broadcast=True)

    # Check expression properties
    expr = result.optimize(fuse=False).expr
    assert not expr._filtered
    assert expr.left._filtered
    assert expr.divisions == expr._divisions()
    assert len(expr.divisions) == 6

    # Check result
    expect = a.compute().merge(b.compute(), left_on=["y"], right_on=["yy"], how="inner")
    assert_eq(result, expect, check_index=False)
