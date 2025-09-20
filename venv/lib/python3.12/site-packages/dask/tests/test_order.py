from __future__ import annotations

import inspect
from collections import defaultdict

import pytest
from packaging.version import Version

import dask
from dask import delayed
from dask._compatibility import WINDOWS
from dask.base import collections_to_expr, key_split, visualize_dsk
from dask.core import get_deps
from dask.order import _connecting_to_roots, diagnostics, ndependencies, order
from dask.utils_test import add, inc


@pytest.fixture(
    params=[
        "abcde",
        "edcba",
    ]
)
def abcde(request):
    return request.param


def f(*args):
    pass


def visualize(dsk, suffix="", **kwargs):
    """Utility to visualize the raw low level graphs in this tests suite. This
    automatically generates a set of visualizations with different metrics and
    writes them out to a file prefixed by the test name and suffixed by the
    measure used."""
    funcname = inspect.stack()[1][3] + suffix
    if hasattr(dsk, "__dask_graph__"):
        dsk = (
            collections_to_expr([dsk], optimize_graph=True).optimize().__dask_graph__()
        )

    node_attrs = {"penwidth": "6"}
    visualize_dsk(dsk, filename=f"{funcname}.pdf", node_attr=node_attrs, **kwargs)
    o = order(dsk, return_stats=True)
    visualize_dsk(
        dsk,
        o=o,
        filename=f"{funcname}-order.pdf",
        color="order",
        node_attr=node_attrs,
        cmap="viridis",
        collapse_outputs=False,
        **kwargs,
    )
    visualize_dsk(
        dsk,
        o=o,
        filename=f"{funcname}-age.pdf",
        color="age",
        node_attr=node_attrs,
        cmap="Reds",
        collapse_outputs=False,
        **kwargs,
    )
    visualize_dsk(
        dsk,
        o=o,
        filename=f"{funcname}-pressure.pdf",
        color="memorypressure",
        node_attr=node_attrs,
        cmap="Reds",
        collapse_outputs=False,
        **kwargs,
    )
    visualize_dsk(
        dsk,
        o=o,
        filename=f"{funcname}-cpath.pdf",
        color="cpath",
        node_attr=node_attrs,
        collapse_outputs=False,
        **kwargs,
    )


def _rechunk_merge_graph():
    return {
        ("transpose", 0, 0, 0): (f, ("concat-groupby", 0, 0, 0)),
        ("transpose", 0, 1, 0): (f, ("concat-groupby", 0, 1, 0)),
        ("transpose", 1, 0, 0): (f, ("concat-groupby", 1, 0, 0)),
        ("transpose", 1, 1, 0): (f, ("concat-groupby", 1, 1, 0)),
        ("groupby-cohort", 0, 0, 0): (f, ("groupby-chunk", 0, 0, 0)),
        ("groupby-cohort", 0, 0, 1): (f, ("groupby-chunk", 0, 0, 1)),
        ("groupby-cohort", 1, 0, 0): (f, ("groupby-chunk", 1, 0, 0)),
        ("groupby-cohort", 1, 0, 1): (f, ("groupby-chunk", 1, 0, 1)),
        ("groupby-cohort-2", 0, 0, 0): (f, ("groupby-chunk-2", 0, 0, 0)),
        ("groupby-cohort-2", 0, 0, 1): (f, ("groupby-chunk-2", 0, 0, 1)),
        ("groupby-cohort-2", 1, 0, 0): (f, ("groupby-chunk-2", 1, 0, 0)),
        ("groupby-cohort-2", 1, 0, 1): (f, ("groupby-chunk-2", 1, 0, 1)),
        ("rechunk-merge", 3, 0, 0): (
            f,
            ("concat-shuffle", 4, 0, 0),
            ("rechunk-split", 12),
        ),
        ("rechunk-merge", 0, 0, 0): (
            f,
            ("rechunk-split", 1),
            ("concat-shuffle", 0, 0, 0),
        ),
        ("rechunk-merge", 3, 1, 0): (
            f,
            ("rechunk-split", 14),
            ("concat-shuffle", 4, 1, 0),
        ),
        ("rechunk-merge", 2, 1, 0): (f, ("rechunk-split", 10), ("rechunk-split", 11)),
        ("rechunk-split", 12): (f, ("concat-shuffle", 3, 0, 0)),
        ("rechunk-merge", 0, 1, 0): (
            f,
            ("rechunk-split", 3),
            ("concat-shuffle", 0, 1, 0),
        ),
        ("rechunk-merge", 1, 0, 0): (f, ("rechunk-split", 4), ("rechunk-split", 5)),
        ("rechunk-merge", 1, 1, 0): (f, ("rechunk-split", 7), ("rechunk-split", 6)),
        ("rechunk-split", 5): (f, ("concat-shuffle", 2, 0, 0)),
        ("rechunk-split", 11): (f, ("concat-shuffle", 3, 1, 0)),
        ("rechunk-merge", 2, 0, 0): (f, ("rechunk-split", 8), ("rechunk-split", 9)),
        ("rechunk-split", 1): (f, ("concat-shuffle", 1, 0, 0)),
        ("rechunk-split", 14): (f, ("concat-shuffle", 3, 1, 0)),
        ("rechunk-split", 4): (f, ("concat-shuffle", 1, 0, 0)),
        ("rechunk-split", 7): (f, ("concat-shuffle", 2, 1, 0)),
        ("rechunk-split", 10): (f, ("concat-shuffle", 2, 1, 0)),
        ("rechunk-split", 6): (f, ("concat-shuffle", 1, 1, 0)),
        ("rechunk-split", 3): (f, ("concat-shuffle", 1, 1, 0)),
        ("rechunk-split", 9): (f, ("concat-shuffle", 3, 0, 0)),
        ("rechunk-split", 8): (f, ("concat-shuffle", 2, 0, 0)),
        ("concat-shuffle", 0, 0, 0): (f, ("shuffle-split", 0), ("shuffle-split", 1)),
        ("concat-shuffle", 0, 1, 0): (
            f,
            ("shuffle-split", 106),
            ("shuffle-split", 107),
        ),
        ("concat-shuffle", 1, 0, 0): (
            f,
            ("shuffle-split", 4665),
            ("shuffle-split", 4664),
        ),
        ("concat-shuffle", 1, 1, 0): (
            f,
            ("shuffle-split", 4770),
            ("shuffle-split", 4771),
        ),
        ("concat-shuffle", 2, 0, 0): (
            f,
            ("shuffle-split", 9328),
            ("shuffle-split", 9329),
            ("shuffle-split", 9330),
        ),
        ("concat-shuffle", 2, 1, 0): (
            f,
            ("shuffle-split", 9487),
            ("shuffle-split", 9488),
            ("shuffle-split", 9489),
        ),
        ("concat-shuffle", 3, 0, 0): (
            f,
            ("shuffle-split", 16324),
            ("shuffle-split", 16325),
        ),
        ("concat-shuffle", 3, 1, 0): (
            f,
            ("shuffle-split", 16430),
            ("shuffle-split", 16431),
        ),
        ("concat-shuffle", 4, 0, 0): (
            f,
            ("shuffle-split", 20989),
            ("shuffle-split", 20988),
        ),
        ("concat-shuffle", 4, 1, 0): (
            f,
            ("shuffle-split", 21094),
            ("shuffle-split", 21095),
        ),
        ("shuffle-split", 9487): (f, ("getitem-2", 2, 1, 0)),
        ("shuffle-split", 9489): (f, ("getitem-2", 14, 1, 0)),
        ("shuffle-split", 106): (f, ("getitem-open", 106)),
        ("shuffle-split", 4664): (f, ("getitem-2", 1, 0, 0)),
        ("shuffle-split", 16431): (f, ("getitem-2", 15, 1, 0)),
        ("shuffle-split", 16324): (f, ("getitem-2", 14, 0, 0)),
        ("shuffle-split", 107): (f, ("getitem-2", 1, 1, 0)),
        ("shuffle-split", 4665): (f, ("getitem-2", 2, 0, 0)),
        ("shuffle-split", 4770): (f, ("getitem-2", 1, 1, 0)),
        ("shuffle-split", 0): (f, ("getitem-open", 0)),
        ("shuffle-split", 9328): (f, ("getitem-2", 2, 0, 0)),
        ("shuffle-split", 9488): (f, ("getitem-open", 9488)),
        ("shuffle-split", 16325): (f, ("getitem-2", 15, 0, 0)),
        ("shuffle-split", 16430): (f, ("getitem-2", 14, 1, 0)),
        ("shuffle-split", 20988): (f, ("getitem-2", 15, 0, 0)),
        ("shuffle-split", 9329): (f, ("getitem-open", 9329)),
        ("shuffle-split", 4771): (f, ("getitem-2", 2, 1, 0)),
        ("shuffle-split", 1): (f, ("getitem-2", 1, 0, 0)),
        ("shuffle-split", 20989): (f, ("getitem-open", 20989)),
        ("shuffle-split", 9330): (f, ("getitem-2", 14, 0, 0)),
        ("shuffle-split", 21094): (f, ("getitem-2", 15, 1, 0)),
        ("shuffle-split", 21095): (f, ("getitem-open", 21095)),
        ("getitem-2", 1, 0, 0): (f, ("open_dataset", 1, 0, 0)),
        ("getitem-2", 14, 0, 0): (f, ("open_dataset", 14, 0, 0)),
        ("getitem-2", 2, 1, 0): (f, ("open_dataset", 2, 1, 0)),
        ("getitem-2", 15, 0, 0): (f, ("open_dataset", 15, 0, 0)),
        ("getitem-2", 15, 1, 0): (f, ("open_dataset", 15, 1, 0)),
        ("getitem-2", 2, 0, 0): (f, ("open_dataset", 2, 0, 0)),
        ("getitem-2", 1, 1, 0): (f, ("open_dataset", 1, 1, 0)),
        ("getitem-2", 14, 1, 0): (f, ("open_dataset", 14, 1, 0)),
        ("groupby-chunk-2", 0, 0, 1): (f, ("rechunk-merge", 2, 0, 0)),
        ("groupby-chunk-2", 0, 0, 0): (f, ("rechunk-merge", 0, 0, 0)),
        ("concat-groupby", 0, 0, 0): (
            f,
            ("groupby-cohort-2", 0, 0, 0),
            ("groupby-cohort-2", 0, 0, 1),
        ),
        ("groupby-chunk", 0, 0, 1): (f, ("rechunk-merge", 3, 0, 0)),
        ("groupby-chunk", 0, 0, 0): (f, ("rechunk-merge", 1, 0, 0)),
        ("concat-groupby", 1, 0, 0): (
            f,
            ("groupby-cohort", 0, 0, 0),
            ("groupby-cohort", 0, 0, 1),
        ),
        ("groupby-chunk", 1, 0, 1): (f, ("rechunk-merge", 3, 1, 0)),
        ("groupby-chunk", 1, 0, 0): (f, ("rechunk-merge", 1, 1, 0)),
        ("concat-groupby", 1, 1, 0): (
            f,
            ("groupby-cohort", 1, 0, 0),
            ("groupby-cohort", 1, 0, 1),
        ),
        ("open_dataset", 14, 1, 0): (f,),
        ("groupby-chunk-2", 1, 0, 0): (f, ("rechunk-merge", 0, 1, 0)),
        ("groupby-chunk-2", 1, 0, 1): (f, ("rechunk-merge", 2, 1, 0)),
        ("concat-groupby", 0, 1, 0): (
            f,
            ("groupby-cohort-2", 1, 0, 1),
            ("groupby-cohort-2", 1, 0, 0),
        ),
        ("getitem-open", 9329): (f,),
        ("open_dataset", 2, 1, 0): (f,),
        ("open_dataset", 15, 1, 0): (f),
        ("getitem-open", 20989): (f,),
        ("getitem-open", 0): (f,),
        ("open_dataset", 1, 0, 0): (f,),
        ("getitem-open", 9488): (f,),
        ("getitem-open", 21095): (f,),
        ("open_dataset", 2, 0, 0): (f,),
        ("getitem-open", 106): (f,),
        ("open_dataset", 1, 1, 0): (f,),
        ("open_dataset", 14, 0, 0): (f),
        ("open_dataset", 15, 0, 0): (f),
    }


@pytest.mark.skipif(WINDOWS, reason="Broken pipe on Windows")
def test_visualize_int_overflow():
    pytest.importorskip("graphviz")
    pytest.importorskip("matplotlib")
    # regression test for https://github.com/dask/dask/pull/11440
    visualize(_rechunk_merge_graph())


def test_ordering_keeps_groups_together(abcde):
    a, b, c, d, e = abcde
    d = {(a, i): (f,) for i in range(4)}
    d.update({(b, 0): (f, (a, 0), (a, 1)), (b, 1): (f, (a, 2), (a, 3))})
    o = order(d)
    assert_topological_sort(d, o)

    assert abs(o[(a, 0)] - o[(a, 1)]) == 1
    assert abs(o[(a, 2)] - o[(a, 3)]) == 1

    d = {(a, i): (f,) for i in range(4)}
    d.update({(b, 0): (f, (a, 0), (a, 2)), (b, 1): (f, (a, 1), (a, 3))})
    o = order(d)
    assert_topological_sort(d, o)

    assert abs(o[(a, 0)] - o[(a, 2)]) == 1
    assert abs(o[(a, 1)] - o[(a, 3)]) == 1


def test_avoid_broker_nodes(abcde):
    r"""

    b0    b1  b2
    |      \  /
    a0      a1

    There are good arguments for both a0 or a1 to run first. Regardless of what
    we choose to run first, we should finish the computation branch before
    moving to the other one
    """
    a, b, c, d, e = abcde
    dsk = {
        (a, 0): (f,),
        (a, 1): (f,),
        (b, 0): (f, (a, 0)),
        (b, 1): (f, (a, 1)),
        (b, 2): (f, (a, 1)),
    }
    o = order(dsk)
    assert_topological_sort(dsk, o)
    assert o[(a, 1)] < o[(b, 0)] or (o[(b, 1)] < o[(a, 0)] and o[(b, 2)] < o[(a, 0)])


@pytest.mark.parametrize("data_root", [True, False])
def test_base_of_reduce_preferred(abcde, data_root):
    r"""
               a3
              /|
            a2 |
           /|  |
         a1 |  |
        /|  |  |
      a0 |  |  |
      |  |  |  |
      b0 b1 b2 b3
        \ \ / /
           c

    We really want to run b0 quickly
    """
    a, b, c, d, e = abcde
    dsk = {(a, i): (f, (a, i - 1), (b, i)) for i in [1, 2, 3]}
    dsk[(a, 0)] = (f, (b, 0))
    dsk.update({(b, i): (f, c, 1) for i in [0, 1, 2, 3]})
    if data_root:
        dsk[c] = 1
    else:
        dsk[c] = (f, 1)

    o = order(dsk)
    assert_topological_sort(dsk, o)

    assert o[(b, 0)] <= 1
    assert o[(b, 1)] <= 3


def test_avoid_upwards_branching(abcde):
    r"""
       a1
       |
       a2
       |
       a3    d1
      /  \  /
    b1    c1
    |     |
    b2    c2
          |
          c3
    """
    a, b, c, d, e = abcde
    dsk = {
        (a, 1): (f, (a, 2)),
        (a, 2): (f, (a, 3)),
        (a, 3): (f, (b, 1), (c, 1)),
        (b, 1): (f, (b, 2)),
        (c, 1): (f, (c, 2)),
        (c, 2): (f, (c, 3)),
        (d, 1): (f, (c, 1)),
        (c, 3): 1,
        (b, 2): 1,
    }

    o = order(dsk)
    assert_topological_sort(dsk, o)

    assert o[(d, 1)] < o[(b, 1)]


def test_avoid_upwards_branching_complex(abcde):
    r"""
         a1
         |
    e2   a2  d2  d3
    |    |    \  /
    e1   a3    d1
     \  /  \  /
      b1    c1
      |     |
      b2    c2
            |
            c3

    Prefer c1 over b1 because c1 will stay in memory less long while b1
    computes
    """
    a, b, c, d, e = abcde
    dsk = {
        (a, 1): (f, (a, 2)),
        (a, 2): (f, (a, 3)),
        (a, 3): (f, (b, 1), (c, 1)),
        (b, 1): (f, (b, 2)),
        (b, 2): (f,),
        (c, 1): (f, (c, 2)),
        (c, 2): (f, (c, 3)),
        (c, 3): (f,),
        (d, 1): (f, (c, 1)),
        (d, 2): (f, (d, 1)),
        (d, 3): (f, (d, 1)),
        (e, 1): (f, (b, 1)),
        (e, 2): (f, (e, 1)),
    }

    o = order(dsk)
    assert_topological_sort(dsk, o)
    assert o[(c, 1)] < o[(b, 1)]
    assert abs(o[(d, 2)] - o[(d, 3)]) == 1


def test_deep_bases_win_over_dependents(abcde):
    r"""
    It's not clear who should run first, e or d

    1.  d is nicer because it exposes parallelism
    2.  e is nicer (hypothetically) because it will be sooner released
        (though in this case we need d to run first regardless)

    Regardless of e or d first, we should run b before c.

            a
          / | \   .
         b  c |
        / \ | /
       e    d
    """
    a, b, c, d, e = abcde
    dsk = {
        a: (f, b, c, d),
        b: (f, d, e),
        c: (f, d),
        d: (f, 1),
        e: (f, 2),
    }

    o = order(dsk)
    assert_topological_sort(dsk, o)
    assert o[b] < o[c]


def test_prefer_deep(abcde):
    """
        c
        |
    e   b
    |   |
    d   a

    Prefer longer chains first so we should start with d
    """
    a, b, c, d, e = abcde
    dsk = {a: 1, b: (f, a), c: (f, b), d: 1, e: (f, d)}

    o = order(dsk)
    assert_topological_sort(dsk, o)
    assert o[a] < o[d]
    assert o[b] < o[d]


def test_stacklimit(abcde):
    dsk = {"x%s" % (i + 1): (inc, "x%s" % i) for i in range(10000)}
    dependencies, dependents = get_deps(dsk)
    ndependencies(dependencies, dependents)


def test_break_ties_by_str(abcde):
    a, b, c, d, e = abcde
    dsk = {("x", i): (inc, i) for i in range(10)}
    x_keys = sorted(dsk)
    dsk["y"] = (f, list(x_keys))

    o = order(dsk)
    assert_topological_sort(dsk, o)
    expected = {"y": 10}
    expected.update({k: i for i, k in enumerate(x_keys[::-1])})

    assert o == expected


def test_order_doesnt_fail_on_mixed_type_keys(abcde):
    order({"x": (inc, 1), ("y", 0): (inc, 2), "z": (add, "x", ("y", 0))})


def test_gh_3055():
    pytest.importorskip("numpy")
    da = pytest.importorskip("dask.array")
    A, B = 20, 99
    orig = x = da.random.normal(size=(A, B), chunks=(1, None))
    for _ in range(2):
        y = (x[:, None, :] * x[:, :, None]).cumsum(axis=0)
        x = x.cumsum(axis=0)
    w = (y * x[:, None]).sum(axis=(1, 2))

    dsk = dict(w.__dask_graph__())
    o = order(dsk)
    assert_topological_sort(dsk, o)
    assert max(diagnostics(dsk, o=o)[1]) <= 8
    L = [o[k] for k in w.__dask_keys__()]
    assert sum(x < len(o) / 2 for x in L) > len(L) / 3  # some complete quickly

    L = [o[k] for kk in orig.__dask_keys__() for k in kk]
    assert sum(x > len(o) / 2 for x in L) > len(L) / 3  # some start later

    assert sorted(L) == L  # operate in order


def test_type_comparisions_ok(abcde):
    a, b, c, d, e = abcde
    dsk = {a: 1, (a, 1): 2, (a, b, 1): 3}
    order(dsk)  # this doesn't err


def test_favor_longest_critical_path(abcde):
    r"""

       a
       |
    d  b  e
     \ | /
       c

    """
    a, b, c, d, e = abcde
    dsk = {c: (f,), d: (f, c), e: (f, c), b: (f, c), a: (f, b)}

    o = order(dsk)
    assert_topological_sort(dsk, o)
    assert o[d] > o[b]
    assert o[e] > o[b]


def test_run_smaller_sections(abcde):
    r"""
            aa
           / |
      b   d  bb dd
     / \ /|  | /
    a   c e  cc

    """
    a, b, c, d, e = abcde
    aa, bb, cc, dd = (x * 2 for x in [a, b, c, d])

    dsk = {
        a: (f,),
        c: (f,),
        e: (f,),
        cc: (f,),
        b: (f, a, c),
        d: (f, c, e),
        bb: (f, cc),
        aa: (f, d, bb),
        dd: (f, cc),
    }
    o = order(dsk)
    assert_topological_sort(dsk, o)
    # The actual order is a bit ambiguous but the optimal pressure is well
    # defined, see also
    # https://github.com/dask/dask/pull/10535/files#r1337528255
    assert max(diagnostics(dsk)[1]) == 3


def test_local_parents_of_reduction(abcde):
    """

            c1
            |
        b1  c2
        |  /|
    a1  b2  c3
    |  /|
    a2  b3
    |
    a3

    Prefer to finish a1 stack before proceeding to b2
    """
    a, b, c, d, e = abcde
    a1, a2, a3 = (a + i for i in "123")
    b1, b2, b3 = (b + i for i in "123")
    c1, c2, c3 = (c + i for i in "123")

    expected = [a3, a2, a1, b3, b2, b1, c3, c2, c1]

    log = []

    def f(x):
        def _(*args):
            log.append(x)

        return _

    dsk = {
        a3: (f(a3),),
        a2: (f(a2), a3),
        a1: (f(a1), a2),
        b3: (f(b3),),
        b2: (f(b2), b3, a2),
        b1: (f(b1), b2),
        c3: (f(c3),),
        c2: (f(c2), c3, b2),
        c1: (f(c1), c2),
    }

    o = order(dsk)
    assert_topological_sort(dsk, o)
    dask.get(dsk, [a1, b1, c1])  # trigger computation

    assert log == expected


def test_nearest_neighbor(abcde):
    r"""

    a1  a2  a3  a4  a5  a6  a7 a8  a9
     \  |  /  \ |  /  \ |  / \ |  /
        b1      b2      b3     b4

    Want to finish off a local group before moving on.
    This is difficult because all groups are connected.
    """
    a, b, c, _, _ = abcde
    a1, a2, a3, a4, a5, a6, a7, a8, a9 = (a + i for i in "123456789")
    b1, b2, b3, b4 = (b + i for i in "1234")

    dsk = {
        b1: (f,),
        b2: (f,),
        b3: (f,),
        b4: (f,),
        a1: (f, b1),
        a2: (f, b1),
        a3: (f, b1, b2),
        a4: (f, b2),
        a5: (f, b2, b3),
        a6: (f, b3),
        a7: (f, b3, b4),
        a8: (f, b4),
        a9: (f, b4),
    }

    o = order(dsk)
    assert_topological_sort(dsk, o)

    assert 3 < sum(o[a + i] < len(o) / 2 for i in "123456789") < 7
    assert 1 < sum(o[b + i] < len(o) / 2 for i in "1234") < 4


def test_string_ordering():
    """Prefer ordering tasks by name first"""
    dsk = {("a", 1): (f,), ("a", 2): (f,), ("a", 3): (f,)}
    o = order(dsk)
    assert_topological_sort(dsk, o)
    assert o == {("a", 1): 0, ("a", 2): 1, ("a", 3): 2} or o == {
        ("a", 1): 2,
        ("a", 2): 1,
        ("a", 3): 0,
    }


@pytest.mark.xfail(reason="see comment", strict=False)
def test_string_ordering_dependents():
    """Prefer ordering tasks by name first even when in dependencies"""
    # XFAIL This is a little too artificial. While we can construct the
    # algorithm in a way that respects the key ordering strictly, there is not
    # necessarily a point to it and we can save ourselves one sorting step with
    # this.
    # See comment in add_to_result
    dsk = {("a", 1): (f, "b"), ("a", 2): (f, "b"), ("a", 3): (f, "b"), "b": (f,)}
    o = order(dsk)
    assert_topological_sort(dsk, o)
    assert o == {"b": 0, ("a", 1): 1, ("a", 2): 2, ("a", 3): 3} or o == {
        "b": 0,
        ("a", 1): 3,
        ("a", 2): 2,
        ("a", 3): 1,
    }


def test_prefer_short_narrow(abcde):
    # See test_prefer_short_ancestor for a fail case.
    a, b, c, _, _ = abcde
    dsk = {
        (a, 0): 0,
        (b, 0): 0,
        (c, 0): 0,
        (c, 1): (f, (c, 0), (a, 0), (b, 0)),
        (a, 1): 1,
        (b, 1): 1,
        (c, 2): (f, (c, 1), (a, 1), (b, 1)),
    }
    o = order(dsk)
    assert_topological_sort(dsk, o)
    assert o[(b, 0)] < o[(b, 1)]
    assert o[(b, 0)] < o[(c, 2)]
    assert o[(c, 1)] < o[(c, 2)]


def test_prefer_short_ancestor(abcde):
    r"""
    From https://github.com/dask/dask-ml/issues/206#issuecomment-395869929

    Two cases, one where chunks of an array are independent, and one where the
    chunks of an array have a shared source. We handled the independent one
    "well" earlier.

    Good:

                    c2
                   / \ \
                  /   \ \
                c1     \ \
              / | \     \ \
            c0  a0 b0   a1 b1

    Bad:

                    c2
                   / \ \
                  /   \ \
                c1     \ \
              / | \     \ \
            c0  a0 b0   a1 b1
                   \ \   / /
                    \ \ / /
                      a-b


    The difference is that all the `a` and `b` tasks now have a common
    ancestor.

    We would like to choose c1 *before* a1, and b1 because

    * we can release a0 and b0 once c1 is done
    * we don't need a1 and b1 to compute c1.
    """
    a, b, c, _, _ = abcde
    ab = a + b

    dsk = {
        ab: 0,
        (a, 0): (f, ab, 0, 0),
        (b, 0): (f, ab, 0, 1),
        (c, 0): 0,
        (c, 1): (f, (c, 0), (a, 0), (b, 0)),
        (a, 1): (f, ab, 1, 0),
        (b, 1): (f, ab, 1, 1),
        (c, 2): (f, (c, 1), (a, 1), (b, 1)),
    }
    o = order(dsk)
    assert_topological_sort(dsk, o)

    assert o[(a, 0)] < o[(a, 1)]
    assert o[(b, 0)] < o[(b, 1)]
    assert o[(b, 0)] < o[(c, 2)]
    assert o[(c, 1)] < o[(c, 2)]
    assert o[(c, 1)] < o[(a, 1)]


def test_map_overlap(abcde):
    r"""
      b1      b3      b5
       |\    / | \  / |
      c1  c2  c3  c4  c5
       |/  | \ | / | \|
      d1  d2  d3  d4  d5
       |       |      |
      e1      e3      e5

    Want to finish b1 before we start on e5
    """
    a, b, c, d, e = abcde
    dsk = {
        (e, 1): (f,),
        (d, 1): (f, (e, 1)),
        (c, 1): (f, (d, 1)),
        (b, 1): (f, (c, 1), (c, 2)),
        (d, 2): (f,),
        (c, 2): (f, (d, 1), (d, 2), (d, 3)),
        (e, 3): (f,),
        (d, 3): (f, (e, 3)),
        (c, 3): (f, (d, 3)),
        (b, 3): (f, (c, 2), (c, 3), (c, 4)),
        (d, 4): (f,),
        (c, 4): (f, (d, 3), (d, 4), (d, 5)),
        (e, 5): (f,),
        (d, 5): (f, (e, 5)),
        (c, 5): (f, (d, 5)),
        (b, 5): (f, (c, 4), (c, 5)),
    }

    o = order(dsk)
    assert_topological_sort(dsk, o)

    assert o[(b, 1)] < o[(e, 5)] or o[(b, 5)] < o[(e, 1)]


def test_use_structure_not_keys(abcde):
    """See https://github.com/dask/dask/issues/5584#issuecomment-554963958

    We were using key names to infer structure, which could result in funny behavior.
    """
    a, b, _, _, _ = abcde
    dsk = {
        (a, 0): (f,),
        (a, 1): (f,),
        (a, 2): (f,),
        (a, 3): (f,),
        (a, 4): (f,),
        (a, 5): (f,),
        (a, 6): (f,),
        (a, 7): (f,),
        (a, 8): (f,),
        (a, 9): (f,),
        (b, 5): (f, (a, 2)),
        (b, 7): (f, (a, 0), (a, 2)),
        (b, 9): (f, (a, 7), (a, 0), (a, 2)),
        (b, 1): (f, (a, 4), (a, 7), (a, 0)),
        (b, 2): (f, (a, 9), (a, 4), (a, 7)),
        (b, 4): (f, (a, 6), (a, 9), (a, 4)),
        (b, 3): (f, (a, 5), (a, 6), (a, 9)),
        (b, 8): (f, (a, 1), (a, 5), (a, 6)),
        (b, 6): (f, (a, 8), (a, 1), (a, 5)),
        (b, 0): (f, (a, 3), (a, 8), (a, 1)),
    }
    o = order(dsk)
    assert_topological_sort(dsk, o)

    assert max(diagnostics(dsk, o=o)[1]) == 3
    As = sorted(val for (letter, _), val in o.items() if letter == a)
    Bs = sorted(val for (letter, _), val in o.items() if letter == b)
    assert Bs[0] in {1, 3}
    if Bs[0] == 3:
        assert As == [0, 1, 2, 4, 6, 8, 10, 12, 14, 16]
        assert Bs == [3, 5, 7, 9, 11, 13, 15, 17, 18, 19]
    else:
        assert As == [0, 2, 4, 6, 8, 10, 12, 14, 16, 18]
        assert Bs == [1, 3, 5, 7, 9, 11, 13, 15, 17, 19]


def test_dont_run_all_dependents_too_early(abcde):
    """From https://github.com/dask/dask-ml/issues/206#issuecomment-395873372"""
    a, b, c, d, e = abcde
    depth = 10
    dsk = {
        (a, 0): (f, 0),
        (b, 0): (f, 1),
        (c, 0): (f, 2),
        (d, 0): (f, (a, 0), (b, 0), (c, 0)),
    }
    for i in range(1, depth):
        dsk[(b, i)] = (f, (b, 0))
        dsk[(c, i)] = (f, (c, 0))
        dsk[(d, i)] = (f, (d, i - 1), (b, i), (c, i))
    o = order(dsk)
    assert_topological_sort(dsk, o)

    expected = [3, 6, 9, 12, 15, 18, 21, 24, 27, 30]
    actual = sorted(v for (letter, num), v in o.items() if letter == d)
    assert expected == actual


def test_many_branches_use_ndependencies(abcde):
    """From https://github.com/dask/dask/pull/5646#issuecomment-562700533

    Sometimes we need larger or wider DAGs to test behavior.  This test
    ensures we choose the branch with more work twice in succession.
    This is important, because ``order`` may search along dependencies
    and then along dependents.

    """
    a, b, c, d, e = abcde
    dd = d + d
    ee = e + e
    dsk = {
        (a, 0): 0,
        (a, 1): (f, (a, 0)),
        (a, 2): (f, (a, 1)),
        (b, 1): (f, (a, 0)),
        (b, 2): (f, (b, 1)),
        (c, 1): (f, (a, 0)),  # most short and thin; should go last
        (d, 1): (f, (a, 0)),
        (d, 2): (f, (d, 1)),
        (dd, 1): (f, (a, 0)),
        (dd, 2): (f, (dd, 1)),
        (dd, 3): (f, (d, 2), (dd, 2)),
        (e, 1): (f, (a, 0)),
        (e, 2): (f, (e, 1)),
        (ee, 1): (f, (a, 0)),
        (ee, 2): (f, (ee, 1)),
        (ee, 3): (f, (e, 2), (ee, 2)),
        (a, 3): (f, (a, 2), (b, 2), (c, 1), (dd, 3), (ee, 3)),
    }
    o = order(dsk)
    assert_topological_sort(dsk, o)
    # run all d's and e's first
    expected = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    actual = sorted(v for (letter, _), v in o.items() if letter in {d, dd, e, ee})
    assert actual == expected
    assert o[(c, 1)] == o[(a, 3)] - 1


def test_order_cycle():
    with pytest.raises(RuntimeError, match="Cycle detected"):
        dask.get({"a": (f, "a")}, "a")  # we encounter this in `get`
    with pytest.raises(RuntimeError, match="Cycle detected"):
        order({"a": (f, "a")})  # trivial self-loop
    with pytest.raises(RuntimeError, match="Cycle detected"):
        order({("a", 0): (f, ("a", 0))})  # non-string
    with pytest.raises(RuntimeError, match="Cycle detected"):
        order({"a": (f, "b"), "b": (f, "c"), "c": (f, "a")})  # non-trivial loop
    with pytest.raises(RuntimeError, match="Cycle detected"):
        order({"a": (f, "b"), "b": (f, "c"), "c": (f, "a", "d"), "d": 1})
    with pytest.raises(RuntimeError, match="Cycle detected"):
        order({"a": (f, "b"), "b": (f, "c"), "c": (f, "a", "d"), "d": (f, "b")})


def test_order_empty():
    assert order({}) == {}


def test_order_with_equal_dependents(abcde):
    """From https://github.com/dask/dask/issues/5859#issuecomment-608422198

    See the visualization of `(maxima, argmax)` example from the above comment.

    This DAG has enough structure to exercise more parts of `order`

    """
    # Lower pressure is better but this is where we are right now. Important is
    # that no variation below should be worse since all variations below should
    # reduce to the same graph when optimized/fused.
    max_pressure = 10
    a, b, c, d, e = abcde
    dsk = {}
    abc = [a, b, c, d]
    for x in abc:
        dsk.update(
            {
                (x, 0): (f, 0),
                (x, 1): (f, (x, 0)),
                (x, 2, 0): (f, (x, 0)),
                (x, 2, 1): (f, (x, 1)),
            }
        )
        for i, y in enumerate(abc):
            dsk.update(
                {
                    (x, 3, i): (f, (x, 2, 0), (y, 2, 1)),  # cross x and y
                    (x, 4, i): (f, (x, 3, i)),
                    (x, 5, i, 0): (f, (x, 4, i)),
                    (x, 5, i, 1): (f, (x, 4, i)),
                    (x, 6, i, 0): (f, (x, 5, i, 0)),
                    (x, 6, i, 1): (f, (x, 5, i, 1)),
                }
            )
    o = order(dsk)
    assert_topological_sort(dsk, o)
    total = 0
    for x in abc:
        for i in range(len(abc)):
            val = abs(o[(x, 6, i, 1)] - o[(x, 6, i, 0)])
            total += val

    assert total <= 32  # ideally, this should be 2 * 16 == 32
    pressure = diagnostics(dsk, o=o)[1]
    assert max(pressure) <= max_pressure
    # Add one to the end of the nine bundles
    dsk2 = dict(dsk)
    for x in abc:
        for i in range(len(abc)):
            dsk2[(x, 7, i, 0)] = (f, (x, 6, i, 0))
    o = order(dsk2)
    assert_topological_sort(dsk2, o)
    total = 0
    for x in abc:
        for i in range(len(abc)):
            val = abs(o[(x, 6, i, 1)] - o[(x, 7, i, 0)])
            total += val
    assert total <= 48  # ideally, this should be 3 * 16 == 48
    pressure = diagnostics(dsk2, o=o)[1]
    assert max(pressure) <= max_pressure

    # Remove one from each of the nine bundles
    dsk3 = dict(dsk)
    for x in abc:
        for i in range(len(abc)):
            del dsk3[(x, 6, i, 1)]
    o = order(dsk3)
    assert_topological_sort(dsk3, o)
    total = 0
    for x in abc:
        for i in range(len(abc)):
            val = abs(o[(x, 5, i, 1)] - o[(x, 6, i, 0)])
            total += val
    assert total <= 32  # ideally, this should be 2 * 16 == 32
    pressure = diagnostics(dsk3, o=o)[1]
    assert max(pressure) <= max_pressure

    # # Remove another one from each of the nine bundles
    dsk4 = dict(dsk3)
    for x in abc:
        for i in range(len(abc)):
            del dsk4[(x, 6, i, 0)]
    o = order(dsk4)
    assert_topological_sort(dsk4, o)
    pressure = diagnostics(dsk4, o=o)[1]
    assert max(pressure) <= max_pressure
    for x in abc:
        for i in range(len(abc)):
            assert abs(o[(x, 5, i, 1)] - o[(x, 5, i, 0)]) <= 2


def test_terminal_node_backtrack():
    r"""
    https://github.com/dask/dask/issues/6745

    We have

    1. A terminal node that depends on the entire graph ('s')
    2. Some shared dependencies near the roots ('a1', 'a4')
    3. But the left and right halves are disconnected, other
       than the terminal node.

                       s
               /   /       \   \
              /   /         \   \
            s00  s10       s01  s11
             |    |         |    |
            b00  b10       b01  b11
            / \  / \       / \ / \
           a0  a1  a2    a3  a4  a5

    Previously we started at 'a', and worked up to 's00'. We'd like to finish
    's00' completely, so we progress to 's' and work through its dependencies.

    Ideally, we would choose 's10', since we've already computed one of its
    (eventual) dependencies: 'a1'. However, all of 's00' through 's11' had
    equal metrics so we fell back to the name tie-breaker and started on
    's01' (via 'a3', a4', 'b01', ...).
    """
    dsk = {
        # left half
        ("a", 0): (0,),
        ("a", 1): (1,),
        ("a", 2): (2,),
        ("b", 0): (f, ("a", 0), ("a", 1)),
        ("b", 1): (f, ("a", 1), ("a", 2)),
        ("store", 0, 0): ("b", 0),
        ("store", 1, 0): ("b", 1),
        # right half
        ("a", 3): (3,),
        ("a", 4): (4,),
        ("a", 5): (5,),
        ("b", 2): (f, ("a", 3), ("a", 4)),
        ("b", 3): (f, ("a", 4), ("a", 5)),
        ("store", 0, 1): ("b", 2),
        ("store", 1, 1): ("b", 3),
        "store": (
            f,
            ("store", 0, 0),
            ("store", 1, 0),
            ("store", 0, 1),
            ("store", 1, 1),
        ),
    }
    o = order(dsk)
    assert_topological_sort(dsk, o)
    assert o[("a", 2)] < o[("a", 3)]


def test_array_store_final_order(tmpdir):
    # https://github.com/dask/dask/issues/6745
    # This essentially tests the same thing as test_terminal_node_backtrack,
    # but with the graph actually generated by da.store.
    pytest.importorskip("numpy")
    da = pytest.importorskip("dask.array")
    zarr = pytest.importorskip("zarr")

    arrays = [da.ones((110, 4), chunks=(100, 2)) for i in range(4)]
    x = da.concatenate(arrays, axis=0).rechunk((100, 2))

    if Version(zarr.__version__) < Version("3.0.0.a0"):
        store = zarr.storage.DirectoryStore(tmpdir)
    else:
        store = zarr.storage.LocalStore(str(tmpdir), read_only=False)

    root = zarr.group(store, overwrite=True)
    dest = root.empty_like(name="dest", data=x, chunks=x.chunksize, overwrite=True)
    d = x.store(dest, lock=False, compute=False)
    o = order(d.dask)
    assert_topological_sort(dict(d.dask), o)
    # Find the lowest store. Dask starts here.
    stores = [k for k in o if isinstance(k, tuple) and k[0].startswith("store-map-")]
    first_store = min(stores, key=lambda k: o[k])
    connected_stores = [k for k in stores if k[-1] == first_store[-1]]
    disconnected_stores = [k for k in stores if k[-1] != first_store[-1]]

    connected_max = max(v for k, v in o.items() if k in connected_stores)
    disconnected_min = min(v for k, v in o.items() if k in disconnected_stores)
    assert connected_max < disconnected_min


def test_eager_to_compute_dependent_to_free_parent():
    r"""https://github.com/dask/dask/pull/7929

    This graph begins with many motifs like the following:

    |      |
    c1    c2
      \ /
       b
       |
       a

    We want to compute c2 and c3 pretty close together, because if we choose to
    compute c1, then we should also compute c2 so we can release b.  Being
    greedy here allows us to release memory sooner and be more globally optimal.
    """
    dsk = {
        "a00": (f, "a06", "a08"),
        "a01": (f, "a28", "a26"),
        "a02": (f, "a24", "a21"),
        "a03": (f, "a22", "a25"),
        "a04": (f, "a29", "a20"),
        "a05": (f, "a23", "a27"),
        "a06": (f, "a04", "a02"),
        "a07": (f, "a00", "a01"),
        "a08": (f, "a05", "a03"),
        "a09": (f, "a43"),
        "a10": (f, "a36"),
        "a11": (f, "a33"),
        "a12": (f, "a47"),
        "a13": (f, "a44"),
        "a14": (f, "a42"),
        "a15": (f, "a37"),
        "a16": (f, "a48"),
        "a17": (f, "a49"),
        "a18": (f, "a35"),
        "a19": (f, "a46"),
        "a20": (f, "a55"),
        "a21": (f, "a53"),
        "a22": (f, "a60"),
        "a23": (f, "a54"),
        "a24": (f, "a59"),
        "a25": (f, "a56"),
        "a26": (f, "a61"),
        "a27": (f, "a52"),
        "a28": (f, "a57"),
        "a29": (f, "a58"),
        "a30": (f, "a19"),
        "a31": (f, "a07"),
        "a32": (f, "a30", "a31"),
        "a33": (f, "a58"),
        "a34": (f, "a11", "a09"),
        "a35": (f, "a60"),
        "a36": (f, "a52"),
        "a37": (f, "a61"),
        "a38": (f, "a14", "a10"),
        "a39": (f, "a38", "a40"),
        "a40": (f, "a18", "a17"),
        "a41": (f, "a34", "a50"),
        "a42": (f, "a54"),
        "a43": (f, "a55"),
        "a44": (f, "a53"),
        "a45": (f, "a16", "a15"),
        "a46": (f, "a51", "a45"),
        "a47": (f, "a59"),
        "a48": (f, "a57"),
        "a49": (f, "a56"),
        "a50": (f, "a12", "a13"),
        "a51": (f, "a41", "a39"),
        "a52": (f, "a62"),
        "a53": (f, "a68"),
        "a54": (f, "a70"),
        "a55": (f, "a67"),
        "a56": (f, "a71"),
        "a57": (f, "a64"),
        "a58": (f, "a65"),
        "a59": (f, "a63"),
        "a60": (f, "a69"),
        "a61": (f, "a66"),
        "a62": (f, f),
        "a63": (f, f),
        "a64": (f, f),
        "a65": (f, f),
        "a66": (f, f),
        "a67": (f, f),
        "a68": (f, f),
        "a69": (f, f),
        "a70": (f, f),
        "a71": (f, f),
    }
    o = order(dsk)
    assert_topological_sort(dsk, o)

    _, pressure = diagnostics(dsk, o=o)
    assert max(pressure) <= 8
    # Visualizing the graph shows that there are two deep, thick branches and
    # two shallow ones. We prefer to run the deep ones first and reduce them as
    # far as possible before starting another computation branch.
    # We don't care about the ordering so deep_roots1 and deep_roots2 are
    # ambiguous. At time of writing, deep_roots1 was executed first.
    shallow_roots = ["a64", "a66"]
    deep_roots1 = ["a68", "a63", "a65", "a67"]
    # These are the final reducers of the first, thick branch
    reducers_1 = ["a06", "a41"]
    deep_roots2 = ["a69", "a70", "a71", "a62"]
    reducers_2 = ["a39", "a00"]
    for deep_roots in [deep_roots1, deep_roots2]:
        assert max(o[r] for r in deep_roots) < min(o[r] for r in shallow_roots)
    # These two are the reduction nodes furthest along before additional roots
    # have to be loaded
    for reducer in [reducers_1, reducers_2]:
        for red in reducer:
            assert o[red] < min(o[r] for r in shallow_roots)


def test_diagnostics(abcde):
    r"""
        a1  b1  c1  d1  e1
        /|\ /|\ /|\ /|  /
       / | X | X | X | /
      /  |/ \|/ \|/ \|/
    a0  b0  c0  d0  e0
    """
    a, b, c, d, e = abcde
    dsk = {
        (a, 0): (f,),
        (b, 0): (f,),
        (c, 0): (f,),
        (d, 0): (f,),
        (e, 0): (f,),
        (a, 1): (f, (a, 0), (b, 0), (c, 0)),
        (b, 1): (f, (b, 0), (c, 0), (d, 0)),
        (c, 1): (f, (c, 0), (d, 0), (e, 0)),
        (d, 1): (f, (d, 0), (e, 0)),
        (e, 1): (f, (e, 0)),
    }
    o = order(dsk)
    assert_topological_sort(dsk, o)
    info, memory_over_time = diagnostics(dsk)
    # this is ambiguous, depending on whether we start from left or right
    assert all(o[key] == val.order for key, val in info.items())
    if o[(e, 1)] == 1:
        assert o[(e, 1)] == 1
        assert o[(d, 1)] == 3
        assert o[(c, 1)] == 5
        assert memory_over_time == [0, 1, 1, 2, 2, 3, 2, 3, 2, 3]
        assert {key: val.age for key, val in info.items()} == {
            (a, 0): 1,
            (b, 0): 3,
            (c, 0): 5,
            (d, 0): 5,
            (e, 0): 5,
            (a, 1): 0,
            (b, 1): 0,
            (c, 1): 0,
            (d, 1): 0,
            (e, 1): 0,
        }
        assert {key: val.num_dependencies_freed for key, val in info.items()} == {
            (a, 0): 0,
            (b, 0): 0,
            (c, 0): 0,
            (d, 0): 0,
            (e, 0): 0,
            (a, 1): 3,
            (b, 1): 1,
            (c, 1): 1,
            (d, 1): 0,
            (e, 1): 0,
        }
        assert {key: val.num_data_when_run for key, val in info.items()} == {
            (a, 0): 2,
            (b, 0): 2,
            (c, 0): 2,
            (d, 0): 1,
            (e, 0): 0,
            (a, 1): 3,
            (b, 1): 3,
            (c, 1): 3,
            (d, 1): 2,
            (e, 1): 1,
        }
        assert {key: val.num_data_when_released for key, val in info.items()} == {
            (a, 0): 3,
            (b, 0): 3,
            (c, 0): 3,
            (d, 0): 3,
            (e, 0): 3,
            (a, 1): 3,
            (b, 1): 3,
            (c, 1): 3,
            (d, 1): 2,
            (e, 1): 1,
        }
    else:
        assert o[(e, 1)] == len(dsk) - 1
        assert o[(d, 1)] == len(dsk) - 2
        assert o[(c, 1)] == len(dsk) - 3
        assert memory_over_time == [0, 1, 2, 3, 2, 3, 2, 3, 2, 1]
        assert {key: val.age for key, val in info.items()} == {
            (a, 0): 3,
            (b, 0): 4,
            (c, 0): 5,
            (d, 0): 4,
            (e, 0): 3,
            (a, 1): 0,
            (b, 1): 0,
            (c, 1): 0,
            (d, 1): 0,
            (e, 1): 0,
        }
        assert {key: val.num_dependencies_freed for key, val in info.items()} == {
            (a, 0): 0,
            (b, 0): 0,
            (c, 0): 0,
            (d, 0): 0,
            (e, 0): 0,
            (a, 1): 1,
            (b, 1): 1,
            (c, 1): 1,
            (d, 1): 1,
            (e, 1): 1,
        }
        assert {key: val.num_data_when_run for key, val in info.items()} == {
            (a, 0): 0,
            (b, 0): 1,
            (c, 0): 2,
            (d, 0): 2,
            (e, 0): 2,
            (a, 1): 3,
            (b, 1): 3,
            (c, 1): 3,
            (d, 1): 2,
            (e, 1): 1,
        }
        assert {key: val.num_data_when_released for key, val in info.items()} == {
            (a, 0): 3,
            (b, 0): 3,
            (c, 0): 3,
            (d, 0): 2,
            (e, 0): 1,
            (a, 1): 3,
            (b, 1): 3,
            (c, 1): 3,
            (d, 1): 2,
            (e, 1): 1,
        }


def test_xarray_like_reduction():
    a, b, c, d, e = list("abcde")

    dsk = {}
    for ix in range(3):
        part = {
            # Part1
            (a, 0, ix): (f,),
            (a, 1, ix): (f,),
            (b, 0, ix): (f, (a, 0, ix)),
            (b, 1, ix): (f, (a, 0, ix), (a, 1, ix)),
            (b, 2, ix): (f, (a, 1, ix)),
            (c, 0, ix): (f, (b, 0, ix)),
            (c, 1, ix): (f, (b, 1, ix)),
            (c, 2, ix): (f, (b, 2, ix)),
        }
        dsk.update(part)
    for ix in range(3):
        dsk.update(
            {
                (d, ix): (f, (c, ix, 0), (c, ix, 1), (c, ix, 2)),
            }
        )
    o = order(dsk)
    assert_topological_sort(dsk, o)
    _, pressure = diagnostics(dsk, o=o)
    assert max(pressure) <= 9


@pytest.mark.parametrize(
    "optimize",
    [
        True,
        False,
    ],
)
def test_array_vs_dataframe(optimize):
    xr = pytest.importorskip("xarray")
    pytest.importorskip("dask.dataframe")

    import dask.array as da

    size = 5000
    ds = xr.Dataset(
        dict(
            anom_u=(
                ["time", "face", "j", "i"],
                da.random.random((size, 1, 987, 1920), chunks=(10, 1, -1, -1)),
            ),
            anom_v=(
                ["time", "face", "j", "i"],
                da.random.random((size, 1, 987, 1920), chunks=(10, 1, -1, -1)),
            ),
        )
    )

    quad = ds**2
    quad["uv"] = ds.anom_u * ds.anom_v
    mean = quad.mean("time")
    diag_array = diagnostics(
        collections_to_expr([mean], optimize_graph=optimize).optimize().__dask_graph__()
    )
    diag_df = diagnostics(
        collections_to_expr([mean.to_dask_dataframe()], optimize_graph=optimize)
        .optimize()
        .__dask_graph__()
    )
    assert max(diag_df[1]) == 15
    assert max(diag_array[1]) == 38
    assert max(diag_array[1]) < 50


def test_anom_mean():
    np = pytest.importorskip("numpy")
    xr = pytest.importorskip("xarray")

    import dask.array as da

    data = da.random.random((200, 1), chunks=(1, -1))

    ngroups = 5
    arr = xr.DataArray(
        data,
        dims=["time", "x"],
        coords={"day": ("time", np.arange(data.shape[0]) % ngroups)},
    )

    clim = arr.groupby("day").mean(dim="time")
    anom = arr.groupby("day") - clim
    anom_mean = anom.mean(dim="time")
    graph = collections_to_expr([anom_mean]).optimize().__dask_graph__()
    _, dependents = get_deps(graph)
    diags, pressure = diagnostics(graph)
    # Encoding the "best" ordering for this graph is tricky. When inspecting the
    # visualization, one sees that there are small, connected tree-like steps at
    # the beginning (tested well below in anom_mean_raw) followed by a
    # concat+transpose per group (see ngroups above). This transpose task fans
    # out into many (20-30) getitem tasks that are tiny and feed into a
    # `mean_chunk` which is the primary reducer in this graph. Therefore we want
    # to run those as quickly as possible.
    # This is difficult to assert on but the pressure is an ok-ish proxy
    assert max(pressure) <= 178
    from collections import defaultdict

    count_dependents = defaultdict(set)
    for k in dict(graph).keys():
        count_dependents[len(dependents[k])].add(k)

    # array-taker has the most dependents, but it's not what we want to look at
    n_splits = sorted(count_dependents)[-2]
    # There is a transpose/stack group that is splitting into many tasks
    # see https://github.com/dask/dask/pull/10660#discussion_r1420571664
    # the name depends on the version of xarray
    assert n_splits > 30  # at time of writing 40
    transpose_tasks = count_dependents[n_splits]
    transpose_metrics = {k: diags[k] for k in transpose_tasks}
    assert len(transpose_metrics) == ngroups, {key_split(k) for k in diags}
    # This is a pretty tightly connected graph overall and we'll have to hold
    # many tasks in memory until this can complete. However, we should ensure
    # that we get to the mean_chunks asap while the transposes are released
    # quickly.
    # If this breaks, I suggest to visually inspect the graph and run the above
    # on a single threaded LocalCluster and verify that the progress is indeed
    # in five steps (i.e. five groups)
    ages_mean_chunks = {k: v.age for k, v in diags.items() if "mean_chunk" in k[0]}
    avg_age_mean_chunks = sum(ages_mean_chunks.values()) / len(ages_mean_chunks)
    max_age_mean_chunks = max(ages_mean_chunks.values())
    ages_tranpose = {k: v.age for k, v in transpose_metrics.items()}
    assert max_age_mean_chunks > 900
    assert avg_age_mean_chunks > 100
    avg_age_transpose = sum(ages_tranpose.values()) / len(ages_tranpose)
    max_age_transpose = max(ages_tranpose.values())
    assert max_age_transpose < 150
    assert avg_age_transpose < 100
    assert sum(pressure) / len(pressure) < 101


def test_anom_mean_raw(abcde):
    a, b, c, d, e = abcde
    g, h = "gh"
    dsk = {
        (d, 0, 0): (f, (a, 0, 0), (b, 1, 0, 0)),
        (d, 1, 0): (f, (a, 1, 0), (b, 1, 1, 0)),
        (d, 2, 0): (f, (a, 2, 0), (b, 1, 2, 0)),
        (d, 3, 0): (f, (a, 3, 0), (b, 1, 3, 0)),
        (d, 4, 0): (f, (a, 4, 0), (b, 1, 4, 0)),
        (a, 0, 0): (f, f, "random_sample", None, (1, 1), [], {}),
        (a, 1, 0): (f, f, "random_sample", None, (1, 1), [], {}),
        (a, 2, 0): (f, f, "random_sample", None, (1, 1), [], {}),
        (a, 3, 0): (f, f, "random_sample", None, (1, 1), [], {}),
        (a, 4, 0): (f, f, "random_sample", None, (1, 1), [], {}),
        (e, 0, 0): (f, (g, 1, 0)),
        (e, 1, 0): (f, (g, 3, 0)),
        (b, 0, 0, 0): (f, (a, 0, 0)),
        (b, 0, 1, 0): (f, (a, 2, 0)),
        (b, 0, 2, 0): (f, (a, 4, 0)),
        (c, 0, 0, 0): (f, (b, 0, 0, 0)),
        (c, 0, 1, 0): (f, (b, 0, 1, 0)),
        (c, 0, 2, 0): (f, (b, 0, 2, 0)),
        (g, 1, 0): (f, [(c, 0, 0, 0), (c, 0, 1, 0), (c, 0, 2, 0)]),
        (b, 2, 0, 0): (f, (a, 1, 0)),
        (b, 2, 1, 0): (f, (a, 3, 0)),
        (c, 1, 0, 0): (f, (b, 2, 0, 0)),
        (c, 1, 1, 0): (f, (b, 2, 1, 0)),
        (g, 3, 0): (f, [(c, 1, 0, 0), (c, 1, 1, 0)]),
        (b, 1, 0, 0): (f, (e, 0, 0)),
        (b, 1, 1, 0): (f, (e, 1, 0)),
        (b, 1, 2, 0): (f, (e, 0, 0)),
        (b, 1, 3, 0): (f, (e, 1, 0)),
        (b, 1, 4, 0): (f, (e, 0, 0)),
        (c, 2, 0, 0): (f, (d, 0, 0)),
        (c, 2, 1, 0): (f, (d, 1, 0)),
        (c, 2, 2, 0): (f, (d, 2, 0)),
        (c, 2, 3, 0): (f, (d, 3, 0)),
        (c, 2, 4, 0): (f, (d, 4, 0)),
        (h, 0, 0): (f, [(c, 2, 0, 0), (c, 2, 1, 0), (c, 2, 2, 0), (c, 2, 3, 0)]),
        (h, 1, 0): (f, [(c, 2, 4, 0)]),
        (g, 2, 0): (f, [(h, 0, 0), (h, 1, 0)]),
    }

    o = order(dsk)
    assert_topological_sort(dsk, o)
    # The left hand computation branch should complete before we start loading
    # more data
    nodes_to_finish_before_loading_more_data = [
        (h, 1, 0),
        (d, 0, 0),
        (d, 2, 0),
        (d, 4, 0),
    ]
    for n in nodes_to_finish_before_loading_more_data:
        assert o[n] < o[(a, 1, 0)]
        assert o[n] < o[(a, 3, 0)]


def test_flaky_array_reduction():
    first = {
        ("mean_agg-aggregate-10d721567ef5a0d6a0e1afae8a87c066", 0, 0, 0): (
            f,
            [
                ("mean_combine-partial-17c7b5c6eed42e203858b3f6dde16003", 0, 0, 0, 0),
                ("mean_combine-partial-17c7b5c6eed42e203858b3f6dde16003", 1, 0, 0, 0),
            ],
        ),
        ("mean_combine-partial-17c7b5c6eed42e203858b3f6dde16003", 0, 0, 0, 0): (
            f,
            [
                ("mean_chunk-98a32cd9f4fadbed908fffb32e0c9679", 0, 0, 0, 0),
                ("mean_chunk-98a32cd9f4fadbed908fffb32e0c9679", 1, 0, 0, 0),
                ("mean_chunk-98a32cd9f4fadbed908fffb32e0c9679", 2, 0, 0, 0),
                ("mean_chunk-98a32cd9f4fadbed908fffb32e0c9679", 3, 0, 0, 0),
            ],
        ),
        ("mean_combine-partial-17c7b5c6eed42e203858b3f6dde16003", 1, 0, 0, 0): (
            "mean_chunk-mean_combine-partial-17c7b5c6eed42e203858b3f6dde16003",
            1,
            0,
            0,
            0,
        ),
        ("mean_chunk-98a32cd9f4fadbed908fffb32e0c9679", 0, 0, 0, 0): (
            f,
            ("random_sample-e16bcfb15a013023c98a21e2f03d66a9", 0, 0, 0, 0),
            ("random_sample-02eaa4a8dbb23fac4db22ad034c401b3", 0, 0, 0, 0),
        ),
        ("mean_chunk-98a32cd9f4fadbed908fffb32e0c9679", 1, 0, 0, 0): (
            f,
            ("random_sample-02eaa4a8dbb23fac4db22ad034c401b3", 1, 0, 0, 0),
            ("random_sample-e16bcfb15a013023c98a21e2f03d66a9", 1, 0, 0, 0),
        ),
        ("mean_chunk-98a32cd9f4fadbed908fffb32e0c9679", 2, 0, 0, 0): (
            f,
            ("random_sample-e16bcfb15a013023c98a21e2f03d66a9", 2, 0, 0, 0),
            ("random_sample-02eaa4a8dbb23fac4db22ad034c401b3", 2, 0, 0, 0),
        ),
        ("mean_chunk-98a32cd9f4fadbed908fffb32e0c9679", 3, 0, 0, 0): (
            f,
            ("random_sample-e16bcfb15a013023c98a21e2f03d66a9", 3, 0, 0, 0),
            ("random_sample-02eaa4a8dbb23fac4db22ad034c401b3", 3, 0, 0, 0),
        ),
        ("mean_agg-aggregate-fdb340546b01334890192fcfa55fa0d9", 0, 0, 0): (
            f,
            [
                ("mean_combine-partial-23adb4747560e6e33afd63c5bb179709", 0, 0, 0, 0),
                ("mean_combine-partial-23adb4747560e6e33afd63c5bb179709", 1, 0, 0, 0),
            ],
        ),
        ("mean_combine-partial-23adb4747560e6e33afd63c5bb179709", 0, 0, 0, 0): (
            f,
            [
                ("mean_chunk-7edba1c5a284fcec88b9efdda6c2135f", 0, 0, 0, 0),
                ("mean_chunk-7edba1c5a284fcec88b9efdda6c2135f", 1, 0, 0, 0),
                ("mean_chunk-7edba1c5a284fcec88b9efdda6c2135f", 2, 0, 0, 0),
                ("mean_chunk-7edba1c5a284fcec88b9efdda6c2135f", 3, 0, 0, 0),
            ],
        ),
        ("mean_combine-partial-23adb4747560e6e33afd63c5bb179709", 1, 0, 0, 0): (
            "mean_chunk-mean_combine-partial-23adb4747560e6e33afd63c5bb179709",
            1,
            0,
            0,
            0,
        ),
        ("mean_chunk-7edba1c5a284fcec88b9efdda6c2135f", 0, 0, 0, 0): (
            f,
            ("random_sample-02eaa4a8dbb23fac4db22ad034c401b3", 0, 0, 0, 0),
            2,
        ),
        ("mean_chunk-7edba1c5a284fcec88b9efdda6c2135f", 1, 0, 0, 0): (
            f,
            ("random_sample-02eaa4a8dbb23fac4db22ad034c401b3", 1, 0, 0, 0),
            2,
        ),
        ("mean_chunk-7edba1c5a284fcec88b9efdda6c2135f", 2, 0, 0, 0): (
            f,
            ("random_sample-02eaa4a8dbb23fac4db22ad034c401b3", 2, 0, 0, 0),
            2,
        ),
        ("mean_chunk-7edba1c5a284fcec88b9efdda6c2135f", 3, 0, 0, 0): (
            f,
            ("random_sample-02eaa4a8dbb23fac4db22ad034c401b3", 3, 0, 0, 0),
            2,
        ),
        ("random_sample-02eaa4a8dbb23fac4db22ad034c401b3", 0, 0, 0, 0): (f, 1),
        ("random_sample-02eaa4a8dbb23fac4db22ad034c401b3", 1, 0, 0, 0): (f, 1),
        ("random_sample-02eaa4a8dbb23fac4db22ad034c401b3", 2, 0, 0, 0): (f, 1),
        ("random_sample-02eaa4a8dbb23fac4db22ad034c401b3", 3, 0, 0, 0): (f, 1),
        ("random_sample-02eaa4a8dbb23fac4db22ad034c401b3", 4, 0, 0, 0): (f, 1),
        ("mean_agg-aggregate-cc19342c8116d616fc6573f5d20b5762", 0, 0, 0): (
            f,
            [
                ("mean_combine-partial-0c98c5a4517f58f8268985e7464daace", 0, 0, 0, 0),
                ("mean_combine-partial-0c98c5a4517f58f8268985e7464daace", 1, 0, 0, 0),
            ],
        ),
        ("mean_combine-partial-0c98c5a4517f58f8268985e7464daace", 0, 0, 0, 0): (
            f,
            [
                ("mean_chunk-540e88b7d9289f6b5461b95a0787af3e", 0, 0, 0, 0),
                ("mean_chunk-540e88b7d9289f6b5461b95a0787af3e", 1, 0, 0, 0),
                ("mean_chunk-540e88b7d9289f6b5461b95a0787af3e", 2, 0, 0, 0),
                ("mean_chunk-540e88b7d9289f6b5461b95a0787af3e", 3, 0, 0, 0),
            ],
        ),
        ("mean_combine-partial-0c98c5a4517f58f8268985e7464daace", 1, 0, 0, 0): (
            "mean_chunk-mean_combine-partial-0c98c5a4517f58f8268985e7464daace",
            1,
            0,
            0,
            0,
        ),
        ("mean_chunk-540e88b7d9289f6b5461b95a0787af3e", 0, 0, 0, 0): (
            f,
            ("random_sample-e16bcfb15a013023c98a21e2f03d66a9", 0, 0, 0, 0),
        ),
        ("mean_chunk-540e88b7d9289f6b5461b95a0787af3e", 1, 0, 0, 0): (
            f,
            ("random_sample-e16bcfb15a013023c98a21e2f03d66a9", 1, 0, 0, 0),
        ),
        ("mean_chunk-540e88b7d9289f6b5461b95a0787af3e", 2, 0, 0, 0): (
            f,
            ("random_sample-e16bcfb15a013023c98a21e2f03d66a9", 2, 0, 0, 0),
        ),
        ("mean_chunk-540e88b7d9289f6b5461b95a0787af3e", 3, 0, 0, 0): (
            f,
            ("random_sample-e16bcfb15a013023c98a21e2f03d66a9", 3, 0, 0, 0),
        ),
        ("random_sample-e16bcfb15a013023c98a21e2f03d66a9", 0, 0, 0, 0): (f, 1),
        ("random_sample-e16bcfb15a013023c98a21e2f03d66a9", 1, 0, 0, 0): (f, 1),
        ("random_sample-e16bcfb15a013023c98a21e2f03d66a9", 2, 0, 0, 0): (f, 1),
        ("random_sample-e16bcfb15a013023c98a21e2f03d66a9", 3, 0, 0, 0): (f, 1),
        ("random_sample-e16bcfb15a013023c98a21e2f03d66a9", 4, 0, 0, 0): (f, 1),
        (
            "mean_chunk-mean_combine-partial-17c7b5c6eed42e203858b3f6dde16003",
            1,
            0,
            0,
            0,
        ): (
            f,
            [
                (
                    f,
                    ("random_sample-e16bcfb15a013023c98a21e2f03d66a9", 4, 0, 0, 0),
                    ("random_sample-02eaa4a8dbb23fac4db22ad034c401b3", 4, 0, 0, 0),
                )
            ],
        ),
        (
            "mean_chunk-mean_combine-partial-0c98c5a4517f58f8268985e7464daace",
            1,
            0,
            0,
            0,
        ): (
            f,
            [(f, ("random_sample-e16bcfb15a013023c98a21e2f03d66a9", 4, 0, 0, 0), 2)],
        ),
        (
            "mean_chunk-mean_combine-partial-23adb4747560e6e33afd63c5bb179709",
            1,
            0,
            0,
            0,
        ): (
            f,
            [(f, ("random_sample-02eaa4a8dbb23fac4db22ad034c401b3", 4, 0, 0, 0), 2)],
        ),
    }

    other = {
        ("mean_agg-aggregate-e79dd3b9757c9fb2ad7ade96f3f6c814", 0, 0, 0): (
            f,
            [
                ("mean_combine-partial-e7d9fd7c132e12007a4b4f62ce443a75", 0, 0, 0, 0),
                ("mean_combine-partial-e7d9fd7c132e12007a4b4f62ce443a75", 1, 0, 0, 0),
            ],
        ),
        ("mean_combine-partial-e7d9fd7c132e12007a4b4f62ce443a75", 0, 0, 0, 0): (
            f,
            [
                ("mean_chunk-0df65d9a6e168673f32082f59f19576a", 0, 0, 0, 0),
                ("mean_chunk-0df65d9a6e168673f32082f59f19576a", 1, 0, 0, 0),
                ("mean_chunk-0df65d9a6e168673f32082f59f19576a", 2, 0, 0, 0),
                ("mean_chunk-0df65d9a6e168673f32082f59f19576a", 3, 0, 0, 0),
            ],
        ),
        ("mean_combine-partial-e7d9fd7c132e12007a4b4f62ce443a75", 1, 0, 0, 0): (
            "mean_chunk-mean_combine-partial-e7d9fd7c132e12007a4b4f62ce443a75",
            1,
            0,
            0,
            0,
        ),
        ("mean_chunk-0df65d9a6e168673f32082f59f19576a", 0, 0, 0, 0): (
            f,
            ("random_sample-a155d5a37ac5e09ede89c98a3bfcadff", 0, 0, 0, 0),
            ("random_sample-241fdbadc062900adc59d1a79c4c41e1", 0, 0, 0, 0),
        ),
        ("mean_chunk-0df65d9a6e168673f32082f59f19576a", 1, 0, 0, 0): (
            f,
            ("random_sample-a155d5a37ac5e09ede89c98a3bfcadff", 1, 0, 0, 0),
            ("random_sample-241fdbadc062900adc59d1a79c4c41e1", 1, 0, 0, 0),
        ),
        ("mean_chunk-0df65d9a6e168673f32082f59f19576a", 2, 0, 0, 0): (
            f,
            ("random_sample-a155d5a37ac5e09ede89c98a3bfcadff", 2, 0, 0, 0),
            ("random_sample-241fdbadc062900adc59d1a79c4c41e1", 2, 0, 0, 0),
        ),
        ("mean_chunk-0df65d9a6e168673f32082f59f19576a", 3, 0, 0, 0): (
            f,
            ("random_sample-a155d5a37ac5e09ede89c98a3bfcadff", 3, 0, 0, 0),
            ("random_sample-241fdbadc062900adc59d1a79c4c41e1", 3, 0, 0, 0),
        ),
        ("mean_agg-aggregate-c7647920facf0e557f947b7a6626b7be", 0, 0, 0): (
            f,
            [
                ("mean_combine-partial-57413f0bb18da78db0f689a096c7fbbf", 0, 0, 0, 0),
                ("mean_combine-partial-57413f0bb18da78db0f689a096c7fbbf", 1, 0, 0, 0),
            ],
        ),
        ("mean_combine-partial-57413f0bb18da78db0f689a096c7fbbf", 0, 0, 0, 0): (
            f,
            [
                ("mean_chunk-d6bd425ea61739f1eaa71762fe3bbbd7", 0, 0, 0, 0),
                ("mean_chunk-d6bd425ea61739f1eaa71762fe3bbbd7", 1, 0, 0, 0),
                ("mean_chunk-d6bd425ea61739f1eaa71762fe3bbbd7", 2, 0, 0, 0),
                ("mean_chunk-d6bd425ea61739f1eaa71762fe3bbbd7", 3, 0, 0, 0),
            ],
        ),
        ("mean_combine-partial-57413f0bb18da78db0f689a096c7fbbf", 1, 0, 0, 0): (
            "mean_chunk-mean_combine-partial-57413f0bb18da78db0f689a096c7fbbf",
            1,
            0,
            0,
            0,
        ),
        ("mean_chunk-d6bd425ea61739f1eaa71762fe3bbbd7", 0, 0, 0, 0): (
            f,
            ("random_sample-241fdbadc062900adc59d1a79c4c41e1", 0, 0, 0, 0),
            2,
        ),
        ("mean_chunk-d6bd425ea61739f1eaa71762fe3bbbd7", 1, 0, 0, 0): (
            f,
            ("random_sample-241fdbadc062900adc59d1a79c4c41e1", 1, 0, 0, 0),
            2,
        ),
        ("mean_chunk-d6bd425ea61739f1eaa71762fe3bbbd7", 2, 0, 0, 0): (
            f,
            ("random_sample-241fdbadc062900adc59d1a79c4c41e1", 2, 0, 0, 0),
            2,
        ),
        ("mean_chunk-d6bd425ea61739f1eaa71762fe3bbbd7", 3, 0, 0, 0): (
            f,
            ("random_sample-241fdbadc062900adc59d1a79c4c41e1", 3, 0, 0, 0),
            2,
        ),
        ("random_sample-241fdbadc062900adc59d1a79c4c41e1", 0, 0, 0, 0): (
            f,
            "random_sample",
            (10, 1, 987, 1920),
            [],
        ),
        ("random_sample-241fdbadc062900adc59d1a79c4c41e1", 1, 0, 0, 0): (
            f,
            "random_sample",
            (10, 1, 987, 1920),
            [],
        ),
        ("random_sample-241fdbadc062900adc59d1a79c4c41e1", 2, 0, 0, 0): (
            f,
            "random_sample",
            (10, 1, 987, 1920),
            [],
        ),
        ("random_sample-241fdbadc062900adc59d1a79c4c41e1", 3, 0, 0, 0): (
            f,
            "random_sample",
            (10, 1, 987, 1920),
            [],
        ),
        ("random_sample-241fdbadc062900adc59d1a79c4c41e1", 4, 0, 0, 0): (
            f,
            "random_sample",
            (10, 1, 987, 1920),
            [],
        ),
        ("mean_agg-aggregate-05071ebaabb68a64c180f6f443c5c8f4", 0, 0, 0): (
            f,
            [
                ("mean_combine-partial-a7c475f79a46af4265b189ffdc000bb3", 0, 0, 0, 0),
                ("mean_combine-partial-a7c475f79a46af4265b189ffdc000bb3", 1, 0, 0, 0),
            ],
        ),
        ("mean_combine-partial-a7c475f79a46af4265b189ffdc000bb3", 0, 0, 0, 0): (
            f,
            [
                ("mean_chunk-fd17feaf0728ea7a89d119d3fd172c75", 0, 0, 0, 0),
                ("mean_chunk-fd17feaf0728ea7a89d119d3fd172c75", 1, 0, 0, 0),
                ("mean_chunk-fd17feaf0728ea7a89d119d3fd172c75", 2, 0, 0, 0),
                ("mean_chunk-fd17feaf0728ea7a89d119d3fd172c75", 3, 0, 0, 0),
            ],
        ),
        ("mean_combine-partial-a7c475f79a46af4265b189ffdc000bb3", 1, 0, 0, 0): (
            "mean_chunk-mean_combine-partial-a7c475f79a46af4265b189ffdc000bb3",
            1,
            0,
            0,
            0,
        ),
        ("mean_chunk-fd17feaf0728ea7a89d119d3fd172c75", 0, 0, 0, 0): (
            f,
            ("random_sample-a155d5a37ac5e09ede89c98a3bfcadff", 0, 0, 0, 0),
            2,
        ),
        ("mean_chunk-fd17feaf0728ea7a89d119d3fd172c75", 1, 0, 0, 0): (
            f,
            ("random_sample-a155d5a37ac5e09ede89c98a3bfcadff", 1, 0, 0, 0),
            2,
        ),
        ("mean_chunk-fd17feaf0728ea7a89d119d3fd172c75", 2, 0, 0, 0): (
            f,
            ("random_sample-a155d5a37ac5e09ede89c98a3bfcadff", 2, 0, 0, 0),
            2,
        ),
        ("mean_chunk-fd17feaf0728ea7a89d119d3fd172c75", 3, 0, 0, 0): (
            f,
            ("random_sample-a155d5a37ac5e09ede89c98a3bfcadff", 3, 0, 0, 0),
            2,
        ),
        ("random_sample-a155d5a37ac5e09ede89c98a3bfcadff", 0, 0, 0, 0): (
            f,
            "random_sample",
            (10, 1, 987, 1920),
            [],
        ),
        ("random_sample-a155d5a37ac5e09ede89c98a3bfcadff", 1, 0, 0, 0): (
            f,
            "random_sample",
            (10, 1, 987, 1920),
            [],
        ),
        ("random_sample-a155d5a37ac5e09ede89c98a3bfcadff", 2, 0, 0, 0): (
            f,
            "random_sample",
            (10, 1, 987, 1920),
            [],
        ),
        ("random_sample-a155d5a37ac5e09ede89c98a3bfcadff", 3, 0, 0, 0): (
            f,
            "random_sample",
            (10, 1, 987, 1920),
            [],
        ),
        ("random_sample-a155d5a37ac5e09ede89c98a3bfcadff", 4, 0, 0, 0): (
            f,
            "random_sample",
            (10, 1, 987, 1920),
            [],
        ),
        (
            "mean_chunk-mean_combine-partial-a7c475f79a46af4265b189ffdc000bb3",
            1,
            0,
            0,
            0,
        ): (
            f,
            [(f, ("random_sample-a155d5a37ac5e09ede89c98a3bfcadff", 4, 0, 0, 0), 2)],
        ),
        (
            "mean_chunk-mean_combine-partial-e7d9fd7c132e12007a4b4f62ce443a75",
            1,
            0,
            0,
            0,
        ): (
            f,
            [
                (
                    f,
                    ("random_sample-a155d5a37ac5e09ede89c98a3bfcadff", 4, 0, 0, 0),
                    ("random_sample-241fdbadc062900adc59d1a79c4c41e1", 4, 0, 0, 0),
                )
            ],
        ),
        (
            "mean_chunk-mean_combine-partial-57413f0bb18da78db0f689a096c7fbbf",
            1,
            0,
            0,
            0,
        ): (
            f,
            [(f, ("random_sample-241fdbadc062900adc59d1a79c4c41e1", 4, 0, 0, 0), 2)],
        ),
    }
    first_pressure = max(diagnostics(first)[1])
    second_pressure = max(diagnostics(other)[1])
    assert first_pressure == second_pressure


def test_flox_reduction(abcde):
    a, b, c, d, e = abcde
    g = "g"
    dsk = {
        (a, 0): (f, 1),
        (a, 1, 0): (f, 1),
        (a, 1, 1): (f, 1),
        (a, 2, 0): (f, 1),
        (a, 2, 1): (f, 1),
        (b, 1, 0): (f, [(f, (a, 2, 0))]),
        (b, 1, 1): (f, [(f, (a, 2, 1))]),
        (b, 2, 1): (f, [(f, (a, 1, 1))]),
        (b, 2, 0): (f, [(f, (a, 1, 0))]),
        (b, 1, 1, 0): (b, 1, 0),
        (b, 1, 1, 1): (b, 1, 1),
        (b, 2, 2, 0): (b, 2, 0),
        (b, 2, 2, 1): (b, 2, 1),
        (c, 1, 0): (f, (b, 2, 2, 0)),
        (c, 1, 1): (f, (b, 2, 2, 1)),
        (c, 2, 0): (f, (b, 1, 1, 0)),
        (c, 2, 1): (f, (b, 1, 1, 1)),
        (d, 1): (f, [(f, (a, 1, 1), (a, 2, 1), (c, 1, 1), (c, 2, 1))]),
        (d, 0): (f, [(f, (a, 1, 0), (a, 2, 0), (c, 1, 0), (c, 2, 0))]),
        (e, 0): (d, 0),
        (e, 1): (d, 1),
        (g, 1, 0): (f, (a, 0), (b, 1, 1, 0)),
        (g, 1, 1): (f, (a, 0), (b, 2, 2, 0)),
        (g, 1, 2): (f, (a, 0), (e, 0)),
        (g, 2, 0): (f, (a, 0), (b, 1, 1, 1)),
        (g, 2, 1): (f, (a, 0), (b, 2, 2, 1)),
        (g, 2, 2): (f, (a, 0), (e, 1)),
    }
    o = order(dsk)
    assert_topological_sort(dsk, o)
    of1 = list(o[(g, 1, ix)] for ix in range(3))
    of2 = list(o[(g, 2, ix)] for ix in range(3))
    assert max(of1) < min(of2) or max(of2) < min(of1)


@pytest.mark.parametrize("optimize", [True, False])
@pytest.mark.parametrize("keep_self", [True, False])
@pytest.mark.parametrize("ndeps", [2, 5])
@pytest.mark.parametrize("n_reducers", [4, 7])
def test_reduce_with_many_common_dependents(optimize, keep_self, ndeps, n_reducers):
    pytest.importorskip("numpy")
    da = pytest.importorskip("dask.array")
    import numpy as np

    def random(**kwargs):
        assert len(kwargs) == ndeps
        return np.random.random((10, 10))

    trivial_deps = {
        f"k{i}": delayed(object(), name=f"object-{i}") for i in range(ndeps)
    }
    x = da.blockwise(
        random,
        "yx",
        new_axes={"y": (10,) * n_reducers, "x": (10,) * n_reducers},
        dtype=float,
        **trivial_deps,
    )
    graph = x.sum(axis=1, split_every=20)
    from dask.order import order

    if keep_self:
        # Keeping self adds a layer that cannot be fused
        dsk = (
            collections_to_expr([x, graph], optimize_graph=optimize)
            .optimize()
            .__dask_graph__()
        )
    else:
        dsk = (
            collections_to_expr([graph], optimize_graph=optimize)
            .optimize()
            .__dask_graph__()
        )
    dependencies, dependents = get_deps(dsk)
    # Verify assumptions
    before = len(dsk)
    before_dsk = dsk.copy()
    o = order(dsk)
    assert_topological_sort(dsk, o)
    assert before == len(o) == len(dsk)
    assert before_dsk == dsk
    # Verify assumptions (specifically that the reducers are sum-aggregate)
    assert ({"object", "sum", "sum-aggregate"}).issubset({key_split(k) for k in o})
    reducers = {k for k in o if key_split(k) == "sum-aggregate"}
    assert (p := max(diagnostics(dsk)[1])) <= n_reducers + ndeps, p
    # With optimize the metrics below change since there are many layers that
    # are being fused but the pressure above should already be a strong
    # indicator if something is wrong
    if optimize:
        for r in reducers:
            prios_deps = []
            for dep in dependencies[r]:
                prios_deps.append(o[dep])
            assert max(prios_deps) - min(prios_deps) == (len(dependencies[r]) - 1) * (
                2 if keep_self else 1
            )


def assert_topological_sort(dsk, order):
    dependencies, dependents = get_deps(dsk)
    num_needed = {k: len(dependencies[k]) for k in dsk}
    for k in sorted(dsk, key=order.__getitem__):
        assert num_needed[k] == 0
        for dep in dependents[k]:
            num_needed[dep] -= 1


def test_doublediff(abcde):
    a, b, c, d, e = abcde
    dsk = {
        (a, 3, 0, 1): (f, 1),
        (a, 3, 1, 1): (f, 1),
        (a, 3, 1, 0): (f, 1),
        (a, 2, 2, 2): (f, 1),
        (a, 4, 0, 1): (f, 1),
        (a, 4, 1, 0): (f, 1),
        (a, 4, 0, 0): (f, 1),
        (a, 5, 0, 0): (f, 1),
        (c, 0, 0, 0): (b, 1, 0, 0),
        (c, 0, 0, 1): (b, 1, 0, 1),
        (c, 0, 0, 2): (b, 1, 0, 2),
        (c, 0, 1, 0): (b, 1, 1, 0),
        (c, 0, 1, 1): (b, 1, 1, 1),
        (c, 0, 1, 2): (b, 1, 1, 2),
        (c, 0, 2, 0): (b, 1, 2, 0),
        (c, 0, 2, 1): (b, 1, 2, 1),
        (c, 1, 0, 0): (a, 5, 0, 0),
        (c, 1, 0, 1): (b, 0, 0, 1),
        (c, 1, 0, 2): (b, 0, 0, 2),
        (c, 1, 1, 0): (b, 0, 1, 0),
        (c, 1, 1, 1): (b, 0, 1, 1),
        (c, 1, 1, 2): (b, 0, 1, 2),
        (c, 1, 2, 0): (b, 0, 2, 0),
        (b, 1, 0, 0): (f, (a, 4, 0, 0)),
        (b, 1, 0, 1): (f, (a, 4, 0, 0)),
        (b, 1, 0, 2): (f, (a, 4, 0, 1)),
        (b, 1, 1, 0): (f, (a, 4, 0, 0)),
        (b, 1, 1, 1): (f, (a, 4, 0, 0)),
        (b, 1, 1, 2): (f, (a, 4, 0, 1)),
        (b, 1, 2, 0): (f, (a, 4, 1, 0)),
        (b, 1, 2, 1): (f, (a, 4, 1, 0)),
        (c, 1, 2, 1): (b, 0, 2, 1),
        (c, 1, 2, 2): (b, 0, 2, 2),
        (d, 0, 0, 0): (f, (c, 1, 0, 0), (c, 0, 0, 0)),
        (d, 0, 0, 1): (f, (c, 1, 0, 1), (c, 0, 0, 1)),
        (d, 0, 0, 2): (f, (c, 1, 0, 2), (c, 0, 0, 2)),
        (d, 0, 1, 0): (f, (c, 1, 1, 0), (c, 0, 1, 0)),
        (d, 0, 1, 1): (f, (c, 1, 1, 1), (c, 0, 1, 1)),
        (d, 0, 1, 2): (f, (c, 1, 1, 2), (c, 0, 1, 2)),
        (d, 0, 2, 0): (f, (c, 1, 2, 0), (c, 0, 2, 0)),
        (d, 0, 2, 1): (f, (c, 1, 2, 1), (c, 0, 2, 1)),
        (d, 0, 2, 2): (f, (c, 1, 2, 2), (a, 2, 2, 2)),
        (b, 0, 0, 1): (f, (a, 3, 0, 1)),
        (b, 0, 0, 2): (f, (a, 3, 0, 1)),
        (b, 0, 1, 0): (f, (a, 3, 1, 0)),
        (b, 0, 1, 1): (f, (a, 3, 1, 1)),
        (b, 0, 1, 2): (f, (a, 3, 1, 1)),
        (b, 0, 2, 0): (f, (a, 3, 1, 0)),
        (b, 0, 2, 1): (f, (a, 3, 1, 1)),
        (b, 0, 2, 2): (f, (a, 3, 1, 1)),
        (e, 0, 0, 0): (f, (d, 0, 0, 0)),
        (e, 0, 0, 1): (f, (d, 0, 0, 1)),
        (e, 0, 0, 2): (f, (d, 0, 0, 2)),
        (e, 0, 1, 0): (f, (d, 0, 1, 0)),
        (e, 0, 1, 1): (f, (d, 0, 1, 1)),
        (e, 0, 1, 2): (f, (d, 0, 1, 2)),
        (e, 0, 2, 0): (f, (d, 0, 2, 0)),
        (e, 0, 2, 1): (f, (d, 0, 2, 1)),
        (e, 0, 2, 2): (f, (d, 0, 2, 2)),
        "END": [
            (e, 0, 0, 0),
            (e, 0, 0, 1),
            (e, 0, 0, 2),
            (e, 0, 1, 0),
            (e, 0, 1, 1),
            (e, 0, 1, 2),
            (e, 0, 2, 0),
            (e, 0, 2, 1),
            (e, 0, 2, 2),
        ],
    }
    _, pressure = diagnostics(dsk)
    assert max(pressure) <= 11, max(pressure)


def test_recursion_depth_long_linear_chains():
    dsk = {"-1": (f, 1)}
    for ix in range(10000):
        dsk[str(ix)] = (f, str(ix - 1))
    assert len(order(dsk)) == len(dsk)


def test_gh_3055_explicit(abcde):
    # This is a subgraph extracted from gh_3055
    # From a critical path perspective, the root a, 2 only has to be loaded
    # towards the very end. However, loading it so late means that we are
    # blocking many possible reductions which is bad for memory pressure
    a, b, c, d, e = abcde
    g = "g"
    dsk = {
        ("root", 0): (f, 1),
        (a, 0): (f, ("root", 0)),
        (a, 1): (f, 1),
        (a, 2): (f, 1),
        (a, 3): (f, 1),
        (a, 4): (f, 1),
        (b, 0, 0): (f, (a, 0)),
        (b, 0, 1): (f, (a, 0)),
        (c, 0, 1): (f, (a, 0)),
        (b, 1, 0): (f, (a, 1)),
        (b, 1, 2): (f, (a, 1)),
        (c, 0, 0): (f, (b, 0), (a, 2), (a, 1)),
        (d, 0, 0): (f, (c, 0, 1), (c, 0, 0)),
        (d, 0, 1): (f, (c, 0, 1), (c, 0, 0)),
        (f, 1, 1): (f, (d, 0, 1)),
        (c, 1, 0): (f, (b, 1, 0), (b, 1, 2)),
        (c, 0, 2): (f, (b, 0, 0), (b, 0, 1)),
        (e, 0): (f, (c, 1, 0), (c, 0, 2)),
        (g, 1): (f, (e, 0), (a, 3)),
        (g, 2): (f, (g, 1), (a, 4), (d, 0, 0)),
    }
    dependencies, dependents = get_deps(dsk)
    con_r, _ = _connecting_to_roots(dependencies, dependents)
    assert len(con_r) == len(dsk)
    assert con_r[(e, 0)] == {("root", 0), (a, 1)}
    o = order(dsk)
    assert_topological_sort(dsk, o)
    assert max(diagnostics(dsk, o=o)[1]) <= 5
    assert o[(e, 0)] < o[(a, 3)] < o[(a, 4)]
    assert o[(a, 2)] < o[(a, 3)] < o[(a, 4)]


def test_order_flox_reduction_2(abcde):
    # https://github.com/dask/dask/issues/10618
    a, b, c, d, e = abcde
    dsk = {
        (a, 0): 0,
        (a, 1): 0,
        (a, 2): 0,
        (b, 0, 0, 0): (f, (a, 0)),
        (b, 0, 0, 1): (f, (a, 1)),
        (b, 0, 0, 2): (f, (a, 2)),
        (b, 0, 1, 0): (f, (a, 0)),
        (b, 0, 1, 1): (f, (a, 1)),
        (b, 0, 1, 2): (f, (a, 2)),
        (b, 1, 0, 0): (f, (a, 0)),
        (b, 1, 0, 1): (f, (a, 1)),
        (b, 1, 0, 2): (f, (a, 2)),
        (b, 1, 1, 0): (f, (a, 0)),
        (b, 1, 1, 1): (f, (a, 1)),
        (b, 1, 1, 2): (f, (a, 2)),
        (c, 0, 0): (f, [(b, 0, 0, 0), (b, 0, 0, 1), (b, 0, 0, 2)]),
        (c, 0, 1): (f, [(b, 0, 1, 0), (b, 0, 1, 1), (b, 0, 1, 2)]),
        (c, 1, 0): (f, [(b, 1, 0, 0), (b, 1, 0, 1), (b, 1, 0, 2)]),
        (c, 1, 1): (f, [(b, 1, 1, 0), (b, 1, 1, 1), (b, 1, 1, 2)]),
        (d, 0, 0): (c, 0, 0),
        (d, 0, 1): (c, 0, 1),
        (d, 1, 0): (c, 1, 0),
        (d, 1, 1): (c, 1, 1),
    }
    o = order(dsk)
    assert_topological_sort(dsk, o)
    final_nodes = sorted(
        [(d, ix, jx) for ix in range(2) for jx in range(2)], key=o.__getitem__
    )
    for ix in range(1, len(final_nodes)):
        # This assumes that all the data tasks are scheduled first.
        # Then, there are exactly four dependencies to load for every final
        # task.
        assert o[final_nodes[ix]] - o[final_nodes[ix - 1]] == 5


def test_xarray_map_reduce_with_slicing():
    # See https://github.com/dask/distributed/pull/8818
    dsk = {
        ("transpose", 0, 0, 0): (f, ("groupby-agg", 0, 0, 0)),
        ("transpose", 0, 0, 1): (f, ("groupby-agg", 0, 0, 1)),
        ("transpose", 0, 1, 0): (f, ("groupby-agg", 0, 1, 0)),
        ("transpose", 0, 1, 1): (f, ("groupby-agg", 0, 1, 1)),
        ("groupby-chunk", 0, 0, 0): (f, ("array", 0), ("getitem", 0, 0, 0)),
        ("groupby-chunk", 0, 0, 1): (f, ("array", 1), ("getitem", 1, 0, 0)),
        ("groupby-chunk", 0, 1, 0): (f, ("array", 0), ("getitem", 0, 0, 1)),
        ("groupby-chunk", 0, 1, 1): (f, ("array", 1), ("getitem", 1, 0, 1)),
        ("groupby-chunk", 1, 0, 0): (f, ("array", 0), ("getitem", 0, 1, 0)),
        ("groupby-chunk", 1, 0, 1): (f, ("array", 1), ("getitem", 1, 1, 0)),
        ("groupby-chunk", 1, 1, 0): (f, ("array", 0), ("getitem", 0, 1, 1)),
        ("groupby-chunk", 1, 1, 1): (f, ("getitem", 1, 1, 1), ("array", 1)),
        ("getitem", 0, 0, 0): (f, ("open_data", 0, 0, 0)),
        ("getitem", 0, 0, 1): (f, ("open_data", 0, 0, 1)),
        ("getitem", 0, 1, 0): (f, ("open_data", 0, 1, 0)),
        ("getitem", 0, 1, 1): (f, ("open_data", 0, 1, 1)),
        ("getitem", 1, 0, 0): (f, ("open_data", 1, 0, 0)),
        ("getitem", 1, 0, 1): (f, ("open_data", 1, 0, 1)),
        ("getitem", 1, 1, 0): (f, ("open_data", 1, 1, 0)),
        ("getitem", 1, 1, 1): (f, ("open_data", 1, 1, 1)),
        "data": "a",
        ("array", 0): "b",
        ("array", 1): "c",
        ("open_data", 0, 1, 1): (f, "data"),
        ("groupby-agg", 0, 1, 0): (
            f,
            ("groupby-chunk", 1, 0, 1),
            ("groupby-chunk", 1, 0, 0),
        ),
        ("groupby-agg", 0, 1, 1): (
            f,
            ("groupby-chunk", 1, 1, 0),
            ("groupby-chunk", 1, 1, 1),
        ),
        ("open_data", 1, 0, 0): (f, "data"),
        ("open_data", 0, 1, 0): (f, "data"),
        ("open_data", 1, 0, 1): (f, "data"),
        ("open_data", 1, 1, 0): (f, "data"),
        ("groupby-agg", 0, 0, 1): (
            f,
            ("groupby-chunk", 0, 1, 0),
            ("groupby-chunk", 0, 1, 1),
        ),
        ("open_data", 0, 0, 1): (f, "data"),
        ("open_data", 1, 1, 1): (f, "data"),
        ("groupby-agg", 0, 0, 0): (
            f,
            ("groupby-chunk", 0, 0, 0),
            ("groupby-chunk", 0, 0, 1),
        ),
        ("open_data", 0, 0, 0): (f, "data"),
    }
    o = order(dsk)

    assert_topological_sort(dsk, o)
    final_nodes = sorted(
        [("transpose", 0, ix, jx) for ix in range(2) for jx in range(2)],
        key=o.__getitem__,
    )
    all_diffs = []
    for ix in range(1, len(final_nodes)):
        # This assumes that all the data tasks are scheduled first.
        # Then, there are exactly four dependencies to load for every final
        # task.
        all_diffs.append(o[final_nodes[ix]] - o[final_nodes[ix - 1]])

    assert set(all_diffs) == {8}

    _, pressure = diagnostics(dsk, o=o)
    assert max(pressure) <= 5


@pytest.mark.parametrize("use_longest_path", [True, False])
def test_xarray_rechunk_map_reduce_cohorts(use_longest_path):
    dsk = _rechunk_merge_graph()
    if use_longest_path:
        # ensure that we run through longes path True and False
        keys = [("open-dataset", i, 0, 0) for i in range(20, 35)]
        dsk.update({("dummy", 0): (f, keys)})
        dsk.update(dict.fromkeys(keys, (f,)))

    o = order(dsk)

    assert_topological_sort(dsk, o)
    _, pressure = diagnostics(dsk, o=o)
    # cut the dummy tasks in the end
    assert max(pressure[:99]) <= 7

    final_nodes = sorted(
        [("transpose", ix, jx, 0) for ix in range(2) for jx in range(2)],
        key=o.__getitem__,
    )
    all_diffs = []
    for ix in range(1, len(final_nodes)):
        all_diffs.append(o[final_nodes[ix]] - o[final_nodes[ix - 1]])

    # We process a big chunk first and then a small side-branch
    # before we repeat this for the next independent branch
    assert all_diffs == [10, 39, 10]


def test_xarray_8414():
    # https://github.com/pydata/xarray/issues/8414#issuecomment-1793860552
    np = pytest.importorskip("numpy")
    xr = pytest.importorskip("xarray")

    ds = xr.Dataset(
        data_vars=dict(
            a=(("x", "y"), np.arange(80).reshape(8, 10)),
            b=(("y"), np.arange(10)),
        )
    ).chunk(x=1)

    def f(ds):
        d = ds.a * ds.b
        return (d * ds.b).sum("y")

    graph = ds.map_blocks(f).__dask_graph__()

    # Two tasks always reduce to one leaf
    # plus a commonly shared root node makes three
    assert (p := max(diagnostics(graph)[1])) == 3, p

    # If the pressure is actually at three, there is not much that could go
    # wrong. Still, this asserts that the dependents of the shared dependency
    # are only loaded only after the reducers have been computed which given a
    # symmetrical graph, means that the steps between the dependents are the
    # same. If the xarray implementation changes or fusing enters, the step size
    # could grow or shrink but should still be the same everywhere
    dsk = dict(graph)
    dependencies, dependents = get_deps(dsk)
    roots = {k for k, v in dependencies.items() if not v}
    shared_roots = {r for r in roots if len(dependents[r]) > 1}
    assert len(shared_roots) == 1
    shared_root = shared_roots.pop()
    o = order(dsk)
    assert_topological_sort(dsk, o)

    previous = None
    step = 0
    assert len(dependents[shared_root]) > 2
    for dep in sorted(dependents[shared_root], key=o.__getitem__):
        if previous is not None:
            if step == 0:
                step = o[dep] - o[shared_root] - 1
                assert step > 1
            else:
                assert o[dep] == previous + step
        previous = o[dep]


def test_connecting_to_roots_single_root():
    dsk = {
        "a": (f, 1),
        "b1": (f, "a"),
        "b2": (f, "a"),
        "c": (f, "b1", "b2"),
    }
    dependencies, dependents = get_deps(dsk)
    connected_roots, max_dependents = _connecting_to_roots(dependencies, dependents)
    assert connected_roots == {k: {"a"} if k != "a" else set() for k in dsk}
    assert len({id(v) for v in connected_roots.values()}) == 2
    assert max_dependents == {
        "a": 2,
        "b1": 2,
        "b2": 2,
        "c": 2,
    }, max_dependents
    connected_roots, max_dependents = _connecting_to_roots(dependents, dependencies)
    assert connected_roots == {k: {"c"} if k != "c" else set() for k in dsk}
    assert len({id(v) for v in connected_roots.values()}) == 2
    assert max_dependents == {
        "a": 2,
        "b1": 2,
        "b2": 2,
        "c": 2,
    }, max_dependents


def test_connecting_to_roots_tree_reduction():
    dsk = {
        "a0": (f, 1),
        "a1": (f, 1),
        "a2": (f, 1),
        "a3": (f, 1),
        "b1": (f, "a0"),
        "b2": (f, "a1"),
        "b3": (f, "a2"),
        "b4": (f, "a3"),
        "c1": (f, "b1", "b2"),
        "c2": (f, "b3", "b4"),
        "d": (f, "c1", "c2"),
    }
    dependencies, dependents = get_deps(dsk)
    connected_roots, max_dependents = _connecting_to_roots(dependencies, dependents)
    assert connected_roots == {
        "a0": set(),
        "a1": set(),
        "a2": set(),
        "a3": set(),
        "b1": {"a0"},
        "b2": {"a1"},
        "b3": {"a2"},
        "b4": {"a3"},
        "c1": {"a1", "a0"},
        "c2": {"a3", "a2"},
        "d": {"a1", "a2", "a3", "a0"},
    }
    assert len({id(v) for v in connected_roots.values() if len(v) == 0}) == 1
    assert all(v == 1 for v in max_dependents.values())

    connected_roots, max_dependents = _connecting_to_roots(dependents, dependencies)
    assert connected_roots.pop("d") == set()
    assert all(v == {"d"} for v in connected_roots.values()), set(
        connected_roots.values()
    )
    assert all(v == 2 for v in max_dependents.values())


def test_connecting_to_roots_long_linear_chains():
    dsk = {
        "a0": (f, 1),
        "a1": (f, 1),
        "b0": (f, "a0"),
        "c0": (f, "b0"),
        "d0": (f, "c0"),
        "d1": (f, "c0"),
        "e1": (f, "d0", "d1"),
        "b1": (f, "a1"),
        "c1": (f, "b1"),
        "d2": (f, "c1"),
        "d3": (f, "c1"),
        "f": (f, "e1", "d2", "d3", "d1"),
    }
    dependencies, dependents = get_deps(dsk)
    connected_roots, _ = _connecting_to_roots(dependencies, dependents)

    id_a0 = id(connected_roots["b0"])
    for k in ["c0", "d0", "d1", "e1"]:
        assert id(connected_roots[k]) == id_a0

    id_a0 = id(connected_roots["b1"])
    for k in ["c1", "d2", "d3"]:
        assert id(connected_roots[k]) == id_a0


def test_connected_roots_deduplication_mem_usage():
    # see https://github.com/dask/dask/issues/11055
    np = pytest.importorskip("numpy")
    da = pytest.importorskip("dask.array")

    def rotate(X, n_iter=2):
        gamma = 1
        n_samples, n_modes = X.shape
        R = np.eye(n_modes)

        for _ in range(n_iter):
            basis = X @ R

            basis2 = basis * basis.conj()
            basis3 = basis2 * basis
            W = np.diag(np.sum(basis2, axis=0))
            alpha = gamma / n_samples

            transformed = X.conj().T @ (basis3 - (alpha * basis @ W))
            U, svals, VT = da.linalg.svd_compressed(
                transformed, k=n_modes, compute=False
            )
            R = U @ VT

        Xrot = X @ R
        return Xrot

    n_modes = 100

    x = da.random.random((10_000, 100), chunks=(5_000, -1))
    u, s, v = da.linalg.svd_compressed(x, k=n_modes, compute=False)

    u_rot = rotate(u)
    dsk = dict(u_rot.dask)
    dependencies, dependents = get_deps(dsk)
    connected_roots, max_dependents = _connecting_to_roots(dependencies, dependents)

    hashes = defaultdict(set)

    for v in connected_roots.values():
        hashes[tuple(sorted(v))].add(id(v))

    assert set(map(len, hashes.values())) == {1}


def test_connecting_to_roots_asym():
    dsk = {
        "a0": (f, 1),
        "a1": (f, 1),
        "a2": (f, 1),
        "a3": (f, 1),
        "a4": (f, 1),
        # Diamond
        "b1": (f, "a0", "a1"),
        "c1": (f, "b1"),
        "c2": (f, "b1"),
        "d1": (f, "c1", "c2"),
        # Multi stage reducers
        "b2": (f, "a2"),
        "b3": (f, "a3"),
        "c3": (f, "b3", "b2"),
        "d2": (f, "b3", "c3", "a4"),
    }
    dependencies, dependents = get_deps(dsk)
    connected_roots, max_dependents = _connecting_to_roots(dependencies, dependents)
    assert connected_roots == {
        "a0": set(),
        "a1": set(),
        "a2": set(),
        "a3": set(),
        "a4": set(),
        "b1": {"a0", "a1"},
        "c1": {"a0", "a1"},
        "c2": {"a0", "a1"},
        "d1": {"a0", "a1"},
        "b2": {"a2"},
        "b3": {"a3"},
        "c3": {"a2", "a3"},
        "d2": {"a2", "a3", "a4"},
    }
    # Max dependents is just pre-computed for performance but it is itself just
    # a derived property of connected_roots
    assert max_dependents == {
        k: max((len(dependents[r]) for r in v), default=len(dependents[k]))
        for k, v in connected_roots.items()
    }


def test_do_not_mutate_input():
    # Internally we may modify the graph but we don't want to mutate the
    # external dsk
    np = pytest.importorskip("numpy")
    dsk = {
        "a": "data",
        "b": (f, 1),
        "c": np.array([[1, 2], [3, 4]]),
        "d": ["a", "b", "c"],
        "e": (f, "d"),
    }
    dependencies, __build_class__ = get_deps(dsk)
    dependencies_copy = dependencies.copy()
    dsk_copy = dsk.copy()
    o = order(dsk, dependencies=dependencies)
    assert_topological_sort(dsk, o)
    assert dsk == dsk_copy
    assert dependencies == dependencies_copy


def test_stackstac():
    # see https://github.com/dask/dask/issues/11363
    # and https://github.com/dask/dask/pull/11367
    final_keys = [
        ("transpose-9c80c6b26d7ef5cfec7274ee9e6b091e", 0, 0, 0),
        ("transpose-9c80c6b26d7ef5cfec7274ee9e6b091e", 0, 1, 0),
    ]
    dsk = {
        ("transpose-9c80c6b26d7ef5cfec7274ee9e6b091e", 0, 0, 0): (
            f,
            ("vindex-merge", 0, 0, 0),
        ),
        ("transpose-9c80c6b26d7ef5cfec7274ee9e6b091e", 0, 1, 0): (
            f,
            ("vindex-merge", 0, 1, 0),
        ),
        ("vindex-merge", 0, 0, 0): (
            f,
            ("vindex-slice", 0, 0, 0),
            ("vindex-slice", 1, 0, 0),
            ("vindex-slice", 2, 0, 0),
            ("vindex-slice", 3, 0, 0),
        ),
        ("vindex-slice", 2, 0, 0): (f, ("getitem-vindex", 2, 0, 0)),
        ("vindex-merge", 0, 1, 0): (
            f,
            ("vindex-slice", 0, 0, 1),
            ("vindex-slice", 1, 0, 1),
            ("vindex-slice", 2, 0, 1),
            ("vindex-slice", 3, 0, 1),
        ),
        ("vindex-slice", 2, 0, 1): (f, ("getitem-vindex", 2, 0, 1)),
        ("vindex-slice", 1, 0, 0): (f, ("getitem-vindex", 1, 0, 0)),
        ("getitem-vindex", 2, 0, 1): (
            f,
            ("shuffle-split", 321),
            ("shuffle-split", 322),
        ),
        ("vindex-slice", 3, 0, 1): (f, ("getitem-vindex", 3, 0, 1)),
        ("vindex-slice", 0, 0, 1): (f, ("getitem-vindex", 0, 0, 1)),
        ("getitem-vindex", 3, 0, 1): (
            f,
            ("shuffle-split", 299),
            ("shuffle-split", 300),
        ),
        ("getitem-vindex", 2, 0, 0): (f, ("concatenate", 0, 1, 9, 1)),
        ("vindex-slice", 3, 0, 0): (f, ("getitem-vindex", 3, 0, 0)),
        ("shuffle-split", 322): (f, ("concatenate", 0, 1, 9, 1)),
        ("vindex-slice", 0, 0, 0): (f, ("getitem-vindex", 0, 0, 0)),
        ("getitem-vindex", 0, 0, 1): (
            f,
            ("shuffle-split", 341),
            ("shuffle-split", 342),
        ),
        ("vindex-slice", 1, 0, 1): (f, ("getitem-vindex", 1, 0, 1)),
        ("shuffle-split", 321): (f, ("concatenate-getitem", 321)),
        ("shuffle-split", 299): (f, ("concatenate-getitem", 299)),
        ("getitem-vindex", 3, 0, 0): (f, ("concatenate", 0, 1, 8, 1)),
        ("getitem-vindex", 0, 0, 0): (f, ("concatenate", 0, 1, 10, 0)),
        ("getitem-vindex", 1, 0, 0): (f, ("concatenate", 0, 1, 9, 0)),
        ("concatenate-getitem", 299): (f, ("fetch-raster", 0, 0, 8, 1)),
        ("concatenate", 0, 1, 9, 1): (f, ("getitem-sub", 0, 1, 9, 1)),
        ("concatenate-getitem", 321): (f, ("fetch-raster", 0, 0, 9, 1)),
        ("fetch-raster", 0, 0, 8, 1): (f, ("asset_table", 0, 0)),
        ("concatenate", 0, 1, 8, 1): (f, ("getitem-sub", 0, 1, 8, 1)),
        ("fetch-raster", 0, 0, 9, 1): (f, ("asset_table", 0, 0)),
        ("concatenate", 0, 1, 9, 0): (f, ("getitem-sub", 0, 1, 9, 0)),
        ("shuffle-split", 300): (f, ("concatenate", 0, 1, 8, 1)),
        ("concatenate", 0, 1, 10, 0): (f, ("getitem-sub", 0, 1, 10, 0)),
        ("asset_table", 0, 0): (f, ("asset-table-data", 0, 0)),
        ("shuffle-split", 342): (f, ("concatenate", 0, 1, 10, 0)),
        ("getitem-vindex", 1, 0, 1): (
            f,
            ("shuffle-split", 319),
            ("shuffle-split", 320),
        ),
        ("getitem-sub", 0, 1, 10, 0): (f, ("fetch-raster", 0, 0, 10, 0)),
        ("getitem-sub", 0, 1, 8, 1): (f, ("fetch-raster", 0, 0, 8, 1)),
        ("shuffle-split", 341): (f, ("concatenate-getitem", 341)),
        ("getitem-sub", 0, 1, 9, 1): (f, ("fetch-raster", 0, 0, 9, 1)),
        ("concatenate-getitem", 341): (f, ("fetch-raster", 0, 0, 10, 0)),
        ("asset-table-data", 0, 0): (f,),
        ("getitem-sub", 0, 1, 9, 0): (f, ("fetch-raster", 0, 0, 9, 0)),
        ("fetch-raster", 0, 0, 10, 0): (f, ("asset_table", 0, 0)),
        ("shuffle-split", 320): (f, ("concatenate", 0, 1, 9, 0)),
        ("shuffle-split", 319): (f, ("concatenate-getitem", 319)),
        ("fetch-raster", 0, 0, 9, 0): (f, ("asset_table", 0, 0)),
        ("concatenate-getitem", 319): (f, ("fetch-raster", 0, 0, 9, 0)),
    }
    o = order(dsk, return_stats=True)
    _, pressure = diagnostics(dsk)
    assert max(pressure) <= 9, pressure
    for k in final_keys:
        # Ensure that we're not processing the entire graph using
        # process_runnables (fractional values) but are using the critical path
        assert o[k].critical_path in {1, 2}


def test_handle_out_of_graph_dependencies():
    from dask._task_spec import Task

    ta = Task("a", f, ())
    tb = Task("b", f, (ta.ref(),))
    tc = Task("c", f, (tb.ref(),))
    dsk = {t.key: t for t in [tb, tc]}
    o = order(dsk)
    assert len(o) == 2
