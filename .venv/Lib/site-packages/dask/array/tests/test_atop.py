from __future__ import annotations

import collections
from operator import add

import numpy as np
import pytest

import dask
import dask.array as da
from dask._task_spec import Task, TaskRef
from dask.array.utils import assert_eq
from dask.blockwise import (
    _BLOCKWISE_DEFAULT_PREFIX,
    Blockwise,
    _unique_dep,
    index_subs,
    optimize_blockwise,
    rewrite_blockwise,
)
from dask.highlevelgraph import HighLevelGraph
from dask.utils_test import dec, hlg_layer_topological, inc

a, b, c, d, e, f, g = "a", "b", "c", "d", "e", "f", "g"
_0, _1, _2, _3, _4, _5, _6, _7, _8, _9 = tuple(
    map(TaskRef, (f"{_BLOCKWISE_DEFAULT_PREFIX}{i}" for i in range(10)))
)
i, j, k = "i", "j", "k"


@pytest.mark.parametrize(
    "inputs,expected",
    [
        # output name, output index, task, input indices
        # Case 0
        [
            [(b, "i", Task(b, inc, _0), [(a, "i")])],
            (
                b,
                "i",
                Task(b, inc, _0),
                [
                    (a, "i"),
                ],
            ),
        ],
        # Case 1
        [
            [
                (b, "i", Task(b, inc, _0), [(a, "i")]),
                (c, "i", Task(c, dec, _0), [(a, "i")]),
                (d, "i", Task(d, add, _0, _1, _2), [(a, "i"), (b, "i"), (c, "i")]),
            ],
            (
                d,
                "i",
                Task.fuse(
                    Task(
                        d,
                        add,
                        _0,
                        TaskRef(_unique_dep(b, "i")),
                        TaskRef(_unique_dep(c, "i")),
                    ),
                    Task(_unique_dep(c, "i"), dec, _0),
                    Task(_unique_dep(b, "i"), inc, _0),
                ),
                [(a, "i")],
            ),
        ],
        # Case 2
        [
            [
                (b, "i", Task(b, inc, _0), [(a, "i")]),
                (c, "j", Task(c, inc, _0), [(b, "j")]),
            ],
            (
                c,
                "j",
                Task.fuse(
                    Task(c, inc, TaskRef(_unique_dep(b, "j"))),
                    Task(_unique_dep(b, "j"), inc, _0),
                ),
                [(a, "j")],
            ),
        ],
        # Case 3
        [
            [
                (b, "i", Task(b, sum, _0), [(a, "ij")]),
                (c, "k", Task(c, inc, _0), [(b, "k")]),
            ],
            (
                c,
                "k",
                Task.fuse(
                    Task(c, inc, TaskRef(_unique_dep(b, "k"))),
                    Task(_unique_dep(b, "k"), sum, _0),
                ),
                [(a, "kA")],
            ),
        ],
        # Case 4
        [
            [
                (c, "i", Task(c, inc, _0), [(a, "i")]),
                (d, "i", Task(d, inc, _0), [(b, "i")]),
                (g, "ij", Task(g, add, _0, _1), [(c, "i"), (d, "j")]),
            ],
            (
                g,
                "ij",
                Task.fuse(
                    Task(
                        g,
                        add,
                        TaskRef(_unique_dep(c, "i")),
                        TaskRef(_unique_dep(d, "j")),
                    ),
                    Task(_unique_dep(c, "i"), inc, _0),
                    Task(_unique_dep(d, "j"), inc, _1),
                ),
                [(a, "i"), (b, "j")],
            ),
        ],
        # Case 5
        [
            [
                (b, "ji", Task(b, np.transpose, _0), [(a, "ij")]),
                (c, "ij", Task(c, add, _0, _1), [(a, "ij"), (b, "ij")]),
            ],
            (
                c,
                "ij",
                Task.fuse(
                    Task(c, add, _0, TaskRef(_unique_dep(b, "ij"))),
                    Task(_unique_dep(b, "ij"), np.transpose, _1),
                ),
                [(a, "ij"), (a, "ji")],
            ),
        ],
        # Case 6
        [
            [
                (c, "i", Task(c, add, _0, _1), [(a, "i"), (b, "i")]),
                (d, "i", Task(d, inc, _0), [(c, "i")]),
            ],
            (
                d,
                "i",
                Task.fuse(
                    Task(d, inc, TaskRef(_unique_dep(c, "i"))),
                    Task(_unique_dep(c, "i"), add, _0, _1),
                ),
                [(a, "i"), (b, "i")],
            ),
        ],
        # Case 7
        [
            [
                (b, "ij", Task(b, np.transpose, _0), [(a, "ji")]),
                (d, "ij", Task(d, np.dot, _0, _1), [(b, "ik"), (c, "kj")]),
            ],
            (
                d,
                "ij",
                Task.fuse(
                    Task(
                        d,
                        np.dot,
                        TaskRef(_unique_dep(b, "ik")),
                        _0,
                    ),
                    Task(_unique_dep(b, "ik"), np.transpose, _1),
                ),
                [(c, "kj"), (a, "ki")],
            ),
        ],
        # Case 8
        [
            [
                (c, "i", Task(c, add, _0, _1), [(a, "i"), (b, "i")]),
                (f, "i", Task(f, add, _0, _1), [(d, "i"), (e, "i")]),
                (g, "i", Task(g, add, _0, _1), [(c, "i"), (f, "i")]),
            ],
            (
                g,
                "i",
                Task.fuse(
                    Task(
                        g,
                        add,
                        TaskRef(_unique_dep(c, "i")),
                        TaskRef(_unique_dep(f, "i")),
                    ),
                    Task(_unique_dep(f, "i"), add, _2, _3),
                    Task(_unique_dep(c, "i"), add, _0, _1),
                ),
                [(a, i), (b, i), (d, i), (e, i)],
            ),
        ],
        # Case 9
        [
            [
                (c, "i", Task(c, add, _0, _1), [(a, "i"), (b, "i")]),
                (f, "i", Task(f, add, _0, _1), [(a, "i"), (e, "i")]),
                (g, "i", Task(g, add, _0, _1), [(c, "i"), (f, "i")]),
            ],
            (
                g,
                "i",
                Task.fuse(
                    Task(
                        g,
                        add,
                        TaskRef(_unique_dep(c, "i")),
                        TaskRef(_unique_dep(f, "i")),
                    ),
                    Task(_unique_dep(f, "i"), add, _0, _2),
                    Task(_unique_dep(c, "i"), add, _0, _1),
                ),
                [(a, "i"), (b, "i"), (e, "i")],
            ),
        ],
        # Case 10
        [
            [
                (b, "i", Task(b, sum, _0), [(a, "ij")]),
                (c, "i", Task(c, inc, _0), [(b, "i")]),
            ],
            (
                c,
                "i",
                Task.fuse(
                    Task(c, inc, TaskRef(_unique_dep(b, "i"))),
                    Task(_unique_dep(b, "i"), sum, _0),
                ),
                [(a, "iA")],
            ),
        ],
        # Case 11
        [
            [
                (c, "i", Task(c, inc, _0), [(b, "i")]),
                (d, "i", Task(d, add, _0, _1, _2), [(a, "i"), (b, "i"), (c, "i")]),
            ],
            (
                d,
                "i",
                Task.fuse(
                    Task(
                        d,
                        add,
                        _0,
                        _1,
                        TaskRef(_unique_dep(c, "i")),
                    ),
                    Task(_unique_dep(c, "i"), inc, _1),
                ),
                [(a, "i"), (b, "i")],
            ),
        ],
        # Case 12
        # Include literals
        [
            [(b, "i", Task(b, add, _0, _1), [(a, "i"), (123, None)])],
            (b, "i", Task(b, add, _0, _1), [(a, "i"), (123, None)]),
        ],
        # Case 13
        [
            [
                (b, "i", Task(b, add, _0, _1), [(a, "i"), (123, None)]),
                (c, "j", Task(c, add, _0, _1), [(b, "j"), (456, None)]),
            ],
            (
                c,
                "j",
                Task.fuse(
                    Task(c, add, TaskRef(_unique_dep(b, "j")), _0),
                    Task(_unique_dep(b, "j"), add, _1, _2),
                ),
                [(456, None), (a, "j"), (123, None)],
            ),
        ],
        # Case 14
        # Literals that compare equal (e.g. 0 and False) aren't deduplicated
        [
            [
                (b, "i", Task(b, add, _0, _1), [(a, "i"), (0, None)]),
                (c, "j", Task(c, add, _0, _1), [(b, "j"), (False, None)]),
            ],
            (
                c,
                "j",
                Task.fuse(
                    Task(c, add, TaskRef(_unique_dep(b, "j")), _0),
                    Task(_unique_dep(b, "j"), add, _1, _2),
                ),
                [(False, None), (a, "j"), (0, None)],
            ),
        ],
        # Case 15
        # Literals are deduplicated
        [
            [
                (b, "i", Task(b, add, _0, _1), [(a, "i"), (123, None)]),
                (c, "j", Task(c, add, _0, _1), [(b, "j"), (123, None)]),
            ],
            (
                c,
                "j",
                Task.fuse(
                    Task(c, add, TaskRef(_unique_dep(b, "j")), _0),
                    Task(_unique_dep(b, "j"), add, _1, _0),
                ),
                [(123, None), (a, "j")],
            ),
        ],
        # Case 16
        # Check cases where two distinct indices are used
        # for the same dependency name, and where the same
        # dependency-index combination is repeated
        # (See: https://github.com/dask/dask/issues/8535)
        [
            [
                (b, "jk", Task(b, add, _0, _1), [(a, "jk"), (2, None)]),
                (c, "ijk", Task(c, add, _0, _1), [(b, "ij"), (b, "jk")]),
                (d, "ijk", Task(d, inc, _0, _1), [(c, "ijk"), (123, None)]),
            ],
            (
                "d",
                "ijk",
                Task.fuse(
                    Task(d, inc, TaskRef(_unique_dep(c, "ijk")), _0),
                    Task(
                        _unique_dep(c, "ijk"),
                        add,
                        TaskRef(_unique_dep(b, "ij")),
                        TaskRef(_unique_dep(b, "jk")),
                    ),
                    Task(_unique_dep(b, "ij"), add, _1, _2),
                    Task(_unique_dep(b, "jk"), add, _3, _2),
                ),
                [(123, None), (a, "ij"), (2, None), (a, "jk")],
            ),
        ],
        # Case 17
        [
            [
                (b, "jk", Task(b, add, _0, _1), [(a, "jk"), (2, None)]),
                (c, "ijk", Task(c, add, _0, _1), [(b, "ij"), (b, "jk")]),
                (d, "ijk", Task(d, add, _0, _1, _2), [(b, "ij"), (c, "ij"), (b, "ij")]),
            ],
            (
                "d",
                "ijk",
                Task.fuse(
                    Task(
                        d,
                        add,
                        TaskRef(_unique_dep(b, "ij")),
                        TaskRef(_unique_dep(c, "ij")),
                        TaskRef(_unique_dep(b, "ij")),
                    ),
                    Task(
                        _unique_dep(c, "ij"),
                        add,
                        TaskRef(_unique_dep(b, "ij")),
                        TaskRef(_unique_dep(b, "jk")),
                    ),
                    Task(_unique_dep(b, "ij"), add, _0, _1),
                    Task(_unique_dep(b, "jk"), add, _2, _1),
                ),
                [(a, "ij"), (2, None), (a, "jk")],
            ),
        ],
    ],
)
def test_rewrite(inputs, expected):
    inputs = [
        Blockwise(
            *inp, numblocks={k: (1,) * len(v) for k, v in inp[-1] if v is not None}
        )
        for inp in inputs
    ]
    result = rewrite_blockwise(inputs)

    result2 = (
        result.output,
        "".join(result.output_indices),
        result.task,
        [
            (name, "".join(ind) if ind is not None else ind)
            for name, ind in result.indices
        ],
    )
    # Assert on the task first to get a more informative error message
    assert result2[2] == expected[2]
    assert result2 == expected


def test_index_subs():
    assert index_subs(tuple("ij"), {"i": "j", "j": "i"}) == tuple("ji")


def test_optimize_blockwise():
    x = da.ones(10, chunks=(5,))
    y = (((x + 1) + 2) + 3) + 4

    dsk = da.optimization.optimize_blockwise(y.dask)

    assert isinstance(dsk, HighLevelGraph)

    assert (
        len([layer for layer in dsk.layers.values() if isinstance(layer, Blockwise)])
        == 1
    )


def test_optimize_blockwise_control_annotations():
    """
    Can we fuse blockwise layers with different, but compatible
    annotations for retries, priority, etc.
    """

    a = da.ones(10, chunks=(5,))
    b = a + 1

    with dask.annotate(retries=5, workers=["a", "b", "c"], allow_other_workers=False):
        c = b + 2

    with dask.annotate(priority=2, workers=["b", "c", "d"], allow_other_workers=True):
        d = c + 3

    with dask.annotate(retries=3, resources={"GPU": 2, "Memory": 10}):
        e = d + 4

    with dask.annotate(priority=4, resources={"GPU": 5, "Memory": 4}):
        f = e + 5

    # This one will not be fused due to the custom annotation, nor will the one below
    with dask.annotate(foo="bar"):
        g = f + 6

    h = g + 6

    dsk = da.optimization.optimize_blockwise(h.dask)

    # The layers and their annotations should be fusable until the custom one
    assert len(dsk.layers) == 3
    layer = hlg_layer_topological(dsk, 0)  # First layer is the fused one
    annotations = layer.annotations

    assert len(annotations) == 5
    assert annotations["priority"] == 4  # max
    assert annotations["retries"] == 5  # max
    assert annotations["allow_other_workers"] is False  # More restrictive
    assert set(annotations["workers"]) == {"b", "c"}  # intersection
    assert annotations["resources"] == {"GPU": 5, "Memory": 10}  # Max of resources

    # If we disable blockwise annotation fusion, we can only fuse the first two layers.
    with dask.config.set({"optimization.annotations.fuse": False}):
        dsk = da.optimization.optimize_blockwise(h.dask)
        assert len(dsk.layers) == 7


def test_optimize_blockwise_custom_annotations():
    a = da.ones(10, chunks=(5,))
    b = a + 1

    with dask.annotate(qux="foo"):
        c = b + 2
        d = c + 3

    with dask.annotate(qux="baz"):
        e = d + 4
        f = e + 5

    g = f + 6

    dsk = da.optimization.optimize_blockwise(g.dask)

    annotations = (
        layer.annotations
        for layer in dsk.layers.values()
        if isinstance(layer, Blockwise)
    )
    annotations = collections.Counter(
        tuple(a.items()) if type(a) is dict else a for a in annotations
    )

    assert len(annotations) == 3
    assert annotations[None] == 2
    assert annotations[(("qux", "baz"),)] == 1
    assert annotations[(("qux", "foo"),)] == 1


def test_blockwise_diamond_fusion():
    x = da.ones(10, chunks=(5,))
    y = ((x + 1) + 2) + 3
    a = y * 2
    b = y * 3
    c = a + b
    d = ((c + 1) + 2) + 3

    dsk = da.optimization.optimize_blockwise(d.dask)
    assert isinstance(dsk, HighLevelGraph)

    assert (
        len([layer for layer in dsk.layers.values() if isinstance(layer, Blockwise)])
        == 1
    )


def test_blockwise_non_blockwise_output():
    x = da.ones(10, chunks=(5,))
    y = ((x + 1) + 2) + 3
    w = y.sum()
    z = ((y * 2) * 3) * 4

    z_top_before = tuple(z.dask.layers[z.name].indices)
    (zz,) = dask.optimize(z)
    z_top_after = tuple(z.dask.layers[z.name].indices)
    assert z_top_before == z_top_after, "z_top mutated"

    dsk = optimize_blockwise(z.dask, keys=list(dask.core.flatten(z.__dask_keys__())))
    assert isinstance(dsk, HighLevelGraph)
    assert (
        len([layer for layer in dsk.layers.values() if isinstance(layer, Blockwise)])
        == 1
    )

    dsk = optimize_blockwise(
        HighLevelGraph.merge(w.dask, z.dask),
        keys=list(dask.core.flatten([w.__dask_keys__(), z.__dask_keys__()])),
    )
    assert isinstance(dsk, HighLevelGraph)
    assert (
        len([layer for layer in z.dask.layers.values() if isinstance(layer, Blockwise)])
        >= 1
    )


def test_top_len():
    x = da.ones(10, chunks=(5,))
    y = x[:, None] * x[None, :]

    d = y.dask.layers[y.name]
    assert len(d) == 4


def test_inner_compute():
    x = da.ones(10, chunks=(5,)) + 1 + 2 + 3
    a = x.sum()
    y = x * 2 * 3 * 4
    b = y.sum()
    z = x * 2 * 3

    dask.compute(x, a, y, b, z)


@pytest.mark.parametrize("name", ["_", "_0", "_1", ".", ".0"])
def test_common_token_names_args(name):
    x = np.array(["a", "bb", "ccc"], dtype=object)
    d = da.from_array(x, chunks=2)

    result = da.blockwise(add, "i", d, "i", name, None, dtype=object)
    expected = x + name

    assert_eq(result, expected)


@pytest.mark.parametrize("name", ["_0", "_1", ".", ".0", "_"])
def test_common_token_names_kwargs(name):
    x = np.array(["a", "bb", "ccc"], dtype=object)
    d = da.from_array(x, chunks=2)

    result = da.blockwise(lambda x, y: x + y, "i", d, "i", y=name, dtype=object)
    expected = x + name

    assert_eq(result, expected)


def test_blockwise_names():
    x = da.ones(5, chunks=(2,))
    y = da.blockwise(add, "i", x, "i", dtype=x.dtype)
    assert y.name.startswith("add")


def test_blockwise_new_axes():
    def f(x):
        return x[:, None] * np.ones((1, 7))

    x = da.ones(5, chunks=2)
    y = da.blockwise(
        f, "aq", x, "a", new_axes={"q": 7}, concatenate=True, dtype=x.dtype
    )
    assert y.chunks == ((2, 2, 1), (7,))
    assert_eq(y, np.ones((5, 7)))

    def f(x):
        return x[None, :] * np.ones((7, 1))

    x = da.ones(5, chunks=2)
    y = da.blockwise(
        f, "qa", x, "a", new_axes={"q": 7}, concatenate=True, dtype=x.dtype
    )
    assert y.chunks == ((7,), (2, 2, 1))
    assert_eq(y, np.ones((7, 5)))

    def f(x):
        y = x.sum(axis=1)
        return y[:, None] * np.ones((1, 5))

    x = da.ones((4, 6), chunks=(2, 2))
    y = da.blockwise(
        f, "aq", x, "ab", new_axes={"q": 5}, concatenate=True, dtype=x.dtype
    )
    assert y.chunks == ((2, 2), (5,))
    assert_eq(y, np.ones((4, 5)) * 6)


def test_blockwise_new_axes_2():
    x = da.ones((2, 2), chunks=(1, 1))

    def func(x):
        return np.stack([x, -x], axis=-1)

    y = da.blockwise(
        func,
        ("x", "y", "sign"),
        x,
        ("x", "y"),
        dtype=x.dtype,
        concatenate=True,
        new_axes={"sign": 2},
    )

    assert_eq(y, y)


@pytest.mark.parametrize("concatenate", [True, False])
def test_blockwise_stacked_new_axes(concatenate):
    def f(x):
        return x[..., None] * np.ones((1, 7))

    x = da.ones(5, chunks=2)
    y = da.blockwise(
        f, "aq", x, "a", new_axes={"q": 7}, concatenate=concatenate, dtype=x.dtype
    )
    z = da.blockwise(
        f, "abq", y, "ab", new_axes={"q": 7}, concatenate=concatenate, dtype=x.dtype
    )
    assert z.chunks == ((2, 2, 1), (7,), (7,))
    assert_eq(z, np.ones((5, 7, 7)))


@pytest.mark.parametrize("concatenate", [True, False])
def test_blockwise_stacked_new_axes_front(concatenate):
    def f(x):
        if isinstance(x, list):
            x = np.concatenate(x)
        return x[None, ...] * np.ones(7)[(slice(None),) + (None,) * x.ndim]

    x = da.ones(5, chunks=2)
    y = da.blockwise(
        f, "qa", x, "a", new_axes={"q": 7}, concatenate=concatenate, dtype=x.dtype
    )
    z = da.blockwise(
        f, "qab", y, "ab", new_axes={"q": 7}, concatenate=concatenate, dtype=x.dtype
    )
    assert z.chunks == ((7,), (7,), (2, 2, 1))
    assert_eq(z, np.ones((7, 7, 5)))

    w = da.blockwise(
        lambda x: x[:, 0, 0], "a", z, "abc", dtype=x.dtype, concatenate=True
    )
    assert w.chunks == ((7,),)
    assert_eq(w, np.ones((7,)))


@pytest.mark.parametrize("concatenate", [True, False])
def test_blockwise_stacked_new_axes_same_dim(concatenate):
    def f(x):
        return x[..., None] * np.ones((1, 7))

    x = da.ones(5, chunks=2)
    y = da.zeros(5, chunks=2)
    a = da.blockwise(
        f, "aq", x, "a", new_axes={"q": 7}, concatenate=concatenate, dtype=x.dtype
    )
    b = da.blockwise(
        f, "aq", y, "a", new_axes={"q": 7}, concatenate=concatenate, dtype=x.dtype
    )
    c = a + b
    assert c.chunks == ((2, 2, 1), (7,))
    assert_eq(c, np.ones((5, 7)))


def test_blockwise_new_axes_chunked():
    def f(x):
        return x[None, :] * 2

    x = da.arange(0, 6, 1, chunks=2, dtype=np.int32)
    y = da.blockwise(f, "qa", x, "a", new_axes={"q": (1, 1)}, dtype=x.dtype)
    assert y.chunks == ((1, 1), (2, 2, 2))
    assert_eq(y, np.array([[0, 2, 4, 6, 8, 10], [0, 2, 4, 6, 8, 10]], np.int32))


def test_blockwise_no_args():
    def f():
        return np.ones((2, 3), np.float32)

    x = da.blockwise(f, "ab", new_axes={"a": 2, "b": (3, 3)}, dtype=np.float32)
    assert x.chunks == ((2,), (3, 3))
    assert_eq(x, np.ones((2, 6), np.float32))


def test_blockwise_no_array_args():
    def f(dtype):
        return np.ones((2, 3), dtype)

    x = da.blockwise(
        f, "ab", np.float32, None, new_axes={"a": 2, "b": (3, 3)}, dtype=np.float32
    )
    assert x.chunks == ((2,), (3, 3))
    assert_eq(x, np.ones((2, 6), np.float32))


def test_blockwise_kwargs():
    def f(a, b=0):
        return a + b

    x = da.ones(5, chunks=(2,))
    y = da.blockwise(f, "i", x, "i", b=10, dtype=x.dtype)
    assert_eq(y, np.ones(5) + 10)


def test_blockwise_chunks():
    x = da.ones((5, 5), chunks=((2, 1, 2), (3, 2)))

    def double(a, axis=0):
        return np.concatenate([a, a], axis=axis)

    y = da.blockwise(
        double,
        "ij",
        x,
        "ij",
        adjust_chunks={"i": lambda n: 2 * n},
        axis=0,
        dtype=x.dtype,
    )
    assert y.chunks == ((4, 2, 4), (3, 2))
    assert_eq(y, np.ones((10, 5)))

    y = da.blockwise(
        double,
        "ij",
        x,
        "ij",
        adjust_chunks={"j": lambda n: 2 * n},
        axis=1,
        dtype=x.dtype,
    )
    assert y.chunks == ((2, 1, 2), (6, 4))
    assert_eq(y, np.ones((5, 10)))

    x = da.ones((10, 10), chunks=(5, 5))
    y = da.blockwise(
        double, "ij", x, "ij", axis=0, adjust_chunks={"i": 10}, dtype=x.dtype
    )
    assert y.chunks == ((10, 10), (5, 5))
    assert_eq(y, np.ones((20, 10)))

    y = da.blockwise(
        double, "ij", x, "ij", axis=0, adjust_chunks={"i": (10, 10)}, dtype=x.dtype
    )
    assert y.chunks == ((10, 10), (5, 5))
    assert_eq(y, np.ones((20, 10)))


def test_blockwise_numpy_arg():
    x = da.arange(10, chunks=(5,))
    y = np.arange(1000)

    x = x.map_blocks(lambda x, y: x, 1.0)
    x = x.map_blocks(lambda x, y: x, "abc")
    x = x.map_blocks(lambda x, y: x, y)
    x = x.map_blocks(lambda x, y: x, "abc")
    x = x.map_blocks(lambda x, y: x, 1.0)
    x = x.map_blocks(lambda x, y, z: x, "abc", np.array(["a", "b"], dtype=object))
    assert_eq(x, np.arange(10))


def test_bag_array_conversion():
    import dask.bag as db

    b = db.range(10, npartitions=1)
    (x,) = b.map_partitions(np.asarray).to_delayed()
    (x,) = (da.from_delayed(a, shape=(10,), dtype=int) for a in [x])
    z = da.concatenate([x])
    assert_eq(z, np.arange(10), check_graph=False)


def test_svd():
    x = da.ones((1, 1), chunks=(1, 1))
    y = x * 2
    u, s, v = da.linalg.svd(y)
    z = y + u
    assert_eq(z, z)


def test_args_delayed():
    x = da.arange(10, chunks=(5,))
    y = dask.delayed(lambda: 100)()

    z = da.blockwise(add, "i", x, "i", y, None, dtype=x.dtype)
    assert_eq(z, np.arange(10) + 100)

    z = da.blockwise(lambda x, y: x + y, "i", x, "i", y=y, dtype=x.dtype)
    assert_eq(z, np.arange(10) + 100)


@pytest.mark.parametrize(
    "tup", [(1, 2), collections.namedtuple("foo", ["a", "b"])(1, 2)]  # type: ignore
)
def test_namedtuple(tup):
    A = da.random.default_rng().random((20, 20), chunks=(10, 10))

    def f(data, x):
        return data

    B = da.blockwise(f, ("d1", "d2"), A, ("d1", "d2"), x=tup, dtype=A.dtype)

    assert_eq(A, B)


def test_validate_top_inputs():
    A = da.random.default_rng().random((20, 20), chunks=(10, 10))

    with pytest.raises(ValueError) as info:
        da.blockwise(inc, "jk", A, "ij", dtype=A.dtype)

    assert "unknown dimension" in str(info.value).lower()
    assert "k" in str(info.value)
    assert "j" not in str(info.value)

    with pytest.raises(ValueError) as info:
        da.blockwise(inc, "ii", A, "ij", dtype=A.dtype)

    assert "repeated" in str(info.value).lower()
    assert "i" in str(info.value)


def test_dont_merge_before_reductions():
    x = da.ones(10, chunks=(5,))
    y = da.blockwise(inc, "i", x, "i", dtype=x.dtype)
    z = da.blockwise(sum, "", y, "i", dtype=y.dtype)
    w = da.blockwise(sum, "", z, "", dtype=y.dtype)

    dsk = optimize_blockwise(w.dask)

    assert len([d for d in dsk.layers.values() if isinstance(d, Blockwise)]) == 2

    z.compute()


def test_atop_legacy():
    x = da.ones(10, chunks=(5,))
    with pytest.warns(
        UserWarning, match="The da.atop function has moved to da.blockwise"
    ):
        y = da.atop(inc, "i", x, "i", dtype=x.dtype)
    z = da.blockwise(inc, "i", x, "i", dtype=x.dtype)
    assert_eq(y, z)
    assert y.name == z.name


def test_non_hlg():
    # Regression test for https://github.com/dask/dask/issues/5850
    a = da.from_array(np.ones(1, np.float64), chunks=(1,))
    a.dask = dict(a.dask)  # Convert from HighLevelGraph to plain dict
    b = da.from_array(np.zeros(1, np.float64), chunks=(1,))
    x = a + b
    assert_eq(x, a)
