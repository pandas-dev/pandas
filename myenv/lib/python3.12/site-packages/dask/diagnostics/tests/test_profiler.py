from __future__ import annotations

import contextlib
import os
import warnings
from operator import add, mul
from timeit import default_timer

import pytest

from dask.diagnostics import CacheProfiler, Profiler, ResourceProfiler
from dask.diagnostics.profile_visualize import BOKEH_VERSION
from dask.threaded import get
from dask.utils import apply, tmpfile
from dask.utils_test import slowadd

try:
    import bokeh
except ImportError:
    bokeh = None
try:
    import psutil
except ImportError:
    psutil = None  # type: ignore


prof = Profiler()

dsk = {"a": 1, "b": 2, "c": (add, "a", "b"), "d": (mul, "a", "b"), "e": (mul, "c", "d")}
dsk2 = {"a": 1, "b": 2, "c": (slowadd, "a", "b")}


def test_profiler():
    with prof:
        in_context_time = default_timer()
        out = get(dsk, "e")
    assert out == 6
    assert prof.start_time < in_context_time < prof.end_time
    prof_data = sorted(prof.results, key=lambda d: d.key)
    keys = [i.key for i in prof_data]
    assert keys == ["c", "d", "e"]
    tasks = [i.task for i in prof_data]
    assert tasks == [(add, "a", "b"), (mul, "a", "b"), (mul, "c", "d")]
    prof.clear()
    assert prof.results == []


def test_profiler_works_under_error():
    div = lambda x, y: x / y
    dsk = {"x": (div, 1, 1), "y": (div, "x", 2), "z": (div, "y", 0)}

    with contextlib.suppress(ZeroDivisionError):
        with prof:
            get(dsk, "z")

    assert all(len(v) == 5 for v in prof.results)
    assert len(prof.results) == 2


def test_two_gets():
    with prof:
        get(dsk, "e")
    n = len(prof.results)

    dsk2 = {"x": (add, 1, 2), "y": (add, "x", "x")}

    with prof:
        get(dsk2, "y")
    m = len(prof.results)

    with prof:
        get(dsk, "e")
        get(dsk2, "y")
        get(dsk, "e")

    assert len(prof.results) == n + m + n


@pytest.mark.skipif("not psutil")
def test_resource_profiler():
    with ResourceProfiler(dt=0.01) as rprof:
        in_context_time = default_timer()
        get(dsk2, "c")
    results = rprof.results
    assert len(results) > 0
    assert all(isinstance(i, tuple) and len(i) == 3 for i in results)
    assert rprof.start_time < in_context_time < rprof.end_time

    # Tracker stopped on exit
    assert not rprof._is_running()

    rprof.clear()
    assert rprof.results == []

    # Close is idempotent
    rprof.close()
    assert not rprof._is_running()

    # Restarts tracker if already closed
    with rprof:
        get(dsk2, "c")
    assert len(rprof.results) > 0


@pytest.mark.skipif("not psutil")
def test_resource_profiler_multiple_gets():
    with ResourceProfiler(dt=0.01) as rprof:
        get(dsk2, "c")
        assert len(rprof.results) == 0
        get(dsk2, "c")
    results = rprof.results
    assert all(isinstance(i, tuple) and len(i) == 3 for i in results)

    rprof.clear()
    rprof.register()
    get(dsk2, "c")
    assert len(rprof.results) > 0
    get(dsk2, "c")
    rprof.unregister()

    results = rprof.results
    assert all(isinstance(i, tuple) and len(i) == 3 for i in results)

    rprof.close()
    assert not rprof._is_running()


def test_cache_profiler():
    with CacheProfiler() as cprof:
        in_context_time = default_timer()
        get(dsk2, "c")
    results = cprof.results
    assert all(isinstance(i, tuple) and len(i) == 5 for i in results)
    assert cprof.start_time < in_context_time < cprof.end_time

    cprof.clear()
    assert cprof.results == []

    tics = [0]

    def nbytes(res):
        tics[0] += 1
        return tics[0]

    with CacheProfiler(nbytes) as cprof:
        get(dsk2, "c")

    results = cprof.results
    assert tics[-1] == len(results)
    assert tics[-1] == results[-1].metric
    assert cprof._metric_name == "nbytes"
    assert CacheProfiler(metric=nbytes, metric_name="foo")._metric_name == "foo"


@pytest.mark.parametrize(
    "profiler",
    [
        Profiler,
        pytest.param(
            lambda: ResourceProfiler(dt=0.01), marks=pytest.mark.skipif("not psutil")
        ),
        CacheProfiler,
    ],
)
def test_register(profiler):
    prof = profiler()
    try:
        prof.register()
        get(dsk2, "c")
        n = len(prof.results)
        assert n > 0
        get(dsk2, "c")
        assert len(prof.results) > n
    finally:
        prof.unregister()


@pytest.mark.skipif("not bokeh")
def test_unquote():
    from dask.diagnostics.profile_visualize import unquote

    t = {"a": 1, "b": 2, "c": 3}
    task_dask = (dict, [["a", 1], ["b", 2], ["c", 3]])
    assert unquote(task_dask) == t

    t = {"a": [1, 2, 3], "b": 2, "c": 3}
    task_dask = (dict, [["a", [1, 2, 3]], ["b", 2], ["c", 3]])
    assert unquote(task_dask) == t

    t = [1, 2, 3]
    task_dask = [1, 2, 3]
    assert unquote(task_dask) == t


@pytest.mark.skipif("not bokeh")
def test_pprint_task():
    from dask.diagnostics.profile_visualize import pprint_task

    keys = {"a", "b", "c", "d", "e"}
    assert pprint_task((add, "a", 1), keys) == "add(_, *)"
    assert pprint_task((add, (add, "a", 1)), keys) == "add(add(_, *))"
    res = "sum([*, _, add(_, *)])"
    assert pprint_task((sum, [1, "b", (add, "a", 1)]), keys) == res
    assert pprint_task((sum, (1, 2, 3, 4, 5, 6, 7)), keys) == "sum(*)"

    assert len(pprint_task((sum, list(keys) * 100), keys)) < 100
    assert pprint_task((sum, list(keys) * 100), keys) == "sum([_, _, _, ...])"
    assert (
        pprint_task((sum, [1, 2, (sum, ["a", 4]), 5, 6] * 100), keys)
        == "sum([*, *, sum([_, *]), ...])"
    )
    assert (
        pprint_task((sum, [1, 2, (sum, ["a", (sum, [1, 2, 3])]), 5, 6]), keys)
        == "sum([*, *, sum([_, sum(...)]), ...])"
    )

    # With kwargs
    def foo(w, x, y=(), z=3):
        return w + x + sum(y) + z

    task = (apply, foo, (tuple, ["a", "b"]), (dict, [["y", ["a", "b"]], ["z", "c"]]))
    assert pprint_task(task, keys) == "foo(_, _, y=[_, _], z=_)"
    task = (apply, foo, (tuple, ["a", "b"]), (dict, [["y", ["a", 1]], ["z", 1]]))
    assert pprint_task(task, keys) == "foo(_, _, y=[_, *], z=*)"


@pytest.mark.skipif("not bokeh")
def test_profiler_plot():
    with prof:
        get(dsk, "e")
    p = prof.visualize(
        width=500,
        height=300,
        tools="hover",
        title="Not the default",
        show=False,
        save=False,
    )
    if BOKEH_VERSION().major < 3:
        assert p.plot_width == 500
        assert p.plot_height == 300
    else:
        assert p.width == 500
        assert p.height == 300
    assert len(p.tools) == 1
    assert isinstance(p.tools[0], bokeh.models.HoverTool)
    assert p.title.text == "Not the default"
    # Test empty, checking for errors
    prof.clear()
    with warnings.catch_warnings(record=True) as record:
        prof.visualize(show=False, save=False)
    assert not record


@pytest.mark.skipif("not bokeh")
@pytest.mark.skipif("not psutil")
def test_resource_profiler_plot():
    with ResourceProfiler(dt=0.01) as rprof:
        get(dsk2, "c")
    p = rprof.visualize(
        width=500,
        height=300,
        tools="hover",
        title="Not the default",
        show=False,
        save=False,
    )
    if BOKEH_VERSION().major < 3:
        assert p.plot_width == 500
        assert p.plot_height == 300
    else:
        assert p.width == 500
        assert p.height == 300
    assert len(p.tools) == 1
    assert isinstance(p.tools[0], bokeh.models.HoverTool)
    assert p.title.text == "Not the default"

    # Test with empty and one point, checking for errors
    rprof.clear()
    for results in [[], [(1.0, 0, 0)]]:
        rprof.results = results
        rprof.start_time = 0.0
        rprof.end_time = 1.0
        with warnings.catch_warnings(record=True) as record:
            p = rprof.visualize(show=False, save=False)
        assert not record
        # Check bounds are valid
        assert p.x_range.start == 0
        assert p.x_range.end == 1
        assert p.y_range.start == 0
        assert p.y_range.end == 100
        assert p.extra_y_ranges["memory"].start == 0
        assert p.extra_y_ranges["memory"].end == 100


@pytest.mark.skipif("not bokeh")
def test_cache_profiler_plot():
    with CacheProfiler(metric_name="non-standard") as cprof:
        get(dsk, "e")
    p = cprof.visualize(
        width=500,
        height=300,
        tools="hover",
        title="Not the default",
        show=False,
        save=False,
    )
    if BOKEH_VERSION().major < 3:
        assert p.plot_width == 500
        assert p.plot_height == 300
    else:
        assert p.width == 500
        assert p.height == 300
    assert len(p.tools) == 1
    assert isinstance(p.tools[0], bokeh.models.HoverTool)
    assert p.title.text == "Not the default"
    assert p.axis[1].axis_label == "Cache Size (non-standard)"
    # Test empty, checking for errors
    cprof.clear()
    with warnings.catch_warnings(record=True) as record:
        cprof.visualize(show=False, save=False)
    assert not record


@pytest.mark.skipif("not bokeh")
def test_cache_profiler_plot_with_invalid_bokeh_kwarg_raises_error():
    with CacheProfiler(metric_name="non-standard") as cprof:
        get(dsk, "e")
    with pytest.raises(AttributeError, match="foo_bar"):
        cprof.visualize(foo_bar="fake")


@pytest.mark.skipif("not bokeh")
@pytest.mark.skipif("not psutil")
def test_plot_multiple():
    from dask.diagnostics.profile_visualize import visualize

    with ResourceProfiler(dt=0.01) as rprof:
        with prof:
            get(dsk2, "c")
    p = visualize(
        [prof, rprof], label_size=50, title="Not the default", show=False, save=False
    )
    # Grid plot layouts changed in Bokeh 3.
    # See https://github.com/dask/dask/issues/9257 for more details
    if BOKEH_VERSION().major < 3:
        figures = [r[0] for r in p.children[1].children]
    else:
        figures = [r[0] for r in p.children]
    assert len(figures) == 2
    assert figures[0].title.text == "Not the default"
    assert figures[0].xaxis[0].axis_label is None
    assert figures[1].title is None
    assert figures[1].xaxis[0].axis_label == "Time (s)"
    # Test empty, checking for errors
    prof.clear()
    rprof.clear()
    visualize([prof, rprof], show=False, save=False)


@pytest.mark.skipif("not bokeh")
def test_saves_file():
    with tmpfile("html") as fn:
        with prof:
            get(dsk, "e")
        # Run just to see that it doesn't error
        prof.visualize(show=False, filename=fn)

        assert os.path.exists(fn)
        with open(fn) as f:
            assert "html" in f.read().lower()


@pytest.mark.skipif("not bokeh")
def test_saves_file_path_deprecated():
    with tmpfile("html") as fn:
        with prof:
            get(dsk, "e")
        # Run just to see that it warns, but still works.
        with pytest.warns(FutureWarning) as record:
            prof.visualize(show=False, file_path=fn)

        assert 1 <= len(record) <= 2
        assert "file_path keyword argument is deprecated" in str(record[-1].message)
        # This additional warning comes from inside `bokeh`. There's a fix upstream
        # https://github.com/bokeh/bokeh/pull/12690 so for now we just ignore it.
        if len(record) == 2:
            assert "`np.bool8` is a deprecated alias for `np.bool_`" in str(
                record[0].message
            )


@pytest.mark.skipif("not bokeh")
def test_get_colors():
    from bokeh.palettes import Blues5, Blues256, Viridis

    from dask.diagnostics.profile_visualize import get_colors

    funcs = list(range(11))
    cmap = get_colors("Blues", funcs)
    assert set(cmap) < set(Blues256)
    assert len(set(cmap)) == 11

    funcs = list(range(5))
    cmap = get_colors("Blues", funcs)
    lk = dict(zip(funcs, Blues5))
    assert cmap == [lk[i] for i in funcs]

    funcs = [0, 1, 0, 1, 0, 1]
    cmap = get_colors("BrBG", funcs)
    assert len(set(cmap)) == 2

    funcs = list(range(100))
    cmap = get_colors("Viridis", funcs)
    assert len(set(cmap)) == 100

    funcs = list(range(300))
    cmap = get_colors("Viridis", funcs)
    assert len(set(cmap)) == len(set(Viridis[256]))
