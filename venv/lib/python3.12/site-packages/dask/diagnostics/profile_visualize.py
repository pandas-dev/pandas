from __future__ import annotations

import random
import warnings
from bisect import bisect_left
from itertools import cycle
from operator import add, itemgetter

from tlz import accumulate, groupby, pluck, unique

from dask._compatibility import import_optional_dependency
from dask.core import istask
from dask.utils import import_required


def BOKEH_VERSION():
    bokeh = import_optional_dependency("bokeh")
    from packaging.version import Version

    return Version(bokeh.__version__)


_BOKEH_MISSING_MSG = "Diagnostics plots require `bokeh` to be installed"


def unquote(expr):
    if istask(expr):
        if expr[0] in (tuple, list, set):
            return expr[0](map(unquote, expr[1]))
        elif (
            expr[0] == dict
            and isinstance(expr[1], list)
            and isinstance(expr[1][0], list)
        ):
            return dict(map(unquote, expr[1]))
    return expr


def get_colors(palette, funcs):
    """Get a dict mapping funcs to colors from palette.

    Parameters
    ----------
    palette : string
        Name of the bokeh palette to use, must be a member of
        bokeh.palettes.all_palettes.
    funcs : iterable
        Iterable of function names
    """
    palettes = import_required("bokeh.palettes", _BOKEH_MISSING_MSG)

    unique_funcs = sorted(unique(funcs))
    n_funcs = len(unique_funcs)
    palette_lookup = palettes.all_palettes[palette]
    keys = sorted(palette_lookup.keys())
    index = keys[min(bisect_left(keys, n_funcs), len(keys) - 1)]
    palette = palette_lookup[index]
    # Some bokeh palettes repeat colors, we want just the unique set
    palette = list(unique(palette))
    if len(palette) > n_funcs:
        # Consistently shuffle palette - prevents just using low-range
        random.Random(42).shuffle(palette)
    color_lookup = dict(zip(unique_funcs, cycle(palette)))
    return [color_lookup[n] for n in funcs]


def visualize(
    profilers, filename="profile.html", show=True, save=None, mode=None, **kwargs
):
    """Visualize the results of profiling in a bokeh plot.

    If multiple profilers are passed in, the plots are stacked vertically.

    Parameters
    ----------
    profilers : profiler or list
        Profiler or list of profilers.
    filename : string, optional
        Name of the plot output file.
    show : boolean, optional
        If True (default), the plot is opened in a browser.
    save : boolean, optional
        If True (default when not in notebook), the plot is saved to disk.
    mode : str, optional
        Mode passed to bokeh.output_file()
    **kwargs
        Other keyword arguments, passed to bokeh.figure. These will override
        all defaults set by visualize.

    Returns
    -------
    The completed bokeh plot object.
    """
    bp = import_required("bokeh.plotting", _BOKEH_MISSING_MSG)
    from bokeh.io import state

    if "file_path" in kwargs:
        warnings.warn(
            "The file_path keyword argument is deprecated "
            "and will be removed in a future release. "
            "Please use filename instead.",
            category=FutureWarning,
            stacklevel=2,
        )
        filename = kwargs.pop("file_path")

    if save is None:
        save = not state.curstate().notebook

    if not isinstance(profilers, list):
        profilers = [profilers]
    figs = [prof._plot(**kwargs) for prof in profilers]
    # Stack the plots
    if len(figs) == 1:
        p = figs[0]
    else:
        top = figs[0]
        for f in figs[1:]:
            f.x_range = top.x_range
            f.title = None
            f.min_border_top = 20
            f.height -= 30
        for f in figs[:-1]:
            f.xaxis.axis_label = None
            f.min_border_bottom = 20
            f.height -= 30
        for f in figs:
            f.min_border_left = 75
            f.min_border_right = 75
        p = bp.gridplot([[f] for f in figs])
    if show:
        bp.show(p)
    if save:
        bp.output_file(filename, mode=mode)
        bp.save(p)
    return p


def plot_tasks(
    results, dsk, start_time, end_time, palette="Viridis", label_size=60, **kwargs
):
    """Visualize the results of profiling in a bokeh plot.

    Parameters
    ----------
    results : sequence
        Output of Profiler.results
    dsk : dict
        The dask graph being profiled.
    start_time : float
        Start time of the profile in seconds
    end_time : float
        End time of the profile in seconds
    palette : string, optional
        Name of the bokeh palette to use, must be a member of
        bokeh.palettes.all_palettes.
    label_size: int (optional)
        Maximum size of output labels in plot, defaults to 60
    **kwargs
        Other keyword arguments, passed to bokeh.figure. These will override
        all defaults set by visualize.

    Returns
    -------
    The completed bokeh plot object.
    """
    bp = import_required("bokeh.plotting", _BOKEH_MISSING_MSG)
    from bokeh.models import HoverTool

    defaults = dict(
        title="Profile Results",
        tools="hover,save,reset,xwheel_zoom,xpan",
        toolbar_location="above",
        width=800,
        height=300,
    )
    # Support plot_width and plot_height for backwards compatibility
    if "plot_width" in kwargs:
        kwargs["width"] = kwargs.pop("plot_width")
    if "plot_height" in kwargs:
        kwargs["height"] = kwargs.pop("plot_height")
    defaults.update(**kwargs)

    if results:
        keys, tasks, starts, ends, ids = zip(*results)

        id_group = groupby(itemgetter(4), results)
        timings = {
            k: [i.end_time - i.start_time for i in v] for (k, v) in id_group.items()
        }
        id_lk = {
            t[0]: n
            for (n, t) in enumerate(
                sorted(timings.items(), key=itemgetter(1), reverse=True)
            )
        }

        p = bp.figure(
            y_range=[str(i) for i in range(len(id_lk))],
            x_range=[0, end_time - start_time],
            **defaults,
        )

        data = {}
        data["width"] = width = [e - s for (s, e) in zip(starts, ends)]
        data["x"] = [w / 2 + s - start_time for (w, s) in zip(width, starts)]
        data["y"] = [id_lk[i] + 1 for i in ids]
        data["function"] = funcs = [repr(i) for i in tasks]
        data["color"] = get_colors(palette, funcs)
        data["key"] = [str(i) for i in keys]

        source = bp.ColumnDataSource(data=data)

        p.rect(
            source=source,
            x="x",
            y="y",
            height=1,
            width="width",
            color="color",
            line_color="gray",
        )
    else:
        p = bp.figure(y_range=[str(i) for i in range(8)], x_range=[0, 10], **defaults)
    p.grid.grid_line_color = None
    p.axis.axis_line_color = None
    p.axis.major_tick_line_color = None
    p.yaxis.axis_label = "Worker ID"
    p.xaxis.axis_label = "Time (s)"

    hover = p.select(HoverTool)
    hover.tooltips = """
    <div>
        <span style="font-size: 14px; font-weight: bold;">Key:</span>&nbsp;
        <span style="font-size: 10px; font-family: Monaco, monospace;">@key</span>
    </div>
    <div>
        <span style="font-size: 14px; font-weight: bold;">Task:</span>&nbsp;
        <span style="font-size: 10px; font-family: Monaco, monospace;">@function</span>
    </div>
    """
    hover.point_policy = "follow_mouse"

    return p


def plot_resources(results, start_time, end_time, palette="Viridis", **kwargs):
    """Plot resource usage in a bokeh plot.

    Parameters
    ----------
    results : sequence
        Output of ResourceProfiler.results
    start_time : float
        Start time of the profile in seconds
    end_time : float
        End time of the profile in seconds
    palette : string, optional
        Name of the bokeh palette to use, must be a member of
        bokeh.palettes.all_palettes.
    **kwargs
        Other keyword arguments, passed to bokeh.figure. These will override
        all defaults set by plot_resources.

    Returns
    -------
    The completed bokeh plot object.
    """
    bp = import_required("bokeh.plotting", _BOKEH_MISSING_MSG)
    from bokeh import palettes
    from bokeh.models import LinearAxis, Range1d

    defaults = dict(
        title="Profile Results",
        tools="save,reset,xwheel_zoom,xpan",
        toolbar_location="above",
        width=800,
        height=300,
    )
    # Support plot_width and plot_height for backwards compatibility
    if "plot_width" in kwargs:
        kwargs["width"] = kwargs.pop("plot_width")
        if BOKEH_VERSION().major >= 3:
            warnings.warn("Use width instead of plot_width with Bokeh >= 3")
    if "plot_height" in kwargs:
        kwargs["height"] = kwargs.pop("plot_height")
        if BOKEH_VERSION().major >= 3:
            warnings.warn("Use height instead of plot_height with Bokeh >= 3")

    # Drop `label_size` to match `plot_cache` and `plot_tasks` kwargs
    if "label_size" in kwargs:
        kwargs.pop("label_size")

    defaults.update(**kwargs)

    if results:
        t, mem, cpu = zip(*results)
        left = start_time
        right = end_time
        t = [i - left for i in t]
        p = bp.figure(
            y_range=fix_bounds(0, max(cpu), 100),
            x_range=fix_bounds(0, right - left, 1),
            **defaults,
        )
    else:
        t = mem = cpu = []
        p = bp.figure(y_range=(0, 100), x_range=(0, 1), **defaults)
    colors = palettes.all_palettes[palette][6]
    p.line(
        t,
        cpu,
        color=colors[0],
        line_width=4,
        legend_label="% CPU",
    )
    p.yaxis.axis_label = "% CPU"
    p.extra_y_ranges = {
        "memory": Range1d(
            *fix_bounds(min(mem) if mem else 0, max(mem) if mem else 100, 100)
        )
    }
    p.line(
        t,
        mem,
        color=colors[2],
        y_range_name="memory",
        line_width=4,
        legend_label="Memory",
    )
    p.add_layout(LinearAxis(y_range_name="memory", axis_label="Memory (MB)"), "right")
    p.xaxis.axis_label = "Time (s)"
    return p


def fix_bounds(start, end, min_span):
    """Adjust end point to ensure span of at least `min_span`"""
    return start, max(end, start + min_span)


def plot_cache(
    results,
    dsk,
    start_time,
    end_time,
    metric_name,
    palette="Viridis",
    label_size=60,
    **kwargs,
):
    """Visualize the results of profiling in a bokeh plot.

    Parameters
    ----------
    results : sequence
        Output of CacheProfiler.results
    dsk : dict
        The dask graph being profiled.
    start_time : float
        Start time of the profile in seconds
    end_time : float
        End time of the profile in seconds
    metric_name : string
        Metric used to measure cache size
    palette : string, optional
        Name of the bokeh palette to use, must be a member of
        bokeh.palettes.all_palettes.
    label_size: int (optional)
        Maximum size of output labels in plot, defaults to 60
    **kwargs
        Other keyword arguments, passed to bokeh.figure. These will override
        all defaults set by visualize.

    Returns
    -------
    The completed bokeh plot object.
    """
    bp = import_required("bokeh.plotting", _BOKEH_MISSING_MSG)
    from bokeh.models import HoverTool

    defaults = dict(
        title="Profile Results",
        tools="hover,save,reset,wheel_zoom,xpan",
        toolbar_location="above",
        width=800,
        height=300,
    )
    # Support plot_width and plot_height for backwards compatibility
    if "plot_width" in kwargs:
        kwargs["width"] = kwargs.pop("plot_width")
        if BOKEH_VERSION().major >= 3:
            warnings.warn("Use width instead of plot_width with Bokeh >= 3")
    if "plot_height" in kwargs:
        kwargs["height"] = kwargs.pop("plot_height")
        if BOKEH_VERSION().major >= 3:
            warnings.warn("Use height instead of plot_height with Bokeh >= 3")
    defaults.update(**kwargs)

    if results:
        starts, ends = list(zip(*results))[3:]
        tics = sorted(unique(starts + ends))
        groups = groupby(lambda d: repr(d), results)
        data = {}
        for k, vals in groups.items():
            cnts = dict.fromkeys(tics, 0)
            for v in vals:
                cnts[v.cache_time] += v.metric
                cnts[v.free_time] -= v.metric
            data[k] = [0] + list(accumulate(add, pluck(1, sorted(cnts.items()))))

        tics = [0] + [i - start_time for i in tics]
        p = bp.figure(x_range=[0, end_time - start_time], **defaults)

        for (key, val), color in zip(data.items(), get_colors(palette, data.keys())):
            p.line(
                "x",
                "y",
                line_color=color,
                line_width=3,
                source=bp.ColumnDataSource(
                    {"x": tics, "y": val, "label": [key for i in val]}
                ),
            )

    else:
        p = bp.figure(y_range=[0, 10], x_range=[0, 10], **defaults)
    p.yaxis.axis_label = f"Cache Size ({metric_name})"
    p.xaxis.axis_label = "Time (s)"

    hover = p.select(HoverTool)
    hover.tooltips = """
    <div>
        <span style="font-size: 14px; font-weight: bold;">Task:</span>&nbsp;
        <span style="font-size: 10px; font-family: Monaco, monospace;">@label</span>
    </div>
    """
    return p
