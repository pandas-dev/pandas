"""
Typing support for Matplotlib

This module contains Type aliases which are useful for Matplotlib and potentially
downstream libraries.

.. warning::
    **Provisional status of typing**

    The ``typing`` module and type stub files are considered provisional and may change
    at any time without a deprecation period.
"""
from collections.abc import Hashable, Sequence
import pathlib
from typing import Any, Literal, TypeAlias, TypeVar
from collections.abc import Callable, Mapping

from . import path
from ._enums import JoinStyle, CapStyle
from .artist import Artist
from .backend_bases import RendererBase
from .markers import MarkerStyle
from .transforms import Bbox, Transform

DataParamType: TypeAlias = Mapping[str, Any] | None
"""The type of the *data* parameter in plotting functions."""

RGBColorType: TypeAlias = tuple[float, float, float] | str
"""Any RGB color specification accepted by Matplotlib."""

RGBAColorType: TypeAlias = (
    str |  # "none" or "#RRGGBBAA"/"#RGBA" hex strings
    tuple[float, float, float, float] |
    # 2 tuple (color, alpha) representations, not infinitely recursive
    # RGBColorType includes the (str, float) tuple, even for RGBA strings
    tuple[RGBColorType, float] |
    # (4-tuple, float) is odd, but accepted as the outer float overriding A of 4-tuple
    tuple[tuple[float, float, float, float], float]
)
"""Any RGBA color specification accepted by Matplotlib."""

ColorType: TypeAlias = RGBColorType | RGBAColorType
"""Any color specification accepted by Matplotlib. See :mpltype:`color`."""

RGBColourType: TypeAlias = RGBColorType
"""Alias of `.RGBColorType`."""

RGBAColourType: TypeAlias = RGBAColorType
"""Alias of `.RGBAColorType`."""

ColourType: TypeAlias = ColorType
"""Alias of `.ColorType`."""

LineStyleType: TypeAlias = (
    Literal["-", "solid", "--", "dashed", "-.", "dashdot", ":", "dotted",
            "", "none", " ", "None"] |
    tuple[float, Sequence[float]]
)
"""
Any line style specification accepted by Matplotlib.
See :doc:`/gallery/lines_bars_and_markers/linestyles`.
"""

DrawStyleType: TypeAlias = Literal["default", "steps", "steps-pre", "steps-mid",
                                   "steps-post"]
"""See :doc:`/gallery/lines_bars_and_markers/step_demo`."""

MarkEveryType: TypeAlias = (
    None |
    int | tuple[int, int] | slice | list[int] |
    float | tuple[float, float] |
    list[bool]
)
"""See :doc:`/gallery/lines_bars_and_markers/markevery_demo`."""

MarkerType: TypeAlias = (
    path.Path | MarkerStyle | str |  # str required for "$...$" marker
    Literal[
        ".", ",", "o", "v", "^", "<", ">",
        "1", "2", "3", "4", "8", "s", "p",
        "P", "*", "h", "H", "+", "x", "X",
        "D", "d", "|", "_", "none", " ",
        0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11
    ] | list[tuple[int, int]] | tuple[int, Literal[0, 1, 2], int]
)
"""
Marker specification. See :doc:`/gallery/lines_bars_and_markers/marker_reference`.
"""

FillStyleType: TypeAlias = Literal["full", "left", "right", "bottom", "top", "none"]
"""Marker fill styles. See :doc:`/gallery/lines_bars_and_markers/marker_reference`."""

JoinStyleType: TypeAlias = JoinStyle | Literal["miter", "round", "bevel"]
"""Line join styles. See :doc:`/gallery/lines_bars_and_markers/joinstyle`."""

CapStyleType: TypeAlias = CapStyle | Literal["butt", "projecting", "round"]
"""Line cap styles. See :doc:`/gallery/lines_bars_and_markers/capstyle`."""

LogLevel: TypeAlias = Literal["NOTSET", "DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"]
"""Literal type for valid logging levels accepted by `matplotlib.set_loglevel()`."""

CoordsBaseType: TypeAlias = (
    str |
    Artist |
    Transform |
    Callable[[RendererBase], Bbox | Transform]
)
CoordsType: TypeAlias = CoordsBaseType | tuple[CoordsBaseType, CoordsBaseType]
"""Annotation coordinate systems. See :doc:`/users/explain/text/annotations`."""

RcStyleType: TypeAlias = (
    str |
    dict[str, Any] |
    pathlib.Path |
    Sequence[str | pathlib.Path | dict[str, Any]]
)
"""
Valid specifiers for styles as used in `matplotlib.style.use` and
`matplotlib.style.context`.
"""

_HT = TypeVar("_HT", bound=Hashable)
HashableList: TypeAlias = list[_HT | "HashableList[_HT]"]
"""A nested list of Hashable values."""

MouseEventType: TypeAlias = Literal[
    "button_press_event",
    "button_release_event",
    "motion_notify_event",
    "scroll_event",
    "figure_enter_event",
    "figure_leave_event",
    "axes_enter_event",
    "axes_leave_event",
]
"""Literal type for valid `.MouseEvent` names."""

KeyEventType: TypeAlias = Literal[
    "key_press_event",
    "key_release_event"
]
"""Literal type for valid `.KeyEvent` names."""

DrawEventType: TypeAlias = Literal["draw_event"]
"""Literal type for valid `.DrawEvent` names."""
PickEventType: TypeAlias = Literal["pick_event"]
"""Literal type for valid `.PickEvent` names."""
ResizeEventType: TypeAlias = Literal["resize_event"]
"""Literal type for valid `.ResizeEvent` names."""
CloseEventType: TypeAlias = Literal["close_event"]
"""Literal type for valid `.CloseEvent` names."""

EventType: TypeAlias = Literal[
    MouseEventType,
    KeyEventType,
    DrawEventType,
    PickEventType,
    ResizeEventType,
    CloseEventType,
]
"""Literal type for all valid events."""

LegendLocType: TypeAlias = (
    Literal[
        # for simplicity, we don't distinguish the between allowed positions for
        # Axes legend and figure legend. It's still better to limit the allowed
        # range to the union of both rather than to accept arbitrary strings
        "upper right", "upper left", "lower left", "lower right",
        "right", "center left", "center right", "lower center", "upper center",
        "center",
        # Axes only
        "best",
        # Figure only
        "outside upper left", "outside upper center", "outside upper right",
        "outside right upper", "outside right center", "outside right lower",
        "outside lower right", "outside lower center", "outside lower left",
        "outside left lower", "outside left center", "outside left upper",
    ] |
    tuple[float, float] |
    int
)
"""
Supported location specifiers for legends.

This is a superset of permissible entries. "best" is only applicable to Axes legends.
All the "outside ..." locations are only applicable to figure legends.
"""

RcKeyType: TypeAlias = Literal[
    "agg.path.chunksize",
    "animation.bitrate",
    "animation.codec",
    "animation.convert_args",
    "animation.convert_path",
    "animation.embed_limit",
    "animation.ffmpeg_args",
    "animation.ffmpeg_path",
    "animation.frame_format",
    "animation.html",
    "animation.writer",
    "axes.autolimit_mode",
    "axes.axisbelow",
    "axes.edgecolor",
    "axes.facecolor",
    "axes.formatter.limits",
    "axes.formatter.min_exponent",
    "axes.formatter.offset_threshold",
    "axes.formatter.use_locale",
    "axes.formatter.use_mathtext",
    "axes.formatter.useoffset",
    "axes.grid",
    "axes.grid.axis",
    "axes.grid.which",
    "axes.labelcolor",
    "axes.labelpad",
    "axes.labelsize",
    "axes.labelweight",
    "axes.linewidth",
    "axes.prop_cycle",
    "axes.spines.bottom",
    "axes.spines.left",
    "axes.spines.right",
    "axes.spines.top",
    "axes.titlecolor",
    "axes.titlelocation",
    "axes.titlepad",
    "axes.titlesize",
    "axes.titleweight",
    "axes.titley",
    "axes.unicode_minus",
    "axes.xmargin",
    "axes.ymargin",
    "axes.zmargin",
    "axes3d.automargin",
    "axes3d.depthshade",
    "axes3d.depthshade_minalpha",
    "axes3d.grid",
    "axes3d.mouserotationstyle",
    "axes3d.trackballborder",
    "axes3d.snap_rotation",
    "axes3d.trackballsize",
    "axes3d.xaxis.panecolor",
    "axes3d.yaxis.panecolor",
    "axes3d.zaxis.panecolor",
    "backend",
    "backend_fallback",
    "boxplot.bootstrap",
    "boxplot.boxprops.color",
    "boxplot.boxprops.linestyle",
    "boxplot.boxprops.linewidth",
    "boxplot.capprops.color",
    "boxplot.capprops.linestyle",
    "boxplot.capprops.linewidth",
    "boxplot.flierprops.color",
    "boxplot.flierprops.linestyle",
    "boxplot.flierprops.linewidth",
    "boxplot.flierprops.marker",
    "boxplot.flierprops.markeredgecolor",
    "boxplot.flierprops.markeredgewidth",
    "boxplot.flierprops.markerfacecolor",
    "boxplot.flierprops.markersize",
    "boxplot.meanline",
    "boxplot.meanprops.color",
    "boxplot.meanprops.linestyle",
    "boxplot.meanprops.linewidth",
    "boxplot.meanprops.marker",
    "boxplot.meanprops.markeredgecolor",
    "boxplot.meanprops.markerfacecolor",
    "boxplot.meanprops.markersize",
    "boxplot.medianprops.color",
    "boxplot.medianprops.linestyle",
    "boxplot.medianprops.linewidth",
    "boxplot.notch",
    "boxplot.patchartist",
    "boxplot.showbox",
    "boxplot.showcaps",
    "boxplot.showfliers",
    "boxplot.showmeans",
    "boxplot.vertical",
    "boxplot.whiskerprops.color",
    "boxplot.whiskerprops.linestyle",
    "boxplot.whiskerprops.linewidth",
    "boxplot.whiskers",
    "contour.algorithm",
    "contour.corner_mask",
    "contour.linewidth",
    "contour.negative_linestyle",
    "date.autoformatter.day",
    "date.autoformatter.hour",
    "date.autoformatter.microsecond",
    "date.autoformatter.minute",
    "date.autoformatter.month",
    "date.autoformatter.second",
    "date.autoformatter.year",
    "date.converter",
    "date.epoch",
    "date.interval_multiples",
    "docstring.hardcopy",
    "errorbar.capsize",
    "errorbar.capthick",
    "errorbar.elinewidth",
    "figure.autolayout",
    "figure.constrained_layout.h_pad",
    "figure.constrained_layout.hspace",
    "figure.constrained_layout.use",
    "figure.constrained_layout.w_pad",
    "figure.constrained_layout.wspace",
    "figure.dpi",
    "figure.edgecolor",
    "figure.facecolor",
    "figure.figsize",
    "figure.frameon",
    "figure.hooks",
    "figure.labelsize",
    "figure.labelweight",
    "figure.max_open_warning",
    "figure.raise_window",
    "figure.subplot.bottom",
    "figure.subplot.hspace",
    "figure.subplot.left",
    "figure.subplot.right",
    "figure.subplot.top",
    "figure.subplot.wspace",
    "figure.titlesize",
    "figure.titleweight",
    "font.cursive",
    "font.enable_last_resort",
    "font.family",
    "font.fantasy",
    "font.monospace",
    "font.sans-serif",
    "font.serif",
    "font.size",
    "font.stretch",
    "font.style",
    "font.variant",
    "font.weight",
    "grid.alpha",
    "grid.color",
    "grid.linestyle",
    "grid.linewidth",
    "grid.major.alpha",
    "grid.major.color",
    "grid.major.linestyle",
    "grid.major.linewidth",
    "grid.minor.alpha",
    "grid.minor.color",
    "grid.minor.linestyle",
    "grid.minor.linewidth",
    "hatch.color",
    "hatch.linewidth",
    "hist.bins",
    "image.aspect",
    "image.cmap",
    "image.composite_image",
    "image.interpolation",
    "image.interpolation_stage",
    "image.lut",
    "image.origin",
    "image.resample",
    "interactive",
    "keymap.back",
    "keymap.copy",
    "keymap.forward",
    "keymap.fullscreen",
    "keymap.grid",
    "keymap.grid_minor",
    "keymap.help",
    "keymap.home",
    "keymap.pan",
    "keymap.quit",
    "keymap.quit_all",
    "keymap.save",
    "keymap.xscale",
    "keymap.yscale",
    "keymap.zoom",
    "legend.borderaxespad",
    "legend.borderpad",
    "legend.columnspacing",
    "legend.edgecolor",
    "legend.facecolor",
    "legend.fancybox",
    "legend.fontsize",
    "legend.framealpha",
    "legend.frameon",
    "legend.handleheight",
    "legend.handlelength",
    "legend.handletextpad",
    "legend.labelcolor",
    "legend.labelspacing",
    "legend.linewidth",
    "legend.loc",
    "legend.markerscale",
    "legend.numpoints",
    "legend.scatterpoints",
    "legend.shadow",
    "legend.title_fontsize",
    "lines.antialiased",
    "lines.color",
    "lines.dash_capstyle",
    "lines.dash_joinstyle",
    "lines.dashdot_pattern",
    "lines.dashed_pattern",
    "lines.dotted_pattern",
    "lines.linestyle",
    "lines.linewidth",
    "lines.marker",
    "lines.markeredgecolor",
    "lines.markeredgewidth",
    "lines.markerfacecolor",
    "lines.markersize",
    "lines.scale_dashes",
    "lines.solid_capstyle",
    "lines.solid_joinstyle",
    "macosx.window_mode",
    "markers.fillstyle",
    "mathtext.bf",
    "mathtext.bfit",
    "mathtext.cal",
    "mathtext.default",
    "mathtext.fallback",
    "mathtext.fontset",
    "mathtext.it",
    "mathtext.rm",
    "mathtext.sf",
    "mathtext.tt",
    "patch.antialiased",
    "patch.edgecolor",
    "patch.facecolor",
    "patch.force_edgecolor",
    "patch.linewidth",
    "path.effects",
    "path.simplify",
    "path.simplify_threshold",
    "path.sketch",
    "path.snap",
    "pcolor.shading",
    "pcolormesh.snap",
    "pdf.compression",
    "pdf.fonttype",
    "pdf.inheritcolor",
    "pdf.use14corefonts",
    "pgf.preamble",
    "pgf.rcfonts",
    "pgf.texsystem",
    "polaraxes.grid",
    "ps.distiller.res",
    "ps.fonttype",
    "ps.papersize",
    "ps.useafm",
    "ps.usedistiller",
    "savefig.bbox",
    "savefig.directory",
    "savefig.dpi",
    "savefig.edgecolor",
    "savefig.facecolor",
    "savefig.format",
    "savefig.orientation",
    "savefig.pad_inches",
    "savefig.transparent",
    "scatter.edgecolors",
    "scatter.marker",
    "svg.fonttype",
    "svg.hashsalt",
    "svg.id",
    "svg.image_inline",
    "text.antialiased",
    "text.color",
    "text.hinting",
    "text.hinting_factor",
    "text.kerning_factor",
    "text.language",
    "text.latex.engine",
    "text.latex.preamble",
    "text.parse_math",
    "text.usetex",
    "timezone",
    "tk.window_focus",
    "toolbar",
    "webagg.address",
    "webagg.open_in_browser",
    "webagg.port",
    "webagg.port_retries",
    "xaxis.labellocation",
    "xtick.alignment",
    "xtick.bottom",
    "xtick.color",
    "xtick.direction",
    "xtick.labelbottom",
    "xtick.labelcolor",
    "xtick.labelsize",
    "xtick.labeltop",
    "xtick.major.bottom",
    "xtick.major.pad",
    "xtick.major.size",
    "xtick.major.top",
    "xtick.major.width",
    "xtick.minor.bottom",
    "xtick.minor.ndivs",
    "xtick.minor.pad",
    "xtick.minor.size",
    "xtick.minor.top",
    "xtick.minor.visible",
    "xtick.minor.width",
    "xtick.top",
    "yaxis.labellocation",
    "ytick.alignment",
    "ytick.color",
    "ytick.direction",
    "ytick.labelcolor",
    "ytick.labelleft",
    "ytick.labelright",
    "ytick.labelsize",
    "ytick.left",
    "ytick.major.left",
    "ytick.major.pad",
    "ytick.major.right",
    "ytick.major.size",
    "ytick.major.width",
    "ytick.minor.left",
    "ytick.minor.ndivs",
    "ytick.minor.pad",
    "ytick.minor.right",
    "ytick.minor.size",
    "ytick.minor.visible",
    "ytick.minor.width",
    "ytick.right",
]
"""Valid specifiers for keys in `matplotlib.rcParams` and `matplotlib.rc_context`."""

RcGroupKeyType: TypeAlias = Literal[
    "agg",
    "agg.path",
    "animation",
    "axes",
    "axes.formatter",
    "axes.grid",
    "axes.spines",
    "axes3d",
    "axes3d.xaxis",
    "axes3d.yaxis",
    "axes3d.zaxis",
    "boxplot",
    "boxplot.boxprops",
    "boxplot.capprops",
    "boxplot.flierprops",
    "boxplot.meanprops",
    "boxplot.medianprops",
    "boxplot.whiskerprops",
    "contour",
    "date",
    "date.autoformatter",
    "docstring",
    "errorbar",
    "figure",
    "figure.constrained_layout",
    "figure.subplot",
    "font",
    "grid",
    "grid.major",
    "grid.minor",
    "hatch",
    "hist",
    "image",
    "keymap",
    "legend",
    "lines",
    "macosx",
    "markers",
    "mathtext",
    "patch",
    "path",
    "pcolor",
    "pcolormesh",
    "pdf",
    "pgf",
    "polaraxes",
    "ps",
    "ps.distiller",
    "savefig",
    "scatter",
    "svg",
    "text",
    "text.latex",
    "tk",
    "webagg",
    "xaxis",
    "xtick",
    "xtick.major",
    "xtick.minor",
    "yaxis",
    "ytick",
    "ytick.major",
    "ytick.minor",
]
"""Literal type for valid groups accepted by `matplotlib.rc`."""
