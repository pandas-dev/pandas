"""
Classes for including text in a figure.
"""

from collections.abc import Sequence
import functools
import itertools
import logging
import math
from numbers import Real
import weakref

import numpy as np

import matplotlib as mpl
from . import _api, artist, cbook, _docstring, colors as mcolors
from .artist import Artist
from .font_manager import FontProperties, fontManager, get_font
from .patches import FancyArrowPatch, FancyBboxPatch, Rectangle
from .textpath import TextPath, TextToPath  # noqa # Logically located here
from .transforms import (
    Affine2D, Bbox, BboxBase, BboxTransformTo, IdentityTransform, Transform)


_log = logging.getLogger(__name__)


@functools.lru_cache(maxsize=128)
def _rotate(theta):
    """
    Return an Affine2D object that rotates by the given angle in radians.
    """
    return Affine2D().rotate(theta)


def _rotate_point(angle, x, y):
    """
    Rotate point (x, y) by rotation angle in degrees
    """
    if angle == 0:
        return (x, y)
    angle_rad = math.radians(angle)
    cos, sin = math.cos(angle_rad), math.sin(angle_rad)
    return (cos * x - sin * y, sin * x + cos * y)


def _get_text_metrics_with_cache(renderer, text, fontprop, ismath, dpi):
    """Call ``renderer.get_text_width_height_descent``, caching the results."""

    # hit the outer cache layer and get the function to compute the metrics
    # for this renderer instance
    get_text_metrics = _get_text_metrics_function(renderer)
    # call the function to compute the metrics and return
    #
    # We pass a copy of the fontprop because FontProperties is both mutable and
    # has a `__hash__` that depends on that mutable state.  This is not ideal
    # as it means the hash of an object is not stable over time which leads to
    # very confusing behavior when used as keys in dictionaries or hashes.
    return get_text_metrics(text, fontprop.copy(), ismath, dpi)


def _get_text_metrics_function(input_renderer, _cache=weakref.WeakKeyDictionary()):
    """
    Helper function to provide a two-layered cache for font metrics


    To get the rendered size of a size of string we need to know:
      - what renderer we are using
      - the current dpi of the renderer
      - the string
      - the font properties
      - is it math text or not

    We do this as a two-layer cache with the outer layer being tied to a
    renderer instance and the inner layer handling everything else.

    The outer layer is implemented as `.WeakKeyDictionary` keyed on the
    renderer.  As long as someone else is holding a hard ref to the renderer
    we will keep the cache alive, but it will be automatically dropped when
    the renderer is garbage collected.

    The inner layer is provided by an lru_cache with a large maximum size (such
    that we expect very few cache misses in actual use cases).  As the
    dpi is mutable on the renderer, we need to explicitly include it as part of
    the cache key on the inner layer even though we do not directly use it (it is
    used in the method call on the renderer).

    This function takes a renderer and returns a function that can be used to
    get the font metrics.

    Parameters
    ----------
    input_renderer : maplotlib.backend_bases.RendererBase
        The renderer to set the cache up for.

    _cache : dict, optional
        We are using the mutable default value to attach the cache to the function.

        In principle you could pass a different dict-like to this function to inject
        a different cache, but please don't.  This is an internal function not meant to
        be reused outside of the narrow context we need it for.

        There is a possible race condition here between threads, we may need to drop the
        mutable default and switch to a threadlocal variable in the future.

    """
    if (_text_metrics := _cache.get(input_renderer, None)) is None:
        # We are going to include this in the closure we put as values in the
        # cache.  Closing over a hard-ref would create an unbreakable reference
        # cycle.
        renderer_ref = weakref.ref(input_renderer)

        # define the function locally to get a new lru_cache per renderer
        @functools.lru_cache(4096)
        # dpi is unused, but participates in cache invalidation (via the renderer).
        def _text_metrics(text, fontprop, ismath, dpi):
            # this should never happen under normal use, but this is a better error to
            # raise than an AttributeError on `None`
            if (local_renderer := renderer_ref()) is None:
                raise RuntimeError(
                    "Trying to get text metrics for a renderer that no longer exists.  "
                    "This should never happen and is evidence of a bug elsewhere."
                    )
            # do the actual method call we need and return the result
            return local_renderer.get_text_width_height_descent(text, fontprop, ismath)

        # stash the function for later use.
        _cache[input_renderer] = _text_metrics

    # return the inner function
    return _text_metrics


@_docstring.interpd
@_api.define_aliases({
    "color": ["c"],
    "fontproperties": ["font", "font_properties"],
    "fontfamily": ["family"],
    "fontname": ["name"],
    "fontsize": ["size"],
    "fontstretch": ["stretch"],
    "fontstyle": ["style"],
    "fontvariant": ["variant"],
    "fontweight": ["weight"],
    "horizontalalignment": ["ha"],
    "verticalalignment": ["va"],
    "multialignment": ["ma"],
})
class Text(Artist):
    """Handle storing and drawing of text in window or data coordinates."""

    zorder = 3
    _charsize_cache = dict()

    def __repr__(self):
        return f"Text({self._x}, {self._y}, {self._text!r})"

    def __init__(self,
                 x=0, y=0, text='', *,
                 color=None,           # defaults to rc params
                 verticalalignment='baseline',
                 horizontalalignment='left',
                 multialignment=None,
                 fontproperties=None,  # defaults to FontProperties()
                 rotation=None,
                 linespacing=None,
                 rotation_mode=None,
                 usetex=None,          # defaults to rcParams['text.usetex']
                 wrap=False,
                 transform_rotates_text=False,
                 parse_math=None,    # defaults to rcParams['text.parse_math']
                 antialiased=None,  # defaults to rcParams['text.antialiased']
                 **kwargs
                 ):
        """
        Create a `.Text` instance at *x*, *y* with string *text*.

        The text is aligned relative to the anchor point (*x*, *y*) according
        to ``horizontalalignment`` (default: 'left') and ``verticalalignment``
        (default: 'baseline'). See also
        :doc:`/gallery/text_labels_and_annotations/text_alignment`.

        While Text accepts the 'label' keyword argument, by default it is not
        added to the handles of a legend.

        Valid keyword arguments are:

        %(Text:kwdoc)s
        """
        super().__init__()
        self._x, self._y = x, y
        self._text = ''
        self._features = None
        self.set_language(None)
        self._reset_visual_defaults(
            text=text,
            color=color,
            fontproperties=fontproperties,
            usetex=usetex,
            parse_math=parse_math,
            wrap=wrap,
            verticalalignment=verticalalignment,
            horizontalalignment=horizontalalignment,
            multialignment=multialignment,
            rotation=rotation,
            transform_rotates_text=transform_rotates_text,
            linespacing=linespacing,
            rotation_mode=rotation_mode,
            antialiased=antialiased
        )
        self.update(kwargs)

    def _reset_visual_defaults(
        self,
        text='',
        color=None,
        fontproperties=None,
        usetex=None,
        parse_math=None,
        wrap=False,
        verticalalignment='baseline',
        horizontalalignment='left',
        multialignment=None,
        rotation=None,
        transform_rotates_text=False,
        linespacing=None,
        rotation_mode=None,
        antialiased=None
    ):
        self.set_text(text)
        self.set_color(mpl._val_or_rc(color, "text.color"))
        self.set_fontproperties(fontproperties)
        self.set_usetex(usetex)
        self.set_parse_math(mpl._val_or_rc(parse_math, 'text.parse_math'))
        self.set_wrap(wrap)
        self.set_verticalalignment(verticalalignment)
        self.set_horizontalalignment(horizontalalignment)
        self._multialignment = multialignment
        self.set_rotation(rotation)
        self._transform_rotates_text = transform_rotates_text
        self._bbox_patch = None  # a FancyBboxPatch instance
        self._renderer = None
        if linespacing is None:
            linespacing = 'normal'  # Maybe use rcParam later.
        self.set_linespacing(linespacing)
        self.set_rotation_mode(rotation_mode)
        self.set_antialiased(mpl._val_or_rc(antialiased, 'text.antialiased'))

    def update(self, kwargs):
        # docstring inherited
        ret = []
        kwargs = cbook.normalize_kwargs(kwargs, Text)
        sentinel = object()  # bbox can be None, so use another sentinel.
        # Update fontproperties first, as it has lowest priority.
        fontproperties = kwargs.pop("fontproperties", sentinel)
        if fontproperties is not sentinel:
            ret.append(self.set_fontproperties(fontproperties))
        # Update bbox last, as it depends on font properties.
        bbox = kwargs.pop("bbox", sentinel)
        ret.extend(super().update(kwargs))
        if bbox is not sentinel:
            ret.append(self.set_bbox(bbox))
        return ret

    def __getstate__(self):
        d = super().__getstate__()
        # remove the cached _renderer (if it exists)
        d['_renderer'] = None
        return d

    def contains(self, mouseevent):
        """
        Return whether the mouse event occurred inside the axis-aligned
        bounding-box of the text.
        """
        if (self._different_canvas(mouseevent) or not self.get_visible()
                or self._renderer is None):
            return False, {}
        # Explicitly use Text.get_window_extent(self) and not
        # self.get_window_extent() so that Annotation.contains does not
        # accidentally cover the entire annotation bounding box.
        bbox = Text.get_window_extent(self)
        inside = (bbox.x0 <= mouseevent.x <= bbox.x1
                  and bbox.y0 <= mouseevent.y <= bbox.y1)
        cattr = {}
        # if the text has a surrounding patch, also check containment for it,
        # and merge the results with the results for the text.
        if self._bbox_patch:
            patch_inside, patch_cattr = self._bbox_patch.contains(mouseevent)
            inside = inside or patch_inside
            cattr["bbox_patch"] = patch_cattr
        return inside, cattr

    def _get_xy_display(self):
        """
        Get the (possibly unit converted) transformed x, y in display coords.
        """
        x, y = self.get_unitless_position()
        return self.get_transform().transform((x, y))

    def _get_multialignment(self):
        if self._multialignment is not None:
            return self._multialignment
        else:
            return self._horizontalalignment

    def _char_index_at(self, x):
        """
        Calculate the index closest to the coordinate x in display space.

        The position of text[index] is assumed to be the sum of the widths
        of all preceding characters text[:index].

        This works only on single line texts.
        """
        if not self._text:
            return 0

        text = self._text

        fontproperties = str(self._fontproperties)
        if fontproperties not in Text._charsize_cache:
            Text._charsize_cache[fontproperties] = dict()

        charsize_cache = Text._charsize_cache[fontproperties]
        for char in set(text):
            if char not in charsize_cache:
                self.set_text(char)
                bb = self.get_window_extent()
                charsize_cache[char] = bb.x1 - bb.x0

        self.set_text(text)
        bb = self.get_window_extent()

        size_accum = np.cumsum([0] + [charsize_cache[x] for x in text])
        std_x = x - bb.x0
        return (np.abs(size_accum - std_x)).argmin()

    def get_rotation(self):
        """Return the text angle in degrees in the range [0, 360)."""
        if self.get_transform_rotates_text():
            return self.get_transform().transform_angles(
                [self._rotation], [self.get_unitless_position()]).item(0) % 360
        else:
            return self._rotation

    def get_transform_rotates_text(self):
        """
        Return whether rotations of the transform affect the text direction.
        """
        return self._transform_rotates_text

    def set_rotation_mode(self, m):
        """
        Set text rotation mode.

        Parameters
        ----------
        m : {None, 'default', 'anchor', 'xtick', 'ytick'}
            If ``"default"``, the text will be first rotated, then aligned according
            to their horizontal and vertical alignments.  If ``"anchor"``, then
            alignment occurs before rotation. "xtick" and "ytick" adjust the
            horizontal/vertical alignment so that the text is visually pointing
            towards its anchor point. This is primarily used for rotated tick
            labels and positions them nicely towards their ticks. Passing
            ``None`` will set the rotation mode to ``"default"``.
        """
        if m is None:
            m = "default"
        else:
            _api.check_in_list(("anchor", "default", "xtick", "ytick"), rotation_mode=m)
        self._rotation_mode = m
        self.stale = True

    def get_rotation_mode(self):
        """Return the text rotation mode."""
        return self._rotation_mode

    def set_antialiased(self, antialiased):
        """
        Set whether to use antialiased rendering.

        Parameters
        ----------
        antialiased : bool

        Notes
        -----
        Antialiasing will be determined by :rc:`text.antialiased`
        and the parameter *antialiased* will have no effect if the text contains
        math expressions.
        """
        self._antialiased = antialiased
        self.stale = True

    def get_antialiased(self):
        """Return whether antialiased rendering is used."""
        return self._antialiased

    def update_from(self, other):
        # docstring inherited
        super().update_from(other)
        self._color = other._color
        self._multialignment = other._multialignment
        self._verticalalignment = other._verticalalignment
        self._horizontalalignment = other._horizontalalignment
        self._fontproperties = other._fontproperties.copy()
        self._usetex = other._usetex
        self._rotation = other._rotation
        self._transform_rotates_text = other._transform_rotates_text
        self._picker = other._picker
        self._linespacing = other._linespacing
        self._antialiased = other._antialiased
        self.stale = True

    def _get_layout(self, renderer):
        """
        Return

        - the rotated, axis-aligned text bbox;
        - a list of ``(line, (width, ascent, descent), xy)`` tuples for each line;
        - a ``(xy, (width, height))` pair of the lower-left corner and size of the
          rotated, *text-aligned* text box (i.e. describing how to draw the
          text-surrounding box).
        """
        thisx, thisy = 0.0, 0.0
        lines = self._get_wrapped_text().split("\n")  # Ensures lines is not empty.

        # Reminder:  The ascent (a) goes from the baseline to the top and the
        # descent (d) from the baseline to the bottom; both are (typically)
        # nonnegative.  The height h is the sum, h = a + d.
        wads = []  # (width, ascents, descents)
        xs = []
        ys = []

        min_ascent = min_descent = line_gap = None
        dpi = self.get_figure(root=True).dpi
        # Determine full vertical extent of font, including ascenders and descenders:
        if not self.get_usetex():
            if hasattr(renderer, '_get_font_height_metrics'):
                # TODO: This is a temporary internal method call (for _backend_pdf_ps to
                # support AFM files) until we design a proper API for the backends.
                min_ascent, min_descent, line_gap = renderer._get_font_height_metrics(
                    self._fontproperties)
            if min_ascent is None:
                font = get_font(fontManager._find_fonts_by_props(self._fontproperties))
                possible = [
                    ('OS/2', 'sTypoLineGap', 'sTypoAscender', 'sTypoDescender'),
                    ('hhea', 'lineGap', 'ascent', 'descent')
                ]
                for table_name, linegap_key, ascent_key, descent_key in possible:
                    table = font.get_sfnt_table(table_name)
                    if table is None:
                        continue
                    # Rescale to font size/DPI if the metrics were available.
                    fontsize = self._fontproperties.get_size_in_points()
                    units_per_em = font.get_sfnt_table('head')['unitsPerEm']
                    scale = 1 / units_per_em * fontsize * dpi / 72
                    line_gap = table[linegap_key] * scale
                    min_ascent = table[ascent_key] * scale
                    min_descent = -table[descent_key] * scale
                    break
        if None in (min_ascent, min_descent):
            # Fallback to font measurement.
            _, h, min_descent = _get_text_metrics_with_cache(
                renderer, "lp", self._fontproperties,
                ismath="TeX" if self.get_usetex() else False,
                dpi=dpi)
            min_ascent = h - min_descent
            line_gap = 0

        # Don't increase text height too much if it's not multiple lines.
        if len(lines) == 1:
            line_gap = 0

        for line in lines:
            clean_line, ismath = self._preprocess_math(line)
            if clean_line:
                w, h, d = _get_text_metrics_with_cache(
                    renderer, clean_line, self._fontproperties,
                    ismath=ismath, dpi=self.get_figure(root=True).dpi)
            else:
                w = h = d = 0

            a = h - d

            if self.get_usetex() or self._linespacing == 'normal':
                # To ensure good linespacing, pretend that the ascent / descent of all
                # lines is at least as large as the measured sizes.
                a = max(a, min_ascent) + line_gap / 2
                d = max(d, min_descent) + line_gap / 2
            else:
                # If using a fixed line spacing, then every line's spacing will be
                # determined by the font metrics of the first available font.
                line_height = self._linespacing * (min_ascent + min_descent)
                leading = line_height - (a + d)
                a += leading / 2
                d += leading / 2

            # Metrics of the last line that are needed later:
            baseline = a - thisy

            thisy -= a

            wads.append((w, a, d))
            xs.append(thisx)  # == 0.
            ys.append(thisy)

            thisy -= d

        # Metrics of the last line that are needed later:
        descent = d

        # Bounding box definition:
        ws = [w for w, a, d in wads]
        width = max(ws)
        xmin = 0
        xmax = width
        ymax = 0
        ymin = ys[-1] - descent  # baseline of last line minus its descent

        # now offset the individual text lines within the box
        malign = self._get_multialignment()
        if malign == 'left':
            offset_layout = [(x, y) for x, y in zip(xs, ys)]
        elif malign == 'center':
            offset_layout = [(x + width / 2 - w / 2, y)
                             for x, y, w in zip(xs, ys, ws)]
        elif malign == 'right':
            offset_layout = [(x + width - w, y)
                             for x, y, w in zip(xs, ys, ws)]

        # the corners of the unrotated bounding box
        corners_horiz = [(xmin, ymin), (xmin, ymax), (xmax, ymax), (xmax, ymin)]
        size_horiz = (xmax - xmin, ymax - ymin)
        # now rotate the bbox
        angle = self.get_rotation()
        rotate = functools.partial(_rotate_point, angle)
        corners_rotated = [rotate(x, y) for x, y in corners_horiz]

        # compute the bounds of the rotated box
        xs, ys = zip(*corners_rotated)
        xmin, xmax = min(xs), max(xs)
        ymin, ymax = min(ys), max(ys)
        width_rot = xmax - xmin
        height_rot = ymax - ymin

        # Now move the box to the target position offset the display
        # bbox by alignment
        halign = self._horizontalalignment
        valign = self._verticalalignment

        rotation_mode = self.get_rotation_mode()
        if rotation_mode != "anchor":
            if rotation_mode == 'xtick':
                halign = self._ha_for_angle(angle)
            elif rotation_mode == 'ytick':
                valign = self._va_for_angle(angle)
            # compute the text location in display coords and the offsets
            # necessary to align the bbox with that location
            offsetx = (
                xmin if halign == "left" else
                xmax if halign == "right" else
                (xmin + xmax) / 2  # halign == "center"
            )
            offsety = (
                ymin if valign == "bottom" else
                ymax if valign == "top" else
                (ymin + ymax) / 2 if valign == "center" else
                ymin + descent if valign == "baseline" else
                ymin + height_rot - baseline / 2  # valign == "center_baseline"
            )
        else:
            xmin1, ymin1 = corners_horiz[0]
            xmax1, ymax1 = corners_horiz[2]
            offsetx = (
                xmin1 if halign == "left" else
                xmax1 if halign == "right" else
                (xmin1 + xmax1) / 2  # halign == "center"
            )
            offsety = (
                ymin1 if valign == "bottom" else
                ymax1 if valign == "top" else
                (ymin1 + ymax1) / 2 if valign == "center" else
                ymax1 - baseline if valign == "baseline" else
                ymax1 - baseline / 2  # valign == "center_baseline"
            )
            offsetx, offsety = rotate(offsetx, offsety)

        xmin -= offsetx
        ymin -= offsety

        bbox_rot = Bbox.from_bounds(xmin, ymin, width_rot, height_rot)

        # now rotate the positions around the first (x, y) position
        xys = [(x - offsetx, y - offsety)
               for x, y in itertools.starmap(rotate, offset_layout)]
        x, y = corners_rotated[0]
        xy_corner = (x - offsetx, y - offsety)

        return bbox_rot, list(zip(lines, wads, xys)), (xy_corner, size_horiz)

    def set_bbox(self, rectprops):
        """
        Draw a box behind/around the text.

        This can be used to set a background and/or a frame around the text.
        It's realized through a `.FancyBboxPatch` behind the text (see also
        `.Text.get_bbox_patch`). The bbox patch is None by default and only
        created when needed.

        Parameters
        ----------
        rectprops : dict with properties for `.FancyBboxPatch` or None
             The default boxstyle is 'square'. The mutation
             scale of the `.patches.FancyBboxPatch` is set to the fontsize.

             Pass ``None`` to remove the bbox patch completely.

        Examples
        --------
        ::

            t.set_bbox(dict(facecolor='red', alpha=0.5))
        """

        if rectprops is not None:
            props = rectprops.copy()
            boxstyle = props.pop("boxstyle", None)
            pad = props.pop("pad", None)
            if boxstyle is None:
                boxstyle = "square"
                if pad is None:
                    pad = 4  # points
                pad /= self.get_size()  # to fraction of font size
            else:
                if pad is None:
                    pad = 0.3
            # boxstyle could be a callable or a string
            if isinstance(boxstyle, str) and "pad" not in boxstyle:
                boxstyle += ",pad=%0.2f" % pad
            self._bbox_patch = FancyBboxPatch(
                (0, 0), 1, 1,
                boxstyle=boxstyle, transform=IdentityTransform(), **props)
        else:
            self._bbox_patch = None

        self._update_clip_properties()

    def get_bbox_patch(self):
        """
        Return the bbox Patch, or None if the `.patches.FancyBboxPatch`
        is not made.

        For more details see `.Text.set_bbox`.
        """
        return self._bbox_patch

    def update_bbox_position_size(self, renderer):
        """
        Update the location and the size of the bbox.

        This method should be used when the position and size of the bbox needs
        to be updated before actually drawing the bbox.
        """
        if self._bbox_patch:
            # don't use self.get_unitless_position here, which refers to text
            # position in Text:
            posx = float(self.convert_xunits(self._x))
            posy = float(self.convert_yunits(self._y))
            posx, posy = self.get_transform().transform((posx, posy))

            _, _, ((x_box, y_box), (w_box, h_box)) = self._get_layout(renderer)
            self._bbox_patch.set_bounds(0., 0., w_box, h_box)
            self._bbox_patch.set_transform(
                Affine2D()
                .rotate_deg(self.get_rotation())
                .translate(posx + x_box, posy + y_box))
            fontsize_in_pixel = renderer.points_to_pixels(self.get_size())
            self._bbox_patch.set_mutation_scale(fontsize_in_pixel)

    def _update_clip_properties(self):
        if self._bbox_patch:
            clipprops = dict(clip_box=self.clipbox,
                             clip_path=self._clippath,
                             clip_on=self._clipon)
            self._bbox_patch.update(clipprops)

    def set_clip_box(self, clipbox):
        # docstring inherited.
        super().set_clip_box(clipbox)
        self._update_clip_properties()

    def set_clip_path(self, path, transform=None):
        # docstring inherited.
        super().set_clip_path(path, transform)
        self._update_clip_properties()

    def set_clip_on(self, b):
        # docstring inherited.
        super().set_clip_on(b)
        self._update_clip_properties()

    def get_wrap(self):
        """Return whether the text can be wrapped."""
        return self._wrap

    def set_wrap(self, wrap):
        """
        Set whether the text can be wrapped.

        Wrapping makes sure the text is confined to the (sub)figure box. It
        does not take into account any other artists.

        Parameters
        ----------
        wrap : bool

        Notes
        -----
        Wrapping does not work together with
        ``savefig(..., bbox_inches='tight')`` (which is also used internally
        by ``%matplotlib inline`` in IPython/Jupyter). The 'tight' setting
        rescales the canvas to accommodate all content and happens before
        wrapping.
        """
        self._wrap = wrap

    def _get_wrap_line_width(self):
        """
        Return the maximum line width for wrapping text based on the current
        orientation.
        """
        x0, y0 = self.get_transform().transform(self.get_position())
        figure_box = self.get_figure().get_window_extent()

        # Calculate available width based on text alignment
        alignment = self.get_horizontalalignment()
        self.set_rotation_mode('anchor')
        angle = self.get_rotation()

        left = self._get_dist_to_box(angle, x0, y0, figure_box)
        right = self._get_dist_to_box(
            (180 + angle) % 360, x0, y0, figure_box)

        if alignment == 'left':
            line_width = left
        elif alignment == 'right':
            line_width = right
        else:
            line_width = 2 * min(left, right)

        return line_width

    def _get_dist_to_box(self, rotation, x0, y0, figure_box):
        """
        Return the distance from the given points to the boundaries of a
        rotated box, in pixels.
        """
        if rotation > 270:
            quad = rotation - 270
            h1 = (y0 - figure_box.y0) / math.cos(math.radians(quad))
            h2 = (figure_box.x1 - x0) / math.cos(math.radians(90 - quad))
        elif rotation > 180:
            quad = rotation - 180
            h1 = (x0 - figure_box.x0) / math.cos(math.radians(quad))
            h2 = (y0 - figure_box.y0) / math.cos(math.radians(90 - quad))
        elif rotation > 90:
            quad = rotation - 90
            h1 = (figure_box.y1 - y0) / math.cos(math.radians(quad))
            h2 = (x0 - figure_box.x0) / math.cos(math.radians(90 - quad))
        else:
            h1 = (figure_box.x1 - x0) / math.cos(math.radians(rotation))
            h2 = (figure_box.y1 - y0) / math.cos(math.radians(90 - rotation))

        return min(h1, h2)

    def _get_rendered_text_width(self, text):
        """
        Return the width of a given text string, in pixels.
        """

        w, h, d = _get_text_metrics_with_cache(
            self._renderer, text, self.get_fontproperties(),
            cbook.is_math_text(text),
            self.get_figure(root=True).dpi)
        return math.ceil(w)

    def _get_wrapped_text(self):
        """
        Return a copy of the text string with new lines added so that the text
        is wrapped relative to the parent figure (if `get_wrap` is True).
        """
        if not self.get_wrap():
            return self.get_text()

        # Not fit to handle breaking up latex syntax correctly, so
        # ignore latex for now.
        if self.get_usetex():
            return self.get_text()

        # Build the line incrementally, for a more accurate measure of length
        line_width = self._get_wrap_line_width()
        wrapped_lines = []

        # New lines in the user's text force a split
        unwrapped_lines = self.get_text().split('\n')

        # Now wrap each individual unwrapped line
        for unwrapped_line in unwrapped_lines:

            sub_words = unwrapped_line.split(' ')
            # Remove items from sub_words as we go, so stop when empty
            while len(sub_words) > 0:
                if len(sub_words) == 1:
                    # Only one word, so just add it to the end
                    wrapped_lines.append(sub_words.pop(0))
                    continue

                for i in range(2, len(sub_words) + 1):
                    # Get width of all words up to and including here
                    line = ' '.join(sub_words[:i])
                    current_width = self._get_rendered_text_width(line)

                    # If all these words are too wide, append all not including
                    # last word
                    if current_width > line_width:
                        wrapped_lines.append(' '.join(sub_words[:i - 1]))
                        sub_words = sub_words[i - 1:]
                        break

                    # Otherwise if all words fit in the width, append them all
                    elif i == len(sub_words):
                        wrapped_lines.append(' '.join(sub_words[:i]))
                        sub_words = []
                        break

        return '\n'.join(wrapped_lines)

    @artist.allow_rasterization
    def draw(self, renderer):
        # docstring inherited

        if renderer is not None:
            self._renderer = renderer
        if not self.get_visible():
            return
        if self.get_text() == '':
            return

        renderer.open_group('text', self.get_gid())

        bbox, info, _ = self._get_layout(renderer)
        trans = self.get_transform()

        # don't use self.get_position here, which refers to text
        # position in Text:
        x, y = self._x, self._y
        if np.ma.is_masked(x):
            x = np.nan
        if np.ma.is_masked(y):
            y = np.nan
        posx = float(self.convert_xunits(x))
        posy = float(self.convert_yunits(y))
        posx, posy = trans.transform((posx, posy))
        if np.isnan(posx) or np.isnan(posy):
            return  # don't throw a warning here
        if not np.isfinite(posx) or not np.isfinite(posy):
            _log.warning("posx and posy should be finite values")
            return
        canvasw, canvash = renderer.get_canvas_width_height()

        # Update the location and size of the bbox
        # (`.patches.FancyBboxPatch`), and draw it.
        if self._bbox_patch:
            self.update_bbox_position_size(renderer)
            self._bbox_patch.draw(renderer)

        gc = renderer.new_gc()
        gc.set_foreground(mcolors.to_rgba(self.get_color()), isRGBA=True)
        gc.set_alpha(self.get_alpha())
        gc.set_url(self._url)
        gc.set_antialiased(self._antialiased)
        gc.set_snap(self.get_snap())
        self._set_gc_clip(gc)

        angle = self.get_rotation()

        for line, wad, (x, y) in info:

            mtext = self if len(info) == 1 else None
            x = x + posx
            y = y + posy
            if renderer.flipy():
                y = canvash - y
            clean_line, ismath = self._preprocess_math(line)

            if self.get_path_effects():
                from matplotlib.patheffects import PathEffectRenderer
                textrenderer = PathEffectRenderer(self.get_path_effects(), renderer)
            else:
                textrenderer = renderer

            if self.get_usetex():
                textrenderer.draw_tex(gc, x, y, clean_line,
                                      self._fontproperties, angle,
                                      mtext=mtext)
            else:
                textrenderer.draw_text(gc, x, y, clean_line,
                                       self._fontproperties, angle,
                                       ismath=ismath, mtext=mtext)

        gc.restore()
        renderer.close_group('text')
        self.stale = False

    def get_color(self):
        """Return the color of the text."""
        return self._color

    def get_fontproperties(self):
        """Return the `.font_manager.FontProperties`."""
        return self._fontproperties

    def get_fontfamily(self):
        """
        Return the list of font families used for font lookup.

        See Also
        --------
        .font_manager.FontProperties.get_family
        """
        return self._fontproperties.get_family()

    def get_fontfeatures(self):
        """
        Return a tuple of font feature tags to enable.
        """
        return self._features

    def get_fontname(self):
        """
        Return the font name as a string.

        See Also
        --------
        .font_manager.FontProperties.get_name
        """
        return self._fontproperties.get_name()

    def get_fontstyle(self):
        """
        Return the font style as a string.

        See Also
        --------
        .font_manager.FontProperties.get_style
        """
        return self._fontproperties.get_style()

    def get_fontsize(self):
        """
        Return the font size as an integer.

        See Also
        --------
        .font_manager.FontProperties.get_size_in_points
        """
        return self._fontproperties.get_size_in_points()

    def get_fontvariant(self):
        """
        Return the font variant as a string.

        See Also
        --------
        .font_manager.FontProperties.get_variant
        """
        return self._fontproperties.get_variant()

    def get_fontweight(self):
        """
        Return the font weight as a string or a number.

        See Also
        --------
        .font_manager.FontProperties.get_weight
        """
        return self._fontproperties.get_weight()

    def get_stretch(self):
        """
        Return the font stretch as a string or a number.

        See Also
        --------
        .font_manager.FontProperties.get_stretch
        """
        return self._fontproperties.get_stretch()

    def get_horizontalalignment(self):
        """
        Return the horizontal alignment as a string.  Will be one of
        'left', 'center' or 'right'.
        """
        return self._horizontalalignment

    def get_unitless_position(self):
        """Return the (x, y) unitless position of the text."""
        # This will get the position with all unit information stripped away.
        # This is here for convenience since it is done in several locations.
        x = float(self.convert_xunits(self._x))
        y = float(self.convert_yunits(self._y))
        return x, y

    def get_position(self):
        """Return the (x, y) position of the text."""
        # This should return the same data (possible unitized) as was
        # specified with 'set_x' and 'set_y'.
        return self._x, self._y

    def get_text(self):
        """Return the text string."""
        return self._text

    def get_verticalalignment(self):
        """
        Return the vertical alignment as a string.  Will be one of
        'top', 'center', 'bottom', 'baseline' or 'center_baseline'.
        """
        return self._verticalalignment

    def get_window_extent(self, renderer=None, dpi=None):
        """
        Return the `.Bbox` bounding the text, in display units.

        In addition to being used internally, this is useful for specifying
        clickable regions in a png file on a web page.

        Parameters
        ----------
        renderer : Renderer, optional
            A renderer is needed to compute the bounding box.  If the artist
            has already been drawn, the renderer is cached; thus, it is only
            necessary to pass this argument when calling `get_window_extent`
            before the first draw.  In practice, it is usually easier to
            trigger a draw first, e.g. by calling
            `~.Figure.draw_without_rendering` or ``plt.show()``.

        dpi : float, optional
            The dpi value for computing the bbox, defaults to
            ``self.get_figure(root=True).dpi`` (*not* the renderer dpi); should be set
            e.g. if to match regions with a figure saved with a custom dpi value.
        """
        if not self.get_visible():
            return Bbox.unit()

        fig = self.get_figure(root=True)
        if dpi is None:
            dpi = fig.dpi
        if self.get_text() == '':
            with cbook._setattr_cm(fig, dpi=dpi):
                tx, ty = self._get_xy_display()
                return Bbox.from_bounds(tx, ty, 0, 0)

        if renderer is not None:
            self._renderer = renderer
        if self._renderer is None:
            self._renderer = fig._get_renderer()
        if self._renderer is None:
            raise RuntimeError(
                "Cannot get window extent of text w/o renderer. You likely "
                "want to call 'figure.draw_without_rendering()' first.")

        with cbook._setattr_cm(fig, dpi=dpi):
            bbox, _, _ = self._get_layout(self._renderer)
            x, y = self.get_unitless_position()
            x, y = self.get_transform().transform((x, y))
            bbox = bbox.translated(x, y)
            return bbox

    def get_tightbbox(self, renderer=None):
        if not self.get_visible() or self.get_text() == "":
            return Bbox.null()
        # Exclude text at data coordinates outside the valid domain of the axes
        # scales (e.g., negative coordinates with a log scale).
        if (self.axes
                and self.get_transform() == self.axes.transData
                and not self.axes._point_in_data_domain(*self.get_unitless_position())):
            return Bbox.null()
        return super().get_tightbbox(renderer)

    def set_backgroundcolor(self, color):
        """
        Set the background color of the text.

        This is realized through the bbox (see `.set_bbox`),
        creating the bbox patch if needed.

        Parameters
        ----------
        color : :mpltype:`color`

        See Also
        --------
        .set_bbox : To change the position of the bounding box
        """
        if self._bbox_patch is None:
            self.set_bbox(dict(facecolor=color, edgecolor=color))
        else:
            self._bbox_patch.update(dict(facecolor=color))

        self._update_clip_properties()
        self.stale = True

    def set_color(self, color):
        """
        Set the foreground color of the text

        Parameters
        ----------
        color : :mpltype:`color`
        """
        # "auto" is only supported by axisartist, but we can just let it error
        # out at draw time for simplicity.
        if not cbook._str_equal(color, "auto"):
            mpl.colors._check_color_like(color=color)
        self._color = color
        self.stale = True

    def set_horizontalalignment(self, align):
        """
        Set the horizontal alignment relative to the anchor point.

        See also :doc:`/gallery/text_labels_and_annotations/text_alignment`.

        Parameters
        ----------
        align : {'left', 'center', 'right'}
        """
        _api.check_in_list(['center', 'right', 'left'], align=align)
        self._horizontalalignment = align
        self.stale = True

    def set_multialignment(self, align):
        """
        Set the text alignment for multiline texts.

        The layout of the bounding box of all the lines is determined by the
        horizontalalignment and verticalalignment properties. This property
        controls the alignment of the text lines within that box.

        Parameters
        ----------
        align : {'left', 'right', 'center'}
        """
        _api.check_in_list(['center', 'right', 'left'], align=align)
        self._multialignment = align
        self.stale = True

    def set_linespacing(self, spacing):
        """
        Set the line spacing.

        Parameters
        ----------
        spacing : 'normal' or float, default: 'normal'
            If 'normal', then the line spacing is automatically determined by font
            metrics for each line individually.

            If a float, then line spacing will be fixed to this multiple of the font
            size for every line.
        """
        if not cbook._str_equal(spacing, 'normal'):
            _api.check_isinstance(Real, spacing=spacing)
        self._linespacing = spacing
        self.stale = True

    def get_linespacing(self):
        """Get the line spacing."""
        return self._linespacing

    def set_fontfamily(self, fontname):
        """
        Set the font family.  Can be either a single string, or a list of
        strings in decreasing priority.  Each string may be either a real font
        name or a generic font class name.  If the latter, the specific font
        names will be looked up in the corresponding rcParams.

        If a `Text` instance is constructed with ``fontfamily=None``, then the
        font is set to :rc:`font.family`, and the
        same is done when `set_fontfamily()` is called on an existing
        `Text` instance.

        Parameters
        ----------
        fontname : {FONTNAME, 'serif', 'sans-serif', 'cursive', 'fantasy', \
'monospace'}

        See Also
        --------
        .font_manager.FontProperties.set_family
        """
        self._fontproperties.set_family(fontname)
        self.stale = True

    def set_fontfeatures(self, features):
        """
        Set the feature tags to enable on the font.

        Parameters
        ----------
        features : list of str, or tuple of str, or None
            A list of feature tags to be used with the associated font. These strings
            are eventually passed to HarfBuzz, and so all `string formats supported by
            hb_feature_from_string()
            <https://harfbuzz.github.io/harfbuzz-hb-common.html#hb-feature-from-string>`__
            are supported. Note though that subranges are not explicitly supported and
            behaviour may change in the future.

            For example, if your desired font includes Stylistic Sets which enable
            various typographic alternates including one that you do not wish to use
            (e.g., Contextual Ligatures), then you can pass the following to enable one
            and not the other::

                fp.set_features([
                    'ss01',  # Use Stylistic Set 1.
                    '-clig',  # But disable Contextural Ligatures.
                ])

            Available font feature tags may be found at
            https://learn.microsoft.com/en-us/typography/opentype/spec/featurelist
        """
        _api.check_isinstance((Sequence, None), features=features)
        if features is not None:
            features = tuple(features)
        self._features = features
        self.stale = True

    def set_fontvariant(self, variant):
        """
        Set the font variant.

        Parameters
        ----------
        variant : {'normal', 'small-caps'}

        See Also
        --------
        .font_manager.FontProperties.set_variant
        """
        self._fontproperties.set_variant(variant)
        self.stale = True

    def set_fontstyle(self, fontstyle):
        """
        Set the font style.

        Parameters
        ----------
        fontstyle : {'normal', 'italic', 'oblique'}

        See Also
        --------
        .font_manager.FontProperties.set_style
        """
        self._fontproperties.set_style(fontstyle)
        self.stale = True

    def set_fontsize(self, fontsize):
        """
        Set the font size.

        Parameters
        ----------
        fontsize : float or {'xx-small', 'x-small', 'small', 'medium', \
'large', 'x-large', 'xx-large'}
            If a float, the fontsize in points. The string values denote sizes
            relative to the default font size.

        See Also
        --------
        .font_manager.FontProperties.set_size
        """
        self._fontproperties.set_size(fontsize)
        self.stale = True

    def get_math_fontfamily(self):
        """
        Return the font family name for math text rendered by Matplotlib.

        The default value is :rc:`mathtext.fontset`.

        See Also
        --------
        set_math_fontfamily
        """
        return self._fontproperties.get_math_fontfamily()

    def set_math_fontfamily(self, fontfamily):
        """
        Set the font family for math text rendered by Matplotlib.

        This does only affect Matplotlib's own math renderer. It has no effect
        when rendering with TeX (``usetex=True``).

        Parameters
        ----------
        fontfamily : str
            The name of the font family.

            Available font families are defined in the
            :ref:`default matplotlibrc file
            <customizing-with-matplotlibrc-files>`.

        See Also
        --------
        get_math_fontfamily
        """
        self._fontproperties.set_math_fontfamily(fontfamily)

    def set_fontweight(self, weight):
        """
        Set the font weight.

        Parameters
        ----------
        weight : {a numeric value in range 0-1000, 'ultralight', 'light', \
'normal', 'regular', 'book', 'medium', 'roman', 'semibold', 'demibold', \
'demi', 'bold', 'heavy', 'extra bold', 'black'}

        See Also
        --------
        .font_manager.FontProperties.set_weight
        """
        self._fontproperties.set_weight(weight)
        self.stale = True

    def set_fontstretch(self, stretch):
        """
        Set the font stretch (horizontal condensation or expansion).

        Parameters
        ----------
        stretch : {a numeric value in range 0-1000, 'ultra-condensed', \
'extra-condensed', 'condensed', 'semi-condensed', 'normal', 'semi-expanded', \
'expanded', 'extra-expanded', 'ultra-expanded'}

        See Also
        --------
        .font_manager.FontProperties.set_stretch
        """
        self._fontproperties.set_stretch(stretch)
        self.stale = True

    def set_position(self, xy):
        """
        Set the (*x*, *y*) position of the text.

        Parameters
        ----------
        xy : (float, float)
        """
        self.set_x(xy[0])
        self.set_y(xy[1])

    def set_x(self, x):
        """
        Set the *x* position of the text.

        Parameters
        ----------
        x : float
        """
        self._x = x
        self.stale = True

    def set_y(self, y):
        """
        Set the *y* position of the text.

        Parameters
        ----------
        y : float
        """
        self._y = y
        self.stale = True

    def set_rotation(self, s):
        """
        Set the rotation of the text.

        Parameters
        ----------
        s : float or {'vertical', 'horizontal'}
            The rotation angle in degrees in mathematically positive direction
            (counterclockwise). 'horizontal' equals 0, 'vertical' equals 90.
        """
        if isinstance(s, Real):
            self._rotation = float(s) % 360
        elif cbook._str_equal(s, 'horizontal') or s is None:
            self._rotation = 0.
        elif cbook._str_equal(s, 'vertical'):
            self._rotation = 90.
        else:
            raise ValueError("rotation must be 'vertical', 'horizontal' or "
                             f"a number, not {s}")
        self.stale = True

    def set_transform_rotates_text(self, t):
        """
        Whether rotations of the transform affect the text direction.

        Parameters
        ----------
        t : bool
        """
        self._transform_rotates_text = t
        self.stale = True

    def set_verticalalignment(self, align):
        """
        Set the vertical alignment relative to the anchor point.

        See also :doc:`/gallery/text_labels_and_annotations/text_alignment`.

        Parameters
        ----------
        align : {'baseline', 'bottom', 'center', 'center_baseline', 'top'}
        """
        _api.check_in_list(
            ['top', 'bottom', 'center', 'baseline', 'center_baseline'],
            align=align)
        self._verticalalignment = align
        self.stale = True

    def set_text(self, s):
        r"""
        Set the text string *s*.

        It may contain newlines (``\n``) or math in LaTeX syntax.

        Parameters
        ----------
        s : object
            Any object gets converted to its `str` representation, except for
            ``None`` which is converted to an empty string.
        """
        s = '' if s is None else str(s)
        if s != self._text:
            self._text = s
            self.stale = True

    def _preprocess_math(self, s):
        """
        Return the string *s* after mathtext preprocessing, and the kind of
        mathtext support needed.

        - If *self* is configured to use TeX, return *s* unchanged except that
          a single space gets escaped, and the flag "TeX".
        - Otherwise, if *s* is mathtext (has an even number of unescaped dollar
          signs) and ``parse_math`` is not set to False, return *s* and the
          flag True.
        - Otherwise, return *s* with dollar signs unescaped, and the flag
          False.
        """
        if self.get_usetex():
            if s == " ":
                s = r"\ "
            return s, "TeX"
        elif not self.get_parse_math():
            return s, False
        elif cbook.is_math_text(s):
            return s, True
        else:
            return s.replace(r"\$", "$"), False

    def set_fontproperties(self, fp):
        """
        Set the font properties that control the text.

        Parameters
        ----------
        fp : `.font_manager.FontProperties` or `str` or `pathlib.Path`
            If a `str`, it is interpreted as a fontconfig pattern parsed by
            `.FontProperties`.  If a `pathlib.Path`, it is interpreted as the
            absolute path to a font file.
        """
        self._fontproperties = FontProperties._from_any(fp).copy()
        self.stale = True

    @_docstring.kwarg_doc("bool, default: :rc:`text.usetex`")
    def set_usetex(self, usetex):
        """
        Parameters
        ----------
        usetex : bool or None
            Whether to render using TeX, ``None`` means to use
            :rc:`text.usetex`.
        """
        self._usetex = bool(mpl._val_or_rc(usetex, 'text.usetex'))
        self.stale = True

    def get_usetex(self):
        """Return whether this `Text` object uses TeX for rendering."""
        return self._usetex

    def set_parse_math(self, parse_math):
        """
        Override switch to disable any mathtext parsing for this `Text`.

        Parameters
        ----------
        parse_math : bool
            If False, this `Text` will never use mathtext.  If True, mathtext
            will be used if there is an even number of unescaped dollar signs.
        """
        self._parse_math = bool(parse_math)

    def get_parse_math(self):
        """Return whether mathtext parsing is considered for this `Text`."""
        return self._parse_math

    def set_fontname(self, fontname):
        """
        Alias for `set_fontfamily`.

        One-way alias only: the getter differs.

        Parameters
        ----------
        fontname : {FONTNAME, 'serif', 'sans-serif', 'cursive', 'fantasy', \
'monospace'}

        See Also
        --------
        .font_manager.FontProperties.set_family

        """
        self.set_fontfamily(fontname)

    def _ha_for_angle(self, angle):
        """
        Determines horizontal alignment ('ha') for rotation_mode "xtick" based on
        the angle of rotation in degrees and the vertical alignment.
        """
        anchor_at_bottom = self.get_verticalalignment() == 'bottom'
        if (angle <= 10 or 85 <= angle <= 95 or 350 <= angle or
                170 <= angle <= 190 or 265 <= angle <= 275):
            return 'center'
        elif 10 < angle < 85 or 190 < angle < 265:
            return 'left' if anchor_at_bottom else 'right'
        return 'right' if anchor_at_bottom else 'left'

    def _va_for_angle(self, angle):
        """
        Determines vertical alignment ('va') for rotation_mode "ytick" based on
        the angle of rotation in degrees and the horizontal alignment.
        """
        anchor_at_left = self.get_horizontalalignment() == 'left'
        if (angle <= 10 or 350 <= angle or 170 <= angle <= 190
                or 80 <= angle <= 100 or 260 <= angle <= 280):
            return 'center'
        elif 190 < angle < 260 or 10 < angle < 80:
            return 'baseline' if anchor_at_left else 'top'
        return 'top' if anchor_at_left else 'baseline'

    def get_language(self):
        """Return the language this Text is in."""
        return self._language

    def set_language(self, language):
        """
        Set the language of the text.

        Parameters
        ----------
        language : str or None
            The language of the text in a format accepted by libraqm, namely `a BCP47
            language code <https://www.w3.org/International/articles/language-tags/>`_.

            If None, then defaults to :rc:`text.language`.
        """
        _api.check_isinstance((Sequence, str, None), language=language)
        language = mpl._val_or_rc(language, 'text.language')

        if not cbook.is_scalar_or_string(language):
            language = tuple(language)
            for val in language:
                if not isinstance(val, tuple) or len(val) != 3:
                    raise TypeError('language must be list of tuple, not {language!r}')
                sublang, start, end = val
                if not isinstance(sublang, str):
                    raise TypeError(
                        'sub-language specification must be str, not {sublang!r}')
                if not isinstance(start, int):
                    raise TypeError('start location must be int, not {start!r}')
                if not isinstance(end, int):
                    raise TypeError('end location must be int, not {end!r}')

        self._language = language
        self.stale = True


class OffsetFrom:
    """Callable helper class for working with `Annotation`."""

    def __init__(self, artist, ref_coord, unit="points"):
        """
        Parameters
        ----------
        artist : `~matplotlib.artist.Artist` or `.BboxBase` or `.Transform`
            The object to compute the offset from.

        ref_coord : (float, float)
            If *artist* is an `.Artist` or `.BboxBase`, this values is
            the location to of the offset origin in fractions of the
            *artist* bounding box.

            If *artist* is a transform, the offset origin is the
            transform applied to this value.

        unit : {'points, 'pixels'}, default: 'points'
            The screen units to use (pixels or points) for the offset input.
        """
        self._artist = artist
        x, y = ref_coord  # Make copy when ref_coord is an array (and check the shape).
        self._ref_coord = x, y
        self.set_unit(unit)

    def set_unit(self, unit):
        """
        Set the unit for input to the transform used by ``__call__``.

        Parameters
        ----------
        unit : {'points', 'pixels'}
        """
        _api.check_in_list(["points", "pixels"], unit=unit)
        self._unit = unit

    def get_unit(self):
        """Return the unit for input to the transform used by ``__call__``."""
        return self._unit

    def __call__(self, renderer):
        """
        Return the offset transform.

        Parameters
        ----------
        renderer : `RendererBase`
            The renderer to use to compute the offset

        Returns
        -------
        `Transform`
            Maps (x, y) in pixel or point units to screen units
            relative to the given artist.
        """
        if isinstance(self._artist, Artist):
            bbox = self._artist.get_window_extent(renderer)
            xf, yf = self._ref_coord
            x = bbox.x0 + bbox.width * xf
            y = bbox.y0 + bbox.height * yf
        elif isinstance(self._artist, BboxBase):
            bbox = self._artist
            xf, yf = self._ref_coord
            x = bbox.x0 + bbox.width * xf
            y = bbox.y0 + bbox.height * yf
        elif isinstance(self._artist, Transform):
            x, y = self._artist.transform(self._ref_coord)
        else:
            _api.check_isinstance((Artist, BboxBase, Transform), artist=self._artist)
        scale = 1 if self._unit == "pixels" else renderer.points_to_pixels(1)
        return Affine2D().scale(scale).translate(x, y)


class _AnnotationBase:
    def __init__(self,
                 xy,
                 xycoords='data',
                 annotation_clip=None):

        x, y = xy  # Make copy when xy is an array (and check the shape).
        self.xy = x, y
        self.xycoords = xycoords
        self.set_annotation_clip(annotation_clip)

        self._draggable = None

    def _get_xy(self, renderer, xy, coords):
        x, y = xy
        xcoord, ycoord = coords if isinstance(coords, tuple) else (coords, coords)
        if xcoord == 'data':
            x = float(self.convert_xunits(x))
        if ycoord == 'data':
            y = float(self.convert_yunits(y))
        return self._get_xy_transform(renderer, coords).transform((x, y))

    def _get_xy_transform(self, renderer, coords):

        if isinstance(coords, tuple):
            xcoord, ycoord = coords
            from matplotlib.transforms import blended_transform_factory
            tr1 = self._get_xy_transform(renderer, xcoord)
            tr2 = self._get_xy_transform(renderer, ycoord)
            return blended_transform_factory(tr1, tr2)
        elif callable(coords):
            tr = coords(renderer)
            if isinstance(tr, BboxBase):
                return BboxTransformTo(tr)
            elif isinstance(tr, Transform):
                return tr
            else:
                raise TypeError(
                    f"xycoords callable must return a BboxBase or Transform, not a "
                    f"{type(tr).__name__}")
        elif isinstance(coords, Artist):
            bbox = coords.get_window_extent(renderer)
            return BboxTransformTo(bbox)
        elif isinstance(coords, BboxBase):
            return BboxTransformTo(coords)
        elif isinstance(coords, Transform):
            return coords
        elif not isinstance(coords, str):
            raise TypeError(
                f"'xycoords' must be an instance of str, tuple[str, str], Artist, "
                f"Transform, or Callable, not a {type(coords).__name__}")

        if coords == 'data':
            return self.axes.transData
        elif coords == 'polar':
            from matplotlib.projections import PolarAxes
            return PolarAxes.PolarTransform() + self.axes.transData

        try:
            bbox_name, unit = coords.split()
        except ValueError:  # i.e. len(coords.split()) != 2.
            raise ValueError(f"{coords!r} is not a valid coordinate") from None

        bbox0, xy0 = None, None

        # if unit is offset-like
        if bbox_name == "figure":
            bbox0 = self.get_figure(root=False).figbbox
        elif bbox_name == "subfigure":
            bbox0 = self.get_figure(root=False).bbox
        elif bbox_name == "axes":
            bbox0 = self.axes.bbox

        # reference x, y in display coordinate
        if bbox0 is not None:
            xy0 = bbox0.p0
        elif bbox_name == "offset":
            xy0 = self._get_position_xy(renderer)
        else:
            raise ValueError(f"{coords!r} is not a valid coordinate")

        if unit == "points":
            tr = Affine2D().scale(
                self.get_figure(root=True).dpi / 72)  # dpi/72 dots per point
        elif unit == "pixels":
            tr = Affine2D()
        elif unit == "fontsize":
            tr = Affine2D().scale(
                self.get_size() * self.get_figure(root=True).dpi / 72)
        elif unit == "fraction":
            tr = Affine2D().scale(*bbox0.size)
        else:
            raise ValueError(f"{unit!r} is not a recognized unit")

        return tr.translate(*xy0)

    def set_annotation_clip(self, b):
        """
        Set the annotation's clipping behavior.

        Parameters
        ----------
        b : bool or None
            - True: The annotation will be clipped when ``self.xy`` is
              outside the Axes.
            - False: The annotation will always be drawn.
            - None: The annotation will be clipped when ``self.xy`` is
              outside the Axes and ``self.xycoords == "data"``.
        """
        self._annotation_clip = b

    def get_annotation_clip(self):
        """
        Return the annotation's clipping behavior.

        See `set_annotation_clip` for the meaning of return values.
        """
        return self._annotation_clip

    def _get_position_xy(self, renderer):
        """Return the pixel position of the annotated point."""
        return self._get_xy(renderer, self.xy, self.xycoords)

    def _check_xy(self, renderer=None):
        """Check whether the annotation at *xy_pixel* should be drawn."""
        if renderer is None:
            renderer = self.get_figure(root=True)._get_renderer()
        b = self.get_annotation_clip()
        if b or (b is None and self.xycoords == "data"):
            # check if self.xy is inside the Axes.
            xy_pixel = self._get_position_xy(renderer)
            return self.axes.contains_point(xy_pixel)
        return True

    def draggable(self, state=None, use_blit=False):
        """
        Set whether the annotation is draggable with the mouse.

        Parameters
        ----------
        state : bool or None
            - True or False: set the draggability.
            - None: toggle the draggability.
        use_blit : bool, default: False
            Use blitting for faster image composition. For details see
            :ref:`func-animation`.

        Returns
        -------
        DraggableAnnotation or None
            If the annotation is draggable, the corresponding
            `.DraggableAnnotation` helper is returned.
        """
        from matplotlib.offsetbox import DraggableAnnotation
        is_draggable = self._draggable is not None

        # if state is None we'll toggle
        if state is None:
            state = not is_draggable

        if state:
            if self._draggable is None:
                self._draggable = DraggableAnnotation(self, use_blit)
        else:
            if self._draggable is not None:
                self._draggable.disconnect()
            self._draggable = None

        return self._draggable


class Annotation(Text, _AnnotationBase):
    """
    An `.Annotation` is a `.Text` that can refer to a specific position *xy*.
    Optionally an arrow pointing from the text to *xy* can be drawn.

    Attributes
    ----------
    xy
        The annotated position.
    xycoords
        The coordinate system for *xy*.
    arrow_patch
        A `.FancyArrowPatch` to point from *xytext* to *xy*.
    """

    def __str__(self):
        return f"Annotation({self.xy[0]:g}, {self.xy[1]:g}, {self._text!r})"

    def __init__(self, text, xy,
                 xytext=None,
                 xycoords='data',
                 textcoords=None,
                 arrowprops=None,
                 annotation_clip=None,
                 **kwargs):
        """
        Annotate the point *xy* with text *text*.

        In the simplest form, the text is placed at *xy*.

        Optionally, the text can be displayed in another position *xytext*.
        An arrow pointing from the text to the annotated point *xy* can then
        be added by defining *arrowprops*.

        Parameters
        ----------
        text : str
            The text of the annotation.

        xy : (float, float)
            The point *(x, y)* to annotate. The coordinate system is determined
            by *xycoords*.

        xytext : (float, float), default: *xy*
            The position *(x, y)* to place the text at. The coordinate system
            is determined by *textcoords*.

        xycoords : single or two-tuple of str or `.Artist` or `.Transform` or \
callable, default: 'data'

            The coordinate system that *xy* is given in. The following types
            of values are supported:

            - One of the following strings:

              ==================== ============================================
              Value                Description
              ==================== ============================================
              'figure points'      Points from the lower left of the figure
              'figure pixels'      Pixels from the lower left of the figure
              'figure fraction'    Fraction of figure from lower left
              'subfigure points'   Points from the lower left of the subfigure
              'subfigure pixels'   Pixels from the lower left of the subfigure
              'subfigure fraction' Fraction of subfigure from lower left
              'axes points'        Points from lower left corner of the Axes
              'axes pixels'        Pixels from lower left corner of the Axes
              'axes fraction'      Fraction of Axes from lower left
              'data'               Use the coordinate system of the object
                                   being annotated (default)
              'polar'              *(theta, r)* if not native 'data'
                                   coordinates
              ==================== ============================================

              Note that 'subfigure pixels' and 'figure pixels' are the same
              for the parent figure, so users who want code that is usable in
              a subfigure can use 'subfigure pixels'.

            - An `.Artist`: *xy* is interpreted as a fraction of the artist's
              `~matplotlib.transforms.Bbox`. E.g. *(0, 0)* would be the lower
              left corner of the bounding box and *(0.5, 1)* would be the
              center top of the bounding box.

            - A `.Transform` to transform *xy* to screen coordinates.

            - A function with one of the following signatures::

                def transform(renderer) -> Bbox
                def transform(renderer) -> Transform

              where *renderer* is a `.RendererBase` subclass.

              The result of the function is interpreted like the `.Artist` and
              `.Transform` cases above.

            - A tuple *(xcoords, ycoords)* specifying separate coordinate
              systems for *x* and *y*. *xcoords* and *ycoords* must each be
              of one of the above described types.

            See :ref:`plotting-guide-annotation` for more details.

        textcoords : single or two-tuple of str or `.Artist` or `.Transform` \
or callable, default: value of *xycoords*
            The coordinate system that *xytext* is given in.

            All *xycoords* values are valid as well as the following strings:

            =================   =================================================
            Value               Description
            =================   =================================================
            'offset points'     Offset, in points, from the *xy* value
            'offset pixels'     Offset, in pixels, from the *xy* value
            'offset fontsize'   Offset, relative to fontsize, from the *xy* value
            =================   =================================================

        arrowprops : dict, optional
            The properties used to draw a `.FancyArrowPatch` arrow between the
            positions *xy* and *xytext*.  Defaults to None, i.e. no arrow is
            drawn.

            For historical reasons there are two different ways to specify
            arrows, "simple" and "fancy":

            **Simple arrow:**

            If *arrowprops* does not contain the key 'arrowstyle' the
            allowed keys are:

            ==========  =================================================
            Key         Description
            ==========  =================================================
            width       The width of the arrow in points
            headwidth   The width of the base of the arrow head in points
            headlength  The length of the arrow head in points
            shrink      Fraction of total length to shrink from both ends
            ?           Any `.FancyArrowPatch` property
            ==========  =================================================

            The arrow is attached to the edge of the text box, the exact
            position (corners or centers) depending on where it's pointing to.

            **Fancy arrow:**

            This is used if 'arrowstyle' is provided in the *arrowprops*.

            Valid keys are the following `.FancyArrowPatch` parameters:

            ===============  ===================================
            Key              Description
            ===============  ===================================
            arrowstyle       The arrow style
            connectionstyle  The connection style
            relpos           See below; default is (0.5, 0.5)
            patchA           Default is bounding box of the text
            patchB           Default is None
            shrinkA          In points. Default is 2 points
            shrinkB          In points. Default is 2 points
            mutation_scale   Default is text size (in points)
            mutation_aspect  Default is 1
            ?                Any `.FancyArrowPatch` property
            ===============  ===================================

            The exact starting point position of the arrow is defined by
            *relpos*. It's a tuple of relative coordinates of the text box,
            where (0, 0) is the lower left corner and (1, 1) is the upper
            right corner. Values <0 and >1 are supported and specify points
            outside the text box. By default (0.5, 0.5), so the starting point
            is centered in the text box.

        annotation_clip : bool or None, default: None
            Whether to clip (i.e. not draw) the annotation when the annotation
            point *xy* is outside the Axes area.

            - If *True*, the annotation will be clipped when *xy* is outside
              the Axes.
            - If *False*, the annotation will always be drawn.
            - If *None*, the annotation will be clipped when *xy* is outside
              the Axes and *xycoords* is 'data'.

        **kwargs
            Additional kwargs are passed to `.Text`.

        Returns
        -------
        `.Annotation`

        See Also
        --------
        :ref:`annotations`

        """
        _AnnotationBase.__init__(self,
                                 xy,
                                 xycoords=xycoords,
                                 annotation_clip=annotation_clip)
        # warn about wonky input data
        if (xytext is None and
                textcoords is not None and
                textcoords != xycoords):
            _api.warn_external("You have used the `textcoords` kwarg, but "
                               "not the `xytext` kwarg.  This can lead to "
                               "surprising results.")

        # clean up textcoords and assign default
        if textcoords is None:
            textcoords = self.xycoords
        self._textcoords = textcoords

        # cleanup xytext defaults
        if xytext is None:
            xytext = self.xy
        x, y = xytext

        self.arrowprops = arrowprops
        if arrowprops is not None:
            arrowprops = arrowprops.copy()
            if "arrowstyle" in arrowprops:
                self._arrow_relpos = arrowprops.pop("relpos", (0.5, 0.5))
            else:
                # modified YAArrow API to be used with FancyArrowPatch
                for key in ['width', 'headwidth', 'headlength', 'shrink']:
                    arrowprops.pop(key, None)
            self.arrow_patch = FancyArrowPatch((0, 0), (1, 1), **arrowprops)
        else:
            self.arrow_patch = None

        # Must come last, as some kwargs may be propagated to arrow_patch.
        Text.__init__(self, x, y, text, **kwargs)

    def contains(self, mouseevent):
        if self._different_canvas(mouseevent):
            return False, {}
        contains, tinfo = Text.contains(self, mouseevent)
        if self.arrow_patch is not None:
            in_patch, _ = self.arrow_patch.contains(mouseevent)
            contains = contains or in_patch
        return contains, tinfo

    @property
    def xycoords(self):
        return self._xycoords

    @xycoords.setter
    def xycoords(self, xycoords):
        def is_offset(s):
            return isinstance(s, str) and s.startswith("offset")

        if (isinstance(xycoords, tuple) and any(map(is_offset, xycoords))
                or is_offset(xycoords)):
            raise ValueError("xycoords cannot be an offset coordinate")
        self._xycoords = xycoords

    @property
    def xyann(self):
        """
        The text position.

        See also *xytext* in `.Annotation`.
        """
        return self.get_position()

    @xyann.setter
    def xyann(self, xytext):
        self.set_position(xytext)

    def get_anncoords(self):
        """
        Return the coordinate system to use for `.Annotation.xyann`.

        See also *xycoords* in `.Annotation`.
        """
        return self._textcoords

    def set_anncoords(self, coords):
        """
        Set the coordinate system to use for `.Annotation.xyann`.

        See also *xycoords* in `.Annotation`.
        """
        self._textcoords = coords

    anncoords = property(get_anncoords, set_anncoords, doc="""
        The coordinate system to use for `.Annotation.xyann`.""")

    def set_figure(self, fig):
        # docstring inherited
        if self.arrow_patch is not None:
            self.arrow_patch.set_figure(fig)
        Artist.set_figure(self, fig)

    def update_positions(self, renderer):
        """
        Update the pixel positions of the annotation text and the arrow patch.
        """
        # generate transformation
        self.set_transform(self._get_xy_transform(renderer, self.anncoords))

        arrowprops = self.arrowprops
        if arrowprops is None:
            return

        bbox = Text.get_window_extent(self, renderer)

        arrow_end = x1, y1 = self._get_position_xy(renderer)  # Annotated pos.

        ms = arrowprops.get("mutation_scale", self.get_size())
        self.arrow_patch.set_mutation_scale(ms)

        if "arrowstyle" not in arrowprops:
            # Approximately simulate the YAArrow.
            shrink = arrowprops.get('shrink', 0.0)
            width = arrowprops.get('width', 4)
            headwidth = arrowprops.get('headwidth', 12)
            headlength = arrowprops.get('headlength', 12)

            # NB: ms is in pts
            stylekw = dict(head_length=headlength / ms,
                           head_width=headwidth / ms,
                           tail_width=width / ms)

            self.arrow_patch.set_arrowstyle('simple', **stylekw)

            # using YAArrow style:
            # pick the corner of the text bbox closest to annotated point.
            xpos = [(bbox.x0, 0), ((bbox.x0 + bbox.x1) / 2, 0.5), (bbox.x1, 1)]
            ypos = [(bbox.y0, 0), ((bbox.y0 + bbox.y1) / 2, 0.5), (bbox.y1, 1)]
            x, relposx = min(xpos, key=lambda v: abs(v[0] - x1))
            y, relposy = min(ypos, key=lambda v: abs(v[0] - y1))
            self._arrow_relpos = (relposx, relposy)
            r = np.hypot(y - y1, x - x1)
            shrink_pts = shrink * r / renderer.points_to_pixels(1)
            self.arrow_patch.shrinkA = self.arrow_patch.shrinkB = shrink_pts

        # adjust the starting point of the arrow relative to the textbox.
        # TODO : Rotation needs to be accounted.
        arrow_begin = bbox.p0 + bbox.size * self._arrow_relpos
        # The arrow is drawn from arrow_begin to arrow_end.  It will be first
        # clipped by patchA and patchB.  Then it will be shrunk by shrinkA and
        # shrinkB (in points).  If patchA is not set, self.bbox_patch is used.
        self.arrow_patch.set_positions(arrow_begin, arrow_end)

        if "patchA" in arrowprops:
            patchA = arrowprops["patchA"]
        elif self._bbox_patch:
            patchA = self._bbox_patch
        elif self.get_text() == "":
            patchA = None
        else:
            pad = renderer.points_to_pixels(4)
            patchA = Rectangle(
                xy=(bbox.x0 - pad / 2, bbox.y0 - pad / 2),
                width=bbox.width + pad, height=bbox.height + pad,
                transform=IdentityTransform(), clip_on=False)
        self.arrow_patch.set_patchA(patchA)

    @artist.allow_rasterization
    def draw(self, renderer):
        # docstring inherited
        if renderer is not None:
            self._renderer = renderer
        if not self.get_visible() or not self._check_xy(renderer):
            return
        # Update text positions before `Text.draw` would, so that the
        # FancyArrowPatch is correctly positioned.
        self.update_positions(renderer)
        self.update_bbox_position_size(renderer)
        if self.arrow_patch is not None:  # FancyArrowPatch
            if (self.arrow_patch.get_figure(root=False) is None and
                    (fig := self.get_figure(root=False)) is not None):
                self.arrow_patch.set_figure(fig)
            self.arrow_patch.draw(renderer)
        # Draw text, including FancyBboxPatch, after FancyArrowPatch.
        # Otherwise, a wedge arrowstyle can land partly on top of the Bbox.
        Text.draw(self, renderer)

    def get_window_extent(self, renderer=None):
        # docstring inherited
        # This block is the same as in Text.get_window_extent, but we need to
        # set the renderer before calling update_positions().
        if not self.get_visible() or not self._check_xy(renderer):
            return Bbox.unit()
        if renderer is not None:
            self._renderer = renderer
        if self._renderer is None:
            self._renderer = self.get_figure(root=True)._get_renderer()
        if self._renderer is None:
            raise RuntimeError('Cannot get window extent without renderer')

        self.update_positions(self._renderer)

        text_bbox = Text.get_window_extent(self)
        bboxes = [text_bbox]

        if self.arrow_patch is not None:
            bboxes.append(self.arrow_patch.get_window_extent())

        return Bbox.union(bboxes)

    def get_tightbbox(self, renderer=None):
        # docstring inherited
        if not self._check_xy(renderer):
            return Bbox.null()
        return super().get_tightbbox(renderer)


_docstring.interpd.register(Annotation=Annotation.__init__.__doc__)
