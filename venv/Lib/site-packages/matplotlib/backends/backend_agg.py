"""
An `Anti-Grain Geometry`_ (AGG) backend.

Features that are implemented:

* capstyles and join styles
* dashes
* linewidth
* lines, rectangles, ellipses
* clipping to a rectangle
* output to RGBA and Pillow-supported image formats
* alpha blending
* DPI scaling properly - everything scales properly (dashes, linewidths, etc)
* draw polygon
* freetype2 w/ ft2font

Still TODO:

* integrate screen dpi w/ ppi and text

.. _Anti-Grain Geometry: http://agg.sourceforge.net/antigrain.com
"""

from contextlib import nullcontext
import math

import numpy as np
from PIL import features

import matplotlib as mpl
from matplotlib import _api, cbook
from matplotlib.backend_bases import (
    _Backend, FigureCanvasBase, FigureManagerBase, RendererBase)
from matplotlib.dviread import Dvi
from matplotlib.font_manager import fontManager as _fontManager, get_font
from matplotlib.ft2font import LoadFlags, RenderMode
from matplotlib.mathtext import MathTextParser
from matplotlib.path import Path
from matplotlib.transforms import Bbox, BboxBase
from matplotlib.backends._backend_agg import RendererAgg as _RendererAgg


def get_hinting_flag():
    mapping = {
        'default': LoadFlags.DEFAULT,
        'no_autohint': LoadFlags.NO_AUTOHINT,
        'force_autohint': LoadFlags.FORCE_AUTOHINT,
        'no_hinting': LoadFlags.NO_HINTING,
        True: LoadFlags.FORCE_AUTOHINT,
        False: LoadFlags.NO_HINTING,
        'either': LoadFlags.DEFAULT,
        'native': LoadFlags.NO_AUTOHINT,
        'auto': LoadFlags.FORCE_AUTOHINT,
        'none': LoadFlags.NO_HINTING,
    }
    return mapping[mpl.rcParams['text.hinting']]


class RendererAgg(RendererBase):
    """
    The renderer handles all the drawing primitives using a graphics
    context instance that controls the colors/styles
    """

    def __init__(self, width, height, dpi):
        super().__init__()

        self.dpi = dpi
        self.width = width
        self.height = height
        self._renderer = _RendererAgg(int(width), int(height), dpi)
        self._filter_renderers = []

        self._update_methods()
        self.mathtext_parser = MathTextParser('path')

        self.bbox = Bbox.from_bounds(0, 0, self.width, self.height)

    def __getstate__(self):
        # We only want to preserve the init keywords of the Renderer.
        # Anything else can be re-created.
        return {'width': self.width, 'height': self.height, 'dpi': self.dpi}

    def __setstate__(self, state):
        self.__init__(state['width'], state['height'], state['dpi'])

    def _update_methods(self):
        self.draw_gouraud_triangles = self._renderer.draw_gouraud_triangles
        self.draw_image = self._renderer.draw_image
        self.draw_markers = self._renderer.draw_markers
        self.draw_path_collection = self._renderer.draw_path_collection
        self.draw_quad_mesh = self._renderer.draw_quad_mesh
        self.copy_from_bbox = self._renderer.copy_from_bbox

    def draw_path(self, gc, path, transform, rgbFace=None):
        # docstring inherited
        nmax = mpl.rcParams['agg.path.chunksize']  # here at least for testing
        npts = path.vertices.shape[0]

        if (npts > nmax > 100 and path.should_simplify and
                rgbFace is None and gc.get_hatch() is None):
            nch = np.ceil(npts / nmax)
            chsize = int(np.ceil(npts / nch))
            i0 = np.arange(0, npts, chsize)
            i1 = np.zeros_like(i0)
            i1[:-1] = i0[1:] - 1
            i1[-1] = npts
            for ii0, ii1 in zip(i0, i1):
                v = path.vertices[ii0:ii1, :]
                c = path.codes
                if c is not None:
                    c = c[ii0:ii1]
                    c[0] = Path.MOVETO  # move to end of last chunk
                p = Path(v, c)
                p.simplify_threshold = path.simplify_threshold
                try:
                    self._renderer.draw_path(gc, p, transform, rgbFace)
                except OverflowError:
                    msg = (
                        "Exceeded cell block limit in Agg.\n\n"
                        "Please reduce the value of "
                        f"rcParams['agg.path.chunksize'] (currently {nmax}) "
                        "or increase the path simplification threshold"
                        "(rcParams['path.simplify_threshold'] = "
                        f"{mpl.rcParams['path.simplify_threshold']:.2f} by "
                        "default and path.simplify_threshold = "
                        f"{path.simplify_threshold:.2f} on the input)."
                    )
                    raise OverflowError(msg) from None
        else:
            try:
                self._renderer.draw_path(gc, path, transform, rgbFace)
            except OverflowError:
                cant_chunk = ''
                if rgbFace is not None:
                    cant_chunk += "- cannot split filled path\n"
                if gc.get_hatch() is not None:
                    cant_chunk += "- cannot split hatched path\n"
                if not path.should_simplify:
                    cant_chunk += "- path.should_simplify is False\n"
                if len(cant_chunk):
                    msg = (
                        "Exceeded cell block limit in Agg, however for the "
                        "following reasons:\n\n"
                        f"{cant_chunk}\n"
                        "we cannot automatically split up this path to draw."
                        "\n\nPlease manually simplify your path."
                    )

                else:
                    inc_threshold = (
                        "or increase the path simplification threshold"
                        "(rcParams['path.simplify_threshold'] = "
                        f"{mpl.rcParams['path.simplify_threshold']} "
                        "by default and path.simplify_threshold "
                        f"= {path.simplify_threshold} "
                        "on the input)."
                        )
                    if nmax > 100:
                        msg = (
                            "Exceeded cell block limit in Agg.  Please reduce "
                            "the value of rcParams['agg.path.chunksize'] "
                            f"(currently {nmax}) {inc_threshold}"
                        )
                    else:
                        msg = (
                            "Exceeded cell block limit in Agg.  Please set "
                            "the value of rcParams['agg.path.chunksize'], "
                            f"(currently {nmax}) to be greater than 100 "
                            + inc_threshold
                        )

                raise OverflowError(msg) from None

    def _draw_text_glyphs_and_boxes(self, gc, x, y, angle, glyphs, boxes):
        # y is downwards.
        cos = math.cos(math.radians(angle))
        sin = math.sin(math.radians(angle))
        load_flags = get_hinting_flag()
        for font, size, glyph_index, slant, extend, dx, dy in glyphs:  # dy is upwards.
            font.set_size(size, self.dpi)
            font._set_transform(
                (0x10000 * np.array([[cos, -sin], [sin, cos]])
                 @ [[extend, extend * slant], [0, 1]]).round().astype(int),
                [round(0x40 * (x + dx * cos - dy * sin)),
                 # FreeType's y is upwards.
                 round(0x40 * (self.height - y + dx * sin + dy * cos))]
            )
            bitmap = font._render_glyph(
                glyph_index, load_flags,
                RenderMode.NORMAL if gc.get_antialiased() else RenderMode.MONO)
            buffer = bitmap.buffer
            if not gc.get_antialiased():
                buffer *= 0xff
            # draw_text_image's y is downwards & the bitmap bottom side.
            self._renderer.draw_text_image(
                buffer,
                bitmap.left, int(self.height) - bitmap.top + buffer.shape[0],
                0, gc)

        rgba = gc.get_rgb()
        if len(rgba) == 3 or gc.get_forced_alpha():
            rgba = rgba[:3] + (gc.get_alpha(),)
        gc1 = self.new_gc()
        gc1.set_linewidth(0)
        gc1.set_snap(gc.get_snap())
        for dx, dy, w, h in boxes:  # dy is upwards.
            if gc1.get_snap() in [None, True]:
                # Prevent thin bars from disappearing by growing symmetrically.
                if w < 1:
                    dx -= (1 - w) / 2
                    w = 1
                if h < 1:
                    dy -= (1 - h) / 2
                    h = 1
            path = Path._create_closed(
                [(dx, dy), (dx + w, dy), (dx + w, dy + h), (dx, dy + h)])
            self._renderer.draw_path(
                gc1, path,
                mpl.transforms.Affine2D()
                .rotate_deg(angle).translate(x, self.height - y),
                rgba)
        gc1.restore()

    def draw_mathtext(self, gc, x, y, s, prop, angle):
        """Draw mathtext using :mod:`matplotlib.mathtext`."""
        parse = self.mathtext_parser.parse(
            s, self.dpi, prop, antialiased=gc.get_antialiased())
        self._draw_text_glyphs_and_boxes(
            gc, x, y, angle,
            ((font, size, glyph_index, 0, 1, dx, dy)
             for font, size, _char, glyph_index, dx, dy in parse.glyphs),
            parse.rects)

    def draw_text(self, gc, x, y, s, prop, angle, ismath=False, mtext=None):
        # docstring inherited
        if ismath:
            return self.draw_mathtext(gc, x, y, s, prop, angle)
        font = self._prepare_font(prop)
        items = font._layout(
            s, flags=get_hinting_flag(),
            features=mtext.get_fontfeatures() if mtext is not None else None,
            language=mtext.get_language() if mtext is not None else None)
        size = prop.get_size_in_points()
        self._draw_text_glyphs_and_boxes(
            gc, x, y, angle,
            ((item.ft_object, size, item.glyph_index, 0, 1, item.x, item.y)
             for item in items),
            [])

    def get_text_width_height_descent(self, s, prop, ismath):
        # docstring inherited

        _api.check_in_list(["TeX", True, False], ismath=ismath)
        if ismath == "TeX":
            return super().get_text_width_height_descent(s, prop, ismath)

        if ismath:
            parse = self.mathtext_parser.parse(s, self.dpi, prop)
            return parse.width, parse.height, parse.depth

        font = self._prepare_font(prop)
        font.set_text(s, 0.0, flags=get_hinting_flag())
        w, h = font.get_width_height()  # width and height of unrotated string
        d = font.get_descent()
        w /= 64.0  # convert from subpixels
        h /= 64.0
        d /= 64.0
        return w, h, d

    def draw_tex(self, gc, x, y, s, prop, angle, *, mtext=None):
        # docstring inherited
        # todo, handle props, angle, origins

        size = prop.get_size_in_points()

        if mpl.rcParams["text.latex.engine"] == "latex+dvipng":
            Z = self.get_texmanager().get_grey(s, size, self.dpi)
            Z = (Z * 0xff).astype(np.uint8)
            w, h, d = self.get_text_width_height_descent(s, prop, ismath="TeX")
            xd = d * math.sin(math.radians(angle))
            yd = d * math.cos(math.radians(angle))
            x = round(x + xd)
            y = round(y + yd)
            self._renderer.draw_text_image(Z, x, y, angle, gc)
            return

        dvifile = self.get_texmanager().make_dvi(s, size)
        with Dvi(dvifile, self.dpi) as dvi:
            page, = dvi

        self._draw_text_glyphs_and_boxes(
            gc, x, y, angle,
            ((get_font(text.font_path), text.font_size, text.index,
              text.font_effects.get('slant', 0), text.font_effects.get('extend', 1),
              text.x, text.y)
             for text in page.text),
            ((box.x, box.y, box.width, box.height) for box in page.boxes))

    def get_canvas_width_height(self):
        # docstring inherited
        return self.width, self.height

    def _prepare_font(self, font_prop):
        """
        Get the `.FT2Font` for *font_prop*, clear its buffer, and set its size.
        """
        font = get_font(_fontManager._find_fonts_by_props(font_prop))
        font.clear()
        size = font_prop.get_size_in_points()
        font.set_size(size, self.dpi)
        return font

    def points_to_pixels(self, points):
        # docstring inherited
        return points * self.dpi / 72

    def buffer_rgba(self):
        return memoryview(self._renderer)

    def tostring_argb(self):
        return np.asarray(self._renderer).take([3, 0, 1, 2], axis=2).tobytes()

    def clear(self):
        self._renderer.clear()

    def option_image_nocomposite(self):
        # docstring inherited

        # It is generally faster to composite each image directly to
        # the Figure, and there's no file size benefit to compositing
        # with the Agg backend
        return True

    def option_scale_image(self):
        # docstring inherited
        return False

    def restore_region(self, region, bbox=None, xy=None):
        """
        Restore the saved region. If bbox (instance of BboxBase, or
        its extents) is given, only the region specified by the bbox
        will be restored. *xy* (a pair of floats) optionally
        specifies the new position (the LLC of the original region,
        not the LLC of the bbox) where the region will be restored.

        >>> region = renderer.copy_from_bbox()
        >>> x1, y1, x2, y2 = region.get_extents()
        >>> renderer.restore_region(region, bbox=(x1+dx, y1, x2, y2),
        ...                         xy=(x1-dx, y1))

        """
        if bbox is not None or xy is not None:
            if bbox is None:
                x1, y1, x2, y2 = region.get_extents()
            elif isinstance(bbox, BboxBase):
                x1, y1, x2, y2 = bbox.extents
            else:
                x1, y1, x2, y2 = bbox

            if xy is None:
                ox, oy = x1, y1
            else:
                ox, oy = xy

            # The incoming data is float, but the _renderer type-checking wants
            # to see integers.
            self._renderer.restore_region(region, int(x1), int(y1),
                                          int(x2), int(y2), int(ox), int(oy))

        else:
            self._renderer.restore_region(region)

    def start_filter(self):
        """
        Start filtering. It simply creates a new canvas (the old one is saved).
        """
        self._filter_renderers.append(self._renderer)
        self._renderer = _RendererAgg(int(self.width), int(self.height),
                                      self.dpi)
        self._update_methods()

    def stop_filter(self, post_processing):
        """
        Save the current canvas as an image and apply post processing.

        The *post_processing* function::

           def post_processing(image, dpi):
             # ny, nx, depth = image.shape
             # image (numpy array) has RGBA channels and has a depth of 4.
             ...
             # create a new_image (numpy array of 4 channels, size can be
             # different). The resulting image may have offsets from
             # lower-left corner of the original image
             return new_image, offset_x, offset_y

        The saved renderer is restored and the returned image from
        post_processing is plotted (using draw_image) on it.
        """
        orig_img = np.asarray(self.buffer_rgba())
        slice_y, slice_x = cbook._get_nonzero_slices(orig_img[..., 3])
        cropped_img = orig_img[slice_y, slice_x]

        self._renderer = self._filter_renderers.pop()
        self._update_methods()

        if cropped_img.size:
            img, ox, oy = post_processing(cropped_img / 255, self.dpi)
            gc = self.new_gc()
            if img.dtype.kind == 'f':
                img = np.asarray(img * 255., np.uint8)
            self._renderer.draw_image(
                gc, slice_x.start + ox, int(self.height) - slice_y.stop + oy,
                img[::-1])


class FigureCanvasAgg(FigureCanvasBase):
    # docstring inherited

    _lastKey = None  # Overwritten per-instance on the first draw.

    def copy_from_bbox(self, bbox):
        renderer = self.get_renderer()
        return renderer.copy_from_bbox(bbox)

    def restore_region(self, region, bbox=None, xy=None):
        renderer = self.get_renderer()
        return renderer.restore_region(region, bbox, xy)

    def draw(self):
        # docstring inherited
        self.renderer = self.get_renderer()
        self.renderer.clear()
        # Acquire a lock on the shared font cache.
        with (self.toolbar._wait_cursor_for_draw_cm() if self.toolbar
              else nullcontext()):
            self.figure.draw(self.renderer)
            # A GUI class may be need to update a window using this draw, so
            # don't forget to call the superclass.
            super().draw()

    def get_renderer(self):
        w, h = self.get_width_height(physical=True)
        key = w, h, self.figure.dpi
        reuse_renderer = (self._lastKey == key)
        if not reuse_renderer:
            self.renderer = RendererAgg(w, h, self.figure.dpi)
            self._lastKey = key
        return self.renderer

    def tostring_argb(self):
        """
        Get the image as ARGB `bytes`.

        `draw` must be called at least once before this function will work and
        to update the renderer for any subsequent changes to the Figure.
        """
        return self.renderer.tostring_argb()

    def buffer_rgba(self):
        """
        Get the image as a `memoryview` to the renderer's buffer.

        `draw` must be called at least once before this function will work and
        to update the renderer for any subsequent changes to the Figure.
        """
        return self.renderer.buffer_rgba()

    def print_raw(self, filename_or_obj, *, metadata=None):
        if metadata is not None:
            raise ValueError("metadata not supported for raw/rgba")
        FigureCanvasAgg.draw(self)
        renderer = self.get_renderer()
        with cbook.open_file_cm(filename_or_obj, "wb") as fh:
            fh.write(renderer.buffer_rgba())

    print_rgba = print_raw

    def _print_pil(self, filename_or_obj, fmt, pil_kwargs, metadata=None):
        """
        Draw the canvas, then save it using `.image.imsave` (to which
        *pil_kwargs* and *metadata* are forwarded).
        """
        FigureCanvasAgg.draw(self)
        mpl.image.imsave(
            filename_or_obj, self.buffer_rgba(), format=fmt, origin="upper",
            dpi=self.figure.dpi, metadata=metadata, pil_kwargs=pil_kwargs)

    def print_png(self, filename_or_obj, *, metadata=None, pil_kwargs=None):
        """
        Write the figure to a PNG file.

        Parameters
        ----------
        filename_or_obj : str or path-like or file-like
            The file to write to.

        metadata : dict, optional
            Metadata in the PNG file as key-value pairs of bytes or latin-1
            encodable strings.
            According to the PNG specification, keys must be shorter than 79
            chars.

            The `PNG specification`_ defines some common keywords that may be
            used as appropriate:

            - Title: Short (one line) title or caption for image.
            - Author: Name of image's creator.
            - Description: Description of image (possibly long).
            - Copyright: Copyright notice.
            - Creation Time: Time of original image creation
              (usually RFC 1123 format).
            - Software: Software used to create the image.
            - Disclaimer: Legal disclaimer.
            - Warning: Warning of nature of content.
            - Source: Device used to create the image.
            - Comment: Miscellaneous comment;
              conversion from other image format.

            Other keywords may be invented for other purposes.

            If 'Software' is not given, an autogenerated value for Matplotlib
            will be used.  This can be removed by setting it to *None*.

            For more details see the `PNG specification`_.

            .. _PNG specification: \
                https://www.w3.org/TR/2003/REC-PNG-20031110/#11keywords

        pil_kwargs : dict, optional
            Keyword arguments passed to `PIL.Image.Image.save`.

            If the 'pnginfo' key is present, it completely overrides
            *metadata*, including the default 'Software' key.
        """
        self._print_pil(filename_or_obj, "png", pil_kwargs, metadata)

    def print_to_buffer(self):
        FigureCanvasAgg.draw(self)
        renderer = self.get_renderer()
        return (bytes(renderer.buffer_rgba()),
                (int(renderer.width), int(renderer.height)))

    # Note that these methods should typically be called via savefig() and
    # print_figure(), and the latter ensures that `self.figure.dpi` already
    # matches the dpi kwarg (if any).

    def print_gif(self, filename_or_obj, *, metadata=None, pil_kwargs=None):
        self._print_pil(filename_or_obj, "gif", pil_kwargs, metadata)

    def print_jpg(self, filename_or_obj, *, metadata=None, pil_kwargs=None):
        # savefig() has already applied savefig.facecolor; we now set it to
        # white to make imsave() blend semi-transparent figures against an
        # assumed white background.
        with mpl.rc_context({"savefig.facecolor": "white"}):
            self._print_pil(filename_or_obj, "jpeg", pil_kwargs, metadata)

    print_jpeg = print_jpg

    def print_tif(self, filename_or_obj, *, metadata=None, pil_kwargs=None):
        self._print_pil(filename_or_obj, "tiff", pil_kwargs, metadata)

    print_tiff = print_tif

    def print_webp(self, filename_or_obj, *, metadata=None, pil_kwargs=None):
        self._print_pil(filename_or_obj, "webp", pil_kwargs, metadata)

    def print_avif(self, filename_or_obj, *, metadata=None, pil_kwargs=None):
        if not features.check("avif"):
            raise RuntimeError(
                "The installed pillow version does not support avif. Full "
                "avif support has been added in pillow 11.3."
            )
        self._print_pil(filename_or_obj, "avif", pil_kwargs, metadata)

    (print_gif.__doc__,
     print_jpg.__doc__,
     print_tif.__doc__,
     print_webp.__doc__,
     print_avif.__doc__) = map(
        """
        Write the figure to a {} file.

        Parameters
        ----------
        filename_or_obj : str or path-like or file-like
            The file to write to.
        metadata : None
            Unused for pillow-based writers. All supported options
            can be passed via *pil_kwargs*.
        pil_kwargs : dict, optional
            Additional keyword arguments that are passed to
            `PIL.Image.Image.save` when saving the figure.
        """.format, ["GIF", "JPEG", "TIFF", "WebP", "AVIF"])


@_Backend.export
class _BackendAgg(_Backend):
    backend_version = 'v2.2'
    FigureCanvas = FigureCanvasAgg
    FigureManager = FigureManagerBase
