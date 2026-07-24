# art3d.py, original mplot3d version by John Porter
# Parts rewritten by Reinier Heeres <reinier@heeres.eu>
# Minor additions by Ben Axelrod <baxelrod@coroware.com>

"""
Module containing 3D artist code and functions to convert 2D
artists into 3D versions which can be added to an Axes3D.
"""

import math

import numpy as np

from contextlib import contextmanager

from matplotlib import (
    _api, artist, cbook, colors as mcolors, lines, text as mtext,
    path as mpath, rcParams)
from matplotlib.collections import (
    Collection, LineCollection, PolyCollection, PatchCollection, PathCollection)
from matplotlib.patches import Patch
from . import proj3d


def _norm_angle(a):
    """Return the given angle normalized to -180 < *a* <= 180 degrees."""
    a = (a + 360) % 360
    if a > 180:
        a = a - 360
    return a


def _norm_text_angle(a):
    """Return the given angle normalized to -90 < *a* <= 90 degrees."""
    a = (a + 180) % 180
    if a > 90:
        a = a - 180
    return a


def get_dir_vector(zdir):
    """
    Return a direction vector.

    Parameters
    ----------
    zdir : {'x', 'y', 'z', None, 3-tuple}
        The direction. Possible values are:

        - 'x': equivalent to (1, 0, 0)
        - 'y': equivalent to (0, 1, 0)
        - 'z': equivalent to (0, 0, 1)
        - *None*: equivalent to (0, 0, 0)
        - an iterable (x, y, z) is converted to an array

    Returns
    -------
    x, y, z : array
        The direction vector.
    """
    if cbook._str_equal(zdir, 'x'):
        return np.array((1, 0, 0))
    elif cbook._str_equal(zdir, 'y'):
        return np.array((0, 1, 0))
    elif cbook._str_equal(zdir, 'z'):
        return np.array((0, 0, 1))
    elif zdir is None:
        return np.array((0, 0, 0))
    elif np.iterable(zdir) and len(zdir) == 3:
        return np.array(zdir)
    else:
        raise ValueError("'x', 'y', 'z', None or vector of length 3 expected")


def _viewlim_mask(xs, ys, zs, axes):
    """
    Return the mask of the points outside the axes view limits.

    Parameters
    ----------
    xs, ys, zs : array-like
        The points to mask. These should be in data coordinates.
    axes : Axes3D
        The axes to use for the view limits.

    Returns
    -------
    mask : np.array
        The mask of the points as a bool array.
    """
    mask = np.logical_or.reduce((xs < axes.xy_viewLim.xmin,
                                 xs > axes.xy_viewLim.xmax,
                                 ys < axes.xy_viewLim.ymin,
                                 ys > axes.xy_viewLim.ymax,
                                 zs < axes.zz_viewLim.xmin,
                                 zs > axes.zz_viewLim.xmax))
    return mask


def _scale_invalid_mask(xs, ys, zs, axes):
    """
    Return the mask of points whose coordinates are invalid for the axis
    scale they live on (e.g. <=0 on a log axis).

    Parameters
    ----------
    xs, ys, zs : array-like
        The points to check, in data coordinates.
    axes : Axes3D
        The axes whose scales are queried.

    Returns
    -------
    mask : np.ndarray
        Boolean array, ``True`` where any of x/y/z is out of its scale's
        valid domain.
    """
    return np.logical_or.reduce((
        np.logical_not(axes.xaxis._scale.val_in_range(xs)),
        np.logical_not(axes.yaxis._scale.val_in_range(ys)),
        np.logical_not(axes.zaxis._scale.val_in_range(zs))))


class Text3D(mtext.Text):
    """
    Text object with 3D position and direction.

    Parameters
    ----------
    x, y, z : float
        The position of the text.
    text : str
        The text string to display.
    zdir : {'x', 'y', 'z', None, 3-tuple}
        The direction of the text. See `.get_dir_vector` for a description of
        the values.
    axlim_clip : bool, default: False
        Whether to hide text outside the axes view limits.

        .. versionadded:: 3.10

    Other Parameters
    ----------------
    **kwargs
         All other parameters are passed on to `~matplotlib.text.Text`.
    """

    def __init__(self, x=0, y=0, z=0, text='', zdir='z', axlim_clip=False,
                 **kwargs):
        if 'rotation' in kwargs:
            _api.warn_external(
                "The `rotation` parameter has not yet been implemented "
                "and is currently ignored."
            )
        if 'rotation_mode' in kwargs:
            _api.warn_external(
                "The `rotation_mode` parameter has not yet been implemented "
                "and is currently ignored."
            )
        mtext.Text.__init__(self, x, y, text, **kwargs)
        self.set_3d_properties(z, zdir, axlim_clip)

    def get_position_3d(self):
        """Return the (x, y, z) position of the text."""
        return self._x, self._y, self._z

    def set_position_3d(self, xyz, zdir=None):
        """
        Set the (*x*, *y*, *z*) position of the text.

        Parameters
        ----------
        xyz : (float, float, float)
            The position in 3D space.
        zdir : {'x', 'y', 'z', None, 3-tuple}
            The direction of the text. If unspecified, the *zdir* will not be
            changed. See `.get_dir_vector` for a description of the values.
        """
        super().set_position(xyz[:2])
        self.set_z(xyz[2])
        if zdir is not None:
            self._dir_vec = get_dir_vector(zdir)

    def set_z(self, z):
        """
        Set the *z* position of the text.

        Parameters
        ----------
        z : float
        """
        self._z = z
        self.stale = True

    def set_3d_properties(self, z=0, zdir='z', axlim_clip=False):
        """
        Set the *z* position and direction of the text.

        Parameters
        ----------
        z : float
            The z-position in 3D space.
        zdir : {'x', 'y', 'z', 3-tuple}
            The direction of the text. Default: 'z'.
            See `.get_dir_vector` for a description of the values.
        axlim_clip : bool, default: False
            Whether to hide text outside the axes view limits.

            .. versionadded:: 3.10
        """
        self._z = z
        self._dir_vec = get_dir_vector(zdir)
        self._axlim_clip = axlim_clip
        self.stale = True

    @artist.allow_rasterization
    def draw(self, renderer):
        mask = _scale_invalid_mask(self._x, self._y, self._z, self.axes)
        if self._axlim_clip:
            mask |= _viewlim_mask(self._x, self._y, self._z, self.axes)
        if np.any(mask):
            pos3d = np.ma.array([self._x, self._y, self._z],
                                mask=mask, dtype=float).filled(np.nan)
        else:
            pos3d = np.array([self._x, self._y, self._z], dtype=float)

        dir_end = pos3d + self._dir_vec
        points = np.asarray([pos3d, dir_end])
        proj = proj3d._scale_proj_transform(
            points[:, 0], points[:, 1], points[:, 2], self.axes)
        dx = proj[0][1] - proj[0][0]
        dy = proj[1][1] - proj[1][0]
        angle = math.degrees(math.atan2(dy, dx))
        with cbook._setattr_cm(self, _x=proj[0][0], _y=proj[1][0],
                               _rotation=_norm_text_angle(angle)):
            mtext.Text.draw(self, renderer)
        self.stale = False

    def get_tightbbox(self, renderer=None):
        # Overwriting the 2d Text behavior which is not valid for 3d.
        # For now, just return None to exclude from layout calculation.
        return None


def text_2d_to_3d(obj, z=0, zdir='z', axlim_clip=False):
    """
    Convert a `.Text` to a `.Text3D` object.

    Parameters
    ----------
    z : float
        The z-position in 3D space.
    zdir : {'x', 'y', 'z', 3-tuple}
        The direction of the text. Default: 'z'.
        See `.get_dir_vector` for a description of the values.
    axlim_clip : bool, default: False
        Whether to hide text outside the axes view limits.

        .. versionadded:: 3.10
    """
    obj.__class__ = Text3D
    obj.set_3d_properties(z, zdir, axlim_clip)


class Line3D(lines.Line2D):
    """
    3D line object.

    .. note:: Use `get_data_3d` to obtain the data associated with the line.
            `~.Line2D.get_data`, `~.Line2D.get_xdata`, and `~.Line2D.get_ydata` return
            the x- and y-coordinates of the projected 2D-line, not the x- and y-data of
            the 3D-line. Similarly, use `set_data_3d` to set the data, not
            `~.Line2D.set_data`, `~.Line2D.set_xdata`, and `~.Line2D.set_ydata`.
    """

    def __init__(self, xs, ys, zs, *args, axlim_clip=False, **kwargs):
        """

        Parameters
        ----------
        xs : array-like
            The x-data to be plotted.
        ys : array-like
            The y-data to be plotted.
        zs : array-like
            The z-data to be plotted.
        *args, **kwargs
            Additional arguments are passed to `~matplotlib.lines.Line2D`.
        """
        super().__init__([], [], *args, **kwargs)
        self.set_data_3d(xs, ys, zs)
        self._axlim_clip = axlim_clip

    def set_3d_properties(self, zs=0, zdir='z', axlim_clip=False):
        """
        Set the *z* position and direction of the line.

        Parameters
        ----------
        zs : float or array of floats
            The location along the *zdir* axis in 3D space to position the
            line.
        zdir : {'x', 'y', 'z'}
            Plane to plot line orthogonal to. Default: 'z'.
            See `.get_dir_vector` for a description of the values.
        axlim_clip : bool, default: False
            Whether to hide lines with an endpoint outside the axes view limits.

            .. versionadded:: 3.10
        """
        xs = self.get_xdata()
        ys = self.get_ydata()
        zs = cbook._to_unmasked_float_array(zs).ravel()
        zs = np.broadcast_to(zs, len(xs))
        self._verts3d = juggle_axes(xs, ys, zs, zdir)
        self._axlim_clip = axlim_clip
        self.stale = True

    def set_data_3d(self, *args):
        """
        Set the x, y and z data

        Parameters
        ----------
        x : array-like
            The x-data to be plotted.
        y : array-like
            The y-data to be plotted.
        z : array-like
            The z-data to be plotted.

        Notes
        -----
        Accepts x, y, z arguments or a single array-like (x, y, z)
        """
        if len(args) == 1:
            args = args[0]
        for name, xyz in zip('xyz', args):
            if not np.iterable(xyz):
                raise RuntimeError(f'{name} must be a sequence')
        self._verts3d = args
        self.stale = True

    def get_data_3d(self):
        """
        Get the current data

        Returns
        -------
        verts3d : length-3 tuple or array-like
            The current data as a tuple or array-like.
        """
        return self._verts3d

    @artist.allow_rasterization
    def draw(self, renderer):
        scale_mask = _scale_invalid_mask(*self._verts3d, self.axes)
        if self._axlim_clip:
            scale_mask |= _viewlim_mask(*self._verts3d, self.axes)
        if np.any(scale_mask):
            mask = np.broadcast_to(
                scale_mask,
                (len(self._verts3d), *self._verts3d[0].shape)
            )
            xs3d, ys3d, zs3d = np.ma.array(self._verts3d,
                                           dtype=float, mask=mask).filled(np.nan)
        else:
            xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs, tis = proj3d._scale_proj_transform_clip(xs3d, ys3d, zs3d, self.axes)
        self.set_data(xs, ys)
        super().draw(renderer)
        self.stale = False


def line_2d_to_3d(line, zs=0, zdir='z', axlim_clip=False):
    """
    Convert a `.Line2D` to a `.Line3D` object.

    Parameters
    ----------
    zs : float
        The location along the *zdir* axis in 3D space to position the line.
    zdir : {'x', 'y', 'z'}
        Plane to plot line orthogonal to. Default: 'z'.
        See `.get_dir_vector` for a description of the values.
    axlim_clip : bool, default: False
        Whether to hide lines with an endpoint outside the axes view limits.

        .. versionadded:: 3.10
    """

    line.__class__ = Line3D
    line.set_3d_properties(zs, zdir, axlim_clip)


def _path_to_3d_segment(path, zs=0, zdir='z'):
    """Convert a path to a 3D segment."""

    zs = np.broadcast_to(zs, len(path))
    pathsegs = path.iter_segments(simplify=False, curves=False)
    seg = [(x, y, z) for (((x, y), code), z) in zip(pathsegs, zs)]
    seg3d = [juggle_axes(x, y, z, zdir) for (x, y, z) in seg]
    return seg3d


def _paths_to_3d_segments(paths, zs=0, zdir='z'):
    """Convert paths from a collection object to 3D segments."""

    if not np.iterable(zs):
        zs = np.broadcast_to(zs, len(paths))
    else:
        if len(zs) != len(paths):
            raise ValueError('Number of z-coordinates does not match paths.')

    segs = [_path_to_3d_segment(path, pathz, zdir)
            for path, pathz in zip(paths, zs)]
    return segs


def _path_to_3d_segment_with_codes(path, zs=0, zdir='z'):
    """Convert a path to a 3D segment with path codes."""

    zs = np.broadcast_to(zs, len(path))
    pathsegs = path.iter_segments(simplify=False, curves=False)
    seg_codes = [((x, y, z), code) for ((x, y), code), z in zip(pathsegs, zs)]
    if seg_codes:
        seg, codes = zip(*seg_codes)
        seg3d = [juggle_axes(x, y, z, zdir) for (x, y, z) in seg]
    else:
        seg3d = []
        codes = []
    return seg3d, list(codes)


def _paths_to_3d_segments_with_codes(paths, zs=0, zdir='z'):
    """
    Convert paths from a collection object to 3D segments with path codes.
    """

    zs = np.broadcast_to(zs, len(paths))
    segments_codes = [_path_to_3d_segment_with_codes(path, pathz, zdir)
                      for path, pathz in zip(paths, zs)]
    if segments_codes:
        segments, codes = zip(*segments_codes)
    else:
        segments, codes = [], []
    return list(segments), list(codes)


class Collection3D(Collection):
    """A collection of 3D paths."""

    def do_3d_projection(self):
        """Project the points according to renderer matrix."""
        vs_list = [vs for vs, _ in self._3dverts_codes]
        masks = [_scale_invalid_mask(*vs.T, self.axes) for vs in vs_list]
        if self._axlim_clip:
            masks = [m | _viewlim_mask(*vs.T, self.axes)
                     for m, vs in zip(masks, vs_list)]
        vs_list = [np.ma.array(vs, mask=np.broadcast_to(m, vs.shape))
                   if np.any(m) else vs
                   for vs, m in zip(vs_list, masks)]
        xyzs_list = [proj3d._scale_proj_transform(
            vs[:, 0], vs[:, 1], vs[:, 2], self.axes) for vs in vs_list]
        self._paths = [mpath.Path(np.ma.column_stack([xs, ys]), cs)
                       for (xs, ys, _), (_, cs) in zip(xyzs_list, self._3dverts_codes)]
        zs = np.concatenate([zs for _, _, zs in xyzs_list])
        return zs.min() if len(zs) else 1e9


def collection_2d_to_3d(col, zs=0, zdir='z', axlim_clip=False):
    """Convert a `.Collection` to a `.Collection3D` object."""
    zs = np.broadcast_to(zs, len(col.get_paths()))
    col._3dverts_codes = [
        (np.column_stack(juggle_axes(
            *np.column_stack([p.vertices, np.broadcast_to(z, len(p.vertices))]).T,
            zdir)),
         p.codes)
        for p, z in zip(col.get_paths(), zs)]
    col.__class__ = cbook._make_class_factory(Collection3D, "{}3D")(type(col))
    col._axlim_clip = axlim_clip


class Line3DCollection(LineCollection):
    """
    A collection of 3D lines.
    """
    def __init__(self, lines, axlim_clip=False, **kwargs):
        super().__init__(lines, **kwargs)
        self._axlim_clip = axlim_clip
        """
        Parameters
        ----------
        lines : list of (N, 3) array-like
            A sequence ``[line0, line1, ...]`` where each line is a (N, 3)-shape
            array-like containing points:: line0 = [(x0, y0, z0), (x1, y1, z1), ...]
            Each line can contain a different number of points.
        linewidths : float or list of float, default: :rc:`lines.linewidth`
            The width of each line in points.
        colors : :mpltype:`color` or list of color, default: :rc:`lines.color`
            A sequence of RGBA tuples (e.g., arbitrary color strings, etc, not
            allowed).
        antialiaseds : bool or list of bool, default: :rc:`lines.antialiased`
            Whether to use antialiasing for each line.
        facecolors : :mpltype:`color` or list of :mpltype:`color`, default: 'none'
            When setting *facecolors*, each line is interpreted as a boundary
            for an area, implicitly closing the path from the last point to the
            first point. The enclosed area is filled with *facecolor*.
            In order to manually specify what should count as the "interior" of
            each line, please use `.PathCollection` instead, where the
            "interior" can be specified by appropriate usage of
            `~.path.Path.CLOSEPOLY`.
        **kwargs : Forwarded to `.Collection`.
        """

    def set_sort_zpos(self, val):
        """Set the position to use for z-sorting."""
        self._sort_zpos = val
        self.stale = True

    def set_segments(self, segments):
        """
        Set 3D segments.
        """
        self._segments3d = segments
        super().set_segments([])

    def do_3d_projection(self):
        """
        Project the points according to renderer matrix.
        """
        segments = self._segments3d

        # Handle ragged inputs, but prefer a faster path for same-length segments
        segment_lengths = [len(segment) for segment in segments]
        ragged = len(set(segment_lengths)) > 1
        if ragged:
            # Branch masked / non-masked for speed
            if any(np.ma.isMA(segment) for segment in segments):
                segments = np.ma.concatenate(segments)
            else:
                segments = np.concatenate(segments)
        else:
            segments = np.asanyarray(segments)

        # Handle empty segments
        if segments.size == 0:
            LineCollection.set_segments(self, [])
            return np.nan

        mask = False
        if np.ma.isMA(segments) and segments.mask is not np.ma.nomask:
            mask = segments.mask

        scale_mask = _scale_invalid_mask(segments[..., 0],
                                         segments[..., 1],
                                         segments[..., 2],
                                         self.axes)
        if np.any(scale_mask):
            mask |= np.broadcast_to(scale_mask[..., np.newaxis],
                                    (*scale_mask.shape, 3))

        if self._axlim_clip:
            viewlim_mask = _viewlim_mask(segments[..., 0],
                                         segments[..., 1],
                                         segments[..., 2],
                                         self.axes)
            if np.any(viewlim_mask):
                # broadcast mask to 3D
                viewlim_mask = np.broadcast_to(viewlim_mask[..., np.newaxis],
                                               (*viewlim_mask.shape, 3))
                mask = mask | viewlim_mask

        xyzs = proj3d._scale_proj_transform_vectors(segments, self.axes)
        if mask is not False:
            xyzs = np.ma.array(xyzs, mask=mask)
        segments_2d = xyzs[..., 0:2]
        if ragged:
            segments_2d = np.split(segments_2d, np.cumsum(segment_lengths[:-1]))
        LineCollection.set_segments(self, segments_2d)

        # FIXME
        if len(xyzs) > 0:
            minz = min(xyzs[..., 2].min(), 1e9)
        else:
            minz = np.nan
        return minz


def line_collection_2d_to_3d(col, zs=0, zdir='z', axlim_clip=False):
    """Convert a `.LineCollection` to a `.Line3DCollection` object."""
    segments3d = _paths_to_3d_segments(col.get_paths(), zs, zdir)
    col.__class__ = Line3DCollection
    col.set_segments(segments3d)
    col._axlim_clip = axlim_clip


class Patch3D(Patch):
    """
    3D patch object.
    """

    def __init__(self, *args, zs=(), zdir='z', axlim_clip=False, **kwargs):
        """
        Parameters
        ----------
        verts :
        zs : float
            The location along the *zdir* axis in 3D space to position the
            patch.
        zdir : {'x', 'y', 'z'}
            Plane to plot patch orthogonal to. Default: 'z'.
            See `.get_dir_vector` for a description of the values.
        axlim_clip : bool, default: False
            Whether to hide patches with a vertex outside the axes view limits.

            .. versionadded:: 3.10
        """
        super().__init__(*args, **kwargs)
        self.set_3d_properties(zs, zdir, axlim_clip)

    def set_3d_properties(self, verts, zs=0, zdir='z', axlim_clip=False):
        """
        Set the *z* position and direction of the patch.

        Parameters
        ----------
        verts :
        zs : float
            The location along the *zdir* axis in 3D space to position the
            patch.
        zdir : {'x', 'y', 'z'}
            Plane to plot patch orthogonal to. Default: 'z'.
            See `.get_dir_vector` for a description of the values.
        axlim_clip : bool, default: False
            Whether to hide patches with a vertex outside the axes view limits.

            .. versionadded:: 3.10
        """
        zs = np.broadcast_to(zs, len(verts))
        self._segment3d = [juggle_axes(x, y, z, zdir)
                           for ((x, y), z) in zip(verts, zs)]
        self._axlim_clip = axlim_clip

    def get_path(self):
        # docstring inherited
        # self._path2d is not initialized until do_3d_projection
        if not hasattr(self, '_path2d'):
            self.axes.M = self.axes.get_proj()
            self.do_3d_projection()
        return self._path2d

    def do_3d_projection(self):
        s = self._segment3d
        xs0, ys0, zs0 = zip(*s)
        mask = _scale_invalid_mask(xs0, ys0, zs0, self.axes)
        if self._axlim_clip:
            mask |= _viewlim_mask(xs0, ys0, zs0, self.axes)
        if np.any(mask):
            xs, ys, zs = np.ma.array(zip(*s),
                                     dtype=float, mask=mask).filled(np.nan)
        else:
            xs, ys, zs = xs0, ys0, zs0
        vxs, vys, vzs, vis = proj3d._scale_proj_transform_clip(xs, ys, zs, self.axes)
        self._path2d = mpath.Path(np.ma.column_stack([vxs, vys]))
        return min(vzs)


class PathPatch3D(Patch3D):
    """
    3D PathPatch object.
    """

    def __init__(self, path, *, zs=(), zdir='z', axlim_clip=False, **kwargs):
        """
        Parameters
        ----------
        path :
        zs : float
            The location along the *zdir* axis in 3D space to position the
            path patch.
        zdir : {'x', 'y', 'z', 3-tuple}
            Plane to plot path patch orthogonal to. Default: 'z'.
            See `.get_dir_vector` for a description of the values.
        axlim_clip : bool, default: False
            Whether to hide path patches with a point outside the axes view limits.

            .. versionadded:: 3.10
        """
        # Not super().__init__!
        Patch.__init__(self, **kwargs)
        self.set_3d_properties(path, zs, zdir, axlim_clip)

    def set_3d_properties(self, path, zs=0, zdir='z', axlim_clip=False):
        """
        Set the *z* position and direction of the path patch.

        Parameters
        ----------
        path :
        zs : float
            The location along the *zdir* axis in 3D space to position the
            path patch.
        zdir : {'x', 'y', 'z', 3-tuple}
            Plane to plot path patch orthogonal to. Default: 'z'.
            See `.get_dir_vector` for a description of the values.
        axlim_clip : bool, default: False
            Whether to hide path patches with a point outside the axes view limits.

            .. versionadded:: 3.10
        """
        Patch3D.set_3d_properties(self, path.vertices, zs=zs, zdir=zdir,
                                  axlim_clip=axlim_clip)
        self._code3d = path.codes

    def do_3d_projection(self):
        s = self._segment3d
        xs0, ys0, zs0 = zip(*s)
        mask = _scale_invalid_mask(xs0, ys0, zs0, self.axes)
        if self._axlim_clip:
            mask |= _viewlim_mask(xs0, ys0, zs0, self.axes)
        if np.any(mask):
            xs, ys, zs = np.ma.array(zip(*s),
                                     dtype=float, mask=mask).filled(np.nan)
        else:
            xs, ys, zs = xs0, ys0, zs0
        vxs, vys, vzs, vis = proj3d._scale_proj_transform_clip(xs, ys, zs, self.axes)
        self._path2d = mpath.Path(np.ma.column_stack([vxs, vys]), self._code3d)
        return min(vzs)


def _get_patch_verts(patch):
    """Return a list of vertices for the path of a patch."""
    trans = patch.get_patch_transform()
    path = patch.get_path()
    polygons = path.to_polygons(trans)
    return polygons[0] if len(polygons) else np.array([])


def patch_2d_to_3d(patch, z=0, zdir='z', axlim_clip=False):
    """Convert a `.Patch` to a `.Patch3D` object."""
    verts = _get_patch_verts(patch)
    patch.__class__ = Patch3D
    patch.set_3d_properties(verts, z, zdir, axlim_clip)


def pathpatch_2d_to_3d(pathpatch, z=0, zdir='z'):
    """Convert a `.PathPatch` to a `.PathPatch3D` object."""
    path = pathpatch.get_path()
    trans = pathpatch.get_patch_transform()

    mpath = trans.transform_path(path)
    pathpatch.__class__ = PathPatch3D
    pathpatch.set_3d_properties(mpath, z, zdir)


class Patch3DCollection(PatchCollection):
    """
    A collection of 3D patches.
    """

    def __init__(
        self,
        *args,
        zs=0,
        zdir="z",
        depthshade=None,
        depthshade_minalpha=None,
        axlim_clip=False,
        **kwargs
    ):
        """
        Create a collection of flat 3D patches with its normal vector
        pointed in *zdir* direction, and located at *zs* on the *zdir*
        axis. 'zs' can be a scalar or an array-like of the same length as
        the number of patches in the collection.

        Constructor arguments are the same as for
        :class:`~matplotlib.collections.PatchCollection`. In addition,
        keywords *zs=0* and *zdir='z'* are available.

        The keyword argument *depthshade* is available to
        indicate whether or not to shade the patches in order to
        give the appearance of depth (default is *True*).
        This is typically desired in scatter plots.

        *depthshade_minalpha* sets the minimum alpha value applied by
        depth-shading.
        """
        if depthshade is None:
            depthshade = rcParams['axes3d.depthshade']
        if depthshade_minalpha is None:
            depthshade_minalpha = rcParams['axes3d.depthshade_minalpha']
        self._depthshade = depthshade
        self._depthshade_minalpha = depthshade_minalpha
        super().__init__(*args, **kwargs)
        self.set_3d_properties(zs, zdir, axlim_clip)

    def get_depthshade(self):
        return self._depthshade

    def set_depthshade(
        self,
        depthshade,
        depthshade_minalpha=None,
    ):
        """
        Set whether depth shading is performed on collection members.

        Parameters
        ----------
        depthshade : bool
            Whether to shade the patches in order to give the appearance of
            depth.
        depthshade_minalpha : float, default: :rc:`axes3d.depthshade_minalpha`
            Sets the minimum alpha value used by depth-shading.

            .. versionadded:: 3.11
        """
        if depthshade_minalpha is None:
            depthshade_minalpha = rcParams['axes3d.depthshade_minalpha']
        self._depthshade = depthshade
        self._depthshade_minalpha = depthshade_minalpha
        self.stale = True

    def set_sort_zpos(self, val):
        """Set the position to use for z-sorting."""
        self._sort_zpos = val
        self.stale = True

    def set_3d_properties(self, zs, zdir, axlim_clip=False):
        """
        Set the *z* positions and direction of the patches.

        Parameters
        ----------
        zs : float or array of floats
            The location or locations to place the patches in the collection
            along the *zdir* axis.
        zdir : {'x', 'y', 'z'}
            Plane to plot patches orthogonal to.
            All patches must have the same direction.
            See `.get_dir_vector` for a description of the values.
        axlim_clip : bool, default: False
            Whether to hide patches with a vertex outside the axes view limits.

            .. versionadded:: 3.10
        """
        # Force the collection to initialize the face and edgecolors
        # just in case it is a scalarmappable with a colormap.
        self.update_scalarmappable()
        offsets = self.get_offsets()
        if len(offsets) > 0:
            xs, ys = offsets.T
        else:
            xs = []
            ys = []
        self._offsets3d = juggle_axes(xs, ys, np.atleast_1d(zs), zdir)
        self._z_markers_idx = slice(-1)
        self._vzs = None
        self._axlim_clip = axlim_clip
        self.stale = True

    def do_3d_projection(self):
        mask = _scale_invalid_mask(*self._offsets3d, self.axes)
        if self._axlim_clip:
            mask |= _viewlim_mask(*self._offsets3d, self.axes)
        if np.any(mask):
            xs, ys, zs = np.ma.array(self._offsets3d, mask=mask)
        else:
            xs, ys, zs = self._offsets3d
        vxs, vys, vzs, vis = proj3d._scale_proj_transform_clip(xs, ys, zs, self.axes)
        self._vzs = vzs
        if np.ma.isMA(vxs):
            super().set_offsets(np.ma.column_stack([vxs, vys]))
        else:
            super().set_offsets(np.column_stack([vxs, vys]))

        if vzs.size > 0:
            return min(vzs)
        else:
            return np.nan

    def _maybe_depth_shade_and_sort_colors(self, color_array):
        # Adjust the color_array alpha values if point depths are defined
        # and depth shading is active
        alpha = self._alpha
        if self._vzs is not None and self._depthshade:
            color_array = _zalpha(
                color_array,
                self._vzs,
                min_alpha=self._depthshade_minalpha,
            )
            if alpha is not None and color_array.shape[1] == 4:  # RGBA, not RGB
                alpha = alpha * color_array[:, 3]

        # Adjust the order of the color_array using the _z_markers_idx,
        # which has been sorted by z-depth
        if len(color_array) > 1:
            color_array = color_array[self._z_markers_idx]
            if np.ndim(alpha) > 0:
                alpha = np.asarray(alpha)[self._z_markers_idx]

        return mcolors.to_rgba_array(color_array, alpha)

    def get_facecolor(self):
        return self._maybe_depth_shade_and_sort_colors(super().get_facecolor())

    def get_edgecolor(self):
        # We need this check here to make sure we do not double-apply the depth
        # based alpha shading when the edge color is "face" which means the
        # edge colour should be identical to the face colour.
        if cbook._str_equal(self._edgecolors, 'face'):
            return self.get_facecolor()
        return self._maybe_depth_shade_and_sort_colors(super().get_edgecolor())


def _get_data_scale(X, Y, Z):
    """
    Estimate the scale of the 3D data for use in depth shading

    Parameters
    ----------
    X, Y, Z : masked arrays
        The data to estimate the scale of.
    """
    # Account for empty datasets. Assume that X Y and Z have the same number
    # of elements.
    if not np.ma.count(X):
        return 0

    # Estimate the scale using the RSS of the ranges of the dimensions
    # Note that we don't use np.ma.ptp() because we otherwise get a build
    # warning about handing empty arrays.
    ptp_x = X.max() - X.min()
    ptp_y = Y.max() - Y.min()
    ptp_z = Z.max() - Z.min()
    return np.sqrt(ptp_x ** 2 + ptp_y ** 2 + ptp_z ** 2)


class Path3DCollection(PathCollection):
    """
    A collection of 3D paths.
    """

    def __init__(
        self,
        *args,
        zs=0,
        zdir="z",
        depthshade=None,
        depthshade_minalpha=None,
        axlim_clip=False,
        **kwargs
    ):
        """
        Create a collection of flat 3D paths with its normal vector
        pointed in *zdir* direction, and located at *zs* on the *zdir*
        axis. 'zs' can be a scalar or an array-like of the same length as
        the number of paths in the collection.

        Constructor arguments are the same as for
        :class:`~matplotlib.collections.PathCollection`. In addition,
        keywords *zs=0* and *zdir='z'* are available.

        Also, the keyword argument *depthshade* is available to
        indicate whether or not to shade the patches in order to
        give the appearance of depth (default is *True*).
        This is typically desired in scatter plots.

        *depthshade_minalpha* sets the minimum alpha value applied by
        depth-shading.
        """
        if depthshade is None:
            depthshade = rcParams['axes3d.depthshade']
        if depthshade_minalpha is None:
            depthshade_minalpha = rcParams['axes3d.depthshade_minalpha']
        self._depthshade = depthshade
        self._depthshade_minalpha = depthshade_minalpha
        self._in_draw = False
        super().__init__(*args, **kwargs)
        self.set_3d_properties(zs, zdir, axlim_clip)
        self._offset_zordered = None

    def draw(self, renderer):
        with self._use_zordered_offset():
            with cbook._setattr_cm(self, _in_draw=True):
                super().draw(renderer)

    def set_sort_zpos(self, val):
        """Set the position to use for z-sorting."""
        self._sort_zpos = val
        self.stale = True

    def set_3d_properties(self, zs, zdir, axlim_clip=False):
        """
        Set the *z* positions and direction of the paths.

        Parameters
        ----------
        zs : float or array of floats
            The location or locations to place the paths in the collection
            along the *zdir* axis.
        zdir : {'x', 'y', 'z'}
            Plane to plot paths orthogonal to.
            All paths must have the same direction.
            See `.get_dir_vector` for a description of the values.
        axlim_clip : bool, default: False
            Whether to hide paths with a vertex outside the axes view limits.

            .. versionadded:: 3.10
        """
        # Force the collection to initialize the face and edgecolors
        # just in case it is a scalarmappable with a colormap.
        self.update_scalarmappable()
        offsets = self.get_offsets()
        if len(offsets) > 0:
            xs, ys = offsets.T
        else:
            xs = []
            ys = []
        self._zdir = zdir
        self._offsets3d = juggle_axes(xs, ys, np.atleast_1d(zs), zdir)
        # In the base draw methods we access the attributes directly which
        # means we cannot resolve the shuffling in the getter methods like
        # we do for the edge and face colors.
        #
        # This means we need to carry around a cache of the unsorted sizes and
        # widths (postfixed with 3d) and in `do_3d_projection` set the
        # depth-sorted version of that data into the private state used by the
        # base collection class in its draw method.
        #
        # Grab the current sizes and linewidths to preserve them.
        self._sizes3d = self._sizes
        self._linewidths3d = np.array(self._linewidths)
        xs, ys, zs = self._offsets3d

        # Sort the points based on z coordinates
        # Performance optimization: Create a sorted index array and reorder
        # points and point properties according to the index array
        self._z_markers_idx = slice(-1)
        self._vzs = None

        self._axlim_clip = axlim_clip
        self.stale = True

    def set_sizes(self, sizes, dpi=72.0):
        super().set_sizes(sizes, dpi)
        if not self._in_draw:
            self._sizes3d = sizes

    def set_linewidth(self, lw):
        super().set_linewidth(lw)
        if not self._in_draw:
            self._linewidths3d = np.array(self._linewidths)

    def get_depthshade(self):
        return self._depthshade

    def set_depthshade(
        self,
        depthshade,
        depthshade_minalpha=None,
    ):
        """
        Set whether depth shading is performed on collection members.

        Parameters
        ----------
        depthshade : bool
            Whether to shade the patches in order to give the appearance of
            depth.
        depthshade_minalpha : float
            Sets the minimum alpha value used by depth-shading.

            .. versionadded:: 3.11
        """
        if depthshade_minalpha is None:
            depthshade_minalpha = rcParams['axes3d.depthshade_minalpha']
        self._depthshade = depthshade
        self._depthshade_minalpha = depthshade_minalpha
        self.stale = True

    def do_3d_projection(self):
        mask = False
        for xyz in self._offsets3d:
            if np.ma.isMA(xyz):
                mask = mask | xyz.mask
        mask = mask | _scale_invalid_mask(*self._offsets3d, self.axes)
        if self._axlim_clip:
            mask = mask | _viewlim_mask(*self._offsets3d, self.axes)
        if np.any(mask):
            mask = np.broadcast_to(mask,
                                   (len(self._offsets3d), *self._offsets3d[0].shape))
            xyzs = np.ma.array(self._offsets3d, mask=mask)
        else:
            xyzs = self._offsets3d
        vxs, vys, vzs, vis = proj3d._scale_proj_transform_clip(*xyzs, self.axes)
        self._data_scale = _get_data_scale(vxs, vys, vzs)
        # Sort the points based on z coordinates
        # Performance optimization: Create a sorted index array and reorder
        # points and point properties according to the index array
        z_markers_idx = self._z_markers_idx = np.ma.argsort(vzs)[::-1]
        self._vzs = vzs

        # we have to special case the sizes because of code in collections.py
        # as the draw method does
        #      self.set_sizes(self._sizes, self.figure.dpi)
        # so we cannot rely on doing the sorting on the way out via get_*

        if len(self._sizes3d) > 1:
            self._sizes = self._sizes3d[z_markers_idx]

        if len(self._linewidths3d) > 1:
            self._linewidths = self._linewidths3d[z_markers_idx]

        PathCollection.set_offsets(self, np.ma.column_stack((vxs, vys)))

        # Re-order items
        vzs = vzs[z_markers_idx]
        vxs = vxs[z_markers_idx]
        vys = vys[z_markers_idx]

        # Store ordered offset for drawing purpose
        self._offset_zordered = np.ma.column_stack((vxs, vys))

        return np.min(vzs) if vzs.size else np.nan

    @contextmanager
    def _use_zordered_offset(self):
        if self._offset_zordered is None:
            # Do nothing
            yield
        else:
            # Swap offset with z-ordered offset
            old_offset = self._offsets
            super().set_offsets(self._offset_zordered)
            try:
                yield
            finally:
                self._offsets = old_offset

    def _maybe_depth_shade_and_sort_colors(self, color_array):
        # Adjust the color_array alpha values if point depths are defined
        # and depth shading is active
        alpha = self._alpha
        if self._vzs is not None and self._depthshade:
            color_array = _zalpha(
                color_array,
                self._vzs,
                min_alpha=self._depthshade_minalpha,
                _data_scale=self._data_scale,
            )
            if alpha is not None and color_array.shape[1] == 4:  # RGBA, not RGB
                alpha = alpha * color_array[:, 3]

        # Adjust the order of the color_array using the _z_markers_idx,
        # which has been sorted by z-depth
        if len(color_array) > 1:
            color_array = color_array[self._z_markers_idx]
            if np.ndim(alpha) > 0:
                alpha = np.asarray(alpha)[self._z_markers_idx]

        return mcolors.to_rgba_array(color_array, alpha)

    def get_facecolor(self):
        return self._maybe_depth_shade_and_sort_colors(super().get_facecolor())

    def get_edgecolor(self):
        # We need this check here to make sure we do not double-apply the depth
        # based alpha shading when the edge color is "face" which means the
        # edge colour should be identical to the face colour.
        if cbook._str_equal(self._edgecolors, 'face'):
            return self.get_facecolor()
        return self._maybe_depth_shade_and_sort_colors(super().get_edgecolor())


def patch_collection_2d_to_3d(
    col,
    zs=0,
    zdir="z",
    depthshade=None,
    axlim_clip=False,
    *args,
    depthshade_minalpha=None,
):
    """
    Convert a `.PatchCollection` into a `.Patch3DCollection` object
    (or a `.PathCollection` into a `.Path3DCollection` object).

    Parameters
    ----------
    col : `~matplotlib.collections.PatchCollection` or \
`~matplotlib.collections.PathCollection`
        The collection to convert.
    zs : float or array of floats
        The location or locations to place the patches in the collection along
        the *zdir* axis. Default: 0.
    zdir : {'x', 'y', 'z'}
        The axis in which to place the patches. Default: "z".
        See `.get_dir_vector` for a description of the values.
    depthshade : bool, default: :rc:`axes3d.depthshade`
        Whether to shade the patches to give a sense of depth.
    axlim_clip : bool, default: False
        Whether to hide patches with a vertex outside the axes view limits.

        .. versionadded:: 3.10

    depthshade_minalpha : float, default: :rc:`axes3d.depthshade_minalpha`
        Sets the minimum alpha value used by depth-shading.

        .. versionadded:: 3.11
    """
    if isinstance(col, PathCollection):
        col.__class__ = Path3DCollection
        col._offset_zordered = None
    elif isinstance(col, PatchCollection):
        col.__class__ = Patch3DCollection
    if depthshade is None:
        depthshade = rcParams['axes3d.depthshade']
    if depthshade_minalpha is None:
        depthshade_minalpha = rcParams['axes3d.depthshade_minalpha']
    col._depthshade = depthshade
    col._depthshade_minalpha = depthshade_minalpha
    col._in_draw = False
    col.set_3d_properties(zs, zdir, axlim_clip)


class Poly3DCollection(PolyCollection):
    """
    A collection of 3D polygons.

    .. note::
        **Filling of 3D polygons**

        There is no simple definition of the enclosed surface of a 3D polygon
        unless the polygon is planar.

        In practice, Matplotlib fills the 2D projection of the polygon. This
        gives a correct filling appearance only for planar polygons. For all
        other polygons, you'll find orientations in which the edges of the
        polygon intersect in the projection. This will lead to an incorrect
        visualization of the 3D area.

        If you need filled areas, it is recommended to create them via
        `~mpl_toolkits.mplot3d.axes3d.Axes3D.plot_trisurf`, which creates a
        triangulation and thus generates consistent surfaces.
    """

    def __init__(self, verts, *args, zsort='average', shade=False,
                 lightsource=None, axlim_clip=False, **kwargs):
        """
        Parameters
        ----------
        verts : list of (N, 3) array-like
            The sequence of polygons [*verts0*, *verts1*, ...] where each
            element *verts_i* defines the vertices of polygon *i* as a 2D
            array-like of shape (N, 3).
        zsort : {'average', 'min', 'max'}, default: 'average'
            The calculation method for the z-order.
            See `~.Poly3DCollection.set_zsort` for details.
        shade : bool, default: False
            Whether to shade *facecolors* and *edgecolors*. When activating
            *shade*, *facecolors* and/or *edgecolors* must be provided.

            .. versionadded:: 3.7

        lightsource : `~matplotlib.colors.LightSource`, optional
            The lightsource to use when *shade* is True.

            .. versionadded:: 3.7

        axlim_clip : bool, default: False
            Whether to hide polygons with a vertex outside the view limits.

            .. versionadded:: 3.10

        *args, **kwargs
            All other parameters are forwarded to `.PolyCollection`.

        Notes
        -----
        Note that this class does a bit of magic with the _facecolors
        and _edgecolors properties.
        """
        if shade:
            normals = _generate_normals(verts)
            facecolors = kwargs.get('facecolors', None)
            if facecolors is not None:
                kwargs['facecolors'] = _shade_colors(
                    facecolors, normals, lightsource
                )

            edgecolors = kwargs.get('edgecolors', None)
            if edgecolors is not None:
                kwargs['edgecolors'] = _shade_colors(
                    edgecolors, normals, lightsource
                )
            if facecolors is None and edgecolors is None:
                raise ValueError(
                    "You must provide facecolors, edgecolors, or both for "
                    "shade to work.")
        super().__init__(verts, *args, **kwargs)
        if isinstance(verts, np.ndarray):
            if verts.ndim != 3:
                raise ValueError('verts must be a list of (N, 3) array-like')
        else:
            if any(len(np.shape(vert)) != 2 for vert in verts):
                raise ValueError('verts must be a list of (N, 3) array-like')
        self.set_zsort(zsort)
        self._codes3d = None
        self._axlim_clip = axlim_clip

    _zsort_functions = {
        'average': np.average,
        'min': np.min,
        'max': np.max,
    }

    def set_zsort(self, zsort):
        """
        Set the calculation method for the z-order.

        Parameters
        ----------
        zsort : {'average', 'min', 'max'}
            The function applied on the z-coordinates of the vertices in the
            viewer's coordinate system, to determine the z-order.
        """
        self._zsortfunc = self._zsort_functions[zsort]
        self._sort_zpos = None
        self.stale = True

    @_api.deprecated("3.10")
    def get_vector(self, segments3d):
        return self._get_vector(segments3d)

    def _get_vector(self, segments3d):
        """
        Optimize points for projection.

        Parameters
        ----------
        segments3d : NumPy array or list of NumPy arrays
            List of vertices of the boundary of every segment. If all paths are
            of equal length and this argument is a NumPy array, then it should
            be of shape (num_faces, num_vertices, 3).
        """
        if isinstance(segments3d, np.ndarray):
            _api.check_shape((None, None, 3), segments3d=segments3d)
            if isinstance(segments3d, np.ma.MaskedArray):
                self._faces = segments3d.data
                self._invalid_vertices = segments3d.mask.any(axis=-1)
            else:
                self._faces = segments3d
                self._invalid_vertices = False
        else:
            # Turn the potentially ragged list into a numpy array for later speedups
            # If it is ragged, set the unused vertices per face as invalid
            num_faces = len(segments3d)
            num_verts = np.fromiter(map(len, segments3d), dtype=np.intp)
            max_verts = num_verts.max(initial=0)
            segments = np.empty((num_faces, max_verts, 3))
            for i, face in enumerate(segments3d):
                segments[i, :len(face)] = face
            self._faces = segments
            self._invalid_vertices = np.arange(max_verts) >= num_verts[:, None]

    def set_verts(self, verts, closed=True):
        """
        Set 3D vertices.

        Parameters
        ----------
        verts : list of (N, 3) array-like
            The sequence of polygons [*verts0*, *verts1*, ...] where each
            element *verts_i* defines the vertices of polygon *i* as a 2D
            array-like of shape (N, 3).
        closed : bool, default: True
            Whether the polygon should be closed by adding a CLOSEPOLY
            connection at the end.
        """
        self._get_vector(verts)
        # 2D verts will be updated at draw time
        super().set_verts([], False)
        self._closed = closed

    def set_verts_and_codes(self, verts, codes):
        """Set 3D vertices with path codes."""
        # set vertices with closed=False to prevent PolyCollection from
        # setting path codes
        self.set_verts(verts, closed=False)
        # and set our own codes instead.
        self._codes3d = codes

    def set_3d_properties(self, axlim_clip=False):
        # Force the collection to initialize the face and edgecolors
        # just in case it is a scalarmappable with a colormap.
        self.update_scalarmappable()
        self._sort_zpos = None
        self.set_zsort('average')
        self._facecolor3d = PolyCollection.get_facecolor(self)
        self._edgecolor3d = PolyCollection.get_edgecolor(self)
        self._alpha3d = PolyCollection.get_alpha(self)
        self.stale = True

    def set_sort_zpos(self, val):
        """Set the position to use for z-sorting."""
        self._sort_zpos = val
        self.stale = True

    def do_3d_projection(self):
        """
        Perform the 3D projection for this object.
        """
        if self._A is not None:
            # force update of color mapping because we re-order them
            # below.  If we do not do this here, the 2D draw will call
            # this, but we will never port the color mapped values back
            # to the 3D versions.
            #
            # We hold the 3D versions in a fixed order (the order the user
            # passed in) and sort the 2D version by view depth.
            self.update_scalarmappable()
            if self._face_is_mapped:
                self._facecolor3d = self._facecolors
            if self._edge_is_mapped:
                self._edgecolor3d = self._edgecolors

        num_faces = len(self._faces)
        mask = self._invalid_vertices | _scale_invalid_mask(
            self._faces[..., 0], self._faces[..., 1],
            self._faces[..., 2], self.axes)
        needs_masking = np.any(mask)

        # Some faces might contain masked vertices, so we want to ignore any
        # errors that those might cause
        with np.errstate(invalid='ignore', divide='ignore'):
            pfaces = proj3d._scale_proj_transform_vectors(self._faces, self.axes)

        if self._axlim_clip:
            viewlim_mask = _viewlim_mask(self._faces[..., 0], self._faces[..., 1],
                                         self._faces[..., 2], self.axes)
            if np.any(viewlim_mask):
                needs_masking = True
                mask = mask | viewlim_mask

        pzs = pfaces[..., 2]
        if needs_masking:
            pzs = np.ma.MaskedArray(pzs, mask=mask)

        # This extra fuss is to re-order face / edge colors
        cface = self._facecolor3d
        cedge = self._edgecolor3d
        if len(cface) != num_faces:
            cface = cface.repeat(num_faces, axis=0)
        if len(cedge) != num_faces:
            if len(cedge) == 0:
                cedge = cface
            else:
                cedge = cedge.repeat(num_faces, axis=0)

        if len(pzs) > 0:
            face_z = self._zsortfunc(pzs, axis=-1)
        else:
            face_z = pzs
        if needs_masking:
            face_z = face_z.data
        face_order = np.argsort(face_z, axis=-1)[::-1]

        if len(pfaces) > 0:
            faces_2d = pfaces[face_order, :, :2]
        else:
            faces_2d = pfaces
        if self._codes3d is not None and len(self._codes3d) > 0:
            if needs_masking:
                segment_mask = ~mask[face_order, :]
                faces_2d = [face[mask, :] for face, mask
                               in zip(faces_2d, segment_mask)]
            codes = [self._codes3d[idx] for idx in face_order]
            PolyCollection.set_verts_and_codes(self, faces_2d, codes)
        else:
            if needs_masking and len(faces_2d) > 0:
                invalid_vertices_2d = np.broadcast_to(
                    mask[face_order, :, None],
                    faces_2d.shape)
                faces_2d = np.ma.MaskedArray(
                        faces_2d, mask=invalid_vertices_2d)
            PolyCollection.set_verts(self, faces_2d, self._closed)

        if len(cface) > 0:
            self._facecolors2d = cface[face_order]
        else:
            self._facecolors2d = cface
        if len(self._edgecolor3d) == len(cface) and len(cedge) > 0:
            self._edgecolors2d = cedge[face_order]
        else:
            self._edgecolors2d = self._edgecolor3d

        # Return zorder value
        if self._sort_zpos is not None:
            zvec = np.array([[0], [0], [self._sort_zpos], [1]])
            ztrans = proj3d._proj_transform_vec(zvec, self.axes.M)
            return ztrans[2][0]
        elif pzs.size > 0:
            # FIXME: Some results still don't look quite right.
            #        In particular, examine contourf3d_demo2.py
            #        with az = -54 and elev = -45.
            return np.min(pzs)
        else:
            return np.nan

    def set_facecolor(self, colors):
        # docstring inherited
        super().set_facecolor(colors)
        self._facecolor3d = PolyCollection.get_facecolor(self)

    def set_edgecolor(self, colors):
        # docstring inherited
        super().set_edgecolor(colors)
        self._edgecolor3d = PolyCollection.get_edgecolor(self)

    def set_alpha(self, alpha):
        # docstring inherited
        artist.Artist.set_alpha(self, alpha)
        try:
            self._facecolor3d = mcolors.to_rgba_array(
                self._facecolor3d, self._alpha)
        except (AttributeError, TypeError, IndexError):
            pass
        try:
            self._edgecolors = mcolors.to_rgba_array(
                    self._edgecolor3d, self._alpha)
        except (AttributeError, TypeError, IndexError):
            pass
        self.stale = True

    def get_facecolor(self):
        # docstring inherited
        # self._facecolors2d is not initialized until do_3d_projection
        if not hasattr(self, '_facecolors2d'):
            self.axes.M = self.axes.get_proj()
            self.do_3d_projection()
        return np.asarray(self._facecolors2d)

    def get_edgecolor(self):
        # docstring inherited
        # self._edgecolors2d is not initialized until do_3d_projection
        if not hasattr(self, '_edgecolors2d'):
            self.axes.M = self.axes.get_proj()
            self.do_3d_projection()
        return np.asarray(self._edgecolors2d)


def poly_collection_2d_to_3d(col, zs=0, zdir='z', axlim_clip=False):
    """
    Convert a `.PolyCollection` into a `.Poly3DCollection` object.

    Parameters
    ----------
    col : `~matplotlib.collections.PolyCollection`
        The collection to convert.
    zs : float or array of floats
        The location or locations to place the polygons in the collection along
        the *zdir* axis. Default: 0.
    zdir : {'x', 'y', 'z'}
        The axis in which to place the patches. Default: 'z'.
        See `.get_dir_vector` for a description of the values.
    """
    segments_3d, codes = _paths_to_3d_segments_with_codes(
            col.get_paths(), zs, zdir)
    col.__class__ = Poly3DCollection
    col.set_verts_and_codes(segments_3d, codes)
    col.set_3d_properties()
    col._axlim_clip = axlim_clip


def juggle_axes(xs, ys, zs, zdir):
    """
    Reorder coordinates so that 2D *xs*, *ys* can be plotted in the plane
    orthogonal to *zdir*. *zdir* is normally 'x', 'y' or 'z'. However, if
    *zdir* starts with a '-' it is interpreted as a compensation for
    `rotate_axes`.
    """
    if zdir == 'x':
        return zs, xs, ys
    elif zdir == 'y':
        return xs, zs, ys
    elif zdir[0] == '-':
        return rotate_axes(xs, ys, zs, zdir)
    else:
        return xs, ys, zs


def rotate_axes(xs, ys, zs, zdir):
    """
    Reorder coordinates so that the axes are rotated with *zdir* along
    the original z axis. Prepending the axis with a '-' does the
    inverse transform, so *zdir* can be 'x', '-x', 'y', '-y', 'z' or '-z'.
    """
    if zdir in ('x', '-y'):
        return ys, zs, xs
    elif zdir in ('-x', 'y'):
        return zs, xs, ys
    else:
        return xs, ys, zs


def _zalpha(
    colors,
    zs,
    min_alpha=0.3,
    _data_scale=None,
):
    """Modify the alpha values of the color list according to z-depth."""

    if len(colors) == 0 or len(zs) == 0:
        return np.zeros((0, 4))

    # Alpha values beyond the range 0-1 inclusive make no sense, so clip them
    min_alpha = np.clip(min_alpha, 0, 1)

    if _data_scale is None or _data_scale == 0:
        # Don't scale the alpha values since we have no valid data scale for reference
        sats = np.ones_like(zs)

    else:
        # Deeper points have an increasingly transparent appearance
        sats = np.clip(1 - (zs - np.min(zs)) / _data_scale, min_alpha, 1)

    rgba = np.broadcast_to(mcolors.to_rgba_array(colors), (len(zs), 4))

    # Change the alpha values of the colors using the generated alpha multipliers
    return np.column_stack([rgba[:, :3], rgba[:, 3] * sats])


def _all_points_on_plane(xs, ys, zs, atol=1e-8):
    """
    Check if all points are on the same plane. Note that NaN values are
    ignored.

    Parameters
    ----------
    xs, ys, zs : array-like
        The x, y, and z coordinates of the points.
    atol : float, default: 1e-8
        The tolerance for the equality check.
    """
    xs, ys, zs = np.asarray(xs), np.asarray(ys), np.asarray(zs)
    points = np.column_stack([xs, ys, zs])
    points = points[~np.isnan(points).any(axis=1)]
    # Check for the case where we have less than 3 unique points
    points = np.unique(points, axis=0)
    if len(points) <= 3:
        return True
    # Calculate the vectors from the first point to all other points
    vs = (points - points[0])[1:]
    vs = vs / np.linalg.norm(vs, axis=1)[:, np.newaxis]
    # Filter out parallel vectors
    vs = np.unique(vs, axis=0)
    if len(vs) <= 2:
        return True
    # Filter out parallel and antiparallel vectors to the first vector
    cross_norms = np.linalg.norm(np.cross(vs[0], vs[1:]), axis=1)
    zero_cross_norms = np.where(np.isclose(cross_norms, 0, atol=atol))[0] + 1
    vs = np.delete(vs, zero_cross_norms, axis=0)
    if len(vs) <= 2:
        return True
    # Calculate the normal vector from the first three points
    n = np.cross(vs[0], vs[1])
    n = n / np.linalg.norm(n)
    # If the dot product of the normal vector and all other vectors is zero,
    # all points are on the same plane
    dots = np.dot(n, vs.transpose())
    return np.allclose(dots, 0, atol=atol)


def _generate_normals(polygons):
    """
    Compute the normals of a list of polygons, one normal per polygon.

    Normals point towards the viewer for a face with its vertices in
    counterclockwise order, following the right hand rule.

    Uses three points equally spaced around the polygon. This method assumes
    that the points are in a plane. Otherwise, more than one shade is required,
    which is not supported.

    Parameters
    ----------
    polygons : list of (M_i, 3) array-like, or (..., M, 3) array-like
        A sequence of polygons to compute normals for, which can have
        varying numbers of vertices. If the polygons all have the same
        number of vertices and array is passed, then the operation will
        be vectorized.

    Returns
    -------
    normals : (..., 3) array
        A normal vector estimated for the polygon.
    """
    if isinstance(polygons, np.ndarray):
        # optimization: polygons all have the same number of points, so can
        # vectorize
        n = polygons.shape[-2]
        i1, i2, i3 = 0, n//3, 2*n//3
        v1 = polygons[..., i1, :] - polygons[..., i2, :]
        v2 = polygons[..., i2, :] - polygons[..., i3, :]
    else:
        # The subtraction doesn't vectorize because polygons is jagged.
        v1 = np.empty((len(polygons), 3))
        v2 = np.empty((len(polygons), 3))
        for poly_i, ps in enumerate(polygons):
            n = len(ps)
            ps = np.asarray(ps)
            i1, i2, i3 = 0, n//3, 2*n//3
            v1[poly_i, :] = ps[i1, :] - ps[i2, :]
            v2[poly_i, :] = ps[i2, :] - ps[i3, :]
    return np.cross(v1, v2)


def _shade_colors(color, normals, lightsource=None):
    """
    Shade *color* using normal vectors given by *normals*,
    assuming a *lightsource* (using default position if not given).
    *color* can also be an array of the same length as *normals*.
    """
    if lightsource is None:
        # chosen for backwards-compatibility
        lightsource = mcolors.LightSource(azdeg=225, altdeg=19.4712)

    with np.errstate(invalid="ignore"):
        shade = ((normals / np.linalg.norm(normals, axis=1, keepdims=True))
                 @ lightsource.direction)
    mask = ~np.isnan(shade)

    if mask.any():
        # convert dot product to allowed shading fractions
        in_norm = mcolors.Normalize(-1, 1)
        out_norm = mcolors.Normalize(0.3, 1).inverse

        def norm(x):
            return out_norm(in_norm(x))

        shade[~mask] = 0

        color = mcolors.to_rgba_array(color)
        # shape of color should be (M, 4) (where M is number of faces)
        # shape of shade should be (M,)
        # colors should have final shape of (M, 4)
        alpha = color[:, 3]
        colors = norm(shade)[:, np.newaxis] * color
        colors[:, 3] = alpha
    else:
        colors = np.asanyarray(color).copy()

    return colors
