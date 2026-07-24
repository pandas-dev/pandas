import numpy as np

from matplotlib import ticker as mticker, _api
from matplotlib.transforms import Bbox, Transform


def _find_line_box_crossings(xys, bbox):
    """
    Find the points where a polyline crosses a bbox, and the crossing angles.

    Parameters
    ----------
    xys : (N, 2) array
        The polyline coordinates.
    bbox : `.Bbox`
        The bounding box.

    Returns
    -------
    list of ((float, float), float)
        Four separate lists of crossings, for the left, right, bottom, and top
        sides of the bbox, respectively.  For each list, the entries are the
        ``((x, y), ccw_angle_in_degrees)`` of the crossing, where an angle of 0
        means that the polyline is moving to the right at the crossing point.

        The entries are computed by linearly interpolating at each crossing
        between the nearest points on either side of the bbox edges.
    """
    crossings = []
    dxys = xys[1:] - xys[:-1]
    for sl in [slice(None), slice(None, None, -1)]:
        us, vs = xys.T[sl]  # "this" coord, "other" coord
        dus, dvs = dxys.T[sl]
        umin, vmin = bbox.min[sl]
        umax, vmax = bbox.max[sl]
        for u0, inside in [(umin, us > umin), (umax, us < umax)]:
            cross = []
            idxs, = (inside[:-1] ^ inside[1:]).nonzero()
            vv = vs[idxs] + (u0 - us[idxs]) * dvs[idxs] / dus[idxs]
            crossings.append([
                ((u0, v)[sl], np.degrees(np.arctan2(*dxy[::-1])))  # ((x, y), theta)
                for v, dxy in zip(vv, dxys[idxs]) if vmin <= v <= vmax])
    return crossings


class ExtremeFinderSimple:
    """
    A helper class to figure out the range of grid lines that need to be drawn.
    """

    def __init__(self, nx, ny):
        """
        Parameters
        ----------
        nx, ny : int
            The number of samples in each direction.
        """
        self.nx = nx
        self.ny = ny

    def __call__(self, transform_xy, x1, y1, x2, y2):
        """
        Compute an approximation of the bounding box obtained by applying
        *transform_xy* to the box delimited by ``(x1, y1, x2, y2)``.

        The intended use is to have ``(x1, y1, x2, y2)`` in axes coordinates,
        and have *transform_xy* be the transform from axes coordinates to data
        coordinates; this method then returns the range of data coordinates
        that span the actual axes.

        The computation is done by sampling ``nx * ny`` equispaced points in
        the ``(x1, y1, x2, y2)`` box and finding the resulting points with
        extremal coordinates; then adding some padding to take into account the
        finite sampling.

        As each sampling step covers a relative range of ``1/nx`` or ``1/ny``,
        the padding is computed by expanding the span covered by the extremal
        coordinates by these fractions.
        """
        tbbox = self._find_transformed_bbox(
            _User2DTransform(transform_xy, None), Bbox.from_extents(x1, y1, x2, y2))
        return tbbox.x0, tbbox.x1, tbbox.y0, tbbox.y1

    def _find_transformed_bbox(self, trans, bbox):
        """
        Compute an approximation of the bounding box obtained by applying
        *trans* to *bbox*.

        See ``__call__`` for details; this method performs similar
        calculations, but using a different representation of the arguments and
        return value.
        """
        grid = np.reshape(np.meshgrid(np.linspace(bbox.x0, bbox.x1, self.nx),
                                      np.linspace(bbox.y0, bbox.y1, self.ny)),
                          (2, -1)).T
        tbbox = Bbox.null()
        tbbox.update_from_data_xy(trans.transform(grid))
        return tbbox.expanded(1 + 2 / self.nx, 1 + 2 / self.ny)


class _User2DTransform(Transform):
    """A transform defined by two user-set functions."""

    input_dims = output_dims = 2

    def __init__(self, forward, backward):
        """
        Parameters
        ----------
        forward, backward : callable
            The forward and backward transforms, taking ``x`` and ``y`` as
            separate arguments and returning ``(tr_x, tr_y)``.
        """
        # The normal Matplotlib convention would be to take and return an
        # (N, 2) array but axisartist uses the transposed version.
        super().__init__()
        self._forward = forward
        self._backward = backward

    def transform_non_affine(self, values):
        # docstring inherited
        return np.transpose(self._forward(*np.transpose(values)))

    def inverted(self):
        # docstring inherited
        return type(self)(self._backward, self._forward)


class GridFinder:
    """
    Internal helper for `~.grid_helper_curvelinear.GridHelperCurveLinear`, with
    the same constructor parameters; should not be directly instantiated.
    """

    def __init__(self,
                 transform,
                 extreme_finder=None,
                 grid_locator1=None,
                 grid_locator2=None,
                 tick_formatter1=None,
                 tick_formatter2=None):
        if extreme_finder is None:
            extreme_finder = ExtremeFinderSimple(20, 20)
        if grid_locator1 is None:
            grid_locator1 = MaxNLocator()
        if grid_locator2 is None:
            grid_locator2 = MaxNLocator()
        if tick_formatter1 is None:
            tick_formatter1 = FormatterPrettyPrint()
        if tick_formatter2 is None:
            tick_formatter2 = FormatterPrettyPrint()
        self.extreme_finder = extreme_finder
        self.grid_locator1 = grid_locator1
        self.grid_locator2 = grid_locator2
        self.tick_formatter1 = tick_formatter1
        self.tick_formatter2 = tick_formatter2
        self.set_transform(transform)

    def _format_ticks(self, idx, direction, factor, levels):
        """
        Helper to support both standard formatters (inheriting from
        `.mticker.Formatter`) and axisartist-specific ones; should be called instead of
        directly calling ``self.tick_formatter1`` and ``self.tick_formatter2``.  This
        method should be considered as a temporary workaround which will be removed in
        the future at the same time as axisartist-specific formatters.
        """
        fmt = _api.getitem_checked(
            {1: self.tick_formatter1, 2: self.tick_formatter2}, idx=idx)
        return (fmt.format_ticks(levels) if isinstance(fmt, mticker.Formatter)
                else fmt(direction, factor, levels))

    def get_grid_info(self, *args, **kwargs):
        """
        Compute positioning information for grid lines and ticks, given the
        axes' data *bbox*.
        """
        params = _api.select_matching_signature(
            [lambda x1, y1, x2, y2: locals(), lambda bbox: locals()], *args, **kwargs)
        if "x1" in params:
            _api.warn_deprecated("3.11", message=(
                "Passing extents as separate arguments to get_grid_info is deprecated "
                "since %(since)s and support will be removed %(removal)s; pass a "
                "single bbox instead."))
            bbox = Bbox.from_extents(
                params["x1"], params["y1"], params["x2"], params["y2"])
        else:
            bbox = params["bbox"]

        tbbox = self.extreme_finder._find_transformed_bbox(
            self.get_transform().inverted(), bbox)

        lon_levs, lon_n, lon_factor = self.grid_locator1(*tbbox.intervalx)
        lat_levs, lat_n, lat_factor = self.grid_locator2(*tbbox.intervaly)

        lon_values = np.asarray(lon_levs[:lon_n]) / lon_factor
        lat_values = np.asarray(lat_levs[:lat_n]) / lat_factor

        lon_lines, lat_lines = self._get_raw_grid_lines(lon_values, lat_values, tbbox)

        bbox_expanded = bbox.expanded(1 + 2e-10, 1 + 2e-10)
        grid_info = {"extremes": tbbox}  # "lon", "lat" keys filled below.

        for idx, lon_or_lat, levs, factor, values, lines in [
                (1, "lon", lon_levs, lon_factor, lon_values, lon_lines),
                (2, "lat", lat_levs, lat_factor, lat_values, lat_lines),
        ]:
            grid_info[lon_or_lat] = gi = {
                "lines": lines,
                "ticks": {"left": [], "right": [], "bottom": [], "top": []},
            }
            for xys, v, level in zip(lines, values, levs):
                all_crossings = _find_line_box_crossings(xys, bbox_expanded)
                for side, crossings in zip(
                        ["left", "right", "bottom", "top"], all_crossings):
                    for crossing in crossings:
                        gi["ticks"][side].append({"level": level, "loc": crossing})
            for side in gi["ticks"]:
                levs = [tick["level"] for tick in gi["ticks"][side]]
                labels = self._format_ticks(idx, side, factor, levs)
                for tick, label in zip(gi["ticks"][side], labels):
                    tick["label"] = label

        return grid_info

    def _get_raw_grid_lines(self, lon_values, lat_values, bbox):
        trans = self.get_transform()
        lons = np.linspace(bbox.x0, bbox.x1, 100)  # for interpolation
        lats = np.linspace(bbox.y0, bbox.y1, 100)
        lon_lines = [trans.transform(np.column_stack([np.full_like(lats, lon), lats]))
                     for lon in lon_values]
        lat_lines = [trans.transform(np.column_stack([lons, np.full_like(lons, lat)]))
                     for lat in lat_values]
        return lon_lines, lat_lines

    def set_transform(self, aux_trans):
        if isinstance(aux_trans, Transform):
            self._aux_transform = aux_trans
        elif len(aux_trans) == 2 and all(map(callable, aux_trans)):
            self._aux_transform = _User2DTransform(*aux_trans)
        else:
            raise TypeError("'aux_trans' must be either a Transform "
                            "instance or a pair of callables")

    def get_transform(self):
        return self._aux_transform

    update_transform = set_transform  # backcompat alias.

    @_api.deprecated("3.11", alternative="grid_finder.get_transform()")
    def transform_xy(self, x, y):
        return self._aux_transform.transform(np.column_stack([x, y])).T

    @_api.deprecated("3.11", alternative="grid_finder.get_transform().inverted()")
    def inv_transform_xy(self, x, y):
        return self._aux_transform.inverted().transform(
            np.column_stack([x, y])).T

    def update(self, **kwargs):
        for k, v in kwargs.items():
            if k in ["extreme_finder",
                     "grid_locator1",
                     "grid_locator2",
                     "tick_formatter1",
                     "tick_formatter2"]:
                setattr(self, k, v)
            else:
                raise ValueError(f"Unknown update property {k!r}")


class MaxNLocator(mticker.MaxNLocator):
    def __init__(self, nbins=10, steps=None,
                 trim=True,
                 integer=False,
                 symmetric=False,
                 prune=None):
        # trim argument has no effect. It has been left for API compatibility
        super().__init__(nbins, steps=steps, integer=integer,
                         symmetric=symmetric, prune=prune)
        self.create_dummy_axis()

    def __call__(self, v1, v2):
        locs = super().tick_values(v1, v2)
        return np.array(locs), len(locs), 1  # 1: factor (see angle_helper)


class FixedLocator:
    def __init__(self, locs):
        self._locs = locs

    def __call__(self, v1, v2):
        v1, v2 = sorted([v1, v2])
        locs = np.array([l for l in self._locs if v1 <= l <= v2])
        return locs, len(locs), 1  # 1: factor (see angle_helper)


# Tick Formatter

class FormatterPrettyPrint:
    def __init__(self, useMathText=True):
        self._fmt = mticker.ScalarFormatter(
            useMathText=useMathText, useOffset=False)
        self._fmt.create_dummy_axis()

    def __call__(self, direction, factor, values):
        return self._fmt.format_ticks(values)


class DictFormatter:
    def __init__(self, format_dict, formatter=None):
        """
        format_dict : dictionary for format strings to be used.
        formatter : fall-back formatter
        """
        super().__init__()
        self._format_dict = format_dict
        self._fallback_formatter = formatter

    def __call__(self, direction, factor, values):
        """
        factor is ignored if value is found in the dictionary
        """
        if self._fallback_formatter:
            fallback_strings = self._fallback_formatter(
                direction, factor, values)
        else:
            fallback_strings = [""] * len(values)
        return [self._format_dict.get(k, v)
                for k, v in zip(values, fallback_strings)]
