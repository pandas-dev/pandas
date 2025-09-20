from __future__ import annotations

import functools
import itertools
import warnings
from collections.abc import Callable, Hashable, Iterable, MutableMapping
from typing import TYPE_CHECKING, Any, Generic, Literal, TypeVar, cast

import numpy as np

from xarray.core.formatting import format_item
from xarray.core.types import HueStyleOptions, T_DataArrayOrSet
from xarray.plot.utils import (
    _LINEWIDTH_RANGE,
    _MARKERSIZE_RANGE,
    _add_legend,
    _determine_guide,
    _get_nice_quiver_magnitude,
    _guess_coords_to_plot,
    _infer_xy_labels,
    _Normalize,
    _parse_size,
    _process_cmap_cbar_kwargs,
    label_from_attrs,
)

if TYPE_CHECKING:
    from matplotlib.axes import Axes
    from matplotlib.cm import ScalarMappable
    from matplotlib.colorbar import Colorbar
    from matplotlib.figure import Figure
    from matplotlib.legend import Legend
    from matplotlib.quiver import QuiverKey
    from matplotlib.text import Annotation

    from xarray.core.dataarray import DataArray


# Overrides axes.labelsize, xtick.major.size, ytick.major.size
# from mpl.rcParams
_FONTSIZE = "small"
# For major ticks on x, y axes
_NTICKS = 5


def _nicetitle(coord, value, maxchar, template):
    """
    Put coord, value in template and truncate at maxchar
    """
    prettyvalue = format_item(value, quote_strings=False)
    title = template.format(coord=coord, value=prettyvalue)

    if len(title) > maxchar:
        title = title[: (maxchar - 3)] + "..."

    return title


T_FacetGrid = TypeVar("T_FacetGrid", bound="FacetGrid")


class FacetGrid(Generic[T_DataArrayOrSet]):
    """
    Initialize the Matplotlib figure and FacetGrid object.

    The :class:`FacetGrid` is an object that links a xarray DataArray to
    a Matplotlib figure with a particular structure.

    In particular, :class:`FacetGrid` is used to draw plots with multiple
    axes, where each axes shows the same relationship conditioned on
    different levels of some dimension. It's possible to condition on up to
    two variables by assigning variables to the rows and columns of the
    grid.

    The general approach to plotting here is called "small multiples",
    where the same kind of plot is repeated multiple times, and the
    specific use of small multiples to display the same relationship
    conditioned on one or more other variables is often called a "trellis
    plot".

    The basic workflow is to initialize the :class:`FacetGrid` object with
    the DataArray and the variable names that are used to structure the grid.
    Then plotting functions can be applied to each subset by calling
    :meth:`FacetGrid.map_dataarray` or :meth:`FacetGrid.map`.

    Attributes
    ----------
    axs : ndarray of matplotlib.axes.Axes
        Array containing axes in corresponding position, as returned from
        :py:func:`matplotlib.pyplot.subplots`.
    col_labels : list of matplotlib.text.Annotation
        Column titles.
    row_labels : list of matplotlib.text.Annotation
        Row titles.
    fig : matplotlib.figure.Figure
        The figure containing all the axes.
    name_dicts : ndarray of dict
        Array containing dictionaries mapping coordinate names to values. ``None`` is
        used as a sentinel value for axes that should remain empty, i.e.,
        sometimes the rightmost grid positions in the bottom row.
    """

    data: T_DataArrayOrSet
    name_dicts: np.ndarray
    fig: Figure
    axs: np.ndarray
    row_names: list[np.ndarray]
    col_names: list[np.ndarray]
    figlegend: Legend | None
    quiverkey: QuiverKey | None
    cbar: Colorbar | None
    _single_group: bool | Hashable
    _nrow: int
    _row_var: Hashable | None
    _ncol: int
    _col_var: Hashable | None
    _col_wrap: int | None
    row_labels: list[Annotation | None]
    col_labels: list[Annotation | None]
    _x_var: None
    _y_var: None
    _hue_var: DataArray | None
    _cmap_extend: Any | None
    _mappables: list[ScalarMappable]
    _finalized: bool

    def __init__(
        self,
        data: T_DataArrayOrSet,
        col: Hashable | None = None,
        row: Hashable | None = None,
        col_wrap: int | None = None,
        sharex: bool = True,
        sharey: bool = True,
        figsize: Iterable[float] | None = None,
        aspect: float = 1,
        size: float = 3,
        subplot_kws: dict[str, Any] | None = None,
    ) -> None:
        """
        Parameters
        ----------
        data : DataArray or Dataset
            DataArray or Dataset to be plotted.
        row, col : str
            Dimension names that define subsets of the data, which will be drawn
            on separate facets in the grid.
        col_wrap : int, optional
            "Wrap" the grid the for the column variable after this number of columns,
            adding rows if ``col_wrap`` is less than the number of facets.
        sharex : bool, optional
            If true, the facets will share *x* axes.
        sharey : bool, optional
            If true, the facets will share *y* axes.
        figsize : Iterable of float or None, optional
            A tuple (width, height) of the figure in inches.
            If set, overrides ``size`` and ``aspect``.
        aspect : scalar, default: 1
            Aspect ratio of each facet, so that ``aspect * size`` gives the
            width of each facet in inches.
        size : scalar, default: 3
            Height (in inches) of each facet. See also: ``aspect``.
        subplot_kws : dict, optional
            Dictionary of keyword arguments for Matplotlib subplots
            (:py:func:`matplotlib.pyplot.subplots`).

        """

        import matplotlib.pyplot as plt

        # Handle corner case of nonunique coordinates
        rep_col = col is not None and not data[col].to_index().is_unique
        rep_row = row is not None and not data[row].to_index().is_unique
        if rep_col or rep_row:
            raise ValueError(
                "Coordinates used for faceting cannot "
                "contain repeated (nonunique) values."
            )

        # single_group is the grouping variable, if there is exactly one
        single_group: bool | Hashable
        if col and row:
            single_group = False
            nrow = len(data[row])
            ncol = len(data[col])
            nfacet = nrow * ncol
            if col_wrap is not None:
                warnings.warn(
                    "Ignoring col_wrap since both col and row were passed", stacklevel=2
                )
        elif row and not col:
            single_group = row
        elif not row and col:
            single_group = col
        else:
            raise ValueError("Pass a coordinate name as an argument for row or col")

        # Compute grid shape
        if single_group:
            nfacet = len(data[single_group])
            if col:
                # idea - could add heuristic for nice shapes like 3x4
                ncol = nfacet
            if row:
                ncol = 1
            if col_wrap is not None:
                # Overrides previous settings
                ncol = col_wrap
            nrow = int(np.ceil(nfacet / ncol))

        # Set the subplot kwargs
        subplot_kws = {} if subplot_kws is None else subplot_kws

        if figsize is None:
            # Calculate the base figure size with extra horizontal space for a
            # colorbar
            cbar_space = 1
            figsize = (ncol * size * aspect + cbar_space, nrow * size)

        fig, axs = plt.subplots(
            nrow,
            ncol,
            sharex=sharex,
            sharey=sharey,
            squeeze=False,
            figsize=figsize,
            subplot_kw=subplot_kws,
        )

        # Set up the lists of names for the row and column facet variables
        col_names = list(data[col].to_numpy()) if col else []
        row_names = list(data[row].to_numpy()) if row else []

        if single_group:
            full: list[dict[Hashable, Any] | None] = [
                {single_group: x} for x in data[single_group].to_numpy()
            ]
            empty: list[dict[Hashable, Any] | None] = [
                None for x in range(nrow * ncol - len(full))
            ]
            name_dict_list = full + empty
        else:
            rowcols = itertools.product(row_names, col_names)
            name_dict_list = [{row: r, col: c} for r, c in rowcols]

        name_dicts = np.array(name_dict_list).reshape(nrow, ncol)

        # Set up the class attributes
        # ---------------------------

        # First the public API
        self.data = data
        self.name_dicts = name_dicts
        self.fig = fig
        self.axs = axs
        self.row_names = row_names
        self.col_names = col_names

        # guides
        self.figlegend = None
        self.quiverkey = None
        self.cbar = None

        # Next the private variables
        self._single_group = single_group
        self._nrow = nrow
        self._row_var = row
        self._ncol = ncol
        self._col_var = col
        self._col_wrap = col_wrap
        self.row_labels = [None] * nrow
        self.col_labels = [None] * ncol
        self._x_var = None
        self._y_var = None
        self._hue_var = None
        self._cmap_extend = None
        self._mappables = []
        self._finalized = False

    @property
    def axes(self) -> np.ndarray:
        warnings.warn(
            (
                "self.axes is deprecated since 2022.11 in order to align with "
                "matplotlibs plt.subplots, use self.axs instead."
            ),
            DeprecationWarning,
            stacklevel=2,
        )
        return self.axs

    @axes.setter
    def axes(self, axs: np.ndarray) -> None:
        warnings.warn(
            (
                "self.axes is deprecated since 2022.11 in order to align with "
                "matplotlibs plt.subplots, use self.axs instead."
            ),
            DeprecationWarning,
            stacklevel=2,
        )
        self.axs = axs

    @property
    def _left_axes(self) -> np.ndarray:
        return self.axs[:, 0]

    @property
    def _bottom_axes(self) -> np.ndarray:
        return self.axs[-1, :]

    def map_dataarray(
        self: T_FacetGrid,
        func: Callable,
        x: Hashable | None,
        y: Hashable | None,
        **kwargs: Any,
    ) -> T_FacetGrid:
        """
        Apply a plotting function to a 2d facet's subset of the data.

        This is more convenient and less general than ``FacetGrid.map``

        Parameters
        ----------
        func : callable
            A plotting function with the same signature as a 2d xarray
            plotting method such as `xarray.plot.imshow`
        x, y : string
            Names of the coordinates to plot on x, y axes
        **kwargs
            additional keyword arguments to func

        Returns
        -------
        self : FacetGrid object

        """

        if kwargs.get("cbar_ax") is not None:
            raise ValueError("cbar_ax not supported by FacetGrid.")

        cmap_params, cbar_kwargs = _process_cmap_cbar_kwargs(
            func, self.data.to_numpy(), **kwargs
        )

        self._cmap_extend = cmap_params.get("extend")

        # Order is important
        func_kwargs = {
            k: v
            for k, v in kwargs.items()
            if k not in {"cmap", "colors", "cbar_kwargs", "levels"}
        }
        func_kwargs.update(cmap_params)
        # to avoid redundant calling, colorbar and labelling is instead handled
        # by `_finalize_grid` at the end
        func_kwargs["add_colorbar"] = False
        if func.__name__ != "surface":
            func_kwargs["add_labels"] = False

        # Get x, y labels for the first subplot
        x, y = _infer_xy_labels(
            darray=self.data.loc[self.name_dicts.flat[0]],
            x=x,
            y=y,
            imshow=func.__name__ == "imshow",
            rgb=kwargs.get("rgb"),
        )

        for d, ax in zip(self.name_dicts.flat, self.axs.flat, strict=True):
            # None is the sentinel value
            if d is not None:
                subset = self.data.loc[d]
                mappable = func(
                    subset, x=x, y=y, ax=ax, **func_kwargs, _is_facetgrid=True
                )
                self._mappables.append(mappable)

        xlabel = label_from_attrs(self.data[x])
        ylabel = label_from_attrs(self.data[y])

        self._finalize_grid(xlabel, ylabel)

        if kwargs.get("add_colorbar", True):
            self.add_colorbar(**cbar_kwargs)

        return self

    def map_plot1d(
        self: T_FacetGrid,
        func: Callable,
        x: Hashable | None,
        y: Hashable | None,
        *,
        z: Hashable | None = None,
        hue: Hashable | None = None,
        markersize: Hashable | None = None,
        linewidth: Hashable | None = None,
        **kwargs: Any,
    ) -> T_FacetGrid:
        """
        Apply a plotting function to a 1d facet's subset of the data.

        This is more convenient and less general than ``FacetGrid.map``

        Parameters
        ----------
        func :
            A plotting function with the same signature as a 1d xarray
            plotting method such as `xarray.plot.scatter`
        x, y :
            Names of the coordinates to plot on x, y axes
        **kwargs
            additional keyword arguments to func

        Returns
        -------
        self : FacetGrid object

        """
        # Copy data to allow converting categoricals to integers and storing
        # them in self.data. It is not possible to copy in the init
        # unfortunately as there are tests that relies on self.data being
        # mutable (test_names_appear_somewhere()). Maybe something to deprecate
        # not sure how much that is used outside these tests.
        self.data = self.data.copy()

        if kwargs.get("cbar_ax") is not None:
            raise ValueError("cbar_ax not supported by FacetGrid.")

        if func.__name__ == "scatter":
            size_ = kwargs.pop("_size", markersize)
            size_r = _MARKERSIZE_RANGE
        else:
            size_ = kwargs.pop("_size", linewidth)
            size_r = _LINEWIDTH_RANGE

        # Guess what coords to use if some of the values in coords_to_plot are None:
        coords_to_plot: MutableMapping[str, Hashable | None] = dict(
            x=x, z=z, hue=hue, size=size_
        )
        coords_to_plot = _guess_coords_to_plot(self.data, coords_to_plot, kwargs)

        # Handle hues:
        hue = coords_to_plot["hue"]
        hueplt = self.data.coords[hue] if hue else None  # TODO: _infer_line_data2 ?
        hueplt_norm = _Normalize(hueplt)
        self._hue_var = hueplt
        cbar_kwargs = kwargs.pop("cbar_kwargs", {})
        if hueplt_norm.data is not None:
            if not hueplt_norm.data_is_numeric:
                # TODO: Ticks seems a little too hardcoded, since it will always
                # show all the values. But maybe it's ok, since plotting hundreds
                # of categorical data isn't that meaningful anyway.
                cbar_kwargs.update(format=hueplt_norm.format, ticks=hueplt_norm.ticks)
                kwargs.update(levels=hueplt_norm.levels)

            cmap_params, cbar_kwargs = _process_cmap_cbar_kwargs(
                func,
                cast("DataArray", hueplt_norm.values).data,
                cbar_kwargs=cbar_kwargs,
                **kwargs,
            )
            self._cmap_extend = cmap_params.get("extend")
        else:
            cmap_params = {}

        # Handle sizes:
        size_ = coords_to_plot["size"]
        sizeplt = self.data.coords[size_] if size_ else None
        sizeplt_norm = _Normalize(data=sizeplt, width=size_r)
        if sizeplt_norm.data is not None:
            self.data[size_] = sizeplt_norm.values

        # Add kwargs that are sent to the plotting function, # order is important ???
        func_kwargs = {
            k: v
            for k, v in kwargs.items()
            if k not in {"cmap", "colors", "cbar_kwargs", "levels"}
        }
        func_kwargs.update(cmap_params)
        # Annotations will be handled later, skip those parts in the plotfunc:
        func_kwargs["add_colorbar"] = False
        func_kwargs["add_legend"] = False
        func_kwargs["add_title"] = False

        add_labels_ = np.zeros(self.axs.shape + (3,), dtype=bool)
        if kwargs.get("z") is not None:
            # 3d plots looks better with all labels. 3d plots can't sharex either so it
            # is easy to get lost while rotating the plots:
            add_labels_[:] = True
        else:
            # Subplots should have labels on the left and bottom edges only:
            add_labels_[-1, :, 0] = True  # x
            add_labels_[:, 0, 1] = True  # y
            # add_labels_[:, :, 2] = True  # z

        # Set up the lists of names for the row and column facet variables:
        if self._single_group:
            full = tuple(
                {self._single_group: x}
                for x in range(self.data[self._single_group].size)
            )
            empty = tuple(None for x in range(self._nrow * self._ncol - len(full)))
            name_d = full + empty
        else:
            rowcols = itertools.product(
                range(self.data[self._row_var].size),
                range(self.data[self._col_var].size),
            )
            name_d = tuple({self._row_var: r, self._col_var: c} for r, c in rowcols)
        name_dicts = np.array(name_d).reshape(self._nrow, self._ncol)

        # Plot the data for each subplot:
        for add_lbls, d, ax in zip(
            add_labels_.reshape((self.axs.size, -1)),
            name_dicts.flat,
            self.axs.flat,
            strict=True,
        ):
            func_kwargs["add_labels"] = add_lbls
            # None is the sentinel value
            if d is not None:
                subset = self.data.isel(d)
                mappable = func(
                    subset,
                    x=x,
                    y=y,
                    ax=ax,
                    hue=hue,
                    _size=size_,
                    **func_kwargs,
                    _is_facetgrid=True,
                )
                self._mappables.append(mappable)

        # Add titles and some touch ups:
        self._finalize_grid()
        self._set_lims()

        add_colorbar, add_legend = _determine_guide(
            hueplt_norm,
            sizeplt_norm,
            kwargs.get("add_colorbar"),
            kwargs.get("add_legend"),
            # kwargs.get("add_guide", None),
            # kwargs.get("hue_style", None),
        )

        if add_legend:
            use_legend_elements = func.__name__ != "hist"
            if use_legend_elements:
                self.add_legend(
                    use_legend_elements=use_legend_elements,
                    hueplt_norm=hueplt_norm if not add_colorbar else _Normalize(None),
                    sizeplt_norm=sizeplt_norm,
                    primitive=self._mappables,
                    legend_ax=self.fig,
                    plotfunc=func.__name__,
                )
            else:
                self.add_legend(use_legend_elements=use_legend_elements)

        if add_colorbar:
            # Colorbar is after legend so it correctly fits the plot:
            if "label" not in cbar_kwargs:
                cbar_kwargs["label"] = label_from_attrs(hueplt_norm.data)

            self.add_colorbar(**cbar_kwargs)

        return self

    def map_dataarray_line(
        self: T_FacetGrid,
        func: Callable,
        x: Hashable | None,
        y: Hashable | None,
        hue: Hashable | None,
        add_legend: bool = True,
        _labels=None,
        **kwargs: Any,
    ) -> T_FacetGrid:
        from xarray.plot.dataarray_plot import _infer_line_data

        for d, ax in zip(self.name_dicts.flat, self.axs.flat, strict=True):
            # None is the sentinel value
            if d is not None:
                subset = self.data.loc[d]
                mappable = func(
                    subset,
                    x=x,
                    y=y,
                    ax=ax,
                    hue=hue,
                    add_legend=False,
                    _labels=False,
                    **kwargs,
                )
                self._mappables.append(mappable)

        xplt, yplt, hueplt, huelabel = _infer_line_data(
            darray=self.data.loc[self.name_dicts.flat[0]], x=x, y=y, hue=hue
        )
        xlabel = label_from_attrs(xplt)
        ylabel = label_from_attrs(yplt)

        self._hue_var = hueplt
        self._finalize_grid(xlabel, ylabel)

        if add_legend and hueplt is not None and huelabel is not None:
            self.add_legend(label=huelabel)

        return self

    def map_dataset(
        self: T_FacetGrid,
        func: Callable,
        x: Hashable | None = None,
        y: Hashable | None = None,
        hue: Hashable | None = None,
        hue_style: HueStyleOptions = None,
        add_guide: bool | None = None,
        **kwargs: Any,
    ) -> T_FacetGrid:
        from xarray.plot.dataset_plot import _infer_meta_data

        kwargs["add_guide"] = False

        if kwargs.get("markersize"):
            kwargs["size_mapping"] = _parse_size(
                self.data[kwargs["markersize"]], kwargs.pop("size_norm", None)
            )

        meta_data = _infer_meta_data(
            self.data, x, y, hue, hue_style, add_guide, funcname=func.__name__
        )
        kwargs["meta_data"] = meta_data

        if hue and meta_data["hue_style"] == "continuous":
            cmap_params, cbar_kwargs = _process_cmap_cbar_kwargs(
                func, self.data[hue].to_numpy(), **kwargs
            )
            kwargs["meta_data"]["cmap_params"] = cmap_params
            kwargs["meta_data"]["cbar_kwargs"] = cbar_kwargs

        kwargs["_is_facetgrid"] = True

        if func.__name__ == "quiver" and "scale" not in kwargs:
            raise ValueError("Please provide scale.")
            # TODO: come up with an algorithm for reasonable scale choice

        for d, ax in zip(self.name_dicts.flat, self.axs.flat, strict=True):
            # None is the sentinel value
            if d is not None:
                subset = self.data.loc[d]
                maybe_mappable = func(
                    ds=subset, x=x, y=y, hue=hue, hue_style=hue_style, ax=ax, **kwargs
                )
                # TODO: this is needed to get legends to work.
                # but maybe_mappable is a list in that case :/
                self._mappables.append(maybe_mappable)

        self._finalize_grid(meta_data["xlabel"], meta_data["ylabel"])

        if hue:
            hue_label = meta_data.pop("hue_label", None)
            self._hue_label = hue_label
            if meta_data["add_legend"]:
                self._hue_var = meta_data["hue"]
                self.add_legend(label=hue_label)
            elif meta_data["add_colorbar"]:
                self.add_colorbar(label=hue_label, **cbar_kwargs)

        if meta_data["add_quiverkey"]:
            self.add_quiverkey(kwargs["u"], kwargs["v"])

        return self

    def _finalize_grid(self, *axlabels: Hashable) -> None:
        """Finalize the annotations and layout."""
        if not self._finalized:
            self.set_axis_labels(*axlabels)
            self.set_titles()
            self.fig.tight_layout()

            for ax, namedict in zip(self.axs.flat, self.name_dicts.flat, strict=True):
                if namedict is None:
                    ax.set_visible(False)

            self._finalized = True

    def _adjust_fig_for_guide(self, guide) -> None:
        # Draw the plot to set the bounding boxes correctly
        if hasattr(self.fig.canvas, "get_renderer"):
            renderer = self.fig.canvas.get_renderer()
        else:
            raise RuntimeError("MPL backend has no renderer")
        self.fig.draw(renderer)

        # Calculate and set the new width of the figure so the legend fits
        guide_width = guide.get_window_extent(renderer).width / self.fig.dpi
        figure_width = self.fig.get_figwidth()
        total_width = figure_width + guide_width
        self.fig.set_figwidth(total_width)

        # Draw the plot again to get the new transformations
        self.fig.draw(renderer)

        # Now calculate how much space we need on the right side
        guide_width = guide.get_window_extent(renderer).width / self.fig.dpi
        space_needed = guide_width / total_width + 0.02
        # margin = .01
        # _space_needed = margin + space_needed
        right = 1 - space_needed

        # Place the subplot axes to give space for the legend
        self.fig.subplots_adjust(right=right)

    def add_legend(
        self,
        *,
        label: str | None = None,
        use_legend_elements: bool = False,
        **kwargs: Any,
    ) -> None:
        if use_legend_elements:
            self.figlegend = _add_legend(**kwargs)
        else:
            assert self._hue_var is not None
            self.figlegend = self.fig.legend(
                handles=self._mappables[-1],
                labels=list(self._hue_var.to_numpy()),
                title=label if label is not None else label_from_attrs(self._hue_var),
                loc=kwargs.pop("loc", "center right"),
                **kwargs,
            )
        self._adjust_fig_for_guide(self.figlegend)

    def add_colorbar(self, **kwargs: Any) -> None:
        """Draw a colorbar."""
        kwargs = kwargs.copy()
        if self._cmap_extend is not None:
            kwargs.setdefault("extend", self._cmap_extend)
        # dont pass extend as kwarg if it is in the mappable
        if hasattr(self._mappables[-1], "extend"):
            kwargs.pop("extend", None)
        if "label" not in kwargs:
            from xarray import DataArray

            assert isinstance(self.data, DataArray)
            kwargs.setdefault("label", label_from_attrs(self.data))
        self.cbar = self.fig.colorbar(
            self._mappables[-1], ax=list(self.axs.flat), **kwargs
        )

    def add_quiverkey(self, u: Hashable, v: Hashable, **kwargs: Any) -> None:
        kwargs = kwargs.copy()

        magnitude = _get_nice_quiver_magnitude(self.data[u], self.data[v])
        units = self.data[u].attrs.get("units", "")
        self.quiverkey = self.axs.flat[-1].quiverkey(
            self._mappables[-1],
            X=0.8,
            Y=0.9,
            U=magnitude,
            label=f"{magnitude}\n{units}",
            labelpos="E",
            coordinates="figure",
        )

        # TODO: does not work because self.quiverkey.get_window_extent(renderer) = 0
        # https://github.com/matplotlib/matplotlib/issues/18530
        # self._adjust_fig_for_guide(self.quiverkey.text)

    def _get_largest_lims(self) -> dict[str, tuple[float, float]]:
        """
        Get largest limits in the facetgrid.

        Returns
        -------
        lims_largest : dict[str, tuple[float, float]]
            Dictionary with the largest limits along each axis.

        Examples
        --------
        >>> ds = xr.tutorial.scatter_example_dataset(seed=42)
        >>> fg = ds.plot.scatter(x="A", y="B", hue="y", row="x", col="w")
        >>> round(fg._get_largest_lims()["x"][0], 3)
        np.float64(-0.334)
        """
        lims_largest: dict[str, tuple[float, float]] = dict(
            x=(np.inf, -np.inf), y=(np.inf, -np.inf), z=(np.inf, -np.inf)
        )
        for axis in ("x", "y", "z"):
            # Find the plot with the largest xlim values:
            lower, upper = lims_largest[axis]
            for ax in self.axs.flat:
                get_lim: Callable[[], tuple[float, float]] | None = getattr(
                    ax, f"get_{axis}lim", None
                )
                if get_lim:
                    lower_new, upper_new = get_lim()
                    lower, upper = (min(lower, lower_new), max(upper, upper_new))
            lims_largest[axis] = (lower, upper)

        return lims_largest

    def _set_lims(
        self,
        x: tuple[float, float] | None = None,
        y: tuple[float, float] | None = None,
        z: tuple[float, float] | None = None,
    ) -> None:
        """
        Set the same limits for all the subplots in the facetgrid.

        Parameters
        ----------
        x : tuple[float, float] or None, optional
            x axis limits.
        y : tuple[float, float] or None, optional
            y axis limits.
        z : tuple[float, float] or None, optional
            z axis limits.

        Examples
        --------
        >>> ds = xr.tutorial.scatter_example_dataset(seed=42)
        >>> fg = ds.plot.scatter(x="A", y="B", hue="y", row="x", col="w")
        >>> fg._set_lims(x=(-0.3, 0.3), y=(0, 2), z=(0, 4))
        >>> fg.axs[0, 0].get_xlim(), fg.axs[0, 0].get_ylim()
        ((np.float64(-0.3), np.float64(0.3)), (np.float64(0.0), np.float64(2.0)))
        """
        lims_largest = self._get_largest_lims()

        # Set limits:
        for ax in self.axs.flat:
            for (axis, data_limit), parameter_limit in zip(
                lims_largest.items(), (x, y, z), strict=True
            ):
                set_lim = getattr(ax, f"set_{axis}lim", None)
                if set_lim:
                    set_lim(data_limit if parameter_limit is None else parameter_limit)

    def set_axis_labels(self, *axlabels: Hashable) -> None:
        """Set axis labels on the left column and bottom row of the grid."""
        from xarray.core.dataarray import DataArray

        for var, axis in zip(axlabels, ["x", "y", "z"], strict=False):
            if var is not None:
                if isinstance(var, DataArray):
                    getattr(self, f"set_{axis}labels")(label_from_attrs(var))
                else:
                    getattr(self, f"set_{axis}labels")(str(var))

    def _set_labels(
        self, axis: str, axes: Iterable, label: str | None = None, **kwargs
    ) -> None:
        if label is None:
            label = label_from_attrs(self.data[getattr(self, f"_{axis}_var")])
        for ax in axes:
            getattr(ax, f"set_{axis}label")(label, **kwargs)

    def set_xlabels(self, label: str | None = None, **kwargs: Any) -> None:
        """Label the x axis on the bottom row of the grid."""
        self._set_labels("x", self._bottom_axes, label, **kwargs)

    def set_ylabels(self, label: str | None = None, **kwargs: Any) -> None:
        """Label the y axis on the left column of the grid."""
        self._set_labels("y", self._left_axes, label, **kwargs)

    def set_zlabels(self, label: str | None = None, **kwargs: Any) -> None:
        """Label the z axis."""
        self._set_labels("z", self._left_axes, label, **kwargs)

    def set_titles(
        self,
        template: str = "{coord} = {value}",
        maxchar: int = 30,
        size=None,
        **kwargs,
    ) -> None:
        """
        Draw titles either above each facet or on the grid margins.

        Parameters
        ----------
        template : str, default: "{coord} = {value}"
            Template for plot titles containing {coord} and {value}
        maxchar : int, default: 30
            Truncate titles at maxchar
        **kwargs : keyword args
            additional arguments to matplotlib.text

        Returns
        -------
        self: FacetGrid object

        """
        import matplotlib as mpl

        if size is None:
            size = mpl.rcParams["axes.labelsize"]

        nicetitle = functools.partial(_nicetitle, maxchar=maxchar, template=template)

        if self._single_group:
            for d, ax in zip(self.name_dicts.flat, self.axs.flat, strict=True):
                # Only label the ones with data
                if d is not None:
                    coord, value = list(d.items()).pop()
                    title = nicetitle(coord, value)
                    ax.set_title(title, size=size, **kwargs)
        else:
            # The row titles on the right edge of the grid
            for index, (ax, row_name, handle) in enumerate(
                zip(self.axs[:, -1], self.row_names, self.row_labels, strict=True)
            ):
                title = nicetitle(coord=self._row_var, value=row_name)
                if not handle:
                    self.row_labels[index] = ax.annotate(
                        title,
                        xy=(1.02, 0.5),
                        xycoords="axes fraction",
                        rotation=270,
                        ha="left",
                        va="center",
                        **kwargs,
                    )
                else:
                    handle.set_text(title)
                    handle.update(kwargs)

            # The column titles on the top row
            for index, (ax, col_name, handle) in enumerate(
                zip(self.axs[0, :], self.col_names, self.col_labels, strict=True)
            ):
                title = nicetitle(coord=self._col_var, value=col_name)
                if not handle:
                    self.col_labels[index] = ax.set_title(title, size=size, **kwargs)
                else:
                    handle.set_text(title)
                    handle.update(kwargs)

    def set_ticks(
        self,
        max_xticks: int = _NTICKS,
        max_yticks: int = _NTICKS,
        fontsize: str | int = _FONTSIZE,
    ) -> None:
        """
        Set and control tick behavior.

        Parameters
        ----------
        max_xticks, max_yticks : int, optional
            Maximum number of labeled ticks to plot on x, y axes
        fontsize : string or int
            Font size as used by matplotlib text

        Returns
        -------
        self : FacetGrid object

        """
        from matplotlib.ticker import MaxNLocator

        # Both are necessary
        x_major_locator = MaxNLocator(nbins=max_xticks)
        y_major_locator = MaxNLocator(nbins=max_yticks)

        for ax in self.axs.flat:
            ax.xaxis.set_major_locator(x_major_locator)
            ax.yaxis.set_major_locator(y_major_locator)
            for tick in itertools.chain(
                ax.xaxis.get_major_ticks(), ax.yaxis.get_major_ticks()
            ):
                tick.label1.set_fontsize(fontsize)

    def map(
        self: T_FacetGrid, func: Callable, *args: Hashable, **kwargs: Any
    ) -> T_FacetGrid:
        """
        Apply a plotting function to each facet's subset of the data.

        Parameters
        ----------
        func : callable
            A plotting function that takes data and keyword arguments. It
            must plot to the currently active matplotlib Axes and take a
            `color` keyword argument. If faceting on the `hue` dimension,
            it must also take a `label` keyword argument.
        *args : Hashable
            Column names in self.data that identify variables with data to
            plot. The data for each variable is passed to `func` in the
            order the variables are specified in the call.
        **kwargs : keyword arguments
            All keyword arguments are passed to the plotting function.

        Returns
        -------
        self : FacetGrid object

        """
        import matplotlib.pyplot as plt

        for ax, namedict in zip(self.axs.flat, self.name_dicts.flat, strict=True):
            if namedict is not None:
                data = self.data.loc[namedict]
                plt.sca(ax)
                innerargs = [data[a].to_numpy() for a in args]
                maybe_mappable = func(*innerargs, **kwargs)
                # TODO: better way to verify that an artist is mappable?
                # https://stackoverflow.com/questions/33023036/is-it-possible-to-detect-if-a-matplotlib-artist-is-a-mappable-suitable-for-use-w#33023522
                if maybe_mappable and hasattr(maybe_mappable, "autoscale_None"):
                    self._mappables.append(maybe_mappable)

        self._finalize_grid(*args[:2])

        return self


def _easy_facetgrid(
    data: T_DataArrayOrSet,
    plotfunc: Callable,
    kind: Literal["line", "dataarray", "dataset", "plot1d"],
    x: Hashable | None = None,
    y: Hashable | None = None,
    row: Hashable | None = None,
    col: Hashable | None = None,
    col_wrap: int | None = None,
    sharex: bool = True,
    sharey: bool = True,
    aspect: float | None = None,
    size: float | None = None,
    subplot_kws: dict[str, Any] | None = None,
    ax: Axes | None = None,
    figsize: Iterable[float] | None = None,
    **kwargs: Any,
) -> FacetGrid[T_DataArrayOrSet]:
    """
    Convenience method to call xarray.plot.FacetGrid from 2d plotting methods

    kwargs are the arguments to 2d plotting method
    """
    if ax is not None:
        raise ValueError("Can't use axes when making faceted plots.")
    if aspect is None:
        aspect = 1
    if size is None:
        size = 3
    elif figsize is not None:
        raise ValueError("cannot provide both `figsize` and `size` arguments")
    if kwargs.get("z") is not None:
        # 3d plots doesn't support sharex, sharey, reset to mpl defaults:
        sharex = False
        sharey = False

    g = FacetGrid(
        data=data,
        col=col,
        row=row,
        col_wrap=col_wrap,
        sharex=sharex,
        sharey=sharey,
        figsize=figsize,
        aspect=aspect,
        size=size,
        subplot_kws=subplot_kws,
    )

    if kind == "line":
        return g.map_dataarray_line(plotfunc, x, y, **kwargs)

    if kind == "dataarray":
        return g.map_dataarray(plotfunc, x, y, **kwargs)

    if kind == "plot1d":
        return g.map_plot1d(plotfunc, x, y, **kwargs)

    if kind == "dataset":
        return g.map_dataset(plotfunc, x, y, **kwargs)

    raise ValueError(
        f"kind must be one of `line`, `dataarray`, `dataset` or `plot1d`, got {kind}"
    )
