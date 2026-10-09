from __future__ import annotations

from typing import (
    TYPE_CHECKING,
    Any,
)

from matplotlib.colorbar import (
    Colorbar,
    make_axes,
    make_axes_gridspec,
)
from matplotlib.colors import (
    BoundaryNorm,
    NoNorm,
)
from matplotlib.ticker import AutoLocator

if TYPE_CHECKING:
    from matplotlib.axes import Axes
    from matplotlib.cm import ScalarMappable
    from matplotlib.colors import Norm


def is_singular_continuous_norm(norm: Norm) -> bool:
    """Return True when a continuous norm has equal, concrete limits."""
    return (
        norm.vmin is not None
        and norm.vmin == norm.vmax
        and not isinstance(norm, (BoundaryNorm, NoNorm))
    )


def _autoscale_unscaled_norm(mappable: ScalarMappable) -> None:
    # An unscaled replacement norm has not seen the data yet. Scale it the
    # same way a draw without a colorbar would, but do not re-enter
    # update_normal through the norm's "changed" callback.
    # No data to scale from. Matplotlib's colorbar skips autoscale in that
    # case; calling it here would raise TypeError.
    if mappable.norm.scaled() or mappable.get_array() is None:
        return
    with mappable.norm.callbacks.blocked():
        mappable.autoscale_None()


class _ConstantColorbar(Colorbar):
    def __init__(self, ax: Axes, mappable: ScalarMappable, **kwargs: Any) -> None:
        extend = kwargs.get("extend")
        if extend is None:
            extend = mappable.cmap.colorbar_extend
            if extend is False:
                extend = getattr(mappable.norm, "extend", "neither")
        super().__init__(
            ax, mappable, **kwargs, **self._get_swatch_kwds(mappable, extend)
        )

    @staticmethod
    def _get_swatch_kwds(
        mappable: ScalarMappable, extend: str | bool | None
    ) -> dict[str, Any]:
        _autoscale_unscaled_norm(mappable)
        if not is_singular_continuous_norm(mappable.norm):
            return {}
        vmin, vmax = mappable.get_clim()
        # GH 64980: pad only the display boundaries of a single-color swatch.
        # Supplying values keeps Colorbar from expanding the shared norm.
        lower, upper = AutoLocator().nonsingular(vmin, vmax)
        boundaries = [lower, upper]
        if extend in ("min", "both"):
            boundaries.insert(0, lower - (upper - lower))
        if extend in ("max", "both"):
            boundaries.append(upper + (upper - lower))
        return {
            "boundaries": boundaries,
            "values": [vmin] * (len(boundaries) - 1),
            "ticks": [vmin],
        }

    def update_normal(self, mappable: ScalarMappable | None = None) -> None:
        if mappable is None:
            mappable = self.mappable
        kwds = self._get_swatch_kwds(mappable, self.extend)
        old_values = self.values
        norm_changed = self.norm is not mappable.norm
        self.values = kwds.get("values")
        self.boundaries = kwds.get("boundaries")
        if (old_values is not None) != bool(kwds):
            # Invalidate the cached norm when the representation changes so
            # the base class also resets the axis scale and tick formatter.
            # The mappable's norm and its limits are left untouched.
            self.norm = None  # type: ignore[assignment]
        super().update_normal(mappable)
        if kwds and (norm_changed or old_values != self.values):
            self.set_ticks(kwds["ticks"])


def make_constant_colorbar(
    mappable: ScalarMappable, ax: Axes, **kwargs: Any
) -> Colorbar:
    # Used only for scatter plots whose continuous norm is initially singular.
    # Figure.colorbar always constructs matplotlib.colorbar.Colorbar, so this
    # copies its placement. update_normal handles later norm changes.
    fig = ax.get_figure(root=False)
    assert fig is not None
    current_ax = fig.gca()
    engine = fig.get_layout_engine()
    if ax.get_subplotspec() is not None and (
        engine is None or engine.colorbar_gridspec
    ):
        cax, kwargs = make_axes_gridspec(ax, **kwargs)
    else:
        cax, kwargs = make_axes(ax, **kwargs)
    fig.sca(current_ax)
    cax.grid(visible=False, which="both", axis="both")
    result = _ConstantColorbar(cax, mappable, **kwargs)
    fig.stale = True
    return result
