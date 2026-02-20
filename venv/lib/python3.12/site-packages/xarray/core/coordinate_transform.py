from __future__ import annotations

from collections.abc import Hashable, Iterable, Mapping
from typing import Any, overload

import numpy as np


class CoordinateTransform:
    """Abstract coordinate transform with dimension & coordinate names.

    .. caution::
        This API is experimental and subject to change. Please report any bugs or surprising
        behaviour you encounter.
    """

    coord_names: tuple[Hashable, ...]
    dims: tuple[str, ...]
    dim_size: dict[str, int]
    dtype: Any

    def __init__(
        self,
        coord_names: Iterable[Hashable],
        dim_size: Mapping[str, int],
        dtype: Any = None,
    ):
        self.coord_names = tuple(coord_names)
        self.dims = tuple(dim_size)
        self.dim_size = dict(dim_size)

        if dtype is None:
            dtype = np.dtype(np.float64)
        self.dtype = dtype

    def forward(self, dim_positions: dict[str, Any]) -> dict[Hashable, Any]:
        """Perform grid -> world coordinate transformation.

        Parameters
        ----------
        dim_positions : dict
            Grid location(s) along each dimension (axis).

        Returns
        -------
        coord_labels : dict
            World coordinate labels.

        """
        # TODO: cache the results in order to avoid re-computing
        # all labels when accessing the values of each coordinate one at a time
        raise NotImplementedError

    def reverse(self, coord_labels: dict[Hashable, Any]) -> dict[str, Any]:
        """Perform world -> grid coordinate reverse transformation.

        Parameters
        ----------
        labels : dict
            World coordinate labels.

        Returns
        -------
        dim_positions : dict
            Grid relative location(s) along each dimension (axis).

        """
        raise NotImplementedError

    @overload
    def equals(self, other: CoordinateTransform) -> bool: ...

    @overload
    def equals(
        self, other: CoordinateTransform, *, exclude: frozenset[Hashable] | None = None
    ) -> bool: ...

    def equals(self, other: CoordinateTransform, **kwargs) -> bool:
        """Check equality with another CoordinateTransform of the same kind.

        Parameters
        ----------
        other : CoordinateTransform
            The other CoordinateTransform object to compare with this object.
        exclude : frozenset of hashable, optional
            Dimensions excluded from checking. It is None by default, (i.e.,
            when this method is not called in the context of alignment). For a
            n-dimensional transform this option allows a CoordinateTransform to
            optionally ignore any dimension in ``exclude`` when comparing
            ``self`` with ``other``. For a 1-dimensional transform this kwarg
            can be safely ignored, as this method is not called when all of the
            transform's dimensions are also excluded from alignment.
        """
        raise NotImplementedError

    def generate_coords(
        self, dims: tuple[str, ...] | None = None
    ) -> dict[Hashable, Any]:
        """Compute all coordinate labels at once."""
        if dims is None:
            dims = self.dims

        positions = np.meshgrid(
            *[np.arange(self.dim_size[d]) for d in dims],
            indexing="ij",
        )
        dim_positions = {dim: positions[i] for i, dim in enumerate(dims)}

        return self.forward(dim_positions)
