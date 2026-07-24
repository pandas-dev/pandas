from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass

from fontTools.designspaceLib import DesignSpaceDocument
from fontTools.ttLib.ttFont import TTFont
from fontTools.varLib.models import (
    VariationModel,
    noRound,
    normalizeValue,
    piecewiseLinearMap,
)

import typing
import warnings

if typing.TYPE_CHECKING:
    from typing import Self

LocationTuple = tuple[tuple[str, float], ...]
"""A hashable location."""


def Location(location: Mapping[str, float]) -> LocationTuple:
    """Create a hashable location from a dictionary-like location."""
    return tuple(sorted(location.items()))


class VariableScalar:
    """A scalar with different values at different points in the designspace."""

    values: dict[LocationTuple, int]
    """The values across various user-locations. Must always include the default
    location by time of building."""

    def __init__(self, location_value=None):
        self.values = {
            Location(location): value
            for location, value in (location_value or {}).items()
        }
        # Deprecated: only used by the add_to_variation_store() backwards-compat
        # shim. New code should use VariableScalarBuilder instead.
        self.axes = []

    def __repr__(self):
        items = []
        for location, value in self.values.items():
            loc = ",".join(
                [
                    f"{ax}={int(coord) if float(coord).is_integer() else coord}"
                    for ax, coord in location
                ]
            )
            items.append("%s:%i" % (loc, value))
        return "(" + (" ".join(items)) + ")"

    @property
    def does_vary(self) -> bool:
        values = list(self.values.values())
        return any(v != values[0] for v in values[1:])

    def add_value(self, location: Mapping[str, float], value: int):
        self.values[Location(location)] = value

    def add_to_variation_store(self, store_builder, model_cache=None, avar=None):
        """Deprecated: use VariableScalarBuilder.add_to_variation_store() instead."""
        warnings.warn(
            "VariableScalar.add_to_variation_store() is deprecated. "
            "Use VariableScalarBuilder.add_to_variation_store() instead.",
            DeprecationWarning,
            stacklevel=2,
        )
        if not self.axes:
            raise ValueError(
                ".axes must be defined on variable scalar before calling "
                "add_to_variation_store()"
            )
        builder = VariableScalarBuilder(
            axis_triples={
                ax.axisTag: (ax.minValue, ax.defaultValue, ax.maxValue)
                for ax in self.axes
            },
            axis_mappings=({} if avar is None else dict(avar.segments)),
            model_cache=model_cache if model_cache is not None else {},
        )
        return builder.add_to_variation_store(self, store_builder)


@dataclass
class VariableScalarBuilder:
    """A helper class for building variable scalars, or otherwise interrogating
    their variation model for interpolation or similar."""

    axis_triples: dict[str, tuple[float, float, float]]
    """Minimum, default, and maximum for each axis in user-coordinates."""
    axis_mappings: dict[str, Mapping[float, float]]
    """Optional mappings from normalized user-coordinates to normalized
    design-coordinates."""

    model_cache: dict[tuple[LocationTuple, ...], VariationModel]
    """We often use the same exact locations (i.e. font sources) for a large
    number of variable scalars. Instead of creating a model for each, cache
    them. Cache by user-location to avoid repeated mapping computations."""

    @classmethod
    def from_ttf(cls, ttf: TTFont) -> Self:
        return cls(
            axis_triples={
                axis.axisTag: (axis.minValue, axis.defaultValue, axis.maxValue)
                for axis in ttf["fvar"].axes
            },
            axis_mappings=(
                {}
                if (avar := ttf.get("avar")) is None
                else {axis: segments for axis, segments in avar.segments.items()}
            ),
            model_cache={},
        )

    @classmethod
    def from_designspace(cls, doc: DesignSpaceDocument) -> Self:
        return cls(
            axis_triples={
                axis.tag: (axis.minimum, axis.default, axis.maximum)
                for axis in doc.axes
            },
            axis_mappings={
                axis.tag: {
                    normalizeValue(
                        user, (axis.minimum, axis.default, axis.maximum)
                    ): normalizeValue(
                        design,
                        (
                            axis.map_forward(axis.minimum),
                            axis.map_forward(axis.default),
                            axis.map_forward(axis.maximum),
                        ),
                    )
                    for user, design in axis.map
                }
                for axis in doc.axes
                if axis.map
            },
            model_cache={},
        )

    def _fully_specify_location(self, location: LocationTuple) -> LocationTuple:
        """Validate and fully-specify a user-space location by filling in
        missing axes with their user-space defaults."""

        full = {}
        for axtag, value in location:
            if axtag not in self.axis_triples:
                raise ValueError("Unknown axis %s in %s" % (axtag, location))
            full[axtag] = value

        for axtag, (_, axis_default, _) in self.axis_triples.items():
            if axtag not in full:
                full[axtag] = axis_default

        return Location(full)

    def _normalize_location(self, location: LocationTuple) -> dict[str, float]:
        """Normalize a user-space location, applying avar mappings if present.

        TODO: This only handles avar1 (per-axis piecewise linear mappings),
        not avar2 (multi-dimensional mappings).
        """

        result = {}
        for axtag, value in location:
            axis_min, axis_default, axis_max = self.axis_triples[axtag]
            normalized = normalizeValue(value, (axis_min, axis_default, axis_max))
            mapping = self.axis_mappings.get(axtag)
            if mapping is not None:
                normalized = piecewiseLinearMap(normalized, mapping)
            result[axtag] = normalized

        return result

    def _full_locations_and_values(
        self, scalar: VariableScalar
    ) -> list[tuple[LocationTuple, int]]:
        """Return a list of (fully-specified user-space location, value) pairs,
        preserving order and length of scalar.values."""

        return [
            (self._fully_specify_location(loc), val)
            for loc, val in scalar.values.items()
        ]

    def default_value(self, scalar: VariableScalar) -> int:
        """Get the default value of a variable scalar."""

        default_loc = Location(
            {tag: default for tag, (_, default, _) in self.axis_triples.items()}
        )
        for location, value in self._full_locations_and_values(scalar):
            if location == default_loc:
                return value

        raise ValueError("Default value could not be found")

    def value_at_location(
        self, scalar: VariableScalar, location: LocationTuple
    ) -> float:
        """Interpolate the value of a scalar from a user-location."""

        location = self._fully_specify_location(location)
        pairs = self._full_locations_and_values(scalar)

        # If user location matches exactly, no axis mapping or variation model needed.
        for loc, val in pairs:
            if loc == location:
                return val

        values = [val for _, val in pairs]
        normalized_location = self._normalize_location(location)

        value = self.model(scalar).interpolateFromMasters(normalized_location, values)
        if value is None:
            raise ValueError("Insufficient number of values to interpolate")

        return value

    def model(self, scalar: VariableScalar) -> VariationModel:
        """Return a variation model based on a scalar's values.

        Variable scalars with the same fully-specified user-locations will use
        the same cached variation model."""

        pairs = self._full_locations_and_values(scalar)
        cache_key = tuple(loc for loc, _ in pairs)

        cached_model = self.model_cache.get(cache_key)
        if cached_model is not None:
            return cached_model

        normalized_locations = [self._normalize_location(loc) for loc, _ in pairs]
        axisOrder = list(self.axis_triples.keys())
        model = self.model_cache[cache_key] = VariationModel(
            normalized_locations, axisOrder=axisOrder
        )

        return model

    def get_deltas_and_supports(self, scalar: VariableScalar):
        """Calculate deltas and supports from this scalar's variation model."""
        values = list(scalar.values.values())
        return self.model(scalar).getDeltasAndSupports(values, round=round)

    def add_to_variation_store(
        self, scalar: VariableScalar, store_builder
    ) -> tuple[int, int]:
        """Serialize this scalar's variation model to a store, returning the
        default value and variation index."""

        deltas, supports = self.get_deltas_and_supports(scalar)
        store_builder.setSupports(supports)
        index = store_builder.storeDeltas(deltas, round=noRound)

        # NOTE: Default value should be an exact integer by construction of
        #       VariableScalar.
        return int(self.default_value(scalar)), index
