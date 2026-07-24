"""
This module provides a facade over the `botocore.model` module.

It defines custom shape and service model classes that extend Botocore's
classes to add additional functionality and properties specific to Moto.

Because Botocore's shape classes are not designed for subclassing, this module
redefines the class hierarchy to ensure that all shape classes inherit from
both the relevant Botocore shape class as well as a new Moto base class.
This allows us to add common functionality to all shape types while still
retaining the original behavior of the Botocore classes without monkey-patching.

This module also defines a modified ShapeResolver to ensure that all shapes created
through the resolver are instances of our custom classes.

Any Moto code that needs to interact with Botocore's model classes should do so
through this facade to avoid direct dependencies on Botocore's internal structure,
which may change over time.
"""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any, cast

from botocore.model import ListShape as BotocoreListShape
from botocore.model import MapShape as BotocoreMapShape
from botocore.model import OperationModel as BotocoreOperationModel
from botocore.model import ServiceModel as BotocoreServiceModel
from botocore.model import Shape as BotocoreShape
from botocore.model import ShapeResolver as BotocoreShapeResolver
from botocore.model import StringShape as BotocoreStringShape
from botocore.model import StructureShape as BotocoreStructureShape
from botocore.utils import CachedProperty, instance_cache


class Shape(BotocoreShape):
    # Custom Moto model properties that we want available in the serialization dict.
    SERIALIZED_ATTRS = BotocoreShape.SERIALIZED_ATTRS + [
        "locationNameForQueryCompatibility"
    ]

    serialization: dict[str, Any]

    @property
    def is_flattened(self) -> bool:
        return self.serialization.get("flattened", False)

    @property
    def is_http_header_trait(self) -> bool:
        return self.serialization.get("location") in ["header", "headers"]

    @property
    def is_not_bound_to_body(self) -> bool:
        return "location" in self.serialization


class StringShape(BotocoreStringShape, Shape):
    pass


class ListShape(BotocoreListShape, Shape):
    pass


class MapShape(BotocoreMapShape, Shape):
    pass


class StructureShape(BotocoreStructureShape, Shape):
    @CachedProperty
    def members(self) -> dict[str, Shape]:  # type: ignore[override]
        return cast(dict[str, Shape], super().members)

    @CachedProperty
    def error_code_aliases(self) -> list[str]:
        if not self.metadata.get("exception", False):
            return []
        error_metadata = self.metadata.get("error", {})
        aliases = error_metadata.get("codeAliases", [])
        return aliases


class ServiceModel(BotocoreServiceModel):
    def __init__(
        self, service_description: Mapping[str, Any], service_name: str | None = None
    ):
        super().__init__(service_description, service_name)
        # Use our custom shape resolver.
        self._shape_resolver = ShapeResolver(service_description.get("shapes", {}))

    @instance_cache
    def operation_model(self, operation_name: str) -> OperationModel:  # type: ignore[misc]
        operation_model = super().operation_model(operation_name)
        model = getattr(operation_model, "_operation_model", {})
        return OperationModel(model, self, operation_name)

    @CachedProperty
    def is_query_compatible(self) -> bool:
        return "awsQueryCompatible" in self.metadata


class OperationModel(BotocoreOperationModel):
    _operation_model: dict[str, Any]
    _service_model: ServiceModel

    @property
    def service_model(self) -> ServiceModel:
        return self._service_model


class ShapeResolver(BotocoreShapeResolver):
    SHAPE_CLASSES = {
        "structure": StructureShape,  # type: ignore[dict-item]
        "list": ListShape,  # type: ignore[dict-item]
        "map": MapShape,  # type: ignore[dict-item]
        "string": StringShape,  # type: ignore[dict-item]
    }

    def get_shape_by_name(
        self, shape_name: str, member_traits: Mapping[str, Any] | None = None
    ) -> Shape:
        shape = super().get_shape_by_name(shape_name, member_traits)
        # The SHAPE_CLASSES dict only allows us to override Shape subclasses,
        # but Botocore will return its own Shape base class for unknown types.
        # If the returned shape is not a part of our hierarchy, we need to
        # convert it to our Shape base class here.
        if not isinstance(shape, Shape):
            shape = Shape(shape_name, getattr(shape, "_shape_model", {}), self)
        return shape
