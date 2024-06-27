"""Event schema objects."""
from __future__ import annotations

import json
from pathlib import Path, PurePath
from typing import Any, Dict, Union

from jsonschema import FormatChecker, validators
from referencing import Registry
from referencing.jsonschema import DRAFT7

try:
    from jsonschema.protocols import Validator
except ImportError:
    Validator = Any  # type:ignore[assignment, misc]

from . import yaml
from .validators import draft7_format_checker, validate_schema


class EventSchemaUnrecognized(Exception):
    """An error for an unrecognized event schema."""


class EventSchemaLoadingError(Exception):
    """An error for an event schema loading error."""


class EventSchemaFileAbsent(Exception):
    """An error for an absent event schema file."""


SchemaType = Union[Dict[str, Any], str, PurePath]


class EventSchema:
    """A validated schema that can be used.

    On instantiation, validate the schema against
    Jupyter Event's metaschema.

    Parameters
    ----------
    schema: dict or str
        JSON schema to validate against Jupyter Events.

    validator_class: jsonschema.validators
        The validator class from jsonschema used to validate instances
        of this event schema. The schema itself will be validated
        against Jupyter Event's metaschema to ensure that
        any schema registered here follows the expected form
        of Jupyter Events.

    registry:
        Registry for nested JSON schema references.
    """

    def __init__(
        self,
        schema: SchemaType,
        validator_class: type[Validator] = validators.Draft7Validator,  # type:ignore[assignment]
        format_checker: FormatChecker = draft7_format_checker,
        registry: Registry[Any] | None = None,
    ):
        """Initialize an event schema."""
        _schema = self._load_schema(schema)
        # Validate the schema against Jupyter Events metaschema.
        validate_schema(_schema)

        if registry is None:
            registry = DRAFT7.create_resource(_schema) @ Registry()

        # Create a validator for this schema
        self._validator = validator_class(_schema, registry=registry, format_checker=format_checker)  # type: ignore[call-arg]
        self._schema = _schema

    def __repr__(self) -> str:
        """A string repr for an event schema."""
        return json.dumps(self._schema, indent=2)

    @staticmethod
    def _ensure_yaml_loaded(schema: SchemaType, was_str: bool = False) -> None:
        """Ensures schema was correctly loaded into a dictionary. Raises
        EventSchemaLoadingError otherwise."""
        if isinstance(schema, dict):
            return

        error_msg = "Could not deserialize schema into a dictionary."

        def intended_as_path(schema: str) -> bool:
            path = Path(schema)
            return path.match("*.yml") or path.match("*.yaml") or path.match("*.json")

        # detect whether the user specified a string but intended a PurePath to
        # generate a more helpful error message
        if was_str and intended_as_path(schema):  # type:ignore[arg-type]
            error_msg += " Paths to schema files must be explicitly wrapped in a Pathlib object."
        else:
            error_msg += " Double check the schema and ensure it is in the proper form."

        raise EventSchemaLoadingError(error_msg)

    @staticmethod
    def _load_schema(schema: SchemaType) -> dict[str, Any]:
        """Load a JSON schema from different sources/data types.

        `schema` could be a dictionary or serialized string representing the
        schema itself or a Pathlib object representing a schema file on disk.

        Returns a dictionary with schema data.
        """

        # if schema is already a dictionary, return it
        if isinstance(schema, dict):
            return schema

        # if schema is PurePath, ensure file exists at path and then load from file
        if isinstance(schema, PurePath):
            if not Path(schema).exists():
                msg = f'Schema file not present at path "{schema}".'
                raise EventSchemaFileAbsent(msg)

            loaded_schema = yaml.load(schema)
            EventSchema._ensure_yaml_loaded(loaded_schema)
            return loaded_schema  # type:ignore[no-any-return]

        # finally, if schema is string, attempt to deserialize and return the output
        if isinstance(schema, str):
            # note the diff b/w load v.s. loads
            loaded_schema = yaml.loads(schema)
            EventSchema._ensure_yaml_loaded(loaded_schema, was_str=True)
            return loaded_schema  # type:ignore[no-any-return]

        msg = f"Expected a dictionary, string, or PurePath, but instead received {schema.__class__.__name__}."  # type:ignore[unreachable]
        raise EventSchemaUnrecognized(msg)

    @property
    def id(self) -> str:
        """Schema $id field."""
        return self._schema["$id"]  # type:ignore[no-any-return]

    @property
    def version(self) -> int:
        """Schema's version."""
        return self._schema["version"]  # type:ignore[no-any-return]

    @property
    def properties(self) -> dict[str, Any]:
        return self._schema["properties"]  # type:ignore[no-any-return]

    def validate(self, data: dict[str, Any]) -> None:
        """Validate an incoming instance of this event schema."""
        self._validator.validate(data)
