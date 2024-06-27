"""Event validators."""
from __future__ import annotations

import pathlib
from typing import Any

import jsonschema
from jsonschema import Draft7Validator, ValidationError
from referencing import Registry
from referencing.jsonschema import DRAFT7

from . import yaml

draft7_format_checker = (
    Draft7Validator.FORMAT_CHECKER
    if hasattr(Draft7Validator, "FORMAT_CHECKER")
    else jsonschema.draft7_format_checker
)


METASCHEMA_PATH = pathlib.Path(__file__).parent.joinpath("schemas")

EVENT_METASCHEMA_FILEPATH = METASCHEMA_PATH.joinpath("event-metaschema.yml")
EVENT_METASCHEMA = yaml.load(EVENT_METASCHEMA_FILEPATH)

EVENT_CORE_SCHEMA_FILEPATH = METASCHEMA_PATH.joinpath("event-core-schema.yml")
EVENT_CORE_SCHEMA = yaml.load(EVENT_CORE_SCHEMA_FILEPATH)

PROPERTY_METASCHEMA_FILEPATH = METASCHEMA_PATH.joinpath("property-metaschema.yml")
PROPERTY_METASCHEMA = yaml.load(PROPERTY_METASCHEMA_FILEPATH)

SCHEMA_STORE = {
    EVENT_METASCHEMA["$id"]: EVENT_METASCHEMA,
    PROPERTY_METASCHEMA["$id"]: PROPERTY_METASCHEMA,
    EVENT_CORE_SCHEMA["$id"]: EVENT_CORE_SCHEMA,
}

resources = [
    DRAFT7.create_resource(each)
    for each in (EVENT_METASCHEMA, PROPERTY_METASCHEMA, EVENT_CORE_SCHEMA)
]
METASCHEMA_REGISTRY: Registry[Any] = resources @ Registry()

JUPYTER_EVENTS_SCHEMA_VALIDATOR = Draft7Validator(
    schema=EVENT_METASCHEMA,
    registry=METASCHEMA_REGISTRY,
    format_checker=draft7_format_checker,
)

JUPYTER_EVENTS_CORE_VALIDATOR = Draft7Validator(
    schema=EVENT_CORE_SCHEMA,
    registry=METASCHEMA_REGISTRY,
    format_checker=draft7_format_checker,
)


def validate_schema(schema: dict[str, Any]) -> None:
    """Validate a schema dict."""
    try:
        # Validate the schema against Jupyter Events metaschema.
        JUPYTER_EVENTS_SCHEMA_VALIDATOR.validate(schema)
    except ValidationError as err:
        reserved_property_msg = " does not match '^(?!__.*)'"
        if reserved_property_msg in str(err):
            idx = str(err).find(reserved_property_msg)
            bad_property = str(err)[:idx].strip()
            msg = (
                f"{bad_property} is an invalid property name because it "
                "starts with `__`. Properties starting with 'dunder' "
                "are reserved as special meta-fields for Jupyter Events to use."
            )
            raise ValidationError(msg) from err
        raise err
