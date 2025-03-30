import json
import pathlib

import jsonschema

HERE = pathlib.Path(__file__).parent
SCHEMA_FILE = HERE / "schema.json"
SCHEMA = json.loads(SCHEMA_FILE.read_text(encoding="utf-8"))
SPEC_VERSION = SCHEMA["definitions"]["current-version"]["enum"][0]


def make_validator(key):
    """make a JSON Schema (Draft 7) validator"""
    schema = {"$ref": "#/definitions/{}".format(key)}
    schema.update(SCHEMA)
    return jsonschema.validators.Draft7Validator(schema)


SERVERS_RESPONSE = make_validator("servers-response")

LANGUAGE_SERVER_SPEC = make_validator("language-server-spec")

LANGUAGE_SERVER_SPEC_MAP = make_validator("language-server-specs-implementation-map")
