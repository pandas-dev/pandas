import json
import string
from typing import Any, Dict

import yaml

from moto.moto_api._internal import mock_random as random


def create_id() -> str:
    size = 10
    chars = list(range(10)) + list(string.ascii_lowercase)
    return "".join(str(random.choice(chars)) for x in range(size))


def deserialize_body(body: str) -> Dict[str, Any]:
    try:
        api_doc = json.loads(body)
    except json.JSONDecodeError:
        api_doc = yaml.safe_load(body)

    if "openapi" in api_doc or "swagger" in api_doc:
        return api_doc

    return {}


def to_path(prop: str) -> str:
    return "/" + prop
