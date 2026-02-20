# mypy: ignore-errors
import base64
import json
from collections import OrderedDict
from collections.abc import MutableMapping
from datetime import datetime, timezone
from typing import Any, Optional

from botocore import xform_name
from botocore.utils import parse_timestamp as botocore_parse_timestamp

from moto.core.model import OperationModel

UNDEFINED = object()  # Sentinel to signal the absence of a field in the input


def parse_timestamp(value: str) -> datetime:
    """Parse a timestamp and return a naive datetime object in UTC.
    This matches Moto's internal representation of timestamps, based
    on moto.core.utils.utcnow().
    """
    parsed = botocore_parse_timestamp(value)
    as_naive_utc = parsed.astimezone(timezone.utc).replace(tzinfo=None)
    return as_naive_utc


def default_blob_parser(value):
    # We don't decode this to a str because it's possible that
    # the blob contains binary data that can't be decoded.
    return base64.b64decode(value)


class RequestParser:
    MAP_TYPE = dict

    def __init__(
        self,
        timestamp_parser=None,
        blob_parser=None,
        map_type=None,
    ):
        if timestamp_parser is None:
            timestamp_parser = parse_timestamp
        self._timestamp_parser = timestamp_parser
        if blob_parser is None:
            blob_parser = default_blob_parser
        self._blob_parser = blob_parser
        if map_type is not None:
            self.MAP_TYPE = map_type

    def parse(self, request_dict, operation_model: OperationModel):
        raise NotImplementedError()


class QueryParser(RequestParser):
    def parse(self, request_dict, operation_model: OperationModel):
        shape = operation_model.input_shape
        if shape is None:
            return {}
        parsed = self._do_parse(request_dict, shape)
        return parsed if parsed is not UNDEFINED else {}

    def _do_parse(self, request_dict, shape):
        parsed = self._parse_shape(shape, request_dict["query_params"])
        return parsed if parsed is not UNDEFINED else {}

    def _parse_shape(self, shape, node, prefix=""):
        handler = getattr(self, f"_handle_{shape.type_name}", self._default_handle)
        return handler(shape, node, prefix)

    def _gonna_recurse(self, query_params, prefix):
        if prefix == "":
            return False
        return not any(param_key.startswith(prefix) for param_key in query_params)

    def _handle_structure(self, shape, query_params, prefix=""):
        if self._gonna_recurse(query_params, prefix):
            return UNDEFINED
        parsed = self.MAP_TYPE()
        members = shape.members
        for member_name in members:
            member_shape = members[member_name]
            member_prefix = self._get_serialized_name(member_shape, member_name)
            if prefix:
                member_prefix = f"{prefix}.{member_prefix}"
            value = self._parse_shape(member_shape, query_params, member_prefix)
            parsed_key = self._parsed_key_name(member_name)
            if value is not UNDEFINED:
                parsed[parsed_key] = value
        return parsed if parsed != {} else UNDEFINED

    def _handle_list(self, shape, node, prefix=""):
        # The query protocol serializes empty lists as an empty string.
        if node.get(prefix, UNDEFINED) == "":
            return []

        if shape.is_flattened:
            list_prefix = prefix
            location_name = self._get_serialized_name(shape.member, None)
            if location_name is not None:
                list_prefix = ".".join(prefix.split(".")[:-1] + [location_name])
        else:
            list_name = self._get_serialized_name(shape.member, "member")
            list_prefix = f"{prefix}.{list_name}"
        parsed_list = []
        i = 1
        while True:
            element_name = f"{list_prefix}.{i}"
            element_shape = shape.member
            value = self._parse_shape(element_shape, node, element_name)
            if value is UNDEFINED:
                break
            parsed_list.append(value)
            i += 1
        return parsed_list if parsed_list != [] else UNDEFINED

    def _handle_map(self, shape, query_params, prefix=""):
        if shape.is_flattened:
            full_prefix = prefix
        else:
            full_prefix = f"{prefix}.entry"
        template = full_prefix + ".{i}.{suffix}"
        key_shape = shape.key
        value_shape = shape.value
        key_suffix = self._get_serialized_name(key_shape, default_name="key")
        value_suffix = self._get_serialized_name(value_shape, "value")
        parsed_map = self.MAP_TYPE()
        i = 1
        while True:
            key_name = template.format(i=i, suffix=key_suffix)
            value_name = template.format(i=i, suffix=value_suffix)
            key = self._parse_shape(key_shape, query_params, key_name)
            value = self._parse_shape(value_shape, query_params, value_name)
            if key is UNDEFINED:
                break
            parsed_map[key] = value
            i += 1
        return parsed_map if parsed_map != {} else UNDEFINED

    def _handle_timestamp(self, shape, query_params, prefix=""):
        value = self._default_handle(shape, query_params, prefix)
        return value if value is UNDEFINED else self._timestamp_parser(value)

    def _handle_blob(self, shape, query_params, prefix=""):
        value = self._default_handle(shape, query_params, prefix)
        if value is UNDEFINED:
            return value
        return self._blob_parser(value)

    def _handle_boolean(self, shape, query_params, prefix=""):
        value = self._default_handle(shape, query_params, prefix)
        try:
            return value.lower() == "true"
        except AttributeError:
            pass
        return UNDEFINED

    def _handle_integer(self, shape, query_params, prefix=""):
        value = self._default_handle(shape, query_params, prefix)
        return value if value is UNDEFINED else int(value)

    def _handle_float(self, shape, query_params, prefix=""):
        value = self._default_handle(shape, query_params, prefix)
        return value if value is UNDEFINED else float(value)

    _handle_double = _handle_float
    _handle_long = _handle_integer

    def _default_handle(self, shape, value, prefix=""):
        default_value = shape.metadata.get("default", UNDEFINED)
        return value.get(prefix, default_value)

    def _get_serialized_name(self, shape, default_name):
        serialized_name = shape.serialization.get("locationNameForQueryCompatibility")
        if serialized_name:
            return serialized_name
        return shape.serialization.get("name", default_name)

    def _parsed_key_name(self, member_name):
        key_name = member_name
        return key_name


class EC2QueryParser(QueryParser):
    def _get_serialized_name(self, shape, default_name):
        # https://smithy.io/2.0/aws/protocols/aws-ec2-query-protocol.html#query-key-resolution
        if "queryName" in shape.serialization:
            return shape.serialization["queryName"]
        if "name" in shape.serialization:
            name = shape.serialization["name"]
            return name[0].upper() + name[1:]
        return default_name

    def _handle_list(self, shape, node, prefix=""):
        parsed_list = []
        i = 1
        while True:
            element_name = f"{prefix}.{i}"
            element_shape = shape.member
            value = self._parse_shape(element_shape, node, element_name)
            if value is UNDEFINED:
                break
            parsed_list.append(value)
            i += 1
        return parsed_list if parsed_list != [] else UNDEFINED


class JSONParser(RequestParser):
    DEFAULT_ENCODING = "utf-8"

    def parse(self, request_dict, operation_model: OperationModel):
        shape = operation_model.input_shape
        parsed = self._do_parse(request_dict, shape)
        return parsed if parsed is not UNDEFINED else {}

    def _do_parse(self, request_dict, shape):
        parsed = self.MAP_TYPE()
        if shape is not None:
            parsed = self._handle_json_body(request_dict["body"], shape)
        return parsed

    def _handle_json_body(self, raw_body, shape):
        parsed_json = self._parse_body_as_json(raw_body)
        return self._parse_shape(shape, parsed_json)

    def _parse_body_as_json(self, body_contents):
        if not body_contents:
            return {}
        body = body_contents.decode(self.DEFAULT_ENCODING)
        original_parsed = json.loads(body)
        return original_parsed

    def _parse_shape(self, shape, node):
        handler = getattr(self, f"_handle_{shape.type_name}", self._default_handle)
        return handler(shape, node)

    def _default_handle(self, _, value):
        return value

    def _handle_blob(self, shape, value):
        value = self._default_handle(shape, value)
        if value is UNDEFINED:
            return value
        return self._blob_parser(value)

    def _handle_float(self, _, value):
        return float(value) if value is not UNDEFINED else value

    def _handle_list(self, shape, node):
        parsed = []
        member_shape = shape.member
        for item in node:
            parsed.append(self._parse_shape(member_shape, item))
        return parsed

    def _handle_map(self, shape, value):
        parsed = self.MAP_TYPE()
        key_shape = shape.key
        value_shape = shape.value
        for key, item_value in value.items():
            actual_key = self._parse_shape(key_shape, key)
            actual_value = self._parse_shape(value_shape, item_value)
            parsed[actual_key] = actual_value
        return parsed

    def _handle_structure(self, shape, value):
        member_shapes = shape.members
        final_parsed = self.MAP_TYPE()
        for member_name in member_shapes:
            member_shape = member_shapes[member_name]
            json_name = member_shape.serialization.get("name", member_name)
            raw_value = value.get(json_name)
            if raw_value is not None:
                final_parsed[member_name] = self._parse_shape(
                    member_shapes[member_name], raw_value
                )
        return final_parsed

    def _handle_timestamp(self, _, value):
        return self._timestamp_parser(value)

    _handle_double = _handle_float


PROTOCOL_PARSERS = {
    "ec2": EC2QueryParser,
    "json": JSONParser,
    "query": QueryParser,
}


class XFormedDict(MutableMapping[str, Any]):
    """
    A Pascal/Snake case-insensitive  ``dict``-like object.

        xfd = XFormedDict()
        xfd['DBInstanceIdentifier'] = 'identifier'
        xfd['DBInstanceIdentifier'] == 'identifier'  # True
        xfd['db_instance_identifier'] == 'identifier'  # True
        list(xfd) == ['db_instance_identifier']  # True

    ***Do Not Import*** This class is deprecated and will be removed in a future release.
    """

    def __init__(
        self, data: Optional[dict[str, Any]] = None, **kwargs: dict[str, Any]
    ) -> None:
        self._xform_cache: dict[str, Any] = {}
        self._store: MutableMapping[str, Any] = OrderedDict()
        if data is None:
            data = {}
        self.update(data, **kwargs)

    def _xformed(self, key: str) -> str:
        return xform_name(key, _xform_cache=self._xform_cache)

    def __setitem__(self, key: str, value: Any) -> None:
        # Use the xformed key for lookups, but store the actual key alongside the value.
        self._store[self._xformed(key)] = (key, value)

    def __getitem__(self, key: str) -> Any:
        return self._store[self._xformed(key)][1]

    def __delitem__(self, key: str) -> None:
        del self._store[self._xformed(key)]

    def __iter__(self) -> Any:
        return (key for key in self._store.keys())

    def __len__(self) -> int:
        return len(self._store)

    def original_dict(self) -> dict[str, Any]:
        original_dict = {}
        for _, keyval in self._store.items():
            key = keyval[0]
            value = keyval[1]
            if isinstance(value, XFormedDict):
                value = value.original_dict()
            if isinstance(value, list):
                value = [
                    v.original_dict() if isinstance(v, XFormedDict) else v
                    for v in value
                ]
            original_dict[key] = value
        return original_dict

    def __eq__(self, other: object) -> bool:
        if isinstance(other, dict):
            other = XFormedDict(other)
        else:
            return NotImplemented
        return dict(self.items()) == dict(other.items())
