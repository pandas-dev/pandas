# mypy: ignore-errors
"""Request parsers for the various AWS protocol specifications."""

from __future__ import annotations

import base64
import json
import re
from collections import OrderedDict
from collections.abc import MutableMapping
from datetime import datetime, timezone
from typing import TYPE_CHECKING, Any, TypedDict
from urllib.parse import unquote
from xml.etree import ElementTree as ETree
from xml.etree.ElementTree import ParseError as XMLParseError

from botocore import xform_name
from botocore.utils import is_json_value_header
from botocore.utils import parse_timestamp as botocore_parse_timestamp

if TYPE_CHECKING:
    from werkzeug.datastructures import Headers, MultiDict

    from moto.core.model import OperationModel


class RequestParserError(Exception):
    pass


def default_timestamp_parser(value: str) -> datetime:
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


def _text_content(func):
    # This decorator hides the difference between
    # an XML node with text or a plain string.  It's used
    # to ensure that scalar processing operates only on text
    # strings, which allows the same scalar handlers to be used
    # for XML nodes from the body and HTTP headers.
    def _get_text_content(self, shape, node_or_string):
        if hasattr(node_or_string, "text"):
            text = node_or_string.text
            if text is None:
                # If an XML node is empty <foo></foo>,
                # we want to parse that as an empty string,
                # not as a null/None value.
                text = ""
        else:
            text = node_or_string
        return func(self, shape, text)

    return _get_text_content


# Sentinel to signal the absence of a field in the input

UNDEFINED = object()


class RequestDict(TypedDict):
    body: str | bytes
    headers: Headers
    method: str
    url_params: dict[str, Any]
    url_path: str
    values: MultiDict


class RequestParser:
    DEFAULT_ENCODING = "utf-8"
    MAP_TYPE = dict

    def __init__(
        self,
        operation_model: OperationModel,
        timestamp_parser=None,
        blob_parser=None,
        map_type=None,
    ):
        self.operation_model = operation_model
        if timestamp_parser is None:
            timestamp_parser = default_timestamp_parser
        self._timestamp_parser = timestamp_parser
        if blob_parser is None:
            blob_parser = default_blob_parser
        self._blob_parser = blob_parser
        if map_type is not None:
            self.MAP_TYPE = map_type

    def parse(self, request_dict: RequestDict) -> dict[str, Any]:
        input_shape = self.operation_model.input_shape
        if input_shape is None:
            return {}
        parsed = self._do_parse(request_dict, input_shape)
        return parsed

    def _do_parse(self, request_dict, shape):
        raise NotImplementedError(f"{self.__class__.__name__}._do_parse")

    def _parse_shape(self, shape, node):
        handler = getattr(self, f"_parse_{shape.type_name}", self._default_handle)
        return handler(shape, node)

    def _parse_list(self, shape, node):
        # Enough implementations share list parsing that it's move up here to the base class.
        parsed = []
        member_shape = shape.member
        for item in node:
            parsed.append(self._parse_shape(member_shape, item))
        return parsed

    def _parse_map(self, shape, value):
        parsed = self.MAP_TYPE()
        key_shape = shape.key
        value_shape = shape.value
        for key, val in value.items():
            actual_key = self._parse_shape(key_shape, key)
            actual_value = self._parse_shape(value_shape, val)
            parsed[actual_key] = actual_value
        return parsed

    def _default_handle(self, shape, value):
        return value


class QueryParser(RequestParser):
    def _do_parse(self, request_dict, shape):
        parsed = self._parse_shape(shape, request_dict["values"])
        return parsed if parsed is not UNDEFINED else {}

    def _parse_shape(self, shape, node, prefix=""):
        handler = getattr(self, f"_parse_{shape.type_name}", self._default_handle)
        return handler(shape, node, prefix)

    def _gonna_recurse(self, query_params, prefix):
        if prefix == "":
            return False
        return not any(param_key.startswith(prefix) for param_key in query_params)

    def _parse_structure(self, shape, query_params, prefix=""):
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

    def _parse_list(self, shape, node, prefix=""):
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

    def _parse_map(self, shape, query_params, prefix=""):
        if self._is_shape_flattened(shape):
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

    def _parse_blob(self, shape, query_params, prefix=""):
        # Blob args must be base64 encoded.
        value = self._default_handle(shape, query_params, prefix)
        if value is UNDEFINED:
            return value
        return self._blob_parser(value)

    def _parse_timestamp(self, shape, query_params, prefix=""):
        value = self._default_handle(shape, query_params, prefix)
        if value is UNDEFINED:
            return value
        return self._timestamp_parser(value)

    def _parse_boolean(self, shape, query_params, prefix=""):
        value = self._default_handle(shape, query_params, prefix)
        if value is True or value is False:
            return value
        try:
            return value.lower() == "true"
        except AttributeError:
            pass
        return UNDEFINED

    def _parse_integer(self, shape, query_params, prefix=""):
        value = self._default_handle(shape, query_params, prefix)
        if value is UNDEFINED:
            return value
        return int(value)

    def _parse_float(self, shape, query_params, prefix=""):
        value = self._default_handle(shape, query_params, prefix)
        if value is UNDEFINED:
            return value
        return float(value)

    _parse_double = _parse_float
    _parse_long = _parse_integer

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

    def _is_shape_flattened(self, shape):
        return shape.serialization.get("flattened")


class EC2QueryParser(QueryParser):
    def _get_serialized_name(self, shape, default_name):
        if "queryName" in shape.serialization:
            return shape.serialization["queryName"]
        elif "name" in shape.serialization:
            # A locationName is always capitalized
            # on input for the ec2 protocol.
            name = shape.serialization["name"]
            return name[0].upper() + name[1:]
        else:
            return default_name

    def _parse_list(self, shape, node, prefix=""):
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


class BaseJSONParser(RequestParser):
    def _parse_structure(self, shape, value):
        if shape.is_document_type:
            final_parsed = value
        else:
            member_shapes = shape.members
            if value is None:
                # If the comes across the wire as "null" (None in python),
                # we should be returning this unchanged, instead of as an
                # empty dict.
                return None
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

    def _parse_blob(self, shape, value):
        return self._blob_parser(value)

    def _parse_timestamp(self, shape, value):
        return self._timestamp_parser(value)

    def _parse_float(self, shape, value):
        if value is UNDEFINED:
            return value
        return float(value)

    _parse_double = _parse_float

    def _parse_integer(self, _, value):
        if value is UNDEFINED:
            return value
        return int(value)

    _parse_byte = _parse_integer
    _parse_short = _parse_integer
    _parse_long = _parse_integer

    def _parse_body_as_json(self, body_contents):
        if not body_contents:
            return {}
        try:
            body = body_contents.decode(self.DEFAULT_ENCODING)
        except (UnicodeDecodeError, AttributeError):
            body = body_contents
        try:
            original_parsed = json.loads(body)
            return original_parsed
        except ValueError:
            # if the body cannot be parsed, include
            # the literal string as the message
            return {"message": body}


class JSONParser(BaseJSONParser):
    def _do_parse(self, request_dict, shape):
        parsed = self.MAP_TYPE()
        if shape is not None:
            parsed = self._parse_json_body(request_dict["body"], shape)
        return parsed

    def _parse_json_body(self, raw_body, shape):
        # The json.loads() gives us the primitive JSON types,
        # but we need to traverse the parsed JSON data to convert
        # to richer types (blobs, timestamps, etc.
        parsed_json = self._parse_body_as_json(raw_body)
        return self._parse_shape(shape, parsed_json)


class BaseRestParser(RequestParser):
    def _do_parse(self, request_dict, shape):
        final_parsed = {}
        self._add_modeled_parse(request_dict, shape, final_parsed)
        return final_parsed

    def _add_modeled_parse(self, request_dict, shape, final_parsed):
        if shape is None:
            return final_parsed
        member_shapes = shape.members
        self._parse_non_payload_attrs(request_dict, shape, member_shapes, final_parsed)
        self._parse_payload(request_dict, shape, member_shapes, final_parsed)

    def _parse_payload(self, response, shape, member_shapes, final_parsed):
        if "payload" in shape.serialization:
            # If a payload is specified in the output shape, then only that
            # shape is used for the body payload.
            payload_member_name = shape.serialization["payload"]
            body_shape = member_shapes[payload_member_name]
            if body_shape.type_name in ["string", "blob"]:
                # This is a stream
                body = response["body"]
                if isinstance(body, bytes):
                    body = body.decode(self.DEFAULT_ENCODING)
                if body != "":
                    final_parsed[payload_member_name] = body
            else:
                original_parsed = self._initial_body_parse(response["body"])
                value = self._parse_shape(body_shape, original_parsed)
                # Payload for empty dict is <foo /> for XML but not for JSON...
                # Need to utilize subclasses here...
                # may have to clean this up with UNDEFINED vs if value...
                # For now do this isinstance hack that only returns {} if body not empty!
                if value or (
                    response["body"] and value == {} and isinstance(self, RestXMLParser)
                ):
                    final_parsed[payload_member_name] = value
        else:
            original_parsed = self._initial_body_parse(response["body"])
            body_parsed = self._parse_shape(shape, original_parsed)
            final_parsed.update(body_parsed)

    def _parse_non_payload_attrs(self, response, shape, member_shapes, final_parsed):
        headers = response["headers"]
        for name in member_shapes:
            member_shape = member_shapes[name]
            location = member_shape.serialization.get("location")
            if location is None:
                continue
            elif location == "headers":
                final_parsed[name] = self._parse_header_map(member_shape, headers)
            elif location == "header":
                header_name = member_shape.serialization.get("name", name)
                if header_name in headers:
                    final_parsed[name] = self._parse_shape(
                        member_shape, headers[header_name]
                    )
            elif location == "uri":
                member_name = member_shape.serialization.get("name", name)
                uri_params = response["url_params"]
                value = uri_params.get(member_name)
                value = unquote(value)
                final_parsed[name] = self._parse_shape(member_shape, value)
            elif location == "querystring":
                qs = response["values"]
                member_name = member_shape.serialization.get("name", name)
                if member_shape.type_name == "list":
                    value = qs.getlist(member_name)
                elif member_shape.type_name == "map":
                    # Maintain querystring keys with multiple values as lists.
                    value = qs.to_dict(flat=False)
                    value = {k: v if len(v) > 1 else v[0] for k, v in value.items()}
                else:
                    value = qs.get(member_name, None)
                if value is None:
                    continue
                final_parsed[name] = self._parse_shape(member_shape, value)

    def _parse_header_map(self, shape, headers):
        # Note that headers are case-insensitive, so we .lower()
        # all header names and header prefixes.
        parsed = self.MAP_TYPE()
        prefix = shape.serialization.get("name", "").lower()
        for header_name in headers.keys():
            if header_name.lower().startswith(prefix):
                # The key name inserted into the parsed hash
                # strips off the prefix.
                name = header_name[len(prefix) :]
                parsed[name] = headers[header_name]
        return parsed

    def _initial_body_parse(self, body_contents):
        # This method should do the initial xml/json parsing of the
        # body.  We still need to walk the parsed body in order
        # to convert types, but this method will do the first round
        # of parsing.
        raise NotImplementedError("_initial_body_parse")

    def _parse_list(self, shape, node):
        location = shape.serialization.get("location")
        if location == "header" and not isinstance(node, list):
            # List in headers may be a comma separated string as per RFC7230
            node = [e.strip() for e in node.split(",")]
        return super()._parse_list(shape, node)


class RestJSONParser(BaseRestParser, BaseJSONParser):
    def _initial_body_parse(self, body_contents):
        return self._parse_body_as_json(body_contents)

    def _parse_string(self, shape, value):
        parsed = value
        if is_json_value_header(shape):
            decoded = base64.b64decode(value).decode(self.DEFAULT_ENCODING)
            parsed = json.loads(decoded)
        return parsed

    # Has to handle text from query string or JSON value
    def _parse_boolean(self, shape, value):
        if value is True or value is False:
            return value
        try:
            return value.lower() == "true"
        except AttributeError:
            pass
        return UNDEFINED


class BaseXMLParser(RequestParser):
    _namespace_re = re.compile("{.*}")

    def _parse_map(self, shape, node):
        parsed = {}
        key_shape = shape.key
        value_shape = shape.value
        key_location_name = key_shape.serialization.get("name") or "key"
        value_location_name = value_shape.serialization.get("name") or "value"
        if shape.serialization.get("flattened") and not isinstance(node, list):
            node = [node]
        for keyval_node in node:
            for single_pair in keyval_node:
                # Within each <entry> there's a <key> and a <value>
                tag_name = self._node_tag(single_pair)
                if tag_name == key_location_name:
                    key_name = self._parse_shape(key_shape, single_pair)
                elif tag_name == value_location_name:
                    val_name = self._parse_shape(value_shape, single_pair)
                else:
                    raise RequestParserError(f"Unknown tag: {tag_name}")
            parsed[key_name] = val_name
        return parsed

    def _node_tag(self, node):
        return self._namespace_re.sub("", node.tag)

    def _parse_list(self, shape, node):
        # When we use _build_name_to_xml_node, repeated elements are aggregated
        # into a list.  However, we can't tell the difference between a scalar
        # value and a single element flattened list.  So before calling the
        # real _parse_list, we know that "node" should actually be a list if
        # it's flattened, and if it's not, then we make it a one element list.
        if shape.serialization.get("flattened") and not isinstance(node, list):
            node = [node]
        return super()._parse_list(shape, node)

    def _parse_structure(self, shape, node):
        parsed = {}
        members = shape.members
        xml_dict = self._build_name_to_xml_node(node)
        for member_name in members:
            member_shape = members[member_name]
            if (
                "location" in member_shape.serialization
                or member_shape.serialization.get("eventheader")
            ):
                # All members with locations have already been handled,
                # so we don't need to parse these members.
                continue
            xml_name = self._member_key_name(member_shape, member_name)
            member_node = xml_dict.get(xml_name)
            if member_node is not None:
                parsed[member_name] = self._parse_shape(member_shape, member_node)
            elif member_shape.serialization.get("xmlAttribute"):
                attribs = {}
                location_name = member_shape.serialization["name"]
                for key, value in node.attrib.items():
                    new_key = self._namespace_re.sub(
                        location_name.split(":")[0] + ":", key
                    )
                    attribs[new_key] = value
                if location_name in attribs:
                    parsed[member_name] = attribs[location_name]
        return parsed

    def _member_key_name(self, shape, member_name):
        # This method is needed because we have to special case flattened list
        # with a serialization name.  If this is the case we use the
        # locationName from the list's member shape as the key name for the
        # surrounding structure.
        if shape.type_name == "list" and shape.serialization.get("flattened"):
            list_member_serialized_name = shape.member.serialization.get("name")
            if list_member_serialized_name is not None:
                return list_member_serialized_name
        serialized_name = shape.serialization.get("name")
        if serialized_name is not None:
            return serialized_name
        return member_name

    def _build_name_to_xml_node(self, parent_node):
        # If the parent node is actually a list. We should not be trying
        # to serialize it to a dictionary. Instead, return the first element
        # in the list.
        if isinstance(parent_node, list):
            return self._build_name_to_xml_node(parent_node[0])
        xml_dict = {}
        for item in parent_node:
            key = self._node_tag(item)
            if key in xml_dict:
                # If the key already exists, the most natural
                # way to handle this is to aggregate repeated
                # keys into a single list.
                # <foo>1</foo><foo>2</foo> -> {'foo': [Node(1), Node(2)]}
                if isinstance(xml_dict[key], list):
                    xml_dict[key].append(item)
                else:
                    # Convert from a scalar to a list.
                    xml_dict[key] = [xml_dict[key], item]
            else:
                xml_dict[key] = item
        return xml_dict

    def _parse_xml_string_to_dom(self, xml_string):
        try:
            parser = ETree.XMLParser(
                target=ETree.TreeBuilder(), encoding=self.DEFAULT_ENCODING
            )
            parser.feed(xml_string)
            root = parser.close()
        except XMLParseError as e:
            raise RequestParserError(
                f"Unable to parse response ({e}), "
                f"invalid XML received. Further retries may succeed:\n{xml_string}"
            )
        return root

    @_text_content
    def _parse_boolean(self, shape, value):
        if value == "true":
            return True
        else:
            return False

    @_text_content
    def _parse_float(self, shape, text):
        return float(text)

    @_text_content
    def _parse_timestamp(self, shape, text):
        return self._timestamp_parser(text)

    @_text_content
    def _parse_integer(self, shape, text):
        return int(text)

    @_text_content
    def _parse_string(self, shape, text):
        parsed = text
        # This if might be duplicated in JSON parser - can we consolidate?
        if is_json_value_header(shape):
            decoded = base64.b64decode(text).decode(self.DEFAULT_ENCODING)
            parsed = json.loads(decoded)
        return parsed

    @_text_content
    def _parse_blob(self, shape, text):
        return self._blob_parser(text)

    _parse_character = _parse_string
    _parse_double = _parse_float
    _parse_long = _parse_integer


class RestXMLParser(BaseRestParser, BaseXMLParser):
    def _initial_body_parse(self, body_contents):
        if not body_contents:
            return ETree.Element("")
        return self._parse_xml_string_to_dom(body_contents)

    def _parse_map(self, shape, node):
        if shape.serialization.get("location") == "querystring":
            return RequestParser._parse_map(self, shape, node)
        return super()._parse_map(shape, node)

    @_text_content
    def _parse_string(self, shape, text):
        text = super()._parse_string(shape, text)
        return text


PROTOCOL_PARSERS = {
    "ec2": EC2QueryParser,
    "json": JSONParser,
    "query": QueryParser,
    "rest-json": RestJSONParser,
    "rest-xml": RestXMLParser,
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
        self, data: dict[str, Any] | None = None, **kwargs: dict[str, Any]
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
