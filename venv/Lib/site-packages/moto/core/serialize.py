# mypy: disable-error-code="misc, override, var-annotated"
"""Response serializers for the various AWS protocol specifications.

There are some similarities among the different protocols with respect
to response serialization, so the code is structured in a way to avoid
code duplication where possible.  The diagram below illustrates the
inheritance hierarchy of the response serializer classes.

                           +--------------------+
                           | ResponseSerializer |
                           +--------------------+
                             ^       ^       ^
         +-------------------+       |       +---------------------+
         |                           |                             |
    +----+--------------+  +---------+----------+  +---------------+----+
    | BaseXMLSerializer |  | BaseRestSerializer |  | BaseJSONSerializer |
    +-------------------+  +--------------------+  +--------------------+
         ^             ^          ^           ^           ^           ^
         |             |          |           |           |           |
         |         +---+----------+----+  +---+-----------+----+      |
         |         | RestXMLSerializer |  | RestJSONSerializer |      |
         |         +-------------------+  +--------------------+      |
         |                                                            |
    +----+------------+                                  +------------+---+
    | QuerySerializer |                                  | JSONSerializer |
    +-----------------+                                  +----------------+
         ^
         |
    +----+----------+
    | EC2Serializer |
    +---------------+

Return Value
============

The response serializers expose a single public method: ``serialize()``
This method takes in any result (dict, object, Exception, etc.) and
returns a serialized ResponseDict of the following form:

    {
        "body": <RESPONSE_BODY>,
        "headers": <RESPONSE_HEADERS>,
        "status_code": <RESPONSE_HTTP_STATUS_CODE>,
    }

The body serialization output is text (Python strings), not bytes, mostly
for ease of inspection while debugging (pretty printing is also supported).
The exception is blob types, which are assumed to be binary and will be
encoded as ``utf-8``.

"""

from __future__ import annotations

import base64
import calendar
import json
from datetime import datetime
from typing import Any, Mapping, MutableMapping, Optional, TypedDict, Union

import xmltodict
from botocore import xform_name
from botocore.compat import formatdate
from botocore.model import (
    ListShape,
    MapShape,
    NoShapeFoundError,
    OperationModel,
    Shape,
    StructureShape,
)
from botocore.utils import is_json_value_header, parse_to_aware_datetime

Serialized = MutableMapping[str, Any]


class ResponseDict(TypedDict):
    body: str
    headers: MutableMapping[str, str]
    status_code: int


class SerializationContext:
    def __init__(self, request_id: Optional[str] = None) -> None:
        self.request_id = request_id or "request-id"


class ErrorShape(StructureShape):
    # Overriding super class property to keep mypy happy...
    @property
    def error_code(self) -> str:
        code = str(super().error_code)
        return code


class ShapeHelpersMixin:
    @staticmethod
    def get_serialized_name(shape: Shape, default_name: str) -> str:
        return shape.serialization.get("name", default_name)

    @staticmethod
    def is_flattened(shape: Shape) -> bool:
        return shape.serialization.get("flattened", False)

    @staticmethod
    def is_http_header_trait(shape: Shape) -> bool:
        return hasattr(shape, "serialization") and shape.serialization.get(
            "location"
        ) in ["header", "headers"]

    @staticmethod
    def is_not_bound_to_body(shape: Shape) -> bool:
        return hasattr(shape, "serialization") and "location" in shape.serialization


class TimestampSerializer:
    TIMESTAMP_FORMAT_ISO8601 = "iso8601"
    TIMESTAMP_FORMAT_RFC822 = "rfc822"
    TIMESTAMP_FORMAT_UNIX = "unixtimestamp"

    ISO8601 = "%Y-%m-%dT%H:%M:%SZ"
    ISO8601_MICRO = "%Y-%m-%dT%H:%M:%S.%fZ"

    def __init__(self, default_format: str) -> None:
        self.default_format = default_format

    def serialize(
        self, serialized: Serialized, value: Any, shape: Shape, key: str
    ) -> None:
        timestamp_format = shape.serialization.get(
            "timestampFormat", self.default_format
        )
        serialized_value = self._convert_timestamp_to_str(value, timestamp_format)
        serialized[key] = serialized_value

    def _timestamp_iso8601(self, value: datetime) -> str:
        if value.microsecond > 0:
            timestamp_format = self.ISO8601_MICRO
        else:
            timestamp_format = self.ISO8601
        return value.strftime(timestamp_format)

    @staticmethod
    def _timestamp_unixtimestamp(value: datetime) -> float:
        return int(calendar.timegm(value.timetuple()))

    def _timestamp_rfc822(self, value: Union[datetime, float]) -> str:
        if isinstance(value, datetime):
            value = self._timestamp_unixtimestamp(value)
        return formatdate(value, usegmt=True)

    def _convert_timestamp_to_str(
        self, value: Union[int, str, datetime], timestamp_format: str
    ) -> str:
        timestamp_format = timestamp_format.lower()
        converter = getattr(self, "_timestamp_%s" % timestamp_format)
        datetime_obj = parse_to_aware_datetime(value)  # type: ignore
        final_value = converter(datetime_obj)
        return final_value


class HeaderSerializer(ShapeHelpersMixin):
    # https://smithy.io/2.0/spec/http-bindings.html#httpheader-serialization-rules
    DEFAULT_ENCODING = "utf-8"
    DEFAULT_TIMESTAMP_FORMAT = TimestampSerializer.TIMESTAMP_FORMAT_RFC822

    def __init__(self, **kwargs: Mapping[str, Any]) -> None:
        super().__init__(**kwargs)
        self._timestamp_serializer = TimestampSerializer(self.DEFAULT_TIMESTAMP_FORMAT)

    def serialize(
        self, serialized: Serialized, value: Any, shape: Shape, key: str
    ) -> None:
        method = getattr(
            self, "_serialize_type_%s" % shape.type_name, self._default_serialize
        )
        method(serialized, value, shape, key)

    @staticmethod
    def _default_serialize(
        serialized: Serialized, value: Any, _: Shape, key: str
    ) -> None:
        serialized[key] = str(value)

    def _serialize_type_boolean(
        self, serialized: Serialized, value: Any, shape: Shape, key: str
    ) -> None:
        boolean_value = "true" if value else "false"
        self._default_serialize(serialized, boolean_value, shape, key)

    def _serialize_type_list(
        self, serialized: Serialized, value: Any, shape: ListShape, key: str
    ) -> None:
        list_value = ",".join(value)
        self._default_serialize(serialized, list_value, shape, key)

    def _serialize_type_map(
        self, serialized: Serialized, value: Any, shape: MapShape, _: str
    ) -> None:
        header_prefix = self.get_serialized_name(shape, "")
        for key, val in value.items():
            full_key = header_prefix + key
            self._default_serialize(serialized, val, shape, full_key)

    def _serialize_type_string(
        self, serialized: Serialized, value: Any, shape: Shape, key: str
    ) -> None:
        string_value = value
        if is_json_value_header(shape):
            json_value = json.dumps(value, separators=(",", ":"))
            string_value = self._base64(json_value)
        self._default_serialize(serialized, string_value, shape, key)

    def _serialize_type_timestamp(
        self, serialized: Serialized, value: Any, shape: Shape, key: str
    ) -> None:
        wrapper = {}
        self._timestamp_serializer.serialize(wrapper, value, shape, "timestamp")
        self._default_serialize(serialized, wrapper["timestamp"], shape, key)

    def _base64(self, value: Union[str, bytes]) -> str:
        if isinstance(value, str):
            value = value.encode(self.DEFAULT_ENCODING)
        return base64.b64encode(value).strip().decode(self.DEFAULT_ENCODING)


class ResponseSerializer(ShapeHelpersMixin):
    CONTENT_TYPE = "text"
    DEFAULT_ENCODING = "utf-8"
    DEFAULT_RESPONSE_CODE = 200
    DEFAULT_ERROR_RESPONSE_CODE = 400
    # From the spec, the default timestamp format if not specified is iso8601.
    DEFAULT_TIMESTAMP_FORMAT = TimestampSerializer.TIMESTAMP_FORMAT_ISO8601
    # Clients can change this to a different MutableMapping (i.e. OrderedDict) if they want.
    # This is used in the compliance test to match the hash ordering used in the tests.
    # NOTE: This is no longer necessary because dicts post 3.6 are ordered:
    # https://stackoverflow.com/questions/39980323/are-dictionaries-ordered-in-python-3-6
    MAP_TYPE = dict

    def __init__(
        self,
        operation_model: OperationModel,
        context: Optional[SerializationContext] = None,
        pretty_print: Optional[bool] = False,
        value_picker: Any = None,
    ) -> None:
        self.operation_model = operation_model
        self.context = context or SerializationContext()
        self.pretty_print = pretty_print
        if value_picker is None:
            value_picker = self._default_value_picker
        self._value_picker = value_picker
        self._timestamp_serializer = TimestampSerializer(self.DEFAULT_TIMESTAMP_FORMAT)

    def _create_default_response(self) -> ResponseDict:
        response_dict: ResponseDict = {
            "body": "",
            "headers": {},
            "status_code": self.DEFAULT_RESPONSE_CODE,
        }
        return response_dict

    def serialize(self, result: Any) -> ResponseDict:
        resp = self._create_default_response()
        if self._is_error_result(result):
            resp = self._serialize_error(resp, result)
        else:
            resp = self._serialize_result(resp, result)
        return resp

    def _serialize_error(
        self,
        resp: ResponseDict,
        error: Exception,
    ) -> ResponseDict:
        error_shape = self._get_error_shape(error)
        serialized_error = self.MAP_TYPE()
        self._serialize_error_metadata(serialized_error, error, error_shape)
        return self._serialized_error_to_response(
            resp, error, error_shape, serialized_error
        )

    def _serialize_result(self, resp: ResponseDict, result: Any) -> ResponseDict:
        output_shape = self.operation_model.output_shape
        serialized_result = self.MAP_TYPE()
        if output_shape is not None:
            assert isinstance(output_shape, StructureShape)  # mypy hint
            self._serialize(serialized_result, result, output_shape, "")
        return self._serialized_result_to_response(
            resp, result, output_shape, serialized_result
        )

    def _serialized_error_to_response(
        self,
        resp: ResponseDict,
        error: Exception,
        shape: ErrorShape,
        serialized_error: MutableMapping[str, Any],
    ) -> ResponseDict:
        raise NotImplementedError("Must be implemented in subclass.")

    def _serialized_result_to_response(
        self,
        resp: ResponseDict,
        result: Any,
        shape: Optional[StructureShape],
        serialized_result: MutableMapping[str, Any],
    ) -> ResponseDict:
        raise NotImplementedError("Must be implemented in subclass.")

    def _serialize_error_metadata(
        self,
        serialized: Serialized,
        error: Exception,
        shape: ErrorShape,
    ) -> None:
        raise NotImplementedError("Must be implemented in subclass.")

    def _serialize_body(self, body: Any) -> str:
        raise NotImplementedError("Must be implemented in subclass.")

    # Some extra utility methods subclasses can use.
    @staticmethod
    def _is_error_result(result: object) -> bool:
        return isinstance(result, Exception)

    def _base64(self, value: Union[str, bytes]) -> str:
        if isinstance(value, str):
            value = value.encode(self.DEFAULT_ENCODING)
        return base64.b64encode(value).strip().decode(self.DEFAULT_ENCODING)

    @staticmethod
    def _default_value_picker(obj: Any, key: str, _: Shape, default: Any = None) -> Any:
        if not hasattr(obj, "__getitem__"):
            return getattr(obj, key, default)

        try:
            return obj[key]
        except (KeyError, IndexError, TypeError, AttributeError):
            return getattr(obj, key, default)

    def _get_value(self, value: Any, key: str, shape: Shape) -> Any:
        return self._value_picker(value, key, shape)

    @staticmethod
    def _get_error_shape_name(error: Exception) -> str:
        shape_name = getattr(error, "code", error.__class__.__name__)
        return shape_name

    def _get_error_shape(self, error: Exception) -> ErrorShape:
        shape_name = self._get_error_shape_name(error)
        try:
            # TODO: there is also an errors array in the operation model,
            # but I think it only includes the possible errors for that
            # operation.  Maybe we try that first, then try all shapes?
            shape = self.operation_model.service_model.shape_for(shape_name)
            # We convert to ErrorShape to keep mypy happy...
            shape = ErrorShape(
                shape_name,
                shape._shape_model,  # type: ignore[attr-defined]
                shape._shape_resolver,  # type: ignore[attr-defined]
            )
        except NoShapeFoundError:
            generic_error_model = {
                "exception": True,
                "type": "structure",
                "members": {},
                "error": {
                    "code": shape_name,
                },
            }
            shape = ErrorShape(shape_name, generic_error_model)
        return shape

    #
    # Default serializers for the various model Shape types.
    # These can be overridden in subclasses to provide protocol-specific implementations.
    #
    def _serialize(
        self, serialized: Serialized, value: Any, shape: Shape, key: str
    ) -> None:
        method = getattr(
            self, "_serialize_type_%s" % shape.type_name, self._default_serialize
        )
        method(serialized, value, shape, key)

    @staticmethod
    def _default_serialize(
        serialized: Serialized, value: Any, _: Shape, key: str
    ) -> None:
        serialized[key] = value

    def _serialize_type_structure(
        self, serialized: Serialized, value: Any, shape: StructureShape, key: str
    ) -> None:
        if value is None:
            return
        if key:
            new_serialized = self.MAP_TYPE()
            serialized[key] = new_serialized
            serialized = new_serialized
        for member_key, member_shape in shape.members.items():
            self._serialize_structure_member(
                serialized, value, member_shape, member_key
            )

    def _serialize_structure_member(
        self, serialized: Serialized, value: Any, shape: Shape, key: str
    ) -> None:
        member_value = self._get_value(value, key, shape)
        if member_value is not None:
            key_name = self.get_serialized_name(shape, key)
            self._serialize(serialized, member_value, shape, key_name)

    def _serialize_type_map(
        self, serialized: Serialized, value: Any, shape: MapShape, key: str
    ) -> None:
        map_list = []
        if self.is_flattened(shape):
            items_name = key
            serialized[items_name] = map_list
        else:
            items_name = "entry"
            serialized[key] = {items_name: map_list}
        key_shape = shape.key
        assert isinstance(key_shape, Shape)
        value_shape = shape.value
        assert isinstance(value_shape, Shape)
        for key in value:
            wrapper = {"__current__": {}}
            key_prefix = self.get_serialized_name(key_shape, "key")
            value_prefix = self.get_serialized_name(value_shape, "value")
            self._serialize(wrapper["__current__"], key, key_shape, key_prefix)
            self._serialize(
                wrapper["__current__"], value[key], value_shape, value_prefix
            )
            map_list.append(wrapper["__current__"])

    def _serialize_type_timestamp(
        self, serialized: Serialized, value: Any, shape: Shape, key: str
    ) -> None:
        value_wrapper = {}
        value_key = "timestamp"
        self._timestamp_serializer.serialize(value_wrapper, value, shape, value_key)
        self._default_serialize(serialized, value_wrapper[value_key], shape, key)

    def _serialize_type_blob(
        self, serialized: Serialized, value: Any, shape: Shape, key: str
    ) -> None:
        blob_value = self._base64(value)
        self._default_serialize(serialized, blob_value, shape, key)


class BaseJSONSerializer(ResponseSerializer):
    APPLICATION_AMZ_JSON = "application/x-amz-json-{version}"
    DEFAULT_TIMESTAMP_FORMAT = "unixtimestamp"

    def _serialized_result_to_response(
        self,
        resp: ResponseDict,
        result: Any,
        shape: Optional[StructureShape],
        serialized_result: MutableMapping[str, Any],
    ) -> ResponseDict:
        resp["body"] = self._serialize_body(serialized_result)
        resp["headers"]["Content-Type"] = self._get_protocol_specific_content_type()
        return resp

    def _serialized_error_to_response(
        self,
        resp: ResponseDict,
        error: Exception,
        shape: ErrorShape,
        serialized_error: MutableMapping[str, Any],
    ) -> ResponseDict:
        resp["body"] = self._serialize_body(serialized_error)
        status_code = shape.metadata.get("error", {}).get(
            "httpStatusCode", self.DEFAULT_ERROR_RESPONSE_CODE
        )
        resp["status_code"] = status_code
        error_code = self._get_protocol_specific_error_code(shape.error_code)
        resp["headers"]["X-Amzn-Errortype"] = error_code
        resp["headers"]["Content-Type"] = self._get_protocol_specific_content_type()
        return resp

    def _get_protocol_specific_content_type(self) -> str:
        content_type = self.CONTENT_TYPE
        service_model = self.operation_model.service_model
        protocol = service_model.protocol
        if protocol == "json":
            json_version = service_model.metadata.get("jsonVersion", "1.0")
            content_type = self.APPLICATION_AMZ_JSON.format(version=json_version)
        return content_type

    def _get_protocol_specific_error_code(
        self,
        error_code: str,
    ) -> str:
        # https://smithy.io/2.0/aws/protocols/aws-json-1_1-protocol.html#operation-error-serialization
        service_metadata = self.operation_model.service_model.metadata
        json_version = service_metadata.get("jsonVersion")
        prefix = service_metadata.get("targetPrefix")
        if json_version == "1.0" and prefix is not None:
            error_code = prefix + "#" + error_code
        return error_code

    def _serialize_error_metadata(
        self,
        serialized: Serialized,
        error: Exception,
        shape: ErrorShape,
    ) -> None:
        error_code = self._get_protocol_specific_error_code(shape.error_code)
        serialized["__type"] = error_code
        message = getattr(error, "message", None) or str(error)
        if shape is not None:
            self._serialize(serialized, error, shape, "")
        if message:
            serialized["Message"] = message

    def _serialize_body(self, body: Mapping[str, Any]) -> str:
        body_encoded = json.dumps(body, indent=4 if self.pretty_print else None)
        return body_encoded

    def _serialize_type_map(
        self, serialized: Serialized, value: Any, shape: MapShape, key: str
    ) -> None:
        map_obj = self.MAP_TYPE()
        serialized[key] = map_obj
        for sub_key, sub_value in value.items():
            assert isinstance(shape.value, Shape)  # mypy hint
            self._serialize(map_obj, sub_value, shape.value, sub_key)

    def _serialize_type_list(
        self, serialized: Serialized, value: Any, shape: ListShape, key: str
    ) -> None:
        list_obj = []
        serialized[key] = list_obj
        for list_item in value:
            wrapper = {}
            # The JSON list serialization is the only case where we aren't
            # setting a key on a dict.  We handle this by using
            # a __current__ key on a wrapper dict to serialize each
            # list item before appending it to the serialized list.
            assert isinstance(shape.member, Shape)  # mypy hint
            self._serialize(wrapper, list_item, shape.member, "__current__")
            if "__current__" in wrapper:
                list_obj.append(wrapper["__current__"])
            else:
                list_obj.append(list_item)

    def _serialize_type_structure(
        self, serialized: Serialized, value: Any, shape: StructureShape, key: str
    ) -> None:
        if shape.is_document_type:
            serialized[key] = value
            return
        super()._serialize_type_structure(serialized, value, shape, key)


class BaseXMLSerializer(ResponseSerializer):
    CONTENT_TYPE = "text/xml"

    def _serialize_namespace_attribute(self, serialized: Serialized) -> None:
        if (
            self.CONTENT_TYPE == "text/xml"
            and "xmlNamespace" in self.operation_model.metadata
        ):
            namespace = self.operation_model.metadata["xmlNamespace"]
            serialized["@xmlns"] = namespace

    def _serialized_error_to_response(
        self,
        resp: ResponseDict,
        error: Exception,
        shape: ErrorShape,
        serialized_error: MutableMapping[str, Any],
    ) -> ResponseDict:
        error_wrapper = {
            "ErrorResponse": {
                "Error": serialized_error,
                "RequestId": self.context.request_id,
            }
        }
        self._serialize_namespace_attribute(error_wrapper["ErrorResponse"])
        resp["body"] = self._serialize_body(error_wrapper)
        status_code = shape.metadata.get("error", {}).get(
            "httpStatusCode", self.DEFAULT_ERROR_RESPONSE_CODE
        )
        resp["status_code"] = status_code
        resp["headers"]["Content-Type"] = self.CONTENT_TYPE
        return resp

    def _serialized_result_to_response(
        self,
        resp: ResponseDict,
        result: Any,
        shape: Optional[StructureShape],
        serialized_result: MutableMapping[str, Any],
    ) -> ResponseDict:
        result_key = f"{self.operation_model.name}Result"
        result_wrapper = {
            result_key: serialized_result,
        }
        self._serialize_namespace_attribute(result_wrapper[result_key])
        resp["body"] = self._serialize_body(result_wrapper)
        resp["headers"]["Content-Type"] = self.CONTENT_TYPE
        return resp

    def _serialize_error_metadata(
        self,
        serialized: Serialized,
        error: Exception,
        shape: ErrorShape,
    ) -> None:
        sender_fault = shape.metadata.get("error", {}).get("senderFault", True)
        serialized["Type"] = "Sender" if sender_fault else "Receiver"
        serialized["Code"] = shape.error_code
        message = getattr(error, "message", None)
        if message is not None:
            serialized["Message"] = message
        # Serialize any error model attributes.
        self._serialize(serialized, error, shape, "")

    def _serialize_body(self, body: Serialized) -> str:
        body_encoded = xmltodict.unparse(
            body,
            full_document=False,
            short_empty_elements=True,
            pretty=self.pretty_print,
        )
        return body_encoded

    #
    # https://smithy.io/2.0/aws/protocols/aws-query-protocol.html#xml-shape-serialization
    #
    def _serialize_type_boolean(
        self, serialized: Serialized, value: Any, shape: Shape, key: str
    ) -> None:
        # We're slightly more permissive here than we should be because the
        # moto backends are not consistent in how they store boolean values.
        # TODO: This should eventually be turned into a strict `is True` check.
        boolean_conditions = [
            (value is True),
            (str(value).lower() == "true"),
        ]
        boolean_value = "true" if any(boolean_conditions) else "false"
        self._default_serialize(serialized, boolean_value, shape, key)

    def _serialize_type_integer(
        self, serialized: Serialized, value: Any, shape: Shape, key: str
    ) -> None:
        integer_value = int(value)
        self._default_serialize(serialized, integer_value, shape, key)

    def _serialize_type_list(
        self, serialized: Serialized, value: Any, shape: ListShape, key: str
    ) -> None:
        assert isinstance(shape.member, Shape)  # mypy hinting
        list_obj = []
        if self.is_flattened(shape):
            items_name = self.get_serialized_name(shape.member, key)
            serialized[items_name] = list_obj
        else:
            items_name = self.get_serialized_name(shape.member, "member")
            serialized[key] = {items_name: list_obj}
        for list_item in value:
            wrapper = {}
            self._serialize(wrapper, list_item, shape.member, "__current__")
            list_obj.append(wrapper["__current__"])
        if not list_obj:
            serialized[key] = ""

    _serialize_type_long = _serialize_type_integer

    def _serialize_type_string(
        self, serialized: Serialized, value: Any, shape: Shape, key: str
    ) -> None:
        string_value = str(value)
        self._default_serialize(serialized, string_value, shape, key)


class BaseRestSerializer(ResponseSerializer):
    EMPTY_BODY: Serialized = ResponseSerializer.MAP_TYPE()
    REQUIRES_EMPTY_BODY = False

    def _serialized_result_to_response(
        self,
        resp: ResponseDict,
        result: Any,
        shape: Optional[StructureShape],
        serialized_result: MutableMapping[str, Any],
    ) -> ResponseDict:
        if "payload" in serialized_result:
            # Payload trumps all and is delivered as-is.
            resp["body"] = serialized_result["payload"]
        else:
            if not serialized_result["body"]:
                if self.REQUIRES_EMPTY_BODY:
                    resp["body"] = self._serialize_body(self.EMPTY_BODY)
            else:
                resp = super()._serialized_result_to_response(
                    resp, result, shape, serialized_result.get("body", {})
                )
        if "headers" in serialized_result:
            resp["headers"].update(serialized_result["headers"])
        resp["headers"]["Content-Type"] = self.CONTENT_TYPE
        return resp

    def _serialize_result(self, resp: ResponseDict, result: Any) -> ResponseDict:
        output_shape = self.operation_model.output_shape
        serialized_result = {
            "body": {},
            "headers": {},
        }
        if output_shape is not None:
            assert isinstance(output_shape, StructureShape)
            self._serialize(serialized_result, result, output_shape, "")
            payload_member = output_shape.serialization.get("payload")
            if payload_member is not None:
                payload_shape = output_shape.members[payload_member]
                payload_value = self._get_value(result, payload_member, payload_shape)
                self._serialize_payload(serialized_result, payload_value, payload_shape)

        return self._serialized_result_to_response(
            resp, result, output_shape, serialized_result
        )

    def _serialize_payload(
        self,
        serialized: Serialized,
        payload: Any,
        payload_shape: Shape,
    ) -> None:
        if payload_shape.type_name in ["blob", "string"]:
            # If it's streaming, then the body is just the value of the payload.
            serialized["payload"] = payload
        else:
            # If there's a payload member, we serialize only that one member to the body.
            serialized["body"] = self.MAP_TYPE()
            self._serialize(serialized["body"], payload, payload_shape, "")

    def _serialize_structure_member(
        self, serialized: Serialized, value: Any, shape: Shape, key: str
    ) -> None:
        if self.is_not_bound_to_body(shape):
            if self.is_http_header_trait(shape) and "headers" in serialized:
                member_value = self._get_value(value, key, shape)
                if member_value is not None:
                    key_name = self.get_serialized_name(shape, key)
                    header_serializer = HeaderSerializer()
                    header_serializer.serialize(
                        serialized["headers"], member_value, shape, key_name
                    )
        elif "body" in serialized:
            if not serialized["body"]:
                serialized["body"] = self.MAP_TYPE()
            # we're at the top-level structure
            super()._serialize_structure_member(serialized["body"], value, shape, key)
        else:
            # we're in nested structure
            super()._serialize_structure_member(serialized, value, shape, key)


class RestXMLSerializer(BaseRestSerializer, BaseXMLSerializer):
    DEFAULT_TIMESTAMP_FORMAT = TimestampSerializer.TIMESTAMP_FORMAT_ISO8601


class RestJSONSerializer(BaseRestSerializer, BaseJSONSerializer):
    CONTENT_TYPE = "application/json"
    REQUIRES_EMPTY_BODY = True


class JSONSerializer(BaseJSONSerializer):
    pass


class QuerySerializer(BaseXMLSerializer):
    def _serialized_error_to_response(
        self,
        resp: ResponseDict,
        error: Exception,
        shape: ErrorShape,
        serialized_error: MutableMapping[str, Any],
    ) -> ResponseDict:
        error_wrapper = {
            "ErrorResponse": {
                "Error": serialized_error,
                "RequestId": self.context.request_id,
            }
        }
        self._serialize_namespace_attribute(error_wrapper["ErrorResponse"])
        resp["body"] = self._serialize_body(error_wrapper)
        status_code = shape.metadata.get("error", {}).get(
            "httpStatusCode", self.DEFAULT_ERROR_RESPONSE_CODE
        )
        resp["status_code"] = status_code
        resp["headers"]["Content-Type"] = self.CONTENT_TYPE
        return resp

    def _serialized_result_to_response(
        self,
        resp: ResponseDict,
        result: Any,
        shape: Optional[StructureShape],
        serialized_result: MutableMapping[str, Any],
    ) -> ResponseDict:
        response_key = f"{self.operation_model.name}Response"
        response_wrapper = {response_key: {}}
        if shape is not None:
            result_key = shape.serialization.get("resultWrapper", f"{shape.name}Result")
            response_wrapper[response_key][result_key] = serialized_result
        response_wrapper[response_key]["ResponseMetadata"] = {
            "RequestId": self.context.request_id
        }
        self._serialize_namespace_attribute(response_wrapper[response_key])
        resp["body"] = self._serialize_body(response_wrapper)
        resp["headers"]["Content-Type"] = self.CONTENT_TYPE
        return resp


class EC2Serializer(QuerySerializer):
    def _serialize_body(self, body: Mapping[str, Any]) -> str:
        body_serialized = xmltodict.unparse(
            body,
            full_document=True,
            pretty=self.pretty_print,
            short_empty_elements=True,
        )
        return body_serialized

    def _serialize_error_metadata(
        self,
        serialized: MutableMapping[str, Any],
        error: Exception,
        shape: ErrorShape,
    ) -> None:
        serialized["Code"] = shape.error_code
        message = getattr(error, "message", None)
        if message is not None:
            serialized["Message"] = message
        # Serialize any error model attributes.
        self._serialize(serialized, error, shape, "")

    def _serialized_error_to_response(
        self,
        resp: ResponseDict,
        error: Exception,
        shape: ErrorShape,
        serialized_error: MutableMapping[str, Any],
    ) -> ResponseDict:
        error_wrapper = {
            "Response": {
                "Errors": [{"Error": serialized_error}],
                "RequestID": self.context.request_id,
            }
        }
        self._serialize_namespace_attribute(error_wrapper["Response"])
        resp["body"] = self._serialize_body(error_wrapper)
        status_code = shape.metadata.get("error", {}).get(
            "httpStatusCode", self.DEFAULT_ERROR_RESPONSE_CODE
        )
        resp["status_code"] = status_code
        resp["headers"]["Content-Type"] = self.CONTENT_TYPE
        return resp

    def _serialized_result_to_response(
        self,
        resp: ResponseDict,
        result: Any,
        shape: StructureShape,
        serialized_result: MutableMapping[str, Any],
    ) -> ResponseDict:
        response_key = f"{self.operation_model.name}Response"
        result_wrapper = {
            response_key: serialized_result,
        }
        result_wrapper[response_key]["requestId"] = self.context.request_id
        self._serialize_namespace_attribute(result_wrapper[response_key])
        resp["body"] = self._serialize_body(result_wrapper)
        resp["headers"]["Content-Type"] = self.CONTENT_TYPE
        return resp


SERIALIZERS = {
    "ec2": EC2Serializer,
    "json": JSONSerializer,
    "query": QuerySerializer,
    "rest-json": RestJSONSerializer,
    "rest-xml": RestXMLSerializer,
}


class XFormedAttributePicker:
    """Can be injected into a ResponseSerializer to aid in plucking AWS model
    attributes specified in `camelCase` or `PascalCase` from Python objects
    with standard `snake_case` attribute names.

    For a model attribute named `DBInstanceIdentifier`, this class will check
    for the following attributes on the provided object:
       * `DBInstanceIdentifier`
       * `db_instance_identifier`
    If the provided object is a class named `DBInstance`, this class will also
    check for the following attribute on the provided object:
       * `identifier`

    Uses ``botocore.xform_name`` to translate the attribute name.
    """

    def __init__(self) -> None:
        self._xform_cache = {}

    def __call__(self, value: Any, key: str, shape: Shape) -> Any:
        return self._get_value(value, key, shape)

    def xform_name(self, name: str) -> str:
        return xform_name(name, _xform_cache=self._xform_cache)

    def _get_value(self, value: Any, key: str, shape: Shape) -> Any:
        new_value = None
        possible_keys = [key, self.xform_name(key)]
        # If a class `Role` has an attribute named `arn`, that will work for a `RoleArn` key.
        if hasattr(value, "__class__"):
            class_name = value.__class__.__name__
            if key.lower().startswith(class_name.lower()):
                short_key = key[len(class_name) :]
                if short_key:  # Will be empty string if class name same as key
                    possible_keys.append(self.xform_name(short_key))
        if shape is not None:
            serialization_key = shape.serialization.get("name", key)
            if serialization_key != key:
                possible_keys.append(serialization_key)
        for key in possible_keys:
            new_value = ResponseSerializer._default_value_picker(value, key, shape)
            if new_value is not None:
                break
        return new_value
