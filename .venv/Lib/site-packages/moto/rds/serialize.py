# mypy: disable-error-code="misc, var-annotated"

from __future__ import annotations

from collections import defaultdict
from datetime import datetime
from functools import lru_cache
from typing import Any, Dict, List, Mapping, MutableMapping, Optional, Tuple, Union

import xmltodict
from botocore import xform_name
from botocore.model import (
    ListShape,
    NoShapeFoundError,
    OperationModel,
    Shape,
    StructureShape,
)
from botocore.utils import parse_to_aware_datetime

from .utils import get_service_model

Serialized = MutableMapping[str, Any]


class ErrorShape(StructureShape):
    pass


class ShapeHelpersMixin:
    @staticmethod
    def get_serialized_name(shape: Shape, default_name: str) -> str:
        return shape.serialization.get("name", default_name)


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
        timestamp_format = self.ISO8601
        if value.microsecond > 0:
            timestamp_format = self.ISO8601_MICRO
        return value.strftime(timestamp_format)

    def _convert_timestamp_to_str(
        self, value: Union[int, str, datetime], timestamp_format: str
    ) -> str:
        timestamp_format = timestamp_format.lower()
        converter = getattr(self, "_timestamp_%s" % timestamp_format)
        datetime_obj = parse_to_aware_datetime(value)  # type: ignore
        final_value = converter(datetime_obj)
        return final_value


class Serializer(ShapeHelpersMixin):  # , BaseSerializer):
    DEFAULT_RESPONSE_CODE = 200
    DEFAULT_ERROR_RESPONSE_CODE = 400
    # Clients can change this to a different MutableMapping
    # (i.e. OrderedDict) if they want.  This is used in the
    # compliance test to match the hash ordering used in the
    # tests.
    # NOTE: This is no longer necessary because dicts post 3.6 are ordered
    # https://stackoverflow.com/questions/39980323/are-dictionaries-ordered-in-python-3-6
    MAP_TYPE = dict
    DEFAULT_ENCODING = "utf-8"

    # From the spec, the default timestamp format if not specified is iso8601.
    DEFAULT_TIMESTAMP_FORMAT = TimestampSerializer.TIMESTAMP_FORMAT_ISO8601

    def __init__(
        self,
        operation_model: OperationModel,
        context: Optional[dict[str, Any]] = None,
        value_picker: Any = None,
        pretty_print: bool = False,
    ) -> None:
        self.operation_model = operation_model
        self.context = context or {"request_id": "request-id"}
        self.pretty_print = pretty_print
        if value_picker is None:
            value_picker = self._default_value_picker
        self._value_picker = value_picker
        self._timestamp_serializer = TimestampSerializer(self.DEFAULT_TIMESTAMP_FORMAT)

    def serialize_to_response(
        self,
        result: Any,
    ) -> Mapping[str, Any]:
        raise NotImplementedError("serialize_to_response")

    def _create_default_response(self) -> Serialized:
        # Creates a boilerplate default request dict that subclasses
        # can use as a starting point.
        serialized = {
            "status_code": self.DEFAULT_RESPONSE_CODE,
            "headers": {},
            # An empty body is represented as an empty string.
            "body": "",
        }
        return serialized

    # Some extra utility methods subclasses can use.

    @staticmethod
    def _is_error_result(result: object) -> bool:
        return isinstance(result, Exception)

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


class ResponseSerializer(Serializer):
    DEFAULT_TIMESTAMP_FORMAT = TimestampSerializer.TIMESTAMP_FORMAT_ISO8601

    CONTENT_TYPE = "text"

    def _encode_body(self, body: Any) -> str:
        raise NotImplementedError("_encode_body")

    @staticmethod
    def _get_error_shape_name(error: Exception) -> str:
        shape_name = getattr(error, "code", error.__class__.__name__)
        return shape_name

    def _get_error_shape(
        self, error: Exception, operation_model: OperationModel
    ) -> ErrorShape:
        shape_name = self._get_error_shape_name(error)
        try:
            # TODO: there is also an errors array in the operation model,
            # but I think it only includes the possible errors for that
            # operation.  Maybe we try that first, then try all shapes?
            shape = operation_model.service_model.shape_for(shape_name)
            # We convert to ErrorShape to keep mypy happy...
            shape = ErrorShape(
                shape_name,
                shape._shape_model,  # type: ignore
                shape._shape_resolver,  # type: ignore
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

    def _serialize_error(
        self,
        serialized: Serialized,
        error: Exception,
        operation_model: OperationModel,
        request_ctx: Mapping[str, Any],
    ) -> Serialized:
        shape = self._get_error_shape(error, operation_model)
        status_code = shape.metadata.get("error", {}).get(
            "httpStatusCode", self.DEFAULT_ERROR_RESPONSE_CODE
        )
        serialized["status_code"] = status_code
        message = getattr(error, "message", None) or str(error)
        error_wrapper, error_body = self._get_error_wrapper(
            operation_model, request_ctx
        )
        self._inject_error_metadata(error_body, error, shape, operation_model)
        if message:
            error_body["Message"] = message
        if shape is not None:
            self._serialize(error_body, error, shape, "")
        serialized["body"] = error_wrapper
        self._inject_error_headers(serialized["headers"], shape, operation_model)
        return serialized

    def _get_error_wrapper(
        self, operation_model: OperationModel, request_ctx: Mapping[str, Any]
    ) -> Tuple[Serialized, Serialized]:
        raise NotImplementedError("_get_error_wrapper")

    def _inject_error_metadata(
        self,
        serialized: Serialized,
        error: Exception,
        shape: ErrorShape,
        operation_model: OperationModel,
    ) -> None:
        raise NotImplementedError("_inject_error_metadata")

    def _inject_error_headers(
        self, serialized: Serialized, shape: ErrorShape, operation_model: OperationModel
    ) -> None:
        pass

    def _inject_response_metadata(
        self,
        serialized: Serialized,
        operation_model: OperationModel,
        request_ctx: Mapping[str, Any],
    ) -> None:
        raise NotImplementedError("_inject_response_metadata")

    def serialize_to_response(
        self,
        result: Any,
    ) -> Serialized:
        serialized = self._create_default_response()
        if self._is_error_result(result):
            serialized = self._serialize_error(
                serialized, result, self.operation_model, self.context
            )
        else:
            response_wrapper, result_wrapper = self._get_response_wrapper(
                self.operation_model, self.context
            )
            output_shape = self.operation_model.output_shape
            if output_shape is not None:
                self._serialize(result_wrapper, result, output_shape, "")

            root_key = list(response_wrapper.keys())[0]
            self._inject_response_metadata(
                response_wrapper[root_key], self.operation_model, self.context
            )

            serialized["body"] = response_wrapper

        serialized["body"] = self._encode_body(serialized["body"])

        serialized["headers"]["Content-Type"] = self.CONTENT_TYPE
        return serialized

    def _serialize(
        self, serialized: Serialized, value: Any, shape: Shape, key: str
    ) -> None:
        method = getattr(
            self, "_serialize_type_%s" % shape.type_name, self._default_serialize
        )
        method(serialized, value, shape, key)

    def _serialize_type_structure(
        self, serialized: Serialized, value: Any, shape: StructureShape, key: str
    ) -> None:
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

    @staticmethod
    def _default_serialize(
        serialized: Serialized, value: Any, _: Shape, key: str
    ) -> None:
        serialized[key] = value

    def _serialize_type_timestamp(
        self, serialized: Serialized, value: Any, shape: Shape, key: str
    ) -> None:
        value_wrapper = {}
        value_key = "timestamp"
        self._timestamp_serializer.serialize(value_wrapper, value, shape, value_key)
        self._default_serialize(serialized, value_wrapper[value_key], shape, key)

    def _get_response_wrapper(
        self, operation_model: OperationModel, request_ctx: Mapping[str, Any]
    ) -> Tuple[Serialized, Serialized]:
        raise NotImplementedError("shouldn't get here...")


class BaseXMLSerializer(ResponseSerializer):
    @staticmethod
    def _serialize_namespace_attribute(
        serialized: Serialized, operation_model: OperationModel
    ) -> None:
        if "xmlNamespace" in operation_model.metadata:
            namespace = operation_model.metadata["xmlNamespace"]
            serialized["@xmlns"] = namespace

    def _get_error_wrapper(
        self, operation_model: OperationModel, request_ctx: Mapping[str, Any]
    ) -> Tuple[Serialized, Serialized]:
        serialized_error = self.MAP_TYPE()
        error_wrapper = {
            "ErrorResponse": {
                "Error": serialized_error,
                "RequestId": request_ctx.get("request_id"),
            }
        }
        self._serialize_namespace_attribute(
            error_wrapper["ErrorResponse"], operation_model
        )
        return error_wrapper, serialized_error

    def _get_response_wrapper(
        self, operation_model: OperationModel, request_ctx: Mapping[str, Any]
    ) -> Tuple[Serialized, Serialized]:
        serialized_result = self.MAP_TYPE()
        root_key = f"{operation_model.name}Response"
        response_wrapper = {root_key: {}}
        output_shape = operation_model.output_shape
        result_key = None
        if output_shape is not None:
            result_key = output_shape.serialization.get("resultWrapper")
        if result_key is not None:
            response_wrapper[root_key][result_key] = serialized_result
        else:
            response_wrapper[root_key] = serialized_result
        self._serialize_namespace_attribute(response_wrapper[root_key], operation_model)
        return response_wrapper, serialized_result

    def _inject_error_metadata(
        self,
        serialized: Serialized,
        error: Exception,
        shape: ErrorShape,
        operation_model: OperationModel,
    ) -> None:
        sender_fault = shape.metadata.get("error", {}).get("senderFault", True)
        serialized["Type"] = "Sender" if sender_fault else "Receiver"
        serialized["Code"] = shape.error_code
        message = getattr(error, "message", None) or str(error)
        if message:
            serialized["Message"] = message
        if shape is not None:
            self._serialize(serialized, error, shape, "")

    def _encode_body(self, body: Serialized) -> str:
        body_encoded = xmltodict.unparse(
            body,
            full_document=False,
            short_empty_elements=True,
            pretty=self.pretty_print,
        )
        return body_encoded

    def _serialize_type_list(
        self, serialized: Serialized, value: Any, shape: ListShape, key: str
    ) -> None:
        list_obj = []
        items_name = self.get_serialized_name(shape.member, "member")
        serialized[key] = {items_name: list_obj}
        for list_item in value:
            wrapper = {}
            self._serialize(wrapper, list_item, shape.member, "__current__")
            list_obj.append(wrapper["__current__"])
        if not list_obj:
            serialized[key] = ""


class QuerySerializer(BaseXMLSerializer):
    CONTENT_TYPE = "text/xml"

    def _inject_response_metadata(
        self,
        serialized: MutableMapping[str, Any],
        operation_model: OperationModel,
        request_ctx: Mapping[str, Any],
    ) -> None:
        serialized["ResponseMetadata"] = {"RequestId": request_ctx.get("request_id")}


SERIALIZERS = {
    "query": QuerySerializer,
}


class XFormedAttributeAccessMixin:
    """Mixin allowing access to "xformed" attributes:

    obj.DBInstanceIdentifier will retrieve the value of obj.db_instance_identifier

    """

    BOTOCORE_MODEL: Optional[str] = None

    _model_attribute_aliases: Dict[str, List[str]] = {}
    _xform_cache: Dict[str, str] = {}

    def __getattr__(self, name: str) -> Any:
        if name in self.model_attributes:
            return self.get_modeled_attribute(name)
        raise AttributeError(f"Attribute '{name}' not found!")

    def get_modeled_attribute(self, attr_name: str) -> Any:
        for attr_alias in self.model_attribute_aliases[attr_name]:
            try:
                return super().__getattribute__(attr_alias)
            except AttributeError:
                pass
        else:
            raise AttributeError

    @property
    def model_attributes(self) -> List[str]:
        return list(self.model_attribute_aliases.keys())

    @property
    def model_attribute_aliases(self) -> Dict[str, List[str]]:
        if not self._model_attribute_aliases:
            self._model_attribute_aliases = self.get_model_attributes_info()
        return self._model_attribute_aliases

    @classmethod
    @lru_cache()
    def get_model_attributes_info(cls) -> Dict[str, List[str]]:
        service_name = cls.__module__.split(".")[1]
        model_name = cls.BOTOCORE_MODEL or cls.__name__
        service_model = get_service_model(service_name)
        model_shape = service_model.shape_for(model_name)
        valid_attributes: Dict[str, List[str]] = defaultdict(list)
        for member_name, member_shape in model_shape.members.items():  # type: ignore[attr-defined]
            aliases = valid_attributes[member_name]
            if member_shape.type_name == "list":
                if member_name != member_shape.name:
                    xformed_name = cls._xform_name(member_shape.name)
                    aliases.append(xformed_name)
            xformed_member_name = cls._xform_name(member_name)
            aliases.append(xformed_member_name)
            if member_name.startswith(model_name):
                short_name = member_name[len(model_name) :]
                xformed_short_name = cls._xform_name(short_name)
                aliases.append(xformed_short_name)
        return valid_attributes

    @classmethod
    def _xform_name(cls, name: str) -> str:
        return xform_name(name, _xform_cache=cls._xform_cache)
