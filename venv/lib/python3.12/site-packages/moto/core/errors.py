"""
This module provides a Shape subclass to encapsulate service model error definitions.

It also provides a mechanism for mapping ServiceException exception classes to the
corresponding error shape defined in the relevant service model.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any
from warnings import warn

from moto.core.model import StructureShape

if TYPE_CHECKING:
    from moto.core.model import ServiceModel

# These are common error codes that are *not* included in the service definitions.
# For example:
# https://docs.aws.amazon.com/emr/latest/APIReference/CommonErrors.html
# https://docs.aws.amazon.com/AmazonRDS/latest/APIReference/CommonErrors.html
# https://docs.aws.amazon.com/AWSSimpleQueueService/latest/APIReference/CommonErrors.html
# TODO: Augment the service definitions with shape models for these errors.
COMMON_ERROR_CODES = [
    "InvalidParameterCombination",
    "InvalidParameterException",
    "InvalidParameterValue",
    "ValidationError",
    "ValidationException",
]


class ErrorShape(StructureShape):
    _shape_model: dict[str, Any]

    @property
    def error_code(self) -> str:
        code = str(super().error_code)
        return code

    @property
    def query_compatible_error_message(self) -> str:
        error_info = self.metadata.get("error", {})
        error_message = error_info.get("messageForQueryCompatibility", "")
        return error_message

    @property
    def is_sender_fault(self) -> bool:
        internal_fault = self._shape_model.get("fault", False)
        error_info = self.metadata.get("error", {})
        sender_fault = error_info.get("senderFault", False)
        return sender_fault or not internal_fault

    @property
    def namespace(self) -> str | None:
        return self.metadata.get("error", {}).get("namespace")


class ErrorLookup:
    def __init__(self, code_to_error: dict[str, ErrorShape]) -> None:
        self._code_to_error = code_to_error

    def from_code(self, code: str) -> ErrorShape | None:
        return self._code_to_error.get(code)


class ErrorLookupFactory:
    def __init__(self) -> None:
        self._error_lut_cache: dict[str, ErrorLookup] = {}

    def for_service(self, service_model: ServiceModel) -> ErrorLookup:
        service_id = service_model.metadata.get("serviceId")
        if service_id not in self._error_lut_cache:
            error_lut = self._create_error_lut(service_model)
            if service_id is None:
                return error_lut
            self._error_lut_cache[service_id] = error_lut
        return self._error_lut_cache[service_id]

    @staticmethod
    def _create_error_lut(service_model: ServiceModel) -> ErrorLookup:
        """We map an error's code, name, and any alias codes to the same ErrorShape."""
        code_to_shape = {}
        for shape in service_model.error_shapes:
            error_shape = ErrorShape(
                shape.name,
                shape._shape_model,  # type: ignore[attr-defined]
                shape._shape_resolver,  # type: ignore[attr-defined]
            )
            code_to_shape[error_shape.name] = error_shape
            code_to_shape[error_shape.error_code] = error_shape
            for error_code in error_shape.error_code_aliases:
                code_to_shape[error_code] = error_shape
        return ErrorLookup(code_to_shape)


def get_exception_service_model(exception: Exception) -> ServiceModel | None:
    from moto.core.utils import get_service_model

    exception_module = exception.__module__
    if not exception_module.startswith("moto"):
        return None
    service = exception_module.split(".")[1]
    service_model = get_service_model(service)
    return service_model


def get_error_model(
    exception: Exception, default_service_model: ServiceModel
) -> ErrorShape:
    possible_service_models = [
        default_service_model,
        get_exception_service_model(exception),
    ]
    services_checked = []
    code = getattr(exception, "code", exception.__class__.__name__)
    error = None
    for service_model in [sm for sm in possible_service_models if sm is not None]:
        if (service_id := service_model.metadata.get("serviceId")) is not None:
            services_checked.append(service_id)
        error_map = ErrorLookupFactory().for_service(service_model)
        error = error_map.from_code(code)
        if error is not None:
            break
    if error is None:
        if services_checked and code not in COMMON_ERROR_CODES:
            warning = f"Exception({exception.__class__.__name__}) with code {code} does not match an eror shape in service models(s): {services_checked}"  # pragma: no cover
            warn(warning, stacklevel=2)  # pragma: no cover
        error = ErrorShape(
            shape_name=exception.__class__.__name__,
            shape_model={
                "exception": True,
                "type": "structure",
                "members": {},
                "error": {
                    "code": code,
                },
            },
        )
    return error
