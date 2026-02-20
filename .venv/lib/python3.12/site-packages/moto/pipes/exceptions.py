"""Exceptions raised by the pipes service."""

import json
from typing import Any, Optional

from moto.core.exceptions import JsonRESTError


class ValidationException(JsonRESTError):
    code = 400

    # TODO: Add field_list parameter in future iterations as it is not being used yet in current codebase
    def __init__(self, message: str):
        super().__init__("ValidationException", message)
        body: dict[str, Any] = {"message": self.message}
        self.description = json.dumps(body)


class ConflictException(JsonRESTError):
    code = 409

    def __init__(
        self,
        message: str,
        resource_id: Optional[str] = None,
        resource_type: Optional[str] = None,
    ):
        super().__init__("ConflictException", message)
        body: dict[str, Any] = {"message": self.message}
        if resource_id is not None:
            body["resourceId"] = resource_id
        if resource_type is not None:
            body["resourceType"] = resource_type
        self.description = json.dumps(body)


class NotFoundException(JsonRESTError):
    code = 404

    def __init__(self, message: str):
        super().__init__("NotFoundException", message)


class InternalException(JsonRESTError):
    code = 500

    def __init__(self, message: str, retry_after_seconds: Optional[int] = None):
        super().__init__("InternalException", message)
        body: dict[str, Any] = {"message": self.message}
        if retry_after_seconds is not None:
            body["retryAfterSeconds"] = retry_after_seconds
        self.description = json.dumps(body)


class ServiceQuotaExceededException(JsonRESTError):
    code = 402

    def __init__(
        self,
        message: str,
        quota_code: Optional[str] = None,
        resource_id: Optional[str] = None,
        resource_type: Optional[str] = None,
        service_code: Optional[str] = None,
    ):
        super().__init__("ServiceQuotaExceededException", message)
        body: dict[str, Any] = {"message": self.message}
        if quota_code is not None:
            body["quotaCode"] = quota_code
        if resource_id is not None:
            body["resourceId"] = resource_id
        if resource_type is not None:
            body["resourceType"] = resource_type
        if service_code is not None:
            body["serviceCode"] = service_code
        self.description = json.dumps(body)


class ThrottlingException(JsonRESTError):
    code = 429

    def __init__(
        self,
        message: str,
        quota_code: Optional[str] = None,
        retry_after_seconds: Optional[int] = None,
        service_code: Optional[str] = None,
    ):
        super().__init__("ThrottlingException", message)
        body: dict[str, Any] = {"message": self.message}
        if quota_code is not None:
            body["quotaCode"] = quota_code
        if retry_after_seconds is not None:
            body["retryAfterSeconds"] = retry_after_seconds
        if service_code is not None:
            body["serviceCode"] = service_code
        self.description = json.dumps(body)


ResourceNotFoundException = NotFoundException
