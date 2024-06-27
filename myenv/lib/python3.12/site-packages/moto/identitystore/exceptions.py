"""Exceptions raised by the identitystore service."""

import json
from typing import Any

from moto.core.exceptions import AWSError

request_id = "178936da-50ad-4d58-8871-22d9979e8658example"


class IdentityStoreError(AWSError):
    def __init__(self, **kwargs: Any):
        super(AWSError, self).__init__(error_type=self.TYPE, message=kwargs["message"])  # type: ignore
        self.description: str = json.dumps(
            {
                "__type": self.error_type,
                "RequestId": request_id,
                "Message": self.message,
                "ResourceType": kwargs.get("resource_type"),
                "Reason": kwargs.get("reason"),
            }
        )


class ResourceNotFoundException(IdentityStoreError):
    TYPE = "ResourceNotFoundException"
    code = 400


class ValidationException(IdentityStoreError):
    TYPE = "ValidationException"
    code = 400


class ConflictException(IdentityStoreError):
    TYPE = "ConflictException"
    code = 400
