"""Exceptions raised by the networkmanager service."""

import json

from moto.core.exceptions import JsonRESTError


class ValidationError(JsonRESTError):
    def __init__(self, message: str):
        super().__init__("ValidationException", message)


class ResourceNotFound(JsonRESTError):
    def __init__(self, resource_id: str):
        super().__init__("NotFoundException", "Resource not found.")
        body = {
            "ResourceId": resource_id,
            "Message": "Resource not found.",
        }
        self.description = json.dumps(body)
