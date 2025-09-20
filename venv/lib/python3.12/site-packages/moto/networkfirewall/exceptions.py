"""Exceptions raised by the networkfirewall service."""

import json

from moto.core.exceptions import JsonRESTError


class ResourceNotFound(JsonRESTError):
    def __init__(self, resource_id: str):
        super().__init__("NotFoundException", "Resource not found.")
        body = {
            "ResourceId": resource_id,
            "Message": "Resource not found.",
        }
        self.description = json.dumps(body)
