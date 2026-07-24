"""Exceptions raised by the clouddirectory service."""

import json

from moto.core.exceptions import JsonRESTError


class ValidationError(JsonRESTError):
    def __init__(self, message: str):
        super().__init__("ValidationException", message)


class InvalidArnException(JsonRESTError):
    def __init__(self, resource_id: str):
        super().__init__("InvalidArnException", "Invalid Arn")
        body = {
            "ResourceId": resource_id,
            "Message": "Invalid Arn",
        }
        self.description = json.dumps(body)


class ResourceNotFoundException(JsonRESTError):
    def __init__(self, resource_id: str):
        super().__init__("ResourceNotFoundException", "Resource not found")
        body = {
            "ResourceId": resource_id,
            "Message": "Resource not found",
        }
        self.description = json.dumps(body)


class SchemaAlreadyPublishedException(JsonRESTError):
    def __init__(self, schema_arn: str):
        super().__init__("SchemaAlreadyPublishedException", "Schema already published")
        body = {
            "SchemaArn": schema_arn,
            "Message": "Schema already published",
        }
        self.description = json.dumps(body)


class SchemaAlreadyExistsException(JsonRESTError):
    def __init__(self, schema_arn: str):
        super().__init__("SchemaAlreadyExistsException", "Schema already exists")
        body = {
            "SchemaArn": schema_arn,
            "Message": "Schema already exists",
        }
        self.description = json.dumps(body)
