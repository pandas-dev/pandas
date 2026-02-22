"""Exceptions raised by the personalize service."""

from moto.core.exceptions import JsonRESTError


class PersonalizeException(JsonRESTError):
    pass


class ResourceNotFoundException(PersonalizeException):
    def __init__(self, arn: str):
        super().__init__(
            "ResourceNotFoundException", f"Resource Arn {arn} does not exist."
        )
