"""Exceptions raised by the bedrockagent service."""

from moto.core.exceptions import JsonRESTError


class AgentsforBedrockClientError(JsonRESTError):
    code = 400


class ResourceNotFoundException(AgentsforBedrockClientError):
    def __init__(self, msg: str):
        super().__init__("ResourceNotFoundException", f"{msg}")


class ConflictException(AgentsforBedrockClientError):
    def __init__(self, msg: str):
        super().__init__("ConflictException", f"{msg}")


class ValidationException(AgentsforBedrockClientError):
    def __init__(self, msg: str):
        super().__init__(
            "ValidationException",
            "Input validation failed. Check your request parameters and retry the request.",
            f"{msg}",
        )
