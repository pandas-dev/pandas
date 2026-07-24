"""Exceptions raised by the bedrock service."""

from moto.core.exceptions import JsonRESTError

# Bedrock.Client.exceptions.ResourceNotFoundException


class BedrockClientError(JsonRESTError):
    code = 400


class ResourceNotFoundException(BedrockClientError):
    def __init__(self, msg: str):
        super().__init__("ResourceNotFoundException", f"{msg}")


class ResourceInUseException(BedrockClientError):
    def __init__(self, msg: str):
        super().__init__("ResourceInUseException", f"{msg}")


class ValidationException(BedrockClientError):
    def __init__(self, msg: str):
        super().__init__(
            "ValidationException",
            "Input validation failed. Check your request parameters and retry the request.",
            f"{msg}",
        )


class TooManyTagsException(BedrockClientError):
    def __init__(self, msg: str):
        super().__init__("TooManyTagsException", f"{msg}")
