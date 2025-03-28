"""Exceptions raised by the Firehose service."""

from moto.core.exceptions import JsonRESTError


class ConcurrentModificationException(JsonRESTError):
    """Existing config has a version ID that does not match given ID."""

    code = 400

    def __init__(self, message: str):
        super().__init__("ConcurrentModificationException", message)


class InvalidArgumentException(JsonRESTError):
    """The specified input parameter has a value that is not valid."""

    code = 400

    def __init__(self, message: str):
        super().__init__("InvalidArgumentException", message)


class LimitExceededException(JsonRESTError):
    """You have already reached the limit for a requested resource."""

    code = 400

    def __init__(self, message: str):
        super().__init__("LimitExceededException", message)


class ResourceInUseException(JsonRESTError):
    """The resource is already in use and not available for this operation."""

    code = 400

    def __init__(self, message: str):
        super().__init__("ResourceInUseException", message)


class ResourceNotFoundException(JsonRESTError):
    """The specified resource could not be found."""

    code = 400

    def __init__(self, message: str):
        super().__init__("ResourceNotFoundException", message)


class ValidationException(JsonRESTError):
    """The tag key or tag value is not valid."""

    code = 400

    def __init__(self, message: str):
        super().__init__("ValidationException", message)
