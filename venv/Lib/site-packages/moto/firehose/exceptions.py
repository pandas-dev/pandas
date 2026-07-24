"""Exceptions raised by the Firehose service."""

from moto.core.exceptions import ServiceException


class FirehoseException(ServiceException):
    pass


class ConcurrentModificationException(FirehoseException):
    """Existing config has a version ID that does not match given ID."""

    def __init__(self, message: str):
        super().__init__("ConcurrentModificationException", message)


class InvalidArgumentException(FirehoseException):
    """The specified input parameter has a value that is not valid."""

    def __init__(self, message: str):
        super().__init__("InvalidArgumentException", message)


class LimitExceededException(FirehoseException):
    """You have already reached the limit for a requested resource."""

    def __init__(self, message: str):
        super().__init__("LimitExceededException", message)


class ResourceInUseException(FirehoseException):
    """The resource is already in use and not available for this operation."""

    def __init__(self, message: str):
        super().__init__("ResourceInUseException", message)


class ResourceNotFoundException(FirehoseException):
    """The specified resource could not be found."""

    def __init__(self, message: str):
        super().__init__("ResourceNotFoundException", message)


class ValidationException(FirehoseException):
    """The tag key or tag value is not valid."""

    def __init__(self, message: str):
        super().__init__("ValidationException", message)
