"""Exceptions raised by the pinpoint service."""

from moto.core.exceptions import JsonRESTError


class PinpointExceptions(JsonRESTError):
    pass


class ApplicationNotFound(PinpointExceptions):
    code = 404

    def __init__(self) -> None:
        super().__init__("NotFoundException", "Application not found")


class EventStreamNotFound(PinpointExceptions):
    code = 404

    def __init__(self) -> None:
        super().__init__("NotFoundException", "Resource not found")
