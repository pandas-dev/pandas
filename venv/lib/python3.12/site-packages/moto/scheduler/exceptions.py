"""Exceptions raised by the scheduler service."""

from moto.core.exceptions import JsonRESTError


class ScheduleExists(JsonRESTError):
    def __init__(self, name: str) -> None:
        super().__init__("ConflictException", f"Schedule {name} already exists.")


class ScheduleNotFound(JsonRESTError):
    code = 404

    def __init__(self, name: str) -> None:
        super().__init__(
            "ResourceNotFoundException", f"Schedule {name} does not exist."
        )


class ScheduleGroupNotFound(JsonRESTError):
    code = 404

    def __init__(self, name: str) -> None:
        super().__init__(
            "ResourceNotFoundException", f"Schedule group {name} does not exist."
        )


class ValidationException(JsonRESTError):
    code = 400

    def __init__(self, message: str):
        super().__init__(error_type="ValidationException", message=message)
