"""Exceptions raised by the scheduler service."""

from moto.core.exceptions import JsonRESTError


class ScheduleExists(JsonRESTError):
    def __init__(self, name: str) -> None:
        super().__init__("ConflictException", f"Schedule {name} already exists.")


class ScheduleNotFound(JsonRESTError):
    def __init__(self, name: str) -> None:
        super().__init__(
            "ResourceNotFoundException", f"Schedule {name} does not exist."
        )


class ScheduleGroupNotFound(JsonRESTError):
    def __init__(self) -> None:
        super().__init__("ResourceNotFoundException", "ScheduleGroup not found")
