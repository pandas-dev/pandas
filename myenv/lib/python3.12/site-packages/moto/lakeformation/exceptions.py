"""Exceptions raised by the lakeformation service."""

from moto.core.exceptions import JsonRESTError


class EntityNotFound(JsonRESTError):
    def __init__(self) -> None:
        super().__init__("EntityNotFoundException", "Entity not found")


class InvalidInput(JsonRESTError):
    def __init__(self, message: str) -> None:
        super().__init__("InvalidInputException", message)


class AlreadyExists(JsonRESTError):
    def __init__(self, message: str) -> None:
        super().__init__("AlreadyExistsException", message)
