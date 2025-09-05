"""Exceptions raised by the backup service."""

from moto.core.exceptions import JsonRESTError


class BackupClientError(JsonRESTError):
    code = 400


class AlreadyExistsException(BackupClientError):
    def __init__(self, msg: str):
        super().__init__("AlreadyExistsException", f"{msg}")


class ResourceNotFoundException(JsonRESTError):
    def __init__(self, msg: str):
        super().__init__("ResourceNotFoundException", f"{msg}")
