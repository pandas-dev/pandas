"""Exceptions raised by the opensearchserverless service."""

from moto.core.exceptions import JsonRESTError


class BackupClientError(JsonRESTError):
    code = 400


class ConflictException(BackupClientError):
    def __init__(self, msg: str):
        super().__init__("ConflictException", f"{msg}")


class ValidationException(BackupClientError):
    def __init__(self, msg: str):
        super().__init__("ValidationException", f"{msg}")


class ResourceNotFoundException(BackupClientError):
    def __init__(self, msg: str):
        super().__init__("ResourceNotFoundException", f"{msg}")
