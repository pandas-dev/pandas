"""Exceptions raised by the s3tables service."""

from moto.core.exceptions import JsonRESTError


class BadRequestException(JsonRESTError):
    code = 400

    def __init__(self, message: str) -> None:
        super().__init__("BadRequestException", message)


class InvalidContinuationToken(BadRequestException):
    msg = "The continuation token is not valid."

    def __init__(self) -> None:
        super().__init__(self.msg)


class InvalidTableBucketName(BadRequestException):
    msg = "The specified bucket name is not valid."

    def __init__(self) -> None:
        super().__init__(self.msg)


class InvalidTableName(BadRequestException):
    template = "1 validation error detected: Value '%s' at 'name' failed to satisfy constraint: Member must satisfy regular expression pattern: [0-9a-z_]*"

    def __init__(self, name: str) -> None:
        super().__init__(self.template.format(name))


class InvalidNamespaceName(BadRequestException):
    msg = "The specified namespace name is not valid."

    def __init__(self) -> None:
        super().__init__(self.msg)


class InvalidMetadataLocation(BadRequestException):
    msg = "The specified metadata location is not valid."

    def __init__(self) -> None:
        super().__init__(self.msg)


class NothingToRename(BadRequestException):
    msg = "Neither a new namespace name nor a new table name is specified."

    def __init__(self) -> None:
        super().__init__(self.msg)


class NotFoundException(JsonRESTError):
    code = 404

    def __init__(self, message: str) -> None:
        super().__init__("NotFoundException", message)


class NamespaceDoesNotExist(NotFoundException):
    msg = "The specified namespace does not exist."

    def __init__(self) -> None:
        super().__init__(self.msg)


class DestinationNamespaceDoesNotExist(NotFoundException):
    msg = "The specified destination namespace does not exist."

    def __init__(self) -> None:
        super().__init__(self.msg)


class TableDoesNotExist(NotFoundException):
    msg = "The specified table does not exist."

    def __init__(self) -> None:
        super().__init__(self.msg)


class ConflictException(JsonRESTError):
    code = 409

    def __init__(self, message: str) -> None:
        super().__init__("ConflictException", message)


class VersionTokenMismatch(ConflictException):
    msg = "Provided version token does not match the table version token."

    def __init__(self) -> None:
        super().__init__(self.msg)


class TableAlreadyExists(ConflictException):
    msg = "A table with an identical name already exists in the namespace."

    def __init__(self) -> None:
        super().__init__(self.msg)
