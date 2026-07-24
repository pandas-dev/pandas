from typing import Any

from moto.core.exceptions import ServiceException


class IAMException(ServiceException):
    pass


class NotFoundException(IAMException):
    code = "NoSuchEntity"


class DeleteConflictException(IAMException):
    code = "DeleteConflict"


class ReportNotPresentException(IAMException):
    code = "ReportNotPresent"


class LimitExceededException(IAMException):
    code = "LimitExceeded"


class MalformedCertificate(IAMException):
    code = "MalformedCertificate"

    def __init__(self, cert: str):
        super().__init__(f"Certificate {cert} is malformed")


class MalformedPolicyDocument(IAMException):
    code = "MalformedPolicyDocument"


class DuplicateTags(IAMException):
    code = "InvalidInput"
    message = (
        "Duplicate tag keys found. Please note that Tag keys are case insensitive."
    )


class ValidationError(IAMException):
    code = "ValidationError"


class TagKeyTooBig(ValidationError):
    def __init__(self, tag: str, param: str = "tags.X.member.key"):
        super().__init__(
            f"1 validation error detected: Value '{tag}' at '{param}' failed to satisfy "
            "constraint: Member must have length less than or equal to 128.",
        )


class TagValueTooBig(ValidationError):
    def __init__(self, tag: str):
        super().__init__(
            f"1 validation error detected: Value '{tag}' at 'tags.X.member.value' failed to satisfy "
            "constraint: Member must have length less than or equal to 256.",
        )


class InvalidTagCharacters(ValidationError):
    def __init__(self, tag: str, param: str = "tags.X.member.key"):
        message = f"1 validation error detected: Value '{tag}' at '{param}' failed to satisfy constraint: Member must satisfy regular expression pattern: [\\p{{L}}\\p{{Z}}\\p{{N}}_.:/=+\\-@]+"
        super().__init__(message)


class TooManyTags(ValidationError):
    def __init__(self, tags: Any, param: str = "tags"):
        super().__init__(
            f"1 validation error detected: Value '{tags}' at '{param}' failed to satisfy "
            "constraint: Member must have length less than or equal to 50.",
        )


class EntityAlreadyExists(IAMException):
    code = "EntityAlreadyExists"


class InvalidInput(IAMException):
    code = "InvalidInput"


class NoSuchEntity(IAMException):
    code = "NoSuchEntity"
