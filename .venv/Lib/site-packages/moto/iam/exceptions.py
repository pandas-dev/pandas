from typing import Any

from moto.core.exceptions import RESTError

XMLNS_IAM = "https://iam.amazonaws.com/doc/2010-05-08/"


class IAMNotFoundException(RESTError):
    code = 404

    def __init__(self, message: str):
        super().__init__(
            "NoSuchEntity", message, xmlns=XMLNS_IAM, template="wrapped_single_error"
        )


class IAMConflictException(RESTError):
    code = 409

    def __init__(self, code: str = "Conflict", message: str = ""):
        super().__init__(code, message)


class IAMReportNotPresentException(RESTError):
    code = 410

    def __init__(self, message: str):
        super().__init__("ReportNotPresent", message)


class IAMLimitExceededException(RESTError):
    code = 400

    def __init__(self, message: str):
        super().__init__("LimitExceeded", message)


class MalformedCertificate(RESTError):
    code = 400

    def __init__(self, cert: str):
        super().__init__("MalformedCertificate", f"Certificate {cert} is malformed")


class MalformedPolicyDocument(RESTError):
    code = 400

    def __init__(self, message: str = ""):
        super().__init__(
            "MalformedPolicyDocument",
            message,
            xmlns=XMLNS_IAM,
            template="wrapped_single_error",
        )


class DuplicateTags(RESTError):
    code = 400

    def __init__(self) -> None:
        super().__init__(
            "InvalidInput",
            "Duplicate tag keys found. Please note that Tag keys are case insensitive.",
        )


class TagKeyTooBig(RESTError):
    code = 400

    def __init__(self, tag: str, param: str = "tags.X.member.key"):
        super().__init__(
            "ValidationError",
            f"1 validation error detected: Value '{tag}' at '{param}' failed to satisfy "
            "constraint: Member must have length less than or equal to 128.",
        )


class TagValueTooBig(RESTError):
    code = 400

    def __init__(self, tag: str):
        super().__init__(
            "ValidationError",
            f"1 validation error detected: Value '{tag}' at 'tags.X.member.value' failed to satisfy "
            "constraint: Member must have length less than or equal to 256.",
        )


class InvalidTagCharacters(RESTError):
    code = 400

    def __init__(self, tag: str, param: str = "tags.X.member.key"):
        message = f"1 validation error detected: Value '{tag}' at '{param}' failed to satisfy constraint: Member must satisfy regular expression pattern: [\\p{{L}}\\p{{Z}}\\p{{N}}_.:/=+\\-@]+"

        super().__init__("ValidationError", message)


class TooManyTags(RESTError):
    code = 400

    def __init__(self, tags: Any, param: str = "tags"):
        super().__init__(
            "ValidationError",
            f"1 validation error detected: Value '{tags}' at '{param}' failed to satisfy "
            "constraint: Member must have length less than or equal to 50.",
        )


class EntityAlreadyExists(RESTError):
    code = 409

    def __init__(self, message: str):
        super().__init__("EntityAlreadyExists", message)


class ValidationError(RESTError):
    code = 400

    def __init__(self, message: str):
        super().__init__("ValidationError", message)


class InvalidInput(RESTError):
    code = 400

    def __init__(self, message: str):
        super().__init__("InvalidInput", message)


class NoSuchEntity(RESTError):
    code = 404

    def __init__(self, message: str):
        super().__init__(
            "NoSuchEntity", message, xmlns=XMLNS_IAM, template="wrapped_single_error"
        )
