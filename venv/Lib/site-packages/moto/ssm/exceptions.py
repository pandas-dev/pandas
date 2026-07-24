from moto.core.exceptions import ServiceException


class SSMException(ServiceException):
    pass


class DoesNotExistException(SSMException):
    def __init__(self, window_id: str):
        super().__init__(
            "DoesNotExistException", f"Maintenance window {window_id} does not exist"
        )


class InvalidFilterKey(SSMException):
    def __init__(self, message: str):
        super().__init__("InvalidFilterKey", message)


class InvalidFilterOption(SSMException):
    def __init__(self, message: str):
        super().__init__("InvalidFilterOption", message)


class InvalidFilterValue(SSMException):
    def __init__(self, message: str):
        super().__init__("InvalidFilterValue", message)


class InvalidResourceId(SSMException):
    def __init__(self) -> None:
        super().__init__("InvalidResourceId", "Invalid Resource Id")


class InvalidResourceType(SSMException):
    def __init__(self) -> None:
        super().__init__("InvalidResourceType", "Invalid Resource Type")


class ParameterNotFound(SSMException):
    def __init__(self, message: str):
        super().__init__("ParameterNotFound", message)


class ParameterVersionNotFound(SSMException):
    def __init__(self, message: str):
        super().__init__("ParameterVersionNotFound", message)


class ParameterVersionLabelLimitExceeded(SSMException):
    def __init__(self, message: str):
        super().__init__("ParameterVersionLabelLimitExceeded", message)


class ValidationException(SSMException):
    def __init__(self, message: str):
        super().__init__("ValidationException", message)


class DocumentAlreadyExists(SSMException):
    def __init__(self, message: str):
        super().__init__("DocumentAlreadyExists", message)


class DocumentPermissionLimit(SSMException):
    def __init__(self, message: str):
        super().__init__("DocumentPermissionLimit", message)


class InvalidPermissionType(SSMException):
    def __init__(self, message: str):
        super().__init__("InvalidPermissionType", message)


class InvalidDocument(SSMException):
    def __init__(self, message: str):
        super().__init__("InvalidDocument", message)


class InvalidDocumentOperation(SSMException):
    def __init__(self, message: str):
        super().__init__("InvalidDocumentOperation", message)


class AccessDeniedException(SSMException):
    def __init__(self, message: str):
        super().__init__("AccessDeniedException", message)


class InvalidDocumentContent(SSMException):
    def __init__(self, message: str):
        super().__init__("InvalidDocumentContent", message)


class InvalidDocumentVersion(SSMException):
    def __init__(self, message: str):
        super().__init__("InvalidDocumentVersion", message)


class DuplicateDocumentVersionName(SSMException):
    def __init__(self, message: str):
        super().__init__("DuplicateDocumentVersionName", message)


class DuplicateDocumentContent(SSMException):
    def __init__(self, message: str):
        super().__init__("DuplicateDocumentContent", message)


class ParameterMaxVersionLimitExceeded(SSMException):
    def __init__(self, message: str):
        super().__init__("ParameterMaxVersionLimitExceeded", message)


class ParameterAlreadyExists(SSMException):
    def __init__(self) -> None:
        super().__init__(
            "ParameterAlreadyExists",
            "The parameter already exists. To overwrite this value, set the overwrite option in the request to true.",
        )


class AlreadyExistsException(SSMException):
    def __init__(self, operating_system: str) -> None:
        super().__init__(
            "AlreadyExistsException",
            f"Patch Group baseline already has a baseline registered for OperatingSystem {operating_system}.",
        )


class BaselineDoesNotExistException(SSMException):
    def __init__(self) -> None:
        super().__init__(
            "DoesNotExistException", "Patch Baseline to be retrieved does not exist."
        )
