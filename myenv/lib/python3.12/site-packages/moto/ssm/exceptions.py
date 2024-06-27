from moto.core.exceptions import JsonRESTError


class InvalidFilterKey(JsonRESTError):
    code = 400

    def __init__(self, message: str):
        super().__init__("InvalidFilterKey", message)


class InvalidFilterOption(JsonRESTError):
    code = 400

    def __init__(self, message: str):
        super().__init__("InvalidFilterOption", message)


class InvalidFilterValue(JsonRESTError):
    code = 400

    def __init__(self, message: str):
        super().__init__("InvalidFilterValue", message)


class InvalidResourceId(JsonRESTError):
    code = 400

    def __init__(self) -> None:
        super().__init__("InvalidResourceId", "Invalid Resource Id")


class InvalidResourceType(JsonRESTError):
    code = 400

    def __init__(self) -> None:
        super().__init__("InvalidResourceType", "Invalid Resource Type")


class ParameterNotFound(JsonRESTError):
    code = 400

    def __init__(self, message: str):
        super().__init__("ParameterNotFound", message)


class ParameterVersionNotFound(JsonRESTError):
    code = 400

    def __init__(self, message: str):
        super().__init__("ParameterVersionNotFound", message)


class ParameterVersionLabelLimitExceeded(JsonRESTError):
    code = 400

    def __init__(self, message: str):
        super().__init__("ParameterVersionLabelLimitExceeded", message)


class ValidationException(JsonRESTError):
    code = 400

    def __init__(self, message: str):
        super().__init__("ValidationException", message)


class DocumentAlreadyExists(JsonRESTError):
    code = 400

    def __init__(self, message: str):
        super().__init__("DocumentAlreadyExists", message)


class DocumentPermissionLimit(JsonRESTError):
    code = 400

    def __init__(self, message: str):
        super().__init__("DocumentPermissionLimit", message)


class InvalidPermissionType(JsonRESTError):
    code = 400

    def __init__(self, message: str):
        super().__init__("InvalidPermissionType", message)


class InvalidDocument(JsonRESTError):
    code = 400

    def __init__(self, message: str):
        super().__init__("InvalidDocument", message)


class InvalidDocumentOperation(JsonRESTError):
    code = 400

    def __init__(self, message: str):
        super().__init__("InvalidDocumentOperation", message)


class AccessDeniedException(JsonRESTError):
    code = 400

    def __init__(self, message: str):
        super().__init__("AccessDeniedException", message)


class InvalidDocumentContent(JsonRESTError):
    code = 400

    def __init__(self, message: str):
        super().__init__("InvalidDocumentContent", message)


class InvalidDocumentVersion(JsonRESTError):
    code = 400

    def __init__(self, message: str):
        super().__init__("InvalidDocumentVersion", message)


class DuplicateDocumentVersionName(JsonRESTError):
    code = 400

    def __init__(self, message: str):
        super().__init__("DuplicateDocumentVersionName", message)


class DuplicateDocumentContent(JsonRESTError):
    code = 400

    def __init__(self, message: str):
        super().__init__("DuplicateDocumentContent", message)


class ParameterMaxVersionLimitExceeded(JsonRESTError):
    code = 400

    def __init__(self, message: str):
        super().__init__("ParameterMaxVersionLimitExceeded", message)


class ParameterAlreadyExists(JsonRESTError):
    code = 400

    def __init__(self) -> None:
        super().__init__(
            "ParameterAlreadyExists",
            "The parameter already exists. To overwrite this value, set the overwrite option in the request to true.",
        )
