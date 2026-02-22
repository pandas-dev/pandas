from typing import Optional

from moto.core.exceptions import ServiceException


class CloudFormationError(ServiceException):
    pass


class UnformattedGetAttTemplateException(Exception):
    description = (
        "Template error: resource {0} does not support attribute type {1} in Fn::GetAtt"
    )
    status_code = 400


class AlreadyExistsException(CloudFormationError):
    code = "AlreadyExistsException"


class ValidationError(CloudFormationError):
    code = "ValidationError"

    def __init__(self, name_or_id: Optional[str] = None, message: Optional[str] = None):
        # FIXME: The "stack does not exist" message should be provided by the caller, not here.
        if message is None:
            message = f"Stack with id {name_or_id} does not exist"
        super().__init__(message)


class MissingParameterError(CloudFormationError):
    code = "ValidationError"

    def __init__(self, parameter_name: str):
        message = f"Missing parameter {parameter_name}"
        super().__init__(message)


class ExportNotFound(CloudFormationError):
    code = "ValidationError"

    def __init__(self, export_name: str):
        message = f"No export named {export_name} found."
        super().__init__(message)


class StackSetNotEmpty(CloudFormationError):
    code = "StackSetNotEmptyException"
    message = "StackSet is not empty"


class StackSetNotFoundException(CloudFormationError):
    code = "StackSetNotFoundException"

    def __init__(self, name: str):
        message = f"StackSet {name} not found"
        super().__init__(message)


class UnsupportedAttribute(CloudFormationError):
    code = "ValidationError"

    def __init__(self, resource: str, attr: str):
        super().__init__(
            f"Template error: resource {resource} does not support attribute type {attr} in Fn::GetAtt"
        )
