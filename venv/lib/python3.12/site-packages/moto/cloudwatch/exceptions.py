from moto.core.exceptions import ServiceException


class CloudWatchException(ServiceException):
    pass


class InvalidFormat(CloudWatchException):
    code = "InvalidFormat"


class InvalidParameterValue(CloudWatchException):
    code = "InvalidParameterValue"


class InvalidParameterCombination(CloudWatchException):
    code = "InvalidParameterCombination"


class ResourceNotFound(CloudWatchException):
    code = "ResourceNotFound"


class ResourceNotFoundException(CloudWatchException):
    code = "ResourceNotFoundException"
    message = "Unknown"


class ValidationError(CloudWatchException):
    code = "ValidationError"


class DashboardInvalidInputError(CloudWatchException):
    code = "InvalidParameterInput"
