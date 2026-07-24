from moto.core.exceptions import ServiceException


class InvalidRequestException(ServiceException):
    code = "InvalidRequestException"


class ResourceNotFoundException(InvalidRequestException):
    error_code = "ResourceNotFound"


class InvalidCluster(InvalidRequestException):
    error_code = "NoSuchCluster"

    def __init__(self, cluster_id: str):
        message = f"Cluster id '{cluster_id}' is not valid."
        super().__init__(message)


class ValidationException(ServiceException):
    code = "ValidationException"
