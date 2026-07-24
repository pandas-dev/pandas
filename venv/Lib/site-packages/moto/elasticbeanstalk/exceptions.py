from moto.core.exceptions import ServiceException


class ElasticBeanstalkException(ServiceException):
    pass


class InvalidParameterValueError(ServiceException):
    code = "InvalidParameterValue"


class ResourceNotFoundException(ServiceException):
    code = "ResourceNotFoundException"
