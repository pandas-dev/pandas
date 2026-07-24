"""Exceptions raised by the emrcontainers service."""

from moto.core.exceptions import ServiceException


class EMRContainersException(ServiceException):
    pass


class ResourceNotFoundException(EMRContainersException):
    code = "ResourceNotFoundException"
