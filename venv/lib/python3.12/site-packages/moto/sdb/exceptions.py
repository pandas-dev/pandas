"""Exceptions raised by the sdb service."""

from moto.core.exceptions import ServiceException


class SDBError(ServiceException):
    box_usage = 0.0055590278


class InvalidParameterError(SDBError):
    code = "InvalidParameterValue"


class InvalidDomainName(InvalidParameterError):
    def __init__(self, domain_name: str):
        super().__init__(f"Value ({domain_name}) for parameter DomainName is invalid. ")


class UnknownDomainName(SDBError):
    code = "NoSuchDomain"
    message = "The specified domain does not exist."
