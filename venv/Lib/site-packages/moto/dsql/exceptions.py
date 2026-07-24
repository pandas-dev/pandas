"""Exceptions raised by the dsql service."""

from moto.core.exceptions import ServiceException


class ResourceNotFoundException(ServiceException):
    code = "ResourceNotFoundException"

    def __init__(self, arn: str, resource_id: str, resource_type: str):
        self.resource_id = resource_id
        self.resource_type = resource_type
        message = (
            f"The resource with ARN {arn} doesn't exist. Verify the ARN and try again."
        )
        super().__init__(message)
