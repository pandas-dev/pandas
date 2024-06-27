"""Exceptions raised by the acmpca service."""

from moto.core.exceptions import JsonRESTError


class ResourceNotFoundException(JsonRESTError):
    def __init__(self, arn: str):
        super().__init__("ResourceNotFoundException", f"Resource {arn} not found")


class InvalidS3ObjectAclInCrlConfiguration(JsonRESTError):
    code = 400

    def __init__(self, value: str):
        super().__init__(
            "InvalidS3ObjectAclInCrlConfiguration",
            f"Invalid value for parameter RevocationConfiguration.CrlConfiguration.S3ObjectAcl, value: {value}, valid values: ['PUBLIC_READ', 'BUCKET_OWNER_FULL_CONTROL']",
        )
