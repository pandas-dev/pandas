"""Exceptions raised by the sdb service."""

from typing import Any

from moto.core.exceptions import RESTError

SDB_ERROR = """<?xml version="1.0"?>
<Response>
    <Errors>
        <Error>
            <Code>{{ error_type }}</Code>
            <Message>{{ message }}</Message>
            <BoxUsage>0.0055590278</BoxUsage>
        </Error>
    </Errors>
    <RequestID>ba3a8c86-dc37-0a45-ef44-c6cf7876a62f</RequestID>
</Response>"""


class InvalidParameterError(RESTError):
    code = 400
    extended_templates = {"sdb_error": SDB_ERROR}
    env = RESTError.extended_environment(extended_templates)

    def __init__(self, **kwargs: Any):
        kwargs.setdefault("template", "sdb_error")
        kwargs["error_type"] = "InvalidParameterValue"
        super().__init__(**kwargs)


class InvalidDomainName(InvalidParameterError):
    code = 400

    def __init__(self, domain_name: str):
        super().__init__(
            message=f"Value ({domain_name}) for parameter DomainName is invalid. "
        )


class UnknownDomainName(RESTError):
    code = 400
    extended_templates = {"sdb_error": SDB_ERROR}
    env = RESTError.extended_environment(extended_templates)

    def __init__(self, **kwargs: Any):
        kwargs.setdefault("template", "sdb_error")
        kwargs["error_type"] = "NoSuchDomain"
        kwargs["message"] = "The specified domain does not exist."
        super().__init__(**kwargs)
