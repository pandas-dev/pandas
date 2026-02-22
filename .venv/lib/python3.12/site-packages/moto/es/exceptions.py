"""Exceptions raised by the ElasticSearch service."""

from moto.core.exceptions import JsonRESTError


class ElasticSearchError(JsonRESTError):
    code = 400


class ResourceNotFound(ElasticSearchError):
    code = 409

    def __init__(self, resource_type: str, resource_name: str):
        msg = f"{resource_type} not found: {resource_name}"
        super().__init__("ResourceNotFoundException", msg)


class InvalidDomainName(ElasticSearchError):
    def __init__(self, domain_name: str):
        msg = f"1 validation error detected: Value '{domain_name}' at 'domainName' failed to satisfy constraint: Member must satisfy regular expression pattern: [a-z][a-z0-9\\-]+"
        super().__init__("ValidationException", msg)


class DomainNotFound(ResourceNotFound):
    def __init__(self, domain_name: str):
        super().__init__("Domain", domain_name)
