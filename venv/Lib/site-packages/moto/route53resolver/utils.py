"""Pagination control model for Route53Resolver."""

from .exceptions import InvalidNextTokenException

PAGINATION_MODEL = {
    "list_resolver_endpoint_ip_addresses": {
        "input_token": "next_token",
        "limit_key": "max_results",
        "limit_default": 100,
        "unique_attribute": "IpId",
    },
    "list_resolver_endpoints": {
        "input_token": "next_token",
        "limit_key": "max_results",
        "limit_default": 100,
        "unique_attribute": "id",
        "fail_on_invalid_token": InvalidNextTokenException,
    },
    "list_resolver_rules": {
        "input_token": "next_token",
        "limit_key": "max_results",
        "limit_default": 100,
        "unique_attribute": "id",
    },
    "list_resolver_rule_associations": {
        "input_token": "next_token",
        "limit_key": "max_results",
        "limit_default": 100,
        "unique_attribute": "id",
    },
    "list_tags_for_resource": {
        "input_token": "next_token",
        "limit_key": "max_results",
        "limit_default": 100,
        "unique_attribute": "Key",
        "fail_on_invalid_token": InvalidNextTokenException,
    },
}
