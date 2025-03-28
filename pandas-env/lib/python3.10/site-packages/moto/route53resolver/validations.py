"""Route53ResolverBackend validations that result in ValidationException.

Note that ValidationExceptions are accumulative.
"""

import re
from typing import Any, Dict, List, Optional, Tuple

from moto.route53resolver.exceptions import RRValidationException


def validate_args(validators: List[Tuple[str, Any]]) -> None:
    """Raise exception if any of the validations fails.

    validators is a list of tuples each containing the following:
        (printable field name, field value)

    The error messages are accumulated before the exception is raised.
    """
    validation_map = {
        "creatorRequestId": validate_creator_request_id,
        "direction": validate_direction,
        "domainName": validate_domain_name,
        "ipAddresses": validate_ip_addresses,
        "ipAddresses.subnetId": validate_subnets,
        "maxResults": validate_max_results,
        "name": validate_name,
        "resolverEndpointId": validate_endpoint_id,
        "resolverRuleAssociationId": validate_rule_association_id,
        "resolverRuleId": validate_rule_id,
        "ruleType": validate_rule_type,
        "securityGroupIds": validate_security_group_ids,
        "targetIps.port": validate_target_port,
        "vPCId": validate_vpc_id,
    }

    err_msgs = []
    # This eventually could be a switch (python 3.10), eliminating the need
    # for the above map and individual functions.
    for fieldname, value in validators:
        msg = validation_map[fieldname](value)  # type: ignore
        if msg:
            err_msgs.append((fieldname, value, msg))
    if err_msgs:
        raise RRValidationException(err_msgs)


def validate_creator_request_id(value: Optional[str]) -> str:
    """Raise exception if the creator_request_id has invalid length."""
    if value and len(value) > 255:
        return "have length less than or equal to 255"
    return ""


def validate_direction(value: Optional[str]) -> str:
    """Raise exception if direction not one of the allowed values."""
    if value and value not in ["INBOUND", "OUTBOUND"]:
        return "satisfy enum value set: [INBOUND, OUTBOUND]"
    return ""


def validate_domain_name(value: str) -> str:
    """Raise exception if the domain_name has invalid length."""
    if len(value) > 256:
        return "have length less than or equal to 256"
    return ""


def validate_endpoint_id(value: Optional[str]) -> str:
    """Raise exception if resolver endpoint id has invalid length."""
    if value and len(value) > 64:
        return "have length less than or equal to 64"
    return ""


def validate_ip_addresses(value: str) -> str:
    """Raise exception if IPs fail to match length constraint."""
    if len(value) > 10:
        return "have length less than or equal to 10"
    return ""


def validate_max_results(value: Optional[int]) -> str:
    """Raise exception if number of endpoints or IPs is too large."""
    if value and value > 100:
        return "have length less than or equal to 100"
    return ""


def validate_name(value: Optional[str]) -> str:
    """Raise exception if name fails to match constraints."""
    if value:
        if len(value) > 64:
            return "have length less than or equal to 64"
        name_pattern = r"^(?!^[0-9]+$)([a-zA-Z0-9-_' ']+)$"
        if not re.match(name_pattern, value):
            return rf"satisfy regular expression pattern: {name_pattern}"
    return ""


def validate_rule_association_id(value: Optional[str]) -> str:
    """Raise exception if resolver rule association id has invalid length."""
    if value and len(value) > 64:
        return "have length less than or equal to 64"
    return ""


def validate_rule_id(value: Optional[str]) -> str:
    """Raise exception if resolver rule id has invalid length."""
    if value and len(value) > 64:
        return "have length less than or equal to 64"
    return ""


def validate_rule_type(value: Optional[str]) -> str:
    """Raise exception if rule_type not one of the allowed values."""
    if value and value not in ["FORWARD", "SYSTEM", "RECURSIVE"]:
        return "satisfy enum value set: [FORWARD, SYSTEM, RECURSIVE]"
    return ""


def validate_security_group_ids(value: List[str]) -> str:
    """Raise exception if IPs fail to match length constraint."""
    # Too many security group IDs is an InvalidParameterException.
    for group_id in value:
        if len(group_id) > 64:
            return (
                "have length less than or equal to 64 and Member must have "
                "length greater than or equal to 1"
            )
    return ""


def validate_subnets(value: List[Dict[str, Any]]) -> str:
    """Raise exception if subnets fail to match length constraint."""
    for subnet_id in [x["SubnetId"] for x in value]:
        if len(subnet_id) > 32:
            return "have length less than or equal to 32"
    return ""


def validate_target_port(value: Optional[Dict[str, int]]) -> str:
    """Raise exception if target port fails to match length constraint."""
    if value and value["Port"] > 65535:
        return "have value less than or equal to 65535"
    return ""


def validate_vpc_id(value: str) -> str:
    """Raise exception if VPC id has invalid length."""
    if len(value) > 64:
        return "have length less than or equal to 64"
    return ""
