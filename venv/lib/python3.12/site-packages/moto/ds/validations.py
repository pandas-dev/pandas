"""DirectoryServiceBackend checks that result in ValidationException.

Note that ValidationExceptions are accumulative.
"""

import re
from typing import Any

from moto.ds.exceptions import DsValidationException


def validate_args(validators: Any) -> None:
    """Raise exception if any of the validations fails.

    validators is a list of tuples each containing the following:
        (printable field name, field value)

    The error messages are accumulated before the exception is raised.
    """
    validation_map = {
        "alias": validate_alias,
        "description": validate_description,
        "directoryId": validate_directory_id,
        "connectSettings.customerDnsIps": validate_dns_ips,
        "edition": validate_edition,
        "name": validate_name,
        "password": validate_password,
        "shortName": validate_short_name,
        "size": validate_size,
        "ssoPassword": validate_sso_password,
        "connectSettings.vpcSettings.subnetIds": validate_subnet_ids,
        "connectSettings.customerUserName": validate_user_name,
        "userName": validate_user_name,
        "vpcSettings.subnetIds": validate_subnet_ids,
        "trustDirection": validate_trust_direction,
        "remoteDomainName": validate_remote_domain_name,
    }
    err_msgs = []
    # This eventually could be a switch (python 3.10), elminating the need
    # for the above map and individual functions.
    for fieldname, value in validators:
        msg = validation_map[fieldname](value)
        if msg:
            err_msgs.append((fieldname, value, msg))
    if err_msgs:
        raise DsValidationException(err_msgs)


def validate_alias(value: str) -> str:
    """Raise exception if alias fails to conform to length and constraints."""
    if len(value) > 62:
        return "have length less than or equal to 62"

    alias_pattern = r"^(?!D-|d-)([\da-zA-Z]+)([-]*[\da-zA-Z])*$"
    if not re.match(alias_pattern, value):
        return rf"satisfy regular expression pattern: {alias_pattern}"
    return ""


def validate_description(value: str) -> str:
    """Raise exception if description exceeds length."""
    if value and len(value) > 128:
        return "have length less than or equal to 128"
    return ""


def validate_directory_id(value: str) -> str:
    """Raise exception if the directory id is invalid."""
    id_pattern = r"^d-[0-9a-f]{10}$"
    if not re.match(id_pattern, value):
        return rf"satisfy regular expression pattern: {id_pattern}"
    return ""


def validate_dns_ips(value: str) -> str:
    """Raise exception if DNS IPs fail to match constraints."""
    dnsip_pattern = (
        r"^(?:(?:25[0-5]|2[0-4][0-9]|[01]?[0-9][0-9]?)\.){3}"
        r"(?:25[0-5]|2[0-4][0-9]|[01]?[0-9][0-9]?)$"
    )
    for dnsip in value:
        if not re.match(dnsip_pattern, dnsip):
            return rf"satisfy regular expression pattern: {dnsip_pattern}"
    return ""


def validate_edition(value: str) -> str:
    """Raise exception if edition not one of the allowed values."""
    if value and value not in ["Enterprise", "Standard"]:
        return "satisfy enum value set: [Enterprise, Standard]"
    return ""


def validate_name(value: str) -> str:
    """Raise exception if name fails to match constraints."""
    name_pattern = r"^([a-zA-Z0-9]+[\.-])+([a-zA-Z0-9])+$"
    if not re.match(name_pattern, value):
        return rf"satisfy regular expression pattern: {name_pattern}"
    return ""


def validate_password(value: str) -> str:
    """Raise exception if password fails to match constraints."""
    passwd_pattern = (
        r"^(?=^.{8,64}$)((?=.*\d)(?=.*[A-Z])(?=.*[a-z])|"
        r"(?=.*\d)(?=.*[^A-Za-z0-9\s])(?=.*[a-z])|"
        r"(?=.*[^A-Za-z0-9\s])(?=.*[A-Z])(?=.*[a-z])|"
        r"(?=.*\d)(?=.*[A-Z])(?=.*[^A-Za-z0-9\s]))^.*$"
    )
    if not re.match(passwd_pattern, value):
        return rf"satisfy regular expression pattern: {passwd_pattern}"
    return ""


def validate_short_name(value: str) -> str:
    """Raise exception if short name fails to match constraints."""
    short_name_pattern = r'^[^\/:*?"<>|.]+[^\/:*?"<>|]*$'
    if value and not re.match(short_name_pattern, value):
        return rf"satisfy regular expression pattern: {short_name_pattern}"
    return ""


def validate_size(value: str) -> str:
    """Raise exception if size fails to match constraints."""
    if value.lower() not in ["small", "large"]:
        return "satisfy enum value set: [Small, Large]"
    return ""


def validate_sso_password(value: str) -> str:
    """Raise exception is SSO password exceeds length."""
    if value and len(value) > 128:
        return "have length less than or equal to 128"
    return ""


def validate_subnet_ids(value: str) -> str:
    """Raise exception is subnet IDs fail to match constraints."""
    subnet_id_pattern = r"^(subnet-[0-9a-f]{8}|subnet-[0-9a-f]{17})$"
    for subnet in value:
        if not re.match(subnet_id_pattern, subnet):
            return rf"satisfy regular expression pattern: {subnet_id_pattern}"
    return ""


def validate_user_name(value: str) -> str:
    """Raise exception is username fails to match constraints."""
    username_pattern = r"^[a-zA-Z0-9._-]+$"
    if value and not re.match(username_pattern, value):
        return rf"satisfy regular expression pattern: {username_pattern}"
    return ""


def validate_trust_direction(value: str) -> str:
    """Raise exception if trust direction fails to match constraints."""
    if value not in ["One-Way: Outgoing", "One-Way: Incoming", "Two-Way"]:
        return "satisfy enum value set: [One-Way: Outgoing, One-Way: Incoming, Two-Way]"
    return ""


def validate_remote_domain_name(value: str) -> str:
    """Raise exception if remote domain name fails to match constraints."""
    domain_name_pattern = r"^([a-zA-Z0-9]+[\\.-])+([a-zA-Z0-9])+[.]?$"
    if not re.match(domain_name_pattern, value):
        return rf"satisfy regular expression pattern: {domain_name_pattern}"
    elif len(value) > 1024:
        return "have length less than or equal to 1024"
    return ""
