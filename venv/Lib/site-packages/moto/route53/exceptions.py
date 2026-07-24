"""Exceptions raised by the Route53 service."""

from moto.core.exceptions import ServiceException


class Route53ClientError(ServiceException):
    """Base class for Route53 errors."""

    pass


class ConflictingDomainExists(Route53ClientError):
    """Domain already exists."""

    def __init__(self, domain_name: str, delegation_set_id: str | None) -> None:
        message = (
            f"Cannot create hosted zone with DelegationSetId DelegationSetId:{delegation_set_id} as the DNSName"
            f"{domain_name} conflicts with existing ones sharing the delegation set"
        )
        super().__init__("ConflictingDomainExists", message)


class InvalidInput(Route53ClientError):
    """Malformed ARN for the CloudWatch log group."""

    def __init__(self, message: str):
        super().__init__("InvalidInput", message)


class UnsupportedCharacter(Route53ClientError):
    """Malformed ARN for the CloudWatch log group."""

    def __init__(self, code: str, char: str):
        super().__init__(
            code,
            f"FATAL problem: UnsupportedCharacter (Value contains unsupported characters) encountered with '{char}'",
        )


class InvalidCloudWatchArn(InvalidInput):
    def __init__(self) -> None:
        message = "The ARN for the CloudWatch Logs log group is invalid"
        super().__init__(message)


class InvalidActionValue(InvalidInput):
    def __init__(self, value: str):
        message = (
            f"Invalid XML ; cvc-enumeration-valid: Value '{value}' is not facet-valid"
            " with respect to enumeration '[CREATE, DELETE, UPSERT]'. It must be a value from the enumeration."
        )
        super().__init__(message)


class InvalidPaginationToken(Route53ClientError):
    """Bad NextToken specified when listing query logging configs."""

    def __init__(self) -> None:
        message = (
            "Route 53 can't get the next page of query logging configurations "
            "because the specified value for NextToken is invalid."
        )
        super().__init__("InvalidPaginationToken", message)


class InvalidVPCId(Route53ClientError):
    """Missing/Invalid VPC ID"""

    def __init__(self) -> None:
        message = "Invalid or missing VPC Id."
        super().__init__("InvalidVPCId", message)


class NoSuchCloudWatchLogsLogGroup(Route53ClientError):
    """CloudWatch LogGroup has a permissions policy, but does not exist."""

    def __init__(self) -> None:
        message = "The specified CloudWatch Logs log group doesn't exist."
        super().__init__("NoSuchCloudWatchLogsLogGroup", message)


class NoSuchHostedZone(Route53ClientError):
    """HostedZone does not exist."""

    def __init__(self, host_zone_id: str):
        message = f"No hosted zone found with ID: {host_zone_id}"
        super().__init__("NoSuchHostedZone", message)


class NoSuchHealthCheck(Route53ClientError):
    """HealthCheck does not exist."""

    def __init__(self, health_check_id: str):
        message = f"A health check with id {health_check_id} does not exist."
        super().__init__("NoSuchHealthCheck", message)


class HostedZoneNotEmpty(Route53ClientError):
    """HostedZone does not exist."""

    def __init__(self) -> None:
        message = (
            "The hosted zone contains resource records that are not SOA or NS records."
        )
        super().__init__("HostedZoneNotEmpty", message)


class PublicZoneVPCAssociation(Route53ClientError):
    """Public hosted zone can't be associated."""

    def __init__(self) -> None:
        message = "You're trying to associate a VPC with a public hosted zone. Amazon Route 53 doesn't support associating a VPC with a public hosted zone."
        super().__init__("PublicZoneVPCAssociation", message)


class LastVPCAssociation(Route53ClientError):
    """Last VPC can't be disassociate."""

    def __init__(self) -> None:
        message = "The VPC that you're trying to disassociate from the private hosted zone is the last VPC that is associated with the hosted zone. Amazon Route 53 doesn't support disassociating the last VPC from a hosted zone."
        super().__init__("LastVPCAssociation", message)


class NoSuchQueryLoggingConfig(Route53ClientError):
    """Query log config does not exist."""

    def __init__(self) -> None:
        message = "The query logging configuration does not exist"
        super().__init__("NoSuchQueryLoggingConfig", message)


class QueryLoggingConfigAlreadyExists(Route53ClientError):
    """Query log config exists for the hosted zone."""

    def __init__(self) -> None:
        message = "A query logging configuration already exists for this hosted zone"
        super().__init__("QueryLoggingConfigAlreadyExists", message)


class InvalidChangeBatch(Route53ClientError):
    def __init__(self) -> None:
        message = "Number of records limit of 1000 exceeded."
        super().__init__("InvalidChangeBatch", message)


class NoSuchDelegationSet(Route53ClientError):
    def __init__(self, delegation_set_id: str):
        super().__init__("NoSuchDelegationSet", delegation_set_id)


class DnsNameInvalidForZone(Route53ClientError):
    def __init__(self, name: str, zone_name: str):
        error_msg = (
            f"""RRSet with DNS name {name} is not permitted in zone {zone_name}"""
        )
        super().__init__("InvalidChangeBatch", error_msg)


class ResourceRecordAlreadyExists(Route53ClientError):
    def __init__(self, name: str, _type: str):
        super().__init__(
            "InvalidChangeBatch",
            f"Tried to create resource record set [name='{name}', type='{_type}'] but it already exists",
        )
