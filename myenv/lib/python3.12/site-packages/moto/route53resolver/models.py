"""Route53ResolverBackend class with methods for supported APIs."""

import re
from collections import defaultdict
from datetime import datetime, timezone
from ipaddress import IPv4Address, ip_address, ip_network
from typing import Any, Dict, List, Optional, Set

from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import BaseModel
from moto.ec2 import ec2_backends
from moto.ec2.exceptions import InvalidSecurityGroupNotFoundError, InvalidSubnetIdError
from moto.moto_api._internal import mock_random
from moto.route53resolver.exceptions import (
    InvalidParameterException,
    InvalidRequestException,
    LimitExceededException,
    ResourceExistsException,
    ResourceInUseException,
    ResourceNotFoundException,
    TagValidationException,
)
from moto.route53resolver.utils import PAGINATION_MODEL
from moto.route53resolver.validations import validate_args
from moto.utilities.paginator import paginate
from moto.utilities.tagging_service import TaggingService
from moto.utilities.utils import get_partition

CAMEL_TO_SNAKE_PATTERN = re.compile(r"(?<!^)(?=[A-Z])")


class ResolverRuleAssociation(BaseModel):  # pylint: disable=too-few-public-methods
    """Representation of a fake Route53 Resolver Rules Association."""

    MAX_TAGS_PER_RESOLVER_ENDPOINT = 200
    MAX_RULE_ASSOCIATIONS_PER_REGION = 2000

    # There are two styles of filter names and either will be transformed
    # into lowercase snake.
    FILTER_NAMES = [
        "name",
        "resolver_rule_id",
        "status",
        "vpc_id",
    ]

    def __init__(
        self,
        region: str,
        resolver_rule_association_id: str,
        resolver_rule_id: str,
        vpc_id: str,
        name: str,
    ):  # pylint: disable=too-many-arguments
        self.region = region
        self.resolver_rule_id = resolver_rule_id
        self.name = name
        self.vpc_id = vpc_id

        # Constructed members.
        self.id = resolver_rule_association_id  # pylint: disable=invalid-name
        self.status = "COMPLETE"
        self.status_message = ""

    def description(self) -> Dict[str, Any]:
        """Return dictionary of relevant info for resolver rule association."""
        return {
            "Id": self.id,
            "ResolverRuleId": self.resolver_rule_id,
            "Name": self.name,
            "VPCId": self.vpc_id,
            "Status": self.status,
            "StatusMessage": self.status_message,
        }


class ResolverRule(BaseModel):  # pylint: disable=too-many-instance-attributes
    """Representation of a fake Route53 Resolver Rule."""

    MAX_TAGS_PER_RESOLVER_RULE = 200
    MAX_RULES_PER_REGION = 1000

    # There are two styles of filter names and either will be transformed
    # into lowercase snake.
    FILTER_NAMES = [
        "creator_request_id",
        "domain_name",
        "name",
        "resolver_endpoint_id",
        "status",
        "rule_type",  # actual filter is "Type"
    ]

    def __init__(
        self,
        account_id: str,
        region: str,
        rule_id: str,
        creator_request_id: str,
        rule_type: str,
        domain_name: str,
        target_ips: Optional[List[Dict[str, Any]]],
        resolver_endpoint_id: Optional[str],
        name: str,
    ):  # pylint: disable=too-many-arguments
        self.account_id = account_id
        self.region = region
        self.creator_request_id = creator_request_id
        self.name = name
        self.rule_id = rule_id
        self.rule_type = rule_type
        self.domain_name = domain_name + "."
        self.target_ips = target_ips
        self.resolver_endpoint_id = resolver_endpoint_id

        # Constructed members.
        self.id = rule_id  # pylint: disable=invalid-name
        self.status = "COMPLETE"

        # The status message should contain a trace Id which is the value
        # of X-Amzn-Trace-Id.  We don't have that info, so a random number
        # of similar format and length will be used.
        self.status_message = (
            f"[Trace id: 1-{mock_random.get_random_hex(8)}-{mock_random.get_random_hex(24)}] "
            f"Successfully created Resolver Rule"
        )
        self.share_status = "NOT_SHARED"
        self.creation_time = datetime.now(timezone.utc).isoformat()
        self.modification_time = datetime.now(timezone.utc).isoformat()

    @property
    def arn(self) -> str:
        return f"arn:{get_partition(self.region)}:route53resolver:{self.region}:{self.account_id}:resolver-rule/{self.id}"

    def description(self) -> Dict[str, Any]:
        """Return a dictionary of relevant info for this resolver rule."""
        return {
            "Id": self.id,
            "CreatorRequestId": self.creator_request_id,
            "Arn": self.arn,
            "DomainName": self.domain_name,
            "Status": self.status,
            "StatusMessage": self.status_message,
            "RuleType": self.rule_type,
            "Name": self.name,
            "TargetIps": self.target_ips,
            "ResolverEndpointId": self.resolver_endpoint_id,
            "OwnerId": self.account_id,
            "ShareStatus": self.share_status,
            "CreationTime": self.creation_time,
            "ModificationTime": self.modification_time,
        }


class ResolverEndpoint(BaseModel):  # pylint: disable=too-many-instance-attributes
    """Representation of a fake Route53 Resolver Endpoint."""

    MAX_TAGS_PER_RESOLVER_ENDPOINT = 200
    MAX_ENDPOINTS_PER_REGION = 4

    # There are two styles of filter names and either will be transformed
    # into lowercase snake.
    FILTER_NAMES = [
        "creator_request_id",
        "direction",
        "host_vpc_id",
        "ip_address_count",
        "name",
        "security_group_ids",
        "status",
    ]

    def __init__(
        self,
        account_id: str,
        region: str,
        endpoint_id: str,
        creator_request_id: str,
        security_group_ids: List[str],
        direction: str,
        ip_addresses: List[Dict[str, Any]],
        name: str,
    ):  # pylint: disable=too-many-arguments
        self.account_id = account_id
        self.region = region
        self.creator_request_id = creator_request_id
        self.name = name
        self.security_group_ids = security_group_ids
        self.direction = direction
        self.ip_addresses = ip_addresses
        self.ec2_backend = ec2_backends[self.account_id][self.region]

        # Constructed members.
        self.id = endpoint_id  # pylint: disable=invalid-name

        # NOTE; This currently doesn't reflect IPv6 addresses.
        self.subnets = self._build_subnet_info()
        self.eni_ids = self.create_eni()
        self.ip_address_count = len(ip_addresses)

        self.host_vpc_id = self._vpc_id_from_subnet()
        self.status = "OPERATIONAL"

        # The status message should contain a trace Id which is the value
        # of X-Amzn-Trace-Id.  We don't have that info, so a random number
        # of similar format and length will be used.
        self.status_message = (
            f"[Trace id: 1-{mock_random.get_random_hex(8)}-{mock_random.get_random_hex(24)}] "
            f"Successfully created Resolver Endpoint"
        )
        self.creation_time = datetime.now(timezone.utc).isoformat()
        self.modification_time = datetime.now(timezone.utc).isoformat()

    @property
    def arn(self) -> str:
        return f"arn:{get_partition(self.region)}:route53resolver:{self.region}:{self.account_id}:resolver-endpoint/{self.id}"

    def _vpc_id_from_subnet(self) -> str:
        """Return VPC Id associated with the subnet.

        The assumption is that all of the subnets are associated with the
        same VPC.  We don't check that assumption, but otherwise the existence
        of the subnets has already been checked.
        """
        first_subnet_id = self.ip_addresses[0]["SubnetId"]
        subnet_info = self.ec2_backend.describe_subnets(subnet_ids=[first_subnet_id])[0]
        return subnet_info.vpc_id

    def _build_subnet_info(self) -> Dict[str, Any]:
        """Create a dict of subnet info, including ip addrs and ENI ids.

        self.subnets[subnet_id][ip_addr1] = eni-id1 ...
        """
        subnets: Dict[str, Any] = defaultdict(dict)
        for entry in self.ip_addresses:
            subnets[entry["SubnetId"]][entry["Ip"]] = (
                f"rni-{mock_random.get_random_hex(17)}"
            )
        return subnets

    def create_eni(self) -> List[str]:
        """Create a VPC ENI for each combo of AZ, subnet and IP."""
        eni_ids = []
        for subnet, ip_info in self.subnets.items():
            for ip_addr, eni_id in ip_info.items():
                eni_info = self.ec2_backend.create_network_interface(
                    description=f"Route 53 Resolver: {self.id}:{eni_id}",
                    group_ids=self.security_group_ids,
                    interface_type="interface",
                    private_ip_address=ip_addr,
                    private_ip_addresses=[
                        {"Primary": True, "PrivateIpAddress": ip_addr}
                    ],
                    subnet=subnet,
                )
                eni_ids.append(eni_info.id)
        return eni_ids

    def delete_eni(self) -> None:
        """Delete the VPC ENI created for the subnet and IP combos."""
        for eni_id in self.eni_ids:
            self.ec2_backend.delete_network_interface(eni_id)

    def description(self) -> Dict[str, Any]:
        """Return a dictionary of relevant info for this resolver endpoint."""
        return {
            "Id": self.id,
            "CreatorRequestId": self.creator_request_id,
            "Arn": self.arn,
            "Name": self.name,
            "SecurityGroupIds": self.security_group_ids,
            "Direction": self.direction,
            "IpAddressCount": self.ip_address_count,
            "HostVPCId": self.host_vpc_id,
            "Status": self.status,
            "StatusMessage": self.status_message,
            "CreationTime": self.creation_time,
            "ModificationTime": self.modification_time,
        }

    def ip_descriptions(self) -> List[Dict[str, Any]]:
        """Return a list of dicts describing resolver endpoint IP addresses."""
        description = []
        for subnet_id, ip_info in self.subnets.items():
            for ip_addr, eni_id in ip_info.items():
                description.append(
                    {
                        "IpId": eni_id,
                        "SubnetId": subnet_id,
                        "Ip": ip_addr,
                        "Status": "ATTACHED",
                        "StatusMessage": "This IP address is operational.",
                        "CreationTime": self.creation_time,
                        "ModificationTime": self.modification_time,
                    }
                )
        return description

    def update_name(self, name: str) -> None:
        """Replace existing name with new name."""
        self.name = name
        self.modification_time = datetime.now(timezone.utc).isoformat()

    def associate_ip_address(self, value: Dict[str, Any]) -> None:
        self.ip_addresses.append(value)
        self.ip_address_count = len(self.ip_addresses)

        eni_id = f"rni-{mock_random.get_random_hex(17)}"
        self.subnets[value["SubnetId"]][value["Ip"]] = eni_id

        eni_info = self.ec2_backend.create_network_interface(
            description=f"Route 53 Resolver: {self.id}:{eni_id}",
            group_ids=self.security_group_ids,
            interface_type="interface",
            private_ip_address=value.get("Ip"),  # type: ignore[arg-type]
            private_ip_addresses=[
                {"Primary": True, "PrivateIpAddress": value.get("Ip")}
            ],
            subnet=value.get("SubnetId"),
        )
        self.eni_ids.append(eni_info.id)

    def disassociate_ip_address(self, value: Dict[str, Any]) -> None:
        if not value.get("Ip") and value.get("IpId"):
            for ip_addr, eni_id in self.subnets[value.get("SubnetId")].items():  # type: ignore
                if value.get("IpId") == eni_id:
                    value["Ip"] = ip_addr
        if value.get("Ip"):
            self.ip_addresses = list(
                filter(lambda i: i["Ip"] != value.get("Ip"), self.ip_addresses)
            )

            if len(self.subnets[value["SubnetId"]]) == 1:
                self.subnets.pop(value["SubnetId"])
            else:
                self.subnets[value["SubnetId"]].pop(value["Ip"])
            for eni_id in self.eni_ids:
                eni_info = self.ec2_backend.get_network_interface(eni_id)
                if eni_info.private_ip_address == value.get("Ip"):
                    self.ec2_backend.delete_network_interface(eni_id)
                    self.eni_ids.remove(eni_id)
            self.ip_address_count = len(self.ip_addresses)


class Route53ResolverBackend(BaseBackend):
    """Implementation of Route53Resolver APIs."""

    def __init__(self, region_name: str, account_id: str):
        super().__init__(region_name, account_id)
        self.resolver_endpoints: Dict[
            str, ResolverEndpoint
        ] = {}  # Key is self-generated ID (endpoint_id)
        self.resolver_rules: Dict[
            str, ResolverRule
        ] = {}  # Key is self-generated ID (rule_id)
        self.resolver_rule_associations: Dict[
            str, ResolverRuleAssociation
        ] = {}  # Key is resolver_rule_association_id)
        self.tagger = TaggingService()

        self.ec2_backend = ec2_backends[self.account_id][self.region_name]

    @staticmethod
    def default_vpc_endpoint_service(
        service_region: str, zones: List[str]
    ) -> List[Dict[str, str]]:
        """List of dicts representing default VPC endpoints for this service."""
        return BaseBackend.default_vpc_endpoint_service_factory(
            service_region, zones, "route53resolver"
        )

    def associate_resolver_rule(
        self, resolver_rule_id: str, name: str, vpc_id: str
    ) -> ResolverRuleAssociation:
        validate_args(
            [("resolverRuleId", resolver_rule_id), ("name", name), ("vPCId", vpc_id)]
        )

        associations = [
            x
            for x in self.resolver_rule_associations.values()
            if x.region == self.region_name
        ]
        if len(associations) > ResolverRuleAssociation.MAX_RULE_ASSOCIATIONS_PER_REGION:
            # This error message was not verified to be the same for AWS.
            raise LimitExceededException(
                f"Account '{self.account_id}' has exceeded 'max-rule-association'"
            )

        if resolver_rule_id not in self.resolver_rules:
            raise ResourceNotFoundException(
                f"Resolver rule with ID '{resolver_rule_id}' does not exist."
            )

        vpcs = self.ec2_backend.describe_vpcs()
        if vpc_id not in [x.id for x in vpcs]:
            raise InvalidParameterException(f"The vpc ID '{vpc_id}' does not exist")

        # Can't duplicate resolver rule, vpc id associations.
        for association in self.resolver_rule_associations.values():
            if (
                resolver_rule_id == association.resolver_rule_id
                and vpc_id == association.vpc_id
            ):
                raise InvalidRequestException(
                    f"Cannot associate rules with same domain name with same "
                    f"VPC. Conflict with resolver rule '{resolver_rule_id}'"
                )

        rule_association_id = f"rslvr-rrassoc-{mock_random.get_random_hex(17)}"
        rule_association = ResolverRuleAssociation(
            self.region_name, rule_association_id, resolver_rule_id, vpc_id, name
        )
        self.resolver_rule_associations[rule_association_id] = rule_association
        return rule_association

    def _verify_subnet_ips(
        self, ip_addresses: List[Dict[str, Any]], initial: bool = True
    ) -> None:
        """
        Perform additional checks on the IPAddresses.

        NOTE: This does not include IPv6 addresses.
        """
        # only initial endpoint creation requires atleast two ip addresses
        if initial:
            if len(ip_addresses) < 2:
                raise InvalidRequestException(
                    "Resolver endpoint needs to have at least 2 IP addresses"
                )

        subnets: Dict[str, Set[str]] = defaultdict(set)
        for subnet_id, ip_addr in [(x["SubnetId"], x["Ip"]) for x in ip_addresses]:
            try:
                subnet_info = self.ec2_backend.describe_subnets(subnet_ids=[subnet_id])[
                    0
                ]
            except InvalidSubnetIdError as exc:
                raise InvalidParameterException(
                    f"The subnet ID '{subnet_id}' does not exist"
                ) from exc

            # IP in IPv4 CIDR range and not reserved?
            if ip_address(ip_addr) in subnet_info.reserved_ips or ip_address(
                ip_addr
            ) not in ip_network(subnet_info.cidr_block):
                raise InvalidRequestException(
                    f"IP address '{ip_addr}' is either not in subnet "
                    f"'{subnet_id}' CIDR range or is reserved"
                )

            if ip_addr in subnets[subnet_id]:
                raise ResourceExistsException(
                    f"The IP address '{ip_addr}' in subnet '{subnet_id}' is already in use"
                )
            subnets[subnet_id].add(ip_addr)

    def _verify_security_group_ids(self, security_group_ids: List[str]) -> None:
        """Perform additional checks on the security groups."""
        if len(security_group_ids) > 10:
            raise InvalidParameterException("Maximum of 10 security groups are allowed")

        for group_id in security_group_ids:
            if not group_id.startswith("sg-"):
                raise InvalidParameterException(
                    f"Malformed security group ID: Invalid id: '{group_id}' "
                    f"(expecting 'sg-...')"
                )
            try:
                self.ec2_backend.describe_security_groups(group_ids=[group_id])
            except InvalidSecurityGroupNotFoundError as exc:
                raise ResourceNotFoundException(
                    f"The security group '{group_id}' does not exist"
                ) from exc

    def create_resolver_endpoint(
        self,
        region: str,
        creator_request_id: str,
        name: str,
        security_group_ids: List[str],
        direction: str,
        ip_addresses: List[Dict[str, Any]],
        tags: List[Dict[str, str]],
    ) -> ResolverEndpoint:  # pylint: disable=too-many-arguments
        """
        Return description for a newly created resolver endpoint.

        NOTE:  IPv6 IPs are currently not being filtered when
        calculating the create_resolver_endpoint() IpAddresses.
        """
        validate_args(
            [
                ("creatorRequestId", creator_request_id),
                ("direction", direction),
                ("ipAddresses", ip_addresses),
                ("name", name),
                ("securityGroupIds", security_group_ids),
                ("ipAddresses.subnetId", ip_addresses),
            ]
        )
        errmsg = self.tagger.validate_tags(
            tags or [], limit=ResolverEndpoint.MAX_TAGS_PER_RESOLVER_ENDPOINT
        )
        if errmsg:
            raise TagValidationException(errmsg)

        endpoints = [x for x in self.resolver_endpoints.values() if x.region == region]
        if len(endpoints) > ResolverEndpoint.MAX_ENDPOINTS_PER_REGION:
            raise LimitExceededException(
                f"Account '{self.account_id}' has exceeded 'max-endpoints'"
            )

        for x in ip_addresses:
            if not x.get("Ip"):
                subnet_info = self.ec2_backend.describe_subnets(
                    subnet_ids=[x["SubnetId"]]
                )[0]
                x["Ip"] = subnet_info.get_available_subnet_ip(self)  # type: ignore[arg-type]

        self._verify_subnet_ips(ip_addresses)
        self._verify_security_group_ids(security_group_ids)
        if creator_request_id in [
            x.creator_request_id for x in self.resolver_endpoints.values()
        ]:
            raise ResourceExistsException(
                f"Resolver endpoint with creator request ID "
                f"'{creator_request_id}' already exists"
            )

        endpoint_id = f"rslvr-{'in' if direction == 'INBOUND' else 'out'}-{mock_random.get_random_hex(17)}"
        resolver_endpoint = ResolverEndpoint(
            self.account_id,
            region,
            endpoint_id,
            creator_request_id,
            security_group_ids,
            direction,
            ip_addresses,
            name,
        )

        self.resolver_endpoints[endpoint_id] = resolver_endpoint
        self.tagger.tag_resource(resolver_endpoint.arn, tags or [])
        return resolver_endpoint

    def create_resolver_rule(
        self,
        region: str,
        creator_request_id: str,
        name: str,
        rule_type: str,
        domain_name: str,
        target_ips: List[Dict[str, Any]],
        resolver_endpoint_id: str,
        tags: List[Dict[str, str]],
    ) -> ResolverRule:  # pylint: disable=too-many-arguments
        """Return description for a newly created resolver rule."""
        validate_args(
            [
                ("creatorRequestId", creator_request_id),
                ("ruleType", rule_type),
                ("domainName", domain_name),
                ("name", name),
                *[("targetIps.port", x) for x in target_ips],
                ("resolverEndpointId", resolver_endpoint_id),
            ]
        )
        errmsg = self.tagger.validate_tags(
            tags or [], limit=ResolverRule.MAX_TAGS_PER_RESOLVER_RULE
        )
        if errmsg:
            raise TagValidationException(errmsg)

        rules = [x for x in self.resolver_rules.values() if x.region == region]
        if len(rules) > ResolverRule.MAX_RULES_PER_REGION:
            # Did not verify that this is the actual error message.
            raise LimitExceededException(
                f"Account '{self.account_id}' has exceeded 'max-rules'"
            )

        # Per the AWS documentation and as seen with the AWS console, target
        # ips are only relevant when the value of Rule is FORWARD.  However,
        # boto3 ignores this condition and so shall we.

        for ip_addr in [x["Ip"] for x in target_ips]:
            try:
                # boto3 fails with an InternalServiceException if IPv6
                # addresses are used, which isn't helpful.
                if not isinstance(ip_address(ip_addr), IPv4Address):
                    raise InvalidParameterException(
                        f"Only IPv4 addresses may be used: '{ip_addr}'"
                    )
            except ValueError as exc:
                raise InvalidParameterException(
                    f"Invalid IP address: '{ip_addr}'"
                ) from exc

        # The boto3 documentation indicates that ResolverEndpoint is
        # optional, as does the AWS documention.  But if resolver_endpoint_id
        # is set to None or an empty string, it results in boto3 raising
        # a ParamValidationError either regarding the type or len of string.
        if resolver_endpoint_id:
            if resolver_endpoint_id not in [
                x.id for x in self.resolver_endpoints.values()
            ]:
                raise ResourceNotFoundException(
                    f"Resolver endpoint with ID '{resolver_endpoint_id}' does not exist."
                )

            if rule_type == "SYSTEM":
                raise InvalidRequestException(
                    "Cannot specify resolver endpoint ID and target IP "
                    "for SYSTEM type resolver rule"
                )

        if creator_request_id in [
            x.creator_request_id for x in self.resolver_rules.values()
        ]:
            raise ResourceExistsException(
                f"Resolver rule with creator request ID "
                f"'{creator_request_id}' already exists"
            )

        rule_id = f"rslvr-rr-{mock_random.get_random_hex(17)}"
        resolver_rule = ResolverRule(
            account_id=self.account_id,
            region=region,
            rule_id=rule_id,
            creator_request_id=creator_request_id,
            rule_type=rule_type,
            domain_name=domain_name,
            target_ips=target_ips,
            resolver_endpoint_id=resolver_endpoint_id,
            name=name,
        )

        self.resolver_rules[rule_id] = resolver_rule
        self.tagger.tag_resource(resolver_rule.arn, tags or [])
        return resolver_rule

    def _validate_resolver_endpoint_id(self, resolver_endpoint_id: str) -> None:
        """Raise an exception if the id is invalid or unknown."""
        validate_args([("resolverEndpointId", resolver_endpoint_id)])
        if resolver_endpoint_id not in self.resolver_endpoints:
            raise ResourceNotFoundException(
                f"Resolver endpoint with ID '{resolver_endpoint_id}' does not exist"
            )

    def delete_resolver_endpoint(self, resolver_endpoint_id: str) -> ResolverEndpoint:
        self._validate_resolver_endpoint_id(resolver_endpoint_id)

        # Can't delete an endpoint if there are rules associated with it.
        rules = [
            x.id
            for x in self.resolver_rules.values()
            if x.resolver_endpoint_id == resolver_endpoint_id
        ]
        if rules:
            raise InvalidRequestException(
                f"Cannot delete resolver endpoint unless its related resolver "
                f"rules are deleted.  The following rules still exist for "
                f"this resolver endpoint:  {','.join(rules)}"
            )

        self.tagger.delete_all_tags_for_resource(resolver_endpoint_id)
        resolver_endpoint = self.resolver_endpoints.pop(resolver_endpoint_id)
        resolver_endpoint.delete_eni()
        resolver_endpoint.status = "DELETING"
        resolver_endpoint.status_message = resolver_endpoint.status_message.replace(
            "Successfully created", "Deleting"
        )
        return resolver_endpoint

    def _validate_resolver_rule_id(self, resolver_rule_id: str) -> None:
        """Raise an exception if the id is invalid or unknown."""
        validate_args([("resolverRuleId", resolver_rule_id)])
        if resolver_rule_id not in self.resolver_rules:
            raise ResourceNotFoundException(
                f"Resolver rule with ID '{resolver_rule_id}' does not exist"
            )

    def delete_resolver_rule(self, resolver_rule_id: str) -> ResolverRule:
        self._validate_resolver_rule_id(resolver_rule_id)

        # Can't delete an rule unless VPC's are disassociated.
        associations = [
            x.id
            for x in self.resolver_rule_associations.values()
            if x.resolver_rule_id == resolver_rule_id
        ]
        if associations:
            raise ResourceInUseException(
                "Please disassociate this resolver rule from VPC first "
                "before deleting"
            )

        self.tagger.delete_all_tags_for_resource(resolver_rule_id)
        resolver_rule = self.resolver_rules.pop(resolver_rule_id)
        resolver_rule.status = "DELETING"
        resolver_rule.status_message = resolver_rule.status_message.replace(
            "Successfully created", "Deleting"
        )
        return resolver_rule

    def disassociate_resolver_rule(
        self, resolver_rule_id: str, vpc_id: str
    ) -> ResolverRuleAssociation:
        validate_args([("resolverRuleId", resolver_rule_id), ("vPCId", vpc_id)])

        # Non-existent rule or vpc ids?
        if resolver_rule_id not in self.resolver_rules:
            raise ResourceNotFoundException(
                f"Resolver rule with ID '{resolver_rule_id}' does not exist"
            )

        # Find the matching association for this rule and vpc.
        rule_association_id = None
        for association in self.resolver_rule_associations.values():
            if (
                resolver_rule_id == association.resolver_rule_id
                and vpc_id == association.vpc_id
            ):
                rule_association_id = association.id
                break
        else:
            raise ResourceNotFoundException(
                f"Resolver Rule Association between Resolver Rule "
                f"'{resolver_rule_id}' and VPC '{vpc_id}' does not exist"
            )

        rule_association = self.resolver_rule_associations.pop(rule_association_id)
        rule_association.status = "DELETING"
        rule_association.status_message = "Deleting Association"
        return rule_association

    def get_resolver_endpoint(self, resolver_endpoint_id: str) -> ResolverEndpoint:
        self._validate_resolver_endpoint_id(resolver_endpoint_id)
        return self.resolver_endpoints[resolver_endpoint_id]

    def get_resolver_rule(self, resolver_rule_id: str) -> ResolverRule:
        """Return info for specified resolver rule."""
        self._validate_resolver_rule_id(resolver_rule_id)
        return self.resolver_rules[resolver_rule_id]

    def get_resolver_rule_association(
        self, resolver_rule_association_id: str
    ) -> ResolverRuleAssociation:
        validate_args([("resolverRuleAssociationId", resolver_rule_association_id)])
        if resolver_rule_association_id not in self.resolver_rule_associations:
            raise ResourceNotFoundException(
                f"ResolverRuleAssociation '{resolver_rule_association_id}' does not Exist"
            )
        return self.resolver_rule_associations[resolver_rule_association_id]

    @paginate(pagination_model=PAGINATION_MODEL)  # type: ignore[misc]
    def list_resolver_endpoint_ip_addresses(
        self, resolver_endpoint_id: str
    ) -> List[Dict[str, Any]]:
        self._validate_resolver_endpoint_id(resolver_endpoint_id)
        endpoint = self.resolver_endpoints[resolver_endpoint_id]
        return endpoint.ip_descriptions()

    @staticmethod
    def _add_field_name_to_filter(filters: List[Dict[str, Any]]) -> None:  # type: ignore[misc]
        """Convert both styles of filter names to lowercase snake format.

        "IP_ADDRESS_COUNT" or "IpAddressCount" will become "ip_address_count".
        However, "HostVPCId" doesn't fit the pattern, so that's treated
        special.
        """
        for rr_filter in filters:
            filter_name = rr_filter["Name"]
            if "_" not in filter_name:
                if "Vpc" in filter_name:
                    filter_name = "WRONG"
                elif filter_name == "HostVPCId":
                    filter_name = "host_vpc_id"
                elif filter_name == "VPCId":
                    filter_name = "vpc_id"
                elif filter_name in ["Type", "TYPE"]:
                    filter_name = "rule_type"
                elif not filter_name.isupper():
                    filter_name = CAMEL_TO_SNAKE_PATTERN.sub("_", filter_name)
            rr_filter["Field"] = filter_name.lower()

    @staticmethod
    def _validate_filters(filters: Any, allowed_filter_names: List[str]) -> None:  # type: ignore[misc]
        """Raise exception if filter names are not as expected."""
        for rr_filter in filters:
            if rr_filter["Field"] not in allowed_filter_names:
                raise InvalidParameterException(
                    f"The filter '{rr_filter['Name']}' is invalid"
                )
            if "Values" not in rr_filter:
                raise InvalidParameterException(
                    f"No values specified for filter {rr_filter['Name']}"
                )

    @staticmethod
    def _matches_all_filters(entity: Any, filters: Any) -> bool:  # type: ignore[misc]
        """Return True if this entity has fields matching all the filters."""
        for rr_filter in filters:
            field_value = getattr(entity, rr_filter["Field"])

            if isinstance(field_value, list):
                if not set(field_value).intersection(rr_filter["Values"]):
                    return False
            elif isinstance(field_value, int):
                if str(field_value) not in rr_filter["Values"]:
                    return False
            elif field_value not in rr_filter["Values"]:
                return False

        return True

    @paginate(pagination_model=PAGINATION_MODEL)
    def list_resolver_endpoints(self, filters: Any) -> List[ResolverEndpoint]:
        if not filters:
            filters = []

        self._add_field_name_to_filter(filters)
        self._validate_filters(filters, ResolverEndpoint.FILTER_NAMES)

        endpoints = []
        for endpoint in sorted(self.resolver_endpoints.values(), key=lambda x: x.name):
            if self._matches_all_filters(endpoint, filters):
                endpoints.append(endpoint)
        return endpoints

    @paginate(pagination_model=PAGINATION_MODEL)
    def list_resolver_rules(self, filters: Any) -> List[ResolverRule]:
        if not filters:
            filters = []

        self._add_field_name_to_filter(filters)
        self._validate_filters(filters, ResolverRule.FILTER_NAMES)

        rules = []
        for rule in sorted(self.resolver_rules.values(), key=lambda x: x.name):
            if self._matches_all_filters(rule, filters):
                rules.append(rule)
        return rules

    @paginate(pagination_model=PAGINATION_MODEL)
    def list_resolver_rule_associations(
        self, filters: Any
    ) -> List[ResolverRuleAssociation]:
        if not filters:
            filters = []

        self._add_field_name_to_filter(filters)
        self._validate_filters(filters, ResolverRuleAssociation.FILTER_NAMES)

        rules = []
        for rule in sorted(
            self.resolver_rule_associations.values(), key=lambda x: x.name
        ):
            if self._matches_all_filters(rule, filters):
                rules.append(rule)
        return rules

    def _matched_arn(self, resource_arn: str) -> None:
        """Given ARN, raise exception if there is no corresponding resource."""
        for resolver_endpoint in self.resolver_endpoints.values():
            if resolver_endpoint.arn == resource_arn:
                return
        for resolver_rule in self.resolver_rules.values():
            if resolver_rule.arn == resource_arn:
                return
        raise ResourceNotFoundException(
            f"Resolver endpoint with ID '{resource_arn}' does not exist"
        )

    @paginate(pagination_model=PAGINATION_MODEL)
    def list_tags_for_resource(self, resource_arn: str) -> List[Dict[str, str]]:
        self._matched_arn(resource_arn)
        return self.tagger.list_tags_for_resource(resource_arn)["Tags"]

    def tag_resource(self, resource_arn: str, tags: List[Dict[str, str]]) -> None:
        self._matched_arn(resource_arn)
        errmsg = self.tagger.validate_tags(
            tags, limit=ResolverEndpoint.MAX_TAGS_PER_RESOLVER_ENDPOINT
        )
        if errmsg:
            raise TagValidationException(errmsg)
        self.tagger.tag_resource(resource_arn, tags)

    def untag_resource(self, resource_arn: str, tag_keys: List[str]) -> None:
        self._matched_arn(resource_arn)
        self.tagger.untag_resource_using_names(resource_arn, tag_keys)

    def update_resolver_endpoint(
        self, resolver_endpoint_id: str, name: str
    ) -> ResolverEndpoint:
        self._validate_resolver_endpoint_id(resolver_endpoint_id)
        validate_args([("name", name)])
        resolver_endpoint = self.resolver_endpoints[resolver_endpoint_id]
        resolver_endpoint.update_name(name)
        return resolver_endpoint

    def associate_resolver_endpoint_ip_address(
        self, resolver_endpoint_id: str, value: Dict[str, Any]
    ) -> ResolverEndpoint:
        self._validate_resolver_endpoint_id(resolver_endpoint_id)
        resolver_endpoint = self.resolver_endpoints[resolver_endpoint_id]

        if not value.get("Ip"):
            subnet_info = self.ec2_backend.describe_subnets(
                subnet_ids=[value.get("SubnetId")]  # type: ignore[list-item]
            )[0]
            value["Ip"] = subnet_info.get_available_subnet_ip(self)  # type: ignore[arg-type]
        self._verify_subnet_ips([value], False)

        resolver_endpoint.associate_ip_address(value)
        return resolver_endpoint

    def disassociate_resolver_endpoint_ip_address(
        self, resolver_endpoint_id: str, value: Dict[str, Any]
    ) -> ResolverEndpoint:
        self._validate_resolver_endpoint_id(resolver_endpoint_id)
        resolver_endpoint = self.resolver_endpoints[resolver_endpoint_id]

        if not (value.get("Ip") or value.get("IpId")):
            raise InvalidRequestException(
                "[RSLVR-00503] Need to specify either the IP ID or both subnet and IP address in order to remove IP address."
            )

        resolver_endpoint.disassociate_ip_address(value)
        return resolver_endpoint


route53resolver_backends = BackendDict(Route53ResolverBackend, "route53resolver")
