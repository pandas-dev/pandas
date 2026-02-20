import copy
import itertools
import json
from collections import defaultdict
from collections.abc import Iterator
from typing import Any, Optional

from moto.core.common_models import BaseModel, CloudFormationModel
from moto.core.utils import aws_api_matches

from ..exceptions import (
    InvalidCIDRSubnetError,
    InvalidGroupIdMalformedError,
    InvalidParameterCombination,
    InvalidParameterValue,
    InvalidPermissionDuplicateError,
    InvalidPermissionNotFoundError,
    InvalidSecurityGroupDuplicateError,
    InvalidSecurityGroupNotFoundError,
    InvalidSecurityGroupRuleIdNotFoundError,
    MissingParameter,
    MissingParameterError,
    MotoNotImplementedError,
    RulesPerSecurityGroupLimitExceededError,
)
from ..utils import (
    is_tag_filter,
    is_valid_cidr,
    is_valid_ipv6_cidr,
    is_valid_security_group_id,
    random_security_group_id,
    random_security_group_rule_id,
    tag_filter_matches,
)
from .core import TaggedEC2Resource


class SecurityGroupRule(TaggedEC2Resource):
    def __init__(
        self,
        ec2_backend: Any,
        ip_protocol: str,
        group_id: str,
        from_port: Optional[str],
        to_port: Optional[str],
        ip_range: Optional[dict[str, str]],
        source_group: Optional[dict[str, str]] = None,
        prefix_list_id: Optional[dict[str, str]] = None,
        is_egress: bool = True,
        tags: Optional[dict[str, str]] = None,
    ):
        self.ec2_backend = ec2_backend
        self.id = random_security_group_rule_id()
        self.ip_protocol = str(ip_protocol) if ip_protocol else None
        self.ip_range = ip_range or {}
        self.source_group = source_group or {}
        self.prefix_list_id = None
        if prefix_list_id is not None:
            self.prefix_list_id = prefix_list_id["PrefixListId"]
        self.from_port = self.to_port = None
        self.is_egress = is_egress
        self.group_id = group_id

        if self.ip_protocol and self.ip_protocol != "-1":
            self.from_port = int(from_port)  # type: ignore[arg-type]
            self.to_port = int(to_port)  # type: ignore[arg-type]

        ip_protocol_keywords = {
            "tcp": "tcp",
            "6": "tcp",
            "udp": "udp",
            "17": "udp",
            "all": "-1",
            "-1": "-1",
            "tCp": "tcp",
            "UDp": "udp",
            "ALL": "-1",
            "icMp": "icmp",
            "1": "icmp",
            "icmp": "icmp",
        }
        proto = (
            ip_protocol_keywords.get(self.ip_protocol.lower())
            if self.ip_protocol
            else None
        )
        self.ip_protocol = proto if proto else self.ip_protocol
        self.add_tags(tags or {})
        self.filters = {
            "group-id": self.filter_group_id,
            "security-group-rule-id": self.filter_id,
        }

    @property
    def description(self) -> Optional[str]:
        return self.ip_range.get("Description")

    @property
    def cidr_ipv4(self) -> Optional[str]:
        return self.ip_range.get("CidrIp", None)

    @property
    def cidr_ipv6(self) -> Optional[str]:
        return self.ip_range.get("CidrIpv6", None)

    @property
    def referenced_group_info(self) -> Optional[dict[str, str]]:
        return self.source_group if self.source_group else None

    @property
    def owner_id(self) -> str:
        return self.ec2_backend.account_id

    def __eq__(self, other: "SecurityGroupRule") -> bool:  # type: ignore[override]
        if self.ip_protocol != other.ip_protocol:
            return False
        if "CidrIp" in self.ip_range and self.ip_range.get(
            "CidrIp"
        ) != other.ip_range.get("CidrIp"):
            return False
        if "CidrIpv6" in self.ip_range and self.ip_range.get(
            "CidrIpv6"
        ) != other.ip_range.get("CidrIpv6"):
            return False
        if self.source_group != other.source_group:
            return False
        if self.prefix_list_id != other.prefix_list_id:
            return False
        if self.ip_protocol != "-1":
            if self.from_port != other.from_port:
                return False
            if self.to_port != other.to_port:
                return False

        return True

    def __deepcopy__(self, memodict: dict[Any, Any]) -> BaseModel:
        memodict = memodict or {}
        cls = self.__class__
        new = cls.__new__(cls)
        memodict[id(self)] = new
        for k, v in self.__dict__.items():
            if k == "ec2_backend":
                setattr(new, k, self.ec2_backend)
            else:
                setattr(new, k, copy.deepcopy(v, memodict))
        return new

    def filter_id(self, values: list[Any]) -> bool:
        for value in values:
            if aws_api_matches(value, self.id):
                return True
        return False

    def filter_group_id(self, values: list[Any]) -> bool:
        for value in values:
            if aws_api_matches(value, self.group_id):
                return True
        return False

    def matches_filter(self, key: str, filter_value: Any) -> Any:
        if is_tag_filter(key):
            return tag_filter_matches(self, key, filter_value)
        else:
            return self.filters[key](filter_value)

    def matches_filters(self, filters: Any) -> bool:
        for key, value in filters.items():
            if not self.matches_filter(key, value):
                return False
        return True


class GroupedSecurityRuleView:
    def __init__(
        self,
        from_port: Optional[int],
        to_port: Optional[int],
        ip_protocol: Optional[str],
    ):
        self.from_port = from_port
        self.to_port = to_port
        self.ip_protocol = ip_protocol
        self.all_ip_ranges: list[dict[str, str]] = []
        self.source_groups: list[dict[str, str]] = []
        self.prefix_list_ids: list[dict[str, str]] = []

    @property
    def user_id_group_pairs(self) -> list[dict[str, str]]:
        return self.source_groups

    @property
    def ip_ranges(self) -> list[dict[str, str]]:
        return [ip_range for ip_range in self.all_ip_ranges if ip_range.get("CidrIp")]

    @property
    def ipv6_ranges(self) -> list[dict[str, str]]:
        return [ip_range for ip_range in self.all_ip_ranges if ip_range.get("CidrIpv6")]


class SecurityGroup(TaggedEC2Resource, CloudFormationModel):
    def __init__(
        self,
        ec2_backend: Any,
        group_id: str,
        name: str,
        description: str,
        vpc_id: Optional[str] = None,
        tags: Optional[dict[str, str]] = None,
        is_default: Optional[bool] = None,
    ):
        self.ec2_backend = ec2_backend
        self.id = group_id
        self.group_id = self.id
        self.name = name
        self.group_name = self.name
        self.description = description
        self.ingress_rules: list[SecurityGroupRule] = []
        self.egress_rules: list[SecurityGroupRule] = []
        self.vpc_id: Optional[str] = vpc_id
        self.owner_id = ec2_backend.account_id
        self.add_tags(tags or {})
        self.is_default = is_default or False
        self.arn = f"arn:aws:ec2:{ec2_backend.region_name}:{ec2_backend.account_id}:security-group/{group_id}"

        # Append default IPv6 egress rule for VPCs with IPv6 support
        if vpc_id:
            vpc = self.ec2_backend.vpcs.get(vpc_id)
            if vpc:
                self.egress_rules.append(
                    SecurityGroupRule(
                        self.ec2_backend,
                        "-1",
                        self.id,
                        None,
                        None,
                        {"CidrIp": "0.0.0.0/0"},
                    )
                )
            if vpc and len(vpc.get_cidr_block_association_set(ipv6=True)) > 0:
                self.egress_rules.append(
                    SecurityGroupRule(
                        self.ec2_backend,
                        "-1",
                        self.id,
                        None,
                        None,
                        {"CidrIpv6": "::/0"},
                    )
                )

        # each filter as a simple function in a mapping

        self.filters = {
            "description": self.filter_description,
            "egress.ip-permission.cidr": self.filter_egress__ip_permission__cidr,
            "egress.ip-permission.from-port": self.filter_egress__ip_permission__from_port,
            "egress.ip-permission.group-id": self.filter_egress__ip_permission__group_id,
            "egress.ip-permission.group-name": self.filter_egress__ip_permission__group_name,
            "egress.ip-permission.ipv6-cidr": self.filter_egress__ip_permission__ipv6_cidr,
            "egress.ip-permission.prefix-list-id": self.filter_egress__ip_permission__prefix_list_id,
            "egress.ip-permission.protocol": self.filter_egress__ip_permission__protocol,
            "egress.ip-permission.to-port": self.filter_egress__ip_permission__to_port,
            "egress.ip-permission.user-id": self.filter_egress__ip_permission__user_id,
            "group-id": self.filter_group_id,
            "group-name": self.filter_group_name,
            "ip-permission.cidr": self.filter_ip_permission__cidr,
            "ip-permission.from-port": self.filter_ip_permission__from_port,
            "ip-permission.group-id": self.filter_ip_permission__group_id,
            "ip-permission.group-name": self.filter_ip_permission__group_name,
            "ip-permission.ipv6-cidr": self.filter_ip_permission__ipv6_cidr,
            "ip-permission.prefix-list-id": self.filter_ip_permission__prefix_list_id,
            "ip-permission.protocol": self.filter_ip_permission__protocol,
            "ip-permission.to-port": self.filter_ip_permission__to_port,
            "ip-permission.user-id": self.filter_ip_permission__user_id,
            "owner-id": self.filter_owner_id,
            "vpc-id": self.filter_vpc_id,
        }

    @staticmethod
    def cloudformation_name_type() -> str:
        return "GroupName"

    @staticmethod
    def cloudformation_type() -> str:
        # https://docs.aws.amazon.com/AWSCloudFormation/latest/UserGuide/aws-resource-ec2-securitygroup.html
        return "AWS::EC2::SecurityGroup"

    @classmethod
    def create_from_cloudformation_json(  # type: ignore[misc]
        cls,
        resource_name: str,
        cloudformation_json: Any,
        account_id: str,
        region_name: str,
        **kwargs: Any,
    ) -> "SecurityGroup":
        from ..models import ec2_backends

        properties = cloudformation_json["Properties"]

        ec2_backend = ec2_backends[account_id][region_name]
        vpc_id = properties.get("VpcId")
        security_group = ec2_backend.create_security_group(
            name=resource_name,
            description=properties.get("GroupDescription"),
            vpc_id=vpc_id,
        )

        for tag in properties.get("Tags", []):
            tag_key = tag["Key"]
            tag_value = tag["Value"]
            security_group.add_tag(tag_key, tag_value)

        for ingress_rule in properties.get("SecurityGroupIngress", []):
            source_group_id = ingress_rule.get("SourceSecurityGroupId")
            source_group_name = ingress_rule.get("SourceSecurityGroupName")
            source_group = {}
            if source_group_id:
                source_group["GroupId"] = source_group_id
            if source_group_name:
                source_group["GroupName"] = source_group_name

            ec2_backend.authorize_security_group_ingress(
                group_name_or_id=security_group.id,
                ip_protocol=ingress_rule["IpProtocol"],
                from_port=ingress_rule["FromPort"],
                to_port=ingress_rule["ToPort"],
                ip_ranges=ingress_rule.get("CidrIp", []),
                source_groups=[source_group] if source_group else [],
                vpc_id=vpc_id,
            )

        return security_group

    @classmethod
    def update_from_cloudformation_json(  # type: ignore[misc]
        cls,
        original_resource: Any,
        new_resource_name: str,
        cloudformation_json: Any,
        account_id: str,
        region_name: str,
    ) -> "SecurityGroup":
        cls._delete_security_group_given_vpc_id(
            original_resource.name, original_resource.vpc_id, account_id, region_name
        )
        return cls.create_from_cloudformation_json(
            new_resource_name, cloudformation_json, account_id, region_name
        )

    @classmethod
    def delete_from_cloudformation_json(  # type: ignore[misc]
        cls,
        resource_name: str,
        cloudformation_json: Any,
        account_id: str,
        region_name: str,
    ) -> None:
        properties = cloudformation_json["Properties"]
        vpc_id = properties.get("VpcId")
        cls._delete_security_group_given_vpc_id(
            resource_name, vpc_id, account_id, region_name
        )

    @classmethod
    def _delete_security_group_given_vpc_id(
        cls, resource_name: str, vpc_id: str, account_id: str, region_name: str
    ) -> None:
        from ..models import ec2_backends

        ec2_backend = ec2_backends[account_id][region_name]
        security_group = ec2_backend.get_security_group_by_name_or_id(
            resource_name, vpc_id
        )
        if security_group:
            security_group.delete(account_id, region_name)

    def delete(
        self,
        account_id: str,
        region_name: str,
    ) -> None:
        """Not exposed as part of the ELB API - used for CloudFormation."""
        self.ec2_backend.delete_security_group(group_id=self.id)

    @property
    def physical_resource_id(self) -> str:
        return self.id

    def filter_description(self, values: list[Any]) -> bool:
        for value in values:
            if aws_api_matches(value, self.description):
                return True
        return False

    def filter_egress__ip_permission__cidr(self, values: list[Any]) -> bool:
        for value in values:
            for rule in self.egress_rules:
                if aws_api_matches(value, rule.ip_range.get("CidrIp", "NONE")):
                    return True
        return False

    def filter_egress__ip_permission__from_port(self, values: list[Any]) -> bool:
        for value in values:
            for rule in self.egress_rules:
                if rule.ip_protocol != "-1" and aws_api_matches(
                    value, str(rule.from_port)
                ):
                    return True
        return False

    def filter_egress__ip_permission__group_id(self, values: list[Any]) -> bool:
        for value in values:
            for rule in self.egress_rules:
                if aws_api_matches(value, rule.source_group.get("GroupId", None)):
                    return True
        return False

    def filter_egress__ip_permission__group_name(self, values: list[Any]) -> bool:
        for value in values:
            for rule in self.egress_rules:
                if aws_api_matches(value, rule.source_group.get("GroupName", None)):
                    return True
        return False

    def filter_egress__ip_permission__ipv6_cidr(self, values: list[Any]) -> bool:
        raise MotoNotImplementedError("egress.ip-permission.ipv6-cidr filter")

    def filter_egress__ip_permission__prefix_list_id(self, values: list[Any]) -> bool:
        raise MotoNotImplementedError("egress.ip-permission.prefix-list-id filter")

    def filter_egress__ip_permission__protocol(self, values: list[Any]) -> bool:
        for value in values:
            for rule in self.egress_rules:
                if aws_api_matches(value, rule.ip_protocol):
                    return True
        return False

    def filter_egress__ip_permission__to_port(self, values: list[Any]) -> bool:
        for value in values:
            for rule in self.egress_rules:
                if aws_api_matches(value, rule.to_port):
                    return True
        return False

    def filter_egress__ip_permission__user_id(self, values: list[Any]) -> bool:
        for value in values:
            for rule in self.egress_rules:
                if aws_api_matches(value, rule.owner_id):
                    return True
        return False

    def filter_group_id(self, values: list[Any]) -> bool:
        for value in values:
            if aws_api_matches(value, self.id):
                return True
        return False

    def filter_group_name(self, values: list[Any]) -> bool:
        for value in values:
            if aws_api_matches(value, self.group_name):
                return True
        return False

    def filter_ip_permission__cidr(self, values: list[Any]) -> bool:
        for value in values:
            for rule in self.ingress_rules:
                if aws_api_matches(value, rule.ip_range.get("CidrIp", "NONE")):
                    return True
        return False

    def filter_ip_permission__from_port(self, values: list[Any]) -> bool:
        for value in values:
            for rule in self.ingress_rules:
                if aws_api_matches(value, rule.from_port):
                    return True
        return False

    def filter_ip_permission__group_id(self, values: list[Any]) -> bool:
        for value in values:
            for rule in self.ingress_rules:
                if aws_api_matches(value, rule.source_group.get("GroupId", None)):
                    return True
        return False

    def filter_ip_permission__group_name(self, values: list[Any]) -> bool:
        for value in values:
            for rule in self.ingress_rules:
                if aws_api_matches(value, rule.source_group.get("GroupName", None)):
                    return True
        return False

    def filter_ip_permission__ipv6_cidr(self, values: list[Any]) -> None:
        raise MotoNotImplementedError("ip-permission.ipv6 filter")

    def filter_ip_permission__prefix_list_id(self, values: list[Any]) -> None:
        raise MotoNotImplementedError("ip-permission.prefix-list-id filter")

    def filter_ip_permission__protocol(self, values: list[Any]) -> bool:
        for value in values:
            for rule in self.ingress_rules:
                if aws_api_matches(value, rule.ip_protocol):
                    return True
        return False

    def filter_ip_permission__to_port(self, values: list[Any]) -> bool:
        for value in values:
            for rule in self.ingress_rules:
                if aws_api_matches(value, rule.to_port):
                    return True
        return False

    def filter_ip_permission__user_id(self, values: list[Any]) -> bool:
        for value in values:
            for rule in self.ingress_rules:
                if aws_api_matches(value, rule.owner_id):
                    return True
        return False

    def filter_owner_id(self, values: list[Any]) -> bool:
        for value in values:
            if aws_api_matches(value, self.owner_id):
                return True
        return False

    def filter_vpc_id(self, values: list[Any]) -> bool:
        for value in values:
            if aws_api_matches(value, self.vpc_id):
                return True
        return False

    def matches_filter(self, key: str, filter_value: Any) -> Any:
        if is_tag_filter(key):
            tag_value = self.get_filter_value(key)
            if isinstance(filter_value, list):
                return tag_filter_matches(self, key, filter_value)
            return tag_value in filter_value
        else:
            return self.filters[key](filter_value)

    def matches_filters(self, filters: Any) -> bool:
        for key, value in filters.items():
            if not self.matches_filter(key, value):
                return False
        return True

    @classmethod
    def has_cfn_attr(cls, attr: str) -> bool:
        return attr in ["GroupId"]

    def get_cfn_attribute(self, attribute_name: str) -> str:
        from moto.cloudformation.exceptions import UnformattedGetAttTemplateException

        if attribute_name == "GroupId":
            return self.id
        raise UnformattedGetAttTemplateException()

    def get_rule(self, rule_id: str) -> Optional[SecurityGroupRule]:
        """Retrieve a security group rule by its ID."""
        for rule in list(itertools.chain(self.egress_rules, self.ingress_rules)):
            if rule.id == rule_id:
                return rule
        return None

    def add_ingress_rule(self, rule: SecurityGroupRule) -> None:
        if rule in self.ingress_rules:
            raise InvalidPermissionDuplicateError()
        self.ingress_rules.append(rule)

    def add_egress_rule(self, rule: SecurityGroupRule) -> None:
        if rule in self.egress_rules:
            raise InvalidPermissionDuplicateError()
        self.egress_rules.append(rule)

    def get_number_of_ingress_rules(self) -> int:
        return len(self.ingress_rules)

    def get_number_of_egress_rules(self) -> int:
        return len(self.egress_rules)

    @property
    def ip_permissions(self) -> list[GroupedSecurityRuleView]:
        return self._flattened_rules(copy.copy(self.ingress_rules))

    @property
    def ip_permissions_egress(self) -> list[GroupedSecurityRuleView]:
        return self._flattened_rules(copy.copy(self.egress_rules))

    def _flattened_rules(
        self, rules: list[SecurityGroupRule]
    ) -> list[GroupedSecurityRuleView]:
        rules_to_return: list[GroupedSecurityRuleView] = []
        for rule in rules:
            for already_added in rules_to_return:
                if (
                    already_added.from_port == rule.from_port
                    and already_added.to_port == rule.to_port
                    and already_added.ip_protocol == rule.ip_protocol
                ):
                    if rule.ip_range:
                        already_added.all_ip_ranges.append(rule.ip_range)
                    if rule.source_group:
                        already_added.source_groups.append(rule.source_group)
                    if rule.prefix_list_id:
                        already_added.prefix_list_ids.append(
                            {"PrefixListId": rule.prefix_list_id}
                        )
                    break
            else:
                view = GroupedSecurityRuleView(
                    rule.from_port, rule.to_port, rule.ip_protocol
                )
                if rule.ip_range:
                    view.all_ip_ranges.append(rule.ip_range)
                if rule.source_group:
                    view.source_groups.append(rule.source_group)
                if rule.prefix_list_id:
                    view.prefix_list_ids.append({"PrefixListId": rule.prefix_list_id})
                rules_to_return.append(view)
        return rules_to_return


class SecurityGroupBackend:
    def __init__(self) -> None:
        # the key in the dict group is the vpc_id or None (non-vpc)
        self.groups: dict[str, dict[str, SecurityGroup]] = defaultdict(dict)

    def create_security_group(
        self,
        name: str,
        description: str,
        vpc_id: Optional[str] = None,
        tags: Optional[dict[str, str]] = None,
        force: bool = False,
        is_default: Optional[bool] = None,
    ) -> SecurityGroup:
        vpc_id = vpc_id or self.default_vpc.id  # type: ignore[attr-defined]
        if not description:
            raise MissingParameterError("GroupDescription")

        group_id = random_security_group_id()
        if not force:
            existing_group = self.get_security_group_by_name_or_id(name, vpc_id)
            if existing_group:
                raise InvalidSecurityGroupDuplicateError(name)
        group = SecurityGroup(
            self,
            group_id,
            name,
            description,
            vpc_id=vpc_id,
            tags=tags,
            is_default=is_default,
        )

        self.groups[vpc_id][group_id] = group
        return group

    def describe_security_groups(
        self,
        group_ids: Optional[list[str]] = None,
        groupnames: Optional[list[str]] = None,
        filters: Any = None,
    ) -> list[SecurityGroup]:
        all_groups = self.groups.copy()
        matches = list(
            itertools.chain(*[x.copy().values() for x in all_groups.values()])
        )
        if group_ids:
            matches = [grp for grp in matches if grp.id in group_ids]
            if len(group_ids) > len(matches):
                unknown_ids = set(group_ids) - set(matches)  # type: ignore[arg-type]
                raise InvalidSecurityGroupNotFoundError(unknown_ids)
        if groupnames:
            matches = [grp for grp in matches if grp.name in groupnames]
            if len(groupnames) > len(matches):
                unknown_names = set(groupnames) - set(matches)  # type: ignore[arg-type]
                raise InvalidSecurityGroupNotFoundError(unknown_names)
        if filters:
            matches = [grp for grp in matches if grp.matches_filters(filters)]

        return matches

    def describe_security_group_rules(
        self,
        sg_rule_ids: list[str],
        filters: dict[str, list[str]],
    ) -> list[SecurityGroupRule]:
        if "group-id" in filters:
            for group_id in filters["group-id"]:
                if not is_valid_security_group_id(group_id):
                    raise InvalidGroupIdMalformedError(group_id)
        rules = [
            copy.copy(rule)
            for vpc_and_group in self.groups.values()
            for group in vpc_and_group.values()
            for rule in list(itertools.chain(group.egress_rules, group.ingress_rules))
        ]
        # Rules that have no from_port or to_port are set to -1 to match AWS behavior for this API call.
        for rule in rules:
            if rule.from_port is None:
                rule.from_port = -1
            if rule.to_port is None:
                rule.to_port = -1

        if sg_rule_ids:
            # If sg_rule_ids is provided, we convert it into a filter.
            if "security-group-rule-id" not in filters:
                filters["security-group-rule-id"] = []
            rule_ids = [rule.id for rule in rules]
            for rule_id in sg_rule_ids:
                if rule_id not in rule_ids:
                    raise InvalidSecurityGroupRuleIdNotFoundError(rule_id)
                filters["security-group-rule-id"].append(rule_id)
        if filters:
            rules = [rule for rule in rules if rule.matches_filters(filters)]
        return rules

    def _delete_security_group(self, vpc_id: Optional[str], group_id: str) -> None:
        vpc_id = vpc_id or self.default_vpc.id  # type: ignore[attr-defined]
        self.groups[vpc_id].pop(group_id)

    def delete_security_group(
        self, name: Optional[str] = None, group_id: Optional[str] = None
    ) -> None:
        if group_id:
            # loop over all the SGs, find the right one
            for vpc_id, groups in self.groups.items():
                if group_id in groups:
                    return self._delete_security_group(vpc_id, group_id)
            raise InvalidSecurityGroupNotFoundError(group_id)
        elif name:
            # Group Name.  Has to be in standard EC2, VPC needs to be
            # identified by group_id
            group = self.get_security_group_by_name_or_id(name)
            if group:
                return self._delete_security_group(None, group.id)
            raise InvalidSecurityGroupNotFoundError(name)

    def get_security_group_from_id(self, group_id: str) -> Optional[SecurityGroup]:
        # 2 levels of chaining necessary since it's a complex structure
        all_groups = itertools.chain.from_iterable(
            [x.copy().values() for x in self.groups.copy().values()]
        )
        for group in all_groups:
            if group.id == group_id:
                return group
        return None

    def get_security_group_from_name(
        self, name: str, vpc_id: Optional[str] = None
    ) -> Optional[SecurityGroup]:
        if vpc_id:
            for group in self.groups[vpc_id].values():
                if group.name == name:
                    return group
        else:
            for vpc_id in self.groups:
                for group in self.groups[vpc_id].values():
                    if group.name == name:
                        return group
        return None

    def get_security_group_by_name_or_id(
        self, group_name_or_id: str, vpc_id: Optional[str] = None
    ) -> Optional[SecurityGroup]:
        # try searching by id, fallbacks to name search
        group = self.get_security_group_from_id(group_name_or_id)
        if group is None:
            group = self.get_security_group_from_name(group_name_or_id, vpc_id)
        return group

    def get_default_security_group(
        self, vpc_id: Optional[str] = None
    ) -> Optional[SecurityGroup]:
        for group in self.groups[vpc_id or self.default_vpc.id].values():  # type: ignore[attr-defined]
            if group.is_default:
                return group
        return None

    def _iterate_security_rules(
        self,
        ip_protocol: str,
        group_id: str,
        from_port: str,
        to_port: str,
        ip_ranges: list[Any],
        source_groups: list[dict[str, Any]],
        prefix_list_ids: list[dict[str, str]],
        is_egress: bool = False,
        tags: Optional[dict[str, str]] = None,
    ) -> Iterator[SecurityGroupRule]:
        for ip_range in ip_ranges:
            yield SecurityGroupRule(
                self,
                ip_protocol,
                group_id,
                from_port,
                to_port,
                ip_range,
                None,
                None,
                is_egress=is_egress,
                tags=tags or {},
            )

        for source_group in source_groups:
            yield SecurityGroupRule(
                self,
                ip_protocol,
                group_id,
                from_port,
                to_port,
                None,
                source_group,
                None,
                is_egress=is_egress,
                tags=tags or {},
            )

        for prefix_list_id in prefix_list_ids:
            yield SecurityGroupRule(
                self,
                ip_protocol,
                group_id,
                from_port,
                to_port,
                None,
                None,
                prefix_list_id,
                is_egress=is_egress,
                tags=tags or {},
            )

    def modify_security_group_rules(
        self,
        group_id: str,
        rules: dict[str, dict[str, Any]],
    ) -> None:
        group = self.get_security_group_by_name_or_id(group_id)
        if group is None:
            raise InvalidSecurityGroupNotFoundError(group_id)

        for rule_id, new_rule in rules.items():
            cidr_ipv4 = new_rule.get("CidrIpv4")
            cidr_ipv6 = new_rule.get("CidrIpv6")
            prefix_list_id = new_rule.get("PrefixListId")
            reference_group_id = new_rule.get("ReferencedGroupId")

            set_param_count = sum(
                item is not None
                for item in [cidr_ipv4, cidr_ipv6, prefix_list_id, reference_group_id]
            )
            if set_param_count > 1:
                raise InvalidParameterCombination(
                    "Only one of cidrIp, cidrIpv6, prefixListId, or referencedGroupId can be specified"
                )
            if set_param_count < 1:
                raise MissingParameter(
                    "The request must contain exactly one of: cidrIp, cidrIpv6, prefixListId, or referencedGroupId"
                )

            existing_rule = group.get_rule(rule_id)
            if existing_rule is None:
                raise InvalidSecurityGroupRuleIdNotFoundError(rule_id)

            ip_protocol = new_rule.get("IpProtocol")
            if ip_protocol is None:
                raise InvalidParameterValue(
                    "Invalid value 'null' for protocol. VPC security group rules must specify protocols explicitly."
                )

            from_port = new_rule.get("FromPort")
            to_port = new_rule.get("ToPort")
            if from_port is None or to_port is None:
                raise InvalidParameterValue(
                    "Invalid value for portRange. Must specify both from and to ports with TCP/UDP."
                )

            if cidr_ipv4 and "CidrIpv6" in existing_rule.ip_range:
                raise InvalidParameterValue(
                    f"Invalid rule type for security group rule '{rule_id}'. You may not specify CidrIpv6 for an existing IPv4 CIDR rule."
                )

            existing_rule.ip_protocol = ip_protocol
            existing_rule.from_port = from_port
            existing_rule.to_port = to_port

            # TODO: Handle cidr_ipv6, prefix_list_id and reference_group_id

            description = new_rule.get("Description")
            if description is not None:
                existing_rule.ip_range["Description"] = description

            if cidr_ipv4 is not None:
                existing_rule.ip_range["CidrIp"] = cidr_ipv4

    def authorize_security_group_ingress(
        self,
        group_name_or_id: str,
        ip_protocol: str,
        from_port: str,
        to_port: str,
        ip_ranges: list[Any],
        sgrule_tags: Optional[dict[str, str]] = None,
        source_groups: Optional[list[dict[str, str]]] = None,
        prefix_list_ids: Optional[list[dict[str, str]]] = None,
        security_rule_ids: Optional[list[str]] = None,
        vpc_id: Optional[str] = None,
    ) -> tuple[list[SecurityGroupRule], SecurityGroup]:
        group = self.get_security_group_by_name_or_id(group_name_or_id, vpc_id)
        if group is None:
            raise InvalidSecurityGroupNotFoundError(group_name_or_id)
        if ip_ranges:
            if isinstance(ip_ranges, str):
                ip_ranges = [{"CidrIp": str(ip_ranges)}]
            elif not isinstance(ip_ranges, list):
                ip_ranges = [json.loads(ip_ranges)]
        if ip_ranges:
            for cidr in ip_ranges:
                if (
                    isinstance(cidr, dict)
                    and not any(
                        [
                            is_valid_cidr(cidr.get("CidrIp", "")),
                            is_valid_ipv6_cidr(cidr.get("CidrIpv6", "")),
                        ]
                    )
                ) or (
                    isinstance(cidr, str)
                    and not any([is_valid_cidr(cidr), is_valid_ipv6_cidr(cidr)])
                ):
                    raise InvalidCIDRSubnetError(cidr=cidr)

        self._verify_group_will_respect_rule_count_limit(
            group, group.get_number_of_ingress_rules(), ip_ranges, source_groups
        )

        _source_groups = self._add_source_group(source_groups, vpc_id)

        rules_added: list[SecurityGroupRule] = []

        for security_rule in self._iterate_security_rules(
            ip_protocol,
            group.group_id,
            from_port,
            to_port,
            ip_ranges,
            _source_groups,
            prefix_list_ids or [],
            is_egress=False,
            tags=sgrule_tags or {},
        ):
            if security_rule in group.ingress_rules:
                raise InvalidPermissionDuplicateError()
            group.add_ingress_rule(security_rule)
            rules_added.append(security_rule)

        return rules_added, group

    def revoke_security_group_ingress(
        self,
        group_name_or_id: str,
        ip_protocol: str,
        from_port: str,
        to_port: str,
        ip_ranges: list[Any],
        source_groups: Optional[list[dict[str, Any]]] = None,
        prefix_list_ids: Optional[list[dict[str, str]]] = None,
        security_rule_ids: Optional[list[str]] = None,
        vpc_id: Optional[str] = None,
    ) -> None:
        group: SecurityGroup = self.get_security_group_by_name_or_id(
            group_name_or_id, vpc_id
        )  # type: ignore[assignment]

        if group is None:
            raise InvalidSecurityGroupNotFoundError(group_name_or_id)

        rules_to_remove: list[str] = []
        has_unknown_rules = False

        if security_rule_ids:
            ingress_rule_ids = [rule.id for rule in group.ingress_rules]
            for rule_id in security_rule_ids:
                if rule_id in ingress_rule_ids:
                    rules_to_remove.append(rule_id)
                else:
                    has_unknown_rules = True
                    break
        else:
            _source_groups = self._add_source_group(source_groups, vpc_id)

            for security_rule in self._iterate_security_rules(
                ip_protocol,
                group.group_id,
                from_port,
                to_port,
                ip_ranges,
                _source_groups,
                prefix_list_ids or [],
                is_egress=False,
            ):
                try:
                    idx = group.ingress_rules.index(security_rule)
                    rules_to_remove.append(group.ingress_rules[idx].id)
                except ValueError:
                    has_unknown_rules = True
                    break

        if has_unknown_rules:
            raise InvalidPermissionNotFoundError()

        group.ingress_rules = [
            rule for rule in group.ingress_rules if rule.id not in rules_to_remove
        ]

    def authorize_security_group_egress(
        self,
        group_name_or_id: str,
        ip_protocol: str,
        from_port: str,
        to_port: str,
        ip_ranges: list[Any],
        sgrule_tags: Optional[dict[str, str]] = None,
        source_groups: Optional[list[dict[str, Any]]] = None,
        prefix_list_ids: Optional[list[dict[str, str]]] = None,
        security_rule_ids: Optional[list[str]] = None,
        vpc_id: Optional[str] = None,
    ) -> tuple[list[SecurityGroupRule], SecurityGroup]:
        group = self.get_security_group_by_name_or_id(group_name_or_id, vpc_id)
        if group is None:
            raise InvalidSecurityGroupNotFoundError(group_name_or_id)
        if ip_ranges and not isinstance(ip_ranges, list):
            if isinstance(ip_ranges, str) and "CidrIp" not in ip_ranges:
                ip_ranges = [{"CidrIp": ip_ranges}]
            else:
                ip_ranges = [json.loads(ip_ranges)]
        if ip_ranges:
            for cidr in ip_ranges:
                if (
                    isinstance(cidr, dict)
                    and not any(
                        [
                            is_valid_cidr(cidr.get("CidrIp", "")),
                            is_valid_ipv6_cidr(cidr.get("CidrIpv6", "")),
                        ]
                    )
                ) or (
                    isinstance(cidr, str)
                    and not any([is_valid_cidr(cidr), is_valid_ipv6_cidr(cidr)])
                ):
                    raise InvalidCIDRSubnetError(cidr=cidr)
        self._verify_group_will_respect_rule_count_limit(
            group,
            group.get_number_of_egress_rules(),
            ip_ranges,
            source_groups,
            egress=True,
        )

        _source_groups = self._add_source_group(source_groups, vpc_id)

        rules_added: list[SecurityGroupRule] = []

        for security_rule in self._iterate_security_rules(
            ip_protocol,
            group.group_id,
            from_port,
            to_port,
            ip_ranges,
            _source_groups,
            prefix_list_ids or [],
            is_egress=True,
            tags=sgrule_tags or {},
        ):
            if security_rule in group.egress_rules:
                raise InvalidPermissionDuplicateError()
            group.add_egress_rule(security_rule)
            rules_added.append(security_rule)

        return rules_added, group

    def revoke_security_group_egress(
        self,
        group_name_or_id: str,
        ip_protocol: str,
        from_port: str,
        to_port: str,
        ip_ranges: list[Any],
        source_groups: Optional[list[dict[str, Any]]] = None,
        prefix_list_ids: Optional[list[dict[str, str]]] = None,
        security_rule_ids: Optional[list[str]] = None,
        vpc_id: Optional[str] = None,
    ) -> None:
        group: SecurityGroup = self.get_security_group_by_name_or_id(
            group_name_or_id, vpc_id
        )  # type: ignore[assignment]

        if group is None:
            raise InvalidSecurityGroupNotFoundError(group_name_or_id)

        rules_to_remove: list[str] = []
        has_unknown_rules = False

        if security_rule_ids:
            egress_rule_ids = [rule.id for rule in group.egress_rules]
            for rule_id in security_rule_ids:
                if rule_id in egress_rule_ids:
                    rules_to_remove.append(rule_id)
                else:
                    has_unknown_rules = True
                    break
        else:
            _source_groups = self._add_source_group(source_groups, vpc_id)

            # I don't believe this is required after changing the default egress rule
            # to be {'CidrIp': '0.0.0.0/0'} instead of just '0.0.0.0/0'
            # Not sure why this would return only the IP if it was 0.0.0.0/0 instead of
            # the ip_range?
            # for ip in ip_ranges:
            #     ip_ranges = [ip.get("CidrIp") if ip.get("CidrIp") == "0.0.0.0/0" else ip]

            if group.vpc_id:
                vpc = self.vpcs.get(group.vpc_id)  # type: ignore[attr-defined]
                if vpc and not len(vpc.get_cidr_block_association_set(ipv6=True)) > 0:
                    for item in ip_ranges.copy():
                        if "CidrIpv6" in item:
                            ip_ranges.remove(item)

            for security_rule in self._iterate_security_rules(
                ip_protocol,
                group.group_id,
                from_port,
                to_port,
                ip_ranges,
                _source_groups,
                prefix_list_ids or [],
                is_egress=True,
            ):
                try:
                    idx = group.egress_rules.index(security_rule)
                    rules_to_remove.append(group.egress_rules[idx].id)
                except ValueError:
                    has_unknown_rules = True
                    break

        if has_unknown_rules:
            raise InvalidPermissionNotFoundError()

        group.egress_rules = [
            rule for rule in group.egress_rules if rule.id not in rules_to_remove
        ]

    def update_security_group_rule_descriptions_ingress(
        self,
        group_name_or_id: str,
        ip_protocol: str,
        from_port: str,
        to_port: str,
        ip_ranges: list[str],
        source_groups: Optional[list[dict[str, Any]]] = None,
        prefix_list_ids: Optional[list[dict[str, str]]] = None,
        security_rule_ids: Optional[list[str]] = None,
        vpc_id: Optional[str] = None,
    ) -> SecurityGroup:
        group = self.get_security_group_by_name_or_id(group_name_or_id, vpc_id)
        if group is None:
            raise InvalidSecurityGroupNotFoundError(group_name_or_id)
        if ip_ranges and not isinstance(ip_ranges, list):
            if isinstance(ip_ranges, str) and "CidrIp" not in ip_ranges:
                ip_ranges = [{"CidrIp": ip_ranges}]
            else:
                ip_ranges = [json.loads(ip_ranges)]
        if ip_ranges:
            for cidr in ip_ranges:
                if (
                    isinstance(cidr, dict)
                    and not any(
                        [
                            is_valid_cidr(cidr.get("CidrIp", "")),
                            is_valid_ipv6_cidr(cidr.get("CidrIpv6", "")),
                        ]
                    )
                ) or (
                    isinstance(cidr, str)
                    and not any([is_valid_cidr(cidr), is_valid_ipv6_cidr(cidr)])
                ):
                    raise InvalidCIDRSubnetError(cidr=cidr)
        _source_groups = self._add_source_group(source_groups, vpc_id)

        for security_rule in self._iterate_security_rules(
            ip_protocol,
            group.group_id,
            from_port,
            to_port,
            ip_ranges,
            _source_groups,
            prefix_list_ids or [],
            is_egress=False,
        ):
            try:
                idx = group.ingress_rules.index(security_rule)
                self._sg_update_description(security_rule, group.ingress_rules[idx])
            except ValueError:
                continue

        return group

    def update_security_group_rule_descriptions_egress(
        self,
        group_name_or_id: str,
        ip_protocol: str,
        from_port: str,
        to_port: str,
        ip_ranges: list[str],
        source_groups: Optional[list[dict[str, Any]]] = None,
        prefix_list_ids: Optional[list[dict[str, str]]] = None,
        security_rule_ids: Optional[list[str]] = None,
        vpc_id: Optional[str] = None,
    ) -> SecurityGroup:
        group = self.get_security_group_by_name_or_id(group_name_or_id, vpc_id)
        if group is None:
            raise InvalidSecurityGroupNotFoundError(group_name_or_id)
        if ip_ranges and not isinstance(ip_ranges, list):
            if isinstance(ip_ranges, str) and "CidrIp" not in ip_ranges:
                ip_ranges = [{"CidrIp": ip_ranges}]
            else:
                ip_ranges = [json.loads(ip_ranges)]
        if ip_ranges:
            for cidr in ip_ranges:
                if (
                    isinstance(cidr, dict)
                    and not any(
                        [
                            is_valid_cidr(cidr.get("CidrIp", "")),
                            is_valid_ipv6_cidr(cidr.get("CidrIpv6", "")),
                        ]
                    )
                ) or (
                    isinstance(cidr, str)
                    and not any([is_valid_cidr(cidr), is_valid_ipv6_cidr(cidr)])
                ):
                    raise InvalidCIDRSubnetError(cidr=cidr)
        _source_groups = self._add_source_group(source_groups, vpc_id)

        for security_rule in self._iterate_security_rules(
            ip_protocol,
            group.group_id,
            from_port,
            to_port,
            ip_ranges,
            _source_groups,
            prefix_list_ids or [],
            is_egress=True,
        ):
            try:
                idx = group.egress_rules.index(security_rule)
                self._sg_update_description(security_rule, group.egress_rules[idx])
            except ValueError:
                continue

        return group

    def _sg_update_description(
        self, security_rule: SecurityGroupRule, rule: SecurityGroupRule
    ) -> None:
        if "Description" in security_rule.ip_range:
            description = security_rule.ip_range["Description"]
            if "CidrIp" in rule.ip_range and rule.ip_range.get(
                "CidrIp"
            ) == security_rule.ip_range.get("CidrIp"):
                rule.ip_range["Description"] = description
            elif "CidrIpv6" in rule.ip_range and rule.ip_range.get(
                "CidrIpv6"
            ) == security_rule.ip_range.get("CidrIpv6"):
                rule.ip_range["Description"] = description

        if "Description" in security_rule.source_group:
            description = security_rule.source_group["Description"]
            if security_rule.source_group.get("GroupId") == rule.source_group.get(
                "GroupId"
            ) or security_rule.source_group.get("GroupName") == rule.source_group.get(
                "GroupName"
            ):
                rule.source_group["Description"] = description

    def _add_source_group(
        self, source_groups: Optional[list[dict[str, Any]]], vpc_id: Optional[str]
    ) -> list[dict[str, Any]]:
        _source_groups = []
        for item in source_groups or []:
            if "OwnerId" not in item:
                item["OwnerId"] = self.account_id  # type: ignore[attr-defined]
            if "UserId" not in item:
                item["UserId"] = self.account_id  # type: ignore[attr-defined]
            # for VPCs
            if "GroupId" in item:
                if not self.get_security_group_by_name_or_id(item["GroupId"], vpc_id):
                    raise InvalidSecurityGroupNotFoundError(item["GroupId"])
            if "GroupName" in item:
                source_group = self.get_security_group_by_name_or_id(
                    item["GroupName"], vpc_id
                )
                if not source_group:
                    raise InvalidSecurityGroupNotFoundError(item["GroupName"])
                else:
                    item["GroupId"] = source_group.id
                    item.pop("GroupName")

            _source_groups.append(item)
        return _source_groups

    def _verify_group_will_respect_rule_count_limit(
        self,
        group: SecurityGroup,
        current_rule_nb: int,
        ip_ranges: list[str],
        source_groups: Optional[list[dict[str, str]]] = None,
        egress: bool = False,
    ) -> None:
        max_nb_rules = 60 if group.vpc_id else 100
        future_group_nb_rules = current_rule_nb
        if ip_ranges:
            future_group_nb_rules += len(ip_ranges)
        if source_groups:
            future_group_nb_rules += len(source_groups)
        if future_group_nb_rules > max_nb_rules:
            raise RulesPerSecurityGroupLimitExceededError


class SecurityGroupIngress(CloudFormationModel):
    def __init__(self, security_group: SecurityGroup, properties: Any):
        self.security_group = security_group
        self.properties = properties

    @staticmethod
    def cloudformation_name_type() -> str:
        return ""

    @staticmethod
    def cloudformation_type() -> str:
        # https://docs.aws.amazon.com/AWSCloudFormation/latest/UserGuide/aws-resource-ec2-securitygroupingress.html
        return "AWS::EC2::SecurityGroupIngress"

    @classmethod
    def create_from_cloudformation_json(  # type: ignore[misc]
        cls,
        resource_name: str,
        cloudformation_json: Any,
        account_id: str,
        region_name: str,
        **kwargs: Any,
    ) -> "SecurityGroupIngress":
        from ..models import ec2_backends

        properties = cloudformation_json["Properties"]

        ec2_backend = ec2_backends[account_id][region_name]
        group_name = properties.get("GroupName")
        group_id = properties.get("GroupId")
        ip_protocol = properties.get("IpProtocol")
        cidr_ip = properties.get("CidrIp")
        cidr_desc = properties.get("Description")
        cidr_ipv6 = properties.get("CidrIpv6")
        from_port = properties.get("FromPort")
        source_security_group_id = properties.get("SourceSecurityGroupId")
        source_security_group_name = properties.get("SourceSecurityGroupName")
        # source_security_owner_id =
        # properties.get("SourceSecurityGroupOwnerId")  # IGNORED AT THE MOMENT
        to_port = properties.get("ToPort")

        assert group_id or group_name
        assert (
            source_security_group_name
            or cidr_ip
            or cidr_ipv6
            or source_security_group_id
        )
        assert ip_protocol

        source_group = {}
        if source_security_group_id:
            source_group["GroupId"] = source_security_group_id
        if source_security_group_name:
            source_group["GroupName"] = source_security_group_name
        if cidr_ip:
            ip_ranges = [{"CidrIp": cidr_ip, "Description": cidr_desc}]
        else:
            ip_ranges = []

        if group_id:
            security_group = ec2_backend.describe_security_groups(group_ids=[group_id])[
                0
            ]
        else:
            security_group = ec2_backend.describe_security_groups(
                groupnames=[group_name]
            )[0]

        ec2_backend.authorize_security_group_ingress(
            group_name_or_id=security_group.id,
            ip_protocol=ip_protocol,
            from_port=from_port,
            to_port=to_port,
            ip_ranges=ip_ranges,
            source_groups=[source_group] if source_group else [],
        )

        return cls(security_group, properties)
