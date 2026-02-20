from typing import Any

from moto.core.responses import ActionResult, EmptyResult

from ._base_response import EC2BaseResponse


def parse_sg_attributes_from_dict(sg_attributes: dict[str, Any]) -> tuple[Any, ...]:
    ip_protocol = sg_attributes.get("IpProtocol")
    from_port = sg_attributes.get("FromPort")
    to_port = sg_attributes.get("ToPort")

    ip_ranges: list[dict[str, Any]] = []
    ip_ranges += sg_attributes.get("IpRanges", [])
    ip_ranges += sg_attributes.get("Ipv6Ranges", [])  # type: ignore[no-redef]

    if "CidrIp" in sg_attributes:
        cidr_ip = sg_attributes.get("CidrIp")  # type: ignore
        ip_ranges.append({"CidrIp": cidr_ip})

    if "CidrIpv6" in sg_attributes:
        cidr_ipv6 = sg_attributes.get("CidrIpv6")  # type: ignore
        ip_ranges.append({"CidrIpv6": cidr_ipv6})

    source_groups: list[dict[str, Any]] = sg_attributes.get("UserIdGroupPairs", [])
    prefix_list_ids: list[dict[str, Any]] = sg_attributes.get("PrefixListIds", [])

    return ip_protocol, from_port, to_port, ip_ranges, source_groups, prefix_list_ids


class SecurityGroups(EC2BaseResponse):
    def _process_rules_from_querystring(self) -> Any:
        group_name_or_id = self._get_param("GroupName") or self._get_param("GroupId")
        security_rule_ids = self._get_param("SecurityGroupRuleIds", [])
        sg_rule_tags = self._parse_tag_specification().get("security-group-rule", {})

        if "IpPermissions" not in self._get_params():
            # Handle single rule syntax
            (
                ip_protocol,
                from_port,
                to_port,
                ip_ranges,
                source_groups,
                prefix_list_ids,
            ) = parse_sg_attributes_from_dict(self._get_params())

            yield {
                "group_name_or_id": group_name_or_id,
                "ip_protocol": ip_protocol,
                "from_port": from_port,
                "to_port": to_port,
                "ip_ranges": ip_ranges,
                "sgrule_tags": sg_rule_tags,
                "source_groups": source_groups,
                "prefix_list_ids": prefix_list_ids,
                "security_rule_ids": security_rule_ids,
            }

        ip_permissions = self._get_param("IpPermissions", [])
        for ip_permission in ip_permissions:
            (
                ip_protocol,
                from_port,
                to_port,
                ip_ranges,
                source_groups,
                prefix_list_ids,
            ) = parse_sg_attributes_from_dict(ip_permission)

            yield {
                "group_name_or_id": group_name_or_id,
                "ip_protocol": ip_protocol,
                "from_port": from_port,
                "to_port": to_port,
                "ip_ranges": ip_ranges,
                "sgrule_tags": sg_rule_tags,
                "source_groups": source_groups,
                "prefix_list_ids": prefix_list_ids,
                "security_rule_ids": security_rule_ids,
            }

    def authorize_security_group_egress(self) -> ActionResult:
        self.error_on_dryrun()

        for kwargs in self._process_rules_from_querystring():
            rules, _ = self.ec2_backend.authorize_security_group_egress(**kwargs)
        result = {"Return": True, "SecurityGroupRules": rules}
        return ActionResult(result)

    def authorize_security_group_ingress(self) -> ActionResult:
        self.error_on_dryrun()

        for kwargs in self._process_rules_from_querystring():
            rules, _ = self.ec2_backend.authorize_security_group_ingress(**kwargs)
        result = {"Return": True, "SecurityGroupRules": rules}
        return ActionResult(result)

    def create_security_group(self) -> ActionResult:
        name = self._get_param("GroupName")
        description = self._get_param("Description")
        vpc_id = self._get_param("VpcId")
        tags = self._parse_tag_specification().get("security-group", {})

        self.error_on_dryrun()

        group = self.ec2_backend.create_security_group(
            name, description, vpc_id=vpc_id, tags=tags
        )
        result = {
            "GroupId": group.id,
            "SecurityGroupArn": group.arn,
        }
        return ActionResult(result)

    def delete_security_group(self) -> ActionResult:
        # TODO this should raise an error if there are instances in the group.
        # See
        # http://docs.aws.amazon.com/AWSEC2/latest/APIReference/ApiReference-query-DeleteSecurityGroup.html

        name = self._get_param("GroupName")
        sg_id = self._get_param("GroupId")

        self.error_on_dryrun()

        if name:
            self.ec2_backend.delete_security_group(name)
        elif sg_id:
            self.ec2_backend.delete_security_group(group_id=sg_id)

        return EmptyResult()

    def describe_security_groups(self) -> ActionResult:
        groupnames = self._get_param("GroupNames", [])
        group_ids = self._get_param("GroupIds", [])
        filters = self._filters_from_querystring()

        groups = self.ec2_backend.describe_security_groups(
            group_ids=group_ids, groupnames=groupnames, filters=filters
        )

        result = {"SecurityGroups": groups}
        return ActionResult(result)

    def describe_security_group_rules(self) -> ActionResult:
        self.error_on_dryrun()
        sg_rule_ids = self._get_param("SecurityGroupRuleIds", [])
        filters = self._filters_from_querystring()
        rules = self.ec2_backend.describe_security_group_rules(sg_rule_ids, filters)
        result = {"SecurityGroupRules": rules}
        return ActionResult(result)

    def revoke_security_group_egress(self) -> ActionResult:
        self.error_on_dryrun()

        for args in self._process_rules_from_querystring():
            # we don't need this parameter to revoke
            del args["sgrule_tags"]
            self.ec2_backend.revoke_security_group_egress(**args)
        return ActionResult({"Return": True})

    def revoke_security_group_ingress(self) -> ActionResult:
        self.error_on_dryrun()

        for args in self._process_rules_from_querystring():
            # we don't need this parameter to revoke
            del args["sgrule_tags"]
            self.ec2_backend.revoke_security_group_ingress(**args)
        return ActionResult({"Return": True})

    def modify_security_group_rules(self) -> ActionResult:
        self.error_on_dryrun()

        rules = {}
        security_group_rules_param = self._get_param("SecurityGroupRules", [])
        for sgr in security_group_rules_param:
            rules[sgr["SecurityGroupRuleId"]] = sgr["SecurityGroupRule"]

        group_id = self._get_param("GroupId")
        self.ec2_backend.modify_security_group_rules(group_id, rules)

        return ActionResult({"Return": True})

    def update_security_group_rule_descriptions_ingress(self) -> ActionResult:
        for args in self._process_rules_from_querystring():
            # we don't need this parameter to revoke
            del args["sgrule_tags"]
            self.ec2_backend.update_security_group_rule_descriptions_ingress(**args)
        return EmptyResult()

    def update_security_group_rule_descriptions_egress(self) -> ActionResult:
        for args in self._process_rules_from_querystring():
            # we don't need this parameter to revoke
            del args["sgrule_tags"]
            self.ec2_backend.update_security_group_rule_descriptions_egress(**args)
        return EmptyResult()
