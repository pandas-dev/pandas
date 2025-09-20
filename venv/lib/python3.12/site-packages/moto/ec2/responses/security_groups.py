from typing import Any, Dict, List, Tuple

from moto.core.responses import ActionResult, EmptyResult

from ._base_response import EC2BaseResponse


def try_parse_int(value: Any, default: Any = None) -> Any:
    try:
        return int(value)
    except (TypeError, ValueError):
        return default


def parse_sg_attributes_from_dict(sg_attributes: Dict[str, Any]) -> Tuple[Any, ...]:
    ip_protocol = sg_attributes.get("IpProtocol", [None])[0]
    from_port = sg_attributes.get("FromPort", [None])[0]
    to_port = sg_attributes.get("ToPort", [None])[0]

    ip_ranges: List[Dict[str, Any]] = []
    ip_ranges_tree = sg_attributes.get("IpRanges") or {}
    for ip_range_idx in sorted(ip_ranges_tree.keys()):
        ip_range = {"CidrIp": ip_ranges_tree[ip_range_idx]["CidrIp"][0]}
        if ip_ranges_tree[ip_range_idx].get("Description"):
            ip_range["Description"] = ip_ranges_tree[ip_range_idx].get("Description")[0]

        ip_ranges.append(ip_range)

    ip_ranges_tree: Dict[str, Any] = sg_attributes.get("Ipv6Ranges") or {}  # type: ignore[no-redef]
    for ip_range_idx in sorted(ip_ranges_tree.keys()):
        ip_range = {"CidrIpv6": ip_ranges_tree[ip_range_idx]["CidrIpv6"][0]}
        if ip_ranges_tree[ip_range_idx].get("Description"):
            ip_range["Description"] = ip_ranges_tree[ip_range_idx].get("Description")[0]

        ip_ranges.append(ip_range)

    if "CidrIp" in sg_attributes:
        cidr_ip = sg_attributes.get("CidrIp")[0]  # type: ignore
        ip_ranges.append({"CidrIp": cidr_ip})

    if "CidrIpv6" in sg_attributes:
        cidr_ipv6 = sg_attributes.get("CidrIpv6")[0]  # type: ignore
        ip_ranges.append({"CidrIpv6": cidr_ipv6})

    source_groups: List[Dict[str, Any]] = []
    groups_tree: Dict[str, Any] = sg_attributes.get("Groups") or {}
    for group_idx in sorted(groups_tree.keys()):
        group_dict = groups_tree[group_idx]
        source_group = {}
        if "GroupId" in group_dict:
            source_group["GroupId"] = group_dict["GroupId"][0]
        if "GroupName" in group_dict:
            source_group["GroupName"] = group_dict["GroupName"][0]
        if "Description" in group_dict:
            source_group["Description"] = group_dict["Description"][0]
        if "OwnerId" in group_dict:
            source_group["OwnerId"] = group_dict["OwnerId"][0]
        source_groups.append(source_group)

    prefix_list_ids: List[Dict[str, Any]] = []
    pl_tree: Dict[str, Any] = sg_attributes.get("PrefixListIds") or {}
    for pl_index in sorted(pl_tree):
        pl_dict = pl_tree.get(pl_index, {})
        pl_item = {}
        if "PrefixListId" in pl_dict:
            pl_item["PrefixListId"] = pl_dict.get("PrefixListId")[0]
        if "Description" in pl_dict:
            pl_item["Description"] = pl_dict.get("Description")[0]
        if pl_item:
            prefix_list_ids.append(pl_item)
    return (ip_protocol, from_port, to_port, ip_ranges, source_groups, prefix_list_ids)


class SecurityGroups(EC2BaseResponse):
    def _process_rules_from_querystring(self) -> Any:
        group_name_or_id = self._get_param("GroupName") or self._get_param("GroupId")
        security_rule_ids = self._get_multi_param("SecurityGroupRuleId")

        querytree: Dict[str, Any] = {}
        for key, value in self.querystring.items():
            key_splitted = key.split(".")
            key_splitted = [try_parse_int(e, e) for e in key_splitted]

            d = querytree
            for subkey in key_splitted[:-1]:
                if subkey not in d:
                    d[subkey] = {}
                d = d[subkey]
            d[key_splitted[-1]] = value

        sg_rule_tags = self._parse_tag_specification().get("security-group-rule", {})

        if "IpPermissions" not in querytree:
            # Handle single rule syntax
            (
                ip_protocol,
                from_port,
                to_port,
                ip_ranges,
                source_groups,
                prefix_list_ids,
            ) = parse_sg_attributes_from_dict(querytree)

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

        ip_permissions = querytree.get("IpPermissions") or {}
        for ip_permission_idx in sorted(ip_permissions.keys()):
            ip_permission = ip_permissions[ip_permission_idx]

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
        description = self._get_param("GroupDescription")
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
        groupnames = self._get_multi_param("GroupName")
        group_ids = self._get_multi_param("GroupId")
        filters = self._filters_from_querystring()

        groups = self.ec2_backend.describe_security_groups(
            group_ids=group_ids, groupnames=groupnames, filters=filters
        )

        result = {"SecurityGroups": groups}
        return ActionResult(result)

    def describe_security_group_rules(self) -> ActionResult:
        self.error_on_dryrun()
        sg_rule_ids = self._get_multi_param("SecurityGroupRuleId")
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
        security_group_rules_param = self._get_params()["SecurityGroupRule"]
        for idx, sgr in security_group_rules_param.items():
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
