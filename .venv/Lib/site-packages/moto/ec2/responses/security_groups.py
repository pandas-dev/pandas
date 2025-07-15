from typing import Any, Dict, List, Tuple

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

    def authorize_security_group_egress(self) -> str:
        self.error_on_dryrun()

        for kwargs in self._process_rules_from_querystring():
            rules, group = self.ec2_backend.authorize_security_group_egress(**kwargs)
        template = self.response_template(AUTHORIZE_SECURITY_GROUP_EGRESS_RESPONSE)
        return template.render(rules=rules, group=group)

    def authorize_security_group_ingress(self) -> str:
        self.error_on_dryrun()

        for kwargs in self._process_rules_from_querystring():
            rules, group = self.ec2_backend.authorize_security_group_ingress(**kwargs)
        template = self.response_template(AUTHORIZE_SECURITY_GROUP_INGRESS_RESPONSE)
        return template.render(rules=rules, group=group)

    def create_security_group(self) -> str:
        name = self._get_param("GroupName")
        description = self._get_param("GroupDescription")
        vpc_id = self._get_param("VpcId")
        tags = self._parse_tag_specification().get("security-group", {})

        self.error_on_dryrun()

        group = self.ec2_backend.create_security_group(
            name, description, vpc_id=vpc_id, tags=tags
        )
        template = self.response_template(CREATE_SECURITY_GROUP_RESPONSE)
        return template.render(group=group)

    def delete_security_group(self) -> str:
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

        return DELETE_GROUP_RESPONSE

    def describe_security_groups(self) -> str:
        groupnames = self._get_multi_param("GroupName")
        group_ids = self._get_multi_param("GroupId")
        filters = self._filters_from_querystring()

        groups = self.ec2_backend.describe_security_groups(
            group_ids=group_ids, groupnames=groupnames, filters=filters
        )

        template = self.response_template(DESCRIBE_SECURITY_GROUPS_RESPONSE)
        return template.render(groups=groups)

    def describe_security_group_rules(self) -> str:
        self.error_on_dryrun()
        sg_rule_ids = self._get_multi_param("SecurityGroupRuleId")
        filters = self._filters_from_querystring()
        rules = self.ec2_backend.describe_security_group_rules(sg_rule_ids, filters)
        template = self.response_template(DESCRIBE_SECURITY_GROUP_RULES_RESPONSE)
        return template.render(rules=rules)

    def revoke_security_group_egress(self) -> str:
        self.error_on_dryrun()

        for args in self._process_rules_from_querystring():
            # we don't need this parameter to revoke
            del args["sgrule_tags"]
            self.ec2_backend.revoke_security_group_egress(**args)
        return REVOKE_SECURITY_GROUP_EGRESS_RESPONSE

    def revoke_security_group_ingress(self) -> str:
        self.error_on_dryrun()

        for args in self._process_rules_from_querystring():
            # we don't need this parameter to revoke
            del args["sgrule_tags"]
            self.ec2_backend.revoke_security_group_ingress(**args)
        return REVOKE_SECURITY_GROUP_INGRESS_RESPONSE

    def modify_security_group_rules(self) -> str:
        self.error_on_dryrun()

        rules = {}
        security_group_rules_param = self._get_params()["SecurityGroupRule"]
        for idx, sgr in security_group_rules_param.items():
            rules[sgr["SecurityGroupRuleId"]] = sgr["SecurityGroupRule"]

        group_id = self._get_param("GroupId")
        self.ec2_backend.modify_security_group_rules(group_id, rules)

        return MODIFY_SECURITY_GROUP_RULES_RESPONSE

    def update_security_group_rule_descriptions_ingress(self) -> str:
        for args in self._process_rules_from_querystring():
            # we don't need this parameter to revoke
            del args["sgrule_tags"]
            self.ec2_backend.update_security_group_rule_descriptions_ingress(**args)
        return UPDATE_SECURITY_GROUP_RULE_DESCRIPTIONS_INGRESS

    def update_security_group_rule_descriptions_egress(self) -> str:
        for args in self._process_rules_from_querystring():
            # we don't need this parameter to revoke
            del args["sgrule_tags"]
            self.ec2_backend.update_security_group_rule_descriptions_egress(**args)
        return UPDATE_SECURITY_GROUP_RULE_DESCRIPTIONS_EGRESS


CREATE_SECURITY_GROUP_RESPONSE = """<CreateSecurityGroupResponse xmlns="http://ec2.amazonaws.com/doc/2013-10-15/">
   <requestId>59dbff89-35bd-4eac-99ed-be587EXAMPLE</requestId>
   <return>true</return>
   <groupId>{{ group.id }}</groupId>
   <securityGroupArn>{{ group.arn }}</securityGroupArn>
   <tagSet>
    {% for tag in group.get_tags() %}
        <item>
        <key>{{ tag.key }}</key>
        <value>{{ tag.value }}</value>
        </item>
    {% endfor %}
    </tagSet>
</CreateSecurityGroupResponse>"""

DESCRIBE_SECURITY_GROUP_RULES_RESPONSE = """
<DescribeSecurityGroupRulesResponse xmlns="http://ec2.amazonaws.com/doc/2016-11-15/">
  <requestId>{{ request_id }}</requestId>
  <securityGroupRuleSet>
    {% for rule in rules %}
        <item>
            <ipProtocol>{{ rule.ip_protocol }}</ipProtocol>
            <groupId>{{ rule.group_id }}</groupId>
            <groupOwnerId>{{ rule.owner_id }}</groupOwnerId>
            <fromPort>{{ rule.from_port if rule.from_port is not none else -1 }}</fromPort>
            <toPort>{{ rule.to_port if rule.to_port is not none else -1 }}</toPort>
            <isEgress>{{ 'true' if rule.is_egress else 'false' }}</isEgress>
            <securityGroupRuleId>{{ rule.id }}</securityGroupRuleId>
            {% if rule.ip_range %}
              {% if rule.ip_range['CidrIp'] %}<cidrIpv4>{{ rule.ip_range['CidrIp'] }}</cidrIpv4>{% endif %}
              {% if rule.ip_range['CidrIpv6'] %}<cidrIpv6>{{ rule.ip_range['CidrIpv6'] }}</cidrIpv6>{% endif %}
              {% if rule.ip_range['Description'] %}<description>{{ rule.ip_range['Description'] }}</description>{% endif %}
            {% endif %}
            {% if rule.source_group %}
              <referencedGroupInfo>
                <groupId>{{ rule.source_group['GroupId'] }}</groupId><userId>{{ rule.source_group['OwnerId'] }}</userId>
              </referencedGroupInfo>
              {% if rule.source_group['Description'] %}<description>{{ rule.source_group['Description'] }}</description>{% endif %}
            {% endif %}
            {% if rule.prefix_list_id %}
              <prefixListId>{{ rule.prefix_list_id['PrefixListId'] }}</prefixListId>
              {% if rule.prefix_list_id['Description'] %}<description>{{ rule.prefix_list_id['Description'] }}</description>{% endif %}
            {% endif %}
            <tagSet>
            {% for tag in rule.get_tags() %}
                <item>
                  <key>{{ tag.key }}</key>
                  <value>{{ tag.value }}</value>
                </item>
            {% endfor %}
            </tagSet> 
        </item>
    {% endfor %}
  </securityGroupRuleSet>
</DescribeSecurityGroupRulesResponse>"""

DELETE_GROUP_RESPONSE = """<DeleteSecurityGroupResponse xmlns="http://ec2.amazonaws.com/doc/2013-10-15/">
  <requestId>59dbff89-35bd-4eac-99ed-be587EXAMPLE</requestId>
  <return>true</return>
</DeleteSecurityGroupResponse>"""

DESCRIBE_SECURITY_GROUPS_RESPONSE = """<DescribeSecurityGroupsResponse xmlns="http://ec2.amazonaws.com/doc/2013-10-15/">
   <requestId>59dbff89-35bd-4eac-99ed-be587EXAMPLE</requestId>
   <securityGroupInfo>
      {% for group in groups %}
          <item>
             <ownerId>{{ group.owner_id }}</ownerId>
             <groupId>{{ group.id }}</groupId>
             <groupName>{{ group.name }}</groupName>
             <groupDescription>{{ group.description }}</groupDescription>
             <securityGroupArn>{{ group.arn }}</securityGroupArn>
             {% if group.vpc_id %}
             <vpcId>{{ group.vpc_id }}</vpcId>
             {% endif %}
             <ipPermissions>
               {% for rule in group.flattened_ingress_rules %}
                    <item>
                       <ipProtocol>{{ rule.ip_protocol }}</ipProtocol>
                       {% if rule.from_port is not none %}
                       <fromPort>{{ rule.from_port }}</fromPort>
                       {% endif %}
                       {% if rule.to_port is not none %}
                       <toPort>{{ rule.to_port }}</toPort>
                       {% endif %}
                       <groups>
                          {% for source_group in rule.source_groups %}
                              <item>
                                 {% if source_group.OwnerId and source_group.OwnerId != "" %}
                                 <userId>{{ source_group.OwnerId }}</userId>
                                 {% endif %}
                                 {% if source_group.GroupId and source_group.GroupId != "" %}
                                 <groupId>{{ source_group.GroupId }}</groupId>
                                 {% endif %}
                                 {% if source_group.GroupName and source_group.GroupName != "" %}
                                 <groupName>{{ source_group.GroupName }}</groupName>
                                 {% endif %}
                                 {% if source_group.Description and source_group.Description != "" %}
                                 <description>{{ source_group.Description }}</description>
                                 {% endif %}
                              </item>
                          {% endfor %}
                       </groups>
                       <ipRanges>
                          {% for ip_range in rule.ip_ranges %}
                             {% if ip_range['CidrIp'] %}
                              <item>
                                    <cidrIp>{{ ip_range['CidrIp'] }}</cidrIp>
                                    {% if ip_range['Description'] %}
                                    <description>{{ ip_range['Description'] }}</description>
                                    {% endif %}
                              </item>
                              {% endif %}
                          {% endfor %}
                       </ipRanges>
                       <ipv6Ranges>
                        {% for ip_range in rule.ip_ranges %}
                            {% if ip_range['CidrIpv6'] %}
                            <item>
                                <cidrIpv6>{{ ip_range['CidrIpv6'] }}</cidrIpv6>
                                {% if ip_range['Description'] %}
                                <description>{{ ip_range['Description'] }}</description>
                                {% endif %}
                            </item>
                            {% endif %}
                        {% endfor %}
                        </ipv6Ranges>
                        <prefixListIds>
                            {% for prefix_list in rule.prefix_list_ids %}
                            <item>
                                <prefixListId>{{ prefix_list.PrefixListId }}</prefixListId>
                                {% if prefix_list.Description %}
                                <description>{{ prefix_list.Description }}</description>
                                {% endif %}
                            </item>
                            {% endfor %}
                       </prefixListIds>
                    </item>
                {% endfor %}
             </ipPermissions>
             <ipPermissionsEgress>
               {% for rule in group.flattened_egress_rules %}
                    <item>
                       <ipProtocol>{{ rule.ip_protocol }}</ipProtocol>
                       {% if rule.from_port is not none %}
                       <fromPort>{{ rule.from_port }}</fromPort>
                       {% endif %}
                       {% if rule.to_port is not none %}
                       <toPort>{{ rule.to_port }}</toPort>
                       {% endif %}
                       <groups>
                          {% for source_group in rule.source_groups %}
                              <item>
                                 {% if source_group.OwnerId and source_group.OwnerId != "" %}
                                 <userId>{{ source_group.OwnerId }}</userId>
                                 {% endif %}
                                 {% if source_group.GroupId and source_group.GroupId != "" %}
                                 <groupId>{{ source_group.GroupId }}</groupId>
                                 {% endif %}
                                 {% if source_group.GroupName and source_group.GroupName != "" %}
                                 <groupName>{{ source_group.GroupName }}</groupName>
                                 {% endif %}
                                 {% if source_group.Description and source_group.Description != "" %}
                                 <description>{{ source_group.Description }}</description>
                                 {% endif %}
                              </item>
                          {% endfor %}
                       </groups>
                       <ipRanges>
                          {% for ip_range in rule.ip_ranges %}
                             {% if ip_range['CidrIp'] %}
                              <item>
                                    <cidrIp>{{ ip_range['CidrIp'] }}</cidrIp>
                                    {% if ip_range['Description'] %}
                                    <description>{{ ip_range['Description'] }}</description>
                                    {% endif %}
                              </item>
                              {% endif %}
                          {% endfor %}
                       </ipRanges>
                       <ipv6Ranges>
                        {% for ip_range in rule.ip_ranges %}
                            {% if ip_range['CidrIpv6'] %}
                            <item>
                                    <cidrIpv6>{{ ip_range['CidrIpv6'] }}</cidrIpv6>
                                    {% if ip_range['Description'] %}
                                    <description>{{ ip_range['Description'] }}</description>
                                    {% endif %}
                            </item>
                            {% endif %}
                        {% endfor %}
                        </ipv6Ranges>
                        <prefixListIds>
                            {% for prefix_list in rule.prefix_list_ids %}
                            <item>
                                <prefixListId>{{ prefix_list.PrefixListId }}</prefixListId>
                                {% if prefix_list.Description %}
                                <description>{{ prefix_list.Description }}</description>
                                {% endif %}
                            </item>
                            {% endfor %}
                        </prefixListIds>
                    </item>
               {% endfor %}
             </ipPermissionsEgress>
             <tagSet>
               {% for tag in group.get_tags() %}
                 <item>
                   <key>{{ tag.key }}</key>
                   <value>{{ tag.value }}</value>
                 </item>
               {% endfor %}
             </tagSet>
          </item>
      {% endfor %}
   </securityGroupInfo>
</DescribeSecurityGroupsResponse>"""

AUTHORIZE_SECURITY_GROUP_INGRESS_RESPONSE = """<AuthorizeSecurityGroupIngressResponse xmlns="http://ec2.amazonaws.com/doc/2016-11-15/">
    <requestId>b1f67202-c2c2-4ba4-8464-c8b1d8f5af7a</requestId>
    <return>true</return>
    <securityGroupRuleSet>
    {% for rule in rules %}
        <item>
            {% if rule.ip_range %}
                {% if rule.ip_range.CidrIp %}
                <cidrIpv4>{{ rule.ip_range.CidrIp }}</cidrIpv4>
                {% endif %}
                {% if rule.ip_range.CidrIpv6 %}
                <cidrIpv6>{{ rule.ip_range.CidrIpv6 }}</cidrIpv6>
                {% endif %}
                <description>{{ rule.ip_range.Description or '' }}</description>
            {% endif %}
            {% if rule.prefix_list_id %}
                <prefixListId>{{ rule.prefix_list_id.PrefixListId }}</prefixListId>
                <description>{{ rule.prefix_list_id.Description or '' }}</description>
            {% endif %}
            {% if rule.source_group %}
                <referencedGroupInfo>
                    {% if rule.source_group.OwnerId and rule.source_group.OwnerId != "" %}
                    <userId>{{ rule.source_group.OwnerId }}</userId>
                    {% endif %}
                    {% if rule.source_group.GroupId and rule.source_group.GroupId != "" %}
                    <groupId>{{ rule.source_group.GroupId }}</groupId>
                    {% endif %}
                    {% if rule.source_group.VpcId and rule.source_group.VpcId != "" %}
                    <vpcId>{{ rule.source_group.VpcId }}</vpcId>
                    {% endif %}
                </referencedGroupInfo>
                <description>{{ rule.source_group.Description or '' }}</description>
            {% endif %}
            {% if rule.from_port is not none %}
            <fromPort>{{ rule.from_port }}</fromPort>
            {% endif %}
            <groupId>{{ group.id }}</groupId>
            <groupOwnerId>{{ rule.owner_id }}</groupOwnerId>
            <ipProtocol>{{ rule.ip_protocol }}</ipProtocol>
            <isEgress>false</isEgress>
            <securityGroupRuleId>{{ rule.id }}</securityGroupRuleId>
            {% if rule.to_port is not none %}
            <toPort>{{ rule.to_port }}</toPort>
            {% endif %}
            <tagSet>
                {% for tag in rule.get_tags() %}
                    <item>
                      <key>{{ tag.key }}</key>
                      <value>{{ tag.value }}</value>
                    </item>
                {% endfor %}
            </tagSet> 
        </item>
    {% endfor %}
    </securityGroupRuleSet>
</AuthorizeSecurityGroupIngressResponse>"""

REVOKE_SECURITY_GROUP_INGRESS_RESPONSE = """<RevokeSecurityGroupIngressResponse xmlns="http://ec2.amazonaws.com/doc/2013-10-15/">
  <requestId>59dbff89-35bd-4eac-99ed-be587EXAMPLE</requestId>
  <return>true</return>
</RevokeSecurityGroupIngressResponse>"""

AUTHORIZE_SECURITY_GROUP_EGRESS_RESPONSE = """<AuthorizeSecurityGroupEgressResponse xmlns="http://ec2.amazonaws.com/doc/2016-11-15/">
    <requestId>b1f67202-c2c2-4ba4-8464-c8b1d8f5af7a</requestId>
    <return>true</return>
    <securityGroupRuleSet>
    {% for rule in rules %}
        <item>
            {% if rule.ip_range %}
                {% if rule.ip_range.CidrIp %}
                <cidrIpv4>{{ rule.ip_range.CidrIp }}</cidrIpv4>
                {% endif %}
                {% if rule.ip_range.CidrIpv6 %}
                <cidrIpv6>{{ rule.ip_range.CidrIpv6 }}</cidrIpv6>
                {% endif %}
                <description>{{ rule.ip_range.Description or '' }}</description>
            {% endif %}
            {% if rule.prefix_list_id %}
                <prefixListId>{{ rule.prefix_list_id.PrefixListId }}</prefixListId>
                <description>{{ rule.prefix_list_id.Description or '' }}</description>
            {% endif %}
            {% if rule.source_group %}
                <referencedGroupInfo>
                    {% if rule.source_group.OwnerId and rule.source_group.OwnerId != "" %}
                    <userId>{{ rule.source_group.OwnerId }}</userId>
                    {% endif %}
                    {% if rule.source_group.GroupId and rule.source_group.GroupId != "" %}
                    <groupId>{{ rule.source_group.GroupId }}</groupId>
                    {% endif %}
                    {% if rule.source_group.VpcId and rule.source_group.VpcId != "" %}
                    <vpcId>{{ rule.source_group.VpcId }}</vpcId>
                    {% endif %}
                </referencedGroupInfo>
                <description>{{ rule.source_group.Description or '' }}</description>
            {% endif %}
            {% if rule.from_port is not none %}
            <fromPort>{{ rule.from_port }}</fromPort>
            {% endif %}
            <groupId>{{ group.id }}</groupId>
            <groupOwnerId>{{ rule.owner_id }}</groupOwnerId>
            <ipProtocol>{{ rule.ip_protocol }}</ipProtocol>
            <isEgress>true</isEgress>
            <securityGroupRuleId>{{ rule.id }}</securityGroupRuleId>
            {% if rule.to_port is not none %}
            <toPort>{{ rule.to_port }}</toPort>
            {% endif %}
            <tagSet>
                {% for tag in rule.get_tags() %}
                    <item>
                      <key>{{ tag.key }}</key>
                      <value>{{ tag.value }}</value>
                    </item>
                {% endfor %}
            </tagSet> 
        </item>
    {% endfor %}
    </securityGroupRuleSet>
</AuthorizeSecurityGroupEgressResponse>"""

REVOKE_SECURITY_GROUP_EGRESS_RESPONSE = """<RevokeSecurityGroupEgressResponse xmlns="http://ec2.amazonaws.com/doc/2013-10-15/">
  <requestId>59dbff89-35bd-4eac-99ed-be587EXAMPLE</requestId>
  <return>true</return>
</RevokeSecurityGroupEgressResponse>"""

MODIFY_SECURITY_GROUP_RULES_RESPONSE = """<?xml version="1.0" encoding="UTF-8"?>
<ModifySecurityGroupRulesResponse xmlns="http://ec2.amazonaws.com/doc/2016-11-15/">
  <requestId>04d8f970-9213-41fa-81d2-50825EXAMPLE</requestId>
  <return>true</return>
</ModifySecurityGroupRulesResponse>
"""

UPDATE_SECURITY_GROUP_RULE_DESCRIPTIONS_INGRESS = """<UpdateSecurityGroupRuleDescriptionsIngressResponse xmlns="http://ec2.amazonaws.com/doc/2013-10-15/">
  <requestId>59dbff89-35bd-4eac-99ed-be587EXAMPLE</requestId>
  <return>true</return>
</UpdateSecurityGroupRuleDescriptionsIngressResponse>"""

UPDATE_SECURITY_GROUP_RULE_DESCRIPTIONS_EGRESS = """<UpdateSecurityGroupRuleDescriptionsEgressResponse xmlns="http://ec2.amazonaws.com/doc/2013-10-15/">
  <requestId>59dbff89-35bd-4eac-99ed-be587EXAMPLE</requestId>
  <return>true</return>
</UpdateSecurityGroupRuleDescriptionsEgressResponse>"""
