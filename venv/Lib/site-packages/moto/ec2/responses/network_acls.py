from moto.core.responses import ActionResult, EmptyResult

from ._base_response import EC2BaseResponse


class NetworkACLs(EC2BaseResponse):
    def create_network_acl(self) -> ActionResult:
        vpc_id = self._get_param("VpcId")
        tags = self._get_param("TagSpecifications", [])
        if tags:
            tags = tags[0].get("Tags")
        network_acl = self.ec2_backend.create_network_acl(vpc_id, tags=tags)
        return ActionResult({"NetworkAcl": network_acl})

    def create_network_acl_entry(self) -> ActionResult:
        network_acl_id = self._get_param("NetworkAclId")
        rule_number = self._get_param("RuleNumber")
        protocol = self._get_param("Protocol")
        rule_action = self._get_param("RuleAction")
        egress = self._get_param("Egress")
        cidr_block = self._get_param("CidrBlock")
        icmp_code = self._get_param("Icmp.Code")
        icmp_type = self._get_param("Icmp.Type")
        port_range_from = self._get_param("PortRange.From")
        port_range_to = self._get_param("PortRange.To")
        ipv6_cidr_block = self._get_param("Ipv6CidrBlock")

        self.ec2_backend.create_network_acl_entry(
            network_acl_id=network_acl_id,
            rule_number=rule_number,
            protocol=protocol,
            rule_action=rule_action,
            egress=egress,
            cidr_block=cidr_block,
            icmp_code=icmp_code,
            icmp_type=icmp_type,
            port_range_from=port_range_from,
            port_range_to=port_range_to,
            ipv6_cidr_block=ipv6_cidr_block,
        )

        return EmptyResult()

    def delete_network_acl(self) -> ActionResult:
        network_acl_id = self._get_param("NetworkAclId")
        self.ec2_backend.delete_network_acl(network_acl_id)
        return EmptyResult()

    def delete_network_acl_entry(self) -> ActionResult:
        network_acl_id = self._get_param("NetworkAclId")
        rule_number = self._get_param("RuleNumber")
        egress = self._get_param("Egress")
        self.ec2_backend.delete_network_acl_entry(network_acl_id, rule_number, egress)
        return EmptyResult()

    def replace_network_acl_entry(self) -> ActionResult:
        network_acl_id = self._get_param("NetworkAclId")
        rule_number = self._get_param("RuleNumber")
        protocol = self._get_param("Protocol")
        rule_action = self._get_param("RuleAction")
        egress = self._get_param("Egress")
        cidr_block = self._get_param("CidrBlock")
        icmp_code = self._get_param("Icmp.Code")
        icmp_type = self._get_param("Icmp.Type")
        port_range_from = self._get_param("PortRange.From")
        port_range_to = self._get_param("PortRange.To")
        ipv6_cidr_block = self._get_param("Ipv6CidrBlock")

        self.ec2_backend.replace_network_acl_entry(
            network_acl_id=network_acl_id,
            rule_number=rule_number,
            protocol=protocol,
            rule_action=rule_action,
            egress=egress,
            cidr_block=cidr_block,
            icmp_code=icmp_code,
            icmp_type=icmp_type,
            port_range_from=port_range_from,
            port_range_to=port_range_to,
            ipv6_cidr_block=ipv6_cidr_block,
        )

        return EmptyResult()

    def describe_network_acls(self) -> ActionResult:
        network_acl_ids = self._get_param("NetworkAclIds", [])
        filters = self._filters_from_querystring()
        network_acls = self.ec2_backend.describe_network_acls(network_acl_ids, filters)
        return ActionResult({"NetworkAcls": network_acls})

    def replace_network_acl_association(self) -> ActionResult:
        association_id = self._get_param("AssociationId")
        network_acl_id = self._get_param("NetworkAclId")

        association = self.ec2_backend.replace_network_acl_association(
            association_id, network_acl_id
        )
        return ActionResult({"NewAssociationId": association.new_association_id})
