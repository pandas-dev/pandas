from ._base_response import EC2BaseResponse


class NetworkACLs(EC2BaseResponse):
    def create_network_acl(self) -> str:
        vpc_id = self._get_param("VpcId")
        tags = self._get_multi_param("TagSpecification")
        if tags:
            tags = tags[0].get("Tag")
        network_acl = self.ec2_backend.create_network_acl(vpc_id, tags=tags)
        template = self.response_template(CREATE_NETWORK_ACL_RESPONSE)
        return template.render(network_acl=network_acl)

    def create_network_acl_entry(self) -> str:
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

        network_acl_entry = self.ec2_backend.create_network_acl_entry(
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

        template = self.response_template(CREATE_NETWORK_ACL_ENTRY_RESPONSE)
        return template.render(network_acl_entry=network_acl_entry)

    def delete_network_acl(self) -> str:
        network_acl_id = self._get_param("NetworkAclId")
        self.ec2_backend.delete_network_acl(network_acl_id)
        template = self.response_template(DELETE_NETWORK_ACL_ASSOCIATION)
        return template.render()

    def delete_network_acl_entry(self) -> str:
        network_acl_id = self._get_param("NetworkAclId")
        rule_number = self._get_param("RuleNumber")
        egress = self._get_param("Egress")
        self.ec2_backend.delete_network_acl_entry(network_acl_id, rule_number, egress)
        template = self.response_template(DELETE_NETWORK_ACL_ENTRY_RESPONSE)
        return template.render()

    def replace_network_acl_entry(self) -> str:
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

        template = self.response_template(REPLACE_NETWORK_ACL_ENTRY_RESPONSE)
        return template.render()

    def describe_network_acls(self) -> str:
        network_acl_ids = self._get_multi_param("NetworkAclId")
        filters = self._filters_from_querystring()
        network_acls = self.ec2_backend.describe_network_acls(network_acl_ids, filters)
        template = self.response_template(DESCRIBE_NETWORK_ACL_RESPONSE)
        return template.render(network_acls=network_acls)

    def replace_network_acl_association(self) -> str:
        association_id = self._get_param("AssociationId")
        network_acl_id = self._get_param("NetworkAclId")

        association = self.ec2_backend.replace_network_acl_association(
            association_id, network_acl_id
        )
        template = self.response_template(REPLACE_NETWORK_ACL_ASSOCIATION)
        return template.render(association=association)


CREATE_NETWORK_ACL_RESPONSE = """
<CreateNetworkAclResponse xmlns="http://ec2.amazonaws.com/doc/2013-10-15/">
   <requestId>59dbff89-35bd-4eac-99ed-be587EXAMPLE</requestId>
   <networkAcl>
      <networkAclId>{{ network_acl.id }}</networkAclId>
      <vpcId>{{ network_acl.vpc_id }}</vpcId>
      <default>false</default>
      <entrySet/>
      <associationSet/>
      <tagSet>
      {% for tag in network_acl.get_tags() %}
        <item>
          <resourceId>{{ tag.resource_id }}</resourceId>
          <resourceType>{{ tag.resource_type }}</resourceType>
          <key>{{ tag.key }}</key>
          <value>{{ tag.value }}</value>
        </item>
      {% endfor %}
      </tagSet>
   </networkAcl>
</CreateNetworkAclResponse>
"""

DESCRIBE_NETWORK_ACL_RESPONSE = """
<DescribeNetworkAclsResponse xmlns="http://ec2.amazonaws.com/doc/2013-10-15/">
   <requestId>59dbff89-35bd-4eac-99ed-be587EXAMPLE</requestId>
   <networkAclSet>
   {% for network_acl in network_acls %}
   <item>
     <networkAclId>{{ network_acl.id }}</networkAclId>
     <vpcId>{{ network_acl.vpc_id }}</vpcId>
     <ownerId>{{ network_acl.owner_id }}</ownerId>
     <default>{{ network_acl.default }}</default>
     <entrySet>
       {% for entry in network_acl.network_acl_entries %}
         <item>
           <ruleNumber>{{ entry.rule_number }}</ruleNumber>
           <protocol>{{ entry.protocol }}</protocol>
           <ruleAction>{{ entry.rule_action }}</ruleAction>
           <egress>{{ entry.egress.lower() }}</egress>
           {% if entry.cidr_block %}<cidrBlock>{{ entry.cidr_block }}</cidrBlock>{% endif %}
           {% if entry.ipv6_cidr_block %}<ipv6CidrBlock>{{ entry.ipv6_cidr_block }}</ipv6CidrBlock>{% endif %}
           {% if entry.icmp_code or entry.icmp_type %}<icmpTypeCode>
             {% if entry.icmp_code %}<code>{{ entry.icmp_code }}</code>{% endif %}
             {% if entry.icmp_type %}<type>{{ entry.icmp_type }}</type>{% endif %}
           </icmpTypeCode>{% endif %}
           {% if entry.port_range_from or entry.port_range_to %}
             <portRange>
               <from>{{ entry.port_range_from }}</from>
               <to>{{ entry.port_range_to }}</to>
             </portRange>
           {% endif %}
         </item>
       {% endfor %}
     </entrySet>
     <associationSet>
       {% for association in network_acl.associations.values() %}
         <item>
           <networkAclAssociationId>{{ association.id }}</networkAclAssociationId>
           <networkAclId>{{ association.network_acl_id }}</networkAclId>
           <subnetId>{{ association.subnet_id }}</subnetId>
         </item>
       {% endfor %}
     </associationSet>
     <tagSet>
      {% for tag in network_acl.get_tags() %}
        <item>
          <resourceId>{{ tag.resource_id }}</resourceId>
          <resourceType>{{ tag.resource_type }}</resourceType>
          <key>{{ tag.key}}</key>
          <value>{{ tag.value }}</value>
        </item>
      {% endfor %}
     </tagSet>
   </item>
   {% endfor %}
 </networkAclSet>
</DescribeNetworkAclsResponse>
"""

CREATE_NETWORK_ACL_ENTRY_RESPONSE = """
<CreateNetworkAclEntryResponse xmlns="http://ec2.amazonaws.com/doc/2013-10-15/">
   <requestId>59dbff89-35bd-4eac-99ed-be587EXAMPLE</requestId>
   <return>true</return>
</CreateNetworkAclEntryResponse>
"""

REPLACE_NETWORK_ACL_ENTRY_RESPONSE = """
<ReplaceNetworkAclEntryResponse xmlns="http://ec2.amazonaws.com/doc/2013-10-15/">
   <requestId>59dbff89-35bd-4eac-99ed-be587EXAMPLE</requestId>
   <return>true</return>
</ReplaceNetworkAclEntryResponse>
"""

REPLACE_NETWORK_ACL_ASSOCIATION = """
<ReplaceNetworkAclAssociationResponse xmlns="http://ec2.amazonaws.com/doc/2013-10-15/">
   <requestId>59dbff89-35bd-4eac-99ed-be587EXAMPLE</requestId>
   <newAssociationId>{{ association.new_association_id }}</newAssociationId>
</ReplaceNetworkAclAssociationResponse>
"""

DELETE_NETWORK_ACL_ASSOCIATION = """
<DeleteNetworkAclResponse xmlns="http://ec2.amazonaws.com/doc/2013-10-15/">
   <requestId>59dbff89-35bd-4eac-99ed-be587EXAMPLE</requestId>
   <return>true</return>
</DeleteNetworkAclResponse>
"""

DELETE_NETWORK_ACL_ENTRY_RESPONSE = """
<DeleteNetworkAclEntryResponse xmlns="http://ec2.amazonaws.com/doc/2013-10-15/">
   <requestId>59dbff89-35bd-4eac-99ed-be587EXAMPLE</requestId>
   <return>true</return>
</DeleteNetworkAclEntryResponse>
"""
