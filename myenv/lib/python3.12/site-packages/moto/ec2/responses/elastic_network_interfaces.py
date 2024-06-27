from moto.ec2.exceptions import InvalidParameterValueErrorUnknownAttribute
from moto.ec2.utils import add_tag_specification, get_attribute_value

from ._base_response import EC2BaseResponse


class ElasticNetworkInterfaces(EC2BaseResponse):
    def create_network_interface(self) -> str:
        subnet_id = self._get_param("SubnetId")
        private_ip_address = self._get_param("PrivateIpAddress")
        private_ip_addresses = self._get_multi_param("PrivateIpAddresses")
        ipv6_addresses = self._get_multi_param("Ipv6Addresses")
        ipv6_address_count = self._get_int_param("Ipv6AddressCount", 0)
        secondary_ips_count = self._get_param("SecondaryPrivateIpAddressCount")
        groups = self._get_multi_param("SecurityGroupId")
        subnet = self.ec2_backend.get_subnet(subnet_id)
        description = self._get_param("Description")
        tags = add_tag_specification(self._get_multi_param("TagSpecification"))

        self.error_on_dryrun()

        eni = self.ec2_backend.create_network_interface(
            subnet,
            private_ip_address,
            private_ip_addresses,
            groups,
            description,
            tags,
            secondary_ips_count=secondary_ips_count,
            ipv6_addresses=ipv6_addresses,
            ipv6_address_count=ipv6_address_count,
        )
        template = self.response_template(CREATE_NETWORK_INTERFACE_RESPONSE)
        return template.render(eni=eni)

    def delete_network_interface(self) -> str:
        eni_id = self._get_param("NetworkInterfaceId")

        self.error_on_dryrun()

        self.ec2_backend.delete_network_interface(eni_id)
        template = self.response_template(DELETE_NETWORK_INTERFACE_RESPONSE)
        return template.render()

    def describe_network_interface_attribute(self) -> str:
        eni_id = self._get_param("NetworkInterfaceId")
        attribute = self._get_param("Attribute")

        self.error_on_dryrun()

        eni = self.ec2_backend.get_all_network_interfaces([eni_id])[0]

        if attribute == "description":
            template = self.response_template(
                DESCRIBE_NETWORK_INTERFACE_ATTRIBUTE_RESPONSE_DESCRIPTION
            )
        elif attribute == "groupSet":
            template = self.response_template(
                DESCRIBE_NETWORK_INTERFACE_ATTRIBUTE_RESPONSE_GROUPSET
            )
        elif attribute == "sourceDestCheck":
            template = self.response_template(
                DESCRIBE_NETWORK_INTERFACE_ATTRIBUTE_RESPONSE_SOURCEDESTCHECK
            )
        elif attribute == "attachment":
            template = self.response_template(
                DESCRIBE_NETWORK_INTERFACE_ATTRIBUTE_RESPONSE_ATTACHMENT
            )
        else:
            raise InvalidParameterValueErrorUnknownAttribute(attribute)
        return template.render(eni=eni)

    def describe_network_interfaces(self) -> str:
        eni_ids = self._get_multi_param("NetworkInterfaceId")
        filters = self._filters_from_querystring()

        self.error_on_dryrun()

        enis = self.ec2_backend.get_all_network_interfaces(eni_ids, filters)
        template = self.response_template(DESCRIBE_NETWORK_INTERFACES_RESPONSE)
        return template.render(enis=enis)

    def attach_network_interface(self) -> str:
        eni_id = self._get_param("NetworkInterfaceId")
        instance_id = self._get_param("InstanceId")
        device_index = self._get_param("DeviceIndex")

        self.error_on_dryrun()

        attachment_id = self.ec2_backend.attach_network_interface(
            eni_id, instance_id, device_index
        )
        template = self.response_template(ATTACH_NETWORK_INTERFACE_RESPONSE)
        return template.render(attachment_id=attachment_id)

    def detach_network_interface(self) -> str:
        attachment_id = self._get_param("AttachmentId")

        self.error_on_dryrun()

        self.ec2_backend.detach_network_interface(attachment_id)
        template = self.response_template(DETACH_NETWORK_INTERFACE_RESPONSE)
        return template.render()

    def modify_network_interface_attribute(self) -> str:
        eni_id = self._get_param("NetworkInterfaceId")
        group_ids = self._get_multi_param("SecurityGroupId")
        source_dest_check = get_attribute_value("SourceDestCheck", self.querystring)
        description = get_attribute_value("Description", self.querystring)

        self.error_on_dryrun()

        self.ec2_backend.modify_network_interface_attribute(
            eni_id, group_ids, source_dest_check, description
        )
        return MODIFY_NETWORK_INTERFACE_ATTRIBUTE_RESPONSE

    def reset_network_interface_attribute(self) -> str:
        self.error_on_dryrun()

        raise NotImplementedError(
            "ElasticNetworkInterfaces(AmazonVPC).reset_network_interface_attribute is not yet implemented"
        )

    def assign_private_ip_addresses(self) -> str:
        eni_id = self._get_param("NetworkInterfaceId")
        secondary_ips_count = self._get_int_param("SecondaryPrivateIpAddressCount", 0)
        private_ip_addresses = self._get_multi_param("PrivateIpAddress")
        eni = self.ec2_backend.assign_private_ip_addresses(
            eni_id, private_ip_addresses, secondary_ips_count
        )
        template = self.response_template(ASSIGN_PRIVATE_IP_ADDRESSES)
        return template.render(eni=eni)

    def unassign_private_ip_addresses(self) -> str:
        eni_id = self._get_param("NetworkInterfaceId")
        private_ip_address = self._get_multi_param("PrivateIpAddress")
        eni = self.ec2_backend.unassign_private_ip_addresses(eni_id, private_ip_address)
        template = self.response_template(UNASSIGN_PRIVATE_IP_ADDRESSES)
        return template.render(eni=eni)

    def assign_ipv6_addresses(self) -> str:
        eni_id = self._get_param("NetworkInterfaceId")
        ipv6_count = self._get_int_param("Ipv6AddressCount", 0)
        ipv6_addresses = self._get_multi_param("Ipv6Addresses")
        eni = self.ec2_backend.assign_ipv6_addresses(eni_id, ipv6_addresses, ipv6_count)
        template = self.response_template(ASSIGN_IPV6_ADDRESSES)
        return template.render(eni=eni)

    def unassign_ipv6_addresses(self) -> str:
        eni_id = self._get_param("NetworkInterfaceId")
        ips = self._get_multi_param("Ipv6Addresses")
        eni = self.ec2_backend.unassign_ipv6_addresses(eni_id, ips)
        template = self.response_template(UNASSIGN_IPV6_ADDRESSES)
        return template.render(eni=eni, unassigned_ips=ips)


ASSIGN_PRIVATE_IP_ADDRESSES = """<AssignPrivateIpAddressesResponse xmlns="http://ec2.amazonaws.com/doc/2016-11-15/">
  <requestId>3fb591ba-558c-48f8-ae6b-c2f9d6d06425</requestId>
  <networkInterfaceId>{{ eni.id }}</networkInterfaceId>
  <assignedPrivateIpAddressesSet>
    {% for address in eni.private_ip_addresses %}
    <item>
        <privateIpAddress>{{ address.PrivateIpAddress }}</privateIpAddress>
    </item>
    {% endfor %}
  </assignedPrivateIpAddressesSet>
  <return>true</return>
</AssignPrivateIpAddressesResponse>"""


UNASSIGN_PRIVATE_IP_ADDRESSES = """<UnAssignPrivateIpAddressesResponse xmlns="http://ec2.amazonaws.com/doc/2016-11-15/">
  <requestId>3fb591ba-558c-48f8-ae6b-c2f9d6d06425</requestId>
  <networkInterfaceId>{{ eni.id }}</networkInterfaceId>
  <assignedPrivateIpAddressesSet>
    {% for address in eni.private_ip_addresses %}
    <item>
        <privateIpAddress>{{ address.PrivateIpAddress }}</privateIpAddress>
    </item>
    {% endfor %}
  </assignedPrivateIpAddressesSet>
  <return>true</return>
</UnAssignPrivateIpAddressesResponse>"""


ASSIGN_IPV6_ADDRESSES = """<AssignIpv6AddressesResponse xmlns="http://ec2.amazonaws.com/doc/2016-11-15/">
    <requestId>c36d17eb-a0ba-4d38-8727-example</requestId>
    <networkInterfaceId>{{ eni.id }}</networkInterfaceId>
    <assignedIpv6Addresses>
        {% for address in eni.ipv6_addresses %}
        <item>{{address}}</item>
        {% endfor %}
    </assignedIpv6Addresses>
</AssignIpv6AddressesResponse>
"""

UNASSIGN_IPV6_ADDRESSES = """<UnassignIpv6AddressesResponse xmlns="http://ec2.amazonaws.com/doc/2016-11-15/">
    <requestId>94d446d7-fc8e-4918-94f9-example</requestId>
    <networkInterfaceId>{{ eni.id }}</networkInterfaceId>
    <unassignedIpv6Addresses>
        {% for address in unassigned_ips %}
        <item>{{address}}</item>
        {% endfor %}
    </unassignedIpv6Addresses>
</UnassignIpv6AddressesResponse>"""


CREATE_NETWORK_INTERFACE_RESPONSE = """
<CreateNetworkInterfaceResponse xmlns="http://ec2.amazonaws.com/doc/2013-10-15/">
    <requestId>2c6021ec-d705-445a-9780-420d0c7ab793</requestId>
    <networkInterface>
        <networkInterfaceId>{{ eni.id }}</networkInterfaceId>
        <subnetId>{{ eni.subnet.id }}</subnetId>
        <vpcId>{{ eni.subnet.vpc_id }}</vpcId>
        <availabilityZone>{{ eni.subnet.availability_zone }}</availabilityZone>
        {% if eni.description %}
        <description>{{ eni.description }}</description>
        {% endif %}
        <ownerId>{{ eni.owner_id }}</ownerId>
        <requesterId>AIDARCSPW2WNREUEN7XFM</requesterId>
        <requesterManaged>False</requesterManaged>
        <status>{{ eni.status }}</status>
        <macAddress>{{ eni.mac_address }}</macAddress>
        {% if eni.private_ip_address %}
          <privateIpAddress>{{ eni.private_ip_address }}</privateIpAddress>
        {% endif %}
        {% if eni.private_dns_name %}
        <privateDnsName>{{ eni.private_dns_name }}</privateDnsName>
        {% endif %}
        <sourceDestCheck>{{ "true" if eni.source_dest_check == True else "false" }}</sourceDestCheck>
        <groupSet>
        {% for group in eni.group_set %}
            <item>
                <groupId>{{ group.id }}</groupId>
                <groupName>{{ group.name }}</groupName>
             </item>
         {% endfor %}
         </groupSet>
        {% if eni.association %}
        <association>
            <publicIp>{{ eni.public_ip }}</publicIp>
            <ipOwnerId>{{ eni.owner_id }}</ipOwnerId>
            <allocationId>{{ eni.association.allocationId }}</allocationId>
            <associationId>{{ eni.association.associationId }}</associationId>
            <natEnabled>true</natEnabled>
        </association>
        {% endif %}
        <tagSet>
          {% for tag in eni.get_tags() %}
              <item>
                  <key>{{ tag.key }}</key>
                  <value>{{ tag.value }}</value>
              </item>
          {% endfor %}
        </tagSet>
        <privateIpAddressesSet>
          {% for address in eni.private_ip_addresses %}
          <item>
            <privateIpAddress>{{ address.PrivateIpAddress }}</privateIpAddress>
            {% if address.privateDnsName %}
            <privateDnsName>{{ address.PrivateDnsName }}</privateDnsName>
            {% endif %}
            <primary>{{ "true" if address.Primary == True else "false" }}</primary>
          </item>
          {% endfor %}
        </privateIpAddressesSet>
        <ipv6AddressesSet>
        {% for address in eni.ipv6_addresses %}
            <item>
                <ipv6Address>{{address}}</ipv6Address>
            </item>
        {% endfor %}
        </ipv6AddressesSet>
        <interfaceType>{{ eni.interface_type }}</interfaceType>
    </networkInterface>
</CreateNetworkInterfaceResponse>
"""

DESCRIBE_NETWORK_INTERFACES_RESPONSE = """<DescribeNetworkInterfacesResponse xmlns="http://ec2.amazonaws.com/doc/2013-10-15/">
    <requestId>ddb0aaf1-8b65-4f0a-94fa-654b18b8a204</requestId>
    <networkInterfaceSet>
    {% for eni in enis %}
        <item>
            <networkInterfaceId>{{ eni.id }}</networkInterfaceId>
            <subnetId>{{ eni.subnet.id }}</subnetId>
            <vpcId>{{ eni.subnet.vpc_id }}</vpcId>
            <availabilityZone>{{ eni.subnet.availability_zone }}</availabilityZone>
            {% if eni.description %}
            <description>{{ eni.description }}</description>
            {% endif %}
            <ownerId>{{ eni.owner_id }}</ownerId>
            <requesterId>AIDARCSPW2WNREUEN7XFM</requesterId>
            <requesterManaged>False</requesterManaged>
            <status>{{ eni.status }}</status>
            <macAddress>{{ eni.mac_address }}</macAddress>
            {% if eni.private_ip_address %}
            <privateIpAddress>{{ eni.private_ip_address }}</privateIpAddress>
            {% endif %}
            {% if eni.private_dns_name %}
            <privateDnsName>{{ eni.private_dns_name }}</privateDnsName>
            {% endif %}
            <sourceDestCheck>{{ "true" if eni.source_dest_check == True else "false" }}</sourceDestCheck>
            <groupSet>
            {% for group in eni.group_set %}
                <item>
                    <groupId>{{ group.id }}</groupId>
                    <groupName>{{ group.name }}</groupName>
                </item>
            {% endfor %}
            </groupSet>
            {% if eni.association %}
            <association>
                <publicIp>{{ eni.public_ip }}</publicIp>
                <ipOwnerId>{{ eni.owner_id }}</ipOwnerId>
                <allocationId>{{ eni.association.allocationId }}</allocationId>
                <associationId>{{ eni.association.associationId }}</associationId>
                <natEnabled>true</natEnabled>
            </association>
            {% endif %}
            {% if eni.attachment_id %}
            <attachment>
                <attachTime>{{ eni.attach_time }}</attachTime>
                <attachmentId>{{ eni.attachment_id }}</attachmentId>
                <deleteOnTermination>{{ eni.delete_on_termination }}</deleteOnTermination>
                <deviceIndex>{{ eni.device_index }}</deviceIndex>
                <networkCardIndex>0</networkCardIndex>
                <instanceId>{{ eni.instance.id }}</instanceId>
                <instanceOwnerId>{{ eni.instance.owner_id }}</instanceOwnerId>
                <status>attached</status>
            </attachment>
            {% endif %}
            <tagSet>
            {% for tag in eni.get_tags() %}
                <item>
                    <key>{{ tag.key }}</key>
                    <value>{{ tag.value }}</value>
                </item>
            {% endfor %}
            </tagSet>
            <privateIpAddressesSet>
            {% for address in eni.private_ip_addresses %}
            <item>
                <privateIpAddress>{{ address.PrivateIpAddress }}</privateIpAddress>
                {% if address.privateDnsName %}
                <privateDnsName>{{ address.PrivateDnsName }}</privateDnsName>
                {% endif %}
                <primary>{{ "true" if address.Primary == True else "false" }}</primary>
            </item>
            {% endfor %}
            </privateIpAddressesSet>
            <ipv6AddressesSet>
                {% for address in eni.ipv6_addresses %}
                <item>
                    <ipv6Address>{{address}}</ipv6Address>
                </item>
                {% endfor %}
            </ipv6AddressesSet>
            <interfaceType>{{ eni.interface_type }}</interfaceType>
        </item>
    {% endfor %}
    </networkInterfaceSet>
</DescribeNetworkInterfacesResponse>"""

ATTACH_NETWORK_INTERFACE_RESPONSE = """<AttachNetworkInterfaceResponse xmlns="http://ec2.amazonaws.com/doc/2013-10-15/">
  <requestId>59dbff89-35bd-4eac-99ed-be587EXAMPLE</requestId>
  <attachmentId>{{ attachment_id }}</attachmentId>
</AttachNetworkInterfaceResponse>"""

DETACH_NETWORK_INTERFACE_RESPONSE = """<DetachNetworkInterfaceResponse xmlns="http://ec2.amazonaws.com/doc/2013-10-15/">
  <requestId>59dbff89-35bd-4eac-99ed-be587EXAMPLE</requestId>
  <return>true</return>
</DetachNetworkInterfaceResponse>"""

MODIFY_NETWORK_INTERFACE_ATTRIBUTE_RESPONSE = """<ModifyNetworkInterfaceAttributeResponse xmlns="http://ec2.amazonaws.com/doc/2013-10-15/">
  <requestId>59dbff89-35bd-4eac-99ed-be587EXAMPLE</requestId>
  <return>true</return>
</ModifyNetworkInterfaceAttributeResponse>"""

DELETE_NETWORK_INTERFACE_RESPONSE = """
<DeleteNetworkInterfaceResponse xmlns="http://ec2.amazonaws.com/doc/2013-10-15/">
    <requestId>34b5b3b4-d0c5-49b9-b5e2-a468ef6adcd8</requestId>
    <return>true</return>
</DeleteNetworkInterfaceResponse>"""

DESCRIBE_NETWORK_INTERFACE_ATTRIBUTE_RESPONSE_DESCRIPTION = """
<DescribeNetworkInterfaceAttributeResponse xmlns="http://ec2.amazonaws.com/doc/2013-10-15/">
    <networkInterfaceId>{{ eni.id }}</networkInterfaceId>
    <description>
        <value>{{ eni.description }}</value>
    </description>
</DescribeNetworkInterfaceAttributeResponse>"""

DESCRIBE_NETWORK_INTERFACE_ATTRIBUTE_RESPONSE_GROUPSET = """
<DescribeNetworkInterfaceAttributeResponse xmlns="http://ec2.amazonaws.com/doc/2013-10-15/">
    <networkInterfaceId>{{ eni.id }}</networkInterfaceId>
    <groupSet>
        {% for group in eni.group_set %}
        <item>
            <groupId>{{ group.id }}</groupId>
            <groupName>{{ group.name }}</groupName>
        </item>
        {% endfor %}
    </groupSet>
</DescribeNetworkInterfaceAttributeResponse>"""

DESCRIBE_NETWORK_INTERFACE_ATTRIBUTE_RESPONSE_SOURCEDESTCHECK = """
<DescribeNetworkInterfaceAttributeResponse xmlns="http://ec2.amazonaws.com/doc/2013-10-15/">
    <networkInterfaceId>{{ eni.id }}</networkInterfaceId>
    <sourceDestCheck>
        <value>{{ "true" if eni.source_dest_check == True else "false" }}</value>
    </sourceDestCheck>
</DescribeNetworkInterfaceAttributeResponse>"""

DESCRIBE_NETWORK_INTERFACE_ATTRIBUTE_RESPONSE_ATTACHMENT = """
<DescribeNetworkInterfaceAttributeResponse xmlns="http://ec2.amazonaws.com/doc/2013-10-15/">
    <networkInterfaceId>{{ eni.id }}</networkInterfaceId>
    {% if eni.attachment_id %}
    <attachment>
        <attachTime>{{ eni.attach_time }}</attachTime>
        <attachmentId>{{ eni.attachment_id }}</attachmentId>
        <deleteOnTermination>{{ eni.delete_on_termination }}</deleteOnTermination>
        <deviceIndex>{{ eni.device_index }}</deviceIndex>
        <networkCardIndex>0</networkCardIndex>
        <instanceId>{{ eni.instance.id }}</instanceId>
        <instanceOwnerId>{{ eni.instance.owner_id }}</instanceOwnerId>
        <status>attached</status>
    </attachment>
    {% endif %}
</DescribeNetworkInterfaceAttributeResponse>"""
