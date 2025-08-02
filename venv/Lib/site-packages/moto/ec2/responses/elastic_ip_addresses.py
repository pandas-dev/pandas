from moto.ec2.utils import add_tag_specification

from ._base_response import EC2BaseResponse


class ElasticIPAddresses(EC2BaseResponse):
    def allocate_address(self) -> str:
        domain = self._get_param("Domain", if_none=None)
        reallocate_address = self._get_param("Address", if_none=None)
        tag_param = self._get_multi_param("TagSpecification")
        tags = add_tag_specification(tag_param)

        self.error_on_dryrun()

        if reallocate_address:
            address = self.ec2_backend.allocate_address(
                domain, address=reallocate_address, tags=tags
            )
        else:
            address = self.ec2_backend.allocate_address(domain, tags=tags)
        template = self.response_template(ALLOCATE_ADDRESS_RESPONSE)
        return template.render(address=address)

    def associate_address(self) -> str:
        instance = eni = None

        if "InstanceId" in self.querystring:
            instance = self.ec2_backend.get_instance(self._get_param("InstanceId"))
        elif "NetworkInterfaceId" in self.querystring:
            eni = self.ec2_backend.get_network_interface(
                self._get_param("NetworkInterfaceId")
            )
        else:
            self.ec2_backend.raise_error(
                "MissingParameter",
                "Invalid request, expect InstanceId/NetworkId parameter.",
            )

        reassociate = False
        if "AllowReassociation" in self.querystring:
            reassociate = self._get_param("AllowReassociation") == "true"

        self.error_on_dryrun()

        if instance or eni:
            if "PublicIp" in self.querystring:
                eip = self.ec2_backend.associate_address(
                    instance=instance,
                    eni=eni,
                    address=self._get_param("PublicIp"),
                    reassociate=reassociate,
                )
            elif "AllocationId" in self.querystring:
                eip = self.ec2_backend.associate_address(
                    instance=instance,
                    eni=eni,
                    allocation_id=self._get_param("AllocationId"),
                    reassociate=reassociate,
                )
            else:
                self.ec2_backend.raise_error(
                    "MissingParameter",
                    "Invalid request, expect PublicIp/AllocationId parameter.",
                )
        else:
            self.ec2_backend.raise_error(
                "MissingParameter",
                "Invalid request, expect either instance or ENI.",
            )

        template = self.response_template(ASSOCIATE_ADDRESS_RESPONSE)
        return template.render(address=eip)

    def describe_addresses(self) -> str:
        self.error_on_dryrun()
        allocation_ids = self._get_multi_param("AllocationId")
        public_ips = self._get_multi_param("PublicIp")
        filters = self._filters_from_querystring()
        addresses = self.ec2_backend.describe_addresses(
            allocation_ids, public_ips, filters
        )
        template = self.response_template(DESCRIBE_ADDRESS_RESPONSE)
        return template.render(addresses=addresses)

    def describe_addresses_attribute(self) -> str:
        self.error_on_dryrun()
        allocation_ids = self._get_multi_param("AllocationId")
        addresses = self.ec2_backend.describe_addresses_attribute(allocation_ids)
        template = self.response_template(DESCRIBE_ADDRESS_ATTRIBUTE_RESPONSE)
        return template.render(addresses=addresses)

    def disassociate_address(self) -> str:
        if (
            "PublicIp" not in self.querystring
            and "AssociationId" not in self.querystring
        ):
            self.ec2_backend.raise_error(
                "MissingParameter",
                "Invalid request, expect PublicIp/AssociationId parameter.",
            )

        self.error_on_dryrun()

        if "PublicIp" in self.querystring:
            self.ec2_backend.disassociate_address(address=self._get_param("PublicIp"))
        elif "AssociationId" in self.querystring:
            self.ec2_backend.disassociate_address(
                association_id=self._get_param("AssociationId")
            )

        return self.response_template(DISASSOCIATE_ADDRESS_RESPONSE).render()

    def release_address(self) -> str:
        if (
            "PublicIp" not in self.querystring
            and "AllocationId" not in self.querystring
        ):
            self.ec2_backend.raise_error(
                "MissingParameter",
                "Invalid request, expect PublicIp/AllocationId parameter.",
            )

        self.error_on_dryrun()

        if "PublicIp" in self.querystring:
            self.ec2_backend.release_address(address=self._get_param("PublicIp"))
        elif "AllocationId" in self.querystring:
            self.ec2_backend.release_address(
                allocation_id=self._get_param("AllocationId")
            )

        return self.response_template(RELEASE_ADDRESS_RESPONSE).render()


ALLOCATE_ADDRESS_RESPONSE = """<AllocateAddressResponse xmlns="http://ec2.amazonaws.com/doc/2013-10-15/">
  <requestId>59dbff89-35bd-4eac-99ed-be587EXAMPLE</requestId>
  <publicIp>{{ address.public_ip }}</publicIp>
  <domain>{{ address.domain }}</domain>
  {% if address.allocation_id %}
    <allocationId>{{ address.allocation_id }}</allocationId>
  {% endif %}
</AllocateAddressResponse>"""

ASSOCIATE_ADDRESS_RESPONSE = """<AssociateAddressResponse xmlns="http://ec2.amazonaws.com/doc/2013-10-15/">
  <requestId>59dbff89-35bd-4eac-99ed-be587EXAMPLE</requestId>
  <return>true</return>
  {% if address.association_id %}
    <associationId>{{ address.association_id }}</associationId>
  {% endif %}
</AssociateAddressResponse>"""

DESCRIBE_ADDRESS_RESPONSE = """<DescribeAddressesResponse xmlns="http://ec2.amazonaws.com/doc/2013-10-15/">
  <requestId>59dbff89-35bd-4eac-99ed-be587EXAMPLE</requestId>
  <addressesSet>
    {% for address in addresses %}
        <item>
          <publicIp>{{ address.public_ip }}</publicIp>
          <domain>{{ address.domain }}</domain>
          {% if address.instance %}
            <instanceId>{{ address.instance.id }}</instanceId>
          {% else %}
            <instanceId/>
          {% endif %}
          {% if address.eni %}
            <networkInterfaceId>{{ address.eni.id }}</networkInterfaceId>
            <privateIpAddress>{{ address.eni.private_ip_address }}</privateIpAddress>
            <networkInterfaceOwnerId>{{ address.eni.owner_id }}</networkInterfaceOwnerId>
          {% else %}
            <networkInterfaceId/>
          {% endif %}
          {% if address.allocation_id %}
            <allocationId>{{ address.allocation_id }}</allocationId>
          {% endif %}
          {% if address.association_id %}
            <associationId>{{ address.association_id }}</associationId>
          {% endif %}
          <tagSet>
          {% for tag in address.get_tags() %}
              <item>
                  <key>{{ tag.key }}</key>
                  <value>{{ tag.value }}</value>
              </item>
          {% endfor %}
          </tagSet>
        </item>
    {% endfor %}
  </addressesSet>
</DescribeAddressesResponse>"""

DESCRIBE_ADDRESS_ATTRIBUTE_RESPONSE = """<DescribeAddressesAttributeResponse xmlns="http://ec2.amazonaws.com/doc/2013-10-15/">
  <requestId>59dbff89-35bd-4eac-99ed-be587EXAMPLE</requestId>
  <addressSet>
    {% for address in addresses %}
        <item>
          <publicIp>{{ address.public_ip }}</publicIp>
          {% if address.allocation_id %}
            <allocationId>{{ address.allocation_id }}</allocationId>
          {% endif %}
          {% if address.ptrRecord %}
            <ptrRecord>{{ address.ptrRecord }}</ptrRecord>
          {% endif %}
          {% if address.ptrRecordUpdate %}
            <ptrRecordUpdate>{{ address.ptrRecordUpdate }}</ptrRecordUpdate>
          {% endif %}
        </item>
    {% endfor %}
  </addressSet>
</DescribeAddressesAttributeResponse>"""

DISASSOCIATE_ADDRESS_RESPONSE = """<DisassociateAddressResponse xmlns="http://ec2.amazonaws.com/doc/2013-10-15/">
  <requestId>59dbff89-35bd-4eac-99ed-be587EXAMPLE</requestId>
  <return>true</return>
</DisassociateAddressResponse>"""

RELEASE_ADDRESS_RESPONSE = """<ReleaseAddressResponse xmlns="http://ec2.amazonaws.com/doc/2013-10-15/">
  <requestId>59dbff89-35bd-4eac-99ed-be587EXAMPLE</requestId>
  <return>true</return>
</ReleaseAddressResponse>"""
