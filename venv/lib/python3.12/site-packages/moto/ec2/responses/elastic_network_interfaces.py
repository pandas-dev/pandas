from moto.core.responses import ActionResult, EmptyResult
from moto.ec2.exceptions import InvalidParameterValueErrorUnknownAttribute
from moto.ec2.utils import add_tag_specification, get_attribute_value

from ._base_response import EC2BaseResponse


class ElasticNetworkInterfaces(EC2BaseResponse):
    def create_network_interface(self) -> ActionResult:
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
        result = {"NetworkInterface": eni}
        return ActionResult(result)

    def delete_network_interface(self) -> ActionResult:
        eni_id = self._get_param("NetworkInterfaceId")

        self.error_on_dryrun()

        self.ec2_backend.delete_network_interface(eni_id)
        return EmptyResult()

    def describe_network_interface_attribute(self) -> ActionResult:
        eni_id = self._get_param("NetworkInterfaceId")
        attribute = self._get_param("Attribute")

        self.error_on_dryrun()

        eni = self.ec2_backend.get_all_network_interfaces([eni_id])[0]

        result = {"NetworkInterfaceId": eni.id}
        if attribute == "description":
            result["Description"] = {"Value": eni.description}
        elif attribute == "groupSet":
            result["Groups"] = eni.group_set
        elif attribute == "sourceDestCheck":
            result["SourceDestCheck"] = {"Value": eni.source_dest_check}
        elif attribute == "attachment":
            result["Attachment"] = eni.attachment
        else:
            raise InvalidParameterValueErrorUnknownAttribute(attribute)
        return ActionResult(result)

    def describe_network_interfaces(self) -> ActionResult:
        eni_ids = self._get_multi_param("NetworkInterfaceId")
        filters = self._filters_from_querystring()

        self.error_on_dryrun()

        enis = self.ec2_backend.get_all_network_interfaces(eni_ids, filters)
        result = {"NetworkInterfaces": enis}
        return ActionResult(result)

    def attach_network_interface(self) -> ActionResult:
        eni_id = self._get_param("NetworkInterfaceId")
        instance_id = self._get_param("InstanceId")
        device_index = self._get_param("DeviceIndex")

        self.error_on_dryrun()

        attachment = self.ec2_backend.attach_network_interface(
            eni_id, instance_id, device_index
        )
        result = {
            "AttachmentId": attachment.attachment_id,
            "NetworkCardIndex": attachment.network_card_index,
        }
        return ActionResult(result)

    def detach_network_interface(self) -> ActionResult:
        attachment_id = self._get_param("AttachmentId")

        self.error_on_dryrun()

        self.ec2_backend.detach_network_interface(attachment_id)
        return EmptyResult()

    def modify_network_interface_attribute(self) -> ActionResult:
        eni_id = self._get_param("NetworkInterfaceId")
        group_ids = self._get_multi_param("SecurityGroupId")
        source_dest_check = get_attribute_value("SourceDestCheck", self.querystring)
        description = get_attribute_value("Description", self.querystring)

        self.error_on_dryrun()

        self.ec2_backend.modify_network_interface_attribute(
            eni_id, group_ids, source_dest_check, description
        )
        return EmptyResult()

    def reset_network_interface_attribute(self) -> str:
        self.error_on_dryrun()

        raise NotImplementedError(
            "ElasticNetworkInterfaces(AmazonVPC).reset_network_interface_attribute is not yet implemented"
        )

    def assign_private_ip_addresses(self) -> ActionResult:
        eni_id = self._get_param("NetworkInterfaceId")
        secondary_ips_count = self._get_int_param("SecondaryPrivateIpAddressCount", 0)
        private_ip_addresses = self._get_multi_param("PrivateIpAddress")
        eni = self.ec2_backend.assign_private_ip_addresses(
            eni_id, private_ip_addresses, secondary_ips_count
        )
        result = {
            "NetworkInterfaceId": eni.id,
            "AssignedPrivateIpAddresses": eni.private_ip_addresses,
        }
        return ActionResult(result)

    def unassign_private_ip_addresses(self) -> ActionResult:
        eni_id = self._get_param("NetworkInterfaceId")
        private_ip_address = self._get_multi_param("PrivateIpAddress")
        self.ec2_backend.unassign_private_ip_addresses(eni_id, private_ip_address)
        return EmptyResult()

    def assign_ipv6_addresses(self) -> ActionResult:
        eni_id = self._get_param("NetworkInterfaceId")
        ipv6_count = self._get_int_param("Ipv6AddressCount", 0)
        ipv6_addresses = self._get_multi_param("Ipv6Addresses")
        eni = self.ec2_backend.assign_ipv6_addresses(eni_id, ipv6_addresses, ipv6_count)
        result = {
            "AssignedIpv6Addresses": eni.ipv6_addresses,
            "NetworkInterfaceId": eni.id,
        }
        return ActionResult(result)

    def unassign_ipv6_addresses(self) -> ActionResult:
        eni_id = self._get_param("NetworkInterfaceId")
        ips = self._get_multi_param("Ipv6Addresses")
        unassigned_ipv6_addresses = self.ec2_backend.unassign_ipv6_addresses(
            eni_id, ips
        )
        result = {
            "NetworkInterfaceId": eni_id,
            "UnassignedIpv6Addresses": unassigned_ipv6_addresses,
        }
        return ActionResult(result)
