from typing import TYPE_CHECKING, Any, Dict, List, Optional, Union

from moto.core.common_models import CloudFormationModel

from ..exceptions import (
    InvalidNetworkAttachmentIdError,
    InvalidNetworkInterfaceIdError,
    LastEniDetachError,
)
from .core import TaggedEC2Resource
from .security_groups import SecurityGroup

if TYPE_CHECKING:
    from .instances import Instance
from ..utils import (
    generate_dns_from_ip,
    generic_filter,
    random_eni_id,
    random_mac_address,
    random_private_ip,
    random_public_ip,
)


class NetworkInterface(TaggedEC2Resource, CloudFormationModel):
    def __init__(
        self,
        ec2_backend: Any,
        subnet: Any,
        private_ip_address: Union[List[str], str],
        private_ip_addresses: Optional[List[Dict[str, Any]]] = None,
        device_index: int = 0,
        public_ip_auto_assign: bool = False,
        group_ids: Optional[List[str]] = None,
        description: Optional[str] = None,
        tags: Optional[Dict[str, str]] = None,
        **kwargs: Any,
    ):
        self.ec2_backend = ec2_backend
        self.id = random_eni_id()
        self.device_index: Optional[int] = device_index
        if isinstance(private_ip_address, list) and private_ip_address:
            private_ip_address = private_ip_address[0]
        self.private_ip_address: Optional[str] = private_ip_address or None  # type: ignore
        self.private_ip_addresses: List[Dict[str, Any]] = private_ip_addresses or []
        self.ipv6_addresses = kwargs.get("ipv6_addresses") or []

        self.subnet = subnet
        if isinstance(subnet, str):
            self.subnet = self.ec2_backend.get_subnet(subnet)
        self.instance: Optional[Instance] = None
        self.attachment_id: Optional[str] = None
        self.attach_time: Optional[str] = None
        self.delete_on_termination = False
        self.description = description
        self.source_dest_check = True

        self.public_ip: Optional[str] = None
        self.public_ip_auto_assign = public_ip_auto_assign
        self.start()
        self.add_tags(tags or {})
        self.status = "available"
        self.mac_address = random_mac_address()
        self.interface_type = "interface"
        # Local set to the ENI. When attached to an instance, @property group_set
        #   returns groups for both self and the attached instance.
        self._group_set = []

        if self.subnet.ipv6_cidr_block_associations:
            association = list(self.subnet.ipv6_cidr_block_associations.values())[0]
            subnet_ipv6_cidr_block = association.get("ipv6CidrBlock")
            if kwargs.get("ipv6_address_count"):
                while len(self.ipv6_addresses) < kwargs["ipv6_address_count"]:
                    ip = random_private_ip(subnet_ipv6_cidr_block, ipv6=True)
                    if ip not in self.ipv6_addresses:
                        self.ipv6_addresses.append(ip)

        if self.private_ip_addresses:
            primary_selected = True if private_ip_address else False
            for item in self.private_ip_addresses.copy():
                if isinstance(item, str):
                    self.private_ip_addresses.remove(item)
                    self.private_ip_addresses.append(
                        {
                            "Primary": True if not primary_selected else False,
                            "PrivateIpAddress": item,
                        }
                    )
                    primary_selected = True

        if not self.private_ip_address:
            if self.private_ip_addresses:
                for private_ip in self.private_ip_addresses:
                    if isinstance(private_ip, dict) and private_ip.get("Primary"):
                        self.private_ip_address = private_ip.get("PrivateIpAddress")
                        break
            if not self.private_ip_addresses:
                self.private_ip_address = random_private_ip(self.subnet.cidr_block)

        if not self.private_ip_addresses:
            self.private_ip_addresses.append(
                {"Primary": True, "PrivateIpAddress": self.private_ip_address}
            )

        secondary_ips = kwargs.get("secondary_ips_count", None)
        if secondary_ips:
            ips = [
                random_private_ip(self.subnet.cidr_block)
                for index in range(0, int(secondary_ips))
            ]
            if ips:
                self.private_ip_addresses.extend(
                    [{"Primary": False, "PrivateIpAddress": ip} for ip in ips]
                )

        if self.subnet:
            vpc = self.ec2_backend.get_vpc(self.subnet.vpc_id)
            if vpc and vpc.enable_dns_hostnames:
                self.private_dns_name = generate_dns_from_ip(self.private_ip_address)
                for address in self.private_ip_addresses:
                    if address.get("Primary", None):
                        address["PrivateDnsName"] = self.private_dns_name

        group = None
        if group_ids:
            for group_id in group_ids:
                group = self.ec2_backend.get_security_group_from_id(group_id)
                if not group:
                    # Create with specific group ID.
                    group = SecurityGroup(
                        self.ec2_backend,
                        group_id,
                        group_id,
                        group_id,
                        vpc_id=subnet.vpc_id,
                    )
                    self.ec2_backend.groups[subnet.vpc_id][group_id] = group
                if group:
                    self._group_set.append(group)
        if not group_ids:
            group = self.ec2_backend.get_default_security_group(vpc.id)
            if group:
                self._group_set.append(group)

    @property
    def owner_id(self) -> str:
        return self.ec2_backend.account_id

    @property
    def association(self) -> Dict[str, Any]:  # type: ignore[misc]
        association: Dict[str, Any] = {}
        if self.public_ip:
            eips = self.ec2_backend.address_by_ip(
                [self.public_ip], fail_if_not_found=False
            )
            eip = eips[0] if len(eips) > 0 else None
            if eip:
                association["allocationId"] = eip.allocation_id or None
                association["associationId"] = eip.association_id or None
        return association

    @staticmethod
    def cloudformation_name_type() -> str:
        return ""

    @staticmethod
    def cloudformation_type() -> str:
        # https://docs.aws.amazon.com/AWSCloudFormation/latest/UserGuide/aws-resource-ec2-networkinterface.html
        return "AWS::EC2::NetworkInterface"

    @classmethod
    def create_from_cloudformation_json(  # type: ignore[misc]
        cls,
        resource_name: str,
        cloudformation_json: Any,
        account_id: str,
        region_name: str,
        **kwargs: Any,
    ) -> "NetworkInterface":
        from ..models import ec2_backends

        properties = cloudformation_json["Properties"]

        security_group_ids = properties.get("SecurityGroups", [])

        ec2_backend = ec2_backends[account_id][region_name]
        subnet_id = properties.get("SubnetId")
        if subnet_id:
            subnet = ec2_backend.get_subnet(subnet_id)
        else:
            subnet = None

        private_ip_address = properties.get("PrivateIpAddress", None)
        description = properties.get("Description", None)

        network_interface = ec2_backend.create_network_interface(
            subnet,
            private_ip_address,
            group_ids=security_group_ids,
            description=description,
        )
        return network_interface

    def stop(self) -> None:
        if self.public_ip_auto_assign:
            self.public_ip = None

    def start(self) -> None:
        self.check_auto_public_ip()

    def check_auto_public_ip(self) -> None:
        if (
            self.public_ip_auto_assign
            and str(self.public_ip_auto_assign).lower() == "true"
        ):
            self.public_ip = random_public_ip()

    @property
    def group_set(self) -> Any:  # type: ignore[misc]
        if self.instance and self.instance.security_groups:
            return set(self._group_set) | set(self.instance.security_groups)
        else:
            return self._group_set

    @classmethod
    def has_cfn_attr(cls, attr: str) -> bool:
        return attr in ["PrimaryPrivateIpAddress", "SecondaryPrivateIpAddresses"]

    def get_cfn_attribute(self, attribute_name: str) -> Any:
        from moto.cloudformation.exceptions import UnformattedGetAttTemplateException

        if attribute_name == "PrimaryPrivateIpAddress":
            return self.private_ip_address
        elif attribute_name == "SecondaryPrivateIpAddresses":
            raise NotImplementedError(
                '"Fn::GetAtt" : [ "{0}" , "SecondaryPrivateIpAddresses" ]"'
            )
        raise UnformattedGetAttTemplateException()

    @property
    def physical_resource_id(self) -> str:
        return self.id

    def get_filter_value(
        self, filter_name: str, method_name: Optional[str] = None
    ) -> Any:
        if filter_name == "network-interface-id":
            return self.id
        elif filter_name in ("addresses.private-ip-address", "private-ip-address"):
            return self.private_ip_address
        elif filter_name == "subnet-id":
            return self.subnet.id
        elif filter_name == "vpc-id":
            return self.subnet.vpc_id
        elif filter_name == "group-id":
            return [group.id for group in self._group_set]
        elif filter_name == "availability-zone":
            return self.subnet.availability_zone
        elif filter_name == "description":
            return self.description
        elif filter_name == "attachment.instance-id":
            return self.instance.id if self.instance else None
        elif filter_name == "attachment.instance-owner-id":
            return self.owner_id
        else:
            return super().get_filter_value(filter_name, "DescribeNetworkInterfaces")


class NetworkInterfaceBackend:
    def __init__(self) -> None:
        self.enis: Dict[str, NetworkInterface] = {}

    def create_network_interface(
        self,
        subnet: Any,
        private_ip_address: Union[str, List[str]],
        private_ip_addresses: Optional[List[Dict[str, Any]]] = None,
        group_ids: Optional[List[str]] = None,
        description: Optional[str] = None,
        tags: Optional[Dict[str, str]] = None,
        **kwargs: Any,
    ) -> NetworkInterface:
        eni = NetworkInterface(
            self,
            subnet,
            private_ip_address,
            private_ip_addresses,
            group_ids=group_ids,
            description=description,
            tags=tags,
            **kwargs,
        )
        self.enis[eni.id] = eni
        return eni

    def get_network_interface(self, eni_id: str) -> NetworkInterface:
        for eni in self.enis.values():
            if eni_id == eni.id:
                return eni
        raise InvalidNetworkInterfaceIdError(eni_id)

    def delete_network_interface(self, eni_id: str) -> None:
        deleted = self.enis.pop(eni_id, None)
        if not deleted:
            raise InvalidNetworkInterfaceIdError(eni_id)

    def describe_network_interfaces(
        self, filters: Any = None
    ) -> List[NetworkInterface]:
        # Note: This is only used in EC2Backend#do_resources_exist
        # Client-calls use #get_all_network_interfaces()
        # We should probably merge these at some point..
        enis = list(self.enis.values())

        if filters:
            for _filter, _filter_value in filters.items():
                if _filter == "network-interface-id":
                    _filter = "id"
                    enis = [
                        eni for eni in enis if getattr(eni, _filter) in _filter_value
                    ]
                else:
                    self.raise_not_implemented_error(  # type: ignore
                        f"The filter '{_filter}' for DescribeNetworkInterfaces"
                    )
        return enis

    def attach_network_interface(
        self, eni_id: str, instance_id: str, device_index: int
    ) -> str:
        eni = self.get_network_interface(eni_id)
        instance = self.get_instance(instance_id)  # type: ignore[attr-defined]
        return instance.attach_eni(eni, device_index)

    def detach_network_interface(self, attachment_id: str) -> None:
        for eni in self.enis.values():
            if eni.attachment_id == attachment_id:
                if eni.instance and len(eni.instance.nics) == 1:
                    raise LastEniDetachError
                eni.instance.detach_eni(eni)  # type: ignore
                return
        raise InvalidNetworkAttachmentIdError(attachment_id)

    def modify_network_interface_attribute(
        self,
        eni_id: str,
        group_ids: List[str],
        source_dest_check: Optional[bool] = None,
        description: Optional[str] = None,
    ) -> None:
        eni = self.get_network_interface(eni_id)
        groups = [self.get_security_group_from_id(group_id) for group_id in group_ids]  # type: ignore[attr-defined]
        if groups:
            eni._group_set = groups
        if source_dest_check in [True, False]:
            eni.source_dest_check = source_dest_check

        if description:
            eni.description = description

    def get_all_network_interfaces(
        self, eni_ids: Optional[List[str]] = None, filters: Any = None
    ) -> List[NetworkInterface]:
        enis = list(self.enis.values())

        if eni_ids:
            enis = [eni for eni in enis if eni.id in eni_ids]
            if len(enis) != len(eni_ids):
                invalid_id = list(
                    set(eni_ids).difference(set([eni.id for eni in enis]))
                )[0]
                raise InvalidNetworkInterfaceIdError(invalid_id)

        return generic_filter(filters, enis)

    def unassign_private_ip_addresses(
        self, eni_id: str, private_ip_address: Optional[List[str]] = None
    ) -> NetworkInterface:
        eni = self.get_network_interface(eni_id)
        if private_ip_address:
            for item in eni.private_ip_addresses.copy():
                if item.get("PrivateIpAddress") in private_ip_address:
                    eni.private_ip_addresses.remove(item)
        return eni

    def assign_private_ip_addresses(
        self,
        eni_id: str,
        private_ip_addresses: Optional[List[str]] = None,
        secondary_ips_count: Optional[int] = None,
    ) -> NetworkInterface:
        eni = self.get_network_interface(eni_id)
        eni_assigned_ips = [
            item.get("PrivateIpAddress") for item in eni.private_ip_addresses
        ]
        if private_ip_addresses:
            eni.private_ip_addresses.extend(
                {"Primary": False, "PrivateIpAddress": ip}
                for ip in private_ip_addresses
                if ip not in eni_assigned_ips
            )
            return eni
        while secondary_ips_count:
            ip = random_private_ip(eni.subnet.cidr_block)
            if ip not in eni_assigned_ips:
                eni.private_ip_addresses.append(
                    {"Primary": False, "PrivateIpAddress": ip}
                )
                secondary_ips_count -= 1
        return eni

    def assign_ipv6_addresses(
        self,
        eni_id: str,
        ipv6_addresses: Optional[List[str]] = None,
        ipv6_count: Optional[int] = None,
    ) -> NetworkInterface:
        eni = self.get_network_interface(eni_id)
        if ipv6_addresses:
            eni.ipv6_addresses.extend(
                (ip for ip in ipv6_addresses if ip not in eni.ipv6_addresses)
            )

        while ipv6_count:
            association = list(eni.subnet.ipv6_cidr_block_associations.values())[0]
            subnet_ipv6_cidr_block = association.get("ipv6CidrBlock")
            ip = random_private_ip(subnet_ipv6_cidr_block, ipv6=True)
            if ip not in eni.ipv6_addresses:
                eni.ipv6_addresses.append(ip)
                ipv6_count -= 1
        return eni

    def unassign_ipv6_addresses(
        self, eni_id: str, ips: Optional[List[str]] = None
    ) -> NetworkInterface:
        eni = self.get_network_interface(eni_id)
        if ips:
            for ip in eni.ipv6_addresses.copy():
                if ip in ips:
                    eni.ipv6_addresses.remove(ip)
        return eni
