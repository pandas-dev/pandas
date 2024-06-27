import ipaddress
import itertools
from collections import defaultdict
from typing import TYPE_CHECKING, Any, Dict, List, Optional, Set, Tuple

from moto.core.common_models import CloudFormationModel

if TYPE_CHECKING:
    from moto.ec2.models.instances import Instance

from moto.ec2.models.availability_zones_and_regions import Zone

from ..exceptions import (
    GenericInvalidParameterValueError,
    InvalidAvailabilityZoneError,
    InvalidCIDRBlockParameterError,
    InvalidParameterValueError,
    InvalidSubnetCidrBlockAssociationID,
    InvalidSubnetConflictError,
    InvalidSubnetIdError,
    InvalidSubnetRangeError,
)
from ..utils import (
    generic_filter,
    random_subnet_id,
    random_subnet_ipv6_cidr_block_association_id,
)
from .availability_zones_and_regions import RegionsAndZonesBackend
from .core import TaggedEC2Resource


class Subnet(TaggedEC2Resource, CloudFormationModel):
    def __init__(
        self,
        ec2_backend: Any,
        subnet_id: str,
        vpc_id: str,
        cidr_block: str,
        ipv6_cidr_block: Optional[str],
        availability_zone: Zone,
        default_for_az: str,
        map_public_ip_on_launch: str,
        assign_ipv6_address_on_creation: bool = False,
    ):
        self.ec2_backend = ec2_backend
        self.id = subnet_id
        self.vpc_id = vpc_id
        self.cidr_block = cidr_block
        self.cidr = ipaddress.IPv4Network(str(self.cidr_block), strict=False)
        self._available_ip_addresses = self.cidr.num_addresses - 5
        self._availability_zone = availability_zone
        self.default_for_az = default_for_az
        self.map_public_ip_on_launch = map_public_ip_on_launch
        self.assign_ipv6_address_on_creation = assign_ipv6_address_on_creation
        self.ipv6_cidr_block_associations: Dict[str, Dict[str, Any]] = {}
        if ipv6_cidr_block:
            self.attach_ipv6_cidr_block_associations(ipv6_cidr_block)

        # Theory is we assign ip's as we go (as 16,777,214 usable IPs in a /8)
        self._subnet_ip_generator = self.cidr.hosts()
        self.reserved_ips = [
            next(self._subnet_ip_generator) for _ in range(0, 3)
        ]  # Reserved by AWS
        self._unused_ips: Set[str] = (
            set()
        )  # if instance is destroyed hold IP here for reuse
        self._subnet_ips: Dict[str, "Instance"] = {}
        self.state = "available"

        # Placeholder for response templates until Ipv6 support implemented.
        self.ipv6_native = False

    @property
    def owner_id(self) -> str:
        return self.ec2_backend.account_id

    @staticmethod
    def cloudformation_name_type() -> str:
        return ""

    @staticmethod
    def cloudformation_type() -> str:
        # https://docs.aws.amazon.com/AWSCloudFormation/latest/UserGuide/aws-resource-ec2-subnet.html
        return "AWS::EC2::Subnet"

    @classmethod
    def create_from_cloudformation_json(  # type: ignore[misc]
        cls,
        resource_name: str,
        cloudformation_json: Any,
        account_id: str,
        region_name: str,
        **kwargs: Any,
    ) -> "Subnet":
        from ..models import ec2_backends

        properties = cloudformation_json["Properties"]

        vpc_id = properties["VpcId"]
        cidr_block = properties["CidrBlock"]
        availability_zone = properties.get("AvailabilityZone")
        ec2_backend = ec2_backends[account_id][region_name]
        subnet = ec2_backend.create_subnet(
            vpc_id=vpc_id, cidr_block=cidr_block, availability_zone=availability_zone
        )
        for tag in properties.get("Tags", []):
            tag_key = tag["Key"]
            tag_value = tag["Value"]
            subnet.add_tag(tag_key, tag_value)

        return subnet

    @classmethod
    def delete_from_cloudformation_json(  # type: ignore[misc]
        cls,
        resource_name: str,
        cloudformation_json: Dict[str, Any],
        account_id: str,
        region_name: str,
    ) -> None:
        from ..models import ec2_backends

        ec2_backends[account_id][region_name].delete_subnet(resource_name)

    @property
    def available_ip_addresses(self) -> str:
        enis = [
            eni
            for eni in self.ec2_backend.get_all_network_interfaces()
            if eni.subnet.id == self.id
        ]
        addresses_taken = []
        for eni in enis:
            if eni.private_ip_addresses:
                addresses_taken.extend(eni.private_ip_addresses)
        return str(self._available_ip_addresses - len(addresses_taken))

    @property
    def availability_zone(self) -> str:
        return self._availability_zone.name

    @property
    def availability_zone_id(self) -> str:
        return self._availability_zone.zone_id

    @property
    def physical_resource_id(self) -> str:
        return self.id

    def get_filter_value(
        self, filter_name: str, method_name: Optional[str] = None
    ) -> Any:
        """
        API Version 2014-10-01 defines the following filters for DescribeSubnets:

        * availabilityZone
        * available-ip-address-count
        * cidrBlock
        * defaultForAz
        * state
        * subnet-id
        * tag:key=value
        * tag-key
        * tag-value
        * vpc-id

        Taken from: http://docs.aws.amazon.com/AWSEC2/latest/APIReference/API_DescribeSubnets.html
        """
        if filter_name in ("cidr", "cidrBlock", "cidr-block"):
            return self.cidr_block
        elif filter_name in ("vpc-id", "vpcId"):
            return self.vpc_id
        elif filter_name == "subnet-id":
            return self.id
        elif filter_name in ("availabilityZone", "availability-zone"):
            return self.availability_zone
        elif filter_name in ("defaultForAz", "default-for-az"):
            return self.default_for_az
        elif filter_name == "state":
            return self.state
        else:
            return super().get_filter_value(filter_name, "DescribeSubnets")

    @classmethod
    def has_cfn_attr(cls, attr: str) -> bool:
        return attr in ["AvailabilityZone"]

    def get_cfn_attribute(self, attribute_name: str) -> None:
        from moto.cloudformation.exceptions import UnformattedGetAttTemplateException

        if attribute_name == "AvailabilityZone":
            raise NotImplementedError('"Fn::GetAtt" : [ "{0}" , "AvailabilityZone" ]"')
        raise UnformattedGetAttTemplateException()

    def get_available_subnet_ip(self, instance: "Instance") -> str:
        try:
            new_ip = str(self._unused_ips.pop())
        except KeyError:
            new_ip_v4 = next(self._subnet_ip_generator)

            # Skips any IP's if they've been manually specified
            while str(new_ip_v4) in self._subnet_ips:
                new_ip_v4 = next(self._subnet_ip_generator)

            if new_ip_v4 == self.cidr.broadcast_address:
                raise StopIteration()  # Broadcast address cant be used obviously
            new_ip = str(new_ip_v4)
        # TODO StopIteration will be raised if no ip's available, not sure how aws handles this.

        self._subnet_ips[new_ip] = instance

        return new_ip

    # EFS calls this method as request_ip(str, MountTarget)
    # So technically it's not just Instances that are stored
    def request_ip(self, ip: str, instance: "Instance") -> None:
        if ipaddress.ip_address(ip) not in self.cidr:
            raise Exception(f"IP does not fall in the subnet CIDR of {self.cidr}")

        if ip in self._subnet_ips:
            raise Exception("IP already in use")
        try:
            self._unused_ips.remove(ip)
        except KeyError:
            pass

        self._subnet_ips[ip] = instance

    def del_subnet_ip(self, ip: str) -> None:
        try:
            del self._subnet_ips[ip]
            self._unused_ips.add(ip)
        except KeyError:
            pass  # Unknown IP

    def attach_ipv6_cidr_block_associations(
        self, ipv6_cidr_block: str
    ) -> Dict[str, str]:
        association = {
            "associationId": random_subnet_ipv6_cidr_block_association_id(),
            "ipv6CidrBlock": ipv6_cidr_block,
            "ipv6CidrBlockState": "associated",
        }
        self.ipv6_cidr_block_associations[association["associationId"]] = association
        return association

    def detach_subnet_cidr_block(self, association_id: str) -> Dict[str, Any]:
        association = self.ipv6_cidr_block_associations.get(association_id)
        association["ipv6CidrBlockState"] = "disassociated"  # type: ignore[index]
        return association  # type: ignore[return-value]


class SubnetRouteTableAssociation(CloudFormationModel):
    def __init__(self, route_table_id: str, subnet_id: str):
        self.route_table_id = route_table_id
        self.subnet_id = subnet_id

    @property
    def physical_resource_id(self) -> str:
        return self.route_table_id

    @staticmethod
    def cloudformation_name_type() -> str:
        return ""

    @staticmethod
    def cloudformation_type() -> str:
        # https://docs.aws.amazon.com/AWSCloudFormation/latest/UserGuide/aws-resource-ec2-subnetroutetableassociation.html
        return "AWS::EC2::SubnetRouteTableAssociation"

    @classmethod
    def create_from_cloudformation_json(  # type: ignore[misc]
        cls,
        resource_name: str,
        cloudformation_json: Any,
        account_id: str,
        region_name: str,
        **kwargs: Any,
    ) -> "SubnetRouteTableAssociation":
        from ..models import ec2_backends

        properties = cloudformation_json["Properties"]

        route_table_id = properties["RouteTableId"]
        subnet_id = properties["SubnetId"]

        ec2_backend = ec2_backends[account_id][region_name]
        subnet_association = ec2_backend.create_subnet_association(
            route_table_id=route_table_id, subnet_id=subnet_id
        )
        return subnet_association


class SubnetBackend:
    def __init__(self) -> None:
        # maps availability zone to dict of (subnet_id, subnet)
        self.subnets: Dict[str, Dict[str, Subnet]] = defaultdict(dict)
        self.subnet_associations: Dict[str, SubnetRouteTableAssociation] = {}

    def get_subnet(self, subnet_id: str) -> Subnet:
        for subnets_per_zone in self.subnets.values():
            if subnet_id in subnets_per_zone:
                return subnets_per_zone[subnet_id]
        raise InvalidSubnetIdError(subnet_id)

    def get_default_subnet(self, availability_zone: str) -> Subnet:
        return [
            subnet
            for subnet in self.describe_subnets(
                filters={"availabilityZone": availability_zone}
            )
            if subnet.default_for_az
        ][0]

    def create_subnet(
        self,
        vpc_id: str,
        cidr_block: str,
        ipv6_cidr_block: Optional[str] = None,
        availability_zone: Optional[str] = None,
        availability_zone_id: Optional[str] = None,
        tags: Optional[Dict[str, Dict[str, str]]] = None,
    ) -> Subnet:
        subnet_id = random_subnet_id()
        # Validate VPC exists and the supplied CIDR block is a subnet of the VPC's
        vpc = self.get_vpc(vpc_id)  # type: ignore[attr-defined]
        vpc_cidr_blocks = [
            ipaddress.IPv4Network(
                str(cidr_block_association["cidr_block"]), strict=False
            )
            for cidr_block_association in vpc.get_cidr_block_association_set()
        ]
        try:
            subnet_cidr_block = ipaddress.IPv4Network(str(cidr_block), strict=False)
        except ValueError:
            raise InvalidCIDRBlockParameterError(cidr_block)

        subnet_in_vpc_cidr_range = False
        for vpc_cidr_block in vpc_cidr_blocks:
            if (
                vpc_cidr_block.network_address <= subnet_cidr_block.network_address
                and vpc_cidr_block.broadcast_address
                >= subnet_cidr_block.broadcast_address
            ):
                subnet_in_vpc_cidr_range = True
                break

        if not subnet_in_vpc_cidr_range:
            raise InvalidSubnetRangeError(cidr_block)

        # The subnet size must use a /64 prefix length.
        if ipv6_cidr_block and "::/64" not in ipv6_cidr_block:
            raise GenericInvalidParameterValueError("ipv6-cidr-block", ipv6_cidr_block)

        for subnet in self.describe_subnets(filters={"vpc-id": vpc_id}):
            if subnet.cidr.overlaps(subnet_cidr_block):
                raise InvalidSubnetConflictError(cidr_block)

        # if this is the first subnet for an availability zone,
        # consider it the default
        default_for_az = str(availability_zone not in self.subnets).lower()
        map_public_ip_on_launch = default_for_az

        if availability_zone is None and not availability_zone_id:
            availability_zone = "us-east-1a"
        try:
            if availability_zone:
                availability_zone_data = next(
                    zone
                    for zones in RegionsAndZonesBackend.zones.values()
                    for zone in zones
                    if zone.name == availability_zone
                )
            elif availability_zone_id:
                availability_zone_data = next(
                    zone
                    for zones in RegionsAndZonesBackend.zones.values()
                    for zone in zones
                    if zone.zone_id == availability_zone_id
                )

        except StopIteration:
            raise InvalidAvailabilityZoneError(
                availability_zone,
                ", ".join(
                    [
                        zone.name
                        for zones in RegionsAndZonesBackend.zones.values()
                        for zone in zones
                    ]
                ),
            )
        subnet = Subnet(
            self,
            subnet_id,
            vpc_id,
            cidr_block,
            ipv6_cidr_block,
            availability_zone_data,
            default_for_az,
            map_public_ip_on_launch,
            assign_ipv6_address_on_creation=False,
        )

        for k, v in tags.get("subnet", {}).items() if tags else []:
            subnet.add_tag(k, v)

        # AWS associates a new subnet with the default Network ACL
        self.associate_default_network_acl_with_subnet(subnet_id, vpc_id)  # type: ignore[attr-defined]
        self.subnets[availability_zone][subnet_id] = subnet  # type: ignore[index]
        return subnet

    def describe_subnets(
        self, subnet_ids: Optional[List[str]] = None, filters: Optional[Any] = None
    ) -> List[Subnet]:
        # Extract a list of all subnets
        matches = list(
            itertools.chain(*[x.copy().values() for x in self.subnets.copy().values()])
        )
        if subnet_ids:
            matches = [sn for sn in matches if sn.id in subnet_ids]
            if len(subnet_ids) > len(matches):
                unknown_ids = set(subnet_ids) - set(matches)  # type: ignore[arg-type]
                raise InvalidSubnetIdError(list(unknown_ids)[0])
        if filters:
            matches = generic_filter(filters, matches)

        return matches

    def delete_subnet(self, subnet_id: str) -> Subnet:
        for subnets in self.subnets.values():
            if subnet_id in subnets:
                return subnets.pop(subnet_id)
        raise InvalidSubnetIdError(subnet_id)

    def modify_subnet_attribute(
        self, subnet_id: str, attr_name: str, attr_value: str
    ) -> None:
        subnet = self.get_subnet(subnet_id)
        if attr_name in ("map_public_ip_on_launch", "assign_ipv6_address_on_creation"):
            setattr(subnet, attr_name, attr_value)
        else:
            raise InvalidParameterValueError(attr_name)

    def get_subnet_from_ipv6_association(self, association_id: str) -> Optional[Subnet]:
        subnet = None
        for s in self.describe_subnets():
            if association_id in s.ipv6_cidr_block_associations:
                subnet = s
        return subnet

    def associate_subnet_cidr_block(
        self, subnet_id: str, ipv6_cidr_block: str
    ) -> Dict[str, str]:
        subnet = self.get_subnet(subnet_id)
        if not subnet:
            raise InvalidSubnetIdError(subnet_id)
        association = subnet.attach_ipv6_cidr_block_associations(ipv6_cidr_block)
        return association

    def disassociate_subnet_cidr_block(
        self, association_id: str
    ) -> Tuple[str, Dict[str, str]]:
        subnet = self.get_subnet_from_ipv6_association(association_id)
        if not subnet:
            raise InvalidSubnetCidrBlockAssociationID(association_id)
        association = subnet.detach_subnet_cidr_block(association_id)
        return subnet.id, association

    def create_subnet_association(
        self, route_table_id: str, subnet_id: str
    ) -> SubnetRouteTableAssociation:
        subnet_association = SubnetRouteTableAssociation(route_table_id, subnet_id)
        self.subnet_associations[f"{route_table_id}:{subnet_id}"] = subnet_association
        return subnet_association
