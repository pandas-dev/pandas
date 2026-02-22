from __future__ import annotations

import ipaddress
import itertools
from collections import defaultdict
from collections.abc import Iterator
from functools import cache
from typing import TYPE_CHECKING, Any, Literal, Optional, cast

from moto.core.common_models import CloudFormationModel

if TYPE_CHECKING:
    from moto.ec2.models.instances import Instance

from moto.ec2.models.availability_zones_and_regions import Zone

from ..exceptions import (
    DefaultSubnetAlreadyExistsInAvailabilityZoneError,
    DefaultVpcDoesNotExistError,
    InvalidAvailabilityZoneError,
    InvalidCIDRBlockParameterError,
    InvalidCidrReservationNotFound,
    InvalidCidrReservationNotWithinSubnetCidr,
    InvalidCidrReservationOverlapExisting,
    InvalidParameterValueError,
    InvalidPrefixReservationType,
    InvalidSubnetCidrBlockAssociationID,
    InvalidSubnetConflictError,
    InvalidSubnetIdError,
    InvalidSubnetRangeError,
)
from ..utils import (
    generic_filter,
    random_subnet_cidr_reservation_id,
    random_subnet_id,
    random_subnet_ipv6_cidr_block_association_id,
)
from .availability_zones_and_regions import RegionsAndZonesBackend
from .core import TaggedEC2Resource


class SubnetCidrReservation(TaggedEC2Resource):
    def __init__(
        self,
        ec2_backend: Any,
        subnet_id: str,
        cidr: ipaddress.IPv4Network | ipaddress.IPv6Network,
        reservation_type: Literal["explicit"] | Literal["prefix"],
        owner_id: str,
        tag_specifications: Optional[list[dict[str, Any]]] = None,
        subnet_cidr_reservation_id: Optional[str] = None,
    ) -> None:
        self.ec2_backend = ec2_backend
        self.id = subnet_cidr_reservation_id or random_subnet_cidr_reservation_id()
        self.subnet_id = subnet_id
        self.cidr = cidr
        self.reservation_type = reservation_type
        self.owner_id = owner_id
        self.tag_specifications = tag_specifications or []

    def get_filter_value(
        self, filter_name: str, method_name: Optional[str] = None
    ) -> Any:
        if filter_name == "reservationType":
            return self.reservation_type
        elif filter_name == "subnet-id":
            return self.subnet_id

        return super().get_filter_value(filter_name, method_name)


class Subnet(TaggedEC2Resource, CloudFormationModel):
    def __init__(
        self,
        ec2_backend: Any,
        subnet_id: str,
        vpc_id: str,
        cidr_block: str,
        ipv6_cidr_block: Optional[str],
        availability_zone: Zone,
        default_for_az: bool,
        map_public_ip_on_launch: bool,
        ipv6_native: bool = False,
    ):
        self.ec2_backend = ec2_backend
        self.id = subnet_id
        self.vpc_id = vpc_id
        self.cidr_block = cidr_block
        self.cidr = (
            ipaddress.IPv4Network(str(self.cidr_block), strict=False)
            if self.cidr_block
            else ipaddress.IPv6Network(ipv6_cidr_block)
        )
        self._available_ip_addresses = self.cidr.num_addresses - 5
        self._availability_zone = availability_zone
        self.default_for_az = default_for_az
        self.map_public_ip_on_launch = map_public_ip_on_launch
        self.assign_ipv6_address_on_creation = ipv6_native
        self.ipv6_cidr_block_associations: dict[str, dict[str, Any]] = {}
        if ipv6_cidr_block:
            self.attach_ipv6_cidr_block_associations(ipv6_cidr_block)

        self._reserved_ips: set[ipaddress.IPv4Address | ipaddress.IPv6Address] = set()

        # Reserved by AWS
        self.aws_reserved_ips = {
            next(iter(self._subnet_ip_generator)) for _ in range(0, 3)
        }

        self._unused_ips: set[str] = (
            set()
        )  # if instance is destroyed hold IP here for reuse
        self._subnet_ips: dict[str, Instance] = {}
        self.state = "available"

        # Placeholder for response templates until Ipv6 support implemented.
        self.ipv6_native = ipv6_native

        self._cidr_reservations: list[SubnetCidrReservation] = []

    @property
    def arn(self) -> str:
        return f"arn:aws:ec2:{self.ec2_backend.region_name}:{self.owner_id}:subnet/{self.id}"

    @property
    def ipv6_cidr_block_association_set(self) -> list[dict[str, str]]:
        association_set = [
            {
                "ipv6CidrBlock": association["ipv6CidrBlock"],
                "associationId": association["associationId"],
                "ipv6CidrBlockState": {
                    "state": association["ipv6CidrBlockState"],
                },
            }
            for association in self.ipv6_cidr_block_associations.values()
            if association["ipv6CidrBlockState"]["State"] == "associated"
        ]
        return association_set

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
    ) -> Subnet:
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
        cloudformation_json: dict[str, Any],
        account_id: str,
        region_name: str,
    ) -> None:
        from ..models import ec2_backends

        ec2_backends[account_id][region_name].delete_subnet(resource_name)

    @property
    def available_ip_address_count(self) -> int:
        enis = [
            eni
            for eni in self.ec2_backend.get_all_network_interfaces()
            if eni.subnet.id == self.id
        ]
        addresses_taken = []
        for eni in enis:
            if eni.private_ip_addresses:
                addresses_taken.extend(eni.private_ip_addresses)
        return self._available_ip_addresses - len(addresses_taken)

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

    def get_available_subnet_ip(self, instance: Instance) -> str:
        try:
            new_ip = str(self._unused_ips.pop())
        except KeyError:
            new_ip_v4 = next(iter(self._subnet_ip_generator))

            # Skips any IP's if they've been manually specified
            while str(new_ip_v4) in self._subnet_ips:
                new_ip_v4 = next(iter(self._subnet_ip_generator))

            if new_ip_v4 == self.cidr.broadcast_address:
                raise StopIteration()  # Broadcast address cant be used obviously
            new_ip = str(new_ip_v4)
        # TODO StopIteration will be raised if no ip's available, not sure how aws handles this.

        self._subnet_ips[new_ip] = instance

        return new_ip

    # EFS calls this method as request_ip(str, MountTarget)
    # So technically it's not just Instances that are stored
    def request_ip(self, ip: str, instance: Instance) -> None:
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
    ) -> dict[str, Any]:
        association = {
            "associationId": random_subnet_ipv6_cidr_block_association_id(),
            "ipv6CidrBlock": ipv6_cidr_block,
            "ipv6CidrBlockState": {"State": "associated"},
        }
        self.ipv6_cidr_block_associations[str(association["associationId"])] = (
            association
        )
        return association

    def detach_subnet_cidr_block(self, association_id: str) -> dict[str, Any]:
        association = self.ipv6_cidr_block_associations.get(association_id)
        assert association is not None
        association["ipv6CidrBlockState"] = {"State": "disassociated"}
        return association

    def create_cidr_reservation(
        self,
        reservation_type: Literal["explicit"] | Literal["prefix"],
        cidr_block: ipaddress.IPv4Network | ipaddress.IPv6Network,
        tag_specifications: Optional[list[dict[str, Any]]] = None,
    ) -> SubnetCidrReservation:
        if not self._cidr_within_subnet_cidr_blocks(cidr_block):
            raise InvalidCidrReservationNotWithinSubnetCidr(str(self.cidr))

        if self._cidr_overlaps_with_existing_reservation(cidr_block):
            raise InvalidCidrReservationOverlapExisting(str(cidr_block))

        reservation = SubnetCidrReservation(
            ec2_backend=self.ec2_backend,
            subnet_id=self.id,
            cidr=cidr_block,
            reservation_type=reservation_type,
            owner_id=self.owner_id,
            tag_specifications=tag_specifications,
            subnet_cidr_reservation_id=random_subnet_cidr_reservation_id(),
        )

        for ip in cidr_block:
            self._reserved_ips.add(ip)

        if reservation_type == "prefix":
            # prefix reservations impact the avail ip count immediately
            self._available_ip_addresses -= cidr_block.num_addresses

        self._cidr_reservations.append(reservation)

        return reservation

    @cache
    def _base_avail_ip_iterator(
        self,
    ) -> (
        Iterator[ipaddress.IPv4Address]
        | list[ipaddress.IPv4Address]
        | Iterator[ipaddress.IPv6Address]
        | list[ipaddress.IPv6Address]
    ):
        return self.cidr.hosts()

    @property
    def _subnet_ip_generator(
        self,
    ) -> Iterator[ipaddress.IPv4Address | ipaddress.IPv6Address]:
        for host in self._base_avail_ip_iterator():
            # Skip any reserved ips
            if host in self._reserved_ips:
                continue
            yield host

    def _cidr_within_subnet_cidr_blocks(
        self, cidr_to_reserve: ipaddress.IPv4Network | ipaddress.IPv6Network
    ) -> bool:
        if isinstance(cidr_to_reserve, ipaddress.IPv4Network):
            if isinstance(self.cidr, ipaddress.IPv4Network):
                return cidr_to_reserve.subnet_of(self.cidr)
            return False

        for _, ipv6_block_association in self.ipv6_cidr_block_associations.items():
            ipv6_block = ipaddress.ip_network(ipv6_block_association["ipv6CidrBlock"])
            if isinstance(
                ipv6_block, ipaddress.IPv6Network
            ) and cidr_to_reserve.subnet_of(ipv6_block):
                return True
        return False

    def _cidr_overlaps_with_existing_reservation(
        self, cidr_to_reserve: ipaddress.IPv4Network | ipaddress.IPv6Network
    ) -> bool:
        cidr_to_reserve_is_ipv4: bool = isinstance(
            cidr_to_reserve, ipaddress.IPv4Network
        )
        for reservation in self._cidr_reservations:
            existing_block = reservation.cidr
            if (
                isinstance(existing_block, ipaddress.IPv4Network)
                and cidr_to_reserve_is_ipv4
                and cidr_to_reserve.overlaps(existing_block)
            ):
                return True
            elif (
                isinstance(existing_block, ipaddress.IPv6Network)
                and not cidr_to_reserve_is_ipv4
                and cidr_to_reserve.overlaps(existing_block)
            ):
                return True

        return False

    def delete_subnet_cidr_reservation(
        self,
        reservation_id: str,
    ) -> Optional[SubnetCidrReservation]:
        for index, reservation in enumerate(self._cidr_reservations):
            if reservation.id == reservation_id:
                self._cidr_reservations.pop(index)

                if reservation.reservation_type == "prefix":
                    for ip in reservation.cidr:
                        if (
                            ip not in self._subnet_ips
                            and ip not in self.aws_reserved_ips
                        ):
                            self._available_ip_addresses += 1
                return reservation

        return None


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
    ) -> SubnetRouteTableAssociation:
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
        self.subnets: dict[str, dict[str, Subnet]] = defaultdict(dict)
        self.subnet_associations: dict[str, SubnetRouteTableAssociation] = {}

    def get_subnet(self, subnet_id: str) -> Subnet:
        for subnets_per_zone in self.subnets.values():
            if subnet_id in subnets_per_zone:
                return subnets_per_zone[subnet_id]
        raise InvalidSubnetIdError(subnet_id)

    def get_default_subnets(self) -> dict[str, Subnet]:
        return {
            subnet.availability_zone: subnet
            for subnet in self.describe_subnets()
            if subnet.default_for_az
        }

    def create_default_subnet(self, availability_zone: str) -> Subnet:
        default_subnets = self.get_default_subnets()
        if availability_zone in default_subnets:
            raise DefaultSubnetAlreadyExistsInAvailabilityZoneError(
                default_subnets[availability_zone].id, availability_zone
            )

        default_vpc = self.get_default_vpc()  # type: ignore[attr-defined]
        if default_vpc is None:
            raise DefaultVpcDoesNotExistError()

        cidr_block = default_vpc.cidr_block
        blocks = list(ipaddress.ip_network(cidr_block).subnets(new_prefix=20))
        subnet = self.create_subnet(
            vpc_id=default_vpc.id,
            cidr_block=str(blocks[len(default_subnets)]),
            availability_zone=availability_zone,
        )
        return subnet

    def create_subnet(
        self,
        vpc_id: str,
        cidr_block: str,
        ipv6_cidr_block: Optional[str] = None,
        availability_zone: Optional[str] = None,
        availability_zone_id: Optional[str] = None,
        ipv6_native: bool = False,
        tags: Optional[dict[str, dict[str, str]]] = None,
    ) -> Subnet:
        subnet_id = random_subnet_id()
        # Validate VPC exists and the supplied CIDR block is a subnet of the VPC's
        vpc = self.get_vpc(vpc_id)  # type: ignore[attr-defined]
        if cidr_block:
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

            for subnet in self.describe_subnets(filters={"vpc-id": vpc_id}):
                if subnet.cidr.overlaps(subnet_cidr_block):
                    raise InvalidSubnetConflictError(cidr_block)

        # if this is the first subnet for an availability zone,
        # consider it the default
        default_for_az = len(self.subnets.get(availability_zone, [])) == 0  # type: ignore[arg-type]
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
            ipv6_native=ipv6_native,
        )

        for k, v in tags.get("subnet", {}).items() if tags else []:
            subnet.add_tag(k, v)

        # AWS associates a new subnet with the default Network ACL
        self.associate_default_network_acl_with_subnet(subnet_id, vpc_id)  # type: ignore[attr-defined]
        self.subnets[availability_zone][subnet_id] = subnet  # type: ignore[index]
        return subnet

    def describe_subnets(
        self, subnet_ids: Optional[list[str]] = None, filters: Optional[Any] = None
    ) -> list[Subnet]:
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
    ) -> dict[str, Any]:
        subnet = self.get_subnet(subnet_id)
        if not subnet:
            raise InvalidSubnetIdError(subnet_id)
        association = subnet.attach_ipv6_cidr_block_associations(ipv6_cidr_block)
        return association

    def disassociate_subnet_cidr_block(
        self, association_id: str
    ) -> tuple[str, dict[str, str]]:
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

    def create_subnet_cidr_reservation(
        self,
        subnet_id: str,
        reservation_type: str,
        cidr: str,
        tags: Optional[list[dict[str, Any]]] = None,
    ) -> SubnetCidrReservation:
        if reservation_type not in {"prefix", "explicit"}:
            raise InvalidPrefixReservationType(reservation_type)

        reservation_type_lit: Literal["prefix", "explicit"] = cast(
            Literal["prefix", "explicit"], reservation_type
        )

        cidr_to_reserve: ipaddress.IPv4Network | ipaddress.IPv6Network
        try:
            cidr_to_reserve = ipaddress.ip_network(cidr)
        except ValueError:
            raise InvalidCIDRBlockParameterError(cidr)

        subnet = self.get_subnet(subnet_id)

        scr = subnet.create_cidr_reservation(
            reservation_type_lit, cidr_to_reserve, tags
        )
        return scr

    def get_subnet_cidr_reservations(
        self,
        subnet_id: str,
        filters: Optional[Any] = None,
    ) -> dict[str, list[SubnetCidrReservation]]:
        matches: list[SubnetCidrReservation] = self.get_subnet(
            subnet_id
        )._cidr_reservations

        if filters:
            matches = generic_filter(filters, matches)

        ipv4_reservations = []
        ipv6_reservations = []

        for reservation in matches:
            if isinstance(reservation.cidr, ipaddress.IPv4Network):
                ipv4_reservations.append(reservation)
            else:
                ipv6_reservations.append(reservation)

        return {
            "SubnetIpv4CidrReservations": ipv4_reservations,
            "SubnetIpv6CidrReservations": ipv6_reservations,
        }

    def delete_subnet_cidr_reservation(
        self,
        reservation_id: str,
    ) -> SubnetCidrReservation:
        for az_subnet in self.subnets.values():
            for subnet in az_subnet.values():
                deleted_reservation = subnet.delete_subnet_cidr_reservation(
                    reservation_id
                )
                if deleted_reservation:
                    return deleted_reservation

        raise InvalidCidrReservationNotFound(reservation_id)
