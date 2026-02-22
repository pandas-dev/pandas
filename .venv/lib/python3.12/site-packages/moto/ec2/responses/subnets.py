from moto.core.responses import ActionResult, EmptyResult
from moto.core.utils import camelcase_to_underscores
from moto.moto_api._internal import mock_random as random

from ..exceptions import InvalidParameterCombination, MissingParameter
from ._base_response import EC2BaseResponse


class Subnets(EC2BaseResponse):
    def create_subnet(self) -> ActionResult:
        vpc_id = self._get_param("VpcId")
        cidr_block = self._get_param("CidrBlock")
        ipv4_ipam_pool_id = self._get_param("Ipv4IpamPoolId")
        ipv6_cidr_block = self._get_param("Ipv6CidrBlock")
        ipv6_ipam_pool_id = self._get_param("Ipv6IpamPoolId")
        availability_zone = self._get_param("AvailabilityZone")
        availability_zone_id = self._get_param("AvailabilityZoneId")
        tags = self._parse_tag_specification("subnet")
        ipv6_native = self._get_bool_param("Ipv6Native", False)

        if ipv6_native and not ipv6_cidr_block and not ipv6_ipam_pool_id:
            raise MissingParameter(
                "Either 'ipv6CidrBlock' or 'ipv6IpamPoolId' should be provided."
            )

        if ipv6_native and (cidr_block or ipv4_ipam_pool_id):
            raise InvalidParameterCombination(
                "When specifying ipv4 parameters, cidrBlock or ipv4IpamPoolId, you cannot set ipv6Native to true."
            )

        if not availability_zone and not availability_zone_id:
            availability_zone = random.choice(
                self.ec2_backend.describe_availability_zones()
            ).name
        subnet = self.ec2_backend.create_subnet(
            vpc_id,
            cidr_block,
            ipv6_cidr_block,
            availability_zone,
            availability_zone_id,
            ipv6_native=ipv6_native,
            tags=tags,
        )
        result = {"Subnet": subnet}
        return ActionResult(result)

    def create_default_subnet(self) -> ActionResult:
        availability_zone = self._get_param("AvailabilityZone")

        subnet = self.ec2_backend.create_default_subnet(availability_zone)
        result = {"Subnet": subnet}
        return ActionResult(result)

    def delete_subnet(self) -> ActionResult:
        subnet_id = self._get_param("SubnetId")
        self.ec2_backend.delete_subnet(subnet_id)
        return EmptyResult()

    def describe_subnets(self) -> ActionResult:
        self.error_on_dryrun()
        subnet_ids = self._get_param("SubnetIds")
        filters = self._filters_from_querystring()
        subnets = self.ec2_backend.describe_subnets(subnet_ids, filters)
        result = {"Subnets": subnets}
        return ActionResult(result)

    def modify_subnet_attribute(self) -> ActionResult:
        subnet_id = self._get_param("SubnetId")

        for attribute in ("MapPublicIpOnLaunch", "AssignIpv6AddressOnCreation"):
            if self._get_param(f"{attribute}.Value") is not None:
                attr_name = camelcase_to_underscores(attribute)
                attr_value = self._get_param(f"{attribute}.Value")
                self.ec2_backend.modify_subnet_attribute(
                    subnet_id, attr_name, attr_value
                )
        return EmptyResult()

    def associate_subnet_cidr_block(self) -> ActionResult:
        ipv6_cidr_block = self._get_param("Ipv6CidrBlock")
        subnet_id = self._get_param("SubnetId")
        association = self.ec2_backend.associate_subnet_cidr_block(
            subnet_id, ipv6_cidr_block
        )
        result = {"Ipv6CidrBlockAssociation": association, "SubnetId": subnet_id}
        return ActionResult(result)

    def disassociate_subnet_cidr_block(self) -> ActionResult:
        association_id = self._get_param("AssociationId")
        subnet_id, association = self.ec2_backend.disassociate_subnet_cidr_block(
            association_id
        )
        result = {"Ipv6CidrBlockAssociation": association, "SubnetId": subnet_id}
        return ActionResult(result)

    def create_subnet_cidr_reservation(self) -> ActionResult:
        self.error_on_dryrun()
        subnet_id = self._get_param("SubnetId")
        reservation_type = self._get_param("ReservationType")
        cidr = self._get_param("Cidr")
        tag_specifications = self._get_param("TagSpecifications", [])

        reservation = self.ec2_backend.create_subnet_cidr_reservation(
            subnet_id, reservation_type, cidr, tag_specifications
        )

        return ActionResult({"SubnetCidrReservation": reservation})

    def get_subnet_cidr_reservations(self) -> ActionResult:
        self.error_on_dryrun()

        subnet_id = self._get_param("SubnetId")
        filters = self._filters_from_querystring()
        cidr_reservations = self.ec2_backend.get_subnet_cidr_reservations(
            subnet_id, filters
        )
        return ActionResult(cidr_reservations)

    def delete_subnet_cidr_reservation(self) -> ActionResult:
        self.error_on_dryrun()

        reservation_id = self._get_param("SubnetCidrReservationId")

        deleted_reservation = self.ec2_backend.delete_subnet_cidr_reservation(
            reservation_id
        )

        return ActionResult({"DeletedSubnetCidrReservation": deleted_reservation})
