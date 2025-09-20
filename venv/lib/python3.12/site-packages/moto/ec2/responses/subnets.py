from moto.core.responses import ActionResult, EmptyResult
from moto.core.utils import camelcase_to_underscores
from moto.moto_api._internal import mock_random as random

from ._base_response import EC2BaseResponse


class Subnets(EC2BaseResponse):
    def create_subnet(self) -> ActionResult:
        vpc_id = self._get_param("VpcId")
        cidr_block = self._get_param("CidrBlock")
        ipv6_cidr_block = self._get_param("Ipv6CidrBlock")
        availability_zone = self._get_param("AvailabilityZone")
        availability_zone_id = self._get_param("AvailabilityZoneId")
        tags = self._parse_tag_specification("subnet")

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
        subnet_ids = self._get_multi_param("SubnetId")
        filters = self._filters_from_querystring()
        subnets = self.ec2_backend.describe_subnets(subnet_ids, filters)
        result = {"Subnets": subnets}
        return ActionResult(result)

    def modify_subnet_attribute(self) -> ActionResult:
        subnet_id = self._get_param("SubnetId")

        for attribute in ("MapPublicIpOnLaunch", "AssignIpv6AddressOnCreation"):
            if self.querystring.get(f"{attribute}.Value"):
                attr_name = camelcase_to_underscores(attribute)
                attr_value = self.querystring[f"{attribute}.Value"][0]
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
