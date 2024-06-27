from typing import Any, Dict, List, Optional

from moto.core.common_models import CloudFormationModel
from moto.core.utils import iso_8601_datetime_with_milliseconds, utcnow

from ..utils import random_nat_gateway_id, random_private_ip
from .core import TaggedEC2Resource


class NatGateway(CloudFormationModel, TaggedEC2Resource):
    def __init__(
        self,
        backend: Any,
        subnet_id: str,
        allocation_id: str,
        tags: Optional[Dict[str, str]] = None,
        connectivity_type: str = "public",
    ):
        # public properties
        self.id = random_nat_gateway_id()
        self.subnet_id = subnet_id
        self.address_set: List[Dict[str, Any]] = []
        self.state = "available"
        self.private_ip = random_private_ip()
        self.connectivity_type = connectivity_type

        # protected properties
        self._created_at = utcnow()
        self.ec2_backend = backend
        # NOTE: this is the core of NAT Gateways creation
        self._eni = self.ec2_backend.create_network_interface(
            backend.get_subnet(self.subnet_id), self.private_ip
        )

        # associate allocation with ENI
        if allocation_id and connectivity_type != "private":
            self.ec2_backend.associate_address(
                eni=self._eni, allocation_id=allocation_id
            )
        self.add_tags(tags or {})
        self.vpc_id = self.ec2_backend.get_subnet(subnet_id).vpc_id

    @property
    def physical_resource_id(self) -> str:
        return self.id

    @property
    def create_time(self) -> str:
        return iso_8601_datetime_with_milliseconds(self._created_at)

    @staticmethod
    def cloudformation_name_type() -> str:
        return ""

    @staticmethod
    def cloudformation_type() -> str:
        # https://docs.aws.amazon.com/AWSCloudFormation/latest/UserGuide/aws-resource-ec2-natgateway.html
        return "AWS::EC2::NatGateway"

    @classmethod
    def create_from_cloudformation_json(  # type: ignore[misc]
        cls,
        resource_name: str,
        cloudformation_json: Any,
        account_id: str,
        region_name: str,
        **kwargs: Any,
    ) -> "NatGateway":
        from ..models import ec2_backends

        ec2_backend = ec2_backends[account_id][region_name]
        nat_gateway = ec2_backend.create_nat_gateway(
            cloudformation_json["Properties"]["SubnetId"],
            cloudformation_json["Properties"]["AllocationId"],
        )
        return nat_gateway


class NatGatewayBackend:
    def __init__(self) -> None:
        self.nat_gateways: Dict[str, NatGateway] = {}

    def describe_nat_gateways(
        self, filters: Any, nat_gateway_ids: Optional[List[str]]
    ) -> List[NatGateway]:
        nat_gateways = list(self.nat_gateways.values())

        if nat_gateway_ids:
            nat_gateways = [item for item in nat_gateways if item.id in nat_gateway_ids]

        if filters is not None:
            if filters.get("nat-gateway-id") is not None:
                nat_gateways = [
                    nat_gateway
                    for nat_gateway in nat_gateways
                    if nat_gateway.id in filters["nat-gateway-id"]
                ]
            if filters.get("vpc-id") is not None:
                nat_gateways = [
                    nat_gateway
                    for nat_gateway in nat_gateways
                    if nat_gateway.vpc_id in filters["vpc-id"]
                ]
            if filters.get("subnet-id") is not None:
                nat_gateways = [
                    nat_gateway
                    for nat_gateway in nat_gateways
                    if nat_gateway.subnet_id in filters["subnet-id"]
                ]
            if filters.get("state") is not None:
                nat_gateways = [
                    nat_gateway
                    for nat_gateway in nat_gateways
                    if nat_gateway.state in filters["state"]
                ]

        return nat_gateways

    def create_nat_gateway(
        self,
        subnet_id: str,
        allocation_id: str,
        tags: Optional[Dict[str, str]] = None,
        connectivity_type: str = "public",
    ) -> NatGateway:
        nat_gateway = NatGateway(
            self, subnet_id, allocation_id, tags, connectivity_type
        )
        address_set: Dict[str, Any] = {}
        if allocation_id:
            eips = self.address_by_allocation([allocation_id])  # type: ignore[attr-defined]
            eip = eips[0] if len(eips) > 0 else None
            if eip:
                address_set["allocationId"] = allocation_id
                address_set["publicIp"] = eip.public_ip or None
                address_set["associationId"] = eip.association_id or None
        address_set["networkInterfaceId"] = nat_gateway._eni.id
        address_set["privateIp"] = nat_gateway._eni.private_ip_address
        nat_gateway.address_set.append(address_set)
        self.nat_gateways[nat_gateway.id] = nat_gateway
        return nat_gateway

    def delete_nat_gateway(self, nat_gateway_id: str) -> NatGateway:
        nat_gw: NatGateway = self.nat_gateways.get(nat_gateway_id)  # type: ignore
        nat_gw.state = "deleted"
        return nat_gw
