from typing import Any, Dict, List, Optional

from moto.utilities.utils import filter_resources

from ..exceptions import InvalidCarrierGatewayID, InvalidVPCIdError
from ..utils import random_carrier_gateway_id
from .core import TaggedEC2Resource


class CarrierGateway(TaggedEC2Resource):
    def __init__(
        self, ec2_backend: Any, vpc_id: str, tags: Optional[Dict[str, str]] = None
    ):
        self.id = random_carrier_gateway_id()
        self.ec2_backend = ec2_backend
        self.vpc_id = vpc_id
        self.state = "available"
        self.add_tags(tags or {})

    @property
    def physical_resource_id(self) -> str:
        return self.id

    @property
    def owner_id(self) -> str:
        return self.ec2_backend.account_id


class CarrierGatewayBackend:
    def __init__(self) -> None:
        self.carrier_gateways: Dict[str, CarrierGateway] = {}

    def create_carrier_gateway(
        self, vpc_id: str, tags: Optional[Dict[str, str]] = None
    ) -> CarrierGateway:
        vpc = self.get_vpc(vpc_id)  # type: ignore[attr-defined]
        if not vpc:
            raise InvalidVPCIdError(vpc_id)
        carrier_gateway = CarrierGateway(self, vpc_id, tags)
        self.carrier_gateways[carrier_gateway.id] = carrier_gateway
        return carrier_gateway

    def delete_carrier_gateway(self, gateway_id: str) -> CarrierGateway:
        if not self.carrier_gateways.get(gateway_id):
            raise InvalidCarrierGatewayID(gateway_id)
        carrier_gateway = self.carrier_gateways.pop(gateway_id)
        carrier_gateway.state = "deleted"
        return carrier_gateway

    def describe_carrier_gateways(
        self, ids: Optional[List[str]] = None, filters: Any = None
    ) -> List[CarrierGateway]:
        carrier_gateways = list(self.carrier_gateways.values())

        if ids:
            carrier_gateways = [
                carrier_gateway
                for carrier_gateway in carrier_gateways
                if carrier_gateway.id in ids
            ]

        attr_pairs = (
            ("carrier-gateway-id", "id"),
            ("state", "state"),
            ("vpc-id", "vpc_id"),
            ("owner-id", "owner_id"),
        )

        result = carrier_gateways
        if filters:
            result = filter_resources(carrier_gateways, filters, attr_pairs)
        return result
