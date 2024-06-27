from typing import Any, Dict, List, Optional

from ..exceptions import InvalidCustomerGatewayIdError
from ..utils import random_customer_gateway_id
from .core import TaggedEC2Resource


class CustomerGateway(TaggedEC2Resource):
    def __init__(
        self,
        ec2_backend: Any,
        gateway_id: str,
        gateway_type: str,
        ip_address: str,
        bgp_asn: str,
        state: str = "available",
        tags: Optional[Dict[str, str]] = None,
    ):
        self.ec2_backend = ec2_backend
        self.id = gateway_id
        self.type = gateway_type or "ipsec.1"
        self.ip_address = ip_address
        self.bgp_asn = bgp_asn
        self.state = state
        self.add_tags(tags or {})
        super().__init__()

    def get_filter_value(
        self, filter_name: str, method_name: Optional[str] = None
    ) -> Any:
        return super().get_filter_value(filter_name, "DescribeCustomerGateways")


class CustomerGatewayBackend:
    def __init__(self) -> None:
        self.customer_gateways: Dict[str, CustomerGateway] = {}

    def create_customer_gateway(
        self,
        gateway_type: str,
        ip_address: str,
        bgp_asn: str,
        tags: Optional[Dict[str, str]] = None,
    ) -> CustomerGateway:
        customer_gateway_id = random_customer_gateway_id()
        customer_gateway = CustomerGateway(
            self, customer_gateway_id, gateway_type, ip_address, bgp_asn, tags=tags
        )
        self.customer_gateways[customer_gateway_id] = customer_gateway
        return customer_gateway

    def describe_customer_gateways(
        self, filters: Any = None, customer_gateway_ids: Optional[List[str]] = None
    ) -> List[CustomerGateway]:
        customer_gateways = list(self.customer_gateways.copy().values())
        if customer_gateway_ids:
            customer_gateways = [
                cg for cg in customer_gateways if cg.id in customer_gateway_ids
            ]

        if filters is not None:
            if filters.get("customer-gateway-id") is not None:
                customer_gateways = [
                    customer_gateway
                    for customer_gateway in customer_gateways
                    if customer_gateway.id in filters["customer-gateway-id"]
                ]
            if filters.get("type") is not None:
                customer_gateways = [
                    customer_gateway
                    for customer_gateway in customer_gateways
                    if customer_gateway.type in filters["type"]
                ]
            if filters.get("bgp-asn") is not None:
                customer_gateways = [
                    customer_gateway
                    for customer_gateway in customer_gateways
                    if customer_gateway.bgp_asn in filters["bgp-asn"]
                ]
            if filters.get("ip-address") is not None:
                customer_gateways = [
                    customer_gateway
                    for customer_gateway in customer_gateways
                    if customer_gateway.ip_address in filters["ip-address"]
                ]
        return customer_gateways

    def get_customer_gateway(self, customer_gateway_id: str) -> CustomerGateway:
        customer_gateway = self.customer_gateways.get(customer_gateway_id)
        if not customer_gateway:
            raise InvalidCustomerGatewayIdError(customer_gateway_id)
        return customer_gateway

    def delete_customer_gateway(self, customer_gateway_id: str) -> None:
        customer_gateway = self.get_customer_gateway(customer_gateway_id)
        customer_gateway.state = "deleted"
