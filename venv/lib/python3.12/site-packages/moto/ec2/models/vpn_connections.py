from typing import Any, Dict, List, Optional

from ..exceptions import InvalidVpnConnectionIdError
from ..utils import generic_filter, random_vpn_connection_id
from .core import TaggedEC2Resource


class VPNConnection(TaggedEC2Resource):
    def __init__(
        self,
        ec2_backend: Any,
        vpn_connection_id: str,
        vpn_conn_type: str,
        customer_gateway_id: str,
        vpn_gateway_id: Optional[str] = None,
        transit_gateway_id: Optional[str] = None,
        tags: Optional[Dict[str, str]] = None,
    ):
        self.ec2_backend = ec2_backend
        self.id = vpn_connection_id
        self.state = "available"
        self.customer_gateway_configuration: Dict[str, str] = {}
        self.type = vpn_conn_type
        self.customer_gateway_id = customer_gateway_id
        self.vpn_gateway_id = vpn_gateway_id
        self.transit_gateway_id = transit_gateway_id
        self.tunnels = None
        self.options = None
        self.static_routes = None
        self.add_tags(tags or {})

    def get_filter_value(
        self, filter_name: str, method_name: Optional[str] = None
    ) -> Any:
        return super().get_filter_value(filter_name, "DescribeVpnConnections")


class VPNConnectionBackend:
    def __init__(self) -> None:
        self.vpn_connections: Dict[str, VPNConnection] = {}

    def create_vpn_connection(
        self,
        vpn_conn_type: str,
        customer_gateway_id: str,
        vpn_gateway_id: Optional[str] = None,
        transit_gateway_id: Optional[str] = None,
        static_routes_only: Optional[bool] = None,
        tags: Optional[Dict[str, str]] = None,
    ) -> VPNConnection:
        vpn_connection_id = random_vpn_connection_id()
        if static_routes_only:
            pass
        vpn_connection = VPNConnection(
            self,
            vpn_connection_id=vpn_connection_id,
            vpn_conn_type=vpn_conn_type,
            customer_gateway_id=customer_gateway_id,
            vpn_gateway_id=vpn_gateway_id,
            transit_gateway_id=transit_gateway_id,
            tags=tags,
        )
        self.vpn_connections[vpn_connection.id] = vpn_connection
        return vpn_connection

    def delete_vpn_connection(self, vpn_connection_id: str) -> VPNConnection:
        if vpn_connection_id in self.vpn_connections:
            self.vpn_connections[vpn_connection_id].state = "deleted"
        else:
            raise InvalidVpnConnectionIdError(vpn_connection_id)
        return self.vpn_connections[vpn_connection_id]

    def describe_vpn_connections(
        self, vpn_connection_ids: Optional[List[str]] = None, filters: Any = None
    ) -> List[VPNConnection]:
        vpn_connections = list(self.vpn_connections.values())

        if vpn_connection_ids:
            vpn_connections = [
                vpn_connection
                for vpn_connection in vpn_connections
                if vpn_connection.id in vpn_connection_ids
            ]
            if len(vpn_connections) != len(vpn_connection_ids):
                invalid_id = list(
                    set(vpn_connection_ids).difference(
                        set([vpn_connection.id for vpn_connection in vpn_connections])
                    )
                )[0]
                raise InvalidVpnConnectionIdError(invalid_id)

        return generic_filter(filters, vpn_connections)
