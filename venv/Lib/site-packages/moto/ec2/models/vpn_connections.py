from typing import Any

from ..exceptions import (
    InvalidParameterValue,
    InvalidTransitGatewayID,
    InvalidVpnConnectionIdError,
)
from ..utils import generic_filter, random_vpn_connection_id
from .core import TaggedEC2Resource


class VPNConnection(TaggedEC2Resource):
    def __init__(
        self,
        ec2_backend: Any,
        vpn_connection_id: str,
        vpn_conn_type: str,
        customer_gateway_id: str,
        vpn_gateway_id: str | None = None,
        transit_gateway_id: str | None = None,
        tags: dict[str, str] | None = None,
    ):
        self.ec2_backend = ec2_backend
        self.id = vpn_connection_id
        self.state = "available"
        self.type = vpn_conn_type
        self.customer_gateway_id = customer_gateway_id
        self.vpn_gateway_id = vpn_gateway_id
        self.transit_gateway_id = transit_gateway_id
        self.tunnels = None
        self.options = None
        self.static_routes = None
        self.add_tags(tags or {})

    @property
    def customer_gateway_configuration(self) -> str:
        vgw_id = self.vpn_gateway_id if self.vpn_gateway_id is not None else ""
        return _CUSTOMER_GATEWAY_CONFIGURATION.format(
            vpn_connection_id=self.id,
            customer_gateway_id=self.customer_gateway_id,
            vpn_gateway_id=vgw_id,
            vpn_connection_type=self.type,
        )

    def get_filter_value(self, filter_name: str, method_name: str | None = None) -> Any:
        return super().get_filter_value(filter_name, "DescribeVpnConnections")


class VPNConnectionBackend:
    def __init__(self) -> None:
        self.vpn_connections: dict[str, VPNConnection] = {}

    def create_vpn_connection(
        self,
        vpn_conn_type: str,
        customer_gateway_id: str,
        vpn_gateway_id: str | None = None,
        transit_gateway_id: str | None = None,
        tags: dict[str, str] | None = None,
    ) -> VPNConnection:
        if vpn_gateway_id and transit_gateway_id:
            # From the docs (parentheses mine):
            # Creates a VPN connection between an existing (virtual private gateway or transit gateway) and a customer gateway
            raise InvalidParameterValue(
                "The request must not contain both parameter vpnGatewayId and transitGatewayId"
            )
        # Validate Gateways exist
        self.get_customer_gateway(customer_gateway_id=customer_gateway_id)  # type: ignore[attr-defined]
        if transit_gateway_id and transit_gateway_id not in self.transit_gateways:  # type: ignore[attr-defined]
            raise InvalidTransitGatewayID(transit_gateway_id)

        vpn_connection_id = random_vpn_connection_id()
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
        self, vpn_connection_ids: list[str] | None = None, filters: Any = None
    ) -> list[VPNConnection]:
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
                        {vpn_connection.id for vpn_connection in vpn_connections}
                    )
                )[0]
                raise InvalidVpnConnectionIdError(invalid_id)

        return generic_filter(filters, vpn_connections)


_IPSEC_TUNNEL = """\
          <ipsec_tunnel>
            <customer_gateway>
            <tunnel_outside_address>
              <ip_address>12.1.2.3</ip_address>
            </tunnel_outside_address>
            <tunnel_inside_address>
              <ip_address>169.254.44.42</ip_address>
              <network_mask>255.255.255.252</network_mask>
              <network_cidr>30</network_cidr>
            </tunnel_inside_address>
            <bgp>
              <asn>65000</asn>
              <hold_time>30</hold_time>
            </bgp>
            </customer_gateway>
            <vpn_gateway>
            <tunnel_outside_address>
              <ip_address>52.2.144.13</ip_address>
            </tunnel_outside_address>
            <tunnel_inside_address>
              <ip_address>169.254.44.41</ip_address>
              <network_mask>255.255.255.252</network_mask>
              <network_cidr>30</network_cidr>
            </tunnel_inside_address>
            <bgp>
              <asn>7224</asn>
              <hold_time>30</hold_time>
            </bgp>
            </vpn_gateway>
            <ike>
            <authentication_protocol>sha1</authentication_protocol>
            <encryption_protocol>aes-128-cbc</encryption_protocol>
            <lifetime>28800</lifetime>
            <perfect_forward_secrecy>group2</perfect_forward_secrecy>
            <mode>main</mode>
            <pre_shared_key>Iw2IAN9XUsQeYUrkMGP3kP59ugFDkfHg</pre_shared_key>
            </ike>
            <ipsec>
            <protocol>esp</protocol>
            <authentication_protocol>hmac-sha1-96</authentication_protocol>
            <encryption_protocol>aes-128-cbc</encryption_protocol>
            <lifetime>3600</lifetime>
            <perfect_forward_secrecy>group2</perfect_forward_secrecy>
            <mode>tunnel</mode>
            <clear_df_bit>true</clear_df_bit>
            <fragmentation_before_encryption>true</fragmentation_before_encryption>
            <tcp_mss_adjustment>1387</tcp_mss_adjustment>
            <dead_peer_detection>
              <interval>10</interval>
              <retries>3</retries>
            </dead_peer_detection>
            </ipsec>
          </ipsec_tunnel>"""

_CUSTOMER_GATEWAY_CONFIGURATION = (
    """\
          <vpn_connection id="{vpn_connection_id}">
          <customer_gateway_id>{customer_gateway_id}</customer_gateway_id>
          <vpn_gateway_id> {vpn_gateway_id} </vpn_gateway_id>
          <vpn_connection_type>{vpn_connection_type}</vpn_connection_type>
"""
    + _IPSEC_TUNNEL
    + "\n"
    + _IPSEC_TUNNEL
    + """
        </vpn_connection>
"""
)
