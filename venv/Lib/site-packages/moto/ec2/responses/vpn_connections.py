from moto.core.responses import ActionResult, EmptyResult
from moto.ec2.utils import add_tag_specification

from ._base_response import EC2BaseResponse


class VPNConnections(EC2BaseResponse):
    def create_vpn_connection(self) -> ActionResult:
        vpn_conn_type = self._get_param("Type")
        cgw_id = self._get_param("CustomerGatewayId")
        vgw_id = self._get_param("VpnGatewayId")
        tgw_id = self._get_param("TransitGatewayId")
        tags = add_tag_specification(self._get_param("TagSpecifications", []))
        vpn_connection = self.ec2_backend.create_vpn_connection(
            vpn_conn_type,
            cgw_id,
            vpn_gateway_id=vgw_id,
            transit_gateway_id=tgw_id,
            tags=tags,
        )
        if vpn_connection.transit_gateway_id:
            self.ec2_backend.create_transit_gateway_vpn_attachment(
                vpn_id=vpn_connection.id, transit_gateway_id=tgw_id
            )
        return ActionResult({"VpnConnection": vpn_connection})

    def delete_vpn_connection(self) -> ActionResult:
        vpn_connection_id = self._get_param("VpnConnectionId")
        vpn_connection = self.ec2_backend.delete_vpn_connection(vpn_connection_id)
        if vpn_connection.transit_gateway_id:
            transit_gateway_attachments = (
                self.ec2_backend.describe_transit_gateway_attachments()
            )
            for attachment in transit_gateway_attachments:
                if attachment.resource_id == vpn_connection.id:
                    attachment.state = "deleted"
        return EmptyResult()

    def describe_vpn_connections(self) -> ActionResult:
        vpn_connection_ids = self._get_param("VpnConnectionIds", [])
        filters = self._filters_from_querystring()
        vpn_connections = self.ec2_backend.describe_vpn_connections(
            vpn_connection_ids=vpn_connection_ids, filters=filters
        )
        return ActionResult({"VpnConnections": vpn_connections})
