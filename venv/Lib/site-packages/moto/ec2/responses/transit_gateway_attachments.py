from moto.core.responses import ActionResult
from moto.ec2.utils import add_tag_specification

from ._base_response import EC2BaseResponse


class TransitGatewayAttachment(EC2BaseResponse):
    def create_transit_gateway_vpc_attachment(self) -> ActionResult:
        options = self._get_param("Options", {})
        subnet_ids = self._get_param("SubnetIds", [])
        transit_gateway_id = self._get_param("TransitGatewayId")
        vpc_id = self._get_param("VpcId")

        tags = self._parse_tag_specification().get("transit-gateway-route-table", {})

        transit_gateway_attachment = (
            self.ec2_backend.create_transit_gateway_vpc_attachment(
                transit_gateway_id=transit_gateway_id,
                tags=tags,
                vpc_id=vpc_id,
                subnet_ids=subnet_ids,
                options=options,
            )
        )
        return ActionResult({"TransitGatewayVpcAttachment": transit_gateway_attachment})

    def describe_transit_gateway_vpc_attachments(self) -> ActionResult:
        transit_gateways_attachment_ids = self._get_param(
            "TransitGatewayAttachmentIds", []
        )
        filters = self._filters_from_querystring()
        transit_gateway_vpc_attachments = (
            self.ec2_backend.describe_transit_gateway_vpc_attachments(
                transit_gateways_attachment_ids=transit_gateways_attachment_ids,
                filters=filters,
            )
        )
        return ActionResult(
            {"TransitGatewayVpcAttachments": transit_gateway_vpc_attachments}
        )

    def modify_transit_gateway_vpc_attachment(self) -> ActionResult:
        add_subnet_ids = self._get_param("AddSubnetIds", [])
        options = self._get_param("Options", {})
        remove_subnet_ids = self._get_param("RemoveSubnetIds", [])
        transit_gateway_attachment_id = self._get_param("TransitGatewayAttachmentId")

        transit_gateway_attachment = (
            self.ec2_backend.modify_transit_gateway_vpc_attachment(
                add_subnet_ids=add_subnet_ids,
                options=options,
                remove_subnet_ids=remove_subnet_ids,
                transit_gateway_attachment_id=transit_gateway_attachment_id,
            )
        )
        return ActionResult({"TransitGatewayVpcAttachment": transit_gateway_attachment})

    def describe_transit_gateway_attachments(self) -> ActionResult:
        transit_gateways_attachment_ids = self._get_param(
            "TransitGatewayAttachmentIds", []
        )
        filters = self._filters_from_querystring()
        transit_gateway_attachments = (
            self.ec2_backend.describe_transit_gateway_attachments(
                transit_gateways_attachment_ids=transit_gateways_attachment_ids,
                filters=filters,
            )
        )
        return ActionResult({"TransitGatewayAttachments": transit_gateway_attachments})

    def delete_transit_gateway_vpc_attachment(self) -> ActionResult:
        transit_gateway_attachment_id = self._get_param("TransitGatewayAttachmentId")
        transit_gateway_attachment = (
            self.ec2_backend.delete_transit_gateway_vpc_attachment(
                transit_gateway_attachment_id=transit_gateway_attachment_id
            )
        )
        return ActionResult({"TransitGatewayVpcAttachment": transit_gateway_attachment})

    def associate_transit_gateway_route_table(self) -> ActionResult:
        transit_gateway_attachment_id = self._get_param("TransitGatewayAttachmentId")
        transit_gateway_route_table_id = self._get_param("TransitGatewayRouteTableId")
        transit_gateway_association = (
            self.ec2_backend.associate_transit_gateway_route_table(
                transit_gateway_attachment_id=transit_gateway_attachment_id,
                transit_gateway_route_table_id=transit_gateway_route_table_id,
            )
        )
        return ActionResult({"Association": transit_gateway_association})

    def disassociate_transit_gateway_route_table(self) -> ActionResult:
        tgw_attach_id = self._get_param("TransitGatewayAttachmentId")
        tgw_rt_id = self._get_param("TransitGatewayRouteTableId")

        tgw_association = self.ec2_backend.disassociate_transit_gateway_route_table(
            tgw_attach_id, tgw_rt_id
        )
        return ActionResult({"Association": tgw_association})

    def enable_transit_gateway_route_table_propagation(self) -> ActionResult:
        transit_gateway_attachment_id = self._get_param("TransitGatewayAttachmentId")
        transit_gateway_route_table_id = self._get_param("TransitGatewayRouteTableId")
        transit_gateway_propagation = (
            self.ec2_backend.enable_transit_gateway_route_table_propagation(
                transit_gateway_attachment_id=transit_gateway_attachment_id,
                transit_gateway_route_table_id=transit_gateway_route_table_id,
            )
        )
        return ActionResult({"Propagation": transit_gateway_propagation})

    def disable_transit_gateway_route_table_propagation(self) -> ActionResult:
        transit_gateway_attachment_id = self._get_param("TransitGatewayAttachmentId")
        transit_gateway_route_table_id = self._get_param("TransitGatewayRouteTableId")
        transit_gateway_propagation = (
            self.ec2_backend.disable_transit_gateway_route_table_propagation(
                transit_gateway_attachment_id=transit_gateway_attachment_id,
                transit_gateway_route_table_id=transit_gateway_route_table_id,
            )
        )
        return ActionResult({"Propagation": transit_gateway_propagation})

    def create_transit_gateway_peering_attachment(self) -> ActionResult:
        peer_account_id = self._get_param("PeerAccountId")
        peer_region = self._get_param("PeerRegion")
        peer_transit_gateway_id = self._get_param("PeerTransitGatewayId")
        transit_gateway_id = self._get_param("TransitGatewayId")
        tags = add_tag_specification(self._get_param("TagSpecifications", []))
        transit_gateway_peering_attachment = (
            self.ec2_backend.create_transit_gateway_peering_attachment(
                transit_gateway_id,
                peer_transit_gateway_id,
                peer_region,
                peer_account_id,
                tags,
            )
        )
        return ActionResult(
            {"TransitGatewayPeeringAttachment": transit_gateway_peering_attachment}
        )

    def describe_transit_gateway_peering_attachments(self) -> ActionResult:
        transit_gateways_attachment_ids = self._get_param(
            "TransitGatewayAttachmentIds", []
        )
        filters = self._filters_from_querystring()
        transit_gateway_peering_attachments = (
            self.ec2_backend.describe_transit_gateway_peering_attachments(
                transit_gateways_attachment_ids=transit_gateways_attachment_ids,
                filters=filters,
            )
        )
        return ActionResult(
            {"TransitGatewayPeeringAttachments": transit_gateway_peering_attachments}
        )

    def accept_transit_gateway_peering_attachment(self) -> ActionResult:
        transit_gateway_attachment_id = self._get_param("TransitGatewayAttachmentId")
        transit_gateway_peering_attachment = (
            self.ec2_backend.accept_transit_gateway_peering_attachment(
                transit_gateway_attachment_id=transit_gateway_attachment_id
            )
        )
        return ActionResult(
            {"TransitGatewayPeeringAttachment": transit_gateway_peering_attachment}
        )

    def delete_transit_gateway_peering_attachment(self) -> ActionResult:
        transit_gateway_attachment_id = self._get_param("TransitGatewayAttachmentId")
        transit_gateway_peering_attachment = (
            self.ec2_backend.delete_transit_gateway_peering_attachment(
                transit_gateway_attachment_id=transit_gateway_attachment_id
            )
        )
        return ActionResult(
            {"TransitGatewayPeeringAttachment": transit_gateway_peering_attachment}
        )

    def reject_transit_gateway_peering_attachment(self) -> ActionResult:
        transit_gateway_attachment_id = self._get_param("TransitGatewayAttachmentId")
        transit_gateway_peering_attachment = (
            self.ec2_backend.reject_transit_gateway_peering_attachment(
                transit_gateway_attachment_id=transit_gateway_attachment_id
            )
        )
        return ActionResult(
            {"TransitGatewayPeeringAttachment": transit_gateway_peering_attachment}
        )
