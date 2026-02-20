from moto.core.responses import ActionResult

from ._base_response import EC2BaseResponse


class TransitGateways(EC2BaseResponse):
    def create_transit_gateway(self) -> ActionResult:
        description = self._get_param("Description") or None
        options = self._get_param("Options", {})
        tags = self._get_param("TagSpecifications", [])
        if tags:
            tags = tags[0].get("Tags")

        transit_gateway = self.ec2_backend.create_transit_gateway(
            description=description, options=options, tags=tags
        )

        return ActionResult({"TransitGateway": transit_gateway})

    def delete_transit_gateway(self) -> ActionResult:
        transit_gateway_id = self._get_param("TransitGatewayId")
        transit_gateway = self.ec2_backend.delete_transit_gateway(transit_gateway_id)
        return ActionResult({"TransitGateway": transit_gateway})

    def describe_transit_gateways(self) -> ActionResult:
        transit_gateway_ids = self._get_param("TransitGatewayIds", [])
        filters = self._filters_from_querystring()
        transit_gateways = self.ec2_backend.describe_transit_gateways(
            filters, transit_gateway_ids
        )
        return ActionResult({"TransitGateways": transit_gateways})

    def modify_transit_gateway(self) -> ActionResult:
        transit_gateway_id = self._get_param("TransitGatewayId")
        description = self._get_param("Description") or None
        options = self._get_param("Options", {})
        transit_gateway = self.ec2_backend.modify_transit_gateway(
            transit_gateway_id=transit_gateway_id,
            description=description,
            options=options,
        )
        return ActionResult({"TransitGateway": transit_gateway})
