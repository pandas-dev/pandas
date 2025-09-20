from moto.core.responses import ActionResult
from moto.utilities.utils import str2bool

from ._base_response import EC2BaseResponse


class TransitGatewayRouteTable(EC2BaseResponse):
    def create_transit_gateway_route_table(self) -> ActionResult:
        transit_gateway_id = self._get_param("TransitGatewayId")
        tags = self._parse_tag_specification().get("transit-gateway-route-table", {})

        transit_gateway_route_table = (
            self.ec2_backend.create_transit_gateway_route_table(
                transit_gateway_id=transit_gateway_id, tags=tags
            )
        )
        return ActionResult({"transitGatewayRouteTable": transit_gateway_route_table})

    def describe_transit_gateway_route_tables(self) -> ActionResult:
        filters = self._filters_from_querystring()
        transit_gateway_route_table_ids = (
            self._get_multi_param("TransitGatewayRouteTableIds") or None
        )
        transit_gateway_route_tables = (
            self.ec2_backend.get_all_transit_gateway_route_tables(
                transit_gateway_route_table_ids, filters
            )
        )
        return ActionResult({"transitGatewayRouteTables": transit_gateway_route_tables})

    def delete_transit_gateway_route_table(self) -> ActionResult:
        transit_gateway_route_table_id = self._get_param("TransitGatewayRouteTableId")
        transit_gateway_route_table = (
            self.ec2_backend.delete_transit_gateway_route_table(
                transit_gateway_route_table_id
            )
        )
        return ActionResult({"transitGatewayRouteTable": transit_gateway_route_table})

    def create_transit_gateway_route(self) -> ActionResult:
        transit_gateway_attachment_id = self._get_param("TransitGatewayAttachmentId")
        destination_cidr_block = self._get_param("DestinationCidrBlock")
        transit_gateway_route_table_id = self._get_param("TransitGatewayRouteTableId")
        blackhole = str2bool(self._get_param("Blackhole"))
        transit_gateways_route_table = self.ec2_backend.create_transit_gateway_route(
            destination_cidr_block=destination_cidr_block,
            transit_gateway_route_table_id=transit_gateway_route_table_id,
            transit_gateway_attachment_id=transit_gateway_attachment_id,
            blackhole=blackhole,
        )
        return ActionResult({"route": transit_gateways_route_table})

    def delete_transit_gateway_route(self) -> ActionResult:
        destination_cidr_block = self._get_param("DestinationCidrBlock")
        transit_gateway_route_table_id = self._get_param("TransitGatewayRouteTableId")
        transit_gateway_route = self.ec2_backend.delete_transit_gateway_route(
            destination_cidr_block=destination_cidr_block,
            transit_gateway_route_table_id=transit_gateway_route_table_id,
        )
        return ActionResult({"route": transit_gateway_route})

    def search_transit_gateway_routes(self) -> ActionResult:
        transit_gateway_route_table_id = self._get_param("TransitGatewayRouteTableId")
        filters = self._filters_from_querystring()
        max_results = self._get_param("MaxResults")
        transit_gateway_routes = self.ec2_backend.search_transit_gateway_routes(
            transit_gateway_route_table_id=transit_gateway_route_table_id,
            filters=filters,
            max_results=max_results,
        )
        return ActionResult(
            {"routeSet": transit_gateway_routes, "additionalRoutesAvailable": False}
        )

    def get_transit_gateway_route_table_associations(self) -> ActionResult:
        transit_gateway_route_table_id = self._get_param("TransitGatewayRouteTableId")
        filters = self._filters_from_querystring()
        transit_gateway_route_table_associations = (
            self.ec2_backend.get_transit_gateway_route_table_associations(
                transit_gateway_route_table_id, filters
            )
        )
        return ActionResult({"associations": transit_gateway_route_table_associations})

    def get_transit_gateway_route_table_propagations(self) -> ActionResult:
        transit_gateway_route_table_id = self._get_param("TransitGatewayRouteTableId")
        filters = self._filters_from_querystring()
        transit_gateway_route_table_propagations = (
            self.ec2_backend.get_transit_gateway_route_table_propagations(
                transit_gateway_route_table_id, filters
            )
        )
        return ActionResult(
            {
                "transitGatewayRouteTablePropagations": transit_gateway_route_table_propagations
            }
        )
