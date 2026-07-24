from moto.core.responses import ActionResult, EmptyResult

from ._base_response import EC2BaseResponse


class RouteTables(EC2BaseResponse):
    def associate_route_table(self) -> ActionResult:
        route_table_id = self._get_param("RouteTableId")
        gateway_id = self._get_param("GatewayId")
        subnet_id = self._get_param("SubnetId")
        association_id = self.ec2_backend.associate_route_table(
            route_table_id, gateway_id, subnet_id
        )
        return ActionResult(
            {
                "AssociationId": association_id,
                "AssociationState": {"State": "associated"},
            }
        )

    def create_route(self) -> ActionResult:
        route_table_id = self._get_param("RouteTableId")
        destination_cidr_block = self._get_param("DestinationCidrBlock")
        destination_ipv6_cidr_block = self._get_param("DestinationIpv6CidrBlock")
        destination_prefix_list_id = self._get_param("DestinationPrefixListId")
        gateway_id = self._get_param("GatewayId")
        instance_id = self._get_param("InstanceId")
        nat_gateway_id = self._get_param("NatGatewayId")
        egress_only_igw_id = self._get_param("EgressOnlyInternetGatewayId")
        transit_gateway_id = self._get_param("TransitGatewayId")
        interface_id = self._get_param("NetworkInterfaceId")
        pcx_id = self._get_param("VpcPeeringConnectionId")
        carrier_gateway_id = self._get_param("CarrierGatewayId")
        vpc_endpoint_id = self._get_param("VpcEndpointId")

        self.ec2_backend.create_route(
            route_table_id,
            destination_cidr_block,
            destination_ipv6_cidr_block,
            destination_prefix_list_id,
            gateway_id=gateway_id,
            instance_id=instance_id,
            nat_gateway_id=nat_gateway_id,
            egress_only_igw_id=egress_only_igw_id,
            transit_gateway_id=transit_gateway_id,
            interface_id=interface_id,
            vpc_peering_connection_id=pcx_id,
            carrier_gateway_id=carrier_gateway_id,
            vpc_endpoint_id=vpc_endpoint_id,
        )

        return ActionResult({"Return": True})

    def create_route_table(self) -> ActionResult:
        vpc_id = self._get_param("VpcId")
        tags = self._get_param("TagSpecifications", [])
        if tags:
            tags = tags[0].get("Tags") or []
        route_table = self.ec2_backend.create_route_table(vpc_id, tags)
        return ActionResult({"RouteTable": route_table})

    def delete_route(self) -> ActionResult:
        route_table_id = self._get_param("RouteTableId")
        destination_cidr_block = self._get_param("DestinationCidrBlock")
        destination_ipv6_cidr_block = self._get_param("DestinationIpv6CidrBlock")
        destination_prefix_list_id = self._get_param("DestinationPrefixListId")
        self.ec2_backend.delete_route(
            route_table_id,
            destination_cidr_block,
            destination_ipv6_cidr_block,
            destination_prefix_list_id,
        )
        return EmptyResult()

    def delete_route_table(self) -> ActionResult:
        route_table_id = self._get_param("RouteTableId")
        self.ec2_backend.delete_route_table(route_table_id)
        return EmptyResult()

    def describe_route_tables(self) -> ActionResult:
        route_table_ids = self._get_param("RouteTableIds", [])
        filters = self._filters_from_querystring()
        route_tables = self.ec2_backend.describe_route_tables(route_table_ids, filters)
        return ActionResult({"RouteTables": route_tables})

    def disassociate_route_table(self) -> ActionResult:
        association_id = self._get_param("AssociationId")
        self.ec2_backend.disassociate_route_table(association_id)
        return EmptyResult()

    def replace_route(self) -> ActionResult:
        route_table_id = self._get_param("RouteTableId")
        destination_cidr_block = self._get_param("DestinationCidrBlock")
        destination_ipv6_cidr_block = self._get_param("DestinationIpv6CidrBlock")
        destination_prefix_list_id = self._get_param("DestinationPrefixListId")
        gateway_id = self._get_param("GatewayId")
        instance_id = self._get_param("InstanceId")
        interface_id = self._get_param("NetworkInterfaceId")
        pcx_id = self._get_param("VpcPeeringConnectionId")
        nat_gateway_id = self._get_param("NatGatewayId")
        egress_only_igw_id = self._get_param("EgressOnlyInternetGatewayId")
        transit_gateway_id = self._get_param("TransitGatewayId")

        self.ec2_backend.replace_route(
            route_table_id,
            destination_cidr_block,
            destination_ipv6_cidr_block,
            destination_prefix_list_id,
            nat_gateway_id,
            egress_only_igw_id,
            transit_gateway_id,
            gateway_id=gateway_id,
            instance_id=instance_id,
            interface_id=interface_id,
            vpc_peering_connection_id=pcx_id,
        )

        return EmptyResult()

    def replace_route_table_association(self) -> ActionResult:
        route_table_id = self._get_param("RouteTableId")
        association_id = self._get_param("AssociationId")
        new_association_id = self.ec2_backend.replace_route_table_association(
            association_id, route_table_id
        )
        return ActionResult(
            {
                "NewAssociationId": new_association_id,
                "AssociationState": {"State": "associated"},
            }
        )
