import ipaddress
from typing import Any, Dict, List, Optional, Set

from moto.core.common_models import CloudFormationModel
from moto.ec2.models.carrier_gateways import CarrierGateway
from moto.ec2.models.elastic_network_interfaces import NetworkInterface
from moto.ec2.models.instances import Instance
from moto.ec2.models.internet_gateways import EgressOnlyInternetGateway
from moto.ec2.models.managed_prefixes import ManagedPrefixList
from moto.ec2.models.nat_gateways import NatGateway
from moto.ec2.models.transit_gateway import TransitGateway
from moto.ec2.models.vpc_peering_connections import VPCPeeringConnection
from moto.ec2.models.vpn_gateway import VpnGateway

from ..exceptions import (
    DependencyViolationError,
    InvalidAssociationIdError,
    InvalidDestinationCIDRBlockParameterError,
    InvalidParameterValueErrorReplaceRoute,
    InvalidRouteError,
    InvalidRouteTableIdError,
    RouteAlreadyExistsError,
    RouteNotSupportedError,
)
from ..utils import (
    EC2_RESOURCE_TO_PREFIX,
    generate_route_id,
    generic_filter,
    random_route_table_id,
    random_subnet_association_id,
)
from .core import TaggedEC2Resource


class RouteTable(TaggedEC2Resource, CloudFormationModel):
    def __init__(
        self, ec2_backend: Any, route_table_id: str, vpc_id: str, main: bool = False
    ):
        self.ec2_backend = ec2_backend
        self.id = route_table_id
        self.vpc_id = vpc_id
        self.main_association_id = random_subnet_association_id() if main else None
        self.associations: Dict[str, str] = {}
        self.routes: Dict[str, Route] = {}

    @property
    def owner_id(self) -> str:
        return self.ec2_backend.account_id

    @staticmethod
    def cloudformation_name_type() -> str:
        return ""

    @staticmethod
    def cloudformation_type() -> str:
        # https://docs.aws.amazon.com/AWSCloudFormation/latest/UserGuide/aws-resource-ec2-routetable.html
        return "AWS::EC2::RouteTable"

    @classmethod
    def create_from_cloudformation_json(  # type: ignore[misc]
        cls,
        resource_name: str,
        cloudformation_json: Any,
        account_id: str,
        region_name: str,
        **kwargs: Any,
    ) -> "RouteTable":
        from ..models import ec2_backends

        properties = cloudformation_json["Properties"]

        vpc_id = properties["VpcId"]
        ec2_backend = ec2_backends[account_id][region_name]
        route_table = ec2_backend.create_route_table(vpc_id=vpc_id)
        return route_table

    @property
    def physical_resource_id(self) -> str:
        return self.id

    def get_filter_value(
        self, filter_name: str, method_name: Optional[str] = None
    ) -> Any:
        if filter_name == "association.main":
            # Note: Boto only supports 'true'.
            # https://github.com/boto/boto/issues/1742
            if self.main_association_id is not None:
                return "true"
            else:
                return "false"
        elif filter_name == "route-table-id":
            return self.id
        elif filter_name == "vpc-id":
            return self.vpc_id
        elif filter_name == "association.route-table-id":
            return self.id
        elif filter_name == "association.route-table-association-id":
            return self.all_associations_ids
        elif filter_name == "association.subnet-id":
            return self.associations.values()
        elif filter_name == "route.destination-cidr-block":
            return [
                route.destination_cidr_block
                for route in self.routes.values()
                if route.destination_cidr_block is not None
            ]
        elif filter_name == "route.gateway-id":
            return [
                route.gateway.id
                for route in self.routes.values()
                if route.gateway is not None
            ]
        elif filter_name == "route.vpc-peering-connection-id":
            return [
                route.vpc_pcx.id
                for route in self.routes.values()
                if route.vpc_pcx is not None
            ]
        elif filter_name == "route.nat-gateway-id":
            return [
                route.nat_gateway.id
                for route in self.routes.values()
                if route.nat_gateway is not None
            ]
        else:
            return super().get_filter_value(filter_name, "DescribeRouteTables")

    @property
    def all_associations_ids(self) -> Set[str]:
        # NOTE(yoctozepto): Doing an explicit copy to not touch the original.
        all_associations = set(self.associations)
        if self.main_association_id is not None:
            all_associations.add(self.main_association_id)
        return all_associations


class Route(CloudFormationModel):
    def __init__(
        self,
        route_table: RouteTable,
        destination_cidr_block: Optional[str],
        destination_ipv6_cidr_block: Optional[str],
        destination_prefix_list: Optional[ManagedPrefixList] = None,
        local: bool = False,
        gateway: Optional[VpnGateway] = None,
        instance: Optional[Instance] = None,
        nat_gateway: Optional[NatGateway] = None,
        egress_only_igw: Optional[EgressOnlyInternetGateway] = None,
        transit_gateway: Optional[TransitGateway] = None,
        interface: Optional[NetworkInterface] = None,
        vpc_pcx: Optional[VPCPeeringConnection] = None,
        carrier_gateway: Optional[CarrierGateway] = None,
        vpc_endpoint_id: Optional[str] = None,
    ):
        self.id = generate_route_id(
            route_table.id,
            destination_cidr_block,
            destination_ipv6_cidr_block,
            destination_prefix_list.id if destination_prefix_list else None,
        )
        self.route_table = route_table
        self.destination_cidr_block = destination_cidr_block
        self.destination_ipv6_cidr_block = destination_ipv6_cidr_block
        self.destination_prefix_list = destination_prefix_list
        self.local = local
        self.gateway = gateway
        self.instance = instance
        self.nat_gateway = nat_gateway
        self.egress_only_igw = egress_only_igw
        self.transit_gateway = transit_gateway
        self.interface = interface
        self.vpc_pcx = vpc_pcx
        self.carrier_gateway = carrier_gateway
        self.vpc_endpoint_id = vpc_endpoint_id

    @property
    def physical_resource_id(self) -> str:
        return self.id

    @staticmethod
    def cloudformation_name_type() -> str:
        return ""

    @staticmethod
    def cloudformation_type() -> str:
        # https://docs.aws.amazon.com/AWSCloudFormation/latest/UserGuide/aws-resource-ec2-route.html
        return "AWS::EC2::Route"

    @classmethod
    def create_from_cloudformation_json(  # type: ignore[misc]
        cls,
        resource_name: str,
        cloudformation_json: Any,
        account_id: str,
        region_name: str,
        **kwargs: str,
    ) -> "Route":
        from ..models import ec2_backends

        properties = cloudformation_json["Properties"]

        gateway_id = properties.get("GatewayId")
        instance_id = properties.get("InstanceId")
        interface_id = properties.get("NetworkInterfaceId")
        nat_gateway_id = properties.get("NatGatewayId")
        egress_only_igw_id = properties.get("EgressOnlyInternetGatewayId")
        transit_gateway_id = properties.get("TransitGatewayId")
        pcx_id = properties.get("VpcPeeringConnectionId")

        route_table_id = properties["RouteTableId"]
        ec2_backend = ec2_backends[account_id][region_name]
        route = ec2_backend.create_route(
            route_table_id=route_table_id,
            destination_cidr_block=properties.get("DestinationCidrBlock"),
            gateway_id=gateway_id,
            instance_id=instance_id,
            nat_gateway_id=nat_gateway_id,
            egress_only_igw_id=egress_only_igw_id,
            transit_gateway_id=transit_gateway_id,
            interface_id=interface_id,
            vpc_peering_connection_id=pcx_id,
        )
        return route


class RouteBackend:
    def __init__(self) -> None:
        self.route_tables: Dict[str, RouteTable] = {}

    def create_route_table(
        self,
        vpc_id: str,
        tags: Optional[List[Dict[str, str]]] = None,
        main: bool = False,
    ) -> RouteTable:
        route_table_id = random_route_table_id()
        vpc = self.get_vpc(vpc_id)  # type: ignore[attr-defined]  # Validate VPC exists
        route_table = RouteTable(self, route_table_id, vpc_id, main=main)
        for tag in tags or []:
            route_table.add_tag(tag["Key"], tag["Value"])
        self.route_tables[route_table_id] = route_table

        # creating default routes for ipv4 cirds
        ipv4_cidrs = vpc.get_cidr_block_association_set(ipv6=False)
        for ipv4_cidr in ipv4_cidrs:
            self.create_route(route_table_id, ipv4_cidr.get("cidr_block"), local=True)

        # creating default routes for ipv6 cidrs
        ipv6_cidrs = vpc.get_cidr_block_association_set(ipv6=True)
        for ipv6_cidr in ipv6_cidrs:
            self.create_route(
                route_table_id,
                destination_cidr_block=None,
                local=True,
                destination_ipv6_cidr_block=ipv6_cidr.get("cidr_block"),
            )

        return route_table

    def get_route_table(self, route_table_id: str) -> RouteTable:
        route_table = self.route_tables.get(route_table_id, None)
        if not route_table:
            raise InvalidRouteTableIdError(route_table_id)
        return route_table

    def describe_route_tables(
        self, route_table_ids: Optional[List[str]] = None, filters: Any = None
    ) -> List[RouteTable]:
        route_tables = list(self.route_tables.values())

        if route_table_ids:
            route_tables = [
                route_table
                for route_table in route_tables
                if route_table.id in route_table_ids
            ]
            if len(route_tables) != len(route_table_ids):
                invalid_id = list(
                    set(route_table_ids).difference(
                        set([route_table.id for route_table in route_tables])
                    )
                )[0]
                raise InvalidRouteTableIdError(invalid_id)

        return generic_filter(filters, route_tables)

    def delete_route_table(self, route_table_id: str) -> None:
        route_table = self.get_route_table(route_table_id)
        if route_table.associations:
            raise DependencyViolationError(
                f"The routeTable '{route_table_id}' has dependencies and cannot be deleted."
            )
        self.route_tables.pop(route_table_id)

    def associate_route_table(
        self,
        route_table_id: str,
        gateway_id: Optional[str] = None,
        subnet_id: Optional[str] = None,
    ) -> str:
        # Idempotent if association already exists.
        route_tables_by_subnet = self.describe_route_tables(
            filters={"association.subnet-id": [subnet_id]}
        )
        if route_tables_by_subnet:
            for association_id, check_subnet_id in route_tables_by_subnet[
                0
            ].associations.items():
                if subnet_id == check_subnet_id:
                    return association_id

        # Association does not yet exist, so create it.
        route_table = self.get_route_table(route_table_id)
        if gateway_id is None:
            self.get_subnet(subnet_id)  # type: ignore[attr-defined]  # Validate subnet exists
            association_id = random_subnet_association_id()
            route_table.associations[association_id] = subnet_id  # type: ignore[assignment]
            return association_id
        else:
            association_id = random_subnet_association_id()
            route_table.associations[association_id] = gateway_id
            return association_id

    def disassociate_route_table(self, association_id: str) -> Optional[str]:
        for route_table in self.route_tables.values():
            if association_id in route_table.associations:
                return route_table.associations.pop(association_id, None)
        raise InvalidAssociationIdError(association_id)

    def replace_route_table_association(
        self, association_id: str, route_table_id: str
    ) -> str:
        # Idempotent if association already exists.
        new_route_table = self.get_route_table(route_table_id)
        if association_id in new_route_table.all_associations_ids:
            return association_id

        # Find route table which currently has the association, error if none.
        route_tables_by_association_id = self.describe_route_tables(
            filters={"association.route-table-association-id": [association_id]}
        )
        if not route_tables_by_association_id:
            raise InvalidAssociationIdError(association_id)
        previous_route_table = route_tables_by_association_id[0]

        # Remove existing association, create new one.
        new_association_id = random_subnet_association_id()
        if previous_route_table.main_association_id == association_id:
            previous_route_table.main_association_id = None
            new_route_table.main_association_id = new_association_id
        else:
            association_target_id = previous_route_table.associations.pop(
                association_id
            )
            new_route_table.associations[new_association_id] = association_target_id
        return new_association_id

    def create_route(
        self,
        route_table_id: str,
        destination_cidr_block: Optional[str],
        destination_ipv6_cidr_block: Optional[str] = None,
        destination_prefix_list_id: Optional[str] = None,
        local: bool = False,
        gateway_id: Optional[str] = None,
        instance_id: Optional[str] = None,
        nat_gateway_id: Optional[str] = None,
        egress_only_igw_id: Optional[str] = None,
        transit_gateway_id: Optional[str] = None,
        interface_id: Optional[str] = None,
        vpc_peering_connection_id: Optional[str] = None,
        carrier_gateway_id: Optional[str] = None,
        vpc_endpoint_id: Optional[str] = None,
    ) -> Route:
        gateway = None
        nat_gateway = None
        transit_gateway = None
        egress_only_igw = None
        interface = None
        destination_prefix_list = None
        carrier_gateway = None

        if vpc_endpoint_id:
            vpce = self.describe_vpc_endpoints(vpc_end_point_ids=[vpc_endpoint_id])  # type: ignore[attr-defined]
            if not vpce[0].endpoint_type == "GatewayLoadBalancer":
                # https://boto3.amazonaws.com/v1/documentation/api/latest/reference/services/ec2/client/create_route.html
                # VpcEndpointId (string) – The ID of a VPC endpoint. Supported for Gateway Load Balancer endpoints only.
                raise RouteNotSupportedError(vpc_endpoint_id)

        route_table = self.get_route_table(route_table_id)

        if interface_id:
            # for validating interface Id whether it is valid or not.
            interface = self.get_network_interface(interface_id)  # type: ignore[attr-defined]

        else:
            if gateway_id:
                if EC2_RESOURCE_TO_PREFIX["vpn-gateway"] in gateway_id:
                    gateway = self.get_vpn_gateway(gateway_id)  # type: ignore[attr-defined]
                elif EC2_RESOURCE_TO_PREFIX["internet-gateway"] in gateway_id:
                    gateway = self.get_internet_gateway(gateway_id)  # type: ignore[attr-defined]
                elif EC2_RESOURCE_TO_PREFIX["vpc-endpoint"] in gateway_id:
                    gateway = self.get_vpc_end_point(gateway_id)  # type: ignore[attr-defined]

            if destination_cidr_block:
                self.__validate_destination_cidr_block(
                    destination_cidr_block, route_table
                )

            if nat_gateway_id is not None:
                nat_gateway = self.nat_gateways.get(nat_gateway_id)  # type: ignore[attr-defined]
            if egress_only_igw_id is not None:
                egress_only_igw = self.get_egress_only_igw(egress_only_igw_id)  # type: ignore[attr-defined]
            if transit_gateway_id is not None:
                transit_gateway = self.transit_gateways.get(transit_gateway_id)  # type: ignore[attr-defined]
            if destination_prefix_list_id is not None:
                destination_prefix_list = self.managed_prefix_lists.get(  # type: ignore[attr-defined]
                    destination_prefix_list_id
                )
            if carrier_gateway_id is not None:
                carrier_gateway = self.carrier_gateways.get(carrier_gateway_id)  # type: ignore[attr-defined]

        route = Route(
            route_table,
            destination_cidr_block,
            destination_ipv6_cidr_block,
            destination_prefix_list,
            local=local,
            gateway=gateway,
            instance=self.get_instance(instance_id) if instance_id else None,  # type: ignore[attr-defined]
            nat_gateway=nat_gateway,
            egress_only_igw=egress_only_igw,
            transit_gateway=transit_gateway,
            interface=interface,
            carrier_gateway=carrier_gateway,
            vpc_endpoint_id=vpc_endpoint_id,
            vpc_pcx=self.get_vpc_peering_connection(vpc_peering_connection_id)  # type: ignore[attr-defined]
            if vpc_peering_connection_id
            else None,
        )
        route_table.routes[route.id] = route
        return route

    def replace_route(
        self,
        route_table_id: str,
        destination_cidr_block: str,
        destination_ipv6_cidr_block: Optional[str] = None,
        destination_prefix_list_id: Optional[str] = None,
        nat_gateway_id: Optional[str] = None,
        egress_only_igw_id: Optional[str] = None,
        transit_gateway_id: Optional[str] = None,
        gateway_id: Optional[str] = None,
        instance_id: Optional[str] = None,
        interface_id: Optional[str] = None,
        vpc_peering_connection_id: Optional[str] = None,
    ) -> Route:
        cidr = destination_cidr_block
        if destination_ipv6_cidr_block:
            cidr = destination_ipv6_cidr_block
        if destination_prefix_list_id:
            cidr = destination_prefix_list_id
        route_table = self.get_route_table(route_table_id)
        route_id = generate_route_id(
            route_table.id, destination_cidr_block, destination_ipv6_cidr_block
        )
        try:
            route = route_table.routes[route_id]
        except KeyError:
            # This should be 'raise InvalidRouteError(route_table_id, cidr)' in
            # line with the delete_route() equivalent, but for some reason AWS
            # returns InvalidParameterValue instead in this case.
            raise InvalidParameterValueErrorReplaceRoute(cidr) from None

        route.gateway = None
        route.nat_gateway = None
        route.egress_only_igw = None
        route.transit_gateway = None
        if gateway_id:
            if EC2_RESOURCE_TO_PREFIX["vpn-gateway"] in gateway_id:
                route.gateway = self.get_vpn_gateway(gateway_id)  # type: ignore[attr-defined]
            elif EC2_RESOURCE_TO_PREFIX["internet-gateway"] in gateway_id:
                route.gateway = self.get_internet_gateway(gateway_id)  # type: ignore[attr-defined]

        if nat_gateway_id is not None:
            route.nat_gateway = self.nat_gateways.get(nat_gateway_id)  # type: ignore[attr-defined]
        if egress_only_igw_id is not None:
            route.egress_only_igw = self.get_egress_only_igw(egress_only_igw_id)  # type: ignore[attr-defined]
        if transit_gateway_id is not None:
            route.transit_gateway = self.transit_gateways.get(transit_gateway_id)  # type: ignore[attr-defined]
        if destination_prefix_list_id is not None:
            route.prefix_list = self.managed_prefix_lists.get(  # type: ignore[attr-defined]
                destination_prefix_list_id
            )

        route.instance = self.get_instance(instance_id) if instance_id else None  # type: ignore[attr-defined]
        route.interface = (
            self.get_network_interface(interface_id) if interface_id else None  # type: ignore[attr-defined]
        )
        route.vpc_pcx = (
            self.get_vpc_peering_connection(vpc_peering_connection_id)  # type: ignore[attr-defined]
            if vpc_peering_connection_id
            else None
        )

        route_table.routes[route.id] = route
        return route

    def delete_route(
        self,
        route_table_id: str,
        destination_cidr_block: str,
        destination_ipv6_cidr_block: Optional[str] = None,
        destination_prefix_list_id: Optional[str] = None,
    ) -> Route:
        cidr = destination_cidr_block
        route_table = self.get_route_table(route_table_id)
        if destination_ipv6_cidr_block:
            cidr = destination_ipv6_cidr_block
        if destination_prefix_list_id:
            cidr = destination_prefix_list_id
        route_id = generate_route_id(route_table_id, cidr)
        deleted = route_table.routes.pop(route_id, None)
        if not deleted:
            raise InvalidRouteError(route_table_id, cidr)
        return deleted

    def __validate_destination_cidr_block(
        self, destination_cidr_block: str, route_table: RouteTable
    ) -> None:
        """
        Utility function to check the destination CIDR block
        Will validate the format and check for overlap with existing routes
        """
        try:
            ip_v4_network = ipaddress.IPv4Network(
                str(destination_cidr_block), strict=False
            )
        except ValueError:
            raise InvalidDestinationCIDRBlockParameterError(destination_cidr_block)

        if not route_table.routes:
            return
        for route in route_table.routes.values():
            if not route.destination_cidr_block:
                continue
            if not route.local and ip_v4_network == ipaddress.IPv4Network(
                str(route.destination_cidr_block)
            ):
                raise RouteAlreadyExistsError(destination_cidr_block)
