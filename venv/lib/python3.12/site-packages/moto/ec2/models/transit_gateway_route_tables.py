import itertools
from datetime import datetime
from typing import Any, Dict, List, Optional

from moto.core.utils import utcnow
from moto.utilities.utils import filter_resources

from ..utils import random_transit_gateway_route_table_id
from .core import TaggedEC2Resource


class RouteTableAssociation:
    def __init__(self) -> None:
        self.resourceId = ""
        self.resourceType = ""
        self.state = ""
        self.transitGatewayAttachmentId = ""

    def clear(self) -> None:
        self.resourceId = ""
        self.resourceType = ""
        self.state = ""
        self.transitGatewayAttachmentId = ""


class RouteTablePropagation(RouteTableAssociation):
    pass


class TransitGatewayRouteTable(TaggedEC2Resource):
    def __init__(
        self,
        backend: Any,
        transit_gateway_id: str,
        tags: Optional[Dict[str, str]] = None,
        default_association_route_table: bool = False,
        default_propagation_route_table: bool = False,
    ):
        self.ec2_backend = backend
        self.id = random_transit_gateway_route_table_id()
        self.transit_gateway_id = transit_gateway_id

        self._created_at = utcnow()

        self.default_association_route_table = default_association_route_table
        self.default_propagation_route_table = default_propagation_route_table
        self.state = "available"
        self.routes: Dict[str, Dict[str, Optional[str]]] = {}
        self.add_tags(tags or {})
        self.route_table_association = RouteTableAssociation()
        self.route_table_propagation: List[RouteTablePropagation] = []

    @property
    def physical_resource_id(self) -> str:
        return self.id

    @property
    def create_time(self) -> datetime:
        return self._created_at

    @property
    def creation_time(self) -> datetime:
        return self.create_time


class TransitGatewayRelations:
    # this class is for TransitGatewayAssociation and TransitGatewayPropagation
    def __init__(
        self,
        backend: Any,
        transit_gateway_attachment_id: str,
        transit_gateway_route_table_id: str,
        state: str,
    ):
        self.ec2_backend = backend
        self.transit_gateway_attachment_id = transit_gateway_attachment_id
        self.transit_gateway_route_table_id = transit_gateway_route_table_id
        self.resource_id = backend.transit_gateway_attachments[
            transit_gateway_attachment_id
        ].resource_id
        self.resource_type = backend.transit_gateway_attachments[
            transit_gateway_attachment_id
        ].resource_type
        self.state = state


class TransitGatewayRouteTableBackend:
    def __init__(self) -> None:
        self.transit_gateways_route_tables: Dict[str, TransitGatewayRouteTable] = {}
        self.transit_gateway_associations: Dict[str, TransitGatewayRelations] = {}
        self.transit_gateway_propagations: Dict[str, TransitGatewayRelations] = {}

    def create_transit_gateway_route_table(
        self,
        transit_gateway_id: str,
        tags: Optional[Dict[str, str]] = None,
        default_association_route_table: bool = False,
        default_propagation_route_table: bool = False,
    ) -> TransitGatewayRouteTable:
        transit_gateways_route_table = TransitGatewayRouteTable(
            self,
            transit_gateway_id=transit_gateway_id,
            tags=tags,
            default_association_route_table=default_association_route_table,
            default_propagation_route_table=default_propagation_route_table,
        )
        self.transit_gateways_route_tables[transit_gateways_route_table.id] = (
            transit_gateways_route_table
        )
        return transit_gateways_route_table

    def get_all_transit_gateway_route_tables(
        self, transit_gateway_route_table_ids: Optional[str] = None, filters: Any = None
    ) -> List[TransitGatewayRouteTable]:
        transit_gateway_route_tables = list(self.transit_gateways_route_tables.values())

        attr_pairs = (
            ("default-association-route-table", "default_association_route_table"),
            ("default-propagation-route-table", "default_propagation_route_table"),
            ("state", "state"),
            ("transit-gateway-id", "transit_gateway_id"),
            ("transit-gateway-route-table-id", "id"),
        )

        if transit_gateway_route_table_ids:
            transit_gateway_route_tables = [
                transit_gateway_route_table
                for transit_gateway_route_table in transit_gateway_route_tables
                if transit_gateway_route_table.id in transit_gateway_route_table_ids
            ]

        result = transit_gateway_route_tables
        if filters:
            result = filter_resources(transit_gateway_route_tables, filters, attr_pairs)
        return result

    def delete_transit_gateway_route_table(
        self, transit_gateway_route_table_id: str
    ) -> TransitGatewayRouteTable:
        transit_gateways_route_table = self.transit_gateways_route_tables[
            transit_gateway_route_table_id
        ]
        transit_gateways_route_table.state = "deleted"
        return transit_gateways_route_table

    def create_transit_gateway_route(
        self,
        transit_gateway_route_table_id: str,
        destination_cidr_block: str,
        transit_gateway_attachment_id: Optional[str] = None,
        blackhole: bool = False,
    ) -> Dict[str, Optional[str]]:
        transit_gateways_route_table = self.transit_gateways_route_tables[
            transit_gateway_route_table_id
        ]
        transit_gateway_attachment = self.transit_gateway_attachments.get(  # type: ignore[attr-defined]
            transit_gateway_attachment_id
        )
        transit_gateways_route_table.routes[destination_cidr_block] = {
            "destinationCidrBlock": destination_cidr_block,
            "prefixListId": None,
            "state": "blackhole" if blackhole else "active",
            "type": "static",
        }

        if transit_gateway_attachment:
            transit_gateway_attachment_dict = {
                "transitGatewayAttachments": {
                    "resourceId": transit_gateway_attachment.resource_id,
                    "resourceType": transit_gateway_attachment.resource_type,
                    "transitGatewayAttachmentId": transit_gateway_attachment_id,
                }
            }
            transit_gateways_route_table.routes[destination_cidr_block].update(
                transit_gateway_attachment_dict  # type: ignore
            )
        return transit_gateways_route_table.routes[destination_cidr_block]

    def delete_transit_gateway_route(
        self, transit_gateway_route_table_id: str, destination_cidr_block: str
    ) -> dict[str, Optional[str]]:
        transit_gateways_route_table = self.transit_gateways_route_tables[
            transit_gateway_route_table_id
        ]
        transit_gateways_route_table.routes[destination_cidr_block]["state"] = "deleted"
        return transit_gateways_route_table.routes.pop(destination_cidr_block)

    def search_transit_gateway_routes(
        self,
        transit_gateway_route_table_id: str,
        filters: Any,
        max_results: Optional[int] = None,
    ) -> list[Dict[str, Optional[str]]]:
        """
        The following filters are currently supported: type, state, route-search.exact-match
        """
        transit_gateway_route_table = self.transit_gateways_route_tables.get(
            transit_gateway_route_table_id
        )
        if not transit_gateway_route_table:
            return []

        attr_pairs = (
            ("type", "type"),
            ("state", "state"),
            ("route-search.exact-match", "destinationCidrBlock"),
        )

        routes = transit_gateway_route_table.routes.copy()
        for key in transit_gateway_route_table.routes:
            for attrs in attr_pairs:
                values = filters.get(attrs[0]) or None
                if values:
                    if routes[key].get(attrs[1]) not in values:
                        routes.pop(key)
                        break
        if max_results:
            routes = routes[: int(max_results)]  # type: ignore
        return list(routes.values())

    def set_route_table_association(
        self, transit_gateway_attachment_id: str, transit_gateway_route_table_id: str
    ) -> None:
        table = self.transit_gateways_route_tables[transit_gateway_route_table_id]
        attachment = self.transit_gateway_attachments[transit_gateway_attachment_id]  # type: ignore[attr-defined]
        table.route_table_association.resourceId = attachment.resource_id
        table.route_table_association.resourceType = attachment.resource_type
        table.route_table_association.state = "associated"
        table.route_table_association.transitGatewayAttachmentId = (
            transit_gateway_attachment_id
        )

    def unset_route_table_association(self, tgw_rt_id: str) -> None:
        tgw_rt = self.transit_gateways_route_tables[tgw_rt_id]
        tgw_rt.route_table_association.clear()

    def set_route_table_propagation(
        self, transit_gateway_attachment_id: str, transit_gateway_route_table_id: str
    ) -> None:
        route_table = self.transit_gateways_route_tables[transit_gateway_route_table_id]
        attchment = self.transit_gateway_attachments[transit_gateway_attachment_id]  # type: ignore[attr-defined]
        propagation = RouteTablePropagation()
        propagation.resourceId = attchment.resource_id
        propagation.resourceType = attchment.resource_type
        propagation.state = "enabled"
        propagation.transitGatewayAttachmentId = transit_gateway_attachment_id
        route_table.route_table_propagation.append(propagation)

    def disable_route_table_propagation(
        self, transit_gateway_attachment_id: str, transit_gateway_route_table_id: str
    ) -> None:
        route_table = self.transit_gateways_route_tables[transit_gateway_route_table_id]
        route_table.route_table_propagation = [
            prop
            for prop in route_table.route_table_propagation
            if prop.transitGatewayAttachmentId != transit_gateway_attachment_id
        ]

    def get_transit_gateway_route_table_associations(
        self, transit_gateway_route_table_id: List[str], filters: Any = None
    ) -> list[RouteTableAssociation]:
        transit_gateway_route_tables = list(self.transit_gateways_route_tables.values())

        if transit_gateway_route_tables:
            transit_gateway_route_tables = [
                transit_gateway_route_table
                for transit_gateway_route_table in transit_gateway_route_tables
                if transit_gateway_route_table.id in transit_gateway_route_table_id
            ]

        attr_pairs = (
            ("resource-id", "resourceId"),
            ("resource-type", "resourceType"),
            (
                "transit-gateway-attachment-id",
                "transitGatewayAttachmentId",
            ),
        )

        result = [
            table.route_table_association for table in transit_gateway_route_tables
        ]
        if filters:
            result = filter_resources(result, filters, attr_pairs)
        return result

    def get_transit_gateway_route_table_propagations(
        self, transit_gateway_route_table_id: str, filters: Any = None
    ) -> List[RouteTablePropagation]:
        transit_gateway_route_tables = list(self.transit_gateways_route_tables.values())

        if transit_gateway_route_tables:
            transit_gateway_route_tables = [
                transit_gateway_route_table
                for transit_gateway_route_table in transit_gateway_route_tables
                if transit_gateway_route_table.id in transit_gateway_route_table_id
            ]

        attr_pairs = (
            ("resource-id", "resourceId"),
            ("resource-type", "resourceType"),
            (
                "transit-gateway-attachment-id",
                "transitGatewayAttachmentId",
            ),
        )

        result = list(
            itertools.chain.from_iterable(
                table.route_table_propagation for table in transit_gateway_route_tables
            )
        )
        if filters:
            result = filter_resources(result, filters, attr_pairs)
        return result

    def associate_transit_gateway_route_table(
        self, transit_gateway_attachment_id: str, transit_gateway_route_table_id: str
    ) -> TransitGatewayRelations:
        transit_gateway_association = TransitGatewayRelations(
            self,
            transit_gateway_attachment_id,
            transit_gateway_route_table_id,
            state="associated",
        )
        self.set_route_table_association(
            transit_gateway_attachment_id, transit_gateway_route_table_id
        )
        self.set_attachment_association(  # type: ignore[attr-defined]
            transit_gateway_attachment_id, transit_gateway_route_table_id
        )
        self.transit_gateway_associations[transit_gateway_attachment_id] = (
            transit_gateway_association
        )

        return transit_gateway_association

    def enable_transit_gateway_route_table_propagation(
        self, transit_gateway_attachment_id: str, transit_gateway_route_table_id: str
    ) -> TransitGatewayRelations:
        transit_gateway_propagation = TransitGatewayRelations(
            self,
            transit_gateway_attachment_id,
            transit_gateway_route_table_id,
            state="enabled",
        )
        self.set_route_table_propagation(
            transit_gateway_attachment_id, transit_gateway_route_table_id
        )
        self.set_attachment_propagation(  # type: ignore[attr-defined]
            transit_gateway_attachment_id, transit_gateway_route_table_id
        )
        self.transit_gateway_propagations[transit_gateway_attachment_id] = (
            transit_gateway_propagation
        )

        return transit_gateway_propagation

    def disable_transit_gateway_route_table_propagation(
        self, transit_gateway_attachment_id: str, transit_gateway_route_table_id: str
    ) -> TransitGatewayRelations:
        self.disable_route_table_propagation(
            transit_gateway_attachment_id=transit_gateway_attachment_id,
            transit_gateway_route_table_id=transit_gateway_route_table_id,
        )
        self.disable_attachment_propagation(  # type: ignore[attr-defined]
            transit_gateway_attachment_id=transit_gateway_attachment_id
        )
        self.transit_gateway_propagations[
            transit_gateway_attachment_id
        ].state = "disabled"
        transit_gateway_propagation = self.transit_gateway_propagations.pop(
            transit_gateway_attachment_id
        )

        return transit_gateway_propagation

    def disassociate_transit_gateway_route_table(
        self, tgw_attach_id: str, tgw_rt_id: str
    ) -> TransitGatewayRelations:
        tgw_association = self.transit_gateway_associations.pop(tgw_attach_id)
        tgw_association.state = "disassociated"

        self.unset_route_table_association(tgw_rt_id)
        self.unset_attachment_association(tgw_attach_id)  # type: ignore[attr-defined]

        return tgw_association
