import weakref
from collections import defaultdict
from typing import Any, Dict, Iterator, List, Optional

from moto.core.utils import iso_8601_datetime_with_milliseconds, utcnow
from moto.utilities.utils import filter_resources, merge_multiple_dicts

from ..exceptions import InvalidParameterValueErrorPeeringAttachment
from ..utils import describe_tag_filter, random_transit_gateway_attachment_id
from .core import TaggedEC2Resource
from .vpc_peering_connections import PeeringConnectionStatus


class TransitGatewayAttachment(TaggedEC2Resource):
    def __init__(
        self,
        backend: Any,
        resource_id: str,
        resource_type: str,
        transit_gateway_id: str,
        tags: Optional[Dict[str, str]] = None,
    ):
        self.ec2_backend = backend
        self.association: Dict[str, str] = {}
        self.propagation: Dict[str, str] = {}
        self.resource_id = resource_id
        self.resource_type = resource_type

        self.id = random_transit_gateway_attachment_id()
        self.transit_gateway_id = transit_gateway_id

        self.state = "available"
        self.add_tags(tags or {})

        self._created_at = utcnow()
        self.resource_owner_id = backend.account_id
        self.transit_gateway_owner_id = backend.account_id
        self.owner_id = backend.account_id

    @property
    def create_time(self) -> str:
        return iso_8601_datetime_with_milliseconds(self._created_at)


class TransitGatewayVpcAttachment(TransitGatewayAttachment):
    DEFAULT_OPTIONS = {
        "ApplianceModeSupport": "disable",
        "DnsSupport": "enable",
        "Ipv6Support": "disable",
    }

    def __init__(
        self,
        backend: Any,
        transit_gateway_id: str,
        vpc_id: str,
        subnet_ids: List[str],
        tags: Optional[Dict[str, str]] = None,
        options: Optional[Dict[str, str]] = None,
    ):
        super().__init__(
            backend=backend,
            transit_gateway_id=transit_gateway_id,
            resource_id=vpc_id,
            resource_type="vpc",
            tags=tags,
        )

        self.vpc_id = vpc_id
        self.subnet_ids = subnet_ids
        self.options = merge_multiple_dicts(self.DEFAULT_OPTIONS, options or {})


class TransitGatewayPeeringAttachment(TransitGatewayAttachment):
    def __init__(
        self,
        backend: Any,
        transit_gateway_id: str,
        peer_transit_gateway_id: str,
        peer_region: str,
        peer_account_id: str,
        tags: Dict[str, str],
        region_name: str,
    ):
        super().__init__(
            backend=backend,
            transit_gateway_id=transit_gateway_id,
            resource_id=peer_transit_gateway_id,
            resource_type="peering",
            tags=tags,
        )

        self.accepter_tgw_info = {
            "ownerId": peer_account_id,
            "region": peer_region,
            "transitGatewayId": peer_transit_gateway_id,
        }
        self.requester_tgw_info = {
            "ownerId": self.owner_id,
            "region": region_name,
            "transitGatewayId": transit_gateway_id,
        }
        self.status = PeeringConnectionStatus(accepter_id=peer_account_id)


class TransitGatewayAttachmentBackend:
    backend_refs = defaultdict(set)  # type: ignore

    def __init__(self) -> None:
        self.transit_gateway_attachments: Dict[str, TransitGatewayAttachment] = {}
        self.backend_refs[self.__class__].add(weakref.ref(self))

    @classmethod
    def _get_peering_attachment_backend_refs(
        cls,
    ) -> Iterator["TransitGatewayAttachmentBackend"]:
        for backend_ref in cls.backend_refs[cls]:
            backend = backend_ref()
            if backend is not None:
                yield backend

    def create_transit_gateway_vpn_attachment(
        self,
        vpn_id: str,
        transit_gateway_id: str,
        tags: Optional[Dict[str, str]] = None,
    ) -> TransitGatewayAttachment:
        transit_gateway_vpn_attachment = TransitGatewayAttachment(
            self,
            resource_id=vpn_id,
            resource_type="vpn",
            transit_gateway_id=transit_gateway_id,
            tags=tags,
        )
        self.transit_gateway_attachments[transit_gateway_vpn_attachment.id] = (
            transit_gateway_vpn_attachment
        )
        return transit_gateway_vpn_attachment

    def create_transit_gateway_vpc_attachment(
        self,
        transit_gateway_id: str,
        vpc_id: str,
        subnet_ids: List[str],
        tags: Optional[Dict[str, str]] = None,
        options: Optional[Dict[str, str]] = None,
    ) -> TransitGatewayVpcAttachment:
        transit_gateway_vpc_attachment = TransitGatewayVpcAttachment(
            self,
            transit_gateway_id=transit_gateway_id,
            tags=tags,
            vpc_id=vpc_id,
            subnet_ids=subnet_ids,
            options=options,
        )
        self.transit_gateway_attachments[transit_gateway_vpc_attachment.id] = (
            transit_gateway_vpc_attachment
        )
        return transit_gateway_vpc_attachment

    def describe_transit_gateway_attachments(
        self,
        transit_gateways_attachment_ids: Optional[List[str]] = None,
        filters: Any = None,
    ) -> List[TransitGatewayAttachment]:
        transit_gateway_attachments = list(self.transit_gateway_attachments.values())

        attr_pairs = (
            ("resource-id", "resource_id"),
            ("resource-type", "resource_type"),
            ("transit-gateway-id", "transit_gateway_id"),
        )

        if (
            not transit_gateways_attachment_ids == []
            and transit_gateways_attachment_ids is not None
        ):
            transit_gateway_attachments = [
                transit_gateways_attachment
                for transit_gateways_attachment in transit_gateway_attachments
                if transit_gateways_attachment.id in transit_gateways_attachment_ids
            ]

        result = transit_gateway_attachments
        if filters:
            result = filter_resources(transit_gateway_attachments, filters, attr_pairs)
        return result

    def describe_transit_gateway_vpc_attachments(
        self,
        transit_gateways_attachment_ids: Optional[List[str]] = None,
        filters: Any = None,
    ) -> List[TransitGatewayAttachment]:
        transit_gateway_attachments = list(self.transit_gateway_attachments.values())

        attr_pairs = (
            ("state", "state"),
            ("transit-gateway-attachment-id", "id"),
            ("transit-gateway-id", "transit_gateway_id"),
            ("vpc-id", "resource_id"),
        )

        if (
            not transit_gateways_attachment_ids == []
            and transit_gateways_attachment_ids is not None
        ):
            transit_gateway_attachments = [
                transit_gateways_attachment
                for transit_gateways_attachment in transit_gateway_attachments
                if transit_gateways_attachment.id in transit_gateways_attachment_ids
            ]

        result = transit_gateway_attachments
        if filters:
            result = filter_resources(transit_gateway_attachments, filters, attr_pairs)
        return result

    def delete_transit_gateway_vpc_attachment(
        self, transit_gateway_attachment_id: str
    ) -> TransitGatewayAttachment:
        transit_gateway_attachment = self.transit_gateway_attachments.pop(
            transit_gateway_attachment_id
        )
        transit_gateway_attachment.state = "deleted"
        try:
            route_table_id = transit_gateway_attachment.propagation.get(
                "transitGatewayRouteTableId"
            )
            route_table = self.transit_gateways_route_tables[route_table_id]  # type: ignore[attr-defined]
            route_table.route_table_propagation = [
                prop
                for prop in route_table.route_table_propagation
                if prop["transitGatewayAttachmentId"] != transit_gateway_attachment_id
            ]
        except (AttributeError, KeyError, IndexError):
            pass
        return transit_gateway_attachment

    def modify_transit_gateway_vpc_attachment(
        self,
        transit_gateway_attachment_id: str,
        add_subnet_ids: Optional[List[str]] = None,
        options: Optional[Dict[str, str]] = None,
        remove_subnet_ids: Optional[List[str]] = None,
    ) -> TransitGatewayAttachment:
        tgw_attachment = self.transit_gateway_attachments[transit_gateway_attachment_id]
        if remove_subnet_ids:
            tgw_attachment.subnet_ids = [  # type: ignore[attr-defined]
                id
                for id in tgw_attachment.subnet_ids  # type: ignore[attr-defined]
                if id not in remove_subnet_ids
            ]

        if options:
            tgw_attachment.options.update(options)  # type: ignore[attr-defined]

        if add_subnet_ids:
            for subnet_id in add_subnet_ids:
                tgw_attachment.subnet_ids.append(subnet_id)  # type: ignore[attr-defined]

        return tgw_attachment

    def set_attachment_association(
        self, transit_gateway_attachment_id: str, transit_gateway_route_table_id: str
    ) -> None:
        self.transit_gateway_attachments[transit_gateway_attachment_id].association = {
            "state": "associated",
            "transitGatewayRouteTableId": transit_gateway_route_table_id,
        }

    def unset_attachment_association(self, tgw_attach_id: str) -> None:
        self.transit_gateway_attachments[tgw_attach_id].association = {}

    def set_attachment_propagation(
        self, transit_gateway_attachment_id: str, transit_gateway_route_table_id: str
    ) -> None:
        self.transit_gateway_attachments[transit_gateway_attachment_id].propagation = {
            "state": "enabled",
            "transitGatewayRouteTableId": transit_gateway_route_table_id,
        }

    def unset_attachment_propagation(self, tgw_attach_id: str) -> None:
        self.transit_gateway_attachments[tgw_attach_id].propagation = {}

    def disable_attachment_propagation(
        self, transit_gateway_attachment_id: str
    ) -> None:
        self.transit_gateway_attachments[transit_gateway_attachment_id].propagation[
            "state"
        ] = "disabled"

    def create_transit_gateway_peering_attachment(
        self,
        transit_gateway_id: str,
        peer_transit_gateway_id: str,
        peer_region: str,
        peer_account_id: str,
        tags: Dict[str, str],
    ) -> TransitGatewayPeeringAttachment:
        transit_gateway_peering_attachment = TransitGatewayPeeringAttachment(
            self,
            transit_gateway_id=transit_gateway_id,
            peer_transit_gateway_id=peer_transit_gateway_id,
            peer_region=peer_region,
            peer_account_id=peer_account_id,
            tags=tags,
            region_name=self.region_name,  # type: ignore[attr-defined]
        )

        self.transit_gateway_attachments[transit_gateway_peering_attachment.id] = (
            transit_gateway_peering_attachment
        )

        # If the peer is not same as the current account or region, create attachment in peer backend
        if self.account_id != peer_account_id or self.region_name != peer_region:  # type: ignore[attr-defined]
            for backend in self._get_peering_attachment_backend_refs():
                if (
                    backend.account_id == peer_account_id  # type: ignore[attr-defined]
                    and backend.region_name == peer_region  # type: ignore[attr-defined]
                ):
                    backend.transit_gateway_attachments[
                        transit_gateway_peering_attachment.id
                    ] = transit_gateway_peering_attachment

        transit_gateway_peering_attachment.status.pending()
        transit_gateway_peering_attachment.state = "pendingAcceptance"
        return transit_gateway_peering_attachment

    def describe_transit_gateway_peering_attachments(
        self,
        transit_gateways_attachment_ids: Optional[List[str]] = None,
        filters: Any = None,
    ) -> List[TransitGatewayAttachment]:
        transit_gateway_attachments = list(self.transit_gateway_attachments.values())

        attr_pairs = (
            ("state", "state"),
            ("transit-gateway-attachment-id", "id"),
            ("local-owner-id", "requester_tgw_info", "ownerId"),
            ("remote-owner-id", "accepter_tgw_info", "ownerId"),
        )

        if transit_gateways_attachment_ids:
            transit_gateway_attachments = [
                transit_gateways_attachment
                for transit_gateways_attachment in transit_gateway_attachments
                if transit_gateways_attachment.id in transit_gateways_attachment_ids
            ]

        if filters:
            transit_gateway_attachments = filter_resources(
                transit_gateway_attachments, filters, attr_pairs
            )
            transit_gateway_attachments = describe_tag_filter(
                filters, transit_gateway_attachments
            )
        return transit_gateway_attachments

    def accept_transit_gateway_peering_attachment(
        self, transit_gateway_attachment_id: str
    ) -> TransitGatewayAttachment:
        transit_gateway_attachment = self.transit_gateway_attachments[
            transit_gateway_attachment_id
        ]

        requester_account_id = transit_gateway_attachment.requester_tgw_info["ownerId"]  # type: ignore[attr-defined]
        requester_region_name = transit_gateway_attachment.requester_tgw_info["region"]  # type: ignore[attr-defined]
        accepter_account_id = transit_gateway_attachment.accepter_tgw_info["ownerId"]  # type: ignore[attr-defined]
        accepter_region_name = transit_gateway_attachment.accepter_tgw_info["region"]  # type: ignore[attr-defined]

        # For cross-account peering, must be accepted by the accepter
        if (
            requester_account_id != accepter_account_id
            and self.account_id != accepter_account_id  # type: ignore[attr-defined]
        ):
            raise InvalidParameterValueErrorPeeringAttachment(
                "accept", transit_gateway_attachment_id
            )

        if (
            requester_region_name != accepter_region_name
            and self.region_name != accepter_region_name  # type: ignore[attr-defined]
        ):
            raise InvalidParameterValueErrorPeeringAttachment(
                "accept", transit_gateway_attachment_id
            )

        transit_gateway_attachment.state = "available"
        # Bit dodgy - we just assume that we act on a TransitGatewayPeeringAttachment
        # We could just as easily have another sub-class of TransitGatewayAttachment on our hands, which does not have a status-attribute
        transit_gateway_attachment.status.accept()  # type: ignore[attr-defined]
        return transit_gateway_attachment

    def reject_transit_gateway_peering_attachment(
        self, transit_gateway_attachment_id: str
    ) -> TransitGatewayAttachment:
        transit_gateway_attachment = self.transit_gateway_attachments[
            transit_gateway_attachment_id
        ]

        requester_account_id = transit_gateway_attachment.requester_tgw_info["ownerId"]  # type: ignore[attr-defined]
        requester_region_name = transit_gateway_attachment.requester_tgw_info["region"]  # type: ignore[attr-defined]
        accepter_account_id = transit_gateway_attachment.accepter_tgw_info["ownerId"]  # type: ignore[attr-defined]
        accepter_region_name = transit_gateway_attachment.requester_tgw_info["region"]  # type: ignore[attr-defined]

        if (
            requester_account_id != accepter_account_id
            and self.account_id != accepter_account_id  # type: ignore[attr-defined]
        ):
            raise InvalidParameterValueErrorPeeringAttachment(
                "reject", transit_gateway_attachment_id
            )

        if (
            requester_region_name != accepter_region_name
            and self.region_name != accepter_region_name  # type: ignore[attr-defined]
        ):
            raise InvalidParameterValueErrorPeeringAttachment(
                "reject", transit_gateway_attachment_id
            )

        transit_gateway_attachment.state = "rejected"
        transit_gateway_attachment.status.reject()  # type: ignore[attr-defined]
        return transit_gateway_attachment

    def delete_transit_gateway_peering_attachment(
        self, transit_gateway_attachment_id: str
    ) -> TransitGatewayAttachment:
        transit_gateway_attachment = self.transit_gateway_attachments[
            transit_gateway_attachment_id
        ]
        transit_gateway_attachment.state = "deleted"
        transit_gateway_attachment.status.deleted(deleter_id=self.account_id)  # type: ignore[attr-defined]
        return transit_gateway_attachment
