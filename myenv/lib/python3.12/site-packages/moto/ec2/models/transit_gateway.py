from typing import Any, Dict, List, Optional

from moto.core.common_models import CloudFormationModel
from moto.core.utils import iso_8601_datetime_with_milliseconds, utcnow
from moto.utilities.utils import filter_resources, merge_multiple_dicts

from ..utils import (
    describe_tag_filter,
    random_transit_gateway_id,
)
from .core import TaggedEC2Resource


class TransitGateway(TaggedEC2Resource, CloudFormationModel):
    DEFAULT_OPTIONS = {
        "AmazonSideAsn": "64512",
        "AssociationDefaultRouteTableId": "tgw-rtb-0d571391e50cf8514",
        "AutoAcceptSharedAttachments": "disable",
        "DefaultRouteTableAssociation": "enable",
        "DefaultRouteTablePropagation": "enable",
        "DnsSupport": "enable",
        "MulticastSupport": "disable",
        "PropagationDefaultRouteTableId": "tgw-rtb-0d571391e50cf8514",
        "TransitGatewayCidrBlocks": None,
        "VpnEcmpSupport": "enable",
    }

    def __init__(
        self,
        backend: Any,
        description: Optional[str],
        options: Optional[Dict[str, str]] = None,
    ):
        self.ec2_backend = backend
        self.id = random_transit_gateway_id()
        self.description = description
        self.state = "available"
        self.options = merge_multiple_dicts(self.DEFAULT_OPTIONS, options or {})
        self._created_at = utcnow()

    @property
    def physical_resource_id(self) -> str:
        return self.id

    @property
    def create_time(self) -> str:
        return iso_8601_datetime_with_milliseconds(self._created_at)

    @property
    def owner_id(self) -> str:
        return self.ec2_backend.account_id

    @property
    def arn(self) -> str:
        return f"arn:{self.ec2_backend.partition}:ec2:{self.ec2_backend.region_name}:{self.ec2_backend.account_id}:transit-gateway/{self.id}"

    @staticmethod
    def cloudformation_name_type() -> str:
        return ""

    @staticmethod
    def cloudformation_type() -> str:
        # https://docs.aws.amazon.com/AWSCloudFormation/latest/UserGuide/aws-resource-ec2-natgateway.html
        return "AWS::EC2::TransitGateway"

    @classmethod
    def create_from_cloudformation_json(  # type: ignore[misc]
        cls,
        resource_name: str,
        cloudformation_json: Any,
        account_id: str,
        region_name: str,
        **kwargs: Any,
    ) -> "TransitGateway":
        from ..models import ec2_backends

        ec2_backend = ec2_backends[account_id][region_name]
        properties = cloudformation_json["Properties"]
        description = properties["Description"]
        options = dict(properties)
        del options["Description"]
        transit_gateway = ec2_backend.create_transit_gateway(
            description=description, options=options
        )

        for tag in properties.get("Tags", []):
            tag_key = tag["Key"]
            tag_value = tag["Value"]
            transit_gateway.add_tag(tag_key, tag_value)

        return transit_gateway


class TransitGatewayBackend:
    def __init__(self) -> None:
        self.transit_gateways: Dict[str, TransitGateway] = {}

    def create_transit_gateway(
        self,
        description: Optional[str],
        options: Optional[Dict[str, str]] = None,
        tags: Optional[List[Dict[str, str]]] = None,
    ) -> TransitGateway:
        transit_gateway = TransitGateway(self, description, options)
        for tag in tags or []:
            transit_gateway.add_tag(tag["Key"], tag["Value"])

        self.transit_gateways[transit_gateway.id] = transit_gateway
        return transit_gateway

    def describe_transit_gateways(
        self, filters: Any, transit_gateway_ids: Optional[List[str]]
    ) -> List[TransitGateway]:
        transit_gateways = list(self.transit_gateways.values())

        if transit_gateway_ids:
            transit_gateways = [
                item for item in transit_gateways if item.id in transit_gateway_ids
            ]

        attr_pairs = (
            ("transit-gateway-id", "id"),
            ("state", "state"),
            ("owner-id", "owner_id"),
        )

        result = transit_gateways
        if filters:
            result = filter_resources(transit_gateways, filters, attr_pairs)
            result = describe_tag_filter(filters, result)

        return result

    def delete_transit_gateway(self, transit_gateway_id: str) -> TransitGateway:
        return self.transit_gateways.pop(transit_gateway_id)

    def modify_transit_gateway(
        self,
        transit_gateway_id: str,
        description: Optional[str] = None,
        options: Optional[Dict[str, str]] = None,
    ) -> TransitGateway:
        transit_gateway = self.transit_gateways[transit_gateway_id]
        if description:
            transit_gateway.description = description
        if options:
            transit_gateway.options.update(options)
        return transit_gateway
