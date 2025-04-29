from typing import Any, Dict, List, Optional

from moto.core.common_models import CloudFormationModel

from ..exceptions import (
    DependencyViolationError,
    GatewayNotAttachedError,
    InvalidGatewayIDError,
    InvalidInternetGatewayIdError,
    InvalidVPCIdError,
    ResourceAlreadyAssociatedError,
)
from ..utils import (
    filter_internet_gateways,
    random_egress_only_internet_gateway_id,
    random_internet_gateway_id,
)
from .core import TaggedEC2Resource
from .vpn_gateway import VPCGatewayAttachment


class EgressOnlyInternetGateway(TaggedEC2Resource):
    def __init__(
        self, ec2_backend: Any, vpc_id: str, tags: Optional[Dict[str, str]] = None
    ):
        self.id = random_egress_only_internet_gateway_id()
        self.ec2_backend = ec2_backend
        self.vpc_id = vpc_id
        self.state = "attached"
        self.add_tags(tags or {})

    @property
    def physical_resource_id(self) -> str:
        return self.id


class EgressOnlyInternetGatewayBackend:
    def __init__(self) -> None:
        self.egress_only_internet_gateways: Dict[str, EgressOnlyInternetGateway] = {}

    def create_egress_only_internet_gateway(
        self, vpc_id: str, tags: Optional[Dict[str, str]] = None
    ) -> EgressOnlyInternetGateway:
        vpc = self.get_vpc(vpc_id)  # type: ignore[attr-defined]
        if not vpc:
            raise InvalidVPCIdError(vpc_id)
        egress_only_igw = EgressOnlyInternetGateway(self, vpc_id, tags)
        self.egress_only_internet_gateways[egress_only_igw.id] = egress_only_igw
        return egress_only_igw

    def describe_egress_only_internet_gateways(
        self, ids: Optional[List[str]] = None
    ) -> List[EgressOnlyInternetGateway]:
        """
        The Filters-argument is not yet supported
        """
        egress_only_igws = list(self.egress_only_internet_gateways.values())

        if ids:
            egress_only_igws = [
                egress_only_igw
                for egress_only_igw in egress_only_igws
                if egress_only_igw.id in ids
            ]
        return egress_only_igws

    def delete_egress_only_internet_gateway(self, gateway_id: str) -> None:
        egress_only_igw = self.egress_only_internet_gateways.get(gateway_id)
        if not egress_only_igw:
            raise InvalidGatewayIDError(gateway_id)
        if egress_only_igw:
            self.egress_only_internet_gateways.pop(gateway_id)

    def get_egress_only_igw(self, gateway_id: str) -> EgressOnlyInternetGateway:
        igw = self.egress_only_internet_gateways.get(gateway_id)
        if not igw:
            raise InvalidGatewayIDError(gateway_id)
        return igw


class InternetGateway(TaggedEC2Resource, CloudFormationModel):
    def __init__(self, ec2_backend: Any):
        self.ec2_backend = ec2_backend
        self.id = random_internet_gateway_id()
        self.vpc = None

    @property
    def owner_id(self) -> str:
        return self.ec2_backend.account_id

    @staticmethod
    def cloudformation_name_type() -> str:
        return ""

    @staticmethod
    def cloudformation_type() -> str:
        # https://docs.aws.amazon.com/AWSCloudFormation/latest/UserGuide/aws-resource-ec2-internetgateway.html
        return "AWS::EC2::InternetGateway"

    @classmethod
    def create_from_cloudformation_json(  # type: ignore[misc]
        cls,
        resource_name: str,
        cloudformation_json: Any,
        account_id: str,
        region_name: str,
        **kwargs: Any,
    ) -> "InternetGateway":
        from ..models import ec2_backends

        ec2_backend = ec2_backends[account_id][region_name]
        return ec2_backend.create_internet_gateway()

    @property
    def physical_resource_id(self) -> str:
        return self.id

    @property
    def attachment_state(self) -> str:
        if self.vpc:
            return "available"
        else:
            return "detached"


class InternetGatewayBackend:
    def __init__(self) -> None:
        self.internet_gateways: Dict[str, InternetGateway] = {}

    def create_internet_gateway(
        self, tags: Optional[List[Dict[str, str]]] = None
    ) -> InternetGateway:
        igw = InternetGateway(self)
        for tag in tags or []:
            igw.add_tag(tag["Key"], tag["Value"])
        self.internet_gateways[igw.id] = igw
        return igw

    def describe_internet_gateways(
        self, internet_gateway_ids: Optional[List[str]] = None, filters: Any = None
    ) -> List[InternetGateway]:
        igws = []
        if internet_gateway_ids is None:
            igws = list(self.internet_gateways.values())
        else:
            for igw_id in internet_gateway_ids:
                if igw_id in self.internet_gateways:
                    igws.append(self.internet_gateways[igw_id])
                else:
                    raise InvalidInternetGatewayIdError(igw_id)
        if filters is not None:
            igws = filter_internet_gateways(igws, filters)
        return igws

    def delete_internet_gateway(self, internet_gateway_id: str) -> None:
        igw = self.get_internet_gateway(internet_gateway_id)
        if igw.vpc:
            raise DependencyViolationError(
                f"{internet_gateway_id} is being utilized by {igw.vpc.id}"
            )
        self.internet_gateways.pop(internet_gateway_id)

    def detach_internet_gateway(self, internet_gateway_id: str, vpc_id: str) -> None:
        igw = self.get_internet_gateway(internet_gateway_id)
        if not igw.vpc or igw.vpc.id != vpc_id:
            raise GatewayNotAttachedError(internet_gateway_id, vpc_id)
        igw.vpc = None

    def attach_internet_gateway(
        self, internet_gateway_id: str, vpc_id: str
    ) -> VPCGatewayAttachment:
        igw = self.get_internet_gateway(internet_gateway_id)
        if igw.vpc:
            raise ResourceAlreadyAssociatedError(internet_gateway_id)
        vpc = self.get_vpc(vpc_id)  # type: ignore[attr-defined]
        igw.vpc = vpc
        return VPCGatewayAttachment(gateway_id=internet_gateway_id, vpc_id=vpc_id)

    def get_internet_gateway(self, internet_gateway_id: str) -> InternetGateway:
        igw_ids = [internet_gateway_id]
        return self.describe_internet_gateways(internet_gateway_ids=igw_ids)[0]
