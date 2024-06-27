import weakref
from collections import defaultdict
from typing import Any, Dict, Iterator, List, Optional

from moto.core.common_models import CloudFormationModel

from ..exceptions import (
    InvalidVPCPeeringConnectionIdError,
    InvalidVPCPeeringConnectionStateTransitionError,
    OperationNotPermitted2,
    OperationNotPermitted3,
    OperationNotPermitted5,
)
from ..utils import random_vpc_peering_connection_id
from .core import TaggedEC2Resource
from .vpcs import VPC


class PeeringConnectionStatus:
    def __init__(
        self, accepter_id: str, code: str = "initiating-request", message: str = ""
    ):
        self.accepter_id = accepter_id
        self.code = code
        self.message = message

    def deleted(self, deleter_id: str) -> None:
        self.code = "deleted"
        self.message = f"Deleted by {deleter_id}"

    def initiating(self) -> None:
        self.code = "initiating-request"
        self.message = f"Initiating Request to {self.accepter_id}"

    def pending(self) -> None:
        self.code = "pending-acceptance"
        self.message = f"Pending Acceptance by {self.accepter_id}"

    def accept(self) -> None:
        self.code = "active"
        self.message = "Active"

    def reject(self) -> None:
        self.code = "rejected"
        self.message = "Inactive"


class VPCPeeringConnection(TaggedEC2Resource, CloudFormationModel):
    DEFAULT_OPTIONS = {
        "AllowEgressFromLocalClassicLinkToRemoteVpc": "false",
        "AllowEgressFromLocalVpcToRemoteClassicLink": "false",
        "AllowDnsResolutionFromRemoteVpc": "false",
    }

    def __init__(
        self,
        backend: Any,
        vpc_pcx_id: str,
        vpc: VPC,
        peer_vpc: VPC,
        tags: Optional[Dict[str, str]] = None,
    ):
        self.id = vpc_pcx_id
        self.ec2_backend = backend
        self.vpc = vpc
        self.peer_vpc = peer_vpc
        self.requester_options = self.DEFAULT_OPTIONS.copy()
        self.accepter_options = self.DEFAULT_OPTIONS.copy()
        self.add_tags(tags or {})
        self._status = PeeringConnectionStatus(accepter_id=peer_vpc.owner_id)

    @staticmethod
    def cloudformation_name_type() -> str:
        return ""

    @staticmethod
    def cloudformation_type() -> str:
        # https://docs.aws.amazon.com/AWSCloudFormation/latest/UserGuide/aws-resource-ec2-vpcpeeringconnection.html
        return "AWS::EC2::VPCPeeringConnection"

    @classmethod
    def create_from_cloudformation_json(  # type: ignore[misc]
        cls,
        resource_name: str,
        cloudformation_json: Any,
        account_id: str,
        region_name: str,
        **kwargs: Any,
    ) -> "VPCPeeringConnection":
        from ..models import ec2_backends

        properties = cloudformation_json["Properties"]

        ec2_backend = ec2_backends[account_id][region_name]
        vpc = ec2_backend.get_vpc(properties["VpcId"])
        peer_vpc = ec2_backend.get_vpc(properties["PeerVpcId"])

        vpc_pcx = ec2_backend.create_vpc_peering_connection(vpc, peer_vpc)

        return vpc_pcx

    @property
    def physical_resource_id(self) -> str:
        return self.id


class VPCPeeringConnectionBackend:
    # for cross region vpc reference
    vpc_pcx_refs = defaultdict(set)  # type: ignore

    def __init__(self) -> None:
        self.vpc_pcxs: Dict[str, VPCPeeringConnection] = {}
        self.vpc_pcx_refs[self.__class__].add(weakref.ref(self))

    @classmethod
    def get_vpc_pcx_refs(cls) -> Iterator[VPCPeeringConnection]:
        for inst_ref in cls.vpc_pcx_refs[cls]:
            inst = inst_ref()
            if inst is not None:
                yield inst

    def create_vpc_peering_connection(
        self, vpc: VPC, peer_vpc: VPC, tags: Optional[Dict[str, str]] = None
    ) -> VPCPeeringConnection:
        vpc_pcx_id = random_vpc_peering_connection_id()
        vpc_pcx = VPCPeeringConnection(self, vpc_pcx_id, vpc, peer_vpc, tags)
        vpc_pcx._status.pending()
        self.vpc_pcxs[vpc_pcx_id] = vpc_pcx
        # insert cross-account/cross-region peering info
        if vpc.owner_id != peer_vpc.owner_id or vpc.region != peer_vpc.region:
            for backend in peer_vpc.ec2_backend.get_vpc_pcx_refs():
                if (
                    backend.account_id == peer_vpc.owner_id
                    and backend.region_name == peer_vpc.region
                ):
                    backend.vpc_pcxs[vpc_pcx_id] = vpc_pcx

        return vpc_pcx

    def describe_vpc_peering_connections(
        self, vpc_peering_ids: Optional[List[str]] = None
    ) -> List[VPCPeeringConnection]:
        all_pcxs = list(self.vpc_pcxs.values())
        if vpc_peering_ids:
            return [pcx for pcx in all_pcxs if pcx.id in vpc_peering_ids]
        return all_pcxs

    def get_vpc_peering_connection(self, vpc_pcx_id: str) -> VPCPeeringConnection:
        if vpc_pcx_id not in self.vpc_pcxs:
            raise InvalidVPCPeeringConnectionIdError(vpc_pcx_id)
        return self.vpc_pcxs[vpc_pcx_id]

    def delete_vpc_peering_connection(self, vpc_pcx_id: str) -> VPCPeeringConnection:
        deleted = self.get_vpc_peering_connection(vpc_pcx_id)
        deleted._status.deleted(deleter_id=self.account_id)  # type: ignore[attr-defined]
        return deleted

    def accept_vpc_peering_connection(self, vpc_pcx_id: str) -> VPCPeeringConnection:
        vpc_pcx = self.get_vpc_peering_connection(vpc_pcx_id)

        # validate cross-account acceptance
        req_account_id = vpc_pcx.vpc.owner_id
        acp_account_id = vpc_pcx.peer_vpc.owner_id
        if req_account_id != acp_account_id and self.account_id != acp_account_id:  # type: ignore[attr-defined]
            raise OperationNotPermitted5(self.account_id, vpc_pcx_id, "accept")  # type: ignore[attr-defined]

        # validate cross-region acceptance
        pcx_req_region = vpc_pcx.vpc.region
        pcx_acp_region = vpc_pcx.peer_vpc.region
        if pcx_req_region != pcx_acp_region and self.region_name == pcx_req_region:  # type: ignore[attr-defined]
            raise OperationNotPermitted2(self.region_name, vpc_pcx.id, pcx_acp_region)  # type: ignore[attr-defined]

        if vpc_pcx._status.code != "pending-acceptance":
            raise InvalidVPCPeeringConnectionStateTransitionError(vpc_pcx.id)
        vpc_pcx._status.accept()
        return vpc_pcx

    def reject_vpc_peering_connection(self, vpc_pcx_id: str) -> VPCPeeringConnection:
        vpc_pcx = self.get_vpc_peering_connection(vpc_pcx_id)

        # validate cross-account rejection
        req_account_id = vpc_pcx.vpc.owner_id
        acp_account_id = vpc_pcx.peer_vpc.owner_id
        if req_account_id != acp_account_id and self.account_id != acp_account_id:  # type: ignore[attr-defined]
            raise OperationNotPermitted5(self.account_id, vpc_pcx_id, "reject")  # type: ignore[attr-defined]

        # validate cross-region acceptance
        pcx_req_region = vpc_pcx.vpc.region
        pcx_acp_region = vpc_pcx.peer_vpc.region
        if pcx_req_region != pcx_acp_region and self.region_name == pcx_req_region:  # type: ignore[attr-defined]
            raise OperationNotPermitted3(self.region_name, vpc_pcx.id, pcx_acp_region)  # type: ignore[attr-defined]

        if vpc_pcx._status.code != "pending-acceptance":
            raise InvalidVPCPeeringConnectionStateTransitionError(vpc_pcx.id)
        vpc_pcx._status.reject()
        return vpc_pcx

    def modify_vpc_peering_connection_options(
        self,
        vpc_pcx_id: str,
        accepter_options: Optional[Dict[str, Any]] = None,
        requester_options: Optional[Dict[str, Any]] = None,
    ) -> None:
        vpc_pcx = self.get_vpc_peering_connection(vpc_pcx_id)
        if not vpc_pcx:
            raise InvalidVPCPeeringConnectionIdError(vpc_pcx_id)
        # TODO: check if actual vpc has this options enabled
        if accepter_options:
            vpc_pcx.accepter_options.update(accepter_options)
        if requester_options:
            vpc_pcx.requester_options.update(requester_options)
