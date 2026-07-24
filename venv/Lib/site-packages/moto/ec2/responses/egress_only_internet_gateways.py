from moto.core.responses import ActionResult
from moto.ec2.utils import add_tag_specification

from ._base_response import EC2BaseResponse


class EgressOnlyInternetGateway(EC2BaseResponse):
    def create_egress_only_internet_gateway(self) -> ActionResult:
        vpc_id = self._get_param("VpcId")
        tag_param = self._get_param("TagSpecifications", [])
        tags = add_tag_specification(tag_param)

        egress_only_igw = self.ec2_backend.create_egress_only_internet_gateway(
            vpc_id=vpc_id, tags=tags
        )
        result = {"EgressOnlyInternetGateway": egress_only_igw}
        return ActionResult(result)

    def describe_egress_only_internet_gateways(self) -> ActionResult:
        egress_only_igw_ids = self._get_param("EgressOnlyInternetGatewayIds", [])
        egress_only_igws = self.ec2_backend.describe_egress_only_internet_gateways(
            egress_only_igw_ids
        )
        result = {"EgressOnlyInternetGateways": egress_only_igws}
        return ActionResult(result)

    def delete_egress_only_internet_gateway(self) -> ActionResult:
        egress_only_igw_id = self._get_param("EgressOnlyInternetGatewayId")
        self.ec2_backend.delete_egress_only_internet_gateway(
            gateway_id=egress_only_igw_id
        )
        return ActionResult({"ReturnCode": True})
