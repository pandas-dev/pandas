from moto.core.responses import ActionResult
from moto.ec2.utils import add_tag_specification

from ._base_response import EC2BaseResponse


class NatGateways(EC2BaseResponse):
    def create_nat_gateway(self) -> ActionResult:
        subnet_id = self._get_param("SubnetId")
        allocation_id = self._get_param("AllocationId")
        connectivity_type = self._get_param("ConnectivityType")
        tags = add_tag_specification(self._get_param("TagSpecifications", []))

        nat_gateway = self.ec2_backend.create_nat_gateway(
            subnet_id=subnet_id,
            allocation_id=allocation_id,
            tags=tags,
            connectivity_type=connectivity_type,
        )
        return ActionResult({"NatGateway": nat_gateway})

    def delete_nat_gateway(self) -> ActionResult:
        nat_gateway_id = self._get_param("NatGatewayId")
        nat_gateway = self.ec2_backend.delete_nat_gateway(nat_gateway_id)
        return ActionResult({"NatGatewayId": nat_gateway.id})

    def describe_nat_gateways(self) -> ActionResult:
        filters = self._filters_from_querystring()
        nat_gateway_ids = self._get_param("NatGatewayIds", [])
        nat_gateways = self.ec2_backend.describe_nat_gateways(filters, nat_gateway_ids)
        return ActionResult({"NatGateways": nat_gateways})
