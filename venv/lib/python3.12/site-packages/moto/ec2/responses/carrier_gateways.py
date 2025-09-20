from moto.core.responses import ActionResult
from moto.ec2.utils import add_tag_specification

from ._base_response import EC2BaseResponse


class CarrierGateway(EC2BaseResponse):
    def create_carrier_gateway(self) -> ActionResult:
        vpc_id = self._get_param("VpcId")
        tag_param = self._get_multi_param("TagSpecification")
        tags = add_tag_specification(tag_param)

        carrier_gateway = self.ec2_backend.create_carrier_gateway(
            vpc_id=vpc_id, tags=tags
        )
        result = {"CarrierGateway": carrier_gateway}
        return ActionResult(result)

    def delete_carrier_gateway(self) -> ActionResult:
        carrier_gateway_id = self._get_param("CarrierGatewayId")

        carrier_gateway = self.ec2_backend.delete_carrier_gateway(carrier_gateway_id)
        result = {"CarrierGateway": carrier_gateway}
        return ActionResult(result)

    def describe_carrier_gateways(self) -> ActionResult:
        carrier_gateway_ids = self._get_multi_param("CarrierGatewayId")
        filters = self._filters_from_querystring()

        carrier_gateways = self.ec2_backend.describe_carrier_gateways(
            carrier_gateway_ids, filters
        )
        result = {"CarrierGateways": carrier_gateways}
        return ActionResult(result)
