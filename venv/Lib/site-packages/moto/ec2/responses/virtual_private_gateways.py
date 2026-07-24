from moto.core.responses import ActionResult, EmptyResult

from ._base_response import EC2BaseResponse


class VirtualPrivateGateways(EC2BaseResponse):
    def attach_vpn_gateway(self) -> ActionResult:
        vpn_gateway_id = self._get_param("VpnGatewayId")
        vpc_id = self._get_param("VpcId")
        attachment = self.ec2_backend.attach_vpn_gateway(vpn_gateway_id, vpc_id)
        return ActionResult({"VpcAttachment": attachment})

    def create_vpn_gateway(self) -> ActionResult:
        gateway_type = self._get_param("Type")
        amazon_side_asn = self._get_param("AmazonSideAsn")
        availability_zone = self._get_param("AvailabilityZone")
        tags = self._parse_tag_specification().get("vpn-gateway", {})
        vpn_gateway = self.ec2_backend.create_vpn_gateway(
            gateway_type=gateway_type,
            amazon_side_asn=amazon_side_asn,
            availability_zone=availability_zone,
            tags=tags,
        )
        return ActionResult({"VpnGateway": vpn_gateway})

    def delete_vpn_gateway(self) -> ActionResult:
        vpn_gateway_id = self._get_param("VpnGatewayId")
        self.ec2_backend.delete_vpn_gateway(vpn_gateway_id)
        return EmptyResult()

    def describe_vpn_gateways(self) -> ActionResult:
        filters = self._filters_from_querystring()
        vpn_gw_ids = self._get_param("VpnGatewayIds", [])
        vpn_gateways = self.ec2_backend.describe_vpn_gateways(filters, vpn_gw_ids)
        return ActionResult({"VpnGateways": vpn_gateways})

    def detach_vpn_gateway(self) -> ActionResult:
        vpn_gateway_id = self._get_param("VpnGatewayId")
        vpc_id = self._get_param("VpcId")
        self.ec2_backend.detach_vpn_gateway(vpn_gateway_id, vpc_id)
        return EmptyResult()
