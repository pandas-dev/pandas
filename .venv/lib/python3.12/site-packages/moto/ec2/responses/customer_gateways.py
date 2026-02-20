from moto.core.responses import ActionResult

from ._base_response import EC2BaseResponse


class CustomerGateways(EC2BaseResponse):
    def create_customer_gateway(self) -> ActionResult:
        gateway_type = self._get_param("Type")
        ip_address = self._get_param("IpAddress") or self._get_param("PublicIp")
        bgp_asn = self._get_param("BgpAsn")
        tags = self._parse_tag_specification().get("customer-gateway", {})
        customer_gateway = self.ec2_backend.create_customer_gateway(
            gateway_type, ip_address=ip_address, bgp_asn=bgp_asn, tags=tags
        )
        result = {"CustomerGateway": customer_gateway}
        return ActionResult(result)

    def delete_customer_gateway(self) -> ActionResult:
        customer_gateway_id = self._get_param("CustomerGatewayId")
        self.ec2_backend.delete_customer_gateway(customer_gateway_id)
        return ActionResult({})

    def describe_customer_gateways(self) -> ActionResult:
        self.error_on_dryrun()
        filters = self._filters_from_querystring()
        customer_gateway_ids = self._get_param("CustomerGatewayIds", [])
        customer_gateways = self.ec2_backend.describe_customer_gateways(
            filters, customer_gateway_ids
        )
        result = {"CustomerGateways": customer_gateways}
        return ActionResult(result)
