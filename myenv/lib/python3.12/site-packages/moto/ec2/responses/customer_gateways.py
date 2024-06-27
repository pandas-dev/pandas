from ._base_response import EC2BaseResponse


class CustomerGateways(EC2BaseResponse):
    def create_customer_gateway(self) -> str:
        gateway_type = self._get_param("Type")
        ip_address = self._get_param("IpAddress") or self._get_param("PublicIp")
        bgp_asn = self._get_param("BgpAsn")
        tags = self._parse_tag_specification().get("customer-gateway", {})
        customer_gateway = self.ec2_backend.create_customer_gateway(
            gateway_type, ip_address=ip_address, bgp_asn=bgp_asn, tags=tags
        )
        template = self.response_template(CREATE_CUSTOMER_GATEWAY_RESPONSE)
        return template.render(customer_gateway=customer_gateway)

    def delete_customer_gateway(self) -> str:
        customer_gateway_id = self._get_param("CustomerGatewayId")
        self.ec2_backend.delete_customer_gateway(customer_gateway_id)
        template = self.response_template(DELETE_CUSTOMER_GATEWAY_RESPONSE)
        return template.render(delete_status="true")

    def describe_customer_gateways(self) -> str:
        self.error_on_dryrun()
        filters = self._filters_from_querystring()
        customer_gateway_ids = self._get_multi_param("CustomerGatewayId")
        customer_gateways = self.ec2_backend.describe_customer_gateways(
            filters, customer_gateway_ids
        )
        template = self.response_template(DESCRIBE_CUSTOMER_GATEWAYS_RESPONSE)
        return template.render(customer_gateways=customer_gateways)


CREATE_CUSTOMER_GATEWAY_RESPONSE = """
<CreateCustomerGatewayResponse xmlns="http://ec2.amazonaws.com/doc/2013-10-15/">
   <requestId>7a62c49f-347e-4fc4-9331-6e8eEXAMPLE</requestId>
   <customerGateway>
      <customerGatewayId>{{ customer_gateway.id }}</customerGatewayId>
      <state>{{ customer_gateway.state }}</state>
      <type>{{ customer_gateway.type }}</type>
      <ipAddress>{{ customer_gateway.ip_address }}</ipAddress>
      <bgpAsn>{{ customer_gateway.bgp_asn }}</bgpAsn>
      <tagSet>
        {% for tag in customer_gateway.get_tags() %}
          <item>
            <key>{{ tag.key }}</key>
              <value>{{ tag.value }}</value>
          </item>
        {% endfor %}
      </tagSet>
   </customerGateway>
</CreateCustomerGatewayResponse>"""

DELETE_CUSTOMER_GATEWAY_RESPONSE = """
<DeleteCustomerGatewayResponse xmlns="http://ec2.amazonaws.com/doc/2013-10-15/">
   <requestId>7a62c49f-347e-4fc4-9331-6e8eEXAMPLE</requestId>
   <return>{{ delete_status }}</return>
</DeleteCustomerGatewayResponse>"""

DESCRIBE_CUSTOMER_GATEWAYS_RESPONSE = """
<DescribeCustomerGatewaysResponse xmlns="http://ec2.amazonaws.com/doc/2014-10- 01/">
  <requestId>7a62c49f-347e-4fc4-9331-6e8eEXAMPLE</requestId>
  <customerGatewaySet>
  {% for customer_gateway in customer_gateways %}
    <item>
       <customerGatewayId>{{ customer_gateway.id }}</customerGatewayId>
       <state>{{ customer_gateway.state }}</state>
       <type>{{ customer_gateway.type }}</type>
       <ipAddress>{{ customer_gateway.ip_address }}</ipAddress>
       <bgpAsn>{{ customer_gateway.bgp_asn }}</bgpAsn>
       <tagSet>
        {% for tag in customer_gateway.get_tags() %}
          <item>
            <key>{{ tag.key }}</key>
            <value>{{ tag.value }}</value>
          </item>
        {% endfor %}
       </tagSet>
    </item>
  {% endfor %}
  </customerGatewaySet>
</DescribeCustomerGatewaysResponse>"""
