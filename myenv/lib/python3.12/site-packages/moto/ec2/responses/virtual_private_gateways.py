from ._base_response import EC2BaseResponse


class VirtualPrivateGateways(EC2BaseResponse):
    def attach_vpn_gateway(self) -> str:
        vpn_gateway_id = self._get_param("VpnGatewayId")
        vpc_id = self._get_param("VpcId")
        attachment = self.ec2_backend.attach_vpn_gateway(vpn_gateway_id, vpc_id)
        template = self.response_template(ATTACH_VPN_GATEWAY_RESPONSE)
        return template.render(attachment=attachment)

    def create_vpn_gateway(self) -> str:
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
        template = self.response_template(CREATE_VPN_GATEWAY_RESPONSE)
        return template.render(vpn_gateway=vpn_gateway)

    def delete_vpn_gateway(self) -> str:
        vpn_gateway_id = self._get_param("VpnGatewayId")
        vpn_gateway = self.ec2_backend.delete_vpn_gateway(vpn_gateway_id)
        template = self.response_template(DELETE_VPN_GATEWAY_RESPONSE)
        return template.render(vpn_gateway=vpn_gateway)

    def describe_vpn_gateways(self) -> str:
        filters = self._filters_from_querystring()
        vpn_gw_ids = self._get_multi_param("VpnGatewayId")
        vpn_gateways = self.ec2_backend.describe_vpn_gateways(filters, vpn_gw_ids)
        template = self.response_template(DESCRIBE_VPN_GATEWAYS_RESPONSE)
        return template.render(vpn_gateways=vpn_gateways)

    def detach_vpn_gateway(self) -> str:
        vpn_gateway_id = self._get_param("VpnGatewayId")
        vpc_id = self._get_param("VpcId")
        attachment = self.ec2_backend.detach_vpn_gateway(vpn_gateway_id, vpc_id)
        template = self.response_template(DETACH_VPN_GATEWAY_RESPONSE)
        return template.render(attachment=attachment)


CREATE_VPN_GATEWAY_RESPONSE = """
<CreateVpnGatewayResponse xmlns="http://ec2.amazonaws.com/doc/2013-10-15/">
  <requestId>7a62c49f-347e-4fc4-9331-6e8eEXAMPLE</requestId>
  <vpnGateway>
    <vpnGatewayId>{{ vpn_gateway.id }}</vpnGatewayId>
    {% if vpn_gateway.amazon_side_asn %}
    <amazonSideAsn>{{ vpn_gateway.amazon_side_asn }}</amazonSideAsn>
    {% endif %}
    <state>{{ vpn_gateway.state }}</state>
    <type>{{ vpn_gateway.type }}</type>
    <availabilityZone>{{ vpn_gateway.availability_zone }}</availabilityZone>
    <attachments/>
    <tagSet>
      {% for tag in vpn_gateway.get_tags() %}
        <item>
          <key>{{ tag.key }}</key>
          <value>{{ tag.value }}</value>
        </item>
      {% endfor %}
    </tagSet>
  </vpnGateway>
</CreateVpnGatewayResponse>"""

DESCRIBE_VPN_GATEWAYS_RESPONSE = """
<DescribeVpnGatewaysResponse xmlns="http://ec2.amazonaws.com/doc/2013-10-15/">
  <requestId>7a62c49f-347e-4fc4-9331-6e8eEXAMPLE</requestId>
  <vpnGatewaySet>
    {% for vpn_gateway in vpn_gateways %}
      <item>
        <vpnGatewayId>{{ vpn_gateway.id }}</vpnGatewayId>
        {% if vpn_gateway.amazon_side_asn %}
        <amazonSideAsn>{{ vpn_gateway.amazon_side_asn }}</amazonSideAsn>
        {% endif %}
        <state>{{ vpn_gateway.state }}</state>
        <type>{{ vpn_gateway.type }}</type>
        <availabilityZone>{{ vpn_gateway.availability_zone }}</availabilityZone>
        <attachments>
          {% for attachment in vpn_gateway.attachments.values() %}
            <item>
              <vpcId>{{ attachment.vpc_id }}</vpcId>
              <state>{{ attachment.state }}</state>
            </item>
          {% endfor %}
        </attachments>
        <tagSet>
          {% for tag in vpn_gateway.get_tags() %}
            <item>
              <key>{{ tag.key }}</key>
              <value>{{ tag.value }}</value>
            </item>
          {% endfor %}
        </tagSet>
      </item>
    {% endfor %}
  </vpnGatewaySet>
</DescribeVpnGatewaysResponse>"""

ATTACH_VPN_GATEWAY_RESPONSE = """
<AttachVpnGatewayResponse xmlns="http://ec2.amazonaws.com/doc/2013-10-15/">
   <requestId>7a62c49f-347e-4fc4-9331-6e8eEXAMPLE</requestId>
   <attachment>
      <vpcId>{{ attachment.vpc_id }}</vpcId>
      <state>{{ attachment.state }}</state>
   </attachment>
</AttachVpnGatewayResponse>"""

DELETE_VPN_GATEWAY_RESPONSE = """
<DeleteVpnGatewayResponse xmlns="http://ec2.amazonaws.com/doc/2013-10-15/">
   <requestId>7a62c49f-347e-4fc4-9331-6e8eEXAMPLE</requestId>
   <return>true</return>
</DeleteVpnGatewayResponse>
"""

DETACH_VPN_GATEWAY_RESPONSE = """
<DetachVpnGatewayResponse xmlns="http://ec2.amazonaws.com/doc/2013-10-15/">
   <requestId>7a62c49f-347e-4fc4-9331-6e8eEXAMPLE</requestId>
   <return>true</return>
</DetachVpnGatewayResponse>
"""
