from moto.ec2.utils import add_tag_specification

from ._base_response import EC2BaseResponse


class CarrierGateway(EC2BaseResponse):
    def create_carrier_gateway(self) -> str:
        vpc_id = self._get_param("VpcId")
        tag_param = self._get_multi_param("TagSpecification")
        tags = add_tag_specification(tag_param)

        carrier_gateway = self.ec2_backend.create_carrier_gateway(
            vpc_id=vpc_id, tags=tags
        )
        template = self.response_template(CREATE_CARRIER_GATEWAY_RESPONSE)
        return template.render(carrier_gateway=carrier_gateway)

    def delete_carrier_gateway(self) -> str:
        carrier_gateway_id = self._get_param("CarrierGatewayId")

        carrier_gateway = self.ec2_backend.delete_carrier_gateway(carrier_gateway_id)
        template = self.response_template(DELETE_CARRIER_GATEWAY_RESPONSE)
        return template.render(carrier_gateway=carrier_gateway)

    def describe_carrier_gateways(self) -> str:
        carrier_gateway_ids = self._get_multi_param("CarrierGatewayId")
        filters = self._filters_from_querystring()

        carrier_gateways = self.ec2_backend.describe_carrier_gateways(
            carrier_gateway_ids, filters
        )
        template = self.response_template(DESCRIBE_CARRIER_GATEWAYS_RESPONSE)
        return template.render(carrier_gateways=carrier_gateways)


CREATE_CARRIER_GATEWAY_RESPONSE = """<CreateCarrierGatewayResponse xmlns="http://ec2.amazonaws.com/doc/2016-11-15/">
    <requestId>c617595f-6c29-4a00-a941-example</requestId>
    <carrierGateway>
        <state>{{ carrier_gateway.state }}</state>
        <vpcId>{{ carrier_gateway.vpc_id }}</vpcId>
        <carrierGatewayId>{{ carrier_gateway.id }}</carrierGatewayId>
        <ownerId>{{ carrier_gateway.owner_id }}</ownerId>
        <tagSet>
        {% for tag in carrier_gateway.get_tags() %}
            <item>
                <key>{{ tag.key }}</key>
                <value>{{ tag.value }}</value>
            </item>
        {% endfor %}
        </tagSet>
    </carrierGateway>
</CreateCarrierGatewayResponse>
"""

DELETE_CARRIER_GATEWAY_RESPONSE = """<DeleteCarrierGatewayResponse xmlns="http://ec2.amazonaws.com/doc/2016-11-15/">
    <requestId>c617595f-6c29-4a00-a941-example</requestId>
    <carrierGateway>
        <state>{{ carrier_gateway.state }}</state>
        <vpcId>{{ carrier_gateway.vpc_id }}</vpcId>
        <carrierGatewayId>{{ carrier_gateway.id }}</carrierGatewayId>
        <ownerId>{{ carrier_gateway.owner_id }}</ownerId>
        <tagSet>
        {% for tag in carrier_gateway.get_tags() %}
            <item>
                <key>{{ tag.key }}</key>
                <value>{{ tag.value }}</value>
            </item>
        {% endfor %}
        </tagSet>
    </carrierGateway>
</DeleteCarrierGatewayResponse>
"""

DESCRIBE_CARRIER_GATEWAYS_RESPONSE = """<DescribeCarrierGatewaysResponse xmlns="http://ec2.amazonaws.com/doc/2016-11-15/">
    <requestId>151283df-f7dc-4317-89b4-01c9888b1d45</requestId>
    <carrierGatewaySet>
    {% for carrier_gateway in carrier_gateways %}
        <item>
            <state>{{ carrier_gateway.state }}</state>
            <vpcId>{{ carrier_gateway.vpc_id }}</vpcId>
            <carrierGatewayId>{{ carrier_gateway.id }}</carrierGatewayId>
            <ownerId>{{ carrier_gateway.owner_id }}</ownerId>
            <tagSet>
            {% for tag in carrier_gateway.get_tags() %}
                <item>
                    <key>{{ tag.key }}</key>
                    <value>{{ tag.value }}</value>
                </item>
            {% endfor %}
            </tagSet>
        </item>
    {% endfor %}
    </carrierGatewaySet>
</DescribeCarrierGatewaysResponse>
"""
