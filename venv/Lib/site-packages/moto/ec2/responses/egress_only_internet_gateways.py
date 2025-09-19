from moto.ec2.utils import add_tag_specification

from ._base_response import EC2BaseResponse


class EgressOnlyInternetGateway(EC2BaseResponse):
    def create_egress_only_internet_gateway(self) -> str:
        vpc_id = self._get_param("VpcId")
        tag_param = self._get_multi_param("TagSpecification")
        tags = add_tag_specification(tag_param)

        egress_only_igw = self.ec2_backend.create_egress_only_internet_gateway(
            vpc_id=vpc_id, tags=tags
        )
        template = self.response_template(CREATE_EGRESS_ONLY_IGW_RESPONSE)
        return template.render(egress_only_igw=egress_only_igw)

    def describe_egress_only_internet_gateways(self) -> str:
        egress_only_igw_ids = self._get_multi_param("EgressOnlyInternetGatewayId")
        egress_only_igws = self.ec2_backend.describe_egress_only_internet_gateways(
            egress_only_igw_ids
        )
        template = self.response_template(DESCRIBE_EGRESS_ONLY_IGW_RESPONSE)
        return template.render(egress_only_igws=egress_only_igws)

    def delete_egress_only_internet_gateway(self) -> str:
        egress_only_igw_id = self._get_param("EgressOnlyInternetGatewayId")
        self.ec2_backend.delete_egress_only_internet_gateway(
            gateway_id=egress_only_igw_id
        )
        template = self.response_template(DELETE_EGRESS_ONLY_IGW_RESPONSE)
        return template.render()


CREATE_EGRESS_ONLY_IGW_RESPONSE = """<CreateEgressOnlyInternetGatewayResponse xmlns="http://ec2.amazonaws.com/doc/2016-11-15/">
    <requestId>c617595f-6c29-4a00-a941-example</requestId>
    <egressOnlyInternetGateway>
        <attachmentSet>
            <item>
                <state>{{ egress_only_igw.state }}</state>
                <vpcId>{{ egress_only_igw.vpc_id }}</vpcId>
            </item>
        </attachmentSet>
        <egressOnlyInternetGatewayId>{{ egress_only_igw.id }}</egressOnlyInternetGatewayId>
        <tagSet>
        {% for tag in egress_only_igw.get_tags() %}
            <item>
                <key>{{ tag.key }}</key>
                <value>{{ tag.value }}</value>
            </item>
        {% endfor %}
        </tagSet>
    </egressOnlyInternetGateway>
</CreateEgressOnlyInternetGatewayResponse>
"""

DESCRIBE_EGRESS_ONLY_IGW_RESPONSE = """<DescribeEgressOnlyInternetGatewaysResponse xmlns="http://ec2.amazonaws.com/doc/2016-11-15/">
    <requestId>ec441b4c-357f-4483-b4a7-example</requestId>
    <egressOnlyInternetGatewaySet>
        {% for egress_only_igw in egress_only_igws %}
        <item>
            <attachmentSet>
                <item>
                    <state>{{ egress_only_igw.state }}</state>
                    <vpcId>{{ egress_only_igw.vpc_id }}</vpcId>
                </item>
            </attachmentSet>
            <egressOnlyInternetGatewayId>{{ egress_only_igw.id }}</egressOnlyInternetGatewayId>
            <tagSet>
            {% for tag in egress_only_igw.get_tags() %}
                <item>
                    <key>{{ tag.key }}</key>
                    <value>{{ tag.value }}</value>
                </item>
            {% endfor %}
            </tagSet>
        </item>
        {% endfor %}
    </egressOnlyInternetGatewaySet>
</DescribeEgressOnlyInternetGatewaysResponse>"""

DELETE_EGRESS_ONLY_IGW_RESPONSE = """<DeleteEgressOnlyInternetGateway xmlns="http://ec2.amazonaws.com/doc/2016-11-15/">
  <requestId>59dbff89-35bd-4eac-99ed-be587EXAMPLE</requestId>
  <returnCode>true</returnCode>
</DeleteEgressOnlyInternetGateway>"""
