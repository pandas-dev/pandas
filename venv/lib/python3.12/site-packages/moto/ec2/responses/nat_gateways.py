from moto.ec2.utils import add_tag_specification

from ._base_response import EC2BaseResponse


class NatGateways(EC2BaseResponse):
    def create_nat_gateway(self) -> str:
        subnet_id = self._get_param("SubnetId")
        allocation_id = self._get_param("AllocationId")
        connectivity_type = self._get_param("ConnectivityType")
        tags = add_tag_specification(self._get_multi_param("TagSpecification"))

        nat_gateway = self.ec2_backend.create_nat_gateway(
            subnet_id=subnet_id,
            allocation_id=allocation_id,
            tags=tags,
            connectivity_type=connectivity_type,
        )
        template = self.response_template(CREATE_NAT_GATEWAY)
        return template.render(nat_gateway=nat_gateway)

    def delete_nat_gateway(self) -> str:
        nat_gateway_id = self._get_param("NatGatewayId")
        nat_gateway = self.ec2_backend.delete_nat_gateway(nat_gateway_id)
        template = self.response_template(DELETE_NAT_GATEWAY_RESPONSE)
        return template.render(nat_gateway=nat_gateway)

    def describe_nat_gateways(self) -> str:
        filters = self._filters_from_querystring()
        nat_gateway_ids = self._get_multi_param("NatGatewayId")
        nat_gateways = self.ec2_backend.describe_nat_gateways(filters, nat_gateway_ids)
        template = self.response_template(DESCRIBE_NAT_GATEWAYS_RESPONSE)
        return template.render(nat_gateways=nat_gateways)


DESCRIBE_NAT_GATEWAYS_RESPONSE = """<DescribeNatGatewaysResponse xmlns="http://ec2.amazonaws.com/doc/2015-10-01/">
    <requestId>bfed02c6-dae9-47c0-86a2-example</requestId>
    <natGatewaySet>
    {% for nat_gateway in nat_gateways %}
         <item>
            <subnetId>{{ nat_gateway.subnet_id }}</subnetId>
            <natGatewayAddressSet>
            {% for address_set in nat_gateway.address_set %}
                <item>
                    {% if address_set.allocationId %}
                    <allocationId>{{ address_set.allocationId }}</allocationId>
                    {% endif %}
                    {% if address_set.privateIp %}
                    <privateIp>{{ address_set.privateIp }}</privateIp>
                    {% endif %}
                    {% if address_set.publicIp %}
                    <publicIp>{{ address_set.publicIp }}</publicIp>
                    {% endif %}
                    {% if address_set.networkInterfaceId %}
                    <networkInterfaceId>{{ address_set.networkInterfaceId }}</networkInterfaceId>
                    {% endif %}
                    {% if address_set.associationId %}
                    <associationId>{{ address_set.associationId }}</associationId>
                    {% endif %}
                </item>
            {% endfor %}
            </natGatewayAddressSet>
            <createTime>{{ nat_gateway.create_time }}</createTime>
            <vpcId>{{ nat_gateway.vpc_id }}</vpcId>
            <natGatewayId>{{ nat_gateway.id }}</natGatewayId>
            <connectivityType>{{ nat_gateway.connectivity_type }}</connectivityType>
            <state>{{ nat_gateway.state }}</state>
            <tagSet>
                {% for tag in nat_gateway.get_tags() %}
                <item>
                    <key>{{ tag.key }}</key>
                    <value>{{ tag.value }}</value>
                </item>
                {% endfor %}
            </tagSet>
        </item>
    {% endfor %}
    </natGatewaySet>
</DescribeNatGatewaysResponse>
"""

CREATE_NAT_GATEWAY = """<CreateNatGatewayResponse xmlns="http://ec2.amazonaws.com/doc/2015-10-01/">
    <requestId>1b74dc5c-bcda-403f-867d-example</requestId>
    <natGateway>
        <subnetId>{{ nat_gateway.subnet_id }}</subnetId>
        <natGatewayAddressSet>
            {% for address_set in nat_gateway.address_set %}
                <item>
                    {% if address_set.allocationId %}
                    <allocationId>{{ address_set.allocationId }}</allocationId>
                    {% endif %}
                    {% if address_set.privateIp %}
                    <privateIp>{{ address_set.privateIp }}</privateIp>
                    {% endif %}
                    {% if address_set.publicIp %}
                    <publicIp>{{ address_set.publicIp }}</publicIp>
                    {% endif %}
                    {% if address_set.networkInterfaceId %}
                    <networkInterfaceId>{{ address_set.networkInterfaceId }}</networkInterfaceId>
                    {% endif %}
                    {% if address_set.associationId %}
                    <associationId>{{ address_set.associationId }}</associationId>
                    {% endif %}
                </item>
            {% endfor %}
        </natGatewayAddressSet>
        <createTime>{{ nat_gateway.create_time }}</createTime>
        <vpcId>{{ nat_gateway.vpc_id }}</vpcId>
        <natGatewayId>{{ nat_gateway.id }}</natGatewayId>
        <connectivityType>{{ nat_gateway.connectivity_type }}</connectivityType>
        <state>{{ nat_gateway.state }}</state>
        <tagSet>
        {% for tag in nat_gateway.get_tags() %}
            <item>
                <key>{{ tag.key }}</key>
                <value>{{ tag.value }}</value>
            </item>
        {% endfor %}
        </tagSet>
    </natGateway>
</CreateNatGatewayResponse>
"""


DELETE_NAT_GATEWAY_RESPONSE = """<DeleteNatGatewayResponse xmlns="http://ec2.amazonaws.com/doc/2015-10-01/">
    <requestId>741fc8ab-6ebe-452b-b92b-example</requestId>
    <natGatewayId>{{ nat_gateway.id }}</natGatewayId>
</DeleteNatGatewayResponse>"""
