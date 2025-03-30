from ._base_response import EC2BaseResponse


class TransitGateways(EC2BaseResponse):
    def create_transit_gateway(self) -> str:
        description = self._get_param("Description") or None
        options = self._get_multi_param_dict("Options")
        tags = self._get_multi_param("TagSpecification")
        if tags:
            tags = tags[0].get("Tag")

        transit_gateway = self.ec2_backend.create_transit_gateway(
            description=description, options=options, tags=tags
        )

        # creating default route table
        transit_gateway_route_table = (
            self.ec2_backend.create_transit_gateway_route_table(
                transit_gateway_id=transit_gateway.id,
                tags={},
                default_association_route_table=True,
                default_propagation_route_table=True,
            )
        )
        transit_gateway.options["AssociationDefaultRouteTableId"] = (
            transit_gateway_route_table.id
        )
        transit_gateway.options["PropagationDefaultRouteTableId"] = (
            transit_gateway_route_table.id
        )

        template = self.response_template(CREATE_TRANSIT_GATEWAY_RESPONSE)
        return template.render(transit_gateway=transit_gateway)

    def delete_transit_gateway(self) -> str:
        transit_gateway_id = self._get_param("TransitGatewayId")
        transit_gateway = self.ec2_backend.delete_transit_gateway(transit_gateway_id)
        template = self.response_template(DELETE_TRANSIT_GATEWAY_RESPONSE)
        return template.render(transit_gateway=transit_gateway)

    def describe_transit_gateways(self) -> str:
        transit_gateway_ids = self._get_multi_param("TransitGatewayIds")
        filters = self._filters_from_querystring()
        transit_gateways = self.ec2_backend.describe_transit_gateways(
            filters, transit_gateway_ids
        )
        template = self.response_template(DESCRIBE_TRANSIT_GATEWAY_RESPONSE)
        return template.render(transit_gateways=transit_gateways)

    def modify_transit_gateway(self) -> str:
        transit_gateway_id = self._get_param("TransitGatewayId")
        description = self._get_param("Description") or None
        options = self._get_multi_param_dict("Options")
        transit_gateway = self.ec2_backend.modify_transit_gateway(
            transit_gateway_id=transit_gateway_id,
            description=description,
            options=options,
        )
        template = self.response_template(MODIFY_TRANSIT_GATEWAY_RESPONSE)
        return template.render(transit_gateway=transit_gateway)


CREATE_TRANSIT_GATEWAY_RESPONSE = """<CreateTransitGatewayResponse xmlns="http://ec2.amazonaws.com/doc/2016-11-15/">
    <requestId>151283df-f7dc-4317-89b4-01c9888b1d45</requestId>
    <transitGateway>
        <transitGatewayId>{{ transit_gateway.id }}</transitGatewayId>
        <transitGatewayArn>{{ transit_gateway.arn }}</transitGatewayArn>
        <ownerId>{{ transit_gateway.owner_id }}</ownerId>
        <description>{{ transit_gateway.description or '' }}</description>
        <createTime>{{ transit_gateway.create_time }}</createTime>
        <state>{{ transit_gateway.state }}</state>
        {% if transit_gateway.options %}
        <options>
            <amazonSideAsn>{{ transit_gateway.options.AmazonSideAsn }}</amazonSideAsn>
            <autoAcceptSharedAttachments>{{ transit_gateway.options.AutoAcceptSharedAttachments }}</autoAcceptSharedAttachments>
            <defaultRouteTableAssociation>{{ transit_gateway.options.DefaultRouteTableAssociation }}</defaultRouteTableAssociation>
            <defaultRouteTablePropagation>{{ transit_gateway.options.DefaultRouteTablePropagation }}</defaultRouteTablePropagation>
            <dnsSupport>{{ transit_gateway.options.DnsSupport }}</dnsSupport>
            <propagationDefaultRouteTableId>{{ transit_gateway.options.PropagationDefaultRouteTableId }}</propagationDefaultRouteTableId>
            <vpnEcmpSupport>{{ transit_gateway.options.VpnEcmpSupport }}</vpnEcmpSupport>
            <transitGatewayCidrBlocks>{{ transit_gateway.options.TransitGatewayCidrBlocks }}</transitGatewayCidrBlocks>
        </options>
        {% endif %}
        <tagSet>
            {% for tag in transit_gateway.get_tags() %}
                <item>
                    <key>{{ tag.key }}</key>
                    <value>{{ tag.value }}</value>
                </item>
            {% endfor %}
        </tagSet>
    </transitGateway>
</CreateTransitGatewayResponse>
"""

DESCRIBE_TRANSIT_GATEWAY_RESPONSE = """<DescribeTransitGatewaysResponse xmlns="http://ec2.amazonaws.com/doc/2016-11-15/">
    <requestId>151283df-f7dc-4317-89b4-01c9888b1d45</requestId>
    <transitGatewaySet>
    {% for transit_gateway in transit_gateways %}
        <item>
            <creationTime>{{ transit_gateway.create_time }}</creationTime>
            <description>{{ transit_gateway.description or '' }}</description>
            {% if transit_gateway.options %}
            <options>
                <amazonSideAsn>{{ transit_gateway.options.AmazonSideAsn }}</amazonSideAsn>
                <associationDefaultRouteTableId>{{ transit_gateway.options.AssociationDefaultRouteTableId }}</associationDefaultRouteTableId>
                <autoAcceptSharedAttachments>{{ transit_gateway.options.AutoAcceptSharedAttachments }}</autoAcceptSharedAttachments>
                <defaultRouteTableAssociation>{{ transit_gateway.options.DefaultRouteTableAssociation }}</defaultRouteTableAssociation>
                <defaultRouteTablePropagation>{{ transit_gateway.options.DefaultRouteTablePropagation }}</defaultRouteTablePropagation>
                <dnsSupport>{{ transit_gateway.options.DnsSupport }}</dnsSupport>
                <propagationDefaultRouteTableId>{{ transit_gateway.options.PropagationDefaultRouteTableId }}</propagationDefaultRouteTableId>
                <vpnEcmpSupport>{{ transit_gateway.options.VpnEcmpSupport }}</vpnEcmpSupport>
                <transitGatewayCidrBlocks>{{ transit_gateway.options.TransitGatewayCidrBlocks }}</transitGatewayCidrBlocks>
            </options>
            {% endif %}
            <ownerId>{{ transit_gateway.owner_id }}</ownerId>
            <state>{{ transit_gateway.state }}</state>
            <tagSet>
                {% for tag in transit_gateway.get_tags() %}
                    <item>
                        <key>{{ tag.key }}</key>
                        <value>{{ tag.value }}</value>
                    </item>
                {% endfor %}
            </tagSet>
            <transitGatewayArn>{{ transit_gateway.arn }}</transitGatewayArn>
            <transitGatewayId>{{ transit_gateway.id }}</transitGatewayId>
        </item>
    {% endfor %}
    </transitGatewaySet>
</DescribeTransitGatewaysResponse>
"""

DELETE_TRANSIT_GATEWAY_RESPONSE = """<DeleteTransitGatewayResponse xmlns="http://ec2.amazonaws.com/doc/2016-11-15/">
    <requestId>151283df-f7dc-4317-89b4-01c9888b1d45</requestId>
    <transitGatewayId>{{ transit_gateway.id }}</transitGatewayId>
</DeleteTransitGatewayResponse>
"""


MODIFY_TRANSIT_GATEWAY_RESPONSE = """<ModifyTransitGatewaysResponse xmlns="http://ec2.amazonaws.com/doc/2016-11-15/">
    <requestId>151283df-f7dc-4317-89b4-01c9888b1d45</requestId>
    <transitGatewaySet>
    <item>
        <creationTime>{{ transit_gateway.create_time }}</creationTime>
        <description>{{ transit_gateway.description or '' }}</description>
        {% if transit_gateway.options %}
        <options>
            <amazonSideAsn>{{ transit_gateway.options.AmazonSideAsn }}</amazonSideAsn>
            <associationDefaultRouteTableId>{{ transit_gateway.options.AssociationDefaultRouteTableId }}</associationDefaultRouteTableId>
            <autoAcceptSharedAttachments>{{ transit_gateway.options.AutoAcceptSharedAttachments }}</autoAcceptSharedAttachments>
            <defaultRouteTableAssociation>{{ transit_gateway.options.DefaultRouteTableAssociation }}</defaultRouteTableAssociation>
            <defaultRouteTablePropagation>{{ transit_gateway.options.DefaultRouteTablePropagation }}</defaultRouteTablePropagation>
            <dnsSupport>{{ transit_gateway.options.DnsSupport }}</dnsSupport>
            <propagationDefaultRouteTableId>{{ transit_gateway.options.PropagationDefaultRouteTableId }}</propagationDefaultRouteTableId>
            <vpnEcmpSupport>{{ transit_gateway.options.VpnEcmpSupport }}</vpnEcmpSupport>
            <transitGatewayCidrBlocks>{{ transit_gateway.options.TransitGatewayCidrBlocks }}</transitGatewayCidrBlocks>
        </options>
        {% endif %}
        <ownerId>{{ transit_gateway.owner_id }}</ownerId>
        <state>{{ transit_gateway.state }}</state>
        <tagSet>
            {% for tag in transit_gateway.get_tags() %}
                <item>
                    <key>{{ tag.key }}</key>
                    <value>{{ tag.value }}</value>
                </item>
            {% endfor %}
        </tagSet>
        <transitGatewayArn>{{ transit_gateway.arn }}</transitGatewayArn>
        <transitGatewayId>{{ transit_gateway.id }}</transitGatewayId>
    </item>
    </transitGatewaySet>
</ModifyTransitGatewaysResponse>
"""
