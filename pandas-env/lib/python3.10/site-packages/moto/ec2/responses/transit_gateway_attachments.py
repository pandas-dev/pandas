from moto.ec2.utils import add_tag_specification

from ._base_response import EC2BaseResponse


class TransitGatewayAttachment(EC2BaseResponse):
    def create_transit_gateway_vpc_attachment(self) -> str:
        options = self._get_multi_param_dict("Options")
        subnet_ids = self._get_multi_param("SubnetIds")
        transit_gateway_id = self._get_param("TransitGatewayId")
        vpc_id = self._get_param("VpcId")

        tags = self._parse_tag_specification().get("transit-gateway-route-table", {})

        transit_gateway_attachment = (
            self.ec2_backend.create_transit_gateway_vpc_attachment(
                transit_gateway_id=transit_gateway_id,
                tags=tags,
                vpc_id=vpc_id,
                subnet_ids=subnet_ids,
                options=options,
            )
        )
        template = self.response_template(CREATE_TRANSIT_GATEWAY_VPC_ATTACHMENT)
        return template.render(transit_gateway_attachment=transit_gateway_attachment)

    def describe_transit_gateway_vpc_attachments(self) -> str:
        transit_gateways_attachment_ids = self._get_multi_param(
            "TransitGatewayAttachmentIds"
        )
        filters = self._filters_from_querystring()
        transit_gateway_vpc_attachments = (
            self.ec2_backend.describe_transit_gateway_vpc_attachments(
                transit_gateways_attachment_ids=transit_gateways_attachment_ids,
                filters=filters,
            )
        )
        template = self.response_template(DESCRIBE_TRANSIT_GATEWAY_VPC_ATTACHMENTS)
        return template.render(
            transit_gateway_vpc_attachments=transit_gateway_vpc_attachments
        )

    def modify_transit_gateway_vpc_attachment(self) -> str:
        add_subnet_ids = self._get_multi_param("AddSubnetIds")
        options = self._get_multi_param_dict("Options")
        remove_subnet_ids = self._get_multi_param("RemoveSubnetIds")
        transit_gateway_attachment_id = self._get_param("TransitGatewayAttachmentId")

        transit_gateway_attachment = (
            self.ec2_backend.modify_transit_gateway_vpc_attachment(
                add_subnet_ids=add_subnet_ids,
                options=options,
                remove_subnet_ids=remove_subnet_ids,
                transit_gateway_attachment_id=transit_gateway_attachment_id,
            )
        )
        template = self.response_template(MODIFY_TRANSIT_GATEWAY_VPC_ATTACHMENTS)
        return template.render(transit_gateway_attachment=transit_gateway_attachment)

    def describe_transit_gateway_attachments(self) -> str:
        transit_gateways_attachment_ids = self._get_multi_param(
            "TransitGatewayAttachmentIds"
        )
        filters = self._filters_from_querystring()
        transit_gateway_attachments = (
            self.ec2_backend.describe_transit_gateway_attachments(
                transit_gateways_attachment_ids=transit_gateways_attachment_ids,
                filters=filters,
            )
        )
        template = self.response_template(DESCRIBE_TRANSIT_GATEWAY_ATTACHMENTS)
        return template.render(transit_gateway_attachments=transit_gateway_attachments)

    def delete_transit_gateway_vpc_attachment(self) -> str:
        transit_gateway_attachment_id = self._get_param("TransitGatewayAttachmentId")
        transit_gateway_attachment = (
            self.ec2_backend.delete_transit_gateway_vpc_attachment(
                transit_gateway_attachment_id=transit_gateway_attachment_id
            )
        )
        template = self.response_template(DELETE_TRANSIT_GATEWAY_VPC_ATTACHMENTS)
        return template.render(transit_gateway_attachment=transit_gateway_attachment)

    def associate_transit_gateway_route_table(self) -> str:
        transit_gateway_attachment_id = self._get_param("TransitGatewayAttachmentId")
        transit_gateway_route_table_id = self._get_param("TransitGatewayRouteTableId")
        transit_gateway_association = (
            self.ec2_backend.associate_transit_gateway_route_table(
                transit_gateway_attachment_id=transit_gateway_attachment_id,
                transit_gateway_route_table_id=transit_gateway_route_table_id,
            )
        )
        template = self.response_template(TRANSIT_GATEWAY_ASSOCIATION)
        return template.render(transit_gateway_association=transit_gateway_association)

    def disassociate_transit_gateway_route_table(self) -> str:
        tgw_attach_id = self._get_param("TransitGatewayAttachmentId")
        tgw_rt_id = self._get_param("TransitGatewayRouteTableId")

        tgw_association = self.ec2_backend.disassociate_transit_gateway_route_table(
            tgw_attach_id, tgw_rt_id
        )
        template = self.response_template(TRANSIT_GATEWAY_DISASSOCIATION)
        return template.render(tgw_association=tgw_association)

    def enable_transit_gateway_route_table_propagation(self) -> str:
        transit_gateway_attachment_id = self._get_param("TransitGatewayAttachmentId")
        transit_gateway_route_table_id = self._get_param("TransitGatewayRouteTableId")
        transit_gateway_propagation = (
            self.ec2_backend.enable_transit_gateway_route_table_propagation(
                transit_gateway_attachment_id=transit_gateway_attachment_id,
                transit_gateway_route_table_id=transit_gateway_route_table_id,
            )
        )
        template = self.response_template(TRANSIT_GATEWAY_PROPAGATION)
        return template.render(transit_gateway_propagation=transit_gateway_propagation)

    def disable_transit_gateway_route_table_propagation(self) -> str:
        transit_gateway_attachment_id = self._get_param("TransitGatewayAttachmentId")
        transit_gateway_route_table_id = self._get_param("TransitGatewayRouteTableId")
        transit_gateway_propagation = (
            self.ec2_backend.disable_transit_gateway_route_table_propagation(
                transit_gateway_attachment_id=transit_gateway_attachment_id,
                transit_gateway_route_table_id=transit_gateway_route_table_id,
            )
        )
        template = self.response_template(TRANSIT_GATEWAY_PROPAGATION)
        return template.render(transit_gateway_propagation=transit_gateway_propagation)

    def create_transit_gateway_peering_attachment(self) -> str:
        peer_account_id = self._get_param("PeerAccountId")
        peer_region = self._get_param("PeerRegion")
        peer_transit_gateway_id = self._get_param("PeerTransitGatewayId")
        transit_gateway_id = self._get_param("TransitGatewayId")
        tags = add_tag_specification(self._get_multi_param("TagSpecification"))
        transit_gateway_peering_attachment = (
            self.ec2_backend.create_transit_gateway_peering_attachment(
                transit_gateway_id,
                peer_transit_gateway_id,
                peer_region,
                peer_account_id,
                tags,
            )
        )
        template = self.response_template(TRANSIT_GATEWAY_PEERING_ATTACHMENT)
        return template.render(
            method_name="CreateTransitGatewayPeeringAttachment",
            transit_gateway_peering_attachment=transit_gateway_peering_attachment,
        )

    def describe_transit_gateway_peering_attachments(self) -> str:
        transit_gateways_attachment_ids = self._get_multi_param(
            "TransitGatewayAttachmentIds"
        )
        filters = self._filters_from_querystring()
        transit_gateway_peering_attachments = (
            self.ec2_backend.describe_transit_gateway_peering_attachments(
                transit_gateways_attachment_ids=transit_gateways_attachment_ids,
                filters=filters,
            )
        )
        template = self.response_template(DESCRIBE_TRANSIT_GATEWAY_PEERING_ATTACHMENTS)
        return template.render(
            transit_gateway_peering_attachments=transit_gateway_peering_attachments
        )

    def accept_transit_gateway_peering_attachment(self) -> str:
        transit_gateway_attachment_id = self._get_param("TransitGatewayAttachmentId")
        transit_gateway_peering_attachment = (
            self.ec2_backend.accept_transit_gateway_peering_attachment(
                transit_gateway_attachment_id=transit_gateway_attachment_id
            )
        )
        template = self.response_template(TRANSIT_GATEWAY_PEERING_ATTACHMENT)
        return template.render(
            method_name="AcceptTransitGatewayPeeringAttachment",
            transit_gateway_peering_attachment=transit_gateway_peering_attachment,
        )

    def delete_transit_gateway_peering_attachment(self) -> str:
        transit_gateway_attachment_id = self._get_param("TransitGatewayAttachmentId")
        transit_gateway_peering_attachment = (
            self.ec2_backend.delete_transit_gateway_peering_attachment(
                transit_gateway_attachment_id=transit_gateway_attachment_id
            )
        )
        template = self.response_template(TRANSIT_GATEWAY_PEERING_ATTACHMENT)
        return template.render(
            method_name="DeleteTransitGatewayPeeringAttachment",
            transit_gateway_peering_attachment=transit_gateway_peering_attachment,
        )

    def reject_transit_gateway_peering_attachment(self) -> str:
        transit_gateway_attachment_id = self._get_param("TransitGatewayAttachmentId")
        transit_gateway_peering_attachment = (
            self.ec2_backend.reject_transit_gateway_peering_attachment(
                transit_gateway_attachment_id=transit_gateway_attachment_id
            )
        )
        template = self.response_template(TRANSIT_GATEWAY_PEERING_ATTACHMENT)
        return template.render(
            method_name="RejectTransitGatewayPeeringAttachment",
            transit_gateway_peering_attachment=transit_gateway_peering_attachment,
        )


CREATE_TRANSIT_GATEWAY_VPC_ATTACHMENT = """<CreateTransitGatewayVpcAttachmentResponse xmlns="http://ec2.amazonaws.com/doc/2016-11-15/">
        <requestId>9b5766ac-2af6-4b92-9a8a-4d74ae46ae79</requestId>
        <transitGatewayVpcAttachment>
            <createTime>{{ transit_gateway_attachment.create_time }}</createTime>
            <options>
                <applianceModeSupport>{{ transit_gateway_attachment.options.ApplianceModeSupport }}</applianceModeSupport>
                <dnsSupport>{{ transit_gateway_attachment.options.DnsSupport }}</dnsSupport>
                <ipv6Support>{{ transit_gateway_attachment.options.Ipv6Support }}</ipv6Support>
            </options>
            <state>{{ transit_gateway_attachment.state }}</state>
            <subnetIds>
            {% for subnet_id in transit_gateway_attachment.subnet_ids %}
                <item>{{ subnet_id }}</item>
            {% endfor %}
            </subnetIds>
            <tagSet>
            {% for tag in transit_gateway_attachment.get_tags() %}
                <item>
                    <key>{{ tag.key }}</key>
                    <value>{{ tag.value }}</value>
                </item>
            {% endfor %}
            </tagSet>
            <transitGatewayAttachmentId>{{ transit_gateway_attachment.id }}</transitGatewayAttachmentId>
            <transitGatewayId>{{ transit_gateway_attachment.transit_gateway_id }}</transitGatewayId>
            <vpcId>{{ transit_gateway_attachment.vpc_id }}</vpcId>
            <vpcOwnerId>{{ transit_gateway_attachment.resource_owner_id }}</vpcOwnerId>
    </transitGatewayVpcAttachment>
</CreateTransitGatewayVpcAttachmentResponse>"""


DESCRIBE_TRANSIT_GATEWAY_ATTACHMENTS = """<DescribeTransitGatewayAttachmentsResponse xmlns="http://ec2.amazonaws.com/doc/2016-11-15/">
    <requestId>92aa7885-74c0-42d1-a846-e59bd07488a7</requestId>
    <transitGatewayAttachments>
        {% for transit_gateway_attachment in transit_gateway_attachments %}
        <item>
            <association>
                <state>associated</state>
                <transitGatewayRouteTableId>tgw-rtb-0b36edb9b88f0d5e3</transitGatewayRouteTableId>
            </association>
            <creationTime>2021-07-18T08:57:21.000Z</creationTime>
            <resourceId>{{ transit_gateway_attachment.resource_id }}</resourceId>
            <resourceOwnerId>{{ transit_gateway_attachment.resource_owner_id }}</resourceOwnerId>
            <resourceType>{{ transit_gateway_attachment.resource_type }}</resourceType>
            <state>{{ transit_gateway_attachment.state }}</state>
            <tagSet>
            {% for tag in transit_gateway_attachment.get_tags() %}
                <item>
                    <key>{{ tag.key }}</key>
                    <value>{{ tag.value }}</value>
                </item>
            {% endfor %}
            </tagSet>
            <transitGatewayAttachmentId>{{ transit_gateway_attachment.id }}</transitGatewayAttachmentId>
            <transitGatewayId>{{ transit_gateway_attachment.transit_gateway_id }}</transitGatewayId>
            <transitGatewayOwnerId>{{ transit_gateway_attachment.resource_owner_id }}</transitGatewayOwnerId>
        </item>
        {% endfor %}
    </transitGatewayAttachments>
</DescribeTransitGatewayAttachmentsResponse>
"""


DESCRIBE_TRANSIT_GATEWAY_VPC_ATTACHMENTS = """<DescribeTransitGatewayVpcAttachmentsResponse xmlns="http://ec2.amazonaws.com/doc/2016-11-15/">
        <requestId>bebc9670-0205-4f28-ad89-049c97e46633</requestId>
        <transitGatewayVpcAttachments>
        {% for transit_gateway_vpc_attachment in transit_gateway_vpc_attachments %}
            <item>
                <creationTime>2021-07-18T08:57:21.000Z</creationTime>
                {% if transit_gateway_vpc_attachment.options %}
                <options>
                    <applianceModeSupport>{{ transit_gateway_vpc_attachment.options.ApplianceModeSupport }}</applianceModeSupport>
                    <dnsSupport>{{ transit_gateway_vpc_attachment.options.DnsSupport }}</dnsSupport>
                    <ipv6Support>{{ transit_gateway_vpc_attachment.options.Ipv6Support }}</ipv6Support>
                </options>
                {% endif %}
                <state>{{ transit_gateway_vpc_attachment.state }}</state>
                <subnetIds>
                {% for id in transit_gateway_vpc_attachment.subnet_ids %}
                    <item>{{ id }}</item>
                {% endfor %}
                </subnetIds>
                <tagSet>
                {% for tag in transit_gateway_vpc_attachment.get_tags() %}
                    <item>
                        <key>{{ tag.key }}</key>
                        <value>{{ tag.value }}</value>
                    </item>
                {% endfor %}
                </tagSet>
                <transitGatewayAttachmentId>{{ transit_gateway_vpc_attachment.id }}</transitGatewayAttachmentId>
                <transitGatewayId>{{ transit_gateway_vpc_attachment.transit_gateway_id }}</transitGatewayId>
                <vpcId>{{ transit_gateway_vpc_attachment.vpc_id }}</vpcId>
                <vpcOwnerId>{{ transit_gateway_vpc_attachment.resource_owner_id }}</vpcOwnerId>
            </item>
        {% endfor %}
    </transitGatewayVpcAttachments>
</DescribeTransitGatewayVpcAttachmentsResponse>
"""


MODIFY_TRANSIT_GATEWAY_VPC_ATTACHMENTS = """<ModifyTransitGatewayVpcAttachmentResponse xmlns="http://ec2.amazonaws.com/doc/2016-11-15/">
        <requestId>9b5766ac-2af6-4b92-9a8a-4d74ae46ae79</requestId>
        <transitGatewayVpcAttachment>
            <createTime>{{ transit_gateway_attachment.create_time }}</createTime>
            <options>
                <applianceModeSupport>{{ transit_gateway_attachment.options.ApplianceModeSupport }}</applianceModeSupport>
                <dnsSupport>{{ transit_gateway_attachment.options.DnsSupport }}</dnsSupport>
                <ipv6Support>{{ transit_gateway_attachment.options.Ipv6Support }}</ipv6Support>
            </options>
            <state>{{ transit_gateway_attachment.state }}</state>
            <subnetIds>
            {% for subnet_id in transit_gateway_attachment.subnet_ids %}
                <item>{{ subnet_id }}</item>
            {% endfor %}
            </subnetIds>
            <tagSet>
            {% for tag in transit_gateway_attachment.get_tags() %}
                <item>
                    <key>{{ tag.key }}</key>
                    <value>{{ tag.value }}</value>
                </item>
            {% endfor %}
            </tagSet>
            <transitGatewayAttachmentId>{{ transit_gateway_attachment.id }}</transitGatewayAttachmentId>
            <transitGatewayId>{{ transit_gateway_attachment.transit_gateway_id }}</transitGatewayId>
            <vpcId>{{ transit_gateway_attachment.vpc_id }}</vpcId>
            <vpcOwnerId>{{ transit_gateway_attachment.resource_owner_id }}</vpcOwnerId>
    </transitGatewayVpcAttachment>
</ModifyTransitGatewayVpcAttachmentResponse>"""


DELETE_TRANSIT_GATEWAY_VPC_ATTACHMENTS = """<DeleteTransitGatewayVpcAttachmentResponse xmlns="http://ec2.amazonaws.com/doc/2016-11-15/">
        <requestId>9b5766ac-2af6-4b92-9a8a-4d74ae46ae79</requestId>
        <transitGatewayVpcAttachment>
            <createTime>{{ transit_gateway_attachment.create_time }}</createTime>
            <options>
                <applianceModeSupport>{{ transit_gateway_attachment.options.ApplianceModeSupport }}</applianceModeSupport>
                <dnsSupport>{{ transit_gateway_attachment.options.DnsSupport }}</dnsSupport>
                <ipv6Support>{{ transit_gateway_attachment.options.Ipv6Support }}</ipv6Support>
            </options>
            <state>{{ transit_gateway_attachment.state }}</state>
            <subnetIds>
            {% for subnet_id in transit_gateway_attachment.subnet_ids %}
                <item>{{ subnet_id }}</item>
            {% endfor %}
            </subnetIds>
            <tagSet>
            {% for tag in transit_gateway_attachment.get_tags() %}
                <item>
                    <key>{{ tag.key }}</key>
                    <value>{{ tag.value }}</value>
                </item>
            {% endfor %}
            </tagSet>
            <transitGatewayAttachmentId>{{ transit_gateway_attachment.id }}</transitGatewayAttachmentId>
            <transitGatewayId>{{ transit_gateway_attachment.transit_gateway_id }}</transitGatewayId>
            <vpcId>{{ transit_gateway_attachment.vpc_id }}</vpcId>
            <vpcOwnerId>{{ transit_gateway_attachment.resource_owner_id }}</vpcOwnerId>
    </transitGatewayVpcAttachment>
</DeleteTransitGatewayVpcAttachmentResponse>"""


TRANSIT_GATEWAY_ASSOCIATION = """<AssociateTransitGatewayRouteTableResponse xmlns="http://ec2.amazonaws.com/doc/2016-11-15/">
    <requestId>86a597cf-93ec-44a3-9559-4641863642a5</requestId>
    <association>
        <resourceId>{{ transit_gateway_association.resource_id }}</resourceId>
        <resourceType>{{ transit_gateway_association.resource_type }}</resourceType>
        <state>{{ transit_gateway_association.state }}</state>
        <transitGatewayAttachmentId>{{ transit_gateway_association.transit_gateway_attachment_id }}</transitGatewayAttachmentId>
        <transitGatewayRouteTableId>{{ transit_gateway_association.transit_gateway_route_table_id }}</transitGatewayRouteTableId>
    </association>
</AssociateTransitGatewayRouteTableResponse>
"""


TRANSIT_GATEWAY_DISASSOCIATION = """<DisassociateTransitGatewayRouteTableResponse xmlns="http://ec2.amazonaws.com/doc/2016-11-15/">
    <requestId>86a597cf-93ec-44a3-9559-4641863642a5</requestId>
    <association>
        <resourceId>{{ tgw_association.resource_id }}</resourceId>
        <resourceType>{{ tgw_association.resource_type }}</resourceType>
        <state>{{ tgw_association.state }}</state>
        <transitGatewayAttachmentId>{{ tgw_association.transit_gateway_attachment_id }}</transitGatewayAttachmentId>
        <transitGatewayRouteTableId>{{ tgw_association.transit_gateway_route_table_id }}</transitGatewayRouteTableId>
    </association>
</DisassociateTransitGatewayRouteTableResponse>
"""


TRANSIT_GATEWAY_PROPAGATION = """<EnableTransitGatewayRouteTablePropagationResponse xmlns="http://ec2.amazonaws.com/doc/2016-11-15/">
    <requestId>c78427d4-e498-46ae-bc14-32841b16bff4</requestId>
    <propagation>
        <resourceId>{{ transit_gateway_propagation.resource_id }}</resourceId>
        <resourceType>{{ transit_gateway_propagation.resource_type }}</resourceType>
        <state>{{ transit_gateway_propagation.state }}</state>
        <transitGatewayAttachmentId>{{ transit_gateway_propagation.transit_gateway_attachment_id }}</transitGatewayAttachmentId>
        <transitGatewayRouteTableId>{{ transit_gateway_propagation.transit_gateway_route_table_id }}</transitGatewayRouteTableId>
    </propagation>
</EnableTransitGatewayRouteTablePropagationResponse>
"""


TRANSIT_GATEWAY_PEERING_ATTACHMENT = """<{{ method_name }} xmlns="http://ec2.amazonaws.com/doc/2016-11-15/">
        <requestId>9b5766ac-2af6-4b92-9a8a-4d74ae46ae79</requestId>
        <transitGatewayPeeringAttachment>
            <createTime>{{ transit_gateway_peering_attachment.create_time }}</createTime>
            <state>{{ transit_gateway_peering_attachment.state }}</state>
            <accepterTgwInfo>
                <ownerId>{{ transit_gateway_peering_attachment.accepter_tgw_info.ownerId or '' }}</ownerId>
                <region>{{ transit_gateway_peering_attachment.accepter_tgw_info.region or '' }}</region>
                <transitGatewayId>{{ transit_gateway_peering_attachment.accepter_tgw_info.transitGatewayId or '' }}</transitGatewayId>
            </accepterTgwInfo>
            <requesterTgwInfo>
                <ownerId>{{ transit_gateway_peering_attachment.requester_tgw_info.ownerId or '' }}</ownerId>
                <region>{{ transit_gateway_peering_attachment.requester_tgw_info.region or '' }}</region>
                <transitGatewayId>{{ transit_gateway_peering_attachment.requester_tgw_info.transitGatewayId or '' }}</transitGatewayId>
            </requesterTgwInfo>
            <status>{{ transit_gateway_peering_attachment.status.code }}</status>
            <tagSet>
            {% for tag in transit_gateway_peering_attachment.get_tags() %}
                <item>
                    <key>{{ tag.key }}</key>
                    <value>{{ tag.value }}</value>
                </item>
            {% endfor %}
            </tagSet>
            <transitGatewayAttachmentId>{{ transit_gateway_peering_attachment.id }}</transitGatewayAttachmentId>
    </transitGatewayPeeringAttachment>
</{{ method_name }}>"""


DESCRIBE_TRANSIT_GATEWAY_PEERING_ATTACHMENTS = """<DescribeTransitGatewayPeeringAttachments xmlns="http://ec2.amazonaws.com/doc/2016-11-15/">
        <requestId>bebc9670-0205-4f28-ad89-049c97e46633</requestId>
        <transitGatewayPeeringAttachments>
        {% for transit_gateway_peering_attachment in transit_gateway_peering_attachments %}
            <item>
                <createTime>{{ transit_gateway_peering_attachment.create_time }}</createTime>
                <state>{{ transit_gateway_peering_attachment.state }}</state>
                {% if transit_gateway_peering_attachment.accepter_tgw_info %}
                <accepterTgwInfo>
                    <ownerId>{{ transit_gateway_peering_attachment.accepter_tgw_info.ownerId or '' }}</ownerId>
                    <region>{{ transit_gateway_peering_attachment.accepter_tgw_info.region or '' }}</region>
                    <transitGatewayId>{{ transit_gateway_peering_attachment.accepter_tgw_info.transitGatewayId or '' }}</transitGatewayId>
                </accepterTgwInfo>
                {% endif %}
                {% if transit_gateway_peering_attachment.requester_tgw_info %}
                <requesterTgwInfo>
                    <ownerId>{{ transit_gateway_peering_attachment.requester_tgw_info.ownerId or '' }}</ownerId>
                    <region>{{ transit_gateway_peering_attachment.requester_tgw_info.region or '' }}</region>
                    <transitGatewayId>{{ transit_gateway_peering_attachment.requester_tgw_info.transitGatewayId or '' }}</transitGatewayId>
                </requesterTgwInfo>
                {% endif %}
                {% if transit_gateway_peering_attachment.status %}
                <status>{{ transit_gateway_peering_attachment.status.code }}</status>
                {% endif %}
                <tagSet>
                {% for tag in transit_gateway_peering_attachment.get_tags() %}
                    <item>
                        <key>{{ tag.key }}</key>
                        <value>{{ tag.value }}</value>
                    </item>
                {% endfor %}
                </tagSet>
                <transitGatewayAttachmentId>{{ transit_gateway_peering_attachment.id }}</transitGatewayAttachmentId>
            </item>
        {% endfor %}
    </transitGatewayPeeringAttachments>
</DescribeTransitGatewayPeeringAttachments>
"""
