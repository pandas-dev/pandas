from typing import List

from moto.core.utils import camelcase_to_underscores
from moto.ec2.utils import add_tag_specification

from ..models.managed_prefixes import ManagedPrefixList
from ._base_response import EC2BaseResponse


class VPCs(EC2BaseResponse):
    def _get_doc_date(self) -> str:
        return (
            "2013-10-15"
            if "Boto/" in self.headers.get("user-agent", "")
            else "2016-11-15"
        )

    def create_default_vpc(self) -> str:
        vpc = self.ec2_backend.create_default_vpc()
        doc_date = self._get_doc_date()
        template = self.response_template(CREATE_VPC_RESPONSE)
        return template.render(vpc=vpc, doc_date=doc_date)

    def create_vpc(self) -> str:
        cidr_block = self._get_param("CidrBlock")
        tags = self._get_multi_param("TagSpecification")
        instance_tenancy = self._get_param("InstanceTenancy", if_none="default")
        amazon_provided_ipv6_cidr_block = self._get_param(
            "AmazonProvidedIpv6CidrBlock"
        ) in ["true", "True"]
        ipv6_cidr_block_network_border_group = self._get_param(
            "Ipv6CidrBlockNetworkBorderGroup"
        )
        # if network group is not specified, use the region of the VPC
        if not ipv6_cidr_block_network_border_group:
            ipv6_cidr_block_network_border_group = self.region
        if tags:
            tags = tags[0].get("Tag")
        vpc = self.ec2_backend.create_vpc(
            cidr_block,
            instance_tenancy,
            amazon_provided_ipv6_cidr_block=amazon_provided_ipv6_cidr_block,
            ipv6_cidr_block_network_border_group=ipv6_cidr_block_network_border_group,
            tags=tags,
        )
        doc_date = self._get_doc_date()
        template = self.response_template(CREATE_VPC_RESPONSE)
        return template.render(vpc=vpc, doc_date=doc_date)

    def delete_vpc(self) -> str:
        vpc_id = self._get_param("VpcId")
        vpc = self.ec2_backend.delete_vpc(vpc_id)
        template = self.response_template(DELETE_VPC_RESPONSE)
        return template.render(vpc=vpc)

    def describe_vpcs(self) -> str:
        self.error_on_dryrun()
        vpc_ids = self._get_multi_param("VpcId")
        filters = self._filters_from_querystring()
        vpcs = self.ec2_backend.describe_vpcs(vpc_ids=vpc_ids, filters=filters)
        doc_date = (
            "2013-10-15"
            if "Boto/" in self.headers.get("user-agent", "")
            else "2016-11-15"
        )
        template = self.response_template(DESCRIBE_VPCS_RESPONSE)
        return template.render(vpcs=vpcs, doc_date=doc_date, region=self.region)

    def modify_vpc_tenancy(self) -> str:
        vpc_id = self._get_param("VpcId")
        tenancy = self._get_param("InstanceTenancy")
        self.ec2_backend.modify_vpc_tenancy(vpc_id, tenancy)
        return self.response_template(MODIFY_VPC_TENANCY_RESPONSE).render()

    def describe_vpc_attribute(self) -> str:
        vpc_id = self._get_param("VpcId")
        attribute = self._get_param("Attribute")
        attr_name = camelcase_to_underscores(attribute)
        value = self.ec2_backend.describe_vpc_attribute(vpc_id, attr_name)
        template = self.response_template(DESCRIBE_VPC_ATTRIBUTE_RESPONSE)
        return template.render(vpc_id=vpc_id, attribute=attribute, value=value)

    def describe_vpc_classic_link_dns_support(self) -> str:
        vpc_ids = self._get_multi_param("VpcIds")
        filters = self._filters_from_querystring()
        vpcs = self.ec2_backend.describe_vpcs(vpc_ids=vpc_ids, filters=filters)
        doc_date = self._get_doc_date()
        template = self.response_template(
            DESCRIBE_VPC_CLASSIC_LINK_DNS_SUPPORT_RESPONSE
        )
        return template.render(vpcs=vpcs, doc_date=doc_date)

    def enable_vpc_classic_link_dns_support(self) -> str:
        vpc_id = self._get_param("VpcId")
        classic_link_dns_supported = (
            self.ec2_backend.enable_vpc_classic_link_dns_support(vpc_id=vpc_id)
        )
        doc_date = self._get_doc_date()
        template = self.response_template(ENABLE_VPC_CLASSIC_LINK_DNS_SUPPORT_RESPONSE)
        return template.render(
            classic_link_dns_supported=classic_link_dns_supported, doc_date=doc_date
        )

    def disable_vpc_classic_link_dns_support(self) -> str:
        vpc_id = self._get_param("VpcId")
        classic_link_dns_supported = (
            self.ec2_backend.disable_vpc_classic_link_dns_support(vpc_id=vpc_id)
        )
        doc_date = self._get_doc_date()
        template = self.response_template(DISABLE_VPC_CLASSIC_LINK_DNS_SUPPORT_RESPONSE)
        return template.render(
            classic_link_dns_supported=classic_link_dns_supported, doc_date=doc_date
        )

    def describe_vpc_classic_link(self) -> str:
        vpc_ids = self._get_multi_param("VpcId")
        filters = self._filters_from_querystring()
        vpcs = self.ec2_backend.describe_vpcs(vpc_ids=vpc_ids, filters=filters)
        doc_date = self._get_doc_date()
        template = self.response_template(DESCRIBE_VPC_CLASSIC_LINK_RESPONSE)
        return template.render(vpcs=vpcs, doc_date=doc_date)

    def enable_vpc_classic_link(self) -> str:
        vpc_id = self._get_param("VpcId")
        classic_link_enabled = self.ec2_backend.enable_vpc_classic_link(vpc_id=vpc_id)
        doc_date = self._get_doc_date()
        template = self.response_template(ENABLE_VPC_CLASSIC_LINK_RESPONSE)
        return template.render(
            classic_link_enabled=classic_link_enabled, doc_date=doc_date
        )

    def disable_vpc_classic_link(self) -> str:
        vpc_id = self._get_param("VpcId")
        classic_link_enabled = self.ec2_backend.disable_vpc_classic_link(vpc_id=vpc_id)
        doc_date = self._get_doc_date()
        template = self.response_template(DISABLE_VPC_CLASSIC_LINK_RESPONSE)
        return template.render(
            classic_link_enabled=classic_link_enabled, doc_date=doc_date
        )

    def modify_vpc_attribute(self) -> str:
        vpc_id = self._get_param("VpcId")
        for attribute in (
            "EnableDnsSupport",
            "EnableDnsHostnames",
            "EnableNetworkAddressUsageMetrics",
        ):
            if self.querystring.get(f"{attribute}.Value"):
                attr_name = camelcase_to_underscores(attribute)
                attr_value = self.querystring[f"{attribute}.Value"][0]
                self.ec2_backend.modify_vpc_attribute(vpc_id, attr_name, attr_value)
                return MODIFY_VPC_ATTRIBUTE_RESPONSE
        return ""

    def associate_vpc_cidr_block(self) -> str:
        vpc_id = self._get_param("VpcId")
        amazon_provided_ipv6_cidr_blocks = self._get_param(
            "AmazonProvidedIpv6CidrBlock"
        )
        # todo test on AWS if can create an association for IPV4 and IPV6 in the same call?
        cidr_block = (
            self._get_param("CidrBlock")
            if not amazon_provided_ipv6_cidr_blocks
            else None
        )
        value = self.ec2_backend.associate_vpc_cidr_block(
            vpc_id, cidr_block, amazon_provided_ipv6_cidr_blocks
        )
        if not amazon_provided_ipv6_cidr_blocks:
            render_template = ASSOCIATE_VPC_CIDR_BLOCK_RESPONSE
        else:
            render_template = IPV6_ASSOCIATE_VPC_CIDR_BLOCK_RESPONSE
        template = self.response_template(render_template)
        return template.render(
            vpc_id=vpc_id,
            value=value,
            cidr_block=value["cidr_block"],
            association_id=value["association_id"],
            cidr_block_state="associating",
        )

    def disassociate_vpc_cidr_block(self) -> str:
        association_id = self._get_param("AssociationId")
        value = self.ec2_backend.disassociate_vpc_cidr_block(association_id)
        if "::" in value.get("cidr_block", ""):
            render_template = IPV6_DISASSOCIATE_VPC_CIDR_BLOCK_RESPONSE
        else:
            render_template = DISASSOCIATE_VPC_CIDR_BLOCK_RESPONSE
        template = self.response_template(render_template)
        return template.render(
            vpc_id=value["vpc_id"],
            cidr_block=value["cidr_block"],
            association_id=value["association_id"],
            cidr_block_state="disassociating",
        )

    def create_vpc_endpoint(self) -> str:
        vpc_id = self._get_param("VpcId")
        service_name = self._get_param("ServiceName")
        route_table_ids = self._get_multi_param("RouteTableId")
        subnet_ids = self._get_multi_param("SubnetId")
        endpoint_type = self._get_param("VpcEndpointType")
        policy_document = self._get_param("PolicyDocument")
        client_token = self._get_param("ClientToken")
        private_dns_enabled = self._get_bool_param("PrivateDnsEnabled", if_none=True)
        security_group_ids = self._get_multi_param("SecurityGroupId")

        tags = add_tag_specification(self._get_multi_param("TagSpecification"))
        vpc_end_point = self.ec2_backend.create_vpc_endpoint(
            vpc_id=vpc_id,
            service_name=service_name,
            endpoint_type=endpoint_type,
            policy_document=policy_document,
            route_table_ids=route_table_ids,
            subnet_ids=subnet_ids,
            client_token=client_token,
            security_group_ids=security_group_ids,
            tags=tags,
            private_dns_enabled=private_dns_enabled,
        )
        template = self.response_template(CREATE_VPC_END_POINT)
        return template.render(vpc_end_point=vpc_end_point)

    def modify_vpc_endpoint(self) -> str:
        vpc_id = self._get_param("VpcEndpointId")
        add_subnets = self._get_multi_param("AddSubnetId")
        add_route_tables = self._get_multi_param("AddRouteTableId")
        remove_route_tables = self._get_multi_param("RemoveRouteTableId")
        policy_doc = self._get_param("PolicyDocument")
        self.ec2_backend.modify_vpc_endpoint(
            vpc_id=vpc_id,
            policy_doc=policy_doc,
            add_subnets=add_subnets,
            add_route_tables=add_route_tables,
            remove_route_tables=remove_route_tables,
        )
        template = self.response_template(MODIFY_VPC_END_POINT)
        return template.render()

    def describe_vpc_endpoint_services(self) -> str:
        vpc_end_point_services = self.ec2_backend.describe_vpc_endpoint_services(
            service_names=self._get_multi_param("ServiceName"),
            filters=self._get_multi_param("Filter"),
            max_results=self._get_int_param("MaxResults"),
            next_token=self._get_param("NextToken"),
            region=self.region,
        )
        template = self.response_template(DESCRIBE_VPC_ENDPOINT_SERVICES_RESPONSE)
        return template.render(vpc_end_points=vpc_end_point_services)

    def describe_vpc_endpoints(self) -> str:
        vpc_end_points_ids = self._get_multi_param("VpcEndpointId")
        filters = self._filters_from_querystring()
        vpc_end_points = self.ec2_backend.describe_vpc_endpoints(
            vpc_end_point_ids=vpc_end_points_ids, filters=filters
        )
        template = self.response_template(DESCRIBE_VPC_ENDPOINT_RESPONSE)
        return template.render(
            vpc_end_points=vpc_end_points, account_id=self.current_account
        )

    def delete_vpc_endpoints(self) -> str:
        vpc_end_points_ids = self._get_multi_param("VpcEndpointId")
        self.ec2_backend.delete_vpc_endpoints(vpce_ids=vpc_end_points_ids)
        return self.response_template(DELETE_VPC_ENDPOINT_RESPONSE).render()

    def create_managed_prefix_list(self) -> str:
        address_family = self._get_param("AddressFamily")
        max_entries = self._get_param("MaxEntries")
        prefix_list_name = self._get_param("PrefixListName")
        entry = self._get_multi_param("Entry")

        tags = self._parse_tag_specification().get("prefix-list", {})

        managed_prefix_list = self.ec2_backend.create_managed_prefix_list(
            address_family=address_family,
            entry=entry,
            max_entries=max_entries,
            prefix_list_name=prefix_list_name,
            tags=tags,
        )
        template = self.response_template(CREATE_MANAGED_PREFIX_LIST)
        return template.render(managed_prefix_list=managed_prefix_list)

    def describe_managed_prefix_lists(self) -> str:
        prefix_list_ids = self._get_multi_param("PrefixListId")
        filters = self._filters_from_querystring()
        managed_prefix_lists = self.ec2_backend.describe_managed_prefix_lists(
            prefix_list_ids=prefix_list_ids, filters=filters
        )
        template = self.response_template(DESCRIBE_MANAGED_PREFIX_LIST)
        return template.render(managed_prefix_lists=managed_prefix_lists)

    def get_managed_prefix_list_entries(self) -> str:
        prefix_list_id = self._get_param("PrefixListId")
        target_version = self._get_param("TargetVersion")
        managed_prefix_list = self.ec2_backend.get_managed_prefix_list_entries(
            prefix_list_id=prefix_list_id
        )
        entries: List[ManagedPrefixList] = []
        if managed_prefix_list:
            entries = (
                list(managed_prefix_list.entries.values())[-1]
                if managed_prefix_list.entries.values()
                else []
            )
            if target_version:
                target_version = int(target_version)
                entries = managed_prefix_list.entries.get(target_version)
        template = self.response_template(GET_MANAGED_PREFIX_LIST_ENTRIES)
        return template.render(entries=entries)

    def delete_managed_prefix_list(self) -> str:
        prefix_list_id = self._get_param("PrefixListId")
        managed_prefix_list = self.ec2_backend.delete_managed_prefix_list(
            prefix_list_id
        )
        template = self.response_template(DELETE_MANAGED_PREFIX_LIST)
        return template.render(managed_prefix_list=managed_prefix_list)

    def describe_prefix_lists(self) -> str:
        prefix_list_ids = self._get_multi_param("PrefixListId")
        filters = self._filters_from_querystring()
        managed_pls = self.ec2_backend.describe_managed_prefix_lists(
            prefix_list_ids=prefix_list_ids, filters=filters
        )
        template = self.response_template(DESCRIBE_PREFIX_LIST)
        return template.render(managed_pls=managed_pls)

    def modify_managed_prefix_list(self) -> str:
        add_entry = self._get_multi_param("AddEntry")
        prefix_list_id = self._get_param("PrefixListId")
        current_version = self._get_param("CurrentVersion")
        prefix_list_name = self._get_param("PrefixListName")
        remove_entry = self._get_multi_param("RemoveEntry")

        current_version = int(current_version) if current_version else None

        managed_prefix_list = self.ec2_backend.modify_managed_prefix_list(
            add_entry=add_entry,
            prefix_list_id=prefix_list_id,
            current_version=current_version,
            prefix_list_name=prefix_list_name,
            remove_entry=remove_entry,
        )
        template = self.response_template(MODIFY_PREFIX_LIST)
        return template.render(managed_prefix_list=managed_prefix_list)


CREATE_VPC_RESPONSE = """
<CreateVpcResponse xmlns="http://ec2.amazonaws.com/doc/{{doc_date}}/">
   <requestId>7a62c49f-347e-4fc4-9331-6e8eEXAMPLE</requestId>
   <vpc>
      <vpcId>{{ vpc.id }}</vpcId>
      <state>pending</state>
      <cidrBlock>{{ vpc.cidr_block }}</cidrBlock>
      {% if doc_date == "2016-11-15" %}
          <cidrBlockAssociationSet>
              {% for assoc in vpc.get_cidr_block_association_set() %}
                <item>
                    <cidrBlock>{{assoc.cidr_block}}</cidrBlock>
                    <associationId>{{ assoc.association_id }}</associationId>
                    <cidrBlockState>
                        <state>{{assoc.cidr_block_state.state}}</state>
                    </cidrBlockState>
                </item>
              {% endfor %}
          </cidrBlockAssociationSet>
          <ipv6CidrBlockAssociationSet>
              {% for assoc in vpc.get_cidr_block_association_set(ipv6=True) %}
                <item>
                    <ipv6CidrBlock>{{assoc.cidr_block}}</ipv6CidrBlock>
                    <associationId>{{ assoc.association_id }}</associationId>
                    <ipv6CidrBlockState>
                        <state>{{assoc.cidr_block_state.state}}</state>
                    </ipv6CidrBlockState>
                </item>
              {% endfor %}
          </ipv6CidrBlockAssociationSet>
        {% endif %}
      <dhcpOptionsId>{% if vpc.dhcp_options %}{{ vpc.dhcp_options.id }}{% else %}default{% endif %}</dhcpOptionsId>
      <instanceTenancy>{{ vpc.instance_tenancy }}</instanceTenancy>
      <ownerId> {{ vpc.owner_id }}</ownerId>
      <tagSet>
        {% for tag in vpc.get_tags() %}
          <item>
            <resourceId>{{ tag.resource_id }}</resourceId>
            <resourceType>{{ tag.resource_type }}</resourceType>
            <key>{{ tag.key }}</key>
            <value>{{ tag.value }}</value>
          </item>
        {% endfor %}
      </tagSet>
   </vpc>
</CreateVpcResponse>"""

DESCRIBE_VPC_CLASSIC_LINK_DNS_SUPPORT_RESPONSE = """
<DescribeVpcClassicLinkDnsSupportResponse xmlns="http://ec2.amazonaws.com/doc/{{doc_date}}/">
  <requestId>7a62c442-3484-4f42-9342-6942EXAMPLE</requestId>
  <vpcs>
    {% for vpc in vpcs %}
      <item>
        <vpcId>{{ vpc.id }}</vpcId>
        <classicLinkDnsSupported>{{ vpc.classic_link_dns_supported }}</classicLinkDnsSupported>
      </item>
    {% endfor %}
  </vpcs>
</DescribeVpcClassicLinkDnsSupportResponse>"""

ENABLE_VPC_CLASSIC_LINK_DNS_SUPPORT_RESPONSE = """
<EnableVpcClassicLinkDnsSupportResponse xmlns="http://ec2.amazonaws.com/doc/{{doc_date}}/">
  <requestId>7a62c442-3484-4f42-9342-6942EXAMPLE</requestId>
  <return>{{ classic_link_dns_supported }}</return>
</EnableVpcClassicLinkDnsSupportResponse>"""

DISABLE_VPC_CLASSIC_LINK_DNS_SUPPORT_RESPONSE = """
<DisableVpcClassicLinkDnsSupportResponse xmlns="http://ec2.amazonaws.com/doc/{{doc_date}}/">
  <requestId>7a62c442-3484-4f42-9342-6942EXAMPLE</requestId>
  <return>{{ classic_link_dns_supported }}</return>
</DisableVpcClassicLinkDnsSupportResponse>"""

DESCRIBE_VPC_CLASSIC_LINK_RESPONSE = """
<DescribeVpcClassicLinkResponse xmlns="http://ec2.amazonaws.com/doc/{{doc_date}}/">
  <requestId>7a62c442-3484-4f42-9342-6942EXAMPLE</requestId>
  <vpcSet>
    {% for vpc in vpcs %}
      <item>
        <vpcId>{{ vpc.id }}</vpcId>
        <classicLinkEnabled>{{ vpc.classic_link_enabled }}</classicLinkEnabled>
      </item>
    {% endfor %}
  </vpcSet>
</DescribeVpcClassicLinkResponse>"""

ENABLE_VPC_CLASSIC_LINK_RESPONSE = """
<EnableVpcClassicLinkResponse xmlns="http://ec2.amazonaws.com/doc/{{doc_date}}/">
  <requestId>7a62c442-3484-4f42-9342-6942EXAMPLE</requestId>
  <return>{{ classic_link_enabled }}</return>
</EnableVpcClassicLinkResponse>"""

DISABLE_VPC_CLASSIC_LINK_RESPONSE = """
<DisableVpcClassicLinkResponse xmlns="http://ec2.amazonaws.com/doc/{{doc_date}}/">
  <requestId>7a62c442-3484-4f42-9342-6942EXAMPLE</requestId>
  <return>{{ classic_link_enabled }}</return>
</DisableVpcClassicLinkResponse>"""

DESCRIBE_VPCS_RESPONSE = """
<DescribeVpcsResponse xmlns="http://ec2.amazonaws.com/doc/{{doc_date}}/">
  <requestId>7a62c442-3484-4f42-9342-6942EXAMPLE</requestId>
  <vpcSet>
    {% for vpc in vpcs %}
      <item>
        <vpcId>{{ vpc.id }}</vpcId>
        <state>{{ vpc.state }}</state>
        <cidrBlock>{{ vpc.cidr_block }}</cidrBlock>
        {% if doc_date == "2016-11-15" %}
            <cidrBlockAssociationSet>
              {% for assoc in vpc.get_cidr_block_association_set() %}
                <item>
                    <cidrBlock>{{assoc.cidr_block}}</cidrBlock>
                    <associationId>{{ assoc.association_id }}</associationId>
                    <cidrBlockState>
                        <state>{{assoc.cidr_block_state.state}}</state>
                    </cidrBlockState>
                </item>
              {% endfor %}
            </cidrBlockAssociationSet>
            <ipv6CidrBlockAssociationSet>
              {% for assoc in vpc.get_cidr_block_association_set(ipv6=True) %}
                <item>
                    <ipv6CidrBlock>{{assoc.cidr_block}}</ipv6CidrBlock>
                    <associationId>{{ assoc.association_id }}</associationId>
                    <ipv6CidrBlockState>
                        <state>{{assoc.cidr_block_state.state}}</state>
                    </ipv6CidrBlockState>
                    <networkBorderGroup>{{ assoc.ipv6_cidr_block_network_border_group }}</networkBorderGroup>
                    <ipv6Pool>{{ assoc.ipv6_pool }}</ipv6Pool>
                </item>
              {% endfor %}
            </ipv6CidrBlockAssociationSet>
        {% endif %}
        <dhcpOptionsId>{% if vpc.dhcp_options %}{{ vpc.dhcp_options.id }}{% else %}default{% endif %}</dhcpOptionsId>
        <instanceTenancy>{{ vpc.instance_tenancy }}</instanceTenancy>
        <isDefault>{{ vpc.is_default }}</isDefault>
        <ownerId> {{ vpc.owner_id }}</ownerId>
        <tagSet>
          {% for tag in vpc.get_tags() %}
            <item>
              <resourceId>{{ tag.resource_id }}</resourceId>
              <resourceType>{{ tag.resource_type }}</resourceType>
              <key>{{ tag.key }}</key>
              <value>{{ tag.value }}</value>
            </item>
          {% endfor %}
        </tagSet>
      </item>
    {% endfor %}
  </vpcSet>
</DescribeVpcsResponse>"""

DELETE_VPC_RESPONSE = """
<DeleteVpcResponse xmlns="http://ec2.amazonaws.com/doc/2016-11-15/">
   <requestId>7a62c49f-347e-4fc4-9331-6e8eEXAMPLE</requestId>
   <return>true</return>
</DeleteVpcResponse>
"""

MODIFY_VPC_TENANCY_RESPONSE = """
<ModifyVpcTenancyResponse xmlns="http://ec2.amazonaws.com/doc/2016-11-15/">
   <requestId>7a62c49f-347e-4fc4-9331-6e8eEXAMPLE</requestId>
   <return>true</return>
</ModifyVpcTenancyResponse>
"""

DESCRIBE_VPC_ATTRIBUTE_RESPONSE = """
<DescribeVpcAttributeResponse xmlns="http://ec2.amazonaws.com/doc/2016-11-15/">
  <requestId>7a62c49f-347e-4fc4-9331-6e8eEXAMPLE</requestId>
  <vpcId>{{ vpc_id }}</vpcId>
  <{{ attribute }}>
    <value>{{ value }}</value>
  </{{ attribute }}>
</DescribeVpcAttributeResponse>"""

MODIFY_VPC_ATTRIBUTE_RESPONSE = """
<ModifyVpcAttributeResponse xmlns="http://ec2.amazonaws.com/doc/2016-11-15/">
  <requestId>7a62c49f-347e-4fc4-9331-6e8eEXAMPLE</requestId>
  <return>true</return>
</ModifyVpcAttributeResponse>"""

ASSOCIATE_VPC_CIDR_BLOCK_RESPONSE = """
<AssociateVpcCidrBlockResponse xmlns="http://ec2.amazonaws.com/doc/2016-11-15/">
    <requestId>7a62c49f-347e-4fc4-9331-6e8eEXAMPLE</requestId>
    <vpcId>{{vpc_id}}</vpcId>
    <cidrBlockAssociation>
        <associationId>{{association_id}}</associationId>
        <cidrBlock>{{cidr_block}}</cidrBlock>
        <cidrBlockState>
            <state>{{cidr_block_state}}</state>
        </cidrBlockState>
    </cidrBlockAssociation>
</AssociateVpcCidrBlockResponse>"""

DISASSOCIATE_VPC_CIDR_BLOCK_RESPONSE = """
<DisassociateVpcCidrBlockResponse xmlns="http://ec2.amazonaws.com/doc/2016-11-15/">
    <requestId>7a62c49f-347e-4fc4-9331-6e8eEXAMPLE</requestId>
    <vpcId>{{vpc_id}}</vpcId>
    <cidrBlockAssociation>
        <associationId>{{association_id}}</associationId>
        <cidrBlock>{{cidr_block}}</cidrBlock>
        <cidrBlockState>
            <state>{{cidr_block_state}}</state>
        </cidrBlockState>
    </cidrBlockAssociation>
</DisassociateVpcCidrBlockResponse>"""

IPV6_ASSOCIATE_VPC_CIDR_BLOCK_RESPONSE = """
<AssociateVpcCidrBlockResponse xmlns="http://ec2.amazonaws.com/doc/2016-11-15/">
    <requestId>33af6c54-1139-4d50-b4f7-15a8example</requestId>
    <vpcId>{{vpc_id}}</vpcId>
    <ipv6CidrBlockAssociation>
        <associationId>{{association_id}}</associationId>
        <ipv6CidrBlock>{{cidr_block}}</ipv6CidrBlock>
        <ipv6CidrBlockState>
            <state>{{cidr_block_state}}</state>
        </ipv6CidrBlockState>
    </ipv6CidrBlockAssociation>
</AssociateVpcCidrBlockResponse>"""

IPV6_DISASSOCIATE_VPC_CIDR_BLOCK_RESPONSE = """
<DisassociateVpcCidrBlockResponse xmlns="http://ec2.amazonaws.com/doc/2016-11-15/">
    <requestId>33af6c54-1139-4d50-b4f7-15a8example</requestId>
    <vpcId>{{vpc_id}}</vpcId>
    <ipv6CidrBlockAssociation>
        <associationId>{{association_id}}</associationId>
        <ipv6CidrBlock>{{cidr_block}}</ipv6CidrBlock>
        <ipv6CidrBlockState>
            <state>{{cidr_block_state}}</state>
        </ipv6CidrBlockState>
    </ipv6CidrBlockAssociation>
</DisassociateVpcCidrBlockResponse>"""

CREATE_VPC_END_POINT = """ <CreateVpcEndpointResponse xmlns="http://monitoring.amazonaws.com/doc/2010-08-01/">
    <vpcEndpoint>
        <policyDocument>{{ vpc_end_point.policy_document }}</policyDocument>
        <state>{{ vpc_end_point.state }}</state>
        <vpcEndpointPolicySupported> false </vpcEndpointPolicySupported>
        <serviceName>{{ vpc_end_point.service_name }}</serviceName>
        <vpcId>{{ vpc_end_point.vpc_id }}</vpcId>
        <vpcEndpointId>{{ vpc_end_point.id }}</vpcEndpointId>
        <routeTableIdSet>
            {% for routeid in vpc_end_point.route_table_ids %}
                <item>{{ routeid }}</item>
            {% endfor %}
        </routeTableIdSet>
        <networkInterfaceIdSet>
            {% for network_interface_id in vpc_end_point.network_interface_ids %}
                <item>{{ network_interface_id }}</item>
            {% endfor %}
        </networkInterfaceIdSet>
        <subnetIdSet>
            {% for subnetId in vpc_end_point.subnet_ids %}
                <item>{{ subnetId }}</item>
            {% endfor %}
        </subnetIdSet>
        <privateDnsEnabled>{{ 'true' if vpc_end_point.private_dns_enabled else 'false' }}</privateDnsEnabled>
        <dnsEntrySet>
        {% if vpc_end_point.dns_entries  %}
            {% for entry in vpc_end_point.dns_entries %}
            <item>
                <hostedZoneId>{{ entry["hosted_zone_id"] }}</hostedZoneId>
                <dnsName>{{ entry["dns_name"] }}</dnsName>
            </item>
            {% endfor %}
        {% endif %}
        </dnsEntrySet>
        <tagSet>
        {% for tag in vpc_end_point.get_tags() %}
            <item>
                <key>{{ tag.key }}</key>
                <value>{{ tag.value }}</value>
            </item>
        {% endfor %}
        </tagSet>
        <creationTimestamp>{{ vpc_end_point.created_at }}</creationTimestamp>
    </vpcEndpoint>
</CreateVpcEndpointResponse>"""

MODIFY_VPC_END_POINT = """<ModifyVpcEndpointResponse xmlns="http://monitoring.amazonaws.com/doc/2010-08-01/">
    <return>true</return>
</ModifyVpcEndpointResponse>"""

DESCRIBE_VPC_ENDPOINT_SERVICES_RESPONSE = """<DescribeVpcEndpointServicesResponse xmlns="http://ec2.amazonaws.com/doc/2016-11-15/">
    <requestId>19a9ff46-7df6-49b8-9726-3df27527089d</requestId>
    <serviceNameSet>
        {% for serviceName in vpc_end_points.serviceNames %}
            <item>{{ serviceName }}</item>
        {% endfor %}
    </serviceNameSet>
    <serviceDetailSet>
        {% for service in vpc_end_points.servicesDetails %}
            <item>
                <acceptanceRequired>{{ 'true' if service.AcceptanceRequired else 'false' }}</acceptanceRequired>
                <availabilityZoneSet>
                    {% for zone in service.AvailabilityZones %}
                        <item>{{ zone }}</item>
                    {% endfor %}
                </availabilityZoneSet>
                <baseEndpointDnsNameSet>
                    {% for endpoint in service.BaseEndpointDnsNames %}
                        <item>{{ endpoint }}</item>
                    {% endfor %}
                </baseEndpointDnsNameSet>
                <managesVpcEndpoints>{{ 'true' if service.ManagesVpcEndpoints else 'false' }}</managesVpcEndpoints>
                <owner>{{ service.Owner }}</owner>
                {% if service.PrivateDnsName is defined %}
                    <privateDnsName>{{ service.PrivateDnsName }}</privateDnsName>
                    <privateDnsNameSet>
                        {% for dns_name in service.PrivateDnsNames %}
                            <item>
                                <privateDnsName>{{ dns_name.PrivateDnsName }}</privateDnsName>
                            </item>
                        {% endfor %}
                    </privateDnsNameSet>
                    <privateDnsNameVerificationState>{{ service.PrivateDnsNameVerificationState }}</privateDnsNameVerificationState>
                {% endif %}
                <serviceId>{{ service.ServiceId }}</serviceId>
                <serviceName>{{ service.ServiceName }}</serviceName>
                <serviceType>
                    {% for service_type in service.ServiceType %}
                        <item>
                            <serviceType>{{ service_type.ServiceType }}</serviceType>
                        </item>
                    {% endfor %}
                </serviceType>
                <tagSet>
                    {% for tag in service.Tags %}
                        {% for key, value in tag.items() %}
                            <item>
                                <key>{{ key }}</key>
                                <value>{{ value }}</value>
                            </item>
                        {% endfor %}
                    {% endfor %}
                </tagSet>
                <vpcEndpointPolicySupported>{{ 'true' if service.VpcEndpointPolicySupported else 'false' }}</vpcEndpointPolicySupported>
            </item>
        {% endfor %}
    </serviceDetailSet>
    {% if vpc_end_points.nextToken|length %}
        <nextToken>{{ vpc_end_points.nextToken }}</nextToken>
    {% endif %}
</DescribeVpcEndpointServicesResponse>"""

DESCRIBE_VPC_ENDPOINT_RESPONSE = """<DescribeVpcEndpointsResponse xmlns="http://ec2.amazonaws.com/doc/2016-11-15/">
    <requestId>19a9ff46-7df6-49b8-9726-3df27527089d</requestId>
    <vpcEndpointSet>
        {% for vpc_end_point in vpc_end_points %}
            <item>
                {% if vpc_end_point.policy_document %}
                    <policyDocument>{{ vpc_end_point.policy_document }}</policyDocument>
                {% endif %}
                <state>{{ vpc_end_point.state }}</state>
                <privateDnsEnabled>{{ 'true' if vpc_end_point.private_dns_enabled else 'false' }}</privateDnsEnabled>
                <serviceName>{{ vpc_end_point.service_name }}</serviceName>
                <vpcId>{{ vpc_end_point.vpc_id }}</vpcId>
                <vpcEndpointId>{{ vpc_end_point.id }}</vpcEndpointId>
                <vpcEndpointType>{{ vpc_end_point.endpoint_type }}</vpcEndpointType>
                {% if vpc_end_point.subnet_ids %}
                    <subnetIdSet>
                        {% for subnet_id in vpc_end_point.subnet_ids %}
                            <item>{{ subnet_id }}</item>
                        {% endfor %}
                    </subnetIdSet>
                {% endif %}
                {% if vpc_end_point.route_table_ids %}
                    <routeTableIdSet>
                        {% for route_table_id in vpc_end_point.route_table_ids %}
                            <item>{{ route_table_id }}</item>
                        {% endfor %}
                    </routeTableIdSet>
                {% endif %}
                {% if vpc_end_point.network_interface_ids %}
                    <networkInterfaceIdSet>
                        {% for network_interface_id in vpc_end_point.network_interface_ids %}
                            <item>{{ network_interface_id }}</item>
                        {% endfor %}
                    </networkInterfaceIdSet>
                {% endif %}
                <dnsEntrySet>
                {% if vpc_end_point.dns_entries  %}
                    {% for entry in vpc_end_point.dns_entries %}
                    <item>
                        <hostedZoneId>{{ entry["hosted_zone_id"] }}</hostedZoneId>
                        <dnsName>{{ entry["dns_name"] }}</dnsName>
                    </item>
                    {% endfor %}
                {% endif %}
                </dnsEntrySet>
                {% if vpc_end_point.security_group_ids %}
                    <groupSet>
                        {% for group_id in vpc_end_point.security_group_ids %}
                            <item>
                                <groupId>{{ group_id }}</groupId>
                                <groupName>TODO</groupName>
                            </item>
                        {% endfor %}
                    </groupSet>
                {% endif %}
                <tagSet>
                {% for tag in vpc_end_point.get_tags() %}
                    <item>
                        <key>{{ tag.key }}</key>
                        <value>{{ tag.value }}</value>
                    </item>
                {% endfor %}
                </tagSet>
                <ownerId>{{ account_id }}</ownerId>
                <creationTimestamp>{{ vpc_end_point.created_at }}</creationTimestamp>
            </item>
        {% endfor %}
    </vpcEndpointSet>
</DescribeVpcEndpointsResponse>"""


DELETE_VPC_ENDPOINT_RESPONSE = """<DeleteVpcEndpointsResponse xmlns="http://ec2.amazonaws.com/doc/2016-11-15/">
    <requestId>19a9ff46-7df6-49b8-9726-3df27527089d</requestId>
    <unsuccessful></unsuccessful>
</DeleteVpcEndpointsResponse>"""


CREATE_MANAGED_PREFIX_LIST = """<CreateManagedPrefixListResponse xmlns="http://monitoring.amazonaws.com/doc/2010-08-01/">
    <prefixList>
        <addressFamily>{{ managed_prefix_list.address_family }}</addressFamily>
        <maxEntries>{{ managed_prefix_list.max_entries }}</maxEntries>
        <ownerId>{{ managed_prefix_list.owner_id }}</ownerId>
        <prefixListArn>{{ managed_prefix_list.prefix_list_arn }}</prefixListArn>
        <prefixListId>{{ managed_prefix_list.id }}</prefixListId>
        <prefixListName>{{ managed_prefix_list.prefix_list_name }}</prefixListName>
        <state>{{ managed_prefix_list.state }}</state>
        <tagSet>
        {% for tag in managed_prefix_list.get_tags() %}
            <item>
                <key>{{ tag.key }}</key>
                <value>{{ tag.value }}</value>
            </item>
        {% endfor %}
        </tagSet>
        <version>{{ managed_prefix_list.version }}</version>
    </prefixList>
</CreateManagedPrefixListResponse>"""


DESCRIBE_MANAGED_PREFIX_LIST = """<DescribeManagedPrefixListsResponse xmlns="http://ec2.amazonaws.com/doc/2016-11-15/">
        <requestId>934214d3-4501-4797-b896-13e8fc7ec256</requestId>
        <prefixListSet>
        {% for managed_prefix_list in managed_prefix_lists %}
            <item>
                <addressFamily>{{ managed_prefix_list.address_family }}</addressFamily>
                {% if managed_prefix_list.max_entries %}
                <maxEntries>{{ managed_prefix_list.max_entries }}</maxEntries>
                {% endif %}
                <ownerId>{{ managed_prefix_list.owner_id }}</ownerId>
                <prefixListArn>{{ managed_prefix_list.prefix_list_arn }}</prefixListArn>
                <prefixListId>{{ managed_prefix_list.id }}</prefixListId>
                <prefixListName>{{ managed_prefix_list.prefix_list_name }}</prefixListName>
                <state>{{ managed_prefix_list.state }}</state>
                <tagSet>
                {% for tag in managed_prefix_list.get_tags() %}
                    <item>
                        <key>{{ tag.key }}</key>
                        <value>{{ tag.value }}</value>
                    </item>
                {% endfor %}
                </tagSet>
                {% if managed_prefix_list.version %}
                <version>{{ managed_prefix_list.version }}</version>
                {% endif %}
            </item>
        {% endfor %}
    </prefixListSet>
</DescribeManagedPrefixListsResponse>
"""


GET_MANAGED_PREFIX_LIST_ENTRIES = """<GetManagedPrefixListEntriesResponse xmlns="http://ec2.amazonaws.com/doc/2016-11-15/">
    <requestId>39a3c79f-846f-4382-a592-example</requestId>
    <entrySet>
    {% for entry in entries %}
        <item>
            <cidr>{{ entry.Cidr or ''}}</cidr>
            <description>{{ entry.Description or ''}}</description>
        </item>
    {% endfor %}
    </entrySet>
</GetManagedPrefixListEntriesResponse>
"""


DELETE_MANAGED_PREFIX_LIST = """<DeleteManagedPrefixListResponse xmlns="http://ec2.amazonaws.com/doc/2016-11-15/">
    <requestId>39a3c79f-846f-4382-a592-example</requestId>
    <prefixList>
        <addressFamily>{{ managed_prefix_list.address_family }}</addressFamily>
        <maxEntries>{{ managed_prefix_list.max_entries or '' }}</maxEntries>
        <ownerId>{{ managed_prefix_list.owner_id }}</ownerId>
        <prefixListArn>{{ managed_prefix_list.prefix_list_arn }}</prefixListArn>
        <prefixListId>{{ managed_prefix_list.id }}</prefixListId>
        <prefixListName>{{ managed_prefix_list.prefix_list_name }}</prefixListName>
        <state>{{ managed_prefix_list.state }}</state>
        <tagSet>
        {% for tag in managed_prefix_list.get_tags() %}
            <item>
                <key>{{ tag.key }}</key>
                <value>{{ tag.value }}</value>
            </item>
        {% endfor %}
        </tagSet>
        <version>{{ managed_prefix_list.version or ''}}</version>
    </prefixList>
</DeleteManagedPrefixListResponse>
"""


DESCRIBE_PREFIX_LIST = """<DescribePrefixListsResponse xmlns="http://ec2.amazonaws.com/doc/2016-11-15/">
    <requestId>8a2ec0e2-6918-4270-ae45-58e61971e97d</requestId>
    <prefixListSet>
    {% for pl in managed_pls %}
    {% if pl.prefix_list_name and pl.prefix_list_name.startswith("com.amazonaws.") %}
        <item>
            <cidrSet>
                {% for entry in pl.entries.1 %}
                <item>{{ entry.Cidr }}</item>
                {% endfor %}
            </cidrSet>
            <prefixListId>{{ pl.id }}</prefixListId>
            <prefixListName>{{ pl.prefix_list_name }}</prefixListName>
        </item>
    {% endif %}
    {% endfor %}
    </prefixListSet>
</DescribePrefixListsResponse>
"""

MODIFY_PREFIX_LIST = """<ModifyManagedPrefixListResponse xmlns="http://ec2.amazonaws.com/doc/2016-11-15/">
    <requestId>602f3752-c348-4b14-81e2-example</requestId>
    <prefixList>
        <addressFamily>{{ managed_prefix_list.address_family }}</addressFamily>
        <maxEntries>{{ managed_prefix_list.max_entries or '' }}</maxEntries>
        <ownerId>{{ managed_prefix_list.owner_id }}</ownerId>
        <prefixListArn>{{ managed_prefix_list.prefix_list_arn }}</prefixListArn>
        <prefixListId>{{ managed_prefix_list.id }}</prefixListId>
        <prefixListName>{{ managed_prefix_list.prefix_list_name }}</prefixListName>
        <state>{{ managed_prefix_list.state }}</state>
        <tagSet>
        {% for tag in managed_prefix_list.get_tags() %}
            <item>
                <key>{{ tag.key }}</key>
                <value>{{ tag.value }}</value>
            </item>
        {% endfor %}
        </tagSet>
        <version>{{ managed_prefix_list.version or ''}}</version>
    </prefixList>
</ModifyManagedPrefixListResponse>
"""
