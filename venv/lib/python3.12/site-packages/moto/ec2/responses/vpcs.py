from moto.core.responses import ActionResult, EmptyResult
from moto.core.utils import camelcase_to_underscores
from moto.ec2.utils import add_tag_specification

from ..models.managed_prefixes import ManagedPrefixList
from ._base_response import EC2BaseResponse


class VPCs(EC2BaseResponse):
    def create_default_vpc(self) -> ActionResult:
        vpc = self.ec2_backend.create_default_vpc()
        result = {"Vpc": vpc}
        return ActionResult(result)

    def create_vpc(self) -> ActionResult:
        cidr_block = self._get_param("CidrBlock")
        tags = self._get_param("TagSpecifications", [])
        instance_tenancy = self._get_param("InstanceTenancy", if_none="default")
        amazon_provided_ipv6_cidr_block = self._get_param(
            "AmazonProvidedIpv6CidrBlock", False
        )
        ipv6_cidr_block_network_border_group = self._get_param(
            "Ipv6CidrBlockNetworkBorderGroup"
        )
        # if network group is not specified, use the region of the VPC
        if not ipv6_cidr_block_network_border_group:
            ipv6_cidr_block_network_border_group = self.region
        if tags:
            tags = tags[0].get("Tags")
        vpc = self.ec2_backend.create_vpc(
            cidr_block,
            instance_tenancy,
            amazon_provided_ipv6_cidr_block=amazon_provided_ipv6_cidr_block,
            ipv6_cidr_block_network_border_group=ipv6_cidr_block_network_border_group,
            tags=tags,
        )
        result = {"Vpc": vpc}
        return ActionResult(result)

    def delete_vpc(self) -> ActionResult:
        vpc_id = self._get_param("VpcId")
        self.ec2_backend.delete_vpc(vpc_id)
        return EmptyResult()

    def describe_vpcs(self) -> ActionResult:
        self.error_on_dryrun()
        vpc_ids = self._get_param("VpcIds", [])
        filters = self._filters_from_querystring()
        vpcs = self.ec2_backend.describe_vpcs(vpc_ids=vpc_ids, filters=filters)
        result = {"Vpcs": vpcs}
        return ActionResult(result)

    def modify_vpc_tenancy(self) -> ActionResult:
        vpc_id = self._get_param("VpcId")
        tenancy = self._get_param("InstanceTenancy")
        self.ec2_backend.modify_vpc_tenancy(vpc_id, tenancy)
        return EmptyResult()

    def describe_vpc_attribute(self) -> ActionResult:
        vpc_id = self._get_param("VpcId")
        attribute = self._get_param("Attribute")
        attr_name = camelcase_to_underscores(attribute)
        value = self.ec2_backend.describe_vpc_attribute(vpc_id, attr_name)
        result = {"VpcId": vpc_id, attribute: {"Value": value}}
        return ActionResult(result)

    def describe_vpc_classic_link_dns_support(self) -> ActionResult:
        vpc_ids = self._get_param("VpcIds")
        filters = self._filters_from_querystring()
        vpcs = self.ec2_backend.describe_vpcs(vpc_ids=vpc_ids, filters=filters)
        result = {
            "Vpcs": [
                {
                    "VpcId": vpc.id,
                    "ClassicLinkDnsSupported": vpc.classic_link_dns_supported,
                }
                for vpc in vpcs
            ]
        }
        return ActionResult(result)

    def enable_vpc_classic_link_dns_support(self) -> ActionResult:
        vpc_id = self._get_param("VpcId")
        ret = self.ec2_backend.enable_vpc_classic_link_dns_support(vpc_id=vpc_id)
        return ActionResult({"Return": ret})

    def disable_vpc_classic_link_dns_support(self) -> ActionResult:
        vpc_id = self._get_param("VpcId")
        ret = self.ec2_backend.disable_vpc_classic_link_dns_support(vpc_id=vpc_id)
        return ActionResult({"Return": ret})

    def describe_vpc_classic_link(self) -> ActionResult:
        vpc_ids = self._get_param("VpcIds")
        filters = self._filters_from_querystring()
        vpcs = self.ec2_backend.describe_vpcs(vpc_ids=vpc_ids, filters=filters)
        result = {
            "Vpcs": [
                {
                    "VpcId": vpc.id,
                    "ClassicLinkEnabled": vpc.classic_link_enabled,
                    "Tags": vpc.tag_list,
                }
                for vpc in vpcs
            ]
        }
        return ActionResult(result)

    def enable_vpc_classic_link(self) -> ActionResult:
        vpc_id = self._get_param("VpcId")
        ret = self.ec2_backend.enable_vpc_classic_link(vpc_id=vpc_id)
        return ActionResult({"Return": ret})

    def disable_vpc_classic_link(self) -> ActionResult:
        vpc_id = self._get_param("VpcId")
        ret = self.ec2_backend.disable_vpc_classic_link(vpc_id=vpc_id)
        return ActionResult({"Return": ret})

    def modify_vpc_attribute(self) -> ActionResult:
        vpc_id = self._get_param("VpcId")
        for attribute in (
            "EnableDnsSupport",
            "EnableDnsHostnames",
            "EnableNetworkAddressUsageMetrics",
        ):
            if self._get_param(f"{attribute}.Value") is not None:
                attr_name = camelcase_to_underscores(attribute)
                attr_value = self._get_param(f"{attribute}.Value")
                self.ec2_backend.modify_vpc_attribute(vpc_id, attr_name, attr_value)
        return EmptyResult()

    def associate_vpc_cidr_block(self) -> ActionResult:
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
        association = self.ec2_backend.associate_vpc_cidr_block(
            vpc_id, cidr_block, amazon_provided_ipv6_cidr_blocks
        )
        result = {"VpcId": vpc_id}
        if not amazon_provided_ipv6_cidr_blocks:
            result["CidrBlockAssociation"] = {
                "AssociationId": association["association_id"],
                "CidrBlock": association["cidr_block"],
                "CidrBlockState": {"State": "associating"},
            }
        else:
            result["Ipv6CidrBlockAssociation"] = {
                "AssociationId": association["association_id"],
                "Ipv6CidrBlock": association["cidr_block"],
                "Ipv6CidrBlockState": {"State": "associating"},
                "NetworkBorderGroup": association.get(
                    "ipv6_cidr_block_network_border_group"
                ),
                "Ipv6Pool": association.get("ipv6_pool"),
            }
        return ActionResult(result)

    def disassociate_vpc_cidr_block(self) -> ActionResult:
        association_id = self._get_param("AssociationId")
        association = self.ec2_backend.disassociate_vpc_cidr_block(association_id)
        result = {"VpcId": association["vpc_id"]}
        if "::" in association.get("cidr_block", ""):
            result["Ipv6CidrBlockAssociation"] = {
                "AssociationId": association["association_id"],
                "Ipv6CidrBlock": association["cidr_block"],
                "Ipv6CidrBlockState": {"State": "disassociating"},
                "NetworkBorderGroup": association.get(
                    "ipv6_cidr_block_network_border_group"
                ),
                "Ipv6Pool": association.get("ipv6_pool"),
            }
        else:
            result["CidrBlockAssociation"] = {
                "AssociationId": association["association_id"],
                "CidrBlock": association["cidr_block"],
                "CidrBlockState": {"State": "disassociating"},
            }
        return ActionResult(result)

    def create_vpc_endpoint(self) -> ActionResult:
        vpc_id = self._get_param("VpcId")
        service_name = self._get_param("ServiceName")
        route_table_ids = self._get_param("RouteTableIds", [])
        subnet_ids = self._get_param("SubnetIds", [])
        endpoint_type = self._get_param("VpcEndpointType")
        policy_document = self._get_param("PolicyDocument")
        client_token = self._get_param("ClientToken")
        private_dns_enabled = self._get_bool_param("PrivateDnsEnabled", if_none=True)
        security_group_ids = self._get_param("SecurityGroupIds", [])

        tags = add_tag_specification(self._get_param("TagSpecifications", []))
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
        result = {"VpcEndpoint": vpc_end_point}
        return ActionResult(result)

    def modify_vpc_endpoint(self) -> ActionResult:
        vpc_id = self._get_param("VpcEndpointId")
        add_subnets = self._get_param("AddSubnetIds", [])
        remove_subnets = self._get_param("RemoveSubnetIds", [])
        add_route_tables = self._get_param("AddRouteTableIds", [])
        remove_route_tables = self._get_param("RemoveRouteTableIds", [])
        policy_doc = self._get_param("PolicyDocument")
        add_security_groups = self._get_param("AddSecurityGroupIds", [])
        remove_security_groups = self._get_param("RemoveSecurityGroupIds", [])
        self.ec2_backend.modify_vpc_endpoint(
            vpc_id=vpc_id,
            policy_doc=policy_doc,
            add_subnets=add_subnets,
            remove_subnets=remove_subnets,
            add_route_tables=add_route_tables,
            remove_route_tables=remove_route_tables,
            add_security_groups=add_security_groups,
            remove_security_groups=remove_security_groups,
        )
        return ActionResult({"Return": True})

    def describe_vpc_endpoint_services(self) -> ActionResult:
        vpc_end_point_services = self.ec2_backend.describe_vpc_endpoint_services(
            service_names=self._get_param("ServiceNames", []),
            filters=self._get_param("Filters", []),
            max_results=self._get_int_param("MaxResults"),
            next_token=self._get_param("NextToken"),
            region=self.region,
        )
        result = {
            "ServiceNames": vpc_end_point_services["serviceNames"],
            "ServiceDetails": vpc_end_point_services["servicesDetails"],
            "NextToken": vpc_end_point_services["nextToken"],
        }
        return ActionResult(result)

    def describe_vpc_endpoints(self) -> ActionResult:
        vpc_end_points_ids = self._get_param("VpcEndpointIds", [])
        filters = self._filters_from_querystring()
        vpc_end_points = self.ec2_backend.describe_vpc_endpoints(
            vpc_end_point_ids=vpc_end_points_ids, filters=filters
        )
        result = {"VpcEndpoints": vpc_end_points}
        return ActionResult(result)

    def delete_vpc_endpoints(self) -> ActionResult:
        vpc_end_points_ids = self._get_param("VpcEndpointIds", [])
        self.ec2_backend.delete_vpc_endpoints(vpce_ids=vpc_end_points_ids)
        return ActionResult({"Unsuccessful": []})

    def create_managed_prefix_list(self) -> ActionResult:
        address_family = self._get_param("AddressFamily")
        max_entries = self._get_param("MaxEntries")
        prefix_list_name = self._get_param("PrefixListName")
        entry = self._get_param("Entries", [])

        tags = self._parse_tag_specification().get("prefix-list", {})

        managed_prefix_list = self.ec2_backend.create_managed_prefix_list(
            address_family=address_family,
            entry=entry,
            max_entries=max_entries,
            prefix_list_name=prefix_list_name,
            tags=tags,
        )
        result = {"PrefixList": managed_prefix_list}
        return ActionResult(result)

    def describe_managed_prefix_lists(self) -> ActionResult:
        prefix_list_ids = self._get_param("PrefixListIds", [])
        filters = self._filters_from_querystring()
        managed_prefix_lists = self.ec2_backend.describe_managed_prefix_lists(
            prefix_list_ids=prefix_list_ids, filters=filters
        )
        result = {"PrefixLists": managed_prefix_lists}
        return ActionResult(result)

    def get_managed_prefix_list_entries(self) -> ActionResult:
        prefix_list_id = self._get_param("PrefixListId")
        target_version = self._get_param("TargetVersion")
        managed_prefix_list = self.ec2_backend.get_managed_prefix_list_entries(
            prefix_list_id=prefix_list_id
        )
        entries: list[ManagedPrefixList] = []
        if managed_prefix_list:
            entries = (
                list(managed_prefix_list.entries.values())[-1]
                if managed_prefix_list.entries.values()
                else []
            )
            if target_version:
                target_version = int(target_version)
                entries = managed_prefix_list.entries.get(target_version)
        result = {"Entries": entries}
        return ActionResult(result)

    def delete_managed_prefix_list(self) -> ActionResult:
        prefix_list_id = self._get_param("PrefixListId")
        managed_prefix_list = self.ec2_backend.delete_managed_prefix_list(
            prefix_list_id
        )
        result = {"PrefixList": managed_prefix_list}
        return ActionResult(result)

    def describe_prefix_lists(self) -> ActionResult:
        prefix_list_ids = self._get_param("PrefixListIds", [])
        filters = self._filters_from_querystring()
        managed_pls = self.ec2_backend.describe_managed_prefix_lists(
            prefix_list_ids=prefix_list_ids, filters=filters
        )
        result = {
            "PrefixLists": [
                pl
                for pl in managed_pls
                if pl.prefix_list_name.startswith("com.amazonaws.")
            ]
        }
        return ActionResult(result)

    def modify_managed_prefix_list(self) -> ActionResult:
        add_entry = self._get_param("AddEntries", [])
        prefix_list_id = self._get_param("PrefixListId")
        current_version = self._get_param("CurrentVersion")
        prefix_list_name = self._get_param("PrefixListName")
        remove_entry = self._get_param("RemoveEntries", [])

        current_version = int(current_version) if current_version else None

        managed_prefix_list = self.ec2_backend.modify_managed_prefix_list(
            add_entry=add_entry,
            prefix_list_id=prefix_list_id,
            current_version=current_version,
            prefix_list_name=prefix_list_name,
            remove_entry=remove_entry,
        )
        result = {"PrefixList": managed_prefix_list}
        return ActionResult(result)
