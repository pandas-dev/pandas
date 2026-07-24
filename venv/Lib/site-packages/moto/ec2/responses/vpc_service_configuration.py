from moto.core.responses import ActionResult

from ..exceptions import NoLoadBalancersProvided
from ._base_response import EC2BaseResponse


class VPCEndpointServiceConfiguration(EC2BaseResponse):
    def create_vpc_endpoint_service_configuration(self) -> ActionResult:
        gateway_lbs = self._get_param("GatewayLoadBalancerArns", [])
        network_lbs = self._get_param("NetworkLoadBalancerArns", [])
        if not gateway_lbs and not network_lbs:
            raise NoLoadBalancersProvided

        tags = self._get_param("TagSpecifications", [])
        if tags:
            tags = tags[0].get("Tags")
        acceptance_required = self._get_param("AcceptanceRequired", True)
        private_dns_name = self._get_param("PrivateDnsName")
        supported_regions = self._get_param("SupportedRegions", [])

        config = self.ec2_backend.create_vpc_endpoint_service_configuration(
            gateway_lbs or network_lbs,
            acceptance_required=acceptance_required,
            private_dns_name=private_dns_name,
            tags=tags,
            supported_regions=supported_regions,
        )
        return ActionResult({"ServiceConfiguration": config})

    def describe_vpc_endpoint_service_configurations(self) -> ActionResult:
        service_ids = self._get_param("ServiceIds", [])
        configs = self.ec2_backend.describe_vpc_endpoint_service_configurations(
            service_ids
        )
        return ActionResult({"ServiceConfigurations": configs})

    def delete_vpc_endpoint_service_configurations(self) -> ActionResult:
        service_ids = self._get_param("ServiceIds", [])
        missing_configs = self.ec2_backend.delete_vpc_endpoint_service_configurations(
            service_ids
        )
        unsuccessful = [
            {
                "ResourceId": m,
                "Error": {
                    "Code": "InvalidVpcEndpointService.NotFound",
                    "Message": f"The VpcEndpointService Id '{m}' does not exist",
                },
            }
            for m in missing_configs
        ]
        return ActionResult({"Unsuccessful": unsuccessful})

    def describe_vpc_endpoint_service_permissions(self) -> ActionResult:
        service_id = self._get_param("ServiceId")
        principals = self.ec2_backend.describe_vpc_endpoint_service_permissions(
            service_id
        )
        allowed_principals = [{"Principal": p} for p in principals]
        return ActionResult({"AllowedPrincipals": allowed_principals})

    def modify_vpc_endpoint_service_configuration(self) -> ActionResult:
        service_id = self._get_param("ServiceId")
        private_dns_name = self._get_param("PrivateDnsName")
        acceptance_required = self._get_param("AcceptanceRequired")
        add_network_lbs = self._get_param("AddNetworkLoadBalancerArns", [])
        remove_network_lbs = self._get_param("RemoveNetworkLoadBalancerArns", [])
        add_gateway_lbs = self._get_param("AddGatewayLoadBalancerArns", [])
        remove_gateway_lbs = self._get_param("RemoveGatewayLoadBalancerArns", [])
        add_supported_regions = self._get_param("AddSupportedRegions", [])
        remove_supported_regions = self._get_param("RemoveSupportedRegions", [])

        self.ec2_backend.modify_vpc_endpoint_service_configuration(
            service_id,
            acceptance_required=acceptance_required,
            private_dns_name=private_dns_name,
            add_network_lbs=add_network_lbs,
            remove_network_lbs=remove_network_lbs,
            add_gateway_lbs=add_gateway_lbs,
            remove_gateway_lbs=remove_gateway_lbs,
            add_supported_regions=add_supported_regions,
            remove_supported_regions=remove_supported_regions,
        )
        return ActionResult({"Return": True})

    def modify_vpc_endpoint_service_permissions(self) -> ActionResult:
        service_id = self._get_param("ServiceId")
        add_principals = self._get_param("AddAllowedPrincipals", [])
        remove_principals = self._get_param("RemoveAllowedPrincipals", [])

        self.ec2_backend.modify_vpc_endpoint_service_permissions(
            service_id, add_principals, remove_principals
        )
        return ActionResult({"ReturnValue": True})
