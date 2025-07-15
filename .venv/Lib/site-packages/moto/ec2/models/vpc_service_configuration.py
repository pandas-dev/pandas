from typing import Any, Dict, List, Optional

from moto.core.common_models import CloudFormationModel
from moto.moto_api._internal import mock_random

from ..exceptions import UnknownVpcEndpointService
from .core import TaggedEC2Resource


class VPCServiceConfiguration(TaggedEC2Resource, CloudFormationModel):
    def __init__(
        self,
        load_balancers: List[Any],
        region: str,
        acceptance_required: bool,
        private_dns_name: str,
        ec2_backend: Any,
    ):
        self.id = f"vpce-svc-{mock_random.get_random_hex(length=8)}"
        self.service_name = f"com.amazonaws.vpce.{region}.{self.id}"
        self.service_state = "Available"

        self.availability_zones = []
        for lb in load_balancers:
            for subnet in lb.subnets:
                self.availability_zones.append(subnet.availability_zone)

        self.gateway_load_balancer_arns = []
        self.network_load_balancer_arns = []
        for lb in load_balancers:
            if lb.loadbalancer_type == "network":
                self.service_type = "Interface"
                self.network_load_balancer_arns.append(lb.arn)
            else:
                self.service_type = "Gateway"
                self.gateway_load_balancer_arns.append(lb.arn)

        self.acceptance_required = acceptance_required
        self.manages_vpc_endpoints = False
        self.private_dns_name = private_dns_name
        self.endpoint_dns_name = f"{self.id}.{region}.vpce.amazonaws.com"

        self.principals: List[str] = []
        self.ec2_backend = ec2_backend

    def to_dict(self) -> Dict[str, Any]:
        return {
            "AcceptanceRequired": self.acceptance_required,
            "AvailabilityZones": self.availability_zones,
            "BaseEndpointDnsNames": [self.endpoint_dns_name],
            "ManagesVpcEndpoints": self.manages_vpc_endpoints,
            "Owner": self.ec2_backend.account_id,
            "PrivateDnsName": self.private_dns_name,
            "PrivateDnsNames": [{"PrivateDnsName": self.private_dns_name}],
            "ServiceId": self.id,
            "ServiceName": self.service_name,
            "ServiceType": [{"ServiceType": self.service_type}],
            "VpcEndpointPolicySupported": True,
        }


class VPCServiceConfigurationBackend:
    def __init__(self) -> None:
        self.configurations: Dict[str, VPCServiceConfiguration] = {}

    @property
    def elbv2_backend(self) -> Any:  # type: ignore[misc]
        from moto.elbv2.models import elbv2_backends

        return elbv2_backends[self.account_id][self.region_name]  # type: ignore[attr-defined]

    def get_vpc_endpoint_service(
        self, resource_id: str
    ) -> Optional[VPCServiceConfiguration]:
        return self.configurations.get(resource_id)

    def create_vpc_endpoint_service_configuration(
        self,
        lb_arns: List[Any],
        acceptance_required: bool,
        private_dns_name: str,
        tags: List[Dict[str, str]],
    ) -> VPCServiceConfiguration:
        lbs = self.elbv2_backend.describe_load_balancers(arns=lb_arns, names=None)
        config = VPCServiceConfiguration(
            load_balancers=lbs,
            region=self.region_name,  # type: ignore[attr-defined]
            acceptance_required=acceptance_required,
            private_dns_name=private_dns_name,
            ec2_backend=self,
        )
        for tag in tags or []:
            config.add_tag(tag["Key"], tag["Value"])

        self.configurations[config.id] = config
        return config

    def describe_vpc_endpoint_service_configurations(
        self, service_ids: Optional[List[str]]
    ) -> List[VPCServiceConfiguration]:
        """
        The Filters, MaxResults, NextToken parameters are not yet implemented
        """
        if service_ids:
            found_configs = []
            for service_id in service_ids:
                if service_id in self.configurations:
                    found_configs.append(self.configurations[service_id])
                else:
                    raise UnknownVpcEndpointService(service_id)
            return found_configs
        return list(self.configurations.values())

    def delete_vpc_endpoint_service_configurations(
        self, service_ids: List[str]
    ) -> List[str]:
        missing = [s for s in service_ids if s not in self.configurations]
        for s in service_ids:
            self.configurations.pop(s, None)
        return missing

    def describe_vpc_endpoint_service_permissions(self, service_id: str) -> List[str]:
        """
        The Filters, MaxResults, NextToken parameters are not yet implemented
        """
        config = self.describe_vpc_endpoint_service_configurations([service_id])[0]
        return config.principals

    def modify_vpc_endpoint_service_permissions(
        self, service_id: str, add_principals: List[str], remove_principals: List[str]
    ) -> None:
        config = self.describe_vpc_endpoint_service_configurations([service_id])[0]
        config.principals += add_principals
        config.principals = [p for p in config.principals if p not in remove_principals]
        config.principals = list(set(config.principals))

    def modify_vpc_endpoint_service_configuration(
        self,
        service_id: str,
        acceptance_required: Optional[str],
        private_dns_name: Optional[str],
        add_network_lbs: List[str],
        remove_network_lbs: List[str],
        add_gateway_lbs: List[str],
        remove_gateway_lbs: List[str],
    ) -> None:
        """
        The following parameters are not yet implemented: RemovePrivateDnsName
        """
        config = self.describe_vpc_endpoint_service_configurations([service_id])[0]
        if private_dns_name is not None:
            config.private_dns_name = private_dns_name
        if acceptance_required is not None:
            config.acceptance_required = str(acceptance_required).lower() == "true"
        for lb in add_network_lbs:
            config.network_load_balancer_arns.append(lb)
        for lb in remove_network_lbs:
            config.network_load_balancer_arns.remove(lb)
        for lb in add_gateway_lbs:
            config.gateway_load_balancer_arns.append(lb)
        for lb in remove_gateway_lbs:
            config.gateway_load_balancer_arns.remove(lb)
