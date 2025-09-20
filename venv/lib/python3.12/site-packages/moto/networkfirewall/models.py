"""NetworkFirewallBackend class with methods for supported APIs."""

from typing import Any, Dict, List, Optional

from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import BaseModel
from moto.es.exceptions import ResourceNotFound
from moto.utilities.paginator import paginate
from moto.utilities.tagging_service import TaggingService

PAGINATION_MODEL = {
    "list_firewalls": {
        "input_token": "next_token",
        "limit_key": "max_results",
        "limit_default": 100,
        "unique_attribute": "arn",
    },
}


class NetworkFirewallModel(BaseModel):
    def __init__(
        self,
        account_id: str,
        region_name: str,
        firewall_name: str,
        firewall_policy_arn: str,
        vpc_id: str,
        subnet_mappings: List[str],
        delete_protection: bool,
        subnet_change_protection: bool,
        firewall_policy_change_protection: bool,
        description: str,
        tags: List[Dict[str, str]],
        encryption_configuration: Dict[str, str],
        enabled_analysis_types: List[str],
    ):
        self.firewall_name = firewall_name
        self.firewall_policy_arn = firewall_policy_arn
        self.vpc_id = vpc_id
        self.subnet_mappings = subnet_mappings
        self.delete_protection: bool = (
            delete_protection if delete_protection is not None else True
        )
        self.subnet_change_protection: bool = (
            subnet_change_protection if subnet_change_protection is not None else True
        )
        self.firewall_policy_change_protection: bool = (
            firewall_policy_change_protection
            if firewall_policy_change_protection is not None
            else True
        )
        self.description = description
        self.tags = tags
        self.encryption_configuration = encryption_configuration
        self.enabled_analysis_types = enabled_analysis_types

        self.arn = f"arn:aws:network-firewall:{region_name}:{account_id}:firewall/{self.firewall_name}"

        self.update_token = "1a2b3c4d-5e6f-7a8b-9c0d-1e2f3a4b5c6d"
        self.firewall_status = {
            "Status": "READY",
            "ConfigurationSyncStateSummary": "IN_SYNC",
        }
        self.logging_configs: Dict[str, List[Dict[str, Any]]] = {}

    def to_dict(self) -> Dict[str, Any]:
        return {
            "FirewallName": self.firewall_name,
            "FirewallArn": self.arn,
            "FirewallPolicyArn": self.firewall_policy_arn,
            "VpcId": self.vpc_id,
            "SubnetMappings": self.subnet_mappings,
            "DeleteProtection": self.delete_protection,
            "SubnetChangeProtection": self.subnet_change_protection,
            "FirewallPolicyChangeProtection": self.firewall_policy_change_protection,
            "Description": self.description,
            "Tags": self.tags,
            "EncryptionConfiguration": self.encryption_configuration,
            "EnabledAnalysisTypes": self.enabled_analysis_types,
        }


class NetworkFirewallBackend(BaseBackend):
    """Implementation of NetworkFirewall APIs."""

    def __init__(self, region_name: str, account_id: str) -> None:
        super().__init__(region_name, account_id)
        self.firewalls: Dict[str, NetworkFirewallModel] = {}
        self.tagger = TaggingService()

    def create_firewall(
        self,
        firewall_name: str,
        firewall_policy_arn: str,
        vpc_id: str,
        subnet_mappings: List[str],
        delete_protection: bool,
        subnet_change_protection: bool,
        firewall_policy_change_protection: bool,
        description: str,
        tags: List[Dict[str, str]],
        encryption_configuration: Dict[str, str],
        enabled_analysis_types: List[str],
    ) -> NetworkFirewallModel:
        firewall = NetworkFirewallModel(
            self.account_id,
            self.region_name,
            firewall_name=firewall_name,
            firewall_policy_arn=firewall_policy_arn,
            vpc_id=vpc_id,
            subnet_mappings=subnet_mappings,
            delete_protection=delete_protection,
            subnet_change_protection=subnet_change_protection,
            firewall_policy_change_protection=firewall_policy_change_protection,
            description=description,
            tags=tags,
            encryption_configuration=encryption_configuration,
            enabled_analysis_types=enabled_analysis_types,
        )
        self.firewalls[firewall.arn] = firewall

        if tags:
            self.tagger.tag_resource(firewall.arn, tags)

        return firewall

    def _get_firewall(
        self, firewall_arn: Optional[str], firewall_name: Optional[str]
    ) -> NetworkFirewallModel:
        if firewall_arn:
            if firewall_arn in self.firewalls:
                return self.firewalls[firewall_arn]
        for firewall in self.firewalls.values():
            if firewall.firewall_name == firewall_name:
                return firewall
        raise ResourceNotFound("NetworkFirewall", str(firewall_arn or firewall_name))

    def describe_logging_configuration(
        self, firewall_arn: str, firewall_name: str
    ) -> NetworkFirewallModel:
        firewall: NetworkFirewallModel = self._get_firewall(firewall_arn, firewall_name)
        return firewall

    def update_logging_configuration(
        self,
        firewall_arn: str,
        firewall_name: str,
        logging_configuration: Dict[str, List[Dict[str, Any]]],
    ) -> NetworkFirewallModel:
        firewall: NetworkFirewallModel = self._get_firewall(firewall_arn, firewall_name)
        firewall.logging_configs = logging_configuration
        return firewall

    @paginate(pagination_model=PAGINATION_MODEL)
    def list_firewalls(self, vpc_ids: List[str]) -> List[NetworkFirewallModel]:
        firewalls = list(self.firewalls.values())
        if vpc_ids:
            firewalls = [fw for fw in firewalls if fw.vpc_id in vpc_ids]
        return firewalls

    def describe_firewall(
        self, firewall_name: str, firewall_arn: str
    ) -> NetworkFirewallModel:
        return self._get_firewall(firewall_arn, firewall_name)


networkfirewall_backends = BackendDict(NetworkFirewallBackend, "network-firewall")
