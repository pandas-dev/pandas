"""Handles incoming networkfirewall requests, invokes methods, returns responses."""

import json

from moto.core.responses import BaseResponse

from .models import NetworkFirewallBackend, networkfirewall_backends


class NetworkFirewallResponse(BaseResponse):
    """Handler for NetworkFirewall requests and responses."""

    def __init__(self) -> None:
        super().__init__(service_name="networkfirewall")

    @property
    def networkfirewall_backend(self) -> NetworkFirewallBackend:
        """Return backend instance specific for this region."""
        return networkfirewall_backends[self.current_account][self.region]

    def create_firewall(self) -> str:
        firewall_name = self._get_param("FirewallName")
        firewall_policy_arn = self._get_param("FirewallPolicyArn")
        vpc_id = self._get_param("VpcId")
        subnet_mappings = self._get_param("SubnetMappings")
        delete_protection = self._get_param("DeleteProtection")
        subnet_change_protection = self._get_param("SubnetChangeProtection")
        firewall_policy_change_protection = self._get_param(
            "FirewallPolicyChangeProtection"
        )
        description = self._get_param("Description")
        tags = self._get_param("Tags")
        encryption_configuration = self._get_param("EncryptionConfiguration")
        enabled_analysis_types = self._get_param("EnabledAnalysisTypes")
        firewall = self.networkfirewall_backend.create_firewall(
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

        return json.dumps(
            dict(Firewall=firewall.to_dict(), FirewallStatus=firewall.firewall_status)
        )

    def describe_logging_configuration(self) -> str:
        firewall_arn = self._get_param("FirewallArn")
        firewall_name = self._get_param("FirewallName")
        firewall = self.networkfirewall_backend.describe_logging_configuration(
            firewall_arn=firewall_arn,
            firewall_name=firewall_name,
        )
        return json.dumps(
            dict(
                FirewallArn=firewall.arn,
                LoggingConfiguration=firewall.logging_configs,
            )
        )

    def update_logging_configuration(self) -> str:
        firewall_arn = self._get_param("FirewallArn")
        firewall_name = self._get_param("FirewallName")
        logging_configuration = self._get_param("LoggingConfiguration")
        firewall = self.networkfirewall_backend.update_logging_configuration(
            firewall_arn=firewall_arn,
            firewall_name=firewall_name,
            logging_configuration=logging_configuration,
        )
        return json.dumps(
            dict(
                FirewallArn=firewall.arn,
                FirewallName=firewall.firewall_name,
                LoggingConfiguration=firewall.logging_configs,
            )
        )

    def list_firewalls(self) -> str:
        next_token = self._get_param("NextToken")
        vpc_ids = self._get_param("VpcIds")
        max_results = self._get_param("MaxResults")
        firewalls, next_token = self.networkfirewall_backend.list_firewalls(
            next_token=next_token,
            vpc_ids=vpc_ids,
            max_results=max_results,
        )
        firewall_list = [fw.to_dict() for fw in firewalls]
        return json.dumps(dict(nextToken=next_token, Firewalls=firewall_list))

    def describe_firewall(self) -> str:
        firewall_name = self._get_param("FirewallName")
        firewall_arn = self._get_param("FirewallArn")
        firewall = self.networkfirewall_backend.describe_firewall(
            firewall_name=firewall_name,
            firewall_arn=firewall_arn,
        )
        return json.dumps(
            dict(
                UpdateToken=firewall.update_token,
                Firewall=firewall.to_dict(),
                FirewallStatus=firewall.firewall_status,
            )
        )
