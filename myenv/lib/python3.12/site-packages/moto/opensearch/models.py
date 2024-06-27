from typing import Any, Dict, List, Optional

from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import BaseModel
from moto.utilities.tagging_service import TaggingService
from moto.utilities.utils import get_partition

from .data import compatible_versions
from .exceptions import EngineTypeNotFoundException, ResourceNotFoundException

default_cluster_config = {
    "InstanceType": "t3.small.search",
    "InstanceCount": 1,
    "DedicatedMasterEnabled": False,
    "ZoneAwarenessEnabled": False,
    "WarmEnabled": False,
    "ColdStorageOptions": {"Enabled": False},
}
default_advanced_security_options = {
    "Enabled": False,
    "InternalUserDatabaseEnabled": False,
    "AnonymousAuthEnabled": False,
}
default_domain_endpoint_options = {
    "EnforceHTTPS": False,
    "TLSSecurityPolicy": "Policy-Min-TLS-1-0-2019-07",
    "CustomEndpointEnabled": False,
}
default_software_update_options = {
    "CurrentVersion": "",
    "NewVersion": "",
    "UpdateAvailable": False,
    "Cancellable": False,
    "UpdateStatus": "COMPLETED",
    "Description": "There is no software update available for this domain.",
    "AutomatedUpdateDate": "1969-12-31T23:00:00-01:00",
    "OptionalDeployment": True,
}
default_advanced_options = {
    "override_main_response_version": "false",
    "rest.action.multi.allow_explicit_index": "true",
}


class OpenSearchDomain(BaseModel):
    def __init__(
        self,
        account_id: str,
        region: str,
        domain_name: str,
        engine_version: str,
        cluster_config: Dict[str, Any],
        ebs_options: Dict[str, Any],
        access_policies: str,
        snapshot_options: Dict[str, int],
        vpc_options: Dict[str, List[str]],
        cognito_options: Dict[str, Any],
        encryption_at_rest_options: Dict[str, Any],
        node_to_node_encryption_options: Dict[str, bool],
        advanced_options: Dict[str, str],
        log_publishing_options: Dict[str, Any],
        domain_endpoint_options: Dict[str, Any],
        advanced_security_options: Dict[str, Any],
        auto_tune_options: Dict[str, Any],
        off_peak_window_options: Dict[str, Any],
        software_update_options: Dict[str, bool],
    ):
        self.domain_id = f"{account_id}/{domain_name}"
        self.domain_name = domain_name
        self.arn = (
            f"arn:{get_partition(region)}:es:{region}:{account_id}:domain/{domain_name}"
        )
        self.engine_version = engine_version or "OpenSearch 2.5"
        self.cluster_config = cluster_config or {}
        self.ebs_options = ebs_options or {"EBSEnabled": False}
        self.access_policies = access_policies or ""
        self.snapshot_options = snapshot_options or {"AutomatedSnapshotStartHour": 0}
        self.vpc_options = vpc_options
        self.cognito_options = cognito_options or {"Enabled": False}
        self.encryption_at_rest_options = encryption_at_rest_options or {
            "Enabled": False
        }
        self.node_to_node_encryption_options = node_to_node_encryption_options or {
            "Enabled": False
        }
        self.advanced_options = advanced_options or default_advanced_options
        self.log_publishing_options = log_publishing_options
        self.domain_endpoint_options = (
            domain_endpoint_options or default_domain_endpoint_options
        )
        self.advanced_security_options = (
            advanced_security_options or default_advanced_security_options
        )
        self.auto_tune_options = auto_tune_options or {"State": "ENABLE_IN_PROGRESS"}
        self.off_peak_windows_options = off_peak_window_options
        self.software_update_options = (
            software_update_options or default_software_update_options
        )

        self.deleted = False
        self.processing = False

        # Defaults
        for key, value in default_cluster_config.items():
            if key not in self.cluster_config:
                self.cluster_config[key] = value

        if self.vpc_options is None:
            self.endpoint: Optional[str] = f"{domain_name}.{region}.es.amazonaws.com"
            self.endpoints: Optional[Dict[str, str]] = None
        else:
            self.endpoint = None
            self.endpoints = {"vpc": f"{domain_name}.{region}.es.amazonaws.com"}

    def delete(self) -> None:
        self.deleted = True
        self.processing = True

    def dct_options(self) -> Dict[str, Any]:
        return {
            "Endpoint": self.endpoint,
            "Endpoints": self.endpoints,
            "EngineVersion": self.engine_version,
            "ClusterConfig": self.cluster_config,
            "EBSOptions": self.ebs_options,
            "AccessPolicies": self.access_policies,
            "SnapshotOptions": self.snapshot_options,
            "VPCOptions": self.vpc_options,
            "CognitoOptions": self.cognito_options,
            "EncryptionAtRestOptions": self.encryption_at_rest_options,
            "NodeToNodeEncryptionOptions": self.node_to_node_encryption_options,
            "AdvancedOptions": self.advanced_options,
            "LogPublishingOptions": self.log_publishing_options,
            "DomainEndpointOptions": self.domain_endpoint_options,
            "AdvancedSecurityOptions": self.advanced_security_options,
            "AutoTuneOptions": self.auto_tune_options,
            "OffPeakWindowsOptions": self.off_peak_windows_options,
            "SoftwareUpdateOptions": self.software_update_options,
        }

    def to_dict(self) -> Dict[str, Any]:
        dct = {
            "DomainId": self.domain_id,
            "DomainName": self.domain_name,
            "ARN": self.arn,
            "Created": True,
            "Deleted": self.deleted,
            "Processing": self.processing,
            "UpgradeProcessing": False,
        }
        for key, value in self.dct_options().items():
            if value is not None:
                dct[key] = value
        return dct

    def to_config_dict(self) -> Dict[str, Any]:
        dct: Dict[str, Any] = dict()
        for key, value in self.dct_options().items():
            if value is not None:
                dct[key] = {"Options": value}
        return dct

    def update(
        self,
        cluster_config: Dict[str, Any],
        ebs_options: Dict[str, Any],
        access_policies: str,
        snapshot_options: Dict[str, int],
        vpc_options: Dict[str, List[str]],
        cognito_options: Dict[str, Any],
        encryption_at_rest_options: Dict[str, Any],
        node_to_node_encryption_options: Dict[str, bool],
        advanced_options: Dict[str, str],
        log_publishing_options: Dict[str, Any],
        domain_endpoint_options: Dict[str, Any],
        advanced_security_options: Dict[str, Any],
        auto_tune_options: Dict[str, Any],
        off_peak_window_options: Dict[str, Any],
        software_update_options: Dict[str, bool],
    ) -> None:
        self.cluster_config = cluster_config or self.cluster_config
        self.ebs_options = ebs_options or self.ebs_options
        self.access_policies = access_policies or self.access_policies
        self.snapshot_options = snapshot_options or self.snapshot_options
        self.vpc_options = vpc_options or self.vpc_options
        self.cognito_options = cognito_options or self.cognito_options
        self.encryption_at_rest_options = (
            encryption_at_rest_options or self.encryption_at_rest_options
        )
        self.node_to_node_encryption_options = (
            node_to_node_encryption_options or self.node_to_node_encryption_options
        )
        self.advanced_options = advanced_options or self.advanced_options
        self.log_publishing_options = (
            log_publishing_options or self.log_publishing_options
        )
        self.domain_endpoint_options = (
            domain_endpoint_options or self.domain_endpoint_options
        )
        self.advanced_security_options = (
            advanced_security_options or self.advanced_security_options
        )
        self.auto_tune_options = auto_tune_options or self.auto_tune_options
        self.off_peak_windows_options = (
            off_peak_window_options or self.off_peak_windows_options
        )
        self.software_update_options = (
            software_update_options or self.software_update_options
        )


class OpenSearchServiceBackend(BaseBackend):
    """Implementation of OpenSearchService APIs."""

    def __init__(self, region_name: str, account_id: str):
        super().__init__(region_name, account_id)
        self.domains: Dict[str, OpenSearchDomain] = dict()
        self.tagger = TaggingService()

    def create_domain(
        self,
        domain_name: str,
        engine_version: str,
        cluster_config: Dict[str, Any],
        ebs_options: Dict[str, Any],
        access_policies: str,
        snapshot_options: Dict[str, Any],
        vpc_options: Dict[str, Any],
        cognito_options: Dict[str, Any],
        encryption_at_rest_options: Dict[str, Any],
        node_to_node_encryption_options: Dict[str, Any],
        advanced_options: Dict[str, Any],
        log_publishing_options: Dict[str, Any],
        domain_endpoint_options: Dict[str, Any],
        advanced_security_options: Dict[str, Any],
        tag_list: List[Dict[str, str]],
        auto_tune_options: Dict[str, Any],
        off_peak_window_options: Dict[str, Any],
        software_update_options: Dict[str, Any],
    ) -> OpenSearchDomain:
        domain = OpenSearchDomain(
            account_id=self.account_id,
            region=self.region_name,
            domain_name=domain_name,
            engine_version=engine_version,
            cluster_config=cluster_config,
            ebs_options=ebs_options,
            access_policies=access_policies,
            snapshot_options=snapshot_options,
            vpc_options=vpc_options,
            cognito_options=cognito_options,
            encryption_at_rest_options=encryption_at_rest_options,
            node_to_node_encryption_options=node_to_node_encryption_options,
            advanced_options=advanced_options,
            log_publishing_options=log_publishing_options,
            domain_endpoint_options=domain_endpoint_options,
            advanced_security_options=advanced_security_options,
            auto_tune_options=auto_tune_options,
            off_peak_window_options=off_peak_window_options,
            software_update_options=software_update_options,
        )
        self.domains[domain_name] = domain
        if tag_list:
            self.add_tags(domain.arn, tag_list)
        return domain

    def get_compatible_versions(self, domain_name: str) -> List[Dict[str, Any]]:
        if domain_name not in self.domains:
            raise ResourceNotFoundException(domain_name)
        return compatible_versions

    def delete_domain(self, domain_name: str) -> OpenSearchDomain:
        if domain_name not in self.domains:
            raise ResourceNotFoundException(domain_name)
        self.domains[domain_name].delete()
        return self.domains.pop(domain_name)

    def describe_domain(self, domain_name: str) -> OpenSearchDomain:
        if domain_name not in self.domains:
            raise ResourceNotFoundException(domain_name)
        return self.domains[domain_name]

    def describe_domain_config(self, domain_name: str) -> OpenSearchDomain:
        return self.describe_domain(domain_name)

    def update_domain_config(
        self,
        domain_name: str,
        cluster_config: Dict[str, Any],
        ebs_options: Dict[str, Any],
        access_policies: str,
        snapshot_options: Dict[str, Any],
        vpc_options: Dict[str, Any],
        cognito_options: Dict[str, Any],
        encryption_at_rest_options: Dict[str, Any],
        node_to_node_encryption_options: Dict[str, Any],
        advanced_options: Dict[str, Any],
        log_publishing_options: Dict[str, Any],
        domain_endpoint_options: Dict[str, Any],
        advanced_security_options: Dict[str, Any],
        auto_tune_options: Dict[str, Any],
        off_peak_window_options: Dict[str, Any],
        software_update_options: Dict[str, Any],
    ) -> OpenSearchDomain:
        domain = self.describe_domain(domain_name)
        domain.update(
            cluster_config=cluster_config,
            ebs_options=ebs_options,
            access_policies=access_policies,
            snapshot_options=snapshot_options,
            vpc_options=vpc_options,
            cognito_options=cognito_options,
            encryption_at_rest_options=encryption_at_rest_options,
            node_to_node_encryption_options=node_to_node_encryption_options,
            advanced_options=advanced_options,
            log_publishing_options=log_publishing_options,
            domain_endpoint_options=domain_endpoint_options,
            advanced_security_options=advanced_security_options,
            auto_tune_options=auto_tune_options,
            off_peak_window_options=off_peak_window_options,
            software_update_options=software_update_options,
        )
        return domain

    def add_tags(self, arn: str, tags: List[Dict[str, str]]) -> None:
        self.tagger.tag_resource(arn, tags)

    def list_tags(self, arn: str) -> List[Dict[str, str]]:
        return self.tagger.list_tags_for_resource(arn)["Tags"]

    def remove_tags(self, arn: str, tag_keys: List[str]) -> None:
        self.tagger.untag_resource_using_names(arn, tag_keys)

    def list_domain_names(self, engine_type: str) -> List[Dict[str, str]]:
        domains = []
        for domain in self.domains.values():
            if engine_type:
                if engine_type in domain.engine_version:
                    domains.append(
                        {
                            "DomainName": domain.domain_name,
                            "EngineType": engine_type.split("_")[0],
                        }
                    )
                else:
                    raise EngineTypeNotFoundException(domain.domain_name)
            else:
                domains.append(
                    {
                        "DomainName": domain.domain_name,
                        "EngineType": domain.engine_version.split("_")[0],
                    }
                )
        return domains


opensearch_backends = BackendDict(OpenSearchServiceBackend, "opensearch")
