import datetime
from typing import Any, Optional

from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import BaseModel
from moto.core.utils import unix_time
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


class EngineVersion(BaseModel):
    def __init__(self, options: str, create_time: datetime.datetime) -> None:
        self.options = options or "OpenSearch_2.5"
        self.create_time = unix_time(create_time)
        self.update_time = self.create_time

    def to_dict(self) -> dict[str, Any]:
        return {
            "Options": self.options,
            "Status": {
                "CreationDate": self.create_time,
                "PendingDeletion": False,
                "State": "Active",
                "UpdateDate": self.update_time,
                "UpdateVersion": 28,
            },
        }


class OpenSearchDomain(BaseModel):
    def __init__(
        self,
        account_id: str,
        region: str,
        domain_name: str,
        engine_version: str,
        cluster_config: dict[str, Any],
        ebs_options: dict[str, Any],
        access_policies: str,
        snapshot_options: dict[str, int],
        vpc_options: dict[str, list[str]],
        cognito_options: dict[str, Any],
        encryption_at_rest_options: dict[str, Any],
        node_to_node_encryption_options: dict[str, bool],
        advanced_options: dict[str, str],
        log_publishing_options: dict[str, Any],
        domain_endpoint_options: dict[str, Any],
        advanced_security_options: dict[str, Any],
        auto_tune_options: dict[str, Any],
        off_peak_window_options: dict[str, Any],
        software_update_options: dict[str, bool],
        is_es: bool,
        elasticsearch_version: Optional[str],
        elasticsearch_cluster_config: Optional[str],
    ):
        # Add creation_date attribute
        self.creation_date = unix_time(datetime.datetime.now())

        self.domain_id = f"{account_id}/{domain_name}"
        self.domain_name = domain_name
        self.arn = (
            f"arn:{get_partition(region)}:es:{region}:{account_id}:domain/{domain_name}"
        )
        self.engine_version = EngineVersion(engine_version, datetime.datetime.now())
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
        if not self.auto_tune_options.get("State"):
            self.auto_tune_options["State"] = "ENABLED"
        # Rename to singular everywhere
        self.off_peak_window_options = off_peak_window_options
        self.software_update_options = (
            software_update_options or default_software_update_options
        )
        self.engine_type = "Elasticsearch" if is_es else "OpenSearch"
        self.is_es = is_es
        self.elasticsearch_version = elasticsearch_version
        self.elasticsearch_cluster_config = elasticsearch_cluster_config

        self.deleted = False
        self.processing = False

        # Defaults
        for key, value in default_cluster_config.items():
            if key not in self.cluster_config:
                self.cluster_config[key] = value

        if self.vpc_options is None:
            self.endpoint: Optional[str] = f"{domain_name}.{region}.es.amazonaws.com"
            self.endpoints: Optional[dict[str, str]] = None
        else:
            self.endpoint = None
            self.endpoints = {"vpc": f"{domain_name}.{region}.es.amazonaws.com"}

    def delete(self) -> None:
        self.deleted = True
        self.processing = True

    def dct_options(self) -> dict[str, Any]:
        return {
            "Endpoint": self.endpoint,
            "Endpoints": self.endpoints,
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
            # Use singular key and attribute
            "OffPeakWindowOptions": self.off_peak_window_options,
            "SoftwareUpdateOptions": self.software_update_options,
            "ElasticsearchVersion": self.elasticsearch_version,
            "ElasticsearchClusterConfig": self.elasticsearch_cluster_config,
        }

    def to_dict(self) -> dict[str, Any]:
        dct = {
            "DomainId": self.domain_id,
            "DomainName": self.domain_name,
            "ARN": self.arn,
            "Created": True,
            "Deleted": self.deleted,
            "EngineVersion": self.engine_version.options,
            "Processing": self.processing,
            "UpgradeProcessing": False,
        }
        for key, value in self.dct_options().items():
            if value is not None:
                dct[key] = value
        return dct

    def _status_block(self) -> dict[str, Any]:
        return {
            "State": "Active",
            "CreationDate": self.creation_date,
            "UpdateDate": self.creation_date,
            "UpdateVersion": 1,
            "PendingDeletion": False,
        }

    def _wrap(self, options: Any) -> dict[str, Any]:
        # Most DomainConfig sections only need {"Options": ...}
        return {"Options": options}

    def to_config_dict(self) -> dict[str, Any]:
        cfg: dict[str, Any] = {}

        # Cluster config section (key differs for ES vs OS)
        cluster_key = "ElasticsearchClusterConfig" if self.is_es else "ClusterConfig"
        cluster_opts = (
            self.elasticsearch_cluster_config
            if (self.is_es and self.elasticsearch_cluster_config)
            else (self.cluster_config or default_cluster_config)
        )
        cfg[cluster_key] = self._wrap(cluster_opts)

        # Version section:
        # - OpenSearch expects Status for EngineVersion (use EngineVersion.to_dict()).
        # - ES expects only Options for ElasticsearchVersion.
        if self.is_es:
            cfg["ElasticsearchVersion"] = self._wrap(
                self.elasticsearch_version or self.engine_version.options
            )
        else:
            cfg["EngineVersion"] = self.engine_version.to_dict()

        # EBSOptions: default to minimal disabled if not provided
        ebs_opts = (
            self.ebs_options if self.ebs_options is not None else {"EBSEnabled": False}
        )
        cfg["EBSOptions"] = self._wrap(ebs_opts)

        # Node-to-node encryption: default minimal
        n2n_opts = (
            self.node_to_node_encryption_options
            if self.node_to_node_encryption_options is not None
            else {"Enabled": False}
        )
        cfg["NodeToNodeEncryptionOptions"] = self._wrap(n2n_opts)

        # Encryption at rest: default minimal
        ear_opts = (
            self.encryption_at_rest_options
            if self.encryption_at_rest_options is not None
            else {"Enabled": False}
        )
        cfg["EncryptionAtRestOptions"] = self._wrap(ear_opts)

        # Access policies: default empty string
        cfg["AccessPolicies"] = self._wrap(self.access_policies or "")

        # Optional passthrough sections
        if self.snapshot_options is not None:
            cfg["SnapshotOptions"] = self._wrap(self.snapshot_options)
        if self.vpc_options is not None:
            cfg["VPCOptions"] = self._wrap(self.vpc_options)
        if self.cognito_options is not None:
            cfg["CognitoOptions"] = self._wrap(self.cognito_options)
        if self.log_publishing_options is not None:
            cfg["LogPublishingOptions"] = self._wrap(self.log_publishing_options)
        if self.auto_tune_options is not None:
            cfg["AutoTuneOptions"] = self._wrap(self.auto_tune_options)
        if self.off_peak_window_options is not None:
            cfg["OffPeakWindowOptions"] = self._wrap(self.off_peak_window_options)

        # Always include with sensible defaults
        cfg["AdvancedOptions"] = self._wrap(
            self.advanced_options or default_advanced_options
        )
        cfg["DomainEndpointOptions"] = self._wrap(
            self.domain_endpoint_options or default_domain_endpoint_options
        )
        cfg["AdvancedSecurityOptions"] = self._wrap(
            self.advanced_security_options or default_advanced_security_options
        )
        cfg["SoftwareUpdateOptions"] = self._wrap(
            self.software_update_options or default_software_update_options
        )

        return cfg

    def update(
        self,
        cluster_config: dict[str, Any],
        ebs_options: dict[str, Any],
        access_policies: str,
        snapshot_options: dict[str, int],
        vpc_options: dict[str, list[str]],
        cognito_options: dict[str, Any],
        encryption_at_rest_options: dict[str, Any],
        node_to_node_encryption_options: dict[str, bool],
        advanced_options: dict[str, str],
        log_publishing_options: dict[str, Any],
        domain_endpoint_options: dict[str, Any],
        advanced_security_options: dict[str, Any],
        auto_tune_options: dict[str, Any],
        off_peak_window_options: dict[str, Any],
        software_update_options: dict[str, bool],
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
        # Fix attribute name (singular)
        self.off_peak_window_options = (
            off_peak_window_options or self.off_peak_window_options
        )
        self.software_update_options = (
            software_update_options or self.software_update_options
        )
        self.engine_version.update_time = unix_time(datetime.datetime.now())


class OpenSearchServiceBackend(BaseBackend):
    """Implementation of OpenSearchService APIs."""

    def __init__(self, region_name: str, account_id: str):
        super().__init__(region_name, account_id)
        self.domains: dict[str, OpenSearchDomain] = {}
        self.tagger = TaggingService()

    def create_domain(
        self,
        domain_name: str,
        engine_version: str,
        cluster_config: dict[str, Any],
        ebs_options: dict[str, Any],
        access_policies: str,
        snapshot_options: dict[str, Any],
        vpc_options: dict[str, Any],
        cognito_options: dict[str, Any],
        encryption_at_rest_options: dict[str, Any],
        node_to_node_encryption_options: dict[str, Any],
        advanced_options: dict[str, Any],
        log_publishing_options: dict[str, Any],
        domain_endpoint_options: dict[str, Any],
        advanced_security_options: dict[str, Any],
        tag_list: list[dict[str, str]],
        auto_tune_options: dict[str, Any],
        off_peak_window_options: dict[str, Any],
        software_update_options: dict[str, Any],
        is_es: bool,
        elasticsearch_version: Optional[str],
        elasticsearch_cluster_config: Optional[str],
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
            is_es=is_es,
            elasticsearch_version=elasticsearch_version,
            elasticsearch_cluster_config=elasticsearch_cluster_config,
        )
        self.domains[domain_name] = domain
        if tag_list:
            self.add_tags(domain.arn, tag_list)
        return domain

    def get_compatible_versions(
        self, domain_name: Optional[str]
    ) -> list[dict[str, Any]]:
        if domain_name and domain_name not in self.domains:
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

    def describe_domain_config(self, domain_name: str) -> dict[str, Any]:
        domain = self.describe_domain(domain_name)
        return domain.to_config_dict()

    def update_domain_config(
        self,
        domain_name: str,
        cluster_config: dict[str, Any],
        ebs_options: dict[str, Any],
        access_policies: str,
        snapshot_options: dict[str, Any],
        vpc_options: dict[str, Any],
        cognito_options: dict[str, Any],
        encryption_at_rest_options: dict[str, Any],
        node_to_node_encryption_options: dict[str, Any],
        advanced_options: dict[str, Any],
        log_publishing_options: dict[str, Any],
        domain_endpoint_options: dict[str, Any],
        advanced_security_options: dict[str, Any],
        auto_tune_options: dict[str, Any],
        off_peak_window_options: dict[str, Any],
        software_update_options: dict[str, Any],
    ) -> "OpenSearchDomain":
        domain = self.domains[domain_name]
        domain.cluster_config = cluster_config or domain.cluster_config
        domain.ebs_options = ebs_options or domain.ebs_options
        domain.access_policies = access_policies or domain.access_policies
        domain.snapshot_options = snapshot_options or domain.snapshot_options
        domain.vpc_options = vpc_options or domain.vpc_options
        domain.cognito_options = cognito_options or domain.cognito_options
        domain.encryption_at_rest_options = (
            encryption_at_rest_options or domain.encryption_at_rest_options
        )
        domain.node_to_node_encryption_options = (
            node_to_node_encryption_options or domain.node_to_node_encryption_options
        )
        domain.advanced_options = advanced_options or domain.advanced_options
        domain.log_publishing_options = (
            log_publishing_options or domain.log_publishing_options
        )
        domain.domain_endpoint_options = (
            domain_endpoint_options or domain.domain_endpoint_options
        )
        domain.advanced_security_options = (
            advanced_security_options or domain.advanced_security_options
        )
        domain.auto_tune_options = auto_tune_options or domain.auto_tune_options
        # Fix attribute name (singular)
        domain.off_peak_window_options = (
            off_peak_window_options or domain.off_peak_window_options
        )
        domain.software_update_options = (
            software_update_options or domain.software_update_options
        )
        return domain

    def add_tags(self, arn: str, tags: list[dict[str, str]]) -> None:
        self.tagger.tag_resource(arn, tags)

    def list_tags(self, arn: str) -> list[dict[str, str]]:
        return self.tagger.list_tags_for_resource(arn)["Tags"]

    def remove_tags(self, arn: str, tag_keys: list[str]) -> None:
        self.tagger.untag_resource_using_names(arn, tag_keys)

    def list_domain_names(self, engine_type: str) -> list[dict[str, str]]:
        domains = []
        if engine_type and engine_type not in ["Elasticsearch", "OpenSearch"]:
            raise EngineTypeNotFoundException(engine_type)
        for domain in self.domains.values():
            if engine_type:
                if engine_type == domain.engine_type:
                    domains.append(
                        {"DomainName": domain.domain_name, "EngineType": engine_type}
                    )
            else:
                domains.append(
                    {"DomainName": domain.domain_name, "EngineType": domain.engine_type}
                )
        return domains

    def describe_domains(self, domain_names: list[str]) -> list[OpenSearchDomain]:
        queried_domains = []
        for domain_name in domain_names:
            if domain_name in self.domains:
                queried_domains.append(self.domains[domain_name])
        return queried_domains


opensearch_backends = BackendDict(OpenSearchServiceBackend, "opensearch")
