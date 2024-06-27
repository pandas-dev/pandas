"""Handles incoming opensearch requests, invokes methods, returns responses."""

import json

from moto.core.responses import BaseResponse

from .models import OpenSearchServiceBackend, opensearch_backends


class OpenSearchServiceResponse(BaseResponse):
    """Handler for OpenSearchService requests and responses."""

    def __init__(self) -> None:
        super().__init__(service_name="opensearch")

    @property
    def opensearch_backend(self) -> OpenSearchServiceBackend:
        """Return backend instance specific for this region."""
        return opensearch_backends[self.current_account][self.region]

    def create_domain(self) -> str:
        domain_name = self._get_param("DomainName")
        engine_version = self._get_param("EngineVersion")
        cluster_config = self._get_param("ClusterConfig")
        ebs_options = self._get_param("EBSOptions")
        access_policies = self._get_param("AccessPolicies")
        snapshot_options = self._get_param("SnapshotOptions")
        vpc_options = self._get_param("VPCOptions")
        cognito_options = self._get_param("CognitoOptions")
        encryption_at_rest_options = self._get_param("EncryptionAtRestOptions")
        node_to_node_encryption_options = self._get_param("NodeToNodeEncryptionOptions")
        advanced_options = self._get_param("AdvancedOptions")
        log_publishing_options = self._get_param("LogPublishingOptions")
        domain_endpoint_options = self._get_param("DomainEndpointOptions")
        advanced_security_options = self._get_param("AdvancedSecurityOptions")
        tag_list = self._get_param("TagList")
        auto_tune_options = self._get_param("AutoTuneOptions")
        off_peak_window_options = self._get_param("OffPeakWindowOptions")
        software_update_options = self._get_param("SoftwareUpdateOptions")
        domain = self.opensearch_backend.create_domain(
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
            tag_list=tag_list,
            auto_tune_options=auto_tune_options,
            off_peak_window_options=off_peak_window_options,
            software_update_options=software_update_options,
        )
        return json.dumps(dict(DomainStatus=domain.to_dict()))

    def get_compatible_versions(self) -> str:
        domain_name = self._get_param("domainName")
        compatible_versions = self.opensearch_backend.get_compatible_versions(
            domain_name=domain_name,
        )
        return json.dumps(dict(CompatibleVersions=compatible_versions))

    def delete_domain(self) -> str:
        domain_name = self._get_param("DomainName")
        domain = self.opensearch_backend.delete_domain(
            domain_name=domain_name,
        )
        return json.dumps(dict(DomainStatus=domain.to_dict()))

    def describe_domain(self) -> str:
        domain_name = self._get_param("DomainName")
        domain = self.opensearch_backend.describe_domain(
            domain_name=domain_name,
        )
        return json.dumps(dict(DomainStatus=domain.to_dict()))

    def describe_domain_config(self) -> str:
        domain_name = self._get_param("DomainName")
        domain = self.opensearch_backend.describe_domain_config(
            domain_name=domain_name,
        )
        return json.dumps(dict(DomainConfig=domain.to_config_dict()))

    def update_domain_config(self) -> str:
        domain_name = self._get_param("DomainName")
        cluster_config = self._get_param("ClusterConfig")
        ebs_options = self._get_param("EBSOptions")
        access_policies = self._get_param("AccessPolicies")
        snapshot_options = self._get_param("SnapshotOptions")
        vpc_options = self._get_param("VPCOptions")
        cognito_options = self._get_param("CognitoOptions")
        encryption_at_rest_options = self._get_param("EncryptionAtRestOptions")
        node_to_node_encryption_options = self._get_param("NodeToNodeEncryptionOptions")
        advanced_options = self._get_param("AdvancedOptions")
        log_publishing_options = self._get_param("LogPublishingOptions")
        domain_endpoint_options = self._get_param("DomainEndpointOptions")
        advanced_security_options = self._get_param("AdvancedSecurityOptions")
        auto_tune_options = self._get_param("AutoTuneOptions")
        off_peak_window_options = self._get_param("OffPeakWindowOptions")
        software_update_options = self._get_param("SoftwareUpdateOptions")
        domain = self.opensearch_backend.update_domain_config(
            domain_name=domain_name,
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
        return json.dumps(dict(DomainConfig=domain.to_config_dict()))

    def list_tags(self) -> str:
        arn = self._get_param("arn")
        tags = self.opensearch_backend.list_tags(arn)
        return json.dumps({"TagList": tags})

    def add_tags(self) -> str:
        arn = self._get_param("ARN")
        tags = self._get_param("TagList")
        self.opensearch_backend.add_tags(arn, tags)
        return "{}"

    def remove_tags(self) -> str:
        arn = self._get_param("ARN")
        tag_keys = self._get_param("TagKeys")
        self.opensearch_backend.remove_tags(arn, tag_keys)
        return "{}"

    def list_domain_names(self) -> str:
        engine_type = self._get_param("engineType")
        domain_names = self.opensearch_backend.list_domain_names(
            engine_type=engine_type,
        )
        return json.dumps(dict(DomainNames=domain_names))
