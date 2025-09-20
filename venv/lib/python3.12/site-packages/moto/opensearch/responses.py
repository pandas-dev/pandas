"""Handles incoming opensearch requests, invokes methods, returns responses."""

import json
import re
from typing import Any

from moto.core.common_types import TYPE_RESPONSE
from moto.core.responses import BaseResponse
from moto.es.exceptions import InvalidDomainName

from .models import OpenSearchServiceBackend, opensearch_backends


class OpenSearchServiceResponse(BaseResponse):
    """Handler for OpenSearchService requests and responses."""

    def __init__(self) -> None:
        super().__init__(service_name="opensearch")

    @property
    def opensearch_backend(self) -> OpenSearchServiceBackend:
        return opensearch_backends[self.current_account][self.region]

    @classmethod
    def list_domains(cls, request: Any, full_url: str, headers: Any) -> TYPE_RESPONSE:  # type: ignore
        response = cls()
        response.setup_class(request, full_url, headers)
        if request.method == "GET":
            return 200, {}, response.list_domain_names()
        if request.method == "POST":
            return 200, {}, response.describe_domains()

    @classmethod
    def domains(cls, request: Any, full_url: str, headers: Any) -> TYPE_RESPONSE:  # type: ignore
        response = cls()
        response.setup_class(request, full_url, headers)
        if request.method == "POST":
            return 200, {}, response.create_domain()

    @classmethod
    def domain(cls, request: Any, full_url: str, headers: Any) -> TYPE_RESPONSE:  # type: ignore
        response = cls()
        response.setup_class(request, full_url, headers)
        if request.method == "DELETE":
            return 200, {}, response.delete_domain()
        if request.method == "GET":
            return 200, {}, response.describe_domain()

    @classmethod
    def tags(cls, request: Any, full_url: str, headers: Any) -> TYPE_RESPONSE:  # type: ignore
        response = cls()
        response.setup_class(request, full_url, headers)
        if request.method == "GET":
            return 200, {}, response.list_tags()
        if request.method == "POST":
            return 200, {}, response.add_tags()

    @classmethod
    def tag_removal(cls, request: Any, full_url: str, headers: Any) -> TYPE_RESPONSE:  # type: ignore
        response = cls()
        response.setup_class(request, full_url, headers)
        if request.method == "POST":
            return 200, {}, response.remove_tags()

    def create_domain(self) -> str:
        domain_name = self._get_param("DomainName")
        if not re.match(r"^[a-z][a-z0-9\-]+$", domain_name):
            raise InvalidDomainName(domain_name)
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
        # ElasticSearch specific options
        is_es = self.parsed_url.path.endswith("/es/domain")
        elasticsearch_version = self._get_param("ElasticsearchVersion")
        elasticsearch_cluster_config = self._get_param("ElasticsearchClusterConfig")
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
            is_es=is_es,
            elasticsearch_version=elasticsearch_version,
            elasticsearch_cluster_config=elasticsearch_cluster_config,
        )
        return json.dumps(dict(DomainStatus=domain.to_dict()))

    def get_compatible_versions(self) -> str:
        domain_name = self._get_param("domainName")
        compatible_versions = self.opensearch_backend.get_compatible_versions(
            domain_name=domain_name,
        )
        return json.dumps(dict(CompatibleVersions=compatible_versions))

    def delete_domain(self) -> str:
        domain_name = self.path.split("/")[-1]
        domain = self.opensearch_backend.delete_domain(
            domain_name=domain_name,
        )
        return json.dumps(dict(DomainStatus=domain.to_dict()))

    def describe_domain(self) -> str:
        domain_name = self.path.split("/")[-1]
        if not re.match(r"^[a-z][a-z0-9\-]+$", domain_name):
            raise InvalidDomainName(domain_name)
        domain = self.opensearch_backend.describe_domain(
            domain_name=domain_name,
        )
        return json.dumps(dict(DomainStatus=domain.to_dict()))

    def describe_domain_config(self) -> str:
        # Supports both body param and URL form (/domain/{name}/config)
        domain_name = self._get_param("DomainName")
        if not domain_name and self.path:
            parts = [p for p in self.path.split("/") if p]
            if len(parts) >= 2 and parts[-1] == "config":
                domain_name = parts[-2]
        config = self.opensearch_backend.describe_domain_config(domain_name=domain_name)  # type: ignore[arg-type]
        return json.dumps({"DomainConfig": config})

    @classmethod
    def describe_es_domain_config(
        cls, request: Any, full_url: str, headers: Any
    ) -> TYPE_RESPONSE:
        response = cls()
        response.setup_class(request, full_url, headers)

        domain_name = request.url.split("/")[-2]
        domain_config = response.opensearch_backend.describe_domain_config(
            domain_name=domain_name,
        )

        return 200, {}, json.dumps(dict(DomainConfig=domain_config))

    def update_domain_config(self) -> str:
        domain_name = self._get_param("DomainName")
        domain = self.opensearch_backend.update_domain_config(
            domain_name=domain_name,
            cluster_config=self._get_param("ClusterConfig"),
            ebs_options=self._get_param("EBSOptions"),
            access_policies=self._get_param("AccessPolicies"),
            snapshot_options=self._get_param("SnapshotOptions"),
            vpc_options=self._get_param("VPCOptions"),
            cognito_options=self._get_param("CognitoOptions"),
            encryption_at_rest_options=self._get_param("EncryptionAtRestOptions"),
            node_to_node_encryption_options=self._get_param(
                "NodeToNodeEncryptionOptions"
            ),
            advanced_options=self._get_param("AdvancedOptions"),
            log_publishing_options=self._get_param("LogPublishingOptions"),
            domain_endpoint_options=self._get_param("DomainEndpointOptions"),
            advanced_security_options=self._get_param("AdvancedSecurityOptions"),
            auto_tune_options=self._get_param("AutoTuneOptions"),
            off_peak_window_options=self._get_param("OffPeakWindowOptions"),
            software_update_options=self._get_param("SoftwareUpdateOptions"),
        )
        return json.dumps({"DomainConfig": domain.to_config_dict()})

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

    def describe_domains(self) -> str:
        domain_names = self._get_param("DomainNames")
        domains = self.opensearch_backend.describe_domains(
            domain_names=domain_names,
        )
        domain_list = [domain.to_dict() for domain in domains]
        return json.dumps({"DomainStatusList": domain_list})
