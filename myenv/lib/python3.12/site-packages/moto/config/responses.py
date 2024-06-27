import json

from moto.core.responses import BaseResponse

from .models import ConfigBackend, config_backends


class ConfigResponse(BaseResponse):
    def __init__(self) -> None:
        super().__init__(service_name="config")

    @property
    def config_backend(self) -> ConfigBackend:
        return config_backends[self.current_account][self.region]

    def put_configuration_recorder(self) -> str:
        self.config_backend.put_configuration_recorder(
            self._get_param("ConfigurationRecorder")
        )
        return ""

    def put_configuration_aggregator(self) -> str:
        aggregator = self.config_backend.put_configuration_aggregator(
            json.loads(self.body)
        )
        schema = {"ConfigurationAggregator": aggregator}
        return json.dumps(schema)

    def describe_configuration_aggregators(self) -> str:
        aggregators = self.config_backend.describe_configuration_aggregators(
            self._get_param("ConfigurationAggregatorNames"),
            self._get_param("NextToken"),
            self._get_param("Limit"),
        )
        return json.dumps(aggregators)

    def delete_configuration_aggregator(self) -> str:
        self.config_backend.delete_configuration_aggregator(
            self._get_param("ConfigurationAggregatorName")
        )
        return ""

    def put_aggregation_authorization(self) -> str:
        agg_auth = self.config_backend.put_aggregation_authorization(
            self._get_param("AuthorizedAccountId"),
            self._get_param("AuthorizedAwsRegion"),
            self._get_param("Tags"),
        )
        schema = {"AggregationAuthorization": agg_auth}
        return json.dumps(schema)

    def describe_aggregation_authorizations(self) -> str:
        authorizations = self.config_backend.describe_aggregation_authorizations(
            self._get_param("NextToken"), self._get_param("Limit")
        )

        return json.dumps(authorizations)

    def delete_aggregation_authorization(self) -> str:
        self.config_backend.delete_aggregation_authorization(
            self._get_param("AuthorizedAccountId"),
            self._get_param("AuthorizedAwsRegion"),
        )

        return ""

    def describe_configuration_recorders(self) -> str:
        recorders = self.config_backend.describe_configuration_recorders(
            self._get_param("ConfigurationRecorderNames")
        )
        schema = {"ConfigurationRecorders": recorders}
        return json.dumps(schema)

    def describe_configuration_recorder_status(self) -> str:
        recorder_statuses = self.config_backend.describe_configuration_recorder_status(
            self._get_param("ConfigurationRecorderNames")
        )
        schema = {"ConfigurationRecordersStatus": recorder_statuses}
        return json.dumps(schema)

    def put_delivery_channel(self) -> str:
        self.config_backend.put_delivery_channel(self._get_param("DeliveryChannel"))
        return ""

    def describe_delivery_channels(self) -> str:
        delivery_channels = self.config_backend.describe_delivery_channels(
            self._get_param("DeliveryChannelNames")
        )
        schema = {"DeliveryChannels": delivery_channels}
        return json.dumps(schema)

    def describe_delivery_channel_status(self) -> str:
        raise NotImplementedError()

    def delete_delivery_channel(self) -> str:
        self.config_backend.delete_delivery_channel(
            self._get_param("DeliveryChannelName")
        )
        return ""

    def delete_configuration_recorder(self) -> str:
        self.config_backend.delete_configuration_recorder(
            self._get_param("ConfigurationRecorderName")
        )
        return ""

    def start_configuration_recorder(self) -> str:
        self.config_backend.start_configuration_recorder(
            self._get_param("ConfigurationRecorderName")
        )
        return ""

    def stop_configuration_recorder(self) -> str:
        self.config_backend.stop_configuration_recorder(
            self._get_param("ConfigurationRecorderName")
        )
        return ""

    def list_discovered_resources(self) -> str:
        schema = self.config_backend.list_discovered_resources(
            resource_type=self._get_param("resourceType"),
            resource_ids=self._get_param("resourceIds"),
            resource_name=self._get_param("resourceName"),
            limit=self._get_param("limit"),
            next_token=self._get_param("nextToken"),
        )
        return json.dumps(schema)

    def list_aggregate_discovered_resources(self) -> str:
        schema = self.config_backend.list_aggregate_discovered_resources(
            self._get_param("ConfigurationAggregatorName"),
            self._get_param("ResourceType"),
            self._get_param("Filters"),
            self._get_param("Limit"),
            self._get_param("NextToken"),
        )
        return json.dumps(schema)

    def list_tags_for_resource(self) -> str:
        schema = self.config_backend.list_tags_for_resource(
            self._get_param("ResourceArn"), self._get_param("Limit")
        )
        return json.dumps(schema)

    def get_resource_config_history(self) -> str:
        schema = self.config_backend.get_resource_config_history(
            self._get_param("resourceType"), self._get_param("resourceId"), self.region
        )
        return json.dumps(schema)

    def batch_get_resource_config(self) -> str:
        schema = self.config_backend.batch_get_resource_config(
            self._get_param("resourceKeys"), self.region
        )
        return json.dumps(schema)

    def batch_get_aggregate_resource_config(self) -> str:
        schema = self.config_backend.batch_get_aggregate_resource_config(
            self._get_param("ConfigurationAggregatorName"),
            self._get_param("ResourceIdentifiers"),
        )
        return json.dumps(schema)

    def put_evaluations(self) -> str:
        evaluations = self.config_backend.put_evaluations(
            self._get_param("Evaluations"),
            self._get_param("ResultToken"),
            self._get_param("TestMode"),
        )
        return json.dumps(evaluations)

    def put_organization_conformance_pack(self) -> str:
        conformance_pack = self.config_backend.put_organization_conformance_pack(
            name=self._get_param("OrganizationConformancePackName"),
            template_s3_uri=self._get_param("TemplateS3Uri"),
            template_body=self._get_param("TemplateBody"),
            delivery_s3_bucket=self._get_param("DeliveryS3Bucket"),
            delivery_s3_key_prefix=self._get_param("DeliveryS3KeyPrefix"),
            input_parameters=self._get_param("ConformancePackInputParameters"),
            excluded_accounts=self._get_param("ExcludedAccounts"),
        )

        return json.dumps(conformance_pack)

    def describe_organization_conformance_packs(self) -> str:
        conformance_packs = self.config_backend.describe_organization_conformance_packs(
            self._get_param("OrganizationConformancePackNames")
        )
        return json.dumps(conformance_packs)

    def describe_organization_conformance_pack_statuses(self) -> str:
        statuses = self.config_backend.describe_organization_conformance_pack_statuses(
            self._get_param("OrganizationConformancePackNames")
        )
        return json.dumps(statuses)

    def get_organization_conformance_pack_detailed_status(self) -> str:
        # 'Filters' parameter is not implemented yet
        statuses = (
            self.config_backend.get_organization_conformance_pack_detailed_status(
                self._get_param("OrganizationConformancePackName")
            )
        )
        return json.dumps(statuses)

    def delete_organization_conformance_pack(self) -> str:
        self.config_backend.delete_organization_conformance_pack(
            self._get_param("OrganizationConformancePackName")
        )
        return ""

    def tag_resource(self) -> str:
        self.config_backend.tag_resource(
            self._get_param("ResourceArn"), self._get_param("Tags")
        )
        return ""

    def untag_resource(self) -> str:
        self.config_backend.untag_resource(
            self._get_param("ResourceArn"), self._get_param("TagKeys")
        )
        return ""

    def put_config_rule(self) -> str:
        self.config_backend.put_config_rule(
            self._get_param("ConfigRule"), self._get_param("Tags")
        )
        return ""

    def describe_config_rules(self) -> str:
        rules = self.config_backend.describe_config_rules(
            self._get_param("ConfigRuleNames"), self._get_param("NextToken")
        )
        return json.dumps(rules)

    def delete_config_rule(self) -> str:
        self.config_backend.delete_config_rule(self._get_param("ConfigRuleName"))
        return ""

    def put_retention_configuration(self) -> str:
        retention_configuration = self.config_backend.put_retention_configuration(
            self._get_param("RetentionPeriodInDays")
        )
        return json.dumps(retention_configuration)

    def describe_retention_configurations(self) -> str:
        retention_configurations = (
            self.config_backend.describe_retention_configurations(
                self._get_param("RetentionConfigurationNames")
            )
        )
        schema = {"RetentionConfigurations": retention_configurations}
        return json.dumps(schema)

    def delete_retention_configuration(self) -> str:
        self.config_backend.delete_retention_configuration(
            self._get_param("RetentionConfigurationName")
        )
        return ""
