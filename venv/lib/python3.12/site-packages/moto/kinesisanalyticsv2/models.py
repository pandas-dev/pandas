from __future__ import annotations

import random
from datetime import datetime
from typing import Any, Dict, List, Optional

from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import BaseModel
from moto.moto_api._internal import mock_random
from moto.utilities.tagging_service import TaggingService

FAKE_VPC_ID = "vpc-0123456789abcdef0"


class Application(BaseModel):
    def __init__(
        self,
        account_id: str,
        region_name: str,
        application_name: str,
        application_description: Optional[str],
        runtime_environment: str,
        service_execution_role: str,
        application_configuration: Optional[Dict[str, Any]],
        cloud_watch_logging_options: Optional[List[Dict[str, str]]],
        application_mode: Optional[str],
    ):
        self.account_id = account_id
        self.region_name = region_name
        self.application_name = application_name
        self.application_description = application_description
        self.runtime_environment = runtime_environment
        self.service_execution_role = service_execution_role
        self.application_mode = application_mode

        self.app_config_description = (
            self._generate_app_config_description(application_configuration)
            if application_configuration
            else None
        )
        self.cloud_watch_logging_description = self._generate_logging_options(
            cloud_watch_logging_options
        )

        self.application_arn = self._generate_arn()
        self.application_status = "STARTING"
        self.application_version_id = 1
        self.creation_date_time = datetime.now().isoformat()
        self.last_updated_date_time = datetime.now().isoformat()
        self.conditional_token = str(mock_random.uuid4()).replace("-", "")

    def _generate_arn(self) -> str:
        return f"arn:aws:kinesisanalytics:{self.region_name}:{self.account_id}:application/{self.application_name}"

    def _generate_logging_options(
        self, cloud_watch_logging_options: Optional[List[Dict[str, str]]]
    ) -> List[Dict[str, str]] | None:
        cloud_watch_logging_option_descriptions = []
        option_id = f"{str(random.randint(1, 100))}.1"

        # Leaving out RoleARN since it is provided only sometimes for backwards
        # compatibility. Current API versions do not have the resource-level
        # role.
        if cloud_watch_logging_options:
            for i in cloud_watch_logging_options:
                cloud_watch_logging_option_descriptions.append(
                    {
                        "CloudWatchLoggingOptionId": option_id,
                        "LogStreamARN": i["LogStreamARN"],
                    }
                )
            return cloud_watch_logging_option_descriptions
        else:
            return None

    # The app_config description does not include
    # - "SqlApplicationConfigurationDescription" (discontinued)
    # - "RunConfigurationDescription" (which requires start_application)
    def _generate_app_config_description(
        self, app_config: Dict[str, Any]
    ) -> Dict[str, Any]:
        # Keys that do not have extra values in the description besides renamed keys
        UPDATABLE_APP_CONFIG_TOP_LEVEL_KEYS = {
            "EnvironmentProperties": "EnvironmentPropertyDescriptions",
            "ApplicationSnapshotConfiguration": "ApplicationSnapshotConfigurationDescription",
            "ApplicationSystemRollbackConfiguration": "ApplicationSystemRollbackConfigurationDescription",
            "ZeppelinApplicationConfiguration": "ZeppelinApplicationConfigurationDescription",
        }

        APP_CONFIG_SUBFIELD_KEYS = {
            "PropertyGroups": "PropertyGroupDescriptions",
            "MonitoringConfiguration": "MonitoringConfigurationDescription",
            "CatalogConfiguration": "CatalogConfigurationDescription",
            "DeployAsApplicationConfiguration": "DeployAsApplicationConfigurationDescription",
            "S3ContentLocation": "S3ContentLocationDescription",
            "CustomArtifactsConfiguration": "CustomArtifactsConfigurationDescription",
            "GlueDataCatalogConfiguration": "GlueDataCatalogConfigurationDescription",
            "MavenReference": "MavenReferenceDescription",
        }

        app_config_description = {}
        if app_config:
            if "FlinkApplicationConfiguration" in app_config:
                app_config_description["FlinkApplicationConfigurationDescription"] = (
                    self.__generate_flink_app_description(app_config)
                )

            for old_key, new_key in UPDATABLE_APP_CONFIG_TOP_LEVEL_KEYS.items():
                if old_key in app_config:
                    app_config_description[new_key] = self.__update_keys(
                        app_config[old_key], APP_CONFIG_SUBFIELD_KEYS
                    )

            app_code_config = app_config.get("ApplicationCodeConfiguration")
            if app_code_config:
                new_key = "ApplicationCodeConfigurationDescription"

                # S3ContentLocation has a different value, so keeping it
                # separate from APP_CONFIG_SUBFIELD_KEYS
                app_code_config_keys = {
                    "S3ContentLocation": "S3ApplicationCodeLocationDescription",
                    "CodeContent": "CodeContentDescription",
                }

                app_config_description[new_key] = self.__update_keys(
                    app_code_config, app_code_config_keys
                )

                if app_code_config["CodeContentType"] == "ZIPFILE":
                    app_config_description[new_key]["CodeContentDescription"][
                        "CodeMD5"
                    ] = "fakechecksum"
                    app_config_description[new_key]["CodeContentDescription"][
                        "CodeSize"
                    ] = 123

            if "VpcConfigurations" in app_config:
                app_config_description["VpcConfigurationDescriptions"] = app_config[
                    "VpcConfigurations"
                ]
                for index, vpc_config in enumerate(
                    app_config_description["VpcConfigurationDescriptions"]
                ):
                    vpc_config["VpcConfigurationId"] = str(index + 1.1)  # type: ignore[index]
                    # FAKE_VPC_ID hardcoded, not a value from the parameters
                    vpc_config["VpcId"] = FAKE_VPC_ID  # type: ignore[index]

        return app_config_description

    def __generate_flink_app_description(
        self, app_config: Dict[str, Any]
    ) -> Dict[str, Any]:
        flink_config_description = {}
        flink_config = app_config.get("FlinkApplicationConfiguration")
        if flink_config and isinstance(flink_config, dict):
            checkpoint_config = flink_config.get("CheckpointConfiguration")
            if checkpoint_config and isinstance(checkpoint_config, dict):
                if checkpoint_config.get("ConfigurationType") == "DEFAULT":
                    flink_config_description["CheckpointConfigurationDescription"] = {
                        "ConfigurationType": "DEFAULT",
                        "CheckpointingEnabled": True,
                        "CheckpointInterval": 60000,
                        "MinPauseBetweenCheckpoints": 5000,
                    }
                elif checkpoint_config.get("ConfigurationType") == "CUSTOM":
                    flink_config_description["CheckpointConfigurationDescription"] = {
                        "ConfigurationType": "CUSTOM",
                        "CheckpointingEnabled": checkpoint_config.get(
                            "CheckpointingEnabled", True
                        ),
                        "CheckpointInterval": checkpoint_config.get(
                            "CheckpointInterval", 60000
                        ),
                        "MinPauseBetweenCheckpoints": checkpoint_config.get(
                            "MinPauseBetweenCheckpoints", 5000
                        ),
                    }

            monitoring_config = flink_config.get("MonitoringConfiguration")
            if monitoring_config and isinstance(monitoring_config, dict):
                if monitoring_config.get("ConfigurationType") == "DEFAULT":
                    flink_config_description["MonitoringConfigurationDescription"] = {
                        "ConfigurationType": "DEFAULT",
                        "MetricsLevel": "APPLICATION",
                        "LogLevel": "INFO",
                    }
                elif monitoring_config.get("ConfigurationType") == "CUSTOM":
                    flink_config_description["MonitoringConfigurationDescription"] = {
                        "ConfigurationType": "CUSTOM",
                        "MetricsLevel": monitoring_config.get(
                            "MetricsLevel", "APPLICATION"
                        ),
                        "LogLevel": monitoring_config.get("LogLevel", "INFO"),
                    }
            parallel_config = flink_config.get("ParallelismConfiguration")
            if monitoring_config and isinstance(parallel_config, dict):
                if parallel_config.get("ConfigurationType") == "DEFAULT":
                    flink_config_description["ParallelismConfigurationDescription"] = {
                        "ConfigurationType": "DEFAULT",
                        "Parallelism": 1,
                        "ParallelismPerKPU": 1,
                        "AutoScalingEnabled": False,
                        "CurrentParallelism": 1,
                    }
                elif parallel_config.get("ConfigurationType") == "CUSTOM":
                    flink_config_description["ParallelismConfigurationDescription"] = {
                        "ConfigurationType": "CUSTOM",
                        "Parallelism": parallel_config.get("Parallelism", 1),
                        "ParallelismPerKPU": parallel_config.get(
                            "ParallelismPerKPU", 1
                        ),
                        "AutoScalingEnabled": parallel_config.get(
                            "AutoScalingEnabled", False
                        ),
                        "CurrentParallelism": parallel_config.get("Parallelism", 1),
                    }
        return flink_config_description

    def __update_keys(self, old_dict: Any, key_map: Dict[str, str]) -> Any:
        if not isinstance(old_dict, dict):
            return old_dict

        updated_dict = {}
        for old_key, value in old_dict.items():
            # Check if the current key is in key_map, else keep old_key
            new_key = key_map.get(old_key, old_key)

            if isinstance(value, dict):
                updated_dict[new_key] = self.__update_keys(value, key_map)
            elif isinstance(value, list):
                updated_dict[new_key] = [
                    self.__update_keys(list_item, key_map) for list_item in value
                ]
            else:
                updated_dict[new_key] = value
        return updated_dict


class KinesisAnalyticsV2Backend(BaseBackend):
    """Implementation of KinesisAnalyticsV2 APIs."""

    def __init__(self, region_name: str, account_id: str) -> None:
        super().__init__(region_name, account_id)
        self.applications: Dict[str, Application] = {}
        self.tagger = TaggingService(
            tag_name="Tags", key_name="Key", value_name="Value"
        )

    def create_application(
        self,
        application_name: str,
        application_description: Optional[str],
        runtime_environment: str,
        service_execution_role: str,
        application_configuration: Optional[Dict[str, Any]],
        cloud_watch_logging_options: Optional[List[Dict[str, str]]],
        tags: Optional[List[Dict[str, str]]],
        application_mode: Optional[str],
    ) -> Dict[str, Any]:
        app = Application(
            account_id=self.account_id,
            region_name=self.region_name,
            application_name=application_name,
            application_description=application_description,
            runtime_environment=runtime_environment,
            service_execution_role=service_execution_role,
            application_configuration=application_configuration,
            cloud_watch_logging_options=cloud_watch_logging_options,
            application_mode=application_mode,
        )

        self.applications[application_name] = app

        if tags:
            self.tag_resource(resource_arn=app.application_arn, tags=tags)
        return {
            "ApplicationARN": app.application_arn,
            "ApplicationDescription": app.application_description,
            "RuntimeEnvironment": app.runtime_environment,
            "ServiceExecutionRole": app.service_execution_role,
            "ApplicationStatus": app.application_status,
            "ApplicationVersionId": app.application_version_id,
            "CreateTimestamp": app.creation_date_time,
            "LastUpdateTimestamp": app.last_updated_date_time,
            "ApplicationConfigurationDescription": app.app_config_description,
            "CloudWatchLoggingOptionDescriptions": app.cloud_watch_logging_description,
            "ApplicationMaintenanceConfigurationDescription": {
                "ApplicationMaintenanceWindowStartTime": "06:00",
                "ApplicationMaintenanceWindowEndTime": "14:00",
            },
            "ApplicationVersionCreateTimestamp": str(app.creation_date_time),
            "ConditionalToken": app.conditional_token,
            "ApplicationMode": app.application_mode,
        }

    def tag_resource(self, resource_arn: str, tags: List[Dict[str, str]]) -> None:
        self.tagger.tag_resource(resource_arn, tags)

    def list_tags_for_resource(self, resource_arn: str) -> List[Dict[str, str]]:
        return self.tagger.list_tags_for_resource(resource_arn)["Tags"]

    def describe_application(
        self,
        application_name: str,
    ) -> Dict[str, Any]:
        app = self.applications[application_name]
        return {
            "ApplicationARN": app.application_arn,
            "ApplicationDescription": app.application_description,
            "RuntimeEnvironment": app.runtime_environment,
            "ServiceExecutionRole": app.service_execution_role,
            "ApplicationStatus": app.application_status,
            "ApplicationVersionId": app.application_version_id,
            "CreateTimestamp": app.creation_date_time,
            "LastUpdateTimestamp": app.last_updated_date_time,
            "ApplicationConfigurationDescription": app.app_config_description,
            "CloudWatchLoggingOptionDescriptions": app.cloud_watch_logging_description,
            "ApplicationMaintenanceConfigurationDescription": {
                "ApplicationMaintenanceWindowStartTime": "06:00",
                "ApplicationMaintenanceWindowEndTime": "14:00",
            },
            "ApplicationVersionCreateTimestamp": str(app.creation_date_time),
            "ConditionalToken": app.conditional_token,
            "ApplicationMode": app.application_mode,
        }

    def list_applications(self) -> List[Dict[str, Any]]:
        application_summaries = [
            {
                "ApplicationName": app.application_name,
                "ApplicationARN": app.application_arn,
                "ApplicationStatus": app.application_status,
                "ApplicationVersionId": app.application_version_id,
                "RuntimeEnvironment": app.runtime_environment,
                "ApplicationMode": app.application_mode,
            }
            for app in self.applications.values()
        ]
        return application_summaries


kinesisanalyticsv2_backends = BackendDict(
    KinesisAnalyticsV2Backend, "kinesisanalyticsv2"
)
