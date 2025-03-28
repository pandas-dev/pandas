import json
import os
import random
import re
import string
from collections import defaultdict
from datetime import datetime
from typing import Any, DefaultDict, Dict, Iterable, List, Optional, Union, cast

from dateutil.tz import tzutc

from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import BaseModel, CloudFormationModel
from moto.core.utils import camelcase_to_underscores
from moto.sagemaker import validators
from moto.utilities.paginator import paginate
from moto.utilities.utils import ARN_PARTITION_REGEX, get_partition

from .exceptions import (
    AWSValidationException,
    ConflictException,
    MissingModel,
    ResourceInUseException,
    ResourceNotFound,
    ValidationError,
)
from .utils import (
    arn_formatter,
    filter_model_cards,
    get_pipeline_execution_from_arn,
    get_pipeline_from_name,
    get_pipeline_name_from_execution_arn,
    load_pipeline_definition_from_s3,
    validate_model_approval_status,
)

PAGINATION_MODEL = {
    "list_experiments": {
        "input_token": "NextToken",
        "limit_key": "MaxResults",
        "limit_default": 100,
        "unique_attribute": "arn",
    },
    "list_trials": {
        "input_token": "NextToken",
        "limit_key": "MaxResults",
        "limit_default": 100,
        "unique_attribute": "arn",
    },
    "list_trial_components": {
        "input_token": "NextToken",
        "limit_key": "MaxResults",
        "limit_default": 100,
        "unique_attribute": "arn",
    },
    "list_tags": {
        "input_token": "NextToken",
        "limit_key": "MaxResults",
        "limit_default": 50,
        "unique_attribute": "Key",
    },
    "list_model_package_groups": {
        "input_token": "next_token",
        "limit_key": "max_results",
        "limit_default": 100,
        "unique_attribute": "arn",
    },
    "list_model_packages": {
        "input_token": "next_token",
        "limit_key": "max_results",
        "limit_default": 100,
        "unique_attribute": "arn",
    },
    "list_notebook_instances": {
        "input_token": "next_token",
        "limit_key": "max_results",
        "limit_default": 100,
        "unique_attribute": "arn",
    },
    "list_clusters": {
        "input_token": "next_token",
        "limit_key": "max_results",
        "limit_default": 100,
        "unique_attribute": "arn",
    },
    "list_cluster_nodes": {
        "input_token": "next_token",
        "limit_key": "max_results",
        "limit_default": 100,
        "unique_attribute": "arn",
    },
    "list_auto_ml_jobs": {
        "input_token": "next_token",
        "limit_key": "max_results",
        "limit_default": 100,
        "unique_attribute": "arn",
    },
    "list_endpoints": {
        "input_token": "next_token",
        "limit_key": "max_results",
        "limit_default": 100,
        "unique_attribute": "arn",
    },
    "list_endpoint_configs": {
        "input_token": "next_token",
        "limit_key": "max_results",
        "limit_default": 100,
        "unique_attribute": "arn",
    },
    "list_compilation_jobs": {
        "input_token": "next_token",
        "limit_key": "max_results",
        "limit_default": 100,
        "unique_attribute": "arn",
    },
    "list_domains": {
        "input_token": "next_token",
        "limit_key": "max_results",
        "limit_default": 100,
        "unique_attribute": "arn",
    },
    "list_model_explainability_job_definitions": {
        "input_token": "next_token",
        "limit_key": "max_results",
        "limit_default": 100,
        "unique_attribute": "arn",
    },
    "list_hyper_parameter_tuning_jobs": {
        "input_token": "next_token",
        "limit_key": "max_results",
        "limit_default": 100,
        "unique_attribute": "arn",
    },
    "list_model_quality_job_definitions": {
        "input_token": "next_token",
        "limit_key": "max_results",
        "limit_default": 100,
        "unique_attribute": "arn",
    },
    "list_model_cards": {
        "input_token": "next_token",
        "limit_key": "max_results",
        "limit_default": 100,
        "unique_attribute": "arn",
    },
    "list_model_card_versions": {
        "input_token": "next_token",
        "limit_key": "max_results",
        "limit_default": 100,
        "unique_attribute": "model_card_arn",
    },
    "list_model_bias_job_definitions": {
        "input_token": "next_token",
        "limit_key": "max_results",
        "limit_default": 100,
        "unique_attribute": "arn",
    },
    "list_data_quality_job_definitions": {
        "input_token": "next_token",
        "limit_key": "max_results",
        "limit_default": 100,
        "unique_attribute": "arn",
    },
}

METRIC_INFO_TYPE = Dict[str, Union[str, int, float, datetime]]
METRIC_STEP_TYPE = Dict[int, METRIC_INFO_TYPE]


class BaseObject(BaseModel):
    def camelCase(self, key: str) -> str:
        words = []
        for word in key.split("_"):
            words.append(word.title())
        return "".join(words)

    def update(self, details_json: str) -> None:
        details = json.loads(details_json)
        for k in details.keys():
            setattr(self, k, details[k])

    def gen_response_object(self) -> Dict[str, Any]:
        response_object: Dict[str, Any] = dict()
        for key, value in self.__dict__.items():
            if "_" in key:
                response_object[self.camelCase(key)] = value
            else:
                response_object[key[0].upper() + key[1:]] = value
        return response_object

    @property
    def response_object(self) -> Dict[str, Any]:  # type: ignore[misc]
        return self.gen_response_object()


class FakePipelineExecution(BaseObject):
    def __init__(
        self,
        pipeline_execution_arn: str,
        pipeline_execution_display_name: str,
        pipeline_parameters: List[Dict[str, str]],
        pipeline_execution_description: str,
        parallelism_configuration: Dict[str, int],
        pipeline_definition: str,
        client_request_token: str,
    ):
        self.arn = pipeline_execution_arn
        self.pipeline_execution_display_name = pipeline_execution_display_name
        self.pipeline_parameters = pipeline_parameters
        self.pipeline_execution_description = pipeline_execution_description
        self.pipeline_execution_status = "Succeeded"
        self.pipeline_execution_failure_reason = None
        self.parallelism_configuration = parallelism_configuration
        self.pipeline_definition_for_execution = pipeline_definition
        self.client_request_token = client_request_token

        now_string = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        self.creation_time = now_string
        self.last_modified_time = now_string
        self.start_time = now_string

        fake_user_profile_name = "fake-user-profile-name"
        fake_domain_id = "fake-domain-id"
        fake_user_profile_arn = arn_formatter(
            "user-profile",
            f"{fake_domain_id}/{fake_user_profile_name}",
            pipeline_execution_arn.split(":")[4],
            pipeline_execution_arn.split(":")[3],
        )
        self.created_by = {
            "UserProfileArn": fake_user_profile_arn,
            "UserProfileName": fake_user_profile_name,
            "DomainId": fake_domain_id,
        }
        self.last_modified_by = {
            "UserProfileArn": fake_user_profile_arn,
            "UserProfileName": fake_user_profile_name,
            "DomainId": fake_domain_id,
        }


class FakePipeline(BaseObject):
    def __init__(
        self,
        pipeline_name: str,
        pipeline_display_name: str,
        pipeline_definition: str,
        pipeline_description: str,
        role_arn: str,
        tags: List[Dict[str, str]],
        account_id: str,
        region_name: str,
        parallelism_configuration: Dict[str, int],
    ):
        self.pipeline_name = pipeline_name
        self.arn = arn_formatter("pipeline", pipeline_name, account_id, region_name)
        self.pipeline_display_name = pipeline_display_name or pipeline_name
        self.pipeline_definition = pipeline_definition
        self.pipeline_description = pipeline_description
        self.pipeline_executions: Dict[str, FakePipelineExecution] = dict()
        self.role_arn = role_arn
        self.tags = tags or []
        self.parallelism_configuration = parallelism_configuration

        now_string = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        self.creation_time = now_string
        self.last_modified_time = now_string
        self.last_execution_time: Optional[str] = None

        self.pipeline_status = "Active"
        fake_user_profile_name = "fake-user-profile-name"
        fake_domain_id = "fake-domain-id"
        fake_user_profile_arn = arn_formatter(
            "user-profile",
            f"{fake_domain_id}/{fake_user_profile_name}",
            account_id,
            region_name,
        )
        self.created_by = {
            "UserProfileArn": fake_user_profile_arn,
            "UserProfileName": fake_user_profile_name,
            "DomainId": fake_domain_id,
        }
        self.last_modified_by = {
            "UserProfileArn": fake_user_profile_arn,
            "UserProfileName": fake_user_profile_name,
            "DomainId": fake_domain_id,
        }


class FakeProcessingJob(BaseObject):
    def __init__(
        self,
        app_specification: Dict[str, Any],
        experiment_config: Dict[str, str],
        network_config: Dict[str, Any],
        processing_inputs: List[Dict[str, Any]],
        processing_job_name: str,
        processing_output_config: Dict[str, Any],
        account_id: str,
        region_name: str,
        role_arn: str,
        tags: List[Dict[str, str]],
        stopping_condition: Dict[str, int],
    ):
        self.processing_job_name = processing_job_name
        self.arn = FakeProcessingJob.arn_formatter(
            processing_job_name, account_id, region_name
        )
        now_string = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        self.creation_time = now_string
        self.last_modified_time = now_string
        self.processing_end_time = now_string
        self.tags = tags or []
        self.role_arn = role_arn
        self.app_specification = app_specification
        self.experiment_config = experiment_config
        self.network_config = network_config
        self.processing_inputs = processing_inputs
        self.processing_job_status = "Completed"
        self.processing_output_config = processing_output_config
        self.stopping_condition = stopping_condition

    @property
    def response_object(self) -> Dict[str, Any]:  # type: ignore[misc]
        response_object = self.gen_response_object()
        response = {
            k: v for k, v in response_object.items() if v is not None and v != [None]
        }
        response["ProcessingJobArn"] = response.pop("Arn")

        return response

    @property
    def response_create(self) -> Dict[str, str]:
        return {"ProcessingJobArn": self.arn}

    @staticmethod
    def arn_formatter(name: str, account_id: str, region: str) -> str:
        return arn_formatter("processing-job", name, account_id, region)


class FakeTrainingJob(BaseObject):
    def __init__(
        self,
        account_id: str,
        region_name: str,
        training_job_name: str,
        hyper_parameters: Dict[str, str],
        algorithm_specification: Dict[str, Any],
        role_arn: str,
        input_data_config: List[Dict[str, Any]],
        output_data_config: Dict[str, str],
        resource_config: Dict[str, Any],
        vpc_config: Dict[str, List[str]],
        stopping_condition: Dict[str, int],
        tags: List[Dict[str, str]],
        enable_network_isolation: bool,
        enable_inter_container_traffic_encryption: bool,
        enable_managed_spot_training: bool,
        checkpoint_config: Dict[str, str],
        debug_hook_config: Dict[str, Any],
        debug_rule_configurations: List[Dict[str, Any]],
        tensor_board_output_config: Dict[str, str],
        experiment_config: Dict[str, str],
    ):
        self.training_job_name = training_job_name
        self.hyper_parameters = hyper_parameters
        self.algorithm_specification = algorithm_specification
        self.role_arn = role_arn
        self.input_data_config = input_data_config
        self.output_data_config = output_data_config
        self.resource_config = resource_config
        self.vpc_config = vpc_config
        self.stopping_condition = stopping_condition
        self.tags = tags or []
        self.enable_network_isolation = enable_network_isolation
        self.enable_inter_container_traffic_encryption = (
            enable_inter_container_traffic_encryption
        )
        self.enable_managed_spot_training = enable_managed_spot_training
        self.checkpoint_config = checkpoint_config
        self.debug_hook_config = debug_hook_config
        self.debug_rule_configurations = debug_rule_configurations
        self.tensor_board_output_config = tensor_board_output_config
        self.experiment_config = experiment_config
        self.arn = FakeTrainingJob.arn_formatter(
            training_job_name, account_id, region_name
        )
        self.creation_time = self.last_modified_time = datetime.now().strftime(
            "%Y-%m-%d %H:%M:%S"
        )
        self.model_artifacts = {
            "S3ModelArtifacts": os.path.join(
                self.output_data_config["S3OutputPath"],
                self.training_job_name,
                "output",
                "model.tar.gz",
            )
        }
        self.training_job_status = "Completed"
        self.secondary_status = "Completed"
        self.algorithm_specification["MetricDefinitions"] = [
            {
                "Name": "test:dcg",
                "Regex": "#quality_metric: host=\\S+, test dcg <score>=(\\S+)",
            }
        ]
        now_string = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        self.creation_time = now_string
        self.last_modified_time = now_string
        self.training_start_time = now_string
        self.training_end_time = now_string
        self.secondary_status_transitions = [
            {
                "Status": "Starting",
                "StartTime": self.creation_time,
                "EndTime": self.creation_time,
                "StatusMessage": "Preparing the instances for training",
            }
        ]
        self.final_metric_data_list = [
            {
                "MetricName": "train:progress",
                "Value": 100.0,
                "Timestamp": self.creation_time,
            }
        ]

    @property
    def response_object(self) -> Dict[str, Any]:  # type: ignore[misc]
        response_object = self.gen_response_object()
        response = {
            k: v for k, v in response_object.items() if v is not None and v != [None]
        }
        response["TrainingJobArn"] = response.pop("Arn")

        return response

    @property
    def response_create(self) -> Dict[str, str]:
        return {"TrainingJobArn": self.arn}

    @staticmethod
    def arn_formatter(name: str, account_id: str, region_name: str) -> str:
        return arn_formatter("training-job", name, account_id, region_name)


class FakeEndpoint(BaseObject, CloudFormationModel):
    def __init__(
        self,
        account_id: str,
        region_name: str,
        endpoint_name: str,
        endpoint_config_name: str,
        production_variants: List[Dict[str, Any]],
        data_capture_config: Dict[str, Any],
        tags: List[Dict[str, str]],
    ):
        self.endpoint_name = endpoint_name
        self.arn = FakeEndpoint.arn_formatter(endpoint_name, account_id, region_name)
        self.endpoint_config_name = endpoint_config_name
        self.production_variants = self._process_production_variants(
            production_variants
        )
        self.data_capture_config = data_capture_config
        self.tags = tags or []
        self.endpoint_status = "InService"
        self.failure_reason = None
        self.creation_time = self.last_modified_time = datetime.now().strftime(
            "%Y-%m-%d %H:%M:%S"
        )

    def _process_production_variants(
        self, production_variants: List[Dict[str, Any]]
    ) -> List[Dict[str, Any]]:
        endpoint_variants = []
        for production_variant in production_variants:
            temp_variant = {}

            # VariantName is the only required param
            temp_variant["VariantName"] = production_variant["VariantName"]

            if production_variant.get("InitialInstanceCount", None):
                temp_variant["CurrentInstanceCount"] = production_variant[
                    "InitialInstanceCount"
                ]
                temp_variant["DesiredInstanceCount"] = production_variant[
                    "InitialInstanceCount"
                ]

            if production_variant.get("InitialVariantWeight", None):
                temp_variant["CurrentWeight"] = production_variant[
                    "InitialVariantWeight"
                ]
                temp_variant["DesiredWeight"] = production_variant[
                    "InitialVariantWeight"
                ]

            if production_variant.get("ServerlessConfig", None):
                temp_variant["CurrentServerlessConfig"] = production_variant[
                    "ServerlessConfig"
                ]
                temp_variant["DesiredServerlessConfig"] = production_variant[
                    "ServerlessConfig"
                ]

            endpoint_variants.append(temp_variant)

        return endpoint_variants

    def summary(self) -> Dict[str, Any]:
        return {
            "EndpointName": self.endpoint_name,
            "EndpointArn": self.arn,
            "CreationTime": self.creation_time,
            "LastModifiedTime": self.last_modified_time,
            "EndpointStatus": self.endpoint_status,
        }

    @property
    def response_object(self) -> Dict[str, Any]:  # type: ignore[misc]
        response_object = self.gen_response_object()
        response = {
            k: v for k, v in response_object.items() if v is not None and v != [None]
        }
        response["EndpointArn"] = response.pop("Arn")

        return response

    @property
    def response_create(self) -> Dict[str, str]:
        return {"EndpointArn": self.arn}

    @staticmethod
    def arn_formatter(endpoint_name: str, account_id: str, region_name: str) -> str:
        return arn_formatter("endpoint", endpoint_name, account_id, region_name)

    @property
    def physical_resource_id(self) -> str:
        return self.arn

    @classmethod
    def has_cfn_attr(cls, attr: str) -> bool:
        return attr in ["EndpointName"]

    def get_cfn_attribute(self, attribute_name: str) -> str:
        # https://docs.aws.amazon.com/AWSCloudFormation/latest/UserGuide/aws-resource-sagemaker-endpoint.html#aws-resource-sagemaker-endpoint-return-values
        from moto.cloudformation.exceptions import UnformattedGetAttTemplateException

        if attribute_name == "EndpointName":
            return self.endpoint_name
        raise UnformattedGetAttTemplateException()

    @staticmethod
    def cloudformation_name_type() -> str:
        return ""

    @staticmethod
    def cloudformation_type() -> str:
        # https://docs.aws.amazon.com/AWSCloudFormation/latest/UserGuide/aws-resource-sagemaker-endpoint.html
        return "AWS::SageMaker::Endpoint"

    @classmethod
    def create_from_cloudformation_json(  # type: ignore[misc]
        cls,
        resource_name: str,
        cloudformation_json: Any,
        account_id: str,
        region_name: str,
        **kwargs: Any,
    ) -> "FakeEndpoint":
        sagemaker_backend = sagemaker_backends[account_id][region_name]

        # Get required properties from provided CloudFormation template
        properties = cloudformation_json["Properties"]
        endpoint_config_name = properties["EndpointConfigName"]

        endpoint = sagemaker_backend.create_endpoint(
            endpoint_name=resource_name,
            endpoint_config_name=endpoint_config_name,
            tags=properties.get("Tags", []),
        )
        return endpoint

    @classmethod
    def update_from_cloudformation_json(  # type: ignore[misc]
        cls,
        original_resource: Any,
        new_resource_name: str,
        cloudformation_json: Any,
        account_id: str,
        region_name: str,
    ) -> "FakeEndpoint":
        # Changes to the Endpoint will not change resource name
        cls.delete_from_cloudformation_json(
            original_resource.arn, cloudformation_json, account_id, region_name
        )
        new_resource = cls.create_from_cloudformation_json(
            original_resource.endpoint_name,
            cloudformation_json,
            account_id,
            region_name,
        )
        return new_resource

    @classmethod
    def delete_from_cloudformation_json(  # type: ignore[misc]
        cls,
        resource_name: str,
        cloudformation_json: Any,
        account_id: str,
        region_name: str,
    ) -> None:
        # Get actual name because resource_name actually provides the ARN
        # since the Physical Resource ID is the ARN despite SageMaker
        # using the name for most of its operations.
        endpoint_name = resource_name.split("/")[-1]

        sagemaker_backends[account_id][region_name].delete_endpoint(endpoint_name)


class FakeEndpointConfig(BaseObject, CloudFormationModel):
    def __init__(
        self,
        account_id: str,
        region_name: str,
        endpoint_config_name: str,
        production_variants: List[Dict[str, Any]],
        data_capture_config: Dict[str, Any],
        tags: List[Dict[str, Any]],
        kms_key_id: str,
    ):
        self.validate_production_variants(production_variants)

        self.endpoint_config_name = endpoint_config_name
        self.endpoint_config_arn = FakeEndpointConfig.arn_formatter(
            endpoint_config_name, account_id, region_name
        )
        self.arn = (self.endpoint_config_arn,)
        self.production_variants = production_variants or []
        self.data_capture_config = data_capture_config or {}
        self.tags = tags or []
        self.kms_key_id = kms_key_id
        self.creation_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    def validate_production_variants(
        self, production_variants: List[Dict[str, Any]]
    ) -> None:
        for production_variant in production_variants:
            if "InstanceType" in production_variant.keys():
                self.validate_instance_type(production_variant["InstanceType"])
            elif "ServerlessConfig" in production_variant.keys():
                self.validate_serverless_config(production_variant["ServerlessConfig"])
            else:
                message = f"Invalid Keys for ProductionVariant: received {production_variant.keys()} but expected it to contain one of {['InstanceType', 'ServerlessConfig']}"
                raise ValidationError(message=message)

    def validate_serverless_config(self, serverless_config: Dict[str, Any]) -> None:
        VALID_SERVERLESS_MEMORY_SIZE = [1024, 2048, 3072, 4096, 5120, 6144]
        if not validators.is_one_of(
            serverless_config["MemorySizeInMB"], VALID_SERVERLESS_MEMORY_SIZE
        ):
            message = f"Value '{serverless_config['MemorySizeInMB']}' at 'MemorySizeInMB' failed to satisfy constraint: Member must satisfy enum value set: {VALID_SERVERLESS_MEMORY_SIZE}"
            raise ValidationError(message=message)

    def validate_instance_type(self, instance_type: str) -> None:
        VALID_INSTANCE_TYPES = [
            "ml.r5d.12xlarge",
            "ml.r5.12xlarge",
            "ml.p2.xlarge",
            "ml.m5.4xlarge",
            "ml.m4.16xlarge",
            "ml.r5d.24xlarge",
            "ml.r5.24xlarge",
            "ml.p3.16xlarge",
            "ml.m5d.xlarge",
            "ml.m5.large",
            "ml.t2.xlarge",
            "ml.p2.16xlarge",
            "ml.m5d.12xlarge",
            "ml.inf1.2xlarge",
            "ml.m5d.24xlarge",
            "ml.c4.2xlarge",
            "ml.c5.2xlarge",
            "ml.c4.4xlarge",
            "ml.inf1.6xlarge",
            "ml.c5d.2xlarge",
            "ml.c5.4xlarge",
            "ml.g4dn.xlarge",
            "ml.g4dn.12xlarge",
            "ml.c5d.4xlarge",
            "ml.g4dn.2xlarge",
            "ml.c4.8xlarge",
            "ml.c4.large",
            "ml.c5d.xlarge",
            "ml.c5.large",
            "ml.g4dn.4xlarge",
            "ml.c5.9xlarge",
            "ml.g4dn.16xlarge",
            "ml.c5d.large",
            "ml.c5.xlarge",
            "ml.c5d.9xlarge",
            "ml.c4.xlarge",
            "ml.inf1.xlarge",
            "ml.g4dn.8xlarge",
            "ml.inf1.24xlarge",
            "ml.m5d.2xlarge",
            "ml.t2.2xlarge",
            "ml.c5d.18xlarge",
            "ml.m5d.4xlarge",
            "ml.t2.medium",
            "ml.c5.18xlarge",
            "ml.r5d.2xlarge",
            "ml.r5.2xlarge",
            "ml.p3.2xlarge",
            "ml.m5d.large",
            "ml.m5.xlarge",
            "ml.m4.10xlarge",
            "ml.t2.large",
            "ml.r5d.4xlarge",
            "ml.r5.4xlarge",
            "ml.m5.12xlarge",
            "ml.m4.xlarge",
            "ml.m5.24xlarge",
            "ml.m4.2xlarge",
            "ml.p2.8xlarge",
            "ml.m5.2xlarge",
            "ml.r5d.xlarge",
            "ml.r5d.large",
            "ml.r5.xlarge",
            "ml.r5.large",
            "ml.p3.8xlarge",
            "ml.m4.4xlarge",
        ]
        if not validators.is_one_of(instance_type, VALID_INSTANCE_TYPES):
            message = f"Value '{instance_type}' at 'instanceType' failed to satisfy constraint: Member must satisfy enum value set: {VALID_INSTANCE_TYPES}"
            raise ValidationError(message=message)

    def summary(self) -> Dict[str, Any]:
        return {
            "EndpointConfigName": self.endpoint_config_name,
            "EndpointConfigArn": self.endpoint_config_arn,
            "CreationTime": self.creation_time,
        }

    @property
    def response_object(self) -> Dict[str, Any]:  # type: ignore[misc]
        response_object = self.gen_response_object()
        return {
            k: v for k, v in response_object.items() if v is not None and v != [None]
        }

    @property
    def response_create(self) -> Dict[str, str]:
        return {"EndpointConfigArn": self.endpoint_config_arn}

    @staticmethod
    def arn_formatter(
        endpoint_config_name: str, account_id: str, region_name: str
    ) -> str:
        return arn_formatter(
            "endpoint-config", endpoint_config_name, account_id, region_name
        )

    @property
    def physical_resource_id(self) -> str:
        return self.endpoint_config_arn

    @classmethod
    def has_cfn_attr(cls, attr: str) -> bool:
        return attr in ["EndpointConfigName"]

    def get_cfn_attribute(self, attribute_name: str) -> str:
        # https://docs.aws.amazon.com/AWSCloudFormation/latest/UserGuide/aws-resource-sagemaker-endpointconfig.html#aws-resource-sagemaker-endpointconfig-return-values
        from moto.cloudformation.exceptions import UnformattedGetAttTemplateException

        if attribute_name == "EndpointConfigName":
            return self.endpoint_config_name
        raise UnformattedGetAttTemplateException()

    @staticmethod
    def cloudformation_name_type() -> str:
        return ""

    @staticmethod
    def cloudformation_type() -> str:
        # https://docs.aws.amazon.com/AWSCloudFormation/latest/UserGuide/aws-resource-sagemaker-endpointconfig.html
        return "AWS::SageMaker::EndpointConfig"

    @classmethod
    def create_from_cloudformation_json(  # type: ignore[misc]
        cls,
        resource_name: str,
        cloudformation_json: Any,
        account_id: str,
        region_name: str,
        **kwargs: Any,
    ) -> "FakeEndpointConfig":
        sagemaker_backend = sagemaker_backends[account_id][region_name]

        # Get required properties from provided CloudFormation template
        properties = cloudformation_json["Properties"]
        production_variants = properties["ProductionVariants"]

        endpoint_config = sagemaker_backend.create_endpoint_config(
            endpoint_config_name=resource_name,
            production_variants=production_variants,
            data_capture_config=properties.get("DataCaptureConfig", {}),
            kms_key_id=properties.get("KmsKeyId"),
            tags=properties.get("Tags", []),
        )
        return endpoint_config

    @classmethod
    def update_from_cloudformation_json(  # type: ignore[misc]
        cls,
        original_resource: Any,
        new_resource_name: str,
        cloudformation_json: Any,
        account_id: str,
        region_name: str,
    ) -> "FakeEndpointConfig":
        # Most changes to the endpoint config will change resource name for EndpointConfigs
        cls.delete_from_cloudformation_json(
            original_resource.endpoint_config_arn,
            cloudformation_json,
            account_id,
            region_name,
        )
        new_resource = cls.create_from_cloudformation_json(
            new_resource_name, cloudformation_json, account_id, region_name
        )
        return new_resource

    @classmethod
    def delete_from_cloudformation_json(  # type: ignore[misc]
        cls,
        resource_name: str,
        cloudformation_json: Any,
        account_id: str,
        region_name: str,
    ) -> None:
        # Get actual name because resource_name actually provides the ARN
        # since the Physical Resource ID is the ARN despite SageMaker
        # using the name for most of its operations.
        endpoint_config_name = resource_name.split("/")[-1]

        sagemaker_backends[account_id][region_name].delete_endpoint_config(
            endpoint_config_name
        )


class FakeTransformJob(BaseObject):
    def __init__(
        self,
        account_id: str,
        region_name: str,
        transform_job_name: str,
        model_name: str,
        max_concurrent_transforms: int,
        model_client_config: Dict[str, int],
        max_payload_in_mb: int,
        batch_strategy: str,
        environment: Dict[str, str],
        transform_input: Dict[str, Union[Dict[str, str], str]],
        transform_output: Dict[str, str],
        data_capture_config: Dict[str, Union[str, bool]],
        transform_resources: Dict[str, Union[str, int]],
        data_processing: Dict[str, str],
        tags: Dict[str, str],
        experiment_config: Dict[str, str],
    ):
        self.transform_job_name = transform_job_name
        self.model_name = model_name
        self.max_concurrent_transforms = max_concurrent_transforms
        self.model_client_config = model_client_config
        self.max_payload_in_mb = max_payload_in_mb
        self.batch_strategy = batch_strategy
        self.environment = environment
        self.transform_input = transform_input
        self.transform_output = transform_output
        self.data_capture_config = data_capture_config
        self.transform_resources = transform_resources
        self.data_processing = data_processing
        self.tags = tags
        self.experiment_config = experiment_config
        self.arn = FakeTransformJob.arn_formatter(
            transform_job_name, account_id, region_name
        )
        self.transform_job_status = "Completed"
        self.failure_reason = ""
        self.labeling_job_arn = ""
        self.auto_ml_job_arn = ""
        now_string = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        self.creation_time = now_string
        self.transform_start_time = now_string
        self.transform_end_time = now_string
        self.last_modified_time = now_string

    # Override title case
    def camelCase(self, key: str) -> str:
        words = []
        for word in key.split("_"):
            if word == "mb":
                words.append("MB")
            else:
                words.append(word.title())
        return "".join(words)

    @property
    def response_object(self) -> Dict[str, Any]:  # type: ignore[misc]
        response_object = self.gen_response_object()
        response = {
            k: v for k, v in response_object.items() if v is not None and v != [None]
        }
        return response

    @property
    def response_create(self) -> Dict[str, str]:
        return {"TransformJobArn": self.arn}

    @staticmethod
    def arn_formatter(name: str, account_id: str, region_name: str) -> str:
        return arn_formatter("transform-job", name, account_id, region_name)


class Model(BaseObject, CloudFormationModel):
    def __init__(
        self,
        account_id: str,
        region_name: str,
        model_name: str,
        execution_role_arn: str,
        primary_container: Dict[str, Any],
        vpc_config: Dict[str, Any],
        containers: Optional[List[Dict[str, Any]]] = None,
        tags: Optional[List[Dict[str, str]]] = None,
    ):
        self.model_name = model_name
        self.creation_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        self.containers = containers or []
        self.tags = tags or []
        self.enable_network_isolation = False
        self.vpc_config = vpc_config
        self.primary_container = primary_container
        self.execution_role_arn = execution_role_arn or "arn:test"
        self.arn = arn_formatter("model", self.model_name, account_id, region_name)

    @property
    def response_object(self) -> Dict[str, Any]:  # type: ignore[misc]
        response_object = self.gen_response_object()
        response = {
            k: v for k, v in response_object.items() if v is not None and v != [None]
        }
        response["ModelArn"] = response.pop("Arn")

        return response

    @property
    def response_create(self) -> Dict[str, str]:
        return {"ModelArn": self.arn}

    @property
    def physical_resource_id(self) -> str:
        return self.arn

    @classmethod
    def has_cfn_attr(cls, attr: str) -> bool:
        return attr in ["ModelName"]

    def get_cfn_attribute(self, attribute_name: str) -> str:
        # https://docs.aws.amazon.com/AWSCloudFormation/latest/UserGuide/aws-resource-sagemaker-model.html#aws-resource-sagemaker-model-return-values
        from moto.cloudformation.exceptions import UnformattedGetAttTemplateException

        if attribute_name == "ModelName":
            return self.model_name
        raise UnformattedGetAttTemplateException()

    @staticmethod
    def cloudformation_name_type() -> str:
        return ""

    @staticmethod
    def cloudformation_type() -> str:
        # https://docs.aws.amazon.com/AWSCloudFormation/latest/UserGuide/aws-resource-sagemaker-model.html
        return "AWS::SageMaker::Model"

    @classmethod
    def create_from_cloudformation_json(  # type: ignore[misc]
        cls,
        resource_name: str,
        cloudformation_json: Any,
        account_id: str,
        region_name: str,
        **kwargs: Any,
    ) -> "Model":
        sagemaker_backend = sagemaker_backends[account_id][region_name]

        # Get required properties from provided CloudFormation template
        properties = cloudformation_json["Properties"]
        execution_role_arn = properties["ExecutionRoleArn"]
        primary_container = properties["PrimaryContainer"]

        model = sagemaker_backend.create_model(
            model_name=resource_name,
            execution_role_arn=execution_role_arn,
            primary_container=primary_container,
            vpc_config=properties.get("VpcConfig", {}),
            containers=properties.get("Containers", []),
            tags=properties.get("Tags", []),
        )
        return model

    @classmethod
    def update_from_cloudformation_json(  # type: ignore[misc]
        cls,
        original_resource: Any,
        new_resource_name: str,
        cloudformation_json: Any,
        account_id: str,
        region_name: str,
    ) -> "Model":
        # Most changes to the model will change resource name for Models
        cls.delete_from_cloudformation_json(
            original_resource.arn, cloudformation_json, account_id, region_name
        )
        new_resource = cls.create_from_cloudformation_json(
            new_resource_name, cloudformation_json, account_id, region_name
        )
        return new_resource

    @classmethod
    def delete_from_cloudformation_json(  # type: ignore[misc]
        cls,
        resource_name: str,
        cloudformation_json: Any,
        account_id: str,
        region_name: str,
    ) -> None:
        # Get actual name because resource_name actually provides the ARN
        # since the Physical Resource ID is the ARN despite SageMaker
        # using the name for most of its operations.
        model_name = resource_name.split("/")[-1]

        sagemaker_backends[account_id][region_name].delete_model(model_name)


class ModelPackageGroup(BaseObject):
    def __init__(
        self,
        model_package_group_name: str,
        model_package_group_description: str,
        account_id: str,
        region_name: str,
        tags: Optional[List[Dict[str, str]]] = None,
    ) -> None:
        model_package_group_arn = arn_formatter(
            region_name=region_name,
            account_id=account_id,
            _type="model-package-group",
            _id=model_package_group_name,
        )
        fake_user_profile_name = "fake-user-profile-name"
        fake_domain_id = "fake-domain-id"
        fake_user_profile_arn = arn_formatter(
            _type="user-profile",
            _id=f"{fake_domain_id}/{fake_user_profile_name}",
            account_id=account_id,
            region_name=region_name,
        )
        datetime_now = datetime.now(tzutc())
        self.model_package_group_name = model_package_group_name
        self.arn = model_package_group_arn
        self.model_package_group_description = model_package_group_description
        self.creation_time = datetime_now
        self.created_by = {
            "UserProfileArn": fake_user_profile_arn,
            "UserProfileName": fake_user_profile_name,
            "DomainId": fake_domain_id,
        }
        self.model_package_group_status = "Completed"
        self.tags = tags

    def gen_response_object(self) -> Dict[str, Any]:
        response_object = super().gen_response_object()
        for k, v in response_object.items():
            if isinstance(v, datetime):
                response_object[k] = v.isoformat()
        response_values = [
            "ModelPackageGroupName",
            "Arn",
            "ModelPackageGroupDescription",
            "CreationTime",
            "ModelPackageGroupStatus",
            "Tags",
        ]
        response = {k: v for k, v in response_object.items() if k in response_values}
        response["ModelPackageGroupArn"] = response.pop("Arn")

        return response


class FakeModelCard(BaseObject):
    def __init__(
        self,
        account_id: str,
        region_name: str,
        model_card_name: str,
        model_card_version: int,
        content: str,
        model_card_status: str,
        security_config: Optional[Dict[str, str]] = None,
        tags: Optional[List[Dict[str, Any]]] = None,
        creation_time: Optional[str] = None,
        last_modified_time: Optional[str] = None,
    ) -> None:
        datetime_now = str(datetime.now(tzutc()))
        self.arn = arn_formatter("model-card", model_card_name, account_id, region_name)
        self.model_card_name = model_card_name
        self.model_card_version = model_card_version
        self.content = content
        self.model_card_status = model_card_status
        self.creation_time = creation_time if creation_time else datetime_now
        self.last_modified_time = (
            last_modified_time if last_modified_time else datetime_now
        )
        self.security_config = security_config
        self.tags = tags

    def describe(self) -> Dict[str, Any]:
        return {
            "ModelCardArn": self.arn,
            "ModelCardName": self.model_card_name,
            "ModelCardVersion": self.model_card_version,
            "Content": self.content,
            "ModelCardStatus": self.model_card_status,
            "SecurityConfig": self.security_config,
            "CreationTime": self.creation_time,
            "CreatedBy": {},
            "LastModifiedTime": self.creation_time,
            "LastModifiedBy": {},
        }

    def summary(self) -> Dict[str, Any]:
        return {
            "ModelCardName": self.model_card_name,
            "ModelCardArn": self.arn,
            "ModelCardStatus": self.model_card_status,
            "CreationTime": self.creation_time,
            "LastModifiedTime": self.last_modified_time,
        }

    def version_summary(self) -> Dict[str, Any]:
        return {
            "ModelCardName": self.model_card_name,
            "ModelCardArn": self.arn,
            "ModelCardStatus": self.model_card_status,
            "ModelCardVersion": self.model_card_version,
            "CreationTime": self.creation_time,
            "LastModifiedTime": self.last_modified_time,
        }


class FeatureGroup(BaseObject):
    def __init__(
        self,
        region_name: str,
        account_id: str,
        feature_group_name: str,
        record_identifier_feature_name: str,
        event_time_feature_name: str,
        feature_definitions: List[Dict[str, str]],
        offline_store_config: Dict[str, Any],
        role_arn: str,
        tags: Optional[List[Dict[str, str]]] = None,
    ) -> None:
        self.feature_group_name = feature_group_name
        self.record_identifier_feature_name = record_identifier_feature_name
        self.event_time_feature_name = event_time_feature_name
        self.feature_definitions = feature_definitions

        table_name = (
            f"{feature_group_name.replace('-','_')}_{int(datetime.now().timestamp())}"
        )
        offline_store_config["DataCatalogConfig"] = {
            "TableName": table_name,
            "Catalog": "AwsDataCatalog",
            "Database": "sagemaker_featurestore",
        }
        offline_store_config["S3StorageConfig"]["ResolvedOutputS3Uri"] = (
            f'{offline_store_config["S3StorageConfig"]["S3Uri"]}/{account_id}/{region_name}/offline-store/{feature_group_name}-{int(datetime.now().timestamp())}/data'
        )

        self.offline_store_config = offline_store_config
        self.role_arn = role_arn

        self.creation_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        self.arn = arn_formatter(
            region_name=region_name,
            account_id=account_id,
            _type="feature-group",
            _id=f"{self.feature_group_name.lower()}",
        )
        self.tags = tags

    def describe(self) -> Dict[str, Any]:
        return {
            "FeatureGroupArn": self.arn,
            "FeatureGroupName": self.feature_group_name,
            "RecordIdentifierFeatureName": self.record_identifier_feature_name,
            "EventTimeFeatureName": self.event_time_feature_name,
            "FeatureDefinitions": self.feature_definitions,
            "CreationTime": self.creation_time,
            "OfflineStoreConfig": self.offline_store_config,
            "RoleArn": self.role_arn,
            "ThroughputConfig": {"ThroughputMode": "OnDemand"},
            "FeatureGroupStatus": "Created",
        }


class ModelPackage(BaseObject):
    def __init__(
        self,
        model_package_name: str,
        model_package_group_name: Optional[str],
        model_package_version: Optional[int],
        model_package_description: Optional[str],
        inference_specification: Any,
        source_algorithm_specification: Any,
        validation_specification: Any,
        certify_for_marketplace: bool,
        model_approval_status: Optional[str],
        metadata_properties: Any,
        model_metrics: Any,
        approval_description: Optional[str],
        customer_metadata_properties: Any,
        drift_check_baselines: Any,
        domain: str,
        task: str,
        sample_payload_url: str,
        additional_inference_specifications: List[Any],
        client_token: str,
        region_name: str,
        account_id: str,
        model_package_type: str,
        tags: Optional[List[Dict[str, str]]] = None,
    ) -> None:
        fake_user_profile_name = "fake-user-profile-name"
        fake_domain_id = "fake-domain-id"
        fake_user_profile_arn = arn_formatter(
            _type="user-profile",
            _id=f"{fake_domain_id}/{fake_user_profile_name}",
            account_id=account_id,
            region_name=region_name,
        )
        model_package_arn = arn_formatter(
            region_name=region_name,
            account_id=account_id,
            _type="model-package",
            _id=f"{model_package_name.lower()}/{model_package_version}"
            if model_package_version
            else model_package_name.lower(),
        )
        datetime_now = datetime.now(tzutc())
        self.model_package_name = model_package_name
        self.model_package_group_name = model_package_group_name
        self.model_package_version = model_package_version
        self.arn = model_package_arn
        self.model_package_description = model_package_description
        self.creation_time = datetime_now
        self.inference_specification = inference_specification
        self.source_algorithm_specification = source_algorithm_specification
        self.validation_specification = validation_specification
        self.model_package_type = model_package_type
        self.model_package_status_details = {
            "ValidationStatuses": [
                {
                    "Name": model_package_arn,
                    "Status": "Completed",
                }
            ],
            "ImageScanStatuses": [
                {
                    "Name": model_package_arn,
                    "Status": "Completed",
                }
            ],
        }
        self.certify_for_marketplace = certify_for_marketplace
        self.model_approval_status: Optional[str] = None
        self.set_model_approval_status(model_approval_status)
        self.created_by = {
            "UserProfileArn": fake_user_profile_arn,
            "UserProfileName": fake_user_profile_name,
            "DomainId": fake_domain_id,
        }
        self.metadata_properties = metadata_properties
        self.model_metrics = model_metrics
        self.last_modified_time: Optional[datetime] = None
        self.approval_description = approval_description
        self.customer_metadata_properties = customer_metadata_properties
        self.drift_check_baselines = drift_check_baselines
        self.domain = domain
        self.task = task
        self.sample_payload_url = sample_payload_url
        self.additional_inference_specifications: Optional[List[Any]] = None
        self.add_additional_inference_specifications(
            additional_inference_specifications
        )
        self.tags = tags
        self.model_package_status = "Completed"
        self.last_modified_by: Optional[Dict[str, str]] = None
        self.client_token = client_token

    def gen_response_object(self) -> Dict[str, Any]:
        response_object = super().gen_response_object()
        for k, v in response_object.items():
            if isinstance(v, datetime):
                response_object[k] = v.isoformat()
        response_values = [
            "ModelPackageName",
            "ModelPackageGroupName",
            "ModelPackageVersion",
            "Arn",
            "ModelPackageDescription",
            "CreationTime",
            "InferenceSpecification",
            "SourceAlgorithmSpecification",
            "ValidationSpecification",
            "ModelPackageStatus",
            "ModelPackageStatusDetails",
            "CertifyForMarketplace",
            "ModelApprovalStatus",
            "CreatedBy",
            "MetadataProperties",
            "ModelMetrics",
            "LastModifiedTime",
            "LastModifiedBy",
            "ApprovalDescription",
            "CustomerMetadataProperties",
            "DriftCheckBaselines",
            "Domain",
            "Task",
            "SamplePayloadUrl",
            "AdditionalInferenceSpecifications",
            "SkipModelValidation",
        ]
        if self.model_package_type == "Versioned":
            del response_object["ModelPackageName"]
        elif self.model_package_type == "Unversioned":
            del response_object["ModelPackageGroupName"]
        response = {
            k: v
            for k, v in response_object.items()
            if k in response_values
            if v is not None
        }
        response["ModelPackageArn"] = response.pop("Arn")

        return response

    def modifications_done(self) -> None:
        self.last_modified_time = datetime.now(tzutc())
        self.last_modified_by = self.created_by

    def set_model_approval_status(self, model_approval_status: Optional[str]) -> None:
        if model_approval_status is not None:
            validate_model_approval_status(model_approval_status)
        self.model_approval_status = model_approval_status

    def remove_customer_metadata_property(
        self, customer_metadata_properties_to_remove: List[str]
    ) -> None:
        if customer_metadata_properties_to_remove is not None:
            for customer_metadata_property in customer_metadata_properties_to_remove:
                self.customer_metadata_properties.pop(customer_metadata_property, None)

    def add_additional_inference_specifications(
        self, additional_inference_specifications_to_add: Optional[List[Any]]
    ) -> None:
        self.validate_additional_inference_specifications(
            additional_inference_specifications_to_add
        )
        if (
            self.additional_inference_specifications is not None
            and additional_inference_specifications_to_add is not None
        ):
            self.additional_inference_specifications.extend(
                additional_inference_specifications_to_add
            )
        else:
            self.additional_inference_specifications = (
                additional_inference_specifications_to_add
            )

    def validate_additional_inference_specifications(
        self, additional_inference_specifications: Optional[List[Dict[str, Any]]]
    ) -> None:
        specifications_to_validate = additional_inference_specifications or []
        for additional_inference_specification in specifications_to_validate:
            if "SupportedTransformInstanceTypes" in additional_inference_specification:
                self.validate_supported_transform_instance_types(
                    additional_inference_specification[
                        "SupportedTransformInstanceTypes"
                    ]
                )
            if (
                "SupportedRealtimeInferenceInstanceTypes"
                in additional_inference_specification
            ):
                self.validate_supported_realtime_inference_instance_types(
                    additional_inference_specification[
                        "SupportedRealtimeInferenceInstanceTypes"
                    ]
                )

    @staticmethod
    def validate_supported_transform_instance_types(instance_types: List[str]) -> None:
        VALID_TRANSFORM_INSTANCE_TYPES = [
            "ml.m4.xlarge",
            "ml.m4.2xlarge",
            "ml.m4.4xlarge",
            "ml.m4.10xlarge",
            "ml.m4.16xlarge",
            "ml.c4.xlarge",
            "ml.c4.2xlarge",
            "ml.c4.4xlarge",
            "ml.c4.8xlarge",
            "ml.p2.xlarge",
            "ml.p2.8xlarge",
            "ml.p2.16xlarge",
            "ml.p3.2xlarge",
            "ml.p3.8xlarge",
            "ml.p3.16xlarge",
            "ml.c5.xlarge",
            "ml.c5.2xlarge",
            "ml.c5.4xlarge",
            "ml.c5.9xlarge",
            "ml.c5.18xlarge",
            "ml.m5.large",
            "ml.m5.xlarge",
            "ml.m5.2xlarge",
            "ml.m5.4xlarge",
            "ml.m5.12xlarge",
            "ml.m5.24xlarge",
            "ml.g4dn.xlarge",
            "ml.g4dn.2xlarge",
            "ml.g4dn.4xlarge",
            "ml.g4dn.8xlarge",
            "ml.g4dn.12xlarge",
            "ml.g4dn.16xlarge",
        ]
        for instance_type in instance_types:
            if not validators.is_one_of(instance_type, VALID_TRANSFORM_INSTANCE_TYPES):
                message = f"Value '{instance_type}' at 'SupportedTransformInstanceTypes' failed to satisfy constraint: Member must satisfy enum value set: {VALID_TRANSFORM_INSTANCE_TYPES}"
                raise ValidationError(message=message)

    @staticmethod
    def validate_supported_realtime_inference_instance_types(
        instance_types: List[str],
    ) -> None:
        VALID_REALTIME_INFERENCE_INSTANCE_TYPES = [
            "ml.t2.medium",
            "ml.t2.large",
            "ml.t2.xlarge",
            "ml.t2.2xlarge",
            "ml.m4.xlarge",
            "ml.m4.2xlarge",
            "ml.m4.4xlarge",
            "ml.m4.10xlarge",
            "ml.m4.16xlarge",
            "ml.m5.large",
            "ml.m5.xlarge",
            "ml.m5.2xlarge",
            "ml.m5.4xlarge",
            "ml.m5.12xlarge",
            "ml.m5.24xlarge",
            "ml.m5d.large",
            "ml.m5d.xlarge",
            "ml.m5d.2xlarge",
            "ml.m5d.4xlarge",
            "ml.m5d.12xlarge",
            "ml.m5d.24xlarge",
            "ml.c4.large",
            "ml.c4.xlarge",
            "ml.c4.2xlarge",
            "ml.c4.4xlarge",
            "ml.c4.8xlarge",
            "ml.p2.xlarge",
            "ml.p2.8xlarge",
            "ml.p2.16xlarge",
            "ml.p3.2xlarge",
            "ml.p3.8xlarge",
            "ml.p3.16xlarge",
            "ml.c5.large",
            "ml.c5.xlarge",
            "ml.c5.2xlarge",
            "ml.c5.4xlarge",
            "ml.c5.9xlarge",
            "ml.c5.18xlarge",
            "ml.c5d.large",
            "ml.c5d.xlarge",
            "ml.c5d.2xlarge",
            "ml.c5d.4xlarge",
            "ml.c5d.9xlarge",
            "ml.c5d.18xlarge",
            "ml.g4dn.xlarge",
            "ml.g4dn.2xlarge",
            "ml.g4dn.4xlarge",
            "ml.g4dn.8xlarge",
            "ml.g4dn.12xlarge",
            "ml.g4dn.16xlarge",
            "ml.r5.large",
            "ml.r5.xlarge",
            "ml.r5.2xlarge",
            "ml.r5.4xlarge",
            "ml.r5.12xlarge",
            "ml.r5.24xlarge",
            "ml.r5d.large",
            "ml.r5d.xlarge",
            "ml.r5d.2xlarge",
            "ml.r5d.4xlarge",
            "ml.r5d.12xlarge",
            "ml.r5d.24xlarge",
            "ml.inf1.xlarge",
            "ml.inf1.2xlarge",
            "ml.inf1.6xlarge",
            "ml.inf1.24xlarge",
            "ml.c6i.large",
            "ml.c6i.xlarge",
            "ml.c6i.2xlarge",
            "ml.c6i.4xlarge",
            "ml.c6i.8xlarge",
            "ml.c6i.12xlarge",
            "ml.c6i.16xlarge",
            "ml.c6i.24xlarge",
            "ml.c6i.32xlarge",
            "ml.g5.xlarge",
            "ml.g5.2xlarge",
            "ml.g5.4xlarge",
            "ml.g5.8xlarge",
            "ml.g5.12xlarge",
            "ml.g5.16xlarge",
            "ml.g5.24xlarge",
            "ml.g5.48xlarge",
            "ml.p4d.24xlarge",
            "ml.c7g.large",
            "ml.c7g.xlarge",
            "ml.c7g.2xlarge",
            "ml.c7g.4xlarge",
            "ml.c7g.8xlarge",
            "ml.c7g.12xlarge",
            "ml.c7g.16xlarge",
            "ml.m6g.large",
            "ml.m6g.xlarge",
            "ml.m6g.2xlarge",
            "ml.m6g.4xlarge",
            "ml.m6g.8xlarge",
            "ml.m6g.12xlarge",
            "ml.m6g.16xlarge",
            "ml.m6gd.large",
            "ml.m6gd.xlarge",
            "ml.m6gd.2xlarge",
            "ml.m6gd.4xlarge",
            "ml.m6gd.8xlarge",
            "ml.m6gd.12xlarge",
            "ml.m6gd.16xlarge",
            "ml.c6g.large",
            "ml.c6g.xlarge",
            "ml.c6g.2xlarge",
            "ml.c6g.4xlarge",
            "ml.c6g.8xlarge",
            "ml.c6g.12xlarge",
            "ml.c6g.16xlarge",
            "ml.c6gd.large",
            "ml.c6gd.xlarge",
            "ml.c6gd.2xlarge",
            "ml.c6gd.4xlarge",
            "ml.c6gd.8xlarge",
            "ml.c6gd.12xlarge",
            "ml.c6gd.16xlarge",
            "ml.c6gn.large",
            "ml.c6gn.xlarge",
            "ml.c6gn.2xlarge",
            "ml.c6gn.4xlarge",
            "ml.c6gn.8xlarge",
            "ml.c6gn.12xlarge",
            "ml.c6gn.16xlarge",
            "ml.r6g.large",
            "ml.r6g.xlarge",
            "ml.r6g.2xlarge",
            "ml.r6g.4xlarge",
            "ml.r6g.8xlarge",
            "ml.r6g.12xlarge",
            "ml.r6g.16xlarge",
            "ml.r6gd.large",
            "ml.r6gd.xlarge",
            "ml.r6gd.2xlarge",
            "ml.r6gd.4xlarge",
            "ml.r6gd.8xlarge",
            "ml.r6gd.12xlarge",
            "ml.r6gd.16xlarge",
            "ml.p4de.24xlarge",
            "ml.trn1.2xlarge",
            "ml.trn1.32xlarge",
            "ml.inf2.xlarge",
            "ml.inf2.8xlarge",
            "ml.inf2.24xlarge",
            "ml.inf2.48xlarge",
            "ml.p5.48xlarge",
        ]
        for instance_type in instance_types:
            if not validators.is_one_of(
                instance_type, VALID_REALTIME_INFERENCE_INSTANCE_TYPES
            ):
                message = f"Value '{instance_type}' at 'SupportedRealtimeInferenceInstanceTypes' failed to satisfy constraint: Member must satisfy enum value set: {VALID_REALTIME_INFERENCE_INSTANCE_TYPES}"
                raise ValidationError(message=message)


class Cluster(BaseObject):
    def __init__(
        self,
        cluster_name: str,
        region_name: str,
        account_id: str,
        instance_groups: List[Dict[str, Any]],
        vpc_config: Dict[str, List[str]],
        tags: Optional[List[Dict[str, str]]] = None,
    ):
        self.region_name = region_name
        self.account_id = account_id
        self.cluster_name = cluster_name
        if cluster_name in sagemaker_backends[account_id][region_name].clusters:
            raise ResourceInUseException(
                message=f"Resource Already Exists: Cluster with name {cluster_name} already exists. Choose a different name."
            )
        self.instance_groups = instance_groups
        for instance_group in instance_groups:
            self.valid_cluster_node_instance_types(instance_group["InstanceType"])
            if not instance_group["LifeCycleConfig"]["SourceS3Uri"].startswith(
                "s3://sagemaker-"
            ):
                raise ValidationError(
                    message=f"Validation Error: SourceS3Uri {instance_group['LifeCycleConfig']['SourceS3Uri']} does not start with 's3://sagemaker'."
                )
        self.vpc_config = vpc_config
        self.tags = tags or []
        self.arn = arn_formatter("cluster", self.cluster_name, account_id, region_name)
        self.status = "InService"
        self.creation_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        self.failure_message = ""
        self.nodes: Dict[str, ClusterNode] = {}
        for instance_group in self.instance_groups:
            instance_group["CurrentCount"] = instance_group["InstanceCount"]
            instance_group["TargetCount"] = instance_group["InstanceCount"]
            del instance_group["InstanceCount"]

    def describe(self) -> Dict[str, Any]:
        return {
            "ClusterArn": self.arn,
            "ClusterName": self.cluster_name,
            "ClusterStatus": self.status,
            "CreationTime": self.creation_time,
            "FailureMessage": self.failure_message,
            "InstanceGroups": self.instance_groups,
            "VpcConfig": self.vpc_config,
        }

    def summary(self) -> Dict[str, Any]:
        return {
            "ClusterArn": self.arn,
            "ClusterName": self.cluster_name,
            "CreationTime": self.creation_time,
            "ClusterStatus": self.status,
        }

    def valid_cluster_node_instance_types(self, instance_type: str) -> None:
        VALID_CLUSTER_INSTANCE_TYPES = [
            "ml.p4d.24xlarge",
            "ml.p4de.24xlarge",
            "ml.p5.48xlarge",
            "ml.trn1.32xlarge",
            "ml.trn1n.32xlarge",
            "ml.g5.xlarge",
            "ml.g5.2xlarge",
            "ml.g5.4xlarge",
            "ml.g5.8xlarge",
            "ml.g5.12xlarge",
            "ml.g5.16xlarge",
            "ml.g5.24xlarge",
            "ml.g5.48xlarge",
            "ml.c5.large",
            "ml.c5.xlarge",
            "ml.c5.2xlarge",
            "ml.c5.4xlarge",
            "ml.c5.9xlarge",
            "ml.c5.12xlarge",
            "ml.c5.18xlarge",
            "ml.c5.24xlarge",
            "ml.c5n.large",
            "ml.c5n.2xlarge",
            "ml.c5n.4xlarge",
            "ml.c5n.9xlarge",
            "ml.c5n.18xlarge",
            "ml.m5.large",
            "ml.m5.xlarge",
            "ml.m5.2xlarge",
            "ml.m5.4xlarge",
            "ml.m5.8xlarge",
            "ml.m5.12xlarge",
            "ml.m5.16xlarge",
            "ml.m5.24xlarge",
            "ml.t3.medium",
            "ml.t3.large",
            "ml.t3.xlarge",
            "ml.t3.2xlarge",
        ]

        if instance_type not in VALID_CLUSTER_INSTANCE_TYPES:
            message = f"Value '{instance_type}' at 'InstanceType' failed to satisfy constraint: Member must satisfy enum value set: {VALID_CLUSTER_INSTANCE_TYPES}"
            raise ValidationError(message=message)


class ClusterNode(BaseObject):
    def __init__(
        self,
        region_name: str,
        account_id: str,
        cluster_name: str,
        instance_group_name: str,
        instance_type: str,
        life_cycle_config: Dict[str, Any],
        execution_role: str,
        node_id: str,
        threads_per_core: Optional[int] = None,
    ):
        self.region_name = region_name
        self.account_id = account_id
        self.cluster_name = cluster_name
        self.instance_group_name = (
            instance_group_name  # probably need to do something with this
        )
        self.instance_id = node_id  # generate instance id
        self.instance_type = instance_type
        self.life_cycle_config = life_cycle_config
        self.execution_role = execution_role
        self.threads_per_core = threads_per_core
        self.status = "Running"
        self.launch_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    def describe(self) -> Dict[str, Any]:
        return {
            "InstanceGroupName": self.instance_group_name,
            "InstanceId": self.instance_id,
            "InstanceStatus": {"Status": self.status, "Message": "message"},
            "InstanceType": self.instance_type,
            "LaunchTime": self.launch_time,
            "LifeCycleConfig": self.life_cycle_config,
            "ThreadsPerCore": self.threads_per_core,
        }

    def summary(self) -> Dict[str, Any]:
        return {
            "InstanceGroupName": self.instance_group_name,
            "InstanceId": self.instance_id,
            "InstanceType": self.instance_type,
            "LaunchTime": self.launch_time,
            "InstanceStatus": {"Status": self.status, "Message": "message"},
        }


class CompilationJob(BaseObject):
    def __init__(
        self,
        compilation_job_name: str,
        role_arn: str,
        region_name: str,
        account_id: str,
        output_config: Dict[str, Any],
        stopping_condition: Dict[str, Any],
        model_package_version_arn: Optional[str],
        input_config: Optional[Dict[str, Any]],
        vpc_config: Optional[Dict[str, Any]],
        tags: Optional[List[Dict[str, str]]],
    ):
        self.compilation_job_name = compilation_job_name
        if (
            compilation_job_name
            in sagemaker_backends[account_id][region_name].compilation_jobs
        ):
            raise ResourceInUseException(
                message=f"Resource Already Exists: Compilation job with name {compilation_job_name} already exists. Choose a different name."
            )
        self.arn = arn_formatter(
            "compilation-job", self.compilation_job_name, account_id, region_name
        )
        self.compilation_job_status = "COMPLETED"
        self.compilation_start_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        self.compilation_end_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        self.stopping_condition = stopping_condition
        self.inference_image = "InferenceImage"
        self.model_package_version_arn = model_package_version_arn
        self.creation_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        self.last_modified_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        self.failure_reason = ""
        self.model_artifacts = {"S3ModelArtifacts": output_config["S3OutputLocation"]}
        self.model_digests = {
            "ArtifactDigest": "786a02f742015903c6c6fd852552d272912f4740e15847618a86e217f71f5419d25e1031afee585313896444934eb04b903a685b1448b755d56f701afe9be2ce"
        }
        self.role_arn = role_arn
        self.input_config = input_config
        if input_config and model_package_version_arn:
            raise ValidationError(
                message="InputConfig and ModelPackageVersionArn cannot be specified at the same time."
            )
        if not input_config and not model_package_version_arn:
            raise ValidationError(
                message="Either InputConfig or ModelPackageVersionArn must be specified."
            )
        self.output_config = output_config
        self.vpc_config = vpc_config
        self.derived_information = {"DerivedDataInputConfig": "DerivedDataInputConfig"}
        self.tags = tags

    def describe(self) -> Dict[str, Any]:
        return {
            "CompilationJobName": self.compilation_job_name,
            "CompilationJobArn": self.arn,
            "CompilationJobStatus": self.compilation_job_status,
            "CompilationStartTime": self.compilation_start_time,
            "CompilationEndTime": self.compilation_end_time,
            "StoppingCondition": self.stopping_condition,
            "InferenceImage": self.inference_image,
            "ModelPackageVersionArn": self.model_package_version_arn,
            "CreationTime": self.creation_time,
            "LastModifiedTime": self.last_modified_time,
            "FailureReason": self.failure_reason,
            "ModelArtifacts": self.model_artifacts,
            "ModelDigests": self.model_digests,
            "RoleArn": self.role_arn,
            "InputConfig": self.input_config,
            "OutputConfig": self.output_config,
            "VpcConfig": self.vpc_config,
            "DerivedInformation": self.derived_information,
        }

    def summary(self) -> Dict[str, Any]:
        summary = {
            "CompilationJobName": self.compilation_job_name,
            "CompilationJobArn": self.arn,
            "CreationTime": self.creation_time,
            "CompilationStartTime": self.compilation_start_time,
            "CompilationEndTime": self.compilation_end_time,
            "LastModifiedTime": self.last_modified_time,
            "CompilationJobStatus": self.compilation_job_status,
        }
        if "TargetDevice" in self.output_config:
            summary["CompilationTargetDevice"] = self.output_config["TargetDevice"]
        else:
            summary["CompilationTargetPlatformOs"] = self.output_config[
                "TargetPlatform"
            ]["Os"]
            summary["CompilationTargetPlatformArch"] = self.output_config[
                "TargetPlatform"
            ]["Arch"]
            summary["CompilationTargetPlatformAccelerator"] = self.output_config[
                "TargetPlatform"
            ]["Accelerator"]
        return summary


class AutoMLJob(BaseObject):
    def __init__(
        self,
        auto_ml_job_name: str,
        auto_ml_job_input_data_config: List[Dict[str, Any]],
        output_data_config: Dict[str, Any],
        auto_ml_problem_type_config: Dict[str, Any],
        role_arn: str,
        region_name: str,
        account_id: str,
        security_config: Optional[Dict[str, Any]],
        auto_ml_job_objective: Optional[Dict[str, Any]],
        model_deploy_config: Optional[Dict[str, Any]],
        data_split_config: Optional[Dict[str, Any]],
        tags: Optional[List[Dict[str, str]]] = None,
    ):
        self.region_name = region_name
        self.account_id = account_id
        self.auto_ml_job_name = auto_ml_job_name
        if auto_ml_job_name in sagemaker_backends[account_id][region_name].auto_ml_jobs:
            raise ResourceInUseException(
                message=f"Resource Already Exists: Auto ML Job with name {auto_ml_job_name} already exists. Choose a different name."
            )
        self.auto_ml_job_input_data_config = auto_ml_job_input_data_config
        self.output_data_config = output_data_config
        self.auto_ml_problem_type_config = auto_ml_problem_type_config
        self.role_arn = role_arn
        self.security_config = security_config
        self.auto_ml_job_objective = auto_ml_job_objective
        self.auto_ml_problem_type_resolved_attributes = {
            "SDK_UNKNOWN_MEMBER": {"name": "UnknownMemberName"}
        }
        if "ImageClassificationJobConfig" in self.auto_ml_problem_type_config:
            self.auto_ml_job_objective = (
                {"MetricName": "Accuracy"}
                if self.auto_ml_job_objective is None
                else self.auto_ml_job_objective
            )
            self.auto_ml_problem_type_config_name = "ImageClassification"
        elif "TextClassificationJobConfig" in self.auto_ml_problem_type_config:
            self.auto_ml_job_objective = (
                {"MetricName": "Accuracy"}
                if self.auto_ml_job_objective is None
                else self.auto_ml_job_objective
            )
            self.auto_ml_problem_type_config_name = "TextClassification"
        elif "TimeSeriesForecastingJobConfig" in self.auto_ml_problem_type_config:
            self.auto_ml_job_objective = (
                {"MetricName": "AverageWeightedQuantileLoss"}
                if self.auto_ml_job_objective is None
                else self.auto_ml_job_objective
            )
            self.auto_ml_problem_type_config_name = "TimeSeriesForecasting"
        elif "TabularJobConfig" in self.auto_ml_problem_type_config:
            self.auto_ml_problem_type_config_name = "Tabular"
            if (
                self.auto_ml_problem_type_config["TabularJobConfig"]["ProblemType"]
                == "BinaryClassification"
            ):
                self.auto_ml_job_objective = (
                    {"MetricName": "F1"}
                    if self.auto_ml_job_objective is None
                    else self.auto_ml_job_objective
                )
                self.auto_ml_problem_type_resolved_attributes = {
                    "TabularResolvedAttributes": {
                        "TabularProblemType": "BinaryClassification"
                    }
                }
            if (
                self.auto_ml_problem_type_config["TabularJobConfig"]["ProblemType"]
                == "MulticlassClassification"
            ):
                self.auto_ml_job_objective = (
                    {"MetricName": "Accuracy"}
                    if self.auto_ml_job_objective is None
                    else self.auto_ml_job_objective
                )
                self.auto_ml_problem_type_resolved_attributes = {
                    "TabularResolvedAttributes": {
                        "TabularProblemType": "MulticlassClassification"
                    }
                }
            if (
                self.auto_ml_problem_type_config["TabularJobConfig"]["ProblemType"]
                == "Regression"
            ):
                self.auto_ml_job_objective = (
                    {"MetricName": "MSE"}
                    if self.auto_ml_job_objective is None
                    else self.auto_ml_job_objective
                )
                self.auto_ml_problem_type_resolved_attributes = {
                    "TabularResolvedAttributes": {"TabularProblemType": "Regression"}
                }
        elif "TextGenerationJobConfig" in self.auto_ml_problem_type_config:
            self.auto_ml_problem_type_config_name = "TextGeneration"
            self.auto_ml_job_objective = (
                {"MetricName": ""}
                if self.auto_ml_job_objective is None
                else self.auto_ml_job_objective
            )
            self.auto_ml_problem_type_resolved_attributes = {
                "TextGenerationResolvedAttributes": {"BaseModelName": "string"}
            }

        self.model_deploy_config = (
            model_deploy_config
            if model_deploy_config
            else {"AutoGenerateEndpointName": False, "EndpointName": "EndpointName"}
        )
        if (
            "AutoGenerateEndpointName" in self.model_deploy_config
            and self.model_deploy_config["AutoGenerateEndpointName"]
            and "EndpointName" in self.model_deploy_config
        ):
            raise ValidationError(
                message="Validation Error: An EndpointName cannot be provided while AutoGenerateEndpoint name is True."
            )
        self.output_data_config = output_data_config
        self.data_split_config = (
            data_split_config if data_split_config else {"ValidationFraction": 0.2}
        )
        self.tags = tags or []
        self.arn = arn_formatter(
            "automl-job", self.auto_ml_job_name, account_id, region_name
        )
        self.creation_time = str(datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        self.end_time = str(datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        self.last_modified_time = str(datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        self.failure_reason = ""
        self.partial_failure_reasons = [{"PartialFailureMessage": ""}]
        self.best_candidate = {
            "CandidateName": "best_candidate",
            "FinalAutoMLJobObjectiveMetric": {
                "Type": "Maximize",
                "MetricName": "Accuracy",
                "Value": 123,
                "StandardMetricName": "Accuracy",
            },
            "ObjectiveStatus": "Succeeded",
            "CandidateSteps": [
                {
                    "CandidateStepType": "AWS::SageMaker::TrainingJob",
                    "CandidateStepArn": arn_formatter(
                        "training-job", "candidate_step_name", account_id, region_name
                    ),
                    "CandidateStepName": "candidate_step_name",
                },
            ],
            "CandidateStatus": "Completed",
            "InferenceContainers": [
                {
                    "Image": "string",
                    "ModelDataUrl": "string",
                    "Environment": {"string": "string"},
                },
            ],
            "CreationTime": str(datetime(2024, 1, 1)),
            "EndTime": str(datetime(2024, 1, 1)),
            "LastModifiedTime": str(datetime(2024, 1, 1)),
            "FailureReason": "string",
            "CandidateProperties": {
                "CandidateArtifactLocations": {
                    "Explainability": "string",
                    "ModelInsights": "string",
                    "BacktestResults": "string",
                },
                "CandidateMetrics": [
                    {
                        "MetricName": "Accuracy",
                        "Value": 123,
                        "Set": "Train",
                        "StandardMetricName": "Accuracy",
                    },
                ],
            },
            "InferenceContainerDefinitions": {
                "string": [
                    {
                        "Image": "string",
                        "ModelDataUrl": "string",
                        "Environment": {"string": "string"},
                    },
                ]
            },
        }
        self.auto_ml_job_status = "InProgress"
        self.auto_ml_job_secondary_status = "Completed"
        self.auto_ml_job_artifacts = {
            "CandidateDefinitionNotebookLocation": "candidate/notebook/location",
            "DataExplorationNotebookLocation": "data/notebook/location",
        }

        self.resolved_attributes = {
            "AutoMLJobObjective": self.auto_ml_job_objective,
            "CompletionCriteria": self.auto_ml_problem_type_config[
                self.auto_ml_problem_type_config_name + "JobConfig"
            ]["CompletionCriteria"],
            "AutoMLProblemTypeResolvedAttributes": self.auto_ml_problem_type_resolved_attributes,
        }

        self.model_deploy_result = {
            "EndpointName": self.model_deploy_config["EndpointName"]
            if self.model_deploy_config
            else "endpoint_name",
        }

    def describe(self) -> Dict[str, Any]:
        return {
            "AutoMLJobName": self.auto_ml_job_name,
            "AutoMLJobArn": self.arn,
            "AutoMLJobInputDataConfig": self.auto_ml_job_input_data_config,
            "OutputDataConfig": self.output_data_config,
            "RoleArn": self.role_arn,
            "AutoMLJobObjective": self.auto_ml_job_objective,
            "AutoMLProblemTypeConfig": self.auto_ml_problem_type_config,
            "AutoMLProblemTypeConfigName": self.auto_ml_problem_type_config_name,
            "CreationTime": self.creation_time,
            "EndTime": self.end_time,
            "LastModifiedTime": self.last_modified_time,
            "FailureReason": self.failure_reason,
            "PartialFailureReasons": self.partial_failure_reasons,
            "BestCandidate": self.best_candidate,
            "AutoMLJobStatus": self.auto_ml_job_status,
            "AutoMLJobSecondaryStatus": self.auto_ml_job_secondary_status,
            "AutoMLJobArtifacts": self.auto_ml_job_artifacts,
            "ResolvedAttributes": self.resolved_attributes,
            "ModelDeployConfig": self.model_deploy_config,
            "ModelDeployResult": self.model_deploy_result,
            "DataSplitConfig": self.data_split_config,
            "SecurityConfig": self.security_config,
        }

    def summary(self) -> Dict[str, Any]:
        return {
            "AutoMLJobName": self.auto_ml_job_name,
            "AutoMLJobArn": self.arn,
            "AutoMLJobStatus": self.auto_ml_job_status,
            "AutoMLJobSecondaryStatus": self.auto_ml_job_secondary_status,
            "CreationTime": self.creation_time,
            "EndTime": self.end_time,
            "LastModifiedTime": self.last_modified_time,
            "FailureReason": self.failure_reason,
            "PartialFailureReasons": self.partial_failure_reasons,
        }


class Domain(BaseObject):
    def __init__(
        self,
        domain_name: str,
        auth_mode: str,
        default_user_settings: Dict[str, Any],
        subnet_ids: List[str],
        vpc_id: str,
        account_id: str,
        region_name: str,
        domain_settings: Optional[Dict[str, Any]],
        tags: Optional[List[Dict[str, str]]],
        app_network_access_type: Optional[str],
        home_efs_file_system_kms_key_id: Optional[str],
        kms_key_id: Optional[str],
        app_security_group_management: Optional[str],
        default_space_settings: Optional[Dict[str, Any]],
    ):
        self.domain_name = domain_name
        if domain_name in sagemaker_backends[account_id][region_name].domains:
            raise ResourceInUseException(
                message=f"Resource Already Exists: Domain with name {domain_name} already exists. Choose a different name."
            )
        self.auth_mode = auth_mode
        self.default_user_settings = default_user_settings
        self.subnet_ids = subnet_ids
        self.vpc_id = vpc_id
        self.account_id = account_id
        self.region_name = region_name
        self.domain_settings = domain_settings
        self.tags = tags
        self.app_network_access_type = (
            app_network_access_type if app_network_access_type else "PublicInternetOnly"
        )
        self.home_efs_file_system_kms_key_id = (
            home_efs_file_system_kms_key_id
            if home_efs_file_system_kms_key_id
            else kms_key_id
        )
        self.kms_key_id = kms_key_id
        self.app_security_group_management = app_security_group_management
        self.default_space_settings = default_space_settings

        self.id = f"d-{domain_name}"
        self.arn = arn_formatter("domain", self.id, account_id, region_name)
        self.home_efs_file_system_id = f"{domain_name}-efs-id"
        self.single_sign_on_managed_application_instance_id = f"{domain_name}-sso-id"
        self.single_sign_on_managed_application_arn = arn_formatter(
            "sso", f"application/{domain_name}/apl-{domain_name}", account_id, ""
        )
        self.status = "InService"
        self.creation_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        self.last_modified_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        self.failure_reason = ""
        self.security_group_id_for_domain_boundary = f"sg-{domain_name}"
        self.url = f"{domain_name}.{region_name}.sagemaker.test.com"

    def describe(self) -> Dict[str, Any]:
        return {
            "DomainArn": self.arn,
            "DomainId": self.id,
            "DomainName": self.domain_name,
            "HomeEfsFileSystemId": self.home_efs_file_system_id,
            "SingleSignOnManagedApplicationInstanceId": self.single_sign_on_managed_application_instance_id,
            "SingleSignOnApplicationArn": self.single_sign_on_managed_application_arn,
            "Status": self.status,
            "CreationTime": self.creation_time,
            "LastModifiedTime": self.last_modified_time,
            "FailureReason": self.failure_reason,
            "SecurityGroupIdForDomainBoundary": self.security_group_id_for_domain_boundary,
            "AuthMode": self.auth_mode,
            "DefaultUserSettings": self.default_user_settings,
            "DomainSetting": self.domain_settings,
            "AppNetworkAccessType": self.app_network_access_type,
            "HomeEfsFileSystemKmsKeyId": self.home_efs_file_system_kms_key_id,
            "SubnetIds": self.subnet_ids,
            "Url": self.url,
            "VpcId": self.vpc_id,
            "KmsKeyId": self.kms_key_id,
            "AppSecurityGroupManagement": self.app_security_group_management,
            "DefaultSpaceSettings": self.default_space_settings,
        }

    def summary(self) -> Dict[str, Any]:
        return {
            "DomainArn": self.arn,
            "DomainId": self.id,
            "DomainName": self.domain_name,
            "Status": self.status,
            "CreationTime": self.creation_time,
            "LastModifiedTime": self.last_modified_time,
            "Url": self.url,
        }


class ModelExplainabilityJobDefinition(BaseObject):
    def __init__(
        self,
        job_definition_name: str,
        model_explainability_baseline_config: Optional[Dict[str, Any]],
        model_explainability_app_specification: Dict[str, Any],
        model_explainability_job_input: Dict[str, Any],
        model_explainability_job_output_config: Dict[str, Any],
        job_resources: Dict[str, Any],
        network_config: Optional[Dict[str, Any]],
        role_arn: str,
        stopping_condition: Optional[Dict[str, Any]],
        region_name: str,
        account_id: str,
        tags: Optional[List[Dict[str, str]]],
    ):
        self.job_definition_name = job_definition_name
        if (
            job_definition_name
            in sagemaker_backends[account_id][
                region_name
            ].model_explainability_job_definitions
        ):
            raise ResourceInUseException(
                message=f"Resource Already Exists: ModelExplainabilityJobDefinition with name {job_definition_name} already exists. Choose a different name."
            )
        self.model_explainability_baseline_config = model_explainability_baseline_config
        self.model_explainability_app_specification = (
            model_explainability_app_specification
        )
        self.model_explainability_job_input = model_explainability_job_input
        self.model_explainability_job_output_config = (
            model_explainability_job_output_config
        )
        self.job_resources = job_resources
        self.network_config = network_config
        self.role_arn = role_arn
        self.stopping_condition = stopping_condition
        self.region_name = region_name
        self.account_id = account_id
        self.tags = tags

        self.arn = arn_formatter(
            "model-explainability-job-definition",
            job_definition_name,
            self.account_id,
            self.region_name,
        )
        self.creation_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        self.endpoint_name = model_explainability_job_input["EndpointInput"][
            "EndpointName"
        ]

    def describe(self) -> Dict[str, Any]:
        return {
            "JobDefinitionArn": self.arn,
            "JobDefinitionName": self.job_definition_name,
            "CreationTime": self.creation_time,
            "ModelExplainabilityBaselineConfig": self.model_explainability_baseline_config,
            "ModelExplainabilityAppSpecification": self.model_explainability_app_specification,
            "ModelExplainabilityJobInput": self.model_explainability_job_input,
            "ModelExplainabilityJobOutputConfig": self.model_explainability_job_output_config,
            "JobResources": self.job_resources,
            "NetworkConfig": self.network_config,
            "RoleArn": self.role_arn,
            "StoppingConditions": self.stopping_condition,
        }

    def summary(self) -> Dict[str, Any]:
        return {
            "MonitoringJobDefinitionName": self.job_definition_name,
            "MonitoringJobDefinitionArn": self.arn,
            "CreationTime": self.creation_time,
            "EndpointName": self.endpoint_name,
        }


class HyperParameterTuningJob(BaseObject):
    def __init__(
        self,
        hyper_parameter_tuning_job_name: str,
        hyper_parameter_tuning_job_config: Dict[str, Any],
        region_name: str,
        account_id: str,
        training_job_definition: Optional[Dict[str, Any]],
        training_job_definitions: Optional[List[Dict[str, Any]]],
        warm_start_config: Optional[Dict[str, Any]],
        tags: Optional[List[Dict[str, str]]],
        autotune: Optional[Dict[str, Any]],
    ):
        self.hyper_parameter_tuning_job_name = hyper_parameter_tuning_job_name
        if (
            hyper_parameter_tuning_job_name
            in sagemaker_backends[account_id][region_name].hyper_parameter_tuning_jobs
        ):
            raise ResourceInUseException(
                message=f"Resource Already Exists: Hyper Parameter Tuning Job with name {hyper_parameter_tuning_job_name} already exists. Choose a different name."
            )
        self.arn = arn_formatter(
            "hyper-parameter-tuning-job",
            self.hyper_parameter_tuning_job_name,
            account_id,
            region_name,
        )
        self.hyper_parameter_tuning_job_config = hyper_parameter_tuning_job_config
        self.region_name = region_name
        self.account_id = account_id
        self.training_job_definition = training_job_definition
        self.training_job_definitions = training_job_definitions
        self.hyper_parameter_tuning_job_status = "Completed"
        self.creation_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        self.last_modified_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        self.hyper_parameter_tuning_end_time = datetime.now().strftime(
            "%Y-%m-%d %H:%M:%S"
        )
        self.training_job_status_counters = {
            "Completed": 1,
            "InProgress": 0,
            "RetryableError": 0,
            "NonRetryableError": 0,
            "Stopped": 0,
        }
        self.objective_status_counters = {
            "Succeeded": 1,
            "Pending": 0,
            "Failed": 0,
        }
        self.best_training_job = {
            "TrainingJobDefinitionName": "string",
            "TrainingJobName": "FakeTrainingJobName",
            "TrainingJobArn": "FakeTrainingJobArn",
            "TuningJobName": "FakeTuningJobName",
            "CreationTime": str(datetime(2024, 1, 1)),
            "TrainingStartTime": str(datetime(2024, 1, 1)),
            "TrainingEndTime": str(datetime(2024, 1, 1)),
            "TrainingJobStatus": "Completed",
            "TunedHyperParameters": {"string": "TunedHyperParameters"},
            "FailureReason": "string",
            "FinalHyperParameterTuningJobObjectiveMetric": {
                "Type": "Maximize",
                "MetricName": "Accuracy",
                "Value": 1,
            },
            "ObjectiveStatus": "Succeeded",
        }
        self.OverallBestTrainingJob = {
            "TrainingJobDefinitionName": "FakeTrainingJobDefinitionName",
            "TrainingJobName": "FakeTrainingJobName",
            "TrainingJobArn": "FakeTrainingJobArn",
            "TuningJobName": "FakeTuningJobName",
            "CreationTime": str(datetime(2024, 1, 1)),
            "TrainingStartTime": str(datetime(2024, 1, 1)),
            "TrainingEndTime": str(datetime(2024, 1, 1)),
            "TrainingJobStatus": "Completed",
            "TunedHyperParameters": {"string": "FakeTunedHyperParameters"},
            "FailureReason": "FakeFailureReason",
            "FinalHyperParameterTuningJobObjectiveMetric": {
                "Type": "Maximize",
                "MetricName": "Acccuracy",
                "Value": 1,
            },
            "ObjectiveStatus": "Succeeded",
        }
        self.warm_start_config = warm_start_config
        self.failure_reason = ""
        self.tuning_job_completion_details = {
            "NumberOfTrainingJobsObjectiveNotImproving": 123,
            "ConvergenceDetectedTime": str(datetime(2024, 1, 1)),
        }
        self.consumed_resources = {"RuntimeInSeconds": 123}
        self.tags = tags
        self.autotune = autotune

    def describe(self) -> Dict[str, Any]:
        return {
            "HyperParameterTuningJobName": self.hyper_parameter_tuning_job_name,
            "HyperParameterTuningJobArn": self.arn,
            "HyperParameterTuningJobConfig": self.hyper_parameter_tuning_job_config,
            "TrainingJobDefinition": self.training_job_definition,
            "TrainingJobDefinitions": self.training_job_definitions,
            "HyperParameterTuningJobStatus": self.hyper_parameter_tuning_job_status,
            "CreationTime": self.creation_time,
            "HyperParameterTuningEndTime": self.hyper_parameter_tuning_end_time,
            "LastModifiedTime": self.last_modified_time,
            "TrainingJobStatusCounters": self.training_job_status_counters,
            "ObjectiveStatusCounters": self.objective_status_counters,
            "BestTrainingJob": self.best_training_job,
            "OverallBestTrainingJob": self.OverallBestTrainingJob,
            "WarmStartConfig": self.warm_start_config,
            "Autotune": self.autotune,
            "FailureReason": self.failure_reason,
            "TuningJobCompletionDetails": self.tuning_job_completion_details,
            "ConsumedResources": self.consumed_resources,
        }

    def summary(self) -> Dict[str, Any]:
        return {
            "HyperParameterTuningJobName": self.hyper_parameter_tuning_job_name,
            "HyperParameterTuningJobArn": self.arn,
            "HyperParameterTuningJobStatus": self.hyper_parameter_tuning_job_status,
            "Strategy": self.hyper_parameter_tuning_job_config["Strategy"],
            "CreationTime": self.creation_time,
            "HyperParameterTuningEndTime": self.hyper_parameter_tuning_end_time,
            "LastModifiedTime": self.last_modified_time,
            "TrainingJobStatusCounters": self.training_job_status_counters,
            "ObjectiveStatusCounters": self.objective_status_counters,
            "ResourceLimits": self.hyper_parameter_tuning_job_config["ResourceLimits"],
        }


class ModelQualityJobDefinition(BaseObject):
    def __init__(
        self,
        job_definition_name: str,
        model_quality_baseline_config: Optional[Dict[str, Any]],
        model_quality_app_specification: Dict[str, Any],
        model_quality_job_input: Dict[str, Any],
        model_quality_job_output_config: Dict[str, Any],
        job_resources: Dict[str, Any],
        network_config: Optional[Dict[str, Any]],
        role_arn: str,
        stopping_condition: Optional[Dict[str, Any]],
        tags: Optional[List[Dict[str, str]]],
        region_name: str,
        account_id: str,
    ):
        self.region_name = region_name
        self.account_id = account_id
        self.job_definition_name = job_definition_name
        if (
            job_definition_name
            in sagemaker_backends[account_id][region_name].model_quality_job_definitions
        ):
            raise ResourceInUseException(
                message=f"Resource Already Exists: Model Quality Job Definition with name {job_definition_name} already exists. Choose a different name."
            )
        self.model_quality_baseline_config = model_quality_baseline_config
        self.model_quality_app_specification = model_quality_app_specification
        self.model_quality_job_input = model_quality_job_input
        self.model_quality_job_output_config = model_quality_job_output_config
        self.job_resources = job_resources
        self.network_config = network_config
        self.role_arn = role_arn
        self.stopping_condition = stopping_condition
        self.tags = tags or []
        self.arn = arn_formatter(
            "model-quality-job-definition",
            self.job_definition_name,
            account_id,
            region_name,
        )
        self.creation_time = str(datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        self.endpoint_name = self.model_quality_job_input["EndpointInput"][
            "EndpointName"
        ]

    def describe(self) -> Dict[str, Any]:
        return {
            "JobDefinitionArn": self.arn,
            "JobDefinitionName": self.job_definition_name,
            "CreationTime": self.creation_time,
            "ModelQualityBaselineConfig": self.model_quality_baseline_config,
            "ModelQualityAppSpecification": self.model_quality_app_specification,
            "ModelQualityJobInput": self.model_quality_job_input,
            "ModelQualityJobOutputConfig": self.model_quality_job_output_config,
            "JobResources": self.job_resources,
            "NetworkConfig": self.network_config,
            "RoleArn": self.role_arn,
            "StoppingCondition": self.stopping_condition,
        }

    def summary(self) -> Dict[str, Any]:
        return {
            "MonitoringJobDefinitionName": self.job_definition_name,
            "MonitoringJobDefinitionArn": self.arn,
            "CreationTime": self.creation_time,
            "EndpointName": self.endpoint_name,
        }


class VpcConfig(BaseObject):
    def __init__(self, security_group_ids: List[str], subnets: List[str]):
        self.security_group_ids = security_group_ids
        self.subnets = subnets

    @property
    def response_object(self) -> Dict[str, List[str]]:
        response_object = self.gen_response_object()
        return {
            k: v for k, v in response_object.items() if v is not None and v != [None]
        }


class Container(BaseObject):
    def __init__(self, **kwargs: Any):
        self.container_hostname = kwargs.get("container_hostname", "localhost")
        self.model_data_url = kwargs.get("data_url", "")
        self.model_package_name = kwargs.get("package_name", "pkg")
        self.image = kwargs.get("image", "")
        self.environment = kwargs.get("environment", {})

    @property
    def response_object(self) -> Dict[str, Any]:  # type: ignore[misc]
        response_object = self.gen_response_object()
        return {
            k: v for k, v in response_object.items() if v is not None and v != [None]
        }


class FakeSagemakerNotebookInstance(CloudFormationModel):
    def __init__(
        self,
        account_id: str,
        region_name: str,
        notebook_instance_name: str,
        instance_type: str,
        role_arn: str,
        subnet_id: Optional[str],
        security_group_ids: Optional[List[str]],
        kms_key_id: Optional[str],
        tags: Optional[List[Dict[str, str]]],
        lifecycle_config_name: Optional[str],
        direct_internet_access: str,
        volume_size_in_gb: int,
        accelerator_types: Optional[List[str]],
        default_code_repository: Optional[str],
        additional_code_repositories: Optional[List[str]],
        root_access: Optional[str],
    ):
        self.validate_volume_size_in_gb(volume_size_in_gb)
        self.validate_instance_type(instance_type)

        self.region_name = region_name
        self.notebook_instance_name = notebook_instance_name
        self.instance_type = instance_type
        self.role_arn = role_arn
        self.subnet_id = subnet_id
        self.security_group_ids = security_group_ids
        self.kms_key_id = kms_key_id
        self.tags = tags or []
        self.lifecycle_config_name = lifecycle_config_name
        self.direct_internet_access = direct_internet_access
        self.volume_size_in_gb = volume_size_in_gb
        self.accelerator_types = accelerator_types
        self.default_code_repository = default_code_repository
        self.additional_code_repositories = additional_code_repositories
        self.root_access = root_access
        self.status = "Pending"
        self.creation_time = self.last_modified_time = datetime.now()
        self.arn = arn_formatter(
            "notebook-instance", notebook_instance_name, account_id, region_name
        )
        self.start()

    def validate_volume_size_in_gb(self, volume_size_in_gb: int) -> None:
        if not validators.is_integer_between(volume_size_in_gb, mn=5, optional=True):
            message = "Invalid range for parameter VolumeSizeInGB, value: {}, valid range: 5-inf"
            raise ValidationError(message=message)

    def validate_instance_type(self, instance_type: str) -> None:
        VALID_INSTANCE_TYPES = [
            "ml.p2.xlarge",
            "ml.m5.4xlarge",
            "ml.m4.16xlarge",
            "ml.t3.xlarge",
            "ml.p3.16xlarge",
            "ml.t2.xlarge",
            "ml.p2.16xlarge",
            "ml.c4.2xlarge",
            "ml.c5.2xlarge",
            "ml.c4.4xlarge",
            "ml.c5d.2xlarge",
            "ml.c5.4xlarge",
            "ml.c5d.4xlarge",
            "ml.c4.8xlarge",
            "ml.c5d.xlarge",
            "ml.c5.9xlarge",
            "ml.c5.xlarge",
            "ml.c5d.9xlarge",
            "ml.c4.xlarge",
            "ml.t2.2xlarge",
            "ml.c5d.18xlarge",
            "ml.t3.2xlarge",
            "ml.t3.medium",
            "ml.t2.medium",
            "ml.c5.18xlarge",
            "ml.p3.2xlarge",
            "ml.m5.xlarge",
            "ml.m4.10xlarge",
            "ml.t2.large",
            "ml.m5.12xlarge",
            "ml.m4.xlarge",
            "ml.t3.large",
            "ml.m5.24xlarge",
            "ml.m4.2xlarge",
            "ml.p2.8xlarge",
            "ml.m5.2xlarge",
            "ml.p3.8xlarge",
            "ml.m4.4xlarge",
        ]
        if not validators.is_one_of(instance_type, VALID_INSTANCE_TYPES):
            message = f"Value '{instance_type}' at 'instanceType' failed to satisfy constraint: Member must satisfy enum value set: {VALID_INSTANCE_TYPES}"
            raise ValidationError(message=message)

    @property
    def url(self) -> str:
        return (
            f"{self.notebook_instance_name}.notebook.{self.region_name}.sagemaker.aws"
        )

    def start(self) -> None:
        self.status = "InService"

    @property
    def is_deletable(self) -> bool:
        return self.status in ["Stopped", "Failed"]

    def stop(self) -> None:
        self.status = "Stopped"

    @property
    def physical_resource_id(self) -> str:
        return self.arn

    @classmethod
    def has_cfn_attr(cls, attr: str) -> bool:
        return attr in ["NotebookInstanceName"]

    def get_cfn_attribute(self, attribute_name: str) -> str:
        # https://docs.aws.amazon.com/AWSCloudFormation/latest/UserGuide/aws-resource-sagemaker-notebookinstance.html#aws-resource-sagemaker-notebookinstance-return-values
        from moto.cloudformation.exceptions import UnformattedGetAttTemplateException

        if attribute_name == "NotebookInstanceName":
            return self.notebook_instance_name
        raise UnformattedGetAttTemplateException()

    @staticmethod
    def cloudformation_name_type() -> str:
        return ""

    @staticmethod
    def cloudformation_type() -> str:
        # https://docs.aws.amazon.com/AWSCloudFormation/latest/UserGuide/aws-resource-sagemaker-notebookinstance.html
        return "AWS::SageMaker::NotebookInstance"

    @classmethod
    def create_from_cloudformation_json(  # type: ignore[misc]
        cls,
        resource_name: str,
        cloudformation_json: Any,
        account_id: str,
        region_name: str,
        **kwargs: Any,
    ) -> "FakeSagemakerNotebookInstance":
        # Get required properties from provided CloudFormation template
        properties = cloudformation_json["Properties"]
        instance_type = properties["InstanceType"]
        role_arn = properties["RoleArn"]

        notebook = sagemaker_backends[account_id][region_name].create_notebook_instance(
            notebook_instance_name=resource_name,
            instance_type=instance_type,
            role_arn=role_arn,
        )
        return notebook

    @classmethod
    def update_from_cloudformation_json(  # type: ignore[misc]
        cls,
        original_resource: Any,
        new_resource_name: str,
        cloudformation_json: Any,
        account_id: str,
        region_name: str,
    ) -> "FakeSagemakerNotebookInstance":
        # Operations keep same resource name so delete old and create new to mimic update
        cls.delete_from_cloudformation_json(
            original_resource.arn, cloudformation_json, account_id, region_name
        )
        new_resource = cls.create_from_cloudformation_json(
            original_resource.notebook_instance_name,
            cloudformation_json,
            account_id,
            region_name,
        )
        return new_resource

    @classmethod
    def delete_from_cloudformation_json(  # type: ignore[misc]
        cls,
        resource_name: str,
        cloudformation_json: Any,
        account_id: str,
        region_name: str,
    ) -> None:
        # Get actual name because resource_name actually provides the ARN
        # since the Physical Resource ID is the ARN despite SageMaker
        # using the name for most of its operations.
        notebook_instance_name = resource_name.split("/")[-1]

        backend = sagemaker_backends[account_id][region_name]
        backend.stop_notebook_instance(notebook_instance_name)
        backend.delete_notebook_instance(notebook_instance_name)

    def to_dict(self) -> Dict[str, Any]:
        return {
            "NotebookInstanceArn": self.arn,
            "NotebookInstanceName": self.notebook_instance_name,
            "NotebookInstanceStatus": self.status,
            "Url": self.url,
            "InstanceType": self.instance_type,
            "SubnetId": self.subnet_id,
            "SecurityGroups": self.security_group_ids,
            "RoleArn": self.role_arn,
            "KmsKeyId": self.kms_key_id,
            # ToDo: NetworkInterfaceId
            "LastModifiedTime": str(self.last_modified_time),
            "CreationTime": str(self.creation_time),
            "NotebookInstanceLifecycleConfigName": self.lifecycle_config_name,
            "DirectInternetAccess": self.direct_internet_access,
            "VolumeSizeInGB": self.volume_size_in_gb,
            "AcceleratorTypes": self.accelerator_types,
            "DefaultCodeRepository": self.default_code_repository,
            "AdditionalCodeRepositories": self.additional_code_repositories,
            "RootAccess": self.root_access,
        }


class FakeSageMakerNotebookInstanceLifecycleConfig(BaseObject, CloudFormationModel):
    def __init__(
        self,
        account_id: str,
        region_name: str,
        notebook_instance_lifecycle_config_name: str,
        on_create: List[Dict[str, str]],
        on_start: List[Dict[str, str]],
    ):
        self.region_name = region_name
        self.notebook_instance_lifecycle_config_name = (
            notebook_instance_lifecycle_config_name
        )
        self.on_create = on_create
        self.on_start = on_start
        self.creation_time = self.last_modified_time = datetime.now().strftime(
            "%Y-%m-%d %H:%M:%S"
        )
        self.arn = FakeSageMakerNotebookInstanceLifecycleConfig.arn_formatter(
            self.notebook_instance_lifecycle_config_name, account_id, region_name
        )

    @staticmethod
    def arn_formatter(name: str, account_id: str, region_name: str) -> str:
        return arn_formatter(
            "notebook-instance-lifecycle-config", name, account_id, region_name
        )

    @property
    def response_object(self) -> Dict[str, Any]:  # type: ignore[misc]
        response_object = self.gen_response_object()
        response = {
            k: v for k, v in response_object.items() if v is not None and v != [None]
        }
        response["NotebookInstanceLifecycleConfigArn"] = response.pop("Arn")

        return response

    @property
    def physical_resource_id(self) -> str:
        return self.arn

    @classmethod
    def has_cfn_attr(cls, attr: str) -> bool:
        return attr in ["NotebookInstanceLifecycleConfigName"]

    def get_cfn_attribute(self, attribute_name: str) -> str:
        # https://docs.aws.amazon.com/AWSCloudFormation/latest/UserGuide/aws-resource-sagemaker-notebookinstancelifecycleconfig.html#aws-resource-sagemaker-notebookinstancelifecycleconfig-return-values
        from moto.cloudformation.exceptions import UnformattedGetAttTemplateException

        if attribute_name == "NotebookInstanceLifecycleConfigName":
            return self.notebook_instance_lifecycle_config_name
        raise UnformattedGetAttTemplateException()

    @staticmethod
    def cloudformation_name_type() -> str:
        return ""

    @staticmethod
    def cloudformation_type() -> str:
        # https://docs.aws.amazon.com/AWSCloudFormation/latest/UserGuide/aws-resource-sagemaker-notebookinstancelifecycleconfig.html
        return "AWS::SageMaker::NotebookInstanceLifecycleConfig"

    @classmethod
    def create_from_cloudformation_json(  # type: ignore[misc]
        cls,
        resource_name: str,
        cloudformation_json: Any,
        account_id: str,
        region_name: str,
        **kwargs: Any,
    ) -> "FakeSageMakerNotebookInstanceLifecycleConfig":
        properties = cloudformation_json["Properties"]

        config = sagemaker_backends[account_id][
            region_name
        ].create_notebook_instance_lifecycle_config(
            notebook_instance_lifecycle_config_name=resource_name,
            on_create=properties.get("OnCreate"),
            on_start=properties.get("OnStart"),
        )
        return config

    @classmethod
    def update_from_cloudformation_json(  # type: ignore[misc]
        cls,
        original_resource: Any,
        new_resource_name: str,
        cloudformation_json: Any,
        account_id: str,
        region_name: str,
    ) -> "FakeSageMakerNotebookInstanceLifecycleConfig":
        # Operations keep same resource name so delete old and create new to mimic update
        cls.delete_from_cloudformation_json(
            original_resource.arn,
            cloudformation_json,
            account_id,
            region_name,
        )
        new_resource = cls.create_from_cloudformation_json(
            original_resource.notebook_instance_lifecycle_config_name,
            cloudformation_json,
            account_id,
            region_name,
        )
        return new_resource

    @classmethod
    def delete_from_cloudformation_json(  # type: ignore[misc]
        cls,
        resource_name: str,
        cloudformation_json: Any,
        account_id: str,
        region_name: str,
    ) -> None:
        # Get actual name because resource_name actually provides the ARN
        # since the Physical Resource ID is the ARN despite SageMaker
        # using the name for most of its operations.
        config_name = resource_name.split("/")[-1]

        backend = sagemaker_backends[account_id][region_name]
        backend.delete_notebook_instance_lifecycle_config(config_name)


class SageMakerModelBackend(BaseBackend):
    def __init__(self, region_name: str, account_id: str):
        super().__init__(region_name, account_id)
        self._models: Dict[str, Model] = {}
        self.notebook_instances: Dict[str, FakeSagemakerNotebookInstance] = {}
        self.endpoint_configs: Dict[str, FakeEndpointConfig] = {}
        self.endpoints: Dict[str, FakeEndpoint] = {}
        self.experiments: Dict[str, FakeExperiment] = {}
        self.pipelines: Dict[str, FakePipeline] = {}
        self.pipeline_executions: Dict[str, FakePipelineExecution] = {}
        self.processing_jobs: Dict[str, FakeProcessingJob] = {}
        self.trials: Dict[str, FakeTrial] = {}
        self.trial_components: Dict[str, FakeTrialComponent] = {}
        self.training_jobs: Dict[str, FakeTrainingJob] = {}
        self.transform_jobs: Dict[str, FakeTransformJob] = {}
        self.notebook_instance_lifecycle_configurations: Dict[
            str, FakeSageMakerNotebookInstanceLifecycleConfig
        ] = {}
        self.model_cards: DefaultDict[str, List[FakeModelCard]] = defaultdict(list)
        self.model_package_groups: Dict[str, ModelPackageGroup] = {}
        self.model_packages: Dict[str, ModelPackage] = {}
        self.model_package_name_mapping: Dict[str, str] = {}
        self.feature_groups: Dict[str, FeatureGroup] = {}
        self.clusters: Dict[str, Cluster] = {}
        self.data_quality_job_definitions: Dict[str, FakeDataQualityJobDefinition] = {}
        self.model_bias_job_definitions: Dict[str, FakeModelBiasJobDefinition] = {}
        self.auto_ml_jobs: Dict[str, AutoMLJob] = {}
        self.compilation_jobs: Dict[str, CompilationJob] = {}
        self.domains: Dict[str, Domain] = {}
        self.model_explainability_job_definitions: Dict[
            str, ModelExplainabilityJobDefinition
        ] = {}
        self.hyper_parameter_tuning_jobs: Dict[str, HyperParameterTuningJob] = {}
        self.model_quality_job_definitions: Dict[str, ModelQualityJobDefinition] = {}

    @staticmethod
    def default_vpc_endpoint_service(
        service_region: str, zones: List[str]
    ) -> List[Dict[str, str]]:
        """Default VPC endpoint services."""
        api_service = BaseBackend.default_vpc_endpoint_service_factory(
            service_region, zones, "api.sagemaker", special_service_name="sagemaker.api"
        )

        notebook_service_id = f"vpce-svc-{BaseBackend.vpce_random_number()}"
        studio_service_id = f"vpce-svc-{BaseBackend.vpce_random_number()}"

        notebook_service = {
            "AcceptanceRequired": False,
            "AvailabilityZones": zones,
            "BaseEndpointDnsNames": [
                f"{notebook_service_id}.{service_region}.vpce.amazonaws.com",
                f"notebook.{service_region}.vpce.sagemaker.aws",
            ],
            "ManagesVpcEndpoints": False,
            "Owner": "amazon",
            "PrivateDnsName": f"*.notebook.{service_region}.sagemaker.aws",
            "PrivateDnsNameVerificationState": "verified",
            "PrivateDnsNames": [
                {"PrivateDnsName": f"*.notebook.{service_region}.sagemaker.aws"}
            ],
            "ServiceId": notebook_service_id,
            "ServiceName": f"aws.sagemaker.{service_region}.notebook",
            "ServiceType": [{"ServiceType": "Interface"}],
            "Tags": [],
            "VpcEndpointPolicySupported": True,
        }
        studio_service = {
            "AcceptanceRequired": False,
            "AvailabilityZones": zones,
            "BaseEndpointDnsNames": [
                f"{studio_service_id}.{service_region}.vpce.amazonaws.com",
                f"studio.{service_region}.vpce.sagemaker.aws",
            ],
            "ManagesVpcEndpoints": False,
            "Owner": "amazon",
            "PrivateDnsName": f"*.studio.{service_region}.sagemaker.aws",
            "PrivateDnsNameVerificationState": "verified",
            "PrivateDnsNames": [
                {"PrivateDnsName": f"*.studio.{service_region}.sagemaker.aws"}
            ],
            "ServiceId": studio_service_id,
            "ServiceName": f"aws.sagemaker.{service_region}.studio",
            "ServiceType": [{"ServiceType": "Interface"}],
            "Tags": [],
            "VpcEndpointPolicySupported": True,
        }
        return api_service + [notebook_service, studio_service]

    def create_model(
        self,
        model_name: str,
        execution_role_arn: str,
        primary_container: Optional[Dict[str, Any]],
        vpc_config: Optional[Dict[str, Any]],
        containers: Optional[List[Dict[str, Any]]],
        tags: Optional[List[Dict[str, str]]],
    ) -> Model:
        model_obj = Model(
            account_id=self.account_id,
            region_name=self.region_name,
            model_name=model_name,
            execution_role_arn=execution_role_arn,
            primary_container=primary_container or {},
            vpc_config=vpc_config or {},
            containers=containers or [],
            tags=tags or [],
        )

        self._models[model_name] = model_obj
        return model_obj

    def describe_model(self, model_name: str) -> Model:
        model = self._models.get(model_name)
        if model:
            return model
        arn = arn_formatter("model", model_name, self.account_id, self.region_name)
        raise ValidationError(message=f"Could not find model '{arn}'.")

    def list_models(self) -> Iterable[Model]:
        return self._models.values()

    def delete_model(self, model_name: str) -> None:
        for model in self._models.values():
            if model.model_name == model_name:
                self._models.pop(model.model_name)
                break
        else:
            raise MissingModel(model=model_name)

    def create_experiment(self, experiment_name: str) -> Dict[str, str]:
        experiment = FakeExperiment(
            account_id=self.account_id,
            region_name=self.region_name,
            experiment_name=experiment_name,
            tags=[],
        )
        self.experiments[experiment_name] = experiment
        return experiment.response_create

    def describe_experiment(self, experiment_name: str) -> Dict[str, Any]:
        experiment_data = self.experiments[experiment_name]
        return {
            "ExperimentName": experiment_data.experiment_name,
            "ExperimentArn": experiment_data.arn,
            "CreationTime": experiment_data.creation_time,
            "LastModifiedTime": experiment_data.last_modified_time,
        }

    def _get_resource_from_arn(self, arn: str) -> Any:
        resources = {
            "model": self._models,
            "notebook-instance": self.notebook_instances,
            "endpoint": self.endpoints,
            "endpoint-config": self.endpoint_configs,
            "training-job": self.training_jobs,
            "transform-job": self.transform_jobs,
            "experiment": self.experiments,
            "experiment-trial": self.trials,
            "experiment-trial-component": self.trial_components,
            "processing-job": self.processing_jobs,
            "pipeline": self.pipelines,
            "model-package-group": self.model_package_groups,
            "cluster": self.clusters,
            "data-quality-job-definition": self.data_quality_job_definitions,
            "model-bias-job-definition": self.model_bias_job_definitions,
            "automl-job": self.auto_ml_jobs,
            "compilation-job": self.compilation_jobs,
            "domain": self.domains,
            "model-explainability-job-definition": self.model_explainability_job_definitions,
            "hyper-parameter-tuning-job": self.hyper_parameter_tuning_jobs,
            "model-quality-job-definition": self.model_quality_job_definitions,
            "model-card": self.model_cards,
        }
        target_resource, target_name = arn.split(":")[-1].split("/")
        try:
            resource = resources.get(target_resource).get(target_name)  # type: ignore
        except KeyError:
            message = f"Could not find {target_resource} with name {target_name}"
            raise ValidationError(message=message)
        if isinstance(resource, list):
            return resource[0]
        return resource

    def add_tags(self, arn: str, tags: List[Dict[str, str]]) -> List[Dict[str, str]]:
        resource = self._get_resource_from_arn(arn)
        resource.tags.extend(tags)
        return resource.tags

    @paginate(pagination_model=PAGINATION_MODEL)
    def list_tags(self, arn: str) -> List[Dict[str, str]]:
        resource = self._get_resource_from_arn(arn)
        return resource.tags

    def delete_tags(self, arn: str, tag_keys: List[str]) -> None:
        resource = self._get_resource_from_arn(arn)
        resource.tags = [tag for tag in resource.tags if tag["Key"] not in tag_keys]

    @paginate(pagination_model=PAGINATION_MODEL)
    def list_experiments(self) -> List["FakeExperiment"]:
        return list(self.experiments.values())

    def search(self, resource: Any = None, search_expression: Any = None) -> Any:
        """
        Only a few SearchExpressions are implemented. Please open a bug report if you find any issues.
        """
        next_index = None

        valid_resources = {
            "Pipeline": self.pipelines.values(),
            "ModelPackageGroup": self.model_package_groups.values(),
            "TrainingJob": self.training_jobs.values(),
            "ExperimentTrialComponent": self.trial_components.values(),
            "FeatureGroup": self.feature_groups.values(),
            "Endpoint": self.endpoints.values(),
            "PipelineExecution": self.pipeline_executions.values(),
            "Project": [],
            "ExperimentTrial": self.trials.values(),
            "Image": [],
            "ImageVersion": [],
            "ModelPackage": self.model_packages.values(),
            "Experiment": self.experiments.values(),
        }

        if resource not in valid_resources:
            raise AWSValidationException(
                f"An error occurred (ValidationException) when calling the Search operation: 1 validation error detected: Value '{resource}' at 'resource' failed to satisfy constraint: Member must satisfy enum value set: {valid_resources}"
            )

        def compare_value(actual: Any, expected: Any, operator: str) -> bool:
            # Defeault:  operator == "Equals"
            if operator == "Contains":
                return expected in actual
            if operator == "NotEquals":
                return expected != actual
            return actual == expected

        def evaluate_search_expression(item: Any) -> bool:
            filters = None
            if search_expression is not None:
                filters = search_expression.get("Filters")

            if filters is not None:
                for f in filters:
                    prop_key = camelcase_to_underscores(f["Name"])
                    if f["Name"].startswith("Tags."):
                        key = f["Name"][5:]
                        value = f["Value"]

                        if f["Operator"] == "Equals":
                            if not [
                                e
                                for e in item.tags
                                if e["Key"] == key and e["Value"] == value
                            ]:
                                return False
                        return True

                    elif f["Name"] == "TrialName":
                        raise AWSValidationException(
                            f"An error occurred (ValidationException) when calling the Search operation: Unknown property name: {f['Name']}"
                        )

                    elif f["Name"] == "Parents.TrialName":
                        trial_name = f["Value"]
                        if getattr(item, "trial_name") != trial_name:
                            return False

                    elif hasattr(item, prop_key):
                        if not compare_value(
                            getattr(item, prop_key), f["Value"], f["Operator"]
                        ):
                            return False
                    else:
                        raise ValidationError(
                            message=f"Unknown property name: {f['Name']}"
                        )

            return True

        result: Dict[str, Any] = {
            "Results": [],
            "NextToken": str(next_index) if next_index is not None else None,
        }
        # ResourceName, ResultName, Resources
        result_names = {
            "ExperimentTrial": "Trial",
            "ExperimentTrialComponent": "TrialComponent",
        }
        resources_found = [
            x for x in valid_resources[resource] if evaluate_search_expression(x)
        ]
        result_name = result_names.get(resource, resource)
        for found in resources_found:
            result["Results"].append({result_name: found.gen_response_object()})

        return result

    def delete_experiment(self, experiment_name: str) -> None:
        try:
            del self.experiments[experiment_name]
        except KeyError:
            arn = FakeTrial.arn_formatter(
                experiment_name, self.account_id, self.region_name
            )
            raise ValidationError(
                message=f"Could not find experiment configuration '{arn}'."
            )

    def create_trial(self, trial_name: str, experiment_name: str) -> Dict[str, str]:
        trial = FakeTrial(
            account_id=self.account_id,
            region_name=self.region_name,
            trial_name=trial_name,
            experiment_name=experiment_name,
            tags=[],
            trial_components=[],
        )
        self.trials[trial_name] = trial
        return trial.response_create

    def describe_trial(self, trial_name: str) -> Dict[str, Any]:
        try:
            return self.trials[trial_name].response_object
        except KeyError:
            arn = FakeTrial.arn_formatter(trial_name, self.account_id, self.region_name)
            raise ValidationError(message=f"Could not find trial '{arn}'.")

    def delete_trial(self, trial_name: str) -> None:
        try:
            del self.trials[trial_name]
        except KeyError:
            arn = FakeTrial.arn_formatter(trial_name, self.account_id, self.region_name)
            raise ValidationError(
                message=f"Could not find trial configuration '{arn}'."
            )

    @paginate(pagination_model=PAGINATION_MODEL)
    def list_trials(
        self,
        experiment_name: Optional[str] = None,
        trial_component_name: Optional[str] = None,
    ) -> List["FakeTrial"]:
        trials_fetched = list(self.trials.values())

        def evaluate_filter_expression(trial_data: FakeTrial) -> bool:
            if experiment_name is not None:
                if trial_data.experiment_name != experiment_name:
                    return False

            if trial_component_name is not None:
                if trial_component_name not in trial_data.trial_components:
                    return False

            return True

        return [
            trial_data
            for trial_data in trials_fetched
            if evaluate_filter_expression(trial_data)
        ]

    def create_trial_component(
        self,
        trial_component_name: str,
        trial_name: str,
        status: Dict[str, str],
        start_time: Optional[datetime],
        end_time: Optional[datetime],
        display_name: Optional[str],
        parameters: Optional[Dict[str, Dict[str, Union[str, float]]]],
        input_artifacts: Optional[Dict[str, Dict[str, str]]],
        output_artifacts: Optional[Dict[str, Dict[str, str]]],
        metadata_properties: Optional[Dict[str, str]],
    ) -> Dict[str, Any]:
        trial_component = FakeTrialComponent(
            account_id=self.account_id,
            region_name=self.region_name,
            display_name=display_name,
            start_time=start_time,
            end_time=end_time,
            parameters=parameters,
            input_artifacts=input_artifacts,
            output_artifacts=output_artifacts,
            metadata_properties=metadata_properties,
            trial_component_name=trial_component_name,
            trial_name=trial_name,
            status=status,
            tags=[],
        )
        self.trial_components[trial_component_name] = trial_component
        return trial_component.response_create

    def delete_trial_component(self, trial_component_name: str) -> None:
        try:
            del self.trial_components[trial_component_name]
        except KeyError:
            arn = FakeTrial.arn_formatter(
                trial_component_name, self.account_id, self.region_name
            )
            raise ValidationError(
                message=f"Could not find trial-component configuration '{arn}'."
            )

    def describe_trial_component(self, trial_component_name: str) -> Dict[str, Any]:
        try:
            return self.trial_components[trial_component_name].response_object
        except KeyError:
            arn = FakeTrialComponent.arn_formatter(
                trial_component_name, self.account_id, self.region_name
            )
            raise ValidationError(message=f"Could not find trial component '{arn}'.")

    def _update_trial_component_details(
        self, trial_component_name: str, details_json: str
    ) -> None:
        self.trial_components[trial_component_name].update(details_json)

    @paginate(pagination_model=PAGINATION_MODEL)
    def list_trial_components(
        self, trial_name: Optional[str] = None
    ) -> List["FakeTrialComponent"]:
        trial_components_fetched = list(self.trial_components.values())

        return [
            trial_component_data
            for trial_component_data in trial_components_fetched
            if trial_name is None or trial_component_data.trial_name == trial_name
        ]

    def associate_trial_component(
        self, trial_name: str, trial_component_name: str
    ) -> Dict[str, str]:
        if trial_name in self.trials.keys():
            self.trials[trial_name].trial_components.extend([trial_component_name])
        else:
            raise ResourceNotFound(
                message=f"Trial 'arn:{get_partition(self.region_name)}:sagemaker:{self.region_name}:{self.account_id}:experiment-trial/{trial_name}' does not exist."
            )

        if trial_component_name in self.trial_components.keys():
            self.trial_components[trial_component_name].trial_name = trial_name

        return {
            "TrialComponentArn": self.trial_components[trial_component_name].arn,
            "TrialArn": self.trials[trial_name].arn,
        }

    def disassociate_trial_component(
        self, trial_name: str, trial_component_name: str
    ) -> Dict[str, str]:
        if trial_component_name in self.trial_components.keys():
            self.trial_components[trial_component_name].trial_name = None

        if trial_name in self.trials.keys():
            self.trials[trial_name].trial_components = list(
                filter(
                    lambda x: x != trial_component_name,
                    self.trials[trial_name].trial_components,
                )
            )

        return {
            "TrialComponentArn": f"arn:{get_partition(self.region_name)}:sagemaker:{self.region_name}:{self.account_id}:experiment-trial-component/{trial_component_name}",
            "TrialArn": f"arn:{get_partition(self.region_name)}:sagemaker:{self.region_name}:{self.account_id}:experiment-trial/{trial_name}",
        }

    def update_trial_component(
        self,
        trial_component_name: str,
        status: Optional[Dict[str, str]],
        display_name: Optional[str],
        start_time: Optional[datetime],
        end_time: Optional[datetime],
        parameters: Optional[Dict[str, Dict[str, Union[str, float]]]],
        parameters_to_remove: Optional[List[str]],
        input_artifacts: Optional[Dict[str, Dict[str, str]]],
        input_artifacts_to_remove: Optional[List[str]],
        output_artifacts: Optional[Dict[str, Dict[str, str]]],
        output_artifacts_to_remove: Optional[List[str]],
    ) -> Dict[str, str]:
        try:
            trial_component = self.trial_components[trial_component_name]
        except KeyError:
            arn = FakeTrialComponent.arn_formatter(
                trial_component_name, self.account_id, self.region_name
            )
            raise ValidationError(message=f"Could not find trial component '{arn}'")

        if status:
            trial_component.status = status
        if display_name:
            trial_component.display_name = display_name
        if start_time:
            trial_component.start_time = start_time
        if end_time:
            trial_component.end_time = end_time
        if parameters:
            trial_component.parameters = parameters
        if input_artifacts:
            trial_component.input_artifacts = input_artifacts
        if output_artifacts:
            trial_component.output_artifacts = output_artifacts

        trial_component.last_modified_time = datetime.now().strftime(
            "%Y-%m-%d %H:%M:%S"
        )

        for parameter_to_remove in parameters_to_remove or []:
            trial_component.parameters.pop(parameter_to_remove)

        for input_artifact_to_remove in input_artifacts_to_remove or []:
            trial_component.input_artifacts.pop(input_artifact_to_remove)

        for output_artifact_to_remove in output_artifacts_to_remove or []:
            trial_component.output_artifacts.pop(output_artifact_to_remove)

        return {
            "TrialComponentArn": FakeTrialComponent.arn_formatter(
                trial_component_name, self.account_id, self.region_name
            )
        }

    def create_notebook_instance(
        self,
        notebook_instance_name: str,
        instance_type: str,
        role_arn: str,
        subnet_id: Optional[str] = None,
        security_group_ids: Optional[List[str]] = None,
        kms_key_id: Optional[str] = None,
        tags: Optional[List[Dict[str, str]]] = None,
        lifecycle_config_name: Optional[str] = None,
        direct_internet_access: str = "Enabled",
        volume_size_in_gb: int = 5,
        accelerator_types: Optional[List[str]] = None,
        default_code_repository: Optional[str] = None,
        additional_code_repositories: Optional[List[str]] = None,
        root_access: Optional[str] = None,
    ) -> FakeSagemakerNotebookInstance:
        self._validate_unique_notebook_instance_name(notebook_instance_name)

        notebook_instance = FakeSagemakerNotebookInstance(
            account_id=self.account_id,
            region_name=self.region_name,
            notebook_instance_name=notebook_instance_name,
            instance_type=instance_type,
            role_arn=role_arn,
            subnet_id=subnet_id,
            security_group_ids=security_group_ids,
            kms_key_id=kms_key_id,
            tags=tags,
            lifecycle_config_name=lifecycle_config_name,
            direct_internet_access=direct_internet_access
            if direct_internet_access is not None
            else "Enabled",
            volume_size_in_gb=volume_size_in_gb if volume_size_in_gb is not None else 5,
            accelerator_types=accelerator_types,
            default_code_repository=default_code_repository,
            additional_code_repositories=additional_code_repositories,
            root_access=root_access,
        )
        self.notebook_instances[notebook_instance_name] = notebook_instance
        return notebook_instance

    def _validate_unique_notebook_instance_name(
        self, notebook_instance_name: str
    ) -> None:
        if notebook_instance_name in self.notebook_instances:
            duplicate_arn = self.notebook_instances[notebook_instance_name].arn
            message = f"Cannot create a duplicate Notebook Instance ({duplicate_arn})"
            raise ValidationError(message=message)

    def get_notebook_instance(
        self, notebook_instance_name: str
    ) -> FakeSagemakerNotebookInstance:
        try:
            return self.notebook_instances[notebook_instance_name]
        except KeyError:
            raise ValidationError(message="RecordNotFound")

    def start_notebook_instance(self, notebook_instance_name: str) -> None:
        notebook_instance = self.get_notebook_instance(notebook_instance_name)
        notebook_instance.start()

    def stop_notebook_instance(self, notebook_instance_name: str) -> None:
        notebook_instance = self.get_notebook_instance(notebook_instance_name)
        notebook_instance.stop()

    def delete_notebook_instance(self, notebook_instance_name: str) -> None:
        notebook_instance = self.get_notebook_instance(notebook_instance_name)
        if not notebook_instance.is_deletable:
            message = f"Status ({notebook_instance.status}) not in ([Stopped, Failed]). Unable to transition to (Deleting) for Notebook Instance ({notebook_instance.arn})"
            raise ValidationError(message=message)
        del self.notebook_instances[notebook_instance_name]

    @paginate(pagination_model=PAGINATION_MODEL)
    def list_notebook_instances(
        self,
        sort_by: str,
        sort_order: str,
        name_contains: Optional[str],
        status: Optional[str],
    ) -> List[FakeSagemakerNotebookInstance]:
        """
        The following parameters are not yet implemented:
        CreationTimeBefore, CreationTimeAfter, LastModifiedTimeBefore, LastModifiedTimeAfter, NotebookInstanceLifecycleConfigNameContains, DefaultCodeRepositoryContains, AdditionalCodeRepositoryEquals
        """
        instances = list(self.notebook_instances.values())
        if name_contains:
            instances = [
                i for i in instances if name_contains in i.notebook_instance_name
            ]
        if status:
            instances = [i for i in instances if i.status == status]
        reverse = sort_order == "Descending"
        if sort_by == "Name":
            instances = sorted(
                instances, key=lambda x: x.notebook_instance_name, reverse=reverse
            )
        if sort_by == "CreationTime":
            instances = sorted(
                instances, key=lambda x: x.creation_time, reverse=reverse
            )
        if sort_by == "Status":
            instances = sorted(instances, key=lambda x: x.status, reverse=reverse)
        return instances

    def create_notebook_instance_lifecycle_config(
        self,
        notebook_instance_lifecycle_config_name: str,
        on_create: List[Dict[str, str]],
        on_start: List[Dict[str, str]],
    ) -> FakeSageMakerNotebookInstanceLifecycleConfig:
        if (
            notebook_instance_lifecycle_config_name
            in self.notebook_instance_lifecycle_configurations
        ):
            arn = FakeSageMakerNotebookInstanceLifecycleConfig.arn_formatter(
                notebook_instance_lifecycle_config_name,
                self.account_id,
                self.region_name,
            )
            message = f"Unable to create Notebook Instance Lifecycle Config {arn}. (Details: Notebook Instance Lifecycle Config already exists.)"
            raise ValidationError(message=message)
        lifecycle_config = FakeSageMakerNotebookInstanceLifecycleConfig(
            account_id=self.account_id,
            region_name=self.region_name,
            notebook_instance_lifecycle_config_name=notebook_instance_lifecycle_config_name,
            on_create=on_create,
            on_start=on_start,
        )
        self.notebook_instance_lifecycle_configurations[
            notebook_instance_lifecycle_config_name
        ] = lifecycle_config
        return lifecycle_config

    def describe_notebook_instance_lifecycle_config(
        self, notebook_instance_lifecycle_config_name: str
    ) -> Dict[str, Any]:
        try:
            return self.notebook_instance_lifecycle_configurations[
                notebook_instance_lifecycle_config_name
            ].response_object
        except KeyError:
            arn = FakeSageMakerNotebookInstanceLifecycleConfig.arn_formatter(
                notebook_instance_lifecycle_config_name,
                self.account_id,
                self.region_name,
            )
            message = f"Unable to describe Notebook Instance Lifecycle Config '{arn}'. (Details: Notebook Instance Lifecycle Config does not exist.)"
            raise ValidationError(message=message)

    def delete_notebook_instance_lifecycle_config(
        self, notebook_instance_lifecycle_config_name: str
    ) -> None:
        try:
            del self.notebook_instance_lifecycle_configurations[
                notebook_instance_lifecycle_config_name
            ]
        except KeyError:
            arn = FakeSageMakerNotebookInstanceLifecycleConfig.arn_formatter(
                notebook_instance_lifecycle_config_name,
                self.account_id,
                self.region_name,
            )
            message = f"Unable to delete Notebook Instance Lifecycle Config '{arn}'. (Details: Notebook Instance Lifecycle Config does not exist.)"
            raise ValidationError(message=message)

    def create_endpoint_config(
        self,
        endpoint_config_name: str,
        production_variants: List[Dict[str, Any]],
        data_capture_config: Dict[str, Any],
        tags: List[Dict[str, str]],
        kms_key_id: str,
    ) -> FakeEndpointConfig:
        endpoint_config = FakeEndpointConfig(
            account_id=self.account_id,
            region_name=self.region_name,
            endpoint_config_name=endpoint_config_name,
            production_variants=production_variants,
            data_capture_config=data_capture_config,
            tags=tags,
            kms_key_id=kms_key_id,
        )
        self.validate_production_variants(production_variants)

        self.endpoint_configs[endpoint_config_name] = endpoint_config
        return endpoint_config

    def validate_production_variants(
        self, production_variants: List[Dict[str, Any]]
    ) -> None:
        for production_variant in production_variants:
            if production_variant["ModelName"] not in self._models:
                arn = arn_formatter(
                    "model",
                    production_variant["ModelName"],
                    self.account_id,
                    self.region_name,
                )
                raise ValidationError(message=f"Could not find model '{arn}'.")

    def describe_endpoint_config(self, endpoint_config_name: str) -> Dict[str, Any]:
        try:
            return self.endpoint_configs[endpoint_config_name].response_object
        except KeyError:
            arn = FakeEndpointConfig.arn_formatter(
                endpoint_config_name, self.account_id, self.region_name
            )
            raise ValidationError(
                message=f"Could not find endpoint configuration '{arn}'."
            )

    def delete_endpoint_config(self, endpoint_config_name: str) -> None:
        try:
            del self.endpoint_configs[endpoint_config_name]
        except KeyError:
            arn = FakeEndpointConfig.arn_formatter(
                endpoint_config_name, self.account_id, self.region_name
            )
            raise ValidationError(
                message=f"Could not find endpoint configuration '{arn}'."
            )

    def create_endpoint(
        self, endpoint_name: str, endpoint_config_name: str, tags: List[Dict[str, str]]
    ) -> FakeEndpoint:
        try:
            endpoint_config = self.describe_endpoint_config(endpoint_config_name)
        except KeyError:
            arn = FakeEndpointConfig.arn_formatter(
                endpoint_config_name, self.account_id, self.region_name
            )
            raise ValidationError(message=f"Could not find endpoint_config '{arn}'.")

        endpoint = FakeEndpoint(
            account_id=self.account_id,
            region_name=self.region_name,
            endpoint_name=endpoint_name,
            endpoint_config_name=endpoint_config_name,
            production_variants=endpoint_config["ProductionVariants"],
            data_capture_config=endpoint_config["DataCaptureConfig"],
            tags=tags,
        )

        self.endpoints[endpoint_name] = endpoint
        return endpoint

    def describe_endpoint(self, endpoint_name: str) -> Dict[str, Any]:
        try:
            return self.endpoints[endpoint_name].response_object
        except KeyError:
            arn = FakeEndpoint.arn_formatter(
                endpoint_name, self.account_id, self.region_name
            )
            raise ValidationError(message=f"Could not find endpoint '{arn}'.")

    def delete_endpoint(self, endpoint_name: str) -> None:
        try:
            del self.endpoints[endpoint_name]
        except KeyError:
            arn = FakeEndpoint.arn_formatter(
                endpoint_name, self.account_id, self.region_name
            )
            raise ValidationError(message=f"Could not find endpoint '{arn}'.")

    def create_processing_job(
        self,
        app_specification: Dict[str, Any],
        experiment_config: Dict[str, str],
        network_config: Dict[str, Any],
        processing_inputs: List[Dict[str, Any]],
        processing_job_name: str,
        processing_output_config: Dict[str, Any],
        role_arn: str,
        tags: List[Dict[str, str]],
        stopping_condition: Dict[str, int],
    ) -> FakeProcessingJob:
        processing_job = FakeProcessingJob(
            app_specification=app_specification,
            experiment_config=experiment_config,
            network_config=network_config,
            processing_inputs=processing_inputs,
            processing_job_name=processing_job_name,
            processing_output_config=processing_output_config,
            account_id=self.account_id,
            region_name=self.region_name,
            role_arn=role_arn,
            stopping_condition=stopping_condition,
            tags=tags,
        )
        self.processing_jobs[processing_job_name] = processing_job
        return processing_job

    def describe_processing_job(self, processing_job_name: str) -> Dict[str, Any]:
        try:
            return self.processing_jobs[processing_job_name].response_object
        except KeyError:
            arn = FakeProcessingJob.arn_formatter(
                processing_job_name, self.account_id, self.region_name
            )
            raise ValidationError(message=f"Could not find processing job '{arn}'.")

    def create_pipeline(
        self,
        pipeline_name: str,
        pipeline_display_name: str,
        pipeline_definition: str,
        pipeline_definition_s3_location: Dict[str, Any],
        pipeline_description: str,
        role_arn: str,
        tags: List[Dict[str, str]],
        parallelism_configuration: Dict[str, int],
    ) -> FakePipeline:
        if not any([pipeline_definition, pipeline_definition_s3_location]):
            raise ValidationError(
                "An error occurred (ValidationException) when calling the CreatePipeline operation: Either "
                "Pipeline Definition or Pipeline Definition S3 location should be provided"
            )
        if all([pipeline_definition, pipeline_definition_s3_location]):
            raise ValidationError(
                "An error occurred (ValidationException) when calling the CreatePipeline operation: "
                "Both Pipeline Definition and Pipeline Definition S3 Location shouldn't be present"
            )

        if pipeline_name in self.pipelines:
            raise ValidationError(
                f"An error occurred (ValidationException) when calling the CreatePipeline operation: Pipeline names "
                f"must be unique within an AWS account and region. Pipeline with name ({pipeline_name}) already exists."
            )

        if pipeline_definition_s3_location:
            pipeline_definition = load_pipeline_definition_from_s3(  # type: ignore
                pipeline_definition_s3_location,
                account_id=self.account_id,
                partition=self.partition,
            )

        pipeline = FakePipeline(
            pipeline_name,
            pipeline_display_name,
            pipeline_definition,
            pipeline_description,
            role_arn,
            tags,
            self.account_id,
            self.region_name,
            parallelism_configuration,
        )

        self.pipelines[pipeline_name] = pipeline
        return pipeline

    def delete_pipeline(self, pipeline_name: str) -> str:
        pipeline = get_pipeline_from_name(self.pipelines, pipeline_name)
        del self.pipelines[pipeline.pipeline_name]
        return pipeline.arn

    def update_pipeline(self, pipeline_name: str, **kwargs: Any) -> str:
        pipeline = get_pipeline_from_name(self.pipelines, pipeline_name)
        if all(
            [
                kwargs.get("pipeline_definition"),
                kwargs.get("pipeline_definition_s3_location"),
            ]
        ):
            raise ValidationError(
                "An error occurred (ValidationException) when calling the UpdatePipeline operation: "
                "Both Pipeline Definition and Pipeline Definition S3 Location shouldn't be present"
            )

        for attr_key, attr_value in kwargs.items():
            if attr_value:
                if attr_key == "pipeline_definition_s3_location":
                    self.pipelines[
                        pipeline_name
                    ].pipeline_definition = load_pipeline_definition_from_s3(  # type: ignore
                        attr_value,
                        self.account_id,
                        partition=self.partition,
                    )
                    continue
                setattr(self.pipelines[pipeline_name], attr_key, attr_value)

        return pipeline.arn

    def start_pipeline_execution(
        self,
        pipeline_name: str,
        pipeline_execution_display_name: str,
        pipeline_parameters: List[Dict[str, Any]],
        pipeline_execution_description: str,
        parallelism_configuration: Dict[str, int],
        client_request_token: str,
    ) -> Dict[str, str]:
        pipeline = get_pipeline_from_name(self.pipelines, pipeline_name)
        execution_id = "".join(
            random.choices(string.ascii_lowercase + string.digits, k=12)
        )
        pipeline_execution_arn = arn_formatter(
            _type="pipeline",
            _id=f"{pipeline.pipeline_name}/execution/{execution_id}",
            account_id=self.account_id,
            region_name=self.region_name,
        )

        fake_pipeline_execution = FakePipelineExecution(
            pipeline_execution_arn=pipeline_execution_arn,
            pipeline_execution_display_name=pipeline_execution_display_name,
            pipeline_parameters=pipeline_parameters,
            pipeline_execution_description=pipeline_execution_description,
            pipeline_definition=pipeline.pipeline_definition,
            parallelism_configuration=parallelism_configuration
            or pipeline.parallelism_configuration,
            client_request_token=client_request_token,
        )

        self.pipelines[pipeline_name].pipeline_executions[pipeline_execution_arn] = (
            fake_pipeline_execution
        )
        self.pipelines[
            pipeline_name
        ].last_execution_time = fake_pipeline_execution.start_time

        return {"PipelineExecutionArn": pipeline_execution_arn}

    def list_pipeline_executions(self, pipeline_name: str) -> Dict[str, Any]:
        pipeline = get_pipeline_from_name(self.pipelines, pipeline_name)
        return {
            "PipelineExecutionSummaries": [
                {
                    "PipelineExecutionArn": arn,
                    "StartTime": pipeline_execution.start_time,
                    "PipelineExecutionStatus": pipeline_execution.pipeline_execution_status,
                    "PipelineExecutionDescription": pipeline_execution.pipeline_execution_description,
                    "PipelineExecutionDisplayName": pipeline_execution.pipeline_execution_display_name,
                    "PipelineExecutionFailureReason": str(
                        pipeline_execution.pipeline_execution_failure_reason
                    ),
                }
                for arn, pipeline_execution in pipeline.pipeline_executions.items()
            ]
        }

    def describe_pipeline_definition_for_execution(
        self, pipeline_execution_arn: str
    ) -> Dict[str, Any]:
        pipeline_execution = get_pipeline_execution_from_arn(
            self.pipelines, pipeline_execution_arn
        )
        return {
            "PipelineDefinition": str(
                pipeline_execution.pipeline_definition_for_execution
            ),
            "CreationTime": pipeline_execution.creation_time,
        }

    def list_pipeline_parameters_for_execution(
        self, pipeline_execution_arn: str
    ) -> Dict[str, Any]:
        pipeline_execution = get_pipeline_execution_from_arn(
            self.pipelines, pipeline_execution_arn
        )
        return {
            "PipelineParameters": pipeline_execution.pipeline_parameters,
        }

    def describe_pipeline_execution(
        self, pipeline_execution_arn: str
    ) -> Dict[str, Any]:
        pipeline_execution = get_pipeline_execution_from_arn(
            self.pipelines, pipeline_execution_arn
        )
        pipeline_name = get_pipeline_name_from_execution_arn(pipeline_execution_arn)
        pipeline = get_pipeline_from_name(self.pipelines, pipeline_name)

        return {
            "PipelineArn": pipeline.arn,
            "PipelineExecutionArn": pipeline_execution.arn,
            "PipelineExecutionDisplayName": pipeline_execution.pipeline_execution_display_name,
            "PipelineExecutionStatus": pipeline_execution.pipeline_execution_status,
            "PipelineExecutionDescription": pipeline_execution.pipeline_execution_description,
            "PipelineExperimentConfig": {},
            "FailureReason": "",
            "CreationTime": pipeline_execution.creation_time,
            "LastModifiedTime": pipeline_execution.last_modified_time,
            "CreatedBy": pipeline_execution.created_by,
            "LastModifiedBy": pipeline_execution.last_modified_by,
            "ParallelismConfiguration": pipeline_execution.parallelism_configuration,
        }

    def describe_pipeline(self, pipeline_name: str) -> Dict[str, Any]:
        pipeline = get_pipeline_from_name(self.pipelines, pipeline_name)
        return {
            "PipelineArn": pipeline.arn,
            "PipelineName": pipeline.pipeline_name,
            "PipelineDisplayName": pipeline.pipeline_display_name,
            "PipelineDescription": pipeline.pipeline_description,
            "PipelineDefinition": pipeline.pipeline_definition,
            "RoleArn": pipeline.role_arn,
            "PipelineStatus": pipeline.pipeline_status,
            "CreationTime": pipeline.creation_time,
            "LastModifiedTime": pipeline.last_modified_time,
            "LastRunTime": pipeline.last_execution_time,
            "CreatedBy": pipeline.created_by,
            "LastModifiedBy": pipeline.last_modified_by,
            "ParallelismConfiguration": pipeline.parallelism_configuration,
        }

    def list_pipelines(
        self,
        pipeline_name_prefix: str,
        created_after: str,
        created_before: str,
        next_token: str,
        max_results: int,
        sort_by: str,
        sort_order: str,
    ) -> Dict[str, Any]:
        if next_token:
            try:
                starting_index = int(next_token)
                if starting_index > len(self.pipelines):
                    raise ValueError  # invalid next_token
            except ValueError:
                raise AWSValidationException('Invalid pagination token because "{0}".')
        else:
            starting_index = 0

        if max_results:
            end_index = max_results + starting_index
            pipelines_fetched: Iterable[FakePipeline] = list(self.pipelines.values())[
                starting_index:end_index
            ]
            if end_index >= len(self.pipelines):
                next_index = None
            else:
                next_index = end_index
        else:
            pipelines_fetched = list(self.pipelines.values())
            next_index = None

        if pipeline_name_prefix is not None:
            pipelines_fetched = filter(
                lambda x: pipeline_name_prefix in x.pipeline_name,
                pipelines_fetched,
            )

        def format_time(x: Any) -> str:
            return (
                x
                if isinstance(x, str)
                else datetime.fromtimestamp(x).strftime("%Y-%m-%d %H:%M:%S")
            )

        if created_after is not None:
            pipelines_fetched = filter(
                lambda x: x.creation_time > format_time(created_after),
                pipelines_fetched,
            )

        if created_before is not None:
            pipelines_fetched = filter(
                lambda x: x.creation_time < format_time(created_before),
                pipelines_fetched,
            )

        sort_key = "pipeline_name" if sort_by == "Name" else "creation_time"
        pipelines_fetched = sorted(
            pipelines_fetched,
            key=lambda pipeline_fetched: getattr(pipeline_fetched, sort_key),
            reverse=sort_order != "Ascending",
        )

        pipeline_summaries = [
            {
                "PipelineArn": pipeline_data.arn,
                "PipelineName": pipeline_data.pipeline_name,
                "PipelineDisplayName": pipeline_data.pipeline_display_name,
                "PipelineDescription": pipeline_data.pipeline_description,
                "RoleArn": pipeline_data.role_arn,
                "CreationTime": pipeline_data.creation_time,
                "LastModifiedTime": pipeline_data.last_modified_time,
                "LastExecutionTime": pipeline_data.last_execution_time,
            }
            for pipeline_data in pipelines_fetched
        ]

        return {
            "PipelineSummaries": pipeline_summaries,
            "NextToken": str(next_index) if next_index is not None else None,
        }

    def list_processing_jobs(
        self,
        next_token: str,
        max_results: int,
        creation_time_after: str,
        creation_time_before: str,
        last_modified_time_after: str,
        last_modified_time_before: str,
        name_contains: str,
        status_equals: str,
    ) -> Dict[str, Any]:
        if next_token:
            try:
                starting_index = int(next_token)
                if starting_index > len(self.processing_jobs):
                    raise ValueError  # invalid next_token
            except ValueError:
                raise AWSValidationException('Invalid pagination token because "{0}".')
        else:
            starting_index = 0

        if max_results:
            end_index = max_results + starting_index
            processing_jobs_fetched: Iterable[FakeProcessingJob] = list(
                self.processing_jobs.values()
            )[starting_index:end_index]
            if end_index >= len(self.processing_jobs):
                next_index = None
            else:
                next_index = end_index
        else:
            processing_jobs_fetched = list(self.processing_jobs.values())
            next_index = None

        if name_contains is not None:
            processing_jobs_fetched = filter(
                lambda x: name_contains in x.processing_job_name,
                processing_jobs_fetched,
            )

        if creation_time_after is not None:
            processing_jobs_fetched = filter(
                lambda x: x.creation_time > creation_time_after, processing_jobs_fetched
            )

        if creation_time_before is not None:
            processing_jobs_fetched = filter(
                lambda x: x.creation_time < creation_time_before,
                processing_jobs_fetched,
            )

        if last_modified_time_after is not None:
            processing_jobs_fetched = filter(
                lambda x: x.last_modified_time > last_modified_time_after,
                processing_jobs_fetched,
            )

        if last_modified_time_before is not None:
            processing_jobs_fetched = filter(
                lambda x: x.last_modified_time < last_modified_time_before,
                processing_jobs_fetched,
            )
        if status_equals is not None:
            processing_jobs_fetched = filter(
                lambda x: x.processing_job_status == status_equals,
                processing_jobs_fetched,
            )

        processing_job_summaries = [
            {
                "ProcessingJobName": processing_job_data.processing_job_name,
                "ProcessingJobArn": processing_job_data.arn,
                "CreationTime": processing_job_data.creation_time,
                "ProcessingEndTime": processing_job_data.processing_end_time,
                "LastModifiedTime": processing_job_data.last_modified_time,
                "ProcessingJobStatus": processing_job_data.processing_job_status,
            }
            for processing_job_data in processing_jobs_fetched
        ]

        return {
            "ProcessingJobSummaries": processing_job_summaries,
            "NextToken": str(next_index) if next_index is not None else None,
        }

    def create_transform_job(
        self,
        transform_job_name: str,
        model_name: str,
        max_concurrent_transforms: int,
        model_client_config: Dict[str, int],
        max_payload_in_mb: int,
        batch_strategy: str,
        environment: Dict[str, str],
        transform_input: Dict[str, Union[Dict[str, str], str]],
        transform_output: Dict[str, str],
        data_capture_config: Dict[str, Union[str, bool]],
        transform_resources: Dict[str, Union[str, int]],
        data_processing: Dict[str, str],
        tags: Dict[str, str],
        experiment_config: Dict[str, str],
    ) -> FakeTransformJob:
        transform_job = FakeTransformJob(
            account_id=self.account_id,
            region_name=self.region_name,
            transform_job_name=transform_job_name,
            model_name=model_name,
            max_concurrent_transforms=max_concurrent_transforms,
            model_client_config=model_client_config,
            max_payload_in_mb=max_payload_in_mb,
            batch_strategy=batch_strategy,
            environment=environment,
            transform_input=transform_input,
            transform_output=transform_output,
            data_capture_config=data_capture_config,
            transform_resources=transform_resources,
            data_processing=data_processing,
            tags=tags,
            experiment_config=experiment_config,
        )
        self.transform_jobs[transform_job_name] = transform_job
        return transform_job

    def list_transform_jobs(
        self,
        next_token: str,
        max_results: int,
        creation_time_after: str,
        creation_time_before: str,
        last_modified_time_after: str,
        last_modified_time_before: str,
        name_contains: str,
        status_equals: str,
    ) -> Dict[str, Any]:
        if next_token:
            try:
                starting_index = int(next_token)
                if starting_index > len(self.transform_jobs):
                    raise ValueError  # invalid next_token
            except ValueError:
                raise AWSValidationException('Invalid pagination token because "{0}".')
        else:
            starting_index = 0

        if max_results:
            end_index = max_results + starting_index
            transform_jobs_fetched: Iterable[FakeTransformJob] = list(
                self.transform_jobs.values()
            )[starting_index:end_index]
            if end_index >= len(self.transform_jobs):
                next_index = None
            else:
                next_index = end_index
        else:
            transform_jobs_fetched = list(self.transform_jobs.values())
            next_index = None

        if name_contains is not None:
            transform_jobs_fetched = filter(
                lambda x: name_contains in x.transform_job_name, transform_jobs_fetched
            )

        if creation_time_after is not None:
            transform_jobs_fetched = filter(
                lambda x: x.creation_time > creation_time_after, transform_jobs_fetched
            )

        if creation_time_before is not None:
            transform_jobs_fetched = filter(
                lambda x: x.creation_time < creation_time_before, transform_jobs_fetched
            )

        if last_modified_time_after is not None:
            transform_jobs_fetched = filter(
                lambda x: x.last_modified_time > last_modified_time_after,
                transform_jobs_fetched,
            )

        if last_modified_time_before is not None:
            transform_jobs_fetched = filter(
                lambda x: x.last_modified_time < last_modified_time_before,
                transform_jobs_fetched,
            )
        if status_equals is not None:
            transform_jobs_fetched = filter(
                lambda x: x.transform_job_status == status_equals,
                transform_jobs_fetched,
            )

        transform_job_summaries = [
            {
                "TransformJobName": transform_job_data.transform_job_name,
                "TransformJobArn": transform_job_data.arn,
                "CreationTime": transform_job_data.creation_time,
                "TransformEndTime": transform_job_data.transform_end_time,
                "LastModifiedTime": transform_job_data.last_modified_time,
                "TransformJobStatus": transform_job_data.transform_job_status,
            }
            for transform_job_data in transform_jobs_fetched
        ]

        return {
            "TransformJobSummaries": transform_job_summaries,
            "NextToken": str(next_index) if next_index is not None else None,
        }

    def describe_transform_job(self, transform_job_name: str) -> Dict[str, Any]:
        try:
            return self.transform_jobs[transform_job_name].response_object
        except KeyError:
            arn = FakeTransformJob.arn_formatter(
                transform_job_name, self.account_id, self.region_name
            )
            message = f"Could not find transform job '{arn}'."
            raise ValidationError(message=message)

    def create_training_job(
        self,
        training_job_name: str,
        hyper_parameters: Dict[str, str],
        algorithm_specification: Dict[str, Any],
        role_arn: str,
        input_data_config: List[Dict[str, Any]],
        output_data_config: Dict[str, str],
        resource_config: Dict[str, Any],
        vpc_config: Dict[str, List[str]],
        stopping_condition: Dict[str, int],
        tags: List[Dict[str, str]],
        enable_network_isolation: bool,
        enable_inter_container_traffic_encryption: bool,
        enable_managed_spot_training: bool,
        checkpoint_config: Dict[str, str],
        debug_hook_config: Dict[str, Any],
        debug_rule_configurations: List[Dict[str, Any]],
        tensor_board_output_config: Dict[str, str],
        experiment_config: Dict[str, str],
    ) -> FakeTrainingJob:
        training_job = FakeTrainingJob(
            account_id=self.account_id,
            region_name=self.region_name,
            training_job_name=training_job_name,
            hyper_parameters=hyper_parameters,
            algorithm_specification=algorithm_specification,
            role_arn=role_arn,
            input_data_config=input_data_config,
            output_data_config=output_data_config,
            resource_config=resource_config,
            vpc_config=vpc_config,
            stopping_condition=stopping_condition,
            tags=tags,
            enable_network_isolation=enable_network_isolation,
            enable_inter_container_traffic_encryption=enable_inter_container_traffic_encryption,
            enable_managed_spot_training=enable_managed_spot_training,
            checkpoint_config=checkpoint_config,
            debug_hook_config=debug_hook_config,
            debug_rule_configurations=debug_rule_configurations,
            tensor_board_output_config=tensor_board_output_config,
            experiment_config=experiment_config,
        )
        self.training_jobs[training_job_name] = training_job
        return training_job

    def describe_training_job(self, training_job_name: str) -> Dict[str, Any]:
        try:
            return self.training_jobs[training_job_name].response_object
        except KeyError:
            arn = FakeTrainingJob.arn_formatter(
                training_job_name, self.account_id, self.region_name
            )
            message = f"Could not find training job '{arn}'."
            raise ValidationError(message=message)

    def list_training_jobs(
        self,
        next_token: str,
        max_results: int,
        creation_time_after: str,
        creation_time_before: str,
        last_modified_time_after: str,
        last_modified_time_before: str,
        name_contains: str,
        status_equals: str,
    ) -> Dict[str, Any]:
        if next_token:
            try:
                starting_index = int(next_token)
                if starting_index > len(self.training_jobs):
                    raise ValueError  # invalid next_token
            except ValueError:
                raise AWSValidationException('Invalid pagination token because "{0}".')
        else:
            starting_index = 0

        if max_results:
            end_index = max_results + starting_index
            training_jobs_fetched: Iterable[FakeTrainingJob] = list(
                self.training_jobs.values()
            )[starting_index:end_index]
            if end_index >= len(self.training_jobs):
                next_index = None
            else:
                next_index = end_index
        else:
            training_jobs_fetched = list(self.training_jobs.values())
            next_index = None

        if name_contains is not None:
            training_jobs_fetched = filter(
                lambda x: name_contains in x.training_job_name, training_jobs_fetched
            )

        if creation_time_after is not None:
            training_jobs_fetched = filter(
                lambda x: x.creation_time > creation_time_after, training_jobs_fetched
            )

        if creation_time_before is not None:
            training_jobs_fetched = filter(
                lambda x: x.creation_time < creation_time_before, training_jobs_fetched
            )

        if last_modified_time_after is not None:
            training_jobs_fetched = filter(
                lambda x: x.last_modified_time > last_modified_time_after,
                training_jobs_fetched,
            )

        if last_modified_time_before is not None:
            training_jobs_fetched = filter(
                lambda x: x.last_modified_time < last_modified_time_before,
                training_jobs_fetched,
            )
        if status_equals is not None:
            training_jobs_fetched = filter(
                lambda x: x.training_job_status == status_equals, training_jobs_fetched
            )

        training_job_summaries = [
            {
                "TrainingJobName": training_job_data.training_job_name,
                "TrainingJobArn": training_job_data.arn,
                "CreationTime": training_job_data.creation_time,
                "TrainingEndTime": training_job_data.training_end_time,
                "LastModifiedTime": training_job_data.last_modified_time,
                "TrainingJobStatus": training_job_data.training_job_status,
            }
            for training_job_data in training_jobs_fetched
        ]

        return {
            "TrainingJobSummaries": training_job_summaries,
            "NextToken": str(next_index) if next_index is not None else None,
        }

    def update_endpoint_weights_and_capacities(
        self, endpoint_name: str, desired_weights_and_capacities: List[Dict[str, Any]]
    ) -> str:
        # Validate inputs
        endpoint = self.endpoints.get(endpoint_name, None)
        if not endpoint:
            arn = FakeEndpoint.arn_formatter(
                endpoint_name, self.account_id, self.region_name
            )
            raise AWSValidationException(f'Could not find endpoint "{arn}".')

        names_checked = []
        for variant_config in desired_weights_and_capacities:
            name = variant_config.get("VariantName")

            if name in names_checked:
                raise AWSValidationException(
                    f'The variant name "{name}" was non-unique within the request.'
                )

            if not any(
                variant["VariantName"] == name
                for variant in endpoint.production_variants
            ):
                raise AWSValidationException(
                    f'The variant name(s) "{name}" is/are not present within endpoint configuration "{endpoint.endpoint_config_name}".'
                )

            names_checked.append(name)

        # Update endpoint variants
        endpoint.endpoint_status = "Updating"

        for variant_config in desired_weights_and_capacities:
            name = variant_config.get("VariantName")
            desired_weight = variant_config.get("DesiredWeight")
            desired_instance_count = variant_config.get("DesiredInstanceCount")

            for variant in endpoint.production_variants:
                if variant.get("VariantName") == name:
                    variant["DesiredWeight"] = desired_weight
                    variant["CurrentWeight"] = desired_weight
                    variant["DesiredInstanceCount"] = desired_instance_count
                    variant["CurrentInstanceCount"] = desired_instance_count
                    break

        endpoint.endpoint_status = "InService"
        return endpoint.arn

    def create_model_package_group(
        self,
        model_package_group_name: str,
        model_package_group_description: str,
        tags: Optional[List[Dict[str, str]]] = None,
    ) -> str:
        self.model_package_groups[model_package_group_name] = ModelPackageGroup(
            model_package_group_name=model_package_group_name,
            model_package_group_description=model_package_group_description,
            account_id=self.account_id,
            region_name=self.region_name,
            tags=tags or [],
        )
        return self.model_package_groups[model_package_group_name].arn

    def _get_versioned_or_not(
        self, model_package_type: Optional[str], model_package_version: Optional[int]
    ) -> bool:
        if model_package_type == "Versioned":
            return model_package_version is not None
        elif model_package_type == "Unversioned" or model_package_type is None:
            return model_package_version is None
        elif model_package_type == "Both":
            return True
        raise ValueError(f"Invalid model package type: {model_package_type}")

    @paginate(pagination_model=PAGINATION_MODEL)
    def list_model_package_groups(
        self,
        creation_time_after: Optional[int],
        creation_time_before: Optional[int],
        name_contains: Optional[str],
        sort_by: Optional[str],
        sort_order: Optional[str],
    ) -> List[ModelPackageGroup]:
        if isinstance(creation_time_before, int):
            creation_time_before_datetime = datetime.fromtimestamp(
                creation_time_before, tz=tzutc()
            )
        if isinstance(creation_time_after, int):
            creation_time_after_datetime = datetime.fromtimestamp(
                creation_time_after, tz=tzutc()
            )
        model_package_group_summary_list = list(
            filter(
                lambda x: (
                    creation_time_after is None
                    or x.creation_time > creation_time_after_datetime
                )
                and (
                    creation_time_before is None
                    or x.creation_time < creation_time_before_datetime
                )
                and (
                    name_contains is None
                    or x.model_package_group_name.find(name_contains) != -1
                ),
                self.model_package_groups.values(),
            )
        )
        model_package_group_summary_list = list(
            sorted(
                model_package_group_summary_list,
                key={
                    "Name": lambda x: x.model_package_group_name,
                    "CreationTime": lambda x: x.creation_time,
                    None: lambda x: x.creation_time,
                }[sort_by],
                reverse=sort_order == "Descending",
            )
        )
        return model_package_group_summary_list

    def describe_model_package_group(
        self, model_package_group_name: str
    ) -> ModelPackageGroup:
        model_package_group = self.model_package_groups.get(model_package_group_name)
        if model_package_group is None:
            model_package_group_arn = arn_formatter(
                region_name=self.region_name,
                account_id=self.account_id,
                _type="model-package-group",
                _id=f"{model_package_group_name}",
            )
            raise ValidationError(
                f"ModelPackageGroup {model_package_group_arn} does not exist."
            )
        return model_package_group

    @paginate(pagination_model=PAGINATION_MODEL)
    def list_model_packages(
        self,
        creation_time_after: Optional[int],
        creation_time_before: Optional[int],
        name_contains: Optional[str],
        model_approval_status: Optional[str],
        model_package_group_name: Optional[str],
        model_package_type: Optional[str],
        sort_by: Optional[str],
        sort_order: Optional[str],
    ) -> List[ModelPackage]:
        if isinstance(creation_time_before, int):
            creation_time_before_datetime = datetime.fromtimestamp(
                creation_time_before, tz=tzutc()
            )
        if isinstance(creation_time_after, int):
            creation_time_after_datetime = datetime.fromtimestamp(
                creation_time_after, tz=tzutc()
            )
        if model_package_group_name is not None:
            model_package_type = "Versioned"
            if re.match(ARN_PARTITION_REGEX, model_package_group_name):
                model_package_group_name = model_package_group_name.split("/")[-1]
        model_package_summary_list = list(
            filter(
                lambda x: (
                    creation_time_after is None
                    or x.creation_time > creation_time_after_datetime
                )
                and (
                    creation_time_before is None
                    or x.creation_time < creation_time_before_datetime
                )
                and (
                    name_contains is None
                    or x.model_package_name.find(name_contains) != -1
                )
                and (
                    model_approval_status is None
                    or x.model_approval_status == model_approval_status
                )
                and (
                    model_package_group_name is None
                    or x.model_package_group_name == model_package_group_name
                )
                and self._get_versioned_or_not(
                    model_package_type, x.model_package_version
                ),
                self.model_packages.values(),
            )
        )
        model_package_summary_list = list(
            sorted(
                model_package_summary_list,
                key={
                    "Name": lambda x: x.model_package_name,
                    "CreationTime": lambda x: x.creation_time,
                    None: lambda x: x.creation_time,
                }[sort_by],
                reverse=sort_order == "Descending",
            )
        )
        return model_package_summary_list

    def describe_model_package(self, model_package_name: str) -> ModelPackage:
        model_package_name_mapped = self.model_package_name_mapping.get(
            model_package_name, model_package_name
        )
        model_package = self.model_packages.get(model_package_name_mapped)
        if model_package is None:
            raise ValidationError(f"Model package {model_package_name} not found")
        return model_package

    def update_model_package(
        self,
        model_package_arn: str,
        model_approval_status: Optional[str],
        approval_description: Optional[str],
        customer_metadata_properties: Optional[Dict[str, str]],
        customer_metadata_properties_to_remove: List[str],
        additional_inference_specifications_to_add: Optional[List[Any]],
    ) -> str:
        model_package_name_mapped = self.model_package_name_mapping.get(
            model_package_arn, model_package_arn
        )
        model_package = self.model_packages.get(model_package_name_mapped)

        if model_package is None:
            raise ValidationError(f"Model package {model_package_arn} not found")

        model_package.set_model_approval_status(model_approval_status)
        model_package.approval_description = approval_description
        model_package.customer_metadata_properties = customer_metadata_properties
        model_package.remove_customer_metadata_property(
            customer_metadata_properties_to_remove
        )
        model_package.add_additional_inference_specifications(
            additional_inference_specifications_to_add
        )
        model_package.modifications_done()

        return model_package.arn

    def create_model_package(
        self,
        model_package_name: Optional[str],
        model_package_group_name: Optional[str],
        model_package_description: Optional[str],
        inference_specification: Any,
        validation_specification: Any,
        source_algorithm_specification: Any,
        certify_for_marketplace: bool,
        tags: Any,
        model_approval_status: Optional[str],
        metadata_properties: Any,
        model_metrics: Any,
        client_token: Any,
        customer_metadata_properties: Any,
        drift_check_baselines: Any,
        domain: Any,
        task: Any,
        sample_payload_url: Any,
        additional_inference_specifications: Any,
    ) -> str:
        model_package_version = None
        if model_package_group_name and model_package_name:
            raise AWSValidationException(
                "An error occurred (ValidationException) when calling the CreateModelPackage operation: Both ModelPackageName and ModelPackageGroupName are provided in the input. Cannot determine which one to use."
            )
        elif not model_package_group_name and not model_package_name:
            raise AWSValidationException(
                "An error ocurred (ValidationException) when calling the CreateModelPackag operation: Missing ARN."
            )
        elif model_package_group_name:
            model_package_type = "Versioned"
            model_package_name = model_package_group_name
            model_packages_for_group = [
                x
                for x in self.model_packages.values()
                if x.model_package_group_name == model_package_group_name
            ]
            if model_package_group_name not in self.model_package_groups:
                raise AWSValidationException(
                    "An error ocurred (ValidationException) when calling the CreateModelPackage operation: Model Package Group does not exist."
                )
            model_package_version = len(model_packages_for_group) + 1
        else:
            model_package_type = "Unversioned"

        model_package = ModelPackage(
            model_package_name=cast(str, model_package_name),
            model_package_group_name=model_package_group_name,
            model_package_description=model_package_description,
            inference_specification=inference_specification,
            validation_specification=validation_specification,
            source_algorithm_specification=source_algorithm_specification,
            certify_for_marketplace=certify_for_marketplace,
            tags=tags,
            model_approval_status=model_approval_status,
            metadata_properties=metadata_properties,
            model_metrics=model_metrics,
            customer_metadata_properties=customer_metadata_properties,
            drift_check_baselines=drift_check_baselines,
            domain=domain,
            task=task,
            sample_payload_url=sample_payload_url,
            additional_inference_specifications=additional_inference_specifications,
            model_package_version=model_package_version,
            approval_description=None,
            region_name=self.region_name,
            account_id=self.account_id,
            client_token=client_token,
            model_package_type=model_package_type,
        )
        self.model_package_name_mapping[model_package.model_package_name] = (
            model_package.arn
        )
        self.model_package_name_mapping[model_package.arn] = model_package.arn
        self.model_packages[model_package.arn] = model_package
        return model_package.arn

    def create_feature_group(
        self,
        feature_group_name: str,
        record_identifier_feature_name: str,
        event_time_feature_name: str,
        feature_definitions: List[Dict[str, str]],
        offline_store_config: Dict[str, Any],
        role_arn: str,
        tags: Any,
    ) -> str:
        feature_group_arn = arn_formatter(
            region_name=self.region_name,
            account_id=self.account_id,
            _type="feature-group",
            _id=f"{feature_group_name.lower()}",
        )
        if feature_group_arn in self.feature_groups:
            raise ResourceInUseException(
                message=f"An error occurred (ResourceInUse) when calling the CreateFeatureGroup operation: Resource Already Exists: FeatureGroup with name {feature_group_name} already exists. Choose a different name.\nInfo: Feature Group '{feature_group_name}' already exists."
            )

        feature_group = FeatureGroup(
            feature_group_name=feature_group_name,
            record_identifier_feature_name=record_identifier_feature_name,
            event_time_feature_name=event_time_feature_name,
            feature_definitions=feature_definitions,
            offline_store_config=offline_store_config,
            role_arn=role_arn,
            region_name=self.region_name,
            account_id=self.account_id,
            tags=tags,
        )
        self.feature_groups[feature_group.arn] = feature_group
        return feature_group.arn

    def describe_feature_group(
        self,
        feature_group_name: str,
    ) -> Dict[str, Any]:
        feature_group_arn = arn_formatter(
            region_name=self.region_name,
            account_id=self.account_id,
            _type="feature-group",
            _id=f"{feature_group_name.lower()}",
        )

        feature_group = self.feature_groups[feature_group_arn]
        return feature_group.describe()

    def create_cluster(
        self,
        cluster_name: str,
        instance_groups: List[Dict[str, Any]],
        vpc_config: Dict[str, List[str]],
        tags: Any,
    ) -> str:
        cluster = Cluster(
            cluster_name=cluster_name,
            region_name=self.region_name,
            account_id=self.account_id,
            instance_groups=instance_groups,
            vpc_config=vpc_config,
            tags=tags,
        )
        self.clusters[cluster_name] = cluster

        # create Cluster Nodes
        for instance_group in instance_groups:
            for i in range(instance_group["TargetCount"]):
                node_id = f"{instance_group['InstanceGroupName']}-{i}"
                fake_cluster_node = ClusterNode(
                    region_name=self.region_name,
                    account_id=self.account_id,
                    cluster_name=cluster_name,
                    instance_group_name=instance_group["InstanceGroupName"],
                    instance_type=instance_group["InstanceType"],
                    life_cycle_config=instance_group["LifeCycleConfig"],
                    execution_role=instance_group["ExecutionRole"],
                    node_id=node_id,
                    threads_per_core=instance_group["ThreadsPerCore"],
                )
                cluster.nodes[node_id] = fake_cluster_node

        return cluster.arn

    def describe_cluster(self, cluster_name: str) -> Dict[str, Any]:
        if cluster_name.startswith(f"arn:{self.partition}:sagemaker:"):
            cluster_name = (cluster_name.split(":")[-1]).split("/")[-1]
        cluster = self.clusters.get(cluster_name)
        if not cluster:
            raise ValidationError(message=f"Could not find cluster '{cluster_name}'.")
        return cluster.describe()

    def delete_cluster(self, cluster_name: str) -> str:
        if cluster_name.startswith(f"arn:{self.partition}:sagemaker:"):
            cluster_name = (cluster_name.split(":")[-1]).split("/")[-1]
        cluster = self.clusters.get(cluster_name)
        if not cluster:
            raise ValidationError(message=f"Could not find cluster '{cluster_name}'.")
        arn = cluster.arn

        del self.clusters[cluster_name]
        return arn

    def describe_cluster_node(self, cluster_name: str, node_id: str) -> Dict[str, Any]:
        if cluster_name.startswith(f"arn:{self.partition}:sagemaker:"):
            cluster_name = (cluster_name.split(":")[-1]).split("/")[-1]
        cluster = self.clusters.get(cluster_name)
        if not cluster:
            raise ValidationError(message=f"Could not find cluster '{cluster_name}'.")
        if node_id in cluster.nodes:
            return cluster.nodes[node_id].describe()
        else:
            raise ValidationError(
                message=f"Could not find node '{node_id}' in cluster '{cluster_name}'."
            )

    @paginate(pagination_model=PAGINATION_MODEL)
    def list_clusters(
        self,
        creation_time_after: Optional[datetime],
        creation_time_before: Optional[datetime],
        name_contains: Optional[str],
        sort_by: Optional[str],
        sort_order: Optional[str],
    ) -> List[Cluster]:
        clusters = list(self.clusters.values())
        if name_contains:
            clusters = [i for i in clusters if name_contains in i.cluster_name]
        if creation_time_before:
            clusters = [
                i for i in clusters if i.creation_time < str(creation_time_before)
            ]
        if creation_time_after:
            clusters = [
                i for i in clusters if i.creation_time > str(creation_time_after)
            ]
        reverse = sort_order == "Descending"
        if sort_by == "Name":
            clusters = sorted(clusters, key=lambda x: x.cluster_name, reverse=reverse)
        if sort_by == "CreationTime" or sort_by is None:
            clusters = sorted(clusters, key=lambda x: x.creation_time, reverse=reverse)
        return clusters

    @paginate(pagination_model=PAGINATION_MODEL)
    def list_cluster_nodes(
        self,
        cluster_name: str,
        creation_time_after: Optional[str],
        creation_time_before: Optional[str],
        instance_group_name_contains: Optional[str],
        sort_by: Optional[str],
        sort_order: Optional[str],
    ) -> List[ClusterNode]:
        if cluster_name.startswith(f"arn:{self.partition}:sagemaker:"):
            cluster_name = (cluster_name.split(":")[-1]).split("/")[-1]
        cluster = self.clusters.get(cluster_name)
        if not cluster:
            raise ValidationError(message=f"Could not find cluster '{cluster_name}'.")
        nodes_list = list(cluster.nodes.values())

        if instance_group_name_contains:
            nodes_list = [
                i
                for i in nodes_list
                if instance_group_name_contains in i.instance_group_name
            ]
        if creation_time_before:
            nodes_list = [
                i for i in nodes_list if i.launch_time < str(creation_time_before)
            ]
        if creation_time_after:
            nodes_list = [
                i for i in nodes_list if i.launch_time > str(creation_time_after)
            ]
        reverse = sort_order == "Descending"
        if sort_by == "Name":
            nodes_list = sorted(
                nodes_list, key=lambda x: x.instance_group_name, reverse=reverse
            )
        if sort_by == "CreationTime" or sort_by is None:
            nodes_list = sorted(
                nodes_list, key=lambda x: x.launch_time, reverse=reverse
            )
        return nodes_list

    def create_model_bias_job_definition(
        self,
        account_id: str,
        job_definition_name: str,
        tags: List[Dict[str, str]] = [],
        role_arn: str = "",
        job_resources: Optional[Dict[str, Any]] = None,
        stopping_condition: Optional[Dict[str, Any]] = None,
        environment: Optional[Dict[str, str]] = None,
        network_config: Optional[Dict[str, Any]] = None,
        model_bias_baseline_config: Optional[Dict[str, Any]] = None,
        model_bias_app_specification: Optional[Dict[str, Any]] = None,
        model_bias_job_input: Optional[Dict[str, Any]] = None,
        model_bias_job_output_config: Optional[Dict[str, Any]] = None,
    ) -> Dict[str, str]:
        job_definition = FakeModelBiasJobDefinition(
            account_id=account_id,
            region_name=self.region_name,
            job_definition_name=job_definition_name,
            tags=tags,
            role_arn=role_arn,
            job_resources=job_resources,
            stopping_condition=stopping_condition,
            environment=environment,
            network_config=network_config,
            model_bias_baseline_config=model_bias_baseline_config,
            model_bias_app_specification=model_bias_app_specification,
            model_bias_job_input=model_bias_job_input,
            model_bias_job_output_config=model_bias_job_output_config,
        )
        self.model_bias_job_definitions[job_definition_name] = job_definition
        return job_definition.response_create

    @paginate(pagination_model=PAGINATION_MODEL)
    def list_model_bias_job_definitions(self) -> List[Dict[str, str]]:
        return [job.summary_object for job in self.model_bias_job_definitions.values()]

    def describe_model_bias_job_definition(
        self, job_definition_name: str
    ) -> Dict[str, Any]:
        job_definition = self.model_bias_job_definitions.get(job_definition_name)
        if job_definition is None:
            raise ResourceNotFound(f"Job definition {job_definition_name} not found")
        return job_definition.response_object

    def delete_model_bias_job_definition(self, job_definition_name: str) -> None:
        if job_definition_name in self.model_bias_job_definitions:
            del self.model_bias_job_definitions[job_definition_name]
        else:
            raise ResourceNotFound(f"Job definition {job_definition_name} not found")

    def create_auto_ml_job_v2(
        self,
        auto_ml_job_name: str,
        auto_ml_job_input_data_config: List[Dict[str, Any]],
        output_data_config: Dict[str, Any],
        auto_ml_problem_type_config: Dict[str, Any],
        role_arn: str,
        tags: Optional[List[Dict[str, str]]],
        security_config: Optional[Dict[str, Any]],
        auto_ml_job_objective: Optional[Dict[str, str]],
        model_deploy_config: Optional[Dict[str, Any]],
        data_split_config: Optional[Dict[str, Any]],
    ) -> str:
        auto_ml_job = AutoMLJob(
            auto_ml_job_name=auto_ml_job_name,
            auto_ml_job_input_data_config=auto_ml_job_input_data_config,
            output_data_config=output_data_config,
            auto_ml_problem_type_config=auto_ml_problem_type_config,
            role_arn=role_arn,
            region_name=self.region_name,
            account_id=self.account_id,
            tags=tags,
            security_config=security_config,
            auto_ml_job_objective=auto_ml_job_objective,
            model_deploy_config=model_deploy_config,
            data_split_config=data_split_config,
        )

        self.auto_ml_jobs[auto_ml_job_name] = auto_ml_job
        return auto_ml_job.arn

    def describe_auto_ml_job_v2(self, auto_ml_job_name: str) -> Dict[str, Any]:
        if auto_ml_job_name not in self.auto_ml_jobs:
            raise ResourceNotFound(
                f"Could not find AutoML job with name {auto_ml_job_name}."
            )
        auto_ml_job = self.auto_ml_jobs[auto_ml_job_name]
        return auto_ml_job.describe()

    @paginate(pagination_model=PAGINATION_MODEL)
    def list_auto_ml_jobs(
        self,
        creation_time_after: Optional[str],
        creation_time_before: Optional[str],
        last_modified_time_after: Optional[str],
        last_modified_time_before: Optional[str],
        name_contains: Optional[str],
        status_equals: Optional[str],
        sort_order: Optional[str],
        sort_by: Optional[str],
    ) -> List[AutoMLJob]:
        auto_ml_jobs = list(self.auto_ml_jobs.values())
        if name_contains:
            auto_ml_jobs = [
                i for i in auto_ml_jobs if name_contains in i.auto_ml_job_name
            ]
        if status_equals:
            auto_ml_jobs = [
                i for i in auto_ml_jobs if status_equals == i.auto_ml_job_status
            ]
        if creation_time_before:
            auto_ml_jobs = [
                i for i in auto_ml_jobs if i.creation_time < str(creation_time_before)
            ]
        if creation_time_after:
            auto_ml_jobs = [
                i for i in auto_ml_jobs if i.creation_time > str(creation_time_after)
            ]
        if last_modified_time_before:
            auto_ml_jobs = [
                i
                for i in auto_ml_jobs
                if i.last_modified_time < str(last_modified_time_before)
            ]
        if last_modified_time_after:
            auto_ml_jobs = [
                i
                for i in auto_ml_jobs
                if i.last_modified_time > str(last_modified_time_after)
            ]
        reverse = sort_order == "Descending"
        if sort_by == "Status":
            auto_ml_jobs = sorted(
                auto_ml_jobs, key=lambda x: x.auto_ml_job_status, reverse=reverse
            )
        if sort_by == "CreationTime":
            auto_ml_jobs = sorted(
                auto_ml_jobs, key=lambda x: x.creation_time, reverse=reverse
            )
        if sort_by == "Name" or sort_by is None:
            auto_ml_jobs = sorted(
                auto_ml_jobs, key=lambda x: x.auto_ml_job_name, reverse=reverse
            )
        return auto_ml_jobs

    def stop_auto_ml_job(self, auto_ml_job_name: str) -> None:
        if auto_ml_job_name not in self.auto_ml_jobs:
            raise ResourceNotFound(
                f"Could not find AutoML job with name {auto_ml_job_name}."
            )
        auto_ml_job = self.auto_ml_jobs[auto_ml_job_name]
        auto_ml_job.auto_ml_job_status = "Stopped"
        auto_ml_job.auto_ml_job_secondary_status = "Stopped"

    @paginate(pagination_model=PAGINATION_MODEL)
    def list_endpoints(
        self,
        sort_by: Optional[str],
        sort_order: Optional[str],
        name_contains: Optional[str],
        creation_time_before: Optional[str],
        creation_time_after: Optional[str],
        last_modified_time_before: Optional[str],
        last_modified_time_after: Optional[str],
        status_equals: Optional[str],
    ) -> List[FakeEndpoint]:
        endpoints = list(self.endpoints.values())
        if name_contains:
            endpoints = [i for i in endpoints if name_contains in i.endpoint_name]
        if status_equals:
            endpoints = [i for i in endpoints if status_equals == i.endpoint_status]
        if creation_time_before:
            endpoints = [
                i for i in endpoints if i.creation_time < str(creation_time_before)
            ]
        if creation_time_after:
            endpoints = [
                i for i in endpoints if i.creation_time > str(creation_time_after)
            ]
        if last_modified_time_before:
            endpoints = [
                i
                for i in endpoints
                if i.last_modified_time < str(last_modified_time_before)
            ]
        if last_modified_time_after:
            endpoints = [
                i
                for i in endpoints
                if i.last_modified_time > str(last_modified_time_after)
            ]
        reverse = sort_order == "Descending"
        if sort_by == "Name":
            endpoints = sorted(
                endpoints, key=lambda x: x.endpoint_name, reverse=reverse
            )
        elif sort_by == "Status":
            endpoints = sorted(
                endpoints, key=lambda x: x.endpoint_status, reverse=reverse
            )
        else:
            endpoints = sorted(
                endpoints, key=lambda x: x.creation_time, reverse=reverse
            )
        return endpoints

    @paginate(pagination_model=PAGINATION_MODEL)
    def list_endpoint_configs(
        self,
        sort_by: Optional[str],
        sort_order: Optional[str],
        name_contains: Optional[str],
        creation_time_before: Optional[str],
        creation_time_after: Optional[str],
    ) -> List[FakeEndpointConfig]:
        endpoint_configs = list(self.endpoint_configs.values())
        if name_contains:
            endpoint_configs = [
                i for i in endpoint_configs if name_contains in i.endpoint_config_name
            ]
        if creation_time_before:
            endpoint_configs = [
                i
                for i in endpoint_configs
                if i.creation_time < str(creation_time_before)
            ]
        if creation_time_after:
            endpoint_configs = [
                i
                for i in endpoint_configs
                if i.creation_time > str(creation_time_after)
            ]
        reverse = sort_order == "Descending"
        if sort_by == "Name":
            endpoint_configs = sorted(
                endpoint_configs, key=lambda x: x.endpoint_config_name, reverse=reverse
            )
        else:
            endpoint_configs = sorted(
                endpoint_configs, key=lambda x: x.creation_time, reverse=reverse
            )
        return endpoint_configs

    def create_compilation_job(
        self,
        compilation_job_name: str,
        role_arn: str,
        output_config: Dict[str, Any],
        stopping_condition: Dict[str, Any],
        model_package_version_arn: Optional[str],
        input_config: Optional[Dict[str, Any]],
        vpc_config: Optional[Dict[str, Any]],
        tags: Optional[List[Dict[str, str]]],
    ) -> str:
        compilation_job = CompilationJob(
            compilation_job_name=compilation_job_name,
            role_arn=role_arn,
            region_name=self.region_name,
            account_id=self.account_id,
            model_package_version_arn=model_package_version_arn,
            input_config=input_config,
            output_config=output_config,
            vpc_config=vpc_config,
            stopping_condition=stopping_condition,
            tags=tags,
        )
        self.compilation_jobs[compilation_job_name] = compilation_job
        return compilation_job.arn

    def describe_compilation_job(self, compilation_job_name: str) -> Dict[str, Any]:
        if compilation_job_name not in self.compilation_jobs:
            raise ResourceNotFound(
                message=f"Could not find compilation job '{compilation_job_name}'."
            )
        compilation_job = self.compilation_jobs[compilation_job_name]
        return compilation_job.describe()

    @paginate(pagination_model=PAGINATION_MODEL)
    def list_compilation_jobs(
        self,
        creation_time_after: Optional[str],
        creation_time_before: Optional[str],
        last_modified_time_after: Optional[str],
        last_modified_time_before: Optional[str],
        name_contains: Optional[str],
        status_equals: Optional[str],
        sort_by: Optional[str],
        sort_order: Optional[str],
    ) -> List[CompilationJob]:
        compilation_jobs = list(self.compilation_jobs.values())
        if name_contains:
            compilation_jobs = [
                i for i in compilation_jobs if name_contains in i.compilation_job_name
            ]
        if creation_time_before:
            compilation_jobs = [
                i
                for i in compilation_jobs
                if i.creation_time < str(creation_time_before)
            ]
        if creation_time_after:
            compilation_jobs = [
                i
                for i in compilation_jobs
                if i.creation_time > str(creation_time_after)
            ]
        if last_modified_time_before:
            compilation_jobs = [
                i
                for i in compilation_jobs
                if i.last_modified_time < str(last_modified_time_before)
            ]
        if creation_time_after:
            compilation_jobs = [
                i
                for i in compilation_jobs
                if i.last_modified_time > str(last_modified_time_after)
            ]
        if status_equals:
            compilation_jobs = [
                i for i in compilation_jobs if i.compilation_job_status == status_equals
            ]
        reverse = sort_order == "Descending"
        if sort_by == "Name":
            compilation_jobs = sorted(
                compilation_jobs, key=lambda x: x.compilation_job_name, reverse=reverse
            )
        if sort_by == "Status":
            compilation_jobs = sorted(
                compilation_jobs,
                key=lambda x: x.compilation_job_status,
                reverse=reverse,
            )
        if sort_by == "CreationTime" or sort_by is None:
            compilation_jobs = sorted(
                compilation_jobs, key=lambda x: x.creation_time, reverse=reverse
            )
        return compilation_jobs

    def delete_compilation_job(self, compilation_job_name: str) -> None:
        if compilation_job_name not in self.compilation_jobs:
            raise ResourceNotFound(
                message=f"Could not find compilation job '{compilation_job_name}'."
            )
        del self.compilation_jobs[compilation_job_name]

    def create_domain(
        self,
        domain_name: str,
        auth_mode: str,
        default_user_settings: Dict[str, Any],
        subnet_ids: List[str],
        vpc_id: str,
        domain_settings: Optional[Dict[str, Any]],
        tags: Optional[List[Dict[str, str]]],
        app_network_access_type: Optional[str],
        home_efs_file_system_kms_key_id: Optional[str],
        kms_key_id: Optional[str],
        app_security_group_management: Optional[str],
        default_space_settings: Optional[Dict[str, Any]],
    ) -> Dict[str, Any]:
        domain = Domain(
            domain_name=domain_name,
            auth_mode=auth_mode,
            default_user_settings=default_user_settings,
            subnet_ids=subnet_ids,
            vpc_id=vpc_id,
            domain_settings=domain_settings,
            tags=tags,
            app_network_access_type=app_network_access_type,
            home_efs_file_system_kms_key_id=home_efs_file_system_kms_key_id,
            kms_key_id=kms_key_id,
            app_security_group_management=app_security_group_management,
            default_space_settings=default_space_settings,
            region_name=self.region_name,
            account_id=self.account_id,
        )
        self.domains[domain.id] = domain
        return {"DomainArn": domain.arn, "Url": domain.url}

    def describe_domain(self, domain_id: str) -> Dict[str, Any]:
        if domain_id not in self.domains:
            raise ValidationError(message=f"Could not find domain '{domain_id}'.")
        return self.domains[domain_id].describe()

    @paginate(pagination_model=PAGINATION_MODEL)
    def list_domains(self) -> List[Domain]:
        return list(self.domains.values())

    def delete_domain(
        self, domain_id: str, retention_policy: Optional[Dict[str, str]]
    ) -> None:
        # 'retention_policy' parameter is not used
        if domain_id not in self.domains:
            raise ValidationError(message=f"Could not find domain '{domain_id}'.")
        del self.domains[domain_id]

    def create_model_explainability_job_definition(
        self,
        job_definition_name: str,
        model_explainability_baseline_config: Optional[Dict[str, Any]],
        model_explainability_app_specification: Dict[str, Any],
        model_explainability_job_input: Dict[str, Any],
        model_explainability_job_output_config: Dict[str, Any],
        job_resources: Dict[str, Any],
        network_config: Optional[Dict[str, Any]],
        role_arn: str,
        stopping_condition: Optional[Dict[str, Any]],
        tags: List[Dict[str, str]],
    ) -> str:
        model_explainability_job_definition = ModelExplainabilityJobDefinition(
            job_definition_name=job_definition_name,
            model_explainability_baseline_config=model_explainability_baseline_config,
            model_explainability_app_specification=model_explainability_app_specification,
            model_explainability_job_input=model_explainability_job_input,
            model_explainability_job_output_config=model_explainability_job_output_config,
            job_resources=job_resources,
            region_name=self.region_name,
            account_id=self.account_id,
            network_config=network_config,
            role_arn=role_arn,
            stopping_condition=stopping_condition,
            tags=tags,
        )
        self.model_explainability_job_definitions[
            model_explainability_job_definition.job_definition_name
        ] = model_explainability_job_definition
        return model_explainability_job_definition.arn

    def describe_model_explainability_job_definition(
        self, job_definition_name: str
    ) -> Dict[str, Any]:
        if job_definition_name not in self.model_explainability_job_definitions:
            raise ResourceNotFound(
                message=f"Could not find model explainability job definition with name '{job_definition_name}'."
            )
        return self.model_explainability_job_definitions[job_definition_name].describe()

    @paginate(pagination_model=PAGINATION_MODEL)
    def list_model_explainability_job_definitions(
        self,
        endpoint_name: Optional[str],
        sort_by: Optional[str],
        sort_order: Optional[str],
        name_contains: Optional[str],
        creation_time_before: Optional[str],
        creation_time_after: Optional[str],
    ) -> List[ModelExplainabilityJobDefinition]:
        model_explainability_job_definitions = list(
            self.model_explainability_job_definitions.values()
        )
        if endpoint_name:
            model_explainability_job_definitions = [
                i
                for i in model_explainability_job_definitions
                if endpoint_name == i.endpoint_name
            ]
        if name_contains:
            model_explainability_job_definitions = [
                i
                for i in model_explainability_job_definitions
                if name_contains in i.job_definition_name
            ]
        if creation_time_before:
            model_explainability_job_definitions = [
                i
                for i in model_explainability_job_definitions
                if i.creation_time < str(creation_time_before)
            ]
        if creation_time_after:
            model_explainability_job_definitions = [
                i
                for i in model_explainability_job_definitions
                if i.creation_time > str(creation_time_after)
            ]
        reverse = sort_order == "Descending"
        if sort_by == "Name":
            model_explainability_job_definitions = sorted(
                model_explainability_job_definitions,
                key=lambda x: x.job_definition_name,
                reverse=reverse,
            )
        if sort_by == "CreationTime" or sort_by is None:
            model_explainability_job_definitions = sorted(
                model_explainability_job_definitions,
                key=lambda x: x.creation_time,
                reverse=reverse,
            )
        return model_explainability_job_definitions

    def delete_model_explainability_job_definition(
        self, job_definition_name: str
    ) -> None:
        if job_definition_name not in self.model_explainability_job_definitions:
            raise ResourceNotFound(
                message=f"Could not find model explainability job definition with name '{job_definition_name}'."
            )
        del self.model_explainability_job_definitions[job_definition_name]

    def create_hyper_parameter_tuning_job(
        self,
        hyper_parameter_tuning_job_name: str,
        hyper_parameter_tuning_job_config: Dict[str, Any],
        training_job_definition: Optional[Dict[str, Any]],
        training_job_definitions: Optional[List[Dict[str, Any]]],
        warm_start_config: Optional[Dict[str, Any]],
        tags: Optional[List[Dict[str, str]]],
        autotune: Optional[Dict[str, Any]],
    ) -> str:
        hyper_parameter_tuning_job = HyperParameterTuningJob(
            hyper_parameter_tuning_job_name=hyper_parameter_tuning_job_name,
            hyper_parameter_tuning_job_config=hyper_parameter_tuning_job_config,
            region_name=self.region_name,
            account_id=self.account_id,
            training_job_definition=training_job_definition,
            training_job_definitions=training_job_definitions,
            warm_start_config=warm_start_config,
            tags=tags,
            autotune=autotune,
        )

        self.hyper_parameter_tuning_jobs[hyper_parameter_tuning_job_name] = (
            hyper_parameter_tuning_job
        )
        return hyper_parameter_tuning_job.arn

    def describe_hyper_parameter_tuning_job(
        self, hyper_parameter_tuning_job_name: str
    ) -> Dict[str, Any]:
        if hyper_parameter_tuning_job_name not in self.hyper_parameter_tuning_jobs:
            raise ResourceNotFound(
                message=f"Could not find hyper parameter tuning job '{hyper_parameter_tuning_job_name}'."
            )
        return self.hyper_parameter_tuning_jobs[
            hyper_parameter_tuning_job_name
        ].describe()

    @paginate(pagination_model=PAGINATION_MODEL)
    def list_hyper_parameter_tuning_jobs(
        self,
        sort_by: Optional[str],
        sort_order: Optional[str],
        name_contains: Optional[str],
        creation_time_after: Optional[str],
        creation_time_before: Optional[str],
        last_modified_time_after: Optional[str],
        last_modified_time_before: Optional[str],
        status_equals: Optional[str],
    ) -> List[HyperParameterTuningJob]:
        hyper_parameter_tuning_jobs = list(self.hyper_parameter_tuning_jobs.values())
        if name_contains:
            hyper_parameter_tuning_jobs = [
                i
                for i in hyper_parameter_tuning_jobs
                if name_contains in i.hyper_parameter_tuning_job_name
            ]
        if status_equals:
            hyper_parameter_tuning_jobs = [
                i
                for i in hyper_parameter_tuning_jobs
                if status_equals == i.hyper_parameter_tuning_job_status
            ]
        if creation_time_before:
            hyper_parameter_tuning_jobs = [
                i
                for i in hyper_parameter_tuning_jobs
                if i.creation_time < str(creation_time_before)
            ]
        if creation_time_after:
            hyper_parameter_tuning_jobs = [
                i
                for i in hyper_parameter_tuning_jobs
                if i.creation_time > str(creation_time_after)
            ]
        if last_modified_time_before:
            hyper_parameter_tuning_jobs = [
                i
                for i in hyper_parameter_tuning_jobs
                if i.last_modified_time < str(last_modified_time_before)
            ]
        if last_modified_time_after:
            hyper_parameter_tuning_jobs = [
                i
                for i in hyper_parameter_tuning_jobs
                if i.last_modified_time > str(last_modified_time_after)
            ]
        reverse = sort_order == "Descending"
        if sort_by == "Name":
            hyper_parameter_tuning_jobs = sorted(
                hyper_parameter_tuning_jobs,
                key=lambda x: x.hyper_parameter_tuning_job_name,
                reverse=reverse,
            )
        elif sort_by == "Status":
            hyper_parameter_tuning_jobs = sorted(
                hyper_parameter_tuning_jobs,
                key=lambda x: x.hyper_parameter_tuning_job_status,
                reverse=reverse,
            )
        else:
            hyper_parameter_tuning_jobs = sorted(
                hyper_parameter_tuning_jobs,
                key=lambda x: x.creation_time,
                reverse=reverse,
            )
        return hyper_parameter_tuning_jobs

    def delete_hyper_parameter_tuning_job(
        self, hyper_parameter_tuning_job_name: str
    ) -> None:
        if hyper_parameter_tuning_job_name not in self.hyper_parameter_tuning_jobs:
            raise ResourceNotFound(
                message=f"Could not find hyper parameter tuning job '{hyper_parameter_tuning_job_name}'."
            )
        del self.hyper_parameter_tuning_jobs[hyper_parameter_tuning_job_name]

    def create_model_quality_job_definition(
        self,
        job_definition_name: str,
        model_quality_baseline_config: Optional[Dict[str, Any]],
        model_quality_app_specification: Dict[str, Any],
        model_quality_job_input: Dict[str, Any],
        model_quality_job_output_config: Dict[str, Any],
        job_resources: Dict[str, Any],
        network_config: Optional[Dict[str, Any]],
        role_arn: str,
        stopping_condition: Optional[Dict[str, Any]],
        tags: Optional[List[Dict[str, str]]],
    ) -> str:
        model_quality_job_definition = ModelQualityJobDefinition(
            job_definition_name=job_definition_name,
            model_quality_baseline_config=model_quality_baseline_config,
            model_quality_app_specification=model_quality_app_specification,
            model_quality_job_input=model_quality_job_input,
            model_quality_job_output_config=model_quality_job_output_config,
            job_resources=job_resources,
            network_config=network_config,
            role_arn=role_arn,
            stopping_condition=stopping_condition,
            region_name=self.region_name,
            account_id=self.account_id,
            tags=tags,
        )
        self.model_quality_job_definitions[job_definition_name] = (
            model_quality_job_definition
        )
        return model_quality_job_definition.arn

    def describe_model_quality_job_definition(
        self, job_definition_name: str
    ) -> Dict[str, Any]:
        if job_definition_name not in self.model_quality_job_definitions:
            raise ResourceNotFound(
                message=f"Could not find model quality job definition '{job_definition_name}'."
            )
        return self.model_quality_job_definitions[job_definition_name].describe()

    @paginate(pagination_model=PAGINATION_MODEL)
    def list_model_quality_job_definitions(
        self,
        endpoint_name: Optional[str],
        sort_by: Optional[str],
        sort_order: Optional[str],
        name_contains: Optional[str],
        creation_time_before: Optional[str],
        creation_time_after: Optional[str],
    ) -> List[ModelQualityJobDefinition]:
        model_quality_job_definitions = list(
            self.model_quality_job_definitions.values()
        )
        if endpoint_name:
            model_quality_job_definitions = [
                i
                for i in model_quality_job_definitions
                if endpoint_name == i.endpoint_name
            ]
        if name_contains:
            model_quality_job_definitions = [
                i
                for i in model_quality_job_definitions
                if name_contains in i.job_definition_name
            ]
        if creation_time_before:
            model_quality_job_definitions = [
                i
                for i in model_quality_job_definitions
                if i.creation_time < str(creation_time_before)
            ]
        if creation_time_after:
            model_quality_job_definitions = [
                i
                for i in model_quality_job_definitions
                if i.creation_time > str(creation_time_after)
            ]
        reverse = sort_order == "Descending"
        if sort_by == "Name":
            model_quality_job_definitions = sorted(
                model_quality_job_definitions,
                key=lambda x: x.job_definition_name,
                reverse=reverse,
            )
        if sort_by == "CreationTime" or sort_by is None:
            model_quality_job_definitions = sorted(
                model_quality_job_definitions,
                key=lambda x: x.creation_time,
                reverse=reverse,
            )
        return model_quality_job_definitions

    def delete_model_quality_job_definition(self, job_definition_name: str) -> None:
        if job_definition_name not in self.model_quality_job_definitions:
            raise ResourceNotFound(
                message=f"Could not find model quality job definition '{job_definition_name}'."
            )
        del self.model_quality_job_definitions[job_definition_name]

    def create_model_card(
        self,
        model_card_name: str,
        security_config: Optional[Dict[str, str]],
        content: str,
        model_card_status: str,
        tags: Optional[List[Dict[str, str]]],
        model_card_version: Optional[int] = None,
        creation_time: Optional[str] = None,
        last_modified_time: Optional[str] = None,
    ) -> str:
        if model_card_name in self.model_cards:
            raise ConflictException(f"Modelcard {model_card_name} already exists")

        if not model_card_version:
            model_card_version = 1

        # implement here
        model_card = FakeModelCard(
            account_id=self.account_id,
            region_name=self.region_name,
            model_card_name=model_card_name,
            model_card_version=model_card_version,
            content=content,
            model_card_status=model_card_status,
            security_config=security_config,
            tags=tags,
        )

        self.model_cards[model_card_name].append(model_card)
        return model_card.arn

    def update_model_card(
        self, model_card_name: str, content: str, model_card_status: str
    ) -> str:
        if model_card_name not in self.model_cards:
            raise ResourceNotFound(f"Modelcard {model_card_name} does not exist.")

        datetime_now = str(datetime.now(tzutc()))

        first_version = self.model_cards[model_card_name][0]
        creation_time = first_version.creation_time

        most_recent_version = self.model_cards[model_card_name][-1]
        next_version = most_recent_version.model_card_version + 1
        security_config = most_recent_version.security_config
        tags = most_recent_version.tags

        model_card = FakeModelCard(
            account_id=self.account_id,
            region_name=self.region_name,
            model_card_name=model_card_name,
            model_card_version=next_version,
            security_config=security_config,
            content=content,
            model_card_status=model_card_status,
            tags=tags,
            creation_time=creation_time,
            last_modified_time=datetime_now,
        )

        self.model_cards[model_card_name].append(model_card)
        return model_card.arn

    @paginate(pagination_model=PAGINATION_MODEL)
    def list_model_cards(
        self,
        creation_time_after: Optional[datetime],
        creation_time_before: Optional[datetime],
        name_contains: Optional[str],
        model_card_status: Optional[str],
        sort_by: Optional[str],
        sort_order: Optional[str],
    ) -> List[FakeModelCard]:
        model_cards = self.model_cards

        return filter_model_cards(
            model_cards,
            creation_time_after,
            creation_time_before,
            name_contains,
            model_card_status,
            sort_by,
            sort_order,
        )

    @paginate(pagination_model=PAGINATION_MODEL)
    def list_model_card_versions(
        self,
        creation_time_after: Optional[datetime],
        creation_time_before: Optional[datetime],
        model_card_name: str,
        model_card_status: Optional[str],
        sort_by: Optional[str],
        sort_order: Optional[str],
    ) -> List[FakeModelCard]:
        if model_card_name not in self.model_cards:
            raise ResourceNotFound(f"Modelcard {model_card_name} does not exist")

        versions = self.model_cards[model_card_name]
        if creation_time_after:
            versions = [
                v for v in versions if v.last_modified_time > str(creation_time_after)
            ]
        if creation_time_before:
            versions = [
                v for v in versions if v.last_modified_time < str(creation_time_before)
            ]
        if model_card_status:
            versions = [v for v in versions if v.model_card_status == model_card_status]

        reverse = sort_order == "Descending"

        return sorted(versions, key=lambda x: x.model_card_version, reverse=reverse)

    def describe_model_card(
        self, model_card_name: str, model_card_version: int
    ) -> Dict[str, Any]:
        if model_card_name not in self.model_cards:
            raise ResourceNotFound(f"Modelcard {model_card_name} does not exist")

        versions = self.model_cards[model_card_name]
        if model_card_version:
            filtered = [
                v for v in versions if v.model_card_version == model_card_version
            ]
            if filtered:
                version = filtered[0]
                return version.describe()
            else:
                raise ResourceNotFound(
                    f"Modelcard with name {model_card_name} and version: {model_card_version} does not exist"
                )
        return versions[-1].describe()

    def delete_model_card(self, model_card_name: str) -> None:
        if model_card_name not in self.model_cards:
            raise ResourceNotFound(f"Modelcard {model_card_name} does not exist")

        del self.model_cards[model_card_name]

    def create_data_quality_job_definition(
        self,
        account_id: str,
        job_definition_name: str,
        tags: List[Dict[str, str]] = [],
        role_arn: str = "",
        job_resources: Optional[Dict[str, Any]] = None,
        stopping_condition: Optional[Dict[str, Any]] = None,
        environment: Optional[Dict[str, str]] = None,
        network_config: Optional[Dict[str, Any]] = None,
        data_quality_baseline_config: Optional[Dict[str, Any]] = None,
        data_quality_app_specification: Optional[Dict[str, Any]] = None,
        data_quality_job_input: Optional[Dict[str, Any]] = None,
        data_quality_job_output_config: Optional[Dict[str, Any]] = None,
    ) -> Dict[str, str]:
        job_definition = FakeDataQualityJobDefinition(
            account_id=account_id,
            region_name=self.region_name,
            job_definition_name=job_definition_name,
            tags=tags,
            role_arn=role_arn,
            job_resources=job_resources,
            stopping_condition=stopping_condition,
            environment=environment,
            network_config=network_config,
            data_quality_baseline_config=data_quality_baseline_config,
            data_quality_app_specification=data_quality_app_specification,
            data_quality_job_input=data_quality_job_input,
            data_quality_job_output_config=data_quality_job_output_config,
        )
        self.data_quality_job_definitions[job_definition_name] = job_definition
        return job_definition.response_create

    @paginate(pagination_model=PAGINATION_MODEL)
    def list_data_quality_job_definitions(self) -> List[Dict[str, str]]:
        return [
            job.summary_object for job in self.data_quality_job_definitions.values()
        ]

    def describe_data_quality_job_definition(
        self, job_definition_name: str
    ) -> Dict[str, Any]:
        job_definition = self.data_quality_job_definitions.get(job_definition_name)
        if job_definition is None:
            raise ResourceNotFound(f"Job definition {job_definition_name} not found")
        return job_definition.response_object

    def delete_data_quality_job_definition(self, job_definition_name: str) -> None:
        if job_definition_name in self.data_quality_job_definitions:
            del self.data_quality_job_definitions[job_definition_name]
        else:
            raise ResourceNotFound(f"Job definition {job_definition_name} not found")


class FakeDataQualityJobDefinition(BaseObject):
    def __init__(
        self,
        account_id: str,
        region_name: str,
        job_definition_name: str,
        tags: List[Dict[str, str]] = [],
        role_arn: str = "",
        job_resources: Optional[Dict[str, Any]] = None,
        stopping_condition: Optional[Dict[str, Any]] = None,
        environment: Optional[Dict[str, str]] = None,
        network_config: Optional[Dict[str, Any]] = None,
        data_quality_baseline_config: Optional[Dict[str, Any]] = None,
        data_quality_app_specification: Optional[Dict[str, Any]] = None,
        data_quality_job_input: Optional[Dict[str, Any]] = None,
        data_quality_job_output_config: Optional[Dict[str, Any]] = None,
    ):
        self.job_definition_name = job_definition_name
        self.arn = FakeDataQualityJobDefinition.arn_formatter(
            job_definition_name, account_id, region_name
        )
        self.tags = tags
        self.role_arn = role_arn
        self.job_resources = job_resources or {}
        self.stopping_condition = stopping_condition or {}
        self.environment = environment or {}
        self.network_config = network_config or {}
        self.data_quality_baseline_config = data_quality_baseline_config or {}
        self.data_quality_app_specification = data_quality_app_specification or {}
        self.data_quality_job_input = data_quality_job_input or {}
        self.data_quality_job_output_config = data_quality_job_output_config or {}
        self.creation_time = self.last_modified_time = datetime.now().strftime(
            "%Y-%m-%d %H:%M:%S"
        )

    @property
    def response_object(self) -> Dict[str, str]:
        response_object = self.gen_response_object()
        response = {
            k: v for k, v in response_object.items() if v is not None and v != [None]
        }
        response["JobDefinitionArn"] = response.pop("Arn")
        return response

    @property
    def response_create(self) -> Dict[str, str]:
        return {"JobDefinitionArn": self.arn}

    @staticmethod
    def arn_formatter(name: str, account_id: str, region: str) -> str:
        return arn_formatter("data-quality-job-definition", name, account_id, region)

    @property
    def summary_object(self) -> Dict[str, str]:
        return {
            "MonitoringJobDefinitionName": self.job_definition_name,
            "MonitoringJobDefinitionArn": self.arn,
            "CreationTime": self.creation_time,
            "EndpointName": "EndpointName",
        }


class FakeExperiment(BaseObject):
    def __init__(
        self,
        account_id: str,
        region_name: str,
        experiment_name: str,
        tags: List[Dict[str, str]],
    ):
        self.experiment_name = experiment_name
        self.arn = arn_formatter("experiment", experiment_name, account_id, region_name)
        self.tags = tags
        self.creation_time = self.last_modified_time = datetime.now().strftime(
            "%Y-%m-%d %H:%M:%S"
        )

    @property
    def response_object(self) -> Dict[str, Any]:  # type: ignore[misc]
        response_object = self.gen_response_object()
        return {
            k: v for k, v in response_object.items() if v is not None and v != [None]
        }

    @property
    def response_create(self) -> Dict[str, str]:
        return {"ExperimentArn": self.arn}


class FakeTrial(BaseObject):
    def __init__(
        self,
        account_id: str,
        region_name: str,
        trial_name: str,
        experiment_name: str,
        tags: List[Dict[str, str]],
        trial_components: List[str],
    ):
        self.trial_name = trial_name
        self.arn = FakeTrial.arn_formatter(trial_name, account_id, region_name)
        self.tags = tags
        self.trial_components = trial_components
        self.experiment_name = experiment_name
        self.creation_time = self.last_modified_time = datetime.now().strftime(
            "%Y-%m-%d %H:%M:%S"
        )

    @property
    def response_object(self) -> Dict[str, Any]:  # type: ignore[misc]
        response_object = self.gen_response_object()
        response = {
            k: v for k, v in response_object.items() if v is not None and v != [None]
        }
        response["TrialArn"] = response.pop("Arn")

        return response

    @property
    def response_create(self) -> Dict[str, str]:
        return {"TrialArn": self.arn}

    @staticmethod
    def arn_formatter(name: str, account_id: str, region: str) -> str:
        return arn_formatter("experiment-trial", name, account_id, region)


class FakeTrialComponent(BaseObject):
    def __init__(
        self,
        account_id: str,
        region_name: str,
        trial_component_name: str,
        display_name: Optional[str],
        start_time: Optional[datetime],
        end_time: Optional[datetime],
        parameters: Optional[Dict[str, Dict[str, Union[str, float]]]],
        input_artifacts: Optional[Dict[str, Dict[str, str]]],
        output_artifacts: Optional[Dict[str, Dict[str, str]]],
        metadata_properties: Optional[Dict[str, str]],
        status: Optional[Dict[str, str]],
        trial_name: Optional[str],
        tags: List[Dict[str, str]],
    ):
        self.trial_component_name = trial_component_name
        self.display_name = (
            display_name if display_name is not None else trial_component_name
        )
        self.arn = FakeTrialComponent.arn_formatter(
            trial_component_name, account_id, region_name
        )
        self.status = status
        self.tags = tags
        self.trial_name = trial_name
        self.start_time = start_time
        self.end_time = end_time
        now_string = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        self.creation_time = self.last_modified_time = now_string
        self.created_by: Dict[str, Union[Dict[str, str], str]] = {}
        self.last_modified_by: Dict[str, Union[Dict[str, str], str]] = {}
        self.parameters = parameters if parameters is not None else {}
        self.input_artifacts = input_artifacts if input_artifacts is not None else {}
        self.output_artifacts = output_artifacts if output_artifacts is not None else {}
        self.metadata_properties = metadata_properties
        self.metrics: Dict[str, Dict[str, Union[str, int, METRIC_STEP_TYPE]]] = {}
        self.sources: List[Dict[str, str]] = []

    @property
    def response_object(self) -> Dict[str, Any]:  # type: ignore[misc]
        response_object = self.gen_response_object()
        response_object["Metrics"] = self.gen_metrics_response_object()
        response = {
            k: v for k, v in response_object.items() if v is not None and v != [None]
        }
        response["TrialComponentArn"] = response.pop("Arn")

        return response

    def gen_metrics_response_object(
        self,
    ) -> List[Dict[str, Union[str, int, float, datetime]]]:
        metrics_names = self.metrics.keys()
        metrics_response_objects = []
        for metrics_name in metrics_names:
            metrics_steps: METRIC_STEP_TYPE = cast(
                METRIC_STEP_TYPE, self.metrics[metrics_name]["Values"]
            )
            max_step = max(list(metrics_steps.keys()))
            metrics_steps_values: List[float] = list(
                map(
                    lambda metric: cast(float, metric["Value"]),
                    list(metrics_steps.values()),
                )
            )
            count = len(metrics_steps_values)
            mean = sum(metrics_steps_values) / count
            std = (
                sum(map(lambda value: (value - mean) ** 2, metrics_steps_values))
                / count
            ) ** 0.5
            timestamp_int: int = cast(int, self.metrics[metrics_name]["Timestamp"])
            metrics_response_object = {
                "MetricName": metrics_name,
                "SourceArn": self.arn,
                "TimeStamp": datetime.fromtimestamp(timestamp_int, tz=tzutc()).strftime(
                    "%Y-%m-%d %H:%M:%S"
                ),
                "Max": max(metrics_steps_values),
                "Min": min(metrics_steps_values),
                "Last": metrics_steps[max_step]["Value"],
                "Count": count,
                "Avg": mean,
                "StdDev": std,
            }
            metrics_response_objects.append(metrics_response_object)
        return metrics_response_objects

    @property
    def response_create(self) -> Dict[str, str]:
        return {"TrialComponentArn": self.arn}

    @staticmethod
    def arn_formatter(
        trial_component_name: str, account_id: str, region_name: str
    ) -> str:
        return arn_formatter(
            "experiment-trial-component", trial_component_name, account_id, region_name
        )


class FakeModelBiasJobDefinition(BaseObject):
    def __init__(
        self,
        account_id: str,
        region_name: str,
        job_definition_name: str,
        tags: List[Dict[str, str]] = [],
        role_arn: str = "",
        job_resources: Optional[Dict[str, Any]] = None,
        stopping_condition: Optional[Dict[str, Any]] = None,
        environment: Optional[Dict[str, str]] = None,
        network_config: Optional[Dict[str, Any]] = None,
        model_bias_baseline_config: Optional[Dict[str, Any]] = None,
        model_bias_app_specification: Optional[Dict[str, Any]] = None,
        model_bias_job_input: Optional[Dict[str, Any]] = None,
        model_bias_job_output_config: Optional[Dict[str, Any]] = None,
    ):
        self.job_definition_name = job_definition_name
        self.arn = FakeModelBiasJobDefinition.arn_formatter(
            job_definition_name, account_id, region_name
        )
        self.tags = tags
        self.role_arn = role_arn
        self.job_resources = job_resources or {}
        self.stopping_condition = stopping_condition or {}
        self.environment = environment or {}
        self.network_config = network_config or {}
        self.model_bias_baseline_config = model_bias_baseline_config or {}
        self.model_bias_app_specification = model_bias_app_specification or {}
        self.model_bias_job_input = model_bias_job_input or {}
        self.model_bias_job_output_config = model_bias_job_output_config or {}
        self.creation_time = self.last_modified_time = datetime.now().strftime(
            "%Y-%m-%d %H:%M:%S"
        )

    @property
    def response_object(self) -> Dict[str, str]:
        response_object = self.gen_response_object()
        response = {
            k: v for k, v in response_object.items() if v is not None and v != [None]
        }
        response["JobDefinitionArn"] = response.pop("Arn")
        return response

    @property
    def response_create(self) -> Dict[str, str]:
        return {"JobDefinitionArn": self.arn}

    @staticmethod
    def arn_formatter(name: str, account_id: str, region: str) -> str:
        return f"arn:{get_partition(region)}:sagemaker:{region}:{account_id}:model-bias-job-definition/{name}"

    @property
    def summary_object(self) -> Dict[str, str]:
        return {
            "MonitoringJobDefinitionName": self.job_definition_name,
            "MonitoringJobDefinitionArn": self.arn,
            "CreationTime": self.creation_time,
            "EndpointName": self.model_bias_job_input.get("EndpointInput", {}).get(
                "EndpointName", "EndpointName"
            ),
        }


sagemaker_backends = BackendDict(SageMakerModelBackend, "sagemaker")
