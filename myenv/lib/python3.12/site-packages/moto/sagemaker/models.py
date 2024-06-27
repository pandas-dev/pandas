import json
import os
import random
import re
import string
from datetime import datetime
from typing import Any, Dict, Iterable, List, Optional, Union, cast

from dateutil.tz import tzutc

from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import BaseModel, CloudFormationModel
from moto.sagemaker import validators
from moto.utilities.paginator import paginate
from moto.utilities.utils import ARN_PARTITION_REGEX, get_partition

from .exceptions import (
    AWSValidationException,
    MissingModel,
    ResourceInUseException,
    ResourceNotFound,
    ValidationError,
)
from .utils import (
    arn_formatter,
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
        "unique_attribute": "experiment_arn",
    },
    "list_trials": {
        "input_token": "NextToken",
        "limit_key": "MaxResults",
        "limit_default": 100,
        "unique_attribute": "trial_arn",
    },
    "list_trial_components": {
        "input_token": "NextToken",
        "limit_key": "MaxResults",
        "limit_default": 100,
        "unique_attribute": "trial_component_arn",
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
        "unique_attribute": "model_package_group_arn",
    },
    "list_model_packages": {
        "input_token": "next_token",
        "limit_key": "max_results",
        "limit_default": 100,
        "unique_attribute": "model_package_arn",
    },
    "list_notebook_instances": {
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
        self.pipeline_execution_arn = pipeline_execution_arn
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
        self.pipeline_arn = arn_formatter(
            "pipeline", pipeline_name, account_id, region_name
        )
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
        self.processing_job_arn = FakeProcessingJob.arn_formatter(
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
        return {
            k: v for k, v in response_object.items() if v is not None and v != [None]
        }

    @property
    def response_create(self) -> Dict[str, str]:
        return {"ProcessingJobArn": self.processing_job_arn}

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
        self.training_job_arn = FakeTrainingJob.arn_formatter(
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
        return {
            k: v for k, v in response_object.items() if v is not None and v != [None]
        }

    @property
    def response_create(self) -> Dict[str, str]:
        return {"TrainingJobArn": self.training_job_arn}

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
        self.endpoint_arn = FakeEndpoint.arn_formatter(
            endpoint_name, account_id, region_name
        )
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

    @property
    def response_object(self) -> Dict[str, Any]:  # type: ignore[misc]
        response_object = self.gen_response_object()
        return {
            k: v for k, v in response_object.items() if v is not None and v != [None]
        }

    @property
    def response_create(self) -> Dict[str, str]:
        return {"EndpointArn": self.endpoint_arn}

    @staticmethod
    def arn_formatter(endpoint_name: str, account_id: str, region_name: str) -> str:
        return arn_formatter("endpoint", endpoint_name, account_id, region_name)

    @property
    def physical_resource_id(self) -> str:
        return self.endpoint_arn

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
            original_resource.endpoint_arn, cloudformation_json, account_id, region_name
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
        self.transform_job_arn = FakeTransformJob.arn_formatter(
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
        return {"TransformJobArn": self.transform_job_arn}

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
        self.model_arn = arn_formatter(
            "model", self.model_name, account_id, region_name
        )

    @property
    def response_object(self) -> Dict[str, Any]:  # type: ignore[misc]
        response_object = self.gen_response_object()
        return {
            k: v for k, v in response_object.items() if v is not None and v != [None]
        }

    @property
    def response_create(self) -> Dict[str, str]:
        return {"ModelArn": self.model_arn}

    @property
    def physical_resource_id(self) -> str:
        return self.model_arn

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
            original_resource.model_arn, cloudformation_json, account_id, region_name
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
        self.model_package_group_arn = model_package_group_arn
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
            "ModelPackageGroupArn",
            "ModelPackageGroupDescription",
            "CreationTime",
            "ModelPackageGroupStatus",
        ]
        return {k: v for k, v in response_object.items() if k in response_values}


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
        self.feature_group_arn = arn_formatter(
            region_name=region_name,
            account_id=account_id,
            _type="feature-group",
            _id=f"{self.feature_group_name.lower()}",
        )
        self.tags = tags

    def describe(self) -> Dict[str, Any]:
        return {
            "FeatureGroupArn": self.feature_group_arn,
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
        self.model_package_arn = model_package_arn
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
            "ModelPackageArn",
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
        return {
            k: v
            for k, v in response_object.items()
            if k in response_values
            if v is not None
        }

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
        self.notebook_instance_lifecycle_config_arn = (
            FakeSageMakerNotebookInstanceLifecycleConfig.arn_formatter(
                self.notebook_instance_lifecycle_config_name, account_id, region_name
            )
        )

    @staticmethod
    def arn_formatter(name: str, account_id: str, region_name: str) -> str:
        return arn_formatter(
            "notebook-instance-lifecycle-configuration", name, account_id, region_name
        )

    @property
    def response_object(self) -> Dict[str, Any]:  # type: ignore[misc]
        response_object = self.gen_response_object()
        return {
            k: v for k, v in response_object.items() if v is not None and v != [None]
        }

    @property
    def physical_resource_id(self) -> str:
        return self.notebook_instance_lifecycle_config_arn

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
            original_resource.notebook_instance_lifecycle_config_arn,
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
        self.model_package_groups: Dict[str, ModelPackageGroup] = {}
        self.model_packages: Dict[str, ModelPackage] = {}
        self.model_package_name_mapping: Dict[str, str] = {}
        self.feature_groups: Dict[str, FeatureGroup] = {}

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
            "ExperimentArn": experiment_data.experiment_arn,
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
        }
        target_resource, target_name = arn.split(":")[-1].split("/")
        try:
            resource = resources.get(target_resource).get(target_name)  # type: ignore
        except KeyError:
            message = f"Could not find {target_resource} with name {target_name}"
            raise ValidationError(message=message)
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
        next_index = None

        valid_resources = [
            "Pipeline",
            "ModelPackageGroup",
            "TrainingJob",
            "ExperimentTrialComponent",
            "FeatureGroup",
            "Endpoint",
            "PipelineExecution",
            "Project",
            "ExperimentTrial",
            "Image",
            "ImageVersion",
            "ModelPackage",
            "Experiment",
        ]

        if resource not in valid_resources:
            raise AWSValidationException(
                f"An error occurred (ValidationException) when calling the Search operation: 1 validation error detected: Value '{resource}' at 'resource' failed to satisfy constraint: Member must satisfy enum value set: {valid_resources}"
            )

        def evaluate_search_expression(item: Any) -> bool:
            filters = None
            if search_expression is not None:
                filters = search_expression.get("Filters")

            if filters is not None:
                for f in filters:
                    if f["Operator"] == "Equals":
                        if f["Name"].startswith("Tags."):
                            key = f["Name"][5:]
                            value = f["Value"]

                            if (
                                len(
                                    [
                                        e
                                        for e in item.tags
                                        if e["Key"] == key and e["Value"] == value
                                    ]
                                )
                                == 0
                            ):
                                return False
                        if f["Name"] == "ExperimentName":
                            experiment_name = f["Value"]

                            if hasattr(item, "experiment_name"):
                                if getattr(item, "experiment_name") != experiment_name:
                                    return False
                            else:
                                raise ValidationError(
                                    message="Unknown property name: ExperimentName"
                                )

                        if f["Name"] == "TrialName":
                            raise AWSValidationException(
                                f"An error occurred (ValidationException) when calling the Search operation: Unknown property name: {f['Name']}"
                            )

                        if f["Name"] == "Parents.TrialName":
                            trial_name = f["Value"]

                            if getattr(item, "trial_name") != trial_name:
                                return False

            return True

        result: Dict[str, Any] = {
            "Results": [],
            "NextToken": str(next_index) if next_index is not None else None,
        }
        if resource == "Experiment":
            experiments_fetched = list(self.experiments.values())

            experiment_summaries = [
                {
                    "ExperimentName": experiment_data.experiment_name,
                    "ExperimentArn": experiment_data.experiment_arn,
                    "CreationTime": experiment_data.creation_time,
                    "LastModifiedTime": experiment_data.last_modified_time,
                }
                for experiment_data in experiments_fetched
                if evaluate_search_expression(experiment_data)
            ]

            for experiment_summary in experiment_summaries:
                result["Results"].append({"Experiment": experiment_summary})

        if resource == "ExperimentTrial":
            trials_fetched = list(self.trials.values())

            trial_summaries = [
                {
                    "TrialName": trial_data.trial_name,
                    "TrialArn": trial_data.trial_arn,
                    "CreationTime": trial_data.creation_time,
                    "LastModifiedTime": trial_data.last_modified_time,
                }
                for trial_data in trials_fetched
                if evaluate_search_expression(trial_data)
            ]

            for trial_summary in trial_summaries:
                result["Results"].append({"Trial": trial_summary})

        if resource == "ExperimentTrialComponent":
            trial_components_fetched = list(self.trial_components.values())

            trial_component_summaries = [
                {
                    "TrialComponentName": trial_component_data.trial_component_name,
                    "TrialComponentArn": trial_component_data.trial_component_arn,
                    "CreationTime": trial_component_data.creation_time,
                    "LastModifiedTime": trial_component_data.last_modified_time,
                }
                for trial_component_data in trial_components_fetched
                if evaluate_search_expression(trial_component_data)
            ]

            for trial_component_summary in trial_component_summaries:
                result["Results"].append({"TrialComponent": trial_component_summary})
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
            "TrialComponentArn": self.trial_components[
                trial_component_name
            ].trial_component_arn,
            "TrialArn": self.trials[trial_name].trial_arn,
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
        return pipeline.pipeline_arn

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

        return pipeline.pipeline_arn

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
                    "PipelineExecutionArn": pipeline_execution_arn,
                    "StartTime": pipeline_execution.start_time,
                    "PipelineExecutionStatus": pipeline_execution.pipeline_execution_status,
                    "PipelineExecutionDescription": pipeline_execution.pipeline_execution_description,
                    "PipelineExecutionDisplayName": pipeline_execution.pipeline_execution_display_name,
                    "PipelineExecutionFailureReason": str(
                        pipeline_execution.pipeline_execution_failure_reason
                    ),
                }
                for pipeline_execution_arn, pipeline_execution in pipeline.pipeline_executions.items()
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
            "PipelineArn": pipeline.pipeline_arn,
            "PipelineExecutionArn": pipeline_execution.pipeline_execution_arn,
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
            "PipelineArn": pipeline.pipeline_arn,
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
                else datetime.fromtimestamp(x).strftime("%Y-%m-%d " "%H:%M:%S")
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
                "PipelineArn": pipeline_data.pipeline_arn,
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
                "ProcessingJobArn": processing_job_data.processing_job_arn,
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
                "TransformJobArn": transform_job_data.transform_job_arn,
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
                "TrainingJobArn": training_job_data.training_job_arn,
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
        return endpoint.endpoint_arn

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
        return self.model_package_groups[
            model_package_group_name
        ].model_package_group_arn

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

        return model_package.model_package_arn

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
            model_package.model_package_arn
        )
        self.model_package_name_mapping[model_package.model_package_arn] = (
            model_package.model_package_arn
        )
        self.model_packages[model_package.model_package_arn] = model_package
        return model_package.model_package_arn

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
        self.feature_groups[feature_group.feature_group_arn] = feature_group
        return feature_group.feature_group_arn

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


class FakeExperiment(BaseObject):
    def __init__(
        self,
        account_id: str,
        region_name: str,
        experiment_name: str,
        tags: List[Dict[str, str]],
    ):
        self.experiment_name = experiment_name
        self.experiment_arn = arn_formatter(
            "experiment", experiment_name, account_id, region_name
        )
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
        return {"ExperimentArn": self.experiment_arn}


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
        self.trial_arn = FakeTrial.arn_formatter(trial_name, account_id, region_name)
        self.tags = tags
        self.trial_components = trial_components
        self.experiment_name = experiment_name
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
        return {"TrialArn": self.trial_arn}

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
        self.trial_component_arn = FakeTrialComponent.arn_formatter(
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
        return {
            k: v for k, v in response_object.items() if v is not None and v != [None]
        }

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
                "SourceArn": self.trial_component_arn,
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
        return {"TrialComponentArn": self.trial_component_arn}

    @staticmethod
    def arn_formatter(
        trial_component_name: str, account_id: str, region_name: str
    ) -> str:
        return arn_formatter(
            "experiment-trial-component", trial_component_name, account_id, region_name
        )


sagemaker_backends = BackendDict(SageMakerModelBackend, "sagemaker")
