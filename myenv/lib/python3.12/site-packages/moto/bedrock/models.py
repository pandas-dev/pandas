"""BedrockBackend class with methods for supported APIs."""

import re
from datetime import datetime
from typing import Any, Dict, List, Optional

from moto.bedrock.exceptions import (
    ResourceInUseException,
    ResourceNotFoundException,
    TooManyTagsException,
    ValidationException,
)
from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import BaseModel
from moto.utilities.paginator import paginate
from moto.utilities.tagging_service import TaggingService
from moto.utilities.utils import get_partition


class ModelCustomizationJob(BaseModel):
    def __init__(
        self,
        job_name: str,
        custom_model_name: str,
        role_arn: str,
        base_model_identifier: str,
        training_data_config: Dict[str, str],
        output_data_config: Dict[str, str],
        hyper_parameters: Dict[str, str],
        region_name: str,
        account_id: str,
        client_request_token: Optional[str],
        customization_type: Optional[str],
        custom_model_kms_key_id: Optional[str],
        job_tags: Optional[List[Dict[str, str]]],
        custom_model_tags: Optional[List[Dict[str, str]]],
        validation_data_config: Optional[Dict[str, Any]],
        vpc_config: Optional[Dict[str, Any]],
    ):
        self.job_name = job_name
        self.custom_model_name = custom_model_name
        self.role_arn = role_arn
        self.client_request_token = client_request_token
        self.base_model_identifier = base_model_identifier
        self.customization_type = customization_type
        self.custom_model_kms_key_id = custom_model_kms_key_id
        self.job_tags = job_tags
        self.custom_model_tags = custom_model_tags
        if "s3Uri" not in training_data_config or not re.match(
            r"s3://.*", training_data_config["s3Uri"]
        ):
            raise ValidationException(
                "Validation error detected: "
                f"Value '{training_data_config}' at 'training_data_config' failed to satisfy constraint: "
                "Member must satisfy regular expression pattern: "
                "s3://.*"
            )
        self.training_data_config = training_data_config
        if validation_data_config:
            if "validators" in validation_data_config:
                for validator in validation_data_config["validators"]:
                    if not re.match(r"s3://.*", validator["s3Uri"]):
                        raise ValidationException(
                            "Validation error detected: "
                            f"Value '{validator}' at 'validation_data_config' failed to satisfy constraint: "
                            "Member must satisfy regular expression pattern: "
                            "s3://.*"
                        )
        self.validation_data_config = validation_data_config
        if "s3Uri" not in output_data_config or not re.match(
            r"s3://.*", output_data_config["s3Uri"]
        ):
            raise ValidationException(
                "Validation error detected: "
                f"Value '{output_data_config}' at 'output_data_config' failed to satisfy constraint: "
                "Member must satisfy regular expression pattern: "
                "s3://.*"
            )
        self.output_data_config = output_data_config
        self.hyper_parameters = hyper_parameters
        self.vpc_config = vpc_config
        self.region_name = region_name
        self.account_id = account_id
        self.job_arn = f"arn:{get_partition(self.region_name)}:bedrock:{self.region_name}:{self.account_id}:model-customization-job/{self.job_name}"
        self.output_model_name = f"{self.custom_model_name}-{self.job_name}"
        self.output_model_arn = f"arn:{get_partition(self.region_name)}:bedrock:{self.region_name}:{self.account_id}:custom-model/{self.output_model_name}"
        self.status = "InProgress"
        self.failure_message = "Failure Message"
        self.creation_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        self.last_modified_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        self.end_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        self.base_model_arn = f"arn:{get_partition(self.region_name)}:bedrock:{self.region_name}::foundation-model/{self.base_model_identifier}"
        self.output_model_kms_key_arn = f"arn:{get_partition(self.region_name)}:kms:{self.region_name}:{self.account_id}:key/{self.output_model_name}-kms-key"
        self.training_metrics = {"trainingLoss": 0.0}  # hard coded
        self.validation_metrics = [{"validationLoss": 0.0}]  # hard coded

    def to_dict(self) -> Dict[str, Any]:
        dct = {
            "baseModelArn": self.base_model_arn,
            "clientRequestToken": self.client_request_token,
            "creationTime": self.creation_time,
            "customizationType": self.customization_type,
            "endTime": self.end_time,
            "failureMessage": self.failure_message,
            "hyperParameters": self.hyper_parameters,
            "jobArn": self.job_arn,
            "jobName": self.job_name,
            "lastModifiedTime": self.last_modified_time,
            "outputDataConfig": self.output_data_config,
            "outputModelArn": self.output_model_arn,
            "outputModelKmsKeyArn": self.output_model_kms_key_arn,
            "outputModelName": self.output_model_name,
            "roleArn": self.role_arn,
            "status": self.status,
            "trainingDataConfig": self.training_data_config,
            "trainingMetrics": self.training_metrics,
            "validationDataConfig": self.validation_data_config,
            "validationMetrics": self.validation_metrics,
            "vpcConfig": self.vpc_config,
        }
        return {k: v for k, v in dct.items() if v}


class CustomModel(BaseModel):
    def __init__(
        self,
        model_name: str,
        job_name: str,
        job_arn: str,
        base_model_arn: str,
        hyper_parameters: Dict[str, str],
        output_data_config: Dict[str, str],
        training_data_config: Dict[str, str],
        training_metrics: Dict[str, float],
        base_model_name: str,
        region_name: str,
        account_id: str,
        customization_type: Optional[str],
        model_kms_key_arn: Optional[str],
        validation_data_config: Optional[Dict[str, Any]],
        validation_metrics: Optional[List[Dict[str, float]]],
    ):
        self.model_name = model_name
        self.job_name = job_name
        self.job_arn = job_arn
        self.base_model_arn = base_model_arn
        self.customization_type = customization_type
        self.model_kms_key_arn = model_kms_key_arn
        self.hyper_parameters = hyper_parameters
        self.training_data_config = training_data_config
        self.validation_data_config = validation_data_config
        self.output_data_config = output_data_config
        self.training_metrics = training_metrics
        self.validation_metrics = validation_metrics
        self.region_name = region_name
        self.account_id = account_id
        self.model_arn = f"arn:{get_partition(self.region_name)}:bedrock:{self.region_name}:{self.account_id}:custom-model/{self.model_name}"
        self.creation_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        self.base_model_name = base_model_name

    def to_dict(self) -> Dict[str, Any]:
        dct = {
            "baseModelArn": self.base_model_arn,
            "creationTime": self.creation_time,
            "customizationType": self.customization_type,
            "hyperParameters": self.hyper_parameters,
            "jobArn": self.job_arn,
            "jobName": self.job_name,
            "modelArn": self.model_arn,
            "modelKmsKeyArn": self.model_kms_key_arn,
            "modelName": self.model_name,
            "outputDataConfig": self.output_data_config,
            "trainingDataConfig": self.training_data_config,
            "trainingMetrics": self.training_metrics,
            "validationDataConfig": self.validation_data_config,
            "validationMetrics": self.validation_metrics,
        }
        return {k: v for k, v in dct.items() if v}


class model_invocation_logging_configuration(BaseModel):
    def __init__(self, logging_config: Dict[str, Any]) -> None:
        self.logging_config = logging_config


class BedrockBackend(BaseBackend):
    """Implementation of Bedrock APIs."""

    PAGINATION_MODEL = {
        "list_model_customization_jobs": {
            "input_token": "next_token",
            "limit_key": "max_results",
            "limit_default": 100,
            "unique_attribute": "jobArn",
        },
        "list_custom_models": {
            "input_token": "next_token",
            "limit_key": "max_results",
            "limit_default": 100,
            "unique_attribute": "modelArn",
        },
    }

    def __init__(self, region_name: str, account_id: str) -> None:
        super().__init__(region_name, account_id)
        self.model_customization_jobs: Dict[str, ModelCustomizationJob] = {}
        self.custom_models: Dict[str, CustomModel] = {}
        self.model_invocation_logging_configuration: Optional[
            model_invocation_logging_configuration
        ] = None
        self.tagger = TaggingService()

    def _list_arns(self) -> List[str]:
        return [job.job_arn for job in self.model_customization_jobs.values()] + [
            model.model_arn for model in self.custom_models.values()
        ]

    def create_model_customization_job(
        self,
        job_name: str,
        custom_model_name: str,
        role_arn: str,
        base_model_identifier: str,
        training_data_config: Dict[str, Any],
        output_data_config: Dict[str, str],
        hyper_parameters: Dict[str, str],
        client_request_token: Optional[str],
        customization_type: Optional[str],
        custom_model_kms_key_id: Optional[str],
        job_tags: Optional[List[Dict[str, str]]],
        custom_model_tags: Optional[List[Dict[str, str]]],
        validation_data_config: Optional[Dict[str, Any]],
        vpc_config: Optional[Dict[str, Any]],
    ) -> str:
        if job_name in self.model_customization_jobs.keys():
            raise ResourceInUseException(
                f"Model customization job {job_name} already exists"
            )
        if custom_model_name in self.custom_models.keys():
            raise ResourceInUseException(
                f"Custom model {custom_model_name} already exists"
            )
        model_customization_job = ModelCustomizationJob(
            job_name,
            custom_model_name,
            role_arn,
            base_model_identifier,
            training_data_config,
            output_data_config,
            hyper_parameters,
            self.region_name,
            self.account_id,
            client_request_token,
            customization_type,
            custom_model_kms_key_id,
            job_tags,
            custom_model_tags,
            validation_data_config,
            vpc_config,
        )
        self.model_customization_jobs[job_name] = model_customization_job
        if job_tags:
            self.tag_resource(model_customization_job.job_arn, job_tags)
        # Create associated custom model
        custom_model = CustomModel(
            custom_model_name,
            job_name,
            model_customization_job.job_arn,
            model_customization_job.base_model_arn,
            model_customization_job.hyper_parameters,
            model_customization_job.output_data_config,
            model_customization_job.training_data_config,
            model_customization_job.training_metrics,
            model_customization_job.base_model_identifier,
            self.region_name,
            self.account_id,
            model_customization_job.customization_type,
            model_customization_job.output_model_kms_key_arn,
            model_customization_job.validation_data_config,
            model_customization_job.validation_metrics,
        )
        self.custom_models[custom_model_name] = custom_model
        if custom_model_tags:
            self.tag_resource(custom_model.model_arn, custom_model_tags)
        return model_customization_job.job_arn

    def get_model_customization_job(self, job_identifier: str) -> ModelCustomizationJob:
        if job_identifier not in self.model_customization_jobs:
            raise ResourceNotFoundException(
                f"Model customization job {job_identifier} not found"
            )
        else:
            return self.model_customization_jobs[job_identifier]

    def stop_model_customization_job(self, job_identifier: str) -> None:
        if job_identifier in self.model_customization_jobs:
            self.model_customization_jobs[job_identifier].status = "Stopped"
        else:
            raise ResourceNotFoundException(
                f"Model customization job {job_identifier} not found"
            )
        return

    @paginate(pagination_model=PAGINATION_MODEL)  # type: ignore
    def list_model_customization_jobs(
        self,
        creation_time_after: Optional[datetime],
        creation_time_before: Optional[datetime],
        status_equals: Optional[str],
        name_contains: Optional[str],
        max_results: Optional[int],
        next_token: Optional[str],
        sort_by: Optional[str],
        sort_order: Optional[str],
    ) -> List[Any]:
        customization_jobs_fetched = list(self.model_customization_jobs.values())

        if name_contains is not None:
            customization_jobs_fetched = list(
                filter(
                    lambda x: name_contains in x.job_name,
                    customization_jobs_fetched,
                )
            )

        if creation_time_after is not None:
            customization_jobs_fetched = list(
                filter(
                    lambda x: x.creation_time > str(creation_time_after),
                    customization_jobs_fetched,
                )
            )

        if creation_time_before is not None:
            customization_jobs_fetched = list(
                filter(
                    lambda x: x.creation_time < str(creation_time_before),
                    customization_jobs_fetched,
                )
            )
        if status_equals is not None:
            customization_jobs_fetched = list(
                filter(
                    lambda x: x.status == status_equals,
                    customization_jobs_fetched,
                )
            )

        if sort_by is not None:
            if sort_by == "CreationTime":
                if sort_order is not None and sort_order == "Ascending":
                    customization_jobs_fetched = sorted(
                        customization_jobs_fetched, key=lambda x: x.creation_time
                    )
                elif sort_order is not None and sort_order == "Descending":
                    customization_jobs_fetched = sorted(
                        customization_jobs_fetched,
                        key=lambda x: x.creation_time,
                        reverse=True,
                    )
                else:
                    raise ValidationException(f"Invalid sort order: {sort_order}")
            else:
                raise ValidationException(f"Invalid sort by field: {sort_by}")

        model_customization_job_summaries = []
        for model in customization_jobs_fetched:
            model_customization_job_summaries.append(
                {
                    "jobArn": model.job_arn,
                    "baseModelArn": model.base_model_arn,
                    "jobName": model.job_name,
                    "status": model.status,
                    "lastModifiedTime": model.last_modified_time,
                    "creationTime": model.creation_time,
                    "endTime": model.end_time,
                    "customModelArn": model.output_model_arn,
                    "customModelName": model.custom_model_name,
                    "customizationType": model.customization_type,
                }
            )
        return model_customization_job_summaries

    def get_model_invocation_logging_configuration(self) -> Optional[Dict[str, Any]]:
        if self.model_invocation_logging_configuration:
            return self.model_invocation_logging_configuration.logging_config
        else:
            return {}

    def put_model_invocation_logging_configuration(
        self, logging_config: Dict[str, Any]
    ) -> None:
        invocation_logging = model_invocation_logging_configuration(logging_config)
        self.model_invocation_logging_configuration = invocation_logging
        return

    def get_custom_model(self, model_identifier: str) -> CustomModel:
        if model_identifier[:3] == "arn":
            for model in self.custom_models.values():
                if model.model_arn == model_identifier:
                    return model
            raise ResourceNotFoundException(
                f"Custom model {model_identifier} not found"
            )
        elif model_identifier in self.custom_models:
            return self.custom_models[model_identifier]
        else:
            raise ResourceNotFoundException(
                f"Custom model {model_identifier} not found"
            )

    def delete_custom_model(self, model_identifier: str) -> None:
        if model_identifier in self.custom_models:
            del self.custom_models[model_identifier]
        else:
            raise ResourceNotFoundException(
                f"Custom model {model_identifier} not found"
            )
        return

    @paginate(pagination_model=PAGINATION_MODEL)  # type: ignore
    def list_custom_models(
        self,
        creation_time_before: Optional[datetime],
        creation_time_after: Optional[datetime],
        name_contains: Optional[str],
        base_model_arn_equals: Optional[str],
        foundation_model_arn_equals: Optional[str],
        max_results: Optional[int],
        next_token: Optional[str],
        sort_by: Optional[str],
        sort_order: Optional[str],
    ) -> List[Any]:
        """
        The foundation_model_arn_equals-argument is not yet supported
        """
        custom_models_fetched = list(self.custom_models.values())

        if name_contains is not None:
            custom_models_fetched = list(
                filter(
                    lambda x: name_contains in x.job_name,
                    custom_models_fetched,
                )
            )

        if creation_time_after is not None:
            custom_models_fetched = list(
                filter(
                    lambda x: x.creation_time > str(creation_time_after),
                    custom_models_fetched,
                )
            )

        if creation_time_before is not None:
            custom_models_fetched = list(
                filter(
                    lambda x: x.creation_time < str(creation_time_before),
                    custom_models_fetched,
                )
            )
        if base_model_arn_equals is not None:
            custom_models_fetched = list(
                filter(
                    lambda x: x.base_model_arn == base_model_arn_equals,
                    custom_models_fetched,
                )
            )

        if sort_by is not None:
            if sort_by == "CreationTime":
                if sort_order is not None and sort_order == "Ascending":
                    custom_models_fetched = sorted(
                        custom_models_fetched, key=lambda x: x.creation_time
                    )
                elif sort_order is not None and sort_order == "Descending":
                    custom_models_fetched = sorted(
                        custom_models_fetched,
                        key=lambda x: x.creation_time,
                        reverse=True,
                    )
                else:
                    raise ValidationException(f"Invalid sort order: {sort_order}")
            else:
                raise ValidationException(f"Invalid sort by field: {sort_by}")
        model_summaries = []
        for model in custom_models_fetched:
            model_summaries.append(
                {
                    "modelArn": model.model_arn,
                    "modelName": model.model_name,
                    "creationTime": model.creation_time,
                    "baseModelArn": model.base_model_arn,
                    "baseModelName": model.base_model_name,
                    "jobArn": model.job_arn,
                    "customizationType": model.customization_type,
                }
            )
        return model_summaries

    def tag_resource(self, resource_arn: str, tags: List[Dict[str, str]]) -> None:
        if resource_arn not in self._list_arns():
            raise ResourceNotFoundException(f"Resource {resource_arn} not found")
        fixed_tags = []
        if len(tags) + len(self.tagger.list_tags_for_resource(resource_arn)) > 50:
            raise TooManyTagsException(
                "Member must have length less than or equal to 50"
            )
        for tag_dict in tags:
            fixed_tags.append({"Key": tag_dict["key"], "Value": tag_dict["value"]})
        self.tagger.tag_resource(resource_arn, fixed_tags)
        return

    def untag_resource(self, resource_arn: str, tag_keys: List[str]) -> None:
        if resource_arn not in self._list_arns():
            raise ResourceNotFoundException(f"Resource {resource_arn} not found")
        self.tagger.untag_resource_using_names(resource_arn, tag_keys)
        return

    def list_tags_for_resource(self, resource_arn: str) -> List[Dict[str, str]]:
        if resource_arn not in self._list_arns():
            raise ResourceNotFoundException(f"Resource {resource_arn} not found")
        tags = self.tagger.list_tags_for_resource(resource_arn)
        fixed_tags = []
        for tag_dict in tags["Tags"]:
            fixed_tags.append({"key": tag_dict["Key"], "value": tag_dict["Value"]})
        return fixed_tags

    def delete_model_invocation_logging_configuration(self) -> None:
        if self.model_invocation_logging_configuration:
            self.model_invocation_logging_configuration.logging_config = {}
        return


bedrock_backends = BackendDict(BedrockBackend, "bedrock")
