"""Handles incoming bedrock requests, invokes methods, returns responses."""

import json
from urllib.parse import unquote

from moto.core.responses import BaseResponse

from .models import BedrockBackend, bedrock_backends


class BedrockResponse(BaseResponse):
    """Handler for Bedrock requests and responses."""

    def __init__(self) -> None:
        super().__init__(service_name="bedrock")

    @property
    def bedrock_backend(self) -> BedrockBackend:
        """Return backend instance specific for this region."""
        return bedrock_backends[self.current_account][self.region]

    def create_model_customization_job(self) -> str:
        params = json.loads(self.body)
        job_name = params.get("jobName")
        custom_model_name = params.get("customModelName")
        role_arn = params.get("roleArn")
        client_request_token = params.get("clientRequestToken")
        base_model_identifier = params.get("baseModelIdentifier")
        customization_type = params.get("customizationType")
        custom_model_kms_key_id = params.get("customModelKmsKeyId")
        job_tags = params.get("jobTags")
        custom_model_tags = params.get("customModelTags")
        training_data_config = params.get("trainingDataConfig")
        validation_data_config = params.get("validationDataConfig")
        output_data_config = params.get("outputDataConfig")
        hyper_parameters = params.get("hyperParameters")
        vpc_config = params.get("vpcConfig")
        job_arn = self.bedrock_backend.create_model_customization_job(
            job_name=job_name,
            custom_model_name=custom_model_name,
            role_arn=role_arn,
            client_request_token=client_request_token,
            base_model_identifier=base_model_identifier,
            customization_type=customization_type,
            custom_model_kms_key_id=custom_model_kms_key_id,
            job_tags=job_tags,
            custom_model_tags=custom_model_tags,
            training_data_config=training_data_config,
            validation_data_config=validation_data_config,
            output_data_config=output_data_config,
            hyper_parameters=hyper_parameters,
            vpc_config=vpc_config,
        )
        return json.dumps(dict(jobArn=job_arn))

    def get_model_customization_job(self) -> str:
        job_identifier = self.path.split("/")[-1]
        model_customization_job = self.bedrock_backend.get_model_customization_job(
            job_identifier=job_identifier
        )
        return json.dumps(dict(model_customization_job.to_dict()))

    def get_model_invocation_logging_configuration(self) -> str:
        logging_config = (
            self.bedrock_backend.get_model_invocation_logging_configuration()
        )
        return json.dumps(dict(loggingConfig=logging_config))

    def put_model_invocation_logging_configuration(self) -> None:
        params = json.loads(self.body)
        logging_config = params.get("loggingConfig")
        self.bedrock_backend.put_model_invocation_logging_configuration(
            logging_config=logging_config
        )
        return

    def tag_resource(self) -> None:
        params = json.loads(self.body)
        resource_arn = params.get("resourceARN")
        tags = params.get("tags")
        self.bedrock_backend.tag_resource(
            resource_arn=resource_arn,
            tags=tags,
        )
        return

    def untag_resource(self) -> str:
        params = json.loads(self.body)
        resource_arn = params.get("resourceARN")
        tag_keys = params.get("tagKeys")
        self.bedrock_backend.untag_resource(
            resource_arn=resource_arn,
            tag_keys=tag_keys,
        )
        return json.dumps(dict())

    def list_tags_for_resource(self) -> str:
        params = json.loads(self.body)
        resource_arn = params.get("resourceARN")
        tags = self.bedrock_backend.list_tags_for_resource(
            resource_arn=resource_arn,
        )
        return json.dumps(dict(tags=tags))

    def get_custom_model(self) -> str:
        model_identifier = unquote(self.path.split("/")[-1])
        custom_model = self.bedrock_backend.get_custom_model(
            model_identifier=model_identifier
        )
        return json.dumps(dict(custom_model.to_dict()))

    def list_custom_models(self) -> str:
        params = self._get_params()
        creation_time_before = params.get("creationTimeBefore")
        creation_time_after = params.get("creationTimeAfter")
        name_contains = params.get("nameContains")
        base_model_arn_equals = params.get("baseModelArnEquals")
        foundation_model_arn_equals = params.get("foundationModelArnEquals")
        max_results = params.get("maxResults")
        next_token = params.get("nextToken")
        sort_by = params.get("sortBy")
        sort_order = params.get("sortOrder")

        max_results = int(max_results) if max_results else None
        model_summaries, next_token = self.bedrock_backend.list_custom_models(
            creation_time_before=creation_time_before,
            creation_time_after=creation_time_after,
            name_contains=name_contains,
            base_model_arn_equals=base_model_arn_equals,
            foundation_model_arn_equals=foundation_model_arn_equals,
            max_results=max_results,
            next_token=next_token,
            sort_by=sort_by,
            sort_order=sort_order,
        )
        return json.dumps(dict(nextToken=next_token, modelSummaries=model_summaries))

    def list_model_customization_jobs(self) -> str:
        params = self._get_params()
        creation_time_after = params.get("creationTimeAfter")
        creation_time_before = params.get("creationTimeBefore")
        status_equals = params.get("statusEquals")
        name_contains = params.get("nameContains")
        max_results = params.get("maxResults")
        next_token = params.get("nextToken")
        sort_by = params.get("sortBy")
        sort_order = params.get("sortOrder")

        max_results = int(max_results) if max_results else None
        (
            model_customization_job_summaries,
            next_token,
        ) = self.bedrock_backend.list_model_customization_jobs(
            creation_time_after=creation_time_after,
            creation_time_before=creation_time_before,
            status_equals=status_equals,
            name_contains=name_contains,
            max_results=max_results,
            next_token=next_token,
            sort_by=sort_by,
            sort_order=sort_order,
        )
        return json.dumps(
            dict(
                nextToken=next_token,
                modelCustomizationJobSummaries=model_customization_job_summaries,
            )
        )

    def delete_custom_model(self) -> str:
        model_identifier = self.path.split("/")[-1]
        self.bedrock_backend.delete_custom_model(
            model_identifier=model_identifier,
        )
        return json.dumps(dict())

    def stop_model_customization_job(self) -> str:
        job_identifier = self.path.split("/")[-2]
        self.bedrock_backend.stop_model_customization_job(
            job_identifier=job_identifier,
        )
        return json.dumps(dict())

    def delete_model_invocation_logging_configuration(self) -> str:
        self.bedrock_backend.delete_model_invocation_logging_configuration()
        return json.dumps(dict())
