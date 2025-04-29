import json
from typing import Any

from moto.core.common_types import TYPE_RESPONSE
from moto.core.responses import BaseResponse
from moto.sagemaker.exceptions import AWSValidationException

from .models import SageMakerModelBackend, sagemaker_backends


def format_enum_error(value: str, attribute: str, allowed: Any) -> str:
    return f"Value '{value}' at '{attribute}' failed to satisfy constraint: Member must satisfy enum value set: {allowed}"


class SageMakerResponse(BaseResponse):
    def __init__(self) -> None:
        super().__init__(service_name="sagemaker")

    @property
    def sagemaker_backend(self) -> SageMakerModelBackend:
        return sagemaker_backends[self.current_account][self.region]

    def describe_model(self) -> str:
        model_name = self._get_param("ModelName")
        model = self.sagemaker_backend.describe_model(model_name)
        return json.dumps(model.response_object)

    def create_model(self) -> str:
        model_name = self._get_param("ModelName")
        execution_role_arn = self._get_param("ExecutionRoleArn")
        primary_container = self._get_param("PrimaryContainer")
        vpc_config = self._get_param("VpcConfig")
        containers = self._get_param("Containers")
        tags = self._get_param("Tags")
        model = self.sagemaker_backend.create_model(
            model_name=model_name,
            execution_role_arn=execution_role_arn,
            primary_container=primary_container,
            vpc_config=vpc_config,
            containers=containers,
            tags=tags,
        )
        return json.dumps(model.response_create)

    def delete_model(self) -> str:
        model_name = self._get_param("ModelName")
        self.sagemaker_backend.delete_model(model_name)
        return "{}"

    def list_models(self) -> str:
        models = self.sagemaker_backend.list_models()
        return json.dumps({"Models": [model.response_object for model in models]})

    def create_notebook_instance(self) -> TYPE_RESPONSE:
        sagemaker_notebook = self.sagemaker_backend.create_notebook_instance(
            notebook_instance_name=self._get_param("NotebookInstanceName"),
            instance_type=self._get_param("InstanceType"),
            subnet_id=self._get_param("SubnetId"),
            security_group_ids=self._get_param("SecurityGroupIds"),
            role_arn=self._get_param("RoleArn"),
            kms_key_id=self._get_param("KmsKeyId"),
            tags=self._get_param("Tags"),
            lifecycle_config_name=self._get_param("LifecycleConfigName"),
            direct_internet_access=self._get_param("DirectInternetAccess"),
            volume_size_in_gb=self._get_param("VolumeSizeInGB"),
            accelerator_types=self._get_param("AcceleratorTypes"),
            default_code_repository=self._get_param("DefaultCodeRepository"),
            additional_code_repositories=self._get_param("AdditionalCodeRepositories"),
            root_access=self._get_param("RootAccess"),
        )
        return 200, {}, json.dumps({"NotebookInstanceArn": sagemaker_notebook.arn})

    def describe_notebook_instance(self) -> str:
        notebook_instance_name = self._get_param("NotebookInstanceName")
        notebook_instance = self.sagemaker_backend.get_notebook_instance(
            notebook_instance_name
        )
        return json.dumps(notebook_instance.to_dict())

    def start_notebook_instance(self) -> TYPE_RESPONSE:
        notebook_instance_name = self._get_param("NotebookInstanceName")
        self.sagemaker_backend.start_notebook_instance(notebook_instance_name)
        return 200, {}, json.dumps("{}")

    def stop_notebook_instance(self) -> TYPE_RESPONSE:
        notebook_instance_name = self._get_param("NotebookInstanceName")
        self.sagemaker_backend.stop_notebook_instance(notebook_instance_name)
        return 200, {}, json.dumps("{}")

    def delete_notebook_instance(self) -> TYPE_RESPONSE:
        notebook_instance_name = self._get_param("NotebookInstanceName")
        self.sagemaker_backend.delete_notebook_instance(notebook_instance_name)
        return 200, {}, json.dumps("{}")

    def list_notebook_instances(self) -> str:
        sort_by = self._get_param("SortBy", "Name")
        sort_order = self._get_param("SortOrder", "Ascending")
        name_contains = self._get_param("NameContains")
        status = self._get_param("StatusEquals")
        max_results = self._get_param("MaxResults")
        next_token = self._get_param("NextToken")
        instances, next_token = self.sagemaker_backend.list_notebook_instances(
            sort_by=sort_by,
            sort_order=sort_order,
            name_contains=name_contains,
            status=status,
            max_results=max_results,
            next_token=next_token,
        )
        return json.dumps(
            {
                "NotebookInstances": [i.to_dict() for i in instances],
                "NextToken": next_token,
            }
        )

    def list_tags(self) -> TYPE_RESPONSE:
        arn = self._get_param("ResourceArn")
        max_results = self._get_param("MaxResults")
        next_token = self._get_param("NextToken")
        paged_results, next_token = self.sagemaker_backend.list_tags(
            arn=arn, MaxResults=max_results, NextToken=next_token
        )
        response = {"Tags": paged_results}
        if next_token:
            response["NextToken"] = next_token
        return 200, {}, json.dumps(response)

    def add_tags(self) -> TYPE_RESPONSE:
        arn = self._get_param("ResourceArn")
        tags = self._get_param("Tags")
        tags = self.sagemaker_backend.add_tags(arn, tags)
        return 200, {}, json.dumps({"Tags": tags})

    def delete_tags(self) -> TYPE_RESPONSE:
        arn = self._get_param("ResourceArn")
        tag_keys = self._get_param("TagKeys")
        self.sagemaker_backend.delete_tags(arn, tag_keys)
        return 200, {}, json.dumps({})

    def create_endpoint_config(self) -> TYPE_RESPONSE:
        endpoint_config = self.sagemaker_backend.create_endpoint_config(
            endpoint_config_name=self._get_param("EndpointConfigName"),
            production_variants=self._get_param("ProductionVariants"),
            data_capture_config=self._get_param("DataCaptureConfig"),
            tags=self._get_param("Tags"),
            kms_key_id=self._get_param("KmsKeyId"),
        )
        return (
            200,
            {},
            json.dumps({"EndpointConfigArn": endpoint_config.endpoint_config_arn}),
        )

    def describe_endpoint_config(self) -> str:
        endpoint_config_name = self._get_param("EndpointConfigName")
        response = self.sagemaker_backend.describe_endpoint_config(endpoint_config_name)
        return json.dumps(response)

    def delete_endpoint_config(self) -> TYPE_RESPONSE:
        endpoint_config_name = self._get_param("EndpointConfigName")
        self.sagemaker_backend.delete_endpoint_config(endpoint_config_name)
        return 200, {}, json.dumps("{}")

    def create_endpoint(self) -> TYPE_RESPONSE:
        endpoint = self.sagemaker_backend.create_endpoint(
            endpoint_name=self._get_param("EndpointName"),
            endpoint_config_name=self._get_param("EndpointConfigName"),
            tags=self._get_param("Tags"),
        )
        return 200, {}, json.dumps({"EndpointArn": endpoint.arn})

    def describe_endpoint(self) -> str:
        endpoint_name = self._get_param("EndpointName")
        response = self.sagemaker_backend.describe_endpoint(endpoint_name)
        return json.dumps(response)

    def delete_endpoint(self) -> TYPE_RESPONSE:
        endpoint_name = self._get_param("EndpointName")
        self.sagemaker_backend.delete_endpoint(endpoint_name)
        return 200, {}, json.dumps("{}")

    def create_processing_job(self) -> TYPE_RESPONSE:
        processing_job = self.sagemaker_backend.create_processing_job(
            app_specification=self._get_param("AppSpecification"),
            experiment_config=self._get_param("ExperimentConfig"),
            network_config=self._get_param("NetworkConfig"),
            processing_inputs=self._get_param("ProcessingInputs"),
            processing_job_name=self._get_param("ProcessingJobName"),
            processing_output_config=self._get_param("ProcessingOutputConfig"),
            role_arn=self._get_param("RoleArn"),
            stopping_condition=self._get_param("StoppingCondition"),
            tags=self._get_param("Tags"),
        )
        response = {"ProcessingJobArn": processing_job.arn}
        return 200, {}, json.dumps(response)

    def describe_processing_job(self) -> str:
        processing_job_name = self._get_param("ProcessingJobName")
        response = self.sagemaker_backend.describe_processing_job(processing_job_name)
        return json.dumps(response)

    def create_transform_job(self) -> TYPE_RESPONSE:
        transform_job = self.sagemaker_backend.create_transform_job(
            transform_job_name=self._get_param("TransformJobName"),
            model_name=self._get_param("ModelName"),
            max_concurrent_transforms=self._get_param("MaxConcurrentTransforms"),
            model_client_config=self._get_param("ModelClientConfig"),
            max_payload_in_mb=self._get_param("MaxPayloadInMB"),
            batch_strategy=self._get_param("BatchStrategy"),
            environment=self._get_param("Environment"),
            transform_input=self._get_param("TransformInput"),
            transform_output=self._get_param("TransformOutput"),
            data_capture_config=self._get_param("DataCaptureConfig"),
            transform_resources=self._get_param("TransformResources"),
            data_processing=self._get_param("DataProcessing"),
            tags=self._get_param("Tags"),
            experiment_config=self._get_param("ExperimentConfig"),
        )
        response = {
            "TransformJobArn": transform_job.arn,
        }
        return 200, {}, json.dumps(response)

    def describe_transform_job(self) -> str:
        transform_job_name = self._get_param("TransformJobName")
        response = self.sagemaker_backend.describe_transform_job(
            transform_job_name=transform_job_name
        )
        return json.dumps(response)

    def create_training_job(self) -> TYPE_RESPONSE:
        training_job = self.sagemaker_backend.create_training_job(
            training_job_name=self._get_param("TrainingJobName"),
            hyper_parameters=self._get_param("HyperParameters"),
            algorithm_specification=self._get_param("AlgorithmSpecification"),
            role_arn=self._get_param("RoleArn"),
            input_data_config=self._get_param("InputDataConfig"),
            output_data_config=self._get_param("OutputDataConfig"),
            resource_config=self._get_param("ResourceConfig"),
            vpc_config=self._get_param("VpcConfig"),
            stopping_condition=self._get_param("StoppingCondition"),
            tags=self._get_param("Tags"),
            enable_network_isolation=self._get_param("EnableNetworkIsolation", False),
            enable_inter_container_traffic_encryption=self._get_param(
                "EnableInterContainerTrafficEncryption", False
            ),
            enable_managed_spot_training=self._get_param(
                "EnableManagedSpotTraining", False
            ),
            checkpoint_config=self._get_param("CheckpointConfig"),
            debug_hook_config=self._get_param("DebugHookConfig"),
            debug_rule_configurations=self._get_param("DebugRuleConfigurations"),
            tensor_board_output_config=self._get_param("TensorBoardOutputConfig"),
            experiment_config=self._get_param("ExperimentConfig"),
        )
        response = {
            "TrainingJobArn": training_job.arn,
        }
        return 200, {}, json.dumps(response)

    def describe_training_job(self) -> str:
        training_job_name = self._get_param("TrainingJobName")
        response = self.sagemaker_backend.describe_training_job(training_job_name)
        return json.dumps(response)

    def create_notebook_instance_lifecycle_config(self) -> TYPE_RESPONSE:
        lifecycle_configuration = (
            self.sagemaker_backend.create_notebook_instance_lifecycle_config(
                notebook_instance_lifecycle_config_name=self._get_param(
                    "NotebookInstanceLifecycleConfigName"
                ),
                on_create=self._get_param("OnCreate"),
                on_start=self._get_param("OnStart"),
            )
        )
        response = {
            "NotebookInstanceLifecycleConfigArn": lifecycle_configuration.arn,
        }
        return 200, {}, json.dumps(response)

    def describe_notebook_instance_lifecycle_config(self) -> str:
        response = self.sagemaker_backend.describe_notebook_instance_lifecycle_config(
            notebook_instance_lifecycle_config_name=self._get_param(
                "NotebookInstanceLifecycleConfigName"
            )
        )
        return json.dumps(response)

    def delete_notebook_instance_lifecycle_config(self) -> TYPE_RESPONSE:
        self.sagemaker_backend.delete_notebook_instance_lifecycle_config(
            notebook_instance_lifecycle_config_name=self._get_param(
                "NotebookInstanceLifecycleConfigName"
            )
        )
        return 200, {}, json.dumps("{}")

    def search(self) -> TYPE_RESPONSE:
        response = self.sagemaker_backend.search(
            resource=self._get_param("Resource"),
            search_expression=self._get_param("SearchExpression"),
        )
        return 200, {}, json.dumps(response)

    def list_experiments(self) -> TYPE_RESPONSE:
        MaxResults = self._get_param("MaxResults")
        NextToken = self._get_param("NextToken")

        paged_results, next_token = self.sagemaker_backend.list_experiments(
            MaxResults=MaxResults, NextToken=NextToken
        )

        experiment_summaries = [
            {
                "ExperimentName": experiment_data.experiment_name,
                "ExperimentArn": experiment_data.arn,
                "CreationTime": experiment_data.creation_time,
                "LastModifiedTime": experiment_data.last_modified_time,
            }
            for experiment_data in paged_results
        ]

        response = {
            "ExperimentSummaries": experiment_summaries,
            "NextToken": next_token,
        }

        return 200, {}, json.dumps(response)

    def delete_experiment(self) -> TYPE_RESPONSE:
        self.sagemaker_backend.delete_experiment(
            experiment_name=self._get_param("ExperimentName")
        )
        return 200, {}, json.dumps({})

    def create_experiment(self) -> TYPE_RESPONSE:
        response = self.sagemaker_backend.create_experiment(
            experiment_name=self._get_param("ExperimentName")
        )
        return 200, {}, json.dumps(response)

    def describe_experiment(self) -> TYPE_RESPONSE:
        response = self.sagemaker_backend.describe_experiment(
            experiment_name=self._get_param("ExperimentName")
        )
        return 200, {}, json.dumps(response)

    def list_trials(self) -> TYPE_RESPONSE:
        MaxResults = self._get_param("MaxResults")
        NextToken = self._get_param("NextToken")

        paged_results, next_token = self.sagemaker_backend.list_trials(
            NextToken=NextToken,
            MaxResults=MaxResults,
            experiment_name=self._get_param("ExperimentName"),
            trial_component_name=self._get_param("TrialComponentName"),
        )

        trial_summaries = [
            {
                "TrialName": trial_data.trial_name,
                "TrialArn": trial_data.arn,
                "CreationTime": trial_data.creation_time,
                "LastModifiedTime": trial_data.last_modified_time,
            }
            for trial_data in paged_results
        ]

        response = {"TrialSummaries": trial_summaries, "NextToken": next_token}

        return 200, {}, json.dumps(response)

    def create_trial(self) -> TYPE_RESPONSE:
        response = self.sagemaker_backend.create_trial(
            trial_name=self._get_param("TrialName"),
            experiment_name=self._get_param("ExperimentName"),
        )
        return 200, {}, json.dumps(response)

    def list_trial_components(self) -> TYPE_RESPONSE:
        MaxResults = self._get_param("MaxResults")
        NextToken = self._get_param("NextToken")

        paged_results, next_token = self.sagemaker_backend.list_trial_components(
            NextToken=NextToken,
            MaxResults=MaxResults,
            trial_name=self._get_param("TrialName"),
        )

        trial_component_summaries = [
            {
                "TrialComponentName": trial_component_data.trial_component_name,
                "TrialComponentArn": trial_component_data.arn,
                "CreationTime": trial_component_data.creation_time,
                "LastModifiedTime": trial_component_data.last_modified_time,
            }
            for trial_component_data in paged_results
        ]

        response = {
            "TrialComponentSummaries": trial_component_summaries,
            "NextToken": next_token,
        }

        return 200, {}, json.dumps(response)

    def create_trial_component(self) -> TYPE_RESPONSE:
        response = self.sagemaker_backend.create_trial_component(
            trial_component_name=self._get_param("TrialComponentName"),
            start_time=self._get_param("StartTime"),
            end_time=self._get_param("EndTime"),
            display_name=self._get_param("DisplayName"),
            parameters=self._get_param("Parameters"),
            input_artifacts=self._get_param("InputArtifacts"),
            output_artifacts=self._get_param("OutputArtifacts"),
            metadata_properties=self._get_param("MetadataProperties"),
            status=self._get_param("Status"),
            trial_name=self._get_param("TrialName"),
        )
        return 200, {}, json.dumps(response)

    def describe_trial(self) -> str:
        trial_name = self._get_param("TrialName")
        response = self.sagemaker_backend.describe_trial(trial_name)
        return json.dumps(response)

    def delete_trial(self) -> TYPE_RESPONSE:
        trial_name = self._get_param("TrialName")
        self.sagemaker_backend.delete_trial(trial_name)
        return 200, {}, json.dumps({})

    def delete_trial_component(self) -> TYPE_RESPONSE:
        trial_component_name = self._get_param("TrialComponentName")
        self.sagemaker_backend.delete_trial_component(trial_component_name)
        return 200, {}, json.dumps({})

    def describe_trial_component(self) -> str:
        trial_component_name = self._get_param("TrialComponentName")
        response = self.sagemaker_backend.describe_trial_component(trial_component_name)
        return json.dumps(response)

    def associate_trial_component(self) -> TYPE_RESPONSE:
        trial_name = self._get_param("TrialName")
        trial_component_name = self._get_param("TrialComponentName")
        response = self.sagemaker_backend.associate_trial_component(
            trial_name, trial_component_name
        )
        return 200, {}, json.dumps(response)

    def disassociate_trial_component(self) -> TYPE_RESPONSE:
        trial_component_name = self._get_param("TrialComponentName")
        trial_name = self._get_param("TrialName")
        response = self.sagemaker_backend.disassociate_trial_component(
            trial_name, trial_component_name
        )
        return 200, {}, json.dumps(response)

    def update_trial_component(self) -> TYPE_RESPONSE:
        response = self.sagemaker_backend.update_trial_component(
            trial_component_name=self._get_param("TrialComponentName"),
            status=self._get_param("Status"),
            display_name=self._get_param("DisplayName"),
            start_time=self._get_param("StartTime"),
            end_time=self._get_param("EndTime"),
            parameters=self._get_param("Parameters"),
            parameters_to_remove=self._get_param("ParametersToRemove"),
            input_artifacts=self._get_param("InputArtifacts"),
            input_artifacts_to_remove=self._get_param("InputArtifactsToRemove"),
            output_artifacts=self._get_param("OutputArtifacts"),
            output_artifacts_to_remove=self._get_param("OutputArtifactsToRemove"),
        )
        return 200, {}, json.dumps(response)

    def describe_pipeline(self) -> TYPE_RESPONSE:
        response = self.sagemaker_backend.describe_pipeline(
            self._get_param("PipelineName")
        )
        return 200, {}, json.dumps(response)

    def start_pipeline_execution(self) -> TYPE_RESPONSE:
        response = self.sagemaker_backend.start_pipeline_execution(
            self._get_param("PipelineName"),
            self._get_param("PipelineExecutionDisplayName"),
            self._get_param("PipelineParameters"),
            self._get_param("PipelineExecutionDescription"),
            self._get_param("ParallelismConfiguration"),
            self._get_param("ClientRequestToken"),
        )
        return 200, {}, json.dumps(response)

    def describe_pipeline_execution(self) -> TYPE_RESPONSE:
        response = self.sagemaker_backend.describe_pipeline_execution(
            self._get_param("PipelineExecutionArn")
        )
        return 200, {}, json.dumps(response)

    def describe_pipeline_definition_for_execution(self) -> TYPE_RESPONSE:
        response = self.sagemaker_backend.describe_pipeline_definition_for_execution(
            self._get_param("PipelineExecutionArn")
        )
        return 200, {}, json.dumps(response)

    def list_pipeline_parameters_for_execution(self) -> TYPE_RESPONSE:
        response = self.sagemaker_backend.list_pipeline_parameters_for_execution(
            self._get_param("PipelineExecutionArn")
        )
        return 200, {}, json.dumps(response)

    def list_pipeline_executions(self) -> TYPE_RESPONSE:
        response = self.sagemaker_backend.list_pipeline_executions(
            self._get_param("PipelineName")
        )
        return 200, {}, json.dumps(response)

    def create_pipeline(self) -> TYPE_RESPONSE:
        pipeline = self.sagemaker_backend.create_pipeline(
            pipeline_name=self._get_param("PipelineName"),
            pipeline_display_name=self._get_param("PipelineDisplayName"),
            pipeline_definition=self._get_param("PipelineDefinition"),
            pipeline_definition_s3_location=self._get_param(
                "PipelineDefinitionS3Location"
            ),
            pipeline_description=self._get_param("PipelineDescription"),
            role_arn=self._get_param("RoleArn"),
            tags=self._get_param("Tags"),
            parallelism_configuration=self._get_param("ParallelismConfiguration"),
        )
        response = {
            "PipelineArn": pipeline.arn,
        }

        return 200, {}, json.dumps(response)

    def delete_pipeline(self) -> TYPE_RESPONSE:
        arn = self.sagemaker_backend.delete_pipeline(
            pipeline_name=self._get_param("PipelineName"),
        )
        response = {"PipelineArn": arn}
        return 200, {}, json.dumps(response)

    def update_pipeline(self) -> TYPE_RESPONSE:
        arn = self.sagemaker_backend.update_pipeline(
            pipeline_name=self._get_param("PipelineName"),
            pipeline_display_name=self._get_param("PipelineDisplayName"),
            pipeline_definition=self._get_param("PipelineDefinition"),
            pipeline_definition_s3_location=self._get_param(
                "PipelineDefinitionS3Location"
            ),
            pipeline_description=self._get_param("PipelineDescription"),
            role_arn=self._get_param("RoleArn"),
            parallelism_configuration=self._get_param("ParallelismConfiguration"),
        )

        response = {"PipelineArn": arn}
        return 200, {}, json.dumps(response)

    def list_pipelines(self) -> TYPE_RESPONSE:
        max_results_range = range(1, 101)
        allowed_sort_by = ("Name", "CreationTime")
        allowed_sort_order = ("Ascending", "Descending")

        pipeline_name_prefix = self._get_param("PipelineNamePrefix")
        created_after = self._get_param("CreatedAfter")
        created_before = self._get_param("CreatedBefore")
        sort_by = self._get_param("SortBy", "CreationTime")
        sort_order = self._get_param("SortOrder", "Descending")
        next_token = self._get_param("NextToken")
        max_results = self._get_param("MaxResults")

        errors = []
        if max_results and max_results not in max_results_range:
            errors.append(
                f"Value '{max_results}' at 'maxResults' failed to satisfy constraint: Member must have value less than or equal to {max_results_range[-1]}"
            )

        if sort_by not in allowed_sort_by:
            errors.append(format_enum_error(sort_by, "SortBy", allowed_sort_by))
        if sort_order not in allowed_sort_order:
            errors.append(
                format_enum_error(sort_order, "SortOrder", allowed_sort_order)
            )
        if errors:
            raise AWSValidationException(
                f"{len(errors)} validation errors detected: {';'.join(errors)}"
            )

        response = self.sagemaker_backend.list_pipelines(
            pipeline_name_prefix=pipeline_name_prefix,
            created_after=created_after,
            created_before=created_before,
            next_token=next_token,
            max_results=max_results,
            sort_by=sort_by,
            sort_order=sort_order,
        )

        return 200, {}, json.dumps(response)

    def list_processing_jobs(self) -> TYPE_RESPONSE:
        max_results_range = range(1, 101)
        allowed_sort_by = ["Name", "CreationTime", "Status"]
        allowed_sort_order = ["Ascending", "Descending"]
        allowed_status_equals = [
            "Completed",
            "Stopped",
            "InProgress",
            "Stopping",
            "Failed",
        ]

        max_results = self._get_int_param("MaxResults")
        sort_by = self._get_param("SortBy", "CreationTime")
        sort_order = self._get_param("SortOrder", "Ascending")
        status_equals = self._get_param("StatusEquals")
        next_token = self._get_param("NextToken")
        errors = []
        if max_results and max_results not in max_results_range:
            errors.append(
                f"Value '{max_results}' at 'maxResults' failed to satisfy constraint: Member must have value less than or equal to {max_results_range[-1]}"
            )

        if sort_by not in allowed_sort_by:
            errors.append(format_enum_error(sort_by, "sortBy", allowed_sort_by))
        if sort_order not in allowed_sort_order:
            errors.append(
                format_enum_error(sort_order, "sortOrder", allowed_sort_order)
            )

        if status_equals and status_equals not in allowed_status_equals:
            errors.append(
                format_enum_error(status_equals, "statusEquals", allowed_status_equals)
            )

        if errors != []:
            raise AWSValidationException(
                f"{len(errors)} validation errors detected: {';'.join(errors)}"
            )

        response = self.sagemaker_backend.list_processing_jobs(
            next_token=next_token,
            max_results=max_results,
            creation_time_after=self._get_param("CreationTimeAfter"),
            creation_time_before=self._get_param("CreationTimeBefore"),
            last_modified_time_after=self._get_param("LastModifiedTimeAfter"),
            last_modified_time_before=self._get_param("LastModifiedTimeBefore"),
            name_contains=self._get_param("NameContains"),
            status_equals=status_equals,
        )
        return 200, {}, json.dumps(response)

    def list_transform_jobs(self) -> TYPE_RESPONSE:
        max_results_range = range(1, 101)
        allowed_sort_by = ["Name", "CreationTime", "Status"]
        allowed_sort_order = ["Ascending", "Descending"]
        allowed_status_equals = [
            "Completed",
            "Stopped",
            "InProgress",
            "Stopping",
            "Failed",
        ]

        max_results = self._get_int_param("MaxResults")
        sort_by = self._get_param("SortBy", "CreationTime")
        sort_order = self._get_param("SortOrder", "Ascending")
        status_equals = self._get_param("StatusEquals")
        next_token = self._get_param("NextToken")
        errors = []
        if max_results and max_results not in max_results_range:
            errors.append(
                f"Value '{max_results}' at 'maxResults' failed to satisfy constraint: Member must have value less than or equal to {max_results_range[-1]}"
            )

        if sort_by not in allowed_sort_by:
            errors.append(format_enum_error(sort_by, "sortBy", allowed_sort_by))
        if sort_order not in allowed_sort_order:
            errors.append(
                format_enum_error(sort_order, "sortOrder", allowed_sort_order)
            )

        if status_equals and status_equals not in allowed_status_equals:
            errors.append(
                format_enum_error(status_equals, "statusEquals", allowed_status_equals)
            )

        if errors != []:
            raise AWSValidationException(
                f"{len(errors)} validation errors detected: {';'.join(errors)}"
            )

        response = self.sagemaker_backend.list_transform_jobs(
            next_token=next_token,
            max_results=max_results,
            creation_time_after=self._get_param("CreationTimeAfter"),
            creation_time_before=self._get_param("CreationTimeBefore"),
            last_modified_time_after=self._get_param("LastModifiedTimeAfter"),
            last_modified_time_before=self._get_param("LastModifiedTimeBefore"),
            name_contains=self._get_param("NameContains"),
            status_equals=status_equals,
        )
        return 200, {}, json.dumps(response)

    def list_training_jobs(self) -> TYPE_RESPONSE:
        max_results_range = range(1, 101)
        allowed_sort_by = ["Name", "CreationTime", "Status"]
        allowed_sort_order = ["Ascending", "Descending"]
        allowed_status_equals = [
            "Completed",
            "Stopped",
            "InProgress",
            "Stopping",
            "Failed",
        ]

        max_results = self._get_int_param("MaxResults")
        sort_by = self._get_param("SortBy", "CreationTime")
        sort_order = self._get_param("SortOrder", "Ascending")
        status_equals = self._get_param("StatusEquals")
        next_token = self._get_param("NextToken")
        errors = []
        if max_results and max_results not in max_results_range:
            errors.append(
                f"Value '{max_results}' at 'maxResults' failed to satisfy constraint: Member must have value less than or equal to {max_results_range[-1]}"
            )

        if sort_by not in allowed_sort_by:
            errors.append(format_enum_error(sort_by, "sortBy", allowed_sort_by))
        if sort_order not in allowed_sort_order:
            errors.append(
                format_enum_error(sort_order, "sortOrder", allowed_sort_order)
            )

        if status_equals and status_equals not in allowed_status_equals:
            errors.append(
                format_enum_error(status_equals, "statusEquals", allowed_status_equals)
            )

        if errors != []:
            raise AWSValidationException(
                f"{len(errors)} validation errors detected: {';'.join(errors)}"
            )

        response = self.sagemaker_backend.list_training_jobs(
            next_token=next_token,
            max_results=max_results,
            creation_time_after=self._get_param("CreationTimeAfter"),
            creation_time_before=self._get_param("CreationTimeBefore"),
            last_modified_time_after=self._get_param("LastModifiedTimeAfter"),
            last_modified_time_before=self._get_param("LastModifiedTimeBefore"),
            name_contains=self._get_param("NameContains"),
            status_equals=status_equals,
        )
        return 200, {}, json.dumps(response)

    def update_endpoint_weights_and_capacities(self) -> TYPE_RESPONSE:
        endpoint_name = self._get_param("EndpointName")
        desired_weights_and_capacities = self._get_param("DesiredWeightsAndCapacities")
        endpoint_arn = self.sagemaker_backend.update_endpoint_weights_and_capacities(
            endpoint_name=endpoint_name,
            desired_weights_and_capacities=desired_weights_and_capacities,
        )
        return 200, {}, json.dumps({"EndpointArn": endpoint_arn})

    def list_model_package_groups(self) -> str:
        creation_time_after = self._get_param("CreationTimeAfter")
        creation_time_before = self._get_param("CreationTimeBefore")
        max_results = self._get_param("MaxResults")
        name_contains = self._get_param("NameContains")
        next_token = self._get_param("NextToken")
        sort_by = self._get_param("SortBy")
        sort_order = self._get_param("SortOrder")
        (
            model_package_group_summary_list,
            next_token,
        ) = self.sagemaker_backend.list_model_package_groups(
            creation_time_after=creation_time_after,
            creation_time_before=creation_time_before,
            max_results=max_results,
            name_contains=name_contains,
            next_token=next_token,
            sort_by=sort_by,
            sort_order=sort_order,
        )
        model_package_group_summary_list_response_object = [
            x.gen_response_object() for x in model_package_group_summary_list
        ]
        return json.dumps(
            dict(
                ModelPackageGroupSummaryList=model_package_group_summary_list_response_object,
                NextToken=next_token,
            )
        )

    def list_model_packages(self) -> str:
        creation_time_after = self._get_param("CreationTimeAfter")
        creation_time_before = self._get_param("CreationTimeBefore")
        max_results = self._get_param("MaxResults")
        name_contains = self._get_param("NameContains")
        model_approval_status = self._get_param("ModelApprovalStatus")
        model_package_group_name = self._get_param("ModelPackageGroupName")
        model_package_type = self._get_param("ModelPackageType", "Unversioned")
        next_token = self._get_param("NextToken")
        sort_by = self._get_param("SortBy")
        sort_order = self._get_param("SortOrder")
        (
            model_package_summary_list,
            next_token,
        ) = self.sagemaker_backend.list_model_packages(
            creation_time_after=creation_time_after,
            creation_time_before=creation_time_before,
            max_results=max_results,
            name_contains=name_contains,
            model_approval_status=model_approval_status,
            model_package_group_name=model_package_group_name,
            model_package_type=model_package_type,
            next_token=next_token,
            sort_by=sort_by,
            sort_order=sort_order,
        )
        model_package_summary_list_response_object = [
            x.gen_response_object() for x in model_package_summary_list
        ]
        return json.dumps(
            dict(
                ModelPackageSummaryList=model_package_summary_list_response_object,
                NextToken=next_token,
            )
        )

    def describe_model_package(self) -> str:
        model_package_name = self._get_param("ModelPackageName")
        model_package = self.sagemaker_backend.describe_model_package(
            model_package_name=model_package_name,
        )
        return json.dumps(
            model_package.gen_response_object(),
        )

    def describe_model_package_group(self) -> str:
        model_package_group_name = self._get_param("ModelPackageGroupName")
        model_package_group = self.sagemaker_backend.describe_model_package_group(
            model_package_group_name=model_package_group_name,
        )
        return json.dumps(
            model_package_group.gen_response_object(),
        )

    def update_model_package(self) -> str:
        model_package_arn = self._get_param("ModelPackageArn")
        model_approval_status = self._get_param("ModelApprovalStatus")
        approval_dexcription = self._get_param("ApprovalDescription")
        customer_metadata_properties = self._get_param("CustomerMetadataProperties")
        customer_metadata_properties_to_remove = self._get_param(
            "CustomerMetadataPropertiesToRemove", []
        )
        additional_inference_specification_to_add = self._get_param(
            "AdditionalInferenceSpecificationsToAdd"
        )
        updated_model_package_arn = self.sagemaker_backend.update_model_package(
            model_package_arn=model_package_arn,
            model_approval_status=model_approval_status,
            approval_description=approval_dexcription,
            customer_metadata_properties=customer_metadata_properties,
            customer_metadata_properties_to_remove=customer_metadata_properties_to_remove,
            additional_inference_specifications_to_add=additional_inference_specification_to_add,
        )
        return json.dumps(dict(ModelPackageArn=updated_model_package_arn))

    def create_model_package(self) -> str:
        model_package_name = self._get_param("ModelPackageName")
        model_package_group_name = self._get_param("ModelPackageGroupName")
        model_package_description = self._get_param("ModelPackageDescription")
        inference_specification = self._get_param("InferenceSpecification")
        validation_specification = self._get_param("ValidationSpecification")
        source_algorithm_specification = self._get_param("SourceAlgorithmSpecification")
        certify_for_marketplace = self._get_param("CertifyForMarketplace", False)
        tags = self._get_param("Tags")
        model_approval_status = self._get_param("ModelApprovalStatus")
        metadata_properties = self._get_param("MetadataProperties")
        model_metrics = self._get_param("ModelMetrics")
        client_token = self._get_param("ClientToken")
        customer_metadata_properties = self._get_param("CustomerMetadataProperties")
        drift_check_baselines = self._get_param("DriftCheckBaselines")
        domain = self._get_param("Domain")
        task = self._get_param("Task")
        sample_payload_url = self._get_param("SamplePayloadUrl")
        additional_inference_specifications = self._get_param(
            "AdditionalInferenceSpecifications", None
        )
        model_package_arn = self.sagemaker_backend.create_model_package(
            model_package_name=model_package_name,
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
            client_token=client_token,
        )
        return json.dumps(dict(ModelPackageArn=model_package_arn))

    def create_model_package_group(self) -> str:
        model_package_group_name = self._get_param("ModelPackageGroupName")
        model_package_group_description = self._get_param(
            "ModelPackageGroupDescription"
        )
        tags = self._get_param("Tags")
        model_package_group_arn = self.sagemaker_backend.create_model_package_group(
            model_package_group_name=model_package_group_name,
            model_package_group_description=model_package_group_description,
            tags=tags,
        )
        return json.dumps(dict(ModelPackageGroupArn=model_package_group_arn))

    def create_feature_group(self) -> str:
        feature_group_arn = self.sagemaker_backend.create_feature_group(
            feature_group_name=self._get_param("FeatureGroupName"),
            record_identifier_feature_name=self._get_param(
                "RecordIdentifierFeatureName"
            ),
            event_time_feature_name=self._get_param("EventTimeFeatureName"),
            feature_definitions=self._get_param("FeatureDefinitions"),
            offline_store_config=self._get_param("OfflineStoreConfig"),
            role_arn=self._get_param("RoleArn"),
            tags=self._get_param("Tags"),
        )
        return json.dumps(dict(FeatureGroupArn=feature_group_arn))

    def describe_feature_group(self) -> str:
        resp = self.sagemaker_backend.describe_feature_group(
            feature_group_name=self._get_param("FeatureGroupName"),
        )
        return json.dumps(resp)

    def create_cluster(self) -> str:
        cluster_name = self._get_param("ClusterName")
        instance_groups = self._get_param("InstanceGroups")
        vpc_config = self._get_param("VpcConfig")
        tags = self._get_param("Tags")
        cluster_arn = self.sagemaker_backend.create_cluster(
            cluster_name=cluster_name,
            instance_groups=instance_groups,
            vpc_config=vpc_config,
            tags=tags,
        )
        return json.dumps(dict(ClusterArn=cluster_arn))

    def delete_cluster(self) -> str:
        cluster_name = self._get_param("ClusterName")
        cluster_arn = self.sagemaker_backend.delete_cluster(
            cluster_name=cluster_name,
        )
        return json.dumps(dict(ClusterArn=cluster_arn))

    def describe_cluster(self) -> str:
        cluster_name = self._get_param("ClusterName")
        cluster_description = self.sagemaker_backend.describe_cluster(
            cluster_name=cluster_name,
        )
        return json.dumps(cluster_description)

    def describe_cluster_node(self) -> str:
        cluster_name = self._get_param("ClusterName")
        node_id = self._get_param("NodeId")
        node_details = self.sagemaker_backend.describe_cluster_node(
            cluster_name=cluster_name,
            node_id=node_id,
        )
        return json.dumps(dict(NodeDetails=node_details))

    def list_clusters(self) -> str:
        creation_time_after = self._get_param("CreationTimeAfter")
        creation_time_before = self._get_param("CreationTimeBefore")
        max_results = self._get_param("MaxResults")
        name_contains = self._get_param("NameContains")
        next_token = self._get_param("NextToken")
        sort_by = self._get_param("SortBy")
        sort_order = self._get_param("SortOrder")
        clusters, next_token = self.sagemaker_backend.list_clusters(
            creation_time_after=creation_time_after,
            creation_time_before=creation_time_before,
            max_results=max_results,
            name_contains=name_contains,
            next_token=next_token,
            sort_by=sort_by,
            sort_order=sort_order,
        )
        cluster_summaries = [cluster.summary() for cluster in clusters]
        return json.dumps(
            dict(NextToken=next_token, ClusterSummaries=cluster_summaries)
        )

    def list_cluster_nodes(self) -> str:
        cluster_name = self._get_param("ClusterName")
        creation_time_after = self._get_param("CreationTimeAfter")
        creation_time_before = self._get_param("CreationTimeBefore")
        instance_group_name_contains = self._get_param("InstanceGroupNameContains")
        max_results = self._get_param("MaxResults")
        next_token = self._get_param("NextToken")
        sort_by = self._get_param("SortBy")
        sort_order = self._get_param("SortOrder")
        cluster_nodes, next_token = self.sagemaker_backend.list_cluster_nodes(
            cluster_name=cluster_name,
            creation_time_after=creation_time_after,
            creation_time_before=creation_time_before,
            instance_group_name_contains=instance_group_name_contains,
            max_results=max_results,
            next_token=next_token,
            sort_by=sort_by,
            sort_order=sort_order,
        )
        cluster_node_summaries = [node.summary() for node in cluster_nodes]
        return json.dumps(
            dict(NextToken=next_token, ClusterNodeSummaries=cluster_node_summaries)
        )

    def create_model_bias_job_definition(self) -> str:
        account_id = self.current_account
        job_definition_name = self._get_param("JobDefinitionName")
        tags = self._get_param("Tags", [])
        role_arn = self._get_param("RoleArn")
        job_resources = self._get_param("JobResources")
        stopping_condition = self._get_param("StoppingCondition")
        environment = self._get_param("Environment", {})
        network_config = self._get_param("NetworkConfig", {})
        model_bias_baseline_config = self._get_param("ModelBiasBaselineConfig", {})
        model_bias_app_specification = self._get_param("ModelBiasAppSpecification")
        model_bias_job_input = self._get_param("ModelBiasJobInput")
        model_bias_job_output_config = self._get_param("ModelBiasJobOutputConfig")

        response = self.sagemaker_backend.create_model_bias_job_definition(
            account_id=account_id,
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
        return json.dumps(response)

    def list_model_bias_job_definitions(self) -> str:
        result, next_token = self.sagemaker_backend.list_model_bias_job_definitions()
        return json.dumps({"JobDefinitionSummaries": result, "NextToken": next_token})

    def describe_model_bias_job_definition(self) -> str:
        job_definition_name = self._get_param("JobDefinitionName")
        job_definition = self.sagemaker_backend.describe_model_bias_job_definition(
            job_definition_name
        )
        return json.dumps(job_definition)

    def delete_model_bias_job_definition(self) -> str:
        job_definition_name = self._get_param("JobDefinitionName")
        self.sagemaker_backend.delete_model_bias_job_definition(job_definition_name)
        return json.dumps({})

    def create_data_quality_job_definition(self) -> str:
        account_id = self.current_account
        job_definition_name = self._get_param("JobDefinitionName")
        tags = self._get_param("Tags", [])
        role_arn = self._get_param("RoleArn")
        job_resources = self._get_param("JobResources")
        stopping_condition = self._get_param("StoppingCondition")
        environment = self._get_param("Environment", {})
        network_config = self._get_param("NetworkConfig", {})
        data_quality_baseline_config = self._get_param("DataQualityBaselineConfig", {})
        data_quality_app_specification = self._get_param("DataQualityAppSpecification")
        data_quality_job_input = self._get_param("DataQualityJobInput")
        data_quality_job_output_config = self._get_param("DataQualityJobOutputConfig")

        response = self.sagemaker_backend.create_data_quality_job_definition(
            account_id=account_id,
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
        return json.dumps(response)

    def list_data_quality_job_definitions(self) -> str:
        result, next_token = self.sagemaker_backend.list_data_quality_job_definitions()
        return json.dumps({"JobDefinitionSummaries": result, "NextToken": next_token})

    def describe_data_quality_job_definition(self) -> str:
        job_definition_name = self._get_param("JobDefinitionName")
        job_definition = self.sagemaker_backend.describe_data_quality_job_definition(
            job_definition_name
        )
        return json.dumps(job_definition)

    def delete_data_quality_job_definition(self) -> str:
        job_definition_name = self._get_param("JobDefinitionName")
        self.sagemaker_backend.delete_data_quality_job_definition(job_definition_name)
        return json.dumps({})

    def create_auto_ml_job_v2(self) -> str:
        auto_ml_job_name = self._get_param("AutoMLJobName")
        auto_ml_job_input_data_config = self._get_param("AutoMLJobInputDataConfig")
        output_data_config = self._get_param("OutputDataConfig")
        auto_ml_problem_type_config = self._get_param("AutoMLProblemTypeConfig")
        role_arn = self._get_param("RoleArn")
        tags = self._get_param("Tags")
        security_config = self._get_param("SecurityConfig")
        auto_ml_job_objective = self._get_param("AutoMLJobObjective")
        model_deploy_config = self._get_param("ModelDeployConfig")
        data_split_config = self._get_param("DataSplitConfig")
        auto_ml_job_arn = self.sagemaker_backend.create_auto_ml_job_v2(
            auto_ml_job_name=auto_ml_job_name,
            auto_ml_job_input_data_config=auto_ml_job_input_data_config,
            output_data_config=output_data_config,
            auto_ml_problem_type_config=auto_ml_problem_type_config,
            role_arn=role_arn,
            tags=tags,
            security_config=security_config,
            auto_ml_job_objective=auto_ml_job_objective,
            model_deploy_config=model_deploy_config,
            data_split_config=data_split_config,
        )
        return json.dumps(dict(AutoMLJobArn=auto_ml_job_arn))

    def describe_auto_ml_job_v2(self) -> str:
        auto_ml_job_name = self._get_param("AutoMLJobName")
        auto_ml_job_description = self.sagemaker_backend.describe_auto_ml_job_v2(
            auto_ml_job_name=auto_ml_job_name,
        )
        return json.dumps(auto_ml_job_description)

    def list_auto_ml_jobs(self) -> str:
        creation_time_after = self._get_param("CreationTimeAfter")
        creation_time_before = self._get_param("CreationTimeBefore")
        last_modified_time_after = self._get_param("LastModifiedTimeAfter")
        last_modified_time_before = self._get_param("LastModifiedTimeBefore")
        name_contains = self._get_param("NameContains")
        status_equals = self._get_param("StatusEquals")
        sort_order = self._get_param("SortOrder")
        sort_by = self._get_param("SortBy")
        max_results = self._get_param("MaxResults")
        next_token = self._get_param("NextToken")
        auto_ml_jobs, next_token = self.sagemaker_backend.list_auto_ml_jobs(
            creation_time_after=creation_time_after,
            creation_time_before=creation_time_before,
            last_modified_time_after=last_modified_time_after,
            last_modified_time_before=last_modified_time_before,
            name_contains=name_contains,
            status_equals=status_equals,
            sort_order=sort_order,
            sort_by=sort_by,
            max_results=max_results,
            next_token=next_token,
        )
        auto_ml_job_summaries = [auto_ml_job.summary() for auto_ml_job in auto_ml_jobs]
        return json.dumps(
            dict(AutoMLJobSummaries=auto_ml_job_summaries, NextToken=next_token)
        )

    def stop_auto_ml_job(self) -> str:
        auto_ml_job_name = self._get_param("AutoMLJobName")
        self.sagemaker_backend.stop_auto_ml_job(
            auto_ml_job_name=auto_ml_job_name,
        )
        return json.dumps(dict())

    def list_endpoints(self) -> str:
        sort_by = self._get_param("SortBy")
        sort_order = self._get_param("SortOrder")
        next_token = self._get_param("NextToken")
        max_results = self._get_param("MaxResults")
        name_contains = self._get_param("NameContains")
        creation_time_before = self._get_param("CreationTimeBefore")
        creation_time_after = self._get_param("CreationTimeAfter")
        last_modified_time_before = self._get_param("LastModifiedTimeBefore")
        last_modified_time_after = self._get_param("LastModifiedTimeAfter")
        status_equals = self._get_param("StatusEquals")
        endpoints, next_token = self.sagemaker_backend.list_endpoints(
            sort_by=sort_by,
            sort_order=sort_order,
            next_token=next_token,
            max_results=max_results,
            name_contains=name_contains,
            creation_time_before=creation_time_before,
            creation_time_after=creation_time_after,
            last_modified_time_before=last_modified_time_before,
            last_modified_time_after=last_modified_time_after,
            status_equals=status_equals,
        )
        endpoint_summaries = [endpoint.summary() for endpoint in endpoints]
        return json.dumps(dict(Endpoints=endpoint_summaries, NextToken=next_token))

    def list_endpoint_configs(self) -> str:
        sort_by = self._get_param("SortBy")
        sort_order = self._get_param("SortOrder")
        next_token = self._get_param("NextToken")
        max_results = self._get_param("MaxResults")
        name_contains = self._get_param("NameContains")
        creation_time_before = self._get_param("CreationTimeBefore")
        creation_time_after = self._get_param("CreationTimeAfter")
        endpoint_configs, next_token = self.sagemaker_backend.list_endpoint_configs(
            sort_by=sort_by,
            sort_order=sort_order,
            next_token=next_token,
            max_results=max_results,
            name_contains=name_contains,
            creation_time_before=creation_time_before,
            creation_time_after=creation_time_after,
        )
        endpoint_summaries = [
            endpoint_config.summary() for endpoint_config in endpoint_configs
        ]
        return json.dumps(
            dict(EndpointConfigs=endpoint_summaries, NextToken=next_token)
        )

    def create_compilation_job(self) -> str:
        compilation_job_name = self._get_param("CompilationJobName")
        role_arn = self._get_param("RoleArn")
        model_package_version_arn = self._get_param("ModelPackageVersionArn")
        input_config = self._get_param("InputConfig")
        output_config = self._get_param("OutputConfig")
        vpc_config = self._get_param("VpcConfig")
        stopping_condition = self._get_param("StoppingCondition")
        tags = self._get_param("Tags")
        compilation_job_arn = self.sagemaker_backend.create_compilation_job(
            compilation_job_name=compilation_job_name,
            role_arn=role_arn,
            model_package_version_arn=model_package_version_arn,
            input_config=input_config,
            output_config=output_config,
            vpc_config=vpc_config,
            stopping_condition=stopping_condition,
            tags=tags,
        )
        return json.dumps(dict(CompilationJobArn=compilation_job_arn))

    def describe_compilation_job(self) -> str:
        compilation_job_name = self._get_param("CompilationJobName")
        compilation_job_description = self.sagemaker_backend.describe_compilation_job(
            compilation_job_name=compilation_job_name,
        )
        return json.dumps(compilation_job_description)

    def list_compilation_jobs(self) -> str:
        next_token = self._get_param("NextToken")
        max_results = self._get_param("MaxResults")
        creation_time_after = self._get_param("CreationTimeAfter")
        creation_time_before = self._get_param("CreationTimeBefore")
        last_modified_time_after = self._get_param("LastModifiedTimeAfter")
        last_modified_time_before = self._get_param("LastModifiedTimeBefore")
        name_contains = self._get_param("NameContains")
        status_equals = self._get_param("StatusEquals")
        sort_by = self._get_param("SortBy")
        sort_order = self._get_param("SortOrder")
        compilation_jobs, next_token = self.sagemaker_backend.list_compilation_jobs(
            next_token=next_token,
            max_results=max_results,
            creation_time_after=creation_time_after,
            creation_time_before=creation_time_before,
            last_modified_time_after=last_modified_time_after,
            last_modified_time_before=last_modified_time_before,
            name_contains=name_contains,
            status_equals=status_equals,
            sort_by=sort_by,
            sort_order=sort_order,
        )
        compilation_job_summaries = [x.summary() for x in compilation_jobs]
        return json.dumps(
            dict(
                CompilationJobSummaries=compilation_job_summaries, NextToken=next_token
            )
        )

    def delete_compilation_job(self) -> str:
        compilation_job_name = self._get_param("CompilationJobName")
        self.sagemaker_backend.delete_compilation_job(
            compilation_job_name=compilation_job_name,
        )
        return json.dumps({})

    def create_domain(self) -> str:
        domain_name = self._get_param("DomainName")
        auth_mode = self._get_param("AuthMode")
        default_user_settings = self._get_param("DefaultUserSettings")
        domain_settings = self._get_param("DomainSettings")
        subnet_ids = self._get_param("SubnetIds")
        vpc_id = self._get_param("VpcId")
        tags = self._get_param("Tags")
        app_network_access_type = self._get_param("AppNetworkAccessType")
        home_efs_file_system_kms_key_id = self._get_param("HomeEfsFileSystemKmsKeyId")
        kms_key_id = self._get_param("KmsKeyId")
        app_security_group_management = self._get_param("AppSecurityGroupManagement")
        default_space_settings = self._get_param("DefaultSpaceSettings")
        resp = self.sagemaker_backend.create_domain(
            domain_name=domain_name,
            auth_mode=auth_mode,
            default_user_settings=default_user_settings,
            domain_settings=domain_settings,
            subnet_ids=subnet_ids,
            vpc_id=vpc_id,
            tags=tags,
            app_network_access_type=app_network_access_type,
            home_efs_file_system_kms_key_id=home_efs_file_system_kms_key_id,
            kms_key_id=kms_key_id,
            app_security_group_management=app_security_group_management,
            default_space_settings=default_space_settings,
        )
        return json.dumps(resp)

    def describe_domain(self) -> str:
        domain_id = self._get_param("DomainId")
        domain_description = self.sagemaker_backend.describe_domain(
            domain_id=domain_id,
        )
        return json.dumps(domain_description)

    def list_domains(self) -> str:
        next_token = self._get_param("NextToken")
        max_results = self._get_param("MaxResults")
        domains, next_token = self.sagemaker_backend.list_domains(
            next_token=next_token,
            max_results=max_results,
        )
        domain_summaries = [domain.summary() for domain in domains]
        return json.dumps(dict(Domains=domain_summaries, NextToken=next_token))

    def delete_domain(self) -> str:
        domain_id = self._get_param("DomainId")
        retention_policy = self._get_param("RetentionPolicy")
        self.sagemaker_backend.delete_domain(
            domain_id=domain_id,
            retention_policy=retention_policy,
        )
        return json.dumps(dict())

    def create_model_explainability_job_definition(self) -> str:
        job_definition_name = self._get_param("JobDefinitionName")
        model_explainability_baseline_config = self._get_param(
            "ModelExplainabilityBaselineConfig"
        )
        model_explainability_app_specification = self._get_param(
            "ModelExplainabilityAppSpecification"
        )
        model_explainability_job_input = self._get_param("ModelExplainabilityJobInput")
        model_explainability_job_output_config = self._get_param(
            "ModelExplainabilityJobOutputConfig"
        )
        job_resources = self._get_param("JobResources")
        network_config = self._get_param("NetworkConfig")
        role_arn = self._get_param("RoleArn")
        stopping_condition = self._get_param("StoppingCondition")
        tags = self._get_param("Tags")
        job_definition_arn = self.sagemaker_backend.create_model_explainability_job_definition(
            job_definition_name=job_definition_name,
            model_explainability_baseline_config=model_explainability_baseline_config,
            model_explainability_app_specification=model_explainability_app_specification,
            model_explainability_job_input=model_explainability_job_input,
            model_explainability_job_output_config=model_explainability_job_output_config,
            job_resources=job_resources,
            network_config=network_config,
            role_arn=role_arn,
            stopping_condition=stopping_condition,
            tags=tags,
        )
        return json.dumps(dict(JobDefinitionArn=job_definition_arn))

    def describe_model_explainability_job_definition(self) -> str:
        job_definition_name = self._get_param("JobDefinitionName")
        description = (
            self.sagemaker_backend.describe_model_explainability_job_definition(
                job_definition_name=job_definition_name,
            )
        )
        return json.dumps(description)

    def list_model_explainability_job_definitions(self) -> str:
        endpoint_name = self._get_param("EndpointName")
        sort_by = self._get_param("SortBy")
        sort_order = self._get_param("SortOrder")
        next_token = self._get_param("NextToken")
        max_results = self._get_param("MaxResults")
        name_contains = self._get_param("NameContains")
        creation_time_before = self._get_param("CreationTimeBefore")
        creation_time_after = self._get_param("CreationTimeAfter")
        job_definitions, next_token = (
            self.sagemaker_backend.list_model_explainability_job_definitions(
                endpoint_name=endpoint_name,
                sort_by=sort_by,
                sort_order=sort_order,
                next_token=next_token,
                max_results=max_results,
                name_contains=name_contains,
                creation_time_before=creation_time_before,
                creation_time_after=creation_time_after,
            )
        )
        job_definition_summaries = [job.summary() for job in job_definitions]
        return json.dumps(
            dict(JobDefinitionSummaries=job_definition_summaries, NextToken=next_token)
        )

    def delete_model_explainability_job_definition(self) -> str:
        job_definition_name = self._get_param("JobDefinitionName")
        self.sagemaker_backend.delete_model_explainability_job_definition(
            job_definition_name=job_definition_name,
        )
        return json.dumps(dict())

    def create_hyper_parameter_tuning_job(self) -> str:
        hyper_parameter_tuning_job_name = self._get_param("HyperParameterTuningJobName")
        hyper_parameter_tuning_job_config = self._get_param(
            "HyperParameterTuningJobConfig"
        )
        training_job_definition = self._get_param("TrainingJobDefinition")
        training_job_definitions = self._get_param("TrainingJobDefinitions")
        warm_start_config = self._get_param("WarmStartConfig")
        tags = self._get_param("Tags")
        autotune = self._get_param("Autotune")
        hyper_parameter_tuning_job_arn = (
            self.sagemaker_backend.create_hyper_parameter_tuning_job(
                hyper_parameter_tuning_job_name=hyper_parameter_tuning_job_name,
                hyper_parameter_tuning_job_config=hyper_parameter_tuning_job_config,
                training_job_definition=training_job_definition,
                training_job_definitions=training_job_definitions,
                warm_start_config=warm_start_config,
                tags=tags,
                autotune=autotune,
            )
        )
        return json.dumps(
            dict(HyperParameterTuningJobArn=hyper_parameter_tuning_job_arn)
        )

    def describe_hyper_parameter_tuning_job(self) -> str:
        hyper_parameter_tuning_job_name = self._get_param("HyperParameterTuningJobName")
        hyper_parameter_tuning_job_description = (
            self.sagemaker_backend.describe_hyper_parameter_tuning_job(
                hyper_parameter_tuning_job_name=hyper_parameter_tuning_job_name,
            )
        )
        return json.dumps(hyper_parameter_tuning_job_description)

    def list_hyper_parameter_tuning_jobs(self) -> str:
        next_token = self._get_param("NextToken")
        max_results = self._get_param("MaxResults")
        sort_by = self._get_param("SortBy")
        sort_order = self._get_param("SortOrder")
        name_contains = self._get_param("NameContains")
        creation_time_after = self._get_param("CreationTimeAfter")
        creation_time_before = self._get_param("CreationTimeBefore")
        last_modified_time_after = self._get_param("LastModifiedTimeAfter")
        last_modified_time_before = self._get_param("LastModifiedTimeBefore")
        status_equals = self._get_param("StatusEquals")
        hyper_parameter_tuning_jobs, next_token = (
            self.sagemaker_backend.list_hyper_parameter_tuning_jobs(
                next_token=next_token,
                max_results=max_results,
                sort_by=sort_by,
                sort_order=sort_order,
                name_contains=name_contains,
                creation_time_after=creation_time_after,
                creation_time_before=creation_time_before,
                last_modified_time_after=last_modified_time_after,
                last_modified_time_before=last_modified_time_before,
                status_equals=status_equals,
            )
        )
        hyper_parameter_tuning_job_summaries = [
            job.summary() for job in hyper_parameter_tuning_jobs
        ]
        return json.dumps(
            dict(
                HyperParameterTuningJobSummaries=hyper_parameter_tuning_job_summaries,
                NextToken=next_token,
            )
        )

    def delete_hyper_parameter_tuning_job(self) -> str:
        hyper_parameter_tuning_job_name = self._get_param("HyperParameterTuningJobName")
        self.sagemaker_backend.delete_hyper_parameter_tuning_job(
            hyper_parameter_tuning_job_name=hyper_parameter_tuning_job_name,
        )
        return json.dumps(dict())

    def create_model_quality_job_definition(self) -> str:
        job_definition_name = self._get_param("JobDefinitionName")
        model_quality_baseline_config = self._get_param("ModelQualityBaselineConfig")
        model_quality_app_specification = self._get_param(
            "ModelQualityAppSpecification"
        )
        model_quality_job_input = self._get_param("ModelQualityJobInput")
        model_quality_job_output_config = self._get_param("ModelQualityJobOutputConfig")
        job_resources = self._get_param("JobResources")
        network_config = self._get_param("NetworkConfig")
        role_arn = self._get_param("RoleArn")
        stopping_condition = self._get_param("StoppingCondition")
        tags = self._get_param("Tags")
        job_definition_arn = self.sagemaker_backend.create_model_quality_job_definition(
            job_definition_name=job_definition_name,
            model_quality_baseline_config=model_quality_baseline_config,
            model_quality_app_specification=model_quality_app_specification,
            model_quality_job_input=model_quality_job_input,
            model_quality_job_output_config=model_quality_job_output_config,
            job_resources=job_resources,
            network_config=network_config,
            role_arn=role_arn,
            stopping_condition=stopping_condition,
            tags=tags,
        )
        return json.dumps(dict(JobDefinitionArn=job_definition_arn))

    def describe_model_quality_job_definition(self) -> str:
        job_definition_name = self._get_param("JobDefinitionName")
        description = self.sagemaker_backend.describe_model_quality_job_definition(
            job_definition_name=job_definition_name,
        )
        return json.dumps(description)

    def list_model_quality_job_definitions(self) -> str:
        endpoint_name = self._get_param("EndpointName")
        sort_by = self._get_param("SortBy")
        sort_order = self._get_param("SortOrder")
        next_token = self._get_param("NextToken")
        max_results = self._get_param("MaxResults")
        name_contains = self._get_param("NameContains")
        creation_time_before = self._get_param("CreationTimeBefore")
        creation_time_after = self._get_param("CreationTimeAfter")
        job_definitions, next_token = (
            self.sagemaker_backend.list_model_quality_job_definitions(
                endpoint_name=endpoint_name,
                sort_by=sort_by,
                sort_order=sort_order,
                next_token=next_token,
                max_results=max_results,
                name_contains=name_contains,
                creation_time_before=creation_time_before,
                creation_time_after=creation_time_after,
            )
        )
        job_definition_summaries = [x.summary() for x in job_definitions]
        return json.dumps(
            dict(JobDefinitionSummaries=job_definition_summaries, NextToken=next_token)
        )

    def delete_model_quality_job_definition(self) -> str:
        job_definition_name = self._get_param("JobDefinitionName")
        self.sagemaker_backend.delete_model_quality_job_definition(
            job_definition_name=job_definition_name,
        )
        return json.dumps(dict())

    def create_model_card(self) -> str:
        model_card_name = self._get_param("ModelCardName")
        security_config = self._get_param("SecurityConfig")
        content = self._get_param("Content")
        model_card_status = self._get_param("ModelCardStatus")
        tags = self._get_param("Tags")
        model_card_arn = self.sagemaker_backend.create_model_card(
            model_card_name=model_card_name,
            security_config=security_config,
            content=content,
            model_card_status=model_card_status,
            tags=tags,
        )
        return json.dumps(dict(ModelCardArn=model_card_arn))

    def list_model_cards(self) -> str:
        creation_time_after = self._get_param("CreationTimeAfter")
        creation_time_before = self._get_param("CreationTimeBefore")
        max_results = self._get_param("MaxResults")
        name_contains = self._get_param("NameContains")
        model_card_status = self._get_param("ModelCardStatus")
        next_token = self._get_param("NextToken")
        sort_by = self._get_param("SortBy")
        sort_order = self._get_param("SortOrder")
        model_cards, next_token = self.sagemaker_backend.list_model_cards(
            creation_time_after=creation_time_after,
            creation_time_before=creation_time_before,
            max_results=max_results,
            name_contains=name_contains,
            model_card_status=model_card_status,
            next_token=next_token,
            sort_by=sort_by,
            sort_order=sort_order,
        )
        model_card_summaries = [model_card.summary() for model_card in model_cards]
        return json.dumps(
            dict(ModelCardSummaries=model_card_summaries, NextToken=next_token)
        )

    def list_model_card_versions(self) -> str:
        creation_time_after = self._get_param("CreationTimeAfter")
        creation_time_before = self._get_param("CreationTimeBefore")
        max_results = self._get_param("MaxResults")
        model_card_name = self._get_param("ModelCardName")
        model_card_status = self._get_param("ModelCardStatus")
        next_token = self._get_param("NextToken")
        sort_by = self._get_param("SortBy")
        sort_order = self._get_param("SortOrder")
        model_card_versions, next_token = (
            self.sagemaker_backend.list_model_card_versions(
                creation_time_after=creation_time_after,
                creation_time_before=creation_time_before,
                max_results=max_results,
                model_card_name=model_card_name,
                model_card_status=model_card_status,
                next_token=next_token,
                sort_by=sort_by,
                sort_order=sort_order,
            )
        )
        model_card_version_summaries = [
            mcv.version_summary() for mcv in model_card_versions
        ]
        return json.dumps(
            dict(
                ModelCardVersionSummaryList=model_card_version_summaries,
                NextToken=next_token,
            )
        )

    def update_model_card(self) -> str:
        model_card_name = self._get_param("ModelCardName")
        content = self._get_param("Content")
        model_card_status = self._get_param("ModelCardStatus")
        model_card_arn = self.sagemaker_backend.update_model_card(
            model_card_name=model_card_name,
            content=content,
            model_card_status=model_card_status,
        )
        return json.dumps(dict(ModelCardArn=model_card_arn))

    def describe_model_card(self) -> str:
        model_card_name = self._get_param("ModelCardName")
        model_card_version = self._get_param("ModelCardVersion")
        model_card_description = self.sagemaker_backend.describe_model_card(
            model_card_name=model_card_name,
            model_card_version=model_card_version,
        )
        return json.dumps(model_card_description)

    def delete_model_card(self) -> str:
        model_card_name = self._get_param("ModelCardName")
        self.sagemaker_backend.delete_model_card(
            model_card_name=model_card_name,
        )
        return json.dumps(dict())
