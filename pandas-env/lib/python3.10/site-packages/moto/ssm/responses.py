import json
from typing import Any, Dict, Tuple, Union

from moto.core.responses import BaseResponse

from .exceptions import ValidationException
from .models import SimpleSystemManagerBackend, ssm_backends


class SimpleSystemManagerResponse(BaseResponse):
    def __init__(self) -> None:
        super().__init__(service_name="ssm")

    @property
    def ssm_backend(self) -> SimpleSystemManagerBackend:
        return ssm_backends[self.current_account][self.region]

    def create_document(self) -> str:
        content = self._get_param("Content")
        requires = self._get_param("Requires")
        attachments = self._get_param("Attachments")
        name = self._get_param("Name")
        version_name = self._get_param("VersionName")
        document_type = self._get_param("DocumentType")
        document_format = self._get_param("DocumentFormat", "JSON")
        target_type = self._get_param("TargetType")
        tags = self._get_param("Tags")

        result = self.ssm_backend.create_document(
            content=content,
            requires=requires,
            attachments=attachments,
            name=name,
            version_name=version_name,
            document_type=document_type,
            document_format=document_format,
            target_type=target_type,
            tags=tags,
        )

        return json.dumps({"DocumentDescription": result})

    def delete_document(self) -> str:
        name = self._get_param("Name")
        document_version = self._get_param("DocumentVersion")
        version_name = self._get_param("VersionName")
        force = self._get_param("Force", False)
        self.ssm_backend.delete_document(
            name=name,
            document_version=document_version,
            version_name=version_name,
            force=force,
        )

        return json.dumps({})

    def get_document(self) -> str:
        name = self._get_param("Name")
        version_name = self._get_param("VersionName")
        document_version = self._get_param("DocumentVersion")
        document_format = self._get_param("DocumentFormat", "JSON")

        document = self.ssm_backend.get_document(
            name=name,
            document_version=document_version,
            document_format=document_format,
            version_name=version_name,
        )

        return json.dumps(document)

    def describe_document(self) -> str:
        name = self._get_param("Name")
        document_version = self._get_param("DocumentVersion")
        version_name = self._get_param("VersionName")

        result = self.ssm_backend.describe_document(
            name=name, document_version=document_version, version_name=version_name
        )

        return json.dumps({"Document": result})

    def update_document(self) -> str:
        content = self._get_param("Content")
        attachments = self._get_param("Attachments")
        name = self._get_param("Name")
        version_name = self._get_param("VersionName")
        document_version = self._get_param("DocumentVersion")
        document_format = self._get_param("DocumentFormat", "JSON")
        target_type = self._get_param("TargetType")

        result = self.ssm_backend.update_document(
            content=content,
            attachments=attachments,
            name=name,
            version_name=version_name,
            document_version=document_version,
            document_format=document_format,
            target_type=target_type,
        )

        return json.dumps({"DocumentDescription": result})

    def update_document_default_version(self) -> str:
        name = self._get_param("Name")
        document_version = self._get_param("DocumentVersion")

        result = self.ssm_backend.update_document_default_version(
            name=name, document_version=document_version
        )
        return json.dumps({"Description": result})

    def list_documents(self) -> str:
        document_filter_list = self._get_param("DocumentFilterList")
        filters = self._get_param("Filters")
        max_results = self._get_param("MaxResults", 10)
        next_token = self._get_param("NextToken", "0")

        documents, token = self.ssm_backend.list_documents(
            document_filter_list=document_filter_list,
            filters=filters,
            max_results=max_results,
            token=next_token,
        )

        return json.dumps({"DocumentIdentifiers": documents, "NextToken": token})

    def describe_document_permission(self) -> str:
        name = self._get_param("Name")

        result = self.ssm_backend.describe_document_permission(name=name)
        return json.dumps(result)

    def modify_document_permission(self) -> str:
        account_ids_to_add = self._get_param("AccountIdsToAdd")
        account_ids_to_remove = self._get_param("AccountIdsToRemove")
        name = self._get_param("Name")
        permission_type = self._get_param("PermissionType")
        shared_document_version = self._get_param("SharedDocumentVersion")

        self.ssm_backend.modify_document_permission(
            name=name,
            account_ids_to_add=account_ids_to_add,
            account_ids_to_remove=account_ids_to_remove,
            shared_document_version=shared_document_version,
            permission_type=permission_type,
        )
        return "{}"

    def delete_parameter(self) -> Union[str, Tuple[str, Dict[str, int]]]:
        name = self._get_param("Name")
        result = self.ssm_backend.delete_parameter(name)
        if result is None:
            error = {
                "__type": "ParameterNotFound",
                "message": f"Parameter {name} not found.",
            }
            return json.dumps(error), dict(status=400)
        return json.dumps({})

    def delete_parameters(self) -> str:
        names = self._get_param("Names")
        result = self.ssm_backend.delete_parameters(names)

        response: Dict[str, Any] = {"DeletedParameters": [], "InvalidParameters": []}

        for name in names:
            if name in result:
                response["DeletedParameters"].append(name)
            else:
                response["InvalidParameters"].append(name)
        return json.dumps(response)

    def get_parameter(self) -> Union[str, Tuple[str, Dict[str, int]]]:
        name = self._get_param("Name")
        with_decryption = self._get_param("WithDecryption")

        if (
            name.startswith("/aws/reference/secretsmanager/")
            and with_decryption is not True
        ):
            raise ValidationException(
                "WithDecryption flag must be True for retrieving a Secret Manager secret."
            )

        result = self.ssm_backend.get_parameter(name)

        if result is None:
            error = {
                "__type": "ParameterNotFound",
                "message": f"Parameter {name} not found.",
            }
            return json.dumps(error), dict(status=400)

        response = {"Parameter": result.response_object(with_decryption, self.region)}
        return json.dumps(response)

    def get_parameters(self) -> str:
        names = self._get_param("Names")
        with_decryption = self._get_param("WithDecryption")

        result = self.ssm_backend.get_parameters(names)

        response: Dict[str, Any] = {"Parameters": [], "InvalidParameters": []}

        for name, parameter in result.items():
            param_data = parameter.response_object(with_decryption, self.region)
            response["Parameters"].append(param_data)

        valid_param_names = [name for name, parameter in result.items()]
        for name in names:
            if name not in valid_param_names:
                response["InvalidParameters"].append(name)
        return json.dumps(response)

    def get_parameters_by_path(self) -> str:
        path = self._get_param("Path")
        with_decryption = self._get_param("WithDecryption")
        recursive = self._get_param("Recursive", False)
        filters = self._get_param("ParameterFilters")
        token = self._get_param("NextToken")
        max_results = self._get_param("MaxResults", 10)

        result, next_token = self.ssm_backend.get_parameters_by_path(
            path,
            recursive,
            filters,
            next_token=token,
            max_results=max_results,
        )

        response: Dict[str, Any] = {"Parameters": [], "NextToken": next_token}

        for parameter in result:
            param_data = parameter.response_object(with_decryption, self.region)
            response["Parameters"].append(param_data)

        return json.dumps(response)

    def describe_parameters(self) -> str:
        page_size = 10
        filters = self._get_param("Filters")
        parameter_filters = self._get_param("ParameterFilters")
        token = self._get_param("NextToken")
        if hasattr(token, "strip"):
            token = token.strip()
        if not token:
            token = "0"
        token = int(token)

        result = self.ssm_backend.describe_parameters(filters, parameter_filters)

        response: Dict[str, Any] = {"Parameters": []}

        end = token + page_size
        for parameter in result[token:]:
            response["Parameters"].append(parameter.describe_response_object(False))

            token += 1
            if len(response["Parameters"]) == page_size:
                response["NextToken"] = str(end)
                break

        return json.dumps(response)

    def put_parameter(self) -> Union[str, Tuple[str, Dict[str, int]]]:
        name = self._get_param("Name")
        description = self._get_param("Description")
        value = self._get_param("Value")
        type_ = self._get_param("Type")
        allowed_pattern = self._get_param("AllowedPattern")
        keyid = self._get_param("KeyId")
        overwrite = self._get_param("Overwrite", False)
        tags = self._get_param("Tags")
        data_type = self._get_param("DataType", "text")
        tier = self._get_param("Tier")
        policies = self._get_param("Policies")

        param = self.ssm_backend.put_parameter(
            name,
            description,
            value,
            type_,
            allowed_pattern,
            keyid,
            overwrite,
            tags,
            data_type,
            tier=tier,
            policies=policies,
        )

        response = {"Version": param.version}
        return json.dumps(response)

    def get_parameter_history(self) -> Union[str, Tuple[str, Dict[str, int]]]:
        name = self._get_param("Name")
        with_decryption = self._get_param("WithDecryption")
        next_token = self._get_param("NextToken")
        max_results = self._get_param("MaxResults", 50)

        result, new_next_token = self.ssm_backend.get_parameter_history(
            name, next_token, max_results
        )

        if result is None:
            error = {
                "__type": "ParameterNotFound",
                "message": f"Parameter {name} not found.",
            }
            return json.dumps(error), dict(status=400)

        response = {
            "Parameters": [
                p_v.describe_response_object(
                    decrypt=with_decryption, include_labels=True
                )
                for p_v in result
            ],
            "NextToken": new_next_token,
        }

        return json.dumps(response)

    def label_parameter_version(self) -> str:
        name = self._get_param("Name")
        version = self._get_param("ParameterVersion")
        labels = self._get_param("Labels")

        invalid_labels, version = self.ssm_backend.label_parameter_version(
            name, version, labels
        )

        response = {"InvalidLabels": invalid_labels, "ParameterVersion": version}
        return json.dumps(response)

    def add_tags_to_resource(self) -> str:
        resource_id = self._get_param("ResourceId")
        resource_type = self._get_param("ResourceType")
        tags = {t["Key"]: t["Value"] for t in self._get_param("Tags")}
        self.ssm_backend.add_tags_to_resource(
            resource_type=resource_type, resource_id=resource_id, tags=tags
        )
        return json.dumps({})

    def remove_tags_from_resource(self) -> str:
        resource_id = self._get_param("ResourceId")
        resource_type = self._get_param("ResourceType")
        keys = self._get_param("TagKeys")
        self.ssm_backend.remove_tags_from_resource(
            resource_type=resource_type, resource_id=resource_id, keys=keys
        )
        return json.dumps({})

    def list_tags_for_resource(self) -> str:
        resource_id = self._get_param("ResourceId")
        resource_type = self._get_param("ResourceType")
        tags = self.ssm_backend.list_tags_for_resource(
            resource_type=resource_type, resource_id=resource_id
        )
        tag_list = [{"Key": k, "Value": v} for (k, v) in tags.items()]
        response = {"TagList": tag_list}
        return json.dumps(response)

    def send_command(self) -> str:
        comment = self._get_param("Comment", "")
        document_name = self._get_param("DocumentName")
        timeout_seconds = self._get_int_param("TimeoutSeconds")
        instance_ids = self._get_param("InstanceIds", [])
        max_concurrency = self._get_param("MaxConcurrency", "50")
        max_errors = self._get_param("MaxErrors", "0")
        notification_config = self._get_param("NotificationConfig")
        output_s3_bucket_name = self._get_param("OutputS3BucketName", "")
        output_s3_key_prefix = self._get_param("OutputS3KeyPrefix", "")
        output_s3_region = self._get_param("OutputS3Region", "")
        parameters = self._get_param("Parameters", {})
        service_role_arn = self._get_param("ServiceRoleArn", "")
        targets = self._get_param("Targets", [])
        command = self.ssm_backend.send_command(
            comment=comment,
            document_name=document_name,
            timeout_seconds=timeout_seconds,
            instance_ids=instance_ids,
            max_concurrency=max_concurrency,
            max_errors=max_errors,
            notification_config=notification_config,
            output_s3_bucket_name=output_s3_bucket_name,
            output_s3_key_prefix=output_s3_key_prefix,
            output_s3_region=output_s3_region,
            parameters=parameters,
            service_role_arn=service_role_arn,
            targets=targets,
        )
        return json.dumps({"Command": command.response_object()})

    def list_commands(self) -> str:
        command_id = self._get_param("CommandId")
        instance_id = self._get_param("InstanceId")
        commands = self.ssm_backend.list_commands(command_id, instance_id)
        response = {"Commands": [command.response_object() for command in commands]}
        return json.dumps(response)

    def get_command_invocation(self) -> str:
        command_id = self._get_param("CommandId")
        instance_id = self._get_param("InstanceId")
        plugin_name = self._get_param("PluginName")
        response = self.ssm_backend.get_command_invocation(
            command_id, instance_id, plugin_name
        )
        return json.dumps(response)

    def create_maintenance_window(self) -> str:
        name = self._get_param("Name")
        desc = self._get_param("Description", None)
        duration = self._get_int_param("Duration")
        cutoff = self._get_int_param("Cutoff")
        schedule = self._get_param("Schedule")
        schedule_timezone = self._get_param("ScheduleTimezone")
        schedule_offset = self._get_int_param("ScheduleOffset")
        start_date = self._get_param("StartDate")
        end_date = self._get_param("EndDate")
        tags = self._get_param("Tags")
        window_id = self.ssm_backend.create_maintenance_window(
            name=name,
            description=desc,
            duration=duration,
            cutoff=cutoff,
            schedule=schedule,
            schedule_timezone=schedule_timezone,
            schedule_offset=schedule_offset,
            start_date=start_date,
            end_date=end_date,
            tags=tags,
        )
        return json.dumps({"WindowId": window_id})

    def get_maintenance_window(self) -> str:
        window_id = self._get_param("WindowId")
        window = self.ssm_backend.get_maintenance_window(window_id)
        return json.dumps(window.to_json())

    def register_target_with_maintenance_window(self) -> str:
        window_target_id = self.ssm_backend.register_target_with_maintenance_window(
            window_id=self._get_param("WindowId"),
            resource_type=self._get_param("ResourceType"),
            targets=self._get_param("Targets"),
            owner_information=self._get_param("OwnerInformation"),
            name=self._get_param("Name"),
            description=self._get_param("Description"),
        )
        return json.dumps({"WindowTargetId": window_target_id})

    def describe_maintenance_window_targets(self) -> str:
        window_id = self._get_param("WindowId")
        filters = self._get_param("Filters", [])
        targets = [
            target.to_json()
            for target in self.ssm_backend.describe_maintenance_window_targets(
                window_id, filters
            )
        ]
        return json.dumps({"Targets": targets})

    def deregister_target_from_maintenance_window(self) -> str:
        window_id = self._get_param("WindowId")
        window_target_id = self._get_param("WindowTargetId")
        self.ssm_backend.deregister_target_from_maintenance_window(
            window_id, window_target_id
        )
        return "{}"

    def describe_maintenance_windows(self) -> str:
        filters = self._get_param("Filters", None)
        windows = [
            window.to_json()
            for window in self.ssm_backend.describe_maintenance_windows(filters)
        ]
        return json.dumps({"WindowIdentities": windows})

    def delete_maintenance_window(self) -> str:
        window_id = self._get_param("WindowId")
        self.ssm_backend.delete_maintenance_window(window_id)
        return "{}"

    def create_patch_baseline(self) -> str:
        baseline_id = self.ssm_backend.create_patch_baseline(
            name=self._get_param("Name"),
            operating_system=self._get_param("OperatingSystem"),
            global_filters=self._get_param("GlobalFilters", {}),
            approval_rules=self._get_param("ApprovalRules", {}),
            approved_patches=self._get_param("ApprovedPatches", []),
            approved_patches_compliance_level=self._get_param(
                "ApprovedPatchesComplianceLevel"
            ),
            approved_patches_enable_non_security=self._get_param(
                "ApprovedPatchesEnableNonSecurity"
            ),
            rejected_patches=self._get_param("RejectedPatches", []),
            rejected_patches_action=self._get_param("RejectedPatchesAction"),
            description=self._get_param("Description"),
            sources=self._get_param("Sources", []),
            tags=self._get_param("Tags", []),
        )
        return json.dumps({"BaselineId": baseline_id})

    def describe_patch_baselines(self) -> str:
        filters = self._get_param("Filters", None)
        baselines = [
            baseline.to_json()
            for baseline in self.ssm_backend.describe_patch_baselines(filters)
        ]
        return json.dumps({"BaselineIdentities": baselines})

    def delete_patch_baseline(self) -> str:
        baseline_id = self._get_param("BaselineId")
        self.ssm_backend.delete_patch_baseline(baseline_id)
        return "{}"

    def register_task_with_maintenance_window(self) -> str:
        window_task_id = self.ssm_backend.register_task_with_maintenance_window(
            window_id=self._get_param("WindowId"),
            targets=self._get_param("Targets"),
            task_arn=self._get_param("TaskArn"),
            service_role_arn=self._get_param("ServiceRoleArn"),
            task_type=self._get_param("TaskType"),
            task_parameters=self._get_param("TaskParameters"),
            task_invocation_parameters=self._get_param("TaskInvocationParameters"),
            priority=self._get_param("Priority"),
            max_concurrency=self._get_param("MaxConcurrency"),
            max_errors=self._get_param("MaxErrors"),
            logging_info=self._get_param("LoggingInfo"),
            name=self._get_param("Name"),
            description=self._get_param("Description"),
            cutoff_behavior=self._get_param("CutoffBehavior"),
            alarm_configurations=self._get_param("AlarmConfigurations"),
        )
        return json.dumps({"WindowTaskId": window_task_id})

    def describe_maintenance_window_tasks(self) -> str:
        window_id = self._get_param("WindowId")
        filters = self._get_param("Filters", [])
        tasks = [
            task.to_json()
            for task in self.ssm_backend.describe_maintenance_window_tasks(
                window_id, filters
            )
        ]
        return json.dumps({"Tasks": tasks})

    def deregister_task_from_maintenance_window(self) -> str:
        window_id = self._get_param("WindowId")
        window_task_id = self._get_param("WindowTaskId")
        self.ssm_backend.deregister_task_from_maintenance_window(
            window_id, window_task_id
        )
        return "{}"
