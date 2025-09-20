import json
import re
from typing import Any, Dict, List, Optional, Tuple, Type

import yaml
from yaml.parser import ParserError
from yaml.scanner import ScannerError

from moto.core.common_models import CloudFormationModel
from moto.core.responses import ActionResult, BaseResponse, EmptyResult
from moto.s3.exceptions import S3ClientError

from .exceptions import MissingParameterError, ValidationError
from .models import (
    Change,
    CloudFormationBackend,
    Stack,
    StackSet,
    cloudformation_backends,
)
from .utils import get_stack_from_s3_url, yaml_tag_constructor


class StackResourceDTO:
    def __init__(self, stack: Stack, resource: Type[CloudFormationModel]) -> None:
        self.stack_id = stack.stack_id
        self.stack_name = stack.name
        # FIXME: None of these attributes are part of CloudFormationModel interface...
        self.logical_resource_id = getattr(resource, "logical_resource_id", None)
        self.physical_resource_id = getattr(resource, "physical_resource_id", None)
        self.resource_type = getattr(resource, "cf_resource_type", None)
        self.stack_status = stack.status
        self.creation_time = stack.creation_time
        self.timestamp = "2010-07-27T22:27:28Z"  # Hardcoded in original XML template.
        self.resource_status = stack.status


class StackSetOperationDTO:
    def __init__(self, operation: Dict[str, Any], stack_set: StackSet) -> None:
        self.execution_role_name = stack_set.execution_role_name
        self.administration_role_arn = stack_set.administration_role_arn
        self.stack_set_id = stack_set.id
        self.creation_timestamp = operation["CreationTimestamp"]
        self.operation_id = operation["OperationId"]
        self.action = operation["Action"]
        self.end_timestamp = operation.get("EndTimestamp", None)
        self.status = operation["Status"]


class ChangeDTO:
    def __init__(self, change: Change) -> None:
        self.type = "Resource"
        self.resource_change = change


def get_template_summary_response_from_template(template_body: str) -> Dict[str, Any]:
    def get_resource_types(template_dict: Dict[str, Any]) -> List[Any]:
        resources = {}
        for key, value in template_dict.items():
            if key == "Resources":
                resources = value

        resource_types = []
        for key, value in resources.items():
            resource_types.append(value["Type"])
        return resource_types

    yaml.add_multi_constructor("", yaml_tag_constructor)

    try:
        template_dict = yaml.load(template_body, Loader=yaml.Loader)
    except (ParserError, ScannerError):
        template_dict = json.loads(template_body)

    resources_types = get_resource_types(template_dict)
    template_dict["ResourceTypes"] = resources_types
    template_dict["Version"] = template_dict["AWSTemplateFormatVersion"]
    parameters = []
    for key, value in template_dict.get("Parameters", {}).items():
        parameter = {
            "ParameterKey": key,
            "Description": value.get("Description", ""),
            "DefaultValue": value.get("Default", None),
            "NoEcho": value.get("NoEcho", False),
            "ParameterType": value.get("Type", "String"),
            "ParameterConstraints": {},
        }
        if "AllowedValues" in value:
            parameter["ParameterConstraints"]["AllowedValues"] = value["AllowedValues"]
        parameters.append(parameter)
    template_dict["Parameters"] = parameters
    return template_dict


def transform_dict(
    data: dict[str, str], key_for_key: str = "Key", key_for_value: str = "Value"
) -> list[dict[str, str]]:
    transformed = [
        {key_for_key: key, key_for_value: value} for key, value in data.items()
    ]
    return transformed


def transform_parameters(data: dict[str, str]) -> list[dict[str, str]]:
    return transform_dict(
        data, key_for_key="ParameterKey", key_for_value="ParameterValue"
    )


class CloudFormationResponse(BaseResponse):
    RESPONSE_KEY_PATH_TO_TRANSFORMER = {
        "DescribeChangeSetOutput.Changes": lambda x: [ChangeDTO(c) for c in x],
        "DescribeChangeSetOutput.Parameters": transform_parameters,
        "DescribeStackSetOutput.StackSet.Parameters": transform_parameters,
        "DescribeStackSetOutput.StackSet.Tags": transform_dict,
        "DescribeStacksOutput.Stacks.Stack.Tags": transform_dict,
    }

    def __init__(self) -> None:
        super().__init__(service_name="cloudformation")
        self.automated_parameter_parsing = True

    @property
    def cloudformation_backend(self) -> CloudFormationBackend:
        return cloudformation_backends[self.current_account][self.region]

    @classmethod
    def cfnresponse(cls, *args: Any, **kwargs: Any) -> Any:  # type: ignore[misc]
        request, full_url, headers = args
        full_url += "&Action=ProcessCfnResponse"
        cf = CloudFormationResponse()
        cf.automated_parameter_parsing = False
        return cf._dispatch(request=request, full_url=full_url, headers=headers)

    def _get_stack_from_s3_url(self, template_url: str) -> str:
        return get_stack_from_s3_url(
            template_url, account_id=self.current_account, partition=self.partition
        )

    def _get_params_from_list(
        self, parameters_list: List[Dict[str, Any]]
    ) -> Dict[str, Any]:
        # Hack dict-comprehension
        return dict(
            [
                (parameter["parameter_key"], parameter["parameter_value"])
                for parameter in parameters_list
            ]
        )

    def _get_param_values(
        self, parameters_list: List[Dict[str, str]], existing_params: Dict[str, str]
    ) -> Dict[str, Any]:
        result = {}
        for parameter in parameters_list:
            if set(parameter.keys()) >= {"parameter_key", "parameter_value"}:
                result[parameter["parameter_key"]] = parameter["parameter_value"]
            elif (
                set(parameter.keys()) >= {"parameter_key", "use_previous_value"}
                and parameter["parameter_key"] in existing_params
            ):
                result[parameter["parameter_key"]] = existing_params[
                    parameter["parameter_key"]
                ]
            else:
                raise MissingParameterError(parameter["parameter_key"])
        return result

    def process_cfn_response(self) -> Tuple[int, Dict[str, int], str]:
        status = self._get_param("Status")
        if status == "SUCCESS":
            stack_id = self._get_param("StackId")
            logical_resource_id = self._get_param("LogicalResourceId")
            outputs = self._get_param("Data")
            stack = self.cloudformation_backend.get_stack(stack_id)
            custom_resource = stack.get_custom_resource(logical_resource_id)
            custom_resource.set_data(outputs)
            stack.verify_readiness()

        return 200, {"status": 200}, json.dumps("{}")

    def create_stack(self) -> ActionResult:
        stack_name = self._get_param("StackName")
        stack_body = self._get_param("TemplateBody")
        template_url = self._get_param("TemplateURL")
        role_arn = self._get_param("RoleARN")
        enable_termination_protection = self._get_param("EnableTerminationProtection")
        timeout_in_mins = self._get_param("TimeoutInMinutes")
        stack_policy_body = self._get_param("StackPolicyBody")
        parameters_list = self._get_param("Parameters", [])
        tags = dict(
            (item["key"], item["value"]) for item in self._get_param("Tags", [])
        )

        parameters = self._get_params_from_list(parameters_list)

        if template_url:
            stack_body = self._get_stack_from_s3_url(template_url)
        stack_notification_arns = self._get_param("NotificationARNs", [])

        stack = self.cloudformation_backend.create_stack(
            name=stack_name,
            template=stack_body,
            parameters=parameters,
            notification_arns=stack_notification_arns,
            tags=tags,
            role_arn=role_arn,
            enable_termination_protection=enable_termination_protection,
            timeout_in_mins=timeout_in_mins,
            stack_policy_body=stack_policy_body,
        )
        result = {"StackId": stack.stack_id}
        return ActionResult(result)

    def validate_template_and_stack_body(self) -> None:
        if (
            self._get_param("TemplateBody") or self._get_param("TemplateURL")
        ) and self._get_bool_param("UsePreviousTemplate", False):
            raise ValidationError(
                message="An error occurred (ValidationError) when calling the CreateChangeSet operation: You cannot specify both usePreviousTemplate and Template Body/Template URL."
            )
        elif (
            not self._get_param("TemplateBody")
            and not self._get_param("TemplateURL")
            and not self._get_bool_param("UsePreviousTemplate", False)
        ):
            raise ValidationError(
                message="An error occurred (ValidationError) when calling the CreateChangeSet operation: Either Template URL or Template Body must be specified."
            )

    def create_change_set(self) -> ActionResult:
        stack_name = self._get_param("StackName")
        change_set_name = self._get_param("ChangeSetName")
        stack_body = self._get_param("TemplateBody")
        template_url = self._get_param("TemplateURL")
        update_or_create = self._get_param("ChangeSetType", "CREATE")
        use_previous_template = self._get_bool_param("UsePreviousTemplate", False)
        if update_or_create == "UPDATE":
            stack = self.cloudformation_backend.get_stack(stack_name)
            self.validate_template_and_stack_body()
            if use_previous_template:
                stack_body = stack.template
        description = self._get_param("Description")
        role_arn = self._get_param("RoleARN")
        parameters_list = self._get_param("Parameters", [])
        tags = dict(
            (item["key"], item["value"]) for item in self._get_param("Tags", [])
        )
        parameters = {
            param["parameter_key"]: (
                stack.stack_parameters[param["parameter_key"]]
                if param.get("use_previous_value", False)
                else param["parameter_value"]
            )
            for param in parameters_list
        }
        if update_or_create == "UPDATE":
            self._validate_different_update(parameters_list, stack_body, stack)

        if template_url:
            stack_body = self._get_stack_from_s3_url(template_url)
        stack_notification_arns = self._get_param("NotificationARNs", [])
        change_set_id, stack_id = self.cloudformation_backend.create_change_set(
            stack_name=stack_name,
            change_set_name=change_set_name,
            template=stack_body,
            parameters=parameters,
            description=description,
            notification_arns=stack_notification_arns,
            tags=tags,
            role_arn=role_arn,
            change_set_type=update_or_create,
        )
        result = {"Id": change_set_id, "StackId": stack_id}
        return ActionResult(result)

    def delete_change_set(self) -> ActionResult:
        change_set_name = self._get_param("ChangeSetName")

        self.cloudformation_backend.delete_change_set(change_set_name=change_set_name)
        return EmptyResult()

    def describe_change_set(self) -> ActionResult:
        change_set_name = self._get_param("ChangeSetName")
        change_set = self.cloudformation_backend.describe_change_set(
            change_set_name=change_set_name
        )
        return ActionResult(change_set)

    def execute_change_set(self) -> ActionResult:
        stack_name = self._get_param("StackName")
        change_set_name = self._get_param("ChangeSetName")
        self.cloudformation_backend.execute_change_set(
            stack_name=stack_name, change_set_name=change_set_name
        )
        return EmptyResult()

    def describe_stacks(self) -> ActionResult:
        stack_name_or_id = self._get_param("StackName")
        token = self._get_param("NextToken")
        stacks = self.cloudformation_backend.describe_stacks(stack_name_or_id)
        stack_ids = [stack.stack_id for stack in stacks]
        if token:
            start = stack_ids.index(token) + 1
        else:
            start = 0
        max_results = 50  # using this to mske testing of paginated stacks more convenient than default 1 MB
        stacks_resp = stacks[start : start + max_results]
        next_token = None
        if len(stacks) > (start + max_results):
            next_token = stacks_resp[-1].stack_id
        result = {"Stacks": stacks_resp, "NextToken": next_token}
        return ActionResult(result)

    def describe_stack_resource(self) -> ActionResult:
        stack_name = self._get_param("StackName")
        logical_resource_id = self._get_param("LogicalResourceId")
        stack, resource = self.cloudformation_backend.describe_stack_resource(
            stack_name, logical_resource_id
        )
        result = {"StackResourceDetail": StackResourceDTO(stack, resource)}
        return ActionResult(result)

    def describe_stack_resources(self) -> ActionResult:
        stack_name = self._get_param("StackName")
        logical_resource_id = self._get_param("LogicalResourceId")
        stack, resources = self.cloudformation_backend.describe_stack_resources(
            stack_name, logical_resource_id
        )

        result = {
            "StackResources": [
                StackResourceDTO(stack, resource) for resource in resources
            ]
        }
        return ActionResult(result)

    def describe_stack_events(self) -> ActionResult:
        stack_name = self._get_param("StackName")
        events = self.cloudformation_backend.describe_stack_events(stack_name)

        result = {"StackEvents": events[::-1]}
        return ActionResult(result)

    def list_change_sets(self) -> ActionResult:
        change_sets = self.cloudformation_backend.list_change_sets()
        result = {"Summaries": change_sets}
        return ActionResult(result)

    def list_stacks(self) -> ActionResult:
        status_filter = self._get_param("StackStatusFilter", [])
        stacks = self.cloudformation_backend.list_stacks(status_filter)
        result = {"StackSummaries": stacks}
        return ActionResult(result)

    def list_stack_resources(self) -> ActionResult:
        stack_name_or_id = self._get_param("StackName")
        resources = self.cloudformation_backend.list_stack_resources(stack_name_or_id)
        # FIXME: This was migrated from the original XML template, including hardcoded values.
        result = {
            "StackResourceSummaries": [
                {
                    "ResourceStatus": "CREATE_COMPLETE",
                    "LogicalResourceId": getattr(resource, "logical_resource_id", ""),
                    "LastUpdateTimestamp": "2011-06-21T20:15:58Z",
                    "PhysicalResourceId": getattr(resource, "physical_resource_id", ""),
                    "ResourceType": getattr(resource, "cf_resource_type", ""),
                }
                for resource in resources
            ]
        }
        return ActionResult(result)

    def get_template(self) -> ActionResult:
        name_or_stack_id = self._get_param("StackName")
        stack_template = self.cloudformation_backend.get_template(name_or_stack_id)
        result = {"TemplateBody": stack_template}
        return ActionResult(result)

    def get_template_summary(self) -> ActionResult:
        stack_name = self._get_param("StackName")
        template_url = self._get_param("TemplateURL")
        stack_body = self._get_param("TemplateBody")

        if stack_name:
            stack = self.cloudformation_backend.get_stack(stack_name)
            if stack.status == "REVIEW_IN_PROGRESS":
                raise ValidationError(
                    message="GetTemplateSummary cannot be called on REVIEW_IN_PROGRESS stacks."
                )
            stack_body = stack.template
        elif template_url:
            stack_body = self._get_stack_from_s3_url(template_url)

        template_summary = get_template_summary_response_from_template(stack_body)
        return ActionResult(template_summary)

    def _validate_different_update(
        self,
        incoming_params: Optional[List[Dict[str, Any]]],
        stack_body: str,
        old_stack: Stack,
    ) -> None:
        if incoming_params and stack_body:
            new_params = self._get_param_values(
                incoming_params, old_stack.stack_parameters
            )
            if (
                old_stack.template == stack_body
                and old_stack.stack_parameters == new_params
            ):
                raise ValidationError(
                    old_stack.name, message="No updates are to be performed."
                )

    def _validate_status(self, stack: Stack) -> None:
        if stack.status == "ROLLBACK_COMPLETE":
            raise ValidationError(
                stack.stack_id,
                message=f"Stack:{stack.stack_id} is in ROLLBACK_COMPLETE state and can not be updated.",
            )

    def update_stack(self) -> ActionResult:
        stack_name = self._get_param("StackName")
        role_arn = self._get_param("RoleARN")
        template_url = self._get_param("TemplateURL")
        stack_body = self._get_param("TemplateBody")
        stack = self.cloudformation_backend.get_stack(stack_name)
        if self._get_bool_param("UsePreviousTemplate", False):
            stack_body = stack.template
        elif not stack_body and template_url:
            stack_body = self._get_stack_from_s3_url(template_url)

        incoming_params = self._get_param("Parameters", [])
        for param in incoming_params:
            if param.get("use_previous_value") and param.get("parameter_value"):
                raise ValidationError(
                    message=f"Invalid input for parameter key {param['parameter_key']}. Cannot specify usePreviousValue as true and non empty value for a parameter"
                )
        # boto3 is supposed to let you clear the tags by passing an empty value, but the request body doesn't
        # end up containing anything we can use to differentiate between passing an empty value versus not
        # passing anything. so until that changes, moto won't be able to clear tags, only update them.
        tags: Optional[Dict[str, str]] = dict(
            (item["key"], item["value"]) for item in self._get_param("Tags", [])
        )
        # so that if we don't pass the parameter, we don't clear all the tags accidentally
        if not tags:
            tags = None

        stack = self.cloudformation_backend.get_stack(stack_name)
        self._validate_different_update(incoming_params, stack_body, stack)
        self._validate_status(stack)

        stack = self.cloudformation_backend.update_stack(
            name=stack_name,
            template=stack_body,
            role_arn=role_arn,
            parameters=incoming_params,
            tags=tags,
        )
        result = {"StackId": stack.stack_id}
        return ActionResult(result)

    def delete_stack(self) -> ActionResult:
        name_or_stack_id = self._get_param("StackName")
        self.cloudformation_backend.delete_stack(name_or_stack_id)
        return EmptyResult()

    def list_exports(self) -> ActionResult:
        token = self._get_param("NextToken")
        exports, next_token = self.cloudformation_backend.list_exports(tokenstr=token)
        result = {"Exports": exports, "NextToken": next_token}
        return ActionResult(result)

    def validate_template(self) -> ActionResult:
        template_body = self._get_param("TemplateBody")
        template_url = self._get_param("TemplateURL")
        if template_url:
            template_body = self._get_stack_from_s3_url(template_url)

        cfn_lint = self.cloudformation_backend.validate_template(template_body)
        if cfn_lint:
            raise ValidationError(cfn_lint[0].message)
        description = ""
        try:
            description = json.loads(template_body)["Description"]
        except (ValueError, KeyError):
            pass
        try:
            description = yaml.load(template_body, Loader=yaml.Loader)["Description"]
        except (ParserError, ScannerError, KeyError):
            pass
        result = {
            "Description": description,
            "Parameters": [],
            "Capabilities": [],
            "DeclaredTransforms": [],
        }
        return ActionResult(result)

    def create_stack_set(self) -> ActionResult:
        stackset_name = self._get_param("StackSetName")
        if not re.match(r"^[a-zA-Z][-a-zA-Z0-9]*$", stackset_name):
            raise ValidationError(
                message=f"1 validation error detected: Value '{stackset_name}' at 'stackSetName' failed to satisfy constraint: Member must satisfy regular expression pattern: [a-zA-Z][-a-zA-Z0-9]*"
            )
        stack_body = self._get_param("TemplateBody")
        template_url = self._get_param("TemplateURL")
        permission_model = self._get_param("PermissionModel")
        parameters_list = self._get_param("Parameters", [])
        admin_role = self._get_param("AdministrationRoleARN")
        exec_role = self._get_param("ExecutionRoleName")
        description = self._get_param("Description")
        tags = dict(
            (item["key"], item["value"]) for item in self._get_param("Tags", [])
        )

        # Copy-Pasta - Hack dict-comprehension
        parameters = dict(
            [
                (parameter["parameter_key"], parameter["parameter_value"])
                for parameter in parameters_list
            ]
        )
        if template_url:
            stack_body = self._get_stack_from_s3_url(template_url)

        stackset = self.cloudformation_backend.create_stack_set(
            name=stackset_name,
            template=stack_body,
            parameters=parameters,
            tags=tags,
            permission_model=permission_model,
            admin_role=admin_role,
            exec_role=exec_role,
            description=description,
        )
        result = {"StackSetId": stackset.id}
        return ActionResult(result)

    def create_stack_instances(self) -> ActionResult:
        stackset_name = self._get_param("StackSetName")
        accounts = self._get_param("Accounts", [])
        regions = self._get_param("Regions", [])
        parameters = self._get_param("ParameterOverrides", [])
        deployment_targets = self._get_params().get("DeploymentTargets")
        if deployment_targets and "OrganizationalUnitIds" in deployment_targets:
            for ou_id in deployment_targets.get("OrganizationalUnitIds", []):
                if not re.match(
                    r"^(ou-[a-z0-9]{4,32}-[a-z0-9]{8,32}|r-[a-z0-9]{4,32})$", ou_id
                ):
                    raise ValidationError(
                        message=f"1 validation error detected: Value '[{ou_id}]' at 'deploymentTargets.organizationalUnitIds' failed to satisfy constraint: Member must satisfy constraint: [Member must have length less than or equal to 68, Member must have length greater than or equal to 6, Member must satisfy regular expression pattern: ^(ou-[a-z0-9]{{4,32}}-[a-z0-9]{{8,32}}|r-[a-z0-9]{{4,32}})$]"
                    )

        operation_id = self.cloudformation_backend.create_stack_instances(
            stackset_name,
            accounts,
            regions,
            parameters,
            deployment_targets=deployment_targets,
        )
        result = {"OperationId": operation_id}
        return ActionResult(result)

    def delete_stack_set(self) -> ActionResult:
        stackset_name = self._get_param("StackSetName")
        self.cloudformation_backend.delete_stack_set(stackset_name)
        return EmptyResult()

    def delete_stack_instances(self) -> ActionResult:
        stackset_name = self._get_param("StackSetName")
        accounts = self._get_param("Accounts", [])
        regions = self._get_param("Regions", [])
        operation = self.cloudformation_backend.delete_stack_instances(
            stackset_name, accounts, regions
        )

        result = {"OperationId": operation.operations[-1]["OperationId"]}
        return ActionResult(result)

    def describe_stack_set(self) -> ActionResult:
        stackset_name = self._get_param("StackSetName")
        stackset = self.cloudformation_backend.describe_stack_set(stackset_name)
        result = {"StackSet": stackset}
        return ActionResult(result)

    def describe_stack_instance(self) -> ActionResult:
        stackset_name = self._get_param("StackSetName")
        account = self._get_param("StackInstanceAccount")
        region = self._get_param("StackInstanceRegion")

        instance = self.cloudformation_backend.describe_stack_instance(
            stackset_name, account, region
        )
        result = {"StackInstance": instance}
        return ActionResult(result)

    def list_stack_sets(self) -> ActionResult:
        stacksets = self.cloudformation_backend.list_stack_sets()
        result = {"Summaries": stacksets}
        return ActionResult(result)

    def list_stack_instances(self) -> ActionResult:
        stackset_name = self._get_param("StackSetName")
        instances = self.cloudformation_backend.list_stack_instances(stackset_name)
        result = {"Summaries": instances}
        return ActionResult(result)

    def list_stack_set_operations(self) -> ActionResult:
        stackset_name = self._get_param("StackSetName")
        operations = self.cloudformation_backend.list_stack_set_operations(
            stackset_name
        )
        result = {"Summaries": operations}
        return ActionResult(result)

    def stop_stack_set_operation(self) -> ActionResult:
        stackset_name = self._get_param("StackSetName")
        operation_id = self._get_param("OperationId")
        self.cloudformation_backend.stop_stack_set_operation(
            stackset_name, operation_id
        )
        return EmptyResult()

    def describe_stack_set_operation(self) -> ActionResult:
        stackset_name = self._get_param("StackSetName")
        operation_id = self._get_param("OperationId")
        stackset, operation = self.cloudformation_backend.describe_stack_set_operation(
            stackset_name, operation_id
        )
        result = {"StackSetOperation": StackSetOperationDTO(operation, stackset)}
        return ActionResult(result)

    def list_stack_set_operation_results(self) -> ActionResult:
        stackset_name = self._get_param("StackSetName")
        operation_id = self._get_param("OperationId")
        operation = self.cloudformation_backend.list_stack_set_operation_results(
            stackset_name, operation_id
        )
        # FIXME: Hardcoded response values come from original XML template.
        result = {
            "Summaries": [
                {
                    "AccountGateResult": {
                        "StatusReason": f"Function not found: arn:aws:lambda:us-west-2:{account}:function:AWSCloudFormationStackSetAccountGate",
                        "Status": "SKIPPED",
                    },
                    "Region": region,
                    "Account": account,
                    "Status": operation["Status"],
                }
                for instance in operation["Instances"]
                for account, region in instance.items()
            ]
        }
        return ActionResult(result)

    def update_stack_set(self) -> ActionResult:
        stackset_name = self._get_param("StackSetName")
        operation_id = self._get_param("OperationId")
        description = self._get_param("Description")
        execution_role = self._get_param("ExecutionRoleName")
        admin_role = self._get_param("AdministrationRoleARN")
        accounts = self._get_param("Accounts", [])
        regions = self._get_param("Regions", [])
        template_body = self._get_param("TemplateBody")
        template_url = self._get_param("TemplateURL")
        if template_url:
            template_body = self._get_stack_from_s3_url(template_url)
        tags = dict(
            (item["key"], item["value"]) for item in self._get_param("Tags", [])
        )
        parameters_list = self._get_param("Parameters", [])

        operation = self.cloudformation_backend.update_stack_set(
            stackset_name=stackset_name,
            template=template_body,
            description=description,
            parameters=parameters_list,
            tags=tags,
            admin_role=admin_role,
            execution_role=execution_role,
            accounts=accounts,
            regions=regions,
            operation_id=operation_id,
        )

        result = {"OperationId": operation["OperationId"]}
        return ActionResult(result)

    def update_stack_instances(self) -> ActionResult:
        stackset_name = self._get_param("StackSetName")
        accounts = self._get_param("Accounts", [])
        regions = self._get_param("Regions", [])
        parameters = self._get_param("ParameterOverrides", [])
        operation = self.cloudformation_backend.update_stack_instances(
            stackset_name, accounts, regions, parameters
        )
        result = {"OperationId": operation["OperationId"]}
        return ActionResult(result)

    def get_stack_policy(self) -> ActionResult:
        stack_name = self._get_param("StackName")
        policy = self.cloudformation_backend.get_stack_policy(stack_name)
        result = {"StackPolicyBody": policy if policy else None}
        return ActionResult(result)

    def set_stack_policy(self) -> ActionResult:
        stack_name = self._get_param("StackName")
        policy_url = self._get_param("StackPolicyURL")
        policy_body = self._get_param("StackPolicyBody")
        if policy_body and policy_url:
            raise ValidationError(
                message="You cannot specify both StackPolicyURL and StackPolicyBody"
            )
        if policy_url:
            try:
                policy_body = self._get_stack_from_s3_url(policy_url)
            except S3ClientError as s3_e:
                raise ValidationError(
                    message=f"S3 error: Access Denied: {s3_e.error_type}"
                )
        self.cloudformation_backend.set_stack_policy(
            stack_name, policy_body=policy_body
        )
        return EmptyResult()
