import json
import re
from typing import Any, Dict, List, Optional, Tuple, Union

import yaml
from yaml.parser import ParserError  # pylint:disable=c-extension-no-member
from yaml.scanner import ScannerError  # pylint:disable=c-extension-no-member

from moto.core.responses import BaseResponse
from moto.s3.exceptions import S3ClientError

from .exceptions import MissingParameterError, ValidationError
from .models import CloudFormationBackend, FakeStack, cloudformation_backends
from .utils import get_stack_from_s3_url, yaml_tag_constructor


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
    template_dict["resourceTypes"] = resources_types
    return template_dict


class CloudFormationResponse(BaseResponse):
    def __init__(self) -> None:
        super().__init__(service_name="cloudformation")

    @property
    def cloudformation_backend(self) -> CloudFormationBackend:
        return cloudformation_backends[self.current_account][self.region]

    @classmethod
    def cfnresponse(cls, *args: Any, **kwargs: Any) -> Any:  # type: ignore[misc]  # pylint: disable=unused-argument
        request, full_url, headers = args
        full_url += "&Action=ProcessCfnResponse"
        return cls.dispatch(request=request, full_url=full_url, headers=headers)

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
            if parameter.keys() >= {"parameter_key", "parameter_value"}:
                result[parameter["parameter_key"]] = parameter["parameter_value"]
            elif (
                parameter.keys() >= {"parameter_key", "use_previous_value"}
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

    def create_stack(self) -> Union[str, Tuple[int, Dict[str, int], str]]:
        stack_name = self._get_param("StackName")
        stack_body = self._get_param("TemplateBody")
        template_url = self._get_param("TemplateURL")
        role_arn = self._get_param("RoleARN")
        enable_termination_protection = self._get_param("EnableTerminationProtection")
        timeout_in_mins = self._get_param("TimeoutInMinutes")
        stack_policy_body = self._get_param("StackPolicyBody")
        parameters_list = self._get_list_prefix("Parameters.member")
        tags = dict(
            (item["key"], item["value"])
            for item in self._get_list_prefix("Tags.member")
        )

        if self.stack_name_exists(new_stack_name=stack_name):
            template = self.response_template(
                CREATE_STACK_NAME_EXISTS_RESPONSE_TEMPLATE
            )
            return 400, {"status": 400}, template.render(name=stack_name)

        parameters = self._get_params_from_list(parameters_list)

        if template_url:
            stack_body = self._get_stack_from_s3_url(template_url)
        stack_notification_arns = self._get_multi_param("NotificationARNs.member")

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
        if self.request_json:
            return json.dumps(
                {
                    "CreateStackResponse": {
                        "CreateStackResult": {"StackId": stack.stack_id}
                    }
                }
            )
        else:
            template = self.response_template(CREATE_STACK_RESPONSE_TEMPLATE)
            return template.render(stack=stack)

    def stack_name_exists(self, new_stack_name: str) -> bool:
        for stack in self.cloudformation_backend.stacks.values():
            if stack.name == new_stack_name:
                return True
        return False

    def create_change_set(self) -> str:
        stack_name = self._get_param("StackName")
        change_set_name = self._get_param("ChangeSetName")
        stack_body = self._get_param("TemplateBody")
        template_url = self._get_param("TemplateURL")
        description = self._get_param("Description")
        role_arn = self._get_param("RoleARN")
        update_or_create = self._get_param("ChangeSetType", "CREATE")
        parameters_list = self._get_list_prefix("Parameters.member")
        tags = dict(
            (item["key"], item["value"])
            for item in self._get_list_prefix("Tags.member")
        )
        parameters = {
            param["parameter_key"]: param["parameter_value"]
            for param in parameters_list
        }
        if template_url:
            stack_body = self._get_stack_from_s3_url(template_url)
        stack_notification_arns = self._get_multi_param("NotificationARNs.member")
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
        if self.request_json:
            return json.dumps(
                {
                    "CreateChangeSetResponse": {
                        "CreateChangeSetResult": {
                            "Id": change_set_id,
                            "StackId": stack_id,
                        }
                    }
                }
            )
        else:
            template = self.response_template(CREATE_CHANGE_SET_RESPONSE_TEMPLATE)
            return template.render(stack_id=stack_id, change_set_id=change_set_id)

    def delete_change_set(self) -> str:
        change_set_name = self._get_param("ChangeSetName")

        self.cloudformation_backend.delete_change_set(change_set_name=change_set_name)
        if self.request_json:
            return json.dumps(
                {"DeleteChangeSetResponse": {"DeleteChangeSetResult": {}}}
            )
        else:
            template = self.response_template(DELETE_CHANGE_SET_RESPONSE_TEMPLATE)
            return template.render()

    def describe_change_set(self) -> str:
        change_set_name = self._get_param("ChangeSetName")
        change_set = self.cloudformation_backend.describe_change_set(
            change_set_name=change_set_name
        )
        template = self.response_template(DESCRIBE_CHANGE_SET_RESPONSE_TEMPLATE)
        return template.render(change_set=change_set)

    def execute_change_set(self) -> str:
        stack_name = self._get_param("StackName")
        change_set_name = self._get_param("ChangeSetName")
        self.cloudformation_backend.execute_change_set(
            stack_name=stack_name, change_set_name=change_set_name
        )
        if self.request_json:
            return json.dumps(
                {"ExecuteChangeSetResponse": {"ExecuteChangeSetResult": {}}}
            )
        else:
            template = self.response_template(EXECUTE_CHANGE_SET_RESPONSE_TEMPLATE)
            return template.render()

    def describe_stacks(self) -> str:
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
        template = self.response_template(DESCRIBE_STACKS_TEMPLATE)
        return template.render(stacks=stacks_resp, next_token=next_token)

    def describe_stack_resource(self) -> str:
        stack_name = self._get_param("StackName")
        logical_resource_id = self._get_param("LogicalResourceId")
        stack, resource = self.cloudformation_backend.describe_stack_resource(
            stack_name, logical_resource_id
        )

        template = self.response_template(DESCRIBE_STACK_RESOURCE_RESPONSE_TEMPLATE)
        return template.render(stack=stack, resource=resource)

    def describe_stack_resources(self) -> str:
        stack_name = self._get_param("StackName")
        stack, resources = self.cloudformation_backend.describe_stack_resources(
            stack_name
        )

        template = self.response_template(DESCRIBE_STACK_RESOURCES_RESPONSE)
        return template.render(stack=stack, resources=resources)

    def describe_stack_events(self) -> str:
        stack_name = self._get_param("StackName")
        events = self.cloudformation_backend.describe_stack_events(stack_name)

        template = self.response_template(DESCRIBE_STACK_EVENTS_RESPONSE)
        return template.render(events=events)

    def list_change_sets(self) -> str:
        change_sets = self.cloudformation_backend.list_change_sets()
        template = self.response_template(LIST_CHANGE_SETS_RESPONSE)
        return template.render(change_sets=change_sets)

    def list_stacks(self) -> str:
        status_filter = self._get_multi_param("StackStatusFilter.member")
        stacks = self.cloudformation_backend.list_stacks(status_filter)
        template = self.response_template(LIST_STACKS_RESPONSE)
        return template.render(stacks=stacks)

    def list_stack_resources(self) -> str:
        stack_name_or_id = self._get_param("StackName")
        resources = self.cloudformation_backend.list_stack_resources(stack_name_or_id)

        template = self.response_template(LIST_STACKS_RESOURCES_RESPONSE)
        return template.render(resources=resources)

    def get_template(self) -> str:
        name_or_stack_id = self.querystring.get("StackName")[0]  # type: ignore[index]
        stack_template = self.cloudformation_backend.get_template(name_or_stack_id)

        if self.request_json:
            return json.dumps(
                {
                    "GetTemplateResponse": {
                        "GetTemplateResult": {
                            "TemplateBody": stack_template,
                            "ResponseMetadata": {
                                "RequestId": "2d06e36c-ac1d-11e0-a958-f9382b6eb86bEXAMPLE"
                            },
                        }
                    }
                }
            )
        else:
            template = self.response_template(GET_TEMPLATE_RESPONSE_TEMPLATE)
            return template.render(stack_template=stack_template)

    def get_template_summary(self) -> str:
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
        template = self.response_template(GET_TEMPLATE_SUMMARY_TEMPLATE)
        return template.render(template_summary=template_summary)

    def _validate_different_update(
        self,
        incoming_params: Optional[List[Dict[str, Any]]],
        stack_body: str,
        old_stack: FakeStack,
    ) -> None:
        if incoming_params and stack_body:
            new_params = self._get_param_values(incoming_params, old_stack.parameters)
            if old_stack.template == stack_body and old_stack.parameters == new_params:
                raise ValidationError(
                    old_stack.name, message="No updates are to be performed."
                )

    def _validate_status(self, stack: FakeStack) -> None:
        if stack.status == "ROLLBACK_COMPLETE":
            raise ValidationError(
                stack.stack_id,
                message=f"Stack:{stack.stack_id} is in ROLLBACK_COMPLETE state and can not be updated.",
            )

    def update_stack(self) -> str:
        stack_name = self._get_param("StackName")
        role_arn = self._get_param("RoleARN")
        template_url = self._get_param("TemplateURL")
        stack_body = self._get_param("TemplateBody")
        stack = self.cloudformation_backend.get_stack(stack_name)
        if self._get_param("UsePreviousTemplate") == "true":
            stack_body = stack.template
        elif not stack_body and template_url:
            stack_body = self._get_stack_from_s3_url(template_url)

        incoming_params = self._get_list_prefix("Parameters.member")
        # boto3 is supposed to let you clear the tags by passing an empty value, but the request body doesn't
        # end up containing anything we can use to differentiate between passing an empty value versus not
        # passing anything. so until that changes, moto won't be able to clear tags, only update them.
        tags: Optional[Dict[str, str]] = dict(
            (item["key"], item["value"])
            for item in self._get_list_prefix("Tags.member")
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
        if self.request_json:
            stack_body = {
                "UpdateStackResponse": {"UpdateStackResult": {"StackId": stack.name}}
            }
            return json.dumps(stack_body)
        else:
            template = self.response_template(UPDATE_STACK_RESPONSE_TEMPLATE)
            return template.render(stack=stack)

    def delete_stack(self) -> str:
        name_or_stack_id = self.querystring.get("StackName")[0]  # type: ignore[index]

        self.cloudformation_backend.delete_stack(name_or_stack_id)
        if self.request_json:
            return json.dumps({"DeleteStackResponse": {"DeleteStackResult": {}}})
        else:
            template = self.response_template(DELETE_STACK_RESPONSE_TEMPLATE)
            return template.render()

    def list_exports(self) -> str:
        token = self._get_param("NextToken")
        exports, next_token = self.cloudformation_backend.list_exports(tokenstr=token)
        template = self.response_template(LIST_EXPORTS_RESPONSE)
        return template.render(exports=exports, next_token=next_token)

    def validate_template(self) -> str:
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
        template = self.response_template(VALIDATE_STACK_RESPONSE_TEMPLATE)
        return template.render(description=description)

    def create_stack_set(self) -> str:
        stackset_name = self._get_param("StackSetName")
        if not re.match(r"^[a-zA-Z][-a-zA-Z0-9]*$", stackset_name):
            raise ValidationError(
                message=f"1 validation error detected: Value '{stackset_name}' at 'stackSetName' failed to satisfy constraint: Member must satisfy regular expression pattern: [a-zA-Z][-a-zA-Z0-9]*"
            )
        stack_body = self._get_param("TemplateBody")
        template_url = self._get_param("TemplateURL")
        permission_model = self._get_param("PermissionModel")
        parameters_list = self._get_list_prefix("Parameters.member")
        admin_role = self._get_param("AdministrationRoleARN")
        exec_role = self._get_param("ExecutionRoleName")
        description = self._get_param("Description")
        tags = dict(
            (item["key"], item["value"])
            for item in self._get_list_prefix("Tags.member")
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
        if self.request_json:
            return json.dumps(
                {
                    "CreateStackSetResponse": {
                        "CreateStackSetResult": {"StackSetId": stackset.id}
                    }
                }
            )
        else:
            template = self.response_template(CREATE_STACK_SET_RESPONSE_TEMPLATE)
            return template.render(stackset=stackset)

    def create_stack_instances(self) -> str:
        stackset_name = self._get_param("StackSetName")
        accounts = self._get_multi_param("Accounts.member")
        regions = self._get_multi_param("Regions.member")
        parameters = self._get_multi_param("ParameterOverrides.member")
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
        template = self.response_template(CREATE_STACK_INSTANCES_TEMPLATE)
        return template.render(operation_id=operation_id)

    def delete_stack_set(self) -> str:
        stackset_name = self._get_param("StackSetName")
        self.cloudformation_backend.delete_stack_set(stackset_name)
        template = self.response_template(DELETE_STACK_SET_RESPONSE_TEMPLATE)
        return template.render()

    def delete_stack_instances(self) -> str:
        stackset_name = self._get_param("StackSetName")
        accounts = self._get_multi_param("Accounts.member")
        regions = self._get_multi_param("Regions.member")
        operation = self.cloudformation_backend.delete_stack_instances(
            stackset_name, accounts, regions
        )

        template = self.response_template(DELETE_STACK_INSTANCES_TEMPLATE)
        return template.render(operation=operation)

    def describe_stack_set(self) -> str:
        stackset_name = self._get_param("StackSetName")
        stackset = self.cloudformation_backend.describe_stack_set(stackset_name)

        if not stackset.execution_role:
            stackset.execution_role = "AWSCloudFormationStackSetExecutionRole"

        template = self.response_template(DESCRIBE_STACK_SET_RESPONSE_TEMPLATE)
        return template.render(stackset=stackset)

    def describe_stack_instance(self) -> str:
        stackset_name = self._get_param("StackSetName")
        account = self._get_param("StackInstanceAccount")
        region = self._get_param("StackInstanceRegion")

        instance = self.cloudformation_backend.describe_stack_instance(
            stackset_name, account, region
        )
        template = self.response_template(DESCRIBE_STACK_INSTANCE_TEMPLATE)
        rendered = template.render(instance=instance)
        return rendered

    def list_stack_sets(self) -> str:
        stacksets = self.cloudformation_backend.list_stack_sets()
        template = self.response_template(LIST_STACK_SETS_TEMPLATE)
        return template.render(stacksets=stacksets)

    def list_stack_instances(self) -> str:
        stackset_name = self._get_param("StackSetName")
        instances = self.cloudformation_backend.list_stack_instances(stackset_name)
        template = self.response_template(LIST_STACK_INSTANCES_TEMPLATE)
        return template.render(instances=instances)

    def list_stack_set_operations(self) -> str:
        stackset_name = self._get_param("StackSetName")
        operations = self.cloudformation_backend.list_stack_set_operations(
            stackset_name
        )
        template = self.response_template(LIST_STACK_SET_OPERATIONS_RESPONSE_TEMPLATE)
        return template.render(operations=operations)

    def stop_stack_set_operation(self) -> str:
        stackset_name = self._get_param("StackSetName")
        operation_id = self._get_param("OperationId")
        self.cloudformation_backend.stop_stack_set_operation(
            stackset_name, operation_id
        )
        template = self.response_template(STOP_STACK_SET_OPERATION_RESPONSE_TEMPLATE)
        return template.render()

    def describe_stack_set_operation(self) -> str:
        stackset_name = self._get_param("StackSetName")
        operation_id = self._get_param("OperationId")
        stackset, operation = self.cloudformation_backend.describe_stack_set_operation(
            stackset_name, operation_id
        )
        template = self.response_template(DESCRIBE_STACKSET_OPERATION_RESPONSE_TEMPLATE)
        return template.render(stackset=stackset, operation=operation)

    def list_stack_set_operation_results(self) -> str:
        stackset_name = self._get_param("StackSetName")
        operation_id = self._get_param("OperationId")
        operation = self.cloudformation_backend.list_stack_set_operation_results(
            stackset_name, operation_id
        )
        template = self.response_template(
            LIST_STACK_SET_OPERATION_RESULTS_RESPONSE_TEMPLATE
        )
        return template.render(operation=operation)

    def update_stack_set(self) -> str:
        stackset_name = self._get_param("StackSetName")
        operation_id = self._get_param("OperationId")
        description = self._get_param("Description")
        execution_role = self._get_param("ExecutionRoleName")
        admin_role = self._get_param("AdministrationRoleARN")
        accounts = self._get_multi_param("Accounts.member")
        regions = self._get_multi_param("Regions.member")
        template_body = self._get_param("TemplateBody")
        template_url = self._get_param("TemplateURL")
        if template_url:
            template_body = self._get_stack_from_s3_url(template_url)
        tags = dict(
            (item["key"], item["value"])
            for item in self._get_list_prefix("Tags.member")
        )
        parameters_list = self._get_list_prefix("Parameters.member")

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

        template = self.response_template(UPDATE_STACK_SET_RESPONSE_TEMPLATE)
        return template.render(operation=operation)

    def update_stack_instances(self) -> str:
        stackset_name = self._get_param("StackSetName")
        accounts = self._get_multi_param("Accounts.member")
        regions = self._get_multi_param("Regions.member")
        parameters = self._get_multi_param("ParameterOverrides.member")
        operation = self.cloudformation_backend.update_stack_instances(
            stackset_name, accounts, regions, parameters
        )
        template = self.response_template(UPDATE_STACK_INSTANCES_RESPONSE_TEMPLATE)
        return template.render(operation=operation)

    def get_stack_policy(self) -> str:
        stack_name = self._get_param("StackName")
        policy = self.cloudformation_backend.get_stack_policy(stack_name)
        template = self.response_template(GET_STACK_POLICY_RESPONSE)
        return template.render(policy=policy)

    def set_stack_policy(self) -> str:
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
        return SET_STACK_POLICY_RESPONSE


VALIDATE_STACK_RESPONSE_TEMPLATE = """<ValidateTemplateResponse>
        <ValidateTemplateResult>
        <Capabilities></Capabilities>
<CapabilitiesReason></CapabilitiesReason>
<DeclaredTransforms></DeclaredTransforms>
<Description>{{ description }}</Description>
<Parameters></Parameters>
</ValidateTemplateResult>
</ValidateTemplateResponse>"""

CREATE_STACK_RESPONSE_TEMPLATE = """<CreateStackResponse>
  <CreateStackResult>
    <StackId>{{ stack.stack_id }}</StackId>
  </CreateStackResult>
  <ResponseMetadata>
    <RequestId>b9b4b068-3a41-11e5-94eb-example</RequestId>
    </ResponseMetadata>
</CreateStackResponse>
"""

CREATE_STACK_NAME_EXISTS_RESPONSE_TEMPLATE = """<ErrorResponse xmlns="http://cloudformation.amazonaws.com/doc/2010-05-15/">
  <Error>
    <Type>Sender</Type>
    <Code>AlreadyExistsException</Code>
    <Message>Stack [{{ name }}] already exists</Message>
  </Error>
  <RequestId>950ff8d7-812a-44b3-bb0c-9b271b954104</RequestId>
</ErrorResponse>"""

UPDATE_STACK_RESPONSE_TEMPLATE = """<UpdateStackResponse xmlns="http://cloudformation.amazonaws.com/doc/2010-05-15/">
  <UpdateStackResult>
    <StackId>{{ stack.stack_id }}</StackId>
  </UpdateStackResult>
  <ResponseMetadata>
    <RequestId>b9b4b068-3a41-11e5-94eb-example</RequestId>
  </ResponseMetadata>
</UpdateStackResponse>
"""

CREATE_CHANGE_SET_RESPONSE_TEMPLATE = """<CreateStackResponse>
  <CreateChangeSetResult>
    <Id>{{change_set_id}}</Id>
    <StackId>{{ stack_id }}</StackId>
  </CreateChangeSetResult>
  <ResponseMetadata>
    <RequestId>{{ request_id }}</RequestId>
  </ResponseMetadata>
</CreateStackResponse>
"""

DELETE_CHANGE_SET_RESPONSE_TEMPLATE = """<DeleteChangeSetResponse>
  <DeleteChangeSetResult>
  </DeleteChangeSetResult>
  <ResponseMetadata>
    <RequestId>3d3200a1-810e-3023-6cc3-example</RequestId>
  </ResponseMetadata>
</DeleteChangeSetResponse>
"""

DESCRIBE_CHANGE_SET_RESPONSE_TEMPLATE = """<DescribeChangeSetResponse>
  <DescribeChangeSetResult>
    <ChangeSetId>{{ change_set.change_set_id }}</ChangeSetId>
    <ChangeSetName>{{ change_set.change_set_name }}</ChangeSetName>
    <StackId>{{ change_set.stack_id }}</StackId>
    <StackName>{{ change_set.stack_name }}</StackName>
    <Description>{{ change_set.description }}</Description>
    <Parameters>
      {% for param_name, param_value in change_set.parameters.items() %}
       <member>
          <ParameterKey>{{ param_name }}</ParameterKey>
          <ParameterValue>{{ param_value }}</ParameterValue>
        </member>
      {% endfor %}
    </Parameters>
    <CreationTime>{{ change_set.creation_time_iso_8601 }}</CreationTime>
    <ExecutionStatus>{{ change_set.execution_status }}</ExecutionStatus>
    <Status>{{ change_set.status }}</Status>
    <StatusReason>{{ change_set.status_reason }}</StatusReason>
    {% if change_set.notification_arns %}
    <NotificationARNs>
      {% for notification_arn in change_set.notification_arns %}
      <member>{{ notification_arn }}</member>
      {% endfor %}
    </NotificationARNs>
    {% else %}
    <NotificationARNs/>
    {% endif %}
    {% if change_set.role_arn %}
    <RoleARN>{{ change_set.role_arn }}</RoleARN>
    {% endif %}
    {% if change_set.changes %}
    <Changes>
      {% for change in change_set.changes %}
      <member>
        <Type>Resource</Type>
        <ResourceChange>
            <Action>{{ change.action }}</Action>
            <LogicalResourceId>{{ change.logical_resource_id }}</LogicalResourceId>
            <ResourceType>{{ change.resource_type }}</ResourceType>
        </ResourceChange>
      </member>
      {% endfor %}
    </Changes>
    {% endif %}
    {% if next_token %}
    <NextToken>{{ next_token }}</NextToken>
    {% endif %}
  </DescribeChangeSetResult>
</DescribeChangeSetResponse>"""

EXECUTE_CHANGE_SET_RESPONSE_TEMPLATE = """<ExecuteChangeSetResponse>
  <ExecuteChangeSetResult>
      <ExecuteChangeSetResult/>
  </ExecuteChangeSetResult>
  <ResponseMetadata>
    <RequestId>{{ request_id }}</RequestId>
  </ResponseMetadata>
</ExecuteChangeSetResponse>
"""

DESCRIBE_STACKS_TEMPLATE = """<DescribeStacksResponse>
  <DescribeStacksResult>
    <Stacks>
      {% for stack in stacks %}
      <member>
        <StackName>{{ stack.name }}</StackName>
        <StackId>{{ stack.stack_id }}</StackId>
        {% if stack.change_set_id %}
        <ChangeSetId>{{ stack.change_set_id }}</stack.timeout_in_minsChangeSetId>
        {% endif %}
        <Description><![CDATA[{{ stack.description }}]]></Description>
        <CreationTime>{{ stack.creation_time_iso_8601 }}</CreationTime>
        <StackStatus>{{ stack.status }}</StackStatus>
        {% if stack.notification_arns %}
        <NotificationARNs>
          {% for notification_arn in stack.notification_arns %}
          <member>{{ notification_arn }}</member>
          {% endfor %}
        </NotificationARNs>
        {% else %}
        <NotificationARNs/>
        {% endif %}
        <DisableRollback>false</DisableRollback>
        {%if stack.stack_outputs %}
        <Outputs>
        {% for output in stack.stack_outputs %}
          <member>
            <OutputKey>{{ output.key }}</OutputKey>
            <OutputValue>{{ output.value }}</OutputValue>
            {% if output.description %}<Description>{{ output.description }}</Description>{% endif %}
          </member>
        {% endfor %}
        </Outputs>
        {% endif %}
        <Parameters>
        {% for param_name, param_value in stack.stack_parameters.items() %}
          <member>
            <ParameterKey>{{ param_name }}</ParameterKey>
            {% if param_name in stack.resource_map.no_echo_parameter_keys %}
                <ParameterValue>****</ParameterValue>
            {% else %}
                <ParameterValue>{{ param_value }}</ParameterValue>
            {% endif %}
          </member>
        {% endfor %}
        </Parameters>
        {% if stack.role_arn %}
        <RoleARN>{{ stack.role_arn }}</RoleARN>
        {% endif %}
        <Tags>
          {% for tag_key, tag_value in stack.tags.items() %}
            <member>
              <Key>{{ tag_key }}</Key>
              <Value>{{ tag_value }}</Value>
            </member>
          {% endfor %}
        </Tags>
        <EnableTerminationProtection>{{ stack.enable_termination_protection }}</EnableTerminationProtection>
        {% if stack.timeout_in_mins %}
        <TimeoutInMinutes>{{ stack.timeout_in_mins }}</TimeoutInMinutes>
        {% endif %}
      </member>
      {% endfor %}
    </Stacks>
    {% if next_token %}
    <NextToken>{{ next_token }}</NextToken>
    {% endif %}
  </DescribeStacksResult>
</DescribeStacksResponse>"""

DESCRIBE_STACK_RESOURCE_RESPONSE_TEMPLATE = """<DescribeStackResourceResponse>
  <DescribeStackResourceResult>
    <StackResourceDetail>
      <StackId>{{ stack.stack_id }}</StackId>
      <StackName>{{ stack.name }}</StackName>
      <LogicalResourceId>{{ resource.logical_resource_id }}</LogicalResourceId>
      <PhysicalResourceId>{{ resource.physical_resource_id }}</PhysicalResourceId>
      <ResourceType>{{ resource.type }}</ResourceType>
      <Timestamp>2010-07-27T22:27:28Z</Timestamp>
      <ResourceStatus>{{ stack.status }}</ResourceStatus>
    </StackResourceDetail>
  </DescribeStackResourceResult>
</DescribeStackResourceResponse>"""

DESCRIBE_STACK_RESOURCES_RESPONSE = """<DescribeStackResourcesResponse>
    <DescribeStackResourcesResult>
      <StackResources>
        {% for resource in resources %}
        <member>
          <StackId>{{ stack.stack_id }}</StackId>
          <StackName>{{ stack.name }}</StackName>
          <LogicalResourceId>{{ resource.logical_resource_id }}</LogicalResourceId>
          <PhysicalResourceId>{{ resource.physical_resource_id }}</PhysicalResourceId>
          <ResourceType>{{ resource.type }}</ResourceType>
          <Timestamp>2010-07-27T22:27:28Z</Timestamp>
          <ResourceStatus>{{ stack.status }}</ResourceStatus>
        </member>
        {% endfor %}
      </StackResources>
    </DescribeStackResourcesResult>
</DescribeStackResourcesResponse>"""

DESCRIBE_STACK_EVENTS_RESPONSE = """<DescribeStackEventsResponse xmlns="http://cloudformation.amazonaws.com/doc/2010-05-15/">
  <DescribeStackEventsResult>
    <StackEvents>
      {% for event in events[::-1] %}
      <member>
        <Timestamp>{{ event.timestamp.strftime('%Y-%m-%dT%H:%M:%S.%fZ') }}</Timestamp>
        <ResourceStatus>{{ event.resource_status }}</ResourceStatus>
        <StackId>{{ event.stack_id }}</StackId>
        <EventId>{{ event.event_id }}</EventId>
        <LogicalResourceId>{{ event.logical_resource_id }}</LogicalResourceId>
        {% if event.resource_status_reason %}<ResourceStatusReason>{{ event.resource_status_reason }}</ResourceStatusReason>{% endif %}
        <StackName>{{ event.stack_name }}</StackName>
        <PhysicalResourceId>{{ event.physical_resource_id }}</PhysicalResourceId>
        {% if event.resource_properties %}<ResourceProperties>{{ event.resource_properties }}</ResourceProperties>{% endif %}
        <ResourceType>{{ event.resource_type }}</ResourceType>
      </member>
      {% endfor %}
    </StackEvents>
  </DescribeStackEventsResult>
  <ResponseMetadata>
    <RequestId>b9b4b068-3a41-11e5-94eb-example</RequestId>
  </ResponseMetadata>
</DescribeStackEventsResponse>"""

LIST_CHANGE_SETS_RESPONSE = """<ListChangeSetsResponse>
 <ListChangeSetsResult>
  <Summaries>
    {% for change_set in change_sets %}
    <member>
        <StackId>{{ change_set.stack_id }}</StackId>
        <StackName>{{ change_set.stack_name }}</StackName>
        <ChangeSetId>{{ change_set.change_set_id }}</ChangeSetId>
        <ChangeSetName>{{ change_set.change_set_name }}</ChangeSetName>
        <ExecutionStatus>{{ change_set.execution_status }}</ExecutionStatus>
        <Status>{{ change_set.status }}</Status>
        <StatusReason>{{ change_set.status_reason }}</StatusReason>
        <CreationTime>2011-05-23T15:47:44Z</CreationTime>
        <Description>{{ change_set.description }}</Description>
    </member>
    {% endfor %}
  </Summaries>
 </ListChangeSetsResult>
</ListChangeSetsResponse>"""

LIST_STACKS_RESPONSE = """<ListStacksResponse>
 <ListStacksResult>
  <StackSummaries>
    {% for stack in stacks %}
    <member>
        <StackId>{{ stack.stack_id }}</StackId>
        <StackStatus>{{ stack.status }}</StackStatus>
        <StackName>{{ stack.name }}</StackName>
        <CreationTime>{{ stack.creation_time_iso_8601 }}</CreationTime>
        <TemplateDescription>{{ stack.description }}</TemplateDescription>
    </member>
    {% endfor %}
  </StackSummaries>
 </ListStacksResult>
</ListStacksResponse>"""

LIST_STACKS_RESOURCES_RESPONSE = """<ListStackResourcesResponse>
  <ListStackResourcesResult>
    <StackResourceSummaries>
      {% for resource in resources %}
      <member>
        <ResourceStatus>CREATE_COMPLETE</ResourceStatus>
        <LogicalResourceId>{{ resource.logical_resource_id }}</LogicalResourceId>
        <LastUpdatedTimestamp>2011-06-21T20:15:58Z</LastUpdatedTimestamp>
        <PhysicalResourceId>{{ resource.physical_resource_id }}</PhysicalResourceId>
        <ResourceType>{{ resource.type }}</ResourceType>
      </member>
      {% endfor %}
    </StackResourceSummaries>
  </ListStackResourcesResult>
  <ResponseMetadata>
    <RequestId>2d06e36c-ac1d-11e0-a958-f9382b6eb86b</RequestId>
  </ResponseMetadata>
</ListStackResourcesResponse>"""

GET_TEMPLATE_RESPONSE_TEMPLATE = """<GetTemplateResponse>
  <GetTemplateResult>
    <TemplateBody>{{ stack_template }}</TemplateBody>
  </GetTemplateResult>
  <ResponseMetadata>
    <RequestId>b9b4b068-3a41-11e5-94eb-example</RequestId>
  </ResponseMetadata>
</GetTemplateResponse>"""

DELETE_STACK_RESPONSE_TEMPLATE = """<DeleteStackResponse>
  <ResponseMetadata>
    <RequestId>5ccc7dcd-744c-11e5-be70-example</RequestId>
  </ResponseMetadata>
</DeleteStackResponse>
"""

LIST_EXPORTS_RESPONSE = """<ListExportsResponse xmlns="http://cloudformation.amazonaws.com/doc/2010-05-15/">
  <ListExportsResult>
    <Exports>
      {% for export in exports %}
      <member>
        <ExportingStackId>{{ export.exporting_stack_id }}</ExportingStackId>
        <Name>{{ export.name }}</Name>
        <Value>{{ export.value }}</Value>
      </member>
      {% endfor %}
    </Exports>
    {% if next_token %}
    <NextToken>{{ next_token }}</NextToken>
    {% endif %}
  </ListExportsResult>
  <ResponseMetadata>
    <RequestId>5ccc7dcd-744c-11e5-be70-example</RequestId>
  </ResponseMetadata>
</ListExportsResponse>"""

CREATE_STACK_SET_RESPONSE_TEMPLATE = """<CreateStackSetResponse xmlns="http://internal.amazon.com/coral/com.amazonaws.maestro.service.v20160713/">
  <CreateStackSetResult>
    <StackSetId>{{ stackset.id }}</StackSetId>
  </CreateStackSetResult>
  <ResponseMetadata>
    <RequestId>f457258c-391d-41d1-861f-example</RequestId>
  </ResponseMetadata>
</CreateStackSetResponse>
"""

DESCRIBE_STACK_SET_RESPONSE_TEMPLATE = """<DescribeStackSetResponse xmlns="http://internal.amazon.com/coral/com.amazonaws.maestro.service.v20160713/">
  <DescribeStackSetResult>
    <StackSet>
      <Capabilities/>
      <StackSetARN>{{ stackset.arn }}</StackSetARN>
      <ExecutionRoleName>{{ stackset.execution_role }}</ExecutionRoleName>
      <AdministrationRoleARN>{{ stackset.admin_role }}</AdministrationRoleARN>
      <StackSetId>{{ stackset.id }}</StackSetId>
      <TemplateBody>{{ stackset.template }}</TemplateBody>
      <StackSetName>{{ stackset.name }}</StackSetName>
        <Parameters>
        {% for param_name, param_value in stackset.parameters.items() %}
          <member>
            <ParameterKey>{{ param_name }}</ParameterKey>
            <ParameterValue>{{ param_value }}</ParameterValue>
          </member>
        {% endfor %}
        </Parameters>
        <Tags>
          {% for tag_key, tag_value in stackset.tags.items() %}
            <member>
              <Key>{{ tag_key }}</Key>
              <Value>{{ tag_value }}</Value>
            </member>
          {% endfor %}
        </Tags>
      <Status>{{ stackset.status }}</Status>
      <PermissionModel>{{ stackset.permission_model }}</PermissionModel>
      {% if stackset.description %}
      <Description>{{ stackset.description }}</Description>
      {% endif %}
      <ManagedExecution><Active>false</Active></ManagedExecution>
    </StackSet>
  </DescribeStackSetResult>
  <ResponseMetadata>
    <RequestId>d8b64e11-5332-46e1-9603-example</RequestId>
  </ResponseMetadata>
</DescribeStackSetResponse>"""

DELETE_STACK_SET_RESPONSE_TEMPLATE = """<DeleteStackSetResponse xmlns="http://internal.amazon.com/coral/com.amazonaws.maestro.service.v20160713/">
  <DeleteStackSetResult/>
  <ResponseMetadata>
    <RequestId>c35ec2d0-d69f-4c4d-9bd7-example</RequestId>
  </ResponseMetadata>
</DeleteStackSetResponse>"""

CREATE_STACK_INSTANCES_TEMPLATE = """<CreateStackInstancesResponse xmlns="http://internal.amazon.com/coral/com.amazonaws.maestro.service.v20160713/">
  <CreateStackInstancesResult>
    <OperationId>{{ operation_id }}</OperationId>
  </CreateStackInstancesResult>
  <ResponseMetadata>
    <RequestId>6b29f7e3-69be-4d32-b374-example</RequestId>
  </ResponseMetadata>
</CreateStackInstancesResponse>
"""

LIST_STACK_INSTANCES_TEMPLATE = """<ListStackInstancesResponse xmlns="http://internal.amazon.com/coral/com.amazonaws.maestro.service.v20160713/">
  <ListStackInstancesResult>
    <Summaries>
    {% for instance in instances %}
      <member>
        <StackId>{{ instance["StackId"] }}</StackId>
        <StackSetId>{{ instance["StackSetId"] }}</StackSetId>
        <Region>{{ instance["Region"] }}</Region>
        <Account>{{ instance["Account"] }}</Account>
        <Status>{{ instance["Status"] }}</Status>
      </member>
    {% endfor %}
    </Summaries>
  </ListStackInstancesResult>
  <ResponseMetadata>
    <RequestId>83c27e73-b498-410f-993c-example</RequestId>
  </ResponseMetadata>
</ListStackInstancesResponse>
"""

DELETE_STACK_INSTANCES_TEMPLATE = """<DeleteStackInstancesResponse xmlns="http://internal.amazon.com/coral/com.amazonaws.maestro.service.v20160713/">
  <DeleteStackInstancesResult>
    <OperationId>{{ operation.OperationId }}</OperationId>
  </DeleteStackInstancesResult>
  <ResponseMetadata>
    <RequestId>e5325090-66f6-4ecd-a531-example</RequestId>
  </ResponseMetadata>
</DeleteStackInstancesResponse>
"""

DESCRIBE_STACK_INSTANCE_TEMPLATE = """<DescribeStackInstanceResponse xmlns="http://internal.amazon.com/coral/com.amazonaws.maestro.service.v20160713/">
  <DescribeStackInstanceResult>
    <StackInstance>
      <StackId>{{ instance["StackId"] }}</StackId>
      <StackSetId>{{ instance["StackSetId"] }}</StackSetId>
    {% if instance["ParameterOverrides"] %}
      <ParameterOverrides>
      {% for override in instance["ParameterOverrides"] %}
      {% if override['ParameterKey'] or override['ParameterValue'] %}
        <member>
          <ParameterKey>{{ override["ParameterKey"] }}</ParameterKey>
          <UsePreviousValue>false</UsePreviousValue>
          <ParameterValue>{{ override["ParameterValue"] }}</ParameterValue>
        </member>
      {% endif %}
      {% endfor %}
      </ParameterOverrides>
    {% else %}
      <ParameterOverrides/>
    {% endif %}
      <Region>{{ instance["Region"] }}</Region>
      <Account>{{ instance["Account"] }}</Account>
      <Status>{{ instance["Status"] }}</Status>
      <StackInstanceStatus>
          <DetailedStatus>{{ instance["StackInstanceStatus"]["DetailedStatus"] }}</DetailedStatus>
      </StackInstanceStatus>
    </StackInstance>
  </DescribeStackInstanceResult>
  <ResponseMetadata>
    <RequestId>c6c7be10-0343-4319-8a25-example</RequestId>
  </ResponseMetadata>
</DescribeStackInstanceResponse>
"""

LIST_STACK_SETS_TEMPLATE = """<ListStackSetsResponse xmlns="http://internal.amazon.com/coral/com.amazonaws.maestro.service.v20160713/">
  <ListStackSetsResult>
    <Summaries>
    {% for stackset in stacksets %}
      <member>
        <StackSetName>{{ stackset.name }}</StackSetName>
        <StackSetId>{{ stackset.id }}</StackSetId>
        <Status>{{ stackset.status }}</Status>
      </member>
    {% endfor %}
    </Summaries>
  </ListStackSetsResult>
  <ResponseMetadata>
    <RequestId>4dcacb73-841e-4ed8-b335-example</RequestId>
  </ResponseMetadata>
</ListStackSetsResponse>
"""

UPDATE_STACK_INSTANCES_RESPONSE_TEMPLATE = """<UpdateStackInstancesResponse xmlns="http://internal.amazon.com/coral/com.amazonaws.maestro.service.v20160713/">
  <UpdateStackInstancesResult>
    <OperationId>{{ operation }}</OperationId>
  </UpdateStackInstancesResult>
  <ResponseMetadata>
    <RequestId>bdbf8e94-19b6-4ce4-af85-example</RequestId>
  </ResponseMetadata>
</UpdateStackInstancesResponse>
"""

UPDATE_STACK_SET_RESPONSE_TEMPLATE = """<UpdateStackSetResponse xmlns="http://internal.amazon.com/coral/com.amazonaws.maestro.service.v20160713/">
  <UpdateStackSetResult>
    <OperationId>{{ operation.OperationId }}</OperationId>
  </UpdateStackSetResult>
  <ResponseMetadata>
    <RequestId>adac907b-17e3-43e6-a254-example</RequestId>
  </ResponseMetadata>
</UpdateStackSetResponse>
"""

LIST_STACK_SET_OPERATIONS_RESPONSE_TEMPLATE = """<ListStackSetOperationsResponse xmlns="http://internal.amazon.com/coral/com.amazonaws.maestro.service.v20160713/">
  <ListStackSetOperationsResult>
    <Summaries>
    {% for operation in operations %}
      <member>
        <CreationTimestamp>{{ operation.CreationTimestamp }}</CreationTimestamp>
        <OperationId>{{ operation.OperationId }}</OperationId>
        <Action>{{ operation.Action }}</Action>
        <EndTimestamp>{{ operation.EndTimestamp }}</EndTimestamp>
        <Status>{{ operation.Status }}</Status>
      </member>
    {% endfor %}
    </Summaries>
  </ListStackSetOperationsResult>
  <ResponseMetadata>
    <RequestId>65b9d9be-08bb-4a43-9a21-example</RequestId>
  </ResponseMetadata>
</ListStackSetOperationsResponse>
"""

STOP_STACK_SET_OPERATION_RESPONSE_TEMPLATE = """<StopStackSetOperationResponse xmlns="http://internal.amazon.com/coral/com.amazonaws.maestro.service.v20160713/">
  <StopStackSetOperationResult/>
  <ResponseMetadata>
    <RequestId>2188554a-07c6-4396-b2c5-example</RequestId>
  </ResponseMetadata>
</StopStackSetOperationResponse>
"""

DESCRIBE_STACKSET_OPERATION_RESPONSE_TEMPLATE = """<DescribeStackSetOperationResponse xmlns="http://internal.amazon.com/coral/com.amazonaws.maestro.service.v20160713/">
  <DescribeStackSetOperationResult>
    <StackSetOperation>
      <ExecutionRoleName>{{ stackset.execution_role }}</ExecutionRoleName>
      <AdministrationRoleARN>{{ stackset.admin_role }}</AdministrationRoleARN>
      <StackSetId>{{ stackset.id }}</StackSetId>
      <CreationTimestamp>{{ operation.CreationTimestamp }}</CreationTimestamp>
      <OperationId>{{ operation.OperationId }}</OperationId>
      <Action>{{ operation.Action }}</Action>
      <OperationPreferences>
        <RegionOrder/>
      </OperationPreferences>
      <EndTimestamp>{{ operation.EndTimestamp }}</EndTimestamp>
      <Status>{{ operation.Status }}</Status>
    </StackSetOperation>
  </DescribeStackSetOperationResult>
  <ResponseMetadata>
    <RequestId>2edc27b6-9ce2-486a-a192-example</RequestId>
  </ResponseMetadata>
</DescribeStackSetOperationResponse>
"""

LIST_STACK_SET_OPERATION_RESULTS_RESPONSE_TEMPLATE = """<ListStackSetOperationResultsResponse xmlns="http://internal.amazon.com/coral/com.amazonaws.maestro.service.v20160713/">
  <ListStackSetOperationResultsResult>
    <Summaries>
    {% for instance in operation.Instances %}
    {% for account, region in instance.items() %}
      <member>
        <AccountGateResult>
          <StatusReason>Function not found: arn:aws:lambda:us-west-2:{{ account }}:function:AWSCloudFormationStackSetAccountGate</StatusReason>
          <Status>SKIPPED</Status>
        </AccountGateResult>
        <Region>{{ region }}</Region>
        <Account>{{ account }}</Account>
        <Status>{{ operation.Status }}</Status>
      </member>
    {% endfor %}
    {% endfor %}
    </Summaries>
  </ListStackSetOperationResultsResult>
  <ResponseMetadata>
    <RequestId>ac05a9ce-5f98-4197-a29b-example</RequestId>
  </ResponseMetadata>
</ListStackSetOperationResultsResponse>
"""

# https://docs.aws.amazon.com/AWSCloudFormation/latest/APIReference/API_GetTemplateSummary.html
# TODO:implement fields: ResourceIdentifierSummaries, Capabilities, CapabilitiesReason
GET_TEMPLATE_SUMMARY_TEMPLATE = """<GetTemplateSummaryResponse xmlns="http://cloudformation.amazonaws.com/doc/2010-05-15/">
  <GetTemplateSummaryResult>
    <Description>{{ template_summary.Description }}</Description>
    {% for resource in template_summary.resourceTypes %}
      <ResourceTypes>
        <member>{{ resource }}</member>
      </ResourceTypes>
    {% endfor %}
    <Parameters>
        {% for k,p in template_summary.get('Parameters',{}).items() %}
        <member>
            <ParameterKey>{{ k }}</ParameterKey> ,
            <Description>{{ p.get('Description', '') }}</Description>,
            {% if p.Default %}
            <DefaultValue>{{ p.Default }}</DefaultValue>
            {% endif %}
            <NoEcho>{{ p.get('NoEcho', False) }}</NoEcho>
            <ParameterType>{{ p.get('Type', 'String') }}</ParameterType>
            <ParameterConstraints>
              {% if p.AllowedValues %}
              <AllowedValues>
                {% for v in p.AllowedValues %}
                <member>{{ v }}</member>
                {% endfor %}
              </AllowedValues>
              {% endif %}
            </ParameterConstraints>
        </member>
        {% endfor %}
    </Parameters>
    <Version>{{ template_summary.AWSTemplateFormatVersion }}</Version>
  </GetTemplateSummaryResult>
  <ResponseMetadata>
    <RequestId>b9b4b068-3a41-11e5-94eb-example</RequestId>
  </ResponseMetadata>
</GetTemplateSummaryResponse>
"""

SET_STACK_POLICY_RESPONSE = """<SetStackPolicyResponse xmlns="http://cloudformation.amazonaws.com/doc/2010-05-15/">
  <ResponseMetadata>
    <RequestId>abe48993-e23f-4167-b703-5b0f1b6aa84f</RequestId>
  </ResponseMetadata>
</SetStackPolicyResponse>"""


GET_STACK_POLICY_RESPONSE = """<GetStackPolicyResponse xmlns="http://cloudformation.amazonaws.com/doc/2010-05-15/">
  <GetStackPolicyResult>
    {% if policy %}
    <StackPolicyBody>{{ policy }}</StackPolicyBody>
    {% endif %}
  </GetStackPolicyResult>
  <ResponseMetadata>
    <RequestId>e9e39eb6-1c05-4f0e-958a-b63f420e0a07</RequestId>
  </ResponseMetadata>
</GetStackPolicyResponse>"""
