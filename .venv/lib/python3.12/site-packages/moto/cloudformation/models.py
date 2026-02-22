from __future__ import annotations

import json
from collections import OrderedDict
from collections.abc import Iterable
from datetime import timedelta
from typing import Any, Optional, Union

import yaml
from yaml.parser import ParserError
from yaml.scanner import ScannerError

from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import BaseModel, CloudFormationModel
from moto.core.utils import iso_8601_datetime_with_milliseconds, utcnow
from moto.moto_api._internal import mock_random
from moto.organizations.models import OrganizationsBackend, organizations_backends
from moto.sns.models import sns_backends
from moto.utilities.utils import get_partition

from .custom_model import CustomModel
from .exceptions import (
    AlreadyExistsException,
    StackSetNotEmpty,
    StackSetNotFoundException,
    ValidationError,
)
from .parsing import Export, OutputMap, ResourceMap
from .utils import (
    generate_changeset_id,
    generate_stack_id,
    generate_stackset_arn,
    generate_stackset_id,
    get_stack_from_s3_url,
    validate_create_change_set,
    validate_template_cfn_lint,
    yaml_tag_constructor,
)


class StackSet(BaseModel):
    def __init__(
        self,
        stackset_id: str,
        account_id: str,
        name: str,
        template: str,
        region: str,
        description: Optional[str],
        parameters: dict[str, str],
        permission_model: str,
        tags: Optional[dict[str, str]],
        admin_role: Optional[str],
        execution_role: Optional[str],
    ):
        self.id = stackset_id
        self.arn = generate_stackset_arn(stackset_id, region, account_id)
        self.name = name
        self.template = template
        self.description = description
        self.parameters = parameters
        self.tags = tags
        self.administration_role_arn = (
            admin_role
            or f"arn:{get_partition(region)}:iam::{account_id}:role/AWSCloudFormationStackSetAdministrationRole"
        )
        self.execution_role_name = (
            execution_role or "AWSCloudFormationStackSetExecutionRole"
        )
        self.status = "ACTIVE"
        self.instances = StackInstances(
            account_id=account_id,
            region=region,
            template=template,
            parameters=parameters,
            stackset_id=self.id,
            stackset_name=self.name,
        )
        self.stack_instances = self.instances.stack_instances
        self.operations: list[dict[str, Any]] = []
        self.permission_model = permission_model or "SELF_MANAGED"

    @property
    def template_body(self) -> str:
        return self.template

    @property
    def managed_execution(self) -> dict[str, Any]:
        return {"Active": False}

    def _create_operation(
        self,
        operation_id: str,
        action: str,
        status: str,
        accounts: Optional[list[str]] = None,
        regions: Optional[list[str]] = None,
    ) -> dict[str, Any]:
        accounts = accounts or []
        regions = regions or []
        operation = {
            "OperationId": operation_id,
            "Action": action,
            "Status": status,
            "CreationTimestamp": utcnow(),
            "EndTimestamp": utcnow() + timedelta(minutes=2),
            "Instances": [
                {account: region} for account in accounts for region in regions
            ],
        }

        self.operations += [operation]
        return operation

    def get_operation(self, operation_id: str) -> dict[str, Any]:
        for operation in self.operations:
            if operation_id == operation["OperationId"]:
                return operation
        raise ValidationError(operation_id)

    def update_operation(self, operation_id: str, status: str) -> str:
        operation = self.get_operation(operation_id)
        operation["Status"] = status
        return operation_id

    def delete(self) -> None:
        self.status = "DELETED"

    def update(
        self,
        template: str,
        description: str,
        parameters: dict[str, str],
        tags: dict[str, str],
        admin_role: str,
        execution_role: str,
        accounts: list[str],
        regions: list[str],
        operation_id: str,
    ) -> dict[str, Any]:
        self.template = template or self.template
        self.description = description if description is not None else self.description
        self.parameters = parameters or self.parameters
        self.tags = tags or self.tags
        self.administration_role_arn = admin_role or self.administration_role_arn
        self.execution_role_name = execution_role or self.execution_role_name

        if accounts and regions:
            self.update_instances(accounts, regions, self.parameters)  # type: ignore[arg-type]

        operation = self._create_operation(
            operation_id=operation_id,
            action="UPDATE",
            status="SUCCEEDED",
            accounts=accounts,
            regions=regions,
        )
        return operation

    def create_stack_instances(
        self,
        accounts: list[str],
        regions: list[str],
        deployment_targets: Optional[dict[str, Any]],
        parameters: list[dict[str, Any]],
    ) -> str:
        if self.permission_model == "SERVICE_MANAGED":
            if not deployment_targets:
                raise ValidationError(
                    message="StackSets with SERVICE_MANAGED permission model can only have OrganizationalUnit as target"
                )
            elif "OrganizationalUnitIds" not in deployment_targets:
                raise ValidationError(message="OrganizationalUnitIds are required")
        if self.permission_model == "SELF_MANAGED":
            if deployment_targets and "OrganizationalUnitIds" in deployment_targets:
                raise ValidationError(
                    message="StackSets with SELF_MANAGED permission model can only have accounts as target"
                )
        operation_id = str(mock_random.uuid4())
        if not parameters:
            parameters = self.parameters  # type: ignore[assignment]

        self.instances.create_instances(
            accounts,
            regions,
            parameters,
            deployment_targets or {},
            permission_model=self.permission_model,
        )
        self._create_operation(
            operation_id=operation_id,
            action="CREATE",
            status="SUCCEEDED",
            accounts=accounts,
            regions=regions,
        )
        return operation_id

    def delete_stack_instances(self, accounts: list[str], regions: list[str]) -> None:
        operation_id = str(mock_random.uuid4())

        self.instances.delete(accounts, regions)

        self._create_operation(
            operation_id=operation_id,
            action="DELETE",
            status="SUCCEEDED",
            accounts=accounts,
            regions=regions,
        )

    def update_instances(
        self, accounts: list[str], regions: list[str], parameters: list[dict[str, Any]]
    ) -> dict[str, Any]:
        operation_id = str(mock_random.uuid4())

        self.instances.update(accounts, regions, parameters)
        operation = self._create_operation(
            operation_id=operation_id,
            action="UPDATE",
            status="SUCCEEDED",
            accounts=accounts,
            regions=regions,
        )
        return operation


class StackInstance(BaseModel):
    def __init__(
        self,
        account_id: str,
        region_name: str,
        stackset_id: str,
        stack_name: str,
        name: str,
        template: str,
        parameters: Optional[list[dict[str, Any]]],
        permission_model: str,
    ):
        self.account_id = account_id
        self.region_name = region_name
        self.stackset_id = stackset_id
        self.stack_name = stack_name
        self.name = name
        self.template = template
        self.parameters = parameters or []
        self.permission_model = permission_model

        # Incoming parameters can be in two formats: {key: value} or [{"": key, "": value}, ..]
        if isinstance(parameters, dict):
            params = parameters
        elif isinstance(parameters, list):
            params = {p["ParameterKey"]: p["ParameterValue"] for p in parameters}

        if permission_model == "SELF_MANAGED":
            self.stack = cloudformation_backends[account_id][region_name].create_stack(
                name=f"StackSet:{name}", template=template, parameters=params
            )
        else:
            stack_id = generate_stack_id(
                "hiddenstackfor" + self.name, self.region_name, self.account_id
            )
            self.stack = Stack(
                stack_id=stack_id,
                name=self.name,
                template=self.template,
                parameters=params,
                account_id=self.account_id,
                region_name=self.region_name,
                notification_arns=[],
                tags=None,
                role_arn=None,
                cross_stack_resources={},
                enable_termination_protection=False,
            )
            self.stack.create_resources()

    def delete(self) -> None:
        if self.permission_model == "SELF_MANAGED":
            cloudformation_backends[self.account_id][self.region_name].delete_stack(
                self.stack.name
            )
        else:
            # Our stack is hidden - we have to delete it manually
            self.stack.delete()

    def to_dict(self) -> dict[str, Any]:
        return {
            "StackId": generate_stack_id(
                self.stack_name, self.region_name, self.account_id
            ),
            "StackSetId": self.stackset_id,
            "Region": self.region_name,
            "Account": self.account_id,
            "Status": "CURRENT",
            "ParameterOverrides": self.parameters,
            "StackInstanceStatus": {"DetailedStatus": "SUCCEEDED"},
        }


class StackInstances(BaseModel):
    def __init__(
        self,
        account_id: str,
        region: str,
        template: str,
        parameters: dict[str, str],
        stackset_id: str,
        stackset_name: str,
    ):
        self.account_id = account_id
        self.partition = get_partition(region)
        self.template = template
        self.parameters = parameters or {}
        self.stackset_id = stackset_id
        self.stack_name = f"StackSet-{stackset_id}"
        self.stackset_name = stackset_name
        self.stack_instances: list[StackInstance] = []

    @property
    def org_backend(self) -> OrganizationsBackend:
        return organizations_backends[self.account_id][self.partition]

    def create_instances(
        self,
        accounts: list[str],
        regions: list[str],
        parameters: Optional[list[dict[str, Any]]],
        deployment_targets: dict[str, Any],
        permission_model: str,
    ) -> list[dict[str, Any]]:
        targets: list[tuple[str, str]] = []
        all_accounts = self.org_backend.accounts
        requested_ous = deployment_targets.get("OrganizationalUnitIds", [])
        child_ous = [
            ou.id for ou in self.org_backend.ou if ou.parent_id in requested_ous
        ]
        for region in regions:
            for account in accounts:
                targets.append((region, account))
            for ou_id in requested_ous + child_ous:
                for acnt in all_accounts:
                    if acnt.parent_id == ou_id:
                        targets.append((region, acnt.id))

        new_instances = []
        for region, account in targets:
            instance = StackInstance(
                account_id=account,
                region_name=region,
                stackset_id=self.stackset_id,
                stack_name=self.stack_name,
                name=self.stackset_name,
                template=self.template,
                parameters=parameters,
                permission_model=permission_model,
            )
            new_instances.append(instance)
        self.stack_instances += new_instances
        return [i.to_dict() for i in new_instances]

    def update(
        self,
        accounts: list[str],
        regions: list[str],
        parameters: Optional[list[dict[str, Any]]],
    ) -> Any:
        for account in accounts:
            for region in regions:
                instance = self.get_instance(account, region)
                instance.parameters = parameters or []

    def delete(self, accounts: list[str], regions: list[str]) -> None:
        to_delete = [
            i
            for i in self.stack_instances
            if i.region_name in regions and i.account_id in accounts
        ]
        for instance in to_delete:
            instance.delete()
            self.stack_instances.remove(instance)

    def get_instance(self, account: str, region: str) -> StackInstance:  # type: ignore[return]
        for i, instance in enumerate(self.stack_instances):
            if instance.region_name == region and instance.account_id == account:
                return self.stack_instances[i]


class Stack(CloudFormationModel):
    class Meta:
        serialization_aliases = {
            "Parameters": "parameter_list",
        }

    def __init__(
        self,
        stack_id: str,
        name: str,
        template: Union[str, dict[str, Any]],
        parameters: dict[str, str],
        account_id: str,
        region_name: str,
        notification_arns: Optional[list[str]] = None,
        tags: Optional[dict[str, str]] = None,
        role_arn: Optional[str] = None,
        cross_stack_resources: Optional[dict[str, Export]] = None,
        enable_termination_protection: Optional[bool] = False,
        timeout_in_mins: Optional[int] = None,
        stack_policy_body: Optional[str] = None,
    ):
        self.stack_id = stack_id
        self.name = name
        self.account_id = account_id
        self.template = template
        self.template_dict: dict[str, Any]
        if template != {}:
            self._parse_template()
            self.description = self.template_dict.get("Description")
        else:
            self.template_dict = {}
            self.description = None
        self.parameters = parameters
        self.region_name = region_name
        self.notification_arns = notification_arns if notification_arns else []
        self.role_arn = role_arn
        self.tags = tags if tags else {}
        self.events: list[Event] = []
        self.timeout_in_minutes = timeout_in_mins
        self.policy = stack_policy_body or ""

        self.cross_stack_resources: dict[str, Export] = cross_stack_resources or {}
        self.enable_termination_protection: bool = (
            enable_termination_protection or False
        )
        self.resource_map = self._create_resource_map()

        self.custom_resources: dict[str, CustomModel] = {}

        self.output_map = self._create_output_map()
        self.creation_time = utcnow()
        self.status = "CREATE_PENDING"
        self.disable_rollback = False

    def has_template(self, other_template: str) -> bool:
        self._parse_template()
        return self.template_dict == self.parse_template(other_template)

    def has_parameters(self, other_parameters: dict[str, Any]) -> bool:
        return self.parameters == other_parameters

    def _create_resource_map(self) -> ResourceMap:
        resource_map = ResourceMap(
            self.stack_id,
            self.name,
            self.parameters,
            self.tags,
            account_id=self.account_id,
            region_name=self.region_name,
            template=self.template_dict,
            cross_stack_resources=self.cross_stack_resources,
        )
        resource_map.load()
        return resource_map

    def _create_output_map(self) -> OutputMap:
        return OutputMap(self.resource_map, self.template_dict, self.stack_id)

    def _add_stack_event(
        self,
        resource_status: str,
        resource_status_reason: Optional[str] = None,
        resource_properties: Optional[str] = None,
    ) -> None:
        event = Event(
            stack_id=self.stack_id,
            stack_name=self.name,
            logical_resource_id=self.name,
            physical_resource_id=self.stack_id,
            resource_type="AWS::CloudFormation::Stack",
            resource_status=resource_status,
            resource_status_reason=resource_status_reason,
            resource_properties=resource_properties,
        )

        event.sendToSns(self.account_id, self.region_name, self.notification_arns)
        self.events.append(event)

    def _parse_template(self) -> None:
        self.template_dict = self.parse_template(self.template)  # type: ignore[arg-type]

    @staticmethod
    def parse_template(template: str) -> dict[str, Any]:  # type: ignore[misc]
        yaml.add_multi_constructor("", yaml_tag_constructor)
        try:
            return yaml.load(template, Loader=yaml.Loader)
        except (ParserError, ScannerError):
            return json.loads(template)

    @property
    def parameter_list(self) -> list[dict[str, Any]]:  # type: ignore[misc]
        parameters = [
            {
                "ParameterKey": k,
                "ParameterValue": v
                if v not in self.resource_map.no_echo_parameter_keys
                else "****",
            }
            for k, v in self.resource_map.resolved_parameters.items()
        ]
        return parameters

    @property
    def stack_parameters(self) -> dict[str, Any]:  # type: ignore[misc]
        return self.resource_map.resolved_parameters

    @property
    def stack_resources(self) -> Iterable[type[CloudFormationModel]]:
        return self.resource_map.values()

    @property
    def outputs(self) -> Optional[list[dict[str, Any]]]:
        def get_export_name(output_value: Any) -> Optional[str]:
            for export in self.exports:
                if output_value == export.value:
                    return export.name
            return None

        outputs = [
            {
                "OutputKey": o.key,
                "OutputValue": o.value,
                "Description": o.description,
                "ExportName": get_export_name(o.value),
            }
            for o in self.output_map.values()
            if o
        ]
        return outputs if outputs else None

    @property
    def exports(self) -> list[Export]:
        return self.output_map.exports

    @property
    def template_description(self) -> Optional[str]:
        return self.template_dict.get("Description")

    def add_custom_resource(self, custom_resource: CustomModel) -> None:
        self.custom_resources[custom_resource.logical_id] = custom_resource

    def get_custom_resource(self, custom_resource: str) -> CustomModel:
        return self.custom_resources[custom_resource]

    def create_resources(self) -> None:
        self.status = "CREATE_IN_PROGRESS"
        all_resources_ready = self.resource_map.create(self.template_dict)
        # Set the description of the stack
        self.description = self.template_dict.get("Description")
        if all_resources_ready:
            self.mark_creation_complete()

    def verify_readiness(self) -> None:
        if self.resource_map.creation_complete():
            self.mark_creation_complete()

    def mark_creation_complete(self) -> None:
        self.status = "CREATE_COMPLETE"
        self._add_stack_event("CREATE_COMPLETE")

    def update(
        self,
        template: str,
        role_arn: Optional[str] = None,
        parameters: Optional[dict[str, Any]] = None,
        tags: Optional[dict[str, str]] = None,
    ) -> None:
        self._add_stack_event(
            "UPDATE_IN_PROGRESS", resource_status_reason="User Initiated"
        )
        self.template = template
        self._parse_template()
        self.resource_map.update(self.template_dict, parameters)
        self.output_map = self._create_output_map()
        self._add_stack_event("UPDATE_COMPLETE")
        self.status = "UPDATE_COMPLETE"
        self.role_arn = role_arn
        if parameters:
            self.parameters = parameters
        # only overwrite tags if passed
        if tags is not None:
            self.tags = tags
            # TODO: update tags in the resource map

    def delete(self) -> None:
        self._add_stack_event(
            "DELETE_IN_PROGRESS", resource_status_reason="User Initiated"
        )
        self.resource_map.delete()
        self._add_stack_event("DELETE_COMPLETE")
        self.status = "DELETE_COMPLETE"

    @staticmethod
    def cloudformation_type() -> str:
        return "AWS::CloudFormation::Stack"

    @classmethod
    def has_cfn_attr(cls, attr: str) -> bool:
        return True

    @property
    def physical_resource_id(self) -> str:
        return self.name

    @classmethod
    def create_from_cloudformation_json(  # type: ignore[misc]
        cls,
        resource_name: str,
        cloudformation_json: dict[str, Any],
        account_id: str,
        region_name: str,
        **kwargs: Any,
    ) -> Stack:
        cf_backend: CloudFormationBackend = cloudformation_backends[account_id][
            region_name
        ]
        properties = cloudformation_json["Properties"]

        template_body = get_stack_from_s3_url(
            properties["TemplateURL"],
            account_id=account_id,
            partition=get_partition(region_name),
        )
        parameters = properties.get("Parameters", {})

        return cf_backend.create_stack(
            name=resource_name, template=template_body, parameters=parameters
        )

    @classmethod
    def update_from_cloudformation_json(  # type: ignore[misc]
        cls,
        original_resource: Any,
        new_resource_name: str,
        cloudformation_json: Any,
        account_id: str,
        region_name: str,
    ) -> Stack:
        cls.delete_from_cloudformation_json(
            original_resource.name, cloudformation_json, account_id, region_name
        )
        return cls.create_from_cloudformation_json(
            new_resource_name, cloudformation_json, account_id, region_name
        )

    @classmethod
    def delete_from_cloudformation_json(  # type: ignore[misc]
        cls,
        resource_name: str,
        cloudformation_json: dict[str, Any],
        account_id: str,
        region_name: str,
    ) -> None:
        cf_backend: CloudFormationBackend = cloudformation_backends[account_id][
            region_name
        ]
        cf_backend.delete_stack(resource_name)


class Change(BaseModel):
    def __init__(self, action: str, logical_resource_id: str, resource_type: str):
        self.action = action
        self.logical_resource_id = logical_resource_id
        self.resource_type = resource_type


class ChangeSet(BaseModel):
    def __init__(
        self,
        change_set_type: str,
        change_set_id: str,
        change_set_name: str,
        stack: Stack,
        template: str,
        parameters: dict[str, str],
        description: str,
        notification_arns: Optional[list[str]] = None,
        tags: Optional[dict[str, str]] = None,
        role_arn: Optional[str] = None,
    ):
        self.change_set_type = change_set_type
        self.change_set_id = change_set_id
        self.change_set_name = change_set_name

        self.stack = stack
        self.stack_id = self.stack.stack_id
        self.stack_name = self.stack.name
        self.notification_arns = notification_arns
        self.description = description
        self.tags = tags
        self.role_arn = role_arn
        self.template = template
        self.parameters = parameters
        self._parse_template()

        self.creation_time = utcnow()
        self.changes = self.diff()

        self.status: Optional[str] = None
        self.execution_status: Optional[str] = None
        self.status_reason: Optional[str] = None

    def _parse_template(self) -> None:
        yaml.add_multi_constructor("", yaml_tag_constructor)
        try:
            self.template_dict = yaml.load(self.template, Loader=yaml.Loader)
        except (ParserError, ScannerError):
            self.template_dict = json.loads(self.template)

    def diff(self) -> list[Change]:
        changes = []
        resources_by_action = self.stack.resource_map.build_change_set_actions(
            self.template_dict
        )
        for action, resources in resources_by_action.items():
            for resource_name, resource in resources.items():
                changes.append(
                    Change(
                        action=action,
                        logical_resource_id=resource_name,
                        resource_type=resource["ResourceType"],
                    )
                )
        return changes

    def apply(self) -> None:
        self.stack.resource_map.update(self.template_dict, self.parameters)


class Event(BaseModel):
    def __init__(
        self,
        stack_id: str,
        stack_name: str,
        logical_resource_id: str,
        physical_resource_id: str,
        resource_type: str,
        resource_status: str,
        resource_status_reason: Optional[str],
        resource_properties: Optional[str],
    ):
        self.stack_id = stack_id
        self.stack_name = stack_name
        self.logical_resource_id = logical_resource_id
        self.physical_resource_id = physical_resource_id
        self.resource_type = resource_type
        self.resource_status = resource_status
        self.resource_status_reason = resource_status_reason
        self.resource_properties = resource_properties
        self.timestamp = utcnow()
        self.event_id = mock_random.uuid4()
        self.client_request_token = None

    def sendToSns(
        self, account_id: str, region: str, sns_topic_arns: list[str]
    ) -> None:
        message = f"""StackId='{self.stack_id}'
Timestamp='{iso_8601_datetime_with_milliseconds(self.timestamp)}'
EventId='{self.event_id}'
LogicalResourceId='{self.logical_resource_id}'
Namespace='{account_id}'
ResourceProperties='{self.resource_properties}'
ResourceStatus='{self.resource_status}'
ResourceStatusReason='{self.resource_status_reason}'
ResourceType='{self.resource_type}'
StackName='{self.stack_name}'
ClientRequestToken='{self.client_request_token}'"""

        for sns_topic_arn in sns_topic_arns:
            sns_backends[account_id][region].publish(
                message, subject="AWS CloudFormation Notification", arn=sns_topic_arn
            )


def filter_stacks(
    all_stacks: list[Stack], status_filter: Optional[list[str]]
) -> list[Stack]:
    filtered_stacks = []
    if not status_filter:
        return all_stacks
    for stack in all_stacks:
        if stack.status in status_filter:
            filtered_stacks.append(stack)
    return filtered_stacks


class CloudFormationBackend(BaseBackend):
    """
    CustomResources are supported when running Moto in ServerMode.
    Because creating these resources involves running a Lambda-function that informs the MotoServer about the status of the resources, the MotoServer has to be reachable for outside connections.
    This means it has to run inside a Docker-container, or be started using `moto_server -h 0.0.0.0`.
    """

    def __init__(self, region_name: str, account_id: str):
        super().__init__(region_name, account_id)
        self.stacks: dict[str, Stack] = OrderedDict()
        self.stacksets: dict[str, StackSet] = OrderedDict()
        self.deleted_stacks: dict[str, Stack] = {}
        self.exports: dict[str, Export] = OrderedDict()
        self.change_sets: dict[str, ChangeSet] = OrderedDict()

    @staticmethod
    def default_vpc_endpoint_service(
        service_region: str, zones: list[str]
    ) -> list[dict[str, str]]:
        """Default VPC endpoint service."""
        return BaseBackend.default_vpc_endpoint_service_factory(
            service_region, zones, "cloudformation", policy_supported=False
        )

    def _resolve_update_parameters(
        self,
        instance: Union[Stack, StackSet],
        incoming_params: list[dict[str, str]],
    ) -> dict[str, str]:
        parameters = {
            parameter["ParameterKey"]: parameter["ParameterValue"]
            for parameter in incoming_params
            if "ParameterValue" in parameter
        }
        previous = {
            parameter["ParameterKey"]: instance.parameters[parameter["ParameterKey"]]
            for parameter in incoming_params
            if parameter.get("UsePreviousValue", False)
        }
        parameters.update(previous)

        return parameters

    def create_stack_set(
        self,
        name: str,
        template: str,
        parameters: dict[str, str],
        tags: dict[str, str],
        permission_model: str,
        admin_role: Optional[str],
        exec_role: Optional[str],
        description: Optional[str],
    ) -> StackSet:
        """
        The following parameters are not yet implemented: StackId, AdministrationRoleARN, AutoDeployment, ExecutionRoleName, CallAs, ClientRequestToken, ManagedExecution
        """
        stackset_id = generate_stackset_id(name)
        new_stackset = StackSet(
            stackset_id=stackset_id,
            account_id=self.account_id,
            name=name,
            region=self.region_name,
            template=template,
            parameters=parameters,
            description=description,
            tags=tags,
            permission_model=permission_model,
            admin_role=admin_role,
            execution_role=exec_role,
        )
        self.stacksets[stackset_id] = new_stackset
        return new_stackset

    def describe_stack_set(self, name: str) -> StackSet:
        stacksets = self.stacksets.keys()
        if name in stacksets and self.stacksets[name].status != "DELETED":
            return self.stacksets[name]
        for stackset in stacksets:
            if (
                self.stacksets[stackset].name == name
                and self.stacksets[stackset].status != "DELETED"
            ):
                return self.stacksets[stackset]
        raise StackSetNotFoundException(name)

    def delete_stack_set(self, name: str) -> None:
        stackset_to_delete: Optional[StackSet] = None
        if name in self.stacksets:
            stackset_to_delete = self.stacksets[name]
        for stackset in self.stacksets.values():
            if stackset.name == name:
                stackset_to_delete = stackset

        if stackset_to_delete is not None:
            if stackset_to_delete.stack_instances:
                raise StackSetNotEmpty()
            # We don't remove StackSets from the list - they still show up when calling list_stack_sets
            stackset_to_delete.delete()

    def list_stack_sets(self) -> Iterable[StackSet]:
        return self.stacksets.values()

    def list_stack_set_operations(self, stackset_name: str) -> list[dict[str, Any]]:
        stackset = self.describe_stack_set(stackset_name)
        return stackset.operations

    def stop_stack_set_operation(self, stackset_name: str, operation_id: str) -> None:
        stackset = self.describe_stack_set(stackset_name)
        stackset.update_operation(operation_id, "STOPPED")

    def describe_stack_set_operation(
        self, stackset_name: str, operation_id: str
    ) -> tuple[StackSet, dict[str, Any]]:
        stackset = self.describe_stack_set(stackset_name)
        operation = stackset.get_operation(operation_id)
        return stackset, operation

    def list_stack_set_operation_results(
        self, stackset_name: str, operation_id: str
    ) -> dict[str, Any]:
        stackset = self.describe_stack_set(stackset_name)
        return stackset.get_operation(operation_id)

    def create_stack_instances(
        self,
        stackset_name: str,
        accounts: list[str],
        regions: list[str],
        parameters: list[dict[str, str]],
        deployment_targets: Optional[dict[str, Any]],
    ) -> str:
        """
        The following parameters are not yet implemented: DeploymentTargets.AccountFilterType, DeploymentTargets.AccountsUrl, OperationPreferences, CallAs
        """
        stackset = self.describe_stack_set(stackset_name)

        operation_id = stackset.create_stack_instances(
            accounts=accounts,
            regions=regions,
            deployment_targets=deployment_targets,
            parameters=parameters,
        )
        return operation_id

    def update_stack_instances(
        self,
        stackset_name: str,
        accounts: list[str],
        regions: list[str],
        parameters: list[dict[str, Any]],
    ) -> dict[str, Any]:
        """
        Calling this will update the parameters, but the actual resources are not updated
        """
        stack_set = self.describe_stack_set(stackset_name)
        return stack_set.update_instances(accounts, regions, parameters)

    def update_stack_set(
        self,
        stackset_name: str,
        template: str,
        description: str,
        parameters: list[dict[str, str]],
        tags: dict[str, str],
        admin_role: str,
        execution_role: str,
        accounts: list[str],
        regions: list[str],
        operation_id: str,
    ) -> dict[str, Any]:
        stackset = self.describe_stack_set(stackset_name)
        resolved_parameters = self._resolve_update_parameters(
            instance=stackset, incoming_params=parameters
        )
        update = stackset.update(
            template=template,
            description=description,
            parameters=resolved_parameters,
            tags=tags,
            admin_role=admin_role,
            execution_role=execution_role,
            accounts=accounts,
            regions=regions,
            operation_id=operation_id,
        )
        return update

    def delete_stack_instances(
        self, stackset_name: str, accounts: list[str], regions: list[str]
    ) -> StackSet:
        """
        The following parameters are not  yet implemented: DeploymentTargets, OperationPreferences, RetainStacks, OperationId, CallAs
        """
        stackset = self.describe_stack_set(stackset_name)
        stackset.delete_stack_instances(accounts, regions)
        return stackset

    def create_stack(
        self,
        name: str,
        template: str,
        parameters: dict[str, Any],
        notification_arns: Optional[list[str]] = None,
        tags: Optional[dict[str, str]] = None,
        role_arn: Optional[str] = None,
        enable_termination_protection: Optional[bool] = False,
        timeout_in_mins: Optional[int] = None,
        stack_policy_body: Optional[str] = None,
    ) -> Stack:
        """
        The functionality behind EnableTerminationProtection is not yet implemented.
        """
        if name in [stack.name for stack in self.stacks.values()]:
            raise AlreadyExistsException(f"Stack {name} already exists")
        stack_id = generate_stack_id(name, self.region_name, self.account_id)
        new_stack = Stack(
            stack_id=stack_id,
            name=name,
            template=template,
            parameters=parameters,
            account_id=self.account_id,
            region_name=self.region_name,
            notification_arns=notification_arns,
            tags=tags,
            role_arn=role_arn,
            cross_stack_resources=self.exports,
            enable_termination_protection=enable_termination_protection,
            timeout_in_mins=timeout_in_mins,
            stack_policy_body=stack_policy_body,
        )
        self.stacks[stack_id] = new_stack
        self._validate_export_uniqueness(new_stack)
        for export in new_stack.exports:
            self.exports[export.name] = export
        new_stack._add_stack_event(
            "CREATE_IN_PROGRESS", resource_status_reason="User Initiated"
        )
        new_stack.create_resources()
        return new_stack

    def create_change_set(
        self,
        stack_name: str,
        change_set_name: str,
        template: str,
        parameters: dict[str, str],
        description: str,
        change_set_type: str,
        notification_arns: Optional[list[str]] = None,
        tags: Optional[dict[str, str]] = None,
        role_arn: Optional[str] = None,
    ) -> tuple[str, str]:
        validate_create_change_set(change_set_name)
        if change_set_type == "UPDATE":
            for stack in self.stacks.values():
                if stack.name == stack_name:
                    break
            else:
                raise ValidationError(stack_name)
        else:
            stack_id = generate_stack_id(stack_name, self.region_name, self.account_id)
            stack = Stack(
                stack_id=stack_id,
                name=stack_name,
                template={},
                parameters=parameters,
                account_id=self.account_id,
                region_name=self.region_name,
                notification_arns=notification_arns,
                tags=tags,
                role_arn=role_arn,
            )
            self.stacks[stack_id] = stack
            stack.status = "REVIEW_IN_PROGRESS"
            stack._add_stack_event(
                "REVIEW_IN_PROGRESS", resource_status_reason="User Initiated"
            )

        change_set_id = generate_changeset_id(
            change_set_name, self.region_name, self.account_id
        )

        new_change_set = ChangeSet(
            change_set_type=change_set_type,
            change_set_id=change_set_id,
            change_set_name=change_set_name,
            stack=stack,
            template=template,
            parameters=parameters,
            description=description,
            notification_arns=notification_arns,
            tags=tags,
            role_arn=role_arn,
        )
        if (
            change_set_type == "UPDATE"
            and stack.has_template(template)
            and stack.has_parameters(parameters)
        ):
            # Nothing has changed - mark it as such
            new_change_set.status = "FAILED"
            new_change_set.execution_status = "UNAVAILABLE"
            new_change_set.status_reason = "The submitted information didn't contain changes. Submit different information to create a change set."
        else:
            new_change_set.status = "CREATE_COMPLETE"
            new_change_set.execution_status = "AVAILABLE"
        self.change_sets[change_set_id] = new_change_set
        return change_set_id, stack.stack_id

    def delete_change_set(self, change_set_name: str) -> None:
        if change_set_name in self.change_sets:
            # This means arn was passed in
            del self.change_sets[change_set_name]
        else:
            for cs in self.change_sets:
                if self.change_sets[cs].change_set_name == change_set_name:
                    to_delete = cs
                    break
            del self.change_sets[to_delete]

    def describe_change_set(self, change_set_name: str) -> Optional[ChangeSet]:
        change_set = None
        if change_set_name in self.change_sets:
            # This means arn was passed in
            change_set = self.change_sets[change_set_name]
        else:
            for cs in self.change_sets:
                if self.change_sets[cs].change_set_name == change_set_name:
                    change_set = self.change_sets[cs]
        if change_set is None:
            raise ValidationError(change_set_name)
        return change_set

    def execute_change_set(
        self, change_set_name: str, stack_name: Optional[str] = None
    ) -> None:
        if change_set_name in self.change_sets:
            # This means arn was passed in
            change_set = self.change_sets[change_set_name]
        else:
            for cs in self.change_sets:
                if self.change_sets[cs].change_set_name == change_set_name:
                    change_set = self.change_sets[cs]

        if change_set is None:
            raise ValidationError(stack_name)

        stack = self.stacks[change_set.stack_id]
        # TODO: handle execution errors and implement rollback
        if change_set.change_set_type == "CREATE":
            stack._add_stack_event(
                "CREATE_IN_PROGRESS", resource_status_reason="User Initiated"
            )
            change_set.apply()
            stack._add_stack_event("CREATE_COMPLETE")
        else:
            stack._add_stack_event("UPDATE_IN_PROGRESS")
            change_set.apply()
            stack._add_stack_event("UPDATE_COMPLETE")

        # set the execution status of the changeset
        change_set.execution_status = "EXECUTE_COMPLETE"

        # set the status of the stack
        stack.status = f"{change_set.change_set_type}_COMPLETE"
        stack.template = change_set.template

    def describe_stacks(self, name_or_stack_id: str) -> list[Stack]:
        stacks = self.stacks.values()
        if name_or_stack_id:
            for stack in stacks:
                if stack.name == name_or_stack_id or stack.stack_id == name_or_stack_id:
                    return [stack]
            if self.deleted_stacks:
                deleted_stacks = self.deleted_stacks.values()
                for stack in deleted_stacks:
                    if stack.stack_id == name_or_stack_id:
                        return [stack]
            raise ValidationError(name_or_stack_id)
        else:
            return list(stacks)

    def describe_stack_instance(
        self, stack_set_name: str, account_id: str, region: str
    ) -> dict[str, Any]:
        stack_set = self.describe_stack_set(stack_set_name)
        return stack_set.instances.get_instance(account_id, region).to_dict()

    def list_stack_instances(self, stackset_name: str) -> list[dict[str, Any]]:
        """
        Pagination is not yet implemented.
        The parameters StackInstanceAccount/StackInstanceRegion are not yet implemented.
        """
        stack_set = self.describe_stack_set(stackset_name)
        return [i.to_dict() for i in stack_set.instances.stack_instances]

    def list_change_sets(self) -> Iterable[ChangeSet]:
        return self.change_sets.values()

    def list_stacks(self, status_filter: Optional[list[str]] = None) -> list[Stack]:
        total_stacks = list(self.stacks.values()) + list(self.deleted_stacks.values())
        return filter_stacks(total_stacks, status_filter)

    def get_stack(self, name_or_stack_id: str) -> Stack:
        all_stacks = dict(self.deleted_stacks, **self.stacks)
        if name_or_stack_id in all_stacks:
            # Lookup by stack id - deleted stacks incldued
            return all_stacks[name_or_stack_id]
        else:
            # Lookup by stack name - undeleted stacks only
            for stack in self.stacks.values():
                if stack.name == name_or_stack_id:
                    return stack
            raise ValidationError(name_or_stack_id)

    def update_stack(
        self,
        name: str,
        template: str,
        role_arn: Optional[str],
        parameters: list[dict[str, Any]],
        tags: Optional[dict[str, str]],
    ) -> Stack:
        stack = self.get_stack(name)
        resolved_parameters = self._resolve_update_parameters(
            instance=stack, incoming_params=parameters
        )
        stack.update(template, role_arn, parameters=resolved_parameters, tags=tags)
        return stack

    def get_stack_policy(self, stack_name: str) -> str:
        try:
            stack = self.get_stack(stack_name)
        except ValidationError:
            raise ValidationError(message=f"Stack: {stack_name} does not exist")
        return stack.policy

    def set_stack_policy(self, stack_name: str, policy_body: str) -> None:
        """
        Note that Moto does no validation/parsing/enforcement of this policy - we simply persist it.
        """
        try:
            stack = self.get_stack(stack_name)
        except ValidationError:
            raise ValidationError(message=f"Stack: {stack_name} does not exist")
        stack.policy = policy_body

    def describe_stack_resource(
        self, stack_name: str, logical_resource_id: str
    ) -> tuple[Stack, type[CloudFormationModel]]:
        stack = self.get_stack(stack_name)

        for stack_resource in stack.stack_resources:
            if stack_resource.logical_resource_id == logical_resource_id:  # type: ignore[attr-defined]
                return stack, stack_resource

        message = (
            f"Resource {logical_resource_id} does not exist for stack {stack_name}"
        )
        raise ValidationError(stack_name, message)

    def describe_stack_resources(
        self, stack_name: str, logical_resource_id: str
    ) -> tuple[Stack, Iterable[type[CloudFormationModel]]]:
        stack = self.get_stack(stack_name)

        if logical_resource_id is not None:
            for res in stack.stack_resources:
                if res.logical_resource_id == logical_resource_id:  # type: ignore[attr-defined]
                    return stack, [res]
            return stack, []

        return stack, stack.stack_resources

    def list_stack_resources(
        self, stack_name_or_id: str
    ) -> Iterable[type[CloudFormationModel]]:
        stack = self.get_stack(stack_name_or_id)
        return stack.stack_resources

    def delete_stack(self, name_or_stack_id: str) -> None:
        if name_or_stack_id in self.stacks:
            # Delete by stack id
            stack = self.stacks.pop(name_or_stack_id)
            export_names = [export.name for export in stack.exports]
            stack.delete()
            self.deleted_stacks[stack.stack_id] = stack
            for export_name in export_names:
                self.exports.pop(export_name)
            self.stacks.pop(name_or_stack_id, None)
        else:
            # Delete by stack name
            for stack in list(self.stacks.values()):
                if stack.name == name_or_stack_id:
                    self.delete_stack(stack.stack_id)

    def list_exports(
        self, tokenstr: Optional[str]
    ) -> tuple[list[Export], Optional[str]]:
        all_exports = list(self.exports.values())
        if tokenstr is None:
            exports = all_exports[0:100]
            next_token = "100" if len(all_exports) > 100 else None
        else:
            token = int(tokenstr)
            exports = all_exports[token : token + 100]
            next_token = str(token + 100) if len(all_exports) > token + 100 else None
        return exports, next_token

    def describe_stack_events(self, stack_name: str) -> list[Event]:
        return self.get_stack(stack_name).events

    def get_template(self, name_or_stack_id: str) -> Union[str, dict[str, Any]]:
        return self.get_stack(name_or_stack_id).template

    def validate_template(self, template: str) -> list[Any]:
        return validate_template_cfn_lint(template)

    def _validate_export_uniqueness(self, stack: Stack) -> None:
        new_stack_export_names = [x.name for x in stack.exports]
        export_names = self.exports.keys()
        if not set(export_names).isdisjoint(new_stack_export_names):
            raise ValidationError(
                stack.stack_id,
                message="Export names must be unique across a given region",
            )


cloudformation_backends = BackendDict(CloudFormationBackend, "cloudformation")
