"""CodeDeployBackend class with methods for supported APIs."""

import uuid
from enum import Enum
from typing import Any, Dict, List, Optional

from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import BaseModel
from moto.core.utils import iso_8601_datetime_with_milliseconds
from moto.utilities.tagging_service import TaggingService

from .exceptions import (
    ApplicationAlreadyExistsException,
    ApplicationDoesNotExistException,
    ApplicationNameRequiredException,
    DeploymentDoesNotExistException,
    DeploymentGroupAlreadyExistsException,
    DeploymentGroupDoesNotExistException,
    DeploymentGroupNameRequiredException,
)


class Application(BaseModel):
    def __init__(
        self, application_name: str, compute_platform: str, tags: List[Dict[str, str]]
    ):
        self.id = str(uuid.uuid4())
        self.application_name = application_name
        self.compute_platform = compute_platform
        self.tags = tags.copy() if tags else []

        # Boto docs mention that the field should be datetime, but AWS API says number
        self.create_time = iso_8601_datetime_with_milliseconds()

        # these GitHub fields need to be set by the user in the console
        # so will be omitting them for now since they are not required and require console
        # self.github_account_name = ""
        # self.linked_to_github = False

    def to_dict(self) -> Dict[str, Any]:
        return {
            "applicationId": self.id,
            "applicationName": self.application_name,
            "createTime": self.create_time,
            "computePlatform": self.compute_platform,
        }


class CodeDeployDefault(str, Enum):
    # https://docs.aws.amazon.com/codedeploy/latest/userguide/deployment-configurations.html
    AllAtOnce = "AllAtOnce"
    HalfAtATime = "HalfAtATime"
    OneAtATime = "OneAtATime"

    def __str__(self) -> str:
        return f"{self.__class__.__name__}.{self.value}"


class AlarmConfiguration(BaseModel):
    def __init__(
        self,
        alarms: Optional[List[Dict[str, Any]]] = None,
        enabled: Optional[bool] = False,
        ignore_poll_alarm_failure: bool = False,
    ):
        self.alarms = alarms or []
        self.enabled = enabled
        self.ignore_poll_alarm_failure = ignore_poll_alarm_failure


class DeploymentGroup(BaseModel):
    def __init__(
        self,
        application: Application,
        deployment_group_name: str,
        deployment_config_name: Optional[str],
        ec2_tag_filters: Optional[List[Any]],
        on_premises_instance_tag_filters: Optional[List[Any]],
        auto_scaling_groups: Optional[List[str]],
        service_role_arn: str,
        trigger_configurations: Optional[List[Any]],
        alarm_configuration: Optional[AlarmConfiguration],
        auto_rollback_configuration: Optional[Dict[str, Any]],
        outdated_instances_strategy: Optional[str],
        deployment_style: Optional[Any],
        blue_green_deployment_configuration: Optional[Any],
        load_balancer_info: Optional[Any],
        ec2_tag_set: Optional[Any],
        ecs_services: Optional[List[Any]],
        on_premises_tag_set: Optional[Any],
        tags: Optional[List[Dict[str, str]]],
        termination_hook_enabled: Optional[bool],
    ):
        self.application = application
        self.deployment_group_name = deployment_group_name
        self.deployment_config_name = deployment_config_name
        self.ec2_tag_filters = ec2_tag_filters or []
        self.on_premises_instance_tag_filters = on_premises_instance_tag_filters or []
        self.auto_scaling_groups = auto_scaling_groups or []
        self.service_role_arn = service_role_arn
        self.trigger_configurations = trigger_configurations or []
        self.alarm_configuration = alarm_configuration
        self.auto_rollback_configuration = auto_rollback_configuration or {}
        self.outdated_instances_strategy = outdated_instances_strategy
        self.deployment_style = deployment_style or {}
        self.blue_green_deployment_configuration = (
            blue_green_deployment_configuration or {}
        )
        self.load_balancer_info = load_balancer_info or {}
        self.ec2_tag_set = ec2_tag_set or {}
        self.ecs_services = ecs_services or []
        self.on_premises_tag_set = on_premises_tag_set or {}
        self.tags = tags or []
        self.termination_hook_enabled = termination_hook_enabled
        self.deployment_group_id = str(uuid.uuid4())

    def to_dict(self) -> Dict[str, Any]:
        return {
            "applicationName": self.application.application_name,
            "deploymentGroupId": self.deployment_group_id,
            "deploymentGroupName": self.deployment_group_name,
            "deploymentConfigName": str(self.deployment_config_name),
            "ec2TagFilters": self.ec2_tag_filters,
            "onPremisesInstanceTagFilters": self.on_premises_instance_tag_filters,
            "autoScalingGroups": self.auto_scaling_groups,
            "serviceRoleArn": self.service_role_arn,
            "targetRevision": {},  # TODO
            "triggerConfigurations": self.trigger_configurations,
            "alarmConfiguration": {},  # TODO
            "autoRollbackConfiguration": self.auto_rollback_configuration,
            "deploymentStyle": self.deployment_style,
            "outdatedInstancesStrategy": self.outdated_instances_strategy,
            "blueGreenDeploymentConfiguration": self.blue_green_deployment_configuration,
            "loadBalancerInfo": self.load_balancer_info,
            "lastSuccessfulDeployment": {},  # TODO
            "lastAttemptedDeployment": {},  # TODO
            "ec2TagSet": self.ec2_tag_set,
            "onPremisesTagSet": self.on_premises_tag_set,
            "computePlatform": self.application.compute_platform,
            "ecsServices": self.ecs_services,
            "terminationHookEnabled": self.termination_hook_enabled,
        }


class DeploymentInfo(BaseModel):
    def __init__(
        self,
        application: Application,
        deployment_group: DeploymentGroup,
        revision: str,
        deployment_config_name: Optional[str],
        description: Optional[str],
        ignore_application_stop_failures: Optional[bool],
        targetInstances: Optional[Dict[str, Any]],
        auto_rollback_configuration: Optional[Dict[str, Any]],
        update_outdated_instances_only: Optional[bool],
        file_exists_behavior: Optional[str],
        override_alarm_configuration: Optional[AlarmConfiguration],
        creator: Optional[str],
    ):
        self.application = application
        self.deployment_group = deployment_group
        self.deployment_id = str(uuid.uuid4())
        self.application_name = application.application_name
        self.deployment_group_name = deployment_group.deployment_group_name
        self.revision = revision
        self.status = "Created"

        # Boto docs mention that the time fields should be datetime, but AWS API says number
        self.create_time = iso_8601_datetime_with_milliseconds()
        self.start_time = None  # iso_8601_datetime_with_milliseconds()
        self.complete_time = None  # iso_8601_datetime_with_milliseconds()

        # summary of deployment status of the instances in the deployment
        self.deployment_overview = {
            "Pending": 0,
            "InProgress": 0,
            "Succeeded": 0,
            "Failed": 0,
            "Skipped": 0,
            "Ready": 0,
        }
        self.description = description

        # the means by which the deployment was created: {user, autoscaling, codeDeployRollback, CodeDeployAutoUpdate}
        self.creator = "user" if not creator else creator

        self.deployment_config_name = deployment_config_name

        self.ignore_application_stop_failures = ignore_application_stop_failures
        self.target_instances = targetInstances
        self.auto_rollback_configuration = auto_rollback_configuration
        self.update_outdated_instances_only = update_outdated_instances_only
        self.instance_termination_wait_time_started = False

        self.additional_deployment_status_info = ""

        self.file_exists_behavior = file_exists_behavior
        self.deployment_status_messages: list[str] = []
        self.external_id = ""
        self.related_deployments: Dict[str, Any] = {}
        self.override_alarm_configuration = override_alarm_configuration

    def to_dict(self) -> Dict[str, Any]:
        return {
            "applicationName": self.application_name,
            "deploymentGroupName": self.deployment_group_name,
            "deploymentConfigName": str(self.deployment_config_name),
            "deploymentId": self.deployment_id,
            "previousRevision": {},  # TODO
            "revision": self.revision,
            "status": self.status,
            "errorInformation": {},  # TODO
            "createTime": self.create_time,
            "startTime": self.start_time,
            "completeTime": self.complete_time,
            "deploymentOverview": self.deployment_overview,
            "description": self.description,
            "creator": self.creator,
            "ignoreApplicationStopFailures": self.ignore_application_stop_failures,
            "autoRollbackConfiguration": self.auto_rollback_configuration,
            "updateOutdatedInstancesOnly": self.update_outdated_instances_only,
            "rollbackInfo": {},  # TODO information about a deployment rollback
            "deploymentStyle": self.deployment_group.deployment_style,
            "targetInstances": self.target_instances,
            "instanceTerminationWaitTimeStarted": self.instance_termination_wait_time_started,  # TODO
            "blueGreenDeploymentConfiguration": self.deployment_group.blue_green_deployment_configuration,
            "loadBalancerInfo": self.deployment_group.load_balancer_info,
            "additionalDeploymentStatusInfo": self.additional_deployment_status_info,  # TODO
            "fileExistsBehavior": self.file_exists_behavior,
            "deploymentStatusMessages": self.deployment_status_messages,  # TODO
            "computePlatform": self.application.compute_platform,
            "externalId": self.external_id,
            "relatedDeployments": self.related_deployments,  # TODO
            "overrideAlarmConfiguration": self.override_alarm_configuration,
        }


class CodeDeployBackend(BaseBackend):
    """Implementation of CodeDeploy APIs."""

    def __init__(self, region_name: str, account_id: str):
        super().__init__(region_name, account_id)
        self.applications: Dict[str, Application] = {}
        self.deployments: Dict[str, DeploymentInfo] = {}
        self.deployment_groups: Dict[str, Dict[str, DeploymentGroup]] = {}
        self.tagger = TaggingService()

    def get_application(self, application_name: str) -> Application:
        if application_name not in self.applications:
            raise ApplicationDoesNotExistException(
                f"The application {application_name} does not exist with the user or AWS account."
            )
        return self.applications[application_name]

    def batch_get_applications(self, application_names: List[str]) -> List[Application]:
        applications_info = []
        for app_name in application_names:
            app_info = self.get_application(app_name)
            applications_info.append(app_info)

        return applications_info

    def get_deployment(self, deployment_id: str) -> DeploymentInfo:
        if deployment_id not in self.deployments:
            raise DeploymentDoesNotExistException(
                f"The deployment {deployment_id} does not exist with the user or AWS account."
            )
        return self.deployments[deployment_id]

    def get_deployment_group(
        self, application_name: str, deployment_group_name: str
    ) -> DeploymentGroup:
        if application_name not in self.applications:
            raise ApplicationDoesNotExistException(
                f"The application {application_name} does not exist with the user or AWS account."
            )

        # application can also exist but just not associated with a deployment group
        if (
            application_name not in self.deployment_groups
            or deployment_group_name not in self.deployment_groups[application_name]
        ):
            raise DeploymentGroupDoesNotExistException(
                f"The deployment group {deployment_group_name} does not exist with the user or AWS account."
            )
        return self.deployment_groups[application_name][deployment_group_name]

    def batch_get_deployments(self, deployment_ids: List[str]) -> List[DeploymentInfo]:
        deployments = []
        for id in deployment_ids:
            if id in self.deployments:
                deployment_info = self.deployments[id]
                deployments.append(deployment_info)

        return deployments

    def create_application(
        self, application_name: str, compute_platform: str, tags: List[Dict[str, str]]
    ) -> str:
        if application_name in self.applications:
            raise ApplicationAlreadyExistsException(
                f"The application {application_name} already exists with the user or AWS account."
            )

        app = Application(application_name, compute_platform, tags)
        self.applications[app.application_name] = app

        if tags:
            app_arn = f"arn:aws:codedeploy:{self.region_name}:{self.account_id}:application:{application_name}"
            self.tagger.tag_resource(app_arn, tags)

        return app.id

    def create_deployment(
        self,
        application_name: str,
        deployment_group_name: str,
        revision: str,
        deployment_config_name: Optional[str] = None,
        description: Optional[str] = None,
        ignore_application_stop_failures: Optional[bool] = None,
        target_instances: Optional[Any] = None,
        auto_rollback_configuration: Optional[Any] = None,
        update_outdated_instances_only: Optional[bool] = None,
        file_exists_behavior: Optional[str] = None,
        override_alarm_configuration: Optional[Any] = None,
    ) -> str:
        if application_name not in self.applications:
            raise ApplicationDoesNotExistException(
                f"The application {application_name} does not exist with the user or AWS account."
            )

        # Deployment Group Name appears to be optional in create_deployment boto3 documents
        # but seems required in most cases, depending on the deployment type
        # https://docs.aws.amazon.com/AWSCloudFormation/latest/UserGuide/aws-resource-codedeploy-deploymentgroup.html
        # assume required for now

        if deployment_group_name:
            if deployment_group_name not in self.deployment_groups.get(
                application_name, {}
            ):
                raise DeploymentGroupDoesNotExistException(
                    "Deployment group name does not exist."
                )
        else:
            raise DeploymentGroupNameRequiredException(
                "Deployment group name is required."
            )

        if not deployment_config_name:
            # get the deployment from the deployment group if config name is not specified
            deployment_config_name = self.deployment_groups[application_name][
                deployment_group_name
            ].deployment_config_name

        deployment = DeploymentInfo(
            self.applications[application_name],
            self.deployment_groups[application_name][deployment_group_name],
            revision,
            deployment_config_name,
            description,
            ignore_application_stop_failures,
            target_instances,
            auto_rollback_configuration,
            update_outdated_instances_only,
            file_exists_behavior,
            override_alarm_configuration,
            "user",
        )

        self.deployments[deployment.deployment_id] = deployment

        deployment_arn = f"arn:aws:codedeploy:{self.region_name}:{self.account_id}:deployment:{deployment.deployment_id}"
        if self.deployment_groups[application_name][deployment_group_name].tags:
            self.tagger.tag_resource(
                deployment_arn,
                self.deployment_groups[application_name][deployment_group_name].tags,
            )

        return deployment.deployment_id

    # TODO support all optional fields
    def create_deployment_group(
        self,
        application_name: str,
        deployment_group_name: str,
        deployment_config_name: Optional[str],
        ec2_tag_filters: Optional[List[Dict[str, str]]],
        on_premises_instance_tag_filters: Optional[List[Dict[str, str]]],
        auto_scaling_groups: Optional[List[str]],
        service_role_arn: str,
        trigger_configurations: Optional[List[Dict[str, Any]]] = None,
        alarm_configuration: Optional[AlarmConfiguration] = None,
        auto_rollback_configuration: Optional[Dict[str, Any]] = None,
        outdated_instances_strategy: Optional[str] = None,
        deployment_style: Optional[Dict[str, str]] = None,
        blue_green_deployment_configuration: Optional[Dict[str, Any]] = None,
        load_balancer_info: Optional[Dict[str, Any]] = None,
        ec2_tag_set: Optional[Dict[str, Any]] = None,
        ecs_services: Optional[List[Dict[str, str]]] = None,
        on_premises_tag_set: Optional[Dict[str, Any]] = None,
        tags: Optional[List[Dict[str, str]]] = None,
        termination_hook_enabled: Optional[bool] = None,
    ) -> str:
        if application_name not in self.applications:
            raise ApplicationDoesNotExistException(
                f"The application {application_name} does not exist with the user or AWS account."
            )

        if deployment_group_name in self.deployment_groups.get(application_name, {}):
            raise DeploymentGroupAlreadyExistsException(
                f"Deployment group {deployment_group_name} already exists."
            )

        # if deployment_config_name is not specified, use the default
        if not deployment_config_name:
            deployment_config_name = CodeDeployDefault.OneAtATime

        dg = DeploymentGroup(
            self.applications[application_name],
            deployment_group_name,
            deployment_config_name,
            ec2_tag_filters,
            on_premises_instance_tag_filters,
            auto_scaling_groups,
            service_role_arn,
            trigger_configurations,
            alarm_configuration,
            auto_rollback_configuration,
            outdated_instances_strategy,
            deployment_style,
            blue_green_deployment_configuration,
            load_balancer_info,
            ec2_tag_set,
            ecs_services,
            on_premises_tag_set,
            tags,
            termination_hook_enabled,
        )

        if application_name not in self.deployment_groups:
            self.deployment_groups[application_name] = {}
        self.deployment_groups[application_name][dg.deployment_group_name] = dg

        if tags:
            dg_arn = f"arn:aws:codedeploy:{self.region_name}:{self.account_id}:deploymentgroup:{application_name}/{deployment_group_name}"
            self.tagger.tag_resource(dg_arn, tags)

        return dg.deployment_group_id

    # TODO: implement pagination
    def list_applications(self) -> List[str]:
        return list(self.applications.keys())

    # TODO: implement pagination and complete filtering
    def list_deployments(
        self,
        application_name: str,
        deployment_group_name: str,
        external_id: str,
        include_only_statuses: List[str],
        create_time_range: Dict[str, Any],
    ) -> List[str]:
        # Ensure if applicationName is specified, then deploymentGroupName must be specified.
        # If deploymentGroupName is specified, application must be specified else error.
        if application_name and not deployment_group_name:
            raise DeploymentGroupNameRequiredException(
                "If applicationName is specified, then deploymentGroupName must be specified."
            )

        if deployment_group_name and not application_name:
            raise ApplicationNameRequiredException(
                "If deploymentGroupName is specified, applicationName must be specified."
            )

        def matches_filters(deployment: DeploymentInfo) -> bool:
            if application_name and deployment.application_name != application_name:
                return False
            if deployment_group_name:
                if application_name not in self.deployment_groups:
                    return False
                if (
                    deployment_group_name
                    not in self.deployment_groups[application_name]
                ):
                    return False
                if deployment.deployment_group_name != deployment_group_name:
                    return False
                if (
                    include_only_statuses
                    and deployment.status not in include_only_statuses
                ):
                    return False
            return True

        return [
            deployment.deployment_id
            for deployment in self.deployments.values()
            if matches_filters(deployment)
        ]

    # TODO: implement pagination
    def list_deployment_groups(
        self, application_name: str, next_token: str
    ) -> List[str]:
        if application_name not in self.deployment_groups:
            return []

        return [
            deployment_group.deployment_group_name
            for deployment_group in self.deployment_groups[application_name].values()
        ]

    def list_tags_for_resource(
        self, resource_arn: str
    ) -> Dict[str, List[Dict[str, str]]]:
        return self.tagger.list_tags_for_resource(resource_arn)

    def tag_resource(
        self, resource_arn: str, tags: List[Dict[str, str]]
    ) -> Dict[str, Any]:
        self.tagger.tag_resource(resource_arn, tags)
        return {}

    def untag_resource(self, resource_arn: str, tag_keys: List[str]) -> Dict[str, Any]:
        self.tagger.untag_resource_using_names(resource_arn, tag_keys)
        return {}


codedeploy_backends = BackendDict(CodeDeployBackend, "codedeploy")
