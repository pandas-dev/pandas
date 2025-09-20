"""Handles incoming codedeploy requests, invokes methods, returns responses."""

import json

from moto.core.responses import BaseResponse

from .models import CodeDeployBackend, codedeploy_backends


class CodeDeployResponse(BaseResponse):
    """Handler for CodeDeploy requests and responses."""

    def __init__(self) -> None:
        super().__init__(service_name="codedeploy")
        self.default_response_headers = {"Content-Type": "application/json"}

    @property
    def codedeploy_backend(self) -> CodeDeployBackend:
        """Return backend instance specific for this region."""
        return codedeploy_backends[self.current_account][self.region]

    def batch_get_applications(self) -> str:
        application_names = self._get_param("applicationNames")
        applications = self.codedeploy_backend.batch_get_applications(
            application_names=application_names,
        )

        applications_info = {
            "applicationsInfo": [app.to_dict() for app in applications]
        }
        return json.dumps(applications_info)

    def get_application(self) -> str:
        application_name = self._get_param("applicationName")
        application = self.codedeploy_backend.get_application(
            application_name=application_name,
        )

        return json.dumps({"application": application.to_dict()})

    def get_deployment(self) -> str:
        deployment_id = self._get_param("deploymentId")
        deployment = self.codedeploy_backend.get_deployment(
            deployment_id=deployment_id,
        )
        return json.dumps({"deploymentInfo": deployment.to_dict()})

    def get_deployment_group(self) -> str:
        application_name = self._get_param("applicationName")
        deployment_group_name = self._get_param("deploymentGroupName")
        deployment_group = self.codedeploy_backend.get_deployment_group(
            application_name=application_name,
            deployment_group_name=deployment_group_name,
        )
        return json.dumps({"deploymentGroupInfo": deployment_group.to_dict()})

    def batch_get_deployments(self) -> str:
        deployment_ids = self._get_param("deploymentIds")
        deployments = self.codedeploy_backend.batch_get_deployments(
            deployment_ids=deployment_ids,
        )

        deployments_info = {
            "deploymentsInfo": [deployment.to_dict() for deployment in deployments]
        }
        return json.dumps(deployments_info)

    def create_application(self) -> str:
        application_name = self._get_param("applicationName")
        compute_platform = self._get_param("computePlatform")
        tags = self._get_param("tags")
        application_id = self.codedeploy_backend.create_application(
            application_name=application_name,
            compute_platform=compute_platform,
            tags=tags,
        )
        return json.dumps(dict(applicationId=application_id))

    def create_deployment(self) -> str:
        application_name = self._get_param("applicationName")
        deployment_group_name = self._get_param("deploymentGroupName")
        revision = self._get_param("revision")
        deployment_config_name = self._get_param("deploymentConfigName")
        description = self._get_param("description")
        ignore_application_stop_failures = self._get_bool_param(
            "ignoreApplicationStopFailures"
        )
        target_instances = self._get_param("targetInstances")
        auto_rollback_configuration = self._get_param("autoRollbackConfiguration")
        update_outdated_instances_only = self._get_param("updateOutdatedInstancesOnly")
        file_exists_behavior = self._get_param("fileExistsBehavior")
        override_alarm_configuration = self._get_param("overrideAlarmConfiguration")
        deployment_id = self.codedeploy_backend.create_deployment(
            application_name=application_name,
            deployment_group_name=deployment_group_name,
            revision=revision,
            deployment_config_name=deployment_config_name,
            description=description,
            ignore_application_stop_failures=ignore_application_stop_failures,
            target_instances=target_instances,
            auto_rollback_configuration=auto_rollback_configuration,
            update_outdated_instances_only=update_outdated_instances_only,
            file_exists_behavior=file_exists_behavior,
            override_alarm_configuration=override_alarm_configuration,
        )
        return json.dumps(dict(deploymentId=deployment_id))

    def create_deployment_group(self) -> str:
        application_name = self._get_param("applicationName")
        deployment_group_name = self._get_param("deploymentGroupName")
        deployment_config_name = self._get_param("deploymentConfigName")
        ec2_tag_filters = self._get_param("ec2TagFilters")
        on_premises_instance_tag_filters = self._get_param(
            "onPremisesInstanceTagFilters"
        )
        auto_scaling_groups = self._get_param("autoScalingGroups")
        service_role_arn = self._get_param("serviceRoleArn")
        trigger_configurations = self._get_param("triggerConfigurations")
        alarm_configuration = self._get_param("alarmConfiguration")
        auto_rollback_configuration = self._get_param("autoRollbackConfiguration")
        outdated_instances_strategy = self._get_param("outdatedInstancesStrategy")
        deployment_style = self._get_param("deploymentStyle")
        blue_green_deployment_configuration = self._get_param(
            "blueGreenDeploymentConfiguration"
        )
        load_balancer_info = self._get_param("loadBalancerInfo")
        ec2_tag_set = self._get_param("ec2TagSet")
        ecs_services = self._get_param("ecsServices")
        on_premises_tag_set = self._get_param("onPremisesTagSet")
        tags = self._get_param("tags")
        termination_hook_enabled = self._get_param("terminationHookEnabled")
        deployment_group_id = self.codedeploy_backend.create_deployment_group(
            application_name=application_name,
            deployment_group_name=deployment_group_name,
            deployment_config_name=deployment_config_name,
            ec2_tag_filters=ec2_tag_filters,
            on_premises_instance_tag_filters=on_premises_instance_tag_filters,
            auto_scaling_groups=auto_scaling_groups,
            service_role_arn=service_role_arn,
            trigger_configurations=trigger_configurations,
            alarm_configuration=alarm_configuration,
            auto_rollback_configuration=auto_rollback_configuration,
            outdated_instances_strategy=outdated_instances_strategy,
            deployment_style=deployment_style,
            blue_green_deployment_configuration=blue_green_deployment_configuration,
            load_balancer_info=load_balancer_info,
            ec2_tag_set=ec2_tag_set,
            ecs_services=ecs_services,
            on_premises_tag_set=on_premises_tag_set,
            tags=tags,
            termination_hook_enabled=termination_hook_enabled,
        )
        return json.dumps(dict(deploymentGroupId=deployment_group_id))

    def list_applications(self) -> str:
        applications = self.codedeploy_backend.list_applications()
        return json.dumps(dict(applications=applications))

    def list_deployments(self) -> str:
        application_name = self._get_param("applicationName")
        deployment_group_name = self._get_param("deploymentGroupName")
        external_id = self._get_param("externalId")
        include_only_statuses = self._get_param("includeOnlyStatuses")
        create_time_range = self._get_param("createTimeRange")
        deployments = self.codedeploy_backend.list_deployments(
            application_name=application_name,
            deployment_group_name=deployment_group_name,
            external_id=external_id,
            include_only_statuses=include_only_statuses,
            create_time_range=create_time_range,
        )
        return json.dumps(dict(deployments=deployments))

    def list_deployment_groups(self) -> str:
        application_name = self._get_param("applicationName")
        next_token = self._get_param("nextToken", "")
        deployment_groups = self.codedeploy_backend.list_deployment_groups(
            application_name=application_name,
            next_token=next_token,
        )
        return json.dumps(
            dict(
                applicationName=application_name,
                deploymentGroups=deployment_groups,
                nextToken=next_token,
            )
        )

    def list_tags_for_resource(self) -> str:
        """Handler for list_tags_for_resource API call."""
        resource_arn = self._get_param("ResourceArn")
        tags_response = self.codedeploy_backend.list_tags_for_resource(resource_arn)
        return json.dumps(tags_response)

    def tag_resource(self) -> str:
        """Handler for tag_resource API call."""
        resource_arn = self._get_param("ResourceArn")
        tags = self._get_param("Tags")
        response = self.codedeploy_backend.tag_resource(resource_arn, tags)
        return json.dumps(response)

    def untag_resource(self) -> str:
        """Handler for untag_resource API call."""
        resource_arn = self._get_param("ResourceArn")
        tag_keys = self._get_param("TagKeys")
        response = self.codedeploy_backend.untag_resource(resource_arn, tag_keys)
        return json.dumps(response)
