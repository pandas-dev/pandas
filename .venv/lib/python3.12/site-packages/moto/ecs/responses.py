from __future__ import annotations

import json

from moto.core.responses import ActionResult, BaseResponse, EmptyResult
from moto.core.utils import camelcase_to_underscores

from .models import EC2ContainerServiceBackend, ecs_backends

DESCRIBE_CLUSTERS_TAGS = "DescribeClustersResponse.clusters.Cluster.tags"
DESCRIBE_SERVICES_TAGS = "DescribeServicesResponse.services.Service.tags"
DESCRIBE_TASKS_TAGS = "DescribeTasksResponse.tasks.Task.tags"
DESCRIBE_TASK_SETS_TAGS = "DescribeTaskSetsResponse.taskSets.TaskSet.tags"


# ContainerInstance.Attribute is a multi-field structure (name, value, targetType, targetId),
# but is currently stored as a dictionary (key=name, value=value) on the Moto ECS backend.
# For now, we use this transformer to convert the attributes dictionary into a data structure
# that matches the AWS spec.
# TODO: model ContainerInstance.Attribute properly on the backend and remove this transformer.
def format_attributes(attributes: dict[str, str] | None) -> list[dict[str, str]] | None:
    if attributes is None:
        return None
    formatted = [{"name": name, "value": value} for name, value in attributes.items()]
    return formatted


class EC2ContainerServiceResponse(BaseResponse):
    RESPONSE_KEY_PATH_TO_TRANSFORMER = {
        "DescribeContainerInstancesResponse.containerInstances.ContainerInstance.attributes": format_attributes,
        "RegisterContainerInstanceResponse.containerInstance.attributes": format_attributes,
    }

    def __init__(self) -> None:
        super().__init__(service_name="ecs")
        self.automated_parameter_parsing = True

    @property
    def ecs_backend(self) -> EC2ContainerServiceBackend:
        return ecs_backends[self.current_account][self.region]

    def create_capacity_provider(self) -> ActionResult:
        name = self._get_param("name")
        asg_provider = self._get_param("autoScalingGroupProvider")
        tags = self._get_param("tags")
        provider = self.ecs_backend.create_capacity_provider(name, asg_provider, tags)
        return ActionResult({"capacityProvider": provider})

    def create_cluster(self) -> ActionResult:
        cluster_name = self._get_param("clusterName")
        tags = self._get_param("tags")
        settings = self._get_param("settings")
        configuration = self._get_param("configuration")
        capacity_providers = self._get_param("capacityProviders")
        default_capacity_provider_strategy = self._get_param(
            "defaultCapacityProviderStrategy"
        )
        service_connect_defaults = self._get_param("serviceConnectDefaults")
        if cluster_name is None:
            cluster_name = "default"
        cluster = self.ecs_backend.create_cluster(
            cluster_name,
            tags,
            settings,
            configuration,
            capacity_providers,
            default_capacity_provider_strategy,
            service_connect_defaults=service_connect_defaults,
        )
        return ActionResult({"cluster": cluster})

    def list_clusters(self) -> ActionResult:
        cluster_arns = self.ecs_backend.list_clusters()
        return ActionResult({"clusterArns": cluster_arns})

    def update_cluster(self) -> ActionResult:
        cluster_name = self._get_param("cluster")
        settings = self._get_param("settings")
        configuration = self._get_param("configuration")
        service_connect_defaults = self._get_param("serviceConnectDefaults")
        cluster = self.ecs_backend.update_cluster(
            cluster_name=cluster_name,
            cluster_settings=settings,
            configuration=configuration,
            service_connect_defaults=service_connect_defaults,
        )
        return ActionResult({"cluster": cluster})

    def put_cluster_capacity_providers(self) -> ActionResult:
        cluster_name = self._get_param("cluster")
        capacity_providers = self._get_param("capacityProviders")
        default_capacity_provider_strategy = self._get_param(
            "defaultCapacityProviderStrategy"
        )
        cluster = self.ecs_backend.put_cluster_capacity_providers(
            cluster_name, capacity_providers, default_capacity_provider_strategy
        )
        return ActionResult({"cluster": cluster})

    def delete_capacity_provider(self) -> ActionResult:
        name = self._get_param("capacityProvider")
        provider = self.ecs_backend.delete_capacity_provider(name)
        return ActionResult({"capacityProvider": provider})

    def update_capacity_provider(self) -> ActionResult:
        name = self._get_param("name")
        asg_provider = self._get_param("autoScalingGroupProvider")
        provider = self.ecs_backend.update_capacity_provider(name, asg_provider)
        return ActionResult({"capacityProvider": provider})

    def describe_capacity_providers(self) -> ActionResult:
        names = self._get_param("capacityProviders")
        providers, failures = self.ecs_backend.describe_capacity_providers(names)
        return ActionResult(
            {
                "capacityProviders": providers,
                "failures": failures,
            }
        )

    def describe_clusters(self) -> ActionResult:
        cluster_names = self._get_param("clusters", [])
        clusters, failures = self.ecs_backend.describe_clusters(cluster_names)
        include = self._get_param("include", [])
        if "TAGS" in include:
            self._include_in_response(DESCRIBE_CLUSTERS_TAGS)
        else:
            self._exclude_from_response(DESCRIBE_CLUSTERS_TAGS)
        result = {
            "clusters": clusters,
            "failures": failures,
        }
        return ActionResult(result)

    def delete_cluster(self) -> ActionResult:
        cluster_str = self._get_param("cluster")
        cluster = self.ecs_backend.delete_cluster(cluster_str)
        return ActionResult({"cluster": cluster})

    def register_task_definition(self) -> ActionResult:
        family = self._get_param("family")
        container_definitions = self._get_param("containerDefinitions")
        volumes = self._get_param("volumes")
        tags = self._get_param("tags")
        network_mode = self._get_param("networkMode")
        placement_constraints = self._get_param("placementConstraints")
        requires_compatibilities = self._get_param("requiresCompatibilities")
        cpu = self._get_param("cpu")
        memory = self._get_param("memory")
        task_role_arn = self._get_param("taskRoleArn")
        execution_role_arn = self._get_param("executionRoleArn")
        proxy_configuration = self._get_param("proxyConfiguration")
        inference_accelerators = self._get_param("inferenceAccelerators")
        runtime_platform = self._get_param("runtimePlatform")
        ipc_mode = self._get_param("ipcMode")
        pid_mode = self._get_param("pidMode")
        ephemeral_storage = self._get_param("ephemeralStorage")

        task_definition = self.ecs_backend.register_task_definition(
            family,
            container_definitions,
            volumes=volumes,
            network_mode=network_mode,
            tags=tags,
            placement_constraints=placement_constraints,
            requires_compatibilities=requires_compatibilities,
            cpu=cpu,
            memory=memory,
            task_role_arn=task_role_arn,
            execution_role_arn=execution_role_arn,
            proxy_configuration=proxy_configuration,
            inference_accelerators=inference_accelerators,
            runtime_platform=runtime_platform,
            ipc_mode=ipc_mode,
            pid_mode=pid_mode,
            ephemeral_storage=ephemeral_storage,
        )
        return ActionResult({"taskDefinition": task_definition})

    def list_task_definitions(self) -> ActionResult:
        family_prefix = self._get_param("familyPrefix")
        status = self._get_param("status", "ACTIVE")
        task_definition_arns = self.ecs_backend.list_task_definitions(
            family_prefix, status
        )
        return ActionResult({"taskDefinitionArns": task_definition_arns})

    def describe_task_definition(self) -> ActionResult:
        task_definition_str = self._get_param("taskDefinition")
        task_definition = self.ecs_backend.describe_task_definition(task_definition_str)
        resp = {"taskDefinition": task_definition, "failures": []}
        if "TAGS" in self._get_param("include", []):
            resp["tags"] = self.ecs_backend.list_tags_for_resource(task_definition.arn)
        return ActionResult(resp)

    def deregister_task_definition(self) -> ActionResult:
        task_definition_str = self._get_param("taskDefinition")
        task_definition = self.ecs_backend.deregister_task_definition(
            task_definition_str
        )
        return ActionResult({"taskDefinition": task_definition})

    def run_task(self) -> ActionResult:
        cluster_str = self._get_param("cluster", "default")
        overrides = self._get_param("overrides")
        task_definition_str = self._get_param("taskDefinition")
        count = self._get_int_param("count", 1)
        started_by = self._get_param("startedBy")
        tags = self._get_param("tags")
        launch_type = self._get_param("launchType")
        network_configuration = self._get_param("networkConfiguration")
        group = self._get_param("group")
        platform_version = self._get_param("platformVersion")
        tasks = self.ecs_backend.run_task(
            cluster_str,
            task_definition_str,
            count,
            overrides,
            started_by,
            tags,
            launch_type,
            network_configuration,
            group,
            platform_version,
        )
        return ActionResult({"tasks": tasks, "failures": []})

    def describe_tasks(self) -> ActionResult:
        cluster = self._get_param("cluster", "default")
        task_ids = self._get_param("tasks", [])
        tasks = self.ecs_backend.describe_tasks(cluster, task_ids)
        include = self._get_param("include", [])
        if "TAGS" in include:
            self._include_in_response(DESCRIBE_TASKS_TAGS)
        else:
            self._exclude_from_response(DESCRIBE_TASKS_TAGS)
        result = {
            "tasks": tasks,
            "failures": [],
        }
        return ActionResult(result)

    def start_task(self) -> ActionResult:
        cluster_str = self._get_param("cluster", "default")
        overrides = self._get_param("overrides")
        task_definition_str = self._get_param("taskDefinition")
        container_instances = self._get_param("containerInstances")
        started_by = self._get_param("startedBy")
        tags = self._get_param("tags")
        group = self._get_param("group")
        tasks = self.ecs_backend.start_task(
            cluster_str,
            task_definition_str,
            container_instances,
            overrides,
            started_by,
            tags,
            group,
        )
        return ActionResult({"tasks": tasks, "failures": []})

    def list_tasks(self) -> ActionResult:
        cluster_str = self._get_param("cluster", "default")
        container_instance = self._get_param("containerInstance")
        family = self._get_param("family")
        started_by = self._get_param("startedBy")
        service_name = self._get_param("serviceName")
        desired_status = self._get_param("desiredStatus")
        tasks = self.ecs_backend.list_tasks(
            cluster_str,
            container_instance,
            family,
            started_by,
            service_name,
            desired_status,
        )
        return ActionResult({"taskArns": [t.task_arn for t in tasks]})

    def stop_task(self) -> ActionResult:
        cluster_str = self._get_param("cluster", "default")
        task = self._get_param("task")
        reason = self._get_param("reason")
        task = self.ecs_backend.stop_task(cluster_str, task, reason)
        return ActionResult({"task": task})

    def create_service(self) -> ActionResult:
        cluster_str = self._get_param("cluster", "default")
        service_name = self._get_param("serviceName")
        task_definition_str = self._get_param("taskDefinition")
        desired_count = self._get_int_param("desiredCount")
        load_balancers = self._get_param("loadBalancers")
        scheduling_strategy = self._get_param("schedulingStrategy")
        service_registries = self._get_param("serviceRegistries")
        tags = self._get_param("tags")
        deployment_controller = self._get_param("deploymentController")
        launch_type = self._get_param("launchType")
        platform_version = self._get_param("platformVersion")
        propagate_tags = self._get_param("propagateTags") or "NONE"
        network_configuration = self._get_param("networkConfiguration")
        role = self._get_param("role")
        service = self.ecs_backend.create_service(
            cluster_str,
            service_name,
            desired_count,
            task_definition_str,
            load_balancers,
            scheduling_strategy,
            tags,
            deployment_controller,
            launch_type,
            network_configuration=network_configuration,
            service_registries=service_registries,
            platform_version=platform_version,
            propagate_tags=propagate_tags,
            role_arn=role,
        )
        return ActionResult({"service": service})

    def list_services(self) -> ActionResult:
        cluster_str = self._get_param("cluster", "default")
        scheduling_strategy = self._get_param("schedulingStrategy")
        launch_type = self._get_param("launchType")
        service_arns = self.ecs_backend.list_services(
            cluster_str, scheduling_strategy, launch_type=launch_type
        )
        return ActionResult({"serviceArns": service_arns})

    def describe_services(self) -> ActionResult:
        cluster_str = self._get_param("cluster", "default")
        service_names = self._get_param("services", [])
        services, failures = self.ecs_backend.describe_services(
            cluster_str, service_names
        )
        include = self._get_param("include", [])
        if "TAGS" in include:
            self._include_in_response(DESCRIBE_SERVICES_TAGS)
        else:
            self._exclude_from_response(DESCRIBE_SERVICES_TAGS)
        result = {
            "services": services,
            "failures": failures,
        }
        return ActionResult(result)

    def update_service(self) -> ActionResult:
        service_properties = self._get_params()
        parsed_props = {}
        for k, v in service_properties.items():
            parsed_props[camelcase_to_underscores(k)] = v

        service = self.ecs_backend.update_service(parsed_props)
        return ActionResult({"service": service})

    def delete_service(self) -> ActionResult:
        service_name = self._get_param("service")
        cluster_name = self._get_param("cluster", "default")
        force = self._get_param("force", False)
        service = self.ecs_backend.delete_service(cluster_name, service_name, force)
        return ActionResult({"service": service})

    def register_container_instance(self) -> ActionResult:
        cluster_str = self._get_param("cluster", "default")
        instance_identity_document_str = self._get_param("instanceIdentityDocument")
        instance_identity_document = json.loads(instance_identity_document_str)
        ec2_instance_id = instance_identity_document["instanceId"]
        container_instance = self.ecs_backend.register_container_instance(
            cluster_str, ec2_instance_id
        )
        return ActionResult({"containerInstance": container_instance})

    def deregister_container_instance(self) -> ActionResult:
        cluster_str = self._get_param("cluster", "default")
        container_instance_str = self._get_param("containerInstance")
        force = self._get_param("force")
        container_instance = self.ecs_backend.deregister_container_instance(
            cluster_str, container_instance_str, force
        )
        return ActionResult({"containerInstance": container_instance})

    def list_container_instances(self) -> ActionResult:
        cluster_str = self._get_param("cluster", "default")
        container_instance_arns = self.ecs_backend.list_container_instances(cluster_str)
        return ActionResult({"containerInstanceArns": container_instance_arns})

    def describe_container_instances(self) -> ActionResult:
        cluster_str = self._get_param("cluster", "default")
        list_container_instance_arns = self._get_param("containerInstances")
        container_instances, failures = self.ecs_backend.describe_container_instances(
            cluster_str, list_container_instance_arns
        )
        return ActionResult(
            {
                "failures": failures,
                "containerInstances": container_instances,
            }
        )

    def update_container_instances_state(self) -> ActionResult:
        cluster_str = self._get_param("cluster", "default")
        list_container_instance_arns = self._get_param("containerInstances")
        status_str = self._get_param("status")
        (
            container_instances,
            failures,
        ) = self.ecs_backend.update_container_instances_state(
            cluster_str, list_container_instance_arns, status_str
        )
        return ActionResult(
            {
                "failures": failures,
                "containerInstances": container_instances,
            }
        )

    def put_attributes(self) -> ActionResult:
        cluster_name = self._get_param("cluster")
        attributes = self._get_param("attributes", [])

        self.ecs_backend.put_attributes(cluster_name, attributes)

        return ActionResult({"attributes": attributes})

    def list_attributes(self) -> ActionResult:
        cluster_name = self._get_param("cluster")
        attr_name = self._get_param("attributeName")
        attr_value = self._get_param("attributeValue")
        target_type = self._get_param("targetType")

        results = self.ecs_backend.list_attributes(
            target_type, cluster_name, attr_name, attr_value
        )
        # Result will be [item will be {0 cluster_name, 1 arn, 2 name, 3 value}]

        formatted_results = []
        for _, arn, name, value in results:
            tmp_result = {"name": name, "targetId": arn}
            if value is not None:
                tmp_result["value"] = value
            formatted_results.append(tmp_result)

        return ActionResult({"attributes": formatted_results})

    def delete_attributes(self) -> ActionResult:
        cluster_name = self._get_param("cluster", "default")
        attributes = self._get_param("attributes")

        self.ecs_backend.delete_attributes(cluster_name, attributes)

        return ActionResult({"attributes": attributes})

    def discover_poll_endpoint(self) -> ActionResult:
        # Here are the arguments, this api is used by the ecs client so obviously no decent
        # documentation. Hence I've responded with valid but useless data
        # cluster_name = self._get_param('cluster')
        # instance = self._get_param('containerInstance')
        return ActionResult(
            {"endpoint": "http://localhost", "telemetryEndpoint": "http://localhost"}
        )

    def list_task_definition_families(self) -> ActionResult:
        family_prefix = self._get_param("familyPrefix")
        results = self.ecs_backend.list_task_definition_families(family_prefix)

        return ActionResult({"families": list(results)})

    def list_tags_for_resource(self) -> ActionResult:
        resource_arn = self._get_param("resourceArn")
        tags = self.ecs_backend.list_tags_for_resource(resource_arn)
        return ActionResult({"tags": tags})

    def tag_resource(self) -> ActionResult:
        resource_arn = self._get_param("resourceArn")
        tags = self._get_param("tags")
        self.ecs_backend.tag_resource(resource_arn, tags)
        return EmptyResult()

    def untag_resource(self) -> ActionResult:
        resource_arn = self._get_param("resourceArn")
        tag_keys = self._get_param("tagKeys")
        self.ecs_backend.untag_resource(resource_arn, tag_keys)
        return EmptyResult()

    def create_task_set(self) -> ActionResult:
        service_str = self._get_param("service")
        cluster_str = self._get_param("cluster", "default")
        task_definition = self._get_param("taskDefinition")
        external_id = self._get_param("externalId")
        network_configuration = self._get_param("networkConfiguration")
        load_balancers = self._get_param("loadBalancers")
        service_registries = self._get_param("serviceRegistries")
        launch_type = self._get_param("launchType")
        capacity_provider_strategy = self._get_param("capacityProviderStrategy")
        platform_version = self._get_param("platformVersion")
        scale = self._get_param("scale")
        client_token = self._get_param("clientToken")
        tags = self._get_param("tags")
        task_set = self.ecs_backend.create_task_set(
            service_str,
            cluster_str,
            task_definition,
            external_id=external_id,
            network_configuration=network_configuration,
            load_balancers=load_balancers,
            service_registries=service_registries,
            launch_type=launch_type,
            capacity_provider_strategy=capacity_provider_strategy,
            platform_version=platform_version,
            scale=scale,
            client_token=client_token,
            tags=tags,
        )
        return ActionResult({"taskSet": task_set})

    def describe_task_sets(self) -> ActionResult:
        cluster_str = self._get_param("cluster", "default")
        service_str = self._get_param("service")
        task_set_ids = self._get_param("taskSets", [])
        task_sets = self.ecs_backend.describe_task_sets(
            cluster_str, service_str, task_set_ids
        )
        include = self._get_param("include", [])
        if "TAGS" in include:
            self._include_in_response(DESCRIBE_TASK_SETS_TAGS)
        else:
            self._exclude_from_response(DESCRIBE_TASK_SETS_TAGS)
        return ActionResult({"taskSets": task_sets})

    def delete_task_set(self) -> ActionResult:
        cluster_str = self._get_param("cluster")
        service_str = self._get_param("service")
        task_set = self._get_param("taskSet")
        task_set = self.ecs_backend.delete_task_set(cluster_str, service_str, task_set)
        return ActionResult({"taskSet": task_set})

    def update_task_set(self) -> ActionResult:
        cluster_str = self._get_param("cluster", "default")
        service_str = self._get_param("service")
        task_set = self._get_param("taskSet")
        scale = self._get_param("scale")

        task_set = self.ecs_backend.update_task_set(
            cluster_str, service_str, task_set, scale
        )
        return ActionResult({"taskSet": task_set})

    def update_service_primary_task_set(self) -> ActionResult:
        cluster_str = self._get_param("cluster", "default")
        service_str = self._get_param("service")
        primary_task_set = self._get_param("primaryTaskSet")

        task_set = self.ecs_backend.update_service_primary_task_set(
            cluster_str, service_str, primary_task_set
        )
        return ActionResult({"taskSet": task_set})

    def put_account_setting(self) -> ActionResult:
        name = self._get_param("name")
        value = self._get_param("value")
        account_setting = self.ecs_backend.put_account_setting(name, value)
        return ActionResult({"setting": account_setting})

    def list_account_settings(self) -> ActionResult:
        name = self._get_param("name")
        value = self._get_param("value")
        account_settings = self.ecs_backend.list_account_settings(name, value)
        return ActionResult({"settings": account_settings})

    def delete_account_setting(self) -> ActionResult:
        name = self._get_param("name")
        self.ecs_backend.delete_account_setting(name)
        return EmptyResult()

    def delete_task_definitions(self) -> ActionResult:
        task_definitions = self._get_param("taskDefinitions")
        deleted_task_definitions, failures = self.ecs_backend.delete_task_definitions(
            task_definitions=task_definitions,
        )
        return ActionResult(
            {
                "taskDefinitions": deleted_task_definitions,
                "failures": failures,
            }
        )
