import re
from copy import copy
from datetime import datetime, timezone
from os import getenv
from typing import Any, Dict, Iterator, List, Optional, Tuple

from moto import settings
from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import BaseModel, CloudFormationModel
from moto.core.exceptions import JsonRESTError
from moto.core.utils import pascal_to_camelcase, remap_nested_keys, unix_time
from moto.ec2 import ec2_backends
from moto.moto_api._internal import mock_random
from moto.moto_api._internal.managed_state_model import ManagedState
from moto.utilities.utils import ARN_PARTITION_REGEX, get_partition

from ..ec2.utils import random_private_ip
from .exceptions import (
    ClusterNotFoundException,
    EcsClientException,
    InvalidParameterException,
    RevisionNotFoundException,
    ServiceNotFoundException,
    TaskDefinitionMemoryError,
    TaskDefinitionMissingPropertyError,
    TaskDefinitionNotFoundException,
    TaskSetNotFoundException,
    UnknownAccountSettingException,
)


class BaseObject(BaseModel):
    def camelCase(self, key: str) -> str:
        words = []
        for i, word in enumerate(key.split("_")):
            if i > 0:
                words.append(word.title())
            else:
                words.append(word)
        return "".join(words)

    def gen_response_object(self) -> Dict[str, Any]:
        response_object = copy(self.__dict__)
        for key, value in self.__dict__.items():
            if key.startswith("_"):
                del response_object[key]
            elif "_" in key:
                response_object[self.camelCase(key)] = value
                del response_object[key]
        return response_object

    @property
    def response_object(self) -> Dict[str, Any]:  # type: ignore[misc]
        return self.gen_response_object()


class AccountSetting(BaseObject):
    def __init__(self, name: str, value: str):
        self.name = name
        self.value = value


class Cluster(BaseObject, CloudFormationModel):
    def __init__(
        self,
        cluster_name: str,
        account_id: str,
        region_name: str,
        cluster_settings: Optional[List[Dict[str, str]]] = None,
        configuration: Optional[Dict[str, Any]] = None,
        capacity_providers: Optional[List[str]] = None,
        default_capacity_provider_strategy: Optional[List[Dict[str, Any]]] = None,
        tags: Optional[List[Dict[str, str]]] = None,
        service_connect_defaults: Optional[Dict[str, str]] = None,
    ):
        self.active_services_count = 0
        self.arn = f"arn:{get_partition(region_name)}:ecs:{region_name}:{account_id}:cluster/{cluster_name}"
        self.name = cluster_name
        self.pending_tasks_count = 0
        self.registered_container_instances_count = 0
        self.running_tasks_count = 0
        self.status = "ACTIVE"
        self.region_name = region_name
        self.settings = cluster_settings or [
            {"name": "containerInsights", "value": "disabled"}
        ]
        self.configuration = configuration
        self.capacity_providers = capacity_providers
        self.default_capacity_provider_strategy = default_capacity_provider_strategy
        self.tags = tags
        self.service_connect_defaults = service_connect_defaults

    @property
    def physical_resource_id(self) -> str:
        return self.name

    @property
    def response_object(self) -> Dict[str, Any]:  # type: ignore[misc]
        response_object = self.gen_response_object()
        response_object["clusterArn"] = self.arn
        response_object["clusterName"] = self.name
        response_object["capacityProviders"] = self.capacity_providers
        response_object["defaultCapacityProviderStrategy"] = (
            self.default_capacity_provider_strategy
        )
        del response_object["arn"], response_object["name"]
        return response_object

    @staticmethod
    def cloudformation_name_type() -> str:
        return "ClusterName"

    @staticmethod
    def cloudformation_type() -> str:
        # https://docs.aws.amazon.com/AWSCloudFormation/latest/UserGuide/aws-resource-ecs-cluster.html
        return "AWS::ECS::Cluster"

    @classmethod
    def create_from_cloudformation_json(  # type: ignore[misc]
        cls,
        resource_name: str,
        cloudformation_json: Any,
        account_id: str,
        region_name: str,
        **kwargs: Any,
    ) -> "Cluster":
        ecs_backend = ecs_backends[account_id][region_name]
        return ecs_backend.create_cluster(
            # ClusterName is optional in CloudFormation, thus create a random
            # name if necessary
            cluster_name=resource_name
        )

    @classmethod
    def update_from_cloudformation_json(  # type: ignore[misc]
        cls,
        original_resource: Any,
        new_resource_name: str,
        cloudformation_json: Any,
        account_id: str,
        region_name: str,
    ) -> "Cluster":
        if original_resource.name != new_resource_name:
            ecs_backend = ecs_backends[account_id][region_name]
            ecs_backend.delete_cluster(original_resource.arn)
            return ecs_backend.create_cluster(
                # ClusterName is optional in CloudFormation, thus create a
                # random name if necessary
                cluster_name=new_resource_name
            )
        else:
            # no-op when nothing changed between old and new resources
            return original_resource

    @classmethod
    def has_cfn_attr(cls, attr: str) -> bool:
        return attr in ["Arn"]

    def get_cfn_attribute(self, attribute_name: str) -> str:
        from moto.cloudformation.exceptions import UnformattedGetAttTemplateException

        if attribute_name == "Arn":
            return self.arn
        raise UnformattedGetAttTemplateException()


class TaskDefinition(BaseObject, CloudFormationModel):
    def __init__(
        self,
        family: str,
        revision: int,
        container_definitions: List[Dict[str, Any]],
        account_id: str,
        region_name: str,
        network_mode: Optional[str] = None,
        volumes: Optional[List[Dict[str, Any]]] = None,
        tags: Optional[List[Dict[str, str]]] = None,
        placement_constraints: Optional[List[Dict[str, str]]] = None,
        requires_compatibilities: Optional[List[str]] = None,
        cpu: Optional[str] = None,
        memory: Optional[str] = None,
        task_role_arn: Optional[str] = None,
        execution_role_arn: Optional[str] = None,
        proxy_configuration: Optional[Dict[str, Any]] = None,
        inference_accelerators: Optional[List[Dict[str, str]]] = None,
        runtime_platform: Optional[Dict[str, str]] = None,
        ipc_mode: Optional[str] = None,
        pid_mode: Optional[str] = None,
        ephemeral_storage: Optional[Dict[str, int]] = None,
    ):
        self.family = family
        self.revision = revision
        self.arn = f"arn:{get_partition(region_name)}:ecs:{region_name}:{account_id}:task-definition/{family}:{revision}"

        default_container_definition = {
            "cpu": 0,
            "portMappings": [],
            "essential": True,
            "environment": [],
            "mountPoints": [],
            "volumesFrom": [],
        }
        self.container_definitions = []
        for container_definition in container_definitions:
            full_definition = default_container_definition.copy()
            full_definition.update(container_definition)
            self.container_definitions.append(full_definition)

        self.tags = tags if tags is not None else []

        if volumes is None:
            self.volumes = []
        else:
            self.volumes = volumes
            for volume in volumes:
                if "efsVolumeConfiguration" in volume:
                    # We should reach into EFS to verify this volume exists
                    efs_config = volume["efsVolumeConfiguration"]
                    if "rootDirectory" not in efs_config:
                        efs_config["rootDirectory"] = "/"

        if not requires_compatibilities or requires_compatibilities == ["EC2"]:
            self.compatibilities = ["EC2"]
        else:
            self.compatibilities = ["EC2", "FARGATE"]

        if network_mode is None and "FARGATE" not in self.compatibilities:
            self.network_mode: Optional[str] = "bridge"
        elif "FARGATE" in self.compatibilities:
            self.network_mode: Optional[str] = "awsvpc"  # type: ignore[no-redef]
        else:
            self.network_mode = network_mode

        if task_role_arn is not None:
            self.task_role_arn = task_role_arn
        if execution_role_arn is not None:
            self.execution_role_arn = execution_role_arn

        self.placement_constraints = (
            placement_constraints if placement_constraints is not None else []
        )

        self.requires_compatibilities = requires_compatibilities
        self.proxy_configuration = proxy_configuration
        self.inference_accelerators = inference_accelerators
        self.runtime_platform = runtime_platform
        self.ipc_mode = ipc_mode
        self.pid_mode = pid_mode
        self.ephemeral_storage = ephemeral_storage

        self.cpu = cpu
        self.memory = memory
        self.status = "ACTIVE"

    @property
    def response_object(self) -> Dict[str, Any]:  # type: ignore[misc]
        response_object = self.gen_response_object()
        response_object["taskDefinitionArn"] = response_object["arn"]
        del response_object["arn"]
        del response_object["tags"]

        if not response_object["requiresCompatibilities"]:
            del response_object["requiresCompatibilities"]
        if not response_object["cpu"]:
            del response_object["cpu"]
        if not response_object["memory"]:
            del response_object["memory"]

        return response_object

    @property
    def physical_resource_id(self) -> str:
        return self.arn

    @staticmethod
    def cloudformation_name_type() -> str:
        return ""

    @staticmethod
    def cloudformation_type() -> str:
        # https://docs.aws.amazon.com/AWSCloudFormation/latest/UserGuide/aws-resource-ecs-taskdefinition.html
        return "AWS::ECS::TaskDefinition"

    @classmethod
    def create_from_cloudformation_json(  # type: ignore[misc]
        cls,
        resource_name: str,
        cloudformation_json: Any,
        account_id: str,
        region_name: str,
        **kwargs: Any,
    ) -> "TaskDefinition":
        properties = cloudformation_json["Properties"]

        family = properties.get(
            "Family", f"task-definition-{int(mock_random.random() * 10**6)}"
        )
        container_definitions = remap_nested_keys(
            properties.get("ContainerDefinitions", []), pascal_to_camelcase
        )
        volumes = remap_nested_keys(properties.get("Volumes", []), pascal_to_camelcase)

        ecs_backend = ecs_backends[account_id][region_name]
        return ecs_backend.register_task_definition(
            family=family, container_definitions=container_definitions, volumes=volumes
        )

    @classmethod
    def update_from_cloudformation_json(  # type: ignore[misc]
        cls,
        original_resource: Any,
        new_resource_name: str,
        cloudformation_json: Any,
        account_id: str,
        region_name: str,
    ) -> "TaskDefinition":
        properties = cloudformation_json["Properties"]
        family = properties.get(
            "Family", f"task-definition-{int(mock_random.random() * 10**6)}"
        )
        container_definitions = properties["ContainerDefinitions"]
        volumes = properties.get("Volumes")
        if (
            original_resource.family != family
            or original_resource.container_definitions != container_definitions
            or original_resource.volumes != volumes
        ):
            # currently TaskRoleArn isn't stored at TaskDefinition
            # instances
            ecs_backend = ecs_backends[account_id][region_name]
            ecs_backend.deregister_task_definition(original_resource.arn)
            return ecs_backend.register_task_definition(
                family=family,
                container_definitions=container_definitions,
                volumes=volumes,
            )
        else:
            # no-op when nothing changed between old and new resources
            return original_resource


class Task(BaseObject, ManagedState):
    def __init__(
        self,
        cluster: Cluster,
        task_definition: TaskDefinition,
        container_instance_arn: Optional[str],
        resource_requirements: Optional[Dict[str, str]],
        backend: "EC2ContainerServiceBackend",
        launch_type: str = "",
        overrides: Optional[Dict[str, Any]] = None,
        started_by: str = "",
        tags: Optional[List[Dict[str, str]]] = None,
        networking_configuration: Optional[Dict[str, Any]] = None,
    ):
        # Configure ManagedState
        # https://docs.aws.amazon.com/AmazonECS/latest/developerguide/task-lifecycle.html
        super().__init__(
            model_name="ecs::task",
            transitions=[
                # We start in RUNNING state in order not to break existing tests.
                # ("PROVISIONING", "PENDING"),
                # ("PENDING", "ACTIVATING"),
                # ("ACTIVATING", "RUNNING"),
                ("RUNNING", "DEACTIVATING"),
                ("DEACTIVATING", "STOPPING"),
                ("STOPPING", "DEPROVISIONING"),
                ("DEPROVISIONING", "STOPPED"),
                # There seems to be race condition, where the waiter expects the task to be in
                # STOPPED state, but it is already in DELETED state.
                # ("STOPPED", "DELETED"),
            ],
        )
        self.id = str(mock_random.uuid4())
        self.cluster_name = cluster.name
        self.cluster_arn = cluster.arn
        self.container_instance_arn = container_instance_arn
        self.desired_status = "RUNNING"
        self.task_definition_arn = task_definition.arn
        self.overrides = overrides or {}
        self.containers: List[Dict[str, Any]] = []
        self.started_by = started_by
        self.tags = tags or []
        self.launch_type = launch_type
        self.stopped_reason = ""
        self.resource_requirements = resource_requirements
        self.region_name = cluster.region_name
        self._account_id = backend.account_id
        self._backend = backend
        self.attachments = []

        if task_definition.network_mode == "awsvpc":
            if not networking_configuration:
                raise InvalidParameterException(
                    "Network Configuration must be provided when networkMode 'awsvpc' is specified."
                )

            self.network_configuration = networking_configuration
            net_conf = networking_configuration["awsvpcConfiguration"]
            ec2_backend = ec2_backends[self._account_id][self.region_name]

            eni = ec2_backend.create_network_interface(
                subnet=net_conf["subnets"][0],
                private_ip_address=random_private_ip(),
                group_ids=net_conf["securityGroups"],
                description="moto ECS",
            )
            eni.status = "in-use"
            eni.device_index = 0

            self.attachments.append(
                {
                    "id": str(mock_random.uuid4()),
                    "type": "ElasticNetworkInterface",
                    "status": "ATTACHED",
                    "details": [
                        {"name": "subnetId", "value": net_conf["subnets"][0]},
                        {"name": "networkInterfaceId", "value": eni.id},
                        {"name": "macAddress", "value": eni.mac_address},
                        {"name": "privateDnsName", "value": eni.private_dns_name},
                        {"name": "privateIPv4Address", "value": eni.private_ip_address},
                    ],
                }
            )

    @property
    def last_status(self) -> Optional[str]:
        return self.status  # managed state

    @last_status.setter
    def last_status(self, value: Optional[str]) -> None:
        self.status = value

    @property
    def task_arn(self) -> str:
        if self._backend.enable_long_arn_for_name(name="taskLongArnFormat"):
            return f"arn:{get_partition(self.region_name)}:ecs:{self.region_name}:{self._account_id}:task/{self.cluster_name}/{self.id}"
        return f"arn:{get_partition(self.region_name)}:ecs:{self.region_name}:{self._account_id}:task/{self.id}"

    def response_object(self, include_tags: bool = True) -> Dict[str, Any]:  # type: ignore
        response_object = self.gen_response_object()
        if not include_tags:
            response_object.pop("tags", None)
        response_object["taskArn"] = self.task_arn
        response_object["lastStatus"] = self.last_status
        return response_object


class CapacityProvider(BaseObject):
    def __init__(
        self,
        account_id: str,
        region_name: str,
        name: str,
        asg_details: Dict[str, Any],
        tags: Optional[List[Dict[str, str]]],
    ):
        self._id = str(mock_random.uuid4())
        self.capacity_provider_arn = f"arn:{get_partition(region_name)}:ecs:{region_name}:{account_id}:capacity-provider/{name}"
        self.name = name
        self.status = "ACTIVE"
        self.auto_scaling_group_provider = self._prepare_asg_provider(asg_details)
        self.tags = tags

        self.update_status: Optional[str] = None

    def _prepare_asg_provider(self, asg_details: Dict[str, Any]) -> Dict[str, Any]:
        if "managedScaling" not in asg_details:
            asg_details["managedScaling"] = {}
        if asg_details["managedScaling"].get("instanceWarmupPeriod") is None:
            asg_details["managedScaling"]["instanceWarmupPeriod"] = 300
        if not asg_details["managedScaling"].get("minimumScalingStepSize"):
            asg_details["managedScaling"]["minimumScalingStepSize"] = 1
        if not asg_details["managedScaling"].get("maximumScalingStepSize"):
            asg_details["managedScaling"]["maximumScalingStepSize"] = 10000
        if not asg_details["managedScaling"].get("targetCapacity"):
            asg_details["managedScaling"]["targetCapacity"] = 100
        if not asg_details["managedScaling"].get("status"):
            asg_details["managedScaling"]["status"] = "DISABLED"
        if "managedTerminationProtection" not in asg_details:
            asg_details["managedTerminationProtection"] = "DISABLED"
        return asg_details

    def update(self, asg_details: Dict[str, Any]) -> None:
        if "managedTerminationProtection" in asg_details:
            self.auto_scaling_group_provider["managedTerminationProtection"] = (
                asg_details["managedTerminationProtection"]
            )
        if "managedScaling" in asg_details:
            scaling_props = [
                "status",
                "targetCapacity",
                "minimumScalingStepSize",
                "maximumScalingStepSize",
                "instanceWarmupPeriod",
            ]
            for prop in scaling_props:
                if prop in asg_details["managedScaling"]:
                    self.auto_scaling_group_provider["managedScaling"][prop] = (
                        asg_details["managedScaling"][prop]
                    )
        self.auto_scaling_group_provider = self._prepare_asg_provider(
            self.auto_scaling_group_provider
        )
        self.update_status = "UPDATE_COMPLETE"


class CapacityProviderFailure(BaseObject):
    def __init__(self, reason: str, name: str, account_id: str, region_name: str):
        self.reason = reason
        self.arn = f"arn:{get_partition(region_name)}:ecs:{region_name}:{account_id}:capacity_provider/{name}"

    @property
    def response_object(self) -> Dict[str, Any]:  # type: ignore[misc]
        response_object = self.gen_response_object()
        response_object["reason"] = self.reason
        response_object["arn"] = self.arn
        return response_object


class Service(BaseObject, CloudFormationModel):
    """Set the environment variable MOTO_ECS_SERVICE_RUNNING to a number of running tasks you want
    the service to transition to, ie if set to 2:

    MOTO_ECS_SERVICE_RUNNING=2

    then describe_services call to return runningCount of the service AND deployment to 2
    """

    def __init__(
        self,
        cluster: Cluster,
        service_name: str,
        desired_count: int,
        backend: "EC2ContainerServiceBackend",
        task_definition: Optional[TaskDefinition] = None,
        load_balancers: Optional[List[Dict[str, Any]]] = None,
        scheduling_strategy: Optional[List[Dict[str, Any]]] = None,
        tags: Optional[List[Dict[str, str]]] = None,
        deployment_controller: Optional[Dict[str, str]] = None,
        launch_type: Optional[str] = None,
        service_registries: Optional[List[Dict[str, Any]]] = None,
        platform_version: Optional[str] = None,
    ):
        self.cluster_name = cluster.name
        self.cluster_arn = cluster.arn
        self.name = service_name
        self.status = "ACTIVE"
        self.task_definition = task_definition.arn if task_definition else None
        self.desired_count = desired_count
        self.task_sets: List[TaskSet] = []
        self.deployment_controller = deployment_controller or {"type": "ECS"}
        self.events: List[Dict[str, Any]] = []
        self.launch_type = launch_type
        self.service_registries = service_registries or []
        self.load_balancers = load_balancers if load_balancers is not None else []
        self.scheduling_strategy = (
            scheduling_strategy if scheduling_strategy is not None else "REPLICA"
        )
        self.platform_version = platform_version
        self.tags = tags if tags is not None else []
        self.region_name = cluster.region_name
        self._account_id = backend.account_id
        self._backend = backend

        try:
            # negative running count not allowed, set to 0 if so
            ecs_running_count = max(int(getenv("MOTO_ECS_SERVICE_RUNNING", 0)), 0)
        except ValueError:
            # Unable to parse value of MOTO_ECS_SERVICE_RUNNING as an integer, set to default 0
            ecs_running_count = 0

        self.running_count = ecs_running_count
        self.pending_count = desired_count - ecs_running_count
        if self.deployment_controller["type"] == "ECS":
            self.deployments = [
                {
                    "createdAt": datetime.now(timezone.utc),
                    "desiredCount": self.desired_count,
                    "id": f"ecs-svc/{mock_random.randint(0, 32**12)}",
                    "launchType": self.launch_type,
                    "pendingCount": self.pending_count,
                    "runningCount": ecs_running_count,
                    "status": "PRIMARY",
                    "taskDefinition": self.task_definition,
                    "updatedAt": datetime.now(timezone.utc),
                }
            ]
        else:
            self.deployments = []

    @property
    def arn(self) -> str:
        if self._backend.enable_long_arn_for_name(name="serviceLongArnFormat"):
            return f"arn:{get_partition(self.region_name)}:ecs:{self.region_name}:{self._account_id}:service/{self.cluster_name}/{self.name}"
        return f"arn:{get_partition(self.region_name)}:ecs:{self.region_name}:{self._account_id}:service/{self.name}"

    @property
    def physical_resource_id(self) -> str:
        return self.arn

    @property
    def response_object(self) -> Dict[str, Any]:  # type: ignore[misc]
        response_object = self.gen_response_object()
        del response_object["name"], response_object["tags"]
        response_object["serviceName"] = self.name
        response_object["serviceArn"] = self.arn
        response_object["schedulingStrategy"] = self.scheduling_strategy
        response_object["platformVersion"] = self.platform_version
        if response_object["deploymentController"]["type"] == "ECS":
            del response_object["deploymentController"]
            del response_object["taskSets"]
        else:
            response_object["taskSets"] = [
                t.response_object for t in response_object["taskSets"]
            ]

        for deployment in response_object["deployments"]:
            if isinstance(deployment["createdAt"], datetime):
                deployment["createdAt"] = unix_time(
                    deployment["createdAt"].replace(tzinfo=None)
                )
            if isinstance(deployment["updatedAt"], datetime):
                deployment["updatedAt"] = unix_time(
                    deployment["updatedAt"].replace(tzinfo=None)
                )

        return response_object

    @staticmethod
    def cloudformation_name_type() -> str:
        return "ServiceName"

    @staticmethod
    def cloudformation_type() -> str:
        # https://docs.aws.amazon.com/AWSCloudFormation/latest/UserGuide/aws-resource-ecs-service.html
        return "AWS::ECS::Service"

    @classmethod
    def create_from_cloudformation_json(  # type: ignore[misc]
        cls,
        resource_name: str,
        cloudformation_json: Any,
        account_id: str,
        region_name: str,
        **kwargs: Any,
    ) -> "Service":
        properties = cloudformation_json["Properties"]
        if isinstance(properties["Cluster"], Cluster):
            cluster = properties["Cluster"].name
        else:
            cluster = properties["Cluster"]
        if isinstance(properties["TaskDefinition"], TaskDefinition):
            task_definition = properties["TaskDefinition"].family
        else:
            task_definition = properties["TaskDefinition"]
        desired_count = properties.get("DesiredCount", None)
        # TODO: LoadBalancers
        # TODO: Role

        ecs_backend = ecs_backends[account_id][region_name]
        return ecs_backend.create_service(
            cluster, resource_name, desired_count, task_definition_str=task_definition
        )

    @classmethod
    def update_from_cloudformation_json(  # type: ignore[misc]
        cls,
        original_resource: Any,
        new_resource_name: str,
        cloudformation_json: Any,
        account_id: str,
        region_name: str,
    ) -> "Service":
        properties = cloudformation_json["Properties"]
        if isinstance(properties["Cluster"], Cluster):
            cluster_name = properties["Cluster"].name
        else:
            cluster_name = properties["Cluster"]
        if isinstance(properties["TaskDefinition"], TaskDefinition):
            task_definition = properties["TaskDefinition"].family
        else:
            task_definition = properties["TaskDefinition"]
        desired_count = properties.get("DesiredCount", None)

        ecs_backend = ecs_backends[account_id][region_name]
        service_name = original_resource.name
        if (
            original_resource.cluster_arn
            != Cluster(cluster_name, account_id, region_name).arn
        ):
            # TODO: LoadBalancers
            # TODO: Role
            ecs_backend.delete_service(
                original_resource.cluster_name, service_name, force=True
            )
            return ecs_backend.create_service(
                cluster_name,
                new_resource_name,
                desired_count,
                task_definition_str=task_definition,
            )
        else:
            return ecs_backend.update_service(
                {
                    "cluster": cluster_name,
                    "service": service_name,
                    "task_definition": task_definition,
                    "desired_count": desired_count,
                }
            )

    @classmethod
    def has_cfn_attr(cls, attr: str) -> bool:
        return attr in ["Name"]

    def get_cfn_attribute(self, attribute_name: str) -> str:
        from moto.cloudformation.exceptions import UnformattedGetAttTemplateException

        if attribute_name == "Name":
            return self.name
        raise UnformattedGetAttTemplateException()


class ContainerInstance(BaseObject):
    def __init__(
        self,
        ec2_instance_id: str,
        account_id: str,
        region_name: str,
        cluster_name: str,
        backend: "EC2ContainerServiceBackend",
    ):
        self.ec2_instance_id = ec2_instance_id
        self.agent_connected = True
        self.status = "ACTIVE"
        self.registered_resources: List[Dict[str, Any]] = [
            {
                "doubleValue": 0.0,
                "integerValue": 4096,
                "longValue": 0,
                "name": "CPU",
                "type": "INTEGER",
            },
            {
                "doubleValue": 0.0,
                "integerValue": 7482,
                "longValue": 0,
                "name": "MEMORY",
                "type": "INTEGER",
            },
            {
                "doubleValue": 0.0,
                "integerValue": 0,
                "longValue": 0,
                "name": "PORTS",
                "stringSetValue": ["22", "2376", "2375", "51678", "51679"],
                "type": "STRINGSET",
            },
            {
                "doubleValue": 0.0,
                "integerValue": 0,
                "longValue": 0,
                "name": "PORTS_UDP",
                "stringSetValue": [],
                "type": "STRINGSET",
            },
        ]
        self.pending_tasks_count = 0
        self.remaining_resources: List[Dict[str, Any]] = [
            {
                "doubleValue": 0.0,
                "integerValue": 4096,
                "longValue": 0,
                "name": "CPU",
                "type": "INTEGER",
            },
            {
                "doubleValue": 0.0,
                "integerValue": 7482,
                "longValue": 0,
                "name": "MEMORY",
                "type": "INTEGER",
            },
            {
                "doubleValue": 0.0,
                "integerValue": 0,
                "longValue": 0,
                "name": "PORTS",
                "stringSetValue": ["22", "2376", "2375", "51678", "51679"],
                "type": "STRINGSET",
            },
            {
                "doubleValue": 0.0,
                "integerValue": 0,
                "longValue": 0,
                "name": "PORTS_UDP",
                "stringSetValue": [],
                "type": "STRINGSET",
            },
        ]
        self.running_tasks_count = 0
        self.version_info = {
            "agentVersion": "1.0.0",
            "agentHash": "4023248",
            "dockerVersion": "DockerVersion: 1.5.0",
        }
        ec2_backend = ec2_backends[account_id][region_name]
        ec2_instance = ec2_backend.get_instance(ec2_instance_id)
        self.attributes = {
            "ecs.ami-id": ec2_instance.image_id,
            "ecs.availability-zone": ec2_instance.placement,
            "ecs.instance-type": ec2_instance.instance_type,
            "ecs.os-type": ec2_instance.platform
            if ec2_instance.platform == "windows"
            else "linux",  # options are windows and linux, linux is default
        }
        self.registered_at = datetime.now(timezone.utc)
        self.region_name = region_name
        self.id = str(mock_random.uuid4())
        self.cluster_name = cluster_name
        self._account_id = backend.account_id
        self._backend = backend

    @property
    def container_instance_arn(self) -> str:
        if self._backend.enable_long_arn_for_name(
            name="containerInstanceLongArnFormat"
        ):
            return f"arn:{get_partition(self.region_name)}:ecs:{self.region_name}:{self._account_id}:container-instance/{self.cluster_name}/{self.id}"
        return f"arn:{get_partition(self.region_name)}:ecs:{self.region_name}:{self._account_id}:container-instance/{self.id}"

    @property
    def response_object(self) -> Dict[str, Any]:  # type: ignore[misc]
        response_object = self.gen_response_object()
        response_object["containerInstanceArn"] = self.container_instance_arn
        response_object["attributes"] = [
            self._format_attribute(name, value)
            for name, value in response_object["attributes"].items()
        ]
        if isinstance(response_object["registeredAt"], datetime):
            response_object["registeredAt"] = unix_time(
                response_object["registeredAt"].replace(tzinfo=None)
            )
        return response_object

    def _format_attribute(self, name: str, value: Optional[str]) -> Dict[str, str]:
        formatted_attr = {"name": name}
        if value is not None:
            formatted_attr["value"] = value
        return formatted_attr


class ClusterFailure(BaseObject):
    def __init__(
        self, reason: str, cluster_name: str, account_id: str, region_name: str
    ):
        self.reason = reason
        self.arn = f"arn:{get_partition(region_name)}:ecs:{region_name}:{account_id}:cluster/{cluster_name}"

    @property
    def response_object(self) -> Dict[str, Any]:  # type: ignore[misc]
        response_object = self.gen_response_object()
        response_object["reason"] = self.reason
        response_object["arn"] = self.arn
        return response_object


class ContainerInstanceFailure(BaseObject):
    def __init__(
        self, reason: str, container_instance_id: str, account_id: str, region_name: str
    ):
        self.reason = reason
        self.arn = f"arn:{get_partition(region_name)}:ecs:{region_name}:{account_id}:container-instance/{container_instance_id}"

    @property
    def response_object(self) -> Dict[str, Any]:  # type: ignore[misc]
        response_object = self.gen_response_object()
        response_object["reason"] = self.reason
        response_object["arn"] = self.arn
        return response_object


class TaskSet(BaseObject):
    def __init__(
        self,
        service: str,
        cluster: str,
        task_definition: str,
        account_id: str,
        region_name: str,
        external_id: Optional[str] = None,
        network_configuration: Optional[Dict[str, Any]] = None,
        load_balancers: Optional[List[Dict[str, Any]]] = None,
        service_registries: Optional[List[Dict[str, Any]]] = None,
        launch_type: Optional[str] = None,
        capacity_provider_strategy: Optional[List[Dict[str, Any]]] = None,
        platform_version: Optional[str] = None,
        scale: Optional[Dict[str, Any]] = None,
        client_token: Optional[str] = None,
        tags: Optional[List[Dict[str, str]]] = None,
    ):
        self.service = service
        self.cluster = cluster
        self.status = "ACTIVE"
        self.task_definition = task_definition or ""
        self.region_name = region_name
        self.external_id = external_id or ""
        self.network_configuration = network_configuration or None
        self.load_balancers = load_balancers or []
        self.service_registries = service_registries or []
        self.launch_type = launch_type
        self.capacity_provider_strategy = capacity_provider_strategy or []
        self.platform_version = platform_version or "LATEST"
        self.scale = scale or {"value": 100.0, "unit": "PERCENT"}
        self.client_token = client_token or ""
        self.tags = tags or []
        self.stabilityStatus = "STEADY_STATE"
        self.createdAt = datetime.now(timezone.utc)
        self.updatedAt = datetime.now(timezone.utc)
        self.stabilityStatusAt = datetime.now(timezone.utc)
        self.id = f"ecs-svc/{mock_random.randint(0, 32**12)}"
        self.service_arn = ""
        self.cluster_arn = ""

        cluster_name = self.cluster.split("/")[-1]
        service_name = self.service.split("/")[-1]
        self.task_set_arn = f"arn:{get_partition(region_name)}:ecs:{region_name}:{account_id}:task-set/{cluster_name}/{service_name}/{self.id}"

    @property
    def response_object(self) -> Dict[str, Any]:  # type: ignore[misc]
        response_object = self.gen_response_object()
        if isinstance(response_object["createdAt"], datetime):
            response_object["createdAt"] = unix_time(
                self.createdAt.replace(tzinfo=None)
            )
        if isinstance(response_object["updatedAt"], datetime):
            response_object["updatedAt"] = unix_time(
                self.updatedAt.replace(tzinfo=None)
            )
        if isinstance(response_object["stabilityStatusAt"], datetime):
            response_object["stabilityStatusAt"] = unix_time(
                self.stabilityStatusAt.replace(tzinfo=None)
            )
        del response_object["service"]
        del response_object["cluster"]
        return response_object


class EC2ContainerServiceBackend(BaseBackend):
    """
    ECS resources use the new ARN format by default.
    Use the following environment variable to revert back to the old/short ARN format:
    `MOTO_ECS_NEW_ARN=false`

    AWS reference: https://aws.amazon.com/blogs/compute/migrating-your-amazon-ecs-deployment-to-the-new-arn-and-resource-id-format-2/
    """

    def __init__(self, region_name: str, account_id: str):
        super().__init__(region_name, account_id)
        self.account_settings: Dict[str, AccountSetting] = dict()
        self.capacity_providers: Dict[str, CapacityProvider] = dict()
        self.clusters: Dict[str, Cluster] = {}
        self.task_definitions: Dict[str, Dict[int, TaskDefinition]] = {}
        self.tasks: Dict[str, Dict[str, Task]] = {}
        self.services: Dict[str, Service] = {}
        self.container_instances: Dict[str, Dict[str, ContainerInstance]] = {}

    def _get_cluster(self, name: str) -> Cluster:
        # short name or full ARN of the cluster
        cluster_name = name.split("/")[-1]

        cluster = self.clusters.get(cluster_name)
        if not cluster:
            raise ClusterNotFoundException

        return cluster

    def create_capacity_provider(
        self,
        name: str,
        asg_details: Dict[str, Any],
        tags: Optional[List[Dict[str, str]]],
    ) -> CapacityProvider:
        capacity_provider = CapacityProvider(
            self.account_id, self.region_name, name, asg_details, tags
        )
        self.capacity_providers[name] = capacity_provider
        return capacity_provider

    def describe_task_definition(self, task_definition_str: str) -> TaskDefinition:
        task_definition_name = task_definition_str.split("/")[-1]
        if ":" in task_definition_name:
            family, rev = task_definition_name.split(":")
            revision = int(rev)
        else:
            family = task_definition_name
            revision = self._get_last_task_definition_revision_id(family)

        if (
            family in self.task_definitions
            and revision in self.task_definitions[family]
        ):
            return self.task_definitions[family][revision]
        else:
            raise TaskDefinitionNotFoundException()

    def create_cluster(
        self,
        cluster_name: str,
        tags: Any = None,
        cluster_settings: Any = None,
        configuration: Optional[Dict[str, Any]] = None,
        capacity_providers: Optional[List[str]] = None,
        default_capacity_provider_strategy: Optional[List[Dict[str, Any]]] = None,
        service_connect_defaults: Optional[Dict[str, str]] = None,
    ) -> Cluster:
        cluster = Cluster(
            cluster_name,
            self.account_id,
            self.region_name,
            cluster_settings,
            configuration,
            capacity_providers,
            default_capacity_provider_strategy,
            tags,
            service_connect_defaults=service_connect_defaults,
        )
        self.clusters[cluster_name] = cluster
        return cluster

    def update_cluster(
        self,
        cluster_name: str,
        cluster_settings: Optional[List[Dict[str, str]]],
        configuration: Optional[Dict[str, Any]],
        service_connect_defaults: Optional[Dict[str, str]],
    ) -> Cluster:
        """
        The serviceConnectDefaults-parameter is not yet implemented
        """
        cluster = self._get_cluster(cluster_name)
        if cluster_settings is not None:
            cluster.settings = cluster_settings
        if configuration is not None:
            cluster.configuration = configuration
        if service_connect_defaults is not None:
            cluster.service_connect_defaults = service_connect_defaults
        return cluster

    def put_cluster_capacity_providers(
        self,
        cluster_name: str,
        capacity_providers: Optional[List[str]],
        default_capacity_provider_strategy: Optional[List[Dict[str, Any]]],
    ) -> Cluster:
        cluster = self._get_cluster(cluster_name)
        if capacity_providers is not None:
            cluster.capacity_providers = capacity_providers
        if default_capacity_provider_strategy is not None:
            cluster.default_capacity_provider_strategy = (
                default_capacity_provider_strategy
            )
        return cluster

    def _get_provider(self, name_or_arn: str) -> Optional[CapacityProvider]:
        for provider in self.capacity_providers.values():
            if (
                provider.name == name_or_arn
                or provider.capacity_provider_arn == name_or_arn
            ):
                return provider
        return None

    def describe_capacity_providers(
        self, names: List[str]
    ) -> Tuple[List[CapacityProvider], List[CapacityProviderFailure]]:
        providers = []
        failures = []
        for name in names:
            provider = self._get_provider(name)
            if provider:
                providers.append(provider)
            else:
                failures.append(
                    CapacityProviderFailure(
                        "MISSING", name, self.account_id, self.region_name
                    )
                )
        return providers, failures

    def delete_capacity_provider(self, name_or_arn: str) -> CapacityProvider:
        provider: CapacityProvider = self._get_provider(name_or_arn)  # type: ignore[assignment]
        self.capacity_providers.pop(provider.name)
        return provider

    def update_capacity_provider(
        self, name_or_arn: str, asg_provider: Dict[str, Any]
    ) -> CapacityProvider:
        provider: CapacityProvider = self._get_provider(name_or_arn)  # type: ignore[assignment]
        provider.update(asg_provider)
        return provider

    def list_clusters(self) -> List[str]:
        """
        maxSize and pagination not implemented
        """
        return [cluster.arn for cluster in self.clusters.values()]

    def describe_clusters(
        self,
        list_clusters_name: Optional[List[str]] = None,
        include: Optional[List[str]] = None,
    ) -> Tuple[List[Dict[str, Any]], List[ClusterFailure]]:
        """
        Only include=TAGS is currently supported.
        """
        list_clusters = []
        failures = []
        if list_clusters_name is None:
            if "default" in self.clusters:
                list_clusters.append(self.clusters["default"].response_object)
        else:
            for cluster_name in list_clusters_name:
                cluster_name = cluster_name.split("/")[-1]
                if cluster_name in self.clusters:
                    list_clusters.append(self.clusters[cluster_name].response_object)
                else:
                    failures.append(
                        ClusterFailure(
                            "MISSING", cluster_name, self.account_id, self.region_name
                        )
                    )

        if not include or "TAGS" not in (include):
            for cluster in list_clusters:
                cluster["tags"] = None

        return list_clusters, failures

    def delete_cluster(self, cluster_str: str) -> Cluster:
        cluster = self._get_cluster(cluster_str)

        # A cluster is not immediately removed - just marked as inactive
        # It is only deleted later on
        # https://boto3.amazonaws.com/v1/documentation/api/latest/reference/services/ecs.html#ECS.Client.delete_cluster
        cluster.status = "INACTIVE"

        return cluster

    def register_task_definition(
        self,
        family: str,
        container_definitions: List[Dict[str, Any]],
        volumes: Optional[List[Dict[str, Any]]] = None,
        network_mode: Optional[str] = None,
        tags: Optional[List[Dict[str, str]]] = None,
        placement_constraints: Optional[List[Dict[str, str]]] = None,
        requires_compatibilities: Optional[List[str]] = None,
        cpu: Optional[str] = None,
        memory: Optional[str] = None,
        task_role_arn: Optional[str] = None,
        execution_role_arn: Optional[str] = None,
        proxy_configuration: Optional[Dict[str, Any]] = None,
        inference_accelerators: Optional[List[Dict[str, str]]] = None,
        runtime_platform: Optional[Dict[str, str]] = None,
        ipc_mode: Optional[str] = None,
        pid_mode: Optional[str] = None,
        ephemeral_storage: Optional[Dict[str, int]] = None,
    ) -> TaskDefinition:
        if requires_compatibilities and "FARGATE" in requires_compatibilities:
            # TODO need more validation for Fargate
            if pid_mode and pid_mode != "task":
                raise EcsClientException(
                    f"Tasks using the Fargate launch type do not support pidMode '{pid_mode}'. The supported value for pidMode is 'task'."
                )
        self._validate_container_defs(
            memory, container_definitions, requires_compatibilities
        )

        if family in self.task_definitions:
            last_id = self._get_last_task_definition_revision_id(family)
            revision = (last_id or 0) + 1
        else:
            self.task_definitions[family] = {}
            revision = 1
        task_definition = TaskDefinition(
            family,
            revision,
            container_definitions,
            self.account_id,
            self.region_name,
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
        self.task_definitions[family][revision] = task_definition

        return task_definition

    @staticmethod
    def _validate_container_defs(  # type: ignore[misc]
        memory: Optional[str],
        container_definitions: List[Dict[str, Any]],
        requires_compatibilities: Optional[List[str]],
    ) -> None:
        # The capitalised keys are passed by Cloudformation
        for cd in container_definitions:
            if "name" not in cd and "Name" not in cd:
                raise TaskDefinitionMissingPropertyError("name")
            if "image" not in cd and "Image" not in cd:
                raise TaskDefinitionMissingPropertyError("image")
            if (
                requires_compatibilities
                and "EC2" in requires_compatibilities
                and ("memory" not in cd and "Memory" not in cd and not memory)
            ):
                raise TaskDefinitionMemoryError(cd["name"])
            if (
                "memory" not in cd
                and "Memory" not in cd
                and "memoryReservation" not in cd
                and "MemoryReservation" not in cd
                and not memory
            ):
                raise TaskDefinitionMemoryError(cd["name"])

    def list_task_definitions(self, family_prefix: str) -> List[str]:
        task_arns = []
        for task_definition_list in self.task_definitions.values():
            task_arns.extend(
                [
                    task_definition.arn
                    for task_definition in task_definition_list.values()
                    if family_prefix is None or task_definition.family == family_prefix
                ]
            )
        return task_arns

    def deregister_task_definition(self, task_definition_str: str) -> TaskDefinition:
        task_definition_name = task_definition_str.split("/")[-1]
        try:
            family, rev = task_definition_name.split(":")
        except ValueError:
            raise RevisionNotFoundException
        try:
            revision = int(rev)
        except ValueError:
            raise InvalidParameterException("Invalid revision number. Number: " + rev)
        if (
            family in self.task_definitions
            and revision in self.task_definitions[family]
        ):
            # https://boto3.amazonaws.com/v1/documentation/api/latest/reference/services/ecs.html#ECS.Client.deregister_task_definition
            # At this time, INACTIVE task definitions remain discoverable in your account indefinitely.
            task_definition = self.task_definitions[family][revision]
            task_definition.status = "INACTIVE"
            return task_definition
        else:
            raise TaskDefinitionNotFoundException

    def run_task(
        self,
        cluster_str: str,
        task_definition_str: str,
        count: int,
        overrides: Optional[Dict[str, Any]],
        started_by: str,
        tags: Optional[List[Dict[str, str]]],
        launch_type: Optional[str],
        networking_configuration: Optional[Dict[str, Any]] = None,
    ) -> List[Task]:
        if launch_type and launch_type not in ["EC2", "FARGATE", "EXTERNAL"]:
            raise InvalidParameterException(
                "launch type should be one of [EC2,FARGATE,EXTERNAL]"
            )

        cluster = self._get_cluster(cluster_str)

        task_definition = self.describe_task_definition(task_definition_str)
        resource_requirements = self._calculate_task_resource_requirements(
            task_definition
        )
        if cluster.name not in self.tasks:
            self.tasks[cluster.name] = {}
        tasks = []
        if launch_type == "FARGATE":
            for _ in range(count):
                task = Task(
                    cluster=cluster,
                    task_definition=task_definition,
                    container_instance_arn=None,
                    resource_requirements=resource_requirements,
                    backend=self,
                    overrides=overrides or {},
                    started_by=started_by or "",
                    tags=tags or [],
                    launch_type=launch_type or "",
                    networking_configuration=networking_configuration,
                )
                tasks.append(task)
                self.tasks[cluster.name][task.task_arn] = task
            return tasks

        container_instances = list(
            self.container_instances.get(cluster.name, {}).keys()
        )
        if not container_instances:
            raise Exception(f"No instances found in cluster {cluster.name}")
        active_container_instances = [
            x
            for x in container_instances
            if self.container_instances[cluster.name][x].status == "ACTIVE"
        ]
        # TODO: return event about unable to place task if not able to place enough tasks to meet count
        placed_count = 0
        for name in active_container_instances:
            container_instance = self.container_instances[cluster.name][name]
            container_instance_arn = container_instance.container_instance_arn
            try_to_place = True
            while try_to_place:
                can_be_placed = self._can_be_placed(
                    container_instance, resource_requirements
                )
                if can_be_placed:
                    task = Task(
                        cluster,
                        task_definition,
                        container_instance_arn,
                        resource_requirements,
                        backend=self,
                        overrides=overrides or {},
                        started_by=started_by or "",
                        tags=tags or [],
                        launch_type=launch_type or "",
                        networking_configuration=networking_configuration,
                    )
                    self.update_container_instance_resources(
                        container_instance, resource_requirements
                    )
                    tasks.append(task)
                    self.tasks[cluster.name][task.task_arn] = task
                    placed_count += 1
                    if placed_count == count:
                        return tasks
                else:
                    try_to_place = False
        return tasks

    @staticmethod
    def _calculate_task_resource_requirements(  # type: ignore[misc]
        task_definition: TaskDefinition,
    ) -> Dict[str, Any]:
        resource_requirements: Dict[str, Any] = {
            "CPU": 0,
            "MEMORY": 0,
            "PORTS": [],
            "PORTS_UDP": [],
        }
        for container_definition in task_definition.container_definitions:
            # cloudformation uses capitalized properties, while boto uses all lower case

            # CPU is optional
            resource_requirements["CPU"] += container_definition.get(
                "cpu", container_definition.get("Cpu", 0)
            )

            # either memory or memory reservation must be provided
            if (
                "Memory" in container_definition
                or "MemoryReservation" in container_definition
            ):
                resource_requirements["MEMORY"] += container_definition.get(
                    "Memory", container_definition.get("MemoryReservation")
                )
            else:
                resource_requirements["MEMORY"] += container_definition.get(
                    "memory", container_definition.get("memoryReservation")
                )

            port_mapping_key = (
                "PortMappings"
                if "PortMappings" in container_definition
                else "portMappings"
            )
            for port_mapping in container_definition.get(port_mapping_key, []):  # type: ignore[attr-defined]
                if "hostPort" in port_mapping:
                    resource_requirements["PORTS"].append(port_mapping.get("hostPort"))
                elif "HostPort" in port_mapping:
                    resource_requirements["PORTS"].append(port_mapping.get("HostPort"))

        return resource_requirements

    @staticmethod
    def _can_be_placed(  # type: ignore[misc]
        container_instance: ContainerInstance,
        task_resource_requirements: Dict[str, Any],
    ) -> bool:
        """

        :param container_instance: The container instance trying to be placed onto
        :param task_resource_requirements: The calculated resource requirements of the task in the form of a dict
        :return: A boolean stating whether the given container instance has enough resources to have the task placed on
        it as well as a description, if it cannot be placed this will describe why.
        """
        # TODO: Implement default and other placement strategies as well as constraints:
        # docs.aws.amazon.com/AmazonECS/latest/developerguide/task-placement.html
        remaining_cpu = 0
        remaining_memory = 0
        reserved_ports: List[str] = []
        for resource in container_instance.remaining_resources:
            if resource.get("name") == "CPU":
                remaining_cpu = resource.get("integerValue")  # type: ignore[assignment]
            elif resource.get("name") == "MEMORY":
                remaining_memory = resource.get("integerValue")  # type: ignore[assignment]
            elif resource.get("name") == "PORTS":
                reserved_ports = resource.get("stringSetValue")  # type: ignore[assignment]
        if task_resource_requirements.get("CPU") > remaining_cpu:  # type: ignore[operator]
            return False
        if task_resource_requirements.get("MEMORY") > remaining_memory:  # type: ignore[operator]
            return False
        ports_needed = task_resource_requirements.get("PORTS")
        for port in ports_needed:  # type: ignore[union-attr]
            if str(port) in reserved_ports:
                return False
        return True

    def start_task(
        self,
        cluster_str: str,
        task_definition_str: str,
        container_instances: List[str],
        overrides: Dict[str, Any],
        started_by: str,
        tags: Optional[List[Dict[str, str]]] = None,
    ) -> List[Task]:
        cluster = self._get_cluster(cluster_str)

        task_definition = self.describe_task_definition(task_definition_str)
        if cluster.name not in self.tasks:
            self.tasks[cluster.name] = {}
        tasks = []
        if not container_instances:
            raise EcsClientException("Container Instances cannot be empty.")

        container_instance_ids = [x.split("/")[-1] for x in container_instances]
        resource_requirements = self._calculate_task_resource_requirements(
            task_definition
        )
        for container_instance_id in container_instance_ids:
            container_instance = self.container_instances[cluster.name][
                container_instance_id
            ]
            task = Task(
                cluster,
                task_definition,
                container_instance.container_instance_arn,
                resource_requirements,
                backend=self,
                overrides=overrides or {},
                started_by=started_by or "",
                tags=tags,
            )
            tasks.append(task)
            self.update_container_instance_resources(
                container_instance, resource_requirements
            )
            self.tasks[cluster.name][task.task_arn] = task
        return tasks

    def describe_tasks(self, cluster_str: str, tasks: Optional[str]) -> List[Task]:
        """
        Only include=TAGS is currently supported.
        """
        self._get_cluster(cluster_str)

        if not tasks:
            raise InvalidParameterException("Tasks cannot be empty.")
        response = []
        for cluster_tasks in self.tasks.values():
            for task_arn, task in cluster_tasks.items():
                task_id = task_arn.split("/")[-1]
                if (
                    task_arn in tasks
                    or task.task_arn in tasks
                    or any(task_id in task for task in tasks)
                ):
                    task.advance()
                    response.append(task)

        return response

    def list_tasks(
        self,
        cluster_str: Optional[str] = None,
        container_instance: Optional[str] = None,
        family: Optional[str] = None,
        started_by: Optional[str] = None,
        service_name: Optional[str] = None,
        desiredStatus: Optional[str] = None,
    ) -> List[Task]:
        filtered_tasks = []
        for tasks in self.tasks.values():
            for task in tasks.values():
                filtered_tasks.append(task)
        if cluster_str:
            cluster = self._get_cluster(cluster_str)

            filtered_tasks = list(
                filter(lambda t: cluster.name in t.cluster_arn, filtered_tasks)
            )

        if container_instance:
            filtered_tasks = list(
                filter(
                    lambda t: container_instance in t.container_instance_arn,  # type: ignore
                    filtered_tasks,
                )
            )

        if family:
            task_definition_arns = self.list_task_definitions(family)
            filtered_tasks = list(
                filter(
                    lambda t: t.task_definition_arn in task_definition_arns,
                    filtered_tasks,
                )
            )

        if started_by:
            filtered_tasks = list(
                filter(lambda t: started_by == t.started_by, filtered_tasks)
            )

        if service_name:
            # TODO: We can't filter on `service_name` until the backend actually
            # launches tasks as part of the service creation process.
            pass

        if desiredStatus:
            filtered_tasks = list(
                filter(lambda t: t.desired_status == desiredStatus, filtered_tasks)
            )

        return filtered_tasks

    def stop_task(self, cluster_str: str, task_str: str, reason: str) -> Task:
        cluster = self._get_cluster(cluster_str)

        task_id = task_str.split("/")[-1]
        tasks = self.tasks.get(cluster.name, None)
        if not tasks:
            raise Exception(f"Cluster {cluster.name} has no registered tasks")
        for task in tasks.keys():
            if task.endswith(task_id):
                container_instance_arn = tasks[task].container_instance_arn
                if container_instance_arn:
                    container_instance = self.container_instances[cluster.name][
                        container_instance_arn.split("/")[-1]
                    ]
                    self.update_container_instance_resources(
                        container_instance,
                        tasks[task].resource_requirements,  # type: ignore[arg-type]
                        removing=True,
                    )
                tasks[task].last_status = "STOPPED"
                tasks[task].desired_status = "STOPPED"
                tasks[task].stopped_reason = reason
                return tasks[task]
        raise Exception(f"Could not find task {task_str} on cluster {cluster.name}")

    def _get_service(self, cluster_str: str, service_str: str) -> Service:
        cluster = self._get_cluster(cluster_str)
        for service in self.services.values():
            if service.cluster_name == cluster.name and (
                service.name == service_str or service.arn == service_str
            ):
                return service
        raise ServiceNotFoundException

    def create_service(
        self,
        cluster_str: str,
        service_name: str,
        desired_count: int,
        task_definition_str: Optional[str] = None,
        load_balancers: Optional[List[Dict[str, Any]]] = None,
        scheduling_strategy: Optional[List[Dict[str, Any]]] = None,
        tags: Optional[List[Dict[str, str]]] = None,
        deployment_controller: Optional[Dict[str, str]] = None,
        launch_type: Optional[str] = None,
        service_registries: Optional[List[Dict[str, Any]]] = None,
        platform_version: Optional[str] = None,
    ) -> Service:
        cluster = self._get_cluster(cluster_str)

        if task_definition_str:
            task_definition = self.describe_task_definition(task_definition_str)
        else:
            task_definition = None
        desired_count = desired_count if desired_count is not None else 0

        launch_type = launch_type if launch_type is not None else "EC2"
        if launch_type not in ["EC2", "FARGATE"]:
            raise EcsClientException("launch type should be one of [EC2,FARGATE]")

        service = Service(
            cluster=cluster,
            service_name=service_name,
            desired_count=desired_count,
            task_definition=task_definition,
            load_balancers=load_balancers,
            scheduling_strategy=scheduling_strategy,
            tags=tags,
            deployment_controller=deployment_controller,
            launch_type=launch_type,
            backend=self,
            service_registries=service_registries,
            platform_version=platform_version,
        )
        cluster_service_pair = f"{cluster.name}:{service_name}"
        self.services[cluster_service_pair] = service

        return service

    def list_services(
        self,
        cluster_str: str,
        scheduling_strategy: Optional[str] = None,
        launch_type: Optional[str] = None,
    ) -> List[str]:
        cluster = self._get_cluster(cluster_str)
        service_arns = []
        for key, service in self.services.items():
            if cluster.name + ":" not in key:
                continue

            if (
                scheduling_strategy is not None
                and service.scheduling_strategy != scheduling_strategy
            ):
                continue

            if launch_type is not None and service.launch_type != launch_type:
                continue

            service_arns.append(service.arn)

        return sorted(service_arns)

    def describe_services(
        self, cluster_str: str, service_names_or_arns: List[str]
    ) -> Tuple[List[Service], List[Dict[str, str]]]:
        cluster = self._get_cluster(cluster_str)

        result = []
        failures = []
        for name_or_arn in service_names_or_arns:
            name = name_or_arn.split("/")[-1]
            cluster_service_pair = f"{cluster.name}:{name}"
            if cluster_service_pair in self.services:
                result.append(self.services[cluster_service_pair])
            else:
                if re.match(ARN_PARTITION_REGEX + ":ecs", name_or_arn):
                    missing_arn = name_or_arn
                else:
                    missing_arn = f"arn:{get_partition(self.region_name)}:ecs:{self.region_name}:{self.account_id}:service/{name}"
                failures.append({"arn": missing_arn, "reason": "MISSING"})

        return result, failures

    def update_service(self, service_properties: Dict[str, Any]) -> Service:
        cluster_str = service_properties.pop("cluster", "default")
        task_definition_str = service_properties.pop("task_definition", None)
        cluster = self._get_cluster(cluster_str)
        service_name = service_properties.pop("service").split("/")[-1]
        cluster_service_pair = f"{cluster.name}:{service_name}"

        if cluster_service_pair in self.services:
            current_service = self.services[cluster_service_pair]
            for prop_name, prop_val in service_properties.items():
                if prop_val is not None:
                    current_service.__setattr__(prop_name, prop_val)
                    if prop_name == "desired_count":
                        current_service.__setattr__("running_count", prop_val)
                        current_service.__setattr__("pending_count", 0)
            if task_definition_str:
                self.describe_task_definition(task_definition_str)
                current_service.task_definition = task_definition_str
            return current_service
        else:
            raise ServiceNotFoundException

    def delete_service(
        self, cluster_name: str, service_name: str, force: bool
    ) -> Service:
        cluster = self._get_cluster(cluster_name)
        service = self._get_service(cluster_name, service_name)

        cluster_service_pair = f"{cluster.name}:{service.name}"

        service = self.services[cluster_service_pair]
        if service.desired_count > 0 and not force:
            raise InvalidParameterException(
                "The service cannot be stopped while it is scaled above 0."
            )
        else:
            # A service is not immediately removed - just marked as inactive
            # It is only deleted later on
            # https://boto3.amazonaws.com/v1/documentation/api/latest/reference/services/ecs.html#ECS.Client.delete_service
            service.status = "INACTIVE"
            service.pending_count = 0
            return service

    def register_container_instance(
        self, cluster_str: str, ec2_instance_id: str
    ) -> ContainerInstance:
        cluster_name = cluster_str.split("/")[-1]
        if cluster_name not in self.clusters:
            raise Exception(f"{cluster_name} is not a cluster")
        container_instance = ContainerInstance(
            ec2_instance_id,
            self.account_id,
            self.region_name,
            cluster_name,
            backend=self,
        )
        if not self.container_instances.get(cluster_name):
            self.container_instances[cluster_name] = {}
        container_instance_id = container_instance.container_instance_arn.split("/")[-1]
        self.container_instances[cluster_name][container_instance_id] = (
            container_instance
        )
        self.clusters[cluster_name].registered_container_instances_count += 1
        return container_instance

    def list_container_instances(self, cluster_str: str) -> List[str]:
        cluster_name = cluster_str.split("/")[-1]
        container_instances_values = self.container_instances.get(
            cluster_name, {}
        ).values()
        container_instances = [
            ci.container_instance_arn for ci in container_instances_values
        ]
        return sorted(container_instances)

    def describe_container_instances(
        self, cluster_str: str, list_container_instance_ids: List[str]
    ) -> Tuple[List[ContainerInstance], List[ContainerInstanceFailure]]:
        cluster = self._get_cluster(cluster_str)

        if not list_container_instance_ids:
            raise EcsClientException("Container Instances cannot be empty.")
        failures = []
        container_instance_objects = []
        for container_instance_id in list_container_instance_ids:
            container_instance_id = container_instance_id.split("/")[-1]
            container_instance = self.container_instances[cluster.name].get(
                container_instance_id, None
            )
            if container_instance is not None:
                container_instance_objects.append(container_instance)
            else:
                failures.append(
                    ContainerInstanceFailure(
                        "MISSING",
                        container_instance_id,
                        self.account_id,
                        self.region_name,
                    )
                )

        return container_instance_objects, failures

    def update_container_instances_state(
        self, cluster_str: str, list_container_instance_ids: List[str], status: str
    ) -> Tuple[List[ContainerInstance], List[ContainerInstanceFailure]]:
        cluster = self._get_cluster(cluster_str)

        status = status.upper()
        if status not in ["ACTIVE", "DRAINING"]:
            raise InvalidParameterException(
                "Container instance status should be one of [ACTIVE, DRAINING]"
            )
        failures = []
        container_instance_objects = []
        list_container_instance_ids = [
            x.split("/")[-1] for x in list_container_instance_ids
        ]
        for container_instance_id in list_container_instance_ids:
            container_instance = self.container_instances[cluster.name].get(
                container_instance_id, None
            )
            if container_instance is not None:
                container_instance.status = status
                container_instance_objects.append(container_instance)
            else:
                failures.append(
                    ContainerInstanceFailure(
                        "MISSING",
                        container_instance_id,
                        self.account_id,
                        self.region_name,
                    )
                )

        return container_instance_objects, failures

    def update_container_instance_resources(
        self,
        container_instance: ContainerInstance,
        task_resources: Dict[str, Any],
        removing: bool = False,
    ) -> None:
        resource_multiplier = 1
        if removing:
            resource_multiplier = -1
        for resource in container_instance.remaining_resources:
            if resource.get("name") == "CPU":
                resource["integerValue"] -= (
                    task_resources.get("CPU") * resource_multiplier  # type: ignore[operator]
                )
            elif resource.get("name") == "MEMORY":
                resource["integerValue"] -= (
                    task_resources.get("MEMORY") * resource_multiplier  # type: ignore[operator]
                )
            elif resource.get("name") == "PORTS":
                for port in task_resources.get("PORTS"):  # type: ignore[union-attr]
                    if removing:
                        resource["stringSetValue"].remove(str(port))
                    else:
                        resource["stringSetValue"].append(str(port))
        container_instance.running_tasks_count += resource_multiplier * 1

    def deregister_container_instance(
        self, cluster_str: str, container_instance_str: str, force: bool
    ) -> ContainerInstance:
        cluster = self._get_cluster(cluster_str)

        container_instance_id = container_instance_str.split("/")[-1]
        container_instance = self.container_instances[cluster.name].get(
            container_instance_id
        )
        if container_instance is None:
            raise Exception("{0} is not a container id in the cluster")
        if not force and container_instance.running_tasks_count > 0:
            raise JsonRESTError(
                error_type="InvalidParameter",
                message="Found running tasks on the instance.",
            )
        # Currently assume that people might want to do something based around deregistered instances
        # with tasks left running on them - but nothing if no tasks were running already
        elif force and container_instance.running_tasks_count > 0:
            if not self.container_instances.get("orphaned"):
                self.container_instances["orphaned"] = {}
            self.container_instances["orphaned"][container_instance_id] = (
                container_instance
            )
        del self.container_instances[cluster.name][container_instance_id]
        self._respond_to_cluster_state_update(cluster_str)
        return container_instance

    def _respond_to_cluster_state_update(self, cluster_str: str) -> None:
        self._get_cluster(cluster_str)

        pass

    def put_attributes(
        self, cluster_name: str, attributes: Optional[List[Dict[str, Any]]] = None
    ) -> None:
        cluster = self._get_cluster(cluster_name)

        if attributes is None:
            raise InvalidParameterException("attributes can not be empty")

        for attr in attributes:
            self._put_attribute(
                cluster.name,
                attr["name"],
                attr.get("value"),
                attr.get("targetId"),
                attr.get("targetType"),
            )

    def _put_attribute(
        self,
        cluster_name: str,
        name: str,
        value: Optional[str] = None,
        target_id: Optional[str] = None,
        target_type: Optional[str] = None,
    ) -> None:
        if target_id is None and target_type is None:
            for instance in self.container_instances[cluster_name].values():
                instance.attributes[name] = value
        elif target_type is None:
            # targetId is full container instance arn
            try:
                arn = target_id.rsplit("/", 1)[-1]  # type: ignore[union-attr]
                self.container_instances[cluster_name][arn].attributes[name] = value
            except KeyError:
                raise JsonRESTError(
                    "TargetNotFoundException", f"Could not find {target_id}"
                )
        else:
            # targetId is container uuid, targetType must be container-instance
            try:
                if target_type != "container-instance":
                    raise JsonRESTError(
                        "TargetNotFoundException", f"Could not find {target_id}"
                    )

                self.container_instances[cluster_name][target_id].attributes[  # type: ignore[index]
                    name
                ] = value
            except KeyError:
                raise JsonRESTError(
                    "TargetNotFoundException", f"Could not find {target_id}"
                )

    def list_attributes(
        self,
        target_type: str,
        cluster_name: Optional[str] = None,
        attr_name: Optional[str] = None,
        attr_value: Optional[str] = None,
    ) -> Any:
        """
        Pagination is not yet implemented
        """
        if target_type != "container-instance":
            raise JsonRESTError(
                "InvalidParameterException", "targetType must be container-instance"
            )

        filters = [lambda x: True]

        # item will be {0 cluster_name, 1 arn, 2 name, 3 value}
        if cluster_name is not None:
            filters.append(lambda item: item[0] == cluster_name)
        if attr_name:
            filters.append(lambda item: item[2] == attr_name)
        if attr_name:
            filters.append(lambda item: item[3] == attr_value)

        all_attrs = []
        for cluster_name, cobj in self.container_instances.items():
            for container_instance in cobj.values():
                for key, value in container_instance.attributes.items():
                    all_attrs.append(
                        (
                            cluster_name,
                            container_instance.container_instance_arn,
                            key,
                            value,
                        )
                    )

        return filter(lambda x: all(f(x) for f in filters), all_attrs)  # type: ignore

    def delete_attributes(
        self, cluster_name: str, attributes: Optional[List[Dict[str, Any]]] = None
    ) -> None:
        cluster = self._get_cluster(cluster_name)

        if attributes is None:
            raise JsonRESTError(
                "InvalidParameterException", "attributes value is required"
            )

        for attr in attributes:
            self._delete_attribute(
                cluster.name,
                attr["name"],
                attr.get("value"),
                attr.get("targetId"),
                attr.get("targetType"),
            )

    def _delete_attribute(
        self,
        cluster_name: str,
        name: str,
        value: Optional[str] = None,
        target_id: Optional[str] = None,
        target_type: Optional[str] = None,
    ) -> None:
        if target_id is None and target_type is None:
            for instance in self.container_instances[cluster_name].values():
                if name in instance.attributes and instance.attributes[name] == value:
                    del instance.attributes[name]
        elif target_type is None:
            # targetId is full container instance arn
            try:
                arn = target_id.rsplit("/", 1)[-1]  # type: ignore[union-attr]
                instance = self.container_instances[cluster_name][arn]
                if name in instance.attributes and instance.attributes[name] == value:
                    del instance.attributes[name]
            except KeyError:
                raise JsonRESTError(
                    "TargetNotFoundException", f"Could not find {target_id}"
                )
        else:
            # targetId is container uuid, targetType must be container-instance
            try:
                if target_type != "container-instance":
                    raise JsonRESTError(
                        "TargetNotFoundException", f"Could not find {target_id}"
                    )

                instance = self.container_instances[cluster_name][target_id]  # type: ignore[index]
                if name in instance.attributes and instance.attributes[name] == value:
                    del instance.attributes[name]
            except KeyError:
                raise JsonRESTError(
                    "TargetNotFoundException", f"Could not find {target_id}"
                )

    def list_task_definition_families(
        self, family_prefix: Optional[str] = None
    ) -> Iterator[str]:
        """
        The Status and pagination parameters are not yet implemented
        """
        for task_fam in self.task_definitions:
            if family_prefix is not None and not task_fam.startswith(family_prefix):
                continue

            yield task_fam

    @staticmethod
    def _parse_resource_arn(resource_arn: str) -> Dict[str, str]:
        regexes = [
            ARN_PARTITION_REGEX
            + ":ecs:(?P<region>[^:]+):(?P<account_id>[^:]+):(?P<service>[^:]+)/(?P<cluster_id>[^:]+)/(?P<service_id>[^:]+)/ecs-svc/(?P<id>.*)$",
            ARN_PARTITION_REGEX
            + ":ecs:(?P<region>[^:]+):(?P<account_id>[^:]+):(?P<service>[^:]+)/(?P<cluster_id>[^:]+)/(?P<id>.*)$",
            ARN_PARTITION_REGEX
            + ":ecs:(?P<region>[^:]+):(?P<account_id>[^:]+):(?P<service>[^:]+)/(?P<id>.*)$",
        ]
        for regex in regexes:
            match = re.match(regex, resource_arn)
            if match:
                return match.groupdict()
        raise JsonRESTError("InvalidParameterException", "The ARN provided is invalid.")

    def _get_resource(self, resource_arn: str, parsed_arn: Dict[str, str]) -> Any:
        if parsed_arn["service"] == "cluster":
            return self._get_cluster(parsed_arn["id"])
        if parsed_arn["service"] == "service":
            for service in self.services.values():
                if service.arn == resource_arn:
                    return service
            raise ServiceNotFoundException
        elif parsed_arn["service"] == "task-set":
            c_id = parsed_arn["cluster_id"]
            s_id = parsed_arn["service_id"]
            services, _ = self.describe_services(
                cluster_str=c_id, service_names_or_arns=[s_id]
            )
            for service in services:
                for task_set in service.task_sets:
                    if task_set.task_set_arn == resource_arn:
                        return task_set
            raise ServiceNotFoundException
        elif parsed_arn["service"] == "task-definition":
            task_def = self.describe_task_definition(
                task_definition_str=parsed_arn["id"]
            )
            return task_def
        elif parsed_arn["service"] == "capacity-provider":
            return self._get_provider(parsed_arn["id"])
        elif parsed_arn["service"] == "task":
            for task in self.list_tasks():
                if task.task_arn == resource_arn:
                    return task
        raise NotImplementedError()

    def list_tags_for_resource(self, resource_arn: str) -> List[Dict[str, str]]:
        """Currently implemented only for task definitions and services"""
        parsed_arn = self._parse_resource_arn(resource_arn)
        resource = self._get_resource(resource_arn, parsed_arn)
        return resource.tags

    def _get_last_task_definition_revision_id(self, family: str) -> int:  # type: ignore[return]
        definitions = self.task_definitions.get(family)
        if definitions:
            return max(definitions.keys())

    def tag_resource(self, resource_arn: str, tags: List[Dict[str, str]]) -> None:
        parsed_arn = self._parse_resource_arn(resource_arn)
        resource = self._get_resource(resource_arn, parsed_arn)
        resource.tags = self._merge_tags(resource.tags or [], tags)

    def _merge_tags(
        self, existing_tags: List[Dict[str, str]], new_tags: List[Dict[str, str]]
    ) -> List[Dict[str, str]]:
        merged_tags = new_tags
        new_keys = self._get_keys(new_tags)
        for existing_tag in existing_tags:
            if existing_tag["key"] not in new_keys:
                merged_tags.append(existing_tag)
        return merged_tags

    @staticmethod
    def _get_keys(tags: List[Dict[str, str]]) -> List[str]:
        return [tag["key"] for tag in tags]

    def untag_resource(self, resource_arn: str, tag_keys: List[str]) -> None:
        parsed_arn = self._parse_resource_arn(resource_arn)
        resource = self._get_resource(resource_arn, parsed_arn)
        resource.tags = [tag for tag in resource.tags if tag["key"] not in tag_keys]

    def create_task_set(
        self,
        service: str,
        cluster_str: str,
        task_definition: str,
        external_id: Optional[str] = None,
        network_configuration: Optional[Dict[str, Any]] = None,
        load_balancers: Optional[List[Dict[str, Any]]] = None,
        service_registries: Optional[List[Dict[str, Any]]] = None,
        launch_type: Optional[str] = None,
        capacity_provider_strategy: Optional[List[Dict[str, Any]]] = None,
        platform_version: Optional[str] = None,
        scale: Optional[Dict[str, Any]] = None,
        client_token: Optional[str] = None,
        tags: Optional[List[Dict[str, str]]] = None,
    ) -> TaskSet:
        launch_type = launch_type if launch_type is not None else "EC2"
        if launch_type not in ["EC2", "FARGATE"]:
            raise EcsClientException("launch type should be one of [EC2,FARGATE]")

        task_set = TaskSet(
            service,
            cluster_str,
            task_definition,
            self.account_id,
            self.region_name,
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

        service_name = service.split("/")[-1]

        cluster_obj = self._get_cluster(cluster_str)
        service_obj = self.services.get(f"{cluster_obj.name}:{service_name}")
        if not service_obj:
            raise ServiceNotFoundException

        task_set.task_definition = self.describe_task_definition(task_definition).arn
        task_set.service_arn = service_obj.arn
        task_set.cluster_arn = cluster_obj.arn

        service_obj.task_sets.append(task_set)
        # TODO: validate load balancers

        return task_set

    def describe_task_sets(
        self, cluster_str: str, service: str, task_sets: Optional[List[str]] = None
    ) -> List[TaskSet]:
        task_sets = task_sets or []

        cluster_obj = self._get_cluster(cluster_str)

        service_name = service.split("/")[-1]
        service_key = f"{cluster_obj.name}:{service_name}"

        service_obj = self.services.get(service_key)
        if not service_obj:
            raise ServiceNotFoundException

        task_set_results = []
        if task_sets:
            for task_set in service_obj.task_sets:
                # Match full ARN
                if task_set.task_set_arn in task_sets:
                    task_set_results.append(task_set)
                # Match partial ARN if only the taskset ID is provided
                elif "/".join(task_set.task_set_arn.split("/")[-2:]) in task_sets:
                    task_set_results.append(task_set)
        else:
            task_set_results = service_obj.task_sets

        return task_set_results

    def delete_task_set(self, cluster: str, service: str, task_set: str) -> TaskSet:
        """
        The Force-parameter is not yet implemented
        """
        cluster_name = cluster.split("/")[-1]
        service_name = service.split("/")[-1]

        service_key = f"{cluster_name}:{service_name}"
        task_set_element = None
        for i, ts in enumerate(self.services[service_key].task_sets):
            if task_set == ts.task_set_arn or task_set == "/".join(
                ts.task_set_arn.split("/")[-2:]
            ):
                task_set_element = i

        if task_set_element is not None:
            deleted_task_set = self.services[service_key].task_sets.pop(
                task_set_element
            )
        else:
            raise TaskSetNotFoundException

        # TODO: add logic for `force` to raise an exception if `PRIMARY` task has not been scaled to 0.

        return deleted_task_set

    def update_task_set(
        self, cluster: str, service: str, task_set: str, scale: Dict[str, Any]
    ) -> TaskSet:
        cluster_name = cluster.split("/")[-1]
        service_name = service.split("/")[-1]
        task_set_obj = self.describe_task_sets(
            cluster_name, service_name, task_sets=[task_set]
        )[0]
        task_set_obj.scale = scale
        return task_set_obj

    def update_service_primary_task_set(
        self, cluster: str, service: str, primary_task_set: str
    ) -> TaskSet:
        """Updates task sets be PRIMARY or ACTIVE for given cluster:service task sets"""
        cluster_name = cluster.split("/")[-1]
        service_name = service.split("/")[-1]
        task_set_obj = self.describe_task_sets(
            cluster_name, service_name, task_sets=[primary_task_set]
        )[0]

        services, _ = self.describe_services(cluster, [service])
        service_obj = services[0]
        service_obj.load_balancers = task_set_obj.load_balancers
        service_obj.task_definition = task_set_obj.task_definition

        for task_set in service_obj.task_sets:
            if task_set.task_set_arn == primary_task_set:
                task_set.status = "PRIMARY"
            else:
                task_set.status = "ACTIVE"
        return task_set_obj

    def list_account_settings(
        self, name: Optional[str] = None, value: Optional[str] = None
    ) -> List[AccountSetting]:
        expected_names = [
            "serviceLongArnFormat",
            "taskLongArnFormat",
            "containerInstanceLongArnFormat",
            "containerLongArnFormat",
            "awsvpcTrunking",
            "containerInsights",
            "dualStackIPv6",
        ]
        if name and name not in expected_names:
            raise UnknownAccountSettingException()
        all_settings = self.account_settings.values()
        return [
            s
            for s in all_settings
            if (not name or s.name == name) and (not value or s.value == value)
        ]

    def put_account_setting(self, name: str, value: str) -> AccountSetting:
        account_setting = AccountSetting(name, value)
        self.account_settings[name] = account_setting
        return account_setting

    def delete_account_setting(self, name: str) -> None:
        self.account_settings.pop(name, None)

    def enable_long_arn_for_name(self, name: str) -> bool:
        account = self.account_settings.get(name, None)
        if account and account.value == "disabled":
            return False
        return settings.ecs_new_arn_format()


ecs_backends = BackendDict(EC2ContainerServiceBackend, "ecs")
