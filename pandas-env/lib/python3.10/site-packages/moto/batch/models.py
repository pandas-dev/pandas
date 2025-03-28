import datetime
import logging
import re
import threading
import time
from itertools import cycle
from sys import platform
from time import sleep
from typing import Any, Dict, List, Optional, Set, Tuple, Union

import dateutil.parser

from moto import settings
from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import BaseModel, CloudFormationModel
from moto.core.utils import unix_time_millis
from moto.ec2.exceptions import InvalidSubnetIdError
from moto.ec2.models import EC2Backend, ec2_backends
from moto.ec2.models.instance_types import INSTANCE_FAMILIES as EC2_INSTANCE_FAMILIES
from moto.ec2.models.instance_types import INSTANCE_TYPES as EC2_INSTANCE_TYPES
from moto.ec2.models.instances import Instance
from moto.ecs.models import EC2ContainerServiceBackend, ecs_backends
from moto.iam.exceptions import IAMNotFoundException
from moto.iam.models import IAMBackend, iam_backends
from moto.logs.models import LogsBackend, logs_backends
from moto.moto_api._internal import mock_random
from moto.moto_api._internal.managed_state_model import ManagedState
from moto.utilities.docker_utilities import DockerModel
from moto.utilities.tagging_service import TaggingService
from moto.utilities.utils import get_partition

from .exceptions import ClientException, InvalidParameterValueException, ValidationError
from .utils import (
    JobStatus,
    lowercase_first_key,
    make_arn_for_compute_env,
    make_arn_for_job,
    make_arn_for_job_queue,
    make_arn_for_task_def,
)

logger = logging.getLogger(__name__)
COMPUTE_ENVIRONMENT_NAME_REGEX = re.compile(
    r"^[A-Za-z0-9][A-Za-z0-9_-]{1,126}[A-Za-z0-9]$"
)
JOB_NAME_REGEX = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_-]{1,127}$")


def datetime2int_milliseconds(date: datetime.datetime) -> int:
    """
    AWS returns timestamps in milliseconds
    We don't use milliseconds timestamps internally,
    this method should be used only in describe() method
    """
    return int(date.timestamp() * 1000)


def datetime2int(date: datetime.datetime) -> int:
    return int(time.mktime(date.timetuple()))


class ComputeEnvironment(CloudFormationModel):
    def __init__(
        self,
        compute_environment_name: str,
        _type: str,
        state: str,
        compute_resources: Dict[str, Any],
        service_role: str,
        account_id: str,
        region_name: str,
    ):
        self.name = compute_environment_name
        self.env_type = _type
        self.state = state
        self.compute_resources = compute_resources
        self.service_role = service_role
        self.arn = make_arn_for_compute_env(
            account_id, compute_environment_name, region_name
        )

        self.instances: List[Instance] = []
        self.ecs_arn = ""
        self.ecs_name = ""

    def add_instance(self, instance: Instance) -> None:
        self.instances.append(instance)

    def set_ecs(self, arn: str, name: str) -> None:
        self.ecs_arn = arn
        self.ecs_name = name

    @property
    def physical_resource_id(self) -> str:
        return self.arn

    @staticmethod
    def cloudformation_name_type() -> str:
        return "ComputeEnvironmentName"

    @staticmethod
    def cloudformation_type() -> str:
        # https://docs.aws.amazon.com/AWSCloudFormation/latest/UserGuide/aws-resource-batch-computeenvironment.html
        return "AWS::Batch::ComputeEnvironment"

    @classmethod
    def create_from_cloudformation_json(  # type: ignore[misc]
        cls,
        resource_name: str,
        cloudformation_json: Dict[str, Any],
        account_id: str,
        region_name: str,
        **kwargs: Any,
    ) -> "ComputeEnvironment":
        backend = batch_backends[account_id][region_name]
        properties = cloudformation_json["Properties"]

        return backend.create_compute_environment(
            resource_name,
            properties["Type"],
            properties.get("State", "ENABLED"),
            lowercase_first_key(properties["ComputeResources"]),
            properties["ServiceRole"],
        )


class JobQueue(CloudFormationModel):
    def __init__(
        self,
        name: str,
        priority: str,
        state: str,
        environments: List[ComputeEnvironment],
        env_order_json: List[Dict[str, Any]],
        schedule_policy: Optional[str],
        backend: "BatchBackend",
        tags: Optional[Dict[str, str]] = None,
    ):
        """
        :param name: Job queue name
        :type name: str
        :param priority: Job queue priority
        :type priority: int
        :param state: Either ENABLED or DISABLED
        :type state: str
        :param environments: Compute Environments
        :type environments: list of ComputeEnvironment
        :param env_order_json: Compute Environments JSON for use when describing
        :type env_order_json: list of dict
        """
        self.name = name
        self.priority = priority
        self.state = state
        self.environments = environments
        self.env_order_json = env_order_json
        self.schedule_policy = schedule_policy
        self.arn = make_arn_for_job_queue(backend.account_id, name, backend.region_name)
        self.status = "VALID"
        self.backend = backend

        if tags:
            backend.tag_resource(self.arn, tags)

        self.jobs: List[Job] = []

    def describe(self) -> Dict[str, Any]:
        return {
            "computeEnvironmentOrder": self.env_order_json,
            "jobQueueArn": self.arn,
            "jobQueueName": self.name,
            "priority": self.priority,
            "schedulingPolicyArn": self.schedule_policy,
            "state": self.state,
            "status": self.status,
            "tags": self.backend.list_tags_for_resource(self.arn),
        }

    @property
    def physical_resource_id(self) -> str:
        return self.arn

    @staticmethod
    def cloudformation_name_type() -> str:
        return "JobQueueName"

    @staticmethod
    def cloudformation_type() -> str:
        # https://docs.aws.amazon.com/AWSCloudFormation/latest/UserGuide/aws-resource-batch-jobqueue.html
        return "AWS::Batch::JobQueue"

    @classmethod
    def create_from_cloudformation_json(  # type: ignore[misc]
        cls,
        resource_name: str,
        cloudformation_json: Dict[str, Any],
        account_id: str,
        region_name: str,
        **kwargs: Any,
    ) -> "JobQueue":
        backend = batch_backends[account_id][region_name]
        properties = cloudformation_json["Properties"]

        # Need to deal with difference case from cloudformation compute_resources, e.g. instanceRole vs InstanceRole
        # Hacky fix to normalise keys, is making me think I want to start spamming cAsEiNsEnSiTiVe dictionaries
        compute_envs = [
            lowercase_first_key(dict_item)
            for dict_item in properties["ComputeEnvironmentOrder"]
        ]

        return backend.create_job_queue(
            queue_name=resource_name,
            priority=properties["Priority"],
            state=properties.get("State", "ENABLED"),
            compute_env_order=compute_envs,
            schedule_policy=None,
        )


class JobDefinition(CloudFormationModel):
    def __init__(
        self,
        name: str,
        parameters: Optional[Dict[str, Any]],
        _type: str,
        container_properties: Dict[str, Any],
        node_properties: Dict[str, Any],
        tags: Dict[str, str],
        retry_strategy: Dict[str, str],
        timeout: Dict[str, int],
        backend: "BatchBackend",
        platform_capabilities: List[str],
        propagate_tags: bool,
        revision: Optional[int] = 0,
    ):
        self.name = name
        self.retry_strategy = retry_strategy
        self.type = _type
        self.revision = revision or 0
        self._region = backend.region_name
        self.container_properties = container_properties
        self.node_properties = node_properties
        self.status = "ACTIVE"
        self.parameters = parameters or {}
        self.timeout = timeout
        self.backend = backend
        self.platform_capabilities = platform_capabilities
        self.propagate_tags = propagate_tags

        if self.container_properties is not None:
            # Set some default values
            default_values: Dict[str, List[Any]] = {
                "command": [],
                "resourceRequirements": [],
                "secrets": [],
                "environment": [],
                "mountPoints": [],
                "ulimits": [],
                "volumes": [],
            }
            for key, val in default_values.items():
                if key not in self.container_properties:
                    self.container_properties[key] = val

            # Set default FARGATE configuration
            if "FARGATE" in (self.platform_capabilities or []):
                if "fargatePlatformConfiguration" not in self.container_properties:
                    self.container_properties["fargatePlatformConfiguration"] = {
                        "platformVersion": "LATEST"
                    }

            # Remove any empty environment variables
            self.container_properties["environment"] = [
                env_var
                for env_var in self.container_properties["environment"]
                if env_var.get("value") != ""
            ]

        self._validate()
        self.revision += 1
        self.arn = make_arn_for_task_def(
            self.backend.account_id, self.name, self.revision, self._region
        )

        tag_list = self._format_tags(tags or {})
        # Validate the tags before proceeding.
        errmsg = self.backend.tagger.validate_tags(tag_list)
        if errmsg:
            raise ValidationError(errmsg)

        self.backend.tagger.tag_resource(self.arn, tag_list)

    def _format_tags(self, tags: Dict[str, str]) -> List[Dict[str, str]]:
        return [{"Key": k, "Value": v} for k, v in tags.items()]

    def _get_resource_requirement(self, req_type: str, default: Any = None) -> Any:
        """
        Get resource requirement from container properties.

        Resource requirements like "memory" and "vcpus" are now specified in
        "resourceRequirements". This function retrieves a resource requirement
        from either container_properties.resourceRequirements (preferred) or
        directly from container_properties (deprecated).

        :param req_type: The type of resource requirement to retrieve.
        :type req_type: ["gpu", "memory", "vcpus"]

        :param default: The default value to return if the resource requirement is not found.
        :type default: any, default=None

        :return: The value of the resource requirement, or None.
        :rtype: any
        """
        resource_reqs = self.container_properties.get("resourceRequirements", [])

        # Filter the resource requirements by the specified type.
        # Note that VCPUS are specified in resourceRequirements without the
        # trailing "s", so we strip that off in the comparison below.
        required_resource = list(
            filter(
                lambda req: req["type"].lower() == req_type.lower().rstrip("s"),
                resource_reqs,
            )
        )

        if required_resource:
            if req_type == "vcpus":
                return float(required_resource[0]["value"])
            elif req_type == "memory":
                return int(required_resource[0]["value"])
            else:
                return required_resource[0]["value"]
        else:
            return self.container_properties.get(req_type, default)

    def _validate(self) -> None:
        # For future use when containers arnt the only thing in batch
        VALID_TYPES = ("container", "multinode")
        if self.type not in VALID_TYPES:
            raise ClientException(f"type must be one of {VALID_TYPES}")

        if not isinstance(self.parameters, dict):
            raise ClientException("parameters must be a string to string map")

        if self.type == "container":
            if "image" not in self.container_properties:
                raise ClientException("containerProperties must contain image")

            memory = self._get_resource_requirement("memory")
            if memory is None:
                raise ClientException("containerProperties must contain memory")
            if memory < 4:
                raise ClientException("container memory limit must be greater than 4")

            vcpus = self._get_resource_requirement("vcpus")
            if vcpus is None:
                raise ClientException("containerProperties must contain vcpus")
            if vcpus <= 0:
                raise ClientException("container vcpus limit must be greater than 0")

    def deregister(self) -> None:
        self.status = "INACTIVE"

    def update(
        self,
        parameters: Optional[Dict[str, Any]],
        _type: str,
        container_properties: Dict[str, Any],
        node_properties: Dict[str, Any],
        retry_strategy: Dict[str, Any],
        tags: Dict[str, str],
        timeout: Dict[str, int],
    ) -> "JobDefinition":
        if self.status != "INACTIVE":
            if parameters is None:
                parameters = self.parameters

            if _type is None:
                _type = self.type

            if container_properties is None:
                container_properties = self.container_properties

            if retry_strategy is None:
                retry_strategy = self.retry_strategy

        return JobDefinition(
            self.name,
            parameters,
            _type,
            container_properties,
            node_properties=node_properties,
            revision=self.revision,
            retry_strategy=retry_strategy,
            tags=tags,
            timeout=timeout,
            backend=self.backend,
            platform_capabilities=self.platform_capabilities,
            propagate_tags=self.propagate_tags,
        )

    def describe(self) -> Dict[str, Any]:
        result = {
            "jobDefinitionArn": self.arn,
            "jobDefinitionName": self.name,
            "parameters": self.parameters,
            "revision": self.revision,
            "status": self.status,
            "type": self.type,
            "tags": self.backend.tagger.get_tag_dict_for_resource(self.arn),
            "platformCapabilities": self.platform_capabilities,
            "retryStrategy": self.retry_strategy,
            "propagateTags": self.propagate_tags,
        }
        if self.container_properties is not None:
            result["containerProperties"] = self.container_properties
        if self.timeout:
            result["timeout"] = self.timeout

        return result

    @property
    def physical_resource_id(self) -> str:
        return self.arn

    @staticmethod
    def cloudformation_name_type() -> str:
        return "JobDefinitionName"

    @staticmethod
    def cloudformation_type() -> str:
        # https://docs.aws.amazon.com/AWSCloudFormation/latest/UserGuide/aws-resource-batch-jobdefinition.html
        return "AWS::Batch::JobDefinition"

    @classmethod
    def create_from_cloudformation_json(  # type: ignore[misc]
        cls,
        resource_name: str,
        cloudformation_json: Dict[str, Any],
        account_id: str,
        region_name: str,
        **kwargs: Any,
    ) -> "JobDefinition":
        backend = batch_backends[account_id][region_name]
        properties = cloudformation_json["Properties"]
        return backend.register_job_definition(
            def_name=resource_name,
            parameters=lowercase_first_key(properties.get("Parameters", {})),
            _type="container",
            tags=lowercase_first_key(properties.get("Tags", {})),
            retry_strategy=lowercase_first_key(properties["RetryStrategy"]),
            container_properties=(
                lowercase_first_key(properties["ContainerProperties"])  # type: ignore[arg-type]
                if "ContainerProperties" in properties
                else None
            ),
            node_properties=(
                lowercase_first_key(properties["NodeProperties"])  # type: ignore[arg-type]
                if "NodeProperties" in properties
                else None
            ),
            timeout=lowercase_first_key(properties.get("timeout", {})),
            platform_capabilities=None,  # type: ignore[arg-type]
            propagate_tags=None,  # type: ignore[arg-type]
        )


class Job(threading.Thread, BaseModel, DockerModel, ManagedState):
    def __init__(
        self,
        name: str,
        job_def: JobDefinition,
        job_queue: JobQueue,
        log_backend: LogsBackend,
        container_overrides: Optional[Dict[str, Any]],
        depends_on: Optional[List[Dict[str, str]]],
        parameters: Optional[Dict[str, str]],
        all_jobs: Dict[str, "Job"],
        timeout: Optional[Dict[str, int]],
        array_properties: Dict[str, Any],
        provided_job_id: Optional[str] = None,
    ):
        threading.Thread.__init__(self)
        DockerModel.__init__(self)
        ManagedState.__init__(
            self,
            "batch::job",
            JobStatus.status_transitions(),
        )

        self.job_name = name
        self.job_id = provided_job_id or str(mock_random.uuid4())
        self.job_definition = job_def
        self.container_overrides: Dict[str, Any] = container_overrides or {}
        self.job_queue = job_queue
        self.job_queue.jobs.append(self)
        self.job_created_at = datetime.datetime.now()
        self.job_started_at = datetime.datetime(1970, 1, 1)
        self.job_stopped_at = datetime.datetime(1970, 1, 1)
        self.job_stopped = False
        self.job_stopped_reason: Optional[str] = None
        self.depends_on = depends_on
        self.parameters = {**self.job_definition.parameters, **(parameters or {})}
        self.timeout = timeout
        self.all_jobs = all_jobs
        self.array_properties: Dict[str, Any] = array_properties

        self.arn = make_arn_for_job(
            job_def.backend.account_id, self.job_id, job_def._region
        )

        self.stop = False
        self.exit_code: Optional[int] = None

        self.daemon = True

        self.name = "MOTO-BATCH-" + self.job_id

        self._log_backend = log_backend
        self._log_group = "/aws/batch/job"
        self._stream_name = f"{self.job_definition.name}/default/{self.job_id}"
        self.log_stream_name: Optional[str] = None

        self.attempts: List[Dict[str, Any]] = []
        self.latest_attempt: Optional[Dict[str, Any]] = None
        self._child_jobs: Optional[List[Job]] = None

    def describe_short(self) -> Dict[str, Any]:
        result = {
            "jobId": self.job_id,
            "jobArn": self.arn,
            "jobName": self.job_name,
            "createdAt": datetime2int_milliseconds(self.job_created_at),
            "status": self.status,
            "jobDefinition": self.job_definition.arn,
        }
        if self.job_stopped_reason is not None:
            result["statusReason"] = self.job_stopped_reason
        if self.status is not None:
            if JobStatus.is_job_already_started(self.status):
                result["startedAt"] = datetime2int_milliseconds(self.job_started_at)
        if self.job_stopped:
            result["stoppedAt"] = datetime2int_milliseconds(self.job_stopped_at)
            if self.exit_code is not None:
                result["container"] = {"exitCode": self.exit_code}
        return result

    def describe(self) -> Dict[str, Any]:
        result = self.describe_short()
        result["jobQueue"] = self.job_queue.arn
        result["dependsOn"] = self.depends_on or []
        result["parameters"] = {**self.job_definition.parameters, **self.parameters}
        if self.job_definition.type == "container":
            result["container"] = self._container_details()
        elif self.job_definition.type == "multinode":
            result["container"] = {
                "logStreamName": self.log_stream_name,
            }
            result["nodeProperties"] = self.job_definition.node_properties
        if self.job_stopped:
            result["stoppedAt"] = datetime2int_milliseconds(self.job_stopped_at)
        if self.timeout:
            result["timeout"] = self.timeout
        result["attempts"] = self.attempts
        if self._child_jobs:
            child_statuses = {
                "STARTING": 0,
                "FAILED": 0,
                "RUNNING": 0,
                "SUCCEEDED": 0,
                "RUNNABLE": 0,
                "SUBMITTED": 0,
                "PENDING": 0,
            }
            for child_job in self._child_jobs:
                if child_job.status is not None:
                    child_statuses[child_job.status] += 1
            result["arrayProperties"] = {
                "statusSummary": child_statuses,
                "size": len(self._child_jobs),
            }
            if len(self._child_jobs) == child_statuses["SUCCEEDED"]:
                self.status = "SUCCEEDED"
                result["status"] = self.status
        return result

    def _container_details(self) -> Dict[str, Any]:
        details = {}
        details["command"] = self._get_container_property("command", [])
        details["privileged"] = self._get_container_property("privileged", False)
        details["readonlyRootFilesystem"] = self._get_container_property(
            "readonlyRootFilesystem", False
        )
        details["ulimits"] = self._get_container_property("ulimits", {})
        details["vcpus"] = self._get_container_property("vcpus", 1)
        details["memory"] = self._get_container_property("memory", 512)
        details["volumes"] = self._get_container_property("volumes", [])
        details["environment"] = self._get_container_property("environment", [])
        if self.log_stream_name:
            details["logStreamName"] = self.log_stream_name
        return details

    def _get_container_property(self, p: str, default: Any) -> Any:
        if p == "environment":
            job_env = self.container_overrides.get(p, default)
            jd_env = self.job_definition.container_properties.get(p, default)

            job_env_dict = {_env["name"]: _env["value"] for _env in job_env}
            jd_env_dict = {_env["name"]: _env["value"] for _env in jd_env}

            for key in jd_env_dict.keys():
                if key not in job_env_dict.keys():
                    job_env.append({"name": key, "value": jd_env_dict[key]})

            job_env.append({"name": "AWS_BATCH_JOB_ID", "value": self.job_id})

            return job_env

        if p in ["vcpus", "memory"]:
            return self.container_overrides.get(
                p, self.job_definition._get_resource_requirement(p, default)
            )

        return self.container_overrides.get(
            p, self.job_definition.container_properties.get(p, default)
        )

    def _get_attempt_duration(self) -> Optional[int]:
        if self.timeout:
            return self.timeout["attemptDurationSeconds"]
        if self.job_definition.timeout:
            return self.job_definition.timeout["attemptDurationSeconds"]
        return None

    def _add_parameters_to_command(self, command: Union[str, List[str]]) -> List[str]:
        if isinstance(command, str):
            command = [command]

        if not self.parameters:
            return command

        return [
            next(
                (
                    command_part.replace(f"Ref::{param}", value)
                    for param, value in self.parameters.items()
                    if f"Ref::{param}" in command_part
                ),
                command_part,
            )
            for command_part in command
        ]

    def run(self) -> None:
        """
        Run the container.

        Logic is as follows:
        Generate container info (eventually from task definition)
        Start container
        Loop whilst not asked to stop and the container is running.
          Get all logs from container between the last time I checked and now.
        Convert logs into cloudwatch format
        Put logs into cloudwatch

        :return:
        """
        try:
            import docker
        except ImportError as err:
            logger.error(f"Failed to run AWS Batch container {self.name}. Error {err}")
            self._mark_stopped(success=False)
            return

        try:
            containers: List[docker.models.containers.Container] = []

            self.advance()
            while self.status == JobStatus.SUBMITTED:
                # Wait until we've moved onto state 'PENDING'
                sleep(0.5)

            # Wait until all dependent jobs have finished
            # If any of the dependent jobs have failed, not even start
            if self.depends_on and not self._wait_for_dependencies():
                return

            container_kwargs = []
            if self.job_definition.container_properties:
                volumes = {
                    v["name"]: v["host"]
                    for v in self._get_container_property("volumes", [])
                }
                container_kwargs.append(
                    {
                        "image": self.job_definition.container_properties.get(
                            "image", "alpine:latest"
                        ),
                        "privileged": self.job_definition.container_properties.get(
                            "privileged", False
                        ),
                        "command": self._add_parameters_to_command(
                            self._get_container_property(
                                "command",
                                '/bin/sh -c "for a in `seq 1 10`; do echo Hello World; sleep 1; done"',
                            )
                        ),
                        "environment": {
                            e["name"]: e["value"]
                            for e in self._get_container_property("environment", [])
                        },
                        "mounts": [
                            docker.types.Mount(
                                m["containerPath"],
                                volumes[m["sourceVolume"]]["sourcePath"],
                                type="bind",
                                read_only=m["readOnly"],
                            )
                            for m in self._get_container_property("mountPoints", [])
                        ],
                        "name": f"{self.job_name}-{self.job_id.replace(':', '-')}",
                    }
                )
            else:
                node_properties = self.job_definition.node_properties
                num_nodes = node_properties["numNodes"]
                node_containers = {}
                for node_range in node_properties["nodeRangeProperties"]:
                    start, sep, end = node_range["targetNodes"].partition(":")
                    if sep == "":
                        start = end = int(start)
                    else:
                        if start == "":
                            start = 0
                        else:
                            start = int(start)
                        if end == "":
                            end = num_nodes - 1
                        else:
                            end = int(end)
                    for i in range(start, end + 1):
                        node_containers[i] = node_range["container"]

                for i in range(num_nodes):
                    spec = node_containers[i]
                    volumes = {v["name"]: v["host"] for v in spec.get("volumes", [])}
                    container_kwargs.append(
                        {
                            "image": spec.get("image", "alpine:latest"),
                            "privileged": spec.get("privileged", False),
                            "command": self._add_parameters_to_command(
                                spec.get(
                                    "command",
                                    '/bin/sh -c "for a in `seq 1 10`; do echo Hello World; sleep 1; done"',
                                )
                            ),
                            "environment": {
                                e["name"]: e["value"]
                                for e in spec.get("environment", [])
                            },
                            "mounts": [
                                docker.types.Mount(
                                    m["containerPath"],
                                    volumes[m["sourceVolume"]]["sourcePath"],
                                    type="bind",
                                    read_only=m["readOnly"],
                                )
                                for m in spec.get("mountPoints", [])
                            ],
                            "name": f"{self.job_name}-{self.job_id}-{i}",
                        }
                    )

            self.advance()
            while self.status == JobStatus.PENDING:
                # Wait until the state is no longer pending, but 'RUNNABLE'
                sleep(0.5)
            # TODO setup ecs container instance

            self.job_started_at = datetime.datetime.now()

            self._start_attempt()

            # add host.docker.internal host on linux to emulate Mac + Windows behavior
            #   for communication with other mock AWS services running on localhost
            extra_hosts = (
                {"host.docker.internal": "host-gateway"}
                if platform == "linux" or platform == "linux2"
                else {}
            )

            network_mode = settings.moto_network_mode()
            network_name = settings.moto_network_name()

            for kwargs in container_kwargs:
                environment = kwargs["environment"]
                environment["MOTO_HOST"] = settings.moto_server_host()
                environment["MOTO_PORT"] = settings.moto_server_port()
                environment["MOTO_HTTP_ENDPOINT"] = (
                    f'{environment["MOTO_HOST"]}:{environment["MOTO_PORT"]}'
                )

                if network_name:
                    kwargs["network"] = network_name
                elif network_mode:
                    kwargs["network_mode"] = network_mode

            log_config = docker.types.LogConfig(type=docker.types.LogConfig.types.JSON)
            self.advance()
            while self.status == JobStatus.RUNNABLE:
                # Wait until the state is no longer runnable, but 'STARTING'
                sleep(0.5)

            self.advance()
            while self.status == JobStatus.STARTING:
                # Wait until the state is no longer runnable, but 'RUNNING'
                sleep(0.5)

            for kwargs in container_kwargs:
                if len(containers) > 0:
                    env = kwargs["environment"]
                    network_settings = containers[0].attrs["NetworkSettings"]
                    networks = network_settings["Networks"]
                    if network_name in networks:
                        ip = networks[network_name]["IPAddress"]
                    else:
                        ip = network_settings["IPAddress"]
                    env["AWS_BATCH_JOB_MAIN_NODE_PRIVATE_IPV4_ADDRESS"] = ip
                container = self.docker_client.containers.run(
                    detach=True,
                    log_config=log_config,
                    extra_hosts=extra_hosts,
                    **kwargs,
                )
                container.reload()
                containers.append(container)

            for i, container in enumerate(containers):
                try:
                    container.reload()

                    max_time = None
                    if self._get_attempt_duration():
                        attempt_duration = self._get_attempt_duration()
                        max_time = self.job_started_at + datetime.timedelta(
                            seconds=attempt_duration  # type: ignore[arg-type]
                        )

                    while container.status == "running" and not self.stop:
                        container.reload()
                        time.sleep(0.5)

                        if max_time and datetime.datetime.now() > max_time:
                            raise Exception(
                                "Job time exceeded the configured attemptDurationSeconds"
                            )

                    # Container should be stopped by this point... unless asked to stop
                    if container.status == "running":
                        container.kill()

                    # Log collection
                    logs_stdout = []
                    logs_stderr = []
                    logs_stderr.extend(
                        container.logs(
                            stdout=False,
                            stderr=True,
                            timestamps=True,
                            since=datetime2int(self.job_started_at),
                        )
                        .decode()
                        .split("\n")
                    )
                    logs_stdout.extend(
                        container.logs(
                            stdout=True,
                            stderr=False,
                            timestamps=True,
                            since=datetime2int(self.job_started_at),
                        )
                        .decode()
                        .split("\n")
                    )

                    # Process logs
                    logs_stdout = [x for x in logs_stdout if len(x) > 0]
                    logs_stderr = [x for x in logs_stderr if len(x) > 0]
                    logs = []
                    for line in logs_stdout + logs_stderr:
                        date, line = line.split(" ", 1)
                        date_obj = (
                            dateutil.parser.parse(date)
                            .astimezone(datetime.timezone.utc)
                            .replace(tzinfo=None)
                        )
                        date = unix_time_millis(date_obj)
                        logs.append({"timestamp": date, "message": line.strip()})
                    logs = sorted(logs, key=lambda log: log["timestamp"])

                    # Send to cloudwatch
                    self.log_stream_name = self._stream_name
                    self._log_backend.ensure_log_group(self._log_group)
                    self._log_backend.ensure_log_stream(
                        self._log_group, self.log_stream_name
                    )
                    self._log_backend.put_log_events(
                        self._log_group, self.log_stream_name, logs
                    )

                    result = container.wait() or {}
                    exit_code = result.get("StatusCode", 0)
                    self.exit_code = exit_code
                    job_failed = self.stop or exit_code > 0
                    if job_failed:
                        self._mark_stopped(success=False)
                        break

                except Exception as err:
                    logger.error(
                        f"Failed to run AWS Batch container {self.name}. Error {err}"
                    )
                    self._mark_stopped(success=False)

            self._mark_stopped(success=True)
        except Exception as err:
            logger.error(f"Failed to run AWS Batch container {self.name}. Error {err}")
            self._mark_stopped(success=False)
        finally:
            for container in containers:
                container.reload()
                if container.status == "running":
                    container.kill()
                container.remove()

    def _mark_stopped(self, success: bool = True) -> None:
        if self.job_stopped:
            return
        # Ensure that job_stopped/job_stopped_at-attributes are set first
        # The describe-method needs them immediately when status is set
        self.job_stopped = True
        self.job_stopped_at = datetime.datetime.now()
        self.status = JobStatus.SUCCEEDED if success else JobStatus.FAILED
        self._stop_attempt()

    def _start_attempt(self) -> None:
        self.latest_attempt = {
            "container": {
                "containerInstanceArn": "TBD",
                "logStreamName": self.log_stream_name,
                "networkInterfaces": [],
                "taskArn": self.job_definition.arn,
            }
        }
        self.latest_attempt["startedAt"] = datetime2int_milliseconds(
            self.job_started_at
        )
        self.attempts.append(self.latest_attempt)

    def _stop_attempt(self) -> None:
        if self.latest_attempt:
            self.latest_attempt["container"]["logStreamName"] = self.log_stream_name
            self.latest_attempt["stoppedAt"] = datetime2int_milliseconds(
                self.job_stopped_at
            )

    def terminate(self, reason: str) -> None:
        if not self.stop:
            self.stop = True
            self.job_stopped_reason = reason

    def _wait_for_dependencies(self) -> bool:
        dependent_ids = [dependency["jobId"] for dependency in self.depends_on]  # type: ignore[union-attr]
        successful_dependencies: Set[str] = set()
        while len(successful_dependencies) != len(dependent_ids):
            for dependent_id in dependent_ids:
                if dependent_id in self.all_jobs:
                    dependent_job = self.all_jobs[dependent_id]
                    if dependent_job.status == JobStatus.SUCCEEDED:
                        successful_dependencies.add(dependent_id)
                    if dependent_job.status == JobStatus.FAILED:
                        logger.error(
                            f"Terminating job {self.name} due to failed dependency {dependent_job.name}"
                        )
                        self._mark_stopped(success=False)
                        return False

            time.sleep(1)
            if self.stop:
                # This job has been cancelled while it was waiting for a dependency
                self._mark_stopped(success=False)
                return False

        return True


class SchedulingPolicy(BaseModel):
    def __init__(
        self,
        account_id: str,
        region: str,
        name: str,
        fairshare_policy: Dict[str, Any],
        backend: "BatchBackend",
        tags: Dict[str, str],
    ):
        self.name = name
        self.arn = f"arn:{get_partition(region)}:batch:{region}:{account_id}:scheduling-policy/{name}"
        self.fairshare_policy = {
            "computeReservation": fairshare_policy.get("computeReservation") or 0,
            "shareDecaySeconds": fairshare_policy.get("shareDecaySeconds") or 0,
            "shareDistribution": fairshare_policy.get("shareDistribution") or [],
        }
        self.backend = backend
        if tags:
            backend.tag_resource(self.arn, tags)

    def to_dict(self, create: bool = False) -> Dict[str, Any]:
        resp: Dict[str, Any] = {"name": self.name, "arn": self.arn}
        if not create:
            resp["fairsharePolicy"] = self.fairshare_policy
            resp["tags"] = self.backend.list_tags_for_resource(self.arn)
        return resp


class BatchBackend(BaseBackend):
    """
    Batch-jobs are executed inside a Docker-container. Everytime the `submit_job`-method is called, a new Docker container is started.
    A job is marked as 'Success' when the Docker-container exits without throwing an error.

    Use `@mock_batch_simple` instead if you do not want to use a Docker-container.
    With this decorator, jobs are simply marked as 'Success' without trying to execute any commands/scripts.
    """

    def __init__(self, region_name: str, account_id: str):
        super().__init__(region_name, account_id)
        self.tagger = TaggingService()

        self._compute_environments: Dict[str, ComputeEnvironment] = {}
        self._job_queues: Dict[str, JobQueue] = {}
        self._job_definitions: Dict[str, JobDefinition] = {}
        self._jobs: Dict[str, Job] = {}
        self._scheduling_policies: Dict[str, SchedulingPolicy] = {}

    @property
    def iam_backend(self) -> IAMBackend:
        """
        :return: IAM Backend
        :rtype: moto.iam.models.IAMBackend
        """
        return iam_backends[self.account_id][self.partition]

    @property
    def ec2_backend(self) -> EC2Backend:
        """
        :return: EC2 Backend
        :rtype: moto.ec2.models.EC2Backend
        """
        return ec2_backends[self.account_id][self.region_name]

    @property
    def ecs_backend(self) -> EC2ContainerServiceBackend:
        """
        :return: ECS Backend
        :rtype: moto.ecs.models.EC2ContainerServiceBackend
        """
        return ecs_backends[self.account_id][self.region_name]

    @property
    def logs_backend(self) -> LogsBackend:
        """
        :return: ECS Backend
        :rtype: moto.logs.models.LogsBackend
        """
        return logs_backends[self.account_id][self.region_name]

    def reset(self) -> None:
        for job in self._jobs.values():
            if job.status not in (JobStatus.FAILED, JobStatus.SUCCEEDED):
                job.stop = True
                # Try to join
                if job.is_alive():
                    job.join(0.2)

        super().reset()

    def get_compute_environment_by_arn(self, arn: str) -> Optional[ComputeEnvironment]:
        return self._compute_environments.get(arn)

    def get_compute_environment_by_name(
        self, name: str
    ) -> Optional[ComputeEnvironment]:
        for comp_env in self._compute_environments.values():
            if comp_env.name == name:
                return comp_env
        return None

    def get_compute_environment(self, identifier: str) -> Optional[ComputeEnvironment]:
        """
        Get compute environment by name or ARN
        :param identifier: Name or ARN
        :type identifier: str

        :return: Compute Environment or None
        :rtype: ComputeEnvironment or None
        """
        return self.get_compute_environment_by_arn(
            identifier
        ) or self.get_compute_environment_by_name(identifier)

    def get_job_queue_by_arn(self, arn: str) -> Optional[JobQueue]:
        return self._job_queues.get(arn)

    def get_job_queue_by_name(self, name: str) -> Optional[JobQueue]:
        for comp_env in self._job_queues.values():
            if comp_env.name == name:
                return comp_env
        return None

    def get_job_queue(self, identifier: str) -> Optional[JobQueue]:
        """
        Get job queue by name or ARN
        :param identifier: Name or ARN
        :type identifier: str

        :return: Job Queue or None
        :rtype: JobQueue or None
        """
        return self.get_job_queue_by_arn(identifier) or self.get_job_queue_by_name(
            identifier
        )

    def get_job_definition_by_arn(self, arn: str) -> Optional[JobDefinition]:
        return self._job_definitions.get(arn)

    def get_job_definition_by_name(self, name: str) -> Optional[JobDefinition]:
        latest_revision = -1
        latest_job = None
        for job_def in self._job_definitions.values():
            if job_def.name == name and job_def.revision > latest_revision:
                latest_job = job_def
                latest_revision = job_def.revision
        return latest_job

    def get_job_definition_by_name_revision(
        self, name: str, revision: str
    ) -> Optional[JobDefinition]:
        for job_def in self._job_definitions.values():
            if job_def.name == name and job_def.revision == int(revision):
                return job_def
        return None

    def get_job_definition(self, identifier: str) -> Optional[JobDefinition]:
        """
        Get job definitions by name or ARN
        :param identifier: Name or ARN
        :type identifier: str

        :return: Job definition or None
        :rtype: JobDefinition or None
        """
        job_def = self.get_job_definition_by_arn(identifier)
        if job_def is None:
            if ":" in identifier:
                job_def = self.get_job_definition_by_name_revision(
                    *identifier.split(":", 1)
                )
            else:
                job_def = self.get_job_definition_by_name(identifier)
        return job_def

    def get_job_definitions(self, identifier: str) -> List[JobDefinition]:
        """
        Get job definitions by name or ARN
        :param identifier: Name or ARN
        :type identifier: str

        :return: Job definition or None
        :rtype: list of JobDefinition
        """
        result = []
        env = self.get_job_definition_by_arn(identifier)
        if env is not None:
            result.append(env)
        else:
            for value in self._job_definitions.values():
                if value.name == identifier:
                    result.append(value)

        return result

    def get_job_by_id(self, identifier: str) -> Optional[Job]:
        try:
            return self._jobs[identifier]
        except KeyError:
            return None

    def describe_compute_environments(
        self, environments: Optional[List[str]] = None
    ) -> List[Dict[str, Any]]:
        """
        Pagination is not yet implemented
        """
        envs = set()
        if environments is not None:
            envs = set(environments)

        result = []
        for arn, environment in self._compute_environments.items():
            # Filter shortcut
            if len(envs) > 0 and arn not in envs and environment.name not in envs:
                continue

            json_part: Dict[str, Any] = {
                "computeEnvironmentArn": arn,
                "computeEnvironmentName": environment.name,
                "ecsClusterArn": environment.ecs_arn,
                "serviceRole": environment.service_role,
                "state": environment.state,
                "type": environment.env_type,
                "status": "VALID",
                "statusReason": "Compute environment is available",
            }
            if environment.env_type == "MANAGED":
                json_part["computeResources"] = environment.compute_resources

            result.append(json_part)

        return result

    def create_compute_environment(
        self,
        compute_environment_name: str,
        _type: str,
        state: str,
        compute_resources: Dict[str, Any],
        service_role: str,
    ) -> ComputeEnvironment:
        # Validate
        if COMPUTE_ENVIRONMENT_NAME_REGEX.match(compute_environment_name) is None:
            raise InvalidParameterValueException(
                "Compute environment name does not match ^[A-Za-z0-9][A-Za-z0-9_-]{1,126}[A-Za-z0-9]$"
            )

        if self.get_compute_environment_by_name(compute_environment_name) is not None:
            raise InvalidParameterValueException(
                f"A compute environment already exists with the name {compute_environment_name}"
            )

        if not service_role:
            raise ClientException(
                f"Error executing request, Exception : ServiceRole is required.,"
                f" RequestId: {mock_random.uuid4()}"
            )

        # Look for IAM role
        try:
            self.iam_backend.get_role_by_arn(service_role)
        except IAMNotFoundException:
            raise InvalidParameterValueException(
                f"Could not find IAM role {service_role}"
            )

        if _type not in ("MANAGED", "UNMANAGED"):
            raise InvalidParameterValueException(
                f"type {_type} must be one of MANAGED | UNMANAGED"
            )

        if state is not None and state not in ("ENABLED", "DISABLED"):
            raise InvalidParameterValueException(
                f"state {state} must be one of ENABLED | DISABLED"
            )

        if compute_resources is None and _type == "MANAGED":
            raise InvalidParameterValueException(
                f"computeResources must be specified when creating a {state} environment"
            )
        elif compute_resources is not None:
            self._validate_compute_resources(compute_resources)

        # By here, all values except SPOT ones have been validated
        new_comp_env = ComputeEnvironment(
            compute_environment_name,
            _type,
            state,
            compute_resources,
            service_role,
            account_id=self.account_id,
            region_name=self.region_name,
        )
        self._compute_environments[new_comp_env.arn] = new_comp_env

        # Ok by this point, everything is legit, so if its Managed then start some instances
        if _type == "MANAGED" and "FARGATE" not in compute_resources["type"]:
            cpus = int(
                compute_resources.get("desiredvCpus", compute_resources["minvCpus"])
            )
            instance_types = compute_resources["instanceTypes"]
            needed_instance_types = self.find_min_instances_to_meet_vcpus(
                instance_types, cpus
            )
            # Create instances

            # Will loop over and over so we get decent subnet coverage
            subnet_cycle = cycle(compute_resources["subnets"])

            for instance_type in needed_instance_types:
                reservation = self.ec2_backend.run_instances(
                    image_id="ami-03cf127a",  # Todo import AMIs
                    count=1,
                    user_data=None,
                    security_group_names=[],
                    instance_type=instance_type,
                    region_name=self.region_name,
                    subnet_id=next(subnet_cycle),
                    key_name=compute_resources.get("ec2KeyPair", "AWS_OWNED"),
                    security_group_ids=compute_resources["securityGroupIds"],
                    is_instance_type_default=False,
                )

                new_comp_env.add_instance(reservation.instances[0])

        # Create ECS cluster
        # Should be of format P2OnDemand_Batch_UUID
        cluster_name = "OnDemand_Batch_" + str(mock_random.uuid4())
        ecs_cluster = self.ecs_backend.create_cluster(cluster_name)
        new_comp_env.set_ecs(ecs_cluster.arn, cluster_name)

        return new_comp_env

    def _validate_compute_resources(self, cr: Dict[str, Any]) -> None:
        """
        Checks contents of sub dictionary for managed clusters

        :param cr: computeResources
        :type cr: dict
        """
        if int(cr["maxvCpus"]) < 0:
            raise InvalidParameterValueException("maxVCpus must be positive")
        if "FARGATE" not in cr["type"]:
            # Most parameters are not applicable to jobs that are running on Fargate resources:
            # non exhaustive list: minvCpus, instanceTypes, imageId, ec2KeyPair, instanceRole, tags
            if "instanceRole" not in cr:
                raise ClientException(
                    "Error executing request, Exception : Instance role is required."
                )
            for profile in self.iam_backend.get_instance_profiles():
                if profile.arn == cr["instanceRole"]:
                    break
            else:
                raise InvalidParameterValueException(
                    f"could not find instanceRole {cr['instanceRole']}"
                )

            if "minvCpus" not in cr:
                raise ClientException(
                    "Error executing request, Exception : Resource minvCpus is required."
                )
            if int(cr["minvCpus"]) < 0:
                raise InvalidParameterValueException("minvCpus must be positive")
            if int(cr["maxvCpus"]) < int(cr["minvCpus"]):
                raise InvalidParameterValueException(
                    "maxVCpus must be greater than minvCpus"
                )

            if len(cr["instanceTypes"]) == 0:
                raise InvalidParameterValueException(
                    "At least 1 instance type must be provided"
                )
            for instance_type in cr["instanceTypes"]:
                if instance_type == "optimal":
                    pass  # Optimal should pick from latest of current gen
                elif (
                    instance_type not in EC2_INSTANCE_TYPES
                    and instance_type not in EC2_INSTANCE_FAMILIES
                ):
                    raise InvalidParameterValueException(
                        f"Instance type {instance_type} does not exist"
                    )

        for sec_id in cr["securityGroupIds"]:
            if self.ec2_backend.get_security_group_from_id(sec_id) is None:
                raise InvalidParameterValueException(
                    f"security group {sec_id} does not exist"
                )
        if len(cr["securityGroupIds"]) == 0:
            raise InvalidParameterValueException(
                "At least 1 security group must be provided"
            )

        for subnet_id in cr["subnets"]:
            try:
                self.ec2_backend.get_subnet(subnet_id)
            except InvalidSubnetIdError:
                raise InvalidParameterValueException(
                    f"subnet {subnet_id} does not exist"
                )
        if len(cr["subnets"]) == 0:
            raise InvalidParameterValueException("At least 1 subnet must be provided")

        if cr["type"] not in {"EC2", "SPOT", "FARGATE", "FARGATE_SPOT"}:
            raise InvalidParameterValueException(
                "computeResources.type must be either EC2 | SPOT | FARGATE | FARGATE_SPOT"
            )

    @staticmethod
    def find_min_instances_to_meet_vcpus(
        instance_types: List[str], target: float
    ) -> List[str]:
        """
        Finds the minimum needed instances to meed a vcpu target

        :param instance_types: Instance types, like ['t2.medium', 't2.small']
        :param target: VCPU target
        :return: List of instance types
        """
        # vcpus = [ (vcpus, instance_type), (vcpus, instance_type), ... ]
        instance_vcpus = []
        instances = []

        for instance_type in instance_types:
            if instance_type == "optimal":
                instance_type = "m4.4xlarge"

            if "." not in instance_type:
                # instance_type can be a family of instance types (c2, t3, etc)
                # We'll just use the first instance_type in this family
                instance_type = [
                    i for i in EC2_INSTANCE_TYPES.keys() if i.startswith(instance_type)
                ][0]
            instance_vcpus.append(
                (
                    EC2_INSTANCE_TYPES[instance_type]["VCpuInfo"]["DefaultVCpus"],
                    instance_type,
                )
            )

        instance_vcpus = sorted(instance_vcpus, key=lambda item: item[0], reverse=True)
        # Loop through,
        #   if biggest instance type smaller than target, and len(instance_types)> 1, then use biggest type
        #   if biggest instance type bigger than target, and len(instance_types)> 1, then remove it and move on

        #   if biggest instance type bigger than target and len(instan_types) == 1 then add instance and finish
        #   if biggest instance type smaller than target and len(instan_types) == 1 then loop adding instances until target == 0
        #   ^^ boils down to keep adding last till target vcpus is negative
        #   #Algorithm ;-) ... Could probably be done better with some quality lambdas
        while target > 0:
            current_vcpu, current_instance = instance_vcpus[0]

            if len(instance_vcpus) > 1:
                if current_vcpu <= target:
                    target -= current_vcpu
                    instances.append(current_instance)
                else:
                    # try next biggest instance
                    instance_vcpus.pop(0)
            else:
                # Were on the last instance
                target -= current_vcpu
                instances.append(current_instance)

        return instances

    def delete_compute_environment(self, compute_environment_name: str) -> None:
        if compute_environment_name is None:
            raise InvalidParameterValueException("Missing computeEnvironment parameter")

        compute_env = self.get_compute_environment(compute_environment_name)

        if compute_env is not None:
            # Pop ComputeEnvironment
            self._compute_environments.pop(compute_env.arn)

            # Delete ECS cluster
            self.ecs_backend.delete_cluster(compute_env.ecs_name)

            if compute_env.env_type == "MANAGED":
                # Delete compute environment
                instance_ids = [instance.id for instance in compute_env.instances]
                if instance_ids:
                    self.ec2_backend.terminate_instances(instance_ids)

    def update_compute_environment(
        self,
        compute_environment_name: str,
        state: Optional[str],
        compute_resources: Optional[Any],
        service_role: Optional[str],
    ) -> Tuple[str, str]:
        # Validate
        compute_env = self.get_compute_environment(compute_environment_name)
        if compute_env is None:
            raise ClientException("Compute environment {0} does not exist")

        # Look for IAM role
        if service_role is not None:
            try:
                self.iam_backend.get_role_by_arn(service_role)
            except IAMNotFoundException:
                raise InvalidParameterValueException(
                    f"Could not find IAM role {service_role}"
                )

            compute_env.service_role = service_role

        if state is not None:
            if state not in ("ENABLED", "DISABLED"):
                raise InvalidParameterValueException(
                    f"state {state} must be one of ENABLED | DISABLED"
                )

            compute_env.state = state

        if compute_resources is not None:
            # TODO Implement resizing of instances based on changing vCpus
            # compute_resources CAN contain desiredvCpus, maxvCpus, minvCpus, and can contain none of them.
            pass

        return compute_env.name, compute_env.arn

    def create_job_queue(
        self,
        queue_name: str,
        priority: str,
        schedule_policy: Optional[str],
        state: str,
        compute_env_order: List[Dict[str, str]],
        tags: Optional[Dict[str, str]] = None,
    ) -> JobQueue:
        for variable, var_name in (
            (queue_name, "jobQueueName"),
            (priority, "priority"),
            (state, "state"),
            (compute_env_order, "computeEnvironmentOrder"),
        ):
            if variable is None:
                raise ClientException(f"{var_name} must be provided")

        if state not in ("ENABLED", "DISABLED"):
            raise ClientException(f"state {state} must be one of ENABLED | DISABLED")
        if self.get_job_queue_by_name(queue_name) is not None:
            raise ClientException(f"Job queue {queue_name} already exists")

        if len(compute_env_order) == 0:
            raise ClientException("At least 1 compute environment must be provided")
        try:
            # orders and extracts computeEnvironment names
            ordered_compute_environments = [
                item["computeEnvironment"]
                for item in sorted(compute_env_order, key=lambda x: x["order"])
            ]
            env_objects = []
            # Check each ARN exists, then make a list of compute env's
            for arn in ordered_compute_environments:
                env = self.get_compute_environment_by_arn(arn)
                if env is None:
                    raise ClientException(f"Compute environment {arn} does not exist")
                env_objects.append(env)
        except Exception:
            raise ClientException("computeEnvironmentOrder is malformed")

        # Create new Job Queue
        queue = JobQueue(
            queue_name,
            priority,
            state,
            env_objects,
            compute_env_order,
            schedule_policy=schedule_policy,
            backend=self,
            tags=tags,
        )
        self._job_queues[queue.arn] = queue

        return queue

    def describe_job_queues(
        self, job_queues: Optional[List[str]] = None
    ) -> List[Dict[str, Any]]:
        """
        Pagination is not yet implemented
        """
        envs = set()
        if job_queues is not None:
            envs = set(job_queues)

        result = []
        for arn, job_queue in self._job_queues.items():
            # Filter shortcut
            if len(envs) > 0 and arn not in envs and job_queue.name not in envs:
                continue

            result.append(job_queue.describe())

        return result

    def update_job_queue(
        self,
        queue_name: str,
        priority: Optional[str],
        state: Optional[str],
        compute_env_order: Optional[List[Dict[str, Any]]],
        schedule_policy: Optional[str],
    ) -> Tuple[str, str]:
        if queue_name is None:
            raise ClientException("jobQueueName must be provided")

        job_queue = self.get_job_queue(queue_name)
        if job_queue is None:
            raise ClientException(f"Job queue {queue_name} does not exist")

        if state is not None:
            if state not in ("ENABLED", "DISABLED"):
                raise ClientException(
                    f"state {state} must be one of ENABLED | DISABLED"
                )

            job_queue.state = state

        if compute_env_order is not None:
            if len(compute_env_order) == 0:
                raise ClientException("At least 1 compute environment must be provided")
            try:
                # orders and extracts computeEnvironment names
                ordered_compute_environments = [
                    item["computeEnvironment"]
                    for item in sorted(compute_env_order, key=lambda x: x["order"])
                ]
                env_objects = []
                # Check each ARN exists, then make a list of compute env's
                for arn in ordered_compute_environments:
                    env = self.get_compute_environment_by_arn(arn)
                    if env is None:
                        raise ClientException(
                            f"Compute environment {arn} does not exist"
                        )
                    env_objects.append(env)
            except Exception:
                raise ClientException("computeEnvironmentOrder is malformed")

            job_queue.env_order_json = compute_env_order
            job_queue.environments = env_objects

        if priority is not None:
            job_queue.priority = priority
        if schedule_policy is not None:
            job_queue.schedule_policy = schedule_policy

        return queue_name, job_queue.arn

    def delete_job_queue(self, queue_name: str) -> None:
        job_queue = self.get_job_queue(queue_name)

        if job_queue is not None:
            del self._job_queues[job_queue.arn]

    def register_job_definition(
        self,
        def_name: str,
        parameters: Dict[str, Any],
        _type: str,
        tags: Dict[str, str],
        retry_strategy: Dict[str, Any],
        container_properties: Dict[str, Any],
        node_properties: Dict[str, Any],
        timeout: Dict[str, int],
        platform_capabilities: List[str],
        propagate_tags: bool,
    ) -> JobDefinition:
        if def_name is None:
            raise ClientException("jobDefinitionName must be provided")

        job_def = self.get_job_definition_by_name(def_name)
        if retry_strategy is not None and "evaluateOnExit" in retry_strategy:
            for strat in retry_strategy["evaluateOnExit"]:
                if "action" in strat:
                    strat["action"] = strat["action"].lower()
        if not tags:
            tags = {}
        if job_def is None:
            job_def = JobDefinition(
                def_name,
                parameters,
                _type,
                container_properties,
                node_properties=node_properties,
                tags=tags,
                retry_strategy=retry_strategy,
                timeout=timeout,
                backend=self,
                platform_capabilities=platform_capabilities,
                propagate_tags=propagate_tags,
            )
        else:
            # Make new jobdef
            job_def = job_def.update(
                parameters,
                _type,
                container_properties,
                node_properties,
                retry_strategy,
                tags,
                timeout,
            )

        self._job_definitions[job_def.arn] = job_def

        return job_def

    def deregister_job_definition(self, def_name: str) -> None:
        job_def = self.get_job_definition_by_arn(def_name)
        if job_def is None and ":" in def_name:
            name, revision = def_name.split(":", 1)
            job_def = self.get_job_definition_by_name_revision(name, revision)

        if job_def is not None:
            self._job_definitions[job_def.arn].deregister()

    def describe_job_definitions(
        self,
        job_def_name: Optional[str] = None,
        job_def_list: Optional[List[str]] = None,
        status: Optional[str] = None,
    ) -> List[JobDefinition]:
        """
        Pagination is not yet implemented
        """
        jobs = []

        # As a job name can reference multiple revisions, we get a list of them
        if job_def_name is not None:
            job_def = self.get_job_definitions(job_def_name)
            if job_def is not None:
                jobs.extend(job_def)
        elif job_def_list is not None:
            for jdn in job_def_list:
                job_def = self.get_job_definitions(jdn)
                if job_def is not None:
                    jobs.extend(job_def)
        else:
            jobs.extend(self._job_definitions.values())

        # Got all the job defs were after, filter then by status
        if status is not None:
            return [job for job in jobs if job.status == status]
        for job in jobs:
            job.describe()
        return jobs

    def submit_job(
        self,
        job_name: str,
        job_def_id: str,
        job_queue: str,
        array_properties: Dict[str, int],
        depends_on: Optional[List[Dict[str, str]]] = None,
        container_overrides: Optional[Dict[str, Any]] = None,
        timeout: Optional[Dict[str, int]] = None,
        parameters: Optional[Dict[str, str]] = None,
    ) -> Tuple[str, str, str]:
        """
        Parameters RetryStrategy and Parameters are not yet implemented.
        """
        if JOB_NAME_REGEX.match(job_name) is None:
            raise ClientException("Job name should match valid pattern")

        # Look for job definition
        job_def = self.get_job_definition(job_def_id)
        if job_def is None:
            raise ClientException(f"Job definition {job_def_id} does not exist")

        queue = self.get_job_queue(job_queue)
        if queue is None:
            raise ClientException(f"Job queue {job_queue} does not exist")

        job = Job(
            job_name,
            job_def,
            queue,
            log_backend=self.logs_backend,
            container_overrides=container_overrides,
            depends_on=depends_on,
            all_jobs=self._jobs,
            timeout=timeout,
            array_properties=array_properties or {},
            parameters=parameters,
        )
        self._jobs[job.job_id] = job

        if "size" in array_properties:
            child_jobs = []
            for array_index in range(array_properties["size"]):
                provided_job_id = f"{job.job_id}:{array_index}"
                child_job = Job(
                    job_name,
                    job_def,
                    queue,
                    log_backend=self.logs_backend,
                    container_overrides=container_overrides,
                    depends_on=depends_on,
                    all_jobs=self._jobs,
                    timeout=timeout,
                    array_properties={"statusSummary": {}, "index": array_index},
                    parameters=parameters,
                    provided_job_id=provided_job_id,
                )
                child_jobs.append(child_job)
                self._jobs[child_job.job_id] = child_job
                child_job.start()

            # The 'parent' job doesn't need to be executed
            # it just needs to keep track of it's children
            job._child_jobs = child_jobs
        else:
            # Here comes the fun
            job.start()
        return job_name, job.job_id, job.arn

    def describe_jobs(self, jobs: Optional[List[str]]) -> List[Dict[str, Any]]:
        job_filter = set()
        if jobs is not None:
            job_filter = set(jobs)

        result = []
        for key, job in self._jobs.items():
            if (
                len(job_filter) > 0
                and key not in job_filter
                and job.arn not in job_filter
            ):
                continue

            result.append(job.describe())

        return result

    def list_jobs(
        self,
        job_queue_name: Optional[str],
        array_job_id: Optional[str],
        job_status: Optional[str] = None,
        filters: Optional[List[Dict[str, Any]]] = None,
    ) -> List[Job]:
        """
        Pagination is not yet implemented
        """
        jobs_to_check = []
        jobs = []

        if job_queue_name:
            if job_queue := self.get_job_queue(job_queue_name):
                jobs_to_check.extend(job_queue.jobs)
            else:
                raise ClientException(f"Job queue {job_queue_name} does not exist")
        if array_job_id:
            if array_job := self.get_job_by_id(array_job_id):
                jobs_to_check.extend(array_job._child_jobs or [])

        if job_status is not None and job_status not in JobStatus.job_statuses():
            raise ClientException(
                "Job status is not one of SUBMITTED | PENDING | RUNNABLE | STARTING | RUNNING | SUCCEEDED | FAILED"
            )

        for job in jobs_to_check:
            if job_status is not None and job.status != job_status:
                continue

            if filters is not None:
                matches = True
                for filt in filters:
                    name = filt["name"]
                    values = filt["values"]
                    if name == "JOB_NAME":
                        if job.job_name not in values:
                            matches = False
                            break
                if not matches:
                    continue

            jobs.append(job)

        return jobs

    def cancel_job(self, job_id: str, reason: str) -> None:
        if job_id == "":
            raise ClientException(
                "'jobId' is a required field (cannot be an empty string)"
            )
        if reason == "":
            raise ClientException(
                "'reason' is a required field (cannot be an empty string)"
            )

        job = self.get_job_by_id(job_id)
        if job is not None:
            if job.status is None:
                return
            if JobStatus.is_job_before_starting(job.status):
                job.terminate(reason)
            # No-Op for jobs that have already started - user has to explicitly terminate those

    def terminate_job(self, job_id: str, reason: str) -> None:
        if job_id == "":
            raise ClientException(
                "'jobId' is a required field (cannot be a empty string)"
            )
        if reason == "":
            raise ClientException(
                "'reason' is a required field (cannot be a empty string)"
            )

        job = self.get_job_by_id(job_id)
        if job is not None:
            job.terminate(reason)

    def tag_resource(self, resource_arn: str, tags: Dict[str, str]) -> None:
        tag_list = self.tagger.convert_dict_to_tags_input(tags or {})
        self.tagger.tag_resource(resource_arn, tag_list)

    def list_tags_for_resource(self, resource_arn: str) -> Dict[str, str]:
        return self.tagger.get_tag_dict_for_resource(resource_arn)

    def untag_resource(self, resource_arn: str, tag_keys: List[str]) -> None:
        self.tagger.untag_resource_using_names(resource_arn, tag_keys)

    def create_scheduling_policy(
        self, name: str, fairshare_policy: Dict[str, Any], tags: Dict[str, str]
    ) -> SchedulingPolicy:
        policy = SchedulingPolicy(
            self.account_id, self.region_name, name, fairshare_policy, self, tags
        )
        self._scheduling_policies[policy.arn] = policy
        return self._scheduling_policies[policy.arn]

    def describe_scheduling_policies(self, arns: List[str]) -> List[SchedulingPolicy]:
        return [pol for arn, pol in self._scheduling_policies.items() if arn in arns]

    def list_scheduling_policies(self) -> List[str]:
        """
        Pagination is not yet implemented
        """
        return list(self._scheduling_policies.keys())

    def delete_scheduling_policy(self, arn: str) -> None:
        self._scheduling_policies.pop(arn, None)

    def update_scheduling_policy(
        self, arn: str, fairshare_policy: Dict[str, Any]
    ) -> None:
        self._scheduling_policies[arn].fairshare_policy = fairshare_policy


batch_backends = BackendDict(BatchBackend, "batch")
