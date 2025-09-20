from __future__ import annotations

from datetime import datetime, timedelta
from functools import cache, cached_property
from typing import Any, Dict, List, Optional, Tuple

from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import BaseModel, CloudFormationModel
from moto.core.utils import utcnow
from moto.ec2.models.instances import Instance as EC2Instance
from moto.emr.exceptions import (
    InvalidCluster,
    InvalidRequestException,
    ResourceNotFoundException,
    ValidationException,
)
from moto.utilities.utils import CamelToUnderscoresWalker, get_partition, load_resource

from .utils import (
    EmrSecurityGroupManager,
    random_cluster_id,
    random_instance_group_id,
    random_step_id,
)

EXAMPLE_AMI_ID = "ami-12c6146b"


class Application(BaseModel):
    def __init__(
        self, name: str, version: str, args: List[str], additional_info: Dict[str, str]
    ):
        self.additional_info = additional_info or {}
        self.args = args or []
        self.name = name or None
        self.version = version or None


class BootstrapAction(BaseModel):
    def __init__(self, name: str, script_bootstrap_action: dict[str, Any]):
        self.name = name
        self.script_bootstrap_action = script_bootstrap_action

    @property
    def script_path(self) -> Optional[str]:
        return self.script_bootstrap_action.get("path")

    @property
    def args(self) -> List[str]:
        return self.script_bootstrap_action.get("args", [])


class Instance(BaseModel):
    def __init__(
        self,
        ec2_instance: EC2Instance,
        instance_group: InstanceGroup,
        instance_fleet_id: Optional[str] = None,
        instance_id: Optional[str] = None,
    ):
        self.id = instance_id or random_instance_group_id()
        self.ec2_instance = ec2_instance
        self.instance_group = instance_group
        self.instance_fleet_id = instance_fleet_id or ""

    @property
    def ec2_instance_id(self) -> str:
        return self.ec2_instance.id

    @property
    def public_dns_name(self) -> str | None:
        return self.ec2_instance.public_dns

    @property
    def public_ip_address(self) -> str | None:
        return self.ec2_instance.public_ip

    @property
    def private_ip_address(self) -> str | None:
        return self.ec2_instance.private_ip

    @property
    def private_dns_name(self) -> str:
        return self.ec2_instance.private_dns

    @property
    def instance_group_id(self) -> str:
        return self.instance_group.id

    @property
    def market(self) -> str:
        return self.instance_group.market

    @property
    def instance_type(self) -> str:
        return self.ec2_instance.instance_type

    @property
    def ebs_volumes(self) -> list[dict[str, str]]:
        ebs_volumes = []
        for volume in self.ec2_instance.block_device_mapping:
            member = {
                "Device": volume,
                "VolumeId": self.ec2_instance.block_device_mapping[volume].volume_id,
            }
            ebs_volumes.append(member)
        return ebs_volumes

    @property
    def status(self) -> dict[str, Any]:
        status_dict = {
            "State": self.instance_group.state,
            "StateChangeReason": {
                "Message": getattr(self, "state_change_reason", None),
            },
            "Timeline": {
                "CreationDateTime": self.instance_group.creation_date_time,
                "EndDateTime": self.instance_group.end_date_time,
                "ReadyDateTime": self.instance_group.ready_date_time,
            },
        }
        return status_dict


class InstanceGroup(CloudFormationModel):
    def __init__(
        self,
        cluster_id: str,
        instance_count: int,
        instance_role: str,
        instance_type: str,
        market: str = "ON_DEMAND",
        name: Optional[str] = None,
        instance_group_id: Optional[str] = None,
        bid_price: Optional[str] = None,
        ebs_configuration: Optional[Dict[str, Any]] = None,
        auto_scaling_policy: Optional[Dict[str, Any]] = None,
    ) -> None:
        self.id = instance_group_id or random_instance_group_id()
        self.cluster_id = cluster_id

        self.bid_price = bid_price
        self.market = market
        if name is None:
            if instance_role == "MASTER":
                name = "master"
            elif instance_role == "CORE":
                name = "slave"
            else:
                name = "Task instance group"
        self.name = name
        self.num_instances = instance_count
        self.role = instance_role
        self.instance_group_type = self.role
        self.instance_type = instance_type
        self.ebs_configuration = ebs_configuration
        self.auto_scaling_policy = auto_scaling_policy
        self.creation_date_time: datetime = utcnow()
        self.start_date_time: datetime = utcnow()
        self.ready_date_time: datetime = utcnow()
        self.end_date_time: Optional[datetime] = None
        self.state: str = "RUNNING"

    def set_instance_count(self, instance_count: int) -> None:
        self.num_instances = instance_count

    @property
    def status(self) -> dict[str, Any]:
        status_dict = {
            "State": self.state,
            "StateChangeReason": {
                "Message": getattr(self, "state_change_reason", None),
                "Code": "USER_REQUEST",
            },
            "Timeline": {
                "CreationDateTime": self.creation_date_time,
                "EndDateTime": self.end_date_time,
                "ReadyDateTime": self.ready_date_time,
            },
        }
        return status_dict

    @property
    def ebs_block_devices(self) -> list[dict[str, Any]]:
        ebs_block_devices = []
        if self.ebs_configuration is not None:
            for ebs_block_device in self.ebs_configuration.get(
                "ebs_block_device_configs", []
            ):
                for i in range(ebs_block_device.get("volumes_per_instance", 0)):
                    volume_spec = ebs_block_device["volume_specification"]
                    member = {
                        "VolumeSpecification": {
                            "VolumeType": volume_spec.get("volume_type"),
                            "Iops": volume_spec.get("iops"),
                            "SizeInGB": volume_spec.get("size_in_gb"),
                        },
                        "Device": f"/dev/sd{i}",
                    }
                    ebs_block_devices.append(member)
        return ebs_block_devices

    @property
    def instance_request_count(self) -> int:
        return self.num_instances

    @property
    def instance_role(self) -> str:
        return self.role

    @property
    def instance_running_count(self) -> int:
        return self.num_instances

    @property
    def requested_instance_count(self) -> int:
        return self.num_instances

    @property
    def running_instance_count(self) -> int:
        return self.num_instances

    @property
    def auto_scaling_policy(self) -> Any:
        return self._auto_scaling_policy

    @auto_scaling_policy.setter
    def auto_scaling_policy(self, value: Any) -> None:
        if value is None:
            self._auto_scaling_policy = value
            return
        self._auto_scaling_policy = CamelToUnderscoresWalker.parse(value)
        self._auto_scaling_policy["status"] = {"state": "ATTACHED"}
        # Transform common ${emr.clusterId} placeholder in any dimensions it occurs in.
        if "rules" in self._auto_scaling_policy:
            for rule in self._auto_scaling_policy["rules"]:
                if (
                    "trigger" in rule
                    and "cloud_watch_alarm_definition" in rule["trigger"]
                    and "dimensions" in rule["trigger"]["cloud_watch_alarm_definition"]
                ):
                    for dimension in rule["trigger"]["cloud_watch_alarm_definition"][
                        "dimensions"
                    ]:
                        if (
                            "value" in dimension
                            and dimension["value"] == "${emr.clusterId}"
                        ):
                            dimension["value"] = self.cluster_id

    @property
    def physical_resource_id(self) -> str:
        return self.id

    @staticmethod
    def cloudformation_type() -> str:
        return "AWS::EMR::InstanceGroupConfig"

    @classmethod
    def create_from_cloudformation_json(
        cls,
        resource_name: str,
        cloudformation_json: Any,
        account_id: str,
        region_name: str,
        **kwargs: Any,
    ) -> InstanceGroup:
        properties = cloudformation_json["Properties"]
        job_flow_id = properties["JobFlowId"]
        ebs_config = properties.get("EbsConfiguration")
        if ebs_config:
            ebs_config = CamelToUnderscoresWalker.parse_dict(ebs_config)
        props = {
            "instance_count": properties.get("InstanceCount"),
            "instance_role": properties.get("InstanceRole"),
            "instance_type": properties.get("InstanceType"),
            "market": properties.get("Market"),
            "bid_price": properties.get("BidPrice"),
            "name": properties.get("Name"),
            "auto_scaling_policy": properties.get("AutoScalingPolicy"),
            "ebs_configuration": ebs_config,
        }

        emr_backend: ElasticMapReduceBackend = emr_backends[account_id][region_name]
        return emr_backend.add_instance_groups(
            cluster_id=job_flow_id, instance_groups=[props]
        )[0]


class Step(BaseModel):
    def __init__(
        self,
        state: str,
        name: str = "",
        hadoop_jar_step: Optional[Dict[str, Any]] = None,
        action_on_failure: str = "TERMINATE_CLUSTER",
    ):
        self.id = random_step_id()
        self.action_on_failure = action_on_failure
        self.name = name
        self.hadoop_jar_step = hadoop_jar_step or {}
        self.jar = self.hadoop_jar_step.get("Jar", "")
        self.args = self.hadoop_jar_step.get("Args", [])
        self.properties = self.hadoop_jar_step.get("Properties", {})
        self.creation_date_time = utcnow()
        self.end_date_time = None
        self.ready_date_time = None
        self.start_date_time: Optional[datetime] = None
        self.state = state

    @property
    def config(self) -> dict[str, Any]:
        config = {
            "ActionOnFailure": self.action_on_failure,
            "HadoopJarStep": self.hadoop_jar_step,
            "Name": self.name,
            "Jar": self.jar,
            "Args": self.args,
            "Properties": {p["Key"]: p["Value"] for p in self.properties},
        }
        return config

    @property
    def execution_status_detail(self) -> dict[str, Any]:
        status_dict = {
            "CreationDateTime": self.creation_date_time,
            "StartDateTime": self.start_date_time,
            "ReadyDateTime": self.ready_date_time,
            "EndDateTime": self.end_date_time,
            "LastStateChangeReason": getattr(self, "last_state_change_reason", None),
            "State": self.state,
        }
        return status_dict

    @property
    def status(self) -> dict[str, Any]:
        status_dict = {
            "State": self.state,
            "StateChangeReason": {
                "Message": getattr(self, "last_state_change_reason", None)
            },
            "Timeline": {
                "CreationDateTime": self.creation_date_time,
                "EndDateTime": self.end_date_time,
                "StartDateTime": self.start_date_time,
            },
        }
        return status_dict

    def start(self) -> None:
        self.start_date_time = utcnow()


class Ec2InstanceAttributes:
    def __init__(self, cluster: Cluster) -> None:
        self.cluster = cluster
        self.additional_master_security_groups = (
            cluster.additional_master_security_groups
        )
        self.additional_slave_security_groups = cluster.additional_slave_security_groups
        self.ec2_availability_zone = cluster.availability_zone
        self.ec2_key_name = cluster.ec2_key_name
        self.ec2_subnet_id = cluster.ec2_subnet_id
        self.iam_instance_profile = cluster.role
        self.emr_managed_master_security_group = cluster.master_security_group
        self.emr_managed_slave_security_group = cluster.slave_security_group
        self.service_access_security_group = cluster.service_access_security_group


class Cluster(CloudFormationModel):
    def __init__(
        self,
        emr_backend: ElasticMapReduceBackend,
        name: str,
        log_uri: str,
        job_flow_role: str,
        service_role: str,
        steps: List[Dict[str, Any]],
        instance_attrs: Dict[str, Any],
        bootstrap_actions: Optional[List[Dict[str, Any]]] = None,
        configurations: Optional[List[Dict[str, Any]]] = None,
        cluster_id: Optional[str] = None,
        visible_to_all_users: bool = False,
        release_label: Optional[str] = None,
        requested_ami_version: Optional[str] = None,
        running_ami_version: Optional[str] = None,
        custom_ami_id: Optional[str] = None,
        step_concurrency_level: int = 1,
        security_configuration: Optional[str] = None,
        kerberos_attributes: Optional[Dict[str, str]] = None,
        auto_scaling_role: Optional[str] = None,
        ebs_root_volume_size: Optional[int] = None,
        ebs_root_volume_iops: Optional[int] = None,
        ebs_root_volume_throughput: Optional[int] = None,
    ):
        self.id = cluster_id or random_cluster_id()
        emr_backend.clusters[self.id] = self
        self.emr_backend = emr_backend

        self.applications: List[Application] = []

        self.bootstrap_actions: List[BootstrapAction] = []
        for bootstrap_action in bootstrap_actions or []:
            self.add_bootstrap_action(bootstrap_action)

        self.configurations = configurations or []

        self._tags: Dict[str, str] = {}

        self.log_uri = log_uri
        self.name = name
        self.normalized_instance_hours = 0

        self.steps: List[Step] = []
        self.add_steps(steps)

        self.set_visibility(visible_to_all_users)

        self.instance_group_ids: List[str] = []
        self.ec2_instances: List[Instance] = []
        self.master_instance_group_id: Optional[str] = None
        self.core_instance_group_id: Optional[str] = None
        if (
            "master_instance_type" in instance_attrs
            and instance_attrs["master_instance_type"]
        ):
            self.emr_backend.add_instance_groups(
                self.id,
                [
                    {
                        "instance_count": 1,
                        "instance_role": "MASTER",
                        "instance_type": instance_attrs["master_instance_type"],
                        "market": "ON_DEMAND",
                        "name": "master",
                    }
                ],
            )
        if (
            "slave_instance_type" in instance_attrs
            and instance_attrs["slave_instance_type"]
        ):
            self.emr_backend.add_instance_groups(
                self.id,
                [
                    {
                        "instance_count": instance_attrs["instance_count"] - 1,
                        "instance_role": "CORE",
                        "instance_type": instance_attrs["slave_instance_type"],
                        "market": "ON_DEMAND",
                        "name": "slave",
                    }
                ],
            )
        self.additional_master_security_groups = instance_attrs.get(
            "additional_master_security_groups"
        )
        self.additional_slave_security_groups = instance_attrs.get(
            "additional_slave_security_groups"
        )
        self.availability_zone = instance_attrs.get("availability_zone")
        self.ec2_key_name = instance_attrs.get("ec2_key_name")
        self.ec2_subnet_id = instance_attrs.get("ec2_subnet_id")
        self.hadoop_version = instance_attrs.get("hadoop_version")
        self.keep_job_flow_alive_when_no_steps = instance_attrs.get(
            "keep_job_flow_alive_when_no_steps"
        )
        self.master_security_group = instance_attrs.get(
            "emr_managed_master_security_group"
        )
        self.service_access_security_group = instance_attrs.get(
            "service_access_security_group"
        )
        self.slave_security_group = instance_attrs.get(
            "emr_managed_slave_security_group"
        )
        self.termination_protected = instance_attrs.get("termination_protected", False)

        self.release_label = release_label
        self.requested_ami_version = requested_ami_version
        self.running_ami_version = running_ami_version
        self.custom_ami_id = custom_ami_id

        self.ebs_root_volume_size = ebs_root_volume_size
        self.ebs_root_volume_iops = ebs_root_volume_iops
        self.ebs_root_volume_throughput = ebs_root_volume_throughput

        self.role = job_flow_role or "EMRJobflowDefault"
        self.service_role = service_role
        self.step_concurrency_level = step_concurrency_level

        self.creation_datetime = utcnow()
        self.start_datetime: Optional[datetime] = None
        self.ready_datetime: Optional[datetime] = None
        self.end_datetime: Optional[datetime] = None
        self.state: Optional[str] = None

        self.start_cluster()
        self.run_bootstrap_actions()
        if self.steps:
            self.steps[0].start()
        self.security_configuration = (
            security_configuration  # ToDo: Raise if doesn't already exist.
        )
        self.kerberos_attributes = kerberos_attributes
        self.auto_scaling_role = auto_scaling_role
        self.supported_products: list[str] = []

    @property
    def arn(self) -> str:
        return f"arn:{get_partition(self.emr_backend.region_name)}:elasticmapreduce:{self.emr_backend.region_name}:{self.emr_backend.account_id}:cluster/{self.id}"

    @property
    def ami_version(self) -> str | None:
        return self.running_ami_version

    @property
    def auto_terminate(self) -> bool:
        return not self.keep_job_flow_alive_when_no_steps

    @property
    def bootstrap_action_detail_list(self) -> list[dict[str, Any]]:
        details = []
        for bootstrap_action in self.bootstrap_actions:
            member = {"BootstrapActionConfig": bootstrap_action}
            details.append(member)
        return details

    @property
    def ec2_instance_attributes(self) -> Ec2InstanceAttributes:
        return Ec2InstanceAttributes(self)

    @property
    def job_flow_instances_detail(self) -> dict[str, Any]:
        instances_detail = {
            "Ec2KeyName": self.ec2_key_name,
            "Ec2SubnetId": self.ec2_subnet_id,
            "HadoopVersion": self.hadoop_version,
            "InstanceCount": self.instance_count,
            "InstanceGroups": self.instance_groups,
            "KeepJobFlowAliveWhenNoSteps": self.keep_job_flow_alive_when_no_steps,
            "MasterInstanceId": getattr(self, "master_instance_id", None),
            "MasterInstanceType": self.master_instance_type,
            "MasterPublicDnsName": self.master_public_dns_name,
            "NormalizedInstanceHours": self.normalized_instance_hours,
            "Placement": {"AvailabilityZone": self.availability_zone},
            "SlaveInstanceType": self.slave_instance_type,
            "TerminationProtected": self.termination_protected,
        }
        return instances_detail

    @property
    def execution_status_detail(self) -> dict[str, Any]:
        status_dict = {
            "State": self.state,
            "CreationDateTime": self.creation_datetime,
            "StartDateTime": self.start_datetime,
            "ReadyDateTime": self.ready_datetime,
            "EndDateTime": self.end_datetime,
            "LastStateChangeReason": getattr(self, "last_state_change_reason", None),
        }
        return status_dict

    @property
    def job_flow_id(self) -> str:
        return self.id

    @property
    def job_flow_role(self) -> str:
        return self.role

    @property
    def master_public_dns_name(self) -> str:
        if self.master_instance_group_id:
            master_instance_group = self.emr_backend.instance_groups[
                self.master_instance_group_id  # type: ignore
            ]

            for ec2_instance in self.ec2_instances:
                if ec2_instance.instance_group.id == master_instance_group.id:
                    dashed_ip, domain = ec2_instance.public_dns_name.split(".", 1)  # type: ignore
                    return ".".join([dashed_ip, self.emr_backend.region_name, domain])

        # TODO: Moto does not create instances when InstanceGroups is not provided to RunJobFlow operation.
        #  In fact it should also do so when InstanceCount/InstanceType are provided.
        #  For now, just return a hardcoded value for such cases.
        return "ec2-184-0-0-1.us-west-1.compute.amazonaws.com"

    @property
    def status(self) -> dict[str, Any]:
        status_dict = {
            "State": self.state,
            "StateChangeReason": {
                "Message": getattr(self, "last_state_change_reason", None)
            },
            "Timeline": {
                "CreationDateTime": self.creation_datetime,
                "EndDateTime": self.end_datetime,
                "ReadyDateTime": self.ready_datetime,
            },
        }
        return status_dict

    @property
    def tags(self) -> list[dict[str, str]]:
        tags = [{"Key": k, "Value": v} for k, v in self._tags.items()]
        return tags

    @property
    def instance_groups(self) -> List[InstanceGroup]:
        return self.emr_backend.get_instance_groups(self.instance_group_ids)

    @property
    def master_instance_type(self) -> str:
        return self.emr_backend.instance_groups[
            self.master_instance_group_id  # type: ignore
        ].instance_type

    @property
    def slave_instance_type(self) -> str:
        return self.emr_backend.instance_groups[
            self.core_instance_group_id  # type: ignore
        ].instance_type

    @property
    def instance_count(self) -> int:
        return sum(group.num_instances for group in self.instance_groups)

    def start_cluster(self) -> None:
        self.state = "STARTING"
        self.start_datetime = utcnow()

    def run_bootstrap_actions(self) -> None:
        self.state = "BOOTSTRAPPING"
        self.ready_datetime = utcnow()
        self.state = "WAITING"
        if not self.steps:
            if not self.keep_job_flow_alive_when_no_steps:
                self.terminate()

    def terminate(self) -> None:
        self.state = "TERMINATING"
        self.end_datetime = utcnow()
        self.state = "TERMINATED"

    def add_applications(self, applications: List[Dict[str, Any]]) -> None:
        self.applications.extend(
            [
                Application(
                    name=app.get("name", ""),
                    version=app.get("version", ""),
                    args=app.get("args", []),
                    additional_info=app.get("additional_info", {}),
                )
                for app in applications
            ]
        )

    def add_bootstrap_action(self, bootstrap_action: Dict[str, Any]) -> None:
        self.bootstrap_actions.append(BootstrapAction(**bootstrap_action))

    def add_instance_group(self, instance_group: InstanceGroup) -> None:
        if instance_group.role == "MASTER":
            if self.master_instance_group_id:
                raise Exception("Cannot add another master instance group")
            self.master_instance_group_id = instance_group.id
            num_master_nodes = instance_group.num_instances
            if num_master_nodes > 1:
                # Cluster is HA
                if num_master_nodes != 3:
                    raise ValidationException(
                        "Master instance group must have exactly 3 instances for HA clusters."
                    )
                self.keep_job_flow_alive_when_no_steps = True
                self.termination_protected = True
        if instance_group.role == "CORE":
            if self.core_instance_group_id:
                raise Exception("Cannot add another core instance group")
            self.core_instance_group_id = instance_group.id
        self.instance_group_ids.append(instance_group.id)

    def add_instance(self, instance: Instance) -> None:
        self.ec2_instances.append(instance)

    def add_steps(self, steps: List[Dict[str, Any]]) -> List[Step]:
        added_steps = []
        for step in steps:
            if self.steps:
                # If we already have other steps, this one is pending
                fake = Step(state="PENDING", **step)
            else:
                fake = Step(state="RUNNING", **step)
            self.steps.append(fake)
            added_steps.append(fake)
        self.state = "RUNNING"
        return added_steps

    def add_tags(self, tags: Dict[str, str]) -> None:
        self._tags.update(tags)

    def remove_tags(self, tag_keys: List[str]) -> None:
        for key in tag_keys:
            self._tags.pop(key, None)

    def set_termination_protection(self, value: bool) -> None:
        self.termination_protected = value

    def set_visibility(self, visibility: bool) -> None:
        self.visible_to_all_users = visibility

    @property
    def physical_resource_id(self) -> str:
        return self.id

    @staticmethod
    def cloudformation_type() -> str:
        return "AWS::EMR::Cluster"

    @classmethod
    def has_cfn_attr(cls, attr: str) -> bool:
        return attr in ["Id"]

    def get_cfn_attribute(self, attribute_name: str) -> str:
        from moto.cloudformation.exceptions import UnformattedGetAttTemplateException

        if attribute_name == "Id":
            return self.id
        raise UnformattedGetAttTemplateException()

    @classmethod
    def create_from_cloudformation_json(
        cls,
        resource_name: str,
        cloudformation_json: Any,
        account_id: str,
        region_name: str,
        **kwargs: Any,
    ) -> Cluster:
        properties = cloudformation_json["Properties"]

        instance_attrs = properties.get("Instances", {})
        instance_attrs["ec2_subnet_id"] = instance_attrs.get("Ec2SubnetId")
        instance_attrs["emr_managed_master_security_group"] = instance_attrs.get(
            "EmrManagedMasterSecurityGroup"
        )
        instance_attrs["emr_managed_slave_security_group"] = instance_attrs.get(
            "EmrManagedSlaveSecurityGroup"
        )
        instance_attrs["service_access_security_group"] = instance_attrs.get(
            "ServiceAccessSecurityGroup"
        )
        instance_attrs["additional_master_security_groups"] = instance_attrs.get(
            "AdditionalMasterSecurityGroups", []
        )
        instance_attrs["additional_slave_security_groups"] = instance_attrs.get(
            "AdditionalSlaveSecurityGroups", []
        )

        emr_backend: ElasticMapReduceBackend = emr_backends[account_id][region_name]
        cluster = emr_backend.run_job_flow(
            name=properties["Name"],
            log_uri=properties.get("LogUri"),
            job_flow_role=properties["JobFlowRole"],
            service_role=properties["ServiceRole"],
            steps=[],
            instance_attrs=instance_attrs,
            kerberos_attributes=properties.get("KerberosAttributes", {}),
            release_label=properties.get("ReleaseLabel"),
            custom_ami_id=properties.get("CustomAmiId"),
        )
        tags = {item["Key"]: item["Value"] for item in properties.get("Tags", [])}
        cluster.add_tags(tags)
        return cluster

    @classmethod
    def delete_from_cloudformation_json(
        cls,
        resource_name: str,
        cloudformation_json: Dict[str, Any],
        account_id: str,
        region_name: str,
    ) -> None:
        emr_backend: ElasticMapReduceBackend = emr_backends[account_id][region_name]

        emr_backend.terminate_job_flows([resource_name])


class SecurityConfiguration(CloudFormationModel):
    def __init__(self, name: str, security_configuration: str):
        self.name = name
        self.security_configuration = security_configuration
        self.creation_date_time = utcnow()

    @property
    def physical_resource_id(self) -> str:
        return self.name

    @staticmethod
    def cloudformation_type() -> str:
        return "AWS::EMR::SecurityConfiguration"

    @classmethod
    def create_from_cloudformation_json(
        cls,
        resource_name: str,
        cloudformation_json: Any,
        account_id: str,
        region_name: str,
        **kwargs: Any,
    ) -> SecurityConfiguration:
        emr_backend: ElasticMapReduceBackend = emr_backends[account_id][region_name]

        properties = cloudformation_json["Properties"]
        return emr_backend.create_security_configuration(
            name=properties.get("Name") or resource_name,
            security_configuration=properties.get("SecurityConfiguration", {}),
        )

    @classmethod
    def delete_from_cloudformation_json(
        cls,
        resource_name: str,
        cloudformation_json: Dict[str, Any],
        account_id: str,
        region_name: str,
    ) -> None:
        emr_backend: ElasticMapReduceBackend = emr_backends[account_id][region_name]

        properties = cloudformation_json["Properties"]
        name = properties.get("Name") or resource_name
        emr_backend.delete_security_configuration(name)


class ElasticMapReduceBackend(BaseBackend):
    def __init__(self, region_name: str, account_id: str):
        super().__init__(region_name, account_id)
        self.clusters: Dict[str, Cluster] = {}
        self.instance_groups: Dict[str, InstanceGroup] = {}
        self.security_configurations: Dict[str, SecurityConfiguration] = {}
        self.block_public_access_configuration: Dict[str, Any] = {}

    @cached_property
    def _release_labels(self) -> list[str]:
        """Returns all available release labels in this region"""
        return load_resource(
            __name__, f"resources/release-labels-{self.region_name}.json"
        )

    @cache
    def _get_instance_type_names(self, release_label: str) -> list[str]:
        """Returns all instance type names that can be used with this release label"""
        return load_resource(__name__, f"resources/instance-types-{release_label}.json")

    @cached_property
    def _instance_types(self) -> dict[str, Any]:
        """Returns a dictionary of {instance-type-name: instance-type-details}"""
        return load_resource(__name__, "resources/instance_types.json")

    @property
    def ec2_backend(self) -> Any:
        """
        :return: EC2 Backend
        :rtype: moto.ec2.models.EC2Backend
        """
        from moto.ec2 import ec2_backends

        return ec2_backends[self.account_id][self.region_name]

    def add_applications(
        self, cluster_id: str, applications: List[Dict[str, Any]]
    ) -> None:
        cluster = self.describe_cluster(cluster_id)
        cluster.add_applications(applications)

    def add_instance_groups(
        self, cluster_id: str, instance_groups: List[Dict[str, Any]]
    ) -> List[InstanceGroup]:
        cluster = self.clusters[cluster_id]
        result_groups = []
        for instance_group in instance_groups:
            group = InstanceGroup(cluster_id=cluster_id, **instance_group)
            self.instance_groups[group.id] = group
            cluster.add_instance_group(group)
            result_groups.append(group)
        return result_groups

    def run_instances(
        self,
        cluster_id: str,
        instances: Dict[str, Any],
        instance_group: InstanceGroup,
    ) -> None:
        cluster = self.clusters[cluster_id]
        instances["is_instance_type_default"] = not instances.get("instance_type")
        response = self.ec2_backend.run_instances(
            EXAMPLE_AMI_ID, instances["instance_count"], "", [], **instances
        )
        for instance in response.instances:
            instance = Instance(ec2_instance=instance, instance_group=instance_group)
            cluster.add_instance(instance)

    def add_job_flow_steps(
        self, job_flow_id: str, steps: List[Dict[str, Any]]
    ) -> List[Step]:
        cluster = self.clusters[job_flow_id]
        return cluster.add_steps(steps)

    def add_tags(self, cluster_id: str, tags: Dict[str, str]) -> None:
        cluster = self.describe_cluster(cluster_id)
        cluster.add_tags(tags)

    def describe_job_flows(
        self,
        job_flow_ids: Optional[List[str]] = None,
        job_flow_states: Optional[List[str]] = None,
        created_after: Optional[datetime] = None,
        created_before: Optional[datetime] = None,
    ) -> List[Cluster]:
        clusters = list(self.clusters.values())

        within_two_month = utcnow() - timedelta(days=60)
        clusters = [c for c in clusters if c.creation_datetime >= within_two_month]

        if job_flow_ids:
            clusters = [c for c in clusters if c.id in job_flow_ids]
        if job_flow_states:
            clusters = [c for c in clusters if c.state in job_flow_states]
        if created_after:
            clusters = [c for c in clusters if c.creation_datetime > created_after]
        if created_before:
            clusters = [c for c in clusters if c.creation_datetime < created_before]

        # Amazon EMR can return a maximum of 512 job flow descriptions
        return sorted(clusters, key=lambda x: x.id)[:512]

    def describe_step(self, cluster_id: str, step_id: str) -> Optional[Step]:
        cluster = self.clusters[cluster_id]
        for step in cluster.steps:
            if step.id == step_id:
                return step
        return None

    def describe_cluster(self, cluster_id: str) -> Cluster:
        if cluster_id in self.clusters:
            return self.clusters[cluster_id]
        raise InvalidCluster(cluster_id)

    def get_instance_groups(self, instance_group_ids: List[str]) -> List[InstanceGroup]:
        return [
            group
            for group_id, group in self.instance_groups.items()
            if group_id in instance_group_ids
        ]

    def list_bootstrap_actions(
        self, cluster_id: str, marker: Optional[str] = None
    ) -> Tuple[List[BootstrapAction], Optional[str]]:
        max_items = 50
        actions = self.clusters[cluster_id].bootstrap_actions
        start_idx = 0 if marker is None else int(marker)
        marker = (
            None
            if len(actions) <= start_idx + max_items
            else str(start_idx + max_items)
        )
        return actions[start_idx : start_idx + max_items], marker

    def list_clusters(
        self,
        cluster_states: Optional[List[str]] = None,
        created_after: Optional[datetime] = None,
        created_before: Optional[datetime] = None,
        marker: Optional[str] = None,
    ) -> Tuple[List[Cluster], Optional[str]]:
        max_items = 50
        clusters = list(self.clusters.values())
        if cluster_states:
            clusters = [c for c in clusters if c.state in cluster_states]
        if created_after:
            clusters = [c for c in clusters if c.creation_datetime > created_after]
        if created_before:
            clusters = [c for c in clusters if c.creation_datetime < created_before]
        clusters = sorted(clusters, key=lambda x: x.id)
        start_idx = 0 if marker is None else int(marker)
        marker = (
            None
            if len(clusters) <= start_idx + max_items
            else str(start_idx + max_items)
        )
        return clusters[start_idx : start_idx + max_items], marker

    def list_instance_groups(
        self, cluster_id: str, marker: Optional[str] = None
    ) -> Tuple[List[InstanceGroup], Optional[str]]:
        max_items = 50
        groups = sorted(self.clusters[cluster_id].instance_groups, key=lambda x: x.id)
        start_idx = 0 if marker is None else int(marker)
        marker = (
            None if len(groups) <= start_idx + max_items else str(start_idx + max_items)
        )
        return groups[start_idx : start_idx + max_items], marker

    def list_instances(
        self,
        cluster_id: str,
        marker: Optional[str] = None,
        instance_group_id: Optional[str] = None,
        instance_group_types: Optional[List[str]] = None,
    ) -> Tuple[List[Instance], Optional[str]]:
        max_items = 50
        groups = sorted(self.clusters[cluster_id].ec2_instances, key=lambda x: x.id)
        start_idx = 0 if marker is None else int(marker)
        marker = (
            None if len(groups) <= start_idx + max_items else str(start_idx + max_items)
        )
        if instance_group_id:
            groups = [g for g in groups if g.instance_group.id == instance_group_id]
        if instance_group_types:
            groups = [
                g for g in groups if g.instance_group.role in instance_group_types
            ]
        return groups[start_idx : start_idx + max_items], marker

    def list_steps(
        self,
        cluster_id: str,
        marker: Optional[str] = None,
        step_ids: Optional[List[str]] = None,
        step_states: Optional[List[str]] = None,
    ) -> Tuple[List[Step], Optional[str]]:
        max_items = 50
        steps = sorted(
            self.clusters[cluster_id].steps,
            key=lambda o: o.creation_date_time,
            reverse=True,
        )
        if step_ids:
            steps = [s for s in steps if s.id in step_ids]
        if step_states:
            steps = [s for s in steps if s.state in step_states]
        start_idx = 0 if marker is None else int(marker)
        marker = (
            None if len(steps) <= start_idx + max_items else str(start_idx + max_items)
        )
        return steps[start_idx : start_idx + max_items], marker

    def modify_cluster(self, cluster_id: str, step_concurrency_level: int) -> Cluster:
        cluster = self.clusters[cluster_id]
        cluster.step_concurrency_level = step_concurrency_level
        return cluster

    def modify_instance_groups(self, instance_groups: List[Dict[str, Any]]) -> None:
        for instance_group in instance_groups:
            group = self.instance_groups[instance_group["instance_group_id"]]
            group.set_instance_count(int(instance_group["instance_count"]))

    def remove_tags(self, cluster_id: str, tag_keys: List[str]) -> None:
        cluster = self.describe_cluster(cluster_id)
        cluster.remove_tags(tag_keys)

    def _manage_security_groups(
        self,
        ec2_subnet_id: str,
        emr_managed_master_security_group: str,
        emr_managed_slave_security_group: str,
        service_access_security_group: str,
        **_: Any,
    ) -> Tuple[str, str, str]:
        default_return_value = (
            emr_managed_master_security_group,
            emr_managed_slave_security_group,
            service_access_security_group,
        )
        if not ec2_subnet_id:
            # TODO: Set up Security Groups in Default VPC.
            return default_return_value

        from moto.ec2.exceptions import InvalidSubnetIdError

        try:
            subnet = self.ec2_backend.get_subnet(ec2_subnet_id)
        except InvalidSubnetIdError:
            return default_return_value

        manager = EmrSecurityGroupManager(self.ec2_backend, subnet.vpc_id)
        master, slave, service = manager.manage_security_groups(
            emr_managed_master_security_group,
            emr_managed_slave_security_group,
            service_access_security_group,
        )
        return master.id, slave.id, service.id

    def run_job_flow(self, **kwargs: Any) -> Cluster:
        (
            kwargs["instance_attrs"]["emr_managed_master_security_group"],
            kwargs["instance_attrs"]["emr_managed_slave_security_group"],
            kwargs["instance_attrs"]["service_access_security_group"],
        ) = self._manage_security_groups(**kwargs["instance_attrs"])
        return Cluster(self, **kwargs)

    def set_visible_to_all_users(
        self, job_flow_ids: List[str], visible_to_all_users: bool
    ) -> None:
        for job_flow_id in job_flow_ids:
            cluster = self.clusters[job_flow_id]
            cluster.set_visibility(visible_to_all_users)

    def set_termination_protection(self, job_flow_ids: List[str], value: bool) -> None:
        for job_flow_id in job_flow_ids:
            cluster = self.clusters[job_flow_id]
            cluster.set_termination_protection(value)

    def terminate_job_flows(self, job_flow_ids: List[str]) -> List[Cluster]:
        clusters_terminated = []
        clusters_protected = []
        for job_flow_id in job_flow_ids:
            cluster = self.clusters[job_flow_id]
            if cluster.termination_protected:
                clusters_protected.append(cluster)
                continue
            cluster.terminate()
            clusters_terminated.append(cluster)
        if clusters_protected:
            raise ValidationException(
                "Could not shut down one or more job flows since they are termination protected."
            )
        return clusters_terminated

    def put_auto_scaling_policy(
        self, instance_group_id: str, auto_scaling_policy: Optional[Dict[str, Any]]
    ) -> Optional[InstanceGroup]:
        instance_groups = self.get_instance_groups(
            instance_group_ids=[instance_group_id]
        )
        if len(instance_groups) == 0:
            return None
        instance_group = instance_groups[0]
        instance_group.auto_scaling_policy = auto_scaling_policy
        return instance_group

    def remove_auto_scaling_policy(self, instance_group_id: str) -> None:
        self.put_auto_scaling_policy(instance_group_id, auto_scaling_policy=None)

    def create_security_configuration(
        self, name: str, security_configuration: str
    ) -> SecurityConfiguration:
        if name in self.security_configurations:
            raise InvalidRequestException(
                message=f"SecurityConfiguration with name '{name}' already exists."
            )
        config = SecurityConfiguration(
            name=name, security_configuration=security_configuration
        )
        self.security_configurations[name] = config
        return config

    def get_security_configuration(self, name: str) -> SecurityConfiguration:
        if name not in self.security_configurations:
            raise ResourceNotFoundException(
                f"Security configuration with name '{name}' does not exist."
            )
        return self.security_configurations[name]

    def delete_security_configuration(self, name: str) -> None:
        if name not in self.security_configurations:
            raise ResourceNotFoundException(
                f"Security configuration with name '{name}' does not exist."
            )
        del self.security_configurations[name]

    def get_block_public_access_configuration(
        self,
    ) -> Dict[str, Any]:  # type ignore[misc]
        return self.block_public_access_configuration

    def put_block_public_access_configuration(
        self,
        block_public_security_group_rules: bool,
        rule_ranges: Optional[List[Dict[str, int]]],
    ) -> None:
        from moto.sts import sts_backends

        sts_backend = sts_backends[self.account_id]["global"]
        _, user_arn, _ = sts_backend.get_caller_identity(
            self.account_id, region=self.region_name
        )
        self.block_public_access_configuration = {
            "block_public_access_configuration": {
                "block_public_security_group_rules": block_public_security_group_rules,
                "permitted_public_security_group_rule_ranges": [
                    {
                        "min_range": rule_range.get("MinRange"),
                        "max_range": rule_range.get("MaxRange"),
                    }
                    for rule_range in rule_ranges or []
                ],
            },
            "block_public_access_configuration_metadata": {
                "creation_date_time": utcnow(),
                "created_by_arn": user_arn,
            },
        }
        return

    def list_release_labels(self) -> list[str]:
        """
        Pagination and Filtering is not yet implemented
        """
        return self._release_labels

    def list_supported_instance_types(self, release_label: str) -> list[dict[str, Any]]:
        """
        Pagination is not yet implemented
        """
        instance_type_names = self._get_instance_type_names(release_label)
        return [
            details
            for name, details in self._instance_types.items()
            if name in instance_type_names
        ]


emr_backends = BackendDict(ElasticMapReduceBackend, "emr")
