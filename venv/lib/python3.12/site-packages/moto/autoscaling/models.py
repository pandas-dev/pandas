from __future__ import annotations

import itertools
import math
from collections import OrderedDict
from datetime import datetime
from typing import Any, Optional, Union

from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import BaseModel, CloudFormationModel
from moto.core.types import Base64EncodedString
from moto.core.utils import utcnow
from moto.ec2 import ec2_backends
from moto.ec2.exceptions import InvalidInstanceIdError
from moto.ec2.models import EC2Backend
from moto.ec2.models.instances import Instance
from moto.ec2.models.launch_templates import LaunchTemplate
from moto.elb.exceptions import LoadBalancerNotFoundError
from moto.elb.models import ELBBackend, elb_backends
from moto.elbv2.models import ELBv2Backend, elbv2_backends
from moto.moto_api._internal import mock_random as random
from moto.packages.boto.ec2.blockdevicemapping import (
    BlockDeviceMapping,
    BlockDeviceType,
)
from moto.utilities.utils import get_partition

from .exceptions import (
    AutoscalingClientError,
    InvalidInstanceError,
    ResourceContentionError,
    ValidationError,
)

# http://docs.aws.amazon.com/AutoScaling/latest/DeveloperGuide/AS_Concepts.html#Cooldown
DEFAULT_COOLDOWN = 300

ASG_NAME_TAG = "aws:autoscaling:groupName"


def make_int(value: Union[None, str, int]) -> Optional[int]:
    return int(value) if value is not None else value


class Activity:
    def __init__(
        self,
        description: str,
        cause: str,
        auto_scaling_group: FakeAutoScalingGroup,
        activity_id: Optional[str] = None,
        start_time: Optional[datetime] = None,
        end_time: Optional[datetime] = None,
        status_code: str = "InProgress",
    ):
        self.activity_id = activity_id or str(random.uuid4())
        self.auto_scaling_group = auto_scaling_group
        self.description = description
        self.cause = cause
        self.start_time = start_time or utcnow()
        self.end_time = end_time or utcnow()
        self.status_code = status_code
        self.progress = 0

    @property
    def auto_scaling_group_name(self) -> str:
        return self.auto_scaling_group.name


class TerminateInstanceActivity(Activity):
    def __init__(self, instance: Instance, original_capacity: int):
        auto_scaling_group = instance.autoscaling_group  # type: ignore[attr-defined]
        desired_capacity = auto_scaling_group.desired_capacity
        should_decrement = desired_capacity < original_capacity
        description = f"Terminating EC2 instance: {instance.id}"
        timestamp = utcnow()
        cause = f"At {timestamp}, instance {instance.id} was taken out of service in response to a user request"
        if should_decrement:
            cause += f", shrinking the capacity from {original_capacity} to {desired_capacity}."
        else:
            cause += "."
        super().__init__(description, cause, auto_scaling_group, start_time=timestamp)


class EnterStandbyActivity(Activity):
    def __init__(self, instance: Instance, original_capacity: Optional[int] = None):
        auto_scaling_group = instance.autoscaling_group  # type: ignore[attr-defined]
        desired_capacity = auto_scaling_group.desired_capacity
        should_decrement = desired_capacity < original_capacity
        description = f"Moving EC2 instance to StandBy: {instance.id}"
        timestamp = utcnow()
        cause = f"At {timestamp}, instance {instance.id} was moved to standby in response to a user request"
        if should_decrement:
            cause += f", shrinking the capacity from {original_capacity} to {desired_capacity}."
        else:
            cause += "."
        super().__init__(description, cause, auto_scaling_group, start_time=timestamp)
        self.progress = 50


class ExitStandbyActivity(Activity):
    def __init__(self, instance: Instance, original_capacity: Optional[int] = None):
        auto_scaling_group = instance.autoscaling_group  # type: ignore[attr-defined]
        desired_capacity = auto_scaling_group.desired_capacity
        description = f"Moving EC2 instance out of StandBy: {instance.id}"
        timestamp = utcnow()
        cause = f"At {timestamp}, instance {instance.id} was moved out of standby in response to a user request, increasing the capacity from {original_capacity} to {desired_capacity}."
        super().__init__(description, cause, auto_scaling_group, start_time=timestamp)
        self.progress = 30
        self.status_code = "PreInService"


class DetachInstanceActivity(Activity):
    def __init__(self, instance: Instance):
        auto_scaling_group = instance.autoscaling_group  # type: ignore[attr-defined]
        description = f"Detaching EC2 instance: {instance.id}"
        timestamp = utcnow()
        cause = f"At {timestamp}, instance {instance.id} was detached in response to a user request."
        super().__init__(description, cause, auto_scaling_group, start_time=timestamp)
        self.progress = 50


class InstanceState:
    def __init__(
        self,
        instance: Instance,
        lifecycle_state: str = "InService",
        health_status: str = "Healthy",
        protected_from_scale_in: Optional[bool] = False,
        autoscaling_group: Optional[FakeAutoScalingGroup] = None,
    ):
        self.instance = instance
        self.lifecycle_state = lifecycle_state
        self.health_status = health_status
        self.protected_from_scale_in = protected_from_scale_in
        if not hasattr(self.instance, "autoscaling_group"):
            self.instance.autoscaling_group = autoscaling_group  # type: ignore[attr-defined]
        self.auto_scaling_group = self.instance.autoscaling_group  # type: ignore[attr-defined]
        self.auto_scaling_group_name = self.auto_scaling_group.name
        self.availability_zone = self.instance.placement  # type: ignore[attr-defined]
        self.instance_id = self.instance.id
        self.instance_type = self.instance.instance_type

    @property
    def launch_template(self) -> Optional[dict[str, Any]]:
        if (
            self.auto_scaling_group is not None
            and self.auto_scaling_group.ec2_launch_template is None
        ):
            return None
        lt = {
            "LaunchTemplateId": self.auto_scaling_group.ec2_launch_template.id,
            "LaunchTemplateName": self.auto_scaling_group.ec2_launch_template.name,
            "Version": self.auto_scaling_group.ec2_launch_template.default_version_number,
        }
        return lt

    @property
    def launch_configuration_name(self) -> Optional[str]:
        return (
            self.auto_scaling_group.launch_configuration_name
            if self.auto_scaling_group is not None
            else None
        )


class LifecycleHook(BaseModel):
    def __init__(
        self,
        name: str,
        as_name: str,
        transition: Optional[str],
        timeout: Optional[int],
        result: Optional[str],
    ):
        self.name = name
        self.auto_scaling_group_name = as_name
        self.lifecycle_transition = transition
        self.heartbeat_timeout = timeout or 3600
        self.default_result = result or "ABANDON"
        # TODO: These were hardcoded in the original XML template, but should be implemented properly.
        self.role_arn = "arn:aws:iam::1234567890:role/my-auto-scaling-role"
        self.notification_target_arn = "arn:aws:sqs:us-east-1:123456789012:my-queue"
        self.global_timeout = 172800


class TargetTrackingConfiguration:
    def __init__(self, data: Optional[dict[str, Any]]) -> None:
        data = data or {}
        customized_metric_spec = data.get("CustomizedMetricSpecification", {})
        if customized_metric_spec:
            if "Dimensions" not in customized_metric_spec:
                customized_metric_spec["Dimensions"] = []
            for metric in customized_metric_spec.get("Metrics", []):
                if "ReturnData" not in metric:
                    metric["ReturnData"] = True
        self.__dict__.update(data)


class FakeScalingPolicy(BaseModel):
    def __init__(
        self,
        name: str,
        policy_type: str,
        metric_aggregation_type: str,
        adjustment_type: str,
        as_name: str,
        min_adjustment_magnitude: str,
        scaling_adjustment: Optional[int],
        cooldown: Optional[int],
        target_tracking_config: dict[str, Any],
        step_adjustments: str,
        estimated_instance_warmup: str,
        predictive_scaling_configuration: str,
        autoscaling_backend: AutoScalingBackend,
    ):
        self.name = name
        self.policy_name = name  # property alias
        self.policy_type = policy_type
        self.metric_aggregation_type = metric_aggregation_type
        self.adjustment_type = adjustment_type
        self.auto_scaling_group_name = as_name
        self.min_adjustment_magnitude = min_adjustment_magnitude
        self.scaling_adjustment = scaling_adjustment
        self.cooldown = None
        if self.policy_type == "SimpleScaling":
            self.cooldown = cooldown if cooldown is not None else DEFAULT_COOLDOWN
        self.target_tracking_configuration: Optional[TargetTrackingConfiguration] = None
        if self.policy_type == "TargetTrackingScaling":
            self.target_tracking_configuration = TargetTrackingConfiguration(
                target_tracking_config
            )
        self.step_adjustments = step_adjustments
        self.estimated_instance_warmup = estimated_instance_warmup
        self.predictive_scaling_configuration = predictive_scaling_configuration
        self.autoscaling_backend = autoscaling_backend

    @property
    def arn(self) -> str:
        return f"arn:{get_partition(self.autoscaling_backend.region_name)}:autoscaling:{self.autoscaling_backend.region_name}:{self.autoscaling_backend.account_id}:scalingPolicy:c322761b-3172-4d56-9a21-0ed9d6161d67:autoScalingGroupName/{self.auto_scaling_group_name}:policyName/{self.name}"

    policy_arn = arn  # property alias

    def execute(self) -> None:
        if self.adjustment_type == "ExactCapacity":
            self.autoscaling_backend.set_desired_capacity(
                self.auto_scaling_group_name, self.scaling_adjustment
            )
        elif self.adjustment_type == "ChangeInCapacity":
            self.autoscaling_backend.change_capacity(
                self.auto_scaling_group_name, self.scaling_adjustment
            )
        elif self.adjustment_type == "PercentChangeInCapacity":
            self.autoscaling_backend.change_capacity_percent(
                self.auto_scaling_group_name, self.scaling_adjustment
            )


class FakeLaunchConfiguration(CloudFormationModel):
    def __init__(
        self,
        name: str,
        image_id: str,
        key_name: Optional[str],
        ramdisk_id: str,
        kernel_id: str,
        security_groups: list[str],
        user_data: Optional[Base64EncodedString],
        instance_type: str,
        instance_monitoring: bool,
        instance_profile_name: Optional[str],
        spot_price: Optional[str],
        ebs_optimized: bool,
        associate_public_ip_address: bool,
        block_device_mapping_dict: list[dict[str, Any]],
        account_id: str,
        region_name: str,
        metadata_options: Optional[str],
        classic_link_vpc_id: Optional[str],
        classic_link_vpc_security_groups: Optional[str],
    ):
        self.name = name
        self.image_id = image_id
        self.key_name = key_name
        self.ramdisk_id = ramdisk_id
        self.kernel_id = kernel_id
        self.security_groups = security_groups if security_groups else []
        self.user_data = user_data
        self.instance_type = instance_type
        self.instance_monitoring_enabled = instance_monitoring
        self.iam_instance_profile = instance_profile_name
        self.spot_price = spot_price
        self.ebs_optimized = ebs_optimized
        self.associate_public_ip_address = associate_public_ip_address
        self.block_device_mapping_dict = block_device_mapping_dict
        self.metadata_options = metadata_options
        self.classic_link_vpc_id = classic_link_vpc_id
        self.classic_link_vpc_security_groups = classic_link_vpc_security_groups
        self.arn = f"arn:{get_partition(region_name)}:autoscaling:{region_name}:{account_id}:launchConfiguration:9dbbbf87-6141-428a-a409-0752edbe6cad:launchConfigurationName/{self.name}"
        self.created_time = utcnow()

    @classmethod
    def create_from_instance(
        cls, name: str, instance: Instance, backend: AutoScalingBackend
    ) -> FakeLaunchConfiguration:
        security_group_names = [sg.name for sg in instance.security_groups]
        config = backend.create_launch_configuration(
            name=name,
            image_id=instance.image_id,
            kernel_id="",
            ramdisk_id="",
            key_name=instance.key_name,
            security_groups=security_group_names,
            user_data=instance.user_data,
            instance_type=instance.instance_type,
            instance_monitoring=False,
            instance_profile_name=None,
            spot_price=None,
            ebs_optimized=instance.ebs_optimized,
            associate_public_ip_address=instance.associate_public_ip,
            # We expect a dictionary in the same format as when the user calls it
            block_device_mappings=instance.block_device_mapping.to_source_dict(),
        )
        return config

    @staticmethod
    def cloudformation_name_type() -> str:
        return "LaunchConfigurationName"

    @staticmethod
    def cloudformation_type() -> str:
        # https://docs.aws.amazon.com/AWSCloudFormation/latest/UserGuide/aws-resource-autoscaling-launchconfiguration.html
        return "AWS::AutoScaling::LaunchConfiguration"

    @classmethod
    def create_from_cloudformation_json(  # type: ignore[misc]
        cls,
        resource_name: str,
        cloudformation_json: Any,
        account_id: str,
        region_name: str,
        **kwargs: Any,
    ) -> FakeLaunchConfiguration:
        properties = cloudformation_json["Properties"]

        instance_profile_name = properties.get("IamInstanceProfile")

        backend = autoscaling_backends[account_id][region_name]
        config = backend.create_launch_configuration(
            name=resource_name,
            image_id=properties.get("ImageId"),
            kernel_id=properties.get("KernelId"),
            ramdisk_id=properties.get("RamdiskId"),
            key_name=properties.get("KeyName"),
            security_groups=properties.get("SecurityGroups"),
            user_data=properties.get("UserData"),
            instance_type=properties.get("InstanceType"),
            instance_monitoring=properties.get("InstanceMonitoring"),
            instance_profile_name=instance_profile_name,
            spot_price=properties.get("SpotPrice"),
            ebs_optimized=properties.get("EbsOptimized"),
            associate_public_ip_address=properties.get("AssociatePublicIpAddress"),
            block_device_mappings=properties.get("BlockDeviceMapping.member"),
        )
        return config

    @classmethod
    def update_from_cloudformation_json(  # type: ignore[misc]
        cls,
        original_resource: Any,
        new_resource_name: str,
        cloudformation_json: Any,
        account_id: str,
        region_name: str,
    ) -> FakeLaunchConfiguration:
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
        cloudformation_json: Any,
        account_id: str,
        region_name: str,
    ) -> None:
        backend = autoscaling_backends[account_id][region_name]
        try:
            backend.delete_launch_configuration(resource_name)
        except KeyError:
            pass

    def delete(self, account_id: str, region_name: str) -> None:
        backend = autoscaling_backends[account_id][region_name]
        backend.delete_launch_configuration(self.name)

    @property
    def physical_resource_id(self) -> str:
        return self.name

    @property
    def block_device_mappings(self) -> list[dict[str, Any]]:
        if not self.block_device_mapping_dict:
            return []
        parsed = self._parse_block_device_mappings()
        value = [
            {
                "VirtualName": mapping.ephemeral_name,
                "DeviceName": mount_point,
                "Ebs": {
                    "SnapshotId": mapping.snapshot_id,
                    "VolumeSize": mapping.size,
                    "VolumeType": mapping.volume_type,
                    "DeleteOnTermination": mapping.delete_on_termination,
                    "Iops": mapping.iops,
                    "Encrypted": mapping.encrypted,
                    "Throughput": mapping.throughput,
                },
                "NoDevice": mapping.no_device,
            }
            for mount_point, mapping in parsed.items()
        ]
        return value

    @property
    def instance_monitoring(self) -> dict[str, bool]:
        return {"Enabled": self.instance_monitoring_enabled}

    def _parse_block_device_mappings(self) -> BlockDeviceMapping:
        block_device_map = BlockDeviceMapping()
        for mapping in self.block_device_mapping_dict:
            block_type = BlockDeviceType()
            mount_point = mapping.get("DeviceName")
            if mapping.get("VirtualName") and "ephemeral" in mapping.get("VirtualName"):  # type: ignore[operator]
                block_type.ephemeral_name = mapping.get("VirtualName")
            elif mapping.get("NoDevice", "false") == "true":
                block_type.no_device = "true"
            else:
                ebs = mapping.get("Ebs", {})
                block_type.volume_type = ebs.get("VolumeType")
                block_type.snapshot_id = ebs.get("SnapshotId")
                block_type.delete_on_termination = ebs.get("DeleteOnTermination")
                block_type.size = ebs.get("VolumeSize")
                block_type.iops = ebs.get("Iops")
                block_type.throughput = ebs.get("Throughput")
                block_type.encrypted = ebs.get("Encrypted")
            block_device_map[mount_point] = block_type
        return block_device_map


class FakeScheduledAction(CloudFormationModel):
    def __init__(
        self,
        autos_caling_group_name: str,
        desired_capacity: Optional[int],
        max_size: Optional[int],
        min_size: Optional[int],
        scheduled_action_name: str,
        start_time: Optional[str],
        end_time: Optional[str],
        recurrence: Optional[str],
        time_zone: Optional[str],
    ):
        self.auto_scaling_group_name = autos_caling_group_name
        self.desired_capacity = desired_capacity
        self.max_size = max_size
        self.min_size = min_size
        self.start_time = start_time
        self.end_time = end_time
        self.recurrence = recurrence
        self.scheduled_action_name = scheduled_action_name
        self.time_zone = time_zone

    @staticmethod
    def cloudformation_name_type() -> str:
        return "ScheduledActionName"

    @staticmethod
    def cloudformation_type() -> str:
        # https://docs.aws.amazon.com/AWSCloudFormation/latest/UserGuide/aws-resource-as-scheduledaction.html
        return "AWS::AutoScaling::ScheduledAction"

    @classmethod
    def create_from_cloudformation_json(  # type: ignore[misc]
        cls,
        resource_name: str,
        cloudformation_json: dict[str, Any],
        account_id: str,
        region_name: str,
        **kwargs: Any,
    ) -> FakeScheduledAction:
        properties = cloudformation_json["Properties"]

        backend = autoscaling_backends[account_id][region_name]

        scheduled_action_name = (
            kwargs["LogicalId"]
            if kwargs.get("LogicalId")
            else "ScheduledScalingAction-{random.randint(0,100)}"
        )

        scheduled_action = backend.put_scheduled_update_group_action(
            name=properties.get("AutoScalingGroupName"),
            desired_capacity=properties.get("DesiredCapacity"),
            max_size=properties.get("MaxSize"),
            min_size=properties.get("MinSize"),
            scheduled_action_name=scheduled_action_name,
            start_time=properties.get("StartTime"),
            end_time=properties.get("EndTime"),
            recurrence=properties.get("Recurrence"),
            timezone=properties.get("TimeZone"),
        )
        return scheduled_action


class FailedScheduledUpdateGroupActionRequest:
    def __init__(
        self,
        *,
        scheduled_action_name: str,
        error_code: Optional[str] = None,
        error_message: Optional[str] = None,
    ) -> None:
        self.scheduled_action_name = scheduled_action_name
        self.error_code = error_code
        self.error_message = error_message


class FakeWarmPool(CloudFormationModel):
    def __init__(
        self,
        max_group_prepared_capacity: Optional[int],
        min_size: Optional[int],
        pool_state: Optional[str],
        instance_reuse_policy: Optional[dict[str, bool]],
    ):
        self.max_group_prepared_capacity = max_group_prepared_capacity
        self.min_size = min_size or 0
        self.pool_state = pool_state or "Stopped"
        self.instance_reuse_policy = instance_reuse_policy


class FakeAutoScalingGroup(CloudFormationModel):
    def __init__(
        self,
        name: str,
        availability_zones: list[str],
        desired_capacity: Optional[int],
        max_size: Optional[int],
        min_size: Optional[int],
        launch_config_name: str,
        launch_template: dict[str, Any],
        vpc_zone_identifier: Optional[str],
        default_cooldown: Optional[int],
        health_check_period: Optional[int],
        health_check_type: Optional[str],
        load_balancers: list[str],
        target_group_arns: list[str],
        placement_group: Optional[str],
        termination_policies: list[str],
        autoscaling_backend: AutoScalingBackend,
        ec2_backend: EC2Backend,
        tags: list[dict[str, str]],
        mixed_instances_policy: Optional[dict[str, Any]],
        capacity_rebalance: bool,
        new_instances_protected_from_scale_in: bool = False,
    ):
        self.autoscaling_backend = autoscaling_backend
        self.ec2_backend = ec2_backend
        self.name = name
        self._id = str(random.uuid4())
        self.region = self.autoscaling_backend.region_name
        self.account_id = self.autoscaling_backend.account_id
        partition = get_partition(self.region)
        self.service_linked_role_arn = f"arn:{partition}:iam::{self.account_id}:role/aws-service-role/autoscaling.amazonaws.com/AWSServiceRoleForAutoScaling"

        self.vpc_zone_identifier: Optional[str] = None
        self._set_azs_and_vpcs(availability_zones, vpc_zone_identifier)

        self.max_size = max_size
        self.min_size = min_size

        self.mixed_instances_policy = mixed_instances_policy
        self.ec2_launch_template: Optional[LaunchTemplate] = None
        # Will be None if self.launch_template is used instead
        self.launch_config: FakeLaunchConfiguration = None  # type: ignore[assignment]

        # Some defaults, if not set
        if (
            self.mixed_instances_policy
            and "InstancesDistribution" not in self.mixed_instances_policy
        ):
            self.mixed_instances_policy["InstancesDistribution"] = {
                "OnDemandAllocationStrategy": "prioritized",
                "OnDemandBaseCapacity": 0,
                "OnDemandPercentageAboveBaseCapacity": 100,
                "SpotAllocationStrategy": "lowest-price",
                "SpotInstancePools": 2,
            }

        self._set_launch_configuration(
            launch_config_name, launch_template, mixed_instances_policy
        )

        self.default_cooldown = (
            default_cooldown if default_cooldown else DEFAULT_COOLDOWN
        )
        self.health_check_grace_period = health_check_period
        self.health_check_type = health_check_type if health_check_type else "EC2"
        self.load_balancer_names = load_balancers
        self.target_group_arns = target_group_arns
        self.placement_group = placement_group
        self.capacity_rebalance = capacity_rebalance
        self.termination_policies = termination_policies or ["Default"]
        self.new_instances_protected_from_scale_in = (
            new_instances_protected_from_scale_in
        )

        self.suspended_processes = []
        self.instance_states: list[InstanceState] = []
        self.tags: list[dict[str, str]] = tags or []
        self.set_desired_capacity(desired_capacity)

        self.metrics: list[str] = []
        self.warm_pool: Optional[FakeWarmPool] = None
        self.created_time = datetime.now().isoformat()

    @property
    def launch_template(self) -> Optional[dict[str, Any]]:
        if self.ec2_launch_template is None:
            return None
        lt = {
            "LaunchTemplateId": self.ec2_launch_template.id,
            "LaunchTemplateName": self.ec2_launch_template.name,
            "Version": self.provided_launch_template_version,
        }
        return lt

    @property
    def enabled_metrics(self) -> list[dict[str, str]]:
        return [{"Metric": metric, "Granularity": "1Minute"} for metric in self.metrics]

    @property
    def suspended_processes(self) -> list[dict[str, str]]:
        return [
            {"ProcessName": process, "SuspensionReason": ""}
            for process in self._suspended_processes
        ]

    @suspended_processes.setter
    def suspended_processes(self, processes: list[str]) -> None:
        self._suspended_processes = processes

    @property
    def tags(self) -> list[dict[str, str]]:
        return self._tags

    @tags.setter
    def tags(self, tags: list[dict[str, Any]]) -> None:
        for tag in tags:
            if "ResourceId" not in tag or not tag["ResourceId"]:
                tag["ResourceId"] = self.name
            if "ResourceType" not in tag or not tag["ResourceType"]:
                tag["ResourceType"] = "auto-scaling-group"
            if "PropagateAtLaunch" not in tag:
                tag["PropagateAtLaunch"] = False
        self._tags = tags

    @property
    def arn(self) -> str:
        return f"arn:{get_partition(self.region)}:autoscaling:{self.region}:{self.account_id}:autoScalingGroup:{self._id}:autoScalingGroupName/{self.name}"

    def active_instances(self) -> list[InstanceState]:
        return [x for x in self.instance_states if x.lifecycle_state == "InService"]

    def _set_azs_and_vpcs(
        self,
        availability_zones: list[str],
        vpc_zone_identifier: Optional[str],
        update: bool = False,
    ) -> None:
        # for updates, if only AZs are provided, they must not clash with
        # the AZs of existing VPCs
        if update and availability_zones and not vpc_zone_identifier:
            vpc_zone_identifier = self.vpc_zone_identifier

        if vpc_zone_identifier:
            # extract azs for vpcs
            subnet_ids = vpc_zone_identifier.split(",")
            subnets = self.autoscaling_backend.ec2_backend.describe_subnets(
                subnet_ids=subnet_ids
            )
            vpc_zones = [subnet.availability_zone for subnet in subnets]

            if availability_zones and set(availability_zones) != set(vpc_zones):
                raise AutoscalingClientError(
                    "ValidationError",
                    "The availability zones of the specified subnets and the Auto Scaling group do not match",
                )
            availability_zones = vpc_zones
        elif not availability_zones:
            if not update:
                raise AutoscalingClientError(
                    "ValidationError",
                    "At least one Availability Zone or VPC Subnet is required.",
                )
            return

        self.availability_zones = availability_zones
        self.vpc_zone_identifier = vpc_zone_identifier

    def _set_launch_configuration(
        self,
        launch_config_name: str,
        launch_template: dict[str, Any],
        mixed_instances_policy: Optional[dict[str, Any]],
    ) -> None:
        if launch_config_name:
            self.launch_config = self.autoscaling_backend.launch_configurations[
                launch_config_name
            ]
            self.launch_configuration_name = launch_config_name

        if launch_template or mixed_instances_policy:
            if launch_template:
                launch_template_id = launch_template.get("LaunchTemplateId")
                launch_template_name = launch_template.get("LaunchTemplateName")
                # If no version is specified, AWS will use '$Default'
                # However, AWS will never show the version if it is not specified
                # (If the user explicitly specifies '$Default', it will be returned)
                self.launch_template_version = (
                    launch_template.get("Version") or "$Default"
                )
                self.provided_launch_template_version = launch_template.get("Version")
            elif mixed_instances_policy:
                spec = mixed_instances_policy["LaunchTemplate"][
                    "LaunchTemplateSpecification"
                ]
                launch_template_id = spec.get("LaunchTemplateId")
                launch_template_name = spec.get("LaunchTemplateName")
                self.launch_template_version = spec.get("Version") or "$Default"

            if not (launch_template_id or launch_template_name) or (
                launch_template_id and launch_template_name
            ):
                raise ValidationError(
                    "Valid requests must contain either launchTemplateId or LaunchTemplateName"
                )

            if launch_template_id:
                self.ec2_launch_template = self.ec2_backend.get_launch_template(
                    launch_template_id
                )
            elif launch_template_name:
                self.ec2_launch_template = self.ec2_backend.get_launch_template_by_name(
                    launch_template_name
                )

            # This logic was in the original XML template and may need to be rethought.
            if (
                self.ec2_launch_template is not None
                and self.mixed_instances_policy is not None
            ):
                try:
                    self.mixed_instances_policy["LaunchTemplate"][
                        "LaunchTemplateSpecification"
                    ] = {
                        "LaunchTemplateId": self.ec2_launch_template.id,
                        "LaunchTemplateName": self.ec2_launch_template.name,
                        "Version": self.launch_template_version,
                    }
                except (AttributeError, KeyError, TypeError):
                    pass

            try:
                if (
                    self.mixed_instances_policy
                    and "Overrides" not in self.mixed_instances_policy["LaunchTemplate"]
                ):
                    self.mixed_instances_policy["LaunchTemplate"]["Overrides"] = []
            except (AttributeError, KeyError, TypeError):
                pass

    @staticmethod
    def cloudformation_name_type() -> str:
        return "AutoScalingGroupName"

    @staticmethod
    def cloudformation_type() -> str:
        # https://docs.aws.amazon.com/AWSCloudFormation/latest/UserGuide/aws-resource-autoscaling-autoscalinggroup.html
        return "AWS::AutoScaling::AutoScalingGroup"

    @classmethod
    def create_from_cloudformation_json(  # type: ignore[misc]
        cls,
        resource_name: str,
        cloudformation_json: dict[str, Any],
        account_id: str,
        region_name: str,
        **kwargs: Any,
    ) -> FakeAutoScalingGroup:
        properties = cloudformation_json["Properties"]

        launch_config_name = properties.get("LaunchConfigurationName")
        launch_template = properties.get("LaunchTemplate", {})
        load_balancer_names = properties.get("LoadBalancerNames", [])
        target_group_arns = properties.get("TargetGroupARNs", [])
        mixed_instances_policy = properties.get("MixedInstancesPolicy", {})

        backend: AutoScalingBackend = autoscaling_backends[account_id][region_name]
        group = backend.create_auto_scaling_group(
            name=resource_name,
            availability_zones=properties.get("AvailabilityZones", []),
            desired_capacity=properties.get("DesiredCapacity"),
            max_size=properties.get("MaxSize"),
            min_size=properties.get("MinSize"),
            launch_config_name=launch_config_name,
            launch_template=launch_template,
            vpc_zone_identifier=(
                ",".join(properties.get("VPCZoneIdentifier", [])) or None
            ),
            default_cooldown=properties.get("Cooldown"),
            health_check_period=properties.get("HealthCheckGracePeriod"),
            health_check_type=properties.get("HealthCheckType"),
            load_balancers=load_balancer_names,
            target_group_arns=target_group_arns,
            placement_group=None,
            termination_policies=properties.get("TerminationPolicies", []),
            tags=properties.get("Tags", []),
            new_instances_protected_from_scale_in=properties.get(
                "NewInstancesProtectedFromScaleIn", False
            ),
            mixed_instances_policy=mixed_instances_policy,
        )
        return group

    @classmethod
    def update_from_cloudformation_json(  # type: ignore[misc]
        cls,
        original_resource: Any,
        new_resource_name: str,
        cloudformation_json: dict[str, Any],
        account_id: str,
        region_name: str,
    ) -> FakeAutoScalingGroup:
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
        backend = autoscaling_backends[account_id][region_name]
        try:
            backend.delete_auto_scaling_group(resource_name)
        except KeyError:
            pass

    def delete(self, account_id: str, region_name: str) -> None:
        backend = autoscaling_backends[account_id][region_name]
        backend.delete_auto_scaling_group(self.name)

    @property
    def physical_resource_id(self) -> str:
        return self.name

    @property
    def image_id(self) -> str:
        if self.ec2_launch_template:
            version = self.ec2_launch_template.get_version(self.launch_template_version)
            return version.image_id

        return self.launch_config.image_id  # type: ignore[union-attr]

    @property
    def instance_type(self) -> str:
        if self.ec2_launch_template:
            version = self.ec2_launch_template.get_version(self.launch_template_version)
            return version.instance_type

        return self.launch_config.instance_type  # type: ignore[union-attr]

    @property
    def user_data(self) -> Optional[Base64EncodedString]:
        if self.ec2_launch_template:
            version = self.ec2_launch_template.get_version(self.launch_template_version)
            return version.user_data

        return self.launch_config.user_data

    @property
    def security_groups(self) -> list[str]:
        if self.ec2_launch_template:
            version = self.ec2_launch_template.get_version(self.launch_template_version)
            return version.security_groups

        return self.launch_config.security_groups  # type: ignore[union-attr]

    @property
    def instance_tags(self) -> dict[str, str]:
        if self.ec2_launch_template:
            version = self.ec2_launch_template.get_version(self.launch_template_version)
            return version.instance_tags
        return {}

    @property
    def instances(self) -> list[InstanceState]:
        return self.instance_states

    @property
    def warm_pool_configuration(self) -> Optional[FakeWarmPool]:
        return self.warm_pool

    def update(
        self,
        availability_zones: list[str],
        desired_capacity: Optional[int],
        max_size: Optional[int],
        min_size: Optional[int],
        launch_config_name: str,
        launch_template: dict[str, Any],
        vpc_zone_identifier: str,
        health_check_period: int,
        health_check_type: str,
        new_instances_protected_from_scale_in: Optional[bool] = None,
        mixed_instances_policy: Optional[dict[str, Any]] = None,
    ) -> None:
        self._set_azs_and_vpcs(availability_zones, vpc_zone_identifier, update=True)

        if max_size is not None:
            self.max_size = max_size
        if min_size is not None:
            self.min_size = min_size

        if desired_capacity is None:
            if min_size is not None and min_size > len(self.instance_states):
                desired_capacity = min_size
            if max_size is not None and max_size < len(self.instance_states):
                desired_capacity = max_size

        self.mixed_instances_policy = mixed_instances_policy
        self._set_launch_configuration(
            launch_config_name,
            launch_template,
            mixed_instances_policy=mixed_instances_policy,
        )

        if health_check_period is not None:
            self.health_check_grace_period = health_check_period
        if health_check_type is not None:
            self.health_check_type = health_check_type
        if new_instances_protected_from_scale_in is not None:
            self.new_instances_protected_from_scale_in = (
                new_instances_protected_from_scale_in
            )

        if desired_capacity is not None:
            self.set_desired_capacity(desired_capacity)

    def set_desired_capacity(self, new_capacity: Optional[int]) -> None:
        if new_capacity is None:
            self.desired_capacity = self.min_size
        else:
            self.desired_capacity = new_capacity

        current_instance_count = len(self.active_instances())

        is_mixed_instances = self.mixed_instances_policy and len(
            self.mixed_instances_policy.get("LaunchTemplate", {}).get("Overrides", [])
        )
        target_capacity = self.desired_capacity or 0

        if is_mixed_instances:
            policy = self.mixed_instances_policy or {}
            # These calculations assume a strategy of "prioritized"
            overrides = policy["LaunchTemplate"]["Overrides"]
            distribution = policy.get("InstancesDistribution", {})

            on_demand_base = int(distribution.get("OnDemandBaseCapacity", 0))
            percent_above_base = int(
                distribution.get("OnDemandPercentageAboveBaseCapacity", 100)
            )

            # When using a "prioritized" strategy, AWS will treat the overrides as priority list when deciding
            # which instances to launch, meaning we always pick the first entry in a mocked environment.
            primary_weight_str = overrides[0].get("WeightedCapacity", "1")
            primary_weight = (
                int(primary_weight_str) if primary_weight_str.isdigit() else 1
            )
            if primary_weight == 0:
                primary_weight = 1

            if on_demand_base >= target_capacity:
                # If the base capacity meets or exceeds desired capacity, the entire desired capacity is fulfilled by On-Demand.
                total_on_demand_capacity = target_capacity
                total_spot_capacity = 0

            else:
                # After fulfilling the OnDemandBase, we need to add more on-demand and spot instances according
                # to the passed percentage.
                above_base_capacity = target_capacity - on_demand_base

                on_demand_above_base_capacity = (
                    int(above_base_capacity * percent_above_base) / 100.0
                )
                spot_above_base_capacity = int(
                    above_base_capacity - on_demand_above_base_capacity
                )

                total_on_demand_capacity = int(
                    on_demand_base + on_demand_above_base_capacity
                )
                total_spot_capacity = spot_above_base_capacity

            on_demand_instances = math.ceil(total_on_demand_capacity / primary_weight)
            spot_instances = math.ceil(total_spot_capacity / primary_weight)

            total_target_instances = on_demand_instances + spot_instances

        else:
            total_target_instances = target_capacity

        instance_count_delta = total_target_instances - current_instance_count

        if instance_count_delta == 0:
            pass
        elif instance_count_delta > 0:
            count_needed = instance_count_delta
            propagated_tags = self.get_propagated_tags()
            self.replace_autoscaling_group_instances(count_needed, propagated_tags)
        else:
            count_to_remove = abs(instance_count_delta)

            instances_to_remove = [  # only remove unprotected
                state
                for state in self.instance_states
                if not state.protected_from_scale_in
            ][:count_to_remove]

            if instances_to_remove:  # just in case not instances to remove
                instance_ids_to_remove = [
                    instance.instance.id for instance in instances_to_remove
                ]
                self.autoscaling_backend.ec2_backend.terminate_instances(
                    instance_ids_to_remove
                )
                self.instance_states = list(
                    set(self.instance_states) - set(instances_to_remove)
                )

        if self.name in self.autoscaling_backend.autoscaling_groups:
            self.autoscaling_backend.update_attached_elbs(self.name)
            self.autoscaling_backend.update_attached_target_groups(self.name)

    def get_propagated_tags(self) -> dict[str, str]:
        propagated_tags = {}
        for tag in self.tags:
            if tag.get("PropagateAtLaunch"):
                propagated_tags[tag["Key"]] = tag["Value"]
        return propagated_tags

    def replace_autoscaling_group_instances(
        self, count_needed: int, propagated_tags: dict[str, str]
    ) -> None:
        propagated_tags[ASG_NAME_TAG] = self.name
        propagated_tags.update(self.instance_tags)

        # VPCZoneIdentifier:
        # A comma-separated list of subnet IDs for a virtual private cloud (VPC) where instances in the Auto Scaling group can be created.
        # We'll create all instances in a single subnet to make things easier
        subnet_id = (
            self.vpc_zone_identifier.split(",")[0] if self.vpc_zone_identifier else None
        )
        associate_public_ip = (
            self.launch_config.associate_public_ip_address
            if self.launch_config
            else None
        )
        launch_template = None
        if self.ec2_launch_template:
            if self.ec2_launch_template.id:
                launch_template = {"LaunchTemplateId": self.ec2_launch_template.id}
            elif self.ec2_launch_template.name:
                launch_template = {"LaunchTemplateName": self.ec2_launch_template.name}
        reservation = self.autoscaling_backend.ec2_backend.run_instances(
            self.image_id,
            count_needed,
            self.user_data,
            self.security_groups,
            instance_type=self.instance_type,
            tags={"instance": propagated_tags},
            placement=random.choice(self.availability_zones),
            launch_config=self.launch_config,
            launch_template=launch_template,
            is_instance_type_default=False,
            associate_public_ip=associate_public_ip,
            subnet_id=subnet_id,
        )
        for instance in reservation.instances:
            instance.autoscaling_group = self
            self.instance_states.append(
                InstanceState(
                    instance,
                    protected_from_scale_in=self.new_instances_protected_from_scale_in,
                )
            )

    def append_target_groups(self, target_group_arns: list[str]) -> None:
        append = [x for x in target_group_arns if x not in self.target_group_arns]
        self.target_group_arns.extend(append)

    def enable_metrics_collection(self, metrics: list[str]) -> None:
        self.metrics = metrics or []

    def put_warm_pool(
        self,
        max_group_prepared_capacity: Optional[int],
        min_size: Optional[int],
        pool_state: Optional[str],
        instance_reuse_policy: Optional[dict[str, bool]],
    ) -> None:
        self.warm_pool = FakeWarmPool(
            max_group_prepared_capacity=max_group_prepared_capacity,
            min_size=min_size,
            pool_state=pool_state,
            instance_reuse_policy=instance_reuse_policy,
        )

    def get_warm_pool(self) -> Optional[FakeWarmPool]:
        return self.warm_pool


class AutoScalingBackend(BaseBackend):
    def __init__(self, region_name: str, account_id: str):
        super().__init__(region_name, account_id)
        self.autoscaling_groups: dict[str, FakeAutoScalingGroup] = OrderedDict()
        self.launch_configurations: dict[str, FakeLaunchConfiguration] = OrderedDict()
        self.scheduled_actions: dict[str, FakeScheduledAction] = OrderedDict()
        self.policies: dict[str, FakeScalingPolicy] = {}
        self.lifecycle_hooks: dict[str, LifecycleHook] = {}
        self.ec2_backend: EC2Backend = ec2_backends[self.account_id][region_name]
        self.elb_backend: ELBBackend = elb_backends[self.account_id][region_name]
        self.elbv2_backend: ELBv2Backend = elbv2_backends[self.account_id][region_name]

    def create_launch_configuration(
        self,
        name: str,
        image_id: str,
        key_name: Optional[str],
        kernel_id: str,
        ramdisk_id: str,
        security_groups: list[str],
        user_data: Optional[Base64EncodedString],
        instance_type: str,
        instance_monitoring: bool,
        instance_profile_name: Optional[str],
        spot_price: Optional[str],
        ebs_optimized: bool,
        associate_public_ip_address: bool,
        block_device_mappings: list[dict[str, Any]],
        instance_id: Optional[str] = None,
        metadata_options: Optional[str] = None,
        classic_link_vpc_id: Optional[str] = None,
        classic_link_vpc_security_groups: Optional[str] = None,
    ) -> FakeLaunchConfiguration:
        valid_requests = [
            instance_id is not None,
            image_id is not None and instance_type is not None,
        ]
        if not any(valid_requests):
            raise ValidationError(
                "Valid requests must contain either the InstanceID parameter or both the ImageId and InstanceType parameters."
            )
        if instance_id is not None:
            # TODO: https://docs.aws.amazon.com/autoscaling/ec2/userguide/create-lc-with-instanceID.html
            pass
        launch_configuration = FakeLaunchConfiguration(
            name=name,
            image_id=image_id,
            key_name=key_name,
            kernel_id=kernel_id,
            ramdisk_id=ramdisk_id,
            security_groups=security_groups,
            user_data=user_data,
            instance_type=instance_type,
            instance_monitoring=instance_monitoring,
            instance_profile_name=instance_profile_name,
            spot_price=spot_price,
            ebs_optimized=ebs_optimized,
            associate_public_ip_address=associate_public_ip_address,
            block_device_mapping_dict=block_device_mappings,
            account_id=self.account_id,
            region_name=self.region_name,
            metadata_options=metadata_options,
            classic_link_vpc_id=classic_link_vpc_id,
            classic_link_vpc_security_groups=classic_link_vpc_security_groups,
        )
        self.launch_configurations[name] = launch_configuration
        return launch_configuration

    def describe_launch_configurations(
        self, names: Optional[list[str]]
    ) -> list[FakeLaunchConfiguration]:
        configurations = self.launch_configurations.values()
        if names:
            return [
                configuration
                for configuration in configurations
                if configuration.name in names
            ]
        else:
            return list(configurations)

    def delete_launch_configuration(self, launch_configuration_name: str) -> None:
        self.launch_configurations.pop(launch_configuration_name, None)

    def put_scheduled_update_group_action(
        self,
        name: str,
        desired_capacity: Union[None, str, int],
        max_size: Union[None, str, int],
        min_size: Union[None, str, int],
        scheduled_action_name: str,
        start_time: Optional[str],
        end_time: Optional[str],
        recurrence: Optional[str],
        timezone: Optional[str],
    ) -> FakeScheduledAction:
        max_size = make_int(max_size)
        min_size = make_int(min_size)
        desired_capacity = make_int(desired_capacity)

        scheduled_action = FakeScheduledAction(
            autos_caling_group_name=name,
            desired_capacity=desired_capacity,
            max_size=max_size,
            min_size=min_size,
            scheduled_action_name=scheduled_action_name,
            start_time=start_time,
            end_time=end_time,
            recurrence=recurrence,
            time_zone=timezone,
        )

        self.scheduled_actions[scheduled_action_name] = scheduled_action
        return scheduled_action

    def batch_put_scheduled_update_group_action(
        self, name: str, actions: list[dict[str, Any]]
    ) -> list[FailedScheduledUpdateGroupActionRequest]:
        result = []
        for action in actions:
            try:
                self.put_scheduled_update_group_action(
                    name=name,
                    desired_capacity=action.get("DesiredCapacity"),
                    max_size=action.get("MaxSize"),
                    min_size=action.get("MinSize"),
                    scheduled_action_name=action["ScheduledActionName"],
                    start_time=action.get("StartTime"),
                    end_time=action.get("EndTime"),
                    recurrence=action.get("Recurrence"),
                    timezone=action.get("TimeZone"),
                )
            except AutoscalingClientError as err:
                result.append(
                    FailedScheduledUpdateGroupActionRequest(
                        scheduled_action_name=action["ScheduledActionName"],
                        error_code=err.code,
                        error_message=err.message,
                    )
                )
        return result

    def describe_scheduled_actions(
        self,
        autoscaling_group_name: Optional[str] = None,
        scheduled_action_names: Optional[list[str]] = None,
    ) -> list[FakeScheduledAction]:
        scheduled_actions = []
        for scheduled_action in self.scheduled_actions.values():
            if (
                not autoscaling_group_name
                or scheduled_action.auto_scaling_group_name == autoscaling_group_name
            ):
                if (
                    not scheduled_action_names
                    or scheduled_action.scheduled_action_name in scheduled_action_names
                ):
                    scheduled_actions.append(scheduled_action)

        return scheduled_actions

    def delete_scheduled_action(
        self, auto_scaling_group_name: str, scheduled_action_name: str
    ) -> None:
        scheduled_action = self.describe_scheduled_actions(
            auto_scaling_group_name, [scheduled_action_name]
        )
        if scheduled_action:
            self.scheduled_actions.pop(scheduled_action_name, None)
        else:
            raise ValidationError("Scheduled action name not found")

    def batch_delete_scheduled_action(
        self, auto_scaling_group_name: str, scheduled_action_names: list[str]
    ) -> list[FailedScheduledUpdateGroupActionRequest]:
        result = []
        for scheduled_action_name in scheduled_action_names:
            try:
                self.delete_scheduled_action(
                    auto_scaling_group_name, scheduled_action_name
                )
            except AutoscalingClientError as err:
                result.append(
                    FailedScheduledUpdateGroupActionRequest(
                        scheduled_action_name=scheduled_action_name,
                        error_code=err.code,
                        error_message=err.message,
                    )
                )
        return result

    def create_auto_scaling_group(
        self,
        name: str,
        availability_zones: list[str],
        desired_capacity: Union[None, str, int],
        max_size: Union[None, str, int],
        min_size: Union[None, str, int],
        launch_config_name: str,
        launch_template: dict[str, Any],
        vpc_zone_identifier: Optional[str],
        default_cooldown: Optional[int],
        health_check_period: Union[None, str, int],
        health_check_type: Optional[str],
        load_balancers: list[str],
        target_group_arns: list[str],
        placement_group: Optional[str],
        termination_policies: list[str],
        tags: list[dict[str, str]],
        capacity_rebalance: bool = False,
        new_instances_protected_from_scale_in: bool = False,
        instance_id: Optional[str] = None,
        mixed_instances_policy: Optional[dict[str, Any]] = None,
    ) -> FakeAutoScalingGroup:
        max_size = make_int(max_size)
        min_size = make_int(min_size)
        desired_capacity = make_int(desired_capacity)
        default_cooldown = make_int(default_cooldown)

        # Verify only a single launch config-like parameter is provided.
        params = [
            launch_config_name,
            launch_template,
            instance_id,
            mixed_instances_policy,
        ]
        num_params = sum([1 for param in params if param])

        if num_params != 1:
            raise ValidationError(
                "Valid requests must contain either LaunchTemplate, LaunchConfigurationName, "
                "InstanceId or MixedInstancesPolicy parameter."
            )

        if instance_id:
            try:
                instance = self.ec2_backend.get_instance(instance_id)
                launch_config_name = name
                FakeLaunchConfiguration.create_from_instance(
                    launch_config_name, instance, self
                )
            except InvalidInstanceIdError:
                raise InvalidInstanceError(instance_id)

        group = FakeAutoScalingGroup(
            name=name,
            availability_zones=availability_zones,
            desired_capacity=desired_capacity,
            max_size=max_size,
            min_size=min_size,
            launch_config_name=launch_config_name,
            launch_template=launch_template,
            vpc_zone_identifier=vpc_zone_identifier,
            default_cooldown=default_cooldown,
            health_check_period=make_int(health_check_period or 300),
            health_check_type=health_check_type,
            load_balancers=load_balancers,
            target_group_arns=target_group_arns,
            placement_group=placement_group,
            termination_policies=termination_policies,
            autoscaling_backend=self,
            ec2_backend=self.ec2_backend,
            tags=tags,
            new_instances_protected_from_scale_in=new_instances_protected_from_scale_in,
            mixed_instances_policy=mixed_instances_policy,
            capacity_rebalance=capacity_rebalance,
        )

        self.autoscaling_groups[name] = group
        self.update_attached_elbs(group.name)
        self.update_attached_target_groups(group.name)
        return group

    def update_auto_scaling_group(
        self,
        name: str,
        availability_zones: list[str],
        desired_capacity: Optional[int],
        max_size: Optional[int],
        min_size: Optional[int],
        launch_config_name: str,
        launch_template: dict[str, Any],
        vpc_zone_identifier: str,
        health_check_period: int,
        health_check_type: str,
        new_instances_protected_from_scale_in: Optional[bool] = None,
        mixed_instances_policy: Optional[dict[str, Any]] = None,
    ) -> FakeAutoScalingGroup:
        """
        The parameter DefaultCooldown, PlacementGroup, TerminationPolicies are not yet implemented
        """
        # Verify only a single launch config-like parameter is provided.
        if launch_config_name and launch_template:
            raise ValidationError(
                "Valid requests must contain either LaunchTemplate, LaunchConfigurationName "
                "or MixedInstancesPolicy parameter."
            )
        if name not in self.autoscaling_groups:
            raise ValidationError("AutoScalingGroup name not found - null")

        group = self.autoscaling_groups[name]
        group.update(
            availability_zones=availability_zones,
            desired_capacity=desired_capacity,
            max_size=max_size,
            min_size=min_size,
            launch_config_name=launch_config_name,
            launch_template=launch_template,
            vpc_zone_identifier=vpc_zone_identifier,
            health_check_period=health_check_period,
            health_check_type=health_check_type,
            new_instances_protected_from_scale_in=new_instances_protected_from_scale_in,
            mixed_instances_policy=mixed_instances_policy,
        )
        return group

    def describe_auto_scaling_groups(
        self, names: list[str], filters: Optional[list[dict[str, str]]] = None
    ) -> list[FakeAutoScalingGroup]:
        groups = list(self.autoscaling_groups.values())

        if filters:
            for f in filters:
                if f["Name"] == "tag-key":
                    groups = [
                        group
                        for group in groups
                        if any(tag["Key"] in f["Values"] for tag in group.tags)
                    ]
                elif f["Name"] == "tag-value":
                    groups = [
                        group
                        for group in groups
                        if any(tag["Value"] in f["Values"] for tag in group.tags)
                    ]
                elif f["Name"].startswith("tag:"):
                    tag_key = f["Name"][4:]
                    groups = [
                        group
                        for group in groups
                        if any(
                            tag["Key"] == tag_key and tag["Value"] in f["Values"]
                            for tag in group.tags
                        )
                    ]

        if names:
            groups = [group for group in groups if group.name in names]

        return groups

    def delete_auto_scaling_group(self, group_name: str) -> None:
        self.set_desired_capacity(group_name, 0)
        self.autoscaling_groups.pop(group_name, None)

    def describe_auto_scaling_instances(
        self, instance_ids: list[str]
    ) -> list[InstanceState]:
        instance_states = []
        for group in self.autoscaling_groups.values():
            instance_states.extend(
                [
                    x
                    for x in group.instance_states
                    if not instance_ids or x.instance.id in instance_ids
                ]
            )
        return instance_states

    def attach_instances(self, group_name: str, instance_ids: list[str]) -> None:
        group = self.autoscaling_groups[group_name]
        original_size = len(group.instance_states)

        if (original_size + len(instance_ids)) > group.max_size:  # type: ignore[operator]
            raise ResourceContentionError
        else:
            group.desired_capacity = original_size + len(instance_ids)
            new_instances = [
                InstanceState(
                    self.ec2_backend.get_instance(x),
                    protected_from_scale_in=group.new_instances_protected_from_scale_in,
                    autoscaling_group=group,
                )
                for x in instance_ids
            ]
            for instance in new_instances:
                self.ec2_backend.create_tags(
                    [instance.instance.id], {ASG_NAME_TAG: group.name}
                )
            group.instance_states.extend(new_instances)
            self.update_attached_elbs(group.name)
            self.update_attached_target_groups(group.name)

    def set_instance_health(self, instance_id: str, health_status: str) -> None:
        """
        The ShouldRespectGracePeriod-parameter is not yet implemented
        """
        instance = self.ec2_backend.get_instance(instance_id)
        instance_state = next(
            instance_state
            for group in self.autoscaling_groups.values()
            for instance_state in group.instance_states
            if instance_state.instance.id == instance.id
        )
        instance_state.health_status = health_status

    def detach_instances(
        self, group_name: str, instance_ids: list[str], should_decrement: bool
    ) -> list[DetachInstanceActivity]:
        group = self.autoscaling_groups[group_name]
        original_size = group.desired_capacity
        activities = []
        detached_instance_states = [
            x for x in group.instance_states if x.instance.id in instance_ids
        ]
        for instance_state in detached_instance_states:
            instance = instance_state.instance
            self.ec2_backend.delete_tags([instance.id], {ASG_NAME_TAG: group.name})
            activity = DetachInstanceActivity(instance)
            activities.append(activity)

        new_instance_state = [
            x for x in group.instance_states if x.instance.id not in instance_ids
        ]
        group.instance_states = new_instance_state

        if should_decrement:
            group.desired_capacity = original_size - len(instance_ids)  # type: ignore[operator]

        group.set_desired_capacity(group.desired_capacity)
        return activities

    def set_desired_capacity(
        self, group_name: str, desired_capacity: Optional[int]
    ) -> None:
        group = self.autoscaling_groups[group_name]
        group.set_desired_capacity(desired_capacity)
        self.update_attached_elbs(group_name)

    def change_capacity(
        self, group_name: str, scaling_adjustment: Optional[int]
    ) -> None:
        group = self.autoscaling_groups[group_name]
        desired_capacity = group.desired_capacity + scaling_adjustment  # type: ignore[operator]
        self.set_desired_capacity(group_name, desired_capacity)

    def change_capacity_percent(
        self, group_name: str, scaling_adjustment: Optional[int]
    ) -> None:
        """http://docs.aws.amazon.com/AutoScaling/latest/DeveloperGuide/as-scale-based-on-demand.html
        If PercentChangeInCapacity returns a value between 0 and 1,
        Auto Scaling will round it off to 1. If the PercentChangeInCapacity
        returns a value greater than 1, Auto Scaling will round it off to the
        lower value. For example, if PercentChangeInCapacity returns 12.5,
        then Auto Scaling will round it off to 12."""
        group = self.autoscaling_groups[group_name]
        percent_change = 1 + (scaling_adjustment / 100.0)  # type: ignore[operator]
        desired_capacity = group.desired_capacity * percent_change  # type: ignore[operator]
        if group.desired_capacity < desired_capacity < group.desired_capacity + 1:  # type: ignore[operator]
            desired_capacity = group.desired_capacity + 1  # type: ignore[operator]
        else:
            desired_capacity = int(desired_capacity)
        self.set_desired_capacity(group_name, desired_capacity)

    def create_lifecycle_hook(
        self,
        name: str,
        as_name: str,
        transition: str,
        timeout: Optional[int],
        result: str,
    ) -> LifecycleHook:
        lifecycle_hook = LifecycleHook(name, as_name, transition, timeout, result)

        self.lifecycle_hooks[f"{as_name}_{name}"] = lifecycle_hook
        return lifecycle_hook

    def describe_lifecycle_hooks(
        self, as_name: str, lifecycle_hook_names: Optional[list[str]] = None
    ) -> list[LifecycleHook]:
        return [
            lifecycle_hook
            for lifecycle_hook in self.lifecycle_hooks.values()
            if (lifecycle_hook.auto_scaling_group_name == as_name)
            and (
                not lifecycle_hook_names or lifecycle_hook.name in lifecycle_hook_names
            )
        ]

    def delete_lifecycle_hook(self, as_name: str, name: str) -> None:
        self.lifecycle_hooks.pop(f"{as_name}_{name}", None)

    def put_scaling_policy(
        self,
        name: str,
        policy_type: str,
        metric_aggregation_type: str,
        adjustment_type: str,
        as_name: str,
        min_adjustment_magnitude: str,
        scaling_adjustment: Optional[int],
        cooldown: Optional[int],
        target_tracking_config: dict[str, Any],
        step_adjustments: str,
        estimated_instance_warmup: str,
        predictive_scaling_configuration: str,
    ) -> FakeScalingPolicy:
        policy = FakeScalingPolicy(
            name,
            policy_type,
            metric_aggregation_type,
            adjustment_type=adjustment_type,
            as_name=as_name,
            min_adjustment_magnitude=min_adjustment_magnitude,
            scaling_adjustment=scaling_adjustment,
            cooldown=cooldown,
            target_tracking_config=target_tracking_config,
            step_adjustments=step_adjustments,
            estimated_instance_warmup=estimated_instance_warmup,
            predictive_scaling_configuration=predictive_scaling_configuration,
            autoscaling_backend=self,
        )

        self.policies[name] = policy
        return policy

    def describe_policies(
        self,
        autoscaling_group_name: Optional[str] = None,
        policy_names: Optional[list[str]] = None,
        policy_types: Optional[list[str]] = None,
    ) -> list[FakeScalingPolicy]:
        return [
            policy
            for policy in self.policies.values()
            if (
                not autoscaling_group_name
                or policy.auto_scaling_group_name == autoscaling_group_name
            )
            and (not policy_names or policy.name in policy_names)
            and (not policy_types or policy.policy_type in policy_types)
        ]

    def delete_policy(self, group_name: str) -> None:
        self.policies.pop(group_name, None)

    def execute_policy(self, group_name: str) -> None:
        policy = self.policies[group_name]
        policy.execute()

    def update_attached_elbs(self, group_name: str) -> None:
        group = self.autoscaling_groups[group_name]
        group_instance_ids = {state.instance.id for state in group.active_instances()}

        # skip this if group.load_balancers is empty
        # otherwise elb_backend.describe_load_balancers returns all available load balancers
        if not group.load_balancer_names:
            return
        try:
            elbs = self.elb_backend.describe_load_balancers(
                names=group.load_balancer_names
            )
        except LoadBalancerNotFoundError:
            # ELBs can be deleted before their autoscaling group
            return

        for elb in elbs:
            elb_instace_ids = set(elb.instance_ids)
            self.elb_backend.register_instances(
                elb.name, group_instance_ids - elb_instace_ids, from_autoscaling=True
            )
            self.elb_backend.deregister_instances(
                elb.name, elb_instace_ids - group_instance_ids, from_autoscaling=True
            )

    def update_attached_target_groups(self, group_name: str) -> None:
        group = self.autoscaling_groups[group_name]
        group_instance_ids = {state.instance.id for state in group.instance_states}

        # no action necessary if target_group_arns is empty
        if not group.target_group_arns:
            return

        target_groups = self.elbv2_backend.describe_target_groups(
            target_group_arns=group.target_group_arns,
            load_balancer_arn=None,
            names=None,
        )

        for target_group in target_groups:
            asg_targets = [
                {"Id": x, "Port": target_group.port} for x in group_instance_ids
            ]
            self.elbv2_backend.register_targets(target_group.arn, (asg_targets))

    def create_or_update_tags(self, tags: list[dict[str, str]]) -> None:
        for tag in tags:
            group_name = tag["ResourceId"]
            group = self.autoscaling_groups[group_name]
            old_tags = group.tags

            new_tags = []
            # if key was in old_tags, update old tag
            for old_tag in old_tags:
                if old_tag["Key"] == tag["Key"]:
                    new_tags.append(tag)
                else:
                    new_tags.append(old_tag)

            # if key was never in old_tag's add it (create tag)
            if not any(new_tag["Key"] == tag["Key"] for new_tag in new_tags):
                new_tags.append(tag)

            group.tags = new_tags

    def delete_tags(self, tags: list[dict[str, str]]) -> None:
        for tag_to_delete in tags:
            group_name = tag_to_delete["ResourceId"]
            key_to_delete = tag_to_delete["Key"]
            group = self.autoscaling_groups[group_name]
            old_tags = group.tags
            group.tags = [x for x in old_tags if x["Key"] != key_to_delete]

    def attach_load_balancers(
        self, group_name: str, load_balancer_names: list[str]
    ) -> None:
        group = self.autoscaling_groups[group_name]
        group.load_balancer_names.extend(
            [x for x in load_balancer_names if x not in group.load_balancer_names]
        )
        self.update_attached_elbs(group_name)

    def describe_load_balancers(self, group_name: str) -> list[str]:
        return self.autoscaling_groups[group_name].load_balancer_names

    def detach_load_balancers(
        self, group_name: str, load_balancer_names: list[str]
    ) -> None:
        group = self.autoscaling_groups[group_name]
        group_instance_ids = {state.instance.id for state in group.instance_states}
        elbs = self.elb_backend.describe_load_balancers(names=group.load_balancer_names)
        for elb in elbs:
            self.elb_backend.deregister_instances(
                elb.name, group_instance_ids, from_autoscaling=True
            )
        group.load_balancer_names = [
            x for x in group.load_balancer_names if x not in load_balancer_names
        ]

    def attach_load_balancer_target_groups(
        self, group_name: str, target_group_arns: list[str]
    ) -> None:
        group = self.autoscaling_groups[group_name]
        group.append_target_groups(target_group_arns)
        self.update_attached_target_groups(group_name)

    def describe_load_balancer_target_groups(self, group_name: str) -> list[str]:
        return self.autoscaling_groups[group_name].target_group_arns

    def detach_load_balancer_target_groups(
        self, group_name: str, target_group_arns: list[str]
    ) -> None:
        group = self.autoscaling_groups[group_name]
        group.target_group_arns = [
            x for x in group.target_group_arns if x not in target_group_arns
        ]
        for target_group in target_group_arns:
            asg_targets = [{"Id": x.instance.id} for x in group.instance_states]
            self.elbv2_backend.deregister_targets(target_group, (asg_targets))

    def suspend_processes(self, group_name: str, scaling_processes: list[str]) -> None:
        all_proc_names = [
            "Launch",
            "Terminate",
            "AddToLoadBalancer",
            "AlarmNotification",
            "AZRebalance",
            "HealthCheck",
            "InstanceRefresh",
            "ReplaceUnhealthy",
            "ScheduledActions",
        ]
        group = self.autoscaling_groups[group_name]
        set_to_add = set(scaling_processes or all_proc_names)
        suspended_processes = [p["ProcessName"] for p in group.suspended_processes]
        group.suspended_processes = list(set(suspended_processes).union(set_to_add))

    def resume_processes(self, group_name: str, scaling_processes: list[str]) -> None:
        group = self.autoscaling_groups[group_name]
        if scaling_processes:
            suspended_processes = [p["ProcessName"] for p in group.suspended_processes]
            group.suspended_processes = list(
                set(suspended_processes).difference(set(scaling_processes))
            )
        else:
            group.suspended_processes = []

    def set_instance_protection(
        self,
        group_name: str,
        instance_ids: list[str],
        protected_from_scale_in: Optional[bool],
    ) -> None:
        group = self.autoscaling_groups[group_name]
        protected_instances = [
            x for x in group.instance_states if x.instance.id in instance_ids
        ]
        for instance in protected_instances:
            instance.protected_from_scale_in = protected_from_scale_in

    def notify_terminate_instances(self, instance_ids: list[str]) -> None:
        for (
            autoscaling_group_name,
            autoscaling_group,
        ) in self.autoscaling_groups.items():
            original_active_instance_count = len(autoscaling_group.active_instances())
            autoscaling_group.instance_states = list(
                filter(
                    lambda i_state: i_state.instance.id not in instance_ids,
                    autoscaling_group.instance_states,
                )
            )
            difference = original_active_instance_count - len(
                autoscaling_group.active_instances()
            )
            if difference > 0:
                autoscaling_group.replace_autoscaling_group_instances(
                    difference, autoscaling_group.get_propagated_tags()
                )
                self.update_attached_elbs(autoscaling_group_name)

    def enter_standby_instances(
        self, group_name: str, instance_ids: list[str], should_decrement: bool
    ) -> list[EnterStandbyActivity]:
        group = self.autoscaling_groups[group_name]
        activities = []
        for instance_state in group.instance_states:
            if instance_state.instance.id in instance_ids:
                instance_state.lifecycle_state = "Standby"
                original_size = group.desired_capacity
                instance = instance_state.instance
                if should_decrement:
                    group.desired_capacity -= 1  # type: ignore[operator]
                activity = EnterStandbyActivity(instance, original_size)
                activities.append(activity)
        group.set_desired_capacity(group.desired_capacity)
        return activities

    def exit_standby_instances(
        self, group_name: str, instance_ids: list[str]
    ) -> list[ExitStandbyActivity]:
        group = self.autoscaling_groups[group_name]
        activities = []
        for instance_state in group.instance_states:
            if instance_state.instance.id in instance_ids:
                instance_state.lifecycle_state = "InService"
                original_size = group.desired_capacity
                instance = instance_state.instance
                group.desired_capacity += 1  # type: ignore[operator]
                activity = ExitStandbyActivity(instance, original_size)
                activities.append(activity)
        group.set_desired_capacity(group.desired_capacity)
        return activities

    def terminate_instance(
        self, instance_id: str, should_decrement: bool
    ) -> TerminateInstanceActivity:
        instance_state = next(
            instance_state
            for group in self.autoscaling_groups.values()
            for instance_state in group.instance_states
            if instance_state.instance.id == instance_id
        )
        instance = instance_state.instance
        group = instance.autoscaling_group  # type: ignore[attr-defined]
        original_size = group.desired_capacity
        self.detach_instances(group.name, [instance.id], should_decrement)
        self.ec2_backend.terminate_instances([instance.id])
        return TerminateInstanceActivity(instance, original_size)

    def describe_tags(self, filters: list[dict[str, str]]) -> list[dict[str, str]]:
        """
        Pagination is not yet implemented.
        """
        resources = self.autoscaling_groups.values()
        tags = list(itertools.chain(*[r.tags for r in resources]))
        for f in filters:
            if f["Name"] == "auto-scaling-group":
                tags = [t for t in tags if t["ResourceId"] in f["Values"]]
            if f["Name"] == "propagate-at-launch":
                values = list(f["Values"])
                tags = [
                    t for t in tags if str(t.get("PropagateAtLaunch", False)) in values
                ]
            if f["Name"] == "key":
                tags = [t for t in tags if t["Key"] in f["Values"]]
            if f["Name"] == "value":
                tags = [t for t in tags if t["Value"] in f["Values"]]
        return tags

    def enable_metrics_collection(self, group_name: str, metrics: list[str]) -> None:
        group = self.describe_auto_scaling_groups([group_name])[0]
        group.enable_metrics_collection(metrics)

    def put_warm_pool(
        self,
        group_name: str,
        max_group_prepared_capacity: Optional[int],
        min_size: Optional[int],
        pool_state: Optional[str],
        instance_reuse_policy: Optional[dict[str, bool]],
    ) -> None:
        group = self.describe_auto_scaling_groups([group_name])[0]
        group.put_warm_pool(
            max_group_prepared_capacity=max_group_prepared_capacity,
            min_size=min_size,
            pool_state=pool_state,
            instance_reuse_policy=instance_reuse_policy,
        )

    def describe_warm_pool(self, group_name: str) -> Optional[FakeWarmPool]:
        """
        Pagination is not yet implemented. Does not create/return any Instances currently.
        """
        group = self.describe_auto_scaling_groups([group_name])[0]
        return group.get_warm_pool()

    def delete_warm_pool(self, group_name: str) -> None:
        group = self.describe_auto_scaling_groups([group_name])[0]
        group.warm_pool = None


autoscaling_backends = BackendDict(AutoScalingBackend, "autoscaling")
