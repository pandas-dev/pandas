from moto.core.common_types import TYPE_RESPONSE
from moto.core.responses import ActionResult, BaseResponse, EmptyResult
from moto.ec2.utils import parse_user_data
from moto.utilities.aws_headers import amz_crc32

from .models import AutoScalingBackend, autoscaling_backends


class AutoScalingResponse(BaseResponse):
    def __init__(self) -> None:
        super().__init__(service_name="autoscaling")
        self.automated_parameter_parsing = True

    @property
    def autoscaling_backend(self) -> AutoScalingBackend:
        return autoscaling_backends[self.current_account][self.region]

    @amz_crc32
    def call_action(self) -> TYPE_RESPONSE:
        return super().call_action()

    def create_launch_configuration(self) -> ActionResult:
        params = self._get_params()
        user_data = parse_user_data(params.get("UserData"))
        self.autoscaling_backend.create_launch_configuration(
            name=params.get("LaunchConfigurationName"),  # type: ignore[arg-type]
            image_id=params.get("ImageId"),  # type: ignore[arg-type]
            key_name=params.get("KeyName"),
            ramdisk_id=params.get("RamdiskId"),  # type: ignore[arg-type]
            kernel_id=params.get("KernelId"),  # type: ignore[arg-type]
            security_groups=self._get_param("SecurityGroups", []),
            user_data=user_data,
            instance_type=params.get("InstanceType"),  # type: ignore[arg-type]
            instance_monitoring=self._get_param("InstanceMonitoring.Enabled", False),
            instance_profile_name=params.get("IamInstanceProfile"),
            spot_price=params.get("SpotPrice"),
            ebs_optimized=self._get_bool_param("EbsOptimized", False),
            associate_public_ip_address=self._get_bool_param(
                "AssociatePublicIpAddress", False
            ),
            block_device_mappings=params.get("BlockDeviceMappings"),  # type: ignore[arg-type]
            instance_id=params.get("InstanceId"),
            metadata_options=params.get("MetadataOptions"),
            classic_link_vpc_id=params.get("ClassicLinkVPCId"),
            classic_link_vpc_security_groups=params.get("ClassicLinkVPCSecurityGroups"),
        )
        return EmptyResult()

    def describe_launch_configurations(self) -> ActionResult:
        names = self._get_param("LaunchConfigurationNames", [])
        all_launch_configurations = (
            self.autoscaling_backend.describe_launch_configurations(names)
        )
        marker = self._get_param("NextToken")
        all_names = [lc.name for lc in all_launch_configurations]
        if marker:
            start = all_names.index(marker) + 1
        else:
            start = 0
        # the default is 100, but using 50 to make testing easier
        max_records = self._get_int_param("MaxRecords") or 50
        launch_configurations_resp = all_launch_configurations[
            start : start + max_records
        ]
        next_token = None
        if len(all_launch_configurations) > start + max_records:
            next_token = launch_configurations_resp[-1].name

        result = {
            "LaunchConfigurations": launch_configurations_resp,
            "NextToken": next_token,
        }
        return ActionResult(result)

    def delete_launch_configuration(self) -> ActionResult:
        launch_configurations_name = self._get_param("LaunchConfigurationName")
        self.autoscaling_backend.delete_launch_configuration(launch_configurations_name)
        return EmptyResult()

    def create_auto_scaling_group(self) -> ActionResult:
        params = self._get_params()
        self.autoscaling_backend.create_auto_scaling_group(
            name=self._get_param("AutoScalingGroupName"),
            availability_zones=self._get_param("AvailabilityZones", []),
            desired_capacity=self._get_int_param("DesiredCapacity"),
            max_size=self._get_int_param("MaxSize"),
            min_size=self._get_int_param("MinSize"),
            instance_id=self._get_param("InstanceId"),
            launch_config_name=self._get_param("LaunchConfigurationName"),
            launch_template=self._get_param("LaunchTemplate", {}),
            mixed_instances_policy=params.get("MixedInstancesPolicy"),
            vpc_zone_identifier=self._get_param("VPCZoneIdentifier"),
            default_cooldown=self._get_int_param("DefaultCooldown"),
            health_check_period=self._get_int_param("HealthCheckGracePeriod"),
            health_check_type=self._get_param("HealthCheckType"),
            load_balancers=self._get_param("LoadBalancerNames", []),
            target_group_arns=self._get_param("TargetGroupARNs", []),
            placement_group=self._get_param("PlacementGroup"),
            termination_policies=self._get_param("TerminationPolicies", []),
            tags=params.get("Tags", []),
            capacity_rebalance=self._get_bool_param("CapacityRebalance", False),
            new_instances_protected_from_scale_in=self._get_bool_param(
                "NewInstancesProtectedFromScaleIn", False
            ),
        )
        return EmptyResult()

    def put_scheduled_update_group_action(self) -> ActionResult:
        self.autoscaling_backend.put_scheduled_update_group_action(
            name=self._get_param("AutoScalingGroupName"),
            desired_capacity=self._get_int_param("DesiredCapacity"),
            max_size=self._get_int_param("MaxSize"),
            min_size=self._get_int_param("MinSize"),
            scheduled_action_name=self._get_param("ScheduledActionName"),
            start_time=self._get_param("StartTime"),
            end_time=self._get_param("EndTime"),
            recurrence=self._get_param("Recurrence"),
            timezone=self._get_param("TimeZone"),
        )
        return EmptyResult()

    def batch_put_scheduled_update_group_action(self) -> ActionResult:
        failed_actions = (
            self.autoscaling_backend.batch_put_scheduled_update_group_action(
                name=self._get_param("AutoScalingGroupName"),
                actions=self._get_param("ScheduledUpdateGroupActions", []),
            )
        )
        result = {"FailedScheduledUpdateGroupActions": failed_actions}
        return ActionResult(result)

    def describe_scheduled_actions(self) -> ActionResult:
        scheduled_actions = self.autoscaling_backend.describe_scheduled_actions(
            autoscaling_group_name=self._get_param("AutoScalingGroupName"),
            scheduled_action_names=self._get_param("ScheduledActionNames", []),
        )
        result = {"ScheduledUpdateGroupActions": scheduled_actions}
        return ActionResult(result)

    def delete_scheduled_action(self) -> ActionResult:
        auto_scaling_group_name = self._get_param("AutoScalingGroupName")
        scheduled_action_name = self._get_param("ScheduledActionName")
        self.autoscaling_backend.delete_scheduled_action(
            auto_scaling_group_name=auto_scaling_group_name,
            scheduled_action_name=scheduled_action_name,
        )
        return EmptyResult()

    def batch_delete_scheduled_action(self) -> ActionResult:
        auto_scaling_group_name = self._get_param("AutoScalingGroupName")
        scheduled_action_names = self._get_param("ScheduledActionNames", [])
        failed_actions = self.autoscaling_backend.batch_delete_scheduled_action(
            auto_scaling_group_name=auto_scaling_group_name,
            scheduled_action_names=scheduled_action_names,
        )
        result = {"FailedScheduledActions": failed_actions}
        return ActionResult(result)

    def describe_scaling_activities(self) -> ActionResult:
        result = {"Activities": []}  # type: ignore[var-annotated]
        return ActionResult(result)

    def attach_instances(self) -> ActionResult:
        group_name = self._get_param("AutoScalingGroupName")
        instance_ids = self._get_param("InstanceIds", [])
        self.autoscaling_backend.attach_instances(group_name, instance_ids)
        return EmptyResult()

    def set_instance_health(self) -> ActionResult:
        instance_id = self._get_param("InstanceId")
        health_status = self._get_param("HealthStatus")
        if health_status not in ["Healthy", "Unhealthy"]:
            raise ValueError("Valid instance health states are: [Healthy, Unhealthy]")
        self.autoscaling_backend.set_instance_health(instance_id, health_status)
        return EmptyResult()

    def detach_instances(self) -> ActionResult:
        group_name = self._get_param("AutoScalingGroupName")
        instance_ids = self._get_param("InstanceIds", [])
        should_decrement = self._get_bool_param("ShouldDecrementDesiredCapacity", False)
        activities = self.autoscaling_backend.detach_instances(
            group_name, instance_ids, should_decrement
        )
        result = {"Activities": activities}
        return ActionResult(result)

    def attach_load_balancer_target_groups(self) -> ActionResult:
        group_name = self._get_param("AutoScalingGroupName")
        target_group_arns = self._get_param("TargetGroupARNs", [])

        self.autoscaling_backend.attach_load_balancer_target_groups(
            group_name, target_group_arns
        )
        return EmptyResult()

    def describe_load_balancer_target_groups(self) -> ActionResult:
        group_name = self._get_param("AutoScalingGroupName")
        target_group_arns = (
            self.autoscaling_backend.describe_load_balancer_target_groups(group_name)
        )
        result = {
            "LoadBalancerTargetGroups": [
                {"LoadBalancerTargetGroupARN": arn, "State": "Added"}
                for arn in target_group_arns
            ]
        }
        return ActionResult(result)

    def detach_load_balancer_target_groups(self) -> ActionResult:
        group_name = self._get_param("AutoScalingGroupName")
        target_group_arns = self._get_param("TargetGroupARNs", [])

        self.autoscaling_backend.detach_load_balancer_target_groups(
            group_name, target_group_arns
        )
        return EmptyResult()

    def describe_auto_scaling_groups(self) -> ActionResult:
        names = self._get_param("AutoScalingGroupNames", [])
        token = self._get_param("NextToken")
        filters = self._get_param("Filters", [])
        all_groups = self.autoscaling_backend.describe_auto_scaling_groups(
            names, filters=filters
        )
        all_names = [group.name for group in all_groups]
        if token:
            start = all_names.index(token) + 1
        else:
            start = 0
        max_records = self._get_int_param("MaxRecords", 50)
        if max_records > 100:
            raise ValueError
        groups = all_groups[start : start + max_records]
        next_token = None
        if max_records and len(all_groups) > start + max_records:
            next_token = groups[-1].name
        result = {"AutoScalingGroups": groups, "NextToken": next_token}
        return ActionResult(result)

    def update_auto_scaling_group(self) -> ActionResult:
        self.autoscaling_backend.update_auto_scaling_group(
            name=self._get_param("AutoScalingGroupName"),
            availability_zones=self._get_param("AvailabilityZones", []),
            desired_capacity=self._get_int_param("DesiredCapacity"),
            max_size=self._get_int_param("MaxSize"),
            min_size=self._get_int_param("MinSize"),
            launch_config_name=self._get_param("LaunchConfigurationName"),
            launch_template=self._get_param("LaunchTemplate", {}),
            vpc_zone_identifier=self._get_param("VPCZoneIdentifier"),
            health_check_period=self._get_int_param("HealthCheckGracePeriod"),
            health_check_type=self._get_param("HealthCheckType"),
            new_instances_protected_from_scale_in=self._get_bool_param(
                "NewInstancesProtectedFromScaleIn", None
            ),
            mixed_instances_policy=self._get_param("MixedInstancesPolicy"),
        )
        return EmptyResult()

    def delete_auto_scaling_group(self) -> ActionResult:
        group_name = self._get_param("AutoScalingGroupName")
        self.autoscaling_backend.delete_auto_scaling_group(group_name)
        return EmptyResult()

    def set_desired_capacity(self) -> ActionResult:
        group_name = self._get_param("AutoScalingGroupName")
        desired_capacity = self._get_int_param("DesiredCapacity")
        self.autoscaling_backend.set_desired_capacity(group_name, desired_capacity)
        return EmptyResult()

    def create_or_update_tags(self) -> ActionResult:
        self.autoscaling_backend.create_or_update_tags(self._get_param("Tags", []))
        return EmptyResult()

    def delete_tags(self) -> ActionResult:
        self.autoscaling_backend.delete_tags(self._get_params().get("Tags", []))
        return EmptyResult()

    def describe_auto_scaling_instances(self) -> ActionResult:
        instance_states = self.autoscaling_backend.describe_auto_scaling_instances(
            instance_ids=self._get_param("InstanceIds", [])
        )
        result = {"AutoScalingInstances": instance_states}
        return ActionResult(result)

    def put_lifecycle_hook(self) -> ActionResult:
        self.autoscaling_backend.create_lifecycle_hook(
            name=self._get_param("LifecycleHookName"),
            as_name=self._get_param("AutoScalingGroupName"),
            transition=self._get_param("LifecycleTransition"),
            timeout=self._get_int_param("HeartbeatTimeout"),
            result=self._get_param("DefaultResult"),
        )
        return EmptyResult()

    def describe_lifecycle_hooks(self) -> ActionResult:
        lifecycle_hooks = self.autoscaling_backend.describe_lifecycle_hooks(
            as_name=self._get_param("AutoScalingGroupName"),
            lifecycle_hook_names=self._get_param("LifecycleHookNames", []),
        )
        result = {"LifecycleHooks": lifecycle_hooks}
        return ActionResult(result)

    def delete_lifecycle_hook(self) -> ActionResult:
        as_name = self._get_param("AutoScalingGroupName")
        name = self._get_param("LifecycleHookName")
        self.autoscaling_backend.delete_lifecycle_hook(as_name, name)
        return EmptyResult()

    def put_scaling_policy(self) -> ActionResult:
        params = self._get_params()
        policy = self.autoscaling_backend.put_scaling_policy(
            name=params.get("PolicyName"),  # type: ignore[arg-type]
            policy_type=params.get("PolicyType", "SimpleScaling"),
            metric_aggregation_type=params.get("MetricAggregationType"),  # type: ignore[arg-type]
            adjustment_type=params.get("AdjustmentType"),  # type: ignore[arg-type]
            as_name=params.get("AutoScalingGroupName"),  # type: ignore[arg-type]
            min_adjustment_magnitude=params.get("MinAdjustmentMagnitude"),  # type: ignore[arg-type]
            scaling_adjustment=self._get_int_param("ScalingAdjustment"),
            cooldown=self._get_int_param("Cooldown"),
            target_tracking_config=params.get("TargetTrackingConfiguration", {}),
            step_adjustments=params.get("StepAdjustments", []),
            estimated_instance_warmup=params.get("EstimatedInstanceWarmup"),  # type: ignore[arg-type]
            predictive_scaling_configuration=params.get(
                "PredictiveScalingConfiguration", {}
            ),
        )
        return ActionResult({"PolicyArn": policy.arn})

    def describe_policies(self) -> ActionResult:
        policies = self.autoscaling_backend.describe_policies(
            autoscaling_group_name=self._get_param("AutoScalingGroupName"),
            policy_names=self._get_param("PolicyNames", []),
            policy_types=self._get_param("PolicyTypes", []),
        )
        result = {"ScalingPolicies": policies}
        return ActionResult(result)

    def delete_policy(self) -> ActionResult:
        group_name = self._get_param("PolicyName")
        self.autoscaling_backend.delete_policy(group_name)
        return EmptyResult()

    def execute_policy(self) -> ActionResult:
        group_name = self._get_param("PolicyName")
        self.autoscaling_backend.execute_policy(group_name)
        return EmptyResult()

    def attach_load_balancers(self) -> ActionResult:
        group_name = self._get_param("AutoScalingGroupName")
        load_balancer_names = self._get_param("LoadBalancerNames", [])
        self.autoscaling_backend.attach_load_balancers(group_name, load_balancer_names)
        return EmptyResult()

    def describe_load_balancers(self) -> ActionResult:
        group_name = self._get_param("AutoScalingGroupName")
        load_balancers = self.autoscaling_backend.describe_load_balancers(group_name)
        result = {
            "LoadBalancers": [
                {"LoadBalancerName": name, "State": "Added"} for name in load_balancers
            ]
        }
        return ActionResult(result)

    def detach_load_balancers(self) -> ActionResult:
        group_name = self._get_param("AutoScalingGroupName")
        load_balancer_names = self._get_param("LoadBalancerNames", [])
        self.autoscaling_backend.detach_load_balancers(group_name, load_balancer_names)
        return EmptyResult()

    def enter_standby(self) -> ActionResult:
        group_name = self._get_param("AutoScalingGroupName")
        instance_ids = self._get_param("InstanceIds", [])
        should_decrement = self._get_bool_param("ShouldDecrementDesiredCapacity")
        activities = self.autoscaling_backend.enter_standby_instances(
            group_name, instance_ids, should_decrement
        )
        result = {"Activities": activities}
        return ActionResult(result)

    def exit_standby(self) -> ActionResult:
        group_name = self._get_param("AutoScalingGroupName")
        instance_ids = self._get_param("InstanceIds", [])
        activities = self.autoscaling_backend.exit_standby_instances(
            group_name, instance_ids
        )
        result = {"Activities": activities}
        return ActionResult(result)

    def suspend_processes(self) -> ActionResult:
        autoscaling_group_name = self._get_param("AutoScalingGroupName")
        scaling_processes = self._get_param("ScalingProcesses", [])
        self.autoscaling_backend.suspend_processes(
            autoscaling_group_name, scaling_processes
        )
        return EmptyResult()

    def resume_processes(self) -> ActionResult:
        autoscaling_group_name = self._get_param("AutoScalingGroupName")
        scaling_processes = self._get_param("ScalingProcesses", [])
        self.autoscaling_backend.resume_processes(
            autoscaling_group_name, scaling_processes
        )
        return EmptyResult()

    def set_instance_protection(self) -> ActionResult:
        group_name = self._get_param("AutoScalingGroupName")
        instance_ids = self._get_param("InstanceIds", [])
        protected_from_scale_in = self._get_bool_param("ProtectedFromScaleIn")
        self.autoscaling_backend.set_instance_protection(
            group_name, instance_ids, protected_from_scale_in
        )
        return EmptyResult()

    def terminate_instance_in_auto_scaling_group(self) -> ActionResult:
        instance_id = self._get_param("InstanceId")
        should_decrement = self._get_bool_param("ShouldDecrementDesiredCapacity", False)
        activity = self.autoscaling_backend.terminate_instance(
            instance_id, should_decrement
        )
        result = {"Activity": activity}
        return ActionResult(result)

    def describe_tags(self) -> ActionResult:
        filters = self._get_param("Filters", [])
        tags = self.autoscaling_backend.describe_tags(filters=filters)
        result = {"Tags": tags, "NextToken": None}
        return ActionResult(result)

    def enable_metrics_collection(self) -> ActionResult:
        group_name = self._get_param("AutoScalingGroupName")
        metrics = self._get_param("Metrics")
        self.autoscaling_backend.enable_metrics_collection(group_name, metrics)  # type: ignore[arg-type]
        return EmptyResult()

    def put_warm_pool(self) -> ActionResult:
        params = self._get_params()
        group_name = params.get("AutoScalingGroupName")
        max_group_prepared_capacity = params.get("MaxGroupPreparedCapacity")
        min_size = params.get("MinSize")
        pool_state = params.get("PoolState")
        instance_reuse_policy = params.get("InstanceReusePolicy")
        self.autoscaling_backend.put_warm_pool(
            group_name=group_name,  # type: ignore[arg-type]
            max_group_prepared_capacity=max_group_prepared_capacity,
            min_size=min_size,
            pool_state=pool_state,
            instance_reuse_policy=instance_reuse_policy,
        )
        return EmptyResult()

    def describe_warm_pool(self) -> ActionResult:
        group_name = self._get_param("AutoScalingGroupName")
        warm_pool = self.autoscaling_backend.describe_warm_pool(group_name=group_name)
        result = {"WarmPoolConfiguration": warm_pool, "Instances": []}  # type: ignore[var-annotated]
        return ActionResult(result)

    def delete_warm_pool(self) -> ActionResult:
        group_name = self._get_param("AutoScalingGroupName")
        self.autoscaling_backend.delete_warm_pool(group_name=group_name)
        return EmptyResult()
