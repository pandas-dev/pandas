from moto.core.common_types import TYPE_RESPONSE
from moto.core.responses import BaseResponse
from moto.core.utils import iso_8601_datetime_with_milliseconds
from moto.utilities.aws_headers import amz_crc32

from .models import AutoScalingBackend, autoscaling_backends


class AutoScalingResponse(BaseResponse):
    def __init__(self) -> None:
        super().__init__(service_name="autoscaling")

    @property
    def autoscaling_backend(self) -> AutoScalingBackend:
        return autoscaling_backends[self.current_account][self.region]

    @amz_crc32
    def call_action(self) -> TYPE_RESPONSE:
        return super().call_action()

    def create_launch_configuration(self) -> str:
        instance_monitoring_string = self._get_param("InstanceMonitoring.Enabled")
        if instance_monitoring_string == "true":
            instance_monitoring = True
        else:
            instance_monitoring = False
        params = self._get_params()
        self.autoscaling_backend.create_launch_configuration(
            name=params.get("LaunchConfigurationName"),  # type: ignore[arg-type]
            image_id=params.get("ImageId"),  # type: ignore[arg-type]
            key_name=params.get("KeyName"),
            ramdisk_id=params.get("RamdiskId"),  # type: ignore[arg-type]
            kernel_id=params.get("KernelId"),  # type: ignore[arg-type]
            security_groups=self._get_multi_param("SecurityGroups.member"),
            user_data=params.get("UserData"),  # type: ignore[arg-type]
            instance_type=params.get("InstanceType"),  # type: ignore[arg-type]
            instance_monitoring=instance_monitoring,
            instance_profile_name=params.get("IamInstanceProfile"),
            spot_price=params.get("SpotPrice"),
            ebs_optimized=params.get("EbsOptimized"),  # type: ignore[arg-type]
            associate_public_ip_address=params.get("AssociatePublicIpAddress"),  # type: ignore[arg-type]
            block_device_mappings=params.get("BlockDeviceMappings"),  # type: ignore[arg-type]
            instance_id=params.get("InstanceId"),
            metadata_options=params.get("MetadataOptions"),
            classic_link_vpc_id=params.get("ClassicLinkVPCId"),
            classic_link_vpc_security_groups=params.get("ClassicLinkVPCSecurityGroups"),
        )
        template = self.response_template(CREATE_LAUNCH_CONFIGURATION_TEMPLATE)
        return template.render()

    def describe_launch_configurations(self) -> str:
        names = self._get_multi_param("LaunchConfigurationNames.member")
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

        template = self.response_template(DESCRIBE_LAUNCH_CONFIGURATIONS_TEMPLATE)
        return template.render(
            launch_configurations=launch_configurations_resp, next_token=next_token
        )

    def delete_launch_configuration(self) -> str:
        launch_configurations_name = self.querystring.get("LaunchConfigurationName")[0]  # type: ignore[index]
        self.autoscaling_backend.delete_launch_configuration(launch_configurations_name)
        template = self.response_template(DELETE_LAUNCH_CONFIGURATION_TEMPLATE)
        return template.render()

    def create_auto_scaling_group(self) -> str:
        params = self._get_params()
        self.autoscaling_backend.create_auto_scaling_group(
            name=self._get_param("AutoScalingGroupName"),
            availability_zones=self._get_multi_param("AvailabilityZones.member"),
            desired_capacity=self._get_int_param("DesiredCapacity"),
            max_size=self._get_int_param("MaxSize"),
            min_size=self._get_int_param("MinSize"),
            instance_id=self._get_param("InstanceId"),
            launch_config_name=self._get_param("LaunchConfigurationName"),
            launch_template=self._get_dict_param("LaunchTemplate."),
            mixed_instances_policy=params.get("MixedInstancesPolicy"),
            vpc_zone_identifier=self._get_param("VPCZoneIdentifier"),
            default_cooldown=self._get_int_param("DefaultCooldown"),
            health_check_period=self._get_int_param("HealthCheckGracePeriod"),
            health_check_type=self._get_param("HealthCheckType"),
            load_balancers=self._get_multi_param("LoadBalancerNames.member"),
            target_group_arns=self._get_multi_param("TargetGroupARNs.member"),
            placement_group=self._get_param("PlacementGroup"),
            termination_policies=self._get_multi_param("TerminationPolicies.member"),
            tags=params.get("Tags", []),
            capacity_rebalance=self._get_bool_param("CapacityRebalance", False),
            new_instances_protected_from_scale_in=self._get_bool_param(
                "NewInstancesProtectedFromScaleIn", False
            ),
        )
        template = self.response_template(CREATE_AUTOSCALING_GROUP_TEMPLATE)
        return template.render()

    def put_scheduled_update_group_action(self) -> str:
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
        template = self.response_template(PUT_SCHEDULED_UPDATE_GROUP_ACTION_TEMPLATE)
        return template.render()

    def batch_put_scheduled_update_group_action(self) -> str:
        failed_actions = (
            self.autoscaling_backend.batch_put_scheduled_update_group_action(
                name=self._get_param("AutoScalingGroupName"),
                actions=self._get_multi_param("ScheduledUpdateGroupActions.member"),
            )
        )
        template = self.response_template(
            BATCH_PUT_SCHEDULED_UPDATE_GROUP_ACTION_TEMPLATE
        )
        return template.render(failed_actions=failed_actions)

    def describe_scheduled_actions(self) -> str:
        scheduled_actions = self.autoscaling_backend.describe_scheduled_actions(
            autoscaling_group_name=self._get_param("AutoScalingGroupName"),
            scheduled_action_names=self._get_multi_param("ScheduledActionNames.member"),
        )
        template = self.response_template(DESCRIBE_SCHEDULED_ACTIONS)
        return template.render(scheduled_actions=scheduled_actions)

    def delete_scheduled_action(self) -> str:
        auto_scaling_group_name = self._get_param("AutoScalingGroupName")
        scheduled_action_name = self._get_param("ScheduledActionName")
        self.autoscaling_backend.delete_scheduled_action(
            auto_scaling_group_name=auto_scaling_group_name,
            scheduled_action_name=scheduled_action_name,
        )
        template = self.response_template(DELETE_SCHEDULED_ACTION_TEMPLATE)
        return template.render()

    def batch_delete_scheduled_action(self) -> str:
        auto_scaling_group_name = self._get_param("AutoScalingGroupName")
        scheduled_action_names = self._get_multi_param("ScheduledActionNames.member")
        failed_actions = self.autoscaling_backend.batch_delete_scheduled_action(
            auto_scaling_group_name=auto_scaling_group_name,
            scheduled_action_names=scheduled_action_names,
        )
        template = self.response_template(BATCH_DELETE_SCHEDULED_ACTION_TEMPLATE)
        return template.render(failed_actions=failed_actions)

    def describe_scaling_activities(self) -> str:
        template = self.response_template(DESCRIBE_SCALING_ACTIVITIES_TEMPLATE)
        return template.render()

    def attach_instances(self) -> str:
        group_name = self._get_param("AutoScalingGroupName")
        instance_ids = self._get_multi_param("InstanceIds.member")
        self.autoscaling_backend.attach_instances(group_name, instance_ids)
        template = self.response_template(ATTACH_INSTANCES_TEMPLATE)
        return template.render()

    def set_instance_health(self) -> str:
        instance_id = self._get_param("InstanceId")
        health_status = self._get_param("HealthStatus")
        if health_status not in ["Healthy", "Unhealthy"]:
            raise ValueError("Valid instance health states are: [Healthy, Unhealthy]")
        self.autoscaling_backend.set_instance_health(instance_id, health_status)
        template = self.response_template(SET_INSTANCE_HEALTH_TEMPLATE)
        return template.render()

    def detach_instances(self) -> str:
        group_name = self._get_param("AutoScalingGroupName")
        instance_ids = self._get_multi_param("InstanceIds.member")
        should_decrement_string = self._get_param("ShouldDecrementDesiredCapacity")
        if should_decrement_string == "true":
            should_decrement = True
        else:
            should_decrement = False
        detached_instances = self.autoscaling_backend.detach_instances(
            group_name, instance_ids, should_decrement
        )
        template = self.response_template(DETACH_INSTANCES_TEMPLATE)
        return template.render(detached_instances=detached_instances)

    def attach_load_balancer_target_groups(self) -> str:
        group_name = self._get_param("AutoScalingGroupName")
        target_group_arns = self._get_multi_param("TargetGroupARNs.member")

        self.autoscaling_backend.attach_load_balancer_target_groups(
            group_name, target_group_arns
        )
        template = self.response_template(ATTACH_LOAD_BALANCER_TARGET_GROUPS_TEMPLATE)
        return template.render()

    def describe_load_balancer_target_groups(self) -> str:
        group_name = self._get_param("AutoScalingGroupName")
        target_group_arns = (
            self.autoscaling_backend.describe_load_balancer_target_groups(group_name)
        )
        template = self.response_template(DESCRIBE_LOAD_BALANCER_TARGET_GROUPS)
        return template.render(target_group_arns=target_group_arns)

    def detach_load_balancer_target_groups(self) -> str:
        group_name = self._get_param("AutoScalingGroupName")
        target_group_arns = self._get_multi_param("TargetGroupARNs.member")

        self.autoscaling_backend.detach_load_balancer_target_groups(
            group_name, target_group_arns
        )
        template = self.response_template(DETACH_LOAD_BALANCER_TARGET_GROUPS_TEMPLATE)
        return template.render()

    def describe_auto_scaling_groups(self) -> str:
        names = self._get_multi_param("AutoScalingGroupNames.member")
        token = self._get_param("NextToken")
        filters = self._get_params().get("Filters", [])
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
        template = self.response_template(DESCRIBE_AUTOSCALING_GROUPS_TEMPLATE)
        return template.render(groups=groups, next_token=next_token)

    def update_auto_scaling_group(self) -> str:
        self.autoscaling_backend.update_auto_scaling_group(
            name=self._get_param("AutoScalingGroupName"),
            availability_zones=self._get_multi_param("AvailabilityZones.member"),
            desired_capacity=self._get_int_param("DesiredCapacity"),
            max_size=self._get_int_param("MaxSize"),
            min_size=self._get_int_param("MinSize"),
            launch_config_name=self._get_param("LaunchConfigurationName"),
            launch_template=self._get_dict_param("LaunchTemplate."),
            vpc_zone_identifier=self._get_param("VPCZoneIdentifier"),
            health_check_period=self._get_int_param("HealthCheckGracePeriod"),
            health_check_type=self._get_param("HealthCheckType"),
            new_instances_protected_from_scale_in=self._get_bool_param(
                "NewInstancesProtectedFromScaleIn", None
            ),
        )
        template = self.response_template(UPDATE_AUTOSCALING_GROUP_TEMPLATE)
        return template.render()

    def delete_auto_scaling_group(self) -> str:
        group_name = self._get_param("AutoScalingGroupName")
        self.autoscaling_backend.delete_auto_scaling_group(group_name)
        template = self.response_template(DELETE_AUTOSCALING_GROUP_TEMPLATE)
        return template.render()

    def set_desired_capacity(self) -> str:
        group_name = self._get_param("AutoScalingGroupName")
        desired_capacity = self._get_int_param("DesiredCapacity")
        self.autoscaling_backend.set_desired_capacity(group_name, desired_capacity)
        template = self.response_template(SET_DESIRED_CAPACITY_TEMPLATE)
        return template.render()

    def create_or_update_tags(self) -> str:
        self.autoscaling_backend.create_or_update_tags(
            self._get_params().get("Tags", [])
        )
        template = self.response_template(UPDATE_AUTOSCALING_GROUP_TEMPLATE)
        return template.render()

    def delete_tags(self) -> str:
        self.autoscaling_backend.delete_tags(self._get_params().get("Tags", []))
        template = self.response_template(UPDATE_AUTOSCALING_GROUP_TEMPLATE)
        return template.render()

    def describe_auto_scaling_instances(self) -> str:
        instance_states = self.autoscaling_backend.describe_auto_scaling_instances(
            instance_ids=self._get_multi_param("InstanceIds.member")
        )
        template = self.response_template(DESCRIBE_AUTOSCALING_INSTANCES_TEMPLATE)
        return template.render(instance_states=instance_states)

    def put_lifecycle_hook(self) -> str:
        lifecycle_hook = self.autoscaling_backend.create_lifecycle_hook(
            name=self._get_param("LifecycleHookName"),
            as_name=self._get_param("AutoScalingGroupName"),
            transition=self._get_param("LifecycleTransition"),
            timeout=self._get_int_param("HeartbeatTimeout"),
            result=self._get_param("DefaultResult"),
        )
        template = self.response_template(CREATE_LIFECYLE_HOOK_TEMPLATE)
        return template.render(lifecycle_hook=lifecycle_hook)

    def describe_lifecycle_hooks(self) -> str:
        lifecycle_hooks = self.autoscaling_backend.describe_lifecycle_hooks(
            as_name=self._get_param("AutoScalingGroupName"),
            lifecycle_hook_names=self._get_multi_param("LifecycleHookNames.member"),
        )
        template = self.response_template(DESCRIBE_LIFECYCLE_HOOKS_TEMPLATE)
        return template.render(lifecycle_hooks=lifecycle_hooks)

    def delete_lifecycle_hook(self) -> str:
        as_name = self._get_param("AutoScalingGroupName")
        name = self._get_param("LifecycleHookName")
        self.autoscaling_backend.delete_lifecycle_hook(as_name, name)
        template = self.response_template(DELETE_LIFECYCLE_HOOK_TEMPLATE)
        return template.render()

    def put_scaling_policy(self) -> str:
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
        template = self.response_template(CREATE_SCALING_POLICY_TEMPLATE)
        return template.render(policy=policy)

    def describe_policies(self) -> str:
        policies = self.autoscaling_backend.describe_policies(
            autoscaling_group_name=self._get_param("AutoScalingGroupName"),
            policy_names=self._get_multi_param("PolicyNames.member"),
            policy_types=self._get_multi_param("PolicyTypes.member"),
        )
        template = self.response_template(DESCRIBE_SCALING_POLICIES_TEMPLATE)
        return template.render(policies=policies)

    def delete_policy(self) -> str:
        group_name = self._get_param("PolicyName")
        self.autoscaling_backend.delete_policy(group_name)
        template = self.response_template(DELETE_POLICY_TEMPLATE)
        return template.render()

    def execute_policy(self) -> str:
        group_name = self._get_param("PolicyName")
        self.autoscaling_backend.execute_policy(group_name)
        template = self.response_template(EXECUTE_POLICY_TEMPLATE)
        return template.render()

    def attach_load_balancers(self) -> str:
        group_name = self._get_param("AutoScalingGroupName")
        load_balancer_names = self._get_multi_param("LoadBalancerNames.member")
        self.autoscaling_backend.attach_load_balancers(group_name, load_balancer_names)
        template = self.response_template(ATTACH_LOAD_BALANCERS_TEMPLATE)
        return template.render()

    def describe_load_balancers(self) -> str:
        group_name = self._get_param("AutoScalingGroupName")
        load_balancers = self.autoscaling_backend.describe_load_balancers(group_name)
        template = self.response_template(DESCRIBE_LOAD_BALANCERS_TEMPLATE)
        return template.render(load_balancers=load_balancers)

    def detach_load_balancers(self) -> str:
        group_name = self._get_param("AutoScalingGroupName")
        load_balancer_names = self._get_multi_param("LoadBalancerNames.member")
        self.autoscaling_backend.detach_load_balancers(group_name, load_balancer_names)
        template = self.response_template(DETACH_LOAD_BALANCERS_TEMPLATE)
        return template.render()

    def enter_standby(self) -> str:
        group_name = self._get_param("AutoScalingGroupName")
        instance_ids = self._get_multi_param("InstanceIds.member")
        should_decrement_string = self._get_param("ShouldDecrementDesiredCapacity")
        if should_decrement_string == "true":
            should_decrement = True
        else:
            should_decrement = False
        (
            standby_instances,
            original_size,
            desired_capacity,
        ) = self.autoscaling_backend.enter_standby_instances(
            group_name, instance_ids, should_decrement
        )
        template = self.response_template(ENTER_STANDBY_TEMPLATE)
        return template.render(
            standby_instances=standby_instances,
            should_decrement=should_decrement,
            original_size=original_size,
            desired_capacity=desired_capacity,
            timestamp=iso_8601_datetime_with_milliseconds(),
        )

    def exit_standby(self) -> str:
        group_name = self._get_param("AutoScalingGroupName")
        instance_ids = self._get_multi_param("InstanceIds.member")
        (
            standby_instances,
            original_size,
            desired_capacity,
        ) = self.autoscaling_backend.exit_standby_instances(group_name, instance_ids)
        template = self.response_template(EXIT_STANDBY_TEMPLATE)
        return template.render(
            standby_instances=standby_instances,
            original_size=original_size,
            desired_capacity=desired_capacity,
            timestamp=iso_8601_datetime_with_milliseconds(),
        )

    def suspend_processes(self) -> str:
        autoscaling_group_name = self._get_param("AutoScalingGroupName")
        scaling_processes = self._get_multi_param("ScalingProcesses.member")
        self.autoscaling_backend.suspend_processes(
            autoscaling_group_name, scaling_processes
        )
        template = self.response_template(SUSPEND_PROCESSES_TEMPLATE)
        return template.render()

    def resume_processes(self) -> str:
        autoscaling_group_name = self._get_param("AutoScalingGroupName")
        scaling_processes = self._get_multi_param("ScalingProcesses.member")
        self.autoscaling_backend.resume_processes(
            autoscaling_group_name, scaling_processes
        )
        template = self.response_template(RESUME_PROCESSES_TEMPLATE)
        return template.render()

    def set_instance_protection(self) -> str:
        group_name = self._get_param("AutoScalingGroupName")
        instance_ids = self._get_multi_param("InstanceIds.member")
        protected_from_scale_in = self._get_bool_param("ProtectedFromScaleIn")
        self.autoscaling_backend.set_instance_protection(
            group_name, instance_ids, protected_from_scale_in
        )
        template = self.response_template(SET_INSTANCE_PROTECTION_TEMPLATE)
        return template.render()

    def terminate_instance_in_auto_scaling_group(self) -> str:
        instance_id = self._get_param("InstanceId")
        should_decrement_string = self._get_param("ShouldDecrementDesiredCapacity")
        if should_decrement_string == "true":
            should_decrement = True
        else:
            should_decrement = False
        (
            instance,
            original_size,
            desired_capacity,
        ) = self.autoscaling_backend.terminate_instance(instance_id, should_decrement)
        template = self.response_template(TERMINATE_INSTANCES_TEMPLATE)
        return template.render(
            instance=instance,
            should_decrement=should_decrement,
            original_size=original_size,
            desired_capacity=desired_capacity,
            timestamp=iso_8601_datetime_with_milliseconds(),
        )

    def describe_tags(self) -> str:
        filters = self._get_params().get("Filters", [])
        tags = self.autoscaling_backend.describe_tags(filters=filters)
        template = self.response_template(DESCRIBE_TAGS_TEMPLATE)
        return template.render(tags=tags, next_token=None)

    def enable_metrics_collection(self) -> str:
        group_name = self._get_param("AutoScalingGroupName")
        metrics = self._get_params().get("Metrics")
        self.autoscaling_backend.enable_metrics_collection(group_name, metrics)  # type: ignore[arg-type]
        template = self.response_template(ENABLE_METRICS_COLLECTION_TEMPLATE)
        return template.render()

    def put_warm_pool(self) -> str:
        params = self._get_params()
        group_name = params.get("AutoScalingGroupName")
        max_capacity = params.get("MaxGroupPreparedCapacity")
        min_size = params.get("MinSize")
        pool_state = params.get("PoolState")
        instance_reuse_policy = params.get("InstanceReusePolicy")
        self.autoscaling_backend.put_warm_pool(
            group_name=group_name,  # type: ignore[arg-type]
            max_capacity=max_capacity,
            min_size=min_size,
            pool_state=pool_state,
            instance_reuse_policy=instance_reuse_policy,
        )
        template = self.response_template(PUT_WARM_POOL_TEMPLATE)
        return template.render()

    def describe_warm_pool(self) -> str:
        group_name = self._get_param("AutoScalingGroupName")
        warm_pool = self.autoscaling_backend.describe_warm_pool(group_name=group_name)
        template = self.response_template(DESCRIBE_WARM_POOL_TEMPLATE)
        return template.render(pool=warm_pool)

    def delete_warm_pool(self) -> str:
        group_name = self._get_param("AutoScalingGroupName")
        self.autoscaling_backend.delete_warm_pool(group_name=group_name)
        template = self.response_template(DELETE_WARM_POOL_TEMPLATE)
        return template.render()


CREATE_LAUNCH_CONFIGURATION_TEMPLATE = """<CreateLaunchConfigurationResponse xmlns="http://autoscaling.amazonaws.com/doc/2011-01-01/">
<ResponseMetadata>
   <RequestId>7c6e177f-f082-11e1-ac58-3714bEXAMPLE</RequestId>
</ResponseMetadata>
</CreateLaunchConfigurationResponse>"""

DESCRIBE_LAUNCH_CONFIGURATIONS_TEMPLATE = """<DescribeLaunchConfigurationsResponse xmlns="http://autoscaling.amazonaws.com/doc/2011-01-01/">
  <DescribeLaunchConfigurationsResult>
    <LaunchConfigurations>
      {% for launch_configuration in launch_configurations %}
        <member>
          <AssociatePublicIpAddress>{{ 'true' if launch_configuration.associate_public_ip_address else 'false' }}</AssociatePublicIpAddress>
          {% if launch_configuration.classic_link_vpc_id %}
          <ClassicLinkVPCId>{{ launch_configuration.classic_link_vpc_id }}</ClassicLinkVPCId>
          {% endif %}
          {% if launch_configuration.classic_link_vpc_security_groups %}
          <ClassicLinkVPCSecurityGroups>
            {% for sg in launch_configuration.classic_link_vpc_security_groups %}
            <member>{{ sg }}</member>
            {% endfor %}
          </ClassicLinkVPCSecurityGroups>
          {% endif %}
          <SecurityGroups>
            {% for security_group in launch_configuration.security_groups %}
              <member>{{ security_group }}</member>
            {% endfor %}
          </SecurityGroups>
          <CreatedTime>2013-01-21T23:04:42.200Z</CreatedTime>
          {% if launch_configuration.kernel_id %}
          <KernelId>{{ launch_configuration.kernel_id }}</KernelId>
          {% else %}
          <KernelId/>
          {% endif %}
          {% if launch_configuration.instance_profile_name %}
            <IamInstanceProfile>{{ launch_configuration.instance_profile_name }}</IamInstanceProfile>
          {% endif %}
          <LaunchConfigurationName>{{ launch_configuration.name }}</LaunchConfigurationName>
          {% if launch_configuration.user_data %}
            <UserData>{{ launch_configuration.user_data }}</UserData>
          {% else %}
            <UserData/>
          {% endif %}
          <InstanceType>{{ launch_configuration.instance_type }}</InstanceType>
          <LaunchConfigurationARN>{{ launch_configuration.arn }}</LaunchConfigurationARN>
          {% if launch_configuration.block_device_mappings %}
            <BlockDeviceMappings>
            {% for mount_point, mapping in launch_configuration.block_device_mappings.items() %}
              <member>
                <DeviceName>{{ mount_point }}</DeviceName>
                {% if mapping.ephemeral_name %}
                <VirtualName>{{ mapping.ephemeral_name }}</VirtualName>
                {% elif mapping.no_device %}
                <NoDevice>true</NoDevice>
                {% else %}
                <Ebs>
                {% if mapping.snapshot_id %}
                  <SnapshotId>{{ mapping.snapshot_id }}</SnapshotId>
                {% endif %}
                {% if mapping.size %}
                  <VolumeSize>{{ mapping.size }}</VolumeSize>
                {% endif %}
                {% if mapping.iops %}
                  <Iops>{{ mapping.iops }}</Iops>
                {% endif %}
                {% if mapping.throughput %}
                  <Throughput>{{ mapping.throughput }}</Throughput>
                {% endif %}
                {% if mapping.delete_on_termination is not none %}
                  <DeleteOnTermination>{{ mapping.delete_on_termination }}</DeleteOnTermination>
                {% endif %}
                {% if mapping.volume_type %}
                  <VolumeType>{{ mapping.volume_type }}</VolumeType>
                {% endif %}
                  {% if mapping.encrypted %}
                  <Encrypted>{{ mapping.encrypted }}</Encrypted>
                  {% endif %}
                </Ebs>
                {% endif %}
              </member>
            {% endfor %}
            </BlockDeviceMappings>
          {% else %}
            <BlockDeviceMappings/>
          {% endif %}
          <ImageId>{{ launch_configuration.image_id }}</ImageId>
          {% if launch_configuration.key_name %}
            <KeyName>{{ launch_configuration.key_name }}</KeyName>
          {% else %}
            <KeyName/>
          {% endif %}
          {% if launch_configuration.ramdisk_id %}
          <RamdiskId>{{ launch_configuration.ramdisk_id }}</RamdiskId>
          {% else %}
          <RamdiskId/>
          {% endif %}
          <EbsOptimized>{{ launch_configuration.ebs_optimized }}</EbsOptimized>
          <InstanceMonitoring>
            <Enabled>{{ launch_configuration.instance_monitoring_enabled }}</Enabled>
          </InstanceMonitoring>
          {% if launch_configuration.spot_price %}
            <SpotPrice>{{ launch_configuration.spot_price }}</SpotPrice>
          {% endif %}
          {% if launch_configuration.metadata_options %}
          <MetadataOptions>
            <HttpTokens>{{ launch_configuration.metadata_options.get("HttpTokens") }}</HttpTokens>
            <HttpPutResponseHopLimit>{{ launch_configuration.metadata_options.get("HttpPutResponseHopLimit") }}</HttpPutResponseHopLimit>
            <HttpEndpoint>{{ launch_configuration.metadata_options.get("HttpEndpoint") }}</HttpEndpoint>
          </MetadataOptions>
          {% endif %}
        </member>
      {% endfor %}
    </LaunchConfigurations>
    {% if next_token %}
    <NextToken>{{ next_token }}</NextToken>
    {% endif %}
  </DescribeLaunchConfigurationsResult>
  <ResponseMetadata>
    <RequestId>d05a22f8-b690-11e2-bf8e-2113fEXAMPLE</RequestId>
  </ResponseMetadata>
</DescribeLaunchConfigurationsResponse>"""

DELETE_LAUNCH_CONFIGURATION_TEMPLATE = """<DeleteLaunchConfigurationResponse xmlns="http://autoscaling.amazonaws.com/doc/2011-01-01/">
  <ResponseMetadata>
    <RequestId>7347261f-97df-11e2-8756-35eEXAMPLE</RequestId>
  </ResponseMetadata>
</DeleteLaunchConfigurationResponse>"""

CREATE_AUTOSCALING_GROUP_TEMPLATE = """<CreateAutoScalingGroupResponse xmlns="http://autoscaling.amazonaws.com/doc/2011-01-01/">
<ResponseMetadata>
<RequestId>8d798a29-f083-11e1-bdfb-cb223EXAMPLE</RequestId>
</ResponseMetadata>
</CreateAutoScalingGroupResponse>"""

PUT_SCHEDULED_UPDATE_GROUP_ACTION_TEMPLATE = """<PutScheduledUpdateGroupActionResponse xmlns="http://autoscaling.amazonaws.com/doc/2011-01-01/">
<ResponseMetadata>
<RequestId></RequestId>
</ResponseMetadata>
</PutScheduledUpdateGroupActionResponse>"""


BATCH_PUT_SCHEDULED_UPDATE_GROUP_ACTION_TEMPLATE = """<BatchPutScheduledUpdateGroupActionResponse xmlns="http://autoscaling.amazonaws.com/doc/2011-01-01/">
  <BatchPutScheduledUpdateGroupActionResult>
    <FailedScheduledUpdateGroupActions>
      {% for failed_action in failed_actions %}
      <member>
        <ScheduledActionName>{{ failed_action.scheduled_action_name }}</ScheduledActionName>
        {% if failed_action.error_code %}
        <ErrorCode>{{ failed_action.error_code }}</ErrorCode>
        {% endif %}
        {% if failed_action.error_message %}
        <ErrorMessage>{{ failed_action.error_message }}</ErrorMessage>
        {% endif %}
      </member>
      {% endfor %}
    </FailedScheduledUpdateGroupActions>
  </BatchPutScheduledUpdateGroupActionResult>
  <ResponseMetadata>
    <RequestId></RequestId>
  </ResponseMetadata>
</BatchPutScheduledUpdateGroupActionResponse>"""

DESCRIBE_SCHEDULED_ACTIONS = """<DescribeScheduledActionsResponse xmlns="http://autoscaling.amazonaws.com/doc/2011-01-01/">
  <DescribeScheduledActionsResult>
    <ScheduledUpdateGroupActions>
      {% for scheduled_action in scheduled_actions %}
      <member>
        <AutoScalingGroupName>{{ scheduled_action.name }}</AutoScalingGroupName>
        <ScheduledActionName>{{ scheduled_action.scheduled_action_name }}</ScheduledActionName>
        {% if scheduled_action.start_time %}
        <StartTime>{{ scheduled_action.start_time }}</StartTime>
        {% endif %}
        {% if scheduled_action.end_time %}
        <EndTime>{{ scheduled_action.end_time }}</EndTime>
        {% endif %}
        {% if scheduled_action.recurrence %}
        <Recurrence>{{ scheduled_action.recurrence }}</Recurrence>
        {% endif %}
        {% if scheduled_action.min_size is not none %}
        <MinSize>{{ scheduled_action.min_size }}</MinSize>
        {% endif %}
        {% if scheduled_action.max_size is not none %}
        <MaxSize>{{ scheduled_action.max_size }}</MaxSize>
        {% endif %}
        {% if scheduled_action.desired_capacity is not none %}
        <DesiredCapacity>{{ scheduled_action.desired_capacity }}</DesiredCapacity>
        {% endif %}
        {% if scheduled_action.timezone %}
        <TimeZone>{{ scheduled_action.timezone }}</TimeZone>
        {% endif %}
      </member>
      {% endfor %}
    </ScheduledUpdateGroupActions>
  </DescribeScheduledActionsResult>
</DescribeScheduledActionsResponse>
"""

DELETE_SCHEDULED_ACTION_TEMPLATE = """<DeleteScheduledActionResponse xmlns="http://autoscaling.amazonaws.com/doc/2011-01-01/">
<ResponseMetadata>
    <RequestId>70a76d42-9665-11e2-9fdf-211deEXAMPLE</RequestId>
  </ResponseMetadata>
</DeleteScheduledActionResponse>"""

BATCH_DELETE_SCHEDULED_ACTION_TEMPLATE = """<BatchDeleteScheduledActionResponse xmlns="http://autoscaling.amazonaws.com/doc/2011-01-01/">
  <BatchDeleteScheduledActionResult>
    <FailedScheduledActions>
      {% for failed_action in failed_actions %}
      <member>
        <ScheduledActionName>{{ failed_action.scheduled_action_name }}</ScheduledActionName>
        {% if failed_action.error_code %}
        <ErrorCode>{{ failed_action.error_code }}</ErrorCode>
        {% endif %}
        {% if failed_action.error_message %}
        <ErrorMessage>{{ failed_action.error_message }}</ErrorMessage>
        {% endif %}
      </member>
      {% endfor %}
    </FailedScheduledActions>
  </BatchDeleteScheduledActionResult>
  <ResponseMetadata>
    <RequestId></RequestId>
  </ResponseMetadata>
</BatchDeleteScheduledActionResponse>"""

ATTACH_LOAD_BALANCER_TARGET_GROUPS_TEMPLATE = """<AttachLoadBalancerTargetGroupsResponse xmlns="http://autoscaling.amazonaws.com/doc/2011-01-01/">
<AttachLoadBalancerTargetGroupsResult>
</AttachLoadBalancerTargetGroupsResult>
<ResponseMetadata>
<RequestId></RequestId>
</ResponseMetadata>
</AttachLoadBalancerTargetGroupsResponse>"""

ATTACH_INSTANCES_TEMPLATE = """<AttachInstancesResponse xmlns="http://autoscaling.amazonaws.com/doc/2011-01-01/">
<AttachInstancesResult>
</AttachInstancesResult>
<ResponseMetadata>
<RequestId></RequestId>
</ResponseMetadata>
</AttachInstancesResponse>"""

DESCRIBE_LOAD_BALANCER_TARGET_GROUPS = """<DescribeLoadBalancerTargetGroupsResponse xmlns="http://autoscaling.amazonaws.com/doc/2011-01-01/">
<DescribeLoadBalancerTargetGroupsResult>
  <LoadBalancerTargetGroups>
  {% for arn in target_group_arns %}
    <member>
      <LoadBalancerTargetGroupARN>{{ arn }}</LoadBalancerTargetGroupARN>
      <State>Added</State>
    </member>
  {% endfor %}
  </LoadBalancerTargetGroups>
</DescribeLoadBalancerTargetGroupsResult>
<ResponseMetadata>
<RequestId></RequestId>
</ResponseMetadata>
</DescribeLoadBalancerTargetGroupsResponse>"""

DETACH_INSTANCES_TEMPLATE = """<DetachInstancesResponse xmlns="http://autoscaling.amazonaws.com/doc/2011-01-01/">
<DetachInstancesResult>
  <Activities>
    {% for instance in detached_instances %}
    <member>
      <ActivityId>5091cb52-547a-47ce-a236-c9ccbc2cb2c9EXAMPLE</ActivityId>
      <AutoScalingGroupName>{{ group_name }}</AutoScalingGroupName>
      <Cause>
      At 2017-10-15T15:55:21Z instance {{ instance.instance.id }} was detached in response to a user request.
      </Cause>
      <Description>Detaching EC2 instance: {{ instance.instance.id }}</Description>
      <StartTime>2017-10-15T15:55:21Z</StartTime>
      <EndTime>2017-10-15T15:55:21Z</EndTime>
      <StatusCode>InProgress</StatusCode>
      <StatusMessage>InProgress</StatusMessage>
      <Progress>50</Progress>
      <Details>details</Details>
    </member>
    {% endfor %}
  </Activities>
</DetachInstancesResult>
<ResponseMetadata>
<RequestId></RequestId>
</ResponseMetadata>
</DetachInstancesResponse>"""

DETACH_LOAD_BALANCER_TARGET_GROUPS_TEMPLATE = """<DetachLoadBalancerTargetGroupsResponse xmlns="http://autoscaling.amazonaws.com/doc/2011-01-01/">
<DetachLoadBalancerTargetGroupsResult>
</DetachLoadBalancerTargetGroupsResult>
<ResponseMetadata>
<RequestId></RequestId>
</ResponseMetadata>
</DetachLoadBalancerTargetGroupsResponse>"""

DESCRIBE_AUTOSCALING_GROUPS_TEMPLATE = """<DescribeAutoScalingGroupsResponse xmlns="http://autoscaling.amazonaws.com/doc/2011-01-01/">
<DescribeAutoScalingGroupsResult>
    <AutoScalingGroups>
      {% for group in groups %}
      <member>
        <Tags>
          {% for tag in group.tags %}
          <member>
            <ResourceType>{{ tag.resource_type or tag.ResourceType }}</ResourceType>
            <ResourceId>{{ tag.resource_id or tag.ResourceId }}</ResourceId>
            <PropagateAtLaunch>{{ tag.propagate_at_launch or tag.PropagateAtLaunch }}</PropagateAtLaunch>
            <Key>{{ tag.key or tag.Key }}</Key>
            <Value>{{ tag.value or tag.Value }}</Value>
          </member>
          {% endfor %}
        </Tags>
        <SuspendedProcesses>
          {% for suspended_process in group.suspended_processes %}
          <member>
            <ProcessName>{{suspended_process}}</ProcessName>
            <SuspensionReason></SuspensionReason>
          </member>
          {% endfor %}
        </SuspendedProcesses>
        <AutoScalingGroupName>{{ group.name }}</AutoScalingGroupName>
        <HealthCheckType>{{ group.health_check_type }}</HealthCheckType>
        <CreatedTime>{{ group.created_time }}</CreatedTime>
        {% if group.launch_config_name %}
        <LaunchConfigurationName>{{ group.launch_config_name }}</LaunchConfigurationName>
        {% elif group.mixed_instances_policy %}
        <MixedInstancesPolicy>
          <LaunchTemplate>
            <LaunchTemplateSpecification>
              <LaunchTemplateId>{{ group.launch_template.id }}</LaunchTemplateId>
              <Version>{{ group.launch_template_version }}</Version>
              <LaunchTemplateName>{{ group.launch_template.name }}</LaunchTemplateName>
            </LaunchTemplateSpecification>
            {% if group.mixed_instances_policy.get("LaunchTemplate", {}).get("Overrides", []) %}
            <Overrides>
              {% for member in group.mixed_instances_policy.get("LaunchTemplate", {}).get("Overrides", []) %}
              <member>
                {% if member.get("InstanceType") %}
                <InstanceType>{{ member.get("InstanceType") }}</InstanceType>
                {% endif %}
                {% if member.get("WeightedCapacity") %}
                <WeightedCapacity>{{ member.get("WeightedCapacity") }}</WeightedCapacity>
                {% endif %}
              </member>
              {% endfor %}
            </Overrides>
            {% endif %}
          </LaunchTemplate>
          {% if group.mixed_instances_policy.get("InstancesDistribution") %}
          <InstancesDistribution>
            {% if group.mixed_instances_policy.get("InstancesDistribution").get("OnDemandAllocationStrategy") %}
            <OnDemandAllocationStrategy>{{ group.mixed_instances_policy.get("InstancesDistribution").get("OnDemandAllocationStrategy") }}</OnDemandAllocationStrategy>
            {% endif %}
            {% if group.mixed_instances_policy.get("InstancesDistribution").get("OnDemandBaseCapacity") %}
            <OnDemandBaseCapacity>{{ group.mixed_instances_policy.get("InstancesDistribution").get("OnDemandBaseCapacity") }}</OnDemandBaseCapacity>
            {% endif %}
            {% if group.mixed_instances_policy.get("InstancesDistribution").get("OnDemandPercentageAboveBaseCapacity") %}
            <OnDemandPercentageAboveBaseCapacity>{{ group.mixed_instances_policy.get("InstancesDistribution").get("OnDemandPercentageAboveBaseCapacity") }}</OnDemandPercentageAboveBaseCapacity>
            {% endif %}
            {% if group.mixed_instances_policy.get("InstancesDistribution").get("SpotAllocationStrategy") %}
            <SpotAllocationStrategy>{{ group.mixed_instances_policy.get("InstancesDistribution").get("SpotAllocationStrategy") }}</SpotAllocationStrategy>
            {% endif %}
            {% if group.mixed_instances_policy.get("InstancesDistribution").get("SpotInstancePools") %}
            <SpotInstancePools>{{ group.mixed_instances_policy.get("InstancesDistribution").get("SpotInstancePools") }}</SpotInstancePools>
            {% endif %}
            {% if group.mixed_instances_policy.get("InstancesDistribution").get("SpotMaxPrice") %}
            <SpotMaxPrice>{{ group.mixed_instances_policy.get("InstancesDistribution").get("SpotMaxPrice") }}</SpotMaxPrice>
            {% endif %}
          </InstancesDistribution>
          {% endif %}
        </MixedInstancesPolicy>
        {% elif group.launch_template %}
        <LaunchTemplate>
          <LaunchTemplateId>{{ group.launch_template.id }}</LaunchTemplateId>
          {% if group.provided_launch_template_version %}}
          <Version>{{ group.provided_launch_template_version }}</Version>
          {% endif %}
          <LaunchTemplateName>{{ group.launch_template.name }}</LaunchTemplateName>
        </LaunchTemplate>
        {% endif %}
        <Instances>
          {% for instance_state in group.instance_states %}
          <member>
            <HealthStatus>{{ instance_state.health_status }}</HealthStatus>
            <AvailabilityZone>{{ instance_state.instance.placement }}</AvailabilityZone>
            <InstanceId>{{ instance_state.instance.id }}</InstanceId>
            <InstanceType>{{ instance_state.instance.instance_type }}</InstanceType>
            {% if group.launch_config_name %}
            <LaunchConfigurationName>{{ group.launch_config_name }}</LaunchConfigurationName>
            {% elif group.launch_template %}
            <LaunchTemplate>
              <LaunchTemplateId>{{ group.launch_template.id }}</LaunchTemplateId>
              <Version>{{ group.launch_template_version }}</Version>
              <LaunchTemplateName>{{ group.launch_template.name }}</LaunchTemplateName>
            </LaunchTemplate>
            {% endif %}
            <LifecycleState>{{ instance_state.lifecycle_state }}</LifecycleState>
            <ProtectedFromScaleIn>{{ instance_state.protected_from_scale_in|string|lower }}</ProtectedFromScaleIn>
          </member>
          {% endfor %}
        </Instances>
        <DesiredCapacity>{{ group.desired_capacity }}</DesiredCapacity>
        <CapacityRebalance>{{ 'true' if group.capacity_rebalance else 'false' }}</CapacityRebalance>
        <AvailabilityZones>
          {% for availability_zone in group.availability_zones %}
          <member>{{ availability_zone }}</member>
          {% endfor %}
        </AvailabilityZones>
        {% if group.load_balancers %}
          <LoadBalancerNames>
          {% for load_balancer in group.load_balancers %}
            <member>{{ load_balancer }}</member>
          {% endfor %}
          </LoadBalancerNames>
        {% else %}
          <LoadBalancerNames/>
        {% endif %}
        {% if group.target_group_arns %}
          <TargetGroupARNs>
          {% for target_group_arn in group.target_group_arns %}
            <member>{{ target_group_arn }}</member>
          {% endfor %}
          </TargetGroupARNs>
        {% else %}
          <TargetGroupARNs/>
        {% endif %}
        <MinSize>{{ group.min_size }}</MinSize>
        {% if group.vpc_zone_identifier %}
          <VPCZoneIdentifier>{{ group.vpc_zone_identifier }}</VPCZoneIdentifier>
        {% else %}
          <VPCZoneIdentifier/>
        {% endif %}
        <HealthCheckGracePeriod>{{ group.health_check_period }}</HealthCheckGracePeriod>
        <DefaultCooldown>{{ group.default_cooldown }}</DefaultCooldown>
        <AutoScalingGroupARN>{{ group.arn }}</AutoScalingGroupARN>
        {% if group.termination_policies %}
        <TerminationPolicies>
          {% for policy in group.termination_policies %}
          <member>{{ policy }}</member>
          {% endfor %}
        </TerminationPolicies>
        {% else %}
        <TerminationPolicies/>
        {% endif %}
        <MaxSize>{{ group.max_size }}</MaxSize>
        {% if group.placement_group %}
        <PlacementGroup>{{ group.placement_group }}</PlacementGroup>
        {% endif %}
        <NewInstancesProtectedFromScaleIn>{{ group.new_instances_protected_from_scale_in|string|lower }}</NewInstancesProtectedFromScaleIn>
        {% if group.metrics %}
        <EnabledMetrics>
          {% for met in group.metrics %}
          <member>
          <Metric>{{ met }}</Metric>
          <Granularity>1Minute</Granularity>
          </member>
          {% endfor %}
        </EnabledMetrics>
        {% endif %}
        <ServiceLinkedRoleARN>{{ group.service_linked_role }}</ServiceLinkedRoleARN>
        {% if group.warm_pool %}
        <WarmPoolConfiguration>
          <MaxGroupPreparedCapacity>{{ group.warm_pool.max_capacity }}</MaxGroupPreparedCapacity>
          <MinSize>{{ group.warm_pool.min_size or 0 }}</MinSize>
          {% if group.warm_pool.pool_state %}
          <PoolState>{{ group.warm_pool.pool_state }}</PoolState>
          {% endif %}
          <InstanceReusePolicy>
            <ReuseOnScaleIn>{{ 'true' if group.warm_pool.instance_reuse_policy["ReuseOnScaleIn"] else 'false' }}</ReuseOnScaleIn>
          </InstanceReusePolicy>
        </WarmPoolConfiguration>
        {% endif %}
      </member>
      {% endfor %}
    </AutoScalingGroups>
    {% if next_token %}
    <NextToken>{{ next_token }}</NextToken>
    {% endif %}
  </DescribeAutoScalingGroupsResult>
  <ResponseMetadata>
    <RequestId>0f02a07d-b677-11e2-9eb0-dd50EXAMPLE</RequestId>
  </ResponseMetadata>
</DescribeAutoScalingGroupsResponse>"""

UPDATE_AUTOSCALING_GROUP_TEMPLATE = """<UpdateAutoScalingGroupResponse xmlns="http://autoscaling.amazonaws.com/doc/2011-01-01/">
  <ResponseMetadata>
    <RequestId>adafead0-ab8a-11e2-ba13-ab0ccEXAMPLE</RequestId>
  </ResponseMetadata>
</UpdateAutoScalingGroupResponse>"""

DELETE_AUTOSCALING_GROUP_TEMPLATE = """<DeleteAutoScalingGroupResponse xmlns="http://autoscaling.amazonaws.com/doc/2011-01-01/">
  <ResponseMetadata>
    <RequestId>70a76d42-9665-11e2-9fdf-211deEXAMPLE</RequestId>
  </ResponseMetadata>
</DeleteAutoScalingGroupResponse>"""

DESCRIBE_SCALING_ACTIVITIES_TEMPLATE = """<DescribeScalingActivitiesResponse xmlns="http://autoscaling.amazonaws.com/doc/2011-01-01/">
<DescribeScalingActivitiesResult>
</DescribeScalingActivitiesResult>
<ResponseMetadata>
<RequestId></RequestId>
</ResponseMetadata>
</DescribeScalingActivitiesResponse>"""

DESCRIBE_AUTOSCALING_INSTANCES_TEMPLATE = """<DescribeAutoScalingInstancesResponse xmlns="http://autoscaling.amazonaws.com/doc/2011-01-01/">
  <DescribeAutoScalingInstancesResult>
    <AutoScalingInstances>
      {% for instance_state in instance_states %}
      <member>
        <HealthStatus>{{ instance_state.health_status }}</HealthStatus>
        <AutoScalingGroupName>{{ instance_state.instance.autoscaling_group.name }}</AutoScalingGroupName>
        <AvailabilityZone>{{ instance_state.instance.placement }}</AvailabilityZone>
        <InstanceId>{{ instance_state.instance.id }}</InstanceId>
        <InstanceType>{{ instance_state.instance.instance_type }}</InstanceType>
        {% if instance_state.instance.autoscaling_group.launch_config_name %}
        <LaunchConfigurationName>{{ instance_state.instance.autoscaling_group.launch_config_name }}</LaunchConfigurationName>
        {% elif instance_state.instance.autoscaling_group.launch_template %}
        <LaunchTemplate>
          <LaunchTemplateId>{{ instance_state.instance.autoscaling_group.launch_template.id }}</LaunchTemplateId>
          <Version>{{ instance_state.instance.autoscaling_group.launch_template_version }}</Version>
          <LaunchTemplateName>{{ instance_state.instance.autoscaling_group.launch_template.name }}</LaunchTemplateName>
        </LaunchTemplate>
        {% endif %}
        <LifecycleState>{{ instance_state.lifecycle_state }}</LifecycleState>
        <ProtectedFromScaleIn>{{ instance_state.protected_from_scale_in|string|lower }}</ProtectedFromScaleIn>
      </member>
      {% endfor %}
    </AutoScalingInstances>
  </DescribeAutoScalingInstancesResult>
  <ResponseMetadata>
    <RequestId>df992dc3-b72f-11e2-81e1-750aa6EXAMPLE</RequestId>
  </ResponseMetadata>
</DescribeAutoScalingInstancesResponse>"""

CREATE_LIFECYLE_HOOK_TEMPLATE = """<PutLifecycleHookResponse xmlns="http://autoscaling.amazonaws.com/doc/2011-01-01/">
  <PutLifecycleHookResult/>
  <ResponseMetadata>
    <RequestId>3cfc6fef-c08b-11e2-a697-2922EXAMPLE</RequestId>
  </ResponseMetadata>
</PutLifecycleHookResponse>"""

DESCRIBE_LIFECYCLE_HOOKS_TEMPLATE = """<DescribeLifecycleHooksResponse xmlns="http://autoscaling.amazonaws.com/doc/2011-01-01/">
  <DescribeLifecycleHooksResult>
    <LifecycleHooks>
      {% for lifecycle_hook in lifecycle_hooks %}
        <member>
          <AutoScalingGroupName>{{ lifecycle_hook.as_name }}</AutoScalingGroupName>
          <RoleARN>arn:aws:iam::1234567890:role/my-auto-scaling-role</RoleARN>
          <LifecycleTransition>{{ lifecycle_hook.transition }}</LifecycleTransition>
          <GlobalTimeout>172800</GlobalTimeout>
          <LifecycleHookName>{{ lifecycle_hook.name }}</LifecycleHookName>
          <HeartbeatTimeout>{{ lifecycle_hook.timeout }}</HeartbeatTimeout>
          <DefaultResult>{{ lifecycle_hook.result }}</DefaultResult>
          <NotificationTargetARN>arn:aws:sqs:us-east-1:123456789012:my-queue</NotificationTargetARN>
        </member>
      {% endfor %}
    </LifecycleHooks>
  </DescribeLifecycleHooksResult>
  <ResponseMetadata>
    <RequestId>ec3bffad-b739-11e2-b38d-15fbEXAMPLE</RequestId>
  </ResponseMetadata>
</DescribeLifecycleHooksResponse>"""

DELETE_LIFECYCLE_HOOK_TEMPLATE = """<DeleteLifecycleHookResponse xmlns="http://autoscaling.amazonaws.com/doc/2011-01-01/">
  <DeleteLifecycleHookResult>
  </DeleteLifecycleHookResult>
  <ResponseMetadata>
    <RequestId>70a76d42-9665-11e2-9fdf-211deEXAMPLE</RequestId>
  </ResponseMetadata>
</DeleteLifecycleHookResponse>"""

CREATE_SCALING_POLICY_TEMPLATE = """<PutScalingPolicyResponse xmlns="http://autoscaling.amazonaws.com/doc/2011-01-01/">
  <PutScalingPolicyResult>
    <PolicyARN>arn:aws:autoscaling:us-east-1:803981987763:scalingPolicy:b0dcf5e8
-02e6-4e31-9719-0675d0dc31ae:autoScalingGroupName/my-test-asg:policyName/my-scal
eout-policy</PolicyARN>
  </PutScalingPolicyResult>
  <ResponseMetadata>
    <RequestId>3cfc6fef-c08b-11e2-a697-2922EXAMPLE</RequestId>
  </ResponseMetadata>
</PutScalingPolicyResponse>"""

DESCRIBE_SCALING_POLICIES_TEMPLATE = """<DescribePoliciesResponse xmlns="http://autoscaling.amazonaws.com/doc/2011-01-01/">
  <DescribePoliciesResult>
    <ScalingPolicies>
      {% for policy in policies %}
      <member>
        <PolicyARN>{{ policy.arn }}</PolicyARN>
        {% if policy.adjustment_type %}
        <AdjustmentType>{{ policy.adjustment_type }}</AdjustmentType>
        {% endif %}
        {% if policy.scaling_adjustment %}
        <ScalingAdjustment>{{ policy.scaling_adjustment }}</ScalingAdjustment>
        {% endif %}
        {% if policy.min_adjustment_magnitude %}
        <MinAdjustmentMagnitude>{{ policy.min_adjustment_magnitude }}</MinAdjustmentMagnitude>
        {% endif %}
        <PolicyName>{{ policy.name }}</PolicyName>
        <PolicyType>{{ policy.policy_type }}</PolicyType>
        <MetricAggregationType>{{ policy.metric_aggregation_type }}</MetricAggregationType>
        <AutoScalingGroupName>{{ policy.as_name }}</AutoScalingGroupName>
        {% if policy.policy_type == 'SimpleScaling' %}
        <Cooldown>{{ policy.cooldown }}</Cooldown>
        {% endif %}
        {% if policy.policy_type == 'TargetTrackingScaling' %}
        <TargetTrackingConfiguration>
            {% if policy.target_tracking_config.get("PredefinedMetricSpecification") %}
            <PredefinedMetricSpecification>
                <PredefinedMetricType>{{ policy.target_tracking_config.get("PredefinedMetricSpecification", {}).get("PredefinedMetricType", "") }}</PredefinedMetricType>
                {% if policy.target_tracking_config.get("PredefinedMetricSpecification", {}).get("ResourceLabel") %}
                <ResourceLabel>{{ policy.target_tracking_config.get("PredefinedMetricSpecification", {}).get("ResourceLabel") }}</ResourceLabel>
                {% endif %}
            </PredefinedMetricSpecification>
            {% endif %}
            {% if policy.target_tracking_config.get("CustomizedMetricSpecification") %}
            <CustomizedMetricSpecification>
              <MetricName>{{ policy.target_tracking_config["CustomizedMetricSpecification"].get("MetricName") }}</MetricName>
              <Namespace>{{ policy.target_tracking_config["CustomizedMetricSpecification"].get("Namespace") }}</Namespace>
              <Dimensions>
                {% for dim in policy.target_tracking_config["CustomizedMetricSpecification"].get("Dimensions", []) %}
                <member>
                  <Name>{{ dim.get("Name") }}</Name>
                  <Value>{{ dim.get("Value") }}</Value>
                </member>
                {% endfor %}
              </Dimensions>
              <Statistic>{{ policy.target_tracking_config["CustomizedMetricSpecification"].get("Statistic") }}</Statistic>
              {% if policy.target_tracking_config["CustomizedMetricSpecification"].get("Unit") %}
              <Unit>{{ policy.target_tracking_config["CustomizedMetricSpecification"].get("Unit") }}</Unit>
              {% endif %}
              {% if policy.target_tracking_config["CustomizedMetricSpecification"].get("Metrics") %}
              <Metrics>
                {% for metric in policy.target_tracking_config["CustomizedMetricSpecification"].get("Metrics", []) %}
                <member>
                  <Id>{{ metric.get("Id") }}</Id>
                  {% if metric.get("MetricStat") is none %}
                  <Expression>{{ metric.get("Expression") }}</Expression>
                  {% endif %}
                  {% if metric.get("Expression") is none %}
                  <MetricStat>
                    <Metric>
                      <Namespace>{{ metric.get("MetricStat", {}).get("Metric", {}).get("Namespace") }}</Namespace>
                      <MetricName>{{ metric.get("MetricStat", {}).get("Metric", {}).get("MetricName") }}</MetricName>
                      <Dimensions>
                      {% for dim in metric.get("MetricStat", {}).get("Metric", {}).get("Dimensions", []) %}
                        <member>
                          <Name>{{ dim.get("Name") }}</Name>
                          <Value>{{ dim.get("Value") }}</Value>
                        </member>
                      {% endfor %}
                      </Dimensions>
                    </Metric>
                    <Stat>{{ metric.get("MetricStat", {}).get("Stat") }}</Stat>
                    <Unit>{{ metric.get("MetricStat", {}).get("Unit") }}</Unit>
                  </MetricStat>
                  {% endif %}
                  <Label>{{ metric.get("Label") }}</Label>
                  <ReturnData>{{ 'true' if metric.get("ReturnData") is none else metric.get("ReturnData") }}</ReturnData>
                </member>
                {% endfor %}
              </Metrics>
              {% endif %}
            </CustomizedMetricSpecification>
            {% endif %}
            <TargetValue>{{ policy.target_tracking_config.get("TargetValue") }}</TargetValue>
        </TargetTrackingConfiguration>
        {% endif %}
        {% if policy.policy_type == 'StepScaling' %}
        <StepAdjustments>
        {% for step in policy.step_adjustments %}
        <member>
            {% if "MetricIntervalLowerBound" in step %}
            <MetricIntervalLowerBound>{{ step.get("MetricIntervalLowerBound") }}</MetricIntervalLowerBound>
            {% endif %}
            {% if "MetricIntervalUpperBound" in step %}
            <MetricIntervalUpperBound>{{ step.get("MetricIntervalUpperBound") }}</MetricIntervalUpperBound>
            {% endif %}
            {% if "ScalingAdjustment" in step %}
            <ScalingAdjustment>{{ step.get("ScalingAdjustment") }}</ScalingAdjustment>
            {% endif %}
        </member>
        {% endfor %}
        </StepAdjustments>
        {% endif %}
        {% if policy.estimated_instance_warmup %}
        <EstimatedInstanceWarmup>{{ policy.estimated_instance_warmup }}</EstimatedInstanceWarmup>
        {% endif %}
        {% if policy.policy_type == 'PredictiveScaling' %}
        <PredictiveScalingConfiguration>
            <MetricSpecifications>
                {% for config in policy.predictive_scaling_configuration.get("MetricSpecifications", []) %}
                <member>
                  <TargetValue>{{ config.get("TargetValue") }}</TargetValue>
                  {% if config.get("PredefinedMetricPairSpecification", {}).get("PredefinedMetricType") %}
                  <PredefinedMetricPairSpecification>
                    <PredefinedMetricType>{{ config.get("PredefinedMetricPairSpecification", {}).get("PredefinedMetricType") }}</PredefinedMetricType>
                    <ResourceLabel>{{ config.get("PredefinedMetricPairSpecification", {}).get("ResourceLabel", "") }}</ResourceLabel>
                  </PredefinedMetricPairSpecification>
                  {% endif %}
                  {% if config.get("PredefinedScalingMetricSpecification", {}).get("PredefinedMetricType") %}
                  <PredefinedScalingMetricSpecification>
                    <PredefinedMetricType>{{ config.get("PredefinedScalingMetricSpecification", {}).get("PredefinedMetricType", "") }}</PredefinedMetricType>
                    <ResourceLabel>{{ config.get("PredefinedScalingMetricSpecification", {}).get("ResourceLabel", "") }}</ResourceLabel>
                  </PredefinedScalingMetricSpecification>
                  {% endif %}
                  {% if config.get("PredefinedLoadMetricSpecification", {}).get("PredefinedMetricType") %}
                  <PredefinedLoadMetricSpecification>
                    <PredefinedMetricType>{{ config.get("PredefinedLoadMetricSpecification", {}).get("PredefinedMetricType", "") }}</PredefinedMetricType>
                    <ResourceLabel>{{ config.get("PredefinedLoadMetricSpecification", {}).get("ResourceLabel", "") }}</ResourceLabel>
                  </PredefinedLoadMetricSpecification>
                  {% endif %}
                  {% if config.get("CustomizedScalingMetricSpecification", {}).get("MetricDataQueries") %}
                  <CustomizedScalingMetricSpecification>
                    <MetricDataQueries>
                    {% for query in config.get("CustomizedScalingMetricSpecification", {}).get("MetricDataQueries", []) %}
                    <member>
                      <Id>{{ query.get("Id") }}</Id>
                      <Expression>{{ query.get("Expression") }}</Expression>
                      <MetricStat>
                        <Metric>
                          <Namespace>{{ query.get("MetricStat", {}).get("Metric", {}).get("Namespace") }}</Namespace>
                          <MetricName>{{ query.get("MetricStat", {}).get("Metric", {}).get("MetricName") }}</MetricName>
                          <Dimensions>
                          {% for dim in query.get("MetricStat", {}).get("Metric", {}).get("Dimensions", []) %}
                            <Name>{{ dim.get("Name") }}</Name>
                            <Value>{{ dim.get("Value") }}</Value>
                          {% endfor %}
                          </Dimensions>
                        </Metric>
                        <Stat>{{ query.get("MetricStat", {}).get("Stat") }}</Stat>
                        <Unit>{{ query.get("MetricStat", {}).get("Unit") }}</Unit>
                      </MetricStat>
                      <Label>{{ query.get("Label") }}</Label>
                      <ReturnData>{{ 'true' if query.get("ReturnData") else 'false' }}</ReturnData>
                    </member>
                    {% endfor %}
                    </MetricDataQueries>
                  </CustomizedScalingMetricSpecification>
                  {% endif %}
                  {% if config.get("CustomizedLoadMetricSpecification", {}).get("MetricDataQueries") %}
                  <CustomizedLoadMetricSpecification>
                    <MetricDataQueries>
                    {% for query in config.get("CustomizedLoadMetricSpecification", {}).get("MetricDataQueries", []) %}
                    <member>
                      <Id>{{ query.get("Id") }}</Id>
                      <Expression>{{ query.get("Expression") }}</Expression>
                      <MetricStat>
                        <Metric>
                          <Namespace>{{ query.get("MetricStat", {}).get("Metric", {}).get("Namespace") }}</Namespace>
                          <MetricName>{{ query.get("MetricStat", {}).get("Metric", {}).get("MetricName") }}</MetricName>
                          <Dimensions>
                          {% for dim in query.get("MetricStat", {}).get("Metric", {}).get("Dimensions", []) %}
                            <Name>{{ dim.get("Name") }}</Name>
                            <Value>{{ dim.get("Value") }}</Value>
                          {% endfor %}
                          </Dimensions>
                        </Metric>
                        <Stat>{{ query.get("MetricStat", {}).get("Stat") }}</Stat>
                        <Unit>{{ query.get("MetricStat", {}).get("Unit") }}</Unit>
                      </MetricStat>
                      <Label>{{ query.get("Label") }}</Label>
                      <ReturnData>{{ 'true' if query.get("ReturnData") else 'false' }}</ReturnData>
                    </member>
                    {% endfor %}
                    </MetricDataQueries>
                  </CustomizedLoadMetricSpecification>
                  {% endif %}
                  {% if config.get("CustomizedCapacityMetricSpecification", {}).get("MetricDataQueries") %}
                  <CustomizedCapacityMetricSpecification>
                    <MetricDataQueries>
                    {% for query in config.get("CustomizedCapacityMetricSpecification", {}).get("MetricDataQueries", []) %}
                    <member>
                      <Id>{{ query.get("Id") }}</Id>
                      <Expression>{{ query.get("Expression") }}</Expression>
                      <MetricStat>
                        <Metric>
                          <Namespace>{{ query.get("MetricStat", {}).get("Metric", {}).get("Namespace") }}</Namespace>
                          <MetricName>{{ query.get("MetricStat", {}).get("Metric", {}).get("MetricName") }}</MetricName>
                          <Dimensions>
                          {% for dim in query.get("MetricStat", {}).get("Metric", {}).get("Dimensions", []) %}
                            <Name>{{ dim.get("Name") }}</Name>
                            <Value>{{ dim.get("Value") }}</Value>
                          {% endfor %}
                          </Dimensions>
                        </Metric>
                        <Stat>{{ query.get("MetricStat", {}).get("Stat") }}</Stat>
                        <Unit>{{ query.get("MetricStat", {}).get("Unit") }}</Unit>
                      </MetricStat>
                      <Label>{{ query.get("Label") }}</Label>
                      <ReturnData>{{ 'true' if query.get("ReturnData") else 'false' }}</ReturnData>
                    </member>
                    {% endfor %}
                    </MetricDataQueries>
                  </CustomizedCapacityMetricSpecification>
                  {% endif %}
                </member>
                {% endfor %}
            </MetricSpecifications>
            {% if "Mode" in policy.predictive_scaling_configuration %}
            <Mode>{{ policy.predictive_scaling_configuration.get("Mode") }}</Mode>
            {% endif %}
            {% if "SchedulingBufferTime" in policy.predictive_scaling_configuration %}
            <SchedulingBufferTime>{{ policy.predictive_scaling_configuration.get("SchedulingBufferTime") }}</SchedulingBufferTime>
            {% endif %}
            {% if "MaxCapacityBreachBehavior" in policy.predictive_scaling_configuration %}
            <MaxCapacityBreachBehavior>{{ policy.predictive_scaling_configuration.get("MaxCapacityBreachBehavior") }}</MaxCapacityBreachBehavior>
            {% endif %}
            {% if "MaxCapacityBuffer" in policy.predictive_scaling_configuration %}
            <MaxCapacityBuffer>{{ policy.predictive_scaling_configuration.get("MaxCapacityBuffer") }}</MaxCapacityBuffer>
            {% endif %}
        </PredictiveScalingConfiguration>
        {% endif %}
        <Alarms/>
      </member>
      {% endfor %}
    </ScalingPolicies>
  </DescribePoliciesResult>
  <ResponseMetadata>
    <RequestId>ec3bffad-b739-11e2-b38d-15fbEXAMPLE</RequestId>
  </ResponseMetadata>
</DescribePoliciesResponse>"""

SET_DESIRED_CAPACITY_TEMPLATE = """<SetDesiredCapacityResponse xmlns="http://autoscaling.amazonaws.com/doc/2011-01-01/">
  <ResponseMetadata>
    <RequestId>9fb7e2db-6998-11e2-a985-57c82EXAMPLE</RequestId>
  </ResponseMetadata>
</SetDesiredCapacityResponse>"""

EXECUTE_POLICY_TEMPLATE = """<ExecuteScalingPolicyResponse xmlns="http://autoscaling.amazonaws.com/doc/2011-01-01/">
  <ResponseMetadata>
    <RequestId>70a76d42-9665-11e2-9fdf-211deEXAMPLE</RequestId>
  </ResponseMetadata>
</ExecuteScalingPolicyResponse>"""

DELETE_POLICY_TEMPLATE = """<DeleteScalingPolicyResponse xmlns="http://autoscaling.amazonaws.com/doc/2011-01-01/">
  <ResponseMetadata>
    <RequestId>70a76d42-9665-11e2-9fdf-211deEXAMPLE</RequestId>
  </ResponseMetadata>
</DeleteScalingPolicyResponse>"""

ATTACH_LOAD_BALANCERS_TEMPLATE = """<AttachLoadBalancersResponse xmlns="http://autoscaling.amazonaws.com/doc/2011-01-01/">
<AttachLoadBalancersResult></AttachLoadBalancersResult>
<ResponseMetadata>
<RequestId></RequestId>
</ResponseMetadata>
</AttachLoadBalancersResponse>"""

DESCRIBE_LOAD_BALANCERS_TEMPLATE = """<DescribeLoadBalancersResponse xmlns="http://autoscaling.amazonaws.com/doc/2011-01-01/">
<DescribeLoadBalancersResult>
  <LoadBalancers>
    {% for load_balancer in load_balancers %}
      <member>
        <LoadBalancerName>{{ load_balancer }}</LoadBalancerName>
        <State>Added</State>
      </member>
    {% endfor %}
  </LoadBalancers>
</DescribeLoadBalancersResult>
<ResponseMetadata>
<RequestId></RequestId>
</ResponseMetadata>
</DescribeLoadBalancersResponse>"""

DETACH_LOAD_BALANCERS_TEMPLATE = """<DetachLoadBalancersResponse xmlns="http://autoscaling.amazonaws.com/doc/2011-01-01/">
<DetachLoadBalancersResult></DetachLoadBalancersResult>
<ResponseMetadata>
<RequestId></RequestId>
</ResponseMetadata>
</DetachLoadBalancersResponse>"""

SUSPEND_PROCESSES_TEMPLATE = """<SuspendProcessesResponse xmlns="http://autoscaling.amazonaws.com/doc/2011-01-01/">
<ResponseMetadata>
   <RequestId>7c6e177f-f082-11e1-ac58-3714bEXAMPLE</RequestId>
</ResponseMetadata>
</SuspendProcessesResponse>"""

RESUME_PROCESSES_TEMPLATE = """<ResumeProcessesResponse xmlns="http://autoscaling.amazonaws.com/doc/2011-01-01/">
<ResponseMetadata>
   <RequestId></RequestId>
</ResponseMetadata>
</ResumeProcessesResponse>"""

SET_INSTANCE_HEALTH_TEMPLATE = """<SetInstanceHealthResponse xmlns="http://autoscaling.amazonaws.com/doc/2011-01-01/">
<SetInstanceHealthResponse></SetInstanceHealthResponse>
<ResponseMetadata>
<RequestId></RequestId>
</ResponseMetadata>
</SetInstanceHealthResponse>"""

SET_INSTANCE_PROTECTION_TEMPLATE = """<SetInstanceProtectionResponse xmlns="http://autoscaling.amazonaws.com/doc/2011-01-01/">
<SetInstanceProtectionResult></SetInstanceProtectionResult>
<ResponseMetadata>
<RequestId></RequestId>
</ResponseMetadata>
</SetInstanceProtectionResponse>"""

ENTER_STANDBY_TEMPLATE = """<EnterStandbyResponse xmlns="http://autoscaling.amazonaws.com/doc/2011-01-01/">
  <EnterStandbyResult>
    <Activities>
      {% for instance in standby_instances %}
      <member>
        <ActivityId>12345678-1234-1234-1234-123456789012</ActivityId>
        <AutoScalingGroupName>{{ group_name }}</AutoScalingGroupName>
        {% if should_decrement %}
        <Cause>At {{ timestamp }} instance {{ instance.instance.id }} was moved to standby in response to a user request, shrinking the capacity from {{ original_size }} to {{ desired_capacity }}.</Cause>
        {% else %}
        <Cause>At {{ timestamp }} instance {{ instance.instance.id }} was moved to standby in response to a user request.</Cause>
        {% endif %}
        <Description>Moving EC2 instance to Standby: {{ instance.instance.id }}</Description>
        <Progress>50</Progress>
        <StartTime>{{ timestamp }}</StartTime>
        <Details>{&quot;Subnet ID&quot;:&quot;??&quot;,&quot;Availability Zone&quot;:&quot;{{ instance.instance.placement }}&quot;}</Details>
        <StatusCode>InProgress</StatusCode>
      </member>
      {% endfor %}
    </Activities>
  </EnterStandbyResult>
  <ResponseMetadata>
    <RequestId>7c6e177f-f082-11e1-ac58-3714bEXAMPLE</RequestId>
  </ResponseMetadata>
</EnterStandbyResponse>"""

EXIT_STANDBY_TEMPLATE = """<ExitStandbyResponse xmlns="http://autoscaling.amazonaws.com/doc/2011-01-01/">
  <ExitStandbyResult>
    <Activities>
      {% for instance in standby_instances %}
      <member>
        <ActivityId>12345678-1234-1234-1234-123456789012</ActivityId>
        <AutoScalingGroupName>{{ group_name }}</AutoScalingGroupName>
        <Description>Moving EC2 instance out of Standby: {{ instance.instance.id }}</Description>
        <Progress>30</Progress>
        <Cause>At {{ timestamp }} instance {{ instance.instance.id }} was moved out of standby in response to a user request, increasing the capacity from {{ original_size }} to {{ desired_capacity }}.</Cause>
        <StartTime>{{ timestamp }}</StartTime>
        <Details>{&quot;Subnet ID&quot;:&quot;??&quot;,&quot;Availability Zone&quot;:&quot;{{ instance.instance.placement }}&quot;}</Details>
        <StatusCode>PreInService</StatusCode>
      </member>
      {% endfor %}
    </Activities>
  </ExitStandbyResult>
  <ResponseMetadata>
    <RequestId>7c6e177f-f082-11e1-ac58-3714bEXAMPLE</RequestId>
  </ResponseMetadata>
</ExitStandbyResponse>"""

TERMINATE_INSTANCES_TEMPLATE = """<TerminateInstanceInAutoScalingGroupResponse xmlns="http://autoscaling.amazonaws.com/doc/2011-01-01/">
  <TerminateInstanceInAutoScalingGroupResult>
    <Activity>
      <ActivityId>35b5c464-0b63-2fc7-1611-467d4a7f2497EXAMPLE</ActivityId>
      <AutoScalingGroupName>{{ group_name }}</AutoScalingGroupName>
      {% if should_decrement %}
      <Cause>At {{ timestamp }} instance {{ instance.instance.id }} was taken out of service in response to a user request, shrinking the capacity from {{ original_size }} to {{ desired_capacity }}.</Cause>
      {% else %}
      <Cause>At {{ timestamp }} instance {{ instance.instance.id }} was taken out of service in response to a user request.</Cause>
      {% endif %}
      <Description>Terminating EC2 instance: {{ instance.instance.id }}</Description>
      <Progress>0</Progress>
      <StartTime>{{ timestamp }}</StartTime>
      <Details>{&quot;Subnet ID&quot;:&quot;??&quot;,&quot;Availability Zone&quot;:&quot;{{ instance.instance.placement }}&quot;}</Details>
      <StatusCode>InProgress</StatusCode>
    </Activity>
  </TerminateInstanceInAutoScalingGroupResult>
  <ResponseMetadata>
    <RequestId>a1ba8fb9-31d6-4d9a-ace1-a7f76749df11EXAMPLE</RequestId>
  </ResponseMetadata>
</TerminateInstanceInAutoScalingGroupResponse>"""

DESCRIBE_TAGS_TEMPLATE = """<DescribeTagsResponse xmlns="http://autoscaling.amazonaws.com/doc/2011-01-01/">
  <ResponseMetadata>
    <RequestId>1549581b-12b7-11e3-895e-1334aEXAMPLE</RequestId>
  </ResponseMetadata>
  <DescribeTagsResult>
    <Tags>
{% for tag in tags %}
      <member>
        <ResourceId>{{ tag.resource_id or tag.ResourceId }}</ResourceId>
        <ResourceType>{{ tag.resource_type or tag.ResourceType }}</ResourceType>
        <Key>{{ tag.key or tag.Key }}</Key>
        <Value>{{ tag.value or tag.Value }}</Value>
        <PropagateAtLaunch>{{ tag.propagate_at_launch or tag.PropagateAtLaunch }}</PropagateAtLaunch>
      </member>
{% endfor %}
    </Tags>
    {% if next_token %}
    <NextToken>{{ next_token }}</NextToken>
    {% endif %}
  </DescribeTagsResult>
</DescribeTagsResponse>"""


ENABLE_METRICS_COLLECTION_TEMPLATE = """<EnableMetricsCollectionResponse xmlns="http://autoscaling.amazonaws.com/doc/2011-01-01/">
<ResponseMetadata>
   <RequestId></RequestId>
</ResponseMetadata>
</EnableMetricsCollectionResponse>"""


PUT_WARM_POOL_TEMPLATE = """<PutWarmPoolResponse xmlns="http://autoscaling.amazonaws.com/doc/2011-01-01/">
<ResponseMetadata>
   <RequestId></RequestId>
</ResponseMetadata>
<PutWarmPoolResult></PutWarmPoolResult>
</PutWarmPoolResponse>"""


DESCRIBE_WARM_POOL_TEMPLATE = """<DescribeWarmPoolResponse xmlns="http://autoscaling.amazonaws.com/doc/2011-01-01/">
<ResponseMetadata>
   <RequestId></RequestId>
</ResponseMetadata>
<DescribeWarmPoolResult>
  {% if pool %}
  <WarmPoolConfiguration>
    {% if pool.max_capacity %}
    <MaxGroupPreparedCapacity>{{ pool.max_capacity }}</MaxGroupPreparedCapacity>
    {% endif %}
    <MinSize>{{ pool.min_size }}</MinSize>
    {% if pool.pool_state %}
    <PoolState>{{ pool.pool_state }}</PoolState>
    {% endif %}
    {% if pool.instance_reuse_policy %}
    <InstanceReusePolicy>
      <ReuseOnScaleIn>{{ 'true' if pool.instance_reuse_policy["ReuseOnScaleIn"] else 'false' }}</ReuseOnScaleIn>
    </InstanceReusePolicy>
    {% endif %}
  </WarmPoolConfiguration>
  {% endif %}
  <Instances>
  </Instances>
</DescribeWarmPoolResult>
</DescribeWarmPoolResponse>"""


DELETE_WARM_POOL_TEMPLATE = """<DeleteWarmPoolResponse xmlns="http://autoscaling.amazonaws.com/doc/2011-01-01/">
<ResponseMetadata>
   <RequestId></RequestId>
</ResponseMetadata>
<DeleteWarmPoolResult></DeleteWarmPoolResult>
</DeleteWarmPoolResponse>"""
