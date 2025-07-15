import re
from typing import Any, Dict, List, Pattern

from moto.core.responses import ActionResult, AWSServiceSpec, BaseResponse, EmptyResult
from moto.core.utils import tags_from_query_string

from .exceptions import ValidationException
from .models import ElasticMapReduceBackend, emr_backends
from .utils import ReleaseLabel, Unflattener, steps_from_query_string


class ElasticMapReduceResponse(BaseResponse):
    # EMR end points are inconsistent in the placement of region name
    # in the URL, so parsing it out needs to be handled differently
    emr_region_regex: List[Pattern[str]] = [
        re.compile(r"elasticmapreduce\.(.+?)\.amazonaws\.com"),
        re.compile(r"(.+?)\.elasticmapreduce\.amazonaws\.com"),
    ]

    aws_service_spec = AWSServiceSpec("data/emr/2009-03-31/service-2.json")

    def __init__(self) -> None:
        super().__init__(service_name="emr")

    def get_region_from_url(self, request: Any, full_url: str) -> str:
        for regex in ElasticMapReduceResponse.emr_region_regex:
            match = regex.search(self.parsed_url.netloc)
            if match:
                return match.group(1)
        return self.default_region

    @property
    def backend(self) -> ElasticMapReduceBackend:
        return emr_backends[self.current_account][self.region]

    def add_instance_groups(self) -> ActionResult:
        jobflow_id = self._get_param("JobFlowId")
        instance_groups = self._get_list_prefix("InstanceGroups.member")
        for item in instance_groups:
            item["instance_count"] = int(item["instance_count"])
            # Adding support to EbsConfiguration
            self._parse_ebs_configuration(item)
            # Adding support for auto_scaling_policy
            Unflattener.unflatten_complex_params(item, "auto_scaling_policy")
        fake_groups = self.backend.add_instance_groups(jobflow_id, instance_groups)
        result = {"InstanceGroups": fake_groups}
        return ActionResult(result)

    def add_job_flow_steps(self) -> ActionResult:
        job_flow_id = self._get_param("JobFlowId")
        steps = self.backend.add_job_flow_steps(
            job_flow_id, steps_from_query_string(self._get_list_prefix("Steps.member"))
        )
        result = {"StepIds": [step.id for step in steps]}
        return ActionResult(result)

    def add_tags(self) -> ActionResult:
        cluster_id = self._get_param("ResourceId")
        tags = tags_from_query_string(self.querystring, prefix="Tags")
        self.backend.add_tags(cluster_id, tags)
        return EmptyResult()

    def create_security_configuration(self) -> ActionResult:
        name = self._get_param("Name")
        security_configuration = self._get_param("SecurityConfiguration")
        security_configuration = self.backend.create_security_configuration(
            name=name, security_configuration=security_configuration
        )
        return ActionResult(security_configuration)

    def describe_security_configuration(self) -> ActionResult:
        name = self._get_param("Name")
        security_configuration = self.backend.get_security_configuration(name=name)
        return ActionResult(security_configuration)

    def delete_security_configuration(self) -> ActionResult:
        name = self._get_param("Name")
        self.backend.delete_security_configuration(name=name)
        return EmptyResult()

    def describe_cluster(self) -> ActionResult:
        cluster_id = self._get_param("ClusterId")
        cluster = self.backend.describe_cluster(cluster_id)
        result = {"Cluster": cluster}
        return ActionResult(result)

    def describe_job_flows(self) -> ActionResult:
        created_after = self._get_param("CreatedAfter")
        created_before = self._get_param("CreatedBefore")
        job_flow_ids = self._get_multi_param("JobFlowIds.member")
        job_flow_states = self._get_multi_param("JobFlowStates.member")
        clusters = self.backend.describe_job_flows(
            job_flow_ids, job_flow_states, created_after, created_before
        )
        result = {"JobFlows": clusters}
        return ActionResult(result)

    def describe_step(self) -> ActionResult:
        cluster_id = self._get_param("ClusterId")
        step_id = self._get_param("StepId")
        step = self.backend.describe_step(cluster_id, step_id)
        result = {"Step": step}
        return ActionResult(result)

    def list_bootstrap_actions(self) -> ActionResult:
        cluster_id = self._get_param("ClusterId")
        marker = self._get_param("Marker")
        bootstrap_actions, marker = self.backend.list_bootstrap_actions(
            cluster_id, marker
        )
        result = {"BootstrapActions": bootstrap_actions, "Marker": marker}
        return ActionResult(result)

    def list_clusters(self) -> ActionResult:
        cluster_states = self._get_multi_param("ClusterStates.member")
        created_after = self._get_param("CreatedAfter")
        created_before = self._get_param("CreatedBefore")
        marker = self._get_param("Marker")
        clusters, marker = self.backend.list_clusters(
            cluster_states, created_after, created_before, marker
        )
        result = {"Clusters": clusters, "Marker": marker}
        return ActionResult(result)

    def list_instance_groups(self) -> ActionResult:
        cluster_id = self._get_param("ClusterId")
        marker = self._get_param("Marker")
        instance_groups, marker = self.backend.list_instance_groups(
            cluster_id, marker=marker
        )
        result = {"InstanceGroups": instance_groups, "Marker": marker}
        return ActionResult(result)

    def list_instances(self) -> ActionResult:
        cluster_id = self._get_param("ClusterId")
        marker = self._get_param("Marker")
        instance_group_id = self._get_param("InstanceGroupId")
        instance_group_types = self._get_param("InstanceGroupTypes")
        instances, marker = self.backend.list_instances(
            cluster_id,
            marker=marker,
            instance_group_id=instance_group_id,
            instance_group_types=instance_group_types,
        )
        result = {"Instances": instances, "Marker": marker}
        return ActionResult(result)

    def list_steps(self) -> ActionResult:
        cluster_id = self._get_param("ClusterId")
        marker = self._get_param("Marker")
        step_ids = self._get_multi_param("StepIds.member")
        step_states = self._get_multi_param("StepStates.member")
        steps, marker = self.backend.list_steps(
            cluster_id, marker=marker, step_ids=step_ids, step_states=step_states
        )
        result = {"Steps": steps, "Marker": marker}
        return ActionResult(result)

    def modify_cluster(self) -> ActionResult:
        cluster_id = self._get_param("ClusterId")
        step_concurrency_level = self._get_param("StepConcurrencyLevel")
        cluster = self.backend.modify_cluster(cluster_id, step_concurrency_level)
        result = {"StepConcurrencyLevel": cluster.step_concurrency_level}
        return ActionResult(result)

    def modify_instance_groups(self) -> ActionResult:
        instance_groups = self._get_list_prefix("InstanceGroups.member")
        for item in instance_groups:
            item["instance_count"] = int(item["instance_count"])
        self.backend.modify_instance_groups(instance_groups)
        return EmptyResult()

    def remove_tags(self) -> ActionResult:
        cluster_id = self._get_param("ResourceId")
        tag_keys = self._get_multi_param("TagKeys.member")
        self.backend.remove_tags(cluster_id, tag_keys)
        return EmptyResult()

    def run_job_flow(self) -> ActionResult:
        instance_attrs = dict(
            master_instance_type=self._get_param("Instances.MasterInstanceType"),
            slave_instance_type=self._get_param("Instances.SlaveInstanceType"),
            instance_count=self._get_int_param("Instances.InstanceCount", 1),
            ec2_key_name=self._get_param("Instances.Ec2KeyName"),
            ec2_subnet_id=self._get_param("Instances.Ec2SubnetId"),
            hadoop_version=self._get_param("Instances.HadoopVersion"),
            availability_zone=self._get_param(
                "Instances.Placement.AvailabilityZone", self.backend.region_name + "a"
            ),
            keep_job_flow_alive_when_no_steps=self._get_bool_param(
                "Instances.KeepJobFlowAliveWhenNoSteps", False
            ),
            termination_protected=self._get_bool_param(
                "Instances.TerminationProtected", False
            ),
            emr_managed_master_security_group=self._get_param(
                "Instances.EmrManagedMasterSecurityGroup"
            ),
            emr_managed_slave_security_group=self._get_param(
                "Instances.EmrManagedSlaveSecurityGroup"
            ),
            service_access_security_group=self._get_param(
                "Instances.ServiceAccessSecurityGroup"
            ),
            additional_master_security_groups=self._get_multi_param(
                "Instances.AdditionalMasterSecurityGroups.member."
            ),
            additional_slave_security_groups=self._get_multi_param(
                "Instances.AdditionalSlaveSecurityGroups.member."
            ),
        )

        kwargs = dict(
            name=self._get_param("Name"),
            log_uri=self._get_param("LogUri"),
            job_flow_role=self._get_param("JobFlowRole"),
            service_role=self._get_param("ServiceRole"),
            auto_scaling_role=self._get_param("AutoScalingRole"),
            steps=steps_from_query_string(self._get_list_prefix("Steps.member")),
            visible_to_all_users=self._get_bool_param("VisibleToAllUsers", False),
            instance_attrs=instance_attrs,
        )

        bootstrap_actions = self._get_list_prefix("BootstrapActions.member")
        if bootstrap_actions:
            for ba in bootstrap_actions:
                args = []
                idx = 1
                keyfmt = "script_bootstrap_action._args.member.{0}"
                key = keyfmt.format(idx)
                while key in ba:
                    args.append(ba.pop(key))
                    idx += 1
                    key = keyfmt.format(idx)
                ba["args"] = args
                ba["script_path"] = ba.pop("script_bootstrap_action._path")
            kwargs["bootstrap_actions"] = bootstrap_actions

        configurations = self._get_list_prefix("Configurations.member")
        if configurations:
            for idx, config in enumerate(configurations, 1):
                for key in list(config.keys()):
                    if key.startswith("properties."):
                        config.pop(key)
                config["properties"] = {}
                map_items = self._get_map_prefix(
                    f"Configurations.member.{idx}.Properties.entry"
                )
                config["properties"] = map_items

            kwargs["configurations"] = configurations

        release_label = self._get_param("ReleaseLabel")
        ami_version = self._get_param("AmiVersion")
        if release_label:
            kwargs["release_label"] = release_label
            if ami_version:
                message = (
                    "Only one AMI version and release label may be specified. "
                    "Provided AMI: {0}, release label: {1}."
                ).format(ami_version, release_label)
                raise ValidationException(message=message)
        else:
            if ami_version:
                kwargs["requested_ami_version"] = ami_version
                kwargs["running_ami_version"] = ami_version
            else:
                kwargs["running_ami_version"] = "1.0.0"

        custom_ami_id = self._get_param("CustomAmiId")
        if custom_ami_id:
            kwargs["custom_ami_id"] = custom_ami_id
            if release_label and (
                ReleaseLabel(release_label) < ReleaseLabel("emr-5.7.0")
            ):
                message = "Custom AMI is not allowed"
                raise ValidationException(message=message)
            elif ami_version:
                message = "Custom AMI is not supported in this version of EMR"
                raise ValidationException(message=message)

        step_concurrency_level = self._get_param("StepConcurrencyLevel")
        if step_concurrency_level:
            kwargs["step_concurrency_level"] = step_concurrency_level

        security_configuration = self._get_param("SecurityConfiguration")
        if security_configuration:
            kwargs["security_configuration"] = security_configuration

        kerberos_attributes: Dict[str, Any] = {}
        kwargs["kerberos_attributes"] = kerberos_attributes

        realm = self._get_param("KerberosAttributes.Realm")
        if realm:
            kerberos_attributes["Realm"] = realm

        kdc_admin_password = self._get_param("KerberosAttributes.KdcAdminPassword")
        if kdc_admin_password:
            kerberos_attributes["KdcAdminPassword"] = kdc_admin_password

        cross_realm_principal_password = self._get_param(
            "KerberosAttributes.CrossRealmTrustPrincipalPassword"
        )
        if cross_realm_principal_password:
            kerberos_attributes["CrossRealmTrustPrincipalPassword"] = (
                cross_realm_principal_password
            )

        ad_domain_join_user = self._get_param("KerberosAttributes.ADDomainJoinUser")
        if ad_domain_join_user:
            kerberos_attributes["ADDomainJoinUser"] = ad_domain_join_user

        ad_domain_join_password = self._get_param(
            "KerberosAttributes.ADDomainJoinPassword"
        )
        if ad_domain_join_password:
            kerberos_attributes["ADDomainJoinPassword"] = ad_domain_join_password

        cluster = self.backend.run_job_flow(**kwargs)

        applications = self._get_list_prefix("Applications.member")
        if applications:
            self.backend.add_applications(cluster.id, applications)
        else:
            self.backend.add_applications(
                cluster.id, [{"Name": "Hadoop", "Version": "0.18"}]
            )

        instance_groups = self._get_list_prefix("Instances.InstanceGroups.member")
        if instance_groups:
            for ig in instance_groups:
                ig["instance_count"] = int(ig["instance_count"])
                # Adding support to EbsConfiguration
                self._parse_ebs_configuration(ig)
                # Adding support for auto_scaling_policy
                Unflattener.unflatten_complex_params(ig, "auto_scaling_policy")
            instance_group_result = self.backend.add_instance_groups(
                cluster.id, instance_groups
            )
            for i in range(0, len(instance_group_result)):
                self.backend.run_instances(
                    cluster.id, instance_groups[i], instance_group_result[i]
                )

        tags = self._get_list_prefix("Tags.member")
        if tags:
            self.backend.add_tags(
                cluster.id, dict((d["key"], d["value"]) for d in tags)
            )
        result = {
            "JobFlowId": cluster.job_flow_id,
            "ClusterArn": cluster.arn,
        }
        return ActionResult(result)

    def _has_key_prefix(self, key_prefix: str, value: Dict[str, Any]) -> bool:
        for key in value:  # iter on both keys and values
            if key.startswith(key_prefix):
                return True
        return False

    def _parse_ebs_configuration(self, instance_group: Dict[str, Any]) -> None:
        key_ebs_config = "ebs_configuration"
        ebs_configuration = dict()
        # Filter only EBS config keys
        for key in instance_group:
            if key.startswith(key_ebs_config):
                ebs_configuration[key] = instance_group[key]

        if len(ebs_configuration) > 0:
            # Key that should be extracted
            ebs_optimized = "ebs_optimized"
            ebs_block_device_configs = "ebs_block_device_configs"
            volume_specification = "volume_specification"
            size_in_gb = "size_in_gb"
            volume_type = "volume_type"
            iops = "iops"
            volumes_per_instance = "volumes_per_instance"

            key_ebs_optimized = f"{key_ebs_config}._{ebs_optimized}"
            # EbsOptimized config
            if key_ebs_optimized in ebs_configuration:
                instance_group.pop(key_ebs_optimized)
                ebs_configuration[ebs_optimized] = ebs_configuration.pop(
                    key_ebs_optimized
                )

            # Ebs Blocks
            ebs_blocks = []
            idx = 1
            keyfmt = f"{key_ebs_config}._{ebs_block_device_configs}.member.{{}}"
            key = keyfmt.format(idx)
            while self._has_key_prefix(key, ebs_configuration):
                vlespc_keyfmt = f"{key}._{volume_specification}._{{}}"
                vol_size = vlespc_keyfmt.format(size_in_gb)
                vol_iops = vlespc_keyfmt.format(iops)
                vol_type = vlespc_keyfmt.format(volume_type)

                ebs_block: Dict[str, Any] = dict()
                ebs_block[volume_specification] = dict()
                if vol_size in ebs_configuration:
                    instance_group.pop(vol_size)
                    ebs_block[volume_specification][size_in_gb] = int(
                        ebs_configuration.pop(vol_size)
                    )
                if vol_iops in ebs_configuration:
                    instance_group.pop(vol_iops)
                    ebs_block[volume_specification][iops] = ebs_configuration.pop(
                        vol_iops
                    )
                if vol_type in ebs_configuration:
                    instance_group.pop(vol_type)
                    ebs_block[volume_specification][volume_type] = (
                        ebs_configuration.pop(vol_type)
                    )

                per_instance = f"{key}._{volumes_per_instance}"
                if per_instance in ebs_configuration:
                    instance_group.pop(per_instance)
                    ebs_block[volumes_per_instance] = int(
                        ebs_configuration.pop(per_instance)
                    )

                if len(ebs_block) > 0:
                    ebs_blocks.append(ebs_block)
                idx += 1
                key = keyfmt.format(idx)

            if len(ebs_blocks) > 0:
                ebs_configuration[ebs_block_device_configs] = ebs_blocks
            instance_group[key_ebs_config] = ebs_configuration

    def set_termination_protection(self) -> ActionResult:
        termination_protection = self._get_bool_param("TerminationProtected")
        job_ids = self._get_multi_param("JobFlowIds.member")
        self.backend.set_termination_protection(job_ids, termination_protection)
        return EmptyResult()

    def set_visible_to_all_users(self) -> ActionResult:
        visible_to_all_users = self._get_bool_param("VisibleToAllUsers", False)
        job_ids = self._get_multi_param("JobFlowIds.member")
        self.backend.set_visible_to_all_users(job_ids, visible_to_all_users)
        return EmptyResult()

    def terminate_job_flows(self) -> ActionResult:
        job_ids = self._get_multi_param("JobFlowIds.member.")
        self.backend.terminate_job_flows(job_ids)
        return EmptyResult()

    def put_auto_scaling_policy(self) -> ActionResult:
        cluster_id = self._get_param("ClusterId")
        cluster = self.backend.describe_cluster(cluster_id)
        instance_group_id = self._get_param("InstanceGroupId")
        auto_scaling_policy = self._get_param("AutoScalingPolicy")
        instance_group = self.backend.put_auto_scaling_policy(
            instance_group_id, auto_scaling_policy
        )
        assert instance_group is not None
        result = {
            "ClusterId": cluster.id,
            "InstanceGroupId": instance_group.id,
            "AutoScalingPolicy": instance_group.auto_scaling_policy,
            "ClusterArn": cluster.arn,
        }
        return ActionResult(result)

    def remove_auto_scaling_policy(self) -> ActionResult:
        instance_group_id = self._get_param("InstanceGroupId")
        self.backend.remove_auto_scaling_policy(instance_group_id)
        return EmptyResult()

    def get_block_public_access_configuration(self) -> ActionResult:
        configuration = self.backend.get_block_public_access_configuration()
        config = configuration.get("block_public_access_configuration") or {}
        metadata = configuration.get("block_public_access_configuration_metadata") or {}
        result = {
            "BlockPublicAccessConfiguration": config,
            "BlockPublicAccessConfigurationMetadata": metadata,
        }
        return ActionResult(result)

    def put_block_public_access_configuration(self) -> ActionResult:
        params = self._get_params()
        block_public_access_configuration = (
            params.get("BlockPublicAccessConfiguration") or {}
        )
        self.backend.put_block_public_access_configuration(
            block_public_security_group_rules=block_public_access_configuration.get(
                "BlockPublicSecurityGroupRules"
            )
            or True,
            rule_ranges=block_public_access_configuration.get(
                "PermittedPublicSecurityGroupRuleRanges"
            ),
        )
        return EmptyResult()
