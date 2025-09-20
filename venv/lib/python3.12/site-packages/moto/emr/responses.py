import re
from typing import Any, List, Pattern

from moto.core.responses import ActionResult, BaseResponse, EmptyResult

from .exceptions import ValidationException
from .models import ElasticMapReduceBackend, emr_backends
from .utils import ReleaseLabel


class ElasticMapReduceResponse(BaseResponse):
    # EMR end points are inconsistent in the placement of region name
    # in the URL, so parsing it out needs to be handled differently
    emr_region_regex: List[Pattern[str]] = [
        re.compile(r"elasticmapreduce\.(.+?)\.amazonaws\.com"),
        re.compile(r"(.+?)\.elasticmapreduce\.amazonaws\.com"),
    ]

    RESPONSE_KEY_PATH_TO_TRANSFORMER = {
        "Properties": lambda x: x.original_dict() if hasattr(x, "original_dict") else x,
    }

    def __init__(self) -> None:
        super().__init__(service_name="emr")
        self.automated_parameter_parsing = True

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
        instance_groups = self._get_param("InstanceGroups", [])
        fake_groups = self.backend.add_instance_groups(jobflow_id, instance_groups)
        result = {"InstanceGroups": fake_groups}
        return ActionResult(result)

    def add_job_flow_steps(self) -> ActionResult:
        job_flow_id = self._get_param("JobFlowId")
        steps = self.backend.add_job_flow_steps(
            job_flow_id, self._get_param("Steps", [])
        )
        result = {"StepIds": [step.id for step in steps]}
        return ActionResult(result)

    def add_tags(self) -> ActionResult:
        cluster_id = self._get_param("ResourceId")
        tags = self._get_param("Tags", [])
        tags = {d["Key"]: d["Value"] for d in tags}
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
        job_flow_ids = self._get_param("JobFlowIds", [])
        job_flow_states = self._get_param("JobFlowStates", [])
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
        cluster_states = self._get_param("ClusterStates", [])
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
        step_ids = self._get_param("StepIds", [])
        step_states = self._get_param("StepStates", [])
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
        instance_groups = self._get_param("InstanceGroups", [])
        self.backend.modify_instance_groups(instance_groups)
        return EmptyResult()

    def remove_tags(self) -> ActionResult:
        cluster_id = self._get_param("ResourceId")
        tag_keys = self._get_param("TagKeys", [])
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
            additional_master_security_groups=self._get_param(
                "Instances.AdditionalMasterSecurityGroups", []
            ),
            additional_slave_security_groups=self._get_param(
                "Instances.AdditionalSlaveSecurityGroups", []
            ),
        )

        kwargs = dict(
            name=self._get_param("Name"),
            log_uri=self._get_param("LogUri"),
            job_flow_role=self._get_param("JobFlowRole"),
            service_role=self._get_param("ServiceRole"),
            auto_scaling_role=self._get_param("AutoScalingRole"),
            steps=self._get_param("Steps", []),
            ebs_root_volume_iops=self._get_param("EbsRootVolumeIops"),
            ebs_root_volume_size=self._get_param("EbsRootVolumeSize"),
            ebs_root_volume_throughput=self._get_param("EbsRootVolumeThroughput"),
            visible_to_all_users=self._get_bool_param("VisibleToAllUsers", False),
            instance_attrs=instance_attrs,
        )

        bootstrap_actions = self._get_param("BootstrapActions", [])
        if bootstrap_actions:
            kwargs["bootstrap_actions"] = bootstrap_actions

        configurations = self._get_param("Configurations", [])
        if configurations:
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
                raise ValidationException(message)
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
                raise ValidationException(message)
            elif ami_version:
                message = "Custom AMI is not supported in this version of EMR"
                raise ValidationException(message)

        step_concurrency_level = self._get_param("StepConcurrencyLevel")
        if step_concurrency_level:
            kwargs["step_concurrency_level"] = step_concurrency_level

        security_configuration = self._get_param("SecurityConfiguration")
        if security_configuration:
            kwargs["security_configuration"] = security_configuration

        kerberos_attributes = self._get_param("KerberosAttributes", {})
        kwargs["kerberos_attributes"] = kerberos_attributes

        cluster = self.backend.run_job_flow(**kwargs)

        applications = self._get_param("Applications", [])
        if applications:
            self.backend.add_applications(cluster.id, applications)
        else:
            self.backend.add_applications(
                cluster.id, [{"Name": "Hadoop", "Version": "0.18"}]
            )

        instance_groups = self._get_param("Instances.InstanceGroups", [])
        if instance_groups:
            instance_group_result = self.backend.add_instance_groups(
                cluster.id, instance_groups
            )
            for i in range(0, len(instance_group_result)):
                self.backend.run_instances(
                    cluster.id, instance_groups[i], instance_group_result[i]
                )

        # TODO: Instances also must be created when `Instances.InstanceType` and `Instances.InstanceCount` are specified in the request.

        tags = self._get_param("Tags", [])
        if tags:
            self.backend.add_tags(
                cluster.id, dict((d["key"], d["value"]) for d in tags)
            )
        result = {
            "JobFlowId": cluster.job_flow_id,
            "ClusterArn": cluster.arn,
        }
        return ActionResult(result)

    def set_termination_protection(self) -> ActionResult:
        termination_protection = self._get_bool_param("TerminationProtected")
        job_ids = self._get_param("JobFlowIds", [])
        self.backend.set_termination_protection(job_ids, termination_protection)
        return EmptyResult()

    def set_visible_to_all_users(self) -> ActionResult:
        visible_to_all_users = self._get_bool_param("VisibleToAllUsers", False)
        job_ids = self._get_param("JobFlowIds", [])
        self.backend.set_visible_to_all_users(job_ids, visible_to_all_users)
        return EmptyResult()

    def terminate_job_flows(self) -> ActionResult:
        job_ids = self._get_param("JobFlowIds")
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

    def list_release_labels(self) -> ActionResult:
        release_labels = self.backend.list_release_labels()
        return ActionResult({"ReleaseLabels": release_labels})

    def list_supported_instance_types(self) -> ActionResult:
        release_label = self._get_param("ReleaseLabel")
        instance_types = self.backend.list_supported_instance_types(release_label)
        return ActionResult({"SupportedInstanceTypes": instance_types})
