import json
import re
from datetime import datetime, timezone
from functools import wraps
from typing import Any, Callable, Dict, List, Pattern

from moto.core.responses import AWSServiceSpec, BaseResponse, xml_to_json_response
from moto.core.utils import tags_from_query_string

from .exceptions import ValidationException
from .models import ElasticMapReduceBackend, emr_backends
from .utils import ReleaseLabel, Unflattener, steps_from_query_string


def generate_boto3_response(
    operation: str,
) -> Callable[
    [Callable[["ElasticMapReduceResponse"], str]],
    Callable[["ElasticMapReduceResponse"], str],
]:
    """The decorator to convert an XML response to JSON, if the request is
    determined to be from boto3. Pass the API action as a parameter.

    """

    def _boto3_request(
        method: Callable[["ElasticMapReduceResponse"], str],
    ) -> Callable[["ElasticMapReduceResponse"], str]:
        @wraps(method)
        def f(self: "ElasticMapReduceResponse") -> str:
            rendered = method(self)
            if "json" in self.headers.get("Content-Type", []):
                self.response_headers.update(
                    {
                        "x-amzn-requestid": "2690d7eb-ed86-11dd-9877-6fad448a8419",
                        "date": datetime.now(timezone.utc).strftime(
                            "%a, %d %b %Y %H:%M:%S %Z"
                        ),
                        "content-type": "application/x-amz-json-1.1",
                    }
                )
                resp = xml_to_json_response(self.aws_service_spec, operation, rendered)
                return "" if resp is None else json.dumps(resp)
            return rendered

        return f

    return _boto3_request


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

    @generate_boto3_response("AddInstanceGroups")
    def add_instance_groups(self) -> str:
        jobflow_id = self._get_param("JobFlowId")
        instance_groups = self._get_list_prefix("InstanceGroups.member")
        for item in instance_groups:
            item["instance_count"] = int(item["instance_count"])
            # Adding support to EbsConfiguration
            self._parse_ebs_configuration(item)
            # Adding support for auto_scaling_policy
            Unflattener.unflatten_complex_params(item, "auto_scaling_policy")
        fake_groups = self.backend.add_instance_groups(jobflow_id, instance_groups)
        template = self.response_template(ADD_INSTANCE_GROUPS_TEMPLATE)
        return template.render(instance_groups=fake_groups)

    @generate_boto3_response("AddJobFlowSteps")
    def add_job_flow_steps(self) -> str:
        job_flow_id = self._get_param("JobFlowId")
        steps = self.backend.add_job_flow_steps(
            job_flow_id, steps_from_query_string(self._get_list_prefix("Steps.member"))
        )
        template = self.response_template(ADD_JOB_FLOW_STEPS_TEMPLATE)
        return template.render(steps=steps)

    @generate_boto3_response("AddTags")
    def add_tags(self) -> str:
        cluster_id = self._get_param("ResourceId")
        tags = tags_from_query_string(self.querystring, prefix="Tags")
        self.backend.add_tags(cluster_id, tags)
        template = self.response_template(ADD_TAGS_TEMPLATE)
        return template.render()

    @generate_boto3_response("CreateSecurityConfiguration")
    def create_security_configuration(self) -> str:
        name = self._get_param("Name")
        security_configuration = self._get_param("SecurityConfiguration")
        resp = self.backend.create_security_configuration(
            name=name, security_configuration=security_configuration
        )
        template = self.response_template(CREATE_SECURITY_CONFIGURATION_TEMPLATE)
        return template.render(name=name, creation_date_time=resp.creation_date_time)

    @generate_boto3_response("DescribeSecurityConfiguration")
    def describe_security_configuration(self) -> str:
        name = self._get_param("Name")
        security_configuration = self.backend.get_security_configuration(name=name)
        template = self.response_template(DESCRIBE_SECURITY_CONFIGURATION_TEMPLATE)
        return template.render(security_configuration=security_configuration)

    @generate_boto3_response("DeleteSecurityConfiguration")
    def delete_security_configuration(self) -> str:
        name = self._get_param("Name")
        self.backend.delete_security_configuration(name=name)
        template = self.response_template(DELETE_SECURITY_CONFIGURATION_TEMPLATE)
        return template.render()

    @generate_boto3_response("DescribeCluster")
    def describe_cluster(self) -> str:
        cluster_id = self._get_param("ClusterId")
        cluster = self.backend.describe_cluster(cluster_id)
        template = self.response_template(DESCRIBE_CLUSTER_TEMPLATE)
        return template.render(cluster=cluster)

    @generate_boto3_response("DescribeJobFlows")
    def describe_job_flows(self) -> str:
        created_after = self._get_param("CreatedAfter")
        created_before = self._get_param("CreatedBefore")
        job_flow_ids = self._get_multi_param("JobFlowIds.member")
        job_flow_states = self._get_multi_param("JobFlowStates.member")
        clusters = self.backend.describe_job_flows(
            job_flow_ids, job_flow_states, created_after, created_before
        )
        template = self.response_template(DESCRIBE_JOB_FLOWS_TEMPLATE)
        return template.render(clusters=clusters)

    @generate_boto3_response("DescribeStep")
    def describe_step(self) -> str:
        cluster_id = self._get_param("ClusterId")
        step_id = self._get_param("StepId")
        step = self.backend.describe_step(cluster_id, step_id)
        template = self.response_template(DESCRIBE_STEP_TEMPLATE)
        return template.render(step=step)

    @generate_boto3_response("ListBootstrapActions")
    def list_bootstrap_actions(self) -> str:
        cluster_id = self._get_param("ClusterId")
        marker = self._get_param("Marker")
        bootstrap_actions, marker = self.backend.list_bootstrap_actions(
            cluster_id, marker
        )
        template = self.response_template(LIST_BOOTSTRAP_ACTIONS_TEMPLATE)
        return template.render(bootstrap_actions=bootstrap_actions, marker=marker)

    @generate_boto3_response("ListClusters")
    def list_clusters(self) -> str:
        cluster_states = self._get_multi_param("ClusterStates.member")
        created_after = self._get_param("CreatedAfter")
        created_before = self._get_param("CreatedBefore")
        marker = self._get_param("Marker")
        clusters, marker = self.backend.list_clusters(
            cluster_states, created_after, created_before, marker
        )
        template = self.response_template(LIST_CLUSTERS_TEMPLATE)
        return template.render(clusters=clusters, marker=marker)

    @generate_boto3_response("ListInstanceGroups")
    def list_instance_groups(self) -> str:
        cluster_id = self._get_param("ClusterId")
        marker = self._get_param("Marker")
        instance_groups, marker = self.backend.list_instance_groups(
            cluster_id, marker=marker
        )
        template = self.response_template(LIST_INSTANCE_GROUPS_TEMPLATE)
        return template.render(instance_groups=instance_groups, marker=marker)

    @generate_boto3_response("ListInstances")
    def list_instances(self) -> str:
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
        template = self.response_template(LIST_INSTANCES_TEMPLATE)
        return template.render(instances=instances, marker=marker)

    @generate_boto3_response("ListSteps")
    def list_steps(self) -> str:
        cluster_id = self._get_param("ClusterId")
        marker = self._get_param("Marker")
        step_ids = self._get_multi_param("StepIds.member")
        step_states = self._get_multi_param("StepStates.member")
        steps, marker = self.backend.list_steps(
            cluster_id, marker=marker, step_ids=step_ids, step_states=step_states
        )
        template = self.response_template(LIST_STEPS_TEMPLATE)
        return template.render(steps=steps, marker=marker)

    @generate_boto3_response("ModifyCluster")
    def modify_cluster(self) -> str:
        cluster_id = self._get_param("ClusterId")
        step_concurrency_level = self._get_param("StepConcurrencyLevel")
        cluster = self.backend.modify_cluster(cluster_id, step_concurrency_level)
        template = self.response_template(MODIFY_CLUSTER_TEMPLATE)
        return template.render(cluster=cluster)

    @generate_boto3_response("ModifyInstanceGroups")
    def modify_instance_groups(self) -> str:
        instance_groups = self._get_list_prefix("InstanceGroups.member")
        for item in instance_groups:
            item["instance_count"] = int(item["instance_count"])
        self.backend.modify_instance_groups(instance_groups)
        template = self.response_template(MODIFY_INSTANCE_GROUPS_TEMPLATE)
        return template.render()

    @generate_boto3_response("RemoveTags")
    def remove_tags(self) -> str:
        cluster_id = self._get_param("ResourceId")
        tag_keys = self._get_multi_param("TagKeys.member")
        self.backend.remove_tags(cluster_id, tag_keys)
        template = self.response_template(REMOVE_TAGS_TEMPLATE)
        return template.render()

    @generate_boto3_response("RunJobFlow")
    def run_job_flow(self) -> str:
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

        template = self.response_template(RUN_JOB_FLOW_TEMPLATE)
        return template.render(cluster=cluster)

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

    @generate_boto3_response("SetTerminationProtection")
    def set_termination_protection(self) -> str:
        termination_protection = self._get_bool_param("TerminationProtected")
        job_ids = self._get_multi_param("JobFlowIds.member")
        self.backend.set_termination_protection(job_ids, termination_protection)
        template = self.response_template(SET_TERMINATION_PROTECTION_TEMPLATE)
        return template.render()

    @generate_boto3_response("SetVisibleToAllUsers")
    def set_visible_to_all_users(self) -> str:
        visible_to_all_users = self._get_param("VisibleToAllUsers")
        job_ids = self._get_multi_param("JobFlowIds.member")
        self.backend.set_visible_to_all_users(job_ids, visible_to_all_users)
        template = self.response_template(SET_VISIBLE_TO_ALL_USERS_TEMPLATE)
        return template.render()

    @generate_boto3_response("TerminateJobFlows")
    def terminate_job_flows(self) -> str:
        job_ids = self._get_multi_param("JobFlowIds.member.")
        self.backend.terminate_job_flows(job_ids)
        template = self.response_template(TERMINATE_JOB_FLOWS_TEMPLATE)
        return template.render()

    @generate_boto3_response("PutAutoScalingPolicy")
    def put_auto_scaling_policy(self) -> str:
        cluster_id = self._get_param("ClusterId")
        cluster = self.backend.describe_cluster(cluster_id)
        instance_group_id = self._get_param("InstanceGroupId")
        auto_scaling_policy = self._get_param("AutoScalingPolicy")
        instance_group = self.backend.put_auto_scaling_policy(
            instance_group_id, auto_scaling_policy
        )
        template = self.response_template(PUT_AUTO_SCALING_POLICY)
        return template.render(
            cluster_id=cluster_id, cluster=cluster, instance_group=instance_group
        )

    @generate_boto3_response("RemoveAutoScalingPolicy")
    def remove_auto_scaling_policy(self) -> str:
        instance_group_id = self._get_param("InstanceGroupId")
        self.backend.remove_auto_scaling_policy(instance_group_id)
        template = self.response_template(REMOVE_AUTO_SCALING_POLICY)
        return template.render()

    @generate_boto3_response("GetBlockPublicAccessConfiguration")
    def get_block_public_access_configuration(self) -> str:
        configuration = self.backend.get_block_public_access_configuration()
        config = configuration.get("block_public_access_configuration") or {}
        metadata = configuration.get("block_public_access_configuration_metadata") or {}
        template = self.response_template(
            GET_BLOCK_PUBLIC_ACCESS_CONFIGURATION_TEMPLATE
        )
        return template.render(
            block_public_security_group_rules=config.get(
                "block_public_security_group_rules"
            ),
            permitted_public_security_group_rule_ranges=config.get(
                "permitted_public_security_group_rule_ranges"
            ),
            creation_date_time=metadata.get("creation_date_time"),
            created_by_arn=metadata.get("created_by_arn"),
        )

    @generate_boto3_response("PutBlockPublicAccessConfiguration")
    def put_block_public_access_configuration(self) -> str:
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
        template = self.response_template(
            PUT_BLOCK_PUBLIC_ACCESS_CONFIGURATION_TEMPLATE
        )
        return template.render()


ADD_INSTANCE_GROUPS_TEMPLATE = """<AddInstanceGroupsResponse xmlns="http://elasticmapreduce.amazonaws.com/doc/2009-03-31">
  <AddInstanceGroupsResult>
    <InstanceGroupIds>
      {% for instance_group in instance_groups %}
      <member>{{ instance_group.id }}</member>
      {% endfor %}
    </InstanceGroupIds>
  </AddInstanceGroupsResult>
  <ResponseMetadata>
    <RequestId>2690d7eb-ed86-11dd-9877-6fad448a8419</RequestId>
  </ResponseMetadata>
</AddInstanceGroupsResponse>"""

ADD_JOB_FLOW_STEPS_TEMPLATE = """<AddJobFlowStepsResponse xmlns="http://elasticmapreduce.amazonaws.com/doc/2009-03-31">
  <AddJobFlowStepsResult>
    <StepIds>
      {% for step in steps %}
      <member>{{ step.id }}</member>
      {% endfor %}
    </StepIds>
  </AddJobFlowStepsResult>
  <ResponseMetadata>
    <RequestId>df6f4f4a-ed85-11dd-9877-6fad448a8419</RequestId>
  </ResponseMetadata>
</AddJobFlowStepsResponse>"""

ADD_TAGS_TEMPLATE = """<AddTagsResponse xmlns="http://elasticmapreduce.amazonaws.com/doc/2009-03-31">
  <ResponseMetadata>
    <RequestId>2690d7eb-ed86-11dd-9877-6fad448a8419</RequestId>
  </ResponseMetadata>
</AddTagsResponse>"""

DESCRIBE_CLUSTER_TEMPLATE = """<DescribeClusterResponse xmlns="http://elasticmapreduce.amazonaws.com/doc/2009-03-31">
  <DescribeClusterResult>
    <Cluster>
      <Applications>
        {% for application in cluster.applications %}
        <member>
          <Name>{{ application.name }}</Name>
          <Version>{{ application.version }}</Version>
        </member>
        {% endfor %}
      </Applications>
      <AutoTerminate>{{ (not cluster.keep_job_flow_alive_when_no_steps)|lower }}</AutoTerminate>
      <Configurations>
        {% for configuration in cluster.configurations %}
        <member>
          <Classification>{{ configuration['classification'] }}</Classification>
          <Properties>
            {% for key, value in configuration['properties'].items() %}
            <entry>
              <key>{{ key }}</key>
              <value>{{ value }}</value>
            </entry>
            {% endfor %}
          </Properties>
        </member>
        {% endfor %}
      </Configurations>
      {% if cluster.custom_ami_id is not none %}
      <CustomAmiId>{{ cluster.custom_ami_id }}</CustomAmiId>
      {% endif %}
      <Ec2InstanceAttributes>
        <AdditionalMasterSecurityGroups>
        {% for each in cluster.additional_master_security_groups %}
          <member>{{ each }}</member>
        {% endfor %}
        </AdditionalMasterSecurityGroups>
        <AdditionalSlaveSecurityGroups>
        {% for each in cluster.additional_slave_security_groups %}
          <member>{{ each }}</member>
        {% endfor %}
        </AdditionalSlaveSecurityGroups>
        <Ec2AvailabilityZone>{{ cluster.availability_zone }}</Ec2AvailabilityZone>
        <Ec2KeyName>{{ cluster.ec2_key_name }}</Ec2KeyName>
        <Ec2SubnetId>{{ cluster.ec2_subnet_id }}</Ec2SubnetId>
        <IamInstanceProfile>{{ cluster.role }}</IamInstanceProfile>
        <EmrManagedMasterSecurityGroup>{{ cluster.master_security_group }}</EmrManagedMasterSecurityGroup>
        <EmrManagedSlaveSecurityGroup>{{ cluster.slave_security_group }}</EmrManagedSlaveSecurityGroup>
        <ServiceAccessSecurityGroup>{{ cluster.service_access_security_group }}</ServiceAccessSecurityGroup>
      </Ec2InstanceAttributes>
      <Id>{{ cluster.id }}</Id>
      <KerberosAttributes>
        {% if 'Realm' in cluster.kerberos_attributes%}
        <Realm>{{ cluster.kerberos_attributes['Realm'] }}</Realm>
        {% endif %}
        {% if 'KdcAdminPassword' in cluster.kerberos_attributes%}
        <KdcAdminPassword>{{ cluster.kerberos_attributes['KdcAdminPassword'] }}</KdcAdminPassword>
        {% endif %}
        {% if 'CrossRealmTrustPrincipalPassword' in cluster.kerberos_attributes%}
        <CrossRealmTrustPrincipalPassword>{{ cluster.kerberos_attributes['CrossRealmTrustPrincipalPassword'] }}</CrossRealmTrustPrincipalPassword>
        {% endif %}
        {% if 'ADDomainJoinUser' in cluster.kerberos_attributes%}
        <ADDomainJoinUser>{{ cluster.kerberos_attributes['ADDomainJoinUser'] }}</ADDomainJoinUser>
        {% endif %}
        {% if 'ADDomainJoinPassword' in cluster.kerberos_attributes%}
        <ADDomainJoinPassword>{{ cluster.kerberos_attributes['ADDomainJoinPassword'] }}</ADDomainJoinPassword>
        {% endif %}
      </KerberosAttributes>
      <LogUri>{{ cluster.log_uri }}</LogUri>
      <MasterPublicDnsName>ec2-184-0-0-1.us-west-1.compute.amazonaws.com</MasterPublicDnsName>
      <Name>{{ cluster.name }}</Name>
      <NormalizedInstanceHours>{{ cluster.normalized_instance_hours }}</NormalizedInstanceHours>
      {% if cluster.release_label is not none %}
      <ReleaseLabel>{{ cluster.release_label }}</ReleaseLabel>
      {% endif %}
      {% if cluster.requested_ami_version is not none %}
      <RequestedAmiVersion>{{ cluster.requested_ami_version }}</RequestedAmiVersion>
      {% endif %}
      {% if cluster.running_ami_version is not none %}
      <RunningAmiVersion>{{ cluster.running_ami_version }}</RunningAmiVersion>
      {% endif %}
      {% if cluster.security_configuration is not none %}
      <SecurityConfiguration>{{ cluster.security_configuration }}</SecurityConfiguration>
      {% endif %}
      <ServiceRole>{{ cluster.service_role }}</ServiceRole>
      <AutoScalingRole>{{ cluster.auto_scaling_role }}</AutoScalingRole>
      <Status>
        <State>{{ cluster.state }}</State>
        <StateChangeReason>
          {% if cluster.last_state_change_reason is not none %}
          <Message>{{ cluster.last_state_change_reason }}</Message>
          {% endif %}
          <Code>USER_REQUEST</Code>
        </StateChangeReason>
        <Timeline>
          <CreationDateTime>{{ cluster.creation_datetime.isoformat() }}</CreationDateTime>
          {% if cluster.end_datetime is not none %}
          <EndDateTime>{{ cluster.end_datetime.isoformat() }}</EndDateTime>
          {% endif %}
          {% if cluster.ready_datetime is not none %}
          <ReadyDateTime>{{ cluster.ready_datetime.isoformat() }}</ReadyDateTime>
          {% endif %}
        </Timeline>
      </Status>
      <Tags>
        {% for tag_key, tag_value in cluster.tags.items() %}
        <member>
          <Key>{{ tag_key }}</Key>
          <Value>{{ tag_value }}</Value>
        </member>
        {% endfor %}
      </Tags>
      <TerminationProtected>{{ cluster.termination_protected|lower }}</TerminationProtected>
      <VisibleToAllUsers>{{ cluster.visible_to_all_users|lower }}</VisibleToAllUsers>
      <StepConcurrencyLevel>{{ cluster.step_concurrency_level }}</StepConcurrencyLevel>
      <ClusterArn>{{ cluster.arn }}</ClusterArn>
    </Cluster>
  </DescribeClusterResult>
  <ResponseMetadata>
    <RequestId>aaaaaaaa-bbbb-cccc-dddd-eeeeeeeeeeee</RequestId>
  </ResponseMetadata>
</DescribeClusterResponse>"""

DESCRIBE_JOB_FLOWS_TEMPLATE = """<DescribeJobFlowsResponse xmlns="http://elasticmapreduce.amazonaws.com/doc/2009-03-31">
  <DescribeJobFlowsResult>
    <JobFlows>
      {% for cluster in clusters %}
      <member>
        {% if cluster.running_ami_version is not none %}
        <AmiVersion>{{ cluster.running_ami_version }}</AmiVersion>
        {% endif %}
        {% if cluster.bootstrap_actions %}
        <BootstrapActions>
          {% for bootstrap_action in cluster.bootstrap_actions %}
          <member>
            <BootstrapActionConfig>
              <Name>{{ bootstrap_action.name }}</Name>
              <ScriptBootstrapAction>
                <Args>
                  {% for arg in bootstrap_action.args %}
                  <member>{{ arg | escape }}</member>
                  {% endfor %}
                </Args>
                <Path>{{ bootstrap_action.script_path | escape }}</Path>
              </ScriptBootstrapAction>
            </BootstrapActionConfig>
          </member>
          {% endfor %}
        </BootstrapActions>
        {% endif %}
        <ExecutionStatusDetail>
          <CreationDateTime>{{ cluster.creation_datetime.isoformat() }}</CreationDateTime>
          {% if cluster.end_datetime is not none %}
          <EndDateTime>{{ cluster.end_datetime.isoformat() }}</EndDateTime>
          {% endif %}
          {% if cluster.last_state_change_reason is not none %}
          <LastStateChangeReason>{{ cluster.last_state_change_reason }}</LastStateChangeReason>
          {% endif %}
          {% if cluster.ready_datetime is not none %}
          <ReadyDateTime>{{ cluster.ready_datetime.isoformat() }}</ReadyDateTime>
          {% endif %}
          {% if cluster.start_datetime is not none %}
          <StartDateTime>{{ cluster.start_datetime.isoformat() }}</StartDateTime>
          {% endif %}
          <State>{{ cluster.state }}</State>
        </ExecutionStatusDetail>
        <Instances>
          {% if cluster.ec2_key_name is not none %}
          <Ec2KeyName>{{ cluster.ec2_key_name }}</Ec2KeyName>
          {% endif %}
          {% if cluster.ec2_subnet_id is not none %}
          <Ec2SubnetId>{{ cluster.ec2_subnet_id }}</Ec2SubnetId>
          {% endif %}
          <HadoopVersion>{{ cluster.hadoop_version }}</HadoopVersion>
          <InstanceCount>{{ cluster.instance_count }}</InstanceCount>
          <InstanceGroups>
            {% for instance_group in cluster.instance_groups %}
            <member>
              {% if instance_group.bid_price is not none %}
              <BidPrice>{{ instance_group.bid_price }}</BidPrice>
              {% endif %}
              <CreationDateTime>{{ instance_group.creation_datetime.isoformat() }}</CreationDateTime>
              {% if instance_group.end_datetime is not none %}
              <EndDateTime>{{ instance_group.end_datetime.isoformat() }}</EndDateTime>
              {% endif %}

              <InstanceGroupId>{{ instance_group.id }}</InstanceGroupId>
              <InstanceRequestCount>{{ instance_group.num_instances }}</InstanceRequestCount>
              <InstanceRole>{{ instance_group.role }}</InstanceRole>
              <InstanceRunningCount>{{ instance_group.num_instances }}</InstanceRunningCount>
              <InstanceType>{{ instance_group.instance_type }}</InstanceType>
              <LastStateChangeReason/>
              <Market>{{ instance_group.market }}</Market>
              <Name>{{ instance_group.name }}</Name>
              {% if instance_group.ready_datetime is not none %}
              <ReadyDateTime>{{ instance_group.ready_datetime.isoformat() }}</ReadyDateTime>
              {% endif %}
              {% if instance_group.start_datetime is not none %}
              <StartDateTime>{{ instance_group.start_datetime.isoformat() }}</StartDateTime>
              {% endif %}
              <State>{{ instance_group.state }}</State>
            </member>
            {% endfor %}
          </InstanceGroups>
          <KeepJobFlowAliveWhenNoSteps>{{ cluster.keep_job_flow_alive_when_no_steps|lower }}</KeepJobFlowAliveWhenNoSteps>
          <MasterInstanceId>{{ cluster.master_instance_id }}</MasterInstanceId>
          <MasterInstanceType>{{ cluster.master_instance_type }}</MasterInstanceType>
          <MasterPublicDnsName>ec2-184-0-0-1.{{ cluster.region }}.compute.amazonaws.com</MasterPublicDnsName>
          <NormalizedInstanceHours>{{ cluster.normalized_instance_hours }}</NormalizedInstanceHours>
          <Placement>
            <AvailabilityZone>{{ cluster.availability_zone }}</AvailabilityZone>
          </Placement>
          <SlaveInstanceType>{{ cluster.slave_instance_type }}</SlaveInstanceType>
          <TerminationProtected>{{ cluster.termination_protected|lower }}</TerminationProtected>
        </Instances>
        <JobFlowId>{{ cluster.id }}</JobFlowId>
        <JobFlowRole>{{ cluster.role }}</JobFlowRole>
        <LogUri>{{ cluster.log_uri }}</LogUri>
        <Name>{{ cluster.name }}</Name>
        <ServiceRole>{{ cluster.service_role }}</ServiceRole>
        <Steps>
          {% for step in cluster.steps %}
          <member>
            <ExecutionStatusDetail>
              <CreationDateTime>{{ step.creation_datetime.isoformat() }}</CreationDateTime>
              {% if step.end_datetime is not none %}
              <EndDateTime>{{ step.end_datetime.isoformat() }}</EndDateTime>
              {% endif %}
              {% if step.last_state_change_reason is not none %}
              <LastStateChangeReason>{{ step.last_state_change_reason }}</LastStateChangeReason>
              {% endif %}
              {% if step.ready_datetime is not none %}
              <ReadyDateTime>{{ step.ready_datetime.isoformat() }}</ReadyDateTime>
              {% endif %}
              {% if step.start_datetime is not none %}
              <StartDateTime>{{ step.start_datetime.isoformat() }}</StartDateTime>
              {% endif %}
              <State>{{ step.state }}</State>
            </ExecutionStatusDetail>
            <StepConfig>
              <ActionOnFailure>{{ step.action_on_failure }}</ActionOnFailure>
              <HadoopJarStep>
                <Jar>{{ step.jar }}</Jar>
                <MainClass>{{ step.main_class }}</MainClass>
                <Args>
                  {% for arg in step.args %}
                  <member>{{ arg | escape }}</member>
                  {% endfor %}
                </Args>
                <Properties>
                {% for key, val in step.properties.items() %}
                  <member>
                    <Key>{{ key }}</Key>
                    <Value>{{ val | escape }}</Value>
                  </member>
                {% endfor %}
                </Properties>
              </HadoopJarStep>
              <Name>{{ step.name | escape }}</Name>
            </StepConfig>
          </member>
          {% endfor %}
        </Steps>
        <SupportedProducts/>
        <VisibleToAllUsers>{{ cluster.visible_to_all_users|lower }}</VisibleToAllUsers>
      </member>
      {% endfor %}
    </JobFlows>
  </DescribeJobFlowsResult>
  <ResponseMetadata>
    <RequestId>9cea3229-ed85-11dd-9877-6fad448a8419</RequestId>
  </ResponseMetadata>
</DescribeJobFlowsResponse>"""

DESCRIBE_STEP_TEMPLATE = """<DescribeStepResponse xmlns="http://elasticmapreduce.amazonaws.com/doc/2009-03-31">
  <DescribeStepResult>
    <Step>
      <ActionOnFailure>{{ step.action_on_failure }}</ActionOnFailure>
      <Config>
        <Args>
          {% for arg in step.args %}
          <member>{{ arg | escape }}</member>
          {% endfor %}
        </Args>
        <Jar>{{ step.jar }}</Jar>
        <MainClass/>
        <Properties>
          {% for key, val in step.properties.items() %}
          <entry>
            <key>{{ key }}</key>
            <value>{{ val | escape }}</value>
          </entry>
          {% endfor %}
        </Properties>
      </Config>
      <Id>{{ step.id }}</Id>
      <Name>{{ step.name | escape }}</Name>
      <Status>
        <FailureDetails>
          <Reason/>
          <Message/>
          <LogFile/>
        </FailureDetails>
        <State>{{ step.state }}</State>
        <StateChangeReason>{{ step.state_change_reason }}</StateChangeReason>
        <Timeline>
          <CreationDateTime>{{ step.creation_datetime.isoformat() }}</CreationDateTime>
          {% if step.end_datetime is not none %}
          <EndDateTime>{{ step.end_datetime.isoformat() }}</EndDateTime>
          {% endif %}
          {% if step.ready_datetime is not none %}
          <StartDateTime>{{ step.start_datetime.isoformat() }}</StartDateTime>
          {% endif %}
        </Timeline>
      </Status>
    </Step>
  </DescribeStepResult>
  <ResponseMetadata>
    <RequestId>df6f4f4a-ed85-11dd-9877-6fad448a8419</RequestId>
  </ResponseMetadata>
</DescribeStepResponse>"""

LIST_BOOTSTRAP_ACTIONS_TEMPLATE = """<ListBootstrapActionsResponse xmlns="http://elasticmapreduce.amazonaws.com/doc/2009-03-31">
  <ListBootstrapActionsResult>
    <BootstrapActions>
      {% for bootstrap_action in bootstrap_actions %}
      <member>
        <Args>
          {% for arg in bootstrap_action.args %}
          <member>{{ arg | escape }}</member>
          {% endfor %}
        </Args>
        <Name>{{ bootstrap_action.name }}</Name>
        <ScriptPath>{{ bootstrap_action.script_path }}</ScriptPath>
      </member>
      {% endfor %}
    </BootstrapActions>
    {% if marker is not none %}
    <Marker>{{ marker }}</Marker>
    {% endif %}
  </ListBootstrapActionsResult>
  <ResponseMetadata>
    <RequestId>df6f4f4a-ed85-11dd-9877-6fad448a8419</RequestId>
  </ResponseMetadata>
</ListBootstrapActionsResponse>"""

LIST_CLUSTERS_TEMPLATE = """<ListClustersResponse xmlns="http://elasticmapreduce.amazonaws.com/doc/2009-03-31">
  <ListClustersResult>
    <Clusters>
      {% for cluster in clusters %}
      <member>
        <Id>{{ cluster.id }}</Id>
        <Name>{{ cluster.name }}</Name>
        <NormalizedInstanceHours>{{ cluster.normalized_instance_hours }}</NormalizedInstanceHours>
        <Status>
          <State>{{ cluster.state }}</State>
          <StateChangeReason>
            <Code>USER_REQUEST</Code>
            {% if cluster.last_state_change_reason is not none %}
            <Message>{{ cluster.last_state_change_reason }}</Message>
            {% endif %}
          </StateChangeReason>
          <Timeline>
            <CreationDateTime>{{ cluster.creation_datetime.isoformat() }}</CreationDateTime>
            {% if cluster.end_datetime is not none %}
            <EndDateTime>{{ cluster.end_datetime.isoformat() }}</EndDateTime>
            {% endif %}
            {% if cluster.ready_datetime is not none %}
            <ReadyDateTime>{{ cluster.ready_datetime.isoformat() }}</ReadyDateTime>
            {% endif %}
          </Timeline>
        </Status>
        <ClusterArn>{{ cluster.arn }}</ClusterArn>
      </member>
      {% endfor %}
    </Clusters>
    {% if marker is not none %}
    <Marker>{{ marker }}</Marker>
    {% endif %}
  </ListClustersResult>
  <ResponseMetadata>
    <RequestId>2690d7eb-ed86-11dd-9877-6fad448a8418</RequestId>
  </ResponseMetadata>
</ListClustersResponse>"""

LIST_INSTANCE_GROUPS_TEMPLATE = """<ListInstanceGroupsResponse xmlns="http://elasticmapreduce.amazonaws.com/doc/2009-03-31">
  <ListInstanceGroupsResult>
    <InstanceGroups>
      {% for instance_group in instance_groups %}
      <member>
        {% if instance_group.bid_price is not none %}
        <BidPrice>{{ instance_group.bid_price }}</BidPrice>
        {% endif %}
        <Configurations/>
        {% if instance_group.ebs_configuration is not none %}
        <EbsBlockDevices>
            {% for ebs_block_device in instance_group.ebs_configuration.ebs_block_device_configs %}
              {% for i in range(ebs_block_device.volumes_per_instance) %}
          <member>
            <VolumeSpecification>
              <VolumeType>{{ebs_block_device.volume_specification.volume_type}}</VolumeType>
              <Iops>{{ebs_block_device.volume_specification.iops}}</Iops>
              <SizeInGB>{{ebs_block_device.volume_specification.size_in_gb}}</SizeInGB>
            </VolumeSpecification>
            <Device>/dev/sd{{i}}</Device>
          </member>
              {% endfor %}
            {% endfor %}
        </EbsBlockDevices>
        {% endif %}
        {% if instance_group.auto_scaling_policy is not none %}
        <AutoScalingPolicy>
            {% if instance_group.auto_scaling_policy.constraints is not none %}
            <Constraints>
                {% if instance_group.auto_scaling_policy.constraints.min_capacity is not none %}
                <MinCapacity>{{instance_group.auto_scaling_policy.constraints.min_capacity}}</MinCapacity>
                {% endif %}
                {% if instance_group.auto_scaling_policy.constraints.max_capacity is not none %}
                <MaxCapacity>{{instance_group.auto_scaling_policy.constraints.max_capacity}}</MaxCapacity>
                {% endif %}
            </Constraints>
            {% endif %}
            {% if instance_group.auto_scaling_policy.rules is not none %}
            <Rules>
                {% for rule in instance_group.auto_scaling_policy.rules %}
                <member>
                    {% if 'name' in rule %}
                    <Name>{{rule['name']}}</Name>
                    {% endif %}
                    {% if 'description' in rule %}
                    <Description>{{rule['description']}}</Description>
                    {% endif %}
                    {% if 'action' in rule %}
                    <Action>
                        {% if 'market' in rule['action'] %}
                        <Market>{{rule['action']['market']}}</Market>
                        {% endif %}
                        {% if 'simple_scaling_policy_configuration' in rule['action'] %}
                        <SimpleScalingPolicyConfiguration>
                            {% if 'adjustment_type' in rule['action']['simple_scaling_policy_configuration'] %}
                            <AdjustmentType>{{rule['action']['simple_scaling_policy_configuration']['adjustment_type']}}</AdjustmentType>
                            {% endif %}
                            {% if 'scaling_adjustment' in rule['action']['simple_scaling_policy_configuration'] %}
                            <ScalingAdjustment>{{rule['action']['simple_scaling_policy_configuration']['scaling_adjustment']}}</ScalingAdjustment>
                            {% endif %}
                            {% if 'cool_down' in rule['action']['simple_scaling_policy_configuration'] %}
                            <CoolDown>{{rule['action']['simple_scaling_policy_configuration']['cool_down']}}</CoolDown>
                            {% endif %}
                        </SimpleScalingPolicyConfiguration>
                        {% endif %}
                    </Action>
                    {% endif %}
                    {% if 'trigger' in rule %}
                    <Trigger>
                        {% if 'cloud_watch_alarm_definition' in rule['trigger'] %}
                        <CloudWatchAlarmDefinition>
                            {% if 'comparison_operator' in rule['trigger']['cloud_watch_alarm_definition'] %}
                            <ComparisonOperator>{{rule['trigger']['cloud_watch_alarm_definition']['comparison_operator']}}</ComparisonOperator>
                            {% endif %}
                            {% if 'evaluation_periods' in rule['trigger']['cloud_watch_alarm_definition'] %}
                            <EvaluationPeriods>{{rule['trigger']['cloud_watch_alarm_definition']['evaluation_periods']}}</EvaluationPeriods>
                            {% endif %}
                            {% if 'metric_name' in rule['trigger']['cloud_watch_alarm_definition'] %}
                            <MetricName>{{rule['trigger']['cloud_watch_alarm_definition']['metric_name']}}</MetricName>
                            {% endif %}
                            {% if 'namespace' in rule['trigger']['cloud_watch_alarm_definition'] %}
                            <Namespace>{{rule['trigger']['cloud_watch_alarm_definition']['namespace']}}</Namespace>
                            {% endif %}
                            {% if 'period' in rule['trigger']['cloud_watch_alarm_definition'] %}
                            <Period>{{rule['trigger']['cloud_watch_alarm_definition']['period']}}</Period>
                            {% endif %}
                            {% if 'statistic' in rule['trigger']['cloud_watch_alarm_definition'] %}
                            <Statistic>{{rule['trigger']['cloud_watch_alarm_definition']['statistic']}}</Statistic>
                            {% endif %}
                            {% if 'threshold' in rule['trigger']['cloud_watch_alarm_definition'] %}
                            <Threshold>{{rule['trigger']['cloud_watch_alarm_definition']['threshold']}}</Threshold>
                            {% endif %}
                            {% if 'unit' in rule['trigger']['cloud_watch_alarm_definition'] %}
                            <Unit>{{rule['trigger']['cloud_watch_alarm_definition']['unit']}}</Unit>
                            {% endif %}
                            {% if 'dimensions' in rule['trigger']['cloud_watch_alarm_definition'] %}
                            <Dimensions>
                                {% for dimension in rule['trigger']['cloud_watch_alarm_definition']['dimensions'] %}
                                <member>
                                    {% if 'key' in dimension %}
                                    <Key>{{dimension['key']}}</Key>
                                    {% endif %}
                                    {% if 'value' in dimension %}
                                    <Value>{{dimension['value']}}</Value>
                                    {% endif %}
                                </member>
                                {% endfor %}
                            </Dimensions>
                            {% endif %}
                        </CloudWatchAlarmDefinition>
                        {% endif %}
                    </Trigger>
                    {% endif %}
                </member>
                {% endfor %}
            </Rules>
            {% endif %}
            {% if instance_group.auto_scaling_policy.status is not none %}
            <Status>
                {% if 'state' in instance_group.auto_scaling_policy.status %}
                <State>{{instance_group.auto_scaling_policy.status['state']}}</State>
                {% endif %}
            </Status>
            {% endif %}
        </AutoScalingPolicy>
        {% endif %}
        {% if instance_group.ebs_optimized is not none %}
        <EbsOptimized>{{ instance_group.ebs_optimized }}</EbsOptimized>
        {% endif %}
        <Id>{{ instance_group.id }}</Id>
        <InstanceGroupType>{{ instance_group.role }}</InstanceGroupType>
        <InstanceType>{{ instance_group.instance_type }}</InstanceType>
        <Market>{{ instance_group.market }}</Market>
        <Name>{{ instance_group.name }}</Name>
        <RequestedInstanceCount>{{ instance_group.num_instances }}</RequestedInstanceCount>
        <RunningInstanceCount>{{ instance_group.num_instances }}</RunningInstanceCount>
        <Status>
          <State>{{ instance_group.state }}</State>
          <StateChangeReason>
            {% if instance_group.state_change_reason is not none %}
            <Message>{{ instance_group.state_change_reason }}</Message>
            {% endif %}
            <Code>USER_REQUEST</Code>
          </StateChangeReason>
          <Timeline>
            <CreationDateTime>{{ instance_group.creation_datetime.isoformat() }}</CreationDateTime>
            {% if instance_group.end_datetime is not none %}
            <EndDateTime>{{ instance_group.end_datetime.isoformat() }}</EndDateTime>
            {% endif %}
            {% if instance_group.ready_datetime is not none %}
            <ReadyDateTime>{{ instance_group.ready_datetime.isoformat() }}</ReadyDateTime>
            {% endif %}
          </Timeline>
        </Status>
      </member>
      {% endfor %}
    </InstanceGroups>
    {% if marker is not none %}
    <Marker>{{ marker }}</Marker>
    {% endif %}
  </ListInstanceGroupsResult>
  <ResponseMetadata>
    <RequestId>8296d8b8-ed85-11dd-9877-6fad448a8419</RequestId>
  </ResponseMetadata>
</ListInstanceGroupsResponse>"""

LIST_INSTANCES_TEMPLATE = """<ListInstancesResponse xmlns="http://elasticmapreduce.amazonaws.com/doc/2009-03-31">
  <ListInstancesResult>
    <Instances>
     {% for instance in instances %}
      <member>
        <Id>{{ instance.id }}</Id>
        <Ec2InstanceId>{{ instance.ec2_instance_id }}</Ec2InstanceId>
        <PublicDnsName>{{ instance.details.public_dns }}</PublicDnsName>
        <PublicIpAddress>{{ instance.details.public_ip }}</PublicIpAddress>
        <PrivateDnsName>{{ instance.details.private_dns }}</PrivateDnsName>
        <PrivateIpAddress>{{ instance.details.private_ip }}</PrivateIpAddress>
        <InstanceGroupId>{{ instance.instance_group.id }}</InstanceGroupId>
        <InstanceFleetId>{{ instance.instance_fleet_id }}</InstanceFleetId>
        <Market>{{ instance.instance_group.market }}</Market>
        <InstanceType>{{ instance.details.instance_type }}</InstanceType>
         <EbsVolumes>
              {% for volume in instance.details.block_device_mapping %}
          <member>
              <Device>{{ volume }}</Device>
              <VolumeId>{{ instance.details.block_device_mapping[volume].volume_id }}</VolumeId>
          </member>
              {% endfor %}
        </EbsVolumes>
       <Status>
          <State>{{ instance.instance_group.state }}</State>
          <StateChangeReason>
            {% if instance.state_change_reason is not none %}
            <Message>{{ instance.state_change_reason }}</Message>
            {% endif %}
          </StateChangeReason>
          <Timeline>
            <CreationDateTime>{{ instance.instance_group.creation_datetime.isoformat() }}</CreationDateTime>
            {% if instance.instance_group.end_datetime is not none %}
            <EndDateTime>{{ instance.instance_group.end_datetime.isoformat() }}</EndDateTime>
            {% endif %}
            {% if instance.instance_group.ready_datetime is not none %}
            <ReadyDateTime>{{ instance.instance_group.ready_datetime.isoformat() }}</ReadyDateTime>
            {% endif %}
          </Timeline>
        </Status>
      </member>
    {% endfor %}
    </Instances>
 </ListInstancesResult>
 <ResponseMetadata>
    <RequestId>4248c46c-71c0-4772-b155-0e992dc30027</RequestId>
  </ResponseMetadata>
</ListInstancesResponse>"""

LIST_STEPS_TEMPLATE = """<ListStepsResponse xmlns="http://elasticmapreduce.amazonaws.com/doc/2009-03-31">
  <ListStepsResult>
    <Steps>
      {% for step in steps %}
      <member>
        <ActionOnFailure>{{ step.action_on_failure }}</ActionOnFailure>
        <Config>
          <Args>
            {% for arg in step.args %}
            <member>{{ arg | escape }}</member>
            {% endfor %}
          </Args>
          <Jar>{{ step.jar | escape }}</Jar>
          <MainClass/>
          <Properties>
            {% for key, val in step.properties.items() %}
            <entry>
              <key>{{ key }}</key>
              <value>{{ val | escape }}</value>
            </entry>
            {% endfor %}
          </Properties>
        </Config>
        <Id>{{ step.id }}</Id>
        <Name>{{ step.name | escape }}</Name>
        <Status>
<!-- does not exist for botocore 1.4.28
          <FailureDetails>
            <Reason/>
            <Message/>
            <LogFile/>
          </FailureDetails>
-->
          <State>{{ step.state }}</State>
          <StateChangeReason>{{ step.state_change_reason }}</StateChangeReason>
          <Timeline>
            <CreationDateTime>{{ step.creation_datetime.isoformat() }}</CreationDateTime>
            {% if step.end_datetime is not none %}
            <EndDateTime>{{ step.end_datetime.isoformat() }}</EndDateTime>
            {% endif %}
            {% if step.start_datetime is not none %}
            <StartDateTime>{{ step.start_datetime.isoformat() }}</StartDateTime>
            {% endif %}
          </Timeline>
        </Status>
      </member>
      {% endfor %}
    </Steps>
    {% if marker is not none %}
    <Marker>{{ marker }}</Marker>
    {% endif %}
  </ListStepsResult>
  <ResponseMetadata>
    <RequestId>df6f4f4a-ed85-11dd-9877-6fad448a8419</RequestId>
  </ResponseMetadata>
</ListStepsResponse>"""

MODIFY_CLUSTER_TEMPLATE = """<ModifyClusterResponse xmlns="http://elasticmapreduce.amazonaws.com/doc/2009-03-31">
  <ModifyClusterResult>
    <StepConcurrencyLevel>{{ cluster.step_concurrency_level }}</StepConcurrencyLevel>
  </ModifyClusterResult>
  <ResponseMetadata>
    <RequestId>0751c837-e78d-4aef-95c9-9c4d29a092ff</RequestId>
  </ResponseMetadata>
</ModifyClusterResponse>
"""

MODIFY_INSTANCE_GROUPS_TEMPLATE = """<ModifyInstanceGroupsResponse xmlns="http://elasticmapreduce.amazonaws.com/doc/2009-03-31">
  <ResponseMetadata>
    <RequestId>2690d7eb-ed86-11dd-9877-6fad448a8419</RequestId>
  </ResponseMetadata>
</ModifyInstanceGroupsResponse>"""

REMOVE_TAGS_TEMPLATE = """<RemoveTagsResponse xmlns="http://elasticmapreduce.amazonaws.com/doc/2009-03-31">
  <ResponseMetadata>
    <RequestId>2690d7eb-ed86-11dd-9877-6fad448a8419</RequestId>
  </ResponseMetadata>
</RemoveTagsResponse>"""

RUN_JOB_FLOW_TEMPLATE = """<RunJobFlowResponse xmlns="http://elasticmapreduce.amazonaws.com/doc/2009-03-31">
  <RunJobFlowResult>
    <JobFlowId>{{ cluster.id }}</JobFlowId>
    <ClusterArn>{{ cluster.arn }}</ClusterArn>
  </RunJobFlowResult>
  <ResponseMetadata>
    <RequestId>8296d8b8-ed85-11dd-9877-6fad448a8419</RequestId>
  </ResponseMetadata>
</RunJobFlowResponse>"""

SET_TERMINATION_PROTECTION_TEMPLATE = """<SetTerminationProtection xmlns="http://elasticmapreduce.amazonaws.com/doc/2009-03-31">
  <ResponseMetadata>
    <RequestId>2690d7eb-ed86-11dd-9877-6fad448a8419</RequestId>
  </ResponseMetadata>
</SetTerminationProtection>"""

SET_VISIBLE_TO_ALL_USERS_TEMPLATE = """<SetVisibleToAllUsersResponse xmlns="http://elasticmapreduce.amazonaws.com/doc/2009-03-31">
  <ResponseMetadata>
    <RequestId>2690d7eb-ed86-11dd-9877-6fad448a8419</RequestId>
  </ResponseMetadata>
</SetVisibleToAllUsersResponse>"""

TERMINATE_JOB_FLOWS_TEMPLATE = """<TerminateJobFlowsResponse xmlns="http://elasticmapreduce.amazonaws.com/doc/2009-03-31">
  <ResponseMetadata>
    <RequestId>2690d7eb-ed86-11dd-9877-6fad448a8419</RequestId>
  </ResponseMetadata>
</TerminateJobFlowsResponse>"""

PUT_AUTO_SCALING_POLICY = """<PutAutoScalingPolicyResponse xmlns="http://elasticmapreduce.amazonaws.com/doc/2009-03-31">
  <PutAutoScalingPolicyResult>
    <ClusterId>{{cluster_id}}</ClusterId>
    <InstanceGroupId>{{instance_group.id}}</InstanceGroupId>
    {% if instance_group.auto_scaling_policy is not none %}
    <AutoScalingPolicy>
        {% if instance_group.auto_scaling_policy.constraints is not none %}
        <Constraints>
            {% if instance_group.auto_scaling_policy.constraints.min_capacity is not none %}
            <MinCapacity>{{instance_group.auto_scaling_policy.constraints.min_capacity}}</MinCapacity>
            {% endif %}
            {% if instance_group.auto_scaling_policy.constraints.max_capacity is not none %}
            <MaxCapacity>{{instance_group.auto_scaling_policy.constraints.max_capacity}}</MaxCapacity>
            {% endif %}
        </Constraints>
        {% endif %}
        {% if instance_group.auto_scaling_policy.rules is not none %}
        <Rules>
            {% for rule in instance_group.auto_scaling_policy.rules %}
            <member>
                {% if 'name' in rule %}
                <Name>{{rule['name']}}</Name>
                {% endif %}
                {% if 'description' in rule %}
                <Description>{{rule['description']}}</Description>
                {% endif %}
                {% if 'action' in rule %}
                <Action>
                    {% if 'market' in rule['action'] %}
                    <Market>{{rule['action']['market']}}</Market>
                    {% endif %}
                    {% if 'simple_scaling_policy_configuration' in rule['action'] %}
                    <SimpleScalingPolicyConfiguration>
                        {% if 'adjustment_type' in rule['action']['simple_scaling_policy_configuration'] %}
                        <AdjustmentType>{{rule['action']['simple_scaling_policy_configuration']['adjustment_type']}}</AdjustmentType>
                        {% endif %}
                        {% if 'scaling_adjustment' in rule['action']['simple_scaling_policy_configuration'] %}
                        <ScalingAdjustment>{{rule['action']['simple_scaling_policy_configuration']['scaling_adjustment']}}</ScalingAdjustment>
                        {% endif %}
                        {% if 'cool_down' in rule['action']['simple_scaling_policy_configuration'] %}
                        <CoolDown>{{rule['action']['simple_scaling_policy_configuration']['cool_down']}}</CoolDown>
                        {% endif %}
                    </SimpleScalingPolicyConfiguration>
                    {% endif %}
                </Action>
                {% endif %}
                {% if 'trigger' in rule %}
                <Trigger>
                    {% if 'cloud_watch_alarm_definition' in rule['trigger'] %}
                    <CloudWatchAlarmDefinition>
                        {% if 'comparison_operator' in rule['trigger']['cloud_watch_alarm_definition'] %}
                        <ComparisonOperator>{{rule['trigger']['cloud_watch_alarm_definition']['comparison_operator']}}</ComparisonOperator>
                        {% endif %}
                        {% if 'evaluation_periods' in rule['trigger']['cloud_watch_alarm_definition'] %}
                        <EvaluationPeriods>{{rule['trigger']['cloud_watch_alarm_definition']['evaluation_periods']}}</EvaluationPeriods>
                        {% endif %}
                        {% if 'metric_name' in rule['trigger']['cloud_watch_alarm_definition'] %}
                        <MetricName>{{rule['trigger']['cloud_watch_alarm_definition']['metric_name']}}</MetricName>
                        {% endif %}
                        {% if 'namespace' in rule['trigger']['cloud_watch_alarm_definition'] %}
                        <Namespace>{{rule['trigger']['cloud_watch_alarm_definition']['namespace']}}</Namespace>
                        {% endif %}
                        {% if 'period' in rule['trigger']['cloud_watch_alarm_definition'] %}
                        <Period>{{rule['trigger']['cloud_watch_alarm_definition']['period']}}</Period>
                        {% endif %}
                        {% if 'statistic' in rule['trigger']['cloud_watch_alarm_definition'] %}
                        <Statistic>{{rule['trigger']['cloud_watch_alarm_definition']['statistic']}}</Statistic>
                        {% endif %}
                        {% if 'threshold' in rule['trigger']['cloud_watch_alarm_definition'] %}
                        <Threshold>{{rule['trigger']['cloud_watch_alarm_definition']['threshold']}}</Threshold>
                        {% endif %}
                        {% if 'unit' in rule['trigger']['cloud_watch_alarm_definition'] %}
                        <Unit>{{rule['trigger']['cloud_watch_alarm_definition']['unit']}}</Unit>
                        {% endif %}
                        {% if 'dimensions' in rule['trigger']['cloud_watch_alarm_definition'] %}
                        <Dimensions>
                            {% for dimension in rule['trigger']['cloud_watch_alarm_definition']['dimensions'] %}
                            <member>
                                {% if 'key' in dimension %}
                                <Key>{{dimension['key']}}</Key>
                                {% endif %}
                                {% if 'value' in dimension %}
                                <Value>{{dimension['value']}}</Value>
                                {% endif %}
                            </member>
                            {% endfor %}
                        </Dimensions>
                        {% endif %}
                    </CloudWatchAlarmDefinition>
                    {% endif %}
                </Trigger>
                {% endif %}
            </member>
            {% endfor %}
        </Rules>
        {% endif %}
        {% if instance_group.auto_scaling_policy.status is not none %}
        <Status>
            {% if 'state' in instance_group.auto_scaling_policy.status %}
            <State>{{instance_group.auto_scaling_policy.status['state']}}</State>
            {% endif %}
        </Status>
        {% endif %}
    </AutoScalingPolicy>
    {% endif %}
    <ClusterArn>{{ cluster.arn }}</ClusterArn>
  </PutAutoScalingPolicyResult>
  <ResponseMetadata>
    <RequestId>d47379d9-b505-49af-9335-a68950d82535</RequestId>
  </ResponseMetadata>
</PutAutoScalingPolicyResponse>"""

REMOVE_AUTO_SCALING_POLICY = """<RemoveAutoScalingPolicyResponse xmlns="http://elasticmapreduce.amazonaws.com/doc/2009-03-31">
  <ResponseMetadata>
    <RequestId>c04a1042-5340-4c0a-a7b5-7779725ce4f7</RequestId>
  </ResponseMetadata>
</RemoveAutoScalingPolicyResponse>"""

CREATE_SECURITY_CONFIGURATION_TEMPLATE = """<CreateSecurityConfigurationResponse xmlns="http://elasticmapreduce.amazonaws.com/doc/2009-03-31">
  <CreateSecurityConfigurationResult>
    <Name>{{name}}</Name>
    <CreationDateTime>{{creation_date_time}}</CreationDateTime>
  </CreateSecurityConfigurationResult>
  <ResponseMetadata>
    <RequestId>2690d7eb-ed86-11dd-9877-6fad448a8419</RequestId>
  </ResponseMetadata>
</CreateSecurityConfigurationResponse>"""

DESCRIBE_SECURITY_CONFIGURATION_TEMPLATE = """<DescribeSecurityConfigurationResponse xmlns="http://elasticmapreduce.amazonaws.com/doc/2009-03-31">
  <DescribeSecurityConfigurationResult>
    <Name>{{security_configuration['name']}}</Name>
    <SecurityConfiguration>{{security_configuration['security_configuration']}}</SecurityConfiguration>
    <CreationDateTime>{{security_configuration['creation_date_time']}}</CreationDateTime>
  </DescribeSecurityConfigurationResult>
  <ResponseMetadata>
    <RequestId>2690d7eb-ed86-11dd-9877-6fad448a8419</RequestId>
  </ResponseMetadata>
</DescribeSecurityConfigurationResponse>"""

DELETE_SECURITY_CONFIGURATION_TEMPLATE = """<DeleteSecurityConfigurationResponse xmlns="http://elasticmapreduce.amazonaws.com/doc/2009-03-31">
  <ResponseMetadata>
    <RequestId>2690d7eb-ed86-11dd-9877-6fad448a8419</RequestId>
  </ResponseMetadata>
</DeleteSecurityConfigurationResponse>"""

PUT_BLOCK_PUBLIC_ACCESS_CONFIGURATION_TEMPLATE = """<PutBlockPublicAccessConfigurationResponse xmlns="http://elasticmapreduce.amazonaws.com/doc/2009-03-31">
  <ResponseMetadata>
    <RequestId>2690d7eb-ed86-11dd-9877-6fad448a8419</RequestId>
  </ResponseMetadata>
</PutBlockPublicAccessConfigurationResponse>"""

GET_BLOCK_PUBLIC_ACCESS_CONFIGURATION_TEMPLATE = """
  <GetBlockPublicAccessConfigurationResponse xmlns="http://elasticmapreduce.amazonaws.com/doc/2009-03-31">
    <GetBlockPublicAccessConfigurationResult>
      <BlockPublicAccessConfiguration>
        <BlockPublicSecurityGroupRules>
          {{block_public_security_group_rules}}
        </BlockPublicSecurityGroupRules>
        <PermittedPublicSecurityGroupRuleRanges>
          {% for rule_range in permitted_public_security_group_rule_ranges %}
            <member>
              <MinRange>{{rule_range['min_range']}}</MinRange>
              <MaxRange>{{rule_range['max_range']}}</MaxRange>      
            </member>
          {% endfor %}
        </PermittedPublicSecurityGroupRuleRanges>
      </BlockPublicAccessConfiguration>
      <BlockPublicAccessConfigurationMetadata>
        <CreationDateTime>{{creation_date_time}}</CreationDateTime>
        <CreatedByArn>{{created_by_arn}}</CreatedByArn>
      </BlockPublicAccessConfigurationMetadata>
    </GetBlockPublicAccessConfigurationResult>
    <ResponseMetadata>
      <RequestId>2690d7eb-ed86-11dd-9877-6fad448a8419</RequestId>
    </ResponseMetadata>
  </GetBlockPublicAccessConfigurationResponse>
"""
