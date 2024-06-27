from collections import defaultdict
from typing import Any, Dict, Iterable, List

from moto.core.common_types import TYPE_RESPONSE
from moto.core.responses import BaseResponse
from moto.ec2.models import ec2_backends
from moto.neptune.responses import (
    CREATE_GLOBAL_CLUSTER_TEMPLATE,
    DELETE_GLOBAL_CLUSTER_TEMPLATE,
    DESCRIBE_GLOBAL_CLUSTERS_TEMPLATE,
    REMOVE_FROM_GLOBAL_CLUSTER_TEMPLATE,
    NeptuneResponse,
)

from .exceptions import DBParameterGroupNotFoundError
from .models import RDSBackend, rds_backends


class RDSResponse(BaseResponse):
    def __init__(self) -> None:
        super().__init__(service_name="rds")
        # Neptune and RDS share a HTTP endpoint RDS is the lucky guy that catches all requests
        # So we have to determine whether we can handle an incoming request here, or whether it needs redirecting to Neptune
        self.neptune = NeptuneResponse()

    @property
    def backend(self) -> RDSBackend:
        return rds_backends[self.current_account][self.region]

    def _dispatch(self, request: Any, full_url: str, headers: Any) -> TYPE_RESPONSE:
        # Because some requests are send through to Neptune, we have to prepare the NeptuneResponse-class
        self.neptune.setup_class(request, full_url, headers)
        return super()._dispatch(request, full_url, headers)

    def __getattribute__(self, name: str) -> Any:
        if name in ["create_db_cluster", "create_global_cluster"]:
            if self._get_param("Engine") == "neptune":
                return object.__getattribute__(self.neptune, name)
        return object.__getattribute__(self, name)

    def _get_db_kwargs(self) -> Dict[str, Any]:
        args = {
            "auto_minor_version_upgrade": self._get_param("AutoMinorVersionUpgrade"),
            "allocated_storage": self._get_int_param("AllocatedStorage"),
            "availability_zone": self._get_param("AvailabilityZone"),
            "backup_retention_period": self._get_param("BackupRetentionPeriod"),
            "copy_tags_to_snapshot": self._get_param("CopyTagsToSnapshot"),
            "db_instance_class": self._get_param("DBInstanceClass"),
            "db_cluster_identifier": self._get_param("DBClusterIdentifier"),
            "db_instance_identifier": self._get_param("DBInstanceIdentifier"),
            "db_name": self._get_param("DBName"),
            "db_parameter_group_name": self._get_param("DBParameterGroupName"),
            "db_snapshot_identifier": self._get_param("DBSnapshotIdentifier"),
            "db_subnet_group_name": self._get_param("DBSubnetGroupName"),
            "engine": self._get_param("Engine"),
            "engine_version": self._get_param("EngineVersion"),
            "enable_cloudwatch_logs_exports": self._get_params().get(
                "EnableCloudwatchLogsExports"
            ),
            "enable_iam_database_authentication": self._get_bool_param(
                "EnableIAMDatabaseAuthentication"
            ),
            "license_model": self._get_param("LicenseModel"),
            "iops": self._get_int_param("Iops"),
            "kms_key_id": self._get_param("KmsKeyId"),
            "master_user_password": self._get_param("MasterUserPassword"),
            "master_username": self._get_param("MasterUsername"),
            "multi_az": self._get_bool_param("MultiAZ"),
            "option_group_name": self._get_param("OptionGroupName"),
            "port": self._get_param("Port"),
            "preferred_backup_window": self._get_param(
                "PreferredBackupWindow", "13:14-13:44"
            ),
            "preferred_maintenance_window": self._get_param(
                "PreferredMaintenanceWindow", "wed:06:38-wed:07:08"
            ).lower(),
            "publicly_accessible": self._get_param("PubliclyAccessible"),
            "account_id": self.current_account,
            "region": self.region,
            "security_groups": self._get_multi_param(
                "DBSecurityGroups.DBSecurityGroupName"
            ),
            "storage_encrypted": self._get_param("StorageEncrypted"),
            "storage_type": self._get_param("StorageType", None),
            "vpc_security_group_ids": self._get_multi_param(
                "VpcSecurityGroupIds.VpcSecurityGroupId"
            ),
            "tags": list(),
            "deletion_protection": self._get_bool_param("DeletionProtection"),
        }
        args["tags"] = self.unpack_list_params("Tags", "Tag")
        return args

    def _get_modify_db_cluster_kwargs(self) -> Dict[str, Any]:
        args = {
            "auto_minor_version_upgrade": self._get_param("AutoMinorVersionUpgrade"),
            "allocated_storage": self._get_int_param("AllocatedStorage"),
            "availability_zone": self._get_param("AvailabilityZone"),
            "backup_retention_period": self._get_param("BackupRetentionPeriod"),
            "copy_tags_to_snapshot": self._get_param("CopyTagsToSnapshot"),
            "db_instance_class": self._get_param("DBInstanceClass"),
            "db_cluster_identifier": self._get_param("DBClusterIdentifier"),
            "new_db_cluster_identifier": self._get_param("NewDBClusterIdentifier"),
            "db_instance_identifier": self._get_param("DBInstanceIdentifier"),
            "db_name": self._get_param("DBName"),
            "db_parameter_group_name": self._get_param("DBParameterGroupName"),
            "db_cluster_parameter_group_name": self._get_param(
                "DBClusterParameterGroupName"
            ),
            "db_snapshot_identifier": self._get_param("DBSnapshotIdentifier"),
            "db_subnet_group_name": self._get_param("DBSubnetGroupName"),
            "engine": self._get_param("Engine"),
            "engine_version": self._get_param("EngineVersion"),
            "enable_cloudwatch_logs_exports": self._get_params().get(
                "CloudwatchLogsExportConfiguration"
            ),
            "enable_iam_database_authentication": self._get_bool_param(
                "EnableIAMDatabaseAuthentication"
            ),
            "enable_http_endpoint": self._get_bool_param("EnableHttpEndpoint"),
            "license_model": self._get_param("LicenseModel"),
            "iops": self._get_int_param("Iops"),
            "kms_key_id": self._get_param("KmsKeyId"),
            "master_user_password": self._get_param("MasterUserPassword"),
            "master_username": self._get_param("MasterUsername"),
            "multi_az": self._get_bool_param("MultiAZ"),
            "option_group_name": self._get_param("OptionGroupName"),
            "port": self._get_param("Port"),
            "preferred_backup_window": self._get_param("PreferredBackupWindow"),
            "preferred_maintenance_window": self._get_param(
                "PreferredMaintenanceWindow"
            ),
            "publicly_accessible": self._get_param("PubliclyAccessible"),
            "account_id": self.current_account,
            "region": self.region,
            "security_groups": self._get_multi_param(
                "DBSecurityGroups.DBSecurityGroupName"
            ),
            "storage_encrypted": self._get_param("StorageEncrypted"),
            "storage_type": self._get_param("StorageType", None),
            "vpc_security_group_ids": self._get_multi_param(
                "VpcSecurityGroupIds.VpcSecurityGroupId"
            ),
            "tags": list(),
            "deletion_protection": self._get_bool_param("DeletionProtection"),
        }
        args["tags"] = self.unpack_list_params("Tags", "Tag")
        return args

    def _get_db_replica_kwargs(self) -> Dict[str, Any]:
        return {
            "auto_minor_version_upgrade": self._get_param("AutoMinorVersionUpgrade"),
            "availability_zone": self._get_param("AvailabilityZone"),
            "db_instance_class": self._get_param("DBInstanceClass"),
            "db_instance_identifier": self._get_param("DBInstanceIdentifier"),
            "db_subnet_group_name": self._get_param("DBSubnetGroupName"),
            "iops": self._get_int_param("Iops"),
            # OptionGroupName
            "port": self._get_param("Port"),
            "publicly_accessible": self._get_param("PubliclyAccessible"),
            "source_db_identifier": self._get_param("SourceDBInstanceIdentifier"),
            "storage_type": self._get_param("StorageType"),
        }

    def _get_option_group_kwargs(self) -> Dict[str, Any]:
        return {
            "major_engine_version": self._get_param("MajorEngineVersion"),
            "description": self._get_param("OptionGroupDescription"),
            "engine_name": self._get_param("EngineName"),
            "name": self._get_param("OptionGroupName"),
        }

    def _get_db_parameter_group_kwargs(self) -> Dict[str, Any]:
        return {
            "description": self._get_param("Description"),
            "family": self._get_param("DBParameterGroupFamily"),
            "name": self._get_param("DBParameterGroupName"),
            "tags": self.unpack_list_params("Tags", "Tag"),
        }

    def _get_db_cluster_kwargs(self) -> Dict[str, Any]:
        params = self._get_params()
        return {
            "availability_zones": self._get_multi_param(
                "AvailabilityZones.AvailabilityZone"
            ),
            "enable_cloudwatch_logs_exports": params.get("EnableCloudwatchLogsExports"),
            "db_name": self._get_param("DatabaseName"),
            "db_cluster_identifier": self._get_param("DBClusterIdentifier"),
            "db_subnet_group_name": self._get_param("DBSubnetGroupName"),
            "deletion_protection": self._get_bool_param("DeletionProtection"),
            "engine": self._get_param("Engine"),
            "engine_version": self._get_param("EngineVersion"),
            "engine_mode": self._get_param("EngineMode"),
            "allocated_storage": self._get_param("AllocatedStorage"),
            "global_cluster_identifier": self._get_param("GlobalClusterIdentifier"),
            "iops": self._get_param("Iops"),
            "storage_encrypted": self._get_param("StorageEncrypted"),
            "enable_global_write_forwarding": self._get_param(
                "EnableGlobalWriteForwarding"
            ),
            "storage_type": self._get_param("StorageType"),
            "kms_key_id": self._get_param("KmsKeyId"),
            "master_username": self._get_param("MasterUsername"),
            "master_user_password": self._get_param("MasterUserPassword"),
            "network_type": self._get_param("NetworkType"),
            "port": self._get_param("Port"),
            "parameter_group": self._get_param("DBClusterParameterGroupName"),
            "region": self.region,
            "db_cluster_instance_class": self._get_param("DBClusterInstanceClass"),
            "enable_http_endpoint": self._get_bool_param("EnableHttpEndpoint"),
            "copy_tags_to_snapshot": self._get_param("CopyTagsToSnapshot"),
            "tags": self.unpack_list_params("Tags", "Tag"),
            "scaling_configuration": self._get_dict_param("ScalingConfiguration."),
            "serverless_v2_scaling_configuration": params.get(
                "ServerlessV2ScalingConfiguration"
            ),
            "replication_source_identifier": self._get_param(
                "ReplicationSourceIdentifier"
            ),
            "vpc_security_group_ids": self.unpack_list_params(
                "VpcSecurityGroupIds", "VpcSecurityGroupId"
            ),
        }

    def _get_export_task_kwargs(self) -> Dict[str, Any]:
        return {
            "export_task_identifier": self._get_param("ExportTaskIdentifier"),
            "source_arn": self._get_param("SourceArn"),
            "s3_bucket_name": self._get_param("S3BucketName"),
            "iam_role_arn": self._get_param("IamRoleArn"),
            "kms_key_id": self._get_param("KmsKeyId"),
            "s3_prefix": self._get_param("S3Prefix"),
            "export_only": self.unpack_list_params("ExportOnly", "member"),
        }

    def _get_event_subscription_kwargs(self) -> Dict[str, Any]:
        return {
            "subscription_name": self._get_param("SubscriptionName"),
            "sns_topic_arn": self._get_param("SnsTopicArn"),
            "source_type": self._get_param("SourceType"),
            "event_categories": self.unpack_list_params(
                "EventCategories", "EventCategory"
            ),
            "source_ids": self.unpack_list_params("SourceIds", "SourceId"),
            "enabled": self._get_param("Enabled"),
            "tags": self.unpack_list_params("Tags", "Tag"),
        }

    def unpack_list_params(self, label: str, child_label: str) -> List[Dict[str, Any]]:
        root = self._get_multi_param_dict(label) or {}
        return root.get(child_label, [])

    def create_db_instance(self) -> str:
        db_kwargs = self._get_db_kwargs()
        database = self.backend.create_db_instance(db_kwargs)
        template = self.response_template(CREATE_DATABASE_TEMPLATE)
        return template.render(database=database)

    def create_db_instance_read_replica(self) -> str:
        db_kwargs = self._get_db_replica_kwargs()

        database = self.backend.create_db_instance_read_replica(db_kwargs)
        template = self.response_template(CREATE_DATABASE_REPLICA_TEMPLATE)
        return template.render(database=database)

    def describe_db_instances(self) -> str:
        db_instance_identifier = self._get_param("DBInstanceIdentifier")
        filters = self._get_multi_param("Filters.Filter.")
        filter_dict = {f["Name"]: f["Values"] for f in filters}
        all_instances = list(
            self.backend.describe_db_instances(
                db_instance_identifier, filters=filter_dict
            )
        )
        marker = self._get_param("Marker")
        all_ids = [instance.db_instance_identifier for instance in all_instances]
        if marker:
            start = all_ids.index(marker) + 1
        else:
            start = 0
        page_size = self._get_int_param(
            "MaxRecords", 50
        )  # the default is 100, but using 50 to make testing easier
        instances_resp = all_instances[start : start + page_size]
        next_marker = None
        if len(all_instances) > start + page_size:
            next_marker = instances_resp[-1].db_instance_identifier

        template = self.response_template(DESCRIBE_DATABASES_TEMPLATE)
        return template.render(databases=instances_resp, marker=next_marker)

    def modify_db_instance(self) -> str:
        db_instance_identifier = self._get_param("DBInstanceIdentifier")
        db_kwargs = self._get_db_kwargs()
        # NOTE modify_db_instance does not support tags
        del db_kwargs["tags"]
        new_db_instance_identifier = self._get_param("NewDBInstanceIdentifier")
        if new_db_instance_identifier:
            db_kwargs["new_db_instance_identifier"] = new_db_instance_identifier
        database = self.backend.modify_db_instance(db_instance_identifier, db_kwargs)
        template = self.response_template(MODIFY_DATABASE_TEMPLATE)
        return template.render(database=database)

    def delete_db_instance(self) -> str:
        db_instance_identifier = self._get_param("DBInstanceIdentifier")
        db_snapshot_name = self._get_param("FinalDBSnapshotIdentifier")
        database = self.backend.delete_db_instance(
            db_instance_identifier, db_snapshot_name
        )
        template = self.response_template(DELETE_DATABASE_TEMPLATE)
        return template.render(database=database)

    def reboot_db_instance(self) -> str:
        db_instance_identifier = self._get_param("DBInstanceIdentifier")
        database = self.backend.reboot_db_instance(db_instance_identifier)
        template = self.response_template(REBOOT_DATABASE_TEMPLATE)
        return template.render(database=database)

    def create_db_snapshot(self) -> str:
        db_instance_identifier = self._get_param("DBInstanceIdentifier")
        db_snapshot_identifier = self._get_param("DBSnapshotIdentifier")
        tags = self.unpack_list_params("Tags", "Tag")
        snapshot = self.backend.create_db_snapshot(
            db_instance_identifier, db_snapshot_identifier, tags=tags
        )
        template = self.response_template(CREATE_SNAPSHOT_TEMPLATE)
        return template.render(snapshot=snapshot)

    def copy_db_snapshot(self) -> str:
        source_snapshot_identifier = self._get_param("SourceDBSnapshotIdentifier")
        target_snapshot_identifier = self._get_param("TargetDBSnapshotIdentifier")
        tags = self.unpack_list_params("Tags", "Tag")
        snapshot = self.backend.copy_db_snapshot(
            source_snapshot_identifier, target_snapshot_identifier, tags
        )
        template = self.response_template(COPY_SNAPSHOT_TEMPLATE)
        return template.render(snapshot=snapshot)

    def describe_db_snapshots(self) -> str:
        db_instance_identifier = self._get_param("DBInstanceIdentifier")
        db_snapshot_identifier = self._get_param("DBSnapshotIdentifier")
        filters = self._get_multi_param("Filters.Filter.")
        filter_dict = {f["Name"]: f["Values"] for f in filters}
        snapshots = self.backend.describe_db_snapshots(
            db_instance_identifier, db_snapshot_identifier, filter_dict
        )
        template = self.response_template(DESCRIBE_SNAPSHOTS_TEMPLATE)
        return template.render(snapshots=snapshots)

    def promote_read_replica(self) -> str:
        db_instance_identifier = self._get_param("DBInstanceIdentifier")
        db_kwargs = self._get_db_kwargs()
        database = self.backend.promote_read_replica(db_kwargs)
        database = self.backend.modify_db_instance(db_instance_identifier, db_kwargs)
        template = self.response_template(PROMOTE_REPLICA_TEMPLATE)
        return template.render(database=database)

    def delete_db_snapshot(self) -> str:
        db_snapshot_identifier = self._get_param("DBSnapshotIdentifier")
        snapshot = self.backend.delete_db_snapshot(db_snapshot_identifier)
        template = self.response_template(DELETE_SNAPSHOT_TEMPLATE)
        return template.render(snapshot=snapshot)

    def restore_db_instance_from_db_snapshot(self) -> str:
        db_snapshot_identifier = self._get_param("DBSnapshotIdentifier")
        db_kwargs = self._get_db_kwargs()
        new_instance = self.backend.restore_db_instance_from_db_snapshot(
            db_snapshot_identifier, db_kwargs
        )
        template = self.response_template(RESTORE_INSTANCE_FROM_SNAPSHOT_TEMPLATE)
        return template.render(database=new_instance)

    def restore_db_instance_to_point_in_time(self) -> str:
        source_db_identifier = self._get_param("SourceDBInstanceIdentifier")
        target_db_identifier = self._get_param("TargetDBInstanceIdentifier")

        db_kwargs = self._get_db_kwargs()
        new_instance = self.backend.restore_db_instance_to_point_in_time(
            source_db_identifier, target_db_identifier, db_kwargs
        )
        template = self.response_template(RESTORE_INSTANCE_TO_POINT_IN_TIME_TEMPLATE)
        return template.render(database=new_instance)

    def list_tags_for_resource(self) -> str:
        arn = self._get_param("ResourceName")
        template = self.response_template(LIST_TAGS_FOR_RESOURCE_TEMPLATE)
        tags = self.backend.list_tags_for_resource(arn)
        return template.render(tags=tags)

    def add_tags_to_resource(self) -> str:
        arn = self._get_param("ResourceName")
        tags = self.unpack_list_params("Tags", "Tag")
        tags = self.backend.add_tags_to_resource(arn, tags)
        template = self.response_template(ADD_TAGS_TO_RESOURCE_TEMPLATE)
        return template.render(tags=tags)

    def remove_tags_from_resource(self) -> str:
        arn = self._get_param("ResourceName")
        tag_keys = self.unpack_list_params("TagKeys", "member")
        self.backend.remove_tags_from_resource(arn, tag_keys)  # type: ignore
        template = self.response_template(REMOVE_TAGS_FROM_RESOURCE_TEMPLATE)
        return template.render()

    def stop_db_instance(self) -> str:
        db_instance_identifier = self._get_param("DBInstanceIdentifier")
        db_snapshot_identifier = self._get_param("DBSnapshotIdentifier")
        database = self.backend.stop_db_instance(
            db_instance_identifier, db_snapshot_identifier
        )
        template = self.response_template(STOP_DATABASE_TEMPLATE)
        return template.render(database=database)

    def start_db_instance(self) -> str:
        db_instance_identifier = self._get_param("DBInstanceIdentifier")
        database = self.backend.start_db_instance(db_instance_identifier)
        template = self.response_template(START_DATABASE_TEMPLATE)
        return template.render(database=database)

    def create_db_security_group(self) -> str:
        group_name = self._get_param("DBSecurityGroupName")
        description = self._get_param("DBSecurityGroupDescription")
        tags = self.unpack_list_params("Tags", "Tag")
        security_group = self.backend.create_db_security_group(
            group_name, description, tags
        )
        template = self.response_template(CREATE_SECURITY_GROUP_TEMPLATE)
        return template.render(security_group=security_group)

    def describe_db_security_groups(self) -> str:
        security_group_name = self._get_param("DBSecurityGroupName")
        security_groups = self.backend.describe_security_groups(security_group_name)
        template = self.response_template(DESCRIBE_SECURITY_GROUPS_TEMPLATE)
        return template.render(security_groups=security_groups)

    def delete_db_security_group(self) -> str:
        security_group_name = self._get_param("DBSecurityGroupName")
        security_group = self.backend.delete_security_group(security_group_name)
        template = self.response_template(DELETE_SECURITY_GROUP_TEMPLATE)
        return template.render(security_group=security_group)

    def authorize_db_security_group_ingress(self) -> str:
        security_group_name = self._get_param("DBSecurityGroupName")
        cidr_ip = self._get_param("CIDRIP")
        security_group = self.backend.authorize_security_group(
            security_group_name, cidr_ip
        )
        template = self.response_template(AUTHORIZE_SECURITY_GROUP_TEMPLATE)
        return template.render(security_group=security_group)

    def create_db_subnet_group(self) -> str:
        subnet_name = self._get_param("DBSubnetGroupName")
        description = self._get_param("DBSubnetGroupDescription")
        subnet_ids = self._get_multi_param("SubnetIds.SubnetIdentifier")
        tags = self.unpack_list_params("Tags", "Tag")
        subnets = [
            ec2_backends[self.current_account][self.region].get_subnet(subnet_id)
            for subnet_id in subnet_ids
        ]
        subnet_group = self.backend.create_subnet_group(
            subnet_name, description, subnets, tags
        )
        template = self.response_template(CREATE_SUBNET_GROUP_TEMPLATE)
        return template.render(subnet_group=subnet_group)

    def describe_db_subnet_groups(self) -> str:
        subnet_name = self._get_param("DBSubnetGroupName")
        subnet_groups = self.backend.describe_db_subnet_groups(subnet_name)
        template = self.response_template(DESCRIBE_SUBNET_GROUPS_TEMPLATE)
        return template.render(subnet_groups=subnet_groups)

    def modify_db_subnet_group(self) -> str:
        subnet_name = self._get_param("DBSubnetGroupName")
        description = self._get_param("DBSubnetGroupDescription")
        subnet_ids = self._get_multi_param("SubnetIds.SubnetIdentifier")
        subnets = [
            ec2_backends[self.current_account][self.region].get_subnet(subnet_id)
            for subnet_id in subnet_ids
        ]
        subnet_group = self.backend.modify_db_subnet_group(
            subnet_name, description, subnets
        )
        template = self.response_template(MODIFY_SUBNET_GROUPS_TEMPLATE)
        return template.render(subnet_group=subnet_group)

    def delete_db_subnet_group(self) -> str:
        subnet_name = self._get_param("DBSubnetGroupName")
        subnet_group = self.backend.delete_subnet_group(subnet_name)
        template = self.response_template(DELETE_SUBNET_GROUP_TEMPLATE)
        return template.render(subnet_group=subnet_group)

    def create_option_group(self) -> str:
        kwargs = self._get_option_group_kwargs()
        option_group = self.backend.create_option_group(kwargs)
        template = self.response_template(CREATE_OPTION_GROUP_TEMPLATE)
        return template.render(option_group=option_group)

    def delete_option_group(self) -> str:
        kwargs = self._get_option_group_kwargs()
        option_group = self.backend.delete_option_group(kwargs["name"])
        template = self.response_template(DELETE_OPTION_GROUP_TEMPLATE)
        return template.render(option_group=option_group)

    def describe_option_groups(self) -> str:
        kwargs = self._get_option_group_kwargs()
        kwargs["max_records"] = self._get_int_param("MaxRecords")
        kwargs["marker"] = self._get_param("Marker")
        option_groups = self.backend.describe_option_groups(kwargs)
        template = self.response_template(DESCRIBE_OPTION_GROUP_TEMPLATE)
        return template.render(option_groups=option_groups)

    def describe_option_group_options(self) -> str:
        engine_name = self._get_param("EngineName")
        major_engine_version = self._get_param("MajorEngineVersion")
        return self.backend.describe_option_group_options(
            engine_name, major_engine_version
        )

    def modify_option_group(self) -> str:
        option_group_name = self._get_param("OptionGroupName")
        count = 1
        options_to_include = []
        # TODO: This can probably be refactored with a single call to super.get_multi_param, but there are not enough tests (yet) to verify this
        while self._get_param(f"OptionsToInclude.member.{count}.OptionName"):
            options_to_include.append(
                {
                    "Port": self._get_param(f"OptionsToInclude.member.{count}.Port"),
                    "OptionName": self._get_param(
                        f"OptionsToInclude.member.{count}.OptionName"
                    ),
                    "DBSecurityGroupMemberships": self._get_param(
                        f"OptionsToInclude.member.{count}.DBSecurityGroupMemberships"
                    ),
                    "OptionSettings": self._get_param(
                        f"OptionsToInclude.member.{count}.OptionSettings"
                    ),
                    "VpcSecurityGroupMemberships": self._get_param(
                        f"OptionsToInclude.member.{count}.VpcSecurityGroupMemberships"
                    ),
                }
            )
            count += 1

        count = 1
        options_to_remove = []
        while self._get_param(f"OptionsToRemove.member.{count}"):
            options_to_remove.append(self._get_param(f"OptionsToRemove.member.{count}"))
            count += 1
        option_group = self.backend.modify_option_group(
            option_group_name, options_to_include, options_to_remove
        )
        template = self.response_template(MODIFY_OPTION_GROUP_TEMPLATE)
        return template.render(option_group=option_group)

    def create_db_parameter_group(self) -> str:
        kwargs = self._get_db_parameter_group_kwargs()
        db_parameter_group = self.backend.create_db_parameter_group(kwargs)
        template = self.response_template(CREATE_DB_PARAMETER_GROUP_TEMPLATE)
        return template.render(db_parameter_group=db_parameter_group)

    def describe_db_parameter_groups(self) -> str:
        kwargs = self._get_db_parameter_group_kwargs()
        kwargs["max_records"] = self._get_int_param("MaxRecords")
        kwargs["marker"] = self._get_param("Marker")
        db_parameter_groups = self.backend.describe_db_parameter_groups(kwargs)
        template = self.response_template(DESCRIBE_DB_PARAMETER_GROUPS_TEMPLATE)
        return template.render(db_parameter_groups=db_parameter_groups)

    def modify_db_parameter_group(self) -> str:
        db_parameter_group_name = self._get_param("DBParameterGroupName")
        db_parameter_group_parameters = self._get_db_parameter_group_parameters()
        db_parameter_group = self.backend.modify_db_parameter_group(
            db_parameter_group_name, db_parameter_group_parameters
        )
        template = self.response_template(MODIFY_DB_PARAMETER_GROUP_TEMPLATE)
        return template.render(db_parameter_group=db_parameter_group)

    def _get_db_parameter_group_parameters(self) -> Iterable[Dict[str, Any]]:
        parameter_group_parameters: Dict[str, Any] = defaultdict(dict)
        for param_name, value in self.querystring.items():
            if not param_name.startswith("Parameters.Parameter"):
                continue

            split_param_name = param_name.split(".")
            param_id = split_param_name[2]
            param_setting = split_param_name[3]

            parameter_group_parameters[param_id][param_setting] = value[0]

        return parameter_group_parameters.values()

    def describe_db_parameters(self) -> str:
        db_parameter_group_name = self._get_param("DBParameterGroupName")
        db_parameter_groups = self.backend.describe_db_parameter_groups(
            {"name": db_parameter_group_name}
        )
        if not db_parameter_groups:
            raise DBParameterGroupNotFoundError(db_parameter_group_name)

        template = self.response_template(DESCRIBE_DB_PARAMETERS_TEMPLATE)
        return template.render(db_parameter_group=db_parameter_groups[0])

    def delete_db_parameter_group(self) -> str:
        kwargs = self._get_db_parameter_group_kwargs()
        db_parameter_group = self.backend.delete_db_parameter_group(kwargs["name"])
        template = self.response_template(DELETE_DB_PARAMETER_GROUP_TEMPLATE)
        return template.render(db_parameter_group=db_parameter_group)

    def describe_db_cluster_parameters(self) -> str:
        db_parameter_group_name = self._get_param("DBParameterGroupName")
        db_parameter_groups = self.backend.describe_db_cluster_parameters()
        if db_parameter_groups is None:
            raise DBParameterGroupNotFoundError(db_parameter_group_name)

        template = self.response_template(DESCRIBE_DB_CLUSTER_PARAMETERS_TEMPLATE)
        return template.render(db_parameter_group=db_parameter_groups)

    def create_db_cluster(self) -> str:
        kwargs = self._get_db_cluster_kwargs()
        cluster = self.backend.create_db_cluster(kwargs)
        template = self.response_template(CREATE_DB_CLUSTER_TEMPLATE)
        return template.render(cluster=cluster)

    def modify_db_cluster(self) -> str:
        kwargs = self._get_modify_db_cluster_kwargs()
        cluster = self.backend.modify_db_cluster(kwargs)
        template = self.response_template(MODIFY_DB_CLUSTER_TEMPLATE)
        return template.render(cluster=cluster)

    def describe_db_clusters(self) -> str:
        _id = self._get_param("DBClusterIdentifier")
        filters = self._get_multi_param("Filters.Filter.")
        filter_dict = {f["Name"]: f["Values"] for f in filters}
        clusters = self.backend.describe_db_clusters(
            cluster_identifier=_id, filters=filter_dict
        )
        template = self.response_template(DESCRIBE_CLUSTERS_TEMPLATE)
        return template.render(clusters=clusters)

    def delete_db_cluster(self) -> str:
        _id = self._get_param("DBClusterIdentifier")
        snapshot_name = self._get_param("FinalDBSnapshotIdentifier")
        cluster = self.backend.delete_db_cluster(
            cluster_identifier=_id, snapshot_name=snapshot_name
        )
        template = self.response_template(DELETE_CLUSTER_TEMPLATE)
        return template.render(cluster=cluster)

    def start_db_cluster(self) -> str:
        _id = self._get_param("DBClusterIdentifier")
        cluster = self.backend.start_db_cluster(cluster_identifier=_id)
        template = self.response_template(START_CLUSTER_TEMPLATE)
        return template.render(cluster=cluster)

    def stop_db_cluster(self) -> str:
        _id = self._get_param("DBClusterIdentifier")
        cluster = self.backend.stop_db_cluster(cluster_identifier=_id)
        template = self.response_template(STOP_CLUSTER_TEMPLATE)
        return template.render(cluster=cluster)

    def create_db_cluster_snapshot(self) -> str:
        db_cluster_identifier = self._get_param("DBClusterIdentifier")
        db_snapshot_identifier = self._get_param("DBClusterSnapshotIdentifier")
        tags = self.unpack_list_params("Tags", "Tag")
        snapshot = self.backend.create_db_cluster_snapshot(
            db_cluster_identifier, db_snapshot_identifier, tags=tags
        )
        template = self.response_template(CREATE_CLUSTER_SNAPSHOT_TEMPLATE)
        return template.render(snapshot=snapshot)

    def copy_db_cluster_snapshot(self) -> str:
        source_snapshot_identifier = self._get_param(
            "SourceDBClusterSnapshotIdentifier"
        )
        target_snapshot_identifier = self._get_param(
            "TargetDBClusterSnapshotIdentifier"
        )
        tags = self.unpack_list_params("Tags", "Tag")
        snapshot = self.backend.copy_db_cluster_snapshot(
            source_snapshot_identifier, target_snapshot_identifier, tags
        )
        template = self.response_template(COPY_CLUSTER_SNAPSHOT_TEMPLATE)
        return template.render(snapshot=snapshot)

    def describe_db_cluster_snapshots(self) -> str:
        db_cluster_identifier = self._get_param("DBClusterIdentifier")
        db_snapshot_identifier = self._get_param("DBClusterSnapshotIdentifier")
        filters = self._get_multi_param("Filters.Filter.")
        filter_dict = {f["Name"]: f["Values"] for f in filters}
        snapshots = self.backend.describe_db_cluster_snapshots(
            db_cluster_identifier, db_snapshot_identifier, filter_dict
        )
        template = self.response_template(DESCRIBE_CLUSTER_SNAPSHOTS_TEMPLATE)
        return template.render(snapshots=snapshots)

    def delete_db_cluster_snapshot(self) -> str:
        db_snapshot_identifier = self._get_param("DBClusterSnapshotIdentifier")
        snapshot = self.backend.delete_db_cluster_snapshot(db_snapshot_identifier)
        template = self.response_template(DELETE_CLUSTER_SNAPSHOT_TEMPLATE)
        return template.render(snapshot=snapshot)

    def restore_db_cluster_from_snapshot(self) -> str:
        db_snapshot_identifier = self._get_param("SnapshotIdentifier")
        db_kwargs = self._get_db_cluster_kwargs()
        new_cluster = self.backend.restore_db_cluster_from_snapshot(
            db_snapshot_identifier, db_kwargs
        )
        template = self.response_template(RESTORE_CLUSTER_FROM_SNAPSHOT_TEMPLATE)
        return template.render(cluster=new_cluster)

    def start_export_task(self) -> str:
        kwargs = self._get_export_task_kwargs()
        export_task = self.backend.start_export_task(kwargs)
        template = self.response_template(START_EXPORT_TASK_TEMPLATE)
        return template.render(task=export_task)

    def cancel_export_task(self) -> str:
        export_task_identifier = self._get_param("ExportTaskIdentifier")
        export_task = self.backend.cancel_export_task(export_task_identifier)
        template = self.response_template(CANCEL_EXPORT_TASK_TEMPLATE)
        return template.render(task=export_task)

    def describe_export_tasks(self) -> str:
        export_task_identifier = self._get_param("ExportTaskIdentifier")
        tasks = self.backend.describe_export_tasks(export_task_identifier)
        template = self.response_template(DESCRIBE_EXPORT_TASKS_TEMPLATE)
        return template.render(tasks=tasks)

    def create_event_subscription(self) -> str:
        kwargs = self._get_event_subscription_kwargs()
        subscription = self.backend.create_event_subscription(kwargs)
        template = self.response_template(CREATE_EVENT_SUBSCRIPTION_TEMPLATE)
        return template.render(subscription=subscription)

    def delete_event_subscription(self) -> str:
        subscription_name = self._get_param("SubscriptionName")
        subscription = self.backend.delete_event_subscription(subscription_name)
        template = self.response_template(DELETE_EVENT_SUBSCRIPTION_TEMPLATE)
        return template.render(subscription=subscription)

    def describe_event_subscriptions(self) -> str:
        subscription_name = self._get_param("SubscriptionName")
        subscriptions = self.backend.describe_event_subscriptions(subscription_name)
        template = self.response_template(DESCRIBE_EVENT_SUBSCRIPTIONS_TEMPLATE)
        return template.render(subscriptions=subscriptions)

    def describe_orderable_db_instance_options(self) -> str:
        engine = self._get_param("Engine")
        engine_version = self._get_param("EngineVersion")
        options = self.backend.describe_orderable_db_instance_options(
            engine, engine_version
        )
        template = self.response_template(DESCRIBE_ORDERABLE_CLUSTER_OPTIONS)
        return template.render(options=options, marker=None)

    def describe_global_clusters(self) -> str:
        clusters = self.backend.describe_global_clusters()
        template = self.response_template(DESCRIBE_GLOBAL_CLUSTERS_TEMPLATE)
        return template.render(clusters=clusters)

    def create_global_cluster(self) -> str:
        params = self._get_params()
        cluster = self.backend.create_global_cluster(
            global_cluster_identifier=params["GlobalClusterIdentifier"],
            source_db_cluster_identifier=params.get("SourceDBClusterIdentifier"),
            engine=params.get("Engine"),
            engine_version=params.get("EngineVersion"),
            storage_encrypted=params.get("StorageEncrypted"),
            deletion_protection=params.get("DeletionProtection"),
        )
        template = self.response_template(CREATE_GLOBAL_CLUSTER_TEMPLATE)
        return template.render(cluster=cluster)

    def delete_global_cluster(self) -> str:
        params = self._get_params()
        cluster = self.backend.delete_global_cluster(
            global_cluster_identifier=params["GlobalClusterIdentifier"],
        )
        template = self.response_template(DELETE_GLOBAL_CLUSTER_TEMPLATE)
        return template.render(cluster=cluster)

    def remove_from_global_cluster(self) -> str:
        params = self._get_params()
        global_cluster = self.backend.remove_from_global_cluster(
            global_cluster_identifier=params["GlobalClusterIdentifier"],
            db_cluster_identifier=params["DbClusterIdentifier"],
        )
        template = self.response_template(REMOVE_FROM_GLOBAL_CLUSTER_TEMPLATE)
        return template.render(cluster=global_cluster)

    def create_db_cluster_parameter_group(self) -> str:
        group_name = self._get_param("DBClusterParameterGroupName")
        family = self._get_param("DBParameterGroupFamily")
        desc = self._get_param("Description")
        db_cluster_parameter_group = self.backend.create_db_cluster_parameter_group(
            group_name=group_name,
            family=family,
            description=desc,
        )
        template = self.response_template(CREATE_DB_CLUSTER_PARAMETER_GROUP_TEMPLATE)
        return template.render(db_cluster_parameter_group=db_cluster_parameter_group)

    def describe_db_cluster_parameter_groups(self) -> str:
        group_name = self._get_param("DBClusterParameterGroupName")
        db_parameter_groups = self.backend.describe_db_cluster_parameter_groups(
            group_name=group_name,
        )
        template = self.response_template(DESCRIBE_DB_CLUSTER_PARAMETER_GROUPS_TEMPLATE)
        return template.render(db_parameter_groups=db_parameter_groups)

    def delete_db_cluster_parameter_group(self) -> str:
        group_name = self._get_param("DBClusterParameterGroupName")
        self.backend.delete_db_cluster_parameter_group(
            group_name=group_name,
        )
        template = self.response_template(DELETE_DB_CLUSTER_PARAMETER_GROUP_TEMPLATE)
        return template.render()

    def promote_read_replica_db_cluster(self) -> str:
        db_cluster_identifier = self._get_param("DBClusterIdentifier")
        cluster = self.backend.promote_read_replica_db_cluster(db_cluster_identifier)
        template = self.response_template(PROMOTE_READ_REPLICA_DB_CLUSTER_TEMPLATE)
        return template.render(cluster=cluster)

    def describe_db_snapshot_attributes(self) -> str:
        params = self._get_params()
        db_snapshot_identifier = params["DBSnapshotIdentifier"]
        db_snapshot_attributes_result = self.backend.describe_db_snapshot_attributes(
            db_snapshot_identifier=db_snapshot_identifier,
        )
        template = self.response_template(DESCRIBE_DB_SNAPSHOT_ATTRIBUTES_TEMPLATE)
        return template.render(
            db_snapshot_attributes_result=db_snapshot_attributes_result,
            db_snapshot_identifier=db_snapshot_identifier,
        )

    def modify_db_snapshot_attribute(self) -> str:
        params = self._get_params()
        db_snapshot_identifier = params["DBSnapshotIdentifier"]
        db_snapshot_attributes_result = self.backend.modify_db_snapshot_attribute(
            db_snapshot_identifier=db_snapshot_identifier,
            attribute_name=params["AttributeName"],
            values_to_add=params.get("ValuesToAdd"),
            values_to_remove=params.get("ValuesToRemove"),
        )
        template = self.response_template(MODIFY_DB_SNAPSHOT_ATTRIBUTE_TEMPLATE)
        return template.render(
            db_snapshot_attributes_result=db_snapshot_attributes_result,
            db_snapshot_identifier=db_snapshot_identifier,
        )

    def describe_db_cluster_snapshot_attributes(self) -> str:
        params = self._get_params()
        db_cluster_snapshot_identifier = params["DBClusterSnapshotIdentifier"]
        db_cluster_snapshot_attributes_result = (
            self.backend.describe_db_cluster_snapshot_attributes(
                db_cluster_snapshot_identifier=db_cluster_snapshot_identifier,
            )
        )
        template = self.response_template(
            DESCRIBE_DB_CLUSTER_SNAPSHOT_ATTRIBUTES_TEMPLATE
        )
        return template.render(
            db_cluster_snapshot_attributes_result=db_cluster_snapshot_attributes_result,
            db_cluster_snapshot_identifier=db_cluster_snapshot_identifier,
        )

    def modify_db_cluster_snapshot_attribute(self) -> str:
        params = self._get_params()
        db_cluster_snapshot_identifier = params["DBClusterSnapshotIdentifier"]
        db_cluster_snapshot_attributes_result = (
            self.backend.modify_db_cluster_snapshot_attribute(
                db_cluster_snapshot_identifier=db_cluster_snapshot_identifier,
                attribute_name=params["AttributeName"],
                values_to_add=params.get("ValuesToAdd"),
                values_to_remove=params.get("ValuesToRemove"),
            )
        )
        template = self.response_template(MODIFY_DB_CLUSTER_SNAPSHOT_ATTRIBUTE_TEMPLATE)
        return template.render(
            db_cluster_snapshot_attributes_result=db_cluster_snapshot_attributes_result,
            db_cluster_snapshot_identifier=db_cluster_snapshot_identifier,
        )

    def describe_db_proxies(self) -> str:
        params = self._get_params()
        db_proxy_name = params.get("DBProxyName")
        # filters = params.get("Filters")
        marker = params.get("Marker")
        db_proxies = self.backend.describe_db_proxies(
            db_proxy_name=db_proxy_name,
            # filters=filters,
        )
        template = self.response_template(DESCRIBE_DB_PROXIES_TEMPLATE)
        rendered = template.render(dbproxies=db_proxies, marker=marker)
        return rendered

    def create_db_proxy(self) -> str:
        params = self._get_params()
        db_proxy_name = params["DBProxyName"]
        engine_family = params["EngineFamily"]
        auth = params["Auth"]
        role_arn = params["RoleArn"]
        vpc_subnet_ids = params["VpcSubnetIds"]
        vpc_security_group_ids = params.get("VpcSecurityGroupIds")
        require_tls = params.get("RequireTLS")
        idle_client_timeout = params.get("IdleClientTimeout")
        debug_logging = params.get("DebugLogging")
        tags = self.unpack_list_params("Tags", "Tag")
        db_proxy = self.backend.create_db_proxy(
            db_proxy_name=db_proxy_name,
            engine_family=engine_family,
            auth=auth,
            role_arn=role_arn,
            vpc_subnet_ids=vpc_subnet_ids,
            vpc_security_group_ids=vpc_security_group_ids,
            require_tls=require_tls,
            idle_client_timeout=idle_client_timeout,
            debug_logging=debug_logging,
            tags=tags,
        )
        template = self.response_template(CREATE_DB_PROXY_TEMPLATE)
        return template.render(dbproxy=db_proxy)


CREATE_DATABASE_TEMPLATE = """<CreateDBInstanceResponse xmlns="http://rds.amazonaws.com/doc/2014-09-01/">
  <CreateDBInstanceResult>
  {{ database.to_xml() }}
  </CreateDBInstanceResult>
  <ResponseMetadata>
    <RequestId>523e3218-afc7-11c3-90f5-f90431260ab4</RequestId>
  </ResponseMetadata>
</CreateDBInstanceResponse>"""

CREATE_DATABASE_REPLICA_TEMPLATE = """<CreateDBInstanceReadReplicaResponse xmlns="http://rds.amazonaws.com/doc/2014-09-01/">
  <CreateDBInstanceReadReplicaResult>
  {{ database.to_xml() }}
  </CreateDBInstanceReadReplicaResult>
  <ResponseMetadata>
    <RequestId>5e60c46d-a844-11e4-bb68-17f36418e58f</RequestId>
  </ResponseMetadata>
</CreateDBInstanceReadReplicaResponse>"""

DESCRIBE_DATABASES_TEMPLATE = """<DescribeDBInstancesResponse xmlns="http://rds.amazonaws.com/doc/2014-09-01/">
  <DescribeDBInstancesResult>
    <DBInstances>
    {%- for database in databases -%}
      {{ database.to_xml() }}
    {%- endfor -%}
    </DBInstances>
    {% if marker %}
    <Marker>{{ marker }}</Marker>
    {% endif %}
  </DescribeDBInstancesResult>
  <ResponseMetadata>
    <RequestId>523e3218-afc7-11c3-90f5-f90431260ab4</RequestId>
  </ResponseMetadata>
</DescribeDBInstancesResponse>"""

MODIFY_DATABASE_TEMPLATE = """<ModifyDBInstanceResponse xmlns="http://rds.amazonaws.com/doc/2014-09-01/">
  <ModifyDBInstanceResult>
  {{ database.to_xml() }}
  </ModifyDBInstanceResult>
  <ResponseMetadata>
    <RequestId>bb58476c-a1a8-11e4-99cf-55e92d4bbada</RequestId>
  </ResponseMetadata>
</ModifyDBInstanceResponse>"""

PROMOTE_REPLICA_TEMPLATE = """<PromoteReadReplicaResponse xmlns="http://rds.amazonaws.com/doc/2014-09-01/">
  <PromoteReadReplicaResult>
  {{ database.to_xml() }}
  </PromoteReadReplicaResult>
  <ResponseMetadata>
    <RequestId>8e8c0d64-be21-11d3-a71c-13dc2f771e41</RequestId>
  </ResponseMetadata>
</PromoteReadReplicaResponse>"""

REBOOT_DATABASE_TEMPLATE = """<RebootDBInstanceResponse xmlns="http://rds.amazonaws.com/doc/2014-09-01/">
  <RebootDBInstanceResult>
  {{ database.to_xml() }}
  </RebootDBInstanceResult>
  <ResponseMetadata>
    <RequestId>d55711cb-a1ab-11e4-99cf-55e92d4bbada</RequestId>
  </ResponseMetadata>
</RebootDBInstanceResponse>"""

START_DATABASE_TEMPLATE = """<StartDBInstanceResponse xmlns="http://rds.amazonaws.com/doc/2014-10-31/">
  <StartDBInstanceResult>
  {{ database.to_xml() }}
  </StartDBInstanceResult>
  <ResponseMetadata>
    <RequestId>523e3218-afc7-11c3-90f5-f90431260ab9</RequestId>
  </ResponseMetadata>
</StartDBInstanceResponse>"""

STOP_DATABASE_TEMPLATE = """<StopDBInstanceResponse xmlns="http://rds.amazonaws.com/doc/2014-10-31/">
  <StopDBInstanceResult>
  {{ database.to_xml() }}
  </StopDBInstanceResult>
  <ResponseMetadata>
    <RequestId>523e3218-afc7-11c3-90f5-f90431260ab8</RequestId>
  </ResponseMetadata>
</StopDBInstanceResponse>"""

DELETE_DATABASE_TEMPLATE = """<DeleteDBInstanceResponse xmlns="http://rds.amazonaws.com/doc/2014-09-01/">
  <DeleteDBInstanceResult>
    {{ database.to_xml() }}
  </DeleteDBInstanceResult>
  <ResponseMetadata>
    <RequestId>7369556f-b70d-11c3-faca-6ba18376ea1b</RequestId>
  </ResponseMetadata>
</DeleteDBInstanceResponse>"""

DELETE_CLUSTER_TEMPLATE = """<DeleteDBClusterResponse xmlns="http://rds.amazonaws.com/doc/2014-09-01/">
  <DeleteDBClusterResult>
    {{ cluster.to_xml() }}
  </DeleteDBClusterResult>
  <ResponseMetadata>
    <RequestId>7369556f-b70d-11c3-faca-6ba18376ea1b</RequestId>
  </ResponseMetadata>
</DeleteDBClusterResponse>"""

RESTORE_INSTANCE_FROM_SNAPSHOT_TEMPLATE = """<RestoreDBInstanceFromDBSnapshotResponse xmlns="http://rds.amazonaws.com/doc/2014-09-01/">
  <RestoreDBInstanceFromDBSnapshotResult>
  {{ database.to_xml() }}
  </RestoreDBInstanceFromDBSnapshotResult>
  <ResponseMetadata>
    <RequestId>523e3218-afc7-11c3-90f5-f90431260ab4</RequestId>
  </ResponseMetadata>
</RestoreDBInstanceFromDBSnapshotResponse>"""


RESTORE_INSTANCE_TO_POINT_IN_TIME_TEMPLATE = """<RestoreDBInstanceToPointInTimeResponse xmlns="http://rds.amazonaws.com/doc/2014-09-01/">
  <RestoreDBInstanceToPointInTimeResult>
  {{ database.to_xml() }}
  </RestoreDBInstanceToPointInTimeResult>
  <ResponseMetadata>
    <RequestId>523e3218-afc7-11c3-90f5-f90431260ab4</RequestId>
  </ResponseMetadata>
</RestoreDBInstanceToPointInTimeResponse>"""

CREATE_SNAPSHOT_TEMPLATE = """<CreateDBSnapshotResponse xmlns="http://rds.amazonaws.com/doc/2014-09-01/">
  <CreateDBSnapshotResult>
  {{ snapshot.to_xml() }}
  </CreateDBSnapshotResult>
  <ResponseMetadata>
    <RequestId>523e3218-afc7-11c3-90f5-f90431260ab4</RequestId>
  </ResponseMetadata>
</CreateDBSnapshotResponse>
"""

COPY_SNAPSHOT_TEMPLATE = """<CopyDBSnapshotResponse xmlns="http://rds.amazonaws.com/doc/2014-09-01/">
  <CopyDBSnapshotResult>
  {{ snapshot.to_xml() }}
  </CopyDBSnapshotResult>
  <ResponseMetadata>
    <RequestId>523e3218-afc7-11c3-90f5-f90431260ab4</RequestId>
  </ResponseMetadata>
</CopyDBSnapshotResponse>
"""

DESCRIBE_SNAPSHOTS_TEMPLATE = """<DescribeDBSnapshotsResponse xmlns="http://rds.amazonaws.com/doc/2014-09-01/">
  <DescribeDBSnapshotsResult>
    <DBSnapshots>
    {%- for snapshot in snapshots -%}
      {{ snapshot.to_xml() }}
    {%- endfor -%}
    </DBSnapshots>
    {% if marker %}
    <Marker>{{ marker }}</Marker>
    {% endif %}
  </DescribeDBSnapshotsResult>
  <ResponseMetadata>
    <RequestId>523e3218-afc7-11c3-90f5-f90431260ab4</RequestId>
  </ResponseMetadata>
</DescribeDBSnapshotsResponse>"""

DELETE_SNAPSHOT_TEMPLATE = """<DeleteDBSnapshotResponse xmlns="http://rds.amazonaws.com/doc/2014-09-01/">
  <DeleteDBSnapshotResult>
  {{ snapshot.to_xml() }}
  </DeleteDBSnapshotResult>
  <ResponseMetadata>
    <RequestId>523e3218-afc7-11c3-90f5-f90431260ab4</RequestId>
  </ResponseMetadata>
</DeleteDBSnapshotResponse>
"""

CREATE_SECURITY_GROUP_TEMPLATE = """<CreateDBSecurityGroupResponse xmlns="http://rds.amazonaws.com/doc/2014-09-01/">
  <CreateDBSecurityGroupResult>
  {{ security_group.to_xml() }}
  </CreateDBSecurityGroupResult>
  <ResponseMetadata>
    <RequestId>462165d0-a77a-11e4-a5fa-75b30c556f97</RequestId>
  </ResponseMetadata>
</CreateDBSecurityGroupResponse>"""

DESCRIBE_SECURITY_GROUPS_TEMPLATE = """<DescribeDBSecurityGroupsResponse xmlns="http://rds.amazonaws.com/doc/2014-09-01/">
  <DescribeDBSecurityGroupsResult>
    <DBSecurityGroups>
    {% for security_group in security_groups %}
      {{ security_group.to_xml() }}
    {% endfor %}
   </DBSecurityGroups>
  </DescribeDBSecurityGroupsResult>
  <ResponseMetadata>
    <RequestId>5df2014e-a779-11e4-bdb0-594def064d0c</RequestId>
  </ResponseMetadata>
</DescribeDBSecurityGroupsResponse>"""

DELETE_SECURITY_GROUP_TEMPLATE = """<DeleteDBSecurityGroupResponse xmlns="http://rds.amazonaws.com/doc/2014-09-01/">
  <ResponseMetadata>
    <RequestId>97e846bd-a77d-11e4-ac58-91351c0f3426</RequestId>
  </ResponseMetadata>
</DeleteDBSecurityGroupResponse>"""

AUTHORIZE_SECURITY_GROUP_TEMPLATE = """<AuthorizeDBSecurityGroupIngressResponse xmlns="http://rds.amazonaws.com/doc/2014-09-01/">
  <AuthorizeDBSecurityGroupIngressResult>
  {{ security_group.to_xml() }}
  </AuthorizeDBSecurityGroupIngressResult>
  <ResponseMetadata>
    <RequestId>75d32fd5-a77e-11e4-8892-b10432f7a87d</RequestId>
  </ResponseMetadata>
</AuthorizeDBSecurityGroupIngressResponse>"""

CREATE_SUBNET_GROUP_TEMPLATE = """<CreateDBSubnetGroupResponse xmlns="http://rds.amazonaws.com/doc/2014-09-01/">
  <CreateDBSubnetGroupResult>
  {{ subnet_group.to_xml() }}
  </CreateDBSubnetGroupResult>
  <ResponseMetadata>
    <RequestId>3a401b3f-bb9e-11d3-f4c6-37db295f7674</RequestId>
  </ResponseMetadata>
</CreateDBSubnetGroupResponse>"""

DESCRIBE_SUBNET_GROUPS_TEMPLATE = """<DescribeDBSubnetGroupsResponse xmlns="http://rds.amazonaws.com/doc/2014-09-01/">
  <DescribeDBSubnetGroupsResult>
    <DBSubnetGroups>
    {% for subnet_group in subnet_groups %}
      {{ subnet_group.to_xml() }}
    {% endfor %}
    </DBSubnetGroups>
  </DescribeDBSubnetGroupsResult>
  <ResponseMetadata>
    <RequestId>b783db3b-b98c-11d3-fbc7-5c0aad74da7c</RequestId>
  </ResponseMetadata>
</DescribeDBSubnetGroupsResponse>"""

MODIFY_SUBNET_GROUPS_TEMPLATE = """<ModifyDBSubnetGroupResponse xmlns="http://rds.amazonaws.com/doc/2014-09-01/">
  <ModifyDBSubnetGroupResult>
    {{ subnet_group.to_xml() }}
  </ModifyDBSubnetGroupResult>
  <ResponseMetadata>
    <RequestId>b783db3b-b98c-11d3-fbc7-5c0aad74da7c</RequestId>
  </ResponseMetadata>
</ModifyDBSubnetGroupResponse>"""

DELETE_SUBNET_GROUP_TEMPLATE = """<DeleteDBSubnetGroupResponse xmlns="http://rds.amazonaws.com/doc/2014-09-01/">
  <ResponseMetadata>
    <RequestId>13785dd5-a7fc-11e4-bb9c-7f371d0859b0</RequestId>
  </ResponseMetadata>
</DeleteDBSubnetGroupResponse>"""

CREATE_OPTION_GROUP_TEMPLATE = """<CreateOptionGroupResponse xmlns="http://rds.amazonaws.com/doc/2014-09-01/">
  <CreateOptionGroupResult>
  {{ option_group.to_xml() }}
  </CreateOptionGroupResult>
  <ResponseMetadata>
    <RequestId>1e38dad4-9f50-11e4-87ea-a31c60ed2e36</RequestId>
  </ResponseMetadata>
</CreateOptionGroupResponse>"""

DELETE_OPTION_GROUP_TEMPLATE = """<DeleteOptionGroupResponse xmlns="http://rds.amazonaws.com/doc/2014-09-01/">
  <ResponseMetadata>
    <RequestId>e2590367-9fa2-11e4-99cf-55e92d41c60e</RequestId>
  </ResponseMetadata>
</DeleteOptionGroupResponse>"""

DESCRIBE_OPTION_GROUP_TEMPLATE = """<DescribeOptionGroupsResponse xmlns="http://rds.amazonaws.com/doc/2014-09-01/">
  <DescribeOptionGroupsResult>
    <OptionGroupsList>
    {%- for option_group in option_groups -%}
      {{ option_group.to_xml() }}
    {%- endfor -%}
    </OptionGroupsList>
  </DescribeOptionGroupsResult>
  <ResponseMetadata>
    <RequestId>4caf445d-9fbc-11e4-87ea-a31c60ed2e36</RequestId>
  </ResponseMetadata>
</DescribeOptionGroupsResponse>"""

DESCRIBE_OPTION_GROUP_OPTIONS_TEMPLATE = """<DescribeOptionGroupOptionsResponse xmlns="http://rds.amazonaws.com/doc/2014-09-01/">
  <DescribeOptionGroupOptionsResult>
    <OptionGroupOptions>
    {%- for option_group_option in option_group_options -%}
      {{ option_group_option.to_xml() }}
    {%- endfor -%}
    </OptionGroupOptions>
  </DescribeOptionGroupOptionsResult>
  <ResponseMetadata>
    <RequestId>457f7bb8-9fbf-11e4-9084-5754f80d5144</RequestId>
  </ResponseMetadata>
</DescribeOptionGroupOptionsResponse>"""

MODIFY_OPTION_GROUP_TEMPLATE = """<ModifyOptionGroupResponse xmlns="http://rds.amazonaws.com/doc/2014-09-01/">
  <ModifyOptionGroupResult>
    {{ option_group.to_xml() }}
  </ModifyOptionGroupResult>
  <ResponseMetadata>
    <RequestId>ce9284a5-a0de-11e4-b984-a11a53e1f328</RequestId>
  </ResponseMetadata>
</ModifyOptionGroupResponse>"""

CREATE_DB_PARAMETER_GROUP_TEMPLATE = """<CreateDBParameterGroupResponse xmlns="http://rds.amazonaws.com/doc/2014-09-01/">
  <CreateDBParameterGroupResult>
    {{ db_parameter_group.to_xml() }}
  </CreateDBParameterGroupResult>
  <ResponseMetadata>
    <RequestId>7805c127-af22-11c3-96ac-6999cc5f7e72</RequestId>
  </ResponseMetadata>
</CreateDBParameterGroupResponse>"""

DESCRIBE_DB_PARAMETER_GROUPS_TEMPLATE = """<DescribeDBParameterGroupsResponse xmlns="http://rds.amazonaws.com/doc/2014-09-01/">
  <DescribeDBParameterGroupsResult>
    <DBParameterGroups>
    {%- for db_parameter_group in db_parameter_groups -%}
      {{ db_parameter_group.to_xml() }}
    {%- endfor -%}
    </DBParameterGroups>
  </DescribeDBParameterGroupsResult>
  <ResponseMetadata>
    <RequestId>b75d527a-b98c-11d3-f272-7cd6cce12cc5</RequestId>
  </ResponseMetadata>
</DescribeDBParameterGroupsResponse>"""

MODIFY_DB_PARAMETER_GROUP_TEMPLATE = """<ModifyDBParameterGroupResponse xmlns="http://rds.amazonaws.com/doc/2014-09-01/">
  <ModifyDBParameterGroupResult>
    <DBParameterGroupName>{{ db_parameter_group.name }}</DBParameterGroupName>
  </ModifyDBParameterGroupResult>
  <ResponseMetadata>
    <RequestId>12d7435e-bba0-11d3-fe11-33d33a9bb7e3</RequestId>
  </ResponseMetadata>
</ModifyDBParameterGroupResponse>"""

DELETE_DB_PARAMETER_GROUP_TEMPLATE = """<DeleteDBParameterGroupResponse xmlns="http://rds.amazonaws.com/doc/2014-09-01/">
  <ResponseMetadata>
    <RequestId>cad6c267-ba25-11d3-fe11-33d33a9bb7e3</RequestId>
  </ResponseMetadata>
</DeleteDBParameterGroupResponse>"""

DESCRIBE_DB_PARAMETERS_TEMPLATE = """<DescribeDBParametersResponse xmlns="http://rds.amazonaws.com/doc/2014-09-01/">
  <DescribeDBParametersResult>
    <Parameters>
      {%- for db_parameter_name, db_parameter in db_parameter_group.parameters.items() -%}
      <Parameter>
        {%- for parameter_name, parameter_value in db_parameter.items() -%}
        <{{ parameter_name }}>{{ parameter_value }}</{{ parameter_name }}>
        {%- endfor -%}
      </Parameter>
      {%- endfor -%}
    </Parameters>
  </DescribeDBParametersResult>
  <ResponseMetadata>
    <RequestId>8c40488f-b9ff-11d3-a15e-7ac49293f4fa</RequestId>
  </ResponseMetadata>
</DescribeDBParametersResponse>
"""

DESCRIBE_DB_CLUSTER_PARAMETERS_TEMPLATE = """<DescribeDBClusterParametersResponse xmlns="http://rds.amazonaws.com/doc/2014-09-01/">
  <DescribeDBClusterParametersResult>
    <Parameters>
      {%- for param in db_parameter_group -%}
      <Parameter>
        {%- for parameter_name, parameter_value in db_parameter.items() -%}
        <{{ parameter_name }}>{{ parameter_value }}</{{ parameter_name }}>
        {%- endfor -%}
      </Parameter>
      {%- endfor -%}
    </Parameters>
  </DescribeDBClusterParametersResult>
  <ResponseMetadata>
    <RequestId>8c40488f-b9ff-11d3-a15e-7ac49293f4fa</RequestId>
  </ResponseMetadata>
</DescribeDBClusterParametersResponse>
"""

LIST_TAGS_FOR_RESOURCE_TEMPLATE = """<ListTagsForResourceResponse xmlns="http://rds.amazonaws.com/doc/2014-10-31/">
  <ListTagsForResourceResult>
    <TagList>
    {%- for tag in tags -%}
      <Tag>
        <Key>{{ tag['Key'] }}</Key>
        <Value>{{ tag['Value'] }}</Value>
      </Tag>
    {%- endfor -%}
    </TagList>
  </ListTagsForResourceResult>
  <ResponseMetadata>
    <RequestId>8c21ba39-a598-11e4-b688-194eaf8658fa</RequestId>
  </ResponseMetadata>
</ListTagsForResourceResponse>"""

ADD_TAGS_TO_RESOURCE_TEMPLATE = """<AddTagsToResourceResponse xmlns="http://rds.amazonaws.com/doc/2014-10-31/">
  <ResponseMetadata>
    <RequestId>b194d9ca-a664-11e4-b688-194eaf8658fa</RequestId>
  </ResponseMetadata>
</AddTagsToResourceResponse>"""

REMOVE_TAGS_FROM_RESOURCE_TEMPLATE = """<RemoveTagsFromResourceResponse xmlns="http://rds.amazonaws.com/doc/2014-10-31/">
  <ResponseMetadata>
    <RequestId>b194d9ca-a664-11e4-b688-194eaf8658fa</RequestId>
  </ResponseMetadata>
</RemoveTagsFromResourceResponse>"""

CREATE_DB_CLUSTER_TEMPLATE = """<CreateDBClusterResponse xmlns="http://rds.amazonaws.com/doc/2014-09-01/">
  <CreateDBClusterResult>
  {{ cluster.to_xml() }}
  </CreateDBClusterResult>
  <ResponseMetadata>
    <RequestId>523e3218-afc7-11c3-90f5-f90431260ab4</RequestId>
  </ResponseMetadata>
</CreateDBClusterResponse>"""

MODIFY_DB_CLUSTER_TEMPLATE = """<ModifyDBClusterResponse xmlns="http://rds.amazonaws.com/doc/2014-10-31/">
  <ModifyDBClusterResult>
  {{ cluster.to_xml() }}
  </ModifyDBClusterResult>
  <ResponseMetadata>
    <RequestId>69673d54-e48e-4ba4-9333-c5a6c1e7526a</RequestId>
  </ResponseMetadata>
</ModifyDBClusterResponse>"""

DESCRIBE_CLUSTERS_TEMPLATE = """<DescribeDBClustersResponse xmlns="http://rds.amazonaws.com/doc/2014-09-01/">
  <DescribeDBClustersResult>
    <DBClusters>
    {%- for cluster in clusters -%}
      {{ cluster.to_xml() }}
    {%- endfor -%}
    </DBClusters>
    {% if marker %}
    <Marker>{{ marker }}</Marker>
    {% endif %}
  </DescribeDBClustersResult>
  <ResponseMetadata>
    <RequestId>523e3218-afc7-11c3-90f5-f90431260ab4</RequestId>
  </ResponseMetadata>
</DescribeDBClustersResponse>"""

START_CLUSTER_TEMPLATE = """<StartDBClusterResponse xmlns="http://rds.amazonaws.com/doc/2014-10-31/">
  <StartDBClusterResult>
  {{ cluster.to_xml() }}
  </StartDBClusterResult>
  <ResponseMetadata>
    <RequestId>523e3218-afc7-11c3-90f5-f90431260ab9</RequestId>
  </ResponseMetadata>
</StartDBClusterResponse>"""

STOP_CLUSTER_TEMPLATE = """<StopDBClusterResponse xmlns="http://rds.amazonaws.com/doc/2014-10-31/">
  <StopDBClusterResult>
  {{ cluster.to_xml() }}
  </StopDBClusterResult>
  <ResponseMetadata>
    <RequestId>523e3218-afc7-11c3-90f5-f90431260ab8</RequestId>
  </ResponseMetadata>
</StopDBClusterResponse>"""

RESTORE_CLUSTER_FROM_SNAPSHOT_TEMPLATE = """<RestoreDBClusterFromDBSnapshotResponse xmlns="http://rds.amazonaws.com/doc/2014-09-01/">
  <RestoreDBClusterFromSnapshotResult>
  {{ cluster.to_xml() }}
  </RestoreDBClusterFromSnapshotResult>
  <ResponseMetadata>
    <RequestId>523e3218-afc7-11c3-90f5-f90431260ab4</RequestId>
  </ResponseMetadata>
</RestoreDBClusterFromDBSnapshotResponse>
"""

CREATE_CLUSTER_SNAPSHOT_TEMPLATE = """<CreateDBClusterSnapshotResponse xmlns="http://rds.amazonaws.com/doc/2014-09-01/">
  <CreateDBClusterSnapshotResult>
  {{ snapshot.to_xml() }}
  </CreateDBClusterSnapshotResult>
  <ResponseMetadata>
    <RequestId>523e3218-afc7-11c3-90f5-f90431260ab4</RequestId>
  </ResponseMetadata>
</CreateDBClusterSnapshotResponse>
"""

COPY_CLUSTER_SNAPSHOT_TEMPLATE = """<CopyDBClusterSnapshotResponse xmlns="http://rds.amazonaws.com/doc/2014-09-01/">
  <CopyDBClusterSnapshotResult>
  {{ snapshot.to_xml() }}
  </CopyDBClusterSnapshotResult>
  <ResponseMetadata>
    <RequestId>523e3218-afc7-11c3-90f5-f90431260ab4</RequestId>
  </ResponseMetadata>
</CopyDBClusterSnapshotResponse>
"""

DESCRIBE_CLUSTER_SNAPSHOTS_TEMPLATE = """<DescribeDBClusterSnapshotsResponse xmlns="http://rds.amazonaws.com/doc/2014-09-01/">
  <DescribeDBClusterSnapshotsResult>
    <DBClusterSnapshots>
    {%- for snapshot in snapshots -%}
      {{ snapshot.to_xml() }}
    {%- endfor -%}
    </DBClusterSnapshots>
    {% if marker %}
    <Marker>{{ marker }}</Marker>
    {% endif %}
  </DescribeDBClusterSnapshotsResult>
  <ResponseMetadata>
    <RequestId>523e3218-afc7-11c3-90f5-f90431260ab4</RequestId>
  </ResponseMetadata>
</DescribeDBClusterSnapshotsResponse>"""

DELETE_CLUSTER_SNAPSHOT_TEMPLATE = """<DeleteDBClusterSnapshotResponse xmlns="http://rds.amazonaws.com/doc/2014-09-01/">
  <DeleteDBClusterSnapshotResult>
  {{ snapshot.to_xml() }}
  </DeleteDBClusterSnapshotResult>
  <ResponseMetadata>
    <RequestId>523e3218-afc7-11c3-90f5-f90431260ab4</RequestId>
  </ResponseMetadata>
</DeleteDBClusterSnapshotResponse>
"""

START_EXPORT_TASK_TEMPLATE = """<StartExportTaskResponse xmlns="http://rds.amazonaws.com/doc/2014-09-01/">
  <StartExportTaskResult>
  {{ task.to_xml() }}
  </StartExportTaskResult>
  <ResponseMetadata>
    <RequestId>523e3218-afc7-11c3-90f5-f90431260ab4</RequestId>
  </ResponseMetadata>
</StartExportTaskResponse>
"""

CANCEL_EXPORT_TASK_TEMPLATE = """<CancelExportTaskResponse xmlns="http://rds.amazonaws.com/doc/2014-09-01/">
  <CancelExportTaskResult>
  {{ task.to_xml() }}
  </CancelExportTaskResult>
  <ResponseMetadata>
    <RequestId>523e3218-afc7-11c3-90f5-f90431260ab4</RequestId>
  </ResponseMetadata>
</CancelExportTaskResponse>
"""

DESCRIBE_EXPORT_TASKS_TEMPLATE = """<DescribeExportTasksResponse xmlns="http://rds.amazonaws.com/doc/2014-09-01/">
  <DescribeExportTasksResult>
    <ExportTasks>
    {%- for task in tasks -%}
      <ExportTask>{{ task.to_xml() }}</ExportTask>
    {%- endfor -%}
    </ExportTasks>
    {% if marker %}
    <Marker>{{ marker }}</Marker>
    {% endif %}
  </DescribeExportTasksResult>
  <ResponseMetadata>
    <RequestId>523e3218-afc7-11c3-90f5-f90431260ab4</RequestId>
  </ResponseMetadata>
</DescribeExportTasksResponse>
"""

CREATE_EVENT_SUBSCRIPTION_TEMPLATE = """<CreateEventSubscriptionResponse xmlns="http://rds.amazonaws.com/doc/2014-09-01/">
  <CreateEventSubscriptionResult>
  {{ subscription.to_xml() }}
  </CreateEventSubscriptionResult>
  <ResponseMetadata>
    <RequestId>523e3218-afc7-11c3-90f5-f90431260ab4</RequestId>
  </ResponseMetadata>
</CreateEventSubscriptionResponse>
"""

DELETE_EVENT_SUBSCRIPTION_TEMPLATE = """<DeleteEventSubscriptionResponse xmlns="http://rds.amazonaws.com/doc/2014-09-01/">
  <DeleteEventSubscriptionResult>
  {{ subscription.to_xml() }}
  </DeleteEventSubscriptionResult>
  <ResponseMetadata>
    <RequestId>523e3218-afc7-11c3-90f5-f90431260ab4</RequestId>
  </ResponseMetadata>
</DeleteEventSubscriptionResponse>
"""

DESCRIBE_EVENT_SUBSCRIPTIONS_TEMPLATE = """<DescribeEventSubscriptionsResponse xmlns="http://rds.amazonaws.com/doc/2014-09-01/">
  <DescribeEventSubscriptionsResult>
    <EventSubscriptionsList>
      {%- for subscription in subscriptions -%}
        {{ subscription.to_xml() }}
      {%- endfor -%}
    </EventSubscriptionsList>
  </DescribeEventSubscriptionsResult>
  <ResponseMetadata>
    <RequestId>523e3218-afc7-11c3-90f5-f90431260ab4</RequestId>
  </ResponseMetadata>
</DescribeEventSubscriptionsResponse>
"""


DESCRIBE_ORDERABLE_CLUSTER_OPTIONS = """<DescribeOrderableDBInstanceOptionsResponse xmlns="http://rds.amazonaws.com/doc/2014-10-31/">
  <DescribeOrderableDBInstanceOptionsResult>
    <OrderableDBInstanceOptions>
    {% for option in options %}
      <OrderableDBInstanceOption>
        <OutpostCapable>false</OutpostCapable>
        <AvailabilityZones>
          {% for zone in option["AvailabilityZones"] %}
          <AvailabilityZone>
            <Name>{{ zone["Name"] }}</Name>
          </AvailabilityZone>
          {% endfor %}
        </AvailabilityZones>
        <SupportsStorageThroughput>{{ option["SupportsStorageThroughput"] }}</SupportsStorageThroughput>
        <SupportedEngineModes>
          <member>provisioned</member>
        </SupportedEngineModes>
        <SupportsGlobalDatabases>{{ option["SupportsGlobalDatabases"] }}</SupportsGlobalDatabases>
        <SupportsClusters>{{ option["SupportsClusters"] }}</SupportsClusters>
        <Engine>{{ option["Engine"] }}</Engine>
        <SupportedActivityStreamModes/>
        <SupportsEnhancedMonitoring>false</SupportsEnhancedMonitoring>
        <EngineVersion>{{ option["EngineVersion"] }}</EngineVersion>
        <ReadReplicaCapable>false</ReadReplicaCapable>
        <Vpc>true</Vpc>
        <DBInstanceClass>{{ option["DBInstanceClass"] }}</DBInstanceClass>
        <SupportsStorageEncryption>{{ option["SupportsStorageEncryption"] }}</SupportsStorageEncryption>
        <SupportsKerberosAuthentication>{{ option["SupportsKerberosAuthentication"] }}</SupportsKerberosAuthentication>
        <SupportedNetworkTypes>
          <member>IPV4</member>
        </SupportedNetworkTypes>
        <AvailableProcessorFeatures/>
        <SupportsPerformanceInsights>{{ option["SupportsPerformanceInsights"] }}</SupportsPerformanceInsights>
        <LicenseModel>{{ option["LicenseModel"] }}</LicenseModel>
        <MultiAZCapable>{{ option["MultiAZCapable"] }}</MultiAZCapable>
        <RequiresCustomProcessorFeatures>{{ option["RequiresCustomProcessorFeatures"] }}</RequiresCustomProcessorFeatures>
        <StorageType>{{ option["StorageType"] }}</StorageType>
        <SupportsIops>{{ option["SupportsIops"] }}</SupportsIops>
        <SupportsIAMDatabaseAuthentication>{{ option["SupportsIAMDatabaseAuthentication"] }}</SupportsIAMDatabaseAuthentication>
      </OrderableDBInstanceOption>
      {% endfor %}
    </OrderableDBInstanceOptions>
    {% if marker %}
    <Marker>{{ marker }}</Marker>
    {% endif %}
  </DescribeOrderableDBInstanceOptionsResult>
  <ResponseMetadata>
    <RequestId>54212dc5-16c4-4eb8-a88e-448691e877ab</RequestId>
  </ResponseMetadata>
</DescribeOrderableDBInstanceOptionsResponse>"""

CREATE_DB_CLUSTER_PARAMETER_GROUP_TEMPLATE = """<CreateDBClusterParameterGroupResponse xmlns="http://rds.amazonaws.com/doc/2014-09-01/">
  <CreateDBClusterParameterGroupResult>
    {{ db_cluster_parameter_group.to_xml() }}
  </CreateDBClusterParameterGroupResult>
  <ResponseMetadata>
    <RequestId>7805c127-af22-11c3-96ac-6999cc5f7e72</RequestId>
  </ResponseMetadata>
</CreateDBClusterParameterGroupResponse>"""


DESCRIBE_DB_CLUSTER_PARAMETER_GROUPS_TEMPLATE = """<DescribeDBClusterParameterGroupsResponse xmlns="http://rds.amazonaws.com/doc/2014-09-01/">
  <DescribeDBClusterParameterGroupsResult>
    <DBClusterParameterGroups>
    {%- for db_parameter_group in db_parameter_groups -%}
      {{ db_parameter_group.to_xml() }}
    {%- endfor -%}
    </DBClusterParameterGroups>
  </DescribeDBClusterParameterGroupsResult>
  <ResponseMetadata>
    <RequestId>b75d527a-b98c-11d3-f272-7cd6cce12cc5</RequestId>
  </ResponseMetadata>
</DescribeDBClusterParameterGroupsResponse>"""

DELETE_DB_CLUSTER_PARAMETER_GROUP_TEMPLATE = """<DeleteDBClusterParameterGroupResponse xmlns="http://rds.amazonaws.com/doc/2014-09-01/">
  <ResponseMetadata>
    <RequestId>cad6c267-ba25-11d3-fe11-33d33a9bb7e3</RequestId>
  </ResponseMetadata>
</DeleteDBClusterParameterGroupResponse>"""

PROMOTE_READ_REPLICA_DB_CLUSTER_TEMPLATE = """<PromoteReadReplicaDBClusterResponse xmlns="http://rds.amazonaws.com/doc/2014-09-01/">
  <PromoteReadReplicaDBClusterResult>
    {{ cluster.to_xml() }}
  </PromoteReadReplicaDBClusterResult>
  <ResponseMetadata>
    <RequestId>7369556f-b70d-11c3-faca-6ba18376ea1b</RequestId>
  </ResponseMetadata>
</PromoteReadReplicaDBClusterResponse>"""

DESCRIBE_DB_SNAPSHOT_ATTRIBUTES_TEMPLATE = """<DescribeDBSnapshotAttributesResponse xmlns="http://rds.amazonaws.com/doc/2014-10-31/">
  <DescribeDBSnapshotAttributesResult>
    <DBSnapshotAttributesResult>
      <DBSnapshotAttributes>
        {%- for attribute in db_snapshot_attributes_result -%}
          <DBSnapshotAttribute>
            <AttributeName>{{ attribute["AttributeName"] }}</AttributeName>
            <AttributeValues>
              {%- for value in attribute["AttributeValues"] -%}
                <AttributeValue>{{ value }}</AttributeValue>
              {%- endfor -%}
            </AttributeValues>
          </DBSnapshotAttribute>
        {%- endfor -%}
      </DBSnapshotAttributes>
      <DBSnapshotIdentifier>{{ db_snapshot_identifier }}</DBSnapshotIdentifier>
    </DBSnapshotAttributesResult>
  </DescribeDBSnapshotAttributesResult>
  <ResponseMetadata>
    <RequestId>1549581b-12b7-11e3-895e-1334a</RequestId>
  </ResponseMetadata>
</DescribeDBSnapshotAttributesResponse>"""

MODIFY_DB_SNAPSHOT_ATTRIBUTE_TEMPLATE = """<ModifyDBSnapshotAttributeResponse xmlns="http://rds.amazonaws.com/doc/2014-10-31/">
  <ModifyDBSnapshotAttributeResult>
    <DBSnapshotAttributesResult>
      <DBSnapshotAttributes>
        {%- for attribute in db_snapshot_attributes_result -%}
          <DBSnapshotAttribute>
            <AttributeName>{{ attribute["AttributeName"] }}</AttributeName>
            <AttributeValues>
              {%- for value in attribute["AttributeValues"] -%}
                <AttributeValue>{{ value }}</AttributeValue>
              {%- endfor -%}
            </AttributeValues>
          </DBSnapshotAttribute>
        {%- endfor -%}
      </DBSnapshotAttributes>
      <DBSnapshotIdentifier>{{ db_snapshot_identifier }}</DBSnapshotIdentifier>
    </DBSnapshotAttributesResult>
  </ModifyDBSnapshotAttributeResult>
  <ResponseMetadata>
    <RequestId>1549581b-12b7-11e3-895e-1334aEXAMPLE</RequestId>
  </ResponseMetadata>
</ModifyDBSnapshotAttributeResponse>"""

MODIFY_DB_CLUSTER_SNAPSHOT_ATTRIBUTE_TEMPLATE = """<ModifyDBClusterSnapshotAttributeResponse xmlns="http://rds.amazonaws.com/doc/2014-10-31/">
  <ModifyDBClusterSnapshotAttributeResult>
    <DBClusterSnapshotAttributesResult>
      <DBClusterSnapshotAttributes>
        {%- for attribute in db_cluster_snapshot_attributes_result -%}
          <DBClusterSnapshotAttribute>
            <AttributeName>{{ attribute["AttributeName"] }}</AttributeName>
            <AttributeValues>
              {%- for value in attribute["AttributeValues"] -%}
                <AttributeValue>{{ value }}</AttributeValue>
              {%- endfor -%}
            </AttributeValues>
          </DBClusterSnapshotAttribute>
        {%- endfor -%}
      </DBClusterSnapshotAttributes>
      <DBClusterSnapshotIdentifier>{{ db_cluster_snapshot_identifier }}</DBClusterSnapshotIdentifier>
    </DBClusterSnapshotAttributesResult>
  </ModifyDBClusterSnapshotAttributeResult>
  <ResponseMetadata>
    <RequestId>1549581b-12b7-11e3-895e-1334a</RequestId>
  </ResponseMetadata>
</ModifyDBClusterSnapshotAttributeResponse>"""

DESCRIBE_DB_CLUSTER_SNAPSHOT_ATTRIBUTES_TEMPLATE = """<DescribeDBClusterSnapshotAttributesResponse xmlns="http://rds.amazonaws.com/doc/2014-10-31/">
  <DescribeDBClusterSnapshotAttributesResult>
    <DBClusterSnapshotAttributesResult>
      <DBClusterSnapshotAttributes>
        {%- for attribute in db_cluster_snapshot_attributes_result -%}
          <DBClusterSnapshotAttribute>
            <AttributeName>{{ attribute["AttributeName"] }}</AttributeName>
            <AttributeValues>
              {%- for value in attribute["AttributeValues"] -%}
                <AttributeValue>{{ value }}</AttributeValue>
              {%- endfor -%}
            </AttributeValues>
          </DBClusterSnapshotAttribute> 
        {%- endfor -%}
      </DBClusterSnapshotAttributes>
      <DBClusterSnapshotIdentifier>{{ db_cluster_snapshot_identifier }}</DBClusterSnapshotIdentifier>
    </DBClusterSnapshotAttributesResult>
  </DescribeDBClusterSnapshotAttributesResult>
  <ResponseMetadata>
    <RequestId>1549581b-12b7-11e3-895e-1334a</RequestId>
  </ResponseMetadata>
</DescribeDBClusterSnapshotAttributesResponse>"""

CREATE_DB_PROXY_TEMPLATE = """<CreateDBProxyResponse xmlns="http://rds.amazonaws.com/doc/2014-10-31/">
  <CreateDBProxyResult>
    <DBProxy>
    {{ dbproxy.to_xml() }}
    </DBProxy>
  </CreateDBProxyResult>
  <ResponseMetadata>
    <RequestId>1549581b-12b7-11e3-895e-1334aEXAMPLE</RequestId>
  </ResponseMetadata>
</CreateDBProxyResponse>"""

DESCRIBE_DB_PROXIES_TEMPLATE = """<DescribeDBProxiesResponse xmlns="http://rds.amazonaws.com/doc/2014-10-31/">
    <DescribeDBProxiesResult>
        <DBProxies>
          {% for dbproxy in dbproxies %}
            <member>
              {{ dbproxy.to_xml() }}
            </member>
          {% endfor %}
        </DBProxies>
    </DescribeDBProxiesResult>
    <ResponseMetadata>
        <RequestId>1549581b-12b7-11e3-895e-1334a</RequestId>
    </ResponseMetadata>
</DescribeDBProxiesResponse>
"""
