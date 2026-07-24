from moto.core.responses import ActionResult, BaseResponse, EmptyResult

from .models import RedshiftBackend, redshift_backends


class RedshiftResponse(BaseResponse):
    RESPONSE_KEY_PATH_TO_TRANSFORMER = {
        "CreateClusterResult.Cluster.ClusterStatus": lambda _: "creating",
        "RestoreFromClusterSnapshotResult.Cluster.ClusterStatus": lambda _: "creating",
    }

    def __init__(self) -> None:
        super().__init__(service_name="redshift")
        self.automated_parameter_parsing = True

    @property
    def redshift_backend(self) -> RedshiftBackend:
        return redshift_backends[self.current_account][self.region]

    def create_cluster(self) -> ActionResult:
        cluster_kwargs = {
            "cluster_identifier": self._get_param("ClusterIdentifier"),
            "node_type": self._get_param("NodeType"),
            "master_username": self._get_param("MasterUsername"),
            "master_user_password": self._get_param("MasterUserPassword"),
            "db_name": self._get_param("DBName"),
            "cluster_type": self._get_param("ClusterType"),
            "cluster_security_groups": self._get_param("ClusterSecurityGroups", []),
            "vpc_security_group_ids": self._get_param("VpcSecurityGroupIds", []),
            "cluster_subnet_group_name": self._get_param("ClusterSubnetGroupName"),
            "availability_zone": self._get_param("AvailabilityZone"),
            "preferred_maintenance_window": self._get_param(
                "PreferredMaintenanceWindow"
            ),
            "cluster_parameter_group_name": self._get_param(
                "ClusterParameterGroupName"
            ),
            "automated_snapshot_retention_period": self._get_int_param(
                "AutomatedSnapshotRetentionPeriod"
            ),
            "port": self._get_int_param("Port"),
            "cluster_version": self._get_param("ClusterVersion"),
            "allow_version_upgrade": self._get_bool_param("AllowVersionUpgrade", True),
            "number_of_nodes": self._get_int_param("NumberOfNodes"),
            "publicly_accessible": self._get_bool_param("PubliclyAccessible", False),
            "encrypted": self._get_bool_param("Encrypted", False),
            "region_name": self.region,
            "tags": self._get_param("Tags", []),
            "iam_roles_arn": self._get_param("IamRoles", []),
            "enhanced_vpc_routing": self._get_bool_param("EnhancedVpcRouting", False),
            "kms_key_id": self._get_param("KmsKeyId"),
        }
        cluster = self.redshift_backend.create_cluster(**cluster_kwargs)
        return ActionResult({"Cluster": cluster})

    def pause_cluster(self) -> ActionResult:
        cluster_id = self._get_param("ClusterIdentifier")
        cluster = self.redshift_backend.pause_cluster(cluster_id)
        return ActionResult({"Cluster": cluster})

    def resume_cluster(self) -> ActionResult:
        cluster_id = self._get_param("ClusterIdentifier")
        cluster = self.redshift_backend.resume_cluster(cluster_id)
        return ActionResult({"Cluster": cluster})

    def restore_from_cluster_snapshot(self) -> ActionResult:
        enhanced_vpc_routing = self._get_bool_param("EnhancedVpcRouting")
        node_type = self._get_param("NodeType")
        number_of_nodes = self._get_int_param("NumberOfNodes")
        restore_kwargs = {
            "snapshot_identifier": self._get_param("SnapshotIdentifier"),
            "cluster_identifier": self._get_param("ClusterIdentifier"),
            "port": self._get_int_param("Port"),
            "availability_zone": self._get_param("AvailabilityZone"),
            "allow_version_upgrade": self._get_bool_param("AllowVersionUpgrade"),
            "cluster_subnet_group_name": self._get_param("ClusterSubnetGroupName"),
            "publicly_accessible": self._get_param("PubliclyAccessible"),
            "cluster_parameter_group_name": self._get_param(
                "ClusterParameterGroupName"
            ),
            "cluster_security_groups": self._get_param("ClusterSecurityGroups", []),
            "vpc_security_group_ids": self._get_param("VpcSecurityGroupIds", []),
            "preferred_maintenance_window": self._get_param(
                "PreferredMaintenanceWindow"
            ),
            "automated_snapshot_retention_period": self._get_int_param(
                "AutomatedSnapshotRetentionPeriod"
            ),
            "region_name": self.region,
            "iam_roles_arn": self._get_param("IamRoles", []),
        }
        if enhanced_vpc_routing is not None:
            restore_kwargs["enhanced_vpc_routing"] = enhanced_vpc_routing
        if node_type is not None:
            restore_kwargs["node_type"] = node_type
        if number_of_nodes is not None:
            restore_kwargs["number_of_nodes"] = number_of_nodes
        cluster = self.redshift_backend.restore_from_cluster_snapshot(**restore_kwargs)
        return ActionResult({"Cluster": cluster})

    def describe_clusters(self) -> ActionResult:
        cluster_identifier = self._get_param("ClusterIdentifier")
        tag_keys = self._get_param("TagKeys")
        clusters = self.redshift_backend.describe_clusters(cluster_identifier, tag_keys)

        return ActionResult({"Clusters": clusters})

    def modify_cluster(self) -> ActionResult:
        request_kwargs = {
            "cluster_identifier": self._get_param("ClusterIdentifier"),
            "new_cluster_identifier": self._get_param("NewClusterIdentifier"),
            "node_type": self._get_param("NodeType"),
            "master_user_password": self._get_param("MasterUserPassword"),
            "cluster_type": self._get_param("ClusterType"),
            "cluster_security_groups": self._get_param("ClusterSecurityGroups", []),
            "vpc_security_group_ids": self._get_param("VpcSecurityGroupIds", []),
            "cluster_subnet_group_name": self._get_param("ClusterSubnetGroupName"),
            "preferred_maintenance_window": self._get_param(
                "PreferredMaintenanceWindow"
            ),
            "cluster_parameter_group_name": self._get_param(
                "ClusterParameterGroupName"
            ),
            "automated_snapshot_retention_period": self._get_int_param(
                "AutomatedSnapshotRetentionPeriod"
            ),
            "cluster_version": self._get_param("ClusterVersion"),
            "allow_version_upgrade": self._get_bool_param("AllowVersionUpgrade"),
            "number_of_nodes": self._get_int_param("NumberOfNodes"),
            "publicly_accessible": self._get_param("PubliclyAccessible"),
            "encrypted": self._get_param("Encrypted"),
            "iam_roles_arn": self._get_param("IamRoles", []),
            "enhanced_vpc_routing": self._get_param("EnhancedVpcRouting"),
        }
        cluster_kwargs = {}
        # We only want parameters that were actually passed in, otherwise
        # we'll stomp all over our cluster metadata with None values.
        for key, value in request_kwargs.items():
            if value is not None and value != []:
                cluster_kwargs[key] = value

        cluster = self.redshift_backend.modify_cluster(**cluster_kwargs)

        return ActionResult({"Cluster": cluster})

    def delete_cluster(self) -> ActionResult:
        request_kwargs = {
            "cluster_identifier": self._get_param("ClusterIdentifier"),
            "final_cluster_snapshot_identifier": self._get_param(
                "FinalClusterSnapshotIdentifier"
            ),
            "skip_final_snapshot": self._get_bool_param("SkipFinalClusterSnapshot"),
        }

        cluster = self.redshift_backend.delete_cluster(**request_kwargs)

        return ActionResult({"Cluster": cluster})

    def create_cluster_subnet_group(self) -> ActionResult:
        cluster_subnet_group_name = self._get_param("ClusterSubnetGroupName")
        description = self._get_param("Description")
        subnet_ids = self._get_param("SubnetIds", [])
        tags = self._get_param("Tags", [])

        subnet_group = self.redshift_backend.create_cluster_subnet_group(
            cluster_subnet_group_name=cluster_subnet_group_name,
            description=description,
            subnet_ids=subnet_ids,
            region_name=self.region,
            tags=tags,
        )

        return ActionResult({"ClusterSubnetGroup": subnet_group})

    def describe_cluster_subnet_groups(self) -> ActionResult:
        subnet_identifier = self._get_param("ClusterSubnetGroupName")
        subnet_groups = self.redshift_backend.describe_cluster_subnet_groups(
            subnet_identifier
        )

        return ActionResult({"ClusterSubnetGroups": subnet_groups})

    def delete_cluster_subnet_group(self) -> ActionResult:
        subnet_identifier = self._get_param("ClusterSubnetGroupName")
        self.redshift_backend.delete_cluster_subnet_group(subnet_identifier)

        return EmptyResult()

    def create_cluster_security_group(self) -> ActionResult:
        cluster_security_group_name = self._get_param("ClusterSecurityGroupName")
        description = self._get_param("Description")
        tags = self._get_param("Tags", [])

        security_group = self.redshift_backend.create_cluster_security_group(
            cluster_security_group_name=cluster_security_group_name,
            description=description,
            tags=tags,
        )

        return ActionResult({"ClusterSecurityGroup": security_group})

    def describe_cluster_security_groups(self) -> ActionResult:
        cluster_security_group_name = self._get_param("ClusterSecurityGroupName")
        security_groups = self.redshift_backend.describe_cluster_security_groups(
            cluster_security_group_name
        )

        return ActionResult({"ClusterSecurityGroups": security_groups})

    def delete_cluster_security_group(self) -> ActionResult:
        security_group_identifier = self._get_param("ClusterSecurityGroupName")
        self.redshift_backend.delete_cluster_security_group(security_group_identifier)

        return EmptyResult()

    def authorize_cluster_security_group_ingress(self) -> ActionResult:
        cluster_security_group_name = self._get_param("ClusterSecurityGroupName")
        cidr_ip = self._get_param("CIDRIP")

        security_group = self.redshift_backend.authorize_cluster_security_group_ingress(
            cluster_security_group_name, cidr_ip
        )
        result = {
            "ClusterSecurityGroup": {
                "ClusterSecurityGroupName": cluster_security_group_name,
                "Description": security_group.description,
                "IPRanges": [
                    {
                        "Status": "authorized",
                        "CIDRIP": cidr_ip,
                        "Tags": security_group.tags,
                    },
                ],
            }
        }
        return ActionResult(result)

    def create_cluster_parameter_group(self) -> ActionResult:
        cluster_parameter_group_name = self._get_param("ParameterGroupName")
        group_family = self._get_param("ParameterGroupFamily")
        description = self._get_param("Description")
        tags = self._get_param("Tags", [])

        parameter_group = self.redshift_backend.create_cluster_parameter_group(
            cluster_parameter_group_name, group_family, description, self.region, tags
        )

        return ActionResult({"ClusterParameterGroup": parameter_group})

    def describe_cluster_parameter_groups(self) -> ActionResult:
        cluster_parameter_group_name = self._get_param("ParameterGroupName")
        parameter_groups = self.redshift_backend.describe_cluster_parameter_groups(
            cluster_parameter_group_name
        )

        return ActionResult({"ParameterGroups": parameter_groups})

    def delete_cluster_parameter_group(self) -> ActionResult:
        cluster_parameter_group_name = self._get_param("ParameterGroupName")
        self.redshift_backend.delete_cluster_parameter_group(
            cluster_parameter_group_name
        )

        return EmptyResult()

    def describe_default_cluster_parameters(self) -> ActionResult:
        family = self._get_param("ParameterGroupFamily")
        params = self.redshift_backend.describe_default_cluster_parameters()
        return ActionResult(
            {
                "DefaultClusterParameters": {
                    "ParameterGroupFamily": family,
                    "Parameters": params,
                }
            }
        )

    def describe_cluster_parameters(self) -> ActionResult:
        group_name = self._get_param("ParameterGroupName")
        params = self.redshift_backend.describe_cluster_parameters(
            parameter_group_name=group_name
        )
        return ActionResult(
            {
                "Parameters": params,
            }
        )

    def create_cluster_snapshot(self) -> ActionResult:
        cluster_identifier = self._get_param("ClusterIdentifier")
        snapshot_identifier = self._get_param("SnapshotIdentifier")
        tags = self._get_param("Tags", [])

        snapshot = self.redshift_backend.create_cluster_snapshot(
            cluster_identifier, snapshot_identifier, self.region, tags
        )
        return ActionResult({"Snapshot": snapshot})

    def describe_cluster_snapshots(self) -> ActionResult:
        cluster_identifier = self._get_param("ClusterIdentifier")
        snapshot_identifier = self._get_param("SnapshotIdentifier")
        snapshot_type = self._get_param("SnapshotType")
        snapshots = self.redshift_backend.describe_cluster_snapshots(
            cluster_identifier, snapshot_identifier, snapshot_type
        )
        return ActionResult({"Snapshots": snapshots})

    def delete_cluster_snapshot(self) -> ActionResult:
        snapshot_identifier = self._get_param("SnapshotIdentifier")
        snapshot = self.redshift_backend.delete_cluster_snapshot(snapshot_identifier)

        return ActionResult({"Snapshot": snapshot})

    def create_snapshot_copy_grant(self) -> ActionResult:
        copy_grant_kwargs = {
            "snapshot_copy_grant_name": self._get_param("SnapshotCopyGrantName"),
            "kms_key_id": self._get_param("KmsKeyId"),
            "region_name": self._get_param("Region"),
        }

        copy_grant = self.redshift_backend.create_snapshot_copy_grant(
            **copy_grant_kwargs
        )
        return ActionResult({"SnapshotCopyGrant": copy_grant})

    def delete_snapshot_copy_grant(self) -> ActionResult:
        copy_grant_kwargs = {
            "snapshot_copy_grant_name": self._get_param("SnapshotCopyGrantName")
        }
        self.redshift_backend.delete_snapshot_copy_grant(**copy_grant_kwargs)
        return EmptyResult()

    def describe_snapshot_copy_grants(self) -> ActionResult:
        copy_grant_kwargs = {
            "snapshot_copy_grant_name": self._get_param("SnapshotCopyGrantName")
        }

        copy_grants = self.redshift_backend.describe_snapshot_copy_grants(
            **copy_grant_kwargs
        )
        return ActionResult({"SnapshotCopyGrants": copy_grants})

    def create_tags(self) -> ActionResult:
        resource_name = self._get_param("ResourceName")
        tags = self._get_param("Tags", [])

        self.redshift_backend.create_tags(resource_name, tags)

        return EmptyResult()

    def describe_tags(self) -> ActionResult:
        resource_name = self._get_param("ResourceName")
        resource_type = self._get_param("ResourceType")

        tagged_resources = self.redshift_backend.describe_tags(
            resource_name, resource_type
        )
        return ActionResult({"TaggedResources": tagged_resources})

    def delete_tags(self) -> ActionResult:
        resource_name = self._get_param("ResourceName")
        tag_keys = self._get_param("TagKeys", [])

        self.redshift_backend.delete_tags(resource_name, tag_keys)

        return EmptyResult()

    def enable_snapshot_copy(self) -> ActionResult:
        snapshot_copy_kwargs = {
            "cluster_identifier": self._get_param("ClusterIdentifier"),
            "destination_region": self._get_param("DestinationRegion"),
            "retention_period": self._get_param("RetentionPeriod", 7),
            "snapshot_copy_grant_name": self._get_param("SnapshotCopyGrantName"),
        }
        cluster = self.redshift_backend.enable_snapshot_copy(**snapshot_copy_kwargs)

        return ActionResult({"Cluster": cluster})

    def disable_snapshot_copy(self) -> ActionResult:
        snapshot_copy_kwargs = {
            "cluster_identifier": self._get_param("ClusterIdentifier")
        }
        cluster = self.redshift_backend.disable_snapshot_copy(**snapshot_copy_kwargs)

        return ActionResult({"Cluster": cluster})

    def modify_snapshot_copy_retention_period(self) -> ActionResult:
        snapshot_copy_kwargs = {
            "cluster_identifier": self._get_param("ClusterIdentifier"),
            "retention_period": self._get_param("RetentionPeriod"),
        }
        cluster = self.redshift_backend.modify_snapshot_copy_retention_period(
            **snapshot_copy_kwargs
        )

        return ActionResult({"Clusters": [cluster]})

    def get_cluster_credentials(self) -> ActionResult:
        cluster_identifier = self._get_param("ClusterIdentifier")
        db_user = self._get_param("DbUser")
        auto_create = self._get_bool_param("AutoCreate", False)
        duration_seconds = self._get_int_param("DurationSeconds", 900)

        cluster_credentials = self.redshift_backend.get_cluster_credentials(
            cluster_identifier, db_user, auto_create, duration_seconds
        )

        return ActionResult(cluster_credentials)

    def enable_logging(self) -> ActionResult:
        cluster_identifier = self._get_param("ClusterIdentifier")
        bucket_name = self._get_param("BucketName")
        s3_key_prefix = self._get_param("S3KeyPrefix")
        log_destination_type = self._get_param("LogDestinationType")
        log_exports = self._get_param("LogExports")
        config = self.redshift_backend.enable_logging(
            cluster_identifier=cluster_identifier,
            bucket_name=bucket_name,
            s3_key_prefix=s3_key_prefix,
            log_destination_type=log_destination_type,
            log_exports=log_exports,
        )
        return ActionResult(config)

    def disable_logging(self) -> ActionResult:
        cluster_identifier = self._get_param("ClusterIdentifier")
        config = self.redshift_backend.disable_logging(
            cluster_identifier=cluster_identifier,
        )
        return ActionResult(config)

    def describe_logging_status(self) -> ActionResult:
        cluster_identifier = self._get_param("ClusterIdentifier")
        config = self.redshift_backend.describe_logging_status(
            cluster_identifier=cluster_identifier,
        )
        return ActionResult(config)
