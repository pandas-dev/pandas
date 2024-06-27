from moto.core.responses import BaseResponse

from .models import NeptuneBackend, neptune_backends


class NeptuneResponse(BaseResponse):
    """Handler for Neptune requests and responses."""

    def __init__(self) -> None:
        super().__init__(service_name="neptune")

    @property
    def neptune_backend(self) -> NeptuneBackend:
        """Return backend instance specific for this region."""
        return neptune_backends[self.current_account][self.region]

    @property
    def global_backend(self) -> NeptuneBackend:
        """Return backend instance of the region that stores Global Clusters"""
        return neptune_backends[self.current_account]["us-east-1"]

    def create_db_cluster(self) -> str:
        params = self._get_params()
        availability_zones = params.get("AvailabilityZones")
        backup_retention_period = params.get("BackupRetentionPeriod")
        character_set_name = params.get("CharacterSetName")
        copy_tags_to_snapshot = params.get("CopyTagsToSnapshot")
        database_name = params.get("DatabaseName")
        db_cluster_identifier = params.get("DBClusterIdentifier")
        db_cluster_parameter_group_name = params.get("DBClusterParameterGroupName")
        vpc_security_group_ids = params.get("VpcSecurityGroupIds")
        db_subnet_group_name = params.get("DBSubnetGroupName")
        engine = params.get("Engine")
        engine_version = params.get("EngineVersion")
        port = params.get("Port")
        master_username = params.get("MasterUsername")
        master_user_password = params.get("MasterUserPassword")
        option_group_name = params.get("OptionGroupName")
        preferred_backup_window = params.get("PreferredBackupWindow")
        preferred_maintenance_window = params.get("PreferredMaintenanceWindow")
        replication_source_identifier = params.get("ReplicationSourceIdentifier")
        tags = (self._get_multi_param_dict("Tags") or {}).get("Tag", [])
        storage_encrypted = params.get("StorageEncrypted", "")
        kms_key_id = params.get("KmsKeyId")
        pre_signed_url = params.get("PreSignedUrl")
        enable_iam_database_authentication = params.get(
            "EnableIAMDatabaseAuthentication"
        )
        enable_cloudwatch_logs_exports = params.get("EnableCloudwatchLogsExports")
        deletion_protection = params.get("DeletionProtection")
        serverless_v2_scaling_configuration = params.get(
            "ServerlessV2ScalingConfiguration"
        )
        global_cluster_identifier = params.get("GlobalClusterIdentifier")
        source_region = params.get("SourceRegion")
        db_cluster = self.neptune_backend.create_db_cluster(
            availability_zones=availability_zones,
            backup_retention_period=backup_retention_period,
            character_set_name=character_set_name,
            copy_tags_to_snapshot=copy_tags_to_snapshot,
            database_name=database_name,
            db_cluster_identifier=db_cluster_identifier,
            db_cluster_parameter_group_name=db_cluster_parameter_group_name,
            vpc_security_group_ids=vpc_security_group_ids,
            db_subnet_group_name=db_subnet_group_name,
            engine=engine,
            engine_version=engine_version,
            port=port,
            master_username=master_username,
            master_user_password=master_user_password,
            option_group_name=option_group_name,
            preferred_backup_window=preferred_backup_window,
            preferred_maintenance_window=preferred_maintenance_window,
            replication_source_identifier=replication_source_identifier,
            tags=tags,
            storage_encrypted=storage_encrypted,
            kms_key_id=kms_key_id,
            pre_signed_url=pre_signed_url,
            enable_iam_database_authentication=enable_iam_database_authentication,
            enable_cloudwatch_logs_exports=enable_cloudwatch_logs_exports,
            deletion_protection=deletion_protection,
            serverless_v2_scaling_configuration=serverless_v2_scaling_configuration,
            global_cluster_identifier=global_cluster_identifier,
            source_region=source_region,
        )
        template = self.response_template(CREATE_DB_CLUSTER_TEMPLATE)
        return template.render(cluster=db_cluster)

    def describe_db_clusters(self) -> str:
        params = self._get_params()
        db_cluster_identifier = params["DBClusterIdentifier"]
        db_clusters = self.neptune_backend.describe_db_clusters(
            db_cluster_identifier=db_cluster_identifier
        )
        template = self.response_template(DESCRIBE_DB_CLUSTERS_TEMPLATE)
        return template.render(db_clusters=db_clusters)

    def describe_global_clusters(self) -> str:
        clusters = self.global_backend.describe_global_clusters()
        template = self.response_template(DESCRIBE_GLOBAL_CLUSTERS_TEMPLATE)
        return template.render(clusters=clusters)

    def create_global_cluster(self) -> str:
        params = self._get_params()
        cluster = self.global_backend.create_global_cluster(
            global_cluster_identifier=params["GlobalClusterIdentifier"],
            engine=params.get("Engine"),
            engine_version=params.get("EngineVersion"),
            storage_encrypted=params.get("StorageEncrypted"),
            deletion_protection=params.get("DeletionProtection"),
        )
        template = self.response_template(CREATE_GLOBAL_CLUSTER_TEMPLATE)
        return template.render(cluster=cluster)

    def delete_global_cluster(self) -> str:
        params = self._get_params()
        cluster = self.global_backend.delete_global_cluster(
            global_cluster_identifier=params["GlobalClusterIdentifier"],
        )
        template = self.response_template(DELETE_GLOBAL_CLUSTER_TEMPLATE)
        return template.render(cluster=cluster)


CREATE_DB_CLUSTER_TEMPLATE = """<CreateDBClusterResponse xmlns="http://rds.amazonaws.com/doc/2014-10-31/">
  <ResponseMetadata>
    <RequestId>1549581b-12b7-11e3-895e-1334aEXAMPLE</RequestId>
  </ResponseMetadata>
  <CreateDBClusterResult>
        {{ cluster.to_xml() }}
  </CreateDBClusterResult>
</CreateDBClusterResponse>"""

DESCRIBE_DB_CLUSTERS_TEMPLATE = """<DescribeDBClustersResponse xmlns="http://rds.amazonaws.com/doc/2014-10-31/">
  <ResponseMetadata>
    <RequestId>1549581b-12b7-11e3-895e-1334aEXAMPLE</RequestId>
  </ResponseMetadata>
  <DescribeDBClustersResult>
    <DBClusters>
{% for cluster in db_clusters %}
        {{ cluster.to_xml() }}
{% endfor %}
    </DBClusters>
  </DescribeDBClustersResult>
</DescribeDBClustersResponse>"""

CREATE_GLOBAL_CLUSTER_TEMPLATE = """<CreateGlobalClusterResponse xmlns="http://rds.amazonaws.com/doc/2014-10-31/">
  <ResponseMetadata>
    <RequestId>1549581b-12b7-11e3-895e-1334aEXAMPLE</RequestId>
  </ResponseMetadata>
  <CreateGlobalClusterResult>
  <GlobalCluster>
    {{ cluster.to_xml() }}
  </GlobalCluster>
  </CreateGlobalClusterResult>
</CreateGlobalClusterResponse>"""

DELETE_GLOBAL_CLUSTER_TEMPLATE = """<DeleteGlobalClusterResponse xmlns="http://rds.amazonaws.com/doc/2014-10-31/">
  <ResponseMetadata>
    <RequestId>1549581b-12b7-11e3-895e-1334aEXAMPLE</RequestId>
  </ResponseMetadata>
  <DeleteGlobalClusterResult>
  <GlobalCluster>
    {{ cluster.to_xml() }}
  </GlobalCluster>
  </DeleteGlobalClusterResult>
</DeleteGlobalClusterResponse>"""

DESCRIBE_GLOBAL_CLUSTERS_TEMPLATE = """<DescribeGlobalClustersResponse xmlns="http://rds.amazonaws.com/doc/2014-10-31/">
  <ResponseMetadata>
    <RequestId>1549581b-12b7-11e3-895e-1334aEXAMPLE</RequestId>
  </ResponseMetadata>
  <DescribeGlobalClustersResult>
    <GlobalClusters>
{% for cluster in clusters %}
    <GlobalClusterMember>
        {{ cluster.to_xml() }}
        </GlobalClusterMember>
{% endfor %}
    </GlobalClusters>
  </DescribeGlobalClustersResult>
</DescribeGlobalClustersResponse>"""

REMOVE_FROM_GLOBAL_CLUSTER_TEMPLATE = """<RemoveFromGlobalClusterResponse xmlns="http://rds.amazonaws.com/doc/2014-10-31/">
  <ResponseMetadata>
    <RequestId>1549581b-12b7-11e3-895e-1334aEXAMPLE</RequestId>
  </ResponseMetadata>
  <RemoveFromGlobalClusterResult>
  {% if cluster %}
  <GlobalCluster>
    {{ cluster.to_xml() }}
  </GlobalCluster>
  {% endif %}
  </RemoveFromGlobalClusterResult>
</RemoveFromGlobalClusterResponse>"""
