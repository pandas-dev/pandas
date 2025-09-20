"""Handles incoming memorydb requests, invokes methods, returns responses."""

import json

from moto.core.responses import BaseResponse

from .models import MemoryDBBackend, memorydb_backends


class MemoryDBResponse(BaseResponse):
    """Handler for MemoryDB requests and responses."""

    def __init__(self) -> None:
        super().__init__(service_name="memorydb")

    @property
    def memorydb_backend(self) -> MemoryDBBackend:
        """Return backend instance specific for this region."""
        return memorydb_backends[self.current_account][self.region]

    def create_cluster(self) -> str:
        params = json.loads(self.body)
        cluster_name = params.get("ClusterName")
        node_type = params.get("NodeType")
        parameter_group_name = params.get("ParameterGroupName")
        description = params.get("Description")
        num_shards = params.get("NumShards")
        num_replicas_per_shard = params.get("NumReplicasPerShard")
        subnet_group_name = params.get("SubnetGroupName")
        security_group_ids = params.get("SecurityGroupIds")
        maintenance_window = params.get("MaintenanceWindow")
        port = params.get("Port")
        sns_topic_arn = params.get("SnsTopicArn")
        tls_enabled = params.get("TLSEnabled")
        kms_key_id = params.get("KmsKeyId")
        snapshot_arns = params.get("SnapshotArns")
        snapshot_name = params.get("SnapshotName")
        snapshot_retention_limit = params.get("SnapshotRetentionLimit")
        tags = params.get("Tags")
        snapshot_window = params.get("SnapshotWindow")
        acl_name = params.get("ACLName")
        engine_version = params.get("EngineVersion")
        auto_minor_version_upgrade = params.get("AutoMinorVersionUpgrade")
        data_tiering = params.get("DataTiering")
        cluster = self.memorydb_backend.create_cluster(
            cluster_name=cluster_name,
            node_type=node_type,
            parameter_group_name=parameter_group_name,
            description=description,
            num_shards=num_shards,
            num_replicas_per_shard=num_replicas_per_shard,
            subnet_group_name=subnet_group_name,
            security_group_ids=security_group_ids,
            maintenance_window=maintenance_window,
            port=port,
            sns_topic_arn=sns_topic_arn,
            tls_enabled=tls_enabled,
            kms_key_id=kms_key_id,
            snapshot_arns=snapshot_arns,
            snapshot_name=snapshot_name,
            snapshot_retention_limit=snapshot_retention_limit,
            tags=tags,
            snapshot_window=snapshot_window,
            acl_name=acl_name,
            engine_version=engine_version,
            auto_minor_version_upgrade=auto_minor_version_upgrade,
            data_tiering=data_tiering,
        )
        return json.dumps(dict(Cluster=cluster.to_dict()))

    def create_subnet_group(self) -> str:
        params = json.loads(self.body)
        subnet_group_name = params.get("SubnetGroupName")
        description = params.get("Description")
        subnet_ids = params.get("SubnetIds")
        tags = params.get("Tags")
        subnet_group = self.memorydb_backend.create_subnet_group(
            subnet_group_name=subnet_group_name,
            description=description,
            subnet_ids=subnet_ids,
            tags=tags,
        )
        return json.dumps(dict(SubnetGroup=subnet_group.to_dict()))

    def create_snapshot(self) -> str:
        params = json.loads(self.body)
        cluster_name = params.get("ClusterName")
        snapshot_name = params.get("SnapshotName")
        kms_key_id = params.get("KmsKeyId")
        tags = params.get("Tags")
        snapshot = self.memorydb_backend.create_snapshot(
            cluster_name=cluster_name,
            snapshot_name=snapshot_name,
            kms_key_id=kms_key_id,
            tags=tags,
        )
        return json.dumps(dict(Snapshot=snapshot.to_dict()))

    def describe_clusters(self) -> str:
        params = json.loads(self.body)
        cluster_name = params.get("ClusterName")
        show_shard_details = params.get("ShowShardDetails")
        clusters = self.memorydb_backend.describe_clusters(
            cluster_name=cluster_name,
        )
        return json.dumps(
            dict(
                Clusters=[
                    cluster.to_desc_dict() if show_shard_details else cluster.to_dict()
                    for cluster in clusters
                ]
            )
        )

    def describe_snapshots(self) -> str:
        params = json.loads(self.body)
        cluster_name = params.get("ClusterName")
        snapshot_name = params.get("SnapshotName")
        source = params.get("Source")
        show_detail = params.get("ShowDetail")
        snapshots = self.memorydb_backend.describe_snapshots(
            cluster_name=cluster_name,
            snapshot_name=snapshot_name,
            source=source,
        )
        return json.dumps(
            dict(
                Snapshots=[
                    snapshot.to_desc_dict() if show_detail else snapshot.to_dict()
                    for snapshot in snapshots
                ]
            )
        )

    def describe_subnet_groups(self) -> str:
        params = json.loads(self.body)
        subnet_group_name = params.get("SubnetGroupName")
        subnet_groups = self.memorydb_backend.describe_subnet_groups(
            subnet_group_name=subnet_group_name,
        )
        return json.dumps(dict(SubnetGroups=[sg.to_dict() for sg in subnet_groups]))

    def list_tags(self) -> str:
        params = json.loads(self.body)
        resource_arn = params.get("ResourceArn")
        tag_list = self.memorydb_backend.list_tags(
            resource_arn=resource_arn,
        )
        return json.dumps(dict(TagList=tag_list))

    def tag_resource(self) -> str:
        params = json.loads(self.body)
        resource_arn = params.get("ResourceArn")
        tags = params.get("Tags")
        tag_list = self.memorydb_backend.tag_resource(
            resource_arn=resource_arn,
            tags=tags,
        )
        return json.dumps(dict(TagList=tag_list))

    def untag_resource(self) -> str:
        params = json.loads(self.body)
        resource_arn = params.get("ResourceArn")
        tag_keys = params.get("TagKeys")
        tag_list = self.memorydb_backend.untag_resource(
            resource_arn=resource_arn,
            tag_keys=tag_keys,
        )
        return json.dumps(dict(TagList=tag_list))

    def update_cluster(self) -> str:
        params = json.loads(self.body)
        cluster_name = params.get("ClusterName")
        description = params.get("Description")
        security_group_ids = params.get("SecurityGroupIds")
        maintenance_window = params.get("MaintenanceWindow")
        sns_topic_arn = params.get("SnsTopicArn")
        sns_topic_status = params.get("SnsTopicStatus")
        parameter_group_name = params.get("ParameterGroupName")
        snapshot_window = params.get("SnapshotWindow")
        snapshot_retention_limit = params.get("SnapshotRetentionLimit")
        node_type = params.get("NodeType")
        engine_version = params.get("EngineVersion")
        replica_configuration = params.get("ReplicaConfiguration")
        shard_configuration = params.get("ShardConfiguration")
        acl_name = params.get("ACLName")
        cluster = self.memorydb_backend.update_cluster(
            cluster_name=cluster_name,
            description=description,
            security_group_ids=security_group_ids,
            maintenance_window=maintenance_window,
            sns_topic_arn=sns_topic_arn,
            sns_topic_status=sns_topic_status,
            parameter_group_name=parameter_group_name,
            snapshot_window=snapshot_window,
            snapshot_retention_limit=snapshot_retention_limit,
            node_type=node_type,
            engine_version=engine_version,
            replica_configuration=replica_configuration,
            shard_configuration=shard_configuration,
            acl_name=acl_name,
        )
        return json.dumps(dict(Cluster=cluster.to_dict()))

    def delete_cluster(self) -> str:
        params = json.loads(self.body)
        cluster_name = params.get("ClusterName")
        final_snapshot_name = params.get("FinalSnapshotName")
        cluster = self.memorydb_backend.delete_cluster(
            cluster_name=cluster_name,
            final_snapshot_name=final_snapshot_name,
        )
        return json.dumps(dict(Cluster=cluster.to_dict()))

    def delete_snapshot(self) -> str:
        params = json.loads(self.body)
        snapshot_name = params.get("SnapshotName")
        snapshot = self.memorydb_backend.delete_snapshot(
            snapshot_name=snapshot_name,
        )
        return json.dumps(dict(Snapshot=snapshot.to_dict()))

    def delete_subnet_group(self) -> str:
        params = json.loads(self.body)
        subnet_group_name = params.get("SubnetGroupName")
        subnet_group = self.memorydb_backend.delete_subnet_group(
            subnet_group_name=subnet_group_name,
        )
        return json.dumps(dict(SubnetGroup=subnet_group.to_dict()))
