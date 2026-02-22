"""MemoryDBBackend class with methods for supported APIs."""

import copy
import random
from datetime import datetime
from typing import Any, Optional

from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import BaseModel
from moto.ec2 import ec2_backends
from moto.utilities.tagging_service import TaggingService

from .exceptions import (
    ClusterAlreadyExistsFault,
    ClusterNotFoundFault,
    InvalidParameterValueException,
    InvalidSubnetError,
    SnapshotAlreadyExistsFault,
    SnapshotNotFoundFault,
    SubnetGroupAlreadyExistsFault,
    SubnetGroupInUseFault,
    SubnetGroupNotFoundFault,
    TagNotFoundFault,
)


class MemoryDBCluster(BaseModel):
    def __init__(
        self,
        cluster_name: str,
        node_type: str,
        parameter_group_name: str,
        description: str,
        num_shards: int,
        num_replicas_per_shard: int,
        subnet_group_name: str,
        vpc_id: str,
        maintenance_window: str,
        port: int,
        sns_topic_arn: str,
        kms_key_id: str,
        snapshot_arns: list[str],
        snapshot_name: str,
        snapshot_retention_limit: int,
        snapshot_window: str,
        acl_name: str,
        engine_version: str,
        region: str,
        account_id: str,
        security_group_ids: list[str],
        auto_minor_version_upgrade: bool,
        data_tiering: bool,
        tls_enabled: bool,
    ):
        self.cluster_name = cluster_name
        self.node_type = node_type
        # Default is set to 'default.memorydb-redis7'.
        self.parameter_group_name = parameter_group_name or "default.memorydb-redis7"
        # Setting it to 'in-sync', other option are 'active' or 'applying'.
        self.parameter_group_status = "in-sync"
        self.description = description
        self.num_shards = num_shards or 1  # Default shards is set to 1
        # Defaults to 1 (i.e. 2 nodes per shard).
        self.num_replicas_per_shard = num_replicas_per_shard or 1
        self.subnet_group_name = subnet_group_name
        self.vpc_id = vpc_id
        self.maintenance_window = maintenance_window or "wed:08:00-wed:09:00"
        self.port = port or 6379  # Default is set to 6379
        self.sns_topic_arn = sns_topic_arn
        self.tls_enabled = tls_enabled if tls_enabled is not None else True
        # Clusters that do not have TLS enabled must use the "open-access" ACL to provide open authentication.
        self.acl_name = "open-access" if tls_enabled is False else acl_name
        self.kms_key_id = kms_key_id
        self.snapshot_arns = snapshot_arns
        self.snapshot_name = snapshot_name
        self.snapshot_retention_limit = snapshot_retention_limit or 0
        self.snapshot_window = snapshot_window or "03:00-04:00"
        self.region = region
        self.engine_version = engine_version
        if engine_version == "7.0":
            self.engine_patch_version = "7.0.7"
        elif engine_version == "6.2":
            self.engine_patch_version = "6.2.6"
        else:
            self.engine_version = "7.1"  # Default is '7.1'.
            self.engine_patch_version = "7.1.1"
        self.auto_minor_version_upgrade = (
            auto_minor_version_upgrade
            if auto_minor_version_upgrade is not None
            else True
        )
        self.data_tiering = "true" if data_tiering else "false"
        # The initial status of the cluster will be set to 'creating'."
        self.status = (
            # Set to 'available', other options are 'creating', 'Updating'.
            "available"
        )
        self.pending_updates: dict[Any, Any] = {}  # TODO
        self.shards = self.get_shard_details()

        self.availability_mode = (
            "SingleAZ" if self.num_replicas_per_shard == 0 else "MultiAZ"
        )
        self.cluster_endpoint = {
            "Address": f"clustercfg.{self.cluster_name}.aoneci.memorydb.{region}.amazonaws.com",
            "Port": self.port,
        }
        self.security_group_ids = security_group_ids or []
        self.security_groups = []
        for sg in self.security_group_ids:
            security_group = {"SecurityGroupId": sg, "Status": "active"}
            self.security_groups.append(security_group)
        self.arn = f"arn:aws:memorydb:{region}:{account_id}:cluster/{self.cluster_name}"
        self.sns_topic_status = "active" if self.sns_topic_arn else ""

    def get_shard_details(self) -> list[dict[str, Any]]:
        shards = []
        for i in range(self.num_shards):
            shard_name = f"{i + 1:04}"
            num_nodes = self.num_replicas_per_shard + 1
            nodes = []
            azs = ["a", "b", "c", "d"]
            for n in range(num_nodes):
                node_name = f"{self.cluster_name}-{shard_name}-{n + 1:03}"
                node = {
                    "Name": node_name,
                    "Status": "available",
                    "AvailabilityZone": f"{self.region}{random.choice(azs)}",
                    "CreateTime": datetime.now().strftime(
                        "%Y-%m-%dT%H:%M:%S.000%f+0000"
                    ),
                    "Endpoint": {
                        "Address": f"{node_name}.{self.cluster_name}.aoneci.memorydb.{self.region}.amazonaws.com",
                        "Port": self.port,
                    },
                }
                nodes.append(node)

            shard = {
                "Name": shard_name,
                # Set to 'available', other options are 'creating', 'modifying' , 'deleting'.
                "Status": "available",
                "Slots": f"0-{str(random.randint(10000, 99999))}",
                "Nodes": nodes,
                "NumberOfNodes": num_nodes,
            }
            shards.append(shard)
        return shards

    def update(
        self,
        description: Optional[str],
        security_group_ids: Optional[list[str]],
        maintenance_window: Optional[str],
        sns_topic_arn: Optional[str],
        sns_topic_status: Optional[str],
        parameter_group_name: Optional[str],
        snapshot_window: Optional[str],
        snapshot_retention_limit: Optional[int],
        node_type: Optional[str],
        engine_version: Optional[str],
        replica_configuration: Optional[dict[str, int]],
        shard_configuration: Optional[dict[str, int]],
        acl_name: Optional[str],
    ) -> None:
        if description is not None:
            self.description = description
        if security_group_ids is not None:
            self.security_group_ids = security_group_ids
        if maintenance_window is not None:
            self.maintenance_window = maintenance_window
        if sns_topic_arn is not None:
            self.sns_topic_arn = sns_topic_arn
        if sns_topic_status is not None:
            self.sns_topic_status = sns_topic_status
        if parameter_group_name is not None:
            self.parameter_group_name = parameter_group_name
        if snapshot_window is not None:
            self.snapshot_window = snapshot_window
        if snapshot_retention_limit is not None:
            self.snapshot_retention_limit = snapshot_retention_limit
        if node_type is not None:
            self.node_type = node_type
        if engine_version is not None:
            self.engine_version = engine_version
        if replica_configuration is not None:
            self.num_replicas_per_shard = replica_configuration["ReplicaCount"]
            self.shards = self.get_shard_details()  # update shards and nodes
        if shard_configuration is not None:
            self.num_shards = shard_configuration["ShardCount"]
            self.shards = self.get_shard_details()  # update shards and nodes
        if acl_name is not None:
            self.acl_name = acl_name

    def to_dict(self) -> dict[str, Any]:
        dct = {
            "Name": self.cluster_name,
            "Description": self.description,
            "Status": self.status,
            "PendingUpdates": self.pending_updates,
            "NumberOfShards": self.num_shards,
            "AvailabilityMode": self.availability_mode,
            "ClusterEndpoint": self.cluster_endpoint,
            "NodeType": self.node_type,
            "EngineVersion": self.engine_version,
            "EnginePatchVersion": self.engine_patch_version,
            "ParameterGroupName": self.parameter_group_name,
            "ParameterGroupStatus": self.parameter_group_status,
            "SecurityGroups": self.security_groups,
            "SubnetGroupName": self.subnet_group_name,
            "KmsKeyId": self.kms_key_id,
            "ARN": self.arn,
            "SnsTopicArn": self.sns_topic_arn,
            "SnsTopicStatus": self.sns_topic_status,
            "MaintenanceWindow": self.maintenance_window,
            "SnapshotWindow": self.snapshot_window,
            "ACLName": self.acl_name,
            "DataTiering": self.data_tiering,
        }
        dct_items = {k: v for k, v in dct.items() if v}
        dct_items["TLSEnabled"] = self.tls_enabled
        dct_items["AutoMinorVersionUpgrade"] = self.auto_minor_version_upgrade
        dct_items["SnapshotRetentionLimit"] = self.snapshot_retention_limit
        return dct_items

    def to_desc_dict(self) -> dict[str, Any]:
        dct = self.to_dict()
        dct["Shards"] = self.shards
        return dct


class MemoryDBSubnetGroup(BaseModel):
    def __init__(
        self,
        region_name: str,
        account_id: str,
        ec2_backend: Any,
        subnet_group_name: str,
        description: str,
        subnet_ids: list[str],
        tags: Optional[list[dict[str, str]]] = None,
    ):
        self.ec2_backend = ec2_backend
        self.subnet_group_name = subnet_group_name
        self.description = description
        self.subnet_ids = subnet_ids
        if not self.subnets:
            raise InvalidSubnetError(subnet_ids)
        self.arn = f"arn:aws:memorydb:{region_name}:{account_id}:subnetgroup/{subnet_group_name}"

    @property
    def subnets(self) -> Any:  # type: ignore[misc]
        return self.ec2_backend.describe_subnets(filters={"subnet-id": self.subnet_ids})

    @property
    def vpc_id(self) -> str:
        return self.subnets[0].vpc_id

    def to_dict(self) -> dict[str, Any]:
        return {
            "Name": self.subnet_group_name,
            "Description": self.description,
            "VpcId": self.vpc_id,
            "Subnets": [
                {
                    "Identifier": subnet.id,
                    "AvailabilityZone": {"Name": subnet.availability_zone},
                }
                for subnet in self.subnets
            ],
            "ARN": self.arn,
        }


class MemoryDBSnapshot(BaseModel):
    def __init__(
        self,
        account_id: str,
        region_name: str,
        cluster: MemoryDBCluster,
        snapshot_name: str,
        kms_key_id: Optional[str],
        tags: Optional[list[dict[str, str]]],
        source: Optional[str],
    ):
        self.cluster = copy.copy(cluster)
        self.cluster_name = self.cluster.cluster_name
        self.snapshot_name = snapshot_name
        self.status = "available"
        self.source = source
        self.kms_key_id = kms_key_id if kms_key_id else cluster.kms_key_id
        self.arn = (
            f"arn:aws:memorydb:{region_name}:{account_id}:snapshot/{snapshot_name}"
        )
        self.vpc_id = self.cluster.vpc_id
        self.shards = []
        for i in self.cluster.shards:
            shard = {
                "Name": i["Name"],
                "Configuration": {
                    "Slots": i["Slots"],
                    "ReplicaCount": self.cluster.num_replicas_per_shard,
                },
                "Size": "11 MB",
                "SnapshotCreationTime": datetime.now().strftime(
                    "%Y-%m-%dT%H:%M:%S.000%f+0000"
                ),
            }
            self.shards.append(shard)

    def to_dict(self) -> dict[str, Any]:
        dct = {
            "Name": self.snapshot_name,
            "Status": self.status,
            "Source": self.source,
            "KmsKeyId": self.kms_key_id,
            "ARN": self.arn,
            "ClusterConfiguration": {
                "Name": self.cluster_name,
                "Description": self.cluster.description,
                "NodeType": self.cluster.node_type,
                "EngineVersion": self.cluster.engine_version,
                "MaintenanceWindow": self.cluster.maintenance_window,
                "TopicArn": self.cluster.sns_topic_arn,
                "Port": self.cluster.port,
                "ParameterGroupName": self.cluster.parameter_group_name,
                "SubnetGroupName": self.cluster.subnet_group_name,
                "VpcId": self.vpc_id,
                "SnapshotRetentionLimit": self.cluster.snapshot_retention_limit,
                "SnapshotWindow": self.cluster.snapshot_window,
                "NumShards": self.cluster.num_shards,
            },
            "DataTiering": self.cluster.data_tiering,
        }
        return {k: v for k, v in dct.items() if v}

    def to_desc_dict(self) -> dict[str, Any]:
        dct = self.to_dict()
        dct["ClusterConfiguration"]["Shards"] = self.shards
        return dct


class MemoryDBBackend(BaseBackend):
    """Implementation of MemoryDB APIs."""

    def __init__(self, region_name: str, account_id: str):
        super().__init__(region_name, account_id)

        self.ec2_backend = ec2_backends[account_id][region_name]
        self.clusters: dict[str, MemoryDBCluster] = {}
        self.subnet_groups: dict[str, MemoryDBSubnetGroup] = {
            "default": MemoryDBSubnetGroup(
                region_name,
                account_id,
                self.ec2_backend,
                "default",
                "Default MemoryDB Subnet Group",
                self.get_default_subnets(),
            )
        }
        self.snapshots: dict[str, MemoryDBSnapshot] = {}
        self.tagger = TaggingService()

    def get_default_subnets(self) -> list[str]:
        default_subnets = self.ec2_backend.describe_subnets(
            filters={"default-for-az": "true"}
        )
        default_subnet_ids = [i.id for i in default_subnets]
        return default_subnet_ids

    def _list_arns(self) -> list[str]:
        cluster_arns = [cluster.arn for cluster in self.clusters.values()]
        snapshot_arns = [snapshot.arn for snapshot in self.snapshots.values()]
        subnet_group_arns = [subnet.arn for subnet in self.subnet_groups.values()]
        return cluster_arns + snapshot_arns + subnet_group_arns

    def create_cluster(
        self,
        cluster_name: str,
        node_type: str,
        parameter_group_name: str,
        description: str,
        subnet_group_name: str,
        security_group_ids: list[str],
        maintenance_window: str,
        port: int,
        sns_topic_arn: str,
        tls_enabled: bool,
        kms_key_id: str,
        snapshot_arns: list[str],
        snapshot_name: str,
        snapshot_retention_limit: int,
        tags: list[dict[str, str]],
        snapshot_window: str,
        acl_name: str,
        engine_version: str,
        auto_minor_version_upgrade: bool,
        data_tiering: bool,
        num_shards: int,
        num_replicas_per_shard: int,
    ) -> MemoryDBCluster:
        if cluster_name in self.clusters:
            raise ClusterAlreadyExistsFault(
                msg="Cluster with specified name already exists."
            )

        subnet_group_name = subnet_group_name or "default"
        subnet_group = self.subnet_groups[subnet_group_name]
        vpc_id = subnet_group.vpc_id
        cluster = MemoryDBCluster(
            cluster_name=cluster_name,
            node_type=node_type,
            parameter_group_name=parameter_group_name,
            description=description,
            num_shards=num_shards,
            num_replicas_per_shard=num_replicas_per_shard,
            subnet_group_name=subnet_group_name,
            vpc_id=vpc_id,
            security_group_ids=security_group_ids,
            maintenance_window=maintenance_window,
            port=port,
            sns_topic_arn=sns_topic_arn,
            tls_enabled=tls_enabled,
            kms_key_id=kms_key_id,
            snapshot_arns=snapshot_arns,
            snapshot_name=snapshot_name,
            snapshot_retention_limit=snapshot_retention_limit,
            snapshot_window=snapshot_window,
            acl_name=acl_name,
            engine_version=engine_version,
            auto_minor_version_upgrade=auto_minor_version_upgrade,
            data_tiering=data_tiering,
            region=self.region_name,
            account_id=self.account_id,
        )
        self.clusters[cluster.cluster_name] = cluster
        self.tag_resource(cluster.arn, tags)
        return cluster

    def create_subnet_group(
        self,
        subnet_group_name: str,
        description: str,
        subnet_ids: list[str],
        tags: Optional[list[dict[str, str]]] = None,
    ) -> MemoryDBSubnetGroup:
        if subnet_group_name in self.subnet_groups:
            raise SubnetGroupAlreadyExistsFault(
                msg=f"Subnet group {subnet_group_name} already exists."
            )
        subnet_group = MemoryDBSubnetGroup(
            self.region_name,
            self.account_id,
            self.ec2_backend,
            subnet_group_name,
            description,
            subnet_ids,
            tags,
        )
        self.subnet_groups[subnet_group_name] = subnet_group
        if tags:
            self.tag_resource(subnet_group.arn, tags)
        return subnet_group

    def create_snapshot(
        self,
        cluster_name: str,
        snapshot_name: str,
        kms_key_id: Optional[str] = None,
        tags: Optional[list[dict[str, str]]] = None,
        source: str = "manual",
    ) -> MemoryDBSnapshot:
        if cluster_name not in self.clusters:
            raise ClusterNotFoundFault(msg=f"Cluster not found: {cluster_name}")
        cluster = self.clusters[cluster_name]
        if snapshot_name in self.snapshots:
            raise SnapshotAlreadyExistsFault(
                msg="Snapshot with specified name already exists."
            )

        snapshot = MemoryDBSnapshot(
            account_id=self.account_id,
            region_name=self.region_name,
            cluster=cluster,
            snapshot_name=snapshot_name,
            kms_key_id=kms_key_id,
            tags=tags,
            source=source,
        )
        self.snapshots[snapshot_name] = snapshot
        if tags:
            self.tag_resource(snapshot.arn, tags)
        return snapshot

    def describe_clusters(
        self, cluster_name: Optional[str] = None
    ) -> list[MemoryDBCluster]:
        if cluster_name:
            if cluster_name in self.clusters:
                cluster = self.clusters[cluster_name]
                return [cluster]
            else:
                raise ClusterNotFoundFault(msg=f"Cluster {cluster_name} not found")
        clusters = list(self.clusters.values())
        return clusters

    def describe_snapshots(
        self,
        cluster_name: Optional[str] = None,
        snapshot_name: Optional[str] = None,
        source: Optional[str] = None,
    ) -> list[MemoryDBSnapshot]:
        sources = ["automated", "manual"] if source is None else [source]

        if cluster_name and snapshot_name:
            for snapshot in list(self.snapshots.values()):
                if (
                    snapshot.cluster_name == cluster_name
                    and snapshot.snapshot_name == snapshot_name
                    and snapshot.source in sources
                ):
                    return [snapshot]
                raise SnapshotNotFoundFault(
                    msg=f"Snapshot with name {snapshot_name} not found"
                )

        if cluster_name:
            snapshots = [
                snapshot
                for snapshot in self.snapshots.values()
                if (snapshot.cluster_name == cluster_name)
                and (snapshot.source in sources)
            ]
            return snapshots

        if snapshot_name:
            snapshots = [
                snapshot
                for snapshot in self.snapshots.values()
                if (snapshot.snapshot_name == snapshot_name)
                and (snapshot.source in sources)
            ]
            if snapshots:
                return snapshots
            raise SnapshotNotFoundFault(
                msg=f"Snapshot with name {snapshot_name} not found"
            )

        snapshots = [
            snapshot
            for snapshot in self.snapshots.values()
            if snapshot.source in sources
        ]
        return snapshots

    def describe_subnet_groups(
        self, subnet_group_name: str
    ) -> list[MemoryDBSubnetGroup]:
        if subnet_group_name:
            if subnet_group_name in self.subnet_groups:
                return [self.subnet_groups[subnet_group_name]]
            raise SubnetGroupNotFoundFault(
                msg=f"Subnet group {subnet_group_name} not found."
            )

        subnet_groups = list(self.subnet_groups.values())
        return subnet_groups

    def list_tags(self, resource_arn: str) -> list[dict[str, str]]:
        if resource_arn not in self._list_arns():
            # Get the resource name from the resource_arn
            resource_name = resource_arn.split("/")[-1]
            if "subnetgroup" in resource_arn:
                raise SubnetGroupNotFoundFault(f"{resource_name} is not present")
            elif "snapshot" in resource_arn:
                raise SnapshotNotFoundFault(f"{resource_name} is not present")
            else:
                raise ClusterNotFoundFault(f"{resource_name} is not present")
        return self.tagger.list_tags_for_resource(arn=resource_arn)["Tags"]

    def tag_resource(
        self, resource_arn: str, tags: list[dict[str, str]]
    ) -> list[dict[str, str]]:
        if resource_arn not in self._list_arns():
            resource_name = resource_arn.split("/")[-1]
            if "subnetgroup" in resource_arn:
                raise SubnetGroupNotFoundFault(f"{resource_name} is not present")
            elif "snapshot" in resource_arn:
                raise SnapshotNotFoundFault(f"{resource_name} is not present")
            else:
                raise ClusterNotFoundFault(f"{resource_name} is not present")
        self.tagger.tag_resource(resource_arn, tags)
        return self.tagger.list_tags_for_resource(arn=resource_arn)["Tags"]

    def untag_resource(
        self, resource_arn: str, tag_keys: list[str]
    ) -> list[dict[str, str]]:
        if resource_arn not in self._list_arns():
            resource_name = resource_arn.split("/")[-1]
            if "subnetgroup" in resource_arn:
                raise SubnetGroupNotFoundFault(f"{resource_name} is not present")
            elif "snapshot" in resource_arn:
                raise SnapshotNotFoundFault(f"{resource_name} is not present")
            else:
                raise ClusterNotFoundFault(f"{resource_name} is not present")
        list_tags = self.list_tags(resource_arn=resource_arn)
        list_keys = [i["Key"] for i in list_tags]
        invalid_keys = [key for key in tag_keys if key not in list_keys]
        if invalid_keys:
            raise TagNotFoundFault(msg=f"These tags are not present : {[invalid_keys]}")
        self.tagger.untag_resource_using_names(resource_arn, tag_keys)
        return self.tagger.list_tags_for_resource(arn=resource_arn)["Tags"]

    def update_cluster(
        self,
        cluster_name: str,
        description: Optional[str],
        security_group_ids: Optional[list[str]],
        maintenance_window: Optional[str],
        sns_topic_arn: Optional[str],
        sns_topic_status: Optional[str],
        parameter_group_name: Optional[str],
        snapshot_window: Optional[str],
        snapshot_retention_limit: Optional[int],
        node_type: Optional[str],
        engine_version: Optional[str],
        replica_configuration: Optional[dict[str, int]],
        shard_configuration: Optional[dict[str, int]],
        acl_name: Optional[str],
    ) -> MemoryDBCluster:
        if cluster_name in self.clusters:
            cluster = self.clusters[cluster_name]
            cluster.update(
                description,
                security_group_ids,
                maintenance_window,
                sns_topic_arn,
                sns_topic_status,
                parameter_group_name,
                snapshot_window,
                snapshot_retention_limit,
                node_type,
                engine_version,
                replica_configuration,
                shard_configuration,
                acl_name,
            )
            return cluster
        raise ClusterNotFoundFault(msg="Cluster not found.")

    def delete_cluster(
        self, cluster_name: str, final_snapshot_name: Optional[str]
    ) -> MemoryDBCluster:
        if cluster_name in self.clusters:
            cluster = self.clusters[cluster_name]
            cluster.status = "deleting"
            if final_snapshot_name is not None:  # create snapshot
                self.create_snapshot(
                    cluster_name=cluster_name,
                    snapshot_name=final_snapshot_name,
                    source="manual",
                )
            return self.clusters.pop(cluster_name)
        raise ClusterNotFoundFault(cluster_name)

    def delete_snapshot(self, snapshot_name: str) -> MemoryDBSnapshot:
        if snapshot_name in self.snapshots:
            snapshot = self.snapshots[snapshot_name]
            snapshot.status = "deleting"
            return self.snapshots.pop(snapshot_name)
        raise SnapshotNotFoundFault(snapshot_name)

    def delete_subnet_group(self, subnet_group_name: str) -> MemoryDBSubnetGroup:
        if subnet_group_name in self.subnet_groups:
            if subnet_group_name == "default":
                raise InvalidParameterValueException(
                    msg="default is reserved and cannot be modified."
                )
            if subnet_group_name in [
                c.subnet_group_name for c in self.clusters.values()
            ]:
                raise SubnetGroupInUseFault(
                    msg=f"Subnet group {subnet_group_name} is currently in use by a cluster."
                )
            return self.subnet_groups.pop(subnet_group_name)
        raise SubnetGroupNotFoundFault(
            msg=f"Subnet group {subnet_group_name} not found."
        )


memorydb_backends = BackendDict(MemoryDBBackend, "memorydb")
