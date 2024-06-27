import copy
import string
from typing import Any, Dict, List, Optional

from jinja2 import Template

from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import BaseModel
from moto.core.utils import iso_8601_datetime_with_milliseconds
from moto.moto_api._internal import mock_random as random
from moto.utilities.utils import get_partition, load_resource

from .exceptions import DBClusterNotFoundError


class GlobalCluster(BaseModel):
    def __init__(
        self,
        account_id: str,
        region_name: str,
        global_cluster_identifier: str,
        engine: Optional[str],
        engine_version: Optional[str],
        storage_encrypted: Optional[str],
        deletion_protection: Optional[str],
    ):
        self.global_cluster_identifier = global_cluster_identifier
        self.global_cluster_resource_id = "cluster-" + random.get_random_hex(8)
        self.global_cluster_arn = f"arn:{get_partition(region_name)}:rds::{account_id}:global-cluster:{global_cluster_identifier}"
        self.engine = engine or "neptune"
        self.engine_version = engine_version or "1.2.0.0"
        self.storage_encrypted = (
            storage_encrypted and storage_encrypted.lower() == "true"
        )
        self.deletion_protection = (
            deletion_protection and deletion_protection.lower() == "true"
        )

    def to_xml(self) -> str:
        template = Template(
            """
          <GlobalClusterIdentifier>{{ cluster.global_cluster_identifier }}</GlobalClusterIdentifier>
          <GlobalClusterResourceId>{{ cluster.global_cluster_resource_id }}</GlobalClusterResourceId>
          <GlobalClusterArn>{{ cluster.global_cluster_arn }}</GlobalClusterArn>
          <Engine>{{ cluster.engine }}</Engine>
          <Status>available</Status>
          <EngineVersion>{{ cluster.engine_version }}</EngineVersion>
          <StorageEncrypted>{{ 'true' if cluster.storage_encrypted else 'false' }}</StorageEncrypted>
          <DeletionProtection>{{ 'true' if cluster.deletion_protection else 'false' }}</DeletionProtection>"""
        )
        return template.render(cluster=self)


class DBCluster(BaseModel):
    def __init__(
        self,
        account_id: str,
        region_name: str,
        db_cluster_identifier: str,
        database_name: Optional[str],
        tags: List[Dict[str, str]],
        storage_encrypted: str,
        parameter_group_name: str,
        engine: str,
        engine_version: str,
        kms_key_id: Optional[str],
        preferred_maintenance_window: Optional[str],
        preferred_backup_window: Optional[str],
        backup_retention_period: Optional[int],
        port: Optional[int],
        serverless_v2_scaling_configuration: Optional[Dict[str, int]],
    ):
        self.account_id = account_id
        self.region_name = region_name
        self.db_cluster_identifier = db_cluster_identifier
        self.resource_id = "cluster-" + random.get_random_hex(8)
        self.tags = tags
        self.storage_encrypted = storage_encrypted.lower() != "false"
        self.db_cluster_parameter_group_name = parameter_group_name
        self.engine = engine
        self.engine_version = engine_version
        self.database_name = database_name
        self.db_subnet_group = "default"
        self.status = "available"
        self.backup_retention_period = backup_retention_period
        self.cluster_create_time = iso_8601_datetime_with_milliseconds()
        self.url_identifier = "".join(
            random.choice(string.ascii_lowercase + string.digits) for _ in range(12)
        )
        self.endpoint = f"{self.db_cluster_identifier}.cluster-{self.url_identifier}.{self.region_name}.neptune.amazonaws.com"
        self.reader_endpoint = f"{self.db_cluster_identifier}.cluster-ro-{self.url_identifier}.{self.region_name}.neptune.amazonaws.com"
        self.resource_id = "cluster-" + "".join(
            random.choice(string.ascii_uppercase + string.digits) for _ in range(26)
        )
        self.hosted_zone_id = "".join(
            random.choice(string.ascii_uppercase + string.digits) for _ in range(14)
        )
        self.kms_key_id = kms_key_id or (
            "default_kms_key_id" if self.storage_encrypted else None
        )
        self.preferred_maintenance_window = preferred_maintenance_window
        self.preferred_backup_window = preferred_backup_window
        self.port = port
        self.availability_zones = [
            f"{self.region_name}a",
            f"{self.region_name}b",
            f"{self.region_name}c",
        ]
        self.serverless_v2_scaling_configuration = serverless_v2_scaling_configuration

    @property
    def db_cluster_arn(self) -> str:
        return f"arn:{get_partition(self.region_name)}:rds:{self.region_name}:{self.account_id}:cluster:{self.db_cluster_identifier}"

    def get_tags(self) -> List[Dict[str, str]]:
        return self.tags

    def add_tags(self, tags: List[Dict[str, str]]) -> List[Dict[str, str]]:
        new_keys = [tag_set["Key"] for tag_set in tags]
        self.tags = [tag_set for tag_set in self.tags if tag_set["Key"] not in new_keys]
        self.tags.extend(tags)
        return self.tags

    def remove_tags(self, tag_keys: List[str]) -> None:
        self.tags = [tag_set for tag_set in self.tags if tag_set["Key"] not in tag_keys]

    def to_xml(self) -> str:
        template = Template(
            """<DBCluster>
            {% if cluster.allocated_storage %}
        <AllocatedStorage>{{ cluster.allocated_storage }}</AllocatedStorage>
        {% endif %}
        <AvailabilityZones>
          {% for zone in cluster.availability_zones %}
          <AvailabilityZone>{{ zone }}</AvailabilityZone>
          {% endfor %}
        </AvailabilityZones>
        {% if cluster.backup_retention_period %}
        <BackupRetentionPeriod>{{ cluster.backup_retention_period }}</BackupRetentionPeriod>
        {% endif %}
        {% if cluster.character_set_name %}
        <CharacterSetName>{{ cluster.character_set_name }}</CharacterSetName>
        {% endif %}
        {% if cluster.database_name %}
        <DatabaseName>{{ cluster.database_name }}</DatabaseName>
        {% endif %}
        <DBClusterIdentifier>{{ cluster.db_cluster_identifier }}</DBClusterIdentifier>
        <DBClusterParameterGroup>{{ cluster.db_cluster_parameter_group_name }}</DBClusterParameterGroup>
        <DBSubnetGroup>{{ cluster.db_subnet_group }}</DBSubnetGroup>
        <Status>{{ cluster.status }}</Status>
        <PercentProgress>{{ cluster.percent_progress }}</PercentProgress>
        {% if cluster.earliest_restorable_time %}
        <EarliestRestorableTime>{{ cluster.earliest_restorable_time }}</EarliestRestorableTime>
        {% endif %}
        <Endpoint>{{ cluster.endpoint }}</Endpoint>
        <ReaderEndpoint>{{ cluster.reader_endpoint }}</ReaderEndpoint>
        <MultiAZ>false</MultiAZ>
        <Engine>{{ cluster.engine }}</Engine>
        <EngineVersion>{{ cluster.engine_version }}</EngineVersion>
        {% if cluster.latest_restorable_time %}
        <LatestRestorableTime>{{ cluster.latest_restorable_time }}</LatestRestorableTime>
        {% endif %}
        {% if cluster.port %}
        <Port>{{ cluster.port }}</Port>
        {% endif %}
        <MasterUsername>{{ cluster.master_username }}</MasterUsername>
        <DBClusterOptionGroupMemberships>
{% for dbclusteroptiongroupmembership in cluster.dbclusteroptiongroupmemberships %}
          <member>
            <DBClusterOptionGroupName>{{ dbclusteroptiongroupmembership.db_cluster_option_group_name }}</DBClusterOptionGroupName>
            <Status>{{ dbclusteroptiongroupmembership.status }}</Status>
          </member>
{% endfor %}
        </DBClusterOptionGroupMemberships>
        <PreferredBackupWindow>{{ cluster.preferred_backup_window }}</PreferredBackupWindow>
        <PreferredMaintenanceWindow>{{ cluster.preferred_maintenance_window }}</PreferredMaintenanceWindow>
        <ReplicationSourceIdentifier>{{ cluster.replication_source_identifier }}</ReplicationSourceIdentifier>
        <ReadReplicaIdentifiers>
{% for readreplicaidentifier in cluster.readreplicaidentifiers %}
          <member/>
{% endfor %}
        </ReadReplicaIdentifiers>
        <DBClusterMembers>
{% for dbclustermember in cluster.dbclustermembers %}
          <member>
            <DBInstanceIdentifier>{{ dbclustermember.db_instance_identifier }}</DBInstanceIdentifier>
            <IsClusterWriter>{{ dbclustermember.is_cluster_writer }}</IsClusterWriter>
            <DBClusterParameterGroupStatus>{{ dbclustermember.db_cluster_parameter_group_status }}</DBClusterParameterGroupStatus>
            <PromotionTier>{{ dbclustermember.promotion_tier }}</PromotionTier>
          </member>
{% endfor %}
        </DBClusterMembers>
        <VpcSecurityGroups>
{% for vpcsecuritygroup in cluster.vpcsecuritygroups %}
          <member>
            <VpcSecurityGroupId>{{ vpcsecuritygroup.vpc_security_group_id }}</VpcSecurityGroupId>
            <Status>{{ vpcsecuritygroup.status }}</Status>
          </member>
{% endfor %}
        </VpcSecurityGroups>
        <HostedZoneId>{{ cluster.hosted_zone_id }}</HostedZoneId>
        <StorageEncrypted>{{ 'true' if cluster.storage_encrypted else 'false'}}</StorageEncrypted>
        <KmsKeyId>{{ cluster.kms_key_id }}</KmsKeyId>
        <DbClusterResourceId>{{ cluster.resource_id }}</DbClusterResourceId>
        <DBClusterArn>{{ cluster.db_cluster_arn }}</DBClusterArn>
        <AssociatedRoles>
{% for associatedrole in cluster.associatedroles %}
          <member>
            <RoleArn>{{ associatedrole.role_arn }}</RoleArn>
            <Status>{{ associatedrole.status }}</Status>
            <FeatureName>{{ associatedrole.feature_name }}</FeatureName>
          </member>
{% endfor %}
        </AssociatedRoles>
        <IAMDatabaseAuthenticationEnabled>false</IAMDatabaseAuthenticationEnabled>
        <CloneGroupId>{{ cluster.clone_group_id }}</CloneGroupId>
        <ClusterCreateTime>{{ cluster.cluster_create_time }}</ClusterCreateTime>
        <CopyTagsToSnapshot>false</CopyTagsToSnapshot>
        <EnabledCloudwatchLogsExports>
{% for enabledcloudwatchlogsexport in cluster.enabledcloudwatchlogsexports %}
          <member/>db_cluster_arn
{% endfor %}
        </EnabledCloudwatchLogsExports>
        <DeletionProtection>false</DeletionProtection>
        <CrossAccountClone>false</CrossAccountClone>
        {% if cluster.automatic_restart_time %}
        <AutomaticRestartTime>{{ cluster.automatic_restart_time }}</AutomaticRestartTime>
        {% endif %}
        {% if cluster.serverless_v2_scaling_configuration %}
        <ServerlessV2ScalingConfiguration>
            <MinCapacity>{{ cluster.serverless_v2_scaling_configuration["MinCapacity"] }}</MinCapacity>
            <MaxCapacity>{{ cluster.serverless_v2_scaling_configuration["MaxCapacity"] }}</MaxCapacity>
        </ServerlessV2ScalingConfiguration>
        {% endif %}
        </DBCluster>"""
        )
        return template.render(cluster=self)


class NeptuneBackend(BaseBackend):
    """Implementation of Neptune APIs."""

    def __init__(self, region_name: str, account_id: str):
        super().__init__(region_name, account_id)
        self.clusters: Dict[str, DBCluster] = dict()
        self.global_clusters: Dict[str, GlobalCluster] = dict()
        self._db_cluster_options: Optional[List[Dict[str, Any]]] = None

    @property
    def global_backend(self) -> "NeptuneBackend":
        return neptune_backends[self.account_id]["us-east-1"]

    @property
    def db_cluster_options(self) -> List[Dict[str, Any]]:  # type: ignore[misc]
        if self._db_cluster_options is None:
            from moto.rds.utils import decode_orderable_db_instance

            decoded_options: List[Dict[str, Any]] = load_resource(
                __name__, "../rds/resources/cluster_options/neptune.json"
            )
            self._db_cluster_options = [
                decode_orderable_db_instance(option) for option in decoded_options
            ]
        return self._db_cluster_options

    def create_db_cluster(self, **kwargs: Any) -> DBCluster:
        cluster = DBCluster(
            account_id=self.account_id,
            region_name=self.region_name,
            db_cluster_identifier=kwargs["db_cluster_identifier"],
            database_name=kwargs.get("database_name"),
            storage_encrypted=kwargs.get("storage_encrypted", True),
            parameter_group_name=kwargs.get("db_cluster_parameter_group_name") or "",
            tags=kwargs.get("tags", []),
            engine=kwargs.get("engine", "neptune"),
            engine_version=kwargs.get("engine_version") or "1.2.0.2",
            kms_key_id=kwargs.get("kms_key_id"),
            preferred_maintenance_window=kwargs.get("preferred_maintenance_window")
            or "none",
            preferred_backup_window=kwargs.get("preferred_backup_window"),
            backup_retention_period=kwargs.get("backup_retention_period") or 1,
            port=kwargs.get("port") or 8192,
            serverless_v2_scaling_configuration=kwargs.get(
                "serverless_v2_scaling_configuration"
            ),
        )
        self.clusters[cluster.db_cluster_identifier] = cluster
        return cluster

    def create_global_cluster(
        self,
        global_cluster_identifier: str,
        engine: Optional[str],
        engine_version: Optional[str],
        storage_encrypted: Optional[str],
        deletion_protection: Optional[str],
    ) -> GlobalCluster:
        cluster = GlobalCluster(
            account_id=self.account_id,
            region_name=self.region_name,
            global_cluster_identifier=global_cluster_identifier,
            engine=engine,
            engine_version=engine_version,
            storage_encrypted=storage_encrypted,
            deletion_protection=deletion_protection,
        )
        self.global_backend.global_clusters[global_cluster_identifier] = cluster
        return cluster

    def delete_global_cluster(self, global_cluster_identifier: str) -> GlobalCluster:
        return self.global_backend.global_clusters.pop(global_cluster_identifier)

    def describe_global_clusters(self) -> List[GlobalCluster]:
        return list(self.global_backend.global_clusters.values())

    def describe_db_clusters(self, db_cluster_identifier: str) -> List[DBCluster]:
        """
        Pagination and the Filters-argument is not yet implemented
        """
        if db_cluster_identifier:
            if db_cluster_identifier not in self.clusters:
                raise DBClusterNotFoundError(db_cluster_identifier)
            return [self.clusters[db_cluster_identifier]]
        return list(self.clusters.values())

    def delete_db_cluster(self, cluster_identifier: str) -> DBCluster:
        """
        The parameters SkipFinalSnapshot and FinalDBSnapshotIdentifier are not yet implemented.
        The DeletionProtection-attribute is not yet enforced
        """
        if cluster_identifier in self.clusters:
            return self.clusters.pop(cluster_identifier)
        raise DBClusterNotFoundError(cluster_identifier)

    def modify_db_cluster(self, kwargs: Any) -> DBCluster:
        cluster_id = kwargs["db_cluster_identifier"]

        cluster = self.clusters[cluster_id]
        del self.clusters[cluster_id]

        kwargs["db_cluster_identifier"] = kwargs.pop("new_db_cluster_identifier")
        for k, v in kwargs.items():
            if v is not None:
                setattr(cluster, k, v)

        cluster_id = kwargs.get("new_db_cluster_identifier", cluster_id)
        self.clusters[cluster_id] = cluster

        initial_state = copy.deepcopy(cluster)  # Return status=creating
        cluster.status = "available"  # Already set the final status in the background
        return initial_state

    def start_db_cluster(self, cluster_identifier: str) -> DBCluster:
        if cluster_identifier not in self.clusters:
            raise DBClusterNotFoundError(cluster_identifier)
        cluster = self.clusters[cluster_identifier]
        temp_state = copy.deepcopy(cluster)
        temp_state.status = "started"
        cluster.status = "available"  # This is the final status - already setting it in the background
        return temp_state

    def describe_orderable_db_instance_options(
        self, engine_version: Optional[str]
    ) -> List[Dict[str, Any]]:
        """
        Only the EngineVersion-parameter is currently implemented.
        """
        if engine_version:
            return [
                option
                for option in self.db_cluster_options
                if option["EngineVersion"] == engine_version
            ]
        return self.db_cluster_options


neptune_backends = BackendDict(NeptuneBackend, "neptune")
