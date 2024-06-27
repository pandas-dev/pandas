import copy
import os
import re
import string
from collections import OrderedDict, defaultdict
from re import compile as re_compile
from typing import TYPE_CHECKING, Any, Dict, Iterable, List, Optional, Tuple, Union

from jinja2 import Template

from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import BaseModel, CloudFormationModel
from moto.core.utils import iso_8601_datetime_with_milliseconds
from moto.ec2.models import ec2_backends
from moto.moto_api._internal import mock_random as random
from moto.neptune.models import NeptuneBackend, neptune_backends
from moto.utilities.utils import ARN_PARTITION_REGEX, get_partition, load_resource

from .exceptions import (
    DBClusterNotFoundError,
    DBClusterParameterGroupNotFoundError,
    DBClusterSnapshotAlreadyExistsError,
    DBClusterSnapshotNotFoundError,
    DBClusterToBeDeletedHasActiveMembers,
    DBInstanceNotFoundError,
    DBParameterGroupNotFoundError,
    DBProxyAlreadyExistsFault,
    DBProxyNotFoundFault,
    DBProxyQuotaExceededFault,
    DBSecurityGroupNotFoundError,
    DBSnapshotAlreadyExistsError,
    DBSnapshotNotFoundError,
    DBSubnetGroupNotFoundError,
    ExportTaskAlreadyExistsError,
    ExportTaskNotFoundError,
    InvalidDBClusterStateFault,
    InvalidDBClusterStateFaultError,
    InvalidDBInstanceEngine,
    InvalidDBInstanceIdentifier,
    InvalidDBInstanceStateError,
    InvalidExportSourceStateError,
    InvalidGlobalClusterStateFault,
    InvalidParameterCombination,
    InvalidParameterValue,
    InvalidSubnet,
    OptionGroupNotFoundFaultError,
    RDSClientError,
    SnapshotQuotaExceededError,
    SubscriptionAlreadyExistError,
    SubscriptionNotFoundError,
)
from .utils import (
    ClusterEngine,
    DbInstanceEngine,
    FilterDef,
    apply_filter,
    merge_filters,
    valid_preferred_maintenance_window,
    validate_filters,
)

if TYPE_CHECKING:
    from moto.ec2.models.subnets import Subnet


def find_cluster(cluster_arn: str) -> "Cluster":
    arn_parts = cluster_arn.split(":")
    region, account = arn_parts[3], arn_parts[4]
    return rds_backends[account][region].describe_db_clusters(cluster_arn)[0]


class GlobalCluster(BaseModel):
    def __init__(
        self,
        account_id: str,
        region_name: str,
        global_cluster_identifier: str,
        engine: str,
        engine_version: Optional[str],
        storage_encrypted: Optional[str],
        deletion_protection: Optional[str],
    ):
        self.global_cluster_identifier = global_cluster_identifier
        self.global_cluster_resource_id = "cluster-" + random.get_random_hex(8)
        self.global_cluster_arn = f"arn:{get_partition(region_name)}:rds::{account_id}:global-cluster:{global_cluster_identifier}"
        self.engine = engine
        self.engine_version = engine_version or "5.7.mysql_aurora.2.11.2"
        self.storage_encrypted = (
            storage_encrypted and storage_encrypted.lower() == "true"
        )
        self.deletion_protection = (
            deletion_protection and deletion_protection.lower() == "true"
        )
        self.members: List[Cluster] = []

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
          <DeletionProtection>{{ 'true' if cluster.deletion_protection else 'false' }}</DeletionProtection>
          <GlobalClusterMembers>
          {% for cluster_member in cluster.members %}
          <GlobalClusterMember>
            <DBClusterArn>{{ cluster_member.db_cluster_arn }}</DBClusterArn>
            <IsWriter>{{ 'true' if cluster_member.is_writer else 'false' }}</IsWriter>
            {% if not cluster_member.is_writer %}<GlobalWriteForwardingStatus>disabled</GlobalWriteForwardingStatus>{% endif %}
                <Readers>
                    {% if cluster_member.is_writer %}
                        {% for reader in cluster.members %}
                            {% if not reader.is_writer %}<Reader>{{ reader.db_cluster_arn }}</Reader>{% endif %}
                        {% endfor %}
                    {% endif %}
                </Readers>
          </GlobalClusterMember>
          {% endfor %}
          </GlobalClusterMembers>
          """
        )
        return template.render(cluster=self)


class Cluster:
    SUPPORTED_FILTERS = {
        "db-cluster-id": FilterDef(
            ["db_cluster_arn", "db_cluster_identifier"], "DB Cluster Identifiers"
        ),
        "engine": FilterDef(["engine"], "Engine Names"),
    }

    def __init__(self, **kwargs: Any):
        self.db_name = kwargs.get("db_name")
        self.db_cluster_identifier = kwargs.get("db_cluster_identifier")
        self.db_cluster_instance_class = kwargs.get("db_cluster_instance_class")
        self.deletion_protection = kwargs.get("deletion_protection")
        self.engine = kwargs.get("engine")
        if self.engine not in ClusterEngine.list_cluster_engines():
            raise InvalidParameterValue(
                (
                    "Engine '{engine}' is not supported "
                    "to satisfy constraint: Member must satisfy enum value set: "
                    "{valid_engines}"
                ).format(
                    engine=self.engine,
                    valid_engines=ClusterEngine.list_cluster_engines(),
                )
            )
        self.engine_version = kwargs.get(
            "engine_version"
        ) or Cluster.default_engine_version(self.engine)
        self.engine_mode = kwargs.get("engine_mode") or "provisioned"
        self.iops = kwargs.get("iops")
        self.kms_key_id = kwargs.get("kms_key_id")
        self.network_type = kwargs.get("network_type") or "IPV4"
        self.status = "active"
        self.account_id = kwargs.get("account_id")
        self.region_name = kwargs["region"]
        self.cluster_create_time = iso_8601_datetime_with_milliseconds()
        self.copy_tags_to_snapshot = kwargs.get("copy_tags_to_snapshot")
        if self.copy_tags_to_snapshot is None:
            self.copy_tags_to_snapshot = True
        self.storage_type = kwargs.get("storage_type")
        if self.storage_type is None:
            self.storage_type = Cluster.default_storage_type(iops=self.iops)
        self.allocated_storage = kwargs.get("allocated_storage")
        if self.allocated_storage is None:
            self.allocated_storage = Cluster.default_allocated_storage(
                engine=self.engine, storage_type=self.storage_type
            )
        self.master_username = kwargs.get("master_username")
        self.global_cluster_identifier = kwargs.get("global_cluster_identifier")
        if not self.master_username and self.global_cluster_identifier:
            pass
        elif not self.master_username:
            raise InvalidParameterValue(
                "The parameter MasterUsername must be provided and must not be blank."
            )
        else:
            self.master_user_password = kwargs.get("master_user_password")  # type: ignore

        self.availability_zones = kwargs.get("availability_zones")
        if not self.availability_zones:
            self.availability_zones = [
                f"{self.region_name}a",
                f"{self.region_name}b",
                f"{self.region_name}c",
            ]
        self.parameter_group = kwargs.get("parameter_group") or "default.aurora8.0"
        self.subnet_group = kwargs.get("db_subnet_group_name") or "default"
        self.status = "creating"
        self.url_identifier = "".join(
            random.choice(string.ascii_lowercase + string.digits) for _ in range(12)
        )
        self.endpoint = f"{self.db_cluster_identifier}.cluster-{self.url_identifier}.{self.region_name}.rds.amazonaws.com"
        self.reader_endpoint = f"{self.db_cluster_identifier}.cluster-ro-{self.url_identifier}.{self.region_name}.rds.amazonaws.com"
        self.port: int = kwargs.get("port")  # type: ignore
        if self.port is None:
            self.port = Cluster.default_port(self.engine)
        self.preferred_backup_window = "01:37-02:07"
        self.preferred_maintenance_window = "wed:02:40-wed:03:10"
        # This should default to the default security group
        self.vpc_security_group_ids: List[str] = kwargs["vpc_security_group_ids"]
        self.hosted_zone_id = "".join(
            random.choice(string.ascii_uppercase + string.digits) for _ in range(14)
        )
        self.resource_id = "cluster-" + "".join(
            random.choice(string.ascii_uppercase + string.digits) for _ in range(26)
        )
        self.tags = kwargs.get("tags", [])
        self.enabled_cloudwatch_logs_exports = (
            kwargs.get("enable_cloudwatch_logs_exports") or []
        )
        self.enable_http_endpoint = kwargs.get("enable_http_endpoint")  # type: ignore
        self.earliest_restorable_time = iso_8601_datetime_with_milliseconds()
        self.scaling_configuration = kwargs.get("scaling_configuration")
        if not self.scaling_configuration and self.engine_mode == "serverless":
            # In AWS, this default configuration only shows up when the Cluster is in a ready state, so a few minutes after creation
            self.scaling_configuration = {
                "min_capacity": 1,
                "max_capacity": 16,
                "auto_pause": True,
                "seconds_until_auto_pause": 300,
                "timeout_action": "RollbackCapacityChange",
                "seconds_before_timeout": 300,
            }
        self.serverless_v2_scaling_configuration = kwargs.get(
            "serverless_v2_scaling_configuration"
        )
        self.cluster_members: List[str] = list()
        self.replication_source_identifier = kwargs.get("replication_source_identifier")
        self.read_replica_identifiers: List[str] = list()
        self.is_writer: bool = False
        self.storage_encrypted = kwargs.get("storage_encrypted", False)
        if self.storage_encrypted:
            self.kms_key_id = kwargs.get("kms_key_id", "default_kms_key_id")
        else:
            self.kms_key_id = kwargs.get("kms_key_id")
        if self.engine == "aurora-mysql" or self.engine == "aurora-postgresql":
            self.global_write_forwarding_requested = kwargs.get(
                "enable_global_write_forwarding"
            )

    @property
    def is_multi_az(self) -> bool:
        return (
            len(self.read_replica_identifiers) > 0
            or self.replication_source_identifier is not None
        )

    @property
    def arn(self) -> str:
        return self.db_cluster_arn

    @property
    def db_cluster_arn(self) -> str:
        return f"arn:{get_partition(self.region_name)}:rds:{self.region_name}:{self.account_id}:cluster:{self.db_cluster_identifier}"

    @property
    def master_user_password(self) -> str:
        return self._master_user_password

    @master_user_password.setter
    def master_user_password(self, val: str) -> None:
        if not val:
            raise InvalidParameterValue(
                "The parameter MasterUserPassword must be provided and must not be blank."
            )
        if len(val) < 8:
            raise InvalidParameterValue(
                "The parameter MasterUserPassword is not a valid password because it is shorter than 8 characters."
            )
        self._master_user_password = val

    @property
    def enable_http_endpoint(self) -> bool:
        return self._enable_http_endpoint

    @enable_http_endpoint.setter
    def enable_http_endpoint(self, val: Optional[bool]) -> None:
        # instead of raising an error on aws rds create-db-cluster commands with
        # incompatible configurations with enable_http_endpoint
        # (e.g. engine_mode is not set to "serverless"), the API
        # automatically sets the enable_http_endpoint parameter to False
        self._enable_http_endpoint = False
        if val is not None:
            if self.engine_mode == "serverless":
                if self.engine == "aurora-mysql" and self.engine_version in [
                    "5.6.10a",
                    "5.6.1",
                    "2.07.1",
                    "5.7.2",
                    "5.7.mysql_aurora.2.07.1",
                    "5.7.mysql_aurora.2.07.2",
                    "5.7.mysql_aurora.2.08.3",
                ]:
                    self._enable_http_endpoint = val
                elif self.engine == "aurora-postgresql" and self.engine_version in [
                    "10.12",
                    "10.14",
                    "10.18",
                    "11.13",
                ]:
                    self._enable_http_endpoint = val
                elif self.engine == "aurora" and self.engine_version in [
                    "5.6.mysql_aurora.1.22.5"
                ]:
                    self._enable_http_endpoint = val

    def get_cfg(self) -> Dict[str, Any]:
        cfg = self.__dict__
        cfg["master_user_password"] = cfg.pop("_master_user_password")
        cfg["enable_http_endpoint"] = cfg.pop("_enable_http_endpoint")
        return cfg

    def to_xml(self) -> str:
        template = Template(
            """<DBCluster>
              <AllocatedStorage>{{ cluster.allocated_storage }}</AllocatedStorage>
              <AvailabilityZones>
              {% for zone in cluster.availability_zones %}
                  <AvailabilityZone>{{ zone }}</AvailabilityZone>
              {% endfor %}
              </AvailabilityZones>
              <BackupRetentionPeriod>1</BackupRetentionPeriod>
              <BacktrackWindow>0</BacktrackWindow>
              <DBInstanceStatus>{{ cluster.status }}</DBInstanceStatus>
              {% if cluster.db_name %}<DatabaseName>{{ cluster.db_name }}</DatabaseName>{% endif %}
              {% if cluster.kms_key_id %}<KmsKeyId>{{ cluster.kms_key_id }}</KmsKeyId>{% endif %}
              {% if cluster.network_type %}<NetworkType>{{ cluster.network_type }}</NetworkType>{% endif %}
              <DBClusterIdentifier>{{ cluster.db_cluster_identifier }}</DBClusterIdentifier>
              <DBClusterParameterGroup>{{ cluster.parameter_group }}</DBClusterParameterGroup>
              <DBSubnetGroup>{{ cluster.subnet_group }}</DBSubnetGroup>
              <ClusterCreateTime>{{ cluster.cluster_create_time }}</ClusterCreateTime>
              <EarliestRestorableTime>{{ cluster.earliest_restorable_time }}</EarliestRestorableTime>
              <Engine>{{ cluster.engine }}</Engine>
              <Status>{{ cluster.status }}</Status>
              <Endpoint>{{ cluster.endpoint }}</Endpoint>
              <ReaderEndpoint>{{ cluster.reader_endpoint }}</ReaderEndpoint>
              <MultiAZ>{{ 'true' if cluster.is_multi_az else 'false' }}</MultiAZ>
              <EngineVersion>{{ cluster.engine_version }}</EngineVersion>
              <Port>{{ cluster.port }}</Port>
              {% if cluster.iops %}
                <Iops>{{ cluster.iops }}</Iops>
                <StorageType>io1</StorageType>
              {% endif %}
              {% if cluster.db_cluster_instance_class %}<DBClusterInstanceClass>{{ cluster.db_cluster_instance_class }}</DBClusterInstanceClass>{% endif %}
              <MasterUsername>{{ cluster.master_username }}</MasterUsername>
              <PreferredBackupWindow>{{ cluster.preferred_backup_window }}</PreferredBackupWindow>
              <PreferredMaintenanceWindow>{{ cluster.preferred_maintenance_window }}</PreferredMaintenanceWindow>
              <ReadReplicaIdentifiers>
              {% for replica_id in cluster.read_replica_identifiers %}
                <ReadReplicaIdentifier>{{ replica_id }}</ReadReplicaIdentifier>
              {% endfor %}
              </ReadReplicaIdentifiers>
              <DBClusterMembers>
              {% for member in cluster.cluster_members %}
              <DBClusterMember>
                <DBInstanceIdentifier>{{ member }}</DBInstanceIdentifier>
                <IsClusterWriter>true</IsClusterWriter>
                <DBClusterParameterGroupStatus>in-sync</DBClusterParameterGroupStatus>
                <PromotionTier>1</PromotionTier>
              </DBClusterMember>
              {% endfor %}
              </DBClusterMembers>
              <VpcSecurityGroups>
              {% for id in cluster.vpc_security_group_ids %}
                  <VpcSecurityGroup>
                      <VpcSecurityGroupId>{{ id }}</VpcSecurityGroupId>
                      <Status>active</Status>
                  </VpcSecurityGroup>
              {% endfor %}
              </VpcSecurityGroups>
              <HostedZoneId>{{ cluster.hosted_zone_id }}</HostedZoneId>
              <StorageEncrypted>{{ 'true' if cluster.storage_encrypted else 'false' }}</StorageEncrypted>
              <GlobalWriteForwardingRequested>{{ cluster.global_write_forwarding_requested }}</GlobalWriteForwardingRequested>
              <DbClusterResourceId>{{ cluster.resource_id }}</DbClusterResourceId>
              <DBClusterArn>{{ cluster.db_cluster_arn }}</DBClusterArn>
              <AssociatedRoles></AssociatedRoles>
              <IAMDatabaseAuthenticationEnabled>false</IAMDatabaseAuthenticationEnabled>
              <EngineMode>{{ cluster.engine_mode }}</EngineMode>
              <DeletionProtection>{{ 'true' if cluster.deletion_protection else 'false' }}</DeletionProtection>
              <HttpEndpointEnabled>{{ 'true' if cluster.enable_http_endpoint else 'false' }}</HttpEndpointEnabled>
              <CopyTagsToSnapshot>{{ cluster.copy_tags_to_snapshot }}</CopyTagsToSnapshot>
              <CrossAccountClone>false</CrossAccountClone>
              <DomainMemberships></DomainMemberships>
              <EnabledCloudwatchLogsExports>
              {% for export in cluster.enabled_cloudwatch_logs_exports %}
              <member>{{ export }}</member>
              {% endfor %}
              </EnabledCloudwatchLogsExports>
              <TagList>
              {%- for tag in cluster.tags -%}
                <Tag>
                  <Key>{{ tag['Key'] }}</Key>
                  <Value>{{ tag['Value'] }}</Value>
                </Tag>
              {%- endfor -%}
              </TagList>
              {% if cluster.scaling_configuration %}
              <ScalingConfigurationInfo>
                {% if "min_capacity" in cluster.scaling_configuration %}<MinCapacity>{{ cluster.scaling_configuration["min_capacity"] }}</MinCapacity>{% endif %}
                {% if "max_capacity" in cluster.scaling_configuration %}<MaxCapacity>{{ cluster.scaling_configuration["max_capacity"] }}</MaxCapacity>{% endif %}
                {% if "auto_pause" in cluster.scaling_configuration %}<AutoPause>{{ cluster.scaling_configuration["auto_pause"] }}</AutoPause>{% endif %}
                {% if "seconds_until_auto_pause" in cluster.scaling_configuration %}<SecondsUntilAutoPause>{{ cluster.scaling_configuration["seconds_until_auto_pause"] }}</SecondsUntilAutoPause>{% endif %}
                {% if "timeout_action" in cluster.scaling_configuration %}<TimeoutAction>{{ cluster.scaling_configuration["timeout_action"] }}</TimeoutAction>{% endif %}
                {% if "seconds_before_timeout" in cluster.scaling_configuration %}<SecondsBeforeTimeout>{{ cluster.scaling_configuration["seconds_before_timeout"] }}</SecondsBeforeTimeout>{% endif %}
              </ScalingConfigurationInfo>
              {% endif %}
              {% if cluster.serverless_v2_scaling_configuration %}
              <ServerlessV2ScalingConfiguration>
                {% if "MinCapacity" in cluster.serverless_v2_scaling_configuration %}<MinCapacity>{{ cluster.serverless_v2_scaling_configuration["MinCapacity"] }}</MinCapacity>{% endif %}
                {% if "MaxCapacity" in cluster.serverless_v2_scaling_configuration %}<MaxCapacity>{{ cluster.serverless_v2_scaling_configuration["MaxCapacity"] }}</MaxCapacity>{% endif %}
              </ServerlessV2ScalingConfiguration>
              {% endif %}
              {%- if cluster.global_cluster_identifier -%}
              <GlobalClusterIdentifier>{{ cluster.global_cluster_identifier }}</GlobalClusterIdentifier>
              {%- endif -%}
              {%- if cluster.replication_source_identifier -%}<ReplicationSourceIdentifier>{{ cluster.replication_source_identifier }}</ReplicationSourceIdentifier>{%- endif -%}
            </DBCluster>"""
        )
        return template.render(cluster=self)

    @staticmethod
    def default_engine_version(engine: str) -> str:
        return {
            "aurora": "5.6.mysql_aurora.1.22.5",
            "aurora-mysql": "5.7.mysql_aurora.2.07.2",
            "aurora-postgresql": "12.7",
            "mysql": "8.0.23",
            "postgres": "13.4",
        }[engine]

    @staticmethod
    def default_port(engine: str) -> int:
        return {
            "aurora": 3306,
            "aurora-mysql": 3306,
            "aurora-postgresql": 5432,
            "mysql": 3306,
            "postgres": 5432,
        }[engine]

    @staticmethod
    def default_storage_type(iops: Any) -> str:  # type: ignore[misc]
        if iops is None:
            return "gp2"
        else:
            return "io1"

    @staticmethod
    def default_allocated_storage(engine: str, storage_type: str) -> int:
        return {
            "aurora": {"gp2": 0, "io1": 0, "standard": 0},
            "aurora-mysql": {"gp2": 20, "io1": 100, "standard": 10},
            "aurora-postgresql": {"gp2": 20, "io1": 100, "standard": 10},
            "mysql": {"gp2": 20, "io1": 100, "standard": 5},
            "postgres": {"gp2": 20, "io1": 100, "standard": 5},
        }[engine][storage_type]

    def get_tags(self) -> List[Dict[str, str]]:
        return self.tags

    def add_tags(self, tags: List[Dict[str, str]]) -> List[Dict[str, str]]:
        new_keys = [tag_set["Key"] for tag_set in tags]
        self.tags = [tag_set for tag_set in self.tags if tag_set["Key"] not in new_keys]
        self.tags.extend(tags)
        return self.tags

    def remove_tags(self, tag_keys: List[str]) -> None:
        self.tags = [tag_set for tag_set in self.tags if tag_set["Key"] not in tag_keys]


class ClusterSnapshot(BaseModel):
    SUPPORTED_FILTERS = {
        "db-cluster-id": FilterDef(
            ["cluster.db_cluster_arn", "cluster.db_cluster_identifier"],
            "DB Cluster Identifiers",
        ),
        "db-cluster-snapshot-id": FilterDef(
            ["snapshot_id"], "DB Cluster Snapshot Identifiers"
        ),
        "snapshot-type": FilterDef(["snapshot_type"], "Snapshot Types"),
        "engine": FilterDef(["cluster.engine"], "Engine Names"),
    }

    def __init__(
        self,
        cluster: Cluster,
        snapshot_id: str,
        snapshot_type: str,
        tags: List[Dict[str, str]],
    ):
        self.cluster = cluster
        self.snapshot_id = snapshot_id
        self.snapshot_type = snapshot_type
        self.tags = tags
        self.status = "available"
        self.created_at = iso_8601_datetime_with_milliseconds()
        self.attributes: List[Dict[str, Any]] = []

    @property
    def arn(self) -> str:
        return self.snapshot_arn

    @property
    def snapshot_arn(self) -> str:
        return f"arn:{get_partition(self.cluster.region_name)}:rds:{self.cluster.region_name}:{self.cluster.account_id}:cluster-snapshot:{self.snapshot_id}"

    def to_xml(self) -> str:
        template = Template(
            """
            <DBClusterSnapshot>
                <DBClusterSnapshotIdentifier>{{ snapshot.snapshot_id }}</DBClusterSnapshotIdentifier>
                <SnapshotCreateTime>{{ snapshot.created_at }}</SnapshotCreateTime>
                <DBClusterIdentifier>{{ cluster.db_cluster_identifier }}</DBClusterIdentifier>
                <ClusterCreateTime>{{ snapshot.created_at }}</ClusterCreateTime>
                <PercentProgress>{{ 100 }}</PercentProgress>
                <AllocatedStorage>{{ cluster.allocated_storage }}</AllocatedStorage>
                <MasterUsername>{{ cluster.master_username }}</MasterUsername>
                <Port>{{ cluster.port }}</Port>
                <Engine>{{ cluster.engine }}</Engine>
                <Status>{{ snapshot.status }}</Status>
                <SnapshotType>{{ snapshot.snapshot_type }}</SnapshotType>
                <DBClusterSnapshotArn>{{ snapshot.snapshot_arn }}</DBClusterSnapshotArn>
                <SourceRegion>{{ cluster.region }}</SourceRegion>
                {% if cluster.iops %}
                <Iops>{{ cluster.iops }}</Iops>
                <StorageType>io1</StorageType>
                {% else %}
                <StorageType>{{ cluster.storage_type }}</StorageType>
                {% endif %}
                <TagList>
                {%- for tag in snapshot.tags -%}
                    <Tag><Key>{{ tag['Key'] }}</Key><Value>{{ tag['Value'] }}</Value></Tag>
                {%- endfor -%}
                </TagList>
                <Timezone></Timezone>
                <LicenseModel>{{ cluster.license_model }}</LicenseModel>
            </DBClusterSnapshot>
            """
        )
        return template.render(snapshot=self, cluster=self.cluster)

    def get_tags(self) -> List[Dict[str, str]]:
        return self.tags

    def add_tags(self, tags: List[Dict[str, str]]) -> List[Dict[str, str]]:
        new_keys = [tag_set["Key"] for tag_set in tags]
        self.tags = [tag_set for tag_set in self.tags if tag_set["Key"] not in new_keys]
        self.tags.extend(tags)
        return self.tags

    def remove_tags(self, tag_keys: List[str]) -> None:
        self.tags = [tag_set for tag_set in self.tags if tag_set["Key"] not in tag_keys]


class Database(CloudFormationModel):
    SUPPORTED_FILTERS = {
        "db-cluster-id": FilterDef(["db_cluster_identifier"], "DB Cluster Identifiers"),
        "db-instance-id": FilterDef(
            ["db_instance_arn", "db_instance_identifier"], "DB Instance Identifiers"
        ),
        "dbi-resource-id": FilterDef(["dbi_resource_id"], "Dbi Resource Ids"),
        "domain": FilterDef(None, ""),
        "engine": FilterDef(["engine"], "Engine Names"),
    }

    default_engine_versions = {
        "MySQL": "5.6.21",
        "mysql": "5.6.21",
        "oracle-se1": "11.2.0.4.v3",
        "oracle-se": "11.2.0.4.v3",
        "oracle-ee": "11.2.0.4.v3",
        "sqlserver-ee": "11.00.2100.60.v1",
        "sqlserver-se": "11.00.2100.60.v1",
        "sqlserver-ex": "11.00.2100.60.v1",
        "sqlserver-web": "11.00.2100.60.v1",
        "postgres": "9.3.3",
    }

    def __init__(self, **kwargs: Any):
        self.status = "available"
        self.is_replica = False
        self.replicas: List[str] = []
        self.account_id: str = kwargs["account_id"]
        self.region_name: str = kwargs["region"]
        self.engine = kwargs.get("engine")
        if self.engine not in DbInstanceEngine.valid_db_instance_engine():
            raise InvalidParameterValue(
                f"Value {self.engine} for parameter Engine is invalid. Reason: engine {self.engine} not supported"
            )
        self.engine_version = kwargs.get("engine_version", None)
        if not self.engine_version and self.engine in self.default_engine_versions:
            self.engine_version = self.default_engine_versions[self.engine]
        self.iops = kwargs.get("iops")
        self.storage_encrypted = kwargs.get("storage_encrypted", False)
        if self.storage_encrypted:
            self.kms_key_id = kwargs.get("kms_key_id", "default_kms_key_id")
        else:
            self.kms_key_id = kwargs.get("kms_key_id")
        self.storage_type = kwargs.get("storage_type")
        if self.storage_type is None:
            self.storage_type = Database.default_storage_type(iops=self.iops)
        self.master_username = kwargs.get("master_username")
        self.master_user_password = kwargs.get("master_user_password")
        self.auto_minor_version_upgrade = kwargs.get("auto_minor_version_upgrade")
        if self.auto_minor_version_upgrade is None:
            self.auto_minor_version_upgrade = True
        self.allocated_storage = kwargs.get("allocated_storage")
        if self.allocated_storage is None:
            self.allocated_storage = Database.default_allocated_storage(
                engine=self.engine, storage_type=self.storage_type
            )
        self.db_cluster_identifier: Optional[str] = kwargs.get("db_cluster_identifier")
        self.db_instance_identifier = kwargs.get("db_instance_identifier")
        self.source_db_identifier: Optional[str] = kwargs.get("source_db_identifier")
        self.db_instance_class = kwargs.get("db_instance_class")
        self.port = kwargs.get("port")
        if self.port is None:
            self.port = Database.default_port(self.engine)
        self.db_instance_identifier = kwargs.get("db_instance_identifier")
        self.db_name = kwargs.get("db_name")
        self.instance_create_time = iso_8601_datetime_with_milliseconds()
        self.publicly_accessible = kwargs.get("publicly_accessible")
        if self.publicly_accessible is None:
            self.publicly_accessible = True
        self.copy_tags_to_snapshot = kwargs.get("copy_tags_to_snapshot")
        if self.copy_tags_to_snapshot is None:
            self.copy_tags_to_snapshot = False
        self.backup_retention_period = kwargs.get("backup_retention_period")
        if self.backup_retention_period is None:
            self.backup_retention_period = 1
        self.availability_zone = kwargs.get("availability_zone")
        if not self.availability_zone:
            self.availability_zone = f"{self.region_name}a"
        self.multi_az = kwargs.get("multi_az")
        self.db_subnet_group_name = kwargs.get("db_subnet_group_name")
        self.db_subnet_group = None
        if self.db_subnet_group_name:
            self.db_subnet_group = rds_backends[self.account_id][
                self.region_name
            ].describe_db_subnet_groups(self.db_subnet_group_name)[0]
        self.security_groups = kwargs.get("security_groups", [])
        self.vpc_security_group_ids = kwargs.get("vpc_security_group_ids", [])
        self.preferred_maintenance_window = kwargs.get("preferred_maintenance_window")
        self.preferred_backup_window = kwargs.get("preferred_backup_window")
        msg = valid_preferred_maintenance_window(
            self.preferred_maintenance_window,
            self.preferred_backup_window,
        )
        if msg:
            raise RDSClientError("InvalidParameterValue", msg)

        self.db_parameter_group_name = kwargs.get("db_parameter_group_name")
        if (
            self.db_parameter_group_name
            and not self.is_default_parameter_group(self.db_parameter_group_name)
            and self.db_parameter_group_name
            not in rds_backends[self.account_id][self.region_name].db_parameter_groups
        ):
            raise DBParameterGroupNotFoundError(self.db_parameter_group_name)

        self.license_model = kwargs.get("license_model", "general-public-license")
        self.option_group_name = kwargs.get("option_group_name", None)
        self.option_group_supplied = self.option_group_name is not None
        if (
            self.option_group_name
            and self.option_group_name
            not in rds_backends[self.account_id][self.region_name].option_groups
        ):
            raise OptionGroupNotFoundFaultError(self.option_group_name)
        self.default_option_groups = {
            "MySQL": "default.mysql5.6",
            "mysql": "default.mysql5.6",
            "postgres": "default.postgres9.3",
        }
        if not self.option_group_name and self.engine in self.default_option_groups:
            self.option_group_name = self.default_option_groups[self.engine]
        self.character_set_name = kwargs.get("character_set_name", None)
        self.enable_iam_database_authentication = kwargs.get(
            "enable_iam_database_authentication", False
        )
        self.dbi_resource_id = "db-M5ENSHXFPU6XHZ4G4ZEI5QIO2U"
        self.tags = kwargs.get("tags", [])
        self.deletion_protection = kwargs.get("deletion_protection", False)
        self.enabled_cloudwatch_logs_exports = (
            kwargs.get("enable_cloudwatch_logs_exports") or []
        )

    @property
    def arn(self) -> str:
        return self.db_instance_arn

    @property
    def db_instance_arn(self) -> str:
        return f"arn:{get_partition(self.region_name)}:rds:{self.region_name}:{self.account_id}:db:{self.db_instance_identifier}"

    @property
    def physical_resource_id(self) -> Optional[str]:
        return self.db_instance_identifier

    def db_parameter_groups(self) -> List["DBParameterGroup"]:
        if not self.db_parameter_group_name or self.is_default_parameter_group(
            self.db_parameter_group_name
        ):
            (
                db_family,
                db_parameter_group_name,
            ) = self.default_db_parameter_group_details()
            description = f"Default parameter group for {db_family}"
            return [
                DBParameterGroup(
                    account_id=self.account_id,
                    name=db_parameter_group_name,
                    family=db_family,
                    description=description,
                    tags=[],
                    region=self.region_name,
                )
            ]
        else:
            backend = rds_backends[self.account_id][self.region_name]
            if self.db_parameter_group_name not in backend.db_parameter_groups:
                raise DBParameterGroupNotFoundError(self.db_parameter_group_name)

            return [backend.db_parameter_groups[self.db_parameter_group_name]]

    def is_default_parameter_group(self, param_group_name: str) -> bool:
        return param_group_name.startswith(f"default.{self.engine.lower()}")  # type: ignore

    def default_db_parameter_group_details(self) -> Tuple[Optional[str], Optional[str]]:
        if not self.engine_version:
            return (None, None)

        minor_engine_version = ".".join(str(self.engine_version).rsplit(".")[:-1])
        db_family = f"{self.engine.lower()}{minor_engine_version}"  # type: ignore

        return db_family, f"default.{db_family}"

    def to_xml(self) -> str:
        template = Template(
            """<DBInstance>
              <AvailabilityZone>{{ database.availability_zone }}</AvailabilityZone>
              <BackupRetentionPeriod>{{ database.backup_retention_period }}</BackupRetentionPeriod>
              <DBInstanceStatus>{{ database.status }}</DBInstanceStatus>
              {% if database.db_name %}<DBName>{{ database.db_name }}</DBName>{% endif %}
              <MultiAZ>{{ 'true' if database.multi_az else 'false' }}</MultiAZ>
              <VpcSecurityGroups>
                {% for vpc_security_group_id in database.vpc_security_group_ids %}
                <VpcSecurityGroupMembership>
                  <Status>active</Status>
                  <VpcSecurityGroupId>{{ vpc_security_group_id }}</VpcSecurityGroupId>
                </VpcSecurityGroupMembership>
                {% endfor %}
              </VpcSecurityGroups>
              {% if database.db_cluster_identifier %}<DBClusterIdentifier>{{ database.db_cluster_identifier }}</DBClusterIdentifier>{% endif %}
              <DBInstanceIdentifier>{{ database.db_instance_identifier }}</DBInstanceIdentifier>
              <DbiResourceId>{{ database.dbi_resource_id }}</DbiResourceId>
              <InstanceCreateTime>{{ database.instance_create_time }}</InstanceCreateTime>
              <PreferredBackupWindow>{{ database.preferred_backup_window }}</PreferredBackupWindow>
              <PreferredMaintenanceWindow>{{ database.preferred_maintenance_window }}</PreferredMaintenanceWindow>
              <ReadReplicaDBInstanceIdentifiers>
                {% for replica_id in database.replicas %}
                    <ReadReplicaDBInstanceIdentifier>{{ replica_id }}</ReadReplicaDBInstanceIdentifier>
                {% endfor %}
              </ReadReplicaDBInstanceIdentifiers>
              <StatusInfos>
                {% if database.is_replica %}
                <DBInstanceStatusInfo>
                    <StatusType>read replication</StatusType>
                    <Status>replicating</Status>
                    <Normal>true</Normal>
                    <Message></Message>
                </DBInstanceStatusInfo>
                {% endif %}
              </StatusInfos>
              <EnabledCloudwatchLogsExports>
              {% for export in database.enabled_cloudwatch_logs_exports %}
              <member>{{ export }}</member>
              {% endfor %}
              </EnabledCloudwatchLogsExports>
              {% if database.is_replica %}
              <ReadReplicaSourceDBInstanceIdentifier>{{ database.source_db_identifier }}</ReadReplicaSourceDBInstanceIdentifier>
              {% endif %}
              <Engine>{{ database.engine }}</Engine>
              <IAMDatabaseAuthenticationEnabled>{{'true' if database.enable_iam_database_authentication else 'false' }}</IAMDatabaseAuthenticationEnabled>
              <LicenseModel>{{ database.license_model }}</LicenseModel>
              <EngineVersion>{{ database.engine_version }}</EngineVersion>
              <OptionGroupMemberships>
                <OptionGroupMembership>
                  <OptionGroupName>{{ database.option_group_name }}</OptionGroupName>
                  <Status>in-sync</Status>
                </OptionGroupMembership>
              </OptionGroupMemberships>
              <DBParameterGroups>
                {% for db_parameter_group in database.db_parameter_groups() %}
                <DBParameterGroup>
                  <ParameterApplyStatus>in-sync</ParameterApplyStatus>
                  <DBParameterGroupName>{{ db_parameter_group.name }}</DBParameterGroupName>
                </DBParameterGroup>
                {% endfor %}
              </DBParameterGroups>
              <DBSecurityGroups>
                {% for security_group in database.security_groups %}
                <DBSecurityGroup>
                  <Status>active</Status>
                  <DBSecurityGroupName>{{ security_group }}</DBSecurityGroupName>
                </DBSecurityGroup>
                {% endfor %}
              </DBSecurityGroups>
              {% if database.db_subnet_group %}
              <DBSubnetGroup>
                <DBSubnetGroupName>{{ database.db_subnet_group.subnet_name }}</DBSubnetGroupName>
                <DBSubnetGroupDescription>{{ database.db_subnet_group.description }}</DBSubnetGroupDescription>
                <SubnetGroupStatus>{{ database.db_subnet_group.status }}</SubnetGroupStatus>
                <Subnets>
                    {% for subnet in database.db_subnet_group.subnets %}
                    <Subnet>
                      <SubnetStatus>Active</SubnetStatus>
                      <SubnetIdentifier>{{ subnet.id }}</SubnetIdentifier>
                      <SubnetAvailabilityZone>
                        <Name>{{ subnet.availability_zone }}</Name>
                        <ProvisionedIopsCapable>false</ProvisionedIopsCapable>
                      </SubnetAvailabilityZone>
                    </Subnet>
                    {% endfor %}
                </Subnets>
                <VpcId>{{ database.db_subnet_group.vpc_id }}</VpcId>
              </DBSubnetGroup>
              {% endif %}
              <PubliclyAccessible>{{ database.publicly_accessible }}</PubliclyAccessible>
              <CopyTagsToSnapshot>{{ database.copy_tags_to_snapshot }}</CopyTagsToSnapshot>
              <AutoMinorVersionUpgrade>{{ database.auto_minor_version_upgrade }}</AutoMinorVersionUpgrade>
              <AllocatedStorage>{{ database.allocated_storage }}</AllocatedStorage>
              <StorageEncrypted>{{ database.storage_encrypted }}</StorageEncrypted>
              {% if database.kms_key_id %}
              <KmsKeyId>{{ database.kms_key_id }}</KmsKeyId>
              {% endif %}
              {% if database.iops %}
              <Iops>{{ database.iops }}</Iops>
              <StorageType>io1</StorageType>
              {% else %}
              <StorageType>{{ database.storage_type }}</StorageType>
              {% endif %}
              <DBInstanceClass>{{ database.db_instance_class }}</DBInstanceClass>
              <MasterUsername>{{ database.master_username }}</MasterUsername>
              <Endpoint>
                <Address>{{ database.address }}</Address>
                <Port>{{ database.port }}</Port>
              </Endpoint>
              <DbInstancePort>{{ database.port }}</DbInstancePort>
              <DBInstanceArn>{{ database.db_instance_arn }}</DBInstanceArn>
              <TagList>
              {%- for tag in database.tags -%}
                <Tag>
                  <Key>{{ tag['Key'] }}</Key>
                  <Value>{{ tag['Value'] }}</Value>
                </Tag>
              {%- endfor -%}
              </TagList>
              <DeletionProtection>{{ 'true' if database.deletion_protection else 'false' }}</DeletionProtection>
            </DBInstance>"""
        )
        return template.render(database=self)

    @property
    def address(self) -> str:
        return f"{self.db_instance_identifier}.aaaaaaaaaa.{self.region_name}.rds.amazonaws.com"

    def add_replica(self, replica: "Database") -> None:
        if self.region_name != replica.region_name:
            # Cross Region replica
            self.replicas.append(replica.db_instance_arn)
        else:
            self.replicas.append(replica.db_instance_identifier)  # type: ignore

    def remove_replica(self, replica: "Database") -> None:
        self.replicas.remove(replica.db_instance_identifier)  # type: ignore

    def set_as_replica(self) -> None:
        self.is_replica = True
        self.replicas = []

    def update(self, db_kwargs: Dict[str, Any]) -> None:
        for key, value in db_kwargs.items():
            if value is not None:
                setattr(self, key, value)

    @classmethod
    def has_cfn_attr(cls, attr: str) -> bool:
        return attr in ["Endpoint.Address", "Endpoint.Port"]

    def get_cfn_attribute(self, attribute_name: str) -> Any:
        # Local import to avoid circular dependency with cloudformation.parsing
        from moto.cloudformation.exceptions import UnformattedGetAttTemplateException

        if attribute_name == "Endpoint.Address":
            return self.address
        elif attribute_name == "Endpoint.Port":
            return self.port
        raise UnformattedGetAttTemplateException()

    @staticmethod
    def default_port(engine: str) -> int:
        return {
            "aurora": 3306,
            "aurora-mysql": 3306,
            "aurora-postgresql": 5432,
            "mysql": 3306,
            "mariadb": 3306,
            "postgres": 5432,
            "oracle-ee": 1521,
            "oracle-se2": 1521,
            "oracle-se1": 1521,
            "oracle-se": 1521,
            "sqlserver-ee": 1433,
            "sqlserver-ex": 1433,
            "sqlserver-se": 1433,
            "sqlserver-web": 1433,
        }[engine]

    @staticmethod
    def default_storage_type(iops: Any) -> str:  # type: ignore[misc]
        if iops is None:
            return "gp2"
        else:
            return "io1"

    @staticmethod
    def default_allocated_storage(engine: str, storage_type: str) -> int:
        return {
            "aurora": {"gp2": 0, "io1": 0, "standard": 0},
            "aurora-mysql": {"gp2": 20, "io1": 100, "standard": 10},
            "aurora-postgresql": {"gp2": 20, "io1": 100, "standard": 10},
            "mysql": {"gp2": 20, "io1": 100, "standard": 5},
            "mariadb": {"gp2": 20, "io1": 100, "standard": 5},
            "postgres": {"gp2": 20, "io1": 100, "standard": 5},
            "oracle-ee": {"gp2": 20, "io1": 100, "standard": 10},
            "oracle-se2": {"gp2": 20, "io1": 100, "standard": 10},
            "oracle-se1": {"gp2": 20, "io1": 100, "standard": 10},
            "oracle-se": {"gp2": 20, "io1": 100, "standard": 10},
            "sqlserver-ee": {"gp2": 200, "io1": 200, "standard": 200},
            "sqlserver-ex": {"gp2": 20, "io1": 100, "standard": 20},
            "sqlserver-se": {"gp2": 200, "io1": 200, "standard": 200},
            "sqlserver-web": {"gp2": 20, "io1": 100, "standard": 20},
        }[engine][storage_type]

    @staticmethod
    def cloudformation_name_type() -> str:
        return "DBInstanceIdentifier"

    @staticmethod
    def cloudformation_type() -> str:
        # https://docs.aws.amazon.com/AWSCloudFormation/latest/UserGuide/aws-resource-rds-dbinstance.html
        return "AWS::RDS::DBInstance"

    @classmethod
    def create_from_cloudformation_json(  # type: ignore[misc]
        cls,
        resource_name: str,
        cloudformation_json: Any,
        account_id: str,
        region_name: str,
        **kwargs: Any,
    ) -> "Database":
        properties = cloudformation_json["Properties"]

        db_security_groups = properties.get("DBSecurityGroups")
        if not db_security_groups:
            db_security_groups = []
        security_groups = [group.group_name for group in db_security_groups]
        db_subnet_group = properties.get("DBSubnetGroupName")
        db_subnet_group_name = db_subnet_group.subnet_name if db_subnet_group else None
        db_kwargs = {
            "auto_minor_version_upgrade": properties.get("AutoMinorVersionUpgrade"),
            "allocated_storage": properties.get("AllocatedStorage"),
            "availability_zone": properties.get("AvailabilityZone"),
            "backup_retention_period": properties.get("BackupRetentionPeriod"),
            "db_instance_class": properties.get("DBInstanceClass"),
            "db_instance_identifier": resource_name.replace("_", "-"),
            "db_name": properties.get("DBName"),
            "preferred_backup_window": properties.get(
                "PreferredBackupWindow", "13:14-13:44"
            ),
            "preferred_maintenance_window": properties.get(
                "PreferredMaintenanceWindow", "wed:06:38-wed:07:08"
            ).lower(),
            "db_subnet_group_name": db_subnet_group_name,
            "engine": properties.get("Engine"),
            "engine_version": properties.get("EngineVersion"),
            "iops": properties.get("Iops"),
            "kms_key_id": properties.get("KmsKeyId"),
            "master_user_password": properties.get("MasterUserPassword"),
            "master_username": properties.get("MasterUsername"),
            "multi_az": properties.get("MultiAZ"),
            "db_parameter_group_name": properties.get("DBParameterGroupName"),
            "port": properties.get("Port", 3306),
            "publicly_accessible": properties.get("PubliclyAccessible"),
            "copy_tags_to_snapshot": properties.get("CopyTagsToSnapshot"),
            "account_id": account_id,
            "region": region_name,
            "security_groups": security_groups,
            "storage_encrypted": properties.get("StorageEncrypted"),
            "storage_type": properties.get("StorageType"),
            "tags": properties.get("Tags"),
            "vpc_security_group_ids": properties.get("VpcSecurityGroupIds", []),
        }

        rds_backend = rds_backends[account_id][region_name]
        source_db_identifier = properties.get("SourceDBInstanceIdentifier")
        if source_db_identifier:
            # Replica
            db_kwargs["source_db_identifier"] = source_db_identifier
            database = rds_backend.create_db_instance_read_replica(db_kwargs)
        else:
            database = rds_backend.create_db_instance(db_kwargs)
        return database

    def to_json(self) -> str:
        template = Template(
            """{
        "AllocatedStorage": 10,
        "AutoMinorVersionUpgrade": "{{ database.auto_minor_version_upgrade }}",
        "AvailabilityZone": "{{ database.availability_zone }}",
        "BackupRetentionPeriod": "{{ database.backup_retention_period }}",
        "CharacterSetName": {%- if database.character_set_name -%}{{ database.character_set_name }}{%- else %} null{%- endif -%},
        "DBInstanceClass": "{{ database.db_instance_class }}",
        {%- if database.db_cluster_identifier -%}"DBClusterIdentifier": "{{ database.db_cluster_identifier }}",{%- endif -%}
        "DBInstanceIdentifier": "{{ database.db_instance_identifier }}",
        "DBInstanceStatus": "{{ database.status }}",
        "DBName": {%- if database.db_name -%}"{{ database.db_name }}"{%- else %} null{%- endif -%},
        {% if database.db_parameter_group_name -%}"DBParameterGroups": {
            "DBParameterGroup": {
            "ParameterApplyStatus": "in-sync",
            "DBParameterGroupName": "{{ database.db_parameter_group_name }}"
          }
        },{%- endif %}
        "DBSecurityGroups": [
          {% for security_group in database.security_groups -%}{%- if loop.index != 1 -%},{%- endif -%}
          {"DBSecurityGroup": {
            "Status": "active",
            "DBSecurityGroupName": "{{ security_group }}"
          }}{% endfor %}
        ],
        {%- if database.db_subnet_group -%}{{ database.db_subnet_group.to_json() }},{%- endif %}
        "Engine": "{{ database.engine }}",
        "EngineVersion": "{{ database.engine_version }}",
        "LatestRestorableTime": null,
        "LicenseModel": "{{ database.license_model }}",
        "MasterUsername": "{{ database.master_username }}",
        "MultiAZ": "{{ database.multi_az }}",{% if database.option_group_name %}
        "OptionGroupMemberships": [{
          "OptionGroupMembership": {
            "OptionGroupName": "{{ database.option_group_name }}",
            "Status": "in-sync"
          }
        }],{%- endif %}
        "PendingModifiedValues": { "MasterUserPassword": "****" },
        "PreferredBackupWindow": "{{ database.preferred_backup_window }}",
        "PreferredMaintenanceWindow": "{{ database.preferred_maintenance_window }}",
        "PubliclyAccessible": "{{ database.publicly_accessible }}",
        "CopyTagsToSnapshot": "{{ database.copy_tags_to_snapshot }}",
        "AllocatedStorage": "{{ database.allocated_storage }}",
        "Endpoint": {
            "Address": "{{ database.address }}",
            "Port": "{{ database.port }}"
        },
        "InstanceCreateTime": "{{ database.instance_create_time }}",
        "Iops": null,
        "ReadReplicaDBInstanceIdentifiers": [{%- for replica in database.replicas -%}
            {%- if not loop.first -%},{%- endif -%}
            "{{ replica }}"
        {%- endfor -%}
        ],
        {%- if database.source_db_identifier -%}
        "ReadReplicaSourceDBInstanceIdentifier": "{{ database.source_db_identifier }}",
        {%- else -%}
        "ReadReplicaSourceDBInstanceIdentifier": null,
        {%- endif -%}
        "SecondaryAvailabilityZone": null,
        "StatusInfos": null,
        "VpcSecurityGroups": [
            {% for vpc_security_group_id in database.vpc_security_group_ids %}
            {
                "Status": "active",
                "VpcSecurityGroupId": "{{ vpc_security_group_id }}"
            }
            {% endfor %}
        ],
        "DBInstanceArn": "{{ database.db_instance_arn }}"
      }"""
        )
        return template.render(database=self)

    def get_tags(self) -> List[Dict[str, str]]:
        return self.tags

    def add_tags(self, tags: List[Dict[str, str]]) -> List[Dict[str, str]]:
        new_keys = [tag_set["Key"] for tag_set in tags]
        self.tags = [tag_set for tag_set in self.tags if tag_set["Key"] not in new_keys]
        self.tags.extend(tags)
        return self.tags

    def remove_tags(self, tag_keys: List[str]) -> None:
        self.tags = [tag_set for tag_set in self.tags if tag_set["Key"] not in tag_keys]

    def delete(self, account_id: str, region_name: str) -> None:
        backend = rds_backends[account_id][region_name]
        backend.delete_db_instance(self.db_instance_identifier)  # type: ignore[arg-type]


class DatabaseSnapshot(BaseModel):
    SUPPORTED_FILTERS = {
        "db-instance-id": FilterDef(
            ["database.db_instance_arn", "database.db_instance_identifier"],
            "DB Instance Identifiers",
        ),
        "db-snapshot-id": FilterDef(["snapshot_id"], "DB Snapshot Identifiers"),
        "dbi-resource-id": FilterDef(["database.dbi_resource_id"], "Dbi Resource Ids"),
        "snapshot-type": FilterDef(["snapshot_type"], "Snapshot Types"),
        "engine": FilterDef(["database.engine"], "Engine Names"),
    }

    def __init__(
        self,
        database: Database,
        snapshot_id: str,
        snapshot_type: str,
        tags: List[Dict[str, str]],
    ):
        self.database = database
        self.snapshot_id = snapshot_id
        self.snapshot_type = snapshot_type
        self.tags = tags
        self.status = "available"
        self.created_at = iso_8601_datetime_with_milliseconds()
        self.attributes: List[Dict[str, Any]] = []

    @property
    def arn(self) -> str:
        return self.snapshot_arn

    @property
    def snapshot_arn(self) -> str:
        return f"arn:{get_partition(self.database.region_name)}:rds:{self.database.region_name}:{self.database.account_id}:snapshot:{self.snapshot_id}"

    def to_xml(self) -> str:
        template = Template(
            """<DBSnapshot>
              <DBSnapshotIdentifier>{{ snapshot.snapshot_id }}</DBSnapshotIdentifier>
              <DBInstanceIdentifier>{{ database.db_instance_identifier }}</DBInstanceIdentifier>
              <DbiResourceId>{{ database.dbi_resource_id }}</DbiResourceId>
              <SnapshotCreateTime>{{ snapshot.created_at }}</SnapshotCreateTime>
              <Engine>{{ database.engine }}</Engine>
              <AllocatedStorage>{{ database.allocated_storage }}</AllocatedStorage>
              <Status>{{ snapshot.status }}</Status>
              <Port>{{ database.port }}</Port>
              <AvailabilityZone>{{ database.availability_zone }}</AvailabilityZone>
              <VpcId>{{ database.db_subnet_group.vpc_id }}</VpcId>
              <InstanceCreateTime>{{ snapshot.created_at }}</InstanceCreateTime>
              {% if database.master_username %}
              <MasterUsername>{{ database.master_username }}</MasterUsername>
              {% endif %}
              <EngineVersion>{{ database.engine_version }}</EngineVersion>
              {% if database.license_model %}
              <LicenseModel>{{ database.license_model }}</LicenseModel>
              {% endif %}
              <SnapshotType>{{ snapshot.snapshot_type }}</SnapshotType>
              {% if database.iops %}
              <Iops>{{ database.iops }}</Iops>
              <StorageType>io1</StorageType>
              {% else %}
              <StorageType>{{ database.storage_type }}</StorageType>
              {% endif %}
              <OptionGroupName>{{ database.option_group_name }}</OptionGroupName>
              <PercentProgress>{{ 100 }}</PercentProgress>
              <SourceRegion>{{ database.region }}</SourceRegion>
              <SourceDBSnapshotIdentifier></SourceDBSnapshotIdentifier>
              <TagList>
              {%- for tag in snapshot.tags -%}
                <Tag><Key>{{ tag['Key'] }}</Key><Value>{{ tag['Value'] }}</Value></Tag>
              {%- endfor -%}
              </TagList>
              <TdeCredentialArn></TdeCredentialArn>
              <Encrypted>{{ database.storage_encrypted }}</Encrypted>
              {% if database.kms_key_id %}
              <KmsKeyId>{{ database.kms_key_id }}</KmsKeyId>
              {% endif %}
              <DBSnapshotArn>{{ snapshot.snapshot_arn }}</DBSnapshotArn>
              <Timezone></Timezone>
              {% if database.enable_iam_database_authentication %}
              <IAMDatabaseAuthenticationEnabled>{{ database.enable_iam_database_authentication|lower }}</IAMDatabaseAuthenticationEnabled>
              {% endif %}
            </DBSnapshot>"""
        )
        return template.render(snapshot=self, database=self.database)

    def get_tags(self) -> List[Dict[str, str]]:
        return self.tags

    def add_tags(self, tags: List[Dict[str, str]]) -> List[Dict[str, str]]:
        new_keys = [tag_set["Key"] for tag_set in tags]
        self.tags = [tag_set for tag_set in self.tags if tag_set["Key"] not in new_keys]
        self.tags.extend(tags)
        return self.tags

    def remove_tags(self, tag_keys: List[str]) -> None:
        self.tags = [tag_set for tag_set in self.tags if tag_set["Key"] not in tag_keys]


class ExportTask(BaseModel):
    def __init__(
        self, snapshot: Union[DatabaseSnapshot, ClusterSnapshot], kwargs: Dict[str, Any]
    ):
        self.snapshot = snapshot

        self.export_task_identifier = kwargs.get("export_task_identifier")
        self.kms_key_id = kwargs.get("kms_key_id", "default_kms_key_id")
        self.source_arn = kwargs.get("source_arn")
        self.iam_role_arn = kwargs.get("iam_role_arn")
        self.s3_bucket_name = kwargs.get("s3_bucket_name")
        self.s3_prefix = kwargs.get("s3_prefix", "")
        self.export_only = kwargs.get("export_only", [])

        self.status = "complete"
        self.created_at = iso_8601_datetime_with_milliseconds()
        self.source_type = (
            "SNAPSHOT" if type(snapshot) is DatabaseSnapshot else "CLUSTER"
        )

    def to_xml(self) -> str:
        template = Template(
            """
            <ExportTaskIdentifier>{{ task.export_task_identifier }}</ExportTaskIdentifier>
            <SourceArn>{{ snapshot.snapshot_arn }}</SourceArn>
            <TaskStartTime>{{ task.created_at }}</TaskStartTime>
            <TaskEndTime>{{ task.created_at }}</TaskEndTime>
            <SnapshotTime>{{ snapshot.created_at }}</SnapshotTime>
            <S3Bucket>{{ task.s3_bucket_name }}</S3Bucket>
            <S3Prefix>{{ task.s3_prefix }}</S3Prefix>
            <IamRoleArn>{{ task.iam_role_arn }}</IamRoleArn>
            <KmsKeyId>{{ task.kms_key_id }}</KmsKeyId>
            {%- if task.export_only -%}
            <ExportOnly>
                {%- for table in task.export_only -%}
                    <member>{{ table }}</member>
                {%- endfor -%}
            </ExportOnly>
            {%- endif -%}
            <Status>{{ task.status }}</Status>
            <PercentProgress>{{ 100 }}</PercentProgress>
            <TotalExtractedDataInGB>{{ 1 }}</TotalExtractedDataInGB>
            <FailureCause></FailureCause>
            <WarningMessage></WarningMessage>
            <SourceType>{{ task.source_type }}</SourceType>
            """
        )
        return template.render(task=self, snapshot=self.snapshot)


class EventSubscription(BaseModel):
    def __init__(self, kwargs: Dict[str, Any]):
        self.subscription_name = kwargs.get("subscription_name")
        self.sns_topic_arn = kwargs.get("sns_topic_arn")
        self.source_type = kwargs.get("source_type")
        self.event_categories = kwargs.get("event_categories", [])
        self.source_ids = kwargs.get("source_ids", [])
        self.enabled = kwargs.get("enabled", True)
        self.tags = kwargs.get("tags", True)

        self.region_name = ""
        self.customer_aws_id = kwargs["account_id"]
        self.status = "active"
        self.created_at = iso_8601_datetime_with_milliseconds()

    @property
    def es_arn(self) -> str:
        return f"arn:{get_partition(self.region_name)}:rds:{self.region_name}:{self.customer_aws_id}:es:{self.subscription_name}"

    def to_xml(self) -> str:
        template = Template(
            """
            <EventSubscription>
              <CustomerAwsId>{{ subscription.customer_aws_id }}</CustomerAwsId>
              <CustSubscriptionId>{{ subscription.subscription_name }}</CustSubscriptionId>
              <SnsTopicArn>{{ subscription.sns_topic_arn }}</SnsTopicArn>
              <SubscriptionCreationTime>{{ subscription.created_at }}</SubscriptionCreationTime>
              <SourceType>{{ subscription.source_type }}</SourceType>
              <SourceIdsList>
                {%- for source_id in subscription.source_ids -%}
                  <SourceId>{{ source_id }}</SourceId>
                {%- endfor -%}
              </SourceIdsList>
              <EventCategoriesList>
                {%- for category in subscription.event_categories -%}
                  <EventCategory>{{ category }}</EventCategory>
                {%- endfor -%}
              </EventCategoriesList>
              <Status>{{ subscription.status }}</Status>
              <Enabled>{{ subscription.enabled }}</Enabled>
              <EventSubscriptionArn>{{ subscription.es_arn }}</EventSubscriptionArn>
              <TagList>
              {%- for tag in subscription.tags -%}
                <Tag><Key>{{ tag['Key'] }}</Key><Value>{{ tag['Value'] }}</Value></Tag>
              {%- endfor -%}
              </TagList>
            </EventSubscription>
            """
        )
        return template.render(subscription=self)

    def get_tags(self) -> List[Dict[str, str]]:
        return self.tags

    def add_tags(self, tags: List[Dict[str, str]]) -> List[Dict[str, str]]:
        new_keys = [tag_set["Key"] for tag_set in tags]
        self.tags = [tag_set for tag_set in self.tags if tag_set["Key"] not in new_keys]
        self.tags.extend(tags)
        return self.tags

    def remove_tags(self, tag_keys: List[str]) -> None:
        self.tags = [tag_set for tag_set in self.tags if tag_set["Key"] not in tag_keys]


class SecurityGroup(CloudFormationModel):
    def __init__(
        self,
        account_id: str,
        group_name: str,
        description: str,
        tags: List[Dict[str, str]],
    ):
        self.group_name = group_name
        self.description = description
        self.status = "authorized"
        self.ip_ranges: List[Any] = []
        self.ec2_security_groups: List[Any] = []
        self.tags = tags
        self.owner_id = account_id
        self.vpc_id = None

    def to_xml(self) -> str:
        template = Template(
            """<DBSecurityGroup>
            <EC2SecurityGroups>
            {% for security_group in security_group.ec2_security_groups %}
                <EC2SecurityGroup>
                    <EC2SecurityGroupId>{{ security_group.id }}</EC2SecurityGroupId>
                    <EC2SecurityGroupName>{{ security_group.name }}</EC2SecurityGroupName>
                    <EC2SecurityGroupOwnerId>{{ security_group.owner_id }}</EC2SecurityGroupOwnerId>
                    <Status>authorized</Status>
                </EC2SecurityGroup>
            {% endfor %}
            </EC2SecurityGroups>

            <DBSecurityGroupDescription>{{ security_group.description }}</DBSecurityGroupDescription>
            <IPRanges>
            {% for ip_range in security_group.ip_ranges %}
                <IPRange>
                    <CIDRIP>{{ ip_range }}</CIDRIP>
                    <Status>authorized</Status>
                </IPRange>
            {% endfor %}
            </IPRanges>
            <OwnerId>{{ security_group.ownder_id }}</OwnerId>
            <DBSecurityGroupName>{{ security_group.group_name }}</DBSecurityGroupName>
        </DBSecurityGroup>"""
        )
        return template.render(security_group=self)

    def to_json(self) -> str:
        template = Template(
            """{
            "DBSecurityGroupDescription": "{{ security_group.description }}",
            "DBSecurityGroupName": "{{ security_group.group_name }}",
            "EC2SecurityGroups": {{ security_group.ec2_security_groups }},
            "IPRanges": [{%- for ip in security_group.ip_ranges -%}
                         {%- if loop.index != 1 -%},{%- endif -%}
                         "{{ ip }}"
                         {%- endfor -%}
                        ],
            "OwnerId": "{{ security_group.owner_id }}",
            "VpcId": "{{ security_group.vpc_id }}"
        }"""
        )
        return template.render(security_group=self)

    def authorize_cidr(self, cidr_ip: str) -> None:
        self.ip_ranges.append(cidr_ip)

    def authorize_security_group(self, security_group: str) -> None:
        self.ec2_security_groups.append(security_group)

    @staticmethod
    def cloudformation_name_type() -> str:
        return ""

    @staticmethod
    def cloudformation_type() -> str:
        # https://docs.aws.amazon.com/AWSCloudFormation/latest/UserGuide/aws-resource-rds-dbsecuritygroup.html
        return "AWS::RDS::DBSecurityGroup"

    @classmethod
    def create_from_cloudformation_json(  # type: ignore[misc]
        cls,
        resource_name: str,
        cloudformation_json: Any,
        account_id: str,
        region_name: str,
        **kwargs: Any,
    ) -> "SecurityGroup":
        properties = cloudformation_json["Properties"]
        group_name = resource_name.lower()
        description = properties["GroupDescription"]
        security_group_ingress_rules = properties.get("DBSecurityGroupIngress", [])
        tags = properties.get("Tags")

        ec2_backend = ec2_backends[account_id][region_name]
        rds_backend = rds_backends[account_id][region_name]
        security_group = rds_backend.create_db_security_group(
            group_name, description, tags
        )
        for security_group_ingress in security_group_ingress_rules:
            for ingress_type, ingress_value in security_group_ingress.items():
                if ingress_type == "CIDRIP":
                    security_group.authorize_cidr(ingress_value)
                elif ingress_type == "EC2SecurityGroupName":
                    subnet = ec2_backend.get_security_group_from_name(ingress_value)
                    security_group.authorize_security_group(subnet)  # type: ignore[arg-type]
                elif ingress_type == "EC2SecurityGroupId":
                    subnet = ec2_backend.get_security_group_from_id(ingress_value)
                    security_group.authorize_security_group(subnet)  # type: ignore[arg-type]
        return security_group

    def get_tags(self) -> List[Dict[str, str]]:
        return self.tags

    def add_tags(self, tags: List[Dict[str, str]]) -> List[Dict[str, str]]:
        new_keys = [tag_set["Key"] for tag_set in tags]
        self.tags = [tag_set for tag_set in self.tags if tag_set["Key"] not in new_keys]
        self.tags.extend(tags)
        return self.tags

    def remove_tags(self, tag_keys: List[str]) -> None:
        self.tags = [tag_set for tag_set in self.tags if tag_set["Key"] not in tag_keys]

    def delete(self, account_id: str, region_name: str) -> None:
        backend = rds_backends[account_id][region_name]
        backend.delete_security_group(self.group_name)


class SubnetGroup(CloudFormationModel):
    def __init__(
        self,
        subnet_name: str,
        description: str,
        subnets: "List[Subnet]",
        tags: List[Dict[str, str]],
        region: str,
        account_id: str,
    ):
        self.subnet_name = subnet_name
        self.description = description
        self.subnets = subnets
        self.status = "Complete"
        self.tags = tags
        self.vpc_id = self.subnets[0].vpc_id
        self.region = region
        self.account_id = account_id

    @property
    def sg_arn(self) -> str:
        return f"arn:{get_partition(self.region)}:rds:{self.region}:{self.account_id}:subgrp:{self.subnet_name}"

    def to_xml(self) -> str:
        template = Template(
            """<DBSubnetGroup>
              <VpcId>{{ subnet_group.vpc_id }}</VpcId>
              <SubnetGroupStatus>{{ subnet_group.status }}</SubnetGroupStatus>
              <DBSubnetGroupDescription>{{ subnet_group.description }}</DBSubnetGroupDescription>
              <DBSubnetGroupName>{{ subnet_group.subnet_name }}</DBSubnetGroupName>
              <DBSubnetGroupArn>{{ subnet_group.sg_arn }}</DBSubnetGroupArn>
              <Subnets>
                {% for subnet in subnet_group.subnets %}
                <Subnet>
                  <SubnetStatus>Active</SubnetStatus>
                  <SubnetIdentifier>{{ subnet.id }}</SubnetIdentifier>
                  <SubnetAvailabilityZone>
                    <Name>{{ subnet.availability_zone }}</Name>
                    <ProvisionedIopsCapable>false</ProvisionedIopsCapable>
                  </SubnetAvailabilityZone>
                </Subnet>
                {% endfor %}
              </Subnets>
            </DBSubnetGroup>"""
        )
        return template.render(subnet_group=self)

    def to_json(self) -> str:
        template = Template(
            """"DBSubnetGroup": {
                "VpcId": "{{ subnet_group.vpc_id }}",
                "SubnetGroupStatus": "{{ subnet_group.status }}",
                "DBSubnetGroupDescription": "{{ subnet_group.description }}",
                "DBSubnetGroupName": "{{ subnet_group.subnet_name }}",
                "Subnets": {
                  "Subnet": [
                    {% for subnet in subnet_group.subnets %}{
                      "SubnetStatus": "Active",
                      "SubnetIdentifier": "{{ subnet.id }}",
                      "SubnetAvailabilityZone": {
                        "Name": "{{ subnet.availability_zone }}",
                        "ProvisionedIopsCapable": "false"
                      }
                    }{%- if not loop.last -%},{%- endif -%}{% endfor %}
                  ]
                }
            }"""
        )
        return template.render(subnet_group=self)

    @staticmethod
    def cloudformation_name_type() -> str:
        return "DBSubnetGroupName"

    @staticmethod
    def cloudformation_type() -> str:
        # https://docs.aws.amazon.com/AWSCloudFormation/latest/UserGuide/aws-resource-rds-dbsubnetgroup.html
        return "AWS::RDS::DBSubnetGroup"

    @classmethod
    def create_from_cloudformation_json(  # type: ignore[misc]
        cls,
        resource_name: str,
        cloudformation_json: Any,
        account_id: str,
        region_name: str,
        **kwargs: Any,
    ) -> "SubnetGroup":
        properties = cloudformation_json["Properties"]

        description = properties["DBSubnetGroupDescription"]
        subnet_ids = properties["SubnetIds"]
        tags = properties.get("Tags")

        ec2_backend = ec2_backends[account_id][region_name]
        subnets = [ec2_backend.get_subnet(subnet_id) for subnet_id in subnet_ids]
        rds_backend = rds_backends[account_id][region_name]
        subnet_group = rds_backend.create_subnet_group(
            resource_name,
            description,
            subnets,
            tags,
        )
        return subnet_group

    def get_tags(self) -> List[Dict[str, str]]:
        return self.tags

    def add_tags(self, tags: List[Dict[str, str]]) -> List[Dict[str, str]]:
        new_keys = [tag_set["Key"] for tag_set in tags]
        self.tags = [tag_set for tag_set in self.tags if tag_set["Key"] not in new_keys]
        self.tags.extend(tags)
        return self.tags

    def remove_tags(self, tag_keys: List[str]) -> None:
        self.tags = [tag_set for tag_set in self.tags if tag_set["Key"] not in tag_keys]

    def delete(self, account_id: str, region_name: str) -> None:
        backend = rds_backends[account_id][region_name]
        backend.delete_subnet_group(self.subnet_name)


class DBProxy(BaseModel):
    def __init__(
        self,
        db_proxy_name: str,
        engine_family: str,
        auth: List[Dict[str, str]],
        role_arn: str,
        vpc_subnet_ids: List[str],
        region_name: str,
        account_id: str,
        vpc_security_group_ids: Optional[List[str]],
        require_tls: Optional[bool] = False,
        idle_client_timeout: Optional[int] = 1800,
        debug_logging: Optional[bool] = False,
        tags: Optional[List[Dict[str, str]]] = None,
    ):
        self.db_proxy_name = db_proxy_name
        self.engine_family = engine_family
        if self.engine_family not in ["MYSQL", "POSTGRESQ", "SQLSERVER"]:
            raise InvalidParameterValue("Provided EngineFamily is not valid.")
        self.auth = auth
        self.role_arn = role_arn
        self.vpc_subnet_ids = vpc_subnet_ids
        self.vpc_security_group_ids = vpc_security_group_ids
        self.require_tls = require_tls
        if idle_client_timeout is None:
            self.idle_client_timeout = 1800
        else:
            if int(idle_client_timeout) < 1:
                self.idle_client_timeout = 1
            elif int(idle_client_timeout) > 28800:
                self.idle_client_timeout = 28800
            else:
                self.idle_client_timeout = idle_client_timeout
        self.debug_logging = debug_logging
        self.created_date = iso_8601_datetime_with_milliseconds()
        self.updated_date = iso_8601_datetime_with_milliseconds()
        if tags is None:
            self.tags = []
        else:
            self.tags = tags
        self.region_name = region_name
        self.account_id = account_id
        self.db_proxy_arn = f"arn:{get_partition(self.region_name)}:rds:{self.region_name}:{self.account_id}:db-proxy:{self.db_proxy_name}"
        self.arn = self.db_proxy_arn
        ec2_backend = ec2_backends[self.account_id][self.region_name]
        subnets = ec2_backend.describe_subnets(subnet_ids=self.vpc_subnet_ids)
        vpcs = []
        for subnet in subnets:
            vpcs.append(subnet.vpc_id)
            if subnet.vpc_id != vpcs[0]:
                raise InvalidSubnet(subnet_identifier=subnet.id)

        self.vpc_id = ec2_backend.describe_subnets(subnet_ids=[self.vpc_subnet_ids[0]])[
            0
        ].vpc_id
        self.status = "availible"
        self.url_identifier = "".join(
            random.choice(string.ascii_lowercase + string.digits) for _ in range(12)
        )
        self.endpoint = f"{self.db_proxy_name}.db-proxy-{self.url_identifier}.{self.region_name}.rds.amazonaws.com"

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
            """
                <RequireTLS>{{ dbproxy.require_tls }}</RequireTLS>
                <VpcSecurityGroupIds>
                {% if dbproxy.VpcSecurityGroupIds %}
                  {% for vpcsecuritygroupid in dbproxy.VpcSecurityGroupIds %}
                    <member>{{ vpcsecuritygroupid }}</member>
                  {% endfor %}
                {% endif %}
                </VpcSecurityGroupIds>
                <Auth>
                  {% for auth in dbproxy.auth %}
                    <member>
                        <UserName>{{ auth["UserName"] }}</UserName>
                        <AuthScheme>{{ auth["AuthScheme"] }}</AuthScheme>
                        <SecretArn>{{ auth["SecretArn"] }}</SecretArn>
                        <IAMAuth>{{ auth["IAMAuth"] }}</IAMAuth>
                        <ClientPasswordAuthType>{{ auth["ClientPasswordAuthType"] }}</ClientPasswordAuthType>
                    </member>
                  {% endfor %}
                </Auth>
                <EngineFamily>{{ dbproxy.engine_family }}</EngineFamily>
                <UpdatedDate>{{ dbproxy.updated_date }}</UpdatedDate>
                <DBProxyName>{{ dbproxy.db_proxy_name }}</DBProxyName>
                <IdleClientTimeout>{{ dbproxy.idle_client_timeout }}</IdleClientTimeout>
                <Endpoint>{{ dbproxy.endpoint }}</Endpoint>
                <CreatedDate>{{ dbproxy.created_date }}</CreatedDate>
                <RoleArn>{{ dbproxy.role_arn }}</RoleArn>
                <DebugLogging>{{ dbproxy.debug_logging }}</DebugLogging>
                <VpcId>{{ dbproxy.vpc_id }}</VpcId>
                <DBProxyArn>{{ dbproxy.db_proxy_arn }}</DBProxyArn>
                <VpcSubnetIds>
                  {% for vpcsubnetid in dbproxy.vpc_subnet_ids %}
                    <member>{{ vpcsubnetid }}</member>
                  {% endfor %}
                </VpcSubnetIds>
                <Status>{{ dbproxy.status }}</Status>
        """
        )
        return template.render(dbproxy=self)


class RDSBackend(BaseBackend):
    def __init__(self, region_name: str, account_id: str):
        super().__init__(region_name, account_id)
        self.arn_regex = re_compile(
            ARN_PARTITION_REGEX
            + r":rds:.*:[0-9]*:(db|cluster|es|og|pg|ri|secgrp|snapshot|cluster-snapshot|subgrp|db-proxy):.*$"
        )
        self.clusters: Dict[str, Cluster] = OrderedDict()
        self.global_clusters: Dict[str, GlobalCluster] = OrderedDict()
        self.databases: Dict[str, Database] = OrderedDict()
        self.database_snapshots: Dict[str, DatabaseSnapshot] = OrderedDict()
        self.cluster_snapshots: Dict[str, ClusterSnapshot] = OrderedDict()
        self.export_tasks: Dict[str, ExportTask] = OrderedDict()
        self.event_subscriptions: Dict[str, EventSubscription] = OrderedDict()
        self.db_parameter_groups: Dict[str, DBParameterGroup] = {}
        self.db_cluster_parameter_groups: Dict[str, DBClusterParameterGroup] = {}
        self.option_groups: Dict[str, OptionGroup] = {}
        self.security_groups: Dict[str, SecurityGroup] = {}
        self.subnet_groups: Dict[str, SubnetGroup] = {}
        self._db_cluster_options: Optional[List[Dict[str, Any]]] = None
        self.db_proxies: Dict[str, DBProxy] = OrderedDict()

    def reset(self) -> None:
        self.neptune.reset()
        super().reset()

    @property
    def neptune(self) -> NeptuneBackend:
        return neptune_backends[self.account_id][self.region_name]

    @property
    def db_cluster_options(self) -> List[Dict[str, Any]]:  # type: ignore
        if self._db_cluster_options is None:
            from moto.rds.utils import decode_orderable_db_instance

            decoded_options = load_resource(
                __name__, "resources/cluster_options/aurora-postgresql.json"
            )
            self._db_cluster_options = [
                decode_orderable_db_instance(option) for option in decoded_options
            ]
        return self._db_cluster_options

    def create_db_instance(self, db_kwargs: Dict[str, Any]) -> Database:
        database_id = db_kwargs["db_instance_identifier"]
        self._validate_db_identifier(database_id)
        database = Database(**db_kwargs)

        cluster_id = database.db_cluster_identifier
        if cluster_id is not None:
            cluster = self.clusters.get(cluster_id)
            if cluster is not None:
                if (
                    cluster.engine in ClusterEngine.serverless_engines()
                    and cluster.engine_mode == "serverless"
                ):
                    raise InvalidParameterValue(
                        "Instances cannot be added to Aurora Serverless clusters."
                    )
                if database.engine != cluster.engine:
                    raise InvalidDBInstanceEngine(
                        str(database.engine), str(cluster.engine)
                    )
                cluster.cluster_members.append(database_id)
        self.databases[database_id] = database
        return database

    def create_auto_snapshot(
        self,
        db_instance_identifier: str,
        db_snapshot_identifier: str,
    ) -> DatabaseSnapshot:
        return self.create_db_snapshot(
            db_instance_identifier, db_snapshot_identifier, snapshot_type="automated"
        )

    def create_db_snapshot(
        self,
        db_instance: Union[str, Database],
        db_snapshot_identifier: str,
        snapshot_type: str = "manual",
        tags: Optional[List[Dict[str, str]]] = None,
    ) -> DatabaseSnapshot:
        if isinstance(db_instance, str):
            database = self.databases.get(db_instance)
            if not database:
                raise DBInstanceNotFoundError(db_instance)
        else:
            database = db_instance

        if db_snapshot_identifier in self.database_snapshots:
            raise DBSnapshotAlreadyExistsError(db_snapshot_identifier)
        if len(self.database_snapshots) >= int(
            os.environ.get("MOTO_RDS_SNAPSHOT_LIMIT", "100")
        ):
            raise SnapshotQuotaExceededError()
        if tags is None:
            tags = list()
        if database.copy_tags_to_snapshot and not tags:
            tags = database.get_tags()
        snapshot = DatabaseSnapshot(
            database, db_snapshot_identifier, snapshot_type, tags
        )
        self.database_snapshots[db_snapshot_identifier] = snapshot
        return snapshot

    def copy_db_snapshot(
        self,
        source_snapshot_identifier: str,
        target_snapshot_identifier: str,
        tags: Optional[List[Dict[str, str]]] = None,
    ) -> DatabaseSnapshot:
        if source_snapshot_identifier not in self.database_snapshots:
            raise DBSnapshotNotFoundError(source_snapshot_identifier)

        source_snapshot = self.database_snapshots[source_snapshot_identifier]
        if tags is None:
            tags = source_snapshot.tags
        else:
            tags = self._merge_tags(source_snapshot.tags, tags)
        return self.create_db_snapshot(
            db_instance=source_snapshot.database,
            db_snapshot_identifier=target_snapshot_identifier,
            tags=tags,
        )

    def delete_db_snapshot(self, db_snapshot_identifier: str) -> DatabaseSnapshot:
        if db_snapshot_identifier not in self.database_snapshots:
            raise DBSnapshotNotFoundError(db_snapshot_identifier)

        return self.database_snapshots.pop(db_snapshot_identifier)

    def promote_read_replica(self, db_kwargs: Dict[str, Any]) -> Database:
        database_id = db_kwargs["db_instance_identifier"]
        database = self.databases[database_id]
        if database.is_replica:
            database.is_replica = False
            database.update(db_kwargs)

        return database

    def create_db_instance_read_replica(self, db_kwargs: Dict[str, Any]) -> Database:
        database_id = db_kwargs["db_instance_identifier"]
        source_database_id = db_kwargs["source_db_identifier"]
        primary = self.find_db_from_id(source_database_id)
        if self.arn_regex.match(source_database_id):
            db_kwargs["region"] = self.region_name

        # Shouldn't really copy here as the instance is duplicated. RDS replicas have different instances.
        replica = copy.copy(primary)
        replica.update(db_kwargs)
        replica.region_name = self.region_name
        replica.set_as_replica()
        self.databases[database_id] = replica
        primary.add_replica(replica)
        return replica

    def describe_db_instances(
        self, db_instance_identifier: Optional[str] = None, filters: Any = None
    ) -> List[Database]:
        databases = self.databases
        if db_instance_identifier:
            filters = merge_filters(
                filters, {"db-instance-id": [db_instance_identifier]}
            )
        if filters:
            databases = self._filter_resources(databases, filters, Database)
        if db_instance_identifier and not databases:
            raise DBInstanceNotFoundError(db_instance_identifier)
        return list(databases.values())

    def describe_db_snapshots(
        self,
        db_instance_identifier: Optional[str],
        db_snapshot_identifier: str,
        filters: Optional[Dict[str, Any]] = None,
    ) -> List[DatabaseSnapshot]:
        snapshots = self.database_snapshots
        if db_instance_identifier:
            filters = merge_filters(
                filters, {"db-instance-id": [db_instance_identifier]}
            )
        if db_snapshot_identifier:
            filters = merge_filters(
                filters, {"db-snapshot-id": [db_snapshot_identifier]}
            )
        if filters:
            snapshots = self._filter_resources(snapshots, filters, DatabaseSnapshot)
        if db_snapshot_identifier and not snapshots and not db_instance_identifier:
            raise DBSnapshotNotFoundError(db_snapshot_identifier)
        return list(snapshots.values())

    def modify_db_instance(
        self, db_instance_identifier: str, db_kwargs: Dict[str, Any]
    ) -> Database:
        database = self.describe_db_instances(db_instance_identifier)[0]
        if "new_db_instance_identifier" in db_kwargs:
            del self.databases[db_instance_identifier]
            db_instance_identifier = db_kwargs["db_instance_identifier"] = (
                db_kwargs.pop("new_db_instance_identifier")
            )
            self.databases[db_instance_identifier] = database
        preferred_backup_window = db_kwargs.get("preferred_backup_window")
        preferred_maintenance_window = db_kwargs.get("preferred_maintenance_window")
        msg = valid_preferred_maintenance_window(
            preferred_maintenance_window, preferred_backup_window
        )
        if msg:
            raise RDSClientError("InvalidParameterValue", msg)
        database.update(db_kwargs)
        return database

    def reboot_db_instance(self, db_instance_identifier: str) -> Database:
        return self.describe_db_instances(db_instance_identifier)[0]

    def restore_db_instance_from_db_snapshot(
        self, from_snapshot_id: str, overrides: Dict[str, Any]
    ) -> Database:
        snapshot = self.describe_db_snapshots(
            db_instance_identifier=None, db_snapshot_identifier=from_snapshot_id
        )[0]
        original_database = snapshot.database
        new_instance_props = copy.deepcopy(original_database.__dict__)
        if not original_database.option_group_supplied:
            # If the option group is not supplied originally, the 'option_group_name' will receive a default value
            # Force this reconstruction, and prevent any validation on the default value
            del new_instance_props["option_group_name"]

        for key, value in overrides.items():
            if value:
                new_instance_props[key] = value

        return self.create_db_instance(new_instance_props)

    def restore_db_instance_to_point_in_time(
        self,
        source_db_identifier: str,
        target_db_identifier: str,
        overrides: Dict[str, Any],
    ) -> Database:
        db_instance = self.describe_db_instances(
            db_instance_identifier=source_db_identifier
        )[0]

        # remove the db subnet group as it cannot be copied
        # and is not used in the restored instance
        source_dict = db_instance.__dict__
        del source_dict["db_subnet_group"]

        new_instance_props = copy.deepcopy(source_dict)
        if not db_instance.option_group_supplied:
            # If the option group is not supplied originally, the 'option_group_name' will receive a default value
            # Force this reconstruction, and prevent any validation on the default value
            del new_instance_props["option_group_name"]

        for key, value in overrides.items():
            if value:
                new_instance_props[key] = value

        # set the new db instance identifier
        new_instance_props["db_instance_identifier"] = target_db_identifier

        return self.create_db_instance(new_instance_props)

    def stop_db_instance(
        self, db_instance_identifier: str, db_snapshot_identifier: Optional[str] = None
    ) -> Database:
        self._validate_db_identifier(db_instance_identifier)
        database = self.describe_db_instances(db_instance_identifier)[0]
        # todo: certain rds types not allowed to be stopped at this time.
        # https://docs.aws.amazon.com/AmazonRDS/latest/UserGuide/USER_StopInstance.html#USER_StopInstance.Limitations
        if database.is_replica or (
            database.multi_az and database.engine.lower().startswith("sqlserver")  # type: ignore
        ):
            # todo: more db types not supported by stop/start instance api
            raise InvalidDBClusterStateFaultError(db_instance_identifier)
        if database.status != "available":
            raise InvalidDBInstanceStateError(db_instance_identifier, "stop")
        if db_snapshot_identifier:
            self.create_auto_snapshot(db_instance_identifier, db_snapshot_identifier)
        database.status = "stopped"
        return database

    def start_db_instance(self, db_instance_identifier: str) -> Database:
        self._validate_db_identifier(db_instance_identifier)
        database = self.describe_db_instances(db_instance_identifier)[0]
        # todo: bunch of different error messages to be generated from this api call
        if database.status != "stopped":
            raise InvalidDBInstanceStateError(db_instance_identifier, "start")
        database.status = "available"
        return database

    def find_db_from_id(self, db_id: str) -> Database:
        if self.arn_regex.match(db_id):
            arn_breakdown = db_id.split(":")
            region = arn_breakdown[3]
            backend = rds_backends[self.account_id][region]
            db_name = arn_breakdown[-1]
        else:
            backend = self
            db_name = db_id

        return backend.describe_db_instances(db_name)[0]

    def delete_db_instance(
        self, db_instance_identifier: str, db_snapshot_name: Optional[str] = None
    ) -> Database:
        self._validate_db_identifier(db_instance_identifier)
        if db_instance_identifier in self.databases:
            if self.databases[db_instance_identifier].deletion_protection:
                raise InvalidParameterValue(
                    "Can't delete Instance with protection enabled"
                )
            if db_snapshot_name:
                self.create_auto_snapshot(db_instance_identifier, db_snapshot_name)
            database = self.databases.pop(db_instance_identifier)
            if database.is_replica:
                primary = self.find_db_from_id(database.source_db_identifier)  # type: ignore
                primary.remove_replica(database)
            if database.db_cluster_identifier in self.clusters:
                self.clusters[database.db_cluster_identifier].cluster_members.remove(
                    db_instance_identifier
                )
            database.status = "deleting"
            return database
        else:
            raise DBInstanceNotFoundError(db_instance_identifier)

    def create_db_security_group(
        self, group_name: str, description: str, tags: List[Dict[str, str]]
    ) -> SecurityGroup:
        security_group = SecurityGroup(self.account_id, group_name, description, tags)
        self.security_groups[group_name] = security_group
        return security_group

    def describe_security_groups(self, security_group_name: str) -> List[SecurityGroup]:
        if security_group_name:
            if security_group_name in self.security_groups:
                return [self.security_groups[security_group_name]]
            else:
                raise DBSecurityGroupNotFoundError(security_group_name)
        return list(self.security_groups.values())

    def delete_security_group(self, security_group_name: str) -> SecurityGroup:
        if security_group_name in self.security_groups:
            return self.security_groups.pop(security_group_name)
        else:
            raise DBSecurityGroupNotFoundError(security_group_name)

    def delete_db_parameter_group(
        self, db_parameter_group_name: str
    ) -> "DBParameterGroup":
        if db_parameter_group_name in self.db_parameter_groups:
            return self.db_parameter_groups.pop(db_parameter_group_name)
        else:
            raise DBParameterGroupNotFoundError(db_parameter_group_name)

    def authorize_security_group(
        self, security_group_name: str, cidr_ip: str
    ) -> SecurityGroup:
        security_group = self.describe_security_groups(security_group_name)[0]
        security_group.authorize_cidr(cidr_ip)
        return security_group

    def create_subnet_group(
        self,
        subnet_name: str,
        description: str,
        subnets: List[Any],
        tags: List[Dict[str, str]],
    ) -> SubnetGroup:
        subnet_group = SubnetGroup(
            subnet_name, description, subnets, tags, self.region_name, self.account_id
        )
        self.subnet_groups[subnet_name] = subnet_group
        return subnet_group

    def describe_db_subnet_groups(self, subnet_group_name: str) -> List[SubnetGroup]:
        if subnet_group_name:
            if subnet_group_name in self.subnet_groups:
                return [self.subnet_groups[subnet_group_name]]
            else:
                raise DBSubnetGroupNotFoundError(subnet_group_name)
        return list(self.subnet_groups.values())

    def modify_db_subnet_group(
        self, subnet_name: str, description: str, subnets: "List[Subnet]"
    ) -> SubnetGroup:
        subnet_group = self.subnet_groups.pop(subnet_name)
        if not subnet_group:
            raise DBSubnetGroupNotFoundError(subnet_name)
        subnet_group.subnet_name = subnet_name
        subnet_group.subnets = subnets
        if description is not None:
            subnet_group.description = description
        return subnet_group

    def delete_subnet_group(self, subnet_name: str) -> SubnetGroup:
        if subnet_name in self.subnet_groups:
            return self.subnet_groups.pop(subnet_name)
        else:
            raise DBSubnetGroupNotFoundError(subnet_name)

    def create_option_group(self, option_group_kwargs: Dict[str, Any]) -> "OptionGroup":
        option_group_id = option_group_kwargs["name"]
        # This list was verified against the AWS Console on 14 Dec 2022
        # Having an automated way (using the CLI) would be nice, but AFAICS that's not possible
        #
        # Some options that are allowed in the CLI, but that do now show up in the Console:
        # - Mysql 5.5
        # - All postgres-versions
        # - oracle-se and oracle-se1 - I could not deduct the available versions
        #   `Cannot find major version 19 for oracle-se`
        #   (The engines do exist, otherwise the error would be `Invalid DB engine`
        valid_option_group_engines = {
            "mariadb": ["10.0", "10.1", "10.2", "10.3", "10.4", "10.5", "10.6"],
            "mysql": ["5.5", "5.6", "5.7", "8.0"],
            "oracle-ee": ["19"],
            "oracle-ee-cdb": ["19", "21"],
            "oracle-se": [],
            "oracle-se1": [],
            "oracle-se2": ["19"],
            "oracle-se2-cdb": ["19", "21"],
            "postgres": ["10", "11", "12", "13"],
            "sqlserver-ee": ["11.00", "12.00", "13.00", "14.00", "15.00"],
            "sqlserver-ex": ["11.00", "12.00", "13.00", "14.00", "15.00"],
            "sqlserver-se": ["11.00", "12.00", "13.00", "14.00", "15.00"],
            "sqlserver-web": ["11.00", "12.00", "13.00", "14.00", "15.00"],
        }
        if option_group_id in self.option_groups:
            raise RDSClientError(
                "OptionGroupAlreadyExistsFault",
                f"An option group named {option_group_id} already exists.",
            )
        if (
            "description" not in option_group_kwargs
            or not option_group_kwargs["description"]
        ):
            raise RDSClientError(
                "InvalidParameterValue",
                "The parameter OptionGroupDescription must be provided and must not be blank.",
            )
        if option_group_kwargs["engine_name"] not in valid_option_group_engines.keys():
            raise RDSClientError(
                "InvalidParameterValue", "Invalid DB engine: non-existent"
            )
        if (
            option_group_kwargs["major_engine_version"]  # type: ignore
            not in valid_option_group_engines[option_group_kwargs["engine_name"]]
        ):
            raise RDSClientError(
                "InvalidParameterCombination",
                f"Cannot find major version {option_group_kwargs['major_engine_version']} for {option_group_kwargs['engine_name']}",
            )
        # AWS also creates default option groups, if they do not yet exist, when creating an option group in the CLI
        # Maybe we should do the same
        # {
        #     "OptionGroupName": "default:postgres-10",
        #     "OptionGroupDescription": "Default option group for postgres 10",
        #     "EngineName": "postgres",
        #     "MajorEngineVersion": "10",
        #     "Options": [],
        #     "AllowsVpcAndNonVpcInstanceMemberships": true,
        #     "OptionGroupArn": "arn:aws:rds:us-east-1:{account}:og:default:postgres-10"
        # }
        # The CLI does not allow deletion of default groups

        # add region and account-id to construct the arn
        option_group_kwargs["region"] = self.region_name
        option_group_kwargs["account_id"] = self.account_id
        option_group = OptionGroup(**option_group_kwargs)
        self.option_groups[option_group_id] = option_group
        return option_group

    def delete_option_group(self, option_group_name: str) -> "OptionGroup":
        if option_group_name in self.option_groups:
            return self.option_groups.pop(option_group_name)
        else:
            raise OptionGroupNotFoundFaultError(option_group_name)

    def describe_option_groups(
        self, option_group_kwargs: Dict[str, Any]
    ) -> List["OptionGroup"]:
        option_group_list = []

        if option_group_kwargs["marker"]:
            marker = option_group_kwargs["marker"]
        else:
            marker = 0
        if option_group_kwargs["max_records"]:
            if (
                option_group_kwargs["max_records"] < 20
                or option_group_kwargs["max_records"] > 100
            ):
                raise RDSClientError(
                    "InvalidParameterValue",
                    "Invalid value for max records. Must be between 20 and 100",
                )
            max_records = option_group_kwargs["max_records"]
        else:
            max_records = 100

        for option_group in self.option_groups.values():
            if (
                option_group_kwargs["name"]
                and option_group.name != option_group_kwargs["name"]
            ):
                continue
            elif (
                option_group_kwargs["engine_name"]
                and option_group.engine_name != option_group_kwargs["engine_name"]
            ):
                continue
            elif (
                option_group_kwargs["major_engine_version"]
                and option_group.major_engine_version
                != option_group_kwargs["major_engine_version"]
            ):
                continue
            else:
                option_group_list.append(option_group)
        if not len(option_group_list):
            raise OptionGroupNotFoundFaultError(option_group_kwargs["name"])
        return option_group_list[marker : max_records + marker]

    @staticmethod
    def describe_option_group_options(
        engine_name: str, major_engine_version: Optional[str] = None
    ) -> str:
        default_option_group_options = {
            "mysql": {
                "5.6": '<DescribeOptionGroupOptionsResponse xmlns="http://rds.amazonaws.com/doc/2014-09-01/">\n  <DescribeOptionGroupOptionsResult>\n    <OptionGroupOptions>\n    \n      <OptionGroupOption><MajorEngineVersion>5.6</MajorEngineVersion><DefaultPort>11211</DefaultPort><PortRequired>True</PortRequired><OptionsDependedOn></OptionsDependedOn><Description>Innodb Memcached for MySQL</Description><Name>MEMCACHED</Name><OptionGroupOptionSettings><OptionGroupOptionSetting><AllowedValues>1-4294967295</AllowedValues><ApplyType>STATIC</ApplyType><DefaultValue>1</DefaultValue><IsModifiable>True</IsModifiable><SettingDescription>Specifies how many memcached read operations (get) to perform before doing a COMMIT to start a new transaction</SettingDescription><SettingName>DAEMON_MEMCACHED_R_BATCH_SIZE</SettingName></OptionGroupOptionSetting><OptionGroupOptionSetting><AllowedValues>1-4294967295</AllowedValues><ApplyType>STATIC</ApplyType><DefaultValue>1</DefaultValue><IsModifiable>True</IsModifiable><SettingDescription>Specifies how many memcached write operations, such as add, set, or incr, to perform before doing a COMMIT to start a new transaction</SettingDescription><SettingName>DAEMON_MEMCACHED_W_BATCH_SIZE</SettingName></OptionGroupOptionSetting><OptionGroupOptionSetting><AllowedValues>1-1073741824</AllowedValues><ApplyType>DYNAMIC</ApplyType><DefaultValue>5</DefaultValue><IsModifiable>True</IsModifiable><SettingDescription>Specifies how often to auto-commit idle connections that use the InnoDB memcached interface.</SettingDescription><SettingName>INNODB_API_BK_COMMIT_INTERVAL</SettingName></OptionGroupOptionSetting><OptionGroupOptionSetting><AllowedValues>0,1</AllowedValues><ApplyType>STATIC</ApplyType><DefaultValue>0</DefaultValue><IsModifiable>True</IsModifiable><SettingDescription>Disables the use of row locks when using the InnoDB memcached interface.</SettingDescription><SettingName>INNODB_API_DISABLE_ROWLOCK</SettingName></OptionGroupOptionSetting><OptionGroupOptionSetting><AllowedValues>0,1</AllowedValues><ApplyType>STATIC</ApplyType><DefaultValue>0</DefaultValue><IsModifiable>True</IsModifiable><SettingDescription>Locks the table used by the InnoDB memcached plugin, so that it cannot be dropped or altered by DDL through the SQL interface.</SettingDescription><SettingName>INNODB_API_ENABLE_MDL</SettingName></OptionGroupOptionSetting><OptionGroupOptionSetting><AllowedValues>0-3</AllowedValues><ApplyType>STATIC</ApplyType><DefaultValue>0</DefaultValue><IsModifiable>True</IsModifiable><SettingDescription>Lets you control the transaction isolation level on queries processed by the memcached interface.</SettingDescription><SettingName>INNODB_API_TRX_LEVEL</SettingName></OptionGroupOptionSetting><OptionGroupOptionSetting><AllowedValues>auto,ascii,binary</AllowedValues><ApplyType>STATIC</ApplyType><DefaultValue>auto</DefaultValue><IsModifiable>True</IsModifiable><SettingDescription>The binding protocol to use which can be either auto, ascii, or binary. The default is auto which means the server automatically negotiates the protocol with the client.</SettingDescription><SettingName>BINDING_PROTOCOL</SettingName></OptionGroupOptionSetting><OptionGroupOptionSetting><AllowedValues>1-2048</AllowedValues><ApplyType>STATIC</ApplyType><DefaultValue>1024</DefaultValue><IsModifiable>True</IsModifiable><SettingDescription>The backlog queue configures how many network connections can be waiting to be processed by memcached</SettingDescription><SettingName>BACKLOG_QUEUE_LIMIT</SettingName></OptionGroupOptionSetting><OptionGroupOptionSetting><AllowedValues>0,1</AllowedValues><ApplyType>STATIC</ApplyType><DefaultValue>0</DefaultValue><IsModifiable>True</IsModifiable><SettingDescription>Disable the use of compare and swap (CAS) which reduces the per-item size by 8 bytes.</SettingDescription><SettingName>CAS_DISABLED</SettingName></OptionGroupOptionSetting><OptionGroupOptionSetting><AllowedValues>1-48</AllowedValues><ApplyType>STATIC</ApplyType><DefaultValue>48</DefaultValue><IsModifiable>True</IsModifiable><SettingDescription>Minimum chunk size in bytes to allocate for the smallest item\'s key, value, and flags. The default is 48 and you can get a significant memory efficiency gain with a lower value.</SettingDescription><SettingName>CHUNK_SIZE</SettingName></OptionGroupOptionSetting><OptionGroupOptionSetting><AllowedValues>1-2</AllowedValues><ApplyType>STATIC</ApplyType><DefaultValue>1.25</DefaultValue><IsModifiable>True</IsModifiable><SettingDescription>Chunk size growth factor that controls the size of each successive chunk with each chunk growing times this amount larger than the previous chunk.</SettingDescription><SettingName>CHUNK_SIZE_GROWTH_FACTOR</SettingName></OptionGroupOptionSetting><OptionGroupOptionSetting><AllowedValues>0,1</AllowedValues><ApplyType>STATIC</ApplyType><DefaultValue>0</DefaultValue><IsModifiable>True</IsModifiable><SettingDescription>If enabled when there is no more memory to store items, memcached will return an error rather than evicting items.</SettingDescription><SettingName>ERROR_ON_MEMORY_EXHAUSTED</SettingName></OptionGroupOptionSetting><OptionGroupOptionSetting><AllowedValues>10-1024</AllowedValues><ApplyType>STATIC</ApplyType><DefaultValue>1024</DefaultValue><IsModifiable>True</IsModifiable><SettingDescription>Maximum number of concurrent connections. Setting this value to anything less than 10 prevents MySQL from starting.</SettingDescription><SettingName>MAX_SIMULTANEOUS_CONNECTIONS</SettingName></OptionGroupOptionSetting><OptionGroupOptionSetting><AllowedValues>v,vv,vvv</AllowedValues><ApplyType>STATIC</ApplyType><DefaultValue>v</DefaultValue><IsModifiable>True</IsModifiable><SettingDescription>Verbose level for memcached.</SettingDescription><SettingName>VERBOSITY</SettingName></OptionGroupOptionSetting></OptionGroupOptionSettings><EngineName>mysql</EngineName></OptionGroupOption>\n    \n    </OptionGroupOptions>\n  </DescribeOptionGroupOptionsResult>\n  <ResponseMetadata>\n    <RequestId>457f7bb8-9fbf-11e4-9084-5754f80d5144</RequestId>\n  </ResponseMetadata>\n</DescribeOptionGroupOptionsResponse>',
                "all": '<DescribeOptionGroupOptionsResponse xmlns="http://rds.amazonaws.com/doc/2014-09-01/">\n  <DescribeOptionGroupOptionsResult>\n    <OptionGroupOptions>\n    \n      <OptionGroupOption><MajorEngineVersion>5.6</MajorEngineVersion><DefaultPort>11211</DefaultPort><PortRequired>True</PortRequired><OptionsDependedOn></OptionsDependedOn><Description>Innodb Memcached for MySQL</Description><Name>MEMCACHED</Name><OptionGroupOptionSettings><OptionGroupOptionSetting><AllowedValues>1-4294967295</AllowedValues><ApplyType>STATIC</ApplyType><DefaultValue>1</DefaultValue><IsModifiable>True</IsModifiable><SettingDescription>Specifies how many memcached read operations (get) to perform before doing a COMMIT to start a new transaction</SettingDescription><SettingName>DAEMON_MEMCACHED_R_BATCH_SIZE</SettingName></OptionGroupOptionSetting><OptionGroupOptionSetting><AllowedValues>1-4294967295</AllowedValues><ApplyType>STATIC</ApplyType><DefaultValue>1</DefaultValue><IsModifiable>True</IsModifiable><SettingDescription>Specifies how many memcached write operations, such as add, set, or incr, to perform before doing a COMMIT to start a new transaction</SettingDescription><SettingName>DAEMON_MEMCACHED_W_BATCH_SIZE</SettingName></OptionGroupOptionSetting><OptionGroupOptionSetting><AllowedValues>1-1073741824</AllowedValues><ApplyType>DYNAMIC</ApplyType><DefaultValue>5</DefaultValue><IsModifiable>True</IsModifiable><SettingDescription>Specifies how often to auto-commit idle connections that use the InnoDB memcached interface.</SettingDescription><SettingName>INNODB_API_BK_COMMIT_INTERVAL</SettingName></OptionGroupOptionSetting><OptionGroupOptionSetting><AllowedValues>0,1</AllowedValues><ApplyType>STATIC</ApplyType><DefaultValue>0</DefaultValue><IsModifiable>True</IsModifiable><SettingDescription>Disables the use of row locks when using the InnoDB memcached interface.</SettingDescription><SettingName>INNODB_API_DISABLE_ROWLOCK</SettingName></OptionGroupOptionSetting><OptionGroupOptionSetting><AllowedValues>0,1</AllowedValues><ApplyType>STATIC</ApplyType><DefaultValue>0</DefaultValue><IsModifiable>True</IsModifiable><SettingDescription>Locks the table used by the InnoDB memcached plugin, so that it cannot be dropped or altered by DDL through the SQL interface.</SettingDescription><SettingName>INNODB_API_ENABLE_MDL</SettingName></OptionGroupOptionSetting><OptionGroupOptionSetting><AllowedValues>0-3</AllowedValues><ApplyType>STATIC</ApplyType><DefaultValue>0</DefaultValue><IsModifiable>True</IsModifiable><SettingDescription>Lets you control the transaction isolation level on queries processed by the memcached interface.</SettingDescription><SettingName>INNODB_API_TRX_LEVEL</SettingName></OptionGroupOptionSetting><OptionGroupOptionSetting><AllowedValues>auto,ascii,binary</AllowedValues><ApplyType>STATIC</ApplyType><DefaultValue>auto</DefaultValue><IsModifiable>True</IsModifiable><SettingDescription>The binding protocol to use which can be either auto, ascii, or binary. The default is auto which means the server automatically negotiates the protocol with the client.</SettingDescription><SettingName>BINDING_PROTOCOL</SettingName></OptionGroupOptionSetting><OptionGroupOptionSetting><AllowedValues>1-2048</AllowedValues><ApplyType>STATIC</ApplyType><DefaultValue>1024</DefaultValue><IsModifiable>True</IsModifiable><SettingDescription>The backlog queue configures how many network connections can be waiting to be processed by memcached</SettingDescription><SettingName>BACKLOG_QUEUE_LIMIT</SettingName></OptionGroupOptionSetting><OptionGroupOptionSetting><AllowedValues>0,1</AllowedValues><ApplyType>STATIC</ApplyType><DefaultValue>0</DefaultValue><IsModifiable>True</IsModifiable><SettingDescription>Disable the use of compare and swap (CAS) which reduces the per-item size by 8 bytes.</SettingDescription><SettingName>CAS_DISABLED</SettingName></OptionGroupOptionSetting><OptionGroupOptionSetting><AllowedValues>1-48</AllowedValues><ApplyType>STATIC</ApplyType><DefaultValue>48</DefaultValue><IsModifiable>True</IsModifiable><SettingDescription>Minimum chunk size in bytes to allocate for the smallest item\'s key, value, and flags. The default is 48 and you can get a significant memory efficiency gain with a lower value.</SettingDescription><SettingName>CHUNK_SIZE</SettingName></OptionGroupOptionSetting><OptionGroupOptionSetting><AllowedValues>1-2</AllowedValues><ApplyType>STATIC</ApplyType><DefaultValue>1.25</DefaultValue><IsModifiable>True</IsModifiable><SettingDescription>Chunk size growth factor that controls the size of each successive chunk with each chunk growing times this amount larger than the previous chunk.</SettingDescription><SettingName>CHUNK_SIZE_GROWTH_FACTOR</SettingName></OptionGroupOptionSetting><OptionGroupOptionSetting><AllowedValues>0,1</AllowedValues><ApplyType>STATIC</ApplyType><DefaultValue>0</DefaultValue><IsModifiable>True</IsModifiable><SettingDescription>If enabled when there is no more memory to store items, memcached will return an error rather than evicting items.</SettingDescription><SettingName>ERROR_ON_MEMORY_EXHAUSTED</SettingName></OptionGroupOptionSetting><OptionGroupOptionSetting><AllowedValues>10-1024</AllowedValues><ApplyType>STATIC</ApplyType><DefaultValue>1024</DefaultValue><IsModifiable>True</IsModifiable><SettingDescription>Maximum number of concurrent connections. Setting this value to anything less than 10 prevents MySQL from starting.</SettingDescription><SettingName>MAX_SIMULTANEOUS_CONNECTIONS</SettingName></OptionGroupOptionSetting><OptionGroupOptionSetting><AllowedValues>v,vv,vvv</AllowedValues><ApplyType>STATIC</ApplyType><DefaultValue>v</DefaultValue><IsModifiable>True</IsModifiable><SettingDescription>Verbose level for memcached.</SettingDescription><SettingName>VERBOSITY</SettingName></OptionGroupOptionSetting></OptionGroupOptionSettings><EngineName>mysql</EngineName></OptionGroupOption>\n    \n    </OptionGroupOptions>\n  </DescribeOptionGroupOptionsResult>\n  <ResponseMetadata>\n    <RequestId>457f7bb8-9fbf-11e4-9084-5754f80d5144</RequestId>\n  </ResponseMetadata>\n</DescribeOptionGroupOptionsResponse>',
            },
            "oracle-ee": {
                "11.2": '<DescribeOptionGroupOptionsResponse xmlns="http://rds.amazonaws.com/doc/2014-09-01/">\n  <DescribeOptionGroupOptionsResult>\n    <OptionGroupOptions>\n    \n      <OptionGroupOption><MajorEngineVersion>11.2</MajorEngineVersion><OptionsDependedOn><OptionName>XMLDB</OptionName></OptionsDependedOn><Description>Oracle Application Express Runtime Environment</Description><Name>APEX</Name><OptionGroupOptionSettings></OptionGroupOptionSettings><EngineName>oracle-ee</EngineName></OptionGroupOption>\n    \n      <OptionGroupOption><MajorEngineVersion>11.2</MajorEngineVersion><OptionsDependedOn><OptionName>APEX</OptionName></OptionsDependedOn><Description>Oracle Application Express Development Environment</Description><Name>APEX-DEV</Name><OptionGroupOptionSettings></OptionGroupOptionSettings><EngineName>oracle-ee</EngineName></OptionGroupOption>\n    \n      <OptionGroupOption><MajorEngineVersion>11.2</MajorEngineVersion><OptionsDependedOn></OptionsDependedOn><Description>Oracle Advanced Security - Native Network Encryption</Description><Name>NATIVE_NETWORK_ENCRYPTION</Name><OptionGroupOptionSettings><OptionGroupOptionSetting><AllowedValues>ACCEPTED,REJECTED,REQUESTED,REQUIRED</AllowedValues><ApplyType>STATIC</ApplyType><DefaultValue>REQUESTED</DefaultValue><IsModifiable>True</IsModifiable><SettingDescription>Specifies the desired encryption behavior</SettingDescription><SettingName>SQLNET.ENCRYPTION_SERVER</SettingName></OptionGroupOptionSetting><OptionGroupOptionSetting><AllowedValues>ACCEPTED,REJECTED,REQUESTED,REQUIRED</AllowedValues><ApplyType>STATIC</ApplyType><DefaultValue>REQUESTED</DefaultValue><IsModifiable>True</IsModifiable><SettingDescription>Specifies the desired data integrity behavior</SettingDescription><SettingName>SQLNET.CRYPTO_CHECKSUM_SERVER</SettingName></OptionGroupOptionSetting><OptionGroupOptionSetting><AllowedValues>RC4_256,AES256,AES192,3DES168,RC4_128,AES128,3DES112,RC4_56,DES,RC4_40,DES40</AllowedValues><ApplyType>STATIC</ApplyType><DefaultValue>RC4_256,AES256,AES192,3DES168,RC4_128,AES128,3DES112,RC4_56,DES,RC4_40,DES40</DefaultValue><IsModifiable>True</IsModifiable><SettingDescription>Specifies list of encryption algorithms in order of intended use</SettingDescription><SettingName>SQLNET.ENCRYPTION_TYPES_SERVER</SettingName></OptionGroupOptionSetting><OptionGroupOptionSetting><AllowedValues>SHA1,MD5</AllowedValues><ApplyType>STATIC</ApplyType><DefaultValue>SHA1,MD5</DefaultValue><IsModifiable>True</IsModifiable><SettingDescription>Specifies list of checksumming algorithms in order of intended use</SettingDescription><SettingName>SQLNET.CRYPTO_CHECKSUM_TYPES_SERVER</SettingName></OptionGroupOptionSetting></OptionGroupOptionSettings><EngineName>oracle-ee</EngineName></OptionGroupOption>\n    \n      <OptionGroupOption><MajorEngineVersion>11.2</MajorEngineVersion><DefaultPort>1158</DefaultPort><PortRequired>True</PortRequired><OptionsDependedOn></OptionsDependedOn><Description>Oracle Enterprise Manager (Database Control only)</Description><Name>OEM</Name><OptionGroupOptionSettings></OptionGroupOptionSettings><EngineName>oracle-ee</EngineName></OptionGroupOption>\n    \n      <OptionGroupOption><MajorEngineVersion>11.2</MajorEngineVersion><OptionsDependedOn></OptionsDependedOn><Description>Oracle Statspack</Description><Name>STATSPACK</Name><OptionGroupOptionSettings></OptionGroupOptionSettings><EngineName>oracle-ee</EngineName></OptionGroupOption>\n    \n      <OptionGroupOption><MajorEngineVersion>11.2</MajorEngineVersion><Persistent>True</Persistent><OptionsDependedOn></OptionsDependedOn><Permanent>True</Permanent><Description>Oracle Advanced Security - Transparent Data Encryption</Description><Name>TDE</Name><OptionGroupOptionSettings></OptionGroupOptionSettings><EngineName>oracle-ee</EngineName></OptionGroupOption>\n    \n      <OptionGroupOption><MajorEngineVersion>11.2</MajorEngineVersion><Persistent>True</Persistent><OptionsDependedOn></OptionsDependedOn><Permanent>True</Permanent><Description>Oracle Advanced Security - TDE with HSM</Description><Name>TDE_HSM</Name><OptionGroupOptionSettings></OptionGroupOptionSettings><EngineName>oracle-ee</EngineName></OptionGroupOption>\n    \n      <OptionGroupOption><MajorEngineVersion>11.2</MajorEngineVersion><Persistent>True</Persistent><OptionsDependedOn></OptionsDependedOn><Permanent>True</Permanent><Description>Change time zone</Description><Name>Timezone</Name><OptionGroupOptionSettings><OptionGroupOptionSetting><AllowedValues>Africa/Cairo,Africa/Casablanca,Africa/Harare,Africa/Monrovia,Africa/Nairobi,Africa/Tripoli,Africa/Windhoek,America/Araguaina,America/Asuncion,America/Bogota,America/Caracas,America/Chihuahua,America/Cuiaba,America/Denver,America/Fortaleza,America/Guatemala,America/Halifax,America/Manaus,America/Matamoros,America/Monterrey,America/Montevideo,America/Phoenix,America/Santiago,America/Tijuana,Asia/Amman,Asia/Ashgabat,Asia/Baghdad,Asia/Baku,Asia/Bangkok,Asia/Beirut,Asia/Calcutta,Asia/Damascus,Asia/Dhaka,Asia/Irkutsk,Asia/Jerusalem,Asia/Kabul,Asia/Karachi,Asia/Kathmandu,Asia/Krasnoyarsk,Asia/Magadan,Asia/Muscat,Asia/Novosibirsk,Asia/Riyadh,Asia/Seoul,Asia/Shanghai,Asia/Singapore,Asia/Taipei,Asia/Tehran,Asia/Tokyo,Asia/Ulaanbaatar,Asia/Vladivostok,Asia/Yakutsk,Asia/Yerevan,Atlantic/Azores,Australia/Adelaide,Australia/Brisbane,Australia/Darwin,Australia/Hobart,Australia/Perth,Australia/Sydney,Brazil/East,Canada/Newfoundland,Canada/Saskatchewan,Europe/Amsterdam,Europe/Athens,Europe/Dublin,Europe/Helsinki,Europe/Istanbul,Europe/Kaliningrad,Europe/Moscow,Europe/Paris,Europe/Prague,Europe/Sarajevo,Pacific/Auckland,Pacific/Fiji,Pacific/Guam,Pacific/Honolulu,Pacific/Samoa,US/Alaska,US/Central,US/Eastern,US/East-Indiana,US/Pacific,UTC</AllowedValues><ApplyType>DYNAMIC</ApplyType><DefaultValue>UTC</DefaultValue><IsModifiable>True</IsModifiable><SettingDescription>Specifies the timezone the user wants to change the system time to</SettingDescription><SettingName>TIME_ZONE</SettingName></OptionGroupOptionSetting></OptionGroupOptionSettings><EngineName>oracle-ee</EngineName></OptionGroupOption>\n    \n      <OptionGroupOption><MajorEngineVersion>11.2</MajorEngineVersion><OptionsDependedOn></OptionsDependedOn><Description>Oracle XMLDB Repository</Description><Name>XMLDB</Name><OptionGroupOptionSettings></OptionGroupOptionSettings><EngineName>oracle-ee</EngineName></OptionGroupOption>\n    \n    </OptionGroupOptions>\n  </DescribeOptionGroupOptionsResult>\n  <ResponseMetadata>\n    <RequestId>457f7bb8-9fbf-11e4-9084-5754f80d5144</RequestId>\n  </ResponseMetadata>\n</DescribeOptionGroupOptionsResponse>',
                "all": '<DescribeOptionGroupOptionsResponse xmlns="http://rds.amazonaws.com/doc/2014-09-01/">\n  <DescribeOptionGroupOptionsResult>\n    <OptionGroupOptions>\n    \n      <OptionGroupOption><MajorEngineVersion>11.2</MajorEngineVersion><OptionsDependedOn><OptionName>XMLDB</OptionName></OptionsDependedOn><Description>Oracle Application Express Runtime Environment</Description><Name>APEX</Name><OptionGroupOptionSettings></OptionGroupOptionSettings><EngineName>oracle-ee</EngineName></OptionGroupOption>\n    \n      <OptionGroupOption><MajorEngineVersion>11.2</MajorEngineVersion><OptionsDependedOn><OptionName>APEX</OptionName></OptionsDependedOn><Description>Oracle Application Express Development Environment</Description><Name>APEX-DEV</Name><OptionGroupOptionSettings></OptionGroupOptionSettings><EngineName>oracle-ee</EngineName></OptionGroupOption>\n    \n      <OptionGroupOption><MajorEngineVersion>11.2</MajorEngineVersion><OptionsDependedOn></OptionsDependedOn><Description>Oracle Advanced Security - Native Network Encryption</Description><Name>NATIVE_NETWORK_ENCRYPTION</Name><OptionGroupOptionSettings><OptionGroupOptionSetting><AllowedValues>ACCEPTED,REJECTED,REQUESTED,REQUIRED</AllowedValues><ApplyType>STATIC</ApplyType><DefaultValue>REQUESTED</DefaultValue><IsModifiable>True</IsModifiable><SettingDescription>Specifies the desired encryption behavior</SettingDescription><SettingName>SQLNET.ENCRYPTION_SERVER</SettingName></OptionGroupOptionSetting><OptionGroupOptionSetting><AllowedValues>ACCEPTED,REJECTED,REQUESTED,REQUIRED</AllowedValues><ApplyType>STATIC</ApplyType><DefaultValue>REQUESTED</DefaultValue><IsModifiable>True</IsModifiable><SettingDescription>Specifies the desired data integrity behavior</SettingDescription><SettingName>SQLNET.CRYPTO_CHECKSUM_SERVER</SettingName></OptionGroupOptionSetting><OptionGroupOptionSetting><AllowedValues>RC4_256,AES256,AES192,3DES168,RC4_128,AES128,3DES112,RC4_56,DES,RC4_40,DES40</AllowedValues><ApplyType>STATIC</ApplyType><DefaultValue>RC4_256,AES256,AES192,3DES168,RC4_128,AES128,3DES112,RC4_56,DES,RC4_40,DES40</DefaultValue><IsModifiable>True</IsModifiable><SettingDescription>Specifies list of encryption algorithms in order of intended use</SettingDescription><SettingName>SQLNET.ENCRYPTION_TYPES_SERVER</SettingName></OptionGroupOptionSetting><OptionGroupOptionSetting><AllowedValues>SHA1,MD5</AllowedValues><ApplyType>STATIC</ApplyType><DefaultValue>SHA1,MD5</DefaultValue><IsModifiable>True</IsModifiable><SettingDescription>Specifies list of checksumming algorithms in order of intended use</SettingDescription><SettingName>SQLNET.CRYPTO_CHECKSUM_TYPES_SERVER</SettingName></OptionGroupOptionSetting></OptionGroupOptionSettings><EngineName>oracle-ee</EngineName></OptionGroupOption>\n    \n      <OptionGroupOption><MajorEngineVersion>11.2</MajorEngineVersion><DefaultPort>1158</DefaultPort><PortRequired>True</PortRequired><OptionsDependedOn></OptionsDependedOn><Description>Oracle Enterprise Manager (Database Control only)</Description><Name>OEM</Name><OptionGroupOptionSettings></OptionGroupOptionSettings><EngineName>oracle-ee</EngineName></OptionGroupOption>\n    \n      <OptionGroupOption><MajorEngineVersion>11.2</MajorEngineVersion><OptionsDependedOn></OptionsDependedOn><Description>Oracle Statspack</Description><Name>STATSPACK</Name><OptionGroupOptionSettings></OptionGroupOptionSettings><EngineName>oracle-ee</EngineName></OptionGroupOption>\n    \n      <OptionGroupOption><MajorEngineVersion>11.2</MajorEngineVersion><Persistent>True</Persistent><OptionsDependedOn></OptionsDependedOn><Permanent>True</Permanent><Description>Oracle Advanced Security - Transparent Data Encryption</Description><Name>TDE</Name><OptionGroupOptionSettings></OptionGroupOptionSettings><EngineName>oracle-ee</EngineName></OptionGroupOption>\n    \n      <OptionGroupOption><MajorEngineVersion>11.2</MajorEngineVersion><Persistent>True</Persistent><OptionsDependedOn></OptionsDependedOn><Permanent>True</Permanent><Description>Oracle Advanced Security - TDE with HSM</Description><Name>TDE_HSM</Name><OptionGroupOptionSettings></OptionGroupOptionSettings><EngineName>oracle-ee</EngineName></OptionGroupOption>\n    \n      <OptionGroupOption><MajorEngineVersion>11.2</MajorEngineVersion><Persistent>True</Persistent><OptionsDependedOn></OptionsDependedOn><Permanent>True</Permanent><Description>Change time zone</Description><Name>Timezone</Name><OptionGroupOptionSettings><OptionGroupOptionSetting><AllowedValues>Africa/Cairo,Africa/Casablanca,Africa/Harare,Africa/Monrovia,Africa/Nairobi,Africa/Tripoli,Africa/Windhoek,America/Araguaina,America/Asuncion,America/Bogota,America/Caracas,America/Chihuahua,America/Cuiaba,America/Denver,America/Fortaleza,America/Guatemala,America/Halifax,America/Manaus,America/Matamoros,America/Monterrey,America/Montevideo,America/Phoenix,America/Santiago,America/Tijuana,Asia/Amman,Asia/Ashgabat,Asia/Baghdad,Asia/Baku,Asia/Bangkok,Asia/Beirut,Asia/Calcutta,Asia/Damascus,Asia/Dhaka,Asia/Irkutsk,Asia/Jerusalem,Asia/Kabul,Asia/Karachi,Asia/Kathmandu,Asia/Krasnoyarsk,Asia/Magadan,Asia/Muscat,Asia/Novosibirsk,Asia/Riyadh,Asia/Seoul,Asia/Shanghai,Asia/Singapore,Asia/Taipei,Asia/Tehran,Asia/Tokyo,Asia/Ulaanbaatar,Asia/Vladivostok,Asia/Yakutsk,Asia/Yerevan,Atlantic/Azores,Australia/Adelaide,Australia/Brisbane,Australia/Darwin,Australia/Hobart,Australia/Perth,Australia/Sydney,Brazil/East,Canada/Newfoundland,Canada/Saskatchewan,Europe/Amsterdam,Europe/Athens,Europe/Dublin,Europe/Helsinki,Europe/Istanbul,Europe/Kaliningrad,Europe/Moscow,Europe/Paris,Europe/Prague,Europe/Sarajevo,Pacific/Auckland,Pacific/Fiji,Pacific/Guam,Pacific/Honolulu,Pacific/Samoa,US/Alaska,US/Central,US/Eastern,US/East-Indiana,US/Pacific,UTC</AllowedValues><ApplyType>DYNAMIC</ApplyType><DefaultValue>UTC</DefaultValue><IsModifiable>True</IsModifiable><SettingDescription>Specifies the timezone the user wants to change the system time to</SettingDescription><SettingName>TIME_ZONE</SettingName></OptionGroupOptionSetting></OptionGroupOptionSettings><EngineName>oracle-ee</EngineName></OptionGroupOption>\n    \n      <OptionGroupOption><MajorEngineVersion>11.2</MajorEngineVersion><OptionsDependedOn></OptionsDependedOn><Description>Oracle XMLDB Repository</Description><Name>XMLDB</Name><OptionGroupOptionSettings></OptionGroupOptionSettings><EngineName>oracle-ee</EngineName></OptionGroupOption>\n    \n    </OptionGroupOptions>\n  </DescribeOptionGroupOptionsResult>\n  <ResponseMetadata>\n    <RequestId>457f7bb8-9fbf-11e4-9084-5754f80d5144</RequestId>\n  </ResponseMetadata>\n</DescribeOptionGroupOptionsResponse>',
            },
            "oracle-sa": {
                "11.2": '<DescribeOptionGroupOptionsResponse xmlns="http://rds.amazonaws.com/doc/2014-09-01/">\n  <DescribeOptionGroupOptionsResult>\n    <OptionGroupOptions>\n    \n      <OptionGroupOption><MajorEngineVersion>11.2</MajorEngineVersion><OptionsDependedOn><OptionName>XMLDB</OptionName></OptionsDependedOn><Description>Oracle Application Express Runtime Environment</Description><Name>APEX</Name><OptionGroupOptionSettings></OptionGroupOptionSettings><EngineName>oracle-ee</EngineName></OptionGroupOption>\n    \n      <OptionGroupOption><MajorEngineVersion>11.2</MajorEngineVersion><OptionsDependedOn><OptionName>APEX</OptionName></OptionsDependedOn><Description>Oracle Application Express Development Environment</Description><Name>APEX-DEV</Name><OptionGroupOptionSettings></OptionGroupOptionSettings><EngineName>oracle-ee</EngineName></OptionGroupOption>\n    \n      <OptionGroupOption><MajorEngineVersion>11.2</MajorEngineVersion><OptionsDependedOn></OptionsDependedOn><Description>Oracle Advanced Security - Native Network Encryption</Description><Name>NATIVE_NETWORK_ENCRYPTION</Name><OptionGroupOptionSettings><OptionGroupOptionSetting><AllowedValues>ACCEPTED,REJECTED,REQUESTED,REQUIRED</AllowedValues><ApplyType>STATIC</ApplyType><DefaultValue>REQUESTED</DefaultValue><IsModifiable>True</IsModifiable><SettingDescription>Specifies the desired encryption behavior</SettingDescription><SettingName>SQLNET.ENCRYPTION_SERVER</SettingName></OptionGroupOptionSetting><OptionGroupOptionSetting><AllowedValues>ACCEPTED,REJECTED,REQUESTED,REQUIRED</AllowedValues><ApplyType>STATIC</ApplyType><DefaultValue>REQUESTED</DefaultValue><IsModifiable>True</IsModifiable><SettingDescription>Specifies the desired data integrity behavior</SettingDescription><SettingName>SQLNET.CRYPTO_CHECKSUM_SERVER</SettingName></OptionGroupOptionSetting><OptionGroupOptionSetting><AllowedValues>RC4_256,AES256,AES192,3DES168,RC4_128,AES128,3DES112,RC4_56,DES,RC4_40,DES40</AllowedValues><ApplyType>STATIC</ApplyType><DefaultValue>RC4_256,AES256,AES192,3DES168,RC4_128,AES128,3DES112,RC4_56,DES,RC4_40,DES40</DefaultValue><IsModifiable>True</IsModifiable><SettingDescription>Specifies list of encryption algorithms in order of intended use</SettingDescription><SettingName>SQLNET.ENCRYPTION_TYPES_SERVER</SettingName></OptionGroupOptionSetting><OptionGroupOptionSetting><AllowedValues>SHA1,MD5</AllowedValues><ApplyType>STATIC</ApplyType><DefaultValue>SHA1,MD5</DefaultValue><IsModifiable>True</IsModifiable><SettingDescription>Specifies list of checksumming algorithms in order of intended use</SettingDescription><SettingName>SQLNET.CRYPTO_CHECKSUM_TYPES_SERVER</SettingName></OptionGroupOptionSetting></OptionGroupOptionSettings><EngineName>oracle-ee</EngineName></OptionGroupOption>\n    \n      <OptionGroupOption><MajorEngineVersion>11.2</MajorEngineVersion><DefaultPort>1158</DefaultPort><PortRequired>True</PortRequired><OptionsDependedOn></OptionsDependedOn><Description>Oracle Enterprise Manager (Database Control only)</Description><Name>OEM</Name><OptionGroupOptionSettings></OptionGroupOptionSettings><EngineName>oracle-ee</EngineName></OptionGroupOption>\n    \n      <OptionGroupOption><MajorEngineVersion>11.2</MajorEngineVersion><OptionsDependedOn></OptionsDependedOn><Description>Oracle Statspack</Description><Name>STATSPACK</Name><OptionGroupOptionSettings></OptionGroupOptionSettings><EngineName>oracle-ee</EngineName></OptionGroupOption>\n    \n      <OptionGroupOption><MajorEngineVersion>11.2</MajorEngineVersion><Persistent>True</Persistent><OptionsDependedOn></OptionsDependedOn><Permanent>True</Permanent><Description>Oracle Advanced Security - Transparent Data Encryption</Description><Name>TDE</Name><OptionGroupOptionSettings></OptionGroupOptionSettings><EngineName>oracle-ee</EngineName></OptionGroupOption>\n    \n      <OptionGroupOption><MajorEngineVersion>11.2</MajorEngineVersion><Persistent>True</Persistent><OptionsDependedOn></OptionsDependedOn><Permanent>True</Permanent><Description>Oracle Advanced Security - TDE with HSM</Description><Name>TDE_HSM</Name><OptionGroupOptionSettings></OptionGroupOptionSettings><EngineName>oracle-ee</EngineName></OptionGroupOption>\n    \n      <OptionGroupOption><MajorEngineVersion>11.2</MajorEngineVersion><Persistent>True</Persistent><OptionsDependedOn></OptionsDependedOn><Permanent>True</Permanent><Description>Change time zone</Description><Name>Timezone</Name><OptionGroupOptionSettings><OptionGroupOptionSetting><AllowedValues>Africa/Cairo,Africa/Casablanca,Africa/Harare,Africa/Monrovia,Africa/Nairobi,Africa/Tripoli,Africa/Windhoek,America/Araguaina,America/Asuncion,America/Bogota,America/Caracas,America/Chihuahua,America/Cuiaba,America/Denver,America/Fortaleza,America/Guatemala,America/Halifax,America/Manaus,America/Matamoros,America/Monterrey,America/Montevideo,America/Phoenix,America/Santiago,America/Tijuana,Asia/Amman,Asia/Ashgabat,Asia/Baghdad,Asia/Baku,Asia/Bangkok,Asia/Beirut,Asia/Calcutta,Asia/Damascus,Asia/Dhaka,Asia/Irkutsk,Asia/Jerusalem,Asia/Kabul,Asia/Karachi,Asia/Kathmandu,Asia/Krasnoyarsk,Asia/Magadan,Asia/Muscat,Asia/Novosibirsk,Asia/Riyadh,Asia/Seoul,Asia/Shanghai,Asia/Singapore,Asia/Taipei,Asia/Tehran,Asia/Tokyo,Asia/Ulaanbaatar,Asia/Vladivostok,Asia/Yakutsk,Asia/Yerevan,Atlantic/Azores,Australia/Adelaide,Australia/Brisbane,Australia/Darwin,Australia/Hobart,Australia/Perth,Australia/Sydney,Brazil/East,Canada/Newfoundland,Canada/Saskatchewan,Europe/Amsterdam,Europe/Athens,Europe/Dublin,Europe/Helsinki,Europe/Istanbul,Europe/Kaliningrad,Europe/Moscow,Europe/Paris,Europe/Prague,Europe/Sarajevo,Pacific/Auckland,Pacific/Fiji,Pacific/Guam,Pacific/Honolulu,Pacific/Samoa,US/Alaska,US/Central,US/Eastern,US/East-Indiana,US/Pacific,UTC</AllowedValues><ApplyType>DYNAMIC</ApplyType><DefaultValue>UTC</DefaultValue><IsModifiable>True</IsModifiable><SettingDescription>Specifies the timezone the user wants to change the system time to</SettingDescription><SettingName>TIME_ZONE</SettingName></OptionGroupOptionSetting></OptionGroupOptionSettings><EngineName>oracle-ee</EngineName></OptionGroupOption>\n    \n      <OptionGroupOption><MajorEngineVersion>11.2</MajorEngineVersion><OptionsDependedOn></OptionsDependedOn><Description>Oracle XMLDB Repository</Description><Name>XMLDB</Name><OptionGroupOptionSettings></OptionGroupOptionSettings><EngineName>oracle-ee</EngineName></OptionGroupOption>\n    \n    </OptionGroupOptions>\n  </DescribeOptionGroupOptionsResult>\n  <ResponseMetadata>\n    <RequestId>457f7bb8-9fbf-11e4-9084-5754f80d5144</RequestId>\n  </ResponseMetadata>\n</DescribeOptionGroupOptionsResponse>',
                "all": '<DescribeOptionGroupOptionsResponse xmlns="http://rds.amazonaws.com/doc/2014-09-01/">\n  <DescribeOptionGroupOptionsResult>\n    <OptionGroupOptions>\n    \n      <OptionGroupOption><MajorEngineVersion>11.2</MajorEngineVersion><OptionsDependedOn><OptionName>XMLDB</OptionName></OptionsDependedOn><Description>Oracle Application Express Runtime Environment</Description><Name>APEX</Name><OptionGroupOptionSettings></OptionGroupOptionSettings><EngineName>oracle-ee</EngineName></OptionGroupOption>\n    \n      <OptionGroupOption><MajorEngineVersion>11.2</MajorEngineVersion><OptionsDependedOn><OptionName>APEX</OptionName></OptionsDependedOn><Description>Oracle Application Express Development Environment</Description><Name>APEX-DEV</Name><OptionGroupOptionSettings></OptionGroupOptionSettings><EngineName>oracle-ee</EngineName></OptionGroupOption>\n    \n      <OptionGroupOption><MajorEngineVersion>11.2</MajorEngineVersion><OptionsDependedOn></OptionsDependedOn><Description>Oracle Advanced Security - Native Network Encryption</Description><Name>NATIVE_NETWORK_ENCRYPTION</Name><OptionGroupOptionSettings><OptionGroupOptionSetting><AllowedValues>ACCEPTED,REJECTED,REQUESTED,REQUIRED</AllowedValues><ApplyType>STATIC</ApplyType><DefaultValue>REQUESTED</DefaultValue><IsModifiable>True</IsModifiable><SettingDescription>Specifies the desired encryption behavior</SettingDescription><SettingName>SQLNET.ENCRYPTION_SERVER</SettingName></OptionGroupOptionSetting><OptionGroupOptionSetting><AllowedValues>ACCEPTED,REJECTED,REQUESTED,REQUIRED</AllowedValues><ApplyType>STATIC</ApplyType><DefaultValue>REQUESTED</DefaultValue><IsModifiable>True</IsModifiable><SettingDescription>Specifies the desired data integrity behavior</SettingDescription><SettingName>SQLNET.CRYPTO_CHECKSUM_SERVER</SettingName></OptionGroupOptionSetting><OptionGroupOptionSetting><AllowedValues>RC4_256,AES256,AES192,3DES168,RC4_128,AES128,3DES112,RC4_56,DES,RC4_40,DES40</AllowedValues><ApplyType>STATIC</ApplyType><DefaultValue>RC4_256,AES256,AES192,3DES168,RC4_128,AES128,3DES112,RC4_56,DES,RC4_40,DES40</DefaultValue><IsModifiable>True</IsModifiable><SettingDescription>Specifies list of encryption algorithms in order of intended use</SettingDescription><SettingName>SQLNET.ENCRYPTION_TYPES_SERVER</SettingName></OptionGroupOptionSetting><OptionGroupOptionSetting><AllowedValues>SHA1,MD5</AllowedValues><ApplyType>STATIC</ApplyType><DefaultValue>SHA1,MD5</DefaultValue><IsModifiable>True</IsModifiable><SettingDescription>Specifies list of checksumming algorithms in order of intended use</SettingDescription><SettingName>SQLNET.CRYPTO_CHECKSUM_TYPES_SERVER</SettingName></OptionGroupOptionSetting></OptionGroupOptionSettings><EngineName>oracle-ee</EngineName></OptionGroupOption>\n    \n      <OptionGroupOption><MajorEngineVersion>11.2</MajorEngineVersion><DefaultPort>1158</DefaultPort><PortRequired>True</PortRequired><OptionsDependedOn></OptionsDependedOn><Description>Oracle Enterprise Manager (Database Control only)</Description><Name>OEM</Name><OptionGroupOptionSettings></OptionGroupOptionSettings><EngineName>oracle-ee</EngineName></OptionGroupOption>\n    \n      <OptionGroupOption><MajorEngineVersion>11.2</MajorEngineVersion><OptionsDependedOn></OptionsDependedOn><Description>Oracle Statspack</Description><Name>STATSPACK</Name><OptionGroupOptionSettings></OptionGroupOptionSettings><EngineName>oracle-ee</EngineName></OptionGroupOption>\n    \n      <OptionGroupOption><MajorEngineVersion>11.2</MajorEngineVersion><Persistent>True</Persistent><OptionsDependedOn></OptionsDependedOn><Permanent>True</Permanent><Description>Oracle Advanced Security - Transparent Data Encryption</Description><Name>TDE</Name><OptionGroupOptionSettings></OptionGroupOptionSettings><EngineName>oracle-ee</EngineName></OptionGroupOption>\n    \n      <OptionGroupOption><MajorEngineVersion>11.2</MajorEngineVersion><Persistent>True</Persistent><OptionsDependedOn></OptionsDependedOn><Permanent>True</Permanent><Description>Oracle Advanced Security - TDE with HSM</Description><Name>TDE_HSM</Name><OptionGroupOptionSettings></OptionGroupOptionSettings><EngineName>oracle-ee</EngineName></OptionGroupOption>\n    \n      <OptionGroupOption><MajorEngineVersion>11.2</MajorEngineVersion><Persistent>True</Persistent><OptionsDependedOn></OptionsDependedOn><Permanent>True</Permanent><Description>Change time zone</Description><Name>Timezone</Name><OptionGroupOptionSettings><OptionGroupOptionSetting><AllowedValues>Africa/Cairo,Africa/Casablanca,Africa/Harare,Africa/Monrovia,Africa/Nairobi,Africa/Tripoli,Africa/Windhoek,America/Araguaina,America/Asuncion,America/Bogota,America/Caracas,America/Chihuahua,America/Cuiaba,America/Denver,America/Fortaleza,America/Guatemala,America/Halifax,America/Manaus,America/Matamoros,America/Monterrey,America/Montevideo,America/Phoenix,America/Santiago,America/Tijuana,Asia/Amman,Asia/Ashgabat,Asia/Baghdad,Asia/Baku,Asia/Bangkok,Asia/Beirut,Asia/Calcutta,Asia/Damascus,Asia/Dhaka,Asia/Irkutsk,Asia/Jerusalem,Asia/Kabul,Asia/Karachi,Asia/Kathmandu,Asia/Krasnoyarsk,Asia/Magadan,Asia/Muscat,Asia/Novosibirsk,Asia/Riyadh,Asia/Seoul,Asia/Shanghai,Asia/Singapore,Asia/Taipei,Asia/Tehran,Asia/Tokyo,Asia/Ulaanbaatar,Asia/Vladivostok,Asia/Yakutsk,Asia/Yerevan,Atlantic/Azores,Australia/Adelaide,Australia/Brisbane,Australia/Darwin,Australia/Hobart,Australia/Perth,Australia/Sydney,Brazil/East,Canada/Newfoundland,Canada/Saskatchewan,Europe/Amsterdam,Europe/Athens,Europe/Dublin,Europe/Helsinki,Europe/Istanbul,Europe/Kaliningrad,Europe/Moscow,Europe/Paris,Europe/Prague,Europe/Sarajevo,Pacific/Auckland,Pacific/Fiji,Pacific/Guam,Pacific/Honolulu,Pacific/Samoa,US/Alaska,US/Central,US/Eastern,US/East-Indiana,US/Pacific,UTC</AllowedValues><ApplyType>DYNAMIC</ApplyType><DefaultValue>UTC</DefaultValue><IsModifiable>True</IsModifiable><SettingDescription>Specifies the timezone the user wants to change the system time to</SettingDescription><SettingName>TIME_ZONE</SettingName></OptionGroupOptionSetting></OptionGroupOptionSettings><EngineName>oracle-ee</EngineName></OptionGroupOption>\n    \n      <OptionGroupOption><MajorEngineVersion>11.2</MajorEngineVersion><OptionsDependedOn></OptionsDependedOn><Description>Oracle XMLDB Repository</Description><Name>XMLDB</Name><OptionGroupOptionSettings></OptionGroupOptionSettings><EngineName>oracle-ee</EngineName></OptionGroupOption>\n    \n    </OptionGroupOptions>\n  </DescribeOptionGroupOptionsResult>\n  <ResponseMetadata>\n    <RequestId>457f7bb8-9fbf-11e4-9084-5754f80d5144</RequestId>\n  </ResponseMetadata>\n</DescribeOptionGroupOptionsResponse>',
            },
            "oracle-sa1": {
                "11.2": '<DescribeOptionGroupOptionsResponse xmlns="http://rds.amazonaws.com/doc/2014-09-01/">\n  <DescribeOptionGroupOptionsResult>\n    <OptionGroupOptions>\n    \n      <OptionGroupOption><MajorEngineVersion>11.2</MajorEngineVersion><OptionsDependedOn><OptionName>XMLDB</OptionName></OptionsDependedOn><Description>Oracle Application Express Runtime Environment</Description><Name>APEX</Name><OptionGroupOptionSettings></OptionGroupOptionSettings><EngineName>oracle-ee</EngineName></OptionGroupOption>\n    \n      <OptionGroupOption><MajorEngineVersion>11.2</MajorEngineVersion><OptionsDependedOn><OptionName>APEX</OptionName></OptionsDependedOn><Description>Oracle Application Express Development Environment</Description><Name>APEX-DEV</Name><OptionGroupOptionSettings></OptionGroupOptionSettings><EngineName>oracle-ee</EngineName></OptionGroupOption>\n    \n      <OptionGroupOption><MajorEngineVersion>11.2</MajorEngineVersion><OptionsDependedOn></OptionsDependedOn><Description>Oracle Advanced Security - Native Network Encryption</Description><Name>NATIVE_NETWORK_ENCRYPTION</Name><OptionGroupOptionSettings><OptionGroupOptionSetting><AllowedValues>ACCEPTED,REJECTED,REQUESTED,REQUIRED</AllowedValues><ApplyType>STATIC</ApplyType><DefaultValue>REQUESTED</DefaultValue><IsModifiable>True</IsModifiable><SettingDescription>Specifies the desired encryption behavior</SettingDescription><SettingName>SQLNET.ENCRYPTION_SERVER</SettingName></OptionGroupOptionSetting><OptionGroupOptionSetting><AllowedValues>ACCEPTED,REJECTED,REQUESTED,REQUIRED</AllowedValues><ApplyType>STATIC</ApplyType><DefaultValue>REQUESTED</DefaultValue><IsModifiable>True</IsModifiable><SettingDescription>Specifies the desired data integrity behavior</SettingDescription><SettingName>SQLNET.CRYPTO_CHECKSUM_SERVER</SettingName></OptionGroupOptionSetting><OptionGroupOptionSetting><AllowedValues>RC4_256,AES256,AES192,3DES168,RC4_128,AES128,3DES112,RC4_56,DES,RC4_40,DES40</AllowedValues><ApplyType>STATIC</ApplyType><DefaultValue>RC4_256,AES256,AES192,3DES168,RC4_128,AES128,3DES112,RC4_56,DES,RC4_40,DES40</DefaultValue><IsModifiable>True</IsModifiable><SettingDescription>Specifies list of encryption algorithms in order of intended use</SettingDescription><SettingName>SQLNET.ENCRYPTION_TYPES_SERVER</SettingName></OptionGroupOptionSetting><OptionGroupOptionSetting><AllowedValues>SHA1,MD5</AllowedValues><ApplyType>STATIC</ApplyType><DefaultValue>SHA1,MD5</DefaultValue><IsModifiable>True</IsModifiable><SettingDescription>Specifies list of checksumming algorithms in order of intended use</SettingDescription><SettingName>SQLNET.CRYPTO_CHECKSUM_TYPES_SERVER</SettingName></OptionGroupOptionSetting></OptionGroupOptionSettings><EngineName>oracle-ee</EngineName></OptionGroupOption>\n    \n      <OptionGroupOption><MajorEngineVersion>11.2</MajorEngineVersion><DefaultPort>1158</DefaultPort><PortRequired>True</PortRequired><OptionsDependedOn></OptionsDependedOn><Description>Oracle Enterprise Manager (Database Control only)</Description><Name>OEM</Name><OptionGroupOptionSettings></OptionGroupOptionSettings><EngineName>oracle-ee</EngineName></OptionGroupOption>\n    \n      <OptionGroupOption><MajorEngineVersion>11.2</MajorEngineVersion><OptionsDependedOn></OptionsDependedOn><Description>Oracle Statspack</Description><Name>STATSPACK</Name><OptionGroupOptionSettings></OptionGroupOptionSettings><EngineName>oracle-ee</EngineName></OptionGroupOption>\n    \n      <OptionGroupOption><MajorEngineVersion>11.2</MajorEngineVersion><Persistent>True</Persistent><OptionsDependedOn></OptionsDependedOn><Permanent>True</Permanent><Description>Oracle Advanced Security - Transparent Data Encryption</Description><Name>TDE</Name><OptionGroupOptionSettings></OptionGroupOptionSettings><EngineName>oracle-ee</EngineName></OptionGroupOption>\n    \n      <OptionGroupOption><MajorEngineVersion>11.2</MajorEngineVersion><Persistent>True</Persistent><OptionsDependedOn></OptionsDependedOn><Permanent>True</Permanent><Description>Oracle Advanced Security - TDE with HSM</Description><Name>TDE_HSM</Name><OptionGroupOptionSettings></OptionGroupOptionSettings><EngineName>oracle-ee</EngineName></OptionGroupOption>\n    \n      <OptionGroupOption><MajorEngineVersion>11.2</MajorEngineVersion><Persistent>True</Persistent><OptionsDependedOn></OptionsDependedOn><Permanent>True</Permanent><Description>Change time zone</Description><Name>Timezone</Name><OptionGroupOptionSettings><OptionGroupOptionSetting><AllowedValues>Africa/Cairo,Africa/Casablanca,Africa/Harare,Africa/Monrovia,Africa/Nairobi,Africa/Tripoli,Africa/Windhoek,America/Araguaina,America/Asuncion,America/Bogota,America/Caracas,America/Chihuahua,America/Cuiaba,America/Denver,America/Fortaleza,America/Guatemala,America/Halifax,America/Manaus,America/Matamoros,America/Monterrey,America/Montevideo,America/Phoenix,America/Santiago,America/Tijuana,Asia/Amman,Asia/Ashgabat,Asia/Baghdad,Asia/Baku,Asia/Bangkok,Asia/Beirut,Asia/Calcutta,Asia/Damascus,Asia/Dhaka,Asia/Irkutsk,Asia/Jerusalem,Asia/Kabul,Asia/Karachi,Asia/Kathmandu,Asia/Krasnoyarsk,Asia/Magadan,Asia/Muscat,Asia/Novosibirsk,Asia/Riyadh,Asia/Seoul,Asia/Shanghai,Asia/Singapore,Asia/Taipei,Asia/Tehran,Asia/Tokyo,Asia/Ulaanbaatar,Asia/Vladivostok,Asia/Yakutsk,Asia/Yerevan,Atlantic/Azores,Australia/Adelaide,Australia/Brisbane,Australia/Darwin,Australia/Hobart,Australia/Perth,Australia/Sydney,Brazil/East,Canada/Newfoundland,Canada/Saskatchewan,Europe/Amsterdam,Europe/Athens,Europe/Dublin,Europe/Helsinki,Europe/Istanbul,Europe/Kaliningrad,Europe/Moscow,Europe/Paris,Europe/Prague,Europe/Sarajevo,Pacific/Auckland,Pacific/Fiji,Pacific/Guam,Pacific/Honolulu,Pacific/Samoa,US/Alaska,US/Central,US/Eastern,US/East-Indiana,US/Pacific,UTC</AllowedValues><ApplyType>DYNAMIC</ApplyType><DefaultValue>UTC</DefaultValue><IsModifiable>True</IsModifiable><SettingDescription>Specifies the timezone the user wants to change the system time to</SettingDescription><SettingName>TIME_ZONE</SettingName></OptionGroupOptionSetting></OptionGroupOptionSettings><EngineName>oracle-ee</EngineName></OptionGroupOption>\n    \n      <OptionGroupOption><MajorEngineVersion>11.2</MajorEngineVersion><OptionsDependedOn></OptionsDependedOn><Description>Oracle XMLDB Repository</Description><Name>XMLDB</Name><OptionGroupOptionSettings></OptionGroupOptionSettings><EngineName>oracle-ee</EngineName></OptionGroupOption>\n    \n    </OptionGroupOptions>\n  </DescribeOptionGroupOptionsResult>\n  <ResponseMetadata>\n    <RequestId>457f7bb8-9fbf-11e4-9084-5754f80d5144</RequestId>\n  </ResponseMetadata>\n</DescribeOptionGroupOptionsResponse>',
                "all": '<DescribeOptionGroupOptionsResponse xmlns="http://rds.amazonaws.com/doc/2014-09-01/">\n  <DescribeOptionGroupOptionsResult>\n    <OptionGroupOptions>\n    \n      <OptionGroupOption><MajorEngineVersion>11.2</MajorEngineVersion><OptionsDependedOn><OptionName>XMLDB</OptionName></OptionsDependedOn><Description>Oracle Application Express Runtime Environment</Description><Name>APEX</Name><OptionGroupOptionSettings></OptionGroupOptionSettings><EngineName>oracle-ee</EngineName></OptionGroupOption>\n    \n      <OptionGroupOption><MajorEngineVersion>11.2</MajorEngineVersion><OptionsDependedOn><OptionName>APEX</OptionName></OptionsDependedOn><Description>Oracle Application Express Development Environment</Description><Name>APEX-DEV</Name><OptionGroupOptionSettings></OptionGroupOptionSettings><EngineName>oracle-ee</EngineName></OptionGroupOption>\n    \n      <OptionGroupOption><MajorEngineVersion>11.2</MajorEngineVersion><OptionsDependedOn></OptionsDependedOn><Description>Oracle Advanced Security - Native Network Encryption</Description><Name>NATIVE_NETWORK_ENCRYPTION</Name><OptionGroupOptionSettings><OptionGroupOptionSetting><AllowedValues>ACCEPTED,REJECTED,REQUESTED,REQUIRED</AllowedValues><ApplyType>STATIC</ApplyType><DefaultValue>REQUESTED</DefaultValue><IsModifiable>True</IsModifiable><SettingDescription>Specifies the desired encryption behavior</SettingDescription><SettingName>SQLNET.ENCRYPTION_SERVER</SettingName></OptionGroupOptionSetting><OptionGroupOptionSetting><AllowedValues>ACCEPTED,REJECTED,REQUESTED,REQUIRED</AllowedValues><ApplyType>STATIC</ApplyType><DefaultValue>REQUESTED</DefaultValue><IsModifiable>True</IsModifiable><SettingDescription>Specifies the desired data integrity behavior</SettingDescription><SettingName>SQLNET.CRYPTO_CHECKSUM_SERVER</SettingName></OptionGroupOptionSetting><OptionGroupOptionSetting><AllowedValues>RC4_256,AES256,AES192,3DES168,RC4_128,AES128,3DES112,RC4_56,DES,RC4_40,DES40</AllowedValues><ApplyType>STATIC</ApplyType><DefaultValue>RC4_256,AES256,AES192,3DES168,RC4_128,AES128,3DES112,RC4_56,DES,RC4_40,DES40</DefaultValue><IsModifiable>True</IsModifiable><SettingDescription>Specifies list of encryption algorithms in order of intended use</SettingDescription><SettingName>SQLNET.ENCRYPTION_TYPES_SERVER</SettingName></OptionGroupOptionSetting><OptionGroupOptionSetting><AllowedValues>SHA1,MD5</AllowedValues><ApplyType>STATIC</ApplyType><DefaultValue>SHA1,MD5</DefaultValue><IsModifiable>True</IsModifiable><SettingDescription>Specifies list of checksumming algorithms in order of intended use</SettingDescription><SettingName>SQLNET.CRYPTO_CHECKSUM_TYPES_SERVER</SettingName></OptionGroupOptionSetting></OptionGroupOptionSettings><EngineName>oracle-ee</EngineName></OptionGroupOption>\n    \n      <OptionGroupOption><MajorEngineVersion>11.2</MajorEngineVersion><DefaultPort>1158</DefaultPort><PortRequired>True</PortRequired><OptionsDependedOn></OptionsDependedOn><Description>Oracle Enterprise Manager (Database Control only)</Description><Name>OEM</Name><OptionGroupOptionSettings></OptionGroupOptionSettings><EngineName>oracle-ee</EngineName></OptionGroupOption>\n    \n      <OptionGroupOption><MajorEngineVersion>11.2</MajorEngineVersion><OptionsDependedOn></OptionsDependedOn><Description>Oracle Statspack</Description><Name>STATSPACK</Name><OptionGroupOptionSettings></OptionGroupOptionSettings><EngineName>oracle-ee</EngineName></OptionGroupOption>\n    \n      <OptionGroupOption><MajorEngineVersion>11.2</MajorEngineVersion><Persistent>True</Persistent><OptionsDependedOn></OptionsDependedOn><Permanent>True</Permanent><Description>Oracle Advanced Security - Transparent Data Encryption</Description><Name>TDE</Name><OptionGroupOptionSettings></OptionGroupOptionSettings><EngineName>oracle-ee</EngineName></OptionGroupOption>\n    \n      <OptionGroupOption><MajorEngineVersion>11.2</MajorEngineVersion><Persistent>True</Persistent><OptionsDependedOn></OptionsDependedOn><Permanent>True</Permanent><Description>Oracle Advanced Security - TDE with HSM</Description><Name>TDE_HSM</Name><OptionGroupOptionSettings></OptionGroupOptionSettings><EngineName>oracle-ee</EngineName></OptionGroupOption>\n    \n      <OptionGroupOption><MajorEngineVersion>11.2</MajorEngineVersion><Persistent>True</Persistent><OptionsDependedOn></OptionsDependedOn><Permanent>True</Permanent><Description>Change time zone</Description><Name>Timezone</Name><OptionGroupOptionSettings><OptionGroupOptionSetting><AllowedValues>Africa/Cairo,Africa/Casablanca,Africa/Harare,Africa/Monrovia,Africa/Nairobi,Africa/Tripoli,Africa/Windhoek,America/Araguaina,America/Asuncion,America/Bogota,America/Caracas,America/Chihuahua,America/Cuiaba,America/Denver,America/Fortaleza,America/Guatemala,America/Halifax,America/Manaus,America/Matamoros,America/Monterrey,America/Montevideo,America/Phoenix,America/Santiago,America/Tijuana,Asia/Amman,Asia/Ashgabat,Asia/Baghdad,Asia/Baku,Asia/Bangkok,Asia/Beirut,Asia/Calcutta,Asia/Damascus,Asia/Dhaka,Asia/Irkutsk,Asia/Jerusalem,Asia/Kabul,Asia/Karachi,Asia/Kathmandu,Asia/Krasnoyarsk,Asia/Magadan,Asia/Muscat,Asia/Novosibirsk,Asia/Riyadh,Asia/Seoul,Asia/Shanghai,Asia/Singapore,Asia/Taipei,Asia/Tehran,Asia/Tokyo,Asia/Ulaanbaatar,Asia/Vladivostok,Asia/Yakutsk,Asia/Yerevan,Atlantic/Azores,Australia/Adelaide,Australia/Brisbane,Australia/Darwin,Australia/Hobart,Australia/Perth,Australia/Sydney,Brazil/East,Canada/Newfoundland,Canada/Saskatchewan,Europe/Amsterdam,Europe/Athens,Europe/Dublin,Europe/Helsinki,Europe/Istanbul,Europe/Kaliningrad,Europe/Moscow,Europe/Paris,Europe/Prague,Europe/Sarajevo,Pacific/Auckland,Pacific/Fiji,Pacific/Guam,Pacific/Honolulu,Pacific/Samoa,US/Alaska,US/Central,US/Eastern,US/East-Indiana,US/Pacific,UTC</AllowedValues><ApplyType>DYNAMIC</ApplyType><DefaultValue>UTC</DefaultValue><IsModifiable>True</IsModifiable><SettingDescription>Specifies the timezone the user wants to change the system time to</SettingDescription><SettingName>TIME_ZONE</SettingName></OptionGroupOptionSetting></OptionGroupOptionSettings><EngineName>oracle-ee</EngineName></OptionGroupOption>\n    \n      <OptionGroupOption><MajorEngineVersion>11.2</MajorEngineVersion><OptionsDependedOn></OptionsDependedOn><Description>Oracle XMLDB Repository</Description><Name>XMLDB</Name><OptionGroupOptionSettings></OptionGroupOptionSettings><EngineName>oracle-ee</EngineName></OptionGroupOption>\n    \n    </OptionGroupOptions>\n  </DescribeOptionGroupOptionsResult>\n  <ResponseMetadata>\n    <RequestId>457f7bb8-9fbf-11e4-9084-5754f80d5144</RequestId>\n  </ResponseMetadata>\n</DescribeOptionGroupOptionsResponse>',
            },
            "sqlserver-ee": {
                "10.50": '<DescribeOptionGroupOptionsResponse xmlns="http://rds.amazonaws.com/doc/2014-09-01/">\n  <DescribeOptionGroupOptionsResult>\n    <OptionGroupOptions>\n    \n      <OptionGroupOption><MajorEngineVersion>10.50</MajorEngineVersion><OptionsDependedOn></OptionsDependedOn><Description>SQLServer Database Mirroring</Description><Name>Mirroring</Name><OptionGroupOptionSettings></OptionGroupOptionSettings><EngineName>sqlserver-ee</EngineName></OptionGroupOption>\n    \n      <OptionGroupOption><MajorEngineVersion>10.50</MajorEngineVersion><Persistent>True</Persistent><OptionsDependedOn></OptionsDependedOn><Description>SQL Server - Transparent Data Encryption</Description><Name>TDE</Name><OptionGroupOptionSettings></OptionGroupOptionSettings><EngineName>sqlserver-ee</EngineName></OptionGroupOption>\n    \n    </OptionGroupOptions>\n  </DescribeOptionGroupOptionsResult>\n  <ResponseMetadata>\n    <RequestId>457f7bb8-9fbf-11e4-9084-5754f80d5144</RequestId>\n  </ResponseMetadata>\n</DescribeOptionGroupOptionsResponse>',
                "11.00": '<DescribeOptionGroupOptionsResponse xmlns="http://rds.amazonaws.com/doc/2014-09-01/">\n  <DescribeOptionGroupOptionsResult>\n    <OptionGroupOptions>\n    \n      <OptionGroupOption><MajorEngineVersion>11.00</MajorEngineVersion><OptionsDependedOn></OptionsDependedOn><Description>SQLServer Database Mirroring</Description><Name>Mirroring</Name><OptionGroupOptionSettings></OptionGroupOptionSettings><EngineName>sqlserver-ee</EngineName></OptionGroupOption>\n    \n      <OptionGroupOption><MajorEngineVersion>11.00</MajorEngineVersion><Persistent>True</Persistent><OptionsDependedOn></OptionsDependedOn><Description>SQL Server - Transparent Data Encryption</Description><Name>TDE</Name><OptionGroupOptionSettings></OptionGroupOptionSettings><EngineName>sqlserver-ee</EngineName></OptionGroupOption>\n    \n    </OptionGroupOptions>\n  </DescribeOptionGroupOptionsResult>\n  <ResponseMetadata>\n    <RequestId>457f7bb8-9fbf-11e4-9084-5754f80d5144</RequestId>\n  </ResponseMetadata>\n</DescribeOptionGroupOptionsResponse>',
                "all": '<DescribeOptionGroupOptionsResponse xmlns="http://rds.amazonaws.com/doc/2014-09-01/">\n  <DescribeOptionGroupOptionsResult>\n    <OptionGroupOptions>\n    \n      <OptionGroupOption><MajorEngineVersion>10.50</MajorEngineVersion><OptionsDependedOn></OptionsDependedOn><Description>SQLServer Database Mirroring</Description><Name>Mirroring</Name><OptionGroupOptionSettings></OptionGroupOptionSettings><EngineName>sqlserver-ee</EngineName></OptionGroupOption>\n    \n      <OptionGroupOption><MajorEngineVersion>10.50</MajorEngineVersion><Persistent>True</Persistent><OptionsDependedOn></OptionsDependedOn><Description>SQL Server - Transparent Data Encryption</Description><Name>TDE</Name><OptionGroupOptionSettings></OptionGroupOptionSettings><EngineName>sqlserver-ee</EngineName></OptionGroupOption>\n    \n      <OptionGroupOption><MajorEngineVersion>11.00</MajorEngineVersion><OptionsDependedOn></OptionsDependedOn><Description>SQLServer Database Mirroring</Description><Name>Mirroring</Name><OptionGroupOptionSettings></OptionGroupOptionSettings><EngineName>sqlserver-ee</EngineName></OptionGroupOption>\n    \n      <OptionGroupOption><MajorEngineVersion>11.00</MajorEngineVersion><Persistent>True</Persistent><OptionsDependedOn></OptionsDependedOn><Description>SQL Server - Transparent Data Encryption</Description><Name>TDE</Name><OptionGroupOptionSettings></OptionGroupOptionSettings><EngineName>sqlserver-ee</EngineName></OptionGroupOption>\n    \n    </OptionGroupOptions>\n  </DescribeOptionGroupOptionsResult>\n  <ResponseMetadata>\n    <RequestId>457f7bb8-9fbf-11e4-9084-5754f80d5144</RequestId>\n  </ResponseMetadata>\n</DescribeOptionGroupOptionsResponse>',
            },
        }

        if engine_name not in default_option_group_options:
            raise RDSClientError(
                "InvalidParameterValue", f"Invalid DB engine: {engine_name}"
            )
        if (
            major_engine_version
            and major_engine_version not in default_option_group_options[engine_name]
        ):
            raise RDSClientError(
                "InvalidParameterCombination",
                f"Cannot find major version {major_engine_version} for {engine_name}",
            )
        if major_engine_version:
            return default_option_group_options[engine_name][major_engine_version]
        return default_option_group_options[engine_name]["all"]

    def modify_option_group(
        self,
        option_group_name: str,
        options_to_include: Optional[List[Dict[str, Any]]] = None,
        options_to_remove: Optional[List[Dict[str, Any]]] = None,
    ) -> "OptionGroup":
        if option_group_name not in self.option_groups:
            raise OptionGroupNotFoundFaultError(option_group_name)
        if not options_to_include and not options_to_remove:
            raise RDSClientError(
                "InvalidParameterValue",
                "At least one option must be added, modified, or removed.",
            )
        if options_to_remove:
            self.option_groups[option_group_name].remove_options(options_to_remove)
        if options_to_include:
            self.option_groups[option_group_name].add_options(options_to_include)
        return self.option_groups[option_group_name]

    def create_db_parameter_group(
        self, db_parameter_group_kwargs: Dict[str, Any]
    ) -> "DBParameterGroup":
        db_parameter_group_id = db_parameter_group_kwargs["name"]
        if db_parameter_group_id in self.db_parameter_groups:
            raise RDSClientError(
                "DBParameterGroupAlreadyExistsFault",
                f"A DB parameter group named {db_parameter_group_id} already exists.",
            )
        if not db_parameter_group_kwargs.get("description"):
            raise RDSClientError(
                "InvalidParameterValue",
                "The parameter Description must be provided and must not be blank.",
            )
        if not db_parameter_group_kwargs.get("family"):
            raise RDSClientError(
                "InvalidParameterValue",
                "The parameter DBParameterGroupName must be provided and must not be blank.",
            )
        db_parameter_group_kwargs["region"] = self.region_name
        db_parameter_group_kwargs["account_id"] = self.account_id
        db_parameter_group = DBParameterGroup(**db_parameter_group_kwargs)
        self.db_parameter_groups[db_parameter_group_id] = db_parameter_group
        return db_parameter_group

    def describe_db_parameter_groups(
        self, db_parameter_group_kwargs: Dict[str, Any]
    ) -> List["DBParameterGroup"]:
        db_parameter_group_list = []

        if db_parameter_group_kwargs.get("marker"):
            marker = db_parameter_group_kwargs["marker"]
        else:
            marker = 0
        if db_parameter_group_kwargs.get("max_records"):
            if (
                db_parameter_group_kwargs["max_records"] < 20
                or db_parameter_group_kwargs["max_records"] > 100
            ):
                raise RDSClientError(
                    "InvalidParameterValue",
                    "Invalid value for max records. Must be between 20 and 100",
                )
            max_records = db_parameter_group_kwargs["max_records"]
        else:
            max_records = 100

        for db_parameter_group in self.db_parameter_groups.values():
            if not db_parameter_group_kwargs.get(
                "name"
            ) or db_parameter_group.name == db_parameter_group_kwargs.get("name"):
                db_parameter_group_list.append(db_parameter_group)
            else:
                continue

        return db_parameter_group_list[marker : max_records + marker]

    def modify_db_parameter_group(
        self,
        db_parameter_group_name: str,
        db_parameter_group_parameters: Iterable[Dict[str, Any]],
    ) -> "DBParameterGroup":
        if db_parameter_group_name not in self.db_parameter_groups:
            raise DBParameterGroupNotFoundError(db_parameter_group_name)

        db_parameter_group = self.db_parameter_groups[db_parameter_group_name]
        db_parameter_group.update_parameters(db_parameter_group_parameters)

        return db_parameter_group

    def describe_db_cluster_parameters(self) -> List[Dict[str, Any]]:
        return []

    def create_db_cluster(self, kwargs: Dict[str, Any]) -> Cluster:
        cluster_id = kwargs["db_cluster_identifier"]
        kwargs["account_id"] = self.account_id
        cluster = Cluster(**kwargs)
        self.clusters[cluster_id] = cluster

        if (
            cluster.global_cluster_identifier
            and cluster.global_cluster_identifier in self.global_clusters
        ):
            global_cluster = self.global_clusters[cluster.global_cluster_identifier]

            # Main DB cluster, does RW on global cluster
            setattr(cluster, "is_writer", True)
            # self.clusters[cluster_id] = cluster
            global_cluster.members.append(cluster)

        # search all backend to check if global cluster named global_cluster_identifier exists
        # anywhere else
        if (
            cluster.global_cluster_identifier
            and cluster.global_cluster_identifier not in self.global_clusters
        ):
            for regional_backend in rds_backends[self.account_id]:
                if (
                    cluster.global_cluster_identifier
                    in rds_backends[self.account_id][regional_backend].global_clusters
                ):
                    global_cluster = rds_backends[self.account_id][
                        regional_backend
                    ].global_clusters[cluster.global_cluster_identifier]
                    global_cluster.members.append(cluster)

        if cluster.replication_source_identifier:
            cluster_identifier = cluster.replication_source_identifier
            original_cluster = find_cluster(cluster_identifier)
            original_cluster.read_replica_identifiers.append(cluster.db_cluster_arn)

        initial_state = copy.deepcopy(cluster)  # Return status=creating
        cluster.status = "available"  # Already set the final status in the background
        return initial_state

    def modify_db_cluster(self, kwargs: Dict[str, Any]) -> Cluster:
        cluster_id = kwargs["db_cluster_identifier"]

        if cluster_id in self.neptune.clusters:
            return self.neptune.modify_db_cluster(kwargs)  # type: ignore

        cluster = self.clusters[cluster_id]
        del self.clusters[cluster_id]

        kwargs["db_cluster_identifier"] = kwargs.pop("new_db_cluster_identifier")
        for k, v in kwargs.items():
            if v is not None:
                setattr(cluster, k, v)

        cwl_exports = kwargs.get("enable_cloudwatch_logs_exports") or {}
        for exp in cwl_exports.get("DisableLogTypes", []):
            cluster.enabled_cloudwatch_logs_exports.remove(exp)
        cluster.enabled_cloudwatch_logs_exports.extend(
            cwl_exports.get("EnableLogTypes", [])
        )

        cluster_id = kwargs.get("new_db_cluster_identifier", cluster_id)
        self.clusters[cluster_id] = cluster

        initial_state = copy.deepcopy(cluster)  # Return status=creating
        cluster.status = "available"  # Already set the final status in the background
        return initial_state

    def promote_read_replica_db_cluster(self, db_cluster_identifier: str) -> Cluster:
        cluster = self.clusters[db_cluster_identifier]
        source_cluster = find_cluster(cluster.replication_source_identifier)  # type: ignore
        source_cluster.read_replica_identifiers.remove(cluster.db_cluster_arn)
        cluster.replication_source_identifier = None
        return cluster

    def create_auto_cluster_snapshot(
        self, db_cluster_identifier: str, db_snapshot_identifier: str
    ) -> ClusterSnapshot:
        return self.create_db_cluster_snapshot(
            db_cluster_identifier, db_snapshot_identifier, snapshot_type="automated"
        )

    def create_db_cluster_snapshot(
        self,
        db_cluster_identifier: str,
        db_snapshot_identifier: str,
        snapshot_type: str = "manual",
        tags: Optional[List[Dict[str, str]]] = None,
    ) -> ClusterSnapshot:
        cluster = self.clusters.get(db_cluster_identifier)
        if cluster is None:
            raise DBClusterNotFoundError(db_cluster_identifier)
        if db_snapshot_identifier in self.cluster_snapshots:
            raise DBClusterSnapshotAlreadyExistsError(db_snapshot_identifier)
        if len(self.cluster_snapshots) >= int(
            os.environ.get("MOTO_RDS_SNAPSHOT_LIMIT", "100")
        ):
            raise SnapshotQuotaExceededError()
        if tags is None:
            tags = list()
        if cluster.copy_tags_to_snapshot:
            tags += cluster.get_tags()
        snapshot = ClusterSnapshot(cluster, db_snapshot_identifier, snapshot_type, tags)
        self.cluster_snapshots[db_snapshot_identifier] = snapshot
        return snapshot

    def copy_db_cluster_snapshot(
        self,
        source_snapshot_identifier: str,
        target_snapshot_identifier: str,
        tags: Optional[List[Dict[str, str]]] = None,
    ) -> ClusterSnapshot:
        if source_snapshot_identifier not in self.cluster_snapshots:
            raise DBClusterSnapshotNotFoundError(source_snapshot_identifier)

        source_snapshot = self.cluster_snapshots[source_snapshot_identifier]
        if tags is None:
            tags = source_snapshot.tags
        else:
            tags = self._merge_tags(source_snapshot.tags, tags)
        return self.create_db_cluster_snapshot(
            db_cluster_identifier=source_snapshot.cluster.db_cluster_identifier,  # type: ignore
            db_snapshot_identifier=target_snapshot_identifier,
            tags=tags,
        )

    def delete_db_cluster_snapshot(
        self, db_snapshot_identifier: str
    ) -> ClusterSnapshot:
        if db_snapshot_identifier not in self.cluster_snapshots:
            raise DBClusterSnapshotNotFoundError(db_snapshot_identifier)

        return self.cluster_snapshots.pop(db_snapshot_identifier)

    def describe_db_clusters(
        self, cluster_identifier: Optional[str] = None, filters: Any = None
    ) -> List[Cluster]:
        clusters = self.clusters
        clusters_neptune = self.neptune.clusters
        if cluster_identifier:
            filters = merge_filters(filters, {"db-cluster-id": [cluster_identifier]})
        if filters:
            clusters = self._filter_resources(clusters, filters, Cluster)
            clusters_neptune = self._filter_resources(
                clusters_neptune, filters, Cluster
            )
        if cluster_identifier and not (clusters or clusters_neptune):
            raise DBClusterNotFoundError(cluster_identifier)
        return list(clusters.values()) + list(clusters_neptune.values())  # type: ignore

    def describe_db_cluster_snapshots(
        self,
        db_cluster_identifier: Optional[str],
        db_snapshot_identifier: str,
        filters: Any = None,
    ) -> List[ClusterSnapshot]:
        snapshots = self.cluster_snapshots
        if db_cluster_identifier:
            filters = merge_filters(filters, {"db-cluster-id": [db_cluster_identifier]})
        if db_snapshot_identifier:
            filters = merge_filters(
                filters, {"db-cluster-snapshot-id": [db_snapshot_identifier]}
            )
        if filters:
            snapshots = self._filter_resources(snapshots, filters, ClusterSnapshot)
        if db_snapshot_identifier and not snapshots and not db_cluster_identifier:
            raise DBClusterSnapshotNotFoundError(db_snapshot_identifier)
        return list(snapshots.values())

    def delete_db_cluster(
        self, cluster_identifier: str, snapshot_name: Optional[str] = None
    ) -> Cluster:
        if cluster_identifier in self.clusters:
            cluster = self.clusters[cluster_identifier]
            if cluster.deletion_protection:
                raise InvalidParameterValue(
                    "Can't delete Cluster with protection enabled"
                )
            if cluster.cluster_members:
                raise DBClusterToBeDeletedHasActiveMembers()
            global_id = cluster.global_cluster_identifier or ""
            if global_id in self.global_clusters:
                self.remove_from_global_cluster(global_id, cluster_identifier)

            if snapshot_name:
                self.create_auto_cluster_snapshot(cluster_identifier, snapshot_name)
            return self.clusters.pop(cluster_identifier)
        if cluster_identifier in self.neptune.clusters:
            return self.neptune.delete_db_cluster(cluster_identifier)  # type: ignore
        raise DBClusterNotFoundError(cluster_identifier)

    def start_db_cluster(self, cluster_identifier: str) -> Cluster:
        if cluster_identifier not in self.clusters:
            return self.neptune.start_db_cluster(cluster_identifier)  # type: ignore
            raise DBClusterNotFoundError(cluster_identifier)
        cluster = self.clusters[cluster_identifier]
        if cluster.status != "stopped":
            raise InvalidDBClusterStateFault(
                "DbCluster cluster-id is not in stopped state."
            )
        temp_state = copy.deepcopy(cluster)
        temp_state.status = "started"
        cluster.status = "available"  # This is the final status - already setting it in the background
        return temp_state

    def restore_db_cluster_from_snapshot(
        self, from_snapshot_id: str, overrides: Dict[str, Any]
    ) -> Cluster:
        snapshot = self.describe_db_cluster_snapshots(
            db_cluster_identifier=None, db_snapshot_identifier=from_snapshot_id
        )[0]
        original_cluster = snapshot.cluster
        new_cluster_props = copy.deepcopy(original_cluster.get_cfg())
        for key, value in overrides.items():
            if value:
                new_cluster_props[key] = value

        return self.create_db_cluster(new_cluster_props)

    def stop_db_cluster(self, cluster_identifier: str) -> "Cluster":
        if cluster_identifier not in self.clusters:
            raise DBClusterNotFoundError(cluster_identifier)
        cluster = self.clusters[cluster_identifier]
        if cluster.status not in ["available"]:
            raise InvalidDBClusterStateFault(
                "DbCluster cluster-id is not in available state."
            )
        previous_state = copy.deepcopy(cluster)
        cluster.status = "stopped"
        return previous_state

    def start_export_task(self, kwargs: Dict[str, Any]) -> ExportTask:
        export_task_id = kwargs["export_task_identifier"]
        source_arn = kwargs["source_arn"]
        snapshot_id = source_arn.split(":")[-1]
        snapshot_type = source_arn.split(":")[-2]

        if export_task_id in self.export_tasks:
            raise ExportTaskAlreadyExistsError(export_task_id)
        if snapshot_type == "snapshot" and snapshot_id not in self.database_snapshots:
            raise DBSnapshotNotFoundError(snapshot_id)
        elif (
            snapshot_type == "cluster-snapshot"
            and snapshot_id not in self.cluster_snapshots
        ):
            raise DBClusterSnapshotNotFoundError(snapshot_id)

        if snapshot_type == "snapshot":
            snapshot: Union[DatabaseSnapshot, ClusterSnapshot] = (
                self.database_snapshots[snapshot_id]
            )
        else:
            snapshot = self.cluster_snapshots[snapshot_id]

        if snapshot.status not in ["available"]:
            raise InvalidExportSourceStateError(snapshot.status)

        export_task = ExportTask(snapshot, kwargs)
        self.export_tasks[export_task_id] = export_task

        return export_task

    def cancel_export_task(self, export_task_identifier: str) -> ExportTask:
        if export_task_identifier in self.export_tasks:
            export_task = self.export_tasks[export_task_identifier]
            export_task.status = "canceled"
            self.export_tasks[export_task_identifier] = export_task
            return export_task
        raise ExportTaskNotFoundError(export_task_identifier)

    def describe_export_tasks(
        self, export_task_identifier: str
    ) -> Iterable[ExportTask]:
        if export_task_identifier:
            if export_task_identifier in self.export_tasks:
                return [self.export_tasks[export_task_identifier]]
            else:
                raise ExportTaskNotFoundError(export_task_identifier)
        return self.export_tasks.values()

    def create_event_subscription(self, kwargs: Any) -> EventSubscription:
        subscription_name = kwargs["subscription_name"]

        if subscription_name in self.event_subscriptions:
            raise SubscriptionAlreadyExistError(subscription_name)

        kwargs["account_id"] = self.account_id
        subscription = EventSubscription(kwargs)
        self.event_subscriptions[subscription_name] = subscription

        return subscription

    def delete_event_subscription(self, subscription_name: str) -> EventSubscription:
        if subscription_name in self.event_subscriptions:
            return self.event_subscriptions.pop(subscription_name)
        raise SubscriptionNotFoundError(subscription_name)

    def describe_event_subscriptions(
        self, subscription_name: str
    ) -> Iterable[EventSubscription]:
        if subscription_name:
            if subscription_name in self.event_subscriptions:
                return [self.event_subscriptions[subscription_name]]
            else:
                raise SubscriptionNotFoundError(subscription_name)
        return self.event_subscriptions.values()

    def list_tags_for_resource(self, arn: str) -> List[Dict[str, str]]:
        if self.arn_regex.match(arn):
            arn_breakdown = arn.split(":")
            resource_type = arn_breakdown[len(arn_breakdown) - 2]
            resource_name = arn_breakdown[len(arn_breakdown) - 1]
            if resource_type == "db":  # Database
                if resource_name in self.databases:
                    return self.databases[resource_name].get_tags()
            elif resource_type == "cluster":  # Cluster
                if resource_name in self.clusters:
                    return self.clusters[resource_name].get_tags()
                if resource_name in self.neptune.clusters:
                    return self.neptune.clusters[resource_name].get_tags()
            elif resource_type == "es":  # Event Subscription
                if resource_name in self.event_subscriptions:
                    return self.event_subscriptions[resource_name].get_tags()
            elif resource_type == "og":  # Option Group
                if resource_name in self.option_groups:
                    return self.option_groups[resource_name].get_tags()
            elif resource_type == "pg":  # Parameter Group
                if resource_name in self.db_parameter_groups:
                    return self.db_parameter_groups[resource_name].get_tags()
            elif resource_type == "ri":  # Reserved DB instance
                # TODO: Complete call to tags on resource type Reserved DB
                # instance
                return []
            elif resource_type == "secgrp":  # DB security group
                if resource_name in self.security_groups:
                    return self.security_groups[resource_name].get_tags()
            elif resource_type == "snapshot":  # DB Snapshot
                if resource_name in self.database_snapshots:
                    return self.database_snapshots[resource_name].get_tags()
            elif resource_type == "cluster-snapshot":  # DB Cluster Snapshot
                if resource_name in self.cluster_snapshots:
                    return self.cluster_snapshots[resource_name].get_tags()
            elif resource_type == "subgrp":  # DB subnet group
                if resource_name in self.subnet_groups:
                    return self.subnet_groups[resource_name].get_tags()
            elif resource_type == "db-proxy":  # DB Proxy
                if resource_name in self.db_proxies:
                    return self.db_proxies[resource_name].get_tags()
        else:
            raise RDSClientError(
                "InvalidParameterValue", f"Invalid resource name: {arn}"
            )
        return []

    def remove_tags_from_resource(self, arn: str, tag_keys: List[str]) -> None:
        if self.arn_regex.match(arn):
            arn_breakdown = arn.split(":")
            resource_type = arn_breakdown[len(arn_breakdown) - 2]
            resource_name = arn_breakdown[len(arn_breakdown) - 1]
            if resource_type == "db":  # Database
                if resource_name in self.databases:
                    self.databases[resource_name].remove_tags(tag_keys)
            elif resource_type == "es":  # Event Subscription
                if resource_name in self.event_subscriptions:
                    self.event_subscriptions[resource_name].remove_tags(tag_keys)
            elif resource_type == "og":  # Option Group
                if resource_name in self.option_groups:
                    self.option_groups[resource_name].remove_tags(tag_keys)
            elif resource_type == "pg":  # Parameter Group
                if resource_name in self.db_parameter_groups:
                    self.db_parameter_groups[resource_name].remove_tags(tag_keys)
            elif resource_type == "ri":  # Reserved DB instance
                return None
            elif resource_type == "secgrp":  # DB security group
                if resource_name in self.security_groups:
                    self.security_groups[resource_name].remove_tags(tag_keys)
            elif resource_type == "snapshot":  # DB Snapshot
                if resource_name in self.database_snapshots:
                    self.database_snapshots[resource_name].remove_tags(tag_keys)
            elif resource_type == "cluster":
                if resource_name in self.clusters:
                    self.clusters[resource_name].remove_tags(tag_keys)
                if resource_name in self.neptune.clusters:
                    self.neptune.clusters[resource_name].remove_tags(tag_keys)
            elif resource_type == "cluster-snapshot":  # DB Cluster Snapshot
                if resource_name in self.cluster_snapshots:
                    self.cluster_snapshots[resource_name].remove_tags(tag_keys)
            elif resource_type == "subgrp":  # DB subnet group
                if resource_name in self.subnet_groups:
                    self.subnet_groups[resource_name].remove_tags(tag_keys)
            elif resource_type == "db-proxy":  # DB Proxy
                if resource_name in self.db_proxies:
                    self.db_proxies[resource_name].remove_tags(tag_keys)
        else:
            raise RDSClientError(
                "InvalidParameterValue", f"Invalid resource name: {arn}"
            )

    def add_tags_to_resource(  # type: ignore[return]
        self, arn: str, tags: List[Dict[str, str]]
    ) -> List[Dict[str, str]]:
        if self.arn_regex.match(arn):
            arn_breakdown = arn.split(":")
            resource_type = arn_breakdown[-2]
            resource_name = arn_breakdown[-1]
            if resource_type == "db":  # Database
                if resource_name in self.databases:
                    return self.databases[resource_name].add_tags(tags)
            elif resource_type == "es":  # Event Subscription
                if resource_name in self.event_subscriptions:
                    return self.event_subscriptions[resource_name].add_tags(tags)
            elif resource_type == "og":  # Option Group
                if resource_name in self.option_groups:
                    return self.option_groups[resource_name].add_tags(tags)
            elif resource_type == "pg":  # Parameter Group
                if resource_name in self.db_parameter_groups:
                    return self.db_parameter_groups[resource_name].add_tags(tags)
            elif resource_type == "ri":  # Reserved DB instance
                return []
            elif resource_type == "secgrp":  # DB security group
                if resource_name in self.security_groups:
                    return self.security_groups[resource_name].add_tags(tags)
            elif resource_type == "snapshot":  # DB Snapshot
                if resource_name in self.database_snapshots:
                    return self.database_snapshots[resource_name].add_tags(tags)
            elif resource_type == "cluster":
                if resource_name in self.clusters:
                    return self.clusters[resource_name].add_tags(tags)
                if resource_name in self.neptune.clusters:
                    return self.neptune.clusters[resource_name].add_tags(tags)
            elif resource_type == "cluster-snapshot":  # DB Cluster Snapshot
                if resource_name in self.cluster_snapshots:
                    return self.cluster_snapshots[resource_name].add_tags(tags)
            elif resource_type == "subgrp":  # DB subnet group
                if resource_name in self.subnet_groups:
                    return self.subnet_groups[resource_name].add_tags(tags)
            elif resource_type == "db-proxy":  # DB Proxy
                if resource_name in self.db_proxies:
                    return self.db_proxies[resource_name].add_tags(tags)
        else:
            raise RDSClientError(
                "InvalidParameterValue", f"Invalid resource name: {arn}"
            )

    @staticmethod
    def _filter_resources(resources: Any, filters: Any, resource_class: Any) -> Any:  # type: ignore[misc]
        try:
            filter_defs = resource_class.SUPPORTED_FILTERS
            validate_filters(filters, filter_defs)
            return apply_filter(resources, filters, filter_defs)
        except KeyError as e:
            # https://stackoverflow.com/questions/24998968/why-does-strkeyerror-add-extra-quotes
            raise InvalidParameterValue(e.args[0])
        except ValueError as e:
            raise InvalidParameterCombination(str(e))

    @staticmethod
    def _merge_tags(  # type: ignore[misc]
        old_tags: List[Dict[str, Any]], new_tags: List[Dict[str, Any]]
    ) -> List[Dict[str, Any]]:
        tags_dict = dict()
        tags_dict.update({d["Key"]: d["Value"] for d in old_tags})
        tags_dict.update({d["Key"]: d["Value"] for d in new_tags})
        return [{"Key": k, "Value": v} for k, v in tags_dict.items()]

    @staticmethod
    def _validate_db_identifier(db_identifier: str) -> None:
        # https://docs.aws.amazon.com/AmazonRDS/latest/APIReference/API_CreateDBInstance.html
        # Constraints:
        # # Must contain from 1 to 63 letters, numbers, or hyphens.
        # # First character must be a letter.
        # # Can't end with a hyphen or contain two consecutive hyphens.
        if re.match(
            "^(?!.*--)([a-zA-Z]?[a-zA-Z0-9-]{0,61}[a-zA-Z0-9])$", db_identifier
        ):
            return
        raise InvalidDBInstanceIdentifier

    def describe_orderable_db_instance_options(
        self, engine: str, engine_version: str
    ) -> List[Dict[str, Any]]:
        """
        Only the Aurora-Postgresql and Neptune-engine is currently implemented
        """
        if engine == "neptune":
            return self.neptune.describe_orderable_db_instance_options(engine_version)
        if engine == "aurora-postgresql":
            if engine_version:
                return [
                    option
                    for option in self.db_cluster_options
                    if option["EngineVersion"] == engine_version
                ]
            return self.db_cluster_options
        return []

    def create_db_cluster_parameter_group(
        self,
        group_name: str,
        family: str,
        description: str,
    ) -> "DBClusterParameterGroup":
        group = DBClusterParameterGroup(
            account_id=self.account_id,
            region=self.region_name,
            name=group_name,
            family=family,
            description=description,
        )
        self.db_cluster_parameter_groups[group_name] = group
        return group

    def describe_db_cluster_parameter_groups(
        self, group_name: str
    ) -> List["DBClusterParameterGroup"]:
        if group_name is not None:
            if group_name not in self.db_cluster_parameter_groups:
                raise DBClusterParameterGroupNotFoundError(group_name)
            return [self.db_cluster_parameter_groups[group_name]]
        return list(self.db_cluster_parameter_groups.values())

    def delete_db_cluster_parameter_group(self, group_name: str) -> None:
        self.db_cluster_parameter_groups.pop(group_name)

    def create_global_cluster(
        self,
        global_cluster_identifier: str,
        source_db_cluster_identifier: Optional[str],
        engine: Optional[str],
        engine_version: Optional[str],
        storage_encrypted: Optional[str],
        deletion_protection: Optional[str],
    ) -> GlobalCluster:
        source_cluster = None
        if source_db_cluster_identifier is not None:
            # validate our source cluster exists
            if not re.match(ARN_PARTITION_REGEX + ":rds", source_db_cluster_identifier):
                raise InvalidParameterValue("Malformed db cluster arn dbci")
            source_cluster = self.describe_db_clusters(
                cluster_identifier=source_db_cluster_identifier
            )[0]
            # We should not specify an engine at the same time, as we'll take it from the source cluster
            if engine is not None:
                raise InvalidParameterCombination(
                    "When creating global cluster from existing db cluster, value for engineName should not be specified since it will be inherited from source cluster"
                )
            engine = source_cluster.engine
            engine_version = source_cluster.engine_version
        elif engine is None:
            raise InvalidParameterValue(
                "When creating standalone global cluster, value for engineName should be specified"
            )
        global_cluster = GlobalCluster(
            account_id=self.account_id,
            region_name=self.region_name,
            global_cluster_identifier=global_cluster_identifier,
            engine=engine,  # type: ignore
            engine_version=engine_version,
            storage_encrypted=storage_encrypted,
            deletion_protection=deletion_protection,
        )
        self.global_clusters[global_cluster_identifier] = global_cluster
        if source_cluster is not None:
            source_cluster.global_cluster_identifier = global_cluster.global_cluster_arn
            global_cluster.members.append(source_cluster)
        return global_cluster

    def describe_global_clusters(self) -> List[GlobalCluster]:
        return (
            list(self.global_clusters.values())
            + self.neptune.describe_global_clusters()  # type: ignore
        )

    def delete_global_cluster(self, global_cluster_identifier: str) -> GlobalCluster:
        try:
            return self.neptune.delete_global_cluster(global_cluster_identifier)  # type: ignore
        except:  # noqa: E722 Do not use bare except
            pass  # It's not a Neptune Global Cluster - assume it's an RDS cluster instead
        global_cluster = self.global_clusters[global_cluster_identifier]
        if global_cluster.members:
            raise InvalidGlobalClusterStateFault(global_cluster.global_cluster_arn)
        return self.global_clusters.pop(global_cluster_identifier)

    def remove_from_global_cluster(
        self, global_cluster_identifier: str, db_cluster_identifier: str
    ) -> Optional[GlobalCluster]:
        try:
            global_cluster = self.global_clusters[global_cluster_identifier]
            cluster = self.describe_db_clusters(
                cluster_identifier=db_cluster_identifier
            )[0]
            global_cluster.members.remove(cluster)
            return global_cluster
        except:  # noqa: E722 Do not use bare except
            pass
        return None

    def describe_db_snapshot_attributes(
        self, db_snapshot_identifier: str
    ) -> List[Dict[str, Any]]:
        snapshot = self.describe_db_snapshots(
            db_instance_identifier=None, db_snapshot_identifier=db_snapshot_identifier
        )[0]
        return snapshot.attributes

    def modify_db_snapshot_attribute(
        self,
        db_snapshot_identifier: str,
        attribute_name: str,
        values_to_add: Optional[Dict[str, Dict[str, str]]] = None,
        values_to_remove: Optional[Dict[str, Dict[str, str]]] = None,
    ) -> List[Dict[str, Any]]:
        snapshot = self.describe_db_snapshots(
            db_instance_identifier=None, db_snapshot_identifier=db_snapshot_identifier
        )[0]
        attribute_present = False
        for attribute in snapshot.attributes:
            if attribute["AttributeName"] == attribute_name:
                attribute_present = True
                if values_to_add:
                    attribute["AttributeValues"] = (
                        values_to_add["AttributeValue"].values()
                        + attribute["AttributeValues"]
                    )
                if values_to_remove:
                    attribute["AttributeValues"] = [
                        i
                        for i in attribute["AttributeValues"]
                        if i not in values_to_remove["AttributeValue"].values()
                    ]
        if not attribute_present and values_to_add:
            snapshot.attributes.append(
                {
                    "AttributeName": attribute_name,
                    "AttributeValues": values_to_add["AttributeValue"].values(),
                }
            )
        return snapshot.attributes

    def describe_db_cluster_snapshot_attributes(
        self, db_cluster_snapshot_identifier: str
    ) -> List[Dict[str, Any]]:
        snapshot = self.describe_db_cluster_snapshots(
            db_cluster_identifier=None,
            db_snapshot_identifier=db_cluster_snapshot_identifier,
        )[0]
        return snapshot.attributes

    def modify_db_cluster_snapshot_attribute(
        self,
        db_cluster_snapshot_identifier: str,
        attribute_name: str,
        values_to_add: Optional[Dict[str, Dict[str, str]]] = None,
        values_to_remove: Optional[Dict[str, Dict[str, str]]] = None,
    ) -> List[Dict[str, Any]]:
        snapshot = self.describe_db_cluster_snapshots(
            db_cluster_identifier=None,
            db_snapshot_identifier=db_cluster_snapshot_identifier,
        )[0]
        attribute_present = False
        for attribute in snapshot.attributes:
            if attribute["AttributeName"] == attribute_name:
                attribute_present = True
                if values_to_add:
                    attribute["AttributeValues"] = (
                        values_to_add["AttributeValue"].values()
                        + attribute["AttributeValues"]
                    )
                if values_to_remove:
                    attribute["AttributeValues"] = [
                        i
                        for i in attribute["AttributeValues"]
                        if i not in values_to_remove["AttributeValue"].values()
                    ]
        if not attribute_present and values_to_add:
            snapshot.attributes.append(
                {
                    "AttributeName": attribute_name,
                    "AttributeValues": values_to_add["AttributeValue"].values(),
                }
            )
        return snapshot.attributes

    def create_db_proxy(
        self,
        db_proxy_name: str,
        engine_family: str,
        auth: List[Dict[str, str]],
        role_arn: str,
        vpc_subnet_ids: List[str],
        vpc_security_group_ids: Optional[List[str]],
        require_tls: Optional[bool],
        idle_client_timeout: Optional[int],
        debug_logging: Optional[bool],
        tags: Optional[List[Dict[str, str]]],
    ) -> DBProxy:
        self._validate_db_identifier(db_proxy_name)
        if db_proxy_name in self.db_proxies:
            raise DBProxyAlreadyExistsFault(db_proxy_name)
        if len(self.db_proxies) >= int(os.environ.get("MOTO_RDS_PROXY_LIMIT", "100")):
            raise DBProxyQuotaExceededFault()
        db_proxy = DBProxy(
            db_proxy_name,
            engine_family,
            auth,
            role_arn,
            vpc_subnet_ids,
            self.region_name,
            self.account_id,
            vpc_security_group_ids,
            require_tls,
            idle_client_timeout,
            debug_logging,
            tags,
        )
        self.db_proxies[db_proxy_name] = db_proxy
        return db_proxy

    def describe_db_proxies(
        self,
        db_proxy_name: Optional[str],
        filters: Optional[List[Dict[str, Any]]] = None,
    ) -> List[DBProxy]:
        """
        The filters-argument is not yet supported
        """
        db_proxies = list(self.db_proxies.values())
        if db_proxy_name and db_proxy_name in self.db_proxies.keys():
            db_proxies = [self.db_proxies[db_proxy_name]]
        if db_proxy_name and db_proxy_name not in self.db_proxies.keys():
            raise DBProxyNotFoundFault(db_proxy_name)
        return db_proxies


class OptionGroup:
    def __init__(
        self,
        name: str,
        engine_name: str,
        major_engine_version: str,
        region: str,
        account_id: str,
        description: Optional[str] = None,
    ):
        self.engine_name = engine_name
        self.major_engine_version = major_engine_version
        self.description = description
        self.name = name
        self.vpc_and_non_vpc_instance_memberships = False
        self.options: Dict[str, Any] = {}
        self.vpcId = "null"
        self.tags: List[Dict[str, str]] = []
        self.arn = f"arn:{get_partition(region)}:rds:{region}:{account_id}:og:{name}"

    def to_json(self) -> str:
        template = Template(
            """{
    "VpcId": null,
    "MajorEngineVersion": "{{ option_group.major_engine_version }}",
    "OptionGroupDescription": "{{ option_group.description }}",
    "AllowsVpcAndNonVpcInstanceMemberships": "{{ option_group.vpc_and_non_vpc_instance_memberships }}",
    "EngineName": "{{ option_group.engine_name }}",
    "Options": [],
    "OptionGroupName": "{{ option_group.name }},
    "OptionGroupArn": "{{ option_group.arn }}
}"""
        )
        return template.render(option_group=self)

    def to_xml(self) -> str:
        template = Template(
            """<OptionGroup>
          <OptionGroupName>{{ option_group.name }}</OptionGroupName>
          <AllowsVpcAndNonVpcInstanceMemberships>{{ option_group.vpc_and_non_vpc_instance_memberships }}</AllowsVpcAndNonVpcInstanceMemberships>
          <MajorEngineVersion>{{ option_group.major_engine_version }}</MajorEngineVersion>
          <EngineName>{{ option_group.engine_name }}</EngineName>
          <OptionGroupDescription>{{ option_group.description }}</OptionGroupDescription>
          <OptionGroupArn>{{ option_group.arn }}</OptionGroupArn>
          <Options/>
        </OptionGroup>"""
        )
        return template.render(option_group=self)

    def remove_options(
        self,
        options_to_remove: Any,  # pylint: disable=unused-argument
    ) -> None:
        # TODO: Check for option in self.options and remove if exists. Raise
        # error otherwise
        return

    def add_options(
        self,
        options_to_add: Any,  # pylint: disable=unused-argument
    ) -> None:
        # TODO: Validate option and add it to self.options. If invalid raise
        # error
        return

    def get_tags(self) -> List[Dict[str, str]]:
        return self.tags

    def add_tags(self, tags: List[Dict[str, str]]) -> List[Dict[str, str]]:
        new_keys = [tag_set["Key"] for tag_set in tags]
        self.tags = [tag_set for tag_set in self.tags if tag_set["Key"] not in new_keys]
        self.tags.extend(tags)
        return self.tags

    def remove_tags(self, tag_keys: List[str]) -> None:
        self.tags = [tag_set for tag_set in self.tags if tag_set["Key"] not in tag_keys]


class DBParameterGroup(CloudFormationModel):
    def __init__(
        self,
        account_id: str,
        name: Optional[str],
        description: str,
        family: Optional[str],
        tags: List[Dict[str, str]],
        region: str,
    ):
        self.name = name
        self.description = description
        self.family = family
        self.tags = tags
        self.parameters: Dict[str, Any] = defaultdict(dict)
        self.arn = f"arn:{get_partition(region)}:rds:{region}:{account_id}:pg:{name}"

    def to_xml(self) -> str:
        template = Template(
            """<DBParameterGroup>
          <DBParameterGroupName>{{ param_group.name }}</DBParameterGroupName>
          <DBParameterGroupFamily>{{ param_group.family }}</DBParameterGroupFamily>
          <Description>{{ param_group.description }}</Description>
          <DBParameterGroupArn>{{ param_group.arn }}</DBParameterGroupArn>
        </DBParameterGroup>"""
        )
        return template.render(param_group=self)

    def get_tags(self) -> List[Dict[str, Any]]:
        return self.tags

    def add_tags(self, tags: List[Dict[str, str]]) -> List[Dict[str, Any]]:
        new_keys = [tag_set["Key"] for tag_set in tags]
        self.tags = [tag_set for tag_set in self.tags if tag_set["Key"] not in new_keys]
        self.tags.extend(tags)
        return self.tags

    def remove_tags(self, tag_keys: List[str]) -> None:
        self.tags = [tag_set for tag_set in self.tags if tag_set["Key"] not in tag_keys]

    def update_parameters(self, new_parameters: Iterable[Dict[str, Any]]) -> None:
        for new_parameter in new_parameters:
            parameter = self.parameters[new_parameter["ParameterName"]]
            parameter.update(new_parameter)

    def delete(self, account_id: str, region_name: str) -> None:
        backend = rds_backends[account_id][region_name]
        backend.delete_db_parameter_group(self.name)  # type: ignore[arg-type]

    @staticmethod
    def cloudformation_name_type() -> str:
        return ""

    @staticmethod
    def cloudformation_type() -> str:
        # https://docs.aws.amazon.com/AWSCloudFormation/latest/UserGuide/aws-resource-rds-dbparametergroup.html
        return "AWS::RDS::DBParameterGroup"

    @classmethod
    def create_from_cloudformation_json(  # type: ignore[misc]
        cls,
        resource_name: str,
        cloudformation_json: Any,
        account_id: str,
        region_name: str,
        **kwargs: Any,
    ) -> "DBParameterGroup":
        properties = cloudformation_json["Properties"]

        db_parameter_group_kwargs = {
            "description": properties["Description"],
            "family": properties["Family"],
            "name": resource_name.lower(),
            "tags": properties.get("Tags"),
        }
        db_parameter_group_parameters = []
        for db_parameter, db_parameter_value in properties.get(
            "Parameters", {}
        ).items():
            db_parameter_group_parameters.append(
                {"ParameterName": db_parameter, "ParameterValue": db_parameter_value}
            )

        rds_backend = rds_backends[account_id][region_name]
        db_parameter_group = rds_backend.create_db_parameter_group(
            db_parameter_group_kwargs
        )
        db_parameter_group.update_parameters(db_parameter_group_parameters)
        return db_parameter_group


class DBClusterParameterGroup(CloudFormationModel):
    def __init__(
        self, account_id: str, region: str, name: str, description: str, family: str
    ):
        self.name = name
        self.description = description
        self.family = family
        self.parameters: Dict[str, Any] = defaultdict(dict)
        self.arn = f"arn:{get_partition(region)}:rds:{region}:{account_id}:cpg:{name}"

    def to_xml(self) -> str:
        template = Template(
            """<DBClusterParameterGroup>
          <DBClusterParameterGroupName>{{ param_group.name }}</DBClusterParameterGroupName>
          <DBParameterGroupFamily>{{ param_group.family }}</DBParameterGroupFamily>
          <Description>{{ param_group.description }}</Description>
          <DBClusterParameterGroupArn>{{ param_group.arn }}</DBClusterParameterGroupArn>
        </DBClusterParameterGroup>"""
        )
        return template.render(param_group=self)


rds_backends = BackendDict(RDSBackend, "rds")
