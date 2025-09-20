from __future__ import annotations

import copy
import json
import math
import os
import re
import string
import uuid
from collections import OrderedDict, defaultdict
from dataclasses import dataclass
from datetime import datetime
from enum import Enum
from functools import lru_cache, partialmethod
from re import compile as re_compile
from typing import (
    TYPE_CHECKING,
    Any,
    Dict,
    Iterable,
    List,
    Literal,
    MutableMapping,
    Optional,
    Protocol,
    Tuple,
    Union,
    overload,
)

from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import BaseModel, CloudFormationModel
from moto.core.utils import unix_time, utcnow
from moto.ec2.models import ec2_backends
from moto.kms.models import KmsBackend, kms_backends
from moto.moto_api._internal import mock_random as random
from moto.secretsmanager.models import FakeSecret, SecretsManagerBackend
from moto.utilities.utils import ARN_PARTITION_REGEX, CaseInsensitiveDict, load_resource

from .exceptions import (
    BlueGreenDeploymentAlreadyExistsFault,
    BlueGreenDeploymentNotFoundFault,
    DBClusterNotFoundError,
    DBClusterParameterGroupNotFoundError,
    DBClusterSnapshotAlreadyExistsError,
    DBClusterSnapshotNotFoundError,
    DBClusterToBeDeletedHasActiveMembers,
    DBInstanceAlreadyExists,
    DBInstanceNotFoundError,
    DBParameterGroupAlreadyExistsError,
    DBParameterGroupNotFoundError,
    DBProxyAlreadyExistsFault,
    DBProxyNotFoundFault,
    DBProxyQuotaExceededFault,
    DBSecurityGroupNotFoundError,
    DBShardGroupAlreadyExistsError,
    DBShardGroupNotFoundFault,
    DBSnapshotAlreadyExistsError,
    DBSnapshotNotFoundFault,
    DBSubnetGroupNotFoundError,
    ExportTaskAlreadyExistsError,
    ExportTaskNotFoundError,
    InvalidBlueGreenDeploymentStateFault,
    InvalidDBClusterSnapshotStateFault,
    InvalidDBClusterStateFault,
    InvalidDBClusterStateFaultError,
    InvalidDBInstanceEngine,
    InvalidDBInstanceIdentifier,
    InvalidDBInstanceStateError,
    InvalidDBInstanceStateFault,
    InvalidDBSnapshotIdentifier,
    InvalidExportSourceStateError,
    InvalidGlobalClusterStateFault,
    InvalidParameterCombination,
    InvalidParameterValue,
    InvalidSubnet,
    KMSKeyNotAccessibleFault,
    OptionGroupNotFoundFaultError,
    RDSClientError,  # TODO: Refactor into specific exceptions
    SharedSnapshotQuotaExceeded,
    SnapshotQuotaExceededFault,
    SourceClusterNotSupportedFault,
    SourceDatabaseNotSupportedFault,
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


def find_cluster(cluster_arn: str) -> DBCluster:
    arn_parts = cluster_arn.split(":")
    region, account = arn_parts[3], arn_parts[4]
    return rds_backends[account][region].describe_db_clusters(cluster_arn)[0]


KMS_ARN_PATTERN = re.compile(
    r"^arn:(aws|aws-us-gov|aws-cn):kms:(?P<region>\w+(?:-\w+)+):(?P<account_id>\d{12}):key\/(?P<key>[A-Za-z0-9]+(?:-[A-Za-z0-9]+)+)$"
)


class SnapshotAttributesMixin:
    ALLOWED_ATTRIBUTE_NAMES = ["restore"]

    attributes: Dict[str, List[str]]

    def __init__(self, **kwargs: Any) -> None:
        super().__init__(**kwargs)
        self.attributes = defaultdict(list)
        for attribute in self.ALLOWED_ATTRIBUTE_NAMES:
            self.attributes[attribute] = []

    def modify_attribute(
        self,
        attribute_name: str,
        values_to_add: Optional[List[str]],
        values_to_remove: Optional[List[str]],
    ) -> None:
        if not values_to_add:
            values_to_add = []
        if not values_to_remove:
            values_to_remove = []
        if attribute_name not in self.ALLOWED_ATTRIBUTE_NAMES:
            raise InvalidParameterValue(f"Invalid snapshot attribute {attribute_name}")
        common_values = set(values_to_add).intersection(values_to_remove)
        if common_values:
            raise InvalidParameterCombination(
                "A value may not appear in both the add list and remove list. "
                + f"{list(common_values)}"
            )
        add = self.attributes[attribute_name] + values_to_add
        new_attribute_values = [value for value in add if value not in values_to_remove]
        if len(new_attribute_values) > int(os.getenv("MAX_SHARED_ACCOUNTS", 20)):
            raise SharedSnapshotQuotaExceeded()
        self.attributes[attribute_name] = new_attribute_values


class ResourceWithEvents(Protocol):
    arn: str
    event_source_type: str

    @property
    def resource_id(self) -> str:
        raise NotImplementedError()


class EventMixin:
    arn: str
    backend: RDSBackend
    event_source_type: str

    def add_event(self, event_type: str) -> None:
        self.backend.add_event(event_type, self)

    @property
    def resource_id(self) -> str:
        raise NotImplementedError(
            "Concrete classes must implement resource_id property."
        )


class TaggingMixin:
    _tags: List[Dict[str, str]] = []

    @property
    def tags(self) -> List[Dict[str, str]]:
        return self._tags

    @tags.setter
    def tags(self, value: Optional[List[Dict[str, str]]]) -> None:
        if value is None:
            value = []
        # Tags may come in as XFormedDict and we want a regular dict.
        coerced = [{"Key": tag["Key"], "Value": tag["Value"]} for tag in value]
        self._tags = coerced

    @property
    def tag_list(self) -> List[Dict[str, str]]:
        return self.tags

    def get_tags(self) -> List[Dict[str, str]]:
        return self.tags

    def add_tags(self, tags: List[Dict[str, str]]) -> List[Dict[str, str]]:
        new_keys = [tag_set["Key"] for tag_set in tags]
        updated_tags = [
            tag_set for tag_set in self.tags if tag_set["Key"] not in new_keys
        ]
        updated_tags.extend(tags)
        self.tags = updated_tags
        return self.tags

    def remove_tags(self, tag_keys: List[str]) -> None:
        self.tags = [tag_set for tag_set in self.tags if tag_set["Key"] not in tag_keys]


class RDSBaseModel(TaggingMixin, BaseModel):
    resource_type: str

    def __init__(self, backend: RDSBackend, **kwargs: Any) -> None:
        super().__init__(**kwargs)
        self.backend = backend
        self.created = utcnow()

    @property
    def resource_id(self) -> str:
        raise NotImplementedError("Subclasses must implement resource_id property.")

    @property
    def region(self) -> str:
        return self.backend.region_name

    @property
    def account_id(self) -> str:
        return self.backend.account_id

    @property
    def partition(self) -> str:
        return self.backend.partition

    @property
    def arn(self) -> str:
        return f"arn:{self.partition}:rds:{self.region}:{self.account_id}:{self.resource_type}:{self.resource_id}"


class DBProxyTarget(RDSBaseModel):
    resource_type = "proxy-target"

    def __init__(
        self,
        backend: RDSBackend,
        resource_id: str,
        endpoint: Optional[str],
        type: str,
    ):
        super().__init__(backend)
        self.endpoint = endpoint
        self.rds_resource_id = resource_id
        self.type = type
        self.port = 5432
        self._registering = True
        # Not implemented yet:
        self.role = None
        self.target_arn = None

    @property
    def registering(self) -> bool:
        if self._registering is True:
            self._registering = False
            return True
        return self._registering

    @property
    def target_health(self) -> Dict[str, str]:
        return {
            "State": "REGISTERING" if self.registering else "AVAILABLE",
        }


class DBProxyTargetGroup(RDSBaseModel):
    resource_type = "target-group"

    def __init__(
        self,
        backend: RDSBackend,
        name: str,
        proxy_name: str,
    ):
        super().__init__(backend)
        self._name = f"prx-tg-{random.get_random_string(length=17, lower_case=True)}"
        self.target_group_name = name
        self.db_proxy_name = proxy_name
        self.targets: List[DBProxyTarget] = []

        self.max_connections = 100
        self.max_idle_connections = 50
        self.borrow_timeout = 120
        self.session_pinning_filters: List[str] = []

        self.created_date = utcnow()
        self.updated_date = utcnow()

        self.status = "available"
        self.is_default = True

    @property
    def resource_id(self) -> str:
        return self._name

    @property
    def target_group_arn(self) -> str:
        return self.arn

    @property
    def connection_pool_config(self) -> Dict[str, Any]:  # type: ignore[misc]
        return {
            "MaxConnectionsPercent": self.max_connections,
            "MaxIdleConnectionsPercent": self.max_idle_connections,
            "ConnectionBorrowTimeout": self.borrow_timeout,
            "SessionPinningFilters": [
                filter_ for filter_ in self.session_pinning_filters
            ],
        }


class GlobalCluster(RDSBaseModel):
    resource_type = "global-cluster"

    def __init__(
        self,
        backend: RDSBackend,
        global_cluster_identifier: str,
        engine: str,
        engine_version: Optional[str],
        storage_encrypted: Optional[bool],
        deletion_protection: Optional[bool],
    ):
        super().__init__(backend)
        self.global_cluster_identifier = global_cluster_identifier
        self.unique_identifier = random.get_random_hex(8)
        self.global_cluster_resource_id = f"cluster-{self.unique_identifier}"
        self.engine = engine
        self.engine_version = engine_version or DBCluster.default_engine_version(
            self.engine
        )
        self.storage_encrypted = storage_encrypted
        if self.storage_encrypted is None:
            self.storage_encrypted = False
        self.deletion_protection = deletion_protection
        if self.deletion_protection is None:
            self.deletion_protection = False
        self.members: List[DBCluster] = []
        self.status = "available"

    @property
    def global_cluster_identifier(self) -> str:
        return self._global_cluster_identifier

    @global_cluster_identifier.setter
    def global_cluster_identifier(self, value: str) -> None:
        self._global_cluster_identifier = value.lower()

    @property
    def resource_id(self) -> str:
        return self.global_cluster_identifier

    @property
    def arn(self) -> str:
        # Global Clusters do not belong to a particular region.
        return super().arn.replace(self.region, "")

    @property
    def endpoint(self) -> str:
        ep = f"{self.global_cluster_identifier}.global-{self.unique_identifier}.global.rds.amazonaws.com"
        return ep

    @property
    def global_cluster_arn(self) -> str:
        return self.arn

    @property
    def readers(self) -> List[str]:
        readers = [
            reader.db_cluster_arn for reader in self.members if not reader.is_writer
        ]
        return readers

    @property
    def global_cluster_members(self) -> List[Dict[str, Any]]:  # type: ignore[misc]
        members: List[Dict[str, Any]] = [
            {
                "DBClusterArn": member.db_cluster_arn,
                "IsWriter": True if member.is_writer else False,
                "DBClusterParameterGroupStatus": "in-sync",
                "PromotionTier": 1,
                # Not sure if this is correct, but current tests assert it being empty for non writers.
                "Readers": [],
            }
            for member in self.members
        ]
        for member in members:
            if member["IsWriter"]:
                member["Readers"] = self.readers
            else:
                member["GlobalWriteForwardingStatus"] = "disabled"
        return members


class MasterUserSecret:
    """Class to manage the master user password for RDS clusters/instances."""

    def __init__(self, resource: DBCluster | DBInstance, kms_key_id: str | None):
        self.resource = resource
        self.kms = resource.backend.kms
        self.secretsmanager = resource.backend.secretsmanager
        self.secret = self._create_secret(kms_key_id)
        self._status = "creating"

    def delete_secret(self) -> None:
        self.secretsmanager.delete_secret(
            self.secret.arn,
            recovery_window_in_days=None,
            force_delete_without_recovery=True,
        )

    def rotate_secret(self) -> None:
        self._status = "rotating"
        self.secretsmanager.rotate_secret(self.secret.arn, rotate_immediately=True)

    @property
    def kms_key_id(self) -> str:
        assert self.secret.kms_key_id is not None
        key = self.kms.describe_key(self.secret.kms_key_id)
        return key.arn

    @property
    def secret_arn(self) -> str:
        return self.secret.arn

    @property
    def secret_status(self) -> str:
        status_to_return = self._status
        # Poor man's state manager - look into using moto's state manager.
        if status_to_return in ["creating", "rotating"]:
            self._status = "active"
        return status_to_return

    def _create_secret(self, kms_key_id: Optional[str]) -> FakeSecret:
        secret = self.secretsmanager.create_managed_secret(
            service_name="rds",
            secret_id=self._generate_secret_name(),
            secret_string=self._generate_secret_string(),
            description=self._generate_secret_description(),
            kms_key_id=kms_key_id,
            tags=self._generate_secret_tags(),
        )
        return secret

    def _generate_secret_name(self) -> str:
        unique_id = str(uuid.uuid4())
        secret_name = f"rds!{self.resource.resource_type}-{unique_id}"
        return secret_name

    def _generate_secret_description(self) -> str:
        resource_type = self.resource.__class__.__name__[2:].lower()
        description = f"The secret associated with the primary RDS DB {resource_type}: {self.resource.arn}"
        return description

    def _generate_secret_string(self) -> str:
        credentials = {"username": "admin", "password": "P@55w0rd!"}
        secret_string = json.dumps(credentials)
        return secret_string

    def _generate_secret_tags(self) -> list[dict[str, str]]:
        resource_type = self.resource.__class__.__name__
        tags = [
            {"Key": f"aws:rds:primary{resource_type}Arn", "Value": self.resource.arn},
        ]
        return tags


class DomainMembership:
    def __init__(
        self,
        domain: Optional[str] = None,
        iam_role_name: Optional[str] = None,
        domain_ou: Optional[str] = None,
        domain_fqdn: Optional[str] = None,
        auth_secret_arn: Optional[str] = None,
        dns_ips: Optional[List[str]] = None,
    ):
        self.domain = domain.split(".")[0] if domain else None
        self.status = "active"
        self.fqdn = domain_fqdn or domain
        self.iam_role_name = iam_role_name or "rds-directory-service-access-role"
        self.ou = domain_ou or None
        if self.ou is None and self.domain:
            self.ou = f"OU={self.domain}OU,DC={self.domain},DC=com"
        self.auth_secret_arn = auth_secret_arn
        self.dns_ips = dns_ips or []


class DBCluster(RDSBaseModel):
    SUPPORTED_FILTERS = {
        "db-cluster-id": FilterDef(
            ["db_cluster_arn", "db_cluster_identifier"],
            "DB Cluster Identifiers",
            case_insensitive=True,
        ),
        "engine": FilterDef(["engine"], "Engine Names"),
    }

    resource_type = "cluster"

    def __init__(
        self,
        backend: RDSBackend,
        db_cluster_identifier: str,
        engine: str,
        allocated_storage: Optional[int] = None,
        auto_minor_version_upgrade: bool = True,
        engine_version: Optional[str] = None,
        master_username: Optional[str] = None,
        master_user_password: Optional[str] = None,
        backup_retention_period: int = 1,
        domain: Optional[str] = None,
        domain_iam_role_name: Optional[str] = None,
        domain_ou: Optional[str] = None,
        domain_fqdn: Optional[str] = None,
        character_set_name: Optional[str] = None,
        copy_tags_to_snapshot: Optional[bool] = False,
        database_name: Optional[str] = None,
        db_cluster_parameter_group_name: Optional[str] = None,
        db_subnet_group_name: Optional[str] = None,
        license_model: str = "general-public-license",
        port: Optional[int] = None,
        preferred_backup_window: str = "01:37-02:07",
        preferred_maintenance_window: str = "wed:02:40-wed:03:10",
        publicly_accessible: bool = False,
        storage_encrypted: bool = False,
        tags: Optional[List[Dict[str, str]]] = None,
        vpc_security_group_ids: Optional[List[str]] = None,
        deletion_protection: Optional[bool] = False,
        kms_key_id: Optional[str] = None,
        manage_master_user_password: Optional[bool] = False,
        master_user_secret_kms_key_id: Optional[str] = None,
        **kwargs: Any,
    ):
        super().__init__(backend)
        self.database_name = database_name
        self.url_identifier = "".join(
            random.choice(string.ascii_lowercase + string.digits) for _ in range(12)
        )
        self.db_cluster_identifier = db_cluster_identifier
        self.db_cluster_instance_class = kwargs.get("db_cluster_instance_class")
        self.auto_minor_version_upgrade = auto_minor_version_upgrade
        self.deletion_protection = deletion_protection
        self.engine = engine
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
        self.engine_version = engine_version or DBCluster.default_engine_version(
            self.engine
        )
        semantic = self.engine_version.split(".")
        option_suffix = semantic[0]
        if len(semantic) > 1:
            option_suffix = option_suffix + "-" + semantic[1]
        self.option_group_name = f"default:{self.engine}-{option_suffix}"
        self.engine_mode = kwargs.get("engine_mode") or "provisioned"
        self.iops = kwargs.get("iops")
        self.network_type = kwargs.get("network_type") or "IPV4"
        self._status = "creating"
        self.cluster_create_time = self.created
        self.copy_tags_to_snapshot = copy_tags_to_snapshot
        self.storage_type = kwargs.get("storage_type")
        if self.storage_type is None:
            self.storage_type = DBCluster.default_storage_type(iops=self.iops)
        self.allocated_storage = (
            allocated_storage
            or DBCluster.default_allocated_storage(
                engine=self.engine, storage_type=self.storage_type
            )
        )
        self.master_username = master_username
        self.character_set_name = character_set_name
        self.domain_memberships: List[DomainMembership] = []
        if domain or domain_iam_role_name or domain_ou or domain_fqdn:
            domain_membership = DomainMembership(
                domain=domain,
                iam_role_name=domain_iam_role_name,
                domain_ou=domain_ou,
                domain_fqdn=domain_fqdn,
                auth_secret_arn=kwargs.get("domain_auth_secret_arn"),
                dns_ips=kwargs.get("domain_dns_ips"),
            )
            self.domain_memberships.append(domain_membership)
        self.global_cluster_identifier = kwargs.get("global_cluster_identifier")
        if (
            not self.master_username
            and self.global_cluster_identifier
            or self.engine == "neptune"
        ):
            pass
        elif not self.master_username:
            raise InvalidParameterValue(
                "The parameter MasterUsername must be provided and must not be blank."
            )
        self.manage_master_user_password = manage_master_user_password
        if self.manage_master_user_password:
            self.master_user_secret = MasterUserSecret(
                self, master_user_secret_kms_key_id
            )
        elif not self.global_cluster_identifier and not self.engine == "neptune":
            self.master_user_password = master_user_password or ""
        self.availability_zones = kwargs.get("availability_zones")
        if not self.availability_zones:
            self.availability_zones = [
                f"{self.region}a",
                f"{self.region}b",
                f"{self.region}c",
            ]
        default_pg = (
            "default.neptune1.3" if self.engine == "neptune" else "default.aurora8.0"
        )
        self.parameter_group = db_cluster_parameter_group_name or default_pg
        self.db_subnet_group_name = db_subnet_group_name or "default"
        self.port = port or DBCluster.default_port(self.engine)
        self.preferred_backup_window = preferred_backup_window
        self.preferred_maintenance_window = preferred_maintenance_window
        self.publicly_accessible = publicly_accessible
        # This should default to the default security group
        self.vpc_security_group_ids = vpc_security_group_ids
        self.hosted_zone_id = "".join(
            random.choice(string.ascii_uppercase + string.digits) for _ in range(14)
        )
        self.db_cluster_resource_id = "cluster-" + "".join(
            random.choice(string.ascii_uppercase + string.digits) for _ in range(26)
        )
        self.tags = tags or []
        self.enabled_cloudwatch_logs_exports = (
            kwargs.get("enable_cloudwatch_logs_exports") or []
        )
        self.enable_http_endpoint = kwargs.get("enable_http_endpoint")  # type: ignore
        self.earliest_restorable_time = utcnow()
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
        self.replication_source_identifier = kwargs.get("replication_source_identifier")
        self.read_replica_identifiers: List[str] = list()
        self.is_writer: bool = False
        self.storage_encrypted = storage_encrypted
        self.kms_key_id = kms_key_id or "default_kms_key_id"
        if self.engine == "aurora-mysql" or self.engine == "aurora-postgresql":
            self._global_write_forwarding_requested = kwargs.get(
                "enable_global_write_forwarding"
            )
        self.backup_retention_period = backup_retention_period

        if backtrack := kwargs.get("backtrack_window"):
            if self.engine == "aurora-mysql":
                # https://docs.aws.amazon.com/AmazonRDS/latest/APIReference/API_CreateDBCluster.html
                if 0 <= backtrack <= 259200:
                    self.backtrack_window: int = backtrack
                else:
                    raise InvalidParameterValue(
                        f"The specified value ({backtrack}) is not a valid Backtrack Window. "
                        "Allowed values are within the range of 0 to 259200"
                    )
            else:
                raise InvalidParameterValue(
                    "Backtrack is not enabled for the postgres engine."
                )
        else:
            self.backtrack_window = 0

        self.iam_auth = kwargs.get("enable_iam_database_authentication", False)
        if self.iam_auth:
            if not self.engine.startswith("aurora-"):
                raise InvalidParameterCombination(
                    "IAM Authentication is currently not supported by Multi-AZ DB clusters."
                )
        self.license_model = license_model

    @property
    def db_cluster_identifier(self) -> str:
        return self._db_cluster_identifier

    @db_cluster_identifier.setter
    def db_cluster_identifier(self, value: str) -> None:
        self._db_cluster_identifier = value.lower()

    @property
    def db_subnet_group(self) -> str:
        # Despite the documentation saying this attribute returns:
        # "Information about the subnet group associated with the DB cluster,
        # including the name, description, and subnets in the subnet group."
        # It just returns the name...
        return self.db_subnet_group_name

    @property
    def resource_id(self) -> str:
        return self.db_cluster_identifier

    @property
    def multi_az(self) -> bool:
        availability_zones = list(
            set([instance.availability_zone for instance in self._members])
        )
        multi_az_conditions = [
            (len(self.read_replica_identifiers) > 0),
            (self.replication_source_identifier is not None),
            (len(availability_zones) > 1),
        ]
        return any(multi_az_conditions)

    @property
    def db_cluster_arn(self) -> str:
        return self.arn

    @property
    def latest_restorable_time(self) -> datetime:
        return utcnow()

    @property
    def master_user_password(self) -> str:
        raise NotImplementedError("Password not retrievable.")

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

    @property
    def http_endpoint_enabled(self) -> bool:
        return True if self.enable_http_endpoint else False

    @property
    def db_cluster_parameter_group(self) -> str:
        return self.parameter_group

    @property
    def status(self) -> str:
        if self._status == "creating":
            self._status = "available"
            return "creating"
        return self._status

    @status.setter
    def status(self, value: str) -> None:
        self._status = value

    @property
    def _members(self) -> List[DBInstanceClustered]:
        return [
            db_instance
            for db_instance in self.backend.databases.values()
            if isinstance(db_instance, DBInstanceClustered)
            and db_instance.db_cluster_identifier == self.resource_id
        ]

    @property
    def writer(self) -> Optional[DBInstanceClustered]:
        return next(
            (
                db_instance
                for db_instance in self._members
                if db_instance.is_cluster_writer
            ),
            None,
        )

    @writer.setter
    def writer(self, db_instance: DBInstanceClustered) -> None:
        db_instance.is_cluster_writer = True

    @property
    def members(self) -> List[DBInstanceClustered]:
        self.designate_writer()
        return self._members

    @property
    def endpoint(self) -> str:
        return f"{self.db_cluster_identifier}.cluster-{self.url_identifier}.{self.region}.rds.amazonaws.com"

    @property
    def reader_endpoint(self) -> str:
        return f"{self.db_cluster_identifier}.cluster-ro-{self.url_identifier}.{self.region}.rds.amazonaws.com"

    def designate_writer(self) -> None:
        if self.writer or not self._members:
            return
        promotion_list = sorted(self._members, key=lambda x: x.promotion_tier)
        self.writer = promotion_list[0]

    @property
    def associated_roles(self) -> List[Dict[str, Any]]:  # type: ignore[misc]
        return []

    @property
    def scaling_configuration_info(self) -> Dict[str, Any]:  # type: ignore[misc]
        configuration = self.scaling_configuration or {}
        info = {
            "MinCapacity": configuration.get("min_capacity"),
            "MaxCapacity": configuration.get("max_capacity"),
            "AutoPause": configuration.get("auto_pause"),
            "SecondsUntilAutoPause": configuration.get("seconds_until_auto_pause"),
            "TimeoutAction": configuration.get("timeout_action"),
            "SecondsBeforeTimeout": configuration.get("seconds_before_timeout"),
        }
        return info

    @property
    def vpc_security_groups(self) -> List[Dict[str, Any]]:  # type: ignore[misc]
        groups = [
            {"VpcSecurityGroupId": sg_id, "Status": "active"}
            for sg_id in self.vpc_security_group_ids
        ]
        return groups

    @property
    def vpc_security_group_ids(self) -> List[str]:
        return self._vpc_security_group_ids

    @vpc_security_group_ids.setter
    def vpc_security_group_ids(
        self, vpc_security_group_ids: Optional[List[str]]
    ) -> None:
        if vpc_security_group_ids is None:
            vpc_security_group_ids = []
        self._vpc_security_group_ids = vpc_security_group_ids

    @property
    def cross_account_clone(self) -> bool:
        return False

    @property
    def global_write_forwarding_requested(self) -> bool:
        # This does not appear to be in the standard response for any clusters
        # Docs say it's only for a secondary cluster in aurora global database...
        return True if self._global_write_forwarding_requested else False

    @property
    def db_cluster_members(self) -> List[Dict[str, Any]]:  # type: ignore[misc]
        members = [
            {
                "DBInstanceIdentifier": member.db_instance_identifier,
                "IsClusterWriter": member.is_cluster_writer,
                "DBClusterParameterGroupStatus": "in-sync",
                "PromotionTier": member.promotion_tier,
            }
            for member in self.members
        ]
        return members

    @property
    def iam_database_authentication_enabled(self) -> bool:
        return True if self.iam_auth else False

    def get_cfg(self) -> Dict[str, Any]:
        cfg = self.__dict__.copy()
        cfg.pop("backend")
        cfg["master_user_password"] = cfg.pop("_master_user_password")
        cfg["enable_http_endpoint"] = cfg.pop("_enable_http_endpoint")
        cfg["vpc_security_group_ids"] = cfg.pop("_vpc_security_group_ids")
        return cfg

    @staticmethod
    def default_engine_version(engine: str) -> str:
        return {
            "aurora": "5.6.mysql_aurora.1.22.5",
            "aurora-mysql": "5.7.mysql_aurora.2.07.2",
            "aurora-postgresql": "12.7",
            "mysql": "8.0.23",
            "neptune": "1.3.2.1",
            "postgres": "13.4",
        }[engine]

    @staticmethod
    def default_port(engine: str) -> int:
        return {
            "aurora": 3306,
            "aurora-mysql": 3306,
            "aurora-postgresql": 5432,
            "mysql": 3306,
            "neptune": 8182,
            "postgres": 5432,
        }[engine]

    @staticmethod
    def default_storage_type(iops: Any) -> str:  # type: ignore[misc]
        return "gp2" if iops is None else "io1"

    @staticmethod
    def default_allocated_storage(engine: str, storage_type: str) -> int:
        return {
            "aurora": {"gp2": 1, "io1": 1, "standard": 1},
            "aurora-mysql": {"gp2": 1, "io1": 1, "standard": 1},
            "aurora-postgresql": {"gp2": 1, "io1": 1, "standard": 1},
            "mysql": {"gp2": 20, "io1": 100, "standard": 5},
            "neptune": {"gp2": 0, "io1": 0, "standard": 0},
            "postgres": {"gp2": 20, "io1": 100, "standard": 5},
        }[engine][storage_type]

    def failover(self, target_instance: DBInstanceClustered) -> None:
        if self.writer is not None:
            self.writer.is_cluster_writer = False
        self.writer = target_instance

    def save_automated_backup(self) -> None:
        time_stamp = utcnow().strftime("%Y-%m-%d-%H-%M")
        snapshot_id = f"rds:{self.db_cluster_identifier}-{time_stamp}"
        self.backend.create_auto_cluster_snapshot(
            self.db_cluster_identifier, snapshot_id
        )


class DBClusterSnapshot(SnapshotAttributesMixin, RDSBaseModel):
    resource_type = "cluster-snapshot"

    SUPPORTED_FILTERS = {
        "db-cluster-id": FilterDef(
            ["db_cluster_arn", "db_cluster_identifier"],
            "DB Cluster Identifiers",
            case_insensitive=True,
        ),
        "db-cluster-snapshot-id": FilterDef(
            ["db_cluster_snapshot_identifier"],
            "DB Cluster Snapshot Identifiers",
            case_insensitive=True,
        ),
        "snapshot-type": FilterDef(["snapshot_type"], "Snapshot Types"),
        "engine": FilterDef(["cluster.engine"], "Engine Names"),
    }

    def __init__(
        self,
        backend: RDSBackend,
        cluster: DBCluster,
        snapshot_id: str,
        snapshot_type: str = "manual",
        tags: Optional[List[Dict[str, str]]] = None,
        kms_key_id: Optional[str] = None,
        **kwargs: Any,
    ):
        super().__init__(backend=backend, **kwargs)
        self.db_cluster_snapshot_identifier = snapshot_id
        self.snapshot_type = snapshot_type
        self.percent_progress = 100
        self.status = "available"
        # If tags are provided at creation, AWS does *not* copy tags from the
        # db_cluster (even if copy_tags_to_snapshot is True).
        if tags is not None:
            self.tags = tags
        elif cluster.copy_tags_to_snapshot:
            self.tags = cluster.tags or []
        else:
            self.tags = []
        self.cluster = copy.copy(cluster)
        self.allocated_storage = self.cluster.allocated_storage
        self.cluster_create_time = self.cluster.created
        self.db_cluster_identifier = self.cluster.db_cluster_identifier
        if kms_key_id is not None:
            self.kms_key_id = self.cluster.kms_key_id = kms_key_id
            self.encrypted = self.cluster.storage_encrypted = True
        else:
            self.kms_key_id = self.cluster.kms_key_id
            self.encrypted = self.cluster.storage_encrypted
        self.engine = self.cluster.engine
        self.engine_version = self.cluster.engine_version
        self.master_username = self.cluster.master_username
        self.port = self.cluster.port
        self.storage_encrypted = self.cluster.storage_encrypted

    @property
    def db_cluster_snapshot_identifier(self) -> str:
        return self._db_cluster_snapshot_identifier

    @db_cluster_snapshot_identifier.setter
    def db_cluster_snapshot_identifier(self, value: str) -> None:
        self._db_cluster_snapshot_identifier = value.lower()

    @property
    def resource_id(self) -> str:
        return self.db_cluster_snapshot_identifier

    @property
    def db_cluster_snapshot_arn(self) -> str:
        return self.arn

    @property
    def snapshot_create_time(self) -> datetime:
        return self.created


class LogFileManager:
    def __init__(self, engine: str) -> None:
        self.log_files = []
        filename = f"error/{engine}.log"
        if engine == "postgres":
            formatted_time = utcnow().strftime("%Y-%m-%d-%H")
            filename = f"error/postgresql.log.{formatted_time}"
        self.log_files.append(DBLogFile(filename))

    @property
    def files(self) -> List[DBLogFile]:
        return self.log_files


class DBLogFile:
    def __init__(self, name: str) -> None:
        self.log_file_name = name
        self.last_written = unix_time()
        self.size = 123


class DBInstance(EventMixin, CloudFormationModel, RDSBaseModel):
    SUPPORTED_FILTERS = {
        "db-cluster-id": FilterDef(
            ["db_cluster_identifier"], "DB Cluster Identifiers", case_insensitive=True
        ),
        "db-instance-id": FilterDef(
            ["db_instance_arn", "db_instance_identifier"],
            "DB Instance Identifiers",
            case_insensitive=True,
        ),
        "dbi-resource-id": FilterDef(["dbi_resource_id"], "Dbi Resource Ids"),
        "domain": FilterDef(None, ""),
        "engine": FilterDef(["engine"], "Engine Names"),
    }

    default_engine_versions = {
        "MySQL": "5.6.21",
        "mysql": "5.6.21",
        "neptune": "1.3.2.1",
        "oracle-ee": "11.2.0.4.v3",
        "oracle-se": "11.2.0.4.v3",
        "oracle-se1": "11.2.0.4.v3",
        "oracle-se2": "11.2.0.4.v3",
        "postgres": "9.3.3",
        "sqlserver-ee": "11.00.2100.60.v1",
        "sqlserver-ex": "11.00.2100.60.v1",
        "sqlserver-se": "11.00.2100.60.v1",
        "sqlserver-web": "11.00.2100.60.v1",
    }
    event_source_type = "db-instance"
    resource_type = "db"

    def __init__(
        self,
        backend: RDSBackend,
        db_instance_identifier: str,
        db_instance_class: str,
        engine: str,
        engine_version: Optional[str] = None,
        port: Optional[int] = None,
        allocated_storage: Optional[int] = None,
        max_allocated_storage: Optional[int] = None,
        backup_retention_period: int = 1,
        character_set_name: Optional[str] = None,
        auto_minor_version_upgrade: bool = True,
        db_name: Optional[str] = None,
        db_security_groups: Optional[List[str]] = None,
        db_subnet_group_name: Optional[str] = None,
        db_cluster_identifier: Optional[str] = None,
        db_parameter_group_name: Optional[str] = None,
        domain: Optional[str] = None,
        domain_iam_role_name: Optional[str] = None,
        domain_ou: Optional[str] = None,
        domain_fqdn: Optional[str] = None,
        copy_tags_to_snapshot: bool = False,
        iops: Optional[str] = None,
        master_username: Optional[str] = None,
        master_user_password: Optional[str] = None,
        multi_az: bool = False,
        license_model: str = "general-public-license",
        preferred_backup_window: str = "13:14-13:44",
        preferred_maintenance_window: str = "wed:06:38-wed:07:08",
        publicly_accessible: Optional[bool] = None,
        source_db_instance_identifier: Optional[str] = None,
        storage_type: Optional[str] = None,
        storage_encrypted: bool = False,
        tags: Optional[List[Dict[str, str]]] = None,
        vpc_security_group_ids: Optional[List[str]] = None,
        deletion_protection: bool = False,
        option_group_name: Optional[str] = None,
        enable_cloudwatch_logs_exports: Optional[List[str]] = None,
        ca_certificate_identifier: str = "rds-ca-default",
        availability_zone: Optional[str] = None,
        manage_master_user_password: Optional[bool] = False,
        master_user_secret_kms_key_id: Optional[str] = None,
        storage_throughput: Optional[int] = None,
        **kwargs: Any,
    ) -> None:
        super().__init__(backend)
        self.status = "available"
        self.is_replica = False
        self.replicas: List[str] = []
        self.engine = engine
        if self.engine not in DbInstanceEngine.valid_db_instance_engine():
            raise InvalidParameterValue(
                f"Value {self.engine} for parameter Engine is invalid. Reason: engine {self.engine} not supported"
            )
        self.log_file_manager = LogFileManager(self.engine)
        self.iops = iops
        self.auto_minor_version_upgrade = auto_minor_version_upgrade
        self.db_instance_identifier = db_instance_identifier
        self.source_db_instance_identifier = source_db_instance_identifier
        self.db_instance_class = db_instance_class
        self.port = port
        if self.port is None:
            self.port = DBInstance.default_port(self.engine)
        self.db_name = db_name
        self.instance_create_time = self.created
        self.publicly_accessible = publicly_accessible
        self.copy_tags_to_snapshot = copy_tags_to_snapshot
        self.availability_zone: str = availability_zone or f"{self.region}a"
        self.multi_az = multi_az
        self.db_subnet_group_name = db_subnet_group_name
        self.db_security_groups = db_security_groups or []
        self.domain_memberships: List[DomainMembership] = []
        if domain or domain_iam_role_name or domain_ou or domain_fqdn:
            domain_membership = DomainMembership(
                domain=domain,
                iam_role_name=domain_iam_role_name,
                domain_ou=domain_ou,
                domain_fqdn=domain_fqdn,
                auth_secret_arn=kwargs.get("domain_auth_secret_arn"),
                dns_ips=kwargs.get("domain_dns_ips"),
            )
            self.domain_memberships.append(domain_membership)
        self.preferred_maintenance_window = preferred_maintenance_window.lower()
        self.db_parameter_group_name = db_parameter_group_name
        if (
            self.db_parameter_group_name
            and not self.is_default_parameter_group(self.db_parameter_group_name)
            and self.db_parameter_group_name
            not in rds_backends[self.account_id][self.region].db_parameter_groups
        ):
            raise DBParameterGroupNotFoundError(self.db_parameter_group_name)
        self.ca_certificate_identifier = ca_certificate_identifier
        self.enable_iam_database_authentication = kwargs.get(
            "enable_iam_database_authentication", False
        )
        self.dbi_resource_id = "db-M5ENSHXFPU6XHZ4G4ZEI5QIO2U"
        self.tags = tags or []
        self.deletion_protection = deletion_protection
        self.enabled_cloudwatch_logs_exports = enable_cloudwatch_logs_exports or []
        self.db_cluster_identifier = db_cluster_identifier
        if self.db_cluster_identifier is None:
            self.vpc_security_group_ids = vpc_security_group_ids or []
            if not self.vpc_security_group_ids:
                ec2_backend = ec2_backends[self.account_id][self.region]
                default_vpc = ec2_backend.default_vpc
                default_sg = ec2_backend.get_default_security_group(default_vpc.id)
                self.vpc_security_group_ids.append(default_sg.id)  # type: ignore
            self.storage_type = storage_type or DBInstance.default_storage_type(
                iops=self.iops
            )
            self.allocated_storage = (
                allocated_storage
                or DBInstance.default_allocated_storage(
                    engine=self.engine, storage_type=self.storage_type
                )
            )
            self.max_allocated_storage = max_allocated_storage or self.allocated_storage
            self.storage_encrypted = storage_encrypted
            if self.storage_encrypted:
                self.kms_key_id = kwargs.get("kms_key_id", "default_kms_key_id")
            else:
                self.kms_key_id = None
            self.backup_retention_period = backup_retention_period
            self.character_set_name = character_set_name
            self.engine_version = (
                engine_version or self.default_engine_versions[self.engine]
            )
            self.license_model = license_model
            self.master_username = master_username
            self.manage_master_user_password = manage_master_user_password
            if self.manage_master_user_password:
                self.master_user_secret = MasterUserSecret(
                    self, master_user_secret_kms_key_id
                )
            else:
                self.master_user_password = master_user_password
            self.preferred_backup_window = preferred_backup_window
            msg = valid_preferred_maintenance_window(
                self.preferred_maintenance_window,
                self.preferred_backup_window,
            )
            if msg:
                raise RDSClientError("InvalidParameterValue", msg)
            self.option_group_supplied = option_group_name is not None
            if (
                option_group_name
                and option_group_name
                not in rds_backends[self.account_id][self.region].option_groups
            ):
                raise OptionGroupNotFoundFaultError(option_group_name)
            assert self.engine and self.engine_version
            semantic = self.engine_version.split(".")
            option_suffix = semantic[0]
            if len(semantic) > 1:
                option_suffix = option_suffix + "-" + semantic[1]
            default_option_group_name = f"default:{self.engine}-{option_suffix}"
            self.option_group_name = option_group_name or default_option_group_name
            self.storage_throughput = storage_throughput

    @property
    def db_instance_identifier(self) -> str:
        return self._db_instance_identifier

    @db_instance_identifier.setter
    def db_instance_identifier(self, value: str) -> None:
        self._db_instance_identifier = value.lower()

    @property
    def db_subnet_group_name(self) -> Optional[str]:
        raise NotImplementedError("write only property")

    @db_subnet_group_name.setter
    def db_subnet_group_name(self, value: Optional[str]) -> None:
        self._db_subnet_group_name = value
        if self._db_subnet_group_name is not None:
            self.db_subnet_group = rds_backends[self.account_id][
                self.region
            ].describe_db_subnet_groups(self._db_subnet_group_name)[0]

    @property
    def allocated_storage(self) -> int:
        return self._allocated_storage

    @allocated_storage.setter
    def allocated_storage(self, value: int) -> None:
        self._allocated_storage = value

    @property
    def backup_retention_period(self) -> int:
        return self._backup_retention_period

    @backup_retention_period.setter
    def backup_retention_period(self, value: int) -> None:
        self._backup_retention_period = value

    @property
    def character_set_name(self) -> Optional[str]:
        return self._character_set_name

    @character_set_name.setter
    def character_set_name(self, value: Optional[str]) -> None:
        self._character_set_name = value

    @property
    def db_name(self) -> Optional[str]:
        return self._db_name

    @db_name.setter
    def db_name(self, value: Optional[str]) -> None:
        self._db_name = value

    @property
    def engine_version(self) -> str:
        return self._engine_version

    @engine_version.setter
    def engine_version(self, value: str) -> None:
        self._engine_version = value

    @property
    def kms_key_id(self) -> Optional[str]:
        return self._kms_key_id

    @kms_key_id.setter
    def kms_key_id(self, value: Optional[str]) -> None:
        self._kms_key_id = value

    @property
    def license_model(self) -> str:
        return self._license_model

    @license_model.setter
    def license_model(self, value: str) -> None:
        self._license_model = value

    @property
    def master_username(self) -> Optional[str]:
        return self._master_username

    @master_username.setter
    def master_username(self, value: Optional[str]) -> None:
        self._master_username = value

    @property
    def master_user_password(self) -> Optional[str]:
        raise NotImplementedError("Password is not retrievable.")

    @master_user_password.setter
    def master_user_password(self, value: Optional[str]) -> None:
        self._master_user_password = value

    @property
    def max_allocated_storage(self) -> Optional[int]:
        if self._max_allocated_storage > self.allocated_storage:
            return self._max_allocated_storage
        return None

    @max_allocated_storage.setter
    def max_allocated_storage(self, value: int) -> None:
        if value < self.allocated_storage:
            raise InvalidParameterCombination(
                "Max storage size must be greater than storage size"
            )
        self._max_allocated_storage = value

    @property
    def option_group_name(self) -> str:
        return self._option_group_name

    @option_group_name.setter
    def option_group_name(self, value: str) -> None:
        self._option_group_name = value

    @property
    def preferred_backup_window(self) -> str:
        return self._preferred_backup_window

    @preferred_backup_window.setter
    def preferred_backup_window(self, value: str) -> None:
        self._preferred_backup_window = value

    @property
    def storage_encrypted(self) -> bool:
        return self._storage_encrypted

    @storage_encrypted.setter
    def storage_encrypted(self, value: bool) -> None:
        self._storage_encrypted = value

    @property
    def storage_type(self) -> str:
        return self._storage_type

    @storage_type.setter
    def storage_type(self, value: str) -> None:
        self._storage_type = value

    @property
    def resource_id(self) -> str:
        return self.db_instance_identifier

    @property
    def db_instance_arn(self) -> str:
        return self.arn

    @property
    def physical_resource_id(self) -> Optional[str]:
        return self.db_instance_identifier

    @property
    def latest_restorable_time(self) -> datetime:
        return utcnow()

    def db_parameter_groups(self) -> List[DBParameterGroup]:
        if not self.db_parameter_group_name or self.is_default_parameter_group(
            self.db_parameter_group_name
        ):
            (
                db_family,
                db_parameter_group_name,
            ) = self.default_db_parameter_group_details()
            description = f"Default parameter group for {db_family}"
            if db_parameter_group_name in self.backend.db_parameter_groups:
                return [self.backend.db_parameter_groups[db_parameter_group_name]]
            default_group = DBParameterGroup(
                backend=self.backend,
                db_parameter_group_name=db_parameter_group_name,
                db_parameter_group_family=db_family,
                description=description,
                tags=[],
            )
            self.backend.db_parameter_groups[db_parameter_group_name] = default_group
            return [default_group]

        return [self.backend.db_parameter_groups[self.db_parameter_group_name]]

    def is_default_parameter_group(self, param_group_name: str) -> bool:
        return param_group_name.startswith(f"default.{self.engine.lower()}")

    def default_db_parameter_group_details(self) -> Tuple[str, str]:
        assert self.engine and self.engine_version
        minor_engine_version = ".".join(str(self.engine_version).rsplit(".")[:-1])
        db_family = f"{self.engine.lower()}{minor_engine_version}"

        return db_family, f"default.{db_family}"

    @property
    def address(self) -> str:
        return (
            f"{self.db_instance_identifier}.aaaaaaaaaa.{self.region}.rds.amazonaws.com"
        )

    @property
    def vpc_security_group_membership_list(self) -> List[Dict[str, Any]]:  # type: ignore[misc]
        groups = [
            {
                "Status": "active",
                "VpcSecurityGroupId": id_,
            }
            for id_ in self.vpc_security_group_ids
        ]
        return groups

    @property
    def db_parameter_group_status_list(self) -> Any:  # type: ignore[misc]
        groups = self.db_parameter_groups()
        for group in groups:
            setattr(group, "ParameterApplyStatus", "in-sync")
        return groups

    @property
    def db_security_group_membership_list(self) -> List[Dict[str, Any]]:  # type: ignore[misc]
        groups = [
            {
                "Status": "active",
                "DBSecurityGroupName": group,
            }
            for group in self.db_security_groups
        ]
        return groups

    @property
    def endpoint(self) -> Dict[str, Any]:  # type: ignore[misc]
        return {
            "Address": self.address,
            "Port": self.port,
        }

    @property
    def option_group_memberships(self) -> List[Dict[str, Any]]:  # type: ignore[misc]
        groups = [
            {
                "OptionGroupName": self.option_group_name,
                "Status": "in-sync",
            }
        ]
        return groups

    @property
    def read_replica_db_instance_identifiers(self) -> List[str]:
        return [replica for replica in self.replicas]

    @property
    def db_instance_port(self) -> Optional[int]:
        return self.port

    @property
    def read_replica_source_db_instance_identifier(self) -> Optional[str]:
        return self.source_db_instance_identifier

    @property
    def iam_database_authentication_enabled(self) -> bool:
        return self.enable_iam_database_authentication

    @property
    def storage_throughput(self) -> Optional[int]:
        return self._storage_throughput

    @storage_throughput.setter
    def storage_throughput(self, value: Optional[int]) -> None:
        self._storage_throughput: Optional[int]
        if value and self.storage_type == "gp3":
            self._storage_throughput = value
        else:
            self._storage_throughput = None

    def add_replica(self, replica: DBInstance) -> None:
        if self.region != replica.region:
            # Cross Region replica
            self.replicas.append(replica.db_instance_arn)
        else:
            self.replicas.append(replica.db_instance_identifier)

    def remove_replica(self, replica: DBInstance) -> None:
        self.replicas.remove(replica.db_instance_identifier)

    def set_as_replica(self) -> None:
        self.is_replica = True
        self.replicas = []

    def update(
        self,
        manage_master_user_password: Optional[bool] = None,
        master_user_secret_kms_key_id: Optional[str] = None,
        rotate_master_user_password: Optional[bool] = None,
        **db_kwargs: Dict[str, Any],
    ) -> None:
        if manage_master_user_password is True:
            self.master_user_secret = MasterUserSecret(
                self, master_user_secret_kms_key_id
            )
            self.manage_master_user_password = True
        elif manage_master_user_password is False and hasattr(
            self, "master_user_secret"
        ):
            self.master_user_secret.delete_secret()
            del self.master_user_secret
            self.manage_master_user_password = False

        if rotate_master_user_password is True:
            self.master_user_secret.rotate_secret()

        for key, value in db_kwargs.items():
            if value is not None:
                setattr(self, key, value)

        cwl_exports = db_kwargs.get("cloudwatch_logs_export_configuration") or {}
        for exp in cwl_exports.get("DisableLogTypes", []):
            self.enabled_cloudwatch_logs_exports.remove(exp)
        self.enabled_cloudwatch_logs_exports.extend(
            cwl_exports.get("EnableLogTypes", [])
        )

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
            "neptune": 8182,
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
        return "gp2" if iops is None else "io1"

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
    ) -> DBInstance:
        properties = cloudformation_json["Properties"]

        db_security_groups = properties.get("DBSecurityGroups")
        if not db_security_groups:
            db_security_groups = []
        security_groups = [group.name for group in db_security_groups]
        db_subnet_group = properties.get("DBSubnetGroupName")
        db_subnet_group_name = db_subnet_group.name if db_subnet_group else None
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
            "db_security_groups": security_groups,
            "storage_encrypted": properties.get("StorageEncrypted"),
            "storage_type": properties.get("StorageType"),
            "tags": properties.get("Tags"),
            "vpc_security_group_ids": properties.get("VpcSecurityGroupIds", []),
        }

        rds_backend = rds_backends[account_id][region_name]
        source_db_identifier = properties.get("SourceDBInstanceIdentifier")
        if source_db_identifier:
            # Replica
            db_kwargs["source_db_instance_identifier"] = source_db_identifier
            database = rds_backend.create_db_instance_read_replica(db_kwargs)
        else:
            database = rds_backend.create_db_instance(db_kwargs)
        return database

    def delete(self, account_id: str, region_name: str) -> None:
        backend = rds_backends[account_id][region_name]
        backend.delete_db_instance(self.db_instance_identifier)

    def save_automated_backup(self) -> None:
        self.add_event("DB_INSTANCE_BACKUP_START")
        time_stamp = utcnow().strftime("%Y-%m-%d-%H-%M")
        snapshot_id = f"rds:{self.db_instance_identifier}-{time_stamp}"
        self.backend.create_auto_snapshot(self.db_instance_identifier, snapshot_id)
        self.add_event("DB_INSTANCE_BACKUP_FINISH")


class DBInstanceClustered(DBInstance):
    def __init__(
        self, db_cluster_identifier: str, promotion_tier: int = 1, **kwargs: Any
    ) -> None:
        super().__init__(db_cluster_identifier=db_cluster_identifier, **kwargs)
        if db_cluster_identifier not in self.backend.clusters:
            raise DBClusterNotFoundError(
                db_cluster_identifier,
                f"The source cluster could not be found or cannot be accessed: {db_cluster_identifier}",
            )
        self.cluster = self.backend.clusters[db_cluster_identifier]
        self.db_cluster_identifier: str = self.cluster.db_cluster_identifier
        self.is_cluster_writer = True if not self.cluster.members else False
        self.promotion_tier = promotion_tier
        self.manage_master_user_password = False

    @property
    def allocated_storage(self) -> int:
        return self.cluster.allocated_storage

    @allocated_storage.setter
    def allocated_storage(self, value: int) -> None:
        raise NotImplementedError("Not valid for clustered db instances.")

    @property
    def backup_retention_period(self) -> int:
        return self.cluster.backup_retention_period

    @backup_retention_period.setter
    def backup_retention_period(self, value: int) -> None:
        raise NotImplementedError("Not valid for clustered db instances.")

    @property
    def character_set_name(self) -> Optional[str]:
        return self.cluster.character_set_name

    @character_set_name.setter
    def character_set_name(self, value: Optional[str]) -> None:
        raise NotImplementedError("Not valid for clustered db instances.")

    # TODO: Need to understand better how this works with Aurora instances.
    # According to the boto3 documentation, `db_name` is valid:
    # "The name of the database to create when the primary instance
    # of the DB cluster is created. If this parameter isn't specified,
    # no database is created in the DB instance."
    # So does that mean the cluster.database_name and the instance.db_name
    # can differ?
    @property
    def db_name(self) -> Optional[str]:
        return self._db_name or self.cluster.database_name

    @db_name.setter
    def db_name(self, value: Optional[str]) -> None:
        self._db_name = value

    @property
    def engine_version(self) -> str:
        return self.cluster.engine_version

    @engine_version.setter
    def engine_version(self, value: str) -> None:
        raise NotImplementedError("Not valid for clustered db instances.")

    @property
    def kms_key_id(self) -> Optional[str]:
        return self.cluster.kms_key_id

    @kms_key_id.setter
    def kms_key_id(self, value: Optional[str]) -> None:
        raise NotImplementedError("Not valid for clustered db instances.")

    @property
    def license_model(self) -> str:
        return self.cluster.license_model

    @license_model.setter
    def license_model(self, value: str) -> None:
        raise NotImplementedError("Not valid for clustered db instances.")

    @property
    def master_username(self) -> Optional[str]:
        return self.cluster.master_username

    @master_username.setter
    def master_username(self, value: Optional[str]) -> None:
        raise NotImplementedError("Not valid for clustered db instances.")

    @property
    def master_user_password(self) -> Optional[str]:
        raise NotImplementedError("Password is not retrievable.")

    @master_user_password.setter
    def master_user_password(self, value: Optional[str]) -> None:
        raise NotImplementedError("Not valid for clustered db instances.")

    @property
    def max_allocated_storage(self) -> Optional[int]:
        return None

    @max_allocated_storage.setter
    def max_allocated_storage(self, value: int) -> None:
        raise NotImplementedError("Not valid for clustered db instances.")

    @property
    def option_group_name(self) -> str:
        return self.cluster.option_group_name

    @option_group_name.setter
    def option_group_name(self, value: str) -> None:
        raise NotImplementedError("Not valid for clustered db instances.")

    @property
    def preferred_backup_window(self) -> str:
        return self.cluster.preferred_backup_window

    @preferred_backup_window.setter
    def preferred_backup_window(self, value: str) -> None:
        raise NotImplementedError("Not valid for clustered db instances.")

    @property
    def storage_encrypted(self) -> bool:
        return self.cluster.storage_encrypted

    @storage_encrypted.setter
    def storage_encrypted(self, value: bool) -> None:
        raise NotImplementedError("Not valid for clustered db instances.")

    @property
    def storage_type(self) -> str:
        return "aurora"

    @storage_type.setter
    def storage_type(self, value: str) -> None:
        raise NotImplementedError("Not valid for clustered db instances.")

    @property
    def vpc_security_group_ids(self) -> List[str]:
        return self.cluster.vpc_security_group_ids

    @vpc_security_group_ids.setter
    def vpc_security_group_ids(self, value: List[str]) -> None:
        raise NotImplementedError("Not valid for clustered db instances.")


class DBSnapshot(EventMixin, SnapshotAttributesMixin, RDSBaseModel):
    event_source_type = "db-snapshot"
    resource_type = "snapshot"

    SUPPORTED_FILTERS = {
        "db-instance-id": FilterDef(
            ["database.db_instance_arn", "database.db_instance_identifier"],
            "DB Instance Identifiers",
            case_insensitive=True,
        ),
        "db-snapshot-id": FilterDef(
            ["db_snapshot_identifier"], "DB Snapshot Identifiers", case_insensitive=True
        ),
        "dbi-resource-id": FilterDef(["database.dbi_resource_id"], "Dbi Resource Ids"),
        "snapshot-type": FilterDef(["snapshot_type"], "Snapshot Types"),
        "engine": FilterDef(["database.engine"], "Engine Names"),
    }

    def __init__(
        self,
        backend: RDSBackend,
        database: DBInstance,
        snapshot_id: str,
        snapshot_type: str = "manual",
        tags: Optional[List[Dict[str, str]]] = None,
        original_created_at: Optional[datetime] = None,
        kms_key_id: Optional[str] = None,
        **kwargs: Any,
    ):
        super().__init__(backend=backend, **kwargs)
        self.database = copy.copy(database)  # TODO: Refactor this out.
        self.db_snapshot_identifier = snapshot_id
        self.snapshot_type = snapshot_type
        self.status = "available"
        self.original_snapshot_create_time = original_created_at or self.created
        # If tags are provided at creation, AWS does *not* copy tags from the
        # db_cluster (even if copy_tags_to_snapshot is True).
        if tags:
            self.tags = tags
        elif database.copy_tags_to_snapshot and database.tags:
            self.tags = database.tags
        else:
            self.tags = []
        # Database attributes are captured at the time the snapshot is taken.
        self.allocated_storage = database.allocated_storage
        self.dbi_resource_id = database.dbi_resource_id
        self.db_instance_identifier = database.db_instance_identifier
        self.engine = database.engine
        self.engine_version = database.engine_version
        self.kms_key_id = kms_key_id or database.kms_key_id
        self.storage_encrypted = (
            self.kms_key_id is not None or database.storage_encrypted
        )
        self.encrypted = self.kms_key_id is not None and self.storage_encrypted
        self.iam_database_authentication_enabled = (
            database.enable_iam_database_authentication
        )
        self.instance_create_time = database.created
        self.master_username = database.master_username
        self.port = database.port

    @property
    def db_snapshot_identifier(self) -> str:
        return self._db_snapshot_identifier

    @db_snapshot_identifier.setter
    def db_snapshot_identifier(self, value: str) -> None:
        self._db_snapshot_identifier = value.lower()

    @property
    def resource_id(self) -> str:
        return self.db_snapshot_identifier

    @property
    def db_snapshot_arn(self) -> str:
        return self.arn

    @property
    def snapshot_create_time(self) -> datetime:
        return self.created


class ExportTask(RDSBaseModel):
    def __init__(
        self,
        backend: RDSBackend,
        snapshot: Union[DBSnapshot, DBClusterSnapshot],
        kwargs: Dict[str, Any],
    ):
        super().__init__(backend)
        self.snapshot = snapshot

        self.export_task_identifier = kwargs.get("export_task_identifier")
        self.kms_key_id = kwargs.get("kms_key_id", "default_kms_key_id")
        self.source_arn = kwargs.get("source_arn")
        self.iam_role_arn = kwargs.get("iam_role_arn")
        self.s3_bucket = kwargs.get("s3_bucket_name")
        self.s3_prefix = kwargs.get("s3_prefix", "")
        self.export_only = kwargs.get("export_only", [])

        self.status = "complete"
        self.created_at = utcnow()
        self.source_type = "SNAPSHOT" if type(snapshot) is DBSnapshot else "CLUSTER"


class EventSubscription(RDSBaseModel):
    resource_type = "es"

    def __init__(self, backend: RDSBackend, subscription_name: str, **kwargs: Any):
        super().__init__(backend)
        self.subscription_name = subscription_name
        self.sns_topic_arn = kwargs.get("sns_topic_arn")
        self.source_type = kwargs.get("source_type")
        self.event_categories = kwargs.get("event_categories", [])
        self.source_ids = kwargs.get("source_ids", [])
        self.enabled = kwargs.get("enabled", False)
        self.tags = kwargs.get("tags", [])
        self.status = "active"
        self.created_at = utcnow()

    @property
    def resource_id(self) -> str:
        return self.subscription_name

    @property
    def cust_subscription_id(self) -> str:
        return self.subscription_name

    @property
    def event_categories_list(self) -> List[str]:
        return self.event_categories

    @property
    def source_ids_list(self) -> List[str]:
        return self.source_ids


class DBSecurityGroup(CloudFormationModel, RDSBaseModel):
    resource_type = "secgrp"

    def __init__(
        self,
        backend: RDSBackend,
        group_name: str,
        description: str,
        tags: List[Dict[str, str]],
    ):
        super().__init__(backend)
        self.name = group_name
        self.description = description
        self.status = "authorized"
        self._ip_ranges: List[Any] = []
        self._ec2_security_groups: List[Any] = []
        self.tags = tags
        self.vpc_id = None

    @property
    def resource_id(self) -> str:
        return self.name

    @property
    def ec2_security_groups(self) -> List[Dict[str, str]]:
        security_groups = [
            {
                "Status": "Active",
                "EC2SecurityGroupName": sg.name,
                "EC2SecurityGroupId": sg.id,
                "EC2SecurityGroupOwnerId": sg.owner_id,
            }
            for sg in self._ec2_security_groups
        ]
        return security_groups

    @property
    def ip_ranges(self) -> List[Dict[str, Any]]:  # type: ignore[misc]
        ranges = [
            {
                "CIDRIP": ip_range,
                "Status": "authorized",
            }
            for ip_range in self._ip_ranges
        ]
        return ranges

    def authorize_cidr(self, cidr_ip: str) -> None:
        self._ip_ranges.append(cidr_ip)

    def authorize_security_group(self, security_group: str) -> None:
        self._ec2_security_groups.append(security_group)

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
    ) -> DBSecurityGroup:
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

    def delete(self, account_id: str, region_name: str) -> None:
        backend = rds_backends[account_id][region_name]
        backend.delete_security_group(self.name)


class DBSubnetGroup(CloudFormationModel, RDSBaseModel):
    resource_type = "subgrp"

    def __init__(
        self,
        backend: RDSBackend,
        subnet_name: str,
        description: str,
        subnets: List[Subnet],
        tags: List[Dict[str, str]],
    ):
        super().__init__(backend)
        self.name = subnet_name
        self.description = description
        self._subnets = subnets
        self.subnet_group_status = "Complete"
        self.tags = tags
        self.vpc_id = self._subnets[0].vpc_id

    @property
    def name(self) -> str:
        return self._name

    @name.setter
    def name(self, value: str) -> None:
        self._name = value.lower()

    @property
    def resource_id(self) -> str:
        return self.name

    @property
    def db_subnet_group_description(self) -> str:
        return self.description

    @property
    def subnets(self) -> List[Dict[str, Any]]:  # type: ignore[misc]
        subnets = [
            {
                "SubnetStatus": "Active",
                "SubnetIdentifier": subnet.id,
                "SubnetAvailabilityZone": {
                    "Name": subnet.availability_zone,
                    "ProvisionedIopsCapable": False,
                },
            }
            for subnet in self._subnets
        ]
        return subnets

    @subnets.setter
    def subnets(self, subnets: List[Subnet]) -> None:
        self._subnets = subnets

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
    ) -> DBSubnetGroup:
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

    def delete(self, account_id: str, region_name: str) -> None:
        backend = rds_backends[account_id][region_name]
        backend.delete_subnet_group(self.name)


class DBProxy(RDSBaseModel):
    resource_type = "db-proxy"

    def __init__(
        self,
        backend: RDSBackend,
        db_proxy_name: str,
        engine_family: str,
        auth: List[Dict[str, str]],
        role_arn: str,
        vpc_subnet_ids: List[str],
        vpc_security_group_ids: Optional[List[str]],
        require_tls: Optional[bool] = False,
        idle_client_timeout: Optional[int] = 1800,
        debug_logging: Optional[bool] = False,
        tags: Optional[List[Dict[str, str]]] = None,
    ):
        super().__init__(backend)
        self.db_proxy_name = db_proxy_name
        self.engine_family = engine_family
        if self.engine_family not in ["MYSQL", "POSTGRESQL", "SQLSERVER"]:
            raise InvalidParameterValue("Provided EngineFamily is not valid.")
        self.auth = auth
        self.role_arn = role_arn
        self.vpc_subnet_ids = vpc_subnet_ids
        self.vpc_security_group_ids = vpc_security_group_ids or []
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
        if self.debug_logging is None:
            self.debug_logging = False
        self.created_date = utcnow()
        self.updated_date = utcnow()
        if tags is None:
            self.tags = []
        else:
            self.tags = tags
        ec2_backend = ec2_backends[self.account_id][self.region]
        subnets = ec2_backend.describe_subnets(subnet_ids=self.vpc_subnet_ids)
        vpcs = []
        for subnet in subnets:
            vpcs.append(subnet.vpc_id)
            if subnet.vpc_id != vpcs[0]:
                raise InvalidSubnet(subnet_identifier=subnet.id)
        if not self.vpc_security_group_ids:
            default_sg = ec2_backend.get_default_security_group(vpcs[0])
            self.vpc_security_group_ids.append(default_sg.id)  # type: ignore

        self.vpc_id = ec2_backend.describe_subnets(subnet_ids=[self.vpc_subnet_ids[0]])[
            0
        ].vpc_id
        self.status = "available"
        self.url_identifier = "".join(
            random.choice(string.ascii_lowercase + string.digits) for _ in range(12)
        )
        self.endpoint = f"{self.db_proxy_name}.db-proxy-{self.url_identifier}.{self.region}.rds.amazonaws.com"

        self.proxy_target_groups = {
            "default": DBProxyTargetGroup(
                backend=self.backend, name="default", proxy_name=db_proxy_name
            )
        }

        self.unique_id = f"prx-{random.get_random_string(17, lower_case=True)}"

    @property
    def resource_id(self) -> str:
        return self.unique_id


class DBInstanceAutomatedBackup:
    def __init__(
        self,
        backend: RDSBackend,
        db_instance_identifier: str,
        automated_snapshots: List[DBSnapshot],
    ) -> None:
        self.backend = backend
        self.db_instance_identifier = db_instance_identifier
        self.automated_snapshots = automated_snapshots

    @property
    def status(self) -> str:
        status = "active"
        if self.db_instance_identifier not in self.backend.databases:
            status = "retained"
        return status


class DBShardGroup(RDSBaseModel):
    resource_type = "shard-group"

    SUPPORTED_FILTERS = {
        "db-shard-group-id": FilterDef(
            ["db_shard_group_identifier"],
            "DB Shard Group Identifier",
        ),
        "db-cluster-id": FilterDef(
            ["db_cluster_identifier"],
            "DB Cluster Identifier",
        ),
    }

    def __init__(
        self,
        backend: RDSBackend,
        db_shard_group_identifier: str,
        db_cluster_identifier: str,
        max_acu: float,
        compute_redundancy: Optional[int] = None,
        min_acu: Optional[float] = None,
        publicly_accessible: Optional[bool] = None,
        tags: Optional[List[Dict[str, str]]] = None,
    ):
        super().__init__(backend)
        self.db_shard_group_identifier = db_shard_group_identifier
        self.db_cluster_identifier = db_cluster_identifier
        self.compute_redundancy = compute_redundancy
        self.max_acu = max_acu
        self.min_acu = min_acu
        self.status = "available"
        self.db_shard_group_resource_id = (
            f"shard-group-{self.db_shard_group_identifier}"
        )
        self.endpoint = f"{self.db_cluster_identifier}.{self.region}.rds.amazonaws.com"
        self.db_shard_group_arn = f"arn:aws:rds:{self.region}:{self.account_id}:shard-group:{self.db_shard_group_identifier}"
        self.publicly_accessible = publicly_accessible
        self.tags = tags or []


class RDSBackend(BaseBackend):
    def __init__(self, region_name: str, account_id: str):
        super().__init__(region_name, account_id)
        self.arn_regex = re_compile(
            ARN_PARTITION_REGEX
            + r":rds:.*:[0-9]*:(db|cluster|es|og|pg|ri|secgrp|snapshot|cluster-snapshot|subgrp|db-proxy):.*$"
        )
        self.clusters: MutableMapping[str, DBCluster] = CaseInsensitiveDict()
        self.global_clusters: MutableMapping[str, GlobalCluster] = CaseInsensitiveDict()
        self.databases: MutableMapping[str, DBInstance] = CaseInsensitiveDict()
        self.database_snapshots: MutableMapping[str, DBSnapshot] = CaseInsensitiveDict()
        self.cluster_snapshots: MutableMapping[str, DBClusterSnapshot] = (
            CaseInsensitiveDict()
        )
        self.export_tasks: Dict[str, ExportTask] = OrderedDict()
        self.event_subscriptions: Dict[str, EventSubscription] = OrderedDict()
        self.db_parameter_groups: MutableMapping[str, DBParameterGroup] = (
            CaseInsensitiveDict()
        )
        self.db_cluster_parameter_groups: MutableMapping[
            str, DBClusterParameterGroup
        ] = CaseInsensitiveDict()
        self.option_groups: MutableMapping[str, OptionGroup] = CaseInsensitiveDict()
        self.security_groups: Dict[str, DBSecurityGroup] = {}
        self.shard_groups: Dict[str, DBShardGroup] = {}
        self.subnet_groups: MutableMapping[str, DBSubnetGroup] = CaseInsensitiveDict()
        self._db_cluster_options: Optional[List[Dict[str, Any]]] = None
        self.db_proxies: Dict[str, DBProxy] = OrderedDict()
        self.events: List[Event] = []
        self.resource_map = {
            DBCluster: self.clusters,
            DBClusterParameterGroup: self.db_cluster_parameter_groups,
            DBClusterSnapshot: self.cluster_snapshots,
            DBInstance: self.databases,
            DBParameterGroup: self.db_parameter_groups,
            DBProxy: self.db_proxies,
            DBSecurityGroup: self.security_groups,
            DBShardGroup: self.shard_groups,
            DBSnapshot: self.database_snapshots,
            DBSubnetGroup: self.subnet_groups,
            EventSubscription: self.event_subscriptions,
            ExportTask: self.export_tasks,
            GlobalCluster: self.global_clusters,
            OptionGroup: self.option_groups,
        }
        self.blue_green_deployments: MutableMapping[str, BlueGreenDeployment] = (
            CaseInsensitiveDict()
        )

    @property
    def kms(self) -> KmsBackend:
        return kms_backends[self.account_id][self.region_name]

    @property
    def secretsmanager(self) -> SecretsManagerBackend:
        return self.get_backend("secretsmanager", self.region_name, self.account_id)

    def _validate_kms_key(self, kms_key_id: str) -> str:
        key = kms_key_id
        kms_backend = self.kms
        if (match := re.fullmatch(KMS_ARN_PATTERN, kms_key_id)) is not None:
            region = match.group("region")
            if region != self.region_name:
                raise KMSKeyNotAccessibleFault(kms_key_id)
            account = match.group("account_id")
            key = match.group("key")
            kms_backend = self.get_backend("kms", region, account)
            assert isinstance(kms_backend, KmsBackend)
        try:
            kms_key = kms_backend.describe_key(key)
        except KeyError:
            raise KMSKeyNotAccessibleFault(kms_key_id)
        validated_key = kms_key.arn
        return validated_key

    @overload
    def get_backend(
        self,
        service: Literal["kms"],
        region: str,
        account_id: Optional[str] = None,
    ) -> KmsBackend: ...

    @overload
    def get_backend(
        self,
        service: Literal["rds"],
        region: str,
        account_id: Optional[str] = None,
    ) -> RDSBackend: ...

    @overload
    def get_backend(
        self,
        service: Literal["secretsmanager"],
        region: str,
        account_id: Optional[str] = None,
    ) -> SecretsManagerBackend: ...

    def get_backend(
        self,
        service: Literal["kms"] | Literal["rds"] | Literal["secretsmanager"],
        region: str,
        account_id: Optional[str] = None,
    ) -> KmsBackend | RDSBackend | SecretsManagerBackend:
        from moto.backends import get_backend as get_moto_backend

        if account_id is None:
            account_id = self.account_id

        return get_moto_backend(service)[account_id][region]

    @overload
    def get_snapshot(
        self,
        identifier: str,
        resource_type: type[DBSnapshot],
        not_found_exception: type[DBSnapshotNotFoundFault],
    ) -> DBSnapshot: ...

    @overload
    def get_snapshot(
        self,
        identifier: str,
        resource_type: type[DBClusterSnapshot],
        not_found_exception: type[DBClusterSnapshotNotFoundError],
    ) -> DBClusterSnapshot: ...

    def get_snapshot(
        self,
        identifier: str,
        resource_type: type[DBSnapshot] | type[DBClusterSnapshot],
        not_found_exception: type[DBSnapshotNotFoundFault]
        | type[DBClusterSnapshotNotFoundError],
    ) -> DBSnapshot | DBClusterSnapshot:
        region = self.region_name
        if identifier.startswith("arn"):
            region = identifier.split(":")[3]
            identifier = identifier.split(":")[-1]
        backend = self.get_backend("rds", region=region)
        snapshots = backend.resource_map[resource_type]
        if identifier not in snapshots:
            raise not_found_exception(identifier)
        return snapshots[identifier]

    def get_db_snapshot(self, identifier: str) -> DBSnapshot:
        return self.get_snapshot(
            identifier,
            resource_type=DBSnapshot,
            not_found_exception=DBSnapshotNotFoundFault,
        )

    def get_db_cluster_snapshot(self, identifier: str) -> DBClusterSnapshot:
        return self.get_snapshot(
            identifier,
            resource_type=DBClusterSnapshot,
            not_found_exception=DBClusterSnapshotNotFoundError,
        )

    def get_shared_snapshots(
        self, resource_type: type[DBSnapshot] | type[DBClusterSnapshot]
    ) -> List[DBSnapshot | DBClusterSnapshot]:
        snapshots_shared = []
        for backend_container in rds_backends.values():
            for backend in backend_container.values():
                if backend.region_name != self.region_name:
                    continue
                snapshots = backend.resource_map[resource_type].values()
                for snapshot in snapshots:
                    if self.account_id in snapshot.attributes["restore"]:
                        snapshots_shared.append(snapshot)
        return snapshots_shared

    get_shared_db_snapshots = partialmethod(
        get_shared_snapshots, resource_type=DBSnapshot
    )
    get_shared_db_cluster_snapshots = partialmethod(
        get_shared_snapshots, resource_type=DBClusterSnapshot
    )

    @lru_cache()
    def db_cluster_options(self, engine) -> List[Dict[str, Any]]:  # type: ignore
        from moto.rds.utils import decode_orderable_db_instance

        decoded_options = load_resource(
            __name__, f"resources/cluster_options/{engine}.json"
        )
        self._db_cluster_options = [
            decode_orderable_db_instance(option) for option in decoded_options
        ]
        return self._db_cluster_options

    def create_db_instance(self, db_kwargs: Dict[str, Any]) -> DBInstance:
        database_id = db_kwargs["db_instance_identifier"]
        if database_id in self.databases:
            raise DBInstanceAlreadyExists()
        self._validate_db_identifier(database_id)
        if db_kwargs.get("db_cluster_identifier") is None:
            database = DBInstance(self, **db_kwargs)
        else:
            database = DBInstanceClustered(backend=self, **db_kwargs)

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
        self.databases[database_id] = database
        database.add_event("DB_INSTANCE_CREATE")
        database.save_automated_backup()
        return database

    def create_auto_snapshot(
        self,
        db_instance_identifier: str,
        db_snapshot_identifier: str,
    ) -> DBSnapshot:
        snapshot = self.create_db_snapshot(
            db_instance_identifier, db_snapshot_identifier, snapshot_type="automated"
        )
        snapshot.add_event("DB_SNAPSHOT_CREATE_AUTOMATED_START")
        snapshot.add_event("DB_SNAPSHOT_CREATE_AUTOMATED_FINISH")
        return snapshot

    def create_db_snapshot(
        self,
        db_instance_identifier: str,
        db_snapshot_identifier: str,
        snapshot_type: str = "manual",
        tags: Optional[List[Dict[str, str]]] = None,
        original_created_at: Optional[datetime] = None,
        kms_key_id: Optional[str] = None,
    ) -> DBSnapshot:
        database = self.databases.get(db_instance_identifier)
        if not database:
            raise DBInstanceNotFoundError(db_instance_identifier)

        if db_snapshot_identifier in self.database_snapshots:
            raise DBSnapshotAlreadyExistsError(db_snapshot_identifier)
        if len(self.database_snapshots) >= int(
            os.environ.get("MOTO_RDS_SNAPSHOT_LIMIT", "100")
        ):
            raise SnapshotQuotaExceededFault()
        snapshot = DBSnapshot(
            self,
            database,
            db_snapshot_identifier,
            snapshot_type,
            tags,
            original_created_at,
            kms_key_id,
        )
        self.database_snapshots[db_snapshot_identifier] = snapshot
        return snapshot

    def copy_db_snapshot(
        self,
        source_db_snapshot_identifier: str,
        target_db_snapshot_identifier: str,
        kms_key_id: Optional[str] = None,
        tags: Optional[List[Dict[str, str]]] = None,
        copy_tags: Optional[bool] = False,
        **_: Any,
    ) -> DBSnapshot:
        if source_db_snapshot_identifier.startswith("arn:aws:rds:"):
            self.extract_snapshot_name_from_arn(source_db_snapshot_identifier)
        if target_db_snapshot_identifier in self.database_snapshots:
            raise DBSnapshotAlreadyExistsError(target_db_snapshot_identifier)

        if len(self.database_snapshots) >= int(
            os.environ.get("MOTO_RDS_SNAPSHOT_LIMIT", "100")
        ):
            raise SnapshotQuotaExceededFault()
        if kms_key_id is not None:
            kms_key_id = self._validate_kms_key(kms_key_id)
        source_db_snapshot = self.get_db_snapshot(source_db_snapshot_identifier)
        # When tags are passed, AWS does NOT copy/merge tags of the
        # source snapshot, even when copy_tags=True is given.
        # But when tags=[], AWS does honor copy_tags=True.
        if copy_tags and not tags:
            tags = source_db_snapshot.tags or []
        target_db_snapshot = DBSnapshot(
            self,
            source_db_snapshot.database,
            target_db_snapshot_identifier,
            tags=tags,
            kms_key_id=kms_key_id,
            original_created_at=source_db_snapshot.original_snapshot_create_time,
        )
        self.database_snapshots[target_db_snapshot_identifier] = target_db_snapshot
        return target_db_snapshot

    def delete_db_snapshot(self, db_snapshot_identifier: str) -> DBSnapshot:
        if db_snapshot_identifier not in self.database_snapshots:
            raise DBSnapshotNotFoundFault(db_snapshot_identifier)

        return self.database_snapshots.pop(db_snapshot_identifier)

    def promote_read_replica(self, db_kwargs: Dict[str, Any]) -> DBInstance:
        database_id = db_kwargs["db_instance_identifier"]
        database = self.databases[database_id]
        if database.is_replica:
            database.is_replica = False
            database.update(**db_kwargs)

        return database

    def create_db_instance_read_replica(self, db_kwargs: Dict[str, Any]) -> DBInstance:
        database_id = db_kwargs["db_instance_identifier"]
        source_database_id = db_kwargs["source_db_instance_identifier"]
        primary = self.find_db_from_id(source_database_id)
        if self.arn_regex.match(source_database_id):
            db_kwargs["backend"] = self

        # Shouldn't really copy here as the instance is duplicated. RDS replicas have different instances.
        replica = copy.copy(primary)
        replica.update(**db_kwargs)
        replica.set_as_replica()
        self.databases[database_id] = replica
        primary.add_replica(replica)
        return replica

    def describe_db_instances(
        self, db_instance_identifier: Optional[str] = None, filters: Any = None
    ) -> List[DBInstance]:
        databases = self.databases
        if db_instance_identifier:
            filters = merge_filters(
                filters, {"db-instance-id": [db_instance_identifier]}
            )
        if filters:
            databases = self._filter_resources(databases, filters, DBInstance)
        if db_instance_identifier and not databases:
            raise DBInstanceNotFoundError(db_instance_identifier)
        return list(databases.values())

    def describe_db_snapshots(
        self,
        db_instance_identifier: Optional[str],
        db_snapshot_identifier: Optional[str] = None,
        snapshot_type: Optional[str] = None,
        filters: Optional[Dict[str, Any]] = None,
    ) -> List[DBSnapshot]:
        if snapshot_type == "shared":
            return self.get_shared_db_snapshots()  # type: ignore[return-value]
        snapshots = self.database_snapshots
        if db_instance_identifier:
            filters = merge_filters(
                filters, {"db-instance-id": [db_instance_identifier]}
            )
        if db_snapshot_identifier:
            filters = merge_filters(
                filters, {"db-snapshot-id": [db_snapshot_identifier]}
            )
        snapshot_types = (
            ["automated", "manual"]
            if (
                snapshot_type is None
                and (filters is not None and "snapshot-type" not in filters)
            )
            else [snapshot_type]
            if snapshot_type is not None
            else []
        )
        if snapshot_types:
            filters = merge_filters(filters, {"snapshot-type": snapshot_types})
        if filters:
            snapshots = self._filter_resources(snapshots, filters, DBSnapshot)
        if db_snapshot_identifier and not snapshots and not db_instance_identifier:
            raise DBSnapshotNotFoundFault(db_snapshot_identifier)
        return list(snapshots.values())

    def modify_db_instance(
        self, db_instance_identifier: str, db_kwargs: Dict[str, Any]
    ) -> DBInstance:
        database = self.describe_db_instances(db_instance_identifier)[0]
        if "new_db_instance_identifier" in db_kwargs:
            del self.databases[db_instance_identifier]
            db_instance_identifier = db_kwargs["db_instance_identifier"] = (
                db_kwargs.pop("new_db_instance_identifier")
            )
            self.databases[db_instance_identifier] = database
        if "db_parameter_group_name" in db_kwargs:
            db_parameter_group_name = db_kwargs["db_parameter_group_name"]
            if db_parameter_group_name not in self.db_parameter_groups:
                raise DBParameterGroupNotFoundError(db_parameter_group_name)
        preferred_backup_window = db_kwargs.get(
            "preferred_backup_window", database.preferred_backup_window
        )
        preferred_maintenance_window = db_kwargs.get(
            "preferred_maintenance_window", database.preferred_maintenance_window
        )
        if preferred_maintenance_window or preferred_backup_window:
            msg = valid_preferred_maintenance_window(
                preferred_maintenance_window, preferred_backup_window
            )
            if msg:
                raise RDSClientError("InvalidParameterValue", msg)

        if db_kwargs.get("rotate_master_user_password") and not db_kwargs.get(
            "apply_immediately"
        ):
            raise RDSClientError(
                "InvalidParameterCombination",
                "You must specify apply immediately when rotating the master user password.",
            )
        database.update(**db_kwargs)
        return database

    def reboot_db_instance(self, db_instance_identifier: str) -> DBInstance:
        return self.describe_db_instances(db_instance_identifier)[0]

    def extract_snapshot_name_from_arn(self, snapshot_arn: str) -> str:
        arn_breakdown = snapshot_arn.split(":")
        region_name, account_id, resource_type, snapshot_name = arn_breakdown[3:7]
        if resource_type != "snapshot":
            raise InvalidParameterValue(
                "The parameter SourceDBSnapshotIdentifier is not a valid identifier. "
                "Identifiers must begin with a letter; must contain only ASCII "
                "letters, digits, and hyphens; and must not end with a hyphen or "
                "contain two consecutive hyphens."
            )
        return snapshot_name

    def restore_db_instance_from_db_snapshot(
        self, from_snapshot_id: str, overrides: Dict[str, Any]
    ) -> DBInstance:
        if from_snapshot_id.startswith("arn:aws:rds:"):
            from_snapshot_id = self.extract_snapshot_name_from_arn(from_snapshot_id)

        snapshot = self.describe_db_snapshots(
            db_instance_identifier=None, db_snapshot_identifier=from_snapshot_id
        )[0]
        original_database = snapshot.database

        if overrides["db_instance_identifier"] in self.databases:
            raise DBInstanceAlreadyExists()

        new_instance_props = {}
        for key, value in original_database.__dict__.items():
            if key.startswith("_"):
                key = key[1:]
            if key not in [
                "backend",
                "db_parameter_group_name",
                "vpc_security_group_ids",
                "max_allocated_storage",
            ]:
                new_instance_props[key] = copy.copy(value)
        new_instance_props["kms_key_id"] = snapshot.kms_key_id
        new_instance_props["storage_encrypted"] = snapshot.encrypted
        if not original_database.option_group_supplied:
            # If the option group is not supplied originally, the 'option_group_name' will receive a default value
            # Force this reconstruction, and prevent any validation on the default value
            del new_instance_props["option_group_name"]
        if "allocated_storage" in overrides:
            if overrides["allocated_storage"] < snapshot.allocated_storage:
                raise InvalidParameterValue(
                    "The allocated storage size can't be less than the source snapshot or backup size."
                )
        for key, value in overrides.items():
            if value:
                new_instance_props[key] = value

        return self.create_db_instance(new_instance_props)

    def restore_db_instance_to_point_in_time(
        self,
        source_db_identifier: str,
        target_db_identifier: str,
        overrides: Dict[str, Any],
    ) -> DBInstance:
        db_instance = self.describe_db_instances(
            db_instance_identifier=source_db_identifier
        )[0]

        new_instance_props = {}
        for key, value in db_instance.__dict__.items():
            if key.startswith("_"):
                key = key[1:]
            # Remove backend / db subnet group as they cannot be copied
            # and are not used in the restored instance.
            if key in ("backend", "db_subnet_group", "max_allocated_storage"):
                continue
            new_instance_props[key] = copy.deepcopy(value)
        new_instance_props["db_name"] = db_instance.db_name

        if not db_instance.option_group_supplied:
            # If the option group is not supplied originally, the 'option_group_name' will receive a default value
            # Force this reconstruction, and prevent any validation on the default value
            del new_instance_props["option_group_name"]
        if "allocated_storage" in overrides:
            if overrides["allocated_storage"] < db_instance.allocated_storage:
                raise InvalidParameterValue(
                    "Allocated storage size can't be less than the source instance size."
                )
        for key, value in overrides.items():
            if value:
                new_instance_props[key] = value

        # set the new db instance identifier
        new_instance_props["db_instance_identifier"] = target_db_identifier

        return self.create_db_instance(new_instance_props)

    def restore_db_cluster_to_point_in_time(
        self,
        db_cluster_identifier: str,
        source_db_cluster_identifier: str,
        restore_type: str = "full-copy",
        restore_to_time: Optional[datetime] = None,
        use_latest_restorable_time: bool = False,
        **overrides: Dict[str, Any],
    ) -> DBCluster:
        db_cluster = self.describe_db_clusters(
            db_cluster_identifier=source_db_cluster_identifier
        )[0]
        new_cluster_props = {}
        for key, value in db_cluster.__dict__.items():
            if key.startswith("_"):
                key = key[1:]
            # Remove backend / db subnet group as they cannot be copied
            # and are not used in the restored instance.
            if key in ("backend", "db_subnet_group", "vpc_security_group_ids"):
                continue
            new_cluster_props[key] = copy.copy(value)
        for key, value in overrides.items():
            new_cluster_props[key] = value
        new_cluster_props["db_cluster_identifier"] = db_cluster_identifier
        return self.create_db_cluster(new_cluster_props)

    def failover_db_cluster(
        self,
        db_cluster_identifier: str,
        target_db_instance_identifier: Optional[str] = None,
    ) -> DBCluster:
        target_instance = None
        if target_db_instance_identifier is not None:
            if target_db_instance_identifier not in self.databases:
                raise InvalidParameterValue(
                    f"Cannot find target instance: {target_db_instance_identifier}"
                )
            target_instance = self.databases[target_db_instance_identifier]
            if target_instance.status != "available":
                raise InvalidDBInstanceStateFault(
                    f"Instance {target_db_instance_identifier} should be in valid state for failover."
                )
        if db_cluster_identifier not in self.clusters:
            raise DBClusterNotFoundError(
                db_cluster_identifier,
                f"The source cluster could not be found or cannot be accessed: {db_cluster_identifier}",
            )
        cluster = self.clusters[db_cluster_identifier]
        if len(cluster.members) < 2:
            raise InvalidDBClusterStateFault(
                f"Database cluster:{cluster.db_cluster_identifier} should have at least two database instances for failover"
            )
        if target_instance is None:
            for instance in cluster.members:
                if not instance.is_cluster_writer and instance.status == "available":
                    target_instance = instance
                    break
            else:
                raise InvalidDBClusterStateFault(
                    f"Database cluster:{db_cluster_identifier} should have at least two database instances in valid states or a valid target instance for failover."
                )
        assert isinstance(target_instance, DBInstanceClustered)
        cluster.failover(target_instance)
        return cluster

    def stop_db_instance(
        self, db_instance_identifier: str, db_snapshot_identifier: Optional[str] = None
    ) -> DBInstance:
        self._validate_db_identifier(db_instance_identifier)
        database = self.describe_db_instances(db_instance_identifier)[0]
        # todo: certain rds types not allowed to be stopped at this time.
        # https://docs.aws.amazon.com/AmazonRDS/latest/UserGuide/USER_StopInstance.html#USER_StopInstance.Limitations
        if database.is_replica or (
            database.multi_az and database.engine.lower().startswith("sqlserver")
        ):
            # todo: more db types not supported by stop/start instance api
            raise InvalidDBClusterStateFaultError(db_instance_identifier)
        if database.status != "available":
            raise InvalidDBInstanceStateError(db_instance_identifier, "stop")
        if db_snapshot_identifier:
            self.create_db_snapshot(db_instance_identifier, db_snapshot_identifier)
        database.status = "stopped"
        return database

    def start_db_instance(self, db_instance_identifier: str) -> DBInstance:
        self._validate_db_identifier(db_instance_identifier)
        database = self.describe_db_instances(db_instance_identifier)[0]
        # todo: bunch of different error messages to be generated from this api call
        if database.status != "stopped":
            raise InvalidDBInstanceStateError(db_instance_identifier, "start")
        database.status = "available"
        return database

    def find_db_from_id(self, db_id: str) -> DBInstance:
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
        self,
        db_instance_identifier: str,
        final_db_snapshot_identifier: Optional[str] = None,
        skip_final_snapshot: Optional[bool] = False,
        delete_automated_backups: Optional[bool] = True,
    ) -> DBInstance:
        self._validate_db_identifier(db_instance_identifier)
        if db_instance_identifier in self.databases:
            if self.databases[db_instance_identifier].deletion_protection:
                raise InvalidParameterCombination(
                    "Cannot delete protected DB Instance, please disable deletion protection and try again."
                )
            if final_db_snapshot_identifier and not skip_final_snapshot:
                self.create_db_snapshot(
                    db_instance_identifier,
                    final_db_snapshot_identifier,
                    snapshot_type="manual",
                )
            database = self.databases.pop(db_instance_identifier)
            if database.manage_master_user_password:
                database.master_user_secret.delete_secret()
            if database.is_replica:
                primary = self.find_db_from_id(database.source_db_instance_identifier)  # type: ignore
                primary.remove_replica(database)
            automated_snapshots = self.describe_db_snapshots(
                db_instance_identifier,
                db_snapshot_identifier=None,
                snapshot_type="automated",
            )
            if delete_automated_backups:
                for snapshot in automated_snapshots:
                    self.delete_db_snapshot(snapshot.db_snapshot_identifier)
            database.status = "deleting"
            return database
        else:
            raise DBInstanceNotFoundError(db_instance_identifier)

    def create_db_security_group(
        self, group_name: str, description: str, tags: List[Dict[str, str]]
    ) -> DBSecurityGroup:
        security_group = DBSecurityGroup(self, group_name, description, tags)
        self.security_groups[group_name] = security_group
        return security_group

    def describe_security_groups(
        self, security_group_name: str
    ) -> List[DBSecurityGroup]:
        if security_group_name:
            if security_group_name in self.security_groups:
                return [self.security_groups[security_group_name]]
            else:
                raise DBSecurityGroupNotFoundError(security_group_name)
        return list(self.security_groups.values())

    def delete_security_group(self, security_group_name: str) -> DBSecurityGroup:
        if security_group_name in self.security_groups:
            return self.security_groups.pop(security_group_name)
        else:
            raise DBSecurityGroupNotFoundError(security_group_name)

    def delete_db_parameter_group(
        self, db_parameter_group_name: str
    ) -> DBParameterGroup:
        if db_parameter_group_name in self.db_parameter_groups:
            return self.db_parameter_groups.pop(db_parameter_group_name)
        else:
            raise DBParameterGroupNotFoundError(db_parameter_group_name)

    def authorize_security_group(
        self, security_group_name: str, cidr_ip: str
    ) -> DBSecurityGroup:
        security_group = self.describe_security_groups(security_group_name)[0]
        security_group.authorize_cidr(cidr_ip)
        return security_group

    def create_subnet_group(
        self,
        subnet_name: str,
        description: str,
        subnets: List[Any],
        tags: List[Dict[str, str]],
    ) -> DBSubnetGroup:
        subnet_group = DBSubnetGroup(self, subnet_name, description, subnets, tags)
        self.subnet_groups[subnet_name] = subnet_group
        return subnet_group

    def describe_db_subnet_groups(self, subnet_group_name: str) -> List[DBSubnetGroup]:
        if subnet_group_name:
            if subnet_group_name in self.subnet_groups:
                return [self.subnet_groups[subnet_group_name]]
            else:
                raise DBSubnetGroupNotFoundError(subnet_group_name)
        return list(self.subnet_groups.values())

    def modify_db_subnet_group(
        self, subnet_name: str, description: str, subnets: List[Subnet]
    ) -> DBSubnetGroup:
        subnet_group = self.subnet_groups.pop(subnet_name)
        if not subnet_group:
            raise DBSubnetGroupNotFoundError(subnet_name)
        subnet_group.name = subnet_name
        subnet_group.subnets = subnets  # type: ignore[assignment]
        if description is not None:
            subnet_group.description = description
        return subnet_group

    def delete_subnet_group(self, subnet_name: str) -> DBSubnetGroup:
        if subnet_name in self.subnet_groups:
            return self.subnet_groups.pop(subnet_name)
        else:
            raise DBSubnetGroupNotFoundError(subnet_name)

    def create_option_group(self, option_group_kwargs: Dict[str, Any]) -> OptionGroup:
        option_group_id = option_group_kwargs["option_group_name"]
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
            "option_group_description" not in option_group_kwargs
            or not option_group_kwargs["option_group_description"]
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
            option_group_kwargs["major_engine_version"]
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

        option_group = OptionGroup(self, **option_group_kwargs)
        self.option_groups[option_group_id] = option_group
        return option_group

    def delete_option_group(self, option_group_name: str) -> OptionGroup:
        if option_group_name in self.option_groups:
            return self.option_groups.pop(option_group_name)
        else:
            raise OptionGroupNotFoundFaultError(option_group_name)

    def describe_option_groups(
        self, option_group_kwargs: Dict[str, Any]
    ) -> List[OptionGroup]:
        option_group_list = []
        for option_group in self.option_groups.values():
            if (
                option_group_kwargs["option_group_name"]
                and option_group.name != option_group_kwargs["option_group_name"]
            ):
                continue
            elif option_group_kwargs.get(
                "engine_name"
            ) and option_group.engine_name != option_group_kwargs.get("engine_name"):
                continue
            elif option_group_kwargs.get(
                "major_engine_version"
            ) and option_group.major_engine_version != option_group_kwargs.get(
                "major_engine_version"
            ):
                continue
            else:
                option_group_list.append(option_group)
        if not len(option_group_list):
            raise OptionGroupNotFoundFaultError(
                option_group_kwargs["option_group_name"]
            )
        return option_group_list

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
        options_to_remove: Optional[List[str]] = None,
    ) -> OptionGroup:
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
    ) -> DBParameterGroup:
        db_parameter_group_id = db_parameter_group_kwargs["db_parameter_group_name"]
        if db_parameter_group_id in self.db_parameter_groups:
            raise RDSClientError(
                "DBParameterGroupAlreadyExists",
                f"A DB parameter group named {db_parameter_group_id} already exists.",
            )
        if not db_parameter_group_kwargs.get("description"):
            raise RDSClientError(
                "InvalidParameterValue",
                "The parameter Description must be provided and must not be blank.",
            )
        if not db_parameter_group_kwargs.get("db_parameter_group_family"):
            raise RDSClientError(
                "InvalidParameterValue",
                "The parameter DBParameterGroupFamily must be provided and must not be blank.",
            )
        db_parameter_group = DBParameterGroup(self, **db_parameter_group_kwargs)
        self.db_parameter_groups[db_parameter_group_id] = db_parameter_group
        return db_parameter_group

    def copy_db_parameter_group(
        self,
        source_db_parameter_group_identifier: str,
        target_db_parameter_group_identifier: str,
        target_db_parameter_group_description: str,
        tags: Optional[List[Dict[str, str]]] = None,
    ) -> DBParameterGroup:
        if source_db_parameter_group_identifier.startswith("arn:aws:rds:"):
            source_db_parameter_group_identifier = (
                source_db_parameter_group_identifier.split(":")[-1]
            )
        if source_db_parameter_group_identifier not in self.db_parameter_groups:
            raise DBParameterGroupNotFoundError(source_db_parameter_group_identifier)
        if target_db_parameter_group_identifier in self.db_parameter_groups:
            raise DBParameterGroupAlreadyExistsError(
                target_db_parameter_group_identifier
            )
        source_db_parameter_group = self.db_parameter_groups[
            source_db_parameter_group_identifier
        ]

        target_db_parameter_group = DBParameterGroup(
            backend=self,
            db_parameter_group_name=target_db_parameter_group_identifier,
            db_parameter_group_family=source_db_parameter_group.family,
            description=target_db_parameter_group_description,
            tags=tags,
        )
        self.db_parameter_groups[target_db_parameter_group_identifier] = (
            target_db_parameter_group
        )

        iterable_source_parameters = [
            {"ParameterName": name, **values}
            for name, values in source_db_parameter_group.parameters.items()
        ]
        self.modify_db_parameter_group(
            db_parameter_group_name=target_db_parameter_group_identifier,
            db_parameter_group_parameters=iterable_source_parameters,
        )
        return target_db_parameter_group

    def describe_db_parameter_groups(
        self, db_parameter_group_kwargs: Dict[str, Any]
    ) -> List[DBParameterGroup]:
        db_parameter_group_list = []
        for db_parameter_group in self.db_parameter_groups.values():
            if not db_parameter_group_kwargs.get(
                "db_parameter_group_name"
            ) or db_parameter_group.name == db_parameter_group_kwargs.get(
                "db_parameter_group_name"
            ):
                db_parameter_group_list.append(db_parameter_group)
            else:
                continue
        return db_parameter_group_list

    def modify_db_parameter_group(
        self,
        db_parameter_group_name: str,
        db_parameter_group_parameters: Iterable[Dict[str, Any]],
    ) -> DBParameterGroup:
        if db_parameter_group_name not in self.db_parameter_groups:
            raise DBParameterGroupNotFoundError(db_parameter_group_name)

        db_parameter_group = self.db_parameter_groups[db_parameter_group_name]
        db_parameter_group.update_parameters(db_parameter_group_parameters)

        return db_parameter_group

    def modify_db_cluster_parameter_group(
        self,
        db_cluster_parameter_group_name: str,
        db_cluster_parameter_group_parameters: Iterable[Dict[str, Any]],
    ) -> DBClusterParameterGroup:
        if db_cluster_parameter_group_name not in self.db_cluster_parameter_groups:
            raise DBClusterParameterGroupNotFoundError(db_cluster_parameter_group_name)

        db_cluster_parameter_group = self.db_cluster_parameter_groups[
            db_cluster_parameter_group_name
        ]
        db_cluster_parameter_group.update_cluster_parameters(
            db_cluster_parameter_group_parameters
        )

        return db_cluster_parameter_group

    def describe_db_cluster_parameters(
        self, db_cluster_parameter_group_name: str
    ) -> List[Dict[str, Any]]:
        if db_cluster_parameter_group_name not in self.db_cluster_parameter_groups:
            raise DBClusterParameterGroupNotFoundError(db_cluster_parameter_group_name)
        db_cluster_parameter_group = self.db_cluster_parameter_groups[
            db_cluster_parameter_group_name
        ]
        parameters = [
            {"ParameterName": name, **values}
            for name, values in db_cluster_parameter_group.parameters.items()
        ]
        return parameters

    def create_db_cluster(self, kwargs: Dict[str, Any]) -> DBCluster:
        cluster_id = kwargs["db_cluster_identifier"]
        cluster = DBCluster(self, **kwargs)
        self.clusters[cluster_id] = cluster
        cluster.save_automated_backup()
        if cluster.global_cluster_identifier:
            for regional_backend in rds_backends[self.account_id]:
                if (
                    cluster.global_cluster_identifier
                    in rds_backends[self.account_id][regional_backend].global_clusters
                ):
                    global_cluster = rds_backends[self.account_id][
                        regional_backend
                    ].global_clusters[cluster.global_cluster_identifier]
                    global_cluster.members.append(cluster)
                    if len(global_cluster.members) == 1:
                        # primary cluster
                        setattr(cluster, "is_writer", True)
                    else:
                        # secondary cluster(s)
                        setattr(cluster, "is_writer", False)

        if cluster.replication_source_identifier:
            cluster_identifier = cluster.replication_source_identifier
            original_cluster = find_cluster(cluster_identifier)
            original_cluster.read_replica_identifiers.append(cluster.db_cluster_arn)

        return cluster

    def modify_db_cluster(self, kwargs: Dict[str, Any]) -> DBCluster:
        cluster_id = kwargs["db_cluster_identifier"]

        cluster = self.clusters[cluster_id]
        del self.clusters[cluster_id]

        if kwargs.get("rotate_master_user_password") and kwargs.get(
            "apply_immediately"
        ):
            cluster.master_user_secret.rotate_secret()
        elif kwargs.get("rotate_master_user_password") and not kwargs.get(
            "apply_immediately"
        ):
            raise RDSClientError(
                "InvalidParameterCombination",
                "You must specify apply immediately when rotating the master user password.",
            )

        if "manage_master_user_password" in kwargs:
            manage_master_user_password = kwargs.pop("manage_master_user_password")
            if manage_master_user_password:
                kms_key_id = kwargs.pop("master_user_secret_kms_key_id", None)
                cluster.master_user_secret = MasterUserSecret(cluster, kms_key_id)
                cluster.manage_master_user_password = True
            elif not manage_master_user_password and hasattr(
                cluster, "master_user_secret"
            ):
                cluster.master_user_secret.delete_secret()
                del cluster.master_user_secret
                cluster.manage_master_user_password = False

        kwargs["db_cluster_identifier"] = kwargs.get("new_db_cluster_identifier", None)
        for k, v in kwargs.items():
            if k == "db_cluster_parameter_group_name":
                k = "parameter_group"
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

        initial_state = copy.copy(cluster)  # Return status=creating

        # Already set the final status in the background
        cluster.status = "available"

        return initial_state

    def promote_read_replica_db_cluster(self, db_cluster_identifier: str) -> DBCluster:
        cluster = self.clusters[db_cluster_identifier]
        source_cluster = find_cluster(cluster.replication_source_identifier)  # type: ignore
        source_cluster.read_replica_identifiers.remove(cluster.db_cluster_arn)
        cluster.replication_source_identifier = None
        return cluster

    def create_auto_cluster_snapshot(
        self, db_cluster_identifier: str, db_snapshot_identifier: str
    ) -> DBClusterSnapshot:
        return self.create_db_cluster_snapshot(
            db_cluster_identifier, db_snapshot_identifier, snapshot_type="automated"
        )

    def create_db_cluster_snapshot(
        self,
        db_cluster_identifier: str,
        db_cluster_snapshot_identifier: str,
        snapshot_type: str = "manual",
        tags: Optional[List[Dict[str, str]]] = None,
    ) -> DBClusterSnapshot:
        if db_cluster_snapshot_identifier in self.cluster_snapshots:
            raise DBClusterSnapshotAlreadyExistsError(db_cluster_snapshot_identifier)
        if db_cluster_identifier not in self.clusters:
            raise DBClusterNotFoundError(db_cluster_identifier)
        if snapshot_type == "manual" and len(self.cluster_snapshots) >= int(
            os.environ.get("MOTO_RDS_SNAPSHOT_LIMIT", "100")
        ):
            raise SnapshotQuotaExceededFault()
        cluster = self.clusters[db_cluster_identifier]
        snapshot = DBClusterSnapshot(
            self, cluster, db_cluster_snapshot_identifier, snapshot_type, tags
        )
        self.cluster_snapshots[db_cluster_snapshot_identifier] = snapshot
        return snapshot

    def copy_db_cluster_snapshot(
        self,
        source_db_cluster_snapshot_identifier: str,
        target_db_cluster_snapshot_identifier: str,
        kms_key_id: Optional[str] = None,
        copy_tags: bool = False,
        tags: Optional[List[Dict[str, str]]] = None,
        **_: Any,
    ) -> DBClusterSnapshot:
        if target_db_cluster_snapshot_identifier in self.cluster_snapshots:
            raise DBClusterSnapshotAlreadyExistsError(
                target_db_cluster_snapshot_identifier
            )

        if len(self.cluster_snapshots) >= int(
            os.environ.get("MOTO_RDS_SNAPSHOT_LIMIT", "100")
        ):
            raise SnapshotQuotaExceededFault()
        if kms_key_id is not None:
            kms_key_id = self._validate_kms_key(kms_key_id)
        source_db_cluster_snapshot = self.get_db_cluster_snapshot(
            source_db_cluster_snapshot_identifier
        )
        # If tags are supplied, copy_tags is ignored (verified against real AWS backend 2025-02-12).
        if copy_tags and not tags:
            tags = source_db_cluster_snapshot.tags or []
        target_db_cluster_snapshot = DBClusterSnapshot(
            self,
            source_db_cluster_snapshot.cluster,
            target_db_cluster_snapshot_identifier,
            tags=tags,
            kms_key_id=kms_key_id,
        )
        self.cluster_snapshots[target_db_cluster_snapshot_identifier] = (
            target_db_cluster_snapshot
        )
        return target_db_cluster_snapshot

    def delete_db_cluster_snapshot(
        self, db_snapshot_identifier: str
    ) -> DBClusterSnapshot:
        if db_snapshot_identifier not in self.cluster_snapshots:
            raise DBClusterSnapshotNotFoundError(db_snapshot_identifier)

        return self.cluster_snapshots.pop(db_snapshot_identifier)

    def describe_db_clusters(
        self, db_cluster_identifier: Optional[str] = None, filters: Any = None
    ) -> List[DBCluster]:
        clusters = self.clusters
        if db_cluster_identifier:
            filters = merge_filters(filters, {"db-cluster-id": [db_cluster_identifier]})
        if filters:
            clusters = self._filter_resources(clusters, filters, DBCluster)
        if db_cluster_identifier and not clusters:
            raise DBClusterNotFoundError(db_cluster_identifier)
        return list(clusters.values())

    def describe_db_cluster_snapshots(
        self,
        db_cluster_identifier: Optional[str],
        db_snapshot_identifier: str,
        snapshot_type: Optional[str] = None,
        filters: Any = None,
    ) -> List[DBClusterSnapshot]:
        if snapshot_type == "shared":
            return self.get_shared_db_cluster_snapshots()  # type: ignore[return-value]
        snapshots = self.cluster_snapshots
        if db_cluster_identifier:
            filters = merge_filters(filters, {"db-cluster-id": [db_cluster_identifier]})
        if db_snapshot_identifier:
            filters = merge_filters(
                filters, {"db-cluster-snapshot-id": [db_snapshot_identifier]}
            )
        snapshot_types = (
            ["automated", "manual"]
            if (
                snapshot_type is None
                and (filters is not None and "snapshot-type" not in filters)
            )
            else [snapshot_type]
            if snapshot_type is not None
            else []
        )
        if snapshot_types:
            filters = merge_filters(filters, {"snapshot-type": snapshot_types})
        if filters:
            snapshots = self._filter_resources(snapshots, filters, DBClusterSnapshot)
        if db_snapshot_identifier and not snapshots and not db_cluster_identifier:
            raise DBClusterSnapshotNotFoundError(db_snapshot_identifier)
        return list(snapshots.values())

    def delete_db_cluster(
        self, cluster_identifier: str, snapshot_name: Optional[str] = None
    ) -> DBCluster:
        if cluster_identifier in self.clusters:
            cluster = self.clusters[cluster_identifier]
            if cluster.deletion_protection:
                raise InvalidParameterCombination(
                    "Cannot delete protected Cluster, please disable deletion protection and try again."
                )
            if cluster.members:
                raise DBClusterToBeDeletedHasActiveMembers()
            global_id = cluster.global_cluster_identifier or ""
            if global_id in self.global_clusters:
                self.remove_from_global_cluster(global_id, cluster_identifier)
            if cluster.manage_master_user_password:
                cluster.master_user_secret.delete_secret()
            if snapshot_name:
                self.create_db_cluster_snapshot(cluster_identifier, snapshot_name)
            return self.clusters.pop(cluster_identifier)
        raise DBClusterNotFoundError(cluster_identifier)

    def start_db_cluster(self, cluster_identifier: str) -> DBCluster:
        if cluster_identifier not in self.clusters:
            raise DBClusterNotFoundError(cluster_identifier)
        cluster = self.clusters[cluster_identifier]
        if cluster.status != "stopped":
            raise InvalidDBClusterStateFault(
                f"DbCluster {cluster_identifier} is in {cluster.status} state "
                "but expected it to be one of "
                "stopped,inaccessible-encryption-credentials-recoverable."
            )
        temp_state = copy.deepcopy(cluster)
        temp_state.status = "started"
        cluster.status = "available"  # This is the final status - already setting it in the background
        return temp_state

    def restore_db_cluster_from_snapshot(
        self, from_snapshot_id: str, overrides: Dict[str, Any]
    ) -> DBCluster:
        snapshot = self.describe_db_cluster_snapshots(
            db_cluster_identifier=None, db_snapshot_identifier=from_snapshot_id
        )[0]
        original_cluster = snapshot.cluster
        new_cluster_props = copy.copy(original_cluster.get_cfg())
        for key, value in overrides.items():
            if value:
                new_cluster_props[key] = value

        return self.create_db_cluster(new_cluster_props)

    def stop_db_cluster(self, cluster_identifier: str) -> DBCluster:
        if cluster_identifier not in self.clusters:
            raise DBClusterNotFoundError(cluster_identifier)
        cluster = self.clusters[cluster_identifier]
        if cluster.status not in ["available"]:
            raise InvalidDBClusterStateFault(
                f"DbCluster {cluster_identifier} is in {cluster.status} state "
                "but expected it to be one of available."
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
            raise DBSnapshotNotFoundFault(snapshot_id)
        elif (
            snapshot_type == "cluster-snapshot"
            and snapshot_id not in self.cluster_snapshots
        ):
            raise DBClusterSnapshotNotFoundError(snapshot_id)

        if snapshot_type == "snapshot":
            snapshot: Union[DBSnapshot, DBClusterSnapshot] = self.database_snapshots[
                snapshot_id
            ]
        else:
            snapshot = self.cluster_snapshots[snapshot_id]

        if snapshot.status not in ["available"]:
            raise InvalidExportSourceStateError(snapshot.status)

        export_task = ExportTask(self, snapshot, kwargs)
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

        subscription = EventSubscription(self, **kwargs)
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

    def _find_resource(self, resource_type: str, resource_name: str) -> Any:
        for resource_class, resources in self.resource_map.items():
            if resource_type == getattr(resource_class, "resource_type", ""):
                if resource_name in resources:
                    return resources[resource_name]
        # The resource_name is the last part of the ARN
        # Usually that's the name - but for DBProxies, the last part of the ARN is a random identifier
        # So we can't just use the dict-keys - we have to manually check the ARN
        if resource_type == "db-proxy":
            for resource in self.db_proxies.values():
                if resource.arn.endswith(resource_name):
                    return resource

    def _get_resource_for_tagging(self, arn: str) -> Any:
        if self.arn_regex.match(arn):
            arn_breakdown = arn.split(":")
            resource_type = arn_breakdown[len(arn_breakdown) - 2]
            resource_name = arn_breakdown[len(arn_breakdown) - 1]
            # FIXME: HACK for automated snapshots
            if resource_type == "rds":
                resource_type = arn_breakdown[-3]
                resource_name = arn_breakdown[-2] + ":" + arn_breakdown[-1]
            resource = self._find_resource(resource_type, resource_name)
            return resource
        raise RDSClientError("InvalidParameterValue", f"Invalid resource name: {arn}")

    def list_tags_for_resource(self, arn: str) -> List[Dict[str, str]]:
        resource = self._get_resource_for_tagging(arn)
        if resource:
            return resource.get_tags()
        return []

    def remove_tags_from_resource(self, arn: str, tag_keys: List[str]) -> None:
        resource = self._get_resource_for_tagging(arn)
        if resource:
            resource.remove_tags(tag_keys)

    def add_tags_to_resource(
        self, arn: str, tags: List[Dict[str, str]]
    ) -> List[Dict[str, str]]:
        resource = self._get_resource_for_tagging(arn)
        if resource:
            return resource.add_tags(tags)
        return []

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
    def _validate_db_identifier(db_identifier: str) -> None:
        # https://docs.aws.amazon.com/AmazonRDS/latest/APIReference/API_CreateDBInstance.html
        # Constraints:
        # # Must contain from 1 to 63 letters, numbers, or hyphens.
        # # First character must be a letter.
        # # Can't end with a hyphen or contain two consecutive hyphens.
        if (
            re.match(
                "^(?!.*--)([a-zA-Z]?[a-zA-Z0-9-]{0,61}[a-zA-Z0-9])$", db_identifier
            )
            and db_identifier[0].isalpha()
        ):
            return
        raise InvalidDBInstanceIdentifier

    @staticmethod
    def validate_db_snapshot_identifier(
        db_snapshot_identifier: str, parameter_name: str
    ) -> None:
        # https://docs.aws.amazon.com/AmazonRDS/latest/APIReference/API_CreateDBSnapshot.html
        # Constraints:
        # # Must contain from 1 to 255 letters, numbers, or hyphens.
        # # First character must be a letter.
        # # Can't end with a hyphen or contain two consecutive hyphens.
        if (
            re.match(
                "^(?!.*--)([a-zA-Z]?[a-zA-Z0-9-]{0,253}[a-zA-Z0-9])$",
                db_snapshot_identifier,
            )
            and db_snapshot_identifier[0].isalpha()
        ):
            return
        raise InvalidDBSnapshotIdentifier(db_snapshot_identifier, parameter_name)

    def describe_orderable_db_instance_options(
        self, engine: str, engine_version: str
    ) -> List[Dict[str, Any]]:
        """
        Only the Aurora-Postgresql and Neptune-engine is currently implemented
        """
        if engine in ["aurora-postgresql", "neptune"]:
            if engine_version:
                return [
                    option
                    for option in self.db_cluster_options(engine)
                    if option["EngineVersion"] == engine_version
                ]
            return self.db_cluster_options(engine)
        return []

    def create_db_cluster_parameter_group(
        self,
        group_name: str,
        family: str,
        description: str,
    ) -> DBClusterParameterGroup:
        group = DBClusterParameterGroup(
            backend=self,
            name=group_name,
            family=family,
            description=description,
        )
        self.db_cluster_parameter_groups[group_name] = group
        return group

    def describe_db_cluster_parameter_groups(
        self, group_name: str
    ) -> List[DBClusterParameterGroup]:
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
        storage_encrypted: Optional[bool],
        deletion_protection: Optional[bool],
    ) -> GlobalCluster:
        source_cluster = None
        if source_db_cluster_identifier is not None:
            # validate our source cluster exists
            if not re.match(ARN_PARTITION_REGEX + ":rds", source_db_cluster_identifier):
                raise InvalidParameterValue("Malformed db cluster arn dbci")
            source_cluster = self.describe_db_clusters(
                db_cluster_identifier=source_db_cluster_identifier
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
            backend=self,
            global_cluster_identifier=global_cluster_identifier,
            engine=engine,
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
        return list(self.global_clusters.values())

    def delete_global_cluster(self, global_cluster_identifier: str) -> GlobalCluster:
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
                db_cluster_identifier=db_cluster_identifier
            )[0]
            global_cluster.members.remove(cluster)
            return global_cluster
        except:  # noqa: E722 Do not use bare except
            pass
        return None

    def describe_db_snapshot_attributes(
        self, db_snapshot_identifier: str
    ) -> Dict[str, List[str]]:
        snapshot = self.describe_db_snapshots(
            db_instance_identifier=None, db_snapshot_identifier=db_snapshot_identifier
        )[0]
        return snapshot.attributes

    def modify_db_snapshot_attribute(
        self,
        db_snapshot_identifier: str,
        attribute_name: str,
        values_to_add: Optional[List[str]] = None,
        values_to_remove: Optional[List[str]] = None,
    ) -> Dict[str, List[str]]:
        snapshot = self.describe_db_snapshots(
            db_instance_identifier=None, db_snapshot_identifier=db_snapshot_identifier
        )[0]
        if db_snapshot_identifier.startswith("rds:"):  # automated snapshot
            raise InvalidParameterValue(
                "The parameter DBSnapshotIdentifier is not a valid identifier."
            )
        snapshot.modify_attribute(attribute_name, values_to_add, values_to_remove)
        return snapshot.attributes

    def describe_db_cluster_snapshot_attributes(
        self, db_cluster_snapshot_identifier: str
    ) -> Dict[str, List[str]]:
        snapshot = self.describe_db_cluster_snapshots(
            db_cluster_identifier=None,
            db_snapshot_identifier=db_cluster_snapshot_identifier,
        )[0]
        return snapshot.attributes

    def modify_db_cluster_snapshot_attribute(
        self,
        db_cluster_snapshot_identifier: str,
        attribute_name: str,
        values_to_add: Optional[List[str]] = None,
        values_to_remove: Optional[List[str]] = None,
    ) -> Dict[str, List[str]]:
        snapshot = self.describe_db_cluster_snapshots(
            db_cluster_identifier=None,
            db_snapshot_identifier=db_cluster_snapshot_identifier,
        )[0]
        if snapshot.snapshot_type != "manual":
            raise InvalidDBClusterSnapshotStateFault(
                "automated snapshots cannot be modified."
            )
        snapshot.modify_attribute(attribute_name, values_to_add, values_to_remove)
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
            self,
            db_proxy_name,
            engine_family,
            auth,
            role_arn,
            vpc_subnet_ids,
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

    def deregister_db_proxy_targets(
        self,
        db_proxy_name: str,
        target_group_name: str,
        db_cluster_identifiers: List[str],
        db_instance_identifiers: List[str],
    ) -> None:
        db_proxy = self.db_proxies[db_proxy_name]
        target_group = db_proxy.proxy_target_groups[target_group_name or "default"]
        target_group.targets = [
            t
            for t in target_group.targets
            if t.rds_resource_id not in db_cluster_identifiers
            and t.rds_resource_id not in db_instance_identifiers
        ]

    def register_db_proxy_targets(
        self,
        db_proxy_name: str,
        target_group_name: str,
        db_cluster_identifiers: List[str],
        db_instance_identifiers: List[str],
    ) -> List[DBProxyTarget]:
        db_proxy = self.db_proxies[db_proxy_name]
        target_group = db_proxy.proxy_target_groups[target_group_name or "default"]
        new_targets = []
        for cluster_id in db_cluster_identifiers:
            cluster = self.clusters[cluster_id]
            target = DBProxyTarget(
                backend=self,
                resource_id=cluster_id,
                endpoint=cluster.endpoint,
                type="TRACKED_CLUSTER",
            )
            new_targets.append(target)
        for instance_id in db_instance_identifiers:
            target = DBProxyTarget(
                backend=self,
                resource_id=instance_id,
                endpoint=None,
                type="RDS_INSTANCE",
            )
            new_targets.append(target)
        target_group.targets.extend(new_targets)
        return new_targets

    def delete_db_proxy(self, proxy_name: str) -> DBProxy:
        return self.db_proxies.pop(proxy_name)

    def describe_db_proxy_targets(self, proxy_name: str) -> List[DBProxyTarget]:
        proxy = self.db_proxies[proxy_name]
        target_group = proxy.proxy_target_groups["default"]
        return target_group.targets

    def describe_db_proxy_target_groups(
        self, proxy_name: str
    ) -> List[DBProxyTargetGroup]:
        proxy = self.db_proxies[proxy_name]
        return list(proxy.proxy_target_groups.values())

    def modify_db_proxy_target_group(
        self, proxy_name: str, config: Dict[str, Any]
    ) -> DBProxyTargetGroup:
        proxy = self.db_proxies[proxy_name]
        target_group = proxy.proxy_target_groups["default"]
        if max_connections := config.get("MaxConnectionsPercent"):
            target_group.max_connections = max_connections
        if max_idle := config.get("MaxIdleConnectionsPercent"):
            target_group.max_idle_connections = max_idle
        else:
            target_group.max_idle_connections = math.floor(
                int(target_group.max_connections) / 2
            )
        target_group.borrow_timeout = config.get(
            "ConnectionBorrowTimeout", target_group.borrow_timeout
        )
        if "SessionPinningFilters" in config:
            target_group.session_pinning_filters = config["SessionPinningFilters"]
        return target_group

    def describe_db_instance_automated_backups(
        self,
        db_instance_identifier: Optional[str] = None,
        **_: Any,
    ) -> List[DBInstanceAutomatedBackup]:
        snapshots = list(self.database_snapshots.values())
        if db_instance_identifier is not None:
            snapshots = [
                snap
                for snap in self.database_snapshots.values()
                if snap.db_instance_identifier == db_instance_identifier
            ]
        snapshots_grouped = defaultdict(list)
        for snapshot in snapshots:
            if snapshot.snapshot_type == "automated":
                snapshots_grouped[snapshot.db_instance_identifier].append(snapshot)
        return [
            DBInstanceAutomatedBackup(self, k, v) for k, v in snapshots_grouped.items()
        ]

    def add_event(self, event_type: str, resource: ResourceWithEvents) -> None:
        event = Event(event_type, resource)
        self.events.append(event)

    def describe_events(
        self,
        source_identifier: Optional[str] = None,
        source_type: Optional[str] = None,
        **_: Any,
    ) -> List[Event]:
        if source_identifier is not None and source_type is None:
            raise InvalidParameterCombination(
                "Cannot specify source identifier without source type"
            )
        events = self.events
        if source_identifier is not None:
            events = [e for e in events if e.source_identifier == source_identifier]
        if source_type is not None:
            events = [e for e in events if e.source_type == source_type]
        return events

    def describe_db_log_files(self, db_instance_identifier: str) -> List[DBLogFile]:
        if db_instance_identifier not in self.databases:
            raise DBInstanceNotFoundError(db_instance_identifier)
        database = self.databases[db_instance_identifier]
        return database.log_file_manager.files

    def create_blue_green_deployment(
        self,
        bg_kwargs: Dict[str, Any],
    ) -> BlueGreenDeployment:
        bg_name = bg_kwargs.get("blue_green_deployment_name", "")

        for existing_bg in self.blue_green_deployments.values():
            if existing_bg.blue_green_deployment_name == bg_name:
                raise BlueGreenDeploymentAlreadyExistsFault(bg_name)
        bg_deployment: BlueGreenDeployment = BlueGreenDeployment(self, **bg_kwargs)
        self.blue_green_deployments[bg_deployment.blue_green_deployment_identifier] = (
            bg_deployment
        )
        # update states only for the response of create command
        bg_copy = copy.copy(bg_deployment)
        bg_copy.set_status(bg_deployment.SupportedStates.PROVISIONING)

        return bg_copy

    def describe_blue_green_deployments(
        self,
        blue_green_deployment_identifier: Optional[str] = None,
        filters: Any = None,
    ) -> List[BlueGreenDeployment]:
        bg_deployments = self.blue_green_deployments
        if blue_green_deployment_identifier:
            filters = merge_filters(
                filters,
                {
                    "blue-green-deployment-identifier": [
                        blue_green_deployment_identifier
                    ]
                },
            )
        if filters:
            bg_deployments = self._filter_resources(
                bg_deployments, filters, BlueGreenDeployment
            )
        if blue_green_deployment_identifier and not bg_deployments:
            raise BlueGreenDeploymentNotFoundFault(blue_green_deployment_identifier)
        return list(bg_deployments.values())

    def switchover_blue_green_deployment(
        self,
        blue_green_deployment_identifier: str,
        switchover_timeout: Optional[int] = 300,
    ) -> BlueGreenDeployment:
        if blue_green_deployment_identifier not in self.blue_green_deployments:
            raise BlueGreenDeploymentNotFoundFault(blue_green_deployment_identifier)

        bg_deployment = self.blue_green_deployments[blue_green_deployment_identifier]

        if bg_deployment.status != "AVAILABLE":
            raise InvalidBlueGreenDeploymentStateFault(blue_green_deployment_identifier)

        bg_deployment_before_switch = copy.copy(bg_deployment)

        bg_deployment.switchover()
        self.blue_green_deployments[blue_green_deployment_identifier] = bg_deployment

        bg_deployment_before_switch.set_status(
            bg_deployment.SupportedStates.SWITCHOVER_IN_PROGRESS
        )

        return bg_deployment_before_switch

    def delete_blue_green_deployment(
        self, blue_green_deployment_identifier: str, delete_target: bool = False
    ) -> BlueGreenDeployment:
        if blue_green_deployment_identifier not in self.blue_green_deployments:
            raise BlueGreenDeploymentNotFoundFault(blue_green_deployment_identifier)
        bg_deployment = self.blue_green_deployments.pop(
            blue_green_deployment_identifier
        )

        if (
            delete_target
            and bg_deployment.status
            != bg_deployment.SupportedStates.SWITCHOVER_COMPLETED.name
        ):
            if self._is_cluster(bg_deployment.target):
                cluster = find_cluster(bg_deployment.target)
                for member in cluster.members:
                    self.delete_db_instance(
                        db_instance_identifier=member.db_instance_identifier
                    )
                self.delete_db_cluster(cluster_identifier=cluster.db_cluster_identifier)
            else:
                instance = self.find_db_from_id(bg_deployment.target)
                self.delete_db_instance(
                    db_instance_identifier=instance.db_instance_identifier
                )

        bg_deployment.set_status(bg_deployment.SupportedStates.DELETING)
        bg_deployment.deletion_time = utcnow()
        return bg_deployment

    def create_db_shard_group(self, kwargs: Dict[str, Any]) -> DBShardGroup:
        db_shard_group = DBShardGroup(self, **kwargs)

        # Validate
        db_shard_group_identifier = kwargs["db_shard_group_identifier"]
        db_cluster_identifier = kwargs["db_cluster_identifier"]
        compute_redundancy = kwargs.get("compute_redundancy")
        max_acu = kwargs["max_acu"]
        min_acu = kwargs.get("min_acu")

        if db_shard_group_identifier in self.shard_groups:
            raise DBShardGroupAlreadyExistsError(db_shard_group_identifier)

        if db_cluster_identifier not in self.clusters:
            raise DBClusterNotFoundError(db_cluster_identifier)

        if compute_redundancy not in (None, 0, 1, 2):
            raise InvalidParameterValue(
                f"Invalid ComputeRedundancy value: '{compute_redundancy}'. "
                "Valid values are 0 (no standby), 1 (1 standby AZ), 2 (2 standby AZs)."
            )

        if "min_acu" in kwargs and min_acu >= max_acu:
            raise InvalidParameterValue("min_acu cannot be larger than max_acu")

        self.shard_groups[db_shard_group.db_shard_group_identifier] = db_shard_group
        return db_shard_group

    def describe_db_shard_groups(
        self,
        db_shard_group_identifier: Optional[str],
        filters: Optional[Dict[str, List[str]]],
    ) -> List[DBShardGroup]:
        shard_groups = self.shard_groups
        if db_shard_group_identifier:
            filters = merge_filters(
                filters, {"db-shard-group-id": [db_shard_group_identifier]}
            )
        if filters:
            shard_groups = self._filter_resources(shard_groups, filters, DBShardGroup)
        if db_shard_group_identifier and not shard_groups:
            raise DBShardGroupNotFoundFault(db_shard_group_identifier)
        return list(shard_groups.values())

    def _is_cluster(self, arn: str) -> bool:
        return arn.split(":")[-2] == "cluster"


class OptionGroup(RDSBaseModel):
    resource_type = "og"

    def __init__(
        self,
        backend: RDSBackend,
        option_group_name: str,
        engine_name: str,
        major_engine_version: str,
        option_group_description: Optional[str] = None,
        tags: Optional[List[Dict[str, str]]] = None,
    ):
        super().__init__(backend)
        self.engine_name = engine_name
        self.major_engine_version = major_engine_version
        self.description = option_group_description
        self.name = option_group_name
        self.vpc_and_non_vpc_instance_memberships = False
        self._options: Dict[str, Any] = {}
        self.vpcId = "null"
        self.tags = tags or []

    @property
    def name(self) -> str:
        return self._name

    @name.setter
    def name(self, value: str) -> None:
        self._name = value.lower()

    @property
    def resource_id(self) -> str:
        return self.name

    @property
    def options(self) -> List[Dict[str, Any]]:  # type: ignore[misc]
        return [
            {
                "OptionName": name,
                "OptionSettings": [
                    {
                        "Name": setting.get("Name"),
                        "Value": setting.get("Value"),
                    }
                    for setting in option_settings
                ],
            }
            for name, option_settings in self._options.items()
        ]

    def remove_options(self, options_to_remove: List[str]) -> None:
        for option in options_to_remove:
            self._options.pop(option, None)

    def add_options(self, options_to_add: List[Dict[str, Any]]) -> None:
        for option in options_to_add:
            self._options[option["OptionName"]] = option.get("OptionSettings", {})


class DBParameterGroup(CloudFormationModel, RDSBaseModel):
    resource_type = "pg"

    def __init__(
        self,
        backend: RDSBackend,
        db_parameter_group_name: str,
        description: str,
        db_parameter_group_family: Optional[str],
        tags: Optional[List[Dict[str, str]]] = None,
    ):
        super().__init__(backend)
        self.name = db_parameter_group_name
        self.description = description
        self.family: str | None = db_parameter_group_family
        self.tags = tags or []
        self.parameters: Dict[str, Any] = defaultdict(dict)

    @property
    def name(self) -> str:
        return self._name

    @name.setter
    def name(self, value: str) -> None:
        self._name = value.lower()

    @property
    def resource_id(self) -> str:
        return self.name

    def update_parameters(self, new_parameters: Iterable[Dict[str, Any]]) -> None:
        for new_parameter in new_parameters:
            parameter = self.parameters[new_parameter["ParameterName"]]
            parameter.update(new_parameter)

    def delete(self, account_id: str, region_name: str) -> None:
        backend = rds_backends[account_id][region_name]
        backend.delete_db_parameter_group(self.name)

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
    ) -> DBParameterGroup:
        properties = cloudformation_json["Properties"]

        db_parameter_group_kwargs = {
            "description": properties["Description"],
            "db_parameter_group_family": properties["Family"],
            "db_parameter_group_name": resource_name.lower(),
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


class DBClusterParameterGroup(CloudFormationModel, RDSBaseModel):
    resource_type = "cluster-pg"

    def __init__(self, backend: RDSBackend, name: str, description: str, family: str):
        super().__init__(backend)
        self.name = name
        self.description = description
        self.db_parameter_group_family = family
        self.parameters: Dict[str, Any] = defaultdict(dict)

    @property
    def name(self) -> str:
        return self._name

    @name.setter
    def name(self, value: str) -> None:
        self._name = value.lower()

    @property
    def resource_id(self) -> str:
        return self.name

    def update_cluster_parameters(
        self, new_parameters: Iterable[Dict[str, Any]]
    ) -> None:
        for new_parameter in new_parameters:
            parameter = self.parameters[new_parameter["ParameterName"]]
            parameter.update(new_parameter)


class Event:
    EVENT_MAP = {
        "DB_INSTANCE_BACKUP_START": {
            "Categories": ["backup"],
            "Message": "Backing up DB instance",
        },
        "DB_INSTANCE_BACKUP_FINISH": {
            "Categories": ["backup"],
            "Message": "Finished DB instance backup",
        },
        "DB_INSTANCE_CREATE": {
            "Categories": ["creation"],
            "Message": "DB instance created",
        },
        "DB_INSTANCE_RESET_MASTER_CREDENTIALS": {
            "Categories": ["configuration change"],
            "Message": "Reset master credentials",
        },
        "DB_SNAPSHOT_CREATE_AUTOMATED_START": {
            "Categories": ["creation"],
            "Message": "Creating automated snapshot",
        },
        "DB_SNAPSHOT_CREATE_AUTOMATED_FINISH": {
            "Categories": ["creation"],
            "Message": "Automated snapshot created",
        },
        "DB_SNAPSHOT_CREATE_MANUAL_START": {
            "Categories": ["creation"],
            "Message": "Creating manual snapshot",
        },
        "DB_SNAPSHOT_CREATE_MANUAL_FINISH": {
            "Categories": ["creation"],
            "Message": "Manual snapshot created",
        },
    }

    def __init__(self, event_type: str, resource: ResourceWithEvents) -> None:
        event_metadata = self.EVENT_MAP[event_type]
        self.source_identifier = resource.resource_id
        self.source_type = resource.event_source_type
        self.message = event_metadata["Message"]
        self.event_categories = event_metadata["Categories"]
        self.source_arn = resource.arn
        self.date = utcnow()


class BlueGreenDeployment(RDSBaseModel):
    class SupportedStates(Enum):
        PROVISIONING = 1
        AVAILABLE = 2
        SWITCHOVER_IN_PROGRESS = 3
        SWITCHOVER_COMPLETED = 4
        INVALID_CONFIGURATION = 5
        SWITCHOVER_FAILED = 6
        DELETING = 7

    SUPPORTED_FILTERS = {
        "blue-green-deployment-identifier": FilterDef(
            ["blue_green_deployment_identifier"],
            "System-generated BlueGreen Deployment Identifier",
        ),
        "blue-green-deployment-name": FilterDef(
            ["blue_green_deployment_name"],
            "User-provided Name for BlueGreen Deployment",
        ),
        "source": FilterDef(["source"], "Source Database Instance or Cluster"),
        "target": FilterDef(["target"], "Target Database Instance or Cluster"),
    }

    def __init__(
        self,
        backend: RDSBackend,
        blue_green_deployment_name: str,
        source: str,
        tags: List[Dict[str, str]] | None = None,
        **kwargs: Any,
    ) -> None:
        super().__init__(backend)
        self.blue_green_deployment_identifier = self._generate_id()
        self.blue_green_deployment_name = blue_green_deployment_name
        self.source = source
        self.target = self._create_green_instance(
            target_engine_version=kwargs.get("target_engine_version", None),
            target_db_parameter_group_name=kwargs.get(
                "target_db_parameter_group_name", None
            ),
            target_db_cluster_parameter_group_name=kwargs.get(
                "target_db_cluster_parameter_group_name", None
            ),
            target_db_instance_class=kwargs.get("target_db_instance_class", None),
            target_iops=kwargs.get("target_iops", None),
            target_storage_type=kwargs.get("target_storage_type", None),
            target_allocated_storage=kwargs.get("target_allocated_storage", None),
            target_storage_throughput=kwargs.get("target_storage_throughput", None),
        )
        self.tasks = [
            BlueGreenDeploymentTask(
                "CREATING_READ_REPLICA_OF_SOURCE",
                BlueGreenDeploymentTask.SupportedStates.COMPLETED.name,
            ),
            BlueGreenDeploymentTask(
                "CONFIGURE_BACKUPS",
                BlueGreenDeploymentTask.SupportedStates.COMPLETED.name,
            ),
        ]
        self.set_status(self.SupportedStates.AVAILABLE)
        self.status_details: str | None = None
        self.create_time = utcnow()
        self.deletion_time: datetime | None = None
        self.tags = tags or []

    @property
    def resource_id(self) -> str:
        return self.blue_green_deployment_identifier

    def set_status(self, value: SupportedStates) -> None:
        self.status = value.name
        self.switchover_details = [
            SwitchoverDetails(self.source, self.target, self.status)
        ]

    def _generate_id(self) -> str:
        return "bgd-" + "".join(
            random.choices(string.ascii_lowercase + string.digits, k=16)
        )

    def _generate_green_id(self, db_identifier: str) -> str:
        return f"{db_identifier}-green-" + "".join(
            random.choices(string.ascii_lowercase, k=6)
        )

    def _create_green_instance(
        self,
        target_engine_version: str | None,
        target_db_parameter_group_name: str | None,
        target_db_cluster_parameter_group_name: str | None,
        target_db_instance_class: str | None,
        target_iops: int | None,
        target_storage_type: str | None,
        target_allocated_storage: int | None,
        target_storage_throughput: int | None,
    ) -> str:
        source_instance: DBInstance | DBCluster = self._get_valid_instance_or_raise(
            self.source
        )
        green_instance: DBInstance | DBCluster

        db_kwargs = {
            "engine": source_instance.engine,
            "engine_version": target_engine_version or source_instance.engine_version,
            "iops": target_iops or source_instance.iops,
            "port": source_instance.port,
            "backup_retention_period": source_instance.backup_retention_period,
            "character_set_name": source_instance.character_set_name,
            "auto_minor_version_upgrade": source_instance.auto_minor_version_upgrade,
            "copy_tags_to_snapshot": source_instance.copy_tags_to_snapshot,
            "master_username": source_instance.master_username,
            "multi_az": source_instance.multi_az,
            "license_model": source_instance.license_model,
            "preferred_backup_window": source_instance.preferred_backup_window,
            "preferred_maintenance_window": source_instance.preferred_maintenance_window,
            "storage_encrypted": source_instance.storage_encrypted,
            "tags": source_instance.tags,
            "deletion_protection": source_instance.deletion_protection,
            "enable_cloudwatch_logs_exports": source_instance.enabled_cloudwatch_logs_exports,
            "master_user_password": source_instance._master_user_password,
            "kms_key_id": source_instance.kms_key_id,
        }

        if isinstance(source_instance, DBInstance):
            db_instance_specific_kwargs = {
                "storage_type": target_storage_type or source_instance.storage_type,
                "allocated_storage": target_allocated_storage
                or source_instance.allocated_storage,
                "db_instance_identifier": self._generate_green_id(
                    source_instance.db_instance_identifier
                ),
                "db_instance_class": target_db_instance_class
                or source_instance.db_instance_class,
                "storage_throughput": target_storage_throughput
                or source_instance.storage_throughput
                or None,
                "max_allocated_storage": source_instance.max_allocated_storage,
                "db_security_groups": source_instance.db_security_groups,
                "db_cluster_identifier": source_instance.db_cluster_identifier,
                "publicly_accessible": source_instance.publicly_accessible,
                "source_db_instance_identifier": source_instance.source_db_instance_identifier,
                "vpc_security_group_ids": source_instance.vpc_security_group_ids,
                "ca_certificate_identifier": source_instance.ca_certificate_identifier,
                "availability_zone": source_instance.availability_zone,
                "db_subnet_group_name": source_instance._db_subnet_group_name,
                "db_name": source_instance.db_name,
            }
            if (
                target_db_parameter_group_name
                or source_instance.db_parameter_group_name
            ):
                db_instance_specific_kwargs["db_parameter_group_name"] = (
                    target_db_parameter_group_name
                    or source_instance.db_parameter_group_name
                )

            green_instance = self.backend.create_db_instance(
                db_kwargs | db_instance_specific_kwargs
            )
            return green_instance.db_instance_arn
        else:
            db_cluster_specific_kwargs = {
                "db_cluster_identifier": self._generate_green_id(
                    source_instance.db_cluster_identifier
                ),
                "allocated_storage": source_instance.allocated_storage,
                "database_name": source_instance.database_name,
                "db_cluster_parameter_group_name": target_db_cluster_parameter_group_name
                or source_instance.db_cluster_parameter_group,
                "vpc_security_group_ids": source_instance._vpc_security_group_ids,
            }

            green_instance = self.backend.create_db_cluster(
                db_kwargs | db_cluster_specific_kwargs
            )

            for cluster_member in source_instance.members:
                cluster_member_specific_kwargs = {
                    "storage_type": cluster_member.storage_type,
                    "allocated_storage": cluster_member.allocated_storage,
                    "db_instance_identifier": self._generate_green_id(
                        cluster_member.db_instance_identifier
                    ),
                    "db_cluster_identifier": green_instance.db_cluster_identifier,
                    "db_instance_class": cluster_member.db_instance_class,
                    "max_allocated_storage": cluster_member.max_allocated_storage,
                    "publicly_accessible": cluster_member.publicly_accessible,
                    "availability_zone": cluster_member.availability_zone,
                    "db_subnet_group_name": cluster_member._db_subnet_group_name,
                }
                if (
                    target_db_parameter_group_name
                    or cluster_member.db_parameter_group_name
                ):
                    cluster_member_specific_kwargs["db_parameter_group_name"] = (
                        target_db_parameter_group_name
                        or cluster_member.db_parameter_group_name
                    )

                self.backend.create_db_instance(
                    db_kwargs | cluster_member_specific_kwargs
                )
            return green_instance.db_cluster_arn

    def _get_valid_instance_or_raise(self, arn: str) -> DBInstance | DBCluster:
        result = re.findall(r"^arn:[A-Za-z][0-9A-Za-z-:._]*", arn)
        if len(result) == 0:
            raise InvalidParameterValue("Provided string is not a valid arn")
        if self.backend._is_cluster(arn):
            cluster = find_cluster(arn)
            if hasattr(cluster, "master_user_secret"):
                raise SourceClusterNotSupportedFault(cluster.arn)
            return cluster
        else:
            instance = self.backend.find_db_from_id(arn)
            if hasattr(instance, "master_user_secret"):
                raise SourceDatabaseNotSupportedFault(instance.arn)
            return instance

    def switchover(self) -> None:
        source = self._get_valid_instance_or_raise(self.source)
        target = self._get_valid_instance_or_raise(self.target)
        source_result: DBCluster | DBInstance
        target_result: DBCluster | DBInstance

        if isinstance(source, DBCluster):
            new_target_name = source.db_cluster_identifier
            source_result = self.backend.modify_db_cluster(
                {
                    "db_cluster_identifier": source.db_cluster_identifier,
                    "new_db_cluster_identifier": f"{source.db_cluster_identifier}-old1",
                }
            )
        else:
            new_target_name = source.db_instance_identifier
            source_result = self.backend.modify_db_instance(
                db_instance_identifier=source.db_instance_identifier,
                db_kwargs={
                    "new_db_instance_identifier": f"{source.db_instance_identifier}-old1"
                },
            )

        if isinstance(target, DBCluster):
            target_result = self.backend.modify_db_cluster(
                {
                    "db_cluster_identifier": target.db_cluster_identifier,
                    "new_db_cluster_identifier": new_target_name,
                }
            )
        else:
            target_result = self.backend.modify_db_instance(
                db_instance_identifier=target.db_instance_identifier,
                db_kwargs={"new_db_instance_identifier": new_target_name},
            )
        self.source = source_result.arn
        self.target = target_result.arn
        self.set_status(self.SupportedStates.SWITCHOVER_COMPLETED)
        self.status_details = "Switchover completed"


@dataclass
class SwitchoverDetails:
    SourceMember: str
    TargetMember: str
    Status: str


@dataclass
class BlueGreenDeploymentTask:
    Name: str
    Status: str

    class SupportedStates(Enum):
        PENDING = 1
        IN_PROGRESS = 2
        COMPLETED = 3
        FAILED = 4


rds_backends = BackendDict(RDSBackend, "rds")
