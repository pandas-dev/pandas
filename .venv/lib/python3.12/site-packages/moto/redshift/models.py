from __future__ import annotations

from collections import OrderedDict
from collections.abc import Iterable
from datetime import datetime, timedelta, timezone
from typing import Any, Optional

from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import BaseModel, CloudFormationModel
from moto.core.utils import utcnow
from moto.ec2 import ec2_backends
from moto.ec2.models.security_groups import SecurityGroup as EC2SecurityGroup
from moto.moto_api._internal import mock_random
from moto.utilities.utils import get_partition

from .exceptions import (
    ClusterAlreadyExistsFaultError,
    ClusterNotFoundError,
    ClusterParameterGroupNotFoundError,
    ClusterSecurityGroupNotFoundError,
    ClusterSnapshotAlreadyExistsError,
    ClusterSnapshotNotFoundError,
    ClusterSubnetGroupNotFoundError,
    InvalidClusterSnapshotStateFaultError,
    InvalidParameterCombinationError,
    InvalidParameterValueError,
    InvalidSubnetError,
    ResourceNotFoundFaultError,
    SnapshotCopyAlreadyDisabledFaultError,
    SnapshotCopyAlreadyEnabledFaultError,
    SnapshotCopyDisabledFaultError,
    SnapshotCopyGrantAlreadyExistsFaultError,
    SnapshotCopyGrantNotFoundFaultError,
    UnknownSnapshotCopyRegionFaultError,
)


class TaggableResourceMixin:
    resource_type = ""

    def __init__(
        self, account_id: str, region_name: str, tags: Optional[list[dict[str, Any]]]
    ):
        self.account_id = account_id
        self.region = region_name
        self.tags = tags or []

    @property
    def resource_id(self) -> str:
        return ""

    @property
    def arn(self) -> str:
        return f"arn:{get_partition(self.region)}:redshift:{self.region}:{self.account_id}:{self.resource_type}:{self.resource_id}"

    def create_tags(self, tags: list[dict[str, str]]) -> list[dict[str, str]]:
        new_keys = [tag_set["Key"] for tag_set in tags]
        self.tags = [tag_set for tag_set in self.tags if tag_set["Key"] not in new_keys]
        self.tags.extend(tags)
        return self.tags

    def delete_tags(self, tag_keys: list[str]) -> list[dict[str, str]]:
        self.tags = [tag_set for tag_set in self.tags if tag_set["Key"] not in tag_keys]
        return self.tags


class Cluster(TaggableResourceMixin, CloudFormationModel):
    resource_type = "cluster"

    def __init__(
        self,
        redshift_backend: RedshiftBackend,
        cluster_identifier: str,
        node_type: str,
        master_username: str,
        master_user_password: str,
        db_name: str,
        cluster_type: str,
        cluster_security_groups: list[str],
        vpc_security_group_ids: list[str],
        cluster_subnet_group_name: str,
        availability_zone: str,
        preferred_maintenance_window: str,
        cluster_parameter_group_name: str,
        automated_snapshot_retention_period: str,
        port: str,
        cluster_version: str,
        allow_version_upgrade: bool,
        number_of_nodes: str,
        publicly_accessible: bool,
        encrypted: bool,
        region_name: str,
        tags: Optional[list[dict[str, str]]] = None,
        iam_roles_arn: Optional[list[str]] = None,
        enhanced_vpc_routing: Optional[bool] = False,
        restored_from_snapshot: bool = False,
        kms_key_id: Optional[str] = None,
    ):
        super().__init__(redshift_backend.account_id, region_name, tags)
        self.redshift_backend = redshift_backend
        self.cluster_identifier = cluster_identifier
        self.unique_cluster_id = mock_random.get_random_hex(length=12)
        self.create_time = utcnow()
        self.status = "available"
        self.node_type = node_type
        self.master_username = master_username
        self.master_user_password = master_user_password
        self.db_name = db_name if db_name else "dev"
        self.vpc_security_group_ids = vpc_security_group_ids
        self.enhanced_vpc_routing = enhanced_vpc_routing
        self.cluster_subnet_group_name = cluster_subnet_group_name
        self.publicly_accessible = publicly_accessible
        self.encrypted = encrypted

        self.allow_version_upgrade = allow_version_upgrade
        self.cluster_version = cluster_version if cluster_version else "1.0"
        self.port = int(port) if port else 5439
        self.automated_snapshot_retention_period = (
            int(automated_snapshot_retention_period)
            if automated_snapshot_retention_period
            else 1
        )
        self.preferred_maintenance_window = (
            preferred_maintenance_window
            if preferred_maintenance_window
            else "Mon:03:00-Mon:03:30"
        )

        if cluster_parameter_group_name:
            self.cluster_parameter_group_name = [cluster_parameter_group_name]
        else:
            self.cluster_parameter_group_name = ["default.redshift-1.0"]

        if cluster_security_groups:
            self.cluster_security_groups = cluster_security_groups
        else:
            self.cluster_security_groups = ["Default"]

        if availability_zone:
            self.availability_zone = availability_zone
        else:
            # This could probably be smarter, but there doesn't appear to be a
            # way to pull AZs for a region in boto
            self.availability_zone = region_name + "a"

        if cluster_type == "single-node":
            self.number_of_nodes = 1
        elif number_of_nodes:
            self.number_of_nodes = int(number_of_nodes)
        else:
            self.number_of_nodes = 1

        self.iam_roles_arn = iam_roles_arn or []
        self.restored_from_snapshot = restored_from_snapshot
        self.kms_key_id = kms_key_id
        self.cluster_snapshot_copy_status: Optional[dict[str, Any]] = None
        self.total_storage_capacity = 0
        self.logging_details = {
            "LoggingEnabled": "false",  # Lower case is required in response so we use string to simplify
            "BucketName": "",
            "S3KeyPrefix": "",
            "LastSuccessfulDeliveryTime": datetime.now(),
            "LastFailureTime": datetime.now(),
            "LastFailureMessage": "",
            "LogDestinationType": "",
            "LogExports": [],
        }

    @staticmethod
    def cloudformation_name_type() -> str:
        return ""

    @staticmethod
    def cloudformation_type() -> str:
        # https://docs.aws.amazon.com/AWSCloudFormation/latest/UserGuide/aws-resource-redshift-cluster.html
        return "AWS::Redshift::Cluster"

    @classmethod
    def create_from_cloudformation_json(
        cls,
        resource_name: str,
        cloudformation_json: Any,
        account_id: str,
        region_name: str,
        **kwargs: Any,
    ) -> Cluster:
        redshift_backend = redshift_backends[account_id][region_name]
        properties = cloudformation_json["Properties"]

        if "ClusterSubnetGroupName" in properties:
            subnet_group_name = properties[
                "ClusterSubnetGroupName"
            ].cluster_subnet_group_name
        else:
            subnet_group_name = None

        cluster = redshift_backend.create_cluster(
            cluster_identifier=resource_name,
            node_type=properties.get("NodeType"),
            master_username=properties.get("MasterUsername"),
            master_user_password=properties.get("MasterUserPassword"),
            db_name=properties.get("DBName"),
            cluster_type=properties.get("ClusterType"),
            cluster_security_groups=properties.get("ClusterSecurityGroups", []),
            vpc_security_group_ids=properties.get("VpcSecurityGroupIds", []),
            cluster_subnet_group_name=subnet_group_name,
            availability_zone=properties.get("AvailabilityZone"),
            preferred_maintenance_window=properties.get("PreferredMaintenanceWindow"),
            cluster_parameter_group_name=properties.get("ClusterParameterGroupName"),
            automated_snapshot_retention_period=properties.get(
                "AutomatedSnapshotRetentionPeriod"
            ),
            port=properties.get("Port"),
            cluster_version=properties.get("ClusterVersion"),
            allow_version_upgrade=properties.get("AllowVersionUpgrade"),
            enhanced_vpc_routing=properties.get("EnhancedVpcRouting"),
            number_of_nodes=properties.get("NumberOfNodes"),
            publicly_accessible=properties.get("PubliclyAccessible"),
            encrypted=properties.get("Encrypted"),
            region_name=region_name,
            kms_key_id=properties.get("KmsKeyId"),
        )
        return cluster

    @classmethod
    def has_cfn_attr(cls, attr: str) -> bool:
        return attr in ["Endpoint.Address", "Endpoint.Port"]

    def get_cfn_attribute(self, attribute_name: str) -> Any:
        from moto.cloudformation.exceptions import UnformattedGetAttTemplateException

        if attribute_name == "Endpoint.Address":
            return self.endpoint
        if attribute_name == "Endpoint.Port":
            return self.port
        raise UnformattedGetAttTemplateException()

    @property
    def address(self) -> str:
        return f"{self.cluster_identifier}.{self.unique_cluster_id}.{self.region}.redshift.amazonaws.com"

    @property
    def security_groups(self) -> list[SecurityGroup]:
        return [
            security_group
            for security_group in self.redshift_backend.describe_cluster_security_groups()
            if security_group.cluster_security_group_name
            in self.cluster_security_groups
        ]

    @property
    def vpc_security_groups(self) -> list[EC2SecurityGroup]:
        return [
            security_group
            for security_group in self.redshift_backend.ec2_backend.describe_security_groups()
            if security_group.id in self.vpc_security_group_ids
        ]

    @property
    def parameter_groups(self) -> list[ParameterGroup]:
        return [
            parameter_group
            for parameter_group in self.redshift_backend.describe_cluster_parameter_groups()
            if parameter_group.name in self.cluster_parameter_group_name
        ]

    @property
    def resource_id(self) -> str:
        return self.cluster_identifier

    def pause(self) -> None:
        self.status = "paused"

    def resume(self) -> None:
        self.status = "available"

    @property
    def vpc_security_group_membership_list(self) -> list[dict[str, str]]:
        return [
            {"VpcSecurityGroupId": group.id, "Status": "active"}
            for group in self.vpc_security_groups
        ]

    @property
    def cluster_parameter_group_status_list(self) -> list[dict[str, str]]:
        return [
            {
                "ParameterGroupName": group.name,
                "ParameterApplyStatus": "in-sync",
            }
            for group in self.parameter_groups
        ]

    @property
    def cluster_security_group_membership_list(self) -> list[dict[str, str]]:
        return [
            {
                "ClusterSecurityGroupName": group.cluster_security_group_name,
                "Status": "active",
            }
            for group in self.security_groups
        ]

    @property
    def endpoint(self) -> dict[str, str | int]:
        return {
            "Address": self.address,
            "Port": self.port,
        }

    @property
    def pending_modified_values(self) -> list[str]:
        return []

    @property
    def iam_roles(self) -> list[dict[str, str]]:
        return [
            {"ApplyStatus": "in-sync", "IamRoleArn": iam_role_arn}
            for iam_role_arn in self.iam_roles_arn
        ]

    @property
    def total_storage_capacity_in_mega_bytes(self) -> int:
        return self.total_storage_capacity

    @property
    def restore_status(self) -> Optional[dict[str, Any]]:
        if not self.restored_from_snapshot:
            return None
        status = {
            "Status": "completed",
            "CurrentRestoreRateInMegaBytesPerSecond": 123.0,
            "SnapshotSizeInMegaBytes": 123,
            "ProgressInMegaBytes": 123,
            "ElapsedTimeInSeconds": 123,
            "EstimatedTimeToCompletionInSeconds": 123,
        }
        return status


class SnapshotCopyGrant(TaggableResourceMixin, BaseModel):
    resource_type = "snapshotcopygrant"

    def __init__(self, snapshot_copy_grant_name: str, kms_key_id: str):
        self.snapshot_copy_grant_name = snapshot_copy_grant_name
        self.kms_key_id = kms_key_id


class SubnetGroup(TaggableResourceMixin, CloudFormationModel):
    resource_type = "subnetgroup"

    def __init__(
        self,
        ec2_backend: Any,
        cluster_subnet_group_name: str,
        description: str,
        subnet_ids: list[str],
        region_name: str,
        tags: Optional[list[dict[str, str]]] = None,
    ):
        super().__init__(ec2_backend.account_id, region_name, tags)
        self.ec2_backend = ec2_backend
        self.cluster_subnet_group_name = cluster_subnet_group_name
        self.description = description
        self.subnet_ids = subnet_ids
        self.status = "Complete"
        if not self.subnets:
            raise InvalidSubnetError(subnet_ids)

    @staticmethod
    def cloudformation_name_type() -> str:
        return ""

    @staticmethod
    def cloudformation_type() -> str:
        # https://docs.aws.amazon.com/AWSCloudFormation/latest/UserGuide/aws-resource-redshift-clustersubnetgroup.html
        return "AWS::Redshift::ClusterSubnetGroup"

    @classmethod
    def create_from_cloudformation_json(
        cls,
        resource_name: str,
        cloudformation_json: Any,
        account_id: str,
        region_name: str,
        **kwargs: Any,
    ) -> SubnetGroup:
        redshift_backend = redshift_backends[account_id][region_name]
        properties = cloudformation_json["Properties"]

        subnet_group = redshift_backend.create_cluster_subnet_group(
            cluster_subnet_group_name=resource_name,
            description=properties.get("Description"),
            subnet_ids=properties.get("SubnetIds", []),
            region_name=region_name,
        )
        return subnet_group

    @property
    def subnets(self) -> Any:
        return self.ec2_backend.describe_subnets(filters={"subnet-id": self.subnet_ids})

    @property
    def vpc_id(self) -> str:
        return self.subnets[0].vpc_id

    @property
    def resource_id(self) -> str:
        return self.cluster_subnet_group_name

    @property
    def subnet_list(self) -> list[dict[str, Any]]:
        return [
            {
                "SubnetStatus": "Active",
                "SubnetIdentifier": subnet.id,
                "SubnetAvailabilityZone": {"Name": subnet.availability_zone},
            }
            for subnet in self.subnets
        ]


class SecurityGroup(TaggableResourceMixin, BaseModel):
    resource_type = "securitygroup"

    def __init__(
        self,
        cluster_security_group_name: str,
        description: str,
        account_id: str,
        region_name: str,
        tags: Optional[list[dict[str, str]]] = None,
    ):
        super().__init__(account_id, region_name, tags)
        self.cluster_security_group_name = cluster_security_group_name
        self.description = description
        self.ingress_rules: list[str] = []
        self.ec2_security_groups: list[str] = []
        self.ip_ranges: list[str] = []

    @property
    def resource_id(self) -> str:
        return self.cluster_security_group_name


class ParameterGroup(TaggableResourceMixin, CloudFormationModel):
    resource_type = "parametergroup"

    def __init__(
        self,
        cluster_parameter_group_name: str,
        group_family: str,
        description: str,
        account_id: str,
        region_name: str,
        tags: Optional[list[dict[str, str]]] = None,
    ):
        super().__init__(account_id, region_name, tags)
        self.name = cluster_parameter_group_name
        self.family = group_family
        self.description = description

    @staticmethod
    def cloudformation_name_type() -> str:
        return ""

    @staticmethod
    def cloudformation_type() -> str:
        # https://docs.aws.amazon.com/AWSCloudFormation/latest/UserGuide/aws-resource-redshift-clusterparametergroup.html
        return "AWS::Redshift::ClusterParameterGroup"

    @classmethod
    def create_from_cloudformation_json(
        cls,
        resource_name: str,
        cloudformation_json: Any,
        account_id: str,
        region_name: str,
        **kwargs: Any,
    ) -> ParameterGroup:
        redshift_backend = redshift_backends[account_id][region_name]
        properties = cloudformation_json["Properties"]

        parameter_group = redshift_backend.create_cluster_parameter_group(
            cluster_parameter_group_name=resource_name,
            description=properties.get("Description"),
            group_family=properties.get("ParameterGroupFamily"),
            region_name=region_name,
        )
        return parameter_group

    @property
    def resource_id(self) -> str:
        return self.name


class Snapshot(TaggableResourceMixin, BaseModel):
    resource_type = "snapshot"

    def __init__(
        self,
        cluster: Any,
        snapshot_identifier: str,
        account_id: str,
        region_name: str,
        tags: Optional[list[dict[str, str]]] = None,
        iam_roles_arn: Optional[list[str]] = None,
        snapshot_type: str = "manual",
    ):
        super().__init__(account_id, region_name, tags)
        self.snapshot_identifier = snapshot_identifier
        self.snapshot_type = snapshot_type
        self.status = "available"
        self.create_time = utcnow()
        self.iam_roles_arn = iam_roles_arn or []
        self.cluster_identifier = cluster.cluster_identifier
        self.port = cluster.port
        self.availability_zone = cluster.availability_zone
        self.master_username = cluster.master_username
        self.master_user_password = cluster.master_user_password
        self.cluster_version = cluster.cluster_version
        self.node_type = cluster.node_type
        self.number_of_nodes = cluster.number_of_nodes
        self.db_name = cluster.db_name
        self.enhanced_vpc_routing = cluster.enhanced_vpc_routing
        self.encrypted = cluster.encrypted

    @property
    def resource_id(self) -> str:
        return f"{self.cluster_identifier}/{self.snapshot_identifier}"

    @property
    def iam_roles(self) -> list[dict[str, str]]:
        return [
            {"ApplyStatus": "in-sync", "IamRoleArn": iam_role_arn}
            for iam_role_arn in self.iam_roles_arn
        ]


class RedshiftBackend(BaseBackend):
    def __init__(self, region_name: str, account_id: str):
        super().__init__(region_name, account_id)
        self.clusters: dict[str, Cluster] = {}
        self.subnet_groups: dict[str, SubnetGroup] = {}
        self.security_groups: dict[str, SecurityGroup] = {
            "Default": SecurityGroup(
                "Default", "Default Redshift Security Group", account_id, region_name
            )
        }
        self.parameter_groups: dict[str, ParameterGroup] = {
            "default.redshift-1.0": ParameterGroup(
                "default.redshift-1.0",
                "redshift-1.0",
                "Default Redshift parameter group",
                self.account_id,
                self.region_name,
            )
        }
        self.ec2_backend = ec2_backends[self.account_id][self.region_name]
        self.snapshots: dict[str, Snapshot] = OrderedDict()
        self.RESOURCE_TYPE_MAP: dict[str, dict[str, TaggableResourceMixin]] = {
            "cluster": self.clusters,  # type: ignore
            "parametergroup": self.parameter_groups,  # type: ignore
            "securitygroup": self.security_groups,  # type: ignore
            "snapshot": self.snapshots,  # type: ignore
            "subnetgroup": self.subnet_groups,  # type: ignore
        }
        self.snapshot_copy_grants: dict[str, SnapshotCopyGrant] = {}
        self.default_params = {
            "auto_analyze": "true",
            "datestyle": "ISO, MDY",
            "enable_user_activity_logging": "false",
            "extra_float_digits": "0",
            "max_concurrency_scaling_clusters": "1",
            "query_group": "default",
            "require_ssl": "false",
            "search_path": "$user, public",
            "statement_timeout": "0",
            "wlm_json_configuration": '[{"auto_wlm":true}]',
            "use_fips_ssl": "false",
        }

    def enable_snapshot_copy(self, **kwargs: Any) -> Cluster:
        cluster_identifier = kwargs["cluster_identifier"]
        cluster = self.clusters[cluster_identifier]
        if cluster.cluster_snapshot_copy_status is None:
            if cluster.encrypted and kwargs["snapshot_copy_grant_name"] is None:
                raise InvalidParameterValueError(
                    "SnapshotCopyGrantName is required for Snapshot Copy on KMS encrypted clusters."
                )
            if kwargs["destination_region"] == self.region_name:
                raise UnknownSnapshotCopyRegionFaultError(
                    f"Invalid region {self.region_name}"
                )
            status = {
                "DestinationRegion": kwargs["destination_region"],
                "RetentionPeriod": kwargs["retention_period"],
                "SnapshotCopyGrantName": kwargs["snapshot_copy_grant_name"],
            }
            cluster.cluster_snapshot_copy_status = status
            return cluster
        raise SnapshotCopyAlreadyEnabledFaultError(cluster_identifier)

    def disable_snapshot_copy(self, **kwargs: Any) -> Cluster:
        cluster_identifier = kwargs["cluster_identifier"]
        cluster = self.clusters[cluster_identifier]
        if cluster.cluster_snapshot_copy_status is not None:
            cluster.cluster_snapshot_copy_status = None
            return cluster
        raise SnapshotCopyAlreadyDisabledFaultError(cluster_identifier)

    def modify_snapshot_copy_retention_period(
        self, cluster_identifier: str, retention_period: str
    ) -> Cluster:
        cluster = self.clusters[cluster_identifier]
        if cluster.cluster_snapshot_copy_status is not None:
            cluster.cluster_snapshot_copy_status["RetentionPeriod"] = retention_period
            return cluster
        else:
            raise SnapshotCopyDisabledFaultError(cluster_identifier)

    def create_cluster(self, **cluster_kwargs: Any) -> Cluster:
        cluster_identifier = cluster_kwargs["cluster_identifier"]
        if cluster_identifier in self.clusters:
            raise ClusterAlreadyExistsFaultError()
        cluster = Cluster(self, **cluster_kwargs)
        self.clusters[cluster_identifier] = cluster
        snapshot_id = f"rs:{cluster_identifier}-{datetime.now(timezone.utc).strftime('%Y-%m-%d-%H-%M')}"
        # Automated snapshots don't copy over the tags
        self.create_cluster_snapshot(
            cluster_identifier,
            snapshot_id,
            cluster.region,
            None,
            snapshot_type="automated",
        )
        return cluster

    def pause_cluster(self, cluster_id: str) -> Cluster:
        if cluster_id not in self.clusters:
            raise ClusterNotFoundError(cluster_identifier=cluster_id)
        self.clusters[cluster_id].pause()
        return self.clusters[cluster_id]

    def resume_cluster(self, cluster_id: str) -> Cluster:
        if cluster_id not in self.clusters:
            raise ClusterNotFoundError(cluster_identifier=cluster_id)
        self.clusters[cluster_id].resume()
        return self.clusters[cluster_id]

    def describe_clusters(
        self,
        cluster_identifier: Optional[str] = None,
        tag_keys: Optional[list[str]] = None,
    ) -> list[Cluster]:
        if cluster_identifier:
            if cluster_identifier in self.clusters:
                return [self.clusters[cluster_identifier]]
            raise ClusterNotFoundError(cluster_identifier)
        elif tag_keys:
            return [
                cluster
                for cluster in self.clusters.values()
                if any(key in [tag["Key"] for tag in cluster.tags] for key in tag_keys)
            ]
        return list(self.clusters.values())

    def modify_cluster(self, **cluster_kwargs: Any) -> Cluster:
        cluster_identifier = cluster_kwargs.pop("cluster_identifier")
        new_cluster_identifier = cluster_kwargs.pop("new_cluster_identifier", None)

        cluster_type = cluster_kwargs.get("cluster_type")
        if cluster_type and cluster_type not in ["multi-node", "single-node"]:
            raise InvalidParameterValueError(
                "Invalid cluster type. Cluster type can be one of multi-node or single-node"
            )
        if cluster_type == "single-node":
            # AWS will always silently override this value for single-node clusters.
            cluster_kwargs["number_of_nodes"] = 1
        elif cluster_type == "multi-node":
            if cluster_kwargs.get("number_of_nodes", 0) < 2:
                raise InvalidParameterCombinationError(
                    "Number of nodes for cluster type multi-node must be greater than or equal to 2"
                )

        cluster = self.describe_clusters(cluster_identifier)[0]

        for key, value in cluster_kwargs.items():
            setattr(cluster, key, value)

        if new_cluster_identifier:
            dic = {
                "cluster_identifier": cluster_identifier,
                "skip_final_snapshot": True,
                "final_cluster_snapshot_identifier": None,
            }
            self.delete_cluster(**dic)
            cluster.cluster_identifier = new_cluster_identifier
            self.clusters[new_cluster_identifier] = cluster

        return cluster

    def delete_automated_snapshots(self, cluster_identifier: str) -> None:
        snapshots = self.describe_cluster_snapshots(
            cluster_identifier=cluster_identifier
        )
        for snapshot in snapshots:
            if snapshot.snapshot_type == "automated":
                self.snapshots.pop(snapshot.snapshot_identifier)

    def delete_cluster(self, **cluster_kwargs: Any) -> Cluster:
        cluster_identifier = cluster_kwargs.pop("cluster_identifier")
        cluster_skip_final_snapshot = cluster_kwargs.pop("skip_final_snapshot")
        cluster_snapshot_identifer = cluster_kwargs.pop(
            "final_cluster_snapshot_identifier"
        )

        if cluster_identifier in self.clusters:
            if (
                cluster_skip_final_snapshot is False
                and cluster_snapshot_identifer is None
            ):
                raise InvalidParameterCombinationError(
                    "FinalClusterSnapshotIdentifier is required unless "
                    "SkipFinalClusterSnapshot is specified."
                )
            if (
                cluster_skip_final_snapshot is False
                and cluster_snapshot_identifer is not None
            ):  # create snapshot
                cluster = self.describe_clusters(cluster_identifier)[0]
                self.create_cluster_snapshot(
                    cluster_identifier,
                    cluster_snapshot_identifer,
                    cluster.region,
                    cluster.tags,
                )
            self.delete_automated_snapshots(cluster_identifier)
            return self.clusters.pop(cluster_identifier)
        raise ClusterNotFoundError(cluster_identifier)

    def create_cluster_subnet_group(
        self,
        cluster_subnet_group_name: str,
        description: str,
        subnet_ids: list[str],
        region_name: str,
        tags: Optional[list[dict[str, str]]] = None,
    ) -> SubnetGroup:
        subnet_group = SubnetGroup(
            self.ec2_backend,
            cluster_subnet_group_name,
            description,
            subnet_ids,
            region_name,
            tags,
        )
        self.subnet_groups[cluster_subnet_group_name] = subnet_group
        return subnet_group

    def describe_cluster_subnet_groups(
        self, subnet_identifier: Optional[str] = None
    ) -> list[SubnetGroup]:
        if subnet_identifier:
            if subnet_identifier in self.subnet_groups:
                return [self.subnet_groups[subnet_identifier]]
            raise ClusterSubnetGroupNotFoundError(subnet_identifier)
        return list(self.subnet_groups.values())

    def delete_cluster_subnet_group(self, subnet_identifier: str) -> SubnetGroup:
        if subnet_identifier in self.subnet_groups:
            return self.subnet_groups.pop(subnet_identifier)
        raise ClusterSubnetGroupNotFoundError(subnet_identifier)

    def create_cluster_security_group(
        self,
        cluster_security_group_name: str,
        description: str,
        tags: Optional[list[dict[str, str]]] = None,
    ) -> SecurityGroup:
        security_group = SecurityGroup(
            cluster_security_group_name,
            description,
            self.account_id,
            self.region_name,
            tags,
        )
        self.security_groups[cluster_security_group_name] = security_group
        return security_group

    def describe_cluster_security_groups(
        self, security_group_name: Optional[str] = None
    ) -> list[SecurityGroup]:
        if security_group_name:
            if security_group_name in self.security_groups:
                return [self.security_groups[security_group_name]]
            raise ClusterSecurityGroupNotFoundError(security_group_name)
        return list(self.security_groups.values())

    def delete_cluster_security_group(
        self, security_group_identifier: str
    ) -> SecurityGroup:
        if security_group_identifier in self.security_groups:
            return self.security_groups.pop(security_group_identifier)
        raise ClusterSecurityGroupNotFoundError(security_group_identifier)

    def authorize_cluster_security_group_ingress(
        self, security_group_name: str, cidr_ip: str
    ) -> SecurityGroup:
        security_group = self.security_groups.get(security_group_name)
        if not security_group:
            raise ClusterSecurityGroupNotFoundError()

        # just adding the cidr_ip as ingress rule for now as there is no security rule
        security_group.ingress_rules.append(cidr_ip)

        return security_group

    def create_cluster_parameter_group(
        self,
        cluster_parameter_group_name: str,
        group_family: str,
        description: str,
        region_name: str,
        tags: Optional[list[dict[str, str]]] = None,
    ) -> ParameterGroup:
        parameter_group = ParameterGroup(
            cluster_parameter_group_name,
            group_family,
            description,
            self.account_id,
            region_name,
            tags,
        )
        self.parameter_groups[cluster_parameter_group_name] = parameter_group

        return parameter_group

    def describe_cluster_parameter_groups(
        self, parameter_group_name: Optional[str] = None
    ) -> list[ParameterGroup]:
        if parameter_group_name:
            if parameter_group_name in self.parameter_groups:
                return [self.parameter_groups[parameter_group_name]]
            raise ClusterParameterGroupNotFoundError(parameter_group_name)
        return list(self.parameter_groups.values())

    def delete_cluster_parameter_group(
        self, parameter_group_name: str
    ) -> ParameterGroup:
        if parameter_group_name in self.parameter_groups:
            return self.parameter_groups.pop(parameter_group_name)
        raise ClusterParameterGroupNotFoundError(parameter_group_name)

    def describe_default_cluster_parameters(self) -> list[dict[str, Any]]:
        return [
            {
                "ParameterName": key,
                "ParameterValue": value,
                "Description": "mock parameter",
                "Source": "engine-default",
                "DataType": "type",
                "ApplyType": "static",
                "IsModifiable": "true",
            }
            for key, value in self.default_params.items()
        ]

    def describe_cluster_parameters(
        self, parameter_group_name: str
    ) -> list[dict[str, Any]]:
        if parameter_group_name not in self.parameter_groups:
            raise ClusterParameterGroupNotFoundError(parameter_group_name)
        return self.describe_default_cluster_parameters()

    def create_cluster_snapshot(
        self,
        cluster_identifier: str,
        snapshot_identifier: str,
        region_name: str,
        tags: Optional[list[dict[str, str]]],
        snapshot_type: str = "manual",
    ) -> Snapshot:
        cluster = self.clusters.get(cluster_identifier)
        if not cluster:
            raise ClusterNotFoundError(cluster_identifier)
        if self.snapshots.get(snapshot_identifier) is not None:
            raise ClusterSnapshotAlreadyExistsError(snapshot_identifier)
        snapshot = Snapshot(
            cluster,
            snapshot_identifier,
            self.account_id,
            region_name,
            tags,
            snapshot_type=snapshot_type,
        )
        self.snapshots[snapshot_identifier] = snapshot
        return snapshot

    def describe_cluster_snapshots(
        self,
        cluster_identifier: Optional[str] = None,
        snapshot_identifier: Optional[str] = None,
        snapshot_type: Optional[str] = None,
    ) -> list[Snapshot]:
        snapshot_types = (
            ["automated", "manual"] if snapshot_type is None else [snapshot_type]
        )
        if cluster_identifier:
            cluster_snapshots = []
            for snapshot in self.snapshots.values():
                if snapshot.cluster_identifier == cluster_identifier:
                    if snapshot.snapshot_type in snapshot_types:
                        cluster_snapshots.append(snapshot)
            if cluster_snapshots:
                return cluster_snapshots

        if snapshot_identifier:
            if snapshot_identifier in self.snapshots:
                if self.snapshots[snapshot_identifier].snapshot_type in snapshot_types:
                    return [self.snapshots[snapshot_identifier]]
            raise ClusterSnapshotNotFoundError(snapshot_identifier)

        return list(self.snapshots.values())

    def delete_cluster_snapshot(self, snapshot_identifier: str) -> Snapshot:
        if snapshot_identifier not in self.snapshots:
            raise ClusterSnapshotNotFoundError(snapshot_identifier)

        snapshot = self.describe_cluster_snapshots(
            snapshot_identifier=snapshot_identifier
        )[0]
        if snapshot.snapshot_type == "automated":
            raise InvalidClusterSnapshotStateFaultError(snapshot_identifier)
        deleted_snapshot = self.snapshots.pop(snapshot_identifier)
        deleted_snapshot.status = "deleted"
        return deleted_snapshot

    def restore_from_cluster_snapshot(self, **kwargs: Any) -> Cluster:
        snapshot_identifier = kwargs.pop("snapshot_identifier")
        snapshot = self.describe_cluster_snapshots(
            snapshot_identifier=snapshot_identifier
        )[0]
        create_kwargs = {
            "node_type": snapshot.node_type,
            "master_username": snapshot.master_username,
            "master_user_password": snapshot.master_user_password,
            "db_name": snapshot.db_name,
            "cluster_type": "multi-node"
            if snapshot.number_of_nodes > 1
            else "single-node",
            "availability_zone": snapshot.availability_zone,
            "port": snapshot.port,
            "cluster_version": snapshot.cluster_version,
            "number_of_nodes": snapshot.number_of_nodes,
            "encrypted": snapshot.encrypted,
            "tags": snapshot.tags,
            "restored_from_snapshot": True,
            "enhanced_vpc_routing": snapshot.enhanced_vpc_routing,
        }
        create_kwargs.update(kwargs)
        return self.create_cluster(**create_kwargs)

    def create_snapshot_copy_grant(self, **kwargs: Any) -> SnapshotCopyGrant:
        snapshot_copy_grant_name = kwargs["snapshot_copy_grant_name"]
        kms_key_id = kwargs["kms_key_id"]
        if snapshot_copy_grant_name not in self.snapshot_copy_grants:
            snapshot_copy_grant = SnapshotCopyGrant(
                snapshot_copy_grant_name, kms_key_id
            )
            self.snapshot_copy_grants[snapshot_copy_grant_name] = snapshot_copy_grant
            return snapshot_copy_grant
        raise SnapshotCopyGrantAlreadyExistsFaultError(snapshot_copy_grant_name)

    def delete_snapshot_copy_grant(self, **kwargs: Any) -> SnapshotCopyGrant:
        snapshot_copy_grant_name = kwargs["snapshot_copy_grant_name"]
        if snapshot_copy_grant_name in self.snapshot_copy_grants:
            return self.snapshot_copy_grants.pop(snapshot_copy_grant_name)
        raise SnapshotCopyGrantNotFoundFaultError(snapshot_copy_grant_name)

    def describe_snapshot_copy_grants(self, **kwargs: Any) -> list[SnapshotCopyGrant]:
        copy_grants = list(self.snapshot_copy_grants.values())
        snapshot_copy_grant_name = kwargs["snapshot_copy_grant_name"]
        if snapshot_copy_grant_name:
            if snapshot_copy_grant_name in self.snapshot_copy_grants:
                return [self.snapshot_copy_grants[snapshot_copy_grant_name]]
            raise SnapshotCopyGrantNotFoundFaultError(snapshot_copy_grant_name)
        return copy_grants

    def _get_resource_from_arn(self, arn: str) -> TaggableResourceMixin:
        try:
            arn_breakdown = arn.split(":")
            resource_type = arn_breakdown[5]
            if resource_type == "snapshot":
                resource_id = arn_breakdown[6].split("/")[1]
            else:
                resource_id = arn_breakdown[6]
        except IndexError:
            resource_type = resource_id = arn
        resources = self.RESOURCE_TYPE_MAP.get(resource_type)
        if resources is None:
            message = (
                "Tagging is not supported for this type of resource: "
                f"'{resource_type}' (the ARN is potentially malformed, "
                "please check the ARN documentation for more information)"
            )
            raise ResourceNotFoundFaultError(message=message)
        try:
            resource = resources[resource_id]
        except KeyError:
            raise ResourceNotFoundFaultError(resource_type, resource_id)
        return resource

    @staticmethod
    def _describe_tags_for_resources(resources: Iterable[Any]) -> list[dict[str, Any]]:
        tagged_resources = []
        for resource in resources:
            for tag in resource.tags:
                data = {
                    "ResourceName": resource.arn,
                    "ResourceType": resource.resource_type,
                    "Tag": {"Key": tag["Key"], "Value": tag["Value"]},
                }
                tagged_resources.append(data)
        return tagged_resources

    def _describe_tags_for_resource_type(
        self, resource_type: str
    ) -> list[dict[str, Any]]:
        resources = self.RESOURCE_TYPE_MAP.get(resource_type)
        if not resources:
            raise ResourceNotFoundFaultError(resource_type=resource_type)
        return self._describe_tags_for_resources(resources.values())

    def _describe_tags_for_resource_name(
        self, resource_name: str
    ) -> list[dict[str, Any]]:
        resource = self._get_resource_from_arn(resource_name)
        return self._describe_tags_for_resources([resource])

    def create_tags(self, resource_name: str, tags: list[dict[str, str]]) -> None:
        resource = self._get_resource_from_arn(resource_name)
        resource.create_tags(tags)

    def describe_tags(
        self, resource_name: str, resource_type: str
    ) -> list[dict[str, Any]]:
        if resource_name and resource_type:
            raise InvalidParameterValueError(
                "You cannot filter a list of resources using an Amazon "
                "Resource Name (ARN) and a resource type together in the "
                "same request. Retry the request using either an ARN or "
                "a resource type, but not both."
            )
        if resource_type:
            return self._describe_tags_for_resource_type(resource_type.lower())
        if resource_name:
            return self._describe_tags_for_resource_name(resource_name)
        # If name and type are not specified, return all tagged resources.
        # TODO: Implement aws marker pagination
        tagged_resources = []
        for resource_type in self.RESOURCE_TYPE_MAP:
            try:
                tagged_resources += self._describe_tags_for_resource_type(resource_type)
            except ResourceNotFoundFaultError:
                pass
        return tagged_resources

    def delete_tags(self, resource_name: str, tag_keys: list[str]) -> None:
        resource = self._get_resource_from_arn(resource_name)
        resource.delete_tags(tag_keys)

    def get_cluster_credentials(
        self,
        cluster_identifier: str,
        db_user: str,
        auto_create: bool,
        duration_seconds: int,
    ) -> dict[str, Any]:
        if duration_seconds < 900 or duration_seconds > 3600:
            raise InvalidParameterValueError(
                "Token duration must be between 900 and 3600 seconds"
            )
        expiration = datetime.now(timezone.utc) + timedelta(seconds=duration_seconds)
        if cluster_identifier in self.clusters:
            user_prefix = "IAM:" if auto_create is False else "IAMA:"
            db_user = user_prefix + db_user
            return {
                "DbUser": db_user,
                "DbPassword": mock_random.get_random_string(32),
                "Expiration": expiration,
            }
        raise ClusterNotFoundError(cluster_identifier)

    def enable_logging(
        self,
        cluster_identifier: str,
        bucket_name: str,
        s3_key_prefix: str,
        log_destination_type: str,
        log_exports: list[str],
    ) -> dict[str, Any]:
        if cluster_identifier not in self.clusters:
            raise ClusterNotFoundError(cluster_identifier)
        cluster = self.clusters[cluster_identifier]
        cluster.logging_details["LoggingEnabled"] = "true"
        cluster.logging_details["BucketName"] = bucket_name
        cluster.logging_details["S3KeyPrefix"] = s3_key_prefix
        cluster.logging_details["LogDestinationType"] = log_destination_type
        cluster.logging_details["LogExports"] = log_exports
        return cluster.logging_details

    def disable_logging(self, cluster_identifier: str) -> dict[str, Any]:
        if cluster_identifier not in self.clusters:
            raise ClusterNotFoundError(cluster_identifier)
        cluster = self.clusters[cluster_identifier]
        cluster.logging_details["LoggingEnabled"] = "false"
        return cluster.logging_details

    def describe_logging_status(self, cluster_identifier: str) -> dict[str, Any]:
        if cluster_identifier not in self.clusters:
            raise ClusterNotFoundError(cluster_identifier)
        cluster = self.clusters[cluster_identifier]
        return cluster.logging_details


redshift_backends = BackendDict(RedshiftBackend, "redshift")
