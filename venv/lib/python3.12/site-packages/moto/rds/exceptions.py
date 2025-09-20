from typing import Optional

from moto.core.exceptions import ServiceException


class RDSClientError(ServiceException):
    pass


class DBInstanceNotFoundError(RDSClientError):
    def __init__(self, database_identifier: str):
        super().__init__(
            "DBInstanceNotFound", f"DBInstance {database_identifier} not found."
        )


class DBInstanceAlreadyExists(RDSClientError):
    def __init__(self) -> None:
        super().__init__("DBInstanceAlreadyExists", "DB instance already exists")


class DBSnapshotNotFoundFault(RDSClientError):
    def __init__(self, snapshot_identifier: str):
        super().__init__(
            "DBSnapshotNotFoundFault", f"DBSnapshot {snapshot_identifier} not found."
        )


class DBSecurityGroupNotFoundError(RDSClientError):
    def __init__(self, security_group_name: str):
        super().__init__(
            "DBSecurityGroupNotFound",
            f"Security Group {security_group_name} not found.",
        )


class DBShardGroupAlreadyExistsError(RDSClientError):
    def __init__(self, shard_group_identifier: str):
        super().__init__(
            "DBShardGroupAlreadyExists",
            f"DB Shard Group {shard_group_identifier} already exists.",
        )


class DBSubnetGroupNotFoundError(RDSClientError):
    def __init__(self, subnet_group_name: str):
        super().__init__(
            "DBSubnetGroupNotFoundFault", f"Subnet Group {subnet_group_name} not found."
        )


class DBParameterGroupNotFoundError(RDSClientError):
    def __init__(self, db_parameter_group_name: str):
        super().__init__(
            "DBParameterGroupNotFound",
            f"DB Parameter Group {db_parameter_group_name} not found.",
        )


class DBParameterGroupAlreadyExistsError(RDSClientError):
    def __init__(self, db_parameter_group_name: str):
        super().__init__(
            "DBParameterGroupAlreadyExists",
            f"DB Parameter Group {db_parameter_group_name} already exists.",
        )


class DBClusterParameterGroupNotFoundError(RDSClientError):
    def __init__(self, group_name: str):
        super().__init__(
            "DBParameterGroupNotFound",
            f"DBClusterParameterGroup not found: {group_name}",
        )


class OptionGroupNotFoundFaultError(RDSClientError):
    def __init__(self, option_group_name: str):
        super().__init__(
            "OptionGroupNotFoundFault",
            f"Specified OptionGroupName: {option_group_name} not found.",
        )


class InvalidDBClusterStateFaultError(RDSClientError):
    def __init__(self, database_identifier: str):
        super().__init__(
            "InvalidDBClusterStateFault",
            f"Invalid DB type, when trying to perform StopDBInstance on {database_identifier}e. See AWS RDS documentation on rds.stop_db_instance",
        )


class DBClusterToBeDeletedHasActiveMembers(RDSClientError):
    def __init__(self) -> None:
        super().__init__(
            "InvalidDBClusterStateFault",
            "Cluster cannot be deleted, it still contains DB instances in non-deleting state.",
        )


class InvalidDBInstanceStateError(RDSClientError):
    def __init__(self, database_identifier: str, istate: str):
        estate = (
            "in available state"
            if istate == "stop"
            else "stopped, it cannot be started"
        )
        super().__init__(
            "InvalidDBInstanceState", f"Instance {database_identifier} is not {estate}."
        )


class InvalidDBInstanceStateFault(RDSClientError):
    def __init__(self, message: str):
        super().__init__("InvalidDBInstanceStateFault", message)


class SnapshotQuotaExceededFault(RDSClientError):
    # This is used for both DBSnapshots and DBClusterSnapshots
    def __init__(self) -> None:
        super().__init__(
            "SnapshotQuotaExceeded",
            "The request cannot be processed because it would exceed the maximum number of snapshots.",
        )


class SharedSnapshotQuotaExceeded(RDSClientError):
    def __init__(self) -> None:
        super().__init__(
            "SharedSnapshotQuotaExceeded",
            "The request cannot be processed because it would exceed the maximum number of snapshots.",
        )


class KMSKeyNotAccessibleFault(RDSClientError):
    fmt = "Specified KMS key [{key_id}] does not exist, is not enabled or you do not have permissions to access it."

    def __init__(self, key_id: str) -> None:
        super().__init__(
            "KMSKeyNotAccessibleFault",
            f"Specified KMS key [{key_id}] does not exist, is not enabled or you do not have permissions to access it.",
        )


class InvalidDBClusterSnapshotStateFault(RDSClientError):
    def __init__(self, message: str):
        super().__init__("InvalidDBClusterSnapshotStateFault", message)


class DBSnapshotAlreadyExistsError(RDSClientError):
    def __init__(self, database_snapshot_identifier: str):
        super().__init__(
            "DBSnapshotAlreadyExists",
            f"Cannot create the snapshot because a snapshot with the identifier {database_snapshot_identifier} already exists.",
        )


class DBShardGroupNotFoundFault(RDSClientError):
    def __init__(self, db_shard_group_identifier: str):
        super().__init__(
            "DBShardGroupNotFoundFault",
            f"DBShardGroup {db_shard_group_identifier} not found.",
        )


class InvalidParameterValue(RDSClientError):
    def __init__(self, message: str):
        super().__init__("InvalidParameterValue", message)


class InvalidParameterCombination(RDSClientError):
    def __init__(self, message: str):
        super().__init__("InvalidParameterCombination", message)


class InvalidDBClusterStateFault(RDSClientError):
    def __init__(self, message: str):
        super().__init__("InvalidDBClusterStateFault", message)


class DBClusterNotFoundError(RDSClientError):
    def __init__(self, cluster_identifier: str, message: Optional[str] = None):
        if message is None:
            message = f"DBCluster {cluster_identifier} not found."
        super().__init__("DBClusterNotFoundFault", message)


class DBClusterSnapshotNotFoundError(RDSClientError):
    def __init__(self, snapshot_identifier: str):
        super().__init__(
            "DBClusterSnapshotNotFoundFault",
            f"DBClusterSnapshot {snapshot_identifier} not found.",
        )


class DBClusterSnapshotAlreadyExistsError(RDSClientError):
    def __init__(self, database_snapshot_identifier: str):
        super().__init__(
            "DBClusterSnapshotAlreadyExistsFault",
            f"Cannot create the snapshot because a snapshot with the identifier {database_snapshot_identifier} already exists.",
        )


class ExportTaskAlreadyExistsError(RDSClientError):
    def __init__(self, export_task_identifier: str):
        super().__init__(
            "ExportTaskAlreadyExists",
            f"Cannot start export task because a task with the identifier {export_task_identifier} already exists.",
        )


class ExportTaskNotFoundError(RDSClientError):
    def __init__(self, export_task_identifier: str):
        super().__init__(
            "ExportTaskNotFound",
            f"Cannot cancel export task because a task with the identifier {export_task_identifier} is not exist.",
        )


class InvalidExportSourceStateError(RDSClientError):
    def __init__(self, status: str):
        super().__init__(
            "InvalidExportSourceStateFault",
            f"Export source should be 'available' but current status is {status}.",
        )


class SubscriptionAlreadyExistError(RDSClientError):
    def __init__(self, subscription_name: str):
        super().__init__(
            "SubscriptionAlreadyExist",
            f"Subscription {subscription_name} already exists.",
        )


class SubscriptionNotFoundError(RDSClientError):
    def __init__(self, subscription_name: str):
        super().__init__(
            "SubscriptionNotFound", f"Subscription {subscription_name} not found."
        )


class InvalidGlobalClusterStateFault(RDSClientError):
    def __init__(self, arn: str):
        super().__init__(
            "InvalidGlobalClusterStateFault", f"Global Cluster {arn} is not empty"
        )


class InvalidDBInstanceIdentifier(InvalidParameterValue):
    def __init__(self) -> None:
        super().__init__(
            "The parameter DBInstanceIdentifier is not a valid identifier. "
            "Identifiers must begin with a letter; must contain only ASCII letters, digits, and hyphens; "
            "and must not end with a hyphen or contain two consecutive hyphens."
        )


class InvalidDBSnapshotIdentifier(InvalidParameterValue):
    def __init__(self, snapshot_identifier: str, parameter_name: str) -> None:
        if snapshot_identifier == "":
            exception_text = f"The parameter {parameter_name} must be provided and must not be blank."
        elif not snapshot_identifier[0].isalpha():
            # On AWS, this error message seems to be triggered when the first character is invalid.
            # The two spaces before the snapshot_identifier are what AWS produces!
            exception_text = f"Invalid snapshot identifier:  {snapshot_identifier}"
        else:
            exception_text = (
                f"The parameter {parameter_name} is not a valid identifier. "
                "Identifiers must begin with a letter; must contain only ASCII letters, digits, and hyphens; "
                "and must not end with a hyphen or contain two consecutive hyphens."
            )
        super().__init__(exception_text)


class InvalidDBInstanceEngine(InvalidParameterCombination):
    def __init__(self, instance_engine: str, cluster_engine: str) -> None:
        super().__init__(
            f"The engine name requested for your DB instance ({instance_engine}) doesn't match "
            f"the engine name of your DB cluster ({cluster_engine})."
        )


class InvalidSubnet(RDSClientError):
    def __init__(self, subnet_identifier: str):
        super().__init__(
            "InvalidSubnet",
            f"The requested subnet {subnet_identifier} is invalid, or multiple subnets were requested that are not all in a common VPC.",
        )


class DBProxyAlreadyExistsFault(RDSClientError):
    def __init__(self, db_proxy_identifier: str):
        super().__init__(
            "DBProxyAlreadyExistsFault",
            f"Cannot create the DBProxy because a DBProxy with the identifier {db_proxy_identifier} already exists.",
        )


class DBProxyQuotaExceededFault(RDSClientError):
    def __init__(self) -> None:
        super().__init__(
            "DBProxyQuotaExceeded",
            "The request cannot be processed because it would exceed the maximum number of DBProxies.",
        )


class DBProxyNotFoundFault(RDSClientError):
    def __init__(self, db_proxy_identifier: str):
        super().__init__(
            "DBProxyNotFoundFault",
            f"The specified proxy name {db_proxy_identifier} doesn't correspond to a proxy owned by your Amazon Web Services account in the specified Amazon Web Services Region.",
        )


class BlueGreenDeploymentAlreadyExistsFault(RDSClientError):
    def __init__(self, bg_name: str):
        super().__init__(
            "BlueGreenDeploymentAlreadyExistsFault",
            f"A blue/green deployment with the specified name {bg_name} already exists.",
        )


class BlueGreenDeploymentNotFoundFault(RDSClientError):
    def __init__(self, bg_identifier: str):
        super().__init__(
            "BlueGreenDeploymentNotFoundFault",
            f"BlueGreenDeploymentIdentifier {bg_identifier} doesn't refer to an existing blue/green deployment.",
        )


class InvalidBlueGreenDeploymentStateFault(RDSClientError):
    def __init__(self, bg_identifier: str):
        super().__init__(
            "InvalidBlueGreenDeploymentStateFault",
            f"The blue/green deployment {bg_identifier} can't be switched over or deleted because there is an invalid configuration in the green environment.",
        )


class SourceDatabaseNotSupportedFault(RDSClientError):
    def __init__(self, source_arn: str):
        super().__init__(
            "SourceDatabaseNotSupportedFault",
            f"The source DB instance {source_arn} isn't supported for a blue/green deployment.",
        )


class SourceClusterNotSupportedFault(RDSClientError):
    def __init__(self, source_arn: str):
        super().__init__(
            "SourceClusterNotSupportedFault",
            f"The source DB cluster {source_arn} isn't supported for a blue/green deployment.",
        )
