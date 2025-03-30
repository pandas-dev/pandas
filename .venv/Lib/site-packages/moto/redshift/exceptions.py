import json
from typing import List, Optional

from moto.core.exceptions import JsonRESTError


class RedshiftClientError(JsonRESTError):
    def __init__(self, code: str, message: str):
        super().__init__(error_type=code, message=message)
        self.description = json.dumps(
            {
                "Error": {"Code": code, "Message": message, "Type": "Sender"},
                "RequestId": "6876f774-7273-11e4-85dc-39e55ca848d1",
            }
        )


class ClusterNotFoundError(RedshiftClientError):
    def __init__(self, cluster_identifier: str):
        super().__init__("ClusterNotFound", f"Cluster {cluster_identifier} not found.")


class ClusterSubnetGroupNotFoundError(RedshiftClientError):
    def __init__(self, subnet_identifier: str):
        super().__init__(
            "ClusterSubnetGroupNotFound", f"Subnet group {subnet_identifier} not found."
        )


class ClusterSecurityGroupNotFoundError(RedshiftClientError):
    def __init__(self, group_identifier: str):
        super().__init__(
            "ClusterSecurityGroupNotFound",
            f"Security group {group_identifier} not found.",
        )


class ClusterParameterGroupNotFoundError(RedshiftClientError):
    def __init__(self, group_identifier: str):
        super().__init__(
            "ClusterParameterGroupNotFound",
            f"Parameter group {group_identifier} not found.",
        )


class InvalidSubnetError(RedshiftClientError):
    def __init__(self, subnet_identifier: List[str]):
        super().__init__("InvalidSubnet", f"Subnet {subnet_identifier} not found.")


class SnapshotCopyGrantAlreadyExistsFaultError(RedshiftClientError):
    def __init__(self, snapshot_copy_grant_name: str):
        super().__init__(
            "SnapshotCopyGrantAlreadyExistsFault",
            "Cannot create the snapshot copy grant because a grant "
            f"with the identifier '{snapshot_copy_grant_name}' already exists",
        )


class SnapshotCopyGrantNotFoundFaultError(RedshiftClientError):
    def __init__(self, snapshot_copy_grant_name: str):
        super().__init__(
            "SnapshotCopyGrantNotFoundFault",
            f"Snapshot copy grant not found: {snapshot_copy_grant_name}",
        )


class ClusterSnapshotNotFoundError(RedshiftClientError):
    def __init__(self, snapshot_identifier: str):
        super().__init__(
            "ClusterSnapshotNotFound", f"Snapshot {snapshot_identifier} not found."
        )


class ClusterSnapshotAlreadyExistsError(RedshiftClientError):
    def __init__(self, snapshot_identifier: str):
        super().__init__(
            "ClusterSnapshotAlreadyExists",
            "Cannot create the snapshot because a snapshot with the "
            f"identifier {snapshot_identifier} already exists",
        )


class InvalidParameterValueError(RedshiftClientError):
    def __init__(self, message: str):
        super().__init__("InvalidParameterValue", message)


class ResourceNotFoundFaultError(RedshiftClientError):
    code = 404

    def __init__(
        self,
        resource_type: Optional[str] = None,
        resource_name: Optional[str] = None,
        message: Optional[str] = None,
    ):
        if resource_type and not resource_name:
            msg = f"resource of type '{resource_type}' not found."
        else:
            msg = f"{resource_type} ({resource_name}) not found."
        if message:
            msg = message
        super().__init__("ResourceNotFoundFault", msg)


class SnapshotCopyDisabledFaultError(RedshiftClientError):
    def __init__(self, cluster_identifier: str):
        super().__init__(
            "SnapshotCopyDisabledFault",
            f"Cannot modify retention period because snapshot copy is disabled on Cluster {cluster_identifier}.",
        )


class SnapshotCopyAlreadyDisabledFaultError(RedshiftClientError):
    def __init__(self, cluster_identifier: str):
        super().__init__(
            "SnapshotCopyAlreadyDisabledFault",
            f"Snapshot Copy is already disabled on Cluster {cluster_identifier}.",
        )


class SnapshotCopyAlreadyEnabledFaultError(RedshiftClientError):
    def __init__(self, cluster_identifier: str):
        super().__init__(
            "SnapshotCopyAlreadyEnabledFault",
            f"Snapshot Copy is already enabled on Cluster {cluster_identifier}.",
        )


class ClusterAlreadyExistsFaultError(RedshiftClientError):
    def __init__(self) -> None:
        super().__init__("ClusterAlreadyExists", "Cluster already exists")


class InvalidParameterCombinationError(RedshiftClientError):
    def __init__(self, message: str):
        super().__init__("InvalidParameterCombination", message)


class UnknownSnapshotCopyRegionFaultError(RedshiftClientError):
    def __init__(self, message: str):
        super().__init__("UnknownSnapshotCopyRegionFault", message)


class ClusterSecurityGroupNotFoundFaultError(RedshiftClientError):
    def __init__(self) -> None:
        super().__init__(
            "ClusterSecurityGroupNotFoundFault",
            "The cluster security group name does not refer to an existing cluster security group.",
        )


class InvalidClusterSnapshotStateFaultError(RedshiftClientError):
    def __init__(self, snapshot_identifier: str):
        super().__init__(
            "InvalidClusterSnapshotStateFault",
            f"Cannot delete the snapshot {snapshot_identifier} because only manual snapshots may be deleted",
        )
