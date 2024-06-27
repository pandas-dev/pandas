import json
from typing import Any, Dict, Tuple, Union

from moto.core.responses import BaseResponse

from .models import EFSBackend, efs_backends

TYPE_RESPONSE = Tuple[str, Dict[str, Union[str, int]]]


class EFSResponse(BaseResponse):
    def __init__(self) -> None:
        super().__init__(service_name="efs")

    @property
    def efs_backend(self) -> EFSBackend:
        return efs_backends[self.current_account][self.region]

    def create_file_system(self) -> TYPE_RESPONSE:
        creation_token = self._get_param("CreationToken")
        performance_mode = self._get_param("PerformanceMode")
        encrypted = self._get_param("Encrypted")
        kms_key_id = self._get_param("KmsKeyId")
        throughput_mode = self._get_param("ThroughputMode")
        provisioned_throughput_in_mibps = self._get_param(
            "ProvisionedThroughputInMibps"
        )
        availability_zone_name = self._get_param("AvailabilityZoneName")
        backup = self._get_param("Backup")
        tags = self._get_param("Tags") or []
        resource = self.efs_backend.create_file_system(
            creation_token=creation_token,
            performance_mode=performance_mode,
            encrypted=encrypted,
            kms_key_id=kms_key_id,
            throughput_mode=throughput_mode,
            provisioned_throughput_in_mibps=provisioned_throughput_in_mibps,
            availability_zone_name=availability_zone_name,
            backup=backup,
            tags=tags,
        )
        return (
            json.dumps(resource.info_json()),
            {"status": 201, "Content-Type": "application/json"},
        )

    def describe_file_systems(self) -> TYPE_RESPONSE:
        max_items = self._get_int_param("MaxItems", 10)
        marker = self._get_param("Marker")
        creation_token = self._get_param("CreationToken")
        file_system_id = self._get_param("FileSystemId")
        next_marker, file_systems = self.efs_backend.describe_file_systems(
            marker=marker,
            max_items=max_items,
            creation_token=creation_token,
            file_system_id=file_system_id,
        )
        resp_json: Dict[str, Any] = {
            "FileSystems": [fs.info_json() for fs in file_systems]
        }
        if marker:
            resp_json["Marker"] = marker
        if next_marker:
            resp_json["NextMarker"] = next_marker
        return json.dumps(resp_json), {"Content-Type": "application/json"}

    def create_mount_target(self) -> TYPE_RESPONSE:
        file_system_id = self._get_param("FileSystemId")
        subnet_id = self._get_param("SubnetId")
        ip_address = self._get_param("IpAddress")
        security_groups = self._get_param("SecurityGroups")
        mount_target = self.efs_backend.create_mount_target(
            file_system_id=file_system_id,
            subnet_id=subnet_id,
            ip_address=ip_address,
            security_groups=security_groups,
        )
        return (
            json.dumps(mount_target.info_json()),
            {"Content-Type": "application/json"},
        )

    def describe_mount_targets(self) -> TYPE_RESPONSE:
        max_items = self._get_int_param("MaxItems", 10)
        marker = self._get_param("Marker")
        file_system_id = self._get_param("FileSystemId")
        mount_target_id = self._get_param("MountTargetId")
        access_point_id = self._get_param("AccessPointId")
        next_marker, mount_targets = self.efs_backend.describe_mount_targets(
            max_items=max_items,
            file_system_id=file_system_id,
            mount_target_id=mount_target_id,
            access_point_id=access_point_id,
            marker=marker,
        )
        resp_json: Dict[str, Any] = {
            "MountTargets": [mt.info_json() for mt in mount_targets]
        }
        if marker:
            resp_json["Marker"] = marker
        if next_marker:
            resp_json["NextMarker"] = next_marker
        return json.dumps(resp_json), {"Content-Type": "application/json"}

    def delete_file_system(self) -> TYPE_RESPONSE:
        file_system_id = self._get_param("FileSystemId")
        self.efs_backend.delete_file_system(file_system_id)
        return json.dumps(dict()), {"status": 204, "Content-Type": "application/json"}

    def delete_mount_target(self) -> TYPE_RESPONSE:
        mount_target_id = self._get_param("MountTargetId")
        self.efs_backend.delete_mount_target(mount_target_id)
        return json.dumps(dict()), {"status": 204, "Content-Type": "application/json"}

    def describe_backup_policy(self) -> TYPE_RESPONSE:
        file_system_id = self._get_param("FileSystemId")
        backup_policy = self.efs_backend.describe_backup_policy(file_system_id)
        resp = {"BackupPolicy": backup_policy}
        return json.dumps(resp), {"Content-Type": "application/json"}

    def put_lifecycle_configuration(self) -> TYPE_RESPONSE:
        file_system_id = self._get_param("FileSystemId")
        policies = self._get_param("LifecyclePolicies")
        self.efs_backend.put_lifecycle_configuration(file_system_id, policies)
        return json.dumps({"LifecyclePolicies": policies}), {
            "Content-Type": "application/json"
        }

    def describe_lifecycle_configuration(self) -> TYPE_RESPONSE:
        file_system_id = self._get_param("FileSystemId")
        policies = self.efs_backend.describe_lifecycle_configuration(file_system_id)
        return json.dumps({"LifecyclePolicies": policies}), {
            "Content-Type": "application/json"
        }

    def describe_mount_target_security_groups(self) -> TYPE_RESPONSE:
        mount_target_id = self._get_param("MountTargetId")
        security_groups = self.efs_backend.describe_mount_target_security_groups(
            mount_target_id
        )
        return json.dumps({"SecurityGroups": security_groups}), {
            "Content-Type": "application/json"
        }

    def modify_mount_target_security_groups(self) -> TYPE_RESPONSE:
        mount_target_id = self._get_param("MountTargetId")
        security_groups = self._get_param("SecurityGroups")
        self.efs_backend.modify_mount_target_security_groups(
            mount_target_id, security_groups
        )
        return "{}", {"Content-Type": "application/json"}

    def create_access_point(self) -> TYPE_RESPONSE:
        client_token = self._get_param("ClientToken")
        tags = self._get_param("Tags") or []
        file_system_id = self._get_param("FileSystemId")
        posix_user = self._get_param("PosixUser")
        root_directory = self._get_param("RootDirectory")
        access_point = self.efs_backend.create_access_point(
            client_token,
            tags=tags,
            file_system_id=file_system_id,
            posix_user=posix_user,
            root_directory=root_directory,
        )
        return json.dumps(access_point.info_json()), {
            "Content-Type": "application/json"
        }

    def describe_access_points(self) -> TYPE_RESPONSE:
        access_point_id = self._get_param("AccessPointId")
        file_system_id = self._get_param("FileSystemId")
        access_points = self.efs_backend.describe_access_points(
            access_point_id, file_system_id
        )
        resp = [ap.info_json() for ap in access_points]
        return json.dumps({"AccessPoints": resp}), {"Content-Type": "application/json"}

    def delete_access_point(self) -> TYPE_RESPONSE:
        access_point_id = self._get_param("AccessPointId")
        self.efs_backend.delete_access_point(access_point_id)
        return "{}", {"Content-Type": "application/json"}

    def list_tags_for_resource(self) -> TYPE_RESPONSE:
        resource_id = self._get_param("ResourceId")
        tags = self.efs_backend.list_tags_for_resource(resource_id)
        return json.dumps({"Tags": tags}), {"Content-Type": "application/json"}

    def tag_resource(self) -> TYPE_RESPONSE:
        resource_id = self._get_param("ResourceId")
        tags = self._get_param("Tags")
        self.efs_backend.tag_resource(resource_id, tags)
        return "{}", {"Content-Type": "application/json"}

    def untag_resource(self) -> TYPE_RESPONSE:
        resource_id = self._get_param("ResourceId")
        tag_keys = self.querystring.get("tagKeys", [])
        self.efs_backend.untag_resource(resource_id, tag_keys)
        return "{}", {"Content-Type": "application/json"}
