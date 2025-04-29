"""FSxBackend class with methods for supported APIs."""

from typing import Any, Dict, List, Optional, Tuple
from uuid import uuid4

from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import BaseModel
from moto.utilities.paginator import paginate

from .utils import FileSystemType

PAGINATION_MODEL = {
    "describe_file_systems": {
        "input_token": "next_token",
        "limit_key": "max_results",
        "limit_default": 2147483647,
        "unique_attribute": "resource_arn",
    }
}


class FileSystem(BaseModel):
    def __init__(
        self,
        account_id: str,
        region_name: str,
        file_system_type: str,
        storage_capacity: int,
        storage_type: str,
        subnet_ids: List[str],
        security_group_ids: List[str],
        tags: Optional[List[Dict[str, str]]],
        kms_key_id: Optional[str],
        windows_configuration: Optional[Dict[str, Any]],
        lustre_configuration: Optional[Dict[str, Any]],
        ontap_configuration: Optional[Dict[str, Any]],
        open_zfs_configuration: Optional[Dict[str, Any]],
    ) -> None:
        self.file_system_id = f"fs-{uuid4().hex[:8]}"
        self.file_system_type = file_system_type
        if self.file_system_type not in FileSystemType.list_values():
            raise ValueError(f"Invalid FileSystemType: {self.file_system_type}")
        self.storage_capacity = storage_capacity
        self.storage_type = storage_type
        self.subnet_ids = subnet_ids
        self.security_group_ids = security_group_ids
        self.dns_name = f"{self.file_system_id}.fsx.{region_name}.amazonaws.com"
        self.kms_key_id = kms_key_id
        self.resource_arn = (
            f"arn:aws:fsx:{region_name}:{account_id}:file-system/{self.file_system_id}"
        )
        self.tags = tags or []
        self.windows_configuration = windows_configuration
        self.lustre_configuration = lustre_configuration
        self.ontap_configuration = ontap_configuration
        self.open_zfs_configuration = open_zfs_configuration

    def to_dict(self) -> Dict[str, Any]:
        dct = {
            "FileSystemId": self.file_system_id,
            "FileSystemType": self.file_system_type,
            "StorageCapacity": self.storage_capacity,
            "StorageType": self.storage_type,
            "SubnetIds": self.subnet_ids,
            "SecurityGroupIds": self.security_group_ids,
            "Tags": self.tags,
            "DNSName": self.dns_name,
            "KmsKeyId": self.kms_key_id,
            "ResourceARN": self.resource_arn,
            "WindowsConfiguration": self.windows_configuration,
            "LustreConfiguration": self.lustre_configuration,
            "OntapConfiguration": self.ontap_configuration,
            "OpenZFSConfiguration": self.open_zfs_configuration,
        }
        return {k: v for k, v in dct.items() if v}


class FSxBackend(BaseBackend):
    """Implementation of FSx APIs."""

    def __init__(self, region_name: str, account_id: str) -> None:
        super().__init__(region_name, account_id)
        self.file_systems: Dict[str, FileSystem] = {}

    def create_file_system(
        self,
        client_request_token: str,
        file_system_type: str,
        storage_capacity: int,
        storage_type: str,
        subnet_ids: List[str],
        security_group_ids: List[str],
        tags: Optional[List[Dict[str, str]]],
        kms_key_id: Optional[str],
        windows_configuration: Optional[Dict[str, Any]],
        lustre_configuration: Optional[Dict[str, Any]],
        ontap_configuration: Optional[Dict[str, Any]],
        file_system_type_version: Optional[str],
        open_zfs_configuration: Optional[Dict[str, Any]],
    ) -> FileSystem:
        file_system = FileSystem(
            account_id=self.account_id,
            region_name=self.region_name,
            file_system_type=file_system_type,
            storage_capacity=storage_capacity,
            storage_type=storage_type,
            subnet_ids=subnet_ids,
            security_group_ids=security_group_ids,
            tags=tags,
            kms_key_id=kms_key_id,
            windows_configuration=windows_configuration,
            ontap_configuration=ontap_configuration,
            open_zfs_configuration=open_zfs_configuration,
            lustre_configuration=lustre_configuration,
        )

        file_system_id = file_system.file_system_id

        self.file_systems[file_system_id] = file_system
        return file_system

    @paginate(pagination_model=PAGINATION_MODEL)
    def describe_file_systems(self, file_system_ids: List[str]) -> List[FileSystem]:
        file_systems = []
        if not file_system_ids:
            file_systems = list(self.file_systems.values())
        else:
            for id in file_system_ids:
                if id in self.file_systems:
                    file_systems.append(self.file_systems[id])
        return file_systems

    def delete_file_system(
        self,
        file_system_id: str,
        client_request_token: str,
        windows_configuration: Optional[Dict[str, Any]],
        lustre_configuration: Optional[Dict[str, Any]],
        open_zfs_configuration: Optional[Dict[str, Any]],
    ) -> Tuple[
        str,
        str,
        Optional[Dict[str, Any]],
        Optional[Dict[str, Any]],
        Optional[Dict[str, Any]],
    ]:
        response_template = {"FinalBackUpId": "", "FinalBackUpTags": []}

        file_system_type = self.file_systems[file_system_id].file_system_type

        lifecycle = "DELETING"
        self.file_systems.pop(file_system_id)

        windows_response = None
        lustre_response = None
        open_zfs_response = None

        if file_system_type == "WINDOWS":
            windows_response = response_template
        elif file_system_type == "LUSTRE":
            lustre_response = response_template
        elif file_system_type == "OPEN_ZFS":
            open_zfs_response = response_template

        return (
            file_system_id,
            lifecycle,
            windows_response,
            lustre_response,
            open_zfs_response,
        )

    def tag_resource(self, resource_arn: str, tags: List[Dict[str, str]]) -> None:
        resource = self._get_resource_from_arn(resource_arn)
        resource.tags.extend(tags)

    def _get_resource_from_arn(self, arn: str) -> FileSystem:
        target_resource, target_name = arn.split(":")[-1].split("/")
        try:
            return self.file_systems[target_name]
        except KeyError:
            message = f"Could not find {target_resource} with name {target_name}"
            raise ValueError(message)

    def untag_resource(self, resource_arn: str, tag_keys: List[str]) -> None:
        resource = self._get_resource_from_arn(resource_arn)
        if tag_keys:
            resource.tags = [tag for tag in resource.tags if tag["Key"] not in tag_keys]


fsx_backends = BackendDict(FSxBackend, "fsx")
