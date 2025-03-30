"""Handles incoming fsx requests, invokes methods, returns responses."""

import json

from moto.core.common_types import TYPE_RESPONSE
from moto.core.responses import BaseResponse

from .models import FSxBackend, fsx_backends


class FSxResponse(BaseResponse):
    """Handler for FSx requests and responses."""

    def __init__(self) -> None:
        super().__init__(service_name="fsx")

    @property
    def fsx_backend(self) -> FSxBackend:
        """Return backend instance specific for this region."""
        return fsx_backends[self.current_account][self.region]

    def create_file_system(self) -> str:
        params = json.loads(self.body)
        client_request_token = params.get("ClientRequestToken")
        file_system_type = params.get("FileSystemType")
        storage_capacity = params.get("StorageCapacity")
        storage_type = params.get("StorageType")
        subnet_ids = params.get("SubnetIds")
        security_group_ids = params.get("SecurityGroupIds")
        tags = params.get("Tags")
        kms_key_id = params.get("KmsKeyId")
        windows_configuration = params.get("WindowsConfiguration")
        lustre_configuration = params.get("LustreConfiguration")
        ontap_configuration = params.get("OntapConfiguration")
        file_system_type_version = params.get("FileSystemTypeVersion")
        open_zfs_configuration = params.get("OpenZFSConfiguration")
        file_system = self.fsx_backend.create_file_system(
            client_request_token=client_request_token,
            file_system_type=file_system_type,
            storage_capacity=storage_capacity,
            storage_type=storage_type,
            subnet_ids=subnet_ids,
            security_group_ids=security_group_ids,
            tags=tags,
            kms_key_id=kms_key_id,
            windows_configuration=windows_configuration,
            lustre_configuration=lustre_configuration,
            ontap_configuration=ontap_configuration,
            file_system_type_version=file_system_type_version,
            open_zfs_configuration=open_zfs_configuration,
        )

        return json.dumps(dict(FileSystem=file_system.to_dict()))

    def describe_file_systems(self) -> str:
        params = json.loads(self.body)
        file_system_ids = params.get("FileSystemIds")
        max_results = params.get("MaxResults")
        next_token = params.get("NextToken")
        file_systems, next_token = self.fsx_backend.describe_file_systems(
            file_system_ids=file_system_ids,
            max_results=max_results,
            next_token=next_token,
        )
        list_file_systems = [file_system.to_dict() for file_system in file_systems]
        return json.dumps(dict(FileSystems=list_file_systems, NextToken=next_token))

    def delete_file_system(self) -> str:
        params = json.loads(self.body)
        file_system_id = params.get("FileSystemId")
        client_request_token = params.get("ClientRequestToken")
        windows_configuration = params.get("WindowsConfiguration")
        lustre_configuration = params.get("LustreConfiguration")
        open_zfs_configuration = params.get("OpenZFSConfiguration")
        (
            file_system_id,
            lifecycle,
            windows_response,
            lustre_response,
            open_zfs_response,
        ) = self.fsx_backend.delete_file_system(
            file_system_id=file_system_id,
            client_request_token=client_request_token,
            windows_configuration=windows_configuration,
            lustre_configuration=lustre_configuration,
            open_zfs_configuration=open_zfs_configuration,
        )

        return json.dumps(
            dict(
                FileSystemId=file_system_id,
                Lifecycle=lifecycle,
                WindowsResponse=windows_response,
                LustreResponse=lustre_response,
                OpenZfsResponse=open_zfs_response,
            )
        )

    def tag_resource(self) -> TYPE_RESPONSE:
        params = json.loads(self.body)
        resource_arn = params.get("ResourceARN")
        tags = params.get("Tags")
        self.fsx_backend.tag_resource(
            resource_arn=resource_arn,
            tags=tags,
        )
        return 200, {}, json.dumps({})

    def untag_resource(self) -> TYPE_RESPONSE:
        params = json.loads(self.body)
        resource_arn = params.get("ResourceARN")
        tag_keys = params.get("TagKeys")
        self.fsx_backend.untag_resource(
            resource_arn=resource_arn,
            tag_keys=tag_keys,
        )
        return 200, {}, json.dumps({})
