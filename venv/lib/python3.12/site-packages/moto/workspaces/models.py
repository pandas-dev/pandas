"""WorkSpacesBackend class with methods for supported APIs."""

import re
from collections.abc import Mapping
from typing import Any, Dict, List, Optional, Tuple

from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import BaseModel
from moto.core.utils import unix_time
from moto.ds import ds_backends
from moto.ds.models import Directory
from moto.ec2 import ec2_backends
from moto.moto_api._internal import mock_random
from moto.utilities.utils import get_partition
from moto.workspaces.exceptions import (
    InvalidParameterValuesException,
    ResourceAlreadyExistsException,
    ResourceNotFoundException,
    ValidationException,
)


class Workspace(BaseModel):
    def __init__(
        self,
        workspace: Dict[str, Any],
        running_mode: str,
        error_code: str,
        error_msg: str,
    ):
        self.workspace_properties: Dict[str, Any]
        self.workspace = workspace
        self.workspace_id = f"ws-{mock_random.get_random_hex(9)}"
        # Create_workspaces operation is asynchronous and returns before the WorkSpaces are created.
        # Initially the 'state' is 'PENDING', but here the 'state' will be set as 'AVAILABLE' since
        # this operation is being mocked.
        self.directory_id = workspace["DirectoryId"]
        self.bundle_id = workspace["BundleId"]
        self.user_name = workspace["UserName"]
        self.state = "AVAILABLE"
        self.error_message = error_msg or ""
        self.error_code = error_code or ""
        self.volume_encryption_key = workspace.get("VolumeEncryptionKey", "")
        self.user_volume_encryption_enabled = workspace.get(
            "UserVolumeEncryptionEnabled", ""
        )
        self.root_volume_encryption_enabled = workspace.get(
            "RootVolumeEncryptionEnabled", ""
        )
        workspace_properties = {"RunningMode": running_mode}
        self.workspace_properties = workspace.get("WorkspaceProperties", "")

        if self.workspace_properties:
            self.workspace_properties["RunningMode"] = running_mode
        else:
            self.workspace_properties = workspace_properties

        self.computer_name = ""  # Workspace Bundle
        self.modification_states: List[
            Dict[str, str]
        ] = []  # modify_workspace_properties
        # create_standy_workspace
        self.related_workspaces: List[Dict[str, str]] = []
        self.data_replication_settings: Dict[str, Any] = {}
        # The properties of the standby WorkSpace related to related_workspaces
        self.standby_workspaces_properties: List[Dict[str, Any]] = []
        self.tags = workspace.get("Tags", [])

    def to_dict_pending(self) -> Dict[str, Any]:
        dct = {
            "WorkspaceId": self.workspace_id,
            "DirectoryId": self.directory_id,
            "UserName": self.user_name,
            "IpAddress": "",  # UnKnown
            "State": self.state,
            "BundleId": self.bundle_id,
            "SubnetId": "",  # UnKnown
            "ErrorMessage": self.error_message,
            "ErrorCode": self.error_code,
            "ComputerName": self.computer_name,
            "VolumeEncryptionKey": self.volume_encryption_key,
            "UserVolumeEncryptionEnabled": self.user_volume_encryption_enabled,
            "RootVolumeEncryptionEnabled": self.root_volume_encryption_enabled,
            "WorkspaceProperties": self.workspace_properties,
            "ModificationStates": self.modification_states,
            "RelatedWorkspaces": self.related_workspaces,
            "DataReplicationSettings": self.data_replication_settings,
            "StandbyWorkspacesProperties": self.standby_workspaces_properties,
        }
        return {k: v for k, v in dct.items() if v}

    def filter_empty_values(self, d: Dict[str, Any]) -> Dict[str, Any]:
        if isinstance(d, Mapping):
            return dict((k, self.filter_empty_values(v)) for k, v in d.items() if v)
        else:
            return d

    def to_dict_failed(self) -> Dict[str, Any]:
        dct = {
            "WorkspaceRequest": {
                "DirectoryId": self.workspace["DirectoryId"],
                "UserName": self.workspace["UserName"],
                "BundleId": self.workspace["BundleId"],
                "SubnetId": "",  # UnKnown
                "VolumeEncryptionKey": self.volume_encryption_key,
                "UserVolumeEncryptionEnabled": self.user_volume_encryption_enabled,
                "RootVolumeEncryptionEnabled": self.root_volume_encryption_enabled,
                "WorkspaceProperties": self.workspace_properties,
                "Tags": self.tags,
            },
            "ErrorCode": self.error_code,
            "ErrorMessage": self.error_message,
        }
        return self.filter_empty_values(dct)


class WorkSpaceDirectory(BaseModel):
    def __init__(
        self,
        account_id: str,
        region: str,
        directory: Directory,
        registration_code: str,
        security_group_id: str,
        subnet_ids: List[str],
        enable_work_docs: bool,
        enable_self_service: bool,
        tenancy: str,
        tags: List[Dict[str, str]],
    ):
        self.account_id = account_id
        self.region = region
        self.directory_id = directory.directory_id
        self.alias = directory.alias
        self.directory_name = directory.name
        self.launch_time = directory.launch_time
        self.registration_code = registration_code
        if directory.directory_type == "ADConnector":
            dir_subnet_ids = directory.connect_settings["SubnetIds"]  # type: ignore[index]
        else:
            dir_subnet_ids = directory.vpc_settings["SubnetIds"]  # type: ignore[index]
        self.subnet_ids = subnet_ids or dir_subnet_ids
        self.dns_ip_addresses = directory.dns_ip_addrs
        self.customer_username = "Administrator"
        self.iam_rold_id = (
            f"arn:{get_partition(region)}:iam::{account_id}:role/workspaces_DefaultRole"
        )
        dir_type = directory.directory_type
        if dir_type == "ADConnector":
            self.directory_type = "AD_CONNECTOR"
        elif dir_type == "SimpleAD":
            self.directory_type = "SIMPLE_AD"
        else:
            self.directory_type = dir_type
        self.workspace_security_group_id = security_group_id
        self.state = "REGISTERED"
        # Default values for workspace_creation_properties
        workspace_creation_properties = {
            "EnableWorkDocs": enable_work_docs,
            "EnableInternetAccess": False,
            "DefaultOu": "",
            "CustomSecurityGroupId": "",
            "UserEnabledAsLocalAdministrator": (
                True if self.customer_username == "Administrator" else False
            ),
            "EnableMaintenanceMode": True,
        }
        # modify creation properites
        self.workspace_creation_properties = workspace_creation_properties
        self.ip_group_ids = ""  # create_ip_group
        # Default values for workspace access properties
        workspace_access_properties = {
            "DeviceTypeWindows": (
                "DENY" if self.directory_type == "AD_CONNECTOR" else "ALLOW"
            ),
            "DeviceTypeOsx": "ALLOW",
            "DeviceTypeWeb": "DENY",
            "DeviceTypeIos": "ALLOW",
            "DeviceTypeAndroid": "ALLOW",
            "DeviceTypeChromeOs": "ALLOW",
            "DeviceTypeZeroClient": (
                "DENY" if self.directory_type == "AD_CONNECTOR" else "ALLOW"
            ),
            "DeviceTypeLinux": "DENY",
        }
        # modify_workspace_access_properties
        self.workspace_access_properties = workspace_access_properties
        self.tenancy = tenancy or "SHARED"

        # Default values for self service permissions
        mode = "DISABLED"
        if enable_self_service:
            mode = "ENABLED"
        self_service_permissions = {
            "RestartWorkspace": "ENABLED",
            "IncreaseVolumeSize": mode,
            "ChangeComputeType": mode,
            "SwitchRunningMode": mode,
            "RebuildWorkspace": mode,
        }
        self.self_service_permissions = self_service_permissions
        # Default values for saml properties
        saml_properties = {
            "Status": "DISABLED",
            "UserAccessUrl": "",
            "RelayStateParameterName": "RelayState",
        }
        self.saml_properties = saml_properties
        # Default values for certificate bases auth properties
        self.certificate_based_auth_properties = {
            "Status": "DISABLED",
        }
        # ModifyCertificateBasedAuthProperties
        self.tags = tags or []
        client_properties = {
            # Log uploading is enabled by default.
            "ReconnectEnabled": "ENABLED",
            "LogUploadEnabled": "ENABLED",  # Remember me is enabled by default
        }
        self.client_properties = client_properties

    def delete_security_group(self) -> None:
        """Delete the given security group."""
        ec2_backends[self.account_id][self.region].delete_security_group(
            group_id=self.workspace_security_group_id
        )

    def to_dict(self) -> Dict[str, Any]:
        dct = {
            "DirectoryId": self.directory_id,
            "Alias": self.alias,
            "DirectoryName": self.directory_name,
            "RegistrationCode": self.registration_code,
            "SubnetIds": self.subnet_ids,
            "DnsIpAddresses": self.dns_ip_addresses,
            "CustomerUserName": self.customer_username,
            "IamRoleId": self.iam_rold_id,
            "DirectoryType": self.directory_type,
            "WorkspaceSecurityGroupId": self.workspace_security_group_id,
            "State": self.state,
            "WorkspaceCreationProperties": self.workspace_creation_properties,
            "ipGroupIds": self.ip_group_ids,
            "WorkspaceAccessProperties": self.workspace_access_properties,
            "Tenancy": self.tenancy,
            "SelfservicePermissions": self.self_service_permissions,
            "SamlProperties": self.saml_properties,
            "CertificateBasedAuthProperties": self.certificate_based_auth_properties,
        }
        return {k: v for k, v in dct.items() if v}


class WorkspaceImage(BaseModel):
    def __init__(
        self,
        name: str,
        description: str,
        tags: List[Dict[str, str]],
        account_id: str,
    ):
        self.image_id = f"wsi-{mock_random.get_random_hex(9)}"
        self.name = name
        self.description = description
        self.operating_system: Dict[str, str] = {}  # Unknown
        # Initially the 'state' is 'PENDING', but here the 'state' will be set as 'AVAILABLE' since
        # this operation is being mocked.
        self.state = "AVAILABLE"
        self.required_tenancy = "DEFAULT"
        self.created = unix_time()
        self.owner_account = account_id
        self.error_code = ""
        self.error_message = ""
        self.image_permissions: List[Dict[str, str]] = []
        self.tags = tags

        # Default updates
        self.updates = {
            "UpdateAvailable": False,
            "Description": "This WorkSpace image does not have updates available",
        }
        self.error_details: List[Dict[str, str]] = []

    def to_dict(self) -> Dict[str, Any]:
        dct = {
            "ImageId": self.image_id,
            "Name": self.name,
            "Description": self.description,
            "OperatingSystem": self.operating_system,
            "State": self.state,
            "RequiredTenancy": self.required_tenancy,
            "Created": self.created,
            "OwnerAccountId": self.owner_account,
        }
        return {k: v for k, v in dct.items() if v}

    def to_desc_dict(self) -> Dict[str, Any]:
        dct = self.to_dict()
        dct_options = {
            "ErrorCode": self.error_code,
            "ErrorMessage": self.error_message,
            "Updates": self.updates,
            "ErrorDetails": self.error_details,
        }
        for key, value in dct_options.items():
            if value is not None:
                dct[key] = value
        return dct


class WorkSpacesBackend(BaseBackend):
    """Implementation of WorkSpaces APIs."""

    # The assumption here is that the limits are the same for all regions.
    DIRECTORIES_LIMIT = 50

    def __init__(self, region_name: str, account_id: str):
        super().__init__(region_name, account_id)
        self.workspaces: Dict[str, Workspace] = dict()
        self.workspace_directories: Dict[str, WorkSpaceDirectory] = dict()
        self.workspace_images: Dict[str, WorkspaceImage] = dict()
        self.directories: List[Directory]

    def validate_directory_id(self, value: str, msg: str) -> None:
        """Raise exception if the directory id is invalid."""
        id_pattern = r"^d-[0-9a-f]{10}$"
        if not re.match(id_pattern, value):
            raise ValidationException(msg)

    def validate_image_id(self, value: str, msg: str) -> None:
        """Raise exception if the image id is invalid."""
        id_pattern = r"^wsi-[0-9a-z]{9}$"
        if not re.match(id_pattern, value):
            raise ValidationException(msg)

    def create_security_group(self, directory_id: str, vpc_id: str) -> str:
        """Create security group for the workspace directory."""
        security_group_info = ec2_backends[self.account_id][
            self.region_name
        ].create_security_group(
            name=f"{directory_id}_workspacesMembers",
            description=("Amazon WorkSpaces Security Group"),
            vpc_id=vpc_id,
        )
        return security_group_info.id

    def create_workspaces(
        self, workspaces: List[Dict[str, Any]]
    ) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]]]:
        failed_requests = []
        pending_requests = []

        for ws in workspaces:
            error_code = ""
            error_msg = ""
            directory_id = ws["DirectoryId"]
            msg = f"The Directory ID {directory_id} in the request is invalid."
            self.validate_directory_id(directory_id, msg)

            # FailedRequests are created if the directory_id is unknown
            if directory_id not in self.workspace_directories:
                error_code = "ResourceNotFound.Directory"
                error_msg = "The specified directory could not be found in the specified region."

            running_mode = "ALWAYS_ON"
            workspace_properties = ws.get("WorkspaceProperties", "")
            if workspace_properties:
                running_mode = workspace_properties.get("RunningMode", running_mode)
                auto_stop_timeout = workspace_properties.get(
                    "RunningModeAutoStopTimeoutInMinutes", ""
                )

                # Requests fail if AutoStopTimeout is given for an AlwaysOn Running mode
                if auto_stop_timeout and running_mode == "ALWAYS_ON":
                    error_code = "AutoStopTimeoutIsNotApplicableForAnAlwaysOnWorkspace"
                    error_msg = "RunningModeAutoStopTimeoutInMinutes is not applicable for WorkSpace with running mode set to ALWAYS_ON."

                # Requests fail if AutoStopTimeout is given for an Manual Running mode
                if auto_stop_timeout and running_mode == "MANUAL":
                    error_code = "AutoStopTimeoutIsNotDefaultForManualWorkspace"

            workspace = Workspace(
                workspace=ws,
                running_mode=running_mode,
                error_code=error_code,
                error_msg=error_msg,
            )
            if error_code:
                failed_requests.append(workspace.to_dict_failed())
            else:
                pending_requests.append(workspace.to_dict_pending())
                self.workspaces[workspace.workspace_id] = workspace

        return failed_requests, pending_requests

    def describe_workspaces(
        self,
        workspace_ids: List[str],
        directory_id: str,
        user_name: str,
        bundle_id: str,
    ) -> List[Workspace]:
        # Pagination not yet implemented

        # Only one of the following are allowed to be specified: BundleId, DirectoryId, WorkSpaceIds.
        if (
            (workspace_ids and directory_id)
            or (directory_id and bundle_id)
            or (workspace_ids and bundle_id)
        ):
            msg = "An invalid number of parameters provided with DescribeWorkspaces. Only one of the following are allowed to be specified: BundleId, DirectoryId, WorkSpaceIds, Filters."
            raise InvalidParameterValuesException(msg)

        # Directory_id parameter is required when Username is given.
        if user_name and not directory_id:
            msg = "The DirectoryId parameter is required when UserName is used."
            raise InvalidParameterValuesException(msg)

        workspaces = list(self.workspaces.values())
        if workspace_ids:
            workspaces = [x for x in workspaces if x.workspace_id in workspace_ids]
        if directory_id:
            workspaces = [x for x in workspaces if x.directory_id == directory_id]
        if directory_id and user_name:
            workspaces = [
                x
                for x in workspaces
                if (x.directory_id == directory_id) and (x.user_name == user_name)
            ]
        if bundle_id:
            workspaces = [x for x in workspaces if x.bundle_id == bundle_id]
        # workspaces = [w.to_dict_pending() for w in workspaces]
        return workspaces

    def register_workspace_directory(
        self,
        directory_id: str,
        subnet_ids: List[str],
        enable_work_docs: bool,
        enable_self_service: bool,
        tenancy: str,
        tags: List[Dict[str, str]],
    ) -> None:
        ran_str = mock_random.get_random_string(length=6)
        registration_code = f"SLiad+{ran_str.upper()}"

        (self.directories, _) = ds_backends[self.account_id][
            self.region_name
        ].describe_directories(directory_ids=[directory_id])
        directory = self.directories[0]

        if directory.directory_type == "ADConnector":
            vpc_id = directory.connect_settings["VpcId"]  # type: ignore[index]
        else:
            vpc_id = directory.vpc_settings["VpcId"]  # type: ignore[index]

        security_group_id = self.create_security_group(directory_id, vpc_id)

        workspace_directory = WorkSpaceDirectory(
            account_id=self.account_id,
            region=self.region_name,
            directory=directory,
            registration_code=registration_code,
            security_group_id=security_group_id,
            subnet_ids=subnet_ids,
            enable_work_docs=enable_work_docs,
            enable_self_service=enable_self_service,
            tenancy=tenancy,
            tags=tags,
        )
        self.workspace_directories[workspace_directory.directory_id] = (
            workspace_directory
        )

    def describe_workspace_directories(
        self, directory_ids: Optional[List[str]] = None
    ) -> List[WorkSpaceDirectory]:
        """Return info on all directories or directories with matching IDs."""
        # Pagination not yet implemented

        workspace_directories = list(self.workspace_directories.values())
        if directory_ids:
            for d in directory_ids:
                msg = "The request is invalid."
                self.validate_directory_id(d, msg)
            workspace_directories = [
                x for x in workspace_directories if x.directory_id in directory_ids
            ]

        return sorted(workspace_directories, key=lambda x: x.launch_time)

    def modify_workspace_creation_properties(
        self, resource_id: str, workspace_creation_properties: Dict[str, Any]
    ) -> None:
        # Raise Exception if Directory doesnot exist.
        if resource_id not in self.workspace_directories:
            raise ValidationException("The request is invalid.")

        res = self.workspace_directories[resource_id]
        res.workspace_creation_properties = workspace_creation_properties

    def create_tags(self, resource_id: str, tags: List[Dict[str, str]]) -> None:
        if resource_id.startswith("d-"):
            ds = self.workspace_directories[resource_id]
            ds.tags.extend(tags)
        if resource_id.startswith("ws-"):
            ws = self.workspaces[resource_id]
            ws.tags.extend(tags)

    def describe_tags(self, resource_id: str) -> List[Dict[str, str]]:
        if resource_id.startswith("d-"):
            ds = self.workspace_directories[resource_id]
            tag_list = ds.tags
        if resource_id.startswith("ws-"):
            ws = self.workspaces[resource_id]
            tag_list = ws.tags
        if resource_id.startswith("wsi-"):
            wsi = self.workspace_images[resource_id]
            tag_list = wsi.tags
        return tag_list

    def describe_client_properties(self, resource_ids: str) -> List[Dict[str, Any]]:
        workspace_directories = list(self.workspace_directories.values())
        workspace_directories = [
            x for x in workspace_directories if x.directory_id in resource_ids
        ]
        client_properties_list = []
        for wd in workspace_directories:
            cpl = {
                "ResourceId": wd.directory_id,
                "ClientProperties": wd.client_properties,
            }
            client_properties_list.append(cpl)
        return client_properties_list

    def modify_client_properties(
        self, resource_id: str, client_properties: Dict[str, str]
    ) -> None:
        res = self.workspace_directories[resource_id]
        res.client_properties = client_properties

    def create_workspace_image(
        self, name: str, description: str, workspace_id: str, tags: List[Dict[str, str]]
    ) -> Dict[str, Any]:
        # Check if workspace exists.
        if workspace_id not in self.workspaces:
            raise ResourceNotFoundException(
                "The specified WorkSpace cannot be found. Confirm that the workspace exists in your AWS account, and try again."
            )
        # Check if image name already exists.
        if name in [x.name for x in self.workspace_images.values()]:
            raise ResourceAlreadyExistsException(
                "A WorkSpace image with the same name exists in the destination Region. Provide a unique destination image name, and try again."
            )

        workspace_image = WorkspaceImage(
            name=name,
            description=description,
            tags=tags,
            account_id=self.account_id,
        )
        self.workspace_images[workspace_image.image_id] = workspace_image
        return workspace_image.to_dict()

    def describe_workspace_images(
        self, image_ids: Optional[List[str]], image_type: Optional[str]
    ) -> List[Dict[str, Any]]:
        # Pagination not yet implemented
        workspace_images = list(self.workspace_images.values())
        if image_type == "OWNED":
            workspace_images = [
                i for i in workspace_images if i.owner_account == self.account_id
            ]
        elif image_type == "SHARED":
            workspace_images = [
                i for i in workspace_images if i.owner_account != self.account_id
            ]
        if image_ids:
            workspace_images = [i for i in workspace_images if i.image_id in image_ids]
        return [w.to_desc_dict() for w in workspace_images]

    def update_workspace_image_permission(
        self, image_id: str, allow_copy_image: bool, shared_account_id: str
    ) -> None:
        shared_account = {"SharedAccountId": shared_account_id}
        res = self.workspace_images[image_id]
        shared_accounts = []
        shared_accounts = res.image_permissions

        if shared_account not in shared_accounts and allow_copy_image:
            shared_accounts.append(shared_account)
        if shared_account in shared_accounts and not allow_copy_image:
            shared_accounts.remove(shared_account)

        res.image_permissions = shared_accounts

    def describe_workspace_image_permissions(
        self, image_id: str
    ) -> Tuple[str, List[Dict[str, str]]]:
        # Pagination not yet implemented

        msg = f"The Image ID {image_id} in the request is invalid"
        self.validate_image_id(image_id, msg)

        image_permissions = []
        if image_id in self.workspace_images:
            res = self.workspace_images[image_id]
            image_permissions = res.image_permissions
        return image_id, image_permissions

    def deregister_workspace_directory(self, directory_id: str) -> None:
        """Deregister Workspace Directory with the matching ID."""
        # self._validate_directory_id(directory_id)
        self.workspace_directories[directory_id].delete_security_group()
        self.workspace_directories.pop(directory_id)

    def modify_selfservice_permissions(
        self, resource_id: str, selfservice_permissions: Dict[str, str]
    ) -> None:
        res = self.workspace_directories[resource_id]
        res.self_service_permissions = selfservice_permissions


workspaces_backends = BackendDict(WorkSpacesBackend, "workspaces")
