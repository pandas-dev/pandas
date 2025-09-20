"""Handles incoming workspaces requests, invokes methods, returns responses."""

import json

from moto.core.responses import BaseResponse

from .models import WorkSpacesBackend, workspaces_backends


class WorkSpacesResponse(BaseResponse):
    """Handler for WorkSpaces requests and responses."""

    def __init__(self) -> None:
        super().__init__(service_name="workspaces")

    @property
    def workspaces_backend(self) -> WorkSpacesBackend:
        """Return backend instance specific for this region."""
        return workspaces_backends[self.current_account][self.region]

    def create_workspaces(self) -> str:
        params = json.loads(self.body)
        workspaces = params.get("Workspaces")
        failed_requests, pending_requests = self.workspaces_backend.create_workspaces(
            workspaces=workspaces,
        )
        return json.dumps(
            dict(FailedRequests=failed_requests, PendingRequests=pending_requests)
        )

    def describe_workspaces(self) -> str:
        params = json.loads(self.body)
        workspace_ids = params.get("WorkspaceIds")
        directory_id = params.get("DirectoryId")
        user_name = params.get("UserName")
        bundle_id = params.get("BundleId")
        workspaces = self.workspaces_backend.describe_workspaces(
            workspace_ids=workspace_ids,
            directory_id=directory_id,
            user_name=user_name,
            bundle_id=bundle_id,
        )
        return json.dumps(dict(Workspaces=[x.to_dict_pending() for x in workspaces]))

    def describe_workspace_directories(self) -> str:
        params = json.loads(self.body)
        directory_ids = params.get("DirectoryIds")
        directories = self.workspaces_backend.describe_workspace_directories(
            directory_ids=directory_ids,
        )
        return json.dumps(dict(Directories=[d.to_dict() for d in directories]))

    def register_workspace_directory(self) -> str:
        params = json.loads(self.body)
        directory_id = params.get("DirectoryId")
        subnet_ids = params.get("SubnetIds")
        enable_work_docs = params.get("EnableWorkDocs")
        enable_self_service = params.get("EnableSelfService")
        tenancy = params.get("Tenancy")
        tags = params.get("Tags")
        self.workspaces_backend.register_workspace_directory(
            directory_id=directory_id,
            subnet_ids=subnet_ids,
            enable_work_docs=enable_work_docs,
            enable_self_service=enable_self_service,
            tenancy=tenancy,
            tags=tags,
        )
        return json.dumps(dict())

    def modify_workspace_creation_properties(self) -> str:
        params = json.loads(self.body)
        resource_id = params.get("ResourceId")
        workspace_creation_properties = params.get("WorkspaceCreationProperties")
        self.workspaces_backend.modify_workspace_creation_properties(
            resource_id=resource_id,
            workspace_creation_properties=workspace_creation_properties,
        )
        return "{}"

    def create_tags(self) -> str:
        params = json.loads(self.body)
        resource_id = params.get("ResourceId")
        tags = params.get("Tags")
        self.workspaces_backend.create_tags(
            resource_id=resource_id,
            tags=tags,
        )
        return "{}"

    def describe_tags(self) -> str:
        params = json.loads(self.body)
        resource_id = params.get("ResourceId")
        tag_list = self.workspaces_backend.describe_tags(
            resource_id=resource_id,
        )
        return json.dumps(dict(TagList=tag_list))

    def describe_client_properties(self) -> str:
        params = json.loads(self.body)
        resource_ids = params.get("ResourceIds")
        client_properties_list = self.workspaces_backend.describe_client_properties(
            resource_ids=resource_ids,
        )
        return json.dumps(dict(ClientPropertiesList=client_properties_list))

    def modify_client_properties(self) -> str:
        params = json.loads(self.body)
        resource_id = params.get("ResourceId")
        client_properties = params.get("ClientProperties")
        self.workspaces_backend.modify_client_properties(
            resource_id=resource_id,
            client_properties=client_properties,
        )
        return "{}"

    def create_workspace_image(self) -> str:
        params = json.loads(self.body)
        name = params.get("Name")
        description = params.get("Description")
        workspace_id = params.get("WorkspaceId")
        tags = params.get("Tags")
        workspace_image = self.workspaces_backend.create_workspace_image(
            name=name,
            description=description,
            workspace_id=workspace_id,
            tags=tags,
        )
        return json.dumps(workspace_image)

    def describe_workspace_images(self) -> str:
        params = json.loads(self.body)
        image_ids = params.get("ImageIds")
        image_type = params.get("ImageType")
        images = self.workspaces_backend.describe_workspace_images(
            image_ids=image_ids,
            image_type=image_type,
        )
        return json.dumps(dict(Images=images))

    def update_workspace_image_permission(self) -> str:
        params = json.loads(self.body)
        image_id = params.get("ImageId")
        allow_copy_image = params.get("AllowCopyImage")
        shared_account_id = params.get("SharedAccountId")
        self.workspaces_backend.update_workspace_image_permission(
            image_id=image_id,
            allow_copy_image=allow_copy_image,
            shared_account_id=shared_account_id,
        )
        return "{}"

    def describe_workspace_image_permissions(self) -> str:
        params = json.loads(self.body)
        image_id = params.get("ImageId")
        (
            image_id,
            image_permissions,
        ) = self.workspaces_backend.describe_workspace_image_permissions(
            image_id=image_id,
        )
        return json.dumps(dict(ImageId=image_id, ImagePermissions=image_permissions))

    def deregister_workspace_directory(self) -> str:
        params = json.loads(self.body)
        directory_id = params.get("DirectoryId")
        self.workspaces_backend.deregister_workspace_directory(
            directory_id=directory_id,
        )
        return "{}"

    def modify_selfservice_permissions(self) -> str:
        params = json.loads(self.body)
        resource_id = params.get("ResourceId")
        selfservice_permissions = params.get("SelfservicePermissions")
        self.workspaces_backend.modify_selfservice_permissions(
            resource_id=resource_id,
            selfservice_permissions=selfservice_permissions,
        )
        return "{}"
