"""Handles Directory Service requests, invokes methods, returns responses."""

import json

from moto.core.exceptions import InvalidToken
from moto.core.responses import BaseResponse
from moto.ds.exceptions import InvalidNextTokenException
from moto.ds.models import DirectoryServiceBackend, ds_backends


class DirectoryServiceResponse(BaseResponse):
    """Handler for DirectoryService requests and responses."""

    def __init__(self) -> None:
        super().__init__(service_name="ds")

    @property
    def ds_backend(self) -> DirectoryServiceBackend:
        """Return backend instance specific for this region."""
        return ds_backends[self.current_account][self.region]

    def connect_directory(self) -> str:
        """Create an AD Connector to connect to a self-managed directory."""
        name = self._get_param("Name")
        short_name = self._get_param("ShortName")
        password = self._get_param("Password")
        description = self._get_param("Description")
        size = self._get_param("Size")
        connect_settings = self._get_param("ConnectSettings")
        tags = self._get_param("Tags", [])
        directory_id = self.ds_backend.connect_directory(
            region=self.region,
            name=name,
            short_name=short_name,
            password=password,
            description=description,
            size=size,
            connect_settings=connect_settings,
            tags=tags,
        )
        return json.dumps({"DirectoryId": directory_id})

    def create_directory(self) -> str:
        """Create a Simple AD directory."""
        name = self._get_param("Name")
        short_name = self._get_param("ShortName")
        password = self._get_param("Password")
        description = self._get_param("Description")
        size = self._get_param("Size")
        vpc_settings = self._get_param("VpcSettings")
        tags = self._get_param("Tags", [])
        directory_id = self.ds_backend.create_directory(
            region=self.region,
            name=name,
            short_name=short_name,
            password=password,
            description=description,
            size=size,
            vpc_settings=vpc_settings,
            tags=tags,
        )
        return json.dumps({"DirectoryId": directory_id})

    def create_alias(self) -> str:
        """Create an alias and assign the alias to the directory."""
        directory_id = self._get_param("DirectoryId")
        alias = self._get_param("Alias")
        response = self.ds_backend.create_alias(directory_id, alias)
        return json.dumps(response)

    def create_microsoft_ad(self) -> str:
        """Create a Microsoft AD directory."""
        name = self._get_param("Name")
        short_name = self._get_param("ShortName")
        password = self._get_param("Password")
        description = self._get_param("Description")
        vpc_settings = self._get_param("VpcSettings")
        edition = self._get_param("Edition")
        tags = self._get_param("Tags", [])
        directory_id = self.ds_backend.create_microsoft_ad(
            region=self.region,
            name=name,
            short_name=short_name,
            password=password,
            description=description,
            vpc_settings=vpc_settings,
            edition=edition,
            tags=tags,
        )
        return json.dumps({"DirectoryId": directory_id})

    def delete_directory(self) -> str:
        """Delete a Directory Service directory."""
        directory_id_arg = self._get_param("DirectoryId")
        directory_id = self.ds_backend.delete_directory(directory_id_arg)
        return json.dumps({"DirectoryId": directory_id})

    def describe_directories(self) -> str:
        """Return directory info for the given IDs or all IDs."""
        directory_ids = self._get_param("DirectoryIds")
        next_token = self._get_param("NextToken")
        limit = self._get_int_param("Limit")
        try:
            (directories, next_token) = self.ds_backend.describe_directories(
                directory_ids, next_token=next_token, limit=limit
            )
        except InvalidToken as exc:
            raise InvalidNextTokenException() from exc

        response = {"DirectoryDescriptions": [x.to_dict() for x in directories]}
        if next_token:
            response["NextToken"] = next_token
        return json.dumps(response)

    def disable_sso(self) -> str:
        """Disable single-sign on for a directory."""
        directory_id = self._get_param("DirectoryId")
        username = self._get_param("UserName")
        password = self._get_param("Password")
        self.ds_backend.disable_sso(directory_id, username, password)
        return ""

    def enable_sso(self) -> str:
        """Enable single-sign on for a directory."""
        directory_id = self._get_param("DirectoryId")
        username = self._get_param("UserName")
        password = self._get_param("Password")
        self.ds_backend.enable_sso(directory_id, username, password)
        return ""

    def get_directory_limits(self) -> str:
        """Return directory limit information for the current region."""
        limits = self.ds_backend.get_directory_limits()
        return json.dumps({"DirectoryLimits": limits})

    def add_tags_to_resource(self) -> str:
        """Add or overwrite on or more tags for specified directory."""
        resource_id = self._get_param("ResourceId")
        tags = self._get_param("Tags")
        self.ds_backend.add_tags_to_resource(resource_id=resource_id, tags=tags)
        return ""

    def remove_tags_from_resource(self) -> str:
        """Removes tags from a directory."""
        resource_id = self._get_param("ResourceId")
        tag_keys = self._get_param("TagKeys")
        self.ds_backend.remove_tags_from_resource(
            resource_id=resource_id, tag_keys=tag_keys
        )
        return ""

    def list_tags_for_resource(self) -> str:
        """Lists all tags on a directory."""
        resource_id = self._get_param("ResourceId")
        next_token = self._get_param("NextToken")
        limit = self._get_param("Limit")
        try:
            tags, next_token = self.ds_backend.list_tags_for_resource(
                resource_id=resource_id, next_token=next_token, limit=limit
            )
        except InvalidToken as exc:
            raise InvalidNextTokenException() from exc

        response = {"Tags": tags}
        if next_token:
            response["NextToken"] = next_token
        return json.dumps(response)
