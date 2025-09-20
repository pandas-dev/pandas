"""Handles incoming workspacesweb requests, invokes methods, returns responses."""

import json
from typing import Any
from urllib.parse import unquote

from moto.core.responses import TYPE_RESPONSE, BaseResponse

from .models import WorkSpacesWebBackend, workspacesweb_backends


class WorkSpacesWebResponse(BaseResponse):
    """Handler for WorkSpacesWeb requests and responses."""

    def __init__(self) -> None:
        super().__init__(service_name="workspaces-web")

    @property
    def workspacesweb_backend(self) -> WorkSpacesWebBackend:
        """Return backend instance specific for this region."""
        return workspacesweb_backends[self.current_account][self.region]

    @staticmethod
    def network_settings(request: Any, full_url: str, headers: Any) -> TYPE_RESPONSE:  # type: ignore[misc]
        handler = WorkSpacesWebResponse()
        handler.setup_class(request, full_url, headers)
        if request.method == "GET":
            return handler.get_network_settings()
        else:
            return handler.delete_network_settings()

    @staticmethod
    def browser_settings(request: Any, full_url: str, headers: Any) -> TYPE_RESPONSE:  # type: ignore[misc]
        handler = WorkSpacesWebResponse()
        handler.setup_class(request, full_url, headers)
        if request.method == "GET":
            return handler.get_browser_settings()
        else:
            return handler.delete_browser_settings()

    @staticmethod
    def user_settings(request: Any, full_url: str, headers: Any) -> TYPE_RESPONSE:  # type: ignore[misc]
        handler = WorkSpacesWebResponse()
        handler.setup_class(request, full_url, headers)
        if request.method == "GET":
            return handler.get_user_settings()
        else:
            return handler.delete_user_settings()

    @staticmethod
    def user_access_logging_settings(  # type: ignore[misc]
        request: Any, full_url: str, headers: Any
    ) -> TYPE_RESPONSE:
        handler = WorkSpacesWebResponse()
        handler.setup_class(request, full_url, headers)
        if request.method == "GET":
            return handler.get_user_access_logging_settings()
        else:
            return handler.delete_user_access_logging_settings()

    @staticmethod
    def portal(request: Any, full_url: str, headers: Any) -> TYPE_RESPONSE:  # type: ignore[misc]
        handler = WorkSpacesWebResponse()
        handler.setup_class(request, full_url, headers)
        if request.method == "GET":
            return handler.get_portal()
        else:
            return handler.delete_portal()

    def create_browser_settings(self) -> str:
        additional_encryption_context = self._get_param("additionalEncryptionContext")
        browser_policy = self._get_param("browserPolicy")
        client_token = self._get_param("clientToken")
        customer_managed_key = self._get_param("customerManagedKey")
        tags = self._get_param("tags")
        browser_settings_arn = self.workspacesweb_backend.create_browser_settings(
            additional_encryption_context=additional_encryption_context,
            browser_policy=browser_policy,
            client_token=client_token,
            customer_managed_key=customer_managed_key,
            tags=tags,
        )
        return json.dumps(dict(browserSettingsArn=browser_settings_arn))

    def create_network_settings(self) -> str:
        security_group_ids = self._get_param("securityGroupIds")
        subnet_ids = self._get_param("subnetIds")
        tags = self._get_param("tags")
        vpc_id = self._get_param("vpcId")
        network_settings_arn = self.workspacesweb_backend.create_network_settings(
            security_group_ids=security_group_ids,
            subnet_ids=subnet_ids,
            tags=tags,
            vpc_id=vpc_id,
        )
        return json.dumps(dict(networkSettingsArn=network_settings_arn))

    def get_network_settings(self) -> TYPE_RESPONSE:
        network_settings_arn = unquote(
            self.parsed_url.path.split("/networkSettings/")[-1]
        )
        network_settings = self.workspacesweb_backend.get_network_settings(
            network_settings_arn=network_settings_arn,
        )
        return 200, {}, json.dumps(dict(networkSettings=network_settings))

    def create_portal(self) -> str:
        additional_encryption_context = self._get_param("additionalEncryptionContext")
        authentication_type = self._get_param("authenticationType")
        client_token = self._get_param("clientToken")
        customer_managed_key = self._get_param("customerManagedKey")
        display_name = self._get_param("displayName")
        instance_type = self._get_param("instanceType")
        max_concurrent_sessions = self._get_param("maxConcurrentSessions")
        tags = self._get_param("tags")
        portal_arn, portal_endpoint = self.workspacesweb_backend.create_portal(
            additional_encryption_context=additional_encryption_context,
            authentication_type=authentication_type,
            client_token=client_token,
            customer_managed_key=customer_managed_key,
            display_name=display_name,
            instance_type=instance_type,
            max_concurrent_sessions=max_concurrent_sessions,
            tags=tags,
        )
        return json.dumps(dict(portalArn=portal_arn, portalEndpoint=portal_endpoint))

    def list_browser_settings(self) -> str:
        browser_settings = self.workspacesweb_backend.list_browser_settings()
        return json.dumps(dict(browserSettings=browser_settings))

    def list_network_settings(self) -> str:
        network_settings = self.workspacesweb_backend.list_network_settings()
        return json.dumps(dict(networkSettings=network_settings))

    def list_portals(self) -> str:
        portals = self.workspacesweb_backend.list_portals()
        return json.dumps(dict(portals=portals))

    def get_browser_settings(self) -> TYPE_RESPONSE:
        browser_settings_arn = unquote(
            self.parsed_url.path.split("/browserSettings/")[-1]
        )
        browser_settings = self.workspacesweb_backend.get_browser_settings(
            browser_settings_arn=browser_settings_arn,
        )
        return 200, {}, json.dumps(dict(browserSettings=browser_settings))

    def delete_browser_settings(self) -> TYPE_RESPONSE:
        browser_settings_arn = unquote(
            self.parsed_url.path.split("/browserSettings/")[-1]
        )
        self.workspacesweb_backend.delete_browser_settings(
            browser_settings_arn=browser_settings_arn
        )
        return 200, {}, "{}"

    def delete_network_settings(self) -> TYPE_RESPONSE:
        network_settings_arn = unquote(
            self.parsed_url.path.split("/networkSettings/")[-1]
        )
        self.workspacesweb_backend.delete_network_settings(
            network_settings_arn=network_settings_arn,
        )
        return 200, {}, "{}"

    def get_portal(self) -> TYPE_RESPONSE:
        portal_arn = unquote(self.parsed_url.path.split("/portals/")[-1])
        portal = self.workspacesweb_backend.get_portal(portal_arn=portal_arn)
        return 200, {}, json.dumps(dict(portal=portal))

    def delete_portal(self) -> TYPE_RESPONSE:
        portal_arn = unquote(self.parsed_url.path.split("/portals/")[-1])
        self.workspacesweb_backend.delete_portal(portal_arn=portal_arn)
        return 200, {}, "{}"

    def associate_browser_settings(self) -> str:
        browser_settings_arn = unquote(self._get_param("browserSettingsArn"))
        portal_arn = unquote(
            self.parsed_url.path.split("/portals/")[-1].split("/browserSettings")[0]
        )
        browser_settings_arn, portal_arn = (
            self.workspacesweb_backend.associate_browser_settings(
                browser_settings_arn=browser_settings_arn,
                portal_arn=portal_arn,
            )
        )
        return json.dumps(
            dict(browserSettingsArn=browser_settings_arn, portalArn=portal_arn)
        )

    def associate_network_settings(self) -> str:
        network_settings_arn = unquote(self._get_param("networkSettingsArn"))
        portal_arn = unquote(
            self.parsed_url.path.split("/portals/")[-1].split("/networkSettings")[0]
        )
        network_settings_arn, portal_arn = (
            self.workspacesweb_backend.associate_network_settings(
                network_settings_arn=network_settings_arn,
                portal_arn=portal_arn,
            )
        )
        return json.dumps(
            dict(networkSettingsArn=network_settings_arn, portalArn=portal_arn)
        )

    def create_user_settings(self) -> str:
        additional_encryption_context = self._get_param("additionalEncryptionContext")
        client_token = self._get_param("clientToken")
        cookie_synchronization_configuration = self._get_param(
            "cookieSynchronizationConfiguration"
        )
        copy_allowed = self._get_param("copyAllowed")
        customer_managed_key = self._get_param("customerManagedKey")
        deep_link_allowed = self._get_param("deepLinkAllowed")
        disconnect_timeout_in_minutes = self._get_param("disconnectTimeoutInMinutes")
        download_allowed = self._get_param("downloadAllowed")
        idle_disconnect_timeout_in_minutes = self._get_param(
            "idleDisconnectTimeoutInMinutes"
        )
        paste_allowed = self._get_param("pasteAllowed")
        print_allowed = self._get_param("printAllowed")
        tags = self._get_param("tags")
        upload_allowed = self._get_param("uploadAllowed")
        user_settings_arn = self.workspacesweb_backend.create_user_settings(
            additional_encryption_context=additional_encryption_context,
            client_token=client_token,
            cookie_synchronization_configuration=cookie_synchronization_configuration,
            copy_allowed=copy_allowed,
            customer_managed_key=customer_managed_key,
            deep_link_allowed=deep_link_allowed,
            disconnect_timeout_in_minutes=disconnect_timeout_in_minutes,
            download_allowed=download_allowed,
            idle_disconnect_timeout_in_minutes=idle_disconnect_timeout_in_minutes,
            paste_allowed=paste_allowed,
            print_allowed=print_allowed,
            tags=tags,
            upload_allowed=upload_allowed,
        )
        return json.dumps(dict(userSettingsArn=user_settings_arn))

    def get_user_settings(self) -> TYPE_RESPONSE:
        user_settings_arn = unquote(self.parsed_url.path.split("/userSettings/")[-1])
        user_settings = self.workspacesweb_backend.get_user_settings(
            user_settings_arn=user_settings_arn,
        )
        return 200, {}, json.dumps(dict(userSettings=user_settings))

    def delete_user_settings(self) -> TYPE_RESPONSE:
        user_settings_arn = unquote(self.parsed_url.path.split("/userSettings/")[-1])
        self.workspacesweb_backend.delete_user_settings(
            user_settings_arn=user_settings_arn,
        )
        return 200, {}, "{}"

    def create_user_access_logging_settings(self) -> str:
        params = self._get_params()
        params = json.loads(list(params.keys())[0])
        client_token = params.get("clientToken")
        kinesis_stream_arn = params.get("kinesisStreamArn")
        tags = params.get("tags")
        user_access_logging_settings_arn = (
            self.workspacesweb_backend.create_user_access_logging_settings(
                client_token=client_token,
                kinesis_stream_arn=kinesis_stream_arn,
                tags=tags,
            )
        )
        return json.dumps(
            dict(userAccessLoggingSettingsArn=user_access_logging_settings_arn)
        )

    def get_user_access_logging_settings(self) -> TYPE_RESPONSE:
        user_access_logging_settings_arn = unquote(
            self.parsed_url.path.split("/userAccessLoggingSettings/")[-1]
        )
        user_access_logging_settings = (
            self.workspacesweb_backend.get_user_access_logging_settings(
                user_access_logging_settings_arn=user_access_logging_settings_arn,
            )
        )
        return (
            200,
            {},
            json.dumps(dict(userAccessLoggingSettings=user_access_logging_settings)),
        )

    def delete_user_access_logging_settings(self) -> TYPE_RESPONSE:
        user_access_logging_settings_arn = unquote(
            self.parsed_url.path.split("/userAccessLoggingSettings/")[-1]
        )
        self.workspacesweb_backend.delete_user_access_logging_settings(
            user_access_logging_settings_arn=user_access_logging_settings_arn,
        )
        return 200, {}, "{}"

    def associate_user_settings(self) -> str:
        user_settings_arn = unquote(self._get_param("userSettingsArn"))
        portal_arn = unquote(
            self.parsed_url.path.split("/portals/")[-1].split("/userSettings")[0]
        )
        user_settings_arn, portal_arn = (
            self.workspacesweb_backend.associate_user_settings(
                user_settings_arn=user_settings_arn,
                portal_arn=portal_arn,
            )
        )
        return json.dumps(dict(userSettingsArn=user_settings_arn, portalArn=portal_arn))

    def associate_user_access_logging_settings(self) -> str:
        user_access_logging_settings_arn = unquote(
            self._get_param("userAccessLoggingSettingsArn")
        )
        portal_arn = unquote(
            self.parsed_url.path.split("/portals/")[-1].split(
                "/userAccessLoggingSettings"
            )[0]
        )
        user_access_logging_settings_arn, portal_arn = (
            self.workspacesweb_backend.associate_user_access_logging_settings(
                user_access_logging_settings_arn=user_access_logging_settings_arn,
                portal_arn=portal_arn,
            )
        )
        return json.dumps(
            dict(
                userAccessLoggingSettingsArn=user_access_logging_settings_arn,
                portalArn=portal_arn,
            )
        )

    def list_user_settings(self) -> str:
        user_settings = self.workspacesweb_backend.list_user_settings()
        return json.dumps(dict(userSettings=user_settings))

    def list_user_access_logging_settings(self) -> str:
        user_access_logging_settings = (
            self.workspacesweb_backend.list_user_access_logging_settings()
        )
        return json.dumps(dict(userAccessLoggingSettings=user_access_logging_settings))

    def tag_resource(self) -> str:
        client_token = self._get_param("clientToken")
        resource_arn = unquote(self._get_param("resourceArn"))
        tags = self._get_param("tags")
        self.workspacesweb_backend.tag_resource(
            client_token=client_token,
            resource_arn=resource_arn,
            tags=tags,
        )
        return json.dumps(dict())

    def untag_resource(self) -> str:
        tagKeys = self.__dict__["data"]["tagKeys"]
        resource_arn = unquote(self._get_param("resourceArn"))
        self.workspacesweb_backend.untag_resource(
            resource_arn=resource_arn,
            tag_keys=tagKeys,
        )
        return json.dumps(dict())

    def list_tags_for_resource(self) -> str:
        resource_arn = unquote(self.parsed_url.path.split("/tags/")[-1])
        tags = self.workspacesweb_backend.list_tags_for_resource(
            resource_arn=resource_arn,
        )
        return json.dumps(dict(tags=tags))
