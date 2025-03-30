"""WorkSpacesWebBackend class with methods for supported APIs."""

import datetime
import uuid
from typing import Any, Dict, List, Optional, Tuple

from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import BaseModel
from moto.utilities.utils import get_partition

from ..utilities.tagging_service import TaggingService


class FakeUserSettings(BaseModel):
    def __init__(
        self,
        additional_encryption_context: Any,
        client_token: str,
        cookie_synchronization_configuration: str,
        copy_allowed: bool,
        customer_managed_key: str,
        deep_link_allowed: bool,
        disconnect_timeout_in_minutes: int,
        download_allowed: bool,
        idle_disconnect_timeout_in_minutes: int,
        paste_allowed: bool,
        print_allowed: bool,
        upload_allowed: bool,
        region_name: str,
        account_id: str,
    ):
        self.user_settings_id = str(uuid.uuid4())
        self.arn = self.arn_formatter(self.user_settings_id, account_id, region_name)
        self.additional_encryption_context = additional_encryption_context
        self.client_token = client_token
        self.cookie_synchronization_configuration = cookie_synchronization_configuration
        self.copy_allowed = copy_allowed if copy_allowed else "Disabled"
        self.customer_managed_key = customer_managed_key
        self.deep_link_allowed = deep_link_allowed if deep_link_allowed else "Disabled"
        self.disconnect_timeout_in_minutes = disconnect_timeout_in_minutes
        self.download_allowed = download_allowed if download_allowed else "Disabled"
        self.idle_disconnect_timeout_in_minutes = idle_disconnect_timeout_in_minutes
        self.paste_allowed = paste_allowed if paste_allowed else "Disabled"
        self.print_allowed = print_allowed if print_allowed else "Disabled"
        self.upload_allowed = upload_allowed if upload_allowed else "Disabled"
        self.associated_portal_arns: List[str] = []

    def arn_formatter(self, _id: str, account_id: str, region_name: str) -> str:
        return f"arn:{get_partition(region_name)}:workspaces-web:{region_name}:{account_id}:user-settings/{_id}"

    def to_dict(self) -> Dict[str, Any]:
        return {
            "associatedPortalArns": self.associated_portal_arns,
            "additionalEncryptionContext": self.additional_encryption_context,
            "clientToken": self.client_token,
            "cookieSynchronizationConfiguration": self.cookie_synchronization_configuration,
            "copyAllowed": self.copy_allowed,
            "customerManagedKey": self.customer_managed_key,
            "deepLinkAllowed": self.deep_link_allowed,
            "disconnectTimeoutInMinutes": self.disconnect_timeout_in_minutes,
            "downloadAllowed": self.download_allowed,
            "idleDisconnectTimeoutInMinutes": self.idle_disconnect_timeout_in_minutes,
            "pasteAllowed": self.paste_allowed,
            "printAllowed": self.print_allowed,
            "uploadAllowed": self.upload_allowed,
            "userSettingsArn": self.arn,
        }


class FakeUserAccessLoggingSettings(BaseModel):
    def __init__(
        self,
        client_token: str,
        kinesis_stream_arn: str,
        region_name: str,
        account_id: str,
    ):
        self.user_access_logging_settings_id = str(uuid.uuid4())
        self.arn = self.arn_formatter(
            self.user_access_logging_settings_id, account_id, region_name
        )
        self.client_token = client_token
        self.kinesis_stream_arn = kinesis_stream_arn
        self.associated_portal_arns: List[str] = []

    def arn_formatter(self, _id: str, account_id: str, region_name: str) -> str:
        return f"arn:{get_partition(region_name)}:workspaces-web:{region_name}:{account_id}:user-access-logging-settings/{_id}"

    def to_dict(self) -> Dict[str, Any]:
        return {
            "associatedPortalArns": self.associated_portal_arns,
            "kinesisStreamArn": self.kinesis_stream_arn,
            "userAccessLoggingSettingsArn": self.arn,
        }


class FakeNetworkSettings(BaseModel):
    def __init__(
        self,
        security_group_ids: List[str],
        subnet_ids: List[str],
        vpc_id: str,
        region_name: str,
        account_id: str,
    ):
        self.network_settings_id = str(uuid.uuid4())
        self.arn = self.arn_formatter(self.network_settings_id, account_id, region_name)
        self.security_group_ids = security_group_ids
        self.subnet_ids = subnet_ids
        self.vpc_id = vpc_id
        self.associated_portal_arns: List[str] = []

    def arn_formatter(self, _id: str, account_id: str, region_name: str) -> str:
        return f"arn:{get_partition(region_name)}:workspaces-web:{region_name}:{account_id}:network-settings/{_id}"

    def to_dict(self) -> Dict[str, Any]:
        return {
            "associatedPortalArns": self.associated_portal_arns,
            "networkSettingsArn": self.arn,
            "securityGroupIds": self.security_group_ids,
            "subnetIds": self.subnet_ids,
            "vpcId": self.vpc_id,
        }


class FakeBrowserSettings(BaseModel):
    def __init__(
        self,
        additional_encryption_context: Any,
        browser_policy: str,
        client_token: str,
        customer_managed_key: str,
        region_name: str,
        account_id: str,
    ):
        self.browser_settings_id = str(uuid.uuid4())
        self.arn = self.arn_formatter(self.browser_settings_id, account_id, region_name)
        self.additional_encryption_context = additional_encryption_context
        self.browser_policy = browser_policy
        self.client_token = client_token
        self.customer_managed_key = customer_managed_key
        self.associated_portal_arns: List[str] = []

    def arn_formatter(self, _id: str, account_id: str, region_name: str) -> str:
        return f"arn:{get_partition(region_name)}:workspaces-web:{region_name}:{account_id}:browser-settings/{_id}"

    def to_dict(self) -> Dict[str, Any]:
        return {
            "associatedPortalArns": self.associated_portal_arns,
            "browserSettingsArn": self.arn,
            "additionalEncryptionContext": self.additional_encryption_context,
            "browserPolicy": self.browser_policy,
            "customerManagedKey": self.customer_managed_key,
        }


class FakePortal(BaseModel):
    def __init__(
        self,
        additional_encryption_context: Any,
        authentication_type: str,
        client_token: str,
        customer_managed_key: str,
        display_name: str,
        instance_type: str,
        max_concurrent_sessions: str,
        region_name: str,
        account_id: str,
    ):
        self.portal_id = str(uuid.uuid4())
        self.arn = self.arn_formatter(self.portal_id, account_id, region_name)
        self.additional_encryption_context = additional_encryption_context
        self.authentication_type = authentication_type
        self.client_token = client_token
        self.customer_managed_key = customer_managed_key
        self.display_name = display_name
        self.instance_type = instance_type
        self.max_concurrent_sessions = max_concurrent_sessions
        self.portal_endpoint = f"{self.portal_id}.portal.aws"
        self.browser_type = "Chrome"
        self.creation_time = datetime.datetime.now().isoformat()
        self.status = "CREATED"
        self.renderer_type = "AppStream"
        self.status_reason = "TestStatusReason"
        self.browser_settings_arn: Optional[str] = None
        self.network_settings_arn: Optional[str] = None
        self.trust_store_arn: Optional[str] = None
        self.ip_access_settings_arn: Optional[str] = None
        self.user_access_logging_settings_arn: Optional[str] = None
        self.user_settings_arn: Optional[str] = None

    def arn_formatter(self, _id: str, account_id: str, region_name: str) -> str:
        return f"arn:{get_partition(region_name)}:workspaces-web:{region_name}:{account_id}:portal/{_id}"

    def to_dict(self) -> Dict[str, Any]:
        return {
            "additionalEncryptionContext": self.additional_encryption_context,
            "authenticationType": self.authentication_type,
            "browserSettingsArn": self.browser_settings_arn,
            "browserType": self.browser_type,
            "creationDate": self.creation_time,
            "customerManagedKey": self.customer_managed_key,
            "displayName": self.display_name,
            "instanceType": self.instance_type,
            "ipAccessSettingsArn": self.ip_access_settings_arn,
            "maxConcurrentSessions": self.max_concurrent_sessions,
            "networkSettingsArn": self.network_settings_arn,
            "portalArn": self.arn,
            "portalEndpoint": self.portal_endpoint,
            "portalStatus": self.status,
            "rendererType": self.renderer_type,
            "statusReason": self.status_reason,
            "trustStoreArn": self.trust_store_arn,
            "userAccessLoggingSettingsArn": self.user_access_logging_settings_arn,
            "userSettingsArn": self.user_settings_arn,
        }


class WorkSpacesWebBackend(BaseBackend):
    """Implementation of WorkSpacesWeb APIs."""

    def __init__(self, region_name: str, account_id: str):
        super().__init__(region_name, account_id)
        self.network_settings: Dict[str, FakeNetworkSettings] = {}
        self.browser_settings: Dict[str, FakeBrowserSettings] = {}
        self.user_settings: Dict[str, FakeUserSettings] = {}
        self.user_access_logging_settings: Dict[str, FakeUserAccessLoggingSettings] = {}
        self.portals: Dict[str, FakePortal] = {}
        self.tagger = TaggingService()

    def create_network_settings(
        self,
        security_group_ids: List[str],
        subnet_ids: List[str],
        tags: Dict[str, str],
        vpc_id: str,
    ) -> str:
        network_settings_object = FakeNetworkSettings(
            security_group_ids,
            subnet_ids,
            vpc_id,
            self.region_name,
            self.account_id,
        )
        self.network_settings[network_settings_object.arn] = network_settings_object
        if tags:
            self.tag_resource("TEMP_CLIENT_TOKEN", network_settings_object.arn, tags)
        return network_settings_object.arn

    def list_network_settings(self) -> List[Dict[str, str]]:
        return [
            {"networkSettingsArn": network_setting.arn, "vpcId": network_setting.vpc_id}
            for network_setting in self.network_settings.values()
        ]

    def get_network_settings(self, network_settings_arn: str) -> Dict[str, Any]:
        return self.network_settings[network_settings_arn].to_dict()

    def delete_network_settings(self, network_settings_arn: str) -> None:
        self.network_settings.pop(network_settings_arn)

    def create_browser_settings(
        self,
        additional_encryption_context: Any,
        browser_policy: str,
        client_token: str,
        customer_managed_key: str,
        tags: Optional[List[Dict[str, str]]] = None,
    ) -> str:
        browser_settings_object = FakeBrowserSettings(
            additional_encryption_context,
            browser_policy,
            client_token,
            customer_managed_key,
            self.region_name,
            self.account_id,
        )
        self.browser_settings[browser_settings_object.arn] = browser_settings_object
        if tags:
            self.tag_resource(client_token, browser_settings_object.arn, tags)
        return browser_settings_object.arn

    def list_browser_settings(self) -> List[Dict[str, str]]:
        return [
            {"browserSettingsArn": browser_setting.arn}
            for browser_setting in self.browser_settings.values()
        ]

    def get_browser_settings(self, browser_settings_arn: str) -> Dict[str, Any]:
        return self.browser_settings[browser_settings_arn].to_dict()

    def delete_browser_settings(self, browser_settings_arn: str) -> None:
        self.browser_settings.pop(browser_settings_arn)

    def create_portal(
        self,
        additional_encryption_context: Any,
        authentication_type: str,
        client_token: str,
        customer_managed_key: str,
        display_name: str,
        instance_type: str,
        max_concurrent_sessions: str,
        tags: Optional[List[Dict[str, str]]] = None,
    ) -> Tuple[str, str]:
        portal_object = FakePortal(
            additional_encryption_context,
            authentication_type,
            client_token,
            customer_managed_key,
            display_name,
            instance_type,
            max_concurrent_sessions,
            self.region_name,
            self.account_id,
        )
        self.portals[portal_object.arn] = portal_object
        if tags:
            self.tag_resource(client_token, portal_object.arn, tags)
        return portal_object.arn, portal_object.portal_endpoint

    def list_portals(self) -> List[Dict[str, Any]]:
        return [
            {
                "authenticationType": portal.authentication_type,
                "browserSettingsArn": portal.browser_settings_arn,
                "browserType": portal.browser_type,
                "creationDate": portal.creation_time,
                "customerManagedKey": portal.customer_managed_key,
                "displayName": portal.display_name,
                "instanceType": portal.instance_type,
                "ipAccessSettingsArn": portal.ip_access_settings_arn,
                "maxConcurrentSessions": portal.max_concurrent_sessions,
                "networkSettingsArn": portal.network_settings_arn,
                "portalArn": portal.arn,
                "portalEndpoint": portal.portal_endpoint,
                "portalStatus": portal.status,
                "rendererType": portal.renderer_type,
                "statusReason": portal.status_reason,
                "trustStoreArn": portal.trust_store_arn,
                "userAccessLoggingSettingsArn": portal.user_access_logging_settings_arn,
                "userSettingsArn": portal.user_settings_arn,
            }
            for portal in self.portals.values()
        ]

    def get_portal(self, portal_arn: str) -> Dict[str, Any]:
        return self.portals[portal_arn].to_dict()

    def delete_portal(self, portal_arn: str) -> None:
        self.portals.pop(portal_arn)

    def associate_browser_settings(
        self, browser_settings_arn: str, portal_arn: str
    ) -> Tuple[str, str]:
        browser_settings_object = self.browser_settings[browser_settings_arn]
        portal_object = self.portals[portal_arn]
        browser_settings_object.associated_portal_arns.append(portal_arn)
        portal_object.browser_settings_arn = browser_settings_arn
        return browser_settings_arn, portal_arn

    def associate_network_settings(
        self, network_settings_arn: str, portal_arn: str
    ) -> Tuple[str, str]:
        network_settings_object = self.network_settings[network_settings_arn]
        portal_object = self.portals[portal_arn]
        network_settings_object.associated_portal_arns.append(portal_arn)
        portal_object.network_settings_arn = network_settings_arn
        return network_settings_arn, portal_arn

    def create_user_settings(
        self,
        additional_encryption_context: Any,
        client_token: Any,
        cookie_synchronization_configuration: Any,
        copy_allowed: Any,
        customer_managed_key: Any,
        deep_link_allowed: Any,
        disconnect_timeout_in_minutes: Any,
        download_allowed: Any,
        idle_disconnect_timeout_in_minutes: Any,
        paste_allowed: Any,
        print_allowed: Any,
        tags: Any,
        upload_allowed: Any,
    ) -> str:
        user_settings_object = FakeUserSettings(
            additional_encryption_context,
            client_token,
            cookie_synchronization_configuration,
            copy_allowed,
            customer_managed_key,
            deep_link_allowed,
            disconnect_timeout_in_minutes,
            download_allowed,
            idle_disconnect_timeout_in_minutes,
            paste_allowed,
            print_allowed,
            upload_allowed,
            self.region_name,
            self.account_id,
        )
        self.user_settings[user_settings_object.arn] = user_settings_object
        if tags:
            self.tag_resource(client_token, user_settings_object.arn, tags)
        return user_settings_object.arn

    def get_user_settings(self, user_settings_arn: str) -> Dict[str, Any]:
        return self.user_settings[user_settings_arn].to_dict()

    def delete_user_settings(self, user_settings_arn: str) -> None:
        self.user_settings.pop(user_settings_arn)

    def create_user_access_logging_settings(
        self, client_token: Any, kinesis_stream_arn: Any, tags: Any
    ) -> str:
        user_access_logging_settings_object = FakeUserAccessLoggingSettings(
            client_token, kinesis_stream_arn, self.region_name, self.account_id
        )
        self.user_access_logging_settings[user_access_logging_settings_object.arn] = (
            user_access_logging_settings_object
        )
        if tags:
            self.tag_resource(
                client_token, user_access_logging_settings_object.arn, tags
            )
        return user_access_logging_settings_object.arn

    def get_user_access_logging_settings(
        self, user_access_logging_settings_arn: str
    ) -> Dict[str, Any]:
        return self.user_access_logging_settings[
            user_access_logging_settings_arn
        ].to_dict()

    def delete_user_access_logging_settings(
        self, user_access_logging_settings_arn: str
    ) -> None:
        self.user_access_logging_settings.pop(user_access_logging_settings_arn)

    def associate_user_settings(
        self, portal_arn: str, user_settings_arn: str
    ) -> Tuple[str, str]:
        user_settings_object = self.user_settings[user_settings_arn]
        portal_object = self.portals[portal_arn]
        user_settings_object.associated_portal_arns.append(portal_arn)
        portal_object.user_settings_arn = user_settings_arn
        return portal_arn, user_settings_arn

    def associate_user_access_logging_settings(
        self, portal_arn: str, user_access_logging_settings_arn: str
    ) -> Tuple[str, str]:
        user_access_logging_settings_object = self.user_access_logging_settings[
            user_access_logging_settings_arn
        ]
        portal_object = self.portals[portal_arn]
        user_access_logging_settings_object.associated_portal_arns.append(portal_arn)
        portal_object.user_access_logging_settings_arn = (
            user_access_logging_settings_arn
        )
        return portal_arn, user_access_logging_settings_arn

    def list_user_settings(self) -> List[Dict[str, str]]:
        return [
            {"userSettingsArn": user_settings.arn}
            for user_settings in self.user_settings.values()
        ]

    def list_user_access_logging_settings(self) -> List[Dict[str, str]]:
        return [
            {"userAccessLoggingSettingsArn": user_access_logging_settings.arn}
            for user_access_logging_settings in self.user_access_logging_settings.values()
        ]

    def tag_resource(self, client_token: str, resource_arn: str, tags: Any) -> None:
        self.tagger.tag_resource(resource_arn, tags)

    def untag_resource(self, resource_arn: str, tag_keys: Any) -> None:
        self.tagger.untag_resource_using_names(resource_arn, tag_keys)

    def list_tags_for_resource(self, resource_arn: str) -> List[Dict[str, str]]:
        tags = self.tagger.get_tag_dict_for_resource(resource_arn)
        Tags = []
        for key, value in tags.items():
            Tags.append({"Key": key, "Value": value})
        return Tags


workspacesweb_backends = BackendDict(WorkSpacesWebBackend, "workspaces-web")
