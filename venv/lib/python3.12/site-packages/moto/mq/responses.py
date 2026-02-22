"""Handles incoming mq requests, invokes methods, returns responses."""

import copy
import json
from urllib.parse import unquote

from moto.core.responses import ActionResult, BaseResponse

from .models import MQBackend, mq_backends


class MQResponse(BaseResponse):
    """Handler for MQ requests and responses."""

    def __init__(self) -> None:
        super().__init__(service_name="mq")

    @property
    def mq_backend(self) -> MQBackend:
        """Return backend instance specific for this region."""
        return mq_backends[self.current_account][self.region]

    def create_broker(self) -> ActionResult:
        params = json.loads(self.body)
        authentication_strategy = params.get("authenticationStrategy")
        auto_minor_version_upgrade = params.get("autoMinorVersionUpgrade")
        broker_name = params.get("brokerName")
        configuration = params.get("configuration")
        deployment_mode = params.get("deploymentMode")
        encryption_options = params.get("encryptionOptions")
        engine_type = params.get("engineType")
        engine_version = params.get("engineVersion")
        host_instance_type = params.get("hostInstanceType")
        ldap_server_metadata = params.get("ldapServerMetadata")
        logs = params.get("logs", {})
        maintenance_window_start_time = params.get("maintenanceWindowStartTime")
        publicly_accessible = params.get("publiclyAccessible")
        security_groups = params.get("securityGroups")
        storage_type = params.get("storageType")
        subnet_ids = params.get("subnetIds", [])
        tags = params.get("tags")
        users = params.get("users", [])
        broker_arn, broker_id = self.mq_backend.create_broker(
            authentication_strategy=authentication_strategy,
            auto_minor_version_upgrade=auto_minor_version_upgrade,
            broker_name=broker_name,
            configuration=configuration,
            deployment_mode=deployment_mode,
            encryption_options=encryption_options,
            engine_type=engine_type,
            engine_version=engine_version,
            host_instance_type=host_instance_type,
            ldap_server_metadata=ldap_server_metadata,
            logs=logs,
            maintenance_window_start_time=maintenance_window_start_time,
            publicly_accessible=publicly_accessible,
            security_groups=security_groups,
            storage_type=storage_type,
            subnet_ids=subnet_ids,
            tags=tags,
            users=users,
        )
        resp = {"brokerArn": broker_arn, "brokerId": broker_id}
        return ActionResult(resp)

    def update_broker(self) -> ActionResult:
        params = json.loads(self.body)
        broker_id = self.path.split("/")[-1]
        authentication_strategy = params.get("authenticationStrategy")
        auto_minor_version_upgrade = params.get("autoMinorVersionUpgrade")
        configuration = params.get("configuration")
        engine_version = params.get("engineVersion")
        host_instance_type = params.get("hostInstanceType")
        ldap_server_metadata = params.get("ldapServerMetadata")
        logs = params.get("logs")
        maintenance_window_start_time = params.get("maintenanceWindowStartTime")
        security_groups = params.get("securityGroups")
        self.mq_backend.update_broker(
            authentication_strategy=authentication_strategy,
            auto_minor_version_upgrade=auto_minor_version_upgrade,
            broker_id=broker_id,
            configuration=configuration,
            engine_version=engine_version,
            host_instance_type=host_instance_type,
            ldap_server_metadata=ldap_server_metadata,
            logs=logs,
            maintenance_window_start_time=maintenance_window_start_time,
            security_groups=security_groups,
        )
        return self.describe_broker()

    def delete_broker(self) -> ActionResult:
        broker_id = self.path.split("/")[-1]
        self.mq_backend.delete_broker(broker_id=broker_id)
        return ActionResult({"BrokerId": broker_id})

    def describe_broker(self) -> ActionResult:
        broker_id = self.path.split("/")[-1]
        broker = self.mq_backend.describe_broker(broker_id=broker_id)
        resp = copy.copy(broker)
        resp.tags = self.mq_backend.list_tags(broker.arn)  # type: ignore[attr-defined]
        return ActionResult(resp)

    def list_brokers(self) -> ActionResult:
        brokers = self.mq_backend.list_brokers()
        return ActionResult({"BrokerSummaries": brokers})

    def create_user(self) -> ActionResult:
        params = json.loads(self.body)
        broker_id = self.path.split("/")[-3]
        username = self.path.split("/")[-1]
        console_access = params.get("consoleAccess", False)
        groups = params.get("groups", [])
        self.mq_backend.create_user(broker_id, username, console_access, groups)
        return ActionResult({})

    def update_user(self) -> ActionResult:
        params = json.loads(self.body)
        broker_id = self.path.split("/")[-3]
        username = self.path.split("/")[-1]
        console_access = params.get("consoleAccess", False)
        groups = params.get("groups", [])
        self.mq_backend.update_user(
            broker_id=broker_id,
            console_access=console_access,
            groups=groups,
            username=username,
        )
        return ActionResult({})

    def describe_user(self) -> ActionResult:
        broker_id = self.path.split("/")[-3]
        username = self.path.split("/")[-1]
        user = self.mq_backend.describe_user(broker_id, username)
        return ActionResult(user)

    def delete_user(self) -> ActionResult:
        broker_id = self.path.split("/")[-3]
        username = self.path.split("/")[-1]
        self.mq_backend.delete_user(broker_id, username)
        return ActionResult({})

    def list_users(self) -> ActionResult:
        broker_id = self.path.split("/")[-2]
        users = self.mq_backend.list_users(broker_id=broker_id)
        resp = {
            "brokerId": broker_id,
            "users": [{"username": u.username} for u in users],
        }
        return ActionResult(resp)

    def create_configuration(self) -> ActionResult:
        params = json.loads(self.body)
        name = params.get("name")
        engine_type = params.get("engineType")
        engine_version = params.get("engineVersion")
        tags = params.get("tags", {})

        config = self.mq_backend.create_configuration(
            name, engine_type, engine_version, tags
        )
        return ActionResult(config)

    def describe_configuration(self) -> ActionResult:
        config_id = self.path.split("/")[-1]
        config = self.mq_backend.describe_configuration(config_id)
        resp = copy.copy(config)
        resp.tags = self.mq_backend.list_tags(config.arn)  # type: ignore[attr-defined]
        return ActionResult(resp)

    def list_configurations(self) -> ActionResult:
        configs = self.mq_backend.list_configurations()
        resp = {"Configurations": configs}
        return ActionResult(resp)

    def update_configuration(self) -> ActionResult:
        config_id = self.path.split("/")[-1]
        params = json.loads(self.body)
        data = params.get("data")
        description = params.get("description")
        config = self.mq_backend.update_configuration(config_id, data, description)
        return ActionResult(config)

    def describe_configuration_revision(self) -> ActionResult:
        revision_id = self.path.split("/")[-1]
        config_id = self.path.split("/")[-3]
        revision = self.mq_backend.describe_configuration_revision(
            config_id, revision_id
        )
        return ActionResult(revision)

    def create_tags(self) -> ActionResult:
        resource_arn = unquote(self.path.split("/")[-1])
        tags = json.loads(self.body).get("tags", {})
        self.mq_backend.create_tags(resource_arn, tags)
        return ActionResult({})

    def delete_tags(self) -> ActionResult:
        resource_arn = unquote(self.path.split("/")[-1])
        tag_keys = self._get_param("tagKeys")
        self.mq_backend.delete_tags(resource_arn, tag_keys)
        return ActionResult({})

    def list_tags(self) -> ActionResult:
        resource_arn = unquote(self.path.split("/")[-1])
        tags = self.mq_backend.list_tags(resource_arn)
        return ActionResult({"Tags": tags})

    def reboot_broker(self) -> ActionResult:
        broker_id = self.path.split("/")[-2]
        self.mq_backend.reboot_broker(broker_id=broker_id)
        return ActionResult({})
