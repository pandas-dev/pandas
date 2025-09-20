"""TransferBackend class with methods for supported APIs."""

from typing import Any, Dict, List, Optional, Tuple

from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.utils import unix_time
from moto.transfer.exceptions import PublicKeyNotFound, ServerNotFound, UserNotFound

from .types import (
    Server,
    ServerDomain,
    ServerEndpointType,
    ServerIdentityProviderType,
    ServerProtocols,
    User,
    UserHomeDirectoryType,
)


class TransferBackend(BaseBackend):
    """Implementation of Transfer APIs."""

    def __init__(self, region_name: str, account_id: str) -> None:
        super().__init__(region_name, account_id)
        self.servers: Dict[str, Server] = {}

    def create_server(
        self,
        certificate: Optional[str],
        domain: Optional[ServerDomain],
        endpoint_details: Optional[Dict[str, Any]],
        endpoint_type: Optional[ServerEndpointType],
        host_key: str,
        identity_provider_details: Optional[Dict[str, Any]],
        identity_provider_type: Optional[ServerIdentityProviderType],
        logging_role: Optional[str],
        post_authentication_login_banner: Optional[str],
        pre_authentication_login_banner: Optional[str],
        protocols: Optional[List[ServerProtocols]],
        protocol_details: Optional[Dict[str, Any]],
        security_policy_name: Optional[str],
        tags: Optional[List[Dict[str, str]]],
        workflow_details: Optional[Dict[str, Any]],
        structured_log_destinations: Optional[List[str]],
        s3_storage_options: Optional[Dict[str, Optional[str]]],
    ) -> str:
        server = Server(
            certificate=certificate,
            domain=domain,
            endpoint_type=endpoint_type,
            host_key_fingerprint=host_key,
            identity_provider_type=identity_provider_type,
            logging_role=logging_role,
            post_authentication_login_banner=post_authentication_login_banner,
            pre_authentication_login_banner=pre_authentication_login_banner,
            protocols=protocols,
            security_policy_name=security_policy_name,
            structured_log_destinations=structured_log_destinations,
            tags=(tags or []),
        )
        if endpoint_details is not None:
            endpoint_details = {
                "address_allocation_ids": endpoint_details.get("AddressAllocationIds"),
                "subnet_ids": endpoint_details.get("SubnetIds"),
                "vpc_endpoint_id": endpoint_details.get("VpcEndpointId"),
                "vpc_id": endpoint_details.get("VpcId"),
                "security_group_ids": endpoint_details.get("SecurityGroupIds"),
            }
            server.endpoint_details = endpoint_details
        if identity_provider_details is not None:
            identity_provider_details = {
                "url": identity_provider_details.get("Url"),
                "invocation_role": identity_provider_details.get("InvocationRole"),
                "directory_id": identity_provider_details.get("DirectoryId"),
                "function": identity_provider_details.get("Function"),
                "sftp_authentication_methods": identity_provider_details.get(
                    "SftpAuthenticationMethods"
                ),
            }
            server.identity_provider_details = identity_provider_details
        if protocol_details is not None:
            protocol_details = {
                "passive_ip": protocol_details.get("PassiveIp"),
                "tls_session_resumption_mode": protocol_details.get(
                    "TlsSessionResumptionMode"
                ),
                "set_stat_option": protocol_details.get("SetStatOption"),
                "as2_transports": protocol_details.get("As2Transports"),
            }
            server.protocol_details = protocol_details
        if s3_storage_options is not None:
            server.s3_storage_options = {
                "directory_listing_optimization": s3_storage_options.get(
                    "DirectoryListingOptimization"
                )
            }
        if workflow_details is not None:
            server.workflow_details = {
                "on_upload": [
                    {
                        "workflow_id": workflow.get("WorkflowId"),
                        "execution_role": workflow.get("ExecutionRole"),
                    }
                    for workflow in (workflow_details.get("OnUpload") or [])
                ],
                "on_partial_upload": [
                    {
                        "workflow_id": workflow.get("WorkflowId"),
                        "execution_role": workflow.get("ExecutionRole"),
                    }
                    for workflow in (workflow_details.get("OnPartialUpload") or [])
                ],
            }
        server_id = server.server_id
        self.servers[server_id] = server
        return server_id

    def describe_server(self, server_id: str) -> Server:
        if server_id not in self.servers:
            ServerNotFound(server_id=server_id)
        server = self.servers[server_id]
        return server

    def delete_server(self, server_id: str) -> None:
        if server_id not in self.servers:
            ServerNotFound(server_id=server_id)
        del self.servers[server_id]
        return

    def create_user(
        self,
        home_directory: Optional[str],
        home_directory_type: Optional[UserHomeDirectoryType],
        home_directory_mappings: Optional[List[Dict[str, Optional[str]]]],
        policy: Optional[str],
        posix_profile: Optional[Dict[str, Any]],
        role: str,
        server_id: str,
        ssh_public_key_body: Optional[str],
        tags: Optional[List[Dict[str, str]]],
        user_name: str,
    ) -> Tuple[str, str]:
        if server_id not in self.servers:
            ServerNotFound(server_id=server_id)
        user = User(
            home_directory=home_directory,
            home_directory_type=home_directory_type,
            policy=policy,
            role=role,
            tags=(tags or []),
            user_name=user_name,
        )
        if home_directory_mappings:
            for mapping in home_directory_mappings:
                user.home_directory_mappings.append(
                    {
                        "entry": mapping.get("Entry"),
                        "target": mapping.get("Target"),
                        "type": mapping.get("Type"),
                    }
                )
        if posix_profile is not None:
            posix_profile = {
                "gid": posix_profile.get("Gid"),
                "uid": posix_profile.get("Uid"),
                "secondary_gids": posix_profile.get("SecondaryGids"),
            }
            user.posix_profile = posix_profile
        if ssh_public_key_body is not None:
            now = unix_time()
            ssh_public_keys = [
                {
                    "date_imported": str(now),
                    "ssh_public_key_body": ssh_public_key_body,
                    "ssh_public_key_id": "mock_ssh_public_key_id_{ssh_public_key_body}_{now}",
                }
            ]
            user.ssh_public_keys = ssh_public_keys
        self.servers[server_id]._users.append(user)
        self.servers[server_id].user_count += 1
        return server_id, user_name

    def describe_user(self, server_id: str, user_name: str) -> Tuple[str, User]:
        if server_id not in self.servers:
            raise ServerNotFound(server_id=server_id)
        for user in self.servers[server_id]._users:
            if user.user_name == user_name:
                return server_id, user
        raise UserNotFound(user_name=user_name, server_id=server_id)

    def delete_user(self, server_id: str, user_name: str) -> None:
        if server_id not in self.servers:
            raise ServerNotFound(server_id=server_id)
        for i, user in enumerate(self.servers[server_id]._users):
            if user.user_name == user_name:
                del self.servers[server_id]._users[i]
                self.servers[server_id].user_count -= 1
                return
        raise UserNotFound(server_id=server_id, user_name=user_name)

    def import_ssh_public_key(
        self, server_id: str, ssh_public_key_body: str, user_name: str
    ) -> Tuple[str, str, str]:
        if server_id not in self.servers:
            raise ServerNotFound(server_id=server_id)
        for user in self.servers[server_id]._users:
            if user.user_name == user_name:
                date_imported = unix_time()
                ssh_public_key_id = (
                    f"{server_id}:{user_name}:public_key:{date_imported}"
                )
                key = {
                    "ssh_public_key_id": ssh_public_key_id,
                    "ssh_public_key_body": ssh_public_key_body,
                    "date_imported": str(date_imported),
                }
                user.ssh_public_keys.append(key)
                return server_id, ssh_public_key_id, user_name
        raise UserNotFound(user_name=user_name, server_id=server_id)

    def delete_ssh_public_key(
        self, server_id: str, ssh_public_key_id: str, user_name: str
    ) -> None:
        if server_id not in self.servers:
            raise ServerNotFound(server_id=server_id)
        for i, user in enumerate(self.servers[server_id]._users):
            if user.user_name == user_name:
                for j, key in enumerate(
                    self.servers[server_id]._users[i].ssh_public_keys
                ):
                    if key["ssh_public_key_id"] == ssh_public_key_id:
                        del user.ssh_public_keys[j]
                        return
                raise PublicKeyNotFound(
                    user_name=user_name,
                    server_id=server_id,
                    ssh_public_key_id=ssh_public_key_id,
                )
        raise UserNotFound(user_name=user_name, server_id=server_id)


transfer_backends = BackendDict(TransferBackend, "transfer")
