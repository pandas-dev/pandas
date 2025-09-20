"""Handles incoming transfer requests, invokes methods, returns responses."""

import json

from moto.core.responses import BaseResponse

from .models import TransferBackend, transfer_backends


class TransferResponse(BaseResponse):
    """Handler for Transfer requests and responses."""

    def __init__(self) -> None:
        super().__init__(service_name="transfer")

    @property
    def transfer_backend(self) -> TransferBackend:
        return transfer_backends[self.current_account][self.region]

    def create_user(self) -> str:
        params = json.loads(self.body)
        server_id, user_name = self.transfer_backend.create_user(
            home_directory=params.get("HomeDirectory"),
            home_directory_type=params.get("HomeDirectoryType"),
            home_directory_mappings=params.get("HomeDirectoryMappings"),
            policy=params.get("Policy"),
            posix_profile=params.get("PosixProfile"),
            role=params.get("Role"),
            server_id=params.get("ServerId"),
            ssh_public_key_body=params.get("SshPublicKeyBody"),
            tags=params.get("Tags"),
            user_name=params.get("UserName"),
        )
        return json.dumps(dict(ServerId=server_id, UserName=user_name))

    def describe_user(self) -> str:
        params = json.loads(self.body)
        server_id, user = self.transfer_backend.describe_user(
            server_id=params.get("ServerId"),
            user_name=params.get("UserName"),
        )
        return json.dumps(dict(ServerId=server_id, User=user.to_dict()))

    def delete_user(self) -> str:
        params = json.loads(self.body)
        self.transfer_backend.delete_user(
            server_id=params.get("ServerId"),
            user_name=params.get("UserName"),
        )
        return json.dumps(dict())

    def import_ssh_public_key(self) -> str:
        params = json.loads(self.body)
        server_id, ssh_public_key_id, user_name = (
            self.transfer_backend.import_ssh_public_key(
                server_id=params.get("ServerId"),
                ssh_public_key_body=params.get("SshPublicKeyBody"),
                user_name=params.get("UserName"),
            )
        )
        return json.dumps(
            dict(
                ServerId=server_id, SshPublicKeyId=ssh_public_key_id, UserName=user_name
            )
        )

    def delete_ssh_public_key(self) -> str:
        params = json.loads(self.body)
        self.transfer_backend.delete_ssh_public_key(
            server_id=params.get("ServerId"),
            ssh_public_key_id=params.get("SshPublicKeyId"),
            user_name=params.get("UserName"),
        )
        return json.dumps(dict())

    def create_server(self) -> str:
        params = json.loads(self.body)
        server_id = self.transfer_backend.create_server(
            certificate=params.get("Certificate"),
            domain=params.get("Domain"),
            endpoint_details=params.get("EndpointDetails"),
            endpoint_type=params.get("EndpointType"),
            host_key=params.get("HostKey"),
            identity_provider_details=params.get("IdentityProviderDetails"),
            identity_provider_type=params.get("IdentityProviderType"),
            logging_role=params.get("LoggingRole"),
            post_authentication_login_banner=params.get(
                "PostAuthenticationLoginBanner"
            ),
            pre_authentication_login_banner=params.get("PreAuthenticationLoginBanner"),
            protocols=params.get("Protocols"),
            protocol_details=params.get("ProtocolDetails"),
            security_policy_name=params.get("SecurityPolicyName"),
            structured_log_destinations=params.get("StructuredLogDestinations"),
            s3_storage_options=params.get("S3StorageOptions"),
            tags=params.get("Tags"),
            workflow_details=params.get("WorkflowDetails"),
        )
        return json.dumps(dict(ServerId=server_id))

    def describe_server(self) -> str:
        params = json.loads(self.body)
        server = self.transfer_backend.describe_server(
            server_id=params.get("ServerId"),
        )
        return json.dumps(dict(Server=server.to_dict()))

    def delete_server(self) -> str:
        params = json.loads(self.body)
        self.transfer_backend.delete_server(
            server_id=params.get("ServerId"),
        )
        return json.dumps(dict())
