from dataclasses import dataclass, field
from datetime import datetime
from enum import Enum
from typing import Any, Dict, List, Literal, Optional, Union

from moto.core.common_models import BaseModel


class UserHomeDirectoryType(str, Enum):
    PATH = "PATH"
    LOGICAL = "LOGICAL"


class UserHomeDirectoryMappingType(str, Enum):
    FILE = "FILE"
    DIRECTORY = "DIRECTORY"


@dataclass
class User(BaseModel):
    home_directory: Optional[str]
    home_directory_type: Optional[UserHomeDirectoryType]
    policy: Optional[str]
    role: str
    user_name: str
    arn: str = field(default="", init=False)
    home_directory_mappings: List[Dict[str, Optional[str]]] = field(
        default_factory=list
    )
    posix_profile: Dict[str, Optional[Union[str, List[str]]]] = field(
        default_factory=dict
    )
    ssh_public_keys: List[Dict[str, str]] = field(default_factory=list)
    tags: List[Dict[str, str]] = field(default_factory=list)

    def __post_init__(self) -> None:
        if self.arn == "":
            self.arn = f"arn:aws:transfer:{self.user_name}:{datetime.now().strftime('%Y%m%d%H%M%S')}"

    def to_dict(self) -> Dict[str, Any]:
        user = {
            "HomeDirectory": self.home_directory,
            "HomeDirectoryType": self.home_directory_type,
            "Policy": self.policy,
            "Role": self.role,
            "UserName": self.user_name,
            "Arn": self.arn,
            "HomeDirectoryMappings": [
                {
                    "Entry": mapping.get("entry"),
                    "Target": mapping.get("target"),
                    "Type": mapping.get("type"),
                }
                for mapping in self.home_directory_mappings
            ],
            "SshPublicKeys": [
                {
                    "DateImported": key.get("date_imported"),
                    "SshPublicKeyBody": key.get("ssh_public_key_body"),
                    "SshPublicKeyId": key.get("ssh_public_key_id"),
                }
                for key in self.ssh_public_keys
            ],
            "PosixProfile": {
                "Uid": self.posix_profile.get("uid"),
                "Gid": self.posix_profile.get("gid"),
                "SecondaryGids": self.posix_profile.get("secondary_gids"),
            },
            "Tags": self.tags,
        }

        return user


class ServerDomain(str, Enum):
    S3 = "S3"
    EFS = "EFS"


class ServerEndpointType(str, Enum):
    PUBLIC = "PUBLIC"
    VPC = "VPC"
    VPC_ENDPOINT = "VPC_ENDPOINT"


class ServerIdentityProviderSftpAuthenticationMethods(str, Enum):
    PASSWORD = "PASSWORD"
    PUBLIC_KEY = "PUBLIC_KEY"
    PUBLIC_KEY_OR_PASSWORD = "PUBLIC_KEY_OR_PASSWORD"
    PUBLIC_KEY_AND_PASSWORD = "PUBLIC_KEY_AND_PASSWORD"


class ServerIdentityProviderType(str, Enum):
    SERVICE_MANAGED = "SERVICE_MANAGED"
    API_GATEWAY = "API_GATEWAY"
    AWS_DIRECTORY_SERVICE = "AWS_DIRECTORY_SERVICE"
    AWS_LAMBDA = "AWS_LAMBDA"


class ServerProtocols(str, Enum):
    SFTP = "SFTP"
    FTP = "FTP"
    FTPS = "FTPS"
    AS2 = "AS2"


class ServerState(str, Enum):
    OFFLINE = "OFFLINE"
    ONLINE = "ONLINE"
    STARTING = "STARTING"
    STOPPING = "STOPPING"
    START_FAILED = "START_FAILED"
    STOP_FAILED = "STOP_FAILED"


AS2_TRANSPORTS_TYPE = List[Literal["HTTP"]]


@dataclass
class Server(BaseModel):
    certificate: Optional[str]
    domain: Optional[ServerDomain]
    endpoint_type: Optional[ServerEndpointType]
    host_key_fingerprint: Optional[str]
    identity_provider_type: Optional[ServerIdentityProviderType]
    logging_role: Optional[str]
    post_authentication_login_banner: Optional[str]
    pre_authentication_login_banner: Optional[str]
    protocols: Optional[List[ServerProtocols]]
    security_policy_name: Optional[str]
    structured_log_destinations: Optional[List[str]]
    arn: str = field(default="", init=False)
    as2_service_managed_egress_ip_addresses: List[str] = field(default_factory=list)
    endpoint_details: Dict[str, str] = field(default_factory=dict)
    identity_provider_details: Dict[str, str] = field(default_factory=dict)
    protocol_details: Dict[str, str] = field(default_factory=dict)
    s3_storage_options: Dict[str, Optional[str]] = field(default_factory=dict)
    server_id: str = field(default="", init=False)
    state: Optional[ServerState] = ServerState.ONLINE
    tags: List[Dict[str, str]] = field(default_factory=list)
    user_count: int = field(default=0)
    workflow_details: Dict[str, List[Dict[str, str]]] = field(default_factory=dict)
    _users: List[User] = field(default_factory=list, repr=False)

    def __post_init__(self) -> None:
        if self.arn == "":
            self.arn = f"arn:aws:transfer:{self.server_id}"
        if self.server_id == "":
            self.server_id = f"{self.identity_provider_type}:{self.server_id}:{datetime.now().strftime('%Y%m%d%H%M%S')}"
        if self.as2_service_managed_egress_ip_addresses == []:
            self.as2_service_managed_egress_ip_addresses.append("0.0.0.0/0")

    def to_dict(self) -> Dict[str, Any]:
        on_upload = []
        on_partial_upload = []
        if self.workflow_details is not None:
            on_upload_workflows = self.workflow_details.get("on_upload")

            if on_upload_workflows is not None:
                for workflow in on_upload_workflows:
                    workflow_id = workflow.get("workflow_id")
                    execution_role = workflow.get("execution_role")
                    if workflow_id and execution_role:
                        on_upload.append(
                            {"WorkflowId": workflow_id, "ExecutionRole": execution_role}
                        )
            on_partial_upload_workflows = self.workflow_details.get("on_partial_upload")
            if on_partial_upload_workflows is not None:
                for workflow in on_partial_upload_workflows:
                    workflow_id = workflow.get("workflow_id")
                    execution_role = workflow.get("execution_role")
                    if workflow_id and execution_role:
                        on_partial_upload.append(
                            {"WorkflowId": workflow_id, "ExecutionRole": execution_role}
                        )
        server = {
            "Certificate": self.certificate,
            "Domain": self.domain,
            "EndpointType": self.endpoint_type,
            "HostKeyFingerprint": self.host_key_fingerprint,
            "IdentityProviderType": self.identity_provider_type,
            "LoggingRole": self.logging_role,
            "PostAuthenticationLoginBanner": self.post_authentication_login_banner,
            "PreAuthenticationLoginBanner": self.pre_authentication_login_banner,
            "Protocols": self.protocols,
            "SecurityPolicyName": self.security_policy_name,
            "StructuredLogDestinations": self.structured_log_destinations,
            "Arn": self.arn,
            "As2ServiceManagedEgressIpAddresses": self.as2_service_managed_egress_ip_addresses,
            "ServerId": self.server_id,
            "State": self.state,
            "Tags": self.tags,
            "UserCount": self.user_count,
            "EndpointDetails": {
                "AddressAllocationIds": self.endpoint_details.get(
                    "address_allocation_ids"
                ),
                "SubnetIds": self.endpoint_details.get("subnet_ids"),
                "VpcEndpointId": self.endpoint_details.get("vpc_endpoint_id"),
                "VpcId": self.endpoint_details.get("vpc_id"),
                "SecurityGroupIds": self.endpoint_details.get("security_group_ids"),
            },
            "IdentityProviderDetails": {
                "Url": self.identity_provider_details.get("url"),
                "InvocationRole": self.identity_provider_details.get("invocation_role"),
                "DirectoryId": self.identity_provider_details.get("directory_id"),
                "Function": self.identity_provider_details.get("function"),
                "SftpAuthenticationMethods": self.identity_provider_details.get(
                    "sftp_authentication_methods"
                ),
            },
            "ProtocolDetails": {
                "PassiveIp": self.protocol_details.get("passive_ip"),
                "TlsSessionResumptionMode": self.protocol_details.get(
                    "tls_session_resumption_mode"
                ),
                "SetStatOption": self.protocol_details.get("set_stat_option"),
                "As2Transports": self.protocol_details.get("as2_transports"),
            },
            "S3StorageOptions": {
                "DirectoryListingOptimization": self.s3_storage_options.get(
                    "directory_listing_optimization"
                )
            },
            "WorkflowDetails": {
                "OnUpload": on_upload,
                "OnPartialUpload": on_partial_upload,
            },
        }
        return server
