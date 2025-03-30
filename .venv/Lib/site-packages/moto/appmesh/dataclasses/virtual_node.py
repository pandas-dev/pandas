from dataclasses import asdict, dataclass, field
from typing import Any, Dict, List, Optional

from moto.appmesh.dataclasses.shared import (
    Duration,
    Metadata,
    MissingField,
    Status,
    Timeout,
)
from moto.appmesh.utils.common import clean_dict


@dataclass
class CertificateFile:
    certificate_chain: str

    def to_dict(self) -> Dict[str, Any]:  # type: ignore[misc]
        return {"certificateChain": self.certificate_chain}


@dataclass
class CertificateFileWithPrivateKey(CertificateFile):
    private_key: str

    def to_dict(self) -> Dict[str, Any]:  # type: ignore[misc]
        return {
            "certificateChain": self.certificate_chain,
            "privateKey": self.private_key,
        }


@dataclass
class SDS:
    secret_name: str

    def to_dict(self) -> Dict[str, Any]:  # type: ignore[misc]
        return {"secretName": self.secret_name}


@dataclass
class Certificate:
    file: Optional[CertificateFileWithPrivateKey]
    sds: Optional[SDS]

    def to_dict(self) -> Dict[str, Any]:  # type: ignore[misc]
        return clean_dict(
            {
                "file": (self.file or MissingField()).to_dict(),
                "sds": (self.sds or MissingField()).to_dict(),
            }
        )


@dataclass
class ListenerCertificateACM:
    certificate_arn: str

    def to_dict(self) -> Dict[str, Any]:  # type: ignore[misc]
        return {"certificateArn": self.certificate_arn}


@dataclass
class TLSListenerCertificate(Certificate):
    acm: Optional[ListenerCertificateACM]

    def to_dict(self) -> Dict[str, Any]:  # type: ignore[misc]
        return clean_dict(
            {
                "acm": (self.acm or MissingField()).to_dict(),
                "file": (self.file or MissingField()).to_dict(),
                "sds": (self.sds or MissingField()).to_dict(),
            }
        )


@dataclass
class Match:
    exact: List[str]

    to_dict = asdict


@dataclass
class SubjectAlternativeNames:
    match: Match

    def to_dict(self) -> Dict[str, Any]:  # type: ignore[misc]
        return {"match": self.match.to_dict()}


@dataclass
class ACM:
    certificate_authority_arns: List[str]

    def to_dict(self) -> Dict[str, Any]:  # type: ignore[misc]
        return {"certificateAuthorityArns": self.certificate_authority_arns}


@dataclass
class Trust:
    file: Optional[CertificateFile]
    sds: Optional[SDS]

    def to_dict(self) -> Dict[str, Any]:  # type: ignore[misc]
        return clean_dict(
            {
                "file": (self.file or MissingField()).to_dict(),
                "sds": (self.sds or MissingField()).to_dict(),
            }
        )


@dataclass
class BackendTrust(Trust):
    acm: Optional[ACM]

    def to_dict(self) -> Dict[str, Any]:  # type: ignore[misc]
        return clean_dict(
            {
                "acm": (self.acm or MissingField()).to_dict(),
                "file": (self.file or MissingField()).to_dict(),
                "sds": (self.sds or MissingField()).to_dict(),
            }
        )


@dataclass
class Validation:
    subject_alternative_names: Optional[SubjectAlternativeNames]

    def to_dict(self) -> Dict[str, Any]:  # type: ignore[misc]
        return clean_dict(
            {
                "subjectAlternativeNames": (
                    self.subject_alternative_names or MissingField()
                ).to_dict()
            }
        )


@dataclass
class TLSBackendValidation(Validation):
    trust: BackendTrust

    def to_dict(self) -> Dict[str, Any]:  # type: ignore[misc]
        return clean_dict(
            {
                "subjectAlternativeNames": (
                    self.subject_alternative_names or MissingField()
                ).to_dict(),
                "trust": self.trust.to_dict(),
            }
        )


@dataclass
class TLSListenerValidation(Validation):
    trust: Trust

    def to_dict(self) -> Dict[str, Any]:  # type: ignore[misc]
        return clean_dict(
            {
                "subjectAlternativeNames": (
                    self.subject_alternative_names or MissingField()
                ).to_dict(),
                "trust": self.trust.to_dict(),
            }
        )


@dataclass
class TLSClientPolicy:
    certificate: Optional[Certificate]
    enforce: Optional[bool]
    ports: Optional[List[int]]
    validation: TLSBackendValidation

    def to_dict(self) -> Dict[str, Any]:  # type: ignore[misc]
        return clean_dict(
            {
                "certificate": (self.certificate or MissingField()).to_dict(),
                "enforce": self.enforce,
                "ports": self.ports,
                "validation": self.validation.to_dict(),
            }
        )


@dataclass
class ClientPolicy:
    tls: Optional[TLSClientPolicy]

    def to_dict(self) -> Dict[str, Any]:  # type: ignore[misc]
        return clean_dict({"tls": (self.tls or MissingField()).to_dict()})


@dataclass
class BackendDefaults:
    client_policy: Optional[ClientPolicy]

    def to_dict(self) -> Dict[str, Any]:  # type: ignore[misc]
        return clean_dict(
            {"clientPolicy": (self.client_policy or MissingField()).to_dict()}
        )


@dataclass
class VirtualService:
    client_policy: Optional[ClientPolicy]
    virtual_service_name: str

    def to_dict(self) -> Dict[str, Any]:  # type: ignore[misc]
        return clean_dict(
            {
                "clientPolicy": (self.client_policy or MissingField()).to_dict(),
                "virtualServiceName": self.virtual_service_name,
            }
        )


@dataclass
class Backend:
    virtual_service: Optional[VirtualService]

    def to_dict(self) -> Dict[str, Any]:  # type: ignore[misc]
        return clean_dict(
            {"virtualService": (self.virtual_service or MissingField()).to_dict()}
        )


@dataclass
class HTTPConnection:
    max_connections: int
    max_pending_requests: Optional[int]

    def to_dict(self) -> Dict[str, Any]:  # type: ignore[misc]
        return clean_dict(
            {
                "maxConnections": self.max_connections,
                "maxPendingRequests": self.max_pending_requests,
            }
        )


@dataclass
class GRPCOrHTTP2Connection:
    max_requests: int

    def to_dict(self) -> Dict[str, Any]:  # type: ignore[misc]
        return {"maxRequests": self.max_requests}


@dataclass
class TCPConnection:
    max_connections: int

    def to_dict(self) -> Dict[str, Any]:  # type: ignore[misc]
        return {"maxConnections": self.max_connections}


@dataclass
class ConnectionPool:
    grpc: Optional[GRPCOrHTTP2Connection]
    http: Optional[HTTPConnection]
    http2: Optional[GRPCOrHTTP2Connection]
    tcp: Optional[TCPConnection]

    def to_dict(self) -> Dict[str, Any]:  # type: ignore[misc]
        return clean_dict(
            {
                "grpc": (self.grpc or MissingField()).to_dict(),
                "http": (self.http or MissingField()).to_dict(),
                "http2": (self.http2 or MissingField()).to_dict(),
                "tcp": (self.tcp or MissingField()).to_dict(),
            }
        )


@dataclass
class HealthCheck:
    healthy_threshold: int
    interval_millis: int
    path: Optional[str]
    port: Optional[int]
    protocol: str
    timeout_millis: int
    unhealthy_threshold: int

    def to_dict(self) -> Dict[str, Any]:  # type: ignore[misc]
        return clean_dict(
            {
                "healthyThreshold": self.healthy_threshold,
                "intervalMillis": self.interval_millis,
                "path": self.path,
                "port": self.port,
                "protocol": self.protocol,
                "timeoutMillis": self.timeout_millis,
                "unhealthyThreshold": self.unhealthy_threshold,
            }
        )


@dataclass
class OutlierDetection:
    base_ejection_duration: Duration
    interval: Duration
    max_ejection_percent: int
    max_server_errors: int

    def to_dict(self) -> Dict[str, Any]:  # type: ignore[misc]
        return {
            "baseEjectionDuration": self.base_ejection_duration.to_dict(),
            "interval": self.interval.to_dict(),
            "maxEjectionPercent": self.max_ejection_percent,
            "maxServerErrors": self.max_server_errors,
        }


@dataclass
class PortMapping:
    port: int
    protocol: str
    to_dict = asdict


@dataclass
class TCPTimeout:
    idle: Duration

    def to_dict(self) -> Dict[str, Any]:  # type: ignore[misc]
        return {"idle": self.idle.to_dict()}


@dataclass
class ProtocolTimeouts:
    grpc: Optional[Timeout]
    http: Optional[Timeout]
    http2: Optional[Timeout]
    tcp: Optional[TCPTimeout]

    def to_dict(self) -> Dict[str, Any]:  # type: ignore[misc]
        return clean_dict(
            {
                "grpc": (self.grpc or MissingField()).to_dict(),
                "http": (self.http or MissingField()).to_dict(),
                "http2": (self.http2 or MissingField()).to_dict(),
                "tcp": (self.tcp or MissingField()).to_dict(),
            }
        )


@dataclass
class ListenerTLS:
    certificate: TLSListenerCertificate
    mode: str
    validation: Optional[TLSListenerValidation]

    def to_dict(self) -> Dict[str, Any]:  # type: ignore[misc]
        return clean_dict(
            {
                "certificate": self.certificate.to_dict(),
                "mode": self.mode,
                "validation": (self.validation or MissingField()).to_dict(),
            }
        )


@dataclass
class Listener:
    connection_pool: Optional[ConnectionPool]
    health_check: Optional[HealthCheck]
    outlier_detection: Optional[OutlierDetection]
    port_mapping: PortMapping
    timeout: Optional[ProtocolTimeouts]
    tls: Optional[ListenerTLS]

    def to_dict(self) -> Dict[str, Any]:  # type: ignore[misc]
        return clean_dict(
            {
                "connectionPool": (self.connection_pool or MissingField()).to_dict(),
                "healthCheck": (self.health_check or MissingField()).to_dict(),
                "outlierDetection": (
                    self.outlier_detection or MissingField()
                ).to_dict(),
                "portMapping": self.port_mapping.to_dict(),
                "timeout": (self.timeout or MissingField()).to_dict(),
                "tls": (self.tls or MissingField()).to_dict(),
            }
        )


@dataclass
class KeyValue:
    key: str
    value: str
    to_dict = asdict


@dataclass
class LoggingFormat:
    json: Optional[List[KeyValue]]
    text: Optional[str]

    def to_dict(self) -> Dict[str, Any]:  # type: ignore[misc]
        return clean_dict(
            {"json": [pair.to_dict() for pair in self.json or []], "text": self.text}
        )


@dataclass
class AccessLogFile:
    format: Optional[LoggingFormat]
    path: str

    def to_dict(self) -> Dict[str, Any]:  # type: ignore[misc]
        return clean_dict(
            {"format": (self.format or MissingField()).to_dict(), "path": self.path}
        )


@dataclass
class AccessLog:
    file: Optional[AccessLogFile]

    def to_dict(self) -> Dict[str, Any]:  # type: ignore[misc]
        return clean_dict({"file": (self.file or MissingField()).to_dict()})


@dataclass
class Logging:
    access_log: Optional[AccessLog]

    def to_dict(self) -> Dict[str, Any]:  # type: ignore[misc]
        return clean_dict({"accessLog": (self.access_log or MissingField()).to_dict()})


@dataclass
class AWSCloudMap:
    attributes: Optional[List[KeyValue]]
    ip_preference: Optional[str]
    namespace_name: str
    service_name: str

    def to_dict(self) -> Dict[str, Any]:  # type: ignore[misc]
        return clean_dict(
            {
                "attributes": [
                    attribute.to_dict() for attribute in self.attributes or []
                ],
                "ipPreference": self.ip_preference,
                "namespaceName": self.namespace_name,
                "serviceName": self.service_name,
            }
        )


@dataclass
class DNS:
    hostname: str
    ip_preference: Optional[str]
    response_type: Optional[str]

    def to_dict(self) -> Dict[str, Any]:  # type: ignore[misc]
        return clean_dict(
            {
                "hostname": self.hostname,
                "ipPreference": self.ip_preference,
                "responseType": self.response_type,
            }
        )


@dataclass
class ServiceDiscovery:
    aws_cloud_map: Optional[AWSCloudMap]
    dns: Optional[DNS]

    def to_dict(self) -> Dict[str, Any]:  # type: ignore[misc]
        return clean_dict(
            {
                "awsCloudMap": (self.aws_cloud_map or MissingField()).to_dict(),
                "dns": (self.dns or MissingField()).to_dict(),
            }
        )


@dataclass
class VirtualNodeSpec:
    backend_defaults: Optional[BackendDefaults]
    backends: Optional[List[Backend]]
    listeners: Optional[List[Listener]]
    logging: Optional[Logging]
    service_discovery: Optional[ServiceDiscovery]

    def to_dict(self) -> Dict[str, Any]:  # type: ignore[misc]
        return clean_dict(
            {
                "backendDefaults": (self.backend_defaults or MissingField()).to_dict(),
                "backends": [backend.to_dict() for backend in self.backends or []],
                "listeners": [listener.to_dict() for listener in self.listeners or []],
                "logging": (self.logging or MissingField()).to_dict(),
                "serviceDiscovery": (
                    self.service_discovery or MissingField()
                ).to_dict(),
            }
        )


@dataclass
class VirtualNodeMetadata(Metadata):
    mesh_name: str = field(default="")
    virtual_node_name: str = field(default="")

    def __post_init__(self) -> None:
        if self.mesh_name == "":
            raise TypeError("__init__ missing 1 required argument: 'mesh_name'")
        if self.mesh_owner == "":
            raise TypeError("__init__ missing 1 required argument: 'route_name'")
        if self.virtual_node_name == "":
            raise TypeError("__init__ missing 1 required argument: 'virtual_node_name'")

    def formatted_for_list_api(self) -> Dict[str, Any]:  # type: ignore
        return {
            "arn": self.arn,
            "createdAt": self.created_at.strftime("%d/%m/%Y, %H:%M:%S"),
            "lastUpdatedAt": self.last_updated_at.strftime("%d/%m/%Y, %H:%M:%S"),
            "meshName": self.mesh_name,
            "meshOwner": self.mesh_owner,
            "resourceOwner": self.resource_owner,
            "version": self.version,
            "virtualNodeName": self.virtual_node_name,
        }

    def formatted_for_crud_apis(self) -> Dict[str, Any]:  # type: ignore
        return {
            "arn": self.arn,
            "createdAt": self.created_at.strftime("%d/%m/%Y, %H:%M:%S"),
            "lastUpdatedAt": self.last_updated_at.strftime("%d/%m/%Y, %H:%M:%S"),
            "meshOwner": self.mesh_owner,
            "resourceOwner": self.resource_owner,
            "uid": self.uid,
            "version": self.version,
        }


@dataclass
class VirtualNode:
    mesh_name: str
    mesh_owner: str
    metadata: VirtualNodeMetadata
    spec: VirtualNodeSpec
    virtual_node_name: str
    status: Status = field(default_factory=lambda: {"status": "ACTIVE"})
    tags: List[Dict[str, str]] = field(default_factory=list)

    def to_dict(self) -> Dict[str, Any]:  # type: ignore[misc]
        return clean_dict(
            {
                "meshName": self.mesh_name,
                "metadata": self.metadata.formatted_for_crud_apis(),
                "spec": self.spec.to_dict(),
                "status": self.status,
                "virtualNodeName": self.virtual_node_name,
            }
        )
