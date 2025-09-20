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
class RouteActionWeightedTarget:
    virtual_node: str
    weight: int
    port: Optional[int] = field(default=None)

    def to_dict(self) -> Dict[str, Any]:  # type: ignore[misc]
        return clean_dict(
            {"port": self.port, "virtualNode": self.virtual_node, "weight": self.weight}
        )


@dataclass
class RouteAction:
    weighted_targets: List[RouteActionWeightedTarget]

    def to_dict(self) -> Dict[str, Any]:  # type: ignore[misc]
        return clean_dict(
            {"weightedTargets": [target.to_dict() for target in self.weighted_targets]}
        )


@dataclass
class Range:
    start: int
    end: int
    to_dict = asdict


@dataclass
class Match:
    exact: Optional[str]
    prefix: Optional[str]
    range: Optional[Range]
    regex: Optional[str]
    suffix: Optional[str]

    def to_dict(self) -> Dict[str, Any]:  # type: ignore[misc]
        return clean_dict(
            {
                "exact": self.exact,
                "prefix": self.prefix,
                "range": (self.range or MissingField()).to_dict(),
                "regex": self.regex,
                "suffix": self.suffix,
            }
        )


@dataclass
class GrpcMetadatum:
    invert: Optional[bool]
    match: Optional[Match]
    name: str

    def to_dict(self) -> Dict[str, Any]:  # type: ignore[misc]
        return clean_dict(
            {
                "invert": self.invert,
                "match": (self.match or MissingField()).to_dict(),
                "name": self.name,
            }
        )


# same object, just different name
HttpRouteMatchHeader = GrpcMetadatum


@dataclass
class RouteMatchPath:
    exact: str
    regex: str
    to_dict = asdict


@dataclass
class QueryParameterMatch:
    exact: str
    to_dict = asdict


@dataclass
class RouteMatchQueryParameter:
    name: str
    match: Optional[QueryParameterMatch] = field(default=None)

    def to_dict(self) -> Dict[str, Any]:  # type: ignore[misc]
        return clean_dict(
            {"match": (self.match or MissingField()).to_dict(), "name": self.name}
        )


@dataclass
class HttpRouteMatch:
    headers: Optional[List[HttpRouteMatchHeader]] = field(default=None)
    method: Optional[str] = field(default=None)
    path: Optional[RouteMatchPath] = field(default=None)
    port: Optional[int] = field(default=None)
    prefix: Optional[str] = field(default=None)
    query_parameters: Optional[List[RouteMatchQueryParameter]] = field(default=None)
    scheme: Optional[str] = field(default=None)

    def to_dict(self) -> Dict[str, Any]:  # type: ignore[misc]
        return clean_dict(
            {
                "headers": [header.to_dict() for header in self.headers or []],
                "method": self.method,
                "path": (self.path or MissingField()).to_dict(),
                "port": self.port,
                "prefix": self.prefix,
                "queryParameters": [
                    param.to_dict() for param in self.query_parameters or []
                ],
                "scheme": self.scheme,
            }
        )


@dataclass
class HttpRouteRetryPolicy:
    max_retries: int
    http_retry_events: Optional[List[str]]
    per_retry_timeout: Duration
    tcp_retry_events: Optional[List[str]]

    def to_dict(self) -> Dict[str, Any]:  # type: ignore[misc]
        return clean_dict(
            {
                "httpRetryEvents": self.http_retry_events or [],
                "maxRetries": self.max_retries,
                "perRetryTimeout": self.per_retry_timeout.to_dict(),
                "tcpRetryEvents": self.tcp_retry_events or [],
            }
        )


@dataclass
class GrpcRouteMatch:
    metadata: Optional[List[GrpcMetadatum]]
    method_name: Optional[str]
    port: Optional[int]
    service_name: Optional[str]

    def to_dict(self) -> Dict[str, Any]:  # type: ignore[misc]
        return clean_dict(
            {
                "metadata": [meta.to_dict() for meta in self.metadata or []],
                "methodName": self.method_name,
                "port": self.port,
                "serviceName": self.service_name,
            }
        )


@dataclass
class GrcpRouteRetryPolicy:
    max_retries: int
    per_retry_timeout: Duration
    grpc_retry_events: Optional[List[str]]
    http_retry_events: Optional[List[str]]
    tcp_retry_events: Optional[List[str]]

    def to_dict(self) -> Dict[str, Any]:  # type: ignore[misc]
        return clean_dict(
            {
                "grpcRetryEvents": self.grpc_retry_events or [],
                "httpRetryEvents": self.http_retry_events or [],
                "maxRetries": self.max_retries,
                "perRetryTimeout": self.per_retry_timeout.to_dict(),
                "tcpRetryEvents": self.tcp_retry_events or [],
            }
        )


@dataclass
class GrpcRoute:
    action: RouteAction
    match: GrpcRouteMatch
    retry_policy: Optional[GrcpRouteRetryPolicy]
    timeout: Optional[Timeout]

    def to_dict(self) -> Dict[str, Any]:  # type: ignore[misc]
        return clean_dict(
            {
                "action": self.action.to_dict(),
                "match": self.match.to_dict(),
                "retryPolicy": (self.retry_policy or MissingField()).to_dict(),
                "timeout": (self.timeout or MissingField()).to_dict(),
            }
        )


@dataclass
class HttpRoute:
    action: RouteAction
    match: HttpRouteMatch
    retry_policy: Optional[HttpRouteRetryPolicy] = field(default=None)
    timeout: Optional[Timeout] = field(default=None)

    def to_dict(self) -> Dict[str, Any]:  # type: ignore[misc]
        return clean_dict(
            {
                "action": self.action.to_dict(),
                "match": self.match.to_dict(),
                "retryPolicy": (self.retry_policy or MissingField()).to_dict(),
                "timeout": (self.timeout or MissingField()).to_dict(),
            }
        )


@dataclass
class TCPRouteMatch:
    port: int
    to_dict = asdict


@dataclass
class TCPRoute:
    action: RouteAction
    match: Optional[TCPRouteMatch]
    timeout: Optional[Timeout]

    def to_dict(self) -> Dict[str, Any]:  # type: ignore[misc]
        return clean_dict(
            {
                "action": self.action.to_dict(),
                "match": (self.match or MissingField()).to_dict(),
                "timeout": (self.timeout or MissingField()).to_dict(),
            }
        )


@dataclass
class RouteSpec:
    priority: Optional[int]
    grpc_route: Optional[GrpcRoute] = field(default=None)
    http_route: Optional[HttpRoute] = field(default=None)
    http2_route: Optional[HttpRoute] = field(default=None)
    tcp_route: Optional[TCPRoute] = field(default=None)

    def to_dict(self) -> Dict[str, Any]:  # type: ignore[misc]
        spec = {
            "grpcRoute": (self.grpc_route or MissingField()).to_dict(),
            "httpRoute": (self.http_route or MissingField()).to_dict(),
            "http2Route": (self.http2_route or MissingField()).to_dict(),
            "priority": self.priority,
            "tcpRoute": (self.tcp_route or MissingField()).to_dict(),
        }
        return clean_dict(spec)


@dataclass
class RouteMetadata(Metadata):
    mesh_name: str = field(default="")
    route_name: str = field(default="")
    virtual_router_name: str = field(default="")

    def __post_init__(self) -> None:
        if self.mesh_name == "":
            raise TypeError("__init__ missing 1 required argument: 'mesh_name'")
        if self.mesh_owner == "":
            raise TypeError("__init__ missing 1 required argument: 'route_name'")
        if self.virtual_router_name == "":
            raise TypeError(
                "__init__ missing 1 required argument: 'virtual_router_name'"
            )

    def formatted_for_list_api(self) -> Dict[str, Any]:  # type: ignore
        return {
            "arn": self.arn,
            "createdAt": self.created_at.strftime("%d/%m/%Y, %H:%M:%S"),
            "lastUpdatedAt": self.last_updated_at.strftime("%d/%m/%Y, %H:%M:%S"),
            "meshName": self.mesh_name,
            "meshOwner": self.mesh_owner,
            "resourceOwner": self.resource_owner,
            "routeName": self.route_name,
            "version": self.version,
            "virtualRouterName": self.virtual_router_name,
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
class Route:
    mesh_name: str
    mesh_owner: str
    metadata: RouteMetadata
    route_name: str
    spec: RouteSpec
    virtual_router_name: str
    status: Status = field(default_factory=lambda: {"status": "ACTIVE"})
    tags: List[Dict[str, str]] = field(default_factory=list)

    def to_dict(self) -> Dict[str, Any]:  # type: ignore[misc]
        return clean_dict(
            {
                "meshName": self.mesh_name,
                "metadata": self.metadata.formatted_for_crud_apis(),
                "routeName": self.route_name,
                "spec": self.spec.to_dict(),
                "status": self.status,
                "tags": self.tags,
                "virtualRouterName": self.virtual_router_name,
            }
        )
