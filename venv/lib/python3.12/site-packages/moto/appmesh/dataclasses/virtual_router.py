from dataclasses import dataclass, field
from typing import Any, Dict, List, Literal, Optional

from moto.appmesh.dataclasses.route import Route
from moto.appmesh.dataclasses.shared import Metadata, Status


@dataclass
class PortMapping:
    port: Optional[int]
    protocol: Optional[str]


@dataclass
class VirtualRouterSpec:
    listeners: List[Dict[Literal["port_mapping"], PortMapping]]


@dataclass
class VirtualRouter:
    mesh_name: str
    metadata: Metadata
    spec: VirtualRouterSpec
    status: Status
    virtual_router_name: str
    routes: Dict[str, Route] = field(default_factory=dict)
    tags: List[Dict[str, str]] = field(default_factory=list)

    def to_dict(self) -> Dict[str, Any]:  # type ignore[misc]
        return {
            "meshName": self.mesh_name,
            "metadata": {
                "arn": self.metadata.arn,
                "createdAt": self.metadata.created_at.strftime("%d/%m/%Y, %H:%M:%S"),
                "lastUpdatedAt": self.metadata.last_updated_at.strftime(
                    "%d/%m/%Y, %H:%M:%S"
                ),
                "meshOwner": self.metadata.mesh_owner,
                "resourceOwner": self.metadata.resource_owner,
                "uid": self.metadata.uid,
                "version": self.metadata.version,
            },
            "spec": {
                "listeners": [
                    {
                        "portMapping": {
                            "port": listener["port_mapping"].port,
                            "protocol": listener["port_mapping"].protocol,
                        }
                    }
                    for listener in self.spec.listeners
                ]
            },
            "status": self.status,
            "virtualRouterName": self.virtual_router_name,
        }
