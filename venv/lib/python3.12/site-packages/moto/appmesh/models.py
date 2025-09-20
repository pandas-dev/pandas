"""AppMeshBackend class with methods for supported APIs."""

from typing import Dict, List, Literal, Optional, Union

from moto.appmesh.dataclasses.mesh import (
    Mesh,
    MeshSpec,
)
from moto.appmesh.dataclasses.route import Route, RouteMetadata, RouteSpec
from moto.appmesh.dataclasses.shared import Metadata
from moto.appmesh.dataclasses.virtual_node import (
    VirtualNode,
    VirtualNodeMetadata,
    VirtualNodeSpec,
)
from moto.appmesh.dataclasses.virtual_router import (
    PortMapping,
    VirtualRouter,
    VirtualRouterSpec,
)
from moto.appmesh.exceptions import (
    MeshNotFoundError,
    MeshOwnerDoesNotMatchError,
    ResourceNotFoundError,
    RouteNameAlreadyTakenError,
    RouteNotFoundError,
    VirtualNodeNameAlreadyTakenError,
    VirtualNodeNotFoundError,
    VirtualRouterNameAlreadyTakenError,
    VirtualRouterNotFoundError,
)
from moto.core.base_backend import BackendDict, BaseBackend
from moto.utilities.paginator import paginate

PAGINATION_MODEL = {
    "list_meshes": {
        "input_token": "next_token",
        "limit_key": "limit",
        "limit_default": 100,
        "unique_attribute": "meshName",
    },
    "list_tags_for_resource": {
        "input_token": "next_token",
        "limit_key": "limit",
        "limit_default": 100,
        "unique_attribute": ["key", "value"],
    },
    "list_virtual_routers": {
        "input_token": "next_token",
        "limit_key": "limit",
        "limit_default": 100,
        "unique_attribute": "virtualRouterName",
    },
    "list_routes": {
        "input_token": "next_token",
        "limit_key": "limit",
        "limit_default": 100,
        "unique_attribute": "route_name",
    },
    "list_virtual_nodes": {
        "input_token": "next_token",
        "limit_key": "limit",
        "limit_default": 100,
        "unique_attribute": "virtual_node_name",
    },
}


class AppMeshBackend(BaseBackend):
    """Implementation of AppMesh APIs."""

    def __init__(self, region_name: str, account_id: str) -> None:
        super().__init__(region_name, account_id)
        self.meshes: Dict[str, Mesh] = dict()

    def _validate_mesh(self, mesh_name: str, mesh_owner: Optional[str]) -> None:
        if mesh_name not in self.meshes:
            raise MeshNotFoundError(mesh_name=mesh_name)
        if (
            mesh_owner is not None
            and self.meshes[mesh_name].metadata.mesh_owner != mesh_owner
        ):
            raise MeshOwnerDoesNotMatchError(mesh_name, mesh_owner)

    def _check_virtual_node_validity(
        self,
        mesh_name: str,
        mesh_owner: Optional[str],
        virtual_node_name: str,
    ) -> None:
        self._validate_mesh(mesh_name=mesh_name, mesh_owner=mesh_owner)
        if virtual_node_name not in self.meshes[mesh_name].virtual_nodes:
            raise VirtualNodeNotFoundError(mesh_name, virtual_node_name)
        return

    def _check_virtual_node_availability(
        self,
        mesh_name: str,
        mesh_owner: Optional[str],
        virtual_node_name: str,
    ) -> None:
        self._validate_mesh(mesh_name=mesh_name, mesh_owner=mesh_owner)
        if virtual_node_name in self.meshes[mesh_name].virtual_nodes:
            raise VirtualNodeNameAlreadyTakenError(
                mesh_name=mesh_name, virtual_node_name=virtual_node_name
            )
        return

    def _check_router_availability(
        self,
        mesh_name: str,
        mesh_owner: Optional[str],
        virtual_router_name: str,
    ) -> None:
        self._validate_mesh(mesh_name=mesh_name, mesh_owner=mesh_owner)
        if virtual_router_name in self.meshes[mesh_name].virtual_routers:
            raise VirtualRouterNameAlreadyTakenError(
                virtual_router_name=virtual_router_name, mesh_name=mesh_name
            )
        return

    def _check_router_validity(
        self,
        mesh_name: str,
        mesh_owner: Optional[str],
        virtual_router_name: str,
    ) -> None:
        self._validate_mesh(mesh_name=mesh_name, mesh_owner=mesh_owner)
        if virtual_router_name not in self.meshes[mesh_name].virtual_routers:
            raise VirtualRouterNotFoundError(
                virtual_router_name=virtual_router_name, mesh_name=mesh_name
            )
        return

    def _check_route_validity(
        self,
        mesh_name: str,
        mesh_owner: Optional[str],
        virtual_router_name: str,
        route_name: str,
    ) -> None:
        self._check_router_validity(
            mesh_name=mesh_name,
            mesh_owner=mesh_owner,
            virtual_router_name=virtual_router_name,
        )
        if (
            route_name
            not in self.meshes[mesh_name].virtual_routers[virtual_router_name].routes
        ):
            raise RouteNotFoundError(
                mesh_name=mesh_name,
                virtual_router_name=virtual_router_name,
                route_name=route_name,
            )
        return

    def _check_route_availability(
        self,
        mesh_name: str,
        mesh_owner: Optional[str],
        virtual_router_name: str,
        route_name: str,
    ) -> None:
        self._check_router_validity(
            mesh_name=mesh_name,
            mesh_owner=mesh_owner,
            virtual_router_name=virtual_router_name,
        )
        if (
            route_name
            in self.meshes[mesh_name].virtual_routers[virtual_router_name].routes
        ):
            raise RouteNameAlreadyTakenError(
                mesh_name=mesh_name,
                virtual_router_name=virtual_router_name,
                route_name=route_name,
            )
        return

    def create_mesh(
        self,
        client_token: Optional[str],
        mesh_name: str,
        egress_filter_type: Optional[str],
        ip_preference: Optional[str],
        tags: Optional[List[Dict[str, str]]],
    ) -> Mesh:
        from moto.sts import sts_backends

        sts_backend = sts_backends[self.account_id]["global"]
        user_id, _, _ = sts_backend.get_caller_identity(
            self.account_id, region=self.region_name
        )

        metadata = Metadata(
            arn=f"arn:aws:appmesh:{self.region_name}:{self.account_id}:{mesh_name}",
            mesh_owner=user_id,
            resource_owner=user_id,
        )
        spec = MeshSpec(
            egress_filter={"type": egress_filter_type},
            service_discovery={"ip_preference": ip_preference},
        )
        mesh = Mesh(
            mesh_name=mesh_name,
            spec=spec,
            status={"status": "ACTIVE"},
            metadata=metadata,
            tags=tags or list(),
        )
        self.meshes[mesh_name] = mesh
        return mesh

    def update_mesh(
        self,
        client_token: Optional[str],
        mesh_name: str,
        egress_filter_type: Optional[str],
        ip_preference: Optional[str],
    ) -> Mesh:
        if mesh_name not in self.meshes:
            raise MeshNotFoundError(mesh_name=mesh_name)
        updated = False
        if egress_filter_type is not None:
            self.meshes[mesh_name].spec.egress_filter["type"] = egress_filter_type
            updated = True

        if ip_preference is not None:
            self.meshes[mesh_name].spec.service_discovery["ip_preference"] = (
                ip_preference
            )
            updated = True

        if updated:
            self.meshes[mesh_name].metadata.update_timestamp()
            self.meshes[mesh_name].metadata.version += 1
        return self.meshes[mesh_name]

    def describe_mesh(self, mesh_name: str, mesh_owner: Optional[str]) -> Mesh:
        self._validate_mesh(mesh_name=mesh_name, mesh_owner=mesh_owner)
        return self.meshes[mesh_name]

    def delete_mesh(self, mesh_name: str) -> Mesh:
        if mesh_name not in self.meshes:
            raise MeshNotFoundError(mesh_name=mesh_name)
        self.meshes[mesh_name].status["status"] = "DELETED"
        mesh = self.meshes[mesh_name]
        del self.meshes[mesh_name]
        return mesh

    @paginate(pagination_model=PAGINATION_MODEL)
    def list_meshes(self) -> List[Dict[str, Union[str, int]]]:
        return [
            {
                "arn": mesh.metadata.arn,
                "createdAt": mesh.metadata.created_at.strftime("%d/%m/%Y, %H:%M:%S"),
                "lastUpdatedAt": mesh.metadata.last_updated_at.strftime(
                    "%d/%m/%Y, %H:%M:%S"
                ),
                "meshName": mesh.mesh_name,
                "meshOwner": mesh.metadata.mesh_owner,
                "resourceOwner": mesh.metadata.resource_owner,
                "version": mesh.metadata.version,
            }
            for mesh in self.meshes.values()
        ]

    def _get_resource_with_arn(
        self, resource_arn: str
    ) -> Union[Mesh, VirtualRouter, Route, VirtualNode]:
        for mesh in self.meshes.values():
            if mesh.metadata.arn == resource_arn:
                return mesh
            for virtual_router in mesh.virtual_routers.values():
                if virtual_router.metadata.arn == resource_arn:
                    return virtual_router
                for route in virtual_router.routes.values():
                    if route.metadata.arn == resource_arn:
                        return route
            for virtual_node in mesh.virtual_nodes.values():
                if virtual_node.metadata.arn == resource_arn:
                    return virtual_node
        raise ResourceNotFoundError(resource_arn)

    @paginate(pagination_model=PAGINATION_MODEL)  # type: ignore
    def list_tags_for_resource(self, resource_arn: str) -> List[Dict[str, str]]:
        return self._get_resource_with_arn(resource_arn=resource_arn).tags

    def tag_resource(self, resource_arn: str, tags: List[Dict[str, str]]) -> None:
        if len(tags) > 0:
            resource = self._get_resource_with_arn(resource_arn=resource_arn)
            resource.tags.extend(tags)
        return

    def describe_virtual_router(
        self, mesh_name: str, mesh_owner: Optional[str], virtual_router_name: str
    ) -> VirtualRouter:
        self._check_router_validity(
            mesh_name=mesh_name,
            mesh_owner=mesh_owner,
            virtual_router_name=virtual_router_name,
        )
        return self.meshes[mesh_name].virtual_routers[virtual_router_name]

    def create_virtual_router(
        self,
        client_token: str,
        mesh_name: str,
        mesh_owner: Optional[str],
        port_mappings: List[PortMapping],
        tags: Optional[List[Dict[str, str]]],
        virtual_router_name: str,
    ) -> VirtualRouter:
        self._check_router_availability(
            mesh_name=mesh_name,
            virtual_router_name=virtual_router_name,
            mesh_owner=mesh_owner,
        )
        owner = mesh_owner or self.meshes[mesh_name].metadata.mesh_owner
        metadata = Metadata(
            mesh_owner=owner,
            resource_owner=owner,
            arn=f"arn:aws:appmesh:{self.region_name}:{self.account_id}:mesh/{mesh_name}/virtualRouter/{virtual_router_name}",
        )
        listeners: List[Dict[Literal["port_mapping"], PortMapping]] = [
            {"port_mapping": port_mapping} for port_mapping in port_mappings
        ]
        spec = VirtualRouterSpec(listeners=listeners)
        virtual_router = VirtualRouter(
            virtual_router_name=virtual_router_name,
            mesh_name=mesh_name,
            metadata=metadata,
            status={"status": "ACTIVE"},
            spec=spec,
            tags=tags or list(),
        )
        self.meshes[mesh_name].virtual_routers[virtual_router_name] = virtual_router
        return virtual_router

    def update_virtual_router(
        self,
        client_token: str,
        mesh_name: str,
        mesh_owner: Optional[str],
        port_mappings: List[PortMapping],
        virtual_router_name: str,
    ) -> VirtualRouter:
        self._check_router_validity(
            mesh_name=mesh_name,
            mesh_owner=mesh_owner,
            virtual_router_name=virtual_router_name,
        )
        listeners: List[Dict[Literal["port_mapping"], PortMapping]] = [
            {"port_mapping": port_mapping} for port_mapping in port_mappings
        ]
        spec = VirtualRouterSpec(listeners=listeners)
        virtual_router = self.meshes[mesh_name].virtual_routers[virtual_router_name]
        virtual_router.spec = spec
        virtual_router.metadata.update_timestamp()
        virtual_router.metadata.version += 1
        return virtual_router

    def delete_virtual_router(
        self, mesh_name: str, mesh_owner: Optional[str], virtual_router_name: str
    ) -> VirtualRouter:
        self._check_router_validity(
            mesh_name=mesh_name,
            mesh_owner=mesh_owner,
            virtual_router_name=virtual_router_name,
        )
        mesh = self.meshes[mesh_name]
        mesh.virtual_routers[virtual_router_name].status["status"] = "DELETED"
        virtual_router = mesh.virtual_routers[virtual_router_name]
        del mesh.virtual_routers[virtual_router_name]
        return virtual_router

    @paginate(pagination_model=PAGINATION_MODEL)
    def list_virtual_routers(
        self, mesh_name: str, mesh_owner: Optional[str]
    ) -> List[Dict[str, Union[str, int]]]:
        self._validate_mesh(mesh_name=mesh_name, mesh_owner=mesh_owner)
        return [
            {
                "arn": virtual_router.metadata.arn,
                "createdAt": virtual_router.metadata.created_at.strftime(
                    "%d/%m/%Y, %H:%M:%S"
                ),
                "lastUpdatedAt": virtual_router.metadata.last_updated_at.strftime(
                    "%d/%m/%Y, %H:%M:%S"
                ),
                "meshName": virtual_router.mesh_name,
                "meshOwner": virtual_router.metadata.mesh_owner,
                "resourceOwner": virtual_router.metadata.resource_owner,
                "version": virtual_router.metadata.version,
                "virtualRouterName": virtual_router.virtual_router_name,
            }
            for virtual_router in self.meshes[mesh_name].virtual_routers.values()
        ]

    def create_route(
        self,
        client_token: Optional[str],
        mesh_name: str,
        mesh_owner: Optional[str],
        route_name: str,
        spec: RouteSpec,
        tags: Optional[List[Dict[str, str]]],
        virtual_router_name: str,
    ) -> Route:
        self._check_route_availability(
            mesh_name=mesh_name,
            mesh_owner=mesh_owner,
            route_name=route_name,
            virtual_router_name=virtual_router_name,
        )
        owner = mesh_owner or self.meshes[mesh_name].metadata.mesh_owner
        metadata = RouteMetadata(
            arn=f"arn:aws:appmesh:{self.region_name}:{self.account_id}:mesh/{mesh_name}/virtualRouter/{virtual_router_name}/route/{route_name}",
            mesh_name=mesh_name,
            mesh_owner=owner,
            resource_owner=owner,
            route_name=route_name,
            virtual_router_name=virtual_router_name,
        )
        route = Route(
            mesh_name=mesh_name,
            mesh_owner=owner,
            metadata=metadata,
            route_name=route_name,
            spec=spec,
            tags=tags or list(),
            virtual_router_name=virtual_router_name,
        )
        self.meshes[mesh_name].virtual_routers[virtual_router_name].routes[
            route_name
        ] = route
        return route

    def describe_route(
        self,
        mesh_name: str,
        mesh_owner: Optional[str],
        route_name: str,
        virtual_router_name: str,
    ) -> Route:
        self._check_route_validity(
            mesh_name=mesh_name,
            mesh_owner=mesh_owner,
            virtual_router_name=virtual_router_name,
            route_name=route_name,
        )
        return (
            self.meshes[mesh_name]
            .virtual_routers[virtual_router_name]
            .routes[route_name]
        )

    def update_route(
        self,
        client_token: Optional[str],
        mesh_name: str,
        mesh_owner: Optional[str],
        route_name: str,
        spec: RouteSpec,
        virtual_router_name: str,
    ) -> Route:
        self._check_route_validity(
            mesh_name=mesh_name,
            mesh_owner=mesh_owner,
            virtual_router_name=virtual_router_name,
            route_name=route_name,
        )
        route = (
            self.meshes[mesh_name]
            .virtual_routers[virtual_router_name]
            .routes[route_name]
        )
        route.spec = spec
        route.metadata.version += 1
        route.metadata.update_timestamp()
        return route

    def delete_route(
        self,
        mesh_name: str,
        mesh_owner: Optional[str],
        route_name: str,
        virtual_router_name: str,
    ) -> Route:
        self._check_route_validity(
            mesh_name=mesh_name,
            mesh_owner=mesh_owner,
            virtual_router_name=virtual_router_name,
            route_name=route_name,
        )
        route = (
            self.meshes[mesh_name]
            .virtual_routers[virtual_router_name]
            .routes[route_name]
        )
        route.status["status"] = "DELETED"
        del (
            self.meshes[mesh_name]
            .virtual_routers[virtual_router_name]
            .routes[route_name]
        )
        return route

    @paginate(pagination_model=PAGINATION_MODEL)
    def list_routes(
        self,
        mesh_name: str,
        mesh_owner: Optional[str],
        virtual_router_name: str,
    ) -> List[RouteMetadata]:
        self._check_router_validity(
            mesh_name=mesh_name,
            mesh_owner=mesh_owner,
            virtual_router_name=virtual_router_name,
        )
        virtual_router = self.meshes[mesh_name].virtual_routers[virtual_router_name]
        return [route.metadata for route in virtual_router.routes.values()]

    def describe_virtual_node(
        self, mesh_name: str, mesh_owner: Optional[str], virtual_node_name: str
    ) -> VirtualNode:
        self._check_virtual_node_validity(
            mesh_name=mesh_name,
            mesh_owner=mesh_owner,
            virtual_node_name=virtual_node_name,
        )
        return self.meshes[mesh_name].virtual_nodes[virtual_node_name]

    def create_virtual_node(
        self,
        client_token: Optional[str],
        mesh_name: str,
        mesh_owner: Optional[str],
        spec: VirtualNodeSpec,
        tags: Optional[List[Dict[str, str]]],
        virtual_node_name: str,
    ) -> VirtualNode:
        self._check_virtual_node_availability(
            mesh_name=mesh_name,
            mesh_owner=mesh_owner,
            virtual_node_name=virtual_node_name,
        )
        owner = mesh_owner or self.meshes[mesh_name].metadata.mesh_owner
        metadata = VirtualNodeMetadata(
            arn=f"arn:aws:appmesh:{self.region_name}:{self.account_id}:mesh/{mesh_name}/virtualNode/{virtual_node_name}",
            mesh_name=mesh_name,
            mesh_owner=owner,
            resource_owner=owner,
            virtual_node_name=virtual_node_name,
        )
        virtual_node = VirtualNode(
            mesh_name=mesh_name,
            mesh_owner=owner,
            metadata=metadata,
            spec=spec,
            tags=tags or list(),
            virtual_node_name=virtual_node_name,
        )
        self.meshes[mesh_name].virtual_nodes[virtual_node_name] = virtual_node
        return virtual_node

    def update_virtual_node(
        self,
        client_token: Optional[str],
        mesh_name: str,
        mesh_owner: Optional[str],
        spec: VirtualNodeSpec,
        virtual_node_name: str,
    ) -> VirtualNode:
        self._check_virtual_node_validity(
            mesh_name=mesh_name,
            mesh_owner=mesh_owner,
            virtual_node_name=virtual_node_name,
        )
        virtual_node = self.meshes[mesh_name].virtual_nodes[virtual_node_name]
        virtual_node.spec = spec
        virtual_node.metadata.version += 1
        virtual_node.metadata.update_timestamp()
        return virtual_node

    def delete_virtual_node(
        self, mesh_name: str, mesh_owner: Optional[str], virtual_node_name: str
    ) -> VirtualNode:
        self._check_virtual_node_validity(
            mesh_name=mesh_name,
            mesh_owner=mesh_owner,
            virtual_node_name=virtual_node_name,
        )
        virtual_node = self.meshes[mesh_name].virtual_nodes[virtual_node_name]
        virtual_node.status["status"] = "DELETED"
        del self.meshes[mesh_name].virtual_nodes[virtual_node_name]
        return virtual_node

    @paginate(pagination_model=PAGINATION_MODEL)
    def list_virtual_nodes(
        self,
        mesh_name: str,
        mesh_owner: Optional[str],
    ) -> List[VirtualNodeMetadata]:
        self._validate_mesh(mesh_name=mesh_name, mesh_owner=mesh_owner)
        virtual_nodes = self.meshes[mesh_name].virtual_nodes
        return [virtual_node.metadata for virtual_node in virtual_nodes.values()]


appmesh_backends = BackendDict(AppMeshBackend, "appmesh")
