"""Handles incoming appmesh requests, invokes methods, returns responses."""

import json

from moto.appmesh.utils.spec_parsing import (
    build_route_spec,
    build_virtual_node_spec,
    port_mappings_from_router_spec,
)
from moto.core.responses import BaseResponse

from .models import AppMeshBackend, appmesh_backends


class AppMeshResponse(BaseResponse):
    """Handler for AppMesh requests and responses."""

    def __init__(self) -> None:
        super().__init__(service_name="appmesh")

    @property
    def appmesh_backend(self) -> AppMeshBackend:
        """Return backend instance specific for this region."""
        return appmesh_backends[self.current_account][self.region]

    def create_mesh(self) -> str:
        params = json.loads(self.body)
        client_token = params.get("clientToken")
        mesh_name = params.get("meshName")
        spec = params.get("spec") or {}
        egress_filter_type = (spec.get("egressFilter") or {}).get("type")
        ip_preference = (spec.get("serviceDiscovery") or {}).get("ipPreference")
        tags = params.get("tags")
        mesh = self.appmesh_backend.create_mesh(
            client_token=client_token,
            mesh_name=mesh_name,
            egress_filter_type=egress_filter_type,
            ip_preference=ip_preference,
            tags=tags,
        )
        return json.dumps(mesh.to_dict())

    def update_mesh(self) -> str:
        params = json.loads(self.body)
        client_token = params.get("clientToken")
        mesh_name = self._get_param("meshName")
        spec = params.get("spec") or {}
        egress_filter_type = (spec.get("egressFilter") or {}).get("type")
        ip_preference = (spec.get("serviceDiscovery") or {}).get("ipPreference")
        mesh = self.appmesh_backend.update_mesh(
            client_token=client_token,
            mesh_name=mesh_name,
            egress_filter_type=egress_filter_type,
            ip_preference=ip_preference,
        )
        return json.dumps(mesh.to_dict())

    def describe_mesh(self) -> str:
        mesh_name = self._get_param(param_name="meshName", if_none="")
        mesh_owner = self._get_param("meshOwner")
        mesh = self.appmesh_backend.describe_mesh(
            mesh_name=mesh_name,
            mesh_owner=mesh_owner,
        )
        return json.dumps(mesh.to_dict())

    def delete_mesh(self) -> str:
        mesh_name = self._get_param("meshName")
        mesh = self.appmesh_backend.delete_mesh(mesh_name=mesh_name)
        return json.dumps(mesh.to_dict())

    def list_meshes(self) -> str:
        params = self._get_params()
        limit = self._get_int_param("limit")
        next_token = params.get("nextToken")
        meshes, next_token = self.appmesh_backend.list_meshes(
            limit=limit,
            next_token=next_token,
        )
        return json.dumps({"meshes": meshes, "nextToken": next_token})

    def list_tags_for_resource(self) -> str:
        params = self._get_params()
        limit = self._get_int_param("limit")
        next_token = params.get("nextToken")
        resource_arn = params.get("resourceArn")
        tags, next_token = self.appmesh_backend.list_tags_for_resource(
            limit=limit,
            next_token=next_token,
            resource_arn=resource_arn,
        )
        return json.dumps({"nextToken": next_token, "tags": tags})

    def tag_resource(self) -> str:
        params = json.loads(self.body)
        resource_arn = self._get_param("resourceArn")
        tags = params.get("tags")
        self.appmesh_backend.tag_resource(
            resource_arn=resource_arn,
            tags=tags,
        )
        return json.dumps({})

    def describe_virtual_router(self) -> str:
        mesh_name = self._get_param("meshName")
        mesh_owner = self._get_param("meshOwner")
        virtual_router_name = self._get_param("virtualRouterName")
        virtual_router = self.appmesh_backend.describe_virtual_router(
            mesh_name=mesh_name,
            mesh_owner=mesh_owner,
            virtual_router_name=virtual_router_name,
        )
        return json.dumps(virtual_router.to_dict())

    def create_virtual_router(self) -> str:
        params = json.loads(self.body)
        client_token = params.get("clientToken")
        mesh_name = self._get_param("meshName")
        mesh_owner = self._get_param("meshOwner")
        port_mappings = port_mappings_from_router_spec(params.get("spec"))
        tags = params.get("tags")
        virtual_router_name = params.get("virtualRouterName")
        virtual_router = self.appmesh_backend.create_virtual_router(
            client_token=client_token,
            mesh_name=mesh_name,
            mesh_owner=mesh_owner,
            port_mappings=port_mappings,
            tags=tags,
            virtual_router_name=virtual_router_name,
        )
        return json.dumps(virtual_router.to_dict())

    def update_virtual_router(self) -> str:
        params = json.loads(self.body)
        client_token = params.get("clientToken")
        mesh_name = self._get_param("meshName")
        mesh_owner = self._get_param("meshOwner")
        port_mappings = port_mappings_from_router_spec(params.get("spec"))
        virtual_router_name = self._get_param("virtualRouterName")
        virtual_router = self.appmesh_backend.update_virtual_router(
            client_token=client_token,
            mesh_name=mesh_name,
            mesh_owner=mesh_owner,
            port_mappings=port_mappings,
            virtual_router_name=virtual_router_name,
        )
        return json.dumps(virtual_router.to_dict())

    def delete_virtual_router(self) -> str:
        mesh_name = self._get_param("meshName")
        mesh_owner = self._get_param("meshOwner")
        virtual_router_name = self._get_param("virtualRouterName")
        virtual_router = self.appmesh_backend.delete_virtual_router(
            mesh_name=mesh_name,
            mesh_owner=mesh_owner,
            virtual_router_name=virtual_router_name,
        )
        return json.dumps(virtual_router.to_dict())

    def list_virtual_routers(self) -> str:
        limit = self._get_int_param("limit")
        mesh_name = self._get_param("meshName")
        mesh_owner = self._get_param("meshOwner")
        next_token = self._get_param("nextToken")
        virtual_routers, next_token = self.appmesh_backend.list_virtual_routers(
            limit=limit,
            mesh_name=mesh_name,
            mesh_owner=mesh_owner,
            next_token=next_token,
        )
        return json.dumps({"nextToken": next_token, "virtualRouters": virtual_routers})

    def create_route(self) -> str:
        params = json.loads(self.body)
        client_token = params.get("clientToken")
        mesh_name = self._get_param("meshName")
        mesh_owner = self._get_param("meshOwner")
        route_name = self._get_param("routeName")
        tags = params.get("tags")
        virtual_router_name = self._get_param("virtualRouterName")
        spec = build_route_spec(params.get("spec") or {})
        route = self.appmesh_backend.create_route(
            client_token=client_token,
            mesh_name=mesh_name,
            mesh_owner=mesh_owner,
            route_name=route_name,
            spec=spec,
            tags=tags,
            virtual_router_name=virtual_router_name,
        )
        return json.dumps(route.to_dict())

    def describe_route(self) -> str:
        mesh_name = self._get_param("meshName")
        mesh_owner = self._get_param("meshOwner")
        route_name = self._get_param("routeName")
        virtual_router_name = self._get_param("virtualRouterName")
        route = self.appmesh_backend.describe_route(
            mesh_name=mesh_name,
            mesh_owner=mesh_owner,
            route_name=route_name,
            virtual_router_name=virtual_router_name,
        )
        return json.dumps(route.to_dict())

    def update_route(self) -> str:
        params = json.loads(self.body)
        client_token = params.get("clientToken")
        mesh_name = self._get_param("meshName")
        mesh_owner = self._get_param("meshOwner")
        route_name = self._get_param("routeName")
        virtual_router_name = self._get_param("virtualRouterName")
        spec = build_route_spec(params.get("spec") or {})
        route = self.appmesh_backend.update_route(
            client_token=client_token,
            mesh_name=mesh_name,
            mesh_owner=mesh_owner,
            route_name=route_name,
            spec=spec,
            virtual_router_name=virtual_router_name,
        )
        return json.dumps(route.to_dict())

    def delete_route(self) -> str:
        mesh_name = self._get_param("meshName")
        mesh_owner = self._get_param("meshOwner")
        route_name = self._get_param("routeName")
        virtual_router_name = self._get_param("virtualRouterName")
        route = self.appmesh_backend.delete_route(
            mesh_name=mesh_name,
            mesh_owner=mesh_owner,
            route_name=route_name,
            virtual_router_name=virtual_router_name,
        )
        return json.dumps(route.to_dict())

    def list_routes(self) -> str:
        limit = self._get_int_param("limit")
        mesh_name = self._get_param("meshName")
        mesh_owner = self._get_param("meshOwner")
        next_token = self._get_param("nextToken")
        virtual_router_name = self._get_param("virtualRouterName")
        routes, next_token = self.appmesh_backend.list_routes(
            limit=limit,
            mesh_name=mesh_name,
            mesh_owner=mesh_owner,
            next_token=next_token,
            virtual_router_name=virtual_router_name,
        )
        return json.dumps(
            {
                "nextToken": next_token,
                "routes": [r.formatted_for_list_api() for r in routes],
            }
        )

    def describe_virtual_node(self) -> str:
        mesh_name = self._get_param("meshName")
        mesh_owner = self._get_param("meshOwner")
        virtual_node_name = self._get_param("virtualNodeName")
        virtual_node = self.appmesh_backend.describe_virtual_node(
            mesh_name=mesh_name,
            mesh_owner=mesh_owner,
            virtual_node_name=virtual_node_name,
        )
        return json.dumps(virtual_node.to_dict())

    def create_virtual_node(self) -> str:
        params = json.loads(self.body)
        client_token = params.get("clientToken")
        mesh_name = self._get_param("meshName")
        mesh_owner = self._get_param("meshOwner")
        spec = build_virtual_node_spec(params.get("spec") or {})
        tags = params.get("tags")
        virtual_node_name = params.get("virtualNodeName")
        virtual_node = self.appmesh_backend.create_virtual_node(
            client_token=client_token,
            mesh_name=mesh_name,
            mesh_owner=mesh_owner,
            spec=spec,
            tags=tags,
            virtual_node_name=virtual_node_name,
        )
        return json.dumps(virtual_node.to_dict())

    def update_virtual_node(self) -> str:
        params = json.loads(self.body)
        client_token = params.get("clientToken")
        mesh_name = self._get_param("meshName")
        mesh_owner = self._get_param("meshOwner")
        spec = build_virtual_node_spec(params.get("spec") or {})
        virtual_node_name = self._get_param("virtualNodeName")
        virtual_node = self.appmesh_backend.update_virtual_node(
            client_token=client_token,
            mesh_name=mesh_name,
            mesh_owner=mesh_owner,
            spec=spec,
            virtual_node_name=virtual_node_name,
        )
        return json.dumps(virtual_node.to_dict())

    def delete_virtual_node(self) -> str:
        mesh_name = self._get_param("meshName")
        mesh_owner = self._get_param("meshOwner")
        virtual_node_name = self._get_param("virtualNodeName")
        virtual_node = self.appmesh_backend.delete_virtual_node(
            mesh_name=mesh_name,
            mesh_owner=mesh_owner,
            virtual_node_name=virtual_node_name,
        )
        return json.dumps(virtual_node.to_dict())

    def list_virtual_nodes(self) -> str:
        limit = self._get_int_param("limit")
        mesh_name = self._get_param("meshName")
        mesh_owner = self._get_param("meshOwner")
        next_token = self._get_param("nextToken")
        virtual_nodes, next_token = self.appmesh_backend.list_virtual_nodes(
            limit=limit,
            mesh_name=mesh_name,
            mesh_owner=mesh_owner,
            next_token=next_token,
        )
        return json.dumps(
            {
                "nextToken": next_token,
                "virtualNodes": [n.formatted_for_list_api() for n in virtual_nodes],
            }
        )
