"""Handles incoming servicediscovery requests, invokes methods, returns responses."""

import json

from moto.core.common_types import TYPE_RESPONSE
from moto.core.responses import BaseResponse

from .models import ServiceDiscoveryBackend, servicediscovery_backends


class ServiceDiscoveryResponse(BaseResponse):
    def __init__(self) -> None:
        super().__init__(service_name="servicediscovery")

    @property
    def servicediscovery_backend(self) -> ServiceDiscoveryBackend:
        """Return backend instance specific for this region."""
        return servicediscovery_backends[self.current_account][self.region]

    def list_namespaces(self) -> TYPE_RESPONSE:
        namespaces = self.servicediscovery_backend.list_namespaces()
        return 200, {}, json.dumps({"Namespaces": [ns.to_json() for ns in namespaces]})

    def create_http_namespace(self) -> str:
        params = json.loads(self.body)
        name = params.get("Name")
        creator_request_id = params.get("CreatorRequestId")
        description = params.get("Description")
        tags = params.get("Tags")
        operation_id = self.servicediscovery_backend.create_http_namespace(
            name=name,
            creator_request_id=creator_request_id,
            description=description,
            tags=tags,
        )
        return json.dumps(dict(OperationId=operation_id))

    def delete_namespace(self) -> str:
        params = json.loads(self.body)
        namespace_id = params.get("Id")
        operation_id = self.servicediscovery_backend.delete_namespace(
            namespace_id=namespace_id
        )
        return json.dumps(dict(OperationId=operation_id))

    def list_operations(self) -> TYPE_RESPONSE:
        operations = self.servicediscovery_backend.list_operations()
        return (
            200,
            {},
            json.dumps({"Operations": [o.to_json(short=True) for o in operations]}),
        )

    def get_operation(self) -> str:
        params = json.loads(self.body)
        operation_id = params.get("OperationId")
        operation = self.servicediscovery_backend.get_operation(
            operation_id=operation_id
        )
        return json.dumps(dict(Operation=operation.to_json()))

    def get_namespace(self) -> str:
        params = json.loads(self.body)
        namespace_id = params.get("Id")
        namespace = self.servicediscovery_backend.get_namespace(
            namespace_id=namespace_id
        )
        return json.dumps(dict(Namespace=namespace.to_json()))

    def tag_resource(self) -> str:
        params = json.loads(self.body)
        resource_arn = params.get("ResourceARN")
        tags = params.get("Tags")
        self.servicediscovery_backend.tag_resource(resource_arn=resource_arn, tags=tags)
        return "{}"

    def untag_resource(self) -> str:
        params = json.loads(self.body)
        resource_arn = params.get("ResourceARN")
        tag_keys = params.get("TagKeys")
        self.servicediscovery_backend.untag_resource(
            resource_arn=resource_arn, tag_keys=tag_keys
        )
        return "{}"

    def list_tags_for_resource(self) -> TYPE_RESPONSE:
        params = json.loads(self.body)
        resource_arn = params.get("ResourceARN")
        tags = self.servicediscovery_backend.list_tags_for_resource(
            resource_arn=resource_arn
        )
        return 200, {}, json.dumps(tags)

    def create_private_dns_namespace(self) -> str:
        params = json.loads(self.body)
        name = params.get("Name")
        creator_request_id = params.get("CreatorRequestId")
        description = params.get("Description")
        vpc = params.get("Vpc")
        tags = params.get("Tags")
        properties = params.get("Properties")
        operation_id = self.servicediscovery_backend.create_private_dns_namespace(
            name=name,
            creator_request_id=creator_request_id,
            description=description,
            vpc=vpc,
            tags=tags,
            properties=properties,
        )
        return json.dumps(dict(OperationId=operation_id))

    def create_public_dns_namespace(self) -> str:
        params = json.loads(self.body)
        name = params.get("Name")
        creator_request_id = params.get("CreatorRequestId")
        description = params.get("Description")
        tags = params.get("Tags")
        properties = params.get("Properties")
        operation_id = self.servicediscovery_backend.create_public_dns_namespace(
            name=name,
            creator_request_id=creator_request_id,
            description=description,
            tags=tags,
            properties=properties,
        )
        return json.dumps(dict(OperationId=operation_id))

    def create_service(self) -> str:
        params = json.loads(self.body)
        name = params.get("Name")
        namespace_id = params.get("NamespaceId")
        creator_request_id = params.get("CreatorRequestId")
        description = params.get("Description")
        dns_config = params.get("DnsConfig")
        health_check_config = params.get("HealthCheckConfig")
        health_check_custom_config = params.get("HealthCheckCustomConfig")
        tags = params.get("Tags")
        service_type = params.get("Type")
        service = self.servicediscovery_backend.create_service(
            name=name,
            namespace_id=namespace_id,
            creator_request_id=creator_request_id,
            description=description,
            dns_config=dns_config,
            health_check_config=health_check_config,
            health_check_custom_config=health_check_custom_config,
            tags=tags,
            service_type=service_type,
        )
        return json.dumps(dict(Service=service.to_json()))

    def get_service(self) -> str:
        params = json.loads(self.body)
        service_id = params.get("Id")
        service = self.servicediscovery_backend.get_service(service_id=service_id)
        return json.dumps(dict(Service=service.to_json()))

    def delete_service(self) -> str:
        params = json.loads(self.body)
        service_id = params.get("Id")
        self.servicediscovery_backend.delete_service(service_id=service_id)
        return "{}"

    def list_services(self) -> str:
        services = self.servicediscovery_backend.list_services()
        return json.dumps(dict(Services=[s.to_json() for s in services]))

    def update_service(self) -> str:
        params = json.loads(self.body)
        service_id = params.get("Id")
        details = params.get("Service")
        operation_id = self.servicediscovery_backend.update_service(
            service_id=service_id, details=details
        )
        return json.dumps(dict(OperationId=operation_id))

    def update_private_dns_namespace(self) -> str:
        params = json.loads(self.body)
        _id = params.get("Id")
        description = params["Namespace"].get("Description")
        properties = params["Namespace"].get("Properties", {}).get("DnsProperties")
        operation_id = self.servicediscovery_backend.update_private_dns_namespace(
            _id=_id,
            description=description,
            properties=properties,
        )
        return json.dumps(dict(OperationId=operation_id))

    def update_public_dns_namespace(self) -> str:
        params = json.loads(self.body)
        _id = params.get("Id")
        description = params["Namespace"].get("Description")
        properties = params["Namespace"].get("Properties", {}).get("DnsProperties")
        operation_id = self.servicediscovery_backend.update_public_dns_namespace(
            _id=_id,
            description=description,
            properties=properties,
        )
        return json.dumps(dict(OperationId=operation_id))
