"""Handles incoming servicediscovery requests, invokes methods, returns responses."""

import json
from typing import Any, Dict, List

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

    def update_http_namespace(self) -> str:
        params = json.loads(self.body)
        _id = params.get("Id")
        updater_request_id = params.get("UpdaterRequestId")
        namespace = params.get("Namespace")
        operation_id = self.servicediscovery_backend.update_http_namespace(
            _id=_id,
            updater_request_id=updater_request_id,
            namespace_dict=namespace,
        )
        return json.dumps(dict(operationId=operation_id))

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

    def register_instance(self) -> str:
        params = json.loads(self.body)
        service_id = params.get("ServiceId")
        instance_id = params.get("InstanceId")
        creator_request_id = params.get("CreatorRequestId")
        attributes = params.get("Attributes")
        operation_id = self.servicediscovery_backend.register_instance(
            service_id=service_id,
            instance_id=instance_id,
            creator_request_id=creator_request_id,
            attributes=attributes,
        )
        return json.dumps(dict(OperationId=operation_id))

    def deregister_instance(self) -> str:
        params = json.loads(self.body)
        service_id = params.get("ServiceId")
        instance_id = params.get("InstanceId")
        operation_id = self.servicediscovery_backend.deregister_instance(
            service_id=service_id,
            instance_id=instance_id,
        )
        return json.dumps(dict(OperationId=operation_id))

    def get_instance(self) -> str:
        params = json.loads(self.body)
        service_id = params.get("ServiceId")
        instance_id = params.get("InstanceId")
        instance = self.servicediscovery_backend.get_instance(
            service_id=service_id,
            instance_id=instance_id,
        )
        return json.dumps(dict(Instance=instance.to_json()))

    def get_instances_health_status(self) -> str:
        params = json.loads(self.body)
        service_id = params.get("ServiceId")
        instances = params.get("Instances")
        max_results = params.get("MaxResults")
        next_token = params.get("NextToken")
        status_records = self.servicediscovery_backend.get_instances_health_status(
            service_id=service_id,
            instances=instances,
        )
        page, new_token = self.servicediscovery_backend.paginate(
            status_records, max_results=max_results, next_token=next_token
        )
        result: Dict[str, Any] = {"Status": {}}
        for record in page:
            result["Status"][record[0]] = record[1]
        if new_token:
            result["NextToken"] = new_token
        return json.dumps(result)

    def update_instance_custom_health_status(self) -> str:
        params = json.loads(self.body)
        service_id = params.get("ServiceId")
        instance_id = params.get("InstanceId")
        status = params.get("Status")
        self.servicediscovery_backend.update_instance_custom_health_status(
            service_id=service_id,
            instance_id=instance_id,
            status=status,
        )
        return "{}"

    def list_instances(self) -> str:
        params = json.loads(self.body)
        service_id = params.get("ServiceId")
        next_token = params.get("NextToken")
        max_results = params.get("MaxResults")
        instances = self.servicediscovery_backend.list_instances(service_id=service_id)
        page, new_token = self.servicediscovery_backend.paginate(
            instances, max_results=max_results, next_token=next_token
        )
        result: Dict[str, Any] = {"Instances": []}
        for instance in page:
            result["Instances"].append(instance.to_json())
        if new_token:
            result["NextToken"] = new_token
        return json.dumps(result)

    def discover_instances(self) -> str:
        params = json.loads(self.body)
        namespace_name = params.get("NamespaceName")
        service_name = params.get("ServiceName")
        max_results = params.get("MaxResults")
        query_parameters = params.get("QueryParameters")
        optional_parameters = params.get("OptionalParameters")
        health_status = params.get("HealthStatus")
        instances, instances_revision = (
            self.servicediscovery_backend.discover_instances(
                namespace_name=namespace_name,
                service_name=service_name,
                query_parameters=query_parameters,
                optional_parameters=optional_parameters,
                health_status=health_status,
            )
        )
        page, new_token = self.servicediscovery_backend.paginate(
            instances, max_results=max_results
        )
        result_instances: List[Dict[str, Any]] = []
        instances_revision_total = 0
        for instance in page:
            result_instances.append(
                {
                    "InstanceId": instance.instance_id,
                    "NamespaceName": namespace_name,
                    "ServiceName": service_name,
                    "Attributes": instance.attributes,
                    "HealthStatus": instance.health_status,
                }
            )
            instances_revision_total += instances_revision[instance.instance_id]
        return json.dumps(
            {
                "Instances": result_instances,
                "InstancesRevision": instances_revision_total,
            }
        )

    def discover_instances_revision(self) -> str:
        params = json.loads(self.body)
        namespace_name = params.get("NamespaceName")
        service_name = params.get("ServiceName")
        instances_revision = self.servicediscovery_backend.discover_instances_revision(
            namespace_name=namespace_name,
            service_name=service_name,
        )
        return json.dumps(dict(InstancesRevision=instances_revision))
