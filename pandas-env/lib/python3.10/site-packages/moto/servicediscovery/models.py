import string
from typing import Any, Dict, Iterable, List, Optional, Tuple

from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import BaseModel
from moto.core.utils import unix_time
from moto.moto_api._internal import mock_random as random
from moto.utilities.tagging_service import TaggingService
from moto.utilities.utils import get_partition

from .exceptions import (
    ConflictingDomainExists,
    CustomHealthNotFound,
    InstanceNotFound,
    InvalidInput,
    NamespaceNotFound,
    OperationNotFound,
    ServiceNotFound,
)


def random_id(size: int) -> str:
    return "".join(
        [random.choice(string.ascii_lowercase + string.digits) for _ in range(size)]
    )


class Namespace(BaseModel):
    def __init__(
        self,
        account_id: str,
        region: str,
        name: str,
        ns_type: str,
        creator_request_id: str,
        description: str,
        dns_properties: Dict[str, Any],
        http_properties: Dict[str, Any],
        vpc: Optional[str] = None,
    ):
        self.id = f"ns-{random_id(20)}"
        self.arn = f"arn:{get_partition(region)}:servicediscovery:{region}:{account_id}:namespace/{self.id}"
        self.name = name
        self.type = ns_type
        self.creator_request_id = creator_request_id
        self.description = description
        self.dns_properties = dns_properties
        self.http_properties = http_properties
        self.vpc = vpc
        self.created = unix_time()
        self.updated = unix_time()

    def to_json(self) -> Dict[str, Any]:
        return {
            "Arn": self.arn,
            "Id": self.id,
            "Name": self.name,
            "Description": self.description,
            "Type": self.type,
            "Properties": {
                "DnsProperties": self.dns_properties,
                "HttpProperties": self.http_properties,
            },
            "CreateDate": self.created,
            "UpdateDate": self.updated,
            "CreatorRequestId": self.creator_request_id,
        }


class Service(BaseModel):
    def __init__(
        self,
        account_id: str,
        region: str,
        name: str,
        namespace_id: str,
        description: str,
        creator_request_id: str,
        dns_config: Dict[str, Any],
        health_check_config: Dict[str, Any],
        health_check_custom_config: Dict[str, int],
        service_type: str,
    ):
        self.id = f"srv-{random_id(8)}"
        self.arn = f"arn:{get_partition(region)}:servicediscovery:{region}:{account_id}:service/{self.id}"
        self.name = name
        self.namespace_id = namespace_id
        self.description = description
        self.creator_request_id = creator_request_id
        self.dns_config: Optional[Dict[str, Any]] = dns_config
        self.health_check_config = health_check_config
        self.health_check_custom_config = health_check_custom_config
        self.service_type = service_type
        self.created = unix_time()
        self.instances: List[ServiceInstance] = []
        self.instances_revision: Dict[str, int] = {}

    def update(self, details: Dict[str, Any]) -> None:
        if "Description" in details:
            self.description = details["Description"]
        if "DnsConfig" in details:
            if self.dns_config is None:
                self.dns_config = {}
            self.dns_config["DnsRecords"] = details["DnsConfig"]["DnsRecords"]
        else:
            # From the docs:
            #    If you omit any existing DnsRecords or HealthCheckConfig configurations from an UpdateService request,
            #    the configurations are deleted from the service.
            self.dns_config = None
        if "HealthCheckConfig" in details:
            self.health_check_config = details["HealthCheckConfig"]

    def to_json(self) -> Dict[str, Any]:
        return {
            "Arn": self.arn,
            "Id": self.id,
            "Name": self.name,
            "NamespaceId": self.namespace_id,
            "CreateDate": self.created,
            "Description": self.description,
            "CreatorRequestId": self.creator_request_id,
            "DnsConfig": self.dns_config,
            "HealthCheckConfig": self.health_check_config,
            "HealthCheckCustomConfig": self.health_check_custom_config,
            "Type": self.service_type,
        }


class ServiceInstance(BaseModel):
    def __init__(
        self,
        service_id: str,
        instance_id: str,
        creator_request_id: Optional[str] = None,
        attributes: Optional[Dict[str, str]] = None,
    ):
        self.service_id = service_id
        self.instance_id = instance_id
        self.attributes = attributes if attributes else {}
        self.creator_request_id = (
            creator_request_id if creator_request_id else random_id(32)
        )
        self.health_status = "HEALTHY"

    def to_json(self) -> Dict[str, Any]:
        return {
            "Id": self.instance_id,
            "CreatorRequestId": self.creator_request_id,
            "Attributes": self.attributes,
        }


class Operation(BaseModel):
    def __init__(self, operation_type: str, targets: Dict[str, str]):
        super().__init__()
        self.id = f"{random_id(32)}-{random_id(8)}"
        self.status = "SUCCESS"
        self.operation_type = operation_type
        self.created = unix_time()
        self.updated = unix_time()
        self.targets = targets

    def to_json(self, short: bool = False) -> Dict[str, Any]:
        if short:
            return {"Id": self.id, "Status": self.status}
        else:
            return {
                "Id": self.id,
                "Status": self.status,
                "Type": self.operation_type,
                "CreateDate": self.created,
                "UpdateDate": self.updated,
                "Targets": self.targets,
            }


class ServiceDiscoveryBackend(BaseBackend):
    """Implementation of ServiceDiscovery APIs."""

    def __init__(self, region_name: str, account_id: str):
        super().__init__(region_name, account_id)
        self.operations: Dict[str, Operation] = dict()
        self.namespaces: Dict[str, Namespace] = dict()
        self.services: Dict[str, Service] = dict()
        self.tagger = TaggingService()

    def list_namespaces(self) -> Iterable[Namespace]:
        """
        Pagination or the Filters-parameter is not yet implemented
        """
        return list(self.namespaces.values())

    def create_http_namespace(
        self,
        name: str,
        creator_request_id: str,
        description: str,
        tags: List[Dict[str, str]],
    ) -> str:
        namespace = Namespace(
            account_id=self.account_id,
            region=self.region_name,
            name=name,
            ns_type="HTTP",
            creator_request_id=creator_request_id,
            description=description,
            dns_properties={"SOA": {}},
            http_properties={"HttpName": name},
        )
        self.namespaces[namespace.id] = namespace
        if tags:
            self.tagger.tag_resource(namespace.arn, tags)
        operation_id = self._create_operation(
            "CREATE_NAMESPACE", targets={"NAMESPACE": namespace.id}
        )
        return operation_id

    def _create_operation(self, op_type: str, targets: Dict[str, str]) -> str:
        operation = Operation(operation_type=op_type, targets=targets)
        self.operations[operation.id] = operation
        return operation.id

    def delete_namespace(self, namespace_id: str) -> str:
        if namespace_id not in self.namespaces:
            raise NamespaceNotFound(namespace_id)
        del self.namespaces[namespace_id]
        operation_id = self._create_operation(
            op_type="DELETE_NAMESPACE", targets={"NAMESPACE": namespace_id}
        )
        return operation_id

    def get_namespace(self, namespace_id: str) -> Namespace:
        if namespace_id not in self.namespaces:
            raise NamespaceNotFound(namespace_id)
        return self.namespaces[namespace_id]

    def list_operations(self) -> Iterable[Operation]:
        """
        Pagination or the Filters-argument is not yet implemented
        """
        # Operations for namespaces will only be listed as long as namespaces exist
        self.operations = {
            op_id: op
            for op_id, op in self.operations.items()
            if op.targets.get("NAMESPACE") in self.namespaces
        }
        return self.operations.values()

    def get_operation(self, operation_id: str) -> Operation:
        if operation_id not in self.operations:
            raise OperationNotFound()
        return self.operations[operation_id]

    def tag_resource(self, resource_arn: str, tags: List[Dict[str, str]]) -> None:
        self.tagger.tag_resource(resource_arn, tags)

    def untag_resource(self, resource_arn: str, tag_keys: List[str]) -> None:
        self.tagger.untag_resource_using_names(resource_arn, tag_keys)

    def list_tags_for_resource(
        self, resource_arn: str
    ) -> Dict[str, List[Dict[str, str]]]:
        return self.tagger.list_tags_for_resource(resource_arn)

    def create_private_dns_namespace(
        self,
        name: str,
        creator_request_id: str,
        description: str,
        vpc: str,
        tags: List[Dict[str, str]],
        properties: Dict[str, Any],
    ) -> str:
        for namespace in self.namespaces.values():
            if namespace.vpc == vpc:
                raise ConflictingDomainExists(vpc)
        dns_properties = (properties or {}).get("DnsProperties", {})
        dns_properties["HostedZoneId"] = "hzi"
        namespace = Namespace(
            account_id=self.account_id,
            region=self.region_name,
            name=name,
            ns_type="DNS_PRIVATE",
            creator_request_id=creator_request_id,
            description=description,
            dns_properties=dns_properties,
            http_properties={},
            vpc=vpc,
        )
        self.namespaces[namespace.id] = namespace
        if tags:
            self.tagger.tag_resource(namespace.arn, tags)
        operation_id = self._create_operation(
            "CREATE_NAMESPACE", targets={"NAMESPACE": namespace.id}
        )
        return operation_id

    def create_public_dns_namespace(
        self,
        name: str,
        creator_request_id: str,
        description: str,
        tags: List[Dict[str, str]],
        properties: Dict[str, Any],
    ) -> str:
        dns_properties = (properties or {}).get("DnsProperties", {})
        dns_properties["HostedZoneId"] = "hzi"
        namespace = Namespace(
            account_id=self.account_id,
            region=self.region_name,
            name=name,
            ns_type="DNS_PUBLIC",
            creator_request_id=creator_request_id,
            description=description,
            dns_properties=dns_properties,
            http_properties={},
        )
        self.namespaces[namespace.id] = namespace
        if tags:
            self.tagger.tag_resource(namespace.arn, tags)
        operation_id = self._create_operation(
            "CREATE_NAMESPACE", targets={"NAMESPACE": namespace.id}
        )
        return operation_id

    def create_service(
        self,
        name: str,
        namespace_id: str,
        creator_request_id: str,
        description: str,
        dns_config: Dict[str, Any],
        health_check_config: Dict[str, Any],
        health_check_custom_config: Dict[str, Any],
        tags: List[Dict[str, str]],
        service_type: str,
    ) -> Service:
        service = Service(
            account_id=self.account_id,
            region=self.region_name,
            name=name,
            namespace_id=namespace_id,
            description=description,
            creator_request_id=creator_request_id,
            dns_config=dns_config,
            health_check_config=health_check_config,
            health_check_custom_config=health_check_custom_config,
            service_type=service_type,
        )
        self.services[service.id] = service
        if tags:
            self.tagger.tag_resource(service.arn, tags)
        return service

    def get_service(self, service_id: str) -> Service:
        if service_id not in self.services:
            raise ServiceNotFound(service_id)
        return self.services[service_id]

    def delete_service(self, service_id: str) -> None:
        self.services.pop(service_id, None)

    def list_services(self) -> Iterable[Service]:
        """
        Pagination or the Filters-argument is not yet implemented
        """
        return self.services.values()

    def update_service(self, service_id: str, details: Dict[str, Any]) -> str:
        service = self.get_service(service_id)
        service.update(details=details)
        operation_id = self._create_operation(
            "UPDATE_SERVICE", targets={"SERVICE": service.id}
        )
        return operation_id

    def update_http_namespace(
        self,
        _id: str,
        namespace_dict: Dict[str, Any],
        updater_request_id: Optional[str] = None,
    ) -> str:
        if "Description" not in namespace_dict:
            raise InvalidInput("Description is required")

        namespace = self.get_namespace(namespace_id=_id)
        if updater_request_id is None:
            # Unused as the operation cannot fail
            updater_request_id = random_id(32)
        namespace.description = namespace_dict["Description"]
        if "Properties" in namespace_dict:
            if "HttpProperties" in namespace_dict["Properties"]:
                namespace.http_properties = namespace_dict["Properties"][
                    "HttpProperties"
                ]
        operation_id = self._create_operation(
            "UPDATE_NAMESPACE", targets={"NAMESPACE": namespace.id}
        )
        return operation_id

    def update_private_dns_namespace(
        self, _id: str, description: str, properties: Dict[str, Any]
    ) -> str:
        namespace = self.get_namespace(namespace_id=_id)
        if description is not None:
            namespace.description = description
        if properties is not None:
            namespace.dns_properties = properties
        operation_id = self._create_operation(
            "UPDATE_NAMESPACE", targets={"NAMESPACE": namespace.id}
        )
        return operation_id

    def update_public_dns_namespace(
        self, _id: str, description: str, properties: Dict[str, Any]
    ) -> str:
        namespace = self.get_namespace(namespace_id=_id)
        if description is not None:
            namespace.description = description
        if properties is not None:
            namespace.dns_properties = properties
        operation_id = self._create_operation(
            "UPDATE_NAMESPACE", targets={"NAMESPACE": namespace.id}
        )
        return operation_id

    def register_instance(
        self,
        service_id: str,
        instance_id: str,
        creator_request_id: str,
        attributes: Dict[str, str],
    ) -> str:
        service = self.get_service(service_id)
        instance = ServiceInstance(
            service_id=service_id,
            instance_id=instance_id,
            creator_request_id=creator_request_id,
            attributes=attributes,
        )
        service.instances.append(instance)
        service.instances_revision[instance_id] = (
            service.instances_revision.get(instance_id, 0) + 1
        )
        operation_id = self._create_operation(
            "REGISTER_INSTANCE", targets={"INSTANCE": instance_id}
        )
        return operation_id

    def deregister_instance(self, service_id: str, instance_id: str) -> str:
        service = self.get_service(service_id)
        i = 0
        while i < len(service.instances):
            instance = service.instances[i]
            if instance.instance_id == instance_id:
                service.instances.remove(instance)
                service.instances_revision[instance_id] = (
                    service.instances_revision.get(instance_id, 0) + 1
                )
                operation_id = self._create_operation(
                    "DEREGISTER_INSTANCE", targets={"INSTANCE": instance_id}
                )
                return operation_id
            i += 1
        raise InstanceNotFound(instance_id)

    def list_instances(self, service_id: str) -> List[ServiceInstance]:
        service = self.get_service(service_id)
        return service.instances

    def get_instance(self, service_id: str, instance_id: str) -> ServiceInstance:
        for instance in self.list_instances(service_id):
            if instance.instance_id == instance_id:
                return instance
        raise InstanceNotFound(instance_id)

    def get_instances_health_status(
        self,
        service_id: str,
        instances: Optional[List[str]] = None,
    ) -> List[Tuple[str, str]]:
        service = self.get_service(service_id)
        status = []
        if instances is None:
            instances = [instance.instance_id for instance in service.instances]
        if not isinstance(instances, list):
            raise InvalidInput("Instances must be a list")
        filtered_instances = [
            instance
            for instance in service.instances
            if instance.instance_id in instances
        ]
        for instance in filtered_instances:
            status.append((instance.instance_id, instance.health_status))
        return status

    def update_instance_custom_health_status(
        self, service_id: str, instance_id: str, status: str
    ) -> None:
        if status not in ["HEALTHY", "UNHEALTHY"]:
            raise CustomHealthNotFound(service_id)
        instance = self.get_instance(service_id, instance_id)
        instance.health_status = status

    def _filter_instances(
        self,
        instances: List[ServiceInstance],
        query_parameters: Optional[Dict[str, str]] = None,
        optional_parameters: Optional[Dict[str, str]] = None,
        health_status: Optional[str] = None,
    ) -> List[ServiceInstance]:
        if query_parameters is None:
            query_parameters = {}
        if optional_parameters is None:
            optional_parameters = {}
        if health_status is None:
            health_status = "ALL"

        filtered_instances = []
        has_healthy = False
        for instance in instances:
            # Filter out instances with mismatching health status
            if (
                health_status not in ["ALL", "HEALTHY_OR_ELSE_ALL"]
                and instance.health_status != health_status
            ):
                continue
            # Record if there is at least one healthy instance for HEALTHY_OR_ELSE_ALL
            if instance.health_status == "HEALTHY":
                has_healthy = True
            # Filter out instances with mismatching query parameters
            matches_query = True
            for param in query_parameters:
                if instance.attributes.get(param) != query_parameters[param]:
                    matches_query = False
                    break
            if not matches_query:
                continue
            # Add instance to the list if it passed all filters
            filtered_instances.append(instance)
        # Handle HEALTHY_OR_ELSE_ALL
        if has_healthy and health_status == "HEALTHY_OR_ELSE_ALL":
            filtered_instances = [
                instance
                for instance in filtered_instances
                if instance.health_status == "HEALTHY"
            ]
        # Filter out instances with mismatching optional parameters
        opt_filtered_instances = []
        for instance in filtered_instances:
            matches_optional = True
            for param in optional_parameters:
                if instance.attributes.get(param) != optional_parameters[param]:
                    matches_optional = False
                    break
            if matches_optional:
                opt_filtered_instances.append(instance)
        # If no instances passed the optional parameters, return the original filtered list
        return opt_filtered_instances if opt_filtered_instances else filtered_instances

    def discover_instances(
        self,
        namespace_name: str,
        service_name: str,
        query_parameters: Optional[Dict[str, str]] = None,
        optional_parameters: Optional[Dict[str, str]] = None,
        health_status: Optional[str] = None,
    ) -> Tuple[List[ServiceInstance], Dict[str, int]]:
        if query_parameters is None:
            query_parameters = {}
        if optional_parameters is None:
            optional_parameters = {}
        if health_status is None:
            health_status = "ALL"
        if health_status not in ["HEALTHY", "UNHEALTHY", "ALL", "HEALTHY_OR_ELSE_ALL"]:
            raise InvalidInput("Invalid health status")
        try:
            namespace = [
                ns for ns in self.list_namespaces() if ns.name == namespace_name
            ][0]
        except IndexError:
            raise NamespaceNotFound(namespace_name)
        try:
            service = [
                srv
                for srv in self.list_services()
                if srv.name == service_name and srv.namespace_id == namespace.id
            ][0]
        except IndexError:
            raise ServiceNotFound(service_name)
        instances = self.list_instances(service.id)
        # Filter instances based on query parameters, optional parameters, and health status
        final_instances = self._filter_instances(
            instances, query_parameters, optional_parameters, health_status
        )
        # Get the revision number for each instance that passed the filters
        instance_revisions = {
            instance.instance_id: service.instances_revision.get(
                instance.instance_id, 0
            )
            for instance in final_instances
        }
        return final_instances, instance_revisions

    def discover_instances_revision(
        self, namespace_name: str, service_name: str
    ) -> int:
        return sum(self.discover_instances(namespace_name, service_name)[1].values())

    def paginate(
        self,
        items: List[Any],
        max_results: Optional[int] = None,
        next_token: Optional[str] = None,
    ) -> Tuple[List[Any], Optional[str]]:
        """
        Paginates a list of items. If called without optional parameters, the entire list is returned as-is.
        """
        # Default to beginning of list
        if next_token is None:
            next_token = "0"
        # Return empty list if next_token is invalid
        if not next_token.isdigit():
            return [], None
        # Default to the entire list
        if max_results is None:
            max_results = len(items)
        new_token = int(next_token) + max_results
        # If the new token overflows the list, return the rest of the list
        if new_token >= len(items):
            return items[int(next_token) :], None
        return items[int(next_token) : new_token], str(new_token)


servicediscovery_backends = BackendDict(ServiceDiscoveryBackend, "servicediscovery")
