import json

from moto.core.common_types import TYPE_RESPONSE
from moto.core.responses import BaseResponse
from moto.core.utils import iso_8601_datetime_with_milliseconds

from .models import GreengrassBackend, greengrass_backends


class GreengrassResponse(BaseResponse):
    def __init__(self) -> None:
        super().__init__(service_name="greengrass")

    @property
    def greengrass_backend(self) -> GreengrassBackend:
        return greengrass_backends[self.current_account][self.region]

    def list_core_definitions(self) -> str:
        res = self.greengrass_backend.list_core_definitions()
        return json.dumps(
            {"Definitions": [core_definition.to_dict() for core_definition in res]}
        )

    def create_core_definition(self) -> TYPE_RESPONSE:
        name = self._get_param("Name")
        initial_version = self._get_param("InitialVersion")
        res = self.greengrass_backend.create_core_definition(
            name=name, initial_version=initial_version
        )
        return 201, {"status": 201}, json.dumps(res.to_dict())

    def get_core_definition(self) -> str:
        core_definition_id = self.path.split("/")[-1]
        res = self.greengrass_backend.get_core_definition(
            core_definition_id=core_definition_id
        )
        return json.dumps(res.to_dict())

    def delete_core_definition(self) -> str:
        core_definition_id = self.path.split("/")[-1]
        self.greengrass_backend.delete_core_definition(
            core_definition_id=core_definition_id
        )
        return json.dumps({})

    def update_core_definition(self) -> str:
        core_definition_id = self.path.split("/")[-1]
        name = self._get_param("Name")
        self.greengrass_backend.update_core_definition(
            core_definition_id=core_definition_id, name=name
        )
        return "{}"

    def create_core_definition_version(self) -> TYPE_RESPONSE:
        core_definition_id = self.path.split("/")[-2]
        cores = self._get_param("Cores")

        res = self.greengrass_backend.create_core_definition_version(
            core_definition_id=core_definition_id, cores=cores
        )
        return 201, {"status": 201}, json.dumps(res.to_dict())

    def list_core_definition_versions(self) -> str:
        core_definition_id = self.path.split("/")[-2]
        res = self.greengrass_backend.list_core_definition_versions(core_definition_id)
        return json.dumps(
            {"Versions": [core_def_ver.to_dict() for core_def_ver in res]}
        )

    def get_core_definition_version(self) -> str:
        core_definition_id = self.path.split("/")[-3]
        core_definition_version_id = self.path.split("/")[-1]
        res = self.greengrass_backend.get_core_definition_version(
            core_definition_id=core_definition_id,
            core_definition_version_id=core_definition_version_id,
        )
        return json.dumps(res.to_dict(include_detail=True))

    def create_device_definition(self) -> TYPE_RESPONSE:
        name = self._get_param("Name")
        initial_version = self._get_param("InitialVersion")
        res = self.greengrass_backend.create_device_definition(
            name=name, initial_version=initial_version
        )
        return 201, {"status": 201}, json.dumps(res.to_dict())

    def list_device_definitions(self) -> str:
        res = self.greengrass_backend.list_device_definitions()
        return json.dumps(
            {"Definitions": [device_definition.to_dict() for device_definition in res]}
        )

    def create_device_definition_version(self) -> TYPE_RESPONSE:
        device_definition_id = self.path.split("/")[-2]
        devices = self._get_param("Devices")

        res = self.greengrass_backend.create_device_definition_version(
            device_definition_id=device_definition_id, devices=devices
        )
        return 201, {"status": 201}, json.dumps(res.to_dict())

    def list_device_definition_versions(self) -> str:
        device_definition_id = self.path.split("/")[-2]
        res = self.greengrass_backend.list_device_definition_versions(
            device_definition_id
        )
        return json.dumps(
            {"Versions": [device_def_ver.to_dict() for device_def_ver in res]}
        )

    def get_device_definition(self) -> str:
        device_definition_id = self.path.split("/")[-1]
        res = self.greengrass_backend.get_device_definition(
            device_definition_id=device_definition_id
        )
        return json.dumps(res.to_dict())

    def delete_device_definition(self) -> str:
        device_definition_id = self.path.split("/")[-1]
        self.greengrass_backend.delete_device_definition(
            device_definition_id=device_definition_id
        )
        return json.dumps({})

    def update_device_definition(self) -> str:
        device_definition_id = self.path.split("/")[-1]
        name = self._get_param("Name")
        self.greengrass_backend.update_device_definition(
            device_definition_id=device_definition_id, name=name
        )
        return "{}"

    def get_device_definition_version(self) -> str:
        device_definition_id = self.path.split("/")[-3]
        device_definition_version_id = self.path.split("/")[-1]
        res = self.greengrass_backend.get_device_definition_version(
            device_definition_id=device_definition_id,
            device_definition_version_id=device_definition_version_id,
        )
        return json.dumps(res.to_dict(include_detail=True))

    def create_resource_definition(self) -> TYPE_RESPONSE:
        initial_version = self._get_param("InitialVersion")
        name = self._get_param("Name")
        res = self.greengrass_backend.create_resource_definition(
            name=name, initial_version=initial_version
        )
        return 201, {"status": 201}, json.dumps(res.to_dict())

    def list_resource_definitions(self) -> str:
        res = self.greengrass_backend.list_resource_definitions()
        return json.dumps({"Definitions": [i.to_dict() for i in res]})

    def get_resource_definition(self) -> str:
        resource_definition_id = self.path.split("/")[-1]
        res = self.greengrass_backend.get_resource_definition(
            resource_definition_id=resource_definition_id
        )
        return json.dumps(res.to_dict())

    def delete_resource_definition(self) -> str:
        resource_definition_id = self.path.split("/")[-1]
        self.greengrass_backend.delete_resource_definition(
            resource_definition_id=resource_definition_id
        )
        return "{}"

    def update_resource_definition(self) -> str:
        resource_definition_id = self.path.split("/")[-1]
        name = self._get_param("Name")
        self.greengrass_backend.update_resource_definition(
            resource_definition_id=resource_definition_id, name=name
        )
        return "{}"

    def create_resource_definition_version(self) -> TYPE_RESPONSE:
        resource_definition_id = self.path.split("/")[-2]
        resources = self._get_param("Resources")

        res = self.greengrass_backend.create_resource_definition_version(
            resource_definition_id=resource_definition_id, resources=resources
        )
        return 201, {"status": 201}, json.dumps(res.to_dict())

    def list_resource_definition_versions(self) -> str:
        resource_device_definition_id = self.path.split("/")[-2]
        res = self.greengrass_backend.list_resource_definition_versions(
            resource_device_definition_id
        )

        return json.dumps(
            {"Versions": [resource_def_ver.to_dict() for resource_def_ver in res]}
        )

    def get_resource_definition_version(self) -> str:
        resource_definition_id = self.path.split("/")[-3]
        resource_definition_version_id = self.path.split("/")[-1]
        res = self.greengrass_backend.get_resource_definition_version(
            resource_definition_id=resource_definition_id,
            resource_definition_version_id=resource_definition_version_id,
        )
        return json.dumps(res.to_dict())

    def create_function_definition(self) -> TYPE_RESPONSE:
        initial_version = self._get_param("InitialVersion")
        name = self._get_param("Name")
        res = self.greengrass_backend.create_function_definition(
            name=name, initial_version=initial_version
        )
        return 201, {"status": 201}, json.dumps(res.to_dict())

    def list_function_definitions(self) -> str:
        res = self.greengrass_backend.list_function_definitions()
        return json.dumps(
            {"Definitions": [func_definition.to_dict() for func_definition in res]}
        )

    def get_function_definition(self) -> str:
        function_definition_id = self.path.split("/")[-1]
        res = self.greengrass_backend.get_function_definition(
            function_definition_id=function_definition_id,
        )
        return json.dumps(res.to_dict())

    def delete_function_definition(self) -> str:
        function_definition_id = self.path.split("/")[-1]
        self.greengrass_backend.delete_function_definition(
            function_definition_id=function_definition_id,
        )
        return "{}"

    def update_function_definition(self) -> str:
        function_definition_id = self.path.split("/")[-1]
        name = self._get_param("Name")
        self.greengrass_backend.update_function_definition(
            function_definition_id=function_definition_id, name=name
        )
        return "{}"

    def create_function_definition_version(self) -> TYPE_RESPONSE:
        default_config = self._get_param("DefaultConfig")
        function_definition_id = self.path.split("/")[-2]
        functions = self._get_param("Functions")

        res = self.greengrass_backend.create_function_definition_version(
            default_config=default_config,
            function_definition_id=function_definition_id,
            functions=functions,
        )
        return 201, {"status": 201}, json.dumps(res.to_dict())

    def list_function_definition_versions(self) -> str:
        function_definition_id = self.path.split("/")[-2]
        res = self.greengrass_backend.list_function_definition_versions(
            function_definition_id=function_definition_id
        )
        versions = [i.to_dict() for i in res.values()]
        return json.dumps({"Versions": versions})

    def get_function_definition_version(self) -> str:
        function_definition_id = self.path.split("/")[-3]
        function_definition_version_id = self.path.split("/")[-1]
        res = self.greengrass_backend.get_function_definition_version(
            function_definition_id=function_definition_id,
            function_definition_version_id=function_definition_version_id,
        )
        return json.dumps(res.to_dict())

    def create_subscription_definition(self) -> TYPE_RESPONSE:
        initial_version = self._get_param("InitialVersion")
        name = self._get_param("Name")
        res = self.greengrass_backend.create_subscription_definition(
            name=name, initial_version=initial_version
        )
        return 201, {"status": 201}, json.dumps(res.to_dict())

    def list_subscription_definitions(self) -> str:
        res = self.greengrass_backend.list_subscription_definitions()
        return json.dumps(
            {
                "Definitions": [
                    subscription_definition.to_dict() for subscription_definition in res
                ]
            }
        )

    def get_subscription_definition(self) -> str:
        subscription_definition_id = self.path.split("/")[-1]
        res = self.greengrass_backend.get_subscription_definition(
            subscription_definition_id=subscription_definition_id
        )
        return json.dumps(res.to_dict())

    def delete_subscription_definition(self) -> str:
        subscription_definition_id = self.path.split("/")[-1]
        self.greengrass_backend.delete_subscription_definition(
            subscription_definition_id=subscription_definition_id
        )
        return "{}"

    def update_subscription_definition(self) -> str:
        subscription_definition_id = self.path.split("/")[-1]
        name = self._get_param("Name")
        self.greengrass_backend.update_subscription_definition(
            subscription_definition_id=subscription_definition_id, name=name
        )
        return "{}"

    def create_subscription_definition_version(self) -> TYPE_RESPONSE:
        subscription_definition_id = self.path.split("/")[-2]
        subscriptions = self._get_param("Subscriptions")
        res = self.greengrass_backend.create_subscription_definition_version(
            subscription_definition_id=subscription_definition_id,
            subscriptions=subscriptions,
        )
        return 201, {"status": 201}, json.dumps(res.to_dict())

    def list_subscription_definition_versions(self) -> str:
        subscription_definition_id = self.path.split("/")[-2]
        res = self.greengrass_backend.list_subscription_definition_versions(
            subscription_definition_id=subscription_definition_id
        )
        versions = [i.to_dict() for i in res.values()]
        return json.dumps({"Versions": versions})

    def get_subscription_definition_version(self) -> str:
        subscription_definition_id = self.path.split("/")[-3]
        subscription_definition_version_id = self.path.split("/")[-1]
        res = self.greengrass_backend.get_subscription_definition_version(
            subscription_definition_id=subscription_definition_id,
            subscription_definition_version_id=subscription_definition_version_id,
        )
        return json.dumps(res.to_dict())

    def create_group(self) -> TYPE_RESPONSE:
        initial_version = self._get_param("InitialVersion")
        name = self._get_param("Name")
        res = self.greengrass_backend.create_group(
            name=name, initial_version=initial_version
        )
        return 201, {"status": 201}, json.dumps(res.to_dict())

    def list_groups(self) -> str:
        res = self.greengrass_backend.list_groups()
        return json.dumps({"Groups": [group.to_dict() for group in res]})

    def get_group(self) -> str:
        group_id = self.path.split("/")[-1]
        res = self.greengrass_backend.get_group(group_id=group_id)
        return json.dumps(res.to_dict())  # type: ignore

    def delete_group(self) -> str:
        group_id = self.path.split("/")[-1]
        self.greengrass_backend.delete_group(group_id=group_id)
        return "{}"

    def update_group(self) -> str:
        group_id = self.path.split("/")[-1]
        name = self._get_param("Name")
        self.greengrass_backend.update_group(group_id=group_id, name=name)
        return "{}"

    def create_group_version(self) -> TYPE_RESPONSE:
        group_id = self.path.split("/")[-2]

        core_definition_version_arn = self._get_param("CoreDefinitionVersionArn")
        device_definition_version_arn = self._get_param("DeviceDefinitionVersionArn")
        function_definition_version_arn = self._get_param(
            "FunctionDefinitionVersionArn"
        )
        resource_definition_version_arn = self._get_param(
            "ResourceDefinitionVersionArn"
        )
        subscription_definition_version_arn = self._get_param(
            "SubscriptionDefinitionVersionArn"
        )

        res = self.greengrass_backend.create_group_version(
            group_id=group_id,
            core_definition_version_arn=core_definition_version_arn,
            device_definition_version_arn=device_definition_version_arn,
            function_definition_version_arn=function_definition_version_arn,
            resource_definition_version_arn=resource_definition_version_arn,
            subscription_definition_version_arn=subscription_definition_version_arn,
        )
        return 201, {"status": 201}, json.dumps(res.to_dict())

    def list_group_versions(self) -> str:
        group_id = self.path.split("/")[-2]
        res = self.greengrass_backend.list_group_versions(group_id=group_id)
        return json.dumps({"Versions": [group_ver.to_dict() for group_ver in res]})

    def get_group_version(self) -> str:
        group_id = self.path.split("/")[-3]
        group_version_id = self.path.split("/")[-1]
        res = self.greengrass_backend.get_group_version(
            group_id=group_id,
            group_version_id=group_version_id,
        )
        return json.dumps(res.to_dict(include_detail=True))

    def create_deployment(self) -> str:
        group_id = self.path.split("/")[-2]
        group_version_id = self._get_param("GroupVersionId")
        deployment_type = self._get_param("DeploymentType")
        deployment_id = self._get_param("DeploymentId")

        res = self.greengrass_backend.create_deployment(
            group_id=group_id,
            group_version_id=group_version_id,
            deployment_type=deployment_type,
            deployment_id=deployment_id,
        )
        return json.dumps(res.to_dict())

    def list_deployments(self) -> str:
        group_id = self.path.split("/")[-2]
        res = self.greengrass_backend.list_deployments(group_id=group_id)

        deployments = [deployment.to_dict(include_detail=True) for deployment in res]

        return json.dumps({"Deployments": deployments})

    def get_deployment_status(self) -> str:
        group_id = self.path.split("/")[-4]
        deployment_id = self.path.split("/")[-2]

        res = self.greengrass_backend.get_deployment_status(
            group_id=group_id,
            deployment_id=deployment_id,
        )
        return json.dumps(res.to_dict())

    def reset_deployments(self) -> str:
        group_id = self.path.split("/")[-3]

        res = self.greengrass_backend.reset_deployments(group_id=group_id)
        return json.dumps(res.to_dict())

    def associate_role_to_group(self) -> str:
        group_id = self.path.split("/")[-2]
        role_arn = self._get_param("RoleArn")
        res = self.greengrass_backend.associate_role_to_group(
            group_id=group_id,
            role_arn=role_arn,
        )
        return json.dumps(res.to_dict())

    def get_associated_role(self) -> str:
        group_id = self.path.split("/")[-2]
        res = self.greengrass_backend.get_associated_role(group_id=group_id)
        return json.dumps(res.to_dict(include_detail=True))

    def disassociate_role_from_group(self) -> str:
        group_id = self.path.split("/")[-2]
        self.greengrass_backend.disassociate_role_from_group(group_id=group_id)
        return json.dumps({"DisassociatedAt": iso_8601_datetime_with_milliseconds()})
