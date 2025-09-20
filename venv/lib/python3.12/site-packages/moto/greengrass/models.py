import json
import re
from collections import OrderedDict
from datetime import datetime
from typing import Any, Dict, Iterable, List, Optional

from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import BaseModel
from moto.core.utils import iso_8601_datetime_with_milliseconds, utcnow
from moto.moto_api._internal import mock_random
from moto.utilities.utils import get_partition

from .exceptions import (
    GreengrassClientError,
    IdNotFoundException,
    InvalidContainerDefinitionException,
    InvalidInputException,
    MissingCoreException,
    ResourceNotFoundException,
    VersionNotFoundException,
)


class FakeCoreDefinition(BaseModel):
    def __init__(self, account_id: str, region_name: str, name: str):
        self.region_name = region_name
        self.name = name
        self.id = str(mock_random.uuid4())
        self.arn = f"arn:{get_partition(region_name)}:greengrass:{region_name}:{account_id}:greengrass/definition/cores/{self.id}"
        self.created_at_datetime = utcnow()
        self.latest_version = ""
        self.latest_version_arn = ""

    def to_dict(self) -> Dict[str, Any]:
        return {
            "Arn": self.arn,
            "CreationTimestamp": iso_8601_datetime_with_milliseconds(
                self.created_at_datetime
            ),
            "Id": self.id,
            "LastUpdatedTimestamp": iso_8601_datetime_with_milliseconds(
                self.created_at_datetime
            ),
            "LatestVersion": self.latest_version,
            "LatestVersionArn": self.latest_version_arn,
            "Name": self.name,
        }


class FakeCoreDefinitionVersion(BaseModel):
    def __init__(
        self,
        account_id: str,
        region_name: str,
        core_definition_id: str,
        definition: Dict[str, Any],
    ):
        self.region_name = region_name
        self.core_definition_id = core_definition_id
        self.definition = definition
        self.version = str(mock_random.uuid4())
        self.arn = f"arn:{get_partition(region_name)}:greengrass:{region_name}:{account_id}:greengrass/definition/cores/{self.core_definition_id}/versions/{self.version}"
        self.created_at_datetime = utcnow()

    def to_dict(self, include_detail: bool = False) -> Dict[str, Any]:
        obj: Dict[str, Any] = {
            "Arn": self.arn,
            "CreationTimestamp": iso_8601_datetime_with_milliseconds(
                self.created_at_datetime
            ),
            "Id": self.core_definition_id,
            "Version": self.version,
        }

        if include_detail:
            obj["Definition"] = self.definition

        return obj


class FakeDeviceDefinition(BaseModel):
    def __init__(
        self,
        account_id: str,
        region_name: str,
        name: str,
        initial_version: Dict[str, Any],
    ):
        self.region_name = region_name
        self.id = str(mock_random.uuid4())
        self.arn = f"arn:{get_partition(region_name)}:greengrass:{region_name}:{account_id}:greengrass/definition/devices/{self.id}"
        self.created_at_datetime = utcnow()
        self.update_at_datetime = utcnow()
        self.latest_version = ""
        self.latest_version_arn = ""
        self.name = name
        self.initial_version = initial_version

    def to_dict(self) -> Dict[str, Any]:
        res = {
            "Arn": self.arn,
            "CreationTimestamp": iso_8601_datetime_with_milliseconds(
                self.created_at_datetime
            ),
            "Id": self.id,
            "LastUpdatedTimestamp": iso_8601_datetime_with_milliseconds(
                self.update_at_datetime
            ),
            "LatestVersion": self.latest_version,
            "LatestVersionArn": self.latest_version_arn,
        }
        if self.name is not None:
            res["Name"] = self.name
        return res


class FakeDeviceDefinitionVersion(BaseModel):
    def __init__(
        self,
        account_id: str,
        region_name: str,
        device_definition_id: str,
        devices: List[Dict[str, Any]],
    ):
        self.region_name = region_name
        self.device_definition_id = device_definition_id
        self.devices = devices
        self.version = str(mock_random.uuid4())
        self.arn = f"arn:{get_partition(region_name)}:greengrass:{region_name}:{account_id}:greengrass/definition/devices/{self.device_definition_id}/versions/{self.version}"
        self.created_at_datetime = utcnow()

    def to_dict(self, include_detail: bool = False) -> Dict[str, Any]:
        obj: Dict[str, Any] = {
            "Arn": self.arn,
            "CreationTimestamp": iso_8601_datetime_with_milliseconds(
                self.created_at_datetime
            ),
            "Id": self.device_definition_id,
            "Version": self.version,
        }

        if include_detail:
            obj["Definition"] = {"Devices": self.devices}

        return obj


class FakeResourceDefinition(BaseModel):
    def __init__(
        self,
        account_id: str,
        region_name: str,
        name: str,
        initial_version: Dict[str, Any],
    ):
        self.region_name = region_name
        self.id = str(mock_random.uuid4())
        self.arn = f"arn:{get_partition(region_name)}:greengrass:{region_name}:{account_id}:greengrass/definition/resources/{self.id}"
        self.created_at_datetime = utcnow()
        self.update_at_datetime = utcnow()
        self.latest_version = ""
        self.latest_version_arn = ""
        self.name = name
        self.initial_version = initial_version

    def to_dict(self) -> Dict[str, Any]:
        return {
            "Arn": self.arn,
            "CreationTimestamp": iso_8601_datetime_with_milliseconds(
                self.created_at_datetime
            ),
            "Id": self.id,
            "LastUpdatedTimestamp": iso_8601_datetime_with_milliseconds(
                self.update_at_datetime
            ),
            "LatestVersion": self.latest_version,
            "LatestVersionArn": self.latest_version_arn,
            "Name": self.name,
        }


class FakeResourceDefinitionVersion(BaseModel):
    def __init__(
        self,
        account_id: str,
        region_name: str,
        resource_definition_id: str,
        resources: List[Dict[str, Any]],
    ):
        self.region_name = region_name
        self.resource_definition_id = resource_definition_id
        self.resources = resources
        self.version = str(mock_random.uuid4())
        self.arn = f"arn:{get_partition(region_name)}:greengrass:{region_name}:{account_id}:greengrass/definition/resources/{self.resource_definition_id}/versions/{self.version}"
        self.created_at_datetime = utcnow()

    def to_dict(self) -> Dict[str, Any]:
        return {
            "Arn": self.arn,
            "CreationTimestamp": iso_8601_datetime_with_milliseconds(
                self.created_at_datetime
            ),
            "Definition": {"Resources": self.resources},
            "Id": self.resource_definition_id,
            "Version": self.version,
        }


class FakeFunctionDefinition(BaseModel):
    def __init__(
        self,
        account_id: str,
        region_name: str,
        name: str,
        initial_version: Dict[str, Any],
    ):
        self.region_name = region_name
        self.id = str(mock_random.uuid4())
        self.arn = f"arn:{get_partition(self.region_name)}:greengrass:{self.region_name}:{account_id}:greengrass/definition/functions/{self.id}"
        self.created_at_datetime = utcnow()
        self.update_at_datetime = utcnow()
        self.latest_version = ""
        self.latest_version_arn = ""
        self.name = name
        self.initial_version = initial_version

    def to_dict(self) -> Dict[str, Any]:
        res = {
            "Arn": self.arn,
            "CreationTimestamp": iso_8601_datetime_with_milliseconds(
                self.created_at_datetime
            ),
            "Id": self.id,
            "LastUpdatedTimestamp": iso_8601_datetime_with_milliseconds(
                self.update_at_datetime
            ),
            "LatestVersion": self.latest_version,
            "LatestVersionArn": self.latest_version_arn,
        }
        if self.name is not None:
            res["Name"] = self.name
        return res


class FakeFunctionDefinitionVersion(BaseModel):
    def __init__(
        self,
        account_id: str,
        region_name: str,
        function_definition_id: str,
        functions: List[Dict[str, Any]],
        default_config: Dict[str, Any],
    ):
        self.region_name = region_name
        self.function_definition_id = function_definition_id
        self.functions = functions
        self.default_config = default_config
        self.version = str(mock_random.uuid4())
        self.arn = f"arn:{get_partition(self.region_name)}:greengrass:{self.region_name}:{account_id}:greengrass/definition/functions/{self.function_definition_id}/versions/{self.version}"
        self.created_at_datetime = utcnow()

    def to_dict(self) -> Dict[str, Any]:
        return {
            "Arn": self.arn,
            "CreationTimestamp": iso_8601_datetime_with_milliseconds(
                self.created_at_datetime
            ),
            "Definition": {"Functions": self.functions},
            "Id": self.function_definition_id,
            "Version": self.version,
        }


class FakeSubscriptionDefinition(BaseModel):
    def __init__(
        self,
        account_id: str,
        region_name: str,
        name: str,
        initial_version: Dict[str, Any],
    ):
        self.region_name = region_name
        self.id = str(mock_random.uuid4())
        self.arn = f"arn:{get_partition(self.region_name)}:greengrass:{self.region_name}:{account_id}:greengrass/definition/subscriptions/{self.id}"
        self.created_at_datetime = utcnow()
        self.update_at_datetime = utcnow()
        self.latest_version = ""
        self.latest_version_arn = ""
        self.name = name
        self.initial_version = initial_version

    def to_dict(self) -> Dict[str, Any]:
        return {
            "Arn": self.arn,
            "CreationTimestamp": iso_8601_datetime_with_milliseconds(
                self.created_at_datetime
            ),
            "Id": self.id,
            "LastUpdatedTimestamp": iso_8601_datetime_with_milliseconds(
                self.update_at_datetime
            ),
            "LatestVersion": self.latest_version,
            "LatestVersionArn": self.latest_version_arn,
            "Name": self.name,
        }


class FakeSubscriptionDefinitionVersion(BaseModel):
    def __init__(
        self,
        account_id: str,
        region_name: str,
        subscription_definition_id: str,
        subscriptions: List[Dict[str, Any]],
    ):
        self.region_name = region_name
        self.subscription_definition_id = subscription_definition_id
        self.subscriptions = subscriptions
        self.version = str(mock_random.uuid4())
        self.arn = f"arn:{get_partition(self.region_name)}:greengrass:{self.region_name}:{account_id}:greengrass/definition/subscriptions/{self.subscription_definition_id}/versions/{self.version}"
        self.created_at_datetime = utcnow()

    def to_dict(self) -> Dict[str, Any]:
        return {
            "Arn": self.arn,
            "CreationTimestamp": iso_8601_datetime_with_milliseconds(
                self.created_at_datetime
            ),
            "Definition": {"Subscriptions": self.subscriptions},
            "Id": self.subscription_definition_id,
            "Version": self.version,
        }


class FakeGroup(BaseModel):
    def __init__(self, account_id: str, region_name: str, name: str):
        self.region_name = region_name
        self.group_id = str(mock_random.uuid4())
        self.name = name
        self.arn = f"arn:{get_partition(self.region_name)}:greengrass:{self.region_name}:{account_id}:greengrass/groups/{self.group_id}"
        self.created_at_datetime = utcnow()
        self.last_updated_datetime = utcnow()
        self.latest_version = ""
        self.latest_version_arn = ""

    def to_dict(self) -> Dict[str, Any]:
        obj = {
            "Arn": self.arn,
            "CreationTimestamp": iso_8601_datetime_with_milliseconds(
                self.created_at_datetime
            ),
            "Id": self.group_id,
            "LastUpdatedTimestamp": iso_8601_datetime_with_milliseconds(
                self.last_updated_datetime
            ),
            "LatestVersion": self.latest_version,
            "LatestVersionArn": self.latest_version_arn,
            "Name": self.name,
        }
        return obj


class FakeGroupVersion(BaseModel):
    def __init__(
        self,
        account_id: str,
        region_name: str,
        group_id: str,
        core_definition_version_arn: Optional[str],
        device_definition_version_arn: Optional[str],
        function_definition_version_arn: Optional[str],
        resource_definition_version_arn: Optional[str],
        subscription_definition_version_arn: Optional[str],
    ):
        self.region_name = region_name
        self.group_id = group_id
        self.version = str(mock_random.uuid4())
        self.arn = f"arn:{get_partition(self.region_name)}:greengrass:{self.region_name}:{account_id}:greengrass/groups/{self.group_id}/versions/{self.version}"
        self.created_at_datetime = utcnow()
        self.core_definition_version_arn = core_definition_version_arn
        self.device_definition_version_arn = device_definition_version_arn
        self.function_definition_version_arn = function_definition_version_arn
        self.resource_definition_version_arn = resource_definition_version_arn
        self.subscription_definition_version_arn = subscription_definition_version_arn

    def to_dict(self, include_detail: bool = False) -> Dict[str, Any]:
        definition = {}
        if self.core_definition_version_arn:
            definition["CoreDefinitionVersionArn"] = self.core_definition_version_arn

        if self.device_definition_version_arn:
            definition["DeviceDefinitionVersionArn"] = (
                self.device_definition_version_arn
            )

        if self.function_definition_version_arn:
            definition["FunctionDefinitionVersionArn"] = (
                self.function_definition_version_arn
            )

        if self.resource_definition_version_arn:
            definition["ResourceDefinitionVersionArn"] = (
                self.resource_definition_version_arn
            )

        if self.subscription_definition_version_arn:
            definition["SubscriptionDefinitionVersionArn"] = (
                self.subscription_definition_version_arn
            )

        obj: Dict[str, Any] = {
            "Arn": self.arn,
            "CreationTimestamp": iso_8601_datetime_with_milliseconds(
                self.created_at_datetime
            ),
            "Id": self.group_id,
            "Version": self.version,
        }

        if include_detail:
            obj["Definition"] = definition

        return obj


class FakeDeployment(BaseModel):
    def __init__(
        self,
        account_id: str,
        region_name: str,
        group_id: str,
        group_arn: str,
        deployment_type: str,
    ):
        self.region_name = region_name
        self.id = str(mock_random.uuid4())
        self.group_id = group_id
        self.group_arn = group_arn
        self.created_at_datetime = utcnow()
        self.update_at_datetime = utcnow()
        self.deployment_status = "InProgress"
        self.deployment_type = deployment_type
        self.arn = f"arn:{get_partition(self.region_name)}:greengrass:{self.region_name}:{account_id}:/greengrass/groups/{self.group_id}/deployments/{self.id}"

    def to_dict(self, include_detail: bool = False) -> Dict[str, Any]:
        obj = {"DeploymentId": self.id, "DeploymentArn": self.arn}

        if include_detail:
            obj["CreatedAt"] = iso_8601_datetime_with_milliseconds(
                self.created_at_datetime
            )
            obj["DeploymentType"] = self.deployment_type
            obj["GroupArn"] = self.group_arn

        return obj


class FakeAssociatedRole(BaseModel):
    def __init__(self, role_arn: str):
        self.role_arn = role_arn
        self.associated_at = utcnow()

    def to_dict(self, include_detail: bool = False) -> Dict[str, Any]:
        obj = {"AssociatedAt": iso_8601_datetime_with_milliseconds(self.associated_at)}
        if include_detail:
            obj["RoleArn"] = self.role_arn

        return obj


class FakeDeploymentStatus(BaseModel):
    def __init__(
        self,
        deployment_type: str,
        updated_at: datetime,
        deployment_status: str = "InProgress",
    ):
        self.deployment_type = deployment_type
        self.update_at_datetime = updated_at
        self.deployment_status = deployment_status

    def to_dict(self) -> Dict[str, Any]:
        return {
            "DeploymentStatus": self.deployment_status,
            "DeploymentType": self.deployment_type,
            "UpdatedAt": iso_8601_datetime_with_milliseconds(self.update_at_datetime),
        }


class GreengrassBackend(BaseBackend):
    def __init__(self, region_name: str, account_id: str):
        super().__init__(region_name, account_id)
        self.groups: Dict[str, FakeGroup] = OrderedDict()
        self.group_role_associations: Dict[str, FakeAssociatedRole] = OrderedDict()
        self.group_versions: Dict[str, Dict[str, FakeGroupVersion]] = OrderedDict()
        self.core_definitions: Dict[str, FakeCoreDefinition] = OrderedDict()
        self.core_definition_versions: Dict[
            str, Dict[str, FakeCoreDefinitionVersion]
        ] = OrderedDict()
        self.device_definitions: Dict[str, FakeDeviceDefinition] = OrderedDict()
        self.device_definition_versions: Dict[
            str, Dict[str, FakeDeviceDefinitionVersion]
        ] = OrderedDict()
        self.function_definitions: Dict[str, FakeFunctionDefinition] = OrderedDict()
        self.function_definition_versions: Dict[
            str, Dict[str, FakeFunctionDefinitionVersion]
        ] = OrderedDict()
        self.resource_definitions: Dict[str, FakeResourceDefinition] = OrderedDict()
        self.resource_definition_versions: Dict[
            str, Dict[str, FakeResourceDefinitionVersion]
        ] = OrderedDict()
        self.subscription_definitions: Dict[str, FakeSubscriptionDefinition] = (
            OrderedDict()
        )
        self.subscription_definition_versions: Dict[
            str, Dict[str, FakeSubscriptionDefinitionVersion]
        ] = OrderedDict()
        self.deployments: Dict[str, FakeDeployment] = OrderedDict()

    def create_core_definition(
        self, name: str, initial_version: Dict[str, Any]
    ) -> FakeCoreDefinition:
        core_definition = FakeCoreDefinition(self.account_id, self.region_name, name)
        self.core_definitions[core_definition.id] = core_definition
        self.create_core_definition_version(
            core_definition.id, initial_version["Cores"]
        )
        return core_definition

    def list_core_definitions(self) -> Iterable[FakeCoreDefinition]:
        return self.core_definitions.values()

    def get_core_definition(self, core_definition_id: str) -> FakeCoreDefinition:
        if core_definition_id not in self.core_definitions:
            raise IdNotFoundException("That Core List Definition does not exist")
        return self.core_definitions[core_definition_id]

    def delete_core_definition(self, core_definition_id: str) -> None:
        if core_definition_id not in self.core_definitions:
            raise IdNotFoundException("That cores definition does not exist.")
        del self.core_definitions[core_definition_id]
        del self.core_definition_versions[core_definition_id]

    def update_core_definition(self, core_definition_id: str, name: str) -> None:
        if name == "":
            raise InvalidContainerDefinitionException(
                "Input does not contain any attributes to be updated"
            )
        if core_definition_id not in self.core_definitions:
            raise IdNotFoundException("That cores definition does not exist.")
        self.core_definitions[core_definition_id].name = name

    def create_core_definition_version(
        self, core_definition_id: str, cores: List[Dict[str, Any]]
    ) -> FakeCoreDefinitionVersion:
        definition = {"Cores": cores}
        core_def_ver = FakeCoreDefinitionVersion(
            self.account_id, self.region_name, core_definition_id, definition
        )
        core_def_vers = self.core_definition_versions.get(
            core_def_ver.core_definition_id, {}
        )
        core_def_vers[core_def_ver.version] = core_def_ver
        self.core_definition_versions[core_def_ver.core_definition_id] = core_def_vers

        self.core_definitions[core_definition_id].latest_version = core_def_ver.version
        self.core_definitions[core_definition_id].latest_version_arn = core_def_ver.arn

        return core_def_ver

    def list_core_definition_versions(
        self, core_definition_id: str
    ) -> Iterable[FakeCoreDefinitionVersion]:
        if core_definition_id not in self.core_definitions:
            raise IdNotFoundException("That cores definition does not exist.")
        return self.core_definition_versions[core_definition_id].values()

    def get_core_definition_version(
        self, core_definition_id: str, core_definition_version_id: str
    ) -> FakeCoreDefinitionVersion:
        if core_definition_id not in self.core_definitions:
            raise IdNotFoundException("That cores definition does not exist.")

        if (
            core_definition_version_id
            not in self.core_definition_versions[core_definition_id]
        ):
            raise VersionNotFoundException(
                f"Version {core_definition_version_id} of Core List Definition {core_definition_id} does not exist."
            )

        return self.core_definition_versions[core_definition_id][
            core_definition_version_id
        ]

    def create_device_definition(
        self, name: str, initial_version: Dict[str, Any]
    ) -> FakeDeviceDefinition:
        device_def = FakeDeviceDefinition(
            self.account_id, self.region_name, name, initial_version
        )
        self.device_definitions[device_def.id] = device_def
        init_ver = device_def.initial_version
        init_device_def = init_ver.get("Devices", {})
        self.create_device_definition_version(device_def.id, init_device_def)

        return device_def

    def list_device_definitions(self) -> Iterable[FakeDeviceDefinition]:
        return self.device_definitions.values()

    def create_device_definition_version(
        self, device_definition_id: str, devices: List[Dict[str, Any]]
    ) -> FakeDeviceDefinitionVersion:
        if device_definition_id not in self.device_definitions:
            raise IdNotFoundException("That devices definition does not exist.")

        device_ver = FakeDeviceDefinitionVersion(
            self.account_id, self.region_name, device_definition_id, devices
        )
        device_vers = self.device_definition_versions.get(
            device_ver.device_definition_id, {}
        )
        device_vers[device_ver.version] = device_ver
        self.device_definition_versions[device_ver.device_definition_id] = device_vers
        self.device_definitions[
            device_definition_id
        ].latest_version = device_ver.version
        self.device_definitions[
            device_definition_id
        ].latest_version_arn = device_ver.arn

        return device_ver

    def list_device_definition_versions(
        self, device_definition_id: str
    ) -> Iterable[FakeDeviceDefinitionVersion]:
        if device_definition_id not in self.device_definitions:
            raise IdNotFoundException("That devices definition does not exist.")
        return self.device_definition_versions[device_definition_id].values()

    def get_device_definition(self, device_definition_id: str) -> FakeDeviceDefinition:
        if device_definition_id not in self.device_definitions:
            raise IdNotFoundException("That Device List Definition does not exist.")
        return self.device_definitions[device_definition_id]

    def delete_device_definition(self, device_definition_id: str) -> None:
        if device_definition_id not in self.device_definitions:
            raise IdNotFoundException("That devices definition does not exist.")
        del self.device_definitions[device_definition_id]
        del self.device_definition_versions[device_definition_id]

    def update_device_definition(self, device_definition_id: str, name: str) -> None:
        if name == "":
            raise InvalidContainerDefinitionException(
                "Input does not contain any attributes to be updated"
            )
        if device_definition_id not in self.device_definitions:
            raise IdNotFoundException("That devices definition does not exist.")
        self.device_definitions[device_definition_id].name = name

    def get_device_definition_version(
        self, device_definition_id: str, device_definition_version_id: str
    ) -> FakeDeviceDefinitionVersion:
        if device_definition_id not in self.device_definitions:
            raise IdNotFoundException("That devices definition does not exist.")

        if (
            device_definition_version_id
            not in self.device_definition_versions[device_definition_id]
        ):
            raise VersionNotFoundException(
                f"Version {device_definition_version_id} of Device List Definition {device_definition_id} does not exist."
            )

        return self.device_definition_versions[device_definition_id][
            device_definition_version_id
        ]

    def create_resource_definition(
        self, name: str, initial_version: Dict[str, Any]
    ) -> FakeResourceDefinition:
        resources = initial_version.get("Resources", [])
        GreengrassBackend._validate_resources(resources)

        resource_def = FakeResourceDefinition(
            self.account_id, self.region_name, name, initial_version
        )
        self.resource_definitions[resource_def.id] = resource_def
        init_ver = resource_def.initial_version
        resources = init_ver.get("Resources", {})
        self.create_resource_definition_version(resource_def.id, resources)

        return resource_def

    def list_resource_definitions(self) -> Iterable[FakeResourceDefinition]:
        return self.resource_definitions.values()

    def get_resource_definition(
        self, resource_definition_id: str
    ) -> FakeResourceDefinition:
        if resource_definition_id not in self.resource_definitions:
            raise IdNotFoundException("That Resource List Definition does not exist.")
        return self.resource_definitions[resource_definition_id]

    def delete_resource_definition(self, resource_definition_id: str) -> None:
        if resource_definition_id not in self.resource_definitions:
            raise IdNotFoundException("That resources definition does not exist.")
        del self.resource_definitions[resource_definition_id]
        del self.resource_definition_versions[resource_definition_id]

    def update_resource_definition(
        self, resource_definition_id: str, name: str
    ) -> None:
        if name == "":
            raise InvalidInputException("Invalid resource name.")
        if resource_definition_id not in self.resource_definitions:
            raise IdNotFoundException("That resources definition does not exist.")
        self.resource_definitions[resource_definition_id].name = name

    def create_resource_definition_version(
        self, resource_definition_id: str, resources: List[Dict[str, Any]]
    ) -> FakeResourceDefinitionVersion:
        if resource_definition_id not in self.resource_definitions:
            raise IdNotFoundException("That resource definition does not exist.")

        GreengrassBackend._validate_resources(resources)

        resource_def_ver = FakeResourceDefinitionVersion(
            self.account_id, self.region_name, resource_definition_id, resources
        )

        resources_ver = self.resource_definition_versions.get(
            resource_def_ver.resource_definition_id, {}
        )
        resources_ver[resource_def_ver.version] = resource_def_ver
        self.resource_definition_versions[resource_def_ver.resource_definition_id] = (
            resources_ver
        )

        self.resource_definitions[
            resource_definition_id
        ].latest_version = resource_def_ver.version

        self.resource_definitions[
            resource_definition_id
        ].latest_version_arn = resource_def_ver.arn

        return resource_def_ver

    def list_resource_definition_versions(
        self, resource_definition_id: str
    ) -> Iterable[FakeResourceDefinitionVersion]:
        if resource_definition_id not in self.resource_definition_versions:
            raise IdNotFoundException("That resources definition does not exist.")

        return self.resource_definition_versions[resource_definition_id].values()

    def get_resource_definition_version(
        self, resource_definition_id: str, resource_definition_version_id: str
    ) -> FakeResourceDefinitionVersion:
        if resource_definition_id not in self.resource_definition_versions:
            raise IdNotFoundException("That resources definition does not exist.")

        if (
            resource_definition_version_id
            not in self.resource_definition_versions[resource_definition_id]
        ):
            raise VersionNotFoundException(
                f"Version {resource_definition_version_id} of Resource List Definition {resource_definition_id} does not exist."
            )

        return self.resource_definition_versions[resource_definition_id][
            resource_definition_version_id
        ]

    @staticmethod
    def _validate_resources(resources: List[Dict[str, Any]]) -> None:  # type: ignore[misc]
        for resource in resources:
            volume_source_path = (
                resource.get("ResourceDataContainer", {})
                .get("LocalVolumeResourceData", {})
                .get("SourcePath", "")
            )
            if volume_source_path == "/sys" or volume_source_path.startswith("/sys/"):
                raise GreengrassClientError(
                    "400",
                    "The resources definition is invalid. (ErrorDetails: [Accessing /sys is prohibited])",
                )

            local_device_resource_data = resource.get("ResourceDataContainer", {}).get(
                "LocalDeviceResourceData", {}
            )
            if local_device_resource_data:
                device_source_path = local_device_resource_data["SourcePath"]
                if not device_source_path.startswith("/dev"):
                    raise GreengrassClientError(
                        "400",
                        f"The resources definition is invalid. (ErrorDetails: [Device resource path should begin with "
                        "/dev"
                        f", but got: {device_source_path}])",
                    )

    def create_function_definition(
        self, name: str, initial_version: Dict[str, Any]
    ) -> FakeFunctionDefinition:
        func_def = FakeFunctionDefinition(
            self.account_id, self.region_name, name, initial_version
        )
        self.function_definitions[func_def.id] = func_def
        init_ver = func_def.initial_version
        init_func_def = init_ver.get("Functions", {})
        init_config = init_ver.get("DefaultConfig", {})
        self.create_function_definition_version(func_def.id, init_func_def, init_config)

        return func_def

    def list_function_definitions(self) -> List[FakeFunctionDefinition]:
        return list(self.function_definitions.values())

    def get_function_definition(
        self, function_definition_id: str
    ) -> FakeFunctionDefinition:
        if function_definition_id not in self.function_definitions:
            raise IdNotFoundException("That Lambda List Definition does not exist.")
        return self.function_definitions[function_definition_id]

    def delete_function_definition(self, function_definition_id: str) -> None:
        if function_definition_id not in self.function_definitions:
            raise IdNotFoundException("That lambdas definition does not exist.")
        del self.function_definitions[function_definition_id]
        del self.function_definition_versions[function_definition_id]

    def update_function_definition(
        self, function_definition_id: str, name: str
    ) -> None:
        if name == "":
            raise InvalidContainerDefinitionException(
                "Input does not contain any attributes to be updated"
            )
        if function_definition_id not in self.function_definitions:
            raise IdNotFoundException("That lambdas definition does not exist.")
        self.function_definitions[function_definition_id].name = name

    def create_function_definition_version(
        self,
        function_definition_id: str,
        functions: List[Dict[str, Any]],
        default_config: Dict[str, Any],
    ) -> FakeFunctionDefinitionVersion:
        if function_definition_id not in self.function_definitions:
            raise IdNotFoundException("That lambdas does not exist.")

        func_ver = FakeFunctionDefinitionVersion(
            self.account_id,
            self.region_name,
            function_definition_id,
            functions,
            default_config,
        )
        func_vers = self.function_definition_versions.get(
            func_ver.function_definition_id, {}
        )
        func_vers[func_ver.version] = func_ver
        self.function_definition_versions[func_ver.function_definition_id] = func_vers
        self.function_definitions[
            function_definition_id
        ].latest_version = func_ver.version
        self.function_definitions[
            function_definition_id
        ].latest_version_arn = func_ver.arn

        return func_ver

    def list_function_definition_versions(
        self, function_definition_id: str
    ) -> Dict[str, FakeFunctionDefinitionVersion]:
        if function_definition_id not in self.function_definition_versions:
            raise IdNotFoundException("That lambdas definition does not exist.")
        return self.function_definition_versions[function_definition_id]

    def get_function_definition_version(
        self, function_definition_id: str, function_definition_version_id: str
    ) -> FakeFunctionDefinitionVersion:
        if function_definition_id not in self.function_definition_versions:
            raise IdNotFoundException("That lambdas definition does not exist.")

        if (
            function_definition_version_id
            not in self.function_definition_versions[function_definition_id]
        ):
            raise IdNotFoundException(
                f"Version {function_definition_version_id} of Lambda List Definition {function_definition_id} does not exist."
            )

        return self.function_definition_versions[function_definition_id][
            function_definition_version_id
        ]

    @staticmethod
    def _is_valid_subscription_target_or_source(target_or_source: str) -> bool:
        if target_or_source in ["cloud", "GGShadowService"]:
            return True

        if re.match(
            r"^arn:aws:iot:[a-zA-Z0-9-]+:[0-9]{12}:thing/[a-zA-Z0-9-]+$",
            target_or_source,
        ):
            return True

        if re.match(
            r"^arn:aws:lambda:[a-zA-Z0-9-]+:[0-9]{12}:function:[a-zA-Z0-9-_]+:[a-zA-Z0-9-_]+$",
            target_or_source,
        ):
            return True

        return False

    @staticmethod
    def _validate_subscription_target_or_source(  # type: ignore[misc]
        subscriptions: List[Dict[str, Any]],
    ) -> None:
        target_errors: List[str] = []
        source_errors: List[str] = []

        for subscription in subscriptions:
            subscription_id = subscription["Id"]
            source = subscription["Source"]
            target = subscription["Target"]

            if not GreengrassBackend._is_valid_subscription_target_or_source(source):
                target_errors.append(
                    f"Subscription source is invalid. ID is '{subscription_id}' and Source is '{source}'"
                )

            if not GreengrassBackend._is_valid_subscription_target_or_source(target):
                target_errors.append(
                    f"Subscription target is invalid. ID is '{subscription_id}' and Target is '{target}'"
                )

        if source_errors:
            error_msg = ", ".join(source_errors)
            raise GreengrassClientError(
                "400",
                f"The subscriptions definition is invalid or corrupted. (ErrorDetails: [{error_msg}])",
            )

        if target_errors:
            error_msg = ", ".join(target_errors)
            raise GreengrassClientError(
                "400",
                f"The subscriptions definition is invalid or corrupted. (ErrorDetails: [{error_msg}])",
            )

    def create_subscription_definition(
        self, name: str, initial_version: Dict[str, Any]
    ) -> FakeSubscriptionDefinition:
        GreengrassBackend._validate_subscription_target_or_source(
            initial_version["Subscriptions"]
        )

        sub_def = FakeSubscriptionDefinition(
            self.account_id, self.region_name, name, initial_version
        )
        self.subscription_definitions[sub_def.id] = sub_def
        init_ver = sub_def.initial_version
        subscriptions = init_ver.get("Subscriptions", {})
        sub_def_ver = self.create_subscription_definition_version(
            sub_def.id, subscriptions
        )

        sub_def.latest_version = sub_def_ver.version
        sub_def.latest_version_arn = sub_def_ver.arn
        return sub_def

    def list_subscription_definitions(self) -> List[FakeSubscriptionDefinition]:
        return list(self.subscription_definitions.values())

    def get_subscription_definition(
        self, subscription_definition_id: str
    ) -> FakeSubscriptionDefinition:
        if subscription_definition_id not in self.subscription_definitions:
            raise IdNotFoundException(
                "That Subscription List Definition does not exist."
            )
        return self.subscription_definitions[subscription_definition_id]

    def delete_subscription_definition(self, subscription_definition_id: str) -> None:
        if subscription_definition_id not in self.subscription_definitions:
            raise IdNotFoundException("That subscriptions definition does not exist.")
        del self.subscription_definitions[subscription_definition_id]
        del self.subscription_definition_versions[subscription_definition_id]

    def update_subscription_definition(
        self, subscription_definition_id: str, name: str
    ) -> None:
        if name == "":
            raise InvalidContainerDefinitionException(
                "Input does not contain any attributes to be updated"
            )
        if subscription_definition_id not in self.subscription_definitions:
            raise IdNotFoundException("That subscriptions definition does not exist.")
        self.subscription_definitions[subscription_definition_id].name = name

    def create_subscription_definition_version(
        self, subscription_definition_id: str, subscriptions: List[Dict[str, Any]]
    ) -> FakeSubscriptionDefinitionVersion:
        GreengrassBackend._validate_subscription_target_or_source(subscriptions)

        if subscription_definition_id not in self.subscription_definitions:
            raise IdNotFoundException("That subscriptions does not exist.")

        sub_def_ver = FakeSubscriptionDefinitionVersion(
            self.account_id, self.region_name, subscription_definition_id, subscriptions
        )

        sub_vers = self.subscription_definition_versions.get(
            subscription_definition_id, {}
        )
        sub_vers[sub_def_ver.version] = sub_def_ver
        self.subscription_definition_versions[subscription_definition_id] = sub_vers

        return sub_def_ver

    def list_subscription_definition_versions(
        self, subscription_definition_id: str
    ) -> Dict[str, FakeSubscriptionDefinitionVersion]:
        if subscription_definition_id not in self.subscription_definition_versions:
            raise IdNotFoundException("That subscriptions definition does not exist.")
        return self.subscription_definition_versions[subscription_definition_id]

    def get_subscription_definition_version(
        self, subscription_definition_id: str, subscription_definition_version_id: str
    ) -> FakeSubscriptionDefinitionVersion:
        if subscription_definition_id not in self.subscription_definitions:
            raise IdNotFoundException("That subscriptions definition does not exist.")

        if (
            subscription_definition_version_id
            not in self.subscription_definition_versions[subscription_definition_id]
        ):
            raise VersionNotFoundException(
                f"Version {subscription_definition_version_id} of Subscription List Definition {subscription_definition_id} does not exist."
            )

        return self.subscription_definition_versions[subscription_definition_id][
            subscription_definition_version_id
        ]

    def create_group(self, name: str, initial_version: Dict[str, Any]) -> FakeGroup:
        group = FakeGroup(self.account_id, self.region_name, name)
        self.groups[group.group_id] = group

        definitions = initial_version or {}
        core_definition_version_arn = definitions.get("CoreDefinitionVersionArn")
        device_definition_version_arn = definitions.get("DeviceDefinitionVersionArn")
        function_definition_version_arn = definitions.get(
            "FunctionDefinitionVersionArn"
        )
        resource_definition_version_arn = definitions.get(
            "ResourceDefinitionVersionArn"
        )
        subscription_definition_version_arn = definitions.get(
            "SubscriptionDefinitionVersionArn"
        )

        self.create_group_version(
            group.group_id,
            core_definition_version_arn=core_definition_version_arn,
            device_definition_version_arn=device_definition_version_arn,
            function_definition_version_arn=function_definition_version_arn,
            resource_definition_version_arn=resource_definition_version_arn,
            subscription_definition_version_arn=subscription_definition_version_arn,
        )

        return group

    def list_groups(self) -> List[FakeGroup]:
        return list(self.groups.values())

    def get_group(self, group_id: str) -> Optional[FakeGroup]:
        if group_id not in self.groups:
            raise IdNotFoundException("That Group Definition does not exist.")
        return self.groups.get(group_id)

    def delete_group(self, group_id: str) -> None:
        if group_id not in self.groups:
            # I don't know why, the error message is different between get_group and delete_group
            raise IdNotFoundException("That group definition does not exist.")
        del self.groups[group_id]
        del self.group_versions[group_id]

    def update_group(self, group_id: str, name: str) -> None:
        if name == "":
            raise InvalidContainerDefinitionException(
                "Input does not contain any attributes to be updated"
            )
        if group_id not in self.groups:
            raise IdNotFoundException("That group definition does not exist.")
        self.groups[group_id].name = name

    def create_group_version(
        self,
        group_id: str,
        core_definition_version_arn: Optional[str],
        device_definition_version_arn: Optional[str],
        function_definition_version_arn: Optional[str],
        resource_definition_version_arn: Optional[str],
        subscription_definition_version_arn: Optional[str],
    ) -> FakeGroupVersion:
        if group_id not in self.groups:
            raise IdNotFoundException("That group does not exist.")

        self._validate_group_version_definitions(
            core_definition_version_arn=core_definition_version_arn,
            device_definition_version_arn=device_definition_version_arn,
            function_definition_version_arn=function_definition_version_arn,
            resource_definition_version_arn=resource_definition_version_arn,
            subscription_definition_version_arn=subscription_definition_version_arn,
        )

        group_ver = FakeGroupVersion(
            self.account_id,
            self.region_name,
            group_id=group_id,
            core_definition_version_arn=core_definition_version_arn,
            device_definition_version_arn=device_definition_version_arn,
            function_definition_version_arn=function_definition_version_arn,
            resource_definition_version_arn=resource_definition_version_arn,
            subscription_definition_version_arn=subscription_definition_version_arn,
        )
        group_vers = self.group_versions.get(group_ver.group_id, {})
        group_vers[group_ver.version] = group_ver
        self.group_versions[group_ver.group_id] = group_vers
        self.groups[group_id].latest_version_arn = group_ver.arn
        self.groups[group_id].latest_version = group_ver.version
        return group_ver

    def _validate_group_version_definitions(
        self,
        core_definition_version_arn: Optional[str] = None,
        device_definition_version_arn: Optional[str] = None,
        function_definition_version_arn: Optional[str] = None,
        resource_definition_version_arn: Optional[str] = None,
        subscription_definition_version_arn: Optional[str] = None,
    ) -> None:
        def _is_valid_def_ver_arn(
            definition_version_arn: Optional[str], kind: str = "cores"
        ) -> bool:
            if definition_version_arn is None:
                return True

            if kind == "cores":
                versions: Any = self.core_definition_versions
            elif kind == "devices":
                versions = self.device_definition_versions
            elif kind == "functions":
                versions = self.function_definition_versions
            elif kind == "resources":
                versions = self.resource_definition_versions
            elif kind == "subscriptions":
                versions = self.subscription_definition_versions
            else:
                raise Exception("invalid args")

            arn_regex = (
                r"^arn:aws:greengrass:[a-zA-Z0-9-]+:[0-9]{12}:greengrass/definition/"
                + kind
                + r"/[a-z0-9-]{36}/versions/[a-z0-9-]{36}$"
            )

            if not re.match(arn_regex, definition_version_arn):
                return False

            definition_id = definition_version_arn.split("/")[-3]

            if definition_id not in versions:
                return False

            definition_version_id = definition_version_arn.split("/")[-1]
            if definition_version_id not in versions[definition_id]:
                return False

            if (
                versions[definition_id][definition_version_id].arn
                != definition_version_arn
            ):
                return False

            return True

        errors = []

        if not _is_valid_def_ver_arn(core_definition_version_arn, kind="cores"):
            errors.append("Cores definition reference does not exist")

        if not _is_valid_def_ver_arn(function_definition_version_arn, kind="functions"):
            errors.append("Lambda definition reference does not exist")

        if not _is_valid_def_ver_arn(resource_definition_version_arn, kind="resources"):
            errors.append("Resource definition reference does not exist")

        if not _is_valid_def_ver_arn(device_definition_version_arn, kind="devices"):
            errors.append("Devices definition reference does not exist")

        if not _is_valid_def_ver_arn(
            subscription_definition_version_arn, kind="subscriptions"
        ):
            errors.append("Subscription definition reference does not exist")

        if errors:
            error_details = ", ".join(errors)
            raise GreengrassClientError(
                "400",
                f"The group is invalid or corrupted. (ErrorDetails: [{error_details}])",
            )

    def list_group_versions(self, group_id: str) -> List[FakeGroupVersion]:
        if group_id not in self.group_versions:
            raise IdNotFoundException("That group definition does not exist.")
        return list(self.group_versions[group_id].values())

    def get_group_version(
        self, group_id: str, group_version_id: str
    ) -> FakeGroupVersion:
        if group_id not in self.group_versions:
            raise IdNotFoundException("That group definition does not exist.")

        if group_version_id not in self.group_versions[group_id]:
            raise VersionNotFoundException(
                f"Version {group_version_id} of Group Definition {group_id} does not exist."
            )

        return self.group_versions[group_id][group_version_id]

    def create_deployment(
        self,
        group_id: str,
        group_version_id: str,
        deployment_type: str,
        deployment_id: Optional[str] = None,
    ) -> FakeDeployment:
        deployment_types = (
            "NewDeployment",
            "Redeployment",
            "ResetDeployment",
            "ForceResetDeployment",
        )
        if deployment_type not in deployment_types:
            raise InvalidInputException(
                f"That deployment type is not valid.  Please specify one of the following types: {{{','.join(deployment_types)}}}."
            )
        if deployment_type == "Redeployment":
            if deployment_id is None:
                raise InvalidInputException(
                    "Your request is missing the following required parameter(s): {DeploymentId}."
                )
            if deployment_id not in self.deployments:
                raise InvalidInputException(
                    f"Deployment ID '{deployment_id}' is invalid."
                )

        if group_id not in self.groups:
            raise ResourceNotFoundException("That group definition does not exist.")

        if group_version_id not in self.group_versions[group_id]:
            raise ResourceNotFoundException(
                f"Version {group_version_id} of Group Definition {group_id} does not exist."
            )

        if (
            self.group_versions[group_id][group_version_id].core_definition_version_arn
            is None
        ):
            err = {
                "ErrorDetails": [
                    {
                        "DetailedErrorCode": "GG-303",
                        "DetailedErrorMessage": "You need a Greengrass Core in this Group before you can deploy.",
                    }
                ]
            }

            raise MissingCoreException(json.dumps(err))
        group_version_arn = self.group_versions[group_id][group_version_id].arn
        deployment = FakeDeployment(
            self.account_id,
            self.region_name,
            group_id,
            group_version_arn,
            deployment_type,
        )
        self.deployments[deployment.id] = deployment
        return deployment

    def list_deployments(self, group_id: str) -> List[FakeDeployment]:
        # ListDeployments API does not check specified group is exists
        return [
            deployment
            for deployment in self.deployments.values()
            if deployment.group_id == group_id
        ]

    def get_deployment_status(
        self, group_id: str, deployment_id: str
    ) -> FakeDeploymentStatus:
        if deployment_id not in self.deployments:
            raise InvalidInputException(f"Deployment '{deployment_id}' does not exist.")

        deployment = self.deployments[deployment_id]

        if deployment.group_id != group_id:
            raise InvalidInputException(f"Deployment '{deployment_id}' does not exist.")

        return FakeDeploymentStatus(
            deployment.deployment_type,
            deployment.update_at_datetime,
            deployment.deployment_status,
        )

    def reset_deployments(self, group_id: str, force: bool = False) -> FakeDeployment:
        if group_id not in self.groups:
            raise ResourceNotFoundException("That Group Definition does not exist.")

        deployment_type = "ForceResetDeployment"
        if not force:
            deployments = list(self.deployments.values())
            reset_error_msg = (
                f"Group id: {group_id} has not been deployed or has already been reset."
            )
            if not deployments:
                raise ResourceNotFoundException(reset_error_msg)
            if deployments[-1].deployment_type not in ["NewDeployment", "Redeployment"]:
                raise ResourceNotFoundException(reset_error_msg)
            deployment_type = "ResetDeployment"

        group = self.groups[group_id]
        deployment = FakeDeployment(
            self.account_id, self.region_name, group_id, group.arn, deployment_type
        )
        self.deployments[deployment.id] = deployment
        return deployment

    def associate_role_to_group(
        self, group_id: str, role_arn: str
    ) -> FakeAssociatedRole:
        # I don't know why, AssociateRoleToGroup does not check specified group is exists
        # So, this API allows any group id such as "a"

        associated_role = FakeAssociatedRole(role_arn)
        self.group_role_associations[group_id] = associated_role
        return associated_role

    def get_associated_role(self, group_id: str) -> FakeAssociatedRole:
        if group_id not in self.group_role_associations:
            raise GreengrassClientError(
                "404", "You need to attach an IAM role to this deployment group."
            )

        return self.group_role_associations[group_id]

    def disassociate_role_from_group(self, group_id: str) -> None:
        if group_id not in self.group_role_associations:
            return
        del self.group_role_associations[group_id]


greengrass_backends = BackendDict(GreengrassBackend, "greengrass")
