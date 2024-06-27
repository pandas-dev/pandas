import base64
import json
import uuid
from datetime import datetime, timedelta
from typing import Any, Dict, List, Optional, Union

from dateutil.tz import tzutc

from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import BaseModel
from moto.moto_api._internal.managed_state_model import ManagedState
from moto.panorama.utils import (
    arn_formatter,
    deep_convert_datetime_to_isoformat,
    generate_package_id,
    hash_device_name,
)
from moto.utilities.paginator import paginate

from .exceptions import (
    ValidationError,
)

PAGINATION_MODEL = {
    "list_devices": {
        "input_token": "next_token",
        "limit_key": "max_results",
        "limit_default": 123,
        "unique_attribute": "device_id",
    },
    "list_nodes": {
        "input_token": "next_token",
        "limit_key": "max_results",
        "limit_default": 123,
        "unique_attribute": "package_id",
    },
}


class BaseObject(BaseModel):
    def camelCase(self, key: str) -> str:
        words = []
        for word in key.split("_"):
            words.append(word.title())
        return "".join(words)

    def update(self, details_json: str) -> None:
        details = json.loads(details_json)
        for k in details.keys():
            setattr(self, k, details[k])

    def gen_response_object(self) -> Dict[str, Any]:
        response_object: Dict[str, Any] = dict()
        for key, value in self.__dict__.items():
            if "_" in key:
                response_object[self.camelCase(key)] = value
            else:
                response_object[key[0].upper() + key[1:]] = value
        return response_object

    def response_object(self) -> Dict[str, Any]:
        return self.gen_response_object()


class Device(BaseObject):
    def __init__(
        self,
        account_id: str,
        region_name: str,
        description: Optional[str],
        name: str,
        network_configuration: Optional[Dict[str, Any]],
        tags: Optional[Dict[str, str]],
    ) -> None:
        # ManagedState is a class that helps us manage the state of a resource.
        # A Panorama Device has a lot of different states that has their own lifecycle.
        # To make all this states in the same object, and avoid changing ManagedState,
        # we use a ManagedState for each state and we manage them in the Device class.
        # Each ManagedState has a name composed with attribute name and device name to make subscription easier.
        self.__device_aggregated_status_manager = ManagedState(
            model_name=f"panorama::device_{name}_aggregated_status",
            transitions=[
                ("NOT-A-STATUS", "AWAITING_PROVISIONING"),
                ("AWAITING_PROVISIONING", "PENDING"),
                ("PENDING", "ONLINE"),
            ],
        )
        self.__device_provisioning_status_manager = ManagedState(
            model_name=f"panorama::device_{name}_provisioning_status",
            transitions=[
                ("NOT-A-STATUS", "AWAITING_PROVISIONING"),
                ("AWAITING_PROVISIONING", "PENDING"),
                ("PENDING", "SUCCEEDED"),
            ],
        )
        self.account_id = account_id
        self.region_name = region_name
        self.description = description
        self.name = name
        self.network_configuration = network_configuration
        self.tags = tags

        self.certificates = base64.b64encode("certificate".encode("utf-8")).decode(
            "utf-8"
        )
        self.arn = arn_formatter("device", self.name, self.account_id, self.region_name)
        self.device_id = f"device-{hash_device_name(name)}"
        self.iot_thing_name = ""

        self.alternate_softwares = [
            {"Version": "0.2.1"},
        ]
        self.brand: str = "AWS_PANORAMA"  # AWS_PANORAMA | LENOVO
        self.created_time = datetime.now(tzutc())
        self.last_updated_time = datetime.now(tzutc())
        self.current_networking_status = {
            "Ethernet0Status": {
                "ConnectionStatus": "CONNECTED",
                "HwAddress": "8C:0F:5F:60:F5:C4",
                "IpAddress": "192.168.1.300/24",
            },
            "Ethernet1Status": {
                "ConnectionStatus": "NOT_CONNECTED",
                "HwAddress": "8C:0F:6F:60:F4:F1",
                "IpAddress": "--",
            },
            "LastUpdatedTime": datetime.now(tzutc()),
            "NtpStatus": {
                "ConnectionStatus": "CONNECTED",
                "IpAddress": "91.224.149.41:123",
                "NtpServerName": "0.pool.ntp.org",
            },
        }
        self.current_software = "6.2.1"
        self.device_connection_status: str = "ONLINE"  # "ONLINE"|"OFFLINE"|"AWAITING_CREDENTIALS"|"NOT_AVAILABLE"|"ERROR"
        self.latest_device_job = {"JobType": "REBOOT", "Status": "COMPLETED"}
        self.latest_software = "6.2.1"
        self.lease_expiration_time = datetime.now(tzutc()) + timedelta(days=5)
        self.serial_number = "GAD81E29013274749"
        self.type: str = "PANORAMA_APPLIANCE"  # "PANORAMA_APPLIANCE_DEVELOPER_KIT", "PANORAMA_APPLIANCE"

    @property
    def device_aggregated_status(self) -> str:
        self.__device_aggregated_status_manager.advance()
        return self.__device_aggregated_status_manager.status  # type: ignore[return-value]

    @property
    def provisioning_status(self) -> str:
        self.__device_provisioning_status_manager.advance()
        return self.__device_provisioning_status_manager.status  # type: ignore[return-value]

    def response_object(self) -> Dict[str, Any]:
        response_object = super().gen_response_object()
        response_object = deep_convert_datetime_to_isoformat(response_object)
        static_response_fields = [
            "AlternateSoftwares",
            "Arn",
            "Brand",
            "CreatedTime",
            "CurrentNetworkingStatus",
            "CurrentSoftware",
            "Description",
            "DeviceConnectionStatus",
            "DeviceId",
            "LatestAlternateSoftware",
            "LatestDeviceJob",
            "LatestSoftware",
            "LeaseExpirationTime",
            "Name",
            "NetworkConfiguration",
            "SerialNumber",
            "Tags",
            "Type",
        ]
        return {
            **{
                k: v
                for k, v in response_object.items()
                if v is not None and k in static_response_fields
            },
            **{
                "DeviceAggregatedStatus": self.device_aggregated_status,
                "ProvisioningStatus": self.provisioning_status,
            },
        }

    def response_listed(self) -> Dict[str, Any]:
        response_object = super().gen_response_object()
        response_object = deep_convert_datetime_to_isoformat(response_object)
        static_response_fields = [
            "Brand",
            "CreatedTime",
            "CurrentSoftware",
            "Description",
            "DeviceId",
            "LastUpdatedTime",
            "LatestDeviceJob",
            "LeaseExpirationTime",
            "Name",
            "Tags",
            "Type",
        ]
        return {
            **{
                k: v
                for k, v in response_object.items()
                if v is not None and k in static_response_fields
            },
            **{
                "DeviceAggregatedStatus": self.device_aggregated_status,
                "ProvisioningStatus": self.provisioning_status,
            },
        }

    @property
    def response_provision(self) -> Dict[str, Union[str, bytes]]:
        return {
            "Arn": self.arn,
            "Certificates": self.certificates,
            "DeviceId": self.device_id,
            "IotThingName": self.iot_thing_name,
            "Status": self.provisioning_status,
        }

    @property
    def response_updated(self) -> Dict[str, str]:
        return {"DeviceId": self.device_id}

    @property
    def response_deleted(self) -> Dict[str, str]:
        return {"DeviceId": self.device_id}


class Package(BaseObject):
    def __init__(
        self,
        category: str,
        description: str,
        name: str,
        account_id: str,
        region_name: str,
        package_name: str,
        package_version: str,
    ):
        self.category = category
        self.description = description
        self.name = name
        now = datetime.now(tzutc())
        self.created_time = now
        self.last_updated_time = now
        self.package_name = package_name
        self.patch_version = generate_package_id(self.package_name)
        self.package_version = package_version
        self.output_package_name = f"{self.package_name}-{self.package_version}-{self.patch_version[:8]}-{self.name}"
        self.owner_account = account_id
        self.package_id = f"package-{hash_device_name(package_name)}"
        self.package_arn = arn_formatter(
            "package", self.package_id, account_id, region_name
        )

    def response_object(self) -> Dict[str, Any]:
        response_object = super().gen_response_object()
        response_object = deep_convert_datetime_to_isoformat(response_object)
        return response_object

    def response_listed(self) -> Dict[str, Any]:
        package_response = self.response_object()
        return package_response


class Node(BaseObject):
    def __init__(
        self,
        job_id: str,
        job_tags: List[Dict[str, Union[str, Dict[str, str]]]],
        node_description: str,
        node_name: str,
        output_package_name: str,
        output_package_version: str,
        template_parameters: Dict[str, str],
        template_type: str,
    ) -> None:
        self.job_id = job_id
        now = datetime.now(tzutc())
        self.created_time = now
        self.last_updated_time = now
        self.job_tags = job_tags
        self.node_description = node_name
        self.status = "PENDING"
        self.node_name = node_name
        self.output_package_name = output_package_name
        self.output_package_version = output_package_version
        self.template_parameters = self.protect_secrets(template_parameters)
        self.template_type = template_type

    def response_object(self) -> Dict[str, Any]:
        response_object = super().gen_response_object()
        response_object = deep_convert_datetime_to_isoformat(response_object)
        return response_object

    def response_created(self) -> Dict[str, str]:
        return {"JobId": self.job_id}

    def response_described(self) -> Dict[str, Any]:
        return self.response_object()

    @staticmethod
    def protect_secrets(template_parameters: Dict[str, str]) -> Dict[str, str]:
        for key in template_parameters.keys():
            if key.lower() == "password":
                template_parameters[key] = "SAVED_AS_SECRET"
            if key.lower() == "username":
                template_parameters[key] = "SAVED_AS_SECRET"
        return template_parameters


class PanoramaBackend(BaseBackend):
    def __init__(self, region_name: str, account_id: str):
        super().__init__(region_name, account_id)
        self.devices_memory: Dict[str, Device] = {}
        self.node_from_template_memory: Dict[str, Node] = {}
        self.nodes_memory: Dict[str, Package] = {}

    def provision_device(
        self,
        description: Optional[str],
        name: str,
        networking_configuration: Optional[Dict[str, Any]],
        tags: Optional[Dict[str, str]],
    ) -> Device:
        device_obj = Device(
            account_id=self.account_id,
            region_name=self.region_name,
            description=description,
            name=name,
            network_configuration=networking_configuration,
            tags=tags,
        )

        self.devices_memory[device_obj.device_id] = device_obj
        return device_obj

    def describe_device(self, device_id: str) -> Device:
        device = self.devices_memory.get(device_id)
        if device is None:
            raise ValidationError(f"Device {device_id} not found")
        return device

    @paginate(pagination_model=PAGINATION_MODEL)
    def list_devices(
        self,
        device_aggregated_status_filter: str,
        name_filter: str,
        sort_by: str,  # "DEVICE_ID", "CREATED_TIME", "NAME", "DEVICE_AGGREGATED_STATUS"
        sort_order: str,  # "ASCENDING", "DESCENDING"
    ) -> List[Device]:
        devices_list = list(
            filter(
                lambda x: (name_filter is None or x.name.startswith(name_filter))
                and (
                    device_aggregated_status_filter is None
                    or x.device_aggregated_status == device_aggregated_status_filter
                ),
                self.devices_memory.values(),
            )
        )
        devices_list = list(
            sorted(
                devices_list,
                key={
                    "DEVICE_ID": lambda x: x.device_id,
                    "CREATED_TIME": lambda x: x.created_time,
                    "NAME": lambda x: x.name,
                    "DEVICE_AGGREGATED_STATUS": lambda x: x.device_aggregated_status,
                    None: lambda x: x.created_time,
                }[sort_by],
                reverse=sort_order == "DESCENDING",
            )
        )
        return devices_list

    def update_device_metadata(self, device_id: str, description: str) -> Device:
        self.devices_memory[device_id].description = description
        return self.devices_memory[device_id]

    def delete_device(self, device_id: str) -> Device:
        return self.devices_memory.pop(device_id)

    def create_node_from_template_job(
        self,
        job_tags: List[Dict[str, Union[str, Dict[str, str]]]],
        node_description: str,
        node_name: str,
        output_package_name: str,
        output_package_version: str,
        template_parameters: Dict[str, str],
        template_type: str,
    ) -> Node:
        job_id = str(uuid.uuid4()).lower()
        self.node_from_template_memory[job_id] = Node(
            job_id=job_id,
            job_tags=job_tags,
            node_description=node_description,
            node_name=node_name,
            output_package_name=output_package_name,
            output_package_version=output_package_version,
            template_parameters=template_parameters,
            template_type=template_type,
        )
        if template_type == "RTSP_CAMERA_STREAM":
            package = Package(
                category="MEDIA_SOURCE",
                description=node_description,
                name=node_name,
                account_id=self.account_id,
                region_name=self.region_name,
                package_name=output_package_name,
                package_version=output_package_version,
            )
            self.nodes_memory[package.package_id] = package

        return self.node_from_template_memory[job_id]

    def describe_node_from_template_job(self, job_id: str) -> Node:
        return self.node_from_template_memory[job_id]

    @paginate(pagination_model=PAGINATION_MODEL)
    def list_nodes(self, category: str) -> List[Package]:
        category_nodes = list(
            filter(
                lambda x: x.category == category,
                self.nodes_memory.values(),
            )
        )
        return category_nodes


panorama_backends = BackendDict(
    PanoramaBackend,
    "panorama",
    False,
    additional_regions=[
        "us-east-1",
        "us-west-2",
        "ca-central-1",
        "eu-west-1",
        "ap-southeast-2",
        "ap-southeast-1",
    ],
)
