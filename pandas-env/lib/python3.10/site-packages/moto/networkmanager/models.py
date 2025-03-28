"""NetworkManagerBackend class with methods for supported APIs."""

from datetime import datetime, timezone
from typing import Any, Dict, List, Optional

from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import BaseModel
from moto.moto_api._internal import mock_random
from moto.utilities.paginator import paginate
from moto.utilities.utils import PARTITION_NAMES

from .exceptions import ResourceNotFound, ValidationError

PAGINATION_MODEL = {
    "describe_global_networks": {
        "input_token": "next_token",
        "limit_key": "max_results",
        "limit_default": 100,
        "unique_attribute": "global_network_arn",
    },
    "list_core_networks": {
        "input_token": "next_token",
        "limit_key": "max_results",
        "limit_default": 100,
        "unique_attribute": "core_network_arn",
    },
    "get_sites": {
        "input_token": "next_token",
        "limit_key": "max_results",
        "limit_default": 100,
        "unique_attribute": "site_arn",
    },
    "get_links": {
        "input_token": "next_token",
        "limit_key": "max_results",
        "limit_default": 100,
        "unique_attribute": "link_arn",
    },
    "get_devices": {
        "input_token": "next_token",
        "limit_key": "max_results",
        "limit_default": 100,
        "unique_attribute": "device_arn",
    },
}


class GlobalNetwork(BaseModel):
    def __init__(
        self,
        account_id: str,
        partition: str,
        description: Optional[str],
        tags: Optional[List[Dict[str, str]]],
    ):
        self.description = description
        self.tags = tags or []
        self.global_network_id = "global-network-" + "".join(
            mock_random.get_random_hex(18)
        )
        self.global_network_arn = f"arn:{partition}:networkmanager:{account_id}:global-network/{self.global_network_id}"
        self.created_at = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%S.%fZ")
        self.state = "PENDING"

    def to_dict(self) -> Dict[str, Any]:
        return {
            "GlobalNetworkId": self.global_network_id,
            "GlobalNetworkArn": self.global_network_arn,
            "Description": self.description,
            "Tags": self.tags,
            "State": self.state,
            "CreatedAt": self.created_at,
        }


class CoreNetwork(BaseModel):
    def __init__(
        self,
        account_id: str,
        partition: str,
        global_network_id: str,
        description: Optional[str],
        tags: Optional[List[Dict[str, str]]],
        policy_document: str,
        client_token: str,
    ):
        self.global_network_id = global_network_id
        self.description = description
        self.tags = tags or []
        self.policy_document = policy_document
        self.client_token = client_token
        self.core_network_id = "core-network-" + "".join(mock_random.get_random_hex(18))
        self.core_network_arn = f"arn:{partition}:networkmanager:{account_id}:core-network/{self.core_network_id}"

        self.created_at = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%S.%fZ")
        self.state = "PENDING"

    def to_dict(self) -> Dict[str, Any]:
        return {
            "CoreNetworkId": self.core_network_id,
            "CoreNetworkArn": self.core_network_arn,
            "GlobalNetworkId": self.global_network_id,
            "Description": self.description,
            "Tags": self.tags,
            "PolicyDocument": self.policy_document,
            "State": self.state,
            "CreatedAt": self.created_at,
        }


class Site(BaseModel):
    def __init__(
        self,
        account_id: str,
        partition: str,
        global_network_id: str,
        description: Optional[str],
        location: Optional[Dict[str, Any]],
        tags: Optional[List[Dict[str, str]]],
    ):
        self.global_network_id = global_network_id
        self.description = description
        self.location = location
        self.tags = tags or []
        self.site_id = "site-" + "".join(mock_random.get_random_hex(18))
        self.site_arn = (
            f"arn:{partition}:networkmanager:{account_id}:site/{self.site_id}"
        )
        self.created_at = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%S.%fZ")
        self.state = "PENDING"

    def to_dict(self) -> Dict[str, Any]:
        return {
            "SiteId": self.site_id,
            "SiteArn": self.site_arn,
            "GlobalNetworkId": self.global_network_id,
            "Description": self.description,
            "Location": self.location,
            "Tags": self.tags,
            "State": self.state,
            "CreatedAt": self.created_at,
        }


class Link(BaseModel):
    def __init__(
        self,
        account_id: str,
        partition: str,
        global_network_id: str,
        description: Optional[str],
        type: Optional[str],
        bandwidth: Dict[str, int],
        provider: Optional[str],
        site_id: str,
        tags: Optional[List[Dict[str, str]]],
    ):
        self.global_network_id = global_network_id
        self.description = description
        self.type = type
        self.bandwidth = bandwidth
        self.provider = provider
        self.site_id = site_id
        self.tags = tags or []
        self.link_id = "link-" + "".join(mock_random.get_random_hex(18))
        self.link_arn = (
            f"arn:{partition}:networkmanager:{account_id}:link/{self.link_id}"
        )
        self.created_at = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%S.%fZ")
        self.state = "PENDING"

    def to_dict(self) -> Dict[str, Any]:
        return {
            "LinkId": self.link_id,
            "LinkArn": self.link_arn,
            "GlobalNetworkId": self.global_network_id,
            "Description": self.description,
            "Type": self.type,
            "Bandwidth": self.bandwidth,
            "Provider": self.provider,
            "SiteId": self.site_id,
            "Tags": self.tags,
            "State": self.state,
            "CreatedAt": self.created_at,
        }


class Device(BaseModel):
    def __init__(
        self,
        account_id: str,
        partition: str,
        global_network_id: str,
        aws_location: Optional[Dict[str, str]],
        description: Optional[str],
        type: Optional[str],
        vendor: Optional[str],
        model: Optional[str],
        serial_number: Optional[str],
        location: Optional[Dict[str, str]],
        site_id: Optional[str],
        tags: Optional[List[Dict[str, str]]],
    ):
        self.global_network_id = global_network_id
        self.aws_location = aws_location
        self.description = description
        self.type = type
        self.vendor = vendor
        self.model = model
        self.serial_number = serial_number
        self.location = location
        self.site_id = site_id
        self.tags = tags or []
        self.device_id = "device-" + "".join(mock_random.get_random_hex(18))
        self.device_arn = (
            f"arn:{partition}:networkmanager:{account_id}:device/{self.device_id}"
        )
        self.created_at = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%S.%fZ")
        self.state = "PENDING"

    def to_dict(self) -> Dict[str, Any]:
        return {
            "DeviceId": self.device_id,
            "DeviceArn": self.device_arn,
            "GlobalNetworkId": self.global_network_id,
            "AWSLocation": self.aws_location,
            "Description": self.description,
            "Type": self.type,
            "Vendor": self.vendor,
            "Model": self.model,
            "SerialNumber": self.serial_number,
            "Location": self.location,
            "SiteId": self.site_id,
            "Tags": self.tags,
            "State": self.state,
            "CreatedAt": self.created_at,
        }


class NetworkManagerBackend(BaseBackend):
    """Implementation of NetworkManager APIs."""

    def __init__(self, region_name: str, account_id: str) -> None:
        super().__init__(region_name, account_id)
        self.global_networks: Dict[str, GlobalNetwork] = {}
        self.core_networks: Dict[str, CoreNetwork] = {}
        self.sites: Dict[str, Dict[str, Site]] = {}
        self.links: Dict[str, Dict[str, Link]] = {}
        self.devices: Dict[str, Dict[str, Device]] = {}

    def _get_resource_from_arn(self, arn: str) -> Any:
        resources_types: Dict[str, Dict[str, Any]] = {
            "core-network": self.core_networks,
            "global-network": self.global_networks,
            "site": self.sites,
            "link": self.links,
            "device": self.devices,
        }
        try:
            target_resource, target_name = arn.split(":")[-1].split("/")
            resources = resources_types.get(target_resource, {})
            # Flatten the nested dictionary stores
            if target_resource not in ["core-network", "global-network"]:
                resources = {
                    k: v
                    for inner_dict in resources.values()
                    for k, v in inner_dict.items()
                }
            return resources[target_name]
        except (KeyError, ValueError, AttributeError):
            raise ResourceNotFound(arn)

    def update_resource_state(self, resource_arn: str, state: str) -> None:
        # Acceptable states: PENDING, AVAILABLE, DELETING, UPDATING
        resource = self._get_resource_from_arn(resource_arn)
        resource.state = state

    def create_global_network(
        self,
        description: Optional[str],
        tags: Optional[List[Dict[str, str]]],
    ) -> GlobalNetwork:
        global_network = GlobalNetwork(
            description=description,
            tags=tags,
            account_id=self.account_id,
            partition=self.partition,
        )
        gnw_id = global_network.global_network_id
        self.global_networks[gnw_id] = global_network
        # Create empty dict for resources
        self.sites[gnw_id] = {}
        self.links[gnw_id] = {}
        self.devices[gnw_id] = {}
        return global_network

    def create_core_network(
        self,
        global_network_id: str,
        description: Optional[str],
        tags: Optional[List[Dict[str, str]]],
        policy_document: str,
        client_token: str,
    ) -> CoreNetwork:
        # check if global network exists
        if global_network_id not in self.global_networks:
            raise ResourceNotFound(global_network_id)

        core_network = CoreNetwork(
            global_network_id=global_network_id,
            description=description,
            tags=tags,
            policy_document=policy_document,
            client_token=client_token,
            account_id=self.account_id,
            partition=self.partition,
        )
        cnw_id = core_network.core_network_id
        self.core_networks[cnw_id] = core_network
        return core_network

    def delete_core_network(self, core_network_id: str) -> CoreNetwork:
        # Check if core network exists
        if core_network_id not in self.core_networks:
            raise ResourceNotFound(core_network_id)
        core_network = self.core_networks.pop(core_network_id)
        core_network.state = "DELETING"
        return core_network

    def tag_resource(self, resource_arn: str, tags: List[Dict[str, Any]]) -> None:
        resource = self._get_resource_from_arn(resource_arn)
        resource.tags.extend(tags)

    def untag_resource(self, resource_arn: str, tag_keys: Optional[List[str]]) -> None:
        resource = self._get_resource_from_arn(resource_arn)
        if tag_keys:
            resource.tags = [tag for tag in resource.tags if tag["Key"] not in tag_keys]

    @paginate(pagination_model=PAGINATION_MODEL)
    def list_core_networks(self) -> List[CoreNetwork]:
        return list(self.core_networks.values())

    def get_core_network(self, core_network_id: str) -> CoreNetwork:
        if core_network_id not in self.core_networks:
            raise ResourceNotFound(core_network_id)
        core_network = self.core_networks[core_network_id]
        return core_network

    @paginate(pagination_model=PAGINATION_MODEL)
    def describe_global_networks(
        self, global_network_ids: List[str]
    ) -> List[GlobalNetwork]:
        queried_global_networks = []
        if not global_network_ids:
            queried_global_networks = list(self.global_networks.values())
        elif isinstance(global_network_ids, str):
            if global_network_ids not in self.global_networks:
                raise ResourceNotFound(global_network_ids)
            queried_global_networks.append(self.global_networks[global_network_ids])
        else:
            for id in global_network_ids:
                if id in self.global_networks:
                    global_network = self.global_networks[id]
                    queried_global_networks.append(global_network)
        return queried_global_networks

    def create_site(
        self,
        global_network_id: str,
        description: Optional[str],
        location: Optional[Dict[str, str]],
        tags: Optional[List[Dict[str, str]]],
    ) -> Site:
        # check if global network exists
        if global_network_id not in self.global_networks:
            raise ResourceNotFound(global_network_id)
        site = Site(
            global_network_id=global_network_id,
            description=description,
            location=location,
            tags=tags,
            account_id=self.account_id,
            partition=self.partition,
        )
        site_id = site.site_id
        self.sites[global_network_id][site_id] = site
        return site

    def delete_site(self, global_network_id: str, site_id: str) -> Site:
        if global_network_id not in self.global_networks:
            raise ResourceNotFound(global_network_id)
        if global_network_id not in self.sites:
            raise ResourceNotFound(site_id)
        gn_sites = self.sites[global_network_id]
        site = gn_sites.pop(site_id)
        site.state = "DELETING"
        return site

    @paginate(pagination_model=PAGINATION_MODEL)
    def get_sites(self, global_network_id: str, site_ids: List[str]) -> List[Site]:
        if global_network_id not in self.global_networks:
            raise ValidationError("Incorrect input.")
        gn_sites = self.sites.get(global_network_id) or {}
        queried = []
        if not site_ids:
            queried = list(gn_sites.values())
        else:
            for id in site_ids:
                if id in gn_sites:
                    q = gn_sites[id]
                    queried.append(q)

        return queried

    def create_link(
        self,
        global_network_id: str,
        description: Optional[str],
        type: Optional[str],
        bandwidth: Dict[str, Any],
        provider: Optional[str],
        site_id: str,
        tags: Optional[List[Dict[str, str]]],
    ) -> Link:
        # check if global network exists
        if global_network_id not in self.global_networks:
            raise ResourceNotFound(global_network_id)
        link = Link(
            global_network_id=global_network_id,
            description=description,
            type=type,
            bandwidth=bandwidth,
            provider=provider,
            site_id=site_id,
            tags=tags,
            account_id=self.account_id,
            partition=self.partition,
        )
        self.links[global_network_id][link.link_id] = link
        return link

    @paginate(pagination_model=PAGINATION_MODEL)
    def get_links(
        self,
        global_network_id: str,
        link_ids: List[str],
        site_id: str,
        type: str,
        provider: str,
    ) -> List[Link]:
        if global_network_id not in self.global_networks:
            raise ValidationError("Incorrect input.")
        # TODO: Implement filtering by site_id, type, provider
        gn_links = self.links.get(global_network_id) or {}
        queried = []
        if not link_ids:
            queried = list(gn_links.values())
        else:
            for id in link_ids:
                if id in gn_links:
                    q = gn_links[id]
                    queried.append(q)
        return queried

    def delete_link(self, global_network_id: str, link_id: str) -> Link:
        try:
            link = self.links[global_network_id].pop(link_id)
        except KeyError:
            raise ResourceNotFound(link_id)
        link.state = "DELETING"
        return link

    def create_device(
        self,
        global_network_id: str,
        aws_location: Optional[Dict[str, str]],
        description: Optional[str],
        type: Optional[str],
        vendor: Optional[str],
        model: Optional[str],
        serial_number: Optional[str],
        location: Optional[Dict[str, str]],
        site_id: Optional[str],
        tags: Optional[List[Dict[str, str]]],
    ) -> Device:
        # check if global network exists
        if global_network_id not in self.global_networks:
            raise ResourceNotFound(global_network_id)
        device = Device(
            global_network_id=global_network_id,
            aws_location=aws_location,
            description=description,
            type=type,
            vendor=vendor,
            model=model,
            serial_number=serial_number,
            location=location,
            site_id=site_id,
            tags=tags,
            account_id=self.account_id,
            partition=self.partition,
        )
        self.devices[global_network_id][device.device_id] = device
        return device

    @paginate(pagination_model=PAGINATION_MODEL)
    def get_devices(
        self, global_network_id: str, device_ids: List[str], site_id: Optional[str]
    ) -> List[Device]:
        if global_network_id not in self.global_networks:
            raise ValidationError("Incorrect input.")
        # TODO: Implement filtering by site_id
        gn_devices = self.devices.get(global_network_id) or {}
        queried = []
        if not device_ids:
            queried = list(gn_devices.values())
        else:
            for id in device_ids:
                if id in gn_devices:
                    q = gn_devices[id]
                    queried.append(q)

        return queried

    def delete_device(self, global_network_id: str, device_id: str) -> Device:
        try:
            device = self.devices[global_network_id].pop(device_id)
        except KeyError:
            raise ResourceNotFound(device_id)
        device.state = "DELETING"
        return device

    def list_tags_for_resource(self, resource_arn: str) -> List[Dict[str, str]]:
        resource = self._get_resource_from_arn(resource_arn)
        tag_list = resource.tags
        return tag_list


networkmanager_backends = BackendDict(
    NetworkManagerBackend,
    "networkmanager",
    use_boto3_regions=False,
    additional_regions=PARTITION_NAMES,
)
