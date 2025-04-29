"""Handles incoming networkmanager requests, invokes methods, returns responses."""

import json
from urllib.parse import unquote

from moto.core.common_types import TYPE_RESPONSE
from moto.core.responses import BaseResponse

from .models import NetworkManagerBackend, networkmanager_backends


class NetworkManagerResponse(BaseResponse):
    """Handler for NetworkManager requests and responses."""

    def __init__(self) -> None:
        super().__init__(service_name="networkmanager")

    @property
    def networkmanager_backend(self) -> NetworkManagerBackend:
        return networkmanager_backends[self.current_account][self.partition]

    def create_global_network(self) -> str:
        params = json.loads(self.body)
        description = params.get("Description")
        tags = params.get("Tags")
        global_network = self.networkmanager_backend.create_global_network(
            description=description,
            tags=tags,
        )
        resp_dict = dict(GlobalNetwork=global_network.to_dict())
        self.networkmanager_backend.update_resource_state(
            global_network.global_network_arn, "AVAILABLE"
        )
        return json.dumps(resp_dict)

    def create_core_network(self) -> str:
        params = json.loads(self.body)
        global_network_id = params.get("GlobalNetworkId")
        description = params.get("Description")
        tags = params.get("Tags")
        policy_document = params.get("PolicyDocument")
        client_token = params.get("ClientToken")
        core_network = self.networkmanager_backend.create_core_network(
            global_network_id=global_network_id,
            description=description,
            tags=tags,
            policy_document=policy_document,
            client_token=client_token,
        )
        resp_dict = dict(CoreNetwork=core_network.to_dict())
        self.networkmanager_backend.update_resource_state(
            core_network.core_network_arn, "AVAILABLE"
        )
        return json.dumps(resp_dict)

    def delete_core_network(self) -> str:
        core_network_id = unquote(self.path.split("/")[-1])
        core_network = self.networkmanager_backend.delete_core_network(
            core_network_id=core_network_id,
        )
        return json.dumps(dict(CoreNetwork=core_network.to_dict()))

    def tag_resource(self) -> TYPE_RESPONSE:
        params = json.loads(self.body)
        tags = params.get("Tags")
        resource_arn = unquote(self.path.split("/tags/")[-1])

        self.networkmanager_backend.tag_resource(
            resource_arn=resource_arn,
            tags=tags,
        )
        return 200, {}, json.dumps({})

    def untag_resource(self) -> TYPE_RESPONSE:
        params = self._get_params()
        tag_keys = params.get("tagKeys")
        resource_arn = unquote(self.path.split("/tags/")[-1])
        self.networkmanager_backend.untag_resource(
            resource_arn=resource_arn,
            tag_keys=tag_keys,
        )
        return 200, {}, json.dumps({})

    def list_core_networks(self) -> str:
        params = self._get_params()
        max_results = params.get("maxResults")
        next_token = params.get("nextToken")
        core_networks, next_token = self.networkmanager_backend.list_core_networks(
            max_results=max_results,
            next_token=next_token,
        )
        list_core_networks = [core_network.to_dict() for core_network in core_networks]
        return json.dumps(dict(CoreNetworks=list_core_networks, NextToken=next_token))

    def get_core_network(self) -> str:
        core_network_id = unquote(self.path.split("/")[-1])
        core_network = self.networkmanager_backend.get_core_network(
            core_network_id=core_network_id,
        )
        return json.dumps(dict(CoreNetwork=core_network.to_dict()))

    def describe_global_networks(self) -> str:
        params = self._get_params()
        global_network_ids = params.get("globalNetworkIds")
        max_results = params.get("maxResults")
        next_token = params.get("nextToken")
        global_networks, next_token = (
            self.networkmanager_backend.describe_global_networks(
                global_network_ids=global_network_ids,
                max_results=max_results,
                next_token=next_token,
            )
        )
        list_global_networks = [
            global_network.to_dict() for global_network in global_networks
        ]
        return json.dumps(
            dict(GlobalNetworks=list_global_networks, nextToken=next_token)
        )

    def create_site(self) -> str:
        params = json.loads(self.body)
        global_network_id = unquote(self.path.split("/")[-2])
        description = params.get("Description")
        location = params.get("Location")
        tags = params.get("Tags")
        site = self.networkmanager_backend.create_site(
            global_network_id=global_network_id,
            description=description,
            location=location,
            tags=tags,
        )
        resp_dict = dict(Site=site.to_dict())
        self.networkmanager_backend.update_resource_state(site.site_arn, "AVAILABLE")
        return json.dumps(resp_dict)

    def delete_site(self) -> str:
        global_network_id = unquote(self.path.split("/")[-3])
        site_id = unquote(self.path.split("/")[-1])
        site = self.networkmanager_backend.delete_site(
            global_network_id=global_network_id,
            site_id=site_id,
        )
        return json.dumps(dict(Site=site.to_dict()))

    def get_sites(self) -> str:
        params = self._get_params()
        global_network_id = unquote(self.path.split("/")[-2])
        site_ids = self.querystring.get("siteIds")
        max_results = params.get("MaxResults")
        next_token = params.get("NextToken")
        sites, next_token = self.networkmanager_backend.get_sites(
            global_network_id=global_network_id,
            site_ids=site_ids,
            max_results=max_results,
            next_token=next_token,
        )
        list_sites = [site.to_dict() for site in sites]
        return json.dumps(dict(Sites=list_sites, nextToken=next_token))

    def create_link(self) -> str:
        params = json.loads(self.body)
        global_network_id = unquote(self.path.split("/")[-2])
        description = params.get("Description")
        type = params.get("Type")
        bandwidth = params.get("Bandwidth")
        provider = params.get("Provider")
        site_id = params.get("SiteId")
        tags = params.get("Tags")
        link = self.networkmanager_backend.create_link(
            global_network_id=global_network_id,
            description=description,
            type=type,
            bandwidth=bandwidth,
            provider=provider,
            site_id=site_id,
            tags=tags,
        )
        resp_dict = dict(Link=link.to_dict())
        self.networkmanager_backend.update_resource_state(link.link_arn, "AVAILABLE")
        return json.dumps(resp_dict)

    def get_links(self) -> str:
        params = self._get_params()
        global_network_id = unquote(self.path.split("/")[-2])
        link_ids = self.querystring.get("linkIds")
        site_id = params.get("SiteId")
        type = params.get("Type")
        provider = params.get("Provider")
        max_results = params.get("MaxResults")
        next_token = params.get("NextToken")
        links, next_token = self.networkmanager_backend.get_links(
            global_network_id=global_network_id,
            link_ids=link_ids,
            site_id=site_id,
            type=type,
            provider=provider,
            max_results=max_results,
            next_token=next_token,
        )
        list_links = [link.to_dict() for link in links]
        return json.dumps(dict(Links=list_links, nextToken=next_token))

    def delete_link(self) -> str:
        global_network_id = unquote(self.path.split("/")[-3])
        link_id = unquote(self.path.split("/")[-1])
        link = self.networkmanager_backend.delete_link(
            global_network_id=global_network_id,
            link_id=link_id,
        )
        return json.dumps(dict(Link=link.to_dict()))

    def create_device(self) -> str:
        params = json.loads(self.body)
        global_network_id = unquote(self.path.split("/")[-2])
        aws_location = params.get("AWSLocation")
        description = params.get("Description")
        type = params.get("Type")
        vendor = params.get("Vendor")
        model = params.get("Model")
        serial_number = params.get("SerialNumber")
        location = params.get("Location")
        site_id = params.get("SiteId")
        tags = params.get("Tags")
        device = self.networkmanager_backend.create_device(
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
        )
        resp_dict = dict(Device=device.to_dict())
        self.networkmanager_backend.update_resource_state(
            device.device_arn, "AVAILABLE"
        )
        return json.dumps(resp_dict)

    def get_devices(self) -> str:
        params = self._get_params()
        global_network_id = unquote(self.path.split("/")[-2])
        device_ids = self.querystring.get("deviceIds")
        site_id = params.get("SiteId")
        max_results = params.get("MaxResults")
        next_token = params.get("NextToken")
        devices, next_token = self.networkmanager_backend.get_devices(
            global_network_id=global_network_id,
            device_ids=device_ids,
            site_id=site_id,
            max_results=max_results,
            next_token=next_token,
        )
        list_devices = [device.to_dict() for device in devices]
        return json.dumps(dict(Devices=list_devices, nextToken=next_token))

    def delete_device(self) -> str:
        global_network_id = unquote(self.path.split("/")[-3])
        device_id = unquote(self.path.split("/")[-1])
        device = self.networkmanager_backend.delete_device(
            global_network_id=global_network_id,
            device_id=device_id,
        )
        return json.dumps(dict(Device=device.to_dict()))

    def list_tags_for_resource(self) -> str:
        resource_arn = unquote(self.path.split("/tags/")[-1])
        tag_list = self.networkmanager_backend.list_tags_for_resource(
            resource_arn=resource_arn,
        )
        return json.dumps(dict(TagList=tag_list))
