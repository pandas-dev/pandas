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
        return json.dumps(dict(GlobalNetwork=global_network.to_dict()))

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
        return json.dumps(dict(CoreNetwork=core_network.to_dict()))

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
