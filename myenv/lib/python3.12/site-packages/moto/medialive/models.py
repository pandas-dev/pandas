from collections import OrderedDict
from typing import Any, Dict, List, Optional

from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import BaseModel
from moto.moto_api._internal import mock_random
from moto.utilities.utils import get_partition


class Input(BaseModel):
    def __init__(self, **kwargs: Any):
        self.arn = kwargs.get("arn")
        self.attached_channels = kwargs.get("attached_channels", [])
        self.destinations = kwargs.get("destinations", [])
        self.input_id = kwargs.get("input_id")
        self.input_class = kwargs.get("input_class", "STANDARD")
        self.input_devices = kwargs.get("input_devices", [])
        self.input_source_type = kwargs.get("input_source_type", "STATIC")
        self.media_connect_flows = kwargs.get("media_connect_flows", [])
        self.name = kwargs.get("name")
        self.role_arn = kwargs.get("role_arn")
        self.security_groups = kwargs.get("security_groups", [])
        self.sources = kwargs.get("sources", [])
        # Possible states: 'CREATING'|'DETACHED'|'ATTACHED'|'DELETING'|'DELETED'
        self.state = kwargs.get("state")
        self.tags = kwargs.get("tags")
        self.input_type = kwargs.get("input_type")

    def to_dict(self) -> Dict[str, Any]:
        return {
            "arn": self.arn,
            "attachedChannels": self.attached_channels,
            "destinations": self.destinations,
            "id": self.input_id,
            "inputClass": self.input_class,
            "inputDevices": self.input_devices,
            "inputSourceType": self.input_source_type,
            "mediaConnectFlows": self.media_connect_flows,
            "name": self.name,
            "roleArn": self.role_arn,
            "securityGroups": self.security_groups,
            "sources": self.sources,
            "state": self.state,
            "tags": self.tags,
            "type": self.input_type,
        }

    def _resolve_transient_states(self) -> None:
        # Resolve transient states before second call
        # (to simulate AWS taking its sweet time with these things)
        if self.state in ["CREATING"]:
            self.state = "DETACHED"  # or ATTACHED
        elif self.state == "DELETING":
            self.state = "DELETED"


class Channel(BaseModel):
    def __init__(self, **kwargs: Any):
        self.arn = kwargs.get("arn")
        self.cdi_input_specification = kwargs.get("cdi_input_specification")
        self.channel_class = kwargs.get("channel_class", "STANDARD")
        self.destinations = kwargs.get("destinations")
        self.egress_endpoints = kwargs.get("egress_endpoints", [])
        self.encoder_settings = kwargs.get("encoder_settings")
        self.channel_id = kwargs.get("channel_id")
        self.input_attachments = kwargs.get("input_attachments")
        self.input_specification = kwargs.get("input_specification")
        self.log_level = kwargs.get("log_level")
        self.name = kwargs.get("name")
        self.pipeline_details = kwargs.get("pipeline_details", [])
        self.role_arn = kwargs.get("role_arn")
        self.state = kwargs.get("state")
        self.tags = kwargs.get("tags")
        self._previous_state = None

    def to_dict(self, exclude: Optional[List[str]] = None) -> Dict[str, Any]:
        data = {
            "arn": self.arn,
            "cdiInputSpecification": self.cdi_input_specification,
            "channelClass": self.channel_class,
            "destinations": self.destinations,
            "egressEndpoints": self.egress_endpoints,
            "encoderSettings": self.encoder_settings,
            "id": self.channel_id,
            "inputAttachments": self.input_attachments,
            "inputSpecification": self.input_specification,
            "logLevel": self.log_level,
            "name": self.name,
            "pipelineDetails": self.pipeline_details,
            "pipelinesRunningCount": 1
            if self.channel_class == "SINGLE_PIPELINE"
            else 2,
            "roleArn": self.role_arn,
            "state": self.state,
            "tags": self.tags,
        }
        if exclude:
            for key in exclude:
                del data[key]
        return data

    def _resolve_transient_states(self) -> None:
        # Resolve transient states before second call
        # (to simulate AWS taking its sweet time with these things)
        if self.state in ["CREATING", "STOPPING"]:
            self.state = "IDLE"
        elif self.state == "STARTING":
            self.state = "RUNNING"
        elif self.state == "DELETING":
            self.state = "DELETED"
        elif self.state == "UPDATING":
            self.state = self._previous_state
            self._previous_state = None


class MediaLiveBackend(BaseBackend):
    def __init__(self, region_name: str, account_id: str):
        super().__init__(region_name, account_id)
        self._channels: Dict[str, Channel] = OrderedDict()
        self._inputs: Dict[str, Input] = OrderedDict()

    def create_channel(
        self,
        cdi_input_specification: Dict[str, Any],
        channel_class: str,
        destinations: List[Dict[str, Any]],
        encoder_settings: Dict[str, Any],
        input_attachments: List[Dict[str, Any]],
        input_specification: Dict[str, str],
        log_level: str,
        name: str,
        role_arn: str,
        tags: Dict[str, str],
    ) -> Channel:
        """
        The RequestID and Reserved parameters are not yet implemented
        """
        channel_id = mock_random.uuid4().hex
        arn = f"arn:{get_partition(self.region_name)}:medialive:channel:{channel_id}"
        channel = Channel(
            arn=arn,
            cdi_input_specification=cdi_input_specification,
            channel_class=channel_class or "STANDARD",
            destinations=destinations,
            egress_endpoints=[],
            encoder_settings=encoder_settings,
            channel_id=channel_id,
            input_attachments=input_attachments,
            input_specification=input_specification,
            log_level=log_level,
            name=name,
            pipeline_details=[],
            role_arn=role_arn,
            state="CREATING",
            tags=tags,
        )
        self._channels[channel_id] = channel
        return channel

    def list_channels(self, max_results: Optional[int]) -> List[Dict[str, Any]]:
        """
        Pagination is not yet implemented
        """
        channels = list(self._channels.values())
        if max_results is not None:
            channels = channels[:max_results]
        return [
            c.to_dict(exclude=["encoderSettings", "pipelineDetails"]) for c in channels
        ]

    def describe_channel(self, channel_id: str) -> Channel:
        channel = self._channels[channel_id]
        channel._resolve_transient_states()
        return channel

    def delete_channel(self, channel_id: str) -> Channel:
        channel = self._channels[channel_id]
        channel.state = "DELETING"
        return channel

    def start_channel(self, channel_id: str) -> Channel:
        channel = self._channels[channel_id]
        channel.state = "STARTING"
        return channel

    def stop_channel(self, channel_id: str) -> Channel:
        channel = self._channels[channel_id]
        channel.state = "STOPPING"
        return channel

    def update_channel(
        self,
        channel_id: str,
        cdi_input_specification: Dict[str, str],
        destinations: List[Dict[str, Any]],
        encoder_settings: Dict[str, Any],
        input_attachments: List[Dict[str, Any]],
        input_specification: Dict[str, str],
        log_level: str,
        name: str,
        role_arn: str,
    ) -> Channel:
        channel = self._channels[channel_id]
        channel.cdi_input_specification = cdi_input_specification
        channel.destinations = destinations
        channel.encoder_settings = encoder_settings
        channel.input_attachments = input_attachments
        channel.input_specification = input_specification
        channel.log_level = log_level
        channel.name = name
        channel.role_arn = role_arn

        channel._resolve_transient_states()
        channel._previous_state = channel.state
        channel.state = "UPDATING"

        return channel

    def create_input(
        self,
        destinations: List[Dict[str, str]],
        input_devices: List[Dict[str, str]],
        input_security_groups: List[str],
        media_connect_flows: List[Dict[str, str]],
        name: str,
        role_arn: str,
        sources: List[Dict[str, str]],
        tags: Dict[str, str],
        input_type: str,
    ) -> Input:
        """
        The VPC and RequestId parameters are not yet implemented
        """
        input_id = mock_random.uuid4().hex
        arn = f"arn:{get_partition(self.region_name)}:medialive:input:{input_id}"
        a_input = Input(
            arn=arn,
            input_id=input_id,
            destinations=destinations,
            input_devices=input_devices,
            input_security_groups=input_security_groups,
            media_connect_flows=media_connect_flows,
            name=name,
            role_arn=role_arn,
            sources=sources,
            tags=tags,
            input_type=input_type,
            state="CREATING",
        )
        self._inputs[input_id] = a_input
        return a_input

    def describe_input(self, input_id: str) -> Input:
        a_input = self._inputs[input_id]
        a_input._resolve_transient_states()
        return a_input

    def list_inputs(self, max_results: Optional[int]) -> List[Dict[str, Any]]:
        """
        Pagination is not yet implemented
        """
        inputs = list(self._inputs.values())
        if max_results is not None:
            inputs = inputs[:max_results]
        return [i.to_dict() for i in inputs]

    def delete_input(self, input_id: str) -> None:
        a_input = self._inputs[input_id]
        a_input.state = "DELETING"

    def update_input(
        self,
        destinations: List[Dict[str, str]],
        input_devices: List[Dict[str, str]],
        input_id: str,
        input_security_groups: List[str],
        media_connect_flows: List[Dict[str, str]],
        name: str,
        role_arn: str,
        sources: List[Dict[str, str]],
    ) -> Input:
        a_input = self._inputs[input_id]
        a_input.destinations = destinations
        a_input.input_devices = input_devices
        a_input.security_groups = input_security_groups
        a_input.media_connect_flows = media_connect_flows
        a_input.name = name
        a_input.role_arn = role_arn
        a_input.sources = sources
        return a_input


medialive_backends = BackendDict(MediaLiveBackend, "medialive")
