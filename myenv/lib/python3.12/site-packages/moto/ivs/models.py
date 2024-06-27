"""IVSBackend class with methods for supported APIs."""

from typing import Any, Dict, List, Optional, Tuple

from moto.core.base_backend import BackendDict, BaseBackend
from moto.ivs.exceptions import ResourceNotFoundException
from moto.moto_api._internal import mock_random
from moto.utilities.paginator import paginate
from moto.utilities.utils import get_partition


class IVSBackend(BaseBackend):
    """Implementation of IVS APIs."""

    PAGINATION_MODEL = {
        "list_channels": {
            "input_token": "next_token",
            "limit_key": "max_results",
            "limit_default": 100,
            "unique_attribute": "arn",
        },
    }

    def __init__(self, region_name: str, account_id: str):
        super().__init__(region_name, account_id)
        self.channels: List[Dict[str, Any]] = []

    def create_channel(
        self,
        authorized: bool,
        insecure_ingest: bool,
        latency_mode: str,
        name: str,
        preset: str,
        recording_configuration_arn: str,
        tags: Dict[str, str],
        channel_type: str,
    ) -> Tuple[Dict[str, Any], Dict[str, Any]]:
        channel_id = mock_random.get_random_string(12)
        channel_arn = f"arn:{get_partition(self.region_name)}:ivs:{self.region_name}:{self.account_id}:channel/{channel_id}"
        channel = {
            "arn": channel_arn,
            "authorized": authorized,
            "ingestEndpoint": "ingest.example.com",
            "insecureIngest": insecure_ingest,
            "latencyMode": latency_mode,
            "name": name,
            "playbackUrl": f"https://playback.example.com/{self.region_name}.{self.account_id}.{channel_id}.m3u8",
            "preset": preset,
            "recordingConfigurationArn": recording_configuration_arn,
            "tags": tags,
            "type": channel_type,
        }
        self.channels.append(channel)
        stream_key_id = mock_random.get_random_string(12)
        stream_key_arn = f"arn:{get_partition(self.region_name)}:ivs:{self.region_name}:{self.account_id}:stream-key/{stream_key_id}"
        stream_key = {
            "arn": stream_key_arn,
            "channelArn": channel_arn,
            "tags": tags,
            "value": f"sk_{self.region_name}_{mock_random.token_urlsafe(32)}",
        }
        return channel, stream_key

    @paginate(pagination_model=PAGINATION_MODEL)  # type: ignore[misc]
    def list_channels(
        self,
        filter_by_name: Optional[str],
        filter_by_recording_configuration_arn: Optional[str],
    ) -> List[Dict[str, Any]]:
        if filter_by_name is not None:
            channels = [
                channel
                for channel in self.channels
                if channel["name"] == filter_by_name
            ]
        elif filter_by_recording_configuration_arn is not None:
            channels = [
                channel
                for channel in self.channels
                if channel["recordingConfigurationArn"]
                == filter_by_recording_configuration_arn
            ]
        else:
            channels = self.channels
        return channels

    def _find_channel(self, arn: str) -> Dict[str, Any]:
        try:
            return next(channel for channel in self.channels if channel["arn"] == arn)
        except StopIteration:
            raise ResourceNotFoundException(f"Resource: {arn} not found")

    def get_channel(self, arn: str) -> Dict[str, Any]:
        return self._find_channel(arn)

    def batch_get_channel(
        self, arns: List[str]
    ) -> Tuple[List[Dict[str, Any]], List[Dict[str, str]]]:
        return [channel for channel in self.channels if channel["arn"] in arns], []

    def update_channel(
        self,
        arn: str,
        authorized: Optional[bool],
        insecure_ingest: Optional[bool],
        latency_mode: Optional[str],
        name: Optional[str],
        preset: Optional[str],
        recording_configuration_arn: Optional[str],
        channel_type: Optional[str],
    ) -> Dict[str, Any]:
        channel = self._find_channel(arn)
        if authorized is not None:
            channel["authorized"] = authorized
        if insecure_ingest is not None:
            channel["insecureIngest"] = insecure_ingest
        if latency_mode is not None:
            channel["latencyMode"] = latency_mode
        if name is not None:
            channel["name"] = name
        if preset is not None:
            channel["preset"] = preset
        if recording_configuration_arn is not None:
            channel["recordingConfigurationArn"] = recording_configuration_arn
        if channel_type is not None:
            channel["type"] = channel_type
        return channel

    def delete_channel(self, arn: str) -> None:
        self._find_channel(arn)
        self.channels = [channel for channel in self.channels if channel["arn"] != arn]


ivs_backends = BackendDict(IVSBackend, "ivs")
