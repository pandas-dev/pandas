"""Handles incoming ivs requests, invokes methods, returns responses."""

import json

from moto.core.responses import BaseResponse

from .models import IVSBackend, ivs_backends


class IVSResponse(BaseResponse):
    """Handler for IVS requests and responses."""

    def __init__(self) -> None:
        super().__init__(service_name="ivs")

    @property
    def ivs_backend(self) -> IVSBackend:
        """Return backend instance specific for this region."""
        return ivs_backends[self.current_account][self.region]

    def create_channel(self) -> str:
        authorized = self._get_param("authorized", False)
        insecure_ingest = self._get_param("insecureIngest", False)
        latency_mode = self._get_param("latencyMode", "LOW")
        name = self._get_param("name")
        preset = self._get_param("preset", "")
        recording_configuration_arn = self._get_param("recordingConfigurationArn", "")
        tags = self._get_param("tags", {})
        channel_type = self._get_param("type", "STANDARD")
        channel, stream_key = self.ivs_backend.create_channel(
            authorized=authorized,
            insecure_ingest=insecure_ingest,
            latency_mode=latency_mode,
            name=name,
            preset=preset,
            recording_configuration_arn=recording_configuration_arn,
            tags=tags,
            channel_type=channel_type,
        )
        return json.dumps(dict(channel=channel, streamKey=stream_key))

    def list_channels(self) -> str:
        filter_by_name = self._get_param("filterByName")
        filter_by_recording_configuration_arn = self._get_param(
            "filterByRecordingConfigurationArn"
        )
        max_results = self._get_param("maxResults")
        next_token = self._get_param("nextToken")
        channels, next_token = self.ivs_backend.list_channels(
            filter_by_name=filter_by_name,
            filter_by_recording_configuration_arn=filter_by_recording_configuration_arn,
            max_results=max_results,
            next_token=next_token,
        )
        return json.dumps(dict(channels=channels, nextToken=next_token))

    def get_channel(self) -> str:
        arn = self._get_param("arn")
        channel = self.ivs_backend.get_channel(
            arn=arn,
        )
        return json.dumps(dict(channel=channel))

    def batch_get_channel(self) -> str:
        arns = self._get_param("arns")
        channels, errors = self.ivs_backend.batch_get_channel(
            arns=arns,
        )
        return json.dumps(dict(channels=channels, errors=errors))

    def update_channel(self) -> str:
        arn = self._get_param("arn")
        authorized = self._get_param("authorized")
        insecure_ingest = self._get_param("insecureIngest")
        latency_mode = self._get_param("latencyMode")
        name = self._get_param("name")
        preset = self._get_param("preset")
        recording_configuration_arn = self._get_param("recordingConfigurationArn")
        channel_type = self._get_param("type")
        channel = self.ivs_backend.update_channel(
            arn=arn,
            authorized=authorized,
            insecure_ingest=insecure_ingest,
            latency_mode=latency_mode,
            name=name,
            preset=preset,
            recording_configuration_arn=recording_configuration_arn,
            channel_type=channel_type,
        )
        return json.dumps(dict(channel=channel))

    def delete_channel(self) -> None:
        arn = self._get_param("arn")
        self.ivs_backend.delete_channel(
            arn=arn,
        )
