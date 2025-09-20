from typing import Tuple

from moto.core.base_backend import BackendDict, BaseBackend
from moto.kinesisvideo.models import KinesisVideoBackend, kinesisvideo_backends
from moto.sts.utils import random_session_token


class KinesisVideoArchivedMediaBackend(BaseBackend):
    @property
    def backend(self) -> KinesisVideoBackend:
        return kinesisvideo_backends[self.account_id][self.region_name]

    def _get_streaming_url(
        self, stream_name: str, stream_arn: str, api_name: str
    ) -> str:
        stream = self.backend._get_stream(stream_name, stream_arn)
        data_endpoint = stream.get_data_endpoint(api_name)
        session_token = random_session_token()
        api_to_relative_path = {
            "GET_HLS_STREAMING_SESSION_URL": "/hls/v1/getHLSMasterPlaylist.m3u8",
            "GET_DASH_STREAMING_SESSION_URL": "/dash/v1/getDASHManifest.mpd",
        }
        relative_path = api_to_relative_path[api_name]
        return f"{data_endpoint}{relative_path}?SessionToken={session_token}"

    def get_hls_streaming_session_url(self, stream_name: str, stream_arn: str) -> str:
        # Ignore option parameters as the format of hls_url doesn't depend on them
        api_name = "GET_HLS_STREAMING_SESSION_URL"
        return self._get_streaming_url(stream_name, stream_arn, api_name)

    def get_dash_streaming_session_url(self, stream_name: str, stream_arn: str) -> str:
        # Ignore option parameters as the format of hls_url doesn't depend on them
        api_name = "GET_DASH_STREAMING_SESSION_URL"
        return self._get_streaming_url(stream_name, stream_arn, api_name)

    def get_clip(self, stream_name: str, stream_arn: str) -> Tuple[str, bytes]:
        self.backend._get_stream(stream_name, stream_arn)
        content_type = "video/mp4"  # Fixed content_type as it depends on input stream
        payload = b"sample-mp4-video"
        return content_type, payload


kinesisvideoarchivedmedia_backends = BackendDict(
    KinesisVideoArchivedMediaBackend, "kinesis-video-archived-media"
)
