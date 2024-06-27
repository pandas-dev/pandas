import json
from typing import Dict, Tuple

from moto.core.responses import BaseResponse

from .models import KinesisVideoArchivedMediaBackend, kinesisvideoarchivedmedia_backends


class KinesisVideoArchivedMediaResponse(BaseResponse):
    def __init__(self) -> None:
        super().__init__(service_name="kinesis-video-archived-media")

    @property
    def kinesisvideoarchivedmedia_backend(self) -> KinesisVideoArchivedMediaBackend:
        return kinesisvideoarchivedmedia_backends[self.current_account][self.region]

    def get_hls_streaming_session_url(self) -> str:
        stream_name = self._get_param("StreamName")
        stream_arn = self._get_param("StreamARN")
        hls_streaming_session_url = (
            self.kinesisvideoarchivedmedia_backend.get_hls_streaming_session_url(
                stream_name=stream_name, stream_arn=stream_arn
            )
        )
        return json.dumps(dict(HLSStreamingSessionURL=hls_streaming_session_url))

    def get_dash_streaming_session_url(self) -> str:
        stream_name = self._get_param("StreamName")
        stream_arn = self._get_param("StreamARN")
        dash_streaming_session_url = (
            self.kinesisvideoarchivedmedia_backend.get_dash_streaming_session_url(
                stream_name=stream_name, stream_arn=stream_arn
            )
        )
        return json.dumps(dict(DASHStreamingSessionURL=dash_streaming_session_url))

    def get_clip(self) -> Tuple[bytes, Dict[str, str]]:
        stream_name = self._get_param("StreamName")
        stream_arn = self._get_param("StreamARN")
        content_type, payload = self.kinesisvideoarchivedmedia_backend.get_clip(
            stream_name=stream_name, stream_arn=stream_arn
        )
        new_headers = {"Content-Type": content_type}
        return payload, new_headers
