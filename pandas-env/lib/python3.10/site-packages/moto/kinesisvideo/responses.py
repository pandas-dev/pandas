import json

from moto.core.responses import BaseResponse

from .models import KinesisVideoBackend, kinesisvideo_backends


class KinesisVideoResponse(BaseResponse):
    def __init__(self) -> None:
        super().__init__(service_name="kinesisvideo")

    @property
    def kinesisvideo_backend(self) -> KinesisVideoBackend:
        return kinesisvideo_backends[self.current_account][self.region]

    def create_stream(self) -> str:
        device_name = self._get_param("DeviceName")
        stream_name = self._get_param("StreamName")
        media_type = self._get_param("MediaType")
        kms_key_id = self._get_param("KmsKeyId")
        data_retention_in_hours = self._get_int_param("DataRetentionInHours")
        tags = self._get_param("Tags")
        stream_arn = self.kinesisvideo_backend.create_stream(
            device_name=device_name,
            stream_name=stream_name,
            media_type=media_type,
            kms_key_id=kms_key_id,
            data_retention_in_hours=data_retention_in_hours,
            tags=tags,
        )
        return json.dumps(dict(StreamARN=stream_arn))

    def describe_stream(self) -> str:
        stream_name = self._get_param("StreamName")
        stream_arn = self._get_param("StreamARN")
        stream_info = self.kinesisvideo_backend.describe_stream(
            stream_name=stream_name, stream_arn=stream_arn
        )
        return json.dumps(dict(StreamInfo=stream_info))

    def list_streams(self) -> str:
        stream_info_list = self.kinesisvideo_backend.list_streams()
        return json.dumps(dict(StreamInfoList=stream_info_list, NextToken=None))

    def delete_stream(self) -> str:
        stream_arn = self._get_param("StreamARN")
        self.kinesisvideo_backend.delete_stream(stream_arn=stream_arn)
        return json.dumps(dict())

    def get_data_endpoint(self) -> str:
        stream_name = self._get_param("StreamName")
        stream_arn = self._get_param("StreamARN")
        api_name = self._get_param("APIName")
        data_endpoint = self.kinesisvideo_backend.get_data_endpoint(
            stream_name=stream_name, stream_arn=stream_arn, api_name=api_name
        )
        return json.dumps(dict(DataEndpoint=data_endpoint))
