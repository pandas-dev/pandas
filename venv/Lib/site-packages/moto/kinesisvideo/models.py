from typing import Any, Dict, List

from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import BaseModel
from moto.core.utils import utcnow
from moto.moto_api._internal import mock_random as random
from moto.utilities.utils import get_partition

from .exceptions import ResourceInUseException, ResourceNotFoundException


class Stream(BaseModel):
    def __init__(
        self,
        account_id: str,
        region_name: str,
        device_name: str,
        stream_name: str,
        media_type: str,
        kms_key_id: str,
        data_retention_in_hours: int,
        tags: Dict[str, str],
    ):
        self.region_name = region_name
        self.stream_name = stream_name
        self.device_name = device_name
        self.media_type = media_type
        self.kms_key_id = kms_key_id
        self.data_retention_in_hours = data_retention_in_hours
        self.tags = tags
        self.status = "ACTIVE"
        self.version = random.get_random_string(include_digits=False, lower_case=True)
        self.creation_time = utcnow()
        stream_arn = f"arn:{get_partition(region_name)}:kinesisvideo:{region_name}:{account_id}:stream/{stream_name}/1598784211076"
        self.data_endpoint_number = random.get_random_hex()
        self.arn = stream_arn

    def get_data_endpoint(self, api_name: str) -> str:
        data_endpoint_prefix = "s-" if api_name in ("PUT_MEDIA", "GET_MEDIA") else "b-"
        return f"https://{data_endpoint_prefix}{self.data_endpoint_number}.kinesisvideo.{self.region_name}.amazonaws.com"

    def to_dict(self) -> Dict[str, Any]:
        return {
            "DeviceName": self.device_name,
            "StreamName": self.stream_name,
            "StreamARN": self.arn,
            "MediaType": self.media_type,
            "KmsKeyId": self.kms_key_id,
            "Version": self.version,
            "Status": self.status,
            "CreationTime": self.creation_time.isoformat(),
            "DataRetentionInHours": self.data_retention_in_hours,
        }


class KinesisVideoBackend(BaseBackend):
    def __init__(self, region_name: str, account_id: str):
        super().__init__(region_name, account_id)
        self.streams: Dict[str, Stream] = {}

    def create_stream(
        self,
        device_name: str,
        stream_name: str,
        media_type: str,
        kms_key_id: str,
        data_retention_in_hours: int,
        tags: Dict[str, str],
    ) -> str:
        streams = [_ for _ in self.streams.values() if _.stream_name == stream_name]
        if len(streams) > 0:
            raise ResourceInUseException(f"The stream {stream_name} already exists.")
        stream = Stream(
            self.account_id,
            self.region_name,
            device_name,
            stream_name,
            media_type,
            kms_key_id,
            data_retention_in_hours,
            tags,
        )
        self.streams[stream.arn] = stream
        return stream.arn

    def _get_stream(self, stream_name: str, stream_arn: str) -> Stream:
        if stream_name:
            streams = [_ for _ in self.streams.values() if _.stream_name == stream_name]
            if len(streams) == 0:
                raise ResourceNotFoundException()
            stream = streams[0]
        elif stream_arn:
            stream = self.streams.get(stream_arn)
            if stream is None:
                raise ResourceNotFoundException()
        return stream

    def describe_stream(self, stream_name: str, stream_arn: str) -> Dict[str, Any]:
        stream = self._get_stream(stream_name, stream_arn)
        return stream.to_dict()

    def list_streams(self) -> List[Dict[str, Any]]:
        """
        Pagination and the StreamNameCondition-parameter are not yet implemented
        """
        return [_.to_dict() for _ in self.streams.values()]

    def delete_stream(self, stream_arn: str) -> None:
        """
        The CurrentVersion-parameter is not yet implemented
        """
        stream = self.streams.get(stream_arn)
        if stream is None:
            raise ResourceNotFoundException()
        del self.streams[stream_arn]

    def get_data_endpoint(
        self, stream_name: str, stream_arn: str, api_name: str
    ) -> str:
        stream = self._get_stream(stream_name, stream_arn)
        return stream.get_data_endpoint(api_name)


kinesisvideo_backends = BackendDict(KinesisVideoBackend, "kinesisvideo")
